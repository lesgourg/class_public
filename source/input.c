/** @file input.c Documented input module.
 *
 * Julien Lesgourgues, 27.08.2010
 * * internal organization of the module structured and improved by Nils Schoeneberg and Matteo Lucca, 07.03.2019
 *
 */

#include "input.h"

/* The input module fills variables belonging to the structures of
   essentially all other modules. Thus we need to include all the
   headers. New in v3.0: These #include fit better here than in
   input.h, to avoid complictaed dependencies slowing down
   compilation. */

#include "quadrature.h"
#include "background.h"
#include "thermodynamics.h"
#include "perturbations.h"
#include "transfer.h"
#include "primordial.h"
#include "harmonic.h"
#include "fourier.h"
#include "lensing.h"
#include "distortions.h"
#include "output.h"

/**
 * Initialize input parameters from external file.
 *
 * @param argc    Input: Number of command line arguments
 * @param argv    Input: Command line argument strings
 * @param ppr     Input: pointer to precision structure
 * @param pba     Input: pointer to background structure
 * @param pth     Input: pointer to thermodynamics structure
 * @param ppt     Input: pointer to perturbation structure
 * @param ptr     Input: pointer to transfer structure
 * @param ppm     Input: pointer to primordial structure
 * @param phr     Input: pointer to harmonic structure
 * @param pfo     Input: pointer to fourier structure
 * @param ple     Input: pointer to lensing structure
 * @param psd     Input: pointer to distorsion structure
 * @param pop     Input: pointer to output structure
 * @param errmsg  Input/Output: Error message
 * @return the error status
 */

int input_init(int argc,
               char **argv,
               struct precision * ppr,
               struct background *pba,
               struct thermodynamics *pth,
               struct perturbations *ppt,
               struct transfer *ptr,
               struct primordial *ppm,
               struct harmonic *phr,
               struct fourier * pfo,
               struct lensing *ple,
               struct distortions *psd,
               struct output *pop,
               ErrorMsg errmsg){

  /** Summary: */

  /** Define local variables */
  struct file_content fc;        // Structure with all parameters

  /** Find and read input file */
  class_call(input_find_file(argc,
                             argv,
                             &fc,
                             errmsg),
             errmsg,
             errmsg);

  /** Initialize all parameters given the input 'file_content' structure.
      If its size is null, all parameters take their default values. */
  class_call(input_read_from_file(&fc,ppr,pba,pth,ppt,ptr,ppm,phr,pfo,ple,psd,pop,
                                  errmsg),
             errmsg,
             errmsg);

  /** Free local struture */
  class_call(parser_free(&fc),
             errmsg,
             errmsg);

  return _SUCCESS_;

}


/**
 * Find and read external file (xxx.ini or xxx.pre) containing the input
 * parameters. All data is stored in the local structure 'file_content'.
 *
 * @param argc    Input: Number of command line arguments
 * @param argv    Input: Command line argument strings
 * @param fc      Output: file_content structure
 * @param errmsg  Input/Output: Error message
 * @return the error status
 */

int input_find_file(int argc,
                    char **argv,
                    struct file_content * fc,
                    ErrorMsg errmsg){

  /** Summary: */

  /** Define local variables */
  struct file_content fc_input;       // Temporary structure with all input parameters
  struct file_content fc_precision;   // Temporary structure with all precision parameters
  struct file_content * pfc_input;    // Pointer to either fc_root or fc_inputroot
  struct file_content fc_setroot;     // Temporary structure for setroot

  int i;
  char extension[5];
  char input_file[_ARGUMENT_LENGTH_MAX_];
  char precision_file[_ARGUMENT_LENGTH_MAX_];

  pfc_input = &fc_input;

  /** Initialize the two file_content structures (for input parameters and
      precision parameters) to some null content. If no arguments are passed,
      they will remain null and inform input_init that all parameters take
      default values. */
  fc->size = 0;
  fc_input.size = 0;
  fc_precision.size = 0;
  input_file[0]='\0';
  precision_file[0]='\0';

  /** If some arguments are passed, identify eventually some 'xxx.ini' and
      'xxx.pre' files, and store their name. */
  if (argc > 1) {
    for (i=1; i<argc; i++) {
      strncpy(extension,(argv[i]+strlen(argv[i])-4),4);
      extension[4]='\0';
      if (strcmp(extension,".ini") == 0) {
        class_test(input_file[0] != '\0',
                   errmsg,
                   "You have passed more than one input file with extension '.ini', choose one.");
        strcpy(input_file,argv[i]);
      }
      else if (strcmp(extension,".pre") == 0) {
        class_test(precision_file[0] != '\0',
                   errmsg,
                   "You have passed more than one precision with extension '.pre', choose one.");
        strcpy(precision_file,argv[i]);
      }
      else {
        fprintf(stdout,"Warning: the file '%s' has an extension different from .ini and .pre, so it has been ignored\n",argv[i]);
      }
    }
  }

  /** If there is an 'xxx.ini' file, read it and store its content. */
  if (input_file[0] != '\0'){
    class_call(parser_read_file(input_file,&fc_input,errmsg),
               errmsg,
               errmsg);

    /* *
     * Set correctly the pfc_input within the set_root function,
     * including the adjusted 'root' field.
     * It is either set to a pointer to fc_input,
     * or when a new filecontent is created with root appended,
     * it is set as the pointer to fc_setroot.
     * */
    class_call(input_set_root(input_file,&pfc_input,&fc_setroot,errmsg),
               errmsg,
               errmsg);
  }

  /** If there is an 'xxx.pre' file, read it and store its content. */
  if (precision_file[0] != '\0'){
    class_call(parser_read_file(precision_file,&fc_precision,errmsg),
               errmsg,
               errmsg);
  }

  /** If one or two files were read, merge their contents in a single 'file_content'
      structure. */
  if ((input_file[0]!='\0') || (precision_file[0]!='\0')){
    class_call(parser_cat(pfc_input,
                          &fc_precision,
                          fc,
                          errmsg),
               errmsg,
               errmsg);
  }

  /** Free local strutures */
  class_call(parser_free(pfc_input),
             errmsg,
             errmsg);
  class_call(parser_free(&fc_precision),
             errmsg,
             errmsg);

  return _SUCCESS_;

}

/**
 * Sets the 'root' variable in the input file content (this will be
 * the beginning of the name of all output files for the current CLASS
 * run)
 *
 * @param input_file      Input: filename of the input file
 * @param ppfc_input      Input/Output: pointer to (pointer to input file structure)
 * @param pfc_setroot     Input: pointer to an allocated temporary file content that will be used here
 * @param errmsg          Input/Output: the error message
 * @return the error status
 */

int input_set_root(char* input_file,
                   struct file_content** ppfc_input,
                   struct file_content * pfc_setroot,
                   ErrorMsg errmsg) {

  /** Define local variables */
  int flag1, filenum, iextens;
  int index_root_in_fc_input = -1;
  int overwrite_root;
  int found_filenum;

  /* The filename of the output root INCLUDING 'output/' */
  FileArg outfname;

  /* Temporary variables for parser reading etc. */
  char tmp_file[_ARGUMENT_LENGTH_MAX_+26]; // 26 is enough to extend the file name [...] with the
                                           // characters "output/[...]%02d_parameters.ini" (as done below)
  struct file_content fc_root;             // Temporary structure with only the root name

  FileArg string1;                         //Is ignored

  int n_extensions = 7;                    //Keep this as the length of the below list
  char* output_extensions[7] = {"cl.dat","pk.dat","tk.dat","parameters.ini","background.dat","thermodynamics.dat","perturbations_k0.dat"};

  /* Shorthand notation */
  struct file_content * pfc = *ppfc_input;


  /** Check whether a root name has been set, and wether overwrite_root is true */
  class_call(parser_read_string(pfc,"root",&string1,&flag1,errmsg),
             errmsg, errmsg);

  /* Set overwrite_root */
  overwrite_root = _FALSE_;
  class_read_flag("overwrite_root",overwrite_root);

  /** If root has not been set, use the default of 'output/<this-filename>' */
  if (flag1 == _FALSE_){
    memcpy(outfname, "output/", 7);
    memcpy(outfname+7, input_file, strlen(input_file)-4);
    outfname[7+strlen(input_file)-4] = '\0';
  }
  /* Check here for the index of the 'root' field in case it was set in fc_input */
  else{
    for(index_root_in_fc_input=0;index_root_in_fc_input<pfc->size;++index_root_in_fc_input){
      if(strcmp(pfc->name[index_root_in_fc_input],"root") == 0){
        strcpy(outfname,pfc->value[index_root_in_fc_input]);
        break;
      }
    }
  }

  /** If we don't want to overwrite the root name, check now for the existence of output for the given root name + N */
  if(overwrite_root == _FALSE_){

    /* Assume files exist, until proven otherwise */
    found_filenum = _TRUE_;
    /** For each 'filenum', test if it exists. Only stop if it has not been found. */
    for (filenum = 0; filenum < _N_FILEROOT_ && found_filenum; filenum++){
      /* No file has been found yet */
      found_filenum = _FALSE_;
      for(iextens = 0; iextens < n_extensions; ++iextens){
        sprintf(tmp_file,"%s%02d_%s", outfname, filenum, output_extensions[iextens]);
        if (file_exists(tmp_file) == _TRUE_){
          /* Found a file, the outer loop is forced to keep searching */
          found_filenum = _TRUE_;
        }
      }
      /* Didn't find a file. This is the correct number. Break the loop. */
      if(found_filenum == _FALSE_){
        break;
      }
    }
    /* If no root was found, add root through the parser routine */
    if(flag1 == _FALSE_){
      class_call(parser_init(&fc_root,
                             1,
                             pfc->filename,
                             errmsg),
                 errmsg,errmsg);
      sprintf(fc_root.name[0],"root");
      sprintf(fc_root.value[0],"%s%02d_",outfname,filenum);
      fc_root.read[0] = _FALSE_;
      class_call(parser_cat(pfc,
                            &fc_root,
                            pfc_setroot,
                            errmsg),
                 errmsg,
                 errmsg);
      class_call(parser_free(pfc),
                 errmsg,
                 errmsg);
      class_call(parser_free(&fc_root),
                 errmsg,
                 errmsg);
      (*ppfc_input) = pfc_setroot;
    }
    /* If root was found, set the index in the fc_input struct */
    else{
      sprintf(pfc->value[index_root_in_fc_input],"%s%02d_",outfname,filenum);
      (*ppfc_input) = pfc;
    }
  }

  /** If we do want to overwrite, just take the given root name */
  else{
    /* If no root was found, add root through the parser routine */
    if(flag1 == _FALSE_){
      class_call(parser_init(&fc_root,
                             1,
                             pfc->filename,
                             errmsg),
                 errmsg,errmsg);
      sprintf(fc_root.name[0],"root");
      sprintf(fc_root.value[0],"%s_",outfname);
      fc_root.read[0] = _FALSE_;
      class_call(parser_cat(pfc,
                            &fc_root,
                            pfc_setroot,
                            errmsg),
                 errmsg,
                 errmsg);
      class_call(parser_free(pfc),
                 errmsg,
                 errmsg);
      class_call(parser_free(&fc_root),
                 errmsg,
                 errmsg);
      (*ppfc_input) = pfc_setroot;
    }
    /* If root was found, set the index in the fc_input struct */
    else{
      sprintf(pfc->value[index_root_in_fc_input],"%s_",outfname);
      (*ppfc_input) = pfc;
    }
  }

  return _SUCCESS_;
}



/**
 * Initialize each parameter, first to its default values, and then
 * from what can be interpreted from the values passed in the input
 * 'file_content' structure. If its size is null, all parameters keep
 * their default values.
 *
 * @param pfc     Input: pointer to local structure
 * @param ppr     Input: pointer to precision structure
 * @param pba     Input: pointer to background structure
 * @param pth     Input: pointer to thermodynamics structure
 * @param ppt     Input: pointer to perturbation structure
 * @param ptr     Input: pointer to transfer structure
 * @param ppm     Input: pointer to primordial structure
 * @param phr     Input: pointer to harmonic structure
 * @param pfo     Input: pointer to fourier structure
 * @param ple     Input: pointer to lensing structure
 * @param psd     Input: pointer to distorsion structure
 * @param pop     Input: pointer to output structure
 * @param errmsg  Input/Output: Error message
 * @return the error status
 */

int input_read_from_file(struct file_content * pfc,
                         struct precision * ppr,
                         struct background *pba,
                         struct thermodynamics *pth,
                         struct perturbations *ppt,
                         struct transfer *ptr,
                         struct primordial *ppm,
                         struct harmonic *phr,
                         struct fourier * pfo,
                         struct lensing *ple,
                         struct distortions *psd,
                         struct output *pop,
                         ErrorMsg errmsg) {


  /** Summary: */

  /** - Define local variables */
  int input_verbose = 0;
  int has_shooting;
  int flag = _FALSE_;
  int i;

  /** Set default values
      Before getting into the assignment of parameters and the shooting, we want
      to already fix our precision parameters. No precision parameter should
      depend on any input parameter  */
  class_call(input_read_precisions(pfc,ppr,pba,pth,ppt,ptr,ppm,phr,pfo,ple,psd,pop,
                                   errmsg),
             errmsg,
             errmsg);

  class_read_int("input_verbose",input_verbose);
  if (input_verbose >0) printf("Reading input parameters\n");

  /** Find out if shooting necessary and, eventually, shoot and initialize
      read parameters */
  class_call(input_shooting(pfc,ppr,pba,pth,ppt,ptr,ppm,phr,pfo,ple,psd,pop,
                            input_verbose,
                            &has_shooting,
                            errmsg),
             errmsg,
             errmsg);

  /** If no shooting is necessary, initialize read parameters without it */
  if(has_shooting == _FALSE_){
    class_call(input_read_parameters(pfc,ppr,pba,pth,ppt,ptr,ppm,phr,pfo,ple,psd,pop,
                                     errmsg),
               errmsg,
               errmsg);
  }

  /** Read the 'write_warnings' flag. This is the correct place to do it,
      since we want it to happen after all the shooting business */
  class_read_flag_or_deprecated("write_warnings","write warnings",flag);
  /* Print warnings for unread parameters */
  if (flag == _TRUE_){
    for (i=0; i<pfc->size; i++) {
      if (pfc->read[i] == _FALSE_)
        fprintf(stdout,"[WARNING: input line not used: '%s=%s']\n",pfc->name[i],pfc->value[i]);
    }
  }

  if (pfo->has_pk_eq == _TRUE_) {

    if (input_verbose > 0) {
      printf(" -> since you want to use Halofit with a non-zero wa_fld and the Pk_equal method,\n");
      printf("    calling background module to extract the effective w(tau), Omega_m(tau) parameters");
      printf("    required by this method\n");
    }
    class_call(input_prepare_pk_eq(ppr,pba,pth,pfo,input_verbose,errmsg),
               errmsg,
               errmsg);
  }

  return _SUCCESS_;

}


/**
 * In CLASS, we call 'shooting' the process of doing preliminary runs
 * of parts of the code in order to find numerically the value of an
 * input variable which cannot be inferred analytically from other
 * input variables passed by the user.
 *
 * A typical example is when the user passes theta_s, the angular
 * scale of the sound horizon at decoupling. This quantity be passed
 * instead of the hubble parameter h, but only if we run CLASS until
 * the thermodynamics module to figure out how h and theta_s relate
 * numerically. The code starts from a guess for h, and runs to find
 * the corresponding theta_s. It adjusts h, shoots again, and repeats
 * this process until it finds some h giving the correct theta_s
 * within some tolerance.
 *
 * This function contains the overall structure to handle these steps.
 *
 * @param pfc               Input/Output: pointer to file content, with input parameters before/after the shooting
 * @param ppr               Input: pointer to precision structure
 * @param pba               Input: pointer to background structure
 * @param pth               Input: pointer to thermodynamics structure
 * @param ppt               Input: pointer to perturbation structure
 * @param ptr               Input: pointer to transfer structure
 * @param ppm               Input: pointer to primordial structure
 * @param phr               Input: pointer to harmonic structure
 * @param pfo               Input: pointer to fourier structure
 * @param ple               Input: pointer to lensing structure
 * @param psd               Input: pointer to distorsion structure
 * @param pop               Input: pointer to output structure
 * @param input_verbose     Input: Verbosity of input
 * @param has_shooting      Output: do we need shooting?
 * @param errmsg            Input/Output: Error message
 * @return the error status
 */

int input_shooting(struct file_content * pfc,
                   struct precision * ppr,
                   struct background *pba,
                   struct thermodynamics *pth,
                   struct perturbations *ppt,
                   struct transfer *ptr,
                   struct primordial *ppm,
                   struct harmonic *phr,
                   struct fourier * pfo,
                   struct lensing *ple,
                   struct distortions *psd,
                   struct output *pop,
                   int input_verbose,
                   int * has_shooting,
                   ErrorMsg errmsg){

  /** Summary: */

  /** Define local variables */
  int flag1;
  double param1;
  double * unknown_parameter;
  int unknown_parameters_size;
  int counter, index_target, i;
  int fevals=0;
  double xzero;
  double *dxdF, *x_inout;
  int target_indices[_NUM_TARGETS_];
  int needs_shooting;
  int shooting_failed=_FALSE_;

  /* array of parameters passed by the user for which we need shooting (= target parameters) */
  char * const target_namestrings[] = {"100*theta_s",
                                       "Omega_dcdmdr",
                                       "omega_dcdmdr",
                                       "Omega_scf",
                                       "Omega_ini_dcdm",
                                       "omega_ini_dcdm",
                                       "sigma8"};

  /* array of corresponding parameters that must be adjusted in order to meet the target (= unknown parameters) */
  char * const unknown_namestrings[] = {"h",                        /* unknown param for target '100*theta_s' */
                                        "Omega_ini_dcdm",           /* unknown param for target 'Omega_dcdmd' */
                                        "Omega_ini_dcdm",           /* unknown param for target 'omega_dcdmdr"' */
                                        "scf_shooting_parameter",   /* unknown param for target 'Omega_scf' */
                                        "Omega_dcdmdr",             /* unknown param for target 'Omega_ini_dcdm' */
                                        "omega_dcdmdr",             /* unknown param for target 'omega_ini_dcdm' */
                                        "A_s"};                     /* unknown param for target 'sigma8' */

  /* for each target, module up to which we need to run CLASS in order
     to compute the targetted quantities (not running the whole code
     each time to saves a lot of time) */
  enum computation_stage target_cs[] = {cs_thermodynamics, /* computation stage for target '100*theta_s' */
                                        cs_background,     /* computation stage for target 'Omega_dcdmdr' */
                                        cs_background,     /* computation stage for target 'omega_dcdmdr' */
                                        cs_background,     /* computation stage for target 'Omega_scf' */
                                        cs_background,     /* computation stage for target 'Omega_ini_dcdm' */
                                        cs_background,     /* computation stage for target 'omega_ini_dcdm' */
                                        cs_nonlinear};       /* computation stage for target 'sigma8' */

  struct fzerofun_workspace fzw;

  *has_shooting=_FALSE_;

  /** Do we need to fix unknown parameters? */
  unknown_parameters_size = 0;
  fzw.required_computation_stage = 0;
  for (index_target = 0; index_target < _NUM_TARGETS_; index_target++){
    class_call(parser_read_double(pfc,target_namestrings[index_target],&param1,&flag1,errmsg),
               errmsg,
               errmsg);
    if (flag1 == _TRUE_){
      /* input_needs_shoting_for_target takes care of the case where, for
         instance, Omega_dcdmdr is set to 0.0, and we don't need shooting */
      class_call(input_needs_shooting_for_target(pfc,
                                                index_target,
                                                param1,
                                                &needs_shooting,
                                                errmsg),
                 errmsg,
                 errmsg);

      if (needs_shooting == _TRUE_){
        target_indices[unknown_parameters_size] = index_target;
        fzw.required_computation_stage = MAX(fzw.required_computation_stage,target_cs[index_target]);
        unknown_parameters_size++;
      }

    }
  }

  /** In the case of unknown parameters, start shooting... */
  if (unknown_parameters_size > 0) {

    /* We need to remember that we shot so we can clean up properly */
    *has_shooting=_TRUE_;

    /* Create file content structure with additional entries */
    class_call(parser_init(&(fzw.fc),
                           pfc->size+unknown_parameters_size,
                           pfc->filename,
                           errmsg),
               errmsg,errmsg);

    /* Copy input file content to the new file content structure: */
    memcpy(fzw.fc.name, pfc->name, pfc->size*sizeof(FileArg));
    memcpy(fzw.fc.value, pfc->value, pfc->size*sizeof(FileArg));
    memcpy(fzw.fc.read, pfc->read, pfc->size*sizeof(short));

    class_alloc(unknown_parameter,
                unknown_parameters_size*sizeof(double),
                errmsg);
    class_alloc(fzw.unknown_parameters_index,
                unknown_parameters_size*sizeof(int),
                errmsg);

    fzw.target_size = unknown_parameters_size;
    class_alloc(fzw.target_name,
                fzw.target_size*sizeof(enum target_names),
                errmsg);
    class_alloc(fzw.target_value,
                fzw.target_size*sizeof(double),
                errmsg);

    /** Go through all cases with unknown parameters */
    for (counter = 0; counter < unknown_parameters_size; counter++){
      index_target = target_indices[counter];
      class_call(parser_read_double(pfc,
                                    target_namestrings[index_target],
                                    &param1,
                                    &flag1,
                                    errmsg),
               errmsg,
               errmsg);

      /* store name of target parameter */
      fzw.target_name[counter] = index_target;
      /* store target value of target parameter */
      fzw.target_value[counter] = param1;
      fzw.unknown_parameters_index[counter]=pfc->size+counter;
      /* substitute the name of the target parameter with the name of the
         corresponding unknown parameter */
      strcpy(fzw.fc.name[fzw.unknown_parameters_index[counter]],unknown_namestrings[index_target]);
    }

    /** If there is only one parameter, we use a more efficient Newton method for 1D cases */
    if (unknown_parameters_size == 1){

      /* We can do 1 dimensional root finding */
      if (input_verbose > 0) {
        fprintf(stdout,
                "Computing unknown input parameter '%s' using input parameter '%s'\n",
                fzw.fc.name[fzw.unknown_parameters_index[0]],
                target_namestrings[fzw.target_name[0]]);
      }

      /* If shooting fails, postpone error to background module to play nice with MontePython. */
      class_call_try(input_find_root(&xzero,
                                     &fevals,
                                     ppr->tol_shooting_deltax_rel,
                                     &fzw,
                                     errmsg),
                     errmsg,
                     pba->shooting_error,
                     shooting_failed=_TRUE_);

      /* Store xzero */
      // This needs to be done with enough accuracy. A standard double has a relative
      // precision of around 1e-16, so 1e-20 should be good enough for the shooting
      sprintf(fzw.fc.value[fzw.unknown_parameters_index[0]],"%.20e",xzero);
      if (input_verbose > 0) {
        fprintf(stdout," -> found '%s = %s'\n",
                fzw.fc.name[fzw.unknown_parameters_index[0]],
                fzw.fc.value[fzw.unknown_parameters_index[0]]);
      }

    }
    /** Otherwise we do multidimensional shooting */
    else{

      /* We need to do multidimensional root finding */
      if (input_verbose > 0) {
        fprintf(stdout,"Computing unknown input parameters\n");
      }

      /* Allocate local variables */
      class_alloc(x_inout,
                  sizeof(double)*unknown_parameters_size,
                  errmsg);
      class_alloc(dxdF,
                  sizeof(double)*unknown_parameters_size,
                  errmsg);

      /* Get the guess for the initial variables */
      class_call(input_get_guess(x_inout, dxdF, &fzw, errmsg),
                 errmsg,
                 errmsg);

      /* Use multi-dimensional Newton method */
      class_call_try(fzero_Newton(input_try_unknown_parameters,
                                  x_inout,
                                  dxdF,
                                  unknown_parameters_size,
                                  ppr->tol_shooting_deltax,
                                  ppr->tol_shooting_deltaF,
                                  &fzw,
                                  &fevals,
                                  errmsg),
                     errmsg,
                     pba->shooting_error,
                     shooting_failed=_TRUE_);

      /* Store xzero */
      // This needs to be done with enough accuracy. A standard double has a relative
      // precision of around 1e-16, so 1e-20 should be good enough for the shooting
      for (counter = 0; counter < unknown_parameters_size; counter++){
        sprintf(fzw.fc.value[fzw.unknown_parameters_index[counter]],
                "%.20e",x_inout[counter]);
        if (input_verbose > 0) {
          fprintf(stdout," -> found '%s = %s'\n",
                  fzw.fc.name[fzw.unknown_parameters_index[counter]],
                  fzw.fc.value[fzw.unknown_parameters_index[counter]]);
        }
      }

      /* Free local variables */
      free(x_inout);
      free(dxdF);
    }

    if (input_verbose > 1) {
      fprintf(stdout,"Shooting completed using %d function evaluations\n",fevals);
    }

    /** Read all parameters from the fc obtained through shooting */
    class_call(input_read_parameters(&(fzw.fc),ppr,pba,pth,ppt,ptr,ppm,phr,pfo,ple,psd,pop,
                                     errmsg),
               errmsg,
               errmsg);

    /** Set status of shooting */
    pba->shooting_failed = shooting_failed;

    /* all parameters read in fzw must be considered as read in pfc. At the same
       time the parameters read before in pfc (like theta_s,...) must still be
       considered as read (hence we could not do a memcopy) */
    for (i=0; i < pfc->size; i ++) {
      if (fzw.fc.read[i] == _TRUE_)
        pfc->read[i] = _TRUE_;
    }

    /* Free tuned pfc */
    parser_free(&(fzw.fc));

    /** Free arrays allocated */
    free(unknown_parameter);
    free(fzw.unknown_parameters_index);
    free(fzw.target_name);
    free(fzw.target_value);
  }

  return _SUCCESS_;

}


/**
 * Related to 'shooting': for each target, check whether it is
 * sufficient to stick to the default value of the unkown parameter
 * (for instance: if the target parameter is a density and the target
 * value is zero, the unkown parameter should remain zero like in the
 * default)
 *
 * @param pfc             Input: pointer to local structure
 * @param target_name     Input: list of possible target names
 * @param target_value    Input: list of possible target values
 * @param needs_shooting  Output: needs shooting?
 * @param errmsg          Input/Output: Error message
 * @return the error status
 */

int input_needs_shooting_for_target(struct file_content * pfc,
                                    enum target_names target_name,
                                    double target_value,
                                    int * needs_shooting,
                                    ErrorMsg errmsg){

  *needs_shooting = _TRUE_;
  switch (target_name){
    case Omega_dcdmdr:
    case omega_dcdmdr:
    case Omega_scf:
    case Omega_ini_dcdm:
    case omega_ini_dcdm:
      /* Check that Omega's or omega's are nonzero: */
      if (target_value == 0.)
        *needs_shooting = _FALSE_;
      break;
    default:
      /* Default is no additional checks */
      *needs_shooting = _TRUE_;
      break;
  }

  return _SUCCESS_;

}

/**
 * Related to 'shooting': Find the root of a one-dimensional
 * function. This function starts from a first guess, then uses a few
 * steps to bracket the root, and then calls another function to
 * actually get the root.
 *
 * @param xzero     Output: root x such that f(x)=0 up to tolerance (f(x) = input_fzerofun_1d)
 * @param fevals    Output: number of iterations (that is, of CLASS runs) needed to find the root
 * @param tol_x_rel Input : Relative tolerance compared to bracket of root that is used to find root.
 * @param pfzw      Input : pointer to workspace containing targets, unkown parameters and other relevant information
 * @param errmsg    Input/Output: Error message
 * @return the error status
*/

int input_find_root(double *xzero,
                    int *fevals,
                    double tol_x_rel,
                    struct fzerofun_workspace *pfzw,
                    ErrorMsg errmsg){

  /** Summary: */

  /** Define local variables */
  double x1, x2, f1, f2, dxdy, dx;
  int iter, iter2;
  int return_function;

  /** Fisrt we do our guess */
  class_call(input_get_guess(&x1, &dxdy, pfzw, errmsg),
             errmsg,
             errmsg);

  class_call(input_fzerofun_1d(x1, pfzw, &f1, errmsg),
             errmsg,
             errmsg);

  (*fevals)++;
  dx = 1.5*f1*dxdy;

  /** Then we do a linear hunt for the boundaries */
  /* Try fifteen times to go above and below the root (i.e. where shooting succeeds) */
  for (iter=1; iter<=15; iter++){
    x2 = x1 - dx;
    /* Try three times to get a 'reasonable' value, i.e. no CLASS error */
    for (iter2=1; iter2 <= 3; iter2++) {
      return_function = input_fzerofun_1d(x2, pfzw, &f2, errmsg);
      (*fevals)++;
      if (return_function ==_SUCCESS_) {
        break;
      }
      else if (iter2 < 3) {
        dx*=0.5;
        x2 = x1-dx;
      }
      else {
        class_stop(errmsg,errmsg);
      }
    }
    if (f1*f2<0.0){
      /* Root has been bracketed */
      break;
    }
    x1 = x2;
    f1 = f2;
  }

  /** Find root using Ridders method (Exchange for bisection if you are old-school) */
  class_call(input_fzero_ridder(input_fzerofun_1d,
                                x1,
                                x2,
                                tol_x_rel*MAX(fabs(x1),fabs(x2)),
                                pfzw,
                                &f1,
                                &f2,
                                xzero,
                                fevals,
                                errmsg),
             errmsg,errmsg);

  return _SUCCESS_;

}

/**
 * Related to 'shooting': defines 1d function of which we want to find
 * the root during the shooting.  The function is simply: "prediction
 * of CLASS for a target parameter y given a parameter x - targeted
 * value of y"
 *
 * @param input         Input: value of x
 * @param pfzw          Input: pointer to workspace containing targets, unkown parameters and other relevant information
 * @param output        Ouput: f(x) = y - y_targeted
 * @param error_message Input/Output: Error message
 * @return the error status
*/

int input_fzerofun_1d(double input,
                      void* pfzw,
                      double *output,
                      ErrorMsg error_message){

  class_call(input_try_unknown_parameters(&input,
                                          1,
                                          pfzw,
                                          output,
                                          error_message),
             error_message,
             error_message);

  return _SUCCESS_;

}


/**
 * Related to 'shooting': using Ridders' method, return the root x of
 * a function f(x) known to lie between x1 and x2, up to some
 * tolerance. Note that this function is very generic and could easily
 * be moved to the tools (and be used in other modules).
 *
 * @param func   Input: function y=f(x), with arguments: x, pointer to y, and another pointer containing several fixed parameters
 * @param x1            Input: lower boundary x1<x
 * @param x2            Input: upper boundary x<x2
 * @param xtol          Input: tolerance: |x- true root|<xtol
 * @param param         Input: fixed parameters passed to f(x)
 * @param Fx1           Input: f(x1)
 * @param Fx2           Input: f(x2)
 * @param xzero         Output: root x
 * @param fevals        Output: number of iterations (that is, of CLASS runs) needed to find the root
 * @param error_message Input/Output: Error message
 * @return the error status
 */

int input_fzero_ridder(int (*func)(double x,
                                   void *param,
                                   double *y,
                                   ErrorMsg error_message),
                       double x1,
                       double x2,
                       double xtol,
                       void *param,
                       double *Fx1,
                       double *Fx2,
                       double *xzero,
                       int *fevals,
                       ErrorMsg error_message){

  /** Summary: */

  /** Define local variables */
  int j,MAXIT=1000;
  double ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew;

  if ((Fx1!=NULL)&&(Fx2!=NULL)){
    fl = *Fx1;
    fh = *Fx2;
  }
  else{
    class_call((*func)(x1, param, &fl, error_message),
               error_message, error_message);
    class_call((*func)(x2, param, &fh, error_message),
               error_message, error_message);

    *fevals = (*fevals)+2;
  }
  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
    xl=x1;
    xh=x2;
    ans=-1.11e11;
    for (j=1;j<=MAXIT;j++) {
      xm=0.5*(xl+xh);
      class_call((*func)(xm, param, &fm, error_message),
                 error_message, error_message);
      *fevals = (*fevals)+1;
      s=sqrt(fm*fm-fl*fh);
      if (s == 0.0){
        *xzero = ans;
        return _SUCCESS_;
      }
      xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
      if (fabs(xnew-ans) <= xtol) {
        *xzero = ans;
        return _SUCCESS_;
      }

      ans=xnew;
      class_call((*func)(ans, param, &fnew, error_message),
                 error_message, error_message);
      *fevals = (*fevals)+1;
      if (fnew == 0.0){
        *xzero = ans;
        return _SUCCESS_;
      }

      if (NRSIGN(fm,fnew) != fm) {
        xl=xm;
        fl=fm;
        xh=ans;
        fh=fnew;
      }
      else if (NRSIGN(fl,fnew) != fl) {
        xh=ans;
        fh=fnew;
      }
      else if (NRSIGN(fh,fnew) != fh) {
        xl=ans;
        fl=fnew;
      }
      else{
        return _FAILURE_;
      }
      if (fabs(xh-xl) <= xtol) {
        *xzero = ans;
        return _SUCCESS_;
      }
    }
    class_stop(error_message,"zriddr exceed maximum iterations");
  }
  else {
    if (fl == 0.0) return x1;
    if (fh == 0.0) return x2;
    class_stop(error_message,"root must be bracketed in zriddr.");
  }
  class_stop(error_message,"Failure in int.");
}


/**
 * Related to 'shooting': we define here a reasonable analytic guess
 * for each unknown parameter as a function of its target
 * parameter. We must also estimate dxdy, i.e. how the unknown
 * parameter responds to the target parameter.  This can simply be
 * estimated as the derivative of the guess formula.
 *
 * @param xguess Output: guess for unkown parameter x given target parameter y
 * @param dxdy   Output: guess for derivative dx/dy
 * @param pfzw   Input : pointer to workspace containing targets, unkown parameters and other relevant information
 * @param errmsg Input/Output: Error message
 * @return the error status
 */

int input_get_guess(double *xguess,
                    double *dxdy,
                    struct fzerofun_workspace * pfzw,
                    ErrorMsg errmsg){

  /** Summary: */

  /** Define local variables */
  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermodynamics th;           /* for thermodynamics */
  struct perturbations pt;         /* for source functions */
  struct transfer tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct harmonic hr;          /* for output spectra */
  struct fourier fo;        /* for non-linear spectra */
  struct lensing le;          /* for lensed spectra */
  struct distortions sd;      /* for spectral distortions */
  struct output op;           /* for output files */
  int i;
  double Omega_M, a_decay, gamma, Omega0_dcdmdr=1.0;
  int index_guess;

  /* Cheat to read only known parameters: */
  pfzw->fc.size -= pfzw->target_size;

  class_call(input_read_precisions(&(pfzw->fc),&pr,&ba,&th,&pt,&tr,&pm,&hr,&fo,&le,&sd,&op,
                                   errmsg),
             errmsg,
             errmsg);
  class_call(input_read_parameters(&(pfzw->fc),&pr,&ba,&th,&pt,&tr,&pm,&hr,&fo,&le,&sd,&op,
                                   errmsg),
             errmsg,
             errmsg);

  pfzw->fc.size += pfzw->target_size;

  /** Estimate dxdy */
  for (index_guess=0; index_guess < pfzw->target_size; index_guess++) {
    switch (pfzw->target_name[index_guess]) {
    case theta_s:
      xguess[index_guess] = 3.54*pow(pfzw->target_value[index_guess],2)-5.455*pfzw->target_value[index_guess]+2.548;
      dxdy[index_guess] = (7.08*pfzw->target_value[index_guess]-5.455);
      /** Update pb to reflect guess */
      ba.h = xguess[index_guess];
      ba.H0 = ba.h *  1.e5 / _c_;
      break;
    case Omega_dcdmdr:
      Omega_M = ba.Omega0_cdm+ba.Omega0_idm_dr+ba.Omega0_dcdmdr+ba.Omega0_b;
      /* *
       * This formula is exact in a Matter + Lambda Universe, but only for Omega_dcdm,
       * not the combined.
       * sqrt_one_minus_M = sqrt(1.0 - Omega_M);
       * xguess[index_guess] = pfzw->target_value[index_guess]*
       *                       exp(2./3.*ba.Gamma_dcdm/ba.H0*
       *                       atanh(sqrt_one_minus_M)/sqrt_one_minus_M);
       * dxdy[index_guess] = 1.0;//exp(2./3.*ba.Gamma_dcdm/ba.H0*atanh(sqrt_one_minus_M)/sqrt_one_minus_M);
       * */
      gamma = ba.Gamma_dcdm/ba.H0;
      if (gamma < 1)
        a_decay = 1.0;
      else
        a_decay = pow(1+(gamma*gamma-1.)/Omega_M,-1./3.);
      xguess[index_guess] = pfzw->target_value[index_guess]/a_decay;
      dxdy[index_guess] = 1./a_decay;
      break;
    case omega_dcdmdr:
      Omega_M = ba.Omega0_cdm+ba.Omega0_idm_dr+ba.Omega0_dcdmdr+ba.Omega0_b;
      gamma = ba.Gamma_dcdm/ba.H0;
      if (gamma < 1)
        a_decay = 1.0;
      else
        a_decay = pow(1+(gamma*gamma-1.)/Omega_M,-1./3.);
      xguess[index_guess] = pfzw->target_value[index_guess]/ba.h/ba.h/a_decay;
      dxdy[index_guess] = 1./a_decay/ba.h/ba.h;
      break;
    case Omega_scf:
      /* *
       * This guess is arbitrary, something nice using WKB should be implemented.
       * Version 2 uses a fit
       * xguess[index_guess] = 1.77835*pow(ba.Omega0_scf,-2./7.);
       * dxdy[index_guess] = -0.5081*pow(ba.Omega0_scf,-9./7.)`;
       * Version 3: use attractor solution
       * */
      if (ba.scf_tuning_index == 0){
        xguess[index_guess] = sqrt(3.0/ba.Omega0_scf);
        dxdy[index_guess] = -0.5*sqrt(3.0)*pow(ba.Omega0_scf,-1.5);
      }
      else{
        /* Default: take the passed value as xguess and set dxdy to 1. */
        xguess[index_guess] = ba.scf_parameters[ba.scf_tuning_index];
        dxdy[index_guess] = 1.;
      }
      break;
    case omega_ini_dcdm:
      Omega0_dcdmdr = 1./(ba.h*ba.h);
    case Omega_ini_dcdm:
      /* This works since correspondence is Omega_ini_dcdm -> Omega_dcdmdr and
         omega_ini_dcdm -> omega_dcdmdr */
      Omega0_dcdmdr *=pfzw->target_value[index_guess];
      Omega_M = ba.Omega0_cdm+ba.Omega0_idm_dr+Omega0_dcdmdr+ba.Omega0_b;
      gamma = ba.Gamma_dcdm/ba.H0;
      if (gamma < 1)
        a_decay = 1.0;
      else
        a_decay = pow(1+(gamma*gamma-1.)/Omega_M,-1./3.);
      xguess[index_guess] = pfzw->target_value[index_guess]*a_decay;
      dxdy[index_guess] = a_decay;
      if (gamma > 100)
        dxdy[index_guess] *= gamma/100;
      break;

    case sigma8:
      /* Assume linear relationship between A_s and sigma8 and fix coefficient
         according to vanilla LambdaCDM. Should be good enough... */
      xguess[index_guess] = 2.43e-9/0.87659*pfzw->target_value[index_guess];
      dxdy[index_guess] = 2.43e-9/0.87659;
      break;
    }
  }

  for (i=0; i<pfzw->fc.size; i++) {
    pfzw->fc.read[i] = _FALSE_;
  }

  /** - Deallocate everything allocated by input_read_parameters */
  background_free_input(&ba);

  return _SUCCESS_;

}

/**
 * Related to 'shooting': when there is one or more targets, call
 * CLASS up to the highest needed computation stage, for a given set
 * of unknown parameters; obtain the corresponding target parameters;
 * and return the vector of each [target - targeted_value].
 *
 * @param unknown_parameter       Input: vector of unkownn parameters x
 * @param unknown_parameters_size Input: size of this vector
 * @param voidpfzw                Input: pointer to workspace containing targets, unkown parameters and other relevant information
 * @param output                  Output: vector of target parameters y
 * @param errmsg                  Input/Output: Error message
 * @return the error status
*/

int input_try_unknown_parameters(double * unknown_parameter,
                                 int unknown_parameters_size,
                                 void * voidpfzw,
                                 double * output,
                                 ErrorMsg errmsg){
  /** Summary */

  /** Define local variables */
  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermodynamics th;           /* for thermodynamics */
  struct perturbations pt;         /* for source functions */
  struct transfer tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct harmonic hr;          /* for output spectra */
  struct fourier fo;        /* for non-linear spectra */
  struct lensing le;          /* for lensed spectra */
  struct distortions sd;      /* for spectral distortions */
  struct output op;           /* for output files */

  int i;
  double rho_dcdm_today, rho_dr_today;
  struct fzerofun_workspace * pfzw;
  int input_verbose;
  int flag;
  int param;
  short compute_sigma8 = _FALSE_;

  pfzw = (struct fzerofun_workspace *) voidpfzw;
  /** Read input parameters */
  // This needs to be done with enough accuracy. A standard double has a relative
  // precision of around 1e-16, so 1e-20 should be good enough for the shooting
  for (i=0; i < unknown_parameters_size; i++) {
    sprintf(pfzw->fc.value[pfzw->unknown_parameters_index[i]],"%.20e",unknown_parameter[i]);
  }

  class_call(input_read_precisions(&(pfzw->fc),&pr,&ba,&th,&pt,&tr,&pm,&hr,&fo,&le,&sd,&op,
                                   errmsg),
             errmsg,
             errmsg);

  class_call(input_read_parameters(&(pfzw->fc),&pr,&ba,&th,&pt,&tr,&pm,&hr,&fo,&le,&sd,&op,
                                   errmsg),
             errmsg,
             errmsg);

  class_call(parser_read_int(&(pfzw->fc),"input_verbose",&param,&flag,errmsg),
             errmsg,
             errmsg);

  if (flag == _TRUE_)
    input_verbose = param;
  else
    input_verbose = 0;

  /** Optimise flags for sigma8 calculation.*/
  for (i=0; i < unknown_parameters_size; i++) {
    if (pfzw->target_name[i] == sigma8) {
      compute_sigma8 = _TRUE_;
    }
  }
  /* Sigma8 depends on linear P(k), so no need to run anything except linear P(k) during shooting */
  if (compute_sigma8 == _TRUE_) {
    pt.k_max_for_pk=10.0; // increased in June 2020 for higher accuracy
    pt.has_pk_matter=_TRUE_;
    pt.has_perturbations = _TRUE_;
    pt.has_cl_cmb_temperature = _FALSE_;
    pt.has_cls = _FALSE_;
    pt.has_cl_cmb_polarization = _FALSE_;
    pt.has_cl_cmb_lensing_potential = _FALSE_;
    pt.has_cl_number_count = _FALSE_;
    pt.has_cl_lensing_potential=_FALSE_;
    pt.has_density_transfers=_FALSE_;
    pt.has_velocity_transfers=_FALSE_;
    fo.has_pk_eq=_FALSE_;
    fo.method=nl_none;
  }

  /** Shoot forward into class up to required stage */
  if (pfzw->required_computation_stage >= cs_background){
    if (input_verbose>2)
      printf("Stage 1: background\n");
    ba.background_verbose = 0;
    class_call(background_init(&pr,&ba), ba.error_message, errmsg);
  }

  if (pfzw->required_computation_stage >= cs_thermodynamics){
   if (input_verbose>2)
     printf("Stage 2: thermodynamics\n");
    pr.thermo_Nz_lin = 10000;
    pr.thermo_Nz_log = 500;
    th.thermodynamics_verbose = 0;
    th.hyrec_verbose = 0;
    class_call_except(thermodynamics_init(&pr,&ba,&th), th.error_message, errmsg, background_free(&ba));
  }

  if (pfzw->required_computation_stage >= cs_perturbations){
       if (input_verbose>2)
         printf("Stage 3: perturbations\n");
    pt.perturbations_verbose = 0;
    class_call_except(perturbations_init(&pr,&ba,&th,&pt), pt.error_message, errmsg, thermodynamics_free(&th);background_free(&ba));
  }

  if (pfzw->required_computation_stage >= cs_primordial){
    if (input_verbose>2)
      printf("Stage 4: primordial\n");
    pm.primordial_verbose = 0;
    class_call_except(primordial_init(&pr,&pt,&pm), pm.error_message, errmsg, perturbations_free(&pt);thermodynamics_free(&th);background_free(&ba));
  }

  if (pfzw->required_computation_stage >= cs_nonlinear){
    if (input_verbose>2)
      printf("Stage 5: nonlinear\n");
    fo.fourier_verbose = 0;
    class_call_except(fourier_init(&pr,&ba,&th,&pt,&pm,&fo), fo.error_message, errmsg, primordial_free(&pm);perturbations_free(&pt);thermodynamics_free(&th);background_free(&ba));
  }

  if (pfzw->required_computation_stage >= cs_transfer){
    if (input_verbose>2)
      printf("Stage 6: transfer\n");
    tr.transfer_verbose = 0;
    class_call_except(transfer_init(&pr,&ba,&th,&pt,&fo,&tr), tr.error_message, errmsg, fourier_free(&fo);primordial_free(&pm);perturbations_free(&pt);thermodynamics_free(&th);background_free(&ba));
  }

  if (pfzw->required_computation_stage >= cs_spectra){
    if (input_verbose>2)
      printf("Stage 7: spectra\n");
    hr.harmonic_verbose = 0;
    class_call_except(harmonic_init(&pr,&ba,&pt,&pm,&fo,&tr,&hr),hr.error_message, errmsg, transfer_free(&tr);fourier_free(&fo);primordial_free(&pm);perturbations_free(&pt);thermodynamics_free(&th);background_free(&ba));
  }

  /** Get the corresponding shoot variable and put into output */
  for (i=0; i < pfzw->target_size; i++) {
    switch (pfzw->target_name[i]) {
    case theta_s:
      output[i] = 100.*th.rs_rec/th.ra_rec-pfzw->target_value[i];
      break;
    case Omega_dcdmdr:
      rho_dcdm_today = ba.background_table[(ba.bt_size-1)*ba.bg_size+ba.index_bg_rho_dcdm];
      if (ba.has_dr == _TRUE_)
        rho_dr_today = ba.background_table[(ba.bt_size-1)*ba.bg_size+ba.index_bg_rho_dr];
      else
        rho_dr_today = 0.;
      output[i] = (rho_dcdm_today+rho_dr_today)/(ba.H0*ba.H0)-pfzw->target_value[i];
      break;
    case omega_dcdmdr:
      rho_dcdm_today = ba.background_table[(ba.bt_size-1)*ba.bg_size+ba.index_bg_rho_dcdm];
      if (ba.has_dr == _TRUE_)
        rho_dr_today = ba.background_table[(ba.bt_size-1)*ba.bg_size+ba.index_bg_rho_dr];
      else
        rho_dr_today = 0.;
      output[i] = (rho_dcdm_today+rho_dr_today)/(ba.H0*ba.H0)-pfzw->target_value[i]/ba.h/ba.h;
      break;
    case Omega_scf:
      /** In case scalar field is used to fill, pba->Omega0_scf is not equal to pfzw->target_value[i].*/
      output[i] = ba.background_table[(ba.bt_size-1)*ba.bg_size+ba.index_bg_rho_scf]/(ba.H0*ba.H0)-ba.Omega0_scf;
      break;
    case Omega_ini_dcdm:
    case omega_ini_dcdm:
      rho_dcdm_today = ba.background_table[(ba.bt_size-1)*ba.bg_size+ba.index_bg_rho_dcdm];
      if (ba.has_dr == _TRUE_)
        rho_dr_today = ba.background_table[(ba.bt_size-1)*ba.bg_size+ba.index_bg_rho_dr];
      else
        rho_dr_today = 0.;
      output[i] = -(rho_dcdm_today+rho_dr_today)/(ba.H0*ba.H0)+ba.Omega0_dcdmdr;
      break;
    case sigma8:
      output[i] = fo.sigma8[fo.index_pk_m]-pfzw->target_value[i];
      break;
    }
  }

  /** Free structures */
  if (pfzw->required_computation_stage >= cs_spectra){
    class_call(harmonic_free(&hr), hr.error_message, errmsg);
  }
  if (pfzw->required_computation_stage >= cs_transfer){
    class_call(transfer_free(&tr), tr.error_message, errmsg);
  }
  if (pfzw->required_computation_stage >= cs_nonlinear){
    class_call(fourier_free(&fo), fo.error_message, errmsg);
  }
  if (pfzw->required_computation_stage >= cs_primordial){
    class_call(primordial_free(&pm), pm.error_message, errmsg);
  }
  if (pfzw->required_computation_stage >= cs_perturbations){
    class_call(perturbations_free(&pt), pt.error_message, errmsg);
  }
  if (pfzw->required_computation_stage >= cs_thermodynamics){
    class_call(thermodynamics_free(&th), th.error_message, errmsg);
  }
  if (pfzw->required_computation_stage >= cs_background){
    class_call(background_free(&ba), ba.error_message, errmsg);
  }

  /** Set filecontent to unread */
  for (i=0; i<pfzw->fc.size; i++) {
    pfzw->fc.read[i] = _FALSE_;
  }

  return _SUCCESS_;

}

/**
 * Initialize the precision parameter structure.
 *
 * All precision parameters used in the other modules are listed here
 * and assigned here a default value.
 *
 * @param pfc     Input: pointer to local structure
 * @param ppr     Input: pointer to precision structure
 * @param pba     Input: pointer to background structure
 * @param pth     Input: pointer to thermodynamics structure
 * @param ppt     Input: pointer to perturbations structure
 * @param ptr     Input: pointer to transfer structure
 * @param ppm     Input: pointer to primordial structure
 * @param phr     Input: pointer to harmonic structure
 * @param pfo     Input: pointer to non-linear structure
 * @param ple     Input: pointer to lensing structure
 * @param pop     Input: pointer to output structure
 * @param psd     Input: pointer to distorsion structure
 * @param errmsg  Input: Error message
 * @return the error status
 */

int input_read_precisions(struct file_content * pfc,
                          struct precision * ppr,
                          struct background *pba,
                          struct thermodynamics *pth,
                          struct perturbations *ppt,
                          struct transfer *ptr,
                          struct primordial *ppm,
                          struct harmonic *phr,
                          struct fourier * pfo,
                          struct lensing *ple,
                          struct distortions *psd,
                          struct output *pop,
                          ErrorMsg errmsg){

  /** Summary: */

  /** - Define local variables */

  /** - Automatic estimate of machine precision */
  ppr->smallest_allowed_variation = DBL_EPSILON;

  class_test(ppr->smallest_allowed_variation < 0,
             ppr->error_message,
             "smallest_allowed_variation = %e < 0",
             ppr->smallest_allowed_variation);

  /* Assign the default precision settings (these very concise lines
     assign all precision parameters thanks to the macros defined in
     macros_precision.h) */
  #define __ASSIGN_DEFAULT_PRECISION__
  #include "precisions.h"
  #undef __ASSIGN_DEFAULT_PRECISION__

  /** Read all precision parameters from input (these very concise
      lines parse all precision parameters thanks to the macros
      defined in macros_precision.h) */
  #define __PARSE_PRECISION_PARAMETER__
  #include "precisions.h"
  #undef __PARSE_PRECISION_PARAMETER__

  return _SUCCESS_;

}


/**
 * If entries are passed in file_content structure, carefully read and
 * interpret each of them, and tune the relevant input parameters
 * accordingly
 *
 * @param pfc     Input: pointer to local structure
 * @param ppr     Input: pointer to precision structure
 * @param pba     Input: pointer to background structure
 * @param pth     Input: pointer to thermodynamics structure
 * @param ppt     Input: pointer to perturbation structure
 * @param ptr     Input: pointer to transfer structure
 * @param ppm     Input: pointer to primordial structure
 * @param phr     Input: pointer to harmonic structure
 * @param pfo     Input: pointer to fourier structure
 * @param ple     Input: pointer to lensing structure
 * @param psd     Input: pointer to distorsion structure
 * @param pop     Input: pointer to output structure
 * @param errmsg  Input: Error message
 * @return the error status
 */

int input_read_parameters(struct file_content * pfc,
                          struct precision * ppr,
                          struct background *pba,
                          struct thermodynamics *pth,
                          struct perturbations *ppt,
                          struct transfer *ptr,
                          struct primordial *ppm,
                          struct harmonic *phr,
                          struct fourier * pfo,
                          struct lensing *ple,
                          struct distortions *psd,
                          struct output *pop,
                          ErrorMsg errmsg){

  /** Summary: */

  /** Define local variables */
  int input_verbose=0;

  /** Set all input parameters to default values */
  class_call(input_default_params(pba,pth,ppt,ptr,ppm,phr,pfo,ple,psd,pop),
             errmsg,
             errmsg);

  /** Read verbose for input structure */
  class_read_int("input_verbose",input_verbose);

  /**
   * Read the general parameters of the
   *  background, thermodynamics, and perturbation structures
   * This function is exclusively for those parameters, NOT
   *  related to any physical species
   * */
  class_call(input_read_parameters_general(pfc,pba,pth,ppt,psd,
                                           errmsg),
             errmsg,
             errmsg);

  /** Read the parameters for each physical species (has to be called after the general read) */
  class_call(input_read_parameters_species(pfc,ppr,pba,pth,ppt,
                                           input_verbose,
                                           errmsg),
             errmsg,
             errmsg);

  /** Read parameters for exotic energy injection quantities */
  class_call(input_read_parameters_injection(pfc,ppr,pth,
                                             errmsg),
             errmsg,
             errmsg);

  /** Read parameters for nonlinear quantities */
  class_call(input_read_parameters_nonlinear(pfc,ppr,pba,pth,ppt,pfo,
                                             input_verbose,
                                             errmsg),
             errmsg,
             errmsg);

  /** Read parameters for primordial quantities */
  class_call(input_read_parameters_primordial(pfc,ppt,ppm,
                                              errmsg),
             errmsg,
             errmsg);

  /** Read parameters for spectra quantities */
  class_call(input_read_parameters_spectra(pfc,ppr,pba,ppm,ppt,ptr,phr,pop,
                                           errmsg),
             errmsg,
             errmsg);

  /** Read parameters for lensing quantities */
  class_call(input_read_parameters_lensing(pfc,ppr,ppt,ptr,ple,
                                           errmsg),
             errmsg,
             errmsg);

  /** Read parameters for distortions quantities */
  class_call(input_read_parameters_distortions(pfc,ppr,psd,
                                               errmsg),
             errmsg,
             errmsg);

  /** Read obsolete parameters */
  class_call(input_read_parameters_additional(pfc,ppr,pba,pth,
                                              errmsg),
             errmsg,
             errmsg);

  /** Read parameters for output quantities */
  class_call(input_read_parameters_output(pfc,pba,pth,ppt,ptr,ppm,phr,pfo,ple,psd,pop,
                                          errmsg),
             errmsg,
             errmsg);

  return _SUCCESS_;

}


/**
 * Read general parameters related to class, including
 *   - background, thermo, and perturbation quantities NOT associated to
 *     any particular species
 *   - calculationary quantities like the gauge/recombination code
 *   - output options
 *
 * @param pfc     Input: pointer to local structure
 * @param pba     Input: pointer to background structure
 * @param pth     Input: pointer to thermodynamics structure
 * @param ppt     Input: pointer to perturbation structure
 * @param psd     Input: pointer to distorsion structure
 * @param errmsg  Input: Error message
 * @return the error status
 */

int input_read_parameters_general(struct file_content * pfc,
                                  struct background * pba,
                                  struct thermodynamics * pth,
                                  struct perturbations * ppt,
                                  struct distortions * psd,
                                  ErrorMsg errmsg){

  /** Summary: */

  /** - Define local variables */
  int flag1,flag2;
  double param1,param2;
  char string1[_ARGUMENT_LENGTH_MAX_];
  char * options_output[33] =  {"tCl","pCl","lCl","nCl","dCl","sCl","mPk","mTk","dTk","vTk","sd",
                                "TCl","PCl","LCl","NCl","DCl","SCl","MPk","MTk","DTk","VTk","Sd",
                                "TCL","PCL","LCL","NCL","DCL","SCL","MPK","MTK","DTK","VTK","SD"};
  char * options_temp_contributions[10] = {"tsw","eisw","lisw","dop","pol","TSW","EISW","LISW","Dop","Pol"};
  char * options_number_count[8] = {"density","dens","rsd","RSD","lensing","lens","gr","GR"};
  char * options_modes[6] = {"s","v","t","S","V","T"};
  char * options_ics[10] = {"ad","bi","cdi","nid","niv","AD","BI","CDI","NID","NIV"};

  /* Set local default values */
  ppt->has_perturbations = _FALSE_;
  ppt->has_cls = _FALSE_;
  psd->has_distortions = _FALSE_;

  /** 1) List of output spectra requested */
  /* Read */
  class_call(parser_read_string(pfc,"output",&string1,&flag1,errmsg),
             errmsg,
             errmsg);
  /* Complete set of parameters */
  if (flag1 == _TRUE_) {
    if ((strstr(string1,"tCl") != NULL) || (strstr(string1,"TCl") != NULL) || (strstr(string1,"TCL") != NULL)) {
      ppt->has_cl_cmb_temperature = _TRUE_;
      ppt->has_perturbations = _TRUE_;
      ppt->has_cls = _TRUE_;
    }
    if ((strstr(string1,"pCl") != NULL) || (strstr(string1,"PCl") != NULL) || (strstr(string1,"PCL") != NULL)) {
      ppt->has_cl_cmb_polarization = _TRUE_;
      ppt->has_perturbations = _TRUE_;
      ppt->has_cls = _TRUE_;
    }
    if ((strstr(string1,"lCl") != NULL) || (strstr(string1,"LCl") != NULL) || (strstr(string1,"LCL") != NULL)) {
      ppt->has_cl_cmb_lensing_potential = _TRUE_;
      ppt->has_perturbations = _TRUE_;
      ppt->has_cls = _TRUE_;
    }
    if ((strstr(string1,"nCl") != NULL) || (strstr(string1,"NCl") != NULL) || (strstr(string1,"NCL") != NULL) ||
        (strstr(string1,"dCl") != NULL) || (strstr(string1,"DCl") != NULL) || (strstr(string1,"DCL") != NULL)) {
      ppt->has_cl_number_count = _TRUE_;
      ppt->has_perturbations = _TRUE_;
      ppt->has_cls = _TRUE_;
    }
    if ((strstr(string1,"sCl") != NULL) || (strstr(string1,"SCl") != NULL) || (strstr(string1,"SCL") != NULL)) {
      ppt->has_cl_lensing_potential=_TRUE_;
      ppt->has_perturbations = _TRUE_;
      ppt->has_cls = _TRUE_;
    }
    if ((strstr(string1,"mPk") != NULL) || (strstr(string1,"MPk") != NULL) || (strstr(string1,"MPK") != NULL)) {
      ppt->has_pk_matter=_TRUE_;
      ppt->has_perturbations = _TRUE_;
    }
    if ((strstr(string1,"mTk") != NULL) || (strstr(string1,"MTk") != NULL) || (strstr(string1,"MTK") != NULL) ||
        (strstr(string1,"dTk") != NULL) || (strstr(string1,"DTk") != NULL) || (strstr(string1,"DTK") != NULL)) {
      ppt->has_density_transfers=_TRUE_;
      ppt->has_perturbations = _TRUE_;
    }
    if ((strstr(string1,"vTk") != NULL) || (strstr(string1,"VTk") != NULL) || (strstr(string1,"VTK") != NULL)) {
      ppt->has_velocity_transfers=_TRUE_;
      ppt->has_perturbations = _TRUE_;
    }
    if ((strstr(string1,"Sd") != NULL) || (strstr(string1,"sd") != NULL) || (strstr(string1,"SD") != NULL)) {
      ppt->has_perturbations = _TRUE_;
      psd->has_distortions=_TRUE_;
      pth->compute_damping_scale=_TRUE_;
    }

    /* The following lines make sure that if perturbations are not computed, idm_dr and idr parameters are still freed */
    if(ppt->has_perturbations == _FALSE_) {

      if (ppt->alpha_idm_dr != NULL)
        free(ppt->alpha_idm_dr);

      if (ppt->beta_idr != NULL)
        free(ppt->beta_idr);
    }

    /* Test */
    class_call(parser_check_options(string1, options_output, 33, &flag1),
               errmsg,
               errmsg);
    class_test(flag1==_FALSE_,
               errmsg, "The options for output are {'tCl','pCl','lCl','nCl','dCl','sCl','mPk','mTk','dTk','vTk','Sd'}, you entered '%s'",string1);
  }

  /** 1.a) Terms contributing to the temperature spectrum */
  if (ppt->has_cl_cmb_temperature == _TRUE_) {
    /* Read */
    class_call(parser_read_string(pfc,"temperature_contributions",&string1,&flag1,errmsg),
               errmsg,
               errmsg);
    /* Compatibility code BEGIN */
    if(flag1 == _FALSE_){
      class_call(parser_read_string(pfc,"temperature contributions",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);
    }
    /* Compatibility code END */
    /* Complete set of parameters */
    if (flag1 == _TRUE_){
      ppt->switch_sw = 0;
      ppt->switch_eisw = 0;
      ppt->switch_lisw = 0;
      ppt->switch_dop = 0;
      ppt->switch_pol = 0;
      if ((strstr(string1,"tsw") != NULL) || (strstr(string1,"TSW") != NULL)){
        ppt->switch_sw = 1;
      }
      if ((strstr(string1,"eisw") != NULL) || (strstr(string1,"EISW") != NULL)){
        ppt->switch_eisw = 1;
      }
      if ((strstr(string1,"lisw") != NULL) || (strstr(string1,"LISW") != NULL)){
        ppt->switch_lisw = 1;
      }
      if ((strstr(string1,"dop") != NULL) || (strstr(string1,"Dop") != NULL)){
        ppt->switch_dop = 1;
      }
      if ((strstr(string1,"pol") != NULL) || (strstr(string1,"Pol") != NULL)){
        ppt->switch_pol = 1;
      }
      /* Test */
      class_call(parser_check_options(string1, options_temp_contributions, 10, &flag1),
                 errmsg,
                 errmsg);
      class_test(flag1==_FALSE_,
                 errmsg, "The options for 'temperature_contributions' are {'tsw','eisw','lisw','dop','pol'}, you entered '%s'",string1);
      class_test((ppt->switch_sw == 0) && (ppt->switch_eisw == 0) && (ppt->switch_lisw == 0) && (ppt->switch_dop == 0) && (ppt->switch_pol == 0),
                 errmsg,
                 "You specified 'temperature_contributions' as '%s'. It has to contain some of {'tsw','eisw','lisw','dop','pol'}.",string1);

      /** 1.a.1) Split value of redshift z at which the isw is considered as late or early */
      /* Read */
      class_read_double("early/late isw redshift",ppt->eisw_lisw_split_z); //Deprecated parameter
      class_read_double("early_late_isw_redshift",ppt->eisw_lisw_split_z);
    }
  }

  /** 1.b) Obsevable number count fluctuation spectrum */
  if (ppt->has_cl_number_count == _TRUE_){
    /* Read */
    class_call(parser_read_string(pfc,"number_count_contributions",&string1,&flag1,errmsg),
               errmsg,
               errmsg);
    /* Compatibility code BEGIN */
    if(flag1 == _FALSE_){
      class_call(parser_read_string(pfc,"number count contributions",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);
    }
    /* Compatibility code END */
    /* Complete set of parameters */
    if (flag1 == _TRUE_) {
      if (strstr(string1,"density") != NULL || strstr(string1,"dens") != NULL){
        ppt->has_nc_density = _TRUE_;
      }
      if (strstr(string1,"rsd") != NULL || strstr(string1,"RSD") != NULL){
        ppt->has_nc_rsd = _TRUE_;
      }
      if (strstr(string1,"lensing") != NULL || strstr(string1,"lens") != NULL){
        ppt->has_nc_lens = _TRUE_;
      }
      if (strstr(string1,"gr") != NULL || strstr(string1,"GR") != NULL){
        ppt->has_nc_gr = _TRUE_;
      }
      /* Test */
      class_call(parser_check_options(string1, options_number_count, 8, &flag1),
                 errmsg,
                 errmsg);
      class_test(flag1==_FALSE_,
                 errmsg, "The options for 'number_count_contributions' are {'density','rsd','lensing','gr'}, you entered '%s'",string1);
      class_test((ppt->has_nc_density == _FALSE_) && (ppt->has_nc_rsd == _FALSE_) && (ppt->has_nc_lens == _FALSE_) && (ppt->has_nc_gr == _FALSE_),
                 errmsg,
                 "You specified 'number_count_contributions' as '%s'. It has to contain some of {'density','rsd','lensing','gr'}.",string1);
    }
    else {
      /* Set default value */
      ppt->has_nc_density = _TRUE_;
    }
  }

  /** 1.c) Transfer function of additional metric fluctuations */
  if (ppt->has_density_transfers == _TRUE_) {
    /* Read */
    class_read_flag_or_deprecated("extra_metric_transfer_functions","extra metric transfer functions",ppt->has_metricpotential_transfers);
  }

  if (ppt->has_perturbations == _TRUE_) {

    /** 2) Perturbed recombination */
    /* Read */
    class_read_flag_or_deprecated("perturbed_recombination","perturbed recombination",ppt->has_perturbed_recombination);

    /** 3) Modes */
    /* Read */
    class_call(parser_read_string(pfc,"modes",&string1,&flag1,errmsg),
               errmsg,
               errmsg);
    /* Complete set of parameters */
    if (flag1 == _TRUE_) {
      /* if no modes are specified, the default is has_scalars=_TRUE_;
         but if they are specified we should reset has_scalars to _FALSE_ before reading */
      ppt->has_scalars=_FALSE_;
      if ((strstr(string1,"s") != NULL) || (strstr(string1,"S") != NULL)){
        ppt->has_scalars=_TRUE_;
      }
      if ((strstr(string1,"v") != NULL) || (strstr(string1,"V") != NULL)){
        ppt->has_vectors=_TRUE_;
      }
      if ((strstr(string1,"t") != NULL) || (strstr(string1,"T") != NULL)){
        ppt->has_tensors=_TRUE_;
      }
      /* Test */
      class_call(parser_check_options(string1, options_modes, 6, &flag1),
                 errmsg,
                 errmsg);
      class_test(flag1==_FALSE_,
                 errmsg, "The options for 'modes' are {'s','v','t'}, you entered '%s'",string1);
      class_test(class_none_of_three(ppt->has_scalars,ppt->has_vectors,ppt->has_tensors),
                 errmsg,
                 "You specified 'modes' as '%s'. It has to contain some of {'s','v','t'}.",string1);
    }
    /* Test */
    if (ppt->has_vectors == _TRUE_){
      class_test((ppt->has_cl_cmb_temperature == _FALSE_) && (ppt->has_cl_cmb_polarization == _FALSE_),
                 errmsg,
                 "Inconsistent input: you asked for vectors, so you should have at least one non-zero tensor source type (temperature or polarization). Please adjust your input.");
    }
    if (ppt->has_tensors == _TRUE_){
      class_test((ppt->has_cl_cmb_temperature == _FALSE_) && (ppt->has_cl_cmb_polarization == _FALSE_),
                 errmsg,
                 "Inconsistent input: you asked for tensors, so you should have at least one non-zero tensor source type (temperature or polarization). Please adjust your input.");
    }

    /** 3.a) List of initial conditions for scalars */
    if (ppt->has_scalars == _TRUE_) {
      /* Read */
      class_call(parser_read_string(pfc,"ic",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);
      /* Complete set of parameters */
      if (flag1 == _TRUE_) {
        /* if no initial conditions are specified, the default is has_ad=_TRUE_;
           but if they are specified we should reset has_ad to _FALSE_ before reading */
        ppt->has_ad=_FALSE_;
        if ((strstr(string1,"ad") != NULL) || (strstr(string1,"AD") != NULL)){
          ppt->has_ad=_TRUE_;
        }
        if ((strstr(string1,"bi") != NULL) || (strstr(string1,"BI") != NULL)){
          ppt->has_bi=_TRUE_;
        }
        if ((strstr(string1,"cdi") != NULL) || (strstr(string1,"CDI") != NULL)){
          ppt->has_cdi=_TRUE_;
        }
        if ((strstr(string1,"nid") != NULL) || (strstr(string1,"NID") != NULL)){
          ppt->has_nid=_TRUE_;
        }
        if ((strstr(string1,"niv") != NULL) || (strstr(string1,"NIV") != NULL)){
          ppt->has_niv=_TRUE_;
        }
        /* Test */
        class_call(parser_check_options(string1, options_ics, 10, &flag1),
                   errmsg,
                   errmsg);
        class_test(flag1==_FALSE_,
                   errmsg, "The options for 'ic' are {'ad','bi','cdi','nid','niv'}, you entered '%s'",string1);
        class_test(ppt->has_ad==_FALSE_ && ppt->has_bi ==_FALSE_ && ppt->has_cdi ==_FALSE_ && ppt->has_nid ==_FALSE_ && ppt->has_niv ==_FALSE_,
                   errmsg,
                   "You specified 'ic' as '%s'. It has to contain some of {'ad','bi','cdi','nid','niv'}.",string1);
      }
    }
    else {
      /* Test */
      class_test(ppt->has_cl_cmb_lensing_potential == _TRUE_,
                 errmsg,
                 "Inconsistency: you want C_l's for cmb lensing potential, but no scalar modes\n");
      class_test(ppt->has_pk_matter == _TRUE_,
                 errmsg,
                 "Inconsistency: you want P(k) of matter, but no scalar modes\n");
    }

    /** 3.b) List of initial conditions for scalars */
    if (ppt->has_tensors == _TRUE_) {
      /* Read */
      class_call(parser_read_string(pfc,"tensor_method",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);
      /* Compatibility code BEGIN */
      if(flag1 == _FALSE_){
        class_call(parser_read_string(pfc,"tensor method",&string1,&flag1,errmsg),
                   errmsg,
                   errmsg);
      }
      /* Compatibility code END */
      /* Complete set of parameters */
      if (flag1 == _TRUE_) {
        if (strstr(string1,"photons") != NULL){
          ppt->tensor_method = tm_photons_only;
        }
        else if (strstr(string1,"massless") != NULL){
          ppt->tensor_method = tm_massless_approximation;
        }
        else if (strstr(string1,"exact") != NULL){
          ppt->tensor_method = tm_exact;
        }
        else{
          class_stop(errmsg,"incomprehensible input '%s' for the field 'tensor_method'",string1);
        }
      }
    }
  }


  /** 4) Gauge */
  /** 4.a) Set gauge */
  /* Read */
  class_call(parser_read_string(pfc,"gauge",&string1,&flag1,errmsg),
             errmsg,
             errmsg);
  /* Complete set of parameters */
  if (flag1 == _TRUE_) {
    if ((strstr(string1,"newtonian") != NULL) || (strstr(string1,"Newtonian") != NULL) || (strstr(string1,"new") != NULL)) {
      ppt->gauge = newtonian;
    }
    else if ((strstr(string1,"synchronous") != NULL) || (strstr(string1,"sync") != NULL) || (strstr(string1,"Synchronous") != NULL)) {
      ppt->gauge = synchronous;
    }
    else{
      class_stop(errmsg,
                 "You specified the gauge as '%s'. It has to be one of {'newtonian','synchronous'}.");
    }
  }

  /** 4.b) Do we want density and velocity transfer functions in Nbody gauge? */
  if ((ppt->has_density_transfers == _TRUE_) || (ppt->has_velocity_transfers == _TRUE_)){

    /* Read */
    class_read_flag_or_deprecated("nbody_gauge_transfer_functions","Nbody gauge transfer functions",ppt->has_Nbody_gauge_transfers);

  }

  /** 5) h in [-] and H_0/c in [1/Mpc = h/2997.9 = h*10^5/c] */
  /* Read */
  class_call(parser_read_double(pfc,"H0",&param1,&flag1,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"h",&param2,&flag2,errmsg),
             errmsg,
             errmsg);
  /* Test */
  class_test((flag1 == _TRUE_) && (flag2 == _TRUE_),
             errmsg,
             "You can only enter one of 'h' or 'H0'.");
  /* Complete set of parameters */
  if (flag1 == _TRUE_){
    pba->H0 = param1*1.e3/_c_;
    pba->h = param1/100.;
  }
  if (flag2 == _TRUE_){
    pba->H0 = param2*1.e5/_c_;
    pba->h = param2;
  }


  /** 6) Primordial helium fraction */
  /* Read */
  class_call(parser_read_string(pfc,"YHe",&string1,&flag1,errmsg),
             errmsg,
             errmsg);
  /* Complete set of parameters */
  if (flag1 == _TRUE_) {
    if ((strstr(string1,"BBN") != NULL) || (strstr(string1,"bbn") != NULL)){
      pth->YHe = _YHE_BBN_;
    }
    else {
      class_read_double("YHe",pth->YHe);
    }
  }


  /** 7) Recombination parameters */
  /* Read */
  class_call(parser_read_string(pfc,"recombination",&string1,&flag1,errmsg),
             errmsg,
             errmsg);
  /* Complete set of parameters */
  if (flag1 == _TRUE_){
    if ((strstr(string1,"RECFAST") != NULL) || (strstr(string1,"recfast") != NULL) || (strstr(string1,"Recfast") != NULL)){
      pth->recombination = recfast;
    }
    else if ((strstr(string1,"HYREC") != NULL) || (strstr(string1,"hyrec") != NULL) || (strstr(string1,"HyRec") != NULL)){
      pth->recombination = hyrec;
    }
    else{
      class_stop(errmsg,
                 "You specified 'recombination' as '%s'. It has to be one of {'recfast','hyrec'}.",string1);
    }
  }

  /** 7.a) Photo-ionization dependence for recfast */
  /* Read */
  if(pth->recombination == recfast){
    class_call(parser_read_string(pfc,"recfast_photoion_dependence",&string1,&flag1,errmsg),
               errmsg,
               errmsg);
    if (flag1 == _TRUE_){
      if((strstr(string1,"Tmat") != NULL) || (strstr(string1,"tmat") != NULL ) || (strstr(string1,"TMAT") !=NULL)){
        pth->recfast_photoion_mode = recfast_photoion_Tmat;
      }
      else if((strstr(string1,"Trad") != NULL) || (strstr(string1,"trad") != NULL ) || (strstr(string1,"TRAD") !=NULL)){
        pth->recfast_photoion_mode = recfast_photoion_Trad;
      }
      else{
        class_stop(errmsg,
                   "You specified 'recfast_photoion_dependence' as '%s'. It has to be one of {'Tmat','Trad'}.",string1);
      }
    }
  }

  /** 8) Reionization parametrization */
  /* Read */
  class_call(parser_read_string(pfc,"reio_parametrization",&string1,&flag1,errmsg),
             errmsg,
             errmsg);
  /* Complete set of parameters */
  if (flag1 == _TRUE_){
    if (strcmp(string1,"reio_none") == 0){
      pth->reio_parametrization = reio_none;
    }
    else if (strcmp(string1,"reio_camb") == 0){
      pth->reio_parametrization = reio_camb;
    }
    else if (strcmp(string1,"reio_bins_tanh") == 0){
      pth->reio_parametrization = reio_bins_tanh;
    }
    else if (strcmp(string1,"reio_half_tanh") == 0){
      pth->reio_parametrization = reio_half_tanh;
    }
    else if (strcmp(string1,"reio_many_tanh") == 0){
      pth->reio_parametrization = reio_many_tanh;
    }
    else if (strcmp(string1,"reio_inter") == 0){
      pth->reio_parametrization = reio_inter;
    }
    else{
      class_stop(errmsg,
                 "You specified 'reio_parametrization' as '%s'. It has to be one of {'reio_none','reio_camb','reio_bins_tanh','reio_half_tanh','reio_many_tanh','reio_inter'}.",string1);
    }
  }

  switch (pth->reio_parametrization) {

  case reio_none:
    /* nothing to be read*/
    break;

  /** 8.a) Reionization parameters if reio_parametrization=reio_camb */
  case reio_camb:
  case reio_half_tanh:
    /* Read */
    class_call(parser_read_double(pfc,"z_reio",&param1,&flag1,errmsg),
               errmsg,
               errmsg);
    class_call(parser_read_double(pfc,"tau_reio",&param2,&flag2,errmsg),
               errmsg,
               errmsg);
    class_read_double("reionization_exponent",pth->reionization_exponent);
    class_read_double("reionization_width",pth->reionization_width);
    class_read_double("helium_fullreio_redshift",pth->helium_fullreio_redshift);
    class_read_double("helium_fullreio_width",pth->helium_fullreio_width);
    /* Test */
    class_test(((flag1 == _TRUE_) && (flag2 == _TRUE_)),
               errmsg,
               "You can only enter one of 'z_reio' or 'tau_reio'.");
    /* Complete set of parameters */
    if (flag1 == _TRUE_){
      pth->z_reio=param1;
      pth->reio_z_or_tau=reio_z;
    }
    if (flag2 == _TRUE_){
      pth->tau_reio=param2;
      pth->reio_z_or_tau=reio_tau;
    }
    break;

    /** 8.b) Reionization parameters if reio_parametrization=reio_bins_tanh */
  case reio_bins_tanh:
    /* Read */
    class_read_int("binned_reio_num",pth->binned_reio_num);
    class_read_list_of_doubles("binned_reio_z",pth->binned_reio_z,pth->binned_reio_num);
    class_read_list_of_doubles("binned_reio_xe",pth->binned_reio_xe,pth->binned_reio_num);
    class_read_double("binned_reio_step_sharpness",pth->binned_reio_step_sharpness);
    break;

    /** 8.c) reionization parameters if reio_parametrization=reio_many_tanh */
  case reio_many_tanh:
    /* Read */
    class_read_int("many_tanh_num",pth->many_tanh_num);
    class_read_list_of_doubles("many_tanh_z",pth->many_tanh_z,pth->many_tanh_num);
    class_read_list_of_doubles("many_tanh_xe",pth->many_tanh_xe,pth->many_tanh_num);
    class_read_double("many_tanh_width",pth->many_tanh_width);
    break;

    /** 8.d) reionization parameters if reio_parametrization=reio_many_tanh */
  case reio_inter:
    /* Read */
    class_read_int("reio_inter_num",pth->reio_inter_num);
    class_read_list_of_doubles("reio_inter_z",pth->reio_inter_z,pth->reio_inter_num);
    class_read_list_of_doubles("reio_inter_xe",pth->reio_inter_xe,pth->reio_inter_num);
    break;

  default:
    class_stop(pth->error_message,"pth->recombination=%d different from all known cases",pth->recombination);
    break;
  }

  /** 9) Damping scale */
  /* Read */
  class_read_flag_or_deprecated("compute_damping_scale","compute damping scale",pth->compute_damping_scale);

  return _SUCCESS_;

}


/**
 * Read the parameters for each physical species
 *
 * @param pfc            Input: pointer to local structure
 * @param ppr            Input: pointer to precision structure
 * @param pba            Input: pointer to background structure
 * @param pth            Input: pointer to thermodynamics structure
 * @param ppt            Input: pointer to perturbation structure
 * @param input_verbose  Input: verbosity of input
 * @param errmsg         Input: Error message
 * @return the error status
 */

int input_read_parameters_species(struct file_content * pfc,
                                  struct precision * ppr,
                                  struct background * pba,
                                  struct thermodynamics * pth,
                                  struct perturbations * ppt,
                                  int input_verbose,
                                  ErrorMsg errmsg){

  /** Summary: */

  /** - Define local variables */
  int flag1, flag2, flag3;
  double param1, param2, param3;
  char string1[_ARGUMENT_LENGTH_MAX_];
  int fileentries;
  int N_ncdm=0, n, entries_read;
  double rho_ncdm;
  double scf_lambda;
  double fnu_factor;
  double Omega_tot;
  double sigma_B; // Stefan-Boltzmann constant
  double stat_f_idr = 7./8.;



  sigma_B = 2.*pow(_PI_,5.)*pow(_k_B_,4.)/15./pow(_h_P_,3.)/pow(_c_,2);  // [W/(m^2 K^4) = Kg/(K^4 s^3)]

  /** 1) Omega_0_g (photons) and T_cmb */
  /* Read */
  class_call(parser_read_double(pfc,"T_cmb",&param1,&flag1,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"Omega_g",&param2,&flag2,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"omega_g",&param3,&flag3,errmsg),
             errmsg,
             errmsg);
  class_test(class_at_least_two_of_three(flag1,flag2,flag3),
             errmsg,
             "You can only enter one of 'T_cmb', 'Omega_g' or 'omega_g'.");
  /* Complete set of parameters
     Note:  Omega0_g = rho_g/rho_c0, each of them expressed in [Kg/(m s^2)]
            rho_g = (4 sigma_B/c) T^4
            rho_c0 = 3 c^2 H_0^2/(8 \pi G) */
  if (class_none_of_three(flag1,flag2,flag3)){
    pba->Omega0_g = (4.*sigma_B/_c_*pow(pba->T_cmb,4.))/(3.*_c_*_c_*1.e10*pba->h*pba->h/_Mpc_over_m_/_Mpc_over_m_/8./_PI_/_G_);
  }
  else {
    if (flag1 == _TRUE_){
      pba->Omega0_g = (4.*sigma_B/_c_*pow(param1,4.))/(3.*_c_*_c_*1.e10*pba->h*pba->h/_Mpc_over_m_/_Mpc_over_m_/8./_PI_/_G_);
      pba->T_cmb=param1;
    }
    if (flag2 == _TRUE_){
      pba->Omega0_g = param2;
      pba->T_cmb = pow(pba->Omega0_g*(3.*_c_*_c_*1.e10*pba->h*pba->h/_Mpc_over_m_/_Mpc_over_m_/8./_PI_/_G_)/(4.*sigma_B/_c_),0.25);
    }
    if (flag3 == _TRUE_){
      pba->Omega0_g = param3/pba->h/pba->h;
      pba->T_cmb = pow(pba->Omega0_g*(3.*_c_*_c_*1.e10*pba->h*pba->h/_Mpc_over_m_/_Mpc_over_m_/8./_PI_/_G_)/(4.*sigma_B/_c_),0.25);
    }
  }


  /** 2) Omega_0_b (baryons) */
  /* Read */
  class_call(parser_read_double(pfc,"Omega_b",&param1,&flag1,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"omega_b",&param2,&flag2,errmsg),
             errmsg,
             errmsg);
  /* Test */
  class_test(((flag1 == _TRUE_) && (flag2 == _TRUE_)),
             errmsg,
             "You can only enter one of 'Omega_b' or 'omega_b'.");
  /* Complete set of parameters */
  if (flag1 == _TRUE_){
    pba->Omega0_b = param1;
  }
  if (flag2 == _TRUE_){
    pba->Omega0_b = param2/pba->h/pba->h;
  }


  /** 3) Omega_0_ur (ultra-relativistic species / massless neutrino) */
  /**
   * We want to keep compatibility with old input files, and as such 'N_eff' is still
   * an allowed parameter name, although it is deprecated and its use is discouraged.
   * */
  /* Read */
  class_call(parser_read_double(pfc,"N_ur",&param1,&flag1,errmsg),
             errmsg,
             errmsg);
  /* Compability code BEGIN */
  class_call(parser_read_double(pfc,"N_eff",&param2,&flag2,errmsg),
               errmsg,
               errmsg);
  class_test((flag1 == _TRUE_) && (flag2 == _TRUE_),
             errmsg,
             "You added both 'N_eff' (deprecated) and 'N_ur'. Please use solely 'N_ur'.");
  if(flag2 == _TRUE_){
    param1 = param2;
    flag1 = _TRUE_;
  }
  /* Compability code END */
  class_call(parser_read_double(pfc,"Omega_ur",&param2,&flag2,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"omega_ur",&param3,&flag3,errmsg),
             errmsg,
             errmsg);
  /* Test */
  class_test(class_at_least_two_of_three(flag1,flag2,flag3),
             errmsg,
             "You can only enter one of 'N_ur', 'Omega_ur' or 'omega_ur'.");
  /* Complete set of parameters assuming as default value N_eff=3.044
      (see 2008.01074 and 2012.02726. This value is more accurate than
      the previous default value of 3.046) */
  if (class_none_of_three(flag1,flag2,flag3)) {
    pba->Omega0_ur = 3.044*7./8.*pow(4./11.,4./3.)*pba->Omega0_g;
  }
  else {
    if (flag1 == _TRUE_) {
      pba->Omega0_ur = param1*7./8.*pow(4./11.,4./3.)*pba->Omega0_g;
    }
    if (flag2 == _TRUE_) {
      pba->Omega0_ur = param2;
    }
    if (flag3 == _TRUE_) {
      pba->Omega0_ur = param3/pba->h/pba->h;
    }
  }

  /** 3.a) Case of non-standard properties */
  /* Read */
  class_call(parser_read_double(pfc,"ceff2_ur",&param1,&flag1,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"cvis2_ur",&param2,&flag2,errmsg),
             errmsg,
             errmsg);
  /* Complete set of parameters */
  if (flag1 == _TRUE_){
    ppt->three_ceff2_ur = 3.*param1;
  }
  if (flag2 == _TRUE_){
    ppt->three_cvis2_ur = 3.*param2;
  }


  /** 4) Omega_0_cdm (CDM) */
  /* Read */
  class_call(parser_read_double(pfc,"Omega_cdm",&param1,&flag1,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"omega_cdm",&param2,&flag2,errmsg),
             errmsg,
             errmsg);
  /* Test */
  class_test(((flag1 == _TRUE_) && (flag2 == _TRUE_)),
             errmsg,
             "You can only enter one of 'Omega_cdm' or 'omega_cdm'.");
  /* Complete set of parameters */
  if (flag1 == _TRUE_){
    pba->Omega0_cdm = param1;
  }
  if (flag2 == _TRUE_){
    pba->Omega0_cdm = param2/pba->h/pba->h;
  }

  if ((ppt->gauge == synchronous) && (pba->Omega0_cdm==0)) pba->Omega0_cdm = ppr->Omega0_cdm_min_synchronous;

  /** 5) Non-cold relics (ncdm) */
  /** 5.a) Number of non-cold relics */
  /* Read */
  class_read_int("N_ncdm",N_ncdm);
  /* Complete set of parameters */
  if (N_ncdm > 0){
    pba->N_ncdm = N_ncdm;
    if (ppt->gauge == synchronous){
      ppr->tol_ncdm = ppr->tol_ncdm_synchronous;
    }
    if (ppt->gauge == newtonian){
      ppr->tol_ncdm = ppr->tol_ncdm_newtonian;
    }

    /** 5.b) Check if filenames for interpolation tables are given */
    /* Read */
    class_read_list_of_integers_or_default("use_ncdm_psd_files",pba->got_files,_FALSE_,N_ncdm);
    /* Complete set of parameters */
    for(n=0,fileentries=0; n<N_ncdm; n++){
      if (pba->got_files[n] == _TRUE_){
        fileentries++;
      }
    }
    if (fileentries > 0) {

      /** 5.b.1) Check if filenames for interpolation tables are given */
      /* Read */
      class_call(parser_read_list_of_strings(pfc,"ncdm_psd_filenames",&entries_read,&(pba->ncdm_psd_files),&flag1,errmsg),
                 errmsg,
                 errmsg);
      /* Test */
      class_test(flag1 == _FALSE_,errmsg,
                 "Entry 'use_ncdm_files' is found, but no corresponding 'ncdm_psd_filenames' were found.");
      class_test(entries_read != fileentries,errmsg,
                 "Number of filenames found (%d) does not match number of _TRUE_ values in use_ncdm_files (%d).",
                 entries_read,fileentries);
    }

    /** 5.c) (optional) p.s.d.-parameters */
    /* Read */
    parser_read_list_of_doubles(pfc,"ncdm_psd_parameters",&entries_read,&(pba->ncdm_psd_parameters),&flag1,errmsg);

    /** 5.d) Mass or Omega of each ncdm species */
    /* Read */
    class_read_list_of_doubles_or_default("m_ncdm",pba->m_ncdm_in_eV,0.0,N_ncdm);
    class_read_list_of_doubles_or_default("Omega_ncdm",pba->Omega0_ncdm,0.0,N_ncdm);
    class_read_list_of_doubles_or_default("omega_ncdm",pba->M_ncdm,0.0,N_ncdm);
    for(n=0; n<N_ncdm; n++){
      if (pba->M_ncdm[n]!=0.0){
        /* Test */
        class_test(pba->Omega0_ncdm[n]!=0,errmsg,
                   "You can only enter one of 'Omega_ncdm' or 'omega_ncdm' for ncdm species %d.",n);
        /* Complete set of parameters */
        pba->Omega0_ncdm[n] = pba->M_ncdm[n]/pba->h/pba->h;
      }
      /* Set default value
         this is the right place for passing the default value of the mass
         (all parameters must have a default value; most of them are defined
         in input_default_params, but the ncdm mass is a bit special and
         there is no better place for setting its default value). We put an
         arbitrary value m << 10^-3 eV, i.e. the ultra-relativistic limit. */
      if ((pba->Omega0_ncdm[n]==0.0) && (pba->m_ncdm_in_eV[n]==0.0)) {
        pba->m_ncdm_in_eV[n]=1.e-5;
      }
    }

    /** 5.e) Temperatures */
    /* Read */
    class_read_list_of_doubles_or_default("T_ncdm",pba->T_ncdm,pba->T_ncdm_default,N_ncdm);

    /** 5.f) Chemical potentials */
    /* Read */
    class_read_list_of_doubles_or_default("ksi_ncdm",pba->ksi_ncdm,pba->ksi_ncdm_default,N_ncdm);

    /** 5.g) Degeneracy of each ncdm species */
    /* Read */
    class_read_list_of_doubles_or_default("deg_ncdm",pba->deg_ncdm,pba->deg_ncdm_default,N_ncdm);

    /** 5.h) Quadrature modes, 0 is qm_auto */
    /* Read */
    class_read_list_of_integers_or_default("Quadrature strategy",pba->ncdm_quadrature_strategy,0,N_ncdm); //Deprecated parameter, still read to keep compatibility
    class_read_list_of_integers_or_default("ncdm_quadrature_strategy",pba->ncdm_quadrature_strategy,0,N_ncdm);

    /** 5.h.1) qmax, if relevant */
    /* Read */
    class_read_list_of_doubles_or_default("Maximum q",pba->ncdm_qmax,15,N_ncdm); //Deprecated parameter, still read to keep compatibility
    class_read_list_of_doubles_or_default("ncdm_maximum_q",pba->ncdm_qmax,15,N_ncdm);

    /** 5.h.2) Number of momentum bins */
    class_read_list_of_integers_or_default("Number of momentum bins",pba->ncdm_input_q_size,150,N_ncdm); //Deprecated parameter, still read to keep compatibility
    class_read_list_of_integers_or_default("ncdm_N_momentum_bins",pba->ncdm_input_q_size,150,N_ncdm);

    /** Last step of 5) (i.e. NCDM) -- Calculate the masses and momenta */
    class_call(background_ncdm_init(ppr,pba),
               pba->error_message,
               errmsg);
    /* Complete set of parameters
     We must calculate M from omega or vice versa if one of them is missing.
     If both are present, we must update the degeneracy parameter to
     reflect the implicit normalization of the distribution function. */
    for (n=0; n < N_ncdm; n++){
      if (pba->m_ncdm_in_eV[n] != 0.0){
        /* Case of only mass or mass and Omega/omega: */
        pba->M_ncdm[n] = pba->m_ncdm_in_eV[n]/_k_B_*_eV_/pba->T_ncdm[n]/pba->T_cmb;
        class_call(background_ncdm_momenta(pba->q_ncdm_bg[n],
                                           pba->w_ncdm_bg[n],
                                           pba->q_size_ncdm_bg[n],
                                           pba->M_ncdm[n],
                                           pba->factor_ncdm[n],
                                           0.,
                                           NULL,
                                           &rho_ncdm,
                                           NULL,
                                           NULL,
                                           NULL),
                   pba->error_message,
                   errmsg);
        if (pba->Omega0_ncdm[n] == 0.0){
          pba->Omega0_ncdm[n] = rho_ncdm/pba->H0/pba->H0;
        }
        else{
          fnu_factor = (pba->H0*pba->H0*pba->Omega0_ncdm[n]/rho_ncdm);
          pba->factor_ncdm[n] *= fnu_factor;
          /* dlnf0dlnq is already computed, but it is independent of any
             normalization of f0. We don't need the factor anymore, but we
             store it nevertheless */
          pba->deg_ncdm[n] *=fnu_factor;
        }
      }
      else{
        /* Case of only Omega/omega: */
        class_call(background_ncdm_M_from_Omega(ppr,pba,n),
                   pba->error_message,
                   errmsg);
        pba->m_ncdm_in_eV[n] = _k_B_/_eV_*pba->T_ncdm[n]*pba->M_ncdm[n]*pba->T_cmb;
      }
      pba->Omega0_ncdm_tot += pba->Omega0_ncdm[n];
    }

  }


  /** 6) Omega_0_k (effective fractional density of curvature) */
  /* Read */
  class_read_double("Omega_k",pba->Omega0_k);
  /* Complete set of parameters */
  pba->K = -pba->Omega0_k*pow(pba->H0,2);
  if (pba->K > 0.){
    pba->sgnK = 1;
  }
  else if (pba->K < 0.){
    pba->sgnK = -1;
  }


  /* 7) ** ADDITIONAL SPECIES ** --> Add your species here */

  /** 7.1) Decaying DM into DR */
  /** 7.1.a) Omega_0_dcdmdr (DCDM, i.e. decaying CDM) */
  /* Read */
  class_call(parser_read_double(pfc,"Omega_dcdmdr",&param1,&flag1,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"omega_dcdmdr",&param2,&flag2,errmsg),
             errmsg,
             errmsg);
  /* Test */
  class_test(((flag1 == _TRUE_) && (flag2 == _TRUE_)),
             errmsg,
             "You can only enter one of 'Omega_dcdmdr' or 'omega_dcdmdr'.");
  /* Complete set of parameters */
  if (flag1 == _TRUE_){
    pba->Omega0_dcdmdr = param1;
  }
  if (flag2 == _TRUE_){
    pba->Omega0_dcdmdr = param2/pba->h/pba->h;
  }

  if (pba->Omega0_dcdmdr > 0) {
    /** 7.1.b) Omega_ini_dcdm or omega_ini_dcdm */
    /* Read */
    class_call(parser_read_double(pfc,"Omega_ini_dcdm",&param1,&flag1,errmsg),
               errmsg,
               errmsg);
    class_call(parser_read_double(pfc,"omega_ini_dcdm",&param2,&flag2,errmsg),
               errmsg,
               errmsg);
    /* Test */
    class_test(((flag1 == _TRUE_) && (flag2 == _TRUE_)),
               errmsg,
               "You can only enter one of 'Omega_ini_dcdm' or 'omega_ini_dcdm'.");
    /* Complete set of parameters */
    if (flag1 == _TRUE_){
      pba->Omega_ini_dcdm = param1;
    }
    if (flag2 == _TRUE_){
      pba->Omega_ini_dcdm = param2/pba->h/pba->h;
    }


    /** 7.1.c) Gamma in same units as H0, i.e. km/(s Mpc)*/
    /* Read */
    class_call(parser_read_double(pfc,"Gamma_dcdm",&param1,&flag1,errmsg),                          // [km/(s Mpc)]
               errmsg,
               errmsg);
    class_call(parser_read_double(pfc,"tau_dcdm",&param2,&flag2,errmsg),                            // [s]
               errmsg,
               errmsg);
    /* Test */
    class_test(((flag1 == _TRUE_) && (flag2 == _TRUE_)),
               errmsg,
               "In input file, you can only enter one of Gamma_dcdm or tau_dcdm, choose one");
    /* Complete set of parameters */
    if (flag1 == _TRUE_){
      pba->Gamma_dcdm = param1*(1.e3/_c_);                                                          // [Mpc]
      pba->tau_dcdm = 1/(param1*1.02e-3)*(1e9*365*24*3600);                                         // [s]
    }//TODO :: fix these factors to be proper (i.e. no 1.02e-3)
    if (flag2 == _TRUE_){
      pba->Gamma_dcdm = 1/(param2/(1e9*365*24*3600))/1.02e-3*(1.e3 / _c_);                          // [Mpc]
      pba->tau_dcdm = param2;                                                                       // [s]
    }
    /* Test */
    class_test(pba->tau_dcdm<0.,
               errmsg,
               "You need to enter a lifetime for the decaying DM 'tau_dcdm > 0.'");
    class_test(pba->Gamma_dcdm<0.,
               errmsg,
               "You need to enter a decay constant for the decaying DM 'Gamma_dcdm > 0.'");
  }

  /** 7.2) Interacting dark matter & dark radiation, ETHOS-parametrization/NADM parametrization, see explanatory.ini */
  /** 7.2.a) Omega_0_idr  */

  /* Read */
  class_call(parser_read_double(pfc,"N_idr",&param1,&flag1,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"N_dg",&param2,&flag2,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"xi_idr",&param3,&flag3,errmsg),
             errmsg,
             errmsg);
  class_test(class_at_least_two_of_three(flag1,flag2,flag3),
             errmsg,
             "In input file, you can only enter one of N_idr, N_dg or xi_idr, choose one");

  /** 7.2.b) stat_f_idr  */
  class_read_double("stat_f_idr",stat_f_idr);

  if (flag1 == _TRUE_) {
    pba->T_idr = pow(param1/stat_f_idr*(7./8.)/pow(11./4.,(4./3.)),(1./4.)) * pba->T_cmb;
    if (input_verbose > 1)
      printf("You passed N_idr = N_dg = %e, this is equivalent to xi_idr = %e in the ETHOS notation. \n", param2, pba->T_idr/pba->T_cmb);
  }
  else if (flag2 == _TRUE_) {
    pba->T_idr = pow(param2/stat_f_idr*(7./8.)/pow(11./4.,(4./3.)),(1./4.)) * pba->T_cmb;
    if (input_verbose > 2)
      printf("You passed N_dg = N_idr = %e, this is equivalent to xi_idr = %e in the ETHOS notation. \n", param2, pba->T_idr/pba->T_cmb);
  }
  else if (flag3 == _TRUE_) {
    pba->T_idr = param3 * pba->T_cmb;
    if (input_verbose > 1)
      printf("You passed xi_idr = %e, this is equivalent to N_idr = N_dg = %e in the NADM notation. \n", param3, stat_f_idr*pow(param3,4.)/(7./8.)*pow(11./4.,(4./3.)));
  }

  pba->Omega0_idr = stat_f_idr*pow(pba->T_idr/pba->T_cmb,4.)*pba->Omega0_g;

  /** - Omega_0_idm_dr (DM interacting with DR) */
  class_call(parser_read_double(pfc,"Omega_idm_dr",&param1,&flag1,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"omega_idm_dr",&param2,&flag2,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"f_idm_dr",&param3,&flag3,errmsg),
             errmsg,
             errmsg);
  class_test(class_at_least_two_of_three(flag1,flag2,flag3),
             errmsg,
             "In input file, you can only enter one of {Omega_idm_dr, omega_idm_dr, f_idm_dr} choose one");

  /* ---> if user passes directly the density of idm_dr */
  if (flag1 == _TRUE_)
    pba->Omega0_idm_dr = param1;
  if (flag2 == _TRUE_)
    pba->Omega0_idm_dr = param2/pba->h/pba->h;

  /* ---> if user passes density of idm_dr as a fraction of the CDM one */
  if (flag3 == _TRUE_) {
    class_test((param3 < 0.) || (param3 > 1.),
               errmsg,
               "The fraction of interacting DM with DR must be between 0 and 1, you asked for f_idm_dr=%e",param3);
    class_test((param3 > 0.) && (pba->Omega0_cdm == 0.),
               errmsg,
               "If you want a fraction of interacting DM with DR, to be consistent, you should not set the fraction of CDM to zero");

    pba->Omega0_idm_dr = param3 * pba->Omega0_cdm;
    /* readjust Omega0_cdm */
    pba->Omega0_cdm -= pba->Omega0_idm_dr;
    /* avoid Omega0_cdm =0 in synchronous gauge */
    if ((ppt->gauge == synchronous) && (pba->Omega0_cdm==0)) {
      pba->Omega0_cdm += ppr->Omega0_cdm_min_synchronous;
      pba->Omega0_idm_dr -= ppr->Omega0_cdm_min_synchronous;
    }
  }

  /* Test */
  if (pba->Omega0_idm_dr > 0.) {

    class_test(pba->Omega0_idr == 0.0,
               errmsg,
               "You have requested interacting DM ith DR, this requires a non-zero density of interacting DR. Please set either N_idr or xi_idr");
    /** 7.2.d) */
    class_read_double_one_of_two("m_idm","m_dm",pth->m_idm_dr);

    /** 7.2.e) */
    class_call(parser_read_double(pfc,"a_idm_dr",&param1,&flag1,errmsg),
               errmsg,
               errmsg);
    class_call(parser_read_double(pfc,"a_dark",&param2,&flag2,errmsg),
               errmsg,
               errmsg);
    class_call(parser_read_double(pfc,"Gamma_0_nadm",&param3,&flag3,errmsg),
               errmsg,
               errmsg);
    class_test(class_at_least_two_of_three(flag1,flag2,flag3),
               errmsg,
               "In input file, you can only enter one of a_idm_dr, a_dark or Gamma_0_nadm, choose one");

    if (flag1 == _TRUE_){
      pth->a_idm_dr = param1;
      if (input_verbose > 1)
        printf("You passed a_idm_dr = a_dark = %e, this is equivalent to Gamma_0_nadm = %e in the NADM notation. \n", param1, param1*(4./3.)*(pba->h*pba->h*pba->Omega0_idr));
    }
    else if (flag2 == _TRUE_){
      pth->a_idm_dr = param2;
      if (input_verbose > 1)
        printf("You passed a_dark = a_idm_dr = %e, this is equivalent to Gamma_0_nadm = %e in the NADM notation. \n", param2, param2*(4./3.)*(pba->h*pba->h*pba->Omega0_idr));
    }
    else if (flag3 == _TRUE_){
      pth->a_idm_dr = param3*(3./4.)/(pba->h*pba->h*pba->Omega0_idr);
      if (input_verbose > 1)
        printf("You passed Gamma_0_nadm = %e, this is equivalent to a_idm_dr = a_dark = %e in the ETHOS notation. \n", param3, pth->a_idm_dr);
    }

    /** 7.2.e.3/4) */
    /* If the user passed Gamma_0_nadm, assume they want nadm parameterisation*/
    if (flag3 == _TRUE_){
      /** Simply set 7.2.e.3/4) */
      pth->nindex_idm_dr = 0;
      ppt->idr_nature = idr_fluid;
      if (input_verbose > 1)
        printf("NADM requested. Defaulting on nindex_idm_dr = %e and idr_nature = fluid \n", pth->nindex_idm_dr);
    }

    /* If the user passed something else, assume they want ETHOS parameterisation*/
    else{

      /** 7.2.e.3) n_index_idm_dr */
      class_read_double_one_of_two("nindex_dark","nindex_idm_dr",pth->nindex_idm_dr);

      /** 7.2.e.4) idr_nature */ class_call(parser_read_string(pfc,"idr_nature",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);

      if (flag1 == _TRUE_) {
        if ((strstr(string1,"free_streaming") != NULL) || (strstr(string1,"Free_Streaming") != NULL) || (strstr(string1,"Free_streaming") != NULL) || (strstr(string1,"FREE_STREAMING") != NULL)) {
          ppt->idr_nature = idr_free_streaming;
        }
        if ((strstr(string1,"fluid") != NULL) || (strstr(string1,"Fluid") != NULL) || (strstr(string1,"FLUID") != NULL)) {
          ppt->idr_nature = idr_fluid;
        }
      }
    }

    /** 7.2.f) Strength of self interactions */
    class_read_double_one_of_two("b_dark","b_idr",pth->b_idr);


  // [NS] :: TODO :: This has to be fixed. For now, we simply always allocate it (exiting the if statement with a '}' and going into an always executed block with '{'.
  }{ // <--------- Bad practice, has to be corrected ASAP

    /** 7.2.g) Read alpha_idm_dr or alpha_dark */
    class_call(parser_read_list_of_doubles(pfc,"alpha_idm_dr",&entries_read,&(ppt->alpha_idm_dr),&flag1,errmsg),
               errmsg,
               errmsg);

    /* try with the other syntax */
    if (flag1 == _FALSE_) {
      class_call(parser_read_list_of_doubles(pfc,"alpha_dark",&entries_read,&(ppt->alpha_idm_dr),&flag1,errmsg),
                 errmsg,
                 errmsg);
    }
    // [NS] :: TODO :: move to perturbations module this allocation somehow? + fix bug that happens when entries_read == l_max_idr-1
    /* Only allocate if perturbations module will be called (otherwise segfaults) */
    if (ppt->has_perturbations) {
      if(flag1 == _TRUE_){
        if(entries_read != (ppr->l_max_idr-1)){
          class_realloc(ppt->alpha_idm_dr,ppt->alpha_idm_dr,(ppr->l_max_idr-1)*sizeof(double),errmsg);
          for(n=entries_read; n<(ppr->l_max_idr-1); n++) ppt->alpha_idm_dr[n] = ppt->alpha_idm_dr[entries_read-1];
        }
      }
      else{
        class_alloc(ppt->alpha_idm_dr,(ppr->l_max_idr-1)*sizeof(double),errmsg);
        for(n=0; n<(ppr->l_max_idr-1); n++) ppt->alpha_idm_dr[n] = 1.5;
      }
    }

    /* 7.2.h) Read beta_idm_dr or beta_dark */

    class_call(parser_read_list_of_doubles(pfc,"beta_idr",&entries_read,&(ppt->beta_idr),&flag1,errmsg),
               errmsg,
               errmsg);

    /* try with the other syntax */
    if (flag1 == _FALSE_) {
      class_call(parser_read_list_of_doubles(pfc,"beta_dark",&entries_read,&(ppt->beta_idr),&flag1,errmsg),
                 errmsg,
                 errmsg);
    }

    // [NS] :: TODO :: move to perturbations module this allocation somehow? + fix bug that happens when entries_read == l_max_idr-1
    /* Only allocate if perturbations module will be called (otherwise segfaults) */
    if (ppt->has_perturbations) {
      if(flag1 == _TRUE_){
        if(entries_read != (ppr->l_max_idr-1)){
          class_realloc(ppt->beta_idr,ppt->beta_idr,(ppr->l_max_idr-1)*sizeof(double),errmsg);
          for(n=entries_read; n<(ppr->l_max_idr-1); n++) ppt->beta_idr[n] = ppt->beta_idr[entries_read-1];
        }
      }
      else{
        class_alloc(ppt->beta_idr,(ppr->l_max_idr-1)*sizeof(double),errmsg);
        for(n=0; n<(ppr->l_max_idr-1); n++) ppt->beta_idr[n] = 1.5;
      }
    }
    else {
      ppt->alpha_idm_dr = NULL;
      ppt->beta_idr = NULL;
    }
  }

  /* ** ADDITIONAL SPECIES ** */


  /* At this point all the species should be set, and used for the budget equation below */
  /** 8) Dark energy
         Omega_0_lambda (cosmological constant), Omega0_fld (dark energy
         fluid), Omega0_scf (scalar field) */
  /* Read */
  class_call(parser_read_double(pfc,"Omega_Lambda",&param1,&flag1,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"Omega_fld",&param2,&flag2,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"Omega_scf",&param3,&flag3,errmsg),
             errmsg,
             errmsg);
  /* Test */
  class_test((flag1 == _TRUE_) && (flag2 == _TRUE_) && ((flag3 == _FALSE_) || (param3 >= 0.)),
             errmsg,
             "'Omega_Lambda' or 'Omega_fld' must be left unspecified, except if 'Omega_scf' is set and < 0.");
  class_test(((flag1 == _FALSE_)||(flag2 == _FALSE_)) && ((flag3 == _TRUE_) && (param3 < 0.)),
             errmsg,
             "You have entered 'Omega_scf' < 0 , so you have to specify both 'Omega_lambda' and 'Omega_fld'.");
  /* Complete set of parameters
     Case of (flag3 == _FALSE_) || (param3 >= 0.) means that either we have not
     read Omega_scf so we are ignoring it (unlike lambda and fld!) OR we have
     read it, but it had a positive value and should not be used for filling.
     We now proceed in two steps:
        1) set each Omega0 and add to the total for each specified component.
        2) go through the components in order {lambda, fld, scf} and fill using
           first unspecified component. */

  /* ** BUDGET EQUATION ** -> Add your species here */
  /* Compute Omega_tot */
  Omega_tot = pba->Omega0_g;
  Omega_tot += pba->Omega0_b;
  Omega_tot += pba->Omega0_ur;
  Omega_tot += pba->Omega0_cdm;
  Omega_tot += pba->Omega0_idm_dr;
  Omega_tot += pba->Omega0_idr;
  Omega_tot += pba->Omega0_ncdm_tot;
  /* Step 1 */
  if (flag1 == _TRUE_){
    pba->Omega0_lambda = param1;
    Omega_tot += pba->Omega0_lambda;
  }
  if (flag2 == _TRUE_){
    pba->Omega0_fld = param2;
    Omega_tot += pba->Omega0_fld;
  }
  if ((flag3 == _TRUE_) && (param3 >= 0.)){
    pba->Omega0_scf = param3;
    Omega_tot += pba->Omega0_scf;
  }
  /* Step 2 */
  if (flag1 == _FALSE_) {
    /* Fill with Lambda */
    pba->Omega0_lambda= 1. - pba->Omega0_k - Omega_tot;
    if (input_verbose > 0){
      printf(" -> matched budget equations by adjusting Omega_Lambda = %g\n",pba->Omega0_lambda);
    }
  }
  else if (flag2 == _FALSE_) {
    /* Fill up with fluid */
    pba->Omega0_fld = 1. - pba->Omega0_k - Omega_tot;
    if (input_verbose > 0){
      printf(" -> matched budget equations by adjusting Omega_fld = %g\n",pba->Omega0_fld);
    }
  }
  else if ((flag3 == _TRUE_) && (param3 < 0.)){
    /* Fill up with scalar field */
    pba->Omega0_scf = 1. - pba->Omega0_k - Omega_tot;
    if (input_verbose > 0){
      printf(" -> matched budget equations by adjusting Omega_scf = %g\n",pba->Omega0_scf);
    }
  }

  /* ** END OF BUDGET EQUATION ** */

  /** 8.a) If Omega fluid is different from 0 */
  if (pba->Omega0_fld != 0.) {
    /** 8.a.1) PPF approximation */
    /* Read */
    class_call(parser_read_string(pfc,"use_ppf",&string1,&flag1,errmsg),
               errmsg,
               errmsg);
    if (flag1 == _TRUE_){
      if(string_begins_with(string1,'y') || string_begins_with(string1,'Y')){
        pba->use_ppf = _TRUE_;
        class_read_double("c_gamma_over_c_fld",pba->c_gamma_over_c_fld);
      }
      else {
        pba->use_ppf = _FALSE_;
      }
    }

    /** 8.a.2) Equation of state */
    /* Read */
    class_call(parser_read_string(pfc,"fluid_equation_of_state",&string1,&flag1,errmsg),
               errmsg,
               errmsg);
    /* Complete set of parameters */
    if (flag1 == _TRUE_) {
      if ((strstr(string1,"CLP") != NULL) || (strstr(string1,"clp") != NULL)) {
        pba->fluid_equation_of_state = CLP;
      }
      else if ((strstr(string1,"EDE") != NULL) || (strstr(string1,"ede") != NULL)) {
        pba->fluid_equation_of_state = EDE;
      }
      else {
        class_stop(errmsg,"incomprehensible input '%s' for the field 'fluid_equation_of_state'",string1);
      }
    }

    if (pba->fluid_equation_of_state == CLP) {
      /** 8.a.2.2) Equation of state of the fluid in 'CLP' case */
      /* Read */
      class_read_double("w0_fld",pba->w0_fld);
      class_read_double("wa_fld",pba->wa_fld);
      class_read_double("cs2_fld",pba->cs2_fld);
    }
    if (pba->fluid_equation_of_state == EDE) {
      /** 8.a.2.3) Equation of state of the fluid in 'EDE' case */
      /* Read */
      class_read_double("w0_fld",pba->w0_fld);
      class_read_double("Omega_EDE",pba->Omega_EDE);
      class_read_double("cs2_fld",pba->cs2_fld);
    }
  }

  /** 8.b) If Omega scalar field (SCF) is different from 0 */
  if (pba->Omega0_scf != 0.){

    /** 8.b.1) Additional SCF parameters */
    /* Read */
    class_call(parser_read_list_of_doubles(pfc,
                                           "scf_parameters",
                                           &(pba->scf_parameters_size),
                                           &(pba->scf_parameters),
                                           &flag1,
                                           errmsg),
               errmsg,errmsg);

    /** 8.b.2) SCF initial conditions from attractor solution */
    /* Read */
    class_call(parser_read_string(pfc,
                                  "attractor_ic_scf",
                                  &string1,
                                  &flag1,
                                  errmsg),
                errmsg,
                errmsg);
    /* Complete set of parameters */
    if (flag1 == _TRUE_){
      if(string_begins_with(string1,'y') || string_begins_with(string1,'Y')){
        pba->attractor_ic_scf = _TRUE_;
      }
      else{
        pba->attractor_ic_scf = _FALSE_;
        /* Test */
        class_test(pba->scf_parameters_size<2,
                   errmsg,
                   "Since you are not using attractor initial conditions, you must specify phi and its derivative phi' as the last two entries in scf_parameters. See explanatory.ini for more details.");
        pba->phi_ini_scf = pba->scf_parameters[pba->scf_parameters_size-2];
        pba->phi_prime_ini_scf = pba->scf_parameters[pba->scf_parameters_size-1];
      }
    }

    /** 8.b.3) SCF tuning parameter */
    /* Read */
    class_read_int("scf_tuning_index",pba->scf_tuning_index);
    /* Test */
    class_test(pba->scf_tuning_index >= pba->scf_parameters_size,
               errmsg,
               "Tuning index 'scf_tuning_index' (%d) is larger than the number of entries (%d) in 'scf_parameters'.",
               pba->scf_tuning_index,pba->scf_parameters_size);

    /** 8.b.4) Shooting parameter */
    /* Read */
    class_read_double("scf_shooting_parameter",pba->scf_parameters[pba->scf_tuning_index]);
    /* Complete set of parameters */
    scf_lambda = pba->scf_parameters[0];
    if ((fabs(scf_lambda) < 3.)&&(pba->background_verbose>1)){
      printf("'scf_lambda' = %e < 3 won't be tracking (for exp quint) unless overwritten by tuning function.",scf_lambda);
    }
  }

  return _SUCCESS_;

}


/**
 * Read the parameters of injection structure
 * (These are all exotic processes of energy injection)
 *
 * @param pfc     Input: pointer to local structure
 * @param ppr     Input: pointer to precision structure
 * @param pth     Input: pointer to thermodynamics structure
 * @param errmsg  Input: Error message
 * @return the error status
 */

int input_read_parameters_injection(struct file_content * pfc,
                                    struct precision * ppr,
                                    struct thermodynamics * pth,
                                    ErrorMsg errmsg){

  /** Summary: */

  /** - Define local variables */
  struct injection* pin = &(pth->in);
  int flag1;
  char string1[_ARGUMENT_LENGTH_MAX_];
  string1[0]='\0';

  /** 1) DM annihilation */
  /** 1.a) Annihilation efficiency */
  /* Read */
  class_read_double("DM_annihilation_efficiency",pin->DM_annihilation_efficiency);
  class_read_double("DM_annihilation_cross_section",pin->DM_annihilation_cross_section);
  class_read_double("DM_annihilation_mass",pin->DM_annihilation_mass);
  class_read_double("DM_annihilation_fraction",pin->DM_annihilation_fraction);

  /* Test consistency of this input */
  class_test(pin->DM_annihilation_efficiency<0,
             errmsg,
             "annihilation efficiency cannot be negative");
  class_test(pin->DM_annihilation_efficiency>1.e-4,
             errmsg,
             "annihilation parameter suspiciously large (%e, while typical bounds are in the range of 1e-7 to 1e-6)",pin->DM_annihilation_efficiency);
  class_test(pin->DM_annihilation_mass<0. || pin->DM_annihilation_cross_section <0,
             errmsg,
             "Both mass and cross section for your dark matter particle must be positive.");
  class_test(pin->DM_annihilation_mass ==0 && pin->DM_annihilation_cross_section >0,
             errmsg,
             "you have annihilation_cross_section > 0 but DM_mass = 0. That is weird, please check your param file and set 'DM_mass' [GeV] to a non-zero value.\n");
  class_test((pin->DM_annihilation_cross_section !=0 || pin->DM_annihilation_mass !=0 || pin->DM_annihilation_fraction !=0) && pin->DM_annihilation_efficiency != 0,
             errmsg,
             "You can only enter one of {'DM_annihilation_cross_section', 'DM_annihilation_mass', 'DM_annihilation_fraction'} or 'annihilation_efficiency'.");
  if ((pin->DM_annihilation_efficiency >0) && (pth->reio_parametrization == reio_none) && (ppr->recfast_Heswitch >= 3) && (pth->recombination==recfast)) {
      printf("Warning: if you have DM annihilation and you use recfast with option recfast_Heswitch >= 3, then the expression for CfHe_t and dy[1] becomes undefined at late times, producing nan's. This is however masked by reionization if you are not in reio_none mode.");
  } //TODO :: check if still occurs !!!

  /* Complete set of parameters */
  if(pin->DM_annihilation_mass > 0 && pin->DM_annihilation_cross_section > 0.){
    pin->DM_annihilation_efficiency = pin->DM_annihilation_cross_section*1.e-6/(pin->DM_annihilation_mass*_eV_*1.e9);
  }

  if(pin->DM_annihilation_efficiency > 0){
    pth->has_exotic_injection = _TRUE_;
  }

  if (pin->DM_annihilation_efficiency > 0.) {
    /** 1.a.1) Model energy fraction absorbed by the gas as a function of redhsift */
    /* Read */
    class_read_double("DM_annihilation_variation",pin->DM_annihilation_variation);
    class_read_double("DM_annihilation_z",pin->DM_annihilation_z);
    class_read_double("DM_annihilation_zmax",pin->DM_annihilation_zmax);
    class_read_double("DM_annihilation_zmin",pin->DM_annihilation_zmin);
    class_read_double("DM_annihilation_f_halo",pin->DM_annihilation_f_halo);
    class_read_double("DM_annihilation_z_halo",pin->DM_annihilation_z_halo);
    /* Test */
    class_test(pin->DM_annihilation_variation>0,
               errmsg,
               "annihilation variation parameter must be negative (decreasing annihilation rate)");
    class_test(pin->DM_annihilation_z<0,
               errmsg,
               "characteristic annihilation redshift cannot be negative");
    class_test(pin->DM_annihilation_zmin<0,
               errmsg,
               "characteristic annihilation redshift cannot be negative");
    class_test(pin->DM_annihilation_zmax<0,
               errmsg,
               "characteristic annihilation redshift cannot be negative");
    class_test(pin->DM_annihilation_f_halo<0,
               errmsg,
               "Parameter for DM annihilation in halos cannot be negative");
    class_test(pin->DM_annihilation_z_halo<0,
               errmsg,
               "Parameter for DM annihilation in halos cannot be negative");
  }

  /** 2) DM decay */
  /** 2.a) Fraction */
  /* Read */
  class_read_double("DM_decay_fraction",pin->DM_decay_fraction);
  if(pin->DM_decay_fraction!=0){
    pth->has_exotic_injection = _TRUE_;
  }
  /* Test */
  class_test(pin->DM_decay_fraction<0,
             errmsg,
             "You need to enter a positive fraction of decaying DM. Please adjust your param file.");

  /** 2.b) Decay width */
  /* Read */
  class_read_double("DM_decay_Gamma",pin->DM_decay_Gamma);


  /** 3) PBH evaporation */
  /** 3.a) Fraction */
  /* Read */
  class_read_double("PBH_evaporation_fraction",pin->PBH_evaporation_fraction);
  if(pin->PBH_evaporation_fraction!=0){
    pth->has_exotic_injection = _TRUE_;
  }
  /* Test */
  class_test(pin->PBH_evaporation_fraction <0.,
             errmsg,
             "You need to enter a positive fraction of evaporating PBH. Please adjust your param file.");

  /** 3.b) Mass */
  /* Read */
  class_read_double("PBH_evaporation_mass",pin->PBH_evaporation_mass);
  /* Test */
  class_test(pin->PBH_evaporation_mass<0.,
             errmsg,
             "You need to enter a positive mass for your PBH.");
  class_test(pin->PBH_evaporation_mass>0. && pin->PBH_evaporation_fraction == 0,
             errmsg,
            "You have 'PBH_evaporation_mass > 0.' but 'PBH_evaporation_fraction = 0'. Please adjust your param file.");
  class_test(pin->PBH_evaporation_fraction>0. && pin->PBH_evaporation_mass == 0.,
             errmsg,
             "You have asked for a fraction of PBH being DM but you have zero mass. Please adjust your param file.");


  /** 4) PBH matter accretion */
  /** 4.a) Fraction */
  /* Read */
  class_read_double("PBH_accretion_fraction",pin->PBH_accretion_fraction);
  if(pin->PBH_accretion_fraction!=0){
    pth->has_exotic_injection = _TRUE_;
  }
  /* Test */
  class_test(pin->PBH_accretion_fraction < 0.,
             errmsg,
             "You need to enter a positive fraction of accreting PBH. Please adjust your param file.");

  /** 4.b) Mass */
  /* Read */
  class_read_double("PBH_accretion_mass",pin->PBH_accretion_mass);
  /* Test */
  class_test(pin->PBH_accretion_mass<0.,
             errmsg,
             "You need to enter a positive mass for your PBH.");
  class_test(pin->PBH_accretion_mass>0. && pin->PBH_accretion_fraction == 0,
             errmsg,
            "You have 'PBH_accretion_mass > 0.' but 'PBH_accretion_fraction = 0'. Please adjust your param file.");
  class_test(pin->PBH_accretion_fraction>0. && pin->PBH_accretion_mass == 0.,
             errmsg,
             "You have asked for a fraction of PBH being DM but you have zero mass. Please adjust your param file.");

  /** 4.c) Recipe */
  /* Read */
  class_call(parser_read_string(pfc,"PBH_accretion_recipe",&string1,&flag1,errmsg),
             errmsg,
             errmsg);
  /* Complete set of parameters */
  if (flag1 == _TRUE_){
    if (strcmp(string1,"spherical_accretion") == 0) {
      pin->PBH_accretion_recipe=spherical_accretion;
    }
    else if (strcmp(string1,"disk_accretion") == 0) {
      pin->PBH_accretion_recipe=disk_accretion;
    }
    else{
      class_stop(errmsg,
                 "You specified 'PBH_accretion_recipe' as '%s'. It has to be one of {'spherical_accretion','disk_accretion'}.",string1);
    }
  }

  /** 4.c.1) Additional parameters specific for spherical accretion */
  if(pin->PBH_accretion_recipe == spherical_accretion){
    /* Read */
    class_read_double("PBH_accretion_relative_velocities",pin->PBH_accretion_relative_velocities);
  }

  /** 4.c.2) Additional parameters specific for disk accretion */
  if(pin->PBH_accretion_recipe == disk_accretion){
    /* Read */
    class_read_double("PBH_accretion_ADAF_delta",pin->PBH_accretion_ADAF_delta);
    class_read_double("PBH_accretion_eigenvalue",pin->PBH_accretion_eigenvalue);
    /* Test */
    class_test(pin->PBH_accretion_ADAF_delta != 1e-3 && pin->PBH_accretion_ADAF_delta != 0.5  && pin->PBH_accretion_ADAF_delta != 0.1,
               errmsg,
               "The parameter 'pth->PBH_ADAF_delta' can currently only be set to 1e-3, 0.1 or 0.5.");
  }


  /** 5) Injection efficiency */
  /* Read */
  class_call(parser_read_string(pfc,"f_eff_type",&string1,&flag1,errmsg),
             errmsg,
             errmsg);
  /* Complete set of parameters */
  if (flag1 == _TRUE_){
    if (strcmp(string1,"on_the_spot") == 0){
      pin->f_eff_type = f_eff_on_the_spot;
      class_read_double("f_eff",pin->f_eff);
    }
    else if (strcmp(string1,"from_file") == 0){
      pin->f_eff_type = f_eff_from_file;
      /* Read */
      class_call(parser_read_string(pfc,"f_eff_file",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);
      /* Test */
      class_test(flag1 == _FALSE_,
                 errmsg,
                 "for the option 'from_file' for 'f_eff_type' the option 'f_eff_file' is required.");
      /* Complete set of parameters */
      strcpy(pin->f_eff_file, string1);
    }
    else{
      class_stop(errmsg,
                 "You specified 'f_eff_type' as '%s'. It has to be one of {'on_the_spot','from_file'}.",string1);
    }
  }

  /** 6) deposition function */
  /* Read */
  class_call(parser_read_string(pfc,"chi_type",&string1,&flag1,errmsg),
               errmsg,
               errmsg);
  /* Complete set of parameters */
  if (flag1 == _TRUE_){
    if (strcmp(string1,"CK_2004") == 0){
      pin->chi_type = chi_CK;
    }
    else if (strcmp(string1,"PF_2005") == 0){
      pin->chi_type = chi_PF;
    }
    else if (strcmp(string1,"Galli_2013_file") == 0){
      pin->chi_type = chi_Galli_file;
    }
    else if (strcmp(string1,"Galli_2013_analytic") == 0){
      pin->chi_type = chi_Galli_analytic;
    }
    else if (strcmp(string1,"heat") == 0){
      pin->chi_type = chi_full_heating;
    }
    else if (strcmp(string1,"from_x_file") == 0){
      pin->chi_type = chi_from_x_file;
    }
    else if (strcmp(string1,"from_z_file") == 0){
      pin->chi_type = chi_from_z_file;
    }
    else{
      class_stop(errmsg,
                   "You specified 'chi_type' as '%s'. It has to be one of {'CK_2004','PF_2005','Galli_2013_file','Galli_2013_analytic','heat','from_x_file','from_z_file'}.",string1);
    }
  }

  if(pin->chi_type == chi_from_x_file || pin->chi_type == chi_from_z_file){
    /** 6.a) External file */
    /* Read */
    class_call(parser_read_string(pfc,"chi_file",&string1,&flag1,errmsg),
               errmsg,
               errmsg);
    /* Test */
    class_test(flag1 == _FALSE_,
               errmsg,
               "for the option 'from_x_file' or 'from_z_file' for 'chi_type' the option 'chi_file' is required.");
    /* Complete set of parameters */
    strcpy(pin->chi_z_file, string1);
    strcpy(pin->chi_x_file, string1);
  }

  return _SUCCESS_;

}


/**
 * Read the parameters of fourier structure.
 *
 * @param pfc            Input: pointer to local structure
 * @param ppr            Input: pointer to precision structure
 * @param pba            Input: pointer to background structure
 * @param pth            Input: pointer to thermodynamics structure
 * @param ppt            Input: pointer to perturbations structure
 * @param pfo            Input: pointer to fourier structure
 * @param input_verbose  Input: verbosity of input
 * @param errmsg         Input: Error message
 * @return the error status
 */

int input_read_parameters_nonlinear(struct file_content * pfc,
                                    struct precision * ppr,
                                    struct background * pba,
                                    struct thermodynamics * pth,
                                    struct perturbations * ppt,
                                    struct fourier * pfo,
                                    int input_verbose,
                                    ErrorMsg errmsg){

  /** Define local variables */
  int flag1,flag2,flag3;
  double param2,param3;
  char string1[_ARGUMENT_LENGTH_MAX_];

  /** 1) Non-linearity */
  /* Read */
  class_call(parser_read_string(pfc,"non_linear",&string1,&flag1,errmsg),
             errmsg,
             errmsg);
  /* Compatibility code BEGIN */
  if(flag1 == _FALSE_){
    class_call(parser_read_string(pfc,"non linear",&string1,&flag1,errmsg),
               errmsg,
               errmsg);
  }
  /* Compatibility code END */

  if (flag1 == _TRUE_) {
    /* Test */
    class_test(ppt->has_perturbations == _FALSE_,
               errmsg,
               "You requested non-linear computation but no perturbations. You must set the 'output' field.");
    /* Complete set of parameters */
    if ((strstr(string1,"halofit") != NULL) || (strstr(string1,"Halofit") != NULL) || (strstr(string1,"HALOFIT") != NULL)) {
      pfo->method=nl_halofit;
      ppt->has_nl_corrections_based_on_delta_m = _TRUE_;
    }
    else if((strstr(string1,"hmcode") != NULL) || (strstr(string1,"HMCODE") != NULL) || (strstr(string1,"HMcode") != NULL) || (strstr(string1,"Hmcode") != NULL)) {
      pfo->method=nl_HMcode;
      ppt->k_max_for_pk = MAX(ppt->k_max_for_pk,MAX(ppr->hmcode_min_k_max,ppr->fourier_min_k_max));
      ppt->has_nl_corrections_based_on_delta_m = _TRUE_;
      class_read_int("extrapolation_method",pfo->extrapolation_method);

      class_call(parser_read_string(pfc,
                                    "feedback model",
                                    &(string1),
                                    &(flag1),
                                    errmsg),
                 errmsg,
                 errmsg);

      if (flag1 == _TRUE_) {

        if (strstr(string1,"emu_dmonly") != NULL) {
          pfo->feedback = nl_emu_dmonly;
        }
        if (strstr(string1,"owls_dmonly") != NULL) {
          pfo->feedback = nl_owls_dmonly;
        }
        if (strstr(string1,"owls_ref") != NULL) {
          pfo->feedback = nl_owls_ref;
        }
        if (strstr(string1,"owls_agn") != NULL) {
          pfo->feedback = nl_owls_agn;
        }
        if (strstr(string1,"owls_dblim") != NULL) {
          pfo->feedback = nl_owls_dblim;
        }
      }

      class_call(parser_read_double(pfc,"eta_0",&param2,&flag2,errmsg),
                 errmsg,
                 errmsg);
      class_call(parser_read_double(pfc,"c_min",&param3,&flag3,errmsg),
                 errmsg,
                 errmsg);

      class_test(((flag1 == _TRUE_) && ((flag2 == _TRUE_) || (flag3 == _TRUE_))),
                 errmsg,
                 "In input file, you cannot enter both a baryonic feedback model and a choice of baryonic feedback parameters, choose one of both methods");

      if ((flag2 == _TRUE_) && (flag3 == _TRUE_)) {
        pfo->feedback = nl_user_defined;
        class_read_double("eta_0", pfo->eta_0);
        class_read_double("c_min", pfo->c_min);
      }
      else if ((flag2 == _TRUE_) && (flag3 == _FALSE_)) {
        pfo->feedback = nl_user_defined;
        class_read_double("eta_0", pfo->eta_0);
        pfo->c_min = (0.98 - pfo->eta_0)/0.12;
      }
      else if ((flag2 == _FALSE_) && (flag3 == _TRUE_)) {
        pfo->feedback = nl_user_defined;
        class_read_double("c_min", pfo->c_min);
        pfo->eta_0 = 0.98 - 0.12*pfo->c_min;
      }

      class_read_double("z_infinity", pfo->z_infinity);
    }
    else if(strstr(string1,"no")!=NULL){
      pfo->method=nl_none;
      ppt->has_nl_corrections_based_on_delta_m = _FALSE_;
    }
    else{
      class_stop(errmsg,
                 "You specified 'non_linear' = '%s'. It has to be one of {'halofit','hmcode','none'}.",string1);
    }
  }

  /** - special steps if we want Halofit with wa_fld non-zero:
      so-called "Pk_equal method" of 0810.0190 and 1601.07230 */

  if (pfo->method == nl_halofit) {

    class_call(parser_read_string(pfc,"pk_eq",&string1,&flag1,errmsg),
               errmsg,
               errmsg);

    if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL))) {

      if ((pba->Omega0_fld != 0.) && (pba->wa_fld != 0.)){

        pfo->has_pk_eq = _TRUE_;
      }
    }
  }

  return _SUCCESS_;
}


/**
 * Perform preliminary steps fur using the method called Pk_equal,
 * described in 0810.0190 and 1601.07230, extending the range of
 * validity of HALOFIT from constant w to (w0,wa) models. In that
 * case, one must compute here some effective values of w0_eff(z_i)
 * and Omega_m_eff(z_i), that will be interpolated later at arbitrary
 * redshift in the non-linear module.
 *
 * Returns table of values [z_i, tau_i, w0_eff_i, Omega_m_eff_i]
 * stored in fourier structure.
 *
 * @param ppr           Input: pointer to precision structure
 * @param pba           Input: pointer to background structure
 * @param pth           Input: pointer to thermodynamics structure
 * @param pfo           Input/Output: pointer to fourier structure
 * @param input_verbose Input: verbosity of this input module
 * @param errmsg        Input/Ouput: error message
 * @return the error status
 */

int input_prepare_pk_eq(struct precision * ppr,
                        struct background *pba,
                        struct thermodynamics *pth,
                        struct fourier *pfo,
                        int input_verbose,
                        ErrorMsg errmsg) {

  /** Summary: */

  /** Define local variables */
  double tau_of_z;
  double delta_tau;
  double error;
  double delta_tau_eq;
  double * pvecback;
  int last_index=0;
  int index_pk_eq_z;
  int index_eq;
  int true_background_verbose;
  int true_thermodynamics_verbose;
  int true_hyrec_verbose;
  double true_w0_fld;
  double true_wa_fld;
  double * z;

  /** Store the true cosmological parameters (w0, wa) somwhere before using temporarily some fake ones in this function */
  true_background_verbose = pba->background_verbose;
  true_thermodynamics_verbose = pth->thermodynamics_verbose;
  true_hyrec_verbose = pth->hyrec_verbose;
  true_w0_fld = pba->w0_fld;
  true_wa_fld = pba->wa_fld;

  /** The fake calls of the background and thermodynamics module will be done in non-verbose mode */
  pba->background_verbose = 0;
  pth->thermodynamics_verbose = 0;
  pth->hyrec_verbose = 0;

  /** Allocate indices and arrays for storing the results */
  pfo->pk_eq_tau_size = ppr->pk_eq_Nzlog;
  class_alloc(pfo->pk_eq_tau,
              pfo->pk_eq_tau_size*sizeof(double),
              errmsg);
  class_alloc(z,
              pfo->pk_eq_tau_size*sizeof(double),
              errmsg);

  index_eq = 0;
  class_define_index(pfo->index_pk_eq_w,_TRUE_,index_eq,1);
  class_define_index(pfo->index_pk_eq_Omega_m,_TRUE_,index_eq,1);
  pfo->pk_eq_size = index_eq;
  class_alloc(pfo->pk_eq_w_and_Omega,
              pfo->pk_eq_tau_size*pfo->pk_eq_size*sizeof(double),
              errmsg);
  class_alloc(pfo->pk_eq_ddw_and_ddOmega,
              pfo->pk_eq_tau_size*pfo->pk_eq_size*sizeof(double),
              errmsg);

  /** Call the background module in order to fill a table of tau_i[z_i] */
  class_call(background_init(ppr,pba), pba->error_message, errmsg);
  for (index_pk_eq_z=0; index_pk_eq_z<pfo->pk_eq_tau_size; index_pk_eq_z++) {
    z[index_pk_eq_z] = exp(log(1.+ppr->pk_eq_z_max)/(pfo->pk_eq_tau_size-1)*index_pk_eq_z)-1.;
    class_call(background_tau_of_z(pba,
                                   z[index_pk_eq_z],
                                   &tau_of_z),
               pba->error_message,
               errmsg);
    pfo->pk_eq_tau[index_pk_eq_z] = tau_of_z;
  }
  class_call(background_free_noinput(pba),
             pba->error_message,
             errmsg);

  /** Loop over z_i values. For each of them, we will call the
     background and thermodynamics module for fake models. The goal is
     to find, for each z_i, and effective w0_eff[z_i] and
     Omega_m_eff{z_i], such that: the true model with (w0,wa) and the
     equivalent model with (w0_eff[z_i],0) have the same conformal
     distance between z_i and z_recombination, namely chi = tau[z_i] -
     tau_rec. It is thus necessary to call both the background and
     thermodynamics module for each fake model and to re-compute
     tau_rec for each of them. Once the eqauivalent model is found we
     compute and store Omega_m_effa(z_i) of the equivalent model */
  for (index_pk_eq_z=0; index_pk_eq_z<pfo->pk_eq_tau_size; index_pk_eq_z++) {

    if (input_verbose > 2)
      printf("    * computing Pk_equal parameters at z=%e\n",z[index_pk_eq_z]);
    /* get chi = (tau[z_i] - tau_rec) in true model */
    pba->w0_fld = true_w0_fld;
    pba->wa_fld = true_wa_fld;
    class_call(background_init(ppr,pba),
               pba->error_message,
               errmsg);
    class_call(thermodynamics_init(ppr,pba,pth),
               pth->error_message,
               errmsg);
    delta_tau = pfo->pk_eq_tau[index_pk_eq_z] - pth->tau_rec;
    /* launch iterations in order to coverge to effective model with wa=0 but the same chi = (tau[z_i] - tau_rec) */
    pba->wa_fld=0.;

    do {
      class_call(background_free_noinput(pba),
                 pba->error_message,
                 errmsg);
      class_call(thermodynamics_free(pth),
                 pth->error_message,
                 errmsg);

      class_call(background_init(ppr,pba),
                 pba->error_message,
                 errmsg);
      class_call(background_tau_of_z(pba,
                                     z[index_pk_eq_z],
                                     &tau_of_z),
                 pba->error_message,
                 errmsg);
      class_call(thermodynamics_init(ppr,pba,pth),
                 pth->error_message,
                 errmsg);

      delta_tau_eq = tau_of_z - pth->tau_rec;

      error = 1.-delta_tau_eq/delta_tau;
      pba->w0_fld = pba->w0_fld*pow(1.+error,10.);

    }
    while(fabs(error) > ppr->pk_eq_tol);

    /* Equivalent model found. Store w0(z) in that model. Find Omega_m(z) in that model and store it. */
    pfo->pk_eq_w_and_Omega[pfo->pk_eq_size*index_pk_eq_z+pfo->index_pk_eq_w] = pba->w0_fld;
    class_alloc(pvecback,
                pba->bg_size*sizeof(double),
                pba->error_message);
    class_call(background_at_tau(pba,
                                 tau_of_z,
                                 long_info,
                                 inter_normal,
                                 &last_index,
                                 pvecback),
               pba->error_message, errmsg);
    pfo->pk_eq_w_and_Omega[pfo->pk_eq_size*index_pk_eq_z+pfo->index_pk_eq_Omega_m] = pvecback[pba->index_bg_Omega_m];
    free(pvecback);

    class_call(background_free_noinput(pba),
               pba->error_message,
               errmsg);
    class_call(thermodynamics_free(pth),
               pth->error_message,
               errmsg);
  }

  /** Restore cosmological parameters (w0, wa) to their true values before main call to CLASS modules */
  pba->background_verbose = true_background_verbose;
  pth->thermodynamics_verbose = true_thermodynamics_verbose;
  pth->hyrec_verbose = true_hyrec_verbose;

  pba->w0_fld = true_w0_fld;
  pba->wa_fld = true_wa_fld;

  /* in verbose mode, report the results */

  if (input_verbose > 1) {
    fprintf(stdout,"    Effective parameters for Pk_equal:\n");

    for (index_pk_eq_z=0; index_pk_eq_z<pfo->pk_eq_tau_size; index_pk_eq_z++) {
      fprintf(stdout,"    * at z=%e, tau=%e w=%e Omega_m=%e\n",
              z[index_pk_eq_z],
              pfo->pk_eq_tau[index_pk_eq_z],
              pfo->pk_eq_w_and_Omega[pfo->pk_eq_size*index_pk_eq_z+pfo->index_pk_eq_w],
              pfo->pk_eq_w_and_Omega[pfo->pk_eq_size*index_pk_eq_z+pfo->index_pk_eq_Omega_m]);
    }
  }
  free(z);

  /** Spline the table for later interpolation */
  class_call(array_spline_table_lines(pfo->pk_eq_tau,
                                      pfo->pk_eq_tau_size,
                                      pfo->pk_eq_w_and_Omega,
                                      pfo->pk_eq_size,
                                      pfo->pk_eq_ddw_and_ddOmega,
                                      _SPLINE_NATURAL_,
                                      errmsg),
             errmsg,
             errmsg);

  return _SUCCESS_;

}


/**
 * Read the parameters of primordial structure.
 *
 * @param pfc     Input: pointer to local structure
 * @param ppt     Input: pointer to perturbations structure
 * @param ppm     Input: pointer to primordial structure
 * @param errmsg  Input: Error message
 * @return the error status
 */

int input_read_parameters_primordial(struct file_content * pfc,
                                     struct perturbations * ppt,
                                     struct primordial * ppm,
                                     ErrorMsg errmsg){

  /** Summary: */

  /** Define local variables */
  int flag1, flag2;
  double param1, param2;
  char string1[_ARGUMENT_LENGTH_MAX_];
  char string2[_ARGUMENT_LENGTH_MAX_];
  double R0,R1,R2,R3,R4;
  double PSR0,PSR1,PSR2,PSR3,PSR4;
  double HSR0,HSR1,HSR2,HSR3,HSR4;
  double k1=0.;
  double k2=0.;
  double prr1=0.;
  double prr2=0.;
  double pii1=0.;
  double pii2=0.;
  double pri1=0.;
  double pri2=0.;
  double n_iso=0.;
  double f_iso=0.;
  double n_cor=0.;
  double c_cor=0.;

  /** 1) Primordial spectrum type */
  /* Read */
  class_call(parser_read_string(pfc,"Pk_ini_type",&string1,&flag1,errmsg),
             errmsg,
             errmsg);
  /* Compatibility code BEGIN */
  if(flag1 == _FALSE_){
    class_call(parser_read_string(pfc,"P_k_ini type",&string1,&flag1,errmsg),
               errmsg,
               errmsg);
  }
  /* Compatibility code END */
  /* Complete set of parameters */
  if (flag1 == _TRUE_) {
    if (strcmp(string1,"analytic_Pk") == 0){
      ppm->primordial_spec_type = analytic_Pk;
    }
    else if (strcmp(string1,"inflation_V") == 0){
      ppm->primordial_spec_type = inflation_V;
    }
    else if (strcmp(string1,"inflation_H") == 0){
      ppm->primordial_spec_type = inflation_H;
    }
    else if (strcmp(string1,"inflation_V_end") == 0){
      ppm->primordial_spec_type = inflation_V_end;
    }
    else if (strcmp(string1,"two_scales") == 0){
      ppm->primordial_spec_type = two_scales;
    }
    else if (strcmp(string1,"external_Pk") == 0){
      ppm->primordial_spec_type = external_Pk;
    }
    else{
      class_stop(errmsg,
                 "You specified 'P_k_ini_type' as '%s'. It has to be one of {'analytic_Pk','inflation_V','inflation_V_end','two_scales','external_Pk'}.",string1);
    }
  }

  /** 1.a) Pivot scale in Mpc-1 */
  /* Read */
  class_read_double("k_pivot",ppm->k_pivot);

  /** 1.b) For type 'analytic_Pk' */
  if (ppm->primordial_spec_type == analytic_Pk) {

    /** 1.b.1) For scalar perturbations */
    if (ppt->has_scalars == _TRUE_) {
      /* Read */
      class_call(parser_read_double(pfc,"A_s",&param1,&flag1,errmsg),
                 errmsg,
                 errmsg);
      class_call(parser_read_double(pfc,"ln10^{10}A_s",&param2,&flag2,errmsg),
                 errmsg,
                 errmsg);
      /* Test */
      class_test((flag1 == _TRUE_) && (flag2 == _TRUE_),
                 errmsg,
                 "You can only enter one of 'A_s' or 'ln10^{10}A_s'.");
      /* Complete set of parameters */
      if (flag1 == _TRUE_){
        ppm->A_s = param1;
      }
      else if (flag2 == _TRUE_){
        ppm->A_s = exp(param2)*1.e-10;
      }

      /** 1.b.1.1) Adiabatic perturbations */
      if (ppt->has_ad == _TRUE_) {
        /* Read */
        class_read_double("n_s",ppm->n_s);
        class_read_double("alpha_s",ppm->alpha_s);
      }

      /** 1.b.1.2) Isocurvature/entropy perturbations */
      /* Read */
      if (ppt->has_bi == _TRUE_) {
        class_read_double("f_bi",ppm->f_bi);
        class_read_double("n_bi",ppm->n_bi);
        class_read_double("alpha_bi",ppm->alpha_bi);
      }
      if (ppt->has_cdi == _TRUE_) {
        class_read_double("f_cdi",ppm->f_cdi);
        class_read_double("n_cdi",ppm->n_cdi);
        class_read_double("alpha_cdi",ppm->alpha_cdi);
      }
      if (ppt->has_nid == _TRUE_) {
        class_read_double("f_nid",ppm->f_nid);
        class_read_double("n_nid",ppm->n_nid);
        class_read_double("alpha_nid",ppm->alpha_nid);
      }
      if (ppt->has_niv == _TRUE_) {
        class_read_double("f_niv",ppm->f_niv);
        class_read_double("n_niv",ppm->n_niv);
        class_read_double("alpha_niv",ppm->alpha_niv);
      }

      /** 1.b.1.3) Cross-correlation between different adiabatic/entropy mode */
      /* Read */
      if ((ppt->has_ad == _TRUE_) && (ppt->has_bi == _TRUE_)) {
        class_read_double_one_of_two("c_ad_bi","c_bi_ad",ppm->c_ad_bi);
        class_read_double_one_of_two("n_ad_bi","n_bi_ad",ppm->n_ad_bi);
        class_read_double_one_of_two("alpha_ad_bi","alpha_bi_ad",ppm->alpha_ad_bi);
      }
      if ((ppt->has_ad == _TRUE_) && (ppt->has_cdi == _TRUE_)) {
        class_read_double_one_of_two("c_ad_cdi","c_cdi_ad",ppm->c_ad_cdi);
        class_read_double_one_of_two("n_ad_cdi","n_cdi_ad",ppm->n_ad_cdi);
        class_read_double_one_of_two("alpha_ad_cdi","alpha_cdi_ad",ppm->alpha_ad_cdi);
      }
      if ((ppt->has_ad == _TRUE_) && (ppt->has_nid == _TRUE_)) {
        class_read_double_one_of_two("c_ad_nid","c_nid_ad",ppm->c_ad_nid);
        class_read_double_one_of_two("n_ad_nid","n_nid_ad",ppm->n_ad_nid);
        class_read_double_one_of_two("alpha_ad_nid","alpha_nid_ad",ppm->alpha_ad_nid);
      }
      if ((ppt->has_ad == _TRUE_) && (ppt->has_niv == _TRUE_)) {
        class_read_double_one_of_two("c_ad_niv","c_niv_ad",ppm->c_ad_niv);
        class_read_double_one_of_two("n_ad_niv","n_niv_ad",ppm->n_ad_niv);
        class_read_double_one_of_two("alpha_ad_niv","alpha_niv_ad",ppm->alpha_ad_niv);
      }
      if ((ppt->has_bi == _TRUE_) && (ppt->has_cdi == _TRUE_)) {
        class_read_double_one_of_two("c_bi_cdi","c_cdi_bi",ppm->c_bi_cdi);
        class_read_double_one_of_two("n_bi_cdi","n_cdi_bi",ppm->n_bi_cdi);
        class_read_double_one_of_two("alpha_bi_cdi","alpha_cdi_bi",ppm->alpha_bi_cdi);
      }
      if ((ppt->has_bi == _TRUE_) && (ppt->has_nid == _TRUE_)) {
        class_read_double_one_of_two("c_bi_nid","c_nid_bi",ppm->c_bi_nid);
        class_read_double_one_of_two("n_bi_nid","n_nid_bi",ppm->n_bi_nid);
        class_read_double_one_of_two("alpha_bi_nid","alpha_nid_bi",ppm->alpha_bi_nid);
      }
      if ((ppt->has_bi == _TRUE_) && (ppt->has_niv == _TRUE_)) {
        class_read_double_one_of_two("c_bi_niv","c_niv_bi",ppm->c_bi_niv);
        class_read_double_one_of_two("n_bi_niv","n_niv_bi",ppm->n_bi_niv);
        class_read_double_one_of_two("alpha_bi_niv","alpha_niv_bi",ppm->alpha_bi_niv);
      }
      if ((ppt->has_cdi == _TRUE_) && (ppt->has_nid == _TRUE_)) {
        class_read_double_one_of_two("c_cdi_nid","c_nid_cdi",ppm->c_cdi_nid);
        class_read_double_one_of_two("n_cdi_nid","n_nid_cdi",ppm->n_cdi_nid);
        class_read_double_one_of_two("alpha_cdi_nid","alpha_nid_cdi",ppm->alpha_cdi_nid);
      }
      if ((ppt->has_cdi == _TRUE_) && (ppt->has_niv == _TRUE_)) {
        class_read_double_one_of_two("c_cdi_niv","c_niv_cdi",ppm->c_cdi_niv);
        class_read_double_one_of_two("n_cdi_niv","n_niv_cdi",ppm->n_cdi_niv);
        class_read_double_one_of_two("alpha_cdi_niv","alpha_niv_cdi",ppm->alpha_cdi_niv);
      }
      if ((ppt->has_nid == _TRUE_) && (ppt->has_niv == _TRUE_)) {
        class_read_double_one_of_two("c_nid_niv","c_niv_nid",ppm->c_nid_niv);
        class_read_double_one_of_two("n_nid_niv","n_niv_nid",ppm->n_nid_niv);
        class_read_double_one_of_two("alpha_nid_niv","alpha_niv_nid",ppm->alpha_nid_niv);
      }
    }

    /** 1.b.2) For tensor perturbations */
    if (ppt->has_tensors == _TRUE_){
      /* Read */
      class_read_double("r",ppm->r);
      if (ppt->has_scalars == _FALSE_){
        class_read_double("A_s",ppm->A_s);
      }
      if (ppm->r <= 0) {
        ppt->has_tensors = _FALSE_;
      }
      else {
        class_call(parser_read_string(pfc,"n_t",&string1,&flag1,errmsg),
                   errmsg,
                   errmsg);
        class_call(parser_read_string(pfc,"alpha_t",&string2,&flag2,errmsg),
                   errmsg,
                   errmsg);
        /* Complete set of parameters */
        if ((flag1 == _TRUE_) && !((strstr(string1,"SCC") != NULL) || (strstr(string1,"scc") != NULL))){
          class_read_double("n_t",ppm->n_t);
        }
        else {
          ppm->n_t = -ppm->r/8.*(2.-ppm->r/8.-ppm->n_s);        // enforce single slow-roll self-consistency condition (order 2 in slow-roll)
        }
        if ((flag2 == _TRUE_) && !((strstr(string2,"SCC") != NULL) || (strstr(string2,"scc") != NULL))) {
          class_read_double("alpha_t",ppm->alpha_t);
        }
        else {
          ppm->alpha_t = ppm->r/8.*(ppm->r/8.+ppm->n_s-1.);     // enforce single slow-roll self-consistency condition (order 2 in slow-roll)
        }
      }
    }
  }

  else if ((ppm->primordial_spec_type == inflation_V) || (ppm->primordial_spec_type == inflation_H)) {

    /** 1.c) For type 'inflation_V' */
    if (ppm->primordial_spec_type == inflation_V) {

      /** 1.c.1) Type of potential */
      /* Read */
      class_call(parser_read_string(pfc,"potential",&string1,&flag1,errmsg),    // only polynomial coded so far:
                 errmsg,                                                        // no need to interpret string1
                 errmsg);

      /** 1.c.2) Coefficients of the Taylor expansion */
      /* Read */
      class_call(parser_read_string(pfc,"PSR_0",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);
      if (flag1 == _TRUE_){
        PSR0=0.;
        PSR1=0.;
        PSR2=0.;
        PSR3=0.;
        PSR4=0.;
        class_read_double("PSR_0",PSR0);
        class_read_double("PSR_1",PSR1);
        class_read_double("PSR_2",PSR2);
        class_read_double("PSR_3",PSR3);
        class_read_double("PSR_4",PSR4);
        /* Test */
        class_test(PSR0 <= 0.,
                   errmsg,
                   "inconsistent parametrization of polynomial inflation potential");
        class_test(PSR1 <= 0.,
                   errmsg,
                   "inconsistent parametrization of polynomial inflation potential");
        /* Complete set of parameters */
        R0 = PSR0;
        R1 = PSR1*16.*_PI_;
        R2 = PSR2*8.*_PI_;
        R3 = PSR3*pow(8.*_PI_,2);
        R4 = PSR4*pow(8.*_PI_,3);
        ppm->V0 = R0*R1*3./128./_PI_;
        ppm->V1 = -sqrt(R1)*ppm->V0;
        ppm->V2 = R2*ppm->V0;
        ppm->V3 = R3*ppm->V0*ppm->V0/ppm->V1;
        ppm->V4 = R4*ppm->V0/R1;
      }
      else {
        /* Read */
        class_call(parser_read_string(pfc,"R_0",&string1,&flag1,errmsg),
                   errmsg,
                   errmsg);
        if (flag1 == _TRUE_) {
          R0=0.;
          R1=0.;
          R2=0.;
          R3=0.;
          R4=0.;
          class_read_double("R_0",R0);
          class_read_double("R_1",R1);
          class_read_double("R_2",R2);
          class_read_double("R_3",R3);
          class_read_double("R_4",R4);
          /* Test */
          class_test(R0 <= 0.,
                     errmsg,
                     "inconsistent parametrization of polynomial inflation potential");
          class_test(R1 <= 0.,
                     errmsg,
                     "inconsistent parametrization of polynomial inflation potential");
          /* Complete set of parameters */
          ppm->V0 = R0*R1*3./128./_PI_;
          ppm->V1 = -sqrt(R1)*ppm->V0;
          ppm->V2 = R2*ppm->V0;
          ppm->V3 = R3*ppm->V0*ppm->V0/ppm->V1;
          ppm->V4 = R4*ppm->V0/R1;
        }
        else {
          /* Read */
          class_read_double("V_0",ppm->V0);
          class_read_double("V_1",ppm->V1);
          class_read_double("V_2",ppm->V2);
          class_read_double("V_3",ppm->V3);
          class_read_double("V_4",ppm->V4);
        }
      }
    }

    /** 1.d) For type 'inflation_H' */
    else {
      /* Read */
      class_call(parser_read_string(pfc,"HSR_0",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);
      if (flag1 == _TRUE_) {
        HSR0=0.;
        HSR1=0.;
        HSR2=0.;
        HSR3=0.;
        HSR4=0.;
        class_read_double("HSR_0",HSR0);
        class_read_double("HSR_1",HSR1);
        class_read_double("HSR_2",HSR2);
        class_read_double("HSR_3",HSR3);
        class_read_double("HSR_4",HSR4);
        /* Complete set of parameters */
        ppm->H0 = sqrt(HSR0*HSR1*_PI_);
        ppm->H1 = -sqrt(4.*_PI_*HSR1)*ppm->H0;
        ppm->H2 = 4.*_PI_*HSR2*ppm->H0;
        ppm->H3 = 4.*_PI_*HSR3*ppm->H0*ppm->H0/ppm->H1;
        ppm->H4 = 4.*_PI_*HSR4*ppm->H0*ppm->H0*ppm->H0/ppm->H1/ppm->H1;

      }
      else {
        /* Read */
        class_read_double("H_0",ppm->H0);
        class_read_double("H_1",ppm->H1);
        class_read_double("H_2",ppm->H2);
        class_read_double("H_3",ppm->H3);
        class_read_double("H_4",ppm->H4);
      }
      /* Test */
      class_test(ppm->H0 <= 0.,
                 errmsg,
                 "inconsistent parametrization of polynomial inflation potential");
    }
  }

  /** 1.e) For type 'inflation_V_end' */
  else if (ppm->primordial_spec_type == inflation_V_end) {

    /** 1.e.1) Value of the field at the minimum of the potential */
    /* Read */
    class_read_double("phi_end",ppm->phi_end);

    /** 1.e.2) Shape of the potential */
    /* Read */
    class_call(parser_read_string(pfc,"full_potential",&string1,&flag1,errmsg),
               errmsg,
               errmsg);
    /* Complete set of parameters */
    if (flag1 == _TRUE_) {
      if (strcmp(string1,"polynomial") == 0){
        ppm->potential = polynomial;
      }
      else if (strcmp(string1,"higgs_inflation") == 0){
        ppm->potential = higgs_inflation;
      }
      else{
        class_stop(errmsg,"You specified 'full_potential' as '%s'. It has to be one of {'polynomial','higgs_inflation'}.",string1);
      }
    }

    /** 1.e.3) Parameters of the potential */
    /* Read */
    class_read_double("Vparam0",ppm->V0);
    class_read_double("Vparam1",ppm->V1);
    class_read_double("Vparam2",ppm->V2);
    class_read_double("Vparam3",ppm->V3);
    class_read_double("Vparam4",ppm->V4);

    /** 1.e.4) How much the scale factor a or the product (aH) increases between
               Hubble crossing for the pivot scale (during inflation) and the
               end of inflation */
    /* Read */
    class_call(parser_read_string(pfc,"ln_aH_ratio",&string1,&flag1,errmsg),
               errmsg,
               errmsg);
    class_call(parser_read_string(pfc,"N_star",&string2,&flag2,errmsg),
               errmsg,
               errmsg);
    /* Test */
    class_test((flag1 == _TRUE_) && (flag2 == _TRUE_),
               errmsg,
               "You can only enter one of 'ln_aH_ratio' or 'N_star'.");
    /* Complete set of parameters */
    if (flag1 == _TRUE_) {
      if ((strstr(string1,"auto") != NULL) || (strstr(string1,"AUTO") != NULL)){
        ppm->phi_pivot_method = ln_aH_ratio_auto;
      }
      else {
        ppm->phi_pivot_method = ln_aH_ratio;
        class_read_double("ln_aH_ratio",ppm->phi_pivot_target);
      }
    }
    if (flag2 == _TRUE_) {
      ppm->phi_pivot_method = N_star;
      class_read_double("N_star",ppm->phi_pivot_target);
    }

     /** 1.e.5) Should the inflation module do its nomral job of numerical
                integration ('numerical') or use analytical slow-roll formulas
                to infer the primordial spectrum from the potential
                ('analytical')? */
    /* Read */
    class_call(parser_read_string(pfc,"inflation_behavior",&string1,&flag1,errmsg),
               errmsg,
               errmsg);
    /* Complete set of parameters */
    if (flag1 == _TRUE_) {
      if (strstr(string1,"numerical") != NULL){
        ppm->behavior = numerical;
      }
      else if (strstr(string1,"analytical") != NULL){
        ppm->behavior = analytical;
      }
      else{
        class_stop(errmsg,"You specified 'inflation_behavior' as '%s'. It has to be one of {'numerical','analytical'}.",string1);
      }
    }
  }

  /** 1.f) For type 'two_scales' */
  else if (ppm->primordial_spec_type == two_scales) {

    /** 1.f.1) Wavenumbers */
    /* Read */
    class_read_double("k1",k1);
    class_read_double("k2",k2);
    /* Test */
    class_test(k1<=0.,errmsg,"enter strictly positive scale k1");
    class_test(k2<=0.,errmsg,"enter strictly positive scale k2");

    if (ppt->has_scalars == _TRUE_){

      /** 1.f.2) Amplitudes for the adiabatic primordial spectrum */
      /* Read */
      class_read_double("P_{RR}^1",prr1);
      class_read_double("P_{RR}^2",prr2);
      /* Test */
      class_test(prr1<=0.,errmsg,"enter strictly positive scale P_{RR}^1");
      class_test(prr2<=0.,errmsg,"enter strictly positive scale P_{RR}^2");
      /* Complete set of parameters */
      ppm->n_s = log(prr2/prr1)/log(k2/k1)+1.;
      ppm->A_s = prr1*exp((ppm->n_s-1.)*log(ppm->k_pivot/k1));

      /** 1.f.3) Isocurvature amplitudes */
      if ((ppt->has_bi == _TRUE_) || (ppt->has_cdi == _TRUE_) || (ppt->has_nid == _TRUE_) || (ppt->has_niv == _TRUE_)){
        /* Read */
        class_read_double("P_{II}^1",pii1);
        class_read_double("P_{II}^2",pii2);
        class_read_double("P_{RI}^1",pri1);
        class_read_double("|P_{RI}^2|",pri2);
        /* Test */
        class_test(pii1 <= 0.,
                   errmsg,
                   "since you request iso modes, you should have P_{ii}^1 strictly positive");
        class_test(pii2 < 0.,
                   errmsg,
                   "since you request iso modes, you should have P_{ii}^2 positive or eventually null");
        class_test(pri2 < 0.,
                   errmsg,
                   "by definition, you should have |P_{ri}^2| positive or eventually null");


        /** 1.f.4) Uncorrelated or anti-correlated? */
        /* Read */
        class_call(parser_read_string(pfc,"special iso",&string1,&flag1,errmsg),
                   errmsg,
                   errmsg);
        /* Complete set of parameters */
        if ((flag1 == _TRUE_) && (strstr(string1,"axion") != NULL)){    // Axion case, only one iso parameter: piir1
          n_iso = 1.;
          n_cor = 0.;
          c_cor = 0.;
        }
        else if ((flag1 == _TRUE_) && (strstr(string1,"anticurvaton") != NULL)){      // Curvaton case, only one iso parameter: piir1
          n_iso = ppm->n_s;
          n_cor = 0.;
          c_cor = 1.;
        }
        else if ((flag1 == _TRUE_) && (strstr(string1,"curvaton") != NULL)){      // inverted-correlation-curvaton case, only one iso parameter: piir1
          n_iso = ppm->n_s;
          n_cor = 0.;
          c_cor = -1.;
        }
        else{          // general case, but if pii2 or pri2=0, the code interprets it as a request for n_iso=n_ad or n_cor=0 respectively
          if (pii2 == 0.){
            n_iso = ppm->n_s;
          }
          else{
            class_test((pii1==0.) || (pii2 == 0.) || (pii1*pii2<0.),errmsg,"should NEVER happen");
            n_iso = log(pii2/pii1)/log(k2/k1)+1.;
          }
          class_test(pri1==0,errmsg,"the general isocurvature case requires a non-zero P_{RI}^1");
          if (pri2 == 0.){
            n_cor = 0.;
          }
          else{
            class_test((pri1==0.) || (pri2 <= 0.) || (pii1*pii2<0),errmsg,"should NEVER happen");
            n_cor = log(pri2/fabs(pri1))/log(k2/k1)-0.5*(ppm->n_s+n_iso-2.);
          }
          c_cor = -pri1/sqrt(pii1*prr1)*exp(n_cor*log(ppm->k_pivot/k1));
          /* Test */
          class_test((pii1*prr1<=0.),errmsg,"should NEVER happen");
          class_test(fabs(pri1)/sqrt(pii1*prr1)>1,errmsg,"too large ad-iso cross-correlation in k1");
          class_test(fabs(pri1)/sqrt(pii1*prr1)*exp(n_cor*log(k2/k1))>1,errmsg,"too large ad-iso cross-correlation in k2");
        }
        /* Complete set of parameters */
        class_test((pii1==0.) || (prr1 == 0.) || (pii1*prr1<0.),errmsg,"should NEVER happen");
        f_iso = sqrt(pii1/prr1)*exp(0.5*(n_iso-ppm->n_s)*log(ppm->k_pivot/k1));
      }
      if (ppt->has_bi == _TRUE_){
        ppm->f_bi = f_iso;
        ppm->n_bi = n_iso;
        ppm->c_ad_bi = c_cor;
        ppm->n_ad_bi = n_cor;
      }
      if (ppt->has_cdi == _TRUE_){
        ppm->f_cdi = f_iso;
        ppm->n_cdi = n_iso;
        ppm->c_ad_cdi = c_cor;
        ppm->n_ad_cdi = n_cor;
      }
      if (ppt->has_nid == _TRUE_){
        ppm->f_nid = f_iso;
        ppm->n_nid = n_iso;
        ppm->c_ad_nid = c_cor;
        ppm->n_ad_nid = n_cor;
      }
      if (ppt->has_niv == _TRUE_){
        ppm->f_niv = f_iso;
        ppm->n_niv = n_iso;
        ppm->c_ad_niv = c_cor;
        ppm->n_ad_niv = n_cor;
      }
    }
    ppm->primordial_spec_type = analytic_Pk;
  }

  /** 1.g) For type 'external_Pk' */
  else if (ppm->primordial_spec_type == external_Pk){

    /** 1.g.1) Command generating the table */
    /* Read */
    class_call(parser_read_string(pfc, "command", &string1, &flag1, errmsg),
               errmsg, errmsg);
    /* Test */
    class_test(strlen(string1) == 0,
               errmsg,
               "You omitted to write a command for the external Pk");
    /* Complete set of parameters */
    ppm->command = (char *) malloc (strlen(string1) + 1);
    strcpy(ppm->command, string1);

    /** 1.g.2) Command generating the table */
    /* Read */
    class_read_double("custom1",ppm->custom1);
    class_read_double("custom2",ppm->custom2);
    class_read_double("custom3",ppm->custom3);
    class_read_double("custom4",ppm->custom4);
    class_read_double("custom5",ppm->custom5);
    class_read_double("custom6",ppm->custom6);
    class_read_double("custom7",ppm->custom7);
    class_read_double("custom8",ppm->custom8);
    class_read_double("custom9",ppm->custom9);
    class_read_double("custom10",ppm->custom10);
  }

  /* Final tests */
  if ((ppm->primordial_spec_type == inflation_V) || (ppm->primordial_spec_type == inflation_H) || (ppm->primordial_spec_type == inflation_V_end)) {
    class_test(ppt->has_scalars == _FALSE_,
               errmsg,
               "inflationary module cannot work if you do not ask for scalar modes");
    class_test(ppt->has_vectors == _TRUE_,
               errmsg,
               "inflationary module cannot work if you ask for vector modes");
    class_test(ppt->has_tensors == _FALSE_,
               errmsg,
               "inflationary module cannot work if you do not ask for tensor modes");
    class_test(ppt->has_bi == _TRUE_ || ppt->has_cdi == _TRUE_ || ppt->has_nid == _TRUE_ || ppt->has_niv == _TRUE_,
               errmsg,
               "inflationary module cannot work if you ask for isocurvature modes");
  }

  return _SUCCESS_;

}


/**
 * Read the parameters of harmonic structure.
 *
 * @param pfc     Input: pointer to local structure
 * @param ppr     Input: pointer to precision structure
 * @param pba     Input: pointer to background structure
 * @param ppm     Input: pointer to primordial structure
 * @param ppt     Input: pointer to perturbations structure
 * @param ptr     Input: pointer to transfer structure
 * @param phr     Input: pointer to harmonic structure
 * @param pop     Input: pointer to output structure
 * @param errmsg  Input: Error message
 * @return the error status
 */

int input_read_parameters_spectra(struct file_content * pfc,
                                  struct precision * ppr,
                                  struct background * pba,
                                  struct primordial * ppm,
                                  struct perturbations * ppt,
                                  struct transfer * ptr,
                                  struct harmonic *phr,
                                  struct output * pop,
                                  ErrorMsg errmsg){

  /** Summary: */

  /** Define local variables */
  int flag1, flag2;
  double param1, param2;
  char string1[_ARGUMENT_LENGTH_MAX_];
  int int1;
  double * pointer1;
  int i;
  double z_max=0.;
  int bin;

  /** 1) Maximum l for CLs */
  /* Read */
  if (ppt->has_cls == _TRUE_) {
    if (ppt->has_scalars == _TRUE_) {
      if ((ppt->has_cl_cmb_temperature == _TRUE_) || (ppt->has_cl_cmb_polarization == _TRUE_) || (ppt->has_cl_cmb_lensing_potential == _TRUE_)){
        class_read_double("l_max_scalars",ppt->l_scalar_max);
      }
      if ((ppt->has_cl_lensing_potential == _TRUE_) || (ppt->has_cl_number_count == _TRUE_)){
        class_read_double("l_max_lss",ppt->l_lss_max);
      }
    }
    if (ppt->has_vectors == _TRUE_){
      class_read_double("l_max_vectors",ppt->l_vector_max);
    }
    if (ppt->has_tensors == _TRUE_) {
      class_read_double("l_max_tensors",ppt->l_tensor_max);
    }
  }


  /** 2) Parameters for the the matter density number count */
  if ((ppt->has_cl_number_count == _TRUE_) || (ppt->has_cl_lensing_potential == _TRUE_)) {

    /** 2.a) Selection functions W(z) of each redshift bin */
    /* Read */
    class_call(parser_read_string(pfc,"selection",&string1,&flag1,errmsg),
               errmsg,
               errmsg);
    /* Complete set of parameters */
    if (flag1 == _TRUE_) {
      if (strstr(string1,"gaussian") != NULL){
        ppt->selection=gaussian;
      }
      else if (strstr(string1,"tophat") != NULL){
        ppt->selection=tophat;
      }
      else if (strstr(string1,"dirac") != NULL){
        ppt->selection=dirac;
      }
      else {
        class_stop(errmsg,"You specified 'selection' as '%s'. It has to be one of {'gaussian','tophat','dirac'}.",string1);
      }
    }

    /* Read */
    class_call(parser_read_list_of_doubles(pfc,"selection_mean",&(int1),&(pointer1),&flag1,errmsg),
               errmsg,
               errmsg);
    if ((flag1 == _TRUE_) && (int1>0)) {
      /* Test */
      class_test(int1 > _SELECTION_NUM_MAX_,errmsg,
                 "you want to compute density Cl's for %d different bins, hence you should increase _SELECTION_NUM_MAX_ in include/perturbations.h to at least this number.",int1);
      /* Complete set of parameters */
      ppt->selection_num = int1;
      for (i=0; i<int1; i++) {
        /* Test */
        class_test((pointer1[i] < 0.) || (pointer1[i] > 1000.),errmsg,
                   "input of selection functions: you asked for a mean redshift equal to %e, is this a mistake?",pointer1[i]);
        /* Complete set of parameters */
        ppt->selection_mean[i] = pointer1[i];
      }
      free(pointer1);
      for (i=1; i<int1; i++) {       // first set all widths to default; correct eventually later
        /* Test */
        class_test(ppt->selection_mean[i]<=ppt->selection_mean[i-1],
                   errmsg,
                   "input of selection functions: the list of mean redshifts must be passed in growing order; you entered %e before %e.",ppt->selection_mean[i-1],ppt->selection_mean[i]);
        /* Complete set of parameters */
        ppt->selection_width[i] = ppt->selection_width[0];
        ptr->selection_bias[i] = ptr->selection_bias[0];
        ptr->selection_magnification_bias[i] = ptr->selection_magnification_bias[0];
      }

      /* Read */
      class_call(parser_read_list_of_doubles(pfc,"selection_width",&int1,&pointer1,&flag1,errmsg),
                 errmsg,
                 errmsg);
      /* Complete set of parameters */
      if ((flag1 == _TRUE_) && (int1>0)) {
        if (int1==1) {
          for (i=0; i<ppt->selection_num; i++) {
            ppt->selection_width[i] = pointer1[0];
          }
        }
        else if (int1==ppt->selection_num) {
          for (i=0; i<int1; i++) {
            ppt->selection_width[i] = pointer1[i];
          }
        }
        else {
          class_stop(errmsg,
                     "In input for selection function, you asked for %d bin centers and %d bin widths; number of bins unclear; you should pass either one bin width (common to all bins) or %d bin widths.",ppt->selection_num,int1,ppt->selection_num);
        }
        free(pointer1);
      }

      /* Read */
      class_call(parser_read_list_of_doubles(pfc,"selection_bias",&int1,&pointer1,&flag1,errmsg),
                 errmsg,
                 errmsg);
      /* Complete set of parameters */
      if ((flag1 == _TRUE_) && (int1>0)) {
        if (int1==1) {
          for (i=0; i<ppt->selection_num; i++) {
            ptr->selection_bias[i] = pointer1[0];
          }
        }
        else if (int1==ppt->selection_num) {
          for (i=0; i<int1; i++) {
            ptr->selection_bias[i] = pointer1[i];
          }
        }
        else {
          class_stop(errmsg,
                     "In input for selection function, you asked for %d bin centers and %d bin biases; number of bins unclear; you should pass either one bin bias (common to all bins) or %d bin biases.",
                     ppt->selection_num,int1,ppt->selection_num);
        }
        free(pointer1);
      }

      /* Read */
      class_call(parser_read_list_of_doubles(pfc,"selection_magnification_bias",&int1,&pointer1,&flag1,errmsg),
                 errmsg,
                 errmsg);
      /* Complete set of parameters */
      if ((flag1 == _TRUE_) && (int1>0)) {
        if (int1==1) {
          for (i=0; i<ppt->selection_num; i++) {
            ptr->selection_magnification_bias[i] = pointer1[0];
          }
        }
        else if (int1==ppt->selection_num) {
          for (i=0; i<int1; i++) {
            ptr->selection_magnification_bias[i] = pointer1[i];
          }
        }
        else {
          class_stop(errmsg,
                     "In input for selection function, you asked for %d bin centers and %d bin biases; number of bins unclear; you should pass either one bin bias (common to all bins) or %d bin biases.",
                     ppt->selection_num,int1,ppt->selection_num);
        }
        free(pointer1);
      }
    }

    /* Read */
    if (ppt->selection_num > 1) {
      class_read_int("non_diagonal",phr->non_diag);
      if ((phr->non_diag<0) || (phr->non_diag>=ppt->selection_num))
        class_stop(errmsg,"Input for non_diagonal is %d, while it is expected to be between 0 and %d\n",
                   phr->non_diag,ppt->selection_num-1);
    }

    /** 2.b) Selection function */
    /* Read */
    class_call(parser_read_string(pfc,"dNdz_selection",&string1,&flag1,errmsg),
               errmsg,
               errmsg);
    /* Complete set of parameters */
    if ((flag1 == _TRUE_)) {
      if ((strstr(string1,"analytic") != NULL)){
        ptr->has_nz_analytic = _TRUE_;
      }
      else{
        ptr->has_nz_file = _TRUE_;
        class_read_string("dNdz_selection",ptr->nz_file_name);
      }
    }

    /** 2.c) Source number counts evolution */
    class_call(parser_read_string(pfc,"dNdz_evolution",&string1,&flag1,errmsg),
               errmsg,
               errmsg);
    /* Complete set of parameters */
    if ((flag1 == _TRUE_)) {
      if ((strstr(string1,"analytic") != NULL)){
        ptr->has_nz_evo_analytic = _TRUE_;
      }
      else{
        ptr->has_nz_evo_file = _TRUE_;
        class_read_string("dNdz_evolution",ptr->nz_evo_file_name);
      }
    }

  }


  /** 3) Power spectrum P(k) */
  if ((ppt->has_pk_matter == _TRUE_) || (ppt->has_density_transfers == _TRUE_) || (ppt->has_velocity_transfers == _TRUE_)){

    /** 3.a) Maximum k in P(k) */
    /* Read */
    class_call(parser_read_double(pfc,"P_k_max_h/Mpc",&param1,&flag1,errmsg),
               errmsg,
               errmsg);
    class_call(parser_read_double(pfc,"P_k_max_1/Mpc",&param2,&flag2,errmsg),
               errmsg,
               errmsg);
    /* Test */
    class_test((flag1 == _TRUE_) && (flag2 == _TRUE_),
               errmsg,
               "You can only enter one of 'P_k_max_h/Mpc' or 'P_k_max_1/Mpc'.");
    /* Complete set of parameters */
    if (flag1 == _TRUE_){
      ppt->k_max_for_pk=param1*pba->h;
    }
    if (flag2 == _TRUE_){
      ppt->k_max_for_pk=param2;
    }

    /** 3.a.1) Maximum k in primordial P(k) */
    /* Read */
    class_call(parser_read_double(pfc,"primordial_P_k_max_h/Mpc",&param1,&flag1,errmsg),
               errmsg,
               errmsg);
    class_call(parser_read_double(pfc,"primordial_P_k_max_1/Mpc",&param2,&flag2,errmsg),
               errmsg,
               errmsg);
    /* Test */
    class_test((flag1 == _TRUE_) && (flag2 == _TRUE_),
               errmsg,
               "You can only enter one of 'primordial_P_k_max_h/Mpc' or 'primordial_P_k_max_1/Mpc'.");
    /* Complete set of parameters */
    if (flag1 == _TRUE_){
      ppm->k_max_for_primordial_pk=param1*pba->h;
      ppm->has_k_max_for_primordial_pk = _TRUE_;
    }
    if (flag2 == _TRUE_){
      ppm->k_max_for_primordial_pk=param2;
      ppm->has_k_max_for_primordial_pk = _TRUE_;
    }

    /** 3.b) Redshift values */
    /* Read */
    class_call(parser_read_list_of_doubles(pfc,"z_pk",&int1,&pointer1,&flag1,errmsg),
               errmsg,
               errmsg);
    /* Test */
    if (flag1 == _TRUE_) {
      class_test(int1 > _Z_PK_NUM_MAX_,
                 errmsg,
                 "you want to write some output for %d different values of z, hence you should increase _Z_PK_NUM_MAX_ in include/output.h to at least this number",
                 int1);
      /* Complete set of parameters */
      pop->z_pk_num = int1;
      for (i=0; i<int1; i++) {
        pop->z_pk[i] = pointer1[i];
      }
      free(pointer1);
    }

  }

  /** 3.c) Maximum redshift */
  if ((ppt->has_pk_matter == _TRUE_) || (ppt->has_density_transfers == _TRUE_) || (ppt->has_velocity_transfers == _TRUE_) || (ppt->has_cl_number_count == _TRUE_) || (ppt->has_cl_lensing_potential == _TRUE_)) {
    /* Read */
    class_call(parser_read_double(pfc,"z_max_pk",&param1,&flag1,errmsg),
               errmsg,
               errmsg);

    /* Complete set of parameters */
    if (flag1==_TRUE_) {
      ppt->z_max_pk = param1;
    }
    /* *
     * If we could not read a z_max, we need to define one.
     * The limit could come from any of the contributions.
     * We test here, which contribution requires the largest z_max
     * */
    else {
      ppt->z_max_pk = 0.;
      /* For the z_pk related quantities, test here the z_pk requirements */
      if ((ppt->has_pk_matter == _TRUE_) || (ppt->has_density_transfers == _TRUE_) || (ppt->has_velocity_transfers == _TRUE_)) {
        for (i=0; i<pop->z_pk_num; i++) {
          ppt->z_max_pk = MAX(ppt->z_max_pk,pop->z_pk[i]);
        }
      }
      /* For the number count / shear related quantities, test the selection function z_max */
      if ((ppt->has_cl_number_count == _TRUE_) || (ppt->has_cl_lensing_potential == _TRUE_)){
        for (bin=0; bin<ppt->selection_num; bin++) {
          /* the few lines below should be consistent with their counterpart in transfer.c, in transfer_selection_times */
          if (ppt->selection==gaussian) {
            z_max = ppt->selection_mean[bin]+ppt->selection_width[bin]*ppr->selection_cut_at_sigma;
          }
          if (ppt->selection==tophat) {
            z_max = ppt->selection_mean[bin]+(1.+ppr->selection_cut_at_sigma*ppr->selection_tophat_edge)*ppt->selection_width[bin];
          }
          if (ppt->selection==dirac) {
            z_max = ppt->selection_mean[bin];
          }
          ppt->z_max_pk = MAX(ppt->z_max_pk,z_max);
        }
      }
      /* Now we have checked all contributions that could change z_max_pk */
    }
  }

  return _SUCCESS_;

}


/**
 * Read the parameters of lensing structure.
 *
 * @param pfc     Input: pointer to local structure
 * @param ppr     Input: pointer to precision structure
 * @param ppt     Input: pointer to perturbations structure
 * @param ptr     Input: pointer to transfer structure
 * @param ple     Input: pointer to lensing structure
 * @param errmsg  Input: Error message
 * @return the error status
 */

int input_read_parameters_lensing(struct file_content * pfc,
                                  struct precision * ppr,
                                  struct perturbations * ppt,
                                  struct transfer * ptr,
                                  struct lensing *ple,
                                  ErrorMsg errmsg){

  /** Summary: */

  /** Define local variables */
  int flag1,flag2;
  double param1,param2;
  char string1[_ARGUMENT_LENGTH_MAX_];

  /** 1) Lensed spectra? */
  /* Read */
  class_call(parser_read_string(pfc,"lensing",&string1,&flag1,errmsg),
             errmsg,
             errmsg);
  /* Complete set of parameters */
  if ((flag1 == _TRUE_) && (string_begins_with(string1,'y') || string_begins_with(string1,'Y'))){
    if ((ppt->has_scalars == _TRUE_) && ((ppt->has_cl_cmb_temperature == _TRUE_) || (ppt->has_cl_cmb_polarization == _TRUE_)) && (ppt->has_cl_cmb_lensing_potential == _TRUE_)){
      ple->has_lensed_cls = _TRUE_;
      /* Slightly increase precision by delta_l_max for more precise lensed Cl's*/
      ppt->l_scalar_max += ppr->delta_l_max;
    }
    else {
      class_stop(errmsg,"you asked for lensed CMB Cls, but this requires a minimal number of options: 'modes' should include 's', 'output' should include 'tCl' and/or 'pCl', and also, importantly, 'lCl', the CMB lensing potential spectrum.");
    }
  }


  /** 2) Should the lensed spectra be rescaled (either with just A_L, or otherwise with amplitude, and tilt and pivot scale in k space) */
  /* Read */
  if ((ppt->has_scalars == _TRUE_) && (ppt->has_cl_cmb_lensing_potential == _TRUE_)) {

    class_call(parser_read_double(pfc,"A_L",&param1,&flag1,errmsg),
               errmsg,
               errmsg);
    class_call(parser_read_double(pfc,"lcmb_rescale",&param2,&flag2,errmsg),
               errmsg,
               errmsg);

    if ((flag1 == _TRUE_) && (flag2 == _TRUE_)) {
      class_stop(errmsg,"You cannot pass both A_l and lcdmb_rescale, choose one");
    }
    else {
      if (flag1 == _TRUE_) {
        ptr->lcmb_rescale = sqrt(param1);
      }
      if (flag2 == _TRUE_) {
        ptr->lcmb_rescale = param2;
      }
    }
    class_read_double("lcmb_tilt",ptr->lcmb_tilt);
    class_read_double("lcmb_pivot",ptr->lcmb_pivot);
  }

  return _SUCCESS_;

}


/**
 * Read free parameters of distortions structure.
 *
 * @param pfc     Input: pointer to local structure
 * @param ppr     Input: pointer to precision structure
 * @param psd     Input: pointer to distortions structure
 * @param errmsg  Input: Error message
 * @return the error status
 */

int input_read_parameters_distortions(struct file_content * pfc,
                                      struct precision * ppr,
                                      struct distortions * psd,
                                      ErrorMsg errmsg){

  /** Summary: */

  /** Define local variables */
  int flag1, flag2;
  char string1[_ARGUMENT_LENGTH_MAX_];
  double param1, param2;
  double updated_nu_max;

  /** 1) Branching ratio approximation */
  /* Read */
  class_call(parser_read_string(pfc,"sd_branching_approx",&string1,&flag1,errmsg),
             errmsg,
             errmsg);
  /* Complete set of parameters */

  if(flag1 == _TRUE_){
    if ( (strstr(string1,"sharp_sharp") != NULL) || (strstr(string1,"sharp sharp") != NULL) ) {
      psd->sd_branching_approx = bra_sharp_sharp;
      psd->sd_PCA_size = 0;
    }
    else if ( (strstr(string1,"sharp_soft") != NULL) || (strstr(string1,"sharp soft") != NULL) ) {
      psd->sd_branching_approx = bra_sharp_soft;
      psd->sd_PCA_size = 0;
    }
    else if ( (strstr(string1,"soft_soft") != NULL) || (strstr(string1,"soft soft") != NULL) ) {
      psd->sd_branching_approx = bra_soft_soft;
      psd->sd_PCA_size = 0;
    }
    else if ( (strstr(string1,"soft_soft_cons") != NULL) || (strstr(string1,"soft soft cons") != NULL) ) {
      psd->sd_branching_approx = bra_soft_soft_cons;
      psd->sd_PCA_size = 0;
    }
    else if ( (strstr(string1,"exact") != NULL) ) {
      psd->sd_branching_approx = bra_exact;
    }
    else{
      class_stop(errmsg,"You specified 'branching_approx' as '%s'. It has to be one of {'sharp_sharp','sharp_soft','soft_soft','soft_soft_cons','exact'}.",string1);
    }
  }

  /* Only read these if 'bra_exact' has been set (could also be set from default) */
  if(psd->sd_branching_approx == bra_exact){

    /** 1.a.1) Number of multipoles in PCA expansion */
    /* Read */
    class_read_int("sd_PCA_size",psd->sd_PCA_size);
    /* Test */
    if(psd->sd_PCA_size < 0 || psd->sd_PCA_size > 6){
      psd->sd_PCA_size = 6;
    }

    /** 1.a.2) Detector name */
    /* Read */
    class_call(parser_read_string(pfc,"sd_detector_name",&string1,&flag1,errmsg),
               errmsg,
               errmsg);
    /* Complete set of parameters */
    if(flag1 == _TRUE_){
      strcpy(psd->sd_detector_name,string1);
      psd->has_user_defined_name = _TRUE_;
    }

    /** 1.a.3) Detector specifics */
    /** 1.a.3.1) From file */
    /* Read */
    class_call(parser_read_string(pfc,"sd_detector_file",&string1,&flag1,errmsg),
               errmsg,
               errmsg);
    /* Complete set of parameters */
    if(flag1 == _TRUE_){
      strcpy(psd->sd_detector_file_name,string1);
      psd->has_detector_file = _TRUE_;
    }

    /** 1.a.3.2) User defined */
    /* Read */
    class_call(parser_read_double(pfc,"sd_detector_nu_min",&param1,&flag1,errmsg),
               errmsg,
               errmsg);
    /* Complete set of parameters */
    if(flag1 == _TRUE_){
      psd->sd_detector_nu_min = param1;
      psd->has_user_defined_detector = _TRUE_;
    }
    /* Read */
    class_call(parser_read_double(pfc,"sd_detector_nu_max",&param1,&flag1,errmsg),
               errmsg,
               errmsg);
    /* Complete set of parameters */
    if(flag1 == _TRUE_){
      psd->sd_detector_nu_max = param1;
      psd->has_user_defined_detector = _TRUE_;
    }
    /* Read */
    class_call(parser_read_double(pfc,"sd_detector_nu_delta",&param1,&flag1,errmsg),
               errmsg,
               errmsg);
     class_call(parser_read_double(pfc,"sd_detector_bin_number",&param2,&flag2,errmsg),
               errmsg,
               errmsg);
    /* Test */
    class_test((flag1 == _TRUE_) && (flag2 == _TRUE_),
               errmsg,
               "You can only enter one of 'sd_detector_nu_delta' or 'sd_detector_bin_number'.",
               psd->sd_detector_nu_delta,psd->sd_detector_bin_number);
    /* Complete set of parameters */
    if(flag1 == _TRUE_){
      psd->sd_detector_nu_delta = param1;
      psd->sd_detector_bin_number = ((int)ceil((psd->sd_detector_nu_max-psd->sd_detector_nu_min)/param1));
      psd->has_user_defined_detector = _TRUE_;
    }
    if(flag2 == _TRUE_){
      psd->sd_detector_nu_delta = (psd->sd_detector_nu_max-psd->sd_detector_nu_min)/param2;
      psd->sd_detector_bin_number = param2;
      psd->has_user_defined_detector = _TRUE_;
   }
    /* Update value of nu_max, given the number of bins */
    updated_nu_max = psd->sd_detector_nu_min+psd->sd_detector_nu_delta*psd->sd_detector_bin_number;
    if(fabs(updated_nu_max-psd->sd_detector_nu_max) > ppr->tol_sd_detector){
      printf(" -> WARNING: The value of 'sd_detector_nu_max' has been updated to %7.3e to accommodate the binning of your detector.\n",updated_nu_max);
      psd->sd_detector_nu_max = updated_nu_max;
    }
    /* Read */
    class_call(parser_read_double(pfc,"sd_detector_delta_Ic",&param1,&flag1,errmsg),
               errmsg,
               errmsg);
    /* Complete set of parameters */
    if(flag1 == _TRUE_){
      psd->sd_detector_delta_Ic = 1.0e-26*param1;
      psd->has_user_defined_detector = _TRUE_;
    }
  }

  /* Final tests */
  class_test(psd->sd_branching_approx != bra_exact && psd->sd_PCA_size > 0,
             errmsg,
             "The PCA expansion is possible only for 'branching_approx = exact'");
  class_test(psd->has_detector_file && psd->has_user_defined_detector,
             errmsg,
             "You can only enter the noise file {'%s'} or the specifications {'%s','%s'/'%s','%s'}.",
             "sd_detector_file","sd_detector_nu_min","sd_detector_nu_max",
             "sd_detector_nu_delta","sd_detector_bin_number","sd_detector_delta_Ic");


  /** 2) Only calculate exotic energy injections and no LCDM processes for spectral distortions ? */
  class_read_flag("sd_only_exotic",psd->include_only_exotic);

  /** 3) Include g distortions? */
  class_read_flag("sd_include_g_distortion",psd->include_g_distortion);


  /** 4) Set g-distortions to zero? */
  class_call(parser_read_double(pfc,"sd_add_y",&psd->sd_add_y,&flag1,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"sd_add_mu",&psd->sd_add_mu,&flag1,errmsg),
             errmsg,
             errmsg);


  /** 5) Include SZ effect from reionization? */
  class_read_flag("include_SZ_effect",psd->has_SZ_effect);

  if(psd->has_SZ_effect == _TRUE_){
    /** 5.a) Type of calculation */
    /* Read */
    class_call(parser_read_string(pfc,"sd_reio_type",&string1,&flag1,errmsg),
               errmsg,
               errmsg);
    /* Complete set of parameters */
    if (flag1 == _TRUE_){
      if (strcmp(string1,"Nozawa_2005") == 0){
        psd->sd_reio_type = sd_reio_Nozawa;
      }
      else if (strcmp(string1,"Chluba_2012") == 0){
        psd->sd_reio_type = sd_reio_Chluba;
      }
      else{
        class_stop(errmsg,
                   "You specified 'sd_reio_type' as '%s'. It has to be one of {'Nozawa_2005','Chluba_2012'}.",string1);
      }
    }
  }

  return _SUCCESS_;

}


/**
 * Read obsolete/additional parameters that are not assigned to a specific structure
 *
 * @param pfc     Input: pointer to local structure
 * @param ppr     Input: pointer to precision structure
 * @param pba     Input: pointer to background structure
 * @param pth     Input: pointer to thermodynamics structure
 * @param errmsg  Input: Error message
 * @return the error status
 */

int input_read_parameters_additional(struct file_content* pfc,
                                     struct precision* ppr,
                                     struct background* pba,
                                     struct thermodynamics* pth,
                                     ErrorMsg errmsg){

  /** Summary: */

  /** Define local variables */
  int flag1;
  double param1;
  char string1[_ARGUMENT_LENGTH_MAX_];

  /**
   * Here we can place all obsolete (deprecated) names for the precision parameters
   * that will still be read as of the current version.
   * There is however, no guarantee that this will be true for future versions as well.
   * The new parameter names should be used preferrably.
   * */
  class_read_double("k_scalar_min_tau0",ppr->k_min_tau0);
  class_read_double("k_scalar_max_tau0_over_l_max",ppr->k_max_tau0_over_l_max);
  class_read_double("k_scalar_step_sub",ppr->k_step_sub);
  class_read_double("k_scalar_step_super",ppr->k_step_super);
  class_read_double("k_scalar_step_transition",ppr->k_step_transition);
  class_read_double("k_scalar_k_per_decade_for_pk",ppr->k_per_decade_for_pk);
  class_read_double("k_scalar_k_per_decade_for_bao",ppr->k_per_decade_for_bao);
  class_read_double("k_scalar_bao_center",ppr->k_bao_center);
  class_read_double("k_scalar_bao_width",ppr->k_bao_width);

  class_read_double("k_step_trans_scalars",ppr->q_linstep);
  class_read_double("k_step_trans_tensors",ppr->q_linstep);
  class_read_double("k_step_trans",ppr->q_linstep);
  class_read_double("q_linstep_trans",ppr->q_linstep);
  class_read_double("q_logstep_trans",ppr->q_logstep_spline);

  /** Here are slgihtly more obsolete parameters, these will not even be read, only give an error message */
  class_call(parser_read_double(pfc,"bias",&param1,&flag1,errmsg),
             errmsg,
             errmsg);
  class_test(flag1 == _TRUE_,
             errmsg,
             "the input parameter 'bias' is obsolete, because you can now pass an independent light-to-mass bias for each bin/selection function. The new input name is 'selection_bias'. It can be set to a single number (common bias for all bins) or as many numbers as the number of bins");

  class_call(parser_read_double(pfc,"s_bias",&param1,&flag1,errmsg),
             errmsg,
             errmsg);
  class_test(flag1 == _TRUE_,
             errmsg,
             "the input parameter 's_bias' is obsolete, because you can now pass an independent magnitude bias for each bin/selection function. The new input name is 'selection_magnitude_bias'. It can be set to a single number (common magnitude bias for all bins) or as many numbers as the number of bins");

  class_call(parser_read_string(pfc,"l_switch_limber_for_cl_density_over_z",&string1,&flag1,errmsg),
             errmsg,
             errmsg);
  class_test(flag1 == _TRUE_,
             errmsg,
             "You passed in input a precision parameter called l_switch_limber_for_cl_density_over_z. This syntax is deprecated since v2.5.0. Please use instead the two precision parameters l_switch_limber_for_nc_local_over_z, l_switch_limber_for_nc_los_over_z, defined in include/common.h, and allowing for better performance.");

  /*
  class_call(parser_read_double(pfc,"annihilation",&param1,&flag1,errmsg),
             errmsg,
             errmsg);
  class_test(flag1 == _TRUE_,
             errmsg,
             "the input parameter 'annihilation' (that is, f_eff <sigma_DM v> / m_DM in units m^3/s/Kg) is obsolete and replaced by 'DM_annihilation_efficiency' (that is, <sigma_DM v> / m_DM in units m^3/s/J). If 'f_eff_type' is set to 'on_the_spot', f_eff is assumed to be one and the old/new parameters only differ by a factor c^2 in m2/s2 = 9.e16. For instance, annihilation=1.e-6 is the same as DM_annihilation_efficiency=9.e-22".
  */


  /** Test additional input parameters related to precision parameters */
  if ((ppr->tight_coupling_approximation == (int)first_order_CLASS) || (ppr->tight_coupling_approximation == (int)second_order_CLASS)) {
    pth->compute_cb2_derivatives = _TRUE_;
  }

  class_test(ppr->ur_fluid_trigger_tau_over_tau_k==ppr->radiation_streaming_trigger_tau_over_tau_k,
             errmsg,
             "please choose different values for precision parameters ur_fluid_trigger_tau_over_tau_k and radiation_streaming_trigger_tau_over_tau_k, in order to avoid switching two approximation schemes at the same time");

  if (pba->N_ncdm>0) {
    class_test(ppr->ncdm_fluid_trigger_tau_over_tau_k==ppr->radiation_streaming_trigger_tau_over_tau_k,
               errmsg,
               "please choose different values for precision parameters ncdm_fluid_trigger_tau_over_tau_k and radiation_streaming_trigger_tau_over_tau_k, in order to avoid switching two approximation schemes at the same time");

    class_test(ppr->ncdm_fluid_trigger_tau_over_tau_k==ppr->ur_fluid_trigger_tau_over_tau_k,
               errmsg,
               "please choose different values for precision parameters ncdm_fluid_trigger_tau_over_tau_k and ur_fluid_trigger_tau_over_tau_k, in order to avoid switching two approximation schemes at the same time");
  }

  if (pba->Omega0_idr != 0.){
    class_test(ppr->idr_streaming_trigger_tau_over_tau_k==ppr->radiation_streaming_trigger_tau_over_tau_k,
               errmsg,
               "please choose different values for precision parameters dark_radiation_trigger_tau_over_tau_k and radiation_streaming_trigger_tau_over_tau_k, in order to avoid switching two approximation schemes at the same time");

    class_test(ppr->idr_streaming_trigger_tau_over_tau_k==ppr->ur_fluid_trigger_tau_over_tau_k,
               errmsg,
               "please choose different values for precision parameters dark_radiation_trigger_tau_over_tau_k and ur_fluid_trigger_tau_over_tau_k, in order to avoid switching two approximation schemes at the same time");

    class_test(ppr->idr_streaming_trigger_tau_over_tau_k==ppr->ncdm_fluid_trigger_tau_over_tau_k,
               errmsg,
               "please choose different values for precision parameters dark_radiation_trigger_tau_over_tau_k and ncdm_fluid_trigger_tau_over_tau_k, in order to avoid switching two approximation schemes at the same time");
  }

  return _SUCCESS_;

}


/**
 * Read the parameters of output structure.
 *
 * @param pfc     Input: pointer to local structure
 * @param pba     Input: pointer to background structure
 * @param pth     Input: pointer to thermodynamics structure
 * @param ppt     Input: pointer to perturbations structure
 * @param ptr     Input: pointer to transfer structure
 * @param ppm     Input: pointer to primordial structure
 * @param phr     Input: pointer to harmonic structure
 * @param pfo     Input: pointer to non-linear structure
 * @param ple     Input: pointer to lensing structure
 * @param psd     Input: pointer to distorsion structure
 * @param pop     Input: pointer to output structure
 * @param errmsg  Input: Error message
 * @return the error status
 */

int input_read_parameters_output(struct file_content * pfc,
                                 struct background *pba,
                                 struct thermodynamics *pth,
                                 struct perturbations *ppt,
                                 struct transfer *ptr,
                                 struct primordial *ppm,
                                 struct harmonic *phr,
                                 struct fourier * pfo,
                                 struct lensing *ple,
                                 struct distortions *psd,
                                 struct output *pop,
                                 ErrorMsg errmsg){

  /** Summary: */

  /** Define local variables */
  int flag1;
  char string1[_ARGUMENT_LENGTH_MAX_];
  int int1;
  double * pointer1;
  int i;
  FILE * param_output;
  FILE * param_unused;
  char param_output_name[_LINE_LENGTH_MAX_];
  char param_unused_name[_LINE_LENGTH_MAX_];

  /** 1) Output for external files */
  /** 1.a) File name */
  /* Read */
  class_call(parser_read_string(pfc,"root",&string1,&flag1,errmsg),
             errmsg,
             errmsg);
  /* Complete set of parameters */
  if (flag1 == _TRUE_){
    class_test(strlen(string1)>_FILENAMESIZE_-32,errmsg,"Root directory name is too long. Please install in other directory, or increase _FILENAMESIZE_ in common.h");
    strcpy(pop->root,string1);
  }

  /** 1.b) Headers */
  /* Read */
  class_read_flag("headers",pop->write_header);

  /** 1.c) Format */
  /* Read */
  class_call(parser_read_string(pfc,"format",&string1,&flag1,errmsg),
             errmsg,
             errmsg);
  /* Complete set of parameters */
  if (flag1 == _TRUE_){
    if ((strstr(string1,"class") != NULL) || (strstr(string1,"CLASS") != NULL)){
      pop->output_format = class_format;
    }
    else if ((strstr(string1,"camb") != NULL) || (strstr(string1,"CAMB") != NULL)){
      pop->output_format = camb_format;
    }
    else{
      class_stop(errmsg,"You specified 'format' as '%s'. It has to be one of {'class','camb'}.",string1);
    }
  }

  /** 1.d) Background quantities */
  /* Read */
  class_read_flag_or_deprecated("write_background","write background",pop->write_background);

  /** 1.e) Thermodynamics quantities */
  /* Read */
  class_read_flag_or_deprecated("write_thermodynamics","write thermodynamics",pop->write_thermodynamics);

  /** 1.f) Table of perturbations for certain wavenumbers k */
  /* Read */
  class_call(parser_read_list_of_doubles(pfc,"k_output_values",&int1,&pointer1,&flag1,errmsg),
             errmsg,
             errmsg);
  if (flag1 == _TRUE_) {
    /* Test */
    class_test(int1 > _MAX_NUMBER_OF_K_FILES_,
               errmsg,
               "you want to write some output for %d different values of k, hence you should increase _MAX_NUMBER_OF_K_FILES_ in include/perturbations.h to at least this number",
               int1);
    /* Complete set of parameters */
    ppt->k_output_values_num = int1;
    for (i=0; i<int1; i++) {
      ppt->k_output_values[i] = pointer1[i];
    }
    free(pointer1);
    qsort (ppt->k_output_values, ppt->k_output_values_num, sizeof(double), compare_doubles);     // Sort the k_array using qsort
    ppt->store_perturbations = _TRUE_;
    pop->write_perturbations = _TRUE_;
  }

  /** 1.g) Primordial spectra */
  /* Read */
  class_read_flag_or_deprecated("write_primordial","write primordial",pop->write_primordial);

  /** 1.h) Exotic energy injection output */
  /* Read */
  class_read_flag_or_deprecated("write_exotic_injection","write exotic injection",pop->write_exotic_injection);

  /** 1.i) Non-injected photon injection */
  /* Read */
  class_read_flag_or_deprecated("write_noninjection","write noninjection",pop->write_noninjection);

  /** 1.k) Spectral Distortions */
  /* Read */
  class_read_flag_or_deprecated("write_distortions","write distortions",pop->write_distortions);

  /** 1.l) Input/precision parameters */
  /* Read */
  flag1 = _FALSE_;
  class_read_flag_or_deprecated("write_parameters","write parameters",flag1);


  /** 2) Verbosity */
  /* Read */
  class_read_int("background_verbose",pba->background_verbose);
  class_read_int("thermodynamics_verbose",pth->thermodynamics_verbose);
  class_read_int("hyrec_verbose",pth->hyrec_verbose);
  class_read_int("perturbations_verbose",ppt->perturbations_verbose);
  class_read_int("transfer_verbose",ptr->transfer_verbose);
  class_read_int("primordial_verbose",ppm->primordial_verbose);
  class_read_int("harmonic_verbose",phr->harmonic_verbose);
  class_read_int("fourier_verbose",pfo->fourier_verbose);
  class_read_int("lensing_verbose",ple->lensing_verbose);
  class_read_int("distortions_verbose",psd->distortions_verbose);
  class_read_int("output_verbose",pop->output_verbose);


  /**
   * This must be the very LAST entry of read_parameters,
   * since it relies on the pfc->read flags being set to _TRUE_ or _FALSE_
   * */
  /* Set the parameters.ini and unused_parameters files */
  if (flag1 == _TRUE_) {
    sprintf(param_output_name,"%s%s",pop->root,"parameters.ini");
    class_open(param_output,param_output_name,"w",errmsg);
    fprintf(param_output,"# List of input/precision parameters actually read\n");
    fprintf(param_output,"# (all other parameters set to default values)\n");
    fprintf(param_output,"# Obtained with CLASS %s (for developers: svn version %s)\n",_VERSION_,_SVN_VERSION_);
    fprintf(param_output,"#\n");
    fprintf(param_output,"# This file can be used as the input file of another run\n");
    fprintf(param_output,"#\n");

    sprintf(param_unused_name,"%s%s",pop->root,"unused_parameters");
    class_open(param_unused,param_unused_name,"w",errmsg);
    fprintf(param_unused,"# List of input/precision parameters passed\n");
    fprintf(param_unused,"# but not used (just for info)\n");
    fprintf(param_unused,"#\n");

    for (i=0; i<pfc->size; i++) {
      if (pfc->read[i] == _TRUE_){
        fprintf(param_output,"%s = %s\n",pfc->name[i],pfc->value[i]);
      }
      else{
        fprintf(param_unused,"%s = %s\n",pfc->name[i],pfc->value[i]);
      }
    }

    fprintf(param_output,"#\n");

    fclose(param_output);
    fclose(param_unused);
  }

  return _SUCCESS_;

}

/**
 * All default parameter values (for input parameters)
 *
 * @param pba Input: pointer to background structure
 * @param pth Input: pointer to thermodynamics structure
 * @param ppt Input: pointer to perturbation structure
 * @param ptr Input: pointer to transfer structure
 * @param ppm Input: pointer to primordial structure
 * @param phr Input: pointer to harmonic structure
 * @param pfo Input: pointer to fourier structure
 * @param ple Input: pointer to lensing structure
 * @param psd     Input: pointer to distorsion structure
 * @param pop Input: pointer to output structure
 * @return the error status
 * @return the error status
 */

int input_default_params(struct background *pba,
                         struct thermodynamics *pth,
                         struct perturbations *ppt,
                         struct transfer *ptr,
                         struct primordial *ppm,
                         struct harmonic *phr,
                         struct fourier * pfo,
                         struct lensing *ple,
                         struct distortions *psd,
                         struct output *pop) {

  /** Summary: */

  /** - Define local variables */
  struct injection* pin = &(pth->in);
  double sigma_B; /* Stefan-Boltzmann constant in \f$ W/m^2/K^4 = Kg/K^4/s^3 \f$*/

  sigma_B = 2. * pow(_PI_,5) * pow(_k_B_,4) / 15. / pow(_h_P_,3) / pow(_c_,2);

  /* 24.02.2021: default parameters.
     If you want to use exactly the Planck 2018 bestfit model, you can always use one of our
     suggested input files, e.g., base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini.
     The default values implemented here are an approximation to this LambdaCDM model
     with the same values for
     omega_b=0.02238280,
     omega_c=0.1201075,
     theta_s=1.041783 (this is not theta_MC nor theta*),
     tau_reio=0.05430842,
     As=2.100549e-09,
     ns=0.9660499,
     but a few differences:
     only massless neutrinos,
     N_eff=3.044 instead of 3.046 (3.044 is a more accurate and updated estimate of the
     effective neutrino number in LCDM),
     non-linear corrections neglected */

  /**
   * Default to input_read_parameters_general
   */

  /** 1) Output spectra */
  ppt->has_cl_cmb_temperature = _FALSE_;
  ppt->has_cl_cmb_polarization = _FALSE_;
  ppt->has_cl_cmb_lensing_potential = _FALSE_;
  ppt->has_cl_number_count = _FALSE_;
  ppt->has_cl_lensing_potential = _FALSE_;
  ppt->has_pk_matter = _FALSE_;
  ppt->has_density_transfers = _FALSE_;
  ppt->has_velocity_transfers = _FALSE_;
  /** 1.a) 'tCl' case */
  ppt->switch_sw = 1;
  ppt->switch_eisw = 1;
  ppt->switch_lisw = 1;
  ppt->switch_dop = 1;
  ppt->switch_pol = 1;
  /** 1.a.1) Split value of redshift z at which the isw is considered as late or early */
  ppt->eisw_lisw_split_z = 120;
  /** 1.b) 'nCl' (or 'dCl') case */
  ppt->has_nc_density = _FALSE_;
  ppt->has_nc_rsd = _FALSE_;
  ppt->has_nc_lens = _FALSE_;
  ppt->has_nc_gr = _FALSE_;
  /** 1.c) 'dTk' (or 'mTk') case */
  ppt->has_metricpotential_transfers = _FALSE_;

  /** 2) Perturbed recombination */
  ppt->has_perturbed_recombination=_FALSE_;
  /** 3) Modes */
  ppt->has_scalars=_TRUE_;
  ppt->has_vectors=_FALSE_;
  ppt->has_tensors=_FALSE_;
  /** 3.a) Initial conditions for scalars */
  ppt->has_ad=_TRUE_;
  ppt->has_bi=_FALSE_;
  ppt->has_cdi=_FALSE_;
  ppt->has_nid=_FALSE_;
  ppt->has_niv=_FALSE_;
  /** 3.b) Initial conditions for tensors */
  ppt->tensor_method = tm_massless_approximation;
  ppt->evolve_tensor_ur = _FALSE_;
  ppt->evolve_tensor_ncdm = _FALSE_;

  /** 4.a) Gauge */
  ppt->gauge=synchronous;
  /** 4.b) N-body gauge */
  ppt->has_Nbody_gauge_transfers = _FALSE_;

  /** 5) Hubble parameter */
  pba->h = 0.67810;
  pba->H0 = pba->h*1.e5/_c_;

  /** 6) Primordial Helium fraction */
  pth->YHe = _YHE_BBN_;

  /** 7) Recombination algorithm */
  pth->recombination=hyrec;
  pth->recfast_photoion_mode=recfast_photoion_Tmat;

  /** 8) Parametrization of reionization */
  pth->reio_parametrization=reio_camb;
  /** 8.a) 'reio_camb' or 'reio_half_tanh' case */
  pth->reio_z_or_tau=reio_z;
  pth->z_reio=7.6711;
  pth->tau_reio=0.05430842;
  pth->reionization_exponent=1.5;
  pth->reionization_width=0.5;
  pth->helium_fullreio_redshift=3.5;
  pth->helium_fullreio_width=0.5;

  /** 8.b) 'reio_bins_tanh' case */
  pth->binned_reio_num=0;
  pth->binned_reio_z=NULL;
  pth->binned_reio_xe=NULL;
  pth->binned_reio_step_sharpness = 0.3;
  /** 8.c) 'reio_many_tanh' case */
  pth->many_tanh_num=0;
  pth->many_tanh_z=NULL;
  pth->many_tanh_xe=NULL;
  pth->many_tanh_width = 0.5;
  /** 8.d) 'reio_inter' case */
  pth->reio_inter_num = 0;
  pth->reio_inter_z = NULL;
  pth->reio_inter_xe = NULL;

  /** 9) Damping scale */
  pth->compute_damping_scale = _FALSE_;

  /**
   * Default to input_read_parameters_species
   */

  /** 1) Photon density */
  pba->T_cmb = 2.7255;
  pba->Omega0_g = (4.*sigma_B/_c_*pow(pba->T_cmb,4.)) / (3.*_c_*_c_*1.e10*pba->h*pba->h/_Mpc_over_m_/_Mpc_over_m_/8./_PI_/_G_);

  /** 2) Baryon density */
  pba->Omega0_b = 0.02238280/pow(pba->h,2);

  /** 3) Ultra-relativistic species / massless neutrino density,
      assuming as default value N_eff=3.044 (see 2008.01074 and
      2012.02726. This value is more accurate than the previous
      default value of 3.046) */
  pba->Omega0_ur = 3.044*7./8.*pow(4./11.,4./3.)*pba->Omega0_g;

  /** 3.a) Effective squared sound speed and viscosity parameter */
  ppt->three_ceff2_ur=1.;
  ppt->three_cvis2_ur=1.;

  /** 4) CDM density */
  pba->Omega0_cdm = 0.1201075/pow(pba->h,2);

  /** 5) ncdm sector */
  /** 5.a) Number of distinct species */
  pba->N_ncdm = 0;
  /** 5.b) List of names of psd files */
  pba->ncdm_psd_files = NULL;
  /** 5.c) Analytic distribution function */
  pba->ncdm_psd_parameters = NULL;
  pba->Omega0_ncdm_tot = 0.;
  /** 5.d) --> See read_parameters_background */
  /** 5.e) ncdm temperature */
  pba->T_ncdm_default = 0.71611; /* this value gives m/omega = 93.14 eV b*/
  pba->T_ncdm = NULL;
  /** 5.f) ncdm chemical potential */
  pba->ksi_ncdm_default = 0.;
  pba->ksi_ncdm = NULL;
  /** 5.g) ncdm degeneracy parameter */
  pba->deg_ncdm_default = 1.;
  pba->deg_ncdm = NULL;
  /** 5.h) --> See read_parameters_background */

  /** 6) Curvature density */
  pba->Omega0_k = 0.;
  pba->K = 0.;
  pba->sgnK = 0;

  /* ** ADDITIONAL SPECIES ** */

  /** 7.1) Decaying CDM into Dark Radiation = dcdm+dr */
  /** 7.1.a) Current fractional density of dcdm+dr */
  pba->Omega0_dcdmdr = 0.0;
  pba->Omega0_dcdm = 0.0;
  /** 7.1.c) Decay constant */
  pba->Gamma_dcdm = 0.0;
  pba->tau_dcdm = 0.0;

  /* ** ADDITIONAL SPECIES ** --> Add your species here */
  /** 7.2) Interacting Dark Matter */
  /** 7.2.a) Current fractional density of idm_dr+idr */
  pba->Omega0_idr = 0.0;
  pba->Omega0_idm_dr = 0.0;
  /** 7.2.b) Current temperature of idm_dr+idr */
  pba->T_idr = 0.0;
  /** 7.2.c) ETHOS parameters of idm_dr+idr */
  pth->a_idm_dr = 0.;
  pth->b_idr = 0.;
  pth->nindex_idm_dr = 4.;
  pth->m_idm_dr = 1.e11;
  /** 7.2.d) Approximation mode of idr */
  ppt->idr_nature=idr_free_streaming;

  /* ** ADDITIONAL SPECIES ** */

  /** 9) Dark energy contributions */
  pba->Omega0_fld = 0.;
  pba->Omega0_scf = 0.;
  pba->Omega0_lambda = 1.-pba->Omega0_k-pba->Omega0_g-pba->Omega0_ur-pba->Omega0_b-pba->Omega0_cdm-pba->Omega0_ncdm_tot-pba->Omega0_dcdmdr+pba->Omega0_idr+pba->Omega0_idm_dr;
  /** 8.a) Omega fluid */
  /** 8.a.1) PPF approximation */
  pba->use_ppf = _TRUE_;
  pba->c_gamma_over_c_fld = 0.4;
  /** 9.a.2) Equation of state */
  pba->fluid_equation_of_state = CLP;
  pba->w0_fld = -1.;
  pba->cs2_fld = 1.;
  /** 9.a.2.1) 'CLP' case */
  pba->wa_fld = 0.;
  /** 9.a.2.2) 'EDE' case */
  pba->Omega_EDE = 0.;
  /** 9.b) Omega scalar field */
  /** 9.b.1) Potential parameters and initial conditions */
  pba->scf_parameters = NULL;
  pba->scf_parameters_size = 0;
  /** 9.b.2) Initial conditions from attractor solution */
  pba->attractor_ic_scf = _TRUE_;
  pba->phi_ini_scf = 1;                // MZ: initial conditions are as multiplicative
  pba->phi_prime_ini_scf = 1;          //     factors of the radiation attractor values
  /** 9.b.3) Tuning parameter */
  pba->scf_tuning_index = 0;
  /** 9.b.4) Shooting parameter */
  pba->shooting_failed = _FALSE_;

  /**
   * Deafult to input_read_parameters_heating
   */
  pth->has_exotic_injection = _FALSE_;

  /** 1) DM annihilation */
  /** 1.a) Energy fraction absorbed by the gas */
  pin->DM_annihilation_efficiency = 0.;
  pin->DM_annihilation_cross_section = 0.;
  pin->DM_annihilation_mass = 0.;
  pin->DM_annihilation_fraction = 0.;
  /** 1.a.1) Redshift dependence */
  pin->DM_annihilation_variation = 0.;
  pin->DM_annihilation_z = 1000.;
  pin->DM_annihilation_zmax = 2500.;
  pin->DM_annihilation_zmin = 30.;
  pin->DM_annihilation_f_halo = 0.;
  pin->DM_annihilation_z_halo = 30.;

  /** 2) DM decay */
  /** 2.a) Fraction */
  pin->DM_decay_fraction = 0.;
  /** 2.b) Decay width */
  pin->DM_decay_Gamma = 0.;

  /** 3) PBH evaporation */
  /** 3.a) Fraction */
  pin->PBH_evaporation_fraction = 0.;
  /** 3.b) Mass */
  pin->PBH_evaporation_mass = 0.;

  /** 4) PBH accretion */
  /** 4.a) Fraction */
  pin->PBH_accretion_fraction = 0.;
  /** 4.b) Mass */
  pin->PBH_accretion_mass = 0.;
  /** 4.c) Recipe */
  pin->PBH_accretion_recipe = disk_accretion;
  /** 4.c.1) Additional parameters for spherical accretion */
  pin->PBH_accretion_relative_velocities = -1.;
  /** 4.c.1) Additional parameters for disk accretion */
  pin->PBH_accretion_eigenvalue = 0.1;
  pin->PBH_accretion_ADAF_delta = 1.e-3;

  /** 5) Injection efficiency */
  pin->f_eff_type = f_eff_on_the_spot;
  pin->f_eff = 1.;
  sprintf(pin->f_eff_file,"/external/heating/example_f_eff_file.dat");

  /** 6) Deposition function */
  pin->chi_type = chi_CK;
  /** 6.1) External file */
  sprintf(pin->chi_z_file,"/external/heating/example_chiz_file.dat");
  sprintf(pin->chi_x_file,"/external/heating/example_chix_file.dat");

  /**
   * Default to input_read_parameters_nonlinear
   */

  /** 1) Non-linearity */
  ppt->has_nl_corrections_based_on_delta_m = _FALSE_;
  pfo->method = nl_none;
  pfo->has_pk_eq = _FALSE_;
  pfo->extrapolation_method = extrap_max_scaled;
  pfo->feedback = nl_emu_dmonly;
  pfo->z_infinity = 10.;

  /**
   * Default to input_read_parameters_primordial
   */

  /** 1) Primordial spectrum type */
  ppm->primordial_spec_type = analytic_Pk;
  /** 1.a) Pivot scale in Mpc-1 */
  ppm->k_pivot = 0.05;

  /** 1.b) For type 'analytic_Pk' */
  /** 1.b.1) For scalar perturbations */
  ppm->A_s = 2.100549e-09;
  /** 1.b.1.1) Adiabatic perturbations */
  ppm->n_s = 0.9660499;
  ppm->alpha_s = 0.;
  /** 1.b.1.2) Isocurvature/entropy perturbations */
  ppm->f_bi = 1.;
  ppm->n_bi = 1.;
  ppm->alpha_bi = 0.;
  ppm->f_cdi = 1.;
  ppm->n_cdi = 1.;
  ppm->alpha_cdi = 0.;
  ppm->f_nid = 1.;
  ppm->n_nid = 1.;
  ppm->alpha_nid = 0.;
  ppm->f_niv = 1.;
  ppm->n_niv = 1.;
  ppm->alpha_niv = 0.;
  /** 1.b.1.3) Cross-correlation between different adiabatic/entropy mode */
  ppm->c_ad_bi = 0.;
  ppm->n_ad_bi = 0.;
  ppm->alpha_ad_bi = 0.;
  ppm->c_ad_cdi = 0.;
  ppm->n_ad_cdi = 0.;
  ppm->alpha_ad_cdi = 0.;
  ppm->c_ad_nid = 0.;
  ppm->n_ad_nid = 0.;
  ppm->alpha_ad_nid = 0.;
  ppm->c_ad_niv = 0.;
  ppm->n_ad_niv = 0.;
  ppm->alpha_ad_niv = 0.;
  ppm->c_bi_cdi = 0.;
  ppm->n_bi_cdi = 0.;
  ppm->alpha_bi_cdi = 0.;
  ppm->c_bi_nid = 0.;
  ppm->n_bi_nid = 0.;
  ppm->alpha_bi_nid = 0.;
  ppm->c_bi_niv = 0.;
  ppm->n_bi_niv = 0.;
  ppm->alpha_bi_niv = 0.;
  ppm->c_cdi_nid = 0.;
  ppm->n_cdi_nid = 0.;
  ppm->alpha_cdi_nid = 0.;
  ppm->c_cdi_niv = 0.;
  ppm->n_cdi_niv = 0.;
  ppm->alpha_cdi_niv = 0.;
  ppm->c_nid_niv = 0.;
  ppm->n_nid_niv = 0.;
  ppm->alpha_nid_niv = 0.;
  /** 1.b.2) For tensor perturbations */
  ppm->r = 1.;
  ppm->n_t = -ppm->r/8.*(2.-ppm->r/8.-ppm->n_s);
  ppm->alpha_t = ppm->r/8.*(ppm->r/8.+ppm->n_s-1.);
  /** 1.c) For type 'inflation_V' */
  /** 1.c.2) Coefficients of the Taylor expansion */
  ppm->V0=1.25e-13;
  ppm->V1=-1.12e-14;
  ppm->V2=-6.95e-14;
  ppm->V3=0.;
  ppm->V4=0.;
  /** 1.d) For type 'inflation_H' */
  ppm->H0=3.69e-6;
  ppm->H1=-5.84e-7;
  ppm->H2=0.;
  ppm->H3=0.;
  ppm->H4=0.;
  /** 1.e) For type 'inflation_V_end' */
  /** 1.e.1) Value of the field at the minimum of the potential */
  ppm->phi_end=0.;
  /** 1.e.2) Shape of the potential */
  ppm->potential=polynomial;
  /** 1.e.4) Increase of scale factor or (aH) between Hubble crossing at pivot
             scale and end of inflation */
  ppm->phi_pivot_method = N_star;
  ppm->phi_pivot_target = 60;
  /** 1.e.5) Nomral numerical integration or analytical slow-roll formulas? */
  ppm->behavior=numerical;
  /** 1.g) For type 'external_Pk' */
  /** 1.g.1) Command generating the table */
  ppm->command=NULL;//"write here your command for the external Pk"

  /** 1.g.2) Parameters to be passed to the command */
  ppm->custom1=0.;
  ppm->custom2=0.;
  ppm->custom3=0.;
  ppm->custom4=0.;
  ppm->custom5=0.;
  ppm->custom6=0.;
  ppm->custom7=0.;
  ppm->custom8=0.;
  ppm->custom9=0.;
  ppm->custom10=0.;

  /**
   * Default to input_read_parameters_spectra
   */

  /** 1) Maximum l for CLs */
  ppt->l_scalar_max=2500;
  ppt->l_vector_max=500;
  ppt->l_tensor_max=500;
  ppt->l_lss_max=300;

  /** 2) Parameters for the the matter density number count */
  /** 2.a) Selection functions W(z) of each redshift bin */
  ppt->selection=gaussian;
  ppt->selection_num=1;
  ppt->selection_mean[0]=1.;
  ppt->selection_width[0]=0.1;
  ptr->selection_bias[0]=1.;
  ptr->selection_magnification_bias[0]=0.;
  phr->non_diag=0;
  /** 2.b) Selection function */
  ptr->has_nz_analytic = _FALSE_;
  ptr->has_nz_file = _FALSE_;
  /** 2.c) Source number counts evolution */
  ptr->has_nz_evo_analytic = _FALSE_;
  ptr->has_nz_evo_file = _FALSE_;

  /** 3) Power spectrum P(k) */
  /** 3.a) Maximum k in P(k) */
  ppt->k_max_for_pk=1.;
  /** 3.a) Maximum k in P(k) primordial */
  ppm->has_k_max_for_primordial_pk = _FALSE_;
  /** 3.b) Redshift values */
  pop->z_pk_num = 1;
  pop->z_pk[0] = 0.;
  /** 3.c) Maximum redshift */
  ppt->z_max_pk=0.;

  /**
   * Default to input_read_parameters_lensing
   */

  /** 1) Lensing */
  ple->has_lensed_cls = _FALSE_;

  /** 2) Should the lensed spectra be rescaled? */
  ptr->lcmb_rescale=1.;
  ptr->lcmb_tilt=0.;
  ptr->lcmb_pivot=0.1;

  /**
   * Default to input_read_parameters_distortions
   */

  /** 1) Branching ratio approximation */
  psd->sd_branching_approx = bra_exact;
  /** 1.a.1) Number of multipoles in PCA expansion */
  psd->sd_PCA_size=2;
  /** 1.a.2) Detector noise file name */
  psd->has_detector_file = _FALSE_;
  /** 1.a.3) Detector name */
  psd->has_user_defined_name = _FALSE_;
  psd->has_user_defined_detector = _FALSE_;
  sprintf(psd->sd_detector_name,"PIXIE");
  /** 1.3.a.1) Detector nu min */
  psd->sd_detector_nu_min = 30.;
  /** 1.3.a.2) Detector nu max */
  psd->sd_detector_nu_max = 1005.;
  /** 1.3.a.3) Detector nu delta/bin number */
  psd->sd_detector_nu_delta = 15.;
  psd->sd_detector_bin_number = 65;
  /** 1.3.a.1) Detector noise */
  psd->sd_detector_delta_Ic = 5.e-26;

  /** 2) Only exotic species? */
  psd->include_only_exotic = _FALSE_;

  /** 3) Include g distortion in total calculation? */
  psd->include_g_distortion = _FALSE_;

  /** 4) Additional y or mu parameters? */
  psd->sd_add_y = 0.;
  psd->sd_add_mu = 0.;

  /** 5) Include SZ effect from reionization? */
  psd->has_SZ_effect = _FALSE_;
  /** 5.a) What type of approximation you want to use for the SZ effect? */
  psd->sd_reio_type = sd_reio_Chluba;

  /**
   * Default to input_read_additional
   */

  pth->compute_cb2_derivatives=_FALSE_;

  /**
   * Default to input_read_parameters_output
   */

  /** 1) Output for external files */
  /** 1.a) File name */
  sprintf(pop->root,"output/");
  /** 1.b) Headers */
  pop->write_header = _TRUE_;
  /** 1.c) Format */
  pop->output_format = class_format;
  /** 1.d) Background quantities */
  pop->write_background = _FALSE_;
  /** 1.e) Thermodynamics quantities */
  pop->write_thermodynamics = _FALSE_;
  /** 1.f) Table of perturbations for certain wavenumbers k */
  ppt->k_output_values_num=0;
  pop->write_perturbations = _FALSE_;
  ppt->store_perturbations = _FALSE_;
  /** 1.g) Primordial spectra */
  pop->write_primordial = _FALSE_;
  /** 1.h) Exotic energy injection function */
  pop->write_exotic_injection = _FALSE_;
  pop->write_noninjection = _FALSE_;
  /** 1.i) Spectral distortions */
  pop->write_distortions = _FALSE_;


  /** 2) Verbosity */
  pba->background_verbose = 0;
  pth->thermodynamics_verbose = 0;
  pth->hyrec_verbose = 0;
  ppt->perturbations_verbose = 0;
  ptr->transfer_verbose = 0;
  ppm->primordial_verbose = 0;
  phr->harmonic_verbose = 0;
  pfo->fourier_verbose = 0;
  ple->lensing_verbose = 0;
  psd->distortions_verbose = 0;
  pop->output_verbose = 0;

  return _SUCCESS_;

}
