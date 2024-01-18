/** @file class.c
 * Julien Lesgourgues, 17.04.2011
 */

/* this main calls CLASS several times in a loop, with different input
   parameters. It illustrates how the code could be interfaced with a
   parameter extraction code. */

/* JL 17.03.2016: implemented here nested openMP loops. The user
   chooses how many instances of CLASS are run in parallel (by setting
   the variable number_of_class_instances below). Each of them uses a
   number of thread such that all cores are used. */

#include "class.h"

int class(
          struct file_content *pfc,
          struct precision * ppr,
          struct background * pba,
          struct thermodynamics * pth,
          struct perturbations * ppt,
          struct primordial * ppm,
          struct fourier * pfo,
          struct transfer * ptr,
          struct harmonic * phr,
          struct lensing * ple,
          struct distortions * psd,
          struct output * pop,
          int l_max,
          double ** cl,
          ErrorMsg errmsg) {

  int l;

  class_call(input_read_from_file(pfc,ppr,pba,pth,ppt,ptr,ppm,phr,pfo,ple,psd,pop,errmsg),
             errmsg,
             errmsg);

  class_call(background_init(ppr,pba),
             pba->error_message,
             errmsg);

  class_call(thermodynamics_init(ppr,pba,pth),
             pth->error_message,
             errmsg);

  class_call(perturbations_init(ppr,pba,pth,ppt),
             ppt->error_message,
             errmsg);

  class_call(primordial_init(ppr,ppt,ppm),
             ppm->error_message,
             errmsg);

  class_call(fourier_init(ppr,pba,pth,ppt,ppm,pfo),
             pfo->error_message,
             errmsg);

  class_call(transfer_init(ppr,pba,pth,ppt,pfo,ptr),
             ptr->error_message,
             errmsg);

  class_call(harmonic_init(ppr,pba,ppt,ppm,pfo,ptr,phr),
             phr->error_message,
             errmsg);

  class_call(lensing_init(ppr,ppt,phr,pfo,ple),
             ple->error_message,
             errmsg);

  /****** write the Cl values in the input array cl[l]  *******/

  for (l=2; l <= l_max; l++) {

    class_call(output_total_cl_at_l(phr,ple,pop,(double)l,cl[l]),
               phr->error_message,
               errmsg);
  }

  /****** all calculations done, now free the structures ******/

  class_call(lensing_free(ple),
             ple->error_message,
             errmsg);

  class_call(harmonic_free(phr),
             phr->error_message,
             errmsg);

  class_call(transfer_free(ptr),
             ptr->error_message,
             errmsg);

  class_call(fourier_free(pfo),
             pfo->error_message,
             errmsg);

  class_call(primordial_free(ppm),
             ppm->error_message,
             errmsg);

  class_call(perturbations_free(ppt),
             ppt->error_message,
             errmsg);

  class_call(thermodynamics_free(pth),
             pth->error_message,
             errmsg);

  class_call(background_free(pba),
             pba->error_message,
             errmsg);

  return _SUCCESS_;

}

int main() {

  /* shared variable that will be common to all CLASS instances */
  int i;
  int l,l_max;
  int num_ct_max=7;
  int num_loops=10;

  struct file_content fc;
  ErrorMsg errmsg_parser;

  int total_number_of_threads;
  int number_of_class_instances;
  int number_of_threads_inside_class;

  int index_ct_tt;
  int index_ct_ee;
  int index_ct_te;

  /* dealing with the openMP part (number of instances, number of
     threads per instance...) */

#ifdef _OPENMP

  /* Determine the number of threads, class instances and nested
     threads */
  total_number_of_threads = omp_get_max_threads();

  /* User-fixed number of CLASS instances to be run in parallel (Total
     number of threads should be dividable by this number) */
  number_of_class_instances=2;

  if ((total_number_of_threads % number_of_class_instances) != 0)
    printf("The total number of threads, %d, is not a mutiple of the requested number of CLASS instances, %d\n",total_number_of_threads,number_of_class_instances);
  number_of_threads_inside_class = total_number_of_threads/number_of_class_instances;

  /* inferred number of threads per instance */
  printf("# Total number of available threads = %d, used to run\n",
         total_number_of_threads);
  printf("# -> %d CLASS executables in parallel\n",
         number_of_class_instances);
  printf("# -> %d threads inside each CLASS executables\n",
         number_of_threads_inside_class);

  /* Turn on nested parallelism */
  omp_set_nested(1);
#endif

  /* choose a value of l_max in C_l's */
  l_max=3000;

  /* all parameters for which we don't want to keep default values
     should be passed to the code through a file_content
     structure. Create such a structure with the size you need: 10 in
     this exemple */
  parser_init(&fc,10,"",errmsg_parser);

  /* assign values to these 9 parameters. Some will be fixed, some
     will be varied in the loop. */
  strcpy(fc.name[0],"output");
  strcpy(fc.value[0],"tCl,pCl,lCl");

  strcpy(fc.name[1],"l_max_scalars");
  sprintf(fc.value[1],"%d",l_max);

  strcpy(fc.name[2],"lensing");
  sprintf(fc.value[2],"yes");

  strcpy(fc.name[3],"H0");
  sprintf(fc.value[3],"%e",72.);

  strcpy(fc.name[4],"omega_b");
  sprintf(fc.value[4],"%e",0.024);

  strcpy(fc.name[5],"omega_cdm");
  sprintf(fc.value[5],"%e",0.05);

  strcpy(fc.name[6],"z_reio");
  sprintf(fc.value[6],"%e",10.);

  strcpy(fc.name[7],"A_s");
  sprintf(fc.value[7],"%e",2.3e-9);

  strcpy(fc.name[8],"n_s");
  sprintf(fc.value[8],"%e",1.);

  strcpy(fc.name[9],"perturbations_verbose");
  sprintf(fc.value[9],"%d",0); // Trick: set to 2 to cross-check actual number of threads per CLASS instance

  /* Create an array of Cl's where all results will be stored for each parameter value in the loop */
  double *** cl;
  cl = malloc(num_loops*sizeof(double**));

  /* Create one thread for each instance of CLASS */
#pragma omp parallel num_threads(number_of_class_instances)
  {

    /* set the number of threads inside each CLASS instance */
#ifdef _OPENMP
    omp_set_num_threads(number_of_threads_inside_class);
#endif

    /* for each thread/instance, create all CLASS input/output
       structures (these variables are being declared insode the
       parallel zone, hence they are openMP private variables) */

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
    ErrorMsg errmsg;            /* for error messages */

    struct file_content fc_local;
    int j,iam;

    /* copy the shared file content into the local file content used by each instance */
    parser_init(&fc_local,fc.size,"",errmsg);
    for (j=0; j < fc.size; j++) {
      strcpy(fc_local.value[j],fc.value[j]);
      strcpy(fc_local.name[j],fc.name[j]);
      fc_local.read[j]=fc.read[j];
    }

    /* loop over (num_loops) values of some parameters: in this exemple, omega_b */
#pragma omp for schedule(static,1)
    for (i=0; i<=num_loops; i++) {

#ifdef _OPENMP
      iam=omp_get_thread_num();
#else
      iam=0;
#endif

      /* assign one value to omega_b */
      sprintf(fc_local.value[4],"%e",0.01+i*0.002);

      printf("# %d\tthread=%d : running with omega_b = %s\n",i,iam,fc_local.value[4]);

      /* allocate the array where the Cl's calculated by one instance
         will be written (we could add another array with P(k), or
         extract other results from the code - here we assume that we
         are interested in the C_l's only */
      cl[i]=malloc((l_max+1)*sizeof(double*));
      for (l=0;l<=l_max;l++)
        cl[i][l]=malloc(num_ct_max*sizeof(double));

      /* calls class and return the C_l's*/
      if (class(&fc_local,&pr,&ba,&th,&pt,&pm,&fo,&tr,&hr,&le,&sd,&op,l_max,cl[i],errmsg) == _FAILURE_) {
        printf("\n\nError in class \n=>%s\n",errmsg);
        //return _FAILURE_;
      }

      /* if this is the first call, extract dynamically the value of indices used in the output */
      if ((i==0) && (iam==0)) {
        index_ct_tt=hr.index_ct_tt;
        index_ct_te=hr.index_ct_te;
        index_ct_ee=hr.index_ct_ee;
      }

    } // end of loop over parameters
  } // end parallel zone

  /* write in file the lensed C_l^TT, C_l^EE, C_l^TE's obtained in all runs */

  FILE * out=fopen("output/test_loops_omp.dat","w");

  for (i=0; i<num_loops; i++) {
    for (l=2;l<=l_max;l++) {
      fprintf(out,"%d  %e  %e  %e\n",
              l,
              l*(l+1)*cl[i][l][index_ct_tt],
              l*(l+1)*cl[i][l][index_ct_ee],
              l*(l+1)*cl[i][l][index_ct_te]);
    }
    fprintf(out,"\n");
  }

  /* free Cl's array */
  for (i=0; i<num_loops; i++) {
    for (l=0;l<l_max;l++) {
      free(cl[i][l]);
    }
    free(cl[i]);
  }
  free(cl);

  return _SUCCESS_;

}
