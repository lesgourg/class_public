/** @file distortions.c Documented module on spectral distortions
 * Matteo Lucca, 31.10.2018
 * Nils Schoeneberg, 18.02.2019
 *
 * When using this module please consider citing:
 * Lucca et al. 2019 (JCAP02(2020)026, arXiv:1910.04619)
 * as well as related pioneering works such as:
 * Chluba & Sunyaev 2012 (MNRAS419(2012)1294-1314, arXiv:1109.6552)
 * Chluba 2013 (MNRAS434(2013)352, arXiv:1304.6120)
 * Clube & Jeong 2014 (MNRAS438(2014)2065–2082, arXiv:1306.5751)
 */

#include "distortions.h"

/**
 * Initialize the distortions structure.
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param ppt        Input: pointer to the perturbations structure
 * @param ppm        Input: pointer to the primordial structure
 * @param psd        Input/Output: pointer to initialized distortions structure
 * @return the error status
 */

int distortions_init(struct precision * ppr,
                     struct background * pba,
                     struct thermodynamics * pth,
                     struct perturbations * ppt,
                     struct primordial * ppm,
                     struct distortions * psd) {

  if(psd->has_distortions == _FALSE_){
    if (psd->distortions_verbose > 0)
      printf("No distortions requested. Distortions module skipped.\n");
    return _SUCCESS_;
  }
  if (psd->distortions_verbose > 0) {
    printf("Computing spectral distortions\n");
  }

  class_test(pth->compute_damping_scale==_FALSE_,
             psd->error_message,
             "Cannot compute spectral distortions without damping scale\n");

  /** Set physical constants */
  class_call(distortions_constants(ppr,pba,pth,psd),
             psd->error_message,
             psd->error_message);

  if(psd->sd_branching_approx == bra_exact){
    /** Set/Check the distortions detector */
    class_call(distortions_set_detector(ppr,psd),
               psd->error_message,
               psd->error_message);
  }

  /** Assign values to all indices in the distortions structure */
  class_call(distortions_indices(psd),
             psd->error_message,
             psd->error_message);

  /** Define z and x arrays */
  class_call(distortions_get_xz_lists(ppr,pba,pth,psd),
             psd->error_message,
             psd->error_message);

  /** Define branching ratios */
  class_call(distortions_compute_branching_ratios(ppr,psd),
             psd->error_message,
             psd->error_message);

  /** Define heating function */
  class_call(distortions_compute_heating_rate(ppr,pba,pth,ppt,ppm,psd),
             psd->error_message,
             psd->error_message);

  /** Define final spectral distortions */
  class_call(distortions_compute_spectral_shapes(ppr,pba,pth,psd),
             psd->error_message,
             psd->error_message);

  return _SUCCESS_;
}

/**
 * Free all memory space allocated by distortions_init()
 *
 * @param psd     Input: pointer to distortions structure (to be freed)
 * @return the error status
 */

int distortions_free(struct distortions * psd) {

  /** Define local variables */
  int index_type;

  if(psd->has_distortions == _TRUE_){
    /** Delete lists */
    free(psd->z);
    free(psd->z_weights);
    free(psd->x);
    free(psd->x_weights);

    /** Delete noise file */
    if(psd->has_detector_file == _TRUE_){
      free(psd->delta_Ic_array);
    }

    /** Delete branching ratios */
    for(index_type=0;index_type<psd->type_size;++index_type){
      free(psd->br_table[index_type]);
    }
    free(psd->br_table);

    /** Delete heating functions */
    free(psd->dQrho_dz_tot);

    /** Delete distortion shapes */
    for(index_type=0;index_type<psd->type_size;++index_type){
      free(psd->sd_shape_table[index_type]);
      free(psd->sd_table[index_type]);
    }
    free(psd->sd_shape_table);
    free(psd->sd_table);

    /** Delete distortion amplitudes */
    free(psd->sd_parameter_table);

    /** Delete total distortion */
    free(psd->DI);
  }

  return _SUCCESS_;
}

/**
 * Calculate physical constant.
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to thermodynamics structure
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */

int distortions_constants(struct precision * ppr,
                          struct background * pba,
                          struct thermodynamics * pth,
                          struct distortions * psd){

  /** Define unit conventions */
  psd->x_to_nu = (_k_B_*pba->T_cmb/_h_P_)/1e9;                    // [GHz]
  psd->DI_units = 2.*pow(_k_B_*pba->T_cmb,3.)/pow(_h_P_*_c_,2.);  // [W/(m^2 Hz sr)]

  /** Define transition redshifts z_muy and z_th */
  psd->z_muy = 5.e4;
  psd->z_th = 1.98e6*
         pow((1.-pth->YHe/2.)/0.8767,-2./5.)*
         pow(pba->Omega0_b*pow(pba->h,2.)/0.02225,-2./5.)*
         pow(pba->T_cmb/2.726,1./5.);

  sprintf(psd->sd_PCA_file_generator,"%s/%s",ppr->sd_external_path,"generate_PCA_files.py");
  sprintf(psd->sd_detector_list_file,"%s/%s",ppr->sd_external_path,"detectors_list.dat");

  return _SUCCESS_;
}

/**
 * Check wether the detector name and the detector properties
 * are a valid combination.
 *
 * There are four options for the user
 *
 * defined_name = true, defined_detector = true
 * Meaning: The user requests a specific detector with specific settings
 * --> Check that the detector exists and has the same settings
 *
 * defined_name = true, defined_detector = false
 * Meaning: The user requests a specific detector
 * --> Check that the detector exists and use the given settings
 *
 * defined_name = false, defined_detector = true
 * Meaning: The user requests specific settings, but does not name their detector
 * --> Check that the settings exists, or create them
 *
 * defined_name = false, defined_detector = false
 * Meaning: The user just wants the default detector and settings
 * --> Just use the default settings and skip this function
 *
 * @param ppr        Input: pointer to precision structure
 * @param psd        Input/Output: pointer to initialized distortions structure
 * @return the error status
 */

int distortions_set_detector(struct precision * ppr,
                             struct distortions * psd){

  /** Local variables */
  FILE* det_list_file;
  char line[_LINE_LENGTH_MAX_];
  DetectorName detector_name;
  DetectorFileName detector_noise_file_name;
  double nu_min,nu_max,nu_delta,delta_Ic;
  int has_detector_noise_file;
  int N_bins;
  char * left;
  int headlines = 0;
  int found_detector;

  has_detector_noise_file = _FALSE_;

  if(psd->has_user_defined_name == _FALSE_){
    /* The user wants the default */
    if(psd->has_user_defined_detector == _FALSE_ && psd->has_detector_file == _FALSE_){
      if(psd->distortions_verbose > 0){
        printf(" -> Using the default (%s) detector\n",psd->sd_detector_name);
      }
      return _SUCCESS_; // Nothing more to do
    }
    /* The user wants a new detector with specified settings, but without name */
    else{
      /* Generate a custom name for this custom detector, so we can check if it has already been defined */
      if(psd->has_detector_file == _TRUE_){
        sprintf(psd->sd_detector_name,
                "Custom__%.80s__Detector",
                psd->sd_detector_file_name);
      }
      else{
        sprintf(psd->sd_detector_name,
                "Custom__%7.2e_%7.2e_%7.2e_%i_%7.2e__Detector",
                psd->sd_detector_nu_min,psd->sd_detector_nu_max,psd->sd_detector_nu_delta,psd->sd_detector_bin_number,psd->sd_detector_delta_Ic);
      }
    }
  }

  /** Open file */
  class_open(det_list_file, psd->sd_detector_list_file, "r",
             psd->error_message);

  found_detector = _FALSE_;
  while (fgets(line,_LINE_LENGTH_MAX_-1,det_list_file) != NULL) {
    headlines++;

    /* Eliminate blank spaces at beginning of line */
    left=line;
    while (left[0]==' ') {
        left++;
    }
    if (left[0] > 39) {
      if(sscanf(line,"%s %lg %lg %lg %i %lg",detector_name,&nu_min,&nu_max,&nu_delta,&N_bins,&delta_Ic) != 6){
        has_detector_noise_file = _TRUE_;
        if(sscanf(line,"%s %s",detector_name,detector_noise_file_name) != 2){
          class_stop(psd->error_message,
                     "Could not read line %i in file '%s'\n",headlines,psd->sd_detector_list_file);
        }
      }
      else{
        has_detector_noise_file = _FALSE_;
      }

      /* Detector has been found */
      if(strcmp(psd->sd_detector_name,detector_name)==0){
        if(psd->distortions_verbose > 0){
          printf(" -> Found detector %s (user defined = %s)\n",detector_name,(psd->has_user_defined_detector?"TRUE":"FALSE"));
        }
        found_detector = _TRUE_;

        if(has_detector_noise_file == _TRUE_){
          if (psd->distortions_verbose > 1){
            printf(" -> Properties:    Noise file name = %s \n",
                        detector_noise_file_name);
          }
          if(psd->has_detector_file == _TRUE_){
            class_test(strcmp(psd->sd_detector_file_name,detector_noise_file_name) != 0,
                       psd->error_message,
                       "Noise file path (sd_detector_file_name) disagrees between stored detector '%s' and input ->  %s (input) vs %s (stored)",
                       detector_name,psd->sd_detector_file_name,detector_noise_file_name);
          }
          class_test(psd->has_user_defined_detector,
                     psd->error_message,
                     "Detector property type disagrees between stored detector '%s' and input  ->  Userdefined (input) vs Noisefile (stored)",
                     detector_name);
          sprintf(psd->sd_detector_file_name, "%s", detector_noise_file_name);
          psd->has_detector_file = has_detector_noise_file;
        }
        else{
          if (psd->distortions_verbose > 1){
            printf(" -> Properties:    nu_min = %lg    nu_max = %lg    delta_nu = %lg    N_bins = %i    delta_Ic = %lg \n",
                        nu_min, nu_max, nu_delta, N_bins, delta_Ic);
          }
          /* If the user has defined the detector, check that their and our definitions agree */
          if(psd->has_user_defined_detector == _TRUE_){
            class_test(fabs(psd->sd_detector_nu_min-nu_min)>ppr->tol_sd_detector,
                       psd->error_message,
                       "Minimal frequency (sd_detector_nu_min) disagrees between stored detector '%s' and input ->  %.10e (input) vs %.10e (stored)",
                       detector_name,psd->sd_detector_nu_min,nu_min);
            class_test(fabs(psd->sd_detector_nu_max-nu_max)>ppr->tol_sd_detector,
                       psd->error_message,
                       "Maximal frequency (sd_detector_nu_min) disagrees between stored detector '%s' and input ->  %.10e (input) vs %.10e (stored)",
                       detector_name,psd->sd_detector_nu_max,nu_max);
            class_test(fabs(psd->sd_detector_nu_delta-nu_delta)>ppr->tol_sd_detector,
                       psd->error_message,
                       "Delta frequency (sd_detector_nu_delta) disagrees between stored detector '%s' and input ->  %.10e (input) vs %.10e (stored)",
                       detector_name,psd->sd_detector_nu_delta,nu_delta);
            class_test(fabs(psd->sd_detector_bin_number-N_bins)>ppr->tol_sd_detector,
                       psd->error_message,
                       "Number of bins (sd_detector_bin_number) disagrees between stored detector '%s' and input ->  %i (input) vs %i (stored)",
                       detector_name,psd->sd_detector_bin_number,N_bins);
            class_test(fabs(psd->sd_detector_delta_Ic-delta_Ic)>ppr->tol_sd_detector,
                       psd->error_message,
                       "Detector accuracy (sd_detector_delta_Ic) disagrees between stored detector '%s' and input ->  %.10e (input) vs %.10e (stored)",
                       detector_name,psd->sd_detector_delta_Ic,delta_Ic);
          }
          class_test(psd->has_detector_file,
                     psd->error_message,
                     "Detector property type disagrees between stored detector '%s' and input  ->  Noisefile (input) vs Userdefined (stored)",
                     detector_name);

          /* In any case, just take the detector definition from the file */
          psd->sd_detector_nu_min = nu_min;
          psd->sd_detector_nu_max = nu_max;
          psd->sd_detector_nu_delta = nu_delta;
          psd->sd_detector_bin_number = N_bins;
          psd->sd_detector_delta_Ic = delta_Ic;
        }
      }
    }
  }

  fclose(det_list_file);

  /* If the detector has not been found, either the user has specified the settings and we create a new one,
   * or the user hasn't specified the settings and we have to stop */
  if(found_detector == _FALSE_){
    if(psd->has_user_defined_detector==_TRUE_ || psd->has_detector_file==_TRUE_){
      if(psd->distortions_verbose > 0){
        printf(" -> Generating detector '%s' \n",psd->sd_detector_name);
      }
      class_call(distortions_generate_detector(ppr,psd),
                 psd->error_message,
                 psd->error_message);
    }
    else{
      class_stop(psd->error_message,
                 "You asked for detector '%s', but it was not in the database '%s'.\nPlease check the name of your detector, or specify its properties if you want to create a new one",
                 psd->sd_detector_name,
                 psd->sd_detector_list_file);
    }
  }

  return _SUCCESS_;
}

/**
 * Evaluate branching ratios, spectral shapes, E and S vectors for a given detector as
 * described in external/distortions/README using generate_PCA_files.py.
 *
 * @param ppr        Input: pointer to precision structure
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */

int distortions_generate_detector(struct precision * ppr,
                                  struct distortions * psd){

  /** Define local variables*/
  int is_success;
  char temporary_string[2*_FILENAMESIZE_+2*_MAX_DETECTOR_NAME_LENGTH_+1024];


  /* Test first whether or not python exists*/
  if(psd->distortions_verbose > 0){
    printf(" -> Testing python\n");
  }
  is_success = system("python --version");
  class_test(is_success == -1,
             psd->error_message,
             "The command 'python --version' failed.\nPlease install a valid version of python.");

  /* Then activate the PCA generator*/
  if(psd->distortions_verbose > 0){
    printf(" -> Executing the PCA generator\n");
  }

  if(psd->has_detector_file == _TRUE_){
    sprintf(temporary_string,"python %s %s %s %s %.10e %.10e %i %i %.10e %.10e %.10e",
            psd->sd_PCA_file_generator,
            psd->sd_detector_name,
            ppr->sd_external_path,
            psd->sd_detector_file_name,
            ppr->sd_z_min,
            ppr->sd_z_max,
            ppr->sd_z_size,
            6,
            psd->z_th,
            psd->DI_units,
            psd->x_to_nu);

  }
  else{
    sprintf(temporary_string,"python %s %s %.10e %.10e %.10e  %i %.10e %.10e %i %.10e %i %.10e %.10e %.10e",
            psd->sd_PCA_file_generator,
            psd->sd_detector_name,
            psd->sd_detector_nu_min,
            psd->sd_detector_nu_max,
            psd->sd_detector_nu_delta,
            psd->sd_detector_bin_number,
            ppr->sd_z_min,
            ppr->sd_z_max,
            ppr->sd_z_size,
            psd->sd_detector_delta_Ic,
            6,
            psd->z_th,
            psd->DI_units,
            psd->x_to_nu);
  }
  is_success = system(temporary_string);
  class_test(is_success == -1,
             psd->error_message,
             "The command 'python %s' failed.\nPlease make sure the file exists.",psd->sd_PCA_file_generator);

  return _SUCCESS_;
}

/**
 * Assign value to each relevant index in vectors of distortions quantities.
 *
 * @param psd     Input: pointer to distortions structure
 * @return the error status
 */

int distortions_indices(struct distortions * psd) {

  /** Define local variables */
  int index_type = 0;

  /** Define indeces for tables - br_table defined in distortions_compute_branching_ratios,
                                - sd_parameter_table and
                                - sd_table defined in distortions_compute_spectral_shapes */
  class_define_index(psd->index_type_g,_TRUE_,index_type,1);
  class_define_index(psd->index_type_y,_TRUE_,index_type,1);
  class_define_index(psd->index_type_mu,_TRUE_,index_type,1);
  class_define_index(psd->index_type_PCA,_TRUE_,index_type,psd->sd_PCA_size);

  psd->type_size = index_type;

  return _SUCCESS_;
}

/**
 * Calculate redshift and frequency vectors and weights for redshift integral.
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param psd        Input/Output: pointer to initialized distortions structure
 * @return the error status
 */

int distortions_get_xz_lists(struct precision * ppr,
                             struct background * pba,
                             struct thermodynamics * pth,
                             struct distortions * psd){

  /** Define local variables */
  int index_z, index_x;

  /** Define and allocate z array */
  psd->z_min = ppr->sd_z_min;
  psd->z_max = ppr->sd_z_max;
  psd->z_size = ppr->sd_z_size;
  psd->z_delta = (log(psd->z_max)-log(psd->z_min))/psd->z_size;

  class_alloc(psd->z,
              psd->z_size*sizeof(double),
              psd->error_message);

  for (index_z = 0; index_z < psd->z_size; index_z++) {
    psd->z[index_z] = exp(log(psd->z_min+1)+psd->z_delta*index_z);
  }

  /** Define and allocate integrating weights for z array */
  class_alloc(psd->z_weights,
              psd->z_size*sizeof(double),
              psd->error_message);
  class_call(array_trapezoidal_weights(
                    psd->z,
                    psd->z_size,
                    psd->z_weights,
                    psd->error_message),
             psd->error_message,
             psd->error_message);

  /** Define and allocate x array */
  if(psd->sd_branching_approx != bra_exact){
    psd->x_min = ppr->sd_x_min;
    psd->x_max = ppr->sd_x_max;
    psd->x_size = ppr->sd_x_size;
    psd->x_delta = (log(psd->x_max)-log(psd->x_min))/psd->x_size;

    class_alloc(psd->x,
                psd->x_size*sizeof(double),
                psd->error_message);

    for (index_x = 0; index_x<psd->x_size; index_x++) {
      psd->x[index_x] = exp(log(psd->x_min)+psd->x_delta*index_x);
    }

  }
  else if(psd->has_detector_file == _FALSE_){
    psd->x_min = psd->sd_detector_nu_min/psd->x_to_nu;
    psd->x_max = psd->sd_detector_nu_max/psd->x_to_nu;
    psd->x_delta = psd->sd_detector_nu_delta/psd->x_to_nu;
    psd->x_size = psd->sd_detector_bin_number+1;

    class_alloc(psd->x,
                psd->x_size*sizeof(double),
                psd->error_message);

    for (index_x = 0; index_x<psd->x_size; index_x++){
      psd->x[index_x] = psd->x_min+psd->x_delta*index_x;
    }
  }

  /** Define and allocate integrating weights for x array */
  class_alloc(psd->x_weights,
              psd->x_size*sizeof(double),
              psd->error_message);
  class_call(array_trapezoidal_weights(
                    psd->x,
                    psd->x_size,
                    psd->x_weights,
                    psd->error_message),
             psd->error_message,
             psd->error_message);

  return _SUCCESS_;
}

/**
 * Calculate branching ratios.
 *
 * Computing the full evolution of the thermal history of the universe is rather time consuming
 * and mathematically challenging. It is therefore not implemented here. However, there are
 * (at least) 5 levels of possible approximatin to evaluate the SD branching ratios (see also
 * Chluba 2016 for useful discussion)
 *    1) Use a sharp transition at z_mu-y and no distortions before z_th ('branching approx'=sharp_sharp)
 *    2) Use a sharp transition at z_mu-y and a soft transition at z_th ('branching approx'=sharp_soft)
 *    3) Use a soft transition at a_mu-y and z_th as described in Chluba 2013 ('branching approx'=soft_soft)
 *       In this case, the user must be aware that energy conservation is violated and no residuals
 *       are taken into consideration.
 *    4) Use a soft transition at a_mu-y and z_th imposing conservation of energy
 *       ('branching approx'=soft_soft_cons)
 *    5) Use a PCA method as described in Chluba & Jeong 2014 ('branching approx'=exact)
 *       In this case, the definition of the BRs is detector dependent and the user has therefore to
 *       specify the detector type and corresponding characteristics.
 *
 * All quantities are stored in the table br_table.
 *
 * @param ppr        Input: pointer to the precision structure
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */

int distortions_compute_branching_ratios(struct precision * ppr,
                                         struct distortions* psd){

  /** Define local variables */
  int index_z,index_type,index_k;
  double f_g, f_y, f_mu;
  double *f_E;
  double bb_vis;
  int last_index = 0;

  /** Allocate space for branching ratios in br_table */
  class_alloc(psd->br_table,
              psd->type_size*sizeof(double*),
              psd->error_message);
  for(index_type=0; index_type<psd->type_size; ++index_type){
    class_alloc(psd->br_table[index_type],
                psd->z_size*sizeof(double),
                psd->error_message);
  }

  /** Calulate branching ratios */
  if(psd->sd_branching_approx != bra_exact){
    for(index_z=0; index_z<psd->z_size; ++index_z){
      bb_vis = exp(-pow(psd->z[index_z]/psd->z_th,2.5));

      /* 1) Calculate branching ratios using sharp_sharp transition */
      if(psd->sd_branching_approx == bra_sharp_sharp){
        if(psd->z[index_z]>psd->z_th){
          f_g = 1.;
          f_y = 0.;
          f_mu = 0.;
        }
        if(psd->z[index_z]<psd->z_th && psd->z[index_z]>psd->z_muy){
          f_g = 0.;
          f_y = 0.;
          f_mu = 1.;
        }
        if(psd->z[index_z]<psd->z_muy){
          f_g = 0.;
          f_y = 1.;
          f_mu = 0.;
        }
      }

      /* 2) Calculate branching ratios using sharp_soft transition */
      if(psd->sd_branching_approx == bra_sharp_soft){
        f_g = 1.-bb_vis;
        if(psd->z[index_z]>psd->z_muy){
          f_y = 0.;
          f_mu = bb_vis;
        }
        if(psd->z[index_z]<psd->z_muy){
          f_y = 1.;
          f_mu = 0.;
        }
      }

      /* 3) Calculate branching ratios unsing soft_soft transitions */
      if(psd->sd_branching_approx == bra_soft_soft){
        f_g = 1.-bb_vis;
        f_y = 1.0/(1.0+pow((1.0+psd->z[index_z])/(6.0e4),2.58));
        f_mu = bb_vis*(1.-exp(-pow((1.0+psd->z[index_z])/(5.8e4),1.88)));
      }

      /* 4) Calculate branching ratios unsing soft_soft_cons transitions */
      if(psd->sd_branching_approx == bra_soft_soft_cons){
        f_g = 1.-bb_vis;
        f_y = 1.0/(1.0+pow((1.0+psd->z[index_z])/(6.0e4),2.58));
        f_mu = bb_vis*(1.-f_y);
      }

      psd->br_table[psd->index_type_g][index_z] = f_g;
      psd->br_table[psd->index_type_y][index_z] = f_y;
      psd->br_table[psd->index_type_mu][index_z] = f_mu;

    }
  }
  else{
    /* 5) Calculate branching ratios according to Chluba & Jeong 2014 */

    /* Read and spline data from file branching_ratios.dat */
    class_call(distortions_read_br_data(ppr,psd),
               psd->error_message,
               psd->error_message);
    class_call(distortions_spline_br_data(psd),
               psd->error_message,
               psd->error_message);

    /* Allocate local variable */
    class_alloc(f_E,
                psd->sd_PCA_size*sizeof(double),
                psd->error_message);

    /* Interpolate over z */
    for(index_z=0; index_z<psd->z_size; ++index_z){
      class_call(distortions_interpolate_br_data(psd,
                                                 psd->z[index_z],
                                                 &f_g,
                                                 &f_y,
                                                 &f_mu,
                                                 f_E,
                                                 &last_index),
                 psd->error_message,
                 psd->error_message);

      /* Store quantities in the table*/
      psd->br_table[psd->index_type_g][index_z] = f_g;
      psd->br_table[psd->index_type_y][index_z] = f_y;
      psd->br_table[psd->index_type_mu][index_z] = f_mu;
      for(index_k=0; index_k<psd->sd_PCA_size; ++index_k){
        psd->br_table[psd->index_type_PCA+index_k][index_z] = f_E[index_k];
      }

    }

    /* Free space allocated in distortions_read_br_data */
    class_call(distortions_free_br_data(psd),
               psd->error_message,
               psd->error_message);
    free(f_E);

  }

  return _SUCCESS_;
}

/**
 * Import heating rates from heating structure.
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param ppt        Input: pointer to the perturbations structure
 * @param ppm        Input: pointer to the primordial structure
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */

int distortions_compute_heating_rate(struct precision* ppr,
                                     struct background* pba,
                                     struct thermodynamics * pth,
                                     struct perturbations * ppt,
                                     struct primordial * ppm,
                                     struct distortions * psd){

  /** Define local variables */
  struct noninjection* pni = &(psd->ni);
  struct injection* pin = &(pth->in);

  int index_z;
  double tau;
  int last_index_back;
  double *pvecback;
  double heat;
  double H, a, rho_g;

  if (psd->include_only_exotic == _FALSE_) {
    /** Update heating table with second order contributions */
    class_call(noninjection_init(ppr,pba,pth,ppt,ppm,pni),
               pni->error_message,
               psd->error_message);
  }

  /** Allocate space for background vector */
  last_index_back = 0;
  class_alloc(pvecback,
              pba->bg_size*sizeof(double),
              psd->error_message);

  /** Allocate space for total heating function */
  class_alloc(psd->dQrho_dz_tot,
              psd->z_size*sizeof(double*),
              psd->error_message);

  /* Loop over z and calculate the heating at each point */
  for(index_z=0; index_z<psd->z_size; ++index_z){

    /** Import quantities from background structure */
    class_call(background_tau_of_z(pba,
                                   psd->z[index_z],
                                   &tau),
               pba->error_message,
               psd->error_message);
    class_call(background_at_tau(pba,
                                 tau,
                                 long_info,
                                 inter_closeby,
                                 &last_index_back,
                                 pvecback),
               pba->error_message,
               psd->error_message);
    H = pvecback[pba->index_bg_H]*_c_/_Mpc_over_m_;               // [1/s]
    a = pvecback[pba->index_bg_a];                                // [-]
    rho_g = pvecback[pba->index_bg_rho_g]*_Jm3_over_Mpc2_;        // [J/m^3]

    heat = 0;

    /** Import heat from non-injection structure */
    if (psd->include_only_exotic == _FALSE_) {
      class_call(noninjection_photon_heating_at_z(pni,
                                                  psd->z[index_z],
                                                  &heat),           // [J/(m^3 s)]
                 pni->error_message,
                 psd->error_message);
    }

    /** Add heat from injection structure */
    if (pth->has_exotic_injection == _TRUE_) {
      class_call(injection_deposition_at_z(pth,
                                           psd->z[index_z]),
                 pin->error_message,
                 psd->error_message);
      heat += pin->pvecdeposition[pin->index_dep_heat];
    }

    /** Calculate total heating rate */
    psd->dQrho_dz_tot[index_z] = heat*a/(H*rho_g);                // [-]
  }

  free(pvecback);

  if (psd->include_only_exotic == _FALSE_) {
    /** Update heating table with second order contributions */
    class_call(noninjection_free(pni),
               pni->error_message,
               psd->error_message);
  }

  return _SUCCESS_;
}

/**
 * Calculate spectral amplitudes and corresponding distortions.
 *
 * The calculation has been done according to Chluba & Jeong 2014 (arxiv:1306.5751).
 * All quantities are stored in the tables sd_parameter_table and sd_table.
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to thermodynamics structure
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */

int distortions_compute_spectral_shapes(struct precision * ppr,
                                        struct background * pba,
                                        struct thermodynamics * pth,
                                        struct distortions * psd){

  /** Define local variables */
  double * S;
  int last_index = 0;
  int index_type, index_x, index_k;
  double sum_S, sum_G;
  double g;
  double y_reio, DI_reio;

  /** Allocate space for spectral distortion amplitude in table sd_parameter_table */
  class_alloc(psd->sd_parameter_table,
              psd->type_size*sizeof(double),
              psd->error_message);

  /** Compute distortion amplitudes corresponding to each branching ratio (g, y and mu) */
  /* Define y, mu, g and mu_k from heating rates */
  for(index_type=0; index_type<psd->type_size; ++index_type){
    class_call(array_trapezoidal_convolution(psd->dQrho_dz_tot,
                                             psd->br_table[index_type],
                                             psd->z_size,
                                             psd->z_weights,
                                             &(psd->sd_parameter_table[index_type]),
                                             psd->error_message),
               psd->error_message,
               psd->error_message);

    if(index_type>=psd->index_type_PCA){
      /* The S_k are not properly normalized, we have to renormalize here */
      psd->sd_parameter_table[index_type] /= (log(1.+psd->z[1])-log(1.+psd->z[0]));
    }
  }

  psd->sd_parameter_table[psd->index_type_g] /= 4.;
  psd->sd_parameter_table[psd->index_type_y] /= 4.;
  psd->sd_parameter_table[psd->index_type_mu] *= 1.401;

  psd->sd_parameter_table[psd->index_type_y] += psd->sd_add_y;
  psd->sd_parameter_table[psd->index_type_mu] += psd->sd_add_mu;

  /** Allocate space for distortions shapes in distortions_table */
  class_alloc(psd->sd_shape_table,
              psd->type_size*sizeof(double*),
              psd->error_message);
  for(index_type=0; index_type<psd->type_size; ++index_type){
    class_alloc(psd->sd_shape_table[index_type],
                psd->x_size*sizeof(double),
                psd->error_message);
  }

  /** Calculate spectral shapes */
  if(psd->sd_branching_approx != bra_exact || psd->sd_PCA_size == 0){
    /* If no PCA analysis is required, the shapes have simple analistical form */
    for(index_x=0; index_x<psd->x_size; ++index_x){
      psd->sd_shape_table[psd->index_type_g][index_x] = pow(psd->x[index_x],4.)*exp(-psd->x[index_x])/
                                                           pow(1.-exp(-psd->x[index_x]),2.);        // [-]
      psd->sd_shape_table[psd->index_type_y][index_x] = psd->sd_shape_table[psd->index_type_g][index_x]*
                                                           (psd->x[index_x]*(1.+exp(-psd->x[index_x]))/
                                                           (1.-exp(-psd->x[index_x]))-4.);          // [-]
      psd->sd_shape_table[psd->index_type_mu][index_x] = psd->sd_shape_table[psd->index_type_g][index_x]*
                                                           (1./2.19229-1./psd->x[index_x]);         // [-]
    }
  }
  else{
    /* If PCA analysis is required, the shapes has to be vectorized. This is done in the external
       file spectral_shapes.dat using generate_PCA_files.py */

    /* Read and spline data from file spectral_shapes.dat */
    class_call(distortions_read_sd_data(ppr,psd),
               psd->error_message,
               psd->error_message);
    class_call(distortions_spline_sd_data(psd),
               psd->error_message,
               psd->error_message);

    /* Allocate local variable */
    class_alloc(S,
                psd->sd_PCA_size*sizeof(double),
                psd->error_message);

    /* Interpolate over z */
    for(index_x=0; index_x<psd->x_size; ++index_x){
      class_call(distortions_interpolate_sd_data(psd,
                                                 psd->x[index_x]*psd->x_to_nu,
                                                 &psd->sd_shape_table[psd->index_type_g][index_x],
                                                 &psd->sd_shape_table[psd->index_type_y][index_x],
                                                 &psd->sd_shape_table[psd->index_type_mu][index_x],
                                                 S,
                                                 &last_index),
                 psd->error_message,
                 psd->error_message);

    for(index_k=0; index_k<psd->sd_PCA_size; ++index_k){
        psd->sd_shape_table[psd->index_type_PCA+index_k][index_x] = S[index_k];
      }
    }

    /* Free allocated space */
    class_call(distortions_free_sd_data(psd),
               psd->error_message,
               psd->error_message);
    free(S);
  }

  /** Compute distortion amplitude for residual parameter epsilon */
  /* For the details of the calculation see Chluba & Jeong 2014, left column of page 6 */
  psd->epsilon = 0.;

  if(psd->sd_branching_approx == bra_exact && psd->sd_PCA_size != 0){
    class_call(array_trapezoidal_integral(psd->sd_shape_table[psd->index_type_g],
                                          psd->x_size,
                                          psd->x_weights,
                                          &(sum_G),
                                          psd->error_message),
               psd->error_message,
               psd->error_message);
    for(index_k=0; index_k<psd->sd_PCA_size; ++index_k){
      class_call(array_trapezoidal_integral(psd->sd_shape_table[psd->index_type_PCA+index_k],
                                            psd->x_size,
                                            psd->x_weights,
                                            &(sum_S),
                                            psd->error_message),
                 psd->error_message,
                 psd->error_message);
      psd->epsilon += (4.*sum_S/sum_G)*psd->sd_parameter_table[psd->index_type_PCA+index_k];
    }
  }

  /** Allocate space for final spectral distortion */
  class_alloc(psd->DI,
              psd->x_size*sizeof(double),
              psd->error_message);

  class_alloc(psd->sd_table,
              psd->type_size*sizeof(double*),
              psd->error_message);
  for(index_type=0; index_type<psd->type_size; ++index_type){
    class_alloc(psd->sd_table[index_type],
                psd->x_size*sizeof(double),
                psd->error_message);
  }

  /** Calculate spectral distortions according to Chluba & Jeong 2014 (arxiv:1306.5751, Eq. (11)) */
  for(index_x=0;index_x<psd->x_size;++index_x){
    psd->DI[index_x] = 0.;

    for(index_type=0;index_type<psd->type_size;++index_type){
      if(index_type==psd->index_type_g){
        if(psd->include_g_distortion == _TRUE_){
          g = psd->sd_parameter_table[psd->index_type_g];
          psd->sd_table[index_type][index_x] = (1.+g)*g*psd->sd_shape_table[psd->index_type_g][index_x]+
                                                 g*g*0.5*psd->sd_shape_table[psd->index_type_mu][index_x];
        }
        else{
          g = 0.;
          psd->sd_table[index_type][index_x] = 0.;
        }
      }
      else{
        psd->sd_table[index_type][index_x] = psd->sd_parameter_table[index_type]*psd->sd_shape_table[index_type][index_x];
      }

      psd->DI[index_x] += psd->sd_table[index_type][index_x];
    }
  }

  /** Include additional sources of distortions */
  /* Superposition of blackbodies */
  //psd->sd_parameter_table[psd->index_type_y] += 2.525e-7;   // CMB Dipole (Chluba & Sunyaev 2004)
  //psd->sd_parameter_table[psd->index_type_y] += 4.59e-13;   // CMB Quadrupole (Chluba & Sunyaev 2004)

  /* Reionization */
  if(psd->has_SZ_effect == _TRUE_){
    for(index_x=0;index_x<psd->x_size;++index_x){
      class_call(distortions_add_effects_reio(pba,pth,psd,
                                              5.e0,
                                              2.e-4,
                                              1./300.,
                                              1./300.,
                                              psd->x[index_x],
                                              &y_reio,
                                              &DI_reio),
                 psd->error_message,
                 psd->error_message);
      psd->DI[index_x] += DI_reio;
    }

    psd->sd_parameter_table[psd->index_type_y] += y_reio;
  }

  /** Compute total heating */
  psd->Drho_over_rho = psd->sd_parameter_table[psd->index_type_y]*4.+
                       psd->sd_parameter_table[psd->index_type_mu]/1.401+
                       psd->epsilon;

  if(psd->include_g_distortion == _TRUE_){
     psd->Drho_over_rho += psd->sd_parameter_table[psd->index_type_g]*4.;
  }

  /** Print found parameters */
  if (psd->distortions_verbose > 1){

    if ( psd->distortions_verbose > 3 && psd->include_g_distortion){
      printf(" -> g-parameter %g (Note, that this does not include contributions from earlier than sd_z_max=%g)\n", psd->sd_parameter_table[psd->index_type_g], ppr->sd_z_max);
    }

    if (psd->sd_parameter_table[psd->index_type_mu] > 9.e-5) {
      printf(" -> mu-parameter = %g. WARNING: The value of your mu-parameter is larger than the FIRAS constraint mu<9e-5.\n", psd->sd_parameter_table[psd->index_type_mu]);
    }
    else{
      printf(" -> mu-parameter = %g\n", psd->sd_parameter_table[psd->index_type_mu]);
    }

    if (psd->sd_parameter_table[psd->index_type_y]>1.5e-5) {
      printf(" -> y-parameter = %g. WARNING: The value of your y-parameter is larger than the FIRAS constraint y<1.5e-5.\n", psd->sd_parameter_table[psd->index_type_y]);
    }
    else{
      printf(" -> y-parameter = %g\n", psd->sd_parameter_table[psd->index_type_y]);
    }

    if(psd->sd_branching_approx == bra_exact && psd->sd_PCA_size != 0){
       if(psd->distortions_verbose > 2){
         for(index_k=0; index_k<psd->sd_PCA_size; ++index_k){
           printf(" -> PCA multipole mu_%d = %g\n", index_k+1, psd->sd_parameter_table[psd->index_type_PCA+index_k]);
         }
       }
       printf(" -> epsilon-parameter = %g\n", psd->epsilon);
    }

    printf(" -> total injected/extracted heat = %g\n", psd->Drho_over_rho);
  }

  return _SUCCESS_;
}

/**
 * Compute relativistic contribution from reionization and structure formation according to
 *        1) Nozawa et al. 2005 (up to order 3 in theta_e) or
 *        2) Chluba et al. 2012 (up to order ? in ?). Note that, for the moment, this appoximation
 *           is only valid for cluster temperatures lower than few KeV.
 *
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to thermodynamics structure
 * @param psd        Input: pointer to the distortions structure
 * @param T_e        Input: electron temperature in keV
 * @param Dtau       Input: optical depth
 * @param beta       Input: peculiar velocity of the cluster
 * @param beta_z     Input: peculiar velocity of the cluster with respect to the line-of-sight
 * @param x          Input: dimensionless frequency
 * @param y_reio     Output: y-parameter
 * @param DI_reio    Output: spectral distortion
 * @return the error status
 */

int distortions_add_effects_reio(struct background * pba,
                                 struct thermodynamics * pth,
                                 struct distortions * psd,
                                 double T_e,
                                 double Dtau,
                                 double beta,
                                 double beta_z,
                                 double x,
                                 double * y_reio,
                                 double * DI_reio){

  /** Define local variables */
  double theta_e, cos_theta, P_1, P_2, x_tilde, S_tilde;
  double G_T, Y_SZ;
  double M_low, M_k, D_low, D_k, Q_low, Q_k;
  int index_k, index_n;
  double Y_0, Y_1, Y_2;
  double B_0, B_1, B_2, B_3;
  double C_0, C_1, C_2, C_3;
  double D_0, D_1, D_2, D_3;
  double DI_tSZ_non_rel, DI_tSZ_rel, DI_tSZ, DI_kSZ;

  /* Compute related quantities */
  theta_e = T_e*1.e3/(_m_e_/_GeV_over_kg_*1.e9);
  cos_theta = beta_z/beta;
  P_1 = cos_theta;
  P_2 = (3.*pow(cos_theta,2.)-1.)/2.;
  x_tilde = x/tanh(x/2.);  // coth=1/tanh
  S_tilde = x/sinh(x/2.);

  G_T = pow(x,4.)*exp(-x)/pow(1.-exp(-x),2.);
  Y_SZ = G_T*(x_tilde-4.);

  double Delta_x[8] = {-1.,
                       x_tilde,
                       -(pow(x_tilde,2.)
                                     +pow(S_tilde,2.)/2.),
                       x_tilde*(pow(x_tilde,2.)
                                     +pow(S_tilde,2.)*2.),
                       -(pow(x_tilde,4.)
                                     +pow(x_tilde,2.)*pow(S_tilde,2.)*11./2.
                                     +pow(S_tilde,4.)),
                       x_tilde*(pow(x_tilde,4.)
                                     +pow(x_tilde,2.)*pow(S_tilde,2.)*13.
                                     +pow(S_tilde,4.)*17./2.),
                       -(pow(x_tilde,6.)
                                     +pow(x_tilde,4.)*pow(S_tilde,2.)*57./2.
                                     +pow(x_tilde,2.)*pow(S_tilde,4.)*45.
                                     +pow(S_tilde,6.)*17./4.),
                       x_tilde*(pow(x_tilde,6.)
                                     +pow(x_tilde,4.)*pow(S_tilde,2.)*60.
                                     +pow(x_tilde,2.)*pow(S_tilde,4.)*192.
                                     +pow(S_tilde,6.)*62.)};

  /** Thermal SZ effect (TSZ) */
 /* Fill coefficient tables from appendix A1 of Chluba et al. 2012 */
  double a[6][3] = {{4., 10.,    15./2.},
                    {1., 47./2., 1023./8.},
                    {0., 42./5., 868./5.},
                    {0., 7./10., 329./5.},
                    {0., 0.,     44./5.},
                    {0., 0.,     11./30.}};

  double Y_k[3] = {0., 0., 0.};

  for(index_k=0; index_k<3; ++index_k){
    Y_k[index_k] = 0.;
    for(index_n=0; index_n<2*index_k+2; ++index_n){
      Y_k[index_k]+=a[index_n][index_k]*Delta_x[index_n];
    }
  }

  /** Non-relativistic TSZ */
  DI_tSZ_non_rel = Dtau*theta_e*G_T*Y_k[0];

  /** Relativistic TSZ */
  DI_tSZ_rel = 0.;
  for(index_k=1; index_k<3; ++index_k){
    DI_tSZ_rel += Dtau*pow(theta_e,index_k+1)*G_T*Y_k[index_k];
  }

  DI_tSZ = DI_tSZ_non_rel+DI_tSZ_rel;

  /** Kinematic SZ effect (kSZ) */
  /* Calculated according to Nozawa et al. 2005 */
  switch(psd->sd_reio_type){
    case sd_reio_Nozawa:
    Y_0 = Y_k[0];
    Y_1 = Y_k[1];
    Y_2 = Y_k[2];

    B_0 = 1.*Y_0/3.;
    B_1 = 5.*Y_0/6.
          +2.*Y_1/3.;
    B_2 = 5.*Y_0/8.
          +3.*Y_1/2.
          +Y_2;
    B_3 = -5.*Y_0/8.
          +5.*Y_1/4.
          +5.*Y_2/2.;

    C_0 = 1.;
    C_1 = 10.
          -47.*x_tilde/5.
          +7.*pow(x_tilde,2.)/5.
          +7.*pow(S_tilde,2.)/10.;
    C_2 = 25.
          -1117.*x_tilde/10.
          +847.*pow(x_tilde,2.)/10.
          -183.*pow(x_tilde,3.)/10.
          +11.*pow(x_tilde,4.)/10.
          +pow(S_tilde,2.)*(847./20.
                            -183.*x_tilde/5.
                            +121.*pow(x_tilde,2.)/20.)
          +11.*pow(S_tilde,4.)/10.;
    C_3 = 75./4.
          -21873.*x_tilde/40.
          +49161.*pow(x_tilde,2.)/40.
          -27519.*pow(x_tilde,3.)/35.
          +6684.*pow(x_tilde,4.)/35.
          -3917.*pow(x_tilde,5.)/210.
          +64.*pow(x_tilde,6.)/105.
          +pow(S_tilde,2.)*(49161./80.
                            -55038.*x_tilde/35.
                            +36762.*pow(x_tilde,2.)/35.
                            -50921.*pow(x_tilde,3.)/210.
                            +608.*pow(x_tilde,4.)/35.)
          +pow(S_tilde,4.)*(6684./35.
                            -66589.*x_tilde/420.
                            +192.*pow(x_tilde,2.)/7.)
          +272.*pow(S_tilde,6.)/105.;

    D_0 = -2./3.
          +11.*x_tilde/30.;
    D_1 = -4.
          +12.*x_tilde
          -6.*pow(x_tilde,2.)
          +19.*pow(x_tilde,3.)/30.
          +pow(S_tilde,2.)*(-3.
                            +19.*x_tilde/15.);
   D_2 = -10.
          +542.*x_tilde/5.
          -843.*pow(x_tilde,2.)/5.
          +10603.*pow(x_tilde,3.)/140.
          -409.*pow(x_tilde,4.)/35.
          +23.*pow(x_tilde,5.)/42.
          +pow(S_tilde,2.)*(-843./10.
                            +10603.*x_tilde/70.
                            -4499.*pow(x_tilde,2.)/70.
                            +299.*pow(x_tilde,3.)/42.)
          +pow(S_tilde,4.)*(-409./35.
                            +391.*x_tilde/84.);
    D_3 = -15./2.
          +4929.*x_tilde/40.
          -39777.*pow(x_tilde,2.)/20.
          +1199897.*pow(x_tilde,3.)/560.
          -4392.*pow(x_tilde,4.)/5.
          +16364.*pow(x_tilde,5.)/105.
          -3764.*pow(x_tilde,6.)/315.
          +101.*pow(x_tilde,7.)/315.
          +pow(S_tilde,2.)*(-39777./40.
                            +1199897.*x_tilde/280.
                            -24156.*pow(x_tilde,2.)/5.
                            +212732.*pow(x_tilde,3.)/105.
                            -35758.*pow(x_tilde,4.)/105.
                            +404.*pow(x_tilde,5.)/21.)
          +pow(S_tilde,4.)*(-4392./5.
                            +139094.*x_tilde/105.
                            -3764.*pow(x_tilde,2.)/7.
                            +6464.*pow(x_tilde,3.)/105.)
          +pow(S_tilde,6.)*(-15997./315.
                            +6262.*x_tilde/315.);

    M_low = G_T*(B_0+theta_e*B_1+pow(theta_e,2.)*B_2+pow(theta_e,3.)*B_3);
    D_low = G_T*(C_0+theta_e*C_1+pow(theta_e,2.)*C_2+pow(theta_e,3.)*C_3);
    Q_low = G_T*(D_0+theta_e*D_1+pow(theta_e,2.)*D_2+pow(theta_e,3.)*D_3);
    break;
  /* Calculated according to Chluba et al. 2012 */
    case sd_reio_Chluba:
    /* Low temperature approximation */
    if(T_e < 10.){
      double d[7][3] = {{-2./5., -1./5.,   407./140.},
                        {-8./5., -24./5., -233./35.},
                        {-2./5., -66./5., -10433./140.},
                        { 0.,    -24./5., -3876./35.},
                        { 0.,    -2./5.,  -1513./35.},
                        { 0.,     0.,     -204./35.},
                        { 0.,     0.,     -17./70.}};

      double q[7][4] = {{-3./5.,   183./70., -429./40.},
                        { 2./5.,  -5./7.,     207./20.},
                        { 1./10.,  115./28.,  1647./80.},
                        { 0.,      12./7.,    44.},
                        { 0.,      1./7.,     19.},
                        { 0.,      0.,        92./35.},
                        { 0.,      0.,        23./210.}};

      M_low = 1./3.*(Y_SZ+G_T);
      for(index_k=0; index_k<3; ++index_k){
        M_k = 0.;
        for(index_n=0; index_n<2*index_k+2; ++index_n){
          M_k += (a[index_n][index_k]-d[index_n][index_k])
                        *(index_n*(index_n+2)*Delta_x[index_n]+
                          (2*index_n+3)*Delta_x[index_n+1]+
                          Delta_x[index_n+2]
                         );
        }
        M_k *= 1./3.;
        M_low += pow(theta_e,index_k+1)*M_k*G_T;
      }

      D_low = G_T;
      for(index_k=0; index_k<3; ++index_k){
        D_k = 0.;
        for(index_n=0; index_n<2*index_k+2; ++index_n){
          D_k += (d[index_n][index_k]-a[index_n][index_k])
                        *(index_n*Delta_x[index_n]+
                          Delta_x[index_n+1]
                         );
        }
        D_low += pow(theta_e,index_k+1)*D_k*G_T;
      }

      Q_low = 11./30.*(x_tilde*G_T);
      for(index_k=0; index_k<3; ++index_k){
        Q_k = 0.;
        for(index_n=0; index_n<2*index_k+2; ++index_n){
          Q_k += (a[index_n][index_k]+q[index_n][index_k]-2.*d[index_n][index_k])
                        *(index_n*(index_n-1)*Delta_x[index_n]+
                          2*index_n*Delta_x[index_n+1]+
                          Delta_x[index_n+2]
                         );
        }
        Q_k *= 1./3.;
        Q_low += pow(theta_e,index_k+1)*Q_k*G_T;
      }
    }
    /* High temperature approximation (not implemented yet) */
    else{
      M_low = 0.;
      D_low = 0.;
      Q_low = 0.;
    }
    break;
    default:
      class_stop(psd->error_message,"Unrecognized sd_reio_type='%i'.",psd->sd_reio_type);
  }

  DI_kSZ = Dtau*beta*(beta*M_low+P_1*D_low+beta*P_2*Q_low);

  /** Total distortion */
  *y_reio = theta_e*Dtau;
  *DI_reio = DI_tSZ+DI_kSZ;

  return _SUCCESS_;
}

/**
 * Reads the external file branching_ratios calculated according to Chluba & Jeong 2014
 *
 * @param ppr        Input: pointer to precision structure
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */

int distortions_read_br_data(struct precision * ppr,
                             struct distortions * psd){

  /** Define local variables */
  int index_k,index_z;
  FILE * infile;
  DetectorFileName br_file;
  char line[_LINE_LENGTH_MAX_];
  char * left;
  int headlines = 0;

  /** Open file */
  sprintf(br_file,"%s/%s_branching_ratios.dat", ppr->sd_external_path, psd->sd_detector_name);
  class_open(infile, br_file, "r", psd->error_message);

  /** Read header */
  psd->br_exact_Nz = 0;
  while (fgets(line,_LINE_LENGTH_MAX_-1,infile) != NULL) {
    headlines++;

    /* Eliminate blank spaces at beginning of line */
    left=line;
    while (left[0]==' ') {
        left++;
    }

    if (left[0] > 39) {
      /** Read number of lines, infer size of arrays and allocate them */
      class_test(sscanf(line,"%d %d", &psd->br_exact_Nz, &psd->E_vec_size) != 2,
                 psd->error_message,
                 "could not header (number of lines, number of multipoles) at line %i in file '%s' \n",headlines,br_file);

      class_alloc(psd->br_exact_z, psd->br_exact_Nz*sizeof(double), psd->error_message);
      class_alloc(psd->f_g_exact, psd->br_exact_Nz*sizeof(double), psd->error_message);
      class_alloc(psd->f_y_exact, psd->br_exact_Nz*sizeof(double), psd->error_message);
      class_alloc(psd->f_mu_exact, psd->br_exact_Nz*sizeof(double), psd->error_message);

      class_alloc(psd->E_vec, psd->br_exact_Nz*psd->E_vec_size*sizeof(double), psd->error_message);
      break;
    }
  }

  /** Read parameters */
  for(index_z=0; index_z<psd->br_exact_Nz; ++index_z){
    class_test(fscanf(infile, "%le",
                      &(psd->br_exact_z[index_z]))!=1,                                // [-]
                      psd->error_message,
                      "Could not read z at line %i in file '%s'",index_z+headlines,br_file);
    class_test(fscanf(infile, "%le",
                      &(psd->f_g_exact[index_z]))!=1,                                 // [-]
                      psd->error_message,
                      "Could not read f_g at line %i in file '%s'",index_z+headlines,br_file);
    class_test(fscanf(infile, "%le",
                      &(psd->f_y_exact[index_z]))!=1,                                 // [-]
                      psd->error_message,
                      "Could not read f_y at line %i in file '%s'",index_z+headlines,br_file);
    class_test(fscanf(infile,"%le",
                      &(psd->f_mu_exact[index_z]))!=1,                                // [-]
                      psd->error_message,
                      "Could not read f_mu at line %i in file '%s'",index_z+headlines,br_file);
    for(index_k=0; index_k<psd->E_vec_size; ++index_k){
      class_test(fscanf(infile,"%le",
                        &(psd->E_vec[index_k*psd->br_exact_Nz+index_z]))!=1,           // [-]
                        psd->error_message,
                        "Could not read E vector at line %i in file '%s'",index_z+headlines,br_file);
    }
  }

  fclose(infile);

  return _SUCCESS_;
}

/**
 * Spline the quantitites read in distortions_read_br_data
 *
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */

int distortions_spline_br_data(struct distortions* psd){

  /** Allocate second derivatives */
  class_alloc(psd->ddf_g_exact,
              psd->br_exact_Nz*sizeof(double),
              psd->error_message);
  class_alloc(psd->ddf_y_exact,
              psd->br_exact_Nz*sizeof(double),
              psd->error_message);
  class_alloc(psd->ddf_mu_exact,
              psd->br_exact_Nz*sizeof(double),
              psd->error_message);
  class_alloc(psd->ddE_vec,
              psd->E_vec_size*psd->br_exact_Nz*sizeof(double),
              psd->error_message);

  /** Spline branching ratios */
  class_call(array_spline_table_columns(psd->br_exact_z,
                                        psd->br_exact_Nz,
                                        psd->f_g_exact,
                                        1,
                                        psd->ddf_g_exact,
                                        _SPLINE_EST_DERIV_,
                                        psd->error_message),
             psd->error_message,
             psd->error_message);
  class_call(array_spline_table_columns(psd->br_exact_z,
                                        psd->br_exact_Nz,
                                        psd->f_y_exact,
                                        1,
                                        psd->ddf_y_exact,
                                        _SPLINE_EST_DERIV_,
                                        psd->error_message),
           psd->error_message,
           psd->error_message);
  class_call(array_spline_table_columns(psd->br_exact_z,
                                        psd->br_exact_Nz,
                                        psd->f_mu_exact,
                                        1,
                                        psd->ddf_mu_exact,
                                        _SPLINE_EST_DERIV_,
                                        psd->error_message),
           psd->error_message,
           psd->error_message);
  class_call(array_spline_table_columns(psd->br_exact_z,
                                        psd->br_exact_Nz,
                                        psd->E_vec,
                                        psd->E_vec_size,
                                        psd->ddE_vec,
                                        _SPLINE_EST_DERIV_,
                                        psd->error_message),
           psd->error_message,
           psd->error_message);

  return _SUCCESS_;
}

/**
 * Interpolate the quantitites splined in distortions_spline_br_data
 *
 * @param psd        Input: pointer to the distortions structure
 * @param z          Input: redshift
 * @param f_g        Output: branching ratio for temperature shift
 * @param f_y        Output: branching ratio for y distortions
 * @param f_mu       Output: branching ratio for mu-distortions
 * @param f_E        Output: branching ratio for residuals (multipole expansion)
 * @param last_index Output: multipole of PCA expansion for f_E
 * @return the error status
 */

int distortions_interpolate_br_data(struct distortions* psd,
                                    double z,
                                    double * f_g,
                                    double * f_y,
                                    double * f_mu,
                                    double * f_E,
                                    int * last_index){

  /** Define local variables */
  int index = *last_index;
  int index_k;
  double h,a,b;

  /** Find z position */
  class_call(array_spline_hunt(psd->br_exact_z,
                               psd->br_exact_Nz,
                               z,
                               &index,
                               &h,&a,&b,
                               psd->error_message),
             psd->error_message,
             psd->error_message);

  /** Evaluate corresponding values for the branching ratios */
  *f_g = 4*array_spline_eval(psd->f_g_exact,
                             psd->ddf_g_exact,
                             index,
                             index+1,
                             h,a,b);
  *f_y = 4*array_spline_eval(psd->f_y_exact,
                             psd->ddf_y_exact,
                             index,
                             index+1,
                             h,a,b);
  *f_mu = 1./1.401*array_spline_eval(psd->f_mu_exact,
                                     psd->ddf_mu_exact,
                                     index,
                                     index+1,
                                     h,a,b);

  for(index_k=0; index_k<psd->sd_PCA_size; ++index_k){
    f_E[index_k] = array_spline_eval(psd->E_vec+index_k*psd->br_exact_Nz,
                                     psd->ddE_vec+index_k*psd->br_exact_Nz,
                                     index,
                                     index+1,
                                     h,a,b);
  }

  *last_index = index;

  return _SUCCESS_;
}

/**
 * Free from distortions_read_br_data and distortions_spline_br_data
 *
 * @param psd     Input: pointer to distortions structure (to be freed)
 * @return the error status
 */

int distortions_free_br_data(struct distortions * psd){

  free(psd->br_exact_z);
  free(psd->f_g_exact);
  free(psd->ddf_g_exact);
  free(psd->f_y_exact);
  free(psd->ddf_y_exact);
  free(psd->f_mu_exact);
  free(psd->ddf_mu_exact);
  free(psd->E_vec);
  free(psd->ddE_vec);

  return _SUCCESS_;
}

/**
 * Reads the external file distortions_shapes calculated according to Chluba & Jeong 2014
 *
 * @param ppr        Input: pointer to precision structure
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */

int distortions_read_sd_data(struct precision * ppr,
                             struct distortions * psd){

  /** Define local variables */
  FILE * infile;
  DetectorFileName sd_file;
  char line[_LINE_LENGTH_MAX_];
  char * left;
  int headlines = 0;
  int index_x,index_k;

  /** Open file */
  sprintf(sd_file,"%s/%s_distortions_shapes.dat",ppr->sd_external_path, psd->sd_detector_name);
  class_open(infile, sd_file, "r", psd->error_message);

  /** Read header */
  psd->PCA_Nnu = 0;
  while (fgets(line,_LINE_LENGTH_MAX_-1,infile) != NULL) {
    headlines++;

    /* Eliminate blank spaces at beginning of line */
    left=line;
    while (left[0]==' ') {
        left++;
    }

    if (left[0] > 39) {
      /** Read number of lines, infer size of arrays and allocate them */
      class_test(sscanf(line, "%d %d", &psd->PCA_Nnu, &psd->S_vec_size) != 2,
                 psd->error_message,
                 "could not header (number of lines, number of multipoles) at line %i in file '%s' \n",headlines,sd_file);

      class_alloc(psd->PCA_nu, psd->PCA_Nnu*sizeof(double), psd->error_message);
      class_alloc(psd->PCA_G_T, psd->PCA_Nnu*sizeof(double), psd->error_message);
      class_alloc(psd->PCA_Y_SZ, psd->PCA_Nnu*sizeof(double), psd->error_message);
      class_alloc(psd->PCA_M_mu, psd->PCA_Nnu*sizeof(double), psd->error_message);

      class_alloc(psd->S_vec, psd->PCA_Nnu*psd->S_vec_size*sizeof(double), psd->error_message);
      break;
    }
  }

  /** Read parameters */
  for(index_x=0; index_x<psd->PCA_Nnu; ++index_x){
    class_test(fscanf(infile,"%le",
                      &(psd->PCA_nu[index_x]))!=1,                                                  // [GHz]
                      psd->error_message,
                      "Could not read z at line %i in file '%s'",index_x+headlines,sd_file);
    class_test(fscanf(infile,"%le",
                      &(psd->PCA_G_T[index_x]))!=1,                                                 // [10^-18 W/(m^2 Hz sr)]
                      psd->error_message,
                      "Could not read f_g at line %i in file '%s'",index_x+headlines,sd_file);
    psd->PCA_G_T[index_x] /= (psd->DI_units*1.e18);                                                 // [-]
    class_test(fscanf(infile,"%le",
                      &(psd->PCA_Y_SZ[index_x]))!=1,                                                // [10^-18 W/(m^2 Hz sr)]
                      psd->error_message,
                      "Could not read f_y at line %i in file '%s'",index_x+headlines,sd_file);
    psd->PCA_Y_SZ[index_x] /= (psd->DI_units*1.e18);                                                // [-]
    class_test(fscanf(infile,"%le",
                      &(psd->PCA_M_mu[index_x]))!=1,                                                // [10^-18 W/(m^2 Hz sr)]
                      psd->error_message,
                      "Could not read f_mu at line %i in file '%s'",index_x+headlines,sd_file);
    psd->PCA_M_mu[index_x] /= (psd->DI_units*1.e18);                                                // [-]
    for(index_k=0; index_k<psd->S_vec_size; ++index_k){
      class_test(fscanf(infile,"%le",
                        &(psd->S_vec[index_k*psd->PCA_Nnu+index_x]))!=1,                            // [10^-18 W/(m^2 Hz sr)]
                        psd->error_message,
                        "Could not read E vector at line %i in file '%s'",index_x+headlines,sd_file);
      psd->S_vec[index_k*psd->PCA_Nnu+index_x] /= (psd->DI_units*1.e18);                            // [-]
    }

  }

  fclose(infile);

  return _SUCCESS_;
}

/**
 * Spline the quantitites read in distortions_read_sd_data
 *
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */

int distortions_spline_sd_data(struct distortions* psd){

  /** Allocate second derivatievs */
  class_alloc(psd->ddPCA_G_T,
              psd->PCA_Nnu*sizeof(double),
              psd->error_message);
  class_alloc(psd->ddPCA_Y_SZ,
              psd->PCA_Nnu*sizeof(double),
              psd->error_message);
  class_alloc(psd->ddPCA_M_mu,
              psd->PCA_Nnu*sizeof(double),
              psd->error_message);
  class_alloc(psd->ddS_vec,
              psd->S_vec_size*psd->PCA_Nnu*sizeof(double),
              psd->error_message);

  /** Spline branching ratios */
  class_call(array_spline_table_columns(psd->PCA_nu,
                                        psd->PCA_Nnu,
                                        psd->PCA_G_T,
                                        1,
                                        psd->ddPCA_G_T,
                                        _SPLINE_EST_DERIV_,
                                        psd->error_message),
             psd->error_message,
             psd->error_message);
  class_call(array_spline_table_columns(psd->PCA_nu,
                                        psd->PCA_Nnu,
                                        psd->PCA_Y_SZ,
                                        1,
                                        psd->ddPCA_Y_SZ,
                                        _SPLINE_EST_DERIV_,
                                        psd->error_message),
           psd->error_message,
           psd->error_message);
  class_call(array_spline_table_columns(psd->PCA_nu,
                                        psd->PCA_Nnu,
                                        psd->PCA_M_mu,
                                        1,
                                        psd->ddPCA_M_mu,
                                        _SPLINE_EST_DERIV_,
                                        psd->error_message),
           psd->error_message,
           psd->error_message);
  class_call(array_spline_table_columns(psd->PCA_nu,
                                        psd->PCA_Nnu,
                                        psd->S_vec,
                                        psd->S_vec_size,
                                        psd->ddS_vec,
                                        _SPLINE_EST_DERIV_,
                                        psd->error_message),
           psd->error_message,
           psd->error_message);

  return _SUCCESS_;
}

/**
 * Interpolate the quantitites splined in distortions_spline_sd_data
 *
 * @param psd        Input: pointer to the distortions structure
 * @param nu         Input: dimnetionless frequency
 * @param G_T        Output: shape of temperature shift
 * @param Y_SZ       Output: shape of y distortions
 * @param M_mu       Output: shape of mu-distortions
 * @param S          Output: shape of residuals (multipole expansion)
 * @param index      Output: multipole of PCA expansion for S
 * @return the error status
 */

int distortions_interpolate_sd_data(struct distortions* psd,
                                    double nu,
                                    double * G_T,
                                    double * Y_SZ,
                                    double * M_mu,
                                    double * S,
                                    int * index){

  /** Define local variables */
  int last_index = *index;
  int index_k;
  double h,a,b;
  double nu_round;

  /** Find z position */
  nu_round = round(nu*pow(10.,3))/pow(10.,3); // The rounding is necessary for the interpolation with the external file
  class_call(array_spline_hunt(psd->PCA_nu,
                               psd->PCA_Nnu,
                               nu_round,
                               &last_index,
                               &h,&a,&b,
                               psd->error_message),
             psd->error_message,
             psd->error_message);

  /** Evaluate corresponding values for the branching ratios */
  *G_T = array_spline_eval(psd->PCA_G_T,
                           psd->ddPCA_G_T,
                           last_index,
                           last_index+1,
                           h,a,b);
  *Y_SZ = array_spline_eval(psd->PCA_Y_SZ,
                            psd->ddPCA_Y_SZ,
                            last_index,
                            last_index+1,
                            h,a,b);
  *M_mu = array_spline_eval(psd->PCA_M_mu,
                            psd->ddPCA_M_mu,
                            last_index,
                            last_index+1,
                            h,a,b);

  for(index_k=0; index_k<psd->sd_PCA_size; ++index_k){
    S[index_k] = array_spline_eval(psd->S_vec+index_k*psd->PCA_Nnu,
                                   psd->ddS_vec+index_k*psd->PCA_Nnu,
                                   last_index,
                                   last_index+1,
                                   h,a,b);
  }

  *index = last_index;

  return _SUCCESS_;
}

/**
 * Free from distortions_read_sd_data and distortions_spline_sd_data
 *
 * @param psd     Input: pointer to distortions structure (in which some fields should be freed)
 * @return the error status
 */

int distortions_free_sd_data(struct distortions * psd){

  free(psd->PCA_nu);
  free(psd->PCA_G_T);
  free(psd->ddPCA_G_T);
  free(psd->PCA_Y_SZ);
  free(psd->ddPCA_Y_SZ);
  free(psd->PCA_M_mu);
  free(psd->ddPCA_M_mu);
  free(psd->S_vec);
  free(psd->ddS_vec);

  return _SUCCESS_;
}

/**
 * Define title of columns in the heat output
 *
 * @param psd     Input: pointer to distortions structure
 * @param titles  Output: title of each column in the output
 */

int distortions_output_heat_titles(struct distortions * psd,
                                   char titles[_MAXTITLESTRINGLENGTH_]){

  class_store_columntitle(titles,"Redshift z",_TRUE_);
  class_store_columntitle(titles,"Heat  [-]",_TRUE_);
  class_store_columntitle(titles,"LHeat [-]",_TRUE_);

  return _SUCCESS_;
}

/**
 * Store data in the heat output
 *
 * @param psd              Input/Output: pointer to distortions structure
 * @param number_of_titles Input: numbert of column in the output
 * @param data             Input: data to be stored
 */

int distortions_output_heat_data(struct distortions * psd,
                                 int number_of_titles,
                                 double * data){
  int storeidx;
  double * dataptr;
  int index_z;

  for (index_z=0; index_z<psd->z_size; index_z++) {
    dataptr = data + index_z*number_of_titles;
    storeidx = 0;
    class_store_double(dataptr, psd->z[index_z], _TRUE_, storeidx);
    class_store_double(dataptr, psd->dQrho_dz_tot[index_z], _TRUE_, storeidx);
    class_store_double(dataptr, psd->dQrho_dz_tot[index_z]*(1.+psd->z[index_z]), _TRUE_, storeidx);
  }

  return _SUCCESS_;
}

/**
 * Define title of columns in the spectral distortion output
 *
 * @param psd     Input: pointer to distortions structure
 * @param titles  Output: title of each column in the output
 */

int distortions_output_sd_titles(struct distortions * psd,
                                 char titles[_MAXTITLESTRINGLENGTH_]){

  char temp_title[256];
  int index_type;
  class_store_columntitle(titles,"Dimensionless frequency x",_TRUE_);
  class_store_columntitle(titles,"Frequency nu [GHz]",_TRUE_);
  class_store_columntitle(titles,"SD_tot",_TRUE_);
  for(index_type=0;index_type<psd->type_size;++index_type){
    if(index_type==psd->index_type_g){
      sprintf(temp_title,"SD[g]");
    }
    if(index_type==psd->index_type_y){
      sprintf(temp_title,"SD[y]");
    }
    if(index_type==psd->index_type_mu){
      sprintf(temp_title,"SD[mu]");
    }
    if(index_type>=psd->index_type_PCA){
      sprintf(temp_title,"SD[e_%i]",(index_type-psd->index_type_PCA));
    }
    class_store_columntitle(titles,temp_title,_TRUE_);
  }

  return _SUCCESS_;
}

/**
 * Store data in the distortion output
 *
 * @param psd              Input/Output: pointer to distortions structure
 * @param number_of_titles Input: numbert of column in the output
 * @param data             Input: data to be stored
 */

int distortions_output_sd_data(struct distortions * psd,
                               int number_of_titles,
                               double * data){
  int index_type;
  int storeidx;
  double * dataptr;
  int index_x;

  for (index_x=0; index_x<psd->x_size; index_x++) {
    dataptr = data + index_x*number_of_titles;
    storeidx = 0;

    class_store_double(dataptr, psd->x[index_x], _TRUE_,storeidx);
    class_store_double(dataptr, psd->x[index_x]*psd->x_to_nu, _TRUE_,storeidx);
    class_store_double(dataptr, psd->DI[index_x]*1.e26*psd->DI_units, _TRUE_,storeidx);
    for(index_type=0;index_type<psd->type_size;++index_type){
      class_store_double(dataptr, psd->sd_table[index_type][index_x]*1.e26*psd->DI_units, _TRUE_,storeidx);
    }
  }

  return _SUCCESS_;
}
