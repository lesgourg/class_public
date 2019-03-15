/** @file distortions.c Documented module on spectral distortions
 * Matteo Lucca, 31.10.2018
 * Nils Schoeneberg, 18.02.2019
 */

#include "distortions.h"

/**
 * Initialize the distortions structure.
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param ppt        Input: pointer to the perturbations structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param ppm        Input: pointer to the primordial structure
 * @param psd        Input/Output: pointer to initialized distortions structure
 * @return the error status
 */
int distortions_init(struct precision * ppr,
                     struct background * pba,
                     struct thermo * pth,
                     struct perturbs * ppt,
                     struct primordial * ppm,
                     struct distortions * psd) {

  /** Define local variables */
  int last_index = 0;
  int index_br;
  int index_x;

  struct heating* phe = &(pth->he);

  if(psd->has_distortions == _FALSE_){
    return _SUCCESS_;
  }
  if (psd->distortions_verbose > 0) {
    printf("Computing spectral distortions\n");
  }

  class_test(pth->compute_damping_scale==_FALSE_,
             psd->error_message,
             "Cannot compute spectral distortions without damping scale\n");

  /** Set physical constants */
  class_call(distortions_constants(pba,pth,psd),
             psd->error_message,
             psd->error_message);

  if(psd->branching_approx == bra_exact){
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
  class_call(distortions_compute_heating_rate(pba,pth,phe,ppt,ppm,psd),
             psd->error_message,
             psd->error_message);

  /** Define final spectral distortions */
  class_call(distortions_compute_spectral_shapes(ppr,pba,psd),
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
  int index_type,index_ht;

  if(psd->has_distortions == _TRUE_){
    /** Delete lists */
    free(psd->z);
    free(psd->x);
    free(psd->z_weights);

    /** Delete branching ratios */
    for(index_type=0;index_type<psd->type_size;++index_type){
      free(psd->br_table[index_type]);
    }
    free(psd->br_table);

    /** Delete heating functions */
    free(psd->dQrho_dz_tot);
    free(psd->dQrho_dz_tot_screened);

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
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to thermodynamics structure
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */
int distortions_constants(struct background * pba,
                          struct thermo * pth,
                          struct distortions * psd){

  /** Define unit conventions */
  psd->x_to_nu = (_k_B_*pba->T_cmb/_h_P_)/1e9;
  psd->DI_units = 2.*pow(_k_B_*pba->T_cmb,3.)/pow(_h_P_*_c_,2.);

  /** Define transition redshifts z_muy and z_th */
  psd->z_muy = 5.e4;
  psd->z_th = 1.98e6*
         pow((1.-pth->YHe/2.)/0.8767,-2./5.)*
         pow(pba->Omega0_b*pow(pba->h,2.)/0.02225,-2./5.)*
         pow(pba->T_cmb/2.726,1./5.);

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
  double nu_min,nu_max,nu_delta,delta_Ic;
  char * left;
  int headlines = 0;
  int i,j;
  int found_detector;

  if(psd->user_defined_name == _FALSE_){
    /* The user wants the default */
    if(psd->user_defined_detector == _FALSE_){
      return _SUCCESS_; //Nothing more to do
    }
    /* The user wants a new detector with specified settings, but without name */
    else{
      /* Generate a custom name for this custom detector, so we can check if it has already been defined */
      sprintf(psd->distortions_detector,"Custom__%7.2e_%7.2e_%7.2e_%7.2e__Detector",psd->nu_min_detector,psd->nu_max_detector,psd->nu_delta_detector,psd->delta_Ic_detector);
    }
  }

  /** Open file */
  class_open(det_list_file, "./external/distortions/detectors_list.dat", "r",
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
      class_test(sscanf(line,"%s %lg %lg %lg %lg",detector_name,&nu_min,&nu_max,&nu_delta,&delta_Ic) != 5,
                 psd->error_message,
                 "Could not read line %i in file '%s'\n",headlines,"./external/distortions/known_detectors.dat");

      /* Detector has been found */
      if(strcmp(psd->distortions_detector,detector_name)==0){
        printf(" -> Found detector %s (user defined = %s)\n",detector_name,(psd->user_defined_detector?"TRUE":"FALSE"));
        found_detector = _TRUE_;

        if (psd->distortions_verbose > 1){
          printf("    Properties:    nu_min = %lg    nu_max = %lg    delta_nu = %lg    delta_Ic = %lg \n",nu_min, nu_max, nu_delta, delta_Ic);
        }
        /* If the user has defined the detector, check that their and our definitions agree */
        if(psd->user_defined_detector){
          class_test(fabs(psd->nu_min_detector-nu_min)>ppr->tol_distortions_detector,
                     psd->error_message,
                     "Minimal frequency (nu_min) disagrees between stored detector '%s' and input ->  %.10e (input) vs %.10e (stored)",detector_name,psd->nu_min_detector,nu_min);
          class_test(fabs(psd->nu_max_detector-nu_max)>ppr->tol_distortions_detector,
                     psd->error_message,
                     "Maximal frequency (nu_max) disagrees between stored detector '%s' and input ->  %.10e (input) vs %.10e (stored)",detector_name,psd->nu_max_detector,nu_max);
          class_test(fabs(psd->nu_delta_detector-nu_delta)>ppr->tol_distortions_detector,
                     psd->error_message,
                     "Delta frequency (nu_delta) disagrees between stored detector '%s' and input ->  %.10e (input) vs %.10e (stored)",detector_name,psd->nu_max_detector,nu_max);
          class_test(fabs(psd->delta_Ic_detector-delta_Ic)>ppr->tol_distortions_detector,
                     psd->error_message,
                     "Detector accuracy (delta_Ic) disagrees between stored detector '%s' and input ->  %.10e (input) vs %.10e (stored)",detector_name,psd->delta_Ic_detector,delta_Ic);
        }

        /* In any case, just take the detector definition from the file */
        psd->nu_min_detector = nu_min;
        psd->nu_max_detector = nu_max;
        psd->nu_delta_detector = nu_delta;
        psd->delta_Ic_detector = delta_Ic;
      }
    }
  }

  /* If the detector has not been found, either the user has specified the settings and we create a new one,
   * or the user hasn't specified the settings and we have to stop */
  if(found_detector == _FALSE_){
    if(psd->user_defined_detector){
      printf(" -> Generating detector '%s' \n",psd->distortions_detector);
      class_call(distortions_generate_detector(ppr,psd),
                 psd->error_message,
                 psd->error_message);
    }
    else{
      class_stop(psd->error_message,
                 "You asked for detector '%s', but it was not in the database '%s'.\nPlease check the name of your detector, or specify its properties if you want to create a new one",
                 psd->distortions_detector,"./external/distortions/known_detectors.dat");
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
  char temporary_string[500];


  /* Test first whether or not python exists*/
  printf(" -> Testing python\n");
  is_success = system("python --version");
  class_test(is_success == -1,
             psd->error_message,
             "The command 'python --version' failed.\nPlease install a valid version of python.");

  /* Then activate the PCA generator*/
  printf(" -> Executing the PCA generator\n");
  sprintf(temporary_string,"python ./external/distortions/generate_PCA_files.py %s %.10e %.10e %.10e %.10e %.10e %i %.10e %i %.10e %.10e %.10e",
          psd->distortions_detector,
          psd->nu_min_detector,
          psd->nu_max_detector,
          psd->nu_delta_detector,
          ppr->distortions_z_min,
          ppr->distortions_z_max,
          ppr->distortions_z_size,
          psd->delta_Ic_detector,
          6,
          psd->z_th,
          psd->DI_units,
          psd->x_to_nu);
  is_success = system(temporary_string);
  class_test(is_success == -1,
             psd->error_message,
             "The command 'python ./external/distortions/generate_PCA_files.py' failed.\nPlease make sure the file exists.");

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
  class_define_index(psd->index_type_PCA,_TRUE_,index_type,psd->N_PCA);

  psd->type_size = index_type;

  return _SUCCESS_;
}


/**
 * Calculate redshift and frequency vectors and weights for redshift integral.
 *
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param psd        Input/Output: pointer to initialized distortions structure
 * @return the error status
 */
int distortions_get_xz_lists(struct precision * ppr,
                             struct background * pba,
                             struct thermo * pth,
                             struct distortions * psd){

  /** Define local variables */
  int index_z, index_x;

  /** Define and allocate z array */
  psd->z_min = ppr->distortions_z_min;
  psd->z_max = ppr->distortions_z_max;
  psd->z_size = ppr->distortions_z_size;
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
  if(psd->branching_approx != bra_exact){
    psd->x_min = ppr->distortions_x_min;
    psd->x_max = ppr->distortions_x_max;
    psd->x_size = ppr->distortions_x_size;
    psd->x_delta = (log(psd->x_max)-log(psd->x_min))/psd->x_size;

    class_alloc(psd->x,
                psd->x_size*sizeof(double),
                psd->error_message);

    for (index_x = 0; index_x<psd->x_size; index_x++) {
      psd->x[index_x] = exp(log(psd->x_min)+psd->x_delta*index_x);
    }

  }
  else{
    psd->x_min = psd->nu_min_detector/psd->x_to_nu;
    psd->x_max = psd->nu_max_detector/psd->x_to_nu;
    psd->x_delta = psd->nu_delta_detector/psd->x_to_nu;
    psd->x_size = (psd->x_max-psd->x_min)/psd->x_delta;

    class_alloc(psd->x,
                psd->x_size*sizeof(double),
                psd->error_message);

    for (index_x = 0; index_x<psd->x_size; index_x++) {
      psd->x[index_x] = psd->x_min+psd->x_delta*index_x;
    }
  }

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
  if(psd->branching_approx != bra_exact){
    for(index_z=0; index_z<psd->z_size; ++index_z){
      bb_vis = exp(-pow(psd->z[index_z]/psd->z_th,2.5));

      /* 1) Calculate branching ratios using sharp_sharp transition */
      if(psd->branching_approx == bra_sharp_sharp){
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
      if(psd->branching_approx == bra_sharp_soft){
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
      if(psd->branching_approx == bra_soft_soft){
        f_g = 1.-bb_vis;
        f_y = 1.0/(1.0+pow((1.0+psd->z[index_z])/(6.0e4),2.58));
        f_mu = bb_vis*(1.-exp(-pow((1.0+psd->z[index_z])/(5.8e4),1.88)));
      }

      /* 4) Calculate branching ratios unsing soft_soft_cons transitions */
      if(psd->branching_approx == bra_soft_soft_cons){
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
                psd->N_PCA*sizeof(double),
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
      for(index_k=0; index_k<psd->N_PCA; ++index_k){
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
 * @param ppt        Input: pointer to the perturbations structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param ppm        Input: pointer to the primordial structure
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */
int distortions_compute_heating_rate(struct background* pba,
                                     struct thermo * pth,
                                     struct heating * phe,
                                     struct perturbs * ppt,
                                     struct primordial * ppm,
                                     struct distortions * psd){

  /** Define local variables */
  int index_z;
  int last_index_back, last_index_thermo;
  double *pvecback, *pvecthermo;
  double tau;
  double x_e, T_b;
  double bb_vis;

  /** Update heating table with second order contributions */
  class_call(heating_add_second_order(pba,
                                      pth,
                                      ppt,
                                      ppm),
             phe->error_message,
             psd->error_message);

  /** Allocate space for background vector */
  last_index_back = 0;
  last_index_thermo = 0;
  class_alloc(pvecback,
              pba->bg_size*sizeof(double),
              psd->error_message);
  class_alloc(pvecthermo,
              pba->bg_size*sizeof(double),
              psd->error_message);

  /** Allocate space for total heating function */
  class_alloc(psd->dQrho_dz_tot,
              psd->z_size*sizeof(double*),
              psd->error_message);
  class_alloc(psd->dQrho_dz_tot_screened,
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
                                 pba->long_info,
                                 pba->inter_closeby,
                                 &last_index_back,
                                 pvecback),
               pba->error_message,
               psd->error_message);

    /** Import quantities from thermodynamics structure */
    class_call(thermodynamics_at_z(pba,
                                   pth,
                                   psd->z[index_z],
                                   pth->inter_normal,
                                   &last_index_thermo,
                                   pvecback,
                                   pvecthermo),
               pth->error_message,
               psd->error_message);

    x_e = pvecthermo[pth->index_th_xe];                                                             // [-]
    T_b = pvecthermo[pth->index_th_Tb];                                                             // [-]

    /* Black body visibility function */
    bb_vis = exp(-pow(psd->z[index_z]/psd->z_th,2.5));

    /* Total heating rate */
    /*
    class_call(heating_at_z(pba,pth,
                            x_e,
                            psd->z[index_z],
                            T_b,
                            pvecback),
               phe->error_message,
               psd->error_message);

     psd->dQrho_dz_tot[index_z] = phe->deposition_table[index_z*phe->dep_size+index_dep_heat]*
                                  pvecback[pba->index_bg_a]/
                                  (pvecback[pba->index_bg_H]*_c_/_Mpc_over_m_)/
                                  (pvecback[pba->index_bg_rho_g]*_GeVcm3_over_Mpc2_*_eV_*1e9*1e6); // [-]
     printf("%g  %g\n", psd->z[index_z],psd->dQrho_dz_tot[index_z]*(1.+psd->z[index_z]));*/
    psd->dQrho_dz_tot[index_z] = 0.;

    psd->dQrho_dz_tot_screened[index_z] = psd->dQrho_dz_tot[index_z]*bb_vis;                        // [-]
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
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */
int distortions_compute_spectral_shapes(struct precision * ppr,
                                        struct background * pba,
                                        struct distortions * psd){

  /** Define local variables */
  double * S;
  int last_index = 0;
  int index_type, index_z, index_x, index_k;
  double * integrand;
  double sum_S, sum_G;
  double g;

  /** Allocate space for spectral distortion amplitude in table sd_parameter_table */
  class_alloc(psd->sd_parameter_table,
              psd->type_size*sizeof(double),
              psd->error_message);

  /** Compute distortion amplitudes corresponding to each branching ratio (g, y and mu) */

  /* Define y, mu, g and mu_k from heating rates */
  for(index_type=0; index_type<psd->type_size; ++index_type){
    class_call(array_trapezoidal_convolution(psd->dQrho_dz_tot_screened,
                                             psd->br_table[index_type],
                                             psd->z_size,
                                             psd->z_weights,
                                             &(psd->sd_parameter_table[index_type]),
                                             psd->error_message),
               psd->error_message,
               psd->error_message);

    if(index_type>=psd->index_type_PCA){
      /* The E_k are not properly normalized, we have to renormalize here */
      psd->sd_parameter_table[index_type]/=(log(1.+psd->z[1])-log(1.+psd->z[0]));
    }
  }

  psd->sd_parameter_table[psd->index_type_g] /= 4.;
  psd->sd_parameter_table[psd->index_type_y] /= 4.;
  psd->sd_parameter_table[psd->index_type_mu] *= 1.401;

  /** Include additional sources of distortions (see also Chluba 2016 for useful discussion) */
  psd->sd_parameter_table[psd->index_type_y] += 2.525e-7;   // CMB Dipole (Chluba & Sunyaev 2004)
  psd->sd_parameter_table[psd->index_type_y] += 4.59e-13;   // CMB Quadrupole (Chluba & Sunyaev 2004)
  psd->sd_parameter_table[psd->index_type_y] += 1.77e-6;    // Reionization and structure formation (Hill et al. 2015)

  /* Print found parameters */
  if (psd->distortions_verbose > 1){

    printf(" -> g-parameter %g\n", psd->sd_parameter_table[psd->index_type_g]);
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
  }

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
  if(psd->branching_approx != bra_exact || psd->N_PCA == 0){
    /* If no PCA analysis is required, the shapes have simple analistical form */
    for(index_x=0; index_x<psd->x_size; ++index_x){
      psd->sd_shape_table[psd->index_type_g][index_x] = pow(psd->x[index_x],4.)*exp(-psd->x[index_x])/
                                                           pow(1.-exp(-psd->x[index_x]),2.);   // [-]
      psd->sd_shape_table[psd->index_type_y][index_x] = psd->sd_shape_table[psd->index_type_g][index_x]*
                                                           (psd->x[index_x]*(1.+exp(-psd->x[index_x]))/
                                                           (1.-exp(-psd->x[index_x]))-4.);     // [-]
      psd->sd_shape_table[psd->index_type_mu][index_x] = psd->sd_shape_table[psd->index_type_g][index_x]*
                                                           (1./2.19229-1./psd->x[index_x]);    // [-]
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
                psd->N_PCA*sizeof(double),
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

      for(index_k=0; index_k<psd->N_PCA; ++index_k){
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

  if(psd->branching_approx == bra_exact && psd->N_PCA != 0){
    for(index_k=0; index_k<psd->N_PCA; ++index_k){
      sum_S = 0.;
      sum_G = 0.;
      for(index_x=0; index_x<psd->x_size; ++index_x){
        sum_S += psd->sd_shape_table[psd->index_type_PCA+index_k][index_x];
        sum_G += psd->sd_shape_table[psd->index_type_g][index_x];
      }
      psd->epsilon += (4.*sum_S/sum_G)*psd->sd_parameter_table[psd->index_type_PCA+index_k];
    }

    /* Print found parameters */
    if (psd->distortions_verbose > 1){
       for(index_k=0; index_k<psd->N_PCA; ++index_k){
         printf(" -> PCA multipole mu_%d = %g\n", index_k+1, psd->sd_parameter_table[psd->index_type_PCA+index_k]);
       }

       printf(" -> epsilon-parameter = %g\n", psd->epsilon);
    }
  }

  /** Compute total heating */
  psd->Drho_over_rho = psd->sd_parameter_table[psd->index_type_g]*4.+
                       psd->sd_parameter_table[psd->index_type_y]*4.+
                       psd->sd_parameter_table[psd->index_type_mu]/1.401+
                       psd->epsilon;

  /* Print found parameter */
  if (psd->distortions_verbose > 1){
    printf(" -> total injected/extracted heat = %g\n", psd->Drho_over_rho);
  }

  /** Allocate space for final spactral distortion */
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
        g = psd->sd_parameter_table[psd->index_type_g];
        psd->sd_table[index_type][index_x] = (1.+g)*g*psd->sd_shape_table[psd->index_type_g][index_x]+
            g*g*0.5*psd->sd_shape_table[psd->index_type_mu][index_x];
      }
      else{
        psd->sd_table[index_type][index_x] = psd->sd_parameter_table[index_type]*psd->sd_shape_table[index_type][index_x];
      }

      psd->DI[index_x] += psd->sd_table[index_type][index_x];
    }
  }

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
  char br_file[500];
  char line[_LINE_LENGTH_MAX_];
  char * left;
  int headlines = 0;

  /** Open file */
  sprintf(br_file,"external/distortions/%s_branching_ratios.dat", psd->distortions_detector);
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
                 "could not header (number of lines, number of columns, number of multipoles) at line %i in file '%s' \n",headlines,br_file);

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

  return _SUCCESS_;

}


/**
 * Spline the quantitites read in distortions_read_br_data
 *
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */
int distortions_spline_br_data(struct distortions* psd){

  /** Allocate second derivatievs */
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
 * @param index      Output: multipole of PCA expansion for f_E
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
                               &h,
                               &a,
                               &b,
                               psd->error_message),
             psd->error_message,
             psd->error_message);

  /** Evaluate corresponding values for the branching ratios */
  *f_g = 4*array_interpolate_spline_hunt(psd->f_g_exact,
                                         psd->ddf_g_exact,
                                         index,
                                         index+1,
                                         h,a,b);
  *f_y = 4*array_interpolate_spline_hunt(psd->f_y_exact,
                                         psd->ddf_y_exact,
                                         index,
                                         index+1,
                                         h,a,b);
  *f_mu = 1./1.401*array_interpolate_spline_hunt(psd->f_mu_exact,
                                                 psd->ddf_mu_exact,
                                                 index,
                                                 index+1,
                                                 h,a,b);

  for(index_k=0; index_k<psd->N_PCA; ++index_k){
    f_E[index_k] = array_interpolate_spline_hunt(psd->E_vec+index_k*psd->br_exact_Nz,
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
  char sd_file[500];
  char line[_LINE_LENGTH_MAX_];
  char * left;
  int headlines = 0;
  int index_x,index_k;

  /** Open file */
  sprintf(sd_file,"external/distortions/%s_distortions_shapes.dat", psd->distortions_detector);
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
                 "could not header (number of lines, number of columns, number of multipoles) at line %i in file '%s' \n",headlines,sd_file);

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
                      &(psd->PCA_nu[index_x]))!=1,                                          // [GHz]
                      psd->error_message,
                      "Could not read z at line %i in file '%s'",index_x+headlines,sd_file);
    class_test(fscanf(infile,"%le",
                      &(psd->PCA_G_T[index_x]))!=1,                                         // [10^-18 W/(m^2 Hz sr)]
                      psd->error_message,
                      "Could not read f_g at line %i in file '%s'",index_x+headlines,sd_file);
    class_test(fscanf(infile,"%le",
                      &(psd->PCA_Y_SZ[index_x]))!=1,                                        // [10^-18 W/(m^2 Hz sr)]
                      psd->error_message,
                      "Could not read f_y at line %i in file '%s'",index_x+headlines,sd_file);
    class_test(fscanf(infile,"%le",
                      &(psd->PCA_M_mu[index_x]))!=1,                                        // [10^-18 W/(m^2 Hz sr)]
                      psd->error_message,
                      "Could not read f_mu at line %i in file '%s'",index_x+headlines,sd_file);
    for(index_k=0; index_k<psd->S_vec_size; ++index_k){
      class_test(fscanf(infile,"%le",
                        &(psd->S_vec[index_k*psd->PCA_Nnu+index_x]))!=1,                       // [10^-18 W/(m^2 Hz sr)]
                        psd->error_message,
                        "Could not read E vector at line %i in file '%s'",index_x+headlines,sd_file);
    }
  }

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

  /** Find z position */
  class_call(array_spline_hunt(psd->PCA_nu,
                               psd->PCA_Nnu,
                               nu,
                               &last_index,
                               &h,
                               &a,
                               &b,
                               psd->error_message),
             psd->error_message,
             psd->error_message);

  /** Evaluate corresponding values for the branching ratios */
  *G_T = array_interpolate_spline_hunt(psd->PCA_G_T,
                                       psd->ddPCA_G_T,
                                       last_index,
                                       last_index+1,
                                       h,a,b);
  *Y_SZ = array_interpolate_spline_hunt(psd->PCA_Y_SZ,
                                        psd->ddPCA_Y_SZ,
                                        last_index,
                                        last_index+1,
                                        h,a,b);
  *M_mu = array_interpolate_spline_hunt(psd->PCA_M_mu,
                                        psd->ddPCA_M_mu,
                                        last_index,
                                        last_index+1,
                                        h,a,b);

  for(index_k=0; index_k<psd->N_PCA; ++index_k){
    S[index_k] = array_interpolate_spline_hunt(psd->S_vec+index_k*psd->PCA_Nnu,
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
 * @param psd     Input: pointer to distortions structure (to be freed)
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
 * Outputs
 */
int heating_output_titles(struct distortions * psd, char titles[_MAXTITLESTRINGLENGTH_]){

  class_store_columntitle(titles,"Redshift z",_TRUE_);
  class_store_columntitle(titles,"Heat  [-]",_TRUE_);
  class_store_columntitle(titles,"LHeat [-]",_TRUE_);
  class_store_columntitle(titles,"EHeat [-]",_TRUE_);

  return _SUCCESS_;
}
int heating_output_data(struct distortions * psd,
                        int number_of_titles,
                        double * data){
  int storeidx;
  double * dataptr;

  for (int index_z=0; index_z<psd->z_size; index_z++) {
    dataptr = data + index_z*number_of_titles;
    storeidx = 0;
    class_store_double(dataptr, psd->z[index_z], _TRUE_, storeidx);
    class_store_double(dataptr, psd->dQrho_dz_tot[index_z], _TRUE_, storeidx);
    class_store_double(dataptr, psd->dQrho_dz_tot[index_z]*(1.+psd->z[index_z]), _TRUE_, storeidx);
    class_store_double(dataptr, psd->dQrho_dz_tot_screened[index_z]*(1.+psd->z[index_z]), _TRUE_, storeidx);
  }

  return _SUCCESS_;
}

int distortions_output_titles(struct distortions * psd, char titles[_MAXTITLESTRINGLENGTH_]){

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
int distortions_output_data(struct distortions * psd,
                            int number_of_titles,
                            double * data){
  int index_type;
  int storeidx;
  double * dataptr;

  for (int index_x=0; index_x<psd->x_size; index_x++) {
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



