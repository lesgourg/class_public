/** @file distortions.c Documented module on spectral distortions
 * Matteo Lucca, 31.10.2018
 * Nils Schoeneberg, 18.02.2019
 */

#include "distortions.h"

/**
 * Initialize the distortions structure.
 *
 * @param pba        Input: pointer to background structure
 * @param ppt        Input: pointer to the perturbations structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param ppm        Input: pointer to the primordial structure
 * @param psd        Input/Output: pointer to initialized distortions structure
 * @return the error status
 */
int distortions_init(struct precision * ppr,
                     struct background * pba,
                     struct perturbs * ppt,
                     struct thermo * pth,
                     struct primordial * ppm,
                     struct distortions * psd) {
  /* Define local variables */
  double * pvecheat;
  double * pvecdist;
  int last_index = 0;
  int index_br;
  int index_x;

  if(psd->has_distortions == _FALSE_){
    return _SUCCESS_;
  }
  if (psd->distortions_verbose > 0) {
    printf("Computing spectral distortions \n");
  }

  psd->distortions_verbose = 2; //TODO :: REMOVE

  class_test(pth->compute_damping_scale==_FALSE_,
             psd->error_message,
             "Cannot compute spectral distortions without damping scale \n");

  /** 
   * Define global quanties, e.g. indeces, z and x arrays and branching ratios
   */ 

  /* Assign values to all indices in the distortions structure */
  class_call(distortions_indices(psd),
             psd->error_message,
             psd->error_message);

  /* Define z and x arrays TODO :: convert to precision parameters the limits and N */
  class_call(distortions_get_xz_lists(psd),
             psd->error_message,
             psd->error_message);

  /* Define visibility function */
  class_call(distortions_visibility_function(pba,pth,psd),
             psd->error_message,
             psd->error_message);

  /* Define branching ratios */
  class_call(distortions_branching_ratios(ppr,psd),
             psd->error_message,
             psd->error_message);

  /* TODO :: what to do with greens function data ?
  class_call(distortions_read_Greens_data(ppr,psd),
             psd->error_message,
             psd->error_message); */

  /**
   * Define heting function
   */

  /* Allocate space for heating output parameters */
  class_alloc(psd->dQrho_dz_tot,
              psd->z_size*sizeof(double),
              psd->error_message);
  /* Get the heating TODO :: get from thermodynamics after exoclass hyrec recfast merge */
  class_alloc(pvecheat,
              psd->ht_size*sizeof(double),
              psd->error_message);
  for (int index_z = 0; index_z < psd->z_size; index_z++) {
    class_call(heating_at_z(ppr,pba,ppt,pth,ppm,psd,
                            psd->z[index_z],
                            pvecheat),
               psd->error_message,
               psd->error_message);
    /* Public quantities */
    psd->dQrho_dz_tot[index_z] = pvecheat[psd->index_ht_dQrho_dz_tot]*psd->bb_visibility_function[index_z];
  }
  /* Free space */
  free(pvecheat);

  /**
   * Define spectral distortion amplitudes
   * This boils down to the integral int ( f(z) d(Q/rho)(z)/dz )
   */

  /* Allocate space for spectral distortion amplitude parameters */
  class_alloc(psd->sd_parameter,
              psd->br_size*sizeof(double),
              psd->error_message);

  for(int index_br=0;index_br<psd->br_size;++index_br){
    class_call(array_trapezoidal_convolution(psd->branching_ratios[index_br],
                                             psd->dQrho_dz_tot,
                                             psd->z_size,
                                             psd->z_weights,
                                             &(psd->sd_parameter[index_br]),
                                             psd->error_message),
               psd->error_message,
               psd->error_message);
  }

  /* Small short-hand notations of the parameters stored in this array */
  psd->g = psd->sd_parameter[psd->index_br_f_g];
  psd->mu = psd->sd_parameter[psd->index_br_f_mu];
  psd->y = psd->sd_parameter[psd->index_br_f_y];
  psd->r = psd->sd_parameter[psd->index_br_f_r];

  /* Include additional sources of distortions (see also Chluba 2016 for useful discussion) */
  psd->y += 2.525e-7;   // CMB Dipole (Chluba & Sunyaev 2004)
  psd->y += 4.59e-13;   // CMB Quadrupole (Chluba & Sunyaev 2004)
  psd->y += 1.77e-6;    // Reionization and structure formation (Hill et al. 2015)

  /* Calculate total heating */
  psd->Drho_over_rho = 4.*psd->g+psd->mu/1.401+4.*psd->y+psd->r;

  /* Print found parameters */
  if (psd->distortions_verbose > 1) { printf("-> total injected/extracted heat = %g\n", psd->Drho_over_rho); }

  if (psd->distortions_verbose > 1) {
    if (psd->mu > 9.e-5) { printf("-> mu-parameter = %g. WARNING: The value of your mu-parameter is larger than the FIRAS constraint mu<9e-5.\n", psd->mu); }
    else{ printf("-> mu-parameter = %g\n", psd->mu); }
  }

  if (psd->distortions_verbose > 1) {
    if (psd->y>1.5e-5) { printf("-> y-parameter = %g. WARNING: The value of your y-parameter is larger than the FIRAS constraint y<1.5e-5.\n", psd->y); }
    else{ printf("-> y-parameter = %g\n", psd->y); }
  }

  if (psd->distortions_verbose > 1) { printf("-> r-parameter = %g\n", psd->r); }

  /**
   * Define final spectral distortions
   */

  /* Allocate space for distortions output parameters */
  class_alloc(psd->DI,
              psd->x_size*sizeof(double),
              psd->error_message);
  /* Get the distortions */
  class_alloc(pvecdist,
              psd->sd_size*sizeof(double),
              psd->error_message);
  for (index_x = 0; index_x<psd->x_size; index_x++) {
    class_call(distortions_at_x(pba,psd,
                                psd->x[index_x],
                                pvecdist),
               psd->error_message,
               psd->error_message);
    /* Public quantities */
    psd->DI[index_x] = pvecdist[psd->index_sd_DI];
  }
  /* Free space */
  free(pvecdist);

  return _SUCCESS_;
}


/**
 * Free all memory space allocated by distortions_init()
 *
 * @param psd     Input: pointer to distortions structure (to be freed)
 * @return the error status
 */
int distortions_free(struct distortions * psd) {
  /* Define local variables */
  int index_br;

  if(psd->has_distortions == _TRUE_){
    /* Free from distortions_get_xz_lists() */
    free(psd->z);
    free(psd->x);
    free(psd->z_weights);

    /* Free from distortions_visibility_function() */
    free(psd->bb_visibility_function);

    /* Free from distortions_branching_ratios() */
    for(index_br=0;index_br<psd->br_size;++index_br){
      free(psd->branching_ratios[index_br]);
    }
    free(psd->branching_ratios);

    /* Free from distortions_init() */
    free(psd->dQrho_dz_tot);
    free(psd->sd_parameter);
    free(psd->DI);
  }

  return _SUCCESS_;
}


/**
 * Assign value to each relevant index in vectors of distortions quantities.
 *
 * @param psd     Input: pointer to distortions structure
 * @return the error status
 */
int distortions_indices(struct distortions * psd) {
  /* Define local variables */
  int index_ht = 0;
  int index_sd = 0;
  int index_br = 0;

  /* Branching ratios from distortions_branching_ratios() */
  class_define_index(psd->index_br_f_g,_TRUE_,index_br,1);
  class_define_index(psd->index_br_f_mu,_TRUE_,index_br,1);
  class_define_index(psd->index_br_f_y,_TRUE_,index_br,1);
  class_define_index(psd->index_br_f_r,_TRUE_,index_br,1);

  psd->br_size = index_br;

  /* Heating functions from heating_at_z() */
  class_define_index(psd->index_ht_dQrho_dz_cool,_TRUE_,index_ht,1);
  class_define_index(psd->index_ht_dQrho_dz_diss,_TRUE_,index_ht,1);
  class_define_index(psd->index_ht_dQrho_dz_ann,_TRUE_,index_ht,1);
  class_define_index(psd->index_ht_dQrho_dz_dec,_TRUE_,index_ht,1);
  class_define_index(psd->index_ht_dQrho_dz_eva_PBH,_TRUE_,index_ht,1);
  class_define_index(psd->index_ht_dQrho_dz_acc_PBH,_TRUE_,index_ht,1);
  class_define_index(psd->index_ht_dQrho_dz_tot,_TRUE_,index_ht,1);

  psd->ht_size = index_ht;

  /* Spectral distortions from distortions_at_x() */
  class_define_index(psd->index_sd_Y,_TRUE_,index_sd,1);
  class_define_index(psd->index_sd_M,_TRUE_,index_sd,1);
  class_define_index(psd->index_sd_G,_TRUE_,index_sd,1);
  class_define_index(psd->index_sd_DI,_TRUE_,index_sd,1);

  psd->sd_size = index_sd;

  return _SUCCESS_;
}


/**
 * Compute redshift and frequency vectors and weights for redshift integral.
 *
 * @param psd     Input: pointer to distortions structure
 * @return the error status
 */
int distortions_get_xz_lists(struct distortions* psd){
  /* Define local variables */
  int index_z;
  int index_x;

  psd->z_min = 1.e3;
  psd->z_max = 5.e6;
  psd->z_size = (int) 550;
  psd->z_delta = (log(psd->z_max)-log(psd->z_min))/psd->z_size;

  psd->x_min = 1.e-2;
  psd->x_max = 5.e1;
  psd->x_size = (int) 550;
  psd->x_delta = (log(psd->x_max)-log(psd->x_min))/psd->x_size;

  class_alloc(psd->z,
              psd->z_size*sizeof(double),
              psd->error_message);
  class_alloc(psd->x,
              psd->x_size*sizeof(double),
              psd->error_message);
  class_alloc(psd->z_weights,
              psd->z_size*sizeof(double),
              psd->error_message);

  for (index_z = 0; index_z < psd->z_size; index_z++) {
    psd->z[index_z] = exp(log(psd->z_min)+psd->z_delta*index_z);
  }
  for (index_x = 0; index_x<psd->x_size; index_x++) {
    psd->x[index_x] = exp(log(psd->x_min)+psd->x_delta*index_x);
  }
  class_call(array_trapezoidal_weights(
                    psd->z,
                    psd->z_size,
                    psd->z_weights,
                    psd->error_message),
             psd->error_message,
             psd->error_message);

  return _SUCCESS_;
}


/**
 * Compute visibility function.
 *
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */
int distortions_visibility_function(struct background * pba,
                                    struct thermo * pth,
                                    struct distortions * psd){
  /* Define local variables */
  int index_z;
  double z,z_th,z_muy;

  class_alloc(psd->bb_visibility_function,
              psd->z_size*sizeof(double),
              psd->error_message);

  for(index_z=0;index_z<psd->z_size;++index_z){
    z = psd->z[index_z];

    z_muy = 5.e4;
    z_th = 1.98e6*
           pow((1.-pth->YHe/2.)/0.8767,-2./5.)*
           pow(pba->Omega0_b*pow(pba->h,2.)/0.02225,-2./5.)*
           pow(pba->T_cmb/2.726,1./5.);

    psd->bb_visibility_function[index_z] = exp(-pow(z/z_th,2.5));
  }

  return _SUCCESS_;
}


/**
 * Calculate all the branching ratios of the Green's functions.
 *
 * Computing the full evolution of the thermal history of the universe is rather time consuming
 * and mathematically challenging. It is therefore not implemented here. However, there are
 * (at least) 5 levels of possible approximatin to evaluate the SD branching ratios (see also
 * Chluba 2016 for useful discussion)
 *    1) Use a sharp transition at z_mu-y and no distortions before z_th
 *    2) Use a sharp transition at z_mu-y and a soft transition at z_th
 *    3) Use a soft transition at a_mu-y and z_th as described in Chluba 2013
 *    4) Use a soft transition at a_mu-y and z_th imposing conservation of energy
 *    5) Use a PCA method as described in Chluba & Jeong 2014 (TODO)
 * The user can select the preferred option with 'branching approx'=sharp_sharp for 1),
 * =sharp_soft for 2), =soft_soft for 3), =soft_soft_cons for 4) and =exact for 5)
 * (default =soft_soft).
 *
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */
int distortions_branching_ratios(struct precision * ppr,
                                 struct distortions* psd){
  /* Define local variables */
  int index_z,index_br;
  double f_g,f_mu,f_y,f_r,f;
  double z_th,z_muy;
  double z;

  class_alloc(psd->branching_ratios,
              psd->br_size*sizeof(double*),
              psd->error_message);
  for(index_br=0;index_br<psd->br_size;++index_br){
    class_alloc(psd->branching_ratios[index_br],
                psd->z_size*sizeof(double),
                psd->error_message);
  }

  for(index_z=0; index_z<psd->z_size; ++index_z){
    z = psd->z[index_z];
    f = psd->bb_visibility_function[index_z];

    /** 1) Calculate branching ratios using sharp_sharp transition */
    if(psd->branching_approx == 1){
      if(z>z_th){
        f_g = 1.;
        f_mu = 0.;
        f_y = 0.;
      }
      if(z<z_th && z>z_muy){
        f_g = 0.;
        f_mu = 1.;
        f_y = 0.;
      }
      if(z<z_muy){
        f_g = 0.;
        f_mu = 0.;
        f_y = 1.;
      }
    }

    /** 2) Calculate branching ratios using sharp_soft transition */
    if(psd->branching_approx == 2){
      f_g = 1.-f;
      if(z>z_muy){
        f_mu = f;
        f_y = 0.;
      }
      if(z<z_muy){
        f_mu = 0.;
        f_y = 1.;
      }
    }

    /** 3) Calculate branching ratios unsing soft_soft transitions */
    if(psd->branching_approx == 3){
      f_g = 1.-f;
      f_mu = f*(1.0-exp(-pow((1.0+z)/(5.8e4),1.88)));
      f_y = 1.0/(1.0+pow((1.0+z)/(6.0e4),2.58));
    }

    /** 4) Calculate branching ratios unsing soft_soft_cons transitions */
    if(psd->branching_approx == 4){
      f_g = 1.-f;
      f_y = 1.0/(1.0+pow((1.0+z)/(6.0e4),2.58));
      f_mu = f*(1.-f_y);

    }

    /** 5) Calculate branching ratios according to Chluba & Jeong 2014 */
    if(psd->branching_approx == 5){
      class_call(distortions_read_BR_exact_data(ppr,psd),
                 psd->error_message,
                 psd->error_message);
      class_alloc(psd->br_exact_table, psd->br_exact_Nz*sizeof(double),psd->error_message);
      class_call(array_spline_table_lines(psd->br_exact_z,
                                          psd->br_exact_Nz,
                                          psd->f_g_exact,
                                          1,
                                          psd->br_exact_table,
                                          _SPLINE_NATURAL_,
                                          psd->error_message),
                 psd->error_message,
                 psd->error_message);
      class_call(array_interpolate_spline(psd->br_exact_z,
                                          psd->br_exact_Nz,
                                          psd->f_g_exact,
                                          psd->br_exact_table,
                                          1,
                                          z,
                                          0,
                                          &f_g,
                                          1,
                                          psd->error_message),
                 psd->error_message,
                 psd->error_message);

      free(psd->br_exact_table);
      class_call(distortions_free_BR_exact_data(psd),
                 psd->error_message,
                 psd->error_message);

    }

    f_r = 1.-f_g-f_mu-f_y;
    
    psd->branching_ratios[psd->index_br_f_g][index_z]=f_g/4.;
    psd->branching_ratios[psd->index_br_f_mu][index_z]=f_mu*1.401;
    psd->branching_ratios[psd->index_br_f_y][index_z]=f_y/4.;
    psd->branching_ratios[psd->index_br_f_r][index_z]=f_r;
  }

  return _SUCCESS_;
}
/**
 * Calculate all redshift dependent quantities needed to compute the spectral distortions, i.e.
 * the branching ratios of the Green's functions and the heating rates.
 *
 * DOES NOT INCLUDE THE VISIBILITY FUNCTION CURRENTLY
 *
 * There are many possible sources of heating (for all details see e.g. Chluba & Sunyaev 2012),
 * some are present even for the standard cosmological model like
 *    1) Adiabatically cooling electrons and barions as described in Chluba & Sunyaev 2012
 *        (see also Khatri, Sunyaev & Chluba 2012 for useful discussion)
 *    2) Dissipation of acoustic waves as described in Chluba, Khatri & Sunyaev 2012 (see also
 *       Chluba 2013 and Diacoumis & Wong 2017 for useful discussion) with two possible
 *       approximations
 *          a) Eq. 42 from Chluba, Khatri & Sunyaev 2012 (approximated to Y_SZ*S_ac) (TODO)
 *          b) Eq. 45 from Chluba, Khatri & Sunyaev 2012
 *       The user can select the preferred option with 'heating approx'=yes/no (default =yes)
 *    3) Cosmological recombination radiation as described in Chluba & Ali-Haimoud 2016  (TODO)
 * while some other are related to new possible physical processes like
 *    4) Annihilating particles (e.g. dark matter) as described in Chluba 2010 and Chluba
 *       & Sunyaev 2012. (see also Chluba 2013 for useful discussion)
 *    5) Decaying relic particles as described in Chluba 2010 and Chluba & Sunyaev 2012
 *       (see also Chluba 2013 for useful discussion)
 *    6) Evaporation of primordial black holes as described in Poulin et al. 2017 (see also
 *       Tashiro & Sugiyama 2008, Carr et al. 2010 and Carr et al. 2016 for useful discussions)
 *    7) Acctretion of matter into primordial black holes both via
 *          a) Spherical accretion as described in Ali-Haimoud & Kamionkowski 2017 (note: "we
 *             find that CMB spectral distortion measurements, both current and upcoming,
 *             do not place any constraints on PBHs.") and
 *          b) Disk accretion as described in Poulin et al. 2017
 *       (see also Carr et al. 2010 for useful discussion)
 *
 * @param pba        Input: pointer to background structure
 * @param ppt        Input: pointer to the perturbations structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param ppm        Input: pointer to the primordial structure
 * @param psd        Input: pointer to the distortions structure
 * @param z          Input: redshift
 * @param pvecheat   Output: vector of heating functions (assumed to be already allocated)
 */
int heating_at_z(struct precision * ppr,
                 struct background* pba,
                 struct perturbs * ppt,
                 struct thermo * pth,
                 struct primordial * ppm,
                 struct distortions * psd,
                 double z,
                 double * pvecheat) {
  /* Definitions */
  double tau;
  int last_index = 0;
  double * pvecback, O_b, O_cdm, h, H, a, t, rho_g, R, T_g0;
  double * pvecthermo, dk, dz_kD, kD, N_e, X_e, Y_He, Y_p, p_ann;
  double z_th, z_muy;
  double alpha_h, tilde_rho_g, theta_g;
  double k_max, k_min, k_size, k_delta, k, pk_primordial_k, * int_dQrho_dz_diss_full, * int_dQrho_dz_diss_approx;
  double g_h;

  /* From z to tau */
  class_call(background_tau_of_z(pba,
                                 z,
                                 &tau),
             pba->error_message,
             psd->error_message);

  /* Import quantities from background */
  class_alloc(pvecback,
              pba->bg_size*sizeof(double),
              psd->error_message);
  class_call(background_at_tau(pba,
                               tau,
                               pba->long_info,
                               pba->inter_closeby,
                               &last_index,
                               pvecback),
             pba->error_message,
             psd->error_message);
  O_b = pba->Omega0_b;                                                                // [-]
  O_cdm = pba->Omega0_cdm;                                                            // [-]
  h = pba->h;                                                                         // [-]
  H = pvecback[pba->index_bg_H];                                                      // [1/Mpc]
  a = pvecback[pba->index_bg_a];                                                      // [-]
  t = pvecback[pba->index_bg_time];                                                   // [Mpc]
  t /= _s_over_Mpc_;                                                                  // [s]
  rho_g = pvecback[pba->index_bg_rho_g];                                              // [1/Mpc^4]
  rho_g /= _GeVcm3_over_Mpc4_;                                                        // [GeV/cm^3]
  R = (3./4.)*pvecback[pba->index_bg_rho_b]/pvecback[pba->index_bg_rho_g];            // [-]
  T_g0 = pba->T_cmb;                                                                  // [K]

  /* Import quantities from thermodynamics */
  class_alloc(pvecthermo,
              pth->tt_size*sizeof(double),
              psd->error_message);
  class_call(thermodynamics_at_z(pba,
                                 pth,
                                 z,
                                 pth->inter_normal,
                                 &last_index,
                                 pvecback,
                                 pvecthermo),
             pth->error_message,
             psd->error_message);
  dk = pvecthermo[pth->index_th_dkappa];                                              // [1/Mpc]
  dz_kD = (1./(H*dk))*(16.0/15.0+pow(R,2.0)/(1.0+R))/(6.0*(1.0+R));                   // [Mpc^2]
  kD = 2.*_PI_/pvecthermo[pth->index_th_r_d];                                         // [1/Mpc]
  N_e = pth->n_e;                                                                     // [1/m^3] (today)
  X_e = pvecthermo[pth->index_th_xe];                                                 // [-]
  Y_He = pth->YHe;                                                                    // [-]
  Y_p = 1.-Y_He;                                                                      // [-]

  /* Free allocated space */
  free(pvecback);
  free(pvecthermo);

  /** 1) Adiabatically cooling electrons and barions */
  tilde_rho_g = rho_g/(_m_e_/_GeV_over_kg_);                                          // [1/cm^3]
  theta_g = (_k_B_*T_g0*(1.+z))/(_m_e_*pow(_c_,2.));                                  // [-]
  alpha_h = (3./2.)*N_e*1.e-6*pow(1.+z,3.)*(1.+Y_He+X_e);                             // [1/cm^3]
  pvecheat[psd->index_ht_dQrho_dz_cool] = -a*alpha_h/tilde_rho_g*theta_g;             // [-]

  /** 2) dissipation of acoustic waves */
  /* a) Full function */
  if (psd->dQrho_dz_diss_approx == _FALSE_){
  }

  /* b) Approximated function */
  if (psd->dQrho_dz_diss_approx == _TRUE_){
    k_max = 5.*kD;
    k_min = 0.12;
    k_size = 500;        /* Found to be reasonable for this particular integral */

    class_alloc(int_dQrho_dz_diss_approx,
                k_size*sizeof(double),
                psd->error_message);

    for (int index_k=0; index_k<k_size; index_k++) {
      k = exp(log(k_min)+(log(k_max)-log(k_min))/(k_size)*index_k);

      /* Import quantities from primordial
       * Note that the the heating caused by dissipation of acustic waves depends on the primordial
       * power spectrum and to analyse the its influence it is enough to change initial parameters. */
      class_call(primordial_spectrum_at_k(ppm,
                                          ppt->index_md_scalars,
                                          linear,
                                          k,
                                          &pk_primordial_k),
                 ppm->error_message,
                 psd->error_message);
      /* Define integrand for approximated function */
      int_dQrho_dz_diss_approx[index_k] = 4.*0.81*
                                          pow(k,2.)*
                                          pk_primordial_k*
                                          dz_kD*exp(-2.*pow(k/kD,2.));                // [-]
    }

    /* Integrate approximate function */
    class_call(simpson_integration(k_size,
                                   int_dQrho_dz_diss_approx,
                                   (log(k_max)-log(k_min))/(k_size),
                                   &pvecheat[psd->index_ht_dQrho_dz_diss],
                                   psd->error_message),
               psd->error_message,
               psd->error_message);

    /* Free space */
    free(int_dQrho_dz_diss_approx);
  }

  /** 3) Cosmological recombination */

  /** 4) Annihilating particles */
  g_h = (1.+Y_He+2*(X_e))/(3.*(1+Y_He));                                              // [-] (TODO: INCOMPLETE)

  pvecheat[psd->index_ht_dQrho_dz_ann] = 0.;
  /*
  if(pth->annihilation != 0.){
      class_call(thermodynamics_DM_annihilation_energy_injection(ppr,
                                                                 pba,
                                                                 pth,
                                                                 z,
                                                                 &pvecheat[psd->index_ht_dQrho_dz_ann],
                                                                 pth->error_message),
                 pth->error_message,
                 psd->error_message);                                               // [J/(m^3 s)]
      pvecheat[psd->index_ht_dQrho_dz_ann] *= a/(H*_s_over_Mpc_)*g_h/
                                              (rho_g/(_eV_over_joules_*1.e-9)*1.e6);

  }
  */

  /** 5) Decaying relic particles */
  pvecheat[psd->index_ht_dQrho_dz_dec] = 0.;
  /*
  if(pth->decay_fraction != 0.){
      class_call(thermodynamics_DM_decay_energy_injection(ppr,
                                                          pba,
                                                          pth,
                                                          z,
                                                          &pvecheat[psd->index_ht_dQrho_dz_dec],
                                                          pth->error_message),
                 pth->error_message,
                 psd->error_message);                                               // [J/(m^3 s)]
      pvecheat[psd->index_ht_dQrho_dz_dec] *= a/(H*_s_over_Mpc_)*g_h/
                                              (rho_g/(_eV_over_joules_*1.e-9)*1.e6);
  }
  */

  /** 6) Evaporation of primordial black holes */
  /* Note: to use this energy injection mechanism you need to set recfast_z_initial=5.1e6 and recfast_Nz0=4000000 in precision.h */
  pvecheat[psd->index_ht_dQrho_dz_eva_PBH] = 0.;
  /*
  if(pth->PBH_evaporating_mass != 0.){
      class_call(thermodynamics_evaporating_pbh_energy_injection(ppr,
                                                                 pba,
                                                                 pth,
                                                                 z,
                                                                 &pvecheat[psd->index_ht_dQrho_dz_eva_PBH],
                                                                 pth->error_message),
                 pth->error_message,
                 psd->error_message);                                               // [J/(m^3 s)]
      pvecheat[psd->index_ht_dQrho_dz_eva_PBH] *= a/(H*_s_over_Mpc_)*g_h/
                                                  (rho_g/(_eV_over_joules_*1.e-9)*1.e6);
  }
  */

  /** 7) Accretion of matter into primordial black holes */
  /* Note: To use this energy injection mechanism you need to set recfast_z_initial=5.1e6 and recfast_Nz0=4000000 in precision.h */
  pvecheat[psd->index_ht_dQrho_dz_acc_PBH] = 0.;
  /*
  if(pth->PBH_accreting_mass != 0.){
      class_call(thermodynamics_accreting_pbh_energy_injection(ppr,
                                                               pba,
                                                               pth,
                                                               z,
                                                               &pvecheat[psd->index_ht_dQrho_dz_acc_PBH],
                                                               pth->error_message),
                 pth->error_message,
                 psd->error_message);                                               // [J/(m^3 s)]
      pvecheat[psd->index_ht_dQrho_dz_acc_PBH] *= a/(H*_s_over_Mpc_)*g_h/
                                                  (rho_g/(_eV_over_joules_*1.e-9)*1.e6);
  }
  */

  /** Total heating rate */
  pvecheat[psd->index_ht_dQrho_dz_tot] = 0.;
  pvecheat[psd->index_ht_dQrho_dz_tot] += pvecheat[psd->index_ht_dQrho_dz_cool];    // [-]
  pvecheat[psd->index_ht_dQrho_dz_tot] += pvecheat[psd->index_ht_dQrho_dz_diss];    // [-]
  pvecheat[psd->index_ht_dQrho_dz_tot] += pvecheat[psd->index_ht_dQrho_dz_ann];     // [-]
  pvecheat[psd->index_ht_dQrho_dz_tot] += pvecheat[psd->index_ht_dQrho_dz_dec];     // [-]
  pvecheat[psd->index_ht_dQrho_dz_tot] += pvecheat[psd->index_ht_dQrho_dz_eva_PBH]; // [-]
  pvecheat[psd->index_ht_dQrho_dz_tot] += pvecheat[psd->index_ht_dQrho_dz_acc_PBH]; // [-]

  return _SUCCESS_;
}


/**
 * Evaluation of the spectral distortions.
 *
 * @param pba        Input: pointer to background structure
 * @param ppt        Input: pointer to the perturbations structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param ppm        Input: pointer to the primordial structure
 * @param psd        Input: pointer to the distortions structure
 * @param z          Input: redshift
 * @param pvecdist   Output: vector of distortions functions (assumed to be already allocated)
 */
int distortions_at_x(struct background* pba,
                     struct distortions * psd,
                     double x,
                     double * pvecdist) {

  /* Calculate spectral distortions */
  pvecdist[psd->index_sd_Y] = 2.*pow(_k_B_*pba->T_cmb, 3.)/pow(_h_P_*_c_,2.)*
                              pow(x,4.)*
                              exp(-x)/pow(1.-exp(-x),2.)*
                              (x*(1.+exp(-x))/(1.-exp(-x))-4.);                     // [W/(m^2 Hz)]
  pvecdist[psd->index_sd_M] = 2.*pow(_k_B_*pba->T_cmb, 3.)/pow(_h_P_*_c_,2.)*
                              pow(x,4.)*
                              exp(-x)/pow(1.-exp(-x),2.)*
                              (1./2.19229-1./x);                                    // [W/(m^2 Hz)]
  pvecdist[psd->index_sd_G] = 2.*pow(_k_B_*pba->T_cmb, 3.)/pow(_h_P_*_c_,2.)*
                              pow(x,4.)*
                              exp(-x)/pow(1.-exp(-x),2.);                           // [W/(m^2 Hz)]

  pvecdist[psd->index_sd_DI] = psd->y*pvecdist[psd->index_sd_Y]+
                               psd->mu*pvecdist[psd->index_sd_M]+
                               psd->g*pvecdist[psd->index_sd_G];                    // [W/(m^2 Hz)]

  return _SUCCESS_;
}


/**
 * Reads the external file Greens_data copied from CosmoTherm by Jens Chluba
 *
 * @param ppr        Input: pointer to precision structure
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */
int distortions_read_Greens_data(struct precision * ppr,
                                 struct distortions * psd){

  FILE * infile;
  char line[_LINE_LENGTH_MAX_];
  char * left;
  int headlines = 0;

  class_open(infile, ppr->Greens_file, "r",
             psd->error_message);

  psd->Greens_Nz = 0;
  psd->Greens_Nx = 0;
  while (fgets(line,_LINE_LENGTH_MAX_-1,infile) != NULL) {
    headlines++;

    /* Eliminate blank spaces at beginning of line */
    left=line;
    while (left[0]==' ') {
        left++;
    }
    if (left[0] > 39) {
      /* read number of lines, infer size of arrays and allocate them */
      class_test(sscanf(line, "%d %d", &psd->Greens_Nz,&psd->Greens_Nx) != 2,
                 psd->error_message,
                 "could not header (number of columns, number of lines) at line %i in file '%s' \n",headlines,ppr->Greens_file);

      class_alloc(psd->Greens_z, psd->Greens_Nz*sizeof(double), psd->error_message);
      class_alloc(psd->Greens_T_ini, psd->Greens_Nz*sizeof(double), psd->error_message);
      class_alloc(psd->Greens_T_last, psd->Greens_Nz*sizeof(double), psd->error_message);
      class_alloc(psd->Greens_rho, psd->Greens_Nz*sizeof(double), psd->error_message);

      class_alloc(psd->Greens_x, psd->Greens_Nx*sizeof(double), psd->error_message);
      class_alloc(psd->Greens_blackbody, psd->Greens_Nx*sizeof(double), psd->error_message);
      class_alloc(psd->Greens_function, psd->Greens_Nz*psd->Greens_Nx*sizeof(double), psd->error_message);
      break;
    }
  }
  for(int i=0;i<psd->Greens_Nz;++i){
    class_test(fscanf(infile,"%le",&(psd->Greens_z[i]))!=1,
                      psd->error_message,
                      "Could not read z values at line %i in file '%s'",headlines+1,ppr->Greens_file);
  }
  for(int i=0;i<psd->Greens_Nz;++i){
    class_test(fscanf(infile,"%le",&(psd->Greens_T_ini[i]))!=1,
                      psd->error_message,
                      "Could not read T_ini values at line %i in file '%s'",headlines+2,ppr->Greens_file);
  }
  for(int i=0;i<psd->Greens_Nz;++i){
    class_test(fscanf(infile,"%le",&(psd->Greens_T_last[i]))!=1,
                      psd->error_message,
                      "Could not read T_last values at line %i in file '%s'",headlines+3,ppr->Greens_file);
  }
  for(int i=0;i<psd->Greens_Nz;++i){
    class_test(fscanf(infile,"%le",&(psd->Greens_rho[i]))!=1,
                      psd->error_message,
                      "Could not read rho values at line %i in file '%s'",headlines+4,ppr->Greens_file);
  }
  for(int i=0;i<psd->Greens_Nx;++i){
    class_test(fscanf(infile,"%le",&(psd->Greens_x[i]))!=1,
                      psd->error_message,
                      "Could not read Greens x at line %i in file '%s'",i+headlines+5,ppr->Greens_file);
    for(int j=0;j<psd->Greens_Nz;++j){
      class_test(fscanf(infile,"%le",&(psd->Greens_function[i*psd->Greens_Nz+j]))!=1,
                        psd->error_message,
                        "Could not read Greens function at line %i in file '%s'",i+headlines+5,ppr->Greens_file);
    }
    class_test(fscanf(infile,"%le",&(psd->Greens_blackbody[i]))!=1,
                      psd->error_message,
                      "Could not read Greens blackbody at line %i in file '%s'",i+headlines+5,ppr->Greens_file);
  }

  return _SUCCESS_;

}


/**
 * Free from distortions_read_Greens_data()
 *
 * @param psd     Input: pointer to distortions structure (to be freed)
 * @return the error status
 */
int distortions_free_Greens_data(struct distortions * psd){
    free(psd->Greens_z);
    free(psd->Greens_x);
    free(psd->Greens_T_ini);
    free(psd->Greens_T_last);
    free(psd->Greens_rho);
    free(psd->Greens_blackbody);
    free(psd->Greens_function);
}



/**
 * Reads the external file branching_ratios_exact calculated according to Chluba & Jeong 2014
 * for PIXIE like configuration (frequency in interval [30,1000] GHz and bin width of 1 GHz)
 *
 * @param ppr        Input: pointer to precision structure
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */
int distortions_read_BR_exact_data(struct precision * ppr,
                                   struct distortions * psd){

  FILE * infile;
  char line[_LINE_LENGTH_MAX_];
  char * left;
  int headlines = 0;

  class_open(infile, ppr->br_exact_file, "r",
             psd->error_message);

  psd->br_exact_Nz = 0;
  psd->br_exact_N_columns = 0;
  while (fgets(line,_LINE_LENGTH_MAX_-1,infile) != NULL) {
    headlines++;

    /* Eliminate blank spaces at beginning of line */
    left=line;
    while (left[0]==' ') {
        left++;
    }
    if (left[0] > 39) {
      /* read number of lines, infer size of arrays and allocate them */
      class_test(sscanf(line, "%d %d %d", &psd->br_exact_Nz,&psd->br_exact_N_columns, &psd->index_e) != 3,
                 psd->error_message,
                 "could not header (number of lines, number of columns, number of multipoles) at line %i in file '%s' \n",headlines,ppr->br_exact_file);

      class_alloc(psd->br_exact_z, psd->br_exact_Nz*sizeof(double), psd->error_message);
      class_alloc(psd->f_g_exact, psd->br_exact_Nz*sizeof(double), psd->error_message);
      class_alloc(psd->f_y_exact, psd->br_exact_Nz*sizeof(double), psd->error_message);
      class_alloc(psd->f_mu_exact, psd->br_exact_Nz*sizeof(double), psd->error_message);

      class_alloc(psd->E_vec, psd->br_exact_Nz*psd->index_e*sizeof(double), psd->error_message);
      break;
    }
  }
  for(int i=0; i<psd->br_exact_Nz; ++i){
    class_test(fscanf(infile,"%le",&(psd->br_exact_z[i]))!=1,
                      psd->error_message,
                      "Could not read x at line %i in file '%s'",i+headlines+5,ppr->br_exact_file);
    class_test(fscanf(infile,"%le",&(psd->f_g_exact[i]))!=1,
                      psd->error_message,
                      "Could not read J_T at line %i in file '%s'",i+headlines+5,ppr->br_exact_file);
    class_test(fscanf(infile,"%le",&(psd->f_y_exact[i]))!=1,
                      psd->error_message,
                      "Could not read J_y at line %i in file '%s'",i+headlines+5,ppr->br_exact_file);
    class_test(fscanf(infile,"%le",&(psd->f_mu_exact[i]))!=1,
                      psd->error_message,
                      "Could not read J_mu at line %i in file '%s'",i+headlines+5,ppr->br_exact_file);
    for(int j=0; j<psd->index_e; ++j){
      class_test(fscanf(infile,"%le",&(psd->E_vec[i*psd->index_e+j]))!=1,
                        psd->error_message,
                        "Could not read S vector at line %i in file '%s'",i+headlines+5,ppr->br_exact_file);
    }
  }

  return _SUCCESS_;

}


/**
 * Free from distortions_read_BR_exact_data()
 *
 * @param psd     Input: pointer to distortions structure (to be freed)
 * @return the error status
 */
int distortions_free_BR_exact_data(struct distortions * psd){
    /* Free from distortions_read_BR_exact_data() */
    free(psd->br_exact_z);
    free(psd->f_g_exact);
    free(psd->f_y_exact);
    free(psd->f_mu_exact);
    free(psd->E_vec);
}


/**
 * Reads the external file PCA_distortions_shapes calculated according to Chluba & Jeong 2014
 * for PIXIE like configuration (frequency in interval [30,1000] GHz and bin width of 1 GHz)
 *
 * @param ppr        Input: pointer to precision structure
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */
int distortions_read_PCA_dist_shapes_data(struct precision * ppr,
                                          struct distortions * psd){

  FILE * infile;
  char line[_LINE_LENGTH_MAX_];
  char * left;
  int headlines = 0;

  class_open(infile, ppr->br_exact_file, "r",
             psd->error_message);

  psd->PCA_Nx = 0;
  psd->PCA_N_columns = 0;
  while (fgets(line,_LINE_LENGTH_MAX_-1,infile) != NULL) {
    headlines++;

    /* Eliminate blank spaces at beginning of line */
    left=line;
    while (left[0]==' ') {
        left++;
    }
    if (left[0] > 39) {
      /* read number of lines, infer size of arrays and allocate them */
      class_test(sscanf(line, "%d %d %d", &psd->PCA_Nx,&psd->PCA_N_columns, &psd->index_s) != 3,
                 psd->error_message,
                 "could not header (number of lines, number of columns, number of multipoles) at line %i in file '%s' \n",headlines,ppr->br_exact_file);

      class_alloc(psd->PCA_x, psd->PCA_Nx*sizeof(double), psd->error_message);
      class_alloc(psd->PCA_J_T, psd->PCA_Nx*sizeof(double), psd->error_message);
      class_alloc(psd->PCA_J_y, psd->PCA_Nx*sizeof(double), psd->error_message);
      class_alloc(psd->PCA_J_mu, psd->PCA_Nx*sizeof(double), psd->error_message);

      class_alloc(psd->S_vec, psd->PCA_Nx*psd->index_s*sizeof(double), psd->error_message);
      break;
    }
  }
  for(int i=0; i<psd->PCA_Nx; ++i){
    class_test(fscanf(infile,"%le",&(psd->PCA_x[i]))!=1,
                      psd->error_message,
                      "Could not read z at line %i in file '%s'",i+headlines+5,ppr->br_exact_file);
    class_test(fscanf(infile,"%le",&(psd->PCA_J_T[i]))!=1,
                      psd->error_message,
                      "Could not read f_g at line %i in file '%s'",i+headlines+5,ppr->br_exact_file);
    class_test(fscanf(infile,"%le",&(psd->PCA_J_y[i]))!=1,
                      psd->error_message,
                      "Could not read f_y at line %i in file '%s'",i+headlines+5,ppr->br_exact_file);
    class_test(fscanf(infile,"%le",&(psd->PCA_J_mu[i]))!=1,
                      psd->error_message,
                      "Could not read f_mu at line %i in file '%s'",i+headlines+5,ppr->br_exact_file);
    for(int j=0; j<psd->index_s; ++j){
      class_test(fscanf(infile,"%le",&(psd->S_vec[i*psd->index_s+j]))!=1,
                        psd->error_message,
                        "Could not read E vector at line %i in file '%s'",i+headlines+5,ppr->br_exact_file);
    }
  }

  return _SUCCESS_;

}


/**
 * Free from distortions_read_PCA_dist_shapes_data()
 *
 * @param psd     Input: pointer to distortions structure (to be freed)
 * @return the error status
 */
int distortions_free_PCA_dist_shapes_data(struct distortions * psd){
    /* Free from distortions_read_BR_exact_data() */
    free(psd->PCA_x);
    free(psd->PCA_J_T);
    free(psd->PCA_J_y);
    free(psd->PCA_J_mu);
    free(psd->S_vec);
}


/**
 * Outputs
 */
int heating_output_titles(char titles[_MAXTITLESTRINGLENGTH_]){
  class_store_columntitle(titles,"z",_TRUE_);
  class_store_columntitle(titles,"Heating function d(Q/rho)/dln(z) [-]",_TRUE_);

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
  }

  return _SUCCESS_;
}

int distortions_output_titles(char titles[_MAXTITLESTRINGLENGTH_]){
  class_store_columntitle(titles,"x",_TRUE_);
  class_store_columntitle(titles,"Spectral distortions DI [10^26 W m^-2 Hz^-1 sr^-1]",_TRUE_);

  return _SUCCESS_;
}
int distortions_output_data(struct distortions * psd,
                            int number_of_titles,
                            double * data){
  int storeidx;
  double * dataptr;

  for (int index_x=0; index_x<psd->x_size; index_x++) {
    dataptr = data + index_x*number_of_titles;
    storeidx = 0;

    class_store_double(dataptr, psd->x[index_x], _TRUE_,storeidx);
    class_store_double(dataptr, psd->DI[index_x]*1.e26, _TRUE_,storeidx);
  }

  return _SUCCESS_;
}



