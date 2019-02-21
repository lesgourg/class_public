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
                     struct perturbs * ppt,
                     struct thermo * pth,
                     struct primordial * ppm,
                     struct distortions * psd) {

  /** Define local variables */
  int last_index = 0;
  int index_br;
  int index_x;

  if(psd->has_distortions == _FALSE_){
    return _SUCCESS_;
  }
  if (psd->distortions_verbose > 0) {
    printf("Computing spectral distortions \n");
  }

  class_test(pth->compute_damping_scale==_FALSE_,
             psd->error_message,
             "Cannot compute spectral distortions without damping scale \n");

  /** Assign values to all indices in the distortions structure */
  class_call(distortions_indices(psd),
             psd->error_message,
             psd->error_message);

  /** Define z and x arrays TODO :: convert to precision parameters the limits and N */
  class_call(distortions_get_xz_lists(pba,pth,psd),
             psd->error_message,
             psd->error_message);

  /** Define branching ratios */
  class_call(distortions_branching_ratios(ppr,psd),
             psd->error_message,
             psd->error_message);

  /** Define heating function */
  class_call(distortions_heating_rate(ppr,pba,ppt,pth,ppm,psd),
             psd->error_message,
             psd->error_message);

  /** Define spectral distortion amplitudes */
  class_call(distortions_amplitudes(psd),
             psd->error_message,
             psd->error_message);

  /** Define final spectral distortions */
  class_call(distortions_spectral_shapes(ppr,pba,psd),
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
  if(psd->has_distortions == _TRUE_){
    free(psd->z);
    free(psd->x);
    free(psd->z_weights);

    free(psd->br_table);
    free(psd->heating_table);
    free(psd->sd_parameter_table);
    free(psd->distortions_table);

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

  /** Define local variables */
  int index_ht = 0;
  int index_type = 0;

  /** Define indeces for table heating_table defined in distortions_heating_rate() */
  class_define_index(psd->index_ht_dQrho_dz_cool,_TRUE_,index_ht,1);
  class_define_index(psd->index_ht_dQrho_dz_diss,_TRUE_,index_ht,1);
  class_define_index(psd->index_ht_dQrho_dz_CRR,_TRUE_,index_ht,1);
  class_define_index(psd->index_ht_dQrho_dz_ann,_TRUE_,index_ht,1);
  class_define_index(psd->index_ht_dQrho_dz_dec,_TRUE_,index_ht,1);
  class_define_index(psd->index_ht_dQrho_dz_eva_PBH,_TRUE_,index_ht,1);
  class_define_index(psd->index_ht_dQrho_dz_acc_PBH,_TRUE_,index_ht,1);
  class_define_index(psd->index_ht_dQrho_dz_tot,_TRUE_,index_ht,1);
  class_define_index(psd->index_ht_dQrho_dlnz_tot_screened,_TRUE_,index_ht,1);

  psd->ht_size = index_ht;

  /** Define indeces for tables - br_table defined in distortions_branching_ratios(),
                                - sd_parameter_table defined in distortions_amplitudes() and
                                - distortions_table defined in distortions_at_x() */
  class_define_index(psd->index_type_g,_TRUE_,index_type,1);
  class_define_index(psd->index_type_y,_TRUE_,index_type,1);
  class_define_index(psd->index_type_mu,_TRUE_,index_type,1);
  class_define_index(psd->index_type_r,_TRUE_,index_type,1);
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
int distortions_get_xz_lists(struct background * pba, 
                             struct thermo * pth, 
                             struct distortions * psd){

  /** Define local variables */
  int index_z, index_x;

  /** Define transition redshifts z_muy and z_th */
  psd->z_muy = 5.e4;
  psd->z_th = 1.98e6*
         pow((1.-pth->YHe/2.)/0.8767,-2./5.)*
         pow(pba->Omega0_b*pow(pba->h,2.)/0.02225,-2./5.)*
         pow(pba->T_cmb/2.726,1./5.);

  /** Define and allocate z array */
  psd->z_min = 1.011e3;
  psd->z_max = 5.e6;
  psd->z_size = (int) 500;
  psd->z_delta = (log(psd->z_max)-log(psd->z_min))/psd->z_size;

  class_alloc(psd->z,
              psd->z_size*sizeof(double),
              psd->error_message);

  for (index_z = 0; index_z < psd->z_size; index_z++) {
    psd->z[index_z] = exp(log(psd->z_min)+psd->z_delta*index_z);
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
  psd->x_min = 1.e-2;
  psd->x_max = 5.e1;
  psd->x_size = (int) 500;
  psd->x_delta = (log(psd->x_max)-log(psd->x_min))/psd->x_size;

  class_alloc(psd->x,
              psd->x_size*sizeof(double),
              psd->error_message);

  for (index_x = 0; index_x<psd->x_size; index_x++) {
    psd->x[index_x] = exp(log(psd->x_min)+psd->x_delta*index_x);
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
 *    4) Use a soft transition at a_mu-y and z_th imposing conservation of energy 
 *       ('branching approx'=soft_soft_cons)
 *    5) Use a PCA method as described in Chluba & Jeong 2014 ('branching approx'=exact)
 *
 * All quantities are stored in the table br_table
 *
 * @param ppr        Input: pointer to the precision structure
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */
int distortions_branching_ratios(struct precision * ppr,
                                 struct distortions* psd){

  /** Define local variables */
  int index_z,index_type,index_e;
  double f_g, f_y, f_mu;
  double * f_E;
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
      psd->br_table[psd->index_type_r][index_z] = 1.-f_g-f_y-f_mu;

    }
  }
  else{
    /* 5) Calculate branching ratios according to Chluba & Jeong 2014
          In this case, read and interpolate precomputed functions (also the multipole expansion of the residual vectors E)
          from external file branching_ratios_exact.dat. The computation has been performed by J. Chluba according
          to Chluba & Jeong 2014. */ 

    /* Read and spline data from file branching_ratios_exact.dat */
    class_call(distortions_read_BR_exact_data(ppr,psd),
               psd->error_message,
               psd->error_message);
    class_call(distortions_spline_BR_exact_data(psd),
               psd->error_message,
               psd->error_message);

    /* Allocate local variable */
    class_alloc(f_E,
                psd->N_PCA*sizeof(double),
                psd->error_message);

    /* Interpolate over z */
    for(index_z=0; index_z<psd->z_size; ++index_z){
      class_call(distortions_interpolate_BR_exact_data(psd,
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
      psd->br_table[psd->index_type_r][index_z] = 1.-f_g-f_y-f_mu;
      for(index_e=0; index_e<psd->N_PCA; ++index_e){
        psd->br_table[psd->index_type_PCA+index_e][index_z] = f_E[index_e];
      }

    }

    /* Free space allocated in distortions_read_BR_exact_data() */
    class_call(distortions_free_BR_exact_data(psd),
               psd->error_message,
               psd->error_message);
    free(f_E);
  }

  return _SUCCESS_;

}


/**
 * Calculate heating rates.
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
 * All quantities are stored in the table heating_table
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param ppt        Input: pointer to the perturbations structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param ppm        Input: pointer to the primordial structure
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */
int distortions_heating_rate(struct precision * ppr,
                             struct background* pba,
                             struct perturbs * ppt,
                             struct thermo * pth,
                             struct primordial * ppm,
                             struct distortions * psd){

  /** Define local variables */
  int index_z;
  double tau;
  int last_index;
  double * pvecback, O_b, h, H, a, t, rho_g, R, T_g0;
  double * pvecthermo, dk, dz_kD, kD, N_e, X_e, Y_He;
  int index_ht;
  double alpha_h, tilde_rho_g, theta_g;
  double k_max, k_min, k_size, k_delta, k, pk_primordial_k, * int_dQrho_dz_diss_full, * int_dQrho_dz_diss_approx;
  double g_h;

  /** Allocate space for heating rates in heating_table */ 
  class_alloc(psd->heating_table,
              psd->ht_size*sizeof(double*),
              psd->error_message);
  for(index_ht=0; index_ht<psd->ht_size; ++index_ht){
    class_alloc(psd->heating_table[index_ht],
                psd->z_size*sizeof(double),
                psd->error_message);
  }

  for(index_z=0; index_z<psd->z_size; ++index_z){
    /* From z to tau */
    class_call(background_tau_of_z(pba,
                                   psd->z[index_z],
                                   &tau),
               pba->error_message,
               psd->error_message);

    /** Import quantities from background */
    last_index = 0;
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
    O_b = pba->Omega0_b;                                                              // [-]
    h = pba->h;                                                                       // [-]
    H = pvecback[pba->index_bg_H];                                                    // [1/Mpc]
    a = pvecback[pba->index_bg_a];                                                    // [-]
    t = pvecback[pba->index_bg_time];                                                 // [Mpc]
    t /= _s_over_Mpc_;                                                                // [s]
    rho_g = pvecback[pba->index_bg_rho_g];                                            // [1/Mpc^4]
    rho_g /= _GeVcm3_over_Mpc4_;                                                      // [GeV/cm^3]
    R = (3./4.)*pvecback[pba->index_bg_rho_b]/pvecback[pba->index_bg_rho_g];          // [-]
    T_g0 = pba->T_cmb;                                                                // [K]

    /** Import quantities from thermodynamics */
    class_alloc(pvecthermo,
                pth->tt_size*sizeof(double),
                psd->error_message);
    class_call(thermodynamics_at_z(pba,
                                   pth,
                                   psd->z[index_z],
                                   pth->inter_normal,
                                   &last_index,
                                   pvecback,
                                   pvecthermo),
               pth->error_message,
               psd->error_message);
    dk = pvecthermo[pth->index_th_dkappa];                                            // [1/Mpc]
    dz_kD = (1./(H*dk))*(16.0/15.0+pow(R,2.0)/(1.0+R))/(6.0*(1.0+R));                 // [Mpc^2]
    kD = 2.*_PI_/pvecthermo[pth->index_th_r_d];                                       // [1/Mpc]
    N_e = pth->n_e;                                                                   // [1/m^3] (today)
    X_e = pvecthermo[pth->index_th_xe];                                               // [-]
    Y_He = pth->YHe;                                                                  // [-]

    /* Free allocated space */
    free(pvecback);
    free(pvecthermo);

    /** Calculate heating rates */

    /* 1) Adiabatically cooling electrons and barions */
    tilde_rho_g = rho_g/(_m_e_/_GeV_over_kg_);                                        // [1/cm^3]
    theta_g = (_k_B_*T_g0*(1.+psd->z[index_z]))/(_m_e_*pow(_c_,2.));                  // [-]
    alpha_h = (3./2.)*N_e*1.e-6*pow(1.+psd->z[index_z],3.)*(1.+Y_He+X_e);             // [1/cm^3]
    psd->heating_table[psd->index_ht_dQrho_dz_cool][index_z] = -a*alpha_h/
                                                               tilde_rho_g*theta_g;   // [-]

    /* 2) dissipation of acoustic waves */
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
                                            dz_kD*exp(-2.*pow(k/kD,2.));              // [-]
      }

      /* Integrate approximate function */
      class_call(simpson_integration(k_size,
                                     int_dQrho_dz_diss_approx,
                                     (log(k_max)-log(k_min))/(k_size),
                                     &psd->heating_table[psd->index_ht_dQrho_dz_diss][index_z],
                                     psd->error_message),
                 psd->error_message,
                 psd->error_message);

      /* Free space */
      free(int_dQrho_dz_diss_approx);
    }

    /* 3) Cosmological recombination */
    psd->heating_table[psd->index_ht_dQrho_dz_CRR][index_z] = 0.;

    /* 4) Annihilating particles */
    g_h = (1.+Y_He+2*(X_e))/(3.*(1+Y_He));                                            // [-] (TODO: INCOMPLETE)
    psd->heating_table[psd->index_ht_dQrho_dz_ann][index_z] = 0.;

    /* 5) Decaying relic particles */
    psd->heating_table[psd->index_ht_dQrho_dz_dec][index_z] = 0.;

    /* 6) Evaporation of primordial black holes */
    psd->heating_table[psd->index_ht_dQrho_dz_eva_PBH][index_z] = 0.;

    /* 7) Accretion of matter into primordial black holes */
    psd->heating_table[psd->index_ht_dQrho_dz_acc_PBH][index_z] = 0.;

    /* Total heating rate */
    psd->heating_table[psd->index_ht_dQrho_dz_tot][index_z] = 
                                        //psd->heating_table[psd->index_ht_dQrho_dz_cool][index_z] +
                                        psd->heating_table[psd->index_ht_dQrho_dz_diss][index_z] +
                                        psd->heating_table[psd->index_ht_dQrho_dz_CRR][index_z] +
                                        psd->heating_table[psd->index_ht_dQrho_dz_ann][index_z] +
                                        psd->heating_table[psd->index_ht_dQrho_dz_dec][index_z] +
                                        psd->heating_table[psd->index_ht_dQrho_dz_eva_PBH][index_z] +
                                        psd->heating_table[psd->index_ht_dQrho_dz_acc_PBH][index_z];
    psd->heating_table[psd->index_ht_dQrho_dlnz_tot_screened][index_z] = 
                                        psd->heating_table[psd->index_ht_dQrho_dz_tot][index_z]*
                                        exp(-pow(psd->z[index_z]/psd->z_th,2.5))*
                                        (1.+psd->z[index_z]);
  }

  return _SUCCESS_;

}


/**
 * Calculate spectral distortions amplitudes.
 *
 * All quantities are stored in the table heating_table
 *
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */
int distortions_amplitudes(struct distortions * psd){

  /** Define local variables */
  int index_type, index_z, index_e;
  double * integrand;

  /** Allocate space for spectral distortion amplitude in table sd_parameter_table */
  class_alloc(psd->sd_parameter_table,
              psd->type_size*sizeof(double),
              psd->error_message);

  /** Compute distortion amplitude corresponding to each branching ratio (g, y, mu and r) */
  class_alloc(integrand,
              psd->z_size*sizeof(double),
              psd->error_message);

  for(index_type=0; index_type<psd->type_size; ++index_type){
    for (index_z=0; index_z<psd->z_size; ++index_z){
      integrand[index_z] = psd->heating_table[psd->index_ht_dQrho_dlnz_tot_screened][index_z]*
                           psd->br_table[index_type][index_z];
    }
    class_call(simpson_integration(psd->z_size,
                                   integrand,
                                   psd->z_delta,
                                   &psd->sd_parameter_table[index_type],
                                   psd->error_message),
               psd->error_message,
               psd->error_message);
  }

  free(integrand);

  psd->sd_parameter_table[psd->index_type_g] /= 4.;
  psd->sd_parameter_table[psd->index_type_y] /= 4.;
  psd->sd_parameter_table[psd->index_type_mu] *= 1.401;

  /** Include additional sources of distortions (see also Chluba 2016 for useful discussion) */
  //psd->sd_parameter_table[psd->index_type_y] += 2.525e-7;   // CMB Dipole (Chluba & Sunyaev 2004)
  //psd->sd_parameter_table[psd->index_type_y] += 4.59e-13;   // CMB Quadrupole (Chluba & Sunyaev 2004)
  //psd->sd_parameter_table[psd->index_type_y] += 1.77e-6;    // Reionization and structure formation (Hill et al. 2015)

  /** Calculate total heating */
  psd->Drho_over_rho = psd->sd_parameter_table[psd->index_type_g]*4.+
                       psd->sd_parameter_table[psd->index_type_y]*4.+
                       psd->sd_parameter_table[psd->index_type_mu]/1.401+
                       psd->sd_parameter_table[psd->index_type_r];

  if(psd->N_PCA != 0){
    /** Calculate mu_k for PCA expansion */
    psd->sd_parameter_table[psd->index_type_PCA+index_e] = 0.;
    for(index_e=0; index_e<psd->N_PCA; ++index_e){
      for(index_z=0; index_z<psd->z_size; ++index_z){
        psd->sd_parameter_table[psd->index_type_PCA+index_e] += psd->br_table[psd->index_type_PCA+index_e][index_z]*
                                                                psd->heating_table[psd->index_type_PCA+index_e][index_z];
      }
    }

  }

  psd->distortions_verbose = 2;

  /* Print found parameters */
  if (psd->distortions_verbose > 1){
    printf("-> total injected/extracted heat = %g\n", psd->Drho_over_rho);

    if (psd->sd_parameter_table[psd->index_type_mu] > 9.e-5) {
      printf("-> mu-parameter = %g. WARNING: The value of your mu-parameter is larger than the FIRAS constraint mu<9e-5.\n", psd->sd_parameter_table[psd->index_type_mu]);
    }
    else{ 
      printf("-> mu-parameter = %g\n", psd->sd_parameter_table[psd->index_type_mu]);
      printf("Chluba 2016 (diss, exact): mu-parameter = %g\n",2.00e-08);
    }

    if (psd->sd_parameter_table[psd->index_type_y]>1.5e-5) {
      printf("-> y-parameter = %g. WARNING: The value of your y-parameter is larger than the FIRAS constraint y<1.5e-5.\n", psd->sd_parameter_table[psd->index_type_y]);
    }
    else{ 
      printf("-> y-parameter = %g\n", psd->sd_parameter_table[psd->index_type_y]);
      printf("Chluba 2016 (diss, exact): y-parameter = %g\n",3.63e-9);
    }

    if(psd->N_PCA == 0){
      printf("-> r-parameter = %g\n", psd->sd_parameter_table[psd->index_type_r]);
    }
    else{
       for(index_e=0; index_e<psd->N_PCA; ++index_e){
         printf("-> PCA multipole mu_%d = %g\n", index_e+1, psd->sd_parameter_table[psd->index_type_PCA+index_e]);
       }
       printf("Chluba 2016 (diss, exact): mu_1 = %g\n",3.81e-08);
       printf("Chluba 2016 (diss, exact): mu_2 = %g\n",-1.19e-09);
    }
  }

  return _SUCCESS_;

}


/**
 * Calculate spectral distortions.
 *
 * @param pba        Input: pointer to background structure
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */
int distortions_spectral_shapes(struct precision * ppr,
                                struct background * pba,
                                struct distortions * psd){

  /** Define local variables */
  double * S;
  int last_index = 0;
  int index_type, index_x, index_s;

  /** Allocate space for distortions shapes in distortions_table */ 
  class_alloc(psd->distortions_table,
              psd->type_size*sizeof(double*),
              psd->error_message);
  for(index_type=0; index_type<psd->type_size; ++index_type){
    class_alloc(psd->distortions_table[index_type],
                psd->x_size*sizeof(double),
                psd->error_message);
  }

  /** Allocate space for final spactral distortion */ 
  class_alloc(psd->DI,
              psd->x_size*sizeof(double),
              psd->error_message);

  if(psd->N_PCA > 0){
    /* Read and spline data from file branching_ratios_exact.dat */
    class_call(distortions_read_PCA_dist_shapes_data(ppr,psd),
               psd->error_message,
               psd->error_message);
    class_call(distortions_spline_PCA_dist_shapes_data(psd),
               psd->error_message,
               psd->error_message);
    /** Calculate spectral distortions */
    for(index_x=0; index_x<psd->x_size; ++index_x){
      /* Allocate loca variable */
      class_alloc(S,
                  psd->N_PCA*sizeof(double),
                  psd->error_message);
      /* Interpolate over z */
      class_call(distortions_interpolate_PCA_dist_shapes_data(psd,
                                                              psd->x[index_x]*(_k_B_*pba->T_cmb/_h_P_),
                                                              &psd->distortions_table[psd->index_type_g][index_x],
                                                              &psd->distortions_table[psd->index_type_y][index_x],
                                                              &psd->distortions_table[psd->index_type_mu][index_x],
                                                              S,
                                                              &last_index),
                 psd->error_message,
                 psd->error_message);

      for(index_s=0; index_s<psd->N_PCA; ++index_s){
        psd->distortions_table[psd->index_type_PCA+index_s][index_x] = S[index_s];
      }

      /* Free allocated space */
      class_call(distortions_free_PCA_dist_shapes_data(psd),
                 psd->error_message,
                 psd->error_message);
      free(S);
    }
  }
  else{

    /** Calculate spectral distortions */
    for(index_x=0; index_x<psd->x_size; ++index_x){
      psd->distortions_table[psd->index_type_g][index_x] = 1.e-18*2.*pow(_k_B_*pba->T_cmb,3.)/pow(_h_P_*_c_,2.)*
                                                           pow(psd->x[index_x],4.)*exp(-psd->x[index_x])/
                                                           pow(1.-exp(-psd->x[index_x]),2.);   // [10^-18 W/(m^2 Hz sr)]
      psd->distortions_table[psd->index_type_y][index_x] = psd->distortions_table[psd->index_type_g][index_x]*
                                                           (psd->x[index_x]*(1.+exp(-psd->x[index_x]))/
                                                           (1.-exp(-psd->x[index_x]))-4.);     // [10^-18 W/(m^2 Hz sr)]
      psd->distortions_table[psd->index_type_mu][index_x] = psd->distortions_table[psd->index_type_g][index_x]*
                                                            (1./2.19229-1./psd->x[index_x]);   // [10^-18 W/(m^2 Hz sr)]
    }
  }

  /** Calculate spectral distortions */
  for(index_x=0; index_x<psd->x_size; ++index_x){
    psd->DI[index_x] = psd->sd_parameter_table[psd->index_type_y]*psd->distortions_table[psd->index_type_g][index_x]+
                       psd->sd_parameter_table[psd->index_type_mu]*psd->distortions_table[psd->index_type_y][index_x]+
                       psd->sd_parameter_table[psd->index_type_g]*psd->distortions_table[psd->index_type_mu][index_x];  // [10^-18 W/(m^2 Hz sr)]
  }

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
  /** Define local variables */
  FILE * infile;
  char line[_LINE_LENGTH_MAX_];
  char * left;
  int headlines = 0;

  /** Open file */
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
      /** Read number of lines, infer size of arrays and allocate them */
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
  /** Read parameters */
  for(int i=0;i<psd->Greens_Nz;++i){
    class_test(fscanf(infile,"%le",&(psd->Greens_z[i]))!=1,                           // [-]
                      psd->error_message,
                      "Could not read z values at line %i in file '%s'",headlines+1,ppr->Greens_file);
  }
  for(int i=0;i<psd->Greens_Nz;++i){
    class_test(fscanf(infile,"%le",&(psd->Greens_T_ini[i]))!=1,                       // [K]
                      psd->error_message,
                      "Could not read T_ini values at line %i in file '%s'",headlines+2,ppr->Greens_file);
  }
  for(int i=0;i<psd->Greens_Nz;++i){
    class_test(fscanf(infile,"%le",&(psd->Greens_T_last[i]))!=1,                      // [K]
                      psd->error_message,
                      "Could not read T_last values at line %i in file '%s'",headlines+3,ppr->Greens_file);
  }
  for(int i=0;i<psd->Greens_Nz;++i){
    class_test(fscanf(infile,"%le",&(psd->Greens_rho[i]))!=1,                         // [??]
                      psd->error_message,
                      "Could not read rho values at line %i in file '%s'",headlines+4,ppr->Greens_file);
  }
  for(int i=0;i<psd->Greens_Nx;++i){
    class_test(fscanf(infile,"%le",&(psd->Greens_x[i]))!=1,                            // [-]
                      psd->error_message,
                      "Could not read Greens x at line %i in file '%s'",i+headlines+5,ppr->Greens_file);
    for(int j=0;j<psd->Greens_Nz;++j){
      class_test(fscanf(infile,"%le",&(psd->Greens_function[i*psd->Greens_Nz+j]))!=1,  // [??]
                        psd->error_message,
                        "Could not read Greens function at line %i in file '%s'",i+headlines+5,ppr->Greens_file);
    }
    class_test(fscanf(infile,"%le",&(psd->Greens_blackbody[i]))!=1,                    // [??]
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

  return _SUCCESS_;

}


/**
 * Reads the external file branching_ratios_exact calculated according to Chluba & Jeong 2014
 * for PIXIE like configuration (frequency in interval [30,1000] GHz and bin width of 1 GHz).
 * The file has been computed by J. Chluba.
 *
 * @param ppr        Input: pointer to precision structure
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */
int distortions_read_BR_exact_data(struct precision * ppr,
                                   struct distortions * psd){

  /** Define local variables */
  int index_e,index_z;
  FILE * infile;
  char line[_LINE_LENGTH_MAX_];
  char * left;
  int headlines = 0;

  /** Open file */
  class_open(infile, ppr->br_exact_file, "r",
             psd->error_message);

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
                 "could not header (number of lines, number of columns, number of multipoles) at line %i in file '%s' \n",headlines,ppr->br_exact_file);

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
                      "Could not read z at line %i in file '%s'",index_z+headlines,ppr->br_exact_file);
    class_test(fscanf(infile, "%le",
                      &(psd->f_g_exact[index_z]))!=1,                                 // [-]
                      psd->error_message,
                      "Could not read f_g at line %i in file '%s'",index_z+headlines,ppr->br_exact_file);
    class_test(fscanf(infile, "%le",
                      &(psd->f_y_exact[index_z]))!=1,                                 // [-]
                      psd->error_message,
                      "Could not read f_y at line %i in file '%s'",index_z+headlines,ppr->br_exact_file);
    class_test(fscanf(infile,"%le",
                      &(psd->f_mu_exact[index_z]))!=1,                                // [-]
                      psd->error_message,
                      "Could not read f_mu at line %i in file '%s'",index_z+headlines,ppr->br_exact_file);
    for(index_e=0; index_e<psd->E_vec_size; ++index_e){
      class_test(fscanf(infile,"%le",
                        &(psd->E_vec[index_z*psd->E_vec_size+index_e]))!=1,           // [-]
                        psd->error_message,
                        "Could not read E vector at line %i in file '%s'",index_z+headlines,ppr->br_exact_file);
    }
  }

  return _SUCCESS_;

}


/**
 * Spline the quantitites read in distortions_read_BR_exact_data()
 *
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */
int distortions_spline_BR_exact_data(struct distortions* psd){

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
 * Interpolate the quantitites splined in distortions_spline_BR_exact_data()
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
int distortions_interpolate_BR_exact_data(struct distortions* psd,
                                          double z,
                                          double * f_g,
                                          double * f_y,
                                          double * f_mu,
                                          double * f_E,
                                          int * last_index){

  /** Define local variables */
  int index = *last_index;
  int index_e;
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
  for(index_e=0;index_e<psd->N_PCA;++index_e){
    f_E[index_e] = array_interpolate_spline_hunt(psd->E_vec,
                                                 psd->ddE_vec,
                                                 index*psd->E_vec_size+index_e,
                                                 (index+1)*psd->E_vec_size+index_e,
                                                 h,a,b);
  }

  *last_index = index;

  return _SUCCESS_;

}


/**
 * Free from distortions_read_BR_exact_data() and distortions_spline_BR_exact_data()
 *
 * @param psd     Input: pointer to distortions structure (to be freed)
 * @return the error status
 */
int distortions_free_BR_exact_data(struct distortions * psd){
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
 * Reads the external file PCA_distortions_shapes calculated according to Chluba & Jeong 2014
 * for PIXIE like configuration (frequency in interval [30,1000] GHz and bin width of 1 GHz).
 * The file has been computed by J. Chluba.
 *
 * @param ppr        Input: pointer to precision structure
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */
int distortions_read_PCA_dist_shapes_data(struct precision * ppr,
                                          struct distortions * psd){
  /** Define local variables */
  FILE * infile;
  char line[_LINE_LENGTH_MAX_];
  char * left;
  int headlines = 0;

  /** Open file */
  class_open(infile, ppr->br_exact_file, "r",
             psd->error_message);

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
                 "could not header (number of lines, number of columns, number of multipoles) at line %i in file '%s' \n",headlines,ppr->br_exact_file);

      class_alloc(psd->PCA_nu, psd->PCA_Nnu*sizeof(double), psd->error_message);
      class_alloc(psd->PCA_G_T, psd->PCA_Nnu*sizeof(double), psd->error_message);
      class_alloc(psd->PCA_Y_SZ, psd->PCA_Nnu*sizeof(double), psd->error_message);
      class_alloc(psd->PCA_M_mu, psd->PCA_Nnu*sizeof(double), psd->error_message);

      class_alloc(psd->S_vec, psd->PCA_Nnu*psd->S_vec_size*sizeof(double), psd->error_message);
      break;
    }
  }
  /** Read parameters */
  for(int i=0; i<psd->PCA_Nnu; ++i){
    class_test(fscanf(infile,"%le",
                      &(psd->PCA_nu[i]))!=1,                                          // [GHz]
                      psd->error_message,
                      "Could not read z at line %i in file '%s'",i+headlines,ppr->br_exact_file);
    class_test(fscanf(infile,"%le",
                      &(psd->PCA_G_T[i]))!=1,                                         // [10^-18 W/(m^2 Hz sr)]
                      psd->error_message,
                      "Could not read f_g at line %i in file '%s'",i+headlines,ppr->br_exact_file);
    class_test(fscanf(infile,"%le",
                      &(psd->PCA_Y_SZ[i]))!=1,                                        // [10^-18 W/(m^2 Hz sr)]
                      psd->error_message,
                      "Could not read f_y at line %i in file '%s'",i+headlines,ppr->br_exact_file);
    class_test(fscanf(infile,"%le",
                      &(psd->PCA_M_mu[i]))!=1,                                        // [10^-18 W/(m^2 Hz sr)]
                      psd->error_message,
                      "Could not read f_mu at line %i in file '%s'",i+headlines,ppr->br_exact_file);
    for(int j=0; j<psd->S_vec_size; ++j){
      class_test(fscanf(infile,"%le",
                        &(psd->S_vec[i*psd->S_vec_size+j]))!=1,                       // [10^-18 W/(m^2 Hz sr)]
                        psd->error_message,
                        "Could not read E vector at line %i in file '%s'",i+headlines,ppr->br_exact_file);
    }
  }

  return _SUCCESS_;

}


/**
 * Spline the quantitites read in distortions_read_PCA_dist_shapes_data()
 *
 * @param psd        Input: pointer to the distortions structure
 * @return the error status
 */
int distortions_spline_PCA_dist_shapes_data(struct distortions* psd){
  /** Allocate second derivatievs */
  class_alloc(psd->ddG_T_PCA,
              psd->PCA_Nnu*sizeof(double),
              psd->error_message);
  class_alloc(psd->ddY_SZ_PCA,
              psd->PCA_Nnu*sizeof(double),
              psd->error_message);
  class_alloc(psd->ddM_mu_PCA,
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
                                        psd->ddG_T_PCA,
                                        _SPLINE_EST_DERIV_,
                                        psd->error_message),
             psd->error_message,
             psd->error_message);
  class_call(array_spline_table_columns(psd->PCA_nu,
                                        psd->PCA_Nnu,
                                        psd->PCA_Y_SZ,
                                        1,
                                        psd->ddY_SZ_PCA,
                                        _SPLINE_EST_DERIV_,
                                        psd->error_message),
           psd->error_message,
           psd->error_message);
  class_call(array_spline_table_columns(psd->PCA_nu,
                                        psd->PCA_Nnu,
                                        psd->PCA_M_mu,
                                        1,
                                        psd->ddM_mu_PCA,
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
 * Interpolate the quantitites splined in distortions_spline_PCA_dist_shapes_data()
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
int distortions_interpolate_PCA_dist_shapes_data(struct distortions* psd,
                                                 double nu,
                                                 double * G_T,
                                                 double * Y_SZ,
                                                 double * M_mu,
                                                 double * S,
                                                 int * index){
  /** Define local variables */
  int last_index = *index;
  int index_s;
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
                                       psd->ddG_T_PCA,
                                       last_index,
                                       last_index+1,
                                       h,a,b);
  *Y_SZ = array_interpolate_spline_hunt(psd->PCA_Y_SZ,
                                        psd->ddY_SZ_PCA,
                                        last_index,
                                        last_index+1,
                                        h,a,b);
  *M_mu = array_interpolate_spline_hunt(psd->PCA_M_mu,
                                        psd->ddM_mu_PCA,
                                        last_index,
                                        last_index+1,
                                        h,a,b);
  for(index_s=0; index_s<psd->N_PCA; ++index_s){
    S[index_s] = array_interpolate_spline_hunt(psd->S_vec,
                                               psd->ddS_vec,
                                               last_index*psd->S_vec_size+index_s,
                                               (last_index+1)*psd->S_vec_size+index_s,
                                               h,a,b);
  }

  *index = last_index;

  return _SUCCESS_;

}


/**
 * Free from distortions_read_PCA_dist_shapes_data() and distortions_spline_PCA_dist_shapes_data()
 *
 * @param psd     Input: pointer to distortions structure (to be freed)
 * @return the error status
 */
int distortions_free_PCA_dist_shapes_data(struct distortions * psd){
  free(psd->PCA_nu);
  free(psd->PCA_G_T);
  free(psd->ddG_T_PCA);
  free(psd->PCA_Y_SZ);
  free(psd->ddY_SZ_PCA);
  free(psd->PCA_M_mu);
  free(psd->ddM_mu_PCA);
  free(psd->S_vec);
  free(psd->ddS_vec);

  return _SUCCESS_;
}


/**
 * Outputs
 */
int heating_output_titles(char titles[_MAXTITLESTRINGLENGTH_]){
  class_store_columntitle(titles,"Redshift z",_TRUE_);
  class_store_columntitle(titles,"Heating function f(z)*d(Q/rho)/dln(z) [-]",_TRUE_);

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
    class_store_double(dataptr, psd->heating_table[psd->index_ht_dQrho_dlnz_tot_screened][index_z], _TRUE_, storeidx);
  }

  return _SUCCESS_;
}

int distortions_output_titles(char titles[_MAXTITLESTRINGLENGTH_]){
  class_store_columntitle(titles,"Dimentionless frequency x",_TRUE_);
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



