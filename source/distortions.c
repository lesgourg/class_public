/** @file distortions.c Documented module on spectral distortions
 * Matteo Lucca, 31.10.2018
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
int distortions_init(
		     struct precision * ppr,
                     struct background * pba,
                     struct perturbs * ppt,
                     struct thermo * pth,
                     struct primordial * ppm,
                     struct distortions * psd) {
    /* Definitions */
    double * pvecheat;
    double * pvecdist;
    int last_index = 0;
    double *dg, * dmu, * dy, y_heat, * dr;

    if(psd->has_distortions == _TRUE_){

        /* Assign values to all indices in the distortions structure */
        class_call(distortions_indices(psd),
                   psd->error_message,
                   psd->error_message);

        /* Set details of z and x arrays as used in distortions.c */
        psd->z_min = 1.e3;
        psd->z_max = 5.e6;
        psd->z_size = (int) 550;
        psd->z_delta = (log(psd->z_max)-log(psd->z_min))/psd->z_size;

        psd->x_min = 1.e-2;
        psd->x_max = 5.e1;
        psd->x_size = (int) 550;
        psd->x_delta = (log(psd->x_max)-log(psd->x_min))/psd->x_size;

        /* Set initial values for distortion quantities */
	psd->y = 0.;

        /* Allocate space for heating output parameters */
        class_alloc(psd->z,
                    psd->z_size*sizeof(double),
                    psd->error_message);
        class_alloc(psd->dQrho_dz_tot,
                    psd->z_size*sizeof(double),
                    psd->error_message);
        /* Allocate space for heating private parameters */
        class_alloc(dg,
                    psd->z_size*sizeof(double),
                    psd->error_message);
        class_alloc(dmu,
                    psd->z_size*sizeof(double),
                    psd->error_message);
        class_alloc(dy,
                    psd->z_size*sizeof(double),
                    psd->error_message);
        class_alloc(dr,
                    psd->z_size*sizeof(double),
                    psd->error_message);
        /* Allocate space for distortions output parameters */
        class_alloc(psd->x,
                    psd->x_size*sizeof(double),
                    psd->error_message);
        class_alloc(psd->DI,
                    psd->x_size*sizeof(double),
                    psd->error_message);

        class_call(read_Greens_data(ppr,psd),
                   psd->error_message,
                   psd->error_message);

        psd->distortions_verbose = 1;

        /** Store quantities from heating_at_z() */ 
        if (psd->distortions_verbose > 0) { printf("Computing heating \n"); }

        for (int index_z = 0; index_z<psd->z_size; index_z++) {
            psd->z[index_z] = exp(log(psd->z_min)+psd->z_delta*index_z);

            class_alloc(pvecheat,
                        psd->ht_size*sizeof(double),
                        psd->error_message);
            class_call(heating_at_z(ppr,pba,ppt,pth,ppm,psd,
                                    psd->z[index_z],
                                    pvecheat),
                       psd->error_message,
                       psd->error_message);

            /* Public quantities */
            psd->dQrho_dz_tot[index_z] = pvecheat[psd->index_ht_dQrho_dz_tot]*(1.+psd->z[index_z]);

            /* Private quantities */
            dg[index_z] = psd->dQrho_dz_tot[index_z]*pvecheat[psd->index_ht_f_g]/4.;
            dmu[index_z] = psd->dQrho_dz_tot[index_z]*pvecheat[psd->index_ht_f_mu]*1.401;
            dy[index_z] = psd->dQrho_dz_tot[index_z]*pvecheat[psd->index_ht_f_y]/4.;
            dr[index_z] = psd->dQrho_dz_tot[index_z]*pvecheat[psd->index_ht_f_r];

            /* Free space */
            free(pvecheat);
        }

        /* Calculate g, mu, y and r parameter from heating function */
        class_call(simpson_integration(psd->z_size,
                                       dg,
                                       psd->z_delta,
                                       &psd->g,
                                       psd->error_message),
                   psd->error_message,
                   psd->error_message);
        free(dg);
        class_call(simpson_integration(psd->z_size,
                                       dmu,
                                       psd->z_delta,
                                       &psd->mu,
                                       psd->error_message),
                   psd->error_message,
                   psd->error_message);
        free(dmu);
        class_call(simpson_integration(psd->z_size,
                                       dy,
                                       psd->z_delta,
                                       &y_heat,
                                       psd->error_message),
                   psd->error_message,
                   psd->error_message);
	psd->y += y_heat;
        free(dy);
        class_call(simpson_integration(psd->z_size,
                                       dr,
                                       psd->z_delta,
                                       &psd->r,
                                       psd->error_message),
                   psd->error_message,
                   psd->error_message);
        free(dr);

 	/* 
         * Include additional sources of distortions:
         *    1) Superposition of blackbodies caused by CMB dipole as described in Chluba & Sunyaev 2004 
         *    2) Reionization and structure formation as described in Hill et al. 2015 
         *        a) without relativistic effects and
         *        b) with relativistic effects (TODO)
         * (see also Chluba 2016 for useful discussion)
         */

        /** 1) Superposition of blackbodies */
        psd->y += 2.525e-7;   // Dipole
        psd->y += 4.59e-13;   // Quadrupole

        /** 2a) Reionization and structure formation without relativistic effects */
        psd->y += 1.77e-6;

        /** 2b) Contribution from relativistic effects */

        /** Calculate total heating */
        psd->Drho_over_rho = 4.*psd->g+psd->mu/1.401+4.*psd->y+psd->r;

        /* 
         * Set terminal output
         */
        if (psd->distortions_verbose > 0) { printf("-> total injected/extracted heat = %g\n", psd->Drho_over_rho); }

        if (psd->distortions_verbose > 0) {
            if (psd->mu > 9.e-5) { printf("-> mu-parameter = %g. WARNING: The value of your mu-parameter is larger than the FIRAS constraint mu<9e-5.\n", psd->mu); }
            else{ printf("-> mu-parameter = %g\n", psd->mu); }
        }

        if (psd->distortions_verbose > 0) {
            if (psd->y>1.5e-5) { printf("-> y-parameter = %g. WARNING: The value of your y-parameter is larger than the FIRAS constraint y<1.5e-5.\n", psd->y); }
            else{ printf("-> y-parameter = %g\n", psd->y); }
        }

        if (psd->distortions_verbose > 0) { printf("-> r-parameter = %g\n", psd->r); }

        /** Store quantities from distortions_at_x() */ 
        if (psd->distortions_verbose > 0) { printf("Computing distortions \n"); }

        for (int index_x = 0; index_x<psd->x_size; index_x++) {
            psd->x[index_x] = exp(log(psd->x_min)+psd->x_delta*index_x);

            class_alloc(pvecdist,
                        psd->sd_size*sizeof(double),
                        psd->error_message);
            class_call(distortions_at_x(pba,ppt,pth,ppm,psd,
                                        psd->x[index_x],
                                        pvecdist),
                       psd->error_message,
                       psd->error_message);

            psd->DI[index_x] = pvecdist[psd->index_sd_DI];

            /* Free space */
            free(pvecdist);
        }
    }

    return _SUCCESS_;
}


/**
 * Free all memory space allocated by distortions_init()
 *
 * @param psd     Input: pointer to distortions structure (to be freed)
 * @return the error status
 */
int distortions_free(
                     struct distortions * psd) {
    if(psd->has_distortions == _TRUE_){
        free(psd->z);
        free(psd->dQrho_dz_tot);

        free(psd->x);
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
int distortions_indices(
                        struct distortions * psd) {
    /* Definitions */
    int index_ht = 0;
    int index_sd = 0;

    /* z-dependent parameters */
    class_define_index(psd->index_ht_dQrho_dz_cool,_TRUE_,index_ht,1);
    class_define_index(psd->index_ht_dQrho_dz_diss,_TRUE_,index_ht,1);
    class_define_index(psd->index_ht_dQrho_dz_ann,_TRUE_,index_ht,1);
    class_define_index(psd->index_ht_dQrho_dz_dec,_TRUE_,index_ht,1);
    class_define_index(psd->index_ht_dQrho_dz_eva_PBH,_TRUE_,index_ht,1);
    class_define_index(psd->index_ht_dQrho_dz_acc_PBH,_TRUE_,index_ht,1);
    class_define_index(psd->index_ht_dQrho_dz_tot,_TRUE_,index_ht,1);

    class_define_index(psd->index_ht_f,_TRUE_,index_ht,1);
    class_define_index(psd->index_ht_f_g,_TRUE_,index_ht,1);
    class_define_index(psd->index_ht_f_mu,_TRUE_,index_ht,1);
    class_define_index(psd->index_ht_f_y,_TRUE_,index_ht,1);
    class_define_index(psd->index_ht_f_r,_TRUE_,index_ht,1);

    psd->ht_size = index_ht;

    /* x-dependent parameters */
    class_define_index(psd->index_sd_Y,_TRUE_,index_sd,1);
    class_define_index(psd->index_sd_M,_TRUE_,index_sd,1);
    class_define_index(psd->index_sd_G,_TRUE_,index_sd,1);
    class_define_index(psd->index_sd_DI,_TRUE_,index_sd,1);

    psd->sd_size = index_sd;

    return _SUCCESS_;
}


/**
 * Calculate all redshift dependent quantities needed to compute the spectral distortions, i.e. 
 * the branching ratios of the Green's functions and the heating rates. 
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
int heating_at_z(
		 struct precision * ppr,
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
    double z_dec, p_dec, tau_dec, t_dec, * pvecback_dec, G_dec;

    /* Assign values to all indices in the distortions structure */
    class_call(distortions_indices(psd),
               psd->error_message,
               psd->error_message);

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
    O_b = pba->Omega0_b;                                                              // [-]
    O_cdm = pba->Omega0_cdm;                                                          // [-]
    h = pba->h;                                                                       // [-]
    H = pvecback[pba->index_bg_H];                                                    // [1/Mpc]
    a = pvecback[pba->index_bg_a];                                                    // [-]
    t = pvecback[pba->index_bg_time];                                                 // [Mpc]
    t /= _s_over_Mpc_;                                                                // [s]
    rho_g = pvecback[pba->index_bg_rho_g];                                            // [1/Mpc^4]
    rho_g /= _GeVcm3_over_Mpc4_;                                                      // [GeV/cm^3]
    R = (3./4.)*pvecback[pba->index_bg_rho_b]/pvecback[pba->index_bg_rho_g];          // [-]
    T_g0 = pba->T_cmb;                                                                // [K]

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

    dk = pvecthermo[pth->index_th_dkappa];                                            // [1/Mpc]
    dz_kD = (1./(H*dk))*(16.0/15.0+pow(R,2.0)/(1.0+R))/(6.0*(1.0+R));                 // [Mpc^2]
    kD = 2.*_PI_/pvecthermo[pth->index_th_r_d];                                       // [1/Mpc]
    N_e = pth->n_e;                                                                   // [1/m^3] (today)
    X_e = pvecthermo[pth->index_th_xe];                                               // [-]
    Y_He = pth->YHe;                                                                  // [-]
    Y_p = 1.-Y_He;                                                                    // [-]

    /* Free allocated space */ 
    free(pvecback);
    free(pvecthermo);

    /** BRANCHING RATIOS */

    z_muy = 5.e4;
    z_th = 1.98e6*
           pow((1.-Y_He/2.)/0.8767,-2./5.)*
           pow(O_b*pow(h,2.)/0.02225,-2./5.)*
           pow(T_g0/2.726,1./5.);
    pvecheat[psd->index_ht_f] = exp(-pow(z/z_th,2.5));

    /** 1) Calculate branching ratios using sharp_sharp transition */ 
    if(psd->branching_approx == 1){
        if(z>z_th){ 
            pvecheat[psd->index_ht_f_g] = 1.; 
            pvecheat[psd->index_ht_f_mu] = 0.;
            pvecheat[psd->index_ht_f_y] = 0.;
        }
        if(z<z_th && z>z_muy){ 
            pvecheat[psd->index_ht_f_g] = 0.; 
            pvecheat[psd->index_ht_f_mu] = 1.;
            pvecheat[psd->index_ht_f_y] = 0.;
        }
        if(z<z_muy){ 
            pvecheat[psd->index_ht_f_g] = 0.; 
            pvecheat[psd->index_ht_f_mu] = 0.;
            pvecheat[psd->index_ht_f_y] = 1.;
        }
        pvecheat[psd->index_ht_f_r] = 1.-pvecheat[psd->index_ht_f_g]
                                        -pvecheat[psd->index_ht_f_mu]
                                        -pvecheat[psd->index_ht_f_y]; 
    }

    /** 2) Calculate branching ratios using sharp_soft transition */ 
    if(psd->branching_approx == 2){
        pvecheat[psd->index_ht_f_g] = 1.-pvecheat[psd->index_ht_f];
        if(z>z_muy){ 
            pvecheat[psd->index_ht_f_mu] = pvecheat[psd->index_ht_f];
            pvecheat[psd->index_ht_f_y] = 0.;
        }
        if(z<z_muy){ 
            pvecheat[psd->index_ht_f_mu] = 0.;
            pvecheat[psd->index_ht_f_y] = 1.;
        }
        pvecheat[psd->index_ht_f_r] = 1.-pvecheat[psd->index_ht_f_g]
                                        -pvecheat[psd->index_ht_f_mu]
                                        -pvecheat[psd->index_ht_f_y]; 
    }

    /** 3) Calculate branching ratios unsing soft_soft transitions */ 
    if(psd->branching_approx == 3){
        pvecheat[psd->index_ht_f_g] = 1.-pvecheat[psd->index_ht_f];
        pvecheat[psd->index_ht_f_mu] = pvecheat[psd->index_ht_f]*
                                       (1.0-exp(-pow((1.0+z)/(5.8e4),1.88)));
        pvecheat[psd->index_ht_f_y] = 1.0/(1.0+pow((1.0+z)/(6.0e4),2.58));
        pvecheat[psd->index_ht_f_r] = 1.-pvecheat[psd->index_ht_f_g]
                                        -pvecheat[psd->index_ht_f_mu]
                                        -pvecheat[psd->index_ht_f_y]; 
    }

    /** 4) Calculate branching ratios unsing soft_soft_cons transitions */ 
    if(psd->branching_approx == 4){
        pvecheat[psd->index_ht_f_g] = 1.-pvecheat[psd->index_ht_f];
        pvecheat[psd->index_ht_f_y] = 1.0/(1.0+pow((1.0+z)/(6.0e4),2.58));
        pvecheat[psd->index_ht_f_mu] = pvecheat[psd->index_ht_f]*(1.-pvecheat[psd->index_ht_f_y]);
        pvecheat[psd->index_ht_f_r] = 1.-pvecheat[psd->index_ht_f_g]
                                        -pvecheat[psd->index_ht_f_mu]
                                        -pvecheat[psd->index_ht_f_y]; 
    }

    /** 5) Calculate branching ratios according to Chluba & Jeong 2014 */ 
    if(psd->branching_approx == 5){
    }


    /** HEATING FUNCTIONS */

    /** 1) Adiabatically cooling electrons and barions */
    tilde_rho_g = rho_g/(_m_e_/_GeV_over_kg_);                                        // [1/cm^3]
    theta_g = (_k_B_*T_g0*(1.+z))/(_m_e_*pow(_c_,2.));                                // [-]
    alpha_h = (3./2.)*N_e*1.e-6*pow(1.+z,3.)*(1.+Y_He+X_e);                           // [1/cm^3]
    pvecheat[psd->index_ht_dQrho_dz_cool] = -a*alpha_h/tilde_rho_g*theta_g*
                                             pvecheat[psd->index_ht_f];               // [-]

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
                                                dz_kD*exp(-2.*pow(k/kD,2.));          // [-]
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
    pvecheat[psd->index_ht_dQrho_dz_diss] *= pvecheat[psd->index_ht_f];               // [-]

    /** 3) Cosmological recombination */

    /** 4) Annihilating particles */
    g_h = (1.+Y_He+2*(X_e))/(3.*(1+Y_He));                                            // [-] (TODO: INCOMPLETE)

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
                                                (rho_g/(_eV_over_joules_*1.e-9)*1.e6)*
                                                pvecheat[psd->index_ht_f];            // [-]

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
                                                (rho_g/(_eV_over_joules_*1.e-9)*1.e6)*
                                                pvecheat[psd->index_ht_f];            // [-]
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
                                                    (rho_g/(_eV_over_joules_*1.e-9)*1.e6)*
                                                    pvecheat[psd->index_ht_f];        // [-]
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
                                                    (rho_g/(_eV_over_joules_*1.e-9)*1.e6)*
                                                    pvecheat[psd->index_ht_f];        // [-]
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
 * Evaluation of the spectral distortions once the value of the g, mu and y-parameter is knwon.
 *
 * @param pba        Input: pointer to background structure
 * @param ppt        Input: pointer to the perturbations structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param ppm        Input: pointer to the primordial structure
 * @param psd        Input: pointer to the distortions structure
 * @param z          Input: redshift
 * @param pvecdist   Output: vector of distortions functions (assumed to be already allocated)
 */
int distortions_at_x(
                     struct background* pba,
                     struct perturbs * ppt,
                     struct thermo * pth,
                     struct primordial * ppm,
                     struct distortions * psd,
                     double x,
                     double * pvecdist) {
    /* Assign values to all indices in the distortions structure */
    class_call(distortions_indices(psd),
               psd->error_message,
               psd->error_message);

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
 */
int read_Greens_data(
		     struct precision * ppr,
		     struct distortions * psd){
    FILE * fA;
    char line[_LINE_LENGTH_MAX_];
    char * left;
    int array_line=0;

    psd->Greens_lines = 0;
    class_open(fA, ppr->Greens_file, "r", 
               psd->error_message);

    char* input_template;
    while (fgets(line,_LINE_LENGTH_MAX_-1,fA) != NULL) {
        /* Eliminate blank spaces at beginning of line */
        left=line;
        while (left[0]==' ') {
            left++;
        }
        /* check that the line is neither blank nor a comment. In ASCII, left[0]>39 means that first non-blank charachter might
           be the beginning of some data (it is not a newline, a #, a %, etc.) */
        if (left[0] > 39) {
            if (psd->Greens_lines == 0) {
                /* read number of lines, infer size of arrays and allocate them */
                class_test(sscanf(line, "%d", &psd->Greens_lines) != 1,
                           psd->error_message,
                           "could not read value of parameters num_lines in file %s\n", ppr->Greens_file);

                class_alloc(psd->Greens_z, psd->Greens_lines*sizeof(double), psd->error_message);
                class_alloc(psd->Greens_T_ini, psd->Greens_lines*sizeof(double), psd->error_message);
                class_alloc(psd->Greens_T_last, psd->Greens_lines*sizeof(double), psd->error_message);
                class_alloc(psd->Greens_rho, psd->Greens_lines*sizeof(double), psd->error_message);

	        class_alloc(input_template, (3*psd->Greens_lines)*sizeof(char), psd->error_message);
	        for(i=0;i<psd->Greens_lines;++i){strcat(input_template,"%d ");}
            }
            else{

                class_test(sscanf(line, "%d", &psd->Greens_lines) != ,
                           psd->error_message,
                           "could not read value of parameters num_lines in file %s\n", ppr->Greens_file);

                psd->Greens_z = (double *)line;
                psd->Greens_T_ini = (double *)line;
                psd->Greens_T_last = (double *)line;
                psd->Greens_rho = (double *)line;

            }


        }
    }
    free(input_template);
                printf("%s\n",line);

return _SUCCESS_;

}


/**
 * Outputs
 */
int heating_output_titles(
                          char titles[_MAXTITLESTRINGLENGTH_]){
    class_store_columntitle(titles,"z",_TRUE_);
    class_store_columntitle(titles,"Heating function d(Q/rho)/dln(z) [-]",_TRUE_);

    return _SUCCESS_;
}
int heating_output_data(
                        struct distortions * psd,
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

int distortions_output_titles(
                              char titles[_MAXTITLESTRINGLENGTH_]){
        class_store_columntitle(titles,"x",_TRUE_);
        class_store_columntitle(titles,"Spectral distortions DI [10^26 W m^-2 Hz^-1 sr^-1]",_TRUE_);

    return _SUCCESS_;
}
int distortions_output_data(
                            struct distortions * psd,
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



