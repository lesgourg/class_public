/** @file noninjection.c Documented non-exotic energy injection module
 *
 * written by Nils Schoeneberg and Matteo Lucca, 27.02.2019
 *
 * The main goal of this module is to calculate the non-injected energy for the
 * photon evolution equation
 * For more details see the description in the README file.
 */
#include "primordial.h"
#include "noninjection.h"

/**
 * Initialize the noninjection structure.
 *
 * @param ppr   Input: pointer to precision structure
 * @param pba   Input: pointer to background structure
 * @param pth   Input: pointer to thermodynamics structure
 * @param ppt   Input: pointer to perturbation structure
 * @param ppm   Input: pointer to primordial structure
 * @param pni   Input: pointer to noninjection structure
 * @return the error status
 */
int noninjection_init(struct precision* ppr,
                      struct background* pba,
                      struct thermodynamics* pth,
                      struct perturbations* ppt,
                      struct primordial* ppm,
                      struct noninjection* pni){

  /** Summary: */

  /** - Define local variables */
  int index_k,index_z;
  int last_index_back, last_index_thermo, last_index_coarse = 0;
  double *pvecback, *pvecthermo;
  double R, dkappa;
  double dEdt;
  double z_wkb;
  double h,a,b;
  double temp_injection;
  double z_coarse;

  /* z-table */
  pni->z_size = pth->tt_size;
  class_alloc(pni->z_table,
              pni->z_size*sizeof(double),
              pni->error_message);
  memcpy(pni->z_table,
         pth->z_table,
         pni->z_size*sizeof(double));
  pni->last_index_z = 0;

  /* Photon non-injected table */
  class_alloc(pni->photon_dep_table,
              pni->z_size*sizeof(double),
              pni->error_message);
  /* Initialize to zero */
  for(index_z = 0; index_z < pni->z_size; ++index_z){
    pni->photon_dep_table[index_z] = 0.;
  }

  /** - Allocate new storage space that is a bit coarser in z as to reduce computational costs */
  pni->logz_max = log(1.+pni->z_table[pni->z_size-1]);
  pni->z_size_coarse = ppr->noninjection_Nz_log; /* This table is only going to be filled after the initial run anyway */
  class_alloc(pni->z_table_coarse,
              pni->z_size_coarse*sizeof(double),
              pni->error_message);
  for(index_z=0;index_z<pni->z_size_coarse-1;++index_z){
    pni->z_table_coarse[index_z] = exp(index_z*pni->logz_max/(pni->z_size_coarse-1))-1.;
  }
  pni->z_table_coarse[pni->z_size_coarse-1] = pni->z_table[pni->z_size-1];

  class_alloc(pni->noninjection_table,
              pni->z_size_coarse*sizeof(double),
              pni->error_message);
  class_alloc(pni->ddnoninjection_table,
              pni->z_size_coarse*sizeof(double),
              pni->error_message);

  /** - Allocate qunatities from primordial structure */
  pni->k_max = ppr->k_max_acc_diss;
  pni->k_min = ppr->k_min_acc_diss;
  pni->k_size = ppr->noninjection_Nk_acc_diss;
  class_alloc(pni->k,
              pni->k_size*sizeof(double),
              pni->error_message);
  class_alloc(pni->k_weights,
              pni->k_size*sizeof(double),
              pni->error_message);
  class_alloc(pni->pk_primordial_k,
              pni->k_size*sizeof(double),
              pni->error_message);

  /** - Import primordial spectrum */
  for (index_k=0; index_k<pni->k_size; index_k++) {
    pni->k[index_k] = exp(log(pni->k_min)+(log(pni->k_max)-log(pni->k_min))/(pni->k_size)*index_k);

    class_call(primordial_spectrum_at_k(ppm,
                                        ppt->index_md_scalars,
                                        linear,
                                        pni->k[index_k],
                                        &pni->pk_primordial_k[index_k]),
               ppm->error_message,
               pni->error_message);
  }

  /** - Prepare wavenumber integration */
  class_call(array_trapezoidal_weights(pni->k,
                                       pni->k_size,
                                       pni->k_weights,
                                       pni->error_message),
             pni->error_message,
             pni->error_message);

  class_alloc(pni->integrand_approx,
              pni->k_size*sizeof(double),
              pni->error_message);

  /** - Allocate backgorund and thermodynamcis vectors */
  last_index_back = 0;
  last_index_thermo = 0;
  class_alloc(pvecback,
              pba->bg_size*sizeof(double),
              pni->error_message);
  class_alloc(pvecthermo,
              pth->tt_size*sizeof(double),
              pni->error_message);

  /** - Calculate WKB approximation ampltidue factor f_nu at early times */
  z_wkb = ppr->z_wkb_acc_diss;
  class_call(background_at_z(pba,
                             z_wkb,
                             long_info,
                             inter_normal,
                             &last_index_back,
                             pvecback),
             pba->error_message,
             pni->error_message);

  pni->f_nu_wkb = (1.-pvecback[pba->index_bg_rho_g]/(pvecback[pba->index_bg_Omega_r]*pvecback[pba->index_bg_rho_crit]));

  /** - Import quantities from other structures */
  /* Background structure */
  pni->H0 = pba->H0*_c_/_Mpc_over_m_;                                                               // [1/s]
  pni->T_g0 = pba->T_cmb;                                                                           // [K]
  pni->Omega0_b = pba->Omega0_b;                                                                    // [-]
  pni->Omega0_cdm = pba->Omega0_cdm;                                                                // [-]
  pni->rho0_cdm = pba->Omega0_cdm*pow(pni->H0,2)*3/8./_PI_/_G_*_c_*_c_;                             // [J/m^3]

  /* Thermodynamics structure */
  pni->fHe = pth->fHe;                                                                              // [-]
  pni->N_e0 = pth->n_e;                                                                             // [1/m^3]

  /** - Loop over z and calculate the heating at each point */
  dEdt = 0.;
  for(index_z=0; index_z<pni->z_size_coarse; ++index_z){

    z_coarse = pni->z_table_coarse[index_z];
    pni->noninjection_table[index_z] = 0.;

    /* Import quantities from background and thermodynamics structure */
    class_call(background_at_z(pba,
                               z_coarse,
                               long_info,
                               inter_closeby,
                               &last_index_back,
                               pvecback),
               pba->error_message,
               pni->error_message);

    pni->H = pvecback[pba->index_bg_H]*_c_/_Mpc_over_m_;                                            // [1/s]
    pni->a = pvecback[pba->index_bg_a];                                                             // [-]
    pni->rho_g = pvecback[pba->index_bg_rho_g]*_Jm3_over_Mpc2_;                                     // [J/m^3]
    R = (3./4.)*pvecback[pba->index_bg_rho_b]/pvecback[pba->index_bg_rho_g];                        // [-]
    pni->nH = pni->N_e0*pow(pni->a,-3);

    class_call(thermodynamics_at_z(pba,
                                   pth,
                                   z_coarse,
                                   inter_normal,
                                   &last_index_thermo,
                                   pvecback,
                                   pvecthermo),
               pth->error_message,
               pni->error_message);

    dkappa = pvecthermo[pth->index_th_dkappa];                                                      // [1/Mpc]
    pni->dkD_dz = 1./(pvecback[pba->index_bg_H]*dkappa)*(16./15.+pow(R,2.)/(1.+R))/(6.*(1.0+R));    // [Mpc^2]
    pni->kD = 2.*_PI_/pvecthermo[pth->index_th_r_d];                                                // [1/Mpc]
    pni->T_b = pvecthermo[pth->index_th_Tb];                                                        // [K]
    pni->T_g = pni->T_g0/pni->a;                                                                    // [K]
    pni->x_e = pvecthermo[pth->index_th_xe];                                                        // [-]
    pni->heat_capacity = (3./2.)*_k_B_*pni->nH*(1.+pni->fHe+pni->x_e);                              // [J/(K m^3)]

    /* Include all non-injected energy that does not need to be deposited (i.e. adiabatic terms as below) */
    /* First order cooling of photons due to adiabatic interaction with baryons */
    class_call(noninjection_rate_adiabatic_cooling(pni,
                                                   z_coarse,
                                                   &dEdt),
               pni->error_message,
               pni->error_message);
    pni->noninjection_table[index_z]+=dEdt;

    /* Second order acoustic dissipation of BAO */
    class_call(noninjection_rate_acoustic_diss(pni,
                                               z_coarse,
                                               &dEdt),
               pni->error_message,
               pni->error_message);
    pni->noninjection_table[index_z]+=dEdt;
  }

  /** - Spline coarse z table in view of interpolation */
  class_call(array_spline_table_columns2(pni->z_table_coarse,
                                         pni->z_size_coarse,
                                         pni->noninjection_table,
                                         1,
                                         pni->ddnoninjection_table,
                                         _SPLINE_EST_DERIV_,
                                         pni->error_message),
             pni->error_message,
             pni->error_message);

  /** - Evaluate the spline coarse table at the fine-grained z samples */
  for(index_z=0;index_z<pni->z_size;++index_z){
    class_call(array_spline_hunt(pni->z_table_coarse,
                                 pni->z_size_coarse,
                                 pni->z_table[index_z],
                                 &last_index_coarse,
                                 &h,&a,&b,
                                 pni->error_message),
               pni->error_message,
               pni->error_message);

    temp_injection = array_spline_eval(pni->noninjection_table,
                                       pni->ddnoninjection_table,
                                       last_index_coarse,
                                       last_index_coarse+1,
                                       h,a,b);
    pni->photon_dep_table[index_z]+=temp_injection;
  }

  /** - Free temporary variables */
  free(pvecback);
  free(pvecthermo);

  return _SUCCESS_;
}

/**
 * Free the noninjection structure.
 *
 * @param pni   Input: pointer to noninjection structure
 * @return the error status
 */
int noninjection_free(struct noninjection* pni){

  free(pni->z_table);
  free(pni->photon_dep_table);

  free(pni->k);
  free(pni->k_weights);
  free(pni->pk_primordial_k);

  free(pni->z_table_coarse);
  free(pni->noninjection_table);
  free(pni->ddnoninjection_table);

  free(pni->integrand_approx);

  return _SUCCESS_;
}


/**
 * Interpolates photon injection from precomputed table at a given value of z.
 *
 * @param pni         Input: pointer to noninjected structure
 * @param z           Input: redshift
 * @param heat        Output: photon heating at given redshift
 * @return the error status
 */
int noninjection_photon_heating_at_z(struct noninjection* pni,
                                     double z,
                                     double* heat){

  /** Define local variables */
  double h,a,b;

  class_call(array_spline_hunt(pni->z_table,
                               pni->z_size,
                               z,
                               &(pni->last_index_z),
                               &h,&a,&b,
                               pni->error_message),
           pni->error_message,
           pni->error_message);

  * heat = ( a*pni->photon_dep_table[pni->last_index_z] + b*pni->photon_dep_table[pni->last_index_z+1] );

  return _SUCCESS_;
}


/**
 * Calculate heating from adiabatic cooling of electrons and baryons.
 *
 * @param pni           Input: pointer to noninjection structure
 * @param z             Input: redshift
 * @param energy_rate   Output: energy rate of non-injection process at given redshift
 * @return the error status
 */
int noninjection_rate_adiabatic_cooling(struct noninjection * pni,
                                        double z,
                                        double * energy_rate){

  /** Calculate heating rates */
  *energy_rate = -pni->heat_capacity*pni->H*pni->T_g;                                               // [J/(m^3 s)]

  return _SUCCESS_;

}


/**
 * Calculate heating from dissipation of acoustic waves.
 *
 * @param pni            Input: pointer to noninjection structure
 * @param z              Input: redshift
 * @param energy_rate    Output: energy rate of non-injection process at given redshift
 * @return the error status
  */
int noninjection_rate_acoustic_diss(struct noninjection * pni,
                                    double z,
                                    double * energy_rate){

  /** Define local variables */
  int index_k;
  double dQrho_dz;
  double A_wkb;

  /** a) Calculate full function */
  // CURRENTLY NOT YET IMPLEMENTED

  /** b) Calculate approximated function */
  if (_TRUE_){

    A_wkb = 1./(1.+4./15.*pni->f_nu_wkb);

    /* Define integrand for approximated function */
    for (index_k=0; index_k<pni->k_size; index_k++) {
      pni->integrand_approx[index_k] = 4.*A_wkb*A_wkb*pow(pni->k[index_k],1.)*
                              pni->pk_primordial_k[index_k]*
                              exp(-2.*pow(pni->k[index_k]/pni->kD,2.))*
                              pni->dkD_dz;
    }

    class_call(array_trapezoidal_integral(pni->integrand_approx,
                                          pni->k_size,
                                          pni->k_weights,
                                          &dQrho_dz,
                                          pni->error_message),
               pni->error_message,
               pni->error_message);

  }

  *energy_rate = dQrho_dz*pni->H*pni->rho_g/pni->a;                                                 // [J/(m^3 s)]

  return _SUCCESS_;
}

/**
 * Outputs
 */
int noninjection_output_titles(struct noninjection * pni, char titles[_MAXTITLESTRINGLENGTH_]){

  class_store_columntitle(titles,"Redshift z",_TRUE_);
  class_store_columntitle(titles,"Photon Heating [-]",_TRUE_);

  return _SUCCESS_;
}

int noninjection_output_data(struct noninjection * pni,
                             int number_of_titles,
                             double * data){
  int storeidx;
  double * dataptr;
  int index_z;

  for (index_z=0; index_z<pni->z_size; index_z++) {
    dataptr = data + index_z*number_of_titles;
    storeidx = 0;
    class_store_double(dataptr, pni->z_table[index_z], _TRUE_, storeidx);
    class_store_double(dataptr, pni->photon_dep_table[index_z], _TRUE_, storeidx);
  }

  return _SUCCESS_;
}
