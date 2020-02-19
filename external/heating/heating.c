/** @file heating.c Documented heating module
 *
 * Initially written by:
 * Nils Schoeneberg and Matteo Lucca, 27.02.2019
 *
 * The main goal of this module is to calculate the deposited energy in form of heating, ionization
 * and Lyman alpha processes. For more details see the description in the README file.
 */


#include "heating.h"
#include "primordial.h"


/**
 * Initialize heating table.
 *
 * @param ppr   Input: pointer to precision structure
 * @param pba   Input: pointer to background structure
 * @param pth   Input: pointer to thermodynamics structure
 * @return the error status
 */
int heating_init(struct precision * ppr,
                 struct background* pba,
                 struct thermo* pth){

  /** Define local variable */
  struct heating* phe = &(pth->he);
  int index_inj, index_dep;

  /** Initialize flags, indeces and parameters */
  phe->has_exotic_injection = _FALSE_;
  phe->has_DM_ann = _FALSE_;
  phe->has_DM_dec = _FALSE_;
  phe->has_PBH_eva = _FALSE_;
  phe->has_PBH_acc = _FALSE_;
  phe->last_index_x_chi = 0;
  phe->last_index_z_chi = 0;
  phe->last_index_z_feff = 0;

  /** Import qunatities from other structures */
  /* Precision structure */
  phe->Nz_size = ppr->thermo_Nz_lin;
  phe->z_initial = ppr->thermo_z_initial;
  phe->z_start_chi_approx = ppr->z_start_chi_approx;
  phe->Nz_PBH = ppr->primordial_black_hole_Nz;
  phe->heating_noninjected_Nz_log = ppr->heating_noninjected_Nz_log;


  /* Background structure */
  phe->H0 = pba->H0*_c_/_Mpc_over_m_;                                                               // [1/s]
  phe->T_g0 = pba->T_cmb;                                                                           // [K]
  phe->Omega0_b = pba->Omega0_b;                                                                    // [-]
  phe->Omega0_cdm = pba->Omega0_cdm;                                                                // [-]
  phe->rho0_cdm = pba->Omega0_cdm*pow(phe->H0,2)*3/8./_PI_/_G_*_c_*_c_;                             // [J/m^3]

  /* Thermodynamics structure */
  phe->fHe = pth->fHe;                                                                              // [-]
  phe->N_e0 = pth->n_e;                                                                             // [1/m^3]

  /** Define redshift tables */
  phe->z_size = pth->tt_size;
  class_alloc(phe->z_table,
              phe->z_size*sizeof(double),
              phe->error_message);
  memcpy(phe->z_table,
         pth->z_table,
         phe->z_size*sizeof(double));

  /** Define additional book-keeping variables for the z table */
  phe->tol_z_table = 1e-10;
  phe->filled_until_index_z = phe->z_size-1;
  phe->filled_until_z = phe->z_table[phe->filled_until_index_z];
  phe->last_index_z_chi = 0;
  phe->last_index_z_feff = 0;
  phe->last_index_z_inj = 0;
  phe->last_index_z = 0;
  phe->index_z_store = 0;

  /** Define indeces of tables */
  phe->to_store = _FALSE_;
  class_call(heating_indices(pth),
             phe->error_message,
             phe->error_message);

  /** Initialize energy injection table */
  /* Allocate space */
  class_alloc(phe->injection_table,
              phe->inj_size*sizeof(double*),
              phe->error_message);
  for(index_inj=0; index_inj<phe->inj_size; ++index_inj){
    class_alloc(phe->injection_table[index_inj],
                phe->z_size*sizeof(double),
                phe->error_message);
  }

  /* Calculate the PBH mass evolution, if needed */
  if(phe->has_PBH_eva == _TRUE_ ){
    class_call(heating_rate_PBH_evaporation_mass_evolution(pba,phe),
               phe->error_message,
               phe->error_message);
  }

  /** Initialize injection efficiency */
  /* Read from external file, if needed */
  if(phe->f_eff_type == f_eff_from_file){
    class_call(heating_read_feff_from_file(ppr,phe,
                                           phe->f_eff_file),
               phe->error_message,
               phe->error_message);
  }

  /** Initialize deposition function */
  /* Allocate space */
  class_alloc(phe->chi,
              phe->dep_size*sizeof(double),
              phe->error_message);

  /* Read from external file, if needed */
  if(phe->chi_type == chi_Galli_file){
    class_call(heating_read_chi_x_from_file(ppr,phe,
                                            ppr->chi_z_Galli),
               phe->error_message,
               phe->error_message);
  }
  else if(phe->chi_type == chi_from_x_file){
    class_call(heating_read_chi_x_from_file(ppr,phe,
                                            phe->chi_x_file),
               phe->error_message,
               phe->error_message);
  }
  else if(phe->chi_type == chi_from_z_file){
    class_call(heating_read_chi_z_from_file(ppr,phe,
                                            phe->chi_z_file),
               phe->error_message,
               phe->error_message);
  }

  /** Initialize energy deposition table */
  /* Allocate space */
  class_alloc(phe->deposition_table,
              phe->dep_size*sizeof(double*),
              phe->error_message);
  for(index_dep=0; index_dep<phe->dep_size; ++index_dep){
    class_alloc(phe->deposition_table[index_dep],
                phe->z_size*sizeof(double),
                phe->error_message);
  }
  class_alloc(phe->photon_dep_table,
              phe->z_size*sizeof(double),
              phe->error_message);

  class_alloc(phe->pvecdeposition,
              phe->dep_size*sizeof(double),
              phe->error_message);

  return _SUCCESS_;
}


/**
 * Initialize indeces of heating table.
 *
 * @param pth   Input: pointer to thermodynamics structure
 * @return the error status
 */
int heating_indices(struct thermo* pth){

  /** Define local variable */
  struct heating* phe = &(pth->he);
  int index_dep,index_inj;


  /* Check energy injection */
  if(phe->DM_annihilation_efficiency!=0 || phe->DM_decay_fraction!=0 || phe->PBH_evaporation_fraction!=0 || phe->PBH_accretion_fraction!=0){
    phe->has_exotic_injection = _TRUE_;
  }
  if(phe->DM_annihilation_efficiency!=0){
    phe->has_DM_ann = _TRUE_;
  }
  if(phe->DM_decay_fraction!=0){
    phe->has_DM_dec = _TRUE_;
  }
  if(phe->PBH_evaporation_fraction!=0){
    phe->has_PBH_eva = _TRUE_;
  }
  if(phe->PBH_accretion_fraction!=0){
    phe->has_PBH_acc = _TRUE_;
  }

  /** Indices for injection table */
  index_inj = 0;
  class_define_index(phe->index_inj_DM_ann  , phe->has_DM_ann  , index_inj, 1);
  class_define_index(phe->index_inj_DM_dec  , phe->has_DM_dec  , index_inj, 1);
  class_define_index(phe->index_inj_PBH_eva , phe->has_PBH_eva , index_inj, 1);
  class_define_index(phe->index_inj_PBH_acc , phe->has_PBH_acc , index_inj, 1);
  class_define_index(phe->index_inj_tot     , _TRUE_           , index_inj, 1);
  phe->inj_size = index_inj;

  /** Indices for deposition (and chi) table */
  index_dep = 0;
  class_define_index(phe->index_dep_heat , _TRUE_, index_dep, 1);
  class_define_index(phe->index_dep_ionH , _TRUE_, index_dep, 1);
  class_define_index(phe->index_dep_ionHe, _TRUE_, index_dep, 1);
  class_define_index(phe->index_dep_lya  , _TRUE_, index_dep, 1);
  class_define_index(phe->index_dep_lowE , _TRUE_, index_dep, 1);
  phe->dep_size = index_dep;

  return _SUCCESS_;
}


/**
 * Free allocated public variables and tables.
 *
 * @param pth   Input: pointer to thermodynamics structure
 * @return the error status
 */
int heating_free(struct thermo* pth){

  /** Define local variables */
  struct heating* phe = &(pth->he);
  int index_inj, index_dep;

  /* Redshift */
  free(phe->z_table);

  /* Energy injection */
  for(index_inj=0;index_inj<phe->inj_size;++index_inj){
    free(phe->injection_table[index_inj]);
  }
  free(phe->injection_table);

  if(phe->has_PBH_eva==_TRUE_){
    free(phe->PBH_table_z);
    free(phe->PBH_table_mass);
    free(phe->PBH_table_mass_dd);
    free(phe->PBH_table_F);
    free(phe->PBH_table_F_dd);
  }

  /* Injection efficiency */
  if(phe->f_eff_type == f_eff_from_file){
    free(phe->feff_table);
  }
  if(phe->chi_type == chi_from_z_file){
    free(phe->chiz_table);
  }
  if(phe->chi_type == chi_from_x_file || phe->chi_type == chi_Galli_file){
    free(phe->chix_table);
  }

  /* Deposition function */
  free(phe->chi);

  /* Energy deposition */
  for(index_dep=0;index_dep<phe->dep_size;++index_dep){
    free(phe->deposition_table[index_dep]);
  }
  free(phe->deposition_table);
  free(phe->photon_dep_table);
  free(phe->pvecdeposition);

  return _SUCCESS_;
}


int heating_noninjected_workspace_init(struct perturbs* ppt, struct primordial* ppm, struct heating* phe){

  struct heating_noninjected_workspace* niws = phe->noninjws;

  int index_k,index_z;
  /** Allocate new storage space that is a bit coarser in z as to reduce computational costs */
  niws->logz_max = log(1.+phe->z_table[phe->z_size-1]);
  niws->z_size_coarse = phe->heating_noninjected_Nz_log; /* This table is only going to be filled after the initial run anyway */
  class_alloc(niws->z_table_coarse,
              niws->z_size_coarse*sizeof(double),
              niws->error_message);
  for(index_z=0;index_z<niws->z_size_coarse-1;++index_z){
    niws->z_table_coarse[index_z] = exp(index_z*niws->logz_max/(niws->z_size_coarse-1))-1.;
  }
  niws->z_table_coarse[niws->z_size_coarse-1] = phe->z_table[phe->z_size-1];

  class_alloc(niws->injected_deposition,
              niws->z_size_coarse*sizeof(double),
              niws->error_message);
  class_alloc(niws->ddinjected_deposition,
              niws->z_size_coarse*sizeof(double),
              niws->error_message);

  /** Allocate qunatities from primordial structure */
  niws->k_max = 1.e6;
  niws->k_min = 0.12;
  niws->k_size = 500;        /* Found to be reasonable for the integral of acoustic dissipation, TODO :: make precision variable */
  class_alloc(niws->k,
              niws->k_size*sizeof(double),
              phe->error_message);
  class_alloc(niws->k_weights,
              niws->k_size*sizeof(double),
              phe->error_message);
  class_alloc(niws->pk_primordial_k,
              niws->k_size*sizeof(double),
              phe->error_message);

  /** Import primordial spectrum */
  for (index_k=0; index_k<niws->k_size; index_k++) {
    niws->k[index_k] = exp(log(niws->k_min)+(log(niws->k_max)-log(niws->k_min))/(niws->k_size)*index_k);
    class_call(primordial_spectrum_at_k(ppm,
                                        ppt->index_md_scalars,
                                        linear,
                                        niws->k[index_k],
                                        &niws->pk_primordial_k[index_k]),
               ppm->error_message,
               niws->error_message);
  }

  class_call(array_trapezoidal_weights(niws->k,
                                       niws->k_size,
                                       niws->k_weights,
                                       niws->error_message),
             niws->error_message,
             niws->error_message);

  if (phe->heating_rate_acoustic_diss_approx == _TRUE_){
    class_alloc(niws->integrand_approx,
                niws->k_size*sizeof(double),
                niws->error_message);
  }

  return _SUCCESS_;
}

int heating_noninjected_workspace_free(struct heating* phe){

  struct heating_noninjected_workspace* niws = phe->noninjws;

  free(niws->k);
  free(niws->k_weights);
  free(niws->pk_primordial_k);

  free(niws->z_table_coarse);
  free(niws->injected_deposition);
  free(niws->ddinjected_deposition);

  if (phe->heating_rate_acoustic_diss_approx == _TRUE_){
    free(niws->integrand_approx);
  }

  return _SUCCESS_;
}
/**
 * Calculate the heating (first order only) at the given redshift.
 * If phe->to_store is set to true, also store the value in
 * a table of heatings, which can later be used to interpolate.
 *
 * @param pba         Input: pointer to background structure
 * @param pth         Input: pointer to thermodynamics structure
 * @param x           Input: freaction of free electrons
 * @param z           Input: redshift
 * @param pvecback    Output: vector of background quantities
 * @return the error status
 */
int heating_calculate_at_z(struct background* pba,
                           struct thermo* pth,
                           double x,
                           double z,
                           double Tmat,
                           double* pvecback){

  /** Define local variables */
  struct heating* phe = &(pth->he);
  int index_dep, iz_store;
  double h,a,b;
  double dEdz_inj;

  /** Redefine input parameters */
  phe->T_b = Tmat;                                                                                  // [K]
  phe->x_e = x;                                                                                     // [-]
  phe->nH = phe->N_e0*pow(1.+z,3);                                                                  // [1/m^3]
  phe->heat_capacity = (3./2.)*_k_B_*phe->nH*(1.+phe->fHe+phe->x_e);                                // [J/(K m^3)]
  dEdz_inj = 0.;

  /** Import varying quantities from background structure */
  phe->H = pvecback[pba->index_bg_H]*_c_/_Mpc_over_m_;                                              // [1/s]
  phe->a = pvecback[pba->index_bg_a];                                                               // [-]
  phe->t = pvecback[pba->index_bg_time]/_s_over_Mpc_;                                               // [s]
  phe->rho_cdm = pvecback[pba->index_bg_rho_cdm]*_Jm3_over_Mpc2_;                                   // [J/m^3]
  phe->rho_g = pvecback[pba->index_bg_rho_g]*_Jm3_over_Mpc2_;                                       // [J/m^3]
  phe->rho_b = pvecback[pba->index_bg_rho_b]*_Jm3_over_Mpc2_;                                       // [J/m^3]

  /** Hunt within the redshift table for the given index of deposition */
  class_call(array_spline_hunt(phe->z_table,
                               phe->z_size,
                               z,
                               &(phe->index_z_store),
                               &h,&a,&b,
                               phe->error_message),
             phe->error_message,
             phe->error_message);

  /** Test if and where the new values should be stored in the injection table */
  /* If this value is important, store it */
  if(phe->to_store){
    /* Calculate where to store the value */
    if(fabs(b-1) < phe->tol_z_table){
      phe->index_z_store = phe->index_z_store+1;
    }
    else if(fabs(b) < phe->tol_z_table){
      phe->index_z_store = phe->index_z_store;
    }
    /* Could not find a matching index in the z table for this z */
    else{
      class_stop(phe->error_message,
                 "Should store z = %.10e, but it was not in the z table (next lower = %.10e , next higher = %.10e )",
                 phe->z_table[phe->index_z_store], phe->z_table[phe->index_z_store+1]);
    }
  }

  /** Get the injected energy that needs to be deposited (i.e. excluding adiabatic terms) */
  class_call(heating_energy_injection_at_z(phe,
                                           z,
                                           &dEdz_inj),
             phe->error_message,
             phe->error_message);

  /** Get the deposition and the efficiency functions */
  class_call(heating_deposition_function_at_z(phe,
                                              x,
                                              z),
             phe->error_message,
             phe->error_message);

  /** Put result into deposition vector */
  for(index_dep = 0; index_dep < phe->dep_size; ++index_dep){
    phe->pvecdeposition[index_dep] = phe->chi[index_dep]*dEdz_inj;
  }

  /** Store z values in table */
  if(phe->to_store){
    for(index_dep = 0; index_dep < phe->dep_size; ++index_dep){
      phe->deposition_table[index_dep][phe->index_z_store] = phe->pvecdeposition[index_dep];
    }
    phe->photon_dep_table[phe->index_z_store] = phe->pvecdeposition[phe->index_dep_heat]; /* All of the heating deposited into baryons also heats the photons */

    class_test(phe->index_z_store < phe->filled_until_index_z-1,
               phe->error_message,
               "Skipping too far ahead in z_table. Check that the heating and thermodynamics module agree in their z sampling.");
    phe->filled_until_index_z = phe->index_z_store;
    phe->filled_until_z = phe->z_table[phe->index_z_store];
  }

  phe->to_store = _FALSE_;

  return _SUCCESS_;
}


/**
 * Calculate energy injection at given redshift.
 *
 * @param phe         Input: pointer to heating structure
 * @param z           Input: redshift
 * @param dEdz_inj    Output: injected energy
 * @return the error status
 */
int heating_energy_injection_at_z(struct heating* phe,
                                  double z,
                                  double* dEdz_inj){

  /** Define local variable */
  double dEdz, rate;
  double h,a,b;
  int index_inj, iz_store;

  /* Initialize local variables */
  dEdz = 0.;

  /** Test if the values are already within the table */
  if(z > phe->filled_until_z == _TRUE_){
    /* If the value is already within the table, just interpolate */
    class_call(array_spline_hunt(phe->z_table,
                                 phe->z_size,
                                 z,
                                 &(phe->last_index_z_inj),
                                 &h,&a,&b,
                                 phe->error_message),
             phe->error_message,
             phe->error_message);
    /* (Linearly) interpolate within the table */
    for(index_inj=0; index_inj<phe->inj_size; ++index_inj){
      dEdz += phe->injection_table[index_inj][phe->last_index_z_inj]*a+phe->injection_table[index_inj][phe->last_index_z_inj+1]*b;
    }

    *dEdz_inj = dEdz;

    return _SUCCESS_;
  }

  /** Exotic energy injection mechanisms */
  if(phe->has_exotic_injection == _TRUE_){

    /* DM annihilation */
    if(phe->has_DM_ann == _TRUE_){
      class_call(heating_rate_DM_annihilation(phe,
                                              z,
                                              &rate),
                 phe->error_message,
                 phe->error_message);
      if(phe->to_store == _TRUE_){
        phe->injection_table[phe->index_inj_DM_ann][phe->index_z_store] = rate;
      }
      dEdz += rate;
    }

    /* DM decay */
    if(phe->has_DM_dec == _TRUE_){
      class_call(heating_rate_DM_decay(phe,
                                       z,
                                       &rate),
                 phe->error_message,
                 phe->error_message);
      if(phe->to_store == _TRUE_){
        phe->injection_table[phe->index_inj_DM_dec][phe->index_z_store] = rate;
      }
      dEdz += rate;
    }

    /* PBH evaporation */
    if(phe->has_PBH_eva == _TRUE_){
      class_call(heating_rate_PBH_evaporation(phe,
                                              z,
                                              &rate),
                 phe->error_message,
                 phe->error_message);
      if(phe->to_store == _TRUE_){
        phe->injection_table[phe->index_inj_PBH_eva][phe->index_z_store] = rate;
      }
      dEdz += rate;
    }

    /* PBH matter acctretion */
    if(phe->has_PBH_acc == _TRUE_){
      class_call(heating_rate_PBH_accretion(phe,
                                            z,
                                            &rate),
                 phe->error_message,
                 phe->error_message);
      if(phe->to_store == _TRUE_){
        phe->injection_table[phe->index_inj_PBH_acc][phe->index_z_store] = rate;
      }
      dEdz += rate;
    }
  }

  /** Total energy injection */
  if(phe->to_store == _TRUE_){
    phe->injection_table[phe->index_inj_tot][phe->index_z_store] = dEdz;

    class_test(phe->index_z_store < phe->filled_until_index_z-1,
               phe->error_message,
               "Skipping too far ahead in z_table. Check that the heating and thermodynamics modules agree in their z sampling.");
  }

  *dEdz_inj = dEdz;

  return _SUCCESS_;
}


/**
 * Calculate deposition function chi and injection efficiency f_eff at given redshift.
 *
 * @param phe         Input: pointer to heating structure
 * @param x           Input: fraction of free electrons
 * @param z           Input: redshift
 * @return the error status
 */
int heating_deposition_function_at_z(struct heating* phe,
                                     double x,
                                     double z){

  /** Define local variables */
  int index_dep;

  if(z > phe->z_start_chi_approx){
    /** In the verz early universe, whole energy goes into heating */
    phe->chi[phe->index_dep_heat]  = 1.;
    phe->chi[phe->index_dep_ionH]  = 0.;
    phe->chi[phe->index_dep_ionHe] = 0.;
    phe->chi[phe->index_dep_lya]   = 0.;
    phe->chi[phe->index_dep_lowE]  = 0.;
  }
  else{
    /** Read the deposition factors for each channel */
    /* Old approximation from Chen and Kamionkowski */
    if(phe->chi_type == chi_CK){
      if(x<1.){
        phe->chi[phe->index_dep_heat]  = (1.+2.*x)/3.;
        phe->chi[phe->index_dep_ionH]  = (1.-x)/3.;
        phe->chi[phe->index_dep_ionHe] = 0.;
        phe->chi[phe->index_dep_lya]   = (1.-x)/3.;
        phe->chi[phe->index_dep_lowE]  = 0.;
      }
      else{
        phe->chi[phe->index_dep_heat]  = 1.;
        phe->chi[phe->index_dep_ionH]  = 0.;
        phe->chi[phe->index_dep_ionHe] = 0.;
        phe->chi[phe->index_dep_lya]   = 0.;
        phe->chi[phe->index_dep_lowE]  = 0.;
      }
    }
    /* Old approximation from Padmanabhan and Finkbeiner */
    else if(phe->chi_type == chi_PF){
      if(x<1.+phe->fHe){
        phe->chi[phe->index_dep_heat]  = (1.+2.*x/(1+phe->fHe))/3.;
        phe->chi[phe->index_dep_ionH]  = (1.-x/(1+phe->fHe))/3.;
        phe->chi[phe->index_dep_ionHe] = 0.;
        phe->chi[phe->index_dep_lya]   = (1.-x/(1+phe->fHe))/3.;
        phe->chi[phe->index_dep_lowE]  = 0.;
      }
      else{
        phe->chi[phe->index_dep_heat]  = 1.;
        phe->chi[phe->index_dep_ionH]  = 0.;
        phe->chi[phe->index_dep_ionHe] = 0.;
        phe->chi[phe->index_dep_lya]   = 0.;
        phe->chi[phe->index_dep_lowE]  = 0.;
      }
    }
    /* Coefficient as revised by Galli et al. 2013 */
    else if(phe->chi_type == chi_Galli_file){
      for(index_dep=0; index_dep<phe->dep_size; ++index_dep){
        class_call(array_interpolate_spline_transposed(phe->chix_table,
                                                       phe->chix_size,
                                                       2*phe->dep_size+1,
                                                       0,
                                                       index_dep+1,
                                                       index_dep+phe->dep_size+1,
                                                       x,
                                                       &phe->last_index_x_chi,
                                                       &(phe->chi[index_dep]),
                                                       phe->error_message),
                   phe->error_message,
                   phe->error_message);
      }
    }
    /* coefficient as revised by Vivian Poulin (analytical interpolation of Galli et al. 2013) */
    else if(phe->chi_type == chi_Galli_analytic){
      if(x<1.){
        phe->chi[phe->index_dep_heat]  = MIN(0.996857*(1.-pow(1.-pow(x,0.300134),1.51035)),1);
        phe->chi[phe->index_dep_ionH]  = 0.369202*pow(1.-pow(x,0.463929),1.70237);
        phe->chi[phe->index_dep_ionHe] = 0.;
        phe->chi[phe->index_dep_lya]   = 0.;
        phe->chi[phe->index_dep_lowE]  = 0.;
      }
      else{
        phe->chi[phe->index_dep_heat]  = 1.;
        phe->chi[phe->index_dep_ionH]  = 0.;
        phe->chi[phe->index_dep_ionHe] = 0.;
        phe->chi[phe->index_dep_lya]   = 0.;
        phe->chi[phe->index_dep_lowE]  = 0.;
      }
    }
    else if(phe->chi_type == chi_full_heating){
      phe->chi[phe->index_dep_heat]  = 1.;
      phe->chi[phe->index_dep_ionH]  = 0.;
      phe->chi[phe->index_dep_ionHe] = 0.;
      phe->chi[phe->index_dep_lya]   = 0.;
      phe->chi[phe->index_dep_lowE]  = 0.;
    }
    /* Read file in ionization fraction */
    else if(phe->chi_type == chi_from_x_file){
      for(index_dep=0; index_dep<phe->dep_size; ++index_dep){
        class_call(array_interpolate_spline_transposed(phe->chix_table,
                                                       phe->chix_size,
                                                       2*phe->dep_size+1,
                                                       0,
                                                       index_dep+1,
                                                       index_dep+phe->dep_size+1,
                                                       x,
                                                       &phe->last_index_x_chi,
                                                       &(phe->chi[index_dep]),
                                                       phe->error_message),
                   phe->error_message,
                   phe->error_message);
      }
    }
    /* Read file in redshift */
    else if(phe->chi_type == chi_from_z_file){
      for(index_dep=0;index_dep<phe->dep_size;++index_dep){
        class_call(array_interpolate_spline_transposed(phe->chiz_table,
                                                       phe->chiz_size,
                                                       2*phe->dep_size+1,
                                                       0,
                                                       index_dep+1,
                                                       index_dep+phe->dep_size+1,
                                                       z,
                                                       &phe->last_index_z_chi,
                                                       &(phe->chi[index_dep]),
                                                       phe->error_message),
                   phe->error_message,
                   phe->error_message);
      }
    }
    else{
      class_stop(phe->error_message,"No valid deposition function has been found found.");
    }
  }

  /** Read the correction factor f_eff */
  /* For the on the spot, we take the user input */
  if(phe->f_eff_type == f_eff_on_the_spot){
    phe->f_eff = 1.;
  }
  /* For the file, read in f_eff from file and multiply */
  else if(phe->f_eff_type == f_eff_from_file){
    class_call(array_interpolate_spline_transposed(phe->feff_table,
                                                   phe->feff_z_size,
                                                   3,
                                                   0,
                                                   1,
                                                   2,
                                                   z,
                                                   &(phe->last_index_z_feff),
                                                   &(phe->f_eff),
                                                   phe->error_message),
               phe->error_message,
               phe->error_message);
  }
  /* Otherwise, something must have gone wrong */
  else{
    class_stop(phe->error_message,
               "Unknown energy deposition mechanism");
  }

  /** Multiply deposition factors with overall correction factor */
  for(index_dep=0; index_dep<phe->dep_size; ++index_dep){
    phe->chi[index_dep] *= phe->f_eff;
  }

  return _SUCCESS_;
}


/**
 * Update heating table with energy injection mechanisms,
 * which do not result from an energy injection which would be redistributed
 * in the primordial plasma
 *
 * @param pba   Input: pointer to background structure
 * @param pth   Input: pointer to thermodynamics structure
 * @param ppt   Input: pointer to perturbation structure
 * @param ppm   Input: pointer to primordial structure
 * @return the error status
 */
int heating_add_noninjected(struct background* pba,
                            struct thermo* pth,
                            struct perturbs* ppt,
                            struct primordial* ppm){

  /** Define local variables */
  struct heating* phe = &(pth->he);
  struct heating_noninjected_workspace* niws;
  int index_z;
  double tau;
  int last_index_back, last_index_thermo, last_index_coarse = 0;
  double *pvecback, *pvecthermo;
  double R, dkappa;
  int index_k;
  double dEdt;
  double z_wkb, tau_wkb;
  double h,a,b;

  /* Store the deposition in here */
  double temp_injection;
  double z_coarse;

  /** Allocate workspace to keep track of the arrays/values used throughout this function */
  class_alloc(phe->noninjws,sizeof(struct heating_noninjected_workspace),phe->error_message);
  class_call(heating_noninjected_workspace_init(ppt,ppm,phe),
             phe->noninjws->error_message,
             phe->error_message);

  /** Define useful shortcut (do AFTER allocating) */
  niws = phe->noninjws;

  /** Allocate backgorund and thermodynamcis vectors */
  last_index_back = 0;
  last_index_thermo = 0;
  class_alloc(pvecback,
              pba->bg_size*sizeof(double),
              phe->error_message);
  class_alloc(pvecthermo,
              pth->tt_size*sizeof(double),
              phe->error_message);

  /** Calculate WKB approximation ampltidue factor f_nu at early times */
  z_wkb = 1.0e6;              /* Found to be reasonable for wkb approximation */
  class_call(background_tau_of_z(pba,
                                 z_wkb,
                                 &tau_wkb),
             pba->error_message,
             phe->error_message);

  class_call(background_at_tau(pba,
                               tau_wkb,
                               pba->long_info,
                               pba->inter_normal,
                               &last_index_back,
                               pvecback),
             pba->error_message,
             phe->error_message);

  phe->f_nu_wkb = (1.-pvecback[pba->index_bg_rho_g]/(pvecback[pba->index_bg_Omega_r]*pvecback[pba->index_bg_rho_crit]));

  dEdt = 0.;
  /* Loop over z and calculate the heating at each point */
  for(index_z=0; index_z<niws->z_size_coarse; ++index_z){

    z_coarse = niws->z_table_coarse[index_z];
    niws->injected_deposition[index_z] = 0.;

    /** Import quantities from background and thermodynamics structure */
    class_call(background_tau_of_z(pba,
                                   z_coarse,
                                   &tau),
               pba->error_message,
               phe->error_message);

    class_call(background_at_tau(pba,
                                 tau,
                                 pba->long_info,
                                 pba->inter_closeby,
                                 &last_index_back,
                                 pvecback),
               pba->error_message,
               phe->error_message);

    phe->H = pvecback[pba->index_bg_H]*_c_/_Mpc_over_m_;                                            // [1/s]
    phe->a = pvecback[pba->index_bg_a];                                                             // [-]
    phe->rho_g = pvecback[pba->index_bg_rho_g]*_Jm3_over_Mpc2_;                                     // [J/m^3]
    R = (3./4.)*pvecback[pba->index_bg_rho_b]/pvecback[pba->index_bg_rho_g];                        // [-]
    phe->nH = phe->N_e0*pow(phe->a,-3);

    class_call(thermodynamics_at_z(pba,
                                   pth,
                                   z_coarse,
                                   pth->inter_normal,
                                   &last_index_thermo,
                                   pvecback,
                                   pvecthermo),
               pth->error_message,
               phe->error_message);

    dkappa = pvecthermo[pth->index_th_dkappa];                                                      // [1/Mpc]
    niws->dkD_dz = 1./(pvecback[pba->index_bg_H]*dkappa)*(16./15.+pow(R,2.)/(1.+R))/(6.*(1.0+R));   // [Mpc^2]
    niws->kD = 2.*_PI_/pvecthermo[pth->index_th_r_d];                                               // [1/Mpc]
    phe->T_b = pvecthermo[pth->index_th_Tb];                                                        // [K]
    phe->T_g = phe->T_g0*pow(phe->a,-1);                                                            // [K]
    phe->x_e = pvecthermo[pth->index_th_xe];                                                        // [-]
    phe->heat_capacity = (3./2.)*_k_B_*phe->nH*(1.+phe->fHe+phe->x_e);                              // [J/(K m^3)]

    /** Injected energy that does not need to be deposited (i.e. adiabatic terms) */
    /* First order cooling of photons due to adiabatic interaction with baryons */
    class_call(heating_rate_adiabatic_cooling(phe,
                                              z_coarse,
                                              &dEdt),
               phe->error_message,
               phe->error_message);
    niws->injected_deposition[index_z]+=dEdt;

    /* Second order acoustic dissipation of BAO */
    class_call(heating_rate_acoustic_diss(phe,
                                          z_coarse,
                                          &dEdt),
               phe->error_message,
               phe->error_message);
    niws->injected_deposition[index_z]+=dEdt;
  }

  class_call(array_spline_table_columns2(niws->z_table_coarse,
                                         niws->z_size_coarse,
                                         niws->injected_deposition,
                                         1,
                                         niws->ddinjected_deposition,
                                         _SPLINE_EST_DERIV_,
                                         phe->error_message),
             phe->error_message,
             phe->error_message);

  for(index_z=0;index_z<phe->z_size;++index_z){
    class_call(array_spline_hunt(niws->z_table_coarse,
                                 niws->z_size_coarse,
                                 phe->z_table[index_z],
                                 &last_index_coarse,
                                 &h,&a,&b,
                                 phe->error_message),
               phe->error_message,
               phe->error_message);

    temp_injection = array_interpolate_spline_hunt(niws->injected_deposition,
                                                   niws->ddinjected_deposition,
                                                   last_index_coarse,
                                                   last_index_coarse+1,
                                                   h,a,b);
    phe->photon_dep_table[index_z]+=temp_injection;
  }

  /* Free allocated space */
  free(pvecback);
  free(pvecthermo);

  class_call(heating_noninjected_workspace_free(phe),
             phe->noninjws->error_message,
             phe->error_message);
  free(phe->noninjws);

  return _SUCCESS_;
}


/**
 * Interpolates heating from precomputed table at a given value of z.
 *
 * @param pth         Input: pointer to thermodynamics structure
 * @param z           Input: redshift
 * @return the error status
 */
int heating_photon_at_z(struct thermo* pth,
                        double z,
                        double* heat){

  /** Define local variables */
  struct heating* phe = &(pth->he);
  int index_dep;
  double h,a,b;

  /** Interpolate at required z in the table */
  class_test(z < phe->filled_until_z,
             phe->error_message,
             "Heating is not yet calculated beyond %.10e (asked for at %.10e)",phe->filled_until_z,z);

  class_call(array_spline_hunt(phe->z_table,
                               phe->z_size,
                               z,
                               &(phe->last_index_z),
                               &h,&a,&b,
                               phe->error_message),
           phe->error_message,
           phe->error_message);

  * heat = ( a*phe->photon_dep_table[phe->last_index_z] + b*phe->photon_dep_table[phe->last_index_z+1] );

  return _SUCCESS_;
}


/**
 * Interpolates heating from precomputed table at a given value of z.
 *
 * @param pth         Input: pointer to thermodynamics structure
 * @param z           Input: redshift
 * @return the error status
 */
int heating_baryon_at_z(struct thermo* pth,
                        double z){

  /** Define local variables */
  struct heating* phe = &(pth->he);
  int index_dep;
  double h,a,b;

  /** Interpolate at required z in the table */
  class_test(z < phe->filled_until_z,
             phe->error_message,
             "Heating is not yet calculated beyond %.10e (asked for at %.10e)",phe->filled_until_z,z);

  class_call(array_spline_hunt(phe->z_table,
                               phe->z_size,
                               z,
                               &(phe->last_index_z),
                               &h,&a,&b,
                               phe->error_message),
           phe->error_message,
           phe->error_message);

  for(index_dep=0; index_dep<phe->dep_size; ++index_dep){
    phe->pvecdeposition[index_dep] = ( a*phe->deposition_table[index_dep][phe->last_index_z]+
                                        b*phe->deposition_table[index_dep][phe->last_index_z+1] );
  }

  return _SUCCESS_;
}


/**
 * Calculate heating from adiabatic cooling of electrons and baryons.
 *
 * @param phe        Input: pointer to heating structure
 * @param ppt        Input: pointer to the perturbations structure
 * @param ppm        Input: pointer to the primordial structure
 * @return the error status
 */
int heating_rate_adiabatic_cooling(struct heating * phe,
                                   double z,
                                   double * energy_rate){

  /** Define local variables */
  double R_g;

  /** Calculate heating rates */
  R_g = (2.*_sigma_/_m_e_/_c_)*(4./3.*phe->rho_g);

  //*energy_rate = R_g*phe->x_e/(1.+phe->x_e+phe->fHe)*(phe->T_b-phe->T_g)*phe->heat_capacity;      // [J/(m^3 s)]
  *energy_rate = -phe->heat_capacity*phe->H*phe->T_g;                                               // [J/(m^3 s)]

  return _SUCCESS_;

}


/**
 * Calculate heating from dissipation of acoustic waves.
 *
 * @param phe            Input: pointer to heating structure
 * @param z              Input: redshift
 * @param energy_rate    Output: energy density injection rate
 * @return the error status
  */
int heating_rate_acoustic_diss(struct heating * phe,
                               double z,
                               double * energy_rate){

  /** Define local variables */
  struct heating_noninjected_workspace* niws = phe->noninjws;
  int index_k;
  double dQrho_dz;
  double A_wkb;

  /** a) Calculate full function */
  if (phe->heating_rate_acoustic_diss_approx == _FALSE_){
    class_stop(phe->error_message,"Full calculation currently not implemented");
  }

  /** b) Calculate approximated function */
  if (phe->heating_rate_acoustic_diss_approx == _TRUE_){

    A_wkb = 1./(1.+4./15.*phe->f_nu_wkb);

    /* Define integrand for approximated function */
    for (index_k=0; index_k<niws->k_size; index_k++) {
      niws->integrand_approx[index_k] = 4.*A_wkb*A_wkb*pow(niws->k[index_k],1.)*
                              niws->pk_primordial_k[index_k]*
                              exp(-2.*pow(niws->k[index_k]/niws->kD,2.))*
                              niws->dkD_dz;
    }

    class_call(array_trapezoidal_integral(niws->integrand_approx,
                                          niws->k_size,
                                          niws->k_weights,
                                          &dQrho_dz,
                                          phe->error_message),
               phe->error_message,
               phe->error_message);

  }

  *energy_rate = dQrho_dz*phe->H*phe->rho_g/phe->a;                                                 // [J/(m^3 s)]

  return _SUCCESS_;
}


/**
 * Calculate heating from DM annihilation.
 *
 * @param phe            Input: pointer to heating structure
 * @param z              Input: redshift
 * @param energy_rate    Output: energy density injection rate
 * @return the error status
 */
int heating_rate_DM_annihilation(struct heating * phe,
                                 double z,
                                 double * energy_rate){

  /** Define local variables */
  double annihilation_at_z, boost_factor;

  /** Calculate chamge in the annihilation efficiency */
  if (z>phe->DM_annihilation_zmax) {
    annihilation_at_z = phe->DM_annihilation_efficiency*
                        exp(-phe->DM_annihilation_variation*pow(log((phe->DM_annihilation_z+1.)/(phe->DM_annihilation_zmax+1.)),2));
  }
  else if (z>phe->DM_annihilation_zmin) {
    annihilation_at_z = phe->DM_annihilation_efficiency*
                        exp(phe->DM_annihilation_variation*(-pow(log((phe->DM_annihilation_z+1.)/(phe->DM_annihilation_zmax+1.)),2)
                                         +pow(log((z+1.)/(phe->DM_annihilation_zmax+1.)),2)));
  }
  else {
    annihilation_at_z = phe->DM_annihilation_efficiency*
                        exp(phe->DM_annihilation_variation*(-pow(log((phe->DM_annihilation_z+1.)/(phe->DM_annihilation_zmax+1.)),2)
                                         +pow(log((phe->DM_annihilation_zmin+1.)/(phe->DM_annihilation_zmax+1.)),2)));
  }

  /** Calculate boost factor due to annihilation in halos */
  if(phe->DM_annihilation_z_halo > 0.){
    boost_factor = phe->DM_annihilation_f_halo * erfc((1+z)/(1+phe->DM_annihilation_z_halo)) / pow(1.+z,3);
  }
  else{
    boost_factor = 0;
  }

  /** Calculate heating rates */
  *energy_rate = pow(phe->rho_cdm,2.)*phe->DM_annihilation_efficiency*(1.+boost_factor);           // [J/(m^3 s)]

  return _SUCCESS_;
}


/**
 * Calculate heating from DM annihilation.
 *
 * @param phe            Input: pointer to heating structure
 * @param z              Input: redshift
 * @param energy_rate    Output: energy density injection rate
 * @return the error status
 */
int heating_rate_DM_decay(struct heating * phe,
                          double z,
                          double * energy_rate){

  /** Calculate heating rates */
  *energy_rate = phe->rho_cdm*phe->DM_decay_fraction*phe->DM_decay_Gamma*
                 exp(-phe->DM_decay_Gamma*phe->t);                                                  // [J/(m^3 s)]

  return _SUCCESS_;
}


/**
 * Determines time evolution of the primordial black hole (PBH) mass.
 * The conventions adopted here are the same as in Stoecker et al. 2018.
 *
 * @param pba            Input: pointer to background structure
 * @param phe            Input: pointer to thermodynamics structure
 * @return the error status
 */
int heating_rate_PBH_evaporation_mass_evolution(struct background * pba,
                                                struct heating * phe){

  /** Define local variables */
  double * pvecback_loop;
  int last_index_back_loop;
  int i_step;
  double current_mass, current_pbh_temperature;
  double f_EM, f_nu, f_q, f_pi, f_bos, f;
  double loop_z, loop_tau, time_now, time_prev, dt, dlnz, lnz_ini;

  /** Set initial parameters */
  current_mass = phe->PBH_evaporation_mass;                                                         // [g]
  phe->PBH_z_evaporation = 0;
  lnz_ini = log(1+phe->z_initial);
  dlnz = lnz_ini/(phe->Nz_PBH-1);
  loop_z = phe->z_initial*1.0001;
  time_prev = 0.;                                                                                   // [s]

  /** Alloate local variables */
  class_alloc(pvecback_loop,
              pba->bg_size*sizeof(double),
              phe->error_message);

  /** Alloate variables for PBH mass evolution */
  class_alloc(phe->PBH_table_z,
              phe->Nz_PBH*sizeof(double),
              phe->error_message);
  class_alloc(phe->PBH_table_mass,
              phe->Nz_PBH*sizeof(double),
              phe->error_message);
  class_alloc(phe->PBH_table_mass_dd,
              phe->Nz_PBH*sizeof(double),
              phe->error_message);
  class_alloc(phe->PBH_table_F,
              phe->Nz_PBH*sizeof(double),
              phe->error_message);
  class_alloc(phe->PBH_table_F_dd,
              phe->Nz_PBH*sizeof(double),
              phe->error_message);

  /** Fill tables with PBH mass evolution */
  /* For the parametrization of F(M) we follow PRD44 (1991) 376 with
   * the additional modification that we dress the "free QCD-particles"
   * (gluons and quarks) with an sigmoid-activation function (in log10-space:
   * Mean at 0.3 GeV and a width of 0.1*"order of magnitude") and the hadrons
   * with (1 - activation) to take the QCD-phase transition into account
   * and to be in agreement with PRD41 (1990) 3052, where the Ansatz is taken
   * that a black hole emmits those particles which appear elementary at the
   * given energy. */
  for(i_step = 0; i_step<phe->Nz_PBH; i_step++) {

    /** Find value of f(M) */
    current_pbh_temperature = 1.06e13/current_mass;                                                 // [GeV]
    phe->PBH_QCD_activation = 1./(1.+exp(-(log(current_pbh_temperature)-log(0.3))/(log(10.)*0.1))); // [-] see Eq. (4.6) of Stoecker et al. 2018

    f_EM = 2.*0.060                                                                                 // gamma
           +4.*0.142*exp(-(current_mass*5.11e-4)/(4.53*1.06e13))                                    // electron
           +4.*0.142*exp(-(current_mass*0.1037)/(4.53*1.06e13))                                     // muon
           +4.*0.142*exp(-(current_mass*1.777)/(4.53*1.06e13));                                     // tau
    f_nu = 6.*0.147;                                                                                // neutrino
    f_q = (12.*0.142*(exp(-(current_mass*2.2e-3)/(4.53*1.06e13))                                    // u
                      +exp(-(current_mass*4.7e-3)/(4.53*1.06e13))                                   // d
                      +exp(-(current_mass*1.28)/(4.53*1.06e13))                                     // c
                      +exp(-(current_mass*9.6e-2)/(4.53*1.06e13))                                   // s
                      +exp(-(current_mass*173.1)/(4.53*1.06e13))                                    // t
                      +exp(-(current_mass*4.18)/(4.53*1.06e13))                                     // b
                     )
           +16.*0.060*exp(-(current_mass*6e-1)/(6.04*1.06e13))                                      // g
          )*phe->PBH_QCD_activation;
    f_pi = (1.*0.267*exp(-(current_mass*1.350e-1)/(2.66*1.06e13))                                   // pi^0
            +2.*0.267*exp(-(current_mass*1.396e-1)/(2.66*1.06e13))                                  // pi^+
           )*(1-phe->PBH_QCD_activation);
    f_bos = 6.*0.060*exp(-(current_mass*80.39)/(6.04*1.06e13))                                      // W
            +3.*0.060*exp(-(current_mass*91.19)/(6.04*1.06e13))                                     // Z
            +1.*0.267*exp(-(current_mass*125.1)/(2.66*1.06e13));                                   // h
    f = f_EM+f_nu+f_q+f_pi+f_bos;

    /** Find current time value */
    class_call(background_tau_of_z(pba,
                                   loop_z,
                                   &loop_tau),
               pba->error_message,
               phe->error_message);
    class_call(background_at_tau(pba,
                                 loop_tau,
                                 pba->long_info,
                                 pba->inter_normal,
                                 &last_index_back_loop,
                                 pvecback_loop),
               pba->error_message,
               phe->error_message);
    time_now = pvecback_loop[pba->index_bg_time]/(_c_/_Mpc_over_m_);                                // [s]
    dt = time_now-time_prev;
    time_prev = time_now;

    if (i_step > 0) {
      //TODO :: check this step
      if (current_mass > 0.5*phe->PBH_evaporation_mass){
        current_mass = current_mass-5.34e25*f*pow(current_mass,-2)*dt;                              // [g]
      }
      else {
        if(phe->PBH_z_evaporation == 0){
          phe->PBH_z_evaporation = loop_z;
        }
        current_mass = 0.;
        f = 0.;
      }
    }

    /** Fill tables */
    phe->PBH_table_z[i_step] = loop_z;
    phe->PBH_table_mass[i_step] = current_mass;                                                     // [g]
    phe->PBH_table_F[i_step] = f;                                                                   // [-]
    loop_z = exp(lnz_ini-dlnz*(i_step+1))-1.;

  }

  /** Free allocated space */
  free(pvecback_loop);

  /** Spline mass and F(M) evolution in z */
  class_call(array_spline_table_lines(phe->PBH_table_z,
                                      phe->Nz_PBH,
                                      phe->PBH_table_mass,
                                      1,
                                      phe->PBH_table_mass_dd,
                                      _SPLINE_NATURAL_,
                                      phe->error_message),
             phe->error_message,
             phe->error_message);
  class_call(array_spline_table_lines(phe->PBH_table_z,
                                      phe->Nz_PBH,
                                      phe->PBH_table_F,
                                      1,
                                      phe->PBH_table_F_dd,
                                      _SPLINE_NATURAL_,
                                      phe->error_message),
             phe->error_message,
             phe->error_message);

  return _SUCCESS_;
}


/**
 * Calculate heating from PBH evaporation
 * The conventions adopted here are the same as in Stoecker et al. 2018.
 *
 * @param pba            Input: pointer to background structure
 * @param pth            Input: pointer to thermodynamics structure
 * @param z              Input: redshift
 * @param energy_rate    Output: energy density injection rate
 * @return the error status
 */
int heating_rate_PBH_evaporation(struct heating * phe,
                                 double z,
                                 double * energy_rate){

  /** Define local variables */
  int last_index_back;
  double mass, f;
  double dMdt, f_em;

  /** Interpolate the PBH mass evolution */
  class_call(array_interpolate_spline(phe->PBH_table_z,
                                      phe->Nz_PBH,
                                      phe->PBH_table_mass,
                                      phe->PBH_table_mass_dd,
                                      1,
                                      z,
                                      &last_index_back,
                                      &mass,                                                        // [g]
                                      1,
                                      phe->error_message),
             phe->error_message,
             phe->error_message);

  class_call(array_interpolate_spline(phe->PBH_table_z,
                                      phe->Nz_PBH,
                                      phe->PBH_table_F,
                                      phe->PBH_table_F_dd,
                                      1,
                                      z,
                                      &last_index_back,
                                      &f,                                                           // [-]
                                      1,
                                      phe->error_message),
             phe->error_message,
             phe->error_message);

  /** Calculate heating rates */
  if(mass <= 0.0001*phe->PBH_evaporation_mass || f <= 0 || z < phe->PBH_z_evaporation){
    *energy_rate = 0.;                                                                                // [J/(m^3 s)]
  }
  else {
    dMdt=5.34e25*f*pow(mass,-2.);                                                                     // [g/s]
    f_em = 0.55*phe->PBH_QCD_activation+(1-phe->PBH_QCD_activation)*(f-6.*0.147)/f;                   // [-]
    *energy_rate = phe->rho_cdm*phe->PBH_evaporation_fraction*f_em*dMdt/phe->PBH_evaporation_mass;    // [J/(m^3 s)]
  }

  return _SUCCESS_;
}


/**
 * Calculate heating from PBH matter accretion
 *
 * @param phe            Input: pointer to heating structure
 * @param z              Input: redshift
 * @param energy_rate    Output: energy density injection rate
 * @return the error status
 */
int heating_rate_PBH_accretion(struct heating * phe,
                               double z,
                               double * energy_rate){

  /** Define local variables */
  double tau, * pvecback;
  int last_index_back;
  double L_ed, M_ed_dot, M_crit, v_B, v_l, v_eff, r_B, t_B;
  double lambda, M_b_dot;
  double Value_min, Value_med, Value_max, a=0, epsilon_0=0.1, epsilon;
  double L_acc;
  double beta_compton_drag, gamma_cooling, tau_cooling;
  double T_ion, T_s, Y_s, theta_s;
  double lambda_1, lambda_2, lambda_ad, lambda_iso, J;

  /** Initialize local variables */
  /* Eddington luminosity */
  L_ed = 4.*_PI_*_G_*(phe->PBH_accretion_mass*_Sun_mass_)*_m_p_/_sigma_*_c_;                        // [W]
  M_ed_dot= 10.*L_ed/pow(_c_,2.);                                                                   // [kg/s]
  M_crit = 0.01*M_ed_dot;                                                                           // [kg/s]

  /* Boldi definitions */
  v_B = sqrt((1.+phe->x_e)*(phe->T_b*_k_B_)/(_m_p_*pow(_c_,2.)))*_c_;                               // [m/s]
  if(phe->PBH_accretion_relative_velocities < 0.){
    v_l = 30.*MIN(1.,(1.+z)/1.e3)*1.e3;                                                             // [m/s]
    if(v_B < v_l){
      v_eff = sqrt(v_B*v_l);                                                                        // [m/s]
    }
    else{
      v_eff = v_B;                                                                                  // [m/s]
    }
  }
  else{
    v_l = phe->PBH_accretion_relative_velocities*1.e3;                                              // [m/s]
    v_eff = pow(v_l*v_l+v_B*v_B,0.5);                                                               // [m/s]
  }
  r_B = _G_*(phe->PBH_accretion_mass*_Sun_mass_)/pow(v_eff,2.);                                     // [m]
  t_B = _G_*(phe->PBH_accretion_mass*_Sun_mass_)/pow(v_eff,3.);                                     // [s]

  /** Disk accretion from Poulin et al. 1707.04206 */
  if(phe->PBH_accretion_recipe == disk_accretion){

    lambda = phe->PBH_accretion_eigenvalue;                                                         // [-]
    M_b_dot = 4.*_PI_*lambda*(phe->rho_b/pow(_c_,2.))*pow(r_B,2.)*v_eff;                            // [kg/s]

    if(phe->PBH_accretion_ADAF_delta == 1e-3){
      Value_min = 7.6e-5;
      Value_med = 4.5e-3;
      Value_max = 7.1e-3;
      if(M_b_dot/M_ed_dot <= Value_min){
        epsilon_0 = 0.065;
        a = 0.71;
      }
      else if(Value_min < M_b_dot/M_ed_dot && M_b_dot/M_ed_dot  <= Value_med){
        epsilon_0 = 0.020;
        a = 0.47;
      }
      else if(Value_med < M_b_dot/M_ed_dot && M_b_dot/M_ed_dot <= Value_max){
        epsilon_0 = 0.26;
        a = 3.67;
      }
      else{
        epsilon_0 = 0.1;
        a = 0.;
      }
    }
    else if(phe->PBH_accretion_ADAF_delta == 0.1){
      Value_min = 9.4e-5;
      Value_med = 5e-3;
      Value_max = 6.6e-3;
      if(M_b_dot/M_ed_dot <= Value_min){
        epsilon_0 = 0.12;
        a = 0.59;
      }
      else if(Value_min < M_b_dot/M_ed_dot && M_b_dot/M_ed_dot  <= Value_med){
        epsilon_0 = 0.026;
        a = 0.27;
      }
      else if(Value_med < M_b_dot/M_ed_dot && M_b_dot/M_ed_dot <= Value_max){
        epsilon_0 = 0.50;
        a = 4.53;
      }
      else{
        epsilon_0 = 0.1;
        a = 0.;
      }
    }
    else if (phe->PBH_accretion_ADAF_delta == 0.5){
      Value_min = 2.9e-5;
      Value_med = 3.3e-3;
      Value_max = 5.3e-3;
      if(M_b_dot/M_ed_dot <= Value_min){
        epsilon_0 = 1.58;
        a = 0.65;
      }
      else if(Value_min < M_b_dot/M_ed_dot && M_b_dot/M_ed_dot  <= Value_med){
        epsilon_0 = 0.055;
        a = 0.076;
      }
      else if(Value_med < M_b_dot/M_ed_dot && M_b_dot/M_ed_dot <= Value_max){
        epsilon_0 = 0.17;
        a = 1.12;
      }
      else{
        epsilon_0 = 0.1;
        a = 0.;
      }
    }

    epsilon = epsilon_0*pow(M_b_dot/M_crit,a);                                                      // [-]
    L_acc = epsilon*M_b_dot*pow(_c_,2.);                                                            // [W]
  }

  /** Spherical accretion from Ali-Haimoud et al. 1612.05644 */
  else if(phe->PBH_accretion_recipe == spherical_accretion){

    beta_compton_drag = 4./3.*phe->x_e*_sigma_*phe->rho_g*t_B/_m_p_/_c_;                            // [-] Eq. (7)
    gamma_cooling = 2.*(_m_p_/_m_e_)/(1+phe->x_e)*beta_compton_drag;                                // [-] Eq. (8)
    tau_cooling = 1.5/(5.+pow(gamma_cooling,2./3.));                                                // [-] Eq. (28)

    lambda_ad = 0.25*pow(3./5.,1.5);                                                                // [-] Eq. (15)
    lambda_iso = 0.25*exp(1.5);                                                                     // [-] Eq. (20)
    lambda_1 = lambda_ad+(lambda_iso-lambda_ad)*
                                pow(pow(gamma_cooling,2.)/(88.+pow(gamma_cooling,2.)),0.22);        // [-] Eq. (27)
    lambda_2 = exp(4.5/(3.+pow(beta_compton_drag,0.75)))/(pow(pow(1.+beta_compton_drag,0.5)+1.,2.));// [-] Eq. (32)
    lambda = lambda_1*lambda_2/lambda_iso;                                                          // [-] Eq. (33)

    T_ion = 1.5e4*_eV_over_Kelvin_;                                                                 // [eV] see line below Eq. (39)
    Y_s = pow((1.+phe->x_e)/2.,2./3.*13.6/T_ion)*
                  tau_cooling/4.*pow(1.-5./2.*tau_cooling,1./3.)*(_m_p_/_m_e_);                     // [-] Eq. (51)
    T_s = _m_e_*pow(_c_,2.)/_k_B_*Y_s*pow(1.+Y_s/0.27,-1./3.);                                      // [K] Eqs. (50) and (47)
    theta_s = (T_s*_k_B_)/(_m_e_*pow(_c_,2.));

    M_b_dot = 4.*_PI_*lambda*(phe->rho_b/pow(_c_,2.))*pow(r_B,2.)*v_eff;                            // [kg/s] Eq. (6)

    if(theta_s > 1.){                                                                               // Eq. (55)
      J = 27./(2.*_PI_)*(log(2.*theta_s*exp(-0.577)+0.08)+4./3.);                                   // [-]
    }
    else{
      J = 4./_PI_*sqrt(2./_PI_)*pow(theta_s,-0.5)*(1+5.5*pow(theta_s,1.25));                        // [-]
    }

    L_acc = 1./137.*(T_s*_k_B_)/(_m_p_*pow(_c_,2.))*J*pow(M_b_dot*_c_*_c_,2.)/L_ed;                 // [W] Eq. (57)
  }

  *energy_rate = phe->rho_cdm*phe->PBH_accretion_fraction/
                      (phe->PBH_accretion_mass*_Sun_mass_)*L_acc/pow(_c_,2.);                       // [J/(m^3 s)]

  return _SUCCESS_;
}


/**
 * Read and interpolate the deposition function from external file.
 *
 * @param ppr   Input: pointer to precision structure
 * @param phe   Input/Output: pointer to heating structure
 * @return the error status
 */
int heating_read_feff_from_file(struct precision* ppr,
                                struct heating* phe,
                                char* f_eff_file){

  /** Define local variables */
  FILE * fA;
  char line[_LINE_LENGTH_MAX_];
  char * left;
  int headlines = 0;
  int index_z;

  phe->feff_z_size = 0;

  /* The file is assumed to contain:
   *    - The number of lines of the file
   *    - The columns ( z, f(z) ) where f(z) represents the "effective" fraction of energy deposited
   *      into the medium  at redshift z, in presence of halo formation. */
  class_open(fA, f_eff_file, "r", phe->error_message);

  while (fgets(line,_LINE_LENGTH_MAX_-1,fA) != NULL) {
    headlines++;

    /* Eliminate blank spaces at beginning of line */
    left=line;
    while (left[0]==' ') {
      left++;
    }

    /* Check that the line is neither blank nor a comment. In ASCII, left[0]>39 means that first non-blank charachter might
       be the beginning of some data (it is not a newline, a #, a %, etc.) */
    if (left[0] > 39) {

      /* If the line contains data, we must interprete it. If num_lines == 0 , the current line must contain
         its value. Otherwise, it must contain (xe , chi_heat, chi_Lya, chi_H, chi_He, chi_lowE). */

      /* Read num_lines, infer size of arrays and allocate them */
      class_test(sscanf(line,"%d",&(phe->feff_z_size)) != 1,
                 phe->error_message,
                 "could not read the initial integer of number of lines in line %i in file '%s' \n",
                 headlines,f_eff_file);

      /* (z, f, ddf)*/
      class_alloc(phe->feff_table,
                  3*phe->feff_z_size*sizeof(double),
                  phe->error_message);
      break;
    }
  }

  for(index_z=0;index_z<phe->feff_z_size;++index_z){
    /* Read coefficients */
    class_test(fscanf(fA,"%lg %lg",
                      &(phe->feff_table[index_z*3+0]),  // z
                      &(phe->feff_table[index_z*3+1])   // f_eff(z)
                     ) != 2,
               phe->error_message,
               "could not read value of parameters coefficients in line %i in file '%s'\n",
               headlines,f_eff_file);
  }

  fclose(fA);

  /* Spline in one dimension */
  class_call(array_spline(phe->feff_table,
                          3,
                          phe->feff_z_size,
                          0,
                          1,
                          2,
                          _SPLINE_NATURAL_,
                          phe->error_message),
             phe->error_message,
             phe->error_message);

  return _SUCCESS_;
}


/**
 * Read and interpolate the branching ratio from external file, if the function
 * in the file is given with respect to redshift.
 *
 * @param ppr   Input: pointer to precision structure
 * @param phe   Input/Output: pointer to heating structure
 * @return the error status
 */
int heating_read_chi_z_from_file(struct precision* ppr,
                                 struct heating* phe,
                                 char* chi_z_file){

  /** Define local variables */
  FILE * fA;
  char line[_LINE_LENGTH_MAX_];
  char * left;
  int headlines = 0;
  int index_z,index_dep;

  phe->chiz_size = 0;

  /* The file is assumed to contain:
   *    - The number of lines of the file
   *    - The columns (xe , chi_heat, chi_Lya, chi_H, chi_He, chi_lowE) where chi_i represents the
   *      branching ratio at redshift z into different heating/ionization channels i */

  class_open(fA, chi_z_file, "r", phe->error_message);

  while (fgets(line,_LINE_LENGTH_MAX_-1,fA) != NULL) {
    headlines++;

    /* Eliminate blank spaces at beginning of line */
    left=line;
    while (left[0]==' ') {
      left++;
    }

    /* Check that the line is neither blank nor a comment. In ASCII, left[0]>39 means that first non-blank charachter might
       be the beginning of some data (it is not a newline, a #, a %, etc.) */
    if (left[0] > 39) {

      /* If the line contains data, we must interprete it. If num_lines == 0 , the current line must contain
         its value. Otherwise, it must contain (xe , chi_heat, chi_Lya, chi_H, chi_He, chi_lowE). */

      /* Read num_lines, infer size of arrays and allocate them */
      class_test(sscanf(line,"%d",&(phe->chiz_size)) != 1,
                 phe->error_message,
                 "could not read the initial integer of number of lines in line %i in file '%s' \n",
                 headlines,chi_z_file);

      /* (z, chi_i)*/
      class_alloc(phe->chiz_table,
                  (2*phe->dep_size+1)*phe->chiz_size*sizeof(double),
                  phe->error_message);
      break;
    }
  }

  for(index_z=0;index_z<phe->chiz_size;++index_z){
    /* Read coefficients */
    class_test(fscanf(fA,"%lg %lg %lg %lg %lg %lg",
                      &(phe->chiz_table[index_z*(2*phe->dep_size+1)+0]), //z
                      &(phe->chiz_table[index_z*(2*phe->dep_size+1)+1+phe->index_dep_heat]), //heat
                      &(phe->chiz_table[index_z*(2*phe->dep_size+1)+1+phe->index_dep_lya]), //lya
                      &(phe->chiz_table[index_z*(2*phe->dep_size+1)+1+phe->index_dep_ionH]), //ionH
                      &(phe->chiz_table[index_z*(2*phe->dep_size+1)+1+phe->index_dep_ionHe]), //ionHe
                      &(phe->chiz_table[index_z*(2*phe->dep_size+1)+1+phe->index_dep_lowE])  //lowE
                     )!= 6,
               phe->error_message,
               "could not read value of parameters coefficients in line %i in file '%s'\n",
               index_z+headlines,chi_z_file);
  }

  fclose(fA);

  /* Spline in one dimension */
  for(index_dep=0;index_dep<phe->dep_size;++index_dep){
    class_call(array_spline(phe->chiz_table,
                            2*phe->dep_size+1,
                            phe->chiz_size,
                            0,
                            1+index_dep,
                            1+index_dep+phe->dep_size,
                            _SPLINE_NATURAL_,
                            phe->error_message),
               phe->error_message,
               phe->error_message);
  }

  return _SUCCESS_;
}


/**
 * Read and interpolate the branching ratio from external file, if the function
 * in the file is given with respect to the fraction of free electrons X_e.
 *
 * @param ppr   Input: pointer to precision structure
 * @param phe   Input/Output: pointer to heating structure
 * @return the error status
 */
int heating_read_chi_x_from_file(struct precision* ppr,
                                 struct heating* phe,
                                 char* chi_x_file){

  /** Define local variables */
  FILE * fA;
  char line[_LINE_LENGTH_MAX_];
  char * left;
  int headlines = 0;
  int index_x,index_dep;

  phe->chix_size = 0;

  /* The file is assumed to contain:
   *    - The number of lines of the file
   *    - The columns (xe , chi_heat, chi_Lya, chi_H, chi_He, chi_lowE) where chi_i represents the
   *      branching ratio at redshift z into different heating/ionization channels i */

  class_open(fA, chi_x_file, "r", phe->error_message);

  while (fgets(line,_LINE_LENGTH_MAX_-1,fA) != NULL) {
    headlines++;

    /* Eliminate blank spaces at beginning of line */
    left=line;
    while (left[0]==' ') {
      left++;
    }

    /* Check that the line is neither blank nor a comment. In ASCII, left[0]>39 means that first non-blank charachter might
       be the beginning of some data (it is not a newline, a #, a %, etc.) */
    if (left[0] > 39) {

      /* If the line contains data, we must interprete it. If num_lines == 0 , the current line must contain
         its value. Otherwise, it must contain (xe , chi_heat, chi_Lya, chi_H, chi_He, chi_lowE). */

      /* Read num_lines, infer size of arrays and allocate them */
      class_test(sscanf(line,"%d",&(phe->chix_size)) != 1,
                 phe->error_message,
                 "could not read the initial integer of number of lines in line %i in file '%s' \n",
                 headlines, chi_x_file);

      /* (z, chi_i)*/
      class_alloc(phe->chix_table,
                  (2*phe->dep_size+1)*phe->chix_size*sizeof(double),
                  phe->error_message);
      break;
    }
  }

  for(index_x = 0; index_x < phe->chix_size;++index_x){
    /* Read coefficients */
    class_test(fscanf(fA,"%lg %lg %lg %lg %lg %lg",
                      &(phe->chix_table[index_x*(2*phe->dep_size+1)+0]), //x
                      &(phe->chix_table[index_x*(2*phe->dep_size+1)+1+phe->index_dep_heat]), //heat
                      &(phe->chix_table[index_x*(2*phe->dep_size+1)+1+phe->index_dep_lya]), //lya
                      &(phe->chix_table[index_x*(2*phe->dep_size+1)+1+phe->index_dep_ionH]), //ionH
                      &(phe->chix_table[index_x*(2*phe->dep_size+1)+1+phe->index_dep_ionHe]), //ionHe
                      &(phe->chix_table[index_x*(2*phe->dep_size+1)+1+phe->index_dep_lowE])  //lowE
                     )!= 6,
               phe->error_message,
               "could not read value of parameters coefficients in line %i in file '%s'\n",
               index_x+headlines,chi_x_file);
  }

  fclose(fA);

  /* Spline in one dimension */
  for(index_dep=0;index_dep<phe->dep_size;++index_dep){
    class_call(array_spline(phe->chix_table,
                            2*phe->dep_size+1,
                            phe->chix_size,
                            0,
                            1+index_dep,
                            1+index_dep+phe->dep_size,
                            _SPLINE_NATURAL_,
                            phe->error_message),
               phe->error_message,
               phe->error_message);
  }

  return _SUCCESS_;
}




/**
 * Outputs
 */
int heating_output_titles(struct heating * phe, char titles[_MAXTITLESTRINGLENGTH_]){

  class_store_columntitle(titles,"Redshift z",_TRUE_);
  class_store_columntitle(titles,"Heat_photon [-]",_TRUE_);
  class_store_columntitle(titles,"Heat_baryon [-]",_TRUE_);
  class_store_columntitle(titles,"IonH_baryon [-]",_TRUE_);
  class_store_columntitle(titles,"IonHe_baryon [-]",_TRUE_);
  class_store_columntitle(titles,"Lya_baryon [-]",_TRUE_);
  class_store_columntitle(titles,"LH_photon [-]",_TRUE_);
  class_store_columntitle(titles,"LH_baryon [-]",_TRUE_);

  return _SUCCESS_;
}

int heating_output_data(struct heating * phe,
                        int number_of_titles,
                        double * data){
  int storeidx;
  double * dataptr;
  int index_z;

  for (index_z=0; index_z<phe->z_size; index_z++) {
    dataptr = data + index_z*number_of_titles;
    storeidx = 0;
    class_store_double(dataptr, phe->z_table[index_z], _TRUE_, storeidx);
    class_store_double(dataptr, phe->photon_dep_table[index_z], _TRUE_, storeidx);
    class_store_double(dataptr, phe->deposition_table[phe->index_dep_heat][index_z], _TRUE_, storeidx);
    class_store_double(dataptr, phe->deposition_table[phe->index_dep_ionH][index_z], _TRUE_, storeidx);
    class_store_double(dataptr, phe->deposition_table[phe->index_dep_ionHe][index_z], _TRUE_, storeidx);
    class_store_double(dataptr, phe->deposition_table[phe->index_dep_lya][index_z], _TRUE_, storeidx);
    class_store_double(dataptr, phe->photon_dep_table[index_z]*(1.+phe->z_table[index_z]), _TRUE_, storeidx);
    class_store_double(dataptr, phe->deposition_table[phe->index_dep_heat][index_z]*(1.+phe->z_table[index_z]), _TRUE_, storeidx);
  }

  return _SUCCESS_;
}
