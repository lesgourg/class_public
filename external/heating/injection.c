/** @file injection.c Documented exotic energy injection module
 *
 * written by Nils Schoeneberg and Matteo Lucca, 27.02.2019
 *
 * The main goal of this module is to calculate the deposited energy in form of heating, ionization
 * and Lyman alpha processes from exotic energy injection processes.
 * For more details see the description in the README file.
 */
#include "injection.h"
#include "thermodynamics.h"

/**
 * Initialize injection structure.
 *
 * @param ppr   Input: pointer to precision structure
 * @param pba   Input: pointer to background structure
 * @param pth   Input: pointer to thermodynamics structure
 * @return the error status
 */
int injection_init(struct precision * ppr,
                   struct background* pba,
                   struct thermodynamics* pth){

  /** Summary: */

  /** - Define local variable */
  struct injection* pin = &(pth->in);
  int index_inj, index_dep;

  /** - Initialize flags, indices and parameters */
  pin->has_DM_ann = _FALSE_;
  pin->has_DM_dec = _FALSE_;
  pin->has_PBH_eva = _FALSE_;
  pin->has_PBH_acc = _FALSE_;
  pin->last_index_x_chi = 0;
  pin->last_index_z_chi = 0;
  pin->last_index_z_feff = 0;

  /** - Import quantities from other structures */
  /* Precision structure */
  pin->Nz_size = ppr->thermo_Nz_lin;
  pin->z_initial = ppr->thermo_z_initial;
  pin->z_start_chi_approx = ppr->z_start_chi_approx;
  pin->Nz_PBH = ppr->primordial_black_hole_Nz;

  /* Background structure */
  pin->H0 = pba->H0*_c_/_Mpc_over_m_;                                                               // [1/s]
  pin->T_g0 = pba->T_cmb;                                                                           // [K]
  pin->Omega0_b = pba->Omega0_b;                                                                    // [-]
  pin->Omega0_cdm = pba->Omega0_cdm;                                                                // [-]
  pin->rho0_cdm = pba->Omega0_cdm*pow(pin->H0,2)*3/8./_PI_/_G_*_c_*_c_;                             // [J/m^3]

  /* Thermodynamics structure */
  pin->fHe = pth->fHe;                                                                              // [-]
  pin->N_e0 = pth->n_e;                                                                             // [1/m^3]

  /** - Define redshift tables */
  pin->z_size = pth->tt_size;
  class_alloc(pin->z_table,
              pin->z_size*sizeof(double),
              pin->error_message);
  memcpy(pin->z_table,
         pth->z_table,
         pin->z_size*sizeof(double));

  /** - Define additional book-keeping variables for the z table */
  pin->tol_z_table = 1e-10;
  pin->filled_until_index_z = pin->z_size-1;
  pin->filled_until_z = pin->z_table[pin->filled_until_index_z];
  pin->last_index_z_chi = 0;
  pin->last_index_z_feff = 0;
  pin->last_index_z_inj = 0;
  pin->last_index_z = 0;
  pin->index_z_store = 0;

  /** - Define indices of tables */
  pin->to_store = _FALSE_;
  class_call(injection_indices(pth),
             pin->error_message,
             pin->error_message);

  /** - Initialize energy injection table */
  /* Allocate space */
  class_alloc(pin->injection_table,
              pin->inj_size*sizeof(double*),
              pin->error_message);
  for(index_inj=0; index_inj<pin->inj_size; ++index_inj){
    class_alloc(pin->injection_table[index_inj],
                pin->z_size*sizeof(double),
                pin->error_message);
  }

  /* Calculate the PBH mass evolution, if needed */
  if(pin->has_PBH_eva == _TRUE_ ){
    class_call(injection_rate_PBH_evaporation_mass_evolution(pba,pin),
               pin->error_message,
               pin->error_message);
  }

  /** - Initialize injection efficiency */
  /* Read from external file, if needed */
  if(pin->f_eff_type == f_eff_from_file){
    class_call(injection_read_feff_from_file(ppr,pin,
                                             pin->f_eff_file),
               pin->error_message,
               pin->error_message);
  }

  /** - Initialize deposition function */
  /* Allocate space */
  class_alloc(pin->chi,
              pin->dep_size*sizeof(double),
              pin->error_message);

  /* Read from external file, if needed */
  if(pin->chi_type == chi_Galli_file){
    class_call(injection_read_chi_x_from_file(ppr,pin,
                                              ppr->chi_z_Galli),
               pin->error_message,
               pin->error_message);
  }
  else if(pin->chi_type == chi_from_x_file){
    class_call(injection_read_chi_x_from_file(ppr,pin,
                                              pin->chi_x_file),
               pin->error_message,
               pin->error_message);
  }
  else if(pin->chi_type == chi_from_z_file){
    class_call(injection_read_chi_z_from_file(ppr,pin,
                                              pin->chi_z_file),
               pin->error_message,
               pin->error_message);
  }

  /** - Initialize energy deposition table */
  /* Allocate space */
  class_alloc(pin->deposition_table,
              pin->dep_size*sizeof(double*),
              pin->error_message);
  for(index_dep=0; index_dep<pin->dep_size; ++index_dep){
    class_alloc(pin->deposition_table[index_dep],
                pin->z_size*sizeof(double),
                pin->error_message);
  }

  class_alloc(pin->pvecdeposition,
              pin->dep_size*sizeof(double),
              pin->error_message);

  return _SUCCESS_;
}


/**
 * Initialize indices of injection table.
 *
 * @param pth   Input: pointer to thermodynamics structure
 * @return the error status
 */
int injection_indices(struct thermodynamics* pth){

  /** - Define local variables */
  struct injection* pin = &(pth->in);
  int index_dep,index_inj;

  /* Check energy injection */
  if(pin->DM_annihilation_efficiency!=0){
    pin->has_DM_ann = _TRUE_;
  }
  if(pin->DM_decay_fraction!=0){
    pin->has_DM_dec = _TRUE_;
  }
  if(pin->PBH_evaporation_fraction!=0){
    pin->has_PBH_eva = _TRUE_;
  }
  if(pin->PBH_accretion_fraction!=0){
    pin->has_PBH_acc = _TRUE_;
  }

  /** - Indices for injection table */
  index_inj = 0;
  class_define_index(pin->index_inj_DM_ann  , pin->has_DM_ann  , index_inj, 1);
  class_define_index(pin->index_inj_DM_dec  , pin->has_DM_dec  , index_inj, 1);
  class_define_index(pin->index_inj_PBH_eva , pin->has_PBH_eva , index_inj, 1);
  class_define_index(pin->index_inj_PBH_acc , pin->has_PBH_acc , index_inj, 1);
  class_define_index(pin->index_inj_tot     , _TRUE_           , index_inj, 1);
  pin->inj_size = index_inj;

  /** - Indices for deposition (and chi) table */
  index_dep = 0;
  class_define_index(pin->index_dep_heat , _TRUE_, index_dep, 1);
  class_define_index(pin->index_dep_ionH , _TRUE_, index_dep, 1);
  class_define_index(pin->index_dep_ionHe, _TRUE_, index_dep, 1);
  class_define_index(pin->index_dep_lya  , _TRUE_, index_dep, 1);
  class_define_index(pin->index_dep_lowE , _TRUE_, index_dep, 1);
  pin->dep_size = index_dep;

  return _SUCCESS_;
}


/**
 * Free allocated public variables and tables.
 *
 * @param pth   Input: pointer to thermodynamics structure
 * @return the error status
 */
int injection_free(struct thermodynamics* pth){

  /** - Define local variables */
  struct injection* pin = &(pth->in);
  int index_inj, index_dep;

  /* Redshift */
  free(pin->z_table);

  /* Energy injection */
  for(index_inj=0;index_inj<pin->inj_size;++index_inj){
    free(pin->injection_table[index_inj]);
  }
  free(pin->injection_table);

  if(pin->has_PBH_eva==_TRUE_){
    free(pin->PBH_table_z);
    free(pin->PBH_table_mass);
    free(pin->PBH_table_mass_dd);
    free(pin->PBH_table_F);
    free(pin->PBH_table_F_dd);
  }

  /* Injection efficiency */
  if(pin->f_eff_type == f_eff_from_file){
    free(pin->feff_table);
  }
  if(pin->chi_type == chi_from_z_file){
    free(pin->chiz_table);
  }
  if(pin->chi_type == chi_from_x_file || pin->chi_type == chi_Galli_file){
    free(pin->chix_table);
  }

  /* Deposition function */
  free(pin->chi);

  /* Energy deposition */
  for(index_dep=0;index_dep<pin->dep_size;++index_dep){
    free(pin->deposition_table[index_dep]);
  }
  free(pin->deposition_table);
  free(pin->pvecdeposition);

  return _SUCCESS_;
}


/**
 * Calculate the injection (first order only) at the given redshift.
 * the results are stored in the deposition vector
 * pth->pin->pvecdeposition, and if pin->to_store is set to true, also
 * in one line of the table pth->pin->deposition_table, which can
 * later be used to interpolate.
 *
 * @param pba         Input: pointer to background structure
 * @param pth         Input/output: pointer to thermodynamics structure
 * @param x           Input: freaction of free electrons
 * @param z           Input: redshift
 * @param pvecback    Input: vector of background quantities
 * @return the error status
 */
int injection_calculate_at_z(struct background* pba,
                             struct thermodynamics* pth,
                             double x,
                             double z,
                             double Tmat,
                             double* pvecback){

  /** - Define local variables */
  struct injection * pin = &(pth->in);
  int index_dep;
  double h,a,b;
  double dEdt_inj;

  /** - Store input parameters in struct */
  pin->T_b = Tmat;                                                                                  // [K]
  pin->x_e = x;                                                                                     // [-]
  pin->nH = pin->N_e0*pow(1.+z,3);                                                                  // [1/m^3]
  pin->heat_capacity = (3./2.)*_k_B_*pin->nH*(1.+pin->fHe+pin->x_e);                                // [J/(K m^3)]
  dEdt_inj = 0.;

  /** - Import varying quantities from background structure, convert to SI units */
  pin->H = pvecback[pba->index_bg_H]*_c_/_Mpc_over_m_;                                              // [1/s]
  pin->a = pvecback[pba->index_bg_a];                                                               // [-]
  pin->t = pvecback[pba->index_bg_time]/_s_over_Mpc_;                                               // [s]
  pin->rho_cdm = pvecback[pba->index_bg_rho_cdm]*_Jm3_over_Mpc2_;                                   // [J/m^3]
  pin->rho_g = pvecback[pba->index_bg_rho_g]*_Jm3_over_Mpc2_;                                       // [J/m^3]
  pin->rho_b = pvecback[pba->index_bg_rho_b]*_Jm3_over_Mpc2_;                                       // [J/m^3]

  /** - Hunt within the redshift table for the given index of deposition */
  class_call(array_spline_hunt(pin->z_table,
                               pin->z_size,
                               z,
                               &(pin->index_z_store),
                               &h,&a,&b,
                               pin->error_message),
             pin->error_message,
             pin->error_message);

  /** - Test if and where the new values should be stored in the injection table */
  /* If this value is important, store it */
  if(pin->to_store == _TRUE_){
    /* Calculate where to store the value */
    if(fabs(b-1) < pin->tol_z_table){
      pin->index_z_store = pin->index_z_store+1;
    }
    else if(fabs(b) < pin->tol_z_table){
      pin->index_z_store = pin->index_z_store;
    }
    /* Could not find a matching index in the z table for this z */
    else{
      class_stop(pin->error_message,
                 "Should store z = %.10e, but it was not in the z table (next lower = %.10e , next higher = %.10e )",
                 pin->z_table[pin->index_z_store], pin->z_table[pin->index_z_store+1]);
    }
  }

  /** - Get the injected energy that needs to be deposited (i.e. excluding adiabatic terms) */
  class_call(injection_energy_injection_at_z(pin,
                                             z,
                                             &dEdt_inj),
             pin->error_message,
             pin->error_message);

  /** - Get the deposition and the efficiency functions */
  class_call(injection_deposition_function_at_z(pin,
                                                x,
                                                z),
             pin->error_message,
             pin->error_message);

  /** - Put result into deposition vector */
  for(index_dep = 0; index_dep < pin->dep_size; ++index_dep){
    pin->pvecdeposition[index_dep] = pin->chi[index_dep]*dEdt_inj;
  }

  /** - Store z values in table */
  if(pin->to_store == _TRUE_){
    for(index_dep = 0; index_dep < pin->dep_size; ++index_dep){
      pin->deposition_table[index_dep][pin->index_z_store] = pin->pvecdeposition[index_dep];
    }

    class_test(pin->index_z_store < pin->filled_until_index_z-1,
               pin->error_message,
               "Skipping too far ahead in z_table. Check that the injection and thermodynamics module agree in their z sampling.");
    pin->filled_until_index_z = pin->index_z_store;
    pin->filled_until_z = pin->z_table[pin->index_z_store];
  }

  pin->to_store = _FALSE_;

  return _SUCCESS_;
}


/**
 * Calculate energy injection at given redshift.
 *
 * @param pin         Input: pointer to injection structure
 * @param z           Input: redshift
 * @param dEdt_inj    Output: injected energy
 * @return the error status
 */
int injection_energy_injection_at_z(struct injection* pin,
                                    double z,
                                    double* dEdt_inj){

  /** - Define local variable */
  double dEdt, rate;
  double h,a,b;
  int index_inj;

  /* Initialize local variables */
  dEdt = 0.;

  /** - Test if the values are already within the table */
  if(z > pin->filled_until_z){
    /* If the value is already within the table, just interpolate */
    class_call(array_spline_hunt(pin->z_table,
                                 pin->z_size,
                                 z,
                                 &(pin->last_index_z_inj),
                                 &h,&a,&b,
                                 pin->error_message),
             pin->error_message,
             pin->error_message);
    /* (Linearly) interpolate within the table */
    for(index_inj=0; index_inj<pin->inj_size; ++index_inj){
      dEdt += pin->injection_table[index_inj][pin->last_index_z_inj]*a+pin->injection_table[index_inj][pin->last_index_z_inj+1]*b;
    }

    *dEdt_inj = dEdt;

  }

  else {
    /** - Exotic energy injection mechanisms */
    /* DM annihilation */
    if(pin->has_DM_ann == _TRUE_){
      class_call(injection_rate_DM_annihilation(pin,
                                                z,
                                                &rate),
                 pin->error_message,
                 pin->error_message);
      if(pin->to_store == _TRUE_){
        pin->injection_table[pin->index_inj_DM_ann][pin->index_z_store] = rate;
      }
      dEdt += rate;
    }

    /* DM decay */
    if(pin->has_DM_dec == _TRUE_){
      class_call(injection_rate_DM_decay(pin,
                                         z,
                                         &rate),
                 pin->error_message,
                 pin->error_message);
      if(pin->to_store == _TRUE_){
        pin->injection_table[pin->index_inj_DM_dec][pin->index_z_store] = rate;
      }
      dEdt += rate;
    }

    /* PBH evaporation */
    if(pin->has_PBH_eva == _TRUE_){
      class_call(injection_rate_PBH_evaporation(pin,
                                                z,
                                                &rate),
                 pin->error_message,
                 pin->error_message);
      if(pin->to_store == _TRUE_){
        pin->injection_table[pin->index_inj_PBH_eva][pin->index_z_store] = rate;
      }
      dEdt += rate;
    }

    /* PBH matter acctretion */
    if(pin->has_PBH_acc == _TRUE_){
      class_call(injection_rate_PBH_accretion(pin,
                                              z,
                                              &rate),
                 pin->error_message,
                 pin->error_message);
      if(pin->to_store == _TRUE_){
        pin->injection_table[pin->index_inj_PBH_acc][pin->index_z_store] = rate;
      }
      dEdt += rate;
    }

    /** Total energy injection */
    if(pin->to_store == _TRUE_){
      pin->injection_table[pin->index_inj_tot][pin->index_z_store] = dEdt;

      class_test(pin->index_z_store < pin->filled_until_index_z-1,
                 pin->error_message,
                 "Skipping too far ahead in z_table. Check that the injection and thermodynamics modules agree in their z sampling.");
    }

    *dEdt_inj = dEdt;

  }

  return _SUCCESS_;
}


/**
 * Calculate deposition function chi and injection efficiency f_eff at given redshift.
 *
 * @param pin         Input: pointer to injection structure
 * @param x           Input: fraction of free electrons
 * @param z           Input: redshift
 * @return the error status
 */
int injection_deposition_function_at_z(struct injection* pin,
                                       double x,
                                       double z){

  /** - Define local variables */
  int index_dep;

  if(z > pin->z_start_chi_approx){
    /** - In the verz early universe, whole energy goes into injection */
    pin->chi[pin->index_dep_heat]  = 1.;
    pin->chi[pin->index_dep_ionH]  = 0.;
    pin->chi[pin->index_dep_ionHe] = 0.;
    pin->chi[pin->index_dep_lya]   = 0.;
    pin->chi[pin->index_dep_lowE]  = 0.;
  }
  else{
    /** - Use the deposition factors for each channel */
    /* Old approximation from Chen and Kamionkowski */
    if(pin->chi_type == chi_CK){
      if(x<1.){
        pin->chi[pin->index_dep_heat]  = (1.+2.*x)/3.;
        pin->chi[pin->index_dep_ionH]  = (1.-x)/3.;
        pin->chi[pin->index_dep_ionHe] = 0.;
        pin->chi[pin->index_dep_lya]   = (1.-x)/3.;
        pin->chi[pin->index_dep_lowE]  = 0.;
      }
      else{
        pin->chi[pin->index_dep_heat]  = 1.;
        pin->chi[pin->index_dep_ionH]  = 0.;
        pin->chi[pin->index_dep_ionHe] = 0.;
        pin->chi[pin->index_dep_lya]   = 0.;
        pin->chi[pin->index_dep_lowE]  = 0.;
      }
    }
    /* Old approximation from Padmanabhan and Finkbeiner */
    else if(pin->chi_type == chi_PF){
      if(x<1.+pin->fHe){
        pin->chi[pin->index_dep_heat]  = (1.+2.*x/(1+pin->fHe))/3.;
        pin->chi[pin->index_dep_ionH]  = (1.-x/(1+pin->fHe))/3.;
        pin->chi[pin->index_dep_ionHe] = 0.;
        pin->chi[pin->index_dep_lya]   = (1.-x/(1+pin->fHe))/3.;
        pin->chi[pin->index_dep_lowE]  = 0.;
      }
      else{
        pin->chi[pin->index_dep_heat]  = 1.;
        pin->chi[pin->index_dep_ionH]  = 0.;
        pin->chi[pin->index_dep_ionHe] = 0.;
        pin->chi[pin->index_dep_lya]   = 0.;
        pin->chi[pin->index_dep_lowE]  = 0.;
      }
    }
    /* Coefficient as revised by Galli et al. 2013 */
    else if(pin->chi_type == chi_Galli_file){
      for(index_dep=0; index_dep<pin->dep_size; ++index_dep){
        class_call(array_interpolate_spline_transposed(pin->chix_table,
                                                       pin->chix_size,
                                                       2*pin->dep_size+1,
                                                       0,
                                                       index_dep+1,
                                                       index_dep+pin->dep_size+1,
                                                       x,
                                                       &pin->last_index_x_chi,
                                                       &(pin->chi[index_dep]),
                                                       pin->error_message),
                   pin->error_message,
                   pin->error_message);
      }
    }
    /* coefficient as revised by Vivian Poulin (analytical interpolation of Galli et al. 2013) */
    else if(pin->chi_type == chi_Galli_analytic){
      if(x<1.){
        pin->chi[pin->index_dep_heat]  = MIN(0.996857*(1.-pow(1.-pow(x,0.300134),1.51035)),1);
        pin->chi[pin->index_dep_ionH]  = 0.369202*pow(1.-pow(x,0.463929),1.70237);
        pin->chi[pin->index_dep_ionHe] = 0.;
        pin->chi[pin->index_dep_lya]   = 0.;
        pin->chi[pin->index_dep_lowE]  = 0.;
      }
      else{
        pin->chi[pin->index_dep_heat]  = 1.;
        pin->chi[pin->index_dep_ionH]  = 0.;
        pin->chi[pin->index_dep_ionHe] = 0.;
        pin->chi[pin->index_dep_lya]   = 0.;
        pin->chi[pin->index_dep_lowE]  = 0.;
      }
    }
    else if(pin->chi_type == chi_full_heating){
      pin->chi[pin->index_dep_heat]  = 1.;
      pin->chi[pin->index_dep_ionH]  = 0.;
      pin->chi[pin->index_dep_ionHe] = 0.;
      pin->chi[pin->index_dep_lya]   = 0.;
      pin->chi[pin->index_dep_lowE]  = 0.;
    }
    /* Read file in ionization fraction */
    else if(pin->chi_type == chi_from_x_file){
      for(index_dep=0; index_dep<pin->dep_size; ++index_dep){
        class_call(array_interpolate_spline_transposed(pin->chix_table,
                                                       pin->chix_size,
                                                       2*pin->dep_size+1,
                                                       0,
                                                       index_dep+1,
                                                       index_dep+pin->dep_size+1,
                                                       x,
                                                       &pin->last_index_x_chi,
                                                       &(pin->chi[index_dep]),
                                                       pin->error_message),
                   pin->error_message,
                   pin->error_message);
      }
    }
    /* Read file in redshift */
    else if(pin->chi_type == chi_from_z_file){
      for(index_dep=0;index_dep<pin->dep_size;++index_dep){
        class_call(array_interpolate_spline_transposed(pin->chiz_table,
                                                       pin->chiz_size,
                                                       2*pin->dep_size+1,
                                                       0,
                                                       index_dep+1,
                                                       index_dep+pin->dep_size+1,
                                                       z,
                                                       &pin->last_index_z_chi,
                                                       &(pin->chi[index_dep]),
                                                       pin->error_message),
                   pin->error_message,
                   pin->error_message);
      }
    }
    else{
      class_stop(pin->error_message,"No valid deposition function has been found found.");
    }
  }

  /** - Read the correction factor f_eff */
  /* For the on the spot, we take the user input */
  if(pin->f_eff_type == f_eff_on_the_spot){
    // pin->f_eff has already been seet by user
  }
  /* For the file, read in f_eff from file and multiply */
  else if(pin->f_eff_type == f_eff_from_file){
    class_call(array_interpolate_spline_transposed(pin->feff_table,
                                                   pin->feff_z_size,
                                                   3,
                                                   0,
                                                   1,
                                                   2,
                                                   z,
                                                   &(pin->last_index_z_feff),
                                                   &(pin->f_eff),
                                                   pin->error_message),
               pin->error_message,
               pin->error_message);
  }
  /* Otherwise, something must have gone wrong */
  else{
    class_stop(pin->error_message,
               "Unknown energy deposition mechanism");
  }

  /** Multiply deposition factors with overall correction factor */
  for(index_dep=0; index_dep<pin->dep_size; ++index_dep){
    pin->chi[index_dep] *= pin->f_eff;
  }

  return _SUCCESS_;
}


/**
 * Interpolates deposition from precomputed table at a given value of z.
 *
 * @param pth         Input: pointer to thermodynamics structure
 * @param z           Input: redshift
 * @return the error status
 */
int injection_deposition_at_z(struct thermodynamics* pth,
                              double z){

  /** - Define local variables */
  struct injection* pin = &(pth->in);
  int index_dep;
  double h,a,b;

  /** - Interpolate at required z in the table */
  class_test(z < pin->filled_until_z,
             pin->error_message,
             "injection is not yet calculated beyond %.10e (asked for at %.10e)",pin->filled_until_z,z);

  class_call(array_spline_hunt(pin->z_table,
                               pin->z_size,
                               z,
                               &(pin->last_index_z),
                               &h,&a,&b,
                               pin->error_message),
           pin->error_message,
           pin->error_message);

  for(index_dep=0; index_dep<pin->dep_size; ++index_dep){
    pin->pvecdeposition[index_dep] = ( a*pin->deposition_table[index_dep][pin->last_index_z]+
                                        b*pin->deposition_table[index_dep][pin->last_index_z+1] );
  }

  return _SUCCESS_;
}




/**
 * Calculate injection from DM annihilation.
 *
 * @param pin            Input: pointer to injection structure
 * @param z              Input: redshift
 * @param energy_rate    Output: energy density injection rate
 * @return the error status
 */
int injection_rate_DM_annihilation(struct injection * pin,
                                   double z,
                                   double * energy_rate){

  /** - Define local variables */
  double annihilation_at_z, boost_factor;

  /** - Calculate change in the annihilation efficiency */
  if (z>pin->DM_annihilation_zmax) {
    annihilation_at_z = pin->DM_annihilation_efficiency*
                        exp(-pin->DM_annihilation_variation*pow(log((pin->DM_annihilation_z+1.)/(pin->DM_annihilation_zmax+1.)),2));
  }
  else if (z>pin->DM_annihilation_zmin) {
    annihilation_at_z = pin->DM_annihilation_efficiency*
                        exp(pin->DM_annihilation_variation*(-pow(log((pin->DM_annihilation_z+1.)/(pin->DM_annihilation_zmax+1.)),2)
                                         +pow(log((z+1.)/(pin->DM_annihilation_zmax+1.)),2)));
  }
  else {
    annihilation_at_z = pin->DM_annihilation_efficiency*
                        exp(pin->DM_annihilation_variation*(-pow(log((pin->DM_annihilation_z+1.)/(pin->DM_annihilation_zmax+1.)),2)
                                         +pow(log((pin->DM_annihilation_zmin+1.)/(pin->DM_annihilation_zmax+1.)),2)));
  }

  /** - Calculate boost factor due to annihilation in halos */
  if(pin->DM_annihilation_z_halo > 0.){
    boost_factor = pin->DM_annihilation_f_halo * erfc((1+z)/(1+pin->DM_annihilation_z_halo)) / pow(1.+z,3);
  }
  else{
    boost_factor = 0;
  }

  /** - Calculate injection rates */
  *energy_rate = pow(pin->rho_cdm,2.)*annihilation_at_z*(1.+boost_factor);           // [J/(m^3 s)]

  return _SUCCESS_;
}


/**
 * Calculate injection from DM decay.
 *
 * @param pin            Input: pointer to injection structure
 * @param z              Input: redshift
 * @param energy_rate    Output: energy density injection rate
 * @return the error status
 */
int injection_rate_DM_decay(struct injection * pin,
                            double z,
                            double * energy_rate){

  /** - Calculate injection rates */
  *energy_rate = pin->rho_cdm*pin->DM_decay_fraction*pin->DM_decay_Gamma*
                 exp(-pin->DM_decay_Gamma*pin->t);                                                  // [J/(m^3 s)]

  return _SUCCESS_;
}


/**
 * Determines time evolution of the primordial black hole (PBH) mass.
 * The conventions adopted here are the same as in Stoecker et al. 2018.
 *
 * @param pba            Input: pointer to background structure
 * @param pin            Input: pointer to thermodynamics structure
 * @return the error status
 */
int injection_rate_PBH_evaporation_mass_evolution(struct background * pba,
                                                  struct injection * pin){

  /** - Define local variables */
  double * pvecback_loop;
  int last_index_back_loop;
  int i_step;
  double current_mass, current_pbh_temperature;
  double f_EM, f_nu, f_q, f_pi, f_bos, f;
  double loop_z, time_now, time_prev, dt, dlnz, lnz_ini;

  /** - Set initial parameters */
  current_mass = pin->PBH_evaporation_mass;                                                         // [g]
  pin->PBH_z_evaporation = 0;
  lnz_ini = log(1+pin->z_initial);
  dlnz = lnz_ini/(pin->Nz_PBH-1);
  loop_z = pin->z_initial*1.0001;
  time_prev = 0.;                                                                                   // [s]

  /** - Alloate local variables */
  class_alloc(pvecback_loop,
              pba->bg_size*sizeof(double),
              pin->error_message);

  /** - Alloate variables for PBH mass evolution */
  class_alloc(pin->PBH_table_z,
              pin->Nz_PBH*sizeof(double),
              pin->error_message);
  class_alloc(pin->PBH_table_mass,
              pin->Nz_PBH*sizeof(double),
              pin->error_message);
  class_alloc(pin->PBH_table_mass_dd,
              pin->Nz_PBH*sizeof(double),
              pin->error_message);
  class_alloc(pin->PBH_table_F,
              pin->Nz_PBH*sizeof(double),
              pin->error_message);
  class_alloc(pin->PBH_table_F_dd,
              pin->Nz_PBH*sizeof(double),
              pin->error_message);

  /** - Fill tables with PBH mass evolution */
  /* For the parametrization of F(M) we follow PRD44 (1991) 376 with
   * the additional modification that we dress the "free QCD-particles"
   * (gluons and quarks) with an sigmoid-activation function (in log10-space:
   * Mean at 0.3 GeV and a width of 0.1*"order of magnitude") and the hadrons
   * with (1 - activation) to take the QCD-phase transition into account
   * and to be in agreement with PRD41 (1990) 3052, where the Ansatz is taken
   * that a black hole emmits those particles which appear elementary at the
   * given energy. */
  for(i_step = 0; i_step<pin->Nz_PBH; i_step++) {

    /** - Find value of f(M) */
    current_pbh_temperature = 1.06e13/current_mass;                                                 // [GeV]
    pin->PBH_QCD_activation = 1./(1.+exp(-(log(current_pbh_temperature)-log(0.3))/(log(10.)*0.1))); // [-] see Eq. (4.6) of Stoecker et al. 2018

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
          )*pin->PBH_QCD_activation;
    f_pi = (1.*0.267*exp(-(current_mass*1.350e-1)/(2.66*1.06e13))                                   // pi^0
            +2.*0.267*exp(-(current_mass*1.396e-1)/(2.66*1.06e13))                                  // pi^+
           )*(1-pin->PBH_QCD_activation);
    f_bos = 6.*0.060*exp(-(current_mass*80.39)/(6.04*1.06e13))                                      // W
            +3.*0.060*exp(-(current_mass*91.19)/(6.04*1.06e13))                                     // Z
            +1.*0.267*exp(-(current_mass*125.1)/(2.66*1.06e13));                                   // h
    f = f_EM+f_nu+f_q+f_pi+f_bos;

    /** - Find current time value */
    class_call(background_at_z(pba,
                               loop_z,
                               long_info,
                               inter_normal,
                               &last_index_back_loop,
                               pvecback_loop),
               pba->error_message,
               pin->error_message);
    time_now = pvecback_loop[pba->index_bg_time]/(_c_/_Mpc_over_m_);                                // [s]
    dt = time_now-time_prev;
    time_prev = time_now;

    if (i_step > 0) {
      //TODO :: check this step
      if (current_mass > 0.5*pin->PBH_evaporation_mass){
        current_mass = current_mass-5.34e25*f*pow(current_mass,-2)*dt;                              // [g]
      }
      else {
        if(pin->PBH_z_evaporation == 0){
          pin->PBH_z_evaporation = loop_z;
        }
        current_mass = 0.;
        f = 0.;
      }
    }

    /** - Fill tables */
    pin->PBH_table_z[i_step] = loop_z;
    pin->PBH_table_mass[i_step] = current_mass;                                                     // [g]
    pin->PBH_table_F[i_step] = f;                                                                   // [-]
    loop_z = exp(lnz_ini-dlnz*(i_step+1))-1.;

  }

  /** - Free local variables */
  free(pvecback_loop);

  /** - Spline mass and F(M) evolution in z */
  class_call(array_spline_table_lines(pin->PBH_table_z,
                                      pin->Nz_PBH,
                                      pin->PBH_table_mass,
                                      1,
                                      pin->PBH_table_mass_dd,
                                      _SPLINE_NATURAL_,
                                      pin->error_message),
             pin->error_message,
             pin->error_message);
  class_call(array_spline_table_lines(pin->PBH_table_z,
                                      pin->Nz_PBH,
                                      pin->PBH_table_F,
                                      1,
                                      pin->PBH_table_F_dd,
                                      _SPLINE_NATURAL_,
                                      pin->error_message),
             pin->error_message,
             pin->error_message);

  return _SUCCESS_;
}


/**
 * Calculate injection from PBH evaporation
 * The conventions adopted here are the same as in Stoecker et al. 2018.
 *
 * @param pba            Input: pointer to background structure
 * @param pth            Input: pointer to thermodynamics structure
 * @param z              Input: redshift
 * @param energy_rate    Output: energy density injection rate
 * @return the error status
 */
int injection_rate_PBH_evaporation(struct injection * pin,
                                   double z,
                                   double * energy_rate){

  /** - Define local variables */
  int last_index_back;
  double mass, f;
  double dMdt, f_em;

  /** - Interpolate the PBH mass evolution */
  class_call(array_interpolate_spline(pin->PBH_table_z,
                                      pin->Nz_PBH,
                                      pin->PBH_table_mass,
                                      pin->PBH_table_mass_dd,
                                      1,
                                      z,
                                      &last_index_back,
                                      &mass,                                                        // [g]
                                      1,
                                      pin->error_message),
             pin->error_message,
             pin->error_message);

  class_call(array_interpolate_spline(pin->PBH_table_z,
                                      pin->Nz_PBH,
                                      pin->PBH_table_F,
                                      pin->PBH_table_F_dd,
                                      1,
                                      z,
                                      &last_index_back,
                                      &f,                                                           // [-]
                                      1,
                                      pin->error_message),
             pin->error_message,
             pin->error_message);

  /** - Calculate injection rates */
  if(mass <= 0.0001*pin->PBH_evaporation_mass || f <= 0 || z < pin->PBH_z_evaporation){
    *energy_rate = 0.;                                                                                // [J/(m^3 s)]
  }
  else {
    dMdt=5.34e25*f*pow(mass,-2.);                                                                     // [g/s]
    f_em = 0.55*pin->PBH_QCD_activation+(1-pin->PBH_QCD_activation)*(f-6.*0.147)/f;                   // [-]
    *energy_rate = pin->rho_cdm*pin->PBH_evaporation_fraction*f_em*dMdt/pin->PBH_evaporation_mass;    // [J/(m^3 s)]
  }

  return _SUCCESS_;
}


/**
 * Calculate injection from PBH matter accretion
 *
 * @param pin            Input: pointer to injection structure
 * @param z              Input: redshift
 * @param energy_rate    Output: energy density injection rate
 * @return the error status
 */
int injection_rate_PBH_accretion(struct injection * pin,
                                 double z,
                                 double * energy_rate){

  /** - Define local variables */
  double L_ed, M_ed_dot, M_crit, v_B, v_l, v_eff, r_B, t_B;
  double lambda, M_b_dot;
  double Value_min, Value_med, Value_max, a=0, epsilon_0=0.1, epsilon;
  double L_acc;
  double beta_compton_drag, gamma_cooling, tau_cooling;
  double T_ion, T_s, Y_s, theta_s;
  double lambda_1, lambda_2, lambda_ad, lambda_iso, J;

  /** - Initialize local variables */
  /* Eddington luminosity */
  L_ed = 4.*_PI_*_G_*(pin->PBH_accretion_mass*_Sun_mass_)*_m_p_/_sigma_*_c_;                        // [W]
  M_ed_dot= 10.*L_ed/pow(_c_,2.);                                                                   // [kg/s]
  M_crit = 0.01*M_ed_dot;                                                                           // [kg/s]

  /* Boldi definitions */
  v_B = sqrt((1.+pin->x_e)*(pin->T_b*_k_B_)/(_m_p_*pow(_c_,2.)))*_c_;                               // [m/s]
  if(pin->PBH_accretion_relative_velocities < 0.){
    v_l = 30.*MIN(1.,(1.+z)/1.e3)*1.e3;                                                             // [m/s]
    if(v_B < v_l){
      v_eff = sqrt(v_B*v_l);                                                                        // [m/s]
    }
    else{
      v_eff = v_B;                                                                                  // [m/s]
    }
  }
  else{
    v_l = pin->PBH_accretion_relative_velocities*1.e3;                                              // [m/s]
    v_eff = pow(v_l*v_l+v_B*v_B,0.5);                                                               // [m/s]
  }
  r_B = _G_*(pin->PBH_accretion_mass*_Sun_mass_)/pow(v_eff,2.);                                     // [m]
  t_B = _G_*(pin->PBH_accretion_mass*_Sun_mass_)/pow(v_eff,3.);                                     // [s]

  switch(pin->PBH_accretion_recipe){
    /** - Disk accretion from Poulin et al. 1707.04206 */
    case disk_accretion:

    lambda = pin->PBH_accretion_eigenvalue;                                                         // [-]
    M_b_dot = 4.*_PI_*lambda*(pin->rho_b/pow(_c_,2.))*pow(r_B,2.)*v_eff;                            // [kg/s]

    if(pin->PBH_accretion_ADAF_delta == 1e-3){
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
    else if(pin->PBH_accretion_ADAF_delta == 0.1){
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
    else if (pin->PBH_accretion_ADAF_delta == 0.5){
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
  break;

    /** - Spherical accretion from Ali-Haimoud et al. 1612.05644 */
    case spherical_accretion:

    beta_compton_drag = 4./3.*pin->x_e*_sigma_*pin->rho_g*t_B/_m_p_/_c_;                            // [-] Eq. (7)
    gamma_cooling = 2.*(_m_p_/_m_e_)/(1+pin->x_e)*beta_compton_drag;                                // [-] Eq. (8)
    tau_cooling = 1.5/(5.+pow(gamma_cooling,2./3.));                                                // [-] Eq. (28)

    lambda_ad = 0.25*pow(3./5.,1.5);                                                                // [-] Eq. (15)
    lambda_iso = 0.25*exp(1.5);                                                                     // [-] Eq. (20)
    lambda_1 = lambda_ad+(lambda_iso-lambda_ad)*
                                pow(pow(gamma_cooling,2.)/(88.+pow(gamma_cooling,2.)),0.22);        // [-] Eq. (27)
    lambda_2 = exp(4.5/(3.+pow(beta_compton_drag,0.75)))/(pow(pow(1.+beta_compton_drag,0.5)+1.,2.));// [-] Eq. (32)
    lambda = lambda_1*lambda_2/lambda_iso;                                                          // [-] Eq. (33)

    T_ion = 1.5e4*_eV_over_Kelvin_;                                                                 // [eV] see line below Eq. (39)
    Y_s = pow((1.+pin->x_e)/2.,2./3.*13.6/T_ion)*
                  tau_cooling/4.*pow(1.-5./2.*tau_cooling,1./3.)*(_m_p_/_m_e_);                     // [-] Eq. (51)
    T_s = _m_e_*pow(_c_,2.)/_k_B_*Y_s*pow(1.+Y_s/0.27,-1./3.);                                      // [K] Eqs. (50) and (47)
    theta_s = (T_s*_k_B_)/(_m_e_*pow(_c_,2.));

    M_b_dot = 4.*_PI_*lambda*(pin->rho_b/pow(_c_,2.))*pow(r_B,2.)*v_eff;                            // [kg/s] Eq. (6)

    if(theta_s > 1.){                                                                               // Eq. (55)
      J = 27./(2.*_PI_)*(log(2.*theta_s*exp(-0.577)+0.08)+4./3.);                                   // [-]
    }
    else{
      J = 4./_PI_*sqrt(2./_PI_)*pow(theta_s,-0.5)*(1+5.5*pow(theta_s,1.25));                        // [-]
    }

    L_acc = 1./137.*(T_s*_k_B_)/(_m_p_*pow(_c_,2.))*J*pow(M_b_dot*_c_*_c_,2.)/L_ed;                 // [W] Eq. (57)
    break;
    default:
      class_stop(pin->error_message,"Invalid PBH_accretion_recipe");
  }

  *energy_rate = pin->rho_cdm*pin->PBH_accretion_fraction/
                      (pin->PBH_accretion_mass*_Sun_mass_)*L_acc/pow(_c_,2.);                       // [J/(m^3 s)]

  return _SUCCESS_;
}


/**
 * Read and interpolate the deposition function from external file.
 *
 * @param ppr   Input: pointer to precision structure
 * @param pin   Input/Output: pointer to injection structure
 * @return the error status
 */
int injection_read_feff_from_file(struct precision* ppr,
                                  struct injection* pin,
                                  char* f_eff_file){

  /** - Define local variables */
  FILE * fA;
  char line[_LINE_LENGTH_MAX_];
  char * left;
  int headlines = 0;
  int index_z;

  pin->feff_z_size = 0;

  /** - Read file header */
  /* The file is assumed to contain:
   *    - The number of lines of the file
   *    - The columns ( z, f(z) ) where f(z) represents the "effective" fraction of energy deposited
   *      into the medium  at redshift z, in presence of halo formation. */
  class_open(fA, f_eff_file, "r", pin->error_message);

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
      class_test(sscanf(line,"%d",&(pin->feff_z_size)) != 1,
                 pin->error_message,
                 "could not read the initial integer of number of lines in line %i in file '%s' \n",
                 headlines,f_eff_file);

      /* (z, f, ddf)*/
      class_alloc(pin->feff_table,
                  3*pin->feff_z_size*sizeof(double),
                  pin->error_message);
      break;
    }
  }

  /** - Read file */
  for(index_z=0;index_z<pin->feff_z_size;++index_z){
    /* Read coefficients */
    class_test(fscanf(fA,"%lg %lg",
                      &(pin->feff_table[index_z*3+0]),  // z
                      &(pin->feff_table[index_z*3+1])   // f_eff(z)
                     ) != 2,
               pin->error_message,
               "could not read value of parameters coefficients in line %i in file '%s'\n",
               headlines,f_eff_file);
  }

  fclose(fA);

  /** - Spline file contents */
  /* Spline in one dimension */
  class_call(array_spline(pin->feff_table,
                          3,
                          pin->feff_z_size,
                          0,
                          1,
                          2,
                          _SPLINE_NATURAL_,
                          pin->error_message),
             pin->error_message,
             pin->error_message);

  return _SUCCESS_;
}


/**
 * Read and interpolate the branching ratio from external file, if the function
 * in the file is given with respect to redshift.
 *
 * @param ppr   Input: pointer to precision structure
 * @param pin   Input/Output: pointer to injection structure
 * @return the error status
 */
int injection_read_chi_z_from_file(struct precision* ppr,
                                   struct injection* pin,
                                   char* chi_z_file){

  /** Define local variables */
  FILE * fA;
  char line[_LINE_LENGTH_MAX_];
  char * left;
  int headlines = 0;
  int index_z,index_dep;

  pin->chiz_size = 0;

  /* The file is assumed to contain:
   *    - The number of lines of the file
   *    - The columns (xe , chi_heat, chi_Lya, chi_H, chi_He, chi_lowE) where chi_i represents the
   *      branching ratio at redshift z into different injection/ionization channels i */

  class_open(fA, chi_z_file, "r", pin->error_message);

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
      class_test(sscanf(line,"%d",&(pin->chiz_size)) != 1,
                 pin->error_message,
                 "could not read the initial integer of number of lines in line %i in file '%s' \n",
                 headlines,chi_z_file);

      /* (z, chi_i)*/
      class_alloc(pin->chiz_table,
                  (2*pin->dep_size+1)*pin->chiz_size*sizeof(double),
                  pin->error_message);
      break;
    }
  }

  for(index_z=0;index_z<pin->chiz_size;++index_z){
    /* Read coefficients */
    class_test(fscanf(fA,"%lg %lg %lg %lg %lg %lg",
                      &(pin->chiz_table[index_z*(2*pin->dep_size+1)+0]), //z
                      &(pin->chiz_table[index_z*(2*pin->dep_size+1)+1+pin->index_dep_heat]), //heat
                      &(pin->chiz_table[index_z*(2*pin->dep_size+1)+1+pin->index_dep_lya]), //lya
                      &(pin->chiz_table[index_z*(2*pin->dep_size+1)+1+pin->index_dep_ionH]), //ionH
                      &(pin->chiz_table[index_z*(2*pin->dep_size+1)+1+pin->index_dep_ionHe]), //ionHe
                      &(pin->chiz_table[index_z*(2*pin->dep_size+1)+1+pin->index_dep_lowE])  //lowE
                     )!= 6,
               pin->error_message,
               "could not read value of parameters coefficients in line %i in file '%s'\n",
               index_z+headlines,chi_z_file);
  }

  fclose(fA);

  /* Spline in one dimension */
  for(index_dep=0;index_dep<pin->dep_size;++index_dep){
    class_call(array_spline(pin->chiz_table,
                            2*pin->dep_size+1,
                            pin->chiz_size,
                            0,
                            1+index_dep,
                            1+index_dep+pin->dep_size,
                            _SPLINE_NATURAL_,
                            pin->error_message),
               pin->error_message,
               pin->error_message);
  }

  return _SUCCESS_;
}


/**
 * Read and interpolate the branching ratio from external file, if the function
 * in the file is given with respect to the fraction of free electrons X_e.
 *
 * @param ppr   Input: pointer to precision structure
 * @param pin   Input/Output: pointer to injection structure
 * @return the error status
 */
int injection_read_chi_x_from_file(struct precision* ppr,
                                   struct injection* pin,
                                   char* chi_x_file){

  /** Define local variables */
  FILE * fA;
  char line[_LINE_LENGTH_MAX_];
  char * left;
  int headlines = 0;
  int index_x,index_dep;

  pin->chix_size = 0;

  /** - Read file header */
  /* The file is assumed to contain:
   *    - The number of lines of the file
   *    - The columns (xe , chi_heat, chi_Lya, chi_H, chi_He, chi_lowE) where chi_i represents the
   *      branching ratio at redshift z into different injection/ionization channels i */

  class_open(fA, chi_x_file, "r", pin->error_message);

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
      class_test(sscanf(line,"%d",&(pin->chix_size)) != 1,
                 pin->error_message,
                 "could not read the initial integer of number of lines in line %i in file '%s' \n",
                 headlines, chi_x_file);

      /* (z, chi_i)*/
      class_alloc(pin->chix_table,
                  (2*pin->dep_size+1)*pin->chix_size*sizeof(double),
                  pin->error_message);
      break;
    }
  }

  /** - Read file */
  for(index_x = 0; index_x < pin->chix_size;++index_x){
    /* Read coefficients */
    class_test(fscanf(fA,"%lg %lg %lg %lg %lg %lg",
                      &(pin->chix_table[index_x*(2*pin->dep_size+1)+0]), //x
                      &(pin->chix_table[index_x*(2*pin->dep_size+1)+1+pin->index_dep_heat]), //heat
                      &(pin->chix_table[index_x*(2*pin->dep_size+1)+1+pin->index_dep_lya]), //lya
                      &(pin->chix_table[index_x*(2*pin->dep_size+1)+1+pin->index_dep_ionH]), //ionH
                      &(pin->chix_table[index_x*(2*pin->dep_size+1)+1+pin->index_dep_ionHe]), //ionHe
                      &(pin->chix_table[index_x*(2*pin->dep_size+1)+1+pin->index_dep_lowE])  //lowE
                     )!= 6,
               pin->error_message,
               "could not read value of parameters coefficients in line %i in file '%s'\n",
               index_x+headlines,chi_x_file);
  }

  fclose(fA);

  /** - Spline file contents */
  /* Spline in one dimension */
  for(index_dep=0;index_dep<pin->dep_size;++index_dep){
    class_call(array_spline(pin->chix_table,
                            2*pin->dep_size+1,
                            pin->chix_size,
                            0,
                            1+index_dep,
                            1+index_dep+pin->dep_size,
                            _SPLINE_NATURAL_,
                            pin->error_message),
               pin->error_message,
               pin->error_message);
  }

  return _SUCCESS_;
}




/**
 * Outputs
 */
int injection_output_titles(struct injection * pin, char titles[_MAXTITLESTRINGLENGTH_]){

  class_store_columntitle(titles,"Redshift z",_TRUE_);
  class_store_columntitle(titles,"Heat [-]",_TRUE_);
  class_store_columntitle(titles,"IonH [-]",_TRUE_);
  class_store_columntitle(titles,"IonHe [-]",_TRUE_);
  class_store_columntitle(titles,"Lya [-]",_TRUE_);

  return _SUCCESS_;
}

int injection_output_data(struct injection * pin,
                          int number_of_titles,
                          double * data){
  int storeidx;
  double * dataptr;
  int index_z;

  for (index_z=0; index_z<pin->z_size; index_z++) {
    dataptr = data + index_z*number_of_titles;
    storeidx = 0;
    class_store_double(dataptr, pin->z_table[index_z], _TRUE_, storeidx);
    class_store_double(dataptr, pin->deposition_table[pin->index_dep_heat][index_z], _TRUE_, storeidx);
    class_store_double(dataptr, pin->deposition_table[pin->index_dep_ionH][index_z], _TRUE_, storeidx);
    class_store_double(dataptr, pin->deposition_table[pin->index_dep_ionHe][index_z], _TRUE_, storeidx);
    class_store_double(dataptr, pin->deposition_table[pin->index_dep_lya][index_z], _TRUE_, storeidx);
  }

  return _SUCCESS_;
}
