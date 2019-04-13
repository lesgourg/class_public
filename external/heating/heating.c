/** @file heating.c Documented heating module
 *
 * Initially written by:
 * Matteo Lucca, 27.02.2019
 * Nils Schoeneberg, 27.02.2019
 *
 * The main goal of this module is to calculate the deposited energy in form of heating, ionization
 * and Lyman alpha processes. For more details see the description in the README file.
 *
 */
#include "heating.h"
#include "primordial.h"

//TODO :: disable branching ratios before z > 2000, and replace with only heating

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

  /** Initialize indeces and parameters */
  phe->last_index_chix = 0;
  phe->last_index_z_chi = 0;
  phe->last_index_z_feff = 0;

  /** Import constant quantities from background structure */
  phe->H0 = pba->H0*_c_/_Mpc_over_m_;                                                               // [1/s]
  phe->T_g0 = pba->T_cmb;                                                                           // [K]
  phe->Omega0_b = pba->Omega0_b;                                                                    // [-]
  phe->rho0_cdm = pba->Omega0_cdm*pow(phe->H0,2)*3/8./_PI_/_G_*_c_*_c_;                             // [J/m^3]

  /** Import constant quantities from thermodynamics structure */
  phe->fHe = pth->fHe;
  //phe->N_e0 = 3.*phe->H0*phe->H0*phe->Omega0_b/(8.*_PI_*_G_*_m_H_)*(1.-phe->Y_He);
  phe->N_e0 = pth->n_e;

  /** Check energy injection */ //TODO :: do properly
  phe->has_exotic_injection = phe->annihilation_efficiency!=0 || phe->decay_fraction!=0;

  /** Check energy injection for DM annihilation */
  class_test((phe->annihilation_efficiency<0),
             phe->error_message,
             "annihilation parameter cannot be negative");

  class_test((phe->annihilation_efficiency>1.e-4),
             phe->error_message,
             "annihilation parameter suspiciously large (%e, while typical bounds are in the range of 1e-7 to 1e-6)",phe->annihilation_efficiency);

  class_test((phe->annihilation_variation>0),
             phe->error_message,
             "annihilation variation parameter must be negative (decreasing annihilation rate)");

  class_test((phe->annihilation_z<0),
             phe->error_message,
             "characteristic annihilation redshift cannot be negative");

  class_test((phe->annihilation_zmin<0),
             phe->error_message,
             "characteristic annihilation redshift cannot be negative");

  class_test((phe->annihilation_zmax<0),
             phe->error_message,
             "characteristic annihilation redshift cannot be negative");

  class_test((phe->annihilation_efficiency>0) && (pba->has_cdm==_FALSE_),
             phe->error_message,
             "CDM annihilation effects require the presence of CDM!");

  class_test((phe->annihilation_f_halo<0),
             phe->error_message,
             "Parameter for DM annihilation in halos cannot be negative");

  class_test((phe->annihilation_z_halo<0),
             phe->error_message,
             "Parameter for DM annihilation in halos cannot be negative");

  if (phe->heating_verbose > 0){
    if ((phe->annihilation_efficiency >0) && (pth->reio_parametrization == reio_none) && (ppr->recfast_Heswitch >= 3) && (pth->recombination==recfast))
      printf("Warning: if you have DM annihilation and you use recfast with option recfast_Heswitch >= 3, then the expression for CfHe_t and dy[1] becomes undefined at late times, producing nan's. This is however masked by reionization if you are not in reio_none mode.");
  } //TODO :: check if still occurs !!!

  phe->has_DM_ann = phe->annihilation_efficiency!=0;

  /** Check energy injection for DM deacy */
  class_test((phe->decay_fraction<0),
             phe->error_message,
             "decay parameter cannot be negative");

  class_test((phe->decay_fraction>0)&&(pba->has_cdm==_FALSE_),
             phe->error_message,
             "CDM decay effects require the presence of CDM!");

  phe->has_DM_dec = phe->decay_fraction != 0;

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

  /** Read file for deposition function */
  if(phe->f_eff_type == f_eff_from_file){
    class_call(heating_read_feff_from_file(ppr,phe),
               phe->error_message,
               phe->error_message);
  }

  /** Read file for branching ratios */
  if(phe->chi_type == chi_from_x_file){
    class_call(heating_read_chi_x_from_file(ppr,phe),
               phe->error_message,
               phe->error_message);
  }
  else  if(phe->chi_type == chi_from_z_file){
    class_call(heating_read_chi_z_from_file(ppr,phe),
               phe->error_message,
               phe->error_message);
  }

  /** Define indeces of tables */
  phe->to_store = _FALSE_;
  class_call(heating_indices(pth),
             phe->error_message,
             phe->error_message);

  /** Allocate tables and pvecs */
  class_alloc(phe->injection_table,
              phe->inj_size*sizeof(double*),
              phe->error_message);
  for(index_inj=0; index_inj<phe->inj_size; ++index_inj){
    class_alloc(phe->injection_table[index_inj],
                phe->z_size*sizeof(double),
                phe->error_message);
  }

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

  class_alloc(phe->chi,
              phe->dep_size*sizeof(double),
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

  /** Indices for injection table */
  index_inj = 0;
  class_define_index(phe->index_inj_cool   , _TRUE_          , index_inj, 1);
  class_define_index(phe->index_inj_diss   , _TRUE_          , index_inj, 1);
  class_define_index(phe->index_inj_DM_ann , phe->has_DM_ann , index_inj, 1);
  class_define_index(phe->index_inj_DM_dec , phe->has_DM_dec , index_inj, 1);
  class_define_index(phe->index_inj_tot    , _TRUE_          , index_inj, 1);
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

  free(phe->z_table);

  for(index_inj=0;index_inj<phe->inj_size;++index_inj){
    free(phe->injection_table[index_inj]);
  }
  free(phe->injection_table);

  for(index_dep=0;index_dep<phe->dep_size;++index_dep){
    free(phe->deposition_table[index_dep]);
  }
  free(phe->deposition_table);
  free(phe->photon_dep_table);
  free(phe->chi);

  free(phe->pvecdeposition);

  if(phe->f_eff_type == f_eff_from_file){
    free(phe->feff_table);
  }
  if(phe->chi_type == chi_from_z_file){
    free(phe->chiz_table);
  }
  if(phe->chi_type == chi_from_x_file){
    free(phe->chix_table);
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
  phe->nH = phe->N_e0 * pow(1.+ z , 3);
  phe->heat_capacity = (3./2.)*_k_B_*phe->nH*(1.+phe->f_He+x);
  dEdz_inj = 0.;

  /** Import varying quantities from background structure */
  phe->H = pvecback[pba->index_bg_H]*_c_/_Mpc_over_m_;                                              // [1/s]
  phe->a = pvecback[pba->index_bg_a];                                                               // [-]
  phe->t = pvecback[pba->index_bg_time]/_s_over_Mpc_;                                               // [s]
  phe->rho_cdm = pvecback[pba->index_bg_rho_cdm]*_Jm3_over_Mpc2_;                                   // [J/m^3]
  phe->rho_g = pvecback[pba->index_bg_rho_g]*_Jm3_over_Mpc2_;

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
  if(z > phe->filled_until_z){
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
  if(phe->has_exotic_injection){

    /* Annihilating Dark Matter */
    if(phe->has_DM_ann){
      class_call(heating_rate_DM_annihilation(phe,
                                              z,
                                              &rate),
                 phe->error_message,
                 phe->error_message);
      if(phe->to_store){
        phe->injection_table[phe->index_inj_DM_ann][phe->index_z_store] = rate;
      }
      dEdz += rate;
    }

    /* Decaying Dark Matter */
    if(phe->has_DM_dec){
      class_call(heating_rate_DM_decay(phe,
                                       z,
                                       &rate),
                 phe->error_message,
                 phe->error_message);
      if(phe->to_store){
        phe->injection_table[phe->index_inj_DM_dec][phe->index_z_store] = rate;
      }
      dEdz += rate;
    }
  }

  /** Total energy injection */
  if(phe->to_store){
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

  /** Read the deposition factors for each channel */
  /* Coefficient as revised by Galli et al. 2013 (in fact it is an interpolation
     by Vivian Poulin of columns 1 and 2 in Table V of Galli et al. 2013) */
  /* Read file in ionization fraction */
  if(phe->chi_type == chi_from_x_file){
    for(index_dep=0; index_dep<phe->dep_size; ++index_dep){
      class_call(array_interpolate_spline_transposed(phe->chix_table,
                                                     phe->chix_size,
                                                     2*phe->dep_size+1,
                                                     0,
                                                     index_dep+1,
                                                     index_dep+phe->dep_size+1,
                                                     x,
                                                     &phe->last_index_chix,
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
  /* Old approximation from Chen and Kamionkowski */
  else if(phe->chi_type == chi_from_SSCK){
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
    /*
      // coefficient as revised by Slatyer et al. 2013 (in fact it is a fit by Vivian Poulin of columns 1 and 2 in Table V of Slatyer et al. 2013)
      if (x < 1.){
        phe->chi[phe->index_dep_ionH] = 0.369202*pow(1.-pow(x,0.463929),1.70237);
      }
      else{
        phe->chi[phe->index_dep_ionH] = 0.;
      }
      // coefficient as revised by Slatyer et al. 2013 (in fact it is a fit by Vivian Poulin of columns 1 and 2 in Table V of Slatyer et al. 2013)
      if (x < 1.)
        chi_heat = MIN(0.996857*(1.-pow(1.-pow(x,0.300134),1.51035)),1);
      else
        chi_heat = 1.;
    */

    /*
     * Recfast ::

      double bHe_ion=branching_ratio_ions_Chen(XHeII/fHe);
      double bHI_ion=branching_ratio_ions_Chen(XHII);
      double b_heat =branching_ratio_heat_Chen(XHII, XHeII/fHe, fHe);
      //
      double fheat=2.0/3.0*const_e/const_kB/(1.0+fHe+XHII+XHeII);

      dXHI_1s_dt+=-bHI_ion/13.6/(1.0+fHe)*dE_dt;
      if(XHeII>0.0) dXHeI_1s_dt+=-bHe_ion/24.6*fHe/(1.0+fHe)*dE_dt;
      drho_dt+=fheat*b_heat*dE_dt/Tg;

     * */
  }
  else if(phe->chi_type == chi_full_heating){
    phe->chi[phe->index_dep_heat]  = 1.;
    phe->chi[phe->index_dep_ionH]  = 0.;
    phe->chi[phe->index_dep_ionHe] = 0.;
    phe->chi[phe->index_dep_lya]   = 0.;
    phe->chi[phe->index_dep_lowE]  = 0.;
  }
  else{
    class_stop(phe->error_message,"No valid deposition function has been found found.");
  }

  /** Read the correction factor f_eff */
  /* For the file, read in f_eff from file and multiply */
  if(phe->f_eff_type == f_eff_from_file){
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
  /* For the on the spot, we take the user input */
  else if(phe->f_eff_type == f_eff_on_the_spot){
    phe->f_eff = 1.;
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
  int index_z;
  double tau;
  int last_index_back, last_index_thermo;
  double *pvecback, *pvecthermo;
  double R, dkappa;
  int index_k;
  double dEdt;
  int index_dep;

  /** Allocate backgorund and thermodynamcis vectors */
  last_index_back = 0;
  last_index_thermo = 0;
  class_alloc(pvecback,
              pba->bg_size*sizeof(double),
              phe->error_message);
  class_alloc(pvecthermo,
              pth->tt_size*sizeof(double),
              phe->error_message);

  /** Allocate qunatities from primordial structure */
  phe->k_max = 1.e6;
  phe->k_min = 0.12;
  phe->k_size = 500;        /* Found to be reasonable for the integral of acoustic dissipation */
  class_alloc(phe->k, //TODO :: make these not part of the struct anymore
              phe->k_size*sizeof(double),
              phe->error_message);
  class_alloc(phe->pk_primordial_k,
              phe->k_size*sizeof(double),
              phe->error_message);

  /** Import primordial spectrum */
  for (index_k=0; index_k<phe->k_size; index_k++) {
    phe->k[index_k] = exp(log(phe->k_min)+(log(phe->k_max)-log(phe->k_min))/(phe->k_size)*index_k);
    class_call(primordial_spectrum_at_k(ppm,
                                        ppt->index_md_scalars,
                                        linear,
                                        phe->k[index_k],
                                        &phe->pk_primordial_k[index_k]),
               ppm->error_message,
               phe->error_message);
  }

  /* Loop over z and calculate the heating at each point */
  for(index_z=0; index_z<phe->z_size; ++index_z){

    /** Import quantities from background and thermodynamics structure */
    class_call(background_tau_of_z(pba,
                                   phe->z_table[index_z],
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
    phe->nH = phe->N_e0 * pow(phe->a, -3);

    class_call(thermodynamics_at_z(pba,
                                   pth,
                                   phe->z_table[index_z],
                                   pth->inter_normal,
                                   &last_index_thermo,
                                   pvecback,
                                   pvecthermo),
               pth->error_message,
               phe->error_message);

    dkappa = pvecthermo[pth->index_th_dkappa];                                                      // [1/Mpc]
    phe->dkD_dz = 1./(pvecback[pba->index_bg_H]*dkappa)*(16./15.+pow(R,2.)/(1.+R))/(6.*(1.0+R));    // [Mpc^2]
    phe->kD = 2.*_PI_/pvecthermo[pth->index_th_r_d];                                                // [1/Mpc]
    phe->T_b = pvecthermo[pth->index_th_Tb];                                                        // [K]
    phe->x_e = pvecthermo[pth->index_th_xe];                                                        // [-]

    /** Injected energy that does not need to be deposited (i.e. adiabatic terms) */
    /* First order cooling of photons due to adiabatic interaction with baryons */
    class_call(heating_rate_adiabatic_cooling(phe,
                                              phe->z_table[index_z],
                                              &dEdt),
               phe->error_message,
               phe->error_message);
    phe->photon_dep_table[index_z] += dEdt;

    /* Second order acoustic dissipation of BAO */
    class_call(heating_rate_acoustic_diss(phe,
                                          phe->z_table[index_z],
                                          &dEdt),
               phe->error_message,
               phe->error_message);
    phe->photon_dep_table[index_z] += dEdt;
  }

  /* Free allocated space */
  free(pvecback);
  free(pvecthermo);
  free(phe->k);
  free(phe->pk_primordial_k);

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
                                 struct heating* phe){

  /** Define local variables */
  FILE * fA;
  char line[_LINE_LENGTH_MAX_];
  char * left;
  int headlines;
  int index_z,index_dep;

  phe->chiz_size = 0;

  /* The file is assumed to contain:
   *    - The number of lines of the file
   *    - The columns (xe , chi_heat, chi_Lya, chi_H, chi_He, chi_lowE) where chi_i represents the
   *      branching ratio at redshift z into different heating/ionization channels i */
  class_test(phe->dep_size != 5,
             phe->error_message,
             "Invalid number of heating/ionization channels for chi(z) file");

  if (phe->chi_type == chi_from_z_file) {
    class_open(fA, phe->chi_z_file, "r", phe->error_message);
  }
  else{
    class_stop(phe->error_message,
               "Unknown chi type option");
  }

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
      if (phe->chiz_size == 0) {

        /* Read num_lines, infer size of arrays and allocate them */
        class_test(sscanf(line,"%d",&(phe->chiz_size)) != 1,
                   phe->error_message,
                   "could not read the initial integer of number of lines in line %i in file '%s' \n",
                   headlines,phe->chi_z_file);

        /* (z, chi_i)*/
        class_alloc(phe->chiz_table,
                    (2*phe->dep_size+1)*phe->chiz_size*sizeof(double),
                    phe->error_message);
      }
      else {
        /* Read coefficients */
        class_test(sscanf(line,"%lg %lg %lg %lg %lg %lg",
                          &(phe->chiz_table[index_z*(2*phe->dep_size+1)+0]), //z
                          &(phe->chiz_table[index_z*(2*phe->dep_size+1)+1]), //heat
                          &(phe->chiz_table[index_z*(2*phe->dep_size+1)+2]), //lya
                          &(phe->chiz_table[index_z*(2*phe->dep_size+1)+3]), //ionH
                          &(phe->chiz_table[index_z*(2*phe->dep_size+1)+4]), //ionHe
                          &(phe->chiz_table[index_z*(2*phe->dep_size+1)+5])  //lowE
                         )!= 6,
                   phe->error_message,
                   "could not read value of parameters coefficients in line %i in file '%s'\n",
                   headlines,phe->chi_z_file);
        index_z++;
      }
    }
  }

  if(phe->chi_type == chi_from_z_file){
    fclose(fA);
  }

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
                                 struct heating* phe){

  /** Define local variables */
  FILE * fA;
  char line[_LINE_LENGTH_MAX_];
  char * left;
  int headlines;
  int index_x,index_dep;

  phe->chix_size = 0;

  /* The file is assumed to contain:
   *    - The number of lines of the file
   *    - The columns (xe , chi_heat, chi_Lya, chi_H, chi_He, chi_lowE) where chi_i represents the
   *      branching ratio at redshift z into different heating/ionization channels i */
  class_test(phe->dep_size != 5,
             phe->error_message,
             "Invalid number of heating/ionization channels for chi(x) file");

  class_open(fA, phe->chi_x_file, "r", phe->error_message);


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
      if (phe->chix_size == 0) {

        /* Read num_lines, infer size of arrays and allocate them */
        class_test(sscanf(line,"%d",&(phe->chix_size)) != 1,
                   phe->error_message,
                   "could not read the initial integer of number of lines in line %i in file '%s' \n",
                   headlines,phe->chi_x_file);

        /* (z, chi_i)*/
        class_alloc(phe->chix_table,
                    (2*phe->dep_size+1)*phe->chix_size*sizeof(double),
                    phe->error_message);
      }
      else {
        /* Read coefficients */
        class_test(sscanf(line,"%lg %lg %lg %lg %lg %lg",
                          &(phe->chiz_table[index_x*(2*phe->dep_size+1)+0]), //x
                          &(phe->chiz_table[index_x*(2*phe->dep_size+1)+1]), //heat
                          &(phe->chiz_table[index_x*(2*phe->dep_size+1)+2]), //lya
                          &(phe->chiz_table[index_x*(2*phe->dep_size+1)+3]), //ionH
                          &(phe->chiz_table[index_x*(2*phe->dep_size+1)+4]), //ionHe
                          &(phe->chiz_table[index_x*(2*phe->dep_size+1)+5])  //lowE
                         )!= 6,
                   phe->error_message,
                   "could not read value of parameters coefficients in line %i in file '%s'\n",
                   headlines,phe->chi_x_file);
        index_x++;
      }
    }
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
 * Read and interpolate the deposition function from external file.
 *
 * @param ppr   Input: pointer to precision structure
 * @param phe   Input/Output: pointer to heating structure
 * @return the error status
 */
int heating_read_feff_from_file(struct precision* ppr,
                                struct heating* phe){

  /** Define local variables */
  FILE * fA;
  char line[_LINE_LENGTH_MAX_];
  char * left;
  int headlines;
  int index_z;

  phe->feff_z_size = 0;

  /* The file is assumed to contain:
   *    - The number of lines of the file
   *    - The columns ( z, f(z) ) where f(z) represents the "effective" fraction of energy deposited
   *      into the medium  at redshift z, in presence of halo formation. */
  class_open(fA,phe->f_eff_file, "r",phe->error_message);

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
      if (phe->feff_z_size == 0) {

        /* Read num_lines, infer size of arrays and allocate them */
        class_test(sscanf(line,"%d",&(phe->feff_z_size)) != 1,
                   phe->error_message,
                   "could not read the initial integer of number of lines in line %i in file '%s' \n",
                   headlines,phe->f_eff_file);

        /* (z, f, ddf)*/
        class_alloc(phe->feff_table,
                    3*phe->feff_z_size*sizeof(double),
                    phe->error_message);
      }
      else {
        /* Read coefficients */
        class_test(sscanf(line,"%lg %lg",
                          &(phe->feff_table[index_z*3+0]),
                          &(phe->feff_table[index_z*3+1]))!= 2,
                   phe->error_message,
                   "could not read value of parameters coefficients in line %i in file '%s'\n",
                   headlines,phe->f_eff_file);
        index_z++;
      }
    }
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
  double T_g, R_g;

  T_g = phe->T_g0*(1.+z);

  /** Calculate heating rates */
  R_g = ( 2. * _sigma_/_m_e_/_c_ ) * ( 4./3. * phe->rho_g );
  double heat_capacity = (3./2.)*_k_B_*phe->nH*(1.+phe->fHe+phe->x_e);
  *energy_rate = R_g * phe->x_e / (1.+phe->x_e+phe->fHe) * (phe->T_b - T_g) * heat_capacity ;
  //*energy_rate = -(3./2.)*phe->N_e0*pow(1.+z,3.)*(1.+phe->f_He+phe->x_e)*phe->H*_k_B_*T_g;          // [J/(m^3 s)]

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
  int index_k;
  double * integrand_full, * integrand_approx;
  double dQrho_dz;

  /** a) Calculate full function */
  if (phe->heating_rate_acoustic_diss_approx == _FALSE_){
  }

  /** b) Calculate approximated function */
  if (phe->heating_rate_acoustic_diss_approx == _TRUE_){

    class_alloc(integrand_approx,
                phe->k_size*sizeof(double),
                phe->error_message);

    /* Define integrand for approximated function */
    for (index_k=0; index_k<phe->k_size; index_k++) {
      integrand_approx[index_k] = 4.*0.81*pow(phe->k[index_k],2.)*
                                  phe->pk_primordial_k[index_k]*
                                  exp(-2.*pow(phe->k[index_k]/phe->kD,2.))*
                                  phe->dkD_dz;                                                      // [-]
    }

    /** Calculate heating rates */
    /* Integrate approximate function */
    class_call(simpson_integration(phe->k_size,
                                   integrand_approx,
                                   (log(phe->k_max)-log(phe->k_min))/(phe->k_size),
                                   &dQrho_dz,                                                       // [-]
                                   phe->error_message),
               phe->error_message,
               phe->error_message);

    /* Free space */
    free(integrand_approx);
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
  if (z>phe->annihilation_zmax) {
    annihilation_at_z = phe->annihilation_efficiency*
                        exp(-phe->annihilation_variation*pow(log((phe->annihilation_z+1.)/(phe->annihilation_zmax+1.)),2));
  }
  else if (z>phe->annihilation_zmin) {
    annihilation_at_z = phe->annihilation_efficiency*
                        exp(phe->annihilation_variation*(-pow(log((phe->annihilation_z+1.)/(phe->annihilation_zmax+1.)),2)
                                         +pow(log((z+1.)/(phe->annihilation_zmax+1.)),2)));
  }
  else {
    annihilation_at_z = phe->annihilation_efficiency*
                        exp(phe->annihilation_variation*(-pow(log((phe->annihilation_z+1.)/(phe->annihilation_zmax+1.)),2)
                                         +pow(log((phe->annihilation_zmin+1.)/(phe->annihilation_zmax+1.)),2)));
  }

  /** Calculate boost factor due to annihilation in halos */
  if(phe->annihilation_z_halo > 0.){
    boost_factor = phe->annihilation_f_halo * erfc((1+z)/(1+phe->annihilation_z_halo)) / pow(1.+z,3);
  }
  else{
    boost_factor = 0;
  }

  /** Calculate heating rates */
  *energy_rate = pow(phe->rho_cdm/_c_,2.)*phe->annihilation_efficiency*(1.+boost_factor);           // [J/(m^3 s)]

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

  /** Define local variables */

  /** Calculate heating rates */
  *energy_rate = phe->rho_cdm*phe->decay_fraction*phe->decay_Gamma*
                 exp(-phe->decay_Gamma*phe->t);                                // [J/(m^3 s)]
  printf("%g  %g  %g\n",z, *energy_rate, phe->decay_Gamma*phe->t);

  return _SUCCESS_;
}



