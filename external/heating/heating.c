#include "heating.h"
#include "primordial.h"

/**
 * The main goal of this module is to calculate deposited energy in form of heating, ionization
 * and Lyman alpha processes.
 *
 * To achieve this goal, we first calculate the injected energy rate for different physical
 * effects. There are many possible sources of energy injection/extraction (for more details
 * see e.g. Chluba & Sunyaev 2012 and Chluba 2016), some are present even for the standard
 * cosmological model like
 *    1) Adiabatically cooling electrons and barions as described in Chluba & Sunyaev 2012
 *       (see also Khatri, Sunyaev & Chluba 2012 for useful discussion)
 *    2) Dissipation of acoustic waves as described in Chluba, Khatri & Sunyaev 2012 (see also
 *       Chluba 2013 and Diacoumis & Wong 2017 for useful discussion) with two possible
 *       approximations
 *          a) Eq. 42 from Chluba, Khatri & Sunyaev 2012 (approximated to Y_SZ*S_ac) (TODO)
 *          b) Eq. 45 from Chluba, Khatri & Sunyaev 2012
 *       The user can select the preferred option with 'heating approx'=yes/no (default =yes)
 *    3) Cosmological recombination radiation as described in Chluba & Ali-Haimoud 2016 (TODO)
 * while some other are related to new possible physical processes like
 *    4) Annihilating particles (e.g. dark matter) as described in Chluba 2010 and Chluba
 *       & Sunyaev 2012. (see also Chluba 2013 for useful discussion)
 *    5) Decaying relic particles as described in Chluba 2010 and Chluba & Sunyaev 2012
 *       (see also Chluba 2013 for useful discussion) (TODO!!)
 *    6) Evaporation of primordial black holes as described in Poulin et al. 2017 (see also
 *       Tashiro & Sugiyama 2008, Carr et al. 2010 and Carr et al. 2016 for useful discussions) (TODO)
 *    7) Acctretion of matter into primordial black holes both via
 *          a) Spherical accretion as described in Ali-Haimoud & Kamionkowski 2017 (TODO) and
 *          b) Disk accretion as described in Poulin et al. 2017 (TODO)
 *       (see also Carr et al. 2010 for useful discussion)
 * Note furthermore that the dissipation of acoustic waves is entrinsically a second order process, i.e.
 * it contribution can be derived by expanding at second order the both the dissipation problem and the
 * energy transfer by Compton scattering (see e.g. Chluba, Khatri & Sunyaev 2012). As such, this process
 * is handled as a second order contribution, separated form the oders.
 *
 * Once the rate of energy injection is known, the program evaluates a so-called deposition function
 * which determines the amount of energy effectively deposited into the different forms, i.e. heating,
 * ionization and Lyman alpha.
 *
 * The module is in part based on ExoCLASS (see Stoecker et al. 2018) which has been mainly developed
 * by Vivian Poulin and Patrick Stöcker. The CLASS implementation has been done by Nils Schoeneberg
 * and Matteo Lucca.
 */

#define deposit_on_the_spot          0
#define deposit_feff_from_file       1

#define chi_from_SSCK      0
#define chi_from_x_file    1
#define chi_from_z_file    2

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

  // TODO :: set in input
  phe->deposit_energy_as = 0;
  phe->chi_type = chi_from_SSCK;
  phe->heating_rate_acoustic_diss_approx = _TRUE_;

  /** Initialize indeces and parameters */
  phe->last_index_chix = 0;
  phe->last_index_z_chi = 0;
  phe->last_index_z_feff = 0;
  phe->f_eff = 1.;

  /** Check energy injection */
  phe->has_dcdm = _FALSE_;
  if(_FALSE_){//if (pba->Omega_ini_dcdm!=0 || pba->Omega0_dcdmdr !=0){
    phe->has_dcdm = _TRUE_;
  }

  phe->has_exotic_injection = phe->annihilation_efficiency!=0 || phe->decay!=0;

  /** Check energy injection parameters for DM annihilation */
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

  /** Check energy injection parameters for DM deacy */
  class_test((phe->decay<0),
             phe->error_message,
             "decay parameter cannot be negative");

  class_test((phe->decay>0)&&(pba->has_cdm==_FALSE_),
             phe->error_message,
             "CDM decay effects require the presence of CDM!");

  phe->has_DM_dec = phe->decay != 0;

  /** Define redshift tables */
  phe->z_size = pth->tt_size;
  class_alloc(phe->z_table,
              phe->z_size*sizeof(double),
              phe->error_message);

  memcpy(phe->z_table,
         pth->z_table,
         phe->z_size*sizeof(double));
  phe->tol_z_table = 1e-10;
  phe->filled_until_index_z_inj = phe->z_size-1;
  phe->filled_until_z_inj = phe->z_table[phe->filled_until_index_z_inj];
  phe->filled_until_index_z_dep = phe->z_size-1;
  phe->filled_until_z_dep = phe->z_table[phe->filled_until_index_z_dep];
  phe->last_index_z_inj = 0;
  phe->last_index_z_dep = 0;

  /** Read file for deposition function */
  if(phe->deposit_energy_as == deposit_feff_from_file){
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
  class_alloc(phe->deposition_table,
              phe->z_size*phe->dep_size*sizeof(double),
              phe->error_message);
  class_alloc(phe->pvecdeposition,
              phe->dep_size*sizeof(double),
              phe->error_message);
  class_alloc(phe->chi_table,
              phe->dep_size*sizeof(double),
              phe->error_message);
  class_alloc(phe->injection_table,
              phe->z_size*phe->inj_size*sizeof(double),
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

  /** Indeces for injection table */
  index_inj = 0;
  class_define_index(phe->index_inj_cool   , _TRUE_          , index_inj, 1);
  class_define_index(phe->index_inj_diss   , _TRUE_          , index_inj, 1);
  class_define_index(phe->index_inj_DM_ann , phe->has_DM_ann , index_inj, 1);
  class_define_index(phe->index_inj_DM_dec , phe->has_DM_dec , index_inj, 1);
  class_define_index(phe->index_inj_tot    , _TRUE_          , index_inj, 1);
  phe->inj_size = index_inj;

  /** Indeces for deposition (and chi) table */
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

  free(phe->z_table);
  free(phe->chi_table);

  free(phe->deposition_table);
  free(phe->injection_table);

  free(phe->pvecdeposition);

  if(phe->deposit_energy_as == deposit_feff_from_file){
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
 *
 * @param pba         Input: pointer to background structure
 * @param pth         Input: pointer to thermodynamics structure
 * @param x           Input: freaction of free electrons
 * @param z           Input: redshift
 * @param pvecback    Output: vector of background quantities
 * @return the error status
 */
int heating_at_z(struct background* pba,
                 struct thermo* pth,
                 double x,
                 double z,
                 double Tmat,
                 double* pvecback){

  /** Define local variables */
  struct heating* phe = &(pth->he);
  int index_z, index_dep, iz_store;
  double h,a,b;
  double dEdz_inj;

  /** Redefine input parameters */
  phe->T_b = Tmat;                                                                                  // [K]
  phe->x_e = x;                                                                                     // [-]
  index_z = 0;
  dEdz_inj = 0.0;

  /** Import quantities from background structure */
  phe->T_g0 = pba->T_cmb;                                                                           // [K]
  phe->Omega_ini_dcdm = pba->Omega_ini_dcdm;                                                        // [-]
  phe->Omega0_dcdmdr = pba->Omega0_dcdmdr;                                                          // [-]
  phe->Gamma_dcdm = pba->Gamma_dcdm;                                                                // [1/s]

  phe->H = pvecback[pba->index_bg_H]*_c_/_Mpc_over_m_;                                              // [1/s]
  phe->a = pvecback[pba->index_bg_a];                                                               // [-]
  phe->rho_cdm = pvecback[pba->index_bg_rho_cdm]*_GeVcm3_over_Mpc2_*_eV_*1e9*1.e6;                  // [J/m^3]
  if(phe->has_dcdm){
    phe->rho_dcdm = pvecback[pba->index_bg_rho_dcdm]*_GeVcm3_over_Mpc2_*_eV_*1e9*1.e6;
  }
  else{
    phe->rho_dcdm = 0.0;
  }

  /** Import quantities from thermodynamics structure */
  phe->Y_He = pth->YHe;                                                                             // [-]
  phe->N_e0 = pth->n_e;                                                                             // [1/m^3]

  /** Hunt within the redshift table for the given index of deposition */
  class_call(array_spline_hunt(phe->z_table,
                               phe->z_size,
                               z,
                               &(phe->last_index_z_dep),
                               &h,&a,&b,
                               phe->error_message),
             phe->error_message,
             phe->error_message);

  /** Test if and where the new values should be stored in the injection table */
  /* If this value is important, store it */
  if(phe->to_store){
    /* Calculate where to store the value */
    if(fabs(b-1) < phe->tol_z_table){
      iz_store = phe->last_index_z_dep+1;
    }
    else if(fabs(b) < phe->tol_z_table){
      iz_store = phe->last_index_z_dep;
    }
    /* Could not find a matching index in the z table for this z */
    else{
      class_stop(phe->error_message,
                 "Should store z = %.10e, but it was not in the z table (next lower = %.10e , next higher = %.10e )",
                 phe->z_table[phe->last_index_z_dep],phe->z_table[phe->last_index_z_dep+1]);
    }
  }

  /** Test if the values are already within the table */
  else if(z > phe->filled_until_z_dep){
    /* (Linearly) interpolate within the table */
    for(index_dep=0; index_dep<phe->dep_size; ++index_dep){
      phe->pvecdeposition[index_dep] = phe->deposition_table[phe->last_index_z_dep*phe->inj_size+index_dep] * a
                                       + phe->injection_table[(phe->last_index_z_dep+1)*phe->inj_size+index_dep] * b;
    }

    return _SUCCESS_;
  }

  /** Step 1 - get the injected energy that needs to be deposited */
  class_call(heating_energy_injection_at_z(phe,
                                           z,
                                           &dEdz_inj),
             phe->error_message,
             phe->error_message);

  /** Step 2 - Now deposit the energy we have injected */
  class_call(heating_deposition_function_at_z(phe,
                                              x,
                                              z),
             phe->error_message,
             phe->error_message);

  /** Step 3 - Put result into deposition vector */
  for(index_dep = 0; index_dep < phe->dep_size; ++index_dep){
    phe->pvecdeposition[index_dep] = phe->chi_table[index_dep]*dEdz_inj;
  }

  /** The output is now successfully stored in the deposition table */
  if(phe->to_store){
    for(index_dep = 0; index_dep < phe->dep_size; ++index_dep){
      phe->deposition_table[iz_store*phe->dep_size+index_dep] = phe->pvecdeposition[index_dep];
    }
    class_test(iz_store < phe->filled_until_index_z_dep-1,
               phe->error_message,
               "Skipping too far ahead in z_table. Check that the heating and thermodynamics module agree in their z sampling.");

    phe->filled_until_index_z_dep = iz_store;
    phe->filled_until_z_dep = phe->z_table[iz_store];
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

  /* Hunt within the table for the given index of injection */
  class_call(array_spline_hunt(phe->z_table,
                               phe->z_size,z,
                               &(phe->last_index_z_inj),
                               &h,&a,&b,
                               phe->error_message),
             phe->error_message,
             phe->error_message);

  /** Test if and where the new values should be stored in the injection table */
  /* If this value is important, store it */
  if(phe->to_store){
    /* Calculate where to store the value*/
    if(fabs(b-1) < phe->tol_z_table){
      iz_store = phe->last_index_z_inj+1;
    }
    else if(fabs(b) < phe->tol_z_table){
      iz_store = phe->last_index_z_inj;
    }
    /* Could not find a matching index in the z table for this z */
    else{
      class_stop(phe->error_message,
                 "Should store z = %.10e, but it was not in the z table (next lower = %.10e , next higher = %.10e )",
                 phe->z_table[phe->last_index_z_inj],phe->z_table[phe->last_index_z_inj+1]);
    }
  }

  /** Test if the values are already within the table */
  else if(z > phe->filled_until_z_inj){
    /* (Linearly) interpolate within the table */
    for(index_inj=0; index_inj<phe->inj_size; ++index_inj){
      dEdz += phe->injection_table[phe->last_index_z_inj*phe->inj_size+index_inj] * a
              + phe->injection_table[(phe->last_index_z_inj+1)*phe->inj_size+index_inj] * b;
    }

    *dEdz_inj = dEdz;

    return _SUCCESS_;
  }

  /** Standard energy injection mechanisms */
  class_call(heating_rate_adiabatic_cooling(phe,z,&rate),
             phe->error_message,
             phe->error_message);
   if(phe->to_store){
      phe->injection_table[iz_store*phe->inj_size+phe->index_inj_cool] = rate;
   }
    dEdz += rate;

  /** Exotic energy injection mechanisms */
  if(phe->has_exotic_injection){

    /* Annihilating Dark Matter */
    if(phe->has_DM_ann){
      class_call(heating_rate_DM_annihilation(phe,z,&rate),
                 phe->error_message,
                 phe->error_message);
      if(phe->to_store){
        phe->injection_table[iz_store*phe->inj_size+phe->index_inj_DM_ann] = rate;
      }
      dEdz += rate;
    }

    /* Decaying Dark Matter */
    if(phe->has_DM_dec){
      class_call(heating_rate_DM_decay(phe,z,&rate),
                 phe->error_message,
                 phe->error_message);
      if(phe->to_store){
        phe->injection_table[iz_store*phe->inj_size+phe->index_inj_DM_dec] = rate;
      }
      dEdz += rate;
    }
  }

  if(phe->to_store){
    phe->injection_table[iz_store*phe->inj_size+phe->index_inj_tot] = dEdz;

    class_test(iz_store < phe->filled_until_index_z_inj-1,
               phe->error_message,
               "Skipping too far ahead in z_table. Check that the heating and thermodynamics module agree in their z sampling.");

    phe->filled_until_index_z_inj = iz_store;
    phe->filled_until_z_inj = phe->z_table[iz_store];
  }

  *dEdz_inj = dEdz;

  return _SUCCESS_;
}


/**
 * Calculate deposition function at given redshift.
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
  double f_eff;

  x = 1.0;                            // TODO :: remove

  /** Step 1 - Read the deposition factors for each channel */
  if (x < 1.){                        // TODO :: why is this a good condition ???!?
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
                                                       &phe->chi_table[index_dep],
                                                       phe->error_message),
                   phe->error_message,
                   phe->error_message);
      }
    }
    /* Read file in redshift */
    if(phe->chi_type == chi_from_z_file){
      for(index_dep=0;index_dep<phe->dep_size;++index_dep){
        class_call(array_interpolate_spline_transposed(phe->chiz_table,
                                                       phe->chiz_size,
                                                       2*phe->dep_size+1,
                                                       0,
                                                       index_dep+1,
                                                       index_dep+phe->dep_size+1,
                                                       z,
                                                       &phe->last_index_z_chi,
                                                       &phe->chi_table[index_dep],
                                                       phe->error_message),
                   phe->error_message,
                   phe->error_message);
      }
    }
    /* Old approximation from Chen and Kamionkowski */
    if(phe->chi_type == chi_from_SSCK){
      phe->chi_table[phe->index_dep_heat]  = (1.+2.*x)/3.;
      phe->chi_table[phe->index_dep_ionH]  = (1.-x)/3.;
      phe->chi_table[phe->index_dep_ionHe] = 0.;
      phe->chi_table[phe->index_dep_lya]   = (1.-x)/3.;
      phe->chi_table[phe->index_dep_lowE]  = 0.;
    }
  }
  else{
    phe->chi_table[phe->index_dep_heat]  = 1.;
    phe->chi_table[phe->index_dep_ionH]  = 0.;
    phe->chi_table[phe->index_dep_ionHe] = 0.;
    phe->chi_table[phe->index_dep_lya]   = 0.;
    phe->chi_table[phe->index_dep_lowE]  = 0.;
  }

  /** Step 2 - Read the correction factor f_eff */
  /* For the file, read in f_eff from file and multiply */
  if(phe->deposit_energy_as == deposit_feff_from_file){
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
  else if(phe->deposit_energy_as == deposit_on_the_spot){
    phe->f_eff = 1.;
  }
  /* Otherwise, something must have gone wrong */
  else{
    class_stop(phe->error_message,
               "Unknown energy deposition mechanism");
  }

  /** Step 3 - Multiply both to get the desired result */
  /* Multiply deposition factors with overall correction factor */
  for(index_dep=0; index_dep<phe->dep_size; ++index_dep){
    phe->chi_table[index_dep] *= f_eff;
  }

  return _SUCCESS_;
}


/**
 * Update heating table with second order energy injection mechanisms.
 *
 * @param pba   Input: pointer to background structure
 * @param pth   Input: pointer to thermodynamics structure
 * @param ppt   Input: pointer to perturbation structure
 * @param ppm   Input: pointer to primordial structure
 * @return the error status
 */
int heating_add_second_order(struct background* pba,
                             struct thermo* pth,
                             struct perturbs* ppt,
                             struct primordial* ppm){

  /** Define local variables */
  struct heating* phe = &(pth->he);
  int index_z;
  double tau;
  int last_index_back, last_index_thermo;
  double *pvecback, *pvecthermo;
  int index_k;
  double rate;

  /** Allocate backgorund and thermodynamcis vectors */
  last_index_back = 0;
  last_index_thermo = 0;
  class_alloc(pvecback,
              pba->bg_size*sizeof(double),
              phe->error_message);
  class_alloc(pvecthermo,
              pth->tt_size*sizeof(double),
              phe->error_message);

  /* Loop over z and calculate the heating at each point */
  for(index_z=0; index_z<phe->z_size; ++index_z){

    /** Import quantities from background structure */
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
    phe->rho_g = pvecback[pba->index_bg_rho_g]*_GeVcm3_over_Mpc2_*_eV_*1e9*1e6;                     // [J/m^3]
    phe->R = (3./4.)*pvecback[pba->index_bg_rho_b]/pvecback[pba->index_bg_rho_g];                   // [-]

    /** Import quantities from thermodynamics structure */
    class_call(thermodynamics_at_z(pba,
                                   pth,
                                   phe->z_table[index_z],
                                   pth->inter_normal,
                                   &last_index_thermo,
                                   pvecback,
                                   pvecthermo),
               pth->error_message,
               phe->error_message);

    phe->dkappa = pvecthermo[pth->index_th_dkappa];                                                 // [1/Mpc]
    phe->dkD_dz = 1./(pvecback[pba->index_bg_H]*phe->dkappa)*
                  (16./15.+pow(phe->R,2.)/(1.+phe->R))/(6.*(1.0+phe->R));                           // [Mpc^2]
    phe->kD = 2.*_PI_/pvecthermo[pth->index_th_r_d];                                                // [1/Mpc]

    /** Import quantities from primordial structure */
    phe->k_max = 5.*phe->kD;
    phe->k_min = 0.12;
    phe->k_size = 500;        /* Found to be reasonable for this particular integral */

    class_alloc(phe->k,
                phe->k_size*sizeof(double),
                phe->error_message);

    class_alloc(phe->pk_primordial_k,
                phe->k_size*sizeof(double),
                phe->error_message);

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

    /** Update injection table with second order energy injection mechanisms */
    class_call(heating_rate_acoustic_diss(phe,
                                          phe->z_table[index_z],
                                          &phe->injection_table[index_z*phe->inj_size+phe->index_inj_diss]),
               phe->error_message,
               phe->error_message);

    /* Free allocated space */
    free(phe->k);
    free(phe->pk_primordial_k);

  }

  /* Free allocated space */
  free(pvecback);
  free(pvecthermo);

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
    class_open(fA, ppr->energy_deposition_chi_z_file, "r", phe->error_message);
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
                   "could not read the initial integer of number of lines in line %i in file '%s' \n",headlines,ppr->energy_deposition_feff_file);

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
                   "could not read value of parameters coefficients in line %i in file '%s'\n",headlines,ppr->energy_deposition_chi_z_file);
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

  class_open(fA, ppr->energy_deposition_chi_x_file, "r", phe->error_message);


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
                   "could not read the initial integer of number of lines in line %i in file '%s' \n",headlines,ppr->energy_deposition_feff_file);

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
                   "could not read value of parameters coefficients in line %i in file '%s'\n",headlines,ppr->energy_deposition_chi_x_file);
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
  class_open(fA,ppr->energy_deposition_feff_file, "r",phe->error_message);

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
                   "could not read the initial integer of number of lines in line %i in file '%s' \n",headlines,ppr->energy_deposition_feff_file);

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
                   "could not read value of parameters coefficients in line %i in file '%s'\n",headlines,ppr->energy_deposition_feff_file);
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
  double alpha_h;

  /** Calculate heating rates */
  alpha_h = (3./2.)*phe->N_e0*pow(1.+z,3.)*(1.+phe->Y_He+phe->x_e);                                 // [1/m^3]

  *energy_rate = -phe->H*alpha_h*_k_B_*phe->T_g0*(1.+z);                                            // [J/(m^3 s)]

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
 * In case of non-minimal cosmology, this function determines the energy rate
 * injected in the IGM at a given redshift z (= on-the-spot annihilation) by
 * DM annihilation.
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

  /** Standard formula for annihilating energy injection */
  *energy_rate = phe->rho_cdm*phe->rho_cdm/_c_/_c_*annihilation_at_z*(1.+boost_factor);  // [J/(m^3 s)]

  return _SUCCESS_;
}


/**
 * In case of non-minimal cosmology, this function determines the energy rate
 * injected in the IGM at a given redshift z (= on-the-spot annihilation) by
 * DM decay.
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

  /* Standard formula for decaying dark matter */
  *energy_rate = phe->rho_dcdm*phe->Gamma_dcdm;                                                     // [J/m^3 * ? * Mpc^(-1)]

  return _SUCCESS_;
}



