#include "heating.h"
/** ENERGY INJECTION FUNCTIONS
 *
 * Developed by Vivian Poulin (added functions for energy repartition from DM annihilations or decays and f_eff),
 *              Patrick Stöcker (20.02.17: added external script to calculate the annihilation coefficients on the fly) and
 *              Matteo Lucca (11.02.19: rewrote section in CLASS style)
 *              Nils Schoeneberg (6.03.19: Added struct and module handling)
 */
#define deposit_on_the_spot          0
#define deposit_feff_from_file       1
#define deposit_from_DarkAges        2
#define deposit_analytical_integral  3

#define chi_from_SSCK      0
#define chi_from_x_file    1
#define chi_from_z_file    2
#define chi_from_DarkAges  3

//TODO :: disable branching ratios before z > 2000, and replace with only heating
int heating_init(struct precision * ppr, struct background* pba, struct thermo* pth){

  struct heating* phe = &(pth->he);
  phe->H0 = pba->H0*_c_/_Mpc_over_m_;
  phe->Omega0_cdm = pba->Omega0_cdm;
  phe->rho_crit0 = phe->H0*phe->H0*3/8./_PI_/_G_*_c_*_c_;
  phe->last_index_bg = 0;
  phe->Gamma_dcdm = pba->Gamma_dcdm;
  phe->has_dcdm = _FALSE_; //pba->Omega_ini_dcdm!=0 || pba->Omega0_dcdmdr !=0;
  phe->chi_type = chi_from_GSVI;
  phe->f_eff = 1.; //TODO :: read from user instead
  phe->has_BH_acc = _FALSE_;
  phe->has_BH_evap = _FALSE_;

  phe->last_index_z_feff = 0;
  phe->last_index_chix = 0;
  phe->last_index_chiz = 0;

  phe->deposit_energy_as = 0; //TODO :: set in input
  /** - check energy injection parameters for annihilation */
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

  class_test((phe->decay<0),
             phe->error_message,
             "decay parameter cannot be negative");

  class_test((phe->decay>0)&&(pba->has_cdm==_FALSE_),
             phe->error_message,
             "CDM decay effects require the presence of CDM!");

  //phe->has_exotic_injection = phe->annihilation!=0 || phe->decay!=0 || phe->PBH_accreting_mass!=0 || phe->PBH_evaporating_mass != 0;
  phe->has_exotic_injection = phe->annihilation_efficiency!=0 || phe->decay!=0;
  phe->has_DM_ann = phe->annihilation_efficiency!=0;
  phe->has_DM_dec = phe->decay != 0;

  phe->z_size = pth->tt_size;
  class_alloc(phe->z_table,
              phe->z_size*sizeof(double),
              phe->error_message);

  memcpy(phe->z_table,pth->z_table,phe->z_size*sizeof(double));
  phe->tol_z_table = 1e-10;
  phe->filled_until_index_z_inj = phe->z_size-1;
  phe->filled_until_z_inj = phe->z_table[phe->filled_until_index_z_inj];
  phe->filled_until_index_z_dep = phe->z_size-1;
  phe->filled_until_z_dep = phe->z_table[phe->filled_until_index_z_dep];
  phe->last_index_z_inj = 0;
  phe->last_index_z_dep = 0;

  if(phe->deposit_energy_as == deposit_feff_from_file){
    class_call(heating_read_feff_from_file(ppr,phe),
               phe->error_message,
               phe->error_message);
  }
  if(phe->chi_type == chi_from_DarkAges || phe->chi_type == chi_from_z_file){
    class_call(heating_read_chi_z(ppr,phe),
               phe->error_message,
               phe->error_message);
  }
  if(phe->chi_type == chi_from_x_file){
    class_call(heating_read_chi_x(ppr,phe),
               phe->error_message,
               phe->error_message);
  }

  phe->to_store = _FALSE_;
  class_call(heating_indices(pth),
             phe->error_message,
             phe->error_message);

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

  phe->decay_fraction *= _c_/_Mpc_over_m_; //TODO :: units
  return _SUCCESS_;
}

int heating_indices(struct thermo* pth){

  struct heating* phe = &(pth->he);
  int index_dep,index_inj;


  /* For injection table */
  index_inj = 0;
  class_define_index(phe->index_inj_BAO      ,_TRUE_            , index_inj, 1);
  class_define_index(phe->index_inj_CRR      ,_TRUE_            , index_inj, 1);
  class_define_index(phe->index_inj_DM_ann   ,phe->has_DM_ann   , index_inj, 1);
  class_define_index(phe->index_inj_DM_dec   ,phe->has_DM_dec   , index_inj, 1);
  class_define_index(phe->index_inj_BH_acc   ,phe->has_BH_acc   , index_inj, 1);
  class_define_index(phe->index_inj_BH_evap  ,phe->has_BH_evap  , index_inj, 1);
  class_define_index(phe->index_inj_tot      ,_TRUE_            , index_inj, 1);
  phe->inj_size = index_inj;


  /* For deposition (and chi) table */
  index_dep = 0;
  class_define_index(phe->index_dep_heat  ,_TRUE_, index_dep, 1);
  class_define_index(phe->index_dep_ionH  ,_TRUE_, index_dep, 1);
  class_define_index(phe->index_dep_ionHe ,_TRUE_, index_dep, 1);
  class_define_index(phe->index_dep_lya   ,_TRUE_, index_dep, 1);
  class_define_index(phe->index_dep_lowE  ,_TRUE_, index_dep, 1);

  phe->dep_size = index_dep;

  return _SUCCESS_;
}

/*
 * At some point, distortions.c will call this function,
 * and the acoustic dissipation contributions will be added to the table of heatings
 * */
int heating_add_second_order_terms(struct thermo* pth, struct perturbs* ppt){

  struct heating* phe = &(pth->he);
  //class_define_index(phe->index_inj_BAO,_TRUE_,phe->ht_size,1);

  return _SUCCESS_;
}


int heating_energy_injection_at_z(struct heating* phe, double z, double* dEdz_inj){

  /** Define local variable */
  double dEdz, rate;
  double h,a,b;
  int index_inj, iz_store;

  /* Initialize local variables */
  dEdz = 0.;

  /* Hunt within the table for the given index of injection */
  class_call(array_spline_hunt(phe->z_table,phe->z_size,z,&(phe->last_index_z_inj),&h,&a,&b,phe->error_message),
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
  else if( z > phe->filled_until_z_inj ){
    /* (Linearly) Interpolate within the table */
    for(index_inj=0;index_inj<phe->inj_size;++index_inj){
      dEdz += phe->injection_table[phe->last_index_z_inj*phe->inj_size+index_inj] * a + phe->injection_table[(phe->last_index_z_inj+1)*phe->inj_size+index_inj] * b;
    }
    *dEdz_inj = dEdz;
    return _SUCCESS_;
  }

  /** Non-exotic energy injection mechanisms */

  // dEdz += non-exotic-stuff (TODO :: put here non-exotic stuff)

  /** exotic energy injection mechanisms */
  if(phe->has_exotic_injection){

    /* Annihilating Dark Matter */
    if(phe->has_DM_ann){
      class_call(heating_DM_annihilation(phe,z,&rate),
                 phe->error_message,
                 phe->error_message);
      if(phe->to_store){phe->injection_table[iz_store*phe->inj_size+phe->index_inj_DM_ann] = rate;}
      dEdz += rate;
    }

    /* Decaying Dark Matter */
    if(phe->has_DM_dec){
      class_call(heating_DM_decay(phe,z,&rate),
                 phe->error_message,
                 phe->error_message);
      if(phe->to_store){phe->injection_table[iz_store*phe->inj_size+phe->index_inj_DM_dec] = rate;}
      dEdz += rate;
    }

    /* Decaying Dark Matter */
    if(phe->has_BH_acc){
      //class_call(heating_BH_accretion(phe,z,&rate),
      //           phe->error_message,
      //           phe->error_message);
      if(phe->to_store){phe->injection_table[iz_store*phe->inj_size+phe->index_inj_BH_acc] = rate;}
      dEdz += rate;
    }

    
    /* Decaying Dark Matter */
    if(phe->has_BH_evap){
      //class_call(heating_BH_evaporation(phe,z,&rate),
      //           phe->error_message,
      //           phe->error_message);
      if(phe->to_store){phe->injection_table[iz_store*phe->inj_size+phe->index_inj_BH_evap] = rate;}
      dEdz += rate;
    }

    // dEdz += exotic-stuff (TODO :: put here more exotic stuff)
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




/* Check if table extends to given z
 *  If yes)
 *   Interpolate from table all types that are known
 *   (i.e. including acous. diss. if already added)
 *  If no)
 *   Calculate heating as required
 **/
int heating_at_z(struct background* pba, struct thermo* pth, double x, double z, double* pvecback){

  /** Define local variables */
  double tau;
  struct heating* phe = &(pth->he);
  int index_z, index_dep, iz_store;
  double h,a,b;
  double dEdz_inj;


  index_z = 0;
  dEdz_inj = 0.0;

  // [J/m^3] //TODO :: fix this
  phe->rho_cdm = pvecback[pba->index_bg_rho_cdm]*_GeVcm3_over_Mpc2_*_eV_*1e9*1e6;

  if(phe->has_dcdm){
    phe->rho_dcdm = pvecback[pba->index_bg_rho_dcdm]*_GeVcm3_over_Mpc2_*_eV_*1e9*1e6;
  }
  else{
    phe->rho_dcdm = 0.0;
  }

  phe->t = pvecback[pba->index_bg_time];



  /* Hunt within the table for the given index of deposition */
  class_call(array_spline_hunt(phe->z_table,phe->z_size,z,&(phe->last_index_z_dep),&h,&a,&b,phe->error_message),
             phe->error_message,
             phe->error_message);

  /** Test if and where the new values should be stored in the injection table */
  /* If this value is important, store it */
  if(phe->to_store){
    /* Calculate where to store the value*/
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
  else if( z > phe->filled_until_z_dep ){
    /* (Linearly) Interpolate within the table */
    for(index_dep=0;index_dep<phe->dep_size;++index_dep){
      phe->pvecdeposition[index_dep] = phe->deposition_table[phe->last_index_z_dep*phe->inj_size+index_dep] * a + phe->injection_table[(phe->last_index_z_dep+1)*phe->inj_size+index_dep] * b;
    }
    return _SUCCESS_;
  }

  /** Step 1 - get the injected energy that needs to be deposited */
  /* In the case of the analytical integral, the energy injection is somewhat special */
  if(phe->deposit_energy_as == deposit_analytical_integral){
    class_call(heating_deposit_analytical_integral(pba,pth,z,&dEdz_inj),
               phe->error_message,
               phe->error_message);
  }
  /* Otherwise get the current injection */
  else{
    class_call(heating_energy_injection_at_z(phe,z,&dEdz_inj),
               phe->error_message,
               phe->error_message);
  }

  /** Step 2 - Now deposit the energy we have injected */
  class_call(heating_deposition_function(phe,x,z),
             phe->error_message,
             phe->error_message);

  /** Step 3 - Put result into deposition vector */
  for(index_dep = 0; index_dep < phe->dep_size; ++index_dep){
    phe->pvecdeposition[index_dep] = phe->chi_table[index_dep] * dEdz_inj;
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


int heating_deposition_function(struct heating* phe, double x, double z){

  int index_dep;
  double f_eff;

  f_eff = 1.; //Default value
  //TODO :: x is uninitialized for first point
  x = 1.0; //TODO :: remove
  /** Step 1 - Read the deposition factors for each channel */
  if (x < 1.){ //TODO :: why is this a good condition ???!?

    /* coefficient as revised by Galli et al. 2013 (in fact it is an interpolation by Vivian Poulin of columns 1 and 2 in Table V of Galli et al. 2013) */
    /* Read file in ionization fraction */
    if(phe->chi_type == chi_from_x_file){
      for(index_dep=0;index_dep<phe->dep_size;++index_dep){
        class_call(array_interpolate_spline_transposed(phe->chix_table,phe->chix_size,
                                                       2*phe->dep_size+1,0,index_dep+1,index_dep+phe->dep_size+1,
                                                       x,&(phe->last_index_chix),phe->chi_table[index_dep],phe->error_message);
      }
    }
    /* Read file in redshift */
    if(phe->chi_type == chi_from_DarkAges || phe->chi_type == chi_from_z_file){
      for(index_dep=0;index_dep<phe->dep_size;++index_dep){
        class_call(array_interpolate_spline_transposed(phe->chiz_table,phe->chiz_size,
                                                       2*phe->dep_size+1,0,index_dep+1,index_dep+phe->dep_size+1,
                                                       z,&(phe->last_index_chiz),phe->chi_table[index_dep],phe->error_message);
      }
    }
    /* old approximation from Chen and Kamionkowski */
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

  /* For the analytical integral, we already integrated the deposition */
  if(phe->deposit_energy_as == deposit_analytical_integral){
    f_eff = 1.;
  }
  /* For the file, read in f_eff from file and multiply */
  else if(phe->deposit_energy_as == deposit_feff_from_file){
    class_call(array_interpolate_spline_transposed(phe->feff_table,phe->feff_z_size,
                                                   3,0,1,2,z,&(phe->last_index_z_feff),
                                                   &(f_eff),phe->error_message),
           phe->error_message,
           phe->error_message);
  }
  /* For the DarkAges, the chi already contain everything */
  else if(phe->deposit_energy_as == deposit_from_DarkAges){
    f_eff = 1.;
  }
  /* For the on the spot, we take the user input */
  else if(phe->deposit_energy_as == deposit_on_the_spot){
    f_eff = phe->f_eff;
  }
  /* Otherwise, something must have gone wrong */
  else{
    class_stop(phe->error_message,
               "Unknown energy deposition mechanism");
  }




  /** Step 3 - Multiply both to get the desired result */

  /* Multiply deposition factors with overall correction factor */
  for(index_dep=0;index_dep<phe->dep_size;++index_dep){
    phe->chi_table[index_dep] *= f_eff;
  }

  return _SUCCESS_;
}

int heating_free(struct thermo* pth){

  struct heating* phe = &(pth->he);

  free(phe->z_table);
  free(phe->chi_table);

  free(phe->deposition_table);
  free(phe->injection_table);

  free(phe->pvecdeposition);

  if(phe->deposit_energy_as == deposit_feff_from_file){
    free(phe->feff_table);
  }
  if(phe->chi_type == chi_from_DarkAges || phe->chi_type == chi_from_z_file){
    free(phe->chiz_table);
  }
  if(phe->chi_type == chi_from_x_file){
    free(phe->chix_table);
  }

  return _SUCCESS_;
}








int heating_DM_annihilation(struct heating * phe,
                            double z,
                            double * energy_rate){
  double boost_factor;

  /* Calculate boost factor due to annihilation in halos */
  if(phe->annihilation_z_halo > 0.){
    boost_factor = phe->annihilation_f_halo * erfc((1+z)/(1+phe->annihilation_z_halo)) / pow(1.+z,3);
  }
  else{
    boost_factor = 0;
  }

  /* Standard formula for annihilating energy injection */
  *energy_rate = phe->rho_cdm*phe->rho_cdm/_c_/_c_ * phe->annihilation_efficiency * (1.+boost_factor);  // [J/(m^3 s)]

  return _SUCCESS_;
}

int heating_DM_decay(struct heating * phe,
                     double z,
                     double * energy_rate){
  double rho_cdm_today, rho_dcdm,decay_factor;

  /* Define energy density of decaying DM, rho_dcdm */
  if(phe->has_dcdm){
    rho_dcdm = phe->rho_dcdm;
  }
  //TODO :: why does this make sense?
  else{
    /* If Omega_dcdm is not given, rho_dcdm = rho_cdm*decay factor */
    if(phe->has_on_the_spot == _FALSE_){
      /* The exponential decay factor is already in the deposition functions */
      decay_factor=1;
    }
    else{
      /* Add the exponential decay factor*/
      decay_factor = exp(-phe->Gamma_dcdm*phe->t);
    }

    /* Surviving rho dcdm*/
    rho_dcdm = phe->rho_cdm*decay_factor;
  }

  /* Standard formula for decaying dark matter*/
  *energy_rate = rho_dcdm*phe->decay_fraction*phe->Gamma_dcdm; //[J/m^3 * ? * Mpc^(-1)]
  //TODO :: figure out units of everything

  return _SUCCESS_;
}


int heating_read_feff_from_file(struct precision* ppr, struct heating* phe){

  /** Define local variables */
  FILE * fA;
  char line[_LINE_LENGTH_MAX_];
  char * left;
  int headlines;
  int index_z;

  phe->feff_z_size = 0;

  /* *
   * The file is assumed to contain:
   * - The number of lines of the file
   * - The columns ( z, f(z) ) where f(z) represents the "effective" fraction of energy deposited into the medium 
   *   at redshift z, in presence of halo formation.
   * */
  class_open(fA,ppr->energy_deposition_feff_file, "r",phe->error_message);

  while (fgets(line,_LINE_LENGTH_MAX_-1,fA) != NULL) {
    headlines++;

    /* eliminate blank spaces at beginning of line */
    left=line;
    while (left[0]==' ') {
      left++;
    }

    /* check that the line is neither blank nor a comment. In ASCII, left[0]>39 means that first non-blank charachter might
       be the beginning of some data (it is not a newline, a #, a %, etc.) */
    if (left[0] > 39) {

      /* if the line contains data, we must interprete it. If num_lines == 0 , the current line must contain
         its value. Otherwise, it must contain (xe , chi_heat, chi_Lya, chi_H, chi_He, chi_lowE). */
      if (phe->feff_z_size == 0) {

        /* read num_lines, infer size of arrays and allocate them */
        class_test(sscanf(line,"%d",&(phe->feff_z_size)) != 1,
                   phe->error_message,
                   "could not read the initial integer of number of lines in line %i in file '%s' \n",headlines,ppr->energy_deposition_feff_file);

        /* (z, f, ddf)*/
        class_alloc(phe->feff_table,
                    3*phe->feff_z_size*sizeof(double),
                    phe->error_message);
      }
      else {
        /* read coefficients */
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

  /* spline in one dimension */
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


#pragma INCOMPLETE

int heating_read_chi_z(struct precision* ppr, struct heating* phe){

  /** Define local variables */
  FILE * fA;
  char line[_LINE_LENGTH_MAX_];
  char * left;
  int headlines;
  int index_z,index_dep;

  /* variables related to the use of DarkAges calculate the annihilation coefficients */
  char command_with_arguments[2*_ARGUMENT_LENGTH_MAX_];
  int status;

  phe->chiz_size = 0;

  /* *
   * The file is assumed to contain:
   * - The number of lines of the file
   * - The columns (xe , chi_heat, chi_Lya, chi_H, chi_He, chi_lowE) where chi_i represents the branching ratio
   *   at redshift z into different heating/ionization channels i
   * */
  class_test(phe->dep_size != 5,
             phe->error_message,
             "Invalid number of heating/ionization channels for chi(z) file");

  if (phe->chi_type == chi_from_z_file) {
    class_open(fA, ppr->energy_deposition_chi_z_file, "r", phe->error_message);
  } 
  else if(phe->chi_type == chi_from_DarkAges){
    /* Write the command */
    sprintf(command_with_arguments, "%s", phe->command_DarkAges);

    if (phe->heating_verbose > 0) {
      printf(" -> Running DarkAges: %s\n", command_with_arguments);
    }

    /* Launch the process and retrieve the output */
    fflush(fA);
    fA = popen(command_with_arguments, "r");
    class_test(fA == NULL, phe->error_message, "The program failed to set the environment for the external command.");
  }
  else{
    class_stop(phe->error_message,
               "Unknown chi type option");
  }

  while (fgets(line,_LINE_LENGTH_MAX_-1,fA) != NULL) {
    headlines++;

    /* eliminate blank spaces at beginning of line */
    left=line;
    while (left[0]==' ') {
      left++;
    }

    /* check that the line is neither blank nor a comment. In ASCII, left[0]>39 means that first non-blank charachter might
       be the beginning of some data (it is not a newline, a #, a %, etc.) */
    if (left[0] > 39) {

      /* if the line contains data, we must interprete it. If num_lines == 0 , the current line must contain
         its value. Otherwise, it must contain (xe , chi_heat, chi_Lya, chi_H, chi_He, chi_lowE). */
      if (phe->chiz_size == 0) {

        /* read num_lines, infer size of arrays and allocate them */
        class_test(sscanf(line,"%d",&(phe->chiz_size)) != 1,
                   phe->error_message,
                   "could not read the initial integer of number of lines in line %i in file '%s' \n",headlines,ppr->energy_deposition_feff_file);

        /* (z, chi_i)*/
        class_alloc(phe->chiz_table,
                    (2*phe->dep_size+1)*phe->chiz_size*sizeof(double),
                    phe->error_message);
      }
      else {
        /* read coefficients */
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

  /* spline in one dimension */
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

int heating_read_chi_x(struct precision* ppr, struct heating* phe){

  /** Define local variables */
  FILE * fA;
  char line[_LINE_LENGTH_MAX_];
  char * left;
  int headlines;
  int index_x,index_dep;

  phe->chix_size = 0;

  /* *
   * The file is assumed to contain:
   * - The number of lines of the file
   * - The columns (xe , chi_heat, chi_Lya, chi_H, chi_He, chi_lowE) where chi_i represents the branching ratio
   *   at redshift z into different heating/ionization channels i
   * */
  class_test(phe->dep_size != 5,
             phe->error_message,
             "Invalid number of heating/ionization channels for chi(x) file");

  class_open(fA, ppr->energy_deposition_chi_x_file, "r", phe->error_message);


  while (fgets(line,_LINE_LENGTH_MAX_-1,fA) != NULL) {
    headlines++;

    /* eliminate blank spaces at beginning of line */
    left=line;
    while (left[0]==' ') {
      left++;
    }

    /* check that the line is neither blank nor a comment. In ASCII, left[0]>39 means that first non-blank charachter might
       be the beginning of some data (it is not a newline, a #, a %, etc.) */
    if (left[0] > 39) {

      /* if the line contains data, we must interprete it. If num_lines == 0 , the current line must contain
         its value. Otherwise, it must contain (xe , chi_heat, chi_Lya, chi_H, chi_He, chi_lowE). */
      if (phe->chix_size == 0) {

        /* read num_lines, infer size of arrays and allocate them */
        class_test(sscanf(line,"%d",&(phe->chix_size)) != 1,
                   phe->error_message,
                   "could not read the initial integer of number of lines in line %i in file '%s' \n",headlines,ppr->energy_deposition_feff_file);

        /* (z, chi_i)*/
        class_alloc(phe->chix_table,
                    (2*phe->dep_size+1)*phe->chix_size*sizeof(double),
                    phe->error_message);
      }
      else {
        /* read coefficients */
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

  /* spline in one dimension */
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

#pragma NOT_YET_PROPERLY_DONE

int heating_deposit_analytical_integral(struct background* pba, struct thermo* pth, double z, double* energy_rate){

  struct heating* phe = &(pth->he);
  double zp,dz;
  double integrand,first_integrand;
  double factor,result;
  double nH0;
  double onthespot;
  double exponent_z,exponent_zp;

  nH0 = 3.*pba->H0*pba->H0*pba->Omega0_b/(8.*_PI_*_G_*_m_H_)*(1.-pth->YHe); // number of hydrogen nuclei today in m^-3

  /* Value from Poulin et al. 1508.01370 */
  /* 
  factor = c sigma_T n_H(0) / (H(0) \sqrt(Omega_m)) (dimensionless)
  factor = _sigma_ * nH0 / pba->H0 * _Mpc_over_m_ / sqrt(pba->Omega0_b+pba->Omega0_cdm);
  exponent_z = 8;
  exponent_zp = 7.5;
  */

  /* Value from Ali-Haimoud & Kamionkowski 1612.05644 */
  factor = 0.1*_sigma_ * nH0 / pba->H0 * _Mpc_over_m_ / sqrt(pba->Omega0_b+pba->Omega0_cdm);
  exponent_z = 7;
  exponent_zp = 6.5;

  /* integral over z'(=zp) with step dz */
  dz=1.;

  /* first point in trapezoidal integral */
  zp = z;
  class_call(heating_energy_injection_at_z(phe,z,&onthespot),
             phe->error_message,
             phe->error_message);
  first_integrand = factor*pow(1+z,exponent_z)/pow(1+zp,exponent_zp)*exp(2./3.*factor*(pow(1+z,1.5)-pow(1+zp,1.5)))*onthespot; 
  // beware: versions before 2.4.3, there were rwrong exponents: 6 and 5.5 instead of 7 and 7.5
  result = 0.5*dz*first_integrand;

  /* other points in trapezoidal integral */
  do{
    zp += dz;
    class_call(heating_energy_injection_at_z(phe,z,&onthespot),
               phe->error_message,
               phe->error_message);
    integrand = factor*pow(1+z,exponent_z)/pow(1+zp,exponent_zp)*exp(2./3.*factor*(pow(1+z,1.5)-pow(1+zp,1.5)))*onthespot; 
    // beware: versions before 2.4.3, there were rwrong exponents: 6 and 5.5 instead of 7 and 7.5
    result += dz*integrand;

  } while (integrand/first_integrand > 0.02);
  if(result < 1e-100) result=0.;

  return _SUCCESS_;
}




#ifdef NOTDEFINED_DEFINITELY_LUL

/**
 * In case of non-minimal cosmology, this function determines time evolution of the primordial black hole (PBH) mass.
 *
 * @param ppr            Input: pointer to precision structure
 * @param pba            Input: pointer to background structure
 * @param pth            Input: pointer to thermodynamics structure
 * @param z              Input: redshift
 * @param energy_rate    Output: energy density injection rate
 * @param error_message  Output: error message
 * @return the error status
 */
int PBH_evaporating_mass_time_evolution(struct precision * ppr,
                                        struct background * pba,
                                        struct thermo * pth,
                                        ErrorMsg error_message){
  double loop_z, loop_tau, current_mass, time_now, time_prev, dt, dz, f;
  double QCD_activation, current_pbh_temperature;
  double * pvecback_loop;
  int  i_step, last_index_back_loop;

  time_prev = 0.;
  pth->PBH_z_evaporation = 0;
  pth->PBH_table_size = ppr->recfast_Nz0;
  dz = ppr->recfast_z_initial/(pth->PBH_table_size);
  loop_z = ppr->recfast_z_initial-dz;
  current_mass = pth->PBH_evaporating_mass;

  class_alloc(pvecback_loop,pba->bg_size*sizeof(double),pba->error_message);
  class_alloc(pth->PBH_table_z,pth->PBH_table_size*sizeof(double),error_message);
  class_alloc(pth->PBH_table_mass,pth->PBH_table_size*sizeof(double),error_message);
  class_alloc(pth->PBH_table_mass_dd,pth->PBH_table_size*sizeof(double),error_message);
  class_alloc(pth->PBH_table_F,pth->PBH_table_size*sizeof(double),error_message);
  class_alloc(pth->PBH_table_F_dd,pth->PBH_table_size*sizeof(double),error_message);

  /* For the parametrization of F(M) we follow PRD44 (1991) 376 with the additional modification that we dress the "free 
     QCD-particles" (gluons and quarks) with an sigmoid-activation function (in log10-space: Mean at 0.3 GeV and a with 
     of 0.1*"order of magnitude") and the hadrons with 1 - activation to take the QCD-phase transition into account
     and to be in agreement with PRD41 (1990) 3052, where the Ansatz is taken that a black hole emmits those particles which appear
     elementary at the given energy. The order of the particles in the following definition of f: photon, neutrino, electron, muon, 
     tau, up, down, charm, strange, top, bottom, W, Z, gluon, Higgs, neutral Pion and charged pion  */
  for(i_step = 0; i_step < pth->PBH_table_size; i_step++) {
    current_pbh_temperature = 1.06e13/current_mass;
    QCD_activation = 1./(1.-exp(-(log(current_pbh_temperature)-log(0.3))/(log(10.)*0.1)));
    f = 2.*0.060
        +6.*0.147
        +4.*0.142*exp(-(current_mass*5.11e-4)/(4.53*1.06e13))
        +4.*0.142*exp(-(current_mass*0.1037)/(4.53*1.06e13))
        +4.*0.142*exp(-(current_mass*1.777)/(4.53*1.06e13))
        +12.*0.142*exp(-(current_mass*2.2e-3)/(4.53*1.06e13))*QCD_activation
        +12.*0.142*exp(-(current_mass*4.7e-3)/(4.53*1.06e13))*QCD_activation
        +12.*0.142*exp(-(current_mass*1.82)/(4.53*1.06e13))*QCD_activation
        +12.*0.142*exp(-(current_mass*9.6e-2)/(4.53*1.06e13))*QCD_activation
        +12.*0.142*exp(-(current_mass*173.1)/(4.53*1.06e13))*QCD_activation
        +12.*0.142*exp(-(current_mass*4.18)/(4.53*1.06e13))*QCD_activation
        +6.*0.060*exp(-(current_mass*80.39)/(6.04*1.06e13))
        +3.*0.060*exp(-(current_mass*91.19)/(6.04*1.06e13))
        +16.*0.060*exp(-(current_mass*6e-1)/(6.04*1.06e13))*QCD_activation
        +1.*0.267*exp(-(current_mass*125.06)/(2.66*1.06e13))
        +1.*0.267*exp(-(current_mass*1.350e-1)/(2.66*1.06e13))*(1-QCD_activation)
        +2.*0.267*exp(-(current_mass*1.396e-1)/(2.66*1.06e13))*(1-QCD_activation);

    class_call(background_tau_of_z(pba,
                                   loop_z,
                                   &loop_tau),
               pba->error_message,
               ppr->error_message);
    class_call(background_at_tau(pba,
                                 loop_tau,
                                 pba->long_info,
                                 pba->inter_normal,
                                 &last_index_back_loop,
                                 pvecback_loop),
               pba->error_message,
               ppr->error_message);
    time_now = pvecback_loop[pba->index_bg_time]/(_c_/_Mpc_over_m_);
    dt = time_now-time_prev;
    time_prev = time_now;

    if (i_step > 0) {
      if (current_mass > 0.5*pth->PBH_evaporating_mass) {
        current_mass = current_mass-5.34e-5*f*pow(current_mass/1e10,-2)*1e10*dt;
      }
      else {
        if(pth->PBH_z_evaporation == 0) pth->PBH_z_evaporation = loop_z;
        current_mass = 0.;
        f = 0.;
      }
    }

    pth->PBH_table_z[i_step] = loop_z;
    pth->PBH_table_mass[i_step] = current_mass;
    pth->PBH_table_F[i_step] = f;
    loop_z = MAX(0,loop_z-dz);

  }

  free(pvecback_loop);

  class_call(array_spline_table_lines(pth->PBH_table_z,
                                      pth->PBH_table_size,
                                      pth->PBH_table_mass,
                                      1,
                                      pth->PBH_table_mass_dd,
                                      _SPLINE_NATURAL_,
                                      error_message),
             pth->error_message,
             pth->error_message);
  class_call(array_spline_table_lines(pth->PBH_table_z,
                                      pth->PBH_table_size,
                                      pth->PBH_table_F,
                                      1,
                                      pth->PBH_table_F_dd,
                                      _SPLINE_NATURAL_,
                                      error_message),
             error_message,
             error_message);

}

/**
 * In case of non-minimal cosmology, this function determines the energy rate injected in the IGM at a given redshift z (= on-the-spot
 * annihilation) by evaporating primordial black holes.
 *
 * @param ppr            Input: pointer to precision structure
 * @param pba            Input: pointer to background structure
 * @param pth            Input: pointer to thermodynamics structure
 * @param z              Input: redshift
 * @param energy_rate    Output: energy density injection rate
 * @param error_message  Output: error message
 * @return the error status
 */
int heating_evaporating_pbh_energy_injection(struct precision * ppr,
                                                    struct background * pba,
                                                    struct thermo * pth,
                                                    double z,
                                                    double * energy_rate,
                                                    ErrorMsg error_message){

  double rho_cdm_today;
  //double tau;
  int last_index_back;
  double f, f_neutrinos, em_branching, pbh_mass;
  double dMdt;

  /* Calculate the PBH-mass evolution at first call of the function */
  if ((pth->PBH_table_is_initialized) == _FALSE_) {
    pth->PBH_table_is_initialized = _TRUE_;
    PBH_evaporating_mass_time_evolution(ppr,pba,pth,error_message);
  }

  rho_cdm_today = pow(pba->H0*_c_/_Mpc_over_m_,2)*3/8./_PI_/_G_*(pba->Omega0_cdm)*_c_*_c_;         // [J/m^3]

  class_test(pth->PBH_table_is_initialized == _FALSE_, error_message, "The PBH table is not initialized");
  class_call(array_interpolate_spline(pth->PBH_table_z,
				      pth->PBH_table_size,
				      pth->PBH_table_mass,
				      pth->PBH_table_mass_dd,
				      1,
				      z,
				      &last_index_back,
				      &(pbh_mass),
				      1,
				      error_message),
	     error_message,
	     error_message);
  class_call(array_interpolate_spline(pth->PBH_table_z,
				      pth->PBH_table_size,
				      pth->PBH_table_F,
				      pth->PBH_table_F_dd,
				      1,
				      z,
				      &last_index_back,
				      &f,
				      1,
				      error_message),
	     error_message,
	     error_message);

  f_neutrinos = 6*0.147;
  em_branching = 1.; // Currently incoporated in the computation of the f(z) functions.
  if(pbh_mass <= 0.0001*pth->PBH_evaporating_mass || f <= 0 || isnan(pbh_mass)==1 || isnan(f)==1 || z < pth->PBH_z_evaporation){
    pbh_mass = 0;
    dMdt = 0;
    f = 0.;
  }
  else {
    dMdt=5.34e-5*f*pow(pbh_mass/1e10,-2)*1e10;
  }

  *energy_rate = rho_cdm_today*pow((1+z),3)*pth->PBH_fraction/pth->PBH_evaporating_mass*em_branching*(dMdt);

  if(isnan(*energy_rate)==1 || *energy_rate < 0){
    *energy_rate=0.;
  }

}


/**
 * In case of non-minimal cosmology, this function determines the energy rate injected in the IGM at a given redshift z (= on-the-spot
 * annihilation) by accetion of matter into primordial black holes.
 *
 * @param ppr            Input: pointer to precision structure
 * @param pba            Input: pointer to background structure
 * @param pth            Input: pointer to thermodynamics structure
 * @param z              Input: redshift
 * @param energy_rate    Output: energy density injection rate
 * @param error_message  Output: error message
 * @return the error status
 */
int heating_accreting_pbh_energy_injection(struct precision * ppr,
                                                  struct background * pba,
                                                  struct thermo * pth,
                                                  double z,
                                                  double * energy_rate,
                                                  ErrorMsg error_message){

  double rho_cdm_today;
  double tau;
  int last_index_back;
  double * pvecback;
  //Parameters related to PBH
  double c_s, v_eff,v_eff_2,v_l, r_B,x_e,beta,beta_eff,beta_hat,x_cr,lambda,n_gas,M_b_dot,M_sun,M_ed_dot,epsilon,L_acc,Integrale,Normalization;
  double m_H, m_dot, m_dot_2, L_acc_2,L_ed,l,l2,M_crit;
  double rho, m_p = 938, m_e = 0.511, T_infinity = 0, rho_infinity = 0, x_e_infinity = 0, P_infinity = 0, rho_cmb = 0, t_B = 0, v_B = 0;
  double lambda_1,lambda_2,lambda_ad,lambda_iso,gamma_cooling,beta_compton_drag, T_s, T_ion, Y_s, J,tau_cooling;
  double Value_min, Value_med, Value_max, a=0, epsilon_0=0.1;
  rho_cdm_today = pow(pba->H0*_c_/_Mpc_over_m_,2)*3/8./_PI_/_G_*(pba->Omega0_cdm)*_c_*_c_;         // [J/m^3]

  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);
  class_call(background_tau_of_z(pba,
                                 z,
                                 &tau),
             pba->error_message,
             pth->error_message);
  class_call(background_at_tau(pba,
                               tau,
                               pba->long_info,
                               pba->inter_normal,
                               &last_index_back,
                               pvecback),
             pba->error_message,
             pth->error_message);

  c_s = 5.7e3*pow(pth->Tm_tmp/2730,0.5);          // [m]
  M_sun = 2e30;                                   // [kg]
  n_gas = 200*1e6*pow((1+z)/1000,3);              // [1/m^3]
  m_H= 1.67e-27;                                  // [kg]
  x_e = pth->xe_tmp;
  T_infinity = pth->Tm_tmp*_eV_over_Kelvin_*1e-6; // [MeV]

  /** Disk accretion from Poulin et al. 1707.04206 */

  if(pth->PBH_accretion_recipe == disk_accretion){
    L_ed = 4*_PI_*_G_*pth->PBH_accreting_mass*M_sun*m_p*1e6/_eV_over_joules_/(_sigma_*_c_);
    M_ed_dot= 10*L_ed/(_c_*_c_);
    M_crit = 0.01*M_ed_dot;
    v_B = sqrt((1+x_e)*T_infinity/m_p)*_c_;

    if(pth->PBH_relative_velocities < 0.){
      v_l = 30*MIN(1,z/1000)*1e3; // [m/s]
      if(v_B < v_l) v_eff = sqrt(v_B*v_l);
      else v_eff = v_B;
    }
    else{
      v_l = pth->PBH_relative_velocities*1e3; // [m/s]
      v_eff = pow(v_l*v_l+v_B*v_B,0.5);
    }

    lambda = pth->PBH_accretion_eigenvalue;
    rho = pvecback[pba->index_bg_rho_b]/pow(_Mpc_over_m_,2)*3/8./_PI_/_G_*_c_*_c_; // [kg/m^3]
    M_b_dot = 4*_PI_*lambda*pow(_G_*pth->PBH_accreting_mass*M_sun,2)*rho*pow(v_eff,-3.);

    if(pth->PBH_ADAF_delta == 1e-3){
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
      epsilon = epsilon_0 * pow(M_b_dot / M_crit,a);
    }
    else if(pth->PBH_ADAF_delta == 0.1){
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
      epsilon = epsilon_0 * pow(M_b_dot / M_crit,a);
    }
    else if (pth->PBH_ADAF_delta == 0.5){
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
      epsilon = epsilon_0 * pow(M_b_dot / M_crit,a);
    }

    L_acc = epsilon*M_b_dot*_c_*_c_;

  }

  /** Spherical accretion from Ali-Haimoud et al. 1612.05644 */

  else if(pth->PBH_accretion_recipe == spherical_accretion){
    rho_cmb = pvecback[pba->index_bg_rho_g]/pow(_Mpc_over_m_,2)*3/8./_PI_/_G_*pow(_c_,4.)* 6.241509e12; // [MeV/m^3]
    // x_e_infinity = 1; // change to 1 for the strong-feedback case
    x_e_infinity = x_e; // change to x_e for the no-feedback case
    v_B = sqrt((1+x_e_infinity)*T_infinity/m_p)*_c_; //sound speed.
    if(pth->PBH_relative_velocities < 0.){
            v_l = 30*MIN(1,z/1000)*1e3; // [m/s]
            if(v_B < v_l) v_eff = sqrt(v_B*v_l);
            else v_eff = v_B;
    }
    else{
            v_l = pth->PBH_relative_velocities*1e3; // [m/s]
            v_eff = pow(v_l*v_l+v_B*v_B,0.5);
    }
    r_B = _G_*pth->PBH_accreting_mass*M_sun*pow(v_eff,-2); // [m]
    t_B = _G_*pth->PBH_accreting_mass*M_sun/pow(v_eff,3); // [s]
    beta_compton_drag = 4./3*x_e_infinity*_sigma_*rho_cmb*t_B/(m_p)*_c_;
    gamma_cooling = 2*m_p/(m_e*(1+x_e_infinity))*beta_compton_drag;
    lambda_iso = 0.25*exp(1.5);
    lambda_ad = 0.25*pow(3./5,1.5);
    lambda_1 = lambda_ad+(lambda_iso-lambda_ad)*pow(gamma_cooling*gamma_cooling/(88+gamma_cooling*gamma_cooling),0.22);
    lambda_2 = exp(4.5/(3+pow(beta_compton_drag,0.75)))*1/(pow(pow(1+beta_compton_drag,0.5)+1,2));
    lambda = lambda_1*lambda_2/lambda_iso;
    rho = pvecback[pba->index_bg_rho_b]/pow(_Mpc_over_m_,2)*3/8./_PI_/_G_*_c_*_c_; // [kg/m^3]
    M_b_dot = 4*_PI_*lambda*rho*r_B*r_B*v_eff; // [kg/s]
    T_ion = 1.5e4*_eV_over_Kelvin_;
    tau_cooling = 1.5/(5+pow(gamma_cooling,2./3));
    Y_s = pow((1+x_e_infinity)/2,2./3*13.6/T_ion)*tau_cooling/4*pow(1-5./2*tau_cooling,1./3)*m_p/m_e;
    T_s = m_e * Y_s*pow(1+Y_s/0.27,-1./3); // [MeV]
    if(T_s/m_e > 1){
      J = 27/(2*_PI_)*(log(2*T_s/(m_e)*exp(-0.577)+0.08)+4./3);
    }
    else{
      J = 4/_PI_*sqrt(2/_PI_)*pow(T_s/m_e,-0.5)*(1+5.5*pow(T_s/m_e,1.25));
    }
    L_ed = 4*_PI_*_G_*pth->PBH_accreting_mass*M_sun*m_p*1e6/_eV_over_joules_/(_sigma_*_c_);
    L_acc = 1./137*T_s/(m_p)*J*pow(M_b_dot*_c_*_c_,2)/L_ed;
  }

  *energy_rate =  (rho_cdm_today/(pth->PBH_accreting_mass*M_sun*_c_*_c_))*pow(1+z,3)*L_acc*pth->PBH_fraction;

  free(pvecback);

}




/**
 * In case of non-minimal cosmology, this function determines the energy rate injected in the IGM at a given redshift z
 * (= on-the-spot annihilation). This energy injection may come e.g. from dark matter annihilation or decay.
 *
 * @param ppr            Input: pointer to precision structure
 * @param pba            Input: pointer to background structure
 * @param ptw            Input: pointer to thermo_workspace structure
 * @param z              Input: redshift
 * @param energy_rate    Output: energy density injection rate
 * @param error_message  Output: error message
 * @return the error status
 */
int thermodynamics_solve_onthespot_energy_injection(struct precision * ppr,
                                                    struct background * pba,
                                                    struct thermo_workspace * ptw,
                                                    double z,
                                                    double * energy_rate,
                                                    ErrorMsg error_message){

  /** Summary: */

  /** Define local variables */
  struct thermo_heating_parameters* pthp = ptw->pthp;

  double annihilation_at_z;
  double rho_cdm_today;
  double u_min;
  double erfc;

  /* redshift-dependent annihilation parameter */

  if (z>pthp->annihilation_zmax) {

    annihilation_at_z = pthp->annihilation*
      exp(-pthp->annihilation_variation*pow(log((pthp->annihilation_z+1.)/(pthp->annihilation_zmax+1.)),2));
  }
  else if (z>pthp->annihilation_zmin) {

    annihilation_at_z = pthp->annihilation*
      exp(pthp->annihilation_variation*(-pow(log((pthp->annihilation_z+1.)/(pthp->annihilation_zmax+1.)),2)
                                         +pow(log((z+1.)/(pthp->annihilation_zmax+1.)),2)));
  }
  else {

    annihilation_at_z = pthp->annihilation*
      exp(pthp->annihilation_variation*(-pow(log((pthp->annihilation_z+1.)/(pthp->annihilation_zmax+1.)),2)
                                         +pow(log((pthp->annihilation_zmin+1.)/(pthp->annihilation_zmax+1.)),2)));
  }

  rho_cdm_today = pow(ptw->SIunit_H0,2)*3/8./_PI_/_G_*pba->Omega0_cdm*_c_*_c_; // in J/m^3

  u_min = (1+z)/(1+pthp->annihilation_z_halo);

  erfc = pow(1.+0.278393*u_min+0.230389*u_min*u_min+0.000972*u_min*u_min*u_min+0.078108*u_min*u_min*u_min*u_min,-4);

  *energy_rate = pow(rho_cdm_today,2)/_c_/_c_*pow((1+z),3)*
    (pow((1.+z),3)*annihilation_at_z+pthp->annihilation_f_halo*erfc)
    +rho_cdm_today*pow((1+z),3)*pthp->decay;  // in J/m^3/s

  return _SUCCESS_;

}


#endif
