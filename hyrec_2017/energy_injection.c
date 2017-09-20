/******************************************************************************************************/
/*                           HYREC: Hydrogen and Helium Recombination Code                            */
/*                      Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                      */
/*                                                                                                    */
/*     energy_injection.c: functions for the energy injection rate by various physical processes      */
/*                                                                                                    */
/*     Written October 2016                                                                          */
/*                                                                                                    */
/******************************************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "include/hyrectools.h"
#include "include/energy_injection.h"


/***************************************************************************************
Total volumic rate of energy *injection*, in eV/cm^3/s due to DM annihilation
in the smooth background and in haloes as in Giesen et al 1209.0247
****************************************************************************************/

double dEdtdV_DM_ann(double z, INJ_PARAMS *params){

  double pann_tot, u_min, erfc;
  double zp1, zp1_ann, zp1_max, zp1_min, var, zp1_halo;

  var       = params->ann_var;
  zp1       = z + 1.;
  zp1_ann   = params->ann_z + 1.;
  zp1_max   = params->ann_zmax + 1.;
  zp1_halo  = params->ann_z_halo + 1.;
  zp1_min   = params->ann_zmin + 1.;

  pann_tot = 0.;
  double Boost_factor = 0.;
  /* Dark matter annihilation in the smooth background */
  // if (params->pann > 0.) {
  //
  //   /* Parametrized variation of pann */
  //   if (zp1 > zp1_max) pann_tot = params->pann *exp(-var *square(log(zp1_ann/zp1_max)));
  //   else if (zp1 > zp1_min) {
  //     pann_tot = params->pann *exp(var*(-square(log(zp1_ann/zp1_max))
	// 			+square(log(zp1/zp1_max))));
  //   }
  //   else {
  //     pann_tot = params->pann *exp(var*(-square(log(zp1_ann/zp1_max))
	// 			+square(log(zp1_min/zp1_max))));
  //   }
  //   pann_tot*=zp1*zp1*zp1;
  // }
  // /* Dark matter annihilation in haloes */
  // if (params->pann_halo > 0.) {
  //   u_min = zp1/zp1_halo;
  //   erfc  = pow(1.+0.278393*u_min+0.230389*u_min*u_min+0.000972*u_min*u_min*u_min+0.078108*u_min*u_min*u_min*u_min,-4);
  //   pann_tot += params->pann_halo *erfc;
  // }



    if(params->ann_f_halo>0.){
      u_min = zp1/zp1_halo;
      erfc  = pow(1.+0.278393*u_min+0.230389*u_min*u_min+0.000972*u_min*u_min*u_min+0.078108*u_min*u_min*u_min*u_min,-4);
      Boost_factor = params->ann_f_halo*erfc/pow(zp1,3);
    }
    else Boost_factor = 0;

    return square(10537.4*params->odmh2)*1e-9*(pow((zp1),6)*params->pann)*(1+Boost_factor);

  // return square(10537.4*params->odmh2) * zp1*zp1*zp1 *1e-9* pann_tot;
  /* the prefactor is 3 H100^2/(8 Pi G) c^2 in eV/cm^3, H100 = 100km/s/Mpc */
  /* pann is given in cm^3/s/GeV, multiply by 1e-9 to get cm^3/s/eV */

}
/******************************Energy Injection DM decay**********************************/
double dEdtdV_DM_decay(double z, INJ_PARAMS *params){
  double result_integrale = 0, energy_rate = 0;
  double rho_cdm_today;
  double _Mpc_over_m_;
  _Mpc_over_m_ = 3.085677581282*pow(10,22);
  rho_cdm_today = params->Omega0_cdm*params->H0*params->H0/(10000)*1.44729366e-9; /* energy density in Kg/m^3 */

  if(params->decay_fraction>0){

     result_integrale = exp(-params->Gamma_dcdm*2*((params->Omega0_b+params->Omega0_cdm)*pow(params->Omega0_r+(params->Omega0_b+params->Omega0_cdm)/(1+z),0.5)
    +2*pow(params->Omega0_r,1.5)*(1+z)-2*params->Omega0_r*pow((1+z)*(params->Omega0_r*(1+z)+(params->Omega0_b+params->Omega0_cdm)),0.5))/(3*pow((params->Omega0_b+params->Omega0_cdm),2)*(1+z)*params->H0));
    energy_rate=rho_cdm_today*pow((1+z),3)*params->decay_fraction*result_integrale*(params->Gamma_dcdm*2.99792458e8/_Mpc_over_m_)/1.e6/1.60217653e-19;
    // if(rho_ini_dcdm!=0)energy_rate=rho_ini_dcdm*pow((1+z),3)*param->decay*result_integrale*(param->Gamma_dcdm*2.99792458e8/_Mpc_over_m_)/1.e6/1.60217653e-19;
    // // else energy_rate=rho_cdm_today*pow((1+z),3)*param->decay/1.e6/1.60217653e-19;
    // else energy_rate=rho_cdm_today*pow((1+z),3)*param->decay*result_integrale*(param->Gamma_dcdm*2.99792458e8/_Mpc_over_m_)/1.e6/1.60217653e-19;
    // // if(z>0)fprintf(stdout, "z = %e energy_rate = %e\n", z, energy_rate*1.e6*1.60217653e-19);
  }
  else energy_rate = 0.;
  // fprintf(stdout, "z %e result_integrale %e energy_rate %e\n",z,result_integrale,energy_rate );


  return energy_rate;


}
/***************************************************************************************
Effect of accreting primordial black holes
Since the accuracy is certainly not at the percent level,
we assume best-fit values for all cosmological parameters and neglect Helium
Throughout, Mpbh is in solar masses, Teff in Kelvins
***************************************************************************************/


/* Dimensionless Compton drag rate */
double beta_pbh(double Mpbh, double z, double xe, double Teff) {
  double a, vB, tB;

  a     = 1./(1.+z);
  // vB    = vbc_av(z,xe,Teff);    /* Bondi speed in cm/s */
  vB    = 9.09e3 * sqrt((1.+xe)*Teff);    /* Bondi speed in cm/s */
  tB    = 1.33e26 *Mpbh/vB/vB/vB;         /* Bondi timescale in sec*/
  // if(isnan(vB)==0)printf("%e %e %e %e %e\n",z,vB,xe,Teff,7.45e-24 *xe/a/a/a/a *tB);
  return 7.45e-24 *xe/a/a/a/a *tB;
}

 /* Dimensionless Compton cooling rate */
double gamma_pbh(double Mpbh, double z, double xe, double Teff) {
  return  3.67e3/(1.+xe) *beta_pbh(Mpbh, z, xe, Teff);
}

/* Dimensionless accretion rate */
double lambda_pbh(double Mpbh, double z, double xe, double Teff) {
  double beta, gamma, lam_ricotti, lam_ad, lam_iso, lam_nodrag;

  beta  = beta_pbh(Mpbh, z, xe, Teff);
  gamma = gamma_pbh(Mpbh, z, xe, Teff);

  lam_ricotti = exp(4.5/(3.+pow(beta, 0.75)))/square(sqrt(1.+beta)+1.);
  /* Fitting formula from Ricotti 2007 for the fully isothermal case */
  lam_ad      = pow(0.6, 1.5)/4.;
  lam_iso     = exp(1.5)/4.;

  lam_nodrag = lam_ad + (lam_iso - lam_ad) * pow(gamma*gamma/(88. + gamma*gamma), 0.22);
  /* Fitting formula for the no-drag case */

  return lam_ricotti *lam_nodrag /lam_iso;
}

/* Accretion rate (in g/s), accounting for Compton drag and Compton cooling.*/
/* This is assuming Omega_b h^2 = 0.022 (sufficient at the level of accuracy we need) */
double Mdot_pbh(double Mpbh, double z, double xe, double Teff) {

  // double vB  = vbc_av(z,xe,Teff);    /* Bondi speed in cm/s */
  double vB  = 9.09e3 * sqrt((1.+xe)*Teff);    /* Bondi speed in cm/s */
  // printf("%e %e \n",z,9.15e22 * Mpbh*Mpbh*cube((1.+z)/vB) *lambda_pbh(Mpbh, z, xe, Teff));
  return 9.15e22 * Mpbh*Mpbh*cube((1.+z)/vB) *lambda_pbh(Mpbh, z, xe, Teff);
}

/* Temperature of the flow near the Shchwarzschild radius divided by m_e c^2 */
/** Changed 01/26/2017: added a parameter "coll_ion":
    if set to 1, assume collisional ionizations
    if set to 0, assume photoionizations.
**/

double TS_over_me_pbh(double Mpbh, double z, double xe, double Teff, int coll_ion) {
  double gamma, tau, YS;

  gamma = gamma_pbh(Mpbh, z, xe, Teff);

  tau = 1.5/(5. + pow(gamma, 2./3.));     /* T/Teff -> tau *rB/r        for r << rB */

  YS = 2./(1.+xe) * tau/4. *pow(1.-2.5*tau,1./3.)*1836.;

  if (coll_ion == 1) YS *= pow((1.+ xe)/2., 8.);

  return YS /pow(1.+YS/0.27, 1./3.);
}

/* Radiative efficiency divided by \dot{m} */
double eps_over_mdot_pbh(double Mpbh, double z, double xe, double Teff, int coll_ion) {
  double X, G;

  X = TS_over_me_pbh(Mpbh, z, xe, Teff, coll_ion);

  /* Fit to the (e-e + e-p) free-free Gaunt factor */
  if (X < 1) G = 4./M_PI * sqrt(2./M_PI/X) *(1.+ 5.5*pow(X, 1.25));
  else       G = 13.5/M_PI *(log(2.*X*0.56146 + 0.08) + 4./3.);

  return X/1836./137. * G;   /* alpha[fine-struct] * TS/me *mp/me * G */
}

/* Luminosity of a single PBH in erg/s*/
double L_pbh(double Mpbh, double z, double xe, double Teff, int coll_ion) {
  double Mdot, mdot, eff;

  Mdot = Mdot_pbh(Mpbh, z, xe, Teff);
  mdot = Mdot / (1.4e17 * Mpbh);        /* Mdot c^2 / L_Eddington */
  eff  = mdot *eps_over_mdot_pbh(Mpbh, z, xe, Teff, coll_ion);
  // printf("%e %e %e\n",z,mdot,eff*Mdot/(1.4e17 * Mpbh));
  return eff * Mdot * 9e20;  /* L = epsilon Mdot c^2 */
}

/* Very approximate value of the rms relative velocity, in cm/s */
double vbc_rms_func(double z) {
  if (z < 1e3) return 3e6 * (1.+z)/1e3;
  else         return 3e6;
}

/* Average of the relative velocity over the distribution of relative velocities */

double vbc_av(double z, double xe, double T) {
  double vB  = 9.09e3 * sqrt((1.+xe)*T);    /* Bondi speed in cm/s */
  if(vB < vbc_rms_func(z)) vB = sqrt(vB*vbc_rms_func(z));
  // printf("%e %e %e %e %e\n",z,vB,9.09e3 * sqrt((1.+xe)*T),xe,T);
  return vB;
}

/* Average of the pbh luminosity (erg/s) over the distribution of relative velocities */

double L_pbh_av(double Mpbh, double z, double xe, double Tgas, int coll_ion) {
  double vbc_rms, vbc_max, vbc, P_vbc, num, denom, x, Teff;
  int i, Nvbc;

  Nvbc = 50; // More than enough at the level of precision we use

  vbc_rms = vbc_rms_func(z);

  vbc_max = 5.*vbc_rms;

  num = denom = 0.;
  for (i = 0; i < Nvbc; i++) {
    vbc    = i*vbc_max/(Nvbc-1.);
    x      = vbc/vbc_rms;
    denom += P_vbc = x*x* exp(-1.5*x*x);

    Teff = Tgas + 1.21e-8 *vbc*vbc/(1.+xe);

    num += L_pbh(Mpbh, z, xe, Teff, coll_ion) * P_vbc;
  }
  // printf(" %e %e %e \n", z,Mdot_pbh(Mpbh, z, xe, Tgas),num/denom/(1.4e17 * Mpbh)/9e20);
  return num/denom;
}

/* Rate of energy *injection* per unit volume (in eV/cm^3/s) due to PBHs */
/* Taking Omega_c h^2 = 0.12                                             */

double dEdtdV_accreting_PBH(double z, double xe, double Tgas, INJ_PARAMS *params) {
  double xe_used = (xe < 1.? xe : 1.); /* Since we are not accounting for Helium */
  double energy_rate = 0;
  // xe_used = feedback + (1-feedback)*xe_used;
  // Old feedback - no feedback model, where xe = 1 throughout was assumed in the strong feedback case
  if (params->fpbh > 0. && params->Mpbh > 0.) {

    if(params->PBH_accretion_recipe == 0){
      // energy_rate = 7.07e-52/Mpbh * cube(1.+z) * fpbh *L_pbh(Mpbh, z, xe_used, Tgas, coll_ion); /** alternative computation using veff instead of the numerical average */
      energy_rate = 7.07e-52/params->Mpbh * cube(1.+z) * params->fpbh *L_pbh_av(params->Mpbh, z, xe_used, Tgas, params->coll_ion);
    }
    else if(params->PBH_accretion_recipe == 1){
      energy_rate = 7.07e-52/params->Mpbh * cube(1.+z) * params->fpbh *L_pbh_ADAF(z, xe_used, Tgas, params);
    }

   }
  return energy_rate;
}

double L_pbh_ADAF(double z, double xe, double Tgas, INJ_PARAMS *params) {
  double L_ed = (1.4e17 * params->Mpbh) * 9e20, M_ed_dot = 10 * L_ed/(9e20), M_crit = 0.01*M_ed_dot, lambda = params->PBH_accretion_eigenvalue; // slightly different definition of M_ed_dot in the ADAF formalism. See e.g. Yuan&Narayan 1406.0586
  double Value_min, Value_med, Value_max, a=0, epsilon_0=0.1,epsilon=0.,Lacc=0.;
  double veff = vbc_av(z,xe,Tgas);
  double Mdot = 9.15e22 * params->Mpbh*params->Mpbh*cube((1.+z)/veff) * lambda;

  if(params->PBH_ADAF_delta == 1e-3){
    Value_min = 7.6e-5;
    Value_med = 4.5e-3;
    Value_max = 7.1e-3;
    if(Mdot/M_ed_dot <= Value_min){
      epsilon_0 = 0.065;
      a = 0.71;
    }
    else if(Value_min < Mdot/M_ed_dot && Mdot/M_ed_dot  <= Value_med){
      epsilon_0 = 0.020;
      a = 0.47;
    }
    else if(Value_med < Mdot/M_ed_dot && Mdot/M_ed_dot <= Value_max){
      epsilon_0 = 0.26;
      a = 3.67;
    }
    else{
      epsilon_0 = 0.1;
      a = 0.;
    }
    epsilon = epsilon_0 * pow(Mdot / M_crit,a);
  }
  else if (params->PBH_ADAF_delta == 0.1){
    Value_min = 9.4e-5;
    Value_med = 5e-3;
    Value_max = 6.6e-3;
    if(Mdot/M_ed_dot <= Value_min){
      epsilon_0 = 0.12;
      a = 0.59;
    }
    else if(Value_min < Mdot/M_ed_dot && Mdot/M_ed_dot  <= Value_med){
      epsilon_0 = 0.026;
      a = 0.27;
    }
    else if(Value_med < Mdot/M_ed_dot && Mdot/M_ed_dot <= Value_max){
      epsilon_0 = 0.50;
      a = 4.53;
    }
    else{
      epsilon_0 = 0.1;
      a = 0.;
    }
    epsilon = epsilon_0 * pow(Mdot / M_crit,a);
  }
  else if (params->PBH_ADAF_delta == 0.5){

    Value_min = 2.9e-5;
    Value_med = 3.3e-3;
    Value_max = 5.3e-3;
    if(Mdot/M_ed_dot <= Value_min){
      epsilon_0 = 1.58;
      a = 0.65;
    }
    else if(Value_min < Mdot/M_ed_dot && Mdot/M_ed_dot  <= Value_med){
      epsilon_0 = 0.055;
      a = 0.076;
    }
    else if(Value_med < Mdot/M_ed_dot && Mdot/M_ed_dot <= Value_max){
      epsilon_0 = 0.17;
      a = 1.12;
    }
    else{
      epsilon_0 = 0.1;
      a = 0.;
    }
    epsilon = epsilon_0 * pow(Mdot / M_crit,a);
  }

  Lacc = epsilon*Mdot*9e20;
  // printf("%e %e %e \n",z,Mdot/(1.4e17 * params->Mpbh),Lacc/L_ed);
  return Lacc;


}
/******************************Energy Injection low mass PBH (evaporation)**********************************/
double dEdVdt_evaporating_PBH(double z, INJ_PARAMS *params){

  double rho_cdm_today;
  //double tau;
  int last_index_back;
  //Parameters related to PBH
  ErrorMsg error_message;
  double f, f_neutrinos, em_branching, pbh_mass;
  double dMdt,energy_rate=0.;
  double zp1 = 1 + z;

  /* Calculate the PBH-mass evolution at first call of the function */
  // if ((params->PBH_table_is_initialized) == _FALSE_) {
  //   params->PBH_table_is_initialized = _TRUE_;
  //   pbh_low_mass_time_evolution(ppr,pba,params,error_message);
  // }
  // class_test(params->PBH_table_is_initialized == _FALSE_, error_message, "The PBH table is not initialized");
  // Need to add a error message in case the PBH table isn't correctly initialised
  /* End of PBH-mass loop */
 if(params->PBH_low_mass > 0){
   class_call(array_interpolate_spline(params->PBH_table_z,
 				      params->PBH_table_size,
 				      params->PBH_table_mass,
 				      params->PBH_table_mass_dd,
 				      1,
 				      z,
 				      &last_index_back,
 				      &(pbh_mass),
 				      1,
 				      error_message),
 	     error_message,
 	     error_message);
   class_call(array_interpolate_spline(params->PBH_table_z,
 				      params->PBH_table_size,
 				      params->PBH_table_F,
 				      params->PBH_table_F_dd,
 				      1,
 				      z,
 				      &last_index_back,
 				      &f,
 				      1,
 				      error_message),
 	     error_message,
 	     error_message);
   // printf("z %e pbhmass %e f %e\n",z,pbh_mass,f);
   // f_neutrinos = 6*0.147;
   // em_branching = (f-f_neutrinos)/f;
   em_branching = 1.; // Currently incoporated in the computation of the f(z) functions.
   // printf("params->PBH_z_evaporation %e\n", params->PBH_z_evaporation);
   if(pbh_mass <= 0.0001*params->PBH_low_mass || f <= 0 || isnan(pbh_mass)==1 || isnan(f)==1 || z < params->PBH_z_evaporation){
     pbh_mass = 0;
     dMdt = 0;
     f = 0.;
   }
   else {
     dMdt=5.34e-5*f*pow(pbh_mass/1e10,-2)*1e10;
   }

    energy_rate = 10537.4*params->odmh2*pow((zp1),3)*params->fpbh/params->PBH_low_mass*em_branching*(dMdt);

   // *energy_rate = rho_cdm_today*pow((1+z),3)*params->PBH_fraction/pbh_mass*em_branching*(dMdt);
   if(isnan(energy_rate)==1 || energy_rate < 0){
     energy_rate=0.;
   }
   // if(pbh_mass>0)fprintf(stdout,"z = %lg | f = %lg | mass = %lg | energy_rate = %lg | params->PBH_low_mass = %lg \n",z,f,pbh_mass,energy_rate,params->PBH_low_mass);

 }
 else energy_rate = 0.;

  return energy_rate;
  // if(pbh_mass>0)fprintf(stdout,"%e %e %e %e \n",z,f,pbh_mass,*energy_rate);
}
// double dEdtdV_stars(double z, INJ_PARAMS *params){
//   double ainv = 1+z;
//   double rho_sfr =params->ap*pow(ainv,params->bp)/(1+pow((ainv)/params->cp,params->dp))*ainv*ainv*ainv*(1+tanh((params->z_start_reio_stars-z1)))/2;//Comoving to physical
//   double dNion_over_dt;
//   double joules_to_ev = 6.24150647996e+18;
//   double MPCcube_to_mcube=pow(3.085677581282e22,3);
//   double stars_xe;
//   dNion_over_dt = params->f_esc*params->Zeta_ion*rho_sfr;
//   double L_x = param->Ex  * param->fx* 2*rho_sfr/(3*MPCcube_to_mcube*kBoltz*nH*H*(1.+xe_in+param->fHe));
//   return dNion_over_dt/(MPCcube_to_mcube);
//
// }
/***********************************************************************************
Total energy *injection* rate per unit volume
Add in your favorite energy injection mechanism here
***********************************************************************************/

double dEdtdV_inj(double z, double xe, double Tgas, INJ_PARAMS *params){

  return   dEdtdV_DM_ann(z, params)
    + dEdtdV_accreting_PBH(z, xe, Tgas, params)
    +dEdVdt_evaporating_PBH(z,params)
    +dEdtdV_DM_decay(z, params);
}


/**********************************************************************************
Energy *deposition* rate per unit volume
Essentially use a very simple ODE solution
**********************************************************************************/

void update_dEdtdV_dep(double z_out, double dlna, double xe, double Tgas,
		       double nH, double H, INJ_PARAMS *params, double *dEdtdV_dep) {

  double inj  = dEdtdV_inj(z_out, xe, Tgas, params);

  if (params->on_the_spot == 1){
    *dEdtdV_dep = inj;
  }

  // Else put in your favorite recipe to translate from injection to deposition
  // Here I assume injected photon spectrum Compton cools at rate dE/dt = - 0.1 n_h c sigma_T E
  // This is valid for E_photon ~ MeV or so.

  else { // 0.1 c sigma_T = 2e-15 (cgs)
    if(params->energy_deposition_treatment == 0)
    *dEdtdV_dep = (*dEdtdV_dep *exp(-7.*dlna) + 2e-15* dlna*nH/H *inj)
                 /(1.+ 2e-15 *dlna*nH/H);
    else if(params->energy_deposition_treatment == 1){
      if(params->energy_repart_functions > 0) hyrec_annihilation_coefficients_interpolate(params,z_out);
      else params->f_eff = 1.;
    *dEdtdV_dep = inj*params->f_eff;
    }

  }

  //  printf("on_the_spot %d *dEdtdV_dep %e inj %e \n", params->on_the_spot,*dEdtdV_dep,inj);


}

/*******************************************************************************
Fraction of energy deposited in the form of heat, ionization and excitations
*******************************************************************************/

int evaluate_chi_heat(INJ_PARAMS *param,double z, double xe){
   if(param->energy_repart_functions==0){
     hyrec_annihilation_coefficients_interpolate(param,z);
   }
   if(param->energy_repart_functions==1){
   hyrec_annihilation_coefficients_interpolate(param,xe);
   }


   /* old approximation from Chen and Kamionkowski */
   if(xe<1){
   if(param->energy_repart_functions==2){
    param->chi_heat = (1.+2.*xe)/3.; // old approximation from Chen and Kamionkowski
   }
   /* coefficient as revised by Slatyer et al. 2013 (in fact it is a fit by Vivian Poulin of columns 1 and 2 in Table V of Slatyer et al. 2013): */
   if(param->energy_repart_functions==3){
    param->chi_heat = 0.996857*(1.-pow(1.-pow(xe,0.300134),1.51035));
   }
  }
   else {
     param->chi_heat = 1.;
   }

   if(param->chi_heat < 0)  param->chi_heat = 0.;
   if(param->chi_heat > 1)  param->chi_heat = 1.;

   return _SUCCESS_;
}

int evaluate_chi_ionisation(INJ_PARAMS *param,double z, double xe){
  if(param->energy_repart_functions==1){
  hyrec_annihilation_coefficients_interpolate(param,xe);
 }
 if(param->energy_repart_functions==0){
   hyrec_annihilation_coefficients_interpolate(param,z);
 }
 /* old approximation from Chen and Kamionkowski */
 if(xe<1){
      if(param->energy_repart_functions==2){
     param->chi_ionH = (1.-xe)/3.;
     param->chi_lya = param->chi_ionH;
     param->chi_ionHe=0;
   }
   /* coefficient as revised by Slatyer et al. 2013 (in fact it is a fit by Vivian Poulin of columns 1 and 2 in Table V of Slatyer et al. 2013): */
   if(param->energy_repart_functions==3){
     param->chi_ionH = 0.369202*pow(1.-pow(xe,0.463929),1.70237);
     param->chi_ionHe =0.0312604*pow(1.-pow(xe,0.200634),0.82247);
     param->chi_lya = 0.335597*pow(1.-pow(xe,0.375314),1.80722);
   }
 }
  else {
    param->chi_ionH = 0.;
    param->chi_ionHe = 0.;
    param->chi_lya = 0.;
  }
   param->chi_ionH=MAX(param->chi_ionH,0.);
   param->chi_ionHe=MAX(param->chi_ionHe,0.);
   param->chi_lya=MAX(param->chi_lya,0.);
   param->chi_ionH=MIN(param->chi_ionH,1.);
   param->chi_ionHe=MIN(param->chi_ionHe,1.);
   param->chi_lya=MIN(param->chi_lya,1.);


   return _SUCCESS_;
}


double chi_heat(double xe) {
  return (1.+2.*xe)/3.; // Approximate value of Chen & Kamionkowski 2004

  // fit by Vivian Poulin of columns 1 and 2 in Table V of Galli et al. 2013
  // overall coefficient slightly changed by YAH so it reaches exactly 1 at xe = 1.
  // return (xe < 1.? 1.-pow(1.-pow(xe,0.300134),1.51035) : 1.);
}

double chi_ion(double xe) {
  return (1.-xe)/3.; // Approximate value of Chen & Kamionkowski 2004

  // fit by Vivian Poulin of columns 1 and 4 in Table V of Galli et al. 2013
  // return 0.369202*pow(1.-pow(xe,0.463929),1.70237);
}

double chi_exc(double xe) {
  return 1. - chi_ion(xe) - chi_heat(xe);
}




/***********************************************************************************************************/
/***********************Return values for the annihilation coefficients after interpolation*****************/
//
//
int hyrec_annihilation_coefficients_interpolate(INJ_PARAMS *inj_params,
                                                          double xe_or_z
                                                          ) {
      int last_index;
      ErrorMsg error_message;
      if(inj_params->energy_repart_functions < 2){
        array_interpolate_spline(inj_params->annihil_coef_xe,
                                  inj_params->annihil_coef_num_lines,
                                  inj_params->annihil_coef_heat,
                                  inj_params->annihil_coef_dd_heat,
                                  1,
                                  xe_or_z,
                                  &last_index,
                                  &(inj_params->chi_heat),
                                  1,
                                  error_message);

        array_interpolate_spline(inj_params->annihil_coef_xe,
                                  inj_params->annihil_coef_num_lines,
                                  inj_params->annihil_coef_lya,
                                  inj_params->annihil_coef_dd_lya,
                                  1,
                                  xe_or_z,
                                  &last_index,
                                  &(inj_params->chi_lya),
                                  1,
                                  error_message);

        array_interpolate_spline(inj_params->annihil_coef_xe,
                                  inj_params->annihil_coef_num_lines,
                                  inj_params->annihil_coef_ionH,
                                  inj_params->annihil_coef_dd_ionH,
                                  1,
                                  xe_or_z,
                                  &last_index,
                                  &(inj_params->chi_ionH),
                                  1,
                                  error_message);

        array_interpolate_spline(inj_params->annihil_coef_xe,
                                  inj_params->annihil_coef_num_lines,
                                  inj_params->annihil_coef_ionHe,
                                  inj_params->annihil_coef_dd_ionHe,
                                  1,
                                  xe_or_z,
                                  &last_index,
                                  &(inj_params->chi_ionHe),
                                  1,
                                  error_message);

        array_interpolate_spline(inj_params->annihil_coef_xe,
                                  inj_params->annihil_coef_num_lines,
                                  inj_params->annihil_coef_lowE,
                                  inj_params->annihil_coef_dd_lowE,
                                  1,
                                  xe_or_z,
                                  &last_index,
                                  &(inj_params->chi_lowE),
                                  1,
                                  error_message);
      }
      if(inj_params->energy_repart_functions > 0){
        array_interpolate_spline(inj_params->annihil_z,
                                inj_params->annihil_f_eff_num_lines,
                                inj_params->annihil_f_eff,
                                inj_params->annihil_dd_f_eff,
                                1,
                                xe_or_z,
                                &last_index,
                                &(inj_params->f_eff),
                                1,
                                error_message);
      }
    // fprintf(stdout,"%e      %e     %e      %e      %e    \n", xe_or_z,inj_params->chi_heat,inj_params->chi_lya, inj_params->chi_ionH,inj_params->chi_ionHe,inj_params->chi_lowE);

        return _SUCCESS_;


}
