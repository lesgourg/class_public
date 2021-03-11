/******************************************************************************************************/
/*                           HYREC-2: Hydrogen and Helium Recombination Code                          */
/*                      Written by Yacine Ali-Haimoud and Chris Hirata (2010-17)                      */
/*                            with contributions from Nanoom Lee (2020)                               */
/*                                                                                                    */
/*     energy_injection.c: functions for the energy injection rate by various physical processes      */
/*                                                                                                    */
/******************************************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "hyrectools.h"
#include "energy_injection.h"


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
  
  /* Dark matter annihilation in the smooth background */
  if (params->pann > 0.) {

    /* Parametrized variation of pann */
    if (zp1 > zp1_max) pann_tot = params->pann *exp(-var *square(log(zp1_ann/zp1_max)));
    else if (zp1 > zp1_min) {
      pann_tot = params->pann *exp(var*(-square(log(zp1_ann/zp1_max))
				+square(log(zp1/zp1_max))));
    }
    else {
      pann_tot = params->pann *exp(var*(-square(log(zp1_ann/zp1_max))
				+square(log(zp1_min/zp1_max))));
    }
	pann_tot = pann_tot*pow(zp1,3.);
  }
  
  /* Dark matter annihilation in haloes */
  if (params->pann_halo > 0.) {
    u_min = zp1/zp1_halo;
    erfc  = pow(1.+0.278393*u_min+0.230389*u_min*u_min+0.000972*u_min*u_min*u_min+0.078108*u_min*u_min*u_min*u_min,-4); 
    pann_tot += params->pann_halo *erfc;
  }
  
  return square(10537.4*params->odmh2) * zp1*zp1*zp1*1e-9* pann_tot
        +10537.4*params->odmh2*pow((1+z),3)*params->decay;
  /* the prefactor is 3 H100^2/(8 Pi G) c^2 in eV/cm^3, H100 = 100km/s/Mpc */
  /* pann is given in cm^3/s/GeV, multiply by 1e-9 to get cm^3/s/eV */
  
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
  vB    = 9.09e3 * sqrt((1.+xe)*Teff);    /* Bondi speed in cm/s */
  tB    = 1.33e26 *Mpbh/vB/vB/vB;         /* Bondi timescale in sec*/

  return 5.60e-24 *xe/a/a/a/a *tB;       
}

 /* Dimensionless Compton cooling rate */
double gamma_pbh(double Mpbh, double z, double xe, double Teff) {
  return  2.*1836./(1.+xe) *beta_pbh(Mpbh, z, xe, Teff);
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
  
  lam_nodrag = lam_ad + (lam_iso - lam_ad) * pow(gamma*gamma/(88. + gamma*gamma), 0.22); /* Fitting formula for the no-drag case */
   
  return lam_ricotti *lam_nodrag /lam_iso;
}

/* Accretion rate (in g/s), accounting for Compton drag and Compton cooling.*/
double Mdot_pbh(double Mpbh, double z, double xe, double Teff) {

  double vB  = 9.09e3 * sqrt((1.+xe)*Teff);    /* Bondi speed in cm/s */
   
  return 8.81e22 * Mpbh*Mpbh*cube((1.+z)/vB) *lambda_pbh(Mpbh, z, xe, Teff);   
}

/* Temperature of the flow near the Shchwarzschild radius divided by m_e c^2 */
double TS_over_me_pbh(double Mpbh, double z, double xe, double Teff) {
  double gamma, tau, omega, Y;
  
  gamma = gamma_pbh(Mpbh, z, xe, Teff);

  tau = 1.5/(5. + pow(gamma, 2./3.));     /* T/Teff -> tau *rB/r        for r << rB */
  omega = sqrt(2. - 5*tau);               /* v/vB -> -omega *sqrt(rB/r) for r << rB */

  Y = pow((1.+ xe)/2., 7.) * tau/2. *pow(omega/4.,2./3.)*1836.;

  return Y /pow(1.+Y/0.27, 1./3.); 
}

/* Radiative efficiency divided by \dot{m} */
double eps_over_mdot_pbh(double Mpbh, double z, double xe, double Teff) {
  double X, G;

  X = TS_over_me_pbh(Mpbh, z, xe, Teff);
  
  /* Fit to the (e-e + e-p) free-free Gaunt factor */
  if (X < 1) G = 4./M_PI * sqrt(2./M_PI/X) *(1.+ 5.5*pow(X, 1.25));
  else       G = 27./2./M_PI *(log(2.*X*0.56146 + 0.08) + 4./3.);

  return X/1836./137. * G;   /* alpha[fine-struct] * TS/mp * G */
}

/* Luminosity of a single PBH in erg/s*/
double L_pbh(double Mpbh, double z, double xe, double Teff) {
  double Mdot, mdot, eff;

  Mdot = Mdot_pbh(Mpbh, z, xe, Teff); 
  mdot = Mdot / (1.4e17 * Mpbh);        /* Mdot c^2 / L_Eddington */ 
  eff  = mdot *eps_over_mdot_pbh(Mpbh, z, xe, Teff); 
  
  return eff * Mdot * 9e20;  /* L = epsilon Mdot c^2 */  
}

/* Very approximate value of the rms relative velocity, in cm/s */
double vbc_rms_func(double z) {
  if (z < 1e3) return 30e5 * (1.+z)/1e3;
  else         return 30e5;
}

/* Average of the pbh luminosity over the distribution of relative velocities */
double L_pbh_av(double Mpbh, double z, double xe, double Tgas) {
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

    Teff = Tgas + 121e-10 *vbc*vbc/(1.+xe); 
    
    num += L_pbh(Mpbh, z, xe, Teff) * P_vbc;
  }

  return num/denom;
}

/* Rate of energy *injection* per unit volume (in erg/cm^3/s) due to PBHs*/
double dEdtdV_pbh(double fpbh, double Mpbh, double z, double xe, double Tgas) { 
  double xe_used = (xe < 1.? xe : 1.); /* Since we are not accounting for Helium */
  if (fpbh > 0.) {
    return 6.12e-52/Mpbh * cube(1.+z) * fpbh *L_pbh_av(Mpbh, z, xe_used, Tgas);
  }
  else return 0.;
}

/***********************************************************************************
Total energy *injection* rate per unit volume
Add in your favorite energy injection mechanism here
***********************************************************************************/

double dEdtdV_inj(double z, double xe, double Tgas, INJ_PARAMS *params){
  return   dEdtdV_DM_ann(z, params)
         + dEdtdV_pbh(params->fpbh, params->Mpbh, z, xe, Tgas);  
}


/**********************************************************************************
Energy *deposition* rate per unit volume
Essentially use a very simple ODE solution
**********************************************************************************/

void update_dEdtdV_dep(double z_out, double dlna, double xe, double Tgas,
		       double nH, double H, INJ_PARAMS *params, double *dEdtdV_dep) {

  double inj  = dEdtdV_inj(z_out, xe, Tgas, params);
  
  if (params->on_the_spot == 1) *dEdtdV_dep = inj;

  // Else put in your favorite recipe to translate from injection to deposition
  // Here I assume injected photon spectrum Compton cools at rate dE/dt = - 0.1 n_h c sigma_T E
  // This is valid for E_photon ~ MeV or so.
  
  else { // c sigma_T = 2e-14 (cgs)
    *dEdtdV_dep = (*dEdtdV_dep *exp(-7.*dlna) + 2e-15* dlna*nH/H *inj)
                 /(1.+ 2e-15 *dlna*nH/H);                              
  }       
}

/*******************************************************************************
Fraction of energy deposited in the form of heat, ionization and excitations
*******************************************************************************/

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


