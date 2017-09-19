/******************************************************************************************************/
/*                           HYREC: Hydrogen and Helium Recombination Code                            */
/*                      Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                      */
/*                                                                                                    */
/*         history.c: functions for numerical integration of the recombination history                */
/*                                                                                                    */
/*         Version: January 2015                                                                      */
/*                                                                                                    */
/*         Revision history:                                                                          */
/*            - written November 2010                                                                 */
/*            - January 2011: changed various switches (notably for post-Saha expansions)             */
/*                             so that they remain valid for arbitrary cosmologies                    */
/*            - November 2011: extended integration down to z = 0 with Peeble's model for z < 20      */
/*                             changed dTm/dlna so it can be called at all times                      */
/*            - May 2012:   - added explicit dependence on fine structure constant and electron mass  */
/*                          - modified call of rec_build_history                                      */
/*                             and improved numerical radiative transfer equations                    */
/*                             so the Lyman-lines spectrum can be extracted                           */
/*                           - split some functions for more clarity                                  */
/*             - October 2012: added some wrapper functions for running CAMB with HyRec               */
/*                             (courtesy of Antony Lewis)                                             */
/*                            - possibility to change fine structure constant/ electron mass          */
/*             -         2015: - added DM annihilation and 21 cm routines.                            */
/*                             - changed cosmological parameters input form                           */
/*			       - possibility to change fine structure constant/ electron mass         */
/*                             - nH0 now in cm^-3 (instead of m^-3 which was only used in helium.c)   */
/******************************************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#include "include/history.h"


#define DXHEII_MAX    1e-5       /* If xHeII - xHeII(Saha) < DXEHII_MAX, use post-Saha expansion for Helium. Lower value = higher accuracy. */
#define DXHII_MAX     3e-4       /* If xHII - xHII(Saha) < DXHII_MAX, use post-Saha expansion for Hydrogen. Switch to ODE integration after that.
                                        IMPORTANT: do not set to a lower value unless using a smaller time-step */
#define XHEII_MIN     1e-6       /* Stop considering Helium recombination once xHeII < XHEII_MIN */
#define DLNT_MAX      5e-4       /* Use the steady-state approximation for Tm as long as 1-Tm/Tr < DLNT_MAX, then switch to ODE integration */

/*************************************************************************************
Hubble expansion rate in sec^-1.
*************************************************************************************/

#ifdef CAMB

/* Use the Hubble rate from CAMB */

extern double dtauda_(double *);

double rec_HubbleRate(REC_COSMOPARAMS *cosmo, double z) {
  double a;

  a = 1./(1.+z);
  /* conversion from d tau/ da in Mpc to H(z) in 1/s */
  return 1./(a*a)/dtauda_(&a) /3.085678e22 * 2.99792458e8;
}

#else

/* Hyrec Hubble expansion rate.
   Note: neutrinos are assumed to be massless. */

double rec_HubbleRate(REC_COSMOPARAMS *cosmo, double z) {
   double a = 1./(1.+z);

   /* Total density parameter, including curvature */
   double rho = cosmo->omh2 /a/a/a    /* Matter */
              + cosmo->okh2 /a/a      /* Curvature */
              + cosmo->odeh2          /* Dark energy */
              + cosmo->orh2 /a/a/a/a; /* Radiation (photons + massless neutrinos) */

   /* Conversion to Hubble rate in sec-1 */
   return( 3.2407792896393e-18 * sqrt(rho) );
}

#endif


/*****************************************************************************************
Matter temperature -- 1st order steady state, from Hirata 2008.
The input and output temperatures are in KELVIN.
Added December 2014: possibility for additional energy deposition dEdtdV in eV/s/cm^3.
******************************************************************************************/

double rec_Tmss(double z, double xe, REC_COSMOPARAMS *cosmo, double dEdtdV) {

  double fsR = cosmo->fsR;
  double meR = cosmo->meR;
  double Tr  = cosmo->T0 *(1.+z);
  double nH  = cosmo->nH0 *cube(1.+z);
  double H   = rec_HubbleRate(cosmo, z);
  if(dEdtdV > 0) evaluate_chi_heat(cosmo->inj_params,z,xe);
  /* Coefficient = 8 sigma_T a_r / (3 m_e c) */
  /* Here Tr, Tm are the actual (not rescaled) temperatures */
  double coeff  = fsR*fsR/meR/meR/meR*4.91466895548409e-22*Tr*Tr*Tr*Tr*xe/(1.+xe+cosmo->fHe)/H;
  // double Tm = Tr/(1.+1./coeff)
  //           + (1.+2.*xe)/3.*dEdtdV/kBoltz /(1.5 *nH*(1.+xe+cosmo->fHe))/H /(1.+coeff);
  double Tm = Tr/(1.+1./coeff)
            + cosmo->inj_params->chi_heat*dEdtdV/kBoltz /(1.5 *nH*(1.+xe+cosmo->fHe))/H /(1.+coeff);

  return Tm;
}

/******************************************************************************************
Matter temperature evolution derivative. Input and output temperatures are in KELVIN.
Added May 2012: when Tm = Tr, return -Tr (needed by CLASS)
Corrected June 2016: in the presence of heating, Tm can cross Tr, so added a criterion
for returning -Tr when quasi-steady state. Note: this is not very "clean", there should
be some flag for quasi-steady-state, will eventually fix.
Added December 2014: possibility of additional energy deposition dEdtdV in eV/s/cm^3.
******************************************************************************************/

double rec_dTmdlna(double z, double xe, double Tm, REC_COSMOPARAMS *cosmo, double dEdtdV) {
  double fsR = cosmo->fsR;
  double meR = cosmo->meR;
  double Tr  = cosmo->T0 *(1.+z);
  double nH  = cosmo->nH0 *cube(1.+z);
  double H   = rec_HubbleRate(cosmo, z);
  if(dEdtdV > 0) evaluate_chi_heat(cosmo->inj_params,z,xe);

  return ( (Tr/Tm-1.<1e-10  && Tr > 3000.)  ? -Tr :
          -2.*Tm + fsR*fsR/meR/meR/meR*4.91466895548409e-22*Tr*Tr*Tr*Tr*xe/(1.+xe+cosmo->fHe)*(Tr-Tm)/H
	   + cosmo->inj_params->chi_heat*dEdtdV /kBoltz /(1.5 *nH*(1.+xe+cosmo->fHe))/H);
   /* Coefficient = 8 sigma_T a_r / (3 m_e c) */
   /* Here Tr, Tm are the actual (not rescaled) temperatures */
}


/************************************************************************
Explicit integrator.
Inputs: deriv: derivative at current time step
        deriv_prev: derivatives at two previous timesteps
        (deriv_prev[0] = previous time, deriv_prev[1] = 2 timesteps ago)
Output: change in x per dt
*************************************************************************/
double hyrec_integrator(double deriv, double deriv_prev[2]) {
  double result;


  // Second-order
  result = 1.25 * deriv - 0.25 *deriv_prev[1];

  // Third-order
  //return (5.*deriv + 8.*deriv_prev[1] - deriv_prev[2])/12.;

  // update derivatives
  deriv_prev[1] = deriv_prev[0];
  deriv_prev[0] = deriv;

  return result;
}

/**********************************************************************************************
Second-order explicit integrator for HeII abundance during HeII->HeI recombination.
Assumes hydrogen ionization is in Saha equilibrium (free electrons being provided by both H and He).
If xHeII is close enough to Saha equilibrium do a post-Saha expansion.
dxHeIIdlna_prev[0] = derivative at previous timestep
dxHeIIdlna_prev[1] = derivative 2 timesteps prior
z = input z
H = Hubble parameter(z)
***********************************************************************************************/

void rec_get_xe_next1_He(REC_COSMOPARAMS *cosmo, double z_in, double *xHeII,
			 double dxHeIIdlna_prev[2], int *post_saha) {

  double H, xH1, xH1_p, xH1_m, xHeIISaha, dxHeIISaha_dlna, DdxHeIIdlna_Dxe, dxHeIIdlna, z_out, Dxe;

  H = rec_HubbleRate(cosmo, z_in);

  xH1        = rec_saha_xH1(*xHeII, cosmo->nH0, cosmo->T0, z_in, cosmo->fsR, cosmo->meR);
  dxHeIIdlna = rec_helium_dxHeIIdlna(xH1, *xHeII, cosmo->nH0, cosmo->T0, cosmo->fHe, H, z_in, cosmo->fsR, cosmo->meR);

    /* Post-Saha approximation during the early phase of HeII->HeI recombination */
    if (*post_saha == 1) {
        z_out     = (1.+z_in)*exp(-DLNA)-1.;
        xHeIISaha = rec_saha_xHeII(cosmo->nH0, cosmo->T0, cosmo->fHe, z_out, cosmo->fsR, cosmo->meR);

        dxHeIISaha_dlna  = (1.+z_out)*(rec_saha_xHeII(cosmo->nH0, cosmo->T0, cosmo->fHe, z_out-0.5, cosmo->fsR, cosmo->meR)
                                      -rec_saha_xHeII(cosmo->nH0, cosmo->T0, cosmo->fHe, z_out+0.5, cosmo->fsR, cosmo->meR));

        Dxe    = 0.01*(cosmo->fHe - xHeIISaha);
        xH1_p = rec_saha_xH1(xHeIISaha+Dxe, cosmo->nH0, cosmo->T0, z_out, cosmo->fsR, cosmo->meR);
        xH1_m = rec_saha_xH1(xHeIISaha-Dxe, cosmo->nH0, cosmo->T0, z_out, cosmo->fsR, cosmo->meR);

	H = rec_HubbleRate(cosmo, z_out);
        DdxHeIIdlna_Dxe  = (rec_helium_dxHeIIdlna(xH1_p, xHeIISaha+Dxe, cosmo->nH0, cosmo->T0,
						  cosmo->fHe, H, z_out, cosmo->fsR, cosmo->meR)
                           -rec_helium_dxHeIIdlna(xH1_m, xHeIISaha-Dxe, cosmo->nH0, cosmo->T0,
						  cosmo->fHe, H, z_out, cosmo->fsR, cosmo->meR))/2./Dxe;

        *xHeII = xHeIISaha + dxHeIISaha_dlna/DdxHeIIdlna_Dxe;

        /* Check that the post-Saha approximation is still valid. If not, switch it off for future iterations */
        if (fabs(*xHeII - xHeIISaha) > DXHEII_MAX) *post_saha = 0;
    }

    /* Otherwise integrate ODE */
    else *xHeII += DLNA * hyrec_integrator(dxHeIIdlna, dxHeIIdlna_prev);

    //printf("%f %E %E %E\n", z_in, dxHeIIdlna, dxHeIIdlna_prev[0], dxHeIIdlna_prev[1]);
}

/***************************************************************************************************
Quasi-equilibrium solution for the hydrogen neutral density at z
Used at early stages when the system is stiff. After that switch to explicit second-order solver.
input:  xH1 at previous time step
output: xH1 at z.
iz_rad is the index for the radiation field at z
***************************************************************************************************/

void rec_xH1_stiff(int model, REC_COSMOPARAMS *cosmo, double z, double xHeII, double *xH1,
		   HYREC_ATOMIC *atomic, RADIATION *rad, unsigned iz_rad,
		   double dEdtdV, int *stiff, int *error){

  double Dxe, dx1s_dlna, dx1s_dlna2;
  double eps, Gamma, xeq, dxeq_dlna, nH, H, T;
  int i;

  nH = cosmo->nH0 *cube(1.+z);
  H  = rec_HubbleRate(cosmo, z);
  T  = cosmo->T0*(1.+z);


  // Stiff approximation: if dx1s_dlna = Gamma (x1s - xeq) with Gamma >> dxeq_dlna/xeq,
  //                      then x_1s \approx xeq + dxeq_dlna/Gamma

  // *****TO DO*****
  // right now I find xeq iteratively (it seems to work very well: very close to Saha equilibrium value)
  // and then approximate dxeq/dlna = (xeq - *xH1)/DLNA, since *xH1 is (nearly) the equilibrium value at the previous time
  // We may be able to improve on this last step and actually compute dxeq/dlna
  // This is what was done in the previous version (with post-Saha, while here we attempt to be more general)
  // and seems to work better.

  eps = 1e-4; // Step in numerical derivatives

  xeq = *xH1;
  for (i = 0; i < 2; i++) { // After 2 iterations it is well converged

    dx1s_dlna2 = -rec_dxHIIdlna(model, xHeII + 1.-xeq*(1.+eps), 1.-xeq*(1.+eps), nH, H, kBoltz*T, kBoltz*T,
			       atomic, rad, iz_rad, z, cosmo->fsR, cosmo->meR, dEdtdV, error, cosmo);

    dx1s_dlna  = -rec_dxHIIdlna(model, xHeII + 1.-xeq*(1.-eps), 1.-xeq*(1.-eps), nH, H, kBoltz*T, kBoltz*T,
			       atomic, rad, iz_rad, z, cosmo->fsR, cosmo->meR, dEdtdV, error, cosmo);

    /* Partial derivative of the time derivative wrt x1s, assuming dEdtdV depends weakly on x1s */
    Gamma = (dx1s_dlna2 - dx1s_dlna)/(2.*eps* xeq);

    dx1s_dlna = -rec_dxHIIdlna(model, xHeII + 1.-xeq, 1.-xeq, nH, H, kBoltz*T, kBoltz*T,
			      atomic, rad, iz_rad, z, cosmo->fsR, cosmo->meR, dEdtdV, error, cosmo);
    xeq -= dx1s_dlna/Gamma;
  }

  dxeq_dlna = (xeq - *xH1)/DLNA;
  *xH1 = xeq + dxeq_dlna/Gamma;

  if (fabs(*xH1/xeq -1) > 1e-2) *stiff = 0;

  if (*xH1 < 0. || *xH1 > 1. || *xH1 != *xH1) {
    fprintf(stderr, "\033[1m\033[31m error\033[22;30m in rec_xH1_stiff: at z = %f, xH1 = %E\n", z, xH1);
    *error = 1;
    return;
  }
  if (*error == 1) {
    fprintf(stderr, "  called from rec_xH1_stiff at z = %f\n", z);
    return;
  }

  /* printf("%f %E\n", z, xeq/rec_saha_xH1(xHeII, cosmo->nH0, cosmo->T0, z, cosmo->fsR, cosmo->meR)-1.); */

  /* if (*stiff == 0) { */
  /*  printf("switched off stiff solver at z = %f\n", z); */
  /*  exit(1); */
  /* } */

}


/*****************************************************************************************************
Second-order integrator (unless stiff system) used to evolve simultaneously HI and HeII.
When hydrogen is close to Saha equilibrium but there is still a significant amount of HeII,
use a post-Saha expansion for hydrogen. The other way around is not treated (simply integrate H and He until
there is almost no HeII left, then integrate H only)
******************************************************************************************************/

void get_rec_next2_HHe(int model, REC_COSMOPARAMS *cosmo, double z_in, double Tm,
                       double *xH1, double *xHeII, HYREC_ATOMIC *atomic, RADIATION *rad, unsigned iz_rad,
		       double dxHIIdlna_prev[2], double dxHeIIdlna_prev[2], double dEdtdV, int *stiff, int *error) {

  double dxHeIIdlna, dxHIIdlna, z_out, xe;
  double nH, H, TR;

  xe = *xHeII + 1.- (*xH1);
  nH = cosmo->nH0 *cube(1.+z_in);
  H  = rec_HubbleRate(cosmo, z_in);
  TR = kBoltz * cosmo->T0 *(1.+z_in);

  /* Evolve HeII by solving ODE */
  dxHeIIdlna  = rec_helium_dxHeIIdlna(*xH1, *xHeII, cosmo->nH0, cosmo->T0, cosmo->fHe, H, z_in, cosmo->fsR, cosmo->meR);
  *xHeII     += DLNA * hyrec_integrator(dxHeIIdlna, dxHeIIdlna_prev);

  /* Compute dxHII/dlna. This also correctly updates the radiation field at z_in,
     which is required even when using the stiff approximation */
  dxHIIdlna = rec_dxHIIdlna(model, xe, 1.-(*xH1), nH, H, kBoltz*Tm, TR, atomic, rad,
			    iz_rad, z_in, cosmo->fsR, cosmo->meR, dEdtdV, error, cosmo);

  /* If system is stiff use the quasi-equilibrium solution */

  if(*stiff == 1){
    z_out = (1.+z_in)*exp(-DLNA)-1.;
    rec_xH1_stiff(model, cosmo, z_out, *xHeII, xH1, atomic, rad, iz_rad+1, dEdtdV, stiff, error);
  }

  /* Otherwise use second-order explicit solver */
  else *xH1 -= DLNA * hyrec_integrator(dxHIIdlna, dxHIIdlna_prev);

  /* Checking for errors */
  if (*xH1 < 0. || *xH1 > 1. || *xH1 != *xH1) {
    fprintf(stderr, "\033[1m\033[31m error\033[22;30m in get_rec_next2_HHe: at z = %f, xH1 = %E\n", z_out, *xH1);
    *error = 1;
    return;
  }
  if (*error == 1) {
    fprintf(stderr, "  called from get_rec_next2_HHe at z = %f\n", z_out);
    return;
  }

}

/*********************************************************************************************************
Second-order integrator (unless stiff system) to evolve hydrogen only, assuming helium has entirely recombined.
Tm is given as an input (to avoid computing it twice) and fixed to quasi-equilibrium value with Tr.
Input : xe [at z_in]
Output: xe [at next timestep]
**********************************************************************************************************/

void rec_get_xe_next1_H(int model, REC_COSMOPARAMS *cosmo, double z_in, double xe_in, double Tm_in,
			double *xe_out, double *Tm_out, HYREC_ATOMIC *atomic, RADIATION *rad, unsigned iz_rad,
			double dxedlna_prev[2], double dEdtdV, int *stiff, int *error) {

  double dxedlna, z_out;
  double nH, H, TR, xH1;

  nH = cosmo->nH0 *cube(1.+z_in);
  H  = rec_HubbleRate(cosmo, z_in);
  TR = kBoltz *cosmo->T0 *(1.+z_in);


  /* Compute dxHII/dlna. This also correctly updates the radiation field at z_in,
     which is required even when using the stiff approximation */
  dxedlna = rec_dxHIIdlna(model, xe_in, xe_in, nH, H, kBoltz*Tm_in, TR, atomic,
			  rad, iz_rad, z_in, cosmo->fsR, cosmo->meR, dEdtdV, error, cosmo);

  z_out   = (1.+z_in)*exp(-DLNA)-1.;

  /* If system is stiff use the quasi-equilibrium solution */
  if (*stiff == 1) {
    xH1 = 1.-xe_in;
    rec_xH1_stiff(model, cosmo, z_out, 0, &xH1, atomic, rad, iz_rad+1, dEdtdV, stiff, error);
    *xe_out = 1.-xH1;
  }

  /* Otherwise use second-order explicit solver */
  else *xe_out = xe_in + DLNA * hyrec_integrator(dxedlna, dxedlna_prev);

  /* Quasi-steady state solution for Tm */
  *Tm_out = rec_Tmss(z_out, *xe_out, cosmo, dEdtdV);

  // Test that the outcome is sensible
  if (*xe_out > 1. || *xe_out < 0. || *xe_out != *xe_out) {
    fprintf(stderr, "\033[1m\033[31m error\033[22;30m in get_rec_next1_H at z = %E, xe = %E\n", z_out, *xe_out);
    *error = 1;
    return;
  }
  if (*error == 1) {
    fprintf(stderr, "  called from get_rec_next1_H at z = %f\n", z_out);
    return;
  }


}


/**********************************************************************************************
Second-order integrator for evolving xe and Tm simultaneously
Used for Hydrogen recombination only
May 2012: added a switch so Peebles model can be used at low redshift.
September 2016: added dEdtdV_dep, the *deposited* energy
(as opposed to dE_dtdV which is the instantaenously injected energy)
***********************************************************************************************/

void rec_get_xe_next2_HTm(int model, REC_COSMOPARAMS *cosmo,
			  double z_in, double xe_in, double Tm_in, double *xe_out, double *Tm_out,
                          HYREC_ATOMIC *atomic, RADIATION *rad, unsigned iz_rad,
			  double dxedlna_prev[2], double dTmdlna_prev[2], double dEdtdV, int *error) {

  double dxedlna, dTmdlna, nH, H, TR;

  nH = cosmo->nH0 *cube(1.+z_in);
  H  = rec_HubbleRate(cosmo, z_in);
  TR = kBoltz *cosmo->T0 *(1.+z_in);

  /*For low redshifts (z < 20 or so) use Peeble's model.
    The precise model does not metter much here as
    1) the free electron fraction is basically zero (~1e-4) in any case and
    2) the universe is going to be reionized around that epoch */

  if (TR/cosmo->fsR/cosmo->fsR/cosmo->meR <= TR_MIN
      || kBoltz*Tm_in/TR <= TM_TR_MIN) model = PEEBLES;

  dxedlna = rec_dxHIIdlna(model, xe_in, xe_in, nH, H, kBoltz*Tm_in, TR, atomic,
			  rad, iz_rad, z_in, cosmo->fsR, cosmo->meR, dEdtdV, error, cosmo);

  dTmdlna = rec_dTmdlna(z_in, xe_in, Tm_in, cosmo, dEdtdV);

  *xe_out = xe_in + DLNA *hyrec_integrator(dxedlna, dxedlna_prev);
  *Tm_out = Tm_in + DLNA *hyrec_integrator(dTmdlna, dTmdlna_prev);

  if (*error == 1) {
    fprintf(stderr, "  called from rec_get_xe_next2_HTm at z = %f\n", z_in);
    return;
  }

}

/****************************************************************************************************
Builds a recombination history with a given model
(model = PEEBLES, RECFAST, EMLA2s2p or FULL)
Added May 2012: The radiation field was added as an input so it can be extracted if desired
Replaced condition iz < param->nz by equivalent z >= 0 (so param.nz not needed anymore)
Added September 2016: follow the deposited energy dEdtdV_dep
****************************************************************************************************/

void rec_build_history(int model, double zstart, double zend,
		       REC_COSMOPARAMS *cosmo, HYREC_ATOMIC *atomic,
		       RADIATION *rad, double *xe_output, double *Tm_output) {

  long iz, iz_rad_0;
  double dxHIIdlna_prev[2], dTmdlna_prev[2], dxHeIIdlna_prev[2];
  double z, Delta_xe, xHeII, xH1, dEdtdV_dep, xe, nH, H;
  int quasi_eq;
  int error = 0;

  // Index at which we start integrating Hydrogen recombination, and corresponding redshift
  iz_rad_0  = (long) floor(1 + log(kBoltz*cosmo->T0/square(cosmo->fsR)/cosmo->meR*(1.+zstart)/TR_MAX)/DLNA);
  rad->z0   = (1.+zstart)*exp(-iz_rad_0 * DLNA) - 1.;


  z = zstart;

  /********* He III -> II Saha phase. Tm = Tr. Stop when xHeIII = 1e-8 *********/
  Delta_xe = cosmo->fHe;   /* Delta_xe = xHeIII here */

  for(iz = 0; z >= 0. && Delta_xe > 1e-8; iz++) {
    z = (1.+zstart)*exp(-DLNA*iz) - 1.;
    xe_output[iz] = rec_xesaha_HeII_III(cosmo->nH0, cosmo->T0, cosmo->fHe, z, &Delta_xe, cosmo->fsR, cosmo->meR);
    Tm_output[iz] = cosmo->T0 * (1.+z);
  }

  /******** He II -> I recombination.
	    Hydrogen in Saha equilibrium.
	    Tm fixed to steady state.
	    Neglect any energy injection.
	    Integrate until TR is low enough that can start integrating hydrogen recombination
	    (this occurs at index izH0 computed in rec_get_cosmoparam).
	    Start with quasi-equilibrium approximation.
  ********/
  dxHeIIdlna_prev[1] = dxHeIIdlna_prev[0] = 0.;

  xHeII    = rec_saha_xHeII(cosmo->nH0, cosmo->T0, cosmo->fHe, z, cosmo->fsR, cosmo->meR);
  quasi_eq = 1;                          /* Start with post-saha expansion */

  for(; iz <= iz_rad_0; iz++) {
    rec_get_xe_next1_He(cosmo, z, &xHeII, dxHeIIdlna_prev, &quasi_eq);
    z             = (1.+zstart)*exp(-DLNA*iz) - 1.;
    xH1           = rec_saha_xH1(xHeII, cosmo->nH0, cosmo->T0, z, cosmo->fsR, cosmo->meR);
    xe_output[iz] = 1.-xH1 + xHeII;
    Tm_output[iz] = rec_Tmss(z, xe_output[iz], cosmo, 0.);

  }

  /******** H II -> I and He II -> I simultaneous recombination (rarely needed but just in case)
	    Tm fixed to steady state.
	    Integrate H and He simultaneously until xHeII < XHEII_MIN
	    Start with post-saha expansion for hydrogen
	    Now account for possible energy injection.
	    Solve for dEdtdV_dep simultaneously;
  ********/

  dxHIIdlna_prev[1] = (xe_output[iz-2] - xe_output[iz-4])/2./DLNA - dxHeIIdlna_prev[1];
  dxHIIdlna_prev[0] = (xe_output[iz-1] - xe_output[iz-3])/2./DLNA - dxHeIIdlna_prev[0];
  quasi_eq          = 1;

  // Initialize energy *deposition*
  dEdtdV_dep = 0.;
  nH = cosmo->nH0*cube(1.+z);
  H  = rec_HubbleRate(cosmo, z);
  update_dEdtdV_dep(z, DLNA, xe_output[iz-1], Tm_output[iz-1], nH, H,
		    cosmo->inj_params, &dEdtdV_dep);
  for(; z >= 0. && xHeII > XHEII_MIN; iz++) {

    get_rec_next2_HHe(model, cosmo, z, Tm_output[iz-1], &xH1, &xHeII, atomic,
		      rad, iz-1-iz_rad_0, dxHIIdlna_prev, dxHeIIdlna_prev, dEdtdV_dep, &quasi_eq, &error);

    xe_output[iz] = 1.-xH1 + xHeII;

    z  = (1.+zstart)*exp(-DLNA*iz) - 1.;
    nH = cosmo->nH0*cube(1.+z);
    H  = rec_HubbleRate(cosmo, z);

    Tm_output[iz] = rec_Tmss(z, xe_output[iz], cosmo, dEdtdV_dep);

    update_dEdtdV_dep(z, DLNA, xe_output[iz], Tm_output[iz], nH, H, cosmo->inj_params, &dEdtdV_dep);

    if (error == 1) exit(1);
  }

  /******** H recombination. Helium assumed entirely neutral.
	    Tm fixed to steady-state until its relative difference from Tr is DLNT_MAX
  ********/

  for (; z >= 0. && fabs(1.-Tm_output[iz-1]/cosmo->T0/(1.+z)) < DLNT_MAX; iz++) {
    rec_get_xe_next1_H(model, cosmo, z, xe_output[iz-1], Tm_output[iz-1], xe_output+iz, Tm_output+iz,
		       atomic, rad, iz-1-iz_rad_0, dxHIIdlna_prev, dEdtdV_dep, &quasi_eq, &error);

    z  = (1.+zstart)*exp(-DLNA*iz) - 1.;
    nH = cosmo->nH0*cube(1.+z);
    H  = rec_HubbleRate(cosmo, z);
    update_dEdtdV_dep(z, DLNA, xe_output[iz], Tm_output[iz], nH, H, cosmo->inj_params, &dEdtdV_dep);

    if (error == 1) exit(1);
  }

  /******** Evolve xe and Tm simultaneously until z = zend
	    Note that the radiative transfer calculation is switched off automatically in the functions
	    rec_get_xe_next1_H and rec_get_xe_next2_HTm when it is no longer relevant.
  ********/

  dTmdlna_prev[1] = (Tm_output[iz-2] - Tm_output[iz-4])/2./DLNA;
  dTmdlna_prev[0] = (Tm_output[iz-1] - Tm_output[iz-3])/2./DLNA;

  for(; z > zend; iz++) {
    rec_get_xe_next2_HTm(model, cosmo, z, xe_output[iz-1], Tm_output[iz-1], xe_output+iz, Tm_output+iz,
			 atomic, rad, iz-1-iz_rad_0, dxHIIdlna_prev, dTmdlna_prev, dEdtdV_dep, &error);
    z  = (1.+zstart)*exp(-DLNA*iz) - 1.;
    nH = cosmo->nH0*cube(1.+z);
    H  = rec_HubbleRate(cosmo, z);
    update_dEdtdV_dep(z, DLNA, xe_output[iz], Tm_output[iz], nH, H, cosmo->inj_params, &dEdtdV_dep);

    if (error == 1) exit(1);
  }
}


/***********************************************************
Function to allocate and initialize HyRec internal tables
***********************************************************/

void hyrec_allocate(HYREC_DATA *data, double zmax, double zmin) {

  // data->zmax = (zmax > 3000.? zmax : 3000.);
  data->zmax = zmax;
  data->zmin = zmin;

  data->atomic = (HYREC_ATOMIC *) malloc(sizeof(HYREC_ATOMIC));
  allocate_and_read_atomic(data->atomic);
  data->cosmo  = (REC_COSMOPARAMS *) malloc(sizeof(REC_COSMOPARAMS));
  data->cosmo->inj_params = (INJ_PARAMS *)  malloc(sizeof(INJ_PARAMS));

  data->Nz = (long int) (log((1.+zmax)/(1.+zmin))/DLNA) + 2;
  data->xe_output = create_1D_array(data->Nz);
  data->Tm_output = create_1D_array(data->Nz);
  data->rad = (RADIATION *) malloc(sizeof(RADIATION));
  // For now assume that radiation field never needed over more than 1 decade in redshift
  // (typically from z ~ 1700 to 800 for recombination history)
  // Will have to adapt for outputting radiation fields at lower z
  allocate_radiation(data->rad, (long int) (log(10.)/DLNA));
}

void hyrec_free(HYREC_DATA *data) {
  free_atomic(data->atomic);
  free(data->cosmo->inj_params);
  free(data->cosmo);
  free(data->xe_output);
  free(data->Tm_output);
  free_radiation(data->rad);
  free(data->rad);
}

/******************************************************************
Compute a recombination history given input cosmological parameters
********************************************************************/

void hyrec_compute(HYREC_DATA *data, int model,
		   double h, double T0, double Omega_b, double Omega_m, double Omega_k, double YHe, double Nnueff,
		   double alphaR, double meR, double pann, double pann_halo, double ann_z, double ann_zmax,
		   double ann_zmin, double ann_var, double ann_z_halo, double Mpbh, double fpbh, int coll_ion, int on_the_spot){

  data->cosmo->T0    = T0;
  data->cosmo->orh2  = 4.48162687719e-7 *T0*T0*T0*T0 *(1. + 0.227107317660239 *Nnueff);
  data->cosmo->omh2  = Omega_m *h*h;
  data->cosmo->obh2  = Omega_b *h*h;
  data->cosmo->okh2  = Omega_k *h*h;
  data->cosmo->odeh2 = (1. -Omega_m -Omega_k - data->cosmo->orh2/h/h)*h*h;

  data->cosmo->nH0 = 11.223846333047e-6 *data->cosmo->obh2*(1.-YHe);  /* number density of hydrogen today in cm-3 */
  data->cosmo->fHe = YHe/(1.-YHe)/3.97153;                       /* fractional abundance of helium by number */
    /* these should depend on fsR and meR, strictly speaking; however, these are corrections to corrections */
  data->cosmo->fsR  = alphaR;
  data->cosmo->meR  = meR;

  /* DM annihilation parameters */
  data->cosmo->inj_params->odmh2      = data->cosmo->omh2 - data->cosmo->obh2;
  data->cosmo->inj_params->pann       = pann;
  data->cosmo->inj_params->pann_halo  = pann_halo;
  data->cosmo->inj_params->ann_z      = ann_z;
  data->cosmo->inj_params->ann_zmax   = ann_zmax;
  data->cosmo->inj_params->ann_zmin   = ann_zmin;
  data->cosmo->inj_params->ann_var    = ann_var;
  data->cosmo->inj_params->ann_z_halo = ann_z_halo;

  /* dark matter PBHs */
  data->cosmo->inj_params->Mpbh = Mpbh;
  data->cosmo->inj_params->fpbh = fpbh;
  data->cosmo->inj_params->coll_ion = coll_ion;

  /* energy deposition recipe */
  data->cosmo->inj_params->on_the_spot = on_the_spot;

  rec_build_history(model, data->zmax, data->zmin, data->cosmo, data->atomic,
		    data->rad, data->xe_output, data->Tm_output);

}

void hyrec_compute_CLASS(HYREC_DATA *data, int model){

  rec_build_history(model, data->zmax, data->zmin, data->cosmo, data->atomic,
		    data->rad, data->xe_output, data->Tm_output);

}
/*****
     Once HYREC_DATA outputs are computed, obtain xe(z) and Tm(z) by interpolation
*****/

double hyrec_xe(double z, HYREC_DATA *rec_data) {
  if (z > rec_data->zmax) return rec_data->xe_output[0];
  if (z < rec_data->zmin) {
    fprintf(stderr, "\033[1m\033[31m error\033[22;30m in hyrec_xe: requesting x_e at z = %f ", z);
    fprintf(stderr, "lower than zmin\n");
    exit(1);
  }
  return rec_interp1d(-log(1.+rec_data->zmax), DLNA, rec_data->xe_output, rec_data->Nz, -log(1.+z));
}

double hyrec_Tm(double z, HYREC_DATA *rec_data) {
  if(z > rec_data->zmax) return rec_data->cosmo->T0*(1.+z);
  if (z < rec_data->zmin) {
    fprintf(stderr, "\033[1m\033[31m error\033[22;30m in hyrec_Tm: requesting x_e at z = %f ", z);
    fprintf(stderr, "lower than zmin\n");
    exit(1);
  }
  return rec_interp1d(-log(1.+rec_data->zmax), DLNA, rec_data->Tm_output, rec_data->Nz, -log(1.+z));
}

double hyrec_dTmdlna(double z, HYREC_DATA *rec_data) {
  double dlna = 1e-2;
  double a = 1./(1.+z);
  double ahi = a*exp( dlna/2.);
  double alo = a*exp(-dlna/2.);

  if (1./alo-1. > rec_data->zmax) return -1.;
  if (1./ahi-1. < rec_data->zmin) return -2.;

  return  (rec_interp1d(-log(1.+rec_data->zmax), DLNA, rec_data->Tm_output, rec_data->Nz, log(ahi))
	  -rec_interp1d(-log(1.+rec_data->zmax), DLNA, rec_data->Tm_output, rec_data->Nz, log(alo)))/dlna;
}
