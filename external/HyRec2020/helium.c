/*************************************************************************************************/
/*                       HYREC-2: Hydrogen and Helium Recombination Code                         */
/*                  Written by Yacine Ali-Haimoud and Chris Hirata (2010-17)                     */
/*                      with contributions from Nanoom Lee (2020)                                */
/*                                                                                               */
/*         helium.c: all functions related to helium recombination (and Saha-equilibria)         */
/*                                                                                               */
/*         Units used in these modules: cm, s, Kelvin                                            */
/*                                                                                               */
/*         Revision history:                                                                     */
/*           - February 2015: changed length units from m to cm, so densities are in cm^-3       */
/*                            like in the rest of the code                                       */
/*           - May 2012: - added explicit dependence on the fine structure constant              */
/*                         and electron mass. Post-Saha functions moved to history.c             */
/*                       - dxHeII/dlna now includes explicit dependence on xH1s                  */
/*                    (in case H and He recombinations overlap so H is not in Saha equilibrium)  */
/*           - January 2011: added input variables so the value of xHeIII and                    */
/*                           post-Saha corrections can be monitored from external routines       */
/*           - written November 2010                                                             */
/*************************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "helium.h"

/************************************************************************************************
Gets xe in He II + III equilibrium
*************************************************************************************************/

double rec_xesaha_HeII_III(REC_COSMOPARAMS *cosmo, double z, double *xHeIII) {

  double Tr, nH, ainv, s;
  double nH0 = cosmo->nH0, Tr0 = cosmo->T0, fHe = cosmo->fHe;
  double fsR = cosmo->fsR, meR = cosmo->meR;
  /* Current conditions */
  Tr = Tr0*(ainv=1+z) /fsR/fsR/meR; /* Rescaled for different alpha, me */
  nH = nH0*ainv*ainv*ainv;

  /* Obtain Saha ratio, xe*xHeIII/xHeII
   * prefactor = (2*pi*m*k/h^2)^1.5. Includes dependence on fine-structure constant and electron mass*/

  s = fsR*fsR*fsR*meR*meR*meR* 2.414194e15*Tr*sqrt(Tr)*exp(-631462.7/Tr)/nH;

  *xHeIII = 2.*s*fHe/(1.+s+fHe)/(1.+sqrt(1.+4.*s*fHe/(1.+s+fHe)/(1.+s+fHe)));

  return(1. + fHe + *xHeIII);
}

/***************************************************************************************************
Gets xHeII (NOT xe) in Saha equilibrium, hydrogen assumed fully ionized.
****************************************************************************************************/

double rec_saha_xHeII(REC_COSMOPARAMS *cosmo, double z) {

  double Tr, nH, ainv, s;
  double nH0 = cosmo->nH0, Tr0 = cosmo->T0, fHe = cosmo->fHe;
  double fsR = cosmo->fsR, meR = cosmo->meR;

  /* Current conditions */
  Tr = Tr0*(ainv=1+z) /fsR/fsR/meR; /* Rescaled for different alpha, me */
  nH = nH0*ainv*ainv*ainv;

  /* Obtain Saha ratio, xe*xHeII/xHeI. prefactor = (2*pi*m*k/h^2)^1.5 */
  s = fsR*fsR*fsR*meR*meR*meR *2.414194e15*Tr*sqrt(Tr)*exp(-285325./Tr)/nH * 4.;

  return 2.*s*fHe/(1.+s)/(1.+sqrt(1.+4.*s*fHe/(1.+s)/(1.+s)));

}

/*******************************************************************************************
Gets xH1s for Saha equilibrium HII + e <-> H(1s) + gamma, where xe = xHII + xHeII
(in case Helium has not full recombined). Added May 2012.
*******************************************************************************************/

double rec_saha_xH1s(REC_COSMOPARAMS *cosmo, double z, double xHeII) {

  double s, Tr, nH, ainv;
  double nH0 = cosmo->nH0, T0 = cosmo->T0;
  double fsR = cosmo->fsR, meR = cosmo->meR;

  Tr = T0 * (ainv=1.+z) /fsR/fsR/meR; /* current (rescaled) radiation temperature */
  nH = nH0*ainv*ainv*ainv;

  /* Obtain Saha ratio, xe xHII/xH1s, with xe = xHII + xHeII and xHII = 1-xH1s
   * Constants:
   * reduced mass Rydberg = 157801.37882 Kelvin (no relativistic corrections)
   * prefactor = (2*pi*m*k/h^2)^1.5  */

  s = fsR*fsR*fsR*meR*meR*meR *2.4127161187130e15*Tr*sqrt(Tr)*exp(-157801.37882/Tr)/nH;
  if      (s == 0.) return 1.;
  else if (s > 1e5) return (1.+xHeII)/s - (xHeII*xHeII + 3.*xHeII + 2.)/s/s;  /* Second-order value for xHII ~ 1 */
  else              return 1.-2./(1.+ xHeII/s + sqrt((1.+ xHeII/s)*(1.+ xHeII/s) +4./s));

}

/*************************************************************************************
 He II->I recombination rate evolution equation.
 Incorporates 2(1)P-->1(1)S escape, with H continuum opacity
              2(1)S-->1(1)S 2 photon decay
              2(3)P-->1(1)S Sobolev
 Assumes Compton temperature equilibrium, no Thomson opacity.
 Added May 2012: Explicit dependence on xH1s, which may differ from its Saha equilibrium
value if hydrogen has already started to recombine (if H and He recombinations overlap).
**************************************************************************************/

double rec_helium_dxHeIIdlna(HYREC_DATA *data, double z, double xH1, double xHeII, double H) {
  double nH0 = data->cosmo->nH0, Tr0 = data->cosmo->T0, fHe = data->cosmo->fHe;
  double fsR = data->cosmo->fsR, meR = data->cosmo->meR;
  int *error = &data->error;
  double Tr, nH, ainv, s, s0;
  double xe, xHeI, y2s, y2p, ydown;
  double etacinv, g2pinc, dnuline, tau2p, pesc, tauc, enh;
  char sub_message[256];
  if (*error == 1) return 0.;

  xe   = xHeII + (1.-xH1);
  xHeI = fHe - xHeII;
  /* First catch potential errors */
  if (xHeI < 0.){
    sprintf(sub_message, "xHeI = %E < 0 at z = %f in rec_helium_dxHeIIdlna.\nYou should try and extend the hydrogen post-saha phase by increasing DXHII_MAX in history.h\n", xHeI, z);
    strcat(data->error_message, sub_message);
    *error = 1;
  }

  /* Current conditions */
  Tr = Tr0*(ainv=1+z) /fsR/fsR/meR;  /* Rescaled temperature to account for different alpha or me */
  nH = nH0*ainv*ainv*ainv;

  /* Saha abundance ratio, and ratio for zero ionization potential */
  s0 = fsR*fsR*fsR*meR*meR*meR *2.414194e15*Tr*sqrt(Tr)/nH * 4.;
  s = s0 * exp(-285325./Tr);

  /* Abundances of excited levels.  y is defined as x_i = y_i * xe * xHeII */
  y2s = exp(46090./Tr)/s0;
  y2p = exp(39101./Tr)/s0 * 3.;

  /* Continuum opacity */
  /* coeff is 1/(sigma lambda) in cm^(-3). Result is in s^{-1}  */

  etacinv = fsR*fsR*fsR *meR*meR*meR* 9.15776e22 * H/nH/xH1;

  /* Incoherent width of 2(1)P incl. ->2s, 3s, 3d, 4s, 4d, 5s, 5d */
  g2pinc = 1.976e6 / (1.-exp(-6989./Tr))
          + 6.03e6 / (exp(19754./Tr)-1.)
          + 1.06e8 / (exp(21539./Tr)-1.)
          + 2.18e6 / (exp(28496./Tr)-1.)
          + 3.37e7 / (exp(29224./Tr)-1.)
          + 1.04e6 / (exp(32414./Tr)-1.)
          + 1.51e7 / (exp(32781./Tr)-1.);

  g2pinc *= fsR*fsR*fsR*fsR*fsR*meR ;

  /* Optical depth, escape prob in x2p */
  tau2p   = 4.277e-8 * nH/H * xHeI /(fsR*meR*meR);
  dnuline = g2pinc * tau2p / (4.*M_PI*M_PI);
  tauc    = dnuline/etacinv;
  enh     = sqrt(1.+M_PI*M_PI*tauc) + 7.74*tauc/(1.+70.*tauc);
  pesc    = enh / tau2p;


  /* Effective increase in escape probability via intercombination line
    * ratio of optical depth to allowed line = 1.023e-7
    * 1-e^-tau23 = absorption prob. in intercom line
    * e^[(E21P-E23P)/T] - e^[-(E21P-E23P)*eta_c] = step in intercom line
    *   relative to N_line in 21P
    * divide by tau2p to get effective increase in escape prob.
    * factor of 0.964525 is phase space factor for intercom vs allowed line -- (584/591)^3
  */
  pesc = pesc + (1.-exp(-1.023e-7*tau2p))*(0.964525*exp(2947./Tr)-enh*exp(-6.14e13*fsR*fsR*meR/etacinv))/tau2p;

  /* Total decay rate */

  ydown = fsR*fsR*fsR*fsR*fsR*fsR*fsR*fsR*meR *50.94*y2s  /*prefactor correcting 2-photon decay rate given alpha, me*/
        + fsR*fsR*fsR*fsR*fsR*meR* 1.7989e9*y2p*pesc;     /*prefactor correcting 1-photon dipole rate */


  return  ydown*(xHeI*s - xHeII*xe)/H;   /* Excitation is obtained by detailed balance */

}

/***********************************************************************************************************/


