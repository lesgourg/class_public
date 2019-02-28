/*************************************************************************************************/
/*                       HYREC: Hydrogen and Helium Recombination Code                           */
/*                  Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                     */
/*                                                                                               */
/*         helium.c: all functions related to helium recombination (and Saha-equilibria)         */
/*                                                                                               */
/*         Units used in these modules: m, s, Kelvin                                             */
/*                                                                                               */
/*         Version: May 2012                                                                     */
/*                                                                                               */
/*         Revision history:                                                                     */
/*           - written November 2010                                                             */
/*           - January 2011: added input variables so the value of xHeIII                        */ 
/*           and post-Saha corrections can be monitored from external routines                   */ 
/*           - May 2012: - added explicit dependence on the fine structure constant              */
/*                         and electron mass. Post-Saha functions moved to history.c             */
/*                       - dxHeII/dlna now includes explicit dependence on xH1s                  */
/*                    (in case H and He recombinations overlap so H is not in Saha equilibrium)  */
/*************************************************************************************************/ 

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/************************************************************************************************
Gets xe in He II + III equilibrium 
*************************************************************************************************/

double rec_xesaha_HeII_III(double nH0, double Tr0, double fHe, double z, double *xHeIII, 
                           double fsR, double meR) {

  double Tr, nH, ainv, s;

  /* Current conditions */
  Tr = Tr0*(ainv=1+z) /fsR/fsR/meR; /* Rescaled for different alpha, me */
  nH = nH0*ainv*ainv*ainv;

  /* Obtain Saha ratio, xe*xHeIII/xHeII
   * prefactor = (2*pi*m*k/h^2)^1.5. Includes dependence on fine-structure constant and electron mass*/

  s = fsR*fsR*fsR*meR*meR*meR* 2.414194e21*Tr*sqrt(Tr)*exp(-631462.7/Tr)/nH;

  *xHeIII = 2.*s*fHe/(1.+s+fHe)/(1.+sqrt(1.+4.*s*fHe/(1.+s+fHe)/(1.+s+fHe)));  

  return(1. + fHe + *xHeIII);
}

/***************************************************************************************************
Gets xHeII (NOT xe) in Saha equilibrium, hydrogen assumed fully ionized. 
****************************************************************************************************/

double rec_saha_xHeII(double nH0, double Tr0, double fHe, double z, double fsR, double meR) {

   double Tr, nH, ainv, s;

   /* Current conditions */
   Tr = Tr0*(ainv=1+z) /fsR/fsR/meR; /* Rescaled for different alpha, me */
   nH = nH0*ainv*ainv*ainv;

   /* Obtain Saha ratio, xe*xHeII/xHeI. prefactor = (2*pi*m*k/h^2)^1.5 */
   s = fsR*fsR*fsR*meR*meR*meR *2.414194e21*Tr*sqrt(Tr)*exp(-285325./Tr)/nH * 4.;
  
  return 2.*s*fHe/(1.+s)/(1.+sqrt(1.+4.*s*fHe/(1.+s)/(1.+s)));  

}

/*******************************************************************************************
Gets xH1s for Saha equilibrium HII + e <-> H(1s) + gamma, where xe = xHII + xHeII 
(in case Helium has not full recombined). Added May 2012.
*******************************************************************************************/

double rec_saha_xH1s(double xHeII, double nH0, double T0, double z, double fsR, double meR) {

   double s, Tr, nH, ainv;

   Tr = T0 * (ainv=1.+z) /fsR/fsR/meR; /* current (rescaled) radiation temperature */
   nH = nH0*ainv*ainv*ainv;

   /* Obtain Saha ratio, xe xHII/xH1s, with xe = xHII + xHeII and xHII = 1-xH1s
    * Constants:
    * reduced mass Rydberg = 157801.37882 Kelvin (no relativistic corrections)
    * prefactor = (2*pi*m*k/h^2)^1.5  */

   s = fsR*fsR*fsR*meR*meR*meR *2.4127161187130e21*Tr*sqrt(Tr)*exp(-157801.37882/Tr)/nH;
  
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

double rec_helium_dxHeIIdlna(double xH1s, double xHeII, double nH0, double Tr0, double fHe, double H, double z,
                             double fsR, double meR) {
    double Tr, nH, ainv, s, s0;
    double xe, xHeI, y2s, y2p, ydown;
    double etacinv, g2pinc, dnuline, tau2p, pesc, tauc, enh;
  

    xe   = xHeII + (1.-xH1s);
    xHeI = fHe - xHeII;

    /* First catch potential errors */
    if (xHeI < 0.){
       fprintf(stderr, "Error in rec_helium_dxHeIIdlna, xHeI = %E < 0 at z = %f .\n", xHeI, z);
       fprintf(stderr, "You should try and extend the hydrogen post-saha phase by increasing DXHII_MAX in history.h\n");
       exit(1);
    }
  
    /* Current conditions */
    Tr = Tr0*(ainv=1+z) /fsR/fsR/meR;  /* Rescaled temperature to account for different alpha or me */
    nH = nH0*ainv*ainv*ainv;

    /* Saha abundance ratio, and ratio for zero ionization potential */
    s0 = fsR*fsR*fsR*meR*meR*meR *2.414194e21*Tr*sqrt(Tr)/nH * 4.;
    s = s0 * exp(-285325./Tr);

    /* Abundances of excited levels.  y is defined as x_i = y_i * xe * xHeII */
    y2s = exp(46090./Tr)/s0;
    y2p = exp(39101./Tr)/s0 * 3.;

    /* Continuum opacity */
    /* coeff is 1/(sigma lambda) in m^(-3). Result is in s^{-1}  */
  
    etacinv = fsR*fsR*fsR *meR*meR*meR* 9.15776e28 * H/nH/xH1s;   

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
    tau2p   = 4.277e-14 * nH/H * xHeI /(fsR*meR*meR);
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
