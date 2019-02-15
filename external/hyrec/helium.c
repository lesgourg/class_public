/*************************************************************************************************/
/*                 HYREC: Hydrogen and Helium Recombination Code                                 */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                              */
/*                                                                                               */
/*         helium.c: all functions related to Hydrogen recombination                             */
/*                                                                                               */
/*         Units used in these modules: m, s, Kelvin                                             */
/*                                                                                               */
/*         Version: January 2011                                                                 */
/*         Revision history:                                                                     */
/*           - written November 2010                                                             */
/*           - January 2011: added input variables so the value of xHeIII                        */ 
/*           and post-Saha corrections can be monitored from external routines                   */ 
/*************************************************************************************************/ 


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/************************************************************************************************
Gets xe in He II + III equilibrium 
*************************************************************************************************/

double rec_sahaHeII(double nH0, double Tr0, double fHe, double z, double *xHeIII) {

  double Tr, nH, ainv, s;

  /* Current conditions */
  Tr = Tr0*(ainv=1+z);
  nH = nH0*ainv*ainv*ainv;

  /* Obtain Saha ratio, xe*xHeIII/xHeII
   * prefactor = (2*pi*m*k/h^2)^1.5 */

  s = 2.414194e21*Tr*sqrt(Tr)*exp(-631462.7/Tr)/nH;

  *xHeIII = 2.*s*fHe/(1.+s+fHe)/(1.+sqrt(1.+4.*s*fHe/(1.+s+fHe)/(1.+s+fHe)));  

  return(1. + fHe + *xHeIII);
}

/**************************************************************************************
Gets xe in He I + II equilibrium 
***************************************************************************************/

double rec_sahaHeI(double nH0, double Tr0, double fHe, double z) {

  double Tr, nH, ainv, s, q;

  /* Current conditions */
  Tr = Tr0*(ainv=1+z);
  nH = nH0*ainv*ainv*ainv;

  /* Obtain Saha ratio, xe*xHeII/xHeI
   * prefactor = (2*pi*m*k/h^2)^1.5 */

  s = 2.414194e21*Tr*sqrt(Tr)*exp(-285325./Tr)/nH * 4.;
  
  q = 2.*s*fHe/(1.+s)/(1.+sqrt(1.+4.*s*fHe/(1.+s)/(1.+s)));  /* q = xHeII */

  return(1.+q);
}

/*************************************************************************************** 
Saha equilibrium for Hydrogen, needed in the evolution equation for helium 
***************************************************************************************/

double rec_saha_xe_H(double nH0, double T0, double z) {

  double xeSaha, s, Tr, nH, ainv;

  Tr = T0 * (ainv=1.+z); /* current radiation temperature */
  nH = nH0*ainv*ainv*ainv;

  /* Obtain Saha ratio, xe^2/(1-xe)
   * Constants:
   * reduced mass Rydberg = 157801.37882 Kelvin (no relativistic corrections)
   * prefactor = (2*pi*m*k/h^2)^1.5  */

  s = 2.4127161187130e21*Tr*sqrt(Tr)*exp(-157801.37882/Tr)/nH;

  xeSaha = 2./(1.+sqrt(1.+4./s));
  return(xeSaha);
}

/*************************************************************************************
 He II->I recombination rate evolution equation.
 Incorporates 2(1)P-->1(1)S escape, with H continuum opacity
              2(1)S-->1(1)S 2 photon decay
              2(3)P-->1(1)S Sobolev
 Assumes Compton temperature equilibrium, no Thomson opacity.
**************************************************************************************/

double rec_helium_dxedt(double xe, double nH0, double Tr0, double fHe, double H, double z) {
  double Tr, nH, ainv, s, s0;
  double xHeI, xHeII, y2s, y2p, xHII;
  double etacinv, g2pinc, dnuline, tau2p, pesc, tauc;
  double xdown, xup, ydown, enh;

  /* Current conditions */
  Tr = Tr0*(ainv=1+z);
  nH = nH0*ainv*ainv*ainv;

  /* Saha abundance ratio, and ratio for zero ionization potential */
  s0 = 2.414194e21*Tr*sqrt(Tr)/nH * 4.;
  s = s0 * exp(-285325./Tr);

  /* Abundances of excited levels.  y is defined as x_i = y_i * xe * xHeII */
  xHII = rec_saha_xe_H(nH0,Tr0,z);
  xHeII = xe - xHII;
  xHeI = fHe - xHeII;
  y2s = exp(46090./Tr)/s0;
  y2p = exp(39101./Tr)/s0 * 3.;

  /* Continuum opacity = H/nH^2 * T^3/2 * exp(-chi/kT) * (2 pi m_e k/h^2)^{3/2} nu/(sigma c) */
  /* coef. is 2.205e50 => log = 115.920                                                      */
  etacinv = H/(nH*nH*xe) * Tr*sqrt(Tr) * exp(115.920-157801.37882/Tr);

  /* Incoherent width of 2(1)P incl. ->2s, 3s, 3d, 4s, 4d, 5s, 5d */
  g2pinc = 1.976e6 / (1.-exp(-6989./Tr))
          + 6.03e6 / (exp(19754./Tr)-1.)
          + 1.06e8 / (exp(21539./Tr)-1.)
          + 2.18e6 / (exp(28496./Tr)-1.)
          + 3.37e7 / (exp(29224./Tr)-1.) 
          + 1.04e6 / (exp(32414./Tr)-1.)
          + 1.51e7 / (exp(32781./Tr)-1.);


  /* Optical depth, escape prob in x2p */
  tau2p = 4.277e-14 * nH/H * xHeI;
  dnuline = g2pinc * tau2p / (4.*M_PI*M_PI);
  tauc = dnuline/etacinv;
  enh = sqrt(1.+M_PI*M_PI*tauc) + 7.74*tauc/(1.+70.*tauc);
  pesc = enh / tau2p;

  
  /* Effective increase in escape probability via intercombination line
    * ratio of optical depth to allowed line = 1.023e-7
    * 1-e^-tau23 = absorption prob. in intercom line
    * e^[(E21P-E23P)/T] - e^[-(E21P-E23P)*eta_c] = step in intercom line
    *   relative to N_line in 21P
    * divide by tau2p to get effective increase in escape prob.
    * factor of 0.964525 is phase space factor for intercom vs allowed line -- (584/591)^3
  */
   pesc = pesc + (1.-exp(-1.023e-7*tau2p))*(0.964525*exp(2947./Tr)-enh*exp(-6.14e13/etacinv))/tau2p;

   /* Net decay rate */
   ydown = 50.94*y2s + 1.7989e9*y2p*pesc;


  /* Excitation via detailed balance */
  xdown = ydown * xHeII * xe;
  xup = ydown * xHeI * s;

  /* Return recombination rate, including derivative of hydrogen Saha */
  return(xup-xdown + H*(1.+z)*(rec_saha_xe_H(nH0,Tr0,z-0.5)-rec_saha_xe_H(nH0,Tr0,z+0.5)) );
}

/******************************************************************************************************
Post-Saha expansion for beginning of Helium II -> Helium I
*******************************************************************************************************/

double xe_PostSahaHe(double nH0, double Tr0, double fHe, double H, double z, double *Delta_xe){
 
   double Tr, nH, s, xeSaha, dxeSahadt, Dxe, DxedotDxe, ainv;

   Tr = Tr0 * (ainv=1.+z); /* current radiation temperature */
   nH = nH0*ainv*ainv*ainv;

   s = 2.414194e21*Tr*sqrt(Tr)*exp(-285325./Tr)/nH * 4.;
   xeSaha = rec_sahaHeI(nH0, Tr0, fHe, z);
   dxeSahadt = - xeSaha*(xeSaha-1.)/(2.*xeSaha +s-1.) * (285325./Tr - 1.5) * H;  /* analytic derivative of xeSaha */
   
   Dxe = 0.01 * (1.+fHe-xeSaha);
   DxedotDxe = (rec_helium_dxedt(xeSaha + Dxe, nH0, Tr0, fHe, H, z)
		- rec_helium_dxedt(xeSaha - Dxe, nH0, Tr0, fHe, H, z))/2./Dxe;
   
   *Delta_xe = dxeSahadt/DxedotDxe;       

   return  xeSaha + *Delta_xe;

}

/********************************************************************************************************/
