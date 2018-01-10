/***********************************************************************
 High-redshift 21 cm module for Hyrec.

Computes the coefficients of the expansion (in Kelvins)
T_b = <T_b> + a \delta_H + b \delta_Tgas 
      + c \delta_H^2 + d \delta_Tgas^2 + e \delta_H \delta_Tgas,
where \delta_H = \delta n_H/n_H, \delta_Tgas = \delta(Tgas)/Tgas.
 
Does not yet account for Wouthuysen-Field coupling.
Written December 2014.
***********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "hyrectools.h"

#define A10 2.86889e-15      /* Spontaneous decay rate of the 21 cm transtion */
#define T10 0.0681683        /* Transition energy in Kelvins */
#define lambda21 21.106114   /* 21cm wavelength */


/* Tables for collisional coefficients.
   H-H and e-H rates from astro-ph/0608032 */

double logT_HH[27] = {0, 0.693, 1.39, 1.79, 2.08, 2.30, 2.71, 3.00, 3.22,
		      3.40, 3.69, 3.91, 4.09, 4.25, 4.38, 4.50, 4.61, 5.30,
		      5.70, 6.21, 6.55, 6.91, 7.60, 8.01, 8.52, 8.85, 9.21};

double logk_HH[27] = {-29.6115, -29.5759, -28.9367, -28.0465, -27.2458, -26.5732, 
                      -25.4227, -24.7518, -24.3241, -24.0282, -23.6457, -23.4027, -23.2316, 
                      -23.1038, -23.006,  -22.9215, -22.8519, -22.4662, -22.2887, -22.0858, 
		      -21.9577, -21.8289, -21.5742, -21.4224, -21.2291, -21.0987, -20.9628};
   
double logT_eH[17] = {0, 0.693, 1.61, 2.30, 3.00, 3.91, 4.61, 5.30, 6.21,
		      6.91, 7.60, 8.01, 8.52, 8.85, 9.21, 9.62, 9.90};

double logk_eH[17] = {-22.1546, -21.8109, -21.3581, -21.0163, -20.6745, -20.2347,
		      -19.9079, -19.5886, -19.1994, -18.9449, -18.7562, -18.6807,
		      -18.6228, -18.6046, -18.5986, -18.6082, -18.6302};

/* Interpolating functions for the collisional couplings */
/* double kappa10_H(double Tgas) { */
/*   if (log(Tgas) <= logT_HH[0]) return exp(logk_HH[0]); */
/*   if (log(Tgas) >= logT_HH[26]) return exp(logk_HH[26]); */
/*   return exp(rec_interpol_G(log(Tgas), logT_HH, logk_HH, 27)); */
/* } */
/* double kappa10_e(double Tgas) { */
/*   if (log(Tgas) <= logT_eH[0]) return exp(logk_eH[0]); */
/*   if (log(Tgas) >= logT_eH[16]) return exp(logk_eH[16]); */
/*   return exp(rec_interpol_G(log(Tgas), logT_eH, logk_eH, 17)); */
/* } */

double kappa10_p(double Tgas) {
  return 0.;
}

/* Simple fitting function as an alternative (works better for perturbations) */
double kappa10_H(double Tgas) {
  return 3.1e-11 *exp(0.357*log(Tgas) - 32./Tgas);
}

double kappa10_e(double Tgas){
  return 0.;
}


/* Spin temperature minus Tcmb */
double Tspin_Tcmb(double xe, double Tgas, double Tcmb, double nh) {
  double n_kappa10; 
  n_kappa10 = nh *((1.-xe) *kappa10_H(Tgas) 
		   + xe *kappa10_e(Tgas) + xe *kappa10_p(Tgas)); 

  return (Tgas - Tcmb)/(1. + A10/n_kappa10 *Tgas/T10);
}

/* Spin temperature */
double Tspin(double xe, double Tgas, double Tcmb, double nh) {
  return Tspin_Tcmb(xe, Tgas, Tcmb, nh) + Tcmb;   
}

/* Optical depth */
double tau21(double xe, double Tgas, double Tcmb, double nh, double H_invsec) {
  return 3.*T10 /(32.* M_PI*(Tcmb + Tspin_Tcmb(xe, Tgas, Tcmb, nh)))
          *lambda21*lambda21*lambda21 *(1.-xe) *nh *A10 /H_invsec;
}

/* 21 cm mean brightness temperature: tau_21*(Ts - Tcmb)/(1+z) */
double T21cm(double a, double xe, double Tgas, double Tcmb, double nh, double H_invsec) {
    return  tau21(xe, Tgas, Tcmb, nh, H_invsec) *Tspin_Tcmb(xe, Tgas, Tcmb, nh)*a;
}

/* d T21/d\delta_H */
double T21cm_dh(double a, double xe, double Tgas, double Tcmb, double nh, double H_invsec) {
   double eps = 1e-3;
   return ( T21cm(a, xe, Tgas, Tcmb, nh*(1.+eps), H_invsec) 
	    -T21cm(a, xe, Tgas, Tcmb, nh*(1.-eps), H_invsec))/(2.*eps);
}

/* 1/2 d^2 T21/(d\delta_H)^2 */
double T21cm_dh2(double a, double xe, double Tgas, double Tcmb, double nh, double H_invsec) {
   double eps = 1e-3;
   return (   T21cm(a, xe, Tgas, Tcmb, nh*(1.+eps), H_invsec) 
	    + T21cm(a, xe, Tgas, Tcmb, nh*(1.-eps), H_invsec)
            - 2.*T21cm(a, xe, Tgas, Tcmb, nh, H_invsec) )/(2.*eps*eps);
}

/* d T21/d\delta_Tgas */
double T21cm_dtgas(double a, double xe, double Tgas, double Tcmb, double nh, double H_invsec) {
   double eps = 1e-3;
   return ( T21cm(a, xe, Tgas*(1.+eps), Tcmb, nh, H_invsec) 
	    -T21cm(a, xe, Tgas*(1.-eps), Tcmb, nh, H_invsec))/(2.*eps);
}

/* 1/2 d^2 T21/(d\delta_Tgas)^2 */
double T21cm_dtgas2(double a, double xe, double Tgas, double Tcmb, double nh, double H_invsec) {
   double eps = 1e-3;
   return (   T21cm(a, xe, Tgas*(1.+eps), Tcmb, nh, H_invsec) 
	    + T21cm(a, xe, Tgas*(1.-eps), Tcmb, nh, H_invsec)
            - 2.*T21cm(a, xe, Tgas, Tcmb, nh, H_invsec) )/(2.*eps*eps);
}

/* d^2 T21/ (d delta_H  d delta_Tgas ) */
double T21cm_dhdtgas(double a, double xe, double Tgas, double Tcmb, double nh, double H_invsec) {
   double eps = 1e-3;
   return (   T21cm(a, xe, Tgas*(1.+eps), Tcmb, nh*(1.+eps), H_invsec) 
	    + T21cm(a, xe, Tgas*(1.-eps), Tcmb, nh*(1.-eps), H_invsec)
            + 2.*T21cm(a, xe, Tgas, Tcmb, nh, H_invsec) 
            - T21cm(a, xe, Tgas, Tcmb, nh*(1.+eps), H_invsec) 
	    - T21cm(a, xe, Tgas, Tcmb, nh*(1.-eps), H_invsec)
            - T21cm(a, xe, Tgas*(1.+eps), Tcmb, nh, H_invsec) 
	    - T21cm(a, xe, Tgas*(1.-eps), Tcmb, nh, H_invsec))/(2.*eps*eps);
}

/* */



