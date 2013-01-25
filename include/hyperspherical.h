/** 
 * definitions for module hyperspherical.c 
 */

#ifndef __HYPERSPHERICAL__
#define __HYPERSPHERICAL__

#include "common.h"


/**
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int HypersphericalClosedGegenbauer(int l, int beta, double y, double *Phi,ErrorMsg errmsg);
  int HypersphericalOpenRecurrence(int l, double beta, double x, double *Phi_l,ErrorMsg errmsg);
  int HypersphericalWKB(int K,int l,double beta,double y, double *Phi,ErrorMsg errmsg);
  int HypersphericalExplicit(int K,int l, double beta,double x, double *Phi,ErrorMsg errmsg);

  double Gegenbauer(int n,int alpha, double x);
  double Hypfactorial(int n);
  int get_CF1(int K,int l,double beta,double x,double *CF,int *isign);
  int HypersphericalForward(int K,int l, double beta,double x, double *Phi_l);
  int HypersphericalBackward(int K,int l,double beta,double x, double *Phi_l);
  double airy_cheb_approx(double z);
  double coef1(double z);
  double coef2(double z);
  double coef3(double z);
  double coef4(double z);
  double cheb(double x, int n, double *A);
  double get_value_at_small_phi(int K,int l,double beta,double Phi);
  int ClosedModY(int l, int beta, double *y);
#ifdef __cplusplus
}
#endif

#endif
