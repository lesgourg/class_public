/** 
 * definitions for module hyperspherical.c 
 */

#ifndef __HYPERSPHERICAL__
#define __HYPERSPHERICAL__

#include "common.h"
#define _HYPER_OVERFLOW_ 1e50
#define _HYPER_SAFETY_ 1e-5
#define _TRIG_PRECISSION_ 1e-7

typedef struct HypersphericalInterpolationStructure{
  int K;                 //Sign of the curvature, (0,-1,1)
  double beta;
  double deltax;         //x-spacing. (xvec is uniformly spaced)
  int trig_order;        //Order of the interpolation formula for SinK and CosK.
  int nl;                //Number of l values
  int *lvec;             //Vector of l values stored
  int nx;                //Number of x-values  
  double *xvec;          //Pointer to x-values
  double *sinK;          //Vector of sin_K(xvec)
  double *cotK;          //Vector of cot_K(xvec)
  double *phivec;        //array of size nl*nx. [y_{l1}(x1) t_{l1}(x2)...]
  double *dphivec;       //Same as phivec, but containing derivatives.
} HyperInterpStruct;


/**
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif
int hyperspherical_HIS_create(int K, 
                              double beta, 
                              int nl, 
                              int *lvec, 
                              double xmin, 
                              double xmax, 
                              double sampling,
                              HyperInterpStruct **ppHIS, 
                              ErrorMsg error_message);
  int hyperspherical_HIS_free(HyperInterpStruct *pHIS);
  int hyperspherical_forwards_recurrence(int K, 
                                         int lmax, 
                                         double beta, 
                                         double x, 
                                         double sinK,
                                         double cotK,
                                         double *sqrtK,
                                         double *PhiL);
  int hyperspherical_backwards_recurrence(int K, 
                                          int lmax, 
                                          double beta, 
                                          double x, 
                                          double sinK,
                                          double cotK,
                                          double *sqrtK,
                                          double *PhiL);
  int hyperspherical_Hermite_interpolation_vector(HyperInterpStruct *pHIS,
                                                  int nxi,
                                                  int lnum,
                                                  double *xinterp,
                                                  double *Phi,
                                                  double *dPhi,
                                                  double *d2Phi,
                                                  double *sinKinterp,
                                                  double *cosKinterp);
                                                  int hyperspherical_get_xmin(HyperInterpStruct *pHIS,
                            double xtol,
                            double phiminabs,
                            double *xmin);
  int hyperspherical_WKB(int K,int l,double beta,double y, double *Phi);
  int get_CF1(int K,int l,double beta, double cotK, double *CF, int *isign);
  int CF1_from_Gegenbauer(int l, int beta, double sinK, double cotK, double *CF);
  double airy_cheb_approx(double z);
  double coef1(double z);
  double coef2(double z);
  double coef3(double z);
  double coef4(double z);
  double cheb(double x, int n, double *A);
  double get_value_at_small_phi(int K,int l,double beta,double Phi);
#ifdef __cplusplus
}
#endif

#endif
