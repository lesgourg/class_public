/**
 * definitions for module hyperspherical.c
 */

#ifndef __HYPERSPHERICAL__
#define __HYPERSPHERICAL__

#include "common.h"
#define _HYPER_OVERFLOW_ 1e200
#define _ONE_OVER_HYPER_OVERFLOW_ 1e-200
#define _HYPER_SAFETY_ 1e-5
#define _TRIG_PRECISSION_ 1e-7
#define _HYPER_BLOCK_ 8
#define _HYPER_CHUNK_ 16
#define _TWO_OVER_THREE_ 0.666666666666666666666666666667e0
#define _HIS_BYTE_ALIGNMENT_ 16

typedef struct HypersphericalInterpolationStructure{
  int K;                 //Sign of the curvature, (0,-1,1)
  double beta;
  double delta_x;         //x-spacing. (xvec is uniformly spaced)
  int trig_order;        //Order of the interpolation formula for SinK and CosK.
  int l_size;                //Number of l values
  int *l;             //Vector of l values stored
  double * chi_at_phimin;     // vector x_min[index-l] below which neglect Bessels
  int x_size;                //Number of x-values
  double *x;          //Pointer to x-values
  double *sinK;          //Vector of sin_K(xvec)
  double *cotK;          //Vector of cot_K(xvec)
  double *phi;        //array of size nl*nx. [y_{l1}(x1) t_{l1}(x2)...]
  double *dphi;       //Same as phivec, but containing derivatives.
} HyperInterpStruct;

struct WKB_parameters{
   int K;
   int l;
   double beta;
   double phiminabs;
};

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
                                int l_WKB,
                                double phiminabs,
                                HyperInterpStruct *pHIS,
                                ErrorMsg error_message);

  int hyperspherical_HIS_free(HyperInterpStruct *pHIS, ErrorMsg error_message);
  int hyperspherical_forwards_recurrence(int K,
                                         int lmax,
                                         double beta,
                                         double x,
                                         double sinK,
                                         double cotK,
                                         double *sqrtK,
                                         double *one_over_sqrtK,
                                         double *PhiL);
int hyperspherical_forwards_recurrence_chunk(int K,
                                             int lmax,
                                             double beta,
                                             double * __restrict__ x,
                                             double * __restrict__ sinK,
                                             double * __restrict__ cotK,
                                             int chunk,
                                             double * __restrict__ sqrtK,
                                             double * __restrict__ one_over_sqrtK,
                                             double * __restrict__ PhiL);
  int hyperspherical_backwards_recurrence(int K,
                                          int lmax,
                                          double beta,
                                          double x,
                                          double sinK,
                                          double cotK,
                                          double *sqrtK,
                                          double *one_over_sqrtK,
                                          double *PhiL);

  int hyperspherical_backwards_recurrence_chunk(int K,
                                                int lmax,
                                                double beta,
                                                double * __restrict__ x,
                                                double * __restrict__ sinK,
                                                double * __restrict__ cotK,
                                                int chunk,
                                                double * __restrict__ sqrtK,
                                                double * __restrict__ one_over_sqrtK,
                                                double * __restrict__ PhiL);
  int hyperspherical_get_xmin(HyperInterpStruct *pHIS,
                              double xtol,
                              double phiminabs,
                              double *xmin);

  int hyperspherical_WKB(int K,int l,double beta,double y, double *Phi);
  int hyperspherical_WKB_vec(int l,
                             double beta,
                             double *sinK_vec,
                             int size_sinK_vec,
                             double *Phi);
  int ClosedModY(int l, int beta, double *y, int * phisign, int * dphisign);
  int get_CF1(int K,int l,double beta, double cotK, double *CF, int *isign);
  int CF1_from_Gegenbauer(int l, int beta, double sinK, double cotK, double *CF);
  double airy_cheb_approx(double z);
  double coef1(double z);
  double coef2(double z);
  double coef3(double z);
  double coef4(double z);
  double cheb(double x, int n, const double A[]);
  double get_value_at_small_phi(int K,int l,double beta,double Phi);

  double PhiWKB_minus_phiminabs(double x, void *param);

  int hyperspherical_get_xmin_from_Airy(int K,
                                        int l,
                                        double beta,
                                        double xtol,
                                        double phiminabs,
                                        double *xmin,
                                        int *fevals
                                        );

  int fzero_ridder(double (*func)(double, void *),
                   double x1,
                   double x2,
                   double xtol,
                   void *param,
                   double *Fx1,
                   double *Fx2,
                   double *xzero,
                   int *fevals
                   );

  int HypersphericalExplicit(int K,int l, double beta,double x, double *Phi);

  int hyperspherical_get_xmin_from_approx(int K,
                                          int l,
                                          double nu,
                                          double ignore1,
                                          double phiminabs,
                                          double *xmin,
                                          int *ignore2);

  size_t hyperspherical_HIS_size(int nl, int nx);
  int hyperspherical_update_pointers(HyperInterpStruct *pHIS_local,
                                     void * HIS_storage_shared);

  int hyperspherical_Hermite_interpolation_vector(HyperInterpStruct *pHIS,
                                                  int nxi,
                                                  int lnum,
                                                  double *xinterp,
                                                  double *Phi,
                                                  double *dPhi,
                                                  double *d2Phi);

  int hyperspherical_Hermite3_interpolation_vector_Phi(HyperInterpStruct *pHIS,int nxi,int lnum,double *xinterp,double *Phi, ErrorMsg error_message);
  int hyperspherical_Hermite3_interpolation_vector_dPhi(HyperInterpStruct *pHIS,int nxi,int lnum,double *xinterp,double *dPhi, ErrorMsg error_message);
  int hyperspherical_Hermite3_interpolation_vector_d2Phi(HyperInterpStruct *pHIS,int nxi,int lnum,double *xinterp,double *d2Phi, ErrorMsg error_message);
  int hyperspherical_Hermite3_interpolation_vector_PhidPhi(HyperInterpStruct *pHIS,int nxi,int lnum,double *xinterp,double *Phi,double *dPhi, ErrorMsg error_message);
  int hyperspherical_Hermite3_interpolation_vector_Phid2Phi(HyperInterpStruct *pHIS,int nxi,int lnum,double *xinterp,double *Phi,double *d2Phi, ErrorMsg error_message);
  int hyperspherical_Hermite3_interpolation_vector_dPhid2Phi(HyperInterpStruct *pHIS,int nxi,int lnum,double *xinterp,double *dPhi,double *d2Phi, ErrorMsg error_message);
  int hyperspherical_Hermite3_interpolation_vector_PhidPhid2Phi(HyperInterpStruct *pHIS,int nxi,int lnum,double *xinterp,double *Phi,double *dPhi,double *d2Phi, ErrorMsg error_message);
  int hyperspherical_Hermite4_interpolation_vector_Phi(HyperInterpStruct *pHIS,int nxi,int lnum,double *xinterp,double *Phi, ErrorMsg error_message);
  int hyperspherical_Hermite4_interpolation_vector_dPhi(HyperInterpStruct *pHIS,int nxi,int lnum,double *xinterp,double *dPhi, ErrorMsg error_message);
  int hyperspherical_Hermite4_interpolation_vector_d2Phi(HyperInterpStruct *pHIS,int nxi,int lnum,double *xinterp,double *d2Phi, ErrorMsg error_message);
  int hyperspherical_Hermite4_interpolation_vector_PhidPhi(HyperInterpStruct *pHIS,int nxi,int lnum,double *xinterp,double *Phi,double *dPhi, ErrorMsg error_message);
  int hyperspherical_Hermite4_interpolation_vector_Phid2Phi(HyperInterpStruct *pHIS,int nxi,int lnum,double *xinterp,double *Phi,double *d2Phi, ErrorMsg error_message);
  int hyperspherical_Hermite4_interpolation_vector_dPhid2Phi(HyperInterpStruct *pHIS,int nxi,int lnum,double *xinterp,double *dPhi,double *d2Phi, ErrorMsg error_message);
  int hyperspherical_Hermite4_interpolation_vector_PhidPhid2Phi(HyperInterpStruct *pHIS,int nxi,int lnum,double *xinterp,double *Phi,double *dPhi,double *d2Phi, ErrorMsg error_message);
  int hyperspherical_Hermite6_interpolation_vector_Phi(HyperInterpStruct *pHIS,int nxi,int lnum,double *xinterp,double *Phi, ErrorMsg error_message);
  int hyperspherical_Hermite6_interpolation_vector_dPhi(HyperInterpStruct *pHIS,int nxi,int lnum,double *xinterp,double *dPhi, ErrorMsg error_message);
  int hyperspherical_Hermite6_interpolation_vector_d2Phi(HyperInterpStruct *pHIS,int nxi,int lnum,double *xinterp,double *d2Phi, ErrorMsg error_message);
  int hyperspherical_Hermite6_interpolation_vector_PhidPhi(HyperInterpStruct *pHIS,int nxi,int lnum,double *xinterp,double *Phi,double *dPhi, ErrorMsg error_message);
  int hyperspherical_Hermite6_interpolation_vector_Phid2Phi(HyperInterpStruct *pHIS,int nxi,int lnum,double *xinterp,double *Phi,double *d2Phi, ErrorMsg error_message);
  int hyperspherical_Hermite6_interpolation_vector_dPhid2Phi(HyperInterpStruct *pHIS,int nxi,int lnum,double *xinterp,double *dPhi,double *d2Phi, ErrorMsg error_message);
  int hyperspherical_Hermite6_interpolation_vector_PhidPhid2Phi(HyperInterpStruct *pHIS,int nxi,int lnum,double *xinterp,double *Phi,double *dPhi,double *d2Phi, ErrorMsg error_message);


#ifdef __cplusplus
}
#endif

#endif
