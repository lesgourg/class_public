/** @file evolver_rkdp45.c Dormand-Prince 5(4) explicit Runge-Kutta evolver
 *
 * An embedded explicit Runge-Kutta pair with adaptive steps, the method
 * of MATLAB's ode45: seven stages, of which the last one is the first
 * stage of the next step (first same as last, FSAL), so that each step
 * costs six evaluations of the derivatives. It has the same arguments
 * as evolver_rk() and evolver_ndf15().
 *
 * Outside of the stiff regimes that tight-coupling and radiation-
 * streaming approximations take care of, the perturbations of the
 * standard cosmological species are not stiff, and there this evolver
 * is faster than ndf15 at the same accuracy. On a stiff system, an
 * explicit method has to take very small steps. Thus, it is used by
 * default for the perturbations only when no species or option that can
 * make them stiff is present (see input_select_perturbations_evolver()).
 *
 * It is the integrator that CLASS++ uses in the same cases, with the
 * same step-size control and the same output:
 *
 * - error control: the largest component of the error estimate relative
 *   to max(|y|, abstol/tolerance) is compared to the tolerance; a step
 *   with a larger error is rejected. The new step is h * max(0.1,
 *   0.8 (tolerance/error)^(1/5)); after two rejections in a row, h is
 *   halved.
 *
 * - output: the values at points x_sampling that fall inside a step are
 *   interpolated. Dormand-Prince's own continuous extension is of order
 *   4 only, one order below the step, so the stored sources would be
 *   systematically less accurate than the solution that the step-size
 *   control protects. Instead, the interpolant is the quintic Hermite
 *   polynomial through the last three step ends, built from the values
 *   and the derivatives there. All of them are already available (the
 *   derivative at a step end is the FSAL stage), so this costs no
 *   evaluation of the derivatives. Only the first step has no history
 *   and uses the continuous extension.
 *
 * - a non-finite value in the solution or in the error estimate stops
 *   the integration with an error message. Otherwise a NaN would be
 *   accepted (NaN > tolerance is false) and an Inf would be rejected
 *   forever, with the step size converging to its lower bound.
 */

#include "evolver_rkdp45.h"

#include <stdint.h>
#include <string.h>

/** Dormand-Prince 5(4) coefficients: nodes c, matrix a (strict lower
    triangle), weights b of the 5th-order solution, error weights e =
    b - bhat, and the continuous extension b_i(theta) = sum_{m=1..4}
    r[i][m-1] theta^m (Shampine's, as in MATLAB's ode45) */

#define _RKDP45_STAGES_ 7

static const double rkdp45_c[_RKDP45_STAGES_] = {0.0, 0.2, 0.3, 0.8, 8.0/9.0, 1.0, 1.0};

static const double rkdp45_a[_RKDP45_STAGES_][_RKDP45_STAGES_] = {
  {0., 0., 0., 0., 0., 0., 0.},
  {0.2, 0., 0., 0., 0., 0., 0.},
  {3.0/40.0, 9.0/40.0, 0., 0., 0., 0., 0.},
  {44.0/45.0, -56.0/15.0, 32.0/9.0, 0., 0., 0., 0.},
  {19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0, 0., 0., 0.},
  {9017.0/3168.0, -355.0/33.0, 46732.0/5247.0, 49.0/176.0, -5103.0/18656.0, 0., 0.},
  {35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0, 0.}};

static const double rkdp45_b[_RKDP45_STAGES_] = {35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0, 0.0};

static const double rkdp45_e[_RKDP45_STAGES_] = {71.0/57600.0, 0.0, -71.0/16695.0, 71.0/1920.0, -17253.0/339200.0, 22.0/525.0, -1.0/40.0};

static const double rkdp45_r[_RKDP45_STAGES_][4] = {
  {1.0, -183.0/64.0, 37.0/12.0, -145.0/128.0},
  {0.0, 0.0, 0.0, 0.0},
  {0.0, 1500.0/371.0, -1000.0/159.0, 1000.0/371.0},
  {0.0, -125.0/32.0, 125.0/12.0, -375.0/64.0},
  {0.0, 9477.0/3392.0, -729.0/106.0, 25515.0/6784.0},
  {0.0, -11.0/7.0, 11.0/3.0, -55.0/28.0},
  {0.0, 1.5, -4.0, 2.5}};

/** step-size control */
#define _RKDP45_SAFETY_ 0.8   /**< safety factor on the step proposed by the error estimate */
#define _RKDP45_FAC_MIN_ 0.1  /**< smallest factor by which a step can shrink */
#define _RKDP45_ABSTOL_ 1.e-15 /**< absolute tolerance, as in evolver_ndf15() */

/**
 * Whether x is NaN or Inf. A bit test rather than isnan() or isinf(),
 * which -ffast-math allows the compiler to fold to false.
 */

static int rkdp45_is_non_finite(double x) {
  uint64_t bits;
  memcpy(&bits,&x,sizeof(bits));
  return ((bits & 0x7ff0000000000000ULL) == 0x7ff0000000000000ULL);
}

/**
 * The integration itself, in a workspace allocated (and freed) by
 * evolver_rkdp45(): 14 vectors of size y_size.
 */

static int evolver_rkdp45_steps(int (*derivs)(double x,
                                              double * y,
                                              double * dy,
                                              void * parameters_and_workspace,
                                              ErrorMsg error_message),
                                double x_ini,
                                double x_end,
                                double * y,
                                int * used_in_output,
                                int neq,
                                void * parameters_and_workspace_for_derivs,
                                double rtol,
                                double * x_sampling,
                                int x_size,
                                int (*output)(double x,
                                              double y[],
                                              double dy[],
                                              int index_x,
                                              void * parameters_and_workspace,
                                              ErrorMsg error_message),
                                int (*print_variables)(double x,
                                                       double y[],
                                                       double dy[],
                                                       void * parameters_and_workspace,
                                                       ErrorMsg error_message),
                                double * workspace,
                                ErrorMsg error_message) {

  const int s = _RKDP45_STAGES_;
  const double threshold = _RKDP45_ABSTOL_/rtol;
  const double pow_grow = 0.2;

  /* stages k_i, i=0..6, then the vectors of size neq */
  double * ki = workspace;
  double * ynew = ki + s*neq;
  double * ytemp = ynew + neq;
  double * err = ytemp + neq;
  double * yinterp = err + neq;
  double * dyinterp = yinterp + neq;
  /* value and derivative at the end of the step before the current one,
     for the Hermite interpolant */
  double * y_prev = dyinterp + neq;
  double * f_prev = y_prev + neq;
  /* derivative at the end of an accepted step: the last (FSAL) stage */
  double * f_end = ki + (s-1)*neq;

  double t,tnew,t_prev=0.,ti,th;
  double h,hnew,hmin,hmax,absh,rh;
  double errmax,errtemp;
  double thpow[5];
  double bi_vec_y[_RKDP45_STAGES_];
  double bi_vec_dy[_RKDP45_STAGES_];
  int tdir,idx,i,j,k,m,kbad;
  short nofailed,nonfinite,have_prev=_FALSE_;

  class_test(x_sampling == NULL,
             error_message,
             "evolver_rkdp45 requires a non-null x_sampling array");
  class_test(rtol <= 0.,
             error_message,
             "evolver_rkdp45 requires tolerance > 0 (got %e)",rtol);

  t = x_ini;

  /** - first stage of the first step: k_0 = f(t, y) */
  class_call((*derivs)(t,y,ki,parameters_and_workspace_for_derivs,error_message),
             error_message,
             error_message);

  /** - initial step, from the scale of the derivatives */
  hmax = fabs(x_end-x_ini)/10.;
  if (x_size > 1)
    absh = MIN(hmax,fabs(x_sampling[1]-x_sampling[0]));
  else
    absh = hmax;
  if (absh == 0.)
    absh = hmax;

  rh = 0.;
  for (k=0; k<neq; k++)
    rh = MAX(rh,fabs(ki[k])/MAX(fabs(y[k]),threshold));
  rh /= _RKDP45_SAFETY_*pow(rtol,pow_grow);
  if (absh*rh > 1.)
    absh = 1./rh;

  tdir = (x_end > x_ini) ? 1 : -1;
  hnew = absh*tdir;
  nofailed = _TRUE_;

  /** - index of the first output point */
  idx = 0;
  if ((t-x_end)*tdir < 0.) {
    for (idx=0; (idx<x_size) && ((x_sampling[idx]-t)*tdir < 0.); idx++);
  }

  while ((t-x_end)*tdir < 0.) {

    h = hnew;
    hmin = 100.*DBL_MIN*fabs(t);
    class_test(fabs(h) < hmin,
               error_message,
               "evolver_rkdp45: step size %e fell below minimum %e at x=%e",fabs(h),hmin,t);
    if (fabs(h) > 0.9*fabs(x_end-t))
      h = x_end-t;

    /** - try a step h: stages 1 to 6 (stage 0 is the last stage of the previous step) */
    for (k=0; k<neq; k++) {
      ynew[k] = y[k]+h*rkdp45_b[0]*ki[k];
      err[k] = h*rkdp45_e[0]*ki[k];
    }
    for (i=1; i<s; i++) {
      for (k=0; k<neq; k++)
        ytemp[k] = y[k];
      for (j=0; j<i; j++) {
        for (k=0; k<neq; k++)
          ytemp[k] += h*rkdp45_a[i][j]*ki[j*neq+k];
      }
      class_call((*derivs)(t+rkdp45_c[i]*h,ytemp,ki+i*neq,parameters_and_workspace_for_derivs,error_message),
                 error_message,
                 error_message);
      for (k=0; k<neq; k++) {
        ynew[k] += h*rkdp45_b[i]*ki[i*neq+k];
        err[k] += h*rkdp45_e[i]*ki[i*neq+k];
      }
    }

    /** - error estimate. Non-finite values are tested component by
        component: the maximum cannot carry a NaN out, since NaN > x is
        false */
    errmax = 0.;
    nonfinite = _FALSE_;
    for (k=0; k<neq; k++) {
      errtemp = fabs(err[k]/MAX(threshold,fabs(ynew[k])));
      if (rkdp45_is_non_finite(ynew[k]) || rkdp45_is_non_finite(err[k]))
        nonfinite = _TRUE_;
      errmax = MAX(errmax,errtemp);
    }

    if (nonfinite == _TRUE_) {
      for (kbad=0; kbad<neq; kbad++) {
        if (rkdp45_is_non_finite(ynew[kbad]) || rkdp45_is_non_finite(err[kbad]))
          break;
      }
      class_stop(error_message,
                 "evolver_rkdp45: non-finite state at x=%e (h=%e): y[%d]=%e, err[%d]=%e. Something upstream produced NaN or Inf.",
                 t,h,kbad,ynew[kbad],kbad,err[kbad]);
    }

    /** - step rejected: shrink by the factor that the error estimate
        suggests, or halve the step after two rejections in a row */
    if (errmax > rtol) {
      if (nofailed == _TRUE_) {
        nofailed = _FALSE_;
        hnew = tdir*MAX(hmin,fabs(h)*MAX(_RKDP45_FAC_MIN_,_RKDP45_SAFETY_*pow(rtol/errmax,pow_grow)));
      }
      else {
        hnew = tdir*MAX(hmin,0.5*fabs(h));
      }
      continue;
    }

    /** - step accepted */
    if (print_variables != NULL) {
      class_call((*print_variables)(t+h,ynew,f_end,parameters_and_workspace_for_derivs,error_message),
                 error_message,
                 error_message);
    }

    hnew = tdir*MAX(hmin,fabs(h)*MAX(_RKDP45_FAC_MIN_,_RKDP45_SAFETY_*pow(rtol/errmax,pow_grow)));
    nofailed = _TRUE_;
    tnew = t+h;

    /** - output at all sampling points within (t, tnew] */
    for (; (idx<x_size) && ((tnew-x_sampling[idx])*tdir >= 0.); idx++) {

      if (tnew == x_sampling[idx]) {
        class_call((*output)(tnew,ynew,f_end,idx,parameters_and_workspace_for_derivs,error_message),
                   error_message,
                   error_message);
        continue;
      }

      ti = x_sampling[idx];

      /** - --> continuous extension of the step */
      th = (ti-t)/h;
      thpow[0] = 1.;
      thpow[1] = th;
      thpow[2] = th*th;
      thpow[3] = thpow[2]*th;
      thpow[4] = thpow[2]*thpow[2];
      for (i=0; i<s; i++) {
        bi_vec_y[i] = 0.;
        bi_vec_dy[i] = 0.;
        for (m=0; m<4; m++) {
          bi_vec_y[i] += rkdp45_r[i][m]*thpow[m+1];
          bi_vec_dy[i] += rkdp45_r[i][m]*(m+1)*thpow[m];
        }
      }
      for (k=0; k<neq; k++) {
        if (used_in_output[k] == _TRUE_) {
          yinterp[k] = y[k];
          dyinterp[k] = 0.;
          for (i=0; i<s; i++) {
            yinterp[k] += h*bi_vec_y[i]*ki[i*neq+k];
            dyinterp[k] += bi_vec_dy[i]*ki[i*neq+k];
          }
        }
      }

      /** - --> quintic Hermite interpolant through (t_prev, t, tnew),
          with value and derivative at each node, when there is a
          previous step:
          h_i(x) = [1 - 2 (x-x_i) L_i'(x_i)] L_i(x)^2 and
          hbar_i(x) = (x-x_i) L_i(x)^2, with L_i the Lagrange basis of
          the three nodes */
      if (have_prev == _TRUE_) {
        double x0 = t_prev, x1 = t, x2 = tnew;
        double d01 = x0-x1, d02 = x0-x2, d12 = x1-x2;
        if ((d01 != 0.) && (d02 != 0.) && (d12 != 0.)) {
          double u0 = ti-x0, u1 = ti-x1, u2 = ti-x2;
          double L0 = (u1*u2)/(d01*d02);
          double L1 = (u0*u2)/(-d01*d12);
          double L2 = (u0*u1)/(d02*d12);
          /* L_i'(x_i) = sum_{j != i} 1/(x_i-x_j) */
          double Lp0 = 1./d01+1./d02;
          double Lp1 = -1./d01+1./d12;
          double Lp2 = -1./d02-1./d12;
          /* derivatives of L_i at ti */
          double dL0 = (u1+u2)/(d01*d02);
          double dL1 = (u0+u2)/(-d01*d12);
          double dL2 = (u0+u1)/(d02*d12);
          double L0s = L0*L0, L1s = L1*L1, L2s = L2*L2;
          double w0 = (1.-2.*u0*Lp0)*L0s;
          double w1 = (1.-2.*u1*Lp1)*L1s;
          double w2 = (1.-2.*u2*Lp2)*L2s;
          double v0 = u0*L0s, v1 = u1*L1s, v2 = u2*L2s;
          double dw0 = -2.*Lp0*L0s+(1.-2.*u0*Lp0)*2.*L0*dL0;
          double dw1 = -2.*Lp1*L1s+(1.-2.*u1*Lp1)*2.*L1*dL1;
          double dw2 = -2.*Lp2*L2s+(1.-2.*u2*Lp2)*2.*L2*dL2;
          double dv0 = L0s+u0*2.*L0*dL0;
          double dv1 = L1s+u1*2.*L1*dL1;
          double dv2 = L2s+u2*2.*L2*dL2;
          for (k=0; k<neq; k++) {
            if (used_in_output[k] == _TRUE_) {
              yinterp[k] = w0*y_prev[k]+w1*y[k]+w2*ynew[k]+v0*f_prev[k]+v1*ki[k]+v2*f_end[k];
              dyinterp[k] = dw0*y_prev[k]+dw1*y[k]+dw2*ynew[k]+dv0*f_prev[k]+dv1*ki[k]+dv2*f_end[k];
            }
          }
        }
      }

      class_call((*output)(ti,yinterp,dyinterp,idx,parameters_and_workspace_for_derivs,error_message),
                 error_message,
                 error_message);
    }

    /** - the start of this step is the history of the next one */
    for (k=0; k<neq; k++) {
      if (used_in_output[k] == _TRUE_) {
        y_prev[k] = y[k];
        f_prev[k] = ki[k];
      }
    }
    t_prev = t;
    have_prev = _TRUE_;

    /** - move on: the last stage is the first one of the next step */
    for (k=0; k<neq; k++) {
      y[k] = ynew[k];
      ki[k] = f_end[k];
    }
    t = tnew;
  }

  return _SUCCESS_;
}

/**
 * Integrate the system dy/dx = derivs(x, y) from x_ini to x_end with
 * the Dormand-Prince 5(4) method, calling output() at each point of
 * x_sampling. The arguments are those of evolver_rk() and
 * evolver_ndf15(); minimum_variation, evaluate_timescale and
 * timestep_over_timescale are not used.
 */

int evolver_rkdp45(int (*derivs)(double x,
                                 double * y,
                                 double * dy,
                                 void * parameters_and_workspace,
                                 ErrorMsg error_message),
                   double x_ini,
                   double x_end,
                   double * y,
                   int * used_in_output,
                   int y_size,
                   void * parameters_and_workspace_for_derivs,
                   double tolerance,
                   double minimum_variation,
                   int (*evaluate_timescale)(double x,
                                             void * parameters_and_workspace,
                                             double * timescale,
                                             ErrorMsg error_message),
                   double timestep_over_timescale,
                   double * x_sampling,
                   int x_size,
                   int (*output)(double x,
                                 double y[],
                                 double dy[],
                                 int index_x,
                                 void * parameters_and_workspace,
                                 ErrorMsg error_message),
                   int (*print_variables)(double x,
                                          double y[],
                                          double dy[],
                                          void * parameters_and_workspace,
                                          ErrorMsg error_message),
                   ErrorMsg error_message) {

  double * workspace;
  int status;

  class_alloc(workspace,(_RKDP45_STAGES_+7)*y_size*sizeof(double),error_message);

  status = evolver_rkdp45_steps(derivs,
                                x_ini,
                                x_end,
                                y,
                                used_in_output,
                                y_size,
                                parameters_and_workspace_for_derivs,
                                tolerance,
                                x_sampling,
                                x_size,
                                output,
                                print_variables,
                                workspace,
                                error_message);

  free(workspace);

  return status;
}
