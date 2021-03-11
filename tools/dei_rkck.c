#include "dei_rkck.h"

/**
 * Initialize the integrator
 *
 */
int initialize_generic_integrator(
				  int n_dim,
				  struct generic_integrator_workspace * pgi){

  /** - Allocate workspace dynamically */

  pgi->n = n_dim;

  class_alloc(pgi->yscal,
	      sizeof(double)*n_dim,
	      pgi->error_message);
  class_alloc(pgi->y,
	      sizeof(double)*n_dim,
	      pgi->error_message);
  class_alloc(pgi->dydx,
	      sizeof(double)*n_dim,
	      pgi->error_message);

  class_alloc(pgi->yerr,
	      sizeof(double)*n_dim,
	      pgi->error_message);
  class_alloc(pgi->ytempo,
	      sizeof(double)*n_dim,
	      pgi->error_message);

  class_alloc(pgi->ak2,
	      sizeof(double)*n_dim,
	      pgi->error_message);
  class_alloc(pgi->ak3,
	      sizeof(double)*n_dim,
	      pgi->error_message);
  class_alloc(pgi->ak4,
	      sizeof(double)*n_dim,
	      pgi->error_message);
  class_alloc(pgi->ak5,
	      sizeof(double)*n_dim,
	      pgi->error_message);
  class_alloc(pgi->ak6,
	      sizeof(double)*n_dim,
	      pgi->error_message);
  class_alloc(pgi->ytemp,
	      sizeof(double)*n_dim,
	      pgi->error_message);

  return _SUCCESS_;
}


/**
 * Free the integrator's memory space
 *
 * Called by background_solve(); thermodynamics_solve_with_recfast(); perturbations_solve().
 */
int cleanup_generic_integrator(struct generic_integrator_workspace * pgi){

  free(pgi->yscal);
  free(pgi->y);
  free(pgi->dydx);

  free(pgi->yerr);
  free(pgi->ytempo);

  free(pgi->ak2);
  free(pgi->ak3);
  free(pgi->ak4);
  free(pgi->ak5);
  free(pgi->ak6);
  free(pgi->ytemp);

  return _SUCCESS_;
}

int generic_integrator(int (*derivs)(double x, double y[], double yprime[], void * parameters_and_workspace, ErrorMsg error_message),
		       double x1,
		       double x2,
		       double ystart[],
		       void * parameters_and_workspace_for_derivs,
		       double eps,
		       double hmin,
		       struct generic_integrator_workspace * pgi)

{
  int nstp,i;
  double x,hnext,hdid,h,h1;

  h1=x2-x1;
  x=x1;
  h=dsign(h1,x2-x1);
  for (i=0;i<pgi->n;i++) pgi->y[i]=ystart[i];
  for (nstp=1;nstp<=_MAXSTP_;nstp++) {
    class_call((*derivs)(x,pgi->y,pgi->dydx,parameters_and_workspace_for_derivs, pgi->error_message),
	       pgi->error_message,
	       pgi->error_message);
    for (i=0;i<pgi->n;i++)
      pgi->yscal[i]=fabs(pgi->y[i])+fabs(pgi->dydx[i]*h)+_TINY_;
    if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
    class_call(rkqs(&x,
		    h,
		    eps,
		    &hdid,
		    &hnext,
		    derivs,
		    parameters_and_workspace_for_derivs,
		    pgi),
	       pgi->error_message,
	       pgi->error_message);
    if ((x-x2)*(x2-x1) >= 0.0) {
      for (i=0;i<pgi->n;i++) ystart[i]=pgi->y[i];
      return _SUCCESS_;
    }
    class_test(fabs(hnext/x1) <= hmin,
	       pgi->error_message,
	       "Step size too small: step:%g, minimum:%g, in interval: [%g:%g]",
	       fabs(hnext/x1),
	       hmin,
	       x1,
	       x2);
    h=hnext;
  }

  class_stop(pgi->error_message,
	     "Too many integration steps needed within interval [%g : %g],\n the system of equations is probably buggy or featuring a discontinuity",x1,x2);

}

int rkqs(double *x, double htry, double eps,
	 double *hdid, double *hnext,
	 int (*derivs)(double, double [], double [], void * parameters_and_workspace, ErrorMsg error_message),
	 void * parameters_and_workspace_for_derivs,
	 struct generic_integrator_workspace * pgi)
{

  int i;
  double errmax,h,htemp,xnew;

  h=htry;
  for (;;) {
    class_call(rkck(*x,h,derivs,parameters_and_workspace_for_derivs,pgi),
	       pgi->error_message,
	       pgi->error_message);
    errmax=0.0;
    for (i=0;i<pgi->n;i++) errmax=MAX(errmax,fabs(pgi->yerr[i]/pgi->yscal[i]));
    errmax /= eps;
    if (errmax <= 1.0) break;
    htemp=_SAFETY_*h*pow(errmax,_PSHRNK_);
    h=(h >= 0.0 ? MAX(htemp,0.1*h) : MIN(htemp,0.1*h));
    xnew=(*x)+h;
    class_test(xnew == *x,
	       pgi->error_message,
	       "stepsize underflow at x=%e",xnew);
  }
  if (errmax > _ERRCON_) *hnext=_SAFETY_*h*pow(errmax,_PGROW_);
  else *hnext=5.0*h;
  *x += (*hdid=h);
  for (i=0;i<pgi->n;i++) pgi->y[i]=pgi->ytemp[i];

  return _SUCCESS_;
}

int rkck(
	 double x,
	 double h,
	 int (*derivs)(double, double [], double [], void * parameters_and_workspace, ErrorMsg error_message),
	 void * parameters_and_workspace_for_derivs,
	 struct generic_integrator_workspace * pgi)
{
  int i;

  for (i=0;i<pgi->n;i++)
    pgi->ytemp[i]=pgi->y[i]+_RKCK_b21_*h*pgi->dydx[i];

  class_call((*derivs)(x+_RKCK_a2_*h,
		       pgi->ytemp,
		       pgi->ak2,
		       parameters_and_workspace_for_derivs,
		       pgi->error_message),
	     pgi->error_message,
	     pgi->error_message);

  for (i=0;i<pgi->n;i++)
    pgi->ytemp[i]=pgi->y[i]+h*(_RKCK_b31_*pgi->dydx[i]+_RKCK_b32_*pgi->ak2[i]);

  class_call((*derivs)(x+_RKCK_a3_*h,
		       pgi->ytemp,
		       pgi->ak3,
		       parameters_and_workspace_for_derivs,
		       pgi->error_message),
	     pgi->error_message,
	     pgi->error_message);

  for (i=0;i<pgi->n;i++)
    pgi->ytemp[i]=pgi->y[i]+h*(_RKCK_b41_*pgi->dydx[i]+_RKCK_b42_*pgi->ak2[i]+_RKCK_b43_*pgi->ak3[i]);

  class_call((*derivs)(x+_RKCK_a4_*h,
		       pgi->ytemp,
		       pgi->ak4,
		       parameters_and_workspace_for_derivs,
		       pgi->error_message),
	     pgi->error_message,
	     pgi->error_message);

  for (i=0;i<pgi->n;i++)
    pgi->ytemp[i]=pgi->y[i]+h*(_RKCK_b51_*pgi->dydx[i]+_RKCK_b52_*pgi->ak2[i]+_RKCK_b53_*pgi->ak3[i]+_RKCK_b54_*pgi->ak4[i]);

  class_call((*derivs)(x+_RKCK_a5_*h,
		       pgi->ytemp,
		       pgi->ak5,parameters_and_workspace_for_derivs,
		       pgi->error_message),
	     pgi->error_message,
	     pgi->error_message);

  for (i=0;i<pgi->n;i++)
    pgi->ytemp[i]=pgi->y[i]+h*(_RKCK_b61_*pgi->dydx[i]+_RKCK_b62_*pgi->ak2[i]+_RKCK_b63_*pgi->ak3[i]+_RKCK_b64_*pgi->ak4[i]+_RKCK_b65_*pgi->ak5[i]);

  class_call((*derivs)(x+_RKCK_a6_*h,
		       pgi->ytemp,
		       pgi->ak6,
		       parameters_and_workspace_for_derivs,
		       pgi->error_message),
	     pgi->error_message,
	     pgi->error_message);

  for (i=0;i<pgi->n;i++)
    pgi->ytemp[i]=pgi->y[i]+h*(_RKCK_c1_*pgi->dydx[i]+_RKCK_c3_*pgi->ak3[i]+_RKCK_c4_*pgi->ak4[i]+_RKCK_c6_*pgi->ak6[i]);

  for (i=0;i<pgi->n;i++)
    pgi->yerr[i]=h*(_RKCK_dc1_*pgi->dydx[i]+_RKCK_dc3_*pgi->ak3[i]+_RKCK_dc4_*pgi->ak4[i]+_RKCK_dc5_*pgi->ak5[i]+_RKCK_dc6_*pgi->ak6[i]);

  return _SUCCESS_;
}


