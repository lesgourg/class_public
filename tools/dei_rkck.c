#include "dei_rkck.h"

#define dsign(a,b) ( (b) > 0. ? (a) : (-(a)) )

char Transmit_Error_Message[2048];

double * yscal;
double * y;
double * dydx;

double * yerr;
double * ytempo;

double * ak2;
double * ak3;
double * ak4;
double * ak5;
double * ak6;
double * ytemp;

int n;

/**
 * Initialize the integrator
 *
 */
int initialize_generic_integrator(int n_dim, Generic_integrator_struct *generic_integrator_in){

  /** - Allocate workspace dynamically */ 

  n = n_dim;

  yscal=malloc(sizeof(double)*n_dim);
  y=malloc(sizeof(double)*n_dim);
  dydx=malloc(sizeof(double)*n_dim);

  yerr=malloc(sizeof(double)*n_dim);
  ytempo=malloc(sizeof(double)*n_dim);

  ak2=malloc(sizeof(double)*n_dim);
  ak3=malloc(sizeof(double)*n_dim);
  ak4=malloc(sizeof(double)*n_dim);
  ak5=malloc(sizeof(double)*n_dim);
  ak6=malloc(sizeof(double)*n_dim);
  ytemp=malloc(sizeof(double)*n_dim);

  return _SUCCESS_;
}


/**
 * Free the integrator's memory space
 *
 * Called by background_solve(); thermodynamics_solve_with_recfast(); perturb_solve().
 */
int cleanup_generic_integrator(Generic_integrator_struct *generic_integrator_in){

  free(yscal);
  free(y);
  free(dydx);

  free(yerr);
  free(ytempo);

  free(ak2);
  free(ak3);
  free(ak4);
  free(ak5);
  free(ak6);
  free(ytemp);

  return _SUCCESS_;
}

/**
 * Integration
 *
 * Called by background_solve(); thermodynamics_solve_with_recfast(); perturb_solve().
 */
int generic_integrator ( void (*derivs)(double x, double y[], double yprime[]),
			 double x_start,
			 double x_end,
			 double y_start[],
			 double tol,
			 Generic_integrator_struct *generic_integrator_in
			 ){

  if (odeint(y_start,x_start,x_end,tol,x_end-x_start,0.,derivs)==_FAILURE_) {
    sprintf(generic_integrator_in->error_message,"Problem in odeint \n=>%s",Transmit_Error_Message);
    return _FAILURE_;
  }

  return _SUCCESS_;

}

int odeint(double ystart[],  
	   double x1, 
	   double x2, 
	   double eps, 
	   double h1,
	   double hmin, 
	   void (*derivs)(double, double [], double []))

{
  int nstp,i;
  double x,hnext,hdid,h;

  x=x1;
  h=dsign(h1,x2-x1);
  for (i=0;i<n;i++) y[i]=ystart[i];
  for (nstp=1;nstp<=_MAXSTP_;nstp++) {
    (*derivs)(x,y,dydx);
    for (i=0;i<n;i++)
      yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+_TINY_;
    if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
    rkqs(&x,h,eps,&hdid,&hnext,derivs);
    if ((x-x2)*(x2-x1) >= 0.0) {
      for (i=0;i<n;i++) ystart[i]=y[i];
      return _SUCCESS_;
    }
    if (fabs(hnext) <= hmin) {
      sprintf(Transmit_Error_Message,"%s(L:%d): Step size too small in odeint\n",__func__,__LINE__);
      return _FAILURE_;
    }
    h=hnext;
  }
  sprintf(Transmit_Error_Message,"%s(L:%d): Too many steps in routine odeint\n",__func__,__LINE__);
  return _FAILURE_;
}

int rkqs(double *x, double htry, double eps,
	 double *hdid, double *hnext,
	 void (*derivs)(double, double [], double []))
{

  int i;
  double errmax,h,htemp,xnew;

  h=htry;
  for (;;) {
    if (rkck(*x,h,derivs)==_FAILURE_) {
      sprintf(Transmit_Error_Message,"%s(L:%d): error in rkck\n",__func__,__LINE__);
      return _FAILURE_;
    }
    errmax=0.0;
    for (i=0;i<n;i++) errmax=max(errmax,fabs(yerr[i]/yscal[i]));
    errmax /= eps;
    if (errmax <= 1.0) break;
    htemp=_SAFETY_*h*pow(errmax,_PSHRNK_);
    h=(h >= 0.0 ? max(htemp,0.1*h) : min(htemp,0.1*h));
    xnew=(*x)+h;
    if (xnew == *x) {
      sprintf(Transmit_Error_Message,"%s(L:%d): stepsize underflow in rkqs\n",__func__,__LINE__);
      return _FAILURE_;
    } 
  }
  if (errmax > _ERRCON_) *hnext=_SAFETY_*h*pow(errmax,_PGROW_);
  else *hnext=5.0*h;
  *x += (*hdid=h);
  for (i=0;i<n;i++) y[i]=ytemp[i];

  return _SUCCESS_;
}

int rkck(
	 double x, 
	 double h,
	 void (*derivs)(double, double [], double []))
{
  int i;

  for (i=0;i<n;i++)
    ytemp[i]=y[i]+_RKCK_b21_*h*dydx[i];
  (*derivs)(x+_RKCK_a2_*h,ytemp,ak2);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(_RKCK_b31_*dydx[i]+_RKCK_b32_*ak2[i]);
  (*derivs)(x+_RKCK_a3_*h,ytemp,ak3);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(_RKCK_b41_*dydx[i]+_RKCK_b42_*ak2[i]+_RKCK_b43_*ak3[i]);
  (*derivs)(x+_RKCK_a4_*h,ytemp,ak4);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(_RKCK_b51_*dydx[i]+_RKCK_b52_*ak2[i]+_RKCK_b53_*ak3[i]+_RKCK_b54_*ak4[i]);
  (*derivs)(x+_RKCK_a5_*h,ytemp,ak5);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(_RKCK_b61_*dydx[i]+_RKCK_b62_*ak2[i]+_RKCK_b63_*ak3[i]+_RKCK_b64_*ak4[i]+_RKCK_b65_*ak5[i]);
  (*derivs)(x+_RKCK_a6_*h,ytemp,ak6);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(_RKCK_c1_*dydx[i]+_RKCK_c3_*ak3[i]+_RKCK_c4_*ak4[i]+_RKCK_c6_*ak6[i]);
  for (i=0;i<n;i++)
    yerr[i]=h*(_RKCK_dc1_*dydx[i]+_RKCK_dc3_*ak3[i]+_RKCK_dc4_*ak4[i]+_RKCK_dc5_*ak5[i]+_RKCK_dc6_*ak6[i]);

  return _SUCCESS_;
}


