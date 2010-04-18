#include "precision.h"

#ifndef __DEI__
#define __DEI__

struct integrator_structure_in     { 

  /** 
    * zone for writing error messages 
    */
  ErrorMsg error_message;
                                     };

typedef struct integrator_structure_in Generic_integrator_struct;

/**************************************************************/

/**
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int initialize_generic_integrator(int n_dim, Generic_integrator_struct *generic_integrator_in);

  int cleanup_generic_integrator(Generic_integrator_struct * generic_integrator_in);

  int generic_integrator ( void (*derivs)(double x, double y[], double yprime[]),
			 double x_start,
			 double x_end,
			 double y_start[],
			 double tol,
			 Generic_integrator_struct *generic_integrator_in
			 );

  int odeint(double ystart[], 
	     double x1, 
	     double x2, 
	     double eps, 
	     double h1,
	     double hmin, 
	     void (*derivs)(double, double [], double []));

  int rkqs(double *x, 
	   double htry, 
	   double eps,
	   double *hdid, 
	   double *hnext,
	   void (*derivs)(double, double [], double []));

  int rkck(double x, 
	   double h,
	   void (*derivs)(double, double [], double []));

#ifdef __cplusplus
}
#endif

/**************************************************************/

#define _MAXSTP_ 10000
#define _TINY_ 1.0e-30
#define _SAFETY_ 0.9
#define _PGROW_ -0.2
#define _PSHRNK_ -0.25
#define _ERRCON_ 1.89e-4

#define _RKCK_a2_ 0.2
#define _RKCK_a3_ 0.3
#define _RKCK_a4_ 0.6
#define _RKCK_a5_ 1.0
#define _RKCK_a6_ 0.875
#define _RKCK_b21_ 0.2
#define _RKCK_b31_ 3.0/40.0
#define _RKCK_b32_ 9.0/40.0
#define _RKCK_b41_ 0.3
#define _RKCK_b42_ -0.9
#define _RKCK_b43_ 1.2
#define _RKCK_b51_ -11.0/54.0
#define _RKCK_b52_ 2.5
#define _RKCK_b53_ -70.0/27.0
#define _RKCK_b54_ 35.0/27.0
#define _RKCK_b61_ 1631.0/55296.0
#define _RKCK_b62_ 175.0/512.0
#define _RKCK_b63_ 575.0/13824.0
#define _RKCK_b64_ 44275.0/110592.0
#define _RKCK_b65_ 253.0/4096.0
#define _RKCK_c1_ 37.0/378.0
#define _RKCK_c3_ 250.0/621.0
#define _RKCK_c4_ 125.0/594.0
#define _RKCK_c6_ 512.0/1771.0
#define _RKCK_dc5_ -277.00/14336.0
#define _RKCK_dc1_ (37.0/378.0-2825.0/27648.)
#define _RKCK_dc3_ (250.0/621.0-18575.0/48384.0)
#define _RKCK_dc4_ (125.0/594.0-13525.0/55296.0)
#define _RKCK_dc6_ (512.0/1771.0-0.25)

#endif
