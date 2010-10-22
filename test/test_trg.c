/** @file test_trg.c 
 * 
 */
 
#include "precision.h"
#include "background.h"
#include "thermodynamics.h"
#include "perturbations.h"
#include "bessel.h"
#include "transfer.h"
#include "primordial.h"
#include "spectra.h"
#include "output.h"
#include "trg.h"

main(int argc, char **argv) {

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions */
  struct bessels bs;          /* for bessel functions */
  struct transfers tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct output op;
  struct spectra sp;          /* for output spectra */
  struct spectra_nl nl; 

  ErrorMsg errmsg;


  /******* just for testing routine array_interpolate_extrapolate_logspline_one_column *****/

/*   int ix,nx; */
/*   int ny; */
/*   int nxstop; */
/*   int i; */
/*   double*x; */
/*   double*y; */
/*   double*ddy; */

/*   int nxx; */
/*   double*xx; */
/*   double*yy; */

/*   nx=100; */
/*   ny=4; */
/*   x=malloc(sizeof(double)*nx); */
/*   y=malloc(sizeof(double)*nx*ny); */
/*   ddy=malloc(sizeof(double)*nx*ny); */
/*   for (ix=0; ix<nx; ix++) { */
/*     x[ix]=pow(1.2,ix); */
/*     y[ix]=exp(sin(log(x[ix]))); */
/*     y[nx+ix]=exp(cos(log(x[ix]))); */
/*     y[2*nx+ix]=exp(2.*sin(log(x[ix]))); */
/*     y[3*nx+ix]=exp(2.*cos(log(x[ix]))); */
/*   } */


/*   for (ix=1; ix<nx-1; ix++) { */
/*     printf("%g %g %g %g %g\n",x[ix],y[2*nx+ix],0.,-(log(y[2*nx+ix+1])-log(y[2*nx+ix-1]))/(log(x[ix+1])-log(x[ix-1])),0.); */
/*   } */
/*   printf("\n\n"); */

/*   nxstop=80; */

/*   array_logspline_table_one_column(x, */
/* 				   nx, */
/* 				   nxstop, */
/* 				   y, */
/* 				   ny, */
/* 				   2, */
/* 				   ddy, */
/* 				   //_SPLINE_EST_DERIV_, */
/* 				   _SPLINE_NATURAL_, */
/* 				   errmsg); */

/*   nxx=1000; */
/*   xx=malloc(sizeof(double)*nxx); */
/*   yy=malloc(sizeof(double)*nxx); */

/*   for (i=0; i <nxx; i++) { */
/*     xx[i]=pow(1.02,i); */
/*     array_interpolate_extrapolate_logspline_one_column(x, */
/* 						       nx, */
/* 						       nxstop, */
/* 						       y, */
/* 						       ny, */
/* 						       2, */
/* 						       ddy, */
/* 						       xx[i], */
/* 						       &(yy[i]), */
/* 						       errmsg); */
/*   } */

/*   for (i=1; i <nxx-1; i++) { */
/*     printf("%g %g %g %g %g\n",xx[i],0.,yy[i],0.,-(log(yy[i+1])-log(yy[i-1]))/(log(xx[i+1])-log(xx[i-1]))); */
/*   } */

/*   return _SUCCESS_; */

  if (input_init_from_arguments(argc, argv,&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }
 
  pt.k_scalar_kmax_for_pk=2000.; 
  /*  pt.k_scalar_kmax_for_pk=600.; */
  pr.k_scalar_k_per_decade_for_pk=10.; 

  pt.has_cl_cmb_temperature = _FALSE_;
  pt.has_cl_cmb_polarization = _FALSE_;
  pt.has_cl_cmb_lensing_potential = _FALSE_;
  pt.has_pk_matter = _TRUE_;

  sp.z_max_pk = 50.; 

  nl.k_max=pt.k_scalar_kmax_for_pk*ba.h-1.;

  if (background_init(&pr,&ba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_init(&pr,&ba,&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (perturb_init(&pr,&ba,&th,&pt) == _FAILURE_) {
    printf("\n\nError in perturb_init \n=>%s\n",pt.error_message);
    return _FAILURE_;
  }

  if (bessel_init(&pr,&bs) == _FAILURE_) {
    printf("\n\nError in bessel_init \n =>%s\n",bs.error_message);
    return _FAILURE_;
  }

  if (transfer_init(&pr,&ba,&th,&pt,&bs,&tr) == _FAILURE_) {
    printf("\n\nError in transfer_init \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  if (primordial_init(&pr,&pt,&pm) == _FAILURE_) {
    printf("\n\nError in transfer_init \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }

  if (spectra_init(&ba,&pt,&tr,&pm,&sp) == _FAILURE_) {
    printf("\n\nError in spectra_init \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }

  if (trg_init(&pr,&ba,&pm,&sp,&nl) == _FAILURE_) {
    printf("\n\nError in trg_init \n=>%s\n",nl.error_message);
    return _FAILURE_;
  }

  /****** done ******/

  
  if (trg_free(&nl) == _FAILURE_) {
    printf("\n\nError in trg_free \n=>%s\n",nl.error_message);
    return _FAILURE_;
  }

  if (spectra_free(&sp) == _FAILURE_) {
    printf("\n\nError in spectra_free \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }

  if (primordial_free(&pm) == _FAILURE_) {
    printf("\n\nError in primordial_free \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }

  if (transfer_free(&tr) == _FAILURE_) {
    printf("\n\nError in transfer_free \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  if (bessel_free(&bs) == _FAILURE_)  {
    printf("\n\nError in bessel_free \n=>%s\n",bs.error_message);
    return _FAILURE_;
  }

  if (perturb_free(&pt) == _FAILURE_) {
    printf("\n\nError in perturb_free \n=>%s\n",pt.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_free(&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_free \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (background_free(&ba) == _FAILURE_) {
    printf("\n\nError in background_free \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }


}
