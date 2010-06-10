/** @file bessel.h Documented includes for Bessel module */

#include "perturbations.h"

#ifndef __BESSEL__
#define __BESSEL__

/**
 * All Bessel functions.
 *
 * Once initialized by bessel_init(), contains table of
 * all Bessel functions (for the moment, only spherical Bessels \f$ j_l(x) \f$).
 */
struct bessels {

  /** @name - parameters defining the exact content of the Bessel table */

  //@{

  int l_max; /**< maximum value of l */

  int l_size; /**< number of multipole values */
  int * l; /**< list of multipole values, l[index_l] */

  double x_step; /**< step dx for sampling Bessel functions \f$ j_l(x) \f$ */
  double x_max; /**< last value of x (always a multiple of x_step!) */
  double j_cut; /**< value of \f$ j_l \f$ below which it is approximated by zero (in the region x << l) */

 //@}

  /** @name - Bessel table, and arrays necessary for reading it */

  //@{

  double * x_min; /**< x_min[index_l] is the minimum value of x for l[index_l], given f_cut; always a multiple of x_step */
  int * x_size; /* x_min[index_l] is the number of x values for l[index_l]; hence x_min[index_l]+x_step*(x_size[index_l]-1) = x_max */
  double ** j; /* (j[index_l])[index_x] is \f$ j_l(x) \f$ for l[index_l] and x=x_min[index_l]+x_step*index_x */ 
  double ** ddj; /* (ddj[index_l])[index_x] is the splined \f$ j_l''(x) \f$ for l[index_l] and x=x_min[index_l]+x_step*index_x */ 

  //@}

  /** @name - flag regulating the amount of information sent to standard output (none if set to zero) */

  //@{

  short bessels_verbose;

  //@}

  /** @name - zone for writing error messages */

  //@{

  ErrorMsg error_message;

  //@}

};

/*************************************************************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int bessel_at_x(
		  struct bessels * pbs,
		  double x,
		  int l,
		  double * j
		  );

  int bessel_init(
		  struct precision * ppr,
		  struct background * pba,
		  struct perturbs * ppt,
		  struct bessels * pbs
		  );

  int bessel_free(
		  struct bessels * pbs
		  );

  int bessel_get_l_list(
			struct precision * ppr,
			struct bessels * pbs
			);

  int bessel_j_for_l(
		     struct precision * ppr,
		     struct bessels * pbs,
		     int index_l,
		     double kmin
		     );

  int bessel_j(
	       struct bessels * pbs,
	       int l,
	       double x,
	       double * jl
	       );
    
#ifdef __cplusplus
}
#endif

/**  
 * @name Constants used for computing Bessel function
 */

//@{

#define _GAMMA1_ 2.6789385347077476336556
#define _GAMMA2_ 1.3541179394264004169452

//@}

#endif
