/** @file bessel.h Documented includes for Bessel module */

#ifndef __BESSEL__
#define __BESSEL__

#include "common.h"
#include "arrays.h"

/**
 * Structure containing everything about spherical bessel functions
 * that other modules need to know.
 *
 * Once initialized by bessel_init(), contains table of
 * spherical Bessel functions \f$ j_l(x) \f$).
 */

struct bessels {

  /** @name - input parameters initialized by user in input module
      (all other quantitites are computed in this module, given these
      parameters and the content of the 'precision' structure) */

  //@{

  int l_max; /**< last l value */

  double x_max; /**< maximum value of x (always multiple of x-step) */

  double x_step; /**< step dx for sampling Bessel functions */

  short bessel_always_recompute; /**< if set to true, Bessels are never read from / written in files */

  short use_pbs; /**< if _TRUE_, try to get bessel function from this module, if _FALSE_, from new hypershperical module. For non-flat models this parameter is forced by input module to be _FALSE_ */

  int get_HIS_from_shared_memory; /**< flag specifying if class should try to get HIS from shared memory */

 //@}

  /** @name - parameters defining uniquely the exact content of the Bessel table
      (hence when reading a file, will compare these values with needed value
      in order to take the decision to recompute or not) */

  //@{

  int l_size; /**< number of multipole values */
  int * l; /**< list of multipole values, l[index_l] */

  double j_cut; /**< value of \f$ j_l \f$ below which it is approximated by zero (in the region x << l) */
  int has_dj;   /**< set to true means j_l'(x) also need to be stored */

 //@}

  /** @name - Bessel table, and arrays necessary for reading it */

  //@{

  int * x_size; /**< x_size[index_l] is the number of x values for l[index_l]; hence *x_min[index_l]+x_step*(x_size[index_l]-1) = x_max */

  int x_size_max;

  double ** buffer; /**< buffer[index_l] is a pointer towards a memory zone containing x_min, all j_l(x) and all j_l''(x) for each value of l */

  double ** x_min; /**< x_min[index_l] is a pointer towards the minimum value of x for l[index_l], given j_cut; always a multiple of x_step */

  double ** j; /* (j[index_l])[index_x] is \f$ j_l(x) \f$ for l[index_l] and x=x_min[index_l]+x_step*index_x */

  double ** ddj; /* (ddj[index_l])[index_x] \f$ j_l''(x) \f$ for l[index_l] and x=x_min[index_l]+x_step*index_x (in view of spline interpolation) */

  double ** dj; /* (dj[index_l])[index_x] is \f$ j_l'(x) \f$ for l[index_l] and x=x_min[index_l]+x_step*index_x */

  double ** dddj; /* (dddj[index_l])[index_x] \f$ j_l'''(x) \f$ for l[index_l] and x=x_min[index_l]+x_step*index_x (in view of spline interpolation) */

  //@}

  /** @name - technical parameters */

  //@{

  short bessels_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */

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
		     int index_l
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
