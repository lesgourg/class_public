/** @file trg.h Document includes for trg module */

#include "spectra.h"
#include "trg.h"

#ifndef __NONLINEAR__
#define __NONLINEAR__

enum non_linear_method {nl_none,nl_trg_linear,nl_trg_one_loop,nl_trg};

/**
 * Structure containing all information on non-linear spectra.
 *
 * Once initialised by nonlinear_init(), contains a table for all two points correlation functions
 * and for all the ai,bj functions (containing the three points correlation functions), for each
 * time and wave-number.
 */

struct nonlinear {

  /** @name - input parameters initialized by user in input module
      (all other quantitites are computed in this module, given these
      parameters and the content of the 'precision', 'background',
      'thermo', 'primordial' and 'spectra' structures) */

  //@{

  enum non_linear_method method;
  
  //@}

  /** @name - table of k values, and related quantities */

  //@{

  double * k;		/**< table containing the values of k used in this module */
  int k_size; 		/**< total number of k values */

  //@}

  /** @name - tables of time values, and related quantities */

  //@{

  double * eta;		/**< table containing eta values used in this module(new time variable defined as log(a/a_ini)) */
  double * z;		/**< table containing z   values used in this module */
  int eta_size;

  //@}

  /** @name - tables of non-linear spectra and their derivatives*/

  //@{

  double * p_dd;
  double * p_dt;
  double * p_tt;

  /** @name - technical parameters */

  //@{

  short nonlinear_verbose;  	/**< amount of information written in standard output */

  ErrorMsg error_message; 	/**< zone for writing error messages */

  //@}
};

/********************************************************************************/

/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int nonlinear_pk_at_z(
			struct nonlinear * pnl,
			enum linear_or_logarithmic mode, 
			double z,
			double * output_tot
			);

  int nonlinear_pk_at_k_and_z(
			      struct nonlinear * pnl,
			      double k,
			      double z,
			      double * pk
			      );

  int nonlinear_init(
		     struct precision *ppr,
		     struct background *pba,
		     struct thermo *pth,
		     struct primordial *ppm,
		     struct spectra *psp,
		     struct nonlinear *pnl
		     );
  
  int nonlinear_free(
	       struct nonlinear *pnl
	       );
  
#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
