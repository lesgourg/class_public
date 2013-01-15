/** @file trg.h Document includes for trg module */

#include "spectra.h"
#include "trg.h"

#ifndef __NONLINEAR__
#define __NONLINEAR__

#define _M_EV_TOO_BIG_FOR_HALOFIT_ 10. /**< above which value of non-CDM mass (in eV) do we stop trusting halofit? */

enum non_linear_method {nl_none,nl_halofit,nl_trg_linear,nl_trg_one_loop,nl_trg};

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
  enum non_linear_ic ic;

  //@}

  /** @name - table of k values, and related quantities */

  //@{

  double * k;		/**< table containing the values of k used in this module */
  int * k_size; 	/**< k_size[index_z] total number of k values at a given redshift z[index_z] */

  //@}

  /** @name - tables of time values, and related quantities */

  //@{

  double * z;		/**< table containing z values used in this module */
  int z_size;

  double * k_nl;        /**< table of non-linear wavenumber at each redshift */

  //@}

  /** @name - tables of non-linear spectra and their second derivatives with respect to redshift */

  //@{

  double * p_density; /* density-density */
  double * p_cross; /* density-velocity */
  double * p_velocity; /* velocity-velocity */

  double * ddp_density;
  double * ddp_cross;
  double * ddp_velocity;

  //@}

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
			double z,
			double * pz_density,
			double * pz_velocity,
			double * pz_cross,
			int * k_size_at_z
			);

  int nonlinear_pk_at_k_and_z(
			      struct nonlinear * pnl,
			      double k,
			      double z,
			      double * pk_density,
			      double * pk_velocity,
			      double * pk_cross,
			      int * k_size_at_z
			      );

  int nonlinear_k_nl_at_z(
			  struct nonlinear * pnl,
			  double z,
			  double * k_nl
			);

  int nonlinear_init(
		     struct precision *ppr,
		     struct background *pba,
		     struct thermo *pth,
		     struct perturbs *ppt,
		     struct bessels * pbs,
		     struct transfers * ptr,
		     struct primordial *ppm,
		     struct spectra *psp,
		     struct nonlinear *pnl
		     );
  
  int nonlinear_free(
	       struct nonlinear *pnl
	       );

  int nonlinear_halofit(
			struct precision *ppr,
			struct background *pba,
			struct primordial *ppm,
			struct spectra *psp,
			struct nonlinear *pnl
			);

  
#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
