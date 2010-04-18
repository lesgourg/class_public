/** @file background.h Documented includes for background module */

#include "precision.h"
#include "growTable.h"
#include "tools_arrays.h"
#include "dei_rkck.h"

#ifndef __BACKGROUND__
#define __BACKGROUND__

/**
 * List of possible formats for the vector of background quantities.
 */
enum format_info {
  short_info, /**< compact format */
  long_info /**< exhaustive, redundent format */ 
};

/**
 * List of possible interpolation modes when calling the background_at_eta() interpolation function
 */
enum interpolation_mode {
  normal,   /**< default mode, OK for any situation */
  closeby /**< when the interpolation variable is arranged in growing order, and the new interpolation point is presumably close to the previous one */ 
};

/**
 * All background parameters and all background evolution.
 *
 * Once initialized by the backgound_init(), contains all necessary
 * information on the background evolution (excepted thermodynamics),
 * and in particular, a table of all background quantities as a
 * function of time and scale factor, used for interpolation in other
 * modules.
 */
struct background
{
  /** @name - background cosmological parameters.
   *
   * The background cosmological parameters listed here form a parameter
   * basis which is directly usable by the background module. Nothing
   * prevents from defining the input cosmological parameters
   * differently, and to pre-process them into this format, using the input
   * module (this might require iterative calls of background_init()
   * e.g. for dark energy or decaying dark matter). */

  //@{

  double H0; /**< \f$ H_0 \f$ : Hubble parameter in \f$ Mpc^{-1} \f$ */
  double Omega0_g; /**< \f$ \Omega_{0 \gamma} \f$ : photons */
  double Omega0_b; /**< \f$ \Omega_{0 b} \f$ : baryons */
  double Omega0_cdm; /**< \f$ \Omega_{0 cdm} \f$ : cold dark matter */
  double Omega0_lambda; /**< \f$ \Omega_{0_\Lambda} \f$ : cosmological constant */
  double Omega0_de; /**< \f$ \Omega_{0 de} \f$ : dark energy fluid with constant \f$ w \f$ */
  double w_de; /**< \f$ w_{DE} \f$ : dark energy equation of state */
  double cs2_de; /**< \f$ c^2_{s~DE} \f$ : dark energy sound speed */
  double Omega0_nur; /**< \f$ \Omega_{0 \nu r} \f$ : ultra-relativistic neutrinos */

  //@}


  /** @name - background cosmological parameters */

  //@{

  struct background_params * params; /**< a cosmo_params structure pointer */

  //@}

  /** @name - related parameters which can be computed only after the background integration */

  //@{

  double age; /**< age in Gyears */
  double conformal_age; /**< conformal age in Mpc */

  //@}

  /** @name - all indices for the vector of background (=bg) quantities */

  //@{

  int index_bg_a;             /**< scale factor */
  int index_bg_H;             /**< Hubble parameter in Mpc^{-1} */
  int index_bg_H_prime;       /**< its derivative w.r.t. conformal time */
  int index_bg_rho_g;         /**< photon density */
  int index_bg_rho_b;         /**< baryon density */
  int index_bg_rho_cdm;       /**< cdm density */
  int index_bg_rho_lambda;    /**< cosmological constant density */
  int index_bg_rho_de;        /**< dark energy fluid with constant w density */
  int index_bg_rho_nur;       /**< relativistic neutrinos/relics density */
  int index_bg_Omega_r;       /**< relativistic density fraction (\f$ \Omega_{\gamma} + \Omega_{\nu r} \f$) */
  int index_bg_rho_crit;      /**< critical density */
  int index_bg_Omega_g;       /**< Omega photons */
  int index_bg_Omega_b;       /**< Omega baryons */
  int index_bg_Omega_cdm;     /**< Omega cdm */
  int index_bg_Omega_lambda;  /**< Omega cosmological constant */
  int index_bg_Omega_de;      /**< Omega dark energy fluid with constant w */
  int index_bg_Omega_nur;     /**< Omega relativistic neutrinos/relics */
  int index_bg_conf_distance; /**< conformal distance (from us) */
  int index_bg_time;          /**< proper (cosmological) time */
  int index_bg_rs;            /**< comoving sound horizon */
  int bg_size_short; /**< size of background vector in the "short format" */
  int bg_size;       /**< size of background vector in the "long format" */

  //@}

  /** @name - background interpolation tables */

  //@{

  int bt_size; /**< number of lines (i.e. time-steps) in the array */
  double * eta_table; /**< values of \f$ \eta \f$ (conformal time) */
  double * z_table; /**< values of \f$ z \f$ (redshift) */
  double * background_table; /**< all other quantities (array of size bg_size*bt_size) **/

  //@}

  /** @name - table of their second derivatives, used for spline interpolation */

  //@{

  double * d2eta_dz2_table; /**< values of \f$ d^2 \eta / dz^2 \f$ (conformal time) */
  double * d2background_deta2_table; /**< values of \f$ d^2 b_i / d\eta^2 \f$ (conformal time) */

  //@}

  /** @name - all indices for the vector of background quantities to be integrated (=bi) 
   *
   * Most background quantities can be immediately inferred from the
   * scale factor. Only few of them require an integration with
   * respect to conformal time (in the minimal case, only one quantity needs to
   * be integrated with time: the scale factor, using the Friedmann
   * equation). These indices refer to the vector of
   * quantities to be integrated with time.
   */

  //@{

  int index_bi_a;    /**< scale factor */
  int index_bi_time; /**< proper (cosmological) time in Mpc */
  int index_bi_rs;   /**< sound horizon */
  int index_bi_eta;  /**< conformal time in Mpc */
  int bi_size;       /**< size of vector of background quantities to be integrated */

  //@}

  /** @name - flags describing the absence or presence of cosmological
      ingredients
      *
      * having one of these flag set to zero allows to skip the
      * corresponding contributions, instead of adding null contributions.
      */


  //@{

  short has_cdm;               /**< presence of cdm? */
  short has_lambda;            /**< presence of cosmological constant? */
  short has_dark_energy_fluid; /**< presence of dark energy fluid with constant w? */
  short has_nur;               /**< presence of ultra-relativistic neutrinos/relics? */

  //@}

  /** @name - flag regulating the amount of information sent to standard output (none if set to zero) */

  //@{

  short background_verbose;

  //@}

  /** @name - zone for writing error messages */

  //@{

  ErrorMsg error_message;

  //@}
};

/**************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int background_at_eta(
			double eta,
			enum format_info return_format,
			enum interpolation_mode intermode,
			int * last_index,
			double * pvecback_local
			);

  int background_functions_of_a(
				double a,
				enum format_info return_format,
				double * pvecback_local
				);

  int background_eta_of_z(
			  double z,
			  double * eta
			  );

  int background_init(
		      struct precision *ppr_input,
		      struct background *pba_output
		      );

  int background_free();

  int background_indices();

  int background_solve();

  int background_initial_conditions(double * pvecback_integration);

  void background_derivs(
			 double z,
			 double * y,
			 double * dy
			 );
      
#ifdef __cplusplus
}
#endif

/**************************************************************/

/**  
 * @name Some limits imposed on cosmological parameter values:
 */

//@{

#define _H0_BIG_ 1./2997.9     /**< maximal \f$ H_0 \f$ in \f$ Mpc^{-1} (h=1.0) \f$ */
#define _H0_SMALL_ 0.3/2997.9  /**< minimal \f$ H_0 \f$ in \f$ Mpc^{-1} (h=0.3) \f$ */
#define _TOLERANCE_ON_CURVATURE_ 1.e-5 /**< if \f$ | \Omega_k | \f$ smaller than this, considered as flat */

//@}

#endif
