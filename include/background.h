/** @file background.h Documented includes for background module */

#ifndef __BACKGROUND__
#define __BACKGROUND__

#include "common.h"
#include "quadrature.h"
#include "growTable.h"
#include "arrays.h"
#include "dei_rkck.h"
#include "parser.h"

/**
 * List of possible formats for the vector of background quantities.
 */

enum format_info {
  short_info,  /**< compact format (when only a, H, H' should be returned) */
  normal_info, /**< normal format (needed when integrating over perturbations: same plus rho_i's and Omega_r) */
  long_info    /**< exhaustive format (same plus proper time, sound horizon, Omega_m, etc.) */ 
};

/**
 * List of possible interpolation modes when calling the background_at_eta() interpolation function
 */

enum interpolation_mode {
  normal,   /**< default mode, OK for any situation */
  closeby /**< when the interpolation variable is arranged in growing order, and the new interpolation point is presumably close to the previous one */ 
};

/**
 * All background parameters and evolution that other modules need to know.
 *
 * Once initialized by the backgound_init(), contains all necessary
 * information on the background evolution (excepted thermodynamics),
 * and in particular, a table of all background quantities as a
 * function of time and scale factor, used for interpolation in other
 * modules.
 */

struct background
{
  /** @name - input parameters initialized by user in input module
   *  (all other quantitites are computed in this module, given these parameters
   *   and the content of the 'precision' structure)
   *
   * The background cosmological parameters listed here form a parameter
   * basis which is directly usable by the background module. Nothing
   * prevents from defining the input cosmological parameters
   * differently, and to pre-process them into this format, using the input
   * module (this might require iterative calls of background_init()
   * e.g. for dark energy or decaying dark matter). */
  
  //@{

  double H0; /**< \f$ H_0 \f$ : Hubble parameter (in fact, [H_0/c]) in \f$ Mpc^{-1} \f$ */

  double Omega0_g; /**< \f$ \Omega_{0 \gamma} \f$ : photons */
  double Tcmb; /**< \f$ T_{cmb} \f$ : current CMB temperature in Kelvins */

  double Omega0_b; /**< \f$ \Omega_{0 b} \f$ : baryons */

  double Omega0_cdm; /**< \f$ \Omega_{0 cdm} \f$ : cold dark matter */

  double Omega0_lambda; /**< \f$ \Omega_{0_\Lambda} \f$ : cosmological constant */

  double Omega0_de; /**< \f$ \Omega_{0 de} \f$ : dark energy fluid with constant \f$ w \f$ */
  double w_de; /**< \f$ w_{DE} \f$ : dark energy equation of state */
  double cs2_de; /**< \f$ c^2_{s~DE} \f$ : dark energy sound speed */

  double Omega0_nur; /**< \f$ \Omega_{0 \nu r} \f$ : ultra-relativistic neutrinos */

  int N_ncdm;      /* Number of distinguishabe ncdm species */
  double *M_ncdm;  /* vector of masses of non-cold relic: m_ncdm1/T_ncdm1 */
  double *T_ncdm;   /* list of 1st parameters in p-s-d of non-cold relics: temperature T_ncdm1/T_gamma */
  double *ksi_ncdm; /* list of 2nd parameters in p-s-d of first non-cold relic: temperature ksi_ncdm1/T_ncdm1 */
  double *deg_ncdm; /* list of degeneracies of ncdm species: 1 for one family of neutrinos (= one neutrino plus its anti-neutrino, total g*=1+1=2 */
  double *Omega0_ncdm; /*list of contributions to Omega0_ncdm */
  double Omega0_ncdm_tot;

  double Omega0_k; /**< \f$ \Omega_{0_k} \f$ : curvature contribution */
  //@}

  /** @name - related parameters */

  //@{

  double h; /** reduced Hubble parameter */
  double age; /**< age in Gyears */
  double conformal_age; /**< conformal age in Mpc */
  double *m_ncdm_in_eV;

  //@}

  /** @name - other background parameters */

  //@{

  double a_today; /** scale factor today (arbitrary and irrelevant for most purposes) */

  //@}

  /** @name - all indices for the vector of background (=bg) quantities stored in table */

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

  int index_bg_rho_ncdm1;
  int index_bg_p_ncdm1;

  int index_bg_Omega_r;       /**< relativistic density fraction (\f$ \Omega_{\gamma} + \Omega_{\nu r} \f$) */
  int index_bg_rho_crit;      /**< critical density */
  int index_bg_Omega_m;       /**< non-relativistic density fraction (\f$ \Omega_b + \Omega_cdm + \Omega_{\nu nr} \f$) */
  int index_bg_conf_distance; /**< conformal distance (from us) */
  int index_bg_time;          /**< proper (cosmological) time */
  int index_bg_rs;            /**< comoving sound horizon */

  int bg_size_short;  /**< size of background vector in the "short format" */
  int bg_size_normal; /**< size of background vector in the "normal format" */
  int bg_size;        /**< size of background vector in the "long format" */

  //@}

  /** @name - background interpolation tables */

  //@{

  int bt_size; /**< number of lines (i.e. time-steps) in the array */
  double * eta_table; /**< vector eta_table[index_eta] with values of \f$ \eta \f$ (conformal time) */
  double * z_table; /**< vector z_table[index_eta] with values of \f$ z \f$ (redshift) */
  double * background_table; /**< table background_table[index_eta*pba->bg_size+pba->index_bg] with all other quantities (array of size bg_size*bt_size) **/

  //@}

  /** @name - table of their second derivatives, used for spline interpolation */

  //@{

  double * d2eta_dz2_table; /**< vector d2eta_dz2_table[index_eta] with values of \f$ d^2 \eta / dz^2 \f$ (conformal time) */
  double * d2background_deta2_table; /**< table d2background_deta2_table[index_eta*pba->bg_size+pba->index_bg] with values of \f$ d^2 b_i / d\eta^2 \f$ (conformal time) */

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
  short has_ncdm;

  //@}

  /** @name - arrays related to sampling and integration of ncdm phase-space ditribution function
   */
  

  //@{

  double ** q_ncdm_bg; /* Pointers to vectors of background sampling in q */
  double ** w_ncdm_bg; /* Pointers to vectors of corresponding weights w */
  double ** q_ncdm;    /* Pointers to vectors of perturbation sampling in q */
  double ** w_ncdm;    /* Pointers to vectors of corresponding weights w */
  double ** dlnf0_dlnq_ncdm; /* Pointers to vectors of logarithmic derivatives of p-d-f */
  int *q_size_ncdm_bg; /* Size of the q_ncdm_bg arrays */
  int *q_size_ncdm;    /* Size of the q_ncdm arrays */
  double *factor_ncdm; /* List of conversion factors for calculating energy density etc.*/

  /** @name - technical parameters */

  //@{

  short background_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */

  //@}
};

/**
 * temporary parameters and workspace passed to the background_derivs function 
 */

struct background_parameters_and_workspace {

  /* structures containing fixed input parameters (indices, ...) */
  struct background * pba; 

  /* workspace */
  double * pvecback;

};

/**
 * temporary parameters and workspace passed to distribution function 
 */

struct background_parameters_for_distributions {

  /* structures containing fixed input parameters (indices, ...) */
  struct background * pba; 

  /* Additional parameters */
  int n_ncdm; /* Current distribution function */

};

/**************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int background_at_eta(
			struct background *pba,
			double eta,
			enum format_info return_format,
			enum interpolation_mode intermode,
			int * last_index,
			double * pvecback
			);

  int background_functions(
			   struct background *pba,
			   double a,
			   enum format_info return_format,
			   double * pvecback
			   );

  int background_eta_of_z(
			  struct background *pba,
			  double z,
			  double * eta
			  );

  int background_init(
		      struct precision *ppr,
		      struct background *pba
		      );

  int background_free(
		      struct background *pba
		      );

  int background_indices(
			 struct background *pba
			 );

  int background_ncdm1_distribution(
				  void *pba,
				  double q,
				  double * f0
				  );

  int background_ncdm1_test_function(
				     void *pba,
				     double q,
				     double * test
				     );

  int background_ncdm_init(
			    struct precision *ppr,
			    struct background *pba
			    );
  

  int background_ncdm_momenta(
                             double * qvec,
                             double * wvec,
                             int qsize,
                             double M,
                             double factor,
                             double z,
                             double * n,
		             double * rho, /* [8piG/3c2] rho in Mpc^-2 */
                             double * p,   /* [8piG/3c2] p in Mpc^-2 */
                             double * drho_dM
                             );

  int background_ncdm_M_from_Omega(
				    struct precision *ppr,
				    struct background *pba,
					int species
				    );

  int background_solve(
		       struct precision *ppr,
		       struct background *pba
		       );

  int background_initial_conditions(
				    struct precision *ppr,
				    struct background *pba,
				    double * pvecback,
				    double * pvecback_integration
				    );

  int background_derivs(
			 double z,
			 double * y,
			 double * dy,
			 void * parameters_and_workspace,
			 ErrorMsg error_message
			 );
      
#ifdef __cplusplus
}
#endif

/**************************************************************/

/**  
 * @name Some conversion factors and fundamental constants needed by background module:
 */

//@{

/*
  #define _Mpc_over_m_ 3.085677581282e22 */ /**< conversion factor from meters to megaparsecs */


/* for testing, set like in CAMB, although less precise: */
#define _Mpc_over_m_ 3.085678e22

#define _Gyr_over_Mpc_ 3.06601394e2 /**< conversion factor from megaparsecs to gigayears 
				       (c=1 units, Julian years of 365.25 days) */
#define _c_ 2.99792458e8 /**< c in m/s */
#define _G_ 6.67428e-11 /**< Newton constant in m^3/Kg/s^2 */
#define _eV_ 1.602176487e-19 /**< 1 eV expressed in J */

/* parameters entering into the Stefan-Boltzmann constant sigma_B */
#define _k_B_ 1.3806504e-23
#define _h_P_ 6.62606896e-34
/* sigma_B = 2pi^5k_B^4/(15h^3c^2) = 5.670400e-8 = Stefan-Boltzmann constant in W/m^2/K^4 = Kg/K^4/s^3 */

//@}

/**  
 * @name Some limits imposed on cosmological parameter values:
 */

//@{

#define _H0_BIG_ 1./2997.9     /**< maximal \f$ H_0 \f$ in \f$ Mpc^{-1} (h=1.0) \f$ */
#define _H0_SMALL_ 0.3/2997.9  /**< minimal \f$ H_0 \f$ in \f$ Mpc^{-1} (h=0.3) \f$ */
#define _TCMB_BIG_ 2.8     /**< maximal \f$ T_{cmb} \f$ in K */
#define _TCMB_SMALL_ 2.7   /**< minimal \f$ T_{cmb}  \f$ in K */
#define _TOLERANCE_ON_CURVATURE_ 1.e-5 /**< if \f$ | \Omega_k | \f$ smaller than this, considered as flat */

//@}

#endif
