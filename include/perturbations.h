/** @file perturbations.h Documented includes for perturbation module */

#ifndef __PERTURBATIONS__
#define __PERTURBATIONS__

#include "thermodynamics.h"

/**  
 * @name Flag values for tight-coupling approximation
 */

//@{

enum tca_flags {tca_off, tca_on};

//@}

/**  
 * @name Flag values for switching on/off the integration of high momenta l>=2 for radiation (photons, massless neutrinos...).  At late time and inside the Hubble radius, theses perturbations vanish according to the free-streaming approximation.
 */

//@{

enum rp_flags {rp_on, rp_off};

//@}

/**
 * All perturbation-related parameters and tables of source functions \f$ S^{X} (k, \eta) \f$.
 *
 * Once initialized by perturb_init(), contains all the necessary
 * information on the perturbations, and in particular, all possible
 * indices and flags, and a table of source functions used for
 * intesrpolation in other modules.
 */
struct perturbs
{
  /** @name - all possible flags stating which perturbations should be computed */

  //@{

  short has_scalars; /**< do we need scalars? */
  short has_vectors; /**< do we need vectors? */
  short has_tensors; /**< do we need tensors? */

  short has_ad;      /**< do we need adiabatic mode? */
  short has_bi;      /**< do we need isocurvature bi mode? */
  short has_cdi;     /**< do we need isocurvature cdi mode? */
  short has_nid;     /**< do we need isocurvature nid mode? */
  short has_niv;     /**< do we need isocurvature niv mode? */

  short has_cl_cmb_temperature;       /**< do we need Cl's for CMB temperature? */
  short has_cl_cmb_polarization;      /**< do we need Cl's for CMB polarization? */
  short has_cl_cmb_lensing_potential; /**< do we need Cl's for CMB lensing potential? */
  short has_pk_matter;                /**< do we need matter Fourier spectrum? */

  short has_source_t;  /**< do we need source for CMB temperature? */
  short has_source_e;  /**< do we need source for CMB E-polarisation? */
  short has_source_b;  /**< do we need source for CMB B-polarisation? */
  short has_source_g;  /**< do we need source for gravitationnal potential? */

  //@}

  /** @name - index running on modes (scalar, vector, tensor) */

  //@{

  int index_md_scalars; /**< index value for scalars */
  int index_md_tensors; /**< index value for tensors */
  int index_md_vectors; /**< index value for vectors */
  int md_size; /**< number of modes included in computation */

  //@}

  /** @name - index running on types (temperature, polarization, lensing, ...) */

  //@{

  int index_tp_t; /**< index value for temperature */
  int index_tp_e; /**< index value for E-polarization */
  int index_tp_b; /**< index value for B-polarization */
  int index_tp_g; /**< index value for gravitationnal potential */
  int * tp_size; /**< number of types tp_size[index_mode] included in computation for each mode */

  //@}

  /** @name - index running on initial conditions (for scalars: ad, cdi, nid, niv; for tensors: only one) */

  //@{

  int index_ic_ad; /**< index value for adiabatic */
  int index_ic_cdi; /**< index value for CDM isocurvature */
  int index_ic_bi; /**< index value for baryon isocurvature */
  int index_ic_nid; /**< index value for neutrino density isocurvature */
  int index_ic_niv; /**< index value for neutrino velocity isocurvature */
  int index_ic_ten; /**< index value for unique possibility for tensors */
  int * ic_size;       /**< for a given mode, ic_size[index_mode] = number of initial conditions included in computation */

  //@}

  /** @name - list of k values for each mode */
  
  //@{

  int * k_size;     /**< k_size[index_mode] = number of values */
  int * k_size_cl;     /**< k_size_cl[index_mode] number of values to take into account in transfer functions for C_l spectra (could be smaller than k_size, e.g. for scalars if extra points needed in P(k) */
  double ** k;      /**< (k[index_mode])[index_k] = list of values */
  double k_scalar_kmax_for_pk; /**< maximum value of k in h/Mpc in P(k) (overseeded by value kmax inferred from k_scalar_oscillations if it is bigger) */

  //@}

  /** @name - list of conformal time values in the source table (for all modes and types) */

  //@{

  int eta_size;     /**< eta_size = number of values */
  double * eta_sampling;      /**< eta_sampling = list of eta values */

  //@}

  /** @name - source functions interpolation table */

  //@{

  double *** sources; /**< Pointer towards the source interpolation table ((sources[index_mode])[index_ic][index_type])[index_eta][index_k] */

  /* ((sources[index_mode])[index_ic][index_type])[index_eta][index_k] = for each mode, initial condition, type, time step and wavenumber, value of source function. ((sources[index_mode])[index_ic][index_type]) defined in the code as an array of size eta_size[index_type]*k_size[index_mode]. (sources[index_mode])  defined in the code as an array (of arrays) of size ic_size[index_mode]*tp_size. sources defined in the code as an array (of arrays of arrays) of size md_size. */

  //@}

  /** @name - flag regulating the amount of information sent to standard output (none if set to zero) */

  //@{

  short perturbations_verbose;

  //@}

  /** @name - zone for writing error messages */

  //@{

  ErrorMsg error_message;

  //@}

};

/**
 * All indices in the vectors describing all perturbations at a given time.
 *
 * For each mode, the vectors of perturbations are
 * different: hence, there will be one such structure created for each
 * mode. This structure is not used in any other modules. 
 */
struct perturb_workspace 
{
  /** @name - all possible useful indices for the vector of perturbed quantities to be integrated over time. "_pt_" stands for "perturbation". */

  //@{

  int index_pt_delta_g;   /**< photon density */
  int index_pt_theta_g;   /**< photon velocity */
  int index_pt_shear_g;   /**< photon shear */
  int index_pt_l3_g;      /**< photon l=3 */
  int l_max_g;            /**< max momentum in Boltzmann hierarchy (at least 3) */
  int index_pt_pol0_g;    /**< photon polarization, l=0 */
  int index_pt_pol1_g;    /**< photon polarization, l=1 */
  int index_pt_pol2_g;    /**< photon polarization, l=2 */
  int index_pt_pol3_g;    /**< photon polarization, l=3 */
  int l_max_pol_g;        /**< max momentum in Boltzmann hierarchy (at least 3) */
  int index_pt_delta_b;   /**< baryon density */
  int index_pt_theta_b;   /**< baryon velocity */
  int index_pt_delta_cdm; /**< cdm density */
  int index_pt_theta_cdm; /**< cdm velocity */
  int index_pt_delta_de; /**< dark energy density */
  int index_pt_theta_de; /**< dark energy velocity */
  int index_pt_delta_nur; /**< density of ultra-relativistic neutrinos/relics */
  int index_pt_theta_nur; /**< velocity of ultra-relativistic neutrinos/relics */
  int index_pt_shear_nur; /**< shear of ultra-relativistic neutrinos/relics */
  int index_pt_l3_nur;    /**< l=3 of ultra-relativistic neutrinos/relics */
  int l_max_nur;          /**< max momentum in Boltzmann hierarchy (at least 3) */
  int index_pt_eta;       /**< synchronous gauge metric perturbation eta*/
  int index_pt_gw;        /**< tensor metric perturbation h (gravitational waves) */
  int index_pt_gwdot;     /**< its time-derivative */
  int pt_size;            /**< size of perturbation vector */

  //@}

  /** @name - all possible useful indices for metric perturbations not to be integrated over time, but to be inferred from Einstein equations. "_mt_" stands for "metric".*/

  //@{

  int index_mt_phi; /**< phi in longitudinal gauge */
  int index_mt_psi; /**< psi in longitudinal gauge */
  int index_mt_phi_prime; /**< (d phi/d conf.time) in longitudinal gauge */
  int index_mt_h_prime; /**< (d h/d conf.time) in synchronous gauge */
  int index_mt_eta_prime; /**< (d \f$ \eta \f$/d conf.time) in synchronous gauge */
  int index_mt_alpha_prime; /**< (d \f$ \alpha \f$/d conf.time) in synchronous gauge, where \f$ \alpha = (h' + 6 \eta') / (2 k^2) \f$ */
  int mt_size; /**< size of metric perturbation vector */

  //@}

  /** @name - all possible useful indices for terms contributing to source function S=S0+S1'+S2'', not to be integrated over time, but to be inferred from "_pt_" and "_mt_" vectors. "_st_" stands for "source term".
   */

  //@{

  int index_st_eta;    /**< conformal time */
  int index_st_S0;     /**< first piece S0 */
  int index_st_S1;     /**< second piece S1 */
  int index_st_S2;     /**< third piece S2 */
  int index_st_dS1;     /**< derivative S1' */
  int index_st_dS2;     /**< derivative S2' */
  int index_st_ddS2;     /**< derivative S2'' */
  int st_size; /**< size of this vector */ 
  
  double * pvecback;
  double * pvecthermo;
  double * pvecperturbations; /**< vector of perturbations to be integrated, used throughout the perturbation module */
  double * pvecderivs; 
  double * pvecmetric;
  double * pvecsource_terms;

  /* table of source terms for each mode, initial condition and wavenumber: (source_terms_table[index_type])[index_eta][index_st] */
  double ** source_term_table;

  int * last_index_back;  /**< the background interpolation function background_at_eta() keeps memory of the last point called through this index */
  int * last_index_thermo; /**< the thermodynamics interpolation function thermodynamics_at_z() keeps memory of the last point called through this index */

  enum tca_flags tca; /**< flag for tight-coupling approximation */
  enum rp_flags rp; /**< flag for free-streaming approximation (switch on/off radiation perturbations) */

};

struct perturb_parameters_and_workspace {

  /* fixed input parameters */
  struct precision * ppr;
  struct background * pba;
  struct thermo * pth;
  struct perturbs * ppt;
  int index_mode;
  double k;

  /* workspace */
  struct perturb_workspace * ppw;
  
};

/*************************************************************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
  extern "C" {
#endif

    int perturb_sources_at_eta(
			       struct perturbs * ppt,
			       int index_mode,
			       int index_ic,
			       int index_k,
			       int index_type,
			       double eta,
			       double * pvecsources
			       );

    int perturb_init(
		     struct precision * ppr,
		     struct background * pba,
		     struct thermo * pth,
		     struct perturbs * ppt
		     );

    int perturb_free(
		     struct perturbs * ppt
		     );

    int perturb_indices_of_perturbs(
				    struct precision * ppr,
				    struct background * pba,
				    struct thermo * pth,
				    struct perturbs * ppt
				    );

    int perturb_timesampling_for_sources(
					 struct precision * ppr,
					 struct background * pba,
					 struct thermo * pth,
					 struct perturbs * ppt
					 );
    
    int perturb_get_k_list(
			   struct precision * ppr,
			   struct background * pba,
			   struct thermo * pth,
			   struct perturbs * ppt,
			   int index_mode);

    int perturb_workspace_init(
			       struct precision * ppr,
			       struct background * pba,
			       struct thermo * pth,
			       struct perturbs * ppt,
			       int index_mode,
			       struct perturb_workspace * ppw
			       );

    int perturb_workspace_free(
			       struct perturbs * ppt,
			       int index_mode,
			       struct perturb_workspace * ppw
			       );

    int perturb_solve(
		      struct precision * ppr,
		      struct background * pba,
		      struct thermo * pth,
		      struct perturbs * ppt,
		      int index_mode,
		      int index_ic,
		      int index_k,
		      struct perturb_workspace * ppw
		      );

    int perturb_initial_conditions(
				   struct precision * ppr,
				   struct background * pba,
				   struct perturbs * ppt,
				   int index_mode,
				   int index_ic,
				   double k,
				   double eta,
				   struct perturb_workspace * ppw
				   );

    int perturb_timescale_and_approximations(
					     struct precision * ppr,
					     struct background * pba,
					     struct thermo * pth,
					     struct perturbs * ppt,
					     int index_mode,
					     double k,
					     double eta,
					     enum interpolation_mode intermode,
					     struct perturb_workspace * ppw,
					     double * timescale
					     );

    int perturb_einstein(
			 struct precision * ppr,
			 struct background * pba,
			 struct perturbs * ppt,
			 int index_mode,
			 double k,
			 double eta,
			 double * y,
			 struct perturb_workspace * ppw
			 );

    int perturb_source_terms(
			     double eta,
			     struct perturb_parameters_and_workspace * pppaw
			     );

    int perturb_sources(
			struct perturbs * ppt,
			int index_mode,
			int index_ic,
			int index_k,
			struct perturb_workspace * ppw
			);

    int perturb_derivs(
			double eta,
			double * y,
			double * dy,
			void * parameters_and_workspace,
			ErrorMsg error_message
			);

#ifdef __cplusplus
  }
#endif

/**************************************************************/

#endif
