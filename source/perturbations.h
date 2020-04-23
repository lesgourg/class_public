/** @file perturbations.h Documented includes for perturbation module */

#ifndef __PERTURBATIONS__
#define __PERTURBATIONS__

#include "thermodynamics.h"
#include "evolver_ndf15.h"
#include "evolver_rkck.h"

#define _scalars_ ((ppt->has_scalars == _TRUE_) && (index_md == index_md_scalars_))
#define _vectors_ ((ppt->has_vectors == _TRUE_) && (index_md == index_md_vectors_))
#define _tensors_ ((ppt->has_tensors == _TRUE_) && (index_md == index_md_tensors_))

//TODO: Remove those when possible!
#define _scalarsEXT_ ((ppt->has_scalars == _TRUE_) && (index_md == perturbations_module_->index_md_scalars_))
#define _vectorsEXT_ ((ppt->has_vectors == _TRUE_) && (index_md == perturbations_module_->index_md_vectors_))
#define _tensorsEXT_ ((ppt->has_tensors == _TRUE_) && (index_md == perturbations_module_->index_md_tensors_))


#define _set_source_(index) sources_[index_md][index_ic*tp_size_[index_md] + index][index_tau*k_size_[index_md] + index_k]

/**
 * flags for various approximation schemes
 * (tca = tight-coupling approximation,
 *  rsa = radiation streaming approximation,
 *  ufa = massless neutrinos / ultra-relativistic relics fluid approximation)
 *
 * CAUTION: must be listed below in chronological order, and cannot be
 * reversible. When integrating equations for a given mode, it is only
 * possible to switch from left to right in the lists below.
 */

//@{

enum tca_flags {tca_on, tca_off};
enum rsa_flags {rsa_off, rsa_on};
enum tca_idm_dr_flags {tca_idm_dr_on, tca_idm_dr_off};
enum rsa_idr_flags {rsa_idr_off, rsa_idr_on};
enum ufa_flags {ufa_off, ufa_on};
enum ncdmfa_flags {ncdmfa_off, ncdmfa_on};

//@}

/**
 * labels for the way in which each approximation scheme is implemented
 */

//@{

enum tca_method {first_order_MB,first_order_CAMB,first_order_CLASS,second_order_CRS,second_order_CLASS,compromise_CLASS};
enum rsa_method {rsa_null,rsa_MD,rsa_MD_with_reio,rsa_none};
enum idr_method {idr_free_streaming,idr_fluid}; /* for the idm-idr case */
enum rsa_idr_method {rsa_idr_none,rsa_idr_MD};  /* for the idm-idr case */
enum ufa_method {ufa_mb,ufa_hu,ufa_CLASS,ufa_none};
enum ncdmfa_method {ncdmfa_mb,ncdmfa_hu,ncdmfa_CLASS,ncdmfa_none};
enum tensor_methods {tm_photons_only,tm_massless_approximation,tm_exact};

//@}

/**
 * List of coded gauges. More gauges can in principle be defined.
 */

//@{

enum possible_gauges {
  newtonian, /**< newtonian (or longitudinal) gauge */
  synchronous /**< synchronous gauge with \f$ \theta_{cdm} = 0 \f$ by convention */
};

//@}

//@{

/**
 * maximum number and types of selection function (for bins of matter density or cosmic shear)
 */
#define _SELECTION_NUM_MAX_ 100
enum selection_type {gaussian,tophat,dirac};

//@}


//@{

/**
 * maximum number of k-values for perturbation output
 */
#define _MAX_NUMBER_OF_K_FILES_ 30

//@}



/**
 * Structure containing everything about perturbations that other
 * modules need to know, in particular tabled values of the source
 * functions \f$ S(k, \tau) \f$ for all requested modes
 * (scalar/vector/tensor), initial conditions, types (temperature,
 * E-polarization, B-polarization, lensing potential, etc), multipole
 * l and wavenumber k.
 *
 */

struct perturbs
{
  /** @name - input parameters initialized by user in input module
   *  (all other quantities are computed in this module, given these
   *  parameters and the content of the 'precision', 'background' and
   *  'thermodynamics' structures) */

  //@{

  short has_perturbations; /**< do we need to compute perturbations at all ? */

  short has_cls; /**< do we need any harmonic space spectrum \f$ C_l \f$ (and hence Bessel functions, transfer functions, ...)? */

  short has_scalars; /**< do we need scalars? */
  short has_vectors; /**< do we need vectors? */
  short has_tensors; /**< do we need tensors? */

  short has_ad;      /**< do we need adiabatic mode? */
  short has_bi;      /**< do we need isocurvature bi mode? */
  short has_cdi;     /**< do we need isocurvature cdi mode? */
  short has_nid;     /**< do we need isocurvature nid mode? */
  short has_niv;     /**< do we need isocurvature niv mode? */

  /* perturbed recombination */
  /** Do we want to consider perturbed temperature and ionization fraction? */
  short has_perturbed_recombination;
  /** Neutrino contribution to tensors */
  enum tensor_methods tensor_method;  /**< way to treat neutrinos in tensor perturbations(neglect, approximate as massless, take exact equations) */

  short has_cl_cmb_temperature;       /**< do we need \f$ C_l \f$'s for CMB temperature? */
  short has_cl_cmb_polarization;      /**< do we need \f$ C_l \f$'s for CMB polarization? */
  short has_cl_cmb_lensing_potential; /**< do we need \f$ C_l \f$'s for CMB lensing potential? */
  short has_cl_lensing_potential;     /**< do we need \f$ C_l \f$'s for galaxy lensing potential? */
  short has_cl_number_count;          /**< do we need \f$ C_l \f$'s for density number count? */
  short has_pk_matter;                /**< do we need matter Fourier spectrum? */
  short has_density_transfers;        /**< do we need to output individual matter density transfer functions? */
  short has_velocity_transfers;       /**< do we need to output individual matter velocity transfer functions? */
  short has_metricpotential_transfers;/**< do we need to output individual transfer functions for scalar metric perturbations? */
  short has_Nbody_gauge_transfers;    /**< should we convert density and velocity transfer functions to Nbody gauge? */

  short has_nl_corrections_based_on_delta_m;  /**< do we want to compute non-linear corrections with an algorithm relying on delta_m (like halofit)? */

  short has_nc_density;  /**< in dCl, do we want density terms ? */
  short has_nc_rsd;      /**< in dCl, do we want redshift space distortion terms ? */
  short has_nc_lens;     /**< in dCl, do we want lensing terms ? */
  short has_nc_gr;       /**< in dCl, do we want gravity terms ? */

  int l_scalar_max; /**< maximum l value for CMB scalars \f$ C_l \f$'s */
  int l_vector_max; /**< maximum l value for CMB vectors \f$ C_l \f$'s */
  int l_tensor_max; /**< maximum l value for CMB tensors \f$ C_l \f$'s */
  int l_lss_max; /**< maximum l value for LSS \f$ C_l \f$'s (density and lensing potential in  bins) */
  double k_max_for_pk; /**< maximum value of k in 1/Mpc in P(k) (if \f$ C_l \f$'s also requested, overseeded by value kmax inferred from l_scalar_max if it is bigger) */

  int selection_num;                            /**< number of selection functions
                                                   (i.e. bins) for matter density \f$ C_l \f$'s */
  enum selection_type selection;                /**< type of selection functions */
  double selection_mean[_SELECTION_NUM_MAX_]; /**< centers of selection functions */
  double selection_width[_SELECTION_NUM_MAX_];  /**< widths of selection functions */

  int switch_sw;   /**< in temperature calculation, do we want to include the intrinsic temperature + Sachs Wolfe term? */
  int switch_eisw; /**< in temperature calculation, do we want to include the early integrated Sachs Wolfe term? */
  int switch_lisw; /**< in temperature calculation, do we want to include the late integrated Sachs Wolfe term? */
  int switch_dop;  /**< in temperature calculation, do we want to include the Doppler term? */
  int switch_pol;  /**< in temperature calculation, do we want to include the polarization-related term? */
  double eisw_lisw_split_z; /**< at which redshift do we define the cut between eisw and lisw ?*/

  int store_perturbations;  /**< Do we want to store perturbations? */
  int k_output_values_num;       /**< Number of perturbation outputs (default=0) */
  double k_output_values[_MAX_NUMBER_OF_K_FILES_];    /**< List of k values where perturbation output is requested. */

  double three_ceff2_ur;/**< 3 x effective squared sound speed for the ultrarelativistic perturbations */
  double three_cvis2_ur;/**< 3 x effective viscosity parameter for the ultrarelativistic perturbations */

  double z_max_pk; /**< when we compute only the matter spectrum / transfer functions, but not the CMB, we are sometimes interested to sample source functions at very high redshift, way before recombination. This z_max_pk will then fix the initial sampling time of the sources. */

  double * alpha_idm_dr; /**< Angular contribution to collisional term at l>=2 for idm_fr-idr */
  double * beta_idr;  /**< Angular contribution to collisional term at l>=2 for idr-idr */

  int idr_nature; /**< Nature of the interacting dark radiation (free streaming or fluid) */

  //@}

  /** @name - gauge in which to perform the calculation */

  //@{

  enum possible_gauges gauge; /**< gauge in which to perform this calculation */

  //@}


  /** @name - list of conformal time values in the source table
      (common to all modes and types) */

  //@{

  double selection_min_of_tau_min; /**< used in presence of selection functions (for matter density, cosmic shear...) */
  double selection_max_of_tau_max; /**< used in presence of selection functions (for matter density, cosmic shear...) */

  double selection_delta_tau; /**< used in presence of selection functions (for matter density, cosmic shear...) */

  double * selection_tau_min; /**< value of conformal time below which W(tau) is considered to vanish for each bin */
  double * selection_tau_max; /**< value of conformal time above which W(tau) is considered to vanish for each bin */
  double * selection_tau; /**< value of conformal time at the center of each bin */
  double * selection_function; /**< selection function W(tau), normalized to \f$ \int W(tau) dtau=1 \f$, stored in selection_function[bin*tau_size_+index_tau] */

  //@}

  /** @name - technical parameters */

  //@{

  short perturbations_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  //@}

};

/**
 * Structure containing the indices and the values of the perturbation
 * variables which are integrated over time (as well as their
 * time-derivatives). For a given wavenumber, the size of these
 * vectors changes when the approximation scheme changes.
 */

struct perturb_vector
{
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
  int index_pt_delta_idm_dr;/**< idm_dr density */
  int index_pt_theta_idm_dr;/**< idm_dr velocity */
  int index_pt_delta_dcdm; /**< dcdm density */
  int index_pt_theta_dcdm; /**< dcdm velocity */
  int index_pt_delta_fld;  /**< dark energy density in true fluid case */
  int index_pt_theta_fld;  /**< dark energy velocity in true fluid case */
  int index_pt_Gamma_fld;  /**< unique dark energy dynamical variable in PPF case */
  int index_pt_phi_scf;  /**< scalar field density */
  int index_pt_phi_prime_scf;  /**< scalar field velocity */
  int index_pt_delta_ur; /**< density of ultra-relativistic neutrinos/relics */
  int index_pt_theta_ur; /**< velocity of ultra-relativistic neutrinos/relics */
  int index_pt_shear_ur; /**< shear of ultra-relativistic neutrinos/relics */
  int index_pt_l3_ur;    /**< l=3 of ultra-relativistic neutrinos/relics */
  int l_max_ur;          /**< max momentum in Boltzmann hierarchy (at least 3) */
  int index_pt_delta_idr; /**< density of interacting dark radiation */
  int index_pt_theta_idr; /**< velocity of interacting dark radiation */
  int index_pt_shear_idr; /**< shear of interacting dark radiation */
  int index_pt_l3_idr;    /**< l=3 of interacting dark radiation */
  int l_max_idr;          /**< max momentum in Boltzmann hierarchy (at least 3) for interacting dark radiation */

/* perturbed recombination */
  int index_pt_perturbed_recombination_delta_temp;		/**< Gas temperature perturbation */
  int index_pt_perturbed_recombination_delta_chi;		/**< Inionization fraction perturbation */

  /** The index to the first Legendre multipole of the DR expansion. Not
      that this is not exactly the usual delta, see Kaplinghat et al.,
      astro-ph/9907388. */
  int index_pt_F0_dr;
  int l_max_dr;          /**< max momentum in Boltzmann hierarchy for dr) */
  int index_pt_psi0_ncdm1; /**< first multipole of perturbation of first ncdm species, Psi_0 */
  int N_ncdm;		/**< number of distinct non-cold-dark-matter (ncdm) species */
  int* l_max_ncdm;	/**< mutipole l at which Boltzmann hierarchy is truncated (for each ncdm species) */
  int* q_size_ncdm;	/**< number of discrete momenta (for each ncdm species) */

  int index_pt_eta;       /**< synchronous gauge metric perturbation eta*/
  int index_pt_phi;	      /**< newtonian gauge metric perturbation phi */
  int index_pt_hv_prime;  /**< vector metric perturbation h_v' in synchronous gauge */
  int index_pt_V;         /**< vector metric perturbation V in Newtonian gauge */

  int index_pt_gw;        /**< tensor metric perturbation h (gravitational waves) */
  int index_pt_gwdot;     /**< its time-derivative */
  int pt_size;            /**< size of perturbation vector */

  double * y;             /**< vector of perturbations to be integrated */
  double * dy;            /**< time-derivative of the same vector */

  int * used_in_sources; /**< boolean array specifying which
                            perturbations enter in the calculation of
                            source functions */

};


/**
 * Workspace containing, among other things, the value at a given time
 * of all background/perturbed quantities, as well as their indices.
 * There will be one such structure created for each mode
 * (scalar/.../tensor) and each thread (in case of parallel computing)
 */

struct perturb_workspace
{

  /** @name - all possible useful indices for those metric
      perturbations which are not integrated over time, but just
      inferred from Einstein equations. "_mt_" stands for "metric".*/

  //@{

  int index_mt_psi;           /**< psi in longitudinal gauge */
  int index_mt_phi_prime;     /**< (d phi/d conf.time) in longitudinal gauge */
  int index_mt_h_prime;       /**< h' (wrt conf. time) in synchronous gauge */
  int index_mt_h_prime_prime; /**< h'' (wrt conf. time) in synchronous gauge */
  int index_mt_eta_prime;     /**< eta' (wrt conf. time) in synchronous gauge */
  int index_mt_alpha;         /**< \f$ \alpha = (h' + 6 \eta') / (2 k^2) \f$ in synchronous gauge */
  int index_mt_alpha_prime;   /**< \f$ \alpha'\f$ wrt conf. time) in synchronous gauge */
  int index_mt_gw_prime_prime;/**< second derivative wrt conformal time of gravitational wave field, often called h */
  int index_mt_V_prime;       /**< derivative of Newtonian gauge vector metric perturbation V */
  int index_mt_hv_prime_prime;/**< Second derivative of Synchronous gauge vector metric perturbation \f$ h_v\f$ */
  int mt_size;                /**< size of metric perturbation vector */

  //@}

  /** @name - value at a given time of all background/perturbed
      quantities
  */

  //@{

  double * pvecback;          /**< background quantities */
  double * pvecthermo;        /**< thermodynamics quantities */
  double * pvecmetric;        /**< metric quantities */
  struct perturb_vector * pv; /**< pointer to vector of integrated
                                 perturbations and their
                                 time-derivatives */

  double delta_rho;		    /**< total density perturbation (gives delta Too) */
  double rho_plus_p_theta;	/**< total (rho+p)*theta perturbation (gives delta Toi) */
  double rho_plus_p_shear;	/**< total (rho+p)*shear (gives delta Tij) */
  double delta_p;		    /**< total pressure perturbation (gives Tii) */

  double rho_plus_p_tot;    /**< total (rho+p) (used to infer theta_tot from rho_plus_p_theta) */

  double gw_source;		    /**< stress-energy source term in Einstein's tensor equations (gives Tij[tensor]) */
  double vector_source_pi;	/**< first stress-energy source term in Einstein's vector equations */
  double vector_source_v;	/**< second stress-energy source term in Einstein's vector equations */

  double tca_shear_g;  /**< photon shear in tight-coupling approximation */
  double tca_slip;     /**< photon-baryon slip in tight-coupling approximation */
  double tca_shear_idm_dr;/**< interacting dark radiation shear in tight coupling appproximation */
  double rsa_delta_g;  /**< photon density in radiation streaming approximation */
  double rsa_theta_g;  /**< photon velocity in radiation streaming approximation */
  double rsa_delta_ur; /**< photon density in radiation streaming approximation */
  double rsa_theta_ur; /**< photon velocity in radiation streaming approximation */
  double rsa_delta_idr; /**< interacting dark radiation density in dark radiation streaming approximation */
  double rsa_theta_idr; /**< interacting dark radiation velocity in dark radiation streaming approximation */

  double * delta_ncdm;	/**< relative density perturbation of each ncdm species */
  double * theta_ncdm;	/**< velocity divergence theta of each ncdm species */
  double * shear_ncdm;	/**< shear for each ncdm species */

  double delta_m;	/**< relative density perturbation of all non-relativistic species */
  double theta_m;	/**< velocity divergence theta of all non-relativistic species */

  double delta_cb;       /**< relative density perturbation of only cdm and baryon */
  double theta_cb;       /**< velocity divergence theta of only cdm and baryon */

  double delta_rho_fld;        /**< density perturbation of fluid, not so trivial in PPF scheme */
  double delta_p_fld;          /**< pressure perturbation of fluid, very non-trivial in PPF scheme */
  double rho_plus_p_theta_fld; /**< velocity divergence of fluid, not so trivial in PPF scheme */
  double S_fld;                /**< S quantity sourcing Gamma_prime evolution in PPF scheme (equivalent to eq. 15 in 0808.3125) */
  double Gamma_prime_fld;      /**< Gamma_prime in PPF scheme (equivalent to eq. 14 in 0808.3125) */

  FILE * perturb_output_file; /**< filepointer to output file*/
  int index_ikout;            /**< index for output k value (when k_output_values is set) */

  //@}

  /** @name - indices useful for searching background/thermo quantities in tables */

  //@{

  short inter_mode;	/**< flag defining the method used for interpolation background/thermo quantities tables */

  int last_index_back;   /**< the background interpolation function background_at_tau() keeps memory of the last point called through this index */
  int last_index_thermo; /**< the thermodynamics interpolation function thermodynamics_at_z() keeps memory of the last point called through this index */

  //@}

  /** @name - approximations used at a given time */

  //@{

  int index_ap_tca; /**< index for tight-coupling approximation */
  int index_ap_rsa; /**< index for radiation streaming approximation */
  int index_ap_tca_idm_dr; /**< index for dark tight-coupling approximation (idm-idr) */
  int index_ap_rsa_idr; /**< index for dark radiation streaming approximation */
  int index_ap_ufa; /**< index for ur fluid approximation */
  int index_ap_ncdmfa; /**< index for ncdm fluid approximation */
  int ap_size;      /**< number of relevant approximations for a given mode */

  int * approx;     /**< array of approximation flags holding at a given time: approx[index_ap] */

  //@}

  /** @name - approximations used at a given time */

  //@{

  int max_l_max;    /**< maximum l_max for any multipole */
  double * s_l;     /**< array of freestreaming coefficients \f$ s_l = \sqrt{1-K*(l^2-1)/k^2} \f$*/

  //@}

};

#endif
