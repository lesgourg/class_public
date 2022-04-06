/** @file perturbations.h Documented includes for perturbation module */

#ifndef __PERTURBATIONS__
#define __PERTURBATIONS__

#include "thermodynamics.h"

#define _scalars_ ((ppt->has_scalars == _TRUE_) && (index_md == ppt->index_md_scalars))
#define _vectors_ ((ppt->has_vectors == _TRUE_) && (index_md == ppt->index_md_vectors))
#define _tensors_ ((ppt->has_tensors == _TRUE_) && (index_md == ppt->index_md_tensors))

#define _set_source_(index) ppt->sources[index_md][index_ic * ppt->tp_size[index_md] + index][index_tau * ppt->k_size[index_md] + index_k]

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

struct perturbations
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

  short evolve_tensor_ur;             /**< will we evolve ur tensor perturbations (either because we have ur species, or we have ncdm species with massless approximation) ? */
  short evolve_tensor_ncdm;             /**< will we evolve ncdm tensor perturbations (if we have ncdm species and we use the exact method) ? */

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
  double k_max_for_pk; /**< maximum value of k in 1/Mpc required for the output of P(k,z) and T(k,z) */

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

  /** @name - useful flags inferred from the ones above */

  //@{

  short has_cmb; /**< do we need CMB-related sources (temperature, polarization) ? */
  short has_lss; /**< do we need LSS-related sources (lensing potential, ...) ? */

  short has_idm_dr; /**< do we have idm-dr interactions? */
  short has_idm_soundspeed; /**< do we need to consider the dark matter sound speed in interaction models? */
  //@}

  /** @name - gauge in which to perform the calculation */

  //@{

  enum possible_gauges gauge; /**< gauge in which to perform this calculation */

  //@}

  /** @name - indices running on modes (scalar, vector, tensor) */

  //@{

  int index_md_scalars; /**< index value for scalars */
  int index_md_tensors; /**< index value for tensors */
  int index_md_vectors; /**< index value for vectors */

  int md_size; /**< number of modes included in computation */

  //@}

  /** @name - indices running on initial conditions (for scalars: ad, cdi, nid, niv; for tensors: only one) */

  //@{

  int index_ic_ad; /**< index value for adiabatic */
  int index_ic_cdi; /**< index value for CDM isocurvature */
  int index_ic_bi; /**< index value for baryon isocurvature */
  int index_ic_nid; /**< index value for neutrino density isocurvature */
  int index_ic_niv; /**< index value for neutrino velocity isocurvature */
  int index_ic_ten; /**< index value for unique possibility for tensors */

  int * ic_size;       /**< for a given mode, ic_size[index_md] = number of initial conditions included in computation */

  //@}

  /** @name - flags and indices running on types (temperature, polarization, lensing, ...) */

  //@{

  short has_source_t;          /**< do we need source for CMB temperature? */
  short has_source_p;          /**< do we need source for CMB polarization? */
  short has_source_delta_m;    /**< do we need source for delta of total matter? */
  short has_source_delta_cb;   /**< do we ALSO need source for delta of ONLY cdm and baryon? */
  short has_source_delta_tot;  /**< do we need source for delta total? */
  short has_source_delta_g;    /**< do we need source for delta of gammas? */
  short has_source_delta_b;    /**< do we need source for delta of baryons? */
  short has_source_delta_cdm;  /**< do we need source for delta of cold dark matter? */
  short has_source_delta_idm;  /**< do we need source for delta of interacting dark matter */
  short has_source_delta_idr;  /**< do we need source for delta of interacting dark radiation? */
  short has_source_delta_dcdm; /**< do we need source for delta of DCDM? */
  short has_source_delta_fld;  /**< do we need source for delta of dark energy? */
  short has_source_delta_scf;  /**< do we need source for delta from scalar field? */
  short has_source_delta_dr;   /**< do we need source for delta of decay radiation? */
  short has_source_delta_ur;   /**< do we need source for delta of ultra-relativistic neutrinos/relics? */
  short has_source_delta_ncdm; /**< do we need source for delta of all non-cold dark matter species (e.g. massive neutrinos)? */
  short has_source_theta_m;    /**< do we need source for theta of total matter? */
  short has_source_theta_cb;   /**< do we ALSO need source for theta of ONLY cdm and baryon? */
  short has_source_theta_tot;  /**< do we need source for theta total? */
  short has_source_theta_g;    /**< do we need source for theta of gammas? */
  short has_source_theta_b;    /**< do we need source for theta of baryons? */
  short has_source_theta_cdm;  /**< do we need source for theta of cold dark matter? */
  short has_source_theta_idm;  /**< do we need source for theta of interacting dark matter */
  short has_source_theta_idr;  /**< do we need source for theta of interacting dark radiation? */
  short has_source_theta_dcdm; /**< do we need source for theta of DCDM? */
  short has_source_theta_fld;  /**< do we need source for theta of dark energy? */
  short has_source_theta_scf;  /**< do we need source for theta of scalar field? */
  short has_source_theta_dr;   /**< do we need source for theta of ultra-relativistic neutrinos/relics? */
  short has_source_theta_ur;   /**< do we need source for theta of ultra-relativistic neutrinos/relics? */
  short has_source_theta_ncdm; /**< do we need source for theta of all non-cold dark matter species (e.g. massive neutrinos)? */
  short has_source_phi;        /**< do we need source for metric fluctuation phi? */
  short has_source_phi_prime;  /**< do we need source for metric fluctuation phi'? */
  short has_source_phi_plus_psi; /**< do we need source for metric fluctuation (phi+psi)? */
  short has_source_psi;        /**< do we need source for metric fluctuation psi? */
  short has_source_h;          /**< do we need source for metric fluctuation h? */
  short has_source_h_prime;    /**< do we need source for metric fluctuation h'? */
  short has_source_eta;        /**< do we need source for metric fluctuation eta? */
  short has_source_eta_prime;  /**< do we need source for metric fluctuation eta'? */
  short has_source_H_T_Nb_prime; /**< do we need source for metric fluctuation H_T_Nb'? */
  short has_source_k2gamma_Nb; /**< do we need source for metric fluctuation gamma in Nbody gauge? */


  /* remember that the temperature source function includes three
     terms that we call 0,1,2 (since the strategy in class v > 1.7 is
     to avoid the integration by part that would reduce the source to
     a single term) */
  int index_tp_t0; /**< index value for temperature (j=0 term) */
  int index_tp_t1; /**< index value for temperature (j=1 term) */
  int index_tp_t2; /**< index value for temperature (j=2 term) */
  int index_tp_p; /**< index value for polarization */
  int index_tp_delta_m; /**< index value for matter density fluctuation */
  int index_tp_delta_cb; /**< index value for delta cb */
  int index_tp_delta_tot; /**< index value for total density fluctuation */
  int index_tp_delta_g;   /**< index value for delta of gammas */
  int index_tp_delta_b;   /**< index value for delta of baryons */
  int index_tp_delta_cdm; /**< index value for delta of cold dark matter */
  int index_tp_delta_idm; /**< index value for delta of interacting dark matter */
  int index_tp_delta_dcdm;/**< index value for delta of DCDM */
  int index_tp_delta_fld;  /**< index value for delta of dark energy */
  int index_tp_delta_scf;  /**< index value for delta of scalar field */
  int index_tp_delta_dr; /**< index value for delta of decay radiation */
  int index_tp_delta_ur; /**< index value for delta of ultra-relativistic neutrinos/relics */
  int index_tp_delta_idr; /**< index value for delta of interacting dark radiation */
  int index_tp_delta_ncdm1; /**< index value for delta of first non-cold dark matter species (e.g. massive neutrinos) */
  int index_tp_perturbed_recombination_delta_temp;		/**< Gas temperature perturbation */
  int index_tp_perturbed_recombination_delta_chi;		/**< Inionization fraction perturbation */

  int index_tp_theta_m;     /**< index value for matter velocity fluctuation */
  int index_tp_theta_cb;    /**< index value for theta cb */
  int index_tp_theta_tot;   /**< index value for total velocity fluctuation */
  int index_tp_theta_g;     /**< index value for theta of gammas */
  int index_tp_theta_b;     /**< index value for theta of baryons */
  int index_tp_theta_cdm;   /**< index value for theta of cold dark matter */
  int index_tp_theta_dcdm;  /**< index value for theta of DCDM */
  int index_tp_theta_fld;   /**< index value for theta of dark energy */
  int index_tp_theta_scf;   /**< index value for theta of scalar field */
  int index_tp_theta_ur;    /**< index value for theta of ultra-relativistic neutrinos/relics */
  int index_tp_theta_idr;   /**< index value for theta of interacting dark radiation */
  int index_tp_theta_idm;   /**< index value for theta of interacting dark matter */
  int index_tp_theta_dr;    /**< index value for F1 of decay radiation */
  int index_tp_theta_ncdm1; /**< index value for theta of first non-cold dark matter species (e.g. massive neutrinos) */

  int index_tp_phi;          /**< index value for metric fluctuation phi */
  int index_tp_phi_prime;    /**< index value for metric fluctuation phi' */
  int index_tp_phi_plus_psi; /**< index value for metric fluctuation phi+psi */
  int index_tp_psi;          /**< index value for metric fluctuation psi */
  int index_tp_h;            /**< index value for metric fluctuation h */
  int index_tp_h_prime;      /**< index value for metric fluctuation h' */
  int index_tp_eta;          /**< index value for metric fluctuation eta */
  int index_tp_eta_prime;    /**< index value for metric fluctuation eta' */
  int index_tp_H_T_Nb_prime; /**< index value for metric fluctuation H_T_Nb' */
  int index_tp_k2gamma_Nb;   /**< index value for metric fluctuation gamma times k^2 in Nbody gauge */

  int * tp_size; /**< number of types tp_size[index_md] included in computation for each mode */

  //@}

  /** @name - list of k values for each mode */

  //@{

  int * k_size_cmb;  /**< k_size_cmb[index_md] number of k values used
                        for CMB calculations, requiring a fine
                        sampling in k-space */

  int * k_size_cl;  /**< k_size_cl[index_md] number of k values used
                       for non-CMB \f$ C_l \f$ calculations, requiring a coarse
                       sampling in k-space. */

  int k_size_pk;    /**< number of k values for the P(k,z) and T(k,z) output, not including possible additional values for non-linear corrections */

  int * k_size;     /**< k_size[index_md] = total number of k values,
                       including those needed for all C_l, P(k),
                       nonlinear corrections */

  double ** k;      /**< k[index_md][index_k] = list of values */

  double k_min;     /**< minimum value (over all modes) */
  double k_max;     /**< maximum value (over all modes) */

  //@}

  /** @name - list of conformal time values in the source table
      (common to all modes and types) */

  //@{

  double * tau_sampling;    /**< array of tau values */
  int tau_size;             /**< number of values in this array */

  double selection_min_of_tau_min; /**< used in presence of selection functions (for matter density, cosmic shear...) */
  double selection_max_of_tau_max; /**< used in presence of selection functions (for matter density, cosmic shear...) */

  double selection_delta_tau; /**< used in presence of selection functions (for matter density, cosmic shear...) */

  double * selection_tau_min; /**< value of conformal time below which W(tau) is considered to vanish for each bin */
  double * selection_tau_max; /**< value of conformal time above which W(tau) is considered to vanish for each bin */
  double * selection_tau; /**< value of conformal time at the center of each bin */
  double * selection_function; /**< selection function W(tau), normalized to \f$ \int W(tau) dtau=1 \f$, stored in selection_function[bin*ppt->tau_size+index_tau] */

  //@}

  /** @name - source functions interpolation table */

  //@{

  double *** sources; /**< Pointer towards the source interpolation table
                         sources[index_md]
                         [index_ic * ppt->tp_size[index_md] + index_tp]
                         [index_tau * ppt->k_size + index_k] */

  //@}

  /** @name - arrays related to the interpolation table for sources at late times, corresponding to z < z_max_pk (used for Fourier transfer function and spectra output) */

  //@{

  double * ln_tau;         /**< log of the arrau tau_sampling, covering only the
                               final time range required for the output of
                               Fourier transfer functions (used for interpolations) */
  int ln_tau_size;         /**< total number of values in this array */
  int index_ln_tau_pk;     /**< first index relevant for output of P(k,z) and T(k,z) */

  double *** late_sources; /**< Pointer towards the source interpolation table
                              late_sources[index_md]
                              [index_ic * ppt->tp_size[index_md] + index_tp]
                              [index_tau * ppt->k_size + index_k]
                              Note that this is not a replication of part of the sources table,
                              it is just poiting towards the same memory zone, at the place where the late_sources actually start */

  double *** ddlate_sources; /**< Pointer towards the splined source interpolation table with second derivatives with respect to time
                                ddlate_sources[index_md]
                                [index_ic * ppt->tp_size[index_md] + index_tp]
                                [index_tau * ppt->k_size + index_k] */

  //@}

  /** @name - arrays storing the evolution of all sources for given k values, passed as k_output_values */

  //@{

  int * index_k_output_values; /**< List of indices corresponding to k-values close to k_output_values for each mode.
                                  index_k_output_values[index_md*k_output_values_num+ik]*/

  char scalar_titles[_MAXTITLESTRINGLENGTH_]; /**< _DELIMITER_ separated string of titles for scalar perturbation output files. */
  char vector_titles[_MAXTITLESTRINGLENGTH_]; /**< _DELIMITER_ separated string of titles for vector perturbation output files. */
  char tensor_titles[_MAXTITLESTRINGLENGTH_]; /**< _DELIMITER_ separated string of titles for tensor perturbation output files. */

  int number_of_scalar_titles; /**< number of titles/columns in scalar perturbation output files */
  int number_of_vector_titles; /**< number of titles/columns in vector perturbation output files*/
  int number_of_tensor_titles; /**< number of titles/columns in tensor perturbation output files*/

  double * scalar_perturbations_data[_MAX_NUMBER_OF_K_FILES_]; /**< Array of double pointers to perturbation output for scalars */
  double * vector_perturbations_data[_MAX_NUMBER_OF_K_FILES_]; /**< Array of double pointers to perturbation output for vectors */
  double * tensor_perturbations_data[_MAX_NUMBER_OF_K_FILES_]; /**< Array of double pointers to perturbation output for tensors */

  int size_scalar_perturbation_data[_MAX_NUMBER_OF_K_FILES_]; /**< Array of sizes of scalar double pointers  */
  int size_vector_perturbation_data[_MAX_NUMBER_OF_K_FILES_]; /**< Array of sizes of vector double pointers  */
  int size_tensor_perturbation_data[_MAX_NUMBER_OF_K_FILES_]; /**< Array of sizes of tensor double pointers  */

  //@}

  /** @name - technical parameters */

  //@{

  short perturbations_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */

  //@}

};

/**
 * Structure containing the indices and the values of the perturbation
 * variables which are integrated over time (as well as their
 * time-derivatives). For a given wavenumber, the size of these
 * vectors changes when the approximation scheme changes.
 */

struct perturbations_vector
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
  int index_pt_delta_idm; /**< idm density */
  int index_pt_theta_idm; /**< idm velocity */
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

struct perturbations_workspace
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
  struct perturbations_vector * pv; /**< pointer to vector of integrated
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
  double tca_shear_idm_dr; /**< interacting dark radiation shear in tight coupling appproximation */
  double rsa_delta_g;  /**< photon density in radiation streaming approximation */
  double rsa_theta_g;  /**< photon velocity in radiation streaming approximation */
  double rsa_delta_ur; /**< photon density in radiation streaming approximation */
  double rsa_theta_ur; /**< photon velocity in radiation streaming approximation */
  double rsa_delta_idr; /**< interacting dark radiation density in dark radiation streaming approximation */
  double rsa_theta_idr; /**< interacting dark radiation velocity in dark radiation streaming approximation */

  double theta_idm; /**< interacting dark matter velocity */
  double theta_idm_prime; /**< derivative of interacting dark matter velocity in regard to conformal time */

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

  FILE * perturbations_output_file; /**< filepointer to output file*/
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

/**
 * Structure pointing towards all what the function that perturbations_derivs
 * needs to know: fixed input parameters and indices contained in the
 * various structures, workspace, etc.
 */

struct perturbations_parameters_and_workspace {

  struct precision * ppr;         /**< pointer to the precision structure */
  struct background * pba;        /**< pointer to the background structure */
  struct thermodynamics * pth;            /**< pointer to the thermodynamics structure */
  struct perturbations * ppt;          /**< pointer to the precision structure */
  int index_md;                   /**< index of mode (scalar/.../vector/tensor) */
  int index_ic;			          /**< index of initial condition (adiabatic/isocurvature(s)/...) */
  int index_k;			          /**< index of wavenumber */
  double k;			              /**< current value of wavenumber in 1/Mpc */
  struct perturbations_workspace * ppw; /**< workspace defined above */

};

/*************************************************************************************************************/
/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int perturbations_sources_at_tau(
                                   struct perturbations * ppt,
                                   int index_md,
                                   int index_ic,
                                   int index_tp,
                                   double tau,
                                   double * psource_at_tau
                                   );

  int perturbations_sources_at_z(
                                 struct background * pba,
                                 struct perturbations * ppt,
                                 int index_md,
                                 int index_ic,
                                 int index_tp,
                                 double z,
                                 double * psource_at_z
                                 );

   int perturbations_sources_at_k_and_z(
                                        struct background * pba,
                                        struct perturbations * ppt,
                                        int index_md,
                                        int index_ic,
                                        int index_tp,
                                        double k,
                                        double z,
                                        double * psource_at_k_and_z
                                        );

  int perturbations_output_data_at_z(
                                     struct background * pba,
                                     struct perturbations * ppt,
                                     enum file_format output_format,
                                     double z,
                                     int number_of_titles,
                                     double *data
                                     );

  int perturbations_output_data_at_index_tau(
                                             struct background * pba,
                                             struct perturbations * ppt,
                                             enum file_format output_format,
                                             int index_tau,
                                             int number_of_titles,
                                             double *data
                                             );

  int perturbations_output_data(
                                struct background * pba,
                                struct perturbations * ppt,
                                enum file_format output_format,
                                double * tkfull,
                                int number_of_titles,
                                double *data
                                );

  int perturbations_output_titles(
                                  struct background *pba,
                                  struct perturbations *ppt,
                                  enum file_format output_format,
                                  char titles[_MAXTITLESTRINGLENGTH_]
                                  );

  int perturbations_output_firstline_and_ic_suffix(
                                                   struct perturbations *ppt,
                                                   int index_ic,
                                                   char first_line[_LINE_LENGTH_MAX_],
                                                   char ic_suffix[_SUFFIXNAMESIZE_]
                                                   );

  int perturbations_init(
                         struct precision * ppr,
                         struct background * pba,
                         struct thermodynamics * pth,
                         struct perturbations * ppt
                         );

  int perturbations_free_input(
                               struct perturbations * ppt
                               );

  int perturbations_free(
                         struct perturbations * ppt
                         );

  int perturbations_indices(
                            struct precision * ppr,
                            struct background * pba,
                            struct thermodynamics * pth,
                            struct perturbations * ppt
                            );

  int perturbations_timesampling_for_sources(
                                             struct precision * ppr,
                                             struct background * pba,
                                             struct thermodynamics * pth,
                                             struct perturbations * ppt
                                             );
  int perturbations_get_k_list(
                               struct precision * ppr,
                               struct background * pba,
                               struct thermodynamics * pth,
                               struct perturbations * ppt
                               );

  int perturbations_workspace_init(
                                   struct precision * ppr,
                                   struct background * pba,
                                   struct thermodynamics * pth,
                                   struct perturbations * ppt,
                                   int index_md,
                                   struct perturbations_workspace * ppw
                                   );

  int perturbations_workspace_free(
                                   struct perturbations * ppt,
                                   int index_md,
                                   struct perturbations_workspace * ppw
                                   );

  int perturbations_solve(
                          struct precision * ppr,
                          struct background * pba,
                          struct thermodynamics * pth,
                          struct perturbations * ppt,
                          int index_md,
                          int index_ic,
                          int index_k,
                          struct perturbations_workspace * ppw
                          );

  int perturbations_prepare_k_output(
                                     struct background * pba,
                                     struct perturbations * ppt
                                     );

  int perturbations_find_approximation_number(
                                              struct precision * ppr,
                                              struct background * pba,
                                              struct thermodynamics * pth,
                                              struct perturbations * ppt,
                                              int index_md,
                                              double k,
                                              struct perturbations_workspace * ppw,
                                              double tau_ini,
                                              double tau_end,
                                              int * interval_number,
                                              int * interval_number_of
                                              );

  int perturbations_find_approximation_switches(
                                                struct precision * ppr,
                                                struct background * pba,
                                                struct thermodynamics * pth,
                                                struct perturbations * ppt,
                                                int index_md,
                                                double k,
                                                struct perturbations_workspace * ppw,
                                                double tau_ini,
                                                double tau_end,
                                                double precision,
                                                int interval_number,
                                                int * interval_number_of,
                                                double * interval_limit,
                                                int ** interval_approx
                                                );

  int perturbations_vector_init(
                                struct precision * ppr,
                                struct background * pba,
                                struct thermodynamics * pth,
                                struct perturbations * ppt,
                                int index_md,
                                int index_ic,
                                double k,
                                double tau,
                                struct perturbations_workspace * ppw,
                                int * pa_old
                                );

  int perturbations_vector_free(
                                struct perturbations_vector * pv
                                );

  int perturbations_initial_conditions(
                                       struct precision * ppr,
                                       struct background * pba,
                                       struct perturbations * ppt,
                                       int index_md,
                                       int index_ic,
                                       double k,
                                       double tau,
                                       struct perturbations_workspace * ppw
                                       );

  int perturbations_approximations(
                                   struct precision * ppr,
                                   struct background * pba,
                                   struct thermodynamics * pth,
                                   struct perturbations * ppt,
                                   int index_md,
                                   double k,
                                   double tau,
                                   struct perturbations_workspace * ppw
                                   );

  int perturbations_timescale(
                              double tau,
                              void * parameters_and_workspace,
                              double * timescale,
                              ErrorMsg error_message
                              );

  int perturbations_einstein(
                             struct precision * ppr,
                             struct background * pba,
                             struct thermodynamics * pth,
                             struct perturbations * ppt,
                             int index_md,
                             double k,
                             double tau,
                             double * y,
                             struct perturbations_workspace * ppw
                             );

  int perturbations_total_stress_energy(
                                        struct precision * ppr,
                                        struct background * pba,
                                        struct thermodynamics * pth,
                                        struct perturbations * ppt,
                                        int index_md,
                                        double k,
                                        double * y,
                                        struct perturbations_workspace * ppw
                                        );

  int perturbations_sources(
                            double tau,
                            double * pvecperturbations,
                            double * pvecderivs,
                            int index_tau,
                            void * parameters_and_workspace,
                            ErrorMsg error_message
                            );

  int perturbations_print_variables(
                                    double tau,
                                    double * y,
                                    double * dy,
                                    void * parameters_and_workspace,
                                    ErrorMsg error_message
                                    );

  int perturbations_derivs(
                           double tau,
                           double * y,
                           double * dy,
                           void * parameters_and_workspace,
                           ErrorMsg error_message
                           );

  int perturbations_tca_slip_and_shear(
                                       double * y,
                                       void * parameters_and_workspace,
                                       ErrorMsg error_message
                                       );

  int perturbations_rsa_delta_and_theta(
                                        struct precision * ppr,
                                        struct background * pba,
                                        struct thermodynamics * pth,
                                        struct perturbations * ppt,
                                        double k,
                                        double * y,
                                        double a_prime_over_a,
                                        double * pvecthermo,
                                        struct perturbations_workspace * ppw,
                                        ErrorMsg error_message
                                        );

  int perturbations_rsa_idr_delta_and_theta(
                                            struct precision * ppr,
                                            struct background * pba,
                                            struct thermodynamics * pth,
                                            struct perturbations * ppt,
                                            double k,
                                            double * y,
                                            double a_prime_over_a,
                                            double * pvecthermo,
                                            struct perturbations_workspace * ppw,
                                            ErrorMsg error_message
                                            );

#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
/* @endcond */
