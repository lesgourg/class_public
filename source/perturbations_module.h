#ifndef PERTURBATIONS_MODULE_H
#define PERTURBATIONS_MODULE_H

#include "input.h"
#include "base_module.h"

class PerturbationsModule : public BaseModule {
public:
  PerturbationsModule(const Input& input, BackgroundModulePtr background_module, ThermodynamicsModulePtr thermodynamics_module);
  ~PerturbationsModule();
  int perturb_output_data(enum file_format output_format, double z, int number_of_titles, double* data) const;
  int perturb_output_titles(enum file_format output_format, char titles[_MAXTITLESTRINGLENGTH_]) const;
  int perturb_output_firstline_and_ic_suffix(int index_ic, char first_line[_LINE_LENGTH_MAX_], FileName ic_suffix) const;

  /** @name - indices running on modes (scalar, vector, tensor) */
  //@{
  int index_md_scalars_; /**< index value for scalars */
  int index_md_tensors_; /**< index value for tensors */
  int index_md_vectors_; /**< index value for vectors */
  int md_size_; /**< number of modes included in computation */
  //@}
  /** @name - indices running on initial conditions (for scalars: ad, cdi, nid, niv; for tensors: only one) */
  //@{
  int index_ic_ad_;  /**< index value for adiabatic */
  int index_ic_cdi_; /**< index value for CDM isocurvature */
  int index_ic_bi_;  /**< index value for baryon isocurvature */
  int index_ic_nid_; /**< index value for neutrino density isocurvature */
  int index_ic_niv_; /**< index value for neutrino velocity isocurvature */
  int index_ic_ten_; /**< index value for unique possibility for tensors */
  int* ic_size_;     /**< for a given mode, ic_size[index_md] = number of initial conditions included in computation */
  //@}
  /** @name - flags and indices running on types (temperature, polarization, lensing, ...) */

  //@{

  /* remember that the temperature source function includes three
     terms that we call 0,1,2 (since the strategy in class v > 1.7 is
     to avoid the integration by part that would reduce the source to
     a single term) */
  int index_tp_t0_; /**< index value for temperature (j=0 term) */
  int index_tp_t1_; /**< index value for temperature (j=1 term) */
  int index_tp_t2_; /**< index value for temperature (j=2 term) */
  int index_tp_p_; /**< index value for polarization */
  int index_tp_delta_tot_; /**< index value for total density fluctuation */
  int index_tp_delta_g_;   /**< index value for delta of gammas */
  int index_tp_delta_b_;   /**< index value for delta of baryons */
  int index_tp_delta_cdm_; /**< index value for delta of cold dark matter */
  int index_tp_delta_dcdm_;/**< index value for delta of DCDM */
  int index_tp_delta_fld_;  /**< index value for delta of dark energy */
  int index_tp_delta_scf_;  /**< index value for delta of scalar field */
  int index_tp_delta_dr_; /**< index value for delta of decay radiation */
  int index_tp_delta_ur_; /**< index value for delta of ultra-relativistic neutrinos/relics */
  int index_tp_delta_idr_; /**< index value for delta of interacting dark radiation */
  int index_tp_delta_idm_dr_;/**< index value for delta of interacting dark matter (with dr)*/
  int index_tp_delta_ncdm1_; /**< index value for delta of first non-cold dark matter species (e.g. massive neutrinos) */
  int index_tp_perturbed_recombination_delta_temp_;    /**< Gas temperature perturbation */
  int index_tp_perturbed_recombination_delta_chi_;    /**< Inionization fraction perturbation */

  int index_tp_theta_m_;     /**< index value for matter velocity fluctuation */
  int index_tp_theta_cb_;    /**< index value for theta cb */
  int index_tp_theta_tot_;   /**< index value for total velocity fluctuation */
  int index_tp_theta_g_;     /**< index value for theta of gammas */
  int index_tp_theta_b_;     /**< index value for theta of baryons */
  int index_tp_theta_cdm_;   /**< index value for theta of cold dark matter */
  int index_tp_theta_dcdm_;  /**< index value for theta of DCDM */
  int index_tp_theta_fld_;   /**< index value for theta of dark energy */
  int index_tp_theta_scf_;   /**< index value for theta of scalar field */
  int index_tp_theta_ur_;    /**< index value for theta of ultra-relativistic neutrinos/relics */
  int index_tp_theta_idr_;   /**< index value for theta of interacting dark radiation */
  int index_tp_theta_idm_dr_;/**< index value for theta of interacting dark matter (with dr)*/
  int index_tp_theta_dr_;    /**< index value for F1 of decay radiation */
  int index_tp_theta_ncdm1_; /**< index value for theta of first non-cold dark matter species (e.g. massive neutrinos) */

  int index_tp_phi_;          /**< index value for metric fluctuation phi */
  int index_tp_phi_prime_;    /**< index value for metric fluctuation phi' */
  int index_tp_phi_plus_psi_; /**< index value for metric fluctuation phi+psi */
  int index_tp_psi_;          /**< index value for metric fluctuation psi */
  int index_tp_h_;            /**< index value for metric fluctuation h */
  int index_tp_h_prime_;      /**< index value for metric fluctuation h' */
  int index_tp_eta_;          /**< index value for metric fluctuation eta */
  int index_tp_eta_prime_;    /**< index value for metric fluctuation eta' */
  int index_tp_H_T_Nb_prime_; /**< index value for metric fluctuation H_T_Nb' */
  int index_tp_k2gamma_Nb_;   /**< index value for metric fluctuation gamma times k^2 in Nbody gauge */

  int index_tp_delta_m_;      /**< index value for matter density fluctuation */
  int index_tp_delta_cb_;     /**< index value for delta cb */
  int* tp_size_;              /**< number of types tp_size[index_md] included in computation for each mode */

  short has_source_t_;           /**< do we need source for CMB temperature? */
  short has_source_p_;           /**< do we need source for CMB polarization? */
  short has_source_delta_m_;     /**< do we need source for delta of total matter? */
  short has_source_delta_cb_;    /**< do we ALSO need source for delta of ONLY cdm and baryon? */
  short has_source_delta_tot_;   /**< do we need source for delta total? */
  short has_source_delta_g_;     /**< do we need source for delta of gammas? */
  short has_source_delta_b_;     /**< do we need source for delta of baryons? */
  short has_source_delta_cdm_;   /**< do we need source for delta of cold dark matter? */
  short has_source_delta_dcdm_;  /**< do we need source for delta of DCDM? */
  short has_source_delta_fld_;   /**< do we need source for delta of dark energy? */
  short has_source_delta_scf_;   /**< do we need source for delta from scalar field? */
  short has_source_delta_dr_;    /**< do we need source for delta of decay radiation? */
  short has_source_delta_ur_;    /**< do we need source for delta of ultra-relativistic neutrinos/relics? */
  short has_source_delta_idr_;   /**< do we need source for delta of interacting dark radiation? */
  short has_source_delta_idm_dr_;/**< do we need source for delta of interacting dark matter (with dr)? */
  short has_source_delta_ncdm_;  /**< do we need source for delta of all non-cold dark matter species (e.g. massive neutrinos)? */
  short has_source_theta_m_;     /**< do we need source for theta of total matter? */
  short has_source_theta_cb_;    /**< do we ALSO need source for theta of ONLY cdm and baryon? */
  short has_source_theta_tot_;   /**< do we need source for theta total? */
  short has_source_theta_g_;     /**< do we need source for theta of gammas? */
  short has_source_theta_b_;     /**< do we need source for theta of baryons? */
  short has_source_theta_cdm_;   /**< do we need source for theta of cold dark matter? */
  short has_source_theta_dcdm_;  /**< do we need source for theta of DCDM? */
  short has_source_theta_fld_;   /**< do we need source for theta of dark energy? */
  short has_source_theta_scf_;   /**< do we need source for theta of scalar field? */
  short has_source_theta_dr_;    /**< do we need source for theta of ultra-relativistic neutrinos/relics? */
  short has_source_theta_ur_;    /**< do we need source for theta of ultra-relativistic neutrinos/relics? */
  short has_source_theta_idr_;   /**< do we need source for theta of interacting dark radiation? */
  short has_source_theta_idm_dr_;/**< do we need source for theta of interacting dark matter (with dr)? */
  short has_source_theta_ncdm_;  /**< do we need source for theta of all non-cold dark matter species (e.g. massive neutrinos)? */
  short has_source_phi_;         /**< do we need source for metric fluctuation phi? */
  short has_source_phi_prime_;   /**< do we need source for metric fluctuation phi'? */
  short has_source_phi_plus_psi_;/**< do we need source for metric fluctuation (phi+psi)? */
  short has_source_psi_;         /**< do we need source for metric fluctuation psi? */
  short has_source_h_;           /**< do we need source for metric fluctuation h? */
  short has_source_h_prime_;     /**< do we need source for metric fluctuation h'? */
  short has_source_eta_;         /**< do we need source for metric fluctuation eta? */
  short has_source_eta_prime_;   /**< do we need source for metric fluctuation eta'? */
  short has_source_H_T_Nb_prime_;/**< do we need source for metric fluctuation H_T_Nb'? */
  short has_source_k2gamma_Nb_;  /**< do we need source for metric fluctuation gamma in Nbody gauge? */

  /** @name - arrays storing the evolution of all sources for given k values, passed as k_output_values */
  //@{
   int* index_k_output_values_;  /**< List of indices corresponding to k-values close to k_output_values for each mode.
                                     index_k_output_values[index_md*k_output_values_num+ik]*/

  char scalar_titles_[_MAXTITLESTRINGLENGTH_]; /**< _DELIMITER_ separated string of titles for scalar perturbation output files. */
  char vector_titles_[_MAXTITLESTRINGLENGTH_]; /**< _DELIMITER_ separated string of titles for vector perturbation output files. */
  char tensor_titles_[_MAXTITLESTRINGLENGTH_]; /**< _DELIMITER_ separated string of titles for tensor perturbation output files. */

  double* scalar_perturbations_data_[_MAX_NUMBER_OF_K_FILES_]; /**< Array of double pointers to perturbation output for scalars */
  double* vector_perturbations_data_[_MAX_NUMBER_OF_K_FILES_]; /**< Array of double pointers to perturbation output for vectors */
  double* tensor_perturbations_data_[_MAX_NUMBER_OF_K_FILES_]; /**< Array of double pointers to perturbation output for tensors */

  int size_scalar_perturbation_data_[_MAX_NUMBER_OF_K_FILES_]; /**< Array of sizes of scalar double pointers  */
  int size_vector_perturbation_data_[_MAX_NUMBER_OF_K_FILES_]; /**< Array of sizes of vector double pointers  */
  int size_tensor_perturbation_data_[_MAX_NUMBER_OF_K_FILES_]; /**< Array of sizes of tensor double pointers  */
  //@}

  /** @name - source functions interpolation table */

  //@{

  double*** sources_; /**< Pointer towards the source interpolation table
                         sources[index_md]
                         [index_ic * perturbations_module_->tp_size_[index_md] + index_tp]
                         [index_tau * perturbations_module_->k_size_ + index_k] */

  //@}
  /** @name - arrays related to the interpolation table for sources at late times, corresponding to z < z_max_pk (used for Fourier transfer function and spectra output) */
  //@{

  double* ln_tau_;      /**< log of the arrau tau_sampling, covering only the final time range required for the output of
                            Fourier transfer functions (used for interpolations) */
  int ln_tau_size_;     /**< number of values in this array */

  //@}

  double* tau_sampling_;    /**< array of tau values */
  int tau_size_;            /**< number of values in this array */

  int* k_size_cl_;  /**< k_size_cl[index_md] number of k values used
                       for non-CMB \f$ C_l \f$ calculations, requiring a coarse
                       sampling in k-space. */
  int* k_size_;     /**< k_size[index_md] = total number of k
                       values, including those needed for P(k) but not
                       for \f$ C_l \f$'s */
  double** k_;      /**< k[index_md][index_k] = list of values */
  double k_min_;    /**< minimum value (over all modes) */
  double k_max_;    /**< maximum value (over all modes) */


private:
  int perturb_sources_at_tau(int index_md, int index_ic, int index_tp, double tau, double* pvecsources) const;
  int perturb_init();
  int perturb_free();
  int perturb_indices_of_perturbs();
  int perturb_timesampling_for_sources();
  int perturb_get_k_list();
  int perturb_workspace_init(int index_md, perturb_workspace* ppw);
  int perturb_workspace_free(int index_md, perturb_workspace* ppw);
  int perturb_solve(int index_md, int index_ic, int index_k, perturb_workspace* ppw);
  int perturb_prepare_k_output();
  int perturb_find_approximation_number(int index_md, double k, perturb_workspace* ppw, double tau_ini, double tau_end, int* interval_number, int* interval_number_of);
  int perturb_find_approximation_switches(int index_md, double k, perturb_workspace* ppw, double tau_ini, double tau_end, double precision, int interval_number, int* interval_number_of, double* interval_limit, int** interval_approx);
  int perturb_vector_init(int index_md, int index_ic, double k, double tau, perturb_workspace* ppw, int* pa_old);
  int perturb_vector_free(struct perturb_vector * pv);
  int perturb_initial_conditions(int index_md, int index_ic, double k, double tau, perturb_workspace* ppw);
  int perturb_approximations(int index_md, double k, double tau, perturb_workspace* ppw);
  int perturb_einstein(int index_md, double k, double tau, double* y, perturb_workspace* ppw);
  int perturb_total_stress_energy(int index_md, double k, double* y, perturb_workspace* ppw);
  int perturb_timescale_member(double tau, void * parameters_and_workspace, double* timescale, ErrorMsg error_message);
  static int perturb_timescale(double tau, void * parameters_and_workspace, double* timescale, ErrorMsg error_message);
  int perturb_sources_member(double tau, double* pvecperturbations, double* pvecderivs, int index_tau, void * parameters_and_workspace, ErrorMsg error_message);
  static int perturb_sources(double tau, double* pvecperturbations, double* pvecderivs, int index_tau, void * parameters_and_workspace, ErrorMsg error_message);
  int perturb_print_variables_member(double tau, double* y, double* dy, void * parameters_and_workspace, ErrorMsg error_message);
  static int perturb_print_variables(double tau, double* y, double* dy, void * parameters_and_workspace, ErrorMsg error_message);
  int perturb_derivs_member(double tau, double* y, double* dy, void * parameters_and_workspace, ErrorMsg error_message);
  static int perturb_derivs(double tau, double* y, double* dy, void * parameters_and_workspace, ErrorMsg error_message);
  int perturb_tca_slip_and_shear(double* y, void * parameters_and_workspace, ErrorMsg error_message);
  int perturb_rsa_delta_and_theta(double k, double* y, double a_prime_over_a, double* pvecthermo, perturb_workspace* ppw);
  int perturb_rsa_idr_delta_and_theta(double k, double* y, double a_prime_over_a, double* pvecthermo, perturb_workspace* ppw);

  BackgroundModulePtr background_module_;
  ThermodynamicsModulePtr thermodynamics_module_;

  short evolve_tensor_ur_;             /**< will we evolve ur tensor perturbations (either because we have ur species, or we have ncdm species with massless approximation) ? */
  short evolve_tensor_ncdm_;             /**< will we evolve ncdm tensor perturbations (if we have ncdm species and we use the exact method) ? */

  /** @name - useful flags  */

  //@{

  short has_cmb_; /**< do we need CMB-related sources (temperature, polarization) ? */
  short has_lss_; /**< do we need LSS-related sources (lensing potential, ...) ? */

  //@}


  double*** late_sources_; /**< Pointer towards the source interpolation table
                                late_sources[index_md]
                                            [index_ic * perturbations_module_->tp_size_[index_md] + index_tp]
                                            [index_tau * perturbations_module_->k_size_ + index_k]
                                Note that this is not a replication of part of the sources table,
                                it is just poiting towards the same memory zone, at the place where the late_sources actually start */

  double*** ddlate_sources_; /**< Pointer towards the splined source interpolation table with second derivatives with respect to time
                              ddlate_sources[index_md]
                                            [index_ic * perturbations_module_->tp_size_[index_md] + index_tp]
                                            [index_tau * perturbations_module_->k_size_ + index_k] */

  //@}
  int number_of_scalar_titles_; /**< number of titles/columns in scalar perturbation output files */
  int number_of_vector_titles_; /**< number of titles/columns in vector perturbation output files*/
  int number_of_tensor_titles_; /**< number of titles/columns in tensor perturbation output files*/

  /** @name - list of k values for each mode */

  //@{

  int* k_size_cmb_;  /**< k_size_cmb[index_md] number of k values used
                        for CMB calculations, requiring a fine
                        sampling in k-space */



  //@}
};

/**
 * Structure pointing towards all what the function that perturb_derivs
 * needs to know: fixed input parameters and indices contained in the
 * various structures, workspace, etc.
 */

struct perturb_parameters_and_workspace {
  perturb_parameters_and_workspace(PerturbationsModule* p_m) : perturbations_module(p_m) {}
  PerturbationsModule* const perturbations_module;
  int index_md;                   /**< index of mode (scalar/.../vector/tensor) */
  int index_ic;                /**< index of initial condition (adiabatic/isocurvature(s)/...) */
  int index_k;                /**< index of wavenumber */
  double k;                    /**< current value of wavenumber in 1/Mpc */
  perturb_workspace* ppw; /**< workspace defined above */

};

#endif //PERTURBATIONS_MODULE_H
