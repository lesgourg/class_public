#ifndef BACKGROUND_MODULE_H
#define BACKGROUND_MODULE_H

#include "input_module.h"
#include "base_module.h"

class BackgroundModule : public BaseModule {
public:
  BackgroundModule(InputModulePtr input_module);
  ~BackgroundModule();
  int background_output_titles(char titles[_MAXTITLESTRINGLENGTH_]) const;
  int background_output_data(int number_of_titles, double* data) const;
  int background_at_tau(double tau, short return_format, short inter_mode, int* last_index, double* pvecback) const;
  int background_tau_of_z(double z, double* tau) const;
  int background_w_fld(double a, double* w_fld, double* dw_over_da_fld, double* integral_fld) const;
  int background_free_noinput() const;
  double dV_scf(double phi) const;

  /** @name - all indices for the vector of background (=bg) quantities stored in table */

  //@{

  int index_bg_a_;             /**< scale factor */
  int index_bg_H_;             /**< Hubble parameter in \f$Mpc^{-1}\f$ */
  int index_bg_H_prime_;       /**< its derivative w.r.t. conformal time */

  /* end of vector in short format, now quantities in normal format */

  int index_bg_rho_g_;         /**< photon density */
  int index_bg_rho_b_;         /**< baryon density */
  int index_bg_rho_cdm_;       /**< cdm density */
  int index_bg_rho_lambda_;    /**< cosmological constant density */
  int index_bg_rho_fld_;       /**< fluid density */
  int index_bg_w_fld_;         /**< fluid equation of state */
  int index_bg_rho_ur_;        /**< relativistic neutrinos/relics density */
  int index_bg_rho_idm_dr_;    /**< density of dark matter interacting with dark radiation */
  int index_bg_rho_idr_;       /**< density of interacting dark radiation */
  int index_bg_rho_dcdm_;      /**< dcdm density */
  int index_bg_rho_dr_;        /**< dr density */

  int index_bg_phi_scf_;       /**< scalar field value */
  int index_bg_phi_prime_scf_; /**< scalar field derivative wrt conformal time */
  int index_bg_V_scf_;         /**< scalar field potential V */
  int index_bg_dV_scf_;        /**< scalar field potential derivative V' */
  int index_bg_ddV_scf_;       /**< scalar field potential second derivative V'' */
  int index_bg_rho_scf_;       /**< scalar field energy density */
  int index_bg_p_scf_;         /**< scalar field pressure */
  int index_bg_p_prime_scf_;         /**< scalar field pressure */

  int index_bg_rho_ncdm1_;     /**< density of first ncdm species (others contiguous) */
  int index_bg_p_ncdm1_;       /**< pressure of first ncdm species (others contiguous) */
  int index_bg_pseudo_p_ncdm1_;/**< another statistical momentum useful in ncdma approximation */

  int index_bg_rho_tot_;       /**< Total density */
  int index_bg_p_tot_;         /**< Total pressure */
  int index_bg_p_tot_prime_;   /**< Conf. time derivative of total pressure */

  int index_bg_Omega_r_;       /**< relativistic density fraction (\f$ \Omega_{\gamma} + \Omega_{\nu r} \f$) */

  /* end of vector in normal format, now quantities in long format */

  int index_bg_rho_crit_;      /**< critical density */
  int index_bg_Omega_m_;       /**< non-relativistic density fraction (\f$ \Omega_b + \Omega_cdm + \Omega_{\nu nr} \f$) */
  int index_bg_conf_distance_; /**< conformal distance (from us) in Mpc */
  int index_bg_ang_distance_;  /**< angular diameter distance in Mpc */
  int index_bg_lum_distance_;  /**< luminosity distance in Mpc */
  int index_bg_time_;          /**< proper (cosmological) time in Mpc */
  int index_bg_rs_;            /**< comoving sound horizon in Mpc */

  int index_bg_D_;             /**< scale independent growth factor D(a) for CDM perturbations */
  int index_bg_f_;             /**< corresponding velocity growth factor [dlnD]/[dln a] */

  int bg_size_short_;  /**< size of background vector in the "short format" */
  int bg_size_normal_; /**< size of background vector in the "normal format" */
  int bg_size_;        /**< size of background vector in the "long format" */

  //@}
  int bt_size_;               /**< number of lines (i.e. time-steps) in the array */
  double* tau_table_;        /**< vector tau_table[index_tau] with values of \f$ \tau \f$ (conformal time) */
  double* background_table_; /**< table background_table[index_tau*bg_size_+pba->index_bg] with all other quantities (array of size bg_size*bt_size) **/

  double conformal_age_; /**< conformal age in Mpc */
  double Neff_; /**< so-called "effective neutrino number", computed at earliest time in interpolation table */
  double a_eq_;      /**< scale factor at radiation/matter equality */
  double H_eq_;      /**< Hubble rate at radiation/matter equality [Mpc^-1] */
  double Omega0_m_;  /**< total non-relativistic matter today */
  double Omega0_de_; /**< total dark energy density today, currently defined as 1 - Omega0_m - Omega0_r - Omega0_k */

  double age_; /**< age in Gyears */
  double Omega0_r_;  /**< total ultra-relativistic radiation today */
  double z_eq_;      /**< redshift at radiation/matter equality */
  double tau_eq_;    /**< conformal time at radiation/matter equality [Mpc] */
  double Omega0_dcdm_; /**< \f$ \Omega_{0 dcdm} \f$: decaying cold dark matter */
  double Omega0_dr_; /**< \f$ \Omega_{0 dr} \f$: decay radiation */

private:
  int background_functions(double* pvecback_B, short return_format, double* pvecback);
  int background_init();
  int background_free();
  int background_indices();
  int background_solve();
  int background_solve_evolver();
  int background_initial_conditions(double* pvecback, double* pvecback_integration);
  int background_find_equality();
  int background_derivs_member(double z, double* y, double* dy, void* parameters_and_workspace, ErrorMsg error_message);
  static int background_derivs(double z, double* y, double* dy, void* parameters_and_workspace, ErrorMsg error_message);
  int background_derivs_loga_member(double loga, double* y, double* dy, void* parameters_and_workspace, ErrorMsg error_message);
  static int background_derivs_loga(double loga, double* y, double* dy, void* parameters_and_workspace, ErrorMsg error_message);
  int background_add_line_to_bg_table_member(double loga, double* y, double* dy, int index_loga, void* parameters_and_workspace, ErrorMsg error_message);
  static int background_add_line_to_bg_table(double loga, double* y, double* dy, int index_loga, void* parameters_and_workspace, ErrorMsg error_message);
  double V_scf(double phi) const;
  double ddV_scf(double phi) const;
  double Q_scf(double phi, double phi_prime);
  double V_e_scf(double phi) const;
  double dV_e_scf(double phi) const;
  double ddV_e_scf(double phi) const;
  double V_p_scf(double phi) const;
  double dV_p_scf(double phi) const;
  double ddV_p_scf(double phi) const;
  int background_output_budget();

  /** @name - all indices for the vector of background quantities to be integrated (=bi)
   *
   * Most background quantities can be immediately inferred from the
   * scale factor. Only few of them require an integration with
   * respect to conformal time (in the minimal case, only one quantity needs to
   * be integrated with time: the scale factor, using the Friedmann
   * equation). These indices refer to the vector of
   * quantities to be integrated with time.
   * {B} quantities are needed by background_functions() while {C} quantities are not.
   */

  //@{

  int index_bi_a_;       /**< {B} scale factor */
  int index_bi_rho_dcdm_;/**< {B} dcdm density */
  int index_bi_rho_dr_;  /**< {B} dr density */
  int index_bi_rho_fld_; /**< {B} fluid density */
  int index_bi_phi_scf_;       /**< {B} scalar field value */
  int index_bi_phi_prime_scf_; /**< {B} scalar field derivative wrt conformal time */

  int index_bi_time_;    /**< {C} proper (cosmological) time in Mpc */
  int index_bi_rs_;      /**< {C} sound horizon */
  int index_bi_tau_;     /**< {C} conformal time in Mpc */
  int index_bi_D_;       /**< {C} scale independent growth factor D(a) for CDM perturbations. */
  int index_bi_D_prime_; /**< {C} D satisfies \f$ [D''(\tau)=-aHD'(\tau)+3/2 a^2 \rho_M D(\tau) \f$ */

  int bi_B_size_;        /**< Number of {B} parameters */
  int bi_size_;          /**< Number of {B}+{C} parameters */

  //@}



  /** @name - background interpolation tables */

  //@{

  double* z_table_;          /**< vector z_table[index_tau] with values of \f$ z \f$ (redshift) */

  //@}

  /** @name - table of their second derivatives, used for spline interpolation */

  //@{

  double* d2tau_dz2_table_; /**< vector d2tau_dz2_table[index_tau] with values of \f$ d^2 \tau / dz^2 \f$ (conformal time) */
  double* d2background_dtau2_table_; /**< table d2background_dtau2_table[index_tau*bg_size_+pba->index_bg] with values of \f$ d^2 b_i / d\tau^2 \f$ (conformal time) */

  //@}
};

/**
 * temporary parameters and workspace passed to the background_derivs function
 */

struct background_parameters_and_workspace {
  background_parameters_and_workspace(BackgroundModule* b_m) : background_module(b_m) {}
  BackgroundModule* const background_module;
  double* pvecback;
};


#endif //BACKGROUND_MODULE_H
