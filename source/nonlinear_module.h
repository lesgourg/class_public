#ifndef NONLINEAR_MODULE_H
#define NONLINEAR_MODULE_H

#include "input_module.h"
#include "base_module.h"

class NonlinearModule : public BaseModule {
public:
  NonlinearModule(InputModulePtr input_module, BackgroundModulePtr background_module, PerturbationsModulePtr perturbations_module, PrimordialModulePtr primordial_module);
  ~NonlinearModule();

  /* external functions (meant to be called from other modules) */
  int nonlinear_pk_at_z(enum linear_or_logarithmic mode, enum pk_outputs pk_output, double z, int index_pk, double* out_pk, double* out_pk_ic) const;
  int nonlinear_pks_at_z(enum linear_or_logarithmic mode, enum pk_outputs pk_output, double z, double* out_pk, double* out_pk_ic, double* out_pk_cb, double* out_pk_cb_ic) const;
  int nonlinear_pk_at_k_and_z(enum pk_outputs pk_output, double k, double z, int index_pk, double* out_pk, double* out_pk_ic) const;
  int nonlinear_pks_at_k_and_z(enum pk_outputs pk_output, double k, double z, double* out_pk, double* out_pk_ic, double* out_pk_cb, double* out_pk_cb_ic) const;
  int nonlinear_pks_at_kvec_and_zvec(enum pk_outputs pk_output, double* kvec, int kvec_size, double* zvec, int zvec_size, double* out_pk, double* out_pk_cb) const;
  int nonlinear_sigmas_at_z(double R, double z, int index_pk, enum out_sigmas sigma_output, double* result) const;
  int nonlinear_pk_tilt_at_k_and_z(enum pk_outputs pk_output, double k, double z, int index_pk, double* pk_tilt) const;
  int nonlinear_k_nl_at_z(double z, double* k_nl, double* k_nl_cb) const;

  // Deprecated:
  int nonlinear_sigma_at_z(double R, double z, int index_pk, double k_per_decade, double* result) const;

  int k_size_;      /**< k_size = total number of k values */
  double* ln_k_;    /**< ln_k[index_k] = list of log(k) values */

  double** nl_corr_density_;   /**< nl_corr_density[index_pk][index_tau * perturbations_module_->k_size_ + index_k] */

  short* is_non_zero_;  /**< for a given mode, is_non_zero[index_md][index_ic1_ic2] is set to true if the pair of initial conditions (index_ic1, index_ic2) are statistically correlated, or to false if they are uncorrelated */
  int ic_size_;         /**< for a given mode, ic_size[index_md] = number of initial conditions included in computation */
  int ic_ic_size_;      /**< for a given mode, ic_ic_size[index_md] = number of pairs of (index_ic1, index_ic2) with index_ic2 >= index_ic1; this number is just N(N+1)/2  where N = ic_size[index_md] */
  short has_pk_m_;  /**< do we want spectra for total matter? */
  short has_pk_cb_; /**< do we want spectra for cdm+baryons? */

  int index_pk_m_;  /**< index of pk for matter (defined only when has_pk_m is TRUE) */
  int index_pk_cb_; /**< index of pk for cold dark matter plus baryons (defined only when has_pk_cb is TRUE */
  int pk_size_;     /**< k_size = total number of pk */
  double* sigma8_;   /**< sigma8[index_pk] */


private:
  /* internal functions */
  int nonlinear_init();
  int nonlinear_free();
  int nonlinear_indices();
  int nonlinear_get_k_list();
  int nonlinear_get_tau_list();
  int nonlinear_get_source(int index_k, int index_ic, int index_tp, int index_tau, double** sources, double* source);
  int nonlinear_pk_linear(int index_pk, int index_tau, int k_size, double* lnpk, double* lnpk_ic);
  int nonlinear_sigmas(double R, double* lnpk_l, double* ddlnpk_l, int k_size, double k_per_decade, enum out_sigmas sigma_output, double* result) const;
  int nonlinear_halofit(int index_pk, double tau, double* pk_nl, double* lnpk_l, double* ddlnpk_l, double* k_nl, short* halofit_found_k_max);
  int nonlinear_halofit_integrate(double* integrand_array, int integrand_size, int ia_size, int index_ia_k, int index_ia_pk, int index_ia_sum, int index_ia_ddsum, double R, enum halofit_integral_type type, double* sum);
  int nonlinear_hmcode(int index_pk, int index_tau, double tau, double*pk_nl, double** lnpk_l, double** ddlnpk_l, double* k_nl, short* halofit_found_k_max, nonlinear_workspace* pnw);
  int nonlinear_hmcode_workspace_init(nonlinear_workspace* pnw);
  int nonlinear_hmcode_workspace_free(nonlinear_workspace* pnw);
  int nonlinear_hmcode_dark_energy_correction(nonlinear_workspace* pnw);
  int nonlinear_hmcode_baryonic_feedback();
  int nonlinear_hmcode_fill_sigtab(int index_tau, double*lnpk_l, double*ddlnpk_l, nonlinear_workspace* pnw);
  int nonlinear_hmcode_fill_growtab(nonlinear_workspace* pnw);
  int nonlinear_hmcode_growint(double a, double w, double wa, double* growth);
  int nonlinear_hmcode_window_nfw(double k, double rv, double c, double* window_nfw);
  int nonlinear_hmcode_halomassfunction(double nu, double* hmf);
  int nonlinear_hmcode_sigma8_at_z(double z, double* sigma_8, double* sigma_8_cb, nonlinear_workspace* pnw);
  int nonlinear_hmcode_sigmadisp_at_z(double z, double* sigma_disp, double* sigma_disp_cb, nonlinear_workspace* pnw);
  int nonlinear_hmcode_sigmadisp100_at_z(double z, double* sigma_disp_100, double* sigma_disp_100_cb, nonlinear_workspace* pnw);
  int nonlinear_hmcode_sigmaprime_at_z(double z, double* sigma_prime, double* sigma_prime_cb, nonlinear_workspace* pnw);
  int prepare_pk_eq();

  BackgroundModulePtr background_module_;
  PerturbationsModulePtr perturbations_module_;
  PrimordialModulePtr primordial_module_;

  /** @name - arrays for the Fourier power spectra P(k,tau) */

  //@{

  short has_pk_matter_; /**< do we need matter Fourier spectrum? */

  double* k_;      /**< k[index_k] = list of k values */

  double* ln_tau_;     /**< log(tau) array, only needed if user wants
                        some output at z>0, instead of only z=0.  This
                        array only covers late times, used for the
                        output of P(k) or T(k), and matching the
                        condition z(tau) < z_max_pk */

  int ln_tau_size_;     /**< number of values in this array */

  double** ln_pk_ic_l_;   /**< Matter power spectrum (linear).
                             Depends on indices index_pk, index_ic1_ic2, index_k, index_tau as:
                             ln_pk_ic_l[index_pk][(index_tau * k_size_ + index_k)* ic_ic_size_ + index_ic1_ic2]
                             where index-pk labels P(k) types (_m = total matter, _cb = baryons+CDM),
                             while index_ic1_ic2 labels ordered pairs (index_ic1, index_ic2) (since
                             the primordial spectrum is symmetric in (index_ic1, index_ic2)).
                             - for diagonal elements (index_ic1 = index_ic2) this arrays contains
                             ln[P(k)] where P(k) is positive by construction.
                             - for non-diagonal elements this arrays contains the k-dependent
                             cosine of the correlation angle, namely
                             P(k)_(index_ic1, index_ic2)/sqrt[P(k)_index_ic1 P(k)_index_ic2]
                             This choice is convenient since the sign of the non-diagonal cross-correlation
                             could be negative. For fully correlated or anti-correlated initial conditions,
                             this non-diagonal element is independent on k, and equal to +1 or -1.
                          */

  double** ddln_pk_ic_l_; /**< second derivative of above array with respect to log(tau), for spline interpolation. So:
                             - for index_ic1 = index_ic, we spline ln[P(k)] vs. ln(k), which is
                             good since this function is usually smooth.
                             - for non-diagonal coefficients, we spline
                             P(k)_(index_ic1, index_ic2)/sqrt[P(k)_index_ic1 P(k)_index_ic2]
                             vs. ln(k), which is fine since this quantity is often assumed to be
                             constant (e.g for fully correlated/anticorrelated initial conditions)
                             or nearly constant, and with arbitrary sign.
                          */

  double** ln_pk_l_;   /**< Total matter power spectrum summed over initial conditions (linear).
                          Only depends on indices index_pk,index_k, index_tau as:
                          ln_pk[index_pk][index_tau * k_size_ + index_k]
                       */

  double** ddln_pk_l_; /**< second derivative of above array with respect to log(tau), for spline interpolation. */

  double** ln_pk_nl_;   /**< Total matter power spectrum summed over initial conditions (nonlinear).
                          Only depends on indices index_pk,index_k, index_tau as:
                          ln_pk[index_pk][index_tau * k_size_ + index_k]
                       */

  double** ddln_pk_nl_; /**< second derivative of above array with respect to log(tau), for spline interpolation. */

  //@}

  /** @name - table non-linear corrections for matter density, sqrt(P_NL(k,z)/P_NL(k,z)) */

  //@{

  int k_size_extra_;/** total number of k values of extrapolated k array (high k)*/

  int tau_size_;    /**< tau_size = number of values */
  double* tau_;    /**< tau[index_tau] = list of time values, covering
                      all the values of the perturbation module */

  double** k_nl_;              /**< wavenumber at which non-linear corrections become important,
                                    defined differently by different non_linear_method's */
  int index_tau_min_nl_;        /**< index of smallest value of tau at which nonlinear corrections have been computed
                                    (so, for tau<tau_min_nl, the array nl_corr_density only contains some factors 1 */

  //@}

  /** @name - information on number of modes and pairs of initial conditions */

  //@{

  int index_md_scalars_; /**< set equal to psp->index_md_scalars
                           (useful since this module only deals with
                           scalars) */

  //@}

  /** @name - information on the type of power spectra (_cb, _m...) */

  //@{


  /* and two redundent but useful indices: */

  int index_pk_total_;      /**< always equal to index_pk_m
                              (always defined, useful e.g. for weak lensing spectrum) */
  int index_pk_cluster_;    /**< equal to index_pk_cb if it exists, otherwise to index_pk_m
                              (always defined, useful e.g. for galaxy clustering spectrum) */

  double c_min_;      /** for HMcode: minimum concentration in Bullock 2001 mass-concentration relation */
  double eta_0_;      /** for HMcode: halo bloating parameter */

  //@}

  /** @name - parameters for the pk_eq method */

  //@{

  int index_pk_eq_w_;                /**< index of w in table pk_eq_w_and_Omega */
  int index_pk_eq_Omega_m_;          /**< index of Omega_m in table pk_eq_w_and_Omega */
  int pk_eq_size_;                   /**< number of indices in table pk_eq_w_and_Omega */
  int pk_eq_tau_size_;               /**< number of times (and raws in table pk_eq_w_and_Omega) */
  double* pk_eq_tau_;                /**< table of time values */
  double* pk_eq_w_and_Omega_;        /**< table of background quantites */
  double* pk_eq_ddw_and_ddOmega_;    /**< table of second derivatives */

  //@}

};

#endif //NONLINEAR_MODULE_H
