/** @file fourier.h Documented includes for trg module */

#include "primordial.h"
#include "trigonometric_integrals.h"

#ifndef __FOURIER__
#define __FOURIER__

#define _M_EV_TOO_BIG_FOR_HALOFIT_ 10. /**< above which value of non-CDM mass (in eV) do we stop trusting halofit? */

#define _M_SUN_ 1.98847e30 /**< Solar mass in Kg */

#define _MAX_NUM_EXTRAPOLATION_ 100000

enum non_linear_method {nl_none,nl_halofit,nl_HMcode};
enum pk_outputs {pk_linear,pk_nonlinear};

enum source_extrapolation {extrap_zero,extrap_only_max,extrap_only_max_units,extrap_max_scaled,extrap_hmcode,extrap_user_defined};

enum halofit_integral_type {halofit_integral_one, halofit_integral_two, halofit_integral_three};

enum hmcode_baryonic_feedback_model {nl_emu_dmonly, nl_owls_dmonly, nl_owls_ref, nl_owls_agn, nl_owls_dblim, nl_user_defined};
enum out_sigmas {out_sigma,out_sigma_prime,out_sigma_disp};

/**
 * Structure containing all information on non-linear spectra.
 *
 * Once initialized by fourier_init(), contains a table for all two points correlation functions
 * and for all the ai,bj functions (containing the three points correlation functions), for each
 * time and wave-number.
 */

struct fourier {

  /** @name - input parameters initialized by user in input module
      (all other quantities are computed in this module, given these
      parameters and the content of the 'precision', 'background',
      'thermo', 'primordial' and 'spectra' structures) */

  //@{

  enum non_linear_method method; /**< method for computing non-linear corrections (none, Halogit, etc.) */

  enum source_extrapolation extrapolation_method; /**< method for analytical extrapolation of sources beyond pre-computed range */

  enum hmcode_baryonic_feedback_model feedback; /** to choose between different baryonic feedback models
                                                in hmcode (dmonly, gas cooling, Agn or supernova feedback) */
  double c_min;      /** for HMcode: minimum concentration in Bullock 2001 mass-concentration relation */
  double eta_0;      /** for HMcode: halo bloating parameter */
  double z_infinity; /** for HMcode: z value at which Dark Energy correction is evaluated needs to be at early times (default */

  short has_pk_eq;  /**< flag: in case wa_fld is defined and non-zero, should we use the pk_eq method? */

  //@}

  /** @name - information on number of modes and pairs of initial conditions */

  //@{

  int index_md_scalars; /**< set equal to phr->index_md_scalars
                           (useful since this module only deals with
                           scalars) */
  int ic_size;         /**< for a given mode, ic_size[index_md] = number of initial conditions included in computation */
  int ic_ic_size;      /**< for a given mode, ic_ic_size[index_md] = number of pairs of (index_ic1, index_ic2) with index_ic2 >= index_ic1; this number is just N(N+1)/2  where N = ic_size[index_md] */
  short * is_non_zero; /**< for a given mode, is_non_zero[index_md][index_ic1_ic2] is set to true if the pair of initial conditions (index_ic1, index_ic2) are statistically correlated, or to false if they are uncorrelated */

  //@}

  /** @name - information on the type of power spectra (_cb, _m...) */

  //@{

  short has_pk_m;  /**< do we want spectra for total matter? */
  short has_pk_cb; /**< do we want spectra for cdm+baryons? */

  int index_pk_m;  /**< index of pk for matter (defined only when has_pk_m is TRUE) */
  int index_pk_cb; /**< index of pk for cold dark matter plus baryons (defined only when has_pk_cb is TRUE */

  /* and two redundent but useful indices: */

  int index_pk_total;      /**< always equal to index_pk_m
                              (always defined, useful e.g. for weak lensing spectrum) */
  int index_pk_cluster;    /**< equal to index_pk_cb if it exists, otherwise to index_pk_m
                              (always defined, useful e.g. for galaxy clustering spectrum) */

  int pk_size;     /**< k_size = total number of pk */

  //@}

  /** @name - arrays for the Fourier power spectra P(k,tau) */

  //@{

  short has_pk_matter; /**< do we need matter Fourier spectrum? */

  int k_size;      /**< k_size = total number of k values */
  double * k;      /**< k[index_k] = list of k values */
  double * ln_k;   /**< ln_k[index_k] = list of log(k) values */

  double * ln_tau;     /**< log(tau) array, only needed if user wants
                        some output at z>0, instead of only z=0.  This
                        array only covers late times, used for the
                        output of P(k) or T(k), and matching the
                        condition z(tau) < z_max_pk */

  int ln_tau_size;     /**< number of values in this array */

  double ** ln_pk_ic_l;   /**< Matter power spectrum (linear).
                             Depends on indices index_pk, index_ic1_ic2, index_k, index_tau as:
                             ln_pk_ic_l[index_pk][(index_tau * pfo->k_size + index_k)* pfo->ic_ic_size + index_ic1_ic2]
                             where index-pk labels P(k) types (m = total matter, cb = baryons+CDM),
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

  double ** ddln_pk_ic_l; /**< second derivative of above array with respect to log(tau), for spline interpolation. So:
                             - for index_ic1 = index_ic, we spline ln[P(k)] vs. ln(k), which is
                             good since this function is usually smooth.
                             - for non-diagonal coefficients, we spline
                             P(k)_(index_ic1, index_ic2)/sqrt[P(k)_index_ic1 P(k)_index_ic2]
                             vs. ln(k), which is fine since this quantity is often assumed to be
                             constant (e.g for fully correlated/anticorrelated initial conditions)
                             or nearly constant, and with arbitrary sign.
                          */

  double ** ln_pk_l;   /**< Total matter power spectrum summed over initial conditions (linear).
                          Only depends on indices index_pk,index_k, index_tau as:
                          ln_pk[index_pk][index_tau * pfo->k_size + index_k]
                       */

  double ** ddln_pk_l; /**< second derivative of above array with respect to log(tau), for spline interpolation. */

  double ** ln_pk_nl;   /**< Total matter power spectrum summed over initial conditions (nonlinear).
                          Only depends on indices index_pk,index_k, index_tau as:
                          ln_pk[index_pk][index_tau * pfo->k_size + index_k]
                       */

  double ** ddln_pk_nl; /**< second derivative of above array with respect to log(tau), for spline interpolation. */

  double * sigma8;   /**< sigma8[index_pk] */

  //@}

  /** @name - table non-linear corrections for matter density, sqrt(P_NL(k,z)/P_NL(k,z)) */

  //@{

  int k_size_extra;/** total number of k values of extrapolated k array (high k)*/

  int tau_size;    /**< tau_size = number of values */
  double * tau;    /**< tau[index_tau] = list of time values, covering
                      all the values of the perturbation module */

  double ** nl_corr_density;   /**< nl_corr_density[index_pk][index_tau * ppt->k_size + index_k] */
  double ** k_nl;              /**< wavenumber at which non-linear corrections become important,
                                    defined differently by different non_linear_method's */
  int index_tau_min_nl;        /**< index of smallest value of tau at which nonlinear corrections have been computed
                                    (so, for tau<tau_min_nl, the array nl_corr_density only contains some factors 1 */

  //@}

  /** @name - parameters for the pk_eq method */

  //@{

  int index_pk_eq_w;                /**< index of w in table pk_eq_w_and_Omega */
  int index_pk_eq_Omega_m;          /**< index of Omega_m in table pk_eq_w_and_Omega */
  int pk_eq_size;                   /**< number of indices in table pk_eq_w_and_Omega */

  int pk_eq_tau_size;               /**< number of times (and raws in table pk_eq_w_and_Omega) */

  double * pk_eq_tau;               /**< table of time values */
  double * pk_eq_w_and_Omega;       /**< table of background quantites */
  double * pk_eq_ddw_and_ddOmega;   /**< table of second derivatives */

  //@}

  /** @name - technical parameters */

  //@{

  short fourier_verbose;  	/**< amount of information written in standard output */

  ErrorMsg error_message; 	/**< zone for writing error messages */

  //@}
};

/**
 * Structure containing variables used only internally in fourier module by various functions.
 *
 */

struct fourier_workspace {

  /** @name - quantitites used by HMcode */

  //@{

  double * rtab; /** List of R values */
  double * stab; /** List of Sigma Values */
  double * ddstab; /** Splined sigma */

  double * growtable;
  double * ztable;
  double * tautable;

  double ** sigma_8;
  double ** sigma_disp;
  double ** sigma_disp_100;
  double ** sigma_prime;

  double dark_energy_correction; /** this is the ratio [g_wcdm(z_infinity)/g_lcdm(z_infinity)]^1.5
                                  * (power comes from Dolag et al. (2004) correction)
                                  * it is 1, if has_fld == _FALSE_ */

  //@}

};

/********************************************************************************/

/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  /* external functions (meant to be called from other modules) */

  int fourier_pk_at_z(
                        struct background * pba,
                        struct fourier *pfo,
                        enum linear_or_logarithmic mode,
                        enum pk_outputs pk_output,
                        double z,
                        int index_pk,
                        double * out_pk,
                        double * out_pk_ic
                        );

  int fourier_pks_at_z(
                         struct background * pba,
                         struct fourier *pfo,
                         enum linear_or_logarithmic mode,
                         enum pk_outputs pk_output,
                         double z,
                         double * out_pk,
                         double * out_pk_ic,
                         double * out_pk_cb,
                         double * out_pk_cb_ic
                         );

  int fourier_pk_at_k_and_z(
                              struct background * pba,
                              struct primordial * ppm,
                              struct fourier *pfo,
                              enum pk_outputs pk_output,
                              double k,
                              double z,
                              int index_pk,
                              double * out_pk,
                              double * out_pk_ic
                              );

  int fourier_pks_at_k_and_z(
                               struct background * pba,
                               struct primordial * ppm,
                               struct fourier *pfo,
                               enum pk_outputs pk_output,
                               double k,
                               double z,
                               double * out_pk,
                               double * out_pk_ic,
                               double * out_pk_cb,
                               double * out_pk_cb_ic
                               );

  int fourier_pks_at_kvec_and_zvec(
                                     struct background * pba,
                                     struct fourier * pfo,
                                     enum pk_outputs pk_output,
                                     double * kvec,
                                     int kvec_size,
                                     double * zvec,
                                     int zvec_size,
                                     double * out_pk,
                                     double * out_pk_cb
                                     );

  int fourier_sigmas_at_z(
                            struct precision * ppr,
                            struct background * pba,
                            struct fourier * pfo,
                            double R,
                            double z,
                            int index_pk,
                            enum out_sigmas sigma_output,
                            double * result
                            );

  int fourier_pk_tilt_at_k_and_z(
                                    struct background * pba,
                                    struct primordial * ppm,
                                    struct fourier * pfo,
                                    enum pk_outputs pk_output,
                                    double k,
                                    double z,
                                    int index_pk,
                                    double * pk_tilt
                                    );

  int fourier_k_nl_at_z(
                          struct background *pba,
                          struct fourier * pfo,
                          double z,
                          double * k_nl,
                          double * k_nl_cb
                          );

  /* internal functions */

  int fourier_init(
                     struct precision *ppr,
                     struct background *pba,
                     struct thermodynamics *pth,
                     struct perturbations *ppt,
                     struct primordial *ppm,
                     struct fourier *pfo
                     );

  int fourier_free(
                     struct fourier *pfo
                     );

  int fourier_indices(
                        struct precision *ppr,
                        struct background *pba,
                        struct perturbations * ppt,
                        struct primordial * ppm,
                        struct fourier * pfo
                        );

  int fourier_get_k_list(
                           struct precision *ppr,
                           struct perturbations * ppt,
                           struct fourier * pfo
                           );

  int fourier_get_tau_list(
                             struct perturbations * ppt,
                             struct fourier * pfo
                             );

  int fourier_get_source(
                           struct background * pba,
                           struct perturbations * ppt,
                           struct fourier * pfo,
                           int index_k,
                           int index_ic,
                           int index_tp,
                           int index_tau,
                           double ** sources,
                           double * source);

  int fourier_pk_linear(
                          struct background *pba,
                          struct perturbations *ppt,
                          struct primordial *ppm,
                          struct fourier *pfo,
                          int index_pk,
                          int index_tau,
                          int k_size,
                          double * lnpk,
                          double * lnpk_ic
                          );

  int fourier_sigmas(
                       struct fourier * pfo,
                       double R,
                       double *lnpk_l,
                       double *ddlnpk_l,
                       int k_size,
                       double k_per_decade,
                       enum out_sigmas sigma_output,
                       double * result
                       );

  int fourier_sigma_at_z(
                           struct background * pba,
                           struct fourier * pfo,
                           double R,
                           double z,
                           int index_pk,
                           double k_per_decade,
                           double * result
                           );

  int fourier_halofit(
                        struct precision *ppr,
                        struct background *pba,
                        struct perturbations *ppt,
                        struct primordial *ppm,
                        struct fourier *pfo,
                        int index_pk,
                        double tau,
                        double *pk_nl,
                        double *lnpk_l,
                        double *ddlnpk_l,
                        double *k_nl,
                        short * halofit_found_k_max
                        );

  int fourier_halofit_integrate(
                                  struct fourier *pfo,
                                  double * integrand_array,
                                  int integrand_size,
                                  int ia_size,
                                  int index_ia_k,
                                  int index_ia_pk,
                                  int index_ia_sum,
                                  int index_ia_ddsum,
                                  double R,
                                  enum halofit_integral_type type,
                                  double * sum
                                  );

  int fourier_hmcode(
                       struct precision *ppr,
                       struct background *pba,
                       struct perturbations *ppt,
                       struct primordial *ppm,
                       struct fourier *pfo,
                       int index_pk,
                       int index_tau,
                       double tau,
                       double *pk_nl,
                       double **lnpk_l,
                       double **ddlnpk_l,
                       double *k_nl,
                       short * halofit_found_k_max,
                       struct fourier_workspace * pnw
                       );

  int fourier_hmcode_workspace_init(
                                      struct precision *ppr,
                                      struct background *pba,
                                      struct fourier *pfo,
                                      struct fourier_workspace * pnw
                                      );

  int fourier_hmcode_workspace_free(
                                      struct fourier *pfo,
                                      struct fourier_workspace * pnw
                                      );

  int fourier_hmcode_dark_energy_correction(
                                              struct precision *ppr,
                                              struct background *pba,
                                              struct fourier *pfo,
                                              struct fourier_workspace * pnw
                                              );

  int fourier_hmcode_baryonic_feedback(
                                         struct fourier *pfo
                                         );


  int fourier_hmcode_fill_sigtab(
                                   struct precision *ppr,
                                   struct background * pba,
                                   struct perturbations *ppt,
                                   struct primordial * ppm,
                                   struct fourier * pfo,
                                   int index_tau,
                                   double *lnpk_l,
                                   double *ddlnpk_l,
                                   struct fourier_workspace * pnw
                                   );

  int fourier_hmcode_fill_growtab(
                                    struct precision *ppr,
                                    struct background * pba,
                                    struct fourier * pfo,
                                    struct fourier_workspace * pnw
                                    );

  int fourier_hmcode_growint(
                               struct precision *ppr,
                               struct background * pba,
                               struct fourier * pfo,
                               double a,
                               double w,
                               double wa,
                               double * growth
                               );

  int fourier_hmcode_window_nfw(
                                  struct fourier * pfo,
                                  double k,
                                  double rv,
                                  double c,
                                  double *window_nfw
                                  );

  int fourier_hmcode_halomassfunction(
                                        double nu,
                                        double *hmf
                                        );

  int fourier_hmcode_sigma8_at_z(
                        struct background *pba,
                        struct fourier * pfo,
                        double z,
                        double * sigma_8,
                        double * sigma_8_cb,
                        struct fourier_workspace * pnw
                        );

  int fourier_hmcode_sigmadisp_at_z(
                        struct background *pba,
                        struct fourier * pfo,
                        double z,
                        double * sigma_disp,
                        double * sigma_disp_cb,
                        struct fourier_workspace * pnw
                        );

  int fourier_hmcode_sigmadisp100_at_z(
                        struct background *pba,
                        struct fourier * pfo,
                        double z,
                        double * sigma_disp_100,
                        double * sigma_disp_100_cb,
                        struct fourier_workspace * pnw
                        );

  int fourier_hmcode_sigmaprime_at_z(
                        struct background *pba,
                        struct fourier * pfo,
                        double z,
                        double * sigma_prime,
                        double * sigma_prime_cb,
                        struct fourier_workspace * pnw
                        );

#ifdef __cplusplus
}
#endif

#endif
/* @endcond */
