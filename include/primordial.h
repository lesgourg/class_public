/** @file primordial.h Documented includes for primordial module. */

#ifndef __PRIMORDIAL__
#define __PRIMORDIAL__

#include "perturbations.h"

/** enum defining how should the primnordial spectrum be compurted */

enum primordial_spectrum_type {
  analytic_Pk,
  two_scales,
  inflation_V,
  inflation_H,
  external_Pk
};

/** enum defining whether the spectrum routine works with linear or logarithmic input/output */

enum linear_or_logarithmic {
  linear,
  logarithmic
};

/** enum defining the type of inflation potential function V(phi) */

enum potential_shape {
  polynomial,
  natural
};

/**
 * Structure containing everything about primordial spectra that other modules need to know.
 *
 * Once initialized by primordial_init(), contains a table of all
 * primordial spectra as a function of wavenumber, mode, and pair of initial conditions.
 */

struct primordial {

  /** @name - input parameters initialized by user in input module
      (all other quantitites are computed in this module, given these parameters
      and the content of the 'precision' and 'perturbs' structures) */

  //@{

  double k_pivot; /**< pivot scale in Mpc-1 */

  enum primordial_spectrum_type primordial_spec_type; /**< type of primordial spectrum (simple analytic from, integration of inflationary perturbations, etc.) */

  /** - parameters describing the case primordial_spec_type = analytic_Pk : amplitudes, tilts, runnings, cross-correlations, ... */

  double A_s;  /**< usual scalar amplitude = curvature power spectrum at pivot scale */
  double n_s;  /**< usual scalar tilt = [curvature power spectrum tilt at pivot scale -1] */
  double alpha_s; /**< usual scalar running */
  double beta_s;  /**< running of running */

  double r;    /**< usual tensor to scalar ratio of power spectra, r=A_T/A_S=P_h/P_R */
  double n_t;  /**< usual tensor tilt = [GW power spectrum tilt at pivot scale] */
  double alpha_t; /**< usual tensor running */

  double f_bi;  /**< baryon isocurvature (BI) entropy-to-curvature ratio S_bi/R */
  double n_bi;  /**< BI tilt */
  double alpha_bi; /**< BI running */

  double f_cdi;  /**< CDM isocurvature (CDI) entropy-to-curvature ratio S_cdi/R */
  double n_cdi;  /**< CDI tilt */
  double alpha_cdi; /**< CDI running */

  double f_nid;  /**< neutrino density isocurvature (NID) entropy-to-curvature ratio S_nid/R */
  double n_nid;  /**< NID tilt */
  double alpha_nid; /**< NID running */

  double f_niv;  /**< neutrino velocity isocurvature (NIV) entropy-to-curvature ratio S_niv/R */
  double n_niv;  /**< NIV tilt */
  double alpha_niv; /**< NIV running */

  double c_ad_bi; /**< ADxBI cross-correlation at pivot scale, from -1 to 1 */
  double n_ad_bi; /**< ADxBI cross-correlation tilt */
  double alpha_ad_bi; /**< ADxBI cross-correlation running */

  double c_ad_cdi; /**< ADxCDI cross-correlation at pivot scale, from -1 to 1 */
  double n_ad_cdi; /**< ADxCDI cross-correlation tilt */
  double alpha_ad_cdi; /**< ADxCDI cross-correlation running */

  double c_ad_nid; /**< ADxNID cross-correlation at pivot scale, from -1 to 1 */
  double n_ad_nid; /**< ADxNID cross-correlation tilt */
  double alpha_ad_nid; /**< ADxNID cross-correlation running */

  double c_ad_niv; /**< ADxNIV cross-correlation at pivot scale, from -1 to 1 */
  double n_ad_niv; /**< ADxNIV cross-correlation tilt */
  double alpha_ad_niv; /**< ADxNIV cross-correlation running */

  double c_bi_cdi; /**< BIxCDI cross-correlation at pivot scale, from -1 to 1 */
  double n_bi_cdi; /**< BIxCDI cross-correlation tilt */
  double alpha_bi_cdi; /**< BIxCDI cross-correlation running */

  double c_bi_nid; /**< BIxNIV cross-correlation at pivot scale, from -1 to 1 */
  double n_bi_nid; /**< BIxNIV cross-correlation tilt */
  double alpha_bi_nid; /**< BIxNIV cross-correlation running */

  double c_bi_niv; /**< BIxNIV cross-correlation at pivot scale, from -1 to 1 */
  double n_bi_niv; /**< BIxNIV cross-correlation tilt */
  double alpha_bi_niv; /**< BIxNIV cross-correlation running */

  double c_cdi_nid; /**< CDIxNID cross-correlation at pivot scale, from -1 to 1 */
  double n_cdi_nid; /**< CDIxNID cross-correlation tilt */
  double alpha_cdi_nid; /**< CDIxNID cross-correlation running */

  double c_cdi_niv; /**< CDIxNIV cross-correlation at pivot scale, from -1 to 1 */
  double n_cdi_niv; /**< CDIxNIV cross-correlation tilt */
  double alpha_cdi_niv; /**< CDIxNIV cross-correlation running */

  double c_nid_niv; /**< NIDxNIV cross-correlation at pivot scale, from -1 to 1 */
  double n_nid_niv; /**< NIDxNIV cross-correlation tilt */
  double alpha_nid_niv; /**< NIDxNIV cross-correlation running */

  /** - parameters describing the case primordial_spec_type = inflation_V */

  enum potential_shape potential;

  double phi_pivot;
  double V0;
  double V1;
  double V2;
  double V3;
  double V4;

  /** - 'external_Pk' mode: command generating the table of Pk and custom parameters to be passed to it */

  char*  command;
  double custom1;
  double custom2;
  double custom3;
  double custom4;
  double custom5;
  double custom6;
  double custom7;
  double custom8;
  double custom9;
  double custom10;

  //@}

  /** @name - pre-computed table of primordial spectra, and related quantities */

  //@{

  int md_size;      /**< number of modes included in computation */

  int * ic_size;    /**< for a given mode, ic_size[index_md] = number of initial conditions included in computation */

  int * ic_ic_size; /**< number of ordered pairs of (index_ic1, index_ic2); this number is just N(N+1)/2  where N = ic_size[index_md] */

  int lnk_size;    /**< number of ln(k) values */

  double * lnk;    /**< list of ln(k) values lnk[index_k] */

  double ** lnpk;  /**< depends on indices index_md, index_ic1, index_ic2, index_k as:
                      lnpk[index_md][index_k*ppm->ic_ic_size[index_md]+index_ic1_ic2]
                      where index_ic1_ic2 labels ordered pairs (index_ic1, index_ic2) (since
                      the primordial spectrum is symmetric in (index_ic1, index_ic2)).
                      - for diagonal elements (index_ic1 = index_ic2) this arrays contains
                      ln[P(k)] where P(k) is positive by construction.
                      - for non-diagonal elements this arrays contains the k-dependent
                      cosine of the correlation angle, namely
                      P(k )_(index_ic1, index_ic2)/sqrt[P(k)_index_ic1 P(k)_index_ic2]
                      This choice is convenient since the sign of the non-diagonal cross-correlation
                      is arbitrary. For fully correlated or anti-correlated initial conditions,
                      this non -diagonal element is independent on k, and equal to +1 or -1.
                   */

  double ** ddlnpk; /**< second derivative of above array, for spline interpolation. So:
                       - for index_ic1 = index_ic, we spline ln[P(k)] vs. ln(k), which is
                       good since this function is usually smooth.
                       - for non-diagonal coefficients, we spline
                       P(k)_(index_ic1, index_ic2)/sqrt[P(k)_index_ic1 P(k)_index_ic2]
                       vs. ln(k), which is fine since this quantity is often assumed to be
                       constant (e.g for fully correlated/anticorrelated initial conditions)
                       or nearly constant, and with arbitrary sign.
                    */

  short ** is_non_zero; /**< is_non_zero[index_md][index_ic1_ic2] set to false if pair
                           (index_ic1, index_ic2) is uncorrelated
                           (ensures more precision and saves time with respect to the option
                           of simply setting P(k)_(index_ic1, index_ic2) to zero) */

  //@}

  //@{

  /** @name - parameters describing the case primordial_spec_type = analytic_Pk : amplitudes, tilts, runnings, cross-correlations, ... */

  double ** amplitude; /**< all amplitudes in matrix form: amplitude[index_md][index_ic1_ic2] */
  double ** tilt;      /**< all tilts in matrix form: tilt[index_md][index_ic1_ic2] */
  double ** running;   /**< all runnings in matrix form: running[index_md][index_ic1_ic2] */

  //@}

  //@{

  /** @name - for the inflation simulator, indices in vector of
      background/perturbation */

  int index_in_a;       /**< scale factor */
  int index_in_phi;     /**< inflaton vev */
  int index_in_dphi;    /**< its time derivative */
  int index_in_ksi_re;  /**< Mukhanov variable (real part) */
  int index_in_ksi_im;  /**< Mukhanov variable (imaginary part) */
  int index_in_dksi_re; /**< Mukhanov variable (real part, time derivative) */
  int index_in_dksi_im; /**< Mukhanov variable (imaginary part, time derivative) */
  int index_in_ah_re;   /**< tensor perturbation (real part) */
  int index_in_ah_im;   /**< tensor perturbation (imaginary part) */
  int index_in_dah_re;  /**< tensor perturbation (real part, time derivative) */
  int index_in_dah_im;  /**< tensor perturbation (imaginary part, time derivative) */
  int in_bg_size;       /**< size of vector of background quantitites only */
  int in_size;          /**< full size of vector */

  //@}

  /** @name - derived parameters */

  //@{

  double phi_min;        /**< in inflationary module, value of phi when k_min=aH */
  double phi_max;        /**< in inflationary module, value of phi when k_max=aH */

  //@}

  /** @name - technical parameters */

  //@{

  short primordial_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  //@}

  ErrorMsg error_message; /**< zone for writing error messages */

};

struct primordial_inflation_parameters_and_workspace {

  struct primordial * ppm;
  double V;
  double dV;
  double ddV;
  double a2V;
  double a2dV;
  double aH;
  double N;
  double a2ddV;
  double zpp_over_z;
  double app_over_a;
  double k;

};


/*************************************************************************************************************/

/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int primordial_spectrum_at_k(
                               struct primordial * ppm,
                               int index_md,
                               enum linear_or_logarithmic mode,
                               double k,
                               double * pk
                               );

  int primordial_init(
                      struct precision  * ppr,
                      struct perturbs   * ppt,
                      struct primordial * ppm
                      );

  int primordial_free(
                      struct primordial * ppm
                      );

  int primordial_indices(
                         struct perturbs   * ppt,
                         struct primordial * ppm
                         );

  int primordial_get_lnk_list(
                              struct primordial * ppm,
                              double kmin,
                              double kmax,
                              double k_per_decade
                              );

  int primordial_analytic_spectrum_init(
                                        struct perturbs   * ppt,
                                        struct primordial * ppm
                                        );

  int primordial_analytic_spectrum(
                                   struct primordial * ppm,
                                   int index_md,
                                   int index_ic1_ic2,
                                   double k,
                                   double * pk
                                   );

  int primordial_inflation_potential(
                                     struct primordial * ppm,
                                     double phi,
                                     double * V,
                                     double * dV,
                                     double * ddV
                                     );

  int primordial_inflation_indices(
                                   struct primordial * ppm
                                   );

  int primordial_inflation_solve_inflation(
                                           struct perturbs * ppt,
                                           struct primordial * ppm,
                                           struct precision * ppr
                                           );

  int primordial_inflation_spectra(
                                   struct perturbs * ppt,
                                   struct primordial * ppm,
                                   struct precision * ppr,
                                   double * y_ini,
                                   double * y,
                                   double * dy
                                   );

  int primordial_inflation_one_k(
                                 struct primordial * ppm,
                                 struct precision * ppr,
                                 double k,
                                 double * y,
                                 double * dy,
                                 double * curvature,
                                 double * tensor
                                 );

  int primordial_inflation_find_attractor(
                                          struct primordial * ppm,
                                          struct precision * ppr,
                                          double phi_0,
                                          double precision,
                                          double * y,
                                          double * dy,
                                          double * H_0,
                                          double * dphidt_0
                                          );

  int primordial_inflation_evolve_background(
                                             struct primordial * ppm,
                                             struct precision * ppr,
                                             double * y,
                                             double * dy,
                                             double phi_stop);

  int primordial_inflation_reach_aH(
                                    struct primordial * ppm,
                                    struct precision * ppr,
                                    double * y,
                                    double * dy,
                                    double aH_stop
                                    );

  int primordial_inflation_check_potential(
                                           struct primordial * ppm,
                                           double phi
                                           );

  int primordial_inflation_get_epsilon(
                                       struct primordial * ppm,
                                       double phi,
                                       double * epsilon
                                       );

  int primordial_inflation_derivs(
                                  double tau,
                                  double * y,
                                  double * dy,
                                  void * parameters_and_workspace,
                                  ErrorMsg error_message
                                  );

  int primordial_external_spectrum_init(
                                        struct perturbs * ppt,
                                        struct primordial * ppm
                                        );

#ifdef __cplusplus
}
#endif

/**************************************************************/

/**
 * @name Some limits imposed on parameter values:
 */

//@{

#define _K_PER_DECADE_PRIMORDIAL_MIN_ 1.

//@}

#endif
