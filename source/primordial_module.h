#ifndef PRIMORDIAL_MODULE_H
#define PRIMORDIAL_MODULE_H

#include "input_module.h"
#include "perturbations_module.h"
#include "base_module.h"

class PrimordialModule : public BaseModule {
public:
  PrimordialModule(InputModulePtr input_module, PerturbationsModulePtr perturbation_module);
  ~PrimordialModule();

  int primordial_spectrum_at_k(int index_md, enum linear_or_logarithmic mode, double k, double* pk) const;
  int primordial_output_titles(char titles[_MAXTITLESTRINGLENGTH_]) const;
  int primordial_output_data(int number_of_titles, double* data) const;

  int* ic_size_;    /**< for a given mode, ic_size[index_md] = number of initial conditions included in computation */
  int* ic_ic_size_; /**< number of ordered pairs of (index_ic1, index_ic2); this number is just N(N+1)/2  where N = ic_size[index_md] */
  short** is_non_zero_; /**< is_non_zero[index_md][index_ic1_ic2] set to false if pair
                          (index_ic1, index_ic2) is uncorrelated
                          (ensures more precision and saves time with respect to the option
                          of simply setting P(k)_(index_ic1, index_ic2) to zero) */
  int lnk_size_;    /**< number of ln(k) values */

  /** @name - derived parameters */
  //@{
  double phi_pivot_;      /**< in inflationary module, value of
                            phi_pivot (set to 0 for inflation_V,
                            inflation_H; found by code for
                            inflation_V_end) */
  double phi_min_;        /**< in inflationary module, value of phi when \f$ k_{min}=aH \f$*/
  double phi_max_;        /**< in inflationary module, value of phi when \f$ k_{max}=aH \f$*/
  double phi_stop_;       /**< in inflationary module, value of phi at the end of inflation */
  //@}
  double A_s_;
  double n_s_;
  double alpha_s_;
  double beta_s_;
  double r_;
  double n_t_;
  double alpha_t_;


private:
  int primordial_init();
  int primordial_free();
  int primordial_indices();
  int primordial_get_lnk_list(double kmin, double kmax, double k_per_decade);
  int primordial_analytic_spectrum_init();
  int primordial_analytic_spectrum(int index_md, int index_ic1_ic2, double k, double* pk) const;
  int primordial_inflation_potential(double phi, double* V, double* dV, double* ddV) const;
  int primordial_inflation_hubble(double phi, double* H, double* dH, double* ddH, double* dddH) const;
  int primordial_inflation_indices();
  int primordial_inflation_solve_inflation();
  int primordial_inflation_analytic_spectra(double* y_ini);
  int primordial_inflation_spectra(double* y_ini);
  int primordial_inflation_one_wavenumber(double* y_ini, int index_k);
  int primordial_inflation_one_k(double k, double* y, double* dy, double* curvature, double* tensor);
  int primordial_inflation_find_attractor(double phi_0, double precision, double* y, double* dy, double* H_0, double* dphidt_0);
  int primordial_inflation_evolve_background(double* y, double* dy, enum target_quantity target, double stop, short check_epsilon, enum integration_direction direction, enum time_definition time);
  int primordial_inflation_check_potential(double phi, double* V, double* dV, double* ddV);
  int primordial_inflation_check_hubble(double phi, double* H, double* dH, double* ddH, double* dddH);
  int primordial_inflation_get_epsilon(double phi, double* epsilon);
  int primordial_inflation_find_phi_pivot(double* y, double* dy);
  int primordial_inflation_derivs_member(double tau, double* y, double* dy, void * parameters_and_workspace, ErrorMsg error_message) const;
  static int primordial_inflation_derivs(double tau, double* y, double* dy, void * parameters_and_workspace, ErrorMsg error_message);
  int primordial_external_spectrum_init();

  PerturbationsModulePtr perturbations_module_;

  /** @name - pre-computed table of primordial spectra, and related quantities */

  //@{
  int md_size_;      /**< number of modes included in computation */
  double* lnk_;    /**< list of ln(k) values lnk[index_k] */
  double** lnpk_;  /**< depends on indices index_md, index_ic1, index_ic2, index_k as:
                      lnpk[index_md][index_k*ic_ic_size_[index_md]+index_ic1_ic2]
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
  double** ddlnpk_; /**< second derivative of above array, for spline interpolation. So:
                       - for index_ic1 = index_ic, we spline ln[P(k)] vs. ln(k), which is
                       good since this function is usually smooth.
                       - for non-diagonal coefficients, we spline
                       P(k)_(index_ic1, index_ic2)/sqrt[P(k)_index_ic1 P(k)_index_ic2]
                       vs. ln(k), which is fine since this quantity is often assumed to be
                       constant (e.g for fully correlated/anticorrelated initial conditions)
                       or nearly constant, and with arbitrary sign.
                    */

  //@}

  //@{
  /** @name - parameters describing the case primordial_spec_type = analytic_Pk : amplitudes, tilts, runnings, cross-correlations, ... */
  double** amplitude_; /**< all amplitudes in matrix form: amplitude[index_md][index_ic1_ic2] */
  double** tilt_;      /**< all tilts in matrix form: tilt[index_md][index_ic1_ic2] */
  double** running_;   /**< all runnings in matrix form: running[index_md][index_ic1_ic2] */
  //@}

  //@{
  /** @name - for the inflation simulator, indices in vector of
      background/perturbation */
  int index_in_a_;       /**< scale factor */
  int index_in_phi_;     /**< inflaton vev */
  int index_in_dphi_;    /**< its time derivative */
  int index_in_ksi_re_;  /**< Mukhanov variable (real part) */
  int index_in_ksi_im_;  /**< Mukhanov variable (imaginary part) */
  int index_in_dksi_re_; /**< Mukhanov variable (real part, time derivative) */
  int index_in_dksi_im_; /**< Mukhanov variable (imaginary part, time derivative) */
  int index_in_ah_re_;   /**< tensor perturbation (real part) */
  int index_in_ah_im_;   /**< tensor perturbation (imaginary part) */
  int index_in_dah_re_;  /**< tensor perturbation (real part, time derivative) */
  int index_in_dah_im_;  /**< tensor perturbation (imaginary part, time derivative) */
  int in_bg_size_;       /**< size of vector of background quantities only */
  int in_size_;          /**< full size of vector */
  //@}

};

#endif //PRIMORDIAL_MODULE_H
