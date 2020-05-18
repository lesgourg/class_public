/** @file primordial.h Documented includes for primordial module. */

#ifndef __PRIMORDIAL__
#define __PRIMORDIAL__

#include "perturbations.h"

/** enum defining how the primordial spectrum should be computed */

enum primordial_spectrum_type {
  analytic_Pk,
  two_scales,
  inflation_V,
  inflation_H,
  inflation_V_end,
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
  natural,
  higgs_inflation
};

/** enum defining which quantity plays the role of a target for evolving inflationary equations */

enum target_quantity {
  _aH_,
  _phi_,
  _end_inflation_,
  _a_
};

/** enum specifying if we want to integrate equations forward or backward in time */

enum integration_direction {
  backward,
  forward
};

/** enum specifying if we want to evolve quantities with conformal or proper time */

enum time_definition {
  conformal,
  proper
};

/** enum specifying how, in the inflation_V_end case, the value of phi_pivot should calculated */

enum phi_pivot_methods {
  N_star,
  ln_aH_ratio,
  ln_aH_ratio_auto
};

/** enum specifying how the inflation module computes the primordial spectrum (default: numerical) */

enum inflation_module_behavior {
  numerical,
  analytical
};

/**
 * Structure containing everything about primordial spectra that other modules need to know.
 *
 * Once initialized by primordial_init(), contains a table of all
 * primordial spectra as a function of wavenumber, mode, and pair of initial conditions.
 */

struct primordial {

  /** @name - input parameters initialized by user in input module
      (all other quantities are computed in this module, given these parameters
      and the content of the 'precision' and 'perturbs' structures) */

  //@{

  double k_pivot; /**< pivot scale in \f$ Mpc^{-1} \f$ */

  enum primordial_spectrum_type primordial_spec_type; /**< type of primordial spectrum (simple analytic from, integration of inflationary perturbations, etc.) */

  /* - parameters describing the case primordial_spec_type = analytic_Pk : amplitudes, tilts, runnings, cross-correlations, ... */

  double A_s;  /**< usual scalar amplitude = curvature power spectrum at pivot scale */
  double n_s;  /**< usual scalar tilt = [curvature power spectrum tilt at pivot scale -1] */
  double alpha_s; /**< usual scalar running */
  
  double r;    /**< usual tensor to scalar ratio of power spectra, \f$ r=A_T/A_S=P_h/P_R \f$*/
  double n_t;  /**< usual tensor tilt = [GW power spectrum tilt at pivot scale] */
  double alpha_t; /**< usual tensor running */

  double f_bi;  /**< baryon isocurvature (BI) entropy-to-curvature ratio \f$ S_{bi}/R \f$*/
  double n_bi;  /**< BI tilt */
  double alpha_bi; /**< BI running */

  double f_cdi;  /**< CDM isocurvature (CDI) entropy-to-curvature ratio \f$ S_{cdi}/R \f$*/
  double n_cdi;  /**< CDI tilt */
  double alpha_cdi; /**< CDI running */

  double f_nid;  /**< neutrino density isocurvature (NID) entropy-to-curvature ratio \f$ S_{nid}/R \f$*/
  double n_nid;  /**< NID tilt */
  double alpha_nid; /**< NID running */

  double f_niv;  /**< neutrino velocity isocurvature (NIV) entropy-to-curvature ratio \f$ S_{niv}/R \f$*/
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

  /** parameters describing the case primordial_spec_type = inflation_V */

  enum potential_shape potential;

  double V0;	/**< one parameter of the function V(phi) */
  double V1;	/**< one parameter of the function V(phi) */
  double V2;	/**< one parameter of the function V(phi) */
  double V3;	/**< one parameter of the function V(phi) */
  double V4;	/**< one parameter of the function V(phi) */

  /* parameters describing the case primordial_spec_type = inflation_H */

  double H0;	/**< one parameter of the function H(phi) */
  double H1;	/**< one parameter of the function H(phi) */
  double H2;	/**< one parameter of the function H(phi) */
  double H3;	/**< one parameter of the function H(phi) */
  double H4;	/**< one parameter of the function H(phi) */

  /* parameters describing inflation_V_end */

  double phi_end;	/**< value of inflaton at the end of inflation */
  enum phi_pivot_methods phi_pivot_method; /**< flag for method used to define and find the pivot scale */
  double phi_pivot_target; /**< For each of the above methods, critical value to be reached between pivot and end of inflation (N_star, [aH]ratio, etc.) */

  /* behavior of the inflation module */
  enum inflation_module_behavior behavior; /**< Specifies if the inflation module computes the primordial spectrum numerically (default) or analytically*/

  /** 'external_Pk' mode: command generating the table of Pk and custom parameters to be passed to it */

  char*  command;  /**< string with the command for calling 'external_Pk' */
  double custom1;  /**< one parameter of the primordial computed in 'external_Pk' */
  double custom2;  /**< one parameter of the primordial computed in 'external_Pk' */
  double custom3;  /**< one parameter of the primordial computed in 'external_Pk' */
  double custom4;  /**< one parameter of the primordial computed in 'external_Pk' */
  double custom5;  /**< one parameter of the primordial computed in 'external_Pk' */
  double custom6;  /**< one parameter of the primordial computed in 'external_Pk' */
  double custom7;  /**< one parameter of the primordial computed in 'external_Pk' */
  double custom8;  /**< one parameter of the primordial computed in 'external_Pk' */
  double custom9;  /**< one parameter of the primordial computed in 'external_Pk' */
  double custom10; /**< one parameter of the primordial computed in 'external_Pk' */

  //@}


  /** @name - technical parameters */

  //@{

  short primordial_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  //@}

};

class PrimordialModule;
struct primordial_inflation_parameters_and_workspace {
  primordial_inflation_parameters_and_workspace(PrimordialModule* p_m) : primordial_module(p_m) {}
  const PrimordialModule* primordial_module;
  double N;
  double a2;

  double V;
  double dV;
  double ddV;
  double aH;

  double H;
  double dH;
  double ddH;
  double dddH;

  double zpp_over_z;
  double app_over_a;

  double k;

  enum integration_direction integrate;
  enum time_definition time;

};

/**
 * @name Some limits imposed on parameter values:
 */

//@{

#define _K_PER_DECADE_PRIMORDIAL_MIN_ 1.

//@}

#endif
/* @endcond */
