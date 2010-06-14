/** @file precision.h Documented includes for precision module */

#ifndef __PRECISION__
#define __PRECISION__

#include "common.h"

/** 
 * List of methods for stopping the transfer function computation 
 * at a given k for each l (saves lots of time). 
 */

enum transfer_cutting {
  tc_none, /**< no transfer cut: for given l, compute transfer functions over full k range (long and usually useless) */ 
  tc_osc, /**< transfer cut with oscillation method: for given l, compute transfer functions until k_max such that oscillations of \f$ \Delta_l(k) \f$ are tiny relatively to largest oscillation */ 
  tc_cl /**< transfer cut with Cl variation method: for given l, compute transfer functions until k_max such that C_l's variation is tiny (C_l being computed approximately and with flat spectrum)  */ 
};

/** 
 * List of coded gauges.
 *
 * More gauges can in principle be defined and coded. 
 */
enum possible_gauges {
  newtonian, /**< newtonian (or longitudinal) gauge */
  synchronous /**< synchronous gauge with \f$ \theta_{cdm} = 0 \f$ by convention */
};


/**
 * All precision parameters. 
 *  
 * Includes integrations
 * steps, flags telling how the computation is to be performed, etc.
 */
struct precision
{

  /** @name - parameters related to the background */
  //@{

  /**
   * default initial value of scale factor in background integration, in
   * units of scale factor today
   */ 
  double a_ini_over_a_today_default; 

  /** 
   * default step d eta in background integration, in units of 
   * conformal Hubble time (\f$ d eta \f$ = back_integration_stepsize / aH )
   */
  double back_integration_stepsize; 

  /**
   * parameter controlling precision of background integration
   */
  double tol_background_integration;

  //@}

  /** @name - parameters related to the thermodynamics */

  //@{

  /** - for recombination */

  double recfast_z_initial; /**< recfast parameter */
  double recfast_z_final; /**< recfast parameter */
  double recfast_x_H0; /**< recfast parameter */
  double recfast_x_He0; /**< recfast parameter */
  double recfast_H_frac; /**< recfast parameter */
  double recfast_fudge; /**< recfast parameter */
  double recfast_fudge_He; /**< recfast 1.4 parameter */
  short recfast_Heswitch; /**< recfast 1.4 parameter */
  double recfast_x_H0_trigger;  /**< recfast parameter */
  double recfast_x_He0_trigger;  /**< recfast parameter */
  int recfast_Nz0; /**< recfast parameter */
  double tol_thermo_integration; /**< recfast parameter */

  /** - for reionization */

  double reionization_z_start_max; /**< maximum redshift at which reionization should start. If not, return an error. */
  double reionization_sampling; /**< control stepsize in z during reionization */
  double reionization_optical_depth_tol; /**< fractional error on optical_depth */
  double reionization_exponent; /**< parameter for CAMB-like parametrization */
  double reionization_width; /**< parameter for CAMB-like parametrization */
  double reionization_start_factor; /**< parameter for CAMB-like parametrization */
  double helium_fullreio_redshift; /**< parameter for CAMB-like parametrization */
  double helium_fullreio_width; /**< parameter for CAMB-like parametrization */
   
  /** - general */

  int thermo_rate_smoothing_radius; /**< plays a minor (almost aesthetic) role in the definition of the variation rate of thermodynamical quantities */

  /**
   * critical values of visibility function g: 
   *
   * (a) when g becomes larger than visibility_threshold_start_sources, start sampling the sources
   */
  double visibility_threshold_start_sources;

  /**
   * (b) when g becomes smaller than
   * visibility_threshold_free_streaming, try to turn off radiation
   * perturbations (actually this condition is overwritten by another
   * one: Omega_r < rad_pert_trigger_Omega_r; the present variable
   * remains just in case one would need it)
   */
  double visibility_threshold_free_streaming;

  //@}

  /** @name - parameters related to the perturbation */

  //@{

  /**
   * gauge in which to perform the calculation
   */
  enum possible_gauges gauge; 

  double k_scalar_min; /**< first mode k_min in units of Hubble radius today */
  double k_scalar_oscillations; /**< number of acoustic oscillations experienced by last mode k_max (to resolve n acoustic peaks, one should set this number to more than n) */
  double k_scalar_step_sub; /**< step in k space, in units of one period of acoustic oscillation at decoupling, for scales inside sound horizon at decoupling */
  double k_scalar_step_super; /**< step in k space, in units of one period of acoustic oscillation at decoupling, for scales above sound horizon at decoupling */  double k_scalar_step_transition; /**< dimensionless number regulaing the transition fro _sub step to _super step. Decrease for more precision. */

  double k_scalar_k_per_decade_for_pk; /**< if values needed between kmax inferred from k_scalar_oscillations and k_scalar_kmax_for_pk, this gives the number of k per decade */

  int k_tensor_number;  /**< number of k values (tensor modes) */
  double k_tensor_min; /**< min k value (tensor modes) */
  double k_tensor_logstep; /**< logstep for k sampling (tensor modes) */

  double k_eta_min; /**< default value of \f$ k \times \eta \f$ for starting integration of given wavenumber */

  double eta_min_over_sampling_min; /**< maximum ratio of eta_min (when perturbations start to be integrated) over the first time at which source are sampled. Value important for the smallest k, otherwise this condition is overwritten by the one on \f$ k \times \eta \f$ */ 

  double k_eta_max; /**< maximum value of \f$ k \times \eta \f$ for which the sources should be computed in a pure CMB run*/

  double l_max_g;     /**< number of momenta in Boltzmann hierarchy for photon temperature (scalar), at least 4 */
  double l_max_pol_g; /**< number of momenta in Boltzmann hierarchy for photon polarisation (scalar), at least 4 */
  double l_max_nur;   /**< number of momenta in Boltzmann hierarchy for relativistic neutrino/relics (scalar), at least 4 */
  double l_max_g_ten;     /**< number of momenta in Boltzmann hierarchy for photon temperature (tensor), at least 4 */
  double l_max_pol_g_ten; /**< number of momenta in Boltzmann hierarchy for photon polarisation (tensor), at least 4 */

  double phi_ini;     /**< initial condition for Bardeen potential for adiabatic */
  double entropy_ini; /**< initial condition for entropy perturbation for isocurvature */ 
  double gw_ini;      /**< initial condition for tensor metric perturbation h */

  /** 
   * default step \f$ d \eta \f$ in perturbation integration, in units of the timescale involved in the equations (usally, the min of \f$ 1/k \f$, \f$ 1/aH \f$, \f$ 1/\dot{\kappa} \f$) 
   */
  double perturb_integration_stepsize;

  /** 
   * default step \f$ d \eta \f$ for sampling the source function, in units of the timescale involved in the sources: \f$ (\dot{\kappa}- \ddot{\kappa}/\dot{\kappa})^{-1} \f$
   */
  double perturb_sampling_stepsize;

  /** 
   * control parameter for the precision of the perturbation integration 
   */
  double tol_perturb_integration;

  /**
   * when to switch off tight-coupling approximation:
   * first condition: \f$ \eta_g/\eta_H \f$ <
   * tight_coupling_trigger_eta_g_over_eta_h 
   */
  double tight_coupling_trigger_eta_g_over_eta_h;

  /**
   * when to switch off tight-coupling approximation:
   * second condition: \f$ \eta_g/\eta_k \equiv k \eta_g \f$ <
   * tight_coupling_trigger_eta_g_over_eta_k
   */
  double tight_coupling_trigger_eta_g_over_eta_k;

  /**
   * when to switch off radiation perturbations, ie when to switch
   * on free-streaming approximation (keep density and theta, set
   * shear and higher momenta of ultrarelativistic particles to zero):
   * first condition: \f$ k/aH \f$ > rad_pert_trigger_k_over_aH
   */
  double rad_pert_trigger_k_over_aH;

  /**
   * when to switch off radiation perturbations, ie when to switch
   * on free-streaming approximation (keep density and theta, set
   * shear and higher momenta of ultrarelativistic particles to zero):
   * second condition: \f$ \Omega_r \f$ < rad_pert_trigger_Omega_r
   */
  double rad_pert_trigger_Omega_r;

  //@}

  /** @name - parameters related to Bessel functions */

  //@{

  int l_linstep; /**< factor for logarithmic spacing of values of l over which transfer functions are sampled */
  double l_logstep; /**< maximum spacing of values of l over which transfer functions are sampled (so, spacing becomes linear instead of logarithmic at some point) */

  double bessel_scalar_x_step; /**< step dx for sampling Bessel functions \f$ j_l(x) \f$ */
  double bessel_scalar_j_cut; /**< value of \f$ j_l \f$ below which it is approximated by zero (in the region \f$ x \ll l \f$) */
  short bessel_always_recompute; /**< if set to true, Bessels are never read from file */

  //@}

  /** @name - parameters related to the primordial spectra */

  //@{

  double k_per_decade_primordial; /**< logarithmic sampling for primordial spectra (number of points per decade in k space) */

  //@}

  /** @name - parameters related to the transfer function */

  //@{

  double k_step_trans; /**< sampling step in k space, in units of \f$ 2\pi/(\eta_0-\eta_{rec}) \f$, which is the typical period of oscillations of \f$ \Delta_l(k) \f$ */

  enum transfer_cutting transfer_cut; /**< flag telling how to cut automatically the transfer function computation at a given \f$ k_{max} \f$ value */

  double transfer_cut_threshold_osc; /**< threshold used for cutting the transfer function computation at a given \f$ k_{max} \f$ value, if transfer_cut = _TC_OSC_ (oscillation method: for given l, compute transfer functions until k_max such that oscillations of \f$ \Delta_l(k) \f$ are tiny relatively to largest oscillation) */

  double transfer_cut_threshold_cl; /**< threshold used for cutting the transfer function computation at a given \f$ k_{max} \f$ value, if transfer_cut = _TC_CL_ (Cl variation method: for given l, compute transfer functions until k_max such that C_l's variation is tiny, C_l being computed approximately and with flat spectrum) */

  double smallest_allowed_variation; /**< machine-dependent, assigned automatically by the code */

  //@}

  /** @name - zone for writing error messages */

  //@{

  ErrorMsg error_message;

  //@}

};

/**************************************************************/

/* Boilerplate for C++ */

#ifdef __cplusplus
extern "C" {
#endif

  int precision_init (struct precision * ppp);

#ifdef __cplusplus
}
#endif


#endif
