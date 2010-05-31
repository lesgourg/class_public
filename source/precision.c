/** @file precision.c Documented common module 
 * Julien Lesgourgues, 18.04.2010    
 *
 *  Initialize the precision
 *  parameter structure, which contains all precision parameters used in
 *  all other modules. 
 */

#include "precision.h"

/** 
 * Initialize the precision parameter structure. 
 * 
 * @param ppr Input/Ouput: a precision_params structure pointer  
 * @return the error status
 *
 * All precision parameters used in the other moduels are assigned
 * here a default value, namely:
 */
int precision_init ( struct precision * ppr ) {

  /** Summary: */

  /**
   * - parameters related to the background
   */

  ppr->a_ini_over_a_today_default = 1.e-9;  /* 1.e-7 unless needs large k_max in P(k) */
  ppr->back_integration_stepsize = 2.e-2;   /* 0.02 */
  ppr->tol_background_integration = 1.e-3;  /* 0.002 */

  /**
   * - parameters related to the thermodynamics
   */

  ppr->recfast_z_initial=1.e4;
  ppr->recfast_z_final=0.;
  ppr->recfast_H_frac=1.e-3; /* from recfast */      
  ppr->recfast_x_H0_trigger=0.99; /* from recfast */   
  ppr->recfast_x_He0_trigger=0.99; /* from recfast */  
  ppr->recfast_fudge=1.14; /* from recfast */
  ppr->recfast_fudge_He=0.86; /* from recfast 1.4 */
  ppr->recfast_Heswitch=6.; /* from recfast 1.4 */ 
  ppr->recfast_Nz0=10000; /* smaller than 6000 gives bug in transfer, need to check why */
  ppr->tol_thermo_integration=1.e-3; /* optimized 9/09/08  */

  ppr->visibility_threshold_start_sources=3.5e-7; /* 3.5e-7 optimized 9/09/08  */
  ppr->visibility_threshold_free_streaming=1.e-5;

  ppr->reionization_z_start_max = 50.;
  ppr->reionization_sampling=1.e-2; /*1.e-2*/
  ppr->reionization_optical_depth_tol=1.e-2;
  ppr->reionization_exponent=1.5;
  ppr->reionization_width=0.5;
  ppr->reionization_start_factor=8.;
  ppr->helium_fullreio_redshift=3.5;
  ppr->helium_fullreio_width=0.5;

  ppr->thermo_rate_smoothing_radius=10;

  /**
   * - parameters related to the perturbations
   */

  ppr->gauge=synchronous;

  ppr->k_scalar_min=0.3; /* 0.3 -> 0.1 */
  ppr->k_scalar_oscillations=7.;  
  ppr->k_scalar_step_sub=0.1;  /* 0.02 -> 0.005 */
  ppr->k_scalar_step_super=0.005;  /* 0.01 -> 0.005 */
  ppr->k_scalar_step_transition=0.4;

  ppr->k_scalar_kmax_for_pk=1.;
  ppr->k_scalar_k_per_decade_for_pk=10.;

  ppr->k_tensor_number=15;
  ppr->k_tensor_min=1.e-4;
  ppr->k_tensor_logstep=2.;

  ppr->k_eta_min=1.e-1; /* 4.5e-6 optimized 9/09/08  */
  ppr->eta_min_over_sampling_min=0.5;
  ppr->k_eta_max=10.; /* 600 */

  ppr->l_max_g=10; /* optimized 9/09/08  */
  ppr->l_max_pol_g=10; /* optimized 9/09/08  */
  ppr->l_max_nur=25;
  ppr->l_max_g_ten=0;
  ppr->l_max_pol_g_ten=0;

  ppr->phi_ini=1.;
  ppr->entropy_ini=1.;
  ppr->gw_ini=1.;

  ppr->perturb_integration_stepsize=0.5; /* 0.5 */ 
  ppr->tol_perturb_integration=1.e-3; 
  ppr->perturb_sampling_stepsize=0.1; /* 0.1 */

  ppr->tight_coupling_trigger_eta_g_over_eta_h=0.006; /* 0.006 */
  ppr->tight_coupling_trigger_eta_g_over_eta_k=1.5e-2; /*1.5e-2*/

  ppr->rad_pert_trigger_k_over_aH = 40.; /* 40 */
  ppr->rad_pert_trigger_Omega_r = 0.1; /* 0.1 */

  /**
   * - parameter related to the Bessel functions
   */

  ppr->l_logstep=1.2 /* 1.4*/;
  ppr->l_linstep=50;
  ppr->l_max = 2500;

  ppr->bessel_scalar_x_step=0.1; /* 1. 1.27 optimized 9/09/08 */
  ppr->bessel_scalar_j_cut=1.e-5; /* 8.1e-5 optimized 9/09/08 */
  ppr->bessel_always_recompute=_FALSE_;

  /**
   * - parameter related to the primordial spectra
   */

  ppr->k_per_decade_primordial = 10.; 

  /**
   * - parameter related to the transfer functions
   */
  
  ppr->k_step_trans=0.15; /* 0.1 sampling step in k space, in units of 2pi/(eta_0-eta_rec), which is the typical period of oscillations of |Delta_l(k)|^2 */

  ppr->transfer_cut=tc_cl;
  ppr->transfer_cut_threshold_osc=0.01; /* 0.01 */
  ppr->transfer_cut_threshold_cl=2.e-6; /* 2.e-6 */

  ppr->l_scalar_max = 2500;
  ppr->l_tensor_max = 1000;

  /**
   * - automatic estimate of machine precision
   */

  get_machine_precision(&(ppr->smallest_allowed_variation));

  class_test(ppr->smallest_allowed_variation < 0,
	     ppr->error_message,
	     "smallest_allowed_variation = %e < 0",ppr->smallest_allowed_variation);

  return _SUCCESS_;

}
  
/** 
 * Computes automatically the machine precision. 
 *
 * @param smallest_allowed_variation a pointer to the smallest allowed variation
 *
 * Returns the smallest
 * allowed variation (minimum epsilon * _TOLVAR_)
 */
int get_machine_precision(double * smallest_allowed_variation) {
  double one, meps, sum;
  
  one = 1.0;
  meps = 1.0;
  do {
    meps /= 2.0;
    sum = one + meps;
  } while (sum != one);
  meps *= 2.0;
  
  *smallest_allowed_variation = meps * _TOLVAR_;

  return _SUCCESS_;

}

