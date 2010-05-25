/** @file thermodynamics.h Documented includes for thermodynamics module */

#include "background.h"

#ifndef __THERMODYNAMICS__
#define __THERMODYNAMICS__

enum reionization_parametrization {
  reio_none, /**< no reionization */
  reio_camb  /**< reionization parameterized like in CAMB */
};

enum reionization_z_or_tau {
  reio_z,  /**< input = redshift */
  reio_tau /**< input = tau */
};



/**
 * All thermodynamics-related parameters and all thermodynamical evolution.
 *
 * Once initialized by thermodynamics_init(), contains all the
 * necessary information on the thermodynamics, and in particular, a
 * table of thermodynamical quantities used for interpolation in
 * other modules
 */
struct thermo 
{
  /** @name - thermodynamics cosmological parameters */
  //@{

  double Tcmb; /**< \f$ T_{cmb} \f$ : CMB temperature */
  double YHe; /**< \f$ Y_{He} \f$ : primordial helium fraction */

  enum reionization_parametrization reio_parametrization;
  enum reionization_z_or_tau reio_z_or_tau;
  double tau_reio;
  double z_reio;

  //@}

  /** @name - all indices for the vector of thermodynamical (=th) quantities */

  //@{

  int index_th_xe;            /**< ionization fraction \f$ x_e \f$ */
  int index_th_dkappa;        /**< Thomson scattering rate \f$ d \kappa / d \eta\f$ (units 1/Mpc) */
  int index_th_ddkappa;       /**< scattering rate derivative \f$ d^2 \kappa / d \eta^2 \f$ */
  int index_th_dddkappa;      /**< scattering rate second derivative \f$ d^3 \kappa / d \eta^3 \f$ */
  int index_th_exp_m_kappa;  /**< \f$ exp^{-\kappa} \f$ */
  int index_th_g;             /**< visibility function \f$ g = (d \kappa / d \eta) * exp^{-\kappa} \f$ */
  int index_th_dg;            /**< visibility function derivative \f$ (d g / d \eta) \f$ */
  int index_th_Tb;            /**< baryon temperature \f$ T_b \f$ */
  int index_th_cb2;           /**< squared baryon sound speed \f$ c_b^2 \f$ */
  int index_th_rate;          /**< maximum variation rate of \f$ exp^{-\kappa}, g and (d g / d \eta) */
  int th_size;                /**< size of thermodynamics vector */ 

  //@}

  /** @name - thermodynamics interpolation tables */

  //@{

  int tt_size; /**< number of lines (redshift steps) in the tables */
  double * z_table; /**< values of redshift (vector of size tt_size) */
  double * thermodynamics_table; /**< all other quantities (array of size th_size*tt_size) */

  //@}

  /** @name - table of their second derivatives, used for spline interpolation */

  //@{

  double * d2thermodynamics_dz2_table; /**< values of \f$ d^2 t_i / dz^2 \f$ (array of size th_size*tt_size) */

  //@}

  /** @name - critical values of the redshift (used for switching on/off approximations in perturbations.c) */

  //@{

  double z_visibility_start_sources; /**< z below which g is non-negligible and sources should be sampled */
  double z_visibility_max; /**< z at which the visibility reaches its maximum (= recombination time) */
  double z_visibility_free_streaming; /**< z below which g is so small that radiation perturbation can be approximated by free-streaming solution */

  //@}


  /** @name - conformal time and sound horizon at recombination (used in transfer.c) */

  //@{

  double eta_rec; /**< conformal time at which the visibility reaches its maximum (= recombination time) */
  double rs_rec; /**< sound horizon at that time */

  //@}

  /** @name - flag regulating the amount of information sent to standard output (none if set to zero) */

  //@{

  short thermodynamics_verbose;

  //@}

  /** @name - zone for writing error messages */

  //@{

  ErrorMsg error_message;

  //@}

};

/**
 * All the recombination history.
 *
 * This structure is used internaly by the thermodynamics module, 
 * but never passed to other modules.
 */ 
struct recombination {

  /** @name - indices of vector of thermodynamics variables related to reionization */

  //@{

  int index_re_z;
  int index_re_xe;
  int index_re_Tb;
  int index_re_cb2;
  int index_re_dkappadeta;
  int re_size;

  //@}

  /** @name - table of the above variables at each redshift, and number of redshits */

  //@{

  double * recombination_table;
  int rt_size;

  //@}

};

/**
 * All the reionization history.
 *
 * This structure is used internaly by the thermodynamics module, 
 * but never passed to other modules.
 */ 
struct reionization {

  /** @name - indices of vector of thermodynamics variables related to reionization */

  //@{

  int index_re_z;
  int index_re_xe;
  int index_re_Tb;
  int index_re_cb2;
  int index_re_dkappadeta;
  int index_re_dkappadz;
  int index_re_d3kappadz3;
  int re_size;

  //@}

  /** @name - table of the above variables at each redshift, and number of redshits */

  //@{

  double * reionization_table;
  int rt_size;

  //@}

  /** @name - reionization optical depth */

  //@{

  double reionization_optical_depth;

  //@}

  /** @name - indices describing parameters used in the definition of the various possible functions x_e(z) */

  //@{

  int index_reio_redshift;
  int index_reio_start;
  int index_reio_xe_before;
  int index_reio_xe_after;
  int index_reio_exponent;
  int index_reio_width;
  int index_helium_fullreio_fraction;
  int index_helium_fullreio_redshift;
  int index_helium_fullreio_width;

  //@}

  /** @name - vector of such parameters, and its size */

  double * reionization_parameters;
  int reio_num_params;

  //@}

  /** @name - index of line in recombination table corresponding to first line of reionization table */

  //@{

  int index_reco_when_reio_start;

  //@}

};

/**************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int thermodynamics_at_z(
			  double z,
			  enum interpolation_mode intermode,
			  int * last_index,
			  double * pvecthermo_local
			  );

  int thermodynamics_init(
			  struct background * pco_input,
			  struct precision * ppp_input,
			  struct thermo * pth_output
			  );

  int thermodynamics_free();

  int thermodynamics_indices(
			     struct thermo * pthermo,
			     struct recombination * preco,
			     struct reionization * preio
			     );

  int thermodynamics_reionization_function(
					   double z,
					   struct reionization * preio,
					   double * xe
					   );

  int thermodynamics_reionization(
				  struct recombination * preco,
				  struct reionization * preio
				  );

  int thermodynamics_reionization_discretize(
					     struct recombination * preco,
					     struct reionization * preio
					     );

  int thermodynamics_get_xe_before_reionization(double z,
						struct recombination * preco,
						double * xe);

  int thermodynamics_recombination(
				   struct recombination * prec
				   );

/*   int thermodynamics_cure_discontinuity( */
/* 					int index_th */
/* 					); */

  int thermodynamics_derivs_with_recfast(
					 double z,
					 double * y,
					 double * dy,
					 void * fixed_parameters,
					 ErrorMsg error_message
					 );

#ifdef __cplusplus
}
#endif

/**************************************************************/

/**  
 * @name Some basic constants needed by RECFAST:
 */

//@{

#define _C_ 2.99792458e8
#define _k_B_ 1.380658e-23
#define _h_P_ 6.6260755e-34
#define _m_e_ 9.1093897e-31
#define _m_H_ 1.673575e-27  /*updated like in recfast 1.4*/
#define _sigma_ 6.6524616e-29
#define _a_ 7.565914e-16
#define _G_ 6.6742e-11  /*updated like in recfast 1.4*/
#define _m_p_ 1.672623e-27
#define _not4_ 3.9715  /*updated like in recfast 1.4*/
#define _Mpc_in_sec_ 1.029272e14

//@}

/**  
 * @name Some specific constants needed by RECFAST:
 */

//@{

#define _RECFAST_INTEG_SIZE_ 3

#define _Lambda_ 8.2245809
#define _Lambda_He_ 51.3
#define _L_H_ion_ 1.096787737e7
#define _L_H_alpha_ 8.225916453e6
#define _L_He1_ion_ 1.98310772e7
#define _L_He2_ion_ 4.389088863e7
#define _L_He_2s_ 1.66277434e7
#define _L_He_2p_ 1.71134891e7
#define _bigH_ 100.0e3/(1.0e6*3.0856775807e16)*2999.7
#define	_A2P_s_		1.798287e9  /*updated like in recfast 1.4*/
#define	_A2P_t_		177.58e0  /*updated like in recfast 1.4*/
#define	_L_He_2Pt_	1.690871466e7  /*updated like in recfast 1.4*/
#define	_L_He_2St_	1.5985597526e7  /*updated like in recfast 1.4*/
#define	_L_He2St_ion_	3.8454693845e6  /*updated like in recfast 1.4*/
#define	_sigma_He_2Ps_	1.436289e-22  /*updated like in recfast 1.4*/
#define	_sigma_He_2Pt_	1.484872e-22  /*updated like in recfast 1.4*/

//@}

/**  
 * @name Some specific constants needed by recfast_derivs:
 */

//@{

#define _a_PPB_ 4.309
#define _b_PPB_ -0.6166
#define _c_PPB_ 0.6703
#define _d_PPB_ 0.5300
#define _T_0_ pow(10.,0.477121) /* from recfast 1.4 */
#define _a_VF_ pow(10.,-16.744)
#define _b_VF_ 0.711
#define _T_1_ pow(10.,5.114)
#define _C2p3P_ 1.69087
#define _C2p1P_ 1.71135
#define _A2p3P_ 233.
#define _g3P_ 3.
#define	_a_trip_ pow(10.,-16.306) /* from recfast 1.4 */
#define	_b_trip_ 0.761 /* from recfast 1.4 */

//@}

/**  
 * @name Some limits imposed on cosmological parameter values:
 */

//@{

#define _TCMB_BIG_ 2.8     /**< maximal \f$ T_{cmb} \f$ in K */
#define _TCMB_SMALL_ 2.7  /**< minimal \f$ T_{cmb}  \f$ in K */
#define _YHE_BIG_ 0.5     /**< maximal \f$ Y_{He} \f$ */
#define _YHE_SMALL_ 0.01  /**< minimal \f$ Y_{He} \f$ */

//@}

#endif
