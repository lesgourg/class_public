/** @file common.h Generic libraries, parameters and functions used in the whole code. */

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "float.h"
#include "svnversion.h"
#ifdef _OPENMP
#include "omp.h"
#endif

#ifndef __COMMON__
#define __COMMON__

#define _VERSION_ "v1.7"

#define _TRUE_ 1 /**< integer associated to true statement */
#define _FALSE_ 0 /**< integer associated to false statement */

#define _SUCCESS_ 0 /**< integer returned after successfull call of a function */
#define _FAILURE_ 1 /**< integer returned after failure in a function */

#define _ERRORMSGSIZE_ 2048 /**< generic error messages are cut beyond this number of characters */
typedef char ErrorMsg[_ERRORMSGSIZE_]; /**< Generic error messages (there is such a field in each structure) */

#define _FILENAMESIZE_ 256 /**< size of the string read in each line of the file (extra characters not taken into account) */
typedef char FileName[_FILENAMESIZE_];

#define _PI_ 3.1415926535897932384626433832795e0 /**< The number pi */

#define _MAX_IT_ 10000/**< default maximum number of iterations in conditional loops (to avoid infinite loops) */

#define _QUADRATURE_MAX_ 250 /**< maximum allowed number of abssices in quadrature integral estimation */

#define _QUADRATURE_MAX_BG_ 800 /**< maximum allowed number of abssices in quadrature integral estimation */

#define _TOLVAR_ 100. /**< The minimum allowed variation is the machine precision times this number */

#define _HUGE_ 1.e99

#define min(a,b) (((a)<(b)) ? (a) : (b) ) /**< the usual "min" function */
#define max(a,b) (((a)<(b)) ? (b) : (a) ) /**< the usual "max" function */
#define sign(a) (((a)>0) ? 1. : -1. )

#define index_symmetric_matrix(i1,i2,N) (((i1)<=(i2)) ? (i2+N*i1-(i1*(i1+1))/2) : (i1+N*i2-(i2*(i2+1))/2)) /**< assigns an index from 0 to [N(N+1)/2-1] to the coefficients M_{i1,i2} of an N*N symmetric matrix; useful for converting a symmetric matrix to a vector, without loosing or double-counting any information */

/* macro for calling function and returning error if it failed */
#define class_call(function,						\
		   error_message_from_function,				\
		   error_message_output)				\
  do {									\
    if (function == _FAILURE_) {					\
      ErrorMsg Transmit_Error_Message;					\
      sprintf(Transmit_Error_Message,"%s(L:%d) : error in %s;\n=>%s",	\
	      __func__,__LINE__,#function,error_message_from_function);	\
      sprintf(error_message_output,"%s",Transmit_Error_Message);	\
      return _FAILURE_;							\
    }									\
  } while(0);

/* macro for calling function and returning error if it failed */
#define class_call_except(function,					\
			  error_message_from_function,			\
			  error_message_output,				\
			  list_of_commands_without_commas)		\
  do {									\
    if (function == _FAILURE_) {					\
      ErrorMsg Transmit_Error_Message;					\
      sprintf(Transmit_Error_Message,"%s(L:%d) : error in %s;\n=>%s",	\
	      __func__,__LINE__,#function,error_message_from_function);	\
      sprintf(error_message_output,"%s",Transmit_Error_Message);	\
      list_of_commands_without_commas;					\
      return _FAILURE_;							\
    }									\
  } while(0);

/* same in parallel region */
#define class_call_parallel(function,					\
			    error_message_from_function,		\
			    error_message_output)			\
  do {									\
    if (abort == _FALSE_) {						\
      if (function == _FAILURE_) {					\
	ErrorMsg Transmit_Error_Message;				\
	sprintf(Transmit_Error_Message,"%s(L:%d) : error in %s;\n=>%s",	\
		__func__,__LINE__,#function,error_message_from_function); \
	sprintf(error_message_output,"%s",Transmit_Error_Message);	\
	abort=_TRUE_;							\
      }									\
    }									\
  } while(0);

/* macro for testing condition and returning error if condition is true;
   args is a variable list of optional arguments, e.g.: args="x=%d",x 
   args cannot be empty, if there is nothing to pass use args="" */
#define class_test(condition,						\
		   error_message_output,				\
		   args...)						\
  do {									\
    if (condition) {							\
      ErrorMsg Transmit_Error_Message;					\
      ErrorMsg Optional_arguments;					\
      sprintf(Transmit_Error_Message,					\
	      "%s(L:%d) : condition (%s) is true",			\
	      __func__,__LINE__,#condition);				\
      sprintf(Optional_arguments,args);					\
      sprintf(error_message_output,"%s; %s",				\
	      Transmit_Error_Message, Optional_arguments);		\
      return _FAILURE_;							\
    }									\
  } while(0);

/* macro for testing condition and returning error if condition is true;
   args is a variable list of optional arguments, e.g.: args="x=%d",x 
   args cannot be empty, if there is nothing to pass use args="" */
#define class_test_except(condition,					\
			  error_message_output,				\
			  list_of_commands_without_commas,		\
			  args...)					\
  do {									\
    if (condition) {							\
      ErrorMsg Transmit_Error_Message;					\
      ErrorMsg Optional_arguments;					\
      sprintf(Transmit_Error_Message,					\
	      "%s(L:%d) : condition (%s) is true",			\
	      __func__,__LINE__,#condition);				\
      sprintf(Optional_arguments,args);					\
      sprintf(error_message_output,"%s; %s",				\
	      Transmit_Error_Message, Optional_arguments);		\
      list_of_commands_without_commas;					\
      return _FAILURE_;							\
    }									\
  } while(0);


/* macro for returning error message;
   args is a variable list of optional arguments, e.g.: args="x=%d",x 
   args cannot be empty, if there is nothing to pass use args="" */
#define class_stop(error_message_output,				\
		   args...)						\
  do {									\
    if (_TRUE_) {							\
      ErrorMsg Transmit_Error_Message;					\
      ErrorMsg Optional_arguments;					\
      sprintf(Transmit_Error_Message,					\
	      "%s(L:%d) : error",					\
	      __func__,__LINE__);					\
      sprintf(Optional_arguments,args);					\
      sprintf(error_message_output,"%s; %s",				\
	      Transmit_Error_Message, Optional_arguments);		\
      return _FAILURE_;							\
    }									\
  } while(0);

/* same in parallel region */
#define class_test_parallel(condition,					\
		   error_message_output,				\
		   args...)						\
  do {									\
    if (abort == _FALSE_) {						\
      if (condition) {							\
	ErrorMsg Transmit_Error_Message;				\
	ErrorMsg Optional_arguments;					\
	sprintf(Transmit_Error_Message,					\
		"%s(L:%d) : condition (%s) is true",			\
		__func__,__LINE__,#condition);				\
	sprintf(Optional_arguments,args);				\
	sprintf(error_message_output,"%s; %s",				\
		Transmit_Error_Message, Optional_arguments);		\
	abort=_TRUE_;							\
      }									\
    }									\
  } while(0);

/* macro for allocating memory and returning error if it failed */
#define class_alloc(pointer,						\
		    size,						\
		    error_message_output)				\
  do {									\
    pointer=malloc(size);						\
    if (pointer == NULL) {						\
      ErrorMsg Transmit_Error_Message;					\
      int size_int;							\
      size_int=size;							\
      sprintf(Transmit_Error_Message,					\
	      "%s(L:%d) : could not allocate %s with size %d",		\
	      __func__,__LINE__,					\
	      #pointer,size_int);					\
      sprintf(error_message_output,"%s",Transmit_Error_Message);	\
      return _FAILURE_;							\
    }									\
  } while(0);

/* same inside parallel structure */
#define class_alloc_parallel(pointer,					\
		    size,						\
		    error_message_output)				\
  do {									\
    if (abort == _FALSE_) {						\
      pointer=malloc(size);						\
      if (pointer == NULL) {						\
	int size_int;							\
	size_int=size;							\
	ErrorMsg Transmit_Error_Message;				\
	sprintf(Transmit_Error_Message,					\
		"%s(L:%d) : could not allocate %s with size %d",	\
		__func__,__LINE__,					\
		#pointer,size_int);					\
	sprintf(error_message_output,"%s",Transmit_Error_Message);	\
	abort=_TRUE_;							\
      }									\
    }									\
    else {								\
      pointer=NULL;							\
    }									\
  } while(0);

/* macro for allocating memory, initializing it with zeros/ and returning error if it failed */
#define class_calloc(pointer,						\
		     number,						\
		     size,						\
		     error_message_output)				\
  do {									\
    pointer=calloc(number,size);					\
    if (pointer == NULL) {						\
      ErrorMsg Transmit_Error_Message;					\
      int size_int;							\
      size_int=number*size;						\
      sprintf(Transmit_Error_Message,					\
	      "%s(L:%d) : could not allocate %s with size %d",		\
	      __func__,__LINE__,					\
	      #pointer,size_int);					\
      sprintf(error_message_output,"%s",Transmit_Error_Message);	\
      return _FAILURE_;							\
    }									\
  } while(0);

/* macro for opening file and returning error if it failed */
#define class_open(pointer,						\
		   filename,						\
  	           mode,						\
		   error_message_output)				\
  do {									\
    pointer=fopen(filename,mode);					\
    if (pointer == NULL) {						\
      ErrorMsg Transmit_Error_Message;					\
      sprintf(Transmit_Error_Message,					\
	      "%s(L:%d) : could not open %s with name %s and mode %s",	\
	      __func__,__LINE__,					\
	      #pointer,filename,#mode);					\
      sprintf(error_message_output,"%s",Transmit_Error_Message);	\
      return _FAILURE_;							\
    }									\
  } while(0);

/* macro for defining indices (usually one, sometimes a block) */
#define class_define_index(index,					\
                           condition,					\
                           running_index,                               \
                           number_of_indices)                           \
  do {									\
    if (condition) {							\
      index = running_index;						\
      running_index += number_of_indices;				\
    }									\
  } while(0);								\

/** parameters related to the precision of the code and to the method of calculation */

/** 
 * List of methods for stopping the transfer function computation 
 * at a given k for each l (saves lots of time). 
 */

enum transfer_cutting {
  tc_none, /**< no transfer cut: for given l, compute transfer functions over full k range (long and usually useless) */ 
  tc_osc, /**< transfer cut with oscillation method: for given l, compute transfer functions until k_max such that oscillations of \f$ \Delta_l(k) \f$ are tiny relatively to largest oscillation */ 
  tc_cl, /**< transfer cut with Cl variation method: for given l, compute transfer functions until k_max such that C_l's variation is tiny (C_l being computed approximately and with flat spectrum)  */ 
  tc_env /**< under development */
};

/**
 * list of evolver types for integrating perturbations over time
 */
enum evolver_type {
  rk, /* Runge-Kutta integrator */
  ndf15 /* stiff integrator */
};

/** 
 * List of ways in which matter power spectrum P(k) can be defined.
 * The standard definition is the first one (delta_m_squared) but
 * alternative definitions can be usfeul in some projects.
 * 
 */
enum pk_def {
  delta_m_squared, /**< normal definition (delta_m includes all non-relativistic species at late times) */
  delta_tot_squared, /**< delta_tot includes all species contributions to (delta rho), and only non-relativistic contributions to rho */
  delta_bc_squared, /**< delta_bc includes contribution of baryons and cdm only to (delta rho) and to rho */
  delta_tot_from_poisson_squared /**< use delta_tot inferred from gravitational potential through Poisson equation */
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
   * default step d tau in background integration, in units of 
   * conformal Hubble time (\f$ d tau \f$ = back_integration_stepsize / aH )
   */
  double back_integration_stepsize; 

  /**
   * parameter controlling precision of background integration
   */
  double tol_background_integration;


  /**
   * parameter controlling how deep inside radiation domination must the
   * initial time be chosen
   */
  double tol_initial_Omega_r;

  /**
   * parameter controlling relative precision of ncdm mass for given
   * ncdm current density
   */
  double tol_M_ncdm;

  /**
   * parameter controlling relative precision of integrals over ncdm
   * phase-space distribution during perturbation calculation
   */
  double tol_ncdm;

  /**
   * parameter controlling relative precision of integrals over ncdm
   * phase-space distribution during background evolution
   */
  double tol_ncdm_bg;

  /**
   * parameter controlling how relativistic must non-cold relics be at
   * initial time
   */
  double tol_ncdm_initial_w;

  //@}

  /** @name - parameters related to the thermodynamics */

  //@{

  /** - for bbn */

  FileName sBBN_file;

  /** - for recombination */

  /* initial and final redshifts in recfast */

  double recfast_z_initial;      /**< initial redshift in recfast */

  /* parameters governing precision of integration */
  
  int recfast_Nz0;               /**< number of integration steps */
  double tol_thermo_integration; /**< precision of each integration step */

  /* He fudge parameters from recfast 1.4 */

  int recfast_Heswitch;           /**< recfast 1.4 parameter */
  double recfast_fudge_He;        /**< recfast 1.4 parameter */

  /* H  fudge parameters from recfast 1.5 (Gaussian fits for extra H physics by Adam Moss) */

  int recfast_Hswitch;            /**< recfast 1.5 switching parameter */
  double recfast_fudge_H;         /**< H fudge factor when recfast_Hswitch set to false (v1.4 fudging) */
  double recfast_delta_fudge_H;   /**< correction to H fudge factor in v1.5 */
  double recfast_AGauss1;         /**< Amplitude of 1st Gaussian */
  double recfast_AGauss2;         /**< Amplitude of 2nd Gaussian */
  double recfast_zGauss1;         /**< ln(1+z) of 1st Gaussian */
  double recfast_zGauss2;         /**< ln(1+z) of 2nd Gaussian */
  double recfast_wGauss1;         /**< Width of 1st Gaussian */
  double recfast_wGauss2;         /**< Width of 2nd Gaussian */

  /* triggers for switching approximations; ranges for doing it smoothly */

  double recfast_z_He_1;              /**< down to which redshift Helium fully ionized */
  double recfast_delta_z_He_1;        /**< z range over which transition is smoothed */

  double recfast_z_He_2;              /**< down to which redshift first Helium recombination 
					   not complete */
  double recfast_delta_z_He_2;        /**< z range over which transition is smoothed */

  double recfast_z_He_3;              /**< down to which redshift Helium singly ionized */
  double recfast_delta_z_He_3;        /**< z range over which transition is smoothed */

  double recfast_x_He0_trigger;       /**< below which Helium ionization fraction start using 
                                           full equation for Helium */
  double recfast_x_He0_trigger2;      /**< a second threshold used in derivative routine */
  double recfast_x_He0_trigger_delta; /**< x_He range over which transition is smoothed */

  double recfast_x_H0_trigger;        /**< below which Helium ionization fraction start using 
                                           full equation for Helium */
  double recfast_x_H0_trigger2;       /**< a second threshold used in derivative routine */
  double recfast_x_H0_trigger_delta;  /**< x_H range over which transition is smoothed */

  double recfast_H_frac;              /**< governs time at which full equation of evolution 
					   for Tmat is used */

  FileName hyrec_Alpha_inf_file;
  FileName hyrec_R_inf_file;
  FileName hyrec_two_photon_tables_file;

  /** - for reionization */

  double reionization_z_start_max; /**< maximum redshift at which reionization should start. If not, return an error. */
  double reionization_sampling; /**< control stepsize in z during reionization */
  double reionization_optical_depth_tol; /**< fractional error on optical_depth */
  double reionization_start_factor; /**< parameter for CAMB-like parametrization */
   
  /** - general */

  int thermo_rate_smoothing_radius; /**< plays a minor (almost aesthetic) role in the definition of the variation rate of thermodynamical quantities */

  //@}

  /** @name - parameters related to the perturbation */

  //@{

  enum evolver_type evolver; /* which type of evolver for integrating perturbations (Runge-Kutta? Stiff?...) */

  enum pk_def pk_definition;

  double k_scalar_min_tau0; /**< number defining k_min for the computation of scalar Cl's and P(k)'s (dimensionless): (k_min tau_0), usually chosen much smaller than one */

  double k_scalar_max_tau0_over_l_max; /**< number defining k_max for the computation of scalar Cl's (dimensionless): (k_max tau_0)/l_max, usually chosen around two */

  double k_scalar_step_sub; /**< step in k space, in units of one period of acoustic oscillation at decoupling, for scales inside sound horizon at decoupling */
  double k_scalar_step_super; /**< step in k space, in units of one period of acoustic oscillation at decoupling, for scales above sound horizon at decoupling */  
  double k_scalar_step_transition; /**< dimensionless number regulating the transition from 'sub' steps to 'super' steps. Decrease for more precision. */

  double k_scalar_k_per_decade_for_pk; /**< if values needed between kmax inferred from k_scalar_oscillations and k_scalar_kmax_for_pk, this gives the number of k per decade outside the BAO region*/

  double k_scalar_k_per_decade_for_bao; /**< if values needed between kmax inferred from k_scalar_oscillations and k_scalar_kmax_for_pk, this gives the number of k per decade inside the BAO region (for finer sampling)*/

  double k_scalar_bao_center; /**< in ln(k) space, the central value of the BAO region where sampling is finer is defined as k_rec times this number (recommended: 3, i.e. finest sampling near 3rd BAO peak) */

  double k_scalar_bao_width; /**< in ln(k) space, width of the BAO region where sampling is finer: this number gives roughly the number of BAO oscillations well resolved on both sides of the central value (recommended: 4, i.e. finest sampling from before first up to 3+4=7th peak) */

  double k_tensor_min_tau0; /**< number defining k_min for the computation of tensor Cl's (dimensionless): (k_min tau_0), usually chosen much smaller than one */

  double k_tensor_max_tau0_over_l_max; /**< number defining k_max for the computation of tensor Cl's (dimensionless): (k_max tau_0)/l_max, usually chosen around two */

  double k_tensor_step_sub; /**< step in k space, in units of one period of oscillation at decoupling, for scales inside horizon at decoupling (tensor modes) */
  double k_tensor_step_super; /**< step in k space, in units of one period of oscillation at decoupling, for scales above horizon at decoupling (tensor modes) */  
  double k_tensor_step_transition; /**< dimensionless number regulaing the transition fro _sub step to _super step. Decrease for more precision. (tensor modes) */

  double start_small_k_at_tau_c_over_tau_h; /**< largest wavelengths start being sampled when universe is sufficiently opaque. This is quantified in terms of the ratio of thermo to hubble time scales, \f$ \tau_c/\tau_H \f$. Start when start_largek_at_tau_c_over_tau_h equals this ratio. Decrease this value to start integrating the wavenumbers earlier in time. */

  double start_large_k_at_tau_h_over_tau_k;  /**< largest wavelengths start being sampled when mode is sufficiently outside Hibble scale. This is quantified in terms of the ratio of hubble time scale to wavenumber time scale, \f$ \tau_h/\tau_k \f$ wich is roughly equal to (k*tau). Start when this ratio equals start_large_k_at_tau_k_over_tau_h. Decrease this value to start integrating the wavenumbers earlier in time. */

  /**
   * when to switch off tight-coupling approximation: first condition:
   * \f$ \tau_c/\tau_H \f$ > tight_coupling_trigger_tau_c_over_tau_h.
   * Decrease this value to switch off earlier in time.  If this
   * number is larger than start_sources_at_tau_c_over_tau_h, the code
   * returns an error, because the source computation requires
   * tight-coupling to be switched off.
   */
  double tight_coupling_trigger_tau_c_over_tau_h;

  /**
   * when to switch off tight-coupling approximation:
   * second condition: \f$ \tau_c/\tau_k \equiv k \tau_c \f$ <
   * tight_coupling_trigger_tau_c_over_tau_k.
   * Decrease this value to switch off earlier in time.
   */
  double tight_coupling_trigger_tau_c_over_tau_k;

  double start_sources_at_tau_c_over_tau_h; /**< sources start being sampled when universe is sufficiently opaque. This is quantified in terms of the ratio of thermo to hubble time scales, \f$ \tau_c/\tau_H \f$. Start when start_sources_at_tau_c_over_tau_h equals this ratio. Decrease this value to start sampling the sources earlier in time. */

  int tight_coupling_approximation;

  int l_max_g;     /**< number of momenta in Boltzmann hierarchy for photon temperature (scalar), at least 4 */
  int l_max_pol_g; /**< number of momenta in Boltzmann hierarchy for photon polarisation (scalar), at least 4 */
  int l_max_ur;   /**< number of momenta in Boltzmann hierarchy for relativistic neutrino/relics (scalar), at least 4 */
  int l_max_ncdm;   /**< number of momenta in Boltzmann hierarchy for relativistic neutrino/relics (scalar), at least 4 */
  int l_max_g_ten;     /**< number of momenta in Boltzmann hierarchy for photon temperature (tensor), at least 4 */
  int l_max_pol_g_ten; /**< number of momenta in Boltzmann hierarchy for photon polarisation (tensor), at least 4 */

  double curvature_ini;     /**< initial condition for curvature for adiabatic */
  double entropy_ini; /**< initial condition for entropy perturbation for isocurvature */ 
  double gw_ini;      /**< initial condition for tensor metric perturbation h */

  /** 
   * default step \f$ d \tau \f$ in perturbation integration, in units of the timescale involved in the equations (usally, the min of \f$ 1/k \f$, \f$ 1/aH \f$, \f$ 1/\dot{\kappa} \f$) 
   */
  double perturb_integration_stepsize;

  /** 
   * default step \f$ d \tau \f$ for sampling the source function, in units of the timescale involved in the sources: \f$ (\dot{\kappa}- \ddot{\kappa}/\dot{\kappa})^{-1} \f$
   */
  double perturb_sampling_stepsize;

  /** 
   * control parameter for the precision of the perturbation integration 
   */
  double tol_perturb_integration;

  /**
   * precision with which the code should determine (by bisection) the
   * times at which sources start being sampled, and at which
   * approximations must be switched on/off (units of Mpc)
   */
  double tol_tau_approx;

  /**
   * method for switching off photon perturbations
   */
  int radiation_streaming_approximation;

  /**
   * when to switch off photon perturbations, ie when to switch
   * on photon free-streaming approximation (keep density and thtau, set
   * shear and higher momenta to zero):
   * first condition: \f$ k \tau \f$ > radiation_streaming_trigger_tau_h_over_tau_k
   */
  double radiation_streaming_trigger_tau_over_tau_k;

  /**
   * when to switch off photon perturbations, ie when to switch
   * on photon free-streaming approximation (keep density and theta, set
   * shear and higher momenta to zero):
   * second condition: 
   */ 
  double radiation_streaming_trigger_tau_c_over_tau;

  int ur_fluid_approximation;

  /**
   * when to switch off ur (massless neutrinos / ultra-relativistic
   * relics) fluid approximation 
   */
  double ur_fluid_trigger_tau_over_tau_k;

  int ncdm_fluid_approximation;

  /**
   * when to switch off ncdm (massive neutrinos / non-cold
   * relics) fluid approximation 
   */
  double ncdm_fluid_trigger_tau_over_tau_k;

  //@}

  /** @name - parameters related to Bessel functions */

  //@{

  int l_linstep; /**< factor for logarithmic spacing of values of l over which bessel and transfer functions are sampled */
  double l_logstep; /**< maximum spacing of values of l over which Bessel and transfer functions are sampled (so, spacing becomes linear instead of logarithmic at some point) */

  double bessel_x_step; /**< step dx for sampling Bessel functions \f$ j_l(x) \f$ */
  double bessel_j_cut; /**< value of \f$ j_l \f$ below which it is approximated by zero (in the region \f$ x \ll l \f$) */
  double bessel_tol_x_min;  /**< precision with which x_min such that j_l(x_min)=j_cut is found (order of magnitude set by k_min) */
  FileName bessel_file_name; /**< name of file where Bessel functions will evnetually be written or read */

  //@}

  /** @name - parameters related to the primordial spectra */

  //@{

  double k_per_decade_primordial; /**< logarithmic sampling for primordial spectra (number of points per decade in k space) */

  double primordial_inflation_ratio_min;
  double primordial_inflation_ratio_max;
  int primordial_inflation_phi_ini_maxit;
  double primordial_inflation_pt_stepsize;
  double primordial_inflation_bg_stepsize;
  double primordial_inflation_tol_integration;
  double primordial_inflation_attractor_precision_pivot;
  double primordial_inflation_attractor_precision_initial;
  int primordial_inflation_attractor_maxit;
  double primordial_inflation_jump_initial;
  double primordial_inflation_tol_curvature;

  //@}

  /** @name - parameters related to the transfer function */

  //@{

  double k_step_trans_scalars; /**< upper bound on linear sampling step in k space, in units of 2pi/tau0 (where tau0 is the conformal time today) */
  double k_step_trans_tensors; /**< upper bound on linear sampling step in k space, in units of 2pi/tau0 (where tau0 is the conformal time today) */

  enum transfer_cutting transfer_cut; /**< flag telling how to cut automatically the transfer function computation at a given \f$ k_{max} \f$ value */

  double transfer_cut_threshold_osc; /**< threshold used for cutting the transfer function computation at a given \f$ k_{max} \f$ value, if transfer_cut = _TC_OSC_ (oscillation method: for given l, compute transfer functions until k_max such that oscillations of \f$ \Delta_l(k) \f$ are tiny relatively to largest oscillation) */

  double transfer_cut_threshold_cl; /**< threshold used for cutting the transfer function computation at a given \f$ k_{max} \f$ value, if transfer_cut = _TC_CL_ (Cl variation method: for given l, compute transfer functions until k_max such that C_l's variation is tiny, C_l being computed approximately and with flat spectrum) */

  /** when to use the Limber approximation for project gravitational potential cl's */
  double l_switch_limber;

  /** when to use the Limber approximation for density cl's (relative to central redshift of each bin) */
  double l_switch_limber_for_cl_density_over_z;

  /** in sigma units, where to cut gaussian selection functions */
  double selection_cut_at_sigma;
  
  /** controls sampling of integral over time when selection functions vary quicker than Bessel functions. Increase for better sampling. */
  double selection_sampling; 

  /** controls sampling of integral over time when selection functions vary slower than Bessel functions. Increase for better sampling */
  double selection_sampling_bessel;

  /** controls how smooth are the edge of top-hat window function (<<1 for very sharp, 0.1 for sharp) */
  double selection_tophat_edge;

  //@}

  /** @name - parameters related to spectra */

  //@{

  /* nothing */

  //@}

  /** @name - parameters related to non-linear computations */

  //@{

  /** parameters relevant for HALOFIT computation */

  double halofit_dz; /* spacing in redshift space defining values of z
			at which HALOFIT will be used. Intermediate
			values will be obtained by
			interpolation. Decrease for more precise
			interpolations, at the expense of increasing
			time spent in nonlinear_init() */ 

  double halofit_min_k_nonlinear; /* value of k in 1/Mpc above
				     which non-linear corrections will
				     be computed */

  double halofit_sigma_precision; /* a smaller value will lead to a
				      more precise halofit result at
				      the highest requested redshift,
				      at the expense of requiring a
				      larger k_max */

  /** parameters relevant for TRG computation */

  int double_escape;      /* number of points to drop at every
			     half-step of the computation (double
			     espace mechanism) */
  double z_ini;           /* starting redshift for computing the
			     non-linear evolution */
  int eta_size;           /* number of steps in the time-variable eta */
  double k_L;             /* scale above which we consider linear
			     theory to be true at any time of the
			     computation. To be more explicit, for k's
			     smaller than this one, the A's function
			     are forced to 0. */
  double k_min;           /* lower bound for computing non-linear
			     matter power spectrum */
  double logstepx_min;    /* precision parameter for computation of
			     the As functions. */
  double logstepk1;       /* various parameters used in the definition
			     of steps in k space in the non-linear
			     calculation */
  double logstepk2;
  double logstepk3;
  double logstepk4;
  double logstepk5;
  double logstepk6;
  double logstepk7;
  double logstepk8;
  double k_growth_factor; /* at which scale k in 1/Mpc units do we
			     want to evaluate the LmabdaCDM growth
			     factor? */

  double k_scalar_max_for_pk_nl; /* max value of k in h/Mpc to be used
				    in the trg module */

  //@}

  /** @name - parameters related to lensing */

  //@{

  int accurate_lensing; /**< switch between Gauss-Legendre quadrature integration and simple quadrature on a subdomain of angles */
  int num_mu_minus_lmax; /**< difference between num_mu and l_max, increase for more precision */
  int delta_l_max; /**<difference between l_max in unlensed and lensed spectra */
  double tol_gauss_legendre; /**< tolerance with which quadrature points are found: must be very small for an accurate integration (if not entered manually, set automatically to match machine precision) */
  //@}

  /** @name - general precision parameters */

  //@{

  double smallest_allowed_variation; /**< machine-dependent, assigned automatically by the code */

  //@}

  /** @name - zone for writing error messages */

  //@{

  ErrorMsg error_message;

  //@}

};



#endif
