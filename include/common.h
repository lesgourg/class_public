/** @file common.h Generic libraries, parameters and functions used in the whole code. */

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "float.h"
#ifdef _OPENMP
#include "omp.h"
#endif

#ifndef __COMMON__
#define __COMMON__

#define _TRUE_ 1 /**< integer associated to true statement */
#define _FALSE_ 0 /**< integer associated to false statement */

#define _SUCCESS_ 0 /**< integer returned after successfull call of a function */
#define _FAILURE_ 1 /**< integer returned after failure in a function */

#define _ERRORMSGSIZE_ 2048 /**< generic error messages are cut beyond this number of characters */
typedef char ErrorMsg[_ERRORMSGSIZE_]; /**< Generic error messages (there is such a field in each structure) */

#define _FILENAMESIZE_ 40 /**< size of the string read in each line of the file (extra characters not taken into account) */
typedef char FileName[_FILENAMESIZE_];

#define _PI_ 3.1415926535897932384626433832795e0 /**< The number pi */

#define _MAX_IT_ 10000/**< default maximum number of iterations in conditional loops (to avoid infinite loops) */

#define _QUADRATURE_MAX_ 150 /**< maximum allowed number of abssices in quadrature integral estimation */

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
#define class_stop(error_message_output)				\
  do {									\
    sprintf(error_message_output,					\
	    "%s(L:%d) : Stop here as requested",			\
	    __func__,__LINE__);						\
    return _FAILURE_;							\
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

/** parameters related to the precision of the code and to the method of calculation */

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
 * list of method for integrating transfer functions
 */
enum transfer_integration {spline,trapezoidal};

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
   * (a) when g becomes smaller than
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

  double l_max_over_k_max_scalars; /**< number defining k_max (in units of Mpc^-1) for scalars: k_max = l_max times this number (order of magnitude given by few Hubble radius in Mpc, i.e. few times 3000) */

  double k_scalar_step_sub; /**< step in k space, in units of one period of acoustic oscillation at decoupling, for scales inside sound horizon at decoupling */
  double k_scalar_step_super; /**< step in k space, in units of one period of acoustic oscillation at decoupling, for scales above sound horizon at decoupling */  
  double k_scalar_step_transition; /**< dimensionless number regulaing the transition fro _sub step to _super step. Decrease for more precision. */

  double k_scalar_k_per_decade_for_pk; /**< if values needed between kmax inferred from k_scalar_oscillations and k_scalar_kmax_for_pk, this gives the number of k per decade */


  double k_tensor_min; /**< first mode k_min in units of Hubble radius today (tensor modes) */

  double l_max_over_k_max_tensors; /**< number defining k_max (in units of Mpc^-1) for tensors: k_max = l_max times this number (order of magnitude given by few Hubble radius in Mpc, i.e. few times 3000) */

  double k_tensor_step_sub; /**< step in k space, in units of one period of oscillation at decoupling, for scales inside horizon at decoupling (tensor modes) */
  double k_tensor_step_super; /**< step in k space, in units of one period of oscillation at decoupling, for scales above horizon at decoupling (tensor modes) */  
  double k_tensor_step_transition; /**< dimensionless number regulaing the transition fro _sub step to _super step. Decrease for more precision. (tensor modes) */

  double start_small_k_at_eta_g_over_eta_h; /**< largest wavelengths start being sampled when universe is sufficiently opaque. This is quantified in terms of the ratio of thermo to hubble time scales, \f$ \eta_g/\eta_H \f$. Start when start_largek_at_eta_g_over_eta_h equals this ratio. Decrease this value to start integrating the wavenumbers earlier in time. */

  double start_large_k_at_eta_g_over_eta_k; /**< largest wavelengths start being sampled when universe is sufficiently opaque. This is quantified in terms of the ratio of thermo to hubble time scales, \f$ \eta_g/\eta_k \f$. Start when start_largek_at_eta_g_over_eta_k equals this ratio. Decrease this value to start integrating the wavenumbers earlier in time. */

  /**
   * when to switch off tight-coupling approximation: first condition:
   * \f$ \eta_g/\eta_H \f$ < tight_coupling_trigger_eta_g_over_eta_h.
   * Decrease this value to switch off earlier in time.  If this
   * number is larger than start_sources_at_eta_g_over_eta_h, the code
   * returns an error, because the source computation requires
   * tight-coupling to be switched off.
   */
  double tight_coupling_trigger_eta_g_over_eta_h;

  /**
   * when to switch off tight-coupling approximation:
   * second condition: \f$ \eta_g/\eta_k \equiv k \eta_g \f$ <
   * tight_coupling_trigger_eta_g_over_eta_k.
   * Decrease this value to switch off earlier in time.
   */
  double tight_coupling_trigger_eta_g_over_eta_k;

  double start_sources_at_eta_g_over_eta_h; /**< sources start being sampled when universe is sufficiently opaque. This is quantified in terms of the ratio of thermo to hubble time scales, \f$ \eta_g/\eta_H \f$. Start when start_sources_at_eta_g_over_eta_h equals this ratio. Decrease this value to start sampling the sources earlier in time. */

  int tight_coupling_approximation;

  double k_eta_max; /**< maximum value of \f$ k \times \eta \f$ for which the sources should be computed in a pure CMB run*/

  int l_max_g;     /**< number of momenta in Boltzmann hierarchy for photon temperature (scalar), at least 4 */
  int l_max_pol_g; /**< number of momenta in Boltzmann hierarchy for photon polarisation (scalar), at least 4 */
  int l_max_nur;   /**< number of momenta in Boltzmann hierarchy for relativistic neutrino/relics (scalar), at least 4 */
  int l_max_g_ten;     /**< number of momenta in Boltzmann hierarchy for photon temperature (tensor), at least 4 */
  int l_max_pol_g_ten; /**< number of momenta in Boltzmann hierarchy for photon polarisation (tensor), at least 4 */

  double curvature_ini;     /**< initial condition for curvature for adiabatic */
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
   * precision with which the code should determine (by bisection) the
   * times at which sources start being sampled, and at which
   * approximations must be switched on/off (units of Mpc)
   */
  double tol_eta_approx;
 
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

  int l_linstep; /**< factor for logarithmic spacing of values of l over which bessel and transfer functions are sampled */
  double l_logstep; /**< maximum spacing of values of l over which Bessel and transfer functions are sampled (so, spacing becomes linear instead of logarithmic at some point) */

  double bessel_x_step; /**< step dx for sampling Bessel functions \f$ j_l(x) \f$ */
  double bessel_j_cut; /**< value of \f$ j_l \f$ below which it is approximated by zero (in the region \f$ x \ll l \f$) */
  double bessel_tol_x_min;  /**< precision with which x_min such that j_l(x_min)=j_cut is found (order of magnitude set by k_min) */
  double bessel_x_max_over_l_max; /**< x_max is defined as this parameter times largest l in Cls computation */ 

  //@}

  /** @name - parameters related to the primordial spectra */

  //@{

  double k_per_decade_primordial; /**< logarithmic sampling for primordial spectra (number of points per decade in k space) */

  //@}

  /** @name - parameters related to the transfer function */

  //@{

  double k_step_trans_scalars; /**< sampling step in k space, in units of \f$ 2\pi/(\eta_0-\eta_{rec}) \f$, which is the typical period of oscillations of \f$ \Delta_l(k) \f$ */
  double k_step_trans_tensors; /**< sampling step in k space, in units of \f$ 2\pi/(\eta_0-\eta_{rec}) \f$, which is the typical period of oscillations of \f$ \Delta_l(k) \f$ */

/*   int k_oversampling_scalars; */
/*   int k_oversampling_tensors; */

  enum transfer_cutting transfer_cut; /**< flag telling how to cut automatically the transfer function computation at a given \f$ k_{max} \f$ value */

  double transfer_cut_threshold_osc; /**< threshold used for cutting the transfer function computation at a given \f$ k_{max} \f$ value, if transfer_cut = _TC_OSC_ (oscillation method: for given l, compute transfer functions until k_max such that oscillations of \f$ \Delta_l(k) \f$ are tiny relatively to largest oscillation) */

  double transfer_cut_threshold_cl; /**< threshold used for cutting the transfer function computation at a given \f$ k_{max} \f$ value, if transfer_cut = _TC_CL_ (Cl variation method: for given l, compute transfer functions until k_max such that C_l's variation is tiny, C_l being computed approximately and with flat spectrum) */

  int l_switch_limber;

  double smallest_allowed_variation; /**< machine-dependent, assigned automatically by the code */

  enum transfer_integration transfer_integrate;

  //@}

  /** @name - parameters related to lensing */

  //@{

  int num_mu_minus_lmax; /**< difference between num_mu and l_max, increase for more precision */
  int delta_l_max; /**<difference between l_max in unlensed and lensed spectra */

  //@}

  /** @name - zone for writing error messages */

  //@{

  ErrorMsg error_message;

  //@}

};



#endif
