/** @file common.h Generic libraries, parameters and functions used in the whole code. */

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "float.h"
#include "svnversion.h"
#include <stdarg.h>

#ifdef _OPENMP
#include "omp.h"
#endif

#ifndef __COMMON__
#define __COMMON__

#define _VERSION_ "v2.7.2"
/* @cond INCLUDE_WITH_DOXYGEN */

#define _TRUE_ 1 /**< integer associated to true statement */
#define _FALSE_ 0 /**< integer associated to false statement */

#define _SUCCESS_ 0 /**< integer returned after successful call of a function */
#define _FAILURE_ 1 /**< integer returned after failure in a function */

#define _ERRORMSGSIZE_ 2048 /**< generic error messages are cut beyond this number of characters */
typedef char ErrorMsg[_ERRORMSGSIZE_]; /**< Generic error messages (there is such a field in each structure) */

#define _FILENAMESIZE_ 256 /**< size of the string read in each line of the file (extra characters not taken into account) */
typedef char FileName[_FILENAMESIZE_];

#define _PI_ 3.1415926535897932384626433832795e0 /**< The number pi */

#define _PIHALF_ 1.57079632679489661923132169164e0 /**< pi divided by 2 */

#define _TWOPI_ 6.283185307179586476925286766559e0 /**< 2 times pi */

#define _SQRT2_ 1.41421356237309504880168872421e0 /** < square root of 2. */

#define _SQRT6_ 2.4494897427831780981972840747059e0 /**< square root of 6. */

#define _SQRT_PI_ 1.77245385090551602729816748334e0 /**< square root of pi. */

#define _MAX_IT_ 10000/**< default maximum number of iterations in conditional loops (to avoid infinite loops) */

#define _QUADRATURE_MAX_ 250 /**< maximum allowed number of abssices in quadrature integral estimation */

#define _QUADRATURE_MAX_BG_ 800 /**< maximum allowed number of abssices in quadrature integral estimation */

#define _TOLVAR_ 100. /**< The minimum allowed variation is the machine precision times this number */

#define _HUGE_ 1.e99

#define _OUTPUTPRECISION_ 12 /**< Number of significant digits in some output files */

#define _COLUMNWIDTH_ 24 /**< Must be at least _OUTPUTPRECISION_+8 for guaranteed fixed width columns */

#define _MAXTITLESTRINGLENGTH_ 8000 /**< Maximum number of characters in title strings */

#define _DELIMITER_ "\t" /**< character used for delimiting titles in the title strings */



#ifndef __CLASSDIR__
#define __CLASSDIR__ "." /**< The directory of CLASS. This is set to the absolute path to the CLASS directory so this is just a failsafe. */
#endif

#define MIN(a,b) (((a)<(b)) ? (a) : (b) ) /**< the usual "min" function */
#define MAX(a,b) (((a)<(b)) ? (b) : (a) ) /**< the usual "max" function */
#define SIGN(a) (((a)>0) ? 1. : -1. )
#define NRSIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define index_symmetric_matrix(i1,i2,N) (((i1)<=(i2)) ? (i2+N*i1-(i1*(i1+1))/2) : (i1+N*i2-(i2*(i2+1))/2)) /**< assigns an index from 0 to [N(N+1)/2-1] to the coefficients M_{i1,i2} of an N*N symmetric matrix; useful for converting a symmetric matrix to a vector, without losing or double-counting any information */
/* @endcond */
// needed because of weird openmp bug on macosx lion...
void class_protect_sprintf(char* dest, char* tpl,...);
void class_protect_fprintf(FILE* dest, char* tpl,...);
void* class_protect_memcpy(void* dest, void* from, size_t sz);

int get_number_of_titles(char * titlestring);

#define class_build_error_string(dest,tmpl,...) {                                                                \
  ErrorMsg FMsg;                                                                                                 \
  class_protect_sprintf(FMsg,tmpl,__VA_ARGS__);                                                                  \
  class_protect_sprintf(dest,"%s(L:%d) :%s",__func__,__LINE__,FMsg);                                             \
}

// Error reporting macros

// Call
#define class_call_message(err_out,extra,err_mess)   \
  class_build_error_string(err_out,"error in %s;\n=>%s",extra,err_mess);

/* macro for calling function and returning error if it failed */
#define class_call_except(function, error_message_from_function, error_message_output,list_of_commands) {        \
  if (function == _FAILURE_) {                                                                                   \
    class_call_message(error_message_output,#function,error_message_from_function);                              \
    list_of_commands;                                                                                            \
    return _FAILURE_;                                                                                            \
  }                                                                                                              \
}

/* macro for trying to call function */
#define class_call_try(function, error_message_from_function, error_message_output,list_of_commands) { \
  if (function == _FAILURE_) {                                                                                   \
    class_call_message(error_message_output,#function,error_message_from_function);                              \
    list_of_commands;                                                                                            \
  }                                                                                                              \
}

/* macro for calling function and returning error if it failed */
#define class_call(function, error_message_from_function, error_message_output)                                  \
  class_call_except(function, error_message_from_function,error_message_output,)

/* same in parallel region */
#define class_call_parallel(function, error_message_from_function, error_message_output) {                       \
  if (abort == _FALSE_) {                                                                                        \
    if (function == _FAILURE_) {                                                                                 \
      class_call_message(error_message_output,#function,error_message_from_function);                            \
      abort=_TRUE_;                                                                                              \
    }                                                                                                            \
  }                                                                                                              \
}




// Alloc
#define class_alloc_message(err_out,extra,sz)                                                                    \
  class_build_error_string(err_out,"could not allocate %s with size %d",extra,sz);

/* macro for allocating memory and returning error if it failed */
#define class_alloc(pointer, size, error_message_output)  {                                                      \
  pointer=malloc(size);                                                                                          \
  if (pointer == NULL) {                                                                                         \
    int size_int;                                                                                                \
    size_int = size;                                                                                             \
    class_alloc_message(error_message_output,#pointer, size_int);                                                \
    return _FAILURE_;                                                                                            \
  }                                                                                                              \
}

/* same inside parallel structure */
#define class_alloc_parallel(pointer, size, error_message_output)  {                                             \
  pointer=NULL;                                                                                                  \
  if (abort == _FALSE_) {                                                                                        \
    pointer=malloc(size);                                                                                        \
    if (pointer == NULL) {                                                                                       \
      int size_int;                                                                                              \
      size_int = size;                                                                                           \
      class_alloc_message(error_message_output,#pointer, size_int);                                              \
      abort=_TRUE_;                                                                                              \
    }                                                                                                            \
  }                                                                                                              \
}

/* macro for allocating memory, initializing it with zeros/ and returning error if it failed */
#define class_calloc(pointer, init,size, error_message_output)  {                                                \
  pointer=calloc(init,size);                                                                                     \
  if (pointer == NULL) {                                                                                         \
    int size_int;                                                                                                \
    size_int = size;                                                                                             \
    class_alloc_message(error_message_output,#pointer, size_int);                                                \
    return _FAILURE_;                                                                                            \
  }                                                                                                              \
}

/* macro for re-allocating memory, returning error if it failed */
#define class_realloc(pointer, newname, size, error_message_output)  {                                          \
    pointer=realloc(newname,size);                                                                               \
  if (pointer == NULL) {                                                                                         \
    int size_int;                                                                                                \
    size_int = size;                                                                                             \
    class_alloc_message(error_message_output,#pointer, size_int);                                                \
    return _FAILURE_;                                                                                            \
  }                                                                                                              \
}

// Testing

#define class_test_message(err_out,extra,args...) {                                                              \
  ErrorMsg Optional_arguments;                                                                                   \
  class_protect_sprintf(Optional_arguments,args);                                                                \
  class_build_error_string(err_out,"condition (%s) is true; %s",extra,Optional_arguments);                       \
}

/* macro for testing condition and returning error if condition is true;
   args is a variable list of optional arguments, e.g.: args="x=%d",x
   args cannot be empty, if there is nothing to pass use args="" */
#define class_test_except(condition, error_message_output,list_of_commands, args...) {                           \
  if (condition) {                                                                                               \
    class_test_message(error_message_output,#condition, args);                                                   \
    list_of_commands;                                                                                            \
    return _FAILURE_;                                                                                            \
  }                                                                                                              \
}

#define class_test(condition, error_message_output, args...) {                                                   \
  if (condition) {                                                                                               \
    class_test_message(error_message_output,#condition, args);                                                   \
    return _FAILURE_;                                                                                            \
  }                                                                                                              \
}

#define class_test_parallel(condition, error_message_output, args...) {                                          \
  if (abort == _FALSE_) {                                                                                        \
    if (condition) {                                                                                             \
      class_test_message(error_message_output,#condition, args);                                                 \
      abort=_TRUE_;                                                                                              \
    }                                                                                                            \
  }                                                                     \
}

/* macro for returning error message;
   args is a variable list of optional arguments, e.g.: args="x=%d",x
   args cannot be empty, if there is nothing to pass use args="" */
#define class_stop(error_message_output,args...) {                                                               \
  ErrorMsg Optional_arguments;                                                                                   \
  class_protect_sprintf(Optional_arguments,args);                                                                \
  class_build_error_string(error_message_output,"error; %s",Optional_arguments);                                 \
  return _FAILURE_;                                                                                              \
}

// IO
/* macro for opening file and returning error if it failed */
#define class_open(pointer, filename,	mode, error_output) {                                                      \
  pointer=fopen(filename,mode);                                                                                  \
  if (pointer == NULL) {                                                                                         \
    class_build_error_string(error_output,"could not open %s with name %s and mode %s",#pointer,filename,#mode); \
    return _FAILURE_;                                                                                            \
  }                                                                                                              \
}

/* macro for defining indices (usually one, sometimes a block) */
#define class_define_index(index,                                       \
                           condition,                                   \
                           running_index,                               \
                           number_of_indices) {                         \
    if (condition) {                                                    \
      index = running_index;                                            \
      running_index += number_of_indices;                               \
    }                                                                   \
  }

/* macros for writing formatted output */
#define class_fprintf_double(file,                                      \
                             output,                                    \
                             condition){                                \
    if (condition == _TRUE_)                                            \
      fprintf(file,"%*.*e ",_COLUMNWIDTH_,_OUTPUTPRECISION_,output);    \
  }

#define class_fprintf_double_or_default(file,                           \
                                        output,                         \
                                        condition,                      \
                                        defaultvalue){                  \
    if (condition == _TRUE_)                                            \
      fprintf(file,"%*.*e ",_COLUMNWIDTH_,_OUTPUTPRECISION_,output);    \
    else                                                                \
      fprintf(file,"%*.*e ",_COLUMNWIDTH_,_OUTPUTPRECISION_,defaultvalue);    \
}

#define class_fprintf_int(file,                                         \
                          output,                                       \
                          condition){                                   \
    if (condition == _TRUE_)                                            \
      fprintf(file,"%*d%*s ",                                           \
              MAX(0,_COLUMNWIDTH_-_OUTPUTPRECISION_-5),                 \
              output, _OUTPUTPRECISION_+5," ");                          \
  }

#define class_fprintf_columntitle(file,                                 \
                                  title,                                \
                                  condition,                            \
                                  colnum){                              \
    if (condition == _TRUE_)                                            \
      fprintf(file,"%*s%2d:%-*s ",                                      \
              MAX(0,MIN(_COLUMNWIDTH_-_OUTPUTPRECISION_-6-3,_COLUMNWIDTH_-((int) strlen(title))-3)), \
              "",colnum++,_OUTPUTPRECISION_+6,title);                   \
  }

#define class_store_columntitle(titlestring,                            \
				title,					\
				condition){				\
    if (condition == _TRUE_){                                           \
      strcat(titlestring,title);                                        \
      strcat(titlestring,_DELIMITER_);                                  \
    }                                                                   \
  }
//,_MAXTITLESTRINGLENGTH_-strlen(titlestring)-1);

#define class_store_double(storage,					\
			   value,					\
			   condition,                                   \
                           dataindex){                                  \
    if (condition == _TRUE_)                                            \
      storage[dataindex++] = value;                                     \
  }

#define class_store_double_or_default(storage,                          \
                                      value,                            \
                                      condition,                        \
                                      dataindex,                        \
                                      defaultvalue){                    \
    if (condition == _TRUE_)                                            \
      storage[dataindex++] = value;                                     \
    else                                                                \
      storage[dataindex++] = defaultvalue;                              \
}

/** parameters related to the precision of the code and to the method of calculation */

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
 * alternative definitions can be useful in some projects.
 *
 */
enum pk_def {
  delta_m_squared, /**< normal definition (delta_m includes all non-relativistic species at late times) */
  delta_tot_squared, /**< delta_tot includes all species contributions to (delta rho), and only non-relativistic contributions to rho */
  delta_bc_squared, /**< delta_bc includes contribution of baryons and cdm only to (delta rho) and to rho */
  delta_tot_from_poisson_squared /**< use delta_tot inferred from gravitational potential through Poisson equation */
};
/**
 * Different ways to present output files
 */

enum file_format {class_format,camb_format};

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
   * conformal Hubble time (\f$ d \tau \f$ = back_integration_stepsize / aH )
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
   * phase-space distribution during perturbation calculation: value
   * to be applied in Newtonian gauge
   */
  double tol_ncdm_newtonian;

  /**
   * parameter controlling relative precision of integrals over ncdm
   * phase-space distribution during perturbation calculation: value
   * to be applied in synchronous gauge
   */
  double tol_ncdm_synchronous;

  /**
   * parameter controlling relative precision of integrals over ncdm
   * phase-space distribution during perturbation calculation: value
   * actually applied in chosen gauge
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

  /**
   * parameter controlling the initial scalar field in background functions
   */
  double safe_phi_scf;

  /**
   * parameter controlling precision with which tau_eq (conformal time
   * at radiation/matter equality) is found (units: Mpc)
   */
  double tol_tau_eq;

  //@}

  /** @name - parameters related to the thermodynamics */

  //@{

  /* - for bbn */
/* @cond INCLUDE_WITH_DOXYGEN */
  FileName sBBN_file;
/* @endcond */
  /* - for recombination */

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

  double recfast_z_He_2;              /**< down to which redshift first Helium recombination not complete */
  double recfast_delta_z_He_2;        /**< z range over which transition is smoothed */

  double recfast_z_He_3;              /**< down to which redshift Helium singly ionized */
  double recfast_delta_z_He_3;        /**< z range over which transition is smoothed */

  double recfast_x_He0_trigger;       /**< value below which recfast uses the full equation for Helium */
  double recfast_x_He0_trigger2;      /**< a second threshold used in derivative routine */
  double recfast_x_He0_trigger_delta; /**< x_He range over which transition is smoothed */

  double recfast_x_H0_trigger;        /**< value below which recfast uses the full equation for Hydrogen */
  double recfast_x_H0_trigger2;       /**< a second threshold used in derivative routine */
  double recfast_x_H0_trigger_delta;  /**< x_H range over which transition is smoothed */

  double recfast_H_frac;              /**< governs time at which full equation of evolution for Tmat is used */
/* @cond INCLUDE_WITH_DOXYGEN */
  FileName hyrec_Alpha_inf_file;
  FileName hyrec_R_inf_file;
  FileName hyrec_two_photon_tables_file;
/* @endcond */
  /* - for reionization */

  double reionization_z_start_max; /**< maximum redshift at which reionization should start. If not, return an error. */
  double reionization_sampling; /**< control stepsize in z during reionization */
  double reionization_optical_depth_tol; /**< fractional error on optical_depth */
  double reionization_start_factor; /**< parameter for CAMB-like parametrization */

  /* - general */

  int thermo_rate_smoothing_radius; /**< plays a minor (almost aesthetic) role in the definition of the variation rate of thermodynamical quantities */

  //@}

  /** @name - parameters related to the perturbation */

  //@{

  enum evolver_type evolver; /**< which type of evolver for integrating perturbations (Runge-Kutta? Stiff?...) */

  double k_min_tau0; /**< number defining k_min for the computation of Cl's and P(k)'s (dimensionless): (k_min tau_0), usually chosen much smaller than one */

  double k_max_tau0_over_l_max; /**< number defining k_max for the computation of Cl's (dimensionless): (k_max tau_0)/l_max, usually chosen around two */

  double k_step_sub; /**< step in k space, in units of one period of acoustic oscillation at decoupling, for scales inside sound horizon at decoupling */
  double k_step_super; /**< step in k space, in units of one period of acoustic oscillation at decoupling, for scales above sound horizon at decoupling */
  double k_step_transition; /**< dimensionless number regulating the transition from 'sub' steps to 'super' steps. Decrease for more precision. */
  double k_step_super_reduction; /**< the step k_step_super is reduced by this amount in the k-->0 limit (below scale of Hubble and/or curvature radius) */

  double k_per_decade_for_pk; /**< if values needed between kmax inferred from k_oscillations and k_kmax_for_pk, this gives the number of k per decade outside the BAO region*/

  double k_per_decade_for_bao; /**< if values needed between kmax inferred from k_oscillations and k_kmax_for_pk, this gives the number of k per decade inside the BAO region (for finer sampling)*/

  double k_bao_center; /**< in ln(k) space, the central value of the BAO region where sampling is finer is defined as k_rec times this number (recommended: 3, i.e. finest sampling near 3rd BAO peak) */

  double k_bao_width; /**< in ln(k) space, width of the BAO region where sampling is finer: this number gives roughly the number of BAO oscillations well resolved on both sides of the central value (recommended: 4, i.e. finest sampling from before first up to 3+4=7th peak) */

  double start_small_k_at_tau_c_over_tau_h; /**< largest wavelengths start being sampled when universe is sufficiently opaque. This is quantified in terms of the ratio of thermo to hubble time scales, \f$ \tau_c/\tau_H \f$. Start when start_largek_at_tau_c_over_tau_h equals this ratio. Decrease this value to start integrating the wavenumbers earlier in time. */

  double start_large_k_at_tau_h_over_tau_k;  /**< largest wavelengths start being sampled when mode is sufficiently outside Hubble scale. This is quantified in terms of the ratio of hubble time scale to wavenumber time scale, \f$ \tau_h/\tau_k \f$ which is roughly equal to (k*tau). Start when this ratio equals start_large_k_at_tau_k_over_tau_h. Decrease this value to start integrating the wavenumbers earlier in time. */

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

  int tight_coupling_approximation; /**< method for tight coupling approximation */

  int l_max_g;     /**< number of momenta in Boltzmann hierarchy for photon temperature (scalar), at least 4 */
  int l_max_pol_g; /**< number of momenta in Boltzmann hierarchy for photon polarization (scalar), at least 4 */
  int l_max_dr;   /**< number of momenta in Boltzmann hierarchy for decay radiation, at least 4 */
  int l_max_ur;   /**< number of momenta in Boltzmann hierarchy for relativistic neutrino/relics (scalar), at least 4 */
  int l_max_ncdm;   /**< number of momenta in Boltzmann hierarchy for relativistic neutrino/relics (scalar), at least 4 */
  int l_max_g_ten;     /**< number of momenta in Boltzmann hierarchy for photon temperature (tensor), at least 4 */
  int l_max_pol_g_ten; /**< number of momenta in Boltzmann hierarchy for photon polarization (tensor), at least 4 */

  double curvature_ini;     /**< initial condition for curvature for adiabatic */
  double entropy_ini; /**< initial condition for entropy perturbation for isocurvature */
  double gw_ini;      /**< initial condition for tensor metric perturbation h */

  /**
   * default step \f$ d \tau \f$ in perturbation integration, in units of the timescale involved in the equations (usually, the min of \f$ 1/k \f$, \f$ 1/aH \f$, \f$ 1/\dot{\kappa} \f$)
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

  int ur_fluid_approximation; /**< method for ultra relativistic fluid approximation */

  /**
   * when to switch off ur (massless neutrinos / ultra-relativistic
   * relics) fluid approximation
   */
  double ur_fluid_trigger_tau_over_tau_k;

  int ncdm_fluid_approximation; /**< method for non-cold dark matter fluid approximation */

  /**
   * when to switch off ncdm (massive neutrinos / non-cold
   * relics) fluid approximation
   */
  double ncdm_fluid_trigger_tau_over_tau_k;

  /**
   * whether CMB source functions can be approximated as zero when
   * visibility function g(tau) is tiny
   */
  double neglect_CMB_sources_below_visibility;

  //@}

  /** @name - parameters related to the primordial spectra */

  //@{

  double k_per_decade_primordial; /**< logarithmic sampling for primordial spectra (number of points per decade in k space) */

  double primordial_inflation_ratio_min; /**< for each k, start following wavenumber when aH = k/primordial_inflation_ratio_min */
  double primordial_inflation_ratio_max; /**< for each k, stop following wavenumber, at the latest, when aH = k/primordial_inflation_ratio_max */
  int primordial_inflation_phi_ini_maxit;      /**< maximum number of iteration when searching a suitable initial field value phi_ini (value reached when no long-enough slow-roll period before the pivot scale) */
  double primordial_inflation_pt_stepsize;     /**< controls the integration timestep for inflaton perturbations */
  double primordial_inflation_bg_stepsize;     /**< controls the integration timestep for inflaton background */
  double primordial_inflation_tol_integration; /**< controls the precision of the ODE integration during inflation */
  double primordial_inflation_attractor_precision_pivot;   /**< targeted precision when searching attractor solution near phi_pivot */
  double primordial_inflation_attractor_precision_initial; /**< targeted precision when searching attractor solution near phi_ini */
  int primordial_inflation_attractor_maxit; /**< maximum number of iteration when searching attractor solution */
  double primordial_inflation_tol_curvature; /**< for each k, stop following wavenumber, at the latest, when curvature perturbation R is stable up to to this tolerance */
  double primordial_inflation_aH_ini_target; /**< control the step size in the search for a suitable initial field value */
  double primordial_inflation_end_dphi; /**< first bracketing width, when trying to bracket the value phi_end at which inflation ends naturally */
  double primordial_inflation_end_logstep; /**< logarithmic step for updating the bracketing width, when trying to bracket the value phi_end at which inflation ends naturally */
  double primordial_inflation_small_epsilon; /**< value of slow-roll parameter epsilon used to define a field value phi_end close to the end of inflation (doesn't need to be exactly at the end): epsilon(phi_end)=small_epsilon (should be smaller than one) */
  double primordial_inflation_small_epsilon_tol; /**< tolerance in the search for phi_end */
  double primordial_inflation_extra_efolds; /**< a small number of efolds, irrelevant at the end, used in the search for the pivot scale (backward from the end of inflation) */

  //@}

  /** @name - parameters related to the transfer function */

  //@{

  int l_linstep; /**< factor for logarithmic spacing of values of l over which bessel and transfer functions are sampled */

  double l_logstep; /**< maximum spacing of values of l over which Bessel and transfer functions are sampled (so, spacing becomes linear instead of logarithmic at some point) */

  /* parameters relevant for bessel functions */
  double hyper_x_min;  /**< flat case: lower bound on the smallest value of x at which we sample \f$ \Phi_l^{\nu}(x)\f$ or \f$ j_l(x)\f$ */
  double hyper_sampling_flat;  /**< flat case: number of sampled points x per approximate wavelength \f$ 2\pi \f$*/
  double hyper_sampling_curved_low_nu;  /**< open/closed cases: number of sampled points x per approximate wavelength \f$ 2\pi/\nu\f$, when \f$ \nu \f$ smaller than hyper_nu_sampling_step */
  double hyper_sampling_curved_high_nu; /**< open/closed cases: number of sampled points x per approximate wavelength \f$ 2\pi/\nu\f$, when \f$ \nu \f$ greater than hyper_nu_sampling_step */
  double hyper_nu_sampling_step;  /**< open/closed cases: value of nu at which sampling changes  */
  double hyper_phi_min_abs;  /**< small value of Bessel function used in calculation of first point x (\f$ \Phi_l^{\nu}(x) \f$ equals hyper_phi_min_abs) */
  double hyper_x_tol;  /**< tolerance parameter used to determine first value of x */
  double hyper_flat_approximation_nu;  /**< value of nu below which the flat approximation is used to compute Bessel function */

  /* parameters relevant for transfer function */

  double q_linstep;         /**< asymptotic linear sampling step in q
                               space, in units of \f$ 2\pi/r_a(\tau_rec) \f$
                               (comoving angular diameter distance to
                               recombination) */

  double q_logstep_spline; /**< initial logarithmic sampling step in q
                                space, in units of \f$ 2\pi/r_a(\tau_{rec})\f$
                                (comoving angular diameter distance to
                                recombination) */

  double q_logstep_open;   /**< in open models, the value of
                                q_logstep_spline must be decreased
                                according to curvature. Increasing
                                this number will make the calculation
                                more accurate for large positive
                                Omega_k */

  double q_logstep_trapzd; /**< initial logarithmic sampling step in q
                                space, in units of \f$ 2\pi/r_a(\tau_{rec}) \f$
                                (comoving angular diameter distance to
                                recombination), in the case of small
                                q's in the closed case, for which one
                                must used trapezoidal integration
                                instead of spline (the number of q's
                                for which this is the case decreases
                                with curvature and vanishes in the
                                flat limit) */

  double q_numstep_transition; /**< number of steps for the transition
                                 from q_logstep_trapzd steps to
                                 q_logstep_spline steps (transition
                                 must be smooth for spline) */

  double transfer_neglect_delta_k_S_t0; /**< for temperature source function T0 of scalar mode, range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero */
  double transfer_neglect_delta_k_S_t1; /**< same for temperature source function T1 of scalar mode */
  double transfer_neglect_delta_k_S_t2; /**< same for temperature source function T2 of scalar mode */
  double transfer_neglect_delta_k_S_e;  /**< same for polarization source function E of scalar mode */
  double transfer_neglect_delta_k_V_t1; /**< same for temperature source function T1 of vector mode */
  double transfer_neglect_delta_k_V_t2; /**< same for temperature source function T2 of vector mode */
  double transfer_neglect_delta_k_V_e;  /**< same for polarization source function E of vector mode */
  double transfer_neglect_delta_k_V_b;  /**< same for polarization source function B of vector mode */
  double transfer_neglect_delta_k_T_t2; /**< same for temperature source function T2 of tensor mode */
  double transfer_neglect_delta_k_T_e;  /**< same for polarization source function E of tensor mode */
  double transfer_neglect_delta_k_T_b;  /**< same for polarization source function B of tensor mode */

  double transfer_neglect_late_source;  /**< value of l below which the CMB source functions can be neglected at late time, excepted when there is a Late ISW contribution */

  /** when to use the Limber approximation for project gravitational potential cl's */
  double l_switch_limber;

  /** when to use the Limber approximation for local number count contributions to cl's (relative to central redshift of each bin) */
  double l_switch_limber_for_nc_local_over_z;

  /** when to use the Limber approximation for number count contributions to cl's integrated along the line-of-sight (relative to central redshift of each bin) */
  double l_switch_limber_for_nc_los_over_z;

  /** in sigma units, where to cut gaussian selection functions */
  double selection_cut_at_sigma;

  /** controls sampling of integral over time when selection functions vary quicker than Bessel functions. Increase for better sampling. */
  double selection_sampling;

  /** controls sampling of integral over time when selection functions vary slower than Bessel functions. Increase for better sampling */
  double selection_sampling_bessel;

  /** controls sampling of integral over time when selection functions vary slower than Bessel functions. This parameter is specific to number counts contributions to Cl integrated along the line of sight. Increase for better sampling */
  double selection_sampling_bessel_los;

  /** controls how smooth are the edge of top-hat window function (<<1 for very sharp, 0.1 for sharp) */
  double selection_tophat_edge;

  //@}

  /** @name - parameters related to non-linear computations */

  //@{

  /** parameters relevant for HALOFIT computation */

  double halofit_min_k_nonlinear; /**< value of k in 1/Mpc below which
				     non-linear corrections will be neglected */

  double halofit_min_k_max; /**< when halofit is used, k_max must be
                               at least equal to this value (otherwise
                               halofit could not find the scale of
                               non-linearity). Calculations are done
                               internally until this k_max, but the
                               output is still controlled by
                               P_k_max_1/Mpc or P_k_max_h/Mpc even if
                               they are smaller */

  double halofit_k_per_decade; /**< halofit needs to evalute integrals
                                  (linear power spectrum times some
                                  kernels). They are sampled using
                                  this logarithmic step size. */

  double halofit_sigma_precision; /**< a smaller value will lead to a
				      more precise halofit result at the *highest*
				      redshift at which halofit can make computations,
				      at the expense of requiring a larger k_max; but
				      this parameter is not relevant for the
				      precision on P_nl(k,z) at other redshifts, so
				      there is normally no need to change it */

  double halofit_tol_sigma; /**< tolerance required on sigma(R) when
                               matching the condition sigma(R_nl)=1,
                               whcih defines the wavenumber of
                               non-linearity, k_nl=1./R_nl */

  double pk_eq_z_max; /**< Maximum z until which the Pk_equal method of 0810.0190 and 1601.07230 is used */

  double pk_eq_tol;   /**< tolerance for finding the equivalent models of the pk_equal method */

  //@}

  /** @name - parameters related to lensing */

  //@{

  int accurate_lensing; /**< switch between Gauss-Legendre quadrature integration and simple quadrature on a subdomain of angles */
  int num_mu_minus_lmax; /**< difference between num_mu and l_max, increase for more precision */
  int delta_l_max; /**< difference between l_max in unlensed and lensed spectra */
  double tol_gauss_legendre; /**< tolerance with which quadrature points are found: must be very small for an accurate integration (if not entered manually, set automatically to match machine precision) */
  //@}

  /** @name - general precision parameters */

  //@{

  double smallest_allowed_variation; /**< machine-dependent, assigned automatically by the code */

  //@}

  /** @name - zone for writing error messages */

  //@{

  ErrorMsg error_message;  /**< zone for writing error messages */

  //@}

};



#endif
