/** @file nonlinear.h Documented includes for trg module */

#include "primordial.h"

#ifndef __NONLINEAR__
#define __NONLINEAR__

#define _M_EV_TOO_BIG_FOR_HALOFIT_ 10. /**< above which value of non-CDM mass (in eV) do we stop trusting halofit? */

enum non_linear_method {nl_none,nl_halofit};
enum halofit_integral_type {halofit_integral_one, halofit_integral_two, halofit_integral_three};
enum halofit_statement {ok, too_small};

/**
 * Structure containing all information on non-linear spectra.
 *
 * Once initialized by nonlinear_init(), contains a table for all two points correlation functions
 * and for all the ai,bj functions (containing the three points correlation functions), for each
 * time and wave-number.
 */

struct nonlinear {

  /** @name - input parameters initialized by user in input module
      (all other quantities are computed in this module, given these
      parameters and the content of the 'precision', 'background',
      'thermo', 'primordial' and 'spectra' structures) */

  //@{

  enum non_linear_method method; /**< method for computing non-linear corrections (none, Halogit, etc.) */

  //@}

  /** @name - table non-linear corrections for matter density, sqrt(P_NL(k,z)/P_NL(k,z)) */

  //@{

  int k_size;      /**< k_size = total number of k values */
  double * k;      /**< k[index_k] = list of k values */
  int tau_size;    /**< tau_size = number of values */
  double * tau;    /**< tau[index_tau] = list of time values */

  double * nl_corr_density;   /**< nl_corr_density[index_tau * ppt->k_size + index_k] */
  double * k_nl;  /**< wavenumber at which non-linear corrections become important, defined differently by different non_linear_method's */
  int index_tau_min_nl; /**< index of smallest value of tau at which nonlinear corrections have been computed (so, for tau<tau_min_nl, the array nl_corr_density only contains some factors 1 */

  //@}

  /** @name - technical parameters */

  //@{

  short nonlinear_verbose;  	/**< amount of information written in standard output */

  ErrorMsg error_message; 	/**< zone for writing error messages */

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

  int nonlinear_k_nl_at_z(
                          struct background *pba,
                          struct nonlinear * pnl,
                          double z,
                          double * k_nl
                          );

  int nonlinear_init(
                     struct precision *ppr,
                     struct background *pba,
                     struct thermo *pth,
                     struct perturbs *ppt,
                     struct primordial *ppm,
                     struct nonlinear *pnl
                     );

  int nonlinear_free(
                     struct nonlinear *pnl
                     );

  int nonlinear_pk_l(struct perturbs *ppt,
                     struct primordial *ppm,
                     struct nonlinear *pnl,
                     int index_tau,
                     double *pk_l,
                     double *lnk,
                     double *lnpk,
                     double *ddlnpk);

  int nonlinear_halofit(
                        struct precision *ppr,
                        struct background *pba,
                        struct primordial *ppm,
                        struct nonlinear *pnl,
                        double tau,
                        double *pk_l,
                        double *pk_nl,
                        double *lnk_l,
                        double *lnpk_l,
                        double *ddlnpk_l,
                        double *k_nl,
                        enum halofit_statement * halofit_found_k_max
                        );

  int nonlinear_halofit_integrate(
                                  struct nonlinear *pnl,
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

#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
/* @endcond */
