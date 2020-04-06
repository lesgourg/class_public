/** @file spectra.h Documented includes for spectra module */

#ifndef __SPECTRA__
#define __SPECTRA__

#include "transfer.h"

/**
 * Structure containing everything about anisotropy and Fourier power spectra that other modules need to know.
 *
 * Once initialized by spectra_init(), contains a table of all
 * \f$ C_l\f$'s and P(k) as a function of multipole/wavenumber,
 * mode (scalar/tensor...), type (for \f$ C_l\f$'s: TT, TE...),
 * and pairs of initial conditions (adiabatic, isocurvatures...).
 */

struct spectra {

  /** @name - input parameters initialized by user in input module
      (all other quantities are computed in this module, given these parameters
      and the content of the 'background', 'perturbs', 'transfers' and
      'primordial' structures) */

  //@{

  double z_max_pk;  /**< maximum value of z at which matter spectrum P(k,z) will be evaluated; keep fixed to zero if P(k) only needed today */

  int non_diag; /**< sets the number of cross-correlation spectra
                   that you want to calculate: 0 means only
                   auto-correlation, 1 means only adjacent bins,
                   and number of bins minus one means all
                   correlations */

  

  /** @name - technical parameters */

  //@{

  struct nonlinear * pnl; /**< a pointer to the nonlinear structure is
                            stored in the spectra structure. This odd,
                            unusual and unelegant feature has been
                            introduced in v2.8 in order to keep in use
                            some deprecated functions spectra_pk_...()
                            that are now pointing at new function
                            nonlinear_pk_...(). In the future, if the
                            deprecated functions are removed, it will
                            be possible to remove also this pointer. */

  short spectra_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */

  //@}
};

/*************************************************************************************************************/
/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif


  /* end deprecated functions */

#ifdef __cplusplus
}
#endif

#endif
/* @endcond */
