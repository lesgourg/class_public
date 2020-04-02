/** @file lensing.h Documented includes for spectra module */

#ifndef __LENSING__
#define __LENSING__

#include "spectra.h"

/**
 * Structure containing everything about lensed spectra that other modules need to know.
 *
 * Once initialized by lensing_init(), contains a table of all lensed
 * \f$ C_l\f$'s for the all modes (scalar/tensor), all types (TT, TE...),
 * and all pairs of initial conditions (adiabatic, isocurvatures...).
 * FOR THE MOMENT, ASSUME ONLY SCALAR & ADIABATIC
 */

struct lensing {

  /** @name - input parameters initialized by user in input module
   *  (all other quantities are computed in this module, given these
   *  parameters and the content of the 'precision', 'background' and
   *  'thermodynamics' structures) */

  //@{

  short has_lensed_cls; /**< do we need to compute lensed \f$ C_l\f$'s at all ? */

  //@}

  /** @name - technical parameters */

  //@{

  short lensing_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

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


#ifdef __cplusplus
}
#endif

#endif
/* @endcond */
