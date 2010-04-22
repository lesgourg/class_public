/** @file spectra.h Documented includes for \f$ C_l^{X}, P(k), ... \f$ module */

#include "primordial.h"

#ifndef __SPECTRA__
#define __SPECTRA__

/**
 * All \f$ C_l^{X} or P(k,z) \f$.
 *
 * Once initialized by cl_init(), contains all multipoles \f$ C_l^{X} or P(k,z) \f$.
 */
struct spectra {
  
  int * l_size; /**< number of multipole values for each requested mode, l_size[index_mode] */
  double ** l; /**< list of multipole values for each requested mode, (l[index_mode])[index_l] */
  
  double ** cl; /**< table of spectrum multipole \f$ C_l^{X} \f$'s for each mode, initial condition and cl_type, (cl[index_mode])[index_ic][index_ct][index_l] */
  double ** ddcl; /**< table of second derivatives in view of spline interpolation */ 

  int index_ct_tt;
  int index_ct_ee;
  int index_ct_te;
  int index_ct_bb;
  int index_ct_pp;
  int index_ct_tp;
  int ct_size; /**< number of ct_type (TT, TE, EE, etc.) ct_size[index_mode]*/

  int k_size;
  int * k;
  double * pk;

  /** @name - flag regulating the amount of information sent to standard output (none if set to zero) */

  //@{

  short spectra_verbose;

  //@}

  ErrorMsg error_message; /**< zone for writing error messages */
};

/*************************************************************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int spectra_cl_at_l(
		      double l,
		      int index_mode,
		      double *cl
		      );

  int spectra_init(
	      struct perturbs * ppt,
	      struct transfers * ptr,
	      struct primordial * ppm,
	      struct spectra * psp
	      );

  int spectra_free();

  int spectra_indices();

#ifdef __cplusplus
}
#endif

#endif
