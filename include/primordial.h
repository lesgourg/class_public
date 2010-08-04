/** @file primordial.h Documented includes for primordial module */

#ifndef __PRIMORDIAL__
#define __PRIMORDIAL__

#include "transfer.h"

enum primordial_spectrum_type {
  analytic_Pk,
  inflation_V,
  inflation_H
};

/**
 * All tables of transfer functions \f$ \Delta_l^{X} (k) \f$.
 *
 * Once initialized by transfer_init(), contains all tables of
 * transfer functions used for interpolation in other modules, for all
 * requested modes (scalar/vector/tensor), initial conditions, type
 * (temperature, polarization, etc), l and k.
 */
struct primordial {

  short has_scalars;
  short has_vectors;
  short has_tensors;

  short has_ad;     
  short has_bi;     
  short has_cdi;    
  short has_nid;    
  short has_niv;    

  int index_md_scalars; /**< index value for scalars */
  int index_md_tensors; /**< index value for tensors */
  int index_md_vectors; /**< index value for vectors */
  int md_size; /**< number of modes included in computation */

  int index_ic_ad; /**< index value for adiabatic */
  int index_ic_cdi; /**< index value for CDM isocurvature */
  int index_ic_bi; /**< index value for baryon isocurvature */
  int index_ic_nid; /**< index value for neutrino density isocurvature */
  int index_ic_niv; /**< index value for neutrino velocity isocurvature */
  int index_ic_ten; /**< index value for unique possibility for tensors */
  int * ic_size;       /**< for a given mode, ic_size[index_mode] = number of initial conditions included in computation */

  enum primordial_spectrum_type primordial_spec_type;

  double k_pivot; /* pivot scale in Mpc-1 */

  double A_s_ad;  /* scalar amplitude (adiabatic) */
  double n_s_ad;  /* scalar tilt (adiabatic) */
  double alpha_s_ad; /* scalar running (adiabatic) */

  double r;  /* tensor to scalar ratio A_T/A_S=P_h/P_R  */
  double n_t;  /* tensor tilt */
  double alpha_t; /* tensor running */

  double * lnk; /* list of ln(k) values lnk[index_k] */
  int lnk_size; /* number of ln(k) values */

  double ** lnpk; /* primordial spectra (lnP[index_mode])[index_ic][index_k] */
  double ** ddlnpk; /* second derivative of lnP for spline interpolation */

  /** @name - flag regulating the amount of information sent to standard output (none if set to zero) */

  //@{

  short primordial_verbose;

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

  int primordial_spectrum_at_k(
			       struct primordial * ppm,
			       int index_mode,
			       double k,
			       double * pk
			       );

  int primordial_init(
		      struct precision  * ppr,
		      struct perturbs   * ppt,
		      struct primordial * ppm
		      );
    
  int primordial_free(
		      struct primordial * ppm
		      );
    
  int primordial_indices(
			 struct perturbs   * ppt,
			 struct primordial * ppm
			 );

  int primordial_analytic_spectrum(
				   struct primordial * ppm,
				   int index_mode,
				   int index_ic,
				   double k,
				   double * pk
				   );

  int primordial_get_lnk_list(
			      struct primordial * ppm,
			      double kmin,
			      double kmax,
			      double k_per_decade
			      );

#ifdef __cplusplus
}
#endif

/**  
 * @name Constants
 */

//@{


//@}

#endif
