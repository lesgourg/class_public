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
  
  int md_size; /**< number of modes included in computation */
  int * ic_size;       /**< for a given mode, ic_size[index_mode] = number of initial conditions included in computation */

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
  int ct_size; /**< number of ct_type (TT, TE, EE, etc.) */


  double z_max_pk; /**< maximum value of z at which matter spectrum P(k,z) will be evaluated; keep fixed to zero if P(k) only needed today */

  int k_size;
  double * k;
  int eta_size; 
  double * eta;
  double * pk;   /**< power spectrum pk[(index_ic * psp->k_size + index_k)*eta_size+index_eta] */
  double * ddpk; /**< table of second derivatives with respect to eta in view of spline interpolation */

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
		      struct spectra * psp,
		      int index_mode,
		      double l,
		      double *cl
		      );

  int spectra_pk_at_z(
		      struct background * pba,
		      struct spectra * psp,
		      int index_mode,
		      double z,
		      double * pk
		      );

  int spectra_pk_at_k_and_z(
			    struct background * pba,
			    struct primordial * ppm,
			    struct spectra * psp,
			    int index_mode,
			    int index_ic,
			    double k,
			    double z,
			    double * pk
			    );

  int spectra_init(
		   struct background * pba,
		   struct perturbs * ppt,
		   struct transfers * ptr,
		   struct primordial * ppm,
		   struct spectra * psp
		   );

  int spectra_free(
		   struct spectra * psp
		   );

  int spectra_indices(
		      struct perturbs * ppt,
		      struct transfers * ptr,
		      struct spectra * psp
		      );

  int spectra_cls(
		  struct perturbs * ppt,
		  struct transfers * ptr,
		  struct primordial * ppm,
		  struct spectra * psp
		  );

  int spectra_compute_cl(
			 struct perturbs * ppt,
			 struct transfers * ptr,
			 struct primordial * ppm,
			 struct spectra * psp,
			 int index_mode,
			 int index_ic,
			 int index_l,
			 int cl_integrand_num_columns,
			 double * cl_integrand,
			 double * primordial_pk,
			 double * transfer
			 );
  
  int spectra_pk(
		 struct background * pba,
		 struct perturbs * ppt,
		 struct primordial * ppm,
		 struct spectra * psp
		 );

#ifdef __cplusplus
}
#endif

#endif
