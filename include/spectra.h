/** @file spectra.h Documented includes for \f$ C_l^{X}, P(k), ... \f$ module */

#ifndef __SPECTRA__
#define __SPECTRA__

#include "primordial.h"

/**
 * All \f$ C_l^{X} or P(k,z) \f$.
 *
 * Once initialized by cl_init(), contains all multipoles \f$ C_l^{X} or P(k,z) \f$.
 */
struct spectra {
  
  int md_size; /**< number of modes included in computation */
  int * ic_size;       /**< for a given mode, ic_size[index_mode] = number of initial conditions included in computation */
  int * ic_ic_size;    /**< number of pairs of (index_ic1, index_ic2) with index_ic2 >= index_ic1; this number is just ic_size[index_mode](ic_size[index_mode]+1)/2  */
  short * * is_non_zero;  /**< is_non_zero[index_mode][index_ic1_ic2] */

  int * l_size; /**< number of multipole values for each requested mode, l_size[index_mode] */
  double ** l; /**< list of multipole values for each requested mode, (l[index_mode])[index_l] */
  
  double ** cl; /**< table of spectrum multipole \f$ C_l^{X} \f$'s for each mode, initial condition and cl_type, cl[index_mode][(index_l * psp->ic_ic_size[index_mode] + index_ic1_ic2) * psp->ct_size + index_ct] */
  double ** ddcl; /**< table of second derivatives w.r.t l in view of spline interpolation */ 

  int * l_max; /**< last multipole (given as an input) at which we trust our C_ls;
		  (l[index_mode][l_size[index_mode]-1] can be larger than l_max[index_mode], 
		  in order to ensure better interpolation with no boundary effects) */

  int l_max_tot; /**< greatest of all l_max[index_mode] */

  int index_ct_tt;
  int index_ct_ee;
  int index_ct_te;
  int index_ct_bb;
  int index_ct_pp;
  int index_ct_tp;
  int ct_size; /**< number of C_l types in the file of total spectra (spectra + tensors if any): TT, TE, EE, BB, phi-phi, T-phi, ... */

  int has_tt;
  int has_ee;
  int has_te;
  int has_bb;
  int has_pp;
  int has_tp;

  double z_max_pk; /**< maximum value of z at which matter spectrum P(k,z) will be evaluated; keep fixed to zero if P(k) only needed today */

  int index_md_scalars;
  int k_size;
  double * k;
  int eta_size; 
  double * eta;
  double * pk;   /**< power spectrum pk[(index_eta * psp->ic_ic_size[index_mode] + index_ic1_ic2) * psp->k_size + index_k] */
  double * ddlnpk; /**< table of second derivatives of ln(P) with respect to ln(eta) in view of spline interpolation */

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
		      double l,
		      double * cl,
		      double * * cl_md,
		      double * * cl_md_ic
		      );

  int spectra_pk_at_z(
		      struct background * pba,
		      struct spectra * psp,
		      double z,
		      double * pk,      /* pk[index_k] (already alloocated) */
		      double * pk_ic    /* pk_ic[index_k][index_ic1_ic2] (already allocated if more than one ic) */
		      );

  int spectra_pk_at_k_and_z(
			    struct background * pba,
			    struct primordial * ppm,
			    struct spectra * psp,
			    double k,
			    double z,
			    double * pk,
			    double * pk_ic   /* pk_ic[index_ic1_ic2] */
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
		      struct primordial * ppm,
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
			 int index_ic1,
			 int index_ic2,
			 int index_l,
			 int cl_integrand_num_columns,
			 double * cl_integrand,
			 double * primordial_pk,
			 double * transfer_ic1,
			 double * transfer_ic2
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
