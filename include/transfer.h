/** @file transfer.h Documented includes for transfer module. */

#ifndef __TRANSFER__
#define __TRANSFER__

#include "bessel.h"
#include "perturbations.h"

/**
 * Structure containing everything about transfer functions in harmonic space \f$ \Delta_l^{X} (k) \f$ that other modules need to know.
 *
 * Once initialized by transfer_init(), contains all tables of
 * transfer functions used for interpolation in other modules, for all
 * requested modes (scalar/vector/tensor), initial conditions, type
 * (temperature, polarization, etc), multipole l and wavenumber k.
 * 
 * The content of this structure is entirely computed in this module,
 * given the content of the 'precision', 'bessels', 'background',
 * 'thermodynamics' and 'perturbation' structures.
 */

struct transfers {

  /** @name - input parameters initialized by user in input module
   *  (all other quantitites are computed in this module, given these
   *  parameters and the content of previous structures) */
  
  //@{

  double lcmb_rescale; /**< normally set to one, can be used
			  excelptionally to rescale by hand the CMB
			  lensing potential */
  double lcmb_tilt;    /**< normally set to zero, can be used
			  excelptionally to tilt by hand the CMB
			  lensing potential */
  double lcmb_pivot;   /**< if lcmb_tilt non-zero, corresponding pivot
			  scale */

  //@}

  /** @name - flag stating whether we need transfer functions at all */

  //@{

  short has_cls; /**< copy of same flag in perturbation structure */

  //@}

  /** @name - number of modes and transfer function types */

  //@{

  int md_size;       /**< number of modes included in computation */

  int index_tt_t;      /**< index for transfer type = temperature */
  int index_tt_e;      /**< index for transfer type = E-polarization */
  int index_tt_b;      /**< index for transfer type = B-polarization */
  int index_tt_lcmb;   /**< index for transfer type = CMB lensing */
  int index_tt_density; /**< index for first bin of transfer type = matter density */
  int index_tt_lensing; /**< index for first bin of transfer type = galaxy lensing */

  int * tt_size;     /**< number of requested transfer types tt_size[index_md] for each mode */

  //@}

  /** @name - number and list of multipoles */

  //@{

  int ** l_size_tt;  /**< number of multipole values for which we effectively compute the transfer function,l_size[index_md][index_tt] */ 

  int * l_size;   /**< number of multipole values for each requested mode, l_size[index_md] */

  int l_size_max; /**< greatest of all l_size[index_md] */

  int * l;        /**< list of multipole values l[index_l] */

  //@}

  /** @name - number and list of wavenumbers */

  //@{

  int * k_size; /**< number of wavenumber values for each requested mode, k_size[index_md] */

  double ** k;  /**< list of wavenumber values for each requested mode, k[index_md][index_k] */

  //@}

  /** @name - transfer functions */

  //@{

  double ** transfer; /**< table of transfer functions for each mode, initial condition, type, multipole and wavenumber, with argument transfer[index_md][((index_ic * ptr->tt_size[index_md] + index_tt) * ptr->l_size[index_md] + index_l) * ptr->k_size[index_md] + index_k] */

  //@}

  /** @name - technical parameters */

  //@{

  short transfer_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */

  //@}
};


/*************************************************************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int transfer_functions_at_k(
			      struct transfers * ptr,
			      int index_md,
			      int index_ic,
			      int index_type,
			      int index_l,
			      double k,
			      double * ptransfer_local
			      );

  int transfer_init(
		    struct precision * ppr,
		    struct background * pba,
		    struct thermo * pth,
		    struct perturbs * ppt,
		    struct bessels * pbs,
		    struct transfers * ptr
		    );
    
  int transfer_free(
		    struct transfers * ptr
		    );

  int transfer_indices_of_transfers(
				    struct precision * ppr,
				    struct perturbs * ppt,
				    struct bessels * pbs,
				    struct transfers * ptr,
				    double tau0
				    );

  int transfer_get_l_list(
			  struct precision * ppr,
			  struct perturbs * ppt,
			  struct bessels * pbs,
			  struct transfers * ptr
			  );

  int transfer_get_k_list(
			  struct precision * ppr,
			  struct perturbs * ppt,
			  struct transfers * ptr,
			  double tau0,
			  int index_md
			  );

  int transfer_get_source_correspondence(
					 struct perturbs * ppt,
					 struct transfers * ptr,
					 int index_md,
					 int * tp_of_tt
					 );
  
  int transfer_source_tau_size(
			       struct precision * ppr,
			       struct background * pba,
			       struct perturbs * ppt,
			       struct transfers * ptr,
			       double tau_rec,
			       double tau0,
			       int index_md,
			       int index_tt,
			       int * tau_size
			       );

  int transfer_interpolate_sources(
				   struct perturbs * ppt,
				   struct transfers * ptr,
				   int index_md,
				   int index_ic,
				   int index_type,
				   double * source_spline,
				   double * interpolated_sources
				   );

  int transfer_sources(
		       struct precision * ppr,
		       struct background * pba,
		       struct perturbs * ppt,
		       struct transfers * ptr,
		       double * interpolated_sources,
		       double tau_rec,
		       int index_md,
		       int index_tt,
		       double * sources,
		       double * tau0_minus_tau,
		       double * delta_tau,
		       double * tau_size_double
		       );

  int transfer_integration_time_steps(
				      struct transfers * ptr,
				      double * tau0_minus_tau,
				      int tau_size,
				      double * delta_tau
				      );

  int transfer_selection_function(
				  struct precision * ppr,
				  struct perturbs * ppt,
				  struct transfers * ptr,
				  int bin,
				  double z,
				  double * selection);

  int transfer_selection_sampling(
				  struct precision * ppr,
				  struct background * pba,
				  struct perturbs * ppt,
				  struct transfers * ptr,
				  int bin,
				  double * tau0_minus_tau,
				  int tau_size);
  
  int transfer_lensing_sampling(
				  struct precision * ppr,
				  struct background * pba,
				  struct perturbs * ppt,
				  struct transfers * ptr,
				  int bin,
				  double tau0,
				  double * tau0_minus_tau,
				  int tau_size);
  
  int transfer_source_resample(
			       struct precision * ppr,
			       struct background * pba,
			       struct perturbs * ppt,
			       struct transfers * ptr,
			       int bin,
			       double * tau0_minus_tau,
			       int tau_size,
			       int index_md,
			       double tau0,
			       double * interpolated_sources,
			       double * sources);

  int transfer_selection_times(
			       struct precision * ppr,
			       struct background * pba,
			       struct perturbs * ppt,
			       struct transfers * ptr,
			       int bin,
			       double * tau_min,
			       double * tau_mean,
			       double * tau_max);
  
  int transfer_selection_compute(
				 struct precision * ppr,
				 struct background * pba,
				 struct perturbs * ppt,
				 struct transfers * ptr,
				 double * selection,
				 double * tau0_minus_tau,
				 double * delta_tau,
				 int tau_size,
				 double * pvecback,
				 double tau0,
				 int bin);

  int transfer_compute_for_each_l(
				  struct precision * ppr,
				  struct perturbs * ppt,
				  struct transfers * ptr,
				  int index_md,
				  int index_ic,
				  int index_tt,
				  int index_l,
				  double l,
				  double x_min_l,
				  double x_step,
				  double * tau0_minus_tau,
				  double * delta_tau,
				  int tau_size,
				  double * sources,
				  double * j_l,
				  double * ddj_l,
				  double * dj_l,
				  double * dddj_l,
				  double k_max_bessel
				  );

  int transfer_use_limber(
			  struct precision * ppr,
			  struct perturbs * ppt,
			  struct transfers * ptr,
			  double k_max_bessel,
			  int index_md,
			  int index_tt,
			  double k,
			  double l,
			  short * use_limber
			  );

  int transfer_integrate(
			 struct transfers * ptr,
			 int tau_size,
			 int index_k,
			 double l,
			 double k,
			 double x_min_l,
			 double x_step,
			 double * tau0_minus_tau,
			 double * delta_tau,
			 double * sources,
			 double *j_l,
			 double *ddj_l,
			 double * trsf
			 );
    
  int transfer_limber(
		      int tau_size,
		      struct transfers * ptr,
		      int index_md,
		      int index_k,
		      double l,
		      double k,
		      double * tau0_minus_tau,
		      double * sources,
		      double * trsf
		      );
  
  int transfer_limber2(
		       int tau_size,
		       struct transfers * ptr,
		       int index_md,
		       int index_k,
		       double l,
		       double k,
		       double * tau0_minus_tau,
		       double * sources,
		       double * trsf
		       );
  
  int transfer_envelop(
		       int tau_size,
		       int index_k,
		       double l,
		       double k,
		       double x_min_l,
		       double x_step,
		       double * tau0_minus_tau,
		       double * delta_tau,
		       double * sources,
		       double *j_l,
		       double *ddj_l,
		       double *dj_l,
		       double *dddj_l,
		       double * trsf
		       );
    
#ifdef __cplusplus
}
#endif

#endif
