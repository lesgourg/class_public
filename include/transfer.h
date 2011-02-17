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

  /** @name - flag stating whether we need transfer functions at all */

  //@{

  short has_cls; /**< copy of same flag in perturbation structure */

  //@}

  /** @name - number of modes and transfer function types */

  //@{

  int md_size;       /**< number of modes included in computation */

  int index_tt_t;    /**< index for transfer type = temperature */
  int index_tt_e;    /**< index for transfer type = E-polarization */
  int index_tt_b;    /**< index for transfer type = B-polarization */
  int index_tt_lcmb; /**< index for transfer type = CMB lensing */

  int * tt_size;     /**< number of requested transfer types tt_size[index_mode] for each mode */

  //@}

  /** @name - number and list of multipoles */

  //@{

  int * l_size; /**< number of multipole values for each requested mode, l_size[index_mode] */

  int ** l;     /**< list of multipole values for each requested mode, l[index_mode][index_l] */

  //@}

  /** @name - number and list of wavenumbers */

  //@{

  int * k_size; /**< number of wavenumber values for each requested mode, k_size[index_mode] */

  double ** k;  /**< list of wavenumber values for each requested mode, k[index_mode][index_k] */

  //@}

  /** @name - transfer functions */

  //@{

  double ** transfer; /**< table of transfer functions for each mode, initial condition, type, multipole and wavenumber, with argument transfer[index_mode][((index_ic * ptr->tt_size[index_mode] + index_tt) * ptr->l_size[index_mode] + index_l) * ptr->k_size[index_mode] + index_k] */

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
			      int index_mode,
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
				    double rs_rec
				    );

  int transfer_get_l_list(
			  struct precision * ppr,
			  struct perturbs * ppt,
			  struct bessels * pbs,
			  struct transfers * ptr,
			  int index_mode
			  );

  int transfer_get_k_list(
			  struct precision * ppr,
			  struct perturbs * ppt,
			  struct transfers * ptr,
			  double rs_rec,
			  int index_mode
			  );

  int transfer_interpolate_sources(
				   struct perturbs * ppt,
				   struct transfers * ptr,
				   double eta0,
				   double eta_rec,
				   int current_index_mode,
				   int current_index_ic,
				   int current_index_type,
				   double * source_spline,
				   double * interpolated_sources
				   );

  int transfer_compute_for_each_l(
				  struct precision * ppr,
				  struct perturbs * ppt,
				  struct transfers * ptr,
				  int index_mode,
				  int index_ic,
				  int index_tt,
				  int index_l,
				  double l,
				  double x_min_l,
				  double x_step,
				  double * eta0_minus_eta,
				  double * delta_eta,
				  double * sources,
				  double * j_l,
				  double * ddj_l
				  );

  int transfer_integrate(
			 int eta_size,
			 int index_k,
			 double l,
			 double k,
			 double x_min_l,
			 double x_step,
			 double * eta0_minus_eta,
			 double * delta_eta,
			 double * sources,
			 double *j_l,
			 double *ddj_l,
			 double * trsf
			 );
    
  int transfer_limber(
		      int eta_size,
		      struct transfers * ptr,
		      int index_mode,
		      int index_k,
		      double l,
		      double k,
		      double * eta0_minus_eta,
		      double * sources,
		      double * trsf
		      );
  
#ifdef __cplusplus
}
#endif

#endif
