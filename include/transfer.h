/** @file transfer.h Documented includes for transfer module */

#include "bessel.h"

#ifndef __TRANSFER__
#define __TRANSFER__

/**
 * All tables of transfer functions \f$ \Delta_l^{X} (k) \f$.
 *
 * Once initialized by transfer_init(), contains all tables of
 * transfer functions used for interpolation in other modules, for all
 * requested modes (scalar/vector/tensor), initial conditions, type
 * (temperature, polarization, etc), l and k.
 */
struct transfers {

  int index_tt_t; /**< index for transfer type = temperature */
  int index_tt_p; /**< index for transfer type = temperature */
  int index_tt_lcmb; /**< index for transfer type = CMB lensing */
  int tt_size;    /**< number of requested transfer types */

  int * l_size; /**< number of multipole values for each requested mode, l_size[index_mode] */
  int ** l; /**< list of multipole values for each requested mode, (l[index_mode])[index_l] */

  int * k_size; /**< number of wavenumber values for each requested mode, k_size[index_mode] */
  double ** k; /**< list of wavenumber values for each requested mode, (k[index_mode])[index_k] */

  double ** transfer; /**< table of transfer functions for each mode, initial condition and type, (transfer[index_mode])[index_ic][index_type][index_l][index_k] */

  /** @name - flag regulating the amount of information sent to standard output (none if set to zero) */

  //@{

  short transfer_verbose;

  //@}

  ErrorMsg error_message; /**< zone for writing error messages */
};

/**
 * Table of integrand of transfer function, together with their second
 * derivatives for spline interpolation:
 */
struct transfer_integrand {

double * trans_int; /* table of integrand \f$ S(k,\eta)*j_l[k(\eta_0-\eta)] \f$ as a function of \f$ \eta \f%, as well as its splined second derivative */

int trans_int_eta; /* index of column for time */
int trans_int_y; /* index of column for integrand */
int trans_int_ddy; /* index of column for secoind derivative of integrand */
int trans_int_col_num; /* it number of columns */
};

/*************************************************************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int transfer_functions_at_k(
			      int index_mode,
			      int index_ic,
			      int index_type,
			      int index_l,
			      double k,
			      double * ptransfer_local
			      );

  int transfer_init(
		    struct background * pba_input,
		    struct thermo * pth_input,
		    struct perturbs * ppt_input,
		    struct bessels * pbs_input,
		    struct precision * ppr_input,
		    struct transfers * ptr_output
		    );
    
  int transfer_free();

  int transfer_indices_of_transfers();

  int transfer_get_l_list_size(
			       int index_mode,
			       int * pl_list_size
			       );

  int transfer_get_l_list(
			  int index_mode,
			  int * pl_list
			  );

  int transfer_get_k_list_size(
			       int index_mode,
			       int * pk_list_size
			       );

  int transfer_get_k_list(
			  int index_mode,
			  double * pk_list
			  );

  int transfer_interpolate_sources(
				   int current_index_mode,
				   int current_index_ic,
				   int current_index_type,
				   int current_index_l,
				   double * source_spline,
				   double * interpolated_sources
				   );

  int transfer_integrate(
			 int current_index_mode,
			 int current_index_ic,
			 int current_index_type,
			 int current_index_l,
			 int current_index_k,
			 double * interpolated_sources,
			 struct transfer_integrand * pti,
			 double * trsf
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
