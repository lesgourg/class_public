/** @file transfer.c Documented transfer module.
 * Julien Lesgourgues, 18.04.2010    
 *
 * This module has two purposes: 
 *
 * - at the beginning, to compute the transfer functions \f$
 *   \Delta_l^{X} (k) \f$, and store them in tables used for
 *   interpolation in other modules.
 *
 * - at any time in the code, to evaluate the transfer functions (for
 *   a given mode, initial condition, type and multipole l) at any
 *   wavenumber k (by interpolating within the interpolation table).
 *
 * Hence the following functions can be called from other modules:
 *
 * -# transfer_init() at the beginning (but after background_init(), thermodynamics_init() and perturb_init())
 * -# transfer_functions_at_k() at any later time 
 * -# transfer_free() at the end, when no more calls to transfer_functions_at_k() are needed
 */

#include "transfer.h"

/** @name - structures used within the transfer module: */

//@{

struct precision * ppr; /**< a precision_params structure pointer for internal use in the perturbation module */
struct perturbs * ppt; /**< a perturbs structure pointer for internal use in the perturbation module */
struct bessels * pbs; /**< a bessels structure pointer for internal use in the perturbation module */
struct transfers * ptr; /**< a transfers structure pointer for internal use in the perturbation module */

//@}

/** @name - instead of having pba and pth as global variables, only
    share these two fields (the only one needed throughout the modules
    from these structures) */

//@{

double eta0; /* conformal age (conformal time today) */
double eta_rec; /* conformal time at recombination */

//@}

/** @name - miscellaneous: */

//@{

char * errmsg; /**< error management pointer */
char Transmit_Error_Message[2048]; /**< contains error message */

//@}

/** Transfer function \f$ \Delta_l^{X} (k) \f$ at a given wavenumber k.
 *
 * For a given mode (scalar, vector, tensor), initial condition, type
 * (temperature, polarization, lensing, etc) and multipole, compute
 * the transfer function for an arbitary value of k by interpolating
 * between pre-computed values of k. Output the result as a
 * two-dimensional vector pvectransfer_local = (k, transfer). This
 * function can be called from whatever module at whatever time,
 * provided that transfer_init() has been called before, and
 * transfer_free() has not been called yet.
 *
 * @param index_mode Input: index of requested mode
 * @param index_ic Input: index of requested initial condition
 * @param index_type Input: index of requested type
 * @param index_l Input: index of requested multipole
 * @param k Input: any wavenumber
 * @param ptransfer_local Output: transfer function
 * @return the error status
 */
int transfer_functions_at_k(
			    int index_mode,
			    int index_ic,
			    int index_type,
			    int index_l,
			    double k,
			    double * ptransfer_local 
			    ) {
  /** Summary: */

  /** - interpolate in pre-computed table using array_interpolate_two() */
  if (array_interpolate_two(
			    ppt->k[index_mode],
			    1,
			    1,
			    ptr->transfer[index_mode] + ( index_ic * ppt->tp_size + index_type) * ptr->l_size[index_mode] + index_l,
			    1,
			    ppt->k_size[index_mode],
			    k,
			    ptransfer_local,
			    1,
			    errmsg) == _FAILURE_) {
    sprintf(ptr->error_message,"%s(L:%d) : error in array_interpolate_two() \n=>%s",__func__,__LINE__,errmsg);
    return _FAILURE_;
  }
  
  return _SUCCESS_;
}

/** Initialize the transfers structure, including transfer functions interpolation table.
 *
 * This function: 
 * - initializes all indices in the transfers structure
 *   and allocate all its arrays using transfer_indices_of_transfers().
 * - for a all requested modes (scalar, vector, tensor),
 *   initial conditions and types (temperature, polarization, lensing,
 *   etc), compute the transfer function \f$ \Delta_l^{X} (k) \f$
 *   following the following steps:
 * -# interpolate sources \f$ S(k, \eta) \f$ to get them at the right values of k 
 *    using transfer_interpolate_sources()
 * -# for each k, resample sources at the right values of \f$ \eta \f$ with 
 *    transfer_resample_sources()
 * -# for each k and l, compute the transfer function by convolving the sources with 
 *    the Bessel functions using transfer_solve()
 * -# store result in the transfer table 
 *    (transfer[index_mode])[index_ic][index_type][index_l][index_k]
 *
 * This function shall be called at the beginning of each run, but
 * only after background_init(), thermodynamics_init() and perturb_init(). It
 * allocates memory spaces which should be freed later with
 * transfer_free().
 *
 * @param pba_input Input : Initialized background structure 
 * @param pth_input Input : Initialized thermodynamics structure 
 * @param ppt_input Input : Initialized perturbation structure
 * @param pbs_input Input : Initialized bessels structure
 * @param ptr_input Input : Parameters describing how the computation is to be performed
 * @param ptr_output Output : Initialized transfers structure
 * @return the error status
 */
int transfer_init(
		  struct background * pba_input,
		  struct thermo * pth_input,
		  struct perturbs * ppt_input,
		  struct bessels * pbs_input,
		  struct precision * ppr_input,
		  struct transfers * ptr_output
		  ) {

  /** Summary: */

  /** - define local variables */

  /* running index for modes */
  int index_mode; 
  /* running index for initial conditions */
  int index_ic; 
  /* running index for wavenumbers */
  int index_k; 
  /* running index for types */
  int index_type; 
  /* running index for multipoles */
  int index_l; 
  /* another index */
  int index;
  /* current wavenumber value */
  double current_k;
  /* result for each transfer function */
  double transfer_function;
  /* flag: for a given l, should the transfer functions stop being computed at next value of k? */
  short cut_transfer;
  /* global maximum of \f$ \Delta_l(k) \f$ as a function of k, used for stopping computation */
  double global_max;
  /* global minimum of \f$ \Delta_l(k) \f$ as a function of k, used for stopping computation */
  double global_min;
  /* last local maximum of \f$ \Delta_l(k) \f$ as a function of k, used for stopping computation */
  double last_local_max;
  /* last local minimum of \f$ S(k, \eta) j_l(k [\eta_0 - \eta]) \f$ as a function of k, used for cutting the integral */
  double last_local_min;
  /* in loop over k, last computed value of transfer function \f$ \Delta_l(k) \f$ */
  double transfer_last;
  /* in loop over k, value of transfer function \f$ \Delta_l(k) \f$ computed two steps before */
  double transfer_last_last;
  /* rough estimate of C_l, used only for cutting the transfer function computation at some k_max */
  double cl;
  /* rough estimate of C_l's variation, used only for cutting the transfer function computation at some k_max */
  double delta_cl;
  /* rough estimate of C_l's fractional variation, used only for cutting the transfer function computation at some k_max */
  double cl_var;
  /* last computed value of C_l's fractional variation, used only for cutting the transfer function computation at some k_max */
  double cl_var_last;
  /* value of C_l's fractional variation computed two steps ago, used only for cutting the transfer function computation at some k_max */
  double cl_var_last_last;
  /* table of source functions interpolated at the right values of k, interpolated_sources[index_k][index_eta] */
  double * interpolated_sources;
  /* array of splines values S''(k,eta) (second derivative with respect to k, not eta!) */
  double * source_spline;
  /* table of integrand of transfer function */
  struct transfer_integrand ti;

  if (ptr_output->has_cls == _FALSE_)
    return _SUCCESS_;

  /** - get conformal age / recombination time from background / thermodynamics structures (only place where these structures are used in this module) */
  eta0 = pba_input->conformal_age;
  eta_rec = pth_input->eta_rec;

  /** - identify the perturbs, precision_params, bessels and transfers structures ppt, ppr, pbs, ptr (used throughout transfer.c as global variables) to the input/output structures of this function (ppt, ppr, pbs are already filled, ptr will be filled by this function) */
  ppt = ppt_input;
  ppr = ppr_input;
  pbs = pbs_input;
  ptr = ptr_output; 

  if (ptr->transfer_verbose > 0)
    printf("Computing transfers\n");

  /** - initialize all indices in the transfers structure and allocate all its arrays using transfer_indices_of_transfers() */
  if (transfer_indices_of_transfers() == _FAILURE_) {
    sprintf(Transmit_Error_Message,"%s(L:%d) : error in transfer_indices_of_transfers() \n=>%s",__func__,__LINE__,ptr->error_message);
    sprintf(ptr->error_message,Transmit_Error_Message);
    return _FAILURE_;
  }

  /** - initialize all indices in the transfer_integrand structure and allocate its array. Fill the eta column. */

  index = 0;
  ti.trans_int_eta = index;
  index++;
  ti.trans_int_y = index;
  index++;
  ti.trans_int_ddy = index;
  index++;
  ti.trans_int_col_num = index;

  ti.trans_int =  malloc(sizeof(double) * ppt->eta_size * ti.trans_int_col_num);
  if (ti.trans_int==NULL) {
    sprintf(ptr->error_message,"%s(L:%d): Cannot allocate ti.trans_int",__func__,__LINE__);
    return _FAILURE_;
  }

  for (index=0; index < ppt->eta_size; index++)
    ti.trans_int[ti.trans_int_col_num*index+ti.trans_int_eta] = ppt->eta_sampling[index];

  /** - loop over all indices of the table of transfer functions. For each mode, initial condition and type: */

  for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {

    /** (a) allocate temporary arrays relevant for this mode */

    interpolated_sources = malloc(ptr->k_size[index_mode]
				  * ppt->eta_size
				  * sizeof(double));
    if (interpolated_sources==NULL) {
      sprintf(ptr->error_message,"%s(L:%d): Cannot allocate interpolated_sources",__func__,__LINE__);
      return _FAILURE_;
    }

    source_spline = malloc(sizeof(double) 
			   * ppt->eta_size 
			   * ppt->k_size[index_mode]);
    if (source_spline==NULL) {
      sprintf(ptr->error_message,"%s(L:%d): Cannot allocate source_spline",__func__,__LINE__);
      return _FAILURE_;
    }

    for (index_ic = 0; index_ic < ppt->ic_size[index_mode]; index_ic++) {

      for (index_type = 0; index_type < ppt->tp_size; index_type++) {

	/** (b) interpolate sources to get them at the right values of k using transfer_interpolate_sources() */

	if (transfer_interpolate_sources(
					 index_mode,
					 index_ic,
					 index_type,
					 index_l,
					 source_spline,
					 interpolated_sources) == _FAILURE_) {
	  sprintf(Transmit_Error_Message,"%s(L:%d) : error in transfer_interpolate_sources()\n=>%s",__func__,__LINE__,ptr->error_message);
	  sprintf(ptr->error_message,Transmit_Error_Message);
	  return _FAILURE_;
	}

	/** (c) loop over l. For each value of l: */

	/***** THIS IS THE LOOP WHICH SHOULD BE PARALLELISED ******/
	for (index_l = 0; index_l < ptr->l_size[index_mode]; index_l++) {

	  if (ptr->transfer_verbose > 1)
	    printf("Compute transfer for l=%d\n",ptr->l[index_mode][index_l]);
	  
	  /** (c.a) if the option of stopping the transfer function computation at some k_max is selected, initialize relevant quantities */
	  
	  if (ppr->transfer_cut == tc_osc) {
	    cut_transfer = _FALSE_;
	    global_max=0.;
	    global_min=0.;
	    last_local_max=0.;
	    last_local_min=0.;
	    transfer_function=0;
	    transfer_last=0;
	    transfer_last_last=0;
	  }

	  if (ppr->transfer_cut == tc_cl) {
	    cut_transfer = _FALSE_;
	    cl=0.;
	    cl_var=1.;
	    cl_var_last=1.;
	    cl_var_last_last=1.;
	  }

	  /** (c.b) loop over k. For each k, if the option of stopping the transfer function computation at some k_max is not selected or if we are below k_max, compute \f$ \Delta_l(k) \f$ with transfer_integrate(); if it is selected and we are above k_max, set the transfer function to zero; if it is selected and we are below k_max, check wether k_max is being reached. */

	  for (index_k = 0; index_k < ptr->k_size[index_mode]; index_k++) {

	    current_k = ptr->k[index_mode][index_k];
	    if (current_k == 0.) {
	      sprintf(Transmit_Error_Message,"%s(L:%d) : k=0, stop to avoid division by zero",__func__,__LINE__);
	      sprintf(ptr->error_message,Transmit_Error_Message);
	      return _FAILURE_;
	    }

	    if (ptr->transfer_verbose > 2)
	      printf("Compute transfer for l=%d k=%e type=%d\n",ptr->l[index_mode][index_l],current_k,index_type);

            /* update previous transfer values in the tc_osc method */
	    if (ppr->transfer_cut == tc_osc) {
	      transfer_last_last = transfer_last;
	      transfer_last = transfer_function;
	    }

            /* update previous relative C_l's variation in the tc_cl method */
	    if (ppr->transfer_cut == tc_cl) {
	      cl_var_last_last = cl_var_last;
	      cl_var_last = cl_var;
	    }

	    /* compute transfer function or set it to zero if above k_max */
	    if ((ppr->transfer_cut == tc_none) || (cut_transfer == _FALSE_)) {
	      if (transfer_integrate(
				     index_mode,
				     index_ic,
				     index_type,
				     index_l,
				     index_k,
				     interpolated_sources,
				     &ti,
				     &transfer_function) == _FAILURE_) {
		sprintf(Transmit_Error_Message,"%s(L:%d) : error callin transfer_integrate()\n=>%s",__func__,__LINE__,ptr->error_message);
		sprintf(ptr->error_message,Transmit_Error_Message);
		return _FAILURE_;
	      }

	    }
	    else {
	      transfer_function = 0.;
	    }
	      
	    /* store transfer function in transfer structure */
	    ptr->transfer[index_mode][((index_ic * ppt->tp_size + index_type)
				       * ptr->l_size[index_mode] + index_l)
				      * ptr->k_size[index_mode] + index_k]
	      = transfer_function;

	    /* in the tc_osc case, update various quantities and check wether k_max is reached */
	    if ((ppr->transfer_cut == tc_osc) && (index_k>=2)) {

	      /* detect local/global maximum of \f$ \Delta_l(k) \f$; if detected, check the criteria for reaching k_max */
	      if ((transfer_last > 0.) && (transfer_function < transfer_last) && (transfer_last > transfer_last_last)) {
		last_local_max = transfer_last;
		if (last_local_max > global_max) {
		  global_max = last_local_max;
		}
		if ((last_local_max-last_local_min) < ppr->transfer_cut_threshold_osc * (global_max-global_min)) {
		  cut_transfer = _TRUE_;
		}

	      }
  
	      /* detect local/global minimum of \f$ \Delta_l(k) \f$; if detected, check the criteria for reaching k_max */ 
	      if ((transfer_last < 0.) && (transfer_function > transfer_last) && (transfer_last < transfer_last_last)) {
		last_local_min = transfer_last;
		if (last_local_min < global_min) {
		  global_min = last_local_min;
		}
		if ((last_local_max-last_local_min) < ppr->transfer_cut_threshold_osc * (global_max-global_min)) {
		  cut_transfer = _TRUE_;
		}
	      }
  	    }

	    /* in the _TC_CUT_ case, update various quantities and check wether k_max is reached */
	    if ((ppr->transfer_cut == tc_cl) && (index_k>=2) && (index_k<ptr->k_size[index_mode]-1) && (transfer_function != 0.)) {

              /* rough estimate of the contribution of the last step to C_l, assuming flat primordial spectrum */
	      delta_cl = transfer_function * transfer_function / current_k * 0.5 * (ptr->k[index_mode][index_k+1] - ptr->k[index_mode][index_k-1]);

	      /* update C_l */
	      cl += delta_cl;

	      /* compute its relative variation */
	      if (cl != 0) {
		cl_var = delta_cl / cl;
	      }
	      else {
		sprintf(ptr->error_message,"%s(L:%d) : cl=0, stop to avoid division by zero",__func__,__LINE__);
		return _FAILURE_;
	      }

	      /* check if k_max is reached */
	      if ((cl_var < cl_var_last) && (cl_var_last > cl_var_last_last)) {
		if (cl_var_last < ppr->transfer_cut_threshold_cl) {
		  cut_transfer = _TRUE_;
		}

	      }
	    }

	    /* end of transfer function computation for given (l,k) */

	  }

	  /* end of loop over k */

	}

	/* end of loop over l */

      }     
      
      /* end of loop over type */

    }

    /* end of loop over initial condition */

    free(interpolated_sources);
    free(source_spline);

  }
  
  /* end of loop over mode */

  free(ti.trans_int);

  return _SUCCESS_;
}

/**
 * Free all memory space allocated by transfer_init().
 *
 * To be called at the end of each run, only when no further calls to
 * transfer_functions_at_k() are needed.
 *
 * @return the error status
 */
int transfer_free() {

  int index_mode;

  for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {
    free(ptr->l[index_mode]);
    free(ptr->k[index_mode]);
    free(ptr->transfer[index_mode]);
  }  

  free(ptr->l_size);
  free(ptr->l);
  free(ptr->k_size);
  free(ptr->k);
  free(ptr->transfer);

  return _SUCCESS_;
  
}

/**
 * Initialize all indices and allocate all arrays in the transfers structure. 
 *
 * Compute list of (k, l) values, allocate and fill corresponding
 * arrays in the transfers structure. Allocate the array of transfer
 * function tables.
 *
 * @return the error status
 */
int transfer_indices_of_transfers() {

  /** Summary: */

  /** - define local variables */

  int index_mode,index_eta;

  /** - allocate arrays of (k, l) values and transfer functions */

  /* number of l values for each mode, l_size[index_mode] */

  ptr->l_size = malloc(ppt->md_size * sizeof(int));
  if (ptr->l_size==NULL) {
    sprintf(ptr->error_message,"%s(L:%d): Cannot allocate ptr->l_size",__func__,__LINE__);
    return _FAILURE_;
  }

  /* list of la values for each mode, l[index_mode] */

  ptr->l = malloc(ppt->md_size * sizeof(int *));
  if (ptr->l==NULL) {
    sprintf(ptr->error_message,"%s(L:%d): Cannot allocate ptr->l",__func__,__LINE__);
    return _FAILURE_;
  }

  /* number of k values for each mode, k_size[index_mode] */

  ptr->k_size = malloc(ppt->md_size * sizeof(int));
  if (ptr->k_size==NULL) {
    sprintf(ptr->error_message,"%s(L:%d): Cannot allocate ptr->k_size",__func__,__LINE__);
    return _FAILURE_;
  }

  /* list of k values for each mode, k[index_mode] */

  ptr->k = malloc(ppt->md_size * sizeof(double *));
  if (ptr->k==NULL) {
    sprintf(ptr->error_message,"%s(L:%d): Cannot allocate ptr->k",__func__,__LINE__);
    return _FAILURE_;
  }

  /* array (of array) of transfer functions for each mode, transfer[index_mode] */

  ptr->transfer = malloc(ppt->md_size * sizeof(double *));
  if (ptr->transfer==NULL) {
    sprintf(ptr->error_message,"%s(L:%d): Cannot allocate ptr->transfer",__func__,__LINE__);
    return _FAILURE_;
  }

  /** - loop over modes (scalar, etc). For each mode: */

  for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {

    /** (a) get number of k values using transfer_get_k_list_size() */
    if (transfer_get_k_list_size(index_mode,&(ptr->k_size[index_mode])) == _FAILURE_) {
      sprintf(Transmit_Error_Message,"%s(L:%d) : error in transfert_get_k_list_size\n=>%s",__func__,__LINE__,ptr->error_message);
      sprintf(ptr->error_message,Transmit_Error_Message);
      return _FAILURE_;
    }

    /** (b) get list of k values using transfer_get_k_list() */
    ptr->k[index_mode] = malloc(ptr->k_size[index_mode]*sizeof(double));
    if (ptr->k[index_mode]==NULL) {
      sprintf(ptr->error_message,"%s(L:%d): Cannot allocate ptr->k[index_mode]",__func__,__LINE__);
      return _FAILURE_;
    }
    if (transfer_get_k_list(index_mode,ptr->k[index_mode]) == _FAILURE_) {
      sprintf(Transmit_Error_Message,"%s(L:%d) : error in transfert_get_k_list()\n=>%s",__func__,__LINE__,ptr->error_message);
      sprintf(ptr->error_message,Transmit_Error_Message);
      return _FAILURE_;
    }

    /** (c) get number of l values using transfer_get_l_list_size() */
    if (transfer_get_l_list_size(index_mode,&(ptr->l_size[index_mode])) == _FAILURE_) {
      sprintf(Transmit_Error_Message,"%s(L:%d) : error in transfer_get_l_list_size()\n=>%s",__func__,__LINE__,ptr->error_message);
      sprintf(ptr->error_message,Transmit_Error_Message);
      return _FAILURE_;
    }

    /** (d) get list of l values using transfer_get_l_list() */
    ptr->l[index_mode] = malloc(ptr->l_size[index_mode]*sizeof(int));
    if (ptr->l[index_mode]==NULL) {
      sprintf(ptr->error_message,"%s(L:%d): Cannot allocate ptr->l[index_mode]",__func__,__LINE__);
      return _FAILURE_;
    }
    if (transfer_get_l_list(index_mode,ptr->l[index_mode]) == _FAILURE_) {
      sprintf(Transmit_Error_Message,"%s(L:%d) : error in transfer_get_l_list()\n=>%s",__func__,__LINE__,ptr->error_message);
      sprintf(ptr->error_message,Transmit_Error_Message);
      return _FAILURE_;
    }

    /** (e) allocate arrays of transfer functions, (ptr->transfer[index_mode])[index_ic][index_type][index_l][index_k] */
    ptr->transfer[index_mode] = malloc(ppt->ic_size[index_mode] * ppt->tp_size * ptr->l_size[index_mode] * ptr->k_size[index_mode] * sizeof(double));
    if (ptr->transfer[index_mode]==NULL) {
      sprintf(ptr->error_message,"%s(L:%d): Cannot allocate ptr->transfer[index_mode]",__func__,__LINE__);
      return _FAILURE_;
    }

  }

  return _SUCCESS_;

}

/**
 * Define number of mutipoles l
 *
 * Define the number of multipoles l for each mode, using information
 * in the precision_params structure.
 *
 * @param index_mode Input: index of requested mode (scalar, tensor, etc) 
 * @param pl_list_size Output: number of multipole 
 * @return the error status
 */
int transfer_get_l_list_size(
			     int index_mode,
			     int * pl_list_size
			     ) {

  int index_l;

  if (ppt->has_scalars && index_mode == ppt->index_md_scalars) {

    if (ppr->l_scalar_max > pbs->l[pbs->l_size-1]) {
      sprintf(Transmit_Error_Message,"%s(L:%d) : For scalar transfer functions, asked for l_max greater than in Bessel table",__func__,__LINE__);
      sprintf(ptr->error_message,Transmit_Error_Message);
      return _FAILURE_;
    }

    index_l=0;
    while((index_l < pbs->l_size-1) && (pbs->l[index_l] <= ppr->l_scalar_max)) {
      index_l++;
    }
    if ((index_l == (pbs->l_size-2)) && (pbs->l[pbs->l_size-1] <= ppr->l_scalar_max)) {
      index_l++;
    }

    *pl_list_size = index_l+1;
     
  }

  if (ppt->has_tensors && index_mode == ppt->index_md_tensors) {

    if (ppr->l_tensor_max > pbs->l[pbs->l_size-1]) {
      sprintf(Transmit_Error_Message,"%s(L:%d) : For transfer transfer functions, asked for l_max greater than in Bessel table",__func__,__LINE__);
      sprintf(ptr->error_message,Transmit_Error_Message);
      return _FAILURE_;
    }

    index_l=0;
    while((index_l < pbs->l_size-1) && (pbs->l[index_l] <= ppr->l_tensor_max)) {
      index_l++;
    }
    if ((index_l == (pbs->l_size-2)) && (pbs->l[pbs->l_size-1] <= ppr->l_tensor_max)) {
      index_l++;
    }

    *pl_list_size = index_l+1;
      
  }
      
}

/**
 * Define list of mutipoles l
 *
 * Define the list of multipoles l for each mode, using information
 * in the precision_params structure.
 *
 * @param index_mode Input: index of requested mode (scalar, tensor, etc) 
 * @param pl_list Output: list of multipole 
 * @return the error status
 */
int transfer_get_l_list(
			int index_mode,
			int * pl_list
			) {

  int index_l;
  
  for (index_l=0; index_l < ptr->l_size[index_mode]; index_l++) {
    pl_list[index_l]=pbs->l[index_l];
  }

  return _SUCCESS_;
  
}

/**
 * Define number of wavenumbers k
 *
 * Define the number of wavenumbers k for each mode, using information
 * in the precision_params and perturbation structure.
 *
 * @param index_mode Input: index of requested mode (scalar, tensor, etc) 
 * @param pk_list_size Output: number of wavenumbers 
 * @return the error status
 */
int transfer_get_k_list_size(
			     int index_mode,
			     int * pk_list_size
			     ) {

  double k_min;
  double k_max_pt;
  double k_step;

  if (ppt->has_scalars && index_mode == ppt->index_md_scalars) {
    k_min = ppt->k[ppt->index_md_scalars][0]; /* first value, inferred from perturbations structure */
    k_max_pt = ppt->k[ppt->index_md_scalars][ppt->k_size[ppt->index_md_scalars]-1]; /* last value, inferred from perturbations structure */
    if ((eta0-eta_rec) != 0.) {
      k_step = 2.*_PI_/(eta0-eta_rec)*ppr->k_step_trans; /* step_size, inferred from precision_params structure */
    }
    else {
      sprintf(ptr->error_message,"%s(L:%d) : (eta0-eta_rec)=0, stop to avoid division by zero",__func__,__LINE__);
      return _FAILURE_;
    }

    if (k_step != 0.) {
      *pk_list_size = (int)((k_max_pt-k_min)/k_step) + 1; /* corresponding number of k values */
          }
    else {
      sprintf(ptr->error_message,"%s(L:%d) : k_step=0, stop to avoid division by zero",__func__,__LINE__);
      return _FAILURE_;
    }

  }

  if (ppt->has_tensors && index_mode == ppt->index_md_tensors)
    *pk_list_size = ppr->k_tensor_number;
  
  return _SUCCESS_;

}

/**
 * Define list of wavenumbers k
 *
 * Define the list of wavenumbers k for each mode, using information
 * in the precision_params and perturbation structures.
 *
 * @param index_mode Input: index of requested mode (scalar, tensor, etc) 
 * @param pk_list Output: list of wavenumbers 
 * @return the error status
 */
int transfer_get_k_list(
			int index_mode,
			double * pk_list
			) {

  double k_min;
  double k_step;
  int index_k;

  if (ppt->has_scalars && index_mode == ppt->index_md_scalars) {

    k_min = ppt->k[ppt->index_md_scalars][0]; /* first value, inferred from perturbations structure */
    if ((eta0-eta_rec) != 0.) {
      k_step = 2.*_PI_/(eta0-eta_rec)*ppr->k_step_trans; /* step_size, inferred from precision_params structure */
    }
    else {
      sprintf(ptr->error_message,"%s(L:%d) : (eta0-eta_rec)=0, stop to avoid division by zero",__func__,__LINE__);
      return _FAILURE_;
    }

    for (index_k = 0; index_k < ptr->k_size[ppt->index_md_scalars]; index_k++) {
      pk_list[index_k] = k_min + index_k * k_step;
    }

    /* k_max = pk_list[ptr->k_size[ppt->index_md_scalars]-1]; */

  }

  if (ppt->has_tensors && index_mode == ppt->index_md_tensors) {
    for (index_k = 0; index_k < ppr->k_tensor_number; index_k++) {
      /* tensor case to be written later */
    }
  }
  
  return _SUCCESS_;
  
}

/**
 * Interpolate sources \f$ S(k, \eta) \f$ using array_spline_table_columns() to get them at the right values of k; get also the second derivative of all source values with respect to \f$ \eta \f$, using again array_spline_table_columns(), in view of later "spline interpolation" of the sources.
 *
 * @param current_index_mode Input : index of mode
 * @param current_index_ic Input : index of initial condition
 * @param current_index_type Input : index of type
 * @param current_index_l Input : index of multipole
 * @param source_spline Input : array of second derivative of sources (filled here but allocated in transfer_init())
 * @param interpolated_sources Output : array of interpolated sources
 * @return the error status
 */
int transfer_interpolate_sources(
				 int current_index_mode,
				 int current_index_ic,
				 int current_index_type,
				 int current_index_l,
				 double * source_spline,
				 double * interpolated_sources) {

  /** Summary: */

  /** - define local variables */

  /* index running on k values in the original source array */
  int index_k, index_eta;

  /* index running on k values in the interpolated source array */
  int index_k_tr;

  /* variables used for spline interpolation algorithm */
  double h, a, b;

  if (array_spline_table_columns(ppt->k[current_index_mode],
				 ppt->k_size[current_index_mode],
				 ppt->sources[current_index_mode][current_index_ic * ppt->tp_size + current_index_type],
				 ppt->eta_size,
				 source_spline,
				 _SPLINE_EST_DERIV_,
				 errmsg
				 ) == _FAILURE_) {
    sprintf(ptr->error_message,"%s(L:%d) : error in array_spline_table_columns()\n=>%s",__func__,__LINE__,errmsg);
    return _FAILURE_;
  }

  /** - interpolate at each k value */

  index_k = 0;
  h = ppt->k[current_index_mode][index_k+1] - ppt->k[current_index_mode][index_k];

  for (index_k_tr = 0; index_k_tr < ptr->k_size[current_index_mode]; index_k_tr++) {

    while (((index_k+1) < ppt->k_size[current_index_mode]) &&
	   (ppt->k[current_index_mode][index_k+1] < 
	    ptr->k[current_index_mode][index_k_tr])) {
      index_k++;
      h = ppt->k[current_index_mode][index_k+1] - ppt->k[current_index_mode][index_k];
    }

    if (h != 0.) {
      b = (ptr->k[current_index_mode][index_k_tr] - ppt->k[current_index_mode][index_k])/h;
    }
    else {
      sprintf(ptr->error_message,"%s(L:%d) : h=0, stop to avoid division by zero",__func__,__LINE__);
      return _FAILURE_;
    }
    a = 1.-b;
    
    for (index_eta = 0; index_eta < ppt->eta_size; index_eta++) {

      interpolated_sources[index_k_tr*ppt->eta_size+index_eta] = 
	a * ppt->sources[current_index_mode]
	[current_index_ic * ppt->tp_size + current_index_type]
	[index_eta*ppt->k_size[current_index_mode]+index_k]
	+ b * ppt->sources[current_index_mode]
	[current_index_ic * ppt->tp_size + current_index_type]
	[index_eta*ppt->k_size[current_index_mode]+index_k+1]
	+ ((a*a*a-a) * source_spline[index_eta*ppt->k_size[current_index_mode]+index_k]
	   +(b*b*b-b) * source_spline[index_eta*ppt->k_size[current_index_mode]+index_k+1])*h*h/6.0;

    }

  }

  /** - fill array of second derivatives with respect to \f$ \eta \f$ (in view of spline interpolation at arbitrary value of \f$ \eta \f$) **/
/*   if (array_spline_table_columns(ppt->eta_sampling, */
/* 				 ppt->eta_size, */
/* 				 interpolated_sources, */
/* 				 ptr->k_size[current_index_mode], */
/* 				 splined_interpolated_sources, */
/* 				 _SPLINE_EST_DERIV_, */
/* 				 errmsg */
/* 				 ) == _FAILURE_) { */
/*     sprintf(ptr->error_message,"problem in array_spline_table_columns \n=>%s",errmsg); */
/*     return _FAILURE_; */
/*   } */

  return _SUCCESS_;

}

/**
 * for given values of l and k, compute the transfer function \f$ \Delta_l(k) \f$ by integrating \f$ S(k, \eta) j_l(k[\eta_0-\eta]) \f$ over \f$ \eta \f$ (in the range of \f$ \eta \f$ values where both functions are non-negligible and actually sampled).
 * 
 * @param current_index_mode Input : index of mode
 * @param current_index_ic Input : index of initial condition
 * @param current_index_type Input : index of type
 * @param current_index_l Input : index of multipole
 * @param current_index_k Input : index of wavenumber
 * @param interpolated_sources Input: array of interpolated sources
 * @param pti Input: pointer towards array of transfer integrand (already allocated and filled with eta values)
 * @param trsf Output: transfer function \f$ \Delta_l(k) \f$ 
 * @return the error status
 */
int transfer_integrate(
		       int current_index_mode,
		       int current_index_ic,
		       int current_index_type,
		       int current_index_l,
		       int current_index_k,
		       double * interpolated_sources,
		       struct transfer_integrand * pti,
		       double * trsf
		       ) {

  /** - define local variables */

  /* maximum value of \f$ \eta \f$ at which \f$ j_l(k[\eta_0-\eta]) \f$ is known, given that \f$ j_l(x) \f$ is sampled above some finite value \f$ x_{\min} \f$ (below which it can be approximated by zero) */  
  double eta_max_bessel;

  /* running value of bessel function */
  double bessel;

  /* index in the source's eta list corresponding to the last point in the overlapping region between sources and bessels */
  int index_eta,index_eta_max;

  /** - find maximum value of \f$ \eta \f$ at which \f$ j_l(k[\eta_0-\eta]) \f$ is known, given that \f$ j_l(x) \f$ is sampled above some finite value \f$ x_{\min} \f$ (below which it can be approximated by zero) */  
  eta_max_bessel = eta0 - pbs->x_min[current_index_l]/ptr->k[current_index_mode][current_index_k]; /* segmentation fault impossible, checked before that current_k != 0 */

  /** - if there is no overlap between the region in which bessels and sources are non-zero, return zero */
  if (eta_max_bessel <= ppt->eta_sampling[0]) {
    *trsf = 0;
    return _SUCCESS_;
  }

  /** - if there is an overlap: */

  /** (a) find index in the source's eta list corresponding to the last point in the overlapping region */ 
  index_eta_max = ppt->eta_size-1;
  while ((ppt->eta_sampling[index_eta_max] > eta_max_bessel) && (index_eta_max > 2))
    index_eta_max--;

  /** (b) the source function can vanish at large $\f k \eta \f$. Check if further points can be eliminated. */
  while ((interpolated_sources[current_index_k * ppt->eta_size + index_eta_max-1] == 0.)  && (index_eta_max > 2)) 
    index_eta_max--;

  /** (c) loop over points: */

  for (index_eta = 0; index_eta <= index_eta_max; index_eta++) {

    if (bessel_at_x(ptr->k[current_index_mode][current_index_k] * (eta0-ppt->eta_sampling[index_eta]),current_index_l,&bessel) == _FAILURE_) {
      sprintf(Transmit_Error_Message,"%s(L:%d) : error when calling bessel_at_x()\n=>%s",__func__,__LINE__,pbs->error_message);
      sprintf(ptr->error_message,Transmit_Error_Message);
      return _FAILURE_;
    }

    pti->trans_int[pti->trans_int_col_num*index_eta+pti->trans_int_y]= 
      interpolated_sources[current_index_k * ppt->eta_size + index_eta]*bessel;

  }

  /** (d) spline the integrand: */

  if (array_spline(pti->trans_int,
		   pti->trans_int_col_num,
		   index_eta_max+1,
		   pti->trans_int_eta,
		   pti->trans_int_y,
		   pti->trans_int_ddy,
		   _SPLINE_EST_DERIV_,
		   errmsg) == _FAILURE_) {
    sprintf(ptr->error_message,"%s(L:%d) : error in array_spline()\n=>%s",__func__,__LINE__,errmsg);
    return _FAILURE_;
  }

  /** (e) integrate: */
  if (array_integrate_all_spline(pti->trans_int,
				 pti->trans_int_col_num,
				 index_eta_max+1,
				 pti->trans_int_eta,
				 pti->trans_int_y,
				 pti->trans_int_ddy,
				 trsf,
				 errmsg) == _FAILURE_) {
    sprintf(ptr->error_message,"%s(L:%d) : error in array_integrate_all_spline()\n=>%s",__func__,__LINE__,errmsg);
    return _FAILURE_;
  }

  /** (f) correct for last piece of integral (up to point where bessel vanishes) */
  *trsf += (eta_max_bessel-ppt->eta_sampling[index_eta_max])
    * pti->trans_int[pti->trans_int_col_num*index_eta_max+pti->trans_int_y]/2.;

  if ((ppt->has_scalars == _TRUE_) && (current_index_mode == ppt->index_md_scalars)) {
    if ((ppt->has_source_p == _TRUE_) && (current_index_type == ppt->index_tp_p)) {
      /* scalar polarization */
      *trsf *= sqrt((pbs->l[current_index_l]+2.) * (pbs->l[current_index_l]+1.) * (pbs->l[current_index_l]) * (pbs->l[current_index_l]-1.)); 
    }
  }
  else {
    if ((ppt->has_tensors == _TRUE_) && (current_index_mode == ppt->index_md_tensors)) {
      /* tensors not coded yet */
    }
  }
  

  return _SUCCESS_;
}
