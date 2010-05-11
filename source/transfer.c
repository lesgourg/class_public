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

#ifdef _OPENMP
#include "omp.h"
#endif
#include "transfer.h"

/** @name - miscellaneous: */

//@{

ErrorMsg Transmit_Error_Message; /**< contains error message */

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
 * @param index_tt Input: index of requested type
 * @param index_l Input: index of requested multipole
 * @param k Input: any wavenumber
 * @param ptransfer_local Output: transfer function
 * @return the error status
 */
int transfer_functions_at_k(
			    struct transfers * ptr,
			    int index_mode,
			    int index_ic,
			    int index_tt,
			    int index_l,
			    double k,
			    double * ptransfer 
			    ) {
  /** Summary: */

  /** - interpolate in pre-computed table using array_interpolate_two() */
  class_call(array_interpolate_two(
				   ptr->k[index_mode],
				   1,
				   0,
				   ptr->transfer[index_mode]
				   +((index_ic * ptr->tt_size + index_tt) * ptr->l_size[index_mode] + index_l)
				   * ptr->k_size[index_mode],
				   1,
				   ptr->k_size[index_mode],
				   k,
				   ptransfer,
				   1,
				   ptr->error_message),
	     ptr->error_message,
	     ptr->error_message);
  
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
 *    (transfer[index_mode])[index_ic][index_tt][index_l][index_k]
 *
 * This function shall be called at the beginning of each run, but
 * only after background_init(), thermodynamics_init() and perturb_init(). It
 * allocates memory spaces which should be freed later with
 * transfer_free().
 *
 * @param pba Input : Initialized background structure 
 * @param pth Input : Initialized thermodynamics structure 
 * @param ppt Input : Initialized perturbation structure
 * @param pbs Input : Initialized bessels structure
 * @param ppr Input : Parameters describing how the computation is to be performed
 * @param ptr Output : Initialized transfers structure
 * @return the error status
 */
int transfer_init(
		  struct precision * ppr,
		  struct background * pba,
		  struct thermo * pth,
		  struct perturbs * ppt,
		  struct bessels * pbs,
		  struct transfers * ptr
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
  int index_tt; 
  /* running index for multipoles */
  int index_l; 
  /* another index */
  int index;
  /* current wavenumber value */
  double current_k;
  /* conformal time today */
  double eta0;
  /* conformal time at recombination */
  double eta_rec;
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
  /* (pointer on) integrand structure */
  struct transfer_integrand *ti;
  /* abort flag, useful for parallel regions */
  int abort;

#ifdef _OPENMP
  /* number of available omp threads */
  int number_of_threads;
  /* table of integrand of transfer function, copied number_of_threads times*/
  struct transfer_integrand *pti;
  /* instrumentation times */
  double tstart, tstop;
#endif

  if ((ppt->has_cl_cmb_temperature == _FALSE_) &&
      (ppt->has_cl_cmb_polarization == _FALSE_) &&
      (ppt->has_cl_cmb_lensing_potential == _FALSE_)) {
    ptr->tt_size=0;
    if (ptr->transfer_verbose > 0)
      printf("No harmonic space transfer functions to compute. Transfer module skipped.\n");
    return _SUCCESS_;
  }

  if (ptr->transfer_verbose > 0)
    printf("Computing transfers\n");

  ptr->md_size = ppt->md_size;

  /** - get conformal age / recombination time from background / thermodynamics structures (only place where these structures are used in this module) */
  eta0 = pba->conformal_age;
  eta_rec = pth->eta_rec;

  /** - initialize all indices in the transfers structure and allocate all its arrays using transfer_indices_of_transfers() */
  class_call(transfer_indices_of_transfers(ppr,ppt,pbs,ptr,eta0,eta_rec),
	     ptr->error_message,
	     ptr->error_message);

  /** - initialize all indices in the transfer_integrand structure and allocate its array. Fill the eta column. */

#ifdef _OPENMP
#pragma omp parallel
  {
    number_of_threads = omp_get_num_threads();
  }
  pti = (struct transfer_integrand*)malloc(number_of_threads * sizeof(struct transfer_integrand));
#endif

  abort=0;
#pragma omp parallel shared(pti) private(ti,index)
  {
#ifdef _OPENMP
    ti = &pti[omp_get_thread_num()];
#else
    ti = (struct transfer_integrand*)malloc(sizeof(struct transfer_integrand));
#endif
    index = 0;
    ti->trans_int_eta = index;
    index++;
    ti->trans_int_y = index;
    index++;
    ti->trans_int_ddy = index;
    index++;
    ti->trans_int_col_num = index;

    class_alloc(ti->trans_int,
		sizeof(double) * ppt->eta_size * ti->trans_int_col_num,
		ptr->error_message);

    if (abort != 0)
      {
	for (index=0; index < ppt->eta_size; index++)
	  ti->trans_int[ti->trans_int_col_num*index+ti->trans_int_eta] = ppt->eta_sampling[index];
      }
  }
  if (abort) return _FAILURE_;

  /** - loop over all indices of the table of transfer functions. For each mode, initial condition and type: */

  for (index_mode = 0; index_mode < ptr->md_size; index_mode++) {

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

      for (index_tt = 0; index_tt < ptr->tt_size; index_tt++) {

	/** (b) interpolate sources to get them at the right values of k using transfer_interpolate_sources() */

	if (transfer_interpolate_sources(ppt,
					 ptr,
					 index_mode,
					 index_ic,
					 index_tt,
					 index_l,
					 eta0,
					 eta_rec,
					 source_spline,
					 interpolated_sources) == _FAILURE_) {
	  sprintf(Transmit_Error_Message,"%s(L:%d) : error in transfer_interpolate_sources()\n=>%s",__func__,__LINE__,ptr->error_message);
	  sprintf(ptr->error_message,"%s",Transmit_Error_Message);
	  return _FAILURE_;
	}

	/** (c) loop over l. For each value of l: */

	/***** THIS IS THE LOOP WHICH SHOULD BE PARALLELISED ******/
	abort=0;
#pragma omp parallel							\
  shared (ptr,ppr,ppt,index_mode,index_ic,index_tt,eta0,eta_rec,	\
	  interpolated_sources,pti,abort)				\
  private (index_l,cut_transfer,global_max,global_min,			\
	   last_local_max,last_local_min,transfer_function,		\
	   transfer_last,transfer_last_last,cl,cl_var,			\
	   cl_var_last,cl_var_last_last,delta_cl,			\
	   index_k,current_k,ti,					\
	   Transmit_Error_Message)

	{
#ifdef _OPENMP
	  ti = &pti[omp_get_thread_num()];
#endif
#ifdef _OPENMP
	tstart = omp_get_wtime();
#endif
#pragma omp for schedule (dynamic)
	  for (index_l = 0; index_l < ptr->l_size[index_mode]; index_l++) {
#pragma omp flush(abort)
	    if (!abort) {
	      
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
		  sprintf(ptr->error_message,"%s",Transmit_Error_Message);
		  abort = 1;
#pragma omp flush (abort)
		}
		
		if (ptr->transfer_verbose > 2)
		  printf("Compute transfer for l=%d k=%e type=%d\n",ptr->l[index_mode][index_l],current_k,index_tt);
		
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
		  if (transfer_integrate(ppt,
					 pbs,
					 ptr,
					 index_mode,
					 index_ic,
					 index_tt,
					 index_l,
					 index_k,
					 eta0,
					 eta_rec,
					 interpolated_sources,
					 ti,
					 &transfer_function) == _FAILURE_) {
		    sprintf(Transmit_Error_Message,"%s(L:%d) : error callin transfer_integrate()\n=>%s",__func__,__LINE__,ptr->error_message);
		    sprintf(ptr->error_message,"%s",Transmit_Error_Message);
		    abort = 1;
#pragma omp flush (abort)
		  }
		  
		}
		else {
		  transfer_function = 0.;
		}
		
		/* store transfer function in transfer structure */
		ptr->transfer[index_mode][((index_ic * ptr->tt_size + index_tt)
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
		    abort = 1;
#pragma omp flush (abort)
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
	  
	  }

	/* end of loop over l */
#ifdef _OPENMP
	  tstop = omp_get_wtime();
	  printf("Time spent in parallel region (loop over ells) = %e\n",tstop-tstart);
#endif
	}

	/* end of parallel region */
	if (abort) return _FAILURE_;

      }     
      
      /* end of loop over type */

    }

    /* end of loop over initial condition */

    free(interpolated_sources);
    free(source_spline);

  }
  
  /* end of loop over mode */
#ifdef _OPENMP
#pragma omp parallel shared(pti) private(ti)
  {
    ti = &pti[omp_get_thread_num()];
    free(ti->trans_int);
  }
  free(pti);
#else
  free(ti->trans_int);
  free(ti);
#endif
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
int transfer_free(
		  struct transfers * ptr
		  ) {

  int index_mode;

  if (ptr->tt_size>0) {
    for (index_mode = 0; index_mode < ptr->md_size; index_mode++) {
      free(ptr->l[index_mode]);
      free(ptr->k[index_mode]);
      free(ptr->transfer[index_mode]);
    }  

    free(ptr->l_size);
    free(ptr->l);
    free(ptr->k_size);
    free(ptr->k);
    free(ptr->transfer);
  }

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
int transfer_indices_of_transfers(
				  struct precision * ppr,
				  struct perturbs * ppt,
				  struct bessels * pbs,
				  struct transfers * ptr,
				  double eta0,
				  double eta_rec
				  ) {

  /** Summary: */

  /** - define local variables */

  int index_mode,index_eta,index_tt;

  /** define indices for transfer types */
  
  index_tt = 0;
  if (ppt->has_cl_cmb_temperature == _TRUE_) {
    ptr->index_tt_t = index_tt;
    index_tt++;
  }
  if (ppt->has_cl_cmb_polarization == _TRUE_) {
    ptr->index_tt_p = index_tt;
    index_tt++;
  }
  if (ppt->has_cl_cmb_lensing_potential == _TRUE_) {
    ptr->index_tt_lcmb = index_tt;
    index_tt++;
  }
  ptr->tt_size=index_tt;

  /** - allocate arrays of (k, l) values and transfer functions */

  /* number of l values for each mode, l_size[index_mode] */

  class_alloc(ptr->l_size,ptr->md_size * sizeof(int),ptr->error_message);

  /* list of la values for each mode, l[index_mode] */

  class_alloc(ptr->l,ptr->md_size * sizeof(int *),ptr->error_message);

  /* number of k values for each mode, k_size[index_mode] */

  class_alloc(ptr->k_size,ptr->md_size * sizeof(int),ptr->error_message);

  /* list of k values for each mode, k[index_mode] */

  class_alloc(ptr->k,ptr->md_size * sizeof(double *),ptr->error_message);

  /* array (of array) of transfer functions for each mode, transfer[index_mode] */

  class_alloc(ptr->transfer,ptr->md_size * sizeof(double *),ptr->error_message);

  /** - loop over modes (scalar, etc). For each mode: */

  for (index_mode = 0; index_mode < ptr->md_size; index_mode++) {

    /** (a) get number of k values using transfer_get_k_list_size() */
    class_call(transfer_get_k_list_size(ppr,ppt,ptr,index_mode,eta0,eta_rec,&(ptr->k_size[index_mode])),
					ptr->error_message,
					ptr->error_message);

    /** (b) get list of k values using transfer_get_k_list() */
    class_alloc(ptr->k[index_mode],ptr->k_size[index_mode]*sizeof(double),ptr->error_message);

    class_call(transfer_get_k_list(ppr,ppt,ptr,index_mode,eta0,eta_rec,ptr->k[index_mode]),
	       ptr->error_message,
	       ptr->error_message);

    /** (c) get number of l values using transfer_get_l_list_size() */
    class_call(transfer_get_l_list_size(ppr,ppt,pbs,ptr,index_mode,&(ptr->l_size[index_mode])),
	       ptr->error_message,
	       ptr->error_message);

    /** (d) get list of l values using transfer_get_l_list() */
    class_alloc(ptr->l[index_mode],ptr->l_size[index_mode]*sizeof(int),ptr->error_message);

    class_call(transfer_get_l_list(ppt,pbs,ptr,index_mode,ptr->l[index_mode]),
	       ptr->error_message,
	       ptr->error_message); 

    /** (e) allocate arrays of transfer functions, (ptr->transfer[index_mode])[index_ic][index_tt][index_l][index_k] */
    class_alloc(ptr->transfer[index_mode],
		ppt->ic_size[index_mode] * ptr->tt_size * ptr->l_size[index_mode] * ptr->k_size[index_mode] * sizeof(double),
		ptr->error_message);

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
			     struct precision * ppr,
			     struct perturbs * ppt,
			     struct bessels * pbs,
			     struct transfers * ptr,
			     int index_mode,
			     int * pl_list_size
			     ) {

  int index_l;

  if (ppt->has_scalars && index_mode == ppt->index_md_scalars) {

    class_test(ppr->l_scalar_max > pbs->l[pbs->l_size-1],
	       ptr->error_message,
	       "For scalar transfer functions, asked for l_max greater than in Bessel table");

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

    class_test(ppr->l_tensor_max > pbs->l[pbs->l_size-1],
	       ptr->error_message,
	       "For tensor transfer functions, asked for l_max greater than in Bessel table");

    index_l=0;
    while((index_l < pbs->l_size-1) && (pbs->l[index_l] <= ppr->l_tensor_max)) {
      index_l++;
    }
    if ((index_l == (pbs->l_size-2)) && (pbs->l[pbs->l_size-1] <= ppr->l_tensor_max)) {
      index_l++;
    }

    *pl_list_size = index_l+1;
      
  }
  return _SUCCESS_;

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
			struct perturbs * ppt,
			struct bessels * pbs,
			struct transfers * ptr,
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
			     struct precision * ppr,
			     struct perturbs * ppt,
			     struct transfers * ptr,
			     int index_mode,
			     double eta0,
			     double eta_rec,
			     int * pk_list_size
			     ) {

  double k_min;
  double k_max_pt;
  double k_step;

  if (ppt->has_scalars && index_mode == ppt->index_md_scalars) {
    k_min = ppt->k[ppt->index_md_scalars][0]; /* first value, inferred from perturbations structure */
    k_max_pt = ppt->k[ppt->index_md_scalars][ppt->k_size_cl[ppt->index_md_scalars]-1]; /* last value, inferred from perturbations structure */

    class_test((eta0-eta_rec) == 0.,
	       ptr->error_message,
	       "stop to avoid division by zero");

    k_step = 2.*_PI_/(eta0-eta_rec)*ppr->k_step_trans; /* step_size, inferred from precision_params structure */

    class_test(k_step == 0.,
	       ptr->error_message,
	       "stop to avoid division by zero");

    *pk_list_size = (int)((k_max_pt-k_min)/k_step) + 1; /* corresponding number of k values */

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
			struct precision * ppr,
			struct perturbs * ppt,
			struct transfers * ptr,
			int index_mode,
			double eta0,
			double eta_rec,
			double * pk_list
			) {

  double k_min;
  double k_step;
  int index_k;

  if (ppt->has_scalars && index_mode == ppt->index_md_scalars) {

    k_min = ppt->k[ppt->index_md_scalars][0]; /* first value, inferred from perturbations structure */

    class_test((eta0-eta_rec) == 0.,
	       ptr->error_message,
	       "stop to avoid division by zero");

    k_step = 2.*_PI_/(eta0-eta_rec)*ppr->k_step_trans; /* step_size, inferred from precision_params structure */

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
 * @param current_index_tt Input : index of type of transfer
 * @param current_index_l Input : index of multipole
 * @param source_spline Input : array of second derivative of sources (filled here but allocated in transfer_init())
 * @param interpolated_sources Output : array of interpolated sources
 * @return the error status
 */
int transfer_interpolate_sources(
				 struct perturbs * ppt,
				 struct transfers * ptr,
				 int current_index_mode,
				 int current_index_ic,
				 int current_index_tt,
				 int current_index_l,
				 double eta0,
				 double eta_rec,
				 double * source_spline,
				 double * interpolated_sources) {

  /** Summary: */

  /** - define local variables */

  /* index running on k values in the original source array */
  int index_k;

  /* index running on time */
  int index_eta;

  /* index running on type of source (not type of transfer) */
  int index_type;

  /* index running on k values in the interpolated source array */
  int index_k_tr;

  /* variables used for spline interpolation algorithm */
  double h, a, b;

  /* which source are we considering? Correspondance between transfer
     type and source type*/

  if ((ppt->has_cl_cmb_temperature == _TRUE_) &&
      (current_index_tt == ptr->index_tt_t)) 
    index_type=ppt->index_tp_t;

  if ((ppt->has_cl_cmb_polarization == _TRUE_) &&
      (current_index_tt == ptr->index_tt_p)) 
    index_type=ppt->index_tp_p;

  if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) &&
      (current_index_tt == ptr->index_tt_lcmb)) 
    index_type=ppt->index_tp_g;

  class_call(array_spline_table_columns(ppt->k[current_index_mode],
					ppt->k_size[current_index_mode],
					ppt->sources[current_index_mode][current_index_ic * ppt->tp_size + index_type],
					ppt->eta_size,
					source_spline,
					_SPLINE_EST_DERIV_,
					ptr->error_message),
	     ptr->error_message,
	     ptr->error_message);

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

    class_test(h == 0.,
	       ptr->error_message,
	       "stop to avoid division by zero");

    b = (ptr->k[current_index_mode][index_k_tr] - ppt->k[current_index_mode][index_k])/h;
    a = 1.-b;
    
    for (index_eta = 0; index_eta < ppt->eta_size; index_eta++) {

      interpolated_sources[index_k_tr*ppt->eta_size+index_eta] = 
	a * ppt->sources[current_index_mode]
	[current_index_ic * ppt->tp_size + index_type]
	[index_eta*ppt->k_size[current_index_mode]+index_k]
	+ b * ppt->sources[current_index_mode]
	[current_index_ic * ppt->tp_size + index_type]
	[index_eta*ppt->k_size[current_index_mode]+index_k+1]
	+ ((a*a*a-a) * source_spline[index_eta*ppt->k_size[current_index_mode]+index_k]
	   +(b*b*b-b) * source_spline[index_eta*ppt->k_size[current_index_mode]+index_k+1])*h*h/6.0;

      /* for lensing, multiply gravitational potential by appropriate window function */

      /* case of cmb lensing */
      if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) &&
	  (current_index_tt == ptr->index_tt_lcmb)) {
	/* lensing source =  4 pi W(eta) psi(k,eta) H(eta-eta_rec) 
	   with 
	   psi = (newtonian) gravitationnal potential  
	   W = 2(eta-eta_rec)/(eta_0-eta)/(eta_0-eta_rec) 
	   H(x) = Heaviside
	   (in eta = eta_0, source = 0 to avoid division by zero (regulated anyway by Bessel)).
	*/
	if ((ppt->eta_sampling[index_eta] > eta_rec) && 
	    ((eta0-ppt->eta_sampling[index_eta]) > 0.)) {
	  interpolated_sources[index_k_tr*ppt->eta_size+index_eta] *=
	    4.*_PI_*(2.*
		     (ppt->eta_sampling[index_eta]-eta_rec)
		     /(eta0-ppt->eta_sampling[index_eta])
		     /(eta0-eta_rec));
	}
	else {
	  interpolated_sources[index_k_tr*ppt->eta_size+index_eta] = 0;
	}

      }

    }

  }

  return _SUCCESS_;

}

/**
 * for given values of l and k, compute the transfer function \f$ \Delta_l(k) \f$ by integrating \f$ S(k, \eta) j_l(k[\eta_0-\eta]) \f$ over \f$ \eta \f$ (in the range of \f$ \eta \f$ values where both functions are non-negligible and actually sampled).
 * 
 * @param current_index_mode Input : index of mode
 * @param current_index_ic Input : index of initial condition
 * @param current_index_tt Input : index of type
 * @param current_index_l Input : index of multipole
 * @param current_index_k Input : index of wavenumber
 * @param interpolated_sources Input: array of interpolated sources
 * @param pti Input: pointer towards array of transfer integrand (already allocated and filled with eta values)
 * @param trsf Output: transfer function \f$ \Delta_l(k) \f$ 
 * @return the error status
 */
int transfer_integrate(
		       struct perturbs * ppt,
		       struct bessels * pbs,
		       struct transfers * ptr,
		       int current_index_mode,
		       int current_index_ic,
		       int current_index_tt,
		       int current_index_l,
		       int current_index_k,
		       double eta0,
		       double eta_rec,
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

    class_call(bessel_at_x(ptr->k[current_index_mode][current_index_k] * (eta0-ppt->eta_sampling[index_eta]),current_index_l,&bessel),
	       ptr->error_message,
	       ptr->error_message);

    pti->trans_int[pti->trans_int_col_num*index_eta+pti->trans_int_y]= 
      interpolated_sources[current_index_k * ppt->eta_size + index_eta]*bessel;

  }

  /** (d) spline the integrand: */

  class_call(array_spline(pti->trans_int,
			  pti->trans_int_col_num,
			  index_eta_max+1,
			  pti->trans_int_eta,
			  pti->trans_int_y,
			  pti->trans_int_ddy,
			  _SPLINE_EST_DERIV_,
			  ptr->error_message),
	     ptr->error_message,
	     ptr->error_message);

  /** (e) integrate: */
  class_call(array_integrate_all_spline(pti->trans_int,
					pti->trans_int_col_num,
					index_eta_max+1,
					pti->trans_int_eta,
					pti->trans_int_y,
					pti->trans_int_ddy,
					trsf,
					ptr->error_message),
	     ptr->error_message,
	     ptr->error_message);

  /** (f) correct for last piece of integral (up to point where bessel vanishes) */
  *trsf += (eta_max_bessel-ppt->eta_sampling[index_eta_max])
    * pti->trans_int[pti->trans_int_col_num*index_eta_max+pti->trans_int_y]/2.;

  if ((ppt->has_scalars == _TRUE_) && (current_index_mode == ppt->index_md_scalars)) {
    if ((ppt->has_cl_cmb_polarization == _TRUE_) && (current_index_tt == ptr->index_tt_p)) {
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
