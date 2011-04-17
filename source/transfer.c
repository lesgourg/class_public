/** @file transfer.c Documented transfer module.
 *
 * Julien Lesgourgues, 26.08.2010    
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
 * -# transfer_init() at the beginning (but after perturb_init() 
 *    and bessel_init())
 *
 * -# transfer_functions_at_k() at any later time 
 *
 * -# transfer_free() at the end, when no more calls to 
 *    transfer_functions_at_k() are needed
 * 
 * Note that in the standard implementation of CLASS, only the pre-computed 
 * values of the transfer functions are used, no interpolation is necessary; 
 * hence the routine transfer_functions_at_k() is actually never called.
 */

#include "transfer.h"

/** 
 * Transfer function \f$ \Delta_l^{X} (k) \f$ at a given wavenumber k.
 *
 * For a given mode (scalar, vector, tensor), initial condition, type
 * (temperature, polarization, lensing, etc) and multipole, compute
 * the transfer function for an arbitary value of k by interpolating
 * between pre-computed values of k. This
 * function can be called from whatever module at whatever time,
 * provided that transfer_init() has been called before, and
 * transfer_free() has not been called yet.
 *
 * @param index_mode Input: index of requested mode
 * @param index_ic   Input: index of requested initial condition
 * @param index_tt   Input: index of requested type
 * @param index_l    Input: index of requested multipole
 * @param k          Input: any wavenumber
 * @param transfer_function Output: transfer function
 * @return the error status
 */
int transfer_functions_at_k(
			    struct transfers * ptr,
			    int index_mode,
			    int index_ic,
			    int index_tt,
			    int index_l,
			    double k,
			    double * transfer_function
			    ) {
  /** Summary: */

  /** - interpolate in pre-computed table using array_interpolate_two() */
  class_call(array_interpolate_two(
				   ptr->k[index_mode],
				   1,
				   0,
				   ptr->transfer[index_mode]
				   +((index_ic * ptr->tt_size[index_mode] + index_tt) * ptr->l_size[index_mode] + index_l)
				   * ptr->k_size[index_mode],
				   1,
				   ptr->k_size[index_mode],
				   k,
				   transfer_function,
				   1,
				   ptr->error_message),
	     ptr->error_message,
	     ptr->error_message);
  
  return _SUCCESS_;
}

/**
 * This routine initializes the transfers structure, (in particular,
 * computes table of transfer functions \f$ \Delta_l^{X} (k) \f$)
 *
 * Main steps: 
 *
 * - initialize all indices in the transfers structure
 *   and allocate all its arrays using transfer_indices_of_transfers().
 *
 * - for a all requested modes (scalar, vector, tensor),
 *   initial conditions and types (temperature, polarization, lensing,
 *   etc), compute the transfer function \f$ \Delta_l^{X} (k) \f$
 *   following the following steps:
 &
 * -# interpolate sources \f$ S(k, \tau) \f$ to get them at the right 
 *    values of k using transfer_interpolate_sources()
 *
 * -# for each l, compute the transfer function by convolving the 
 *    sources with the Bessel functions using transfer_compute_for_each_l()
 *    (this step is parallelized). Store result in the transfer table 
 *    transfer[index_mode][((index_ic * ptr->tt_size[index_mode] + index_tt) * ptr->l_size[index_mode] + index_l) * ptr->k_size[index_mode] + index_k]
 *
 * @param ppr Input : pointer to precision structure 
 * @param pba Input : pointer to background structure 
 * @param pth Input : pointer to thermodynamics structure 
 * @param ppt Input : pointer to perturbation structure
 * @param pbs Input : pointer to bessels structure
 * @param ptr Output: pointer to initialized transfers structure
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
  /* running index for types */
  int index_tt; 
  /* running index for multipoles */
  int index_l; 
  /* running index for conformal time */
  int index_tau;

  /* conformal time today */
  double tau0;
  /* conformal time at recombination */
  double tau_rec;

  /* array of source derivatives S''(k,tau) 
     (second derivative with respect to k, not tau!), 
     used to interpolate sources at the right values of k,
     source_spline[index_tau*ppt->k_size[index_mode]+index_k] */
  double * source_spline;

  /* table of source functions interpolated at the right values of k, 
     interpolated_sources[index_k_tr*ppt->tau_size+index_tau] */
  double * interpolated_sources;

  /* we deal with workspaces which a contiguous memory zone containing various
     fields used by the integration routine */

  /* - pointer used to assign adresses to the various workspace fields */
  double * address_in_workspace;

  /* - first workspace field: list of tau0-tau values, tau0_minus_tau[index_tau] */
  double * tau0_minus_tau;

  /* - second workspace field: list of delta_tau values, delta_tau[index_tau] */
  double * delta_tau;

  /* - third workspace field, identical to above interpolated sources:
     sources[index_k_tr*ppt->tau_size+index_tau] */
  double * sources;

  /* - fourth workspace field, containing just a double: value of x at
     which bessel functions become non-negligible for a given l (
     infered from bessel module) */
  double * x_min_l;

  /* - fifth workspace field containing the list of j_l(x) values */
  double * j_l;

  /* - sixth workspace field containing the list of j_l''(x) values */
  double * ddj_l;

  /* - optional: seventh workspace field containing the list of j_l(x) values */
  double * dj_l;

  /* - optional: eigth workspace field containing the list of j_l''(x) values */
  double * dddj_l;

  /* no more workspace fields */

  /* pointer on one workspace per thread (only one if no openmp) */
  double ** pw;

  /* number of available omp threads (remains always one if no openmp) */
  int number_of_threads=1;

  /* index of each thread (reminas always zero if no openmp) */
  int thread=0;

  /* number of values of x for a given l infered from bessel module */
  int x_size_l;

  /** - number of functions stored in bessel structure: at least j_l,
      j_l'', plus eventually j_l', j_l''' */
  int num_j;

  /* This code can be optionally compiled with the openmp option for parallel computation.
     Inside parallel regions, the use of the command "return" is forbidden.
     For error management, instead of "return _FAILURE_", we will set the variable below
     to "abort = _TRUE_". This will lead to a "return _FAILURE_" jus after leaving the 
     parallel region. */
  int abort;

#ifdef _OPENMP

  /* instrumentation times */
  double tstart, tstop, tspent;

#endif

  /** check whether any spectrum in harmonic space (i.e., any C_l's) is actually requested */

  if (ppt->has_cls == _FALSE_) {
    ptr->has_cls = _FALSE_;
    if (ptr->transfer_verbose > 0)
      printf("No harmonic space transfer functions to compute. Transfer module skipped.\n");
    return _SUCCESS_;
  }
  else
    ptr->has_cls = _TRUE_;

  if (ptr->transfer_verbose > 0)
    printf("Computing transfers\n");

  /** get number of modes (scalars, tensors...) */

  ptr->md_size = ppt->md_size;

  /** - get conformal age / recombination time 
        from background / thermodynamics structures 
	(only place where these structures are used in this module) */
  tau0 = pba->conformal_age;
  tau_rec = pth->tau_rec;

  /** - initialize all indices in the transfers structure and 
        allocate all its arrays using transfer_indices_of_transfers() */
  class_call(transfer_indices_of_transfers(ppr,ppt,pbs,ptr,pth->rs_rec),
	     ptr->error_message,
	     ptr->error_message);

  /** - how many functions stored in bessel structure: just j_l, j_l'', or also j_l', j_l''' ? */
  if (ppr->transfer_cut == tc_env) {
    num_j = 4;
  }
  else {
    num_j = 2;
  }

  /* find number of threads */

#ifdef _OPENMP
#pragma omp parallel
  {
    number_of_threads = omp_get_num_threads();
  }

#endif

  /* allocate the pointer to one workspace per thread */
  class_alloc(pw,number_of_threads*sizeof(double*),ptr->error_message);

  /** - loop over all modes. For each mode: */ 

  for (index_mode = 0; index_mode < ptr->md_size; index_mode++) {

    /** (a) allocate temporary arrays relevant for this mode */

    /** (a.1.) array of sources interpolated at correct values of k */

    class_alloc(interpolated_sources,
		ptr->k_size[index_mode]*ppt->tau_size*sizeof(double),
		ptr->error_message);
    
    /** (a.2.) second derivative of sources in the previous source arrays, useful to
       obtain the one above */

    class_alloc(source_spline,
		ppt->k_size[index_mode]*ppt->tau_size*sizeof(double),
		ptr->error_message);

    /** (a.3.) workspace, allocated in a parallel zone since in openmp
       version there is one workspace per thread */
    
    /* initialize error management flag */
    abort = _FALSE_;
    
    /* beginning of parallel region */
    
#pragma omp parallel							\
  shared(ptr,index_mode,ppt,pbs,pw)					\
  private(thread,address_in_workspace,tau0_minus_tau,delta_tau,index_tau)
    {
      
#ifdef _OPENMP
      thread = omp_get_thread_num();
#endif


      class_alloc_parallel(pw[thread],
			   ((2+ptr->k_size[index_mode])*ppt->tau_size
			    +1+num_j*pbs->x_size_max)*
			   sizeof(double),
			   ptr->error_message);
      
      /* -- define the address of the first two fields in the workspace, tau0_minus_tau and delta_tau */

      address_in_workspace = pw[thread];
	  
      tau0_minus_tau = address_in_workspace;
      address_in_workspace += ppt->tau_size;
      
      delta_tau  = address_in_workspace;

      /* -- fill these two fields since their content does not depend
	 on ic, type and l */

      for (index_tau=0; index_tau < ppt->tau_size; index_tau++) {
	tau0_minus_tau[index_tau] = tau0 - ppt->tau_sampling[index_tau];
      }

      delta_tau[0] = ppt->tau_sampling[1]-ppt->tau_sampling[0];
      
      for (index_tau=1; index_tau < ppt->tau_size-1; index_tau++)
	delta_tau[index_tau] = ppt->tau_sampling[index_tau+1]-ppt->tau_sampling[index_tau-1];
      
      delta_tau[ppt->tau_size-1] = ppt->tau_sampling[ppt->tau_size-1]-ppt->tau_sampling[ppt->tau_size-2];
      
    } /* end of parallel region */
    
    if (abort == _TRUE_) return _FAILURE_;
    
    /** (b) now loop over initial conditions and types: For each of them: */

    for (index_ic = 0; index_ic < ppt->ic_size[index_mode]; index_ic++) {

      for (index_tt = 0; index_tt < ptr->tt_size[index_mode]; index_tt++) {
	
        /** (b.1) interpolate sources to get them at the right values of k 
                using transfer_interpolate_sources() */

	if (ptr->transfer_verbose>2)
	  printf("In %s: Interpolate sources for one mode/ic/type.\n",
		 __func__);

	class_call(transfer_interpolate_sources(ppt,
						ptr,
						tau0,
						tau_rec,
						index_mode,
						index_ic,
						index_tt,
						source_spline,
						interpolated_sources),
		   ptr->error_message,
		   ptr->error_message);

	/** (b.2) store the sources in the workspace and define all
	    fields in this workspace (in a parallel zone, since in
	    open mp version there is one workspace per thread). The
	    (parallelized) loop over l values will take place right
	    after that in the same parallel zone. */

#ifdef _OPENMP
	if (ptr->transfer_verbose>1)
	  printf("In %s: Split transfer function computation between %d threads for one mode/ic/type:\n",
		 __func__,number_of_threads);
#endif
	
	/* initialize error management flag */

	abort = _FALSE_;

	/* beginning of parallel region */

#pragma omp parallel							\
  shared (pw,ptr,ppr,ppt,index_mode,index_ic,index_tt,			\
	  interpolated_sources,abort,num_j)				\
  private (thread,index_l,tstart,tstop,tspent,address_in_workspace,tau0_minus_tau,delta_tau,sources,j_l,ddj_l,x_size_l,x_min_l)
	
	{
	  
#ifdef _OPENMP
	  thread = omp_get_thread_num();
	  tspent = 0.;
#endif

	  /* define address of each field in the workspace */
	  address_in_workspace = pw[thread];
	  
	  tau0_minus_tau = address_in_workspace;
	  address_in_workspace += ppt->tau_size;
	  
	  delta_tau  = address_in_workspace;
	  address_in_workspace += ppt->tau_size;

	  sources = address_in_workspace;
	  address_in_workspace += ptr->k_size[index_mode]*ppt->tau_size;
	  
	  x_min_l = address_in_workspace;
	  address_in_workspace += 1;

	  j_l = address_in_workspace;

	  /* the address of the next fields, ddj_l, etc., will be defined
	     within the l loop, since they depend on l */

	  /* copy the interpolated sources in the workspace */
	  
	  memcpy(sources,
		 interpolated_sources,
		 ptr->k_size[index_mode]*ppt->tau_size*sizeof(double));
	  
#pragma omp for schedule (dynamic)

	  for (index_l = 0; index_l < ptr->l_size[index_mode]; index_l++) {

#ifdef _OPENMP
	    tstart = omp_get_wtime();
#endif

	    /* for this value of l, retrieve the number of x values
	       from the bessel structure */
	    x_size_l=pbs->x_size[index_l];

	    /* copy from the bessel structure all bessel function
	       related information for this particular value of l:
	       x_min, the j_l(x)s, the j_l''(x)s. They are all
	       allocated in a contiguous memory zone, so we can use a
	       single memcpy. */
	    memcpy(x_min_l,pbs->x_min[index_l],(1+num_j*x_size_l)*sizeof(double));
	
	    /* now define the address of the ddj_l field (and
	       eventually additional fields in the workspace) */
	    ddj_l = j_l + x_size_l;

	    if (ppr->transfer_cut == tc_env) {
	      dj_l = ddj_l + x_size_l;
	      dddj_l = dj_l + x_size_l;
	    }

	    /* check that the computation will never need values of
	       j_l(x) with x > x_max (should never happen, since x_max
	       is chosen to be greater than tau0*k_max in bessel
	       module) */

	    class_test_parallel((int)((tau0_minus_tau[0] * ptr->k[index_mode][ptr->k_size[index_mode]-1] - (*x_min_l))/pbs->x_step)+1 >= x_size_l,
				ptr->error_message,
				"Increase x_max in bessel functions! The computation needs index_x up to %d while x_size[%d]=%d for x=%e\n",
				(int)((tau0_minus_tau[0] * ptr->k[index_mode][ptr->k_size[index_mode]-1] - (*x_min_l))/pbs->x_step)+1,
				index_l,
				x_size_l,
				tau0_minus_tau[0] * ptr->k[index_mode][ptr->k_size[index_mode]-1]);

	    /* compute the transfer function for this l */
	    class_call_parallel(transfer_compute_for_each_l(ppr,
							    ppt,
							    ptr,
							    index_mode,
							    index_ic,
							    index_tt,
							    index_l,
							    (double)ptr->l[index_l],
							    *x_min_l,
							    pbs->x_step,
							    tau0_minus_tau,
							    delta_tau,
							    sources,
							    j_l,
							    ddj_l,
							    dj_l,
							    dddj_l),
				ptr->error_message,
				ptr->error_message);

#ifdef _OPENMP
	    tstop = omp_get_wtime();

	    tspent += tstop-tstart;
#endif

#pragma omp flush(abort)

	  } /* end of loop over l */

#ifdef _OPENMP
	  if (ptr->transfer_verbose>1)
	    printf("In %s: time spent in parallel region (loop over l's) = %e s for thread %d\n",
		   __func__,tspent,omp_get_thread_num());
#endif

	} /* end of parallel region */
	
	if (abort == _TRUE_) return _FAILURE_;
	
      } /* end of loop over type */
      
    } /* end of loop over initial condition */

    free(interpolated_sources);
    free(source_spline);    

  } /* end of loop over mode */

#pragma omp parallel shared(pw) private(thread)
  {
#ifdef _OPENMP
    thread = omp_get_thread_num();
#endif
    free(pw[thread]);
  }

  free(pw);

  return _SUCCESS_;
}

/**
 * This routine frees all the memory space allocated by transfer_init().
 *
 * To be called at the end of each run, only when no further calls to
 * transfer_functions_at_k() are needed.
 *
 * @param ptr Input: pointer to transfers structure (which fields must be freed)
 * @return the error status
 */

int transfer_free(
		  struct transfers * ptr
		  ) {

  int index_mode;

  if (ptr->has_cls == _TRUE_) {

    for (index_mode = 0; index_mode < ptr->md_size; index_mode++) {
      free(ptr->k[index_mode]);
      free(ptr->transfer[index_mode]);
    }  
   
    free(ptr->tt_size);
    free(ptr->l_size);
    free(ptr->l);
    free(ptr->k_size);
    free(ptr->k);
    free(ptr->transfer);
    
  }

  return _SUCCESS_;
  
}

/**
 * This routine defines all indices and allocates all tables 
 * in the transfers structure 
 *
 * Compute list of (k, l) values, allocate and fill corresponding
 * arrays in the transfers structure. Allocate the array of transfer
 * function tables.
 *
 * @param ppr Input : pointer to precision structure 
 * @param ppt Input : pointer to perturbation structure
 * @param pbs Input : pointer to bessels structure
 * @param ptr Input/Output: pointer to transfer structure
 * @param rs_rec  Input : comoving distance to recombination
 * @return the error status
 */

int transfer_indices_of_transfers(
				  struct precision * ppr,
				  struct perturbs * ppt,
				  struct bessels * pbs,
				  struct transfers * ptr,
				  double rs_rec
				  ) {

  /** Summary: */

  /** - define local variables */

  int index_mode,index_tt,index_tt_common;

  /** define indices for transfer types */

  class_alloc(ptr->tt_size,ptr->md_size * sizeof(int),ptr->error_message);

  /** - type indices common to scalars and tensors */

  index_tt = 0;

  if (ppt->has_cl_cmb_temperature == _TRUE_) {
    ptr->index_tt_t = index_tt;
    index_tt++;
  }

  if (ppt->has_cl_cmb_polarization == _TRUE_) {
    ptr->index_tt_e = index_tt;
    index_tt++;
  }

  index_tt_common=index_tt;

  /** - type indices for scalars */

  if (ppt->has_scalars == _TRUE_) {
  
    index_tt = index_tt_common;

    if (ppt->has_cl_cmb_lensing_potential == _TRUE_) {
      ptr->index_tt_lcmb = index_tt;
      index_tt++;
    }

    ptr->tt_size[ppt->index_md_scalars]=index_tt;

  }

  /** - type indices for tensors */

  if (ppt->has_tensors == _TRUE_) {
  
    index_tt = index_tt_common;

    if (ppt->has_cl_cmb_polarization == _TRUE_) {
      ptr->index_tt_b = index_tt;
      index_tt++;
    }

    ptr->tt_size[ppt->index_md_tensors]=index_tt;

  }

  /** - allocate arrays of (k, l) values and transfer functions */

  /* number of l values for each mode, l_size[index_mode] */

  class_alloc(ptr->l_size,ptr->md_size * sizeof(int),ptr->error_message);

  /* number of k values for each mode, k_size[index_mode] */

  class_alloc(ptr->k_size,ptr->md_size * sizeof(int),ptr->error_message);

  /* list of k values for each mode, k[index_mode] */

  class_alloc(ptr->k,ptr->md_size * sizeof(double *),ptr->error_message);

  /* array (of array) of transfer functions for each mode, transfer[index_mode] */

  class_alloc(ptr->transfer,ptr->md_size * sizeof(double *),ptr->error_message);

  /** get l values using transfer_get_l_list() */
  class_call(transfer_get_l_list(ppr,ppt,pbs,ptr),
	     ptr->error_message,
	     ptr->error_message);
  
  /** - loop over modes (scalar, etc). For each mode: */
  
  for (index_mode = 0; index_mode < ptr->md_size; index_mode++) {

    /** (a) get k values using transfer_get_k_list() */
    class_call(transfer_get_k_list(ppr,ppt,ptr,rs_rec,index_mode),
	       ptr->error_message,
	       ptr->error_message);

    /** (b) allocate arrays of transfer functions, (ptr->transfer[index_mode])[index_ic][index_tt][index_l][index_k] */
    class_alloc(ptr->transfer[index_mode],
		ppt->ic_size[index_mode] * ptr->tt_size[index_mode] * ptr->l_size[index_mode] * ptr->k_size[index_mode] * sizeof(double),
		ptr->error_message);
    
  }

  return _SUCCESS_;

}

/**
 * This routine defines the number and values of mutipoles l for all modes.
 *
 * @param ppr  Input : pointer to precision structure
 * @param ppt  Input : pointer to perturbation structure
 * @param pbs  Input : pointer to bessels structure
 * @param ptr  Input/Output : pointer to transfers structure containing l's
 * @return the error status
 */

int transfer_get_l_list(
			struct precision * ppr,
			struct perturbs * ppt,
			struct bessels * pbs,
			struct transfers * ptr
			) {

  int index_l;
  int l_max=0;
  int index_mode;

  ptr->l_size_max=0;

  for (index_mode=0; index_mode < ppt->md_size; index_mode++) {

    if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {
      l_max = ppt->l_scalar_max;
    }
    
    if ((ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors)) {
      l_max = ppt->l_tensor_max;
    }

    class_test(l_max > pbs->l[pbs->l_size-1],
	       ptr->error_message,
	       "For mode %d, asked for l_max=%d greater than in Bessel table where l_max=%d",
	       index_mode,
	       l_max,
	       pbs->l[pbs->l_size-1]);

    index_l=0;

    while((index_l < pbs->l_size-1) && (pbs->l[index_l] <= l_max)) {
      index_l++;
    }
  
    ptr->l_size[index_mode] = index_l+1;

    ptr->l_size_max = max(ptr->l_size_max,ptr->l_size[index_mode]);

  }
     
  class_alloc(ptr->l,ptr->l_size_max*sizeof(int),ptr->error_message);
  
  for (index_l=0; index_l < ptr->l_size_max; index_l++) {
    ptr->l[index_l]=pbs->l[index_l];
  }
  
  return _SUCCESS_;

}

/**
 * This routine defines the number and values of wavenumbers k for
 * each mode (different in perturbation module and transfer module:
 * here we impose an upper bound on the linear step. So, typically,
 * for small k, the sampling is identical to that in the perturbation
 * module, while at high k it is denser and source functions are
 * interpolated).
 *
 * @param ppr     Input : pointer to precision structure
 * @param ppt     Input : pointer to perturbation structure
 * @param ptr     Input/Output : pointer to transfers structure containing k's
 * @param rs_rec  Input : comoving distance to recombination
 * @param index_mode Input: index of requested mode (scalar, tensor, etc) 
 * @return the error status
 */

int transfer_get_k_list(
			struct precision * ppr,
			struct perturbs * ppt,
			struct transfers * ptr,
			double rs_rec,
			int index_mode
			) {

  int index_k_pt;
  int index_k_tr;
  double k,k_min,k_max,k_step_max=0.;

  /* find k_step_max, the maximum value of the step */

  if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

    k_step_max = 2.*_PI_/rs_rec*ppr->k_step_trans_scalars; /* step_size, inferred from precision_params structure */

  }

  if ((ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors)) {

    k_step_max = 2.*_PI_/rs_rec*ppr->k_step_trans_tensors; /* step_size, inferred from precision_params structure */

  }

  class_test(k_step_max == 0.,
	     ptr->error_message,
	     "stop to avoid infinite loop");

  /* first and last value in perturbation module */

  k_min = ppt->k[index_mode][0]; /* first value, inferred from perturbations structure */

  k_max = ppt->k[index_mode][ppt->k_size_cl[index_mode]-1]; /* last value, inferred from perturbations structure */

  /* first, count the number of necessary values */

  index_k_pt = 0;
  index_k_tr = 0;

  /* - first point */

  k = k_min;
  index_k_pt++;
  index_k_tr++;

  /* - points taken from perturbation module if step smkall enough */

  while ((index_k_pt < ppt->k_size[index_mode]) && ((ppt->k[index_mode][index_k_pt] -k) < k_step_max)) {
      k = ppt->k[index_mode][index_k_pt];
      index_k_pt++;
      index_k_tr++;
  }

  /* - then, points spaced linearily with step k_step_max */

  while (k < k_max) {
      k += k_step_max;
      index_k_tr++;
  }

  /* - get number of points and allocate list */

  if (k > k_max)
    ptr->k_size[index_mode]=index_k_tr-1;
  else
    ptr->k_size[index_mode]=index_k_tr;

  class_alloc(ptr->k[index_mode],
	      ptr->k_size[index_mode]*sizeof(double),
	      ptr->error_message);

  /* repeat exactly the same steps, but now filling the list */

  index_k_pt = 0;
  index_k_tr = 0;

  ptr->k[index_mode][0] = k_min;
  index_k_pt++;
  index_k_tr++;

  while ((index_k_pt < ppt->k_size[index_mode]) && ((ppt->k[index_mode][index_k_pt] -k) < k_step_max)) {
      k = ppt->k[index_mode][index_k_pt];
      ptr->k[index_mode][index_k_tr] = k;
      index_k_pt++;
      index_k_tr++;
  }

  while ((index_k_tr < ptr->k_size[index_mode]) && (k < k_max)) {
      k += k_step_max;
      ptr->k[index_mode][index_k_tr] = k;
      index_k_tr++;
  }

  /* consistency check */

  class_test(ptr->k[index_mode][ptr->k_size[index_mode]-1] > k_max,
	     ptr->error_message,
	     "bug in k list calculation, k_max larger in transfer than in perturb, should never happen");

  return _SUCCESS_; 

}

/**
 * This routine interpolates sources \f$ S(k, \tau) \f$ for each mode, 
 * initial condition and type, to get them at the right values of k,
 * using the spline interpolation method. 
 *
 * Important: some physics enters here. This is indeed the most
 * efficient place for mutiplying the sources by some factors
 * in order to transform the 'raw source functions' with a given
 * 'source types' into an 'observable source funtion' with a given
 * 'transfer type'.
 *
 * E.g.: here we can multiply the gravitational potential source 
 * by one or several window functions, transforming it into one or 
 * several lensing source functions.
 * 
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param tau0                  Input : conformal time today
 * @param tau_rec               Input : conformal time at recombination
 * @param index_mode            Input : index of mode
 * @param index_ic              Input : index of initial condition
 * @param index_tt              Input : index of type of transfer
 * @param source_spline         Output: array of second derivative of sources (filled here but allocated in transfer_init() to avoid numerous reallocation)
 * @param interpolated_sources  Output: array of interpolated sources (filled here but allocated in transfer_init() to avoid numerous reallocation)
 * @return the error status
 */

int transfer_interpolate_sources(
				 struct perturbs * ppt,
				 struct transfers * ptr,
				 double tau0,
				 double tau_rec,
				 int index_mode,
				 int index_ic,
				 int index_tt,
				 double * source_spline, /* array with argument source_spline[index_tau*ppt->k_size[index_mode]+index_k] (must be allocated) */
				 double * interpolated_sources /* array with argument interpolated_sources[index_k_tr*ppt->tau_size+index_tau] (must be allocated) */
				 ) {

  /** Summary: */

  /** - define local variables */

  /* index running on k values in the original source array */
  int index_k;

  /* index running on time */
  int index_tau;

  /* index running on type of source (not type of transfer) */
  int index_type=0;

  /* index running on k values in the interpolated source array */
  int index_k_tr;

  /* variables used for spline interpolation algorithm */
  double h, a, b;

  /** - which source are we considering? 
        Define correspondence between transfer types and source types */

  if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

    if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t)) 
      index_type=ppt->index_tp_t;

    if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_e)) 
      index_type=ppt->index_tp_e;

    if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (index_tt == ptr->index_tt_lcmb)) 
      index_type=ppt->index_tp_g;
    
  }

  if ((ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors)) {

    if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t)) 
      index_type=ppt->index_tp_t;

    if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_e)) 
      index_type=ppt->index_tp_e;

    if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_b))
      index_type=ppt->index_tp_b;

  }
	  
  /** - find second derivative of original sources with respect to k
        in view of spline interpolation */

  class_call(array_spline_table_columns(ppt->k[index_mode],
					ppt->k_size[index_mode],
					ppt->sources[index_mode][index_ic * ppt->tp_size[index_mode] + index_type],
					ppt->tau_size,
					source_spline,
					_SPLINE_EST_DERIV_,
					ptr->error_message),
	     ptr->error_message,
	     ptr->error_message);

  /** - interpolate at each k value using the usual 
        spline interpolation algorithm. 
        Eventually mutiply the result by a factor accounting for the 
        difference between 'raw source functions' in the perturbation module
        and 'observable source functions' in the transfer module
        (e.g. gravitational potential source -> lensing source) */

  index_k = 0;
  h = ppt->k[index_mode][index_k+1] - ppt->k[index_mode][index_k];

  for (index_k_tr = 0; index_k_tr < ptr->k_size[index_mode]; index_k_tr++) {

    while (((index_k+1) < ppt->k_size[index_mode]) &&
	   (ppt->k[index_mode][index_k+1] < 
	    ptr->k[index_mode][index_k_tr])) {
      index_k++;
      h = ppt->k[index_mode][index_k+1] - ppt->k[index_mode][index_k];
    }

    class_test(h==0.,
	       ptr->error_message,
	       "stop to avoid division by zero");

    b = (ptr->k[index_mode][index_k_tr] - ppt->k[index_mode][index_k])/h;
    a = 1.-b;
    
    for (index_tau = 0; index_tau < ppt->tau_size; index_tau++) {

      /**   a) interpolate for each value of conformal time */

      interpolated_sources[index_k_tr*ppt->tau_size+index_tau] = 
	a * ppt->sources[index_mode]
	[index_ic * ppt->tp_size[index_mode] + index_type]
	[index_tau*ppt->k_size[index_mode]+index_k]
	+ b * ppt->sources[index_mode]
	[index_ic * ppt->tp_size[index_mode] + index_type]
	[index_tau*ppt->k_size[index_mode]+index_k+1]
	+ ((a*a*a-a) * source_spline[index_tau*ppt->k_size[index_mode]+index_k]
	   +(b*b*b-b) * source_spline[index_tau*ppt->k_size[index_mode]+index_k+1])*h*h/6.0;

      /**   b) case of cmb lensing: multiply gravitational potential 
               by appropriate window function */

      if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

	if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (index_tt == ptr->index_tt_lcmb)) {

	  /* lensing source =  W(tau) psi(k,tau) H(tau-tau_rec) 
	     with 
	     psi = (newtonian) gravitationnal potential  
	     W = 2(tau-tau_rec)/(tau_0-tau)/(tau_0-tau_rec) 
	     H(x) = Heaviside
	     (in tau = tau_0, set source = 0 to avoid division by zero;
              regulated anyway by Bessel).
	  */
	  if ((ppt->tau_sampling[index_tau] > tau_rec) && 
	      ((tau0-ppt->tau_sampling[index_tau]) > 0.)) {
	    interpolated_sources[index_k_tr*ppt->tau_size+index_tau] *=
	      2.*(ppt->tau_sampling[index_tau]-tau_rec)
	      /(tau0-ppt->tau_sampling[index_tau])
	      /(tau0-tau_rec);
	  }
	  else {
	    interpolated_sources[index_k_tr*ppt->tau_size+index_tau] = 0;
	  }
	}
      }
    }
  }

  return _SUCCESS_;

}

/**
 * This routine computes the transfer functions \f$ \Delta_l^{X} (k) \f$)
 * as a function of wavenumber k for each mode, initial condition,
 * type and multipole l passed in input. 
 *
 * For a given value of k, the transfer function is infered from 
 * the source function (passed in input in the array interpolated_sources)
 * and from Bessel functions (passed in input in the bessels structure),
 * either by convolving them along tau, or by a Limber appoximation.
 * This elementary task is distributed either to transfer_integrate()
 * or to transfer_limber(). The task of this routine is mainly to
 * loop over k values, and to decide at which k_max the calculation can
 * be stopped, according to some approximation scheme designed to find a 
 * compromise between execution time and precision. The approximation scheme
 * is defined by parameters in bthe precision structure.
 * 
 * @param ppr                   Input : pointer to precision structure 
 * @param ppt                   Input : pointer to perturbation structure
 * @param pbs                   Input : pointer to bessels structure 
 * @param ptr                   Input/output : pointer to transfers structure (result stored there)
 * @param tau0                  Input : conformal time today
 * @param tau_rec               Input : conformal time at recombination
 * @param index_mode            Input : index of mode
 * @param index_ic              Input : index of initial condition
 * @param index_tt              Input : index of type of transfer
 * @param index_l               Input : index of multipole
 * @param interpolated_sources  Input : array containing the sources
 * @param ptw                   Input : pointer to transfer_workspace structure (allocated in transfer_init() to avoid numerous reallocation) 
 * @return the error status
 */

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
				double * tau0_minus_tau,
				double * delta_tau,
				double * sources,
				double * j_l,
				double * ddj_l,
				double * dj_l,
				double * dddj_l
				){

  /** Summary: */

  /** - define local variables */

  /* running index for wavenumbers */	
  int index_k;
  /* current wavenumber value */
  double k;
  /* flag: for a given l, should the transfer functions stop being computed at next value of k? */
  short cut_transfer;
  /* global maximum of \f$ \Delta_l(k) \f$ as a function of k, used for stopping computation */
  double global_max;
  /* global minimum of \f$ \Delta_l(k) \f$ as a function of k, used for stopping computation */
  double global_min;
  /* last local maximum of \f$ \Delta_l(k) \f$ as a function of k, used for stopping computation */
  double last_local_max;
  /* last local minimum of \f$ S(k, \tau) j_l(k [\tau_0 - \tau]) \f$ as a function of k, used for cutting the integral */
  double last_local_min;
  /* value of transfer function */
  double transfer_function;
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

  /* quantitites for multiplying some transfer function by a factor */
  short multiply_by_factor;
  double extra_factor=1.;

  if (ptr->transfer_verbose > 2)
    printf("Compute transfer for l=%d\n",(int)l);

  /** - if the option of stopping the transfer function computation at some k_max is selected, initialize relevant quantities */
	      
  if (ppr->transfer_cut == tc_osc) {
    cut_transfer = _FALSE_;
    global_max=0.;
    global_min=0.;
    last_local_max=0.;
    last_local_min=0.;
    transfer_function=0.;
    transfer_last=0.;
    transfer_last_last=0;
  }
	      
  if (ppr->transfer_cut == tc_cl) {
    cut_transfer = _FALSE_;
    cl=0.;
    cl_var=1.;
    cl_var_last=1.;
    cl_var_last_last=1.;
  }

  /** - is there an extra factor in the transfer function,
      besides the integral over S*j_l(x)? */

  multiply_by_factor=_FALSE_;

  /* for scalar (E-)polarization, multiply by 
     square root of  (l+2)(l+1)l(l-1) */
  
  if ((((ppt->has_scalars == _TRUE_) && 
	(index_mode == ppt->index_md_scalars)) && 
       ((ppt->has_cl_cmb_polarization == _TRUE_) && 
	(index_tt == ptr->index_tt_e)))) {
    
    multiply_by_factor=_TRUE_;
    extra_factor=sqrt((l+2.) * (l+1.) * l * (l-1.));
  }
  
  /* for tensor temperature, multiply by 
     square root of (l+2)(l+1)l(l-1)/2 */
  
  if ((((ppt->has_tensors == _TRUE_) && 
	(index_mode == ppt->index_md_tensors)) &&
       ((ppt->has_cl_cmb_temperature == _TRUE_) && 
	(index_tt == ptr->index_tt_t)))) {
    
    multiply_by_factor=_TRUE_;
    extra_factor=sqrt((l+2.) * (l+1.) * l * (l-1.));
  }
	      
  /** - loop over k. For each k:
        (a) if the option of stopping the transfer function computation 
            at some k_max is not selected or if we are below k_max, 
            compute \f$ \Delta_l(k) \f$ either with transfer_integrate() 
            or with transfer_limber(); store the result in the table
            in transfers structure;
        (b) if it is selected and we are above k_max, set the transfer function 
            to zero; 
        (c) if it is selected and we are below k_max, check wether k_max is 
            being reached. 
  */

  for (index_k = 0; index_k < ptr->k_size[index_mode]; index_k++) {

    k = ptr->k[index_mode][index_k];

    if (ptr->transfer_verbose > 3)
      printf("Compute transfer for l=%d k=%e type=%d\n",(int)l,k,index_tt);
		
    /* update previous transfer values in the tc_osc method */
    if (ppr->transfer_cut == tc_osc) {
      transfer_last_last = transfer_last;
      transfer_last = transfer_function;
    }
		
    /* update previous relative pseudo-C_l variation in the tc_cl method */
    if (ppr->transfer_cut == tc_cl) {
      cl_var_last_last = cl_var_last;
      cl_var_last = cl_var;
    }
		
    /* below k_max: compute transfer function */
    if ((ppr->transfer_cut == tc_none) || (cut_transfer == _FALSE_)) {

      /* criterium for chosing between integration and Limber 
	 must be implemented here */

      if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (index_tt == ptr->index_tt_lcmb) && (l>ppr->l_switch_limber)) {
	
      	class_call(transfer_limber(ppt->tau_size,
				   ptr,
				   index_mode,
				   index_k,
				   l,
				   k,
				   tau0_minus_tau,
				   sources,
				   &transfer_function),
		   ptr->error_message,
		   ptr->error_message); 

      }
      else {
	
	if ((ppr->transfer_cut == tc_env) &&
	    (((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t)) || 
	     ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_t)))) {
	  
	  class_call(transfer_envelop(ppt->tau_size,
				      index_k,
				      l,
				      k,
				      x_min_l,
				      x_step,
				      tau0_minus_tau,
				      delta_tau,
				      sources,
				      j_l,
				      ddj_l,
				      dj_l,
				      dddj_l,
				      &transfer_function),
		     ptr->error_message,
		     ptr->error_message);
	}
	else {

	  class_call(transfer_integrate(ppt->tau_size,
					index_k,
					l,
					k,
					x_min_l,
					x_step,
					tau0_minus_tau,
					delta_tau,
					sources,
					j_l,
					ddj_l,
					&transfer_function),
		     ptr->error_message,
		     ptr->error_message);
	}
      }
    }

    /* above k_max: set transfer function to zero */
    else {
      transfer_function = 0.;
    }

    /* eventually multiply by extra factor */
    if (multiply_by_factor == _TRUE_) transfer_function *= extra_factor;

    /* store transfer function in transfer structure */
    ptr->transfer[index_mode][((index_ic * ptr->tt_size[index_mode] + index_tt)
			       * ptr->l_size[index_mode] + index_l)
			      * ptr->k_size[index_mode] + index_k]
      = transfer_function;
    
    /* in the tc_osc case, update various quantities and 
       check wether k_max is reached */
    if ((ppr->transfer_cut == tc_osc) && (index_k>=2)) {
      
      /* detect local/global maximum of \f$ \Delta_l(k) \f$; 
	 if detected, check the criteria for reaching k_max */
      if ((transfer_last > 0.) && (transfer_function < transfer_last) && (transfer_last > transfer_last_last)) {
	last_local_max = transfer_last;
	if (last_local_max > global_max) {
	  global_max = last_local_max;
	}
	if ((last_local_max-last_local_min) < ppr->transfer_cut_threshold_osc * (global_max-global_min)) {
	  cut_transfer = _TRUE_;
	}
	
      }
		  
      /* detect local/global minimum of \f$ \Delta_l(k) \f$; 
	 if detected, check the criteria for reaching k_max */ 
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
		
    /* in the tc_cl case, update various quantities and check wether k_max is reached */
    if ((ppr->transfer_cut == tc_cl) && (index_k>=2) && (index_k<ptr->k_size[index_mode]-1) && (transfer_function != 0.)) {
		  
      /* rough estimate of the contribution of the last step to pseudo-C_l, 
	 assuming flat primordial spectrum */
      delta_cl = transfer_function * transfer_function / k * 0.5 * (ptr->k[index_mode][index_k+1] - ptr->k[index_mode][index_k-1]);
      
      /* update pseudo-C_l */
      cl += delta_cl;
      
      class_test(cl == 0.,
		 ptr->error_message,
		 "stop to avoid division by zero");

      /* compute its relative variation */
      cl_var = delta_cl / cl;
      
      /* check if k_max is reached */
      if ((cl_var < cl_var_last) && (cl_var_last > cl_var_last_last)) {
	if (cl_var_last < ppr->transfer_cut_threshold_cl) {
	  cut_transfer = _TRUE_;
	}

      }
    }

  } /* end of loop over k */

  return _SUCCESS_;

}

/**
 * This routine computes the transfer functions \f$ \Delta_l^{X} (k) \f$)
 * for each mode, initial condition, type, multipole l and wavenumber k,
 * by convolving  the source function (passed in input in the array 
 * interpolated_sources) with Bessel functions (passed in input in the 
 * bessels structure).
 *
 * @param ppt                   Input : pointer to perturbation structure
 * @param pbs                   Input : pointer to bessels structure 
 * @param ptr                   Input : pointer to transfers structure
 * @param tau0                  Input : conformal time today
 * @param tau_rec               Input : conformal time at recombination
 * @param index_mode            Input : index of mode
 * @param index_tt              Input : index of type
 * @param index_l               Input : index of multipole
 * @param index_k               Input : index of wavenumber
 * @param interpolated_sources  Input: array of interpolated sources
 * @param ptw                   Input : pointer to transfer_workspace structure (allocated in transfer_init() to avoid numerous reallocation) 
 * @param trsf                  Output: transfer function \f$ \Delta_l(k) \f$ 
 * @return the error status
 */

int transfer_integrate(
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
		       ) {

  /** Summary: */

  /** - define local variables */

  /* minimum value of \f$ (\tau0-\tau) \f$ at which \f$ j_l(k[\tau_0-\tau]) \f$ is known, given that \f$ j_l(x) \f$ is sampled above some finite value \f$ x_{\min} \f$ (below which it can be approximated by zero) */  
  double tau0_minus_tau_min_bessel;

  /* running value of bessel function */
/*   double bessel; */
  int index_x;
  double x;
  double a;
  
  double transfer;

  /* index in the source's tau list corresponding to the last point in the overlapping region between sources and bessels */
  int index_tau,index_tau_max;

  /** - find minimum value of (tau0-tau) at which \f$ j_l(k[\tau_0-\tau]) \f$ is known, given that \f$ j_l(x) \f$ is sampled above some finite value \f$ x_{\min} \f$ (below which it can be approximated by zero) */  
  tau0_minus_tau_min_bessel = x_min_l/k; /* segmentation fault impossible, checked before that k != 0 */

  /** - if there is no overlap between the region in which bessels and sources are non-zero, return zero */
  if (tau0_minus_tau_min_bessel >= tau0_minus_tau[0]) {
    *trsf = 0.;
    return _SUCCESS_;
  }

  /** - if there is an overlap: */

  /** (a) find index in the source's tau list corresponding to the last point in the overlapping region. After this step, index_tau_max can be as small as zero, but not negative. */ 
  index_tau_max = tau_size-1;
  while (tau0_minus_tau[index_tau_max] < tau0_minus_tau_min_bessel)
    index_tau_max--;

  /** (b) the source function can vanish at large $\f \tau \f$. Check if further points can be eliminated. After this step and if we did not return a null transfer function, index_tau_max can be as small as zero, but not negative. */
  while (sources[index_k * tau_size + index_tau_max] == 0.) { 
    index_tau_max--;
    if (index_tau_max < 0) {
      *trsf = 0.;
      return _SUCCESS_;
    }
  }

  /** (c) integrate with trapezoidal method */
    
  /* for bessel function interpolation, we could call the subroutine bessel_at_x; however we perform operations directly here in order to speed up the code */
  
  x = k * tau0_minus_tau[index_tau_max];

  index_x = (int)((x-x_min_l)/x_step);
  
  a = (x_min_l+x_step*(index_x+1) - x)/x_step;
  
  transfer = sources[index_k * tau_size + index_tau_max] /* source */
    * (a * j_l[index_x] +                                /* bessel (cubic spline interpolation) */
       (1.-a) * (j_l[index_x+1]
		 - a * ((a+1.) * ddj_l[index_x]
			+(2.-a) * ddj_l[index_x+1]) 
		 * x_step * x_step / 6.0) );
  
  if (index_tau_max > 0)
    transfer *= (tau0_minus_tau[index_tau_max-1]-tau0_minus_tau_min_bessel);
  else
    transfer *= (tau0_minus_tau[index_tau_max]-tau0_minus_tau_min_bessel);

  for (index_tau=0; index_tau<index_tau_max; index_tau++) {
    
    /* for bessel function interpolation, we could call the subroutine bessel_at_x; however we perform operations directly here in order to speed up the code */
    
    x = k * tau0_minus_tau[index_tau];
    
    index_x = (int)((x-x_min_l)/x_step);
    
    a = (x_min_l+x_step*(index_x+1) - x)/x_step;
    
    transfer += sources[index_k * tau_size + index_tau] /* source */
      * (a * j_l[index_x]                               /* bessel (cubic spline interpolation) */
	 + (1.-a) * ( j_l[index_x+1]
		      - a * ((a+1.) * ddj_l[index_x]
			     +(2.-a) * ddj_l[index_x+1]) 
		      * x_step * x_step / 6.0)) 
      * delta_tau[index_tau];                           /* dtau */
  }
  
  *trsf = 0.5*transfer; /* correct for factor 1/2 from trapezoidal rule */

  return _SUCCESS_;
}

/**
 * This routine computes the transfer functions \f$ \Delta_l^{X} (k) \f$)
 * for each mode, initial condition, type, multipole l and wavenumber k,
 * by using the Limber approximation, i.e by evaluating the source function 
 * (passed in input in the array interpolated_sources) at a single value of
 * tau (the Bessel function being approximated as a Dirac distribution)
 *
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param tau0                  Input : conformal time today
 * @param index_mode            Input : index of mode
 * @param index_tt              Input : index of type
 * @param index_l               Input : index of multipole
 * @param index_k               Input : index of wavenumber
 * @param interpolated_sources  Input: array of interpolated sources
 * @param trsf                  Output: transfer function \f$ \Delta_l(k) \f$ 
 * @return the error status
 */

int transfer_limber(
		    int tau_size,
		    struct transfers * ptr,
		    int index_mode,
		    int index_k,
		    double l,
		    double k,
		    double * tau0_minus_tau,
		    double * sources, /* array with argument interpolated_sources[index_k_tr*ppt->tau_size+index_tau] */
		    double * trsf
		    ){

  /** Summary: */

  /** - define local variables */

  /* conformal time at which source must be computed */
  double tau0_minus_tau_limber;
  double transfer;

  /** - get k, l and infer tau such that k(tau0-tau)=l+1/2; 
      check that tau is in appropriate range */

  tau0_minus_tau_limber = (l+0.5)/k;

  if (tau0_minus_tau_limber > tau0_minus_tau[0]) {
    *trsf = 0.;
    return _SUCCESS_;
  }

  /** - get source at this value tau */

  class_call(array_interpolate_two_arrays_one_column(
						     tau0_minus_tau,
						     sources,
						     ptr->k_size[index_mode],
						     index_k,
						     tau_size,
						     tau0_minus_tau_limber,
						     &transfer,
						     ptr->error_message),
	     ptr->error_message,
	     ptr->error_message);

  /** - get transfer = source * sqrt(pi/(2l+1))/k */

  transfer *= sqrt(_PI_/(2.*l+1.))/k;
  
  *trsf = transfer;

  return _SUCCESS_;
  
}

/**
 * This routine computes the envelop of transfer functions \f$
 * \Delta_l^{X} (k) \f$) for each mode, initial condition, type,
 * multipole l and wavenumber k, by convolving the source function
 * (passed in input in the array interpolated_sources) with Bessel
 * functions (passed in input in the bessels structure).
 *
 * @param ppt                   Input : pointer to perturbation structure
 * @param pbs                   Input : pointer to bessels structure 
 * @param ptr                   Input : pointer to transfers structure
 * @param tau0                  Input : conformal time today
 * @param tau_rec               Input : conformal time at recombination
 * @param index_mode            Input : index of mode
 * @param index_tt              Input : index of type
 * @param index_l               Input : index of multipole
 * @param index_k               Input : index of wavenumber
 * @param interpolated_sources  Input: array of interpolated sources
 * @param ptw                   Input : pointer to transfer_workspace structure (allocated in transfer_init() to avoid numerous reallocation) 
 * @param trsf                  Output: transfer function \f$ \Delta_l(k) \f$ 
 * @return the error status
 */

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
		     ) {

  /** Summary: */

  /** - define local variables */

  /* minimum value of \f$ (\tau0-\tau) \f$ at which \f$ j_l(k[\tau_0-\tau]) \f$ is known, given that \f$ j_l(x) \f$ is sampled above some finite value \f$ x_{\min} \f$ (below which it can be approximated by zero) */  
  double tau0_minus_tau_min_bessel;

  /* running value of bessel function */
/*   double bessel; */
  int index_x;
  double x;
  double a;
  
  double transfer, transfer_shifted, source_times_dtau;

  /* index in the source's tau list corresponding to the last point in the overlapping region between sources and bessels */
  int index_tau,index_tau_max;

  /** - find minimum value of (tau0-tau) at which \f$ j_l(k[\tau_0-\tau]) \f$ is known, given that \f$ j_l(x) \f$ is sampled above some finite value \f$ x_{\min} \f$ (below which it can be approximated by zero) */  
  tau0_minus_tau_min_bessel = x_min_l/k; /* segmentation fault impossible, checked before that k != 0 */

  /** - if there is no overlap between the region in which bessels and sources are non-zero, return zero */
  if (tau0_minus_tau_min_bessel >= tau0_minus_tau[0]) {
    *trsf = 0.;
    return _SUCCESS_;
  }

  /** - if there is an overlap: */

  /** (a) find index in the source's tau list corresponding to the last point in the overlapping region. After this step, index_tau_max can be as small as zero, but not negative. */ 
  index_tau_max = tau_size-1;
  while (tau0_minus_tau[index_tau_max] < tau0_minus_tau_min_bessel)
    index_tau_max--;

  /** (b) the source function can vanish at large $\f \tau \f$. Check if further points can be eliminated. After this step and if we did not return a null transfer function, index_tau_max can be as small as zero, but not negative. */
  while (sources[index_k * tau_size + index_tau_max] == 0.) { 
    index_tau_max--;
    if (index_tau_max < 0) {
      *trsf = 0.;
      return _SUCCESS_;
    }
  }

  /** (c) integrate with trapezoidal method */
    
  /* for bessel function interpolation, we could call the subroutine bessel_at_x; however we perform operations directly here in order to speed up the code */
  
  x = k * tau0_minus_tau[index_tau_max];

  index_x = (int)((x-x_min_l)/x_step);
  
  a = (x_min_l+x_step*(index_x+1) - x)/x_step;

  source_times_dtau = sources[index_k * tau_size + index_tau_max];
  
  if (index_tau_max > 0) {
    source_times_dtau *= (tau0_minus_tau[index_tau_max-1]-tau0_minus_tau_min_bessel);
   }
  else {
    source_times_dtau *= (tau0_minus_tau[index_tau_max]-tau0_minus_tau_min_bessel);
  }

  /* source times j_l (cubic spline interpolation) */
  transfer = source_times_dtau
    * (a * j_l[index_x] +
       (1.-a) * (j_l[index_x+1]
		 - a * ((a+1.) * ddj_l[index_x]
			+(2.-a) * ddj_l[index_x+1]) 
		 * x_step * x_step / 6.0) );
  
  /* source times j_l' (cubic spline interpolation) */
  transfer_shifted = source_times_dtau
    * (a * dj_l[index_x] +
       (1.-a) * (dj_l[index_x+1]
		 - a * ((a+1.) * dddj_l[index_x]
			+(2.-a) * dddj_l[index_x+1]) 
		 * x_step * x_step / 6.0) );
  
  for (index_tau=0; index_tau<index_tau_max; index_tau++) {
    
    /* for bessel function interpolation, we could call the subroutine bessel_at_x; however we perform operations directly here in order to speed up the code */
    
    x = k * tau0_minus_tau[index_tau];
    
    index_x = (int)((x-x_min_l)/x_step);
    
    a = (x_min_l+x_step*(index_x+1) - x)/x_step;
    
  source_times_dtau = sources[index_k * tau_size + index_tau] * delta_tau[index_tau];
    
  /* source times j_l (cubic spline interpolation) */
  transfer += source_times_dtau * 
    (a * j_l[index_x] +
     (1.-a) * (j_l[index_x+1]
	       - a * ((a+1.) * ddj_l[index_x]
		      +(2.-a) * ddj_l[index_x+1]) 
	       * x_step * x_step / 6.0));
  
  /* source times j_l' (cubic spline interpolation) */
  transfer_shifted += source_times_dtau
    * (a * dj_l[index_x] +
       (1.-a) * (dj_l[index_x+1]
		 - a * ((a+1.) * dddj_l[index_x]
			+(2.-a) * dddj_l[index_x+1]) 
		 * x_step * x_step / 6.0));
  }
  
  *trsf = 0.5*sqrt(transfer*transfer + 1.3*transfer_shifted*transfer_shifted); /* correct for factor 1/2 from trapezoidal rule */

  return _SUCCESS_;
}
