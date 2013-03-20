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
 * @param index_md Input: index of requested mode
 * @param index_ic   Input: index of requested initial condition
 * @param index_tt   Input: index of requested type
 * @param index_l    Input: index of requested multipole
 * @param k          Input: any wavenumber
 * @param transfer_function Output: transfer function
 * @return the error status
 */
int transfer_functions_at_k(
			    struct transfers * ptr,
			    int index_md,
			    int index_ic,
			    int index_tt,
			    int index_l,
			    double k,
			    double * transfer_function
			    ) {
  /** Summary: */

  /** - interpolate in pre-computed table using array_interpolate_two() */
  class_call(array_interpolate_two(
				   ptr->k[index_md],
				   1,
				   0,
				   ptr->transfer[index_md]
				   +((index_ic * ptr->tt_size[index_md] + index_tt) * ptr->l_size[index_md] + index_l)
				   * ptr->k_size[index_md],
				   1,
				   ptr->k_size[index_md],
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
 *    transfer[index_md][((index_ic * ptr->tt_size[index_md] + index_tt) * ptr->l_size[index_md] + index_l) * ptr->k_size[index_md] + index_k]
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
  int index_md; 
  /* running index for initial conditions */
  int index_ic; 
  /* running index for types */
  int index_tt; 
  /* running index for multipoles */
  int index_l; 

  /* conformal time today */
  double tau0;
  /* conformal time at recombination */
  double tau_rec;

  /* number of sampling times for each transfer sources */
  int tau_size_tt;
  /* maximum number of sampling times for transfer sources */
  int tau_size_max;

  /* array of source derivatives S''(k,tau) 
     (second derivative with respect to k, not tau!), 
     used to interpolate sources at the right values of k,
     source_spline[index_tau*ppt->k_size[index_md]+index_k] */
  double * source_spline;

  /* table of source functions interpolated at the right values of k, 
     interpolated_sources[index_k_tr*ppt->tau_size+index_tau] */
  double * interpolated_sources;

  /* we deal with workspaces, i.e. with contiguous memory zones (one
     per thread) containing various fields used by the integration
     routine */

  /* - pointer used to assign adresses to the various workspace fields */
  double * address_in_workspace;

  /* - first workspace field: list of tau0-tau values, tau0_minus_tau[index_tau] */
  double * tau0_minus_tau;

  /* - second workspace field: list of delta_tau values, delta_tau[index_tau] */
  double * delta_tau;

  /* - third workspace field, containing just a double: number of time values */
  double * tau_size;

  /* - fourth workspace field, identical to above interpolated sources:
     sources[index_k_tr*ppt->tau_size+index_tau] */
  double * sources;

  /* - fifth workspace field, containing just a double: value of x at
     which bessel functions become non-negligible for a given l (
     infered from bessel module) */
  double * x_min_l;

  /* - sixth workspace field containing the list of j_l(x) values */
  double * j_l;

  /* - seventh workspace field containing the list of j_l''(x) values */
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

  /** - for a given l, maximum value of k such that we can convolve
      the source with Bessel functions j_l(x) without reaching x_max */
  double k_max_bessel;

  /** - array with the correspondance between the index of sources in
      the perturbation module and in the transfer module */
  int * tp_of_tt;

  /* a value of index_type */
  int previous_type;

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
  class_call(transfer_indices_of_transfers(ppr,ppt,pbs,ptr,tau0),
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

  for (index_md = 0; index_md < ptr->md_size; index_md++) {

    /** (a) allocate temporary arrays relevant for this mode */

    /** (a.0.) find correspondence between sources in the perturbation and transfer modules */

    class_alloc(tp_of_tt,ptr->tt_size[index_md]*sizeof(int),ptr->error_message);

    class_call(transfer_get_source_correspondence(ppt,ptr,index_md,tp_of_tt),
	       ptr->error_message,
	       ptr->error_message);

    /** (a.1.) arrays that will contain the sources interpolated at
	correct values of k, and their second derivatives with respect
	to k (for spline interpolation) */

    class_alloc(interpolated_sources,
		ptr->k_size[index_md]*ppt->tau_size*sizeof(double),
		ptr->error_message);
     
    class_alloc(source_spline,
		ppt->k_size[index_md]*ppt->tau_size*sizeof(double),
		ptr->error_message);

    /** (a.2.) maximum number of sampled times in the transfer
	sources: needs to be known here, in order to allocate a large
	enough workspace */

    tau_size_max = 0;
    
    for (index_tt = 0; index_tt < ptr->tt_size[index_md]; index_tt++) {

      class_call(transfer_source_tau_size(ppr,
					  pba,
					  ppt,
					  ptr,
					  tau_rec,
					  tau0,
					  index_md,
					  index_tt,
					  &tau_size_tt),
		 ptr->error_message,
		 ptr->error_message);
      
      tau_size_max = max(tau_size_max,tau_size_tt);
    }

    /** (a.3.) workspace, allocated in a parallel zone since in openmp
	version there is one workspace per thread */
    
    /* initialize error management flag */
    abort = _FALSE_;
    
    /* beginning of parallel region */
    
#pragma omp parallel				\
  shared(ptr,index_md,ppt,pbs,pw)		\
  private(thread)
    {
      
#ifdef _OPENMP
      thread = omp_get_thread_num();
#endif


      class_alloc_parallel(pw[thread],
			   ((2+ptr->k_size[index_md])*tau_size_max
			    +2+num_j*pbs->x_size_max)*
			   sizeof(double),
			   ptr->error_message);
            
    } /* end of parallel region */
    
    if (abort == _TRUE_) return _FAILURE_;
    
    /** (b) now loop over initial conditions and types: For each of them: */

    for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++) {
      
      /* initialize the previous type index */
      previous_type=-1;

      for (index_tt = 0; index_tt < ptr->tt_size[index_md]; index_tt++) {

        /* (b.1) check if we must now deal with a new source with a
	   new index ppt->index_type. If yes, interpolate it at the
	   right values of k. */

	if (tp_of_tt[index_tt] != previous_type) {
	  
	  if (ptr->transfer_verbose>2)
	    printf("In %s: Interpolate source for one mode/ic/type.\n",
		   __func__);
	  
	  class_call(transfer_interpolate_sources(ppt,
						  ptr,
						  index_md,
						  index_ic,
						  tp_of_tt[index_tt],
						  source_spline,
						  interpolated_sources),
		     ptr->error_message,
		     ptr->error_message);
	}
	
	previous_type = tp_of_tt[index_tt];

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
  shared (pw,ptr,ppr,ppt,index_md,index_ic,index_tt,			\
	  interpolated_sources,abort,num_j,tau0,tau_rec,tau_size_max)	\
  private (thread,index_l,tstart,tstop,tspent,address_in_workspace,tau0_minus_tau,delta_tau,sources,j_l,ddj_l,x_size_l,x_min_l,k_max_bessel,tau_size)
	
	{
	  
#ifdef _OPENMP
	  thread = omp_get_thread_num();
	  tspent = 0.;
#endif
	  
	  /* define address of each field in the workspace */
	  address_in_workspace = pw[thread];
	  
	  tau0_minus_tau = address_in_workspace;
	  address_in_workspace += tau_size_max;
	  
	  delta_tau  = address_in_workspace;
	  address_in_workspace += tau_size_max;
	  
	  tau_size = address_in_workspace;
	  address_in_workspace += 1;

	  sources = address_in_workspace;
	  address_in_workspace += ptr->k_size[index_md]*tau_size_max;
	  
	  x_min_l = address_in_workspace;
	  address_in_workspace += 1;
	  
	  j_l = address_in_workspace;
      
	  /* the address of the next fields, ddj_l, etc., will be defined
	     within the l loop, since they depend on l */

	  /* the code makes a distinction between "perturbation
	     sources" (e.g. gravitational potential) and "transfer
	     sources" (e.g. total density fluctuations, obtained
	     through the Poisson equation, and observed with a given
	     selection function).

	     The next routine computes the transfer source given the
	     interpolated perturbation source, and copies it in the
	     workspace. */

	  class_call_parallel(transfer_sources(ppr,
					       pba,
					       ppt,
					       ptr,
					       interpolated_sources,
					       tau_rec,
					       index_md,
					       index_tt,
					       sources,
					       tau0_minus_tau,
					       delta_tau,
					       tau_size),
			      ptr->error_message,
			      ptr->error_message);

#pragma omp for schedule (dynamic)

	  for (index_l = 0; index_l < ptr->l_size[index_md]; index_l++) {

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

	    /* for a given l, maximum value of k such that we can convolve
	       the source with Bessel functions j_l(x) without reaching x_max */
	    
	    k_max_bessel = ((*x_min_l)+(x_size_l-1)*pbs->x_step)/tau0_minus_tau[0];

	    /* compute the transfer function for this l */
	    class_call_parallel(transfer_compute_for_each_l(ppr,
							    ppt,
							    ptr,
							    index_md,
							    index_ic,
							    index_tt,
							    index_l,
							    (double)ptr->l[index_l],
							    *x_min_l,
							    pbs->x_step,
							    tau0_minus_tau,
							    delta_tau,
							    (int)(*tau_size),
							    sources,
							    j_l,
							    ddj_l,
							    dj_l,
							    dddj_l,
							    k_max_bessel),
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
    free(tp_of_tt);

#pragma omp parallel shared(pw) private(thread)
    {
#ifdef _OPENMP
      thread = omp_get_thread_num();
#endif
      free(pw[thread]);
    }
    
  } /* end of loop over mode */

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

  int index_md;

  if (ptr->has_cls == _TRUE_) {

    for (index_md = 0; index_md < ptr->md_size; index_md++) {
      free(ptr->l_size_tt[index_md]);
      free(ptr->k[index_md]);
      free(ptr->transfer[index_md]);
    }  
   
    free(ptr->tt_size);
    free(ptr->l_size_tt);
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
				  double tau0
				  ) {

  /** Summary: */

  /** - define local variables */

  int index_md,index_tt,index_tt_common;

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

    if (ppt->has_cl_density == _TRUE_) {
      ptr->index_tt_density = index_tt;
      index_tt+=ppt->selection_num;
    }

    if (ppt->has_cl_lensing_potential == _TRUE_) {
      ptr->index_tt_lensing = index_tt;
      index_tt+=ppt->selection_num;
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

  /* number of l values for each mode and type,
     l_size_tt[index_md][index_tt], and maximized for each mode,
     l_size[index_md] */

  class_alloc(ptr->l_size,ptr->md_size * sizeof(int),ptr->error_message);

  class_alloc(ptr->l_size_tt,ptr->md_size * sizeof(int *),ptr->error_message);

  for (index_md = 0; index_md < ptr->md_size; index_md++) {
    class_alloc(ptr->l_size_tt[index_md],ptr->tt_size[index_md] * sizeof(int),ptr->error_message);
  }

  /* number of k values for each mode, k_size[index_md] */

  class_alloc(ptr->k_size,ptr->md_size * sizeof(int),ptr->error_message);

  /* list of k values for each mode, k[index_md] */

  class_alloc(ptr->k,ptr->md_size * sizeof(double *),ptr->error_message);

  /* array (of array) of transfer functions for each mode, transfer[index_md] */

  class_alloc(ptr->transfer,ptr->md_size * sizeof(double *),ptr->error_message);

  /** get l values using transfer_get_l_list() */
  class_call(transfer_get_l_list(ppr,ppt,pbs,ptr),
	     ptr->error_message,
	     ptr->error_message);
  
  /** - loop over modes (scalar, etc). For each mode: */
  
  for (index_md = 0; index_md < ptr->md_size; index_md++) {

    /** (a) get k values using transfer_get_k_list() */
    class_call(transfer_get_k_list(ppr,ppt,ptr,tau0,index_md),
	       ptr->error_message,
	       ptr->error_message);

    /** (b) allocate arrays of transfer functions, (ptr->transfer[index_md])[index_ic][index_tt][index_l][index_k] */
    class_alloc(ptr->transfer[index_md],
		ppt->ic_size[index_md] * ptr->tt_size[index_md] * ptr->l_size[index_md] * ptr->k_size[index_md] * sizeof(double),
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
  int index_md;
  int index_tt;

  /* allocate and fill l array (taken directly from Bessel module) */

  ptr->l_size_max = pbs->l_size;

  class_alloc(ptr->l,ptr->l_size_max*sizeof(int),ptr->error_message);
  
  for (index_l=0; index_l < ptr->l_size_max; index_l++)
    ptr->l[index_l]=pbs->l[index_l];

  /* for each mode and type, find relevant size of l array,
     l_size_tt[index_md][index_tt] (since for some modes and types
     l_max can be smaller). Also, maximize this size for each mode to
     find l_size[index_md]. */

  for (index_md=0; index_md < ppt->md_size; index_md++) {

    ptr->l_size[index_md] = 0;

    for (index_tt=0;index_tt<ptr->tt_size[index_md];index_tt++) {

      if _scalars_ {

	if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t))
	  l_max=ppt->l_scalar_max;

	if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_e))
	  l_max=ppt->l_scalar_max;

	if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (index_tt == ptr->index_tt_lcmb))
	  l_max=ppt->l_scalar_max;

	if ((ppt->has_cl_density == _TRUE_) && (index_tt >= ptr->index_tt_density) && (index_tt < ptr->index_tt_density+ppt->selection_num))
	  l_max=ppt->l_lss_max;
	
	if ((ppt->has_cl_lensing_potential == _TRUE_) && (index_tt >= ptr->index_tt_lensing) && (index_tt < ptr->index_tt_lensing+ppt->selection_num))
	  l_max=ppt->l_lss_max;
	
      }
      
      if _tensors_
	l_max = ppt->l_tensor_max;

      class_test(l_max > ptr->l[ptr->l_size_max-1],
		 ptr->error_message,
		 "For mode %d, type %d, asked for l_max=%d greater than in Bessel table where l_max=%d",
		 index_md,
		 index_tt,
		 l_max,
		 ptr->l[ptr->l_size_max-1]);

      index_l=0;
      while (ptr->l[index_l] < l_max) index_l++;
      ptr->l_size_tt[index_md][index_tt]=index_l+1;

      if (ptr->l_size_tt[index_md][index_tt] < ptr->l_size_max)
	ptr->l_size_tt[index_md][index_tt]++;
      if (ptr->l_size_tt[index_md][index_tt] < ptr->l_size_max)
	ptr->l_size_tt[index_md][index_tt]++;

      ptr->l_size[index_md] = max(ptr->l_size[index_md],ptr->l_size_tt[index_md][index_tt]);

    }
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
 * @param index_md Input: index of requested mode (scalar, tensor, etc) 
 * @return the error status
 */

int transfer_get_k_list(
			struct precision * ppr,
			struct perturbs * ppt,
			struct transfers * ptr,
			double tau0,
			int index_md
			) {

  int index_k_pt;
  int index_k_tr;
  double k,k_min,k_max,k_step_max=0.;

  /* find k_step_max, the maximum value of the step */

  if _scalars_ {

    k_step_max = 2.*_PI_/tau0*ppr->k_step_trans_scalars;

  }

  if _tensors_ {

    k_step_max = 2.*_PI_/tau0*ppr->k_step_trans_tensors;

  }

  class_test(k_step_max == 0.,
	     ptr->error_message,
	     "stop to avoid infinite loop");

  /* first and last value in perturbation module */

  k_min = ppt->k[index_md][0]; /* first value, inferred from perturbations structure */

  k_max = ppt->k[index_md][ppt->k_size_cl[index_md]-1]; /* last value, inferred from perturbations structure */

  /* first, count the number of necessary values */

  index_k_pt = 0;
  index_k_tr = 0;

  /* - first point */

  k = k_min;
  index_k_pt++;
  index_k_tr++;

  /* - points taken from perturbation module if step small enough */

  while ((index_k_pt < ppt->k_size[index_md]) && ((ppt->k[index_md][index_k_pt] -k) < k_step_max)) {
    k = ppt->k[index_md][index_k_pt];
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
    ptr->k_size[index_md]=index_k_tr-1;
  else
    ptr->k_size[index_md]=index_k_tr;

  class_alloc(ptr->k[index_md],
	      ptr->k_size[index_md]*sizeof(double),
	      ptr->error_message);

  /* repeat exactly the same steps, but now filling the list */

  index_k_pt = 0;
  index_k_tr = 0;

  ptr->k[index_md][0] = k_min;
  k = k_min;
  index_k_pt++;
  index_k_tr++;

  while ((index_k_pt < ppt->k_size[index_md]) && ((ppt->k[index_md][index_k_pt] -k) < k_step_max)) {
    k = ppt->k[index_md][index_k_pt];
    ptr->k[index_md][index_k_tr] = k;
    index_k_pt++;
    index_k_tr++;
  }

  while ((index_k_tr < ptr->k_size[index_md]) && (k < k_max)) {
    k += k_step_max;
    ptr->k[index_md][index_k_tr] = k;
    index_k_tr++;
  }

  /* consistency check */

  class_test(ptr->k[index_md][ptr->k_size[index_md]-1] > k_max,
	     ptr->error_message,
	     "bug in k list calculation, k_max larger in transfer than in perturb, should never happen");

  return _SUCCESS_; 

}

/**
 * This routine defines the correspondence between the sources in the
 * perturbation and transfer module.
 *
 * @param ppt  Input : pointer to perturbation structure
 * @param ptr  Input : pointer to transfers structure containing l's
 * @param index_md : Input: index of mode (scalar, tensor...)
 * @param tp_of_tt : Input/Output: array with the correspondance (allocated before, filled here)
 * @return the error status
 */

int transfer_get_source_correspondence(
				       struct perturbs * ppt,
				       struct transfers * ptr,
				       int index_md,
				       int * tp_of_tt
				       ) {

  /* running index on transfer types */
  int index_tt;

  /** - which source are we considering? Define correspondence
      between transfer types and source types */
  
  for (index_tt=0; index_tt<ptr->tt_size[index_md]; index_tt++) {
    
    if _scalars_ {
      
      if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t)) 
	tp_of_tt[index_tt]=ppt->index_tp_t;
      
      if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_e)) 
	tp_of_tt[index_tt]=ppt->index_tp_e;
      
      if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (index_tt == ptr->index_tt_lcmb)) 
	tp_of_tt[index_tt]=ppt->index_tp_g;
      
      if ((ppt->has_cl_density == _TRUE_) && (index_tt >= ptr->index_tt_density) && (index_tt < ptr->index_tt_density+ppt->selection_num))
	tp_of_tt[index_tt]=ppt->index_tp_g;

      if ((ppt->has_cl_lensing_potential == _TRUE_) && (index_tt >= ptr->index_tt_lensing) && (index_tt < ptr->index_tt_lensing+ppt->selection_num))
	tp_of_tt[index_tt]=ppt->index_tp_g;
      
    }
    
    if _tensors_ {
      
      if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t)) 
	tp_of_tt[index_tt]=ppt->index_tp_t;
      
      if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_e)) 
	tp_of_tt[index_tt]=ppt->index_tp_e;
      
      if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_b))
	tp_of_tt[index_tt]=ppt->index_tp_b;
    }

  }
  
  return _SUCCESS_;
  
}

/**
 * the code makes a distinction between "perturbation sources"
 * (e.g. gravitational potential) and "transfer sources" (e.g. total
 * density fluctuations, obtained through the Poisson equation, and
 * observed with a given selection function).
 *
 * This routine computes the number of sampled time values for each type
 * of transfer sources.
 *
 * @param ppr                   Input : pointer to precision structure
 * @param pba                   Input : pointer to background structure
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param tau_rec               Input : recombination time
 * @param tau0                  Input : time today
 * @param index_md            Input : index of the mode (scalar, tensor)
 * @param index_tt              Input : index of transfer type
 * @param tau_size              Output: pointer to number of smapled times
 * @return the error status
 */

int transfer_source_tau_size(
			     struct precision * ppr,
			     struct background * pba,
			     struct perturbs * ppt,
			     struct transfers * ptr,
			     double tau_rec,
			     double tau0,
			     int index_md,
			     int index_tt,
			     int * tau_size) {

  /* values of conformal time */
  double tau_min,tau_mean,tau_max;
  
  /* minimum value of index_tt */
  int index_tau_min;
  
  /* value of l at which limber approximation is switched on */
  int l_limber;
  
  /* scalar mode */
  if _scalars_ {

    /* scalar temperature */
    if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t))
      *tau_size = ppt->tau_size;

    /* scalar polarisation */
    if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_e))
      *tau_size = ppt->tau_size;

    /* cmb lensing potential */
    if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (index_tt == ptr->index_tt_lcmb)) {

      /* find times before recombination, that will be thrown away */
      index_tau_min=0;
      while (ppt->tau_sampling[index_tau_min]<=tau_rec) index_tau_min++;
      
      /* infer number of time steps after removing early times */
      *tau_size = ppt->tau_size-index_tau_min;
    }

    /* density Cl's */
    if ((ppt->has_cl_density == _TRUE_) && (index_tt >= ptr->index_tt_density) && (index_tt < ptr->index_tt_density+ppt->selection_num)) {

      /* time interval for this bin */
      class_call(transfer_selection_times(ppr,
					  pba,
					  ppt,
					  ptr,
					  index_tt-ptr->index_tt_density,
					  &tau_min,
					  &tau_mean,
					  &tau_max),
		 ptr->error_message,
		 ptr->error_message);
      
      /* case selection=dirac */
      if (tau_min == tau_max) {
	*tau_size = 1;
      }
      /* other cases (gaussian, top-hat...) */
      else {

	/* check that selection function well sampled */
	*tau_size = (int)ppr->selection_sampling;

	/* value of l at which the code switches to Limber approximation
	   (necessary for next step) */
	l_limber=ppr->l_switch_limber_for_cl_density_over_z*ppt->selection_mean[index_tt-ptr->index_tt_density];
	      
	/* check that bessel well sampled, if not define finer sampling
	   overwriting the previous one.
	   One Bessel oscillations corresponds to [Delta tau]=2pi/k.
	   This is minimal for largest relevant k_max, 
	   namely k_max=l_limber/(tau0-tau_mean).
	   We need to cut the interval (tau_max-tau_min) in pieces of size
           [Delta tau]=2pi/k_max. This gives the number below.
	*/
	*tau_size=max(*tau_size,(int)((tau_max-tau_min)/((tau0-tau_mean)/l_limber))*ppr->selection_sampling_bessel);
      }
    }

    /* galaxy lensing Cl's, differs from density Cl's since the source
       function will spead from the selection function region up to
       tau0 */
    if ((ppt->has_cl_lensing_potential == _TRUE_) && (index_tt >= ptr->index_tt_lensing) && (index_tt < ptr->index_tt_lensing+ppt->selection_num)) {

      /* time interval for this bin */
      class_call(transfer_selection_times(ppr,
					  pba,
					  ppt,
					  ptr,
					  index_tt-ptr->index_tt_lensing,
					  &tau_min,
					  &tau_mean,
					  &tau_max),
		 ptr->error_message,
		 ptr->error_message);

      /* check that selection function well sampled */
      *tau_size = (int)ppr->selection_sampling;
      
      /* value of l at which the code switches to Limber approximation
	 (necessary for next step) */
      l_limber=ppr->l_switch_limber_for_cl_density_over_z*ppt->selection_mean[index_tt-ptr->index_tt_lensing];
      
      /* check that bessel well sampled, if not define finer sampling
	 overwriting the previous one.
	 One Bessel oscillations corresponds to [Delta tau]=2pi/k.
	 This is minimal for largest relevant k_max, 
	 namely k_max=l_limber/((tau0-tau_mean)/2).
	 We need to cut the interval (tau_0-tau_min) in pieces of size
	 [Delta tau]=2pi/k_max. This gives the number below. 
      */
      *tau_size=max(*tau_size,(int)((tau0-tau_min)/((tau0-tau_mean)/2./l_limber))*ppr->selection_sampling_bessel);

    }
  }

  /* tensor mode */
  if _tensors_ {

    /* for all tensor types */
    *tau_size = ppt->tau_size;
  }

  return _SUCCESS_;
}

/**
 * This routine interpolates sources \f$ S(k, \tau) \f$ for each mode,
 * initial condition and type (of perturbation module), to get them at
 * the right values of k, using the spline interpolation method.
 *
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param index_md            Input : index of mode
 * @param index_ic              Input : index of initial condition
 * @param index_type            Input : index of type of source (in perturbation module)
 * @param source_spline         Output: array of second derivative of sources (filled here but allocated in transfer_init() to avoid numerous reallocation)
 * @param interpolated_sources  Output: array of interpolated sources (filled here but allocated in transfer_init() to avoid numerous reallocation)
 * @return the error status
 */

int transfer_interpolate_sources(
				 struct perturbs * ppt,
				 struct transfers * ptr,
				 int index_md,
				 int index_ic,
				 int index_type,
				 double * source_spline, /* array with argument source_spline[index_tau*ppt->k_size[index_md]+index_k] (must be allocated) */
				 double * interpolated_sources /* array with argument interpolated_sources[index_k_tr*ppt->tau_size+index_tau] (must be allocated) */
				 ) {

  /** Summary: */

  /** - define local variables */

  /* index running on k values in the original source array */
  int index_k;

  /* index running on time */
  int index_tau;

  /* index running on k values in the interpolated source array */
  int index_k_tr;

  /* variables used for spline interpolation algorithm */
  double h, a, b;

  /** - find second derivative of original sources with respect to k
      in view of spline interpolation */

  class_call(array_spline_table_columns(ppt->k[index_md],
					ppt->k_size[index_md],
					ppt->sources[index_md][index_ic * ppt->tp_size[index_md] + index_type],
					ppt->tau_size,
					source_spline,
					_SPLINE_EST_DERIV_,
					ptr->error_message),
	     ptr->error_message,
	     ptr->error_message);

  /** - interpolate at each k value using the usual 
      spline interpolation algorithm. */
  
  index_k = 0;
  h = ppt->k[index_md][index_k+1] - ppt->k[index_md][index_k];
  
  for (index_k_tr = 0; index_k_tr < ptr->k_size[index_md]; index_k_tr++) {
    
    while (((index_k+1) < ppt->k_size[index_md]) &&
	   (ppt->k[index_md][index_k+1] < 
	    ptr->k[index_md][index_k_tr])) {
      index_k++;
      h = ppt->k[index_md][index_k+1] - ppt->k[index_md][index_k];
    }
    
    class_test(h==0.,
	       ptr->error_message,
	       "stop to avoid division by zero");
    
    b = (ptr->k[index_md][index_k_tr] - ppt->k[index_md][index_k])/h;
    a = 1.-b;
    
    for (index_tau = 0; index_tau < ppt->tau_size; index_tau++) {
      
      interpolated_sources[index_k_tr*ppt->tau_size+index_tau] = 
	a * ppt->sources[index_md]
	[index_ic * ppt->tp_size[index_md] + index_type]
	[index_tau*ppt->k_size[index_md]+index_k]
	+ b * ppt->sources[index_md]
	[index_ic * ppt->tp_size[index_md] + index_type]
	[index_tau*ppt->k_size[index_md]+index_k+1]
	+ ((a*a*a-a) * source_spline[index_tau*ppt->k_size[index_md]+index_k]
	   +(b*b*b-b) * source_spline[index_tau*ppt->k_size[index_md]+index_k+1])*h*h/6.0;
      
    }
  }

  return _SUCCESS_;
  
}

/**
 * the code makes a distinction between "perturbation sources"
 * (e.g. gravitational potential) and "transfer sources" (e.g. total
 * density fluctuations, obtained through the Poisson equation, and
 * observed with a given selection function).
 *
 * This routine computes the transfer source given the interpolated
 * perturbation source, and copies it in the workspace.
 *
 * @param ppr                   Input : pointer to precision structure
 * @param pba                   Input : pointer to background structure
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param interpolated_sources  Input : interpolated perturbation source
 * @param tau_rec               Input : recombination time
 * @param index_md            Input : index of mode
 * @param index_tt              Input : index of type of (transfer) source
 * @param sources               Output: transfer source
 * @param tau0_minus_tau        Output: values of (tau0-tau) at which source are sample
 * @param delta_tau             Output: corresponding values dtau, used later for trapezoidal integration
 * @param tau_size_double       Output: pointer to size of previous two arrays, converted to double
 * @return the error status
 */

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
		     )  {

  /** Summary: */
  
  /** - define local variables */

  /* index running on time */
  int index_tau;

  /* index running on k values in the interpolated source array */
  int index_k_tr;

  /* bin for computation of cl_density */  
  int bin;

  /* number of tau values */
  int tau_size;

  /* minimum tau index kept in transfer sources */
  int index_tau_min;

  /* for calling background_at_eta */
  int last_index;
  double * pvecback = NULL;

  /* conformal time */
  double tau, tau0;

  /* rescaling factor depending on the background at a given time */
  double rescaling;

  /* flag: is there any difference between the perturbation and transfer source? */
  short redefine_source;

  /* array of selection function values at different times */
  double * selection;
  
  /* array of time sampling for lensing source selection function */
  double * tau0_minus_tau_lensing_sources;
  
  /* differential array of time sampling for lensing source selection function */
  double * delta_tau_lensing_sources;

  /* index running on time in previous two arrays */
  int index_tau_sources;

  /* number of time values in previous two arrays */
  int tau_sources_size;

  /* in which cases are perturbation and transfer sources are different?
     I.e., in which case do we need to mutiply the sources by some
     background and/or window function, and eventually to resample it,
     or redfine its time limits? */

  redefine_source = _FALSE_;
  
  if _scalars_ {

    /* cmb lensing potential */
    if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (index_tt == ptr->index_tt_lcmb))
      redefine_source = _TRUE_;
  
    /* density Cl's */
    if ((ppt->has_cl_density == _TRUE_) && (index_tt >= ptr->index_tt_density) && (index_tt < ptr->index_tt_density+ppt->selection_num))
      redefine_source = _TRUE_;

    /* galaxy lensing potential */
    if ((ppt->has_cl_lensing_potential == _TRUE_) && (index_tt >= ptr->index_tt_lensing) && (index_tt < ptr->index_tt_lensing+ppt->selection_num))
      redefine_source = _TRUE_;

  }

  /* conformal time today */
  tau0 = pba->conformal_age;

  /* case where we need to redefine by a window function (or any
     function of the background and of k) */
  if (redefine_source == _TRUE_) {
    
    class_call(transfer_source_tau_size(ppr,
					pba,
					ppt,
					ptr,
					tau_rec,
					tau0,
					index_md,
					index_tt,
					&tau_size),
	       ptr->error_message,
	       ptr->error_message);

    if _scalars_ {

      /* lensing source: throw away times before recombuination, and multiply psi by window function */

      if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (index_tt == ptr->index_tt_lcmb)) {
      
        /* first time step after removing early times */
	index_tau_min =  ppt->tau_size - tau_size;

        /* loop over time and rescale */
	for (index_tau = index_tau_min; index_tau < ppt->tau_size; index_tau++) {

	  /* conformal time */
	  tau = ppt->tau_sampling[index_tau];

	  /* lensing source =  - 2 W(tau) psi(k,tau) Heaviside(tau-tau_rec) 
	     with 
	     psi = (newtonian) gravitationnal potential  
	     W = (tau-tau_rec)/(tau_0-tau)/(tau_0-tau_rec) 
	     H(x) = Heaviside
	     (in tau = tau_0, set source = 0 to avoid division by zero;
	     regulated anyway by Bessel).
	  */
		  
	  if (index_tau == ppt->tau_size-1) {
	    rescaling=0.;
	  }
	  else {
	    rescaling = -2.*(tau-tau_rec)/(tau0-tau)/(tau0-tau_rec);
	  }
	  
          /* copy from input array to output array */
	  for (index_k_tr = 0; index_k_tr < ptr->k_size[index_md]; index_k_tr++) { 
	    sources[index_k_tr*tau_size+(index_tau-index_tau_min)] = 
	      interpolated_sources[index_k_tr*ppt->tau_size+index_tau]
	      * rescaling
	      * ptr->lcmb_rescale
	      * pow(ptr->k[index_md][index_k_tr]/ptr->lcmb_pivot,ptr->lcmb_tilt);
	  }
	  
          /* store value of (tau0-tau) */
	  tau0_minus_tau[index_tau-index_tau_min] = tau0 - tau;
	  
	}

        /* infer values of delta_tau from values of (tau0-tau) */
	class_call(transfer_integration_time_steps(ptr,
						   tau0_minus_tau,
						   tau_size,
						   delta_tau),
		   ptr->error_message,
		   ptr->error_message);
	
      }
    
      /* density source: redefine the time sampling, multiply by
	 coefficient of Poisson equation, and multiply by selection
	 function */

      if ((ppt->has_cl_density == _TRUE_) && (index_tt >= ptr->index_tt_density) && (index_tt < ptr->index_tt_density+ppt->selection_num)) {
      
        /* bin number associated to particular redshift bin and selection function */
	bin=index_tt-ptr->index_tt_density;

        /* allocate temporary arrays for storing sources and for calling background */
	class_alloc(selection,tau_size*sizeof(double),ptr->error_message);
	class_alloc(pvecback,pba->bg_size*sizeof(double),ptr->error_message); 

        /* redefine the time sampling */
	class_call(transfer_selection_sampling(ppr,
					       pba,
					       ppt,
					       ptr,
					       bin,
					       tau0_minus_tau,
					       tau_size),
		   ptr->error_message,
		   ptr->error_message);
	
        /* resample the source at those times */
	class_call(transfer_source_resample(ppr,
					    pba,
					    ppt,
					    ptr,
					    bin,
					    tau0_minus_tau,
					    tau_size,
					    index_md,
					    tau0,
					    interpolated_sources,
					    sources),
		   ptr->error_message,
		   ptr->error_message);

	/* infer values of delta_tau from values of (tau0-tau) */
	class_call(transfer_integration_time_steps(ptr,
						   tau0_minus_tau,
						   tau_size,
						   delta_tau),
		   ptr->error_message,
		   ptr->error_message);

	/* compute values of selection function at sampled values of tau */
	class_call(transfer_selection_compute(ppr,
					      pba,
					      ppt,
					      ptr,
					      selection,
					      tau0_minus_tau,
					      delta_tau,
					      tau_size,
					      pvecback,
					      tau0,
					      bin),
		   ptr->error_message,
		   ptr->error_message);

	/* loop over time and rescale */
	for (index_tau = 0; index_tau < tau_size; index_tau++) {
	
	  /* conformal time */
	  tau = tau0 - tau0_minus_tau[index_tau];

	  /* corresponding background quantities */
	  class_call(background_at_tau(pba,
				       tau,
				       pba->long_info,
				       pba->inter_normal,
				       &last_index,
				       pvecback),
		     pba->error_message,
		     ptr->error_message);

	  /* matter density source =  - [- (dz/dtau) W(z)] * 2/(3 Omega_m(tau) H^2(tau)) * (k/a)^2 psi(k,tau)
	     =  - W(tau) * 2/(3 Omega_m(tau) H^2(tau)) * (k/a)^2 psi(k,tau)
	     with 
	     psi = (newtonian) gravitationnal potential  
	     W(z) = redshift space selection function = dN/dz
	     W(tau) = same wrt conformal time = dN/dtau
	     (in tau = tau_0, set source = 0 to avoid division by zero;
	     regulated anyway by Bessel).
	  */

	  rescaling = selection[index_tau]
	    *(-2.)/3./pvecback[pba->index_bg_Omega_m]/pvecback[pba->index_bg_H]
	    /pvecback[pba->index_bg_H]/pow(pvecback[pba->index_bg_a],2);
	  
	  for (index_k_tr = 0; index_k_tr < ptr->k_size[index_md]; index_k_tr++) {
	    sources[index_k_tr*tau_size+index_tau] *= rescaling*pow(ptr->k[index_md][index_k_tr],2);
	  }
	}

	/*
	  if (bin == 2) {
	  for (index_k_tr = 0; index_k_tr < ptr->k_size[index_md]; index_k_tr++) {
	  for (index_tau = 0; index_tau < tau_size; index_tau++) {
	  fprintf(stdout,"%e %e %e\n",
	  tau0-tau0_minus_tau[index_tau],
	  ptr->k[index_md][index_k_tr],
	  sources[index_k_tr*tau_size+index_tau]);
	  }
	  fprintf(stdout,"\n\n");
	  }
	  }
	*/

	/* deallocate temporary arrays */
	free(pvecback);
	free(selection);
      }

      /* lensing potential: eliminate early times, and multiply by selection
	 function */

      if ((ppt->has_cl_lensing_potential == _TRUE_) && (index_tt >= ptr->index_tt_lensing) && (index_tt < ptr->index_tt_lensing+ppt->selection_num)) {
	
        /* bin number associated to particular redshift bin and selection function */
	bin=index_tt-ptr->index_tt_lensing;

        /* allocate temporary arrays for storing sources and for calling background */
	class_alloc(pvecback,
		    pba->bg_size*sizeof(double),
		    ptr->error_message); 

	/* dirac case */
	if (ppt->selection == dirac) {
	  tau_sources_size=1;
	}
	/* other cases (gaussian, tophat...) */
	else {
	  tau_sources_size=ppr->selection_sampling;
	}

	class_alloc(selection,
		    tau_sources_size*sizeof(double),
		    ptr->error_message);	  
	
	class_alloc(tau0_minus_tau_lensing_sources,
		    tau_sources_size*sizeof(double),
		    ptr->error_message);
	
	class_alloc(delta_tau_lensing_sources,
		    tau_sources_size*sizeof(double),
		    ptr->error_message);
	
	/* time sampling for source selection function */
	class_call(transfer_selection_sampling(ppr,
					       pba,
					       ppt,
					       ptr,
					       bin,
					       tau0_minus_tau_lensing_sources,
					       tau_sources_size),
		   ptr->error_message,
		   ptr->error_message);
	
	/* infer values of delta_tau from values of (tau0-tau) */
	class_call(transfer_integration_time_steps(ptr,
						   tau0_minus_tau_lensing_sources,
						   tau_sources_size,
						   delta_tau_lensing_sources),
		   ptr->error_message,
		   ptr->error_message);
	
	/* compute values of selection function at sampled values of tau */
	class_call(transfer_selection_compute(ppr,
					      pba,
					      ppt,
					      ptr,
					      selection,
					      tau0_minus_tau_lensing_sources,
					      delta_tau_lensing_sources,
					      tau_sources_size,
					      pvecback,
					      tau0,
					      bin),
		   ptr->error_message,
		   ptr->error_message);

	/* redefine the time sampling */
	class_call(transfer_lensing_sampling(ppr,
					     pba,
					     ppt,
					     ptr,
					     bin,
					     tau0, 
					     tau0_minus_tau,
					     tau_size),
		   ptr->error_message,
		   ptr->error_message);
	
        /* resample the source at those times */
	class_call(transfer_source_resample(ppr,
					    pba,
					    ppt,
					    ptr,
					    bin,
					    tau0_minus_tau,
					    tau_size,
					    index_md,
					    tau0,
					    interpolated_sources,
					    sources),
		   ptr->error_message,
		   ptr->error_message);

	/* infer values of delta_tau from values of (tau0-tau) */
	class_call(transfer_integration_time_steps(ptr,
						   tau0_minus_tau,
						   tau_size,
						   delta_tau),
		   ptr->error_message,
		   ptr->error_message);

        /* loop over time and rescale */
	for (index_tau = 0; index_tau < tau_size; index_tau++) {

	  /* lensing source =  - 2 W(tau) psi(k,tau) Heaviside(tau-tau_rec) 
	     with 
	     psi = (newtonian) gravitationnal potential  
	     W = (tau-tau_rec)/(tau_0-tau)/(tau_0-tau_rec) 
	     H(x) = Heaviside
	     (in tau = tau_0, set source = 0 to avoid division by zero;
	     regulated anyway by Bessel).
	  */
		  
	  if (index_tau == tau_size-1) {
	    rescaling=0.;
	  }
	  else {

	    rescaling = 0.;

	    for (index_tau_sources=0; 
		 index_tau_sources < tau_sources_size; 
		 index_tau_sources++) {
	     
	      /* condition for excluding from the sum the sources located in z=zero */
	      if ((tau0_minus_tau_lensing_sources[index_tau_sources] > 0.) && (tau0_minus_tau_lensing_sources[index_tau_sources]-tau0_minus_tau[index_tau] > 0.)) {

		rescaling +=  
		  -2.*(tau0_minus_tau_lensing_sources[index_tau_sources]-tau0_minus_tau[index_tau])
		  /tau0_minus_tau[index_tau]
		  /tau0_minus_tau_lensing_sources[index_tau_sources]
		  * selection[index_tau_sources]
		  * delta_tau_lensing_sources[index_tau_sources];
	      }
	    }
	      
	    rescaling /= 2.;

	  }
	  
          /* copy from input array to output array */
	  for (index_k_tr = 0; index_k_tr < ptr->k_size[index_md]; index_k_tr++) { 
	    sources[index_k_tr*tau_size+index_tau] *= rescaling;
	  }
	  
	}

	/* deallocate temporary arrays */
	free(pvecback);
	free(selection);
	free(tau0_minus_tau_lensing_sources);
	free(delta_tau_lensing_sources);
	
      }
    }
  }

  /* case where we do not need to redefine */

  else {

    /* number of sampled time values */
    tau_size = ppt->tau_size;
    
    /* plain copy from input array to output array */
    memcpy(sources,
	   interpolated_sources,
	   ptr->k_size[index_md]*ppt->tau_size*sizeof(double));
        
    /* store values of (tau0-tau) */
    for (index_tau=0; index_tau < ppt->tau_size; index_tau++) {
      tau0_minus_tau[index_tau] = tau0 - ppt->tau_sampling[index_tau];
    }

    /* infer values of delta_tau from values of (tau0-tau) */
    class_call(transfer_integration_time_steps(ptr,
					       tau0_minus_tau,
					       tau_size,
					       delta_tau),
	       ptr->error_message,
	       ptr->error_message);
  }
  
  /* return tau_size value that will be stored in the workspace (the
     workspace wants a double) */

  *tau_size_double = (double)tau_size;

  return _SUCCESS_;
    
}

/**
 * infer delta_tau array from tau0_minus_tau array (will be used in
 * transfer_integrate for a fast trapezoidal integration
 *
 * @param ptr                   Input: pointer to transfers structure
 * @param tau0_minus_tau        Input: values of (tau0-tau) at which source are sample
 * @param tau_size              Input: size of previous array
 * @param delta_tau             Output: corresponding values delta_tau
 * @return the error status
 */

int transfer_integration_time_steps(
				    struct transfers * ptr,
				    double * tau0_minus_tau,
				    int tau_size,
				    double * delta_tau
				    ) {

  /* running index on time */
  int index_tau;

  /* case selection = dirac */
  if (tau_size == 1) {
    delta_tau[0] = 2.; // factor 2 corrects for dibision by two occuring by convention in all our trapezoidal integrations 
  }
  /* other cases */
  else {

    /* first value */
    delta_tau[0] = tau0_minus_tau[0]-tau0_minus_tau[1];
    
    /* intermediate values */
    for (index_tau=1; index_tau < tau_size-1; index_tau++)
      delta_tau[index_tau] = tau0_minus_tau[index_tau-1]-tau0_minus_tau[index_tau+1];
    
    /* last value */
    delta_tau[tau_size-1] = tau0_minus_tau[tau_size-2]-tau0_minus_tau[tau_size-1];
    
  }

  return _SUCCESS_;
  
}

/**
 * arbitrarily normalized selection function dN/dz(z,bin)
 *
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param bin                   Input : redshift bin number
 * @param z                     Input : one value of redshift
 * @param selection             Output: pointer to selection function
 * @return the error status
 */

int transfer_selection_function(
				struct precision * ppr,
				struct perturbs * ppt,
				struct transfers * ptr,
				int bin,
				double z,
				double * selection) {

  double x;

  /* trivial dirac case */
  if (ppt->selection==dirac) {

    *selection=1.;

    return _SUCCESS_;
  }

  /* difference between z and the bin center (we can take the absolute
     value as long as all selection functions are symmetric around
     x=0) */
  x=fabs(z-ppt->selection_mean[bin]);

  /* gaussian case (the function is anyway normalized later
     automatically, but could not resist to normalize it already
     here) */
  if (ppt->selection==gaussian) {

    *selection = exp(-0.5*pow(x/ppt->selection_width[bin],2))
      /ppt->selection_width[bin]/sqrt(2.*_PI_);

    return _SUCCESS_;
  }
  
  /* top-hat case, with smoothed edges. The problem with sharp edges
     is that the final result will be affected by random
     noise. Indeed, the values of k at which the transfer functions
     Delta_l(k) are sampled will never coicide with the actual edges
     of the true transfer function (computed with or even without the
     Limber approximation). Hence the integral Cl=\int dk
     Delta_l(k)**2 (...) will be unprecise and will fluctuate randomly
     with the resolution along k. With smooth edges, the problemeis
     sloved, and the final Cls become mildly dependent on the
     resolution along k. */
  if (ppt->selection==tophat) {
    
    /* selection function, centered on z=mean (i.e. on x=0), equal to
       one around x=0, with tanh step centered on x=width, of width
       delta x = 0.1*width
    */
    *selection=(1.-tanh((x-ppt->selection_width[bin])/(ppr->selection_tophat_edge*ppt->selection_width[bin])))/2.;
  
    return _SUCCESS_;
  }

  /* get here only if selection type was recognized */
  class_stop(ptr->error_message,
	     "invalid choice of selection function");

  return _SUCCESS_;
}

/**
 * for sources that need to be mutiplied by a selection function,
 * redefine a finer time sampling in a small range
 *
 * @param ppr                   Input : pointer to precision structure
 * @param pba                   Input : pointer to background structure
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param bin                   Input : redshift bin number
 * @param tau0_minus_tau        Output: values of (tau0-tau) at which source are sample
 * @param tau_size              Output: pointer to size of previous array
 * @param index_md            Input : index of mode
 * @param tau0                  Input : time today
 * @param interpolated_sources  Input : interpolated perturbation source
 * @param sources               Output: resampled transfer source
 * @return the error status
 */

int transfer_selection_sampling(
				struct precision * ppr,
				struct background * pba,
				struct perturbs * ppt,
				struct transfers * ptr,
				int bin,
				double * tau0_minus_tau,
				int tau_size) {
  
  /* running index on time */
  int index_tau;

  /* minimum and maximal value of time in new sampled interval */
  double tau_min,tau_mean,tau_max;

  /* time interval for this bin */
  class_call(transfer_selection_times(ppr,
				      pba,
				      ppt,
				      ptr,
				      bin,
				      &tau_min,
				      &tau_mean,
				      &tau_max),
	     ptr->error_message,
	     ptr->error_message);

  /* case selection == dirac */
  if (tau_min == tau_max) {
    class_test(tau_size !=1,
	       ptr->error_message,
	       "for Dirac selection function tau_size should be 1, not %d",tau_size); 
    tau0_minus_tau[0] = pba->conformal_age - tau_mean;
  }
  /* for other cases (gaussian, tophat...) define new sampled values
     of (tau0-tau) with even spacing */
  else {
    for (index_tau=0; index_tau<tau_size; index_tau++) {
      tau0_minus_tau[index_tau]=pba->conformal_age-tau_min-((double)index_tau)/((double)tau_size-1.)*(tau_max-tau_min);
    }
  }

  return _SUCCESS_;

}

/**
 * for lensing sources that need to be convolved with a selection
 * function, redefine the sampling within the range extending from the
 * tau_min of the selection function up to tau0
 *
 * @param ppr                   Input : pointer to precision structure
 * @param pba                   Input : pointer to background structure
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param bin                   Input : redshift bin number
 * @param tau0_minus_tau        Output: values of (tau0-tau) at which source are sample
 * @param tau_size              Output: pointer to size of previous array
 * @param index_md            Input : index of mode
 * @param tau0                  Input : time today
 * @param interpolated_sources  Input : interpolated perturbation source
 * @param sources               Output: resampled transfer source
 * @return the error status
 */

int transfer_lensing_sampling(
			      struct precision * ppr,
			      struct background * pba,
			      struct perturbs * ppt,
			      struct transfers * ptr,
			      int bin,
			      double tau0,
			      double * tau0_minus_tau,
			      int tau_size) {
  
  /* running index on time */
  int index_tau;

  /* minimum and maximal value of time in new sampled interval */
  double tau_min,tau_mean,tau_max;

  /* time interval for this bin */
  class_call(transfer_selection_times(ppr,
				      pba,
				      ppt,
				      ptr,
				      bin,
				      &tau_min,
				      &tau_mean,
				      &tau_max),
	     ptr->error_message,
	     ptr->error_message);

  for (index_tau=0; index_tau<tau_size; index_tau++) {
    tau0_minus_tau[index_tau]=pba->conformal_age-tau_min-((double)index_tau)/((double)tau_size-1.)*(tau0-tau_min);
  }

  return _SUCCESS_;

}


/**
 * for sources that need to be mutiplied by a selection function,
 * redefine a finer time sampling in a small range, and resample the
 * perturbation sources at the new value by linear interpolation
 *
 * @param ppr                   Input : pointer to precision structure
 * @param pba                   Input : pointer to background structure
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param bin                   Input : redshift bin number
 * @param tau0_minus_tau        Output: values of (tau0-tau) at which source are sample
 * @param tau_size              Output: pointer to size of previous array
 * @param index_md            Input : index of mode
 * @param tau0                  Input : time today
 * @param interpolated_sources  Input : interpolated perturbation source
 * @param sources               Output: resampled transfer source
 * @return the error status
 */

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
			     double * sources) {
  
  /* running index on time */
  int index_tau;

  /* running index on wavenumbers */
  int index_k_tr;

  /* array of values of source */
  double * source_at_tau;

  /* array of source values for a given time and for all k's */
  class_alloc(source_at_tau,
	      ptr->k_size[index_md]*sizeof(double),
	      ptr->error_message);

  /* interpolate the sources linearily at the new time values */
  for (index_tau=0; index_tau<tau_size; index_tau++) {

    class_call(array_interpolate_two(ppt->tau_sampling,
				     1,
				     0,
				     interpolated_sources, 
				     // this array is indexed as interpolated_sources[index_k_tr*ppt->tau_size+index_tau]
				     ptr->k_size[index_md],
				     ppt->tau_size,
				     tau0-tau0_minus_tau[index_tau],
				     source_at_tau,
				     ptr->k_size[index_md],
				     ptr->error_message),
	       ptr->error_message,
	       ptr->error_message);

    /* for each k, copy the new values in the output sources array */
    for (index_k_tr=0;index_k_tr<ptr->k_size[index_md];index_k_tr++) {
      sources[index_k_tr * tau_size + index_tau] = source_at_tau[index_k_tr];
    }
  } 

  /* deallocate the temporary array */
  free(source_at_tau);

  return _SUCCESS_;

}

/**
 * for each selection function, compute the min, mean and max values
 * of conformal time (associated to the min, mean and max values of
 * redshift specified by the user)
 *
 * @param ppr                   Input : pointer to precision structure
 * @param pba                   Input : pointer to background structure
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param bin                   Input : redshift bin number
 * @param tau_min               Output: smallest time in the selection interval
 * @param tau_mean              Output: time corresponding to z_mean
 * @param tau_max               Output: largest time in the selection interval
 * @return the error status
 */

int transfer_selection_times(
			     struct precision * ppr,
			     struct background * pba,
			     struct perturbs * ppt,
			     struct transfers * ptr,
			     int bin,
			     double * tau_min,
			     double * tau_mean,
			     double * tau_max) {

  /* a value of redshift */
  double z=0.;

  /* lower edge of time interval for this bin */
  
  if (ppt->selection==gaussian) {
    z = ppt->selection_mean[bin]+ppt->selection_width[bin]*ppr->selection_cut_at_sigma;
  }
  if (ppt->selection==tophat) {
    z = ppt->selection_mean[bin]+(1.+ppr->selection_cut_at_sigma*ppr->selection_tophat_edge)*ppt->selection_width[bin];
  }
  if (ppt->selection==dirac) {
    z = ppt->selection_mean[bin];
  }
  
  class_call(background_tau_of_z(pba,
				 z,
				 tau_min),
	     pba->error_message,
	     ppt->error_message);
  
  /* higher edge of time interval for this bin */
  
  if (ppt->selection==gaussian) {
    z = max(ppt->selection_mean[bin]-ppt->selection_width[bin]*ppr->selection_cut_at_sigma,0.);
  }
  if (ppt->selection==tophat) {
    z = max(ppt->selection_mean[bin]-(1.+ppr->selection_cut_at_sigma*ppr->selection_tophat_edge)*ppt->selection_width[bin],0.);
  }
  if (ppt->selection==dirac) {
    z = ppt->selection_mean[bin];
  }  

  class_call(background_tau_of_z(pba,
				 z,
				 tau_max),
	     pba->error_message,
	     ppt->error_message);

  /* central value of time interval for this bin */
  
  z = max(ppt->selection_mean[bin],0.);

  class_call(background_tau_of_z(pba,
				 z,
				 tau_mean),
	     pba->error_message,
	     ppt->error_message);

  return _SUCCESS_;

}

/**
 * compute and normalise selection function for a set of time values 
 *
 *
 * @param pba                   Input : pointer to background structure
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param selection             Output: normalized selection function
 * @param tau0_minus_tau        Input : values of (tau0-tau) at which source are sample
 * @param delta_tau             Input : corresponding values dtau, used for trapezoidal integration
 * @param tau_size              Input : size of previous two arrays
 * @param pvecback              Input : allocated array of background values
 * @param tau_0                 Input : time today
 * @param bin                   Input : redshift bin number
 * @return the error status
 */

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
			       int bin) {

  /* running index over time */
  int index_tau;

  /* running value of time */
  double tau;

  /* used for normalizing the selection to one */
  double norm;

  /* used for calling background_at_tau() */
  int last_index;

  /* runnign value of redshift */
  double z;

  /* loop over time */
  for (index_tau = 0; index_tau < tau_size; index_tau++) {
    
    /* running value of time */
    tau = tau0 - tau0_minus_tau[index_tau];
    
    /* get background quantitites at this time */
    class_call(background_at_tau(pba,
				 tau,
				 pba->long_info,
				 pba->inter_normal,
				 &last_index,
				 pvecback),
	       pba->error_message,
	       ptr->error_message);
    
    /* infer redhsift */
    z = pba->a_today/pvecback[pba->index_bg_a]-1.;

    /* get corresponding dN/dz(z,bin) */
    class_call(transfer_selection_function(ppr,
					   ppt,
					   ptr,
					   bin,
					   z,
					   &(selection[index_tau])),
	       ptr->error_message,
	       ptr->error_message);
    
    /* get corresponding dN/dtau = dN/dz * dz/dtau = dN/dz * H */
    selection[index_tau] *= pvecback[pba->index_bg_H];
    
  }

  /* compute norm = \int W(tau) dtau */
  norm = 0;
  for (index_tau = 0; index_tau <tau_size; index_tau++) {
    norm += selection[index_tau]*delta_tau[index_tau];
  }
  norm /= 2.; /* correct for factor 1/2 from trapezoidal rule */
  
  /* divide W by norm so that \int W(tau) dtau = 1 */
  for (index_tau = 0; index_tau < tau_size; index_tau++) {
    selection[index_tau]/=norm;
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
 * @param index_md            Input : index of mode
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

  /* whether to use the Limber approximation */
  short use_limber;
  
  /* whether to use a cutting scheme */
  enum transfer_cutting use_cut;

  if (index_l >= ptr->l_size_tt[index_md][index_tt]) { 

    for (index_k = 0; index_k < ptr->k_size[index_md]; index_k++) {
      ptr->transfer[index_md][((index_ic * ptr->tt_size[index_md] + index_tt)
				 * ptr->l_size[index_md] + index_l)
				* ptr->k_size[index_md] + index_k] = 0.;
    }
    return _SUCCESS_;
  }

  if (ptr->transfer_verbose > 2)
    printf("Compute transfer for l=%d\n",(int)l);

  /** use cutting scheme for CMB **/
  use_cut = tc_none;

  if (_scalars_ && ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t)))
    use_cut = ppr->transfer_cut;

  if (_scalars_ && ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_e)))
    use_cut = ppr->transfer_cut;

  if _tensors_
    use_cut = ppr->transfer_cut;

  /** - if the option of stopping the transfer function computation at some k_max is selected, initialize relevant quantities */
	      
  if (use_cut == tc_osc) {
    cut_transfer = _FALSE_;
    global_max=0.;
    global_min=0.;
    last_local_max=0.;
    last_local_min=0.;
    transfer_function=0.;
    transfer_last=0.;
    transfer_last_last=0;
  }
	      
  if (use_cut == tc_cl) {
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
  
  if (_scalars_ && ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_e))) {
    
    multiply_by_factor=_TRUE_;
    extra_factor=sqrt((l+2.) * (l+1.) * l * (l-1.));
  }
  
  /* for tensor temperature, multiply by 
     square root of (l+2)(l+1)l(l-1)/2 */
  
  if (_tensors_ && ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t))) {
    
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

  for (index_k = 0; index_k < ptr->k_size[index_md]; index_k++) {

    k = ptr->k[index_md][index_k];

    if (ptr->transfer_verbose > 3)
      printf("Compute transfer for l=%d k=%e type=%d\n",(int)l,k,index_tt);
		
    /* update previous transfer values in the tc_osc method */
    if (use_cut == tc_osc) {
      transfer_last_last = transfer_last;
      transfer_last = transfer_function;
    }
		
    /* update previous relative pseudo-C_l variation in the tc_cl method */
    if (use_cut == tc_cl) {
      cl_var_last_last = cl_var_last;
      cl_var_last = cl_var;
    }
		
    /* below k_max: compute transfer function */
    if ((use_cut == tc_none) || (cut_transfer == _FALSE_)) {

      class_call(transfer_use_limber(ppr,
				     ppt,
				     ptr,
				     k_max_bessel,
				     index_md,
				     index_tt,
				     k,
				     l,
				     &use_limber),
		 ptr->error_message,
		 ptr->error_message);

      if (use_limber == _TRUE_) {

	class_call(transfer_limber(tau_size,
				   ptr,
				   index_md,
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
	
	if ((use_cut == tc_env) &&
	    (((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t)) || 
	     ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_t)))) {
	  
	  class_call(transfer_envelop(tau_size,
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

	  class_call(transfer_integrate(ptr,
					tau_size,
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
    ptr->transfer[index_md][((index_ic * ptr->tt_size[index_md] + index_tt)
			       * ptr->l_size[index_md] + index_l)
			      * ptr->k_size[index_md] + index_k]
      = transfer_function;
    
    /* in the tc_osc case, update various quantities and 
       check wether k_max is reached */
    if ((use_cut == tc_osc) && (index_k>=2)) {
      
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
    if ((use_cut == tc_cl) && (index_k>=2) && (index_k<ptr->k_size[index_md]-1) && (transfer_function != 0.)) {
		  
      /* rough estimate of the contribution of the last step to pseudo-C_l, 
	 assuming flat primordial spectrum */
      delta_cl = transfer_function * transfer_function / k * 0.5 * (ptr->k[index_md][index_k+1] - ptr->k[index_md][index_k-1]);
      
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

int transfer_use_limber(
			struct precision * ppr,
			struct perturbs * ppt,
			struct transfers * ptr,
			double k_max_bessel,
			int index_md,
			int index_tt,
			double k,
			double l,
			short * use_limber) {


  /* criterium for chosing between integration and Limber 
     must be implemented here */
  
  *use_limber = _FALSE_;

  if (k>k_max_bessel) {
    *use_limber = _TRUE_;
  }
  else {
    
    if _scalars_ {
      
      if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (index_tt == ptr->index_tt_lcmb) && (l>ppr->l_switch_limber))
	*use_limber = _TRUE_;
      
      if ((ppt->has_cl_density == _TRUE_) && (index_tt >= ptr->index_tt_density) && (index_tt < ptr->index_tt_density+ppt->selection_num) && (l>=ppr->l_switch_limber_for_cl_density_over_z*ppt->selection_mean[index_tt-ptr->index_tt_density]))

	if (ppt->selection != dirac)
	  *use_limber = _TRUE_;

      if ((ppt->has_cl_lensing_potential == _TRUE_) && (index_tt >= ptr->index_tt_lensing) && (index_tt < ptr->index_tt_lensing+ppt->selection_num) && (l>=ppr->l_switch_limber_for_cl_density_over_z*ppt->selection_mean[index_tt-ptr->index_tt_lensing]))
	*use_limber = _TRUE_;
	
    }
  }

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
 * @param index_md            Input : index of mode
 * @param index_tt              Input : index of type
 * @param index_l               Input : index of multipole
 * @param index_k               Input : index of wavenumber
 * @param interpolated_sources  Input: array of interpolated sources
 * @param ptw                   Input : pointer to transfer_workspace structure (allocated in transfer_init() to avoid numerous reallocation) 
 * @param trsf                  Output: transfer function \f$ \Delta_l(k) \f$ 
 * @return the error status
 */

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

  //double Phi;

  /** - find minimum value of (tau0-tau) at which \f$ j_l(k[\tau_0-\tau]) \f$ is known, given that \f$ j_l(x) \f$ is sampled above some finite value \f$ x_{\min} \f$ (below which it can be approximated by zero) */  
  tau0_minus_tau_min_bessel = x_min_l/k; /* segmentation fault impossible, checked before that k != 0 */

  /** - if there is no overlap between the region in which bessels and sources are non-zero, return zero */
  if (tau0_minus_tau_min_bessel >= tau0_minus_tau[0]) {
    *trsf = 0.;
    return _SUCCESS_;
  }

  /** - if there is an overlap: */

  /** -> trivial case: the source is a Dirac function and is sampled in only one point */
  if (tau_size == 1) {

    x = k * tau0_minus_tau[0];

    index_x = (int)((x-x_min_l)/x_step);
  
    a = (x_min_l+x_step*(index_x+1) - x)/x_step;

    *trsf = sources[index_k]                         /* source */
      * (a * j_l[index_x]                               /* bessel */
	 + (1.-a) * ( j_l[index_x+1]
		      - a * ((a+1.) * ddj_l[index_x]
			     +(2.-a) * ddj_l[index_x+1]) 
		      * x_step * x_step / 6.0)); 

    return _SUCCESS_;  
    
  }

  /** -> other cases */

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
    
  /* Note on the trapezoidal method used here:
 
     Take a function y(x) sampled in n points [x_0, ..., x_{n-1}].
     The integral Sigma=int_{x_0}^{x_{n-1}} y(x) dx can be written as:
     
     sigma = sum_0^{n-2} [0.5 * (y_i+y_{i+1}) * (x_{i+1}-x_i)]
     = 0.5 * sum_0^{n-1} [y_i * delta_i]
     
     with delta_0     = x_1 - x_0
     delta_i     = x_{i+1} - x_{i-1}
     delta_{n-1} = x_{n-1} - x_{n-2}

     We will use the second expression (we have already defined
     delta_tau as the above delta_i).

     Suppose that we want to truncate the integral at some x_trunc <
     x_{n-1}, knowing that y(x_trunc)=0. We are in this case, with
     x-trunc corresponding to x_min of the bessel function. Let i_max
     be the last index such that y is non-zero: 

     x_{i_max} < x_trunc < x_{i_max+1}.

     We must adapt the formula above, with the last point being
     x_trunc, and knowing that y(x_trunc)=0:

     sigma = 0.5 * sum_0^{i_trunc-1} [y_i * delta_i]
     + 0.5 * y_{i_trunc} * (x_trunc - x_{i_max-1}) 
     + 0.      
	   
     Below we willuse exactly this expression, strating form the last term 
     [y_{i_trunc} * (x_trunc - x_{i_max-1})],
     then adding all the terms
     [y_i * delta_i],
     and finally multiplying by 0.5

     There is just one exception to the formula: the case when
     x_0<x_trunc<x_1, so that i_max=0. Then
      
     sigma = 0.5 * x_0 * (x_trunc-x_0)

     This exception is taken into account below.
   	   
  */

  /* Edge of the integral */
  
  x = k * tau0_minus_tau[index_tau_max];

  index_x = (int)((x-x_min_l)/x_step);
  
  a = (x_min_l+x_step*(index_x+1) - x)/x_step;

  /*  
      if (l==10. && index_k==60) {
      printf("%e %e %e %e %e %e %e %e %e %e\n",
      k,
      tau0_minus_tau[index_tau_max],
      sources[index_k * tau_size + index_tau_max],
      (a*j_l[index_x]+(1.-a) *( j_l[index_x+1]
      - a * ((a+1.) * ddj_l[index_x]
					 
      +(2.-a) * ddj_l[index_x+1]) 
      * x_step * x_step / 6.0)),
      (tau0_minus_tau[index_tau_max-1]-tau0_minus_tau_min_bessel),
      a,
      j_l[index_x],
      j_l[index_x+1],
      ddj_l[index_x],
      ddj_l[index_x+1]
      );
      }
  */

  transfer = sources[index_k * tau_size + index_tau_max] /* source */
    * (a * j_l[index_x] +                                /* bessel */
       (1.-a) * (j_l[index_x+1]
		 - a * ((a+1.) * ddj_l[index_x]
			+(2.-a) * ddj_l[index_x+1]) 
		 * x_step * x_step / 6.0) );

  /* (for bessel function cubic spline interpolation, we could have called the
     subroutine bessel_at_x; however we perform operations directly here
     in order to speed up the code) */

  // case with no truncation: normal edge:
  if (index_tau_max == tau_size-1) {
    transfer *=  delta_tau[index_tau_max];
  }
  // case with truncation at tau0_minus_tau_min_bessel:
  else {
    if (index_tau_max > 0)
      // case with several points in the integral
      transfer *= (tau0_minus_tau[index_tau_max-1]-tau0_minus_tau_min_bessel);
    else
      // case with only one non-zero point y(x_0) in the integral
      transfer *= (tau0_minus_tau[index_tau_max]-tau0_minus_tau_min_bessel);
  }

  /* rest of the integral */

  for (index_tau=0; index_tau<index_tau_max; index_tau++) {
    
    x = k * tau0_minus_tau[index_tau];
    
    index_x = (int)((x-x_min_l)/x_step);
    
    a = (x_min_l+x_step*(index_x+1) - x)/x_step;

    transfer += sources[index_k * tau_size + index_tau] /* source */
      * (a * j_l[index_x]                               /* bessel */
	 + (1.-a) * ( j_l[index_x+1]
		      - a * ((a+1.) * ddj_l[index_x]
			     +(2.-a) * ddj_l[index_x+1]) 
		      * x_step * x_step / 6.0)) 
      * delta_tau[index_tau];                           /* dtau */

    /* (for bessel function cubic spline interpolation, we could have called the
       subroutine bessel_at_x; however we perform operations directly here
       in order to speed up the code) */
    
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
* @param index_md            Input : index of mode
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
		    int index_md,
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
  int index_tau;

  /* interpolated source and its derivatives at this value */
  double S, dS, ddS;

  /** - get k, l and infer tau such that k(tau0-tau)=l+1/2; 
      check that tau is in appropriate range */

  tau0_minus_tau_limber = (l+0.5)/k;

  if ((tau0_minus_tau_limber > tau0_minus_tau[0]) || 
      (tau0_minus_tau_limber < tau0_minus_tau[tau_size-1])) {
    *trsf = 0.;
    return _SUCCESS_;
  }

  /** - find  bracketing indices. 
      index_tau must be at least 1 (so that index_tau-1 is at least 0) 
      and at most tau_size-2 (so that index_tau+1 is at most tau_size-1). 
  */
  index_tau=1;
  while ((tau0_minus_tau[index_tau] > tau0_minus_tau_limber) && (index_tau<tau_size-2))
    index_tau++;
  
  /** - interpolate by fitting a polynomial of order two; get source
      and its first two derivatives. Note that we are not
      interpolating S, but the product S*(tau0-tau). Indeed this
      product is regular in tau=tau0, while S alone diverges for
      lensing. */

  /* the case where the last of the three point is the edge (tau0=tau) must be treated separately, see below */
  if (index_tau < tau_size-2) { 

    class_call(array_interpolate_parabola(tau0_minus_tau[index_tau-1],
					  tau0_minus_tau[index_tau],
					  tau0_minus_tau[index_tau+1],
					  tau0_minus_tau_limber,
					  sources[index_k*tau_size+index_tau-1]*tau0_minus_tau[index_tau-1],
					  sources[index_k*tau_size+index_tau]*tau0_minus_tau[index_tau],
					  sources[index_k*tau_size+index_tau+1]*tau0_minus_tau[index_tau+1],
					  &S,
					  &dS,
					  &ddS,
					  ptr->error_message),
	       ptr->error_message,
	       ptr->error_message);

  }

  /* in this case, we have stored a zero for sources[index_k*tau_size+index_tau+1]. But we can use in very good approximation the fact that S*(tau0-tau) is constant near tau=tau0 and replace sources[index_k*tau_size+index_tau+1]*tau0_minus_tau[index_tau+1] by sources[index_k*tau_size+index_tau]*tau0_minus_tau[index_tau] */
  else {

    class_call(array_interpolate_parabola(tau0_minus_tau[index_tau-1],
					  tau0_minus_tau[index_tau],
					  tau0_minus_tau[index_tau+1],
					  tau0_minus_tau_limber,
					  sources[index_k*tau_size+index_tau-1]*tau0_minus_tau[index_tau-1],
					  sources[index_k*tau_size+index_tau]*tau0_minus_tau[index_tau],
					  sources[index_k*tau_size+index_tau]*tau0_minus_tau[index_tau],
					  &S,
					  &dS,
					  &ddS,
					  ptr->error_message),
	       ptr->error_message,
	       ptr->error_message);
  }

  /** - get transfer = source * sqrt(pi/(2l+1))/k 
      = source*[tau0-tau] * sqrt(pi/(2l+1))/(l+1/2) 
  */

  *trsf = sqrt(_PI_/(2.*l+1.))*S/(l+0.5);

  return _SUCCESS_;
  
}

/**
 * This routine computes the transfer functions \f$ \Delta_l^{X} (k)
 * \f$) for each mode, initial condition, type, multipole l and
 * wavenumber k, by using the Limber approximation at ordet two, i.e
 * as a function of the source function and its first two derivatives
 * at a single value of tau
 *
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param tau0                  Input : conformal time today
 * @param index_md            Input : index of mode
 * @param index_tt              Input : index of type
 * @param index_l               Input : index of multipole
 * @param index_k               Input : index of wavenumber
 * @param interpolated_sources  Input: array of interpolated sources
 * @param trsf                  Output: transfer function \f$ \Delta_l(k) \f$ 
 * @return the error status
 */

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
		     ){

  /** Summary: */

  /** - define local variables */

  /* conformal time at which source must be computed */
  double tau0_minus_tau_limber;
  int index_tau;

  /* interpolated source and its derivatives */
  double S, dS, ddS;

  /** - get k, l and infer tau such that k(tau0-tau)=l+1/2; 
      check that tau is in appropriate range */

  tau0_minus_tau_limber = (l+0.5)/k;

  if ((tau0_minus_tau_limber > tau0_minus_tau[0]) || 
      (tau0_minus_tau_limber < tau0_minus_tau[tau_size-1])) {
    *trsf = 0.;
    return _SUCCESS_;
  }

  /** - find  bracketing indices */
  index_tau=0;
  while ((tau0_minus_tau[index_tau] > tau0_minus_tau_limber) && (index_tau<tau_size-2))
    index_tau++;

  /** - interpolate by fitting a polynomial of order two; get source
      and its first two derivatives */
  class_call(array_interpolate_parabola(tau0_minus_tau[index_tau-1],
					tau0_minus_tau[index_tau],
					tau0_minus_tau[index_tau+1],
					tau0_minus_tau_limber,
					sources[index_k*tau_size+index_tau-1],
					sources[index_k*tau_size+index_tau],
					sources[index_k*tau_size+index_tau+1],
					&S,
					&dS,
					&ddS,
					ptr->error_message),
	     ptr->error_message,
	     ptr->error_message);
	     

  /** - get transfer from 2nd order Limber approx (infered from 0809.5112 [astro-ph]) */

  *trsf = sqrt(_PI_/(2.*l+1.))/k*((1.-3./2./(2.*l+1.)/(2.*l+1.))*S+dS/k/(2.*l+1.)-0.5*ddS/k/k);

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
 * @param index_md            Input : index of mode
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

