/** @file perturbations.c Documented perturbation module
 *
 * Julien Lesgourgues, 23.09.2010    
 *
 * Deals with the perturbation evolution.
 * This module has two purposes: 
 *
 * - at the beginning, to initialize the perturbations, i.e. to
 * integrate the perturbation equations, and store temporarily the terms
 * contributing to the source functions as a function of conformal
 * time. Then, to perform a few manipulations of these terms in order to
 * infer the actual source functions \f$ S^{X} (k, \eta) \f$, and to
 * store them as a function of conformal time inside an interpolation
 * table.
 *
 * - at any time in the code, to evaluate the source functions at a
 * given conformal time (by interpolating within the interpolation
 * table).
 *
 * Hence the following functions can be called from other modules:
 *
 * -# perturb_init() at the beginning (but after background_init() and thermodynamics_init())  
 * -# perturb_sources_at_eta() at any later time
 * -# perturb_free() at the end, when no more calls to perturb_sources_at_eta() are needed
 */

#include "perturbations.h"

/** 
 * Source function \f$ S^{X} (k, \eta) \f$ at a given conformal time eta.
 *
 * Evaluate all source functions at given conformal time eta by reading
 * the pre-computed table and interpolating.
 *
 * @param ppt        Input : pointer to perturbation structure containing interpolation tables
 * @param index_mode Input : index of requested mode
 * @param index_ic   Input : index of requested initial condition
 * @param index_k    Input : index of requested wavenumber
 * @param index_type Input : index of requested source function type
 * @param eta        Input : any value of conformal time
 * @param psource    Output: vector (assumed to be already allocated) of source functions
 * @return the error status
 */

int perturb_sources_at_eta(
			   struct perturbs * ppt,
			   int index_mode,
			   int index_ic,
			   int index_k,
			   int index_type,
			   double eta,
			   double * psource
			   ) {

  /** Summary: */

  /** - interpolate in pre-computed table contained in ppt */
  class_call(array_interpolate_two(&(ppt->sources[index_mode]
				     [index_ic * ppt->tp_size[index_mode] + index_type]
				     [index_k * ppt->eta_size]),
				   1,
				   0,
				   ppt->eta_sampling,
				   1,
				   ppt->eta_size,
				   eta,
				   psource,
				   1,
				   ppt->error_message),
	     ppt->error_message,
	     ppt->error_message);

  return _SUCCESS_;
}

/** 
 * Initialize the perturbs structure, and in particular the table of source functions.
 * 
 * Main steps:
 *
 * - given the values of the flags describing which kind of
 *    perturbations should be considered (modes: scalar/vector/tensor,
 *    initial conditions, type of source functions needed...),
 *    initialize indices and wavenumber list
 *
 * - define the time sampling for the output source functions
 *
 * - for each mode (scalar/vector/tensor): initialize the indices of
 *    relevant perturbations, integrate the differential system,
 *    compute and store the source functions.
 *
 * @param ppr Input : pointer to precision structure
 * @param pba Input : pointer to background strucutre
 * @param pth Input : pointer to thermodynamics structure
 * @param ppt Output: Initialized perturbation structure
 * @return the error status
 */

int perturb_init(
		 struct precision * ppr,
		 struct background * pba,
		 struct thermo * pth,
		 struct perturbs * ppt
		 ) {
  
  /** Summary: */

  /** - define local variables */

  /* running index for modes */
  int index_mode; 
  /* running index for initial conditions */
  int index_ic; 
  /* running index for wavenumbers */
  int index_k; 
  /* struct perturb_workspace pw; */
  struct perturb_workspace * ppw;

  /* This code can be optionally compiled with the openmp option for parallel computation.
     Inside parallel regions, the use of the command "return" is forbidden.
     For error management, instead of "return _FAILURE_", we will set the variable below
     to "abort = _TRUE_". This will lead to a "return _FAILURE_" jus after leaving the 
     parallel region. */
  int abort;

#ifdef _OPENMP
  /* number of available omp threads */
  int number_of_threads;
  /* instrumentation times */
  double tstart, tstop;
  /* pointer to one "ppw" per thread */
  struct perturb_workspace ** pppw;
#endif

  /** - preliminary checks */

  if (ppt->has_perturbations == _FALSE_) {
    if (ppt->perturbations_verbose > 0)
      printf("No sources requested. Perturbation module skipped.\n");
    return _SUCCESS_;
  }
  else {
    if (ppt->perturbations_verbose > 0)
      printf("Computing sources\n");
  }

  class_test(ppt->has_vectors == _TRUE_,
	     ppt->error_message,
	     "Vectors not coded yet");

  if ((ppt->has_bi == _TRUE_) || (ppt->has_cdi == _TRUE_) || (ppt->has_nid == _TRUE_) || (ppt->has_niv == _TRUE_)) {
    printf("Warning: so far, isocurvature initial condition only implemented at first order in (k eta), not very precise...\n");
  }

  if ((ppt->has_cl_cmb_temperature == _TRUE_) && (ppt->has_cl_cmb_polarization == _TRUE_) &&
      (ppt->has_tensors == _TRUE_)) {
      printf("Warning: our C_l^TE for tensors has a minus sign with respect to CAMB 2008. Mistake in one of the two codes? To be checked.\n");
  }

  if (ppt->has_cl_cmb_lensing_potential == _TRUE_) {
      printf("Warning: so far, C_l^phiphi and C_l^Tphi computed using Limber for all l's, not very precise...\n");
  }

  class_test(ppr->eta_min_over_sampling_min >= 1,
	     ppt->error_message,
	     "you have ppr->eta_min_over_sampling_min=%g, should always be < 1 to guarantee that integration starts early enough and that the first value of the source can be trusted"); 

  /** - initialize all indices and lists in perturbs structure using perturb_indices_of_perturbs() */

  class_call(perturb_indices_of_perturbs(ppr,
					 pba,
					 pth,
					 ppt),
	     ppt->error_message,
	     ppt->error_message);

  /** - define the common time sampling for all sources using perturb_timesampling_for_sources() */

  class_call(perturb_timesampling_for_sources(ppr,
					      pba,
					      pth,
					      ppt),
	     ppt->error_message,
	     ppt->error_message);


  /** - create an array of workspaces in multi-thread case */

#ifdef _OPENMP

#pragma omp parallel 
  {
    number_of_threads = omp_get_num_threads();
  }
  class_alloc(pppw,number_of_threads * sizeof(struct perturb_workspace *),ppt->error_message);
  
#endif

  /** - loop over modes (scalar, tensors, etc). For each mode: */

  for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {

    abort = _FALSE_;

#pragma omp parallel				\
  shared(pppw,ppr,pba,pth,ppt,index_mode,abort)	\
  private(ppw)

    {

      /** create a workspace (one per thread in multi-thread case) */

      class_alloc_parallel(ppw,sizeof(struct perturb_workspace),ppt->error_message);


#ifdef _OPENMP
      pppw[omp_get_thread_num()]=ppw;
#endif


      /** (a) initialize indices of vectors of perturbations with perturb_indices_of_current_vectors() */

      class_call_parallel(perturb_workspace_init(ppr,
						 pba,
						 pth,
						 ppt,
						 index_mode,
					 	 ppw),
			  ppt->error_message,
			  ppt->error_message);

    } /* end of parallel region */

    if (abort == _TRUE_) return _FAILURE_;

    /** (c) loop over initial conditions and wavenumbers; for each of them, evolve perturbations and compute source functions with perturb_solve() */

    for (index_ic = 0; index_ic < ppt->ic_size[index_mode]; index_ic++) {

      abort = _FALSE_;

#pragma omp parallel						\
  shared(pppw,ppr,pba,pth,ppt,index_mode,index_ic,abort)	\
  private(index_k,ppw,tstart,tstop)

      {

#ifdef _OPENMP
	ppw=pppw[omp_get_thread_num()];
	tstart = omp_get_wtime();
#endif
	
#pragma omp for schedule (dynamic)

	for (index_k = 0; index_k < ppt->k_size[index_mode]; index_k++) {
	  
	  if ((ppt->perturbations_verbose > 2) && (abort == _FALSE_))
	    printf("evolving mode k=%e /Mpc\n",(ppt->k[index_mode])[index_k]);
	  
	  class_call_parallel(perturb_solve(ppr,
					    pba,
					    pth,
					    ppt,
					    index_mode,
					    index_ic,
					    index_k,
					    ppw),
			      ppt->error_message,
			      ppt->error_message);

#pragma omp flush(abort)

	} /* end of loop over wavenumbers */

#ifdef _OPENMP
	tstop = omp_get_wtime();
	if (ppt->perturbations_verbose>1)
	  printf("In %s: time spent in parallel region (loop over k's) = %e s for thread %d\n",
		 __func__,tstop-tstart,omp_get_thread_num());
#endif

      } /* end of parallel region */

      if (abort == _TRUE_) return _FAILURE_;

    } /* end of loop over initial conditions */
    
    abort = _FALSE_;

#pragma omp parallel				\
  shared(pppw,ppt,index_mode,abort)		\
  private(ppw)

    {

#ifdef _OPENMP
      ppw=pppw[omp_get_thread_num()];
#endif
      
      class_call_parallel(perturb_workspace_free(ppt,index_mode,ppw),
			  ppt->error_message,
			  ppt->error_message);

    } /* end of parallel region */

    if (abort == _TRUE_) return _FAILURE_;

  } /* end loop over modes */    

#ifdef _OPENMP
  free(pppw);
#endif

  return _SUCCESS_;
}

/**
 * Free all memory space allocated by perturb_init().
 * 
 * To be called at the end of each run, only when no further calls to
 * perturb_sources_at_eta() are needed.
 *
 * @param ppt Input: perturbation structure to be freed
 * @return the error status
 */

int perturb_free(
		 struct perturbs * ppt
		 ) {

  int index_mode,index_ic,index_k,index_type;

  if (ppt->has_perturbations == _TRUE_) {

    for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {
    
      for (index_ic = 0; index_ic < ppt->ic_size[index_mode]; index_ic++) {

	for (index_type = 0; index_type < ppt->tp_size[index_mode]; index_type++) {

	  free(ppt->sources[index_mode][index_ic*ppt->tp_size[index_mode]+index_type]);

	}

      }

      free(ppt->k[index_mode]);

      free(ppt->sources[index_mode]);

    }
    
    free(ppt->eta_sampling);
	 
    free(ppt->tp_size);

    free(ppt->ic_size);

    free(ppt->k_size);

    free(ppt->k_size_cl);

    free(ppt->k);

    free(ppt->sources);

  }

  return _SUCCESS_;

}

/** 
 * Initialize all indices and allocate most arrays in perturbs structure.
 *
 * @param ppr Input : pointer to precision structure
 * @param pba Input : pointer to background strucutre
 * @param pth Input : pointer to thermodynamics structure
 * @param ppt Input/Output: Initialized perturbation structure
 * @return the error status
 */

int perturb_indices_of_perturbs(
				struct precision * ppr,
				struct background * pba,
				struct thermo * pth,
				struct perturbs * ppt
				) {

  /** Summary: */

  /** - define local variables */

  int index_type, index_mode, index_ic, index_type_common;
  int k_list_size,k_list_cl_size;

  /** - count modes (scalar, vector, tensor) and assign corresponding indices */

  index_mode = 0;

  if (ppt->has_scalars == _TRUE_) {
    ppt->index_md_scalars = index_mode;
    index_mode++;      
  }
  if (ppt->has_vectors == _TRUE_) {
    ppt->index_md_vectors = index_mode;
    index_mode++;      
  }
  if (ppt->has_tensors == _TRUE_) {
    ppt->index_md_tensors = index_mode;
    index_mode++;      
  }

  ppt->md_size = index_mode;

  class_test(index_mode == 0,
	     ppt->error_message,
	     "you should have at least one out of {scalars, vectors, tensors} !!!");

  /** - allocate array of number of types for each mode, ppt->tp_size[index_mode] */

  class_alloc(ppt->tp_size,ppt->md_size*sizeof(int),ppt->error_message);

  /** - allocate array of number of initial conditions for each mode, ppt->ic_size[index_mode] */

  class_alloc(ppt->ic_size,ppt->md_size*sizeof(int),ppt->error_message);

  /** - allocate array of number of wavenumbers for each mode, ppt->k_size[index_mode] */

  class_alloc(ppt->k_size,ppt->md_size * sizeof(int),ppt->error_message);
  class_alloc(ppt->k_size_cl,ppt->md_size * sizeof(int),ppt->error_message);

  /** - allocate array of lists of wavenumbers for each mode, ppt->k[index_mode] */

  class_alloc(ppt->k,ppt->md_size * sizeof(double *),ppt->error_message);

  /** - allocate array of arrays of source functions for each mode, ppt->source[index_mode] */

  class_alloc(ppt->sources,ppt->md_size * sizeof(double *),ppt->error_message);

  /** - count source types that all modes have in common (temperature, polarization, ...), and assign corresponding indices */

  index_type = 0;
  ppt->has_cmb = _FALSE_; /* initialization (will eventually be set to true later) */
  ppt->has_lss = _FALSE_; /* initialization (will eventually be set to true later) */

  if (ppt->has_cl_cmb_temperature == _TRUE_) {
    ppt->has_source_t = _TRUE_;
    ppt->has_cmb = _TRUE_;
    ppt->index_tp_t = index_type; 
    index_type++;
  }
  else 
    ppt->has_source_t = _FALSE_;

  if (ppt->has_cl_cmb_polarization == _TRUE_) {
    ppt->has_source_e = _TRUE_;
    ppt->has_cmb = _TRUE_;
    ppt->index_tp_e = index_type; 
    index_type++;
  }
  else
    ppt->has_source_e = _FALSE_;

  index_type_common = index_type;

  /** - loop over modes. Initialize flags and indices which are specific to each mode. */

  for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {

    index_ic = 0;
    index_type = index_type_common;

    /** (a) scalars */

    if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

      /** -- count source types specific to scalars (gravitational potential, ...) and assign corresponding indices */

      if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) || (ppt->has_pk_matter == _TRUE_)) { 
	ppt->has_lss = _TRUE_;
        ppt->has_source_g = _TRUE_;
	ppt->index_tp_g = index_type; 
	index_type++;
      }
      else
	ppt->has_source_g = _FALSE_;

      ppt->tp_size[index_mode] = index_type;

      class_test(index_type == 0,
		 ppt->error_message,
		 "inconsistent input: you asked for scalars, so you should have at least one non-zero scalar source type (temperature, polarisation, lensing/gravitational potential, ...). Please adjust your input.");

    /** -- count scalar initial conditions (for scalars: ad, cdi, nid, niv; for tensors: only one) and assign corresponding indices */

      if (ppt->has_ad == _TRUE_) {
	ppt->index_ic_ad = index_ic;
	index_ic++;
      }
      if (ppt->has_bi == _TRUE_) {
	ppt->index_ic_bi = index_ic;
	index_ic++;
      }
      if (ppt->has_cdi == _TRUE_) {
	ppt->index_ic_cdi = index_ic;
	index_ic++;
      }
      if (ppt->has_nid == _TRUE_) {
	ppt->index_ic_nid = index_ic;
	index_ic++;
      }
      if (ppt->has_niv == _TRUE_) {
	ppt->index_ic_niv = index_ic;
	index_ic++;
      }
      ppt->ic_size[index_mode] = index_ic;

      class_test(index_ic == 0,
		 ppt->error_message,
		 "you should have at least one adiabatic or isocurvature initial condition...} !!!");

    }
    
    if ((ppt->has_vectors == _TRUE_) && (index_mode == ppt->index_md_vectors)) {

      /* vectors not treated yet */
    }

    /** (b) tensors */

    if ((ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors)) {

      /** -- count source types specific to tensors (B-polarization, ...) and assign corresponding indices */

      if (ppt->has_cl_cmb_polarization == _TRUE_) {
	ppt->has_source_b = _TRUE_;
	ppt->has_cmb = _TRUE_;
	ppt->index_tp_b = index_type; 
	index_type++;
      }
      else
	ppt->has_source_b = _FALSE_;

      ppt->tp_size[index_mode] = index_type;

      class_test(index_type == 0,
		 ppt->error_message,
		 "inconsistent input: you asked for tensors, so you should have at least one non-zero tensor source type (temperature or polarisation). Please adjust your input.");

      /** -- only one initial condition for tensors*/

      ppt->index_ic_ten = index_ic;
      index_ic++;
      ppt->ic_size[index_mode] = index_ic;

    }

    /** (c) for each mode, count values of k with perturb_get_k_list_size() */

    class_call(perturb_get_k_list(ppr,
				  pba,
				  pth,
				  ppt,
				  index_mode),
	       ppt->error_message,
	       ppt->error_message);

    /** (d) for each mode, allocate array of arrays of source functions for each initial conditions and wavenumber, (ppt->source[index_mode])[index_ic][index_type] */
    
    class_alloc(ppt->sources[index_mode],
		ppt->ic_size[index_mode] * ppt->tp_size[index_mode] * sizeof(double *),
		ppt->error_message);

  }
  
  return _SUCCESS_;

}

/**
 * Define time sampling for source functions. 
 *
 * For each type, compute the list of values of eta at which sources
 * will be sampled.  Knowing the number of eta values, allocate all
 * arrays of source functions.
 *
 * @param ppr Input : pointer to precision structure
 * @param pba Input : pointer to background strucutre
 * @param pth Input : pointer to thermodynamics structure
 * @param ppt Input/Output: Initialized perturbation structure
 * @return the error status
 */ 

int perturb_timesampling_for_sources(
				     struct precision * ppr,
				     struct background * pba,
				     struct thermo * pth,
				     struct perturbs * ppt
				     ) {

  /** Summary: */

  /** - define local variables */

  int counter;
  int index_mode;
  int index_type;
  int index_ic;
  int last_index_back;
  int last_index_thermo;

  double eta;
  double timescale_source;
  double rate_thermo;
  double rate_isw_squared;
  double a_prime_over_a;
  double a_primeprime_over_a;
  double eta_visibility_start_sources;
  double * pvecback;
  double * pvecthermo;

  /** - allocate background/thermodynamics vectors */

  class_alloc(pvecback,pba->bg_size_short*sizeof(double),ppt->error_message);  
  class_alloc(pvecthermo,pth->th_size*sizeof(double),ppt->error_message);

  /** - compute conformal time corresponding to opaque universe
      (starting point for source sampling) */

  class_call(background_eta_of_z(pba,pth->z_visibility_start_sources,&eta_visibility_start_sources),
	     pba->error_message,
	     ppt->error_message);

  /** - first, just count the number of sampling points in order to allocate the array containing all values: */

  /** (a) if CMB requested, first sampling point = when the universe stops being opaque; otherwise,
      start sampling gravitational potential at recombination */
  if (ppt->has_cmb == _TRUE_)
    eta = eta_visibility_start_sources;
  else 
    eta = pth->eta_rec;
  
  counter = 1;

  /** (b) next sampling point = previous + ppr->perturb_sampling_stepsize * timescale_source, where:
      - if CMB requested:
      timescale_source1 = \f$ |g/\dot{g}| = |\dot{\kappa}-\ddot{\kappa}/\dot{\kappa}|^{-1} \f$;
      timescale_source2 = \f$ |2\ddot{a}/a-(\dot{a}/a)^2|^{-1/2} \f$ (to sample correctly the late ISW effect; and 
      timescale_source=1/(1/timescale_source1+1/timescale_source2); repeat till today.
      - if CMB not requested:
      timescale_source = 1/aH; repeat till today.
  */

  while (eta < pba->conformal_age) {

    class_call(background_at_eta(pba,
				 eta, 
				 short_info, 
				 normal, 
				 &last_index_back, 
				 pvecback),
	       pba->error_message,
	       ppt->error_message);

    class_call(thermodynamics_at_z(pba,
				   pth,
				   1./pvecback[pba->index_bg_a]-1.,  /* redshift z=1/a-1 */
				   normal,
				   &last_index_thermo,
				   pvecback,
				   pvecthermo),
	       pth->error_message,
	       ppt->error_message);

    if (ppt->has_cmb == _TRUE_) {

      /* variation rate of thermodynamics variables */
      rate_thermo = pvecthermo[pth->index_th_rate];
      
      /* variation rate of metric due to late ISW effect (important at late times) */
      a_prime_over_a = pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a];
      a_primeprime_over_a = pvecback[pba->index_bg_H_prime] * pvecback[pba->index_bg_a]
	+ 2. * a_prime_over_a * a_prime_over_a;
      rate_isw_squared = fabs(2.*a_primeprime_over_a-a_prime_over_a*a_prime_over_a);
	
      /* compute rate */
      timescale_source = sqrt(rate_thermo*rate_thermo+rate_isw_squared);
    }
    else {
      /* variation rate given by Hubble time */
      a_prime_over_a = pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a];

      timescale_source = a_prime_over_a;
    }

    /* check it is non-zero */
    class_test(timescale_source == 0.,
	       ppt->error_message,
	       "null evolution rate, integration is diverging");

    /* compute inverse rate */
    timescale_source = 1./timescale_source;

    class_test(ppr->perturb_sampling_stepsize*timescale_source < ppr->smallest_allowed_variation,
	       ppt->error_message,
	       "integration step =%e < machine precision : leads either to numerical error or infinite loop",ppr->perturb_sampling_stepsize*timescale_source);

    eta = eta + ppr->perturb_sampling_stepsize*timescale_source; 
    counter++;

  }

  /** - infer total number of time steps, ppt->eta_size */
  ppt->eta_size = counter;

  /** - allocate array of time steps, ppt->eta_sampling[index_eta] */
  class_alloc(ppt->eta_sampling,ppt->eta_size * sizeof(double),ppt->error_message);

  /** - repeat the same steps, now filling the array with each eta value: */

  /** (a) first sampling point = when the universe stops being opaque */
  if (ppt->has_cmb == _TRUE_)
    eta = eta_visibility_start_sources;
  else
    eta = pth->eta_rec;

  counter = 0;
  ppt->eta_sampling[counter]=eta;

  /** (b) next sampling point = previous + ppr->perturb_sampling_stepsize * timescale_source, where
      timescale_source1 = \f$ |g/\dot{g}| = |\dot{\kappa}-\ddot{\kappa}/\dot{\kappa}|^{-1} \f$;
      timescale_source2 = \f$ |2\ddot{a}/a-(\dot{a}/a)^2|^{-1/2} \f$ (to smaple correctly the late ISW effect; and 
      timescale_source=1/(1/timescale_source1+1/timescale_source2); repeat till today
      - if CMB not requested:
      timescale_source = 1/aH; repeat till today.  */

  while (eta < pba->conformal_age) {
    
    class_call(background_at_eta(pba,
				 eta, 
				 short_info, 
				 normal, 
				 &last_index_back, 
				 pvecback),
	       pba->error_message,
	       ppt->error_message);

    class_call(thermodynamics_at_z(pba,
				   pth,
				   1./pvecback[pba->index_bg_a]-1.,  /* redshift z=1/a-1 */
				   normal,
				   &last_index_thermo,
				   pvecback,
				   pvecthermo),
	       pth->error_message,
	       ppt->error_message);

    if (ppt->has_cmb == _TRUE_) {

      /* variation rate of thermodynamics variables */
      rate_thermo = pvecthermo[pth->index_th_rate];

      /* variation rate of metric due to late ISW effect (important at late times) */
      a_prime_over_a = pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a];
      a_primeprime_over_a = pvecback[pba->index_bg_H_prime] * pvecback[pba->index_bg_a]
	+ 2. * a_prime_over_a * a_prime_over_a;
      rate_isw_squared = fabs(2.*a_primeprime_over_a-a_prime_over_a*a_prime_over_a);

      /* compute rate */
      timescale_source = sqrt(rate_thermo*rate_thermo+rate_isw_squared);
    }
    else {
      a_prime_over_a = pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a];
      timescale_source = a_prime_over_a;
    }
    /* check it is non-zero */
    class_test(timescale_source == 0.,
	       ppt->error_message,
	       "null evolution rate, integration is diverging");

    /* compute inverse rate */
    timescale_source = 1./timescale_source;

    class_test(ppr->perturb_sampling_stepsize*timescale_source < ppr->smallest_allowed_variation,
	       ppt->error_message,
	       "integration step =%e < machine precision : leads either to numerical error or infinite loop",ppr->perturb_sampling_stepsize*timescale_source);

    eta = eta + ppr->perturb_sampling_stepsize*timescale_source; 
    counter++;
    ppt->eta_sampling[counter]=eta;

  }

  /** - last sampling point = exactly today */
  ppt->eta_sampling[counter] = pba->conformal_age;

  free(pvecback);
  free(pvecthermo);

  /** - loop over modes, initial conditions and types. For each of
        them, allocate array of source functions. */
  
  for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {
    for (index_ic = 0; index_ic < ppt->ic_size[index_mode]; index_ic++) {
      for (index_type = 0; index_type < ppt->tp_size[index_mode]; index_type++) {

	class_alloc(ppt->sources[index_mode][index_ic*ppt->tp_size[index_mode]+index_type],
		    ppt->k_size[index_mode] * ppt->eta_size * sizeof(double),
		    ppt->error_message);

      }
    }
  }

  return _SUCCESS_;
}

/**
 * Define the number of comoving wavenumbers using the information
 * passed in the precision structure.
 *
 * @param ppr        Input : pointer to precision structure
 * @param pba        Input : pointer to background strucutre
 * @param pth        Input : pointer to thermodynamics structure
 * @param ppt        Input : pointer to perturbation structure
 * @param index_mode Input: index describing the mode (scalar, tensor, etc.)
 * @return the error status
 */

int perturb_get_k_list(
		       struct precision * ppr,
		       struct background * pba,
		       struct thermo * pth,
		       struct perturbs * ppt,
		       int index_mode
		       ) {
  int index_k;
  double k,k_next,k_rec,step,k_max_cl;

  /** Summary: */

  /** - get number of wavenumbers for scalar mode */

  if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

    class_test(ppr->k_scalar_step_transition == 0.,
	       ppt->error_message,
	       "stop to avoid division by zero");
    
    class_test(pth->rs_rec == 0.,
	       ppt->error_message,
	       "stop to avoid division by zero");

    k_rec = 2. * _PI_ / pth->rs_rec; /* comoving scale corresping to sound horizon at recombination */

    index_k=0;
    k = ppr->k_scalar_min * pba->H0;
    index_k=1;

    if (ppt->has_cls == _TRUE_)
      k_max_cl = ppt->l_scalar_max/ppr->l_max_over_k_max_scalars;
    else
      k_max_cl = 0.;

    while (k < k_max_cl) {
      step = ppr->k_scalar_step_super 
	+ 0.5 * (tanh((k-k_rec)/k_rec/ppr->k_scalar_step_transition)+1.) * (ppr->k_scalar_step_sub-ppr->k_scalar_step_super);

      class_test(step * k_rec < ppr->smallest_allowed_variation,
		 ppt->error_message,
		 "k step =%e < machine precision : leads either to numerical error or infinite loop",step * k_rec);

      k_next=k + step * k_rec;
      index_k++;
      k=k_next;
    }
    ppt->k_size_cl[index_mode] = index_k;

    if (ppt->has_pk_matter == _TRUE_) {
      if (k < ppt->k_scalar_kmax_for_pk*pba->h) {

	index_k += (int)((log(ppt->k_scalar_kmax_for_pk*pba->h/k)/log(10.))*ppr->k_scalar_k_per_decade_for_pk)+1;
	
      }
    }
    ppt->k_size[index_mode] = index_k;

    class_alloc(ppt->k[index_mode],ppt->k_size[index_mode]*sizeof(double),ppt->error_message);

    /** - repeat the same steps, now filling the array */

    index_k=0;
    ppt->k[index_mode][index_k] = ppr->k_scalar_min * pba->H0;
    index_k++;
    while (index_k < ppt->k_size_cl[index_mode]) {
      step = ppr->k_scalar_step_super 
	+ 0.5 * (tanh((ppt->k[index_mode][index_k-1]-k_rec)/k_rec/ppr->k_scalar_step_transition)+1.) * (ppr->k_scalar_step_sub-ppr->k_scalar_step_super);

      class_test(step * k_rec < ppr->smallest_allowed_variation,
		 ppt->error_message,
		 "k step =%e < machine precision : leads either to numerical error or infinite loop",step * k_rec);
      ppt->k[index_mode][index_k]=ppt->k[index_mode][index_k-1] + step * k_rec;
      index_k++;
    }

    while (index_k < ppt->k_size[index_mode]) {
      
      ppt->k[index_mode][index_k] = ppt->k[index_mode][index_k-1] 
	* exp(log(ppt->k_scalar_kmax_for_pk*pba->h/ppt->k[index_mode][ppt->k_size_cl[index_mode]-1])/(ppt->k_size[index_mode]-ppt->k_size_cl[index_mode]));
      index_k++;

    }

  }

  /** - get number of wavenumbers for tensor mode */

  if ((ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors)) {

    class_test(ppr->k_tensor_step_transition == 0.,
	       ppt->error_message,
	       "stop to avoid division by zero");

    class_test(pth->rs_rec == 0.,
	       ppt->error_message,
	       "stop to avoid division by zero");

    k_rec = 2. * _PI_ / pth->eta_rec; /* comoving scale corresping to causal horizon at recombination 
					 (roughly, sqrt(3) bigger than sound horizon) */

    index_k=0;
    k = ppr->k_tensor_min * pba->H0;
    index_k=1;

    while (k < ppt->l_tensor_max/ppr->l_max_over_k_max_tensors) {
      step = ppr->k_tensor_step_super 
	+ 0.5 * (tanh((k-k_rec)/k_rec/ppr->k_tensor_step_transition)+1.) * (ppr->k_tensor_step_sub-ppr->k_tensor_step_super);

      class_test(step * k_rec < ppr->smallest_allowed_variation,
		 ppt->error_message,
		 "k step =%e < machine precision : leads either to numerical error or infinite loop",step * k_rec);

      k_next=k + step * k_rec;
      index_k++;
      k=k_next;
    }

    ppt->k_size_cl[index_mode] = index_k;
    ppt->k_size[index_mode] = index_k;

    class_alloc(ppt->k[index_mode],ppt->k_size[index_mode]*sizeof(double),ppt->error_message);

    /** - repeat the same steps, now filling the array */

    index_k=0;
    ppt->k[index_mode][index_k] = ppr->k_tensor_min * pba->H0;
    index_k++;
    while (index_k < ppt->k_size_cl[index_mode]) {
      step = ppr->k_tensor_step_super 
	+ 0.5 * (tanh((ppt->k[index_mode][index_k-1]-k_rec)/k_rec/ppr->k_tensor_step_transition)+1.) * (ppr->k_tensor_step_sub-ppr->k_tensor_step_super);

      class_test(step * k_rec < ppr->smallest_allowed_variation,
		 ppt->error_message,
		 "k step =%e < machine precision : leads either to numerical error or infinite loop",step * k_rec);
      ppt->k[index_mode][index_k]=ppt->k[index_mode][index_k-1] + step * k_rec;
      index_k++;
    }

  }

  /* vectors not coded yet */

  return _SUCCESS_;

}

/**
 * Initialize a perturb_workspace structure. All fields are allocated
 * here, with the exception of the perturb_vector '->pv' field, which
 * is allocated separately in perturb_vector_init. We allocate one
 * such perturb_workspace structure per thread and per mode
 * (scalar/../tensor). Then, for each thread, all initial conditions
 * and wavenumbers will use the same workspace.
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param ppt        Input: pointer to the perturbation structure
 * @param index_mode Input: index of mode under consideration (scalar/.../tensor)
 * @param ppw        Input/Output: pointer to perturb_workspace structure which fields are allocated or filled here
 * @return the error status
 */

int perturb_workspace_init(
			   struct precision * ppr,
			   struct background * pba,
			   struct thermo * pth,
			   struct perturbs * ppt,
			   int index_mode,
			   struct perturb_workspace * ppw
			   ) {

  /** Summary: */

  /** - define local variables */

  int index_mt;
  int index_st;
  int index_type;
  int index_eta;
  int number_of_sources;

  /** - define indices of metric perturbations obeying to constraint
        equations (this can be done once and for all, because the
        vector of metric perturbations is the same whatever the
        approximation scheme, unlike the vector of quantities to
        be integrated, which is allocated separately in
        perturb_vector_init) */

  if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

    index_mt = 0;

    /* newtonian gauge */

    if (ppr->gauge == newtonian) {
      ppw->index_mt_phi = index_mt; /* phi */
      index_mt++;
      ppw->index_mt_psi = index_mt; /* psi */
      index_mt++;
      ppw->index_mt_phi_prime = index_mt; /* phi' */
      index_mt++;
    }
      
    /* synchronous gauge (note that eta is counted in the vector of
        quantities to be integrated, while here we only consider
        quantities obeying to constraint equations) */

    if (ppr->gauge == synchronous) {
      ppw->index_mt_h_prime = index_mt; /* h' */
      index_mt++;
      ppw->index_mt_eta_prime = index_mt; /* eta' */
      index_mt++;
      ppw->index_mt_alpha_prime = index_mt; /* alpha' (with alpha = (h' + 6 eta') / (2 k**2) ) */
      index_mt++;

    }     

  }

  if ((ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors)) {

    index_mt = 0;

  }

  ppw->mt_size = index_mt;

  /** - define indices in the vector of source terms */

  index_st = 0;

  ppw->index_st_eta = index_st;
  index_st++;

  ppw->index_st_S0 = index_st;
  index_st++;

  ppw->index_st_S1 = index_st;
  index_st++;

  ppw->index_st_S2 = index_st;
  index_st++;

  ppw->index_st_dS1 = index_st;
  index_st++;

  ppw->index_st_dS2 = index_st;
  index_st++;

  ppw->index_st_ddS2 = index_st;
  index_st++;

  ppw->st_size = index_st;

  /** - allocate some workspace in which we will store temporarily the
        values of background, thermodynamics, metric and source
        quantities at a given time */

  class_alloc(ppw->pvecback,pba->bg_size_normal*sizeof(double),ppt->error_message);
  class_alloc(ppw->pvecthermo,pth->th_size*sizeof(double),ppt->error_message);
  class_alloc(ppw->pvecmetric,ppw->mt_size*sizeof(double),ppt->error_message);

  /** - allocate array of source terms array for the current mode,
        initial condition and wavenumber:
        source_term_table[index_type][index_eta*ppw->st_size+index_st] */

  class_alloc(ppw->source_term_table,
	      ppt->tp_size[index_mode] * sizeof(double *),
	      ppt->error_message);

  for (index_type = 0; index_type < ppt->tp_size[index_mode]; index_type++) {
    class_alloc(ppw->source_term_table[index_type], 
		ppt->eta_size*ppw->st_size*sizeof(double),
		ppt->error_message);
  }

  /** - allocate the perturb_approximation structure */

  class_alloc(ppw->pa,
	      sizeof(struct perturb_approximations),
	      ppt->error_message);

  /** - For definitness, initialize approximation flags to arbitrary
      values (correct values are overwritten during the first call to
      perturb_timescale_and_approximations) */

  ppw->pa->fsa=fsa_off;
  ppw->pa->tca=tca_on;

  return _SUCCESS_;
}

/**
 * Free the perturb_workspace structure (with the exception of the
 * perturb_vector '->pv' field, which is freed separately in
 * perturb_vector_free).
 *
 * @param ppt        Input: pointer to the perturbation structure
 * @param index_mode Input: index of mode under consideration (scalar/.../tensor)
 * @param ppw        Input: pointer to perturb_workspace structure to be freed
 * @return the error status
 */

int perturb_workspace_free (
			    struct perturbs * ppt,
			    int index_mode,
			    struct perturb_workspace * ppw
			    ) {

  int index_type;
  int index_eta;

  free(ppw->pvecback);
  free(ppw->pvecthermo);
  free(ppw->pvecmetric);
  for (index_type = 0; index_type < ppt->tp_size[index_mode]; index_type++) {
    free(ppw->source_term_table[index_type]);
  }
  free(ppw->source_term_table);
  free(ppw->pa);
  free(ppw);

  return _SUCCESS_;
}

/**
 * Solve the perturbation evolution for a given mode, initial
 * condition and wavenumber, and compute the corresponding source
 * functions.
 *
 * For a given mode, initial condition and wavenumber, this function
 * initializes all perturbations using perturb_vector_init(), and
 * integrates them over time. Whenever a "source sampling time" is
 * passed, the source terms are computed and stored temporarily (for
 * each type) using perturb_source_terms(). Whenever the approximation
 * scheme changes, the list of perturbations to be integrated changes
 * thanks to perturb_vector_init(), and the integration goes on.
 * Finally, the actual source functions are computed using the source
 * terms, and stored in the source table using perturb_sources().
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param ppt        Input/Output: pointer to the perturbation structure (output source functions S(k,eta) written here)
 * @param index_mode Input: index of mode under consideration (scalar/.../tensor)
 * @param index_ic   Input: index of initial condition under consideration (ad, iso...)
 * @param index_k    Input: index of wavenumber
 * @param ppw        Input: pointer to perturb_workspace structure containing index values and workspaces
 * @return the error status
 */

int perturb_solve(
		  struct precision * ppr,
		  struct background * pba,
		  struct thermo * pth,
		  struct perturbs * ppt,
		  int index_mode,
		  int index_ic,
		  int index_k,
		  struct perturb_workspace * ppw
		  ) {

  /** Summary: */

  /** - define local variables */ 

  /* contains all quantities relevant for the integration algorithm */
  struct generic_integrator_workspace gi;

  /* contains all fixed parameters, indices and workspaces used by the perturb_derivs function */
  struct perturb_parameters_and_workspace ppaw;

  /* conformal time */
  double eta;

  /* maximum value of conformal time for current wavenumber */
  double etamax;

  int index_eta;
  int eta_actual_size;

  /* size of each step of integration */
  double timestep;
  
  /* running index for the source term vector */
  int index_st;

  /* running index over types (temperature, etc) */
  int index_type;
  
  /* next index in the list of discrete eta values for which sources must be computed (for each type, temperature, polarization, lensing, etc) */
  int next_index_eta;

  /* fourier mode */
  double k;

  /* store the approximation scheme at previous step, in order to see when approximations should be switched on/off */
  struct perturb_approximations * pa_old;

  /* number of changes in the approximation schemes between the
     starting and ending point of one step */
  int num_changing_approximations;

  double timescale;

  /*  time at which  recombination can be considered as over; check that we always integrate at least till this time */
  double eta_visibility_free_streaming;

  /* initialize indices relevant for back/thermo tables search */
  ppw->last_index_back=0;
  ppw->last_index_thermo=0;
  ppw->intermode = normal;

  /* get wavenumber value */
  k = (ppt->k[index_mode])[index_k];

  class_test(k == 0.,
	     ppt->error_message,
	     "stop to avoid division by zero");

  /** - compute conformal time corresponding to end of efficient recombination using background_eta_of_z() */

  class_call(background_eta_of_z(pba,pth->z_visibility_free_streaming,&eta_visibility_free_streaming),
	     pba->error_message,
	     ppt->error_message);

  /** - compute maximum value of eta for which sources are calculated for this wavenumber */

  /* by default, today */
  etamax = pba->conformal_age;
  eta_actual_size = ppt->eta_size;

  /* eventually stop earlier, when k*eta=k_eta_max, but not before the end of recombination */
  if (ppt->has_lss == _FALSE_) {
    if ((ppr->k_eta_max/k < pba->conformal_age) && (ppr->k_eta_max/k > eta_visibility_free_streaming))
      etamax= ppr->k_eta_max/k;
    if ((ppr->k_eta_max/k < eta_visibility_free_streaming) && (eta_visibility_free_streaming < pba->conformal_age))
      etamax = eta_visibility_free_streaming;

    while (ppt->eta_sampling[eta_actual_size-1] > etamax)
      eta_actual_size--;

    class_test(eta_actual_size < 1,
	       ppt->error_message,
	       "did not consider this case yet");

  }

  /** - compute minimum value of eta for which sources are calculated for this wavenumber */

  /* by default, (k eta)_min / k */
  eta = ppr->k_eta_min / k;

  /* force this time to be before the first point defined in the timesampling_for_sources routine */
  if (eta > ppt->eta_sampling[0]*ppr->eta_min_over_sampling_min)
    eta = ppt->eta_sampling[0]*ppr->eta_min_over_sampling_min;
  
  /** - prepare eveything for the integration over time: */

  /** (a) fill the structure containing all fixed parameters, indices
      and workspaces needed by perturb_derivs */

  ppaw.ppr = ppr;
  ppaw.pba = pba;
  ppaw.pth = pth;
  ppaw.ppt = ppt;
  ppaw.index_mode = index_mode;
  ppaw.k = k;
  ppaw.ppw = ppw;

  /** (b) compute background/thermodynamical quantities, and smallest
      relevant time scale in the system using
      perturb_timescale_and_approximations() */

  class_call(perturb_timescale_and_approximations(eta,
						  &ppaw,
						  &timescale,
						  &num_changing_approximations,
						  ppt->error_message),
	     ppt->error_message,
	     ppt->error_message);

  ppaw.ppw->intermode = closeby;

  /** (c) given these background/thermodynamical quantities, fill the
      vector of perturbation variables with appropriate initial
      conditions using perturb_vector_init(), which will call
      perturb_initial_conditions() */

  class_call(perturb_vector_init(ppr,
				 pba,
				 pth,
				 ppt,
				 index_mode,
				 index_ic,
				 k,
				 eta,
				 ppw,
				 NULL),
	     ppt->error_message,
	     ppt->error_message);

  /** (e) initialize 'previous approximation' flags to initial approximations */

  class_alloc(pa_old,
	      sizeof(struct perturb_approximations),
	      ppt->error_message);

  class_call(perturb_copy_approximations(ppw->pa,pa_old),
	     ppt->error_message,
	     ppt->error_message);


  while (eta < ppt->eta_sampling[eta_actual_size-1]) {

    class_call(generic_evolver(perturb_derivs,
			       &eta,
			       ppw->pv->y,
			       ppw->pv->pt_size,
			       &ppaw,
			       ppr->tol_perturb_integration,
			       ppr->smallest_allowed_variation,
			       timescale,
			       perturb_timescale_and_approximations,
			       ppr->perturb_integration_stepsize,
			       ppt->eta_sampling,
			       eta_actual_size,
			       perturb_source_terms,
			       ppt->error_message),
	       ppt->error_message,
	       ppt->error_message);
    
    if (eta < ppt->eta_sampling[eta_actual_size-1]) {

      /* for testing purposes */
      if (ppt->perturbations_verbose > 3) {
	
	class_test(ppw->pvecthermo[pth->index_th_dkappa] == 0.,
		   ppt->error_message,
		   "you have dkappak=0, stop to avoid division by zero");
	
	class_test(ppw->pvecback[pba->index_bg_H]*ppw->pvecback[pba->index_bg_a] == 0.,
		   ppt->error_message,
		   "you have aH=0, stop to avoid division by zero");
	
	if (pa_old->tca==tca_on && ppw->pa->tca==tca_off)
	  printf("Turn off tight-coupling at eta=%e, with k*eta_c=%e and eta_c/eta_h=%e\n",
		 eta,
		 k/ppw->pvecthermo[pth->index_th_dkappa],
		 (ppw->pvecback[pba->index_bg_H]*ppw->pvecback[pba->index_bg_a])/ppw->pvecthermo[pth->index_th_dkappa]);
	
	if (pa_old->tca==tca_off && ppw->pa->tca==tca_on)
	  printf("Turn on tight-coupling  at eta=%e, with k*eta_c=%e\n",
		 eta,
		 k/ppw->pvecthermo[pth->index_th_dkappa]);
	
	if (pa_old->fsa==fsa_off && ppw->pa->fsa==fsa_on)
	  printf("Turn on free-streaming  at eta=%e, with k/aH   =%e and Omega_r    =%e\n",eta,k/(ppw->pvecback[pba->index_bg_H]*ppw->pvecback[pba->index_bg_a]),ppw->pvecback[pba->index_bg_Omega_r]);
	
	if (pa_old->fsa==fsa_on && ppw->pa->fsa==fsa_off)
	  printf("Turn off free-streaming at eta=%e\n",eta);
	
      } 
      /* end of testing zone */
      
      class_call(perturb_vector_init(ppr,
				     pba,
				     pth,
				     ppt,
				     index_mode,
				     index_ic,
				     k,
				     eta,
				     ppw,
				     pa_old),
		 ppt->error_message,
		 ppt->error_message);
      
      class_call(perturb_copy_approximations(ppw->pa,pa_old),
		 ppt->error_message,
		 ppt->error_message);
      
    } 

  }

    /* trivial case (ppt->eta_sampling[next_index_eta] > etamax): just fill source array with zeros and loop till eta_today */

  for (index_eta = eta_actual_size; index_eta < ppt->eta_size; index_eta++) {
    for (index_type = 0; index_type < ppt->tp_size[index_mode]; index_type++) {
      for (index_st = 0; index_st < ppw->st_size; index_st++) {
	ppw->source_term_table[index_type][index_eta*ppw->st_size+index_st] = 0.;
      }
    }
  }

  /** - free other quantitites allocated at the beginning of the routine */

  class_call(perturb_vector_free(ppw->pv),
	     ppt->error_message,
	     ppt->error_message);

  free(pa_old);

  /** - infer source functions from source terms using perturb_sources() */

  class_call(perturb_sources(ppt,
			     index_mode,
			     index_ic,
			     index_k,
			     ppw),
	     ppt->error_message,
	     ppt->error_message);

  return _SUCCESS_;
}
  
/**
 * Initialize the field '->pv' of a perturb_workspace structure, which
 * is a perturb_vector structure. This structure contains indices and
 * values of all quantitites which need to be integrated with respect
 * to time (and only them: quantitites fixed analytically or obeying a
 * constraint equations are NOT included in this vector). This routine
 * distinguishes between two cases:
 *
 * -> the input pa_old is set to the NULL pointer:
 *
 * This happens when we start integrating over a new wavenumber and we
 * want to set initial conditions for the perturbations. Then, it is
 * assumed that ppw->pv is not yet alloacted. This routine allocates
 * it, defines all indices, and then fill the vector ppw->pv->y with
 * the initial conditions defined in perturb_initial_conditions.
 *
 * -> the input pa_old is not set to the NULL pointer and describes
 * some set of approximations:
 *
 * This happens when we need to change of approximation scheme while
 * integrating over a given wavenumber. The new approximation
 * described by ppw->pa is then different from pa_old. Then, this
 * routine allocates a new vector with a new size and new index
 * values; it fills this vector with initial condtions taken from the
 * previous vector passed as an input in ppw->pv, and eventually with
 * some analytic approximations for the new variables appearing at
 * this time; then the new vector comes in replacement of the old one,
 * which is freed.
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param ppt        Input: pointer to the perturbation structure
 * @param index_mode Input: index of mode under consideration (scalar/.../tensor)
 * @param index_ic   Input: index of initial condition under consideration (ad, iso...)
 * @param k          Input: wavenumber
 * @param eta        Input: conformal time
 * @param ppw        Input/Output: workspace containing input the approximation scheme, the background/thermodynamics/metric quantitites, and eventually the previous vector y; and in output the new vector y.
 * @param pa_old     Input: NULL is we need to set y to initial conditions for a new wavnumber; points towards a perturb_approximations if we want to switch of approximation.
 * @return the error status
 */

int perturb_vector_init(
			struct precision * ppr,
			struct background * pba,
			struct thermo * pth,
			struct perturbs * ppt,
			int index_mode,
			int index_ic,
			double k,
			double eta,
			struct perturb_workspace * ppw, /* ppw->pv unallocated if pa_old = NULL, allocated and filled otherwise */
			struct perturb_approximations * pa_old 
			) {

  /** Summary: */

  /** - define local variables */

  struct perturb_vector * ppv;

  int index_pt;
  int l;
  double h_plus_six_eta_prime;

  /** - allocate a new perturb_vector structure to which ppw->pv will point at the end of the routine */

  class_alloc(ppv,sizeof(struct perturb_vector),ppt->error_message);

  /** - defines all indices in this new vector (depends on approximation scheme, described by the input structure ppw->pa) */

  index_pt = 0;

  if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

    /* reject inconsistent values of the number of mutipoles in photon temperature hierachy */
    class_test(ppr->l_max_g < 4,
	       ppt->error_message,
	       "ppr->l_max_g should be at least 4, i.e. we must integrate at least over photon density, velocity, shear, third and fourth momentum");

    /* reject inconsistent values of the number of mutipoles in photon polarization hierachy */
    class_test(ppr->l_max_pol_g < 4,
	       ppt->error_message,
	       "ppr->l_max_pol_g should be at least 4");

    if (ppw->pa->fsa == fsa_off) { /* if free-streaming approximation is off */

      /* photons */

      ppv->index_pt_delta_g = index_pt; /* photon density */
      index_pt++;

      ppv->index_pt_theta_g = index_pt; /* photon velocity */
      index_pt++;

      if (ppw->pa->tca == tca_off) { /* if tight-coupling approximation is off */

	ppv->index_pt_shear_g = index_pt; /* photon shear */
	index_pt++;

	ppv->index_pt_l3_g = index_pt; /* photon l=3 */
	index_pt++;

	ppv->l_max_g = ppr->l_max_g; /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3) */
	index_pt += (ppv->l_max_g-3);

	ppv->index_pt_pol0_g = index_pt; /* photon polarization, l=0 */
	index_pt++;

	ppv->index_pt_pol1_g = index_pt; /* photon polarization, l=1 */
	index_pt++;

	ppv->index_pt_pol2_g = index_pt; /* photon polarization, l=2 */
	index_pt++;

	ppv->index_pt_pol3_g = index_pt; /* photon polarization, l=3 */
	index_pt++;

	ppv->l_max_pol_g = ppr->l_max_pol_g; /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3) */
	index_pt += (ppv->l_max_pol_g-3); 

      }
    }

    /* baryons */

    ppv->index_pt_delta_b = index_pt;  /* baryon density */
    index_pt++;
    
    ppv->index_pt_theta_b = index_pt;  /* baryon velocity */
    index_pt++;

    /* cdm */

    if (pba->has_cdm == _TRUE_) {       

      ppv->index_pt_delta_cdm = index_pt; /* cdm density */
      index_pt++;

      if (ppr->gauge == newtonian) {
	ppv->index_pt_theta_cdm = index_pt; /* cdm velocity */
	index_pt++;
      }
 
    }

    /* dark energy */    

    if (pba->has_dark_energy_fluid == _TRUE_) {       
      
      ppv->index_pt_delta_de = index_pt; /* dark energy density */
      index_pt++;

      ppv->index_pt_theta_de = index_pt; /* dark energy velocity */
      index_pt++;
      
    }
    
    if (ppw->pa->fsa == fsa_off) { /* if free-streaming approximation is off */

      /* ultra relativistic neutrinos */

      if (pba->has_nur == _TRUE_) {

	/* reject inconsistent values of the number of mutipoles in ultra relativistic neutrino hierachy */
	class_test(ppr->l_max_nur < 4,
		   ppt->error_message,
		   "ppr->l_max_nur should be at least 4, i.e. we must integrate at least over neutrino/relic density, velocity, shear, third and fourth momentum");
	
	ppv->index_pt_delta_nur = index_pt; /* density of ultra-relativistic neutrinos/relics */
	index_pt++;
	
	ppv->index_pt_theta_nur = index_pt; /* velocity of ultra-relativistic neutrinos/relics */
	index_pt++;
	
	ppv->index_pt_shear_nur = index_pt; /* shear of ultra-relativistic neutrinos/relics */
	index_pt++;
	
	ppv->index_pt_l3_nur = index_pt; /* l=3 of ultra-relativistic neutrinos/relics */
	index_pt++;
	
	ppv->l_max_nur = ppr->l_max_nur; /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3) */
	index_pt += (ppv->l_max_nur-3);
	
      }
     
    }

    /* metric (only quantitites to be integrated, not those obeying constraint equations) */
 
    if (ppr->gauge == synchronous) {
      ppv->index_pt_eta = index_pt; /* metric perturbation eta of synchronous gauge */
      index_pt++;
    }
  }

  if ((ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors)) {

    /* reject inconsistent values of the number of mutipoles in photon temperature hierachy */
    class_test(ppr->l_max_g_ten < 4,
	       ppt->error_message,
	       "ppr->l_max_g_ten should be at least 4, i.e. we must integrate at least over photon density, velocity, shear, third momentum");
 
    /* reject inconsistent values of the number of mutipoles in photon polarization hierachy */
    class_test(ppr->l_max_pol_g_ten < 2,
	       ppt->error_message,
	       "ppr->l_max_pol_g_ten should be at least 2");

    ppv->index_pt_delta_g = index_pt; /* photon density */
    index_pt++;

    ppv->index_pt_theta_g = index_pt; /* photon velocity */
    index_pt++;

    ppv->index_pt_shear_g = index_pt; /* photon shear */
    index_pt++;

    ppv->index_pt_l3_g = index_pt; /* photon l=3 */
    index_pt++;

    ppv->l_max_g = ppr->l_max_g_ten; /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2) */
    index_pt += (ppv->l_max_g-3);
      
    ppv->index_pt_pol0_g = index_pt; /* photon polarization, l=0 */
    index_pt++;
      
    ppv->index_pt_pol1_g = index_pt; /* photon polarization, l=1 */
    index_pt++;

    ppv->index_pt_pol2_g = index_pt; /* photon polarization, l=2 */
    index_pt++;

    ppv->index_pt_pol3_g = index_pt; /* photon polarization, l=3 */
    index_pt++;

    ppv->l_max_pol_g = ppr->l_max_pol_g_ten; /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2) */
    index_pt += (ppv->l_max_pol_g-3); 

    /** (b) metric perturbation h is a propagating degree of freedom, so h and hdot are included
	in the vector of ordinary perturbations, no in that of metric perturbations */

    ppv->index_pt_gw = index_pt;     /* tensor metric perturbation h (gravitational waves) */
    index_pt++;

    ppv->index_pt_gwdot = index_pt; /* its time-derivative */
    index_pt++;

  }

  ppv->pt_size = index_pt;

  /** - allocate vectors for storing the values of all these
        quantities and their time-derivatives at a given time */

  class_alloc(ppv->y,ppv->pt_size*sizeof(double),ppt->error_message);
  class_alloc(ppv->dy,ppv->pt_size*sizeof(double),ppt->error_message);

  /** - case of setting initial conditions for a new wavenumber */

  if (pa_old == NULL) {

    if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

      /** (a) check that current approximation scheme is consistent
	  with initial conditions */

      class_test(ppw->pa->fsa == fsa_on,
		 ppt->error_message,
		 "scalar initial conditions assume free-streaming approximation turned off");
      
      class_test(ppw->pa->tca == tca_off,
		 ppt->error_message,
		 "scalar initial conditions assume tight-coupling approximation turned on");
      
    }

    /** (b) let ppw->pv points towards the perturb_vector structure
	that we just created */

    ppw->pv = ppv;

   /** (c) fill the vector ppw->pv->y with appropriate initial conditions */

    class_call(perturb_initial_conditions(ppr,
					  pba,
					  ppt,
					  index_mode,
					  index_ic,
					  k,
					  eta,
					  ppw),
	       ppt->error_message,
	       ppt->error_message);
    
  }
    
  /** - case of switching approximation while a wavenumber is being integrated */

  else {

    /** (a) for the scalar mode: */

    if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

      /** -- check that the change of approximation scheme makes
	  sense (note: before calling this routine there is already a
	  check that we wish to change only one approximation flag at
	  a time) */

      class_test((pa_old->tca == tca_off) && (ppw->pa->tca == tca_on),
		 ppt->error_message,
		 "at eta=%g: the tight-coupling approximation can be switched off, not on",eta);

      /* -- case of switching off tight coupling
	 approximation. Provide correct initial conditions to new set
	 of variables */

      if ((pa_old->tca == tca_on) && (ppw->pa->tca == tca_off)) {

	ppv->y[ppv->index_pt_delta_g] =
	  ppw->pv->y[ppw->pv->index_pt_delta_g];

	ppv->y[ppv->index_pt_theta_g] =
	  ppw->pv->y[ppw->pv->index_pt_theta_g];

	if (ppr->gauge == newtonian)
	  ppv->y[ppv->index_pt_shear_g] = (8./3.*ppw->pv->y[ppw->pv->index_pt_theta_g])
	    /9./ppw->pvecthermo[pth->index_th_dkappa]; /* tight-coupling approximation for sshear_g */;
	
	if (ppr->gauge == synchronous)
	  ppv->y[ppv->index_pt_shear_g] = (8./3.*ppw->pv->y[ppw->pv->index_pt_theta_g] + 
					   4./3.*(ppw->pvecmetric[ppw->index_mt_h_prime] + 6. * ppw->pvecmetric[ppw->index_mt_eta_prime]))
	    /9./ppw->pvecthermo[pth->index_th_dkappa]; /* tight-coupling approximation for shear_g */  
	
	ppv->y[ppv->index_pt_l3_g] = 6./7.*k/ppw->pvecthermo[pth->index_th_dkappa]*
	  ppv->y[ppv->index_pt_shear_g]; /* tight-coupling approximation for l=3 */

	for (l = 4; l <= ppv->l_max_g; l++) /* photon additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3) */
	  ppv->y[ppv->index_pt_delta_g+l] = 0.;

	ppv->y[ppv->index_pt_pol0_g] = 2.5*ppv->y[ppv->index_pt_shear_g]; /* tight-coupling approximation for polarization, l=0 */
	ppv->y[ppv->index_pt_pol1_g] = 7./12.*ppv->y[ppv->index_pt_l3_g]; /* tight-coupling approximation for polarization, l=1 */
	ppv->y[ppv->index_pt_pol2_g] = 0.5*ppv->y[ppv->index_pt_shear_g]; /* tight-coupling approximation for polarization, l=2 */
	ppv->y[ppv->index_pt_pol3_g] = 0.25*ppv->y[ppv->index_pt_l3_g];   /* tight-coupling approximation for polarization, l=3 */
	
	for (l = 4; l <= ppv->l_max_pol_g; l++) /* photon additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3) */
	  ppv->y[ppv->index_pt_pol0_g+l] = 0.;

	ppv->y[ppv->index_pt_delta_b] =
	  ppw->pv->y[ppw->pv->index_pt_delta_b];
	
	ppv->y[ppv->index_pt_theta_b] =
	  ppw->pv->y[ppw->pv->index_pt_theta_b];

	if (pba->has_cdm == _TRUE_) {   

	  ppv->y[ppv->index_pt_delta_cdm] =
	    ppw->pv->y[ppw->pv->index_pt_delta_cdm];
	  
	  if (ppr->gauge == newtonian) {
	  ppv->y[ppv->index_pt_theta_cdm] =
	    ppw->pv->y[ppw->pv->index_pt_theta_cdm];
	  }
	}

	if (pba->has_dark_energy_fluid == _TRUE_) {  

	  ppv->y[ppv->index_pt_delta_de] =
	    ppw->pv->y[ppw->pv->index_pt_delta_de];
	  
	  ppv->y[ppv->index_pt_theta_de] =
	    ppw->pv->y[ppw->pv->index_pt_theta_de];
	}

	if (pba->has_nur == _TRUE_) {

	  ppv->y[ppv->index_pt_delta_nur] =
	    ppw->pv->y[ppw->pv->index_pt_delta_nur];

	  ppv->y[ppv->index_pt_theta_nur] =
	    ppw->pv->y[ppw->pv->index_pt_theta_nur];

	  ppv->y[ppv->index_pt_shear_nur] =
	    ppw->pv->y[ppw->pv->index_pt_shear_nur];

	  ppv->y[ppv->index_pt_l3_nur] =
	    ppw->pv->y[ppw->pv->index_pt_l3_nur];

	  for (l=4; l <= ppv->l_max_nur; l++)
	    ppv->y[ppv->index_pt_delta_nur+l] = 
	      ppw->pv->y[ppw->pv->index_pt_delta_nur+l];

	}

	if (ppr->gauge == synchronous)
	  ppv->y[ppv->index_pt_eta] =
	    ppw->pv->y[ppw->pv->index_pt_eta];

      }

      /* -- case of switching on free streaming
	 approximation. Provide correct initial conditions to new set
	 of variables */

      if ((pa_old->fsa == fsa_off) && (ppw->pa->fsa == fsa_on)) {

	ppv->y[ppv->index_pt_delta_b] =
	  ppw->pv->y[ppw->pv->index_pt_delta_b];
	
	ppv->y[ppv->index_pt_theta_b] =
	  ppw->pv->y[ppw->pv->index_pt_theta_b];

	if (pba->has_cdm == _TRUE_) {   

	  ppv->y[ppv->index_pt_delta_cdm] =
	    ppw->pv->y[ppw->pv->index_pt_delta_cdm];
	  
	  if (ppr->gauge != synchronous) {
	  ppv->y[ppv->index_pt_theta_cdm] =
	    ppw->pv->y[ppw->pv->index_pt_theta_cdm];
	  }
	}

	if (pba->has_dark_energy_fluid == _TRUE_) {  

	  ppv->y[ppv->index_pt_delta_de] =
	    ppw->pv->y[ppw->pv->index_pt_delta_de];
	  
	  ppv->y[ppv->index_pt_theta_de] =
	    ppw->pv->y[ppw->pv->index_pt_theta_de];
	}

	if (ppr->gauge == synchronous)
	  ppv->y[ppv->index_pt_eta] =
	    ppw->pv->y[ppw->pv->index_pt_eta];

      }

      /* -- case of switching off free streaming
	 approximation. Provide correct initial conditions to new set
	 of variables */

      if ((pa_old->fsa == fsa_on) && (ppw->pa->fsa == fsa_off)) {

	if (ppr->gauge == newtonian) {
	  ppv->y[ppv->index_pt_delta_g] = 0.;
	  ppv->y[ppv->index_pt_theta_g] = 0.;
	}

	if (ppr->gauge == synchronous) {
	  ppv->y[ppv->index_pt_delta_g] = -4.*ppw->pvecmetric[ppw->index_mt_alpha_prime];
	  ppv->y[ppv->index_pt_theta_g] = -0.5*ppw->pvecmetric[ppw->index_mt_h_prime];
	}

	ppv->y[ppv->index_pt_shear_g] = 0.;
	ppv->y[ppv->index_pt_l3_g] = 0.;
	for (l = 4; l <= ppv->l_max_g; l++)
	  ppv->y[ppv->index_pt_delta_g+l] = 0.;
	ppv->y[ppv->index_pt_pol0_g] = 0.;
	ppv->y[ppv->index_pt_pol1_g] = 0.;
	ppv->y[ppv->index_pt_pol2_g] = 0.;
	ppv->y[ppv->index_pt_pol3_g] = 0.;
	for (l = 4; l <= ppv->l_max_pol_g; l++)
	  ppv->y[ppv->index_pt_pol0_g+l] = 0.;

	ppv->y[ppv->index_pt_delta_b] =
	  ppw->pv->y[ppw->pv->index_pt_delta_b];
	
	ppv->y[ppv->index_pt_theta_b] =
	  ppw->pv->y[ppw->pv->index_pt_theta_b];

	if (pba->has_cdm == _TRUE_) {   

	  ppv->y[ppv->index_pt_delta_cdm] =
	    ppw->pv->y[ppw->pv->index_pt_delta_cdm];
	  
	  if (ppr->gauge != synchronous) {
	  ppv->y[ppv->index_pt_theta_cdm] =
	    ppw->pv->y[ppw->pv->index_pt_theta_cdm];
	  }
	}

	if (pba->has_dark_energy_fluid == _TRUE_) {  

	  ppv->y[ppv->index_pt_delta_de] =
	    ppw->pv->y[ppw->pv->index_pt_delta_de];
	  
	  ppv->y[ppv->index_pt_theta_de] =
	    ppw->pv->y[ppw->pv->index_pt_theta_de];
	}

	if (pba->has_nur == _TRUE_) {

	  if (ppr->gauge == newtonian) {
	    ppv->y[ppv->index_pt_delta_nur] = 0.;
	    ppv->y[ppv->index_pt_theta_nur] = 0.;
	  }

	  if (ppr->gauge == synchronous) {
	    ppv->y[ppv->index_pt_delta_nur] = -4.*ppw->pvecmetric[ppw->index_mt_alpha_prime];
	    ppv->y[ppv->index_pt_theta_nur] = -0.5*ppw->pvecmetric[ppw->index_mt_h_prime];
	  }

	  ppv->y[ppv->index_pt_shear_nur] = 0.;
	  for (l=4; l <= ppv->l_max_nur; l++)
	    ppv->y[ppv->index_pt_delta_nur+l] = 0.;

	}	

	if (ppr->gauge == synchronous)
	  ppv->y[ppv->index_pt_eta] =
	    ppw->pv->y[ppw->pv->index_pt_eta];

      }

    }

    /** (b) for the tensor mode */

    if ((ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors)) {

      class_test(0 ==0,
		 ppt->error_message,
		 "so far, no approximation scheme defined for tensors, so we should never arrive here");
      
    }

    /** (c) free the previous vector of perturbations */
    
    class_call(perturb_vector_free(ppw->pv),
	       ppt->error_message,
	       ppt->error_message);

    /** (d) let ppw->pv points towards the perturb_vector structure
	that we just created */

    ppw->pv = ppv;

  }

  return _SUCCESS_;
}

/**
 * Free the perturb_vector structure.
 *
 * @param pv        Input: pointer to perturb_vector structure to be freed
 * @return the error status
 */

int perturb_vector_free(
			struct perturb_vector * pv
			) {

  free(pv->y);
  free(pv->dy);
  free(pv);
  
  return _SUCCESS_;
}

/**
 * Set the fields of one approximation structure equal to those of another one. 
 *
 * @param pa_input   Input : pointer to input perturb_approximations structure
 * @param pa_copy    Input : pointer to perturb_approximations (already allocated) which fields must become equal to the other ones
 * @param int           Output: number of differences between the two
 * @return the error status
 */

int perturb_copy_approximations(
				   struct perturb_approximations * pa_input,
				   struct perturb_approximations * pa_copy /* already allocated */
				   ) {

  pa_copy->tca = pa_input->tca;
  pa_copy->fsa = pa_input->fsa;

  return _SUCCESS_;

}

/**
 * For each mode, wavenumber and initial condition, this function initializes all values
 * in the vector of perturbed variables. 
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param ppt        Input: pointer to the perturbation structure
 * @param index_mode Input: index of mode under consideration (scalar/.../tensor)
 * @param index_ic   Input: index of initial condition under consideration (ad, iso...)
 * @param k          Input: wavenumber
 * @param eta        Input: conformal time
 * @param ppw        Input/Output: workspace containing input the approximation scheme, the background/thermodynamics/metric quantitites, and eventually the previous vector y; and in output the new vector y.
 * @return the error status
 */
int perturb_initial_conditions(struct precision * ppr,
			       struct background * pba,
			       struct perturbs * ppt,
			       int index_mode,
			       int index_ic,
			       double k,
			       double eta,
			       struct perturb_workspace * ppw
			       ) {
  /** Summary: */

  /** - define local variables */

  /* multipole l */
  int l;

  /** - first set everything to zero */

  ppw->pv->y[ppw->pv->index_pt_delta_g] = 0.; /* photon density */
  ppw->pv->y[ppw->pv->index_pt_theta_g] = 0.; /* photon velocity */

  /* additional perturbations relevant only for the scalar mode */

  if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

    ppw->pv->y[ppw->pv->index_pt_delta_b] = 0.;  /* baryon density */
    ppw->pv->y[ppw->pv->index_pt_theta_b] = 0.;  /* baryon velocity */

    if (pba->has_cdm == _TRUE_) {       
      ppw->pv->y[ppw->pv->index_pt_delta_cdm] = 0.; /* cdm density */
      if (ppr->gauge == newtonian) 
	ppw->pv->y[ppw->pv->index_pt_theta_cdm] = 0.; /* cdm velocity */
    }
    
    if (pba->has_dark_energy_fluid == _TRUE_) {        
      ppw->pv->y[ppw->pv->index_pt_delta_de] = 0.; /* dark energy density */   
      ppw->pv->y[ppw->pv->index_pt_theta_de] = 0.; /* dark energy velocity */ 
    } 
    
    if (pba->has_nur == _TRUE_) {
      ppw->pv->y[ppw->pv->index_pt_delta_nur] = 0; /* density of ultra-relativistic neutrinos/relics */
      ppw->pv->y[ppw->pv->index_pt_theta_nur] = 0; /* velocity of ultra-relativistic neutrinos/relics */
      ppw->pv->y[ppw->pv->index_pt_shear_nur] = 0.; /* shear of ultra-relativistic neutrinos/relics */
      ppw->pv->y[ppw->pv->index_pt_l3_nur] = 0.; /* l=3 of ultra-relativistic neutrinos/relics */
      for (l=4; l <= ppw->pv->l_max_nur; l++)
	ppw->pv->y[ppw->pv->index_pt_delta_nur+l] = 0.;  /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2) */
    }
    
    if (ppr->gauge == synchronous)
      ppw->pv->y[ppw->pv->index_pt_eta] = 0; /* metric perturbation eta */ 

  }

  /* additional perturbations relevant only for the tensor mode */

  if ((ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors)) {

    ppw->pv->y[ppw->pv->index_pt_shear_g] = 0.; /* photon shear */
    ppw->pv->y[ppw->pv->index_pt_l3_g] = 0.; /* photon shear */
    for (l=4; l <= ppw->pv->l_max_g; l++) ppw->pv->y[ppw->pv->index_pt_delta_g+l] = 0.;  /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3) */
    
    ppw->pv->y[ppw->pv->index_pt_pol0_g] = 0.; /* photon polarization, l=0 */
    ppw->pv->y[ppw->pv->index_pt_pol1_g] = 0.; /* photon polarization, l=1 */
    ppw->pv->y[ppw->pv->index_pt_pol2_g] = 0.; /* photon polarization, l=2 */
    ppw->pv->y[ppw->pv->index_pt_pol3_g] = 0.; /* photon polarization, l=3 */
    for (l=4; l <= ppw->pv->l_max_pol_g; l++) ppw->pv->y[ppw->pv->index_pt_pol0_g+l] = 0.;  /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3) */
    
    ppw->pv->y[ppw->pv->index_pt_gw] = 0.;     /* tensor metric perturbation h (gravitational waves) */
    ppw->pv->y[ppw->pv->index_pt_gwdot] = 0.;  /* its time-derivative */
  }

  /** - initial conditions for scalars */

  if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

    /** (a) adiabatic */ 

    if ((ppt->has_ad == _TRUE_) && (index_ic == ppt->index_ic_ad)) {

      /* relevant background quantities */

      /* 8piG/3 rho_r(t_i) */
      double rho_r = ppw->pvecback[pba->index_bg_rho_g];

      /* 8piG/3 rho_m(t_i) */
      double rho_m = ppw->pvecback[pba->index_bg_rho_b];

      /* 8piG/3 rho_nu(t_i) (all neutrinos/relics being relativistic at that time) */
      double rho_nu = 0.;

      if (pba->has_cdm == _TRUE_) {
	rho_m += ppw->pvecback[pba->index_bg_rho_cdm];
      }

      if (pba->has_nur == _TRUE_) {
	rho_r += ppw->pvecback[pba->index_bg_rho_nur];
	rho_nu += ppw->pvecback[pba->index_bg_rho_nur];
      }
      
      class_test(rho_r == 0.,
		 ppt->error_message,
		 "stop to avoid division by zero");

      /* f_nu = Omega_nu(t_i) / Omega_r(t_i) */
      double fracnu = rho_nu/rho_r;

      /* alpha = 5 eta / (15 + 4 f_nu) */ 
      double alpha = 5.*eta/(15.+4.*fracnu);

      /* omega = Omega_m(t_i) a(t_i) H(t_i) / sqrt(Omega_r(t_i))
	 = (8piG/3 rho_m(t_i)) a(t_i) / sqrt(8piG/3 rho_r(t_i))  in Mpc-1 */
      double om = ppw->pvecback[pba->index_bg_a]*rho_m/sqrt(rho_r);

      /* newtonian gauge */
      if (ppr->gauge == newtonian) {

	class_test(eta == 0.,
		   ppt->error_message,
		   "stop to avoid division by zero");

	ppw->pv->y[ppw->pv->index_pt_delta_g] = -4.*alpha/eta - k*k*eta*eta/3.; /* photon density */
	ppw->pv->y[ppw->pv->index_pt_theta_g] = alpha*k*k - pow(k*eta,3.)*k/36.; /* photon velocity */

	ppw->pv->y[ppw->pv->index_pt_delta_b] = 3./4.*ppw->pv->y[ppw->pv->index_pt_delta_g]; /* baryon density */
	ppw->pv->y[ppw->pv->index_pt_theta_b] = ppw->pv->y[ppw->pv->index_pt_theta_g]; /* baryon velocity */
      
	if (pba->has_cdm == _TRUE_) {       
	  ppw->pv->y[ppw->pv->index_pt_delta_cdm] = 3./4.*ppw->pv->y[ppw->pv->index_pt_delta_g]; /* cdm density */
	  ppw->pv->y[ppw->pv->index_pt_theta_cdm] = alpha*k*k; /* cdm velocity */
	}
	
 	if (pba->has_dark_energy_fluid == _TRUE_) {        
 	  ppw->pv->y[ppw->pv->index_pt_delta_de] = 0.; /* dark energy density (TO BE WRITTEN) */
 	  ppw->pv->y[ppw->pv->index_pt_theta_de] = 0.; /* dark energy velocity (TO BE WRITTEN) */
 	} 
	
	if (pba->has_nur == _TRUE_) {
	  ppw->pv->y[ppw->pv->index_pt_delta_nur] = ppw->pv->y[ppw->pv->index_pt_delta_g]; /* density of ultra-relativistic neutrinos/relics */
	  ppw->pv->y[ppw->pv->index_pt_theta_nur] = alpha*k*k - pow(k*eta,3.)*k/36. * (23.+4.*fracnu)/(15.+4.*fracnu); /* velocity of ultra-relativistic neutrinos/relics */
	  ppw->pv->y[ppw->pv->index_pt_shear_nur] = k*k*eta*eta*2./3./(12.+fracnu); /* shear of ultra-relativistic neutrinos/relics */
	}

      }

      /* synchronous gauge */
      if (ppr->gauge == synchronous) {

	ppw->pv->y[ppw->pv->index_pt_delta_g] = - k*k*eta*eta/3. * (1.-om*eta/5.); /* photon density */
	/* ppw->pv->y[ppw->pv->index_pt_theta_g] = - k*k*k*k*eta*eta*eta/36.; /\* photon velocity *\/ */
	ppw->pv->y[ppw->pv->index_pt_theta_g] = k*k*eta/9.*ppw->pv->y[ppw->pv->index_pt_delta_g]; /* photon velocity */

	ppw->pv->y[ppw->pv->index_pt_delta_b] = 3./4.*ppw->pv->y[ppw->pv->index_pt_delta_g]; /* baryon density */
	ppw->pv->y[ppw->pv->index_pt_theta_b] = ppw->pv->y[ppw->pv->index_pt_theta_g]; /* baryon velocity */
      
	if (pba->has_cdm == _TRUE_) {       
	  ppw->pv->y[ppw->pv->index_pt_delta_cdm] = 3./4.*ppw->pv->y[ppw->pv->index_pt_delta_g]; /* cdm density */
          /* by convention, cdm velocity velocity vanishes in this implementation of the synchronous gauge */
	}

 	if (pba->has_dark_energy_fluid == _TRUE_) {        
 	  ppw->pv->y[ppw->pv->index_pt_delta_de] = 0; /* dark energy density (TO BE WRITTEN) */
	  ppw->pv->y[ppw->pv->index_pt_theta_de] = 0; /* dark energy velocity (TO BE WRITTEN) */
 	} 

	if (pba->has_nur == _TRUE_) {
	  ppw->pv->y[ppw->pv->index_pt_delta_nur] = ppw->pv->y[ppw->pv->index_pt_delta_g]; /* density of ultra-relativistic neutrinos/relics */
	  ppw->pv->y[ppw->pv->index_pt_theta_nur] = - pow(k*eta,3.)*k/36. * (23.+4.*fracnu)/(15.+4.*fracnu); /* velocity of ultra-relativistic neutrinos/relics */
	  ppw->pv->y[ppw->pv->index_pt_shear_nur] = k*k*eta*eta*2./3./(12.+fracnu); /* shear of ultra-relativistic neutrinos/relics */
	}    

	ppw->pv->y[ppw->pv->index_pt_eta] = 1.-(5.+4.*fracnu)/12./(15.+4.*fracnu)*k*k*eta*eta; /* metric perturbation eta */

      }

    }

    /** (b) Cold dark matter Isocurvature */ 

    if ((ppt->has_cdi == _TRUE_) && (index_ic == ppt->index_ic_cdi)) { 
      
      class_test(pba->has_cdm == _FALSE_,
		 ppt->error_message,
		 "not consistent to ask for CDI in absence of CDM!");

      ppw->pv->y[ppw->pv->index_pt_delta_cdm] = ppr->entropy_ini;

    }

    /** (c) Baryon Isocurvature */ 

    if ((ppt->has_bi == _TRUE_) && (index_ic == ppt->index_ic_bi)) {

      ppw->pv->y[ppw->pv->index_pt_delta_b] = ppr->entropy_ini;

    }

    /** (d) Neutrino density Isocurvature */ 

    if ((ppt->has_nid == _TRUE_) && (index_ic == ppt->index_ic_nid)) {

      class_test(pba->has_nur == _FALSE_,
		 ppt->error_message,
		 "not consistent to ask for NID in absence of neutrinos!");

      ppw->pv->y[ppw->pv->index_pt_delta_nur] = ppr->entropy_ini;

    }
     
    /** (e) Neutrino velocity Isocurvature */ 

    if ((ppt->has_niv == _TRUE_) && (index_ic == ppt->index_ic_niv)) {

      class_test(pba->has_nur == _FALSE_,
		 ppt->error_message,
		 "not consistent to ask for NIV in absence of neutrinos!");

      ppw->pv->y[ppw->pv->index_pt_theta_nur] = k*ppr->entropy_ini;

    }

  }

  /** - initial conditions for tensors */

  if ((ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors)) {

    if (index_ic == ppt->index_ic_ten) {
      ppw->pv->y[ppw->pv->index_pt_gw] = ppr->gw_ini;
    }

  }

/*   printf("%e %e ",k,eta); */
/*   for (l=0;l<ppw->pv->pt_size;l++) { */
/*     printf("%e \n",ppw->pv->y[l]); */
/*   } */

  return _SUCCESS_;
}

/**
 * Evaluate background/thermodynamics at \f$ \eta \f$, infer useful flags / time scales for integrating perturbations.
 *
 * Evaluate background quantities at \f$ \eta \f$, as well as thermodynamics for scalar mode; infer useful flags and time scales for integrating the perturbations:
 * - check whether tight-coupling approximation is needed.
 * - check whether radiation (photons, massless neutrinos...) perturbations are needed.
 * - choose step of integration: step = ppr->perturb_integration_stepsize * min_time_scale, where min_time_scale = smallest time scale involved in the equations. There are three time scales to compare:
 * -# that of recombination, \f$ \eta_g = 1/\kappa' \f$
 * -# Hubble time scale, \f$ \eta_h = a/a' \f$
 * -# Fourier mode, \f$ \eta_k = 1/k \f$
 *
 * So, in general, min_time_scale = \f$ \min(\eta_g, \eta_b, \eta_h, \eta_k) \f$.
 *
 * However, if \f$ \eta_g \ll \eta_h \f$ and \f$ \eta_g
 * \ll \eta_k \f$, we can use the tight-coupling regime for photons
 * and write equations in such way that the time scale \f$
 * \eta_g \f$ becomes irrelevant (no effective mass term in \f$
 * 1/\eta_g \f$).  Then, the smallest
 * scale in the equations is only \f$ \min(\eta_h, \eta_k) \f$.
 * In practise, it is sufficient to use only the condition \f$ \eta_g \ll \eta_h \f$.
 * 
 * Also, if \f$ \rho_{matter} \gg \rho_{radiation} \f$ and \f$ k \gg
 * aH \f$, we can switch off radiation perturbations (i.e. switch on
 * the free-streaming approximation) and then the smallest scale is
 * simply \f$ \eta_h \f$.
 *
 * @param eta                       Input: conformal time
 * @param parameters_and_workspace  Input/Output: in output contains the updated background/thermo quantitites, as well as the approximations
 * @param timescale  Ouput: smallest relevant time scale in the differential system of perturbations 
 * @param intermode  Input: interpolation mode (normal or growing_closeby)
 * @return the error status
 */

int perturb_timescale_and_approximations(
					double eta,
					void * parameters_and_workspace,
					double * timescale,
					int * num_of_changing_approximations,
					ErrorMsg error_message
					) {
  /** Summary: */

  /** - define local variables */

  /* (a) time scale of Fourier mode, \f$ \eta_k = 1/k \f$ */  
  double eta_k;
  /* (b) time scale of expansion, \f$ \eta_h = a/a' \f$ */
  double eta_h;
  /* (c) time scale of recombination, \f$ \eta_{\gamma} = 1/\kappa' \f$ */
  double eta_g;

  struct perturb_parameters_and_workspace * pppaw;
  struct precision * ppr;
  struct background * pba;
  struct thermo * pth;
  struct perturbs * ppt;
  int index_mode;
  double k;
  struct perturb_workspace * ppw;
  struct perturb_approximations * pa;
  double * pvecback;
  double * pvecthermo;

  pppaw = parameters_and_workspace;
  k = pppaw->k;
  index_mode = pppaw->index_mode;
  ppr = pppaw->ppr;
  pba = pppaw->pba;
  pth = pppaw->pth;
  ppt = pppaw->ppt;
  ppw = pppaw->ppw;
  pa = ppw->pa;
  pvecback = ppw->pvecback;
  pvecthermo = ppw->pvecthermo;

  class_test(k == 0.,
	     ppt->error_message,
	     "stop to avoid division by zero");

  * num_of_changing_approximations = 0;

  /** - compute Fourier mode time scale = \f$ \eta_k = 1/k \f$ */

  eta_k = 1./k;

  /** - evaluate background quantities with background_at_eta() and
        Hubble time scale \f$ \eta_h = a/a' \f$ */

  class_call(background_at_eta(pba,eta, normal_info, ppw->intermode, &(ppw->last_index_back), pvecback),
	     pba->error_message,
	     error_message);

  class_test(pvecback[pba->index_bg_H]*pvecback[pba->index_bg_a] == 0.,
	     error_message,
	     "aH=0, stop to avoid division by zero");

  eta_h = 1./(pvecback[pba->index_bg_H]*pvecback[pba->index_bg_a]);

  /** - for scalars modes: */

  if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {
    
    /** (a) if \f$ k \gg aH \f$ and \f$ \Omega_r \ll 1 \f$ and
        efficient recombination finished, switch off radiation
        perturbations and keep Hubble time scale as the relevant
        one. Otherwise take \f$ \min(\eta_k, \eta_h). \f$ */

    if ((eta_h/eta_k > ppr->rad_pert_trigger_k_over_aH) && 
	//      (eta > eta_visibility_free_streaming) && /* optionally this line could be restored, to check that this does not happen before recombination is completed) */
	(pvecback[pba->index_bg_Omega_r] < ppr->rad_pert_trigger_Omega_r)) {
      
      if (pa->fsa == fsa_off) (*num_of_changing_approximations)++;
      pa->fsa = fsa_on;
      *timescale = eta_h;
    }
    else {
      if (pa->fsa == fsa_on) (*num_of_changing_approximations)++;
      pa->fsa = fsa_off;
      if (eta_k < eta_h) 
	*timescale = eta_k;
      else  
	*timescale = eta_h;
    }

    /** (b) evaluate thermodynamical quantities with thermodynamics_at_z() */

    class_call(thermodynamics_at_z(pba,
				   pth,
				   1./pvecback[pba->index_bg_a]-1.,  /* redshift z=1/a-1 */
				   ppw->intermode,
				   &(ppw->last_index_thermo),
				   pvecback,
				   pvecthermo),
	       pth->error_message,
	       error_message);
    
    /** (b.1.) if \f$ \kappa'=0 \f$, recombination is finished; check that tight-coupling approximation is off */

    if (pvecthermo[pth->index_th_dkappa] == 0.) {
      if (pa->tca == tca_on) (*num_of_changing_approximations)++;
      pa->tca = tca_off;
      class_test(pvecthermo[pth->index_th_dkappa] == 0.,
		 error_message,
		 "This test is just for debugging, can be erased\n");
    }

    /** (b.2.) if \f$ \kappa' \neq 0 \f$, recombination is not finished: */

    else {

      /** (b.2.a) compute recombination time scale for photons, \f$ \eta_{\gamma} = 1/ \kappa' \f$ */
      eta_g = 1./pvecthermo[pth->index_th_dkappa];

      /** (b.2.b) if \f$ \eta_g \ll max(\eta_h, \eta_k) \f$, turn on the tight-coupling approximation for photons; otherwise turn off  and define relevant time scale as min time scale between \f$ \eta_k \f$, \f$ \eta_h \f$ and \f$ \eta_g \f$ */
      if ((eta_g/eta_h < ppr->tight_coupling_trigger_eta_g_over_eta_h) 
	  && (eta_g/eta_k < ppr->tight_coupling_trigger_eta_g_over_eta_k)) {
	if (pa->tca == tca_off) (*num_of_changing_approximations)++;
	pa->tca = tca_on;
      }
      else {
	if (pa->tca == tca_on) (*num_of_changing_approximations)++;
	pa->tca = tca_off;
	if (eta_g < *timescale) *timescale=eta_g;
      }
    }
  }

  /** - for tensor modes: tight-coupling approximation is off, time scale remains \f$ min (\eta_k , \eta_h) \f$. */

  if ((ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors)) {
    if (pa->tca == tca_on) (*num_of_changing_approximations)++;
    pa->tca = tca_off;
    if (pa->fsa == fsa_on) (*num_of_changing_approximations)++;
    pa->fsa = fsa_off;
    if (eta_k < eta_h) 
      *timescale = eta_k;
    else  
      *timescale = eta_h;
  }
  
  /* vectors not coded yet */
  
  return _SUCCESS_;
}

/**
 * Compute metric perturbations (those not integrated over time) using Einstein equations
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param ppt        Input: pointer to the perturbation structure
 * @param index_mode Input: index of mode under consideration (scalar/.../tensor)
 * @param k          Input: wavenumber
 * @param eta        Input: conformal time
 * @param y          Input: vector of perturbations (those integrated over time) (already allocated)
 * @param ppw        Input/Output: in output contains the updated metric perturbations
 * @return the error status
 */

int perturb_einstein(
		     struct precision * ppr,
		     struct background * pba,
		     struct perturbs * ppt,
		     int index_mode,
		     double k,
		     double eta,
		     double * y,
		     struct perturb_workspace * ppw
		     ) {
  /** Summary: */

  /** - define local variables */

  double k2,a,a2,a_prime_over_a;
  double delta_rho,delta_theta,delta_shear;
  double delta_g, theta_g, shear_g;
  double delta_nur, theta_nur, shear_nur;

  /** - wavenumber and scale factor related quantities */ 

  k2 = k*k;
  a = ppw->pvecback[pba->index_bg_a];
  a2 = a * a;
  a_prime_over_a = ppw->pvecback[pba->index_bg_H]*a;

  /** - for scalar modes: */  

  if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

    /** (a) deal with approximation schemes */
    
    if (ppw->pa->tca == tca_off) {
      if (ppw->pa->fsa == fsa_off) {
	delta_g = y[ppw->pv->index_pt_delta_g];
	theta_g = y[ppw->pv->index_pt_theta_g];
	shear_g = y[ppw->pv->index_pt_shear_g];
	if (pba->has_nur == _TRUE_) {
	  delta_nur = y[ppw->pv->index_pt_delta_nur];
	  theta_nur = y[ppw->pv->index_pt_theta_nur];
	  shear_nur = y[ppw->pv->index_pt_shear_nur];
	}
      }
      else {
	delta_g = 0.;
	theta_g = 0.; /* the free-streaming approximation is much
			 better if we use an approximation for theta_g
			 instead of 0. We will correct for this effect
			 below, after the computation of delta_rho. */
	shear_g = 0.;
	if (pba->has_nur == _TRUE_) {
	  delta_nur = 0.;
	  theta_nur = 0.; /* the free-streaming approximation is much
			 better if we use an approximation for
			 theta_nur instead of 0. We will correct for
			 this effect below, after the computation of
			 delta_rho. */
	  shear_nur = 0.;
	}
      }
    }
    else {
      delta_g = y[ppw->pv->index_pt_delta_g];
      theta_g = y[ppw->pv->index_pt_theta_g];
      shear_g = 0.; /* in the tight-coupling approximation, shear_g is
		       a function of h' and eta'; but h' and eta' are
		       calculated below as a function of delta_g and
		       theta_g, not as a function of shear_g itself;
		       hence it is consistent to set shear_g to zero
		       here, and to its tight-coupling value in
		       perturb_derivs. In the synchronous gauge,
		       neglecting shear_g here results in an
		       inaccurate calculation of alpha' below. However
		       this alpha' is only used for computing the
		       sources. Since sources are not sampled during
		       the tight-coupling regime, we are perfectely
		       safe at least in the synchronous gauge. */
      if (pba->has_nur == _TRUE_) {
	delta_nur = y[ppw->pv->index_pt_delta_nur];
	theta_nur = y[ppw->pv->index_pt_theta_nur];
	shear_nur = y[ppw->pv->index_pt_shear_nur];
      }
    }

    /** (b) compute the total density, velocity and shear perturbations */
 
    /* photon and baryon contribution */
    delta_rho = ppw->pvecback[pba->index_bg_rho_g]*delta_g
      + ppw->pvecback[pba->index_bg_rho_b]*y[ppw->pv->index_pt_delta_b];
    delta_theta = 4./3.*ppw->pvecback[pba->index_bg_rho_g]*theta_g
      + ppw->pvecback[pba->index_bg_rho_b]*y[ppw->pv->index_pt_theta_b];
    delta_shear = 4./3.*ppw->pvecback[pba->index_bg_rho_g]*shear_g;

    /* cdm contribution */
    if (pba->has_cdm == _TRUE_) {
      delta_rho = delta_rho + ppw->pvecback[pba->index_bg_rho_cdm]*y[ppw->pv->index_pt_delta_cdm];
      if (ppr->gauge == newtonian)
	delta_theta = delta_theta + ppw->pvecback[pba->index_bg_rho_cdm]*y[ppw->pv->index_pt_theta_cdm];
    }
    
    /* dark energy fluid contribution */
    if (pba->has_dark_energy_fluid == _TRUE_) { 
      delta_rho = delta_rho + ppw->pvecback[pba->index_bg_rho_de]*y[ppw->pv->index_pt_delta_de]; 
      delta_theta = delta_theta + ppw->pvecback[pba->index_bg_rho_de]*y[ppw->pv->index_pt_theta_de];
    } 

    /* ultra-relativistic neutrino/relics contribution */
    if (pba->has_nur == _TRUE_) {
      delta_rho = delta_rho + ppw->pvecback[pba->index_bg_rho_nur]*delta_nur;
      delta_theta = delta_theta + 4./3.*ppw->pvecback[pba->index_bg_rho_nur]*theta_nur;
      delta_shear = delta_shear + 4./3.*ppw->pvecback[pba->index_bg_rho_nur]*shear_nur;
    }

    /** (c) eventually correct for the free-streaming velocities */
    if (ppw->pa->fsa == fsa_on) {

      if (ppr->gauge == newtonian) { /*TBD*/ }

      /* theta_g = theta_nur = -h_prime/2 */
      if (ppr->gauge == synchronous) {
	theta_g = -0.5 * ( k2 * y[ppw->pv->index_pt_eta] + 1.5 * a2 * delta_rho)/(0.5*a_prime_over_a);
	delta_theta += 4./3.*ppw->pvecback[pba->index_bg_rho_g]*theta_g;
	if (pba->has_nur == _TRUE_) {
	  theta_nur = theta_g;
	  delta_theta += 4./3.*ppw->pvecback[pba->index_bg_rho_nur]*theta_nur;
	}
      }
    }

    /** (d) infer metric perturbations from Einstein equations */

    /* newtonian gauge */
    if (ppr->gauge == newtonian) {
      ppw->pvecmetric[ppw->index_mt_phi] = -1.5 * (a2/k2/k2) * (k2 * delta_rho + 3.*a_prime_over_a * delta_theta); /* phi */
      ppw->pvecmetric[ppw->index_mt_psi] = ppw->pvecmetric[ppw->index_mt_phi] - 4.5 * (a2/k2) * delta_shear;  /* psi */
      ppw->pvecmetric[ppw->index_mt_phi_prime] = - a_prime_over_a * ppw->pvecmetric[ppw->index_mt_psi] + 1.5 * (a2/k2) * delta_theta; /* phi' */
    }

    /* synchronous gauge */
    if (ppr->gauge == synchronous) {
      ppw->pvecmetric[ppw->index_mt_eta_prime] = 1.5 * (a2/k2) * delta_theta;  /* eta' */
      ppw->pvecmetric[ppw->index_mt_h_prime] = 
	( k2 * y[ppw->pv->index_pt_eta] + 1.5 * a2 * delta_rho)/(0.5*a_prime_over_a);  /* h' */
      ppw->pvecmetric[ppw->index_mt_alpha_prime] = 
	- 4.5 * (a2/k2) * delta_shear + y[ppw->pv->index_pt_eta] - 2.*a_prime_over_a*
	(ppw->pvecmetric[ppw->index_mt_h_prime] + 6.*ppw->pvecmetric[ppw->index_mt_eta_prime])
	/ 2./ k2; /* alpha' = (h''+6eta'')/2k2 *

      /* getting phi here is an option */
      /* phi=y[ppw->pv->index_pt_eta]-0.5 * (a_prime_over_a/k2) * (h_plus_six_eta_prime); */   /* phi from gauge transformation (from synchronous to newtonian) */

    }
  }

  /* nothing to be done for tensors: only one propagating degree of freedom, no constraint equation */

  return _SUCCESS_;

}

/**
 * Compute the terms contributing to the source functions
 *
 * Compute the terms contributing to the source functions (for all
 * requested types) and store results in global variable
 * pvecsourceterms. The source functions can be decomposed as 
 * \f$ S = S_0 + S_1' + S_2'' \f$.  
 * This function computes \f$ ( S_0, S_1', S_2') \f$ for temperature
 * and \f$ ( S_0, S_1', S_2'') \f$ for other quantitites.
 *
 * @param eta   Input: conformal time  
 * @param pppaw Input/Output: in input, all parameters needed by perturb_derivs, in ourput, source terms
 * @return the error status
 */

int perturb_source_terms(
			 double eta,
			 double * pvecperturbations,
			 double * pvecderivs,
			 int index_eta,
			 void * parameters_and_workspace,
			 ErrorMsg error_message
			 ) {
  /** Summary: */

  /** - define local variables */

  double k2,a_prime_over_a,a_primeprime_over_a,R;
  double Pi,Pi_prime,Psi,Psi_prime;
  double x;
  int index_type;

  struct perturb_parameters_and_workspace * pppaw;
  struct precision * ppr;
  struct background * pba;
  struct thermo * pth;
  struct perturbs * ppt;
  int index_mode;
  double k;
  struct perturb_workspace * ppw;
  struct perturb_approximations * pa;
  double * pvecback;
  double * pvecthermo;
  double * pvecmetric;
  double ** source_term_table;

  double delta_g;

  /** - rename structure fields (just to avoid heavy notations) */

  pppaw = parameters_and_workspace;
  ppr = pppaw->ppr;
  pba = pppaw->pba;
  pth = pppaw->pth;
  ppt = pppaw->ppt;
  index_mode = pppaw->index_mode;
  k = pppaw->k;
  ppw = pppaw->ppw;

  pvecback = ppw->pvecback;
  pvecthermo = ppw->pvecthermo;
  pvecmetric = ppw->pvecmetric;
  source_term_table = ppw->source_term_table;
  pa = ppw->pa;

  class_test(pa->tca == tca_on,
	     ppt->error_message,
	     "source calculation assume tight-coupling approximation turned off");

  x = k * (pba->conformal_age-eta);

  /** - get background/thermo quantities in this point */

  class_call(background_at_eta(pba,
			       eta, 
			       normal_info, 
			       closeby, 
			       &(ppw->last_index_back), 
			       pvecback),
	     pba->error_message,
	     error_message);

  class_call(thermodynamics_at_z(pba,
				 pth,
				 1./pvecback[pba->index_bg_a]-1.,  /* redshift z=1/a-1 */
				 closeby,
				 &(ppw->last_index_thermo),
				 pvecback,
				 pvecthermo),
	     pth->error_message,
	     error_message);

  /* scalars */
  if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

    class_call(perturb_einstein(ppr,
				pba,
				ppt,
				index_mode,
				k,
				eta,
				pvecperturbations,
				ppw),
	       ppt->error_message,
	       error_message);

    if (pa->fsa == fsa_on) {
      delta_g = -4.*pvecmetric[ppw->index_mt_alpha_prime];
      Pi = 0.;
      Pi_prime = 0.;
    }
    else {
      delta_g = pvecperturbations[ppw->pv->index_pt_delta_g];
      Pi = pvecperturbations[ppw->pv->index_pt_pol0_g] + pvecperturbations[ppw->pv->index_pt_pol2_g] + 2.*pvecperturbations[ppw->pv->index_pt_shear_g];
      Pi_prime = pvecderivs[ppw->pv->index_pt_pol0_g] + pvecderivs[ppw->pv->index_pt_pol2_g] + 2.*pvecderivs[ppw->pv->index_pt_shear_g];
    }

    /** - compute \f$ k^2 \f$, \f$ \Pi = G_{\gamma 0} + G_{\gamma 2} + F_{\gamma 2} \f$, \f$ e^{- \kappa} \f$ */ 
      
    k2 = k * k;
    
    a_prime_over_a = pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a];

    a_primeprime_over_a = pvecback[pba->index_bg_H_prime] * pvecback[pba->index_bg_a] 
      + 2. * a_prime_over_a * a_prime_over_a;

    R = 4./3. * pvecback[pba->index_bg_rho_g]/pvecback[pba->index_bg_rho_b];

    /** - compute \f$ k^2 \f$, \f$ \Pi = G_{\gamma 0} + G_{\gamma 2} + F_{\gamma 2} \f$, \f$ e^{- \kappa} \f$ */ 

    /** - for each type and each mode, compute S0, S1, S2 */
    for (index_type = 0; index_type < ppt->tp_size[index_mode]; index_type++) {

      source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_eta] = eta;

      /* scalar temperature */
      if ((ppt->has_source_t == _TRUE_) && (index_type == ppt->index_tp_t)) {

        /* check that visibility is non-zero (otherwise source = 0) */
	if (pvecthermo[pth->index_th_g] != 0.) {

          /* newtonian gauge */
	  if (ppr->gauge == newtonian) {

            /* S0 */
	    source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_S0] =
	      pvecthermo[pth->index_th_exp_m_kappa] * pvecmetric[ppw->index_mt_phi_prime]
	      + pvecthermo[pth->index_th_g] / 4. * (delta_g + Pi / 4.);
	    
            /* S1 */
	    source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_S1] =
	      pvecthermo[pth->index_th_exp_m_kappa] * pvecmetric[ppw->index_mt_psi]
	      + pvecthermo[pth->index_th_g] * pvecperturbations[ppw->pv->index_pt_theta_b] / k2;

	    /* S2 */
	    source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_S2] =
	      3./16. * pvecthermo[pth->index_th_g] * Pi / k2;

	  }

          /* synchronous gauge */
	  if (ppr->gauge == synchronous) {

	    /* S0 */
	    source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_S0] =
	      pvecthermo[pth->index_th_exp_m_kappa] * pvecmetric[ppw->index_mt_eta_prime]
	      + pvecthermo[pth->index_th_g] / 4. * (delta_g + Pi / 4.);
	  
	    /* S1 */
	    /* source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_S1] = */
	    /* 	      pvecthermo[pth->index_th_g] * pvecperturbations[ppw->pv->index_pt_theta_b] / k2; */

	    /* dS1 */
	    source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_dS1] =
	      pvecthermo[pth->index_th_dg] * pvecperturbations[ppw->pv->index_pt_theta_b] / k2
	      + pvecthermo[pth->index_th_g] * pvecderivs[ppw->pv->index_pt_theta_b] / k2;

	    /* S2 */
	    /* source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_S2] = */
	    /* 	      pvecthermo[pth->index_th_exp_m_kappa] * (pvecmetric[ppw->index_mt_h_prime] + 6. * pvecmetric[ppw->index_mt_eta_prime])/2./k2  */
	    /* 	      + 3./16. * pvecthermo[pth->index_th_g] * Pi / k2; */

            /* dS2 */
	    source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_dS2] =
	      pvecthermo[pth->index_th_g] 
	      * (pvecmetric[ppw->index_mt_h_prime] + 6. * pvecmetric[ppw->index_mt_eta_prime])/2./k2
	      + pvecthermo[pth->index_th_exp_m_kappa] * (pvecmetric[ppw->index_mt_alpha_prime])
	      + 3./16. * pvecthermo[pth->index_th_dg] * Pi / k2
	      + 3./16. * pvecthermo[pth->index_th_g] * Pi_prime / k2;

	    /* 	    Pi_prime = -k*pvecperturbations[ppw->pv->index_pt_pol0_g+1] */
	    /* 	      -pvecthermo[pth->index_th_dkappa] * (pvecperturbations[ppw->pv->index_pt_pol0_g]-Pi/2.) */
	    /* 	      +k/5. * (2.*pvecperturbations[ppw->pv->index_pt_pol2_g-1]-3.*pvecperturbations[ppw->pv->index_pt_pol2_g+1]) */
	    /* 	      -pvecthermo[pth->index_th_dkappa] * (pvecperturbations[ppw->pv->index_pt_pol2_g]-Pi/10.) */
	    /*               +8./15.*pvecperturbations[ppw->pv->index_pt_theta_g] */
	    /* 	      -3./5.*k*pvecperturbations[ppw->pv->index_pt_shear_g+1] */
	    /* 	      +4./15.*pvecmetric[ppw->index_mt_h_prime] */
	    /* 	      +8./5.*pvecmetric[ppw->index_mt_eta_prime] */
	    /* 	      -pvecthermo[pth->index_th_dkappa] * (2.*pvecperturbations[ppw->pv->index_pt_shear_g]-1./10.*Pi); */

	    /* 	    Pi_prime = -3./10.*pvecthermo[pth->index_th_dkappa]*Pi */
	    /* 	      -3./5.*k * (pvecperturbations[ppw->pv->index_pt_pol1_g]+ */
	    /* 				pvecperturbations[ppw->pv->index_pt_pol2_g+1]+ */
	    /* 				pvecperturbations[ppw->pv->index_pt_shear_g+1]) */
	    /* 	      +8./15.*pvecperturbations[ppw->pv->index_pt_theta_g] */
	    /* 	      +4./15. * (pvecmetric[ppw->index_mt_h_prime] + 6. * pvecmetric[ppw->index_mt_eta_prime]); */

	  }
	}

	else {

	  source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_S0] = 0.;
	  source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_dS1] = 0.;
	  source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_dS2] = 0.;

	}
      }

      /* scalar polarization */
      if ((ppt->has_source_e == _TRUE_) && (index_type == ppt->index_tp_e)) {

	if (x > 0.) {
	  /* in all gauges */
	  source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_S0] =
	    + 3./16. * pvecthermo[pth->index_th_g] * Pi /x/x;  
	}
	else {
	  source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_S0] = 0.;
	}
	source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_dS1] = 0.;
	source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_ddS2] = 0.;

      }

      /* gravitational potential */
      if ((ppt->has_source_g == _TRUE_) && (index_type == ppt->index_tp_g)) {
      
	/* newtonian gauge */
	if (ppr->gauge == newtonian) {
	  source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_S0] = 
	    pvecmetric[ppw->index_mt_phi];
	  source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_dS1] = 0.;
	  source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_ddS2] = 0.;
	}

	/* synchronous gauge */
	if (ppr->gauge == synchronous) {
	  source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_S0] = 
	    (a_prime_over_a * (pvecmetric[ppw->index_mt_h_prime] + 6. * pvecmetric[ppw->index_mt_eta_prime])/2./k2 + pvecmetric[ppw->index_mt_alpha_prime]);
	  source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_dS1] = 0.;
	  source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_ddS2] = 0.;

	  /* testing zone */
/* 	  if ((k>1.e-3) && (k<2.e-3)) */
/* 	    printf("%g %g %g %g %g %g\n", */
/* 		   k, */
/* 		   pvecback[pba->index_bg_a],    */
/* 		   a_prime_over_a*pvecperturbations[ppw->pv->index_pt_theta_b], */
/* 		   pvecthermo[pth->index_th_cb2]*k2*pvecperturbations[ppw->pv->index_pt_delta_b], */
/* 		   k2*pvecperturbations[ppw->pv->index_pt_eta], */
/* 		   0.5*a_prime_over_a*pvecmetric[ppw->index_mt_h_prime] */
/* 		   ); */
	 

	}
      }
    }

  }

  /* tensors */
  if ((ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors)) {

    Psi=pvecperturbations[ppw->pv->index_pt_delta_g]/40.
      +pvecperturbations[ppw->pv->index_pt_shear_g]*2./35.
      +pvecperturbations[ppw->pv->index_pt_delta_g+4]/210.
      -pvecperturbations[ppw->pv->index_pt_pol0_g]*3./5. 
      +pvecperturbations[ppw->pv->index_pt_pol2_g]*6./35.
      -pvecperturbations[ppw->pv->index_pt_pol0_g+4]/210.;

    Psi_prime=pvecderivs[ppw->pv->index_pt_delta_g]/40.
      +pvecderivs[ppw->pv->index_pt_shear_g]*2./35.
      +pvecderivs[ppw->pv->index_pt_delta_g+4]/210.
      -pvecderivs[ppw->pv->index_pt_pol0_g]*3./5.  
      +pvecderivs[ppw->pv->index_pt_pol2_g]*6./35. 
      -pvecderivs[ppw->pv->index_pt_pol0_g+4]/210.;

    /** - for each type and each mode, compute S0, S1, S2 */
    for (index_type = 0; index_type < ppt->tp_size[index_mode]; index_type++) {

      source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_eta] = eta;

      /* tensor temperature */
      if ((ppt->has_source_t == _TRUE_) && (index_type == ppt->index_tp_t)) {

	if (x > 0.) {
	  source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_S0] = 
	    (-pvecperturbations[ppw->pv->index_pt_gwdot]*pvecthermo[pth->index_th_exp_m_kappa]
	    +pvecthermo[pth->index_th_g]*Psi)/x/x;
	}
	else {
	  source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_S0] = 0.;
	}

	source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_dS1] = 0.;
	source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_ddS2] = 0.;

      }

      /* tensor polarization */
      if ((ppt->has_source_e == _TRUE_) && (index_type == ppt->index_tp_e)) {

	if (x > 0.) {

	  source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_S0] = 
	    (1.-2./x/x)*pvecthermo[pth->index_th_g]*Psi;
	  source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_dS1] = 
	    -4./x*(pvecthermo[pth->index_th_g]*(Psi/x+Psi_prime/k)+pvecthermo[pth->index_th_dg]*Psi/k);
	  source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_dS2] = 
	    -(pvecthermo[pth->index_th_g]*Psi_prime+pvecthermo[pth->index_th_dg]*Psi)/k/k;

	}

	else {

	  source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_S0] = 0.;
	  source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_dS1] = 0.;
	  source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_dS2] = 0.;

	}

      }

      if ((ppt->has_source_b == _TRUE_) && (index_type == ppt->index_tp_b)) {

	if (x > 0.) {
	
	  source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_S0] = 
 	    -pvecthermo[pth->index_th_g]*(4.*Psi/x+2.*Psi_prime/k)-2.*pvecthermo[pth->index_th_dg]*Psi/k;
	}
	
	else {
	  source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_S0] = 0.;
	}
	
	source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_dS1] = 0.;
	source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_ddS2] = 0.;
      }
    }
  }

  return _SUCCESS_;

}

/**
 * Infer source functions from source terms.
 * 
 * For a given mode, initial condition, wavenumber and source type,
 * take the array of source terms \f$ (S_0, S_1, S_2) \f$ at all times
 * (stored in source_terms_array) and infer the source function \f$ S = S_0 + S_1' +
 * S_2'' \f$ by numerical differentiation
 *
 * @param source_terms_array Input/Output: array of source terms
 * @return the error status
 */
int perturb_sources(
		    struct perturbs * ppt,
		    int index_mode,
		    int index_ic,
		    int index_k,
		    struct perturb_workspace * ppw
		    ) {
  /** Summary: */

  /** - define local variables */

  int index_eta,index_type;

  for (index_type = 0; index_type < ppt->tp_size[index_mode]; index_type++) {

    /** - for scalar temperature, infer \f$ S_2'' \f$ from \f$ S_2' \f$ at each time with array_derive1_order2_table_line_to_line() */

    if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {
      if ((ppt->has_source_t == _TRUE_) && (index_type == ppt->index_tp_t)) {

	/* before computing numerical derivatives, slice out the end of the table if filled with zeros */
	index_eta = ppt->eta_size-1;
	while ((ppw->source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_dS2] == 0.) && (index_eta > 0))
	  index_eta--;
	
	/* numerical derivative */
	class_call(array_derive1_order2_table_line_to_line(
							   ppt->eta_sampling,
							   index_eta+1,
							   ppw->source_term_table[index_type],
							   ppw->st_size,
							   ppw->index_st_dS2,
							   ppw->index_st_ddS2,
							   ppt->error_message),
		   ppt->error_message,
		   ppt->error_message);
      }
    }

    /** - for tensor E-polarization, infer \f$ S_2'' \f$ from \f$ S_2' \f$ at each time with array_derive1_order2_table_line_to_line() */

    if ((ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors)) {
      if ((ppt->has_source_e == _TRUE_) && (index_type == ppt->index_tp_e)) {

	/* before computing numerical derivatives, slice out the end of the table if filled with zeros */
	index_eta = ppt->eta_size-1;
	while ((ppw->source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_dS2] == 0.) && (index_eta > 0))
	  index_eta--;
	
	/* numerical derivative */
	class_call(array_derive1_order2_table_line_to_line(
							   ppt->eta_sampling,
							   index_eta+1,
							   ppw->source_term_table[index_type],
							   ppw->st_size,
							   ppw->index_st_dS2,
							   ppw->index_st_ddS2,
							   ppt->error_message),
		   ppt->error_message,
		   ppt->error_message);
      }
    }

    /** - for each time, sum up \f$ S = S_0 + S_1' + S_2'' \f$ and store in array ((sources[index_mode])[index_ic][index_type])[index_eta][index_k] */

    for (index_eta = 0; index_eta < ppt->eta_size; index_eta++) {

      ppt->sources[index_mode]
	[index_ic * ppt->tp_size[index_mode] + index_type]
	[index_eta * ppt->k_size[index_mode] + index_k] = 
	ppw->source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_S0]
	+ppw->source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_dS1]
	+ppw->source_term_table[index_type][index_eta * ppw->st_size + ppw->index_st_ddS2];

    }

  } /* end of loop over types */

  return _SUCCESS_;
}

/**
 * Compute derivative of all perturbations to be integrated 
 *
 * For each mode (scalar/vector/tensor) and each wavenumber k, this
 * function computes the derivative of all values in the vector of
 * perturbed variables to be integrated.
 *
 * This is one of the few functions in the code which are passed to the generic_integrator() routine. 
 * Since generic_integrator() should work with functions passed from various modules, the format of the arguments
 * is a bit special:
 * - fixed parameters and workspaces are passed through a generic pointer. 
 *   generic_integrator() doesn't know what the content of this pointer.
 * - the error management is a bit special: errors are not written as usual to pth->error_message, but to a generic 
 *   error_message passed in the list of arguments.
 *
 * @param eta Input: conformal time
 * @param y Input: vector of perturbations
 * @param dy Ouput: vector of its derivatives (already allocated)
 * @param parameters_and_workspace Input/Output: in input, fixed parameters (e.g. indices); in output, background and thermo quantities evaluated at eta.
 * @param error_message Output : error message
 */
int perturb_derivs(double eta,       /**< Input : conformal time */
		   double * y,       /**< Input : vector of perturbations */
		   double * dy, /**< Output : derivative of vector of perturbations */
		   void * parameters_and_workspace,
		   ErrorMsg error_message
		   ) {
  /** Summary: */

  /** - define local variables */

  /* multipole */
  int l;
  
  /* squared wavenumber */
  double k2;  

  /* scale factor and related quantities */ 
  double a,a2,a_prime_over_a,a_primeprime_over_a,z;

  /* background density ratios */
  double R,fracnu;

  /* useful terms for tight-coupling approximation */
  double slip,Pi;

  /* useful term for tensors */
  double Psi;

  /* useful combination of synchronous metric perturbations \f$ (h' + 6 \eta') \f$ */
  double h_plus_six_eta_prime;

  struct perturb_parameters_and_workspace * pppaw;

  double k;
  int index_mode;
  struct precision * ppr;
  struct background * pba;
  struct thermo * pth;
  struct perturbs * ppt;
  struct perturb_workspace * ppw;
  struct perturb_approximations * pa;
  double * pvecback;
  double * pvecthermo;
  double * pvecmetric;

  double theta_g;
  double shear_g;

  /** - rename structure fields (just to avoid heavy notations) */

  pppaw = parameters_and_workspace;
  k = pppaw->k;
  index_mode = pppaw->index_mode;
  ppr = pppaw->ppr;
  pba = pppaw->pba;
  pth = pppaw->pth;
  ppt = pppaw->ppt;
  ppw = pppaw->ppw;
  pa = ppw->pa;
  pvecback = ppw->pvecback;
  pvecthermo = ppw->pvecthermo;
  pvecmetric = ppw->pvecmetric;

  k2 = k*k;

  /** - get background/thermo quantities in this point */

  class_call(background_at_eta(pba,
			       eta, 
			       normal_info, 
			       closeby, 
			       &(ppw->last_index_back), 
			       pvecback),
	     pba->error_message,
	     error_message);

  class_call(thermodynamics_at_z(pba,
				 pth,
				 1./pvecback[pba->index_bg_a]-1.,  /* redshift z=1/a-1 */
				 closeby,
				 &(ppw->last_index_thermo),
				 pvecback,
				 pvecthermo),
	     pth->error_message,
	     error_message);

  /** - compute related background quantities */
  k2 = k*k;
  a = pvecback[pba->index_bg_a];
  a2 = a * a;
  a_prime_over_a = pvecback[pba->index_bg_H] * a;
  a_primeprime_over_a = pvecback[pba->index_bg_H_prime] * a + 2. * a_prime_over_a * a_prime_over_a;
  z = 1./a-1.;
  R = 4./3. * pvecback[pba->index_bg_rho_g]/pvecback[pba->index_bg_rho_b];
  fracnu = pvecback[pba->index_bg_rho_nur] / (pvecback[pba->index_bg_rho_g] + pvecback[pba->index_bg_rho_nur]);

  /** - for scalar mode: */
  if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

    /** (a) get metric perturbations with perturb_einstein() */
    class_call(perturb_einstein(ppr,
				pba,
				ppt,
				index_mode,
				k,
				eta,
				y,
				ppw),
	       ppt->error_message,
	       error_message);

    /* compute metric-related quantities */
    h_plus_six_eta_prime = pvecmetric[ppw->index_mt_h_prime] + 6. * pvecmetric[ppw->index_mt_eta_prime]
;

    /*******************************************************/

     /** (b) if some approximation schemes are turned on, enforce a few y[] values */ 

     /** (b.b) fix by hand photon/neutrino density/velocity/shear when free-streaming approximation is on (exclusive of photon tight-coupling being on) */ 
/*     if (pa->fsa == fsa_on) { */

       /* analytic free-streaming solution (during matter domination and inside Hubble) */ 
/*       if (ppr->gauge == newtonian) { */
 	/* Newtonian gauge : */

/*       } */
/*       if (ppr->gauge == synchronous) { */
 	/* Synchronous gauge : */ 
/* 	delta_g = -4.*pvecmetric[ppw->index_mt_alpha_prime]; */
/* 	theta_g = -0.5*pvecmetric[ppw->index_mt_h_prime]; */
/* 	shear_g = 0.; */
/* 	if (pba->has_nur == _TRUE_) { */
/* 	  delta_nur = -4.*pvecmetric[ppw->index_mt_alpha_prime]; */
/* 	  theta_nur = -0.5*pvecmetric[ppw->index_mt_h_prime]; */
/* 	  shear_nur = 0.; */
/* 	} */
/*       } */
/*     } */

    /************************************************************/

    /* order:  (delta_g)     if fsa_off
                delta_b
                theta_b
	       (theta_g)     if fsa_off
	      ((shear_g))    if fsa_off and tca_off
	      ((other's_g))  if fsa_off abd tca_off
	      other species */


    /** (c) Photon temperature density (does not depend on tca) */

    if (pa->fsa == fsa_off) {

      if (ppr->gauge == newtonian)
	  dy[ppw->pv->index_pt_delta_g] = /* photon density */
	    -4./3.*y[ppw->pv->index_pt_theta_g]+4.*pvecmetric[ppw->index_mt_phi_prime];
	  
      if (ppr->gauge == synchronous)
	dy[ppw->pv->index_pt_delta_g] = /* photon density */
	  -4./3.*y[ppw->pv->index_pt_theta_g] - 2./3.*pvecmetric[ppw->index_mt_h_prime];
      
    }

    /** (c) baryon density (does not depend on tca) */

    if (ppr->gauge == newtonian)
      dy[ppw->pv->index_pt_delta_b] = /* baryon density */
	-y[ppw->pv->index_pt_theta_b] + 3.*pvecmetric[ppw->index_mt_phi_prime];
    
    if (ppr->gauge == synchronous)
      dy[ppw->pv->index_pt_delta_b] = /* baryon density */
	-y[ppw->pv->index_pt_theta_b] - 0.5*pvecmetric[ppw->index_mt_h_prime];
    
    /** (d) Baryon velocity (depends on tight-coupling approximation) */

    /** (d.1) Baryon velocity if tight-coupling and free-streaming are off */

    if (pa->tca == tca_off) {

      if (pa->fsa == fsa_off) {

	theta_g = y[ppw->pv->index_pt_theta_g];

      }

      else {

	if (ppr->gauge == newtonian)
	  theta_g = 0.; /* TBD */

	if (ppr->gauge == synchronous)
	  theta_g = -0.5*pvecmetric[ppw->index_mt_h_prime];

      }
      
      if (ppr->gauge == newtonian)
	/* Newtonian gauge : */
	dy[ppw->pv->index_pt_theta_b] = /* baryon velocity */
	  - a_prime_over_a*y[ppw->pv->index_pt_theta_b] 
	  + k2*pvecmetric[ppw->index_mt_psi] 
	  + pvecthermo[pth->index_th_cb2]*k2*y[ppw->pv->index_pt_delta_b]
	  + R*pvecthermo[pth->index_th_dkappa] * (theta_g-y[ppw->pv->index_pt_theta_b]);

      if (ppr->gauge == synchronous)
	/* Synchronous gauge : */
	dy[ppw->pv->index_pt_theta_b] = /* baryon velocity */
	  - a_prime_over_a*y[ppw->pv->index_pt_theta_b] 
	  + pvecthermo[pth->index_th_cb2]*k2*y[ppw->pv->index_pt_delta_b]
	  + R*pvecthermo[pth->index_th_dkappa] * (theta_g-y[ppw->pv->index_pt_theta_b]);
    
    }
  
    /** (d.2) Baryon velocity if baryon tight-coupling is on */
    
    else {

      if (ppr->gauge == newtonian) {
	/* Newtonian gauge : */
	slip=(2.*R/(1.+R)*a_prime_over_a+pvecthermo[pth->index_th_ddkappa]/pvecthermo[pth->index_th_dkappa]) /* tight-coupling (theta_b-theta_g)' */
	  *(y[ppw->pv->index_pt_theta_b]-y[ppw->pv->index_pt_theta_g])
	  +(-a_primeprime_over_a*y[ppw->pv->index_pt_theta_b]
	    -a_prime_over_a*k2*(y[ppw->pv->index_pt_delta_g]/2.+pvecmetric[ppw->index_mt_psi])
	    +k2*(pvecthermo[pth->index_th_cb2]*dy[ppw->pv->index_pt_delta_b]
		 -dy[ppw->pv->index_pt_delta_g]/4.)
	    )/pvecthermo[pth->index_th_dkappa]/(1.+R);

	shear_g=(8./3.*y[ppw->pv->index_pt_theta_g])/9./pvecthermo[pth->index_th_dkappa]; /* tight-coupling shear_g */

	dy[ppw->pv->index_pt_theta_b] = /* tight-coupling baryon velocity */
	  (-a_prime_over_a*y[ppw->pv->index_pt_theta_b]
	   +pvecthermo[pth->index_th_cb2]*k2*y[ppw->pv->index_pt_delta_b]
	   +k2*R*(y[ppw->pv->index_pt_delta_g]/4.-shear_g)
	   +R*slip)/(1.+R)
	  +k2*pvecmetric[ppw->index_mt_psi];
      }

      if (ppr->gauge == synchronous) {
	/* Synchronous gauge : */
	slip=(2.*R/(1.+R)*a_prime_over_a+pvecthermo[pth->index_th_ddkappa]/pvecthermo[pth->index_th_dkappa]) /* tight-coupling (theta_b-theta_g)' */
	  *(y[ppw->pv->index_pt_theta_b]-y[ppw->pv->index_pt_theta_g])
	  +(-a_primeprime_over_a*y[ppw->pv->index_pt_theta_b]
	    -a_prime_over_a*k2*(y[ppw->pv->index_pt_delta_g]/2.)
	    +k2*(pvecthermo[pth->index_th_cb2]*dy[ppw->pv->index_pt_delta_b]
		 -dy[ppw->pv->index_pt_delta_g]/4.)
	    )/pvecthermo[pth->index_th_dkappa]/(1.+R);

	/* tight-coupling (theta_b-theta_g)' */
	/* 	slip=2.*R/(1.+R)*a_prime_over_a*(y[ppw->pv->index_pt_theta_b]-y[ppw->pv->index_pt_theta_g]) */
	/* 	  +(-a_primeprime_over_a*y[ppw->pv->index_pt_theta_b] */
	/* 	    -a_prime_over_a*k2*y[ppw->pv->index_pt_delta_g]/2. */
	/* 	    +k2*(pvecthermo[pth->index_th_cb2]*dy[ppw->pv->index_pt_delta_b] */
	/* 		 -dy[ppw->pv->index_pt_delta_g]/4.) */
	/* 	    )/pvecthermo[pth->index_th_dkappa]/(1.+R); */

	/* for testing */
	/*printf("%e %e\n",1./a-1.,pvecthermo[pth->index_th_ddkappa]/pvecthermo[pth->index_th_dkappa]);*/

	shear_g=(8./3.*y[ppw->pv->index_pt_theta_g]+4./3.*h_plus_six_eta_prime)/9./pvecthermo[pth->index_th_dkappa]; /* tight-coupling shear_g */ 

	dy[ppw->pv->index_pt_theta_b] = /* tight-coupling baryon velocity */
	  (-a_prime_over_a*y[ppw->pv->index_pt_theta_b]
	   +pvecthermo[pth->index_th_cb2]*k2*y[ppw->pv->index_pt_delta_b]
	   +k2*R*(y[ppw->pv->index_pt_delta_g]/4.-shear_g)
	   +R*slip)/(1.+R);
      }
    }

    /** (e) Photon temperature higher momenta and photon polarisation (depend on tight-coupling approximation) : */

    if (pa->fsa == fsa_off) {

      /** (e.1) if photon tight-coupling is off: */ 
      if (pa->tca == tca_off) {

	/** (e.1.a) define \f$ \Pi = G_{\gamma 0} + G_{\gamma 2} + F_{\gamma 2} \f$ */
	Pi = y[ppw->pv->index_pt_pol0_g] + y[ppw->pv->index_pt_pol2_g] + 2.*y[ppw->pv->index_pt_shear_g];

	/** (e.1.b) Photon velocity and shear */ 

	if (ppr->gauge == newtonian) {
	  /* Newtonian gauge : */
	  dy[ppw->pv->index_pt_theta_g] = /* photon velocity */
	    k2*(y[ppw->pv->index_pt_delta_g]/4.
		-y[ppw->pv->index_pt_shear_g]+pvecmetric[ppw->index_mt_psi])
	    +pvecthermo[pth->index_th_dkappa]*(y[ppw->pv->index_pt_theta_b]
					       -y[ppw->pv->index_pt_theta_g]);

	  dy[ppw->pv->index_pt_shear_g] = /* photon shear */
	    0.5*(8./15.*y[ppw->pv->index_pt_theta_g]
		 -3./5.*k*y[ppw->pv->index_pt_l3_g]
		 -pvecthermo[pth->index_th_dkappa]*(2.*y[ppw->pv->index_pt_shear_g]-1./10.*Pi));
	}
      
	if (ppr->gauge == synchronous) {
	  /* Synchronous gauge : */
	  dy[ppw->pv->index_pt_theta_g] = /* photon velocity */
	    k2*(y[ppw->pv->index_pt_delta_g]/4.
		-y[ppw->pv->index_pt_shear_g])
	    + pvecthermo[pth->index_th_dkappa]*(y[ppw->pv->index_pt_theta_b]
						-y[ppw->pv->index_pt_theta_g]);
	  
	  dy[ppw->pv->index_pt_shear_g] = /* photon shear */
	    0.5*(8./15.*y[ppw->pv->index_pt_theta_g]
		 -3./5.*k*y[ppw->pv->index_pt_l3_g]
		 +4./15.*pvecmetric[ppw->index_mt_h_prime]+8./5.*pvecmetric[ppw->index_mt_eta_prime]
		 -pvecthermo[pth->index_th_dkappa]*(2.*y[ppw->pv->index_pt_shear_g]-1./10.*Pi));
	}
      
	/** (e.1.c) Photon temperature higher momenta (l >=3), gauge-independent */ 

	l = 3; /* photon l=3 (special case because F_gamma2=2*shear !!) */
	dy[ppw->pv->index_pt_l3_g] =
	  k/(2.*l+1.)*(l*2.*y[ppw->pv->index_pt_shear_g]-(l+1.)*y[ppw->pv->index_pt_l3_g+1])
	  - pvecthermo[pth->index_th_dkappa]*y[ppw->pv->index_pt_l3_g];

	for (l = 4; l < ppw->pv->l_max_g; l++) { /* photon additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3) */
	  dy[ppw->pv->index_pt_delta_g+l] =
	    k/(2.*l+1)*(l*y[ppw->pv->index_pt_delta_g+l-1]-(l+1.)*y[ppw->pv->index_pt_delta_g+l+1])
	    - pvecthermo[pth->index_th_dkappa]*y[ppw->pv->index_pt_delta_g+l];
	}

	l = ppw->pv->l_max_g; /* l=lmax */
	dy[ppw->pv->index_pt_delta_g+ppw->pv->l_max_g] = /* last photon term */
	  k*y[ppw->pv->index_pt_delta_g+ppw->pv->l_max_g-1]
	  -(1.+l)/eta*y[ppw->pv->index_pt_delta_g+ppw->pv->l_max_g]
	  - pvecthermo[pth->index_th_dkappa]*y[ppw->pv->index_pt_delta_g+ppw->pv->l_max_g];
	
	/** (e.1.d) Photon polarisation */

	dy[ppw->pv->index_pt_pol0_g] = /* photon polarization, l=0 */
	  -k*y[ppw->pv->index_pt_pol0_g+1]
	  -pvecthermo[pth->index_th_dkappa]*(y[ppw->pv->index_pt_pol0_g]-Pi/2.);
	dy[ppw->pv->index_pt_pol1_g] = /* photon polarization, l=1 */
	  k/3.*(y[ppw->pv->index_pt_pol1_g-1]-2.*y[ppw->pv->index_pt_pol1_g+1])
	  -pvecthermo[pth->index_th_dkappa]*y[ppw->pv->index_pt_pol1_g];
	dy[ppw->pv->index_pt_pol2_g] = /* photon polarization, l=2 */
	  k/5.*(2.*y[ppw->pv->index_pt_pol2_g-1]-3.*y[ppw->pv->index_pt_pol2_g+1])
	  -pvecthermo[pth->index_th_dkappa]*(y[ppw->pv->index_pt_pol2_g]-Pi/10.);

	for (l=3; l < ppw->pv->l_max_pol_g; l++)
	  dy[ppw->pv->index_pt_pol0_g+l] =  /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2) */
	    k/(2.*l+1)*(l*y[ppw->pv->index_pt_pol0_g+l-1]-(l+1.)*y[ppw->pv->index_pt_pol0_g+l+1])
	    -pvecthermo[pth->index_th_dkappa]*y[ppw->pv->index_pt_pol0_g+l];

	l = ppw->pv->l_max_pol_g;
	dy[ppw->pv->index_pt_pol0_g+l] =  /* l=lmax */
	  k*y[ppw->pv->index_pt_pol0_g+l-1]-(l+1)/eta*y[ppw->pv->index_pt_pol0_g+l]
	  -pvecthermo[pth->index_th_dkappa]*y[ppw->pv->index_pt_pol0_g+l];

      }

      /** (e.2) if photon tight-coupling is on: */
 
      else {

	/** (e.2.a) photon velocity */

	if (ppr->gauge == newtonian)
	  /* Newtonian gauge : */
	  dy[ppw->pv->index_pt_theta_g] = /* tight-coupling photon velocity */
	    -(dy[ppw->pv->index_pt_theta_b]+a_prime_over_a*y[ppw->pv->index_pt_theta_b]-pvecthermo[pth->index_th_cb2]*k2*y[ppw->pv->index_pt_delta_b])/R
	    +k2*(0.25*y[ppw->pv->index_pt_delta_g]-shear_g)+(1.+R)/R*k2*pvecmetric[ppw->index_mt_psi];
	
	if (ppr->gauge == synchronous)
	  /* Synchronous gauge : */
	  
	  dy[ppw->pv->index_pt_theta_g] = /* tight-coupling photon velocity */
	    -(dy[ppw->pv->index_pt_theta_b]+a_prime_over_a*y[ppw->pv->index_pt_theta_b]-pvecthermo[pth->index_th_cb2]*k2*y[ppw->pv->index_pt_delta_b])/R
	    +k2*(y[ppw->pv->index_pt_delta_g]/4.-shear_g);
      
	/*       	  dy[ppw->pv->index_pt_theta_g] = /\* photon velocity (used here in CAMB, leading to more instabilities) *\/ */
	/*       	    k2*(y[ppw->pv->index_pt_delta_g]/4.-y[ppw->pv->index_pt_shear_g]) */
	/*       	    + pvecthermo[pth->index_th_dkappa]*(y[ppw->pv->index_pt_theta_b] */
	/*       						 -y[ppw->pv->index_pt_theta_g]); */

      }
    }

    /** (f) cdm */

    if (pba->has_cdm == _TRUE_) {  

      if (ppr->gauge == newtonian) {
	/* Newtonian gauge : */
	dy[ppw->pv->index_pt_delta_cdm] = /* cdm density */
	  -y[ppw->pv->index_pt_theta_cdm]+3.*pvecmetric[ppw->index_mt_phi_prime];
	dy[ppw->pv->index_pt_theta_cdm] = /* cdm velocity */
	  - a_prime_over_a*y[ppw->pv->index_pt_theta_cdm] + k2*pvecmetric[ppw->index_mt_psi];
      }

      if (ppr->gauge == synchronous) {
	/* Synchronous gauge : */
	dy[ppw->pv->index_pt_delta_cdm] = /* cdm density */
	  -0.5*pvecmetric[ppw->index_mt_h_prime];
      }

    }
    
    /** (g) dark energy fluid */
    
    if (pba->has_dark_energy_fluid == _TRUE_) {  

      double ache_prime = a_prime_over_a;
      double cs2 = 1.;

      if (ppr->gauge == newtonian) {
	/* Newtonian gauge : */
	dy[ppw->pv->index_pt_delta_de] = /* dark energy density */
	  (-3*(1+ pba->w_de )*ache_prime-3*pvecback[pba->index_bg_H]*(cs2- pba->w_de )*(y[ppw->pv->index_pt_delta_de]/pvecback[pba->index_bg_rho_de]+3*pvecback[pba->index_bg_H]*(1+ pba->w_de )*y[ppw->pv->index_pt_theta_de]/k)-(1+ pba->w_de )*k*y[ppw->pv->index_pt_theta_de])/pvecback[pba->index_bg_rho_de]; // 0;

	dy[ppw->pv->index_pt_theta_de] = /* dark energy velocity */
	  (k*cs2*y[ppw->pv->index_pt_delta_de])/(pvecback[pba->index_bg_rho_de]*(1+ pba->w_de ))-pvecback[pba->index_bg_H]*(1-3*cs2)*y[ppw->pv->index_pt_theta_de]+k*pvecmetric[ppw->index_mt_psi]; // 0;
      }

      if (ppr->gauge == synchronous) {
	/* Synchronous gauge : */
	dy[ppw->pv->index_pt_delta_de] = /* dark energy density */
	  (-3*(1+ pba->w_de )*ache_prime-3*pvecback[pba->index_bg_H]*(cs2- pba->w_de )*(y[ppw->pv->index_pt_delta_de]/pvecback[pba->index_bg_rho_de]+3*pvecback[pba->index_bg_H]*(1+ pba->w_de )*y[ppw->pv->index_pt_theta_de]/k)-(1+ pba->w_de )*k*y[ppw->pv->index_pt_theta_de])/pvecback[pba->index_bg_rho_de]; // 0;

	dy[ppw->pv->index_pt_theta_de] = /* dark energy velocity */
	  (k*cs2*y[ppw->pv->index_pt_delta_de])/(pvecback[pba->index_bg_rho_de]*(1+ pba->w_de ))-pvecback[pba->index_bg_H]*(1-3*cs2)*y[ppw->pv->index_pt_theta_de]; // 0;
      }
      
    }  
    
    /** (h) ultra-relativistic neutrino/relics density, velocity, shear, etc. */
    
    if(pa->fsa == fsa_off) {

      if (pba->has_nur == _TRUE_) {
      
	if (ppr->gauge == newtonian) {

	  /* Newtonian gauge : */
	  dy[ppw->pv->index_pt_delta_nur] = /* density of ultra-relativistic neutrinos/relics */
	    -4./3.*y[ppw->pv->index_pt_theta_nur] + 4.*pvecmetric[ppw->index_mt_phi_prime];
	  
	  dy[ppw->pv->index_pt_theta_nur] = /* velocity of ultra-relativistic neutrinos/relics */
	    k2*(y[ppw->pv->index_pt_delta_nur]/4.
		-y[ppw->pv->index_pt_shear_nur]+pvecmetric[ppw->index_mt_psi]);
	  
	  dy[ppw->pv->index_pt_shear_nur] = /* shear of ultra-relativistic neutrinos/relics */
	    0.5*(8./15.*y[ppw->pv->index_pt_theta_nur]
		 -3./5.*k*y[ppw->pv->index_pt_shear_nur+1]);
	  
	}

	if (ppr->gauge == synchronous) {

	  /* Synchronous gauge : */
	  dy[ppw->pv->index_pt_delta_nur] = /* density of ultra-relativistic neutrinos/relics */
	    -4./3.*y[ppw->pv->index_pt_theta_nur] - 2./3.*pvecmetric[ppw->index_mt_h_prime];
	  
	  dy[ppw->pv->index_pt_theta_nur] = /* velocity of ultra-relativistic neutrinos/relics */
	    k2*(y[ppw->pv->index_pt_delta_nur]/4.
		-y[ppw->pv->index_pt_shear_nur]);
	  
	  dy[ppw->pv->index_pt_shear_nur] = /* shear of ultra-relativistic neutrinos/relics */
	    0.5*(8./15.*y[ppw->pv->index_pt_theta_nur]
		 -3./5.*k*y[ppw->pv->index_pt_shear_nur+1]
		 +4./15.*pvecmetric[ppw->index_mt_h_prime]+8./5.*pvecmetric[ppw->index_mt_eta_prime]);
	}

	l = 3;
	dy[ppw->pv->index_pt_l3_nur] = /* l=3 of ultra-relativistic neutrinos/relics (special case because F_gamma2=2*shear !!) */
	  k/(2.*l+1.)*(l*2.*y[ppw->pv->index_pt_shear_nur]-(l+1.)*y[ppw->pv->index_pt_l3_nur+1]);

	for (l = 4; l < ppw->pv->l_max_nur; l++) {
	  dy[ppw->pv->index_pt_delta_nur+l] = /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3) */
	    k/(2.*l+1)*(l*y[ppw->pv->index_pt_delta_nur+l-1]-(l+1.)*y[ppw->pv->index_pt_delta_nur+l+1]);
	}

	l = ppw->pv->l_max_nur; /* l=lmax */
	dy[ppw->pv->index_pt_delta_nur+ppw->pv->l_max_nur] = /* last term of ultra-relativistic neutrinos/relics */
	  k/(2.*l+1)*(l*y[ppw->pv->index_pt_delta_nur+ppw->pv->l_max_nur-1]-(l+1.)*
		      ((2.*l+1)/k/eta*y[ppw->pv->index_pt_delta_nur+ppw->pv->l_max_nur]-y[ppw->pv->index_pt_delta_nur+ppw->pv->l_max_nur-1]));

      }
    }

    /** (j) metric */

    if (ppr->gauge == synchronous) {
      /* Synchronous gauge */
      dy[ppw->pv->index_pt_eta] = pvecmetric[ppw->index_mt_eta_prime];
    }


      /* for testing, will be useful for improving tight-coupling approximation */
      /*       if ((index_k == 0) && (eta > 800)) */
      /* 	printf("%e %e %e %e %e %e %e %e %e %e %e\n", */
      /* 	       eta, */
      /* 	       y[ppw->pv->index_pt_delta_g], */
      /* 	       y[ppw->pv->index_pt_theta_g], */
      /* 	       y[ppw->pv->index_pt_shear_g], */
      /* 	       y[ppw->pv->index_pt_l3_g], */
      /* 	       y[ppw->pv->index_pt_pol0_g], */
      /* 	       y[ppw->pv->index_pt_pol1_g], */
      /* 	       y[ppw->pv->index_pt_pol2_g], */
      /* 	       y[ppw->pv->index_pt_pol3_g], */
      /* 	       y[ppw->pv->index_pt_delta_b], */
      /* 	       y[ppw->pv->index_pt_theta_b]); */

  }

  /** - tensor mode */

  if ((ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors)) {
      
    Psi = 
      y[ppw->pv->index_pt_delta_g]/40.
      +2.*y[ppw->pv->index_pt_shear_g]/35.
      +y[ppw->pv->index_pt_delta_g+4]/210.
      -3.*y[ppw->pv->index_pt_pol0_g]/5. 
      +6.*y[ppw->pv->index_pt_pol2_g]/35.
      -y[ppw->pv->index_pt_pol0_g+4]/210.;

    /* photon density (4*F_0) */
    dy[ppw->pv->index_pt_delta_g] = 
      -4./3.*y[ppw->pv->index_pt_theta_g]
      -4.*y[ppw->pv->index_pt_gwdot]
      -pvecthermo[pth->index_th_dkappa]*(y[ppw->pv->index_pt_delta_g]-4.*Psi);

    /* photon velocity ((3k/4)*F_1) */
    dy[ppw->pv->index_pt_theta_g] = 
      k2*(y[ppw->pv->index_pt_delta_g]/4.-y[ppw->pv->index_pt_shear_g])
      -pvecthermo[pth->index_th_dkappa]*y[ppw->pv->index_pt_theta_g];

    /* photon shear (0.5*F_2) */
    dy[ppw->pv->index_pt_shear_g] =	
      0.5*(8./15.*y[ppw->pv->index_pt_theta_g]
	   -3./5.*k*y[ppw->pv->index_pt_shear_g+1])
      -pvecthermo[pth->index_th_dkappa]*y[ppw->pv->index_pt_shear_g];

    /* photon l=3 */
    dy[ppw->pv->index_pt_l3_g] = 
      k/7.*(6.*y[ppw->pv->index_pt_shear_g]-4.*y[ppw->pv->index_pt_l3_g+1])
      -pvecthermo[pth->index_th_dkappa]*y[ppw->pv->index_pt_l3_g];

    /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3,4) */
    for (l=4; l < ppw->pv->l_max_g; l++)
      dy[ppw->pv->index_pt_delta_g+l] = 
	k/(2.*l+1.)*(l*y[ppw->pv->index_pt_delta_g+l-1]-(l+1.)*y[ppw->pv->index_pt_delta_g+l+1])
	-pvecthermo[pth->index_th_dkappa]*y[ppw->pv->index_pt_delta_g+l];  

    /* l=lmax */
    l = ppw->pv->l_max_g;
    dy[ppw->pv->index_pt_delta_g+l] = 
      k*y[ppw->pv->index_pt_delta_g+l-1]
      -(1.+l)/eta*y[ppw->pv->index_pt_delta_g+l]
      - pvecthermo[pth->index_th_dkappa]*y[ppw->pv->index_pt_delta_g+l];

    /* photon polarization, l=0 */
    dy[ppw->pv->index_pt_pol0_g] = 
      -k*y[ppw->pv->index_pt_pol0_g+1]
      -pvecthermo[pth->index_th_dkappa]*(y[ppw->pv->index_pt_pol0_g]+Psi); 
    
    /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3,4) */
    for (l=1; l < ppw->pv->l_max_pol_g; l++)
      dy[ppw->pv->index_pt_pol0_g+l] = 
	k/(2.*l+1.)*(l*y[ppw->pv->index_pt_pol0_g+l-1]-(l+1.)*y[ppw->pv->index_pt_pol0_g+l+1])
	-pvecthermo[pth->index_th_dkappa]*y[ppw->pv->index_pt_pol0_g+l];

    /* l=lmax */
    l = ppw->pv->l_max_pol_g;
    dy[ppw->pv->index_pt_pol0_g+l] = 
      k*y[ppw->pv->index_pt_pol0_g+l-1]
      -(l+1.)/eta*y[ppw->pv->index_pt_pol0_g+l]
      -pvecthermo[pth->index_th_dkappa]*y[ppw->pv->index_pt_pol0_g+l];

      /* tensor metric perturbation h (gravitational waves) */
    dy[ppw->pv->index_pt_gw] = y[ppw->pv->index_pt_gwdot];     

    /* its time-derivative */
    dy[ppw->pv->index_pt_gwdot] = -2.*a_prime_over_a*y[ppw->pv->index_pt_gwdot]-k2*y[ppw->pv->index_pt_gw];

/*     fprintf(stderr, */
/* 	    "%g %g %g %g %g %g %g\n", */
/* 	    y[ppw->pv->index_pt_gw], */
/* 	    y[ppw->pv->index_pt_gwdot], */
/* 	    y[ppw->pv->index_pt_delta_g], */
/* 	    y[ppw->pv->index_pt_theta_g], */
/* 	    y[ppw->pv->index_pt_shear_g], */
/* 	    y[ppw->pv->index_pt_l3_g], */
/* 	    y[ppw->pv->index_pt_l3_g+1]); */

/*     class_test(0==1,error_message,"stop here\n"); */
	    
  }

  /*     printf("Leaves derivs with:\n"); */
  /*     printf("gamma : %e %e %e %e %e %e \n",dy[0],dy[1],dy[2],dy[3],dy[4],dy[5]); */
  /*     printf("b     : %e %e \n",dy[6],dy[7]); */
  /*     printf("cdm   : %e \n",dy[8]); */
  /*     printf("dark energy : %e %e \n",dy[9],dy[10]); */
  /*     printf("nu    : %e %e %e %e %e %e \n",dy[10],dy[11],dy[12],dy[13],dy[14],dy[15]); */
  /*     printf("eta   : %e \n",dy[16]); */
  /*     printf("h     : %e \n",pvecmetric[ppw->index_mt_h_prime]); */

  return _SUCCESS_;
}

