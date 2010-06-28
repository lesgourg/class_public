/** @file perturbations.c Documented perturbation module
 * Julien Lesgourgues, 18.04.2010    
 *
 * Deals with the perturbation evolution.
 * This module has two purposes: 
 *
 * - at the beginning, to initialize the perturbations, i.e. to
 integrate the perturbation equations, and store temporarily the terms
 contributing to the source functions as a function of conformal
 time. Then, to perform a few manipulations of these terms in order to
 infer the actual source functions \f$ S^{X} (k, \eta) \f$, and to
 store them as a function of conformal time inside an interpolation
 table.
 *
 * - at any time in the code, to evaluate the source functions at a
 given conformal time (by interpolating within the interpolation
 table).
 *
 * Hence the following functions can be called from other modules:
 *
 * -# perturb_init() at the beginning (but after background_init() and thermodynamics_init())  
 * -# perturb_sources_at_eta() at any later time
 * -# perturb_free() at the end, when no more calls to perturb_sources_at_eta() are needed
 */

#include "perturbations.h"
#ifdef _OPENMP
#include "omp.h"
#endif

/** 
 * Source function \f$ S^{X} (k, \eta) \f$ at given conformal time eta.
 *
 * Evaluate source function at given conformal time eta by reading
 * the pre-computed table ant interpolating.  This function can be
 * called from whatever module at whatever time, provided that
 * perturb_init() has been called before, and perturb_free() has not
 * been called yet.
 *
 * @param index_mode Input: index of requested mode
 * @param index_ic Input: index of requested initial condition
 * @param index_k Input: index of requested wavenumber
 * @param index_type Input: index of requested type
 * @param eta Input: any value of conformal time
 * @param psource_local Output: vector (assumed to be already allocated) of source functions
 * @return the error status
 */
int perturb_sources_at_eta(
			   struct perturbs * ppt,
			   int index_mode,
			   int index_ic,
			   int index_k,
			   int index_type,
			   double eta,
			   double * psource_local
			   ) {

  /** Summary: */

  /** - interpolate in pre-computed table contained in ppt using array_interpolate_two() */
  class_call(array_interpolate_two(&(ppt->sources[index_mode]
				     [index_ic * ppt->tp_size + index_type]
				     [index_k * ppt->eta_size]),
				   1,
				   0,
				   ppt->eta_sampling,
				   1,
				   ppt->eta_size,
				   eta,
				   psource_local,
				   1,
				   ppt->error_message),
	     ppt->error_message,
	     ppt->error_message);

  return _SUCCESS_;
}

/** 
 * Initialize the perturbs structure, including source functions interpolation table.
 * 
 * Initialize all fields in the structure perturbs, in particular:
 *
 * - given the values of the flags describing which kind of perturbations should be considered (modes: scalar/vector/tensor, initial conditions, type of source functions needed...), initialize indices and lists in the perturbs structure using perturb_indices_of_perturbs()
 *
 * - define the time sampling for the output source functions using perturb_timesampling_for_sources()
 *
 * - for each mode (scalar/vector/tensor): (1) initialize the indices of perturbation vectors with perturb_indices_of_current_vectors(); (2) integrate the perturbations and compute the source functions for each initial condition and wavenumber using perturb_solve()
 *
 * This function shall be called at the beginning of each run, but
 * only after background_init() and thermodynamics_init(). It
 * allocates memory spaces which should be freed later with
 * perturb_free().
 *
 * @param pba_input Input : Initialized background structure
 * @param pth_input Input : Initialized thermodynamics structure
 * @param ppr_input Input : Parameters describing how the computation is to be performed
 * @param ppt_output Output : Initialized perturbation structure
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

#pragma omp parallel 
  {
    number_of_threads = omp_get_num_threads();
  }
  class_alloc(pppw,number_of_threads * sizeof(struct perturb_workspace *),ppt->error_message);
  
#endif

  class_test(ppt->has_vectors == _TRUE_,
	     ppt->error_message,
	     "Vectors not coded yet");

  if (ppt->has_tensors == _TRUE_) printf("!!! TENSORS UNDER DEBUGGING PHASE !!!\n");

  if (ppt->has_bi == _TRUE_ || ppt->has_cdi == _TRUE_ || ppt->has_nid == _TRUE_ || ppt->has_niv == _TRUE_) {
    printf("Warning: isocurvature initial condition implemented, but not tested yet\n");
  }

  /** - decide which types of sources must be computed */
  if (ppt->has_cl_cmb_temperature == _TRUE_)
    ppt->has_source_t=_TRUE_;
  else 
    ppt->has_source_t= _FALSE_;
  
  if (ppt->has_cl_cmb_polarization == _TRUE_)
    ppt->has_source_p=_TRUE_;
  else
    ppt->has_source_p=_FALSE_;
  
  if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) ||
      (ppt->has_pk_matter == _TRUE_))
    ppt->has_source_g=_TRUE_;
  else
    ppt->has_source_g=_FALSE_;

  if (((ppt->has_source_t == _FALSE_) && (ppt->has_source_p == _FALSE_)) && (ppt->has_source_g == _FALSE_)) {
    if (ppt->perturbations_verbose > 0)
      printf("No sources requested. Perturbation module skipped.\n");
    return _SUCCESS_;
  }
  else {
    if (ppt->perturbations_verbose > 0)
      printf("Computing sources\n");
  }

  /** - initialize all indices and lists in perturbs structure using perturb_indices_of_perturbs() */
  class_call(perturb_indices_of_perturbs(ppr,
					 pba,
					 pth,
					 ppt),
	     ppt->error_message,
	     ppt->error_message);

  /** - define time sampling for sources using perturb_timesampling_for_sources() */
  class_call(perturb_timesampling_for_sources(ppr,
					      pba,
					      pth,
					      ppt),
	     ppt->error_message,
	     ppt->error_message);

  /** - loop over modes (scalar, tensors, etc). For each mode: */

  for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {

    abort = _FALSE_;

#pragma omp parallel				\
  shared(pppw,ppr,pba,pth,ppt,index_mode,abort)	\
  private(ppw)

    {

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
  shared(pppw,ppt,abort)			\
  private(ppw)

    {

#ifdef _OPENMP
      ppw=pppw[omp_get_thread_num()];
#endif
      
      class_call_parallel(perturb_workspace_free(ppt,ppw),
			  ppt->error_message,
			  ppt->error_message);

    } /* end of parallel region */


    if (abort == _TRUE_) return _FAILURE_;

  } /* end loop over modes */    

  return _SUCCESS_;
}

/**
 * Free all memory space allocated by perturb_init().
 * 
 * To be called at the end of each run, only when no further calls to
 * perturb_sources_at_eta() are needed.
 *
 * @return the error status
 */
int perturb_free(
		 struct perturbs * ppt
		 ) {

  int index_mode,index_ic,index_k,index_type;

  if (((ppt->has_source_t == _TRUE_) || (ppt->has_source_p == _TRUE_)) || (ppt->has_source_g == _TRUE_)) {

    for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {
    
      for (index_ic = 0; index_ic < ppt->ic_size[index_mode]; index_ic++) {

	for (index_type = 0; index_type < ppt->tp_size; index_type++) {

	  free(ppt->sources[index_mode][index_ic*ppt->tp_size+index_type]);

	}

      }

      free(ppt->k[index_mode]);

      free(ppt->sources[index_mode]);

    }
    
    free(ppt->eta_sampling);
	 
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
 * Initialize all indices in perturbs structure, which represent all
 * necessary indices for defining and reading the table of source
 * functions.  Allocate most of the arrays in the perturbs structure.
 * Fill the array of wavenumbers using perturb_get_k_list_size() and
 * perturb_get_k_list().
 *
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

  int index_type, index_mode, index_ic;
  int k_list_size,k_list_cl_size;

  /** - count types (eta, temperature, polarization, lensing, ...) and assign corresponding indices */
  
  index_type = 0;

  if (ppt->has_source_t == _TRUE_) {
    ppt->index_tp_t = index_type; 
    index_type++;
  }
  if (ppt->has_source_p == _TRUE_) {
    ppt->index_tp_p = index_type; 
    index_type++;
  }
  if (ppt->has_source_g == _TRUE_) {
    ppt->index_tp_g = index_type; 
    index_type++;
  }

  ppt->tp_size = index_type;

  class_test(index_type == 0,
	     ppt->error_message,
	     "you should have at least one out of {temperature, polarisation, lensing...} !!!");
  
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

  /** - allocate array of number of initial conditions for each mode, ppt->ic_size[index_mode] */

  class_alloc(ppt->ic_size,ppt->md_size*sizeof(int),ppt->error_message);

  /** - allocate array of number of wavenumbers for each mode, ppt->k_size[index_mode] */

  class_alloc(ppt->k_size,ppt->md_size * sizeof(int),ppt->error_message);
  class_alloc(ppt->k_size_cl,ppt->md_size * sizeof(int),ppt->error_message);

  /** - allocate array of lists of wavenumbers for each mode, ppt->k[index_mode] */

  class_alloc(ppt->k,ppt->md_size * sizeof(double *),ppt->error_message);

  /** - allocate array of arrays of source functions for each mode, ppt->source[index_mode] */

  class_alloc(ppt->sources,ppt->md_size * sizeof(double *),ppt->error_message);

  /** - loop over modes. For each mode: */

  for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {

    /** (a) count initial conditions (for scalars: ad, cdi, nid, niv; for tensors: only one) and assign corresponding indices */

    index_ic = 0;

    if ((ppt->has_scalars) && (index_mode == ppt->index_md_scalars)) {

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

    }
    
    if ((ppt->has_vectors) && (index_mode == ppt->index_md_vectors)) {

      /* vectors not treated yet */
    }

    if ((ppt->has_tensors) && (index_mode == ppt->index_md_tensors)) {

      ppt->index_ic_ten = index_ic;
      index_ic++;
      ppt->ic_size[index_mode] = index_ic;

    }

    /** (b) for each mode, count values of k with perturb_get_k_list_size() */

    class_call(perturb_get_k_list(ppr,
				  pba,
				  pth,
				  ppt,
				  index_mode),
	       ppt->error_message,
	       ppt->error_message);

    /** (e) allocate array of arrays of source functions for each initial conditions and wavenumber, (ppt->source[index_mode])[index_ic][index_type] */
    class_alloc(ppt->sources[index_mode],
		ppt->ic_size[index_mode] * ppt->tp_size * sizeof(double *),
		ppt->error_message);

  }
  
  return _SUCCESS_;

}

/**
 * Define time sampling for sources. 
 *
 * For each type, compute the list of values of eta at which sources
 * will be sampled.  Knowing the number of eta values, allocate all
 * arrays of source functions.
 *
 * Calls background_eta_of_z().
 *
 * Called by perturb_init().
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

  int index_mode, index_type, index_ic;

  /* time and time scale */
  double eta, timescale_source, rate_thermo, rate_isw_squared;

  /* intermediate background quantities */
  double a_prime_over_a,a_primeprime_over_a;

  /* just for calling timescale_and_approximations */
  double timescale;

  double * pvecback;
  double * pvecthermo;

  int last_index_back;
  int last_index_thermo;

  double eta_visibility_start_sources;

  class_alloc(pvecback,pba->bg_size_short*sizeof(double),ppt->error_message);  
  class_alloc(pvecthermo,pth->th_size*sizeof(double),ppt->error_message);

  /** (a) compute conformal time corresponding to opaque universe (starting point for source sampling) using background_eta_of_z() */
  class_call(background_eta_of_z(pba,pth->z_visibility_start_sources,&eta_visibility_start_sources),
	     pba->error_message,
	     ppt->error_message);

  /** (c) first, just count the number of sampling points in order to allocate the array containing all values: */

  /** (c.a) if CMB requested, first sampling point = when the universe stops being opaque; otherwise,
      start sampling gravitational potential at recombination */
  if ((ppt->has_source_t == _TRUE_) || (ppt->has_source_p == _TRUE_)) {
    eta = eta_visibility_start_sources;
  }
  else {
    eta = pth->eta_rec;
  }

  counter = 1;

  /** (c.b) next sampling point = previous + ppr->perturb_sampling_stepsize * timescale_source, where:
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

    if ((ppt->has_source_t == _TRUE_) || (ppt->has_source_p == _TRUE_)) {

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

  /** (e) infer total number of time steps, ppt->eta_size */
  ppt->eta_size = counter;

  /** (f) allocate array of time steps, ppt->eta_sampling[index_eta] */
  class_alloc(ppt->eta_sampling,ppt->eta_size * sizeof(double),ppt->error_message);

  /** (g) repeat the loop, now filling the array with each eta value: */

  /** (g.a) first sampling point = when the universe stops being opaque */
  if ((ppt->has_source_t == _TRUE_) || (ppt->has_source_p == _TRUE_)) {
    eta = eta_visibility_start_sources;
  }
  else {
    eta = pth->eta_rec;
  }
  counter = 0;
  ppt->eta_sampling[counter]=eta;

  /** (g.b) next sampling point = previous + ppr->perturb_sampling_stepsize * timescale_source, where
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

    if ((ppt->has_source_t == _TRUE_) || (ppt->has_source_p == _TRUE_)) {

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

  /** (g.c) last sampling point = exactly today */
  ppt->eta_sampling[counter] = pba->conformal_age;

  free(pvecback);
  free(pvecthermo);

  /** - Loop over modes, initial conditions and types. For each of them, allocate array of source functions, ((ppt->source[index_mode])[index_ic][index_type])[index_k][index_eta] */
  
  for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {
    for (index_ic = 0; index_ic < ppt->ic_size[index_mode]; index_ic++) {
      for (index_type = 0; index_type < ppt->tp_size; index_type++) {

	class_alloc(ppt->sources[index_mode][index_ic*ppt->tp_size+index_type],
		    ppt->k_size[index_mode] * ppt->eta_size * sizeof(double),
		    ppt->error_message);

      }
    }
  }

  return _SUCCESS_;
}

/**
 * Define the number of comoving wavenumbers using the information passed in the precision structure.
 *
 * @param index_mode Input: index describing the mode (scalar, tensor, etc.)
 * @param pk_list_size Output: number of wavenumbers 
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
  double k,k_next,k_rec,step;

  /** Summary: */

  /** - get number of wavenumbers for scalar mode */
  if ((ppt->has_scalars) && (index_mode == ppt->index_md_scalars)) {

    class_test(ppr->k_scalar_step_transition == 0.,
	       ppt->error_message,
	       "stop to avoid division by zero");
    
    class_test(pth->rs_rec == 0.,
	       ppt->error_message,
	       "stop to avoid division by zero");

    k_rec = 2. * _PI_ / pth->rs_rec;

    index_k=0;
    k = ppr->k_scalar_min * pba->H0;
    index_k=1;
    while (k < ppr->k_scalar_oscillations*k_rec) {
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

    if (k < ppt->k_scalar_kmax_for_pk*pba->h) {

      index_k += (int)((log(ppt->k_scalar_kmax_for_pk*pba->h/k)/log(10.))*ppr->k_scalar_k_per_decade_for_pk)+1;

    }

    ppt->k_size[index_mode] = index_k;

    class_alloc(ppt->k[index_mode],ppt->k_size[index_mode]*sizeof(double),ppt->error_message);

    /* now fill the array */

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
  if ((ppt->has_tensors) && (index_mode == ppt->index_md_tensors)) {
    ppt->k_size[index_mode] = ppr->k_tensor_number;
    ppt->k_size_cl[index_mode] = index_k;

    class_alloc(ppt->k[index_mode],ppt->k_size[index_mode]*sizeof(double),ppt->error_message);

    for (index_k = 0; index_k < ppr->k_tensor_number; index_k++) {
      ppt->k[index_mode][index_k] = ppr->k_tensor_min * exp(index_k * log(ppr->k_tensor_logstep));
    }

  }

  /* vectors not coded yet */


  return _SUCCESS_;

}

/**
 * Initialize all indices in vectors of perturbations for the current mode (scalar/vector/tensor).
 *
 * For the current mode (scalar/vector/tensor) with index
 * current_index_mode, this function initializes the list of indices
 * that will be used in the vectors of perturbation variables:
 * perturbations to be integrated, additional metric perturbations and
 * source terms.
 *
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

  /* - define local variables */
  int index_pt;
  int index_mt;
  int index_st;
  int index_type;
  int number_of_sources;

  /** - for scalar mode: */

  if ((ppt->has_scalars) && (index_mode == ppt->index_md_scalars)) {

    /** (a) count and assign values to indices in the vector of perturbations to be integrated */

    index_pt = 0;

    /* reject inconsistent values of the number of mutipoles in photon temperature hierachy */
    class_test(ppr->l_max_g < 4,
	       ppt->error_message,
	       "ppr->l_max_g should be at least 4, i.e. we must integrate at least over photon density, velocity, shear, third and fourth momentum");

    /* reject inconsistent values of the number of mutipoles in photon polarization hierachy */
    class_test(ppr->l_max_pol_g < 4,
	       ppt->error_message,
	       "ppr->l_max_pol_g should be at least 4");

    /* photons */

    ppw->index_pt_delta_g = index_pt; /* photon density */
    index_pt++;

    ppw->index_pt_theta_g = index_pt; /* photon velocity */
    index_pt++;

    ppw->index_pt_shear_g = index_pt; /* photon shear */
    index_pt++;

    ppw->index_pt_l3_g = index_pt; /* photon l=3 */
    index_pt++;

    ppw->l_max_g = ppr->l_max_g; /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3) */
    index_pt += (ppw->l_max_g-3);

    ppw->index_pt_pol0_g = index_pt; /* photon polarization, l=0 */
    index_pt++;

    ppw->index_pt_pol1_g = index_pt; /* photon polarization, l=1 */
    index_pt++;

    ppw->index_pt_pol2_g = index_pt; /* photon polarization, l=2 */
    index_pt++;

    ppw->index_pt_pol3_g = index_pt; /* photon polarization, l=3 */
    index_pt++;

    ppw->l_max_pol_g = ppr->l_max_pol_g; /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3) */
    index_pt += (ppw->l_max_pol_g-3); 

    /* baryons */

    ppw->index_pt_delta_b = index_pt;  /* baryon density */
    index_pt++;
    
    ppw->index_pt_theta_b = index_pt;  /* baryon velocity */
    index_pt++;

    /* cdm */

    if (pba->has_cdm == _TRUE_) {       

      ppw->index_pt_delta_cdm = index_pt; /* cdm density */
      index_pt++;

      if (ppr->gauge == newtonian) {
	ppw->index_pt_theta_cdm = index_pt; /* cdm velocity */
	index_pt++;
      }
 
    }

    /* dark energy */    
    if (pba->has_dark_energy_fluid == _TRUE_) {       
      
      ppw->index_pt_delta_de = index_pt; /* dark energy density */
      index_pt++;

      ppw->index_pt_theta_de = index_pt; /* dark energy velocity */
      index_pt++;
      
    }
    
    /* ultra relativistic neutrinos */
    if (pba->has_nur == _TRUE_) {

      /* reject inconsistent values of the number of mutipoles in ultra relativistic neutrino hierachy */
      class_test(ppr->l_max_nur < 4,
		 ppt->error_message,
		 "ppr->l_max_nur should be at least 4, i.e. we must integrate at least over neutrino/relic density, velocity, shear, third and fourth momentum");

      ppw->index_pt_delta_nur = index_pt; /* density of ultra-relativistic neutrinos/relics */
      index_pt++;

      ppw->index_pt_theta_nur = index_pt; /* velocity of ultra-relativistic neutrinos/relics */
      index_pt++;

      ppw->index_pt_shear_nur = index_pt; /* shear of ultra-relativistic neutrinos/relics */
      index_pt++;

      ppw->index_pt_l3_nur = index_pt; /* l=3 of ultra-relativistic neutrinos/relics */
      index_pt++;

      ppw->l_max_nur = ppr->l_max_nur; /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3) */
      index_pt += (ppw->l_max_nur-3);

    }

    if (ppr->gauge == synchronous) {
      ppw->index_pt_eta = index_pt; /* metric perturbation eta of synchronous gauge */
      index_pt++;
    }

    /** (b) count and assign values to indices of vector of additional metric perturbations */

    index_mt = 0;

    /* longitudinal/newtonian gauge */
    if (ppr->gauge == newtonian) {
      ppw->index_mt_phi = index_mt; /* phi */
      index_mt++;
      ppw->index_mt_psi = index_mt; /* psi */
      index_mt++;
      ppw->index_mt_phi_prime = index_mt; /* phi' */
      index_mt++;
    }
      
    /* synchronous gauge */
    if (ppr->gauge == synchronous) {
      ppw->index_mt_h_prime = index_mt; /* h */
      index_mt++;
      ppw->index_mt_eta_prime = index_mt; /* eta' */
      index_mt++;

      /* computing also alpha' (with alpha = (h' + 6 eta') / (2 k**2) ) is an option */
      ppw->index_mt_alpha_prime = index_mt; /* alpha' */
      index_mt++;

    }     

  }

  /** - for tensor mode: */

  if ((ppt->has_tensors) && (index_mode == ppt->index_md_tensors)) {


    /** (a) count and assign values to indices in the vector of perturbations to be integrated */

    index_pt = 0;

    /* reject inconsistent values of the number of mutipoles in photon temperature hierachy */
    class_test(ppr->l_max_g_ten < 5,
	       ppt->error_message,
	       "ppr->l_max_g_ten should be at least 5, i.e. we must integrate at least over photon density, velocity, shear, third, fourth, fifth momentum");
 
    /* reject inconsistent values of the number of mutipoles in photon polarization hierachy */
    class_test(ppr->l_max_pol_g_ten < 5,
	       ppt->error_message,
	       "ppr->l_max_pol_g_ten should be at least 5");

    ppw->index_pt_delta_g = index_pt; /* photon density */
    index_pt++;

    ppw->index_pt_theta_g = index_pt; /* photon velocity */
    index_pt++;

    ppw->index_pt_shear_g = index_pt; /* photon shear */
    index_pt++;

    ppw->index_pt_l3_g = index_pt; /* photon l=3 */
    index_pt++;

    ppw->l_max_g = ppr->l_max_g_ten; /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2) */
    index_pt += (ppw->l_max_g-3);
      
    ppw->index_pt_pol0_g = index_pt; /* photon polarization, l=0 */
    index_pt++;
      
    ppw->index_pt_pol1_g = index_pt; /* photon polarization, l=1 */
    index_pt++;

    ppw->index_pt_pol2_g = index_pt; /* photon polarization, l=2 */
    index_pt++;

    ppw->index_pt_pol3_g = index_pt; /* photon polarization, l=3 */
    index_pt++;

    ppw->l_max_pol_g = ppr->l_max_pol_g_ten; /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2) */
    index_pt += (ppw->l_max_pol_g-3); 

    /** (b) metric perturbation h is a propagating degree of freedom, so h and hdot are included
	in the vector of ordinary perturbations, no in that of metric perturbations */

    index_mt = 0;

    ppw->index_pt_gw = index_pt;     /* tensor metric perturbation h (gravitational waves) */
    index_pt++;

    ppw->index_pt_gwdot = index_pt; /* its time-derivative */
    index_pt++;

  }

  /** - store total size of each of these two vectors */

  ppw->pt_size = index_pt;
  ppw->mt_size = index_mt;

  /** - count and assign values to indices in the vector of source terms */

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

  class_alloc(ppw->pvecperturbations,ppw->pt_size * sizeof(double),ppt->error_message);
  class_alloc(ppw->pvecderivs,ppw->pt_size * sizeof(double),ppt->error_message);
  class_alloc(ppw->pvecback,pba->bg_size_short*sizeof(double),ppt->error_message);
  class_alloc(ppw->pvecthermo,pth->th_size*sizeof(double),ppt->error_message);
  class_alloc(ppw->pvecmetric,ppw->mt_size*sizeof(double),ppt->error_message);
  class_alloc(ppw->pvecsource_terms,ppw->st_size * ppt->tp_size * sizeof(double),ppt->error_message);
  class_alloc(ppw->last_index_back,sizeof(int),ppt->error_message);
  class_alloc(ppw->last_index_thermo,sizeof(int),ppt->error_message);

  /** - allocate source terms array for the current mode, initial condition and wavenumber: (source_terms_table[index_type])[index_eta][index_st] */
  class_alloc(ppw->source_term_table,
	      ppt->tp_size * sizeof(double *),
	      ppt->error_message);

  for (index_type = 0; index_type < ppt->tp_size; index_type++) {
    class_alloc(ppw->source_term_table[index_type], 
		ppt->eta_size*ppw->st_size*sizeof(double),
		ppt->error_message);
  }

  return _SUCCESS_;
}

int perturb_workspace_free (
			    struct perturbs * ppt,
			    struct perturb_workspace * ppw
			    ) {

  int index_type;

  free(ppw->pvecperturbations);
  free(ppw->pvecderivs);
  free(ppw->pvecback);
  free(ppw->pvecthermo);
  free(ppw->pvecmetric);
  free(ppw->pvecsource_terms);
  free(ppw->last_index_back);
  free(ppw->last_index_thermo);
  for (index_type = 0; index_type < ppt->tp_size; index_type++) {
    free(ppw->source_term_table[index_type]);
  }
  free(ppw->source_term_table);
  free(ppw);

  return _SUCCESS_;
}


/**
 * Solve the perturbation evolution for a given mode, initial condition and wavenumber, and compute the corresponding source functions.
 *
 * For a given mode, initial condition and wavenumber, this function
 * initializes all perturbations using perturb_initial_conditions(),
 * and integrates them over time. Whenever a "source sampling time" is
 * passed, the source terms are computed and stored temporarily (for
 * each type) using perturb_source_terms(). Finally, the actual
 * source functions are computed using the source terms, and stored
 * in the source table using perturb_sources().
 *
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

  /* contains all fixed parameters which should be passed to thermodynamics_derivs_with_recfast */
  struct perturb_parameters_and_workspace ppaw;

  /* conformal time */
  double eta;

  /* maximum value of conformal time for current wavenumber */
  double etamax;

  /* size of each step of integration */
  double timestep;

  /* smallest relevant timescale in the differential system (for computing step size) */
  double timescale;
  
  /* running index for the source term vector */
  int index_st;

  /* running index over types (temperature, etc) */
  int index_type;
  
  /* next index in the list of discrete eta values for which sources must be computed (for each type, temperature, polarization, lensing, etc) */
  int next_index_eta;

  /* fourier mode */
  double k;

  /* for testing purpose only */
  enum tca_flags old_tca;
  enum rp_flags old_rp;

  double eta_visibility_free_streaming;

  old_tca = tca_on;
  old_rp = rp_on;

  *(ppw->last_index_back)=0;
  *(ppw->last_index_thermo)=0;

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
  /* eventually stop earlier, when k*eta=k_eta_max, but not before the end of recombination */
  if (ppt->has_source_g == _FALSE_) {
    if ((ppr->k_eta_max/k < pba->conformal_age) && (ppr->k_eta_max/k > eta_visibility_free_streaming))
      etamax= ppr->k_eta_max/k;
    if ((ppr->k_eta_max/k < eta_visibility_free_streaming) && (eta_visibility_free_streaming < pba->conformal_age))
      etamax = eta_visibility_free_streaming;
  }

  /** - initialize generic integrator with initialize_generic_integrator() */ 

  /* Size of vector to integrate is ppw->pt_size (the number of dynamical variables in the differential system of perturbations */
  class_call(initialize_generic_integrator(ppw->pt_size, &gi),
	     gi.error_message,
	     ppt->error_message);
  
  /** - deal with the first point (initial value of conformal time eta): */

  /** (a) compute initial conformal time using information in the precision structure, and check whether it is smaller than the first sampling time */
  eta = ppr->k_eta_min / k;
  if (eta > ppt->eta_sampling[0]*ppr->eta_min_over_sampling_min)
    eta = ppt->eta_sampling[0]*ppr->eta_min_over_sampling_min;

  /** (b) compute background quantities and smallest relevant time scale in the system using perturb_timescale_and_approximations() */
  class_call(perturb_timescale_and_approximations(ppr,
						  pba,
						  pth,
						  ppt,
						  index_mode,
						  k,
						  eta,
						  normal,
						  ppw,
						  &timescale),
	     ppt->error_message,
	     ppt->error_message);

  /** (c) fill the vector of perturbation variables with appropriate initial conditions using perturb_initial_conditions() */
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

  /** (d) first time step = this time scale times the precision parameter */
  timestep = ppr->perturb_integration_stepsize*timescale;

  /** fill the structure containing fixed parameters and workspaces for perturb_derivs */
  ppaw.ppr = ppr;
  ppaw.pba = pba;
  ppaw.pth = pth;
  ppaw.ppt = ppt;
  ppaw.index_mode = index_mode;
  ppaw.k = k;
  ppaw.ppw = ppw;

  /** (f) if this is the first sampling point, compute source terms with pvecsource_terms(), store it and define next point as eta_sampling[1]; otherwise, just define next point as eta_sampling[0] */
  if (eta == ppt->eta_sampling[0]) {

    next_index_eta = 0;

    /* compute source terms at eta using perturb_source_terms() */
    class_call(perturb_source_terms(eta,&ppaw),
	       ppt->error_message,
	       ppt->error_message);
    
    /* store source terms in the array (ppw->source_term_table[index_type]) for each type */

    for (index_type = 0; index_type < ppt->tp_size; index_type++) {
      for (index_st = 0; index_st < ppw->st_size; index_st++) {
	ppw->source_term_table[index_type][next_index_eta*ppw->st_size+index_st] = 
	  ppw->pvecsource_terms[index_type * ppw->st_size + index_st];
      }
    }    
    
    next_index_eta = 1;
  }
  else {
    next_index_eta = 0;
  }

  /** - start loop over sampling values eta_sampling[next_index_eta]: */

  while (next_index_eta < ppt->eta_size) {

    if (ppt->eta_sampling[next_index_eta] <= etamax) {

      /** (a) perform intermediate integration step before reaching the next eta_sampling[next_index_eta] value: */

      /* counter=0; */

      while (eta + 2.*timestep < ppt->eta_sampling[next_index_eta]) {

	/** (a.1) integrate perturbations till eta + timestep using generic_integrator() */
	class_call(generic_integrator(perturb_derivs,
				      eta,
				      eta+timestep,
				      ppw->pvecperturbations,
				      &ppaw,
				      ppr->tol_perturb_integration,
				      ppr->smallest_allowed_variation,
				      &gi),
		   gi.error_message,
		   ppt->error_message);
      
	/** (a.2) define new time value eta = eta + timestep */
	eta = eta + timestep;

	class_test(timestep < ppr->smallest_allowed_variation,
		   ppt->error_message,
		   "integration step =%e < machine precision : leads either to numerical error or infinite loop",timestep);

	/** (a.3) compute background quantities and smallest relevant time scale in the system using perturb_timescale_and_approximations() */
	class_call(perturb_timescale_and_approximations(ppr,
							pba,
							pth,
							ppt,
							index_mode,
							k,
							eta,
							closeby,
							ppw,
							&timescale),
		   ppt->error_message,
		   ppt->error_message);

	/* for testing purposes */
	if (ppt->perturbations_verbose > 3) {

	  class_test(ppw->pvecthermo[pth->index_th_dkappa] == 0.,
		     ppt->error_message,
		     "you have dkappak=0, stop to avoid division by zero");

	  class_test(ppw->pvecback[pba->index_bg_H]*ppw->pvecback[pba->index_bg_a] == 0.,
		     ppt->error_message,
		     "you have aH=0, stop to avoid division by zero");

	  if (old_tca==tca_on && ppw->tca==tca_off)
	    printf("Turn off tight-coupling at eta=%e, with k*eta_c=%e and eta_c/eta_h=%e\n",
		   eta,
		   k/ppw->pvecthermo[pth->index_th_dkappa],
		   (ppw->pvecback[pba->index_bg_H]*ppw->pvecback[pba->index_bg_a])/ppw->pvecthermo[pth->index_th_dkappa]);

	  if (old_tca==tca_off && ppw->tca==tca_on)
	    printf("Turn on tight-coupling  at eta=%e, with k*eta_c=%e\n",
		   eta,
		   k/ppw->pvecthermo[pth->index_th_dkappa]);

	  if (old_rp==rp_on && ppw->rp==rp_off)
	    printf("Turn on free-streaming  at eta=%e, with k/aH   =%e and Omega_r    =%e\n",eta,k/(ppw->pvecback[pba->index_bg_H]*ppw->pvecback[pba->index_bg_a]),ppw->pvecback[pba->index_bg_Omega_r]);

	  if (old_rp==rp_off && ppw->rp==rp_on)
	    printf("Turn off free-streaming at eta=%e\n",eta);

	  old_tca=ppw->tca;

	  old_rp=ppw->rp;

	} /* end of testing zone */
      

	/** (a.4) new time step = this time scale times the precision parameter */
	timestep = ppr->perturb_integration_stepsize*timescale;

      }

      /** (b) perform integration step till exactly the next eta_sampling[next_index_eta] value: */    

      /** (b.1) integrate perturbations over current step using generic_integrator() */
      class_call(generic_integrator(perturb_derivs,
				    eta,
				    ppt->eta_sampling[next_index_eta],
				    ppw->pvecperturbations,
				    &ppaw,
				    ppr->tol_perturb_integration,
				    ppr->smallest_allowed_variation,
				    &gi),
		 gi.error_message,
		 ppt->error_message);

      /** (a.2) define new time value eta = eta_sampling[next_index_eta] */
      eta = ppt->eta_sampling[next_index_eta];

      /** (a.3) compute background quantities and smallest relevant time scale in the system using perturb_timescale_and_approximations() */
      class_call(perturb_timescale_and_approximations(ppr,
						      pba,
						      pth,
						      ppt,
						      index_mode,
						      k,
						      eta,
						      closeby,
						      ppw,
						      &timescale),
		 ppt->error_message,
		 ppt->error_message);

      /* for testing purpose */     
      if (ppt->perturbations_verbose > 2) {

	class_test(ppw->pvecthermo[pth->index_th_dkappa] == 0.,
		   ppt->error_message,
		   "you have dkappak=0, stop to avoid division by zero");
	
	class_test(ppw->pvecback[pba->index_bg_H]*ppw->pvecback[pba->index_bg_a] == 0.,
		   ppt->error_message,
		   "you have aH=0, stop to avoid division by zero");

	if (old_tca==tca_on && ppw->tca==tca_off)
	  printf("Turn off tight-coupling at eta=%e, with k*eta_c=%e and eta_c/eta_h=%e\n",
		 eta,
		 k/ppw->pvecthermo[pth->index_th_dkappa],
		 (ppw->pvecback[pba->index_bg_H]*ppw->pvecback[pba->index_bg_a])/ppw->pvecthermo[pth->index_th_dkappa]);

	if (old_tca==tca_off && ppw->tca==tca_on)
	  printf("Turn on tight-coupling  at eta=%e, with k*eta_c=%e\n",
		 eta,
		 k/ppw->pvecthermo[pth->index_th_dkappa]);

	if (old_rp==rp_on && ppw->rp==rp_off)
	  printf("Turn on free-streaming  at eta=%e, with k/aH   =%e and Omega_r    =%e\n",eta,k/(ppw->pvecback[pba->index_bg_H]*ppw->pvecback[pba->index_bg_a]),ppw->pvecback[pba->index_bg_Omega_r]);

	if (old_rp==rp_off && ppw->rp==rp_on)
	  printf("Turn off free-streaming at eta=%e\n",eta);

	old_tca=ppw->tca;

	old_rp=ppw->rp;

      } /* end of testing zone */

      /** (a.4) new time step = this time scale times the precision parameter */
      timestep = ppr->perturb_integration_stepsize*timescale;

      class_test(timestep < ppr->smallest_allowed_variation,
		 ppt->error_message,
		 "integration step =%e < machine precision : leads either to numerical error or infinite loop",timestep);

      /** (a.6) compute source terms at eta using perturb_source_terms() */
      class_call(perturb_source_terms(eta,&ppaw),
		 ppt->error_message,
		 ppt->error_message);

      /** (a.7) store source terms in the array (ppw->source_term_table[index_type]) for each type */

      for (index_type = 0; index_type < ppt->tp_size; index_type++) {
	for (index_st = 0; index_st < ppw->st_size; index_st++) {
	  ppw->source_term_table[index_type][next_index_eta*ppw->st_size+index_st] = 
	    ppw->pvecsource_terms[index_type * ppw->st_size + index_st];
	}
      }    
    }
    else {
      for (index_type = 0; index_type < ppt->tp_size; index_type++) {
	for (index_st = 0; index_st < ppw->st_size; index_st++) {
	  ppw->source_term_table[index_type][next_index_eta*ppw->st_size+index_st] = 0.;
	}
      }
    }


    /** (a.7) increment index next_index_eta of eta_sampling[next_index_eta] */ 
    next_index_eta++;

  } /* end of loop over time */

  /** - infer source functions from source terms using perturb_sources() */
      
  class_call(perturb_sources(ppt,
			     index_mode,
			     index_ic,
			     index_k,
			     ppw),
	     ppt->error_message,
	     ppt->error_message);
  
  /** - clean up generic integrator with cleanup_generic_integrator() */

  class_call(cleanup_generic_integrator(&gi),
	     gi.error_message,
	     ppt->error_message);

  return _SUCCESS_;
}
  
/**
 * For each mode, wavenumber and initial condition, this function initializes all values
 * in the vector of perturbed variables. The mode, wavenumber and initial condition information is passed through the global variables index_mode, current_index_k, current_index_ic.
 *
 * @param eta Input: conformal time
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

  ppw->pvecperturbations[ppw->index_pt_delta_g] = 0.; /* photon density */
  ppw->pvecperturbations[ppw->index_pt_theta_g] = 0.; /* photon velocity */
  ppw->pvecperturbations[ppw->index_pt_shear_g] = 0.; /* photon shear */
  ppw->pvecperturbations[ppw->index_pt_l3_g] = 0.; /* photon shear */
  for (l=4; l <= ppw->l_max_g; l++) ppw->pvecperturbations[ppw->index_pt_delta_g+l] = 0.;  /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3) */

  ppw->pvecperturbations[ppw->index_pt_pol0_g] = 0.; /* photon polarization, l=0 */
  ppw->pvecperturbations[ppw->index_pt_pol1_g] = 0.; /* photon polarization, l=1 */
  ppw->pvecperturbations[ppw->index_pt_pol2_g] = 0.; /* photon polarization, l=2 */
  ppw->pvecperturbations[ppw->index_pt_pol3_g] = 0.; /* photon polarization, l=3 */
  for (l=4; l <= ppw->l_max_pol_g; l++) ppw->pvecperturbations[ppw->index_pt_pol0_g+l] = 0.;  /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3) */

  /* additional perturbations relevant only for the scalar mode */

  if ((ppt->has_scalars && index_mode) == (ppt->index_md_scalars)) {

    ppw->pvecperturbations[ppw->index_pt_delta_b] = 0.;  /* baryon density */
    ppw->pvecperturbations[ppw->index_pt_theta_b] = 0.;  /* baryon velocity */

    if (pba->has_cdm == _TRUE_) {       
      ppw->pvecperturbations[ppw->index_pt_delta_cdm] = 0.; /* cdm density */
      if (ppr->gauge == newtonian) 
	ppw->pvecperturbations[ppw->index_pt_theta_cdm] = 0.; /* cdm velocity */
    }
    
    if (pba->has_dark_energy_fluid == _TRUE_) {        
      ppw->pvecperturbations[ppw->index_pt_delta_de] = 0.; /* dark energy density */   
      ppw->pvecperturbations[ppw->index_pt_theta_de] = 0.; /* dark energy velocity */ 
    } 
    
    if (pba->has_nur == _TRUE_) {
      ppw->pvecperturbations[ppw->index_pt_delta_nur] = 0; /* density of ultra-relativistic neutrinos/relics */
      ppw->pvecperturbations[ppw->index_pt_theta_nur] = 0; /* velocity of ultra-relativistic neutrinos/relics */
      ppw->pvecperturbations[ppw->index_pt_shear_nur] = 0.; /* shear of ultra-relativistic neutrinos/relics */
      ppw->pvecperturbations[ppw->index_pt_l3_nur] = 0.; /* l=3 of ultra-relativistic neutrinos/relics */
      for (l=4; l <= ppw->l_max_nur; l++)
	ppw->pvecperturbations[ppw->index_pt_delta_nur+l] = 0.;  /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2) */
    }
    
    if (ppr->gauge == synchronous)
      ppw->pvecperturbations[ppw->index_pt_eta] = 0; /* metric perturbation eta */ 

  }

  /* additional perturbations relevant only for the tensor mode */

  if ((ppt->has_tensors) && (index_mode == ppt->index_md_tensors)) {
    ppw->pvecperturbations[ppw->index_pt_gw] = 0.;     /* tensor metric perturbation h (gravitational waves) */
    ppw->pvecperturbations[ppw->index_pt_gwdot] = 0.;  /* its time-derivative */
  }

  /** - initial conditions for scalars */

  if ((ppt->has_scalars) && (index_mode == ppt->index_md_scalars)) {

    /** (a) adiabatic */ 

    if ((ppt->has_ad) && (index_ic == ppt->index_ic_ad)) {

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

	ppw->pvecperturbations[ppw->index_pt_delta_g] = -4.*alpha/eta - k*k*eta*eta/3.; /* photon density */
	ppw->pvecperturbations[ppw->index_pt_theta_g] = alpha*k*k - pow(k*eta,3.)*k/36.; /* photon velocity */

	ppw->pvecperturbations[ppw->index_pt_delta_b] = 3./4.*ppw->pvecperturbations[ppw->index_pt_delta_g]; /* baryon density */
	ppw->pvecperturbations[ppw->index_pt_theta_b] = ppw->pvecperturbations[ppw->index_pt_theta_g]; /* baryon velocity */
      
	if (pba->has_cdm == _TRUE_) {       
	  ppw->pvecperturbations[ppw->index_pt_delta_cdm] = 3./4.*ppw->pvecperturbations[ppw->index_pt_delta_g]; /* cdm density */
	  ppw->pvecperturbations[ppw->index_pt_theta_cdm] = alpha*k*k; /* cdm velocity */
	}
	
 	if (pba->has_dark_energy_fluid == _TRUE_) {        
 	  ppw->pvecperturbations[ppw->index_pt_delta_de] = 0.; /* dark energy density (TO BE WRITTEN) */
 	  ppw->pvecperturbations[ppw->index_pt_theta_de] = 0.; /* dark energy velocity (TO BE WRITTEN) */
 	} 
	
	if (pba->has_nur == _TRUE_) {
	  ppw->pvecperturbations[ppw->index_pt_delta_nur] = ppw->pvecperturbations[ppw->index_pt_delta_g]; /* density of ultra-relativistic neutrinos/relics */
	  ppw->pvecperturbations[ppw->index_pt_theta_nur] = alpha*k*k - pow(k*eta,3.)*k/36. * (23.+4.*fracnu)/(15.+4.*fracnu); /* velocity of ultra-relativistic neutrinos/relics */
	  ppw->pvecperturbations[ppw->index_pt_shear_nur] = k*k*eta*eta*2./3./(12.+fracnu); /* shear of ultra-relativistic neutrinos/relics */
	}

      }

      /* synchronous gauge */
      if (ppr->gauge == synchronous) {

	ppw->pvecperturbations[ppw->index_pt_delta_g] = - k*k*eta*eta/3. * (1.-om*eta/5.); /* photon density */
	/* ppw->pvecperturbations[ppw->index_pt_theta_g] = - k*k*k*k*eta*eta*eta/36.; /\* photon velocity *\/ */
	ppw->pvecperturbations[ppw->index_pt_theta_g] = k*k*eta/9.*ppw->pvecperturbations[ppw->index_pt_delta_g]; /* photon velocity */

	ppw->pvecperturbations[ppw->index_pt_delta_b] = 3./4.*ppw->pvecperturbations[ppw->index_pt_delta_g]; /* baryon density */
	ppw->pvecperturbations[ppw->index_pt_theta_b] = ppw->pvecperturbations[ppw->index_pt_theta_g]; /* baryon velocity */
      
	if (pba->has_cdm == _TRUE_) {       
	  ppw->pvecperturbations[ppw->index_pt_delta_cdm] = 3./4.*ppw->pvecperturbations[ppw->index_pt_delta_g]; /* cdm density */
          /* by convention, cdm velocity velocity vanishes in this implementation of the synchronous gauge */
	}

 	if (pba->has_dark_energy_fluid == _TRUE_) {        
 	  ppw->pvecperturbations[ppw->index_pt_delta_de] = 0; /* dark energy density (TO BE WRITTEN) */
	  ppw->pvecperturbations[ppw->index_pt_theta_de] = 0; /* dark energy velocity (TO BE WRITTEN) */
 	} 

	if (pba->has_nur == _TRUE_) {
	  ppw->pvecperturbations[ppw->index_pt_delta_nur] = ppw->pvecperturbations[ppw->index_pt_delta_g]; /* density of ultra-relativistic neutrinos/relics */
	  ppw->pvecperturbations[ppw->index_pt_theta_nur] = - pow(k*eta,3.)*k/36. * (23.+4.*fracnu)/(15.+4.*fracnu); /* velocity of ultra-relativistic neutrinos/relics */
	  ppw->pvecperturbations[ppw->index_pt_shear_nur] = k*k*eta*eta*2./3./(12.+fracnu); /* shear of ultra-relativistic neutrinos/relics */
	}    

	ppw->pvecperturbations[ppw->index_pt_eta] = 1.-(5.+4.*fracnu)/12./(15.+4.*fracnu)*k*k*eta*eta; /* metric perturbation eta */

      }

    }

    /** (b) Cold dark matter Isocurvature */ 

    if ((ppt->has_cdi) && (index_ic == ppt->index_ic_cdi)) { 
      
      class_test(pba->has_cdm == _FALSE_,
		 ppt->error_message,
		 "not consistent to ask for CDI in absence of CDM!");

      ppw->pvecperturbations[ppw->index_pt_delta_cdm] = ppr->entropy_ini;

    }

    /** (c) Baryon Isocurvature */ 

    if ((ppt->has_bi) && (index_ic == ppt->index_ic_bi)) {

      ppw->pvecperturbations[ppw->index_pt_delta_b] = ppr->entropy_ini;

    }

    /** (d) Neutrino density Isocurvature */ 

    if ((ppt->has_nid) && (index_ic == ppt->index_ic_nid)) {

      class_test(pba->has_nur == _FALSE_,
		 ppt->error_message,
		 "not consistent to ask for NID in absence of neutrinos!");

      ppw->pvecperturbations[ppw->index_pt_delta_nur] = ppr->entropy_ini;

    }
     
    /** (e) Neutrino velocity Isocurvature */ 

    if ((ppt->has_niv) && (index_ic == ppt->index_ic_niv)) {

      class_test(pba->has_nur == _FALSE_,
		 ppt->error_message,
		 "not consistent to ask for NIV in absence of neutrinos!");

      ppw->pvecperturbations[ppw->index_pt_theta_nur] = ppr->entropy_ini;

    }

  }

  /** - initial conditions for tensors */

  if ((ppt->has_tensors) && (index_mode == ppt->index_md_tensors)) {

    if (index_ic == ppt->index_ic_ten) {
      ppw->pvecperturbations[ppw->index_pt_gw] = ppr->gw_ini; 
    }

  }

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
 * @param eta Input: conformal time  
 * @param intermode Input: interpolation mode (normal or growing_closeby)
 * @param last_index_back Input/Ouput: index of the previous/current point in the background interpolation array (relevant for closeby mode only) 
 * @param last_index_thermo Input/Ouput: index of the previous/current point in the thermodynamics interpolation array (relevant for closeby mode only) 
 * @param tca_local Output: flag: whether the tight-coupling approximation holds
 * @param rp_local Ouput: flag: whether radiation perturbation need to be integrated
 * @param timescale Ouput: smallest relevant time scale in the differential system of perturbations 
 * @return the error status
 */
int perturb_timescale_and_approximations(
					 struct precision * ppr,
					 struct background * pba,
					 struct thermo * pth,
					 struct perturbs * ppt,
					 int index_mode,
					 double k,
					 double eta,
					 enum interpolation_mode intermode,
					 struct perturb_workspace * ppw,
					 double * timescale
					 ) {
  /** Summary: */

  /** - define local variables */

  /** (a) time scale of Fourier mode, \f$ \eta_k = 1/k \f$ */  
  double eta_k;
  /** (b) time scale of expansion, \f$ \eta_h = a/a' \f$ */
  double eta_h;
  /** (c) time scale of recombination, \f$ \eta_{\gamma} = 1/\kappa' \f$ */
  double eta_g;

  class_test(k == 0.,
	     ppt->error_message,
	     "stop to avoid division by zero");

  /** - compute Fourier mode time scale = \f$ \eta_k = 1/k \f$ */
  eta_k = 1./k;

  /** - evaluate background quantities with background_at_eta() and Hubble time scale \f$ \eta_h = a/a' \f$ */
  class_call(background_at_eta(pba,eta, short_info, intermode, ppw->last_index_back, ppw->pvecback),
	     pba->error_message,
	     ppt->error_message);

  class_test(ppw->pvecback[pba->index_bg_H]*ppw->pvecback[pba->index_bg_a] == 0.,
	     ppt->error_message,
	     "aH=0, stop to avoid division by zero");

  eta_h = 1./(ppw->pvecback[pba->index_bg_H]*ppw->pvecback[pba->index_bg_a]);

  /** - if \f$ k \gg aH \f$ and \f$ \Omega_r \ll 1 \f$ and efficient recombination finished, switch off radiation perturbations and keep Hubble time scale as the relevant one. Otherwise take \f$ \min(\eta_k, \eta_h). \f$ */
  if ((eta_h/eta_k > ppr->rad_pert_trigger_k_over_aH) && 
      //      (eta > eta_visibility_free_streaming) && /* optionally this line could be restored, to check that this does not happen before recombination is completed) */
      (ppw->pvecback[pba->index_bg_Omega_r] < ppr->rad_pert_trigger_Omega_r)) {
    ppw->rp = rp_off;
    *timescale = eta_h;
  }
  else {
    ppw->rp = rp_on;
    if (eta_k < eta_h) 
      *timescale = eta_k;
    else  
      *timescale = eta_h;
  }

  /** - for scalars modes: */
  if ((ppt->has_scalars) && (index_mode == ppt->index_md_scalars)) {

    /** (a) evaluate thermodynamical quantities with thermodynamics_at_z() */
    class_call(thermodynamics_at_z(pba,
				   pth,
				   1./ppw->pvecback[pba->index_bg_a]-1.,  /* redshift z=1/a-1 */
				   intermode,
				   ppw->last_index_thermo,
				   ppw->pvecback,
				   ppw->pvecthermo),
	       pth->error_message,
	       ppt->error_message);
    
    /** (b.1.) if \f$ \kappa'=0 \f$, recombination is finished; check that tight-coupling approximation is off */
    if (ppw->pvecthermo[pth->index_th_dkappa] == 0.) {
      ppw->tca = tca_off;
    }

    /** (b.2.) if \f$ \kappa' \neq 0 \f$, recombination is not finished: */
    else {

      /** (b.2.a) compute recombination time scale for photons, \f$ \eta_{\gamma} = 1/ \kappa' \f$ */
      eta_g = 1./ppw->pvecthermo[pth->index_th_dkappa];

      /** (b.2.b) if \f$ \eta_g \ll max(\eta_h, \eta_k) \f$, turn on the tight-coupling approximation for photons; otherwise turn off  and define relevant time scale as min time scale between \f$ \eta_k \f$, \f$ \eta_h \f$ and \f$ \eta_g \f$ */
      if ((eta_g/eta_h < ppr->tight_coupling_trigger_eta_g_over_eta_h) 
	  && (eta_g/eta_k < ppr->tight_coupling_trigger_eta_g_over_eta_k)) {
	ppw->tca = tca_on;
      }
      else {
	ppw->tca = tca_off;
	if (eta_g < *timescale) *timescale=eta_g;
      }
    }
    
  }

  /** - for tensor modes: tight-coupling approximation is off, time scale remains \f$ min (\eta_k , \eta_h) \f$. */
  if ((ppt->has_tensors) && (index_mode == ppt->index_md_tensors)) {
    ppw->tca = tca_off;
  }

  /* vectors not coded yet */

  return _SUCCESS_;
}

/**
 * Compute metric perturbations (those not integrated over time) using Einstein equations
 *
 * At a given time \f$ \eta \f$ and given the input values of perturbations (those
 * which are integrated iver time), 
 * compute the metric perturbations using constraint equations
 * provided by Einstein equations.
 *
 * @param eta Input: conformal time  
 * @return the error status
 */
int perturb_einstein(
		     struct precision * ppr,
		     struct background * pba,
		     struct perturbs * ppt,
		     int index_mode,
		     double k,
		     double eta,       /**< Input : conformal time */
		     double * y,
		     struct perturb_workspace * ppw
		     ) {
  /** Summary: */

  /** - define local variables */

  double k2,a,a2,a_prime_over_a;
  double delta_rho,delta_theta,delta_shear;
  double alpha;
  double * pvecback;
  double * pvecmetric;

  pvecback = ppw->pvecback;
  pvecmetric = ppw->pvecmetric;

  /** - compute \f$ k^2 \f$, \f$ a \f$, \f$ a^2 \f$, \f$ a'/a \f$ */ 

  k2 = k*k;
  a = pvecback[pba->index_bg_a];
  a2 = a * a;
  a_prime_over_a = pvecback[pba->index_bg_H]*a;

  /** - for scalar modes: */  
  if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {
    
    /** (a) compute the total \f$ \delta \rho, \delta \theta, \delta \sigma \f$ */
 
    /* photon and baryon contribution */
    delta_rho = pvecback[pba->index_bg_rho_g]*y[ppw->index_pt_delta_g]
      + pvecback[pba->index_bg_rho_b]*y[ppw->index_pt_delta_b];
    delta_theta = 4./3.*pvecback[pba->index_bg_rho_g]*y[ppw->index_pt_theta_g]
      + pvecback[pba->index_bg_rho_b]*y[ppw->index_pt_theta_b];
    delta_shear = 4./3.*pvecback[pba->index_bg_rho_g]*y[ppw->index_pt_shear_g];

    /* cdm contribution */
    if (pba->has_cdm == _TRUE_) {
      delta_rho = delta_rho + pvecback[pba->index_bg_rho_cdm]*y[ppw->index_pt_delta_cdm];
      if (ppr->gauge == newtonian)
	delta_theta = delta_theta + pvecback[pba->index_bg_rho_cdm]*y[ppw->index_pt_theta_cdm];
    }
    
    /* dark energy fluid contribution */
    if (pba->has_dark_energy_fluid == _TRUE_) { 
      delta_rho = delta_rho + pvecback[pba->index_bg_rho_de]*y[ppw->index_pt_delta_de]; 
      delta_theta = delta_theta + pvecback[pba->index_bg_rho_de]*y[ppw->index_pt_theta_de];
    } 

    /* ultra-relativistic neutrino/relics contribution */
    if (pba->has_nur == _TRUE_) {
      delta_rho = delta_rho + pvecback[pba->index_bg_rho_nur]*y[ppw->index_pt_delta_nur];
      delta_theta = delta_theta + 4./3.*pvecback[pba->index_bg_rho_nur]*y[ppw->index_pt_theta_nur];
      delta_shear = delta_shear + 4./3.*pvecback[pba->index_bg_rho_nur]*y[ppw->index_pt_shear_nur];
    }
    
    /** (b) compute metric perturbations */

    /* newtonian gauge */
    if (ppr->gauge == newtonian) {
      pvecmetric[ppw->index_mt_phi] = -1.5 * (a2/k2/k2) * (k2 * delta_rho + 3.*a_prime_over_a * delta_theta); /* phi */
      pvecmetric[ppw->index_mt_psi] = pvecmetric[ppw->index_mt_phi] - 4.5 * (a2/k2) * delta_shear;  /* psi */
      pvecmetric[ppw->index_mt_phi_prime] = - a_prime_over_a * pvecmetric[ppw->index_mt_psi] + 1.5 * (a2/k2) * delta_theta; /* phi' */
    }

    /* synchronous gauge */
    if (ppr->gauge == synchronous) {
      pvecmetric[ppw->index_mt_eta_prime] = 1.5 * (a2/k2) * delta_theta;  /* eta' */
      pvecmetric[ppw->index_mt_h_prime] = 
	( k2 * y[ppw->index_pt_eta] + 1.5 * a2 * delta_rho)/(0.5*a_prime_over_a);  /* h' */
      pvecmetric[ppw->index_mt_alpha_prime] = 
	- 4.5 * (a2/k2) * delta_shear + y[ppw->index_pt_eta] - 2.*a_prime_over_a*
	(pvecmetric[ppw->index_mt_h_prime] + 6.*pvecmetric[ppw->index_mt_eta_prime])
	/ 2./ k2; /* alpha' = (h''+6eta'')/2k2 */

      /* getting phi here is an option */
      /* phi=y[ppw->index_pt_eta]-0.5 * (a_prime_over_a/k2) * (h_plus_six_eta_prime); */   /* phi from gauge transformation (from synchronous to newtonian) */

    }

    return _SUCCESS_;

  }

  if ((ppt->has_tensors == _TRUE_) || (index_mode == ppt->index_md_tensors)) {

    /* nothing to be done for tensors: only one propagating degree of freedom, no constraint equation */

    return _SUCCESS_;

  }

  class_test(0 == 0,
	     ppt->error_message,
	     "strange value of index_mode: neither scalar not tensor. Vectors not coded yet!");

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
 * @param eta Input: conformal time  
 * @return the error status
 */
int perturb_source_terms(
			 double eta,
			 struct perturb_parameters_and_workspace * pppaw
			 ) {
  /** Summary: */

  /** - define local variables */

  double k2,a_prime_over_a,a_primeprime_over_a,R;
  double Pi,Pi_prime;
  double x2;
  int index_type;

  struct precision * ppr;
  struct background * pba;
  struct thermo * pth;
  struct perturbs * ppt;
  int index_mode;
  double k;
  struct perturb_workspace * ppw;
  double * pvecback;
  double * pvecthermo;
  double * pvecmetric;
  double * pvecperturbations;
  double * pvecderivs; 
  double * pvecsource_terms;

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
  pvecperturbations = ppw->pvecperturbations;
  pvecderivs = ppw->pvecderivs;
  pvecsource_terms = ppw->pvecsource_terms;

  class_call(perturb_derivs(eta,
			    pvecperturbations,
			    pvecderivs,
			    pppaw,
			    ppt->error_message),
	     ppt->error_message,
	     ppt->error_message);

  k2 = k * k;

  a_prime_over_a = pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a];

  a_primeprime_over_a = pvecback[pba->index_bg_H_prime] * pvecback[pba->index_bg_a] 
    + 2. * a_prime_over_a * a_prime_over_a;

  R = 4./3. * pvecback[pba->index_bg_rho_g]/pvecback[pba->index_bg_rho_b];

  /** - compute \f$ k^2 \f$, \f$ \Pi = G_{\gamma 0} + G_{\gamma 2} + F_{\gamma 2} \f$, \f$ e^{- \kappa} \f$ */ 

  Pi = pvecperturbations[ppw->index_pt_pol0_g] + pvecperturbations[ppw->index_pt_pol2_g] + 2.*pvecperturbations[ppw->index_pt_shear_g];

  Pi_prime = pvecderivs[ppw->index_pt_pol0_g] + pvecderivs[ppw->index_pt_pol2_g] + 2.*pvecderivs[ppw->index_pt_shear_g];


  /** - for each type and each mode, compute S0, S1, S2 */
  for (index_type = 0; index_type < ppt->tp_size; index_type++) {

    pvecsource_terms[index_type * ppw->st_size + ppw->index_st_eta] = eta;

    /* temperature */
    if ((ppt->has_source_t == _TRUE_) && (index_type == ppt->index_tp_t)) {

      /* scalar temperature */
      if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

        /* check that visibility is non-zero (otherwise source = 0) */
	if (pvecthermo[pth->index_th_g] != 0.) {

          /* newtonian gauge */
	  if (ppr->gauge == newtonian) {

            /* S0 */
	    pvecsource_terms[index_type * ppw->st_size + ppw->index_st_S0] =
	      pvecthermo[pth->index_th_exp_m_kappa] * pvecmetric[ppw->index_mt_phi_prime]
	      + pvecthermo[pth->index_th_g] / 4. * (pvecperturbations[ppw->index_pt_delta_g] + Pi / 4.);
	    
            /* S1 */
	    pvecsource_terms[index_type * ppw->st_size + ppw->index_st_S1] =
	      pvecthermo[pth->index_th_exp_m_kappa] * pvecmetric[ppw->index_mt_psi]
	      + pvecthermo[pth->index_th_g] * pvecperturbations[ppw->index_pt_theta_b] / k2;

	    /* S2 */
	    pvecsource_terms[index_type * ppw->st_size + ppw->index_st_S2] =
	      3./16. * pvecthermo[pth->index_th_g] * Pi / k2;

	  }

          /* synchronous gauge */
	  if (ppr->gauge == synchronous) {

	    /* S0 */
	    pvecsource_terms[index_type * ppw->st_size + ppw->index_st_S0] =
	      pvecthermo[pth->index_th_exp_m_kappa] * pvecmetric[ppw->index_mt_eta_prime]
	      + pvecthermo[pth->index_th_g] / 4. * (pvecperturbations[ppw->index_pt_delta_g] + Pi / 4.);
	  
	    /* S1 */
	    /* pvecsource_terms[index_type * ppw->st_size + ppw->index_st_S1] = */
	    /* 	      pvecthermo[pth->index_th_g] * pvecperturbations[ppw->index_pt_theta_b] / k2; */

	    /* dS1 */
	    pvecsource_terms[index_type * ppw->st_size + ppw->index_st_dS1] =
	      pvecthermo[pth->index_th_dg] * pvecperturbations[ppw->index_pt_theta_b] / k2
	      + pvecthermo[pth->index_th_g] * pvecderivs[ppw->index_pt_theta_b] / k2;

	    /* S2 */
	    /* pvecsource_terms[index_type * ppw->st_size + ppw->index_st_S2] = */
	    /* 	      pvecthermo[pth->index_th_exp_m_kappa] * (pvecmetric[ppw->index_mt_h_prime] + 6. * pvecmetric[ppw->index_mt_eta_prime])/2./k2  */
	    /* 	      + 3./16. * pvecthermo[pth->index_th_g] * Pi / k2; */

            /* dS2 */
	    pvecsource_terms[index_type * ppw->st_size + ppw->index_st_dS2] =
	      pvecthermo[pth->index_th_g] 
	      * (pvecmetric[ppw->index_mt_h_prime] + 6. * pvecmetric[ppw->index_mt_eta_prime])/2./k2
	      + pvecthermo[pth->index_th_exp_m_kappa] * (pvecmetric[ppw->index_mt_alpha_prime])
	      + 3./16. * pvecthermo[pth->index_th_dg] * Pi / k2
	      + 3./16. * pvecthermo[pth->index_th_g] * Pi_prime / k2;

	    /* 	    Pi_prime = -k*pvecperturbations[ppw->index_pt_pol0_g+1] */
	    /* 	      -pvecthermo[pth->index_th_dkappa] * (pvecperturbations[ppw->index_pt_pol0_g]-Pi/2.) */
	    /* 	      +k/5. * (2.*pvecperturbations[ppw->index_pt_pol2_g-1]-3.*pvecperturbations[ppw->index_pt_pol2_g+1]) */
	    /* 	      -pvecthermo[pth->index_th_dkappa] * (pvecperturbations[ppw->index_pt_pol2_g]-Pi/10.) */
	    /*               +8./15.*pvecperturbations[ppw->index_pt_theta_g] */
	    /* 	      -3./5.*k*pvecperturbations[ppw->index_pt_shear_g+1] */
	    /* 	      +4./15.*pvecmetric[ppw->index_mt_h_prime] */
	    /* 	      +8./5.*pvecmetric[ppw->index_mt_eta_prime] */
	    /* 	      -pvecthermo[pth->index_th_dkappa] * (2.*pvecperturbations[ppw->index_pt_shear_g]-1./10.*Pi); */

	    /* 	    Pi_prime = -3./10.*pvecthermo[pth->index_th_dkappa]*Pi */
	    /* 	      -3./5.*k * (pvecperturbations[ppw->index_pt_pol1_g]+ */
	    /* 				pvecperturbations[ppw->index_pt_pol2_g+1]+ */
	    /* 				pvecperturbations[ppw->index_pt_shear_g+1]) */
	    /* 	      +8./15.*pvecperturbations[ppw->index_pt_theta_g] */
	    /* 	      +4./15. * (pvecmetric[ppw->index_mt_h_prime] + 6. * pvecmetric[ppw->index_mt_eta_prime]); */

	  }
	}

	else {

	  pvecsource_terms[index_type * ppw->st_size + ppw->index_st_S0] = 0.;
	  pvecsource_terms[index_type * ppw->st_size + ppw->index_st_dS1] = 0.;
	  pvecsource_terms[index_type * ppw->st_size + ppw->index_st_dS2] = 0.;

	}
      }

      /* temperature tensors */
      if ((ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors)) {

	/* --------------------------------- */

      }

    }


    /* polarization */
    if ((ppt->has_source_p == _TRUE_) && (index_type == ppt->index_tp_p)) {

      /* polarization scalars */
      if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

	if (pvecthermo[pth->index_th_g] != 0.) {

	  x2 = k2 * (pba->conformal_age-eta) * (pba->conformal_age-eta);

	  if (x2 != 0.) {

	    /* in all gauges */
	    
	    pvecsource_terms[index_type * ppw->st_size + ppw->index_st_S0] =
	      + 3./16. * pvecthermo[pth->index_th_g] * Pi / x2;   /* /x2; */ 
	
	    /* 	    pvecsource_terms[index_type * ppw->st_size + ppw->index_st_S1] = 0.; */
	    pvecsource_terms[index_type * ppw->st_size + ppw->index_st_dS1] = 0.;


	    /* 	    pvecsource_terms[index_type * ppw->st_size + ppw->index_st_S2] = 0.; */
	    /*  	    pvecsource_terms[index_type * ppw->st_size + ppw->index_st_dS2] = 0.; */
	    pvecsource_terms[index_type * ppw->st_size + ppw->index_st_ddS2] = 0.; 

	  }
	  else {
	    pvecsource_terms[index_type * ppw->st_size + ppw->index_st_S0] = 0.;
	    pvecsource_terms[index_type * ppw->st_size + ppw->index_st_dS1] = 0.;
	    pvecsource_terms[index_type * ppw->st_size + ppw->index_st_ddS2] = 0.;
	  }
	}
	else {
	  pvecsource_terms[index_type * ppw->st_size + ppw->index_st_S0] = 0.;
	  pvecsource_terms[index_type * ppw->st_size + ppw->index_st_dS1] = 0.;
	  pvecsource_terms[index_type * ppw->st_size + ppw->index_st_ddS2] = 0.;
	}
      }

      /* polarization tensors */
      if ((ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors)) {

	/* -----------------  */

      }
    }

    /* gravitational potential (scalars only) */
    if ((ppt->has_source_g == _TRUE_) && (index_type == ppt->index_tp_g)) {
      
      if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

	/* newtonian gauge */
	if (ppr->gauge == newtonian) {
	  pvecsource_terms[index_type * ppw->st_size + ppw->index_st_S0] = 
	    pvecmetric[ppw->index_mt_phi];
	  pvecsource_terms[index_type * ppw->st_size + ppw->index_st_dS1] = 0.;
	  pvecsource_terms[index_type * ppw->st_size + ppw->index_st_ddS2] = 0.;
	}

	/* synchronous gauge */
	if (ppr->gauge == synchronous) {
	  pvecsource_terms[index_type * ppw->st_size + ppw->index_st_S0] = 
	    (a_prime_over_a * (pvecmetric[ppw->index_mt_h_prime] + 6. * pvecmetric[ppw->index_mt_eta_prime])/2./k2 + pvecmetric[ppw->index_mt_alpha_prime]);
	  pvecsource_terms[index_type * ppw->st_size + ppw->index_st_dS1] = 0.;
	  pvecsource_terms[index_type * ppw->st_size + ppw->index_st_ddS2] = 0.;
	}
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
  
  for (index_type = 0; index_type < ppt->tp_size; index_type++) {

    /** - for scalar temperature, infer \f$ S_2'' \f$ from \f$ S_2' \f$ at each time with array_derive1_order2_table_line_to_line() */

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

    /** - for each time, sum up \f$ S = S_0 + S_1' + S_2'' \f$ and store in array ((sources[index_mode])[index_ic][index_type])[index_eta][index_k] */

    for (index_eta = 0; index_eta < ppt->eta_size; index_eta++) {

      ppt->sources[index_mode]
	[index_ic * ppt->tp_size + index_type]
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
  double * pvecback;
  double * pvecthermo;
  double * pvecmetric;

  pppaw = parameters_and_workspace;

  k = pppaw->k;
  index_mode = pppaw->index_mode;
  ppr = pppaw->ppr;
  pba = pppaw->pba;
  pth = pppaw->pth;
  ppt = pppaw->ppt;

  ppw = pppaw->ppw;
  pvecback = ppw->pvecback;
  pvecthermo = ppw->pvecthermo;
  pvecmetric = ppw->pvecmetric;
  

  k2 = k*k;

  class_call(background_at_eta(pba,
			       eta, 
			       short_info, 
			       closeby, 
			       ppw->last_index_back, 
			       pvecback),
	     pba->error_message,
	     error_message);

  class_call(thermodynamics_at_z(pba,
				 pth,
				 1./pvecback[pba->index_bg_a]-1.,  /* redshift z=1/a-1 */
				 closeby,
				 ppw->last_index_thermo,
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
  if (ppt->has_scalars && index_mode == ppt->index_md_scalars) {

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
    h_plus_six_eta_prime = pvecmetric[ppw->index_mt_h_prime] + 6. * pvecmetric[ppw->index_mt_eta_prime];

    /** (b) if some approximation schemes are turned on, enforce a few y[] values */

    /** (b.a) photon shear and higher photon momenta if photon tight-coupling is on */
    if (ppw->tca == tca_on) {
      if (ppr->gauge == newtonian)
	y[ppw->index_pt_shear_g]=(8./3.*y[ppw->index_pt_theta_g])/9./pvecthermo[pth->index_th_dkappa]; /* tight-coupling shear_g */
      if (ppr->gauge == synchronous)
	y[ppw->index_pt_shear_g]=(8./3.*y[ppw->index_pt_theta_g]+4./3.*h_plus_six_eta_prime)/9./pvecthermo[pth->index_th_dkappa]; /* tight-coupling shear_g */      
      
      /* gauge-independent approximation for higher momenta */
      y[ppw->index_pt_l3_g]=6./7.*k/pvecthermo[pth->index_th_dkappa]*y[ppw->index_pt_shear_g];
      y[ppw->index_pt_pol0_g]=2.5*y[ppw->index_pt_shear_g];
      y[ppw->index_pt_pol1_g]=7./12.*y[ppw->index_pt_l3_g];
      y[ppw->index_pt_pol2_g]=0.5*y[ppw->index_pt_shear_g];
      y[ppw->index_pt_pol3_g]=0.25*y[ppw->index_pt_l3_g];
    }

    /** (b.b) radiation temperature/polarization momenta when radiation perturbations are off (exclusive of photon tight-coupling being on) */
    if (ppw->rp == rp_off) {

      /* analytic free-streaming solution (during matter domination and inside Hubble) */
      if (ppr->gauge == newtonian) {
	/* Newtonian gauge : */

      }
      if (ppr->gauge == synchronous) {
	/* Synchronous gauge : */
	y[ppw->index_pt_delta_g] = -4.*pvecmetric[ppw->index_mt_alpha_prime];
	y[ppw->index_pt_theta_g] = -0.5*pvecmetric[ppw->index_mt_h_prime];
	if (pba->has_nur == _TRUE_) {
	  y[ppw->index_pt_delta_nur] = -4.*pvecmetric[ppw->index_mt_alpha_prime];
	  y[ppw->index_pt_theta_nur] = -0.5*pvecmetric[ppw->index_mt_h_prime];
	}
      }

      for (l = 2; l < ppw->l_max_g; l++) {
	y[ppw->index_pt_delta_g+l] = 0.;
      }
      
      for (l = 0; l < ppw->l_max_pol_g; l++) {
	y[ppw->index_pt_pol0_g+l] = 0.;
      }

      if (pba->has_nur == _TRUE_) {
	for (l = 2; l < ppw->l_max_nur; l++) {
	  y[ppw->index_pt_delta_nur+l] = 0.;
	}
      }

    }

    /** (c) Photon temperature density and baryon density (do not depend on tight-coupling approximation) */
    if (ppr->gauge == newtonian) {
      /* Newtonian gauge : */
      if (ppw->rp == rp_on) {
	dy[ppw->index_pt_delta_g] = /* photon density */
	  -4./3.*y[ppw->index_pt_theta_g]+4.*pvecmetric[ppw->index_mt_phi_prime];
      }
      else {
	dy[ppw->index_pt_delta_g] = 0.; /* photon density when switched off */
      }
      dy[ppw->index_pt_delta_b] = /* baryon density */
	-y[ppw->index_pt_theta_b] + 3.*pvecmetric[ppw->index_mt_phi_prime];
    }
    
    if (ppr->gauge == synchronous) {
      /* Synchronous gauge : */
      if (ppw->rp == rp_on) {
	dy[ppw->index_pt_delta_g] = /* photon density */
	  -4./3.*y[ppw->index_pt_theta_g] - 2./3.*pvecmetric[ppw->index_mt_h_prime];
      }
      else {
	dy[ppw->index_pt_delta_g] = 0.; /* photon density when switched off */
      }
      dy[ppw->index_pt_delta_b] = /* baryon density */
	-y[ppw->index_pt_theta_b] - 0.5*pvecmetric[ppw->index_mt_h_prime];
    }

    /** (d) Baryon velocity (depend on tight-coupling approximation) */

    /** (d.1) Baryon velocity if baryon tight-coupling is off */

    if (ppw->tca == tca_off) {

      if (ppr->gauge == newtonian)
	/* Newtonian gauge : */
	dy[ppw->index_pt_theta_b] = /* baryon velocity */
	  - a_prime_over_a*y[ppw->index_pt_theta_b] 
	  + k2*pvecmetric[ppw->index_mt_psi] 
	  + pvecthermo[pth->index_th_cb2]*k2*y[ppw->index_pt_delta_b]
	  + R*pvecthermo[pth->index_th_dkappa] * (y[ppw->index_pt_theta_g]
						  -y[ppw->index_pt_theta_b]);
      if (ppr->gauge == synchronous)
	/* Synchronous gauge : */
	dy[ppw->index_pt_theta_b] = /* baryon velocity */
	  - a_prime_over_a*y[ppw->index_pt_theta_b] 
	  + pvecthermo[pth->index_th_cb2]*k2*y[ppw->index_pt_delta_b]
	  + R*pvecthermo[pth->index_th_dkappa] * (y[ppw->index_pt_theta_g]
						  -y[ppw->index_pt_theta_b]);
    }

    /** (d.2) Baryon velocity if baryon tight-coupling is on */

    else {

      if (ppr->gauge == newtonian) {
	/* Newtonian gauge : */
	slip=(2.*R/(1.+R)*a_prime_over_a+pvecthermo[pth->index_th_ddkappa]/pvecthermo[pth->index_th_dkappa]) /* tight-coupling (theta_b-theta_g)' */
	  *(y[ppw->index_pt_theta_b]-y[ppw->index_pt_theta_g])
	  +(-a_primeprime_over_a*y[ppw->index_pt_theta_b]
	    -a_prime_over_a*k2*(y[ppw->index_pt_delta_g]/2.+pvecmetric[ppw->index_mt_psi])
	    +k2*(pvecthermo[pth->index_th_cb2]*dy[ppw->index_pt_delta_b]
		 -dy[ppw->index_pt_delta_g]/4.)
	    )/pvecthermo[pth->index_th_dkappa]/(1.+R);

	dy[ppw->index_pt_theta_b] = /* tight-coupling baryon velocity */
	  (-a_prime_over_a*y[ppw->index_pt_theta_b]
	   +pvecthermo[pth->index_th_cb2]*k2*y[ppw->index_pt_delta_b]
	   +k2*R*(y[ppw->index_pt_delta_g]/4.-y[ppw->index_pt_shear_g])
	   +R*slip)/(1.+R)
	  +k2*pvecmetric[ppw->index_mt_psi];
      }

      if (ppr->gauge == synchronous) {
	/* Synchronous gauge : */

	slip=(2.*R/(1.+R)*a_prime_over_a+pvecthermo[pth->index_th_ddkappa]/pvecthermo[pth->index_th_dkappa]) /* tight-coupling (theta_b-theta_g)' */
	  *(y[ppw->index_pt_theta_b]-y[ppw->index_pt_theta_g])
	  +(-a_primeprime_over_a*y[ppw->index_pt_theta_b]
	    -a_prime_over_a*k2*(y[ppw->index_pt_delta_g]/2.)
	    +k2*(pvecthermo[pth->index_th_cb2]*dy[ppw->index_pt_delta_b]
		 -dy[ppw->index_pt_delta_g]/4.)
	    )/pvecthermo[pth->index_th_dkappa]/(1.+R);

	/* tight-coupling (theta_b-theta_g)' */
	/* 	slip=2.*R/(1.+R)*a_prime_over_a*(y[ppw->index_pt_theta_b]-y[ppw->index_pt_theta_g]) */
	/* 	  +(-a_primeprime_over_a*y[ppw->index_pt_theta_b] */
	/* 	    -a_prime_over_a*k2*y[ppw->index_pt_delta_g]/2. */
	/* 	    +k2*(pvecthermo[pth->index_th_cb2]*dy[ppw->index_pt_delta_b] */
	/* 		 -dy[ppw->index_pt_delta_g]/4.) */
	/* 	    )/pvecthermo[pth->index_th_dkappa]/(1.+R); */

	/* for testing */
	/*printf("%e %e\n",1./a-1.,pvecthermo[pth->index_th_ddkappa]/pvecthermo[pth->index_th_dkappa]);*/

	dy[ppw->index_pt_theta_b] = /* tight-coupling baryon velocity */
	  (-a_prime_over_a*y[ppw->index_pt_theta_b]
	   +pvecthermo[pth->index_th_cb2]*k2*y[ppw->index_pt_delta_b]
	   +k2*R*(y[ppw->index_pt_delta_g]/4.-y[ppw->index_pt_shear_g])
	   +R*slip)/(1.+R);
      }
    }


    /** (e) Photon temperature higher momenta and photon polarisation (depend on tight-coupling approximation) : */

    /** (e.1) if photon tight-coupling is off: */ 
    if (ppw->tca == tca_off) {

      /** (e.1.a) define \f$ \Pi = G_{\gamma 0} + G_{\gamma 2} + F_{\gamma 2} \f$ */
      Pi = y[ppw->index_pt_pol0_g] + y[ppw->index_pt_pol2_g] + 2.*y[ppw->index_pt_shear_g];

      /** (e.1.b) Photon velocity and shear */ 

      if (ppr->gauge == newtonian) {
	/* Newtonian gauge : */
	if (ppw->rp == rp_on) {
	  dy[ppw->index_pt_theta_g] = /* photon velocity */
	    k2*(y[ppw->index_pt_delta_g]/4.
		-y[ppw->index_pt_shear_g]+pvecmetric[ppw->index_mt_psi])
	    +pvecthermo[pth->index_th_dkappa]*(y[ppw->index_pt_theta_b]
					       -y[ppw->index_pt_theta_g]);
	  dy[ppw->index_pt_shear_g] = /* photon shear */
	    0.5*(8./15.*y[ppw->index_pt_theta_g]
		 -3./5.*k*y[ppw->index_pt_shear_g+1]
		 -pvecthermo[pth->index_th_dkappa]*(2.*y[ppw->index_pt_shear_g]-1./10.*Pi));
	}
	else {
	  dy[ppw->index_pt_theta_g] = 0.; /* photon velocity when switched off */
	  dy[ppw->index_pt_shear_g] = 0.; /* photon shear when switched off */
	}

      }
      
      if (ppr->gauge == synchronous) {
	/* Synchronous gauge : */
	if (ppw->rp == rp_on) {
	  dy[ppw->index_pt_theta_g] = /* photon velocity */
	    k2*(y[ppw->index_pt_delta_g]/4.
		-y[ppw->index_pt_shear_g])
	    + pvecthermo[pth->index_th_dkappa]*(y[ppw->index_pt_theta_b]
						-y[ppw->index_pt_theta_g]);

	  dy[ppw->index_pt_shear_g] = /* photon shear */
	    0.5*(8./15.*y[ppw->index_pt_theta_g]
		 -3./5.*k*y[ppw->index_pt_shear_g+1]
		 +4./15.*pvecmetric[ppw->index_mt_h_prime]+8./5.*pvecmetric[ppw->index_mt_eta_prime]
		 -pvecthermo[pth->index_th_dkappa]*(2.*y[ppw->index_pt_shear_g]-1./10.*Pi));
	}
	else {
	  dy[ppw->index_pt_theta_g] = 0.; /* photon velocity when switched off */
	  dy[ppw->index_pt_shear_g] = 0.; /* photon shear when switched off */
	}	  

      }
      
      /** (e.1.c) Photon temperature higher momenta (l >=3), gauge-independent */ 

      if (ppw->rp == rp_on) {

	l = 3; /* photon l=3 (special case because F_gamma2=2*shear !!) */
	dy[ppw->index_pt_l3_g] =
	  k/(2.*l+1.)*(l*2.*y[ppw->index_pt_shear_g]-(l+1.)*y[ppw->index_pt_l3_g+1])
	  - pvecthermo[pth->index_th_dkappa]*y[ppw->index_pt_l3_g];

	for (l = 4; l < ppw->l_max_g; l++) { /* photon additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3) */
	  dy[ppw->index_pt_delta_g+l] =
	    k/(2.*l+1)*(l*y[ppw->index_pt_delta_g+l-1]-(l+1.)*y[ppw->index_pt_delta_g+l+1])
	    - pvecthermo[pth->index_th_dkappa]*y[ppw->index_pt_delta_g+l];
	}

	l = ppw->l_max_g; /* l=lmax */
	dy[ppw->index_pt_delta_g+ppw->l_max_g] = /* last photon term */
	  k*y[ppw->index_pt_delta_g+ppw->l_max_g-1]
	  -(1.+l)/eta*y[ppw->index_pt_delta_g+ppw->l_max_g]
	  - pvecthermo[pth->index_th_dkappa]*y[ppw->index_pt_delta_g+ppw->l_max_g];

      }
      else {
	for (l = 3; l < ppw->l_max_g; l++) {
	  dy[ppw->index_pt_delta_g+l] = 0.; /* photon T high momenta when switched off */
	}

      }

      /** (e.1.d) Photon polarisation */

      if (ppw->rp == rp_on) {

	dy[ppw->index_pt_pol0_g] = /* photon polarization, l=0 */
	  -k*y[ppw->index_pt_pol0_g+1]
	  -pvecthermo[pth->index_th_dkappa]*(y[ppw->index_pt_pol0_g]-Pi/2.);
	dy[ppw->index_pt_pol1_g] = /* photon polarization, l=1 */
	  k/3.*(y[ppw->index_pt_pol1_g-1]-2.*y[ppw->index_pt_pol1_g+1])
	  -pvecthermo[pth->index_th_dkappa]*y[ppw->index_pt_pol1_g];
	dy[ppw->index_pt_pol2_g] = /* photon polarization, l=2 */
	  k/5.*(2.*y[ppw->index_pt_pol2_g-1]-3.*y[ppw->index_pt_pol2_g+1])
	  -pvecthermo[pth->index_th_dkappa]*(y[ppw->index_pt_pol2_g]-Pi/10.);

	for (l=3; l < ppw->l_max_pol_g; l++)
	  dy[ppw->index_pt_pol0_g+l] =  /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2) */
	    k/(2.*l+1)*(l*y[ppw->index_pt_pol0_g+l-1]-(l+1.)*y[ppw->index_pt_pol0_g+l+1])
	    -pvecthermo[pth->index_th_dkappa]*y[ppw->index_pt_pol0_g+l];

	l = ppw->l_max_pol_g;
	dy[ppw->index_pt_pol0_g+l] =  /* l=lmax */
	  k*y[ppw->index_pt_pol0_g+l-1]-(l+1)/eta*y[ppw->index_pt_pol0_g+l]
	  -pvecthermo[pth->index_th_dkappa]*y[ppw->index_pt_pol0_g+l];

      }
      else {

	for (l = 0; l < ppw->l_max_pol_g; l++) {
	  dy[ppw->index_pt_pol0_g+l] = 0.; /* photon polarization when switched off */
	}
      }

    }

    /** (e.2) if photon tight-coupling is on: */
 
    else {

      /** (e.2.a) photon velocity */

      if (ppr->gauge == newtonian)
	/* Newtonian gauge : */
	dy[ppw->index_pt_theta_g] = /* tight-coupling photon velocity */
	  -(dy[ppw->index_pt_theta_b]+a_prime_over_a*y[ppw->index_pt_theta_b]-pvecthermo[pth->index_th_cb2]*k2*y[ppw->index_pt_delta_b])/R
	  +k2*(0.25*y[ppw->index_pt_delta_g]-y[ppw->index_pt_shear_g])+(1.+R)/R*k2*pvecmetric[ppw->index_mt_psi];

      if (ppr->gauge == synchronous)
	/* Synchronous gauge : */

	dy[ppw->index_pt_theta_g] = /* tight-coupling photon velocity */
	  -(dy[ppw->index_pt_theta_b]+a_prime_over_a*y[ppw->index_pt_theta_b]-pvecthermo[pth->index_th_cb2]*k2*y[ppw->index_pt_delta_b])/R
	  +k2*(y[ppw->index_pt_delta_g]/4.-y[ppw->index_pt_shear_g]);

      /*       	  dy[ppw->index_pt_theta_g] = /\* photon velocity (used here in CAMB, leading to more instabilities) *\/ */
      /*       	    k2*(y[ppw->index_pt_delta_g]/4.-y[ppw->index_pt_shear_g]) */
      /*       	    + pvecthermo[pth->index_th_dkappa]*(y[ppw->index_pt_theta_b] */
      /*       						 -y[ppw->index_pt_theta_g]); */

      /** (e.2.b) all other photon perturbations vanish */

      for (l = 2; l <= ppw->l_max_g; l++)
	dy[ppw->index_pt_delta_g+l]=0.; /* tight-coupling: all photon temperature terms excepted l=0,1 have zero derivs */

      for (l = 0; l <= ppw->l_max_pol_g; l++)
	dy[ppw->index_pt_pol0_g+l]=0.; /* tight-coupling: all photon polarisation terms have zero derivs */

    }

    /** (f) cdm */

    if (pba->has_cdm == _TRUE_) {  

      if (ppr->gauge == newtonian) {
	/* Newtonian gauge : */
	dy[ppw->index_pt_delta_cdm] = /* cdm density */
	  -y[ppw->index_pt_theta_cdm]+3.*pvecmetric[ppw->index_mt_phi_prime];
	dy[ppw->index_pt_theta_cdm] = /* cdm velocity */
	  - a_prime_over_a*y[ppw->index_pt_theta_cdm] + k2*pvecmetric[ppw->index_mt_psi];
      }

      if (ppr->gauge == synchronous) {
	/* Synchronous gauge : */
	dy[ppw->index_pt_delta_cdm] = /* cdm density */
	  -0.5*pvecmetric[ppw->index_mt_h_prime];
      }

    }
    
    /** (g) dark energy fluid */
    
    if (pba->has_dark_energy_fluid == _TRUE_) {  

      double ache_prime = a_prime_over_a;
      double cs2 = 1.;

      if (ppr->gauge == newtonian) {
	/* Newtonian gauge : */
	dy[ppw->index_pt_delta_de] = /* dark energy density */
	  (-3*(1+ pba->w_de )*ache_prime-3*pvecback[pba->index_bg_H]*(cs2- pba->w_de )*(y[ppw->index_pt_delta_de]/pvecback[pba->index_bg_rho_de]+3*pvecback[pba->index_bg_H]*(1+ pba->w_de )*y[ppw->index_pt_theta_de]/k)-(1+ pba->w_de )*k*y[ppw->index_pt_theta_de])/pvecback[pba->index_bg_rho_de]; // 0;

	dy[ppw->index_pt_theta_de] = /* dark energy velocity */
	  (k*cs2*y[ppw->index_pt_delta_de])/(pvecback[pba->index_bg_rho_de]*(1+ pba->w_de ))-pvecback[pba->index_bg_H]*(1-3*cs2)*y[ppw->index_pt_theta_de]+k*pvecmetric[ppw->index_mt_psi]; // 0;
      }

      if (ppr->gauge == synchronous) {
	/* Synchronous gauge : */
	dy[ppw->index_pt_delta_de] = /* dark energy density */
	  (-3*(1+ pba->w_de )*ache_prime-3*pvecback[pba->index_bg_H]*(cs2- pba->w_de )*(y[ppw->index_pt_delta_de]/pvecback[pba->index_bg_rho_de]+3*pvecback[pba->index_bg_H]*(1+ pba->w_de )*y[ppw->index_pt_theta_de]/k)-(1+ pba->w_de )*k*y[ppw->index_pt_theta_de])/pvecback[pba->index_bg_rho_de]; // 0;

	dy[ppw->index_pt_theta_de] = /* dark energy velocity */
	  (k*cs2*y[ppw->index_pt_delta_de])/(pvecback[pba->index_bg_rho_de]*(1+ pba->w_de ))-pvecback[pba->index_bg_H]*(1-3*cs2)*y[ppw->index_pt_theta_de]; // 0;
      }
      
    }  
    
    /** (h) ultra-relativistic neutrino/relics density, velocity, shear, etc. */
    
    if (pba->has_nur == _TRUE_) {
      
      if (ppr->gauge == newtonian) {

	if(ppw->rp == rp_on) {
	  /* Newtonian gauge : */
	  dy[ppw->index_pt_delta_nur] = /* density of ultra-relativistic neutrinos/relics */
	    -4./3.*y[ppw->index_pt_theta_nur] + 4.*pvecmetric[ppw->index_mt_phi_prime];

	  dy[ppw->index_pt_theta_nur] = /* velocity of ultra-relativistic neutrinos/relics */
	    k2*(y[ppw->index_pt_delta_nur]/4.
		-y[ppw->index_pt_shear_nur]+pvecmetric[ppw->index_mt_psi]);

	  dy[ppw->index_pt_shear_nur] = /* shear of ultra-relativistic neutrinos/relics */
	    0.5*(8./15.*y[ppw->index_pt_theta_nur]
		 -3./5.*k*y[ppw->index_pt_shear_nur+1]);
	}
	else {
	  dy[ppw->index_pt_delta_nur] = 0.; /* density of ultra-relativistic neutrinos/relics when switched off*/
	  dy[ppw->index_pt_theta_nur] = 0.; /* velocity of ultra-relativistic neutrinos/relics when switched off*/
	  dy[ppw->index_pt_shear_nur] = 0.; /* shear of ultra-relativistic neutrinos/relics when switched off*/
	}

      }

      if (ppr->gauge == synchronous) {

	if(ppw->rp == rp_on) {
	  /* Synchronous gauge : */
	  dy[ppw->index_pt_delta_nur] = /* density of ultra-relativistic neutrinos/relics */
	    -4./3.*y[ppw->index_pt_theta_nur] - 2./3.*pvecmetric[ppw->index_mt_h_prime];
	
	  dy[ppw->index_pt_theta_nur] = /* velocity of ultra-relativistic neutrinos/relics */
	    k2*(y[ppw->index_pt_delta_nur]/4.
		-y[ppw->index_pt_shear_nur]);

	  dy[ppw->index_pt_shear_nur] = /* shear of ultra-relativistic neutrinos/relics */
	    0.5*(8./15.*y[ppw->index_pt_theta_nur]
		 -3./5.*k*y[ppw->index_pt_shear_nur+1]
		 +4./15.*pvecmetric[ppw->index_mt_h_prime]+8./5.*pvecmetric[ppw->index_mt_eta_prime]);
	}
	else {
	  dy[ppw->index_pt_delta_nur] = 0.; /* density of ultra-relativistic neutrinos/relics when switched off*/
	  dy[ppw->index_pt_theta_nur] = 0.; /* velocity of ultra-relativistic neutrinos/relics when switched off*/
	  dy[ppw->index_pt_shear_nur] = 0.; /* shear of ultra-relativistic neutrinos/relics when switched off*/
	}

      }

      /** (i) ultra-relativistic neutrino/relics higher momenta l>=3 if radiation perturbations are on */

      if(ppw->rp == rp_on) {
	l = 3;
	dy[ppw->index_pt_l3_nur] = /* l=3 of ultra-relativistic neutrinos/relics (special case because F_gamma2=2*shear !!) */
	  k/(2.*l+1.)*(l*2.*y[ppw->index_pt_shear_nur]-(l+1.)*y[ppw->index_pt_l3_nur+1]);

	for (l = 4; l < ppw->l_max_nur; l++) {
	  dy[ppw->index_pt_delta_nur+l] = /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3) */
	    k/(2.*l+1)*(l*y[ppw->index_pt_delta_nur+l-1]-(l+1.)*y[ppw->index_pt_delta_nur+l+1]);
	}

	l = ppw->l_max_nur; /* l=lmax */
	dy[ppw->index_pt_delta_nur+ppw->l_max_nur] = /* last term of ultra-relativistic neutrinos/relics */
	  k/(2.*l+1)*(l*y[ppw->index_pt_delta_nur+ppw->l_max_nur-1]-(l+1.)*
		      ((2.*l+1)/k/eta*y[ppw->index_pt_delta_nur+ppw->l_max_nur]-y[ppw->index_pt_delta_nur+ppw->l_max_nur-1]));
      }

      else{
	for (l = 3; l < ppw->l_max_nur; l++) {
	  dy[ppw->index_pt_delta_nur+l] = 0.;
	}
      }
     
    }

    /** (j) metric */

    if (ppr->gauge == synchronous) {
      /* Synchronous gauge */
      dy[ppw->index_pt_eta] = pvecmetric[ppw->index_mt_eta_prime];

      /* for testing, will be useful for improving tight-coupling approximation */
      /*       if ((index_k == 0) && (eta > 800)) */
      /* 	printf("%e %e %e %e %e %e %e %e %e %e %e\n", */
      /* 	       eta, */
      /* 	       y[ppw->index_pt_delta_g], */
      /* 	       y[ppw->index_pt_theta_g], */
      /* 	       y[ppw->index_pt_shear_g], */
      /* 	       y[ppw->index_pt_l3_g], */
      /* 	       y[ppw->index_pt_pol0_g], */
      /* 	       y[ppw->index_pt_pol1_g], */
      /* 	       y[ppw->index_pt_pol2_g], */
      /* 	       y[ppw->index_pt_pol3_g], */
      /* 	       y[ppw->index_pt_delta_b], */
      /* 	       y[ppw->index_pt_theta_b]); */
	       
    }

  }

  /** - tensor mode */

  if (ppt->has_tensors && index_mode == ppt->index_md_tensors) {

    /* photon density */
    dy[ppw->index_pt_delta_g] = 
      -4./3.*y[ppw->index_pt_theta_g]
      -4.*y[ppw->index_pt_gwdot]
      -pvecthermo[pth->index_th_dkappa]*(9./10.*y[ppw->index_pt_theta_g]+12./5.*y[ppw->index_pt_pol0_g]);

    /* photon velocity */
    dy[ppw->index_pt_theta_g] = 
      k2*(y[ppw->index_pt_delta_g]/4.-y[ppw->index_pt_shear_g])
      -pvecthermo[pth->index_th_dkappa]*y[ppw->index_pt_theta_g];

    /* photon shear */
    dy[ppw->index_pt_shear_g] =	
      0.5*(8./15.*y[ppw->index_pt_theta_g]
	   -3./5.*k*y[ppw->index_pt_shear_g+1]
	   -pvecthermo[pth->index_th_dkappa]*(68./35.*y[ppw->index_pt_shear_g]-6./35.*y[ppw->index_pt_pol2_g]));

    /* photon l=3 */
    dy[ppw->index_pt_l3_g] = 
      k/7.*(6.*y[ppw->index_pt_shear_g]-4.*y[ppw->index_pt_l3_g+1])
      -pvecthermo[pth->index_th_dkappa]*y[ppw->index_pt_l3_g];
      
    /* photon l=4 */
    l=4;
    dy[ppw->index_pt_delta_g+l] = 
      k/(2.*l+1.)*(l*y[ppw->index_pt_delta_g+l-1]-(l+1.)*y[ppw->index_pt_delta_g+l+1])
      -pvecthermo[pth->index_th_dkappa]*(209./210.*y[ppw->index_pt_l3_g+1]+1./210.*y[ppw->index_pt_pol2_g+2]);


    /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3,4) */
    for (l=5; l < ppw->l_max_g; l++)
      dy[ppw->index_pt_delta_g+l] = 
	k/(2.*l+1.)*(l*y[ppw->index_pt_delta_g+l-1]-(l+1.)*y[ppw->index_pt_delta_g+l+1])
	-pvecthermo[pth->index_th_dkappa]*y[ppw->index_pt_delta_g+l];  

    /* l=lmax */
    l = ppw->l_max_g;
    dy[ppw->index_pt_delta_g+l] = 
      k*y[ppw->index_pt_delta_g+l-1]
      -(1.+l)/eta*y[ppw->index_pt_delta_g+l]
      - pvecthermo[pth->index_th_dkappa]*y[ppw->index_pt_delta_g+l];

    /* photon polarization, l=0 */
    dy[ppw->index_pt_pol0_g] = 
      -k*y[ppw->index_pt_pol0_g+1]
      -pvecthermo[pth->index_th_dkappa]*(2./5.*y[ppw->index_pt_pol0_g]+1./10.*y[ppw->index_pt_delta_g]); 
    
    /* photon polarization, l=1 */
    dy[ppw->index_pt_pol1_g] = 
      k/3.*(y[ppw->index_pt_pol1_g-1]-2.*y[ppw->index_pt_pol1_g+1])
      -pvecthermo[pth->index_th_dkappa]*y[ppw->index_pt_pol1_g]; 

    /* photon polarization, l=2 */
    dy[ppw->index_pt_pol2_g] = 
      k/5.*(2.*y[ppw->index_pt_pol2_g-1]-3.*y[ppw->index_pt_pol2_g+1])
      -pvecthermo[pth->index_th_dkappa]*(41./35.*y[ppw->index_pt_pol2_g]+2./35.*y[ppw->index_pt_shear_g]); 

    /* photon polarization, l=3 */
    dy[ppw->index_pt_pol3_g] = 
      k/7.*(3.*y[ppw->index_pt_pol3_g-1]-4.*y[ppw->index_pt_pol3_g+1])
      -pvecthermo[pth->index_th_dkappa]*y[ppw->index_pt_pol3_g]; 

    /* photon polarization, l=4 */
    l=4;
    dy[ppw->index_pt_pol0_g+l] = 
      k/(2.*l+1.)*(l*y[ppw->index_pt_pol0_g+l-1]-(l+1.)*y[ppw->index_pt_pol0_g+l+1])
      -pvecthermo[pth->index_th_dkappa]*(209./210.*y[ppw->index_pt_pol0_g+l]+1./210.*y[ppw->index_pt_delta_g+l]); 

    /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3,4) */
    for (l=5; l < ppw->l_max_pol_g; l++)
      dy[ppw->index_pt_pol0_g+l] = 
	k/(2.*l+1.)*(l*y[ppw->index_pt_pol0_g+l-1]-(l+1.)*y[ppw->index_pt_pol0_g+l+1])
	-pvecthermo[pth->index_th_dkappa]*y[ppw->index_pt_pol0_g+l];

    /* l=lmax */
    l = ppw->l_max_pol_g;
    dy[ppw->index_pt_pol0_g+l] = 
      k*y[ppw->index_pt_pol0_g+l-1]
      -(l+1.)/eta*y[ppw->index_pt_pol0_g+l]
      -pvecthermo[pth->index_th_dkappa]*y[ppw->index_pt_pol0_g+l];

      /* tensor metric perturbation h (gravitational waves) */
    dy[ppw->index_pt_gw] = y[ppw->index_pt_gwdot];     

    /* its time-derivative */
    dy[ppw->index_pt_gwdot] = -2.*a_prime_over_a*y[ppw->index_pt_gwdot]-k2*y[ppw->index_pt_gw];
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

