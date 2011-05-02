/** @file perturbations.c Documented perturbation module
 *
 * Julien Lesgourgues, 23.09.2010    
 *
 * Deals with the perturbation evolution.
 * This mdule has two purposes: 
 *
 * - at the beginning, to initialize the perturbations, i.e. to
 * integrate the perturbation equations, and store temporarily the terms
 * contributing to the source functions as a function of conformal
 * time. Then, to perform a few manipulations of these terms in order to
 * infer the actual source functions \f$ S^{X} (k, \tau) \f$, and to
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
 * -# perturb_sources_at_tau() at any later time
 * -# perturb_free() at the end, when no more calls to perturb_sources_at_tau() are needed
 */

#include "perturbations.h"

/** 
 * Source function \f$ S^{X} (k, \tau) \f$ at a given conformal time tau.
 *
 * Evaluate all source functions at given conformal time tau by reading
 * the pre-computed table and interpolating.
 *
 * @param ppt        Input : pointer to perturbation structure containing interpolation tables
 * @param index_mode Input : index of requested mode
 * @param index_ic   Input : index of requested initial condition
 * @param index_k    Input : index of requested wavenumber
 * @param index_type Input : index of requested source function type
 * @param tau        Input : any value of conformal time
 * @param psource    Output: vector (assumed to be already allocated) of source functions
 * @return the error status
 */

int perturb_sources_at_tau(
			   struct perturbs * ppt,
			   int index_mode,
			   int index_ic,
			   int index_k,
			   int index_type,
			   double tau,
			   double * psource
			   ) {

  /** Summary: */

  /** - interpolate in pre-computed table contained in ppt */
  class_call(array_interpolate_two(&(ppt->sources[index_mode]
				     [index_ic * ppt->tp_size[index_mode] + index_type]
				     [index_k * ppt->tau_size]),
				   1,
				   0,
				   ppt->tau_sampling,
				   1,
				   ppt->tau_size,
				   tau,
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
  /* pointer to one struct perturb_workspace per thread (one if no openmp) */
  struct perturb_workspace ** pppw;
  /* number of threads (always one if no openmp) */
  int number_of_threads=1;
  /* index of the thread (always 0 if no openmp) */
  int thread=0;

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

  class_test((ppr->gauge == synchronous) && (pba->has_cdm == _FALSE_),
	     ppt->error_message,
	     "In the present synchronous gauge, it is not self-consistent to assume no CDM: the later is used to define the initial timelike hypersurface");

  class_test ((ppr->tight_coupling_approximation < first_order_MB) ||
	      (ppr->tight_coupling_approximation > compromise_CLASS),
	      ppt->error_message,
	      "your tight_coupling_approximation is set to %d, out of range defined in perturbations.h",ppr->tight_coupling_approximation);

  class_test ((ppr->radiation_streaming_approximation < rsa_null) ||
	      (ppr->radiation_streaming_approximation > rsa_none),
	      ppt->error_message,
	      "your radiation_streaming_approximation is set to %d, out of range defined in perturbations.h",ppr->radiation_streaming_approximation);
  
  if (pba->has_ur == _TRUE_) {

    class_test ((ppr->ur_fluid_approximation < ufa_mb) ||
		(ppr->ur_fluid_approximation > ufa_none),
		ppt->error_message,
		"your ur_fluid_approximation is set to %d, out of range defined in perturbations.h",ppr->ur_fluid_approximation);
  }

  if (pba->has_ncdm == _TRUE_) {

    class_test ((ppr->ncdm_fluid_approximation < ncdmfa_mb) ||
		(ppr->ncdm_fluid_approximation > ncdmfa_none),
		ppt->error_message,
		"your ncdm_fluid_approximation is set to %d, out of range defined in perturbations.h",ppr->ncdm_fluid_approximation);
  }

  class_test(ppt->has_vectors == _TRUE_,
	     ppt->error_message,
	     "Vectors not coded yet");

  if ((ppt->has_bi == _TRUE_) || (ppt->has_cdi == _TRUE_) || (ppt->has_nid == _TRUE_) || (ppt->has_niv == _TRUE_)) {
    printf("Warning: so far, isocurvature initial conditions not tested as thougoughly as adiabatic ones\n");
  }

  if (ppt->has_tensors == _TRUE_)
    printf("Warning: don't trust polarized tensors from version %s, we are still otpimizing them\n",_VERSION_);

  if ((ppt->has_cl_cmb_temperature == _TRUE_) && (ppt->has_cl_cmb_polarization == _TRUE_) &&
      (ppt->has_tensors == _TRUE_)) {
    printf("Warning: our C_l^TE for tensors has a minus sign with respect to CAMB 2008.\n");
  }

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
#endif

  class_alloc(pppw,number_of_threads * sizeof(struct perturb_workspace *),ppt->error_message);

  /** - loop over modes (scalar, tensors, etc). For each mode: */

  for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {

    abort = _FALSE_;

#pragma omp parallel				\
  shared(pppw,ppr,pba,pth,ppt,index_mode,abort)	\
  private(thread)

    {

#ifdef _OPENMP
      thread=omp_get_thread_num();
#endif

      /** create a workspace (one per thread in multi-thread case) */

      class_alloc_parallel(pppw[thread],sizeof(struct perturb_workspace),ppt->error_message);

      /** (a) initialize indices of vectors of perturbations with perturb_indices_of_current_vectors() */

      class_call_parallel(perturb_workspace_init(ppr,
						 pba,
						 pth,
						 ppt,
						 index_mode,
					 	 pppw[thread]),
			  ppt->error_message,
			  ppt->error_message);

    } /* end of parallel region */

    if (abort == _TRUE_) return _FAILURE_;

    /** (c) loop over initial conditions and wavenumbers; for each of them, evolve perturbations and compute source functions with perturb_solve() */

    for (index_ic = 0; index_ic < ppt->ic_size[index_mode]; index_ic++) {

      abort = _FALSE_;

#pragma omp parallel						\
  shared(pppw,ppr,pba,pth,ppt,index_mode,index_ic,abort)	\
  private(index_k,thread,tstart,tstop,tspent)

      {

#ifdef _OPENMP
	thread=omp_get_thread_num();
	tspent=0.;
#endif
	
#pragma omp for schedule (dynamic)

        /* integrating background is slightly more optimal for parallel runs */
	//for (index_k = 0; index_k < ppt->k_size[index_mode]; index_k++) {
	for (index_k = ppt->k_size[index_mode]-1; index_k >=0; index_k--) {  

	  if ((ppt->perturbations_verbose > 2) && (abort == _FALSE_))
	    printf("evolving mode k=%e /Mpc\n",(ppt->k[index_mode])[index_k]);
	  
#ifdef _OPENMP
	  tstart = omp_get_wtime();
#endif

	  class_call_parallel(perturb_solve(ppr,
					    pba,
					    pth,
					    ppt,
					    index_mode,
					    index_ic,
					    index_k,
					    pppw[thread]),
			      ppt->error_message,
			      ppt->error_message);

#ifdef _OPENMP
	  tstop = omp_get_wtime();

	  tspent += tstop-tstart;
#endif

#pragma omp flush(abort)

	} /* end of loop over wavenumbers */

#ifdef _OPENMP
	if (ppt->perturbations_verbose>1)
	  printf("In %s: time spent in parallel region (loop over k's) = %e s for thread %d\n",
		 __func__,tspent,omp_get_thread_num());
#endif

      } /* end of parallel region */

      if (abort == _TRUE_) return _FAILURE_;

    } /* end of loop over initial conditions */
    
    abort = _FALSE_;

#pragma omp parallel				\
  shared(pppw,ppt,index_mode,abort)		\
  private(thread)

    {

#ifdef _OPENMP
      thread=omp_get_thread_num();
#endif
      
      class_call_parallel(perturb_workspace_free(ppt,index_mode,pppw[thread]),
			  ppt->error_message,
			  ppt->error_message);

    } /* end of parallel region */

    if (abort == _TRUE_) return _FAILURE_;

  } /* end loop over modes */    

  free(pppw);

  return _SUCCESS_;
}

/**
 * Free all memory space allocated by perturb_init().
 * 
 * To be called at the end of each run, only when no further calls to
 * perturb_sources_at_tau() are needed.
 *
 * @param ppt Input: perturbation structure to be freed
 * @return the error status
 */

int perturb_free(
		 struct perturbs * ppt
		 ) {

  int index_mode,index_ic,index_type;

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
    
    free(ppt->tau_sampling);
	 
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

      if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) ||
	  ((ppt->has_pk_matter == _TRUE_) && (ppr->pk_definition == delta_tot_from_poisson_squared))) { 
	ppt->has_lss = _TRUE_;
	ppt->has_source_g = _TRUE_;
	ppt->index_tp_g = index_type; 
	index_type++;
      }
      else
	ppt->has_source_g = _FALSE_;

      if ((ppt->has_pk_matter == _TRUE_) && (ppr->pk_definition != delta_tot_from_poisson_squared)) {
	ppt->has_lss = _TRUE_;
        ppt->has_source_delta_pk = _TRUE_;
	ppt->index_tp_delta_pk = index_type; 
	index_type++;
      }
      else
	ppt->has_source_delta_pk = _FALSE_;
      
      if (ppt->has_matter_transfers == _TRUE_) {
	ppt->has_lss = _TRUE_;
	ppt->has_source_delta_g = _TRUE_;
	ppt->index_tp_delta_g = index_type;
	index_type++;
	ppt->has_source_delta_b = _TRUE_;
	ppt->index_tp_delta_b = index_type;
	index_type++;
	if (pba->has_cdm == _TRUE_) {
	  ppt->has_source_delta_cdm = _TRUE_;
	  ppt->index_tp_delta_cdm = index_type;
	  index_type++;
	}
	else ppt->has_source_delta_cdm = _FALSE_;
	if (pba->has_fld == _TRUE_) {
	  ppt->has_source_delta_fld = _TRUE_;
	  ppt->index_tp_delta_fld = index_type;
	  index_type++;
	}
	else ppt->has_source_delta_fld = _FALSE_;
	if (pba->has_ur == _TRUE_) {
	  ppt->has_source_delta_ur = _TRUE_;
	  ppt->index_tp_delta_ur = index_type;
	  index_type++;
	}
	else ppt->has_source_delta_ur = _FALSE_;
	if (pba->has_ncdm == _TRUE_) {
	  ppt->has_source_delta_ncdm = _TRUE_;
	  ppt->index_tp_delta_ncdm1 = index_type;
	  index_type+=pba->N_ncdm;
	}
	else ppt->has_source_delta_ncdm = _FALSE_;
      }
      else {
	ppt->has_source_delta_g = _FALSE_;
	ppt->has_source_delta_b = _FALSE_;
	ppt->has_source_delta_cdm = _FALSE_;
	ppt->has_source_delta_fld = _FALSE_;
	ppt->has_source_delta_ur = _FALSE_;
	ppt->has_source_delta_ncdm = _FALSE_;
      }

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

    /** (c) for each mode, define k values with perturb_get_k_list() */

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
 * For each type, compute the list of values of tau at which sources
 * will be sampled.  Knowing the number of tau values, allocate all
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
  int first_index_back;
  int first_index_thermo;

  double tau;
  double tau_ini;
  double tau_lower;
  double tau_upper;
  double tau_mid;

  double timescale_source;
  double rate_thermo;
  double rate_isw_squared;
  double a_prime_over_a;
  double a_primeprime_over_a;
  double * pvecback;
  double * pvecthermo;

  /** - allocate background/thermodynamics vectors */

  class_alloc(pvecback,pba->bg_size_short*sizeof(double),ppt->error_message);  
  class_alloc(pvecthermo,pth->th_size*sizeof(double),ppt->error_message);

  /** - first, just count the number of sampling points in order to allocate the array containing all values: */

  /** (a) if CMB requested, first sampling point = when the universe
      stops being opaque; otherwise, start sampling gravitational
      potential at recombination */
  if (ppt->has_cmb == _TRUE_) {

    /* using bisection, search time tau such that the ratio of thermo
       to Hubble time scales tau_c/tau_h=aH/kappa' is equal to
       start_sources_at_tau_c_over_tau_h */
    
    tau_lower = pth->tau_ini;

    class_call(background_at_tau(pba,
				 tau_lower, 
				 short_info, 
				 normal, 
				 &first_index_back, 
				 pvecback),
	       pba->error_message,
	       ppt->error_message);
    
    class_call(thermodynamics_at_z(pba,
				   pth,
				   1./pvecback[pba->index_bg_a]-1.,  /* redshift z=1/a-1 */
				   normal,
				   &first_index_thermo,
				   pvecback,
				   pvecthermo),
	       pth->error_message,
	       ppt->error_message);
    
    class_test(pvecback[pba->index_bg_a]*
	       pvecback[pba->index_bg_H]/
	       pvecthermo[pth->index_th_dkappa] > 
	       ppr->start_sources_at_tau_c_over_tau_h,
	       ppt->error_message,
	       "your choice of initial time for computing sources is inappropriate: it corresponds to an earlier time than the one at which the integration of thermodynamical variables started (tau=%g). You should increase either 'start_sources_at_tau_c_over_tau_h' or 'recfast_z_initial'\n",
	       tau_lower);
    
    
    tau_upper = pth->tau_rec;
    
    class_call(background_at_tau(pba,
				 tau_upper, 
				 short_info, 
				 normal, 
				 &first_index_back, 
				 pvecback),
	       pba->error_message,
	       ppt->error_message);
    
    class_call(thermodynamics_at_z(pba,
				   pth,
				   1./pvecback[pba->index_bg_a]-1.,  /* redshift z=1/a-1 */
				   normal,
				   &first_index_thermo,
				   pvecback,
				   pvecthermo),
	       pth->error_message,
	       ppt->error_message);

    class_test(pvecback[pba->index_bg_a]*
	       pvecback[pba->index_bg_H]/
	       pvecthermo[pth->index_th_dkappa] <
	       ppr->start_sources_at_tau_c_over_tau_h,
	       ppt->error_message,
	       "your choice of initial time for computing sources is inappropriate: it corresponds to a time after recombination. You should decrease 'start_sources_at_tau_c_over_tau_h'\n");
    
    tau_mid = 0.5*(tau_lower + tau_upper);
    
    while (tau_upper - tau_lower > ppr->tol_tau_approx) {

      class_call(background_at_tau(pba,
				   tau_mid, 
				   short_info, 
				   normal, 
				   &first_index_back, 
				   pvecback),
		 pba->error_message,
		 ppt->error_message);
      
      class_call(thermodynamics_at_z(pba,
				     pth,
				     1./pvecback[pba->index_bg_a]-1.,  /* redshift z=1/a-1 */
				     normal,
				     &first_index_thermo,
				     pvecback,
				     pvecthermo),
		 pth->error_message,
		 ppt->error_message);
      
      
      if (pvecback[pba->index_bg_a]*
	  pvecback[pba->index_bg_H]/
	  pvecthermo[pth->index_th_dkappa] > 
	  ppr->start_sources_at_tau_c_over_tau_h)
	
	tau_upper = tau_mid;
      else
	tau_lower = tau_mid;

      tau_mid = 0.5*(tau_lower + tau_upper);
      
    }

    tau_ini = tau_mid;

  }
  else {

    /* case when CMB not requested: start at recombination time */
    tau_ini = pth->tau_rec;

    /* set values of first_index_back/thermo */
    class_call(background_at_tau(pba,
				 tau_ini, 
				 short_info, 
				 normal, 
				 &first_index_back, 
				 pvecback),
	       pba->error_message,
	       ppt->error_message);
    
    class_call(thermodynamics_at_z(pba,
				   pth,
				   1./pvecback[pba->index_bg_a]-1.,  /* redshift z=1/a-1 */
				   normal,
				   &first_index_thermo,
				   pvecback,
				   pvecthermo),
	       pth->error_message,
	       ppt->error_message);
  }    

  
  counter = 1;

  /** (b) next sampling point = previous + ppr->perturb_sampling_stepsize * timescale_source, where:
      - if CMB requested:
      timescale_source1 = \f$ |g/\dot{g}| = |\dot{\kappa}-\ddot{\kappa}/\dot{\kappa}|^{-1} \f$;
      timescale_source2 = \f$ |2\ddot{a}/a-(\dot{a}/a)^2|^{-1/2} \f$ (to sample correctly the late ISW effect; and 
      timescale_source=1/(1/timescale_source1+1/timescale_source2); repeat till today.
      - if CMB not requested:
      timescale_source = 1/aH; repeat till today.
  */

  last_index_back = first_index_back; 
  last_index_thermo = first_index_thermo;
  tau = tau_ini;

  while (tau < pba->conformal_age) {

    class_call(background_at_tau(pba,
				 tau, 
				 short_info, 
				 closeby, 
				 &last_index_back, 
				 pvecback),
	       pba->error_message,
	       ppt->error_message);

    class_call(thermodynamics_at_z(pba,
				   pth,
				   1./pvecback[pba->index_bg_a]-1.,  /* redshift z=1/a-1 */
				   closeby,
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

    class_test(fabs(ppr->perturb_sampling_stepsize*timescale_source/tau) < ppr->smallest_allowed_variation,
	       ppt->error_message,
	       "integration step =%e < machine precision : leads either to numerical error or infinite loop",ppr->perturb_sampling_stepsize*timescale_source);

    tau = tau + ppr->perturb_sampling_stepsize*timescale_source; 
    counter++;

  }

  /** - infer total number of time steps, ppt->tau_size */
  ppt->tau_size = counter;

  /** - allocate array of time steps, ppt->tau_sampling[index_tau] */
  class_alloc(ppt->tau_sampling,ppt->tau_size * sizeof(double),ppt->error_message);

  /** - repeat the same steps, now filling the array with each tau value: */

  /** (a) first sampling point = when the universe stops being opaque */

  counter = 0;
  ppt->tau_sampling[counter]=tau_ini;

  /** (b) next sampling point = previous + ppr->perturb_sampling_stepsize * timescale_source, where
      timescale_source1 = \f$ |g/\dot{g}| = |\dot{\kappa}-\ddot{\kappa}/\dot{\kappa}|^{-1} \f$;
      timescale_source2 = \f$ |2\ddot{a}/a-(\dot{a}/a)^2|^{-1/2} \f$ (to smaple correctly the late ISW effect; and 
      timescale_source=1/(1/timescale_source1+1/timescale_source2); repeat till today
      - if CMB not requested:
      timescale_source = 1/aH; repeat till today.  */

  last_index_back = first_index_back; 
  last_index_thermo = first_index_thermo;
  tau = tau_ini;

  while (tau < pba->conformal_age) {
    
    class_call(background_at_tau(pba,
				 tau, 
				 short_info, 
				 closeby, 
				 &last_index_back, 
				 pvecback),
	       pba->error_message,
	       ppt->error_message);

    class_call(thermodynamics_at_z(pba,
				   pth,
				   1./pvecback[pba->index_bg_a]-1.,  /* redshift z=1/a-1 */
				   closeby,
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

    class_test(fabs(ppr->perturb_sampling_stepsize*timescale_source/tau) < ppr->smallest_allowed_variation,
	       ppt->error_message,
	       "integration step =%e < machine precision : leads either to numerical error or infinite loop",ppr->perturb_sampling_stepsize*timescale_source);

    tau = tau + ppr->perturb_sampling_stepsize*timescale_source; 
    counter++;
    ppt->tau_sampling[counter]=tau;

  }

  /** - last sampling point = exactly today */
  ppt->tau_sampling[counter] = pba->conformal_age;

  free(pvecback);
  free(pvecthermo);

  /** - loop over modes, initial conditions and types. For each of
      them, allocate array of source functions. */
  
  for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {
    for (index_ic = 0; index_ic < ppt->ic_size[index_mode]; index_ic++) {
      for (index_type = 0; index_type < ppt->tp_size[index_mode]; index_type++) {

	class_alloc(ppt->sources[index_mode][index_ic*ppt->tp_size[index_mode]+index_type],
		    ppt->k_size[index_mode] * ppt->tau_size * sizeof(double),
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

    k_rec = 2. * _PI_ / pth->rs_rec; /* comoving scale corresponding to sound horizon at recombination */

    index_k=0;
    k=ppr->k_scalar_min_tau0/pba->conformal_age;
    index_k=1;

    if (ppt->has_cls == _TRUE_) {
      k_max_cl = ppr->k_scalar_max_tau0_over_l_max
	*ppt->l_scalar_max
	/pba->conformal_age;
    }
    else
      k_max_cl = 0.;

    while (k < k_max_cl) {
      step = ppr->k_scalar_step_super 
	+ 0.5 * (tanh((k-k_rec)/k_rec/ppr->k_scalar_step_transition)+1.) * (ppr->k_scalar_step_sub-ppr->k_scalar_step_super);

      class_test(step * k_rec / k < ppr->smallest_allowed_variation,
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
    ppt->k[index_mode][index_k] = ppr->k_scalar_min_tau0/pba->conformal_age;;

    /*     printf("%d %e %g\n",index_k,ppt->k[index_mode][index_k],1.); */

    index_k++;

    while (index_k < ppt->k_size_cl[index_mode]) {
      step = ppr->k_scalar_step_super 
	+ 0.5 * (tanh((ppt->k[index_mode][index_k-1]-k_rec)/k_rec/ppr->k_scalar_step_transition)+1.) * (ppr->k_scalar_step_sub-ppr->k_scalar_step_super);

      class_test(step * k_rec / ppt->k[index_mode][index_k-1] < ppr->smallest_allowed_variation,
		 ppt->error_message,
		 "k step =%e < machine precision : leads either to numerical error or infinite loop",step * k_rec);
      ppt->k[index_mode][index_k]=ppt->k[index_mode][index_k-1] + step * k_rec;
      /*  printf("%d %e %g\n",index_k,ppt->k[index_mode][index_k],1.); */
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

    k_rec = 2. * _PI_ / pth->tau_rec; /* comoving scale corresping to causal horizon at recombination 
					 (roughly, sqrt(3) bigger than sound horizon) */

    index_k=0;
    k = ppr->k_tensor_min_tau0/pba->conformal_age;
    index_k=1;

    if (ppt->has_cls == _TRUE_) {
      k_max_cl = ppr->k_tensor_max_tau0_over_l_max
	*ppt->l_tensor_max
	/pba->conformal_age;
    }
    else
      k_max_cl = 0.;

    while (k < k_max_cl) {
      step = ppr->k_tensor_step_super 
	+ 0.5 * (tanh((k-k_rec)/k_rec/ppr->k_tensor_step_transition)+1.) * (ppr->k_tensor_step_sub-ppr->k_tensor_step_super);

      class_test(step * k_rec / k < ppr->smallest_allowed_variation,
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
    ppt->k[index_mode][index_k] = ppr->k_tensor_min_tau0/pba->conformal_age;
    index_k++;
    while (index_k < ppt->k_size_cl[index_mode]) {
      step = ppr->k_tensor_step_super 
	+ 0.5 * (tanh((ppt->k[index_mode][index_k-1]-k_rec)/k_rec/ppr->k_tensor_step_transition)+1.) * (ppr->k_tensor_step_sub-ppr->k_tensor_step_super);

      class_test(step * k_rec / ppt->k[index_mode][index_k-1] < ppr->smallest_allowed_variation,
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

  int index_mt=0;
  int index_st;
  int index_type;
  int index_ap;

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
      
    /* synchronous gauge (note that tau is counted in the vector of
       quantities to be integrated, while here we only consider
       quantities obeying to constraint equations) */

    if (ppr->gauge == synchronous) {
      ppw->index_mt_h_prime = index_mt; /* h' */
      index_mt++;
      ppw->index_mt_h_prime_prime = index_mt; /* h'' */
      index_mt++;
      ppw->index_mt_eta_prime = index_mt; /* tau' */
      index_mt++;
      ppw->index_mt_alpha_prime = index_mt; /* alpha' (with alpha = (h' + 6 tau') / (2 k**2) ) */
      index_mt++;

    }     

  }

  if ((ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors)) {

    index_mt = 0;

  }

  ppw->mt_size = index_mt;

  /** - define indices in the vector of source terms */

  index_st = 0;

  ppw->index_st_tau = index_st;
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
      source_term_table[index_type][index_tau*ppw->st_size+index_st] */

  class_alloc(ppw->source_term_table,
	      ppt->tp_size[index_mode] * sizeof(double *),
	      ppt->error_message);

  for (index_type = 0; index_type < ppt->tp_size[index_mode]; index_type++) {
    class_alloc(ppw->source_term_table[index_type], 
		ppt->tau_size*ppw->st_size*sizeof(double),
		ppt->error_message);
  }

  /** - count number of approximation, initialize their indices, and allocate their flags */
  index_ap=0;

  if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

    index_ap=0;

    ppw->index_ap_tca=index_ap;
    index_ap++;

    ppw->index_ap_rsa=index_ap;
    index_ap++;

    if (pba->has_ur == _TRUE_) {
      ppw->index_ap_ufa=index_ap;
      index_ap++;
    }

    if (pba->has_ncdm == _TRUE_) {
      ppw->index_ap_ncdmfa=index_ap;
      index_ap++;
    }
  }

  ppw->ap_size=index_ap;

  if (ppw->ap_size > 0)
    class_alloc(ppw->approx,ppw->ap_size*sizeof(int),ppt->error_message);

  /** - For definitness, initialize approximation flags to arbitrary
      values (correct values are overwritten in
      pertub_find_approximation_switches) */

  if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

    ppw->approx[ppw->index_ap_tca]=(int)tca_on;
    ppw->approx[ppw->index_ap_rsa]=(int)rsa_off;
    if (pba->has_ur == _TRUE_) {
      ppw->approx[ppw->index_ap_ufa]=(int)ufa_off;
    }
    if (pba->has_ncdm == _TRUE_) {
      ppw->approx[ppw->index_ap_ncdmfa]=(int)ncdmfa_off;
    }
  }

  /** - allocate fields where some of the perturbations are stored */

  if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

    if ((ppt->has_matter_transfers == _TRUE_) || (ppt->has_source_delta_pk == _TRUE_)) {

      class_alloc(ppw->delta_ncdm,pba->N_ncdm*sizeof(double),ppt->error_message);

    }

  }

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

  free(ppw->pvecback);
  free(ppw->pvecthermo);
  free(ppw->pvecmetric);
  for (index_type = 0; index_type < ppt->tp_size[index_mode]; index_type++) {
    free(ppw->source_term_table[index_type]);
  }
  free(ppw->source_term_table);
  if (ppw->ap_size > 0)
    free(ppw->approx);

  if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {
    if (ppt->has_matter_transfers == _TRUE_) {
      free(ppw->delta_ncdm);
    }
  }

  free(ppw);

  return _SUCCESS_;
}

/**
 * Solve the perturbation evolution for a given mode, initial
 * condition and wavenumber, and compute the corresponding source
 * functions.
 *
 * For a given mode, initial condition and wavenumber, this function
 * finds the time ranges over witch the perturbations can be described
 * within a given approximation. For each such range, it initializes
 * (or redistribute) perturbations using perturb_vector_init(), and
 * integrates over time. Whenever a "source sampling time" is passed,
 * the source terms are computed and stored temporarily (for each
 * type) using perturb_source_terms().  Finally, the actual source
 * functions are computed using the source terms, and stored in the
 * source table using perturb_sources().
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param ppt        Input/Output: pointer to the perturbation structure (output source functions S(k,tau) written here)
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

  /* contains all fixed parameters, indices and workspaces used by the perturb_derivs function */
  struct perturb_parameters_and_workspace ppaw;

  /* conformal time */
  double tau,tau_lower,tau_upper,tau_mid;

  /* maximum value of conformal time for current wavenumber */
  double taumax;

  /* index running over time */
  int index_tau;

  /* number of values in the tau_sampling array that should be considered for a given mode */
  int tau_actual_size;

  /* running index for the source term vector */
  int index_st;

  /* running index over types (temperature, etc) */
  int index_type;
  
  /* fourier mode */
  double k;

  /* number of time intervals where the approximation scheme is uniform */
  int interval_number;

  /* index running over such time intervals */
  int index_interval;

  /* number of time intervals where each particular approximation is uniform */
  int * interval_number_of;

  /* edge of intervals where approximation scheme is uniform: tau_ini, tau_switch_1, ..., tau_end */
  double * interval_limit;

  /* array of approximation scheme within each interval: interval_approx[index_interval][index_ap] */
  int ** interval_approx;

  /* index running over approximations */
  int index_ap;

  /* approximation scheme within previous interval: previous_approx[index_ap] */
  int * previous_approx;

  int n_ncdm,is_early_enough;

  /* function pointer to ODE evolver and names of possible evolvers */

  extern int evolver_rk();
  extern int evolver_ndf15(); 	
  int (*generic_evolver)();
  
  /** - initialize indices relevant for back/thermo tables search */
  ppw->last_index_back=0;
  ppw->last_index_thermo=0;
  ppw->intermode = normal;

  /** - get wavenumber value */
  k = (ppt->k[index_mode])[index_k];

  class_test(k == 0.,
	     ppt->error_message,
	     "stop to avoid division by zero");

  /** - compute maximum value of tau for which sources are calculated for this wavenumber */

  /* by default, today */
  taumax = pba->conformal_age;
  tau_actual_size = ppt->tau_size;

  /* eventually stop earlier, when k*tau=k_tau_max, but not before the end of recombination */
  /*   if (ppt->has_lss == _FALSE_) { */

  /* revisit this issue */

  /*     while (ppt->tau_sampling[tau_actual_size-1] > taumax) */
  /*       tau_actual_size--; */

  /*     class_test(tau_actual_size < 1, */
  /* 	       ppt->error_message, */
  /* 	       "did not consider this case yet"); */

  /*   } */

  /** - using bisection, compute minimum value of tau for which this
      wavenumber is integrated */

  /* will be at least the first time in the background table */  
  tau_lower = pba->tau_table[0];
    
  class_call(background_at_tau(pba,
			       tau_lower, 
			       normal_info, 
			       normal, 
			       &(ppw->last_index_back), 
			       ppw->pvecback),
	     pba->error_message,
	     ppt->error_message);
  
  class_call(thermodynamics_at_z(pba,
				 pth,
				 1./ppw->pvecback[pba->index_bg_a]-1.,
				 normal,
				 &(ppw->last_index_thermo),
				 ppw->pvecback,
				 ppw->pvecthermo),
	     pth->error_message,
	     ppt->error_message);
  
  /* check that this initial time is indeed OK given imposed
     conditions on kappa' and on k/aH */

  class_test(ppw->pvecback[pba->index_bg_a]*
	     ppw->pvecback[pba->index_bg_H]/
	     ppw->pvecthermo[pth->index_th_dkappa] >
	     ppr->start_small_k_at_tau_c_over_tau_h, ppt->error_message, "your choice of initial time for integrating wavenumbers is inappropriate: it corresponds to a time before that at which the background has been integrated. You should increase 'start_small_k_at_tau_c_over_tau_h' up to at least %g, or decrease 'a_ini_over_a_today_default'\n", 
	     ppw->pvecback[pba->index_bg_a]*
	     ppw->pvecback[pba->index_bg_H]/
	     ppw->pvecthermo[pth->index_th_dkappa]);
  
  class_test(k/ppw->pvecback[pba->index_bg_a]/ppw->pvecback[pba->index_bg_H] >
	     ppr->start_large_k_at_tau_h_over_tau_k,
             ppt->error_message,
             "your choice of initial time for integrating wavenumbers is inappropriate: it corresponds to a time before that at which the background has been integrated. You should increase 'start_large_k_at_tau_h_over_tau_k' up to at least %g, or decrease 'a_ini_over_a_today_default'\n",
	     ppt->k[index_mode][ppt->k_size[index_mode]-1]/ppw->pvecback[pba->index_bg_a]/ ppw->pvecback[pba->index_bg_H]);
  
  if (pba->has_ncdm == _TRUE_) {
    for (n_ncdm=0; n_ncdm < pba->N_ncdm; n_ncdm++) {
      class_test(fabs(ppw->pvecback[pba->index_bg_p_ncdm1+n_ncdm]/ppw->pvecback[pba->index_bg_rho_ncdm1+n_ncdm]-1./3.)>ppr->tol_ncdm_initial_w,
		 ppt->error_message,     
		 "your choice of initial time for integrating wavenumbers is inappropriate: it corresponds to a time at which the ncdm species number %d is not ultra-relativistic anymore, with w=%g, p=%g and rho=%g\n",
		 n_ncdm,
		 ppw->pvecback[pba->index_bg_p_ncdm1+n_ncdm]/ppw->pvecback[pba->index_bg_rho_ncdm1+n_ncdm],
		 ppw->pvecback[pba->index_bg_p_ncdm1+n_ncdm],
		 ppw->pvecback[pba->index_bg_rho_ncdm1+n_ncdm]);
    }
  }
  
  /* is at most the time at which sources must be sampled */
  tau_upper = ppt->tau_sampling[0];

  /* start bisection */  
  tau_mid = 0.5*(tau_lower + tau_upper);
    
  while ((tau_upper - tau_lower)/tau_lower > ppr->tol_tau_approx) {

    is_early_enough = _TRUE_;

    class_call(background_at_tau(pba,
				 tau_mid, 
				 normal_info, 
				 normal, 
				 &(ppw->last_index_back), 
				 ppw->pvecback),
	       pba->error_message,
	       ppt->error_message);

    /* if there are non-cold relics, check that they are relativistic enough */
    if (pba->has_ncdm == _TRUE_) {
      for (n_ncdm=0; n_ncdm < pba->N_ncdm; n_ncdm++) {
	if (fabs(ppw->pvecback[pba->index_bg_p_ncdm1+n_ncdm]/ppw->pvecback[pba->index_bg_rho_ncdm1+n_ncdm]-1./3.) > ppr->tol_ncdm_initial_w)
	  is_early_enough = _FALSE_;
      }
    }

    /* also check that the two conditions on (aH/kappa') and (aH/k) are fulfilled */
    if (is_early_enough == _TRUE_) {

      class_call(thermodynamics_at_z(pba,
				     pth,
				     1./ppw->pvecback[pba->index_bg_a]-1.,  /* redshift z=1/a-1 */
				     normal,
				     &(ppw->last_index_thermo),
				     ppw->pvecback,
				     ppw->pvecthermo),
		 pth->error_message,
		 ppt->error_message);
    
      if ((ppw->pvecback[pba->index_bg_a]*
	   ppw->pvecback[pba->index_bg_H]/
	   ppw->pvecthermo[pth->index_th_dkappa] >
	   ppr->start_small_k_at_tau_c_over_tau_h) ||
	  (k/ppw->pvecback[pba->index_bg_a]/ppw->pvecback[pba->index_bg_H] >
	   ppr->start_large_k_at_tau_h_over_tau_k))
	
	is_early_enough = _FALSE_;
    }

    if (is_early_enough == _TRUE_)
      tau_lower = tau_mid;
    else
      tau_upper = tau_mid;
    
    tau_mid = 0.5*(tau_lower + tau_upper);
    
  }
  
  tau = tau_mid;

  /** - find the number of intervals over which approximation scheme is constant */

  class_alloc(interval_number_of,ppw->ap_size*sizeof(int),ppt->error_message);

  ppw->intermode = normal;

  class_call(perturb_find_approximation_number(ppr,
					       pba,
					       pth,
					       ppt,
					       index_mode,
					       k,
					       ppw,
					       tau,
					       ppt->tau_sampling[tau_actual_size-1],
					       &interval_number,
					       interval_number_of),
	     ppt->error_message,
	     ppt->error_message);

  class_alloc(interval_limit,(interval_number+1)*sizeof(double),ppt->error_message);

  class_alloc(interval_approx,interval_number*sizeof(int*),ppt->error_message);

  for (index_interval=0; index_interval<interval_number; index_interval++)
    class_alloc(interval_approx[index_interval],ppw->ap_size*sizeof(int),ppt->error_message);

  class_call(perturb_find_approximation_switches(ppr,
						 pba,
						 pth,
						 ppt,
						 index_mode,
						 k,
						 ppw,
						 tau,
						 ppt->tau_sampling[tau_actual_size-1],
						 ppr->tol_tau_approx,
						 interval_number,
						 interval_number_of,
						 interval_limit,
						 interval_approx),
	     ppt->error_message,
	     ppt->error_message);

  free(interval_number_of);

  /** - fill the structure containing all fixed parameters, indices
      and workspaces needed by perturb_derivs */

  ppaw.ppr = ppr;
  ppaw.pba = pba;
  ppaw.pth = pth;
  ppaw.ppt = ppt;
  ppaw.index_mode = index_mode;
  ppaw.k = k;
  ppaw.ppw = ppw;
  ppaw.ppw->intermode = closeby;

  /** - loop over intervals over which approximatiomn scheme is uniform. For each interval: */

  for (index_interval=0; index_interval<interval_number; index_interval++) {

    /** (a) fix the approximation scheme */

    for (index_ap=0; index_ap<ppw->ap_size; index_ap++)
      ppw->approx[index_ap]=interval_approx[index_interval][index_ap];

    /** (b) get the previous approximation scheme. If the current
	interval starts from the initial time tau_ini, the previous
	approximation is set to be a NULL pointer, so that the
	function perturb_vector_init() knows that perturbations must
	be initialized */

    if (index_interval==0) {
      previous_approx=NULL;
    }
    else {
      previous_approx=interval_approx[index_interval-1];
    }

    /** (c) define the vector of perturbations to be integrated
	over. If the current interval starts from the initial time
	tau_ini, fill the vector with initial conditions for each
	mode. If it starts from an approximation switching point,
	redistribute correctly the perturbations from the previous to
	the new vector of perturbations. */

    class_call(perturb_vector_init(ppr,
				   pba,
				   pth,
				   ppt,
				   index_mode,
				   index_ic,
				   k,
				   interval_limit[index_interval],
				   ppw,
				   previous_approx),
	       ppt->error_message,
	       ppt->error_message);

    /** (d) integrate the perturbations over the current interval. */

    if(ppr->evolver == rk){
      generic_evolver = evolver_rk;
    }
    else{
      generic_evolver = evolver_ndf15;
    }

    class_call(generic_evolver(perturb_derivs,
			       interval_limit[index_interval],
			       interval_limit[index_interval+1],
			       ppw->pv->y,
			       ppw->pv->used_in_sources,
			       ppw->pv->pt_size,
			       &ppaw,
			       ppr->tol_perturb_integration,
			       ppr->smallest_allowed_variation,
			       perturb_timescale,
			       ppr->perturb_integration_stepsize,
			       ppt->tau_sampling,
			       tau_actual_size,
			       perturb_source_terms,
			       /* next entry = 'NULL' in standard
				  runs; = 'perturb_print_variables' if
				  you want to output perturbations at
				  each integration step for testing
				  purposes (you'll get many more
				  points if you use the runge-kutta
				  integrator, set evolver=rk) */
			       NULL,
			       //perturb_print_variables,
			       ppt->error_message),
	       ppt->error_message,
	       ppt->error_message);
          
  } 

  /** fill the source terms array with zeros for all times between
      then last integrated time tau_max and tau_today. */

  for (index_tau = tau_actual_size; index_tau < ppt->tau_size; index_tau++) {
    for (index_type = 0; index_type < ppt->tp_size[index_mode]; index_type++) {
      for (index_st = 0; index_st < ppw->st_size; index_st++) {
	ppw->source_term_table[index_type][index_tau*ppw->st_size+index_st] = 0.;
      }
    }
  }

  /** - free quantitites allocated at the beginning of the routine */

  class_call(perturb_vector_free(ppw->pv),
	     ppt->error_message,
	     ppt->error_message);

  for (index_interval=0; index_interval<interval_number; index_interval++)
    free(interval_approx[index_interval]);

  free(interval_approx);

  free(interval_limit);

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
 * For a given mode and wavenumber, find the number of interval of
 * times bewteen tau_ini and tau_end such that the approximation
 * scheme (and the number of perturbation equations) is uniform.
 *
 * @param ppr                Input: pointer to precision structure
 * @param pba                Input: pointer to background structure
 * @param pth                Input: pointer to the thermodynamics structure
 * @param ppt                Input: pointer to the perturbation structure
 * @param index_mode         Input: index of mode under consideration (scalar/.../tensor)
 * @param index_k            Input: index of wavenumber
 * @param ppw                Input: pointer to perturb_workspace structure containing index values and workspaces
 * @param tau_ini            Input: initial time of the perturbation integration
 * @param tau_end            Input: final time of the perturbation integration
 * @param interval_number    Output: total number of intervals
 * @param interval_number_of Output: number of intervals with respect to each particular approximation
 * @return the error status
 */

int perturb_find_approximation_number(
				      struct precision * ppr,
				      struct background * pba,
				      struct thermo * pth,
				      struct perturbs * ppt,
				      int index_mode,
				      double k,
				      struct perturb_workspace * ppw,
				      double tau_ini,
				      double tau_end,
				      int * interval_number,
				      int * interval_number_of /* interval_number_of[index_ap] (already allocated) */
				      ){
  
  /* index running over approximations */
  int index_ap;

  /* value of a given approximation at tau_ini and tau_end */
  int flag_ini,flag_end;

  /** - fix default number of intervals to one (if no approximation switch) */  

  *interval_number=1; 

  /** - loop over each approximation and add the number of approximation switching times */

  for (index_ap=0; index_ap<ppw->ap_size; index_ap++) {

    class_call(perturb_approximations(ppr,
				      pba,
				      pth,
				      ppt,
				      index_mode,
				      k,
				      tau_ini,
				      ppw),
	       ppt->error_message,
	       ppt->error_message);
    
    flag_ini = ppw->approx[index_ap];
    
    class_call(perturb_approximations(ppr,
				      pba,
				      pth,
				      ppt,
				      index_mode,
				      k,
				      tau_end,
				      ppw),
	       ppt->error_message,
	       ppt->error_message);
    
    flag_end = ppw->approx[index_ap];
    
    class_test(flag_end<flag_ini,
	       ppt->error_message,
	       "For each approximation scheme, the declaration of approximation labels in the enumeration must follow chronological order, e.g: enum approx_flags {flag1, flag2, flag3} with flag1 being the initial one and flag3 the final one");
    
    *interval_number += flag_end-flag_ini;

    interval_number_of[index_ap] = flag_end-flag_ini+1;
  }
  
  return _SUCCESS_;

}

/**
 * For a given mode and wavenumber, find the values of time at which
 * the approximation changes.
 *
 * @param ppr                Input: pointer to precision structure
 * @param pba                Input: pointer to background structure
 * @param pth                Input: pointer to the thermodynamics structure
 * @param ppt                Input: pointer to the perturbation structure
 * @param index_mode         Input: index of mode under consideration (scalar/.../tensor)
 * @param index_k            Input: index of wavenumber
 * @param ppw                Input: pointer to perturb_workspace structure containing index values and workspaces
 * @param tau_ini            Input: initial time of the perturbation integration
 * @param tau_end            Input: final time of the perturbation integration
 * @param interval_number    Input: total number of intervals
 * @param interval_number_of Input: number of intervals with respect to each particular approximation
 * @param interval_limit     Output: value of time at the boundary of the intervals: tau_ini, tau_switch1, ..., tau_end 
 * @param interval_approx    Output: value of approximations in each interval
 * @return the error status
 */

int perturb_find_approximation_switches(
					struct precision * ppr,
					struct background * pba,
					struct thermo * pth,
					struct perturbs * ppt,
					int index_mode,
					double k,
					struct perturb_workspace * ppw,
					double tau_ini,
					double tau_end,
					double precision,
					int interval_number,
					int * interval_number_of,
					double * interval_limit, /* interval_limit[index_interval] (already allocated) */
					int ** interval_approx   /* interval_approx[index_interval][index_ap] (already allocated) */
					){

  int index_ap;
  int index_switch;
  int index_switch_tot;
  int num_switch; 
  double tau_min,lower_bound,upper_bound;
  double mid=0;
  double * unsorted_tau_switch;
  double next_tau_switch;
  int flag_ini;

  /** - write in output arrays the initial time and approximation */

  interval_limit[0]=tau_ini;

  class_call(perturb_approximations(ppr,
				    pba,
				    pth,
				    ppt,
				    index_mode,
				    k,
				    tau_ini,
				    ppw),
	     ppt->error_message,
	     ppt->error_message);

  for (index_ap=0; index_ap<ppw->ap_size; index_ap++)
    interval_approx[0][index_ap]=ppw->approx[index_ap];
  
  /** - if there are no approximation switches, just write final time and return */

  if (interval_number == 1) {

    interval_limit[1]=tau_end;

  }

  /** - if there are switches, consider approximations one after each
      other.  Find switching time by bisection. Store all switches in
      arbitrary order in array unsorted_tau_switch[] */

  else {

    class_alloc(unsorted_tau_switch,(interval_number-1)*sizeof(double),ppt->error_message);

    index_switch_tot=0;

    for (index_ap=0; index_ap<ppw->ap_size; index_ap++) {

      if (interval_number_of[index_ap] > 1) {
	
	num_switch = interval_number_of[index_ap]-1;

	tau_min = tau_ini;
	
	flag_ini = interval_approx[0][index_ap];
	
	for (index_switch=0; index_switch<num_switch; index_switch++) {
	  
	  lower_bound=tau_min;
	  upper_bound=tau_end;
	  mid = 0.5*(lower_bound+upper_bound);

	  while (upper_bound - lower_bound > precision) {
	  	    
	    class_call(perturb_approximations(ppr,
					      pba,
					      pth,
					      ppt,
					      index_mode,
					      k,
					      mid,
					      ppw),
		       ppt->error_message,
		       ppt->error_message);

	    if (ppw->approx[index_ap] > flag_ini+index_switch) {
	      upper_bound=mid;
	    }
	    else {
	      lower_bound=mid;
	    }

	    mid = 0.5*(lower_bound+upper_bound);

	  }

	  unsorted_tau_switch[index_switch_tot]=mid;
	  index_switch_tot++;

	  tau_min=mid;

	}
      }
    }

    class_test(index_switch_tot != (interval_number-1),
	       ppt->error_message,
	       "bug in approximation switch search routine: should have %d = %d",
	       index_switch_tot,interval_number-1);
    
    /** - now sort interval limits in correct order */
    
    index_switch_tot=1;
    
    while (index_switch_tot < interval_number) {
      
      next_tau_switch=tau_end;
      for (index_switch=0; index_switch<interval_number-1; index_switch++) {
	if ((unsorted_tau_switch[index_switch] > interval_limit[index_switch_tot-1]) &&
	    (unsorted_tau_switch[index_switch] < next_tau_switch)) {
	  next_tau_switch=unsorted_tau_switch[index_switch];
	}
      }
      interval_limit[index_switch_tot]=next_tau_switch;
      index_switch_tot++;
    }
    
    interval_limit[index_switch_tot]=tau_end;
    
    class_test(index_switch_tot != interval_number,
	       ppt->error_message,
	       "most probably two approximation switching time were found to be equal, which cannot be handled\n");
    
    /** - store each approximation in chronological order */

    for (index_switch=1; index_switch<interval_number; index_switch++) {
      
      class_call(perturb_approximations(ppr,
					pba,
					pth,
					ppt,
					index_mode,
					k,
					0.5*(interval_limit[index_switch]+interval_limit[index_switch+1]),
					ppw),

		 ppt->error_message,
		 ppt->error_message);
      
      for (index_ap=0; index_ap<ppw->ap_size; index_ap++)
	interval_approx[index_switch][index_ap]=ppw->approx[index_ap];

      if (ppt->perturbations_verbose>2) {

	if ((interval_approx[index_switch-1][ppw->index_ap_tca]==(int)tca_on) && 
	    (interval_approx[index_switch][ppw->index_ap_tca]==(int)tca_off))
	  fprintf(stdout,"Mode k=%e: will switch off tight-coupling approximation at tau=%e\n",k,interval_limit[index_switch]);

	if ((interval_approx[index_switch-1][ppw->index_ap_rsa]==(int)rsa_off) && 
	    (interval_approx[index_switch][ppw->index_ap_rsa]==(int)rsa_on))
	  fprintf(stdout,"Mode k=%e: will switch on radiation streaming approximation at tau=%e\n",k,interval_limit[index_switch]);

	if (pba->has_ur == _TRUE_) {
	  if ((interval_approx[index_switch-1][ppw->index_ap_ufa]==(int)ufa_off) && 
	      (interval_approx[index_switch][ppw->index_ap_ufa]==(int)ufa_on)) {
	    fprintf(stdout,"Mode k=%e: will switch on ur fluid approximation at tau=%e\n",k,interval_limit[index_switch]);
	  }
	}
	if (pba->has_ncdm == _TRUE_) {
	  if ((interval_approx[index_switch-1][ppw->index_ap_ncdmfa]==(int)ncdmfa_off) && 
	      (interval_approx[index_switch][ppw->index_ap_ncdmfa]==(int)ncdmfa_on)) {
	    fprintf(stdout,"Mode k=%e: will switch on ncdm fluid approximation at tau=%e\n",k,interval_limit[index_switch]);
	  }
	}
      }
    }
    
    free(unsorted_tau_switch);

  }

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
 * @param tau        Input: conformal time
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
			double tau,
			struct perturb_workspace * ppw, /* ppw->pv unallocated if pa_old = NULL, allocated and filled otherwise */
			int * pa_old 
			) {

  /** Summary: */

  /** - define local variables */

  struct perturb_vector * ppv;

  int index_pt;
  int l;
  int n_ncdm,index_q,ncdm_l_size;
  double rho_plus_p_ncdm,q,q2,epsilon,a,factor;

  /** - allocate a new perturb_vector structure to which ppw->pv will point at the end of the routine */

  class_alloc(ppv,sizeof(struct perturb_vector),ppt->error_message);

  /** - initialize pointers to NULL (they will be allocated later if
    needed), relevant for perturb_vector_free() */
  ppv->l_max_ncdm = NULL;
  ppv->q_size_ncdm = NULL;

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

    if (ppw->approx[ppw->index_ap_rsa] == (int)rsa_off) { /* if radiation streaming approximation is off */

      /* photons */

      ppv->index_pt_delta_g = index_pt; /* photon density */
      index_pt++;

      ppv->index_pt_theta_g = index_pt; /* photon velocity */
      index_pt++;

      if (ppw->approx[ppw->index_ap_tca] == (int)tca_off) { /* if tight-coupling approximation is off */

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

    /* fluid */    

    if (pba->has_fld == _TRUE_) {       
      
      ppv->index_pt_delta_fld = index_pt; /* fluid density */
      index_pt++;

      ppv->index_pt_theta_fld = index_pt; /* fluid velocity */
      index_pt++;
      
    }

    /* ultra relativistic neutrinos */

    if (pba->has_ur == _TRUE_) {

      /* reject inconsistent values of the number of mutipoles in ultra relativistic neutrino hierachy */
      class_test(ppr->l_max_ur < 4,
		 ppt->error_message,
		 "ppr->l_max_ur should be at least 4, i.e. we must integrate at least over neutrino/relic density, velocity, shear, third and fourth momentum");
      
      if (ppw->approx[ppw->index_ap_rsa] == (int)rsa_off) { /* if radiation streaming approximation is off */

	ppv->index_pt_delta_ur = index_pt; /* density of ultra-relativistic neutrinos/relics */
	index_pt++;
	
	ppv->index_pt_theta_ur = index_pt; /* velocity of ultra-relativistic neutrinos/relics */
	index_pt++;
	
	ppv->index_pt_shear_ur = index_pt; /* shear of ultra-relativistic neutrinos/relics */
	index_pt++;

	if (ppw->approx[ppw->index_ap_ufa] == (int)ufa_off) { /* if neutrino free-streaming approximation is off */
	
	  ppv->index_pt_l3_ur = index_pt; /* l=3 of ultra-relativistic neutrinos/relics */
	  index_pt++;
	
	  ppv->l_max_ur = ppr->l_max_ur; /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3) */
	  index_pt += (ppv->l_max_ur-3);
	}
      }
      
    }

    /* non-cold dark matter */

    if (pba->has_ncdm == _TRUE_) {
      ppv->index_pt_psi0_ncdm1 = index_pt; /* density of ultra-relativistic neutrinos/relics */
      ppv->N_ncdm = pba->N_ncdm;
      class_alloc(ppv->l_max_ncdm,ppv->N_ncdm*sizeof(double),ppt->error_message);
      class_alloc(ppv->q_size_ncdm,ppv->N_ncdm*sizeof(double),ppt->error_message);
            
      for(n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++){
	// Set value of ppv->l_max_ncdm:
	if(ppw->approx[ppw->index_ap_ncdmfa] == (int)ncdmfa_off){
	  /* reject inconsistent values of the number of mutipoles in ultra relativistic neutrino hierachy */
	  class_test(ppr->l_max_ncdm < 4,
		     ppt->error_message,
		     "ppr->l_max_ncdm=%d should be at least 4, i.e. we must integrate at least over first four momenta of non-cold dark matter perturbed phase-space distribution",n_ncdm);
	  //Copy value from precision parameter:
	  ppv->l_max_ncdm[n_ncdm] = ppr->l_max_ncdm;
	  ppv->q_size_ncdm[n_ncdm] = pba->q_size_ncdm[n_ncdm];
	}
	else{
	  // In the fluid approximaation, hierarcy is cut at lmax = 2 and q dependance is integrated out:
	  ppv->l_max_ncdm[n_ncdm] = 2;
	  ppv->q_size_ncdm[n_ncdm] = 1;
	}
	index_pt += (ppv->l_max_ncdm[n_ncdm]+1)*ppv->q_size_ncdm[n_ncdm];
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

  class_calloc(ppv->y,ppv->pt_size,sizeof(double),ppt->error_message);
  class_alloc(ppv->dy,ppv->pt_size*sizeof(double),ppt->error_message);
  class_alloc(ppv->used_in_sources,ppv->pt_size*sizeof(int),ppt->error_message);

  /** - specify which perturbations are needed in the evaluation of source terms */

  /* take all of them by default */
  for (index_pt=0; index_pt<ppv->pt_size; index_pt++)
    ppv->used_in_sources[index_pt] = _TRUE_;

  /* indicate which ones are not needed (this is just for saving time,
     omitting perturbations in this list will not change the
     results!) */

  if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

    if (ppw->approx[ppw->index_ap_rsa] == (int)rsa_off) {
      
      if (ppw->approx[ppw->index_ap_tca] == (int)tca_off) {
	
	/* we don't need temperature multipoles above l=2 (but they are
	   defined only when rsa and tca are off) */
	
	for (index_pt=ppv->index_pt_l3_g; index_pt <= ppv->index_pt_delta_g+ppv->l_max_g; index_pt++)
	  ppv->used_in_sources[index_pt]=_FALSE_;
	
	/* for polarisation, we only need l=0,2 (but l =1,3, ... are
	   defined only when rsa and tca are off) */
	
	ppv->used_in_sources[ppv->index_pt_pol1_g]=_FALSE_;
	
	for (index_pt=ppv->index_pt_pol3_g; index_pt <= ppv->index_pt_pol0_g+ppv->l_max_pol_g; index_pt++)
	  ppv->used_in_sources[index_pt]=_FALSE_;
	
      }
      
    }

    if (pba->has_ur == _TRUE_) {
      
      if (ppw->approx[ppw->index_ap_rsa] == (int)rsa_off) {
	
	if (ppw->approx[ppw->index_ap_ufa] == (int)ufa_off) {
	  
	  /* we don't need ur multipoles above l=2 (but they are
	     defined only when rsa and ufa are off) */
	  
	  for (index_pt=ppv->index_pt_l3_ur; index_pt <= ppv->index_pt_delta_ur+ppv->l_max_ur; index_pt++)
	    ppv->used_in_sources[index_pt]=_FALSE_;
	  
	}
      }
    }
    
    if (pba->has_ncdm == _TRUE_) {
      
      /* we don't need ncdm multipoles above l=2 (but they are
	 defined only when ncdmfa is off) */
      
      index_pt = ppv->index_pt_psi0_ncdm1;
      for(n_ncdm = 0; n_ncdm < ppv-> N_ncdm; n_ncdm++){
	for(index_q=0; index_q < ppv->q_size_ncdm[n_ncdm]; index_q++){ 
	  for(l=0; l<=ppv->l_max_ncdm[n_ncdm]; l++){
	    if (l>2) ppv->used_in_sources[index_pt]=_FALSE_;
	    index_pt++;
	  }
	}
      }
    }
  }

  if ((ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors)) {

    /* we don't need temperature multipoles above except l=0,2,4 */
	
    ppv->used_in_sources[ppv->index_pt_theta_g]=_FALSE_;
    ppv->used_in_sources[ppv->index_pt_l3_g]=_FALSE_;

    for (index_pt=ppv->index_pt_delta_g+5; index_pt <= ppv->index_pt_delta_g+ppv->l_max_g; index_pt++) 
      ppv->used_in_sources[index_pt]=_FALSE_;
	
    /* same for polarisation, we only need l=0,2,4 */
	
    ppv->used_in_sources[ppv->index_pt_pol1_g]=_FALSE_;
    ppv->used_in_sources[ppv->index_pt_pol3_g]=_FALSE_;
	
    for (index_pt=ppv->index_pt_pol0_g+5; index_pt <= ppv->index_pt_pol0_g+ppv->l_max_pol_g; index_pt++)
      ppv->used_in_sources[index_pt]=_FALSE_;
    
  }

  /** - case of setting initial conditions for a new wavenumber */

  if (pa_old == NULL) {

    if (ppt->perturbations_verbose>2)
      fprintf(stdout,"Mode k=%e: initializing vector at tau=%e\n",k,tau);

    if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

      /** (a) check that current approximation scheme is consistent
	  with initial conditions */

      class_test(ppw->approx[ppw->index_ap_rsa] == (int)rsa_on,
		 ppt->error_message,
		 "scalar initial conditions assume radiation streaming approximation turned off");
      

      if (pba->has_ur == _TRUE_) {

	class_test(ppw->approx[ppw->index_ap_ufa] == (int)ufa_on,
		   ppt->error_message,
		   "scalar initial conditions assume ur fluid approximation turned off");
	
      }
      
      if (pba->has_ncdm == _TRUE_) {

	class_test(ppw->approx[ppw->index_ap_ncdmfa] == (int)ncdmfa_on,
		   ppt->error_message,
		   "scalar initial conditions assume ncdm fluid approximation turned off");
	
      }
      
      class_test(ppw->approx[ppw->index_ap_tca] == (int)tca_off,
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
					  tau,
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

      class_test((pa_old[ppw->index_ap_tca] == (int)tca_off) && (ppw->approx[ppw->index_ap_tca] == (int)tca_on),
		 ppt->error_message,
		 "at tau=%g: the tight-coupling approximation can be switched off, not on",tau);

      /** -- some variables (b, cdm, fld, ...) are not affected by
             any approximation. They need to be reconducted whatever
             the approximation switching is. We treat them here. Below
             we will treat other variables case by case. */

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
      
      if (pba->has_fld == _TRUE_) {  
	
	ppv->y[ppv->index_pt_delta_fld] =
	  ppw->pv->y[ppw->pv->index_pt_delta_fld];
	
	ppv->y[ppv->index_pt_theta_fld] =
	  ppw->pv->y[ppw->pv->index_pt_theta_fld];
      }
      
      if (ppr->gauge == synchronous)
	ppv->y[ppv->index_pt_eta] =
	  ppw->pv->y[ppw->pv->index_pt_eta];

      /* -- case of switching off tight coupling
	 approximation. Provide correct initial conditions to new set
	 of variables */

      if ((pa_old[ppw->index_ap_tca] == (int)tca_on) && (ppw->approx[ppw->index_ap_tca] == (int)tca_off)) {

	if (ppt->perturbations_verbose>2)
	  fprintf(stdout,"Mode k=%e: switch off tight-coupling approximation at tau=%e\n",k,tau);
	
	ppv->y[ppv->index_pt_delta_g] =
	  ppw->pv->y[ppw->pv->index_pt_delta_g];

	ppv->y[ppv->index_pt_theta_g] =
	  ppw->pv->y[ppw->pv->index_pt_theta_g];
	  
	/* tight-coupling approximation for shear_g (previously
	   computed in perturb_derivs: perturb_derivs is always
	   called at the end of generic_evolver, in order to update
	   all quantities in ppw to the time at which the
	   approximation is switched off) */
	ppv->y[ppv->index_pt_shear_g] = ppw->tca_shear_g;
	
	ppv->y[ppv->index_pt_l3_g] = 6./7.*k/ppw->pvecthermo[pth->index_th_dkappa]*
	  ppv->y[ppv->index_pt_shear_g]; /* second-order tight-coupling approximation for l=3 */

	ppv->y[ppv->index_pt_pol0_g] = 2.5*ppv->y[ppv->index_pt_shear_g]; /* first-order tight-coupling approximation for polarization, l=0 */
	ppv->y[ppv->index_pt_pol1_g] = 7./12.*ppv->y[ppv->index_pt_l3_g]; /* second-order tight-coupling approximation for polarization, l=1 */
	ppv->y[ppv->index_pt_pol2_g] = 0.5*ppv->y[ppv->index_pt_shear_g]; /* first-order tight-coupling approximation for polarization, l=2 */
	ppv->y[ppv->index_pt_pol3_g] = 0.25*ppv->y[ppv->index_pt_l3_g];   /* second-order tight-coupling approximation for polarization, l=3 */
	
	if (pba->has_ur == _TRUE_) {

	  ppv->y[ppv->index_pt_delta_ur] =
	    ppw->pv->y[ppw->pv->index_pt_delta_ur];
	    
	  ppv->y[ppv->index_pt_theta_ur] =
	    ppw->pv->y[ppw->pv->index_pt_theta_ur];

	  ppv->y[ppv->index_pt_shear_ur] =
	    ppw->pv->y[ppw->pv->index_pt_shear_ur];
	    
	  if (ppw->approx[ppw->index_ap_ufa] == (int)ufa_off) {

	    ppv->y[ppv->index_pt_l3_ur] =
	      ppw->pv->y[ppw->pv->index_pt_l3_ur];

	    for (l=4; l <= ppv->l_max_ur; l++)
	      ppv->y[ppv->index_pt_delta_ur+l] = 
		ppw->pv->y[ppw->pv->index_pt_delta_ur+l];

	  }
	}

	if (pba->has_ncdm == _TRUE_) {
	  index_pt = 0;
	  for(n_ncdm = 0; n_ncdm < ppv->N_ncdm; n_ncdm++){
	    for(index_q=0; index_q < ppv->q_size_ncdm[n_ncdm]; index_q++){
	      for(l=0; l<=ppv->l_max_ncdm[n_ncdm];l++){
		// This is correct with or without ncdmfa, since ppv->lmax_ncdm is set accordingly.
		ppv->y[ppv->index_pt_psi0_ncdm1+index_pt] = 
		  ppw->pv->y[ppw->pv->index_pt_psi0_ncdm1+index_pt];
		index_pt++;
	      }
	    }
	  }
	}
      }

      /* -- case of switching on radiation streaming
	 approximation. Provide correct initial conditions to new set
	 of variables */

      if ((pa_old[ppw->index_ap_rsa] == (int)rsa_off) && (ppw->approx[ppw->index_ap_rsa] == (int)rsa_on)) {

	if (ppt->perturbations_verbose>2)
	  fprintf(stdout,"Mode k=%e: switch on radiation streaming approximation at tau=%e with Omega_r=%g\n",k,tau,ppw->pvecback[pba->index_bg_Omega_r]);

	if (pba->has_ncdm == _TRUE_) {
	  index_pt = 0;
	  for(n_ncdm = 0; n_ncdm < ppv->N_ncdm; n_ncdm++){
	    for(index_q=0; index_q < ppv->q_size_ncdm[n_ncdm]; index_q++){ 
	      for(l=0; l<=ppv->l_max_ncdm[n_ncdm]; l++){
		ppv->y[ppv->index_pt_psi0_ncdm1+index_pt] = 
		  ppw->pv->y[ppw->pv->index_pt_psi0_ncdm1+index_pt];
		index_pt++;
	      }
	    }
	  }
	}	
      }

      /* -- case of switching on ur fluid
	 approximation. Provide correct initial conditions to new set
	 of variables */

      if (pba->has_ur == _TRUE_) {

	if ((pa_old[ppw->index_ap_ufa] == (int)ufa_off) && (ppw->approx[ppw->index_ap_ufa] == (int)ufa_on)) {

	  if (ppt->perturbations_verbose>2)
	    fprintf(stdout,"Mode k=%e: switch on ur fluid approximation at tau=%e\n",k,tau);

	  if (ppw->approx[ppw->index_ap_rsa] == (int)rsa_off) {

	    ppv->y[ppv->index_pt_delta_g] =
	      ppw->pv->y[ppw->pv->index_pt_delta_g];

	    ppv->y[ppv->index_pt_theta_g] =
	      ppw->pv->y[ppw->pv->index_pt_theta_g];
	  }

	  if ((ppw->approx[ppw->index_ap_tca] == (int)tca_off) && (ppw->approx[ppw->index_ap_rsa] == (int)rsa_off)) {

	    ppv->y[ppv->index_pt_shear_g] =
	      ppw->pv->y[ppw->pv->index_pt_shear_g];

	    ppv->y[ppv->index_pt_l3_g] =
	      ppw->pv->y[ppw->pv->index_pt_l3_g];

	    for (l = 4; l <= ppw->pv->l_max_g; l++) { 

	      ppv->y[ppv->index_pt_delta_g+l] =
		ppw->pv->y[ppw->pv->index_pt_delta_g+l];
	    }
	  	  
	    ppv->y[ppv->index_pt_pol0_g] =
	      ppw->pv->y[ppw->pv->index_pt_pol0_g];

	    ppv->y[ppv->index_pt_pol1_g] =
	      ppw->pv->y[ppw->pv->index_pt_pol1_g];

	    ppv->y[ppv->index_pt_pol2_g] =
	      ppw->pv->y[ppw->pv->index_pt_pol2_g];

	    ppv->y[ppv->index_pt_pol3_g] =
	      ppw->pv->y[ppw->pv->index_pt_pol3_g];

	    for (l = 4; l <= ppw->pv->l_max_pol_g; l++) { 

	      ppv->y[ppv->index_pt_pol0_g+l] =
		ppw->pv->y[ppw->pv->index_pt_pol0_g+l];
	    }

	  }

	  if (ppw->approx[ppw->index_ap_rsa] == (int)rsa_off) {

	    ppv->y[ppv->index_pt_delta_ur] =
	      ppw->pv->y[ppw->pv->index_pt_delta_ur];
	    
	    ppv->y[ppv->index_pt_theta_ur] =
	      ppw->pv->y[ppw->pv->index_pt_theta_ur];
	    
	    ppv->y[ppv->index_pt_shear_ur] =
	      ppw->pv->y[ppw->pv->index_pt_shear_ur];
	  }

	  if (pba->has_ncdm == _TRUE_) {
	    index_pt = 0;
	    for(n_ncdm = 0; n_ncdm < ppv->N_ncdm; n_ncdm++){
	      for(index_q=0; index_q < ppv->q_size_ncdm[n_ncdm]; index_q++){ 
		for(l=0; l<=ppv->l_max_ncdm[n_ncdm]; l++){
		  /** This is correct even when ncdmfa == off, since ppv->l_max_ncdm and 
		      ppv->q_size_ncdm is updated.*/
		  ppv->y[ppv->index_pt_psi0_ncdm1+index_pt] = 
		    ppw->pv->y[ppw->pv->index_pt_psi0_ncdm1+index_pt];
		  index_pt++;
		}
	      }
	    }
	  }
	}
      }

      /* -- case of switching on ncdm fluid
	 approximation. Provide correct initial conditions to new set
	 of variables */

      if (pba->has_ncdm == _TRUE_) {

	if ((pa_old[ppw->index_ap_ncdmfa] == (int)ncdmfa_off) && (ppw->approx[ppw->index_ap_ncdmfa] == (int)ncdmfa_on)) {
	  
	  if (ppt->perturbations_verbose>2)
	    fprintf(stdout,"Mode k=%e: switch on ncdm fluid approximation at tau=%e\n",k,tau);
	  
	  if (ppw->approx[ppw->index_ap_rsa] == (int)rsa_off) {
	    
	    ppv->y[ppv->index_pt_delta_g] =
	      ppw->pv->y[ppw->pv->index_pt_delta_g];
	    
	    ppv->y[ppv->index_pt_theta_g] =
	      ppw->pv->y[ppw->pv->index_pt_theta_g];
	  }
	  
	  if ((ppw->approx[ppw->index_ap_tca] == (int)tca_off) && (ppw->approx[ppw->index_ap_rsa] == (int)rsa_off)) {
	    
	    ppv->y[ppv->index_pt_shear_g] =
	      ppw->pv->y[ppw->pv->index_pt_shear_g];
	    
	    ppv->y[ppv->index_pt_l3_g] =
	      ppw->pv->y[ppw->pv->index_pt_l3_g];
	    
	    for (l = 4; l <= ppw->pv->l_max_g; l++) { 
	      
	      ppv->y[ppv->index_pt_delta_g+l] =
		ppw->pv->y[ppw->pv->index_pt_delta_g+l];
	    }
	    
	    ppv->y[ppv->index_pt_pol0_g] =
	      ppw->pv->y[ppw->pv->index_pt_pol0_g];

	    ppv->y[ppv->index_pt_pol1_g] =
	      ppw->pv->y[ppw->pv->index_pt_pol1_g];
	    
	    ppv->y[ppv->index_pt_pol2_g] =
	      ppw->pv->y[ppw->pv->index_pt_pol2_g];
	    
	    ppv->y[ppv->index_pt_pol3_g] =
	      ppw->pv->y[ppw->pv->index_pt_pol3_g];
	    
	    for (l = 4; l <= ppw->pv->l_max_pol_g; l++) { 
	      
	      ppv->y[ppv->index_pt_pol0_g+l] =
		ppw->pv->y[ppw->pv->index_pt_pol0_g+l];
	    }
	    
	  }
	  
	  if (pba->has_ur == _TRUE_) {

	    if (ppw->approx[ppw->index_ap_rsa] == (int)rsa_off) {
	      
	      
	      ppv->y[ppv->index_pt_delta_ur] =
		ppw->pv->y[ppw->pv->index_pt_delta_ur];
	      
	      ppv->y[ppv->index_pt_theta_ur] =
		ppw->pv->y[ppw->pv->index_pt_theta_ur];
	      
	      ppv->y[ppv->index_pt_shear_ur] =
		ppw->pv->y[ppw->pv->index_pt_shear_ur];
	      
	      if (ppw->approx[ppw->index_ap_ufa] == (int)ufa_off) {
		
		ppv->y[ppv->index_pt_l3_ur] =
		  ppw->pv->y[ppw->pv->index_pt_l3_ur];
		
		for (l=4; l <= ppv->l_max_ur; l++)
		  ppv->y[ppv->index_pt_delta_ur+l] = 
		    ppw->pv->y[ppw->pv->index_pt_delta_ur+l];
		
	      }
	    }
	  }
	  
	  a = ppw->pvecback[pba->index_bg_a];
	  index_pt = ppw->pv->index_pt_psi0_ncdm1;
	  for(n_ncdm = 0; n_ncdm < ppv->N_ncdm; n_ncdm++){
	    // We are in the fluid approximation, so ncdm_l_size is always 3.
	    ncdm_l_size = ppv->l_max_ncdm[n_ncdm]+1;
	    rho_plus_p_ncdm = ppw->pvecback[pba->index_bg_rho_ncdm1+n_ncdm]+
	      ppw->pvecback[pba->index_bg_p_ncdm1+n_ncdm];
	    for(l=0; l<=2; l++){
	      ppv->y[ppv->index_pt_psi0_ncdm1+ncdm_l_size*n_ncdm+l] = 0.0;
	    }
	    factor = pba->factor_ncdm[n_ncdm]*pow(pba->a_today/a,4);
	    for(index_q=0; index_q < ppw->pv->q_size_ncdm[n_ncdm]; index_q++){
	      // Integrate over distributions:
	      q = pba->q_ncdm[n_ncdm][index_q];
	      q2 = q*q;
	      epsilon = sqrt(q2+a*a*pba->M_ncdm[n_ncdm]*pba->M_ncdm[n_ncdm]);
	      ppv->y[ppv->index_pt_psi0_ncdm1+ncdm_l_size*n_ncdm] +=
		pba->w_ncdm[n_ncdm][index_q]*q2*epsilon*
		ppw->pv->y[index_pt];
	      
	      ppv->y[ppv->index_pt_psi0_ncdm1+ncdm_l_size*n_ncdm+1] +=
		pba->w_ncdm[n_ncdm][index_q]*q2*q*  
		ppw->pv->y[index_pt+1];
	      
	      ppv->y[ppv->index_pt_psi0_ncdm1+ncdm_l_size*n_ncdm+2] +=
		pba->w_ncdm[n_ncdm][index_q]*q2*q2/epsilon*
		ppw->pv->y[index_pt+2];
	      
	      //Jump to next momentum bin in ppw->pv->y:
	      index_pt += (ppw->pv->l_max_ncdm[n_ncdm]+1);
	    }
	    ppv->y[ppv->index_pt_psi0_ncdm1+ncdm_l_size*n_ncdm] *=factor/ppw->pvecback[pba->index_bg_rho_ncdm1+n_ncdm];
	    ppv->y[ppv->index_pt_psi0_ncdm1+ncdm_l_size*n_ncdm+1] *=k*factor/rho_plus_p_ncdm;
	    ppv->y[ppv->index_pt_psi0_ncdm1+ncdm_l_size*n_ncdm+2] *=2.0/3.0*factor/rho_plus_p_ncdm;
	  }
	}
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

  if (pv->l_max_ncdm != NULL) free(pv->l_max_ncdm);
  if (pv->q_size_ncdm != NULL) free(pv->q_size_ncdm);
  free(pv->y);
  free(pv->dy);
  free(pv->used_in_sources);
  free(pv);
  
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
 * @param tau        Input: conformal time
 * @param ppw        Input/Output: workspace containing input the approximation scheme, the background/thermodynamics/metric quantitites, and eventually the previous vector y; and in output the new vector y.
 * @return the error status
 */
int perturb_initial_conditions(struct precision * ppr,
			       struct background * pba,
			       struct perturbs * ppt,
			       int index_mode,
			       int index_ic,
			       double k,
			       double tau,
			       struct perturb_workspace * ppw
			       ) {
  /** Summary: */

  /** - assuming that everything has already been set to zero, write non-zero initial conditions: */

  double a;
  double delta_ur,theta_ur,shear_ur,l3_ur;
  double q,epsilon;
  int index_q,n_ncdm,idx;

  class_call(background_at_tau(pba,
			       tau, 
			       normal_info, 
			       normal, 
			       &(ppw->last_index_back), 
			       ppw->pvecback),
	     pba->error_message,
	     ppt->error_message);
  
  a = ppw->pvecback[pba->index_bg_a];
  
  /** - for scalars */

  if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

    /** (a) adiabatic */ 

    if ((ppt->has_ad == _TRUE_) && (index_ic == ppt->index_ic_ad)) {

      /* relevant background quantities */

      /* 8piG/3 rho_r(t_i) */
      double rho_r = ppw->pvecback[pba->index_bg_rho_g];

      /* 8piG/3 rho_m(t_i) */
      double rho_m = ppw->pvecback[pba->index_bg_rho_b];

      /* 8piG/3 rho_nu(t_i) (all neutrinos and collisionless relics being relativistic at that time) */
      double rho_nu = 0.;

      if (pba->has_cdm == _TRUE_) {
	rho_m += ppw->pvecback[pba->index_bg_rho_cdm];
      }

      if (pba->has_ur == _TRUE_) {
	rho_r += ppw->pvecback[pba->index_bg_rho_ur];
	rho_nu += ppw->pvecback[pba->index_bg_rho_ur];
      }
      
      if (pba->has_ncdm == _TRUE_) {
	for(n_ncdm=0; n_ncdm<pba->N_ncdm; n_ncdm++){
	  rho_r += ppw->pvecback[pba->index_bg_rho_ncdm1 + n_ncdm];
	  rho_nu += ppw->pvecback[pba->index_bg_rho_ncdm1 + n_ncdm];
	}
      }
      
      class_test(rho_r == 0.,
		 ppt->error_message,
		 "stop to avoid division by zero");

      /* f_nu = Omega_nu(t_i) / Omega_r(t_i) */
      double fracnu = rho_nu/rho_r;

      /* f_nu = Omega_b(t_i) / Omega_m(t_i) */
      double fracb = ppw->pvecback[pba->index_bg_rho_b]/rho_m;

      /* alpha = 5 tau / (15 + 4 f_nu) */ 
      double alpha = 5.*tau/(15.+4.*fracnu);

      /* omega = Omega_m(t_i) a(t_i) H(t_i) / sqrt(Omega_r(t_i))
	 = Omega_m(t_0) a(t_0) H(t_0) / sqrt(Omega_r(t_0)) assuming rho_m in a-3 and rho_r in a^-4
	 = (8piG/3 rho_m(t_i)) a(t_i) / sqrt(8piG/3 rho_r(t_i))  in Mpc-1 */
      double om = a*rho_m/sqrt(rho_r);

      double ktau_two=k*k*tau*tau;
      double ktau_three=k*tau*ktau_two;

      /* newtonian gauge */
      if (ppr->gauge == newtonian) {

	class_test(tau == 0.,
		   ppt->error_message,
		   "stop to avoid division by zero");

	ppw->pv->y[ppw->pv->index_pt_delta_g] = -4.*alpha/tau - k*k*tau*tau/3.; /* photon density */
	ppw->pv->y[ppw->pv->index_pt_theta_g] = alpha*k*k - pow(k*tau,3.)*k/36.; /* photon velocity */

	ppw->pv->y[ppw->pv->index_pt_delta_b] = 3./4.*ppw->pv->y[ppw->pv->index_pt_delta_g]; /* baryon density */
	ppw->pv->y[ppw->pv->index_pt_theta_b] = ppw->pv->y[ppw->pv->index_pt_theta_g]; /* baryon velocity */
      
	if (pba->has_cdm == _TRUE_) {       
	  ppw->pv->y[ppw->pv->index_pt_delta_cdm] = 3./4.*ppw->pv->y[ppw->pv->index_pt_delta_g]; /* cdm density */
	  ppw->pv->y[ppw->pv->index_pt_theta_cdm] = alpha*k*k; /* cdm velocity */
	}
	
 	if (pba->has_fld == _TRUE_) {        
 	  ppw->pv->y[ppw->pv->index_pt_delta_fld] = 0.; /* fluid density (TO BE WRITTEN) */
 	  ppw->pv->y[ppw->pv->index_pt_theta_fld] = 0.; /* fluid velocity (TO BE WRITTEN) */
 	} 
	
	if (pba->has_ur == _TRUE_) {
	  ppw->pv->y[ppw->pv->index_pt_delta_ur] = ppw->pv->y[ppw->pv->index_pt_delta_g]; /* density of ultra-relativistic neutrinos/relics */
	  ppw->pv->y[ppw->pv->index_pt_theta_ur] = alpha*k*k - pow(k*tau,3.)*k/36. * (23.+4.*fracnu)/(15.+4.*fracnu); /* velocity of ultra-relativistic neutrinos/relics */
	  ppw->pv->y[ppw->pv->index_pt_shear_ur] = k*k*tau*tau*2./3./(12.+fracnu); /* shear of ultra-relativistic neutrinos/relics */
	}

      }

      /* synchronous gauge */
      if (ppr->gauge == synchronous) {


	/* valid at leading order in (k*tau), and order zero in tight-coupling */
	/* identical to first order terms in CRS, excpet for normalization (when ppr->curvature_ini=1, tau=1: leads to factor 1/2 difference between CRS formulas with beta1=0) */
        /* identical to CAMB when om set to zero in theta_g, theta_ur, shear_ur, tau */

	/* photon density */
	ppw->pv->y[ppw->pv->index_pt_delta_g] = - ktau_two/3. * (1.-om*tau/5.) 
	  * ppr->curvature_ini; 

	/* photon velocity */
	ppw->pv->y[ppw->pv->index_pt_theta_g] = - k*ktau_three/36. * (1.-3.*(1.+5.*fracb-fracnu)/20./(1.-fracnu)*om*tau) 
	  * ppr->curvature_ini;

	/* tighly-coupled baryons */
	ppw->pv->y[ppw->pv->index_pt_delta_b] = 3./4.*ppw->pv->y[ppw->pv->index_pt_delta_g]; /* baryon density */
	ppw->pv->y[ppw->pv->index_pt_theta_b] = ppw->pv->y[ppw->pv->index_pt_theta_g]; /* baryon velocity */
      
	if (pba->has_cdm == _TRUE_) {       
	  ppw->pv->y[ppw->pv->index_pt_delta_cdm] = 3./4.*ppw->pv->y[ppw->pv->index_pt_delta_g]; /* cdm density */
          /* by convention, cdm velocity velocity vanishes in this implementation of the synchronous gauge */
	}

	/* fluid density */
 	if (pba->has_fld == _TRUE_) {        
 	  ppw->pv->y[ppw->pv->index_pt_delta_fld] = - ktau_two/4.*(1.+pba->w_fld)*(4.-3.*pba->cs2_fld)/(4.-6.*pba->w_fld+3.*pba->cs2_fld) * ppr->curvature_ini; /* from 1004.5509 */
	  ppw->pv->y[ppw->pv->index_pt_theta_fld] = - k*ktau_three/4.*pba->cs2_fld/(4.-6.*pba->w_fld+3.*pba->cs2_fld) * ppr->curvature_ini; /* from 1004.5509 */
 	} 

	/* relativistic relics */
	if ((pba->has_ur == _TRUE_) || (pba->has_ncdm == _TRUE_)) {
	
	  delta_ur = ppw->pv->y[ppw->pv->index_pt_delta_g]; /* density of ultra-relativistic neutrinos/relics */

	  theta_ur = - k*ktau_three/36./(4.*fracnu+15.) * (4.*fracnu+23.-3.*(8.*fracnu*fracnu+50.*fracnu+275.)/20./(2.*fracnu+15.)*tau*om) 
	    * ppr->curvature_ini; /* velocity of ultra-relativistic neutrinos/relics */

	  shear_ur = 2.*ktau_two/(45.+12.*fracnu) * (1.+(4.*fracnu-5.)/4./(2.*fracnu+15.)*tau*om)
	    * ppr->curvature_ini; /* shear of ultra-relativistic neutrinos/relics */
	  
	  l3_ur = ktau_three*2./7./(12.*fracnu+45.)* ppr->curvature_ini;
	}

	if (pba->has_ur == _TRUE_) {

	  ppw->pv->y[ppw->pv->index_pt_delta_ur] = delta_ur;

	  ppw->pv->y[ppw->pv->index_pt_theta_ur] = theta_ur;

	  ppw->pv->y[ppw->pv->index_pt_shear_ur] = shear_ur;

	  ppw->pv->y[ppw->pv->index_pt_l3_ur] = l3_ur;

	}    

	if (pba->has_ncdm == _TRUE_) {
	  idx = ppw->pv->index_pt_psi0_ncdm1;
	  for (n_ncdm=0; n_ncdm < pba->N_ncdm; n_ncdm++){

	    for (index_q=0; index_q < ppw->pv->q_size_ncdm[n_ncdm]; index_q++) {
	    
	      q = pba->q_ncdm[n_ncdm][index_q];

	      epsilon = sqrt(q*q+a*a*pba->M_ncdm[n_ncdm]*pba->M_ncdm[n_ncdm]);

	      ppw->pv->y[idx] = -0.25 * delta_ur * pba->dlnf0_dlnq_ncdm[n_ncdm][index_q];

	      ppw->pv->y[idx+1] =  -epsilon/3./q/k*theta_ur* pba->dlnf0_dlnq_ncdm[n_ncdm][index_q];

	      ppw->pv->y[idx+2] = -0.5 * shear_ur * pba->dlnf0_dlnq_ncdm[n_ncdm][index_q];

	      ppw->pv->y[idx+3] = -0.25 * l3_ur * pba->dlnf0_dlnq_ncdm[n_ncdm][index_q];
		
	      //Jump to next momentum bin:
	      idx += (ppw->pv->l_max_ncdm[n_ncdm]+1);

	    }
	  }
	}

	/* metric perturbation eta */
	ppw->pv->y[ppw->pv->index_pt_eta] = ppr->curvature_ini *
	  (1.-ktau_two/12./(15.+4.*fracnu)*(5.+4.*fracnu - (16.*fracnu*fracnu+280.*fracnu+325)/10./(2.*fracnu+15.)*tau*om));

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

      class_test(pba->has_ur == _FALSE_,
		 ppt->error_message,
		 "not consistent to ask for NID in absence of neutrinos!");

      ppw->pv->y[ppw->pv->index_pt_delta_ur] = ppr->entropy_ini;

    }
     
    /** (e) Neutrino velocity Isocurvature */ 

    if ((ppt->has_niv == _TRUE_) && (index_ic == ppt->index_ic_niv)) {

      class_test(pba->has_ur == _FALSE_,
		 ppt->error_message,
		 "not consistent to ask for NIV in absence of neutrinos!");

      ppw->pv->y[ppw->pv->index_pt_theta_ur] = k*ppr->entropy_ini;

    }

  }

  /** - for tensors */

  if ((ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors)) {

    if (index_ic == ppt->index_ic_ten) {
      ppw->pv->y[ppw->pv->index_pt_gw] = ppr->gw_ini;
    }

  }

  return _SUCCESS_;
}

/**
 * Evaluate background/thermodynamics at \f$ \tau \f$, infer useful flags / time scales for integrating perturbations.
 *
 * Evaluate background quantities at \f$ \tau \f$, as well as thermodynamics for scalar mode; infer useful flags and time scales for integrating the perturbations:
 * - check whether tight-coupling approximation is needed.
 * - check whether radiation (photons, massless neutrinos...) perturbations are needed.
 * - choose step of integration: step = ppr->perturb_integration_stepsize * min_time_scale, where min_time_scale = smallest time scale involved in the equations. There are three time scales to compare:
 * -# that of recombination, \f$ \tau_c = 1/\kappa' \f$
 * -# Hubble time scale, \f$ \tau_h = a/a' \f$
 * -# Fourier mode, \f$ \tau_k = 1/k \f$
 *
 * So, in general, min_time_scale = \f$ \min(\tau_c, \tau_b, \tau_h, \tau_k) \f$.
 *
 * However, if \f$ \tau_c \ll \tau_h \f$ and \f$ \tau_c
 * \ll \tau_k \f$, we can use the tight-coupling regime for photons
 * and write equations in such way that the time scale \f$
 * \tau_c \f$ becomes irrelevant (no effective mass term in \f$
 * 1/\tau_c \f$).  Then, the smallest
 * scale in the equations is only \f$ \min(\tau_h, \tau_k) \f$.
 * In practise, it is sufficient to use only the condition \f$ \tau_c \ll \tau_h \f$.
 * 
 * Also, if \f$ \rho_{matter} \gg \rho_{radiation} \f$ and \f$ k \gg
 * aH \f$, we can switch off radiation perturbations (i.e. switch on
 * the free-streaming approximation) and then the smallest scale is
 * simply \f$ \tau_h \f$.
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to thermodynamics structure
 * @param ppt        Input: pointer to the perturbation structure
 * @param index_mode Input: index of mode under consideration (scalar/.../tensor)
 * @param k          Input: wavenumber
 * @param tau        Input: conformal time
 * @param ppw        Input/Output: in output contains the approximation to be used at this time
 * @return the error status
 */

int perturb_approximations(
			   struct precision * ppr,
			   struct background * pba,
			   struct thermo * pth,
			   struct perturbs * ppt,
			   int index_mode,
			   double k,
			   double tau,
			   struct perturb_workspace * ppw
			   ) {
  /** Summary: */

  /** - define local variables */

  /* (a) time scale of Fourier mode, \f$ \tau_k = 1/k \f$ */  
  double tau_k;
  /* (b) time scale of expansion, \f$ \tau_h = a/a' \f$ */
  double tau_h;
  /* (c) time scale of recombination, \f$ \tau_{\gamma} = 1/\kappa' \f$ */
  double tau_c;

  /** - compute Fourier mode time scale = \f$ \tau_k = 1/k \f$ */

  class_test(k == 0.,
	     ppt->error_message,
	     "stop to avoid division by zero");

  tau_k = 1./k;

  /** - evaluate background quantities with background_at_tau() and
      Hubble time scale \f$ \tau_h = a/a' \f$ */

  class_call(background_at_tau(pba,tau, normal_info, ppw->intermode, &(ppw->last_index_back), ppw->pvecback),
	     pba->error_message,
	     ppt->error_message);

  class_test(ppw->pvecback[pba->index_bg_H]*ppw->pvecback[pba->index_bg_a] == 0.,
	     ppt->error_message,
	     "aH=0, stop to avoid division by zero");

  tau_h = 1./(ppw->pvecback[pba->index_bg_H]*ppw->pvecback[pba->index_bg_a]);

  /** - for scalars modes: */

  if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {
    
    /** (a) evaluate thermodynamical quantities with thermodynamics_at_z() */

    class_call(thermodynamics_at_z(pba,
				   pth,
				   1./ppw->pvecback[pba->index_bg_a]-1.,  /* redshift z=1/a-1 */
				   ppw->intermode,
				   &(ppw->last_index_thermo),
				   ppw->pvecback,
				   ppw->pvecthermo),
	       pth->error_message,
	       ppt->error_message);
    
    /** (b.1) if \f$ \kappa'=0 \f$, recombination is finished; tight-coupling approximation must be off */

    if (ppw->pvecthermo[pth->index_th_dkappa] == 0.) {

      ppw->approx[ppw->index_ap_tca] = (int)tca_off;
      
    }

    /** (b.2) if \f$ \kappa' \neq 0 \f$, recombination is not finished: check tight-coupling approximation */

    else {

      /** (b.2.a) compute recombination time scale for photons, \f$ \tau_{\gamma} = 1/ \kappa' \f$ */
      tau_c = 1./ppw->pvecthermo[pth->index_th_dkappa];

      /** (b.2.b) check whether tight-coupling approximation should be on */
      if ((tau_c/tau_h < ppr->tight_coupling_trigger_tau_c_over_tau_h) &&
	  (tau_c/tau_k < ppr->tight_coupling_trigger_tau_c_over_tau_k)) {
	ppw->approx[ppw->index_ap_tca] = (int)tca_on;
      }
      else {
	ppw->approx[ppw->index_ap_tca] = (int)tca_off;
      }

    }

    /* (c) free-streaming approximations */

    if ((tau_h/tau_k > ppr->radiation_streaming_trigger_tau_h_over_tau_k) &&
	(tau > pth->tau_free_streaming) &&
	(ppr->radiation_streaming_approximation != rsa_none)) {

      ppw->approx[ppw->index_ap_rsa] = (int)rsa_on;
    }
    else {
      ppw->approx[ppw->index_ap_rsa] = (int)rsa_off;
    }
   
    if (pba->has_ur == _TRUE_) {

      if ((tau_h/tau_k > ppr->ur_fluid_trigger_tau_h_over_tau_k) &&
	  (ppr->ur_fluid_approximation != ufa_none)) {
	
	ppw->approx[ppw->index_ap_ufa] = (int)ufa_on;
      }
      else {
	ppw->approx[ppw->index_ap_ufa] = (int)ufa_off;
      }  
    }

    if (pba->has_ncdm == _TRUE_) {

      if ((tau_h/tau_k > ppr->ncdm_fluid_trigger_tau_h_over_tau_k) &&
	  (ppr->ncdm_fluid_approximation != ncdmfa_none)) {
	
	ppw->approx[ppw->index_ap_ncdmfa] = (int)ncdmfa_on;
      }
      else {
	ppw->approx[ppw->index_ap_ncdmfa] = (int)ncdmfa_off;
      }  
    }
  }

  /* no approximation implemented so far for tensors */

  return _SUCCESS_;
}

/**
 * Compute typical timescale over which the perturbation equation
 * vary. Some integrators (e.g. Runge-Kunta) benefit from calling this
 * routine at each step in order to adapt the next step.
 *
 * This is one of the few functions in the code which are passed to the generic_integrator() routine. 
 * Since generic_integrator() should work with functions passed from various modules, the format of the arguments
 * is a bit special:
 * - fixed parameters and workspaces are passed through a generic pointer. 
 *   generic_integrator() doesn't know the content of this pointer.
 * - the error management is a bit special: errors are not written as usual to pth->error_message, but to a generic 
 *   error_message passed in the list of arguments.
 *
 * @param tau                      Input : conformal time
 * @param parameters_and_workspace Input : fixed parameters (e.g. indices), workspace, approximation used, etc.
 * @param timescale                Output: perturbation variation timescale (given the apprtoximation used)  
 * @param error_message            Output: error message
 */

int perturb_timescale(
		      double tau,
		      void * parameters_and_workspace,
		      double * timescale,
		      ErrorMsg error_message
		      ) {
  /** Summary: */

  /** - define local variables */

  /* (a) time scale of Fourier mode, \f$ \tau_k = 1/k \f$ */  
  double tau_k;
  /* (b) time scale of expansion, \f$ \tau_h = a/a' \f$ */
  double tau_h;
  /* (c) time scale of recombination, \f$ \tau_{\gamma} = 1/\kappa' \f$ */
  double tau_c;

  /* various pointers allowing to extract the fields of the
     parameter_and_workspace input structure */
  struct perturb_parameters_and_workspace * pppaw;
  struct background * pba;
  struct thermo * pth;
  struct perturbs * ppt;
  struct perturb_workspace * ppw;
  double * pvecback;
  double * pvecthermo;

  /** - extract the fields of the parameter_and_workspace input structure */
  pppaw = parameters_and_workspace;
  pba = pppaw->pba;
  pth = pppaw->pth;
  ppt = pppaw->ppt;
  ppw = pppaw->ppw;
  pvecback = ppw->pvecback;
  pvecthermo = ppw->pvecthermo;

  /** - compute Fourier mode time scale = \f$ \tau_k = 1/k \f$ */

  class_test(pppaw->k == 0.,
	     ppt->error_message,
	     "stop to avoid division by zero");

  tau_k = 1./pppaw->k;

  /** - evaluate background quantities with background_at_tau() and
      Hubble time scale \f$ \tau_h = a/a' \f$ */

  class_call(background_at_tau(pba,tau, normal_info, ppw->intermode, &(ppw->last_index_back), pvecback),
	     pba->error_message,
	     error_message);

  class_test(pvecback[pba->index_bg_H]*pvecback[pba->index_bg_a] == 0.,
	     error_message,
	     "aH=0, stop to avoid division by zero");

  tau_h = 1./(pvecback[pba->index_bg_H]*pvecback[pba->index_bg_a]);

  /** - for scalars modes: */

  if ((ppt->has_scalars == _TRUE_) && (pppaw->index_mode == ppt->index_md_scalars)) {

    *timescale = tau_h;

    if ((ppw->approx[ppw->index_ap_rsa] == (int)rsa_off) || (pba->has_ncdm == _TRUE_))
      *timescale = min(tau_k,*timescale);

    if (ppw->approx[ppw->index_ap_tca] == (int)tca_off) {

      class_call(thermodynamics_at_z(pba,
				     pth,
				     1./pvecback[pba->index_bg_a]-1.,  /* redshift z=1/a-1 */
				     ppw->intermode,
				     &(ppw->last_index_thermo),
				     pvecback,
				     pvecthermo),
		 pth->error_message,
		 error_message);

      if (pvecthermo[pth->index_th_dkappa] != 0.) {

	/** (b.2.a) compute recombination time scale for photons, \f$ \tau_{\gamma} = 1/ \kappa' \f$ */

	tau_c = 1./pvecthermo[pth->index_th_dkappa];

	*timescale = min(tau_c,*timescale);

      }
    }

  }

  /** - for tensor modes: tight-coupling approximation is off, time scale remains \f$ min (\tau_k , \tau_h) \f$. */

  if ((ppt->has_tensors == _TRUE_) && (pppaw->index_mode == ppt->index_md_tensors)) {

    *timescale = min(tau_k,tau_h);

  }
  
  /** - vectors not coded yet */
  
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
 * @param tau        Input: conformal time
 * @param y          Input: vector of perturbations (those integrated over time) (already allocated)
 * @param ppw        Input/Output: in output contains the updated metric perturbations
 * @return the error status
 */

int perturb_einstein(
		     struct precision * ppr,
		     struct background * pba,
		     struct thermo * pth,
		     struct perturbs * ppt,
		     int index_mode,
		     double k,
		     double tau,
		     double * y,
		     struct perturb_workspace * ppw
		     ) {
  /** Summary: */

  /** - define local variables */

  double k2,a,a2,a_prime_over_a;
  double delta_rho,rho_plus_p_theta,rho_plus_p_shear,delta_p;
  double delta_g=0.;
  double theta_g=0.; 
  double shear_g=0.;
  double delta_ur=0.;
  double theta_ur=0.; 
  double shear_ur=0.; 
  double rho_delta_ncdm=0.;
  double rho_plus_p_theta_ncdm=0.; 
  double rho_plus_p_shear_ncdm=0.;
  double delta_p_ncdm=0.;
  double alpha;
  double factor;
  double rho_plus_p_ncdm;
  int index_q,n_ncdm,idx;
  double epsilon,q,q2,cg2_ncdm,w_ncdm,rho_ncdm_bg,p_ncdm_bg,pseudo_p_ncdm;
  double rho_pk,delta_rho_pk;

  /** - wavenumber and scale factor related quantities */ 

  k2 = k*k;
  a = ppw->pvecback[pba->index_bg_a];
  a2 = a * a;
  a_prime_over_a = ppw->pvecback[pba->index_bg_H]*a;

  /** - for scalar modes: */  

  if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

    /** (a) deal with approximation schemes */
    
    /** (a.1) photons */

    if (ppw->approx[ppw->index_ap_tca] == (int)tca_off) {

      if (ppw->approx[ppw->index_ap_rsa] == (int)rsa_off) {

	/** (a.1.1) no approximation */

	delta_g = y[ppw->pv->index_pt_delta_g];
	theta_g = y[ppw->pv->index_pt_theta_g];
	shear_g = y[ppw->pv->index_pt_shear_g];

      }
      else {

	/** (a.1.2) radiation streaming approximation */

	delta_g = 0.; /* actual free streaming approximation imposed after evaluation of 1st einstein equation */
	theta_g = 0.; /* actual free streaming approximation imposed after evaluation of 1st einstein equation */
	shear_g = 0.; /* shear always neglected in free streaming approximatio */
      }
    }
    else {

      /** (a.1.3) tight coupling approximation */

      delta_g = y[ppw->pv->index_pt_delta_g];
      theta_g = y[ppw->pv->index_pt_theta_g];
      
      shear_g = 0.; /* in the tight-coupling approximation, the
		       expression of shear_g (at first-order in a
		       tight-coupling expansion) is a function of h'
		       and eta'; but h' and eta' are calculated below
		       as a function of delta_g and theta_g.  Hence,
		       we set shear_g temporarily to zero, and set it
		       to the right first-order value below, just
		       before using the Einstein equation for the
		       shear. */
    }

    /** (a.2) ur */

    if (pba->has_ur == _TRUE_) {

      if (ppw->approx[ppw->index_ap_rsa] == (int)rsa_off) {

	delta_ur = y[ppw->pv->index_pt_delta_ur];
	theta_ur = y[ppw->pv->index_pt_theta_ur];
	shear_ur = y[ppw->pv->index_pt_shear_ur];

      }

      else {

	delta_ur = 0.; /* actual free streaming approximation imposed after evaluation of 1st einstein equation */
	theta_ur = 0.; /* actual free streaming approximation imposed after evaluation of 1st einstein equation */
	shear_ur = 0.; /* shear always neglected in free streaming approximatio */

      }

    }

    /** (b) compute the total density, velocity and shear perturbations */
 
    /* photon and baryon contribution */
    delta_rho = ppw->pvecback[pba->index_bg_rho_g]*delta_g
      + ppw->pvecback[pba->index_bg_rho_b]*y[ppw->pv->index_pt_delta_b];
    rho_plus_p_theta = 4./3.*ppw->pvecback[pba->index_bg_rho_g]*theta_g
      + ppw->pvecback[pba->index_bg_rho_b]*y[ppw->pv->index_pt_theta_b];
    rho_plus_p_shear = 4./3.*ppw->pvecback[pba->index_bg_rho_g]*shear_g;
    delta_p = 1./3.*ppw->pvecback[pba->index_bg_rho_g]*delta_g
      + ppw->pvecthermo[pth->index_th_cb2]*ppw->pvecback[pba->index_bg_rho_b]*y[ppw->pv->index_pt_delta_b];

    /* cdm contribution */
    if (pba->has_cdm == _TRUE_) {
      delta_rho = delta_rho + ppw->pvecback[pba->index_bg_rho_cdm]*y[ppw->pv->index_pt_delta_cdm];      
      if (ppr->gauge == newtonian)
	rho_plus_p_theta = rho_plus_p_theta + ppw->pvecback[pba->index_bg_rho_cdm]*y[ppw->pv->index_pt_theta_cdm];
    }
    
    /* fluid contribution */
    if (pba->has_fld == _TRUE_) { 
      delta_rho = delta_rho + ppw->pvecback[pba->index_bg_rho_fld]*y[ppw->pv->index_pt_delta_fld]; 
      rho_plus_p_theta = rho_plus_p_theta + ppw->pvecback[pba->index_bg_rho_fld]*y[ppw->pv->index_pt_theta_fld];
      delta_p = delta_p + pba->cs2_fld * ppw->pvecback[pba->index_bg_rho_fld]*y[ppw->pv->index_pt_delta_fld]; 
    } 

    /* ultra-relativistic neutrino/relics contribution */
    if (pba->has_ur == _TRUE_) {
      delta_rho = delta_rho + ppw->pvecback[pba->index_bg_rho_ur]*delta_ur;
      rho_plus_p_theta = rho_plus_p_theta + 4./3.*ppw->pvecback[pba->index_bg_rho_ur]*theta_ur;
      rho_plus_p_shear = rho_plus_p_shear + 4./3.*ppw->pvecback[pba->index_bg_rho_ur]*shear_ur;
      delta_p += 1./3.*ppw->pvecback[pba->index_bg_rho_ur]*delta_ur;
    }

    /* non-cold dark matter contribution */
    if (pba->has_ncdm == _TRUE_) {
      idx = ppw->pv->index_pt_psi0_ncdm1;
      if(ppw->approx[ppw->index_ap_ncdmfa] == (int)ncdmfa_on){
	// The perturbations are evolved integrated:
	for(n_ncdm=0; n_ncdm < pba->N_ncdm; n_ncdm++){
	  rho_ncdm_bg = ppw->pvecback[pba->index_bg_rho_ncdm1+n_ncdm];
	  p_ncdm_bg = ppw->pvecback[pba->index_bg_p_ncdm1+n_ncdm];
	  pseudo_p_ncdm = ppw->pvecback[pba->index_bg_pseudo_p_ncdm1+n_ncdm];

	  rho_plus_p_ncdm = rho_ncdm_bg + p_ncdm_bg;
	  w_ncdm = p_ncdm_bg/rho_ncdm_bg;
	  cg2_ncdm = w_ncdm*(1.0-1.0/(3.0+3.0*w_ncdm)*(3.0*w_ncdm-2.0+pseudo_p_ncdm/p_ncdm_bg));
	  if ((ppt->has_source_delta_ncdm == _TRUE_) || (ppt->has_source_delta_pk == _TRUE_)) {
	    ppw->delta_ncdm[n_ncdm] = y[idx];
	    ppw->theta_ncdm1 = y[ppw->pv->index_pt_psi0_ncdm1+1];
	    ppw->shear_ncdm1 = y[ppw->pv->index_pt_psi0_ncdm1+2];
	  }

	  delta_rho += rho_ncdm_bg*y[idx];
	  rho_plus_p_theta += rho_plus_p_ncdm*y[idx+1];
	  rho_plus_p_shear += rho_plus_p_ncdm*y[idx+2]; 
	  delta_p += cg2_ncdm*rho_ncdm_bg*y[idx];
	  idx += ppw->pv->l_max_ncdm[n_ncdm]+1;
        }
      }
      else{
	// We must integrate to find perturbations:
	for(n_ncdm=0; n_ncdm < pba->N_ncdm; n_ncdm++){
	  rho_delta_ncdm = 0.0;
	  rho_plus_p_theta_ncdm = 0.0;
	  rho_plus_p_shear_ncdm = 0.0;
	  delta_p_ncdm = 0.0;
	  factor = pba->factor_ncdm[n_ncdm]*pow(pba->a_today/a,4);
	 
	  for (index_q=0; index_q < ppw->pv->q_size_ncdm[n_ncdm]; index_q ++) {
	    
	    q = pba->q_ncdm[n_ncdm][index_q];
	    q2 = q*q;
	    epsilon = sqrt(q2+pba->M_ncdm[n_ncdm]*pba->M_ncdm[n_ncdm]*a2);

	    rho_delta_ncdm += q2*epsilon*pba->w_ncdm[n_ncdm][index_q]*y[idx];
	    rho_plus_p_theta_ncdm += q2*q*pba->w_ncdm[n_ncdm][index_q]*y[idx+1];
	    rho_plus_p_shear_ncdm += q2*q2/epsilon*pba->w_ncdm[n_ncdm][index_q]*y[idx+2];
	    delta_p_ncdm += q2*q2/epsilon*pba->w_ncdm[n_ncdm][index_q]*y[idx];

	    //Jump to next momentum bin:
	    idx+=(ppw->pv->l_max_ncdm[n_ncdm]+1);
	  }

	  rho_delta_ncdm *= factor;
	  rho_plus_p_theta_ncdm *= k*factor;
	  rho_plus_p_shear_ncdm *= 2.0/3.0*factor;
	  delta_p_ncdm *= factor/3.;

	  if ((ppt->has_source_delta_ncdm == _TRUE_) || (ppt->has_source_delta_pk == _TRUE_)) {
	    ppw->delta_ncdm[n_ncdm] = rho_delta_ncdm/ppw->pvecback[pba->index_bg_rho_ncdm1+n_ncdm];
	    ppw->theta_ncdm1 = rho_plus_p_theta_ncdm/
	      (ppw->pvecback[pba->index_bg_rho_ncdm1]+ppw->pvecback[pba->index_bg_p_ncdm1]);
	    ppw->shear_ncdm1 = rho_plus_p_shear_ncdm/
	      (ppw->pvecback[pba->index_bg_rho_ncdm1]+ppw->pvecback[pba->index_bg_p_ncdm1]);
	  }

	  delta_rho += rho_delta_ncdm;
	  rho_plus_p_theta += rho_plus_p_theta_ncdm;
	  rho_plus_p_shear += rho_plus_p_shear_ncdm;
	  delta_p += delta_p_ncdm;
	}
      }
    }
    

    /* store delta_pk (for corresponding source function). Since the
       matter power spectrum is usually defined in such way to include
       only non-relativistic components, the sum over each species
       contribution to delta_rho_pk and rho_pk must be done
       'manually'. Only if the P(k) is defined from the total matter
       overdensity we can use the delta_rho computed above.  */

    if (ppt->has_source_delta_pk == _TRUE_) {

      /* do the sum over species contributing to delta_pk */

      if ((ppr->pk_definition == delta_m_squared) ||
	  (ppr->pk_definition == delta_bc_squared)){

	/* include baryons and cold dark matter */

	delta_rho_pk = ppw->pvecback[pba->index_bg_rho_b]*y[ppw->pv->index_pt_delta_b];
	rho_pk = ppw->pvecback[pba->index_bg_rho_b];

	if (pba->has_cdm == _TRUE_) {
	  delta_rho_pk += ppw->pvecback[pba->index_bg_rho_cdm]*y[ppw->pv->index_pt_delta_cdm];      
	  rho_pk += ppw->pvecback[pba->index_bg_rho_cdm];
	}
	 
	/* include any other species non-relativistic today (like ncdm
	   species) */

	if (ppr->pk_definition == delta_m_squared) {
 
	  if (pba->has_ncdm == _TRUE_) {

	    for(n_ncdm=0; n_ncdm < pba->N_ncdm; n_ncdm++){

	      delta_rho_pk += ppw->pvecback[pba->index_bg_rho_ncdm1+n_ncdm]*ppw->delta_ncdm[n_ncdm];
	      rho_pk += ppw->pvecback[pba->index_bg_rho_ncdm1+n_ncdm];
	    }
	  }
	}

	/* infer delta_pk */

	ppw->delta_pk = delta_rho_pk/rho_pk;
      }
      
      /* alternative: take directly into account all species in
	 delta_rho, and matter species in rho (final result differs
	 from delta_m_squared only for modes close to Hubble scale) */

      if (ppr->pk_definition == delta_tot_squared) {

	rho_pk = ppw->pvecback[pba->index_bg_rho_b];
	if (pba->has_cdm == _TRUE_)
	  rho_pk += ppw->pvecback[pba->index_bg_rho_cdm];
	if (pba->has_ncdm == _TRUE_) {
	  for(n_ncdm=0; n_ncdm < pba->N_ncdm; n_ncdm++){
	    rho_pk += ppw->pvecback[pba->index_bg_rho_ncdm1+n_ncdm];
	  }
	}
	ppw->delta_pk = delta_rho/rho_pk;
      }
      
    }


    /** (c) infer metric perturbations from Einstein equations */

    /* newtonian gauge */
    if (ppr->gauge == newtonian) {
      ppw->pvecmetric[ppw->index_mt_phi] = -1.5 * (a2/k2/k2) * (k2 * delta_rho + 3.*a_prime_over_a * rho_plus_p_theta); /* phi */
      ppw->pvecmetric[ppw->index_mt_psi] = ppw->pvecmetric[ppw->index_mt_phi] - 4.5 * (a2/k2) * rho_plus_p_shear;  /* psi */
      ppw->pvecmetric[ppw->index_mt_phi_prime] = - a_prime_over_a * ppw->pvecmetric[ppw->index_mt_psi] + 1.5 * (a2/k2) * rho_plus_p_theta; /* phi' */
    }

    /* synchronous gauge */
    if (ppr->gauge == synchronous) {

      /* first equation involving total density fluctuation */
      ppw->pvecmetric[ppw->index_mt_h_prime] = 
	( k2 * y[ppw->pv->index_pt_eta] + 1.5 * a2 * delta_rho)/(0.5*a_prime_over_a);  /* h' */

      /* eventually, infer streaming approximation for gamma and
	 correct the total velocity */

      if (ppw->approx[ppw->index_ap_rsa] == (int)rsa_on) {

	if (ppr->radiation_streaming_approximation == rsa_null) {
	  ppw->rsa_delta_g = 0.;
	  ppw->rsa_theta_g = 0.;
	}
	else {
	  
	  ppw->rsa_delta_g = 4./k2*(a_prime_over_a*ppw->pvecmetric[ppw->index_mt_h_prime]
				    -k2*y[ppw->pv->index_pt_eta]);
	  ppw->rsa_theta_g = -0.5*ppw->pvecmetric[ppw->index_mt_h_prime];
	}
	    
	if (ppr->radiation_streaming_approximation == rsa_MD_with_reio) {
	  
	  ppw->rsa_delta_g += 
	    -4./k2*ppw->pvecthermo[pth->index_th_dkappa]*(y[ppw->pv->index_pt_theta_b]+0.5*ppw->pvecmetric[ppw->index_mt_h_prime]);
	  
	  ppw->rsa_theta_g += 
	    3./k2*(ppw->pvecthermo[pth->index_th_ddkappa]*
		   (y[ppw->pv->index_pt_theta_b]
		    +0.5*ppw->pvecmetric[ppw->index_mt_h_prime])
		   +ppw->pvecthermo[pth->index_th_dkappa]*
		   (-a_prime_over_a*y[ppw->pv->index_pt_theta_b]
		    + ppw->pvecthermo[pth->index_th_cb2]*k2*y[ppw->pv->index_pt_delta_b]
		    -a_prime_over_a*ppw->pvecmetric[ppw->index_mt_h_prime]
		    +k2*y[ppw->pv->index_pt_eta]));
	}	    

	rho_plus_p_theta += 4./3.*ppw->pvecback[pba->index_bg_rho_g]*ppw->rsa_theta_g;

	if (pba->has_ur == _TRUE_) {

	  if (ppr->radiation_streaming_approximation == rsa_null) {
	    ppw->rsa_delta_ur = 0.;
	    ppw->rsa_theta_ur = 0.;
	  }
	  else {
	    ppw->rsa_delta_ur = 4./k2*(a_prime_over_a*ppw->pvecmetric[ppw->index_mt_h_prime]
					-k2*y[ppw->pv->index_pt_eta]);
	    ppw->rsa_theta_ur = -0.5*ppw->pvecmetric[ppw->index_mt_h_prime];
	  }

	  rho_plus_p_theta += 4./3.*ppw->pvecback[pba->index_bg_rho_ur]*ppw->rsa_theta_ur;

	}
      }
      
      /* second equation involving total velocity */
      ppw->pvecmetric[ppw->index_mt_eta_prime] = 1.5 * (a2/k2) * rho_plus_p_theta;  /* eta' */

      /* third equation involving total pressure */
      ppw->pvecmetric[ppw->index_mt_h_prime_prime] = 
	-2.*a_prime_over_a*ppw->pvecmetric[ppw->index_mt_h_prime]
	+2.*k2*y[ppw->pv->index_pt_eta]
	-9.*a2*delta_p;

      /* intermediate quantity: alpha = (h'+6eta')/2k^2 */
      alpha = (ppw->pvecmetric[ppw->index_mt_h_prime] + 6.*ppw->pvecmetric[ppw->index_mt_eta_prime])/2./k2;

      /* eventually, infer first-order tight-coupling approximation for photon
	 shear, then correct the total shear */
      if (ppw->approx[ppw->index_ap_tca] == (int)tca_on) {
	
	shear_g = 16./45./ppw->pvecthermo[pth->index_th_dkappa]*(y[ppw->pv->index_pt_theta_g]+k2*alpha);
		
	rho_plus_p_shear += 4./3.*ppw->pvecback[pba->index_bg_rho_g]*shear_g;
	
      }
      
      /* fourth equation involving total shear */
      ppw->pvecmetric[ppw->index_mt_alpha_prime] = 
	- 4.5 * (a2/k2) * rho_plus_p_shear + y[ppw->pv->index_pt_eta] - 2.*a_prime_over_a*alpha; /* alpha' = (h''+6eta'')/2k2 */

      /* getting phi here is an option */
      /* phi=y[ppw->pv->index_pt_eta]-0.5 * (a_prime_over_a/k2) * (h_plus_six_eta_prime); */   
      /* phi from gauge transformation (from synchronous to newtonian) */

    }
  }

  /* nothing to be done for tensors: only one propagating degree of
     freedom, no constraint equation */

  return _SUCCESS_;

}

/**
 * Compute the terms contributing to the source functions.
 *
 * pvecsourceterms. The source functions can be decomposed as \f$ S =
 * S_0 + S_1' + S_2'' \f$.  This function computes \f$ ( S_0, S_1',
 * S_2') \f$ for temperature and \f$ ( S_0, S_1', S_2'') \f$ for other
 * quantitites.
 *
 * This is one of the few functions in the code which are passed to
 * the generic_integrator() routine. Since generic_integrator()
 * should work with functions passed from various modules, the format
 * of the arguments is a bit special: 
 *
 * - fixed parameters and workspaces are passed through a generic
 * pointer.  generic_integrator() doesn't know the content of this
 * pointer.
 *
 * - the error management is a bit special: errors are not written as
 * usual to pth->error_message, but to a generic error_message passed
 * in the list of arguments.
 *
 * @param tau                      Input: conformal time 
 * @param y                        Input: vector of perturbations
 * @param dy                       Input: vector of time derivative of perturbations
 * @param index_tau                Input: index in the array tau_sampling, specifies where the result should be stored in the source_term_table
 * @param parameters_and_workspace Input/Output: in input, all parameters needed by perturb_derivs, in ourput, source terms
 * @param error_message            Output: error message
 * @return the error status
 */

int perturb_source_terms(
			 double tau,
			 double * y,
			 double * dy,
			 int index_tau,
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

  /* class_test(ppw->approx[ppw->index_ap_tca] == (int)tca_on, */
  /* 	     ppt->error_message, */
  /* 	     "source calculation assume tight-coupling approximation turned off"); */

  x = k * (pba->conformal_age-tau);

  /** - get background/thermo quantities in this point */

  class_call(background_at_tau(pba,
			       tau, 
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

    /** - compute useful background quantities \f$ */ 
      
    k2 = k * k;
    
    a_prime_over_a = pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a];

    a_primeprime_over_a = pvecback[pba->index_bg_H_prime] * pvecback[pba->index_bg_a] 
      + 2. * a_prime_over_a * a_prime_over_a;

    R = 4./3. * pvecback[pba->index_bg_rho_g]/pvecback[pba->index_bg_rho_b];

    /** - compute metric perturbations */

    class_call(perturb_einstein(ppr,
				pba,
				pth,
				ppt,
				index_mode,
				k,
				tau,
				y,
				ppw),
	       ppt->error_message,
	       error_message);

    /** - compute quantities depending on approximation schemes */

    if (ppw->approx[ppw->index_ap_rsa] == (int)rsa_on) {
      delta_g = ppw->rsa_delta_g;
      Pi = 0.;
      Pi_prime = 0.;
    }
    else {
      delta_g = y[ppw->pv->index_pt_delta_g];

      if (ppw->approx[ppw->index_ap_tca] == (int)tca_on) {
	Pi = 5.*ppw->tca_shear_g; /* (2.5+0.5+2)shear_g */
	Pi_prime = 5.*ppw->tca_shear_g_prime; /* (2.5+0.5+2)shear_g_prime */
      }
      else {
	Pi = y[ppw->pv->index_pt_pol0_g] + y[ppw->pv->index_pt_pol2_g] + 2.*y[ppw->pv->index_pt_shear_g];
	Pi_prime = dy[ppw->pv->index_pt_pol0_g] + dy[ppw->pv->index_pt_pol2_g] + 2.*dy[ppw->pv->index_pt_shear_g];
      }
    }

    /** - for each type and each mode, compute S0, S1, S2 */
    for (index_type = 0; index_type < ppt->tp_size[index_mode]; index_type++) {

      source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_tau] = tau;

      /* scalar temperature */
      if ((ppt->has_source_t == _TRUE_) && (index_type == ppt->index_tp_t)) {

	/* check that visibility is non-zero (otherwise source = 0) */
	if (pvecthermo[pth->index_th_g] != 0.) {

	  /* newtonian gauge */
	  if (ppr->gauge == newtonian) {

	    /* S0 */
	    source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S0] =
	      pvecthermo[pth->index_th_exp_m_kappa] * pvecmetric[ppw->index_mt_phi_prime]
	      + pvecthermo[pth->index_th_g] / 4. * (delta_g + Pi / 4.);
	    
	    /* S1 */
	    source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S1] =
	      pvecthermo[pth->index_th_exp_m_kappa] * pvecmetric[ppw->index_mt_psi]
	      + pvecthermo[pth->index_th_g] * y[ppw->pv->index_pt_theta_b] / k2;

	    /* S2 */
	    source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S2] =
	      3./16. * pvecthermo[pth->index_th_g] * Pi / k2;

	  }

	  /* synchronous gauge */
	  if (ppr->gauge == synchronous) {

	    /* S0 */
	    source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S0] =
	      pvecthermo[pth->index_th_exp_m_kappa] * pvecmetric[ppw->index_mt_eta_prime]
	      + pvecthermo[pth->index_th_g] / 4. * (delta_g + Pi / 4.);
	  
	    /* S1 */
	    /* source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S1] = */
	    /* 	      pvecthermo[pth->index_th_g] * y[ppw->pv->index_pt_theta_b] / k2; */

	    /* dS1 */
	    source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_dS1] =
	      pvecthermo[pth->index_th_dg] * y[ppw->pv->index_pt_theta_b] / k2
	      + pvecthermo[pth->index_th_g] * dy[ppw->pv->index_pt_theta_b] / k2;

	    /* S2 */
	    /* source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S2] = */
	    /* 	      pvecthermo[pth->index_th_exp_m_kappa] * (pvecmetric[ppw->index_mt_h_prime] + 6. * pvecmetric[ppw->index_mt_eta_prime])/2./k2  */
	    /* 	      + 3./16. * pvecthermo[pth->index_th_g] * Pi / k2; */

	    /* dS2 */
	    source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_dS2] =
	      pvecthermo[pth->index_th_g] 
	      * (pvecmetric[ppw->index_mt_h_prime] + 6. * pvecmetric[ppw->index_mt_eta_prime])/2./k2
	      + pvecthermo[pth->index_th_exp_m_kappa] * (pvecmetric[ppw->index_mt_alpha_prime])
	      + 3./16. * pvecthermo[pth->index_th_dg] * Pi / k2
	      + 3./16. * pvecthermo[pth->index_th_g] * Pi_prime / k2;

	    /* 	    Pi_prime = -k*y[ppw->pv->index_pt_pol0_g+1] */
	    /* 	      -pvecthermo[pth->index_th_dkappa] * (y[ppw->pv->index_pt_pol0_g]-Pi/2.) */
	    /* 	      +k/5. * (2.*y[ppw->pv->index_pt_pol2_g-1]-3.*y[ppw->pv->index_pt_pol2_g+1]) */
	    /* 	      -pvecthermo[pth->index_th_dkappa] * (y[ppw->pv->index_pt_pol2_g]-Pi/10.) */
	    /*               +8./15.*y[ppw->pv->index_pt_theta_g] */
	    /* 	      -3./5.*k*y[ppw->pv->index_pt_shear_g+1] */
	    /* 	      +4./15.*pvecmetric[ppw->index_mt_h_prime] */
	    /* 	      +8./5.*pvecmetric[ppw->index_mt_eta_prime] */
	    /* 	      -pvecthermo[pth->index_th_dkappa] * (2.*y[ppw->pv->index_pt_shear_g]-1./10.*Pi); */

	    /* 	    Pi_prime = -3./10.*pvecthermo[pth->index_th_dkappa]*Pi */
	    /* 	      -3./5.*k * (y[ppw->pv->index_pt_pol1_g]+ */
	    /* 				y[ppw->pv->index_pt_pol2_g+1]+ */
	    /* 				y[ppw->pv->index_pt_shear_g+1]) */
	    /* 	      +8./15.*y[ppw->pv->index_pt_theta_g] */
	    /* 	      +4./15. * (pvecmetric[ppw->index_mt_h_prime] + 6. * pvecmetric[ppw->index_mt_eta_prime]); */

	  }
	}

	else {

	  source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S0] = 0.;
	  source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_dS1] = 0.;
	  source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_dS2] = 0.;

	}
      }

      /* scalar polarization */
      if ((ppt->has_source_e == _TRUE_) && (index_type == ppt->index_tp_e)) {

	if (x > 0.) {
	  /* in all gauges */
	  source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S0] =
	    + 3./16. * pvecthermo[pth->index_th_g] * Pi /x/x;  
	}
	else {
	  source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S0] = 0.;
	}
	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_dS1] = 0.;
	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_ddS2] = 0.;

      }

      /* gravitational potential */
      if ((ppt->has_source_g == _TRUE_) && (index_type == ppt->index_tp_g)) {
      
	/* newtonian gauge */
	if (ppr->gauge == newtonian) {
	  source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S0] = 
	    pvecmetric[ppw->index_mt_phi];
	  source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_dS1] = 0.;
	  source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_ddS2] = 0.;
	}

	/* synchronous gauge */
	if (ppr->gauge == synchronous) {
	  source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S0] = 
	    (a_prime_over_a * (pvecmetric[ppw->index_mt_h_prime] + 6. * pvecmetric[ppw->index_mt_eta_prime])/2./k2 + pvecmetric[ppw->index_mt_alpha_prime]);
	  source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_dS1] = 0.;
	  source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_ddS2] = 0.;

	  /* testing zone */
	  /* 	  if ((k>1.e-3) && (k<2.e-3)) */
	  /* 	    printf("%g %g %g %g %g %g\n", */
	  /* 		   k, */
	  /* 		   pvecback[pba->index_bg_a],    */
	  /* 		   a_prime_over_a*y[ppw->pv->index_pt_theta_b], */
	  /* 		   pvecthermo[pth->index_th_cb2]*k2*y[ppw->pv->index_pt_delta_b], */
	  /* 		   k2*y[ppw->pv->index_pt_eta], */
	  /* 		   0.5*a_prime_over_a*pvecmetric[ppw->index_mt_h_prime] */
	  /* 		   ); */
	 

	}

      }

      /* total matter overdensity */

      if ((ppt->has_source_delta_pk == _TRUE_) && (index_type == ppt->index_tp_delta_pk)) {

	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S0] = ppw->delta_pk;
	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_dS1] = 0.;
	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_ddS2] = 0.;

      }

      /* delta_g */
      if ((ppt->has_source_delta_g == _TRUE_) && (index_type == ppt->index_tp_delta_g)) {

	if (ppw->approx[ppw->index_ap_rsa]==(int)rsa_off) {
	  source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S0] = y[ppw->pv->index_pt_delta_g];
	}
	else {
	  source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S0] = ppw->rsa_delta_g;
	}
	
	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S0] = delta_g;
	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_dS1] = 0.;
	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_ddS2] = 0.;
      }

      /* delta_baryon */
      if ((ppt->has_source_delta_b == _TRUE_) && (index_type == ppt->index_tp_delta_b)) {
	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S0] = 
	  y[ppw->pv->index_pt_delta_b]; 
	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_dS1] = 0.;
	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_ddS2] = 0.;
      }

      /* delta_cdm */
      if ((ppt->has_source_delta_cdm == _TRUE_) && (index_type == ppt->index_tp_delta_cdm)) {
	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S0] = 
	  y[ppw->pv->index_pt_delta_cdm]; 
	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_dS1] = 0.;
	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_ddS2] = 0.;
      }

      /* delta_fld */
      if ((ppt->has_source_delta_fld == _TRUE_) && (index_type == ppt->index_tp_delta_fld)) {
	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S0] = 
	  y[ppw->pv->index_pt_delta_fld]; 
	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_dS1] = 0.;
	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_ddS2] = 0.;
      }

      /* delta_ur */
      if ((ppt->has_source_delta_ur == _TRUE_) && (index_type == ppt->index_tp_delta_ur)) {
	if (ppw->approx[ppw->index_ap_rsa]==(int)rsa_off) {
	  source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S0] =
	    y[ppw->pv->index_pt_delta_ur];
	}
	else {
	  source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S0] =
	    ppw->rsa_delta_ur;
	}
	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_dS1] = 0.;
	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_ddS2] = 0.;
      }
      
      /* delta_ncdm1 */
      if ((ppt->has_source_delta_ncdm == _TRUE_) && (index_type >= ppt->index_tp_delta_ncdm1) && (index_type < ppt->index_tp_delta_ncdm1+pba->N_ncdm)) {
	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S0] = 
	  ppw->delta_ncdm[index_type - ppt->index_tp_delta_ncdm1];
	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_dS1] = 0.;
	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_ddS2] = 0.;
      }
    }
  }

  /* tensors */
  if ((ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors)) {

    Psi=y[ppw->pv->index_pt_delta_g]/40.
      +y[ppw->pv->index_pt_shear_g]*2./35.
      +y[ppw->pv->index_pt_delta_g+4]/210.
      -y[ppw->pv->index_pt_pol0_g]*3./5. 
      +y[ppw->pv->index_pt_pol2_g]*6./35.
      -y[ppw->pv->index_pt_pol0_g+4]/210.;

    Psi_prime=dy[ppw->pv->index_pt_delta_g]/40.
      +dy[ppw->pv->index_pt_shear_g]*2./35.
      +dy[ppw->pv->index_pt_delta_g+4]/210.
      -dy[ppw->pv->index_pt_pol0_g]*3./5.  
      +dy[ppw->pv->index_pt_pol2_g]*6./35. 
      -dy[ppw->pv->index_pt_pol0_g+4]/210.;

    /** - for each type and each mode, compute S0, S1, S2 */
    for (index_type = 0; index_type < ppt->tp_size[index_mode]; index_type++) {

      source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_tau] = tau;

      /* tensor temperature */
      if ((ppt->has_source_t == _TRUE_) && (index_type == ppt->index_tp_t)) {

	if (x > 0.) {
	  source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S0] = 
	    (-y[ppw->pv->index_pt_gwdot]*pvecthermo[pth->index_th_exp_m_kappa]
	     +pvecthermo[pth->index_th_g]*Psi)/x/x;
	}
	else {
	  source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S0] = 0.;
	}

	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_dS1] = 0.;
	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_ddS2] = 0.;

      }

      /* tensor polarization */
      if ((ppt->has_source_e == _TRUE_) && (index_type == ppt->index_tp_e)) {

	if (x > 0.) {

	  source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S0] = 
	    (1.-2./x/x)*pvecthermo[pth->index_th_g]*Psi;
	  source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_dS1] = 
	    -4./x*(pvecthermo[pth->index_th_g]*(Psi/x+Psi_prime/k)+pvecthermo[pth->index_th_dg]*Psi/k);
	  source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_dS2] = 
	    -(pvecthermo[pth->index_th_g]*Psi_prime+pvecthermo[pth->index_th_dg]*Psi)/k/k;

	}

	else {

	  source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S0] = 0.;
	  source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_dS1] = 0.;
	  source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_dS2] = 0.;

	}

      }

      if ((ppt->has_source_b == _TRUE_) && (index_type == ppt->index_tp_b)) {

	if (x > 0.) {
	
	  source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S0] = 
	    -pvecthermo[pth->index_th_g]*(4.*Psi/x+2.*Psi_prime/k)-2.*pvecthermo[pth->index_th_dg]*Psi/k;
	}
	
	else {
	  source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S0] = 0.;
	}
	
	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_dS1] = 0.;
	source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_ddS2] = 0.;
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
 * @param ppt        Input : pointer to perturbation structure containing interpolation tables
 * @param index_mode Input : index of requested mode
 * @param index_ic   Input : index of requested initial condition
 * @param index_k    Input : index of requested wavenumber
 * @param ppw        Input/Output : pointer to workspace (contains source terms in input, sources in output)
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

  int index_tau,index_type;

  for (index_type = 0; index_type < ppt->tp_size[index_mode]; index_type++) {

    /** - for scalar temperature, infer \f$ S_2'' \f$ from \f$ S_2' \f$ at each time with array_derive1_order2_table_line_to_line() */

    if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {
      if ((ppt->has_source_t == _TRUE_) && (index_type == ppt->index_tp_t)) {

	/* before computing numerical derivatives, slice out the end of the table if filled with zeros */
	index_tau = ppt->tau_size-1;
	while ((ppw->source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_dS2] == 0.) && (index_tau > 0))
	  index_tau--;
	
	/* numerical derivative */
	class_call(array_derive1_order2_table_line_to_line(
							   ppt->tau_sampling,
							   index_tau+1,
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
	index_tau = ppt->tau_size-1;
	while ((ppw->source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_dS2] == 0.) && (index_tau > 0))
	  index_tau--;
	
	/* numerical derivative */
	class_call(array_derive1_order2_table_line_to_line(
							   ppt->tau_sampling,
							   index_tau+1,
							   ppw->source_term_table[index_type],
							   ppw->st_size,
							   ppw->index_st_dS2,
							   ppw->index_st_ddS2,
							   ppt->error_message),
		   ppt->error_message,
		   ppt->error_message);
      }
    }

    /** - for each time, sum up \f$ S = S_0 + S_1' + S_2'' \f$ and store in array ((sources[index_mode])[index_ic][index_type])[index_tau][index_k] */

    for (index_tau = 0; index_tau < ppt->tau_size; index_tau++) {

      ppt->sources[index_mode]
	[index_ic * ppt->tp_size[index_mode] + index_type]
	[index_tau * ppt->k_size[index_mode] + index_k] = 
	ppw->source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_S0]
	+ppw->source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_dS1]
	+ppw->source_term_table[index_type][index_tau * ppw->st_size + ppw->index_st_ddS2];

    }

  } /* end of loop over types */

  return _SUCCESS_;
}

/**
 * When testing the code or a cosmological model, it can be useful to
 * output perturbations at each step of integration (and not just the
 * delta's at each source sampling point, which is acheived simply by
 * asking for matter transfer functions). Then this function can be
 * passed to the generic_evolver routine.
 *
 * By default, instead of passing this function to generic_evolver,
 * one passes a null pointer. Then this function is just not used.
 *
 * @param tau                      Input: conformal time
 * @param y                        Input: vector of perturbations
 * @param dy                       Input: vector of its derivatives (already allocated)
 * @param parameters_and_workspace Input: fixed parameters (e.g. indices)
 * @param error_message            Output : error message
 *
 */

int perturb_print_variables(double tau,
			    double * y,
			    double * dy,
			    void * parameters_and_workspace,
			    ErrorMsg error_message
			    ) {
  
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
  

  /** - rename structure fields (just to avoid heavy notations) */

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

  /** - print whatever you want for whatever mode of your choice */

  double delta_g,theta_g,shear_g,l3_g,pol0_g,pol1_g,pol2_g,pol3_g;
  double delta_ur=0.,theta_ur=0.,shear_ur=0.;
  //double delta_ncdm=0.,theta_ncdm=0.,shear_ncdm=0.;

  if (pba->has_ur == _TRUE_) {
    if (ppw->approx[ppw->index_ap_rsa]==(int)rsa_off) {
      delta_ur = y[ppw->pv->index_pt_delta_ur];
      theta_ur = y[ppw->pv->index_pt_theta_ur];
      shear_ur = y[ppw->pv->index_pt_shear_ur];
    }
    else {
      delta_ur = ppw->rsa_delta_ur;
      theta_ur = ppw->rsa_theta_ur;
      shear_ur = 0.;
    }
  }
  
  /*  if((k>0.1)&&(k<0.12)){ */
  /*     printf("%g %g %g %g %g %g %g %g\n",k,tau, */
  /* 	   ppw->delta_ncdm[0],ppw->theta_ncdm1,ppw->shear_ncdm1, */
  /* 	   delta_ur,theta_ur,shear_ur); */
  /*   } */
  

  if ((k>=0.069) && (k<0.071)) {

    if (ppw->approx[ppw->index_ap_rsa]==(int)rsa_off) {
      delta_g = y[ppw->pv->index_pt_delta_g];
      theta_g = y[ppw->pv->index_pt_theta_g];
    }
    else {
      delta_g = ppw->rsa_delta_g;
      theta_g = ppw->rsa_theta_g;
    }
    
    if (ppw->approx[ppw->index_ap_rsa]==(int)rsa_off) {
      if (ppw->approx[ppw->index_ap_tca]==(int)tca_on) {
	shear_g = ppw->tca_shear_g;
	l3_g = 6./7.*k/ppw->pvecthermo[pth->index_th_dkappa]*ppw->tca_shear_g;
	pol0_g = 2.5*ppw->tca_shear_g;
	pol1_g = 7./12.*6./7.*k/ppw->pvecthermo[pth->index_th_dkappa]*ppw->tca_shear_g;
	pol2_g = 0.5*ppw->tca_shear_g;
	pol3_g = 0.25*6./7.*k/ppw->pvecthermo[pth->index_th_dkappa]*ppw->tca_shear_g;
      }
      else {
	shear_g = y[ppw->pv->index_pt_shear_g];
	l3_g = y[ppw->pv->index_pt_l3_g];
	pol0_g = y[ppw->pv->index_pt_pol0_g];
	pol1_g = y[ppw->pv->index_pt_pol1_g];
	pol2_g = y[ppw->pv->index_pt_pol2_g];
	pol3_g = y[ppw->pv->index_pt_pol3_g];
      }
    }
    else {
      shear_g = 0;
      l3_g = 0;
      pol0_g = 0;
      pol1_g = 0;
      pol2_g = 0;
      pol3_g = 0.;
    }
    
    if (pba->has_ur == _TRUE_) {
      if (ppw->approx[ppw->index_ap_rsa]==(int)rsa_off) {
	delta_ur = y[ppw->pv->index_pt_delta_ur];
	theta_ur = y[ppw->pv->index_pt_theta_ur];
	shear_ur = y[ppw->pv->index_pt_shear_ur];
      }
      else {
	delta_ur = ppw->rsa_delta_ur;
	theta_ur = ppw->rsa_theta_ur;
	shear_ur = 0.;
      }
    }

    fprintf(stdout,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
	    k,
	    tau,
	    delta_g,
	    theta_g,
	    shear_g,
	    l3_g,
	    pol0_g,
	    pol1_g,
	    pol2_g,
	    pol3_g,
	    y[ppw->pv->index_pt_delta_b],
	    y[ppw->pv->index_pt_theta_b],
	    delta_ur,
	    theta_ur,
	    shear_ur,
	    ppw->delta_ncdm[0],
	    ppw->theta_ncdm1,
	    ppw->shear_ncdm1,
	    y[ppw->pv->index_pt_eta],
	    pvecmetric[ppw->index_mt_eta_prime],
	    pvecmetric[ppw->index_mt_h_prime]);

  }
 
  /** - print whatever you want for whatever mode of your choice */
  /*   qsiz = ppw->pv->q_size_ncdm1; */
  /*   if(((k>0.004)&&(k<0.005))||((k>0.0148)&&(k<0.0152))){ */
  /*     fprintf(stdout,"%g %g ",pvecback[pba->index_bg_a],k); */
  /*     for(l=0; l<=ppw->pv->l_max_ncdm1; l++){ */
  /*       for(index_q=0; index_q < qsiz; index_q++){ */
  /*       index_q = 3; */
  /*       fprintf(stdout,"%g ",y[ppw->pv->index_pt_psi0_ncdm1+l*qsiz+index_q]); */
  /*       } */
  /*     } */
  /*     fprintf(stdout,"\n"); */
  /* } */


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
 * @param tau                      Input: conformal time
 * @param y                        Input: vector of perturbations
 * @param dy                       Ouput: vector of its derivatives (already allocated)
 * @param parameters_and_workspace Input/Output: in input, fixed parameters (e.g. indices); in output, background and thermo quantities evaluated at tau.
 * @param error_message            Output : error message
 */

int perturb_derivs(double tau,
		   double * y,
		   double * dy,
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
  double R;//fracnu;

  /* useful terms for tight-coupling approximation */
  double slip=0.;
  double Pi;
  double tau_c=0.,dtau_c=0.;
  double theta_prime,shear_g_prime=0.,theta_prime_prime;
  double g0,g0_prime,g0_prime_prime;
  double F=0.,F_prime=0.,F_prime_prime=0.;

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
  double * pvecback;
  double * pvecthermo;
  double * pvecmetric;

  /* short-cut notations for the perturbations */
  double delta_g=0.,theta_g=0.,shear_g=0.;
  double delta_b,theta_b;
  double Delta;
  double cb2;

  /* For use with non-cold Dark Matter: */
  int index_q,n_ncdm,idx;
  double q,epsilon,dlnf0_dlnq,qk_div_epsilon;
  double rho_ncdm_bg,p_ncdm_bg,pseudo_p_ncdm,w_ncdm,ca2_ncdm,ceff2_ncdm=0.,cvis2_ncdm=0.;

  /** - rename structure fields (just to avoid heavy notations) */

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

  /** - get background/thermo quantities in this point */

  class_call(background_at_tau(pba,
			       tau, 
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
  
  /** - for scalar mode: */
  if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

    /* short-cut notations for the scalar perturbations */
    if (ppw->approx[ppw->index_ap_rsa] == (int)rsa_off) {
      delta_g = y[ppw->pv->index_pt_delta_g];
      theta_g = y[ppw->pv->index_pt_theta_g];
      if (ppw->approx[ppw->index_ap_tca] == (int)tca_off) {
	shear_g = y[ppw->pv->index_pt_shear_g];
      }
    }
    delta_b = y[ppw->pv->index_pt_delta_b];
    theta_b = y[ppw->pv->index_pt_theta_b];
    cb2 = pvecthermo[pth->index_th_cb2];
    
    /* short-cut notations used only in tight-coupling approximation */
    if (ppw->approx[ppw->index_ap_tca] == (int)tca_on) {
      tau_c = 1./pvecthermo[pth->index_th_dkappa]; /* inverse of opacity */
      dtau_c = -pvecthermo[pth->index_th_ddkappa]*tau_c*tau_c; /* its first derivative wrt conformal time */
      F = tau_c/(1+R); /* F = tau_c/(1+R) */
      if (ppr->tight_coupling_approximation >= (int)second_order_CLASS) {
	F_prime = dtau_c/(1+R)+tau_c*a_prime_over_a*R/(1+R)/(1+R); /*F' needed by second_order_CLASS and compromise_CLASS */
	if (ppr->tight_coupling_approximation == (int)second_order_CLASS) {
	  F_prime_prime =(- pvecthermo[pth->index_th_dddkappa]*tau_c*tau_c /* F'' needed by second_order_CLASS only */
			  + 2.*pvecthermo[pth->index_th_ddkappa]*pvecthermo[pth->index_th_ddkappa]*tau_c*tau_c*tau_c)/(1+R)
	    +2.*dtau_c*a_prime_over_a*R/(1+R)/(1+R)
	    +tau_c*((a_primeprime_over_a-2.*a_prime_over_a*a_prime_over_a)+2.*a_prime_over_a*a_prime_over_a*R/(1+R))*R/(1+R)/(1+R);
	}
      }
    }

    /** (a) get metric perturbations with perturb_einstein() */
    class_call(perturb_einstein(ppr,
				pba,
				pth,
				ppt,
				index_mode,
				k,
				tau,
				y,
				ppw),
	       ppt->error_message,
	       error_message);

    /* compute metric-related quantities */
    h_plus_six_eta_prime = pvecmetric[ppw->index_mt_h_prime] + 6. * pvecmetric[ppw->index_mt_eta_prime];

    /** (b) if some approximation schemes are turned on, enforce a few y[] values computed in perturb_einstein */

    /* free-streaming photon velocity */
    if (ppw->approx[ppw->index_ap_rsa] == (int)rsa_on) 
      theta_g = ppw->rsa_theta_g;

    /** (c) Photon temperature density (does not depend on tca) */

    if (ppw->approx[ppw->index_ap_rsa] == (int)rsa_off) {

      if (ppr->gauge == newtonian)
	dy[ppw->pv->index_pt_delta_g] = /* photon density */
	  -4./3.*theta_g+4.*pvecmetric[ppw->index_mt_phi_prime];
	  
      if (ppr->gauge == synchronous)
	dy[ppw->pv->index_pt_delta_g] = /* photon density */
	  -4./3.*theta_g - 2./3.*pvecmetric[ppw->index_mt_h_prime];
      
    }

    /** (c) baryon density (does not depend on tca) */

    if (ppr->gauge == newtonian)
      dy[ppw->pv->index_pt_delta_b] = /* baryon density */
	-theta_b + 3.*pvecmetric[ppw->index_mt_phi_prime];
    
    if (ppr->gauge == synchronous)
      dy[ppw->pv->index_pt_delta_b] = /* baryon density */
	-theta_b - 0.5*pvecmetric[ppw->index_mt_h_prime];
    
    /** (d) Baryon velocity (depends on tight-coupling approximation) */

    /** (d.1) Baryon velocity if tight-coupling is off */

    if (ppw->approx[ppw->index_ap_tca] == (int)tca_off) {
      
      if (ppr->gauge == newtonian)
	/* Newtonian gauge : */
	dy[ppw->pv->index_pt_theta_b] = /* baryon velocity */
	  - a_prime_over_a*theta_b 
	  + k2*pvecmetric[ppw->index_mt_psi] 
	  + cb2*k2*delta_b
	  + R*pvecthermo[pth->index_th_dkappa]*(theta_g-theta_b);

      if (ppr->gauge == synchronous)
	/* Synchronous gauge : */
	dy[ppw->pv->index_pt_theta_b] = /* baryon velocity */
	  - a_prime_over_a*theta_b
	  + cb2*k2*delta_b
	  + R*pvecthermo[pth->index_th_dkappa]*(theta_g-theta_b);
    
    }
  
    /** (d.2) Baryon velocity if baryon tight-coupling is on */
    
    else {

      if (ppr->gauge == newtonian) {
	/* Newtonian gauge : */
	slip=(2.*R/(1.+R)*a_prime_over_a+pvecthermo[pth->index_th_ddkappa]/pvecthermo[pth->index_th_dkappa]) /* tight-coupling (theta_b-theta_g)' */
	  *(theta_b-theta_g)
	  +(-a_primeprime_over_a*theta_b
	    -a_prime_over_a*k2*(y[ppw->pv->index_pt_delta_g]/2.+pvecmetric[ppw->index_mt_psi])
	    +k2*(cb2*dy[ppw->pv->index_pt_delta_b]
		 -dy[ppw->pv->index_pt_delta_g]/4.)
	    )/pvecthermo[pth->index_th_dkappa]/(1.+R);

	shear_g=(8./3.*theta_g)*2./15./pvecthermo[pth->index_th_dkappa]; /* tight-coupling shear_g (Ma & Bertschinger give 1/9 instead of 2/15 becasue they didn't include consistently the contribution of G_gamma0 and G_gamma2, which are of the same order as sigma_g. This was already consistently included in CAMB)*/

	dy[ppw->pv->index_pt_theta_b] = /* tight-coupling baryon velocity */
	  (-a_prime_over_a*theta_b
	   +cb2*k2*delta_b
	   +k2*R*(y[ppw->pv->index_pt_delta_g]/4.-shear_g)
	   +R*slip)/(1.+R)
	  +k2*pvecmetric[ppw->index_mt_psi];
      }

      if (ppr->gauge == synchronous) {

	/* Synchronous gauge : */

	/* tight-coupling slip = (theta_b-theta_g)' */

	/* like Ma & Bertschinger */
	if (ppr->tight_coupling_approximation == (int)first_order_MB) {

	  slip=2.*R/(1.+R)*a_prime_over_a*(theta_b-theta_g)
	    +F*(-a_primeprime_over_a*theta_b
		+k2*(-a_prime_over_a*delta_g/2.
		     +cb2*dy[ppw->pv->index_pt_delta_b]
		     -dy[ppw->pv->index_pt_delta_g]/4.));
	  
	}

	/* relax assumption dkappa~a^-2 (like in CAMB) */
	if ((ppr->tight_coupling_approximation == (int)first_order_CAMB) || (ppr->tight_coupling_approximation == (int)compromise_CLASS)) {

	  slip=(dtau_c/tau_c-2.*a_prime_over_a/(1.+R))*(theta_b-theta_g)
	    +F*(-a_primeprime_over_a*theta_b
		+k2*(-a_prime_over_a*delta_g/2.
		     +cb2*dy[ppw->pv->index_pt_delta_b]
		     -dy[ppw->pv->index_pt_delta_g]/4.));
	  
	}

	/* also relax assumption cb2~a^-1 */
	if ((ppr->tight_coupling_approximation == (int)first_order_CLASS) || (ppr->tight_coupling_approximation == (int)second_order_CLASS)){
	  
	  slip=(dtau_c/tau_c-2.*a_prime_over_a/(1.+R))*(theta_b-theta_g)
	    +F*(-a_primeprime_over_a*theta_b
		+k2*(-a_prime_over_a*delta_g/2.
		     +pvecthermo[pth->index_th_dcb2]*y[ppw->pv->index_pt_delta_b]
		     +cb2*dy[ppw->pv->index_pt_delta_b]
		     -dy[ppw->pv->index_pt_delta_g]/4.));
	}

	/* shear_g at first order in tight-coupling */
	shear_g=8./45.*tau_c*(2.*theta_g+h_plus_six_eta_prime);
	/* (Ma & Bertschinger give (1/9)*(4/3) instead of (2/15)*(4/3)
	   becasue they didn't include the contribution of G_gamma0
	   and G_gamma2, which are of the same order as sigma_g. This
	   was already consistently included in CAMB) */ 

	/* zero order for theta_b' = theta_g' */
	theta_prime = (-a_prime_over_a*theta_b+k2*(cb2*delta_b+R/4.*delta_g))/(1.+R);
	
	/* shear_g_prime at first order in tight-coupling */
	shear_g_prime=16./45.*(tau_c*(theta_prime+k2*pvecmetric[ppw->index_mt_alpha_prime])
			       +dtau_c*(theta_g+0.5*h_plus_six_eta_prime));

	/* 2nd order as in CRS*/
	if (ppr->tight_coupling_approximation == (int)second_order_CRS) {

	  /* infer Delta from h'' using Einstein equation */

	  Delta = 2*k2*y[ppw->pv->index_pt_eta]
	    -2*a_prime_over_a*pvecmetric[ppw->index_mt_h_prime]
	    -pvecmetric[ppw->index_mt_h_prime_prime];

	  /* monster expression for slip at second-order in tight-coupling */
	  slip=(-2./(1.+R)*a_prime_over_a-pvecthermo[pth->index_th_ddkappa]/pvecthermo[pth->index_th_dkappa])*(theta_b-theta_g)
	    +(-a_primeprime_over_a*theta_b
	      -k2*a_prime_over_a*(delta_g/2.-2.*shear_g)
	      +k2*(cb2*dy[ppw->pv->index_pt_delta_b]
		   -dy[ppw->pv->index_pt_delta_g]/4.
		   +shear_g_prime)
	      )/pvecthermo[pth->index_th_dkappa]/(1.+R)
	    -2.*R*(3.*a_prime_over_a*a_prime_over_a*cb2+(1.+R)*(a_primeprime_over_a-a_prime_over_a*a_prime_over_a)-3.*a_prime_over_a*a_prime_over_a)
	    /(1.+R)/(1.+R)/(1.+R)*(theta_b-theta_g)/pvecthermo[pth->index_th_dkappa]
	    +(
	      a_primeprime_over_a*a_prime_over_a*((2.-3.*cb2)*R-2.)*theta_b/(1.+R)
	      +a_prime_over_a*k2*(1.-3.*cb2)*theta_b/3./(1.+R)
	      +a_primeprime_over_a*k2*cb2*delta_b/(1.+R)
	      +k2*k2*(3.*cb2-1.)*cb2*delta_b/3./(1.+R)
	      +k2*k2*R*(3.*cb2-1.)*delta_g/12./(1.+R)
	      +a_primeprime_over_a*k2*(2.+3.*R)*delta_g/4./(1.+R)
	      +a_prime_over_a*a_prime_over_a*k2*((2.-3.*cb2)*R-1.)*delta_g/2./(1.+R)
	      +a_prime_over_a*k2*cb2*(1.+(3.*cb2-2.)*R)*dy[ppw->pv->index_pt_delta_b]/(1.+R)
	      +a_prime_over_a*k2*(2.+(5.-3.*cb2)*R)*dy[ppw->pv->index_pt_delta_g]/4./(1.+R)
	      +a_prime_over_a*(1.-3.*cb2)*k2*h_plus_six_eta_prime/3.
	      +k2*k2*(3.*cb2-1.)*y[ppw->pv->index_pt_eta]/3.
	      +2.*a_prime_over_a*k2*(3.*cb2-1.)*pvecmetric[ppw->index_mt_eta_prime]
	      +k2*(1.-3.*cb2)*Delta/6.
	      )/pvecthermo[pth->index_th_dkappa]/pvecthermo[pth->index_th_dkappa]/(1.+R)/(1.+R)
	    -(4.*a_primeprime_over_a*theta_b-4.*k2*cb2*dy[ppw->pv->index_pt_delta_b]+2.*a_prime_over_a*k2*delta_g+k2*dy[ppw->pv->index_pt_delta_g])/2./(1.+R)/(1.+R)*pvecthermo[pth->index_th_ddkappa]/pvecthermo[pth->index_th_dkappa]/pvecthermo[pth->index_th_dkappa]/pvecthermo[pth->index_th_dkappa]
	    +4.*a_prime_over_a*R/(1.+R)/(1.+R)*pvecthermo[pth->index_th_ddkappa]/pvecthermo[pth->index_th_dkappa]/pvecthermo[pth->index_th_dkappa]*(theta_b-theta_g);

	  /* second-order correction to shear */
	  shear_g = (1.-11./6.*dtau_c)*shear_g-11./6.*tau_c*16./45.*tau_c*(theta_prime+k2*pvecmetric[ppw->index_mt_alpha_prime]); 

	}

	/* 2nd order like in CLASS paper */
	if (ppr->tight_coupling_approximation == (int)second_order_CLASS) {

	  /* zero order for theta_b'' = theta_g'' */
	  theta_prime_prime = ((R-1.)*a_prime_over_a*theta_prime-(a_primeprime_over_a-a_prime_over_a*a_prime_over_a)*theta_b
			       +k2*(pvecthermo[pth->index_th_dcb2]*delta_b+cb2*dy[ppw->pv->index_pt_delta_b]-a_prime_over_a*R/4.*delta_g+R/4.*dy[ppw->pv->index_pt_delta_g]))/(1.+R);

	  /* zero-order quantities g0, g0', go'' */
	  g0 = -a_prime_over_a*theta_b + k2*(cb2*delta_b-delta_g/4.);
	  g0_prime = -a_prime_over_a*theta_prime-(a_primeprime_over_a-a_prime_over_a*a_prime_over_a)*theta_b+k2*(pvecthermo[pth->index_th_dcb2]*delta_b+(1./3.-cb2)*(theta_b+0.5*pvecmetric[ppw->index_mt_h_prime]));
	  g0_prime_prime = -a_prime_over_a*theta_prime_prime-2.*(a_primeprime_over_a-a_prime_over_a*a_prime_over_a)*theta_prime
	    -(2.*a_prime_over_a*a_prime_over_a*a_prime_over_a-3.*a_primeprime_over_a*a_prime_over_a)*theta_b
	    +k2*(pvecthermo[pth->index_th_ddcb2]*delta_b-2.*pvecthermo[pth->index_th_dcb2]*(theta_b+0.5*pvecmetric[ppw->index_mt_h_prime])+(1./3.-cb2)*(theta_prime+0.5*pvecmetric[ppw->index_mt_h_prime_prime]));

	  /* slip at second order */
	  slip = (1.-2*a_prime_over_a*F)*slip + F*k2*(2.*a_prime_over_a*shear_g+shear_g_prime)
	    -F*(F_prime_prime*g0+2.*F_prime*g0_prime+F*g0_prime_prime);
	  
	  /* second-order correction to shear */
	  shear_g = (1.-11./6.*dtau_c)*shear_g-11./6.*tau_c*16./45.*tau_c*(theta_prime+k2*pvecmetric[ppw->index_mt_alpha_prime]); 

	}

	/* add only the most important 2nd order terms */
	if (ppr->tight_coupling_approximation == (int)compromise_CLASS) {
	  
	  /* slip at second order (only leading second-order terms) */
	  slip = (1.-2.*a_prime_over_a*F)*slip + F*k2*(2.*a_prime_over_a*shear_g+shear_g_prime-(1./3.-cb2)*(F*theta_prime+2.*F_prime*theta_b));
		
	  /* second-order correction to shear */
	  shear_g = (1.-11./6.*dtau_c)*shear_g-11./6.*tau_c*16./45.*tau_c*(theta_prime+k2*pvecmetric[ppw->index_mt_alpha_prime]); 

	}

	dy[ppw->pv->index_pt_theta_b] = /* tight-coupling baryon velocity */
	  (-a_prime_over_a*theta_b+k2*(cb2*delta_b+R*(delta_g/4.-shear_g))+R*slip)/(1.+R);
      }

      ppw->tca_shear_g = shear_g;
      ppw->tca_shear_g_prime = shear_g_prime;

    }

    /** (e) Photon temperature higher momenta and photon polarisation (depend on tight-coupling approximation) : */

    if (ppw->approx[ppw->index_ap_rsa] == (int)rsa_off) {

      /** (e.1) if photon tight-coupling is off: */ 
      if (ppw->approx[ppw->index_ap_tca] == (int)tca_off) {

	/** (e.1.a) define \f$ \Pi = G_{\gamma 0} + G_{\gamma 2} + F_{\gamma 2} \f$ */
	Pi = y[ppw->pv->index_pt_pol0_g] + y[ppw->pv->index_pt_pol2_g] + 2.*shear_g;

	/** (e.1.b) Photon velocity and shear */ 

	if (ppr->gauge == newtonian) {
	  /* Newtonian gauge : */
	  dy[ppw->pv->index_pt_theta_g] = /* photon velocity */
	    k2*(delta_g/4.
		-shear_g+pvecmetric[ppw->index_mt_psi])
	    +pvecthermo[pth->index_th_dkappa]*(theta_b
					       -theta_g);

	  dy[ppw->pv->index_pt_shear_g] = /* photon shear */
	    0.5*(8./15.*theta_g
		 -3./5.*k*y[ppw->pv->index_pt_l3_g]
		 -pvecthermo[pth->index_th_dkappa]*(2.*shear_g-1./10.*Pi));
	}
      
	if (ppr->gauge == synchronous) {
	  /* Synchronous gauge : */
	  dy[ppw->pv->index_pt_theta_g] = /* photon velocity */
	    k2*(delta_g/4.
		-shear_g)
	    + pvecthermo[pth->index_th_dkappa]*(theta_b
						-theta_g);
	  
	  dy[ppw->pv->index_pt_shear_g] = /* photon shear */
	    0.5*(8./15.*theta_g
		 -3./5.*k*y[ppw->pv->index_pt_l3_g]
		 +4./15.*pvecmetric[ppw->index_mt_h_prime]+8./5.*pvecmetric[ppw->index_mt_eta_prime]
		 -pvecthermo[pth->index_th_dkappa]*(2.*shear_g-1./10.*Pi));
	}
      
	/** (e.1.c) Photon temperature higher momenta (l >=3), gauge-independent */ 

	l = 3; /* photon l=3 (special case because F_gamma2=2*shear !!) */
	dy[ppw->pv->index_pt_l3_g] =
	  k/(2.*l+1.)*(l*2.*shear_g-(l+1.)*y[ppw->pv->index_pt_l3_g+1])
	  - pvecthermo[pth->index_th_dkappa]*y[ppw->pv->index_pt_l3_g];

	for (l = 4; l < ppw->pv->l_max_g; l++) { /* photon additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3) */
	  dy[ppw->pv->index_pt_delta_g+l] =
	    k/(2.*l+1)*(l*y[ppw->pv->index_pt_delta_g+l-1]-(l+1.)*y[ppw->pv->index_pt_delta_g+l+1])
	    - pvecthermo[pth->index_th_dkappa]*y[ppw->pv->index_pt_delta_g+l];
	}

	l = ppw->pv->l_max_g; /* l=lmax */
	dy[ppw->pv->index_pt_delta_g+ppw->pv->l_max_g] = /* last photon term */
	  k*y[ppw->pv->index_pt_delta_g+ppw->pv->l_max_g-1]
	  -(1.+l)/tau*y[ppw->pv->index_pt_delta_g+ppw->pv->l_max_g]
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
	  k*y[ppw->pv->index_pt_pol0_g+l-1]-(l+1)/tau*y[ppw->pv->index_pt_pol0_g+l]
	  -pvecthermo[pth->index_th_dkappa]*y[ppw->pv->index_pt_pol0_g+l];

      }

      /** (e.2) if photon tight-coupling is on: */
 
      else {

	/** (e.2.a) photon velocity */

	if (ppr->gauge == newtonian)
	  /* Newtonian gauge : */
	  dy[ppw->pv->index_pt_theta_g] = /* tight-coupling photon velocity */
	    -(dy[ppw->pv->index_pt_theta_b]+a_prime_over_a*theta_b-cb2*k2*delta_b)/R
	    +k2*(0.25*delta_g-shear_g)+(1.+R)/R*k2*pvecmetric[ppw->index_mt_psi];
	
	if (ppr->gauge == synchronous)
	  /* Synchronous gauge : */
	  
	  dy[ppw->pv->index_pt_theta_g] = /* tight-coupling photon velocity */
	    -(dy[ppw->pv->index_pt_theta_b]+a_prime_over_a*theta_b-cb2*k2*delta_b)/R
	    +k2*(delta_g/4.-shear_g);
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
    
    /** (g) fluid */
    
    if (pba->has_fld == _TRUE_) {  

      if (ppr->gauge == newtonian) {
	/* Newtonian gauge : */
	dy[ppw->pv->index_pt_delta_fld] = /* fluid density */
	  (-3*(1+ pba->w_fld )*a_prime_over_a-3*pvecback[pba->index_bg_H]*(pba->cs2_fld- pba->w_fld )*(y[ppw->pv->index_pt_delta_fld]/pvecback[pba->index_bg_rho_fld]+3*pvecback[pba->index_bg_H]*(1+ pba->w_fld )*y[ppw->pv->index_pt_theta_fld]/k)-(1+ pba->w_fld )*k*y[ppw->pv->index_pt_theta_fld])/pvecback[pba->index_bg_rho_fld]; // 0;

	dy[ppw->pv->index_pt_theta_fld] = /* fluid velocity */
	  (k*pba->cs2_fld*y[ppw->pv->index_pt_delta_fld])/(pvecback[pba->index_bg_rho_fld]*(1+ pba->w_fld ))-pvecback[pba->index_bg_H]*(1-3*pba->cs2_fld)*y[ppw->pv->index_pt_theta_fld]+k*pvecmetric[ppw->index_mt_psi]; // 0;
      }

      if (ppr->gauge == synchronous) {
	/* Synchronous gauge : */
	dy[ppw->pv->index_pt_delta_fld] = /* fluid density */
	  -(1+ pba->w_fld )*(y[ppw->pv->index_pt_theta_fld]+0.5*pvecmetric[ppw->index_mt_h_prime])
	  -3.*(pba->cs2_fld-pba->w_fld)*a_prime_over_a*y[ppw->pv->index_pt_delta_fld]
	  -9.*(1+ pba->w_fld )*(pba->cs2_fld-pba->w_fld)*a_prime_over_a*a_prime_over_a*y[ppw->pv->index_pt_theta_fld]/k2;

	dy[ppw->pv->index_pt_theta_fld] = /* fluid velocity */
	  -(1.-3.*pba->cs2_fld)*a_prime_over_a*y[ppw->pv->index_pt_theta_fld]
	  +pba->cs2_fld*k2/(1.+pba->w_fld)*y[ppw->pv->index_pt_delta_fld];
      }
      
    }  
    
    /** (h) ultra-relativistic neutrino/relics density, velocity, shear, etc. */
    
    if (pba->has_ur == _TRUE_) {
      
      if (ppw->approx[ppw->index_ap_rsa] == (int)rsa_off) {

	if (ppr->gauge == newtonian) {
	  
	  /* Newtonian gauge : */
	  dy[ppw->pv->index_pt_delta_ur] = /* density of ultra-relativistic neutrinos/relics */
	    -4./3.*y[ppw->pv->index_pt_theta_ur] + 4.*pvecmetric[ppw->index_mt_phi_prime];
	  
	  dy[ppw->pv->index_pt_theta_ur] = /* velocity of ultra-relativistic neutrinos/relics */
	    k2*(y[ppw->pv->index_pt_delta_ur]/4.
		-y[ppw->pv->index_pt_shear_ur]+pvecmetric[ppw->index_mt_psi]);
	}
	
	if (ppr->gauge == synchronous) {
	  
	  /* Synchronous gauge : */
	  dy[ppw->pv->index_pt_delta_ur] = /* density of ultra-relativistic neutrinos/relics */
	    -4./3.*y[ppw->pv->index_pt_theta_ur] - 2./3.*pvecmetric[ppw->index_mt_h_prime];
	  
	  dy[ppw->pv->index_pt_theta_ur] = /* velocity of ultra-relativistic neutrinos/relics */
	    k2*(y[ppw->pv->index_pt_delta_ur]/4.
		-y[ppw->pv->index_pt_shear_ur]);
	}
      
	if(ppw->approx[ppw->index_ap_ufa] == (int)ufa_off) {
	  
	  if (ppr->gauge == newtonian) {
	    dy[ppw->pv->index_pt_shear_ur] = /* shear of ultra-relativistic neutrinos/relics */
	      0.5*(8./15.*y[ppw->pv->index_pt_theta_ur]
		   -3./5.*k*y[ppw->pv->index_pt_shear_ur+1]);
	  }
	  
	  if (ppr->gauge == synchronous) {
	    dy[ppw->pv->index_pt_shear_ur] = /* shear of ultra-relativistic neutrinos/relics */
	      0.5*(8./15.*y[ppw->pv->index_pt_theta_ur]
		   -3./5.*k*y[ppw->pv->index_pt_shear_ur+1]
		   +4./15.*pvecmetric[ppw->index_mt_h_prime]+8./5.*pvecmetric[ppw->index_mt_eta_prime]);
	  }
	  
	  l = 3;
	  dy[ppw->pv->index_pt_l3_ur] = /* l=3 of ultra-relativistic neutrinos/relics (special case because F_gamma2=2*shear !!) */
	    k/(2.*l+1.)*(l*2.*y[ppw->pv->index_pt_shear_ur]-(l+1.)*y[ppw->pv->index_pt_l3_ur+1]);
	  
	  for (l = 4; l < ppw->pv->l_max_ur; l++) {
	    dy[ppw->pv->index_pt_delta_ur+l] = /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3) */
	      k/(2.*l+1)*(l*y[ppw->pv->index_pt_delta_ur+l-1]-(l+1.)*y[ppw->pv->index_pt_delta_ur+l+1]);
	  }
	  
	  l = ppw->pv->l_max_ur; /* l=lmax */
	  dy[ppw->pv->index_pt_delta_ur+ppw->pv->l_max_ur] = /* last term of ultra-relativistic neutrinos/relics */
	    k*y[ppw->pv->index_pt_delta_ur+ppw->pv->l_max_ur-1]
	    -(1.+l)/tau*y[ppw->pv->index_pt_delta_ur+ppw->pv->l_max_ur];
	  
	}
	
	else {
	  
	  if (ppr->gauge == newtonian) { }
	  
	  if (ppr->gauge == synchronous) {

	    /* shear of ultra-relativistic neutrinos/relics in fluid approach */
	    if (ppr->ur_fluid_approximation == ufa_mb) {
	      
	      dy[ppw->pv->index_pt_shear_ur] =
		-3./tau*y[ppw->pv->index_pt_shear_ur]
		+2./3.*y[ppw->pv->index_pt_theta_ur]
		+1./3.*(pvecmetric[ppw->index_mt_h_prime]+6.*pvecmetric[ppw->index_mt_eta_prime]);
	      
	    }

	    if (ppr->ur_fluid_approximation == ufa_hu) {

	      dy[ppw->pv->index_pt_shear_ur] =
		-3.*a_prime_over_a*y[ppw->pv->index_pt_shear_ur]
		+2./3.*y[ppw->pv->index_pt_theta_ur]
		+1./3.*(pvecmetric[ppw->index_mt_h_prime]+6.*pvecmetric[ppw->index_mt_eta_prime]);
	      
	    }

	    if (ppr->ur_fluid_approximation == ufa_CLASS) {

	      dy[ppw->pv->index_pt_shear_ur] = 
		-3./tau*y[ppw->pv->index_pt_shear_ur]
		+2./3.*y[ppw->pv->index_pt_theta_ur]
		+1./3.*pvecmetric[ppw->index_mt_h_prime];
	      
	    }
	  }
	}
      }
    }

    /** (h) non-cold dark matter (massive neutrinos, WDM, etc.) */

    if (pba->has_ncdm == _TRUE_) {
      idx = ppw->pv->index_pt_psi0_ncdm1;
      if(ppw->approx[ppw->index_ap_ncdmfa] == (int)ncdmfa_on){
	// Use fluid equations:
	for (n_ncdm=0; n_ncdm<ppw->pv->N_ncdm; n_ncdm++){
	  rho_ncdm_bg = ppw->pvecback[pba->index_bg_rho_ncdm1+n_ncdm];
	  p_ncdm_bg = ppw->pvecback[pba->index_bg_p_ncdm1+n_ncdm];
	  pseudo_p_ncdm = ppw->pvecback[pba->index_bg_pseudo_p_ncdm1+n_ncdm];

	  /* w is p/rho */
	  w_ncdm = p_ncdm_bg/rho_ncdm_bg;

	  /* c_a is the adiabatic sound speed */
	  ca2_ncdm = w_ncdm/3.0/(1.0+w_ncdm)*(5.0-pseudo_p_ncdm/p_ncdm_bg);

	  /* c_eff is (delta p / delta rho) in the gauge under
	     consideration (not in the gauge comoving with the
	     fluid) */
	  
	  /* c_vis is introduced in order to close the system */
	  
	  if (ppr->gauge == newtonian) { 
	  
	  }
	  if (ppr->gauge == synchronous) {
	    
	    /* different ansatz for sound speed c)_eff and viscosity speed c_vis */
	    if (ppr->ncdm_fluid_approximation == ncdmfa_mb) {
	       ceff2_ncdm = ca2_ncdm;
	       cvis2_ncdm = 3.*w_ncdm*ca2_ncdm;
	    }
	    if (ppr->ncdm_fluid_approximation == ncdmfa_hu) {
	       ceff2_ncdm = ca2_ncdm;
	       cvis2_ncdm = w_ncdm;
	    }
	    if (ppr->ncdm_fluid_approximation == ncdmfa_CLASS) {
	       ceff2_ncdm = ca2_ncdm;
	       cvis2_ncdm = 3.*w_ncdm*ca2_ncdm;
	    }

	    /* exact continuity equation */
	    dy[idx] = -(1.0+w_ncdm)*(y[idx+1]+0.5*pvecmetric[ppw->index_mt_h_prime])-
	      3.0*a_prime_over_a*(ceff2_ncdm-w_ncdm)*y[idx];
	    
	    /* exact euler equation */
	    dy[idx+1] = -a_prime_over_a*(1.0-3.0*ca2_ncdm)*y[idx+1]+
	      ceff2_ncdm/(1.0+w_ncdm)*k2*y[idx]-k2*y[idx+2];
	 
	    /* different ansatz for approximate shear derivative */
	    if (ppr->ncdm_fluid_approximation == ncdmfa_mb) {

	      dy[idx+2] = -3.0*(a_prime_over_a*(2./3.-ca2_ncdm-pseudo_p_ncdm/p_ncdm_bg/3.)+1./tau)*y[idx+2]
		+4.0/3.0*cvis2_ncdm/(1.0+w_ncdm)*(2.0*y[idx+1]+pvecmetric[ppw->index_mt_h_prime]+6.0*pvecmetric[ppw->index_mt_eta_prime]);
	    
	    }
	
	    if (ppr->ncdm_fluid_approximation == ncdmfa_hu) {

	      dy[idx+2] = -3.0*a_prime_over_a*ca2_ncdm/w_ncdm*y[idx+2]
		+4.0/3.0*cvis2_ncdm/(1.0+w_ncdm)*(2.0*y[idx+1]+pvecmetric[ppw->index_mt_h_prime]+6.0*pvecmetric[ppw->index_mt_eta_prime]);
	    
	    }
	    
	    if (ppr->ncdm_fluid_approximation == ncdmfa_CLASS) {

	      dy[idx+2] = -3.0*(a_prime_over_a*(2./3.-ca2_ncdm-pseudo_p_ncdm/p_ncdm_bg/3.)+1./tau)*y[idx+2]
		+4.0/3.0*cvis2_ncdm/(1.0+w_ncdm)*(2.0*y[idx+1]+pvecmetric[ppw->index_mt_h_prime]);
	    
	    }
	  }
	  // Jump to next species:
	  idx += ppw->pv->l_max_ncdm[n_ncdm]+1;
	}
      }
      else{
	// Use equations for the Boltzmann hierarchy on momentum grid:
	for (n_ncdm=0; n_ncdm<ppw->pv->N_ncdm; n_ncdm++){
	  for (index_q=0; index_q < ppw->pv->q_size_ncdm[n_ncdm]; index_q++){
	    dlnf0_dlnq = pba->dlnf0_dlnq_ncdm[n_ncdm][index_q];
	    q = pba->q_ncdm[n_ncdm][index_q];
	    epsilon = sqrt(q*q+a2*pba->M_ncdm[n_ncdm]*pba->M_ncdm[n_ncdm]);
	    qk_div_epsilon = k*q/epsilon;
			
	    if (ppr->gauge == synchronous) {
	      dy[idx] = -qk_div_epsilon*y[idx+1]+pvecmetric[ppw->index_mt_h_prime]*dlnf0_dlnq/6.;
					
	      dy[idx+1] = qk_div_epsilon/3.*(y[idx]-2.*y[idx+2]);
					
	      dy[idx+2] = qk_div_epsilon/5.0*(2*y[idx+1]-3.*y[idx+3])-h_plus_six_eta_prime/15.*dlnf0_dlnq;
	    }

	    if (ppr->gauge == newtonian){
	      dy[idx] = -qk_div_epsilon*y[idx+1]-pvecmetric[ppw->index_mt_phi_prime]*dlnf0_dlnq;
				
	      dy[idx+1] = qk_div_epsilon/3.0*(y[idx] - 2*y[idx+2])-
		epsilon*k/(3*q)*pvecmetric[ppw->index_mt_psi]*dlnf0_dlnq;
				
	      dy[idx+2] = qk_div_epsilon/5.0*(2*y[idx+1]-3*y[idx+3]);
				
	    }
						
	    for(l=3; l<ppw->pv->l_max_ncdm[n_ncdm]; l++){
	      dy[idx+l] = qk_div_epsilon/(2.*l+1.0)*(l*y[idx+(l-1)]-(l+1.)*y[idx+(l+1)]);
	    }
			
	    // Truncation as in Ma and Bertschinger. We have l = lmax-1;
	    //l is now l_max
	    dy[idx+l] = qk_div_epsilon*y[idx+l-1]-(1.+l)/tau*y[idx+l]; 
			
	    // Jump to next momentum bin:
	    idx += (ppw->pv->l_max_ncdm[n_ncdm]+1);
	  }
	}
      }
    }

    /** (j) metric */

    if (ppr->gauge == synchronous) {
      /* Synchronous gauge */
      dy[ppw->pv->index_pt_eta] = pvecmetric[ppw->index_mt_eta_prime];
    }


    /* for testing, will be useful for improving tight-coupling approximation */
    /*       if ((index_k == 0) && (tau > 800)) */
    /* 	printf("%e %e %e %e %e %e %e %e %e %e %e\n", */
    /* 	       tau, */
    /* 	       delta_g, */
    /* 	       theta_g, */
    /* 	       shear_g, */
    /* 	       y[ppw->pv->index_pt_l3_g], */
    /* 	       y[ppw->pv->index_pt_pol0_g], */
    /* 	       y[ppw->pv->index_pt_pol1_g], */
    /* 	       y[ppw->pv->index_pt_pol2_g], */
    /* 	       y[ppw->pv->index_pt_pol3_g], */
    /* 	       delta_b, */
    /* 	       theta_b); */

  }

  /** - tensor mode */

  if ((ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors)) {
      
    /* short-cut notations for the tensor perturbations */
    delta_g = y[ppw->pv->index_pt_delta_g];
    theta_g = y[ppw->pv->index_pt_theta_g];
    shear_g = y[ppw->pv->index_pt_shear_g];
    
    Psi = 
      delta_g/40.
      +2.*shear_g/35.
      +y[ppw->pv->index_pt_delta_g+4]/210.
      -3.*y[ppw->pv->index_pt_pol0_g]/5. 
      +6.*y[ppw->pv->index_pt_pol2_g]/35.
      -y[ppw->pv->index_pt_pol0_g+4]/210.;

    /* photon density (4*F_0) */
    dy[ppw->pv->index_pt_delta_g] = 
      -4./3.*theta_g
      -4.*y[ppw->pv->index_pt_gwdot]
      -pvecthermo[pth->index_th_dkappa]*(delta_g-4.*Psi);

    /* photon velocity ((3k/4)*F_1) */
    dy[ppw->pv->index_pt_theta_g] = 
      k2*(delta_g/4.-shear_g)
      -pvecthermo[pth->index_th_dkappa]*theta_g;

    /* photon shear (0.5*F_2) */
    dy[ppw->pv->index_pt_shear_g] =	
      0.5*(8./15.*theta_g
	   -3./5.*k*y[ppw->pv->index_pt_shear_g+1])
      -pvecthermo[pth->index_th_dkappa]*shear_g;

    /* photon l=3 */
    dy[ppw->pv->index_pt_l3_g] = 
      k/7.*(6.*shear_g-4.*y[ppw->pv->index_pt_l3_g+1])
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
      -(1.+l)/tau*y[ppw->pv->index_pt_delta_g+l]
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
      -(l+1.)/tau*y[ppw->pv->index_pt_pol0_g+l]
      -pvecthermo[pth->index_th_dkappa]*y[ppw->pv->index_pt_pol0_g+l];

    /* tensor metric perturbation h (gravitational waves) */
    dy[ppw->pv->index_pt_gw] = y[ppw->pv->index_pt_gwdot];     

    /* its time-derivative */
    dy[ppw->pv->index_pt_gwdot] = -2.*a_prime_over_a*y[ppw->pv->index_pt_gwdot]-k2*y[ppw->pv->index_pt_gw];

    /*     fprintf(stderr, */
    /* 	    "%g %g %g %g %g %g %g\n", */
    /* 	    y[ppw->pv->index_pt_gw], */
    /* 	    y[ppw->pv->index_pt_gwdot], */
    /* 	    delta_g, */
    /* 	    theta_g, */
    /* 	    shear_g, */
    /* 	    y[ppw->pv->index_pt_l3_g], */
    /* 	    y[ppw->pv->index_pt_l3_g+1]); */

    /*     class_test(0==1,error_message,"stop here\n"); */
	    
  }

  /*     printf("Leaves derivs with:\n"); */
  /*     printf("gamma : %e %e %e %e %e %e \n",dy[0],dy[1],dy[2],dy[3],dy[4],dy[5]); */
  /*     printf("b     : %e %e \n",dy[6],dy[7]); */
  /*     printf("cdm   : %e \n",dy[8]); */
  /*     printf("fluid : %e %e \n",dy[9],dy[10]); */
  /*     printf("nu    : %e %e %e %e %e %e \n",dy[10],dy[11],dy[12],dy[13],dy[14],dy[15]); */
  /*     printf("eta   : %e \n",dy[16]); */
  /*     printf("h     : %e \n",pvecmetric[ppw->index_mt_h_prime]); */

  return _SUCCESS_;
}



