/** @file transfer.c Documented transfer module.
 *
 * Julien Lesgourgues, 28.07.2013    
 *
 * This module has two purposes:
 *
 * - at the beginning, to compute the transfer functions \f$
 *   \Delta_l^{X} (q) \f$, and store them in tables used for
 *   interpolation in other modules.
 *
 * - at any time in the code, to evaluate the transfer functions (for
 *   a given mode, initial condition, type and multipole l) at any
 *   wavenumber q (by interpolating within the interpolation table).
 *
 * Hence the following functions can be called from other modules:
 *
 * -# transfer_init() at the beginning (but after perturb_init() 
 *    and bessel_init())
 *
 * -# transfer_functions_at_q() at any later time 
 *
 * -# transfer_free() at the end, when no more calls to 
 *    transfer_functions_at_q() are needed
 * 
 * Note that in the standard implementation of CLASS, only the pre-computed 
 * values of the transfer functions are used, no interpolation is necessary; 
 * hence the routine transfer_functions_at_q() is actually never called.
 */

#include "transfer.h"

/** 
 * Transfer function \f$ \Delta_l^{X} (q) \f$ at a given wavenumber q.
 *
 * For a given mode (scalar, vector, tensor), initial condition, type
 * (temperature, polarization, lensing, etc) and multipole, computes
 * the transfer function for an arbitary value of q by interpolating
 * between pre-computed values of q. This
 * function can be called from whatever module at whatever time,
 * provided that transfer_init() has been called before, and
 * transfer_free() has not been called yet.
 *
 * Wavenumbers are called q in this module and k in the perturbation
 * module. In flat universes k=q. In non-flat universes q and k differ
 * through q2 = k2 + K(1+m), where m=0,1,2 for scalar, vector,
 * tensor. q should be used throughout the transfer module, excepted
 * when interpolating or manipulating the source functions S(k,tau)
 * calculated in the perturbation module: for a given value of q, this
 * should be done at the corresponding k(q).
 *
 * @param index_md Input: index of requested mode
 * @param index_ic   Input: index of requested initial condition
 * @param index_tt   Input: index of requested type
 * @param index_l    Input: index of requested multipole
 * @param k          Input: any wavenumber
 * @param transfer_function Output: transfer function
 * @return the error status
 */

int transfer_functions_at_q(
                            struct transfers * ptr,
                            int index_md,
                            int index_ic,
                            int index_tt,
                            int index_l,
                            double q,
                            double * transfer_function
                            ) {
  /** Summary: */

  /** - interpolate in pre-computed table using array_interpolate_two() */
  class_call(array_interpolate_two(
                                   ptr->q,
                                   1,
                                   0,
                                   ptr->transfer[index_md]
                                   +((index_ic * ptr->tt_size[index_md] + index_tt) * ptr->l_size[index_md] + index_l)
                                   * ptr->q_size,
                                   1,
                                   ptr->q_size,
                                   q,
                                   transfer_function,
                                   1,
                                   ptr->error_message),
             ptr->error_message,
             ptr->error_message);
  
  return _SUCCESS_;
}

/**
 * This routine initializes the transfers structure, (in particular,
 * computes table of transfer functions \f$ \Delta_l^{X} (q) \f$)
 *
 * Main steps: 
 *
 * - initialize all indices in the transfers structure
 *   and allocate all its arrays using transfer_indices_of_transfers().
 *
 * - for each thread (in case of parallel run), initialize the fields of a memory zone called the transfer_workspace with transfer_workspace_init()
 *
 * - loop over q values. For each q, compute the bessel functions if needed with transfer_update_HIS(), and defer the calculation of all transfer functions to transfer_compute_for_each_q()
 * - for each thread, free the the workspace with transfer_workspace_free()
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

  /* running index for wavenumbers */
  int index_q;

  /* conformal time today */
  double tau0;
  /* conformal time at recombination */
  double tau_rec;

  /* maximum number of sampling times for transfer sources */
  int tau_size_max;

  /* array of source derivatives S''(k,tau) 
     (second derivative with respect to k, not tau!), 
     used to interpolate sources at the right values of k,
     sources_spline[index_md][index_ic * ppt->tp_size[index_md] + index_tp][index_tau * ppt->k_size + index_k]
  */
  double *** sources_spline;

  /* pointer on workspace (one per thread if openmp) */
  struct transfer_workspace * ptw;

  /** - array with the correspondance between the index of sources in
      the perturbation module and in the transfer module,
      tp_of_tt[index_md][index_tt]
  */
  int ** tp_of_tt;

  /* This code can be optionally compiled with the openmp option for parallel computation.
     Inside parallel regions, the use of the command "return" is forbidden.
     For error management, instead of "return _FAILURE_", we will set the variable below
     to "abort = _TRUE_". This will lead to a "return _FAILURE_" jus after leaving the 
     parallel region. */
  int abort;

#ifdef _OPENMP

  /* number of a given thread */
  int thread;

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
    fprintf(stdout,"Computing transfers\n");

  /** get number of modes (scalars, tensors...) */

  ptr->md_size = ppt->md_size;

  /** - get conformal age / recombination time 
      from background / thermodynamics structures 
      (only place where these structures are used in this module) */

  tau0 = pba->conformal_age;
  tau_rec = pth->tau_rec;

  /** - initialize all indices in the transfers structure and 
      allocate all its arrays using transfer_indices_of_transfers() */

  class_call(transfer_indices_of_transfers(ppr,ppt,pbs,ptr,tau0,pba->K,pba->sgnK),
             ptr->error_message,
             ptr->error_message);

  /** - spline all the sources passed by the perturbation module with respect to k (in order to interpolate later at a given value of k) */

  class_alloc(sources_spline,
              ptr->md_size*sizeof(double*),
              ptr->error_message);

  class_call(transfer_perturbation_source_spline(ppt,ptr,sources_spline),
             ptr->error_message,
             ptr->error_message);

  /** - allocate and fill array describing the correspondence between perturbation types and transfer types */

  class_alloc(tp_of_tt,
              ptr->md_size*sizeof(int*),
              ptr->error_message);

  class_call(transfer_get_source_correspondence(ppt,ptr,tp_of_tt),
             ptr->error_message,
             ptr->error_message);

  /** - evaluate maximum number of sampled times in the transfer
      sources: needs to be known here, in order to allocate a large
      enough workspace */

  class_call(transfer_source_tau_size_max(ppr,pba,ppt,ptr,tau_rec,tau0,&tau_size_max),
             ptr->error_message,
             ptr->error_message);

  /** (a.3.) workspace, allocated in a parallel zone since in openmp
      version there is one workspace per thread */
  
  /* initialize error management flag */
  abort = _FALSE_;
  
  /* beginning of parallel region */
  
#pragma omp parallel                                                    \
  shared(tau_size_max,pbs,ptr,ppr,pba,ppt,tp_of_tt,tau_rec,sources_spline,abort) \
  private(ptw,thread,index_q,tstart,tstop,tspent)
  {
    
#ifdef _OPENMP
    thread = omp_get_thread_num();
    tspent = 0.;
#endif
       
    /* allocate workspace */

    class_call_parallel(transfer_workspace_init(ptr,
                                                &ptw,
                                                ppt->tau_size,
                                                tau_size_max,
                                                pbs->use_pbs,
                                                pba->K,
                                                pba->sgnK,
                                                pbs->get_HIS_from_shared_memory),
                        ptr->error_message,
                        ptr->error_message);
    
    /** - loop over all wavenumbers (parallelised). For each wavenumber: */
    
#pragma omp for schedule (dynamic)
    
    for (index_q = 0; index_q < ptr->q_size; index_q++) {
      
#ifdef _OPENMP
      tstart = omp_get_wtime();
#endif
      
      /* Update interpolation structure: */
      class_call_parallel(transfer_update_HIS(ppr,
                                              pbs,
                                              ptr,
                                              ptw,
                                              index_q,
                                              tau0),
                          ptr->error_message,
                          ptr->error_message);

      class_call_parallel(transfer_compute_for_each_q(ppr,
                                                      pba,
                                                      ppt,
                                                      ptr,
                                                      tp_of_tt,
                                                      index_q,
                                                      tau_size_max,
                                                      tau_rec,
                                                      sources_spline,
                                                      ptw),
                          ptr->error_message,
                          ptr->error_message);
      
#ifdef _OPENMP
      tstop = omp_get_wtime();
      
      tspent += tstop-tstart;
#endif

#pragma omp flush(abort)
      
    } /* end of loop over wavenumber */
    
    /* free workspace allocated inside parallel zone */
    transfer_workspace_free(ptr,ptw);

#ifdef _OPENMP
    if (ptr->transfer_verbose>1)
      printf("In %s: time spent in parallel region (loop over k's) = %e s for thread %d\n",
             __func__,tspent,omp_get_thread_num());
#endif
    
  } /* end of parallel region */

  if (abort == _TRUE_) return _FAILURE_;
  
  /* finally, free arrays allocated outside parallel zone */

  class_call(transfer_perturbation_source_spline_free(ppt,ptr,sources_spline),
             ptr->error_message,
             ptr->error_message);

  class_call(transfer_free_source_correspondence(ptr,tp_of_tt),
             ptr->error_message,
             ptr->error_message);

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
      free(ptr->transfer[index_md]);
      free(ptr->k[index_md]);
    }  
   
    free(ptr->tt_size);
    free(ptr->l_size_tt);
    free(ptr->l_size);
    free(ptr->l);
    free(ptr->q);
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
                                  double tau0,
                                  double K,
                                  int sgnK
                                  ) {

  /** Summary: */

  /** - define local variables */

  int index_md,index_tt,index_tt_common;

  /** define indices for transfer types */

  class_alloc(ptr->tt_size,ptr->md_size * sizeof(int),ptr->error_message);

  /** - type indices common to scalars and tensors */

  index_tt = 0;

  if (ppt->has_cl_cmb_temperature == _TRUE_) {
    ptr->index_tt_t2 = index_tt;
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

    if (ppt->has_cl_cmb_temperature == _TRUE_) {
      ptr->index_tt_t0 = index_tt;
      index_tt++;
      ptr->index_tt_t1 = index_tt;
      index_tt++;
    }

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

  /** - type indices for vectors */

  if (ppt->has_vectors == _TRUE_) {
  
    index_tt = index_tt_common;

    if (ppt->has_cl_cmb_temperature == _TRUE_) {
      ptr->index_tt_t1 = index_tt;
      index_tt++;
    }

    if (ppt->has_cl_cmb_polarization == _TRUE_) {
      ptr->index_tt_b = index_tt;
      index_tt++;
    }

    ptr->tt_size[ppt->index_md_vectors]=index_tt;

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

  /* array (of array) of transfer functions for each mode, transfer[index_md] */

  class_alloc(ptr->transfer,ptr->md_size * sizeof(double *),ptr->error_message);

  /** get q values using transfer_get_q_list() */
  class_call(transfer_get_q_list(ppr,ppt,ptr,tau0,K,sgnK),
             ptr->error_message,
             ptr->error_message);
  if (pbs->get_HIS_from_shared_memory==_TRUE_){
    /** Read q_list from generated file. We should add some tests here. */
    free(ptr->q);
    FILE * q_file;
    int lSize;
    q_file = fopen("q_list_binary.bin","rb");
    // obtain file size:
    fseek (q_file , 0 , SEEK_END);
    lSize = ftell (q_file);
    ptr->q_size = lSize/sizeof(double);
    rewind (q_file);
    class_alloc(ptr->q,sizeof(double)*ptr->q_size,ptr->error_message);
    fread (ptr->q,sizeof(double),ptr->q_size,q_file);
    fclose(q_file);
  }

  /** get k values using transfer_get_k_list() */
  class_call(transfer_get_k_list(ppt,ptr,K),
             ptr->error_message,
             ptr->error_message);

  /** get l values using transfer_get_l_list() */
  class_call(transfer_get_l_list(ppr,ppt,pbs,ptr),
             ptr->error_message,
             ptr->error_message);

  /** - loop over modes (scalar, etc). For each mode: */
  
  for (index_md = 0; index_md < ptr->md_size; index_md++) {

    /** allocate arrays of transfer functions, (ptr->transfer[index_md])[index_ic][index_tt][index_l][index_k] */
    class_alloc(ptr->transfer[index_md],
                ppt->ic_size[index_md] * ptr->tt_size[index_md] * ptr->l_size[index_md] * ptr->q_size * sizeof(double),
                ptr->error_message);
    
  }

  return _SUCCESS_;

}

int transfer_perturbation_source_spline(
                                        struct perturbs * ppt,
                                        struct transfers * ptr,
                                        double *** sources_spline
                                        ) {
  int index_md;
  int index_ic;
  int index_tp;

  for (index_md = 0; index_md < ptr->md_size; index_md++) {

    class_alloc(sources_spline[index_md],
                ppt->ic_size[index_md]*ppt->tp_size[index_md]*sizeof(double*),
                ptr->error_message);

    for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++) {
      
      for (index_tp = 0; index_tp < ppt->tp_size[index_md]; index_tp++) {
      
        class_alloc(sources_spline[index_md][index_ic * ppt->tp_size[index_md] + index_tp],
                    ppt->k_size*ppt->tau_size*sizeof(double),
                    ptr->error_message);

        class_call(array_spline_table_columns(ppt->k,
                                              ppt->k_size,
                                              ppt->sources[index_md][index_ic * ppt->tp_size[index_md] + index_tp],
                                              ppt->tau_size,
                                              sources_spline[index_md][index_ic * ppt->tp_size[index_md] + index_tp],
                                              _SPLINE_EST_DERIV_,
                                              ptr->error_message),
                   ptr->error_message,
                   ptr->error_message);
        
      }
    }
  }

  return _SUCCESS_;

}

int transfer_perturbation_source_spline_free(
                                             struct perturbs * ppt,
                                             struct transfers * ptr,
                                             double *** sources_spline
                                             ) {
  int index_md;
  int index_ic;
  int index_tp;

  for (index_md = 0; index_md < ptr->md_size; index_md++) {
    for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++) {
      for (index_tp = 0; index_tp < ppt->tp_size[index_md]; index_tp++) {
        free(sources_spline[index_md][index_ic * ppt->tp_size[index_md] + index_tp]);
      }
    }
    free(sources_spline[index_md]);
  }
  free(sources_spline);

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
  int increment,current_l;

  /* allocate and fill l array */

  if (pbs->use_pbs == _TRUE_) {
    
    /* take directly from Bessel module) */
    
    ptr->l_size_max = pbs->l_size;
    
    class_alloc(ptr->l,ptr->l_size_max*sizeof(int),ptr->error_message);
    
    for (index_l=0; index_l < ptr->l_size_max; index_l++)
      ptr->l[index_l]=pbs->l[index_l];
  }

  else {

    /** - start from l = 2 and increase with logarithmic step */

    index_l = 0;
    current_l = 2;
    increment = max((int)(current_l * (ppr->l_logstep-1.)),1);
    
    while (((current_l+increment) < pbs->l_max) && 
           (increment < ppr->l_linstep)) {
      
      index_l ++;
      current_l += increment;
      increment = max((int)(current_l * (ppr->l_logstep-1.)),1);

    }

    /** - when the logarithmic step becomes larger than some linear step, 
        stick to this linear step till l_max */

    increment = ppr->l_linstep;

    while ((current_l+increment) <= pbs->l_max) {

      index_l ++;
      current_l += increment;

    }

    /** - last value set to exactly l_max */

    if (current_l != pbs->l_max) {

      index_l ++;
      current_l = pbs->l_max;

    } 

    ptr->l_size_max = index_l+1;

    /** - so far we just counted the number of values. Now repeat the
        whole thing but fill array with values. */

    class_alloc(ptr->l,ptr->l_size_max*sizeof(int),ptr->error_message);

    index_l = 0;
    ptr->l[0] = 2;
    increment = max((int)(ptr->l[0] * (ppr->l_logstep-1.)),1);

    while (((ptr->l[index_l]+increment) < pbs->l_max) && 
           (increment < ppr->l_linstep)) {
      
      index_l ++;
      ptr->l[index_l]=ptr->l[index_l-1]+increment;
      increment = max((int)(ptr->l[index_l] * (ppr->l_logstep-1.)),1);
 
    }

    increment = ppr->l_linstep;

    while ((ptr->l[index_l]+increment) <= pbs->l_max) {

      index_l ++;
      ptr->l[index_l]=ptr->l[index_l-1]+increment;
 
    }

    if (ptr->l[index_l] != pbs->l_max) {

      index_l ++;
      ptr->l[index_l]= pbs->l_max;
       
    }

  }

  /* for each mode and type, find relevant size of l array,
     l_size_tt[index_md][index_tt] (since for some modes and types
     l_max can be smaller). Also, maximize this size for each mode to
     find l_size[index_md]. */

  for (index_md=0; index_md < ppt->md_size; index_md++) {

    ptr->l_size[index_md] = 0;

    for (index_tt=0;index_tt<ptr->tt_size[index_md];index_tt++) {

      if _scalars_ {

          if ((ppt->has_cl_cmb_temperature == _TRUE_) && 
              ((index_tt == ptr->index_tt_t0) || (index_tt == ptr->index_tt_t1) || (index_tt == ptr->index_tt_t2)))
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
 * This routine defines the number and values of wavenumbers q for
 * each mode (different in perturbation module and transfer module:
 * here we impose an upper bound on the linear step. So, typically,
 * for small q, the sampling is identical to that in the perturbation
 * module, while at high q it is denser and source functions are
 * interpolated).
 *
 * @param ppr     Input : pointer to precision structure
 * @param ppt     Input : pointer to perturbation structure
 * @param ptr     Input/Output : pointer to transfers structure containing q's
 * @param rs_rec  Input : comoving distance to recombination
 * @param index_md Input: index of requested mode (scalar, tensor, etc) 
 * @return the error status
 */

int transfer_get_q_list(
                        struct precision * ppr,
                        struct perturbs * ppt,
                        struct transfers * ptr,
                        double tau0,
                        double K,
                        int sgnK
                        ) {

  int index_k;
  int index_q;
  double q_min=0.,q_max,q_step_max=0.,k_max;
  int nu, nu_min, nu_proposed, nu_proposed_following_ppt, nu_proposed_following_step_max;
  int q_size_max;

  /* find q_step_max, the maximum value of the step */

  q_step_max = 2.*_PI_/tau0*ppr->k_step_trans;

  class_test(q_step_max == 0.,
             ptr->error_message,
             "stop to avoid infinite loop");

  /* first deal with case K=0 (flat) and K<0 (open). The case K>0 (closed) is very different, and is dealt with separately below. */

  if (sgnK <= 0) {

    /* first and last value */

    if (sgnK == 0) {
      q_min = ppt->k[0];
      q_max = ppt->k[ppt->k_size_cl-1]; 
      K=0;
    }
    else {
      q_min = sqrt(ppt->k[0]*ppt->k[0]+K);
      k_max = ppt->k[ppt->k_size_cl-1];
      q_max = sqrt(k_max*k_max+K);
      if (ppt->has_vectors == _TRUE_) 
        q_max = min(q_max,sqrt(k_max*k_max+2.*K));
      if (ppt->has_tensors == _TRUE_) 
        q_max = min(q_max,sqrt(k_max*k_max+3.*K));
    }

    /* conservative estimate of maximum size of the list (will be reduced later with realloc) */

    q_size_max = 2+ppt->k_size_cl+(int)((q_max-q_min)/q_step_max);

    class_alloc(ptr->q,
                q_size_max*sizeof(double),
                ptr->error_message);

    /* - first point */
    
    index_q = 0;

    ptr->q[index_q] = q_min; 

    index_q++;

    /* - points taken from perturbation module if step small enough */

    while ((index_q < ppt->k_size_cl) && ((sqrt(ppt->k[index_q]*ppt->k[index_q]+K) - ptr->q[index_q-1]) < q_step_max)) {
      
      class_test(index_q >= q_size_max,ptr->error_message,"buggy q-list definition");
      ptr->q[index_q] = sqrt(ppt->k[index_q]*ppt->k[index_q]+K);
      index_q++;

    }

    /* - then, points spaced linearily with step q_step_max */

    while (ptr->q[index_q-1] < q_max) {
      
      class_test(index_q >= q_size_max,ptr->error_message,"buggy q-list definition");
      ptr->q[index_q] = ptr->q[index_q-1] + q_step_max;
      index_q++;
      
    }
    
    /* - get number of valid points in order to re-allocate list */
    
    if (ptr->q[index_q-1] > q_max)
      ptr->q_size=index_q-1;
    else
      ptr->q_size=index_q;
    
  }

  /* deal with K>0 (closed) case */

  else {

    /* first and last value */
    
    nu_min = 3;
    q_min = nu_min * sqrt(K);
    q_max = ppt->k[ppt->k_size_cl-1]; 
    
    /* conservative estimate of maximum size of the list (will be reduced later with realloc) */

    q_size_max = 2+(int)((q_max-q_min)/sqrt(K));

    class_alloc(ptr->q,
                q_size_max*sizeof(double),
                ptr->error_message);
    
    /* - first point */
    
    index_q = 0;
    index_k = 0;
    
    ptr->q[index_q] = q_min; 
    nu = nu_min;

    index_q++;

    while (ptr->q[index_q-1] < q_max) {

      /* find first point in ppt->k list corresponding to a q that is higher than the last q value */

      while ((index_k < ppt->k_size_cl) && (ppt->k[index_k]*ppt->k[index_k]+K < ptr->q[index_q-1]*ptr->q[index_q-1])) index_k++;

      /* now we can propose as a next point q_propose_following_ppt, which is bigger than the last q value */

      nu_proposed_following_ppt = (int)(sqrt(ppt->k[index_k]*ppt->k[index_k]+K)/sqrt(K));

      /* but if spacing in ppt become too sparse, we should instead refer to the maximum stepsize */

      nu_proposed_following_step_max = (int)((ptr->q[index_q-1]+q_step_max)/sqrt(K));

      nu_proposed = min(nu_proposed_following_ppt,nu_proposed_following_step_max);

      if (nu_proposed <= nu) nu_proposed=nu+1;
      
      ptr->q[index_q] = nu_proposed*sqrt(K);
      nu = nu_proposed;

      index_q++;
    }

    /* - get number of valid points in order to re-allocate list */

    if (ptr->q[index_q-1] > q_max)
      ptr->q_size=index_q-1;
    else
      ptr->q_size=index_q;
    
  }

  /* check size of q_list and realloc the array to the correct size */

  class_test(ptr->q_size<2,ptr->error_message,"buggy q-list definition");

  class_realloc(ptr->q,
                ptr->q,
                ptr->q_size*sizeof(double),
                ptr->error_message);

  /* consistency checks */

  class_test(ptr->q[0] <= 0.,
             ptr->error_message,
             "bug in q list calculation, q_min=%e, should always be strictly positive",ptr->q[0]);

  if (sgnK == 1) {
    class_test(ptr->q[0] < 3.*sqrt(K),
               ptr->error_message,
               "bug in q list calculation, q_min=%e, should be greater or equal to 3sqrt(K)=%e in positivevly curved universe",ptr->q[0],3.*sqrt(K));
  }    

  for (index_q=1; index_q<ptr->q_size; index_q++) {
    class_test(ptr->q[index_q] <= ptr->q[index_q-1],
               ptr->error_message,
               "bug in q list calculation, q values should be in strictly growing order"); 
  }

  FILE *fid=fopen("q_from_class.dat","w");
  int i;
  for (i=0; i<ptr->q_size; i++){
    fprintf(fid,"%.16e ",ptr->q[i]);
  }
  fclose(fid);


  return _SUCCESS_; 

}

/**
 * This routine defines the number and values of wavenumbers q for
 * each mode (different in perturbation module and transfer module:
 * here we impose an upper bound on the linear step. So, typically,
 * for small q, the sampling is identical to that in the perturbation
 * module, while at high q it is denser and source functions are
 * interpolated).
 *
 * @param ppr     Input : pointer to precision structure
 * @param ppt     Input : pointer to perturbation structure
 * @param ptr     Input/Output : pointer to transfers structure containing q's
 * @param rs_rec  Input : comoving distance to recombination
 * @param index_md Input: index of requested mode (scalar, tensor, etc) 
 * @return the error status
 */

int transfer_get_q_list2(
                        struct precision * ppr,
                        struct perturbs * ppt,
                        double tau0,
                        double K,
                        int sgnK,
                        double **qin,
                        int *q_size,
                        ErrorMsg error_message
                        ) {

  int index_k;
  int index_q;
  double q_min=0.,q_max,q_step_max=0.,k_max;
  int nu, nu_min, nu_proposed, nu_proposed_following_ppt, nu_proposed_following_step_max;
  int q_size_max;
  double * q;
  /* find q_step_max, the maximum value of the step */

  q_step_max = 2.*_PI_/tau0*ppr->k_step_trans;

  class_test(q_step_max == 0.,
             error_message,
             "stop to avoid infinite loop");

  /* first deal with case K=0 (flat) and K<0 (open). The case K>0 (closed) is very different, and is dealt with separately below. */

  if (sgnK <= 0) {

    /* first and last value */

    if (sgnK == 0) {
      q_min = ppt->k[0];
      q_max = ppt->k[ppt->k_size_cl-1]; 
      K=0;
    }
    else {
      q_min = sqrt(ppt->k[0]*ppt->k[0]+K);
      k_max = ppt->k[ppt->k_size_cl-1];
      q_max = sqrt(k_max*k_max+K);
      if (ppt->has_vectors == _TRUE_) 
        q_max = min(q_max,sqrt(k_max*k_max+2.*K));
      if (ppt->has_tensors == _TRUE_) 
        q_max = min(q_max,sqrt(k_max*k_max+3.*K));
    }

    /* conservative estimate of maximum size of the list (will be reduced later with realloc) */

    q_size_max = 2+ppt->k_size_cl+(int)((q_max-q_min)/q_step_max);

    class_alloc(q,
                q_size_max*sizeof(double),
                error_message);

    /* - first point */
    
    index_q = 0;

    q[index_q] = q_min; 

    index_q++;

    /* - points taken from perturbation module if step small enough */

    while ((index_q < ppt->k_size_cl) && ((sqrt(ppt->k[index_q]*ppt->k[index_q]+K) - q[index_q-1]) < q_step_max)) {
      
      class_test(index_q >= q_size_max,error_message,"buggy q-list definition");
      q[index_q] = sqrt(ppt->k[index_q]*ppt->k[index_q]+K);
      index_q++;

    }

    /* - then, points spaced linearily with step q_step_max */

    while (q[index_q-1] < q_max) {
      
      class_test(index_q >= q_size_max,error_message,"buggy q-list definition");
      q[index_q] = q[index_q-1] + q_step_max;
      index_q++;
      
    }
    
    /* - get number of valid points in order to re-allocate list */
    
    if (q[index_q-1] > q_max)
      *q_size=index_q-1;
    else
      *q_size=index_q;
    
  }

  /* deal with K>0 (closed) case */

  else {

    /* first and last value */
    
    nu_min = 3;
    q_min = nu_min * sqrt(K);
    q_max = ppt->k[ppt->k_size_cl-1]; 
    
    /* conservative estimate of maximum size of the list (will be reduced later with realloc) */

    q_size_max = 2+(int)((q_max-q_min)/sqrt(K));

    class_alloc(q,
                q_size_max*sizeof(double),
                error_message);
    
    /* - first point */
    
    index_q = 0;
    index_k = 0;
    
    q[index_q] = q_min; 
    nu = nu_min;

    index_q++;

    while (q[index_q-1] < q_max) {

      /* find first point in ppt->k list corresponding to a q that is higher than the last q value */

      while ((index_k < ppt->k_size_cl) && (ppt->k[index_k]*ppt->k[index_k]+K < q[index_q-1]*q[index_q-1])) index_k++;

      /* now we can propose as a next point q_propose_following_ppt, which is bigger than the last q value */

      nu_proposed_following_ppt = (int)(sqrt(ppt->k[index_k]*ppt->k[index_k]+K)/sqrt(K));

      /* but if spacing in ppt become too sparse, we should instead refer to the maximum stepsize */

      nu_proposed_following_step_max = (int)((q[index_q-1]+q_step_max)/sqrt(K));

      nu_proposed = min(nu_proposed_following_ppt,nu_proposed_following_step_max);

      if (nu_proposed <= nu) nu_proposed=nu+1;
      
      q[index_q] = nu_proposed*sqrt(K);
      nu = nu_proposed;

      index_q++;
    }

    /* - get number of valid points in order to re-allocate list */

    if (q[index_q-1] > q_max)
      *q_size=index_q-1;
    else
      *q_size=index_q;
    
  }

  /* check size of q_list and realloc the array to the correct size */

  class_test((*q_size)<2,error_message,"buggy q-list definition");

  class_realloc(q,
                q,
                (*q_size)*sizeof(double),
                error_message);
  /* consistency checks */

  class_test(q[0] <= 0.,
             error_message,
             "bug in q list calculation, q_min=%e, should always be strictly positive",q[0]);

  if (sgnK == 1) {
    class_test(q[0] < 3.*sqrt(K),
               error_message,
               "bug in q list calculation, q_min=%e, should be greater or equal to 3sqrt(K)=%e in positivevly curved universe",q[0],3.*sqrt(K));
  }    

  for (index_q=1; index_q<(*q_size); index_q++) {
    class_test(q[index_q] <= q[index_q-1],
               error_message,
               "bug in q list calculation, q values should be in strictly growing order"); 
  }

  *qin = q;
  return _SUCCESS_; 

}



int transfer_get_k_list(
                        struct perturbs * ppt,
                        struct transfers * ptr,
                        double K
                        ) {

  int index_md;
  int index_q;
  double m=0.;

  class_alloc(ptr->k,ptr->md_size*sizeof(double*),ptr->error_message);

  for (index_md = 0; index_md <  ptr->md_size; index_md++) {

    class_alloc(ptr->k[index_md],ptr->q_size*sizeof(double),ptr->error_message);

    if _scalars_ {
        m=0.;
      }
    if _vectors_ {
        m=1.;
      }
    if _tensors_ {
        m=2.;
      }

    for (index_q=0; index_q < ptr->q_size; index_q++) {
      ptr->k[index_md][index_q] = sqrt(ptr->q[index_q]*ptr->q[index_q]-K*(m+1.));
    }

    class_test(ptr->k[index_md][0] < ppt->k[0],
               ptr->error_message,
               "bug in k_list calculation: in perturbation module k_min=%e, in transfer module k_min[mode=%d]=%e, interpolation impossible",
               ppt->k[0],
               index_md,
               ptr->k[index_md][0]);

    class_test(ptr->k[index_md][ptr->q_size-1] > ppt->k[ppt->k_size_cl-1],
               ptr->error_message,
               "bug in k_list calculation: in perturbation module k_max=%e, in transfer module k_max[mode=%d]=%e, interpolation impossible",
               ppt->k[ppt->k_size_cl],
               index_md,
               ptr->k[index_md][ptr->q_size-1]);
      

  }

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
                                       int ** tp_of_tt
                                       ) {
  
  /* running index on modes */
  int index_md;

  /* running index on transfer types */
  int index_tt;

  /** - which source are we considering? Define correspondence
      between transfer types and source types */
  
  for (index_md = 0; index_md < ptr->md_size; index_md++) {

    class_alloc(tp_of_tt[index_md],ptr->tt_size[index_md]*sizeof(int),ptr->error_message);

    for (index_tt=0; index_tt<ptr->tt_size[index_md]; index_tt++) {
    
      if _scalars_ {
      
          if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t0)) 
            tp_of_tt[index_md][index_tt]=ppt->index_tp_t0;
          
          if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t1)) 
            tp_of_tt[index_md][index_tt]=ppt->index_tp_t1;
          
          if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t2)) 
            tp_of_tt[index_md][index_tt]=ppt->index_tp_t2;
          
          if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_e)) 
            tp_of_tt[index_md][index_tt]=ppt->index_tp_p;
          
          if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (index_tt == ptr->index_tt_lcmb)) 
            tp_of_tt[index_md][index_tt]=ppt->index_tp_g;
          
          if ((ppt->has_cl_density == _TRUE_) && (index_tt >= ptr->index_tt_density) && (index_tt < ptr->index_tt_density+ppt->selection_num))
            tp_of_tt[index_md][index_tt]=ppt->index_tp_g;
          
          if ((ppt->has_cl_lensing_potential == _TRUE_) && (index_tt >= ptr->index_tt_lensing) && (index_tt < ptr->index_tt_lensing+ppt->selection_num))
            tp_of_tt[index_md][index_tt]=ppt->index_tp_g;
          
        }
    
      if _vectors_ {
          
          if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t1)) 
            tp_of_tt[index_md][index_tt]=ppt->index_tp_t1;
          
          if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t2)) 
            tp_of_tt[index_md][index_tt]=ppt->index_tp_t2;
          
          if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_e)) 
            tp_of_tt[index_md][index_tt]=ppt->index_tp_p;
          
          if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_b))
            tp_of_tt[index_md][index_tt]=ppt->index_tp_p;
        }
      
      if _tensors_ {
          
          if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t2)) 
            tp_of_tt[index_md][index_tt]=ppt->index_tp_t2;
          
          if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_e)) 
            tp_of_tt[index_md][index_tt]=ppt->index_tp_p;
          
          if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_b))
            tp_of_tt[index_md][index_tt]=ppt->index_tp_p;
        }
    }
  }
  
  return _SUCCESS_;
  
}

int transfer_free_source_correspondence(
                                        struct transfers * ptr,
                                        int ** tp_of_tt
                                        ) {

  int index_md;

  for (index_md = 0; index_md < ptr->md_size; index_md++) {
    free(tp_of_tt[index_md]);
  }
  free(tp_of_tt);
  
  return _SUCCESS_;

}

int transfer_source_tau_size_max(
                                 struct precision * ppr,
                                 struct background * pba,
                                 struct perturbs * ppt,
                                 struct transfers * ptr,
                                 double tau_rec,
                                 double tau0,
                                 int * tau_size_max
                                 ) {

  int index_md;
  int index_tt;
  int tau_size_tt;

  *tau_size_max = 0;

  for (index_md = 0; index_md < ptr->md_size; index_md++) {
    
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
       
      *tau_size_max = max(*tau_size_max,tau_size_tt);
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
      if ((ppt->has_cl_cmb_temperature == _TRUE_) && 
          ((index_tt == ptr->index_tt_t0) || (index_tt == ptr->index_tt_t1) || (index_tt == ptr->index_tt_t2)))
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

int transfer_compute_for_each_q(
                                struct precision * ppr,
                                struct background * pba,
                                struct perturbs * ppt,
                                struct transfers * ptr,
                                int ** tp_of_tt,
                                int index_q,
                                int tau_size_max,
                                double tau_rec,
                                double *** sources_spline,
                                struct transfer_workspace * ptw
                                ) {
  
  /** Summary: */

  /** - define local variables */
  
  /* running index for modes */
  int index_md; 
  /* running index for initial conditions */
  int index_ic; 
  /* running index for transfer types */
  int index_tt; 
  /* running index for multipoles */
  int index_l; 

  /* we deal with workspaces, i.e. with contiguous memory zones (one
     per thread) containing various fields used by the integration
     routine */

  /* - first workspace field: perturbation source interpolated from perturbation stucture */
  double * interpolated_sources;

  /* - second workspace field: list of tau0-tau values, tau0_minus_tau[index_tau] */
  double * tau0_minus_tau;

  /* - third workspace field: list of trapezoidal weights for integration in tau0_minus_tau */
  double * w_trapz;

  /* - fourth workspace field, containing just a double: number of time values */
  int * tau_size;

  /* - fifth workspace field, identical to above interpolated sources:
     sources[index_tau] */
  double * sources;

  /** - for a given l, maximum value of k such that we can convolve
      the source with Bessel functions j_l(x) without reaching x_max */
  double q_max_bessel;

  /* a value of index_type */
  int previous_type;

  double l;

  short neglect;

  /** store the sources in the workspace and define all
      fields in this workspace */
  interpolated_sources = ptw->interpolated_sources;
  tau0_minus_tau = ptw->tau0_minus_tau;
  w_trapz  = ptw->w_trapz;
  tau_size = &(ptw->tau_size);
  sources = ptw->sources;

  /** - loop over all modes. For each mode: */ 
    
  for (index_md = 0; index_md < ptr->md_size; index_md++) {
    
    /** - loop over initial conditions. For each of them: */
    
    for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++) {
      
      /* initialize the previous type index */
      previous_type=-1;
      
      /** - loop over types. For each of them: */

      for (index_tt = 0; index_tt < ptr->tt_size[index_md]; index_tt++) {
        
        /** check if we must now deal with a new source with a
            new index ppt->index_type. If yes, interpolate it at the
            right values of k. */
        
        if (tp_of_tt[index_md][index_tt] != previous_type) {
          
          class_call(transfer_interpolate_sources(ppt,
                                                  ptr,
                                                  index_q,
                                                  index_md,
                                                  index_ic,
                                                  tp_of_tt[index_md][index_tt],
                                                  sources_spline[index_md][index_ic * ppt->tp_size[index_md] + tp_of_tt[index_md][index_tt]],
                                                  interpolated_sources),
                     ptr->error_message,
                     ptr->error_message);
        }
        
        previous_type = tp_of_tt[index_md][index_tt];
        
        /* the code makes a distinction between "perturbation
           sources" (e.g. gravitational potential) and "transfer
           sources" (e.g. total density fluctuations, obtained
           through the Poisson equation, and observed with a given
           selection function).
           
           The next routine computes the transfer source given the
           interpolated perturbation source, and copies it in the
           workspace. */
        
        class_call(transfer_sources(ppr,
                                    pba,
                                    ppt,
                                    ptr,
                                    interpolated_sources,
                                    tau_rec,
                                    index_q,
                                    index_md,
                                    index_tt,
                                    sources,
                                    tau0_minus_tau,
                                    w_trapz,
                                    tau_size),
                   ptr->error_message,
                   ptr->error_message);

        /* now that the array of times tau0_minus_tau is known, we can
           infer the arry of radial coordinates r(tau0_minus_tau) as well as a
           few other quantities related by trigonometric functions */

        class_call(transfer_radial_coordinates(ptr,ptw,index_md,index_q),
                   ptr->error_message,
                   ptr->error_message);

        for (index_l = 0; index_l < ptr->l_size[index_md]; index_l++) {
          
          l = (double)ptr->l[index_l];

          /* neglect transfer function when l is much smaller than k*tau0 */
          class_call(transfer_can_be_neglected(ppr,
                                               pba,
                                               ppt,
                                               ptr,
                                               index_md,
                                               index_ic,
                                               index_tt,
                                               ptr->q[index_q],
                                               l,
                                               &neglect),
                     ptr->error_message,
                     ptr->error_message);

          /* for K>0 (closed), transfer functions only defined for l<nu */
          if ((ptw->sgnK == 1) && (ptr->l[index_l] >= (int)(ptr->q[index_q]/sqrt(ptw->K)+0.2))) {
            neglect = _TRUE_;
          }
          /* This would maybe go into transfer_can_be_neglected later: */
          if ((ptw->sgnK != 0) && (index_l>=ptw->HIS.l_size)){
            neglect = _TRUE_;
          }
          if (neglect == _TRUE_) {

            ptr->transfer[index_md][((index_ic * ptr->tt_size[index_md] + index_tt)
                                     * ptr->l_size[index_md] + index_l)
                                    * ptr->q_size + index_q] = 0.;
          }
          else {

            /* for a given l, maximum value of k such that we can
               convolve the source with Bessel functions j_l(x)
               without reaching x_max (this is relevant in the flat
               case when the bessels are compiuted with the old bessel
               module. otherwise this condition is guaranteed by the
               choice of proper xmax when computing bessels) */
            if (ptw->sgnK == 0) {
              q_max_bessel = ptw->HIS.x[ptw->HIS.x_size-1]/tau0_minus_tau[0];
            }
            else {
              q_max_bessel = ptr->q[ptr->q_size-1];
            }

            /* compute the transfer function for this l */
            class_call(transfer_compute_for_each_l(
                                                   ptw,
                                                   ppr,
                                                   ppt,
                                                   ptr,
                                                   index_q,
                                                   index_md,
                                                   index_ic,
                                                   index_tt,
                                                   index_l,
                                                   l,
                                                   q_max_bessel
                                                   ),
                       ptr->error_message,
                       ptr->error_message);
          }

        } /* end of loop over l */

      } /* end of loop over type */
        
    } /* end of loop over initial condition */

  } /* end of loop over mode */
    
  return _SUCCESS_;

}

int transfer_radial_coordinates(
                                struct transfers * ptr,
                                struct transfer_workspace * ptw,
                                int index_md,
                                int index_q
                                ) {

  int index_tau;
  double sqrt_absK=0.;

  switch (ptw->sgnK){
  case 1:
    sqrt_absK = sqrt(ptw->K);
    for (index_tau=0; index_tau < ptw->tau_size; index_tau++) {
      ptw->chi[index_tau] = sqrt_absK*ptw->tau0_minus_tau[index_tau];
      ptw->cscKgen[index_tau] = sqrt_absK/ptr->k[index_md][index_q]/sin(ptw->chi[index_tau]);
      ptw->cotKgen[index_tau] = ptw->cscKgen[index_tau]*cos(ptw->chi[index_tau]);
    }
    break;
  case 0:
    for (index_tau=0; index_tau < ptw->tau_size; index_tau++) {
      ptw->chi[index_tau] = ptr->k[index_md][index_q] * ptw->tau0_minus_tau[index_tau];
      ptw->cscKgen[index_tau] = 1.0/ptw->chi[index_tau];
      ptw->cotKgen[index_tau] = 1.0/ptw->chi[index_tau];
    }
    break;
  case -1:
    sqrt_absK = sqrt(-ptw->K);
    for (index_tau=0; index_tau < ptw->tau_size; index_tau++) {
      ptw->chi[index_tau] = sqrt_absK*ptw->tau0_minus_tau[index_tau];
      ptw->cscKgen[index_tau] = sqrt_absK/ptr->k[index_md][index_q]/sinh(ptw->chi[index_tau]);
      ptw->cotKgen[index_tau] = ptw->cscKgen[index_tau]*cosh(ptw->chi[index_tau]);
    }
    break;
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
                                 int index_q,
                                 int index_md,
                                 int index_ic,
                                 int index_type,
                                 double * source_spline, /* array with argument source_spline[index_tau*ppt->k_size[index_md]+index_k] (must be allocated) */
                                 double * interpolated_sources /* array with argument interpolated_sources[index_q*ppt->tau_size+index_tau] (must be allocated) */
                                 ) {

  /** Summary: */

  /** - define local variables */

  /* index running on k values in the original source array */
  int index_k;

  /* index running on time */
  int index_tau;

  /* variables used for spline interpolation algorithm */
  double h, a, b;

  /** - find second derivative of original sources with respect to k
      in view of spline interpolation */
  /*
    class_call(array_spline_table_columns(ppt->k,
    ppt->k_size,
    ppt->sources[index_md][index_ic * ppt->tp_size[index_md] + index_type],
    ppt->tau_size,
    source_spline,
    _SPLINE_EST_DERIV_,
    ptr->error_message),
    ptr->error_message,
    ptr->error_message);
  */

  /** - interpolate at each k value using the usual 
      spline interpolation algorithm. */
  
  index_k = 0;
  h = ppt->k[index_k+1] - ppt->k[index_k];
    
  while (((index_k+1) < ppt->k_size) &&
         (ppt->k[index_k+1] < 
          ptr->k[index_md][index_q])) {
    index_k++;
    h = ppt->k[index_k+1] - ppt->k[index_k];
  }
    
  class_test(h==0.,
             ptr->error_message,
             "stop to avoid division by zero");
    
  b = (ptr->k[index_md][index_q] - ppt->k[index_k])/h;
  a = 1.-b;
    
  for (index_tau = 0; index_tau < ppt->tau_size; index_tau++) {
      
    interpolated_sources[index_tau] = 
      a * ppt->sources[index_md]
      [index_ic * ppt->tp_size[index_md] + index_type]
      [index_tau*ppt->k_size+index_k]
      + b * ppt->sources[index_md]
      [index_ic * ppt->tp_size[index_md] + index_type]
      [index_tau*ppt->k_size+index_k+1]
      + ((a*a*a-a) * source_spline[index_tau*ppt->k_size+index_k]
         +(b*b*b-b) * source_spline[index_tau*ppt->k_size+index_k+1])*h*h/6.0;
      
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
 * @param w_trapz               Output: trapezoidal weights for integration in (tau0-tau)
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
                     int index_q,
                     int index_md,
                     int index_tt,
                     double * sources,
                     double * tau0_minus_tau,
                     double * w_trapz,
                     int * tau_size_out
                     )  {

  /** Summary: */
  
  /** - define local variables */

  /* index running on time */
  int index_tau;

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
  
  /* trapezoidal weights for lensing source selection function */
  double * w_trapz_lensing_sources;

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
            sources[index_tau-index_tau_min] = 
              interpolated_sources[index_tau]
              * rescaling
              * ptr->lcmb_rescale
              * pow(ptr->k[index_md][index_q]/ptr->lcmb_pivot,ptr->lcmb_tilt);
          
            /* store value of (tau0-tau) */
            tau0_minus_tau[index_tau-index_tau_min] = tau0 - tau;
          
          }

          /* Compute trapezoidal weights for integration in (tau0-tau) */
          class_call(array_trapezoidal_weights(tau0_minus_tau,
                                               tau_size,
                                               w_trapz,
                                               ptr->error_message),
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

          /* Compute trapezoidal weights for integration in (tau0-tau) */
          class_call(array_trapezoidal_weights(tau0_minus_tau,
                                               tau_size,
                                               w_trapz,
                                               ptr->error_message),
                     ptr->error_message,
                     ptr->error_message);

          /* compute values of selection function at sampled values of tau */
          class_call(transfer_selection_compute(ppr,
                                                pba,
                                                ppt,
                                                ptr,
                                                selection,
                                                tau0_minus_tau,
                                                w_trapz,
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
          
            sources[index_tau] *= rescaling*pow(ptr->k[index_md][index_q],2);
          }

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
        
          class_alloc(w_trapz_lensing_sources,
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
        
          /* Compute trapezoidal weights for integration in (tau0-tau) */
          class_call(array_trapezoidal_weights(tau0_minus_tau_lensing_sources,
                                               tau_sources_size,
                                               w_trapz_lensing_sources,
                                               ptr->error_message),
                     ptr->error_message,
                     ptr->error_message);
        
          /* compute values of selection function at sampled values of tau */
          class_call(transfer_selection_compute(ppr,
                                                pba,
                                                ppt,
                                                ptr,
                                                selection,
                                                tau0_minus_tau_lensing_sources,
                                                w_trapz_lensing_sources,
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

          /* Compute trapezoidal weights for integration in (tau0-tau) */
          class_call(array_trapezoidal_weights(tau0_minus_tau,
                                               tau_size,
                                               w_trapz,
                                               ptr->error_message),
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
                    * w_trapz_lensing_sources[index_tau_sources];
                }
              }
              
              rescaling /= 2.;

            }
          
            /* copy from input array to output array */
            sources[index_tau] *= rescaling;
          
          }

          /* deallocate temporary arrays */
          free(pvecback);
          free(selection);
          free(tau0_minus_tau_lensing_sources);
          free(w_trapz_lensing_sources);
        
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
           ppt->tau_size*sizeof(double));
        
    /* store values of (tau0-tau) */
    for (index_tau=0; index_tau < ppt->tau_size; index_tau++) {
      tau0_minus_tau[index_tau] = tau0 - ppt->tau_sampling[index_tau];
    }

    /* Compute trapezoidal weights for integration in (tau0-tau) */
    class_call(array_trapezoidal_weights(tau0_minus_tau,
                                         tau_size,
                                         w_trapz,
                                         ptr->error_message),
               ptr->error_message,
               ptr->error_message);
  }
  
  /* return tau_size value that will be stored in the workspace (the
     workspace wants a double) */

  *tau_size_out = tau_size;

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

  /* array of values of source */
  double * source_at_tau;

  /* array of source values for a given time and for all k's */
  class_alloc(source_at_tau,
              sizeof(double),
              ptr->error_message);

  /* interpolate the sources linearily at the new time values */
  for (index_tau=0; index_tau<tau_size; index_tau++) {

    class_call(array_interpolate_two(ppt->tau_sampling,
                                     1,
                                     0,
                                     interpolated_sources, 
                                     1,
                                     ppt->tau_size,
                                     tau0-tau0_minus_tau[index_tau],
                                     source_at_tau,
                                     1,
                                     ptr->error_message),
               ptr->error_message,
               ptr->error_message);

    /* copy the new values in the output sources array */
    sources[index_tau] = source_at_tau[0];
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
 * @param w_trapz               Input : trapezoidal weights for (tau0-tau) integration
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
                               double * w_trapz,
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
  class_call(array_trapezoidal_integral(selection, 
                                        tau_size, 
                                        w_trapz,
                                        &norm,
                                        ptr->error_message),
             ptr->error_message,
             ptr->error_message);
  
  /* divide W by norm so that \int W(tau) dtau = 1 */
  for (index_tau = 0; index_tau < tau_size; index_tau++) {
    selection[index_tau]/=norm;
  }

  return _SUCCESS_;
}

/**
 * This routine computes the transfer functions \f$ \Delta_l^{X} (k) \f$)
 * as a function of wavenumber k for a given mode, initial condition,
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
                                struct transfer_workspace * ptw,
                                struct precision * ppr,
                                struct perturbs * ppt,
                                struct transfers * ptr,
                                int index_q,
                                int index_md,
                                int index_ic,
                                int index_tt,
                                int index_l,
                                double l,
                                double q_max_bessel
                                ){

  /** Summary: */

  /** - define local variables */

  /* current wavenumber value */
  double q,k;

  /* value of transfer function */
  double transfer_function;

  /* whether to use the Limber approximation */
  short use_limber;
  
  /* return zero tranbsfer function if l is above l_max */ 
  if (index_l >= ptr->l_size_tt[index_md][index_tt]) { 
    
    ptr->transfer[index_md][((index_ic * ptr->tt_size[index_md] + index_tt)
                             * ptr->l_size[index_md] + index_l)
                            * ptr->q_size + index_q] = 0.;
    return _SUCCESS_;
  }
  
  if (ptr->transfer_verbose > 2)
    printf("Compute transfer for l=%d\n",(int)l);
  
  q = ptr->q[index_q];
  k = ptr->k[index_md][index_q];

  if (ptr->transfer_verbose > 3)
    printf("Compute transfer for l=%d k=%e type=%d\n",(int)l,q,index_tt);
                
  class_call(transfer_use_limber(ppr,
                                 ppt,
                                 ptr,
                                 q_max_bessel,
                                 index_md,
                                 index_tt,
                                 q,
                                 l,
                                 &use_limber),
             ptr->error_message,
             ptr->error_message);
  
  if (use_limber == _TRUE_) {
    
    class_call(transfer_limber(ptw->tau_size,
                               ptr,
                               index_md,
                               index_q,
                               l,
                               q,
                               ptw->tau0_minus_tau,
                               ptw->sources,
                               &transfer_function),
               ptr->error_message,
               ptr->error_message); 
    
  }
  else {
    class_call(transfer_integrate(
                                  ppt,
                                  ptr,
                                  ptw,
                                  index_q,
                                  index_md,
                                  index_tt,
                                  l,
                                  index_l,
                                  k,
                                  &transfer_function
                                  ),
               ptr->error_message,
               ptr->error_message);
  }

  /* store transfer function in transfer structure */
  ptr->transfer[index_md][((index_ic * ptr->tt_size[index_md] + index_tt)
                           * ptr->l_size[index_md] + index_l)
                          * ptr->q_size + index_q]
    = transfer_function;
  
  return _SUCCESS_;

}

int transfer_use_limber(
                        struct precision * ppr,
                        struct perturbs * ppt,
                        struct transfers * ptr,
                        double q_max_bessel,
                        int index_md,
                        int index_tt,
                        double q,
                        double l,
                        short * use_limber) {


  /* criterium for chosing between integration and Limber 
     must be implemented here */
  
  *use_limber = _FALSE_;

  if (q>q_max_bessel) {
    *use_limber = _TRUE_;
  }
  else {
    
    if _scalars_ {
      
        //TBC: in principle the Limber condition should be adapted to account for curvature effects

        if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (index_tt == ptr->index_tt_lcmb) && (l>ppr->l_switch_limber)) {
          *use_limber = _TRUE_;
        }
        else if ((ppt->has_cl_density == _TRUE_) && (index_tt >= ptr->index_tt_density) && (index_tt < ptr->index_tt_density+ppt->selection_num) && (l>=ppr->l_switch_limber_for_cl_density_over_z*ppt->selection_mean[index_tt-ptr->index_tt_density])) {
          if (ppt->selection != dirac) *use_limber = _TRUE_;
        }

        else if ((ppt->has_cl_lensing_potential == _TRUE_) && (index_tt >= ptr->index_tt_lensing) && (index_tt < ptr->index_tt_lensing+ppt->selection_num) && (l>=ppr->l_switch_limber_for_cl_density_over_z*ppt->selection_mean[index_tt-ptr->index_tt_lensing])) {
          *use_limber = _TRUE_;
        }
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
 * @param index_q               Input : index of wavenumber
 * @param interpolated_sources  Input: array of interpolated sources
 * @param ptw                   Input : pointer to transfer_workspace structure (allocated in transfer_init() to avoid numerous reallocation) 
 * @param trsf                  Output: transfer function \f$ \Delta_l(k) \f$ 
 * @return the error status
 */

int transfer_integrate(
                       struct perturbs * ppt,
                       struct transfers * ptr,
                       struct transfer_workspace *ptw,
                       int index_q,
                       int index_md,
                       int index_tt,
                       double l,
                       int index_l,
                       double k,
                       double * trsf
                       ) {

  /** Summary: */

  /** - define local variables */

  double * tau0_minus_tau = ptw->tau0_minus_tau;
  double * w_trapz = ptw->w_trapz;
  double * sources = ptw->sources;

  /* minimum value of \f$ (\tau0-\tau) \f$ at which \f$ j_l(k[\tau_0-\tau]) \f$ is known, given that \f$ j_l(x) \f$ is sampled above some finite value \f$ x_{\min} \f$ (below which it can be approximated by zero) */  
  double tau0_minus_tau_min_bessel;
    
  /* index in the source's tau list corresponding to the last point in the overlapping region between sources and bessels. Also the index of possible Bessel truncation. */
  int index_tau_max, index_tau_max_Bessel;

  double bessel, *radial_function;

  radial_function_type radial_type; 

  /** - find minimum value of (tau0-tau) at which \f$ j_l(k[\tau_0-\tau]) \f$ is known, given that \f$ j_l(x) \f$ is sampled above some finite value \f$ x_{\min} \f$ (below which it can be approximated by zero) */  
  //tau0_minus_tau_min_bessel = x_min_l/k; /* segmentation fault impossible, checked before that k != 0 */
  //printf("index_l=%d\n",index_l);
  if (ptw->sgnK==0){
    tau0_minus_tau_min_bessel = ptw->chi_at_phiminabs[index_l]/k; /* segmentation fault impossible, checked before that k != 0 */
  }
  else{
    tau0_minus_tau_min_bessel = ptw->chi_at_phiminabs[index_l]/sqrt(ptw->sgnK*ptw->K); 
  }
  /** - if there is no overlap between the region in which bessels and sources are non-zero, return zero */
  if (tau0_minus_tau_min_bessel >= tau0_minus_tau[0]) {
    *trsf = 0.;
    return _SUCCESS_;
  }

  /** - if there is an overlap: */

  /** Select radial function type: */
  class_call(transfer_select_radial_function(
                                             ppt,
                                             ptr,
                                             index_md,
                                             index_tt,
                                             &radial_type),
             ptr->error_message,
             ptr->error_message);
  
  /** -> trivial case: the source is a Dirac function and is sampled in only one point */
  if (ptw->tau_size == 1) {

    class_call(transfer_radial_function(
                                        ptw,
                                        ppt,
                                        ptr,
                                        k,
                                        index_l,
                                        1,
                                        &bessel,
                                        radial_type
                                        ),
               ptr->error_message,
               ptr->error_message);
    
    *trsf = sources[0] * bessel;
    return _SUCCESS_;      
  }

  /** -> other cases */

  /** (a) find index in the source's tau list corresponding to the last point in the overlapping region. After this step, index_tau_max can be as small as zero, but not negative. */ 
  index_tau_max = ptw->tau_size-1;
  while (tau0_minus_tau[index_tau_max] < tau0_minus_tau_min_bessel)
    index_tau_max--;
  /* Set index so we know if the truncation of the convolution integral is due to Bessel and not
     due to the source. */
  index_tau_max_Bessel = index_tau_max;

  /** (b) the source function can vanish at large $\f \tau \f$. Check if further points can be eliminated. After this step and if we did not return a null transfer function, index_tau_max can be as small as zero, but not negative. */
  while (sources[index_tau_max] == 0.) { 
    index_tau_max--;
    if (index_tau_max < 0) {
      *trsf = 0.;
      return _SUCCESS_;
    }
  }

  /** Compute the radial function: */
  radial_function = malloc(sizeof(double)*(index_tau_max+1));
  class_call(transfer_radial_function(
                                      ptw,
                                      ppt,
                                      ptr,
                                      k,
                                      index_l,
                                      index_tau_max+1,
                                      radial_function,
                                      radial_type
                                      ),
             ptr->error_message,
             ptr->error_message);

  /** Now we do most of the convolution integral: */
  class_call(array_trapezoidal_convolution(sources,
                                           radial_function,
                                           index_tau_max+1,
                                           w_trapz,
                                           trsf,
                                           ptr->error_message),
             ptr->error_message,
             ptr->error_message);

  /** This integral is correct for the case where no truncation has
      occured. If it has been truncated at som index_tau_max because
      f[index_tau_max+1]==0, it is still correct. The 'mistake' in using
      the wrong weight w_trapz[index_tau_max] is exactly compensated by the
      triangle we miss. However, for the Bessel cut off, we must subtract the
      wrong triangle and add the correct triangle */
  if ((index_tau_max!=(ptw->tau_size-1))&&(index_tau_max==index_tau_max_Bessel)){
    //Bessel truncation
    *trsf -= 0.5*(tau0_minus_tau[index_tau_max+1]-tau0_minus_tau_min_bessel)*
      radial_function[index_tau_max]*sources[index_tau_max];
  }
  

  free(radial_function);
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
 * @param index_q               Input : index of wavenumber
 * @param interpolated_sources  Input: array of interpolated sources
 * @param trsf                  Output: transfer function \f$ \Delta_l(k) \f$ 
 * @return the error status
 */

int transfer_limber(
                    int tau_size,
                    struct transfers * ptr,
                    int index_md,
                    int index_q,
                    double l,
                    double k,
                    double * tau0_minus_tau,
                    double * sources, /* array with argument interpolated_sources[index_q*ppt->tau_size+index_tau] */
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

  tau0_minus_tau_limber = (l+0.5)/k; //TBC: to be updated to include curvature effects

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
                                          sources[index_tau-1]*tau0_minus_tau[index_tau-1],
                                          sources[index_tau]*tau0_minus_tau[index_tau],
                                          sources[index_tau+1]*tau0_minus_tau[index_tau+1],
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
                                          sources[index_tau-1]*tau0_minus_tau[index_tau-1],
                                          sources[index_tau]*tau0_minus_tau[index_tau],
                                          sources[index_tau]*tau0_minus_tau[index_tau],
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

  tau0_minus_tau_limber = (l+0.5)/k;  //TBC: to be updated to include curvature effects

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
                                        sources[index_tau-1],
                                        sources[index_tau],
                                        sources[index_tau+1],
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

int transfer_can_be_neglected(
                              struct precision * ppr,
                              struct background * pba,
                              struct perturbs * ppt,
                              struct transfers * ptr,
                              int index_md,
                              int index_ic,
                              int index_tt,
                              double k,
                              double l,
                              short * neglect) {

  *neglect = _FALSE_;

  /* implement here some conditions, e.g.:
     if (k*l>ppr->trans_kl_max) *neglect = _TRUE_;
  */      

  if _scalars_ {

      if ((index_tt == ptr->index_tt_t0) && (l < (k-ppr->transfer_neglect_delta_k_S_t0)*pba->conformal_age)) *neglect = _TRUE_;
      else if ((index_tt == ptr->index_tt_t1) && (l < (k-ppr->transfer_neglect_delta_k_S_t1)*pba->conformal_age)) *neglect = _TRUE_;
      else if ((index_tt == ptr->index_tt_t2) && (l < (k-ppr->transfer_neglect_delta_k_S_t2)*pba->conformal_age)) *neglect = _TRUE_;
      else if ((index_tt == ptr->index_tt_e) && (l < (k-ppr->transfer_neglect_delta_k_S_e)*pba->conformal_age)) *neglect = _TRUE_;
      else if ((index_tt == ptr->index_tt_lcmb) && (l < (k-ppr->transfer_neglect_delta_k_S_lcmb)*pba->conformal_age)) *neglect = _TRUE_;
      
    }
  
  else if _vectors_ {
      
      if ((index_tt == ptr->index_tt_t1) && (l < (k-ppr->transfer_neglect_delta_k_V_t1)*pba->conformal_age)) *neglect = _TRUE_;
      else if ((index_tt == ptr->index_tt_t2) && (l < (k-ppr->transfer_neglect_delta_k_V_t2)*pba->conformal_age)) *neglect = _TRUE_;
      else if ((index_tt == ptr->index_tt_e) && (l < (k-ppr->transfer_neglect_delta_k_V_e)*pba->conformal_age)) *neglect = _TRUE_;
      else if ((index_tt == ptr->index_tt_b) && (l < (k-ppr->transfer_neglect_delta_k_V_b)*pba->conformal_age)) *neglect = _TRUE_;

    }

  else if _tensors_ {

      if ((index_tt == ptr->index_tt_t2) && (l < (k-ppr->transfer_neglect_delta_k_T_t2)*pba->conformal_age)) *neglect = _TRUE_;
      else if ((index_tt == ptr->index_tt_e) && (l < (k-ppr->transfer_neglect_delta_k_T_e)*pba->conformal_age)) *neglect = _TRUE_;
      else if ((index_tt == ptr->index_tt_b) && (l < (k-ppr->transfer_neglect_delta_k_T_b)*pba->conformal_age)) *neglect = _TRUE_;

    }

  return _SUCCESS_;

}

int transfer_radial_function(
                             struct transfer_workspace * ptw,
                             struct perturbs * ppt,
                             struct transfers * ptr,
                             double k,
                             int index_l,
                             int x_size,
                             double * radial_function,
                             radial_function_type radial_type
                             ){
  
  HyperInterpStruct * pHIS = &(ptw->HIS);
  double *chi = ptw->chi;
  double *cscKgen = ptw->cscKgen;
  double *cotKgen = ptw->cotKgen;
  int j;
  double *Phi, *dPhi, *d2Phi, *chireverse;
  double K=0.,k2=1.0;
  double sqrt_absK_over_k;
  double absK_over_k2;
  double factor, s0, s2, ssqrt3, si, ssqrt2, ssqrt2i;
  double l = (double)ptr->l[index_l];
  
  K = ptw->K;
  k2 = k*k;

  if (ptw->sgnK==0){
    /* This is the choice consistent with chi=k*(tau0-tau) and nu=1 */
    sqrt_absK_over_k = 1.0;
  }
  else {
    K=ptw->K;
    sqrt_absK_over_k = sqrt(ptw->sgnK*K)/k;
  }
  absK_over_k2 =sqrt_absK_over_k*sqrt_absK_over_k;

  Phi = malloc(sizeof(double)*x_size);
  dPhi = malloc(sizeof(double)*x_size);
  d2Phi = malloc(sizeof(double)*x_size);
  chireverse = malloc(sizeof(double)*x_size);
  //Reverse chi
  for (j=0; j<x_size; j++)
    chireverse[j] = chi[x_size-1-j];

  /** Debug region:
      hyperspherical_Hermite_interpolation_vector(pHIS, x_size, index_l, chireverse, Phi, dPhi, d2Phi, NULL, NULL);
      FILE * debugbessel = fopen("debugbessel.dat","a");
      fprintf(debugbessel,"#l = %d\n",ptw->HIS.l[index_l]);
      for (j=0; j<x_size; j++){
      fprintf(debugbessel,"%.16e %.16e %.16e %.16e\n",chireverse[j], Phi[j], dPhi[j], d2Phi[j]);
      }  
      fprintf(debugbessel,"#end of block\n");
      fclose(debugbessel);
  */

  switch (radial_type){
  case SCALAR_TEMPERATURE_0:
    hyperspherical_Hermite_interpolation_vector(pHIS, x_size, index_l, chireverse, Phi, NULL, NULL, NULL, NULL);
    for (j=0; j<x_size; j++)
      radial_function[x_size-1-j] = Phi[j];
    break;
  case SCALAR_TEMPERATURE_1:
    hyperspherical_Hermite_interpolation_vector(pHIS, x_size, index_l, chireverse, NULL, dPhi, NULL, NULL, NULL);
    for (j=0; j<x_size; j++)
      radial_function[x_size-1-j] = sqrt_absK_over_k*dPhi[j];
    break;
  case SCALAR_TEMPERATURE_2:
    hyperspherical_Hermite_interpolation_vector(pHIS, x_size, index_l, chireverse, Phi, NULL, d2Phi, NULL, NULL);
    s2 = sqrt(1.0-3.0*K/k2);
    factor = 1.0/(2.0*s2);
    for (j=0; j<x_size; j++)
      radial_function[x_size-1-j] = factor*(3*absK_over_k2*d2Phi[j]+Phi[j]);
    break;
  case SCALAR_POLARISATION_E:
    hyperspherical_Hermite_interpolation_vector(pHIS, x_size, index_l, chireverse, Phi, NULL, NULL, NULL, NULL);
    s2 = sqrt(1.0-3.0*K/k2);
    factor = sqrt(3.0/8.0*(l+2.0)*(l+1.0)*l*(l-1.0))/s2;
    for (j=0; j<x_size; j++)
      radial_function[x_size-1-j] = factor*cscKgen[x_size-1-j]*cscKgen[x_size-1-j]*Phi[j];
    break;
  case VECTOR_TEMPERATURE_1:
    hyperspherical_Hermite_interpolation_vector(pHIS, x_size, index_l, chireverse, Phi, NULL, NULL, NULL, NULL);
    s0 = sqrt(1.0+K/k2);
    factor = sqrt(0.5*l*(l+1))/s0;
    for (j=0; j<x_size; j++)
      radial_function[x_size-1-j] = factor*cscKgen[x_size-1-j]*Phi[j];
    break;
  case VECTOR_TEMPERATURE_2:
    hyperspherical_Hermite_interpolation_vector(pHIS, x_size, index_l, chireverse, Phi, dPhi, NULL, NULL, NULL);
    s0 = sqrt(1.0+K/k2);
    ssqrt3 = sqrt(1.0-2.0*K/k2);
    factor = sqrt(1.5*l*(l+1))/s0/ssqrt3;
    for (j=0; j<x_size; j++)
      radial_function[x_size-1-j] = factor*cscKgen[x_size-1-j]*(sqrt_absK_over_k*dPhi[j]-cotKgen[j]*Phi[j]);
    break;
  case VECTOR_POLARISATION_E:
    hyperspherical_Hermite_interpolation_vector(pHIS, x_size, index_l, chireverse, Phi, dPhi, NULL, NULL, NULL);
    s0 = sqrt(1.0+K/k2);
    ssqrt3 = sqrt(1.0-2.0*K/k2);
    factor = 0.5*sqrt((l-1.0)*(l+2.0))/s0/ssqrt3;
    for (j=0; j<x_size; j++)
      radial_function[x_size-1-j] = factor*cscKgen[x_size-1-j]*(cotKgen[j]*Phi[j]+sqrt_absK_over_k*dPhi[j]);
    break;
  case VECTOR_POLARISATION_B:
    hyperspherical_Hermite_interpolation_vector(pHIS, x_size, index_l, chireverse, Phi, NULL, NULL, NULL, NULL);
    s0 = sqrt(1.0+K/k2);
    ssqrt3 = sqrt(1.0-2.0*K/k2);
    si = sqrt(1.0+2.0*K/k2);
    factor = 0.5*sqrt((l-1.0)*(l+2.0))*si/s0/ssqrt3;
    for (j=0; j<x_size; j++)
      radial_function[x_size-1-j] = factor*cscKgen[x_size-1-j]*Phi[j];
    break;
  case TENSOR_TEMPERATURE_2:
    hyperspherical_Hermite_interpolation_vector(pHIS, x_size, index_l, chireverse, Phi, NULL, NULL, NULL, NULL);
    ssqrt2 = sqrt(1.0-1.0*K/k2);
    si = sqrt(1.0+2.0*K/k2);
    factor = sqrt(3.0/8.0*(l+2.0)*(l+1.0)*l*(l-1.0))/si/ssqrt2;
    for (j=0; j<x_size; j++) {
      radial_function[x_size-1-j] = factor*cscKgen[x_size-1-j]*cscKgen[x_size-1-j]*Phi[j];
      /*
      if (ptr->l[index_l]<10) {
        class_call(HypersphericalExplicit(ptw->sgnK,(int)ptr->l[index_l],pHIS->beta,chireverse[j],&onePhi),
                   ptr->error_message,
                   ptr->error_message);
        radial_function[x_size-1-j] = factor*cscKgen[x_size-1-j]*cscKgen[x_size-1-j]*onePhi;

      }
      */
    }
    break;
  case TENSOR_POLARISATION_E:
    hyperspherical_Hermite_interpolation_vector(pHIS, x_size, index_l, chireverse, Phi, dPhi, d2Phi, NULL, NULL);
    ssqrt2 = sqrt(1.0-1.0*K/k2);
    si = sqrt(1.0+2.0*K/k2);
    factor = 0.25/si/ssqrt2;
    for (j=0; j<x_size; j++)
      radial_function[x_size-1-j] = factor*(absK_over_k2*d2Phi[j]+4.0*cotKgen[x_size-1-j]*sqrt_absK_over_k*dPhi[j]-
                                        (1.0+4*K/k2-2.0*cotKgen[x_size-1-j]*cotKgen[x_size-1-j])*Phi[j]);
    break;
  case TENSOR_POLARISATION_B:
    hyperspherical_Hermite_interpolation_vector(pHIS, x_size, index_l, chireverse, Phi, dPhi, NULL, NULL, NULL);
    ssqrt2i = sqrt(1.0+3.0*K/k2);
    ssqrt2 = sqrt(1.0-1.0*K/k2);
    si = sqrt(1.0+2.0*K/k2);
    factor = 0.5*ssqrt2i/ssqrt2/si;
    for (j=0; j<x_size; j++)
      radial_function[x_size-1-j] = factor*(sqrt_absK_over_k*dPhi[j]+2.0*cotKgen[x_size-1-j]*Phi[j]);
    break;
  }
  free(Phi);
  free(dPhi);
  free(d2Phi);
  free(chireverse);

  return _SUCCESS_;
}

int transfer_select_radial_function(
                                    struct perturbs * ppt,
                                    struct transfers * ptr,
                                    int index_md,
                                    int index_tt,
                                    radial_function_type * radial_type
                                    ) {

  /* generic case leading to generic bessel function (it applies also to all nonCMB types: lcmb, density, lensing) */ 
  *radial_type = SCALAR_TEMPERATURE_0; 

  /* other specific cases */
  if _scalars_ {

      if (ppt->has_cl_cmb_temperature == _TRUE_) {
         
        if (index_tt == ptr->index_tt_t0) {
          *radial_type = SCALAR_TEMPERATURE_0;
        }       
        if (index_tt == ptr->index_tt_t1) {
          *radial_type = SCALAR_TEMPERATURE_1;
        }
        if (index_tt == ptr->index_tt_t2) {
          *radial_type = SCALAR_TEMPERATURE_2;
        }
        
      }
              
      if (ppt->has_cl_cmb_polarization == _TRUE_) {
        
        if (index_tt == ptr->index_tt_e) {
          *radial_type = SCALAR_POLARISATION_E;
        }
        
      }
    }
          
  if _vectors_ {
      
      if (ppt->has_cl_cmb_temperature == _TRUE_) {
        
        if (index_tt == ptr->index_tt_t1) {
          *radial_type = VECTOR_TEMPERATURE_1;
        }
        if (index_tt == ptr->index_tt_t2) {
          *radial_type = VECTOR_TEMPERATURE_2;
        }
      }
      
      if (ppt->has_cl_cmb_polarization == _TRUE_) {
                
        if (index_tt == ptr->index_tt_e) {
          *radial_type = VECTOR_POLARISATION_E;
        }
        if (index_tt == ptr->index_tt_b) {
          *radial_type = VECTOR_POLARISATION_B;
        }
      
      }
    }
          
  if _tensors_ {

      if (ppt->has_cl_cmb_temperature == _TRUE_) {

        if (index_tt == ptr->index_tt_t2) {
          *radial_type = TENSOR_TEMPERATURE_2;
        }
      }

      if (ppt->has_cl_cmb_polarization == _TRUE_) {

        if (index_tt == ptr->index_tt_e) {
          *radial_type = TENSOR_POLARISATION_E;
        }
        if (index_tt == ptr->index_tt_b) {
          *radial_type = TENSOR_POLARISATION_B;
        }

      }
    }

  return _SUCCESS_;

} 

int transfer_init_HIS_from_bessel(
                                  struct bessels * pbs,
                                  struct transfers * ptr,
                                  HyperInterpStruct *pHIS
                                  ) {

  /* The purpose of this function is to populate a flat Hyperspherical
     Interpolation Structure from the Bessel structure. The Bessel
     structure must have been computed on an l-independent x-grid.
  */
  int x_size, index_l, index_x;
  double xmin, x;
  
  /* Check the Bessels structure for compatibility: */
  x_size=pbs->x_size[0];
  xmin = *(pbs->x_min[0]);
  for (index_l=1; index_l<pbs->l_size; index_l++){
    class_test(x_size!=pbs->x_size[index_l],
               ptr->error_message,
               "pbs->x_size[%d]=%d, not equal to pbs->x_size[0]=%d.\n",
               index_l,pbs->x_size[index_l],x_size);
    class_test(fabs(xmin-*(pbs->x_min[index_l]))>100*DBL_EPSILON,
               ptr->error_message,
               "*(pbs->xmin[%d])=%.16e, should have been %.16e.\n",
               index_l,*(pbs->x_min[index_l]),xmin);
  }
  /* Set and copy all fields of HIS: */
  pHIS->K = 0;
  pHIS->beta = 1;
  pHIS->delta_x = pbs->x_step;
  pHIS->x_size = x_size;
  pHIS->l_size = pbs->l_size;
  pHIS->trig_order=5;
  /* Allocate storage in pHIS: */
  class_alloc(pHIS->l,hyperspherical_HIS_size(pHIS->l_size, pHIS->x_size),ptr->error_message);
  /* Set pointers: */
  hyperspherical_update_pointers(pHIS, pHIS->l);
  
  for(index_x=0; index_x<x_size; index_x++){
    x = xmin+index_x*pbs->x_step;
    pHIS->x[index_x] = x;
    pHIS->sinK[index_x] = x;
    pHIS->cotK[index_x] = 1.0/x;
  }
  memcpy(pHIS->l, pbs->l, pbs->l_size*sizeof(int));
  for (index_l=0; index_l<pbs->l_size; index_l++)
    memcpy(pHIS->phi+index_l*x_size,pbs->j[index_l],sizeof(double)*x_size);
  for (index_l=0; index_l<pbs->l_size; index_l++)
    memcpy(pHIS->dphi+index_l*x_size,pbs->dj[index_l],sizeof(double)*x_size);
  return _SUCCESS_;
}


int transfer_workspace_init(
                            struct transfers * ptr,
                            struct transfer_workspace **ptw,
                            int perturb_tau_size,
                            int tau_size_max,
                            int get_HIS_from_pbs,
                            double K,
                            int sgnK,
                            int get_HIS_from_shared_memory){

  class_calloc(*ptw,1,sizeof(struct transfer_workspace),ptr->error_message);

  (*ptw)->tau_size_max = tau_size_max;
  (*ptw)->l_size = ptr->l_size_max;
  (*ptw)->HIS_allocated=_FALSE_;
  (*ptw)->HIS_attached=_FALSE_;
  (*ptw)->get_HIS_from_pbs=get_HIS_from_pbs;
  (*ptw)->get_HIS_from_shared_memory=get_HIS_from_shared_memory;

  (*ptw)->K = K;
  (*ptw)->sgnK = sgnK;

  class_alloc((*ptw)->chi_at_phiminabs,(*ptw)->l_size*sizeof(double),ptr->error_message);
  class_alloc((*ptw)->interpolated_sources,perturb_tau_size*sizeof(double),ptr->error_message);
  class_alloc((*ptw)->sources,tau_size_max*sizeof(double),ptr->error_message);
  class_alloc((*ptw)->tau0_minus_tau,tau_size_max*sizeof(double),ptr->error_message);
  class_alloc((*ptw)->w_trapz,tau_size_max*sizeof(double),ptr->error_message);
  class_alloc((*ptw)->chi,tau_size_max*sizeof(double),ptr->error_message);
  class_alloc((*ptw)->cscKgen,tau_size_max*sizeof(double),ptr->error_message);
  class_alloc((*ptw)->cotKgen,tau_size_max*sizeof(double),ptr->error_message);
 
  return _SUCCESS_;
}

int transfer_workspace_free(
                            struct transfers * ptr,
                            struct transfer_workspace *ptw
                            ) {

  HyperInterpStruct *pHIS_shared;

  if (ptw->HIS_allocated==_TRUE_){
    //Free HIS structure:
    hyperspherical_HIS_free(&(ptw->HIS));
  }
  if (ptw->HIS_attached==_TRUE_){
    //detach HIS from shared memory zone:
    pHIS_shared = (HyperInterpStruct *)ptw->HIS.l;
    shmdt((pHIS_shared-1));
  }
  free(ptw->chi_at_phiminabs);
  free(ptw->interpolated_sources);
  free(ptw->sources);
  free(ptw->tau0_minus_tau);
  free(ptw->w_trapz);
  free(ptw->chi);
  free(ptw->cscKgen);
  free(ptw->cotKgen);
  
  free(ptw);
  return _SUCCESS_;
}

int transfer_update_HIS( 
                        struct precision * ppr,
                        struct bessels * pbs,
                        struct transfers * ptr,
                        struct transfer_workspace * ptw,
                        int index_q,
                        double tau0
                         ) {

  double nu,new_nu;
  int int_nu;
  double xmin, xmax, sampling, phiminabs, xtol;
  double sqrt_absK;
  int l_size_max, l_WKB;
  int index_l_left,index_l_right;
  int shmid;
  HyperInterpStruct *pHIS_shared;

  if (ptw->sgnK==0)
    index_q=0;


  if ((ptw->get_HIS_from_shared_memory==_TRUE_)&&
      (ptr->initialise_HIS_cache==_FALSE_)){

    if ((ptw->HIS_attached==_TRUE_) && (ptw->sgnK==0)) return _SUCCESS_;
      
    /** Try to connect to shared memory region containing HIS structure: */
    shmid = shmget(_SHARED_MEMORY_KEYS_START_+index_q, 
                   sizeof(HyperInterpStruct),
                   S_IRUSR | S_IWUSR);
    if (shmid>0){

      if (ptw->HIS_attached==_TRUE_){
        //detach HIS from shared memory zone:
        pHIS_shared = (HyperInterpStruct *)ptw->HIS.l;
        shmdt(pHIS_shared-1);
      }

      /** Now try to attach the region: */
      pHIS_shared = shmat (shmid, 0, 0);
      ptw->HIS_attached=_TRUE_;
      if (pHIS_shared>0){
        /** Eventually do a few more checks here, perhaps macros. */
        //fprintf(stderr,"Connected to HIS struct!\n");
        /** Copy HIS structure from shared HIS to ptw->HIS: */
        memcpy(&(ptw->HIS),pHIS_shared,sizeof(HyperInterpStruct));
        /** Update pointers in ptw->HIS:*/
        hyperspherical_update_pointers(&(ptw->HIS),pHIS_shared+1);
        /** Find x_at_phimin: */
        phiminabs = ppr->hyper_phi_min_abs;
        xtol = ppr->hyper_x_tol;
        
        class_call(hyperspherical_get_xmin(&(ptw->HIS),
                                           xtol,
                                           phiminabs/1000.,
                                           ptw->chi_at_phiminabs),
                   ptr->error_message, 
                   ptr->error_message);
        return _SUCCESS_;
      }
    }
    else {
      fprintf(stderr,"Unable to connect to shared memory\n");
    }
  }

    
  if (ptw->HIS_allocated==_TRUE_){

    if (ptw->sgnK == 0) {
      /* nothing to be done */
      return _SUCCESS_;
    }
    else {
      /** free HIS */
      class_call(hyperspherical_HIS_free(&(ptw->HIS)),
                 ptr->error_message,
                 ptr->error_message);
      ptw->HIS_allocated = _FALSE_;
    }
  }

  if (ptw->get_HIS_from_pbs==_TRUE_){
    class_call(transfer_init_HIS_from_bessel(pbs,
                                             ptr,
                                             &(ptw->HIS)),
               ptr->error_message,
               ptr->error_message);
    ptw->HIS_allocated = _TRUE_;
  }
  else{
    //These number should be set from input structures in the future:
    //printf("Creating interpolation structure...\n");
    xmin = ppr->hyper_x_min; //Will be changed to _HYPER_SAFETY_ by routine
    xmin = max(xmin,_HYPER_SAFETY_);
    if (ptw->sgnK == 0) {
      xmax = ptr->q[ptr->q_size-1]*tau0; // x_max = k_max * tau0
      nu=1.;
      sampling = ppr->hyper_sampling_flat;
      l_size_max = ptr->l_size_max;
      sqrt_absK=0.0;
    }
    else {
      sqrt_absK = sqrt(ptw->sgnK*ptw->K);
    
      xmax = sqrt_absK*tau0;
      nu = ptr->q[index_q]/sqrt_absK;
    
      if (ptw->sgnK == 1) {
        xmax = min(xmax,_PI_/2.0-_HYPER_SAFETY_); //We only need solution on [0;pi/2]
      
        int_nu = (int)(nu+0.2);
        new_nu = (double)int_nu;
        class_test(nu-new_nu > 1.e-6,
                   ptr->error_message,
                   "problem in q list definition in closed case for index_q=%d, nu=%e, nu-int(nu)=%e",index_q,nu,nu-new_nu);
        nu = new_nu;
        
      }
    }
    sampling = ppr->hyper_sampling_curved;
    
  
    /* find the highest value of l such that x_nonzero < xmax = sqrt(|K|) tau0. That will be l_max. */
    l_size_max = ptr->l_size_max;
    if (ptw->sgnK == 1)
      while ((double)ptr->l[l_size_max-1] >= nu)
        l_size_max--;
    
    if (ptw->sgnK == -1){
      xtol = ppr->hyper_x_tol;
      phiminabs = ppr->hyper_phi_min_abs;
  
      /** First try to find lmax using fast approximation: */
      index_l_left=0;
      index_l_right=l_size_max-1;
      class_call(transfer_get_lmax(hyperspherical_get_xmin_from_approx,
                                   ptw->sgnK,
                                   nu,
                                   ptr->l,
                                   l_size_max,
                                   phiminabs,
                                   xmax,
                                   xtol,
                                   &index_l_left,
                                   &index_l_right,
                                   ptr->error_message),
                 ptr->error_message,
                 ptr->error_message);
  
      /** Now use WKB approximation to eventually modify borders: */
      class_call(transfer_get_lmax(hyperspherical_get_xmin_from_Airy,
                                   ptw->sgnK,
                                   nu,
                                   ptr->l,
                                   l_size_max,
                                   phiminabs,
                                   xmax,
                                   xtol,
                                   &index_l_left,
                                   &index_l_right,
                                   ptr->error_message),
                 ptr->error_message,
                 ptr->error_message);
      l_size_max = index_l_right+1;
    }
  
    while ((ptr->l[l_size_max]> 750./nu ) && (ptr->l[l_size_max] > 45 ))
      l_size_max--;

    class_test(nu <= 0.,
               ptr->error_message,
               "nu=%e when index_q=%d, q=%e, K=%e, sqrt(|K|)=%e; instead nu should always be strictly positive",
               nu,index_q,ptr->q[index_q],ptw->K,sqrt_absK);
    

    //fprintf(stderr,"%d %d %d %e\n",ptr->l_size_max,l_size_max,ptr->l[l_size_max-1],nu);
    l_WKB = 46;
    class_call(hyperspherical_HIS_create(ptw->sgnK,
                                         nu,
                                         l_size_max,
                                         ptr->l,
                                         xmin,
                                         xmax,
                                         sampling,
                                         l_WKB,
                                         &(ptw->HIS),
                                         ptr->error_message),
               ptr->error_message,
               ptr->error_message);
    ptw->HIS_allocated = _TRUE_; 
  }

  if (ptr->initialise_HIS_cache==_TRUE_){
    //Try to allocate shared memory segment and copy HIS structure:
    shmid = shmget(_SHARED_MEMORY_KEYS_START_+index_q, 
                   sizeof(HyperInterpStruct)+
                   hyperspherical_HIS_size(ptw->HIS.l_size, ptw->HIS.x_size), 
                   IPC_CREAT | IPC_EXCL | S_IRUSR | S_IWUSR);
    //printf("shmid = %d\n",shmid);
    if (shmid<0){
      fprintf(stderr,"Could not allocate HIS\n");
      //printf("errno = %d. (EEXIST=%d)\n",errno,EEXIST);
    }
    else{
      /*
      printf("Allocating HIS structure for key:%08X with size %dkB\n",
             _SHARED_MEMORY_KEYS_START_+index_q,
             (int)(hyperspherical_HIS_size(ptw->HIS.l_size, ptw->HIS.x_size)/1024));
      */
      pHIS_shared = shmat (shmid, 0, 0);
      //Memcopy: First structure, then storage
      memcpy(pHIS_shared,&(ptw->HIS),sizeof(HyperInterpStruct));
      memcpy((pHIS_shared+1), ptw->HIS.l,
             hyperspherical_HIS_size(ptw->HIS.l_size,ptw->HIS.x_size));
      shmdt(pHIS_shared);
    }
  }
  

  //For each l, find lowest x such that |phi(x)| (or |j_l(x)| = phiminabs.
  phiminabs = ppr->hyper_phi_min_abs;
  xtol = ppr->hyper_x_tol;

  class_call(hyperspherical_get_xmin(&(ptw->HIS),
                                     xtol,
                                     phiminabs/1000.,
                                     ptw->chi_at_phiminabs),
             ptr->error_message, 
             ptr->error_message);
  
  return _SUCCESS_;
}

int transfer_get_lmax(int (*get_xmin_generic)(int sgnK,
                                              int l,
                                              double nu,
                                              double xtol,
                                              double phiminabs,
                                              double *x_nonzero,
                                              int *fevals),
                      int sgnK,
                      double nu,
                      int *lvec,
                      int lsize,
                      double phiminabs,
                      double xmax,
                      double xtol,
                      int *index_l_left,
                      int *index_l_right,
                      ErrorMsg error_message){
  double x_nonzero;
  int fevals=0, index_l_mid;
  int multiplier;
  int right_boundary_checked = _FALSE_;
  int hil=0,hir=0,bini=0;
  class_call(get_xmin_generic(sgnK,
                              lvec[0],
                              nu,
                              xtol,
                              phiminabs,
                              &x_nonzero,
                              &fevals),
             error_message,
             error_message);
  if (x_nonzero >= xmax){
    //printf("None relevant\n");
    //x at left boundary is already larger than xmax. 
    *index_l_right = max(lsize-1,1);
    return _SUCCESS_;
  }
  class_call(get_xmin_generic(sgnK,
                              lvec[lsize-1],
                              nu,
                              xtol,
                              phiminabs,
                              &x_nonzero,
                              &fevals),
             error_message,
             error_message);
  
  if (x_nonzero < xmax){
    //All Bessels are relevant
    //printf("All relevant\n");
    *index_l_left = max(0,(lsize-2));
    return _SUCCESS_;
  }
  /** Hunt for left boundary: */
  for (multiplier=1; ;multiplier *= 5){
    hil++;
    class_call(get_xmin_generic(sgnK,
                                lvec[*index_l_left],
                                nu,
                                xtol,
                                phiminabs,
                                &x_nonzero,
                                &fevals),
               error_message,
               error_message);
    //printf("Hunt left, iter = %d, x_nonzero=%g\n",hil,x_nonzero);
    if (x_nonzero <= xmax){
      //Boundary found
      break;
    }
    else{
      //We can use current index_l_left as index_l_right:
      *index_l_right = *index_l_left;
      right_boundary_checked = _TRUE_;
    }
    //Update index_l_left:
    *index_l_left = (*index_l_left)-multiplier;
    if (*index_l_left<=0){
      *index_l_left = 0;
      break;
    }
  }
  /** If not found, hunt for right boundary: */
  if (right_boundary_checked == _FALSE_){
    for (multiplier=1; ;multiplier *= 5){
      hir++;
      //printf("right iteration %d,index_l_right:%d\n",hir,*index_l_right);
      class_call(get_xmin_generic(sgnK,
                                  lvec[*index_l_right],
                                  nu,
                                  xtol,
                                  phiminabs,
                                  &x_nonzero,
                                  &fevals),
                 error_message,
                 error_message);
      if (x_nonzero >= xmax){
        //Boundary found
        break;
      }
      else{
        //We can use current index_l_right as index_l_left:
        *index_l_left = *index_l_right;
       }
      //Update index_l_right:
      *index_l_right = (*index_l_right)+multiplier;
      if (*index_l_right>=(lsize-1)){
        *index_l_right = lsize-1;
        break;
      }
    }
  }
  //  int fevalshunt=fevals;
  fevals=0;
  //Do binary search
  //  printf("Do binary search in get_lmax. \n");
  //printf("Region: [%d, %d]\n",*index_l_left,*index_l_right);
  while (((*index_l_right) - (*index_l_left)) > 1) {
    bini++;
    index_l_mid= (int)(0.5*((*index_l_right)+(*index_l_left)));
    //printf("left:%d, mid=%d, right=%d\n",*index_l_left,index_l_mid,*index_l_right);
    class_call(get_xmin_generic(sgnK,
                                lvec[index_l_mid],
                                nu,
                                xtol,
                                phiminabs,
                                &x_nonzero,
                                &fevals),
               error_message,
               error_message);
    if (x_nonzero < xmax) 
      *index_l_left=index_l_mid;
    else 
      *index_l_right=index_l_mid;
  }
  //printf("Done\n");
  /**  printf("Hunt left iter=%d, hunt right iter=%d (fevals: %d). For binary seach: %d (fevals: %d)\n",
         hil,hir,fevalshunt,bini,fevals);
  */
  return _SUCCESS_;
}
