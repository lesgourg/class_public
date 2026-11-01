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
 * -# transfer_init() at the beginning (but after perturbations_init()
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
#include "parallel.h"

/**
 * Transfer function \f$ \Delta_l^{X} (q) \f$ at a given wavenumber q.
 *
 * For a given mode (scalar, vector, tensor), initial condition, type
 * (temperature, polarization, lensing, etc) and multipole, computes
 * the transfer function for an arbitrary value of q by interpolating
 * between pre-computed values of q. This
 * function can be called from whatever module at whatever time,
 * provided that transfer_init() has been called before, and
 * transfer_free() has not been called yet.
 *
 * Wavenumbers are called q in this module and k in the perturbation
 * module. In flat universes k=q. In non-flat universes q and k differ
 * through \f$ q2 = k2 + K(1+m)\f$, where m=0,1,2 for scalar, vector,
 * tensor. q should be used throughout the transfer module, excepted
 * when interpolating or manipulating the source functions S(k,tau)
 * calculated in the perturbation module: for a given value of q, this
 * should be done at the corresponding k(q).
 *
 * @param ptr        Input: pointer to transfer structure
 * @param index_md   Input: index of requested mode
 * @param index_ic   Input: index of requested initial condition
 * @param index_tt   Input: index of requested type
 * @param index_l    Input: index of requested multipole
 * @param q          Input: any wavenumber
 * @param transfer_function Output: transfer function
 * @return the error status
 */

int transfer_functions_at_q(
                            struct transfer * ptr,
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
 * This routine initializes the transfer structure, (in particular,
 * computes table of transfer functions \f$ \Delta_l^{X} (q) \f$)
 *
 * Main steps:
 *
 * - initialize all indices in the transfer structure, define the
 *   properties of each transfer type, and allocate all its arrays
 *   using transfer_indices().
 *
 * - define once for all the time sampling of the line-of-sight
 *   integral of each transfer type with transfer_time_samplings(),
 *   and the window function by which each perturbation source must be
 *   multiplied with transfer_precompute_windows()
 *
 * - loop over q values. For each q, initialize the fields of a memory zone called the transfer_workspace with transfer_workspace_init(), compute the Bessel functions if needed with transfer_update_HIS(), defer the calculation of all transfer functions to transfer_compute_for_each_q(), and free the the workspace with transfer_workspace_free()
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param pth Input: pointer to thermodynamics structure
 * @param ppt Input: pointer to perturbation structure
 * @param pfo Input: pointer to fourier structure
 * @param ptr Output: pointer to initialized transfer structure
 * @return the error status
 */

int transfer_init(
                  struct precision * ppr,
                  struct background * pba,
                  struct thermodynamics * pth,
                  struct perturbations * ppt,
                  struct fourier * pfo,
                  struct transfer * ptr
                  ) {

  /** Summary: */

  /** - define local variables */

  /* running index for wavenumbers */
  int index_q;

  /* conformal time today */
  double tau0;
  /* conformal time at recombination */
  double tau_rec;
  /* order of magnitude of the oscillation period of transfer functions */
  double q_period;

  /* array of sources S(k,tau), just taken from perturbation module,
     or transformed if non-linear corrections are needed
     sources[index_md][index_ic * ppt->tp_size[index_md] + index_tp][index_tau * ppt->k_size[index_md] + index_k]
  */
  double *** sources;

  /* array of source derivatives S''(k,tau)
     (second derivative with respect to k, not tau!),
     used to interpolate sources at the right values of k,
     sources_spline[index_md][index_ic * ppt->tp_size[index_md] + index_tp][index_tau * ppt->k_size[index_md] + index_k]
  */
  double *** sources_spline;

  /* structure containing the flat spherical bessel functions */

  HyperInterpStruct BIS;
  double xmax;


  /** - check whether any spectrum in harmonic space (i.e., any \f$C_l\f$'s) is actually requested */

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

  /** - check whether we will need the full Limber q-sampling scheme */

  if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (ppt->want_lcmb_full_limber == _TRUE_)) {
    ptr->do_lcmb_full_limber = _TRUE_;
  }
  else {
    ptr->do_lcmb_full_limber = _FALSE_;
  }

  /** - get number of modes (scalars, tensors...) */

  ptr->md_size = ppt->md_size;

  /** - get conformal age / recombination time
      from background / thermodynamics structures
      (only place where these structures are used in this module) */

  tau0 = pba->conformal_age;
  tau_rec = pth->tau_rec;

  /** - correspondence between k and l depend on angular diameter
      distance, i.e. on curvature. */

  ptr->angular_rescaling = pth->angular_rescaling;

  /** - order of magnitude of the oscillation period of transfer functions */

  q_period = 2.*_PI_/(tau0-tau_rec)*ptr->angular_rescaling;

  /** - initialize all indices in the transfer structure, define the
      properties of each transfer type, and allocate all its arrays
      using transfer_indices() */

  class_call(transfer_indices(ppr,ppt,ptr,q_period,pba->K,pba->sgnK),
             ptr->error_message,
             ptr->error_message);

  /** - copy sources to a local array sources (in fact, only the pointers are copied, not the data), and eventually apply non-linear corrections to the sources */

  class_alloc(sources,
              ptr->md_size*sizeof(double**),
              ptr->error_message);

  class_call(transfer_perturbation_copy_sources_and_nl_corrections(ppt,pfo,ptr,sources),
             ptr->error_message,
             ptr->error_message);

  /** - spline all the sources passed by the perturbation module with respect to k (in order to interpolate later at a given value of k) */

  class_alloc(sources_spline,
              ptr->md_size*sizeof(double**),
              ptr->error_message);

  class_call(transfer_perturbation_source_spline(ppt,ptr,sources,sources_spline),
             ptr->error_message,
             ptr->error_message);

  /** - compute flat spherical bessel functions */

  xmax = ptr->q[ptr->q_size-1]*tau0;
  // If the use the flat approximation for K<0, we need to increase xmax artificially
  if (pba->sgnK == -1 && ptr->index_q_flat_approximation < ptr->q_size)
    xmax *= sqrt(pba->sgnK*pba->K)/ptr->q[ptr->q_size-1]*(ptr->l[ptr->l_size_max-1]+1)/asinh((ptr->l[ptr->l_size_max-1]+1)/ptr->q[ptr->q_size-1]*sqrt(pba->sgnK*pba->K))*1.01;

  class_call(hyperspherical_HIS_create(0,
                                       1.,
                                       ptr->l_size_max,
                                       ptr->l,
                                       ppr->hyper_x_min,
                                       xmax,
                                       ppr->hyper_sampling_flat,
                                       ptr->l[ptr->l_size_max-1]+1,
                                       ppr->hyper_phi_min_abs,
                                       &BIS,
                                       ptr->error_message),
             ptr->error_message,
             ptr->error_message);

  /** - eventually read the selection and evolution functions */

  class_call(transfer_global_selection_read(ptr),
             ptr->error_message,
             ptr->error_message);

  /** - eventually read the selection functions of each bin in a file */

  ptr->selection_file_z_size = 0;

  if (((ppt->has_cl_lensing_potential == _TRUE_) || (ppt->has_cl_number_count == _TRUE_)) && (ppt->selection == tabulated)) {
    class_call(transfer_tabulated_selection_read(ppt,ptr),
               ptr->error_message,
               ptr->error_message);
  }

  /** - find the conformal times delimiting the selection function of each redshift bin */

  if ((ppt->has_cl_lensing_potential == _TRUE_) || (ppt->has_cl_number_count == _TRUE_)) {
    class_call(transfer_selection_times(ppr,pba,ppt,ptr),
               ptr->error_message,
               ptr->error_message);
  }

  /** - define the time sampling of the line-of-sight integral of each transfer type */

  class_call(transfer_time_samplings(ppr,pba,ppt,ptr,tau_rec),
             ptr->error_message,
             ptr->error_message);

  /** - compute the window function W(tau) by which the perturbation source of each transfer type must be multiplied */

  class_call(transfer_precompute_windows(ppr,pba,ppt,ptr,tau_rec),
             ptr->error_message,
             ptr->error_message);

  class_setup_parallel();

  /** - loop over all wavenumbers (parallelized).*/
  /* For each wavenumber: */
  for (index_q = 0; index_q < MAX(ptr->q_size,ptr->q_size_limber); index_q++) {
 class_run_parallel(with_arguments(pba,pth,ppt,ptr,ppr,index_q,tau_rec,sources,sources_spline,tau0,&BIS),

      /* compute the transfer functions in the normal case (not with the
         full Limber q-sampling scheme) */

      struct transfer_workspace tw;
      struct transfer_workspace * ptw = &tw;

      class_call(transfer_workspace_init(ptr,
                                         ppt,
                                         ptw,
                                         pba->K,
                                         pba->sgnK,
                                         tau0-pth->tau_cut,
                                         &BIS),
                 ptr->error_message,
                 ptr->error_message);

      if (index_q < ptr->q_size) {

        if (ptr->transfer_verbose > 2)
        printf("Compute transfer for wavenumber [%d/%zu]\n",index_q,ptr->q_size-1);

        /* Update interpolation structure: */
        class_call(transfer_update_HIS(ppr,
                                       ptr,
                                       ptw,
                                       index_q,
                                       tau0),
                   ptr->error_message,
                   ptr->error_message);

        class_call(transfer_compute_for_each_q(ppr,
                                               pba,
                                               ppt,
                                               ptr,
                                               index_q,
                                               tau_rec,
                                               sources,
                                               sources_spline,
                                               ptw,
                                               _FALSE_),
                   ptr->error_message,
                   ptr->error_message);
      }

      /* compute the transfer functions in the full Limber q-sampling
       scheme (if this scheme is not needed, ptr->q_size_limber=0 and the
       condition is never met) */

      if (index_q < ptr->q_size_limber) {

        class_call(transfer_compute_for_each_q(ppr,
                                               pba,
                                               ppt,
                                               ptr,
                                               index_q,
                                               tau_rec,
                                               sources,
                                               sources_spline,
                                               ptw,
                                               _TRUE_),
                            ptr->error_message,
                            ptr->error_message);
      }

      class_call(transfer_workspace_free(ptr,ptw),
                 ptr->error_message,
                 ptr->error_message);
      return _SUCCESS_;
    );
  } /* end of loop over wavenumber */

  class_finish_parallel();

  /** - finally, free arrays allocated outside parallel zone */

  class_call(transfer_perturbation_sources_spline_free(ppt,ptr,sources_spline),
             ptr->error_message,
             ptr->error_message);

  class_call(transfer_perturbation_sources_free(ppt,pfo,ptr,sources),
             ptr->error_message,
             ptr->error_message);

  class_call(hyperspherical_HIS_free(&BIS,ptr->error_message),
             ptr->error_message,
             ptr->error_message);

  ptr->is_allocated = _TRUE_;
  return _SUCCESS_;
}

/**
 * This routine frees all the memory space allocated by transfer_init().
 *
 * To be called at the end of each run, only when no further calls to
 * transfer_functions_at_k() are needed.
 *
 * @param ptr Input: pointer to transfer structure (which fields must be freed)
 * @return the error status
 */

int transfer_free(
                  struct transfer * ptr
                  ) {

  int index_md;
  int index_ts;
  int index_tt;

  if (ptr->has_cls == _TRUE_) {

    for (index_md = 0; index_md < ptr->md_size; index_md++) {
      free(ptr->l_size_tt[index_md]);
      free(ptr->transfer[index_md]);
      free(ptr->k[index_md]);
      if (ptr->do_lcmb_full_limber == _TRUE_) {
        free(ptr->transfer_limber[index_md]);
        free(ptr->k_limber[index_md]);
      }

      for (index_ts = 0; index_ts < ptr->ts_size[index_md]; index_ts++) {
        free(ptr->tau0_minus_tau[index_md][index_ts]);
        free(ptr->w_trapz[index_md][index_ts]);
        free(ptr->index_tau_pt[index_md][index_ts]);
        free(ptr->weight_tau_pt[index_md][index_ts]);
      }
      free(ptr->tau0_minus_tau[index_md]);
      free(ptr->w_trapz[index_md]);
      free(ptr->index_tau_pt[index_md]);
      free(ptr->weight_tau_pt[index_md]);
      free(ptr->q_max_sampled[index_md]);
      free(ptr->tau_size_ts[index_md]);

      for (index_tt = 0; index_tt < ptr->tt_size[index_md]; index_tt++) {
        free(ptr->window[index_md][index_tt]);
      }
      free(ptr->window[index_md]);

      free(ptr->tp_of_tt[index_md]);
      free(ptr->bin_of_tt[index_md]);
      free(ptr->radial_type[index_md]);
      free(ptr->window_type[index_md]);
      free(ptr->k_power[index_md]);
      free(ptr->k_pivot[index_md]);
      free(ptr->l_limber[index_md]);
      free(ptr->delta_k_neglect[index_md]);
      free(ptr->has_late_source_cut[index_md]);
      free(ptr->has_full_limber[index_md]);
      free(ptr->ts_of_tt[index_md]);
    }

    free(ptr->tt_size);
    free(ptr->l_size_tt);
    free(ptr->l_size);
    free(ptr->l);
    free(ptr->q);
    free(ptr->k);
    free(ptr->transfer);
    if (ptr->do_lcmb_full_limber == _TRUE_) {
      free(ptr->q_limber);
      free(ptr->k_limber);
      free(ptr->transfer_limber);
    }

    free(ptr->tau0_minus_tau);
    free(ptr->w_trapz);
    free(ptr->index_tau_pt);
    free(ptr->weight_tau_pt);
    free(ptr->q_max_sampled);
    free(ptr->tau_size_ts);
    free(ptr->ts_size);
    free(ptr->window);

    free(ptr->tp_of_tt);
    free(ptr->bin_of_tt);
    free(ptr->radial_type);
    free(ptr->window_type);
    free(ptr->k_power);
    free(ptr->k_pivot);
    free(ptr->l_limber);
    free(ptr->delta_k_neglect);
    free(ptr->has_late_source_cut);
    free(ptr->has_full_limber);
    free(ptr->ts_of_tt);

    if (ptr->nz_size > 0) {
      free(ptr->nz_z);
      free(ptr->nz_nz);
      free(ptr->nz_ddnz);
    }

    if (ptr->nz_evo_size > 0) {
      free(ptr->nz_evo_z);
      free(ptr->nz_evo_nz);
      free(ptr->nz_evo_dlog_nz);
      free(ptr->nz_evo_dd_dlog_nz);
    }

    if (ptr->selection_file_z_size > 0) {
      free(ptr->selection_file_z);
      free(ptr->selection_file_dNdz);
      free(ptr->selection_file_dddNdz);
    }
  }

  ptr->is_allocated = _FALSE_;

  return _SUCCESS_;

}

/**
 * This routine defines all indices and allocates all tables
 * in the transfer structure
 *
 * Compute list of (k, l) values, allocate and fill corresponding
 * arrays in the transfer structure. Allocate the array of transfer
 * function tables.
 *
 * @param ppr      Input: pointer to precision structure
 * @param ppt      Input: pointer to perturbation structure
 * @param ptr      Input/Output: pointer to transfer structure
 * @param q_period Input: order of magnitude of the oscillation period of transfer functions
 * @param K        Input: spatial curvature (in absolute value)
 * @param sgnK     Input: spatial curvature sign (open/closed/flat)
 * @return the error status
 */

int transfer_indices(
                     struct precision * ppr,
                     struct perturbations * ppt,
                     struct transfer * ptr,
                     double q_period,
                     double K,
                     int sgnK
                     ) {

  /** Summary: */

  /** - define local variables */

  int index_md,index_tt,index_tt_common;

  /** - define indices for transfer types */

  class_alloc(ptr->tt_size,ptr->md_size * sizeof(int),ptr->error_message);

  /** - type indices common to scalars and tensors */

  index_tt = 0;

  class_define_index(ptr->index_tt_t2,ppt->has_cl_cmb_temperature, index_tt,1);
  class_define_index(ptr->index_tt_e, ppt->has_cl_cmb_polarization,index_tt,1);

  index_tt_common=index_tt;

  /** - type indices for scalars */

  if (ppt->has_scalars == _TRUE_) {

    index_tt = index_tt_common;

    class_define_index(ptr->index_tt_t0,     ppt->has_cl_cmb_temperature,      index_tt,1);
    class_define_index(ptr->index_tt_t1,     ppt->has_cl_cmb_temperature,      index_tt,1);
    class_define_index(ptr->index_tt_lcmb,   ppt->has_cl_cmb_lensing_potential,index_tt,1);
    class_define_index(ptr->index_tt_density,ppt->has_nc_density,              index_tt,ppt->selection_num);
    class_define_index(ptr->index_tt_rsd,    ppt->has_nc_rsd,                  index_tt,ppt->selection_num);
    class_define_index(ptr->index_tt_d0,     ppt->has_nc_rsd,                  index_tt,ppt->selection_num);
    class_define_index(ptr->index_tt_d1,     ppt->has_nc_rsd,                  index_tt,ppt->selection_num);
    class_define_index(ptr->index_tt_nc_lens,ppt->has_nc_lens,                 index_tt,ppt->selection_num);
    class_define_index(ptr->index_tt_nc_g1,  ppt->has_nc_gr,                   index_tt,ppt->selection_num);
    class_define_index(ptr->index_tt_nc_g2,  ppt->has_nc_gr,                   index_tt,ppt->selection_num);
    class_define_index(ptr->index_tt_nc_g3,  ppt->has_nc_gr,                   index_tt,ppt->selection_num);
    class_define_index(ptr->index_tt_nc_g4,  ppt->has_nc_gr,                   index_tt,ppt->selection_num);
    class_define_index(ptr->index_tt_nc_g5,  ppt->has_nc_gr,                   index_tt,ppt->selection_num);
    class_define_index(ptr->index_tt_lensing,ppt->has_cl_lensing_potential,    index_tt,ppt->selection_num);

    ptr->tt_size[ppt->index_md_scalars]=index_tt;

  }

  /** - type indices for vectors */

  if (ppt->has_vectors == _TRUE_) {

    index_tt = index_tt_common;

    /* For vectors we do not use index_tt_t1 and index_tt_b because this creates problems
       (as these are not indices which are common to the three types of modes, but only to two of them).
       Hence we define separate indices index_tt_vector_t1 and index_tt_vector_b */
    class_define_index(ptr->index_tt_vector_t1,ppt->has_cl_cmb_temperature, index_tt,1);
    class_define_index(ptr->index_tt_vector_b, ppt->has_cl_cmb_polarization,index_tt,1);

    ptr->tt_size[ppt->index_md_vectors]=index_tt;

  }

  /** - type indices for tensors */

  if (ppt->has_tensors == _TRUE_) {

    index_tt = index_tt_common;

    class_define_index(ptr->index_tt_b, ppt->has_cl_cmb_polarization,index_tt,1);

    ptr->tt_size[ppt->index_md_tensors]=index_tt;

  }

  /** - define the properties of each transfer type using transfer_type_properties() */

  class_call(transfer_type_properties(ppr,ppt,ptr),
             ptr->error_message,
             ptr->error_message);

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
  if (ptr->do_lcmb_full_limber == _TRUE_) {
    class_alloc(ptr->transfer_limber,ptr->md_size * sizeof(double *),ptr->error_message);
  }

  /** - get q values using transfer_get_q_list() */

  class_call(transfer_get_q_list(ppr,ppt,ptr,q_period,K,sgnK),
             ptr->error_message,
             ptr->error_message);

  /** - get q values of the full Limber q-sampling scheme using transfer_get_q_limber_list() */

  if (ptr->do_lcmb_full_limber == _TRUE_) {
    class_call(transfer_get_q_limber_list(ppr,ppt,ptr,K,sgnK),
               ptr->error_message,
               ptr->error_message);
  }
  else {
    ptr->q_size_limber=0;
  }

  /** - get k values using transfer_get_k_list() */

  class_call(transfer_get_k_list(ppt,ptr,K),
             ptr->error_message,
             ptr->error_message);

  /* for testing, it can be useful to print the q list in a file: */

  /*
    FILE * out=fopen("output/q","w");
    int index_q;

    for (index_q=0; index_q < ptr->q_size; index_q++) {

    fprintf(out,"%d %e %e %e %e\n",
    index_q,
    ptr->q[index_q],
    ptr->k[0][index_q],
    ptr->q[index_q]/sqrt(sgnK*K),
    ptr->q[index_q+1]-ptr->q[index_q]);

    }

    fclose(out);
  */

  /** - get l values using transfer_get_l_list() */
  class_call(transfer_get_l_list(ppr,ppt,ptr),
             ptr->error_message,
             ptr->error_message);

  /** - loop over modes (scalar, etc). For each mode: */

  for (index_md = 0; index_md < ptr->md_size; index_md++) {

    /** - allocate arrays of transfer functions, (ptr->transfer[index_md])[index_ic][index_tt][index_l][index_k] */
    class_alloc(ptr->transfer[index_md],
                ppt->ic_size[index_md] * ptr->tt_size[index_md] * ptr->l_size[index_md] * ptr->q_size * sizeof(double),
                ptr->error_message);

    if (ptr->do_lcmb_full_limber == _TRUE_) {
      class_alloc(ptr->transfer_limber[index_md],
                  ppt->ic_size[index_md] * ptr->tt_size[index_md] * ptr->l_size[index_md] * ptr->q_size_limber * sizeof(double),
                  ptr->error_message);
    }

  }

  return _SUCCESS_;

}

/**
 * This routine prepares the perturbation sources used by the transfer
 * module. If a non-linear method is used, the scalar sources of the
 * matter (or cdm+baryon) density and velocity and of the metric
 * perturbations are multiplied by the non-linear correction factor
 * sqrt(P_NL/P_L) of the matter (or cdm+baryon) power spectrum,
 * pfo->nl_corr_density, and stored in newly allocated arrays. All
 * other sources are not copied: the array only points to the
 * sources of the perturbation module.
 *
 * @param ppt     Input: pointer to perturbation structure
 * @param pfo     Input: pointer to fourier structure
 * @param ptr     Input: pointer to transfer structure
 * @param sources Output: array of sources, sources[index_md][index_ic*ppt->tp_size[index_md]+index_tp][index_tau*ppt->k_size[index_md]+index_k] (first level allocated by the caller, other levels allocated here)
 * @return the error status
 */

int transfer_perturbation_copy_sources_and_nl_corrections(
                                                          struct perturbations * ppt,
                                                          struct fourier * pfo,
                                                          struct transfer * ptr,
                                                          double *** sources
                                                          ) {
  int index_md;
  int index_ic;
  int index_tp;
  int index_k;
  int index_tau;

  for (index_md = 0; index_md < ptr->md_size; index_md++) {

    class_alloc(sources[index_md],
                ppt->ic_size[index_md]*ppt->tp_size[index_md]*sizeof(double*),
                ptr->error_message);

    for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++) {

      for (index_tp = 0; index_tp < ppt->tp_size[index_md]; index_tp++) {

        if ((pfo->method != nl_none) && (_scalars_) &&
            (((ppt->has_source_delta_m == _TRUE_) && (index_tp == ppt->index_tp_delta_m)) ||
             ((ppt->has_source_delta_cb == _TRUE_) && (index_tp == ppt->index_tp_delta_cb)) ||
             ((ppt->has_source_theta_m == _TRUE_) && (index_tp == ppt->index_tp_theta_m)) ||
             ((ppt->has_source_theta_cb == _TRUE_) && (index_tp == ppt->index_tp_theta_cb)) ||
             ((ppt->has_source_phi == _TRUE_) && (index_tp == ppt->index_tp_phi)) ||
             ((ppt->has_source_phi_prime == _TRUE_) && (index_tp == ppt->index_tp_phi_prime)) ||
             ((ppt->has_source_phi_plus_psi == _TRUE_) && (index_tp == ppt->index_tp_phi_plus_psi)) ||
             ((ppt->has_source_psi == _TRUE_) && (index_tp == ppt->index_tp_psi)))) {

          class_alloc(sources[index_md][index_ic * ppt->tp_size[index_md] + index_tp],
                      ppt->k_size[index_md]*ppt->tau_size*sizeof(double),
                      ptr->error_message);

          for (index_tau=0; index_tau<ppt->tau_size; index_tau++) {
            for (index_k=0; index_k<ppt->k_size[index_md]; index_k++) {
              if (((ppt->has_source_delta_cb == _TRUE_) && (index_tp == ppt->index_tp_delta_cb)) ||
                  ((ppt->has_source_theta_cb == _TRUE_) && (index_tp == ppt->index_tp_theta_cb))){
                sources[index_md]
                  [index_ic * ppt->tp_size[index_md] + index_tp]
                  [index_tau * ppt->k_size[index_md] + index_k] =
                  ppt->sources[index_md]
                  [index_ic * ppt->tp_size[index_md] + index_tp]
                  [index_tau * ppt->k_size[index_md] + index_k]
                  * pfo->nl_corr_density[pfo->index_pk_cb][index_tau * ppt->k_size[index_md] + index_k];
              }
              else{
                sources[index_md]
                  [index_ic * ppt->tp_size[index_md] + index_tp]
                  [index_tau * ppt->k_size[index_md] + index_k] =
                  ppt->sources[index_md]
                  [index_ic * ppt->tp_size[index_md] + index_tp]
                  [index_tau * ppt->k_size[index_md] + index_k]
                  * pfo->nl_corr_density[pfo->index_pk_m][index_tau * ppt->k_size[index_md] + index_k];
              }
            }
          }
        }
        else {
          sources[index_md][index_ic * ppt->tp_size[index_md] + index_tp] =
            ppt->sources[index_md][index_ic * ppt->tp_size[index_md] + index_tp];
        }
      }
    }
  }

  return _SUCCESS_;

}

/**
 * This routine computes the second derivatives of the perturbation
 * sources with respect to k, for spline interpolation of the sources
 * at the wavenumbers of the transfer module in
 * transfer_interpolate_sources().
 *
 * @param ppt            Input: pointer to perturbation structure
 * @param ptr            Input: pointer to transfer structure
 * @param sources        Input: array of sources, as filled by transfer_perturbation_copy_sources_and_nl_corrections()
 * @param sources_spline Output: array of second derivatives of the sources with respect to k, same indices (first level allocated by the caller, other levels allocated here)
 * @return the error status
 */

int transfer_perturbation_source_spline(
                                        struct perturbations * ppt,
                                        struct transfer * ptr,
                                        double *** sources,
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
                    ppt->k_size[index_md]*ppt->tau_size*sizeof(double),
                    ptr->error_message);

        class_call(array_spline_table_columns2(ppt->k[index_md],
                                               ppt->k_size[index_md],
                                               sources[index_md][index_ic * ppt->tp_size[index_md] + index_tp],
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

/**
 * This routine frees the array of perturbation sources allocated by
 * transfer_perturbation_copy_sources_and_nl_corrections(), i.e. only
 * the sources that were copied and corrected for non-linear effects
 * (the other ones belong to the perturbation module), and the array
 * itself.
 *
 * @param ppt     Input: pointer to perturbation structure
 * @param pfo     Input: pointer to fourier structure
 * @param ptr     Input: pointer to transfer structure
 * @param sources Input: array of sources to be freed
 * @return the error status
 */

int transfer_perturbation_sources_free(
                                       struct perturbations * ppt,
                                       struct fourier * pfo,
                                       struct transfer * ptr,
                                       double *** sources
                                       ) {
  int index_md;
  int index_ic;
  int index_tp;

  for (index_md = 0; index_md < ptr->md_size; index_md++) {
    for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++) {
      for (index_tp = 0; index_tp < ppt->tp_size[index_md]; index_tp++) {
        if ((pfo->method != nl_none) && (_scalars_) &&
            (((ppt->has_source_delta_m == _TRUE_) && (index_tp == ppt->index_tp_delta_m)) ||
             ((ppt->has_source_theta_m == _TRUE_) && (index_tp == ppt->index_tp_theta_m)) ||
             ((ppt->has_source_delta_cb == _TRUE_) && (index_tp == ppt->index_tp_delta_cb)) ||
             ((ppt->has_source_theta_cb == _TRUE_) && (index_tp == ppt->index_tp_theta_cb)) ||
             ((ppt->has_source_phi == _TRUE_) && (index_tp == ppt->index_tp_phi)) ||
             ((ppt->has_source_phi_prime == _TRUE_) && (index_tp == ppt->index_tp_phi_prime)) ||
             ((ppt->has_source_phi_plus_psi == _TRUE_) && (index_tp == ppt->index_tp_phi_plus_psi)) ||
             ((ppt->has_source_psi == _TRUE_) && (index_tp == ppt->index_tp_psi)))) {

          free(sources[index_md][index_ic * ppt->tp_size[index_md] + index_tp]);
        }
      }
    }
    free(sources[index_md]);
  }
  free(sources);

  return _SUCCESS_;
}

/**
 * This routine frees the array of second derivatives of the
 * perturbation sources allocated by
 * transfer_perturbation_source_spline().
 *
 * @param ppt            Input: pointer to perturbation structure
 * @param ptr            Input: pointer to transfer structure
 * @param sources_spline Input: array of second derivatives to be freed
 * @return the error status
 */

int transfer_perturbation_sources_spline_free(
                                              struct perturbations * ppt,
                                              struct transfer * ptr,
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
 * This routine defines the number and values of multipoles l for all modes.
 *
 * @param ppr  Input: pointer to precision structure
 * @param ppt  Input: pointer to perturbation structure
 * @param ptr  Input/Output: pointer to transfer structure containing l's
 * @return the error status
 */

int transfer_get_l_list(
                        struct precision * ppr,
                        struct perturbations * ppt,
                        struct transfer * ptr
                        ) {

  int index_l;
  int l_max=0;
  int index_md;
  int index_tt;
  int increment,current_l;
  /** Summary: */
  /*
    fprintf(stderr,"rescaling %e logstep %e linstep %e\n",
    ptr->angular_rescaling,
    pow(ppr->l_logstep,ptr->angular_rescaling),
    ppr->l_linstep*ptr->angular_rescaling);
  */

  /* check that largests need value of l_max */

  if (ppt->has_cls == _TRUE_) {

    if (ppt->has_scalars == _TRUE_) {

      if ((ppt->has_cl_cmb_temperature == _TRUE_) ||
          (ppt->has_cl_cmb_polarization == _TRUE_) ||
          (ppt->has_cl_cmb_lensing_potential == _TRUE_))
        l_max=MAX(ppt->l_scalar_max,l_max);

      if ((ppt->has_cl_lensing_potential == _TRUE_) ||
          (ppt->has_cl_number_count == _TRUE_))
        l_max=MAX(ppt->l_lss_max,l_max);
    }

    if (ppt->has_vectors == _TRUE_)
      l_max=MAX(ppt->l_vector_max,l_max);

    if (ppt->has_tensors == _TRUE_)
      l_max=MAX(ppt->l_tensor_max,l_max);

  }

  /** - allocate and fill l array */

  /** - start from l = 2 and increase with logarithmic step */

  index_l = 0;
  current_l = 2;
  increment = MAX((int)(current_l * (pow(ppr->l_logstep,ptr->angular_rescaling)-1.)),1);

  while (((current_l+increment) < l_max) &&
         (increment < ppr->l_linstep*ptr->angular_rescaling)) {

    index_l ++;
    current_l += increment;
    increment = MAX((int)(current_l * (pow(ppr->l_logstep,ptr->angular_rescaling)-1.)),1);

  }

  /** - when the logarithmic step becomes larger than some linear step,
      stick to this linear step till l_max */

  increment = ppr->l_linstep*ptr->angular_rescaling;

  while ((current_l+increment) <= l_max) {

    index_l ++;
    current_l += increment;

  }

  /** - last value set to exactly l_max */

  if (current_l != l_max) {

    index_l ++;
    current_l = l_max;

  }

  ptr->l_size_max = index_l+1;

  /** - so far we just counted the number of values. Now repeat the
      whole thing but fill array with values. */

  class_alloc(ptr->l,ptr->l_size_max*sizeof(int),ptr->error_message);

  index_l = 0;
  ptr->l[0] = 2;
  increment = MAX((int)(ptr->l[0] * (pow(ppr->l_logstep,ptr->angular_rescaling)-1.)),1);

  while (((ptr->l[index_l]+increment) < l_max) &&
         (increment < ppr->l_linstep*ptr->angular_rescaling)) {

    index_l ++;
    ptr->l[index_l]=ptr->l[index_l-1]+increment;
    increment = MAX((int)(ptr->l[index_l] * (pow(ppr->l_logstep,ptr->angular_rescaling)-1.)),1);

  }

  increment = ppr->l_linstep*ptr->angular_rescaling;

  while ((ptr->l[index_l]+increment) <= l_max) {

    index_l ++;
    ptr->l[index_l]=ptr->l[index_l-1]+increment;

  }

  if (ptr->l[index_l] != l_max) {

    index_l ++;
    ptr->l[index_l]= l_max;

  }

  /* for each mode and type, find relevant size of l array,
     l_size_tt[index_md][index_tt] (since for some modes and types
     l_max can be smaller). Also, maximize this size for each mode to
     find l_size[index_md]. */

  for (index_md=0; index_md < ppt->md_size; index_md++) {

    ptr->l_size[index_md] = 0;

    for (index_tt=0;index_tt<ptr->tt_size[index_md];index_tt++) {

      if (_scalars_) {

        /* CMB types (temperature, polarization, CMB lensing) */
        if (ptr->bin_of_tt[index_md][index_tt] < 0)
          l_max=ppt->l_scalar_max;

        /* types associated to a redshift bin (number count and galaxy lensing contributions) */
        else
          l_max=ppt->l_lss_max;

      }

      if (_vectors_) {
        l_max = ppt->l_vector_max;
      }

      if (_tensors_) {
        l_max = ppt->l_tensor_max;
      }

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

      ptr->l_size[index_md] = MAX(ptr->l_size[index_md],ptr->l_size_tt[index_md][index_tt]);

    }
  }

  return _SUCCESS_;

}

/**
 * This routine defines the number and values of wavenumbers q for
 * each mode (goes smoothly from logarithmic step for small q's to
 * linear step for large q's).
 *
 * @param ppr     Input: pointer to precision structure
 * @param ppt     Input: pointer to perturbation structure
 * @param ptr     Input/Output: pointer to transfer structure containing q's
 * @param q_period Input: order of magnitude of the oscillation period of transfer functions
 * @param K        Input: spatial curvature (in absolute value)
 * @param sgnK     Input: spatial curvature sign (open/closed/flat)
 * @return the error status
 */

int transfer_get_q_list(
                        struct precision * ppr,
                        struct perturbations * ppt,
                        struct transfer * ptr,
                        double q_period,
                        double K,
                        int sgnK
                        ) {

  int index_q;
  double q,q_min=0.,q_max=0.,q_step,k_max;
  int nu, nu_min, nu_proposed;
  int q_size_max;
  double q_approximation;
  double last_step=0.;
  int last_index=0;
  double q_logstep_spline;
  double q_logstep_trapzd;
  double q_threshold;
  int index_md;

  /* first and last value in flat case*/

  if (sgnK == 0) {
    q_min = ppt->k_min;

    q_max = 0.;
    for (index_md=0; index_md<ppt->md_size; index_md++) {
      q_max = MAX(q_max,ppt->k[index_md][ppt->k_size_cl[index_md]-1]);
    }

    K=0;
  }

  /* first and last value in open case*/

  else if (sgnK == -1) {
    q_min = sqrt(ppt->k_min*ppt->k_min+K);

    k_max = 0.;
    for (index_md=0; index_md<ppt->md_size; index_md++) {
      k_max = MAX(k_max,ppt->k[index_md][ppt->k_size_cl[index_md]-1]);
    }

    q_max = sqrt(k_max*k_max+K);
    if (ppt->has_vectors == _TRUE_)
      q_max = MIN(q_max,sqrt(k_max*k_max+2.*K));
    if (ppt->has_tensors == _TRUE_)
      q_max = MIN(q_max,sqrt(k_max*k_max+3.*K));
  }

  /* first and last value in closed case*/

  else if (sgnK == 1) {
    nu_min = 3;
    q_min = nu_min * sqrt(K);

    q_max = 0.;
    for (index_md=0; index_md<ppt->md_size; index_md++) {
      q_max = MAX(q_max,ppt->k[index_md][ppt->k_size_cl[index_md]-1]);
    }
  }

  /* adjust the parameter governing the log step size to curvature */

  q_logstep_spline = ppr->q_logstep_spline/pow(ptr->angular_rescaling,ppr->q_logstep_open);
  q_logstep_trapzd = ppr->q_logstep_trapzd;

  /* slightly conservative estimate of number of values */

  // One can divide the q range into 3 regions:
  // For q_threshold = q_linstep/q_logstep, we can show explicitly that
  // If q<q_threshold , then delta q > q * (1+q_logstep/2)
  // If q< N * q_threshold, then delta q > q + N/(N+1) * q_linstep
  // It makes sense to combine the first range, with the second range for N=1,10,100,...
  // In practice, we don't need percent level, so N=10 is sufficient
  // The total formula then reads
  // N_steps = ln(MIN(q_max,q_threshold)/q_min)/ln(1+b/2) + 18/q_logstep + 11/10 * (q_max-MIN(10 q_threshold, q_max))/(a)
  // The additional MIN/MAX just catch cases where q_max < q_threshold
  // The contribution from the term 18/q_logstep accounts for the range [q_threshold, 10 q_threshold], and in prinicple can be ignored when q_max < q_threshold, BUT
  // since usually q_logstep is of order unity, we neglect this detail, simplifying our lives
  // The +1 is just to always round up
  if (sgnK == 1) {

    q_approximation = MIN(ppr->hyper_flat_approximation_nu,(q_max/sqrt(K)))*sqrt(K);

    /* max contribution from integer nu values */
    q_threshold = ppr->q_linstep/q_logstep_trapzd;

    q_step = 1.+0.5*q_period*q_logstep_trapzd;
    q_size_max = (int)(log(MIN(q_approximation,q_threshold)/q_min)/log(q_step)+1);

    // Either linear q step OR AT LEAST delta nu = 1 until we reach the approximation/q_max
    q_step = MAX(q_period*ppr->q_linstep, 1*sqrt(K));
    q_size_max += (int)(1.1*(q_approximation-MIN(q_min,10*q_threshold))/q_step+18*q_threshold/q_step+1);

    /* max contribution from non-integer nu values */
    q_threshold = ppr->q_linstep/q_logstep_spline;

    q_step = 1.+0.5*q_period*q_logstep_spline;
    q_size_max += (int)(log(MIN(q_max, q_threshold)/q_approximation)/log(q_step)+1);

    q_step = q_period*ppr->q_linstep;
    q_size_max += (int)(1.1*(q_max-MIN(q_max, 10*q_threshold))/q_step+18*q_threshold/q_step+1);

    /* Make a final maximum with a very large number -- 2^25 (~3*10^7)
       If we have that many points, we have a problem anyways
       BUT maybe the integer steps in nu will 'save' us in this case
       In any case, let's try! */
    q_size_max = MIN(q_size_max, 1<<25);
  }
  else {

    /* max contribution from non-integer nu values */
    q_threshold = ppr->q_linstep/q_logstep_spline;

    q_step = 1.+0.5*q_period*q_logstep_spline;
    q_size_max = (int)(log(MIN(q_max, q_threshold)/q_min)/log(q_step)+1);

    q_step = q_period*ppr->q_linstep;
    q_size_max += (int)(1.1*(q_max-MIN(q_max, 10*q_threshold))/q_step+18*q_threshold/q_step+1);

  }

  /* create array with this conservative size estimate. The exact size
     will be readjusted below, after filling the array. */

  class_alloc(ptr->q,
              q_size_max*sizeof(double),
              ptr->error_message);

  /* assign the first value before starting the loop */

  index_q = 0;
  ptr->q[index_q] = q_min;
  nu = 3;
  index_q++;

  /* loop over the values */

  while (ptr->q[index_q-1] < q_max) {

    class_test(index_q >= q_size_max,ptr->error_message,"buggy q-list definition (q_size=%i)",q_size_max);

    /* step size formula in flat/open case. Step goes gradually from
       logarithmic to linear:

       - in the small q limit, it is logarithmic with: (delta q / q) =
       q_period * q_logstep_spline

       - in the large q limit, it is linear with: (delta q) = q_period
       * ppr->q_linstep
       */

    if (sgnK<=0) {

      q = ptr->q[index_q-1]
        + q_period * ppr->q_linstep * ptr->q[index_q-1]
        / (ptr->q[index_q-1] + ppr->q_linstep/q_logstep_spline);

    }

    /* step size formula in closed case. Same thing excepted that:

       - in the small q limit, the logarithmic step is reduced, being
       given by q_logstep_trapzd, and values are rounded to integer
       values of nu=q/sqrt(K). This happens as long as
       nu<nu_flat_approximation

       - for nu>nu_flat_approximation, the step gradually catches up
       the same expression as in the flat/open case, and there is no
       need to round up to integer nu's.
    */

    else {

      if (nu < (int)ppr->hyper_flat_approximation_nu) {

        q = ptr->q[index_q-1]
          + q_period * ppr->q_linstep * ptr->q[index_q-1]
          / (ptr->q[index_q-1] + ppr->q_linstep/q_logstep_trapzd);

        nu_proposed = (int)(q/sqrt(K));
        if (nu_proposed <= nu+1)
          nu = nu+1;
        else
          nu = nu_proposed;

        q = nu*sqrt(K);
        last_step = q - ptr->q[index_q-1];
        last_index = index_q+1;
      }
      else {

        q_step = q_period * ppr->q_linstep * ptr->q[index_q-1] / (ptr->q[index_q-1] + ppr->q_linstep/q_logstep_spline);

        if (index_q-last_index < (int)ppr->q_numstep_transition)
          q = ptr->q[index_q-1] + (1-(double)(index_q-last_index)/ppr->q_numstep_transition) * last_step + (double)(index_q-last_index)/ppr->q_numstep_transition * q_step;
        else
          q = ptr->q[index_q-1] + q_step;
      }
    }

    ptr->q[index_q] = q;
    index_q++;
  }

  /* infer total number of values (also checking if we overshot the last point) */

  if (ptr->q[index_q-1] > q_max)
    ptr->q_size=index_q-1;
  else
    ptr->q_size=index_q;

  class_test(ptr->q_size<2,ptr->error_message,"buggy q-list definition");

  /* now, readjust array size */

  class_realloc(ptr->q,
                ptr->q_size*sizeof(double),
                ptr->error_message);

  /* in curved universe, check at which index the flat rescaling
     approximation will start being used */

  if (sgnK != 0) {

    q_approximation = ppr->hyper_flat_approximation_nu * sqrt(sgnK*K);
    /* If q_approximation >= q_max, we can skip any computation and conclude we cannot do an approximation */
    if(q_approximation < ptr->q[ptr->q_size-1]){
      /* Otherwise, just find the first index at which the approximation would break */
      for (ptr->index_q_flat_approximation=0;
           ptr->index_q_flat_approximation < ptr->q_size-1;
           ptr->index_q_flat_approximation++) {
        if (ptr->q[ptr->index_q_flat_approximation] > q_approximation) break;
      }
      if (ptr->transfer_verbose > 1)
        printf("Flat bessel approximation spares hyperspherical bessel computations for %zu wavenumbers over a total of %zu\n",
               ptr->q_size-ptr->index_q_flat_approximation,ptr->q_size);
    }
    else{
      ptr->index_q_flat_approximation = ptr->q_size;
      if (ptr->transfer_verbose > 1)
        printf("Cannot use flat bessel approximation, since q_max = %.5e > q_start_approximation = %.5e\n", ptr->q[ptr->q_size-1], q_approximation);
    }
  }

  return _SUCCESS_;

}

/**
 * This routine defines the number and values of wavenumbers q_limber for
 * each mode (logarithmic step only).
 *
 * @param ppr     Input: pointer to precision structure
 * @param ppt     Input: pointer to perturbation structure
 * @param ptr     Input/Output: pointer to transfer structure containing q's
 * @param K        Input: spatial curvature (in absolute value)
 * @param sgnK     Input: spatial curvature sign (open/closed/flat)
 * @return the error status
 */

int transfer_get_q_limber_list(
                               struct precision * ppr,
                               struct perturbations * ppt,
                               struct transfer * ptr,
                               double K,
                               int sgnK
                               ) {

  int index_q;
  double q_min=0.,q_max=0.,k_max;
  int nu_min;
  int index_md;

  /* first and last value in flat case*/

  if (sgnK == 0) {
    q_min = ppt->k_min;

    q_max = 0.;
    for (index_md=0; index_md<ppt->md_size; index_md++) {
      q_max = MAX(q_max,ppt->k[index_md][ppt->k_size[index_md]-1]);
    }

    K=0;
  }

  /* first and last value in open case*/

  else if (sgnK == -1) {
    q_min = sqrt(ppt->k_min*ppt->k_min+K);

    k_max = 0.;
    for (index_md=0; index_md<ppt->md_size; index_md++) {
      k_max = MAX(k_max,ppt->k[index_md][ppt->k_size[index_md]-1]);
    }

    q_max = sqrt(k_max*k_max+K);
    if (ppt->has_vectors == _TRUE_)
      q_max = MIN(q_max,sqrt(k_max*k_max+2.*K));
    if (ppt->has_tensors == _TRUE_)
      q_max = MIN(q_max,sqrt(k_max*k_max+3.*K));
  }

  /* first and last value in closed case*/

  else if (sgnK == 1) {
    nu_min = 3;
    q_min = nu_min * sqrt(K);

    q_max = 0.;
    for (index_md=0; index_md<ppt->md_size; index_md++) {
      q_max = MAX(q_max,ppt->k[index_md][ppt->k_size[index_md]-1]);
    }
  }

  /* number of values */
  ptr->q_size_limber = (int)(log(q_max/q_min)/log(ppr->q_logstep_limber))+1;

  class_alloc(ptr->q_limber,
              ptr->q_size_limber*sizeof(double),
              ptr->error_message);

  /* assign the first value before starting the loop */

  ptr->q_limber[0] = q_min;
  for (index_q = 1; index_q < ptr->q_size_limber; index_q++) {
    ptr->q_limber[index_q] = ppr->q_logstep_limber *  ptr->q_limber[index_q-1];
  }

  return _SUCCESS_;

}

/**
 * This routine infers from the q values a list of corresponding k
 * values for each mode.
 *
 * @param ppt     Input: pointer to perturbation structure
 * @param ptr     Input/Output: pointer to transfer structure containing q's
 * @param K       Input: spatial curvature
 * @return the error status
 */

int transfer_get_k_list(
                        struct perturbations * ppt,
                        struct transfer * ptr,
                        double K
                        ) {

  int index_md;
  int index_q;
  double m=0.;

  class_alloc(ptr->k,ptr->md_size*sizeof(double*),ptr->error_message);
  if (ptr->do_lcmb_full_limber == _TRUE_) {
    class_alloc(ptr->k_limber,ptr->md_size*sizeof(double*),ptr->error_message);
  }

  for (index_md = 0; index_md <  ptr->md_size; index_md++) {

    class_alloc(ptr->k[index_md],ptr->q_size*sizeof(double),ptr->error_message);
    if (ptr->do_lcmb_full_limber == _TRUE_) {
      class_alloc(ptr->k_limber[index_md],ptr->q_size_limber*sizeof(double),ptr->error_message);
    }

    if (_scalars_) {
      m=0.;
    }
    if (_vectors_) {
      m=1.;
    }
    if (_tensors_) {
      m=2.;
    }

    for (index_q=0; index_q < ptr->q_size; index_q++) {
      ptr->k[index_md][index_q] = sqrt(ptr->q[index_q]*ptr->q[index_q]-K*(m+1.));
    }
    if (ptr->do_lcmb_full_limber == _TRUE_) {
      for (index_q=0; index_q < ptr->q_size_limber; index_q++) {
        ptr->k_limber[index_md][index_q] = sqrt(ptr->q_limber[index_q]*ptr->q_limber[index_q]-K*(m+1.));
      }
    }

    /* check consistency of the first value of ptr->k */
    if (ptr->k[index_md][0] < ppt->k[index_md][0]){
      /* If ptr->k[index_md][0] < ppt->k[index_md][0] at the level of rounding,
         adjust first value of k_list to avoid interpolation errors: */
      if ((ppt->k[index_md][0]-ptr->k[index_md][0]) < 10.*DBL_EPSILON){
        ptr->k[index_md][0] = ppt->k[index_md][0];
      }
      else{
        class_stop(ptr->error_message,
                   "bug in k_list calculation: in perturbation module k_min=%e, in transfer module k_min[mode=%d]=%e, interpolation impossible",
                   ppt->k[0][0],
                   index_md,
                   ptr->k[index_md][0]);
      }
    }

    /* check consistency of the first value of ptr->k_limber */
    if (ptr->do_lcmb_full_limber == _TRUE_) {
      if (ptr->k_limber[index_md][0] < ppt->k[index_md][0]){
        /* If ptr->k_limber[index_md][0] < ppt->k[index_md][0] at the level of rounding,
           adjust first value of k_list to avoid interpolation errors: */
        if ((ppt->k[index_md][0]-ptr->k_limber[index_md][0]) < 10.*DBL_EPSILON){
          ptr->k_limber[index_md][0] = ppt->k[index_md][0];
        }
        else{
          class_stop(ptr->error_message,
                     "bug in k_list calculation: in perturbation module k_min=%e, in transfer module k_min[mode=%d]=%e, interpolation impossible",
                     ppt->k[0][0],
                     index_md,
                     ptr->k_limber[index_md][0]);
        }
      }
    }

    /* check consistency of the last value of ptr->k compared to the one of ppt->k across all modes (hence the 0 index) */
    class_test(ptr->k[index_md][ptr->q_size-1] > ppt->k[0][ppt->k_size_cl[0]-1],
               ptr->error_message,
               "bug in k_list calculation: in perturbation module k_max=%e, in transfer module k_max[mode=%d]=%e, interpolation impossible",
               ppt->k[0][ppt->k_size_cl[0]-1],
               index_md,
               ptr->k[index_md][ptr->q_size-1]);

    /* check consistency of the last value of ptr->k_limber compared to the one of ppt->k across all modes (hence the 0 index) */
    if (ptr->do_lcmb_full_limber == _TRUE_) {
      class_test(ptr->k_limber[index_md][ptr->q_size_limber-1] > ppt->k[0][ppt->k_size[0]-1],
                 ptr->error_message,
                 "bug in k_list calculation: in perturbation module k_max=%e, in transfer module k_max[mode=%d]=%e, interpolation impossible",
                 ppt->k[0][ppt->k_size[0]-1],
                 index_md,
                 ptr->k_limber[index_md][ptr->q_size_limber-1]);
    }

  }

  return _SUCCESS_;

}

/**
 * This routine defines once for all the properties of each transfer
 * type X (index_tt) for each mode: the perturbation source Y(X)
 * (index_tp) entering its line-of-sight integral, the redshift bin
 * (for number count and galaxy lensing types), the radial function
 * convolved with the source, the type of window function W(tau)
 * multiplying the source, a possible factor depending on the
 * wavenumber, and the conditions under which the transfer function is
 * computed with the Limber approximation or neglected.
 *
 * All the distinctions between transfer types are made here, so that
 * the rest of the module can treat all types in the same way.
 *
 * @param ppr  Input: pointer to precision structure
 * @param ppt  Input: pointer to perturbation structure
 * @param ptr  Input/Output: pointer to transfer structure (indices index_tt_... defined before, properties filled here)
 * @return the error status
 */

int transfer_type_properties(
                             struct precision * ppr,
                             struct perturbations * ppt,
                             struct transfer * ptr
                             ) {

  /** Summary: */

  /** - define local variables */

  /* running index on modes */
  int index_md;

  /* running index on transfer types */
  int index_tt;

  /* redshift bin */
  int bin;

  /* is the type a number count contribution? (for those the Limber approximation is never used with Dirac selection functions) */
  short is_number_count;

  /** - allocate arrays of properties */

  class_alloc(ptr->tp_of_tt,ptr->md_size*sizeof(int*),ptr->error_message);
  class_alloc(ptr->bin_of_tt,ptr->md_size*sizeof(int*),ptr->error_message);
  class_alloc(ptr->radial_type,ptr->md_size*sizeof(radial_function_type*),ptr->error_message);
  class_alloc(ptr->window_type,ptr->md_size*sizeof(enum transfer_window_type*),ptr->error_message);
  class_alloc(ptr->k_power,ptr->md_size*sizeof(double*),ptr->error_message);
  class_alloc(ptr->k_pivot,ptr->md_size*sizeof(double*),ptr->error_message);
  class_alloc(ptr->l_limber,ptr->md_size*sizeof(double*),ptr->error_message);
  class_alloc(ptr->delta_k_neglect,ptr->md_size*sizeof(double*),ptr->error_message);
  class_alloc(ptr->has_late_source_cut,ptr->md_size*sizeof(short*),ptr->error_message);
  class_alloc(ptr->has_full_limber,ptr->md_size*sizeof(short*),ptr->error_message);

  for (index_md = 0; index_md < ptr->md_size; index_md++) {

    class_alloc(ptr->tp_of_tt[index_md],ptr->tt_size[index_md]*sizeof(int),ptr->error_message);
    class_alloc(ptr->bin_of_tt[index_md],ptr->tt_size[index_md]*sizeof(int),ptr->error_message);
    class_alloc(ptr->radial_type[index_md],ptr->tt_size[index_md]*sizeof(radial_function_type),ptr->error_message);
    class_alloc(ptr->window_type[index_md],ptr->tt_size[index_md]*sizeof(enum transfer_window_type),ptr->error_message);
    class_alloc(ptr->k_power[index_md],ptr->tt_size[index_md]*sizeof(double),ptr->error_message);
    class_alloc(ptr->k_pivot[index_md],ptr->tt_size[index_md]*sizeof(double),ptr->error_message);
    class_alloc(ptr->l_limber[index_md],ptr->tt_size[index_md]*sizeof(double),ptr->error_message);
    class_alloc(ptr->delta_k_neglect[index_md],ptr->tt_size[index_md]*sizeof(double),ptr->error_message);
    class_alloc(ptr->has_late_source_cut[index_md],ptr->tt_size[index_md]*sizeof(short),ptr->error_message);
    class_alloc(ptr->has_full_limber[index_md],ptr->tt_size[index_md]*sizeof(short),ptr->error_message);

    for (index_tt = 0; index_tt < ptr->tt_size[index_md]; index_tt++) {

      /** - start from default properties: W(tau)=1, no factor
          depending on k, no Limber approximation, never neglected */

      ptr->tp_of_tt[index_md][index_tt] = -1;
      ptr->bin_of_tt[index_md][index_tt] = -1;
      ptr->radial_type[index_md][index_tt] = SCALAR_TEMPERATURE_0;
      ptr->window_type[index_md][index_tt] = window_none;
      ptr->k_power[index_md][index_tt] = 0.;
      ptr->k_pivot[index_md][index_tt] = 1.;
      ptr->l_limber[index_md][index_tt] = _HUGE_;
      ptr->delta_k_neglect[index_md][index_tt] = _HUGE_;
      ptr->has_late_source_cut[index_md][index_tt] = _FALSE_;
      ptr->has_full_limber[index_md][index_tt] = _FALSE_;
      is_number_count = _FALSE_;

      /** - scalar modes */

      if (_scalars_) {

        /** - --> CMB temperature and polarisation: the transfer function
            can be neglected for l < (q-delta_k_neglect) r_a(tau_rec);
            sources at late times can be neglected at high l,
            excepted for T0, which contains the late ISW effect */

        if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t0)) {
          ptr->tp_of_tt[index_md][index_tt] = ppt->index_tp_t0;
          ptr->radial_type[index_md][index_tt] = SCALAR_TEMPERATURE_0;
          ptr->delta_k_neglect[index_md][index_tt] = ppr->transfer_neglect_delta_k_S_t0;
        }

        if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t1)) {
          ptr->tp_of_tt[index_md][index_tt] = ppt->index_tp_t1;
          ptr->radial_type[index_md][index_tt] = SCALAR_TEMPERATURE_1;
          ptr->delta_k_neglect[index_md][index_tt] = ppr->transfer_neglect_delta_k_S_t1;
          ptr->has_late_source_cut[index_md][index_tt] = _TRUE_;
        }

        if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t2)) {
          ptr->tp_of_tt[index_md][index_tt] = ppt->index_tp_t2;
          ptr->radial_type[index_md][index_tt] = SCALAR_TEMPERATURE_2;
          ptr->delta_k_neglect[index_md][index_tt] = ppr->transfer_neglect_delta_k_S_t2;
          ptr->has_late_source_cut[index_md][index_tt] = _TRUE_;
        }

        if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_e)) {
          ptr->tp_of_tt[index_md][index_tt] = ppt->index_tp_p;
          ptr->radial_type[index_md][index_tt] = SCALAR_POLARISATION_E;
          ptr->delta_k_neglect[index_md][index_tt] = ppr->transfer_neglect_delta_k_S_e;
          ptr->has_late_source_cut[index_md][index_tt] = _TRUE_;
        }

        /** - --> CMB lensing potential: lensing kernel with support
            [tau_rec, tau0], possibly rescaled by hand by a factor
            lcmb_rescale (k/lcmb_pivot)^lcmb_tilt. The Limber
            approximation is used for l > l_switch_limber. This is
            the only type computed also in the full Limber q-sampling scheme. */

        if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (index_tt == ptr->index_tt_lcmb)) {
          ptr->tp_of_tt[index_md][index_tt] = ppt->index_tp_phi_plus_psi;
          ptr->window_type[index_md][index_tt] = window_cmb_lensing;
          ptr->k_power[index_md][index_tt] = ptr->lcmb_tilt;
          ptr->k_pivot[index_md][index_tt] = ptr->lcmb_pivot;
          ptr->l_limber[index_md][index_tt] = floor(ppr->l_switch_limber)+1.;
          ptr->has_full_limber[index_md][index_tt] = ptr->do_lcmb_full_limber;
        }

        /** - --> number count contributions which are local in
            each redshift bin (see CLASSgal paper arXiv:1307.1459):
            density, redshift-space distortion, Doppler terms d0 and
            d1, gravity terms g1, g2, g3. Their window is the selection
            function of the bin, times a function of the
            background. The Limber approximation is used for l >=
            l_switch_limber_for_nc_local_over_z * z_mean */

        if (_index_tt_in_range_(ptr->index_tt_density, ppt->selection_num, ppt->has_nc_density)) {
          ptr->bin_of_tt[index_md][index_tt] = index_tt - ptr->index_tt_density;
          /* use here delta_cb rather than delta_m if density number counts calculated only for cold dark matter + baryon */
          /* (this important comment is referenced in a WARNING message in perturbations.c) */
          ptr->tp_of_tt[index_md][index_tt] = ppt->index_tp_delta_m;
          ptr->window_type[index_md][index_tt] = window_selection;
          ptr->l_limber[index_md][index_tt] = ppr->l_switch_limber_for_nc_local_over_z;
          is_number_count = _TRUE_;
        }

        if (_index_tt_in_range_(ptr->index_tt_rsd, ppt->selection_num, ppt->has_nc_rsd)) {
          ptr->bin_of_tt[index_md][index_tt] = index_tt - ptr->index_tt_rsd;
          ptr->tp_of_tt[index_md][index_tt] = ppt->index_tp_theta_m;
          ptr->radial_type[index_md][index_tt] = NC_RSD;
          ptr->window_type[index_md][index_tt] = window_selection;
          ptr->l_limber[index_md][index_tt] = ppr->l_switch_limber_for_nc_local_over_z;
          is_number_count = _TRUE_;
        }

        if (_index_tt_in_range_(ptr->index_tt_d0, ppt->selection_num, ppt->has_nc_rsd)) {
          ptr->bin_of_tt[index_md][index_tt] = index_tt - ptr->index_tt_d0;
          ptr->tp_of_tt[index_md][index_tt] = ppt->index_tp_theta_m;
          ptr->window_type[index_md][index_tt] = window_selection;
          ptr->k_power[index_md][index_tt] = -2.; // Factor from original ClassGAL paper ( arXiv 1307.1459 )
          ptr->l_limber[index_md][index_tt] = ppr->l_switch_limber_for_nc_local_over_z;
          is_number_count = _TRUE_;
        }

        if (_index_tt_in_range_(ptr->index_tt_d1, ppt->selection_num, ppt->has_nc_rsd)) {
          ptr->bin_of_tt[index_md][index_tt] = index_tt - ptr->index_tt_d1;
          ptr->tp_of_tt[index_md][index_tt] = ppt->index_tp_theta_m;
          ptr->radial_type[index_md][index_tt] = SCALAR_TEMPERATURE_1;
          ptr->window_type[index_md][index_tt] = window_selection;
          ptr->k_power[index_md][index_tt] = -1.; // Factor from original ClassGAL paper ( arXiv 1307.1459 )
          ptr->l_limber[index_md][index_tt] = ppr->l_switch_limber_for_nc_local_over_z;
          is_number_count = _TRUE_;
        }

        if (_index_tt_in_range_(ptr->index_tt_nc_g1, ppt->selection_num, ppt->has_nc_gr)) {
          ptr->bin_of_tt[index_md][index_tt] = index_tt - ptr->index_tt_nc_g1;
          ptr->tp_of_tt[index_md][index_tt] = ppt->index_tp_psi;
          ptr->window_type[index_md][index_tt] = window_selection;
          ptr->l_limber[index_md][index_tt] = ppr->l_switch_limber_for_nc_local_over_z;
          is_number_count = _TRUE_;
        }

        if (_index_tt_in_range_(ptr->index_tt_nc_g2, ppt->selection_num, ppt->has_nc_gr)) {
          ptr->bin_of_tt[index_md][index_tt] = index_tt - ptr->index_tt_nc_g2;
          ptr->tp_of_tt[index_md][index_tt] = ppt->index_tp_phi;
          ptr->window_type[index_md][index_tt] = window_selection;
          ptr->l_limber[index_md][index_tt] = ppr->l_switch_limber_for_nc_local_over_z;
          is_number_count = _TRUE_;
        }

        if (_index_tt_in_range_(ptr->index_tt_nc_g3, ppt->selection_num, ppt->has_nc_gr)) {
          ptr->bin_of_tt[index_md][index_tt] = index_tt - ptr->index_tt_nc_g3;
          ptr->tp_of_tt[index_md][index_tt] = ppt->index_tp_phi_prime;
          ptr->window_type[index_md][index_tt] = window_selection;
          ptr->l_limber[index_md][index_tt] = ppr->l_switch_limber_for_nc_local_over_z;
          is_number_count = _TRUE_;
        }

        /** - --> number count contributions integrated along the
            line of sight in each bin (lensing, gravity terms g4 and
            g5) and galaxy lensing potential in each bin. Their window
            is the selection function of the bin convolved with a
            kernel, with support [tau_min, tau0]. The Limber
            approximation is used for l >=
            l_switch_limber_for_nc_los_over_z * z_mean (and for l >=
            l_switch_limber_for_nc_local_over_z * z_mean for g5,
            which involves the derivative of the Bessel function) */

        if (_index_tt_in_range_(ptr->index_tt_nc_lens, ppt->selection_num, ppt->has_nc_lens)) {
          ptr->bin_of_tt[index_md][index_tt] = index_tt - ptr->index_tt_nc_lens;
          ptr->tp_of_tt[index_md][index_tt] = ppt->index_tp_phi_plus_psi;
          ptr->window_type[index_md][index_tt] = window_integrated;
          ptr->l_limber[index_md][index_tt] = ppr->l_switch_limber_for_nc_los_over_z;
          is_number_count = _TRUE_;
        }

        if (_index_tt_in_range_(ptr->index_tt_nc_g4, ppt->selection_num, ppt->has_nc_gr)) {
          ptr->bin_of_tt[index_md][index_tt] = index_tt - ptr->index_tt_nc_g4;
          ptr->tp_of_tt[index_md][index_tt] = ppt->index_tp_phi_plus_psi;
          ptr->window_type[index_md][index_tt] = window_integrated;
          ptr->l_limber[index_md][index_tt] = ppr->l_switch_limber_for_nc_los_over_z;
          is_number_count = _TRUE_;
        }

        if (_index_tt_in_range_(ptr->index_tt_nc_g5, ppt->selection_num, ppt->has_nc_gr)) {
          ptr->bin_of_tt[index_md][index_tt] = index_tt - ptr->index_tt_nc_g5;
          ptr->tp_of_tt[index_md][index_tt] = ppt->index_tp_phi_plus_psi;
          ptr->radial_type[index_md][index_tt] = SCALAR_TEMPERATURE_1;
          ptr->window_type[index_md][index_tt] = window_integrated;
          ptr->k_power[index_md][index_tt] = 1.; // Factor from chi derivative of d/dchi j_ell(k*chi)= d/d(kchi) j_ell(k chi) * k = k * j_ell'(kchi)
          ptr->l_limber[index_md][index_tt] = ppr->l_switch_limber_for_nc_local_over_z;
          is_number_count = _TRUE_;
        }

        if (_index_tt_in_range_(ptr->index_tt_lensing, ppt->selection_num, ppt->has_cl_lensing_potential)) {
          ptr->bin_of_tt[index_md][index_tt] = index_tt - ptr->index_tt_lensing;
          ptr->tp_of_tt[index_md][index_tt] = ppt->index_tp_phi_plus_psi;
          ptr->window_type[index_md][index_tt] = window_integrated;
          ptr->l_limber[index_md][index_tt] = ppr->l_switch_limber_for_nc_los_over_z;
        }

        /** - --> for all types associated to a redshift bin, the
            Limber threshold is proportional to the mean redshift of
            the bin. For number counts with Dirac selection functions,
            the Limber approximation is never used. */

        bin = ptr->bin_of_tt[index_md][index_tt];

        if (bin >= 0) {
          ptr->l_limber[index_md][index_tt] *= ppt->selection_mean[bin];
          if ((ppt->selection == dirac) && (is_number_count == _TRUE_))
            ptr->l_limber[index_md][index_tt] = _HUGE_;
        }
      }

      /** - vector modes */

      if (_vectors_) {

        if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_vector_t1)) {
          ptr->tp_of_tt[index_md][index_tt] = ppt->index_tp_vector_t1;
          ptr->radial_type[index_md][index_tt] = VECTOR_TEMPERATURE_1;
          ptr->delta_k_neglect[index_md][index_tt] = ppr->transfer_neglect_delta_k_V_t1;
          ptr->has_late_source_cut[index_md][index_tt] = _TRUE_;
        }

        if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t2)) {
          ptr->tp_of_tt[index_md][index_tt] = ppt->index_tp_t2;
          ptr->radial_type[index_md][index_tt] = VECTOR_TEMPERATURE_2;
          ptr->delta_k_neglect[index_md][index_tt] = ppr->transfer_neglect_delta_k_V_t2;
          ptr->has_late_source_cut[index_md][index_tt] = _TRUE_;
        }

        if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_e)) {
          ptr->tp_of_tt[index_md][index_tt] = ppt->index_tp_p;
          ptr->radial_type[index_md][index_tt] = VECTOR_POLARISATION_E;
          ptr->delta_k_neglect[index_md][index_tt] = ppr->transfer_neglect_delta_k_V_e;
          ptr->has_late_source_cut[index_md][index_tt] = _TRUE_;
        }

        if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_vector_b)) {
          ptr->tp_of_tt[index_md][index_tt] = ppt->index_tp_p;
          ptr->radial_type[index_md][index_tt] = VECTOR_POLARISATION_B;
          ptr->delta_k_neglect[index_md][index_tt] = ppr->transfer_neglect_delta_k_V_b;
          ptr->has_late_source_cut[index_md][index_tt] = _TRUE_;
        }
      }

      /** - tensor modes (the late time sources of temperature can never be neglected) */

      if (_tensors_) {

        if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t2)) {
          ptr->tp_of_tt[index_md][index_tt] = ppt->index_tp_t2;
          ptr->radial_type[index_md][index_tt] = TENSOR_TEMPERATURE_2;
          ptr->delta_k_neglect[index_md][index_tt] = ppr->transfer_neglect_delta_k_T_t2;
        }

        if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_e)) {
          ptr->tp_of_tt[index_md][index_tt] = ppt->index_tp_p;
          ptr->radial_type[index_md][index_tt] = TENSOR_POLARISATION_E;
          ptr->delta_k_neglect[index_md][index_tt] = ppr->transfer_neglect_delta_k_T_e;
          ptr->has_late_source_cut[index_md][index_tt] = _TRUE_;
        }

        if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_b)) {
          ptr->tp_of_tt[index_md][index_tt] = ppt->index_tp_p;
          ptr->radial_type[index_md][index_tt] = TENSOR_POLARISATION_B;
          ptr->delta_k_neglect[index_md][index_tt] = ppr->transfer_neglect_delta_k_T_b;
          ptr->has_late_source_cut[index_md][index_tt] = _TRUE_;
        }
      }

      class_test(ptr->tp_of_tt[index_md][index_tt] < 0,
                 ptr->error_message,
                 "no perturbation source associated to transfer type %d of mode %d",
                 index_tt,index_md);
    }
  }

  return _SUCCESS_;

}

/**
 * This routine defines the time sampling of the line-of-sight
 * integral of a given transfer type, that is, a type of sampling, an
 * interval [tau_min, tau_max] and a number of points.
 *
 * The time sampling must catch both the time dependence of the
 * transfer source W(tau) S(k, tau) and the oscillations of the Bessel
 * functions. For CMB types, the time sampling of the perturbation
 * module is used. For number count and galaxy lensing types, the
 * transfer source is sampled with evenly spaced times. The number of
 * points is then at least ppr->selection_sampling, and large enough
 * to sample the oscillations of the Bessel function \f$ j_l(k (\tau_0-\tau)) \f$ (with
 * period \f$ \Delta \tau = 2\pi/k \f$) for the largest relevant wavenumber.
 *
 * @param ppr           Input: pointer to precision structure
 * @param ppt           Input: pointer to perturbation structure
 * @param ptr           Input: pointer to transfer structure
 * @param tau_rec       Input: recombination time
 * @param tau0          Input: time today
 * @param index_md      Input: index of the mode
 * @param index_tt      Input: index of transfer type
 * @param sampling_type Output: type of time sampling
 * @param tau_min       Output: smallest sampled time (for sampling_from_perturbations and sampling_logarithmic: sampled times are strictly larger than tau_min)
 * @param tau_max       Output: largest sampled time (for sampling_logarithmic: largest time sampled with logarithmic steps, before tau0)
 * @param log_step      Output: for sampling_logarithmic, relative step Delta tau / (tau0-tau)
 * @param tau_size      Output: number of sampled times
 * @return the error status
 */

int transfer_type_sampling(
                           struct precision * ppr,
                           struct perturbations * ppt,
                           struct transfer * ptr,
                           double tau_rec,
                           double tau0,
                           int index_md,
                           int index_tt,
                           enum transfer_sampling_type * sampling_type,
                           double * tau_min,
                           double * tau_max,
                           double * log_step,
                           int * tau_size
                           ) {

  /* redshift bin of the type */
  int bin = ptr->bin_of_tt[index_md][index_tt];

  /* conformal time at the mean redshift of the bin */
  double tau_mean;

  /* value of l at which limber approximation is switched on */
  int l_limber;

  *log_step = 0.;

  switch (ptr->window_type[index_md][index_tt]) {

    /** - CMB temperature and polarisation: all times of the perturbation module */

  case window_none:

    *sampling_type = sampling_from_perturbations;
    *tau_min = 0.;
    *tau_max = tau0;
    *tau_size = ppt->tau_size;
    break;

    /** - CMB lensing potential: the times of the perturbation module
        after recombination, as long as they are separated by less
        than Delta tau = (tau0-tau)/l_max*ppr->cmb_lensing_sampling_factor;
        then this logarithmic time step, until tau = (1-ppr->cmb_lensing_sampling_fraction)*tau0;
        then tau0 */

  case window_cmb_lensing:

    *sampling_type = sampling_logarithmic;
    *tau_min = tau_rec;
    *tau_max = (1.-ppr->cmb_lensing_sampling_fraction)*tau0;
    *log_step = 1./MIN(ppr->l_switch_limber,ppt->l_scalar_max)*ppr->cmb_lensing_sampling_factor;

    class_test(ppr->cmb_lensing_sampling_fraction*(*log_step) < 10.*ppr->smallest_allowed_variation,
               ptr->error_message,
               "cmb_lensing_sampling_fraction=%e is too small\n and/or cmb_lensing_sampling_factor=%e is too small\n and/or l_switch_limber=%g is too big\n and/or l_scalar_max=%d is too big,\n such that the logarithmic step used to sample the CMB lensing source function\n would become of the same order of magnitude as machine precision...",
               ppr->cmb_lensing_sampling_fraction,
               ppr->cmb_lensing_sampling_factor,
               ppr->l_switch_limber,
               ppt->l_scalar_max);

    /* count the number of sampled times */
    class_call(transfer_logarithmic_sampling(ppt,
                                             tau0,
                                             *tau_min,
                                             *tau_max,
                                             *log_step,
                                             NULL,
                                             tau_size,
                                             ptr->error_message),
               ptr->error_message,
               ptr->error_message);
    break;

    /** - local number count contributions: the support of the
        selection function of the bin, with one point for a Dirac
        selection function */

  case window_selection:

    *sampling_type = sampling_uniform;
    *tau_min = ptr->selection_tau_min[bin];
    *tau_max = ptr->selection_tau_max[bin];
    tau_mean = ptr->selection_tau_mean[bin];

    /* case selection=dirac */
    if (*tau_min == *tau_max) {
      *tau_size = 1;
    }
    /* other cases (gaussian, top-hat...) */
    else {

      /* check that selection function well sampled */
      *tau_size = (int)ppr->selection_sampling;

      /* value of l at which the code switches to Limber approximation
         (necessary for next step) */
      l_limber=ppr->l_switch_limber_for_nc_local_over_z*ppt->selection_mean[bin];

      /* check that bessel well sampled, if not define finer sampling
         overwriting the previous one.
         One Bessel oscillations corresponds to [Delta tau]=2pi/k.
         This is minimal for largest relevant k_max,
         namely k_max=l_limber/(tau0-tau_mean).
         We need to cut the interval (tau_max-tau_min) in pieces of size
         [Delta tau]=2pi/k_max. This gives the number below.
      */
      *tau_size=MAX(*tau_size,(int)((*tau_max-*tau_min)/((tau0-tau_mean)/MIN(l_limber,ppt->l_lss_max)))*ppr->selection_sampling_bessel);
    }
    break;

    /** - number count contributions integrated along the line of
        sight and galaxy lensing: from the beginning of the selection
        function of the bin till today */

  case window_integrated:

    *sampling_type = sampling_uniform;
    *tau_min = ptr->selection_tau_min[bin];
    *tau_max = tau0;
    tau_mean = ptr->selection_tau_mean[bin];

    /* check that selection function well sampled */
    *tau_size = (int)ppr->selection_sampling;

    /* value of l at which the code switches to Limber approximation
       (necessary for next step) */
    if (_index_tt_in_range_(ptr->index_tt_nc_g5, ppt->selection_num, ppt->has_nc_gr)) {
      /* Even if G5 is integrated along the line-of-sight, we do not apply the same Limber criteria as for the other integrated terms, because here we have the derivative of the Bessel.  */
      l_limber=ppr->l_switch_limber_for_nc_local_over_z*ppt->selection_mean[bin];
      *tau_size=MAX(*tau_size,(int)((tau0-*tau_min)/((tau0-tau_mean)/2./MIN(l_limber,ppt->l_lss_max)))*ppr->selection_sampling_bessel);
    }
    else {
      l_limber=ppr->l_switch_limber_for_nc_los_over_z*ppt->selection_mean[bin];
      /* check that bessel well sampled, if not define finer sampling
         overwriting the previous one.
         One Bessel oscillations corresponds to [Delta tau]=2pi/k.
         This is minimal for largest relevant k_max,
         namely k_max=l_limber/((tau0-tau_mean)/2).
         We need to cut the interval (tau_0-tau_min) in pieces of size
         [Delta tau]=2pi/k_max. This gives the number below.
      */
      *tau_size=MAX(*tau_size,(int)((tau0-*tau_min)/((tau0-tau_mean)/2./MIN(l_limber,ppt->l_lss_max)))*ppr->selection_sampling_bessel_los);
    }
    break;
  }

  return _SUCCESS_;
}

/**
 * This routine fills tau_size evenly spaced values of (tau0-tau)
 * between (tau0-tau_min) and (tau0-tau_max), in decreasing order
 * (just (tau0-tau_min) if tau_size=1)
 *
 * @param tau0           Input: time today
 * @param tau_min        Input: smallest time
 * @param tau_max        Input: largest time
 * @param tau_size       Input: number of values
 * @param tau0_minus_tau Output: values of (tau0-tau) (allocated before)
 * @return the error status
 */

int transfer_uniform_sampling(
                              double tau0,
                              double tau_min,
                              double tau_max,
                              int tau_size,
                              double * tau0_minus_tau
                              ) {

  int index_tau;

  if (tau_size == 1) {
    tau0_minus_tau[0] = tau0-tau_min;
  }
  else {
    for (index_tau=0; index_tau<tau_size-1; index_tau++) {
      tau0_minus_tau[index_tau]=tau0-tau_min-((double)index_tau)/((double)tau_size-1.)*(tau_max-tau_min);
    }
    tau0_minus_tau[tau_size-1]=tau0-tau_max;
  }

  return _SUCCESS_;
}

/**
 * This routine defines a time sampling with logarithmic steps at late
 * times: it keeps the times of the perturbation module larger than
 * tau_min, as long as they are separated by less than Delta tau =
 * (tau0-tau)*log_step; then it uses this logarithmic step until
 * tau=tau_max; and it adds a last value tau=tau0. If tau0_minus_tau
 * is NULL, the routine only counts the number of values.
 *
 * @param ppt            Input: pointer to perturbation structure
 * @param tau0           Input: time today
 * @param tau_min        Input: sampled times are larger than tau_min
 * @param tau_max        Input: largest time sampled with logarithmic steps
 * @param log_step       Input: logarithmic step Delta tau/(tau0-tau)
 * @param tau0_minus_tau Output: values of (tau0-tau) in decreasing order (allocated before), or NULL
 * @param tau_size       Output: number of values
 * @param error_message  Output: error message
 * @return the error status
 */

int transfer_logarithmic_sampling(
                                  struct perturbations * ppt,
                                  double tau0,
                                  double tau_min,
                                  double tau_max,
                                  double log_step,
                                  double * tau0_minus_tau,
                                  int * tau_size,
                                  ErrorMsg error_message
                                  ) {

  /* index of time in ppt->tau_sampling array */
  int index_tau_sampling=0;

  /* index of time in the new time sampling */
  int index_tau=0;

  /* running value of conformal time, and last value of (tau0-tau) */
  double tau, tau0_minus_tau_last;

  /* throw away times before tau_min */
  while (ppt->tau_sampling[index_tau_sampling] <= tau_min) index_tau_sampling++;

  class_test((index_tau_sampling < 1) || (index_tau_sampling >= ppt->tau_size),
             error_message,
             "index_tau_sampling=%d out of bounds (ppt->tau_size=%d)",
             index_tau_sampling,
             ppt->tau_size);

  /* use the time steps from the perturbation module as
     long as they are smaller than (tau0-tau)*log_step */
  while ((index_tau_sampling < ppt->tau_size) &&
         ((ppt->tau_sampling[index_tau_sampling]-ppt->tau_sampling[index_tau_sampling-1]) < (tau0-ppt->tau_sampling[index_tau_sampling-1])*log_step)) {
    if (tau0_minus_tau != NULL)
      tau0_minus_tau[index_tau] = tau0 - ppt->tau_sampling[index_tau_sampling];
    index_tau++;
    index_tau_sampling++;
  }

  /* if we have not reached tau0 already, use the logarithmic step until tau = tau_max */
  if (index_tau_sampling < ppt->tau_size) {

    tau = ppt->tau_sampling[index_tau_sampling-1];
    tau0_minus_tau_last = tau0 - tau;

    while (tau < tau_max) {
      tau += (tau0-tau)*log_step;
      tau0_minus_tau_last = tau0 - tau;
      if (tau0_minus_tau != NULL)
        tau0_minus_tau[index_tau] = tau0_minus_tau_last;
      index_tau++;
    }

    /* add the last value given by tau = tau0 ... */
    if (tau0_minus_tau_last > 0.) {
      if (tau0_minus_tau != NULL)
        tau0_minus_tau[index_tau] = 0.;
      index_tau++;
    }
    /* ... unless we overpassed it. Then, we must go back */
    else {
      if (tau0_minus_tau != NULL)
        tau0_minus_tau[index_tau-1] = 0.;
    }
  }

  *tau_size = index_tau;

  return _SUCCESS_;
}

/**
 * This routine defines once for all the time sampling of the
 * line-of-sight integral of each transfer type, using
 * transfer_type_sampling(). Types with an identical time sampling
 * share the same arrays. For each time sampling, it stores the values
 * of (tau0-tau), the weights of the trapezoidal integration, and the
 * indices and weights needed for resampling the perturbation sources
 * by linear interpolation.
 *
 * @param ppr     Input: pointer to precision structure
 * @param pba     Input: pointer to background structure
 * @param ppt     Input: pointer to perturbation structure
 * @param ptr     Input/Output: pointer to transfer structure
 * @param tau_rec Input: recombination time
 * @return the error status
 */

int transfer_time_samplings(
                            struct precision * ppr,
                            struct background * pba,
                            struct perturbations * ppt,
                            struct transfer * ptr,
                            double tau_rec
                            ) {

  /** Summary: */

  /** - define local variables */

  int index_md;
  int index_tt;
  int index_tt_other;
  int index_ts;
  int index_tau;

  /* conformal time today, and a running value of conformal time */
  double tau0 = pba->conformal_age;
  double tau;

  /* description of the time sampling of each type */
  enum transfer_sampling_type * sampling_type;
  double * tau_min;
  double * tau_max;
  double * log_step;
  int * tau_size;

  /* index of the time in ppt->tau_sampling just before a sampled time */
  int index_tau_inf;

  /* number of sampled times, for checking consistency */
  int tau_size_check;

  /* sampled values for a given time sampling */
  double * tau0_minus_tau;
  int * index_tau_pt;
  double * weight_tau_pt;

  class_test(ppr->selection_sampling_bessel_min <= 0.,
             ptr->error_message,
             "selection_sampling_bessel_min=%e should be strictly positive (it is a minimum number of sampled points per radian of the Bessel functions)",
             ppr->selection_sampling_bessel_min);

  /** - allocate arrays */

  class_alloc(ptr->ts_size,ptr->md_size*sizeof(int),ptr->error_message);
  class_alloc(ptr->ts_of_tt,ptr->md_size*sizeof(int*),ptr->error_message);
  class_alloc(ptr->tau_size_ts,ptr->md_size*sizeof(int*),ptr->error_message);
  class_alloc(ptr->tau0_minus_tau,ptr->md_size*sizeof(double**),ptr->error_message);
  class_alloc(ptr->w_trapz,ptr->md_size*sizeof(double**),ptr->error_message);
  class_alloc(ptr->index_tau_pt,ptr->md_size*sizeof(int**),ptr->error_message);
  class_alloc(ptr->weight_tau_pt,ptr->md_size*sizeof(double**),ptr->error_message);
  class_alloc(ptr->q_max_sampled,ptr->md_size*sizeof(double*),ptr->error_message);

  ptr->tau_size_max = 0;

  for (index_md = 0; index_md < ptr->md_size; index_md++) {

    class_alloc(ptr->ts_of_tt[index_md],ptr->tt_size[index_md]*sizeof(int),ptr->error_message);
    class_alloc(sampling_type,ptr->tt_size[index_md]*sizeof(enum transfer_sampling_type),ptr->error_message);
    class_alloc(tau_min,ptr->tt_size[index_md]*sizeof(double),ptr->error_message);
    class_alloc(tau_max,ptr->tt_size[index_md]*sizeof(double),ptr->error_message);
    class_alloc(log_step,ptr->tt_size[index_md]*sizeof(double),ptr->error_message);
    class_alloc(tau_size,ptr->tt_size[index_md]*sizeof(int),ptr->error_message);

    /** - get the time sampling of each type, and give it a new
        index index_ts unless it is identical to the one of a
        previous type */

    ptr->ts_size[index_md] = 0;

    for (index_tt = 0; index_tt < ptr->tt_size[index_md]; index_tt++) {

      class_call(transfer_type_sampling(ppr,
                                        ppt,
                                        ptr,
                                        tau_rec,
                                        tau0,
                                        index_md,
                                        index_tt,
                                        &(sampling_type[index_tt]),
                                        &(tau_min[index_tt]),
                                        &(tau_max[index_tt]),
                                        &(log_step[index_tt]),
                                        &(tau_size[index_tt])),
                 ptr->error_message,
                 ptr->error_message);

      ptr->ts_of_tt[index_md][index_tt] = -1;

      for (index_tt_other = 0; index_tt_other < index_tt; index_tt_other++) {
        if ((sampling_type[index_tt_other] == sampling_type[index_tt]) &&
            (tau_min[index_tt_other] == tau_min[index_tt]) &&
            (tau_max[index_tt_other] == tau_max[index_tt]) &&
            (log_step[index_tt_other] == log_step[index_tt]) &&
            (tau_size[index_tt_other] == tau_size[index_tt])) {
          ptr->ts_of_tt[index_md][index_tt] = ptr->ts_of_tt[index_md][index_tt_other];
          break;
        }
      }

      if (ptr->ts_of_tt[index_md][index_tt] == -1) {
        ptr->ts_of_tt[index_md][index_tt] = ptr->ts_size[index_md];
        ptr->ts_size[index_md]++;
      }
    }

    class_alloc(ptr->tau_size_ts[index_md],ptr->ts_size[index_md]*sizeof(int),ptr->error_message);
    class_alloc(ptr->tau0_minus_tau[index_md],ptr->ts_size[index_md]*sizeof(double*),ptr->error_message);
    class_alloc(ptr->w_trapz[index_md],ptr->ts_size[index_md]*sizeof(double*),ptr->error_message);
    class_alloc(ptr->index_tau_pt[index_md],ptr->ts_size[index_md]*sizeof(int*),ptr->error_message);
    class_alloc(ptr->weight_tau_pt[index_md],ptr->ts_size[index_md]*sizeof(double*),ptr->error_message);
    class_alloc(ptr->q_max_sampled[index_md],ptr->ts_size[index_md]*sizeof(double),ptr->error_message);

    /** - fill each time sampling (using the description of the first type that uses it) */

    for (index_tt = 0; index_tt < ptr->tt_size[index_md]; index_tt++) {

      index_ts = ptr->ts_of_tt[index_md][index_tt];

      /* skip if this time sampling has been filled already */
      for (index_tt_other = 0; index_tt_other < index_tt; index_tt_other++) {
        if (ptr->ts_of_tt[index_md][index_tt_other] == index_ts)
          break;
      }
      if (index_tt_other < index_tt)
        continue;

      ptr->tau_size_ts[index_md][index_ts] = tau_size[index_tt];
      ptr->tau_size_max = MAX(ptr->tau_size_max,tau_size[index_tt]);

      class_alloc(ptr->tau0_minus_tau[index_md][index_ts],tau_size[index_tt]*sizeof(double),ptr->error_message);
      class_alloc(ptr->w_trapz[index_md][index_ts],tau_size[index_tt]*sizeof(double),ptr->error_message);
      class_alloc(ptr->index_tau_pt[index_md][index_ts],tau_size[index_tt]*sizeof(int),ptr->error_message);
      class_alloc(ptr->weight_tau_pt[index_md][index_ts],tau_size[index_tt]*sizeof(double),ptr->error_message);

      tau0_minus_tau = ptr->tau0_minus_tau[index_md][index_ts];
      index_tau_pt = ptr->index_tau_pt[index_md][index_ts];
      weight_tau_pt = ptr->weight_tau_pt[index_md][index_ts];

      switch (sampling_type[index_tt]) {

        /** - --> subset of the times of the perturbation module: the
            sources do not need to be resampled, the resampling
            weights are zero. (For the last time, we use the
            equivalent left index ppt->tau_size-2 and weight one, so
            that the next index always exists) */

      case sampling_from_perturbations:

        for (index_tau = 0; index_tau < tau_size[index_tt]; index_tau++) {

          index_tau_pt[index_tau] = ppt->tau_size-tau_size[index_tt]+index_tau;
          weight_tau_pt[index_tau] = 0.;
          tau0_minus_tau[index_tau] = tau0 - ppt->tau_sampling[index_tau_pt[index_tau]];

          if (index_tau_pt[index_tau] == ppt->tau_size-1) {
            index_tau_pt[index_tau] = ppt->tau_size-2;
            weight_tau_pt[index_tau] = 1.;
          }
        }

        ptr->q_max_sampled[index_md][index_ts] = _HUGE_;
        break;

        /** - --> evenly spaced times, or times with logarithmic
            steps at late times: find for each of them the bracketing
            times in the perturbation module by bisection, and the
            corresponding weight for linear interpolation */

      case sampling_uniform:
      case sampling_logarithmic:

        if (sampling_type[index_tt] == sampling_uniform) {
          class_call(transfer_uniform_sampling(tau0,
                                               tau_min[index_tt],
                                               tau_max[index_tt],
                                               tau_size[index_tt],
                                               tau0_minus_tau),
                     ptr->error_message,
                     ptr->error_message);
        }
        else {
          class_call(transfer_logarithmic_sampling(ppt,
                                                   tau0,
                                                   tau_min[index_tt],
                                                   tau_max[index_tt],
                                                   log_step[index_tt],
                                                   tau0_minus_tau,
                                                   &tau_size_check,
                                                   ptr->error_message),
                     ptr->error_message,
                     ptr->error_message);

          class_test(tau_size_check != tau_size[index_tt],
                     ptr->error_message,
                     "inconsistent number of sampled times: %d instead of %d",
                     tau_size_check,tau_size[index_tt]);
        }

        for (index_tau = 0; index_tau < tau_size[index_tt]; index_tau++) {

          tau = tau0 - tau0_minus_tau[index_tau];

          class_test((tau < ppt->tau_sampling[0]) || (tau > ppt->tau_sampling[ppt->tau_size-1]),
                     ptr->error_message,
                     "transfer source needed at tau=%e, outside the range [%e, %e] in which perturbation sources are sampled",
                     tau,ppt->tau_sampling[0],ppt->tau_sampling[ppt->tau_size-1]);

          /* find index_tau_inf such that ppt->tau_sampling[index_tau_inf] <= tau < ppt->tau_sampling[index_tau_inf+1]
             (or index_tau_inf = ppt->tau_size-2 if tau is the last time) */
          class_call(array_search_bisect(ppt->tau_size,
                                         ppt->tau_sampling,
                                         tau,
                                         &index_tau_inf,
                                         ptr->error_message),
                     ptr->error_message,
                     ptr->error_message);

          index_tau_pt[index_tau] = index_tau_inf;
          weight_tau_pt[index_tau] = (tau-ppt->tau_sampling[index_tau_inf])/(ppt->tau_sampling[index_tau_inf+1]-ppt->tau_sampling[index_tau_inf]);
        }

        /* for evenly spaced times: largest wavenumber for which the
           Bessel functions j_l(q(tau0-tau)) are sampled with at least
           ppr->selection_sampling_bessel_min points per radian (no
           limit for a single time, or for the logarithmic sampling,
           designed to resolve the Bessel functions down to tau0).
           transfer_use_limber() switches to the Limber approximation
           for q > q_max_sampled: the value _HUGE_ means that this
           never happens. */
        if ((sampling_type[index_tt] == sampling_uniform) && (tau_size[index_tt] > 1)) {
          ptr->q_max_sampled[index_md][index_ts] = (tau_size[index_tt]-1.)/(tau_max[index_tt]-tau_min[index_tt])/ppr->selection_sampling_bessel_min;
        }
        else {
          ptr->q_max_sampled[index_md][index_ts] = _HUGE_;
        }
        break;
      }

      /** - --> trapezoidal integration weights */

      class_call(array_trapezoidal_mweights(tau0_minus_tau,
                                            tau_size[index_tt],
                                            ptr->w_trapz[index_md][index_ts],
                                            ptr->error_message),
                 ptr->error_message,
                 ptr->error_message);
    }

    free(sampling_type);
    free(tau_min);
    free(tau_max);
    free(log_step);
    free(tau_size);
  }

  return _SUCCESS_;
}

/**
 * This routine computes the transfer functions \f$ \Delta_l^{X} (q) \f$
 * for a given wavenumber q, for all modes, initial conditions, types
 * and multipoles l.
 *
 * For each mode, it first interpolates all the needed perturbation
 * sources at the corresponding k(q). Then, for each time sampling, it
 * builds the transfer sources of all types sharing this sampling with
 * transfer_sources(). For each l, it chooses the method for computing
 * the transfer function of each of these types and initial conditions
 * with transfer_select_method(), interpolates the needed Bessel
 * functions once for all of them with transfer_bessel_interpolate(),
 * and computes the transfer functions with transfer_limber() or
 * transfer_integrate().
 *
 * @param ppr                 Input: pointer to precision structure
 * @param pba                 Input: pointer to background structure
 * @param ppt                 Input: pointer to perturbation structure
 * @param ptr                 Input/Output: pointer to transfer structure (result stored there)
 * @param index_q             Input: index of wavenumber in ptr->q or ptr->q_limber
 * @param tau_rec             Input: recombination time
 * @param pert_sources        Input: perturbation sources
 * @param pert_sources_spline Input: second derivative of perturbation sources with respect to k
 * @param ptw                 Input: pointer to workspace
 * @param use_full_limber     Input: whether to use the full Limber q-sampling scheme (then index_q refers to ptr->q_limber)
 * @return the error status
 */

int transfer_compute_for_each_q(
                                struct precision * ppr,
                                struct background * pba,
                                struct perturbations * ppt,
                                struct transfer * ptr,
                                int index_q,
                                double tau_rec,
                                double *** pert_sources,
                                double *** pert_sources_spline,
                                struct transfer_workspace * ptw,
                                short use_full_limber
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
  /* running index for time sampling schemes */
  int index_ts;
  /* index of a perturbation source */
  int index_tp;
  /* index running over the transfer types (tt) that share a common time sampling scheme (ts) */
  int index_tt_of_ts;
  /* index in the flattened array of initial conditions (ic) and of transfer types (tt) that share a common time sampling scheme (ts) */
  int index_ic_tt_of_ts;

  /* current time sampling */
  int tau_size;
  double * tau0_minus_tau;
  double * w_trapz;

  /* transfer source of the current type */
  double * sources;

  /* table where results are stored (ptr->transfer or ptr->transfer_limber) and its number of wavenumbers */
  double * transfer;
  size_t q_size;

  /** - for a given l, maximum value of k such that we can convolve
      the source with Bessel functions j_l(x) without reaching x_max */
  double q_max_bessel=0.;

  double l;

  double q,k,k_max;

  double transfer_function;

  /* minimum value of (tau0-tau) at which the Bessel function is non-negligible, for the current l */
  double tau0_minus_tau_min_bessel=0.;

  /* number of times at which the Bessel functions are needed for the current l */
  int bessel_size;

  /* is the current time sampling used by a type computed in the full Limber q-sampling scheme? */
  short has_full_limber;

  /* which Bessel functions and derivatives are needed for the current l, and for a given type */
  short need_Phi, need_dPhi, need_d2Phi;
  short need_Phi_tt, need_dPhi_tt, need_d2Phi_tt;

  /** - loop over all modes. For each mode */

  for (index_md = 0; index_md < ptr->md_size; index_md++) {

    if (use_full_limber == _FALSE_) {
      q = ptr->q[index_q];
      k = ptr->k[index_md][index_q];
      k_max = ppt->k[index_md][ppt->k_size_cl[index_md]-1];
      transfer = ptr->transfer[index_md];
      q_size = ptr->q_size;
    }
    else {
      q = ptr->q_limber[index_q];
      k = ptr->k_limber[index_md][index_q];
      k_max = ppt->k[index_md][ppt->k_size[index_md]-1];
      transfer = ptr->transfer_limber[index_md];
      q_size = ptr->q_size_limber;
    }

    /** - if we reached q_max for this mode, there is nothing to be done */

    if (k > k_max) {
      for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++) {
        for (index_tt = 0; index_tt < ptr->tt_size[index_md]; index_tt++) {
          for (index_l = 0; index_l < ptr->l_size[index_md]; index_l++) {
            transfer[((index_ic * ptr->tt_size[index_md] + index_tt) * ptr->l_size[index_md] + index_l) * q_size + index_q] = 0.;
          }
        }
      }
      continue;
    }

    /** - interpolate at the current value of k all the perturbation
        sources needed by this mode, for each initial condition */

    for (index_ic = 0; index_ic < ptw->ic_size_max; index_ic++)
      for (index_tp = 0; index_tp < ptw->tp_size_max; index_tp++)
        ptw->is_interpolated[index_ic*ptw->tp_size_max+index_tp] = _FALSE_;

    for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++) {
      for (index_tt = 0; index_tt < ptr->tt_size[index_md]; index_tt++) {

        /* in the full Limber q-sampling scheme, only the types with has_full_limber are computed */
        if ((use_full_limber == _TRUE_) && (ptr->has_full_limber[index_md][index_tt] == _FALSE_))
          continue;

        index_tp = ptr->tp_of_tt[index_md][index_tt];

        if (ptw->is_interpolated[index_ic*ptw->tp_size_max+index_tp] == _FALSE_) {

          if (ptw->interpolated_sources[index_ic*ptw->tp_size_max+index_tp] == NULL) {
            class_alloc(ptw->interpolated_sources[index_ic*ptw->tp_size_max+index_tp],
                        ppt->tau_size*sizeof(double),
                        ptr->error_message);
          }

          class_call(transfer_interpolate_sources(ppt,
                                                  ptr,
                                                  k,
                                                  index_md,
                                                  index_ic,
                                                  index_tp,
                                                  pert_sources[index_md][index_ic * ppt->tp_size[index_md] + index_tp],
                                                  pert_sources_spline[index_md][index_ic * ppt->tp_size[index_md] + index_tp],
                                                  ptw->interpolated_sources[index_ic*ptw->tp_size_max+index_tp]),
                     ptr->error_message,
                     ptr->error_message);

          ptw->is_interpolated[index_ic*ptw->tp_size_max+index_tp] = _TRUE_;
        }
      }
    }

    /** - loop over time sampling schemes */

    for (index_ts = 0; index_ts < ptr->ts_size[index_md]; index_ts++) {

      /** - --> if we are in the full Limber q-sampling scheme
          (use_full_limber == _TRUE_), only the types for which the
          Limber approximation must be computed on this q-grid
          (has_full_limber == _TRUE_, i.e. the CMB lensing potential
          of scalar modes) are computed; the transfer functions of
          the other types are set to zero. Hence, if no type using
          the current time sampling scheme has has_full_limber ==
          _TRUE_, we set the transfer functions of all its types to
          zero and skip to the next time sampling scheme */

      has_full_limber = _FALSE_;
      for (index_tt = 0; index_tt < ptr->tt_size[index_md]; index_tt++) {
        if ((ptr->ts_of_tt[index_md][index_tt] == index_ts) && (ptr->has_full_limber[index_md][index_tt] == _TRUE_))
          has_full_limber = _TRUE_;
      }

      if ((use_full_limber == _TRUE_) && (has_full_limber == _FALSE_)) {
        for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++) {
          for (index_tt = 0; index_tt < ptr->tt_size[index_md]; index_tt++) {
            if (ptr->ts_of_tt[index_md][index_tt] != index_ts)
              continue;
            for (index_l = 0; index_l < ptr->l_size[index_md]; index_l++) {
              transfer[((index_ic * ptr->tt_size[index_md] + index_tt) * ptr->l_size[index_md] + index_l) * q_size + index_q] = 0.;
            }
          }
        }
        continue;
      }

      tau_size = ptr->tau_size_ts[index_md][index_ts];
      tau0_minus_tau = ptr->tau0_minus_tau[index_md][index_ts];
      w_trapz = ptr->w_trapz[index_md][index_ts];

      /** - --> infer the array of radial coordinates r(tau0_minus_tau)
          as well as a few other quantities related by trigonometric
          functions */

      class_call(transfer_radial_coordinates(ptw,tau0_minus_tau,tau_size,k),
                 ptr->error_message,
                 ptr->error_message);

      /* for a given l, maximum value of k such that we can
         convolve the source with Bessel functions j_l(x)
         without reaching x_max (this is relevant in the flat
         case when the bessels are computed with the old bessel
         module. otherwise this condition is guaranteed by the
         choice of proper xmax when computing bessels) */
      if (ptw->sgnK == 0) {
        q_max_bessel = ptw->pBIS->x[ptw->pBIS->x_size-1]/tau0_minus_tau[0];
      }
      else {
        q_max_bessel = ptr->q[ptr->q_size-1];
      }

      /** - --> compute the transfer sources W(tau) S(k,tau) (k/k_pivot)^k_power
          of all types using this time sampling, for each initial
          condition (in the full Limber q-sampling scheme, only those of the
          types with has_full_limber are needed) */

      for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++) {
        index_tt_of_ts = 0;
        for (index_tt = 0; index_tt < ptr->tt_size[index_md]; index_tt++) {
          if (ptr->ts_of_tt[index_md][index_tt] != index_ts)
            continue;
          if ((use_full_limber == _FALSE_) || (ptr->has_full_limber[index_md][index_tt] == _TRUE_)) {
            class_call(transfer_sources(ptr,
                                        ptw->interpolated_sources[index_ic*ptw->tp_size_max+ptr->tp_of_tt[index_md][index_tt]],
                                        k,
                                        index_md,
                                        index_tt,
                                        ptw->sources + (index_ic * ptw->tt_per_ts_max + index_tt_of_ts) * ptr->tau_size_max),
                       ptr->error_message,
                       ptr->error_message);
          }
          index_tt_of_ts++;
        }
      }

      /** - --> loop over l */

      for (index_l = 0; index_l < ptr->l_size[index_md]; index_l++) {

        l = (double)ptr->l[index_l];

        /** - ---> find the minimum value of (tau0-tau) at which the
            Bessel function is non-negligible (not needed in the full
            Limber q-sampling scheme) */

        if (use_full_limber == _FALSE_) {
          class_call(transfer_bessel_cut(ptr,
                                         ptw,
                                         index_q,
                                         index_l,
                                         l,
                                         k,
                                         &tau0_minus_tau_min_bessel),
                     ptr->error_message,
                     ptr->error_message);
        }

        /** - ---> for each initial condition and type using this time
            sampling, choose the method for computing the transfer
            function, and find which Bessel functions (and derivatives)
            are needed, and at how many times */

        bessel_size = 0;
        need_Phi = _FALSE_;
        need_dPhi = _FALSE_;
        need_d2Phi = _FALSE_;

        for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++) {
          index_tt_of_ts = 0;
          for (index_tt = 0; index_tt < ptr->tt_size[index_md]; index_tt++) {
            if (ptr->ts_of_tt[index_md][index_tt] != index_ts)
              continue;

            index_ic_tt_of_ts = index_ic * ptw->tt_per_ts_max + index_tt_of_ts;

            class_call(transfer_select_method(ppr,
                                              ptr,
                                              ptw,
                                              index_q,
                                              index_md,
                                              index_tt,
                                              index_l,
                                              l,
                                              q,
                                              (pba->conformal_age-tau_rec)*ptr->angular_rescaling,
                                              q_max_bessel,
                                              tau0_minus_tau_min_bessel,
                                              tau0_minus_tau,
                                              tau_size,
                                              ptw->sources + index_ic_tt_of_ts * ptr->tau_size_max,
                                              use_full_limber,
                                              &(ptw->method[index_ic_tt_of_ts]),
                                              &(ptw->index_tau_max[index_ic_tt_of_ts]),
                                              &(ptw->bessel_truncation[index_ic_tt_of_ts])),
                       ptr->error_message,
                       ptr->error_message);

            if (ptw->method[index_ic_tt_of_ts] == transfer_line_of_sight) {
              bessel_size = MAX(bessel_size,ptw->index_tau_max[index_ic_tt_of_ts]+1);
              class_call(transfer_radial_function_needs(ptr->radial_type[index_md][index_tt],
                                                        &need_Phi_tt,
                                                        &need_dPhi_tt,
                                                        &need_d2Phi_tt),
                         ptr->error_message,
                         ptr->error_message);
              need_Phi = need_Phi || need_Phi_tt;
              need_dPhi = need_dPhi || need_dPhi_tt;
              need_d2Phi = need_d2Phi || need_d2Phi_tt;
            }
            index_tt_of_ts++;
          }
        }

        /** - ---> interpolate these Bessel functions once for all types and initial conditions */

        if (bessel_size > 0) {
          class_call(transfer_bessel_interpolate(ptw,
                                                 ptr,
                                                 index_q,
                                                 index_l,
                                                 bessel_size,
                                                 need_Phi,
                                                 need_dPhi,
                                                 need_d2Phi),
                     ptr->error_message,
                     ptr->error_message);
        }

        /** - ---> compute and store the transfer functions */

        for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++) {
          index_tt_of_ts = 0;
          for (index_tt = 0; index_tt < ptr->tt_size[index_md]; index_tt++) {
            if (ptr->ts_of_tt[index_md][index_tt] != index_ts)
              continue;

            index_ic_tt_of_ts = index_ic * ptw->tt_per_ts_max + index_tt_of_ts;
            sources = ptw->sources + index_ic_tt_of_ts * ptr->tau_size_max;

            switch (ptw->method[index_ic_tt_of_ts]) {

            case transfer_set_zero:
              transfer_function = 0.;
              break;

            case transfer_limber_approximation:
              class_call(transfer_limber(ptr,
                                         ptw,
                                         l,
                                         q,
                                         tau0_minus_tau,
                                         tau_size,
                                         sources,
                                         ptr->radial_type[index_md][index_tt],
                                         &transfer_function),
                         ptr->error_message,
                         ptr->error_message);
              break;

            case transfer_line_of_sight:
              class_call(transfer_integrate(ptr,
                                            ptw,
                                            index_md,
                                            index_tt,
                                            index_l,
                                            k,
                                            tau0_minus_tau,
                                            w_trapz,
                                            sources,
                                            ptw->index_tau_max[index_ic_tt_of_ts],
                                            ptw->bessel_truncation[index_ic_tt_of_ts],
                                            tau0_minus_tau_min_bessel,
                                            &transfer_function),
                         ptr->error_message,
                         ptr->error_message);
              break;
            }

            transfer[((index_ic * ptr->tt_size[index_md] + index_tt) * ptr->l_size[index_md] + index_l) * q_size + index_q] = transfer_function;

            index_tt_of_ts++;
          }
        }

      } /* end of loop over l */

    } /* end of loop over time sampling schemes */

  } /* end of loop over modes */

  return _SUCCESS_;

}

/**
 * This routine computes the radial coordinates (argument of the
 * Bessel functions) corresponding to a time sampling, and some
 * trigonometric functions of them.
 *
 * @param ptw            Input/Output: pointer to workspace (results stored there)
 * @param tau0_minus_tau Input: sampled values of (tau0-tau)
 * @param tau_size       Input: number of sampled values
 * @param k              Input: wavenumber
 * @return the error status
 */

int transfer_radial_coordinates(
                                struct transfer_workspace * ptw,
                                double * tau0_minus_tau,
                                int tau_size,
                                double k
                                ) {

  int index_tau;
  double sqrt_absK=0.;

  switch (ptw->sgnK){
  case 1:
    sqrt_absK = sqrt(ptw->K);
    for (index_tau=0; index_tau < tau_size; index_tau++) {
      ptw->chi[index_tau] = sqrt_absK*tau0_minus_tau[index_tau];
      ptw->cscKgen[index_tau] = sqrt_absK/k/sin(ptw->chi[index_tau]);
      ptw->cotKgen[index_tau] = ptw->cscKgen[index_tau]*cos(ptw->chi[index_tau]);
    }
    break;
  case 0:
    for (index_tau=0; index_tau < tau_size; index_tau++) {
      ptw->chi[index_tau] = k * tau0_minus_tau[index_tau];
      ptw->cscKgen[index_tau] = 1.0/ptw->chi[index_tau];
      ptw->cotKgen[index_tau] = 1.0/ptw->chi[index_tau];
    }
    break;
  case -1:
    sqrt_absK = sqrt(-ptw->K);
    for (index_tau=0; index_tau < tau_size; index_tau++) {
      ptw->chi[index_tau] = sqrt_absK*tau0_minus_tau[index_tau];
      ptw->cscKgen[index_tau] = sqrt_absK/k/sinh(ptw->chi[index_tau]);
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
 * @param ppt                   Input: pointer to perturbation structure
 * @param ptr                   Input: pointer to transfer structure
 * @param k                     Input: wavenumber at which to interpolate
 * @param index_md              Input: index of mode
 * @param index_ic              Input: index of initial condition
 * @param index_type            Input: index of type of source (in perturbation module)
 * @param pert_source           Input: array of sources
 * @param pert_source_spline    Input: array of second derivative of sources
 * @param interpolated_sources  Output: array of interpolated sources (filled here but allocated in transfer_init() to avoid numerous reallocation)
 * @return the error status
 */

int transfer_interpolate_sources(
                                 struct perturbations * ppt,
                                 struct transfer * ptr,
                                 double k,
                                 int index_md,
                                 int index_ic,
                                 int index_type,
                                 double * pert_source,       /* array with argument pert_source[index_tau*ppt->k_size[index_md]+index_k] (must be allocated) */
                                 double * pert_source_spline, /* array with argument pert_source_spline[index_tau*ppt->k_size[index_md]+index_k] (must be allocated) */
                                 double * interpolated_sources /* array with argument interpolated_sources[index_tau] (must be allocated) */
                                 ) {

  /** Summary: */

  /** - define local variables */

  /* index running on k values in the original source array */
  int index_k;

  /* index running on time */
  int index_tau;

  /* variables used for spline interpolation algorithm */
  double h, a, b;

  /** - interpolate at each k value using the usual
      spline interpolation algorithm. */

  index_k = 0;
  h = ppt->k[index_md][index_k+1] - ppt->k[index_md][index_k];

  while (((index_k+1) < ppt->k_size[index_md]) &&
         (ppt->k[index_md][index_k+1] < k)) {
    index_k++;
    h = ppt->k[index_md][index_k+1] - ppt->k[index_md][index_k];
  }

  class_test(h==0.,
             ptr->error_message,
             "stop to avoid division by zero");

  b = (k - ppt->k[index_md][index_k])/h;
  a = 1.-b;

  for (index_tau = 0; index_tau < ppt->tau_size; index_tau++) {

    interpolated_sources[index_tau] =
      a * pert_source[index_tau*ppt->k_size[index_md]+index_k]
      + b * pert_source[index_tau*ppt->k_size[index_md]+index_k+1]
      + ((a*a*a-a) * pert_source_spline[index_tau*ppt->k_size[index_md]+index_k]
         +(b*b*b-b) * pert_source_spline[index_tau*ppt->k_size[index_md]+index_k+1])*h*h/6.0;

  }

  return _SUCCESS_;

}

/**
 * The code makes a distinction between "perturbation sources"
 * (e.g. gravitational potential) and "transfer sources" (e.g. total
 * density fluctuations, obtained through the Poisson equation, and
 * observed with a given selection function).
 *
 * This routine computes the transfer source of a given type, given
 * the perturbation source interpolated at the right value of k:
 * it resamples the perturbation source at the times of the time
 * sampling of this type, and multiplies it by the window function
 * W(tau) of the type and by the factor (k/k_pivot)^k_power. All these
 * quantities have been precomputed by transfer_time_samplings() and
 * transfer_precompute_windows().
 *
 * @param ptr                   Input: pointer to transfer structure
 * @param interpolated_sources  Input: perturbation source interpolated at k, sampled like in the perturbation module
 * @param k                     Input: wavenumber
 * @param index_md              Input: index of mode
 * @param index_tt              Input: index of type of (transfer) source
 * @param sources               Output: transfer source, sampled like the type
 * @return the error status
 */

int transfer_sources(
                     struct transfer * ptr,
                     double * interpolated_sources,
                     double k,
                     int index_md,
                     int index_tt,
                     double * sources
                     )  {

  /* index running on time */
  int index_tau;

  /* time sampling of this type */
  int index_ts = ptr->ts_of_tt[index_md][index_tt];
  int tau_size = ptr->tau_size_ts[index_md][index_ts];
  int * index_tau_pt = ptr->index_tau_pt[index_md][index_ts];
  double * weight_tau_pt = ptr->weight_tau_pt[index_md][index_ts];

  /* window function of this type */
  double * window = ptr->window[index_md][index_tt];

  /* factor depending on the wavenumber */
  double k_factor = 1.;

  if (ptr->k_power[index_md][index_tt] != 0.)
    k_factor = pow(k/ptr->k_pivot[index_md][index_tt],ptr->k_power[index_md][index_tt]);

  /** - resample the perturbation source by linear interpolation, and multiply it by the window function and k factor */

  for (index_tau = 0; index_tau < tau_size; index_tau++) {
    sources[index_tau] =
      (interpolated_sources[index_tau_pt[index_tau]] * (1.-weight_tau_pt[index_tau])
       + weight_tau_pt[index_tau] * interpolated_sources[index_tau_pt[index_tau]+1])
      * window[index_tau]
      * k_factor;
  }

  return _SUCCESS_;

}

/**
 * Arbitrarily normalized selection function dN/dz(z,bin)
 *
 * @param ppr                   Input: pointer to precision structure
 * @param ppt                   Input: pointer to perturbation structure
 * @param ptr                   Input: pointer to transfer structure
 * @param bin                   Input: redshift bin number
 * @param z                     Input: one value of redshift
 * @param selection             Output: pointer to selection function
 * @return the error status
 */

int transfer_selection_function(
                                struct precision * ppr,
                                struct perturbations * ppt,
                                struct transfer * ptr,
                                int bin,
                                double z,
                                double * selection) {

  double x;
  double dNdz;
  double dln_dNdz_dz;
  int last_index;

  /* trivial dirac case */
  if (ppt->selection==dirac) {

    *selection=1.;

    return _SUCCESS_;
  }

  /* tabulated case: spline interpolation of the selection function
     of the bin, set to zero outside of the tabulated range (and
     wherever the spline would become negative) */
  if (ppt->selection==tabulated) {

    if ((z < ptr->selection_file_z[0]) || (z > ptr->selection_file_z[ptr->selection_file_z_size-1])) {
      *selection = 0.;
    }
    else {
      class_call(array_interpolate_spline(ptr->selection_file_z,
                                          ptr->selection_file_z_size,
                                          ptr->selection_file_dNdz+bin*ptr->selection_file_z_size,
                                          ptr->selection_file_dddNdz+bin*ptr->selection_file_z_size,
                                          1,
                                          z,
                                          &last_index,
                                          selection,
                                          1,
                                          ptr->error_message),
                 ptr->error_message,
                 ptr->error_message);

      *selection = MAX(*selection,0.);
    }

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

    if ((ptr->has_nz_file == _TRUE_) || (ptr->has_nz_analytic == _TRUE_)) {

      if (ptr->has_nz_file == _TRUE_) {

        class_test((z<ptr->nz_z[0]) || (z>ptr->nz_z[ptr->nz_size-1]),
                   ptr->error_message,
                   "Your input file for the selection function only covers the redshift range [%f : %f]. However, your input for the selection function requires z=%f",
                   ptr->nz_z[0],
                   ptr->nz_z[ptr->nz_size-1],
                   z);

        class_call(array_interpolate_spline(
                                            ptr->nz_z,
                                            ptr->nz_size,
                                            ptr->nz_nz,
                                            ptr->nz_ddnz,
                                            1,
                                            z,
                                            &last_index,
                                            &dNdz,
                                            1,
                                            ptr->error_message),
                   ptr->error_message,
                   ptr->error_message);
      }
      else {

        class_call(transfer_dNdz_analytic(ptr,
                                          z,
                                          &dNdz,
                                          &dln_dNdz_dz),
                   ptr->error_message,
                   ptr->error_message);
      }

      *selection *= dNdz;
    }

    return _SUCCESS_;
  }

  /* top-hat case, with smoothed edges. The problem with sharp edges
     is that the final result will be affected by random
     noise. Indeed, the values of k at which the transfer functions
     Delta_l(k) are sampled will never coincide with the actual edges
     of the true transfer function (computed with or even without the
     Limber approximation). Hence the integral Cl=\int dk
     Delta_l(k)**2 (...) will be imprecise and will fluctuate randomly
     with the resolution along k. With smooth edges, the problem is
     solved, and the final Cls become mildly dependent on the
     resolution along k. */

  if (ppt->selection==tophat) {

    /* selection function, centered on z=mean (i.e. on x=0), equal to
       one around x=0, with tanh step centered on x=width, of width
       delta x = 0.1*width
    */
    *selection=(1.-tanh((x-ppt->selection_width[bin])/(ppr->selection_tophat_edge*ppt->selection_width[bin])))/2.;

    if ((ptr->has_nz_file == _TRUE_) || (ptr->has_nz_analytic == _TRUE_)) {

      if (ptr->has_nz_file == _TRUE_) {

        class_call(array_interpolate_spline(
                                            ptr->nz_z,
                                            ptr->nz_size,
                                            ptr->nz_nz,
                                            ptr->nz_ddnz,
                                            1,
                                            z,
                                            &last_index,
                                            &dNdz,
                                            1,
                                            ptr->error_message),
                   ptr->error_message,
                   ptr->error_message);
      }
      else {

        class_call(transfer_dNdz_analytic(ptr,
                                          z,
                                          &dNdz,
                                          &dln_dNdz_dz),
                   ptr->error_message,
                   ptr->error_message);
      }

      *selection *= dNdz;
    }

    return _SUCCESS_;
  }

  /* get here only if selection type was recognized */
  class_stop(ptr->error_message,
             "invalid choice of selection function");

  return _SUCCESS_;
}

/**
 * Analytic form for dNdz distribution, from arXiv:1004.4640
 *
 * @param ptr          Input: pointer to transfer structure
 * @param z            Input: redshift
 * @param dNdz         Output: density per redshift, dN/dZ
 * @param dln_dNdz_dz  Output: dln(dN/dz)/dz, used optionally for the source evolution
 * @return the error status
 */

int transfer_dNdz_analytic(
                           struct transfer * ptr,
                           double z,
                           double * dNdz,
                           double * dln_dNdz_dz) {

  /* Implement here your favorite analytic ansatz for the selection
     function. Typical function for photometric sample: dN/dz =
     (z/z0)^alpha exp[-(z/z0)^beta]. Then: dln(dN/dz)/dz = (alpha -
     beta*(z/z0)^beta)/z. In principle, one is free to use different
     ansatz for the selection function and the evolution
     function. Since the selection function uses only dN/dz, while the
     evolution uses only dln(dN/dz)/dz, it is possible to use
     different functions for dN/dz and dln(dN/dz)/dz */

  double z0,alpha,beta;

  //Euclid IST dNdz, do not change this!
  z0 = 0.9/pow(2.,1./2.);
  alpha = 2.0;
  beta = 1.5;

  *dNdz = pow(z/z0,alpha) * exp(-pow(z/z0,beta));

  *dln_dNdz_dz = (alpha - pow(z/z0,beta)*beta)/z;

  return _SUCCESS_;

}

/**
 * For each selection function, compute the min, mean and max values
 * of conformal time (associated to the min, mean and max values of
 * redshift specified by the user), and store them in the transfer
 * structure
 *
 * @param ppr                   Input: pointer to precision structure
 * @param pba                   Input: pointer to background structure
 * @param ppt                   Input: pointer to perturbation structure
 * @param ptr                   Input/Output: pointer to transfer structure
 * @return the error status
 */

int transfer_selection_times(
                             struct precision * ppr,
                             struct background * pba,
                             struct perturbations * ppt,
                             struct transfer * ptr
                             ) {

  /* redshift bin */
  int bin;

  /* a value of redshift */
  double z=0.;

  for (bin = 0; bin < ppt->selection_num; bin++) {

    /* lower edge of time interval for this bin */
    /* the few lines below should be consistent with their counterpart in input.c */
    if (ppt->selection==gaussian) {
      z = ppt->selection_mean[bin]+ppt->selection_width[bin]*ppr->selection_cut_at_sigma;
    }
    if (ppt->selection==tophat) {
      z = ppt->selection_mean[bin]+(1.+ppr->selection_cut_at_sigma*ppr->selection_tophat_edge)*ppt->selection_width[bin];
    }
    if (ppt->selection==dirac) {
      z = ppt->selection_mean[bin];
    }
    if (ppt->selection==tabulated) {
      z = ptr->selection_z_max[bin];
    }

    class_call(background_tau_of_z(pba,
                                   z,
                                   &(ptr->selection_tau_min[bin])),
               pba->error_message,
               ptr->error_message);

    /* higher edge of time interval for this bin */

    if (ppt->selection==gaussian) {
      z = MAX(ppt->selection_mean[bin]-ppt->selection_width[bin]*ppr->selection_cut_at_sigma,0.);
    }
    if (ppt->selection==tophat) {
      z = MAX(ppt->selection_mean[bin]-(1.+ppr->selection_cut_at_sigma*ppr->selection_tophat_edge)*ppt->selection_width[bin],0.);
    }
    if (ppt->selection==dirac) {
      z = ppt->selection_mean[bin];
    }
    if (ppt->selection==tabulated) {
      z = ptr->selection_z_min[bin];
    }

    class_call(background_tau_of_z(pba,
                                   z,
                                   &(ptr->selection_tau_max[bin])),
               pba->error_message,
               ptr->error_message);

    /* central value of time interval for this bin */

    z = MAX(ppt->selection_mean[bin],0.);

    class_call(background_tau_of_z(pba,
                                   z,
                                   &(ptr->selection_tau_mean[bin])),
               pba->error_message,
               ptr->error_message);
  }

  return _SUCCESS_;

}

/**
 * Compute and normalize selection function for a set of time values
 *
 * @param ppr                   Input: pointer to precision structure
 * @param pba                   Input: pointer to background structure
 * @param ppt                   Input: pointer to perturbation structure
 * @param ptr                   Input: pointer to transfer structure
 * @param selection             Output: normalized selection function
 * @param tau0_minus_tau        Input: values of (tau0-tau) at which source are sample
 * @param w_trapz               Input: trapezoidal weights for integration over tau
 * @param tau_size              Input: size of previous two arrays
 * @param pvecback              Input: allocated array of background values
 * @param tau0                  Input: time today
 * @param bin                   Input: redshift bin number
 * @return the error status
 */

int transfer_selection_compute(
                               struct precision * ppr,
                               struct background * pba,
                               struct perturbations * ppt,
                               struct transfer * ptr,
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

  /* running value of redshift */
  double z;

  if (tau_size > 1) {

    /* loop over time */
    for (index_tau = 0; index_tau < tau_size; index_tau++) {

      /* running value of time */
      tau = tau0 - tau0_minus_tau[index_tau];

      /* get background quantities at this time */
      class_call(background_at_tau(pba,
                                   tau,
                                   long_info,
                                   inter_normal,
                                   &last_index,
                                   pvecback),
                 pba->error_message,
                 ptr->error_message);

      /* infer redshift (remember that a in the code is in fact a/a_0) */
      z = 1./pvecback[pba->index_bg_a]-1.;

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

  }

  /* trivial case: dirac distribution */
  else {
    selection[0] = 1.;
  }

  return _SUCCESS_;
}

/**
 * This routine chooses the method for computing the transfer function
 * \f$ \Delta_l^{X} (q) \f$ of a given mode, type and initial condition,
 * at a given multipole l and wavenumber q: the transfer function can
 * be zero (or neglected), computed with the Limber approximation, or
 * computed by convolving the transfer source with a radial function
 * along the line of sight. In the latter case, the routine also finds
 * the last sampled time contributing to the integral.
 *
 * @param ppr                       Input: pointer to precision structure
 * @param ptr                       Input: pointer to transfer structure
 * @param ptw                       Input: pointer to transfer_workspace structure
 * @param index_q                   Input: index of wavenumber
 * @param index_md                  Input: index of mode
 * @param index_tt                  Input: index of type of transfer
 * @param index_l                   Input: index of multipole
 * @param l                         Input: multipole
 * @param q                         Input: wavenumber
 * @param ra_rec                    Input: angular diameter distance to recombination
 * @param q_max_bessel              Input: maximum value of argument q at which Bessel functions are computed
 * @param tau0_minus_tau_min_bessel Input: minimum value of (tau0-tau) at which the Bessel function is non-negligible
 * @param tau0_minus_tau            Input: time sampling of the transfer source
 * @param tau_size                  Input: number of sampled times
 * @param sources                   Input: transfer source
 * @param use_full_limber           Input: whether we are in the full Limber q-sampling scheme
 * @param method                    Output: method for computing the transfer function
 * @param index_tau_max             Output: for transfer_line_of_sight, index of the last sampled time contributing to the integral
 * @param bessel_truncation         Output: for transfer_line_of_sight, whether the integral is truncated by the cut in the Bessel function
 * @return the error status
 */

int transfer_select_method(
                           struct precision * ppr,
                           struct transfer * ptr,
                           struct transfer_workspace * ptw,
                           int index_q,
                           int index_md,
                           int index_tt,
                           int index_l,
                           double l,
                           double q,
                           double ra_rec,
                           double q_max_bessel,
                           double tau0_minus_tau_min_bessel,
                           double * tau0_minus_tau,
                           int tau_size,
                           double * sources,
                           short use_full_limber,
                           enum transfer_integration_method * method,
                           int * index_tau_max,
                           short * bessel_truncation
                           ){

  /** Summary: */

  /** - define local variables */

  short neglect;
  short use_limber;
  short neglect_late_source;

  /* index of the last time before the cut in the Bessel function */
  int index_tau_max_Bessel;

  *method = transfer_set_zero;
  *index_tau_max = -1;
  *bessel_truncation = _FALSE_;

  /** - return zero transfer function if l is above l_max */
  if (index_l >= ptr->l_size_tt[index_md][index_tt])
    return _SUCCESS_;

  /** - in the full Limber q-sampling scheme, only the types with has_full_limber (the CMB lensing potential) are computed */
  if ((use_full_limber == _TRUE_) && (ptr->has_full_limber[index_md][index_tt] == _FALSE_))
    return _SUCCESS_;

  /** - neglect transfer function when l is much smaller than k*tau0 */
  class_call(transfer_can_be_neglected(ptr,
                                       index_md,
                                       index_tt,
                                       ra_rec,
                                       q,
                                       l,
                                       &neglect),
             ptr->error_message,
             ptr->error_message);
  if (neglect == _TRUE_)
    return _SUCCESS_;

  /** - for K>0 (closed), transfer functions only defined for l<nu */
  if ((ptw->sgnK == 1) && (ptr->l[index_l] >= (int)(q/sqrt(ptw->K)+0.2)))
    return _SUCCESS_;

  /** - in the full Limber q-sampling scheme, always use the Limber approximation */
  if (use_full_limber == _TRUE_) {
    *method = transfer_limber_approximation;
    return _SUCCESS_;
  }

  /** - in curved space, hyperspherical Bessel functions are not computed above some l */
  if ((ptw->sgnK != 0) && (index_q < ptr->index_q_flat_approximation) && (index_l>=ptw->HIS.l_size))
    return _SUCCESS_;

  if (ptr->transfer_verbose > 3)
    printf("Compute transfer for l=%d type=%d\n",(int)l,index_tt);

  /** - decide whether to use the Limber approximation */
  class_call(transfer_use_limber(ptr,
                                 q_max_bessel,
                                 index_md,
                                 index_tt,
                                 q,
                                 l,
                                 &use_limber),
             ptr->error_message,
             ptr->error_message);
  if (use_limber == _TRUE_) {
    *method = transfer_limber_approximation;
    return _SUCCESS_;
  }

  /** - otherwise, line-of-sight integral. If there is no overlap
      between the region in which bessels and sources are non-zero,
      return zero */
  if (tau0_minus_tau_min_bessel >= tau0_minus_tau[0])
    return _SUCCESS_;

  /** - find index in the source's tau list corresponding to the last
      point in the overlapping region. After this step,
      index_tau_max can be as small as zero, but not negative. */
  *index_tau_max = tau_size-1;
  while (tau0_minus_tau[*index_tau_max] < tau0_minus_tau_min_bessel)
    (*index_tau_max)--;
  /* Set index so we know if the truncation of the convolution integral is due to Bessel and not
     due to the source. */
  index_tau_max_Bessel = *index_tau_max;

  /** - the source function can vanish at large \f$ \tau \f$. Check
      if further points can be eliminated. */
  while (sources[*index_tau_max] == 0.) {
    (*index_tau_max)--;
    if (*index_tau_max < 0)
      return _SUCCESS_;
  }

  /** - neglect late time CMB sources when l is above threshold */
  class_call(transfer_late_source_can_be_neglected(ppr,
                                                   ptr,
                                                   index_md,
                                                   index_tt,
                                                   l,
                                                   &neglect_late_source),
             ptr->error_message,
             ptr->error_message);

  if (neglect_late_source == _TRUE_) {
    while (tau0_minus_tau[*index_tau_max] < ptw->tau0_minus_tau_cut) {
      (*index_tau_max)--;
      if (*index_tau_max < 0)
        return _SUCCESS_;
    }
  }

  *method = transfer_line_of_sight;
  *bessel_truncation = ((*index_tau_max != tau_size-1) && (*index_tau_max == index_tau_max_Bessel)) ? _TRUE_ : _FALSE_;

  return _SUCCESS_;

}

/**
 * This routine decides whether the transfer function of a given type
 * and multipole l is computed with the Limber approximation.
 *
 * @param ptr            Input: pointer to transfer structure
 * @param q_max_bessel   Input: maximum value of argument q at which Bessel functions are computed
 * @param index_md       Input: index of mode
 * @param index_tt       Input: index of type
 * @param q              Input: wavenumber
 * @param l              Input: multipole
 * @param use_limber     Output: whether the Limber approximation should be used
 * @return the error status
 */

int transfer_use_limber(
                        struct transfer * ptr,
                        double q_max_bessel,
                        int index_md,
                        int index_tt,
                        double q,
                        double l,
                        short * use_limber) {

  /* criteria for choosing between integration and Limber
     must be implemented here */

  //TBC: in principle the Limber condition should be adapted to account for curvature effects

  /* use the Limber approximation: above the largest q at which
     Bessel functions are tabulated, above the largest q at which
     the time sampling of the type resolves the oscillations of the
     Bessel functions (to avoid aliasing), and above the l
     threshold of the type */

  if ((q > q_max_bessel) ||
      (q > ptr->q_max_sampled[index_md][ptr->ts_of_tt[index_md][index_tt]]) ||
      (l >= ptr->l_limber[index_md][index_tt])) {
    *use_limber = _TRUE_;
  }
  else {
    *use_limber = _FALSE_;
  }

  return _SUCCESS_;
}

/**
 * This routine finds the minimum value of \f$ (\tau_0-\tau) \f$ at
 * which \f$ j_l(k[\tau_0-\tau]) \f$ is known, given that \f$ j_l(x)
 * \f$ is sampled above some finite value \f$ x_{\min} \f$ (below which
 * it can be approximated by zero)
 *
 * @param ptr                       Input: pointer to transfer structure
 * @param ptw                       Input: pointer to transfer_workspace structure
 * @param index_q                   Input: index of wavenumber
 * @param index_l                   Input: index of multipole
 * @param l                         Input: multipole
 * @param k                         Input: wavenumber
 * @param tau0_minus_tau_min_bessel Output: minimum value of (tau0-tau)
 * @return the error status
 */

int transfer_bessel_cut(
                        struct transfer * ptr,
                        struct transfer_workspace * ptw,
                        int index_q,
                        int index_l,
                        double l,
                        double k,
                        double * tau0_minus_tau_min_bessel
                        ) {

  double x_turning_point;

  if (ptw->sgnK==0){
    *tau0_minus_tau_min_bessel = ptw->pBIS->chi_at_phimin[index_l]/k; /* segmentation fault impossible, checked before that k != 0 */
  }
  else{

    if (index_q < ptr->index_q_flat_approximation) {

      /* the hyperspherical Bessel functions are only computed for
         index_l < ptw->HIS.l_size; above, the transfer functions are
         set to zero by transfer_select_method() */
      if (index_l < ptw->HIS.l_size)
        *tau0_minus_tau_min_bessel = ptw->HIS.chi_at_phimin[index_l]/sqrt(ptw->sgnK*ptw->K);
      else
        *tau0_minus_tau_min_bessel = _HUGE_;

    }
    else {

      *tau0_minus_tau_min_bessel = ptw->pBIS->chi_at_phimin[index_l]/sqrt(ptw->sgnK*ptw->K);

      if (ptw->sgnK == 1) {
        x_turning_point = asin(sqrt(l*(l+1.))/ptr->q[index_q]*sqrt(ptw->sgnK*ptw->K));
        *tau0_minus_tau_min_bessel *= x_turning_point/sqrt(l*(l+1.));
      }
      else {
        x_turning_point = asinh(sqrt(l*(l+1.))/ptr->q[index_q]*sqrt(ptw->sgnK*ptw->K));
        *tau0_minus_tau_min_bessel *= x_turning_point/sqrt(l*(l+1.));
      }
    }
  }

  return _SUCCESS_;
}

/**
 * This routine computes the transfer functions \f$ \Delta_l^{X} (k) \f$)
 * for a given mode, type, multipole l and wavenumber k,
 * by convolving the transfer source (passed in input in the array
 * sources) with the radial function of the type, built from the
 * Bessel functions interpolated before in the workspace by
 * transfer_bessel_interpolate().
 *
 * @param ptr                       Input: pointer to transfer structure
 * @param ptw                       Input: pointer to transfer_workspace structure
 * @param index_md                  Input: index of mode
 * @param index_tt                  Input: index of type
 * @param index_l                   Input: index of multipole
 * @param k                         Input: wavenumber
 * @param tau0_minus_tau            Input: time sampling of the transfer source
 * @param w_trapz                   Input: corresponding trapezoidal weights
 * @param sources                   Input: transfer source
 * @param index_tau_max             Input: index of the last sampled time contributing to the integral
 * @param bessel_truncation         Input: whether the integral is truncated by the cut in the Bessel function
 * @param tau0_minus_tau_min_bessel Input: minimum value of (tau0-tau) at which the Bessel function is non-negligible
 * @param trsf                      Output: transfer function \f$ \Delta_l(k) \f$
 * @return the error status
 */

int transfer_integrate(
                       struct transfer * ptr,
                       struct transfer_workspace *ptw,
                       int index_md,
                       int index_tt,
                       int index_l,
                       double k,
                       double * tau0_minus_tau,
                       double * w_trapz,
                       double * sources,
                       int index_tau_max,
                       short bessel_truncation,
                       double tau0_minus_tau_min_bessel,
                       double * trsf
                       ) {

  double * radial_function = ptw->radial_function;

  /** - Compute the radial function: */
  class_call(transfer_radial_function(ptw,
                                      ptr,
                                      k,
                                      index_l,
                                      index_tau_max+1,
                                      radial_function,
                                      ptr->radial_type[index_md][index_tt]),
             ptr->error_message,
             ptr->error_message);

  /** - Now we do most of the convolution integral: */
  class_call(array_trapezoidal_convolution(sources,
                                           radial_function,
                                           index_tau_max+1,
                                           w_trapz,
                                           trsf,
                                           ptr->error_message),
             ptr->error_message,
             ptr->error_message);

  /** - This integral is correct for the case where no truncation has
      occurred. If it has been truncated at some index_tau_max because
      f[index_tau_max+1]==0, it is still correct. The 'mistake' in using
      the wrong weight w_trapz[index_tau_max] is exactly compensated by the
      triangle we miss. However, for the Bessel cut off, we must subtract the
      wrong triangle and add the correct triangle. */
  if (bessel_truncation == _TRUE_) {
    //Bessel truncation
    *trsf -= 0.5*(tau0_minus_tau[index_tau_max+1]-tau0_minus_tau_min_bessel)*
      radial_function[index_tau_max]*sources[index_tau_max];
  }

  return _SUCCESS_;
}

/**
 * This routine computes the transfer functions \f$ \Delta_l^{X} (k) \f$)
 * for a given mode, type, multipole l and wavenumber k,
 * by using the Limber approximation, i.e by evaluating the transfer source
 * (passed in input in the array sources) at a single value of
 * tau (the Bessel function being approximated as a Dirac distribution).
 *
 *
 * @param ptr            Input: pointer to transfer structure
 * @param ptw            Input: pointer to transfer workspace structure
 * @param l              Input: multipole
 * @param q              Input: wavenumber
 * @param tau0_minus_tau Input: time sampling of the transfer source
 * @param tau_size       Input: number of sampled times
 * @param sources        Input: transfer source
 * @param radial_type    Input: type of radial (Bessel) functions to convolve with
 * @param trsf           Output: transfer function \f$ \Delta_l(k) \f$
 * @return the error status
 */

int transfer_limber(
                    struct transfer * ptr,
                    struct transfer_workspace * ptw,
                    double l,
                    double q,
                    double * tau0_minus_tau,
                    int tau_size,
                    double * sources,
                    radial_function_type radial_type,
                    double * trsf
                    ){

  /** Summary: */

  /** - define local variables */

  /* interpolated source and its derivatives at this value */
  double S, Sp, Sm;

  double x_limber=0.;
  double tau0_minus_tau_limber=0.;
  double IPhiFlat = 0.;

  if (radial_type == SCALAR_TEMPERATURE_0) {

    /** - get k, l and infer tau such that k(tau0-tau)=l+1/2;
        check that tau is in appropriate range */

    if (ptw->sgnK == 0) {
      tau0_minus_tau_limber = (l+0.5)/q;
    }
    else if (ptw->sgnK == 1) {
      x_limber = asin(sqrt(l*(l+1.))/q*sqrt(ptw->K));
      tau0_minus_tau_limber = x_limber/sqrt(ptw->K);
    }
    else if (ptw->sgnK == -1) {
      x_limber = asinh((l+0.5)/q*sqrt(-ptw->K));
      tau0_minus_tau_limber = x_limber/sqrt(-ptw->K);
    }

    if ((tau0_minus_tau_limber > tau0_minus_tau[0]) ||
        (tau0_minus_tau_limber < tau0_minus_tau[tau_size-1])) {
      *trsf = 0.;
      return _SUCCESS_;
    }

    class_call(transfer_limber_interpolate(ptr,
                                           tau0_minus_tau,
                                           sources,
                                           tau_size,
                                           tau0_minus_tau_limber,
                                           &S),
               ptr->error_message,
               ptr->error_message);

    /** - get transfer = source * \f$ \sqrt{\pi/(2l+1)}/q \f$
        = source*[tau0-tau] * \f$ \sqrt{\pi/(2l+1)}/(l+1/2)\f$
    */

    IPhiFlat = sqrt(_PI_/(2.*l))*(1.-0.25/l+1./32./(l*l));

    *trsf = IPhiFlat*S;

    if (ptw->sgnK == 0) {
      *trsf /= (l+0.5);
    }
    else {
      *trsf *= pow(1.-ptw->K*l*l/q/q,-1./4.)/(tau0_minus_tau_limber*q);
    }

  }

  else if (radial_type == SCALAR_TEMPERATURE_1) {

    if (((l+1.5)/q > tau0_minus_tau[0]) ||
        ((l-0.5)/q < tau0_minus_tau[tau_size-1])) {
      *trsf = 0.;
      return _SUCCESS_;
    }

    class_call(transfer_limber_interpolate(ptr,
                                           tau0_minus_tau,
                                           sources,
                                           tau_size,
                                           (l+1.5)/q,
                                           &Sp),
               ptr->error_message,
               ptr->error_message);

    class_call(transfer_limber_interpolate(ptr,
                                           tau0_minus_tau,
                                           sources,
                                           tau_size,
                                           (l-0.5)/q,
                                           &Sm),
               ptr->error_message,
               ptr->error_message);

    *trsf =
      -sqrt(_PI_/(2.*l+3.))*Sp/(l+1.5) * (l+1.)/(2.*l+1)
      +sqrt(_PI_/(2.*l-1.))*Sm/(l-0.5) * l/(2.*l+1.);

  }

  else if (radial_type == NC_RSD) {

    if (((l+2.5)/q > tau0_minus_tau[0]) ||
        ((l-1.5)/q < tau0_minus_tau[tau_size-1])) {
      *trsf = 0.;
      return _SUCCESS_;
    }

    class_call(transfer_limber_interpolate(ptr,
                                           tau0_minus_tau,
                                           sources,
                                           tau_size,
                                           (l+2.5)/q,
                                           &Sp),
               ptr->error_message,
               ptr->error_message);

    class_call(transfer_limber_interpolate(ptr,
                                           tau0_minus_tau,
                                           sources,
                                           tau_size,
                                           (l-1.5)/q,
                                           &Sm),
               ptr->error_message,
               ptr->error_message);

    class_call(transfer_limber_interpolate(ptr,
                                           tau0_minus_tau,
                                           sources,
                                           tau_size,
                                           (l+0.5)/q,
                                           &S),
               ptr->error_message,
               ptr->error_message);

    *trsf =
      sqrt(_PI_/(2.*l+5.))*Sp/(l+2.5) * l*(l+2.)/(2.*l+1.)/(2.*l+3.)
      -sqrt(_PI_/(2.*l+1.))*S/(l+0.5) * l/(2.*l+1.)*(l/(2.*l-1.)+(l+1.)/(2.*l+3.))
      +sqrt(_PI_/(2.*l-3.))*Sm/(l-1.5) * l*(l-1.)/(2.*l+1.)/(2.*l-1.);

  }

  else {

    class_stop(ptr->error_message,
               "Limber approximation has not been coded for the radial_type of index %d\n",
               radial_type);

  }

  return _SUCCESS_;

}

/**
 * This routine interpolates the transfer source at the time given by
 * the Limber approximation (tau0-tau = (l+1/2)/q in the flat case),
 * for transfer_limber(). It uses a polynomial of order two through three
 * neighboring sampled times, fitted to the product of the source by
 * (tau0-tau) rather than to the source itself, because this product
 * remains regular at tau=tau0 for the lensing sources.
 *
 * @param ptr                   Input: pointer to transfer structure
 * @param tau0_minus_tau        Input: time sampling of the transfer source, in decreasing order
 * @param sources               Input: transfer source
 * @param tau_size              Input: number of sampled times
 * @param tau0_minus_tau_limber Input: value of (tau0-tau) at which the source is interpolated
 * @param S                     Output: interpolated source
 * @return the error status
 */

int transfer_limber_interpolate(
                                struct transfer * ptr,
                                double * tau0_minus_tau,
                                double * sources,
                                int tau_size,
                                double tau0_minus_tau_limber,
                                double * S
                                ){

  int index_tau;
  double dS,ddS;

  /** - find the central point index_tau of the three points used for
      the interpolation below. index_tau must be at least 1 (so that
      index_tau-1 is at least 0) and at most tau_size-2 (so that
      index_tau+1 is at most tau_size-1): it is the smallest index in
      this range such that tau0_minus_tau[index_tau] <
      tau0_minus_tau_limber, or tau_size-2 if there is none. */

  if ((tau_size <= 3) || (tau0_minus_tau[1] < tau0_minus_tau_limber)) {
    index_tau = 1;
  }
  else if (tau0_minus_tau[tau_size-2] >= tau0_minus_tau_limber) {
    index_tau = tau_size-2;
  }
  else {
    /* here tau0_minus_tau[1] >= tau0_minus_tau_limber > tau0_minus_tau[tau_size-2]: since tau0_minus_tau is decreasing,
       array_search_bisect() returns index_tau-1, such that tau0_minus_tau[index_tau-1] >= tau0_minus_tau_limber > tau0_minus_tau[index_tau] */
    class_call(array_search_bisect(tau_size,
                                   tau0_minus_tau,
                                   tau0_minus_tau_limber,
                                   &index_tau,
                                   ptr->error_message),
               ptr->error_message,
               ptr->error_message);
    index_tau++;
  }

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
                                          S,
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
                                          S,
                                          &dS,
                                          &ddS,
                                          ptr->error_message),
               ptr->error_message,
               ptr->error_message);
  }

  return _SUCCESS_;

}

/**
 * This routine computes the transfer functions \f$ \Delta_l^{X} (k)
 * \f$) for each mode, initial condition, type, multipole l and
 * wavenumber k, by using the Limber approximation at order two, i.e
 * as a function of the source function and its first two derivatives
 * at a single value of tau
 *
 * @param tau_size        Input: size of conformal time array
 * @param ptr             Input: pointer to transfer structure
 * @param index_md        Input: index of mode
 * @param index_k         Input: index of wavenumber
 * @param l               Input: multipole
 * @param k               Input: wavenumber
 * @param tau0_minus_tau  Input: array of values of (tau_today - tau)
 * @param sources         Input: source functions
 * @param radial_type     Input: type of radial (Bessel) functions to convolve with
 * @param trsf            Output: transfer function \f$ \Delta_l(k) \f$
 * @return the error status
 */

int transfer_limber2(
                     int tau_size,
                     struct transfer * ptr,
                     int index_md,
                     int index_k,
                     double l,
                     double k,
                     double * tau0_minus_tau,
                     double * sources,
                     radial_function_type radial_type,
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


  /** - get transfer from 2nd order Limber approx (inferred from 0809.5112 [astro-ph]) */

  *trsf = sqrt(_PI_/(2.*l+1.))/k*((1.-3./2./(2.*l+1.)/(2.*l+1.))*S+dS/k/(2.*l+1.)-0.5*ddS/k/k);

  return _SUCCESS_;

}

/**
 * This routine decides whether the transfer function of a given type
 * can be neglected (set to zero) at a given wavenumber and multipole,
 * because l is much smaller than q times the angular diameter
 * distance to recombination (for the types where this is relevant,
 * see transfer_type_properties())
 *
 * @param ptr      Input: pointer to transfer structure
 * @param index_md Input: index of mode
 * @param index_tt Input: index of type
 * @param ra_rec   Input: angular diameter distance to recombination
 * @param q        Input: wavenumber
 * @param l        Input: multipole
 * @param neglect  Output: whether the transfer function can be neglected
 * @return the error status
 */

int transfer_can_be_neglected(
                              struct transfer * ptr,
                              int index_md,
                              int index_tt,
                              double ra_rec,
                              double q,
                              double l,
                              short * neglect) {

  if (l < (q-ptr->delta_k_neglect[index_md][index_tt])*ra_rec)
    *neglect = _TRUE_;
  else
    *neglect = _FALSE_;

  return _SUCCESS_;

}

/**
 * This routine decides whether the source of a given type can be
 * neglected at late times (tau > tau_cut) for a given multipole. This
 * is the case above l = ppr->transfer_neglect_late_source for CMB
 * types without late ISW contribution (see
 * transfer_type_properties())
 *
 * @param ppr      Input: pointer to precision structure
 * @param ptr      Input: pointer to transfer structure
 * @param index_md Input: index of mode
 * @param index_tt Input: index of type
 * @param l        Input: multipole
 * @param neglect  Output: whether the late source can be neglected
 * @return the error status
 */

int transfer_late_source_can_be_neglected(
                                          struct precision * ppr,
                                          struct transfer * ptr,
                                          int index_md,
                                          int index_tt,
                                          double l,
                                          short * neglect) {

  if ((ptr->has_late_source_cut[index_md][index_tt] == _TRUE_) &&
      (l > ppr->transfer_neglect_late_source*ptr->angular_rescaling))
    *neglect = _TRUE_;
  else
    *neglect = _FALSE_;

  return _SUCCESS_;

}

/**
 * This routine tells which Bessel functions and derivatives are
 * needed for building the radial function of a given type in
 * transfer_radial_function().
 *
 * @param radial_type Input: type of radial function
 * @param need_Phi    Output: is the Bessel function needed?
 * @param need_dPhi   Output: is its first derivative needed?
 * @param need_d2Phi  Output: is its second derivative needed?
 * @return the error status
 */

int transfer_radial_function_needs(
                                   radial_function_type radial_type,
                                   short * need_Phi,
                                   short * need_dPhi,
                                   short * need_d2Phi
                                   ) {

  *need_Phi = _FALSE_;
  *need_dPhi = _FALSE_;
  *need_d2Phi = _FALSE_;

  switch (radial_type){
  case SCALAR_TEMPERATURE_0:
  case SCALAR_POLARISATION_E:
  case VECTOR_TEMPERATURE_1:
  case VECTOR_POLARISATION_B:
  case TENSOR_TEMPERATURE_2:
    *need_Phi = _TRUE_;
    break;
  case SCALAR_TEMPERATURE_1:
    *need_dPhi = _TRUE_;
    break;
  case SCALAR_TEMPERATURE_2:
    *need_Phi = _TRUE_;
    *need_d2Phi = _TRUE_;
    break;
  case VECTOR_TEMPERATURE_2:
  case VECTOR_POLARISATION_E:
  case TENSOR_POLARISATION_B:
    *need_Phi = _TRUE_;
    *need_dPhi = _TRUE_;
    break;
  case TENSOR_POLARISATION_E:
    *need_Phi = _TRUE_;
    *need_dPhi = _TRUE_;
    *need_d2Phi = _TRUE_;
    break;
  case NC_RSD:
    *need_d2Phi = _TRUE_;
    break;
  }

  return _SUCCESS_;
}

/**
 * This routine interpolates the Bessel functions (hyperspherical
 * Bessel functions in non-flat space) and the requested derivatives
 * at the first x_size radial coordinates of the current time sampling
 * (stored in the workspace by transfer_radial_coordinates()). The
 * results are stored in the workspace, and can then be used for
 * building the radial functions of all the types sharing this time
 * sampling with transfer_radial_function().
 *
 * @param ptw         Input/Output: pointer to workspace
 * @param ptr         Input: pointer to transfer structure
 * @param index_q     Input: index of wavenumber
 * @param index_l     Input: index of multipole
 * @param x_size      Input: number of values
 * @param need_Phi    Input: is the Bessel function needed?
 * @param need_dPhi   Input: is its first derivative needed?
 * @param need_d2Phi  Input: is its second derivative needed?
 * @return the error status
 */

int transfer_bessel_interpolate(
                                struct transfer_workspace * ptw,
                                struct transfer * ptr,
                                int index_q,
                                int index_l,
                                int x_size,
                                short need_Phi,
                                short need_dPhi,
                                short need_d2Phi
                                ){

  HyperInterpStruct * pHIS;
  double *chi = ptw->chi;
  double *Phi = ptw->Phi;
  double *dPhi = ptw->dPhi;
  double *d2Phi = ptw->d2Phi;
  double *chireverse = ptw->chireverse;
  double *rescale_function = ptw->rescale_function;
  int j;
  double K=0.;
  double nu=0., chi_tp=0.;
  double rescale_argument;
  double rescale_amplitude;
  int (*interpolate_Phi)(HyperInterpStruct*, int, int, double*, double*, char*);
  int (*interpolate_dPhi)(HyperInterpStruct*, int, int, double*, double*, char*);
  int (*interpolate_d2Phi)(HyperInterpStruct*, int, int, double*, double*, char*);
  int (*interpolate_PhidPhi)(HyperInterpStruct*, int, int, double*, double*, double*, char*);
  int (*interpolate_Phid2Phi)(HyperInterpStruct*, int, int, double*, double*, double*, char*);
  int (*interpolate_dPhid2Phi)(HyperInterpStruct*, int, int, double*, double*, double*, char*);
  int (*interpolate_PhidPhid2Phi)(HyperInterpStruct*, int, int, double*, double*, double*, double*, char*);
  enum Hermite_Interpolation_Order HIorder;

  K = ptw->K;

  if (ptw->sgnK == 0) {
    pHIS = ptw->pBIS;
    rescale_argument = 1.;
    rescale_amplitude = 1.;
    HIorder = HERMITE4;
  }
  else if (index_q < ptr->index_q_flat_approximation) {
    pHIS = &(ptw->HIS);
    rescale_argument = 1.;
    rescale_amplitude = 1.;
    HIorder = HERMITE6;
  }
  else {
    pHIS = ptw->pBIS;
    if (ptw->sgnK == 1){
      nu = ptr->q[index_q]/sqrt(K);
      chi_tp = asin(sqrt(ptr->l[index_l]*(ptr->l[index_l]+1.))/nu);
    }
    else{
      nu = ptr->q[index_q]/sqrt(-K);
      chi_tp = asinh(sqrt(ptr->l[index_l]*(ptr->l[index_l]+1.))/nu);
    }
    rescale_argument = sqrt(ptr->l[index_l]*(ptr->l[index_l]+1.))/chi_tp;
    rescale_amplitude = pow(1.-K*ptr->l[index_l]*(ptr->l[index_l]+1.)/ptr->q[index_q]/ptr->q[index_q],-1./12.);
    HIorder = HERMITE4;
  }

  ptw->bessel_size = x_size;
  ptw->rescale_argument = rescale_argument;

  switch (HIorder){
  case HERMITE3:
    interpolate_Phi = hyperspherical_Hermite3_interpolation_vector_Phi;
    interpolate_dPhi = hyperspherical_Hermite3_interpolation_vector_dPhi;
    interpolate_d2Phi = hyperspherical_Hermite3_interpolation_vector_d2Phi;
    interpolate_PhidPhi = hyperspherical_Hermite3_interpolation_vector_PhidPhi;
    interpolate_Phid2Phi = hyperspherical_Hermite3_interpolation_vector_Phid2Phi;
    interpolate_dPhid2Phi = hyperspherical_Hermite3_interpolation_vector_dPhid2Phi;
    interpolate_PhidPhid2Phi = hyperspherical_Hermite3_interpolation_vector_PhidPhid2Phi;
    break;
  case HERMITE4:
    interpolate_Phi = hyperspherical_Hermite4_interpolation_vector_Phi;
    interpolate_dPhi = hyperspherical_Hermite4_interpolation_vector_dPhi;
    interpolate_d2Phi = hyperspherical_Hermite4_interpolation_vector_d2Phi;
    interpolate_PhidPhi = hyperspherical_Hermite4_interpolation_vector_PhidPhi;
    interpolate_Phid2Phi = hyperspherical_Hermite4_interpolation_vector_Phid2Phi;
    interpolate_dPhid2Phi = hyperspherical_Hermite4_interpolation_vector_dPhid2Phi;
    interpolate_PhidPhid2Phi = hyperspherical_Hermite4_interpolation_vector_PhidPhid2Phi;
    break;
  case HERMITE6:
    interpolate_Phi = hyperspherical_Hermite6_interpolation_vector_Phi;
    interpolate_dPhi = hyperspherical_Hermite6_interpolation_vector_dPhi;
    interpolate_d2Phi = hyperspherical_Hermite6_interpolation_vector_d2Phi;
    interpolate_PhidPhi = hyperspherical_Hermite6_interpolation_vector_PhidPhi;
    interpolate_Phid2Phi = hyperspherical_Hermite6_interpolation_vector_Phid2Phi;
    interpolate_dPhid2Phi = hyperspherical_Hermite6_interpolation_vector_dPhid2Phi;
    interpolate_PhidPhid2Phi = hyperspherical_Hermite6_interpolation_vector_PhidPhid2Phi;
    break;
  }

  //Reverse chi
  for (j=0; j<x_size; j++) {
    chireverse[j] = chi[x_size-1-j]*rescale_argument;
    if (rescale_amplitude == 1.) {
      rescale_function[j] = 1.;
    }
    else {
      if (ptw->sgnK == 1) {
        rescale_function[j] =
          MIN(
              rescale_amplitude
              * (1
                 + 0.34 * atan(ptr->l[index_l]/nu) * (chireverse[j]/rescale_argument-chi_tp)
                 + 2.00 * pow(atan(ptr->l[index_l]/nu) * (chireverse[j]/rescale_argument-chi_tp),2)),
              chireverse[j]/rescale_argument/sin(chireverse[j]/rescale_argument)
              );
      }
      else {
        rescale_function[j] =
          MAX(
              rescale_amplitude
              * (1
                 - 0.38 * atan(ptr->l[index_l]/nu) * (chireverse[j]/rescale_argument-chi_tp)
                 + 0.40 * pow(atan(ptr->l[index_l]/nu) * (chireverse[j]/rescale_argument-chi_tp),2)),
              chireverse[j]/rescale_argument/sinh(chireverse[j]/rescale_argument)
              );
      }
    }
  }

  /*
    class_test(pHIS->x[0] > chireverse[0],
    ptr->error_message,
    "Bessels need to be interpolated at %e, outside the range in which they have been computed (>%e). Decrease their x_min.",
    chireverse[0],
    pHIS->x[0]);
  */

  class_test((pHIS->x[pHIS->x_size-1] < chireverse[x_size-1]) && (ptw->sgnK != 1),
             ptr->error_message,
             "Bessels need to be interpolated at %e, outside the range in which they have been computed (<%e). Increase their x_max.",
             chireverse[x_size-1],
             pHIS->x[pHIS->x_size-1]
             );

  /* interpolate only the needed functions */

  if (need_Phi == _TRUE_) {
    if (need_dPhi == _TRUE_) {
      if (need_d2Phi == _TRUE_) {
        class_call(interpolate_PhidPhid2Phi(pHIS, x_size, index_l, chireverse, Phi, dPhi, d2Phi, ptr->error_message),
                   ptr->error_message, ptr->error_message);
      }
      else {
        class_call(interpolate_PhidPhi(pHIS, x_size, index_l, chireverse, Phi, dPhi, ptr->error_message),
                   ptr->error_message, ptr->error_message);
      }
    }
    else {
      if (need_d2Phi == _TRUE_) {
        class_call(interpolate_Phid2Phi(pHIS, x_size, index_l, chireverse, Phi, d2Phi, ptr->error_message),
                   ptr->error_message, ptr->error_message);
      }
      else {
        class_call(interpolate_Phi(pHIS, x_size, index_l, chireverse, Phi, ptr->error_message),
                   ptr->error_message, ptr->error_message);
      }
    }
  }
  else {
    if (need_dPhi == _TRUE_) {
      if (need_d2Phi == _TRUE_) {
        class_call(interpolate_dPhid2Phi(pHIS, x_size, index_l, chireverse, dPhi, d2Phi, ptr->error_message),
                   ptr->error_message, ptr->error_message);
      }
      else {
        class_call(interpolate_dPhi(pHIS, x_size, index_l, chireverse, dPhi, ptr->error_message),
                   ptr->error_message, ptr->error_message);
      }
    }
    else {
      if (need_d2Phi == _TRUE_) {
        class_call(interpolate_d2Phi(pHIS, x_size, index_l, chireverse, d2Phi, ptr->error_message),
                   ptr->error_message, ptr->error_message);
      }
    }
  }

  return _SUCCESS_;
}

/**
 * This routine computes the radial function (Bessel function or
 * combination of its derivatives) convolved with the transfer source
 * of a given type, at the first x_size times of the current time
 * sampling, using the Bessel functions interpolated before by
 * transfer_bessel_interpolate() (for at least x_size times).
 *
 * @param ptw             Input: pointer to workspace (radial coordinates and Bessel functions)
 * @param ptr             Input: pointer to transfer structure
 * @param k               Input: wavenumber
 * @param index_l         Input: index of multipole
 * @param x_size          Input: number of values
 * @param radial_function Output: radial function, in the order of the time sampling
 * @param radial_type     Input: type of radial function
 * @return the error status
 */

int transfer_radial_function(
                             struct transfer_workspace * ptw,
                             struct transfer * ptr,
                             double k,
                             int index_l,
                             int x_size,
                             double * radial_function,
                             radial_function_type radial_type
                             ){

  double *cscKgen = ptw->cscKgen;
  double *cotKgen = ptw->cotKgen;
  /* Bessel functions are stored in reverse order: the value at index_tau of the time sampling is at index n-index_tau */
  double *Phi = ptw->Phi;
  double *dPhi = ptw->dPhi;
  double *d2Phi = ptw->d2Phi;
  double *rescale_function = ptw->rescale_function;
  int n = ptw->bessel_size-1;
  int i;
  double K=0.,k2=1.0;
  double sqrt_absK_over_k;
  double absK_over_k2;
  double factor, s0, s2, ssqrt3, si, ssqrt2, ssqrt2i;
  double l = (double)ptr->l[index_l];
  double rescale_argument = ptw->rescale_argument;

  class_test(x_size > ptw->bessel_size,
             ptr->error_message,
             "radial function needed at %d times, but Bessel functions interpolated at only %d times",
             x_size,ptw->bessel_size);

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

  switch (radial_type){
  case SCALAR_TEMPERATURE_0:
    for (i=0; i<x_size; i++)
      radial_function[i] = Phi[n-i]*rescale_function[n-i];
    break;
  case SCALAR_TEMPERATURE_1:
    for (i=0; i<x_size; i++)
      radial_function[i] = sqrt_absK_over_k*dPhi[n-i]*rescale_argument*rescale_function[n-i];
    break;
  case SCALAR_TEMPERATURE_2:
    s2 = sqrt(1.0-3.0*K/k2);
    factor = 1.0/(2.0*s2);
    for (i=0; i<x_size; i++)
      radial_function[i] = factor*(3*absK_over_k2*d2Phi[n-i]*rescale_argument*rescale_argument+Phi[n-i])*rescale_function[n-i];
    break;
  case SCALAR_POLARISATION_E:
    s2 = sqrt(1.0-3.0*K/k2);
    factor = sqrt(3.0/8.0*(l+2.0)*(l+1.0)*l*(l-1.0))/s2;
    for (i=0; i<x_size; i++)
      radial_function[i] = factor*cscKgen[i]*cscKgen[i]*Phi[n-i]*rescale_function[n-i];
    break;
  case VECTOR_TEMPERATURE_1:
    s0 = sqrt(1.0+K/k2);
    factor = sqrt(0.5*l*(l+1))/s0;
    for (i=0; i<x_size; i++)
      radial_function[i] = factor*cscKgen[i]*Phi[n-i]*rescale_function[n-i];
    break;
  case VECTOR_TEMPERATURE_2:
    s0 = sqrt(1.0+K/k2);
    ssqrt3 = sqrt(1.0-2.0*K/k2);
    factor = sqrt(1.5*l*(l+1))/s0/ssqrt3;
    for (i=0; i<x_size; i++)
      radial_function[i] = factor*cscKgen[i]*(sqrt_absK_over_k*dPhi[n-i]*rescale_argument-cotKgen[i]*Phi[n-i])*rescale_function[n-i];
    break;
  case VECTOR_POLARISATION_E:
    s0 = sqrt(1.0+K/k2);
    ssqrt3 = sqrt(1.0-2.0*K/k2);
    factor = 0.5*sqrt((l-1.0)*(l+2.0))/s0/ssqrt3;
    for (i=0; i<x_size; i++)
      radial_function[i] = factor*cscKgen[i]*(cotKgen[i]*Phi[n-i]+sqrt_absK_over_k*dPhi[n-i]*rescale_argument)*rescale_function[n-i];
    break;
  case VECTOR_POLARISATION_B:
    s0 = sqrt(1.0+K/k2);
    ssqrt3 = sqrt(1.0-2.0*K/k2);
    si = sqrt(1.0+2.0*K/k2);
    factor = 0.5*sqrt((l-1.0)*(l+2.0))*si/s0/ssqrt3;
    for (i=0; i<x_size; i++)
      radial_function[i] = factor*cscKgen[i]*Phi[n-i]*rescale_function[n-i];
    break;
  case TENSOR_TEMPERATURE_2:
    ssqrt2 = sqrt(1.0-1.0*K/k2);
    si = sqrt(1.0+2.0*K/k2);
    factor = sqrt(3.0/8.0*(l+2.0)*(l+1.0)*l*(l-1.0))/si/ssqrt2;
    for (i=0; i<x_size; i++)
      radial_function[i] = factor*cscKgen[i]*cscKgen[i]*Phi[n-i]*rescale_function[n-i];
    break;
  case TENSOR_POLARISATION_E:
    ssqrt2 = sqrt(1.0-1.0*K/k2);
    si = sqrt(1.0+2.0*K/k2);
    factor = 0.25/si/ssqrt2;
    for (i=0; i<x_size; i++)
      radial_function[i] = factor*(absK_over_k2*d2Phi[n-i]*rescale_argument*rescale_argument
                                   +4.0*cotKgen[i]*sqrt_absK_over_k*dPhi[n-i]*rescale_argument
                                   -(1.0+4*K/k2-2.0*cotKgen[i]*cotKgen[i])*Phi[n-i])*rescale_function[n-i];
    break;
  case TENSOR_POLARISATION_B:
    ssqrt2i = sqrt(1.0+3.0*K/k2);
    ssqrt2 = sqrt(1.0-1.0*K/k2);
    si = sqrt(1.0+2.0*K/k2);
    factor = 0.5*ssqrt2i/ssqrt2/si;
    for (i=0; i<x_size; i++)
      radial_function[i] = factor*(sqrt_absK_over_k*dPhi[n-i]*rescale_argument+2.0*cotKgen[i]*Phi[n-i])*rescale_function[n-i];
    break;
  case NC_RSD:
    factor = 1.0;
    for (i=0; i<x_size; i++)
      radial_function[i] = factor*absK_over_k2*d2Phi[n-i]*rescale_argument*rescale_argument*rescale_function[n-i];
    // Note: in previous line there was a missing factor absK_over_k2 until version 2.4.3. Credits Francesco Montanari.
    break;
  }

  return _SUCCESS_;
}

/**
 * This routine reads the optional files describing the global
 * selection function dN/dz (multiplying the selection function of
 * each bin, input dNdz_selection) and the evolution of the number of
 * sources (input dNdz_evolution), and splines them for later
 * interpolation. Each file contains the redshift in the first column
 * and dN/dz in the second one. For the evolution, the code splines
 * d ln(dN/dz)/dz, computed by finite differences.
 *
 * @param ptr Input/Output: pointer to transfer structure
 * @return the error status
 */

int transfer_global_selection_read(
                                   struct transfer * ptr
                                   ) {

  /* for reading selection function */
  FILE * input_file;
  int row,status;
  double tmp1,tmp2;

  ptr->nz_size = 0;

  if (ptr->has_nz_file == _TRUE_) {

    input_file = fopen(ptr->nz_file_name,"r");
    class_test(input_file == NULL,
               ptr->error_message,
               "Could not open file %s!",ptr->nz_file_name);

    /* Find size of table */
    for (row=0,status=2; status==2; row++){
      status = fscanf(input_file,"%lf %lf",&tmp1,&tmp2);
    }
    rewind(input_file);
    ptr->nz_size = row-1;

    /* Allocate room for interpolation table */
    class_alloc(ptr->nz_z,sizeof(double)*ptr->nz_size,ptr->error_message);
    class_alloc(ptr->nz_nz,sizeof(double)*ptr->nz_size,ptr->error_message);
    class_alloc(ptr->nz_ddnz,sizeof(double)*ptr->nz_size,ptr->error_message);

    for (row=0; row<ptr->nz_size; row++){
      status = fscanf(input_file,"%lf %lf",
                      &ptr->nz_z[row],&ptr->nz_nz[row]);
      //printf("%d: (z,dNdz) = (%g,%g)\n",row,ptr->nz_z[row],ptr->nz_nz[row]);
    }
    fclose(input_file);

    /* Call spline interpolation: */
    class_call(array_spline_table_lines(ptr->nz_z,
                                        ptr->nz_size,
                                        ptr->nz_nz,
                                        1,
                                        ptr->nz_ddnz,
                                        _SPLINE_EST_DERIV_,
                                        ptr->error_message),
               ptr->error_message,
               ptr->error_message);
  }

  ptr->nz_evo_size = 0;

  if (ptr->has_nz_evo_file == _TRUE_) {

    input_file = fopen(ptr->nz_evo_file_name,"r");
    class_test(input_file == NULL,
               ptr->error_message,
               "Could not open file %s!",ptr->nz_evo_file_name);

    /* Find size of table */
    for (row=0,status=2; status==2; row++){
      status = fscanf(input_file,"%lf %lf",&tmp1,&tmp2);
    }
    rewind(input_file);
    ptr->nz_evo_size = row-1;

    /* Allocate room for interpolation table */
    class_alloc(ptr->nz_evo_z,sizeof(double)*ptr->nz_evo_size,ptr->error_message);
    class_alloc(ptr->nz_evo_nz,sizeof(double)*ptr->nz_evo_size,ptr->error_message);
    class_alloc(ptr->nz_evo_dlog_nz,sizeof(double)*ptr->nz_evo_size,ptr->error_message);
    class_alloc(ptr->nz_evo_dd_dlog_nz,sizeof(double)*ptr->nz_evo_size,ptr->error_message);

    for (row=0; row<ptr->nz_evo_size; row++){
      status = fscanf(input_file,"%lf %lf",
                      &ptr->nz_evo_z[row],&ptr->nz_evo_nz[row]);
    }
    fclose(input_file);

    /* infer dlog(dN/dz)/dz from dN/dz */
    ptr->nz_evo_dlog_nz[0] =
      (log(ptr->nz_evo_nz[1])-log(ptr->nz_evo_nz[0]))
      /(ptr->nz_evo_z[1]-ptr->nz_evo_z[0]);
    for (row=1; row<ptr->nz_evo_size-1; row++){
      ptr->nz_evo_dlog_nz[row] =
        (log(ptr->nz_evo_nz[row+1])-log(ptr->nz_evo_nz[row-1]))
        /(ptr->nz_evo_z[row+1]-ptr->nz_evo_z[row-1]);
    }
    ptr->nz_evo_dlog_nz[ptr->nz_evo_size-1] =
      (log(ptr->nz_evo_nz[ptr->nz_evo_size-1])-log(ptr->nz_evo_nz[ptr->nz_evo_size-2]))
      /(ptr->nz_evo_z[ptr->nz_evo_size-1]-ptr->nz_evo_z[ptr->nz_evo_size-2]);

    /* to test that the file is read:
       for (row=0; row<ptr->nz_evo_size; row++){
       fprintf(stdout,"%d: (z,dNdz,dlndNdzdz) = (%g,%g,%g)\n",row,ptr->nz_evo_z[row],ptr->nz_evo_nz[row],ptr->nz_evo_dlog_nz[row]);
       }
    */

    /* Call spline interpolation: */
    class_call(array_spline_table_lines(ptr->nz_evo_z,
                                        ptr->nz_evo_size,
                                        ptr->nz_evo_dlog_nz,
                                        1,
                                        ptr->nz_evo_dd_dlog_nz,
                                        _SPLINE_EST_DERIV_,
                                        ptr->error_message),
               ptr->error_message,
               ptr->error_message);
  }

  return _SUCCESS_;

};

/**
 * Read a file containing the selection function dN/dz of several
 * redshift bins (for selection=tabulated): the redshift z in the first
 * column, and the selection function of each bin (arbitrarily
 * normalized) in each of the next columns. Lines starting with '#' are
 * ignored. Also computes for each bin the mean redshift and the range
 * in which the selection function is non-zero.
 *
 * This function is called both by the input module (to infer the number
 * of bins, their mean redshift and their redshift range) and by the
 * transfer module (to store the tabulated selection functions).
 *
 * @param file_name     Input: name of the file
 * @param z_size_out    Output: number of redshift values
 * @param bin_num_out   Output: number of bins
 * @param z_out         Output: redshift values (allocated here)
 * @param dNdz_out      Output: selection function of each bin, dNdz[bin*z_size+index_z] (allocated here)
 * @param z_mean        Output: mean redshift of each bin, weighted by dN/dz (array of size at least _SELECTION_NUM_MAX_)
 * @param z_min         Output: smallest redshift at which the selection function of each bin is non-zero (same size)
 * @param z_max         Output: largest redshift at which the selection function of each bin is non-zero (same size)
 * @param error_message Output: error message
 * @return the error status
 */

int transfer_selection_file_read(
                                 char * file_name,
                                 int * z_size_out,
                                 int * bin_num_out,
                                 double ** z_out,
                                 double ** dNdz_out,
                                 double * z_mean,
                                 double * z_min,
                                 double * z_max,
                                 ErrorMsg error_message
                                 ) {

  FILE * input_file;
  /* current line of the file */
  char line[_SELECTION_FILE_LINE_LENGTH_MAX_];
  /* position in the line of the first character not read yet */
  char * left;
  /* position in the line just after the number read by strtod() */
  char * end;
  /* number of numbers in the current line, and in the first line with numbers */
  int value_num;
  int column_num=0;
  /* local copies of the outputs: number of redshift values and of
     bins, redshift values and selection function of each bin */
  int z_size;
  int bin_num;
  double * z;
  double * dNdz;
  /* selection function of the current bin */
  double * dNdz_bin;
  int index_z;
  int bin;
  int index_first,index_last;
  double dz,norm,mean;

  /** - first pass: count the number of lines with values, and check
      that they all have the same number of columns (the file is
      closed before returning in case of error) */

  class_open(input_file,file_name,"r",error_message);

  z_size = 0;

  while (fgets(line,_SELECTION_FILE_LINE_LENGTH_MAX_,input_file) != NULL) {

    /* fgets() stops before the end of a line that does not fit in the buffer */
    class_test_except((strchr(line,'\n') == NULL) && (feof(input_file) == 0),
                      error_message,
                      fclose(input_file),
                      "In file %s, a line is longer than %d characters",
                      file_name,_SELECTION_FILE_LINE_LENGTH_MAX_-1);

    /* skip the white space at the beginning of the line; then skip
       the line if it is a comment (starting with '#') or if it is
       empty */
    left=line;
    while ((left[0]==' ') || (left[0]=='\t'))
      left++;
    if ((left[0]=='#') || (left[0]=='\n') || (left[0]=='\r') || (left[0]=='\0'))
      continue;

    /* count the numbers in the line: each call to strtod() reads the
       number starting at 'left' (after skipping white space), and
       sets 'end' just after it, or at 'left' if there is no number
       left to read */
    value_num=0;
    strtod(left,&end);
    while (end != left) {
      value_num++;
      left = end;
      strtod(left,&end);
    }

    /* all lines must have as many numbers as the first one */
    if (z_size == 0)
      column_num = value_num;

    class_test_except(value_num != column_num,
                      error_message,
                      fclose(input_file),
                      "In file %s, the line '%s' contains %d numbers, while the first line contained %d (a redshift and the selection function of each bin)",
                      file_name,line,value_num,column_num);

    z_size++;
  }

  fclose(input_file);

  class_test(column_num < 2,
             error_message,
             "File %s should contain at least two columns (redshift and selection function of at least one bin)",
             file_name);

  class_test(z_size < 2,
             error_message,
             "File %s should contain at least two lines",
             file_name);

  bin_num = column_num-1;

  class_test(bin_num > _SELECTION_NUM_MAX_,
             error_message,
             "File %s contains the selection functions of %d bins, hence you should increase _SELECTION_NUM_MAX_ in include/perturbations.h to at least this number",
             file_name,bin_num);

  class_alloc(z,   z_size*sizeof(double),         error_message);
  class_alloc(dNdz,bin_num*z_size*sizeof(double), error_message);

  *z_size_out = z_size;
  *bin_num_out = bin_num;
  *z_out = z;
  *dNdz_out = dNdz;

  /** - second pass: read the values (the file is closed before
      returning in case of error) */

  class_open(input_file,file_name,"r",error_message);

  index_z=0;

  while ((index_z < z_size) && (fgets(line,_SELECTION_FILE_LINE_LENGTH_MAX_,input_file) != NULL)) {

    /* skip white space, comments and empty lines as in the first pass */
    left=line;
    while ((left[0]==' ') || (left[0]=='\t'))
      left++;
    if ((left[0]=='#') || (left[0]=='\n') || (left[0]=='\r') || (left[0]=='\0'))
      continue;

    /* read the redshift, and then the selection function of each bin */
    z[index_z] = strtod(left,&end);
    left = end;
    for (bin=0; bin<bin_num; bin++) {
      dNdz[bin*z_size+index_z] = strtod(left,&end);
      left = end;
      class_test_except(dNdz[bin*z_size+index_z] < 0.,
                        error_message,
                        fclose(input_file),
                        "In file %s, the selection function of bin %d is negative at z=%e",
                        file_name,bin,z[index_z]);
    }

    class_test_except(z[index_z] < 0.,
                      error_message,
                      fclose(input_file),
                      "In file %s, the redshift z=%e is negative",
                      file_name,z[index_z]);

    class_test_except((index_z > 0) && (z[index_z] <= z[index_z-1]),
                      error_message,
                      fclose(input_file),
                      "In file %s, the redshifts must be in strictly growing order, but z=%e comes after z=%e",
                      file_name,z[index_z],z[index_z-1]);

    index_z++;
  }

  fclose(input_file);

  class_test(index_z != z_size,
             error_message,
             "File %s changed while being read",
             file_name);

  /** - for each bin, find the redshift range in which the selection
      function is non-zero (including one more point on each side,
      since the interpolated function is non-zero between the last
      zero and the first non-zero value), and the mean redshift */

  for (bin=0; bin<bin_num; bin++) {

    dNdz_bin = dNdz + bin*z_size;

    index_first=0;
    while ((index_first < z_size) && (dNdz_bin[index_first] == 0.))
      index_first++;

    class_test(index_first == z_size,
               error_message,
               "In file %s, the selection function of bin %d vanishes everywhere",
               file_name,bin);

    index_last=z_size-1;
    while (dNdz_bin[index_last] == 0.)
      index_last--;

    z_min[bin] = z[MAX(index_first-1,0)];
    z_max[bin] = z[MIN(index_last+1,z_size-1)];

    /* mean redshift = int z dN/dz dz / int dN/dz dz (trapezoidal rule) */
    norm = 0.;
    mean = 0.;
    for (index_z=0; index_z<z_size-1; index_z++) {
      dz = z[index_z+1]-z[index_z];
      norm += 0.5*(dNdz_bin[index_z]+dNdz_bin[index_z+1])*dz;
      mean += 0.5*(z[index_z]*dNdz_bin[index_z]+z[index_z+1]*dNdz_bin[index_z+1])*dz;
    }
    z_mean[bin] = mean/norm;
  }

  return _SUCCESS_;
}

/**
 * For selection=tabulated, read the selection functions of all bins in
 * the file ptr->selection_file_name, and spline them in order to
 * interpolate them in transfer_selection_function()
 *
 * @param ppt Input: pointer to perturbation structure
 * @param ptr Input/Output: pointer to transfer structure
 * @return the error status
 */

int transfer_tabulated_selection_read(
                                      struct perturbations * ppt,
                                      struct transfer * ptr
                                      ) {

  int bin_num;
  int bin;
  double z_mean[_SELECTION_NUM_MAX_];
  double z_min[_SELECTION_NUM_MAX_];
  double z_max[_SELECTION_NUM_MAX_];

  class_call(transfer_selection_file_read(ptr->selection_file_name,
                                          &(ptr->selection_file_z_size),
                                          &bin_num,
                                          &(ptr->selection_file_z),
                                          &(ptr->selection_file_dNdz),
                                          z_mean,
                                          z_min,
                                          z_max,
                                          ptr->error_message),
             ptr->error_message,
             ptr->error_message);

  class_test(bin_num != ppt->selection_num,
             ptr->error_message,
             "File %s contains %d selection functions, while %d bins were expected",
             ptr->selection_file_name,bin_num,ppt->selection_num);

  class_alloc(ptr->selection_file_dddNdz,bin_num*ptr->selection_file_z_size*sizeof(double),ptr->error_message);

  for (bin=0; bin<bin_num; bin++) {
    class_call(array_spline_table_lines(ptr->selection_file_z,
                                        ptr->selection_file_z_size,
                                        ptr->selection_file_dNdz+bin*ptr->selection_file_z_size,
                                        1,
                                        ptr->selection_file_dddNdz+bin*ptr->selection_file_z_size,
                                        _SPLINE_EST_DERIV_,
                                        ptr->error_message),
               ptr->error_message,
               ptr->error_message);
  }

  return _SUCCESS_;
}

/**
 * This routine allocates and initializes the workspace used for
 * computing the transfer functions at a given wavenumber
 *
 * @param ptr                Input: pointer to transfer structure
 * @param ppt                Input: pointer to perturbation structure
 * @param ptw                Output: pointer to workspace
 * @param K                  Input: curvature parameter
 * @param sgnK               Input: sign of curvature
 * @param tau0_minus_tau_cut Input: critical value of (tau0-tau) in time cut approximation
 * @param pBIS               Input: pointer to structure containing the flat spherical Bessel functions
 * @return the error status
 */

int transfer_workspace_init(
                            struct transfer * ptr,
                            struct perturbations * ppt,
                            struct transfer_workspace *ptw,
                            double K,
                            int sgnK,
                            double tau0_minus_tau_cut,
                            HyperInterpStruct * pBIS){

  int index_md;
  int index_ts;
  int index_tt;
  int index_tt_of_ts;

  ptw->l_size = ptr->l_size_max;
  ptw->HIS_allocated=_FALSE_;
  ptw->pBIS = pBIS;
  ptw->K = K;
  ptw->sgnK = sgnK;
  ptw->tau0_minus_tau_cut = tau0_minus_tau_cut;

  /** - find the dimensions of the arrays of interpolated sources and of transfer sources */

  ptw->ic_size_max = 0;
  ptw->tp_size_max = 0;
  ptw->tt_per_ts_max = 0;

  for (index_md = 0; index_md < ptr->md_size; index_md++) {

    ptw->ic_size_max = MAX(ptw->ic_size_max,ppt->ic_size[index_md]);
    ptw->tp_size_max = MAX(ptw->tp_size_max,ppt->tp_size[index_md]);

    for (index_ts = 0; index_ts < ptr->ts_size[index_md]; index_ts++) {
      index_tt_of_ts = 0;
      for (index_tt = 0; index_tt < ptr->tt_size[index_md]; index_tt++) {
        if (ptr->ts_of_tt[index_md][index_tt] == index_ts)
          index_tt_of_ts++;
      }
      ptw->tt_per_ts_max = MAX(ptw->tt_per_ts_max,index_tt_of_ts);
    }
  }

  /** - allocate arrays (the interpolated sources will only be allocated when needed) */

  class_calloc(ptw->interpolated_sources,ptw->ic_size_max*ptw->tp_size_max,sizeof(double*),ptr->error_message);
  class_alloc(ptw->is_interpolated,   ptw->ic_size_max*ptw->tp_size_max*sizeof(short),                              ptr->error_message);
  class_alloc(ptw->sources,           ptw->ic_size_max*ptw->tt_per_ts_max*ptr->tau_size_max*sizeof(double),         ptr->error_message);
  class_alloc(ptw->chi,               ptr->tau_size_max*sizeof(double),                                             ptr->error_message);
  class_alloc(ptw->cscKgen,           ptr->tau_size_max*sizeof(double),                                             ptr->error_message);
  class_alloc(ptw->cotKgen,           ptr->tau_size_max*sizeof(double),                                             ptr->error_message);
  class_alloc(ptw->radial_function,   ptr->tau_size_max*sizeof(double),                                             ptr->error_message);
  class_alloc(ptw->Phi,               ptr->tau_size_max*sizeof(double),                                             ptr->error_message);
  class_alloc(ptw->dPhi,              ptr->tau_size_max*sizeof(double),                                             ptr->error_message);
  class_alloc(ptw->d2Phi,             ptr->tau_size_max*sizeof(double),                                             ptr->error_message);
  class_alloc(ptw->chireverse,        ptr->tau_size_max*sizeof(double),                                             ptr->error_message);
  class_alloc(ptw->rescale_function,  ptr->tau_size_max*sizeof(double),                                             ptr->error_message);
  class_alloc(ptw->method,            ptw->ic_size_max*ptw->tt_per_ts_max*sizeof(enum transfer_integration_method), ptr->error_message);
  class_alloc(ptw->index_tau_max,     ptw->ic_size_max*ptw->tt_per_ts_max*sizeof(int),                              ptr->error_message);
  class_alloc(ptw->bessel_truncation, ptw->ic_size_max*ptw->tt_per_ts_max*sizeof(short),                            ptr->error_message);
  ptw->bessel_size = 0;

  return _SUCCESS_;
}

/**
 * This routine frees the workspace
 *
 * @param ptr Input: pointer to transfer structure
 * @param ptw Input: pointer to workspace
 * @return the error status
 */

int transfer_workspace_free(
                            struct transfer * ptr,
                            struct transfer_workspace *ptw
                            ) {

  int index;

  if (ptw->HIS_allocated==_TRUE_){
    //Free HIS structure:
    class_call(hyperspherical_HIS_free(&(ptw->HIS),ptr->error_message),
               ptr->error_message,
               ptr->error_message);
  }

  for (index = 0; index < ptw->ic_size_max*ptw->tp_size_max; index++) {
    if (ptw->interpolated_sources[index] != NULL)
      free(ptw->interpolated_sources[index]);
  }
  free(ptw->interpolated_sources);
  free(ptw->is_interpolated);
  free(ptw->sources);
  free(ptw->chi);
  free(ptw->cscKgen);
  free(ptw->cotKgen);
  free(ptw->radial_function);
  free(ptw->Phi);
  free(ptw->dPhi);
  free(ptw->d2Phi);
  free(ptw->chireverse);
  free(ptw->rescale_function);
  free(ptw->method);
  free(ptw->index_tau_max);
  free(ptw->bessel_truncation);

  return _SUCCESS_;
}

/**
 * In a non-flat universe, this routine computes the hyperspherical
 * Bessel functions Phi_l^nu(x) needed at the wavenumber of index
 * index_q, i.e. for nu = q/sqrt(|K|), and stores them in the
 * hyperspherical interpolation structure ptw->HIS of the workspace
 * (after freeing the previous one). In the open case, the largest
 * multipole for which they are needed is found with
 * transfer_get_lmax(). For wavenumbers above
 * ptr->index_q_flat_approximation, and in a flat universe, nothing is
 * done: the flat Bessel functions of ptw->pBIS are used instead.
 *
 * @param ppr     Input: pointer to precision structure
 * @param ptr     Input: pointer to transfer structure
 * @param ptw     Input/Output: pointer to workspace (ptw->HIS updated)
 * @param index_q Input: index of wavenumber
 * @param tau0    Input: conformal age
 * @return the error status
 */

int transfer_update_HIS(
                        struct precision * ppr,
                        struct transfer * ptr,
                        struct transfer_workspace * ptw,
                        int index_q,
                        double tau0
                        ) {

  double nu,new_nu;
  int int_nu;
  double xmin, xmax, sampling, phiminabs, xtol;
  double sqrt_absK;
  int l_size_max;
  int index_l_left,index_l_right;

  if (ptw->HIS_allocated == _TRUE_) {
    class_call(hyperspherical_HIS_free(&(ptw->HIS),ptr->error_message),
               ptr->error_message,
               ptr->error_message);
    ptw->HIS_allocated = _FALSE_;
  }

  if ((ptw->sgnK!=0) && (index_q < ptr->index_q_flat_approximation)) {

    xmin = ppr->hyper_x_min;

    sqrt_absK = sqrt(ptw->sgnK*ptw->K);

    xmax = sqrt_absK*tau0;
    nu = ptr->q[index_q]/sqrt_absK;

    if (ptw->sgnK == 1) {
      xmax = MIN(xmax,_PI_/2.0-ppr->hyper_x_min); //We only need solution on [0;pi/2]

      int_nu = (int)(nu+0.2);
      new_nu = (double)int_nu;
      class_test(nu-new_nu > 1.e-6,
                 ptr->error_message,
                 "problem in q list definition in closed case for index_q=%d, nu=%e, nu-int(nu)=%e",index_q,nu,nu-new_nu);
      nu = new_nu;

    }

    if (nu > ppr->hyper_nu_sampling_step)
      sampling = ppr->hyper_sampling_curved_high_nu;
    else
      sampling = ppr->hyper_sampling_curved_low_nu;

    /* find the highest value of l such that x_nonzero < xmax = sqrt(|K|) tau0. That will be l_max. */
    l_size_max = ptr->l_size_max;
    if (ptw->sgnK == 1)
      while ((double)ptr->l[l_size_max-1] >= nu)
        l_size_max--;

    if (ptw->sgnK == -1){
      xtol = ppr->hyper_x_tol;
      phiminabs = ppr->hyper_phi_min_abs;

      /* First try to find lmax using fast approximation: */
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

      /* Now use WKB approximation to eventually modify borders: */
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

    class_test(nu <= 0.,
               ptr->error_message,
               "nu=%e when index_q=%d, q=%e, K=%e, sqrt(|K|)=%e; instead nu should always be strictly positive",
               nu,index_q,ptr->q[index_q],ptw->K,sqrt_absK);

    class_call(hyperspherical_HIS_create(ptw->sgnK,
                                         nu,
                                         l_size_max,
                                         ptr->l,
                                         xmin,
                                         xmax,
                                         sampling,
                                         ptr->l[l_size_max-1]+1,
                                         ppr->hyper_phi_min_abs,
                                         &(ptw->HIS),
                                         ptr->error_message),
               ptr->error_message,
               ptr->error_message);

    ptw->HIS_allocated = _TRUE_;

  }

  return _SUCCESS_;
}

/**
 * In a non-flat universe, for a given wavenumber nu = q/sqrt(|K|),
 * find the largest multipole for which the hyperspherical Bessel
 * function Phi_l^nu(x) is non-negligible somewhere in the range [0,
 * xmax] where it is needed. Phi_l^nu(x) is negligible for x <
 * x_nonzero(l), the smallest argument at which |Phi_l^nu(x)| reaches
 * phiminabs, which grows with l. Above this multipole, there is no
 * need to compute and store the Bessel functions.
 *
 * The result is the pair of neighboring indices index_l_left,
 * index_l_right = index_l_left+1 such that x_nonzero(lvec[index_l_left])
 * < xmax <= x_nonzero(lvec[index_l_right]). On input, index_l_left and
 * index_l_right are a first guess of this bracket (for instance the
 * result of a previous call with a cheaper approximation for
 * x_nonzero). The bracket is first enlarged until it contains the
 * boundary (searching with steps growing by a factor 5), and then
 * narrowed down by bisection.
 *
 * @param get_xmin_generic Input: function computing x_nonzero(l) for given sgnK, l, nu, xtol, phiminabs (e.g. hyperspherical_get_xmin_from_approx() or hyperspherical_get_xmin_from_Airy())
 * @param sgnK             Input: sign of the curvature
 * @param nu               Input: wavenumber divided by sqrt(|K|)
 * @param lvec             Input: list of multipoles
 * @param lsize            Input: number of multipoles in lvec to be considered
 * @param phiminabs        Input: threshold below which the Bessel functions are negligible
 * @param xmax             Input: largest argument at which the Bessel functions are needed
 * @param xtol             Input: tolerance on x_nonzero
 * @param index_l_left     Input/Output: index in lvec of the largest multipole with x_nonzero < xmax
 * @param index_l_right    Input/Output: index in lvec of the next multipole
 * @param error_message    Output: error message
 * @return the error status
 */

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

  /* smallest argument at which the Bessel function of a given multipole is non-negligible */
  double x_nonzero;
  /* number of evaluations of x_nonzero (required by the signature of get_xmin_generic, not used here) */
  int fevals=0;
  int index_l_mid;
  /* step in the search for the boundaries of the bracket */
  int multiplier;
  /* has index_l_right been checked to be a valid right boundary? */
  int right_boundary_checked = _FALSE_;

  /** - if the Bessel function of the first multipole is already
      negligible in [0, xmax], return index_l_right = lsize-1 (all
      multipoles are then kept by the caller) */

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
    *index_l_right = MAX(lsize-1,1);
    return _SUCCESS_;
  }

  /** - if the Bessel function of the last multipole is still
      non-negligible in [0, xmax], all multipoles are needed */

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
    *index_l_left = MAX(0,(lsize-2));
    return _SUCCESS_;
  }

  /** - otherwise, hunt for a valid left boundary of the bracket (a
      multipole with x_nonzero <= xmax) by moving index_l_left down
      with growing steps. Each multipole found to be negligible on
      the way is a valid right boundary. */

  for (multiplier=1; ;multiplier *= 5){
    class_call(get_xmin_generic(sgnK,
                                lvec[*index_l_left],
                                nu,
                                xtol,
                                phiminabs,
                                &x_nonzero,
                                &fevals),
               error_message,
               error_message);
    if (x_nonzero <= xmax){
      /* left boundary found */
      break;
    }
    else{
      /* the current index_l_left is a valid right boundary */
      *index_l_right = *index_l_left;
      right_boundary_checked = _TRUE_;
    }
    *index_l_left = (*index_l_left)-multiplier;
    if (*index_l_left<=0){
      *index_l_left = 0;
      break;
    }
  }

  /** - if the input index_l_left was already a valid left boundary,
      the input index_l_right has not been checked yet: hunt for a
      valid right boundary (a multipole with x_nonzero >= xmax) by
      moving index_l_right up with growing steps. Each multipole found
      to be non-negligible on the way is a valid left boundary. */

  if (right_boundary_checked == _FALSE_){
    for (multiplier=1; ;multiplier *= 5){
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
        /* right boundary found */
        break;
      }
      else{
        /* the current index_l_right is a valid left boundary */
        *index_l_left = *index_l_right;
      }
      *index_l_right = (*index_l_right)+multiplier;
      if (*index_l_right>=(lsize-1)){
        *index_l_right = lsize-1;
        break;
      }
    }
  }

  /** - narrow down the bracket by bisection until index_l_right = index_l_left+1 */

  while (((*index_l_right) - (*index_l_left)) > 1) {
    index_l_mid= (int)(0.5*((*index_l_right)+(*index_l_left)));
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

  return _SUCCESS_;
}

/**
 * Here we precompute once for all the window functions W(tau) by
 * which the perturbation source S(k,tau) of each transfer type must be
 * multiplied in the line-of-sight integral, on the time sampling of
 * this type defined in transfer_time_samplings():
 *
 * - for CMB temperature and polarisation types, W=1.
 *
 * - for the CMB lensing potential, W is the lensing kernel, with
 *   support [tau_rec, tau0], possibly rescaled by hand by the factor
 *   lcmb_rescale.
 *
 * - for each type of nCl/dCl/sCl we combine the selection function
 *   with the corresponding prefactor (e.g. 1/aH), and, if required,
 *   we also integrate for integrated (lensed) contributions (In the
 *   original ClassGAL paper these would be labeled g4,g5, and lens)
 *
 * All factors of k are added later, in transfer_sources().
 *
 * @param ppr                   Input: pointer to precision structure
 * @param pba                   Input: pointer to background structure
 * @param ppt                   Input: pointer to perturbation structure
 * @param ptr                   Input/Output: pointer to transfer structure
 * @param tau_rec               Input: recombination time
 * @return the error status
 */

int transfer_precompute_windows(
                                struct precision * ppr,
                                struct background * pba,
                                struct perturbations * ppt,
                                struct transfer * ptr,
                                double tau_rec
                                ){
  /** Summary: */

  /** - define local variables */

  /* running indices on modes, types, time samplings */
  int index_md;
  int index_tt;
  int index_ts;

  /* time sampling of the current type */
  double * tau0_minus_tau;
  double * w_trapz;

  /* window function of the current type */
  double * window;

  /* index running on time */
  int index_tau;

  /* bin for computation of cl_density */
  int bin=0;

  /* number of tau values */
  int tau_size;

  /* for calling background_at_eta */
  int last_index;
  double * pvecback = NULL;

  /* conformal time */
  double tau, tau0;

  /* geometrical quantities */
  double sinKgen_source=0.;
  double sinKgen_source_to_lens=0.;
  double cotKgen_source=0.;
  double cscKgen_lens=0.;

  /* rescaling factor depending on the background at a given time */
  double rescaling=0.;

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

  /* source evolution factor */
  double f_evo = 0.;

  /* allocate temporary arrays for storing selections, and for calling background */
  class_alloc(selection,MAX(ptr->tau_size_max,(int)ppr->selection_sampling)*sizeof(double),ptr->error_message);
  class_alloc(pvecback,pba->bg_size*sizeof(double),ptr->error_message);

  /* conformal time today */
  tau0 = pba->conformal_age;

  class_alloc(ptr->window,ptr->md_size*sizeof(double**),ptr->error_message);

  for (index_md = 0; index_md < ptr->md_size; index_md++) {

    class_alloc(ptr->window[index_md],ptr->tt_size[index_md]*sizeof(double*),ptr->error_message);

    /* Loop through different types to be precomputed */
    for (index_tt = 0; index_tt < ptr->tt_size[index_md]; index_tt++) {

      /* time sampling of this type */
      index_ts = ptr->ts_of_tt[index_md][index_tt];
      tau_size = ptr->tau_size_ts[index_md][index_ts];
      tau0_minus_tau = ptr->tau0_minus_tau[index_md][index_ts];
      w_trapz = ptr->w_trapz[index_md][index_ts];

      /* redshift bin of this type (if any) */
      bin = ptr->bin_of_tt[index_md][index_tt];

      class_alloc(ptr->window[index_md][index_tt],tau_size*sizeof(double),ptr->error_message);
      window = ptr->window[index_md][index_tt];

      switch (ptr->window_type[index_md][index_tt]) {

        /** - CMB temperature and polarisation: no window */

      case window_none:

        for (index_tau = 0; index_tau < tau_size; index_tau++) {
          window[index_tau] = 1.;
        }
        break;

        /** - CMB lensing potential */

      case window_cmb_lensing:

        for (index_tau = 0; index_tau < tau_size; index_tau++) {

          /* conformal time */
          tau = tau0 - tau0_minus_tau[index_tau];

          /* lensing source =  - W(tau) (phi(k,tau) + psi(k,tau)) Heaviside(tau-tau_rec)
             with
             psi,phi = metric perturbation in newtonian gauge (phi+psi = Phi_A-Phi_H of Bardeen)
             W = (tau-tau_rec)/(tau_0-tau)/(tau_0-tau_rec)
             H(x) = Heaviside
             (in tau = tau_0, set source = 0 to avoid division by zero;
             regulated anyway by Bessel).
          */

          if (index_tau == tau_size-1) {
            rescaling=0.;
          }
          else {
            switch (pba->sgnK){
            case 1:
              rescaling = sqrt(pba->K)
                *sin((tau_rec-tau)*sqrt(pba->K))
                /sin((tau0-tau)*sqrt(pba->K))
                /sin((tau0-tau_rec)*sqrt(pba->K));
              break;
            case 0:
              rescaling = (tau_rec-tau)/(tau0-tau)/(tau0-tau_rec);
              break;
            case -1:
              rescaling = sqrt(-pba->K)
                *sinh((tau_rec-tau)*sqrt(-pba->K))
                /sinh((tau0-tau)*sqrt(-pba->K))
                /sinh((tau0-tau_rec)*sqrt(-pba->K));
              break;
            }
            // Note: until 2.4.3 there was a bug here: the curvature effects had been omitted.
          }

          window[index_tau] = rescaling * ptr->lcmb_rescale;
        }
        break;

        /** - non-integrated (local) number count contributions */

      case window_selection:

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

          /* geometrical quantity */
          switch (pba->sgnK){
          case 1:
            cotKgen_source = sqrt(pba->K)
              *cos(tau0_minus_tau[index_tau]*sqrt(pba->K))
              /sin(tau0_minus_tau[index_tau]*sqrt(pba->K));
            break;
          case 0:
            cotKgen_source = 1./(tau0_minus_tau[index_tau]);
            break;
          case -1:
            cotKgen_source = sqrt(-pba->K)
              *cosh(tau0_minus_tau[index_tau]*sqrt(-pba->K))
              /sinh(tau0_minus_tau[index_tau]*sqrt(-pba->K));
            break;
          }

          /* corresponding background quantities */
          class_call(background_at_tau(pba,
                                       tau,
                                       long_info,
                                       inter_normal,
                                       &last_index,
                                       pvecback),
                     pba->error_message,
                     ptr->error_message);

          /* Source evolution, used by nCl doppler and nCl gravity terms */

          if ((_index_tt_in_range_(ptr->index_tt_d0,      ppt->selection_num, ppt->has_nc_rsd)) ||
              (_index_tt_in_range_(ptr->index_tt_d1,      ppt->selection_num, ppt->has_nc_rsd)) ||
              (_index_tt_in_range_(ptr->index_tt_nc_g2,   ppt->selection_num, ppt->has_nc_gr)))
            {
              class_call(transfer_f_evo(pba,ptr,pvecback,last_index,cotKgen_source,&f_evo),
                         ptr->error_message,
                         ptr->error_message);
              /* Error in old CLASS 2.6.3 : Number count evolution did not respect curvature */
            }

          /* matter density source =  [- (dz/dtau) W(z)] * delta_m(k,tau)
             = W(tau) delta_m(k,tau)
             with
             delta_m = total matter perturbation (defined in gauge-independent way, see arXiv 1307.1459)
             W(z) = redshift space selection function = dN/dz
             W(tau) = same wrt conformal time = dN/dtau
             (in tau = tau_0, set source = 0 to avoid division by zero;
             regulated anyway by Bessel).
          */

          if (_index_tt_in_range_(ptr->index_tt_density, ppt->selection_num, ppt->has_nc_density))
            rescaling = ptr->selection_bias[bin]*selection[index_tau];

          /* redshift space distortion source = - [- (dz/dtau) W(z)] * (k/H) * theta(k,tau) */

          if (_index_tt_in_range_(ptr->index_tt_rsd,     ppt->selection_num, ppt->has_nc_rsd))
            rescaling = selection[index_tau]/pvecback[pba->index_bg_H]/pvecback[pba->index_bg_a];

          if (_index_tt_in_range_(ptr->index_tt_d0,      ppt->selection_num, ppt->has_nc_rsd))
            rescaling = (f_evo-3.)*selection[index_tau]*pvecback[pba->index_bg_H]*pvecback[pba->index_bg_a];

          if (_index_tt_in_range_(ptr->index_tt_d1,      ppt->selection_num, ppt->has_nc_rsd))

            rescaling = selection[index_tau]*(1.
                                              +pvecback[pba->index_bg_H_prime]
                                              /pvecback[pba->index_bg_a]
                                              /pvecback[pba->index_bg_H]
                                              /pvecback[pba->index_bg_H]
                                              +(2.-5.*ptr->selection_magnification_bias[bin])
                                              // /tau0_minus_tau[index_tau] // in flat space
                                              *cotKgen_source  // in general case
                                              /pvecback[pba->index_bg_a]
                                              /pvecback[pba->index_bg_H]
                                              +5.*ptr->selection_magnification_bias[bin]
                                              -f_evo
                                              );

          if (_index_tt_in_range_(ptr->index_tt_nc_g1,   ppt->selection_num, ppt->has_nc_gr))

            rescaling = selection[index_tau];

          if (_index_tt_in_range_(ptr->index_tt_nc_g2,   ppt->selection_num, ppt->has_nc_gr))

            rescaling = -selection[index_tau]*(3.
                                               +pvecback[pba->index_bg_H_prime]
                                               /pvecback[pba->index_bg_a]
                                               /pvecback[pba->index_bg_H]
                                               /pvecback[pba->index_bg_H]
                                               +(2.-5.*ptr->selection_magnification_bias[bin])
                                               // /tau0_minus_tau[index_tau]  // in flat space
                                               *cotKgen_source  // in general case
                                               /pvecback[pba->index_bg_a]
                                               /pvecback[pba->index_bg_H]
                                               -f_evo
                                               );

          if (_index_tt_in_range_(ptr->index_tt_nc_g3,   ppt->selection_num, ppt->has_nc_gr))
            rescaling = selection[index_tau]/pvecback[pba->index_bg_a]/pvecback[pba->index_bg_H];

          /* finally store in array */
          window[index_tau] = rescaling;
        }
        break;

        /** - integrated contributions */

      case window_integrated:

        /* dirac case */
        if (ppt->selection == dirac) {
          tau_sources_size=1;
        }
        /* other cases (gaussian, tophat...) */
        else {
          tau_sources_size=ppr->selection_sampling;
        }

        class_alloc(tau0_minus_tau_lensing_sources,
                    tau_sources_size*sizeof(double),
                    ptr->error_message);

        class_alloc(w_trapz_lensing_sources,
                    tau_sources_size*sizeof(double),
                    ptr->error_message);

        /* time sampling for source selection function */
        class_call(transfer_uniform_sampling(tau0,
                                             ptr->selection_tau_min[bin],
                                             ptr->selection_tau_max[bin],
                                             tau_sources_size,
                                             tau0_minus_tau_lensing_sources),
                   ptr->error_message,
                   ptr->error_message);

        /* Compute trapezoidal weights for integration over tau */
        class_call(array_trapezoidal_mweights(tau0_minus_tau_lensing_sources,
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

        /* loop over time and rescale */
        for (index_tau = 0; index_tau < tau_size; index_tau++) {

          /* lensing source =  - W(tau) (phi(k,tau) + psi(k,tau)) Heaviside(tau-tau_rec)
             with
             psi,phi = metric perturbation in newtonian gauge (phi+psi = Phi_A-Phi_H of Bardeen)
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

              switch (pba->sgnK){
              case 1:
                sinKgen_source = sin(tau0_minus_tau_lensing_sources[index_tau_sources]*sqrt(pba->K))/sqrt(pba->K);
                sinKgen_source_to_lens = sin((tau0_minus_tau[index_tau]-tau0_minus_tau_lensing_sources[index_tau_sources])*sqrt(pba->K))/sqrt(pba->K);
                cotKgen_source = cos(tau0_minus_tau_lensing_sources[index_tau_sources]*sqrt(pba->K))/sinKgen_source;
                cscKgen_lens = sqrt(pba->K)/sin(sqrt(pba->K)*tau0_minus_tau[index_tau]);
                break;
              case 0:
                sinKgen_source = tau0_minus_tau_lensing_sources[index_tau_sources];
                sinKgen_source_to_lens = (tau0_minus_tau[index_tau]-tau0_minus_tau_lensing_sources[index_tau_sources]);
                cotKgen_source = 1./(tau0_minus_tau_lensing_sources[index_tau_sources]);
                cscKgen_lens = 1./(tau0_minus_tau[index_tau]);
                break;
              case -1:
                sinKgen_source = sinh(tau0_minus_tau_lensing_sources[index_tau_sources]*sqrt(-pba->K))/sqrt(-pba->K);
                sinKgen_source_to_lens = sinh((tau0_minus_tau[index_tau]-tau0_minus_tau_lensing_sources[index_tau_sources])*sqrt(-pba->K))/sqrt(-pba->K);
                cotKgen_source = cosh(tau0_minus_tau_lensing_sources[index_tau_sources]*sqrt(-pba->K))/sinKgen_source;
                cscKgen_lens = sqrt(-pba->K)/sinh(sqrt(-pba->K)*tau0_minus_tau[index_tau]);
                break;
              }

              /* condition for excluding from the sum the sources located in z=zero */
              if ((tau0_minus_tau_lensing_sources[index_tau_sources] > 0.) && (tau0_minus_tau_lensing_sources[index_tau_sources]-tau0_minus_tau[index_tau] > 0.)) {

                if (_index_tt_in_range_(ptr->index_tt_lensing, ppt->selection_num, ppt->has_cl_lensing_potential)) {

                  rescaling +=
                    //  *(tau0_minus_tau[index_tau]-tau0_minus_tau_lensing_sources[index_tau_sources])
                    //  /tau0_minus_tau[index_tau]
                    //  /tau0_minus_tau_lensing_sources[index_tau_sources]
                    sinKgen_source_to_lens
                    *cscKgen_lens
                    /sinKgen_source
                    * selection[index_tau_sources]
                    * w_trapz_lensing_sources[index_tau_sources];
                }

                if (_index_tt_in_range_(ptr->index_tt_nc_lens, ppt->selection_num, ppt->has_nc_lens)) {

                  rescaling -=
                    (2.-5.*ptr->selection_magnification_bias[bin])/2.
                    //  *(tau0_minus_tau[index_tau]-tau0_minus_tau_lensing_sources[index_tau_sources])
                    //  /tau0_minus_tau[index_tau]
                    //  /tau0_minus_tau_lensing_sources[index_tau_sources]
                    *sinKgen_source_to_lens
                    *cscKgen_lens
                    /sinKgen_source
                    * selection[index_tau_sources]
                    * w_trapz_lensing_sources[index_tau_sources];
                }

                if (_index_tt_in_range_(ptr->index_tt_nc_g4, ppt->selection_num, ppt->has_nc_gr)) {

                  rescaling +=
                    (2.-5.*ptr->selection_magnification_bias[bin])
                    // /tau0_minus_tau_lensing_sources[index_tau_sources]
                    * cotKgen_source
                    * selection[index_tau_sources]
                    * w_trapz_lensing_sources[index_tau_sources];

                }

                if (_index_tt_in_range_(ptr->index_tt_nc_g5, ppt->selection_num, ppt->has_nc_gr)) {

                  /* background quantities at time tau_lensing_source */

                  class_call(background_at_tau(pba,
                                               tau0-tau0_minus_tau_lensing_sources[index_tau_sources],
                                               long_info,
                                               inter_normal,
                                               &last_index,
                                               pvecback),
                             pba->error_message,
                             ptr->error_message);

                  /* Source evolution at time tau_lensing_source */

                  class_call(transfer_f_evo(pba,ptr,pvecback,last_index,cotKgen_source,&f_evo),
                             ptr->error_message,
                             ptr->error_message);


                  rescaling +=
                    (1.
                     + pvecback[pba->index_bg_H_prime]
                     /pvecback[pba->index_bg_a]
                     /pvecback[pba->index_bg_H]
                     /pvecback[pba->index_bg_H]
                     + (2.-5.*ptr->selection_magnification_bias[bin])
                     //  /tau0_minus_tau_lensing_sources[index_tau_sources]
                     * cotKgen_source
                     /pvecback[pba->index_bg_a]
                     /pvecback[pba->index_bg_H]
                     + 5.*ptr->selection_magnification_bias[bin]
                     - f_evo)
                    * selection[index_tau_sources]
                    * w_trapz_lensing_sources[index_tau_sources];
                }
              }
            }
          }

          /* Finally store integrated result for later use */
          window[index_tau] = rescaling;
        }

        /* deallocate temporary arrays */
        free(tau0_minus_tau_lensing_sources);
        free(w_trapz_lensing_sources);
        break;
      }
    }
  }

  /* deallocate temporary arrays */
  free(selection);
  free(pvecback);
  return _SUCCESS_;
}

/**
 * This routine computes the factor f_evo by which the selection
 * function is multiplied in the Doppler (d0, d1) and gravitational
 * (g2) contributions to the number counts, at a given time:
 *
 * f_evo = 2/(a H) cotKgen + H'/(a H^2) - (1/a) d ln(dN/dz)/dz
 *
 * where H' is the derivative of H with respect to conformal time and
 * dN/dz is the evolution of the number of sources (from the file
 * passed in dNdz_evolution or from an analytic form). Without such
 * an evolution, f_evo is set to zero.
 *
 * @param pba        Input: pointer to background structure
 * @param ptr        Input: pointer to transfer structure
 * @param pvecback   Input: background quantities at the current time
 * @param last_index Input: index in a table, only used locally by the spline interpolation (its input value does not matter)
 * @param cotKgen    Input: sqrt(|K|) cot(sqrt(|K|)(tau0-tau)) (with cot replaced by coth if K<0, and 1/(tau0-tau) if K=0) at the current time
 * @param f_evo      Output: value of f_evo
 * @return the error status
 */

int transfer_f_evo(
                   struct background* pba,
                   struct transfer * ptr,
                   double* pvecback,
                   int last_index,
                   double cotKgen, /* Should be FILLED with values of corresponding time */
                   double* f_evo
                   ){
  /* Allocate temporary variables for calculation of f_evo */
  double z;
  double dNdz;
  double dln_dNdz_dz;
  double temp_f_evo;

  if ((ptr->has_nz_evo_file == _TRUE_) || (ptr->has_nz_evo_analytic == _TRUE_)){

    temp_f_evo = 2./pvecback[pba->index_bg_H]/pvecback[pba->index_bg_a]*cotKgen
      + pvecback[pba->index_bg_H_prime]/pvecback[pba->index_bg_H]/pvecback[pba->index_bg_H]/pvecback[pba->index_bg_a];

    /* z = a_0/a-1 (remember that a in the code is in fact a/a_0) */
    z = 1./pvecback[pba->index_bg_a]-1.;

    if (ptr->has_nz_evo_file ==_TRUE_) {

      class_test((z<ptr->nz_evo_z[0]) || (z>ptr->nz_evo_z[ptr->nz_evo_size-1]),
                 ptr->error_message,
                 "Your input file for the selection function only covers the redshift range [%f : %f]. However, your input for the selection function requires z=%f",
                 ptr->nz_evo_z[0],
                 ptr->nz_evo_z[ptr->nz_evo_size-1],
                 z);


      class_call(array_interpolate_spline(
                                          ptr->nz_evo_z,
                                          ptr->nz_evo_size,
                                          ptr->nz_evo_dlog_nz,
                                          ptr->nz_evo_dd_dlog_nz,
                                          1,
                                          z,
                                          &last_index,
                                          &dln_dNdz_dz,
                                          1,
                                          ptr->error_message),
                 ptr->error_message,
                 ptr->error_message);

    }
    else {

      class_call(transfer_dNdz_analytic(ptr,
                                        z,
                                        &dNdz,
                                        &dln_dNdz_dz),
                 ptr->error_message,
                 ptr->error_message);
    }

    temp_f_evo -= dln_dNdz_dz/pvecback[pba->index_bg_a];
  }
  else {
    temp_f_evo = 0.;
  }

  /* after obtaining f_evo, store it in output */
  *f_evo = temp_f_evo;

  return _SUCCESS_;
}
