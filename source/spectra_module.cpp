/** @file spectra.c Documented spectra module
 *
 * Julien Lesgourgues, 1.11.2019
 *
 * This module computes the harmonic power spectra \f$ C_l^{X} \f$'s
 * given the transfer functions and the primordial spectra.
 *
 * The following functions can be called from other modules:
 *
 * -# spectra_init() at the beginning (but after transfer_init())
 * -# spectra_cl_at_l() at any time for computing individual \f$ C_l \f$'s at any l
 * -# spectra_free() at the end
 */

#include "perturbations_module.h"
#include "primordial_module.h"
#include "nonlinear_module.h"
#include "transfer_module.h"
#include "spectra_module.h"
#include "thread_pool.h"

SpectraModule::SpectraModule(InputModulePtr input_module, PerturbationsModulePtr perturbations_module, PrimordialModulePtr primordial_module, NonlinearModulePtr nonlinear_module, TransferModulePtr transfer_module)
: BaseModule(std::move(input_module))
, perturbations_module_(std::move(perturbations_module))
, primordial_module_(std::move(primordial_module))
, nonlinear_module_(std::move(nonlinear_module))
, transfer_module_(std::move(transfer_module)) {
  if (spectra_init() != _SUCCESS_) {
    throw std::runtime_error(error_message_);
  }
}

SpectraModule::~SpectraModule() {
  spectra_free();
}

std::map<std::string, int> SpectraModule::cl_output_index_map() const {

  std::map<std::string, int> index_map;
  if (has_tt_) index_map["tt"] = index_ct_tt_;
  if (has_ee_) index_map["ee"] = index_ct_ee_;
  if (has_te_) index_map["te"] = index_ct_te_;
  if (has_bb_) index_map["bb"] = index_ct_bb_;
  if (has_pp_) index_map["pp"] = index_ct_pp_;
  if (has_tp_) index_map["tp"] = index_ct_tp_;
  if (has_ep_) index_map["ep"] = index_ct_ep_;

  int index = static_cast<int>(index_map.size());
  if (has_dd_ == _TRUE_) {
    for (int index_d1 = 0; index_d1 < d_size_; index_d1++) {
      for (int index_d2 = index_d1; index_d2 <= MIN(index_d1 + psp->non_diag, d_size_ - 1); index_d2++) {
        std::string key = "dens[" + std::to_string(index_d1 + 1) + "]-dens[" + std::to_string(index_d2 + 1) + "]";
        index_map[key] = index++;
      }
    }
  }
  if (has_td_ == _TRUE_) {
    for (int index_d1 = 0; index_d1 < d_size_; index_d1++) {
      std::string key = "T-dens[" + std::to_string(index_d1 + 1) + "]";
      index_map[key] = index++;
    }
  }
  if (has_pd_ == _TRUE_) {
    for (int index_d1 = 0; index_d1 < d_size_; index_d1++) {
      std::string key = "phi-dens[" + std::to_string(index_d1 + 1) + "]";
      index_map[key] = index++;
    }
  }
  if (has_ll_ == _TRUE_) {
    for (int index_d1 = 0; index_d1 < d_size_; index_d1++) {
      for (int index_d2 = index_d1; index_d2 <= MIN(index_d1 + psp->non_diag, d_size_ - 1); index_d2++){
        std::string key = "lens[" + std::to_string(index_d1 + 1) + "]-lens[" + std::to_string(index_d2 + 1) + "]";
        index_map[key] = index++;
      }
    }
  }
  if (has_tl_ == _TRUE_) {
    for (int index_d1 = 0; index_d1 < d_size_; index_d1++) {
      std::string key = "T-lens[" + std::to_string(index_d1 + 1) + "]";
      index_map[key] = index++;
    }
  }
  if (has_dl_ == _TRUE_) {
    for (int index_d1 = 0; index_d1 < d_size_; index_d1++) {
      for (int index_d2 = MAX(index_d1-psp->non_diag, 0); index_d2 <= MIN(index_d1 + psp->non_diag, d_size_ - 1); index_d2++) {
        std::string key = "dens[" + std::to_string(index_d1 + 1) + "]-lens[" + std::to_string(index_d2 + 1) + "]";
        index_map[key] = index++;
      }
    }
  }
  return index_map;
}

void SpectraModule::cl_output_no_copy(int lmax, std::vector<double*>& output_pointers) const {

  ThrowRuntimeErrorIf((lmax > l_max_tot_) || (lmax < 0), "Error: lmax = %d is outside the allowed range [0, %d]\n", lmax, l_max_tot_);
  ThrowRuntimeErrorIf(output_pointers.size() != ct_size_, "Error: Size of input vector (%d) does not match ct_size = %d\n", output_pointers.size(), ct_size_);

 cl_output_index_map();

  double* cl_md[md_size_];
  double cl_md_data[md_size_][ct_size_];
  for (int index_md = 0; index_md < md_size_; index_md++) {
    cl_md[index_md] = &(cl_md_data[index_md][0]);
  }

  int cl_md_ic_size = 0;
  if (md_size_ > 1) {
    for (int index_md = 0; index_md < md_size_; index_md++) {
      if (perturbations_module_->ic_size_[index_md] > 1) {
        cl_md_ic_size += ic_ic_size_[index_md]*ct_size_;
      }
    }
  }
  std::vector<double> cl_md_ic_data(cl_md_ic_size, 0.0);
  double* cl_md_ic[md_size_];
  int cl_md_ic_index = 0;
  for (int index_md = 0; index_md < md_size_; index_md++) {
    cl_md_ic[index_md] = &cl_md_ic_data[cl_md_ic_index];
    if (perturbations_module_->ic_size_[index_md] > 1) {
      cl_md_ic_index += ic_ic_size_[index_md]*ct_size_;
    }
  }

  double cl_tot[ct_size_];
  for (int l = 0; l <= lmax; l++) {
    if (l < 2) {
      for (auto& output_pointer : output_pointers) {
        output_pointer[l] = 0.0;
      }
    }
    else {
      int status = spectra_cl_at_l(l, cl_tot, cl_md, cl_md_ic);
      ThrowRuntimeErrorIf(status != _SUCCESS_, "Error in SpectraModule::cl_output: %s", error_message_);
      for (int index_ct = 0; index_ct < ct_size_; ++index_ct) {
        output_pointers[index_ct][l] = cl_tot[index_ct];
      }
    }
  }
}

std::map<std::string, std::vector<double>> SpectraModule::cl_output(int lmax) const {
  ThrowRuntimeErrorIf(ppt->has_cls == _FALSE_, "Error: Cls have not been computed! lmax = %d\n", lmax);
  ThrowRuntimeErrorIf((lmax > l_max_tot_) || (lmax < 0), "Error: lmax = %d is outside the allowed range [0, %d]\n", lmax, l_max_tot_);
  std::map<std::string, int> index_map = cl_output_index_map();

  double* cl_md[md_size_];
  double cl_md_data[md_size_][ct_size_];
  for (int index_md = 0; index_md < md_size_; index_md++) {
    cl_md[index_md] = &(cl_md_data[index_md][0]);
  }

  int cl_md_ic_size = 0;
  if (md_size_ > 0) {
    for (int index_md = 0; index_md < md_size_; index_md++) {
      if (perturbations_module_->ic_size_[index_md] > 1) {
        cl_md_ic_size += ic_ic_size_[index_md]*ct_size_;
      }
    }
  }
  double* cl_md_ic[md_size_];
  std::vector<double> cl_md_ic_data;
  if (cl_md_ic_size > 0) {
    cl_md_ic_data.resize(cl_md_ic_size);
    int cl_md_ic_index = 0;
    for (int index_md = 0; index_md < md_size_; index_md++) {
      cl_md_ic[index_md] = &cl_md_ic_data[cl_md_ic_index];
      if (perturbations_module_->ic_size_[index_md] > 1) {
        cl_md_ic_index += ic_ic_size_[index_md]*ct_size_;
      }
    }
  }

  std::vector<double> data_vectors[ct_size_];
  for (int i = 0; i < ct_size_; ++i) {
    data_vectors[i] = std::vector<double>(lmax + 1, 0.0);
  }
  double cl_tot[ct_size_];
  for (int l = 2; l <= lmax; l++) {
    int status = spectra_cl_at_l(l, cl_tot, cl_md, cl_md_ic);
    ThrowRuntimeErrorIf(status != _SUCCESS_, "Error in SpectraModule::cl_output: %s", error_message_);
    for (int i = 0; i < ct_size_; ++i) {
      data_vectors[i][l] = cl_tot[i];
    }
  }
  // Now move vectors into map. We could have created the vectors inside the map directly, but that would
  // lead to many unnecessary map-lookups in the l-loop above.
  std::map<std::string, std::vector<double>> output;
  for (const auto& element : index_map) {
    output[element.first] = std::move(data_vectors[element.second]);
  }

  return output;
}


/**
 * Anisotropy power spectra \f$ C_l\f$'s for all types, modes and initial conditions.
 *
 * This routine evaluates all the \f$C_l\f$'s at a given value of l by
 * interpolating in the pre-computed table. When relevant, it also
 * sums over all initial conditions for each mode, and over all modes.
 *
 * This function can be
 * called from whatever module at whatever time, provided that
 * spectra_init() has been called before, and spectra_free() has not
 * been called yet.
 *
 * @param psp        Input: pointer to spectra structure (containing pre-computed table)
 * @param l          Input: multipole number
 * @param cl_tot     Output: total \f$C_l\f$'s for all types (TT, TE, EE, etc..)
 * @param cl_md      Output: \f$C_l\f$'s for all types (TT, TE, EE, etc..) decomposed mode by mode (scalar, tensor, ...) when relevant
 * @param cl_md_ic   Output: \f$C_l\f$'s for all types (TT, TE, EE, etc..) decomposed by pairs of initial conditions (adiabatic, isocurvatures) for each mode (usually, only for the scalar mode) when relevant
 * @return the error status
 */

int SpectraModule::spectra_cl_at_l(double l,
                    double * cl_tot,    /* array with argument cl_tot[index_ct] (must be already allocated) */
                    double * * cl_md,   /* array with argument cl_md[index_md][index_ct] (must be already allocated only if several modes) */
                    double * * cl_md_ic /* array with argument cl_md_ic[index_md][index_ic1_ic2*ct_size_+index_ct] (must be already allocated for a given mode only if several ic's) */
                    ) const {

  /** Summary: */

  /** - define local variables */

  int last_index;
  int index_md;
  int index_ic1,index_ic2,index_ic1_ic2;
  int index_ct;

  /** - (a) treat case in which there is only one mode and one initial condition.
      Then, only cl_tot needs to be filled. */

  if ((md_size_ == 1) && (ic_size_[0] == 1)) {
    index_md = 0;
    if ((int)l <= l_[l_size_[index_md]-1]) {

      /* interpolate at l */
      class_call(array_interpolate_spline(l_,
                                          l_size_[index_md],
                                          cl_[index_md],
                                          ddcl_[index_md],
                                          ct_size_,
                                          l,
                                          &last_index,
                                          cl_tot,
                                          ct_size_,
                                          error_message_),
                 error_message_,
                 error_message_);

      /* set to zero for the types such that l<l_max */
      for (index_ct = 0; index_ct < ct_size_; index_ct++)
        if ((int)l > l_max_ct_[index_md][index_ct])
          cl_tot[index_ct]=0.;
    }
    else {
      for (index_ct=0; index_ct < ct_size_; index_ct++)
        cl_tot[index_ct]=0.;
    }
  }

  /** - (b) treat case in which there is only one mode
      with several initial condition.
      Fill cl_md_ic[index_md=0] and sum it to get cl_tot. */

  if ((md_size_ == 1) && (ic_size_[0] > 1)) {
    index_md = 0;
    for (index_ct = 0; index_ct < ct_size_; index_ct++)
      cl_tot[index_ct]=0.;
    for (index_ic1 = 0; index_ic1 < ic_size_[index_md]; index_ic1++) {
      for (index_ic2 = index_ic1; index_ic2 < ic_size_[index_md]; index_ic2++) {
        index_ic1_ic2 = index_symmetric_matrix(index_ic1, index_ic2, ic_size_[index_md]);
        if (((int)l <= l_[l_size_[index_md] - 1]) &&
            (is_non_zero_[index_md][index_ic1_ic2] == _TRUE_)) {

          class_call(array_interpolate_spline(l_,
                                              l_size_[index_md],
                                              cl_[index_md],
                                              ddcl_[index_md],
                                              ic_ic_size_[index_md]*ct_size_,
                                              l,
                                              &last_index,
                                              cl_md_ic[index_md],
                                              ic_ic_size_[index_md]*ct_size_,
                                              error_message_),
                     error_message_,
                     error_message_);

          for (index_ct = 0; index_ct < ct_size_; index_ct++)
            if ((int)l > l_max_ct_[index_md][index_ct])
              cl_md_ic[index_md][index_ic1_ic2*ct_size_ + index_ct] = 0.;
        }
        else {
          for (index_ct = 0; index_ct < ct_size_; index_ct++)
            cl_md_ic[index_md][index_ic1_ic2*ct_size_ + index_ct] = 0.;
        }

        /* compute cl_tot by summing over cl_md_ic */
        for (index_ct = 0; index_ct < ct_size_; index_ct++) {
          if (index_ic1 == index_ic2)
            cl_tot[index_ct] += cl_md_ic[index_md][index_ic1_ic2*ct_size_ + index_ct];
          else
            cl_tot[index_ct] += 2.*cl_md_ic[index_md][index_ic1_ic2*ct_size_ + index_ct];
        }
      }
    }
  }

  /** - (c) loop over modes */

  if (md_size_ > 1) {

    for (index_ct = 0; index_ct < ct_size_; index_ct++)
      cl_tot[index_ct]=0.;

    for (index_md = 0; index_md < md_size_; index_md++) {

      /** - --> (c.1.) treat case in which the mode under consideration
          has only one initial condition.
          Fill cl_md[index_md]. */

      if (ic_size_[index_md] == 1) {
        if ((int)l <= l_[l_size_[index_md] - 1]) {

          class_call(array_interpolate_spline(l_,
                                              l_size_[index_md],
                                              cl_[index_md],
                                              ddcl_[index_md],
                                              ct_size_,
                                              l,
                                              &last_index,
                                              cl_md[index_md],
                                              ct_size_,
                                              error_message_),
                     error_message_,
                     error_message_);

          for (index_ct = 0; index_ct < ct_size_; index_ct++)
            if ((int)l > l_max_ct_[index_md][index_ct])
              cl_md[index_md][index_ct]=0.;
        }
        else {
          for (index_ct = 0; index_ct < ct_size_; index_ct++)
            cl_md[index_md][index_ct]=0.;
        }
      }

      /** - --> (c.2.) treat case in which the mode under consideration
          has several initial conditions.
          Fill cl_md_ic[index_md] and sum it to get cl_md[index_md] */

      if (ic_size_[index_md] > 1) {

        if ((int)l <= l_[l_size_[index_md] - 1]) {

          /* interpolate all ic and ct */
          class_call(array_interpolate_spline(l_,
                                              l_size_[index_md],
                                              cl_[index_md],
                                              ddcl_[index_md],
                                              ic_ic_size_[index_md]*ct_size_,
                                              l,
                                              &last_index,
                                              cl_md_ic[index_md],
                                              ic_ic_size_[index_md]*ct_size_,
                                              error_message_),
                     error_message_,
                     error_message_);

          /* set to zero some of the components */
          for (index_ic1 = 0; index_ic1 < ic_size_[index_md]; index_ic1++) {
            for (index_ic2 = index_ic1; index_ic2 < ic_size_[index_md]; index_ic2++) {
              index_ic1_ic2 = index_symmetric_matrix(index_ic1, index_ic2, ic_size_[index_md]);
              for (index_ct = 0; index_ct < ct_size_; index_ct++) {

                if (((int)l > l_max_ct_[index_md][index_ct]) || (is_non_zero_[index_md][index_ic1_ic2] == _FALSE_))
                  cl_md_ic[index_md][index_ic1_ic2*ct_size_ + index_ct] = 0.;
              }
            }
          }
        }
        /* if l was too big, set anyway all components to zero */
        else {
          for (index_ic1 = 0; index_ic1 < ic_size_[index_md]; index_ic1++) {
            for (index_ic2 = index_ic1; index_ic2 < ic_size_[index_md]; index_ic2++) {
              index_ic1_ic2 = index_symmetric_matrix(index_ic1, index_ic2, ic_size_[index_md]);
              for (index_ct = 0; index_ct<ct_size_; index_ct++) {
                cl_md_ic[index_md][index_ic1_ic2*ct_size_ + index_ct] = 0.;
              }
            }
          }
        }

        /* sum up all ic for each mode */

        for (index_ct = 0; index_ct < ct_size_; index_ct++) {

          cl_md[index_md][index_ct]=0.;

          for (index_ic1 = 0; index_ic1 < ic_size_[index_md]; index_ic1++) {
            for (index_ic2 = index_ic1; index_ic2 < ic_size_[index_md]; index_ic2++) {
              index_ic1_ic2 = index_symmetric_matrix(index_ic1, index_ic2, ic_size_[index_md]);

              if (index_ic1 == index_ic2)
                cl_md[index_md][index_ct] += cl_md_ic[index_md][index_ic1_ic2*ct_size_ + index_ct];
              else
                cl_md[index_md][index_ct] += 2.*cl_md_ic[index_md][index_ic1_ic2*ct_size_ + index_ct];
            }
          }
        }
      }

      /** - --> (c.3.) add contribution of cl_md[index_md] to cl_tot */

      for (index_ct = 0; index_ct < ct_size_; index_ct++)
        cl_tot[index_ct]+=cl_md[index_md][index_ct];
    }
  }

  return _SUCCESS_;

}

/**
 * This routine initializes the spectra structure (in particular,
 * computes table of anisotropy and Fourier spectra \f$ C_l^{X}, P(k), ... \f$)
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure (will provide H, Omega_m at redshift of interest)
 * @param ppt Input: pointer to perturbation structure
 * @param ptr Input: pointer to transfer structure
 * @param ppm Input: pointer to primordial structure
 * @param pnl Input: pointer to nonlinear structure
 * @param psp Output: pointer to initialized spectra structure
 * @return the error status
 */

int SpectraModule::spectra_init() {

  /** Summary: */

  /** - check that we really want to compute at least one spectrum */

  if (ppt->has_cls == _FALSE_) {
    md_size_ = 0;
    if (psp->spectra_verbose > 0)
      printf("No spectra requested. Spectra module skipped.\n");
    return _SUCCESS_;
  }
  else {
    if (psp->spectra_verbose > 0)
      printf("Computing unlensed harmonic spectra\n");
  }

  /** - initialize indices and allocate some of the arrays in the
      spectra structure */

  class_call(spectra_indices(),
             error_message_,
             error_message_);

  /** - deal with \f$ C_l\f$'s, if any */

  if (ppt->has_cls == _TRUE_) {

    class_call(spectra_cls(),
               error_message_,
               error_message_);

  }
  else {
    ct_size_ = 0;
  }

  return _SUCCESS_;
}

/**
 * This routine frees all the memory space allocated by spectra_init().
 *
 * To be called at the end of each run, only when no further calls to
 * spectra_cls_at_l(), spectra_pk_at_z(), spectra_pk_at_k_and_z() are needed.
 *
 * @param psp Input: pointer to spectra structure (which fields must be freed)
 * @return the error status
 */

int SpectraModule::spectra_free() {

  if (ppt->has_cls == _FALSE_) {
    return _SUCCESS_;
  }

  int index_md;

  if (md_size_ > 0) {
    if (ct_size_ > 0) {

      for (index_md = 0; index_md < md_size_; index_md++) {
        free(l_max_ct_[index_md]);
        free(cl_[index_md]);
        free(ddcl_[index_md]);
      }
      free(l_);
      free(l_size_);
      free(l_max_ct_);
      free(l_max_);
      free(cl_);
      free(ddcl_);
    }
  }

  for (index_md=0; index_md < md_size_; index_md++)
    free(is_non_zero_[index_md]);

  free(is_non_zero_);
  free(ic_size_);
  free(ic_ic_size_);

  return _SUCCESS_;

}

/**
 * This routine defines indices and allocates tables in the spectra structure
 *
 * @param pba  Input: pointer to background structure
 * @param ppt  Input: pointer to perturbation structure
 * @param ptr  Input: pointer to transfers structure
 * @param ppm  Input: pointer to primordial structure
 * @param psp  Input/output: pointer to spectra structure
 * @return the error status
 */

int SpectraModule::spectra_indices() {

  int index_ct;
  int index_md;
  int index_ic1_ic2;

  md_size_ = perturbations_module_->md_size_;
  if (ppt->has_scalars == _TRUE_)
    index_md_scalars_ = perturbations_module_->index_md_scalars_;

  class_alloc(ic_size_,
              sizeof(int)*md_size_,
              error_message_);

  class_alloc(ic_ic_size_,
              sizeof(int)*md_size_,
              error_message_);

  class_alloc(is_non_zero_,
              sizeof(short *)*md_size_,
              error_message_);

  for (index_md = 0; index_md < md_size_; index_md++) {
    ic_size_[index_md] = primordial_module_->ic_size_[index_md];
    ic_ic_size_[index_md] = primordial_module_->ic_ic_size_[index_md];
    class_alloc(is_non_zero_[index_md],
                sizeof(short)*ic_ic_size_[index_md],
                error_message_);
    for (index_ic1_ic2=0; index_ic1_ic2 < ic_ic_size_[index_md]; index_ic1_ic2++)
      is_non_zero_[index_md][index_ic1_ic2] = primordial_module_->is_non_zero_[index_md][index_ic1_ic2];
  }

  if (ppt->has_cls == _TRUE_) {

    /* types of C_l's relevant for both scalars and tensors: TT, EE, TE */

    index_ct=0;

    if (ppt->has_cl_cmb_temperature == _TRUE_) {
      has_tt_ = _TRUE_;
      index_ct_tt_ = index_ct;
      index_ct++;
    }
    else {
      has_tt_ = _FALSE_;
    }

    if (ppt->has_cl_cmb_polarization == _TRUE_) {
      has_ee_ = _TRUE_;
      index_ct_ee_ = index_ct;
      index_ct++;
    }
    else {
      has_ee_ = _FALSE_;
    }

    if ((ppt->has_cl_cmb_temperature == _TRUE_) &&
        (ppt->has_cl_cmb_polarization == _TRUE_)) {
      has_te_ = _TRUE_;
      index_ct_te_ = index_ct;
      index_ct++;
    }
    else {
      has_te_ = _FALSE_;
    }

    if (ppt->has_cl_cmb_polarization == _TRUE_) {
      has_bb_ = _TRUE_;
      index_ct_bb_ = index_ct;
      index_ct++;
    }
    else {
      has_bb_ = _FALSE_;
    }

    /* types of C_l's relevant only for scalars: phi-phi, T-phi, E-phi, d-d, T-d */

    if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
      has_pp_ = _TRUE_;
      index_ct_pp_ = index_ct;
      index_ct++;
    }
    else {
      has_pp_ = _FALSE_;
    }

    if ((ppt->has_cl_cmb_temperature == _TRUE_) && (ppt->has_cl_cmb_lensing_potential == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
      has_tp_ = _TRUE_;
      index_ct_tp_ = index_ct;
      index_ct++;
    }
    else {
      has_tp_ = _FALSE_;
    }

    ct_size_ = index_ct;

    if ((ppt->has_cl_cmb_polarization == _TRUE_) && (ppt->has_cl_cmb_lensing_potential == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
      has_ep_ = _TRUE_;
      index_ct_ep_ = index_ct;
      index_ct++;
    }
    else {
      has_ep_ = _FALSE_;
    }

    if ((ppt->has_scalars == _TRUE_) &&
        ((ppt->has_cl_number_count == _TRUE_) || (ppt->has_cl_lensing_potential == _TRUE_)))
      d_size_ = ppt->selection_num;
    else
      d_size_ = 0;

    if ((ppt->has_cl_number_count == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
      has_dd_ = _TRUE_;
      index_ct_dd_=index_ct;
      index_ct += (d_size_*(d_size_ + 1) - (d_size_ - psp->non_diag)*(d_size_ - 1 - psp->non_diag))/2;
    }
    else {
      has_dd_ = _FALSE_;
    }

    /* the computation of C_l^Td would require a very good sampling of
       transfer functions over a wide range, and a huge computation
       time. In the current version, we prefer to switch it off, rather
       than either slowing down the code considerably, or producing
       very inaccurate spectra.

       if ((ppt->has_cl_cmb_temperature == _TRUE_) && (ppt->has_cl_number_count == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
       has_td_ = _TRUE_;
       index_ct_td_ = index_ct;
       index_ct += d_size_;
       }
       else {
       has_td_ = _FALSE_;
       }
    */
    has_td_ = _FALSE_;

    if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (ppt->has_cl_number_count == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
      has_pd_ = _TRUE_;
      index_ct_pd_=index_ct;
      index_ct += d_size_;
    }
    else {
      has_pd_ = _FALSE_;
    }

    has_td_ = _FALSE_;

    if ((ppt->has_cl_lensing_potential == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
      has_ll_ = _TRUE_;
      index_ct_ll_ = index_ct;
      index_ct += (d_size_*(d_size_ + 1) - (d_size_ - psp->non_diag)*(d_size_ - 1 - psp->non_diag))/2;
    }
    else {
      has_ll_ = _FALSE_;
    }

    /* the computation of C_l^Tl would require a very good sampling of
       transfer functions over a wide range, and a huge computation
       time. In the current version, we prefer to switch it off, rather
       than either slowing down the code considerably, or producing
       very inaccurate spectra.

       if ((ppt->has_cl_cmb_temperature == _TRUE_) && (ppt->has_cl_lensing_potential == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
       has_tl_ = _TRUE_;
       index_ct_tl_ = index_ct;
       index_ct += d_size_;
       }
       else {
       has_tl_ = _FALSE_;
       }
    */
    has_tl_ = _FALSE_;

    if ((ppt->has_cl_number_count == _TRUE_) && (ppt->has_cl_lensing_potential == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
      has_dl_ = _TRUE_;
      index_ct_dl_ = index_ct;
      index_ct += d_size_*d_size_ - (d_size_ - psp->non_diag)*(d_size_ - 1 - psp->non_diag);
    }
    else {
      has_dl_ = _FALSE_;
    }

    ct_size_ = index_ct;

    /* infer from input quantities the l_max for each mode and type,
       l_max_ct[index_md][index_type].  Maximize it over index_ct, and
       then over index_md. */

    class_alloc(l_max_, sizeof(int*)*md_size_, error_message_);
    class_alloc(l_max_ct_, sizeof(int*)*md_size_, error_message_);
    for (index_md = 0; index_md < md_size_; index_md++) {
      class_calloc(l_max_ct_[index_md], ct_size_, sizeof(int), error_message_);
    }

    if (ppt->has_scalars == _TRUE_) {

      /* spectra computed up to l_scalar_max */

      if (has_tt_ == _TRUE_) l_max_ct_[perturbations_module_->index_md_scalars_][index_ct_tt_] = ppt->l_scalar_max;
      if (has_ee_ == _TRUE_) l_max_ct_[perturbations_module_->index_md_scalars_][index_ct_ee_] = ppt->l_scalar_max;
      if (has_te_ == _TRUE_) l_max_ct_[perturbations_module_->index_md_scalars_][index_ct_te_] = ppt->l_scalar_max;
      if (has_pp_ == _TRUE_) l_max_ct_[perturbations_module_->index_md_scalars_][index_ct_pp_] = ppt->l_scalar_max;
      if (has_tp_ == _TRUE_) l_max_ct_[perturbations_module_->index_md_scalars_][index_ct_tp_] = ppt->l_scalar_max;
      if (has_ep_ == _TRUE_) l_max_ct_[perturbations_module_->index_md_scalars_][index_ct_ep_] = ppt->l_scalar_max;

      /* spectra computed up to l_lss_max */

      if (has_dd_ == _TRUE_)
        for (index_ct = index_ct_dd_;
             index_ct < index_ct_dd_ + (d_size_*(d_size_ + 1) - (d_size_ - psp->non_diag)*(d_size_ - 1 - psp->non_diag))/2;
             index_ct++)
          l_max_ct_[perturbations_module_->index_md_scalars_][index_ct] = ppt->l_lss_max;

      if (has_td_ == _TRUE_)
        for (index_ct = index_ct_td_;
             index_ct < index_ct_td_ + d_size_;
             index_ct++)
          l_max_ct_[perturbations_module_->index_md_scalars_][index_ct] = MIN(ppt->l_scalar_max, ppt->l_lss_max);

      if (has_pd_ == _TRUE_)
        for (index_ct = index_ct_pd_;
             index_ct < index_ct_pd_ + d_size_;
             index_ct++)
          l_max_ct_[perturbations_module_->index_md_scalars_][index_ct] = MIN(ppt->l_scalar_max, ppt->l_lss_max);

      if (has_ll_ == _TRUE_)
        for (index_ct = index_ct_ll_;
             index_ct < index_ct_ll_ + (d_size_*(d_size_ + 1) - (d_size_ - psp->non_diag)*(d_size_ - 1-psp->non_diag))/2;
             index_ct++)
          l_max_ct_[perturbations_module_->index_md_scalars_][index_ct] = ppt->l_lss_max;

      if (has_tl_ == _TRUE_)
        for (index_ct = index_ct_tl_;
             index_ct < index_ct_tl_ + d_size_;
             index_ct++)
          l_max_ct_[perturbations_module_->index_md_scalars_][index_ct] = MIN(ppt->l_scalar_max, ppt->l_lss_max);

      if (has_dl_ == _TRUE_)
        for (index_ct = index_ct_dl_;
             index_ct < index_ct_dl_ + (d_size_*d_size_ - (d_size_ - psp->non_diag)*(d_size_ - 1 - psp->non_diag));
             index_ct++)
          l_max_ct_[perturbations_module_->index_md_scalars_][index_ct] = ppt->l_lss_max;

    }
    if (ppt->has_tensors == _TRUE_) {

      /* spectra computed up to l_tensor_max */

      if (has_tt_ == _TRUE_) l_max_ct_[perturbations_module_->index_md_tensors_][index_ct_tt_] = ppt->l_tensor_max;
      if (has_ee_ == _TRUE_) l_max_ct_[perturbations_module_->index_md_tensors_][index_ct_ee_] = ppt->l_tensor_max;
      if (has_te_ == _TRUE_) l_max_ct_[perturbations_module_->index_md_tensors_][index_ct_te_] = ppt->l_tensor_max;
      if (has_bb_ == _TRUE_) l_max_ct_[perturbations_module_->index_md_tensors_][index_ct_bb_] = ppt->l_tensor_max;
    }

    /* maximizations */
    l_max_tot_ = 0.;
    for (index_md = 0; index_md < md_size_; index_md++) {
      l_max_[index_md] = 0.;
      for (index_ct = 0.; index_ct < ct_size_; index_ct++)
        l_max_[index_md] = MAX(l_max_[index_md], l_max_ct_[index_md][index_ct]);
      l_max_tot_ = MAX(l_max_tot_, l_max_[index_md]);
    }
  }

  return _SUCCESS_;

}

/**
 * This routine computes a table of values for all harmonic spectra \f$ C_l \f$'s,
 * given the transfer functions and primordial spectra.
 *
 * @param pba Input: pointer to background structure
 * @param ppt Input: pointer to perturbation structure
 * @param ptr Input: pointer to transfers structure
 * @param ppm Input: pointer to primordial structure
 * @param psp Input/Output: pointer to spectra structure
 * @return the error status
 */

int SpectraModule::spectra_cls() {

  /** Summary: */

  /** - define local variables */

  int index_md;
  int index_ic1,index_ic2,index_ic1_ic2;
  int index_l;
  int index_ct;
  int cl_integrand_num_columns;

  /** - allocate pointers to arrays where results will be stored */

  class_alloc(l_size_, sizeof(int)*md_size_, error_message_);
  class_alloc(cl_, sizeof(double *)*md_size_, error_message_);
  class_alloc(ddcl_, sizeof(double *)*md_size_, error_message_);

  l_size_max_ = transfer_module_->l_size_max_;
  class_alloc(l_, sizeof(double)*l_size_max_,error_message_);

  /** - store values of l */
  for (index_l = 0; index_l < l_size_max_; index_l++) {
    l_[index_l] = (double)transfer_module_->l_[index_l];
  }

  Tools::TaskSystem task_system(pba->number_of_threads);
  std::vector<std::future<int>> future_output;

  /** - loop over modes (scalar, tensors, etc). For each mode: */

  for (index_md = 0; index_md < md_size_; index_md++) {

    /** - --> (a) store number of l values for this mode */

    l_size_[index_md] = transfer_module_->l_size_[index_md];

    /** - --> (b) allocate arrays where results will be stored */

    class_alloc(cl_[index_md],sizeof(double)*l_size_[index_md]*ct_size_*ic_ic_size_[index_md],error_message_);
    class_alloc(ddcl_[index_md], sizeof(double)*l_size_[index_md]*ct_size_*ic_ic_size_[index_md], error_message_);
    cl_integrand_num_columns = 1 + ct_size_*2; /* one for k, ct_size_ for each type, ct_size_ for each second derivative of each type */

    /** - --> (c) loop over initial conditions */

    for (index_ic1 = 0; index_ic1 < ic_size_[index_md]; index_ic1++) {
      for (index_ic2 = index_ic1; index_ic2 < ic_size_[index_md]; index_ic2++) {
        index_ic1_ic2 = index_symmetric_matrix(index_ic1, index_ic2, ic_size_[index_md]);

        /* non-diagonal coefficients should be computed only if non-zero correlation */
        if (is_non_zero_[index_md][index_ic1_ic2] == _TRUE_) {

          future_output.push_back(task_system.AsyncTask([this, index_md, cl_integrand_num_columns, index_ic1, index_ic2] () {
            double * cl_integrand; /* array with argument cl_integrand[index_k*cl_integrand_num_columns+1+psp->index_ct] */
            double * transfer_ic1; /* array with argument transfer_ic1[index_tt] */
            double * transfer_ic2; /* idem */
            double * primordial_pk;  /* array with argument primordial_pk[index_ic_ic]*/


            class_alloc(cl_integrand, transfer_module_->q_size_*cl_integrand_num_columns*sizeof(double), error_message_);
            class_alloc(primordial_pk, ic_ic_size_[index_md]*sizeof(double), error_message_);
            class_alloc(transfer_ic1, transfer_module_->tt_size_[index_md]*sizeof(double), error_message_);
            class_alloc(transfer_ic2, transfer_module_->tt_size_[index_md]*sizeof(double), error_message_);

            /** - ---> loop over l values defined in the transfer module.
                For each l, compute the \f$ C_l\f$'s for all types (TT, TE, ...)
                by convolving primordial spectra with transfer  functions.
                This elementary task is assigned to spectra_compute_cl() */

            for (int index_l = 0; index_l < transfer_module_->l_size_[index_md]; index_l++) {

#pragma omp flush(abort)

              class_call(spectra_compute_cl(
                                            index_md,
                                            index_ic1,
                                            index_ic2,
                                            index_l,
                                            cl_integrand_num_columns,
                                            cl_integrand,
                                            primordial_pk,
                                            transfer_ic1,
                                            transfer_ic2),
                         error_message_,
                         error_message_);

            } /* end of loop over l */

            free(cl_integrand);

            free(primordial_pk);

            free(transfer_ic1);

            free(transfer_ic2);

            return _SUCCESS_;

          }));
        }
        else {

          /* set non-diagonal coefficients to zero if pair of ic's uncorrelated */

          for (index_l = 0; index_l < transfer_module_->l_size_[index_md]; index_l++) {
            for (index_ct = 0; index_ct < ct_size_; index_ct++) {
              cl_[index_md][(index_l*ic_ic_size_[index_md] + index_ic1_ic2)*ct_size_ + index_ct] = 0.;
            }
          }
        }
      }
    }

    for (std::future<int>& future : future_output) {
        if (future.get() != _SUCCESS_) return _FAILURE_;
    }
    future_output.clear();

    /** - --> (d) now that for a given mode, all possible \f$ C_l\f$'s have been computed,
        compute second derivative of the array in which they are stored,
        in view of spline interpolation. */

    class_call(array_spline_table_lines(l_,
                                        l_size_[index_md],
                                        cl_[index_md],
                                        ic_ic_size_[index_md]*ct_size_,
                                        ddcl_[index_md],
                                        _SPLINE_EST_DERIV_,
                                        error_message_),
               error_message_,
               error_message_);
  }

  return _SUCCESS_;

}

/**
 * This routine computes the \f$ C_l\f$'s for a given mode, pair of initial conditions
 * and multipole, but for all types (TT, TE...), by convolving the
 * transfer functions with the primordial spectra.
 *
 * @param pba           Input: pointer to background structure
 * @param ppt           Input: pointer to perturbation structure
 * @param ptr           Input: pointer to transfers structure
 * @param ppm           Input: pointer to primordial structure
 * @param psp           Input/Output: pointer to spectra structure (result stored here)
 * @param index_md      Input: index of mode under consideration
 * @param index_ic1     Input: index of first initial condition in the correlator
 * @param index_ic2     Input: index of second initial condition in the correlator
 * @param index_l       Input: index of multipole under consideration
 * @param cl_integrand_num_columns Input: number of columns in cl_integrand
 * @param cl_integrand  Input: an allocated workspace
 * @param primordial_pk Input: table of primordial spectrum values
 * @param transfer_ic1  Input: table of transfer function values for first initial condition
 * @param transfer_ic2  Input: table of transfer function values for second initial condition
 * @return the error status
 */

int SpectraModule::spectra_compute_cl(int index_md,
                       int index_ic1,
                       int index_ic2,
                       int index_l,
                       int cl_integrand_num_columns,
                       double * cl_integrand,
                       double * primordial_pk,
                       double * transfer_ic1,
                       double * transfer_ic2
                       ) {

  int index_q;
  int index_tt;
  int index_ct;
  int index_d1,index_d2;
  double k;
  double clvalue;
  int index_ic1_ic2;
  double transfer_ic1_temp=0.;
  double transfer_ic2_temp=0.;
  double * transfer_ic1_nc=NULL;
  double * transfer_ic2_nc=NULL;
  double factor;
  int index_q_spline=0;

  index_ic1_ic2 = index_symmetric_matrix(index_ic1, index_ic2, ic_size_[index_md]);

  if (ppt->has_cl_number_count == _TRUE_) {
    class_alloc(transfer_ic1_nc, d_size_*sizeof(double), error_message_);
    class_alloc(transfer_ic2_nc, d_size_*sizeof(double), error_message_);
  }

  for (index_q = 0; index_q < transfer_module_->q_size_; index_q++) {

    k = transfer_module_->k_[index_md][index_q];

    cl_integrand[index_q*cl_integrand_num_columns+0] = k;

    class_call(primordial_module_->primordial_spectrum_at_k(index_md, linear, k, primordial_pk),
               primordial_module_->error_message_,
               error_message_);

    /* above routine checks that k>0: no possible division by zero below */

    for (index_tt = 0; index_tt < transfer_module_->tt_size_[index_md]; index_tt++) {

      transfer_ic1[index_tt] =
        transfer_module_->transfer_[index_md]
        [((index_ic1*transfer_module_->tt_size_[index_md] + index_tt)
          *transfer_module_->l_size_[index_md] + index_l)
         *transfer_module_->q_size_ + index_q];

      if (index_ic1 == index_ic2) {
        transfer_ic2[index_tt] = transfer_ic1[index_tt];
      }
      else {
        transfer_ic2[index_tt] = transfer_module_->transfer_[index_md]
          [((index_ic2*transfer_module_->tt_size_[index_md] + index_tt)
            *transfer_module_->l_size_[index_md] + index_l)
           *transfer_module_->q_size_ + index_q];
      }
    }

    /* define combinations of transfer functions */

    if (ppt->has_cl_cmb_temperature == _TRUE_) {

      if (_scalarsEXT_) {

        transfer_ic1_temp = transfer_ic1[transfer_module_->index_tt_t0_] + transfer_ic1[transfer_module_->index_tt_t1_] + transfer_ic1[transfer_module_->index_tt_t2_];
        transfer_ic2_temp = transfer_ic2[transfer_module_->index_tt_t0_] + transfer_ic2[transfer_module_->index_tt_t1_] + transfer_ic2[transfer_module_->index_tt_t2_];

      }

      if (_vectorsEXT_) {

        transfer_ic1_temp = transfer_ic1[transfer_module_->index_tt_t1_] + transfer_ic1[transfer_module_->index_tt_t2_];
        transfer_ic2_temp = transfer_ic2[transfer_module_->index_tt_t1_] + transfer_ic2[transfer_module_->index_tt_t2_];

      }

      if (_tensorsEXT_) {

        transfer_ic1_temp = transfer_ic1[transfer_module_->index_tt_t2_];
        transfer_ic2_temp = transfer_ic2[transfer_module_->index_tt_t2_];

      }
    }

    if ((_scalarsEXT_) && (ppt->has_cl_number_count == _TRUE_)) {

      for (index_d1 = 0; index_d1 < d_size_; index_d1++) {

        transfer_ic1_nc[index_d1] = 0.;
        transfer_ic2_nc[index_d1] = 0.;

        if (ppt->has_nc_density == _TRUE_) {
          transfer_ic1_nc[index_d1] += transfer_ic1[transfer_module_->index_tt_density_ + index_d1];
          transfer_ic2_nc[index_d1] += transfer_ic2[transfer_module_->index_tt_density_ + index_d1];
        }

        if (ppt->has_nc_rsd     == _TRUE_) {
          transfer_ic1_nc[index_d1]
            += transfer_ic1[transfer_module_->index_tt_rsd_ + index_d1]
            + transfer_ic1[transfer_module_->index_tt_d0_ + index_d1]
            + transfer_ic1[transfer_module_->index_tt_d1_ + index_d1];
          transfer_ic2_nc[index_d1]
            += transfer_ic2[transfer_module_->index_tt_rsd_ + index_d1]
            + transfer_ic2[transfer_module_->index_tt_d0_ + index_d1]
            + transfer_ic2[transfer_module_->index_tt_d1_ + index_d1];
        }

        if (ppt->has_nc_lens == _TRUE_) {
          transfer_ic1_nc[index_d1] += l_[index_l]*(l_[index_l] + 1.)*transfer_ic1[transfer_module_->index_tt_nc_lens_ + index_d1];
          transfer_ic2_nc[index_d1] += l_[index_l]*(l_[index_l] + 1.)*transfer_ic2[transfer_module_->index_tt_nc_lens_ + index_d1];
        }

        if (ppt->has_nc_gr == _TRUE_) {
          transfer_ic1_nc[index_d1]
            += transfer_ic1[transfer_module_->index_tt_nc_g1_ + index_d1]
            + transfer_ic1[transfer_module_->index_tt_nc_g2_ + index_d1]
            + transfer_ic1[transfer_module_->index_tt_nc_g3_ + index_d1]
            + transfer_ic1[transfer_module_->index_tt_nc_g4_ + index_d1]
            + transfer_ic1[transfer_module_->index_tt_nc_g5_ + index_d1];
          transfer_ic2_nc[index_d1]
            += transfer_ic2[transfer_module_->index_tt_nc_g1_ + index_d1]
            + transfer_ic2[transfer_module_->index_tt_nc_g2_ + index_d1]
            + transfer_ic2[transfer_module_->index_tt_nc_g3_ + index_d1]
            + transfer_ic2[transfer_module_->index_tt_nc_g4_ + index_d1]
            + transfer_ic2[transfer_module_->index_tt_nc_g5_ + index_d1];
        }

      }
    }

    /* integrand of Cl's */

    /* note: we must integrate

       C_l = int [4 pi dk/k calP(k) Delta1_l(q) Delta2_l(q)]

       where calP(k) is the dimensionless
       power spectrum equal to a constant in the scale-invariant case,
       and to P(k) = A_s k^(ns-1) otherwise and q=sqrt(k2+K) (scalars)
       or sqrt(k2+2K) (vectors) or sqrt(k2+3K) (tensors)

       In the literature, people often rewrite the integral in terms
       of q and absorb the Jacobian of the change of variables in a redefinition of the primodial
       spectrum. Let us illustrate this for scalars:

       dk/k = kdk/k2 = qdq/k2 = dq/q * (q/k)^2 = dq/q * [q2/(q2-K)] = q2dq * 1/[q(q2-K)]

       This factor 1/[q(q2-K)] is commonly absorbed in the definition of calP. Then one would have

       C_l = int [4 pi q2 dq {A_s k^(ns-1)/[q(q2-K)]} Delta1_l(q) Delta2_l(q)]

       Sometimes in the literature, the factor (k2-3K)=(q2-4K) present
       in the initial conditions of scalar transfer functions (if
       normalized to curvature R=1) is also absorbed in the definition
       of the power spectrum. Then the curvature power spectrum reads

       calP = (q2-4K)/[q(q2-K)] * (k/k)^ns

       In CLASS we prefer to define calP = (k/k)^ns like in the flat
       case, to have the factor (q2-4K) in the initialk conditions,
       and the factor 1/[q(q2-K)] doesn't need to be there since we
       integrate over dk/k.

       For tensors, the change of variable described above gives a slightly different result:

       dk/k = kdk/k2 = qdq/k2 = dq/q * (q/k)^2 = dq/q * [q2/(q2-3K)] = q2dq * 1/[q(q2-3K)]

       But for tensors there are extra curvature-related correction factors to
       take into account. See the comments in the perturbation module,
       related to initial conditions for tensors.

    */

    factor = 4. * _PI_ / k;

    if (has_tt_ == _TRUE_)
      cl_integrand[index_q*cl_integrand_num_columns + 1 + index_ct_tt_] =
        primordial_pk[index_ic1_ic2]
        * transfer_ic1_temp
        * transfer_ic2_temp
        * factor;

    if (has_ee_ == _TRUE_)
      cl_integrand[index_q*cl_integrand_num_columns + 1 + index_ct_ee_] =
        primordial_pk[index_ic1_ic2]
        * transfer_ic1[transfer_module_->index_tt_e_]
        * transfer_ic2[transfer_module_->index_tt_e_]
        * factor;

    if (has_te_ == _TRUE_)
      cl_integrand[index_q*cl_integrand_num_columns + 1 + index_ct_te_] =
        primordial_pk[index_ic1_ic2]
        *0.5*(transfer_ic1_temp*transfer_ic2[transfer_module_->index_tt_e_] +
              transfer_ic1[transfer_module_->index_tt_e_]*transfer_ic2_temp)
        * factor;

    if (_tensorsEXT_ && (has_bb_ == _TRUE_))
      cl_integrand[index_q*cl_integrand_num_columns + 1 + index_ct_bb_] =
        primordial_pk[index_ic1_ic2]
        * transfer_ic1[transfer_module_->index_tt_b_]
        * transfer_ic2[transfer_module_->index_tt_b_]
        * factor;

    if (_scalarsEXT_ && (has_pp_ == _TRUE_))
      cl_integrand[index_q*cl_integrand_num_columns + 1 + index_ct_pp_] =
        primordial_pk[index_ic1_ic2]
        * transfer_ic1[transfer_module_->index_tt_lcmb_]
        * transfer_ic2[transfer_module_->index_tt_lcmb_]
        * factor;

    if (_scalarsEXT_ && (has_tp_ == _TRUE_))
      cl_integrand[index_q*cl_integrand_num_columns + 1 + index_ct_tp_] =
        primordial_pk[index_ic1_ic2]
        *0.5*(transfer_ic1_temp*transfer_ic2[transfer_module_->index_tt_lcmb_] +
              transfer_ic1[transfer_module_->index_tt_lcmb_]*transfer_ic2_temp)
        * factor;

    if (_scalarsEXT_ && (has_ep_ == _TRUE_))
      cl_integrand[index_q*cl_integrand_num_columns + 1 + index_ct_ep_] =
        primordial_pk[index_ic1_ic2]
        *0.5*(transfer_ic1[transfer_module_->index_tt_e_]*transfer_ic2[transfer_module_->index_tt_lcmb_] +
              transfer_ic1[transfer_module_->index_tt_lcmb_]*transfer_ic2[transfer_module_->index_tt_e_])
        * factor;

    if (_scalarsEXT_ && (has_dd_ == _TRUE_)) {
      index_ct=0;
      for (index_d1 = 0; index_d1 < d_size_; index_d1++) {
        for (index_d2 = index_d1; index_d2 <= MIN(index_d1 + psp->non_diag, d_size_ - 1); index_d2++) {
          cl_integrand[index_q*cl_integrand_num_columns + 1 + index_ct_dd_ + index_ct] =
            primordial_pk[index_ic1_ic2]
            * transfer_ic1_nc[index_d1]
            * transfer_ic2_nc[index_d2]
            * factor;
          index_ct++;
        }
      }
    }

    if (_scalarsEXT_ && (has_td_ == _TRUE_)) {
      for (index_d1 = 0; index_d1 < d_size_; index_d1++) {
        cl_integrand[index_q*cl_integrand_num_columns + 1 + index_ct_td_ + index_d1] =
          primordial_pk[index_ic1_ic2]
          * 0.5*(transfer_ic1_temp * transfer_ic2_nc[index_d1] +
                 transfer_ic1_nc[index_d1] * transfer_ic2_temp)
          * factor;
      }
    }

    if (_scalarsEXT_ && (has_pd_ == _TRUE_)) {
      for (index_d1 = 0; index_d1 < d_size_; index_d1++) {
        cl_integrand[index_q*cl_integrand_num_columns + 1 + index_ct_pd_ + index_d1]=
          primordial_pk[index_ic1_ic2]
          *0.5*(transfer_ic1[transfer_module_->index_tt_lcmb_]*transfer_ic2_nc[index_d1] +
                transfer_ic1_nc[index_d1]*transfer_ic2[transfer_module_->index_tt_lcmb_])
          * factor;
      }
    }

    if (_scalarsEXT_ && (has_ll_ == _TRUE_)) {
      index_ct=0;
      for (index_d1 = 0; index_d1 < d_size_; index_d1++) {
        for (index_d2 = index_d1; index_d2 <= MIN(index_d1 + psp->non_diag, d_size_ - 1); index_d2++) {
          cl_integrand[index_q*cl_integrand_num_columns + 1 + index_ct_ll_ + index_ct] =
            primordial_pk[index_ic1_ic2]
            * transfer_ic1[transfer_module_->index_tt_lensing_ + index_d1]
            * transfer_ic2[transfer_module_->index_tt_lensing_ + index_d2]
            * factor;
          index_ct++;
        }
      }
    }

    if (_scalarsEXT_ && (has_tl_ == _TRUE_)) {
      for (index_d1 = 0; index_d1 < d_size_; index_d1++) {
        cl_integrand[index_q*cl_integrand_num_columns + 1 + index_ct_tl_ + index_d1] =
          primordial_pk[index_ic1_ic2]
          *0.5*(transfer_ic1_temp*transfer_ic2[transfer_module_->index_tt_lensing_ + index_d1] +
                transfer_ic1[transfer_module_->index_tt_lensing_ + index_d1]*transfer_ic2_temp)
          * factor;
      }
    }

    if (_scalarsEXT_ && (has_dl_ == _TRUE_)) {
      index_ct=0;
      for (index_d1=0; index_d1<d_size_; index_d1++) {
        for (index_d2 = MAX(index_d1 - psp->non_diag, 0); index_d2 <= MIN(index_d1 + psp->non_diag, d_size_ - 1); index_d2++) {
          cl_integrand[index_q*cl_integrand_num_columns + 1 + index_ct_dl_ + index_ct] =
            primordial_pk[index_ic1_ic2]
            *transfer_ic1_nc[index_d1]*transfer_ic2[transfer_module_->index_tt_lensing_ + index_d2]
            * factor;
          index_ct++;
        }
      }
    }
  }

  for (index_ct = 0; index_ct < ct_size_; index_ct++) {

    /* treat null spectra (C_l^BB of scalars, C_l^pp of tensors, etc. */

    if ((_scalarsEXT_ && (has_bb_ == _TRUE_) && (index_ct == index_ct_bb_)) ||
        (_tensorsEXT_ && (has_pp_ == _TRUE_) && (index_ct == index_ct_pp_)) ||
        (_tensorsEXT_ && (has_tp_ == _TRUE_) && (index_ct == index_ct_tp_)) ||
        (_tensorsEXT_ && (has_ep_ == _TRUE_) && (index_ct == index_ct_ep_)) ||
        (_tensorsEXT_ && (has_dd_ == _TRUE_) && (index_ct == index_ct_dd_)) ||
        (_tensorsEXT_ && (has_td_ == _TRUE_) && (index_ct == index_ct_td_)) ||
        (_tensorsEXT_ && (has_pd_ == _TRUE_) && (index_ct == index_ct_pd_)) ||
        (_tensorsEXT_ && (has_ll_ == _TRUE_) && (index_ct == index_ct_ll_)) ||
        (_tensorsEXT_ && (has_tl_ == _TRUE_) && (index_ct == index_ct_tl_)) ||
        (_tensorsEXT_ && (has_dl_ == _TRUE_) && (index_ct == index_ct_dl_))
        ) {

      cl_[index_md][(index_l*ic_ic_size_[index_md] + index_ic1_ic2)*ct_size_ + index_ct] = 0.;

    }
    /* for non-zero spectra, integrate over q */
    else {

      /* spline the integrand over the whole range of k's */

      class_call(array_spline(cl_integrand,
                              cl_integrand_num_columns,
                              transfer_module_->q_size_,
                              0,
                              1+index_ct,
                              1 + ct_size_ + index_ct,
                              _SPLINE_EST_DERIV_,
                              error_message_),
                 error_message_,
                 error_message_);

      /* Technical point: we will now do a spline integral over the
         whole range of k's, excepted in the closed (K>0) case. In
         that case, it is a bad idea to spline over the values of k
         corresponding to nu<nu_flat_approximation. In this region, nu
         values are integer values, so the steps dq and dk have some
         discrete jumps. This makes the spline routine less accurate
         than a trapezoidal integral with finer sampling. So, in the
         closed case, we set index_q_spline to
         transfer_module_->index_q_flat_approximation_, to tell the integration
         routine that below this index, it should treat the integral
         as a trapezoidal one. For testing, one is free to set
         index_q_spline to 0, to enforce spline integration
         everywhere, or to (transfer_module_->q_size_-1), to enforce trapezoidal
         integration everywhere. */

      if (pba->sgnK == 1) {
        index_q_spline = transfer_module_->index_q_flat_approximation_;
      }

      class_call(array_integrate_all_trapzd_or_spline(cl_integrand,
                                                      cl_integrand_num_columns,
                                                      transfer_module_->q_size_,
                                                      index_q_spline,
                                                      0,
                                                      1+index_ct,
                                                      1 + ct_size_ + index_ct,
                                                      &clvalue,
                                                      error_message_),
                 error_message_,
                 error_message_);

      /* in the closed case, instead of an integral, we have a
         discrete sum. In practice, this does not matter: the previous
         routine does give a correct approximation of the discrete
         sum, both in the trapezoidal and spline regions. The only
         error comes from the first point: the previous routine
         assumes a weight for the first point which is too small
         compared to what it would be in the an actual discrete
         sum. The line below correct this problem in an exact way.
      */

      if (pba->sgnK == 1) {
        clvalue += cl_integrand[1 + index_ct]*transfer_module_->q_[0]/transfer_module_->k_[0][0]*sqrt(pba->K)/2.;
      }

      /* we have the correct C_l now. We can store it in the transfer structure. */

      cl_[index_md][(index_l*ic_ic_size_[index_md] + index_ic1_ic2)*ct_size_ + index_ct] = clvalue;

    }
  }

  if (ppt->has_cl_number_count == _TRUE_) {
    free(transfer_ic1_nc);
    free(transfer_ic2_nc);
  }

  return _SUCCESS_;

}

  /* deprecated functions (since v2.8) */

/**
 * Matter power spectrum for arbitrary redshift and for all initial conditions.
 *
 * This function is deprecated since v2.8. Try using nonlinear_pk_at_z() instead.
 *
 * @param pba           Input: pointer to background structure (used for converting z into tau)
 * @param psp           Input: pointer to spectra structure (containing pre-computed table)
 * @param mode          Input: linear or logarithmic
 * @param z             Input: redshift
 * @param output_tot    Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$ (linear mode), or its logarithms (logarithmic mode)
 * @param output_ic     Output: for each pair of initial conditions, matter power spectra P(k) in \f$ Mpc^3 \f$ (linear mode), or their logarithms and cross-correlation angles (logarithmic mode)
 * @param output_cb_tot Output: CDM+baryon power spectrum P_cb(k) in \f$ Mpc^3 \f$ (linear mode), or its logarithms (logarithmic mode)
 * @param output_cb_ic  Output: for each pair of initial conditions, CDM+baryon power spectra P_cb(k) in \f$ Mpc^3 \f$ (linear mode), or their logarithms and cross-correlation angles (logarithmic mode)
 * @return the error status
 */

int SpectraModule::spectra_pk_at_z(enum linear_or_logarithmic mode,
                    double z,
                    double * output_tot,    /* array with argument output_tot[index_k] (must be already allocated) */
                    double * output_ic,     /* array with argument output_tot[index_k * ic_ic_size_[index_md] + index_ic1_ic2] (must be already allocated only if more than one initial condition) */
                    double * output_cb_tot, /* same as output_tot for the baryon+CDM only */
                    double * output_cb_ic   /* same as output_ic  for the baryon+CDM only */ ) {

  fprintf(stderr," -> [WARNING:] You are calling the function spectra_pk_at_z() which is deprecated since v2.8. Try using nonlinear_pk_at_z() instead.\n");

  class_call(nonlinear_module_->nonlinear_pks_at_z(
                                                  mode,
                                                  pk_linear,
                                                  z,
                                                  output_tot,
                                                  output_ic,
                                                  output_cb_tot,
                                                  output_cb_ic
                                                  ),
             nonlinear_module_->error_message_,
             error_message_);

  return _SUCCESS_;

}

/**
 * Matter power spectrum for arbitrary wavenumber, redshift and initial condition.
 *
 * This function is deprecated since v2.8. Try using nonlinear_pk_linear_at_k_and_z() instead.
 *
 * @param pba        Input: pointer to background structure (used for converting z into tau)
 * @param ppm        Input: pointer to primordial structure (used only in the case 0 < k < kmin)
 * @param psp        Input: pointer to spectra structure (containing pre-computed table)
 * @param k          Input: wavenumber in 1/Mpc
 * @param z          Input: redshift
 * @param pk_tot     Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
 * @param pk_ic      Output: for each pair of initial conditions, matter power spectra P(k) in \f$ Mpc^3\f$
 * @param pk_cb_tot  Output: b+CDM power spectrum P(k) in \f$ Mpc^3 \f$
 * @param pk_cb_ic   Output: for each pair of initial conditions, b+CDM power spectra P(k) in \f$ Mpc^3\f$
 * @return the error status
 */

int SpectraModule::spectra_pk_at_k_and_z(
                          double k,
                          double z,
                          double * pk_tot,    /* pointer to a single number (must be already allocated) */
                          double * pk_ic,     /* array of argument pk_ic[index_ic1_ic2]
                                                 (must be already allocated only if several initial conditions) */
                          double * pk_cb_tot, /* same as pk_tot for baryon+CDM part only */
                          double * pk_cb_ic   /* same as pk_ic  for baryon+CDM part only */) {

  fprintf(stderr," -> [WARNING:] You are calling the function spectra_pk_at_k_and_z() which is deprecated since v2.8. Try using nonlinear_pk_linear_at_k_and_z() instead.\n");

  class_call(nonlinear_module_->nonlinear_pks_at_k_and_z(pk_linear,
                                                        k,
                                                        z,
                                                        pk_tot,
                                                        pk_ic,
                                                        pk_cb_tot,
                                                        pk_cb_ic),
             nonlinear_module_->error_message_,
             error_message_);

  return _SUCCESS_;
}

/**
 * Non-linear total matter power spectrum for arbitrary redshift.
 *
 * This function is deprecated since v2.8. Try using nonlinear_pk_at_z() instead.
 *
 * @param pba           Input: pointer to background structure (used for converting z into tau)
 * @param psp           Input: pointer to spectra structure (containing pre-computed table)
 * @param mode          Input: linear or logarithmic
 * @param z             Input: redshift
 * @param output_tot    Output: total matter power spectrum P(k) in \f$ Mpc^3\f$ (linear mode), or its logarithms (logarithmic mode)
 * @param output_cb_tot Output: b+CDM power spectrum P(k) in \f$ Mpc^3\f$ (linear mode), or its logarithms (logarithmic mode)
 * @return the error status
 */

int SpectraModule::spectra_pk_nl_at_z(
                                      enum linear_or_logarithmic mode,
                                      double z,
                                      double * output_tot,   /* array with argument output_tot[index_k] (must be already allocated) */
                                      double * output_cb_tot
                                      ) {

  fprintf(stderr," -> [WARNING:] You are calling the function spectra_pk_nl_at_z() which is deprecated since v2.8. Try using nonlinear_pk_at_z() instead.\n");

  class_call(nonlinear_module_->nonlinear_pks_at_z(mode,
                                                  pk_nonlinear,
                                                  z,
                                                  output_tot,
                                                  NULL,
                                                  output_cb_tot,
                                                  NULL
                                                  ),
             nonlinear_module_->error_message_,
             error_message_);

  return _SUCCESS_;

}

/**
 * Non-linear total matter power spectrum for arbitrary wavenumber and redshift.
 *
 * This function is deprecated since v2.8. Try using nonlinear_pk_at_k_and_z() instead.
 *
 * @param pba        Input: pointer to background structure (used for converting z into tau)
 * @param ppm        Input: pointer to primordial structure (used only in the case 0 < k < kmin)
 * @param psp        Input: pointer to spectra structure (containing pre-computed table)
 * @param k          Input: wavenumber in 1/Mpc
 * @param z          Input: redshift
 * @param pk_tot     Output: total matter power spectrum P(k) in \f$ Mpc^3\f$
 * @param pk_cb_tot  Output: b+CDM power spectrum P(k) in \f$ Mpc^3\f$
 * @return the error status
 */

int SpectraModule::spectra_pk_nl_at_k_and_z(
                             double k,
                             double z,
                             double * pk_tot,   /* pointer to a single number (must be already allocated) */
                             double * pk_cb_tot /* same as pk_tot for baryon+CDM only */
                             ) {

  fprintf(stderr," -> [WARNING:] You are calling the function spectra_pk_nl_at_k_and_z() which is deprecated since v2.8. Try using nonlinear_pk_at_k_and_z() instead.\n");

  class_call(nonlinear_module_->nonlinear_pks_at_k_and_z(pk_nonlinear, k, z, pk_tot, NULL, pk_cb_tot, NULL),
             nonlinear_module_->error_message_,
             error_message_);

  return _SUCCESS_;

}

/**
 * Return the P(k,z) for a grid of (k_i,z_j) passed in input,
 * for all available pk types (_m, _cb),
 * either linear or nonlinear depending on input.
 *
 * This function is deprecated since v2.8. Try using nonlinear_pks_at_kvec_and_zvec() instead.
 *
 * @param pba            Input: pointer to background structure
 * @param psp            Input: pointer to spectra structure
 * @param kvec           Input: array of wavenumbers in ascending order (in 1/Mpc)
 * @param kvec_size      Input: size of array of wavenumbers
 * @param zvec           Input: array of redshifts in arbitrary order
 * @param zvec_size      Input: size of array of redshifts
 * @param pk_tot_out     Output: P(k_i,z_j) for total matter (if available) in Mpc**3
 * @param pk_cb_tot_out  Output: P_cb(k_i,z_j) for cdm+baryons (if available) in Mpc**3
 * @param nonlinear      Input: _TRUE_ or _FALSE_ (to output nonlinear or linear P(k,z))
 * @return the error status
 */

int SpectraModule::spectra_fast_pk_at_kvec_and_zvec(
                                     double * kvec,
                                     int kvec_size,
                                     double * zvec,
                                     int zvec_size,
                                     double * pk_tot_out, // pk_tot_out[index_zvec*kvec_size+index_kvec],
                                                          // already allocated
                                                          //(or NULL if user knows there is no _m output)
                                     double * pk_cb_tot_out, // idem
                                     int nonlinear
                                     ) {
  enum pk_outputs pk_output;

  fprintf(stderr," -> [WARNING:] You are calling the function spectra_fast_pks_at_kvec_and_zvec() which is deprecated since v2.8. Try using nonlinear_pk_at_kvec_and_zvec() instead.\n");

  if (nonlinear == _TRUE_)
    pk_output = pk_nonlinear;
  else
    pk_output = pk_linear;

  class_call(nonlinear_module_->nonlinear_pks_at_kvec_and_zvec(pk_output,
                                                              kvec,
                                                              kvec_size,
                                                              zvec,
                                                              zvec_size,
                                                              pk_tot_out,
                                                              pk_cb_tot_out),
             nonlinear_module_->error_message_,
             error_message_);

  return _SUCCESS_;
}

/**
 * This routine computes sigma(R) given P(k) for total matter power
 * spectrum (does not check that k_max is large enough)
 *
 * This function is deprecated since v2.8. Try using nonlinear_sigmas_at_z() instead.
 *
 * @param pba   Input: pointer to background structure
 * @param ppm   Input: pointer to primordial structure
 * @param psp   Input: pointer to spectra structure
 * @param R     Input: radius in Mpc
 * @param z     Input: redshift
 * @param sigma Output: variance in a sphere of radius R (dimensionless)
 * @return the error status
 */

int SpectraModule::spectra_sigma(double R, double z, double * sigma) {

  fprintf(stderr," -> [WARNING:] You are calling the function spectra_sigma() which is deprecated since v2.8. Try using nonlinear_sigmas_at_z() instead.\n");

  if (nonlinear_module_->has_pk_m_) {

    class_call(nonlinear_module_->nonlinear_sigma_at_z(R,
                                                      z,
                                                      nonlinear_module_->index_pk_m_,
                                                      80., // hardcoded, yes, but the function is deprecated...
                                                      sigma),
               nonlinear_module_->error_message_,
               error_message_);

  }

  return _SUCCESS_;
}

/**
 * This routine computes sigma(R) given P(k) for baryon+cdm power
 * spectrum (does not check that k_max is large enough)
 *
 * This function is deprecated since v2.8. Try using nonlinear_sigmas_at_z() instead.
 *
 * @param pba      Input: pointer to background structure
 * @param ppm      Input: pointer to primordial structure
 * @param psp      Input: pointer to spectra structure
 * @param R        Input: radius in Mpc
 * @param z        Input: redshift
 * @param sigma_cb Output: variance in a sphere of radius R (dimensionless)
 * @return the error status
 */

int SpectraModule::spectra_sigma_cb(double R, double z, double * sigma_cb) {

  fprintf(stderr," -> [WARNING:] You are calling the function spectra_sigma_cb() which is deprecated since v2.8. Try using nonlinear_sigmas_at_z() instead.\n");

  if (nonlinear_module_->has_pk_cb_) {

    class_call(nonlinear_module_->nonlinear_sigma_at_z(R,
                                                      z,
                                                      nonlinear_module_->index_pk_cb_,
                                                      80., // hardcoded, yes, but the function is deprecated...
                                                      sigma_cb),
               nonlinear_module_->error_message_,
               error_message_);
  }

  return _SUCCESS_;
}

  /* deprecated functions (since v2.1) */

/**
 * Obsolete function, superseeded by perturb_sources_at_tau()
 * (at the time of the switch, this function was anyway never used anywhere)
 *
 * @param pba        Input: pointer to background structure (used for converting z into tau)
 * @param psp        Input: pointer to spectra structure (containing pre-computed table)
 * @param z          Input: redshift
 * @param output     Output: matter transfer functions
 * @return the error status
 */

int SpectraModule::spectra_tk_at_z(
                                   double z,
                                   double * output /* array with argument output[(index_k*ic_size_[index_md]+index_ic)*psp->tr_size+index_tr] (must be already allocated) */
                                   ) {

  class_stop(error_message_,
             "The function spectra_tk_at_z() is obsolete, use instead perturb_sources_at_tau(), it does the same");

  return _SUCCESS_;

}

/**
 * Obsolete function, superseeded by perturb_sources_at_tau()
 * (at the time of the switch, this function was anyway never used anywhere)
 *
 * @param pba        Input: pointer to background structure (used for converting z into tau)
 * @param psp        Input: pointer to spectra structure (containing pre-computed table)
 * @param k          Input: wavenumber in 1/Mpc
 * @param z          Input: redshift
 * @param output     Output: matter transfer functions
 * @return the error status
 */

int SpectraModule::spectra_tk_at_k_and_z(
                                         double k,
                                         double z,
                                         double * output  /* array with argument output[index_ic*psp->tr_size+index_tr] (must be already allocated) */
                                         ) {

  class_stop(error_message_,
             "The function spectra_tk_at_k_and_z() is obsolete, use instead perturb_sources_at_tau(), it does the same provided that you interpolate its output at some wavenumber k");

  return _SUCCESS_;

}

  /* end deprecated functions */
