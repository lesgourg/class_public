/** @file primordial.c Documented primordial module.
 *
 * Julien Lesgourgues, 24.08.2010
 *
 * This module computes the primordial spectra. Can be used in different modes:
 * simple parametric form, evolving inflaton perturbations, etc. So far only
 * the mode corresponding to a simple analytic form in terms of amplitudes, tilts
 * and runnings has been developped.
 *
 * The following functions can be called from other modules:
 *
 * -# primordial_init() at the beginning (anytime after perturb_init() and before spectra_init())
 * -# primordial_spectrum_at_k() at any time for computing P(k) at any k
 * -# primordial_free() at the end
 */

#include "primordial.h"

/**
 * Primordial spectra for arbitrary argument and for all initial conditions.
 *
 * This routine evaluates the primordial spectrum at a given value of k by
 * interpolating in the pre-computed table.
 *
 * When k is not in the pre-computed range but the spectrum can be found
 * analytically, finds it. Otherwise returns an error.
 *
 * Can be called in two modes: linear or logarithmic.
 *
 * - linear: takes k, returns P(k)
 *
 * - logarithmic: takes ln(k), return ln(P(k))
 *
 * One little subtlety: in case of several correlated initial conditions,
 * the cross-correlation spectrum can be negative. Then, in logarithmic mode,
 * the non-diagonal elements contain the cross-correlation angle P_12/sqrt(P_11 P_22)
 * (from -1 to 1) instead of ln(P_12)
 *
 * This function can be
 * called from whatever module at whatever time, provided that
 * primordial_init() has been called before, and primordial_free() has not
 * been called yet.
 *
 * @param ppm        Input: pointer to primordial structure containing tabulated primordial spectrum
 * @param index_md Input: index of mode (scalar, tensor, ...)
 * @param mode       Input: linear or logarithmic
 * @param k          Input: wavenumber in 1/Mpc (linear mode) or its logarithm (logarithmic mode)
 * @param pk         Ouput: for each pair of initial conditions, primordial spectra P(k) in Mpc**3 (linear mode), or their logarithms and cross-correlation angles (logarithmic mode)
 * @return the error status
 */

int primordial_spectrum_at_k(
                             struct primordial * ppm,
                             int index_md,
                             enum linear_or_logarithmic mode,
                             double input,
                             double * output /* array with argument output[index_ic1_ic2] (must be already allocated) */
                             ) {

  /** Summary: */

  /** - define local variables */

  int index_ic1,index_ic2,index_ic1_ic2;
  double lnk;
  int last_index;

  /** - infer ln(k) from input. In linear mode, reject negative value of input k value. */

  if (mode == linear) {
    class_test(input<=0.,
               ppm->error_message,
               "k = %e",input);
    lnk=log(input);
  }
  else {
    lnk = input;
  }

  /** - if ln(k) is not in the interpolation range, return an error, unless
      we are in the case of a analytic spectrum, for which a direct computation is possible */

  if ((lnk > ppm->lnk[ppm->lnk_size-1]) || (lnk < ppm->lnk[0])) {

    class_test(ppm->primordial_spec_type != analytic_Pk,
               ppm->error_message,
               "k=%e out of range [%e : %e]",exp(lnk),exp(ppm->lnk[0]),exp(ppm->lnk[ppm->lnk_size-1]));

    /* direct computation */

    for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++) {
      for (index_ic2 = index_ic1; index_ic2 < ppm->ic_size[index_md]; index_ic2++) {

        index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,ppm->ic_size[index_md]);

        if (ppm->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {

          class_call(primordial_analytic_spectrum(ppm,
                                                  index_md,
                                                  index_ic1_ic2,
                                                  exp(lnk),
                                                  &(output[index_ic1_ic2])),
                     ppm->error_message,
                     ppm->error_message);
        }
        else {
          output[index_ic1_ic2] = 0.;
        }
      }
    }

    /* if mode==linear, output is already in the correct format. Otherwise, apply necessary transformation. */

    if (mode == logarithmic) {

      for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++) {
        index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic1,ppm->ic_size[index_md]);
        output[index_ic1_ic2] = log(output[index_ic1_ic2]);
      }
      for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++) {
        for (index_ic2 = index_ic1+1; index_ic2 < ppm->ic_size[index_md]; index_ic2++) {
          index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,ppm->ic_size[index_md]);
          if (ppm->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {
            output[index_ic1_ic2] /= sqrt(output[index_symmetric_matrix(index_ic1,index_ic1,ppm->ic_size[index_md])]*
                                          output[index_symmetric_matrix(index_ic2,index_ic2,ppm->ic_size[index_md])]);
          }
        }
      }
    }
  }

  /** - otherwise, interpolate in the pre-computed table: */

  else {

    class_call(array_interpolate_spline(
                                        ppm->lnk,
                                        ppm->lnk_size,
                                        ppm->lnpk[index_md],
                                        ppm->ddlnpk[index_md],
                                        ppm->ic_ic_size[index_md],
                                        lnk,
                                        &last_index,
                                        output,
                                        ppm->ic_ic_size[index_md],
                                        ppm->error_message),
               ppm->error_message,
               ppm->error_message);

    /* if mode==logarithmic, output is already in the correct format. Otherwise, apply necessary transformation. */

    if (mode == linear) {

      for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++) {
        index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic1,ppm->ic_size[index_md]);
        output[index_ic1_ic2]=exp(output[index_ic1_ic2]);
      }
      for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++) {
        for (index_ic2 = index_ic1+1; index_ic2 < ppm->ic_size[index_md]; index_ic2++) {
          index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,ppm->ic_size[index_md]);
          if (ppm->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {
            output[index_ic1_ic2] *= sqrt(output[index_symmetric_matrix(index_ic1,index_ic1,ppm->ic_size[index_md])]*
                                          output[index_symmetric_matrix(index_ic2,index_ic2,ppm->ic_size[index_md])]);
          }
          else {
            output[index_ic1_ic2] = 0.;
          }
        }
      }
    }
  }

  return _SUCCESS_;

}

/**
 * This routine initializes the primordial structure (in particular, compute table of primordial spectrum values)
 *
 * @param ppr Input : pointer to precision structure (defines method and precision for all computations)
 * @param ppt Input : pointer to perturbation structure (useful for knowing k_min, k_max, etc.)
 * @param ppm Output: pointer to initialized primordial structure
 * @return the error status
 */

int primordial_init(
                    struct precision  * ppr,
                    struct perturbs   * ppt,
                    struct primordial * ppm
                    ) {

  /** Summary: */

  /** - define local variables */

  double k,k_min,k_max;
  int index_md,index_ic1,index_ic2,index_ic1_ic2,index_k;
  double pk,pk1,pk2;
  double dlnk,lnpk_pivot,lnpk_minus,lnpk_plus,lnpk_minusminus,lnpk_plusplus;
  /* uncomment if you use optional test below
     (for correlated isocurvature modes) */
  //double cos_delta_k;

  /** - check that we really need to compute the primordial spectra */

  if (ppt->has_perturbations == _FALSE_) {
    ppm->lnk_size=0;
    if (ppm->primordial_verbose > 0)
      printf("No perturbations requested. Primordial module skipped.\n");
    return _SUCCESS_;
  }
  else {
    if (ppm->primordial_verbose > 0)
      printf("Computing primordial spectra");
  }

  /** - get kmin and kmax from perturbation structure. Test that they make sense. */

  k_min = ppt->k[0]; /* first value, inferred from perturbations structure */
  k_max = ppt->k[ppt->k_size-1]; /* last value, inferred from perturbations structure */

  class_test(k_min <= 0.,
             ppm->error_message,
             "k_min negative or null: stop to avoid segmentation fault");

  class_test(k_max <= 0.,
             ppm->error_message,
             "k_max negative or null: stop to avoid segmentation fault");

  class_test(ppm->k_pivot <= 0.,
             ppm->error_message,
             "k_pivot negative or null: stop to avoid segmentation fault");

  class_test(ppr->k_per_decade_primordial <= 0.,
             ppm->error_message,
             "k_per_decade_primordial negative or null: stop to avoid segmentation fault");

  class_test(ppr->k_per_decade_primordial <= _K_PER_DECADE_PRIMORDIAL_MIN_,
             ppm->error_message,
             "k_per_decade_primordial = %e: you ask for such a sparse sampling of the primordial spectrum that this is probably a mistake",
             ppr->k_per_decade_primordial);

  /** - allocate and fill values of lnk's */

  class_call(primordial_get_lnk_list(ppm,
                                     k_min,
                                     k_max,
                                     ppr->k_per_decade_primordial
                                     ),
             ppm->error_message,
             ppm->error_message);

  /** - define indices and allocate tables in primordial structure */

  class_call(primordial_indices(ppt,
                                ppm),
             ppm->error_message,
             ppm->error_message);

  /** - deal with case of analytic primordial spectra (with amplitudes, tilts, runnings etc.) */

  if (ppm->primordial_spec_type == analytic_Pk) {

    if (ppm->primordial_verbose > 0)
      printf(" (analytic spectrum)\n");

    class_call_except(primordial_analytic_spectrum_init(ppt,
                                                        ppm),
                      ppm->error_message,
                      ppm->error_message,
                      primordial_free(ppm));

    for (index_k = 0; index_k < ppm->lnk_size; index_k++) {

      k=exp(ppm->lnk[index_k]);

      for (index_md = 0; index_md < ppt->md_size; index_md++) {
        for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++) {
          for (index_ic2 = index_ic1; index_ic2 < ppm->ic_size[index_md]; index_ic2++) {

            index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,ppm->ic_size[index_md]);

            if (ppm->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {

              class_call(primordial_analytic_spectrum(ppm,
                                                      index_md,
                                                      index_ic1_ic2,
                                                      k,
                                                      &pk),
                         ppm->error_message,
                         ppm->error_message);

              if (index_ic1 == index_ic2) {

                /* diagonal coefficients: ln[P(k)] */

                ppm->lnpk[index_md][index_k*ppm->ic_ic_size[index_md]+index_ic1_ic2] = log(pk);
              }
              else {

                /* non-diagonal coefficients: cosDelta(k) = P(k)_12/sqrt[P(k)_1 P(k)_2] */

                class_call(primordial_analytic_spectrum(ppm,
                                                        index_md,
                                                        index_symmetric_matrix(index_ic1,index_ic1,ppm->ic_size[index_md]),
                                                        k,
                                                        &pk1),
                           ppm->error_message,
                           ppm->error_message);

                class_call(primordial_analytic_spectrum(ppm,
                                                        index_md,
                                                        index_symmetric_matrix(index_ic2,index_ic2,ppm->ic_size[index_md]),
                                                        k,
                                                        &pk2),
                           ppm->error_message,
                           ppm->error_message);

                /* either return an error if correlation is too large... */
                /*
                  cos_delta_k = pk/sqrt(pk1*pk2);
                  class_test_except((cos_delta_k < -1.) || (cos_delta_k > 1.),
                  ppm->error_message,
                  primordial_free(ppm),
                  "correlation angle between IC's takes unphysical values");

                  ppm->lnpk[index_md][index_k*ppm->ic_ic_size[index_md]+index_ic1_ic2] = cos_delta_k;
                */

                /* ... or enforce definite positive correlation matrix */

                if (pk > sqrt(pk1*pk2))
                  ppm->lnpk[index_md][index_k*ppm->ic_ic_size[index_md]+index_ic1_ic2] = 1.;
                else if (pk < -sqrt(pk1*pk2))
                  ppm->lnpk[index_md][index_k*ppm->ic_ic_size[index_md]+index_ic1_ic2] = -1.;
                else
                  ppm->lnpk[index_md][index_k*ppm->ic_ic_size[index_md]+index_ic1_ic2] = pk/sqrt(pk1*pk2);


              }
            }
            else {

              /* non-diagonal coefficients when ic's are uncorrelated */

              ppm->lnpk[index_md][index_k*ppm->ic_ic_size[index_md]+index_ic1_ic2] = 0.;
            }
          }
        }
      }
    }
  }

  /** - deal with case of inflation with given V(phi) */

  else if (ppm->primordial_spec_type == inflation_V) {

    class_test(ppt->has_scalars == _FALSE_,
               ppm->error_message,
               "inflationary module cannot work if you do not ask for scalar modes");

    class_test(ppt->has_vectors == _TRUE_,
               ppm->error_message,
               "inflationary module cannot work if you ask for vector modes");

    class_test(ppt->has_tensors == _FALSE_,
               ppm->error_message,
               "inflationary module cannot work if you do not ask for tensor modes");

    class_test(ppt->has_bi == _TRUE_ || ppt->has_cdi == _TRUE_ || ppt->has_nid == _TRUE_ || ppt->has_niv == _TRUE_,
               ppm->error_message,
               "inflationary module cannot work if you ask for isocurvature modes");

    class_call(primordial_inflation_indices(ppm),
               ppm->error_message,
               ppm->error_message);

    if (ppm->primordial_verbose > 0)
      printf(" (simulating inflation)\n");

    class_call_except(primordial_inflation_solve_inflation(ppt,ppm,ppr),
                      ppm->error_message,
                      ppm->error_message,
                      primordial_free(ppm));

  }

  /** - deal with the case of external calculation of Pk */

  else if (ppm->primordial_spec_type == external_Pk) {

    class_test(ppt->has_scalars == _FALSE_,
               ppm->error_message,
               "external Pk module cannot work if you do not ask for scalar modes");

    class_test(ppt->has_vectors == _TRUE_,
               ppm->error_message,
               "external Pk module cannot work if you ask for vector modes");

    class_test(ppt->has_bi == _TRUE_ || ppt->has_cdi == _TRUE_ || ppt->has_nid == _TRUE_ || ppt->has_niv == _TRUE_,
               ppm->error_message,
               "external Pk module cannot work if you ask for isocurvature modes (but that could be implemented easily in the future!)");

    if (ppm->primordial_verbose > 0)
      printf(" (Pk calculated externally)\n");

    class_call_except(primordial_external_spectrum_init(ppt,ppm),
                      ppm->error_message,
                      ppm->error_message,
                      primordial_free(ppm));
  }

  else {

    class_test(0==0,
               ppm->error_message,
               "only analytic, external and inflation_V primordial spectrum coded yet");

  }

  /** - compute second derivative of each lnpk versus lnk  with spline, in view of interpolation */

  for (index_md = 0; index_md < ppm->md_size; index_md++) {

    class_call(array_spline_table_lines(ppm->lnk,
                                        ppm->lnk_size,
                                        ppm->lnpk[index_md],
                                        ppm->ic_ic_size[index_md],
                                        ppm->ddlnpk[index_md],
                                        _SPLINE_EST_DERIV_,
                                        ppm->error_message),
               ppm->error_message,
               ppm->error_message);

  }

  /** derive spectral parameters from numerically computed spectra
      (not used by the rest of the code, but useful to keep in memory for several types of investigations) */

  if (ppm->primordial_spec_type != analytic_Pk) {

    dlnk = log(10.)/ppr->k_per_decade_primordial;

    if (ppt->has_scalars == _TRUE_) {

      class_call(primordial_spectrum_at_k(ppm,
                                          ppt->index_md_scalars,
                                          logarithmic,
                                          log(ppm->k_pivot),
                                          &lnpk_pivot),
                 ppm->error_message,
                 ppm->error_message);

      class_call(primordial_spectrum_at_k(ppm,
                                          ppt->index_md_scalars,
                                          logarithmic,
                                          log(ppm->k_pivot)+dlnk,

                                          &lnpk_plus),
                 ppm->error_message,
                 ppm->error_message);

      class_call(primordial_spectrum_at_k(ppm,
                                          ppt->index_md_scalars,
                                          logarithmic,
                                          log(ppm->k_pivot)-dlnk,
                                          &lnpk_minus),
                 ppm->error_message,
                 ppm->error_message);

      ppm->A_s = exp(lnpk_pivot);
      ppm->n_s = (lnpk_plus-lnpk_minus)/(2.*dlnk)+1.;
      ppm->alpha_s = (lnpk_plus-2.*lnpk_pivot+lnpk_minus)/pow(dlnk,2);

      /** expression for n_s comes from:

          ns_2 = (lnpk_plus-lnpk_pivot)/(dlnk)+1.
          ns_1 = (lnpk_pivot-lnpk_minus)/(dlnk)+1.
          alpha_s = dns/dlnk
          = (ns_2-ns_1)/dlnk
          = (lnpk_plus-lnpk_pivot-lnpk_pivot+lnpk_minus)/(dlnk)/(dlnk)

      **/

      class_call(primordial_spectrum_at_k(ppm,
                                          ppt->index_md_scalars,
                                          logarithmic,
                                          log(ppm->k_pivot)+2.*dlnk,

                                          &lnpk_plusplus),
                 ppm->error_message,
                 ppm->error_message);

      class_call(primordial_spectrum_at_k(ppm,
                                          ppt->index_md_scalars,
                                          logarithmic,
                                          log(ppm->k_pivot)-2.*dlnk,
                                          &lnpk_minusminus),
                 ppm->error_message,
                 ppm->error_message);

      /** expression for beta_s:

          ppm->beta_s = (alpha_plus-alpha_minus)/dlnk
          = (lnpk_plusplus-2.*lnpk_plus+lnpk_pivot - (lnpk_pivot-2.*lnpk_minus+lnpk_minusminus)/pow(dlnk,3);

          This simplifies into:

      **/

      ppm->beta_s = (lnpk_plusplus-2.*lnpk_plus+2.*lnpk_minus-lnpk_minusminus)/pow(dlnk,3);

      if (ppm->primordial_verbose > 0)
        printf(" -> A_s=%g  n_s=%g  alpha_s=%g\n",ppm->A_s,ppm->n_s,ppm->alpha_s);

    }

    if (ppt->has_tensors == _TRUE_) {

      class_call(primordial_spectrum_at_k(ppm,
                                          ppt->index_md_tensors,
                                          logarithmic,
                                          log(ppm->k_pivot),
                                          &lnpk_pivot),
                 ppm->error_message,
                 ppm->error_message);

      class_call(primordial_spectrum_at_k(ppm,
                                          ppt->index_md_tensors,
                                          logarithmic,
                                          log(ppm->k_pivot)+dlnk,
                                          &lnpk_plus),
                 ppm->error_message,
                 ppm->error_message);

      class_call(primordial_spectrum_at_k(ppm,
                                          ppt->index_md_tensors,
                                          logarithmic,
                                          log(ppm->k_pivot)-dlnk,
                                          &lnpk_minus),
                 ppm->error_message,
                 ppm->error_message);

      ppm->r = exp(lnpk_pivot)/ppm->A_s;
      ppm->n_t = (lnpk_plus-lnpk_minus)/(2.*dlnk);
      ppm->alpha_t = (lnpk_plus-2.*lnpk_pivot+lnpk_minus)/pow(dlnk,2);

      if (ppm->primordial_verbose > 0)
        printf(" -> r=%g  n_r=%g  alpha_r=%g\n",ppm->r,ppm->n_t,ppm->alpha_t);

    }

  }

  return _SUCCESS_;

}

/**
 * This routine frees all the memory space allocated by primordial_init().
 *
 * To be called at the end of each run.
 *
 * @param ppm Input: pointer to primordial structure (which fields must be freed)
 * @return the error status
 */

int primordial_free(
                    struct primordial * ppm
                    ) {

  int index_md;

  if (ppm->lnk_size > 0) {

    if (ppm->primordial_spec_type == analytic_Pk) {
      for (index_md = 0; index_md < ppm->md_size; index_md++) {
        free(ppm->amplitude[index_md]);
        free(ppm->tilt[index_md]);
        free(ppm->running[index_md]);
      }
      free(ppm->amplitude);
      free(ppm->tilt);
      free(ppm->running);
    }
    else if (ppm->primordial_spec_type == external_Pk) {
      free(ppm->command);
    }

    for (index_md = 0; index_md < ppm->md_size; index_md++) {
      free(ppm->lnpk[index_md]);
      free(ppm->ddlnpk[index_md]);
      free(ppm->is_non_zero[index_md]);
    }

    free(ppm->lnpk);
    free(ppm->ddlnpk);
    free(ppm->is_non_zero);
    free(ppm->ic_size);
    free(ppm->ic_ic_size);

    free(ppm->lnk);

  }

  return _SUCCESS_;
}

/**
 * This routine defines indices and allocates tables in the primordial structure
 *
 * @param ppt  Input : pointer to perturbation structure
 * @param ppm  Input/output: pointer to primordial structure
 * @return the error status
 */

int primordial_indices(
                       struct perturbs   * ppt,
                       struct primordial * ppm
                       ) {

  int index_md;

  ppm->md_size = ppt->md_size;

  class_alloc(ppm->lnpk,ppt->md_size*sizeof(double*),ppm->error_message);

  class_alloc(ppm->ddlnpk,ppt->md_size*sizeof(double*),ppm->error_message);

  class_alloc(ppm->ic_size,ppt->md_size*sizeof(int*),ppm->error_message);

  class_alloc(ppm->ic_ic_size,ppt->md_size*sizeof(int*),ppm->error_message);

  class_alloc(ppm->is_non_zero,ppm->md_size*sizeof(short *),ppm->error_message);

  for (index_md = 0; index_md < ppt->md_size; index_md++) {

    ppm->ic_size[index_md] = ppt->ic_size[index_md];

    ppm->ic_ic_size[index_md] = (ppm->ic_size[index_md]*(ppm->ic_size[index_md]+1))/2;

    class_alloc(ppm->lnpk[index_md],
                ppm->lnk_size*ppm->ic_ic_size[index_md]*sizeof(double),
                ppm->error_message);

    class_alloc(ppm->ddlnpk[index_md],
                ppm->lnk_size*ppm->ic_ic_size[index_md]*sizeof(double),
                ppm->error_message);

    class_alloc(ppm->is_non_zero[index_md],
                ppm->ic_ic_size[index_md]*sizeof(short),
                ppm->error_message);


  }

  return _SUCCESS_;

}

/**
 * This routine allocates and fills the list of wavenumbers k
 *
 *
 * @param ppm  Input/output: pointer to primordial structure
 * @param kmin Input : first value
 * @param kmax Input : last value that we should encompass
 * @param k_per_decade Input : number of k per decade
 * @return the error status
 */

int primordial_get_lnk_list(
                            struct primordial * ppm,
                            double kmin,
                            double kmax,
                            double k_per_decade
                            ) {

  int i;

  class_test((kmin <= 0.) || (kmax <= kmin),
             ppm->error_message,
             "inconsistent values of kmin=%e, kmax=%e",kmin,kmax);

  ppm->lnk_size = (int)(log(kmax/kmin)/log(10.)*k_per_decade) + 2;

  class_alloc(ppm->lnk,ppm->lnk_size*sizeof(double),ppm->error_message);

  for (i=0; i<ppm->lnk_size; i++)
    ppm->lnk[i]=log(kmin)+i*log(10.)/k_per_decade;

  return _SUCCESS_;

}

/**
 * This routine interprets and stores in a condensed form the input parameters
 * in the case of a simple analytic spectra with amplitudes, tilts, runnings,
 * in such way that later on, the spectrum can be obtained by a quick call to
 * the routine primordial_analytic_spectrum(()
 *
 * @param ppt  Input : pointer to perturbation structure
 * @param ppm  Input/output: pointer to primordial structure
 * @return the error status
 */

int primordial_analytic_spectrum_init(
                                      struct perturbs   * ppt,
                                      struct primordial * ppm
                                      ) {

  int index_md,index_ic1,index_ic2;
  int index_ic1_ic2,index_ic1_ic1,index_ic2_ic2;
  double one_amplitude=0.;
  double one_tilt=0.;
  double one_running=0.;
  double one_correlation=0.;

  class_alloc(ppm->amplitude,
              ppm->md_size*sizeof(double *),
              ppm->error_message);

  class_alloc(ppm->tilt,
              ppm->md_size*sizeof(double *),
              ppm->error_message);

  class_alloc(ppm->running,
              ppm->md_size*sizeof(double *),
              ppm->error_message);

  for (index_md = 0; index_md < ppm->md_size; index_md++) {

    class_alloc(ppm->amplitude[index_md],
                ppm->ic_ic_size[index_md]*sizeof(double),
                ppm->error_message);

    class_alloc(ppm->tilt[index_md],
                ppm->ic_ic_size[index_md]*sizeof(double),
                ppm->error_message);

    class_alloc(ppm->running[index_md],
                ppm->ic_ic_size[index_md]*sizeof(double),
                ppm->error_message);

  }

  for (index_md = 0; index_md < ppm->md_size; index_md++) {

    /* diagonal coefficients */

    for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++) {

      if (_scalars_) {

        if ((ppt->has_ad == _TRUE_) && (index_ic1 == ppt->index_ic_ad)) {
          one_amplitude = ppm->A_s;
          one_tilt = ppm->n_s;
          one_running = ppm->alpha_s;
        }

        if ((ppt->has_bi == _TRUE_) && (index_ic1 == ppt->index_ic_bi)) {
          one_amplitude = ppm->A_s*ppm->f_bi*ppm->f_bi;
          one_tilt = ppm->n_bi;
          one_running = ppm->alpha_bi;
        }

        if ((ppt->has_cdi == _TRUE_) && (index_ic1 == ppt->index_ic_cdi)) {
          one_amplitude = ppm->A_s*ppm->f_cdi*ppm->f_cdi;
          one_tilt = ppm->n_cdi;
          one_running = ppm->alpha_cdi;
        }

        if ((ppt->has_nid == _TRUE_) && (index_ic1 == ppt->index_ic_nid)) {
          one_amplitude = ppm->A_s*ppm->f_nid*ppm->f_nid;
          one_tilt = ppm->n_nid;
          one_running = ppm->alpha_nid;
        }

        if ((ppt->has_niv == _TRUE_) && (index_ic1 == ppt->index_ic_niv)) {
          one_amplitude = ppm->A_s*ppm->f_niv*ppm->f_niv;
          one_tilt = ppm->n_niv;
          one_running = ppm->alpha_niv;
        }
      }

      if (_tensors_) {

        if (index_ic1 == ppt->index_ic_ten) {
          one_amplitude = ppm->A_s*ppm->r;
          one_tilt = ppm->n_t+1.; /* +1 to match usual definition of n_t (equivalent to n_s-1) */
          one_running = ppm->alpha_t;
        }
      }

      class_test(one_amplitude <= 0.,
                 ppm->error_message,
                 "inconsistent input for primordial amplitude: %g for index_md=%d, index_ic=%d\n",
                 one_amplitude,index_md,index_ic1);

      index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic1,ppm->ic_size[index_md]);

      ppm->is_non_zero[index_md][index_ic1_ic2] = _TRUE_;
      ppm->amplitude[index_md][index_ic1_ic2] = one_amplitude;
      ppm->tilt[index_md][index_ic1_ic2] = one_tilt;
      ppm->running[index_md][index_ic1_ic2] = one_running;
    }

    /* non-diagonal coefficients */

    for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++) {
      for (index_ic2 = index_ic1+1; index_ic2 < ppm->ic_size[index_md]; index_ic2++) {

        if (_scalars_) {

          if ((ppt->has_ad == _TRUE_) && (ppt->has_bi == _TRUE_) &&
              (((index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_bi)) ||
               ((index_ic1 == ppt->index_ic_ad) && (index_ic1 == ppt->index_ic_bi)))) {
            one_correlation = ppm->c_ad_bi;
            one_tilt = ppm->n_ad_bi;
            one_running = ppm->alpha_ad_bi;
          }

          if ((ppt->has_ad == _TRUE_) && (ppt->has_cdi == _TRUE_) &&
              (((index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_cdi)) ||
               ((index_ic2 == ppt->index_ic_ad) && (index_ic1 == ppt->index_ic_cdi)))) {
            one_correlation = ppm->c_ad_cdi;
            one_tilt = ppm->n_ad_cdi;
            one_running = ppm->alpha_ad_cdi;
          }

          if ((ppt->has_ad == _TRUE_) && (ppt->has_nid == _TRUE_) &&
              (((index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_nid)) ||
               ((index_ic2 == ppt->index_ic_ad) && (index_ic1 == ppt->index_ic_nid)))) {
            one_correlation = ppm->c_ad_nid;
            one_tilt = ppm->n_ad_nid;
            one_running = ppm->alpha_ad_nid;
          }

          if ((ppt->has_ad == _TRUE_) && (ppt->has_niv == _TRUE_) &&
              (((index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_niv)) ||
               ((index_ic2 == ppt->index_ic_ad) && (index_ic1 == ppt->index_ic_niv)))) {
            one_correlation = ppm->c_ad_niv;
            one_tilt = ppm->n_ad_niv;
            one_running = ppm->alpha_ad_niv;
          }

          if ((ppt->has_bi == _TRUE_) && (ppt->has_cdi == _TRUE_) &&
              (((index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_cdi)) ||
               ((index_ic2 == ppt->index_ic_bi) && (index_ic1 == ppt->index_ic_cdi)))) {
            one_correlation = ppm->c_bi_cdi;
            one_tilt = ppm->n_bi_cdi;
            one_running = ppm->alpha_bi_cdi;
          }

          if ((ppt->has_bi == _TRUE_) && (ppt->has_nid == _TRUE_) &&
              (((index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_nid)) ||
               ((index_ic2 == ppt->index_ic_bi) && (index_ic1 == ppt->index_ic_nid)))) {
            one_correlation = ppm->c_bi_nid;
            one_tilt = ppm->n_bi_nid;
            one_running = ppm->alpha_bi_nid;
          }

          if ((ppt->has_bi == _TRUE_) && (ppt->has_niv == _TRUE_) &&
              (((index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_niv)) ||
               ((index_ic2 == ppt->index_ic_bi) && (index_ic1 == ppt->index_ic_niv)))) {
            one_correlation = ppm->c_bi_niv;
            one_tilt = ppm->n_bi_niv;
            one_running = ppm->alpha_bi_niv;
          }

          if ((ppt->has_cdi == _TRUE_) && (ppt->has_nid == _TRUE_) &&
              (((index_ic1 == ppt->index_ic_cdi) && (index_ic2 == ppt->index_ic_nid)) ||
               ((index_ic2 == ppt->index_ic_cdi) && (index_ic1 == ppt->index_ic_nid)))) {
            one_correlation = ppm->c_cdi_nid;
            one_tilt = ppm->n_cdi_nid;
            one_running = ppm->alpha_cdi_nid;
          }

          if ((ppt->has_cdi == _TRUE_) && (ppt->has_niv == _TRUE_) &&
              (((index_ic1 == ppt->index_ic_cdi) && (index_ic2 == ppt->index_ic_niv)) ||
               ((index_ic2 == ppt->index_ic_cdi) && (index_ic1 == ppt->index_ic_niv)))) {
            one_correlation = ppm->c_cdi_niv;
            one_tilt = ppm->n_cdi_niv;
            one_running = ppm->alpha_cdi_niv;
          }

          if ((ppt->has_nid == _TRUE_) && (ppt->has_niv == _TRUE_) &&
              (((index_ic1 == ppt->index_ic_nid) && (index_ic2 == ppt->index_ic_niv)) ||
               ((index_ic2 == ppt->index_ic_nid) && (index_ic1 == ppt->index_ic_niv)))) {
            one_correlation = ppm->c_nid_niv;
            one_tilt = ppm->n_nid_niv;
            one_running = ppm->alpha_nid_niv;
          }

        }

        class_test((one_correlation < -1) || (one_correlation > 1),
                   ppm->error_message,
                   "inconsistent input for isocurvature cross-correlation\n");

        index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,ppm->ic_size[index_md]);
        index_ic1_ic1 = index_symmetric_matrix(index_ic1,index_ic1,ppm->ic_size[index_md]);
        index_ic2_ic2 = index_symmetric_matrix(index_ic2,index_ic2,ppm->ic_size[index_md]);

        if (one_correlation == 0.) {
          ppm->is_non_zero[index_md][index_ic1_ic2] = _FALSE_;
          ppm->amplitude[index_md][index_ic1_ic2] = 0.;
          ppm->tilt[index_md][index_ic1_ic2] = 0.;
          ppm->running[index_md][index_ic1_ic2] = 0.;
        }
        else {
          ppm->is_non_zero[index_md][index_ic1_ic2] = _TRUE_;
          ppm->amplitude[index_md][index_ic1_ic2] =
            sqrt(ppm->amplitude[index_md][index_ic1_ic1]*
                 ppm->amplitude[index_md][index_ic2_ic2])*
            one_correlation;
          ppm->tilt[index_md][index_ic1_ic2] =
            0.5*(ppm->tilt[index_md][index_ic1_ic1]
                 +ppm->tilt[index_md][index_ic2_ic2])
            + one_tilt;
          ppm->running[index_md][index_ic1_ic2] =
            0.5*(ppm->running[index_md][index_ic1_ic1]
                 +ppm->running[index_md][index_ic2_ic2])
            + one_running;
        }
      }
    }
  }

  return _SUCCESS_;

}

/**
 * This routine returns the primordial spectrum in the simple analytic case with
 * amplitudes, tilts, runnings, for each mode (scalar/tensor...),
 * pair of initial conditions, and wavenumber.
 *
 * @param ppm            Input/output: pointer to primordial structure
 * @param index_md     Input: index of mode (scalar, tensor, ...)
 * @param index_ic1_ic2  Input: pair of initial conditions (ic1, ic2)
 * @param k              Input: wavenumber in same units as pivot scale, i.e. in 1/Mpc
 * @param pk             Output: primordial power spectrum A (k/k_pivot)^(n+...)
 * @return the error status
 */

int primordial_analytic_spectrum(
                                 struct primordial * ppm,
                                 int index_md,
                                 int index_ic1_ic2,
                                 double k,
                                 double * pk
                                 ) {

  if (ppm->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {
    *pk = ppm->amplitude[index_md][index_ic1_ic2]
      *exp((ppm->tilt[index_md][index_ic1_ic2]-1.)*log(k/ppm->k_pivot)
           + 0.5 * ppm->running[index_md][index_ic1_ic2] * pow(log(k/ppm->k_pivot), 2.));

  }
  else {
    *pk = 0.;
  }

  return _SUCCESS_;

}

/**
 * This routine encodes the inflaton scalar potential
 *
 * @param ppm            Input: pointer to primordial structure
 * @param phi            Input: background inflaton field value in units of Mp
 * @param V              Output: inflaton potential in units of MP^4
 * @param dV             Output: first derivative of inflaton potential wrt the field
 * @param ddV            Output: second derivative of inflaton potential wrt the field
 * @return the error status
 */

int primordial_inflation_potential(
                                   struct primordial * ppm,
                                   double phi,
                                   double * V,
                                   double * dV,
                                   double * ddV
                                   ) {

  /* V(phi)=polynomial in (phi-phi*) */
  if (ppm->potential == polynomial) {

    *V   = ppm->V0+(phi-ppm->phi_pivot)*ppm->V1+pow((phi-ppm->phi_pivot),2)/2.*ppm->V2+pow((phi-ppm->phi_pivot),3)/6.*ppm->V3+pow((phi-ppm->phi_pivot),4)/24.*ppm->V4;
    *dV  = ppm->V1+(phi-ppm->phi_pivot)*ppm->V2+pow((phi-ppm->phi_pivot),2)/2.*ppm->V3+pow((phi-ppm->phi_pivot),3)/6.*ppm->V4;
    *ddV = ppm->V2+(phi-ppm->phi_pivot)*ppm->V3+pow((phi-ppm->phi_pivot),2)/2.*ppm->V4;

  }

  /* V(phi) = Lambda^4(1+cos(phi/f)) = V0 (1+cos(phi/V1)) */
  if (ppm->potential == natural) {

    *V   = ppm->V0*(1.+cos(phi/ppm->V1));
    *dV  = -ppm->V0/ppm->V1*sin(phi/ppm->V1);
    *ddV = -ppm->V0/ppm->V1/ppm->V1*cos(phi/ppm->V1);

  }

  /* code here other shapes */

  return _SUCCESS_;
}

/**
 * This routine defines indices used by the inflation simulator
 *
 * @param ppm  Input/output: pointer to primordial structure
 * @return the error status
 */
int primordial_inflation_indices(
                                 struct primordial * ppm
                                 ) {

  int index_in;

  index_in = 0;

  /* indices for background quantitites */
  ppm->index_in_a = index_in;
  index_in ++;
  ppm->index_in_phi = index_in;
  index_in ++;
  ppm->index_in_dphi = index_in;
  index_in ++;

  /* size of background vector */
  ppm->in_bg_size = index_in;

  /* indices for perturbations */
  ppm->index_in_ksi_re = index_in;
  index_in ++;
  ppm->index_in_ksi_im = index_in;
  index_in ++;
  ppm->index_in_dksi_re = index_in;
  index_in ++;
  ppm->index_in_dksi_im = index_in;
  index_in ++;
  ppm->index_in_ah_re = index_in;
  index_in ++;
  ppm->index_in_ah_im = index_in;
  index_in ++;
  ppm->index_in_dah_re = index_in;
  index_in ++;
  ppm->index_in_dah_im = index_in;
  index_in ++;

  /* size of perturbation vector */
  ppm->in_size = index_in;

  return _SUCCESS_;
}

/**
 * Main routine of inflation simulator. Its goal is to check the
 * background evolution before and after the pivot value
 * phi=phi_pivot, and then, if this evolution is suitable, to call the
 * routine primordial_inflation_spectra().
 *
 * @param ppt  Input: pointer to perturbation structure
 * @param ppm  Input/output: pointer to primordial structure
 * @param ppr  Input: pointer to precision structure
 * @return the error status
 */

int primordial_inflation_solve_inflation(
                                         struct perturbs * ppt,
                                         struct primordial * ppm,
                                         struct precision *ppr
                                         ) {

  double * y;
  double * y_ini;
  double * dy;
  double a_pivot,a_try;
  double H_pivot,H_try;
  double phi_try;
  double dphidt_pivot,dphidt_try;
  double aH_ini;
  double k_max;
  int counter;
  double V,dV,ddV;

  //  fprintf(stdout,"Expected slow-roll A_s: %g\n",128.*_PI_/3.*pow(ppm->V0,3)/pow(ppm->V1,2));
  //  fprintf(stdout,"Expected slow-roll T/S: %g\n",pow(ppm->V1/ppm->V0,2)/_PI_);
  //  fprintf(stdout,"Expected slow-roll A_T: %g\n",pow(ppm->V1/ppm->V0,2)/_PI_*128.*_PI_/3.*pow(ppm->V0,3)/pow(ppm->V1,2));
  //  fprintf(stdout,"Expected slow-roll n_s: %g\n",1.-6./16./_PI_*pow(ppm->V1/ppm->V0,2)+2./8./_PI_*(ppm->V2/ppm->V0));
  //  fprintf(stdout,"Expected slow-roll n_t: %g\n",-2./16./_PI_*pow(ppm->V1/ppm->V0,2));

  /* allocate vectors for background/perturbed quantitites */
  class_alloc(y,ppm->in_size*sizeof(double),ppm->error_message);
  class_alloc(y_ini,ppm->in_size*sizeof(double),ppm->error_message);
  class_alloc(dy,ppm->in_size*sizeof(double),ppm->error_message);

  /* check positivity and negative slope of potential in field pivot value */
  class_call_except(primordial_inflation_check_potential(ppm,ppm->phi_pivot),
                    ppm->error_message,
                    ppm->error_message,
                    free(y);free(y_ini);free(dy));

  if (ppm->primordial_verbose > 1)
    printf(" (search attractor at pivot)\n");

  /* find value of phi_dot and H for field's pivot value, assuming slow-roll
     attractor solution has been reached. If no solution, code will
     stop there. */
  class_call_except(primordial_inflation_find_attractor(ppm,
                                                        ppr,
                                                        ppm->phi_pivot,
                                                        ppr->primordial_inflation_attractor_precision_pivot,
                                                        y,
                                                        dy,
                                                        &H_pivot,
                                                        &dphidt_pivot),
                    ppm->error_message,
                    ppm->error_message,
                    free(y);free(y_ini);free(dy));

  /* find a_pivot, value of scale factor when k_pivot crosses horizon while phi=phi_pivot */
  a_pivot = ppm->k_pivot/H_pivot;

  /* integrate background solution starting from phi_pivot and until
     k_max>>aH. This ensure that the inflationary model considered
     here is valid and that the primordial spectrum can be
     computed. Otherwise, if slow-roll brakes too early, model is not
     suitable and run stops. */
  k_max = exp(ppm->lnk[ppm->lnk_size-1]);
  y[ppm->index_in_a] = a_pivot;
  y[ppm->index_in_phi] = ppm->phi_pivot;
  y[ppm->index_in_dphi] = a_pivot*dphidt_pivot;

  if (ppm->primordial_verbose > 1)
    printf(" (check inflation duration after pivot)\n");

  class_call_except(primordial_inflation_reach_aH(ppm,
                                                  ppr,
                                                  y,
                                                  dy,
                                                  k_max/ppr->primordial_inflation_ratio_max),
                    ppm->error_message,
                    ppm->error_message,
                    free(y);free(y_ini);free(dy));

  /* we need to do the opposite: to check that there is an initial
     time such that k_min << (aH)_ini. One such time is found by
     iterations. If no solution exist (no long-enough slow-roll period
     before the pivot scale), the run stops. */
  aH_ini = exp(ppm->lnk[0])/ppr->primordial_inflation_ratio_min;

  a_try = a_pivot;
  H_try = H_pivot;
  phi_try = ppm->phi_pivot;
  counter = 0;

  if (ppm->primordial_verbose > 1)
    printf(" (check inflation duration before pivot, with phi_pivot=%e)\n",phi_try);

  while ((a_try*H_try) >= aH_ini) {

    counter ++;

    class_test_except(counter >= ppr->primordial_inflation_phi_ini_maxit,
                      ppm->error_message,
                      free(y);free(y_ini);free(dy),
                      "when searching for an initial value of phi just before observable inflation takes place, could not converge after %d iterations. The potential does not allow eough inflationary e-folds before reaching the pivot scale",
                      counter);

    class_call_except(primordial_inflation_potential(ppm,phi_try,&V,&dV,&ddV),
                      ppm->error_message,
                      ppm->error_message,
                      free(y);free(y_ini);free(dy));

    phi_try += ppr->primordial_inflation_jump_initial*log(a_try*H_try/aH_ini)*dV/V/8./_PI_;

    printf(" (--> search attractor at phi_try=%e)\n",phi_try);

    class_call_except(primordial_inflation_find_attractor(ppm,
                                                          ppr,
                                                          phi_try,
                                                          ppr->primordial_inflation_attractor_precision_initial,
                                                          y,
                                                          dy,
                                                          &H_try,
                                                          &dphidt_try),
                      ppm->error_message,
                      ppm->error_message,
                      free(y);free(y_ini);free(dy));

    y[ppm->index_in_a] = 1.;
    y[ppm->index_in_phi] = phi_try;
    y[ppm->index_in_dphi] = y[ppm->index_in_a]*dphidt_try;

    if (ppm->primordial_verbose > 1)
      printf(" (--> compute e-folds from phi_try=%e to phi_pivot=%e with dphi/dt_try=%e)\n",phi_try,ppm->phi_pivot,dphidt_try);

    class_call_except(primordial_inflation_evolve_background(ppm,
                                                             ppr,
                                                             y,
                                                             dy,
                                                             ppm->phi_pivot),
                      ppm->error_message,
                      ppm->error_message,
                      free(y);free(y_ini);free(dy));

    a_try = a_pivot/y[ppm->index_in_a];

    if (ppm->primordial_verbose > 1)
      printf(" (--> found %f e-folds\n",-log(a_try));

  }

  /* we found an initial time labeled 'try' with a_try < a_ini, and
     such that an inflationary attractor solution exists. We
     initialize background quantitites at this time. */
  y_ini[ppm->index_in_a] = a_try;
  y_ini[ppm->index_in_phi] = phi_try;
  y_ini[ppm->index_in_dphi] = a_try*dphidt_try;

  if (ppm->primordial_verbose > 1)
    printf(" (compute spectrum)\n");

  /* statting from this time, we run the routine which takes care of computing the primordial spectrum. */
  class_call_except(primordial_inflation_spectra(ppt,
                                                 ppm,
                                                 ppr,
                                                 y_ini,
                                                 y,
                                                 dy),
                    ppm->error_message,
                    ppm->error_message,
                    free(y);free(y_ini);free(dy));

  /* before ending, we want to compute and store the values of phi correspondig to k=aH for k_min and k_max */

  y[ppm->index_in_a] = y_ini[ppm->index_in_a];
  y[ppm->index_in_phi] = y_ini[ppm->index_in_phi];
  y[ppm->index_in_dphi] = y_ini[ppm->index_in_dphi];

  class_call_except(primordial_inflation_reach_aH(ppm,ppr,y,dy,exp(ppm->lnk[0])),
                    ppm->error_message,
                    ppm->error_message,
                    free(y);free(y_ini);free(dy));

  ppm->phi_min=y[ppm->index_in_phi];

  class_call_except(primordial_inflation_reach_aH(ppm,ppr,y,dy,exp(ppm->lnk[ppm->lnk_size])),
                    ppm->error_message,
                    ppm->error_message,
                    free(y);free(y_ini);free(dy));

  ppm->phi_max=y[ppm->index_in_phi];

  if (ppm->primordial_verbose > 1)
    printf(" (observable power spectrum goes from %e to %e)\n",ppm->phi_min,ppm->phi_max);

  /* we are done, we can de-allocate */

  free(y);
  free(y_ini);
  free(dy);

  return _SUCCESS_;
}

/**
 * Routine coordinating the computation of the primordial
 * spectrum. For each wavenumber it calls primordial_inflation_one_k() to
 * integrate the perturbation equations, and then it stores the result
 * for the scalar/tensor spectra.
 *
 * @param ppt   Input: pointer to perturbation structure
 * @param ppm   Input/output: pointer to primordial structure
 * @param ppr   Input: pointer to precision structure
 * @param y_ini Input: initial conditions for the vector of background/perturbations, already allocated and filled
 * @param y     Input: running vector of background/perturbations, already allocated
 * @param dy    Input: running vector of background/perturbation derivatives, already allocated
 * @return the error status
 */

int primordial_inflation_spectra(
                                 struct perturbs * ppt,
                                 struct primordial * ppm,
                                 struct precision * ppr,
                                 double * y_ini,
                                 double * y,
                                 double * dy
                                 ) {
  double aH,k;
  int index_k;
  double V,dV,ddV;
  double curvature,tensors;

  /* check positivity and negative slope of potential in initial field value */
  class_call(primordial_inflation_check_potential(ppm,y_ini[ppm->index_in_phi]),
             ppm->error_message,
             ppm->error_message);

  /* get scalar potential in initial field value */
  class_call(primordial_inflation_potential(ppm,y_ini[ppm->index_in_phi],&V,&dV,&ddV),
             ppm->error_message,
             ppm->error_message);

  /* get initial aH from Friedmann equation */
  aH = sqrt((8*_PI_/3.)*(0.5*y_ini[ppm->index_in_dphi]*y_ini[ppm->index_in_dphi]+y_ini[ppm->index_in_a]*y_ini[ppm->index_in_a]*V));

  class_test(aH >= exp(ppm->lnk[0])/ppr->primordial_inflation_ratio_min,
             ppm->error_message,
             "at initial time, a_k_min > a*H*ratio_min");

  /* loop over Fourier wavenumbers */
  for (index_k=0; index_k < ppm->lnk_size; index_k++) {

    k = exp(ppm->lnk[index_k]);

    /* initialize the background part of the running vector */
    y[ppm->index_in_a] = y_ini[ppm->index_in_a];
    y[ppm->index_in_phi] = y_ini[ppm->index_in_phi];
    y[ppm->index_in_dphi] = y_ini[ppm->index_in_dphi];

    /* evolve the background until the relevant initial time for integrating perturbations */
    class_call(primordial_inflation_reach_aH(ppm,ppr,y,dy,k/ppr->primordial_inflation_ratio_min),
               ppm->error_message,
               ppm->error_message);

    /* evolve the background/perturbation equations from this time and until some time fater Horizon crossing */
    class_call(primordial_inflation_one_k(ppm,ppr,k,y,dy,&curvature,&tensors),
               ppm->error_message,
               ppm->error_message);

    class_test(curvature<=0.,
               ppm->error_message,
               "negative curvature spectrum");

    class_test(tensors<=0.,
               ppm->error_message,
               "negative tensor spectrum");

    /* store the obtained result for curvatute and tensor perturbations */
    ppm->lnpk[ppt->index_md_scalars][index_k] = log(curvature);
    ppm->lnpk[ppt->index_md_tensors][index_k] = log(tensors);

    /* fprintf(stderr,"%e %e %e\n", */
    /* 	    ppm->lnk[index_k], */
    /* 	    ppm->lnpk[ppt->index_md_scalars][index_k], */
    /* 	    ppm->lnpk[ppt->index_md_tensors][index_k]); */

  }

  ppm->is_non_zero[ppt->index_md_scalars][ppt->index_ic_ad] = _TRUE_;
  ppm->is_non_zero[ppt->index_md_tensors][ppt->index_ic_ten] = _TRUE_;

  return _SUCCESS_;

}

/**
 * Routine integrating the background plus perturbation equations for
 * or each wavenumber, and returning the scalr and tensor spectrum.
 *
 * @param ppm   Input: pointer to primordial structure
 * @param ppr   Input: pointer to precision structure
 * @param k     Input: Fourier wavenumber
 * @param y     Input: running vector of background/perturbations, already allocated and initialized
 * @param dy    Input: running vector of background/perturbation derivatives, already allocated
 * @param curvature  Output: curvature perturbation
 * @param tensor     Output: tensor perturbation
 * @return the error status
 */

int primordial_inflation_one_k(
                               struct primordial * ppm,
                               struct precision * ppr,
                               double k,
                               double * y,
                               double * dy,
                               double * curvature,
                               double * tensor
                               ) {

  double tau_start,tau_end,dtau;
  double z,ksi2,ah2;
  double aH;
  double curvature_old;
  double curvature_new;
  double dlnPdN;

  struct primordial_inflation_parameters_and_workspace pipaw;
  struct generic_integrator_workspace gi;

  /* initialize the generic integrator (same integrator already used
     in background, thermodynamics and perturbation modules) */
  pipaw.ppm = ppm;
  pipaw.N = ppm->in_size;
  pipaw.k = k;

  class_call(initialize_generic_integrator(pipaw.N,&gi),
             gi.error_message,
             ppm->error_message);

  /* initial conditions for the perturbations, Bunch-Davies vacuum */
  y[ppm->index_in_ksi_re]=1./sqrt(2.*k);
  y[ppm->index_in_ksi_im]=0.;
  y[ppm->index_in_dksi_re]=0.;
  y[ppm->index_in_dksi_im]=-k*y[ppm->index_in_ksi_re];

  y[ppm->index_in_ah_re]=1./sqrt(2.*k);
  y[ppm->index_in_ah_im]=0.;
  y[ppm->index_in_dah_re]=0.;
  y[ppm->index_in_dah_im]=-k*y[ppm->index_in_ah_re];

  /* intialize variable used for deciding when to stop the calculation (= when the curvature remains stable) */
  curvature_new=1.e10;

  /* intialize conformal time to arbitrary value (here, only variations
     of tau matter: the equations that we integrate do not depend
     explicitely on time) */
  tau_end = 0;

  /* compute derivative of initial vector and infer first value of adaptative time-step */
  class_call(primordial_inflation_derivs(tau_end,
                                         y,
                                         dy,
                                         &pipaw,
                                         ppm->error_message),
             ppm->error_message,
             ppm->error_message);

  dtau = ppr->primordial_inflation_pt_stepsize*2.*_PI_/MAX(sqrt(fabs(dy[ppm->index_in_dksi_re]/y[ppm->index_in_ksi_re])),k);

  /* loop over time */
  do {

    /* new time interval [tau_start, tau_end] over which equations will be integrated */
    tau_start = tau_end;

    tau_end = tau_start + dtau;

    class_test(dtau/tau_start < ppr->smallest_allowed_variation,
               ppm->error_message,
               "integration step: relative change in time =%e < machine precision : leads either to numerical error or infinite loop",dtau/tau_start);

    /* evolve the system */
    class_call(generic_integrator(primordial_inflation_derivs,
                                  tau_start,
                                  tau_end,
                                  y,
                                  &pipaw,
                                  ppr->primordial_inflation_tol_integration,
                                  ppr->smallest_allowed_variation,
                                  &gi),
               gi.error_message,
               ppm->error_message);

    /* compute derivatives at tau_end, useful to infer new time step and spectra */
    class_call(primordial_inflation_derivs(tau_end,
                                           y,
                                           dy,
                                           &pipaw,
                                           ppm->error_message),
               ppm->error_message,
               ppm->error_message);

    /* new time step */
    dtau = ppr->primordial_inflation_pt_stepsize*2.*_PI_/MAX(sqrt(fabs(dy[ppm->index_in_dksi_re]/y[ppm->index_in_ksi_re])),k);

    /* new aH */
    aH = dy[ppm->index_in_a]/y[ppm->index_in_a];

    /* store previous value of curvature (at tau_start) */
    curvature_old =  curvature_new;

    /* new curvature */
    z = y[ppm->index_in_a]*y[ppm->index_in_dphi]/aH;
    ksi2 = y[ppm->index_in_ksi_re]*y[ppm->index_in_ksi_re]+y[ppm->index_in_ksi_im]*y[ppm->index_in_ksi_im];
    curvature_new = k*k*k/2./_PI_/_PI_*ksi2/z/z;

    /* variation of curvature with time (dimensionless) */
    dlnPdN = (curvature_new-curvature_old)/dtau*y[ppm->index_in_a]/dy[ppm->index_in_a]/curvature_new;

    /* stop when (k >> aH) AND curvature is stable */
  } while ((k/aH >= ppr->primordial_inflation_ratio_max) || (fabs(dlnPdN) > ppr->primordial_inflation_tol_curvature));

  /* clean the generic integrator */
  class_call(cleanup_generic_integrator(&gi),
             gi.error_message,
             ppm->error_message);

  /* store final value of curvature for this wavenumber */
  *curvature = curvature_new;

  /* stor final value of tensor perturbation for this wavenumber */
  ah2 = y[ppm->index_in_ah_re]*y[ppm->index_in_ah_re]+y[ppm->index_in_ah_im]*y[ppm->index_in_ah_im];
  *tensor = 32.*k*k*k/_PI_*ah2/y[ppm->index_in_a]/y[ppm->index_in_a];

  //fprintf(stdout,"%g %g %g %g %g\n",k,*curvature,*tensor,*tensor/(*curvature),dlnPdN);

  return _SUCCESS_;
}

/**
 * Routine searching for the inflationary attractor solution, by
 * iterations, with a given tolerance. If no solution found within
 * tolerance, returns error message.
 *
 * @param ppm       Input: pointer to primordial structure
 * @param ppr       Input: pointer to precision structure
 * @param phi_0     Input: field value at which we wish to find the solution
 * @param precision Input: tolerance on output values (if too large, an attractor will always considered to be found)
 * @param y         Input: running vector of background variables, already allocated and initialized
 * @param dy        Input: running vector of background derivatives, already allocated
 * @param H0        Output: Hubble value at phi_0 for attractor solution
 * @param dphidt_0  Output: field derivative value at phi_0 for attractor solution
 * @return the error status
 */

int primordial_inflation_find_attractor(
                                        struct primordial * ppm,
                                        struct precision * ppr,
                                        double phi_0,
                                        double precision,
                                        double * y,
                                        double * dy,
                                        double * H_0,
                                        double * dphidt_0
                                        ) {

  double V_0,dV_0,ddV_0;
  double V=0.,dV=0.,ddV=0.;
  double a;
  double dphidt,dphidt_0new,dphidt_0old,phi;
  int counter;

  class_call(primordial_inflation_potential(ppm,phi_0,&V_0,&dV_0,&ddV_0),
             ppm->error_message,
             ppm->error_message);

  dphidt_0new = -dV_0/3./sqrt((8.*_PI_/3.)*V_0);
  phi = phi_0;
  counter = 0;

  dphidt_0old = dphidt_0new/(precision+2.);

  while (fabs(dphidt_0new/dphidt_0old-1.) >= precision) {

    counter ++;
    class_test(counter >= ppr->primordial_inflation_attractor_maxit,
               ppm->error_message,
               "could not converge after %d iterations: there exists no attractor solution near phi=%g. Potential probably too steep in this region, or precision parameter primordial_inflation_attractor_precision=%g too small",
               counter,
               phi_0,
               precision);

    dphidt_0old = dphidt_0new;

    phi=phi+dV_0/V_0/16./_PI_;

    class_call(primordial_inflation_check_potential(ppm,phi),
               ppm->error_message,
               ppm->error_message);

    class_call(primordial_inflation_potential(ppm,phi,&V,&dV,&ddV),
               ppm->error_message,
               ppm->error_message);

    a = 1.;
    dphidt = -dV/3./sqrt((8.*_PI_/3.)*V);
    y[ppm->index_in_a]=a;
    y[ppm->index_in_phi]=phi;
    y[ppm->index_in_dphi]=a*dphidt;

    class_call(primordial_inflation_evolve_background(ppm,ppr,y,dy,phi_0),
               ppm->error_message,
               ppm->error_message);

    dphidt_0new = y[ppm->index_in_dphi]/y[ppm->index_in_a];

  }

  *dphidt_0 = dphidt_0new;
  *H_0 = sqrt((8.*_PI_/3.)*(0.5*dphidt_0new*dphidt_0new+V_0));

  // fprintf(stderr,"attractor found with %g %g\n",*dphidt_0,*H_0);

  return _SUCCESS_;
}

/**
 * Routine integrating background equations only, from initial field
 * value (stored in y) to final value (passed as phi_stop). In output,
 * y contrains the final background values.
 *
 * @param ppm       Input: pointer to primordial structure
 * @param ppr       Input: pointer to precision structure
 * @param y         Input/output: running vector of background variables, already allocated and initialized
 * @param dy        Input: running vector of background derivatives, already allocated
 * @param phi_stop  Input: final field value
 * @return the error status
 */

int primordial_inflation_evolve_background(
                                           struct primordial * ppm,
                                           struct precision * ppr,
                                           double * y,
                                           double * dy,
                                           double phi_stop) {

  double tau_start,tau_end,dtau;
  double aH;
  double epsilon,epsilon_old;
  struct primordial_inflation_parameters_and_workspace pipaw;
  struct generic_integrator_workspace gi;

  pipaw.ppm = ppm;
  pipaw.N = ppm->in_bg_size;

  class_call(initialize_generic_integrator(pipaw.N,&gi),
             gi.error_message,
             ppm->error_message);

  class_call(primordial_inflation_get_epsilon(ppm,
                                              y[ppm->index_in_phi],
                                              &epsilon),
             ppm->error_message,
             ppm->error_message);

  tau_end = 0;

  class_call(primordial_inflation_derivs(tau_end,
                                         y,
                                         dy,
                                         &pipaw,
                                         ppm->error_message),
             ppm->error_message,
             ppm->error_message);

  aH = dy[ppm->index_in_a]/y[ppm->index_in_a];
  dtau = ppr->primordial_inflation_bg_stepsize*MIN(1./aH,fabs(y[ppm->index_in_dphi]/dy[ppm->index_in_dphi]));

  while (y[ppm->index_in_phi] <= (phi_stop-y[ppm->index_in_dphi]*dtau)) {

    class_call_except(primordial_inflation_check_potential(ppm,
                                                           y[ppm->index_in_phi]),
                      ppm->error_message,
                      ppm->error_message,
                      cleanup_generic_integrator(&gi));

    tau_start = tau_end;

    class_call_except(primordial_inflation_derivs(tau_start,
                                                  y,
                                                  dy,
                                                  &pipaw,
                                                  ppm->error_message),
                      ppm->error_message,
                      ppm->error_message,
                      cleanup_generic_integrator(&gi));

    aH = dy[ppm->index_in_a]/y[ppm->index_in_a];
    dtau = ppr->primordial_inflation_bg_stepsize*MIN(1./aH,fabs(y[ppm->index_in_dphi]/dy[ppm->index_in_dphi]));

    tau_end = tau_start + dtau;

    class_test_except(dtau/tau_start < ppr->smallest_allowed_variation,
                      ppm->error_message,
                      cleanup_generic_integrator(&gi),
                      "integration step: relative change in time =%e < machine precision : leads either to numerical error or infinite loop",
                      dtau/tau_start);

    class_call_except(generic_integrator(primordial_inflation_derivs,
                                         tau_start,
                                         tau_end,
                                         y,
                                         &pipaw,
                                         ppr->primordial_inflation_tol_integration,
                                         ppr->smallest_allowed_variation,
                                         &gi),
                      gi.error_message,
                      ppm->error_message,
                      cleanup_generic_integrator(&gi));

    epsilon_old = epsilon;

    class_call_except(primordial_inflation_get_epsilon(ppm,
                                                       y[ppm->index_in_phi],
                                                       &epsilon),
                      ppm->error_message,
                      ppm->error_message,
                      cleanup_generic_integrator(&gi));

    class_test_except((epsilon > 1) && (epsilon_old <= 1),
                      ppm->error_message,
                      cleanup_generic_integrator(&gi),
                      "Inflaton evolution crosses the border from epsilon<1 to epsilon>1 at phi=%g. Inflation disrupted during the observable e-folds",
                      y[ppm->index_in_phi]);



  }

  class_call(cleanup_generic_integrator(&gi),
             gi.error_message,
             ppm->error_message);

  class_call(primordial_inflation_derivs(tau_end,
                                         y,
                                         dy,
                                         &pipaw,
                                         ppm->error_message),
             ppm->error_message,
             ppm->error_message);

  dtau = (phi_stop-y[ppm->index_in_phi])/dy[ppm->index_in_dphi];
  y[ppm->index_in_a] += dy[ppm->index_in_a]*dtau;
  y[ppm->index_in_phi] += dy[ppm->index_in_phi]*dtau;
  y[ppm->index_in_dphi] += dy[ppm->index_in_dphi]*dtau;

  return _SUCCESS_;
}

/**
 * Routine integrating background equations only, from initial AH
 * value (stored in y) to final value (passed as aH_stop). In output,
 * y contrains the final background values.
 *
 * @param ppm       Input: pointer to primordial structure
 * @param ppr       Input: pointer to precision structure
 * @param y         Input/output: running vector of background variables, already allocated and initialized
 * @param dy        Input: running vector of background derivatives, already allocated
 * @param aH_stop  Input: final aH value
 * @return the error status
 */

int primordial_inflation_reach_aH(
                                  struct primordial * ppm,
                                  struct precision * ppr,
                                  double * y,
                                  double * dy,
                                  double aH_stop
                                  ) {

  double tau_start,tau_end,dtau;
  double aH;
  struct primordial_inflation_parameters_and_workspace pipaw;
  struct generic_integrator_workspace gi;

  pipaw.ppm = ppm;
  pipaw.N = ppm->in_bg_size;

  class_call(initialize_generic_integrator(pipaw.N,&gi),
             gi.error_message,
             ppm->error_message);

  tau_end = 0;

  class_call(primordial_inflation_derivs(tau_end,
                                         y,
                                         dy,
                                         &pipaw,
                                         ppm->error_message),
             ppm->error_message,
             ppm->error_message);

  while (dy[ppm->index_in_a]/y[ppm->index_in_a] < aH_stop) {

    class_call(primordial_inflation_check_potential(ppm,
                                                    y[ppm->index_in_phi]),
               ppm->error_message,
               ppm->error_message);

    tau_start = tau_end;

    class_call(primordial_inflation_derivs(tau_start,
                                           y,
                                           dy,
                                           &pipaw,
                                           ppm->error_message),
               ppm->error_message,
               ppm->error_message);

    aH = dy[ppm->index_in_a]/y[ppm->index_in_a];
    dtau = ppr->primordial_inflation_bg_stepsize*MIN(1./aH,fabs(y[ppm->index_in_dphi]/dy[ppm->index_in_dphi]));

    tau_end = tau_start + dtau;

    class_test(dtau/tau_start < ppr->smallest_allowed_variation,
               ppm->error_message,
               "integration step: relative change in time =%e < machine precision : leads either to numerical error or infinite loop",dtau/tau_start);

    class_call(generic_integrator(primordial_inflation_derivs,
                                  tau_start,
                                  tau_end,
                                  y,
                                  &pipaw,
                                  ppr->primordial_inflation_tol_integration,
                                  ppr->smallest_allowed_variation,
                                  &gi),
               gi.error_message,
               ppm->error_message);

  }

  class_call(cleanup_generic_integrator(&gi),
             gi.error_message,
             ppm->error_message);

  return _SUCCESS_;
}

/**
 * Routine checking positivity and negative slope of potential. The
 * negative slope is an arbitrary choice. Currently the code can only
 * deal with monotonic variations of the inflaton during inflation. So
 * the slope had to be always negative or always positive... we took
 * the first option.
 *
 * @param ppm       Input: pointer to primordial structure
 * @param phi       Input: field value where to perform the check
 * @return the error status
 */

int primordial_inflation_check_potential(
                                         struct primordial * ppm,
                                         double phi
                                         ) {

  double V=0.,dV=0.,ddV=0.;

  class_call(primordial_inflation_potential(ppm,phi,&V,&dV,&ddV),
             ppm->error_message,
             ppm->error_message);

  class_test(V <= 0.,
             ppm->error_message,
             "This potential becomes negative at phi=%g, before the end of observable inflation. It  cannot be treated by this code",
             phi);

  class_test(dV >= 0.,
             ppm->error_message,
             "All the code is written for the case dV/dphi<0. Here, in phi=%g, we have dV/dphi=%g. This potential cannot be treated by this code",
             phi,dV);

  return _SUCCESS_;
}

/**
 * Routine computing the first slow-roll parameter epsilon
 *
 * @param ppm       Input: pointer to primordial structure
 * @param phi       Input: field value where to compute epsilon
 * @param epsilon   Ouput: result
 * @return the error status
 */

int primordial_inflation_get_epsilon(
                                     struct primordial * ppm,
                                     double phi,
                                     double * epsilon
                                     ) {

  double V=0.,dV=0.,ddV=0.;

  class_call(primordial_inflation_potential(ppm,phi,&V,&dV,&ddV),
             ppm->error_message,
             ppm->error_message);

  *epsilon = 1./16./_PI_*pow(dV/V,2);

  //*eta = 1./8./pi*(ddV/V)

  return _SUCCESS_;

}

/**
 * Routine creturning derivative of system of background/perturbation
 * variables. Like other routines used by the generic integrator
 * (background_derivs, thermodynamics_derivs, perturb_derivs), this
 * routine has a generci list of arguments, and a slightly different error management, with the error
 * message returned directly in an ErrMsg field.
 *
 * @param tau                      Input: time (not used explicitely inside the routine, but requested by the generic integrator)
 * @param y                        Input/output: running vector of background variables, already allocated and initialized
 * @param dy                       Input: running vector of background derivatives, already allocated
 * @param parameters_and_workspace Input: all necessary input variables apart from y
 * @return the error status
 */

int primordial_inflation_derivs(
                                double tau,
                                double * y,
                                double * dy,
                                void * parameters_and_workspace,
                                ErrorMsg error_message
                                ) {

  struct primordial_inflation_parameters_and_workspace * ppipaw;
  struct primordial * ppm;

  ppipaw = parameters_and_workspace;
  ppm = ppipaw->ppm;

  class_call(primordial_inflation_potential(ppm,
                                            y[ppm->index_in_phi],
                                            &(ppipaw->V),
                                            &(ppipaw->dV),
                                            &(ppipaw->ddV)),
             ppm->error_message,
             ppm->error_message);

  // BACKGROUND

  // a**2 V
  ppipaw->a2V=y[ppm->index_in_a]*y[ppm->index_in_a]*ppipaw->V;
  // a**2 dV/dphi
  ppipaw->a2dV=y[ppm->index_in_a]*y[ppm->index_in_a]*ppipaw->dV;
  // a H = a'/a
  ppipaw->aH = sqrt((8*_PI_/3.)*(0.5*y[ppm->index_in_dphi]*y[ppm->index_in_dphi]+ppipaw->a2V));

  // 1: a
  dy[ppm->index_in_a]=y[ppm->index_in_a]*ppipaw->aH;
  // 2: phi
  dy[ppm->index_in_phi]=y[ppm->index_in_dphi];
  // 3: dphi/dtau
  dy[ppm->index_in_dphi]=-2.*ppipaw->aH*y[ppm->index_in_dphi]-ppipaw->a2dV;

  if (ppipaw->N == ppm->in_bg_size)
    return _SUCCESS_;

  // PERTURBATIONS

  // a**2 d2V/dphi2
  ppipaw->a2ddV=y[ppm->index_in_a]*y[ppm->index_in_a]*ppipaw->ddV;
  // z''/z
  ppipaw->zpp_over_z=2*ppipaw->aH*ppipaw->aH - ppipaw->a2ddV - 4.*_PI_*(7.*y[ppm->index_in_dphi]*y[ppm->index_in_dphi]+4.*y[ppm->index_in_dphi]/ppipaw->aH*ppipaw->a2dV)
    +32.*_PI_*_PI_*pow(y[ppm->index_in_dphi],4)/pow(ppipaw->aH,2);
  // a''/a
  ppipaw->app_over_a=2.*ppipaw->aH*ppipaw->aH - 4.*_PI_*y[ppm->index_in_dphi]*y[ppm->index_in_dphi];

  // SCALARS
  // 4: ksi_re
  dy[ppm->index_in_ksi_re]=y[ppm->index_in_dksi_re];
  // 5: ksi_im
  dy[ppm->index_in_ksi_im]=y[ppm->index_in_dksi_im];
  // 6: d ksi_re / dtau
  dy[ppm->index_in_dksi_re]=-(ppipaw->k*ppipaw->k-ppipaw->zpp_over_z)*y[ppm->index_in_ksi_re];
  // 7: d ksi_im / dtau
  dy[ppm->index_in_dksi_im]=-(ppipaw->k*ppipaw->k-ppipaw->zpp_over_z)*y[ppm->index_in_ksi_im];

  // TENSORS
  // 8: ah_re
  dy[ppm->index_in_ah_re]=y[ppm->index_in_dah_re];
  // 9: ah_im
  dy[ppm->index_in_ah_im]=y[ppm->index_in_dah_im];
  // 10: d ah_re / dtau
  dy[ppm->index_in_dah_re]=-(ppipaw->k*ppipaw->k-ppipaw->app_over_a)*y[ppm->index_in_ah_re];
  // 11: d ah_im / dtau
  dy[ppm->index_in_dah_im]=-(ppipaw->k*ppipaw->k-ppipaw->app_over_a)*y[ppm->index_in_ah_im];

  return _SUCCESS_;

}

/**
 * This routine reads the primordial spectrum from an external command,
 * and stores the tabulated values.
 * The sampling of the k's given by the external command is preserved.
 *
 * Author: Jesus Torrado (torradocacho@lorentz.leidenuniv.nl)
 * Date:   2013-12-20
 *
 * @param ppm  Input/output: pointer to perturbs structure
 * @param ppm  Input/output: pointer to primordial structure
 * @return the error status
 */

int primordial_external_spectrum_init(
                                      struct perturbs * ppt,
                                      struct primordial * ppm
                                      ) {

  char arguments[_ARGUMENT_LENGTH_MAX_];
  char line[_LINE_LENGTH_MAX_];
  char command_with_arguments[2*_ARGUMENT_LENGTH_MAX_];
  FILE *process;
  int n_data_guess, n_data = 0;
  double *k = NULL, *pks = NULL, *pkt = NULL, *tmp = NULL;
  double this_k, this_pks, this_pkt;
  int status;
  int index_k;

  /** 1. Initialization */
  /* Prepare the data (with some initial size) */
  n_data_guess = 100;
  k   = (double *)malloc(n_data_guess*sizeof(double));
  pks = (double *)malloc(n_data_guess*sizeof(double));
  if (ppt->has_tensors == _TRUE_)
    pkt = (double *)malloc(n_data_guess*sizeof(double));
  /* Prepare the command */
  /* If the command is just a "cat", no arguments need to be passed */
  if(strncmp("cat ", ppm->command, 4) == 0) {
    sprintf(arguments, " ");
  }
  /* otherwise pass the list of arguments */
  else {
    sprintf(arguments, " %g %g %g %g %g %g %g %g %g %g",
            ppm->custom1, ppm->custom2, ppm->custom3, ppm->custom4, ppm->custom5,
            ppm->custom6, ppm->custom7, ppm->custom8, ppm->custom9, ppm->custom10);
  }
  /* write the actual command in a string */
  sprintf(command_with_arguments, "%s %s", ppm->command, arguments);
  if (ppm->primordial_verbose > 0)
    printf(" -> running: %s\n",command_with_arguments);

  /** 2. Launch the command and retrieve the output */
  /* Launch the process */
  process = popen(command_with_arguments, "r");
  class_test(process == NULL,
             ppm->error_message,
             "The program failed to set the environment for the external command. Maybe you ran out of memory.");
  /* Read output and store it */
  while (fgets(line, sizeof(line)-1, process) != NULL) {
    if (ppt->has_tensors == _TRUE_) {
      sscanf(line, "%lf %lf %lf", &this_k, &this_pks, &this_pkt);
    }
    else {
      sscanf(line, "%lf %lf", &this_k, &this_pks);
    }
    /* Standard technique in C: if too many data, double the size of the vectors */
    /* (it is faster and safer that reallocating every new line) */
    if((n_data+1) > n_data_guess) {
      n_data_guess *= 2;
      tmp = (double *)realloc(k,   n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 ppm->error_message,
                 "Error allocating memory to read the external spectrum.\n");
      k = tmp;
      tmp = (double *)realloc(pks, n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 ppm->error_message,
                 "Error allocating memory to read the external spectrum.\n");
      pks = tmp;
      if (ppt->has_tensors == _TRUE_) {
        tmp = (double *)realloc(pkt, n_data_guess*sizeof(double));
        class_test(tmp == NULL,
                   ppm->error_message,
                   "Error allocating memory to read the external spectrum.\n");
        pkt = tmp;
      };
    };
    /* Store */
    k  [n_data]   = this_k;
    pks[n_data]   = this_pks;
    if (ppt->has_tensors == _TRUE_) {
      pkt[n_data] = this_pkt;
    }
    n_data++;
    /* Check ascending order of the k's */
    if(n_data>1) {
      class_test(k[n_data-1] <= k[n_data-2],
                 ppm->error_message,
                 "The k's are not strictly sorted in ascending order, "
                 "as it is required for the calculation of the splines.\n");
    }
  }
  /* Close the process */
  status = pclose(process);
  class_test(status != 0.,
             ppm->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");
  /* Test limits of the k's */
  class_test(k[1] > ppt->k[0],
             ppm->error_message,
             "Your table for the primordial spectrum does not have "
             "at least 2 points before the minimum value of k: %e . "
             "The splines interpolation would not be safe.",ppt->k[0]);
  class_test(k[n_data-2] < ppt->k[ppt->k_size-1],
             ppm->error_message,
             "Your table for the primordial spectrum does not have "
             "at least 2 points after the maximum value of k: %e . "
             "The splines interpolation would not be safe.",ppt->k[ppt->k_size-1]);

  /** 3. Store the read results into CLASS structures */
  ppm->lnk_size = n_data;
  /** Make room */
  class_realloc(ppm->lnk,
                ppm->lnk,
                ppm->lnk_size*sizeof(double),
                ppm->error_message);
  class_realloc(ppm->lnpk[ppt->index_md_scalars],
                ppm->lnpk[ppt->index_md_scalars],
                ppm->lnk_size*sizeof(double),
                ppm->error_message);
  class_realloc(ppm->ddlnpk[ppt->index_md_scalars],
                ppm->ddlnpk[ppt->index_md_scalars],
                ppm->lnk_size*sizeof(double),
                ppm->error_message);
  if (ppt->has_tensors == _TRUE_) {
    class_realloc(ppm->lnpk[ppt->index_md_tensors],
                  ppm->lnpk[ppt->index_md_tensors],
                  ppm->lnk_size*sizeof(double),
                  ppm->error_message);
    class_realloc(ppm->ddlnpk[ppt->index_md_tensors],
                  ppm->ddlnpk[ppt->index_md_tensors],
                  ppm->lnk_size*sizeof(double),
                  ppm->error_message);
  };
  /** Store them */
  for (index_k=0; index_k<ppm->lnk_size; index_k++) {
    ppm->lnk[index_k] = log(k[index_k]);
    ppm->lnpk[ppt->index_md_scalars][index_k] = log(pks[index_k]);
    if (ppt->has_tensors == _TRUE_)
      ppm->lnpk[ppt->index_md_tensors][index_k] = log(pkt[index_k]);
    /* DEBUG (with tensors)
       fprintf(stderr,"Storing[%d(+1) of %d]: \n k = %g == %g\n pks = %g == %g\n pkt = %g == %g\n",
       index_k, n_data,
       ppm->lnk[index_k], log(k[index_k]),
       ppm->lnpk[ppt->index_md_scalars][index_k], log(pks[index_k]),
       ppm->lnpk[ppt->index_md_tensors][index_k], log(pkt[index_k]));
    */
  };
  /** Release the memory used locally */
  free(k);
  free(pks);
  if (ppt->has_tensors == _TRUE_)
    free(pkt);
  /** Tell CLASS that the are scalar (and tensor) modes */
  ppm->is_non_zero[ppt->index_md_scalars][ppt->index_ic_ad] = _TRUE_;
  if (ppt->has_tensors == _TRUE_)
    ppm->is_non_zero[ppt->index_md_tensors][ppt->index_ic_ten] = _TRUE_;

  return _SUCCESS_;
}
