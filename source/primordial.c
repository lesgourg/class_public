/** @file primordial.c Documented primordial module.
 *
 * Julien Lesgourgues, 24.08.2010
 *
 * This module computes the primordial spectra. It can be used in different modes:
 * simple parametric form, evolving inflaton perturbations, etc. So far only
 * the mode corresponding to a simple analytic form in terms of amplitudes, tilts
 * and runnings has been developed.
 *
 * The following functions can be called from other modules:
 *
 * -# primordial_init() at the beginning (anytime after perturbations_init() and before harmonic_init())
 * -# primordial_spectrum_at_k() at any time for computing P(k) at any k
 * -# primordial_free() at the end
 */

#include "primordial.h"
#include "parallel.h"

/**
 * Primordial spectra for arbitrary argument and for all initial conditions.
 *
 * This routine evaluates the primordial spectrum at a given value of k by
 * interpolating in the pre-computed table.
 *
 * When k is not in the pre-computed range but the spectrum can be found
 * analytically, it finds it. Otherwise returns an error.
 *
 * Can be called in two modes; linear or logarithmic:
 *
 * - linear: takes k, returns P(k)
 *
 * - logarithmic: takes ln(k), return ln(P(k))
 *
 * One little subtlety: in case of several correlated initial conditions,
 * the cross-correlation spectrum can be negative. Then, in logarithmic mode,
 * the non-diagonal elements contain the cross-correlation angle \f$ P_{12}/\sqrt{P_{11} P_{22}}\f$
 * (from -1 to 1) instead of \f$\ln{P_{12}}\f$
 *
 * This function can be
 * called from whatever module at whatever time, provided that
 * primordial_init() has been called before, and primordial_free() has not
 * been called yet.
 *
 * @param ppm        Input: pointer to primordial structure containing tabulated primordial spectrum
 * @param index_md   Input: index of mode (scalar, tensor, ...)
 * @param mode       Input: linear or logarithmic
 * @param input      Input: wavenumber in 1/Mpc (linear mode) or its logarithm (logarithmic mode)
 * @param output     Output: for each pair of initial conditions, primordial spectra P(k) in \f$Mpc^3\f$ (linear mode), or their logarithms and cross-correlation angles (logarithmic mode)
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

  /** - otherwise, interpolate in the pre-computed table */

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
 * This routine initializes the primordial structure (in particular, it computes table of primordial spectrum values)
 *
 * @param ppr Input: pointer to precision structure (defines method and precision for all computations)
 * @param ppt Input: pointer to perturbation structure (useful for knowing k_min, k_max, etc.)
 * @param ppm Output: pointer to initialized primordial structure
 * @return the error status
 */

int primordial_init(
                    struct precision  * ppr,
                    struct perturbations   * ppt,
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

  k_min = ppt->k_min; /* first value, inferred from perturbations structure */
  if (ppm->has_k_max_for_primordial_pk == _TRUE_){
    k_max = ppm->k_max_for_primordial_pk; /* last value, user-defined (i.e. if specified in .ini file) */
  }
  else{
    k_max = ppt->k_max; /* last value, inferred from perturbations structure */
  }

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

  /** - allocate and fill values of \f$ \ln{k}\f$'s */

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

  /** - deal with case of analytic primordial spectra (with amplitudes, tilts, runnings, etc.) */

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

  /** - deal with case of inflation with given \f$V(\phi)\f$ or \f$H(\phi)\f$ */

  else if ((ppm->primordial_spec_type == inflation_V) || (ppm->primordial_spec_type == inflation_H) || (ppm->primordial_spec_type == inflation_V_end)) {

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

  /** - deal with the case of external calculation of \f$ P_k \f$*/

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
               "primordial spectrum type not recognized");

  }

  /** - compute second derivative of each \f$ \ln{P_k} \f$ versus lnk with spline, in view of interpolation */

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

  /** - derive spectral parameters from numerically computed spectra
      (not used by the rest of the code, but useful to keep in memory for several types of investigation) */

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

      /** - expression for alpha_s comes from:

          `ns_2 = (lnpk_plus-lnpk_pivot)/(dlnk)+1`

          `ns_1 = (lnpk_pivot-lnpk_minus)/(dlnk)+1`

          `alpha_s = dns/dlnk = (ns_2-ns_1)/dlnk = (lnpk_plus-lnpk_pivot-lnpk_pivot+lnpk_minus)/(dlnk)/(dlnk)`

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

      /** - expression for beta_s:

          `ppm->beta_s = (alpha_plus-alpha_minus)/dlnk = (lnpk_plusplus-2.*lnpk_plus+lnpk_pivot -
          (lnpk_pivot-2.*lnpk_minus+lnpk_minusminus)/pow(dlnk,3)`
      **/

      /* Simplification of the beta_s expression: */

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
        printf(" -> r=%g  n_t=%g  alpha_t=%g\n",ppm->r,ppm->n_t,ppm->alpha_t);

    }

  }

  ppm->is_allocated = _TRUE_;

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

  ppm->is_allocated = _FALSE_;

  return _SUCCESS_;
}

/**
 * This routine defines indices and allocates tables in the primordial structure
 *
 * @param ppt  Input: pointer to perturbation structure
 * @param ppm  Input/output: pointer to primordial structure
 * @return the error status
 */

int primordial_indices(
                       struct perturbations   * ppt,
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
 * @param kmin Input: first value
 * @param kmax Input: last value that we should encompass
 * @param k_per_decade Input: number of k per decade
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
 * @param ppt  Input: pointer to perturbation structure
 * @param ppm  Input/output: pointer to primordial structure
 * @return the error status
 */

int primordial_analytic_spectrum_init(
                                      struct perturbations   * ppt,
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
 * @param V              Output: inflaton potential in units of \f$ Mp^4\f$
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

  double e,de,dde,mu,dmu,ddmu,l,dl,ddl,p,dp,ddp;

  switch (ppm->potential) {

    /* V(phi)=polynomial in phi */
  case polynomial:

    *V   = ppm->V0+phi*ppm->V1+pow(phi,2)/2.*ppm->V2+pow(phi,3)/6.*ppm->V3+pow(phi,4)/24.*ppm->V4;
    *dV  = ppm->V1+phi*ppm->V2+pow(phi,2)/2.*ppm->V3+pow(phi,3)/6.*ppm->V4;
    *ddV = ppm->V2+phi*ppm->V3+pow(phi,2)/2.*ppm->V4;
    break;

    /* V(phi) = Lambda^4(1+cos(phi/f)) = V0 (1+cos(phi/V1)) */
  case natural:

    *V   = ppm->V0*(1.+cos(phi/ppm->V1));
    *dV  = -ppm->V0/ppm->V1*sin(phi/ppm->V1);
    *ddV = -ppm->V0/ppm->V1/ppm->V1*cos(phi/ppm->V1);
    break;

    /* Higgs inflation from arXiv:1403.6078 */
  case higgs_inflation:

    // correspondence with 1403.6078:
    // V0 = b
    // V1 = ksi
    // V2 = kappa
    // V3 = delta_lambda
    // mu = bar(mu)/M_P
    // phi = -chi/M_P

    e = exp(2./sqrt(6.)*sqrt(8.*_PI_)*phi);
    de = 2./sqrt(6.)*sqrt(8.*_PI_)*e;
    dde = 2./3. * 8.*_PI_ * e;

    mu = pow(1.-e,0.5);
    dmu = -0.5*de*pow(1.-e,-0.5);
    ddmu = -0.5*dde*pow(1.-e,-0.5)-0.25*de*de*pow(1.-e,-1.5);

    l = log(mu/ppm->V2);
    dl = dmu/mu;
    ddl = ddmu/mu - dl*dl;

    p = 1./16. + ppm->V3/ppm->V0 + l*l;
    dp = 2.*dl*l;
    ddp = 2.*ddl*l+2.*dl*dl;

    *V = ppm->V0/4./pow(8.*_PI_,2)/ppm->V1/ppm->V1*p*pow(mu,4);
    *dV = ppm->V0/4./pow(8.*_PI_,2)/ppm->V1/ppm->V1*(dp*pow(mu,4)+4.*p*dmu*pow(mu,3));
    *ddV = ppm->V0/4./pow(8.*_PI_,2)/ppm->V1/ppm->V1*(ddp*pow(mu,4)+8.*dp*dmu*pow(mu,3)+4.*p*ddmu*pow(mu,3)+12.*p*pow(dmu*mu,2));

    //fprintf(stderr,"%e  %e  %e\n",*V,p,mu);

    break;

    /* code here other shapes */

  default:
    class_stop(ppm->error_message,"ppm->potential=%d different from all known cases",ppm->potential);
    break;
  }

  return _SUCCESS_;
}

/**
 * This routine encodes the function \f$ H(\phi)\f$
 *
 * @param ppm            Input: pointer to primordial structure
 * @param phi            Input: background inflaton field value in units of Mp
 * @param H              Output: Hubble parameters in units of Mp
 * @param dH             Output: \f$ dH / d\phi \f$
 * @param ddH            Output: \f$ d^2H / d\phi^2 \f$
 * @param dddH           Output: \f$ d^3H / d\phi^3 \f$
 * @return the error status
 */

int primordial_inflation_hubble(
                                struct primordial * ppm,
                                double phi,
                                double * H,
                                double * dH,
                                double * ddH,
                                double * dddH
                                ) {

  *H =    ppm->H0 + phi*ppm->H1 + pow(phi,2)/2.*ppm->H2 + pow(phi,3)/6.*ppm->H3 + pow(phi,4)/24.*ppm->H4;
  *dH =   ppm->H1 + phi*ppm->H2 + pow(phi,2)/2.*ppm->H3 + pow(phi,3)/6.*ppm->H4;
  *ddH =  ppm->H2 + phi*ppm->H3 + pow(phi,2)/2.*ppm->H4;
  *dddH = ppm->H3 + phi*ppm->H4;

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

  /* indices for background quantities */
  ppm->index_in_a = index_in;
  index_in ++;
  ppm->index_in_phi = index_in;
  index_in ++;
  if ((ppm->primordial_spec_type == inflation_V) || (ppm->primordial_spec_type == inflation_V_end)) {
    ppm->index_in_dphi = index_in;
    index_in ++;
  }

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
                                         struct perturbations * ppt,
                                         struct primordial * ppm,
                                         struct precision *ppr
                                         ) {
  /** Summary: */
  /** - define local variables */
  double * y;
  double * y_ini;
  double * dy;
  double a_pivot;
  double a_try;
  double H_pivot;
  double H_try;
  double phi_try;
  double dphidt_pivot;
  double dphidt_try;
  double aH_ini,aH_end;
  double k_max,k_min;
  int counter;
  double dH,ddH,dddH;

  /** - allocate vectors for background/perturbed quantities */
  class_alloc(y,ppm->in_size*sizeof(double),ppm->error_message);
  class_alloc(y_ini,ppm->in_size*sizeof(double),ppm->error_message);
  class_alloc(dy,ppm->in_size*sizeof(double),ppm->error_message);

  /** - eventually, needs first to find phi_pivot */
  if (ppm->primordial_spec_type == inflation_V_end) {

    class_call(primordial_inflation_find_phi_pivot(ppm,ppr,y,dy),
               ppm->error_message,
               ppm->error_message);

  }
  else {
    ppm->phi_pivot = 0.;
  }

  // uncomment these lines if for checking, you want first-order slow-roll predictions
  /*
    if (ppm->primordial_verbose>0) {
    if ((ppm->primordial_spec_type == inflation_V) || (ppm->primordial_spec_type == inflation_V_end)) {
    double V,dV,ddV;
    class_call(primordial_inflation_check_potential(ppm,ppm->phi_pivot,&V,&dV,&ddV),
    ppm->error_message,
    ppm->error_message);
    fprintf(stdout," -> 1st-order slow-roll prediction for A_s: %g\n",128.*_PI_/3.*pow(V,3)/pow(dV,2));
    fprintf(stdout," -> 1st-order slow-roll prediction for T/S: %g\n",pow(dV/V,2)/_PI_);
    fprintf(stdout," -> 1st-order slow-roll prediction for A_T: %g\n",pow(dV/V,2)/_PI_*128.*_PI_/3.*pow(V,3)/pow(dV,2));
    fprintf(stdout," -> 1st-order slow-roll prediction for n_s: %g\n",1.-6./16./_PI_*pow(dV/V,2)+2./8./_PI_*(ddV/V));
    fprintf(stdout," -> 1st-order slow-roll prediction for n_t: %g\n",-2./16./_PI_*pow(dV/V,2));
    }
    }
  */

  /** - compute H_pivot at phi_pivot */
  switch (ppm->primordial_spec_type) {

  case inflation_V:
  case inflation_V_end:

    /** - check positivity and negative slope of potential in field pivot
        value, and find value of phi_dot and H for field's pivot value,
        assuming slow-roll attractor solution has been reached. If no
        solution, code will stop there. */

    if (ppm->primordial_verbose > 1)
      printf(" (search attractor at pivot)\n");

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
    break;

  case inflation_H:

    /** - check positivity and negative slope of \f$ H(\phi)\f$ in field pivot
        value, and get H_pivot */

    class_call_except(primordial_inflation_check_hubble(ppm,
                                                        ppm->phi_pivot,
                                                        &H_pivot,
                                                        &dH,
                                                        &ddH,
                                                        &dddH),
                      ppm->error_message,
                      ppm->error_message,
                      free(y);free(y_ini);free(dy));
    break;

  default:
    free(y);free(y_ini);free(dy);
    class_stop(ppm->error_message,"ppm->primordial_spec_type=%d different from possible relevant cases",ppm->primordial_spec_type);
    break;
  }

  /** - find a_pivot, value of scale factor when k_pivot crosses horizon while phi=phi_pivot */

  a_pivot = ppm->k_pivot/H_pivot;

  /** - integrate background solution starting from phi_pivot and until
      k_max>>aH. This ensures that the inflationary model considered
      here is valid and that the primordial spectrum can be
      computed. Otherwise, if slow-roll brakes too early, model is not
      suitable and run stops. */

  if (ppm->primordial_verbose > 1)
    printf(" (check inflation duration after phi_pivot=%e)\n",ppm->phi_pivot);

  k_max = exp(ppm->lnk[ppm->lnk_size-1]);
  aH_end = k_max/ppr->primordial_inflation_ratio_max;

  y[ppm->index_in_a] = a_pivot;
  y[ppm->index_in_phi] = ppm->phi_pivot;
  if ((ppm->primordial_spec_type == inflation_V) || (ppm->primordial_spec_type == inflation_V_end))
    y[ppm->index_in_dphi] = a_pivot*dphidt_pivot;

  class_call_except(primordial_inflation_evolve_background(ppm,
                                                           ppr,
                                                           y,
                                                           dy,
                                                           _aH_,
                                                           aH_end,
                                                           _TRUE_,
                                                           forward,
                                                           conformal),
                    ppm->error_message,
                    ppm->error_message,
                    free(y);free(y_ini);free(dy));

  /* we need to do the opposite: to check that there is an initial
     time such that k_min << (aH)_ini. A guess is made by integrating
     backward in time. This can be done exactly for inflation_H, or
     only approximately for inflation_V (using the first-order
     approximation to the attractor inflationary solution). However
     this approximation is irrelevant because nevertheless, later on,
     we compute the attractor solution at the initial time with high
     accuracy, and then we integrate the background equations forward
     in time. Hence the approximation made here introduces zero
     mistake on the final result. It is just a way to find quickly a
     reasonable initial phi value. In the inflation_V case, if the
     exact forward integration reveals that the guess was not good
     (i.e. does not correspond to "early enough"), we iterate over
     sequences of backward/forward integration, until a correct time is
     found. For potential such that no solution exists (no long-enough
     slow-roll period before the pivot scale), the run stops. */

  if (ppm->primordial_verbose > 1)
    printf(" (check inflation duration before pivot)\n");

  k_min = exp(ppm->lnk[0]);
  aH_ini = k_min/ppr->primordial_inflation_ratio_min;

  switch (ppm->primordial_spec_type) {

  case inflation_V:
  case inflation_V_end:

    counter = 0;

    y[ppm->index_in_a] = a_pivot;
    y[ppm->index_in_phi] = ppm->phi_pivot;

    do {

      /* counter to avoid infinite loop */
      counter ++;

      class_test_except(counter >= ppr->primordial_inflation_phi_ini_maxit,
                        ppm->error_message,
                        free(y);free(y_ini);free(dy),
                        "when searching for an initial value of phi just before observable inflation takes place, could not converge after %d iterations. The potential does not allow eough inflationary e-folds before reaching the pivot scale",
                        counter);

      /* try to find a value phi_try such that
         aH=aH_ini*(ppr->primordial_inflation_aH_ini_target) (default:
         aH_ini*0.9). But this is using the approximate backward
         solution. So, anyway, we will check using the exact forward
         solution that at this phi_try, we really have aH < aH_ini; if
         this is not the case, we will iterate until a correct phi_try
         is found. */

      class_call_except(primordial_inflation_evolve_background(ppm,
                                                               ppr,
                                                               y,
                                                               dy,
                                                               _aH_,
                                                               aH_ini*ppr->primordial_inflation_aH_ini_target,
                                                               _TRUE_,
                                                               backward,
                                                               conformal),
                        ppm->error_message,
                        ppm->error_message,
                        free(y);free(y_ini);free(dy));

      phi_try = y[ppm->index_in_phi];

      /* in inflation_V case, find the accurate attractor solution for
         phi_ini', and then the correct value of a_ini, and finally of
         dphi/dtau_ini */

      /* find dphi/dt_ini (unlike dphi/dtau_ini, this does not depend on normalization of a) */
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

      /* we need to normalize a properly so that a=a_pivot when
         phi=phi_pivot. To do so, we evolve starting arbitrarily from
         a_ini=1, and then we rescale a_ini appropriately. */

      y[ppm->index_in_a] = 1.;
      y[ppm->index_in_phi] = phi_try;
      y[ppm->index_in_dphi] = y[ppm->index_in_a]*dphidt_try; // dphi/dtau = a dphi/dt

      class_call_except(primordial_inflation_evolve_background(ppm,
                                                               ppr,
                                                               y,
                                                               dy,
                                                               _phi_,
                                                               ppm->phi_pivot,
                                                               _TRUE_,
                                                               forward,
                                                               conformal),
                        ppm->error_message,
                        ppm->error_message,
                        free(y);free(y_ini);free(dy));

      /* now impose the correct a_ini */
      a_try = a_pivot/y[ppm->index_in_a];

      /* in case another iteration will be needed, set a new starting point for the routine primordial_inflation_evolve_background(...,backward) */
      y[ppm->index_in_a] = a_try;
      y[ppm->index_in_phi] = phi_try;

    } while (a_try*H_try > aH_ini);

    y_ini[ppm->index_in_a] = a_try;
    y_ini[ppm->index_in_phi] = phi_try;
    y_ini[ppm->index_in_dphi] = y_ini[ppm->index_in_a]*dphidt_try; // dphi/dtau = a dphi/dt

    break;

  case inflation_H:

    y[ppm->index_in_a] = a_pivot;
    y[ppm->index_in_phi] = ppm->phi_pivot;

    class_call_except(primordial_inflation_evolve_background(ppm,
                                                             ppr,
                                                             y,
                                                             dy,
                                                             _aH_,
                                                             aH_ini,
                                                             _TRUE_,
                                                             backward,
                                                             conformal),
                      ppm->error_message,
                      ppm->error_message,
                      free(y);free(y_ini);free(dy));

    y_ini[ppm->index_in_a] = y[ppm->index_in_a];
    y_ini[ppm->index_in_phi] = y[ppm->index_in_phi];

    break;

  default:
    free(y);free(y_ini);free(dy);
    class_stop(ppm->error_message,"ppm->primordial_spec_type=%d different from possible relevant cases",ppm->primordial_spec_type);
    break;
  }

  /** - starting from this time, i.e. from y_ini[ ], we run the routine
      which takes care of computing the primordial spectrum. */

  if (ppm->primordial_verbose > 1)
    printf(" (compute spectrum)\n");

  if (ppm->behavior == numerical) {

    class_call_except(primordial_inflation_spectra(ppt,
                                                   ppm,
                                                   ppr,
                                                   y_ini),
                      ppm->error_message,
                      ppm->error_message,
                      free(y);free(y_ini);free(dy));
  }
  else if (ppm->behavior == analytical) {

    class_call_except(primordial_inflation_analytic_spectra(ppt,
                                                            ppm,
                                                            ppr,
                                                            y_ini),
                      ppm->error_message,
                      ppm->error_message,
                      free(y);free(y_ini);free(dy));
  }
  else {
    class_stop(ppm->error_message,"Uncomprehensible value of the flag ppm->behavior=%d",ppm->behavior);
  }

  /** - before ending, we want to compute and store the values of \f$ \phi \f$
      corresponding to k=aH for k_min and k_max */

  y[ppm->index_in_a] = y_ini[ppm->index_in_a];
  y[ppm->index_in_phi] = y_ini[ppm->index_in_phi];
  if ((ppm->primordial_spec_type == inflation_V) || (ppm->primordial_spec_type == inflation_V_end))
    y[ppm->index_in_dphi] = y_ini[ppm->index_in_dphi];

  class_call_except(primordial_inflation_evolve_background(ppm,
                                                           ppr,
                                                           y,
                                                           dy,
                                                           _aH_,
                                                           k_min,
                                                           _FALSE_,
                                                           forward,
                                                           conformal),
                    ppm->error_message,
                    ppm->error_message,
                    free(y);free(y_ini);free(dy));

  ppm->phi_min=y[ppm->index_in_phi];

  class_call_except(primordial_inflation_evolve_background(ppm,
                                                           ppr,
                                                           y,
                                                           dy,
                                                           _aH_,
                                                           k_max,
                                                           _FALSE_,
                                                           forward,
                                                           conformal),
                    ppm->error_message,
                    ppm->error_message,
                    free(y);free(y_ini);free(dy));

  ppm->phi_max=y[ppm->index_in_phi];

  if (ppm->primordial_verbose > 1)
    printf(" (observable power spectrum goes from %e to %e)\n",
           ppm->phi_min,
           ppm->phi_max);

  /** - finally, we can de-allocate */

  free(y);
  free(y_ini);
  free(dy);

  return _SUCCESS_;
}

/**
 * Routine for the computation of an analytic apporoximation to the
 * the primordial spectrum. In general, should be used only for
 * comparing with exact numerical computation performed by
 * primordial_inflation_spectra().
 *
 * @param ppt   Input: pointer to perturbation structure
 * @param ppm   Input/output: pointer to primordial structure
 * @param ppr   Input: pointer to precision structure
 * @param y_ini Input: initial conditions for the vector of background/perturbations, already allocated and filled
 * @return the error status
 */

int primordial_inflation_analytic_spectra(
                                          struct perturbations * ppt,
                                          struct primordial * ppm,
                                          struct precision * ppr,
                                          double * y_ini
                                          ) {
  double * y;
  double * dy;
  int index_k;
  double k,phi_k;
  double curvature,tensors;
  double V,dV,ddV;

  /** Summary */
  /** - allocate vectors for background/perturbed quantities */
  class_alloc(y,ppm->in_size*sizeof(double),ppm->error_message);
  class_alloc(dy,ppm->in_size*sizeof(double),ppm->error_message);

  /** - initialize the background part of the running vector */
  y[ppm->index_in_a] = y_ini[ppm->index_in_a];
  y[ppm->index_in_phi] = y_ini[ppm->index_in_phi];
  if ((ppm->primordial_spec_type == inflation_V) || (ppm->primordial_spec_type == inflation_V_end))
    y[ppm->index_in_dphi] = y_ini[ppm->index_in_dphi];

  /** - loop over Fourier wavenumbers */
  for (index_k=0; index_k < ppm->lnk_size; index_k++) {

    k = exp(ppm->lnk[index_k]);

    /* evolve background until k=aH is reached */
    class_call(primordial_inflation_evolve_background(ppm,
                                                      ppr,
                                                      y,
                                                      dy,
                                                      _aH_,
                                                      k,
                                                      _FALSE_,
                                                      forward,
                                                      conformal),
               ppm->error_message,
               ppm->error_message);

    /** - read value of phi at time when k=aH */
    phi_k = y[ppm->index_in_phi];

    /** - get potential (and its derivatives) at this value */
    class_call(primordial_inflation_check_potential(ppm,phi_k,&V,&dV,&ddV),
               ppm->error_message,
               ppm->error_message);

    /** - calculate the analytic slow-roll formula for the spectra */
    curvature = 128.*_PI_/3.*pow(V,3)/pow(dV,2);
    tensors = pow(dV/V,2)/_PI_*128.*_PI_/3.*pow(V,3)/pow(dV,2);

    /** - store the obtained result for curvature and tensor perturbations */
    ppm->lnpk[ppt->index_md_scalars][index_k] = log(curvature);
    ppm->lnpk[ppt->index_md_tensors][index_k] = log(tensors);
  }

  ppm->is_non_zero[ppt->index_md_scalars][ppt->index_ic_ad] = _TRUE_;
  ppm->is_non_zero[ppt->index_md_tensors][ppt->index_ic_ten] = _TRUE_;

  return _SUCCESS_;
}

/**
 * Routine with a loop over wavenumbers for the computation of the primordial
 * spectrum. For each wavenumber it calls primordial_inflation_one_wavenumber()
 *
 * @param ppt   Input: pointer to perturbation structure
 * @param ppm   Input/output: pointer to primordial structure
 * @param ppr   Input: pointer to precision structure
 * @param y_ini Input: initial conditions for the vector of background/perturbations, already allocated and filled
 * @return the error status
 */

int primordial_inflation_spectra(
                                 struct perturbations * ppt,
                                 struct primordial * ppm,
                                 struct precision * ppr,
                                 double * y_ini
                                 ) {
  int index_k;

  class_setup_parallel();
  /* loop over Fourier wavenumbers */
  for (index_k=0; index_k < ppm->lnk_size; index_k++) {

    class_run_parallel(with_arguments(ppt,ppm,ppr,y_ini,index_k),

    class_call(primordial_inflation_one_wavenumber(ppt,ppm,ppr,y_ini,index_k),
               ppm->error_message,
               ppm->error_message);
    return _SUCCESS_;
    );
  }

  class_finish_parallel();

  ppm->is_non_zero[ppt->index_md_scalars][ppt->index_ic_ad] = _TRUE_;
  ppm->is_non_zero[ppt->index_md_tensors][ppt->index_ic_ten] = _TRUE_;

  return _SUCCESS_;

}

/**
 * Routine coordinating the computation of the primordial
 * spectrum for one wavenumber. It calls primordial_inflation_one_k() to
 * integrate the perturbation equations, and then it stores the result
 * for the scalar/tensor spectra.
 *
 * @param ppt     Input: pointer to perturbation structure
 * @param ppm     Input/output: pointer to primordial structure
 * @param ppr     Input: pointer to precision structure
 * @param y_ini   Input: initial conditions for the vector of background/perturbations, already allocated and filled
 * @param index_k Input: index of wavenumber to be considered
 * @return the error status
 */

int primordial_inflation_one_wavenumber(
                                        struct perturbations * ppt,
                                        struct primordial * ppm,
                                        struct precision * ppr,
                                        double * y_ini,
                                        int index_k
                                        ) {
  double k;
  double curvature,tensors;
  double * y;
  double * dy;

  k = exp(ppm->lnk[index_k]);

  /** Summary */
  /** - allocate vectors for background/perturbed quantities */
  class_alloc(y,ppm->in_size*sizeof(double),ppm->error_message);
  class_alloc(dy,ppm->in_size*sizeof(double),ppm->error_message);

  /** - initialize the background part of the running vector */
  y[ppm->index_in_a] = y_ini[ppm->index_in_a];
  y[ppm->index_in_phi] = y_ini[ppm->index_in_phi];
  if ((ppm->primordial_spec_type == inflation_V) || (ppm->primordial_spec_type == inflation_V_end))
    y[ppm->index_in_dphi] = y_ini[ppm->index_in_dphi];

  /** - evolve the background until the relevant initial time for
      integrating perturbations */
  class_call(primordial_inflation_evolve_background(ppm,
                                                    ppr,
                                                    y,
                                                    dy,
                                                    _aH_,
                                                    k/ppr->primordial_inflation_ratio_min,
                                                    _FALSE_,
                                                    forward,
                                                    conformal),
             ppm->error_message,
             ppm->error_message);

  /** - evolve the background/perturbation equations from this time and
      until some time after Horizon crossing */
  class_call(primordial_inflation_one_k(ppm,
                                        ppr,
                                        k,
                                        y,
                                        dy,
                                        &curvature,
                                        &tensors),
             ppm->error_message,
             ppm->error_message);

  free(y);
  free(dy);

  class_test(curvature<=0.,
             ppm->error_message,
             "negative curvature spectrum");

  class_test(tensors<=0.,
             ppm->error_message,
             "negative tensor spectrum");

  /** - store the obtained result for curvature and tensor perturbations */
  ppm->lnpk[ppt->index_md_scalars][index_k] = log(curvature);
  ppm->lnpk[ppt->index_md_tensors][index_k] = log(tensors);

  /* uncomment if you want to print here the spectra for testing */
  /* fprintf(stderr,"%e %e %e\n", */
  /* 	    ppm->lnk[index_k], */
  /* 	    ppm->lnpk[ppt->index_md_scalars][index_k], */
  /* 	    ppm->lnpk[ppt->index_md_tensors][index_k]); */

  return _SUCCESS_;
}

/**
 * Routine integrating the background plus perturbation equations for
 * each wavenumber, and returning the scalar and tensor spectrum.
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

  /** Summary: */

  /** - define local variables */
  double tau_start,tau_end,dtau;
  double z,ksi2,ah2;
  double aH;
  double curvature_old;
  double curvature_new;
  double dlnPdN;

  struct primordial_inflation_parameters_and_workspace pipaw;
  struct generic_integrator_workspace gi;

  /** - initialize the generic integrator (same integrator already used
      in background, thermodynamics and perturbation modules) */

  pipaw.ppm = ppm;
  pipaw.N = ppm->in_size;
  pipaw.integrate = forward;
  pipaw.time = conformal;
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

  /** - initialize variable used for deciding when to stop the calculation (= when the curvature remains stable) */
  curvature_new = _HUGE_;

  /** - initialize conformal time to arbitrary value (here, only variations
      of tau matter: the equations that we integrate do not depend
      explicitly on time) */
  tau_end = 0;

  /** - compute derivative of initial vector and infer first value of adaptive time-step */
  class_call(primordial_inflation_derivs(tau_end,
                                         y,
                                         dy,
                                         &pipaw,
                                         ppm->error_message),
             ppm->error_message,
             ppm->error_message);

  dtau = ppr->primordial_inflation_pt_stepsize*2.*_PI_
    /MAX(sqrt(fabs(dy[ppm->index_in_dksi_re]/y[ppm->index_in_ksi_re])),k);

  /** - loop over time */
  do {

    /*  new time interval [tau_start, tau_end] over which equations will be integrated */
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
    dtau = ppr->primordial_inflation_pt_stepsize*2.*_PI_
      /MAX(sqrt(fabs(dy[ppm->index_in_dksi_re]/y[ppm->index_in_ksi_re])),k);

    /* new aH */
    aH = dy[ppm->index_in_a]/y[ppm->index_in_a];

    /* store previous value of curvature (at tau_start) */
    curvature_old =  curvature_new;

    /* new curvature */
    z = y[ppm->index_in_a]*dy[ppm->index_in_phi]/aH;
    ksi2 = y[ppm->index_in_ksi_re]*y[ppm->index_in_ksi_re]+y[ppm->index_in_ksi_im]*y[ppm->index_in_ksi_im];
    curvature_new = k*k*k/2./_PI_/_PI_*ksi2/z/z;

    /* variation of curvature with time (dimensionless) */
    dlnPdN = (curvature_new-curvature_old)/dtau*y[ppm->index_in_a]/dy[ppm->index_in_a]/curvature_new;

    /* stop when (k >> aH) AND curvature is stable */
  } while ((k/aH >= ppr->primordial_inflation_ratio_max) || (fabs(dlnPdN) > ppr->primordial_inflation_tol_curvature));

  /** - clean the generic integrator */
  class_call(cleanup_generic_integrator(&gi),
             gi.error_message,
             ppm->error_message);

  /** - store final value of curvature for this wavenumber */
  *curvature = curvature_new;

  /** - store final value of tensor perturbation for this wavenumber */
  ah2 = y[ppm->index_in_ah_re]*y[ppm->index_in_ah_re]+y[ppm->index_in_ah_im]*y[ppm->index_in_ah_im];
  *tensor = 32.*k*k*k/_PI_*ah2/y[ppm->index_in_a]/y[ppm->index_in_a];

  //fprintf(stdout,"%g %g %g %g %g\n",k,*curvature,*tensor,*tensor/(*curvature),dlnPdN);

  return _SUCCESS_;
}

/**
 * Routine searching for the inflationary attractor solution at a
 * given phi_0, by iterations, with a given tolerance. If no solution
 * found within tolerance, returns error message. The principle is the
 * following. The code starts integrating the background equations
 * from various values of phi, corresponding to earlier and earlier
 * value before phi_0, and separated by a small arbitrary step size,
 * corresponding roughly to 1 e-fold of inflation. Each time, the
 * integration starts with the initial condition \f$ \phi=-V'/3H\f$ (slow-roll
 * prediction). If the found value of \f$\phi'\f$ in phi_0 is stable (up to
 * the parameter "precision"), the code considers that there is an
 * attractor, and stops iterating. If this process does not converge,
 * it returns an error message.
 *
 * @param ppm       Input: pointer to primordial structure
 * @param ppr       Input: pointer to precision structure
 * @param phi_0     Input: field value at which we wish to find the solution
 * @param precision Input: tolerance on output values (if too large, an attractor will always considered to be found)
 * @param y         Input: running vector of background variables, already allocated and initialized
 * @param dy        Input: running vector of background derivatives, already allocated
 * @param H_0       Output: Hubble value at phi_0 for attractor solution
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

  /* we want a series of value of phi' in phi_0, obtained by
     integrating the system from earlier and earlier time. The first
     value iof the series is the slow-roll prediction phi' =
     -V'/3H. The following lines compute this value and initialize relevant quantities. */

  class_call(primordial_inflation_check_potential(ppm,phi_0,&V_0,&dV_0,&ddV_0),
             ppm->error_message,
             ppm->error_message);

  dphidt_0new = -dV_0/3./sqrt((8.*_PI_/3.)*V_0);
  phi = phi_0;
  counter = 0;

  dphidt_0old = dphidt_0new/(precision+2.); // this silly value just
                                            // ensures that the loop
                                            // below will be executed
                                            // at least once.

  /* loop over different values of phi, from which the background
     equations are integrated until phi_0 */

  while (fabs(dphidt_0new/dphidt_0old-1.) >= precision) {

    counter ++;
    class_test(counter >= ppr->primordial_inflation_attractor_maxit,
               ppm->error_message,
               "could not converge after %d iterations: there exists no attractor solution near phi=%g. Potential probably too steep in this region, or precision parameter primordial_inflation_attractor_precision=%g too small",
               counter,
               phi_0,
               precision);

    dphidt_0old = dphidt_0new;

    /* take one step in phi, corresponding roughly to adding one more
       e-fold of inflation */

    phi=phi+dV_0/V_0/16./_PI_;

    /* fix the initial phi' to the slow-roll prediction in that point,
       and initialize other relevant quantities */

    class_call(primordial_inflation_check_potential(ppm,phi,&V,&dV,&ddV),
               ppm->error_message,
               ppm->error_message);

    a = 1.;
    dphidt = -dV/3./sqrt((8.*_PI_/3.)*V);
    y[ppm->index_in_a]=a;
    y[ppm->index_in_phi]=phi;
    y[ppm->index_in_dphi]=a*dphidt;

    /* evolve the background equations until phi_0 is reached */

    class_call(primordial_inflation_evolve_background(ppm,
                                                      ppr,
                                                      y,
                                                      dy,
                                                      _phi_,
                                                      phi_0,
                                                      _TRUE_,
                                                      forward,
                                                      conformal),
               ppm->error_message,
               ppm->error_message);

    /* compute phi' in phi_0, this is the new point in the series
       which convergence we want to check */

    dphidt_0new = y[ppm->index_in_dphi]/y[ppm->index_in_a];

  }

  /* if we have converged and found the attractor, we take the last
     value of phi' in phi_0 to be the correct one for the attractor
     solution */

  *dphidt_0 = dphidt_0new;
  *H_0 = sqrt((8.*_PI_/3.)*(0.5*dphidt_0new*dphidt_0new+V_0));

  if (ppm->primordial_verbose > 1) {
    printf(" (attractor found in phi=%g with phi'=%g, H=%g)\n",phi_0,*dphidt_0,*H_0);
  }

  return _SUCCESS_;
}

/**
 * Routine integrating background equations only, from initial values
 * stored in y, to a final value (if target = _aH_, until aH =
 * aH_stop; if target = _phi_, till phi = phi_stop; if target =
 * _end_inflation_, until \f$ d^2a/dt^2 = 0\f$ (here t = proper time)). In
 * output, y contains the final background values. In addition, if
 * check_epsilon is true, the routine controls at each step that the
 * expansion is accelerated and that inflation holds (wepsilon>1),
 * otherwise it returns an error. Thanks to the last argument, it is
 * also possible to specify whether the integration should be carried
 * forward or backward in time. For the inflation_H case, only a 1st
 * order differential equation is involved, so the forward and
 * backward case can be done exactly without problems. For the
 * inflation_V case, the equation of motion is 2nd order. What the
 * module will do in the backward case is to search for an approximate
 * solution, corresponding to the (first-order) attractor inflationary
 * solution. This approximate backward solution is used in order to
 * estimate some initial times, but the approximation made here will
 * never impact the final result: the module is written in such a way
 * that after using this approximation, the code always computes (and
 * relies on) the exact forward solution.
 *
 * @param ppm           Input: pointer to primordial structure
 * @param ppr           Input: pointer to precision structure
 * @param y             Input/output: running vector of background variables, already allocated and initialized
 * @param dy            Input: running vector of background derivatives, already allocated
 * @param target        Input: whether the goal is to reach a given aH or \f$ \phi \f$
 * @param stop          Input: the target value of either aH or \f$ \phi \f$
 * @param check_epsilon Input: whether we should impose inflation (epsilon>1) at each step
 * @param direction     Input: whether we should integrate forward or backward in time
 * @param time          Input: definition of time (proper or conformal)
 * @return the error status
 */

int primordial_inflation_evolve_background(
                                           struct primordial * ppm,
                                           struct precision * ppr,
                                           double * y,
                                           double * dy,
                                           enum target_quantity target,
                                           double stop,
                                           short check_epsilon,
                                           enum integration_direction direction,
                                           enum time_definition time
                                           ) {

  struct primordial_inflation_parameters_and_workspace pipaw;
  struct generic_integrator_workspace gi;
  double tau_start,tau_end,dtau=0.;
  double H,dH,ddH,dddH;
  double epsilon,epsilon_old;
  double quantity=0.;
  double V,dV,ddV;
  double sign_dtau=0.;

  pipaw.ppm = ppm;

  pipaw.N = ppm->in_bg_size;

  if ((direction == backward) && ((ppm->primordial_spec_type == inflation_V) || (ppm->primordial_spec_type == inflation_V_end))) {
    // -1 to remove the differential equation for phi', since we stick to the attractor
    pipaw.N -= 1;
  }

  pipaw.integrate = direction;
  pipaw.time = time;

  switch (direction) {
  case forward:
    sign_dtau = 1.;
    break;
  case backward:
    sign_dtau = -1.;
    break;
  }

  class_call(initialize_generic_integrator(pipaw.N,&gi),
             gi.error_message,
             ppm->error_message);

  /* at starting point, compute eventually epsilon */

  if (check_epsilon == _TRUE_) {

    class_call(primordial_inflation_get_epsilon(ppm,
                                                y[ppm->index_in_phi],
                                                &epsilon),
               ppm->error_message,
               ppm->error_message);
  }

  /* at starting point, compute the stepsize dtau */

  tau_end = 0;

  class_call(primordial_inflation_derivs(tau_end,
                                         y,
                                         dy,
                                         &pipaw,
                                         ppm->error_message),
             ppm->error_message,
             ppm->error_message);

  // compute timestep (if time = conformal, dtau is the conformal time step,
  // if time = proper, dtau is in fact dt, the proper time step)
  if ((direction == forward) && ((ppm->primordial_spec_type == inflation_V) || (ppm->primordial_spec_type == inflation_V_end))) {
    dtau = ppr->primordial_inflation_bg_stepsize
      *MIN(y[ppm->index_in_a]/dy[ppm->index_in_a],fabs(y[ppm->index_in_dphi]/dy[ppm->index_in_dphi]));
  }
  else {
    // minus sign for backward in time
    dtau = sign_dtau * ppr->primordial_inflation_bg_stepsize*y[ppm->index_in_a]/dy[ppm->index_in_a];
  }

  /* expected value of target quantity after the next step */
  switch (target) {
  case _aH_:
    // next (approximate) value of aH after next step
    // (a+[da/dx]*dx) H = aH (1 + [da/dx] / a dx)
    // where dtau can be conformal or proper time
    quantity = dy[ppm->index_in_a] * (1.+ dy[ppm->index_in_a]/y[ppm->index_in_a] * dtau);
    if (time == conformal) quantity /= y[ppm->index_in_a];
    break;
  case _phi_:
    // next (approximate) value of phi after next step
    quantity = y[ppm->index_in_phi]+dy[ppm->index_in_phi]*dtau;
    break;
  case _end_inflation_:
    // in this case, the goal is to reach d2a/dt2 = 0 (end of accelerated expansion)
    stop = 0.;
    // current value of quantity = - d2a/dt2 /a = [- (a'/a)^2 + 3/2 8pi/3 phi'^2]/a^2
    quantity = -pow(dy[ppm->index_in_a]/y[ppm->index_in_a],2) + 4*_PI_ *  y[ppm->index_in_dphi] * y[ppm->index_in_dphi];
    if (time == conformal) quantity /= pow(y[ppm->index_in_a],2);

    // check that we are in the right case
    class_test(ppm->primordial_spec_type != inflation_V_end,
               ppm->error_message,
               "the target _end_inflation_ is only coded to work with inflation_V_end (but could be generalized if needed)");
    break;
  case _a_:
    // next (approximate) value of a after next step
    quantity = y[ppm->index_in_a]+dy[ppm->index_in_a]*dtau;
    break;
  }

  /* loop over time steps, checking that there will be no overshooting */

  while (sign_dtau*(quantity - stop) < 0.) {

    /* check that V(phi) or H(phi) do not take forbidden values
       (negative or positive derivative) */

    if ((ppm->primordial_spec_type == inflation_V) || (ppm->primordial_spec_type == inflation_V_end)) {
      class_call(primordial_inflation_check_potential(ppm,
                                                      y[ppm->index_in_phi],
                                                      &V,
                                                      &dV,
                                                      &ddV),
                 ppm->error_message,
                 ppm->error_message);
    }
    else {
      class_call(primordial_inflation_check_hubble(ppm,
                                                   y[ppm->index_in_phi],
                                                   &H,
                                                   &dH,
                                                   &ddH,
                                                   &dddH),
                 ppm->error_message,
                 ppm->error_message);
    }

    /* take one time step */

    tau_start = tau_end;

    tau_end = tau_start + dtau;

    // mind the fabs(...) below (works for both forward and backward integration)
    class_test(fabs(dtau/tau_start) < ppr->smallest_allowed_variation,
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

    /* eventually, check that epsilon is not becoming greater than one */

    if (check_epsilon == _TRUE_) {

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

    /* recompute new value of next conformal time step */

    class_call(primordial_inflation_derivs(tau_end,
                                           y,
                                           dy,
                                           &pipaw,
                                           ppm->error_message),
               ppm->error_message,
               ppm->error_message);

    // compute timestep (if time = conformal, dtau is the conformal time step,
    // if time = proper, dtau is in fact dt, the proper time step)
    if ((direction == forward) && ((ppm->primordial_spec_type == inflation_V) || (ppm->primordial_spec_type == inflation_V_end))) {
      dtau = ppr->primordial_inflation_bg_stepsize
        *MIN(y[ppm->index_in_a]/dy[ppm->index_in_a],fabs(y[ppm->index_in_dphi]/dy[ppm->index_in_dphi]));
    }
    else {
      // minus sign for backward in time
      dtau = sign_dtau * ppr->primordial_inflation_bg_stepsize*y[ppm->index_in_a]/dy[ppm->index_in_a];
    }

    /* expected value of target quantity after the next step */

    switch (target) {
    case _aH_:
      // next (approximate) value of aH after next step
      // (a+[da/dx]*dx) H = aH (1 + [da/dx] / a dx)
      // where dtau can be conformal or proper time
      quantity = dy[ppm->index_in_a] * (1.+ dy[ppm->index_in_a]/y[ppm->index_in_a] * dtau);
      if (time == conformal) quantity /= y[ppm->index_in_a];
      break;
    case _phi_:
      // next (approximate) value of phi after next step
      quantity = y[ppm->index_in_phi]+dy[ppm->index_in_phi]*dtau;
      break;
    case _end_inflation_:
      // current value of quantity = - d2a/dt2 /a = [- (a'/a)^2 + 3/2 8pi/3 phi'^2]/a^2
      quantity = -pow(dy[ppm->index_in_a]/y[ppm->index_in_a],2) + 4*_PI_ *  y[ppm->index_in_dphi] * y[ppm->index_in_dphi];
      if (time == conformal) quantity /= pow(y[ppm->index_in_a],2);
      break;
    case _a_:
      // next (approximate) value of a after next step
      quantity = y[ppm->index_in_a]+dy[ppm->index_in_a]*dtau;
      break;
    }

  }

  /* won't use the integrator anymore */

  class_call(cleanup_generic_integrator(&gi),
             gi.error_message,
             ppm->error_message);

  /* Perform one last step with a simple trapezoidal integral. This
     will bring exactly phi or a forward to phi_stop or a_stop, or
     approximately aH forward to aH_stop, or approximately [-d2a/dt2
     /a] backward to zero. */

  switch (target) {
  case _aH_:
    switch (time){
    case proper:
      dtau = (stop/dy[ppm->index_in_a]-1.)/dy[ppm->index_in_a];
      break;
    case conformal:
      dtau = (stop/(dy[ppm->index_in_a]/y[ppm->index_in_a])-1.)/(dy[ppm->index_in_a]/y[ppm->index_in_a]);
      break;
    }
    break;
  case _phi_:
    dtau = (stop-y[ppm->index_in_phi])/dy[ppm->index_in_phi];
    break;
  case _end_inflation_:
    class_call(primordial_inflation_check_potential(ppm,y[ppm->index_in_phi],&V,&dV,&ddV),
               ppm->error_message,
               ppm->error_message);
    // We can easily pull back quantity=-d2a/dt2 /a by noticing that
    // d(quantity)/dtau = 8piG phi' phi'' / a^2 (exact relation!)
    // or
    // d(quantity)/dtau = 8piG phi^dot (a phi^dot)^dot = 8piG phi^dot (a^dot phi^dot+ a phi^dotdot)
    // By taking the step dtau = - quantity / [d(quantity)/dtau] we nearly reach quantity=0 (end of inflation), up to very good approximation
    switch (time){
    case proper:
      dtau = -quantity/(8.*_PI_*dy[ppm->index_in_phi]*(dy[ppm->index_in_a]*dy[ppm->index_in_phi]+y[ppm->index_in_a]*dy[ppm->index_in_dphi]));

      break;
    case conformal:
      dtau = -quantity/(8.*_PI_/y[ppm->index_in_a]/y[ppm->index_in_a]*dy[ppm->index_in_phi]*dy[ppm->index_in_dphi]);
      break;
    }
    break;
  case _a_:
    dtau = (stop-y[ppm->index_in_a])/dy[ppm->index_in_a];
    break;
  }

  y[ppm->index_in_a] += dy[ppm->index_in_a]*dtau;
  y[ppm->index_in_phi] += dy[ppm->index_in_phi]*dtau;
  if ((direction == forward) && ((ppm->primordial_spec_type == inflation_V)||(ppm->primordial_spec_type == inflation_V_end)))
    y[ppm->index_in_dphi] += dy[ppm->index_in_dphi]*dtau;

  // this last step updates also the dy[]
  class_call(primordial_inflation_derivs(tau_end,
                                         y,
                                         dy,
                                         &pipaw,
                                         ppm->error_message),
             ppm->error_message,
             ppm->error_message);

  // uncomment if you want to test that the routine really reached the point at which d2a/dt2=0
  /*
    if (target == _end_inflation_) {
    class_call(primordial_inflation_derivs(tau_end,
    y,
    dy,
    &pipaw,
    ppm->error_message),
    ppm->error_message,
    ppm->error_message);

    aH = dy[ppm->index_in_a]/y[ppm->index_in_a];
    quantity = (-aH*aH + 4*_PI_ *  y[ppm->index_in_dphi] * y[ppm->index_in_dphi])/y[ppm->index_in_a]/y[ppm->index_in_a];
    if (ppm->primordial_verbose>1)
    printf(" (-d2a/dt2 /a = %e)\n",quantity);
    }
  */

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
 * @param V              Output: inflaton potential in units of \f$ Mp^4\f$
 * @param dV             Output: first derivative of inflaton potential wrt the field
 * @param ddV            Output: second derivative of inflaton potential wrt the field
 * @return the error status
 */

int primordial_inflation_check_potential(
                                         struct primordial * ppm,
                                         double phi,
                                         double * V,
                                         double * dV,
                                         double * ddV
                                         ) {

  class_call(primordial_inflation_potential(ppm,phi,V,dV,ddV),
             ppm->error_message,
             ppm->error_message);

  class_test(*V <= 0.,
             ppm->error_message,
             "This potential becomes negative at phi=%g, before the end of observable inflation. It  cannot be treated by this code",
             phi);

  class_test(*dV >= 0.,
             ppm->error_message,
             "All the code is written for the case dV/dphi<0. Here, in phi=%g, we have dV/dphi=%g. This potential cannot be treated by this code",
             phi,*dV);

  return _SUCCESS_;
}

/**
 * Routine checking positivity and negative slope of \f$ H(\phi)\f$. The
 * negative slope is an arbitrary choice. Currently the code can only
 * deal with monotonic variations of the inflaton during
 * inflation. And H can only decrease with time. So the slope \f$ dH/d\phi\f$
 * has to be always negative or always positive... we took the first
 * option: phi increases, H decreases.
 *
 * @param ppm       Input: pointer to primordial structure
 * @param phi       Input: field value where to perform the check
 * @param H         Output: Hubble parameters in units of Mp
 * @param dH        Output: \f$ dH / d\phi \f$
 * @param ddH       Output: \f$ d^2H / d\phi^2 \f$
 * @param dddH      Output: \f$ d^3H / d\phi^3 \f$
 * @return the error status
 */

int primordial_inflation_check_hubble(
                                      struct primordial * ppm,
                                      double phi,
                                      double * H,
                                      double * dH,
                                      double * ddH,
                                      double * dddH
                                      ) {

  class_call(primordial_inflation_hubble(ppm,
                                         phi,
                                         H,dH,ddH,dddH),
             ppm->error_message,
             ppm->error_message);

  class_test(*H < 0.,
             ppm->error_message,
             "this H(phi) is not physical. H = %e",
             *H);

  class_test(*dH > 0.,
             ppm->error_message,
             "this H(phi) is not decreasing with growing phi. dH/dphi = %e",
             *dH);

  return _SUCCESS_;

}

/**
 * Routine computing the first slow-roll parameter epsilon
 *
 * @param ppm       Input: pointer to primordial structure
 * @param phi       Input: field value where to compute epsilon
 * @param epsilon   Output: result
 * @return the error status
 */

int primordial_inflation_get_epsilon(
                                     struct primordial * ppm,
                                     double phi,
                                     double * epsilon
                                     ) {

  double V,dV,ddV;
  double H,dH,ddH,dddH;

  switch (ppm->primordial_spec_type) {

  case inflation_V:
  case inflation_V_end:

    class_call(primordial_inflation_potential(ppm,
                                              phi,
                                              &V,&dV,&ddV),
               ppm->error_message,
               ppm->error_message);

    *epsilon = 1./16./_PI_*pow(dV/V,2);
    //*eta = 1./8./pi*(ddV/V)
    break;

  case inflation_H:

    class_call(primordial_inflation_hubble(ppm,
                                           phi,
                                           &H,&dH,&ddH,&dddH),
               ppm->error_message,
               ppm->error_message);

    *epsilon = 1./4./_PI_*pow(dH/H,2);
    break;

  default:
    class_stop(ppm->error_message,"ppm->primordial_spec_type=%d different from possible relevant cases",ppm->primordial_spec_type);
    break;
  }

  return _SUCCESS_;
}

/**
 * Routine searching phi_pivot when a given amount of inflation is requested.
 *
 * @param ppm       Input/output: pointer to primordial structure
 * @param ppr       Input: pointer to precision structure
 * @param y         Input: running vector of background variables, already allocated and initialized
 * @param dy        Input: running vector of background derivatives, already allocated
 * @return the error status
 */

int primordial_inflation_find_phi_pivot(
                                        struct primordial * ppm,
                                        struct precision * ppr,
                                        double * y,
                                        double * dy
                                        ) {
  /** Summary: */

  /** - define local variables */
  double epsilon,dphi;
  double phi_try,H_try,dphidt_try,ratio_try=0.;
  double phi_left,phi_right,phi_mid;
  double phi_small_epsilon,phi_stop;
  double dphidt_small_epsilon;
  double H_small_epsilon;
  double aH_ratio_after_small_epsilon=0.;
  double a_ratio_after_small_epsilon=0.;
  double target=0.;
  double a_pivot,aH_pivot;

  double rho_end;
  double h;
  double H0;
  double rho_c0;
  double sigma_B;
  double Omega_g0;
  double Omega_r0;

  /** - check whether in vicinity of phi_end, inflation is still ongoing */

  class_call(primordial_inflation_get_epsilon(ppm,ppm->phi_end-ppr->primordial_inflation_end_dphi,&epsilon),
             ppm->error_message,
             ppm->error_message);

  /** - case in which epsilon>1: hence we must find the value phi_stop <
      phi_end where inflation ends up naturally */

  if (epsilon > 1.) {

    // assume that inflation ends up naturally

    /** - --> find latest value of the field such that epsilon = primordial_inflation_small_epsilon (default: 0.1) */

    /** - --> bracketing right-hand value is phi_end (but the potential will not be evaluated exactly there, only closeby */
    phi_right = ppm->phi_end;

    /** - --> bracketing left-hand value is found by iterating with logarithmic step until epsilon < primordial_inflation_small_epsilon */
    dphi = ppr->primordial_inflation_end_dphi;
    do {
      dphi *= ppr->primordial_inflation_end_logstep;
      class_call(primordial_inflation_get_epsilon(ppm,ppm->phi_end-dphi,&epsilon),
                 ppm->error_message,
                 ppm->error_message);
    } while (epsilon > ppr->primordial_inflation_small_epsilon);
    phi_left = ppm->phi_end-dphi;

    /** - --> find value such that epsilon = primordial_inflation_small_epsilon by bisection */
    do {
      phi_mid = 0.5*(phi_left+phi_right);
      class_call(primordial_inflation_get_epsilon(ppm,phi_mid,&epsilon),
                 ppm->error_message,
                 ppm->error_message);
      if (epsilon < ppr->primordial_inflation_small_epsilon) phi_left=phi_mid;
      else phi_right=phi_mid;
    } while (fabs(epsilon-ppr->primordial_inflation_small_epsilon) > ppr->primordial_inflation_small_epsilon_tol);

    /** - --> value found and stored as phi_small_epsilon */
    phi_small_epsilon = phi_mid;

    /** - --> find inflationary attractor in phi_small_epsilon (should exist since epsilon<<1 there) */
    class_call(primordial_inflation_find_attractor(ppm,
                                                   ppr,
                                                   phi_small_epsilon,
                                                   ppr->primordial_inflation_attractor_precision_initial,
                                                   y,
                                                   dy,
                                                   &H_small_epsilon,
                                                   &dphidt_small_epsilon),
               ppm->error_message,
               ppm->error_message);

    /** - --> compute amount of inflation between this phi_small_epsilon and the end of inflation */
    y[ppm->index_in_a]=1.;
    y[ppm->index_in_phi]= phi_small_epsilon;
    y[ppm->index_in_dphi]=y[ppm->index_in_a]*dphidt_small_epsilon;

    class_call(primordial_inflation_evolve_background(ppm,
                                                      ppr,
                                                      y,
                                                      dy,
                                                      _end_inflation_,
                                                      0.,
                                                      _FALSE_,
                                                      forward,
                                                      conformal),
               ppm->error_message,
               ppm->error_message);

    // we have used here conformal time, so aH = dy[a]/y[a]
    aH_ratio_after_small_epsilon = dy[ppm->index_in_a]/y[ppm->index_in_a]/H_small_epsilon;
    a_ratio_after_small_epsilon = y[ppm->index_in_a];

    switch (ppm->phi_pivot_method) {

    case ln_aH_ratio_auto:

      /* get the target value of ln_aH_ratio */

      rho_end = 2./8./_PI_*pow(dy[ppm->index_in_a]/y[ppm->index_in_a],2);
      rho_end = 8*_PI_/3.*rho_end/(_G_*_h_P_/pow(_c_,3))*pow(_Mpc_over_m_,2);
      h = 0.7;
      H0 = h * 1.e5 / _c_;
      rho_c0 = pow(H0,2);

      sigma_B = 2. * pow(_PI_,5) * pow(_k_B_,4) / 15. / pow(_h_P_,3) / pow(_c_,2);
      Omega_g0 = (4.*sigma_B/_c_*pow(2.726,4.)) / (3.*_c_*_c_*1.e10*h*h/_Mpc_over_m_/_Mpc_over_m_/8./_PI_/_G_);
      Omega_r0 = 3.044*7./8.*pow(4./11.,4./3.)*Omega_g0;

      target = log(H0/0.05*pow(Omega_r0,0.5)*pow(2./100.,1./12.)*pow(rho_end/rho_c0,0.25));

      //fprintf(stderr,"auto: log(aH_end/aH_*)=%e\n",target);
      break;

    case ln_aH_ratio:

      target = ppm->phi_pivot_target;
      //fprintf(stderr,"fixed: log(aH_end/aH_*)=%e\n",target);
      break;

    case N_star:

      target = ppm->phi_pivot_target;
      //fprintf(stderr,"fixed: log(a_end/a_*)=%e\n",target);
      break;
    }

    /** - --> by starting from phi_small_epsilon and integrating an approximate
        solution backward in time, try to estimate roughly a value close
        to phi_pivot but a bit smaller. This is done by trying to reach
        an amount of inflation equal to the requested one, minus the
        amount after phi_small_epsilon, and plus
        primordial_inflation_extra_efolds efolds (default: two). Note
        that it is not aggressive to require two extra e-folds of
        inflation before the pivot, since the calculation of the spectrum
        in the observable range will require even more. */

    y[ppm->index_in_a]=1.;
    y[ppm->index_in_phi]= phi_small_epsilon;

    switch (ppm->phi_pivot_method) {

    case ln_aH_ratio_auto:
    case ln_aH_ratio:

      class_call(primordial_inflation_evolve_background(ppm,
                                                        ppr,
                                                        y,
                                                        dy,
                                                        _aH_,
                                                        H_small_epsilon/exp(target+ppr->primordial_inflation_extra_efolds)*aH_ratio_after_small_epsilon,
                                                        _TRUE_,
                                                        backward,
                                                        conformal),
                 ppm->error_message,
                 ppm->error_message);
      break;

    case N_star:

      class_call(primordial_inflation_evolve_background(ppm,
                                                        ppr,
                                                        y,
                                                        dy,
                                                        _a_,
                                                        1./exp(target+ppr->primordial_inflation_extra_efolds)*a_ratio_after_small_epsilon,
                                                        _TRUE_,
                                                        backward,
                                                        conformal),
                 ppm->error_message,
                 ppm->error_message);
      break;
    }

    /* we now have a value phi_try believed to be close to and slightly smaller than phi_pivot */

    phi_try = y[ppm->index_in_phi];

    /** - --> find attractor in phi_try */

    class_call(primordial_inflation_find_attractor(ppm,
                                                   ppr,
                                                   phi_try,
                                                   ppr->primordial_inflation_attractor_precision_initial,
                                                   y,
                                                   dy,
                                                   &H_try,
                                                   &dphidt_try),
               ppm->error_message,
               ppm->error_message);

    /** - --> check the total amount of inflation between phi_try and the end of inflation */

    y[ppm->index_in_a]=1.;
    y[ppm->index_in_phi]= phi_try;
    y[ppm->index_in_dphi]= dphidt_try;

    class_call(primordial_inflation_evolve_background(ppm,
                                                      ppr,
                                                      y,
                                                      dy,
                                                      _end_inflation_,
                                                      0.,
                                                      _FALSE_,
                                                      forward,
                                                      proper),
               ppm->error_message,
               ppm->error_message);

    switch (ppm->phi_pivot_method) {

    case ln_aH_ratio_auto:
    case ln_aH_ratio:

      // aH_ratio (we have used here proper time, so aH = dy[a])
      ratio_try = dy[ppm->index_in_a]/H_try;
      break;

    case N_star:

      // a_ratio
      ratio_try = y[ppm->index_in_a];
      break;
    }

    class_test(log(ratio_try) < target,
               ppm->error_message,
               "phi_try not small enough, log(aH_stop/aH_try) or log(a_stop/a_try) (depending on what you asked) is equal to %e instead of requested %e; must write here a loop to deal automatically with this situation (by decreasing phi_try iteratively), or must increase precision parameter primordial_inflation_extra_efolds",
               log(ratio_try),
               target);

    phi_stop = y[1];

    if (ppm->primordial_verbose > 1)
      printf(" (inflation stops in phi_stop = %e)\n",phi_stop);

    /** - --> go back to phi_try, and now find phi_pivot such that the amount
        of inflation between phi_pivot and the end of inflation is
        exactly the one requested. */
    y[ppm->index_in_a]=1.;
    y[ppm->index_in_phi]= phi_try;
    y[ppm->index_in_dphi]= dphidt_try;

    switch (ppm->phi_pivot_method) {

    case ln_aH_ratio_auto:
    case ln_aH_ratio:

      class_call(primordial_inflation_evolve_background(ppm,
                                                        ppr,
                                                        y,
                                                        dy,
                                                        _aH_,
                                                        H_try*ratio_try/exp(target),
                                                        _FALSE_,
                                                        forward,
                                                        proper),
                 ppm->error_message,
                 ppm->error_message);
      break;

    case N_star:

      class_call(primordial_inflation_evolve_background(ppm,
                                                        ppr,
                                                        y,
                                                        dy,
                                                        _a_,
                                                        ratio_try/exp(target),
                                                        _FALSE_,
                                                        forward,
                                                        proper),
                 ppm->error_message,
                 ppm->error_message);
      break;
    }

    ppm->phi_pivot = y[1];

    if (ppm->primordial_verbose > 1) {

      printf(" (reached phi_pivot=%e)\n",ppm->phi_pivot);

      /* - --> In verbose mode, check that phi_pivot is correct. Done by
         restarting from phi_pivot and going again till the end of
         inflation. */

      aH_pivot = dy[0];
      a_pivot = y[0];
      class_call(primordial_inflation_evolve_background(ppm,
                                                        ppr,
                                                        y,
                                                        dy,
                                                        _end_inflation_,
                                                        0.,
                                                        _FALSE_,
                                                        forward,
                                                        proper),
                 ppm->error_message,
                 ppm->error_message);
      printf(" (from phi_pivot till the end, ln(aH_2/aH_1) = %e, ln(a_2/a_1) = %e)\n",log(dy[0]/aH_pivot),log(y[0]/a_pivot));
    }


  }
  /** - case in which epsilon<1: */
  else {

    /** - --> find inflationary attractor in phi_small_epsilon (should exist since epsilon<1 there) */
    class_call(primordial_inflation_find_attractor(ppm,
                                                   ppr,
                                                   ppm->phi_end,
                                                   ppr->primordial_inflation_attractor_precision_initial,
                                                   y,
                                                   dy,
                                                   &H_small_epsilon,
                                                   &dphidt_small_epsilon),
               ppm->error_message,
               ppm->error_message);

    /** - --> by starting from phi_end and integrating an approximate
        solution backward in time, try to estimate roughly a value close
        to phi_pivot but a bit smaller. This is done by trying to reach
        an amount of inflation equal to the requested one, minus the
        amount after phi_small_epsilon, and plus
        primordial_inflation_extra_efolds efolds (default: two). Note
        that it is not aggressive to require two extra e-folds of
        inflation before the pivot, since the calculation of the spectrum
        in the observable range will require even more. */

    y[ppm->index_in_a]=1.;
    y[ppm->index_in_phi]= ppm->phi_end;

    switch (ppm->phi_pivot_method) {

    case ln_aH_ratio_auto:
    case ln_aH_ratio:

      class_call(primordial_inflation_evolve_background(ppm,
                                                        ppr,
                                                        y,
                                                        dy,
                                                        _aH_,
                                                        H_small_epsilon/exp(target+ppr->primordial_inflation_extra_efolds)*aH_ratio_after_small_epsilon,
                                                        _TRUE_,
                                                        backward,
                                                        conformal),
                 ppm->error_message,
                 ppm->error_message);
      break;

    case N_star:

      class_call(primordial_inflation_evolve_background(ppm,
                                                        ppr,
                                                        y,
                                                        dy,
                                                        _a_,
                                                        1./exp(target+ppr->primordial_inflation_extra_efolds)*a_ratio_after_small_epsilon,
                                                        _TRUE_,
                                                        backward,
                                                        conformal),
                 ppm->error_message,
                 ppm->error_message);
      break;
    }

    /** - --> we now have a value phi_try believed to be close to and slightly smaller than phi_pivot */

    phi_try = y[ppm->index_in_phi];

    /** - --> find attractor in phi_try */

    class_call(primordial_inflation_find_attractor(ppm,
                                                   ppr,
                                                   phi_try,
                                                   ppr->primordial_inflation_attractor_precision_initial,
                                                   y,
                                                   dy,
                                                   &H_try,
                                                   &dphidt_try),
               ppm->error_message,
               ppm->error_message);

    /** - --> check the total amount of inflation between phi_try and the end of inflation */

    y[ppm->index_in_a]=1.;
    y[ppm->index_in_phi]= phi_try;
    y[ppm->index_in_dphi]= dphidt_try;

    class_call(primordial_inflation_evolve_background(ppm,
                                                      ppr,
                                                      y,
                                                      dy,
                                                      _phi_,
                                                      ppm->phi_end,
                                                      _FALSE_,
                                                      forward,
                                                      proper),
               ppm->error_message,
               ppm->error_message);

    switch (ppm->phi_pivot_method) {

    case ln_aH_ratio_auto:
    case ln_aH_ratio:

      // aH_ratio (we have used here proper time, so aH = dy[a])
      ratio_try = dy[ppm->index_in_a]/H_try;
      break;

    case N_star:

      // a_ratio
      ratio_try = y[ppm->index_in_a];
      break;
    }

    class_test(log(ratio_try) < target,
               ppm->error_message,
               "phi_try not small enough, log(aH_stop/aH_try) or log(a_stop/a_try) (depending on what you asked) is equal to %e instead of requested %e; must write here a loop to deal automatically with this situation (by decreasing phi_try iteratively), or must increase precision parameter primordial_inflation_extra_efolds",
               log(ratio_try),
               target);

    phi_stop = y[1];

    if (ppm->primordial_verbose > 1)
      printf(" (inflation stops in phi_stop = %e)\n",phi_stop);

    /** - --> go back to phi_try, and now find phi_pivot such that the amount
        of inflation between phi_pivot and the end of inflation is
        exactly the one requested. */
    y[ppm->index_in_a]=1.;
    y[ppm->index_in_phi]= phi_try;
    y[ppm->index_in_dphi]= dphidt_try;

    switch (ppm->phi_pivot_method) {

    case ln_aH_ratio_auto:
    case ln_aH_ratio:

      class_call(primordial_inflation_evolve_background(ppm,
                                                        ppr,
                                                        y,
                                                        dy,
                                                        _aH_,
                                                        H_try*ratio_try/exp(target),
                                                        _FALSE_,
                                                        forward,
                                                        proper),
                 ppm->error_message,
                 ppm->error_message);
      break;

    case N_star:

      class_call(primordial_inflation_evolve_background(ppm,
                                                        ppr,
                                                        y,
                                                        dy,
                                                        _a_,
                                                        ratio_try/exp(target),
                                                        _FALSE_,
                                                        forward,
                                                        proper),
                 ppm->error_message,
                 ppm->error_message);
      break;
    }

    ppm->phi_pivot = y[1];

    if (ppm->primordial_verbose > 1) {

      printf(" (reached phi_pivot=%e)\n",ppm->phi_pivot);

      /** - --> In verbose mode, check that phi_pivot is correct. Done by
          restarting from phi_pivot and going again till the end of
          inflation. */

      aH_pivot = dy[0];
      a_pivot = y[0];
      class_call(primordial_inflation_evolve_background(ppm,
                                                        ppr,
                                                        y,
                                                        dy,
                                                        _phi_,
                                                        ppm->phi_end,
                                                        _FALSE_,
                                                        forward,
                                                        proper),
                 ppm->error_message,
                 ppm->error_message);
      printf(" (from phi_pivot till the end, ln(aH_2/aH_1) = %e, ln(a_2/a_1) = %e)\n",log(dy[0]/aH_pivot),log(y[0]/a_pivot));
    }

  }

  return _SUCCESS_;
}

/**
 * Routine returning derivative of system of background/perturbation
 * variables. Like other routines used by the generic integrator
 * (background_derivs, thermodynamics_derivs, perturbations_derivs), this
 * routine has a generic list of arguments, and a slightly different
 * error management, with the error message returned directly in an
 * ErrMsg field.
 *
 * @param tau                      Input: time (not used explicitly inside the routine, but requested by the generic integrator)
 * @param y                        Input/output: running vector of background variables, already allocated and initialized
 * @param dy                       Input: running vector of background derivatives, already allocated
 * @param parameters_and_workspace Input: all necessary input variables apart from y
 * @param error_message            Output: error message
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

  ppipaw = (struct primordial_inflation_parameters_and_workspace *)parameters_and_workspace;
  ppm = ppipaw->ppm;

  // a2
  ppipaw->a2=y[ppm->index_in_a]*y[ppm->index_in_a];

  // BACKGROUND

  switch (ppm->primordial_spec_type) {

  case inflation_V:
  case inflation_V_end:

    class_call(primordial_inflation_potential(ppm,
                                              y[ppm->index_in_phi],
                                              &(ppipaw->V),
                                              &(ppipaw->dV),
                                              &(ppipaw->ddV)),
               ppm->error_message,
               ppm->error_message);

    switch (ppipaw->integrate) {

    case forward:

      switch (ppipaw->time) {

      case conformal:

        // a H = a'/a
        ppipaw->aH = sqrt((8*_PI_/3.)*(0.5*y[ppm->index_in_dphi]*y[ppm->index_in_dphi]+ppipaw->a2*ppipaw->V));
        // 1: a
        dy[ppm->index_in_a]=y[ppm->index_in_a]*ppipaw->aH;
        // 2: phi
        dy[ppm->index_in_phi]=y[ppm->index_in_dphi];
        // 3: dphi/dtau
        dy[ppm->index_in_dphi]=-2.*ppipaw->aH*y[ppm->index_in_dphi]-ppipaw->a2*ppipaw->dV;
        break;

      case proper:

        // a H = adot
        ppipaw->aH = y[ppm->index_in_a]*sqrt((8*_PI_/3.)*(0.5*y[ppm->index_in_dphi]*y[ppm->index_in_dphi]+ppipaw->V));
        // 1: a
        dy[ppm->index_in_a]=ppipaw->aH;
        // 2: phi
        dy[ppm->index_in_phi]=y[ppm->index_in_dphi];
        // 3: dphi/dt
        dy[ppm->index_in_dphi]=-3.*ppipaw->aH/y[ppm->index_in_a]*y[ppm->index_in_dphi]-ppipaw->dV;
        break;
      }

      // z''/z (assumes that conformal time is requested)
      ppipaw->zpp_over_z=
        2*ppipaw->aH*ppipaw->aH
        - ppipaw->a2*ppipaw->ddV
        - 4.*_PI_*(7.*y[ppm->index_in_dphi]*y[ppm->index_in_dphi]
                   +4.*y[ppm->index_in_dphi]/ppipaw->aH*ppipaw->a2*ppipaw->dV)
        +32.*_PI_*_PI_*pow(y[ppm->index_in_dphi],4)/pow(ppipaw->aH,2);

      // a''/a (assumes that conformal time is requested)
      ppipaw->app_over_a=2.*ppipaw->aH*ppipaw->aH - 4.*_PI_*y[ppm->index_in_dphi]*y[ppm->index_in_dphi];

      break;

      // For backward integration of approximate slow-roll solution:
      // Neglect kinetic energy of the field phi'^2/(2a^2) w.r.t. potential energy V
      // Neglect phi'' w.r.t 2aHphi', reducing 2nd order Klein-Gordon to approximate 1st-order
    case backward:

      switch (ppipaw->time) {

      case conformal:

        // a H = a'/a
        ppipaw->aH = sqrt((8*_PI_/3.)*ppipaw->a2*ppipaw->V);
        // 1: a
        dy[ppm->index_in_a]=y[ppm->index_in_a]*ppipaw->aH;
        // 2: phi
        dy[ppm->index_in_phi]= -ppipaw->a2*ppipaw->dV/3./ppipaw->aH;
        break;

      case proper:

        // a H = da/dt
        ppipaw->aH = y[ppm->index_in_a]*sqrt((8*_PI_/3.)*ppipaw->V);
        // 1: a
        dy[ppm->index_in_a]=ppipaw->aH;
        // 2: phi
        dy[ppm->index_in_phi]= -ppipaw->dV/3./ppipaw->aH*y[ppm->index_in_a];
        break;
      }

      break;
    }

    break;

  case inflation_H:

    class_call(primordial_inflation_hubble(ppm,
                                           y[ppm->index_in_phi],
                                           &(ppipaw->H),
                                           &(ppipaw->dH),
                                           &(ppipaw->ddH),
                                           &(ppipaw->dddH)),
               ppm->error_message,
               ppm->error_message);

    switch (ppipaw->time) {

    case conformal:

      // 1: a
      dy[ppm->index_in_a]=ppipaw->a2*ppipaw->H;
      // 2: phi
      dy[ppm->index_in_phi]=-1./4./_PI_*y[ppm->index_in_a]*ppipaw->dH;
      break;

    case proper:

      // 1: a
      dy[ppm->index_in_a]=y[ppm->index_in_a]*ppipaw->H;
      // 2: phi
      dy[ppm->index_in_phi]=-1./4./_PI_*ppipaw->dH;
      break;
    }

    // z''/z (assumes that conformal time is requested)
    ppipaw->zpp_over_z =
      2.               *ppipaw->a2*ppipaw->H*ppipaw->H
      -3./4./_PI_      *ppipaw->a2*ppipaw->H*ppipaw->ddH
      +1./16./_PI_/_PI_*ppipaw->a2*ppipaw->ddH*ppipaw->ddH
      +1./16./_PI_/_PI_*ppipaw->a2*ppipaw->dH*ppipaw->dddH
      -1./4./_PI_/_PI_ *ppipaw->a2*ppipaw->dH*ppipaw->dH*ppipaw->ddH/ppipaw->H
      +1./2./_PI_      *ppipaw->a2*ppipaw->dH*ppipaw->dH
      +1./8./_PI_/_PI_ *ppipaw->a2*ppipaw->dH*ppipaw->dH*ppipaw->dH*ppipaw->dH/ppipaw->H/ppipaw->H;

    // a''/a (assumes that conformal time is requested)
    ppipaw->app_over_a = 2.*ppipaw->a2*ppipaw->H*ppipaw->H
      -4.*_PI_*dy[ppm->index_in_phi]*dy[ppm->index_in_phi];

    break;

  default:
    class_stop(ppm->error_message,"ppm->primordial_spec_type=%d different from possible relevant cases",ppm->primordial_spec_type);
    break;

  }

  if (ppipaw->N <= ppm->in_bg_size) // mind the <= instead of ==, necessary because for backward integration 1 equation is removed
    return _SUCCESS_;

  // PERTURBATIONS

  class_test(ppipaw->time == proper,
             ppm->error_message,
             "For inflaton perturbations, only conformal time is coded.");

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
 * @param ppt  Input/output: pointer to perturbation structure
 * @param ppm  Input/output: pointer to primordial structure
 * @return the error status
 */

int primordial_external_spectrum_init(
                                      struct perturbations * ppt,
                                      struct primordial * ppm
                                      ) {
  /** Summary: */

  char arguments[_ARGUMENT_LENGTH_MAX_];
  char line[_LINE_LENGTH_MAX_];
  char command_with_arguments[2*_ARGUMENT_LENGTH_MAX_];
  FILE *process;
  int n_data_guess, n_data = 0;
  double *k = NULL, *pks = NULL, *pkt = NULL;
  double this_k, this_pks, this_pkt;
  int status;
  int index_k;

  /** - Initialization */
  /* Prepare the data (with some initial size) */
  n_data_guess = 100;
  k   = (double *)malloc(n_data_guess*sizeof(double));
  pks = (double *)malloc(n_data_guess*sizeof(double));
  if (ppt->has_tensors == _TRUE_)
    pkt = (double *)malloc(n_data_guess*sizeof(double));
  /* Prepare the command */
  /* If the command is just a "cat", no arguments need to be passed */
  if (strncmp("cat ", ppm->command, 4) == 0) {
    class_sprintf(arguments, " ");
  }
  /* otherwise pass the list of arguments */
  else {
    class_sprintf(arguments, " %g %g %g %g %g %g %g %g %g %g",
            ppm->custom1, ppm->custom2, ppm->custom3, ppm->custom4, ppm->custom5,
            ppm->custom6, ppm->custom7, ppm->custom8, ppm->custom9, ppm->custom10);
  }
  /* write the actual command in a string */
  class_sprintf(command_with_arguments, "%s %s", ppm->command, arguments);
  if (ppm->primordial_verbose > 0)
    printf(" -> running: %s\n",command_with_arguments);

  /** - Launch the command and retrieve the output */
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
    if ((n_data+1) > n_data_guess) {
      n_data_guess *= 2;
      class_realloc(k, n_data_guess*sizeof(double), ppm->error_message);
      class_realloc(pks, n_data_guess*sizeof(double), ppm->error_message);
      if (ppt->has_tensors == _TRUE_) {
        class_realloc(pkt, n_data_guess*sizeof(double), ppm->error_message);
      }
    }
    /* Store */
    k  [n_data]   = this_k;
    pks[n_data]   = this_pks;
    if (ppt->has_tensors == _TRUE_) {
      pkt[n_data] = this_pkt;
    }
    n_data++;
    /* Check ascending order of the k's */
    if (n_data>1) {
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
  class_test(k[1] > ppt->k_min,
             ppm->error_message,
             "Your table for the primordial spectrum does not have "
             "at least 2 points before the minimum value of k: %e . "
             "The splines interpolation would not be safe.",ppt->k_min);
  class_test(k[n_data-2] < ppt->k_max,
             ppm->error_message,
             "Your table for the primordial spectrum does not have "
             "at least 2 points after the maximum value of k: %e . "
             "The splines interpolation would not be safe.",ppt->k_max);

  /** - Store the read results into CLASS structures */
  ppm->lnk_size = n_data;
  /** - Make room */
  class_realloc(ppm->lnk,
                ppm->lnk_size*sizeof(double),
                ppm->error_message);
  class_realloc(ppm->lnpk[ppt->index_md_scalars],
                ppm->lnk_size*sizeof(double),
                ppm->error_message);
  class_realloc(ppm->ddlnpk[ppt->index_md_scalars],
                ppm->lnk_size*sizeof(double),
                ppm->error_message);
  if (ppt->has_tensors == _TRUE_) {
    class_realloc(ppm->lnpk[ppt->index_md_tensors],
                  ppm->lnk_size*sizeof(double),
                  ppm->error_message);
    class_realloc(ppm->ddlnpk[ppt->index_md_tensors],
                  ppm->lnk_size*sizeof(double),
                  ppm->error_message);
  };
  /** - Store values */
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
  /** - Release the memory used locally */
  free(k);
  free(pks);
  if (ppt->has_tensors == _TRUE_)
    free(pkt);
  /** - Tell CLASS that there are scalar (and tensor) modes */
  ppm->is_non_zero[ppt->index_md_scalars][ppt->index_ic_ad] = _TRUE_;
  if (ppt->has_tensors == _TRUE_)
    ppm->is_non_zero[ppt->index_md_tensors][ppt->index_ic_ten] = _TRUE_;

  return _SUCCESS_;
}

int primordial_output_titles(struct perturbations * ppt,
                             struct primordial * ppm,
                             char titles[_MAXTITLESTRINGLENGTH_]
                             ){
  class_store_columntitle(titles,"k [1/Mpc]",_TRUE_);
  class_store_columntitle(titles,"P_scalar(k)",_TRUE_);
  class_store_columntitle(titles,"P_tensor(k)",ppt->has_tensors);

  return _SUCCESS_;

}

int primordial_output_data(struct perturbations * ppt,
                           struct primordial * ppm,
                           int number_of_titles,
                           double *data){

  int index_k, storeidx;
  double *dataptr;

  for (index_k=0; index_k<ppm->lnk_size; index_k++) {
    dataptr = data + index_k*number_of_titles;
    storeidx = 0;

    class_store_double(dataptr, exp(ppm->lnk[index_k]), _TRUE_,storeidx);
    class_store_double(dataptr, exp(ppm->lnpk[ppt->index_md_scalars][index_k]), _TRUE_,storeidx);
    class_store_double(dataptr, exp(ppm->lnpk[ppt->index_md_tensors][index_k]), ppt->has_tensors,storeidx);
  }


  return _SUCCESS_;

}
