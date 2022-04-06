/** @file fourier.c Documented fourier module
 *
 * Julien Lesgourgues, 6.03.2014
 *
 * New module replacing an older one present up to version 2.0 The new
 * module is located in a better place in the main, allowing it to
 * compute non-linear correction to \f$ C_l\f$'s and not just \f$ P(k)\f$. It will
 * also be easier to generalize to new methods.  The old implementation
 * of one-loop calculations and TRG calculations has been dropped from
 * this version, they can still be found in older versions.
 *
 */

#include "fourier.h"

/**
 * Return the P(k,z) for a given redshift z and pk type (_m, _cb)
 * (linear if pk_output = pk_linear, nonlinear if pk_output = pk_nonlinear)
 *
 * In the linear case, if there are several initial conditions *and* the
 * input pointer out_pk_ic is not set to NULL, the function also
 * returns the decomposition into different IC contributions.
 *
 * Hints on input index_pk:
 *
 * a. if you want the total matter spectrum P_m(k,z), pass in input
 *    pfo->index_pk_total
 *    (this index is always defined)
 *
 * b. if you want the power spectrum relevant for galaxy or halos,
 *    given by P_cb if there is non-cold-dark-matter (e.g. massive neutrinos)
 *    and to P_m otherwise, pass in input
 *    pfo->index_pk_cluster
 *    (this index is always defined)
 *
 * c. there is another possible syntax (use it only if you know what you are doing):
 *    if pfo->has_pk_m == _TRUE_ you may pass pfo->index_pk_m to get P_m
 *    if pfo->has_pk_cb == _TRUE_ you may pass pfo->index_pk_cb to get P_cb
 *
 * Output format:
 *
 * 1. if mode = logarithmic (most straightforward for the code):
 *     out_pk = ln(P(k))
 *     out_pk_ic[diagonal] = ln(P_ic(k))
 *     out_pk_ic[non-diagonal] = cos(correlation angle icxic)
 *
 * 2. if mode = linear (a conversion is done internally in this function)
 *     out_pk = P(k)
 *     out_pk_ic[diagonal] = P_ic(k)
 *     out_pk_ic[non-diagonal] = P_icxic(k)
 *
 * @param pba         Input: pointer to background structure
 * @param pfo         Input: pointer to fourier structure
 * @param mode        Input: linear or logarithmic
 * @param pk_output   Input: linear or nonlinear
 * @param z           Input: redshift
 * @param index_pk    Input: index of pk type (_m, _cb)
 * @param out_pk      Output: P(k) returned as out_pk_l[index_k]
 * @param out_pk_ic   Output:  P_ic(k) returned as  out_pk_ic[index_k * pfo->ic_ic_size + index_ic1_ic2]
 * @return the error status
 */

int fourier_pk_at_z(
                    struct background * pba,
                    struct fourier *pfo,
                    enum linear_or_logarithmic mode,
                    enum pk_outputs pk_output,
                    double z,
                    int index_pk,
                    double * out_pk, // array out_pk[index_k]
                    double * out_pk_ic // array out_pk_ic[index_k * pfo->ic_ic_size + index_ic1_ic2]
                    ) {
  double tau;
  double ln_tau;
  int index_k;
  int index_ic1;
  int index_ic2;
  int index_ic1_ic1;
  int index_ic2_ic2;
  int index_ic1_ic2;
  int last_index;
  short do_ic = _FALSE_;

  /** - check whether we need the decomposition into contributions from each initial condition */

  if ((pk_output == pk_linear) && (pfo->ic_size > 1) && (out_pk_ic != NULL))
    do_ic = _TRUE_;

  /** - case z=0 requiring no interpolation in z */
  if (z == 0) {

    for (index_k=0; index_k<pfo->k_size; index_k++) {

      if (pk_output == pk_linear) {
        out_pk[index_k] = pfo->ln_pk_l[index_pk][(pfo->ln_tau_size-1)*pfo->k_size+index_k];

        if (do_ic == _TRUE_) {
          for (index_ic1_ic2 = 0; index_ic1_ic2 < pfo->ic_ic_size; index_ic1_ic2++) {
            out_pk_ic[index_k * pfo->ic_ic_size + index_ic1_ic2] =
              pfo->ln_pk_ic_l[index_pk][((pfo->ln_tau_size-1)*pfo->k_size+index_k)*pfo->ic_ic_size+index_ic1_ic2];
          }
        }
      }
      else {
        out_pk[index_k] = pfo->ln_pk_nl[index_pk][(pfo->ln_tau_size-1)*pfo->k_size+index_k];
      }
    }
  }

  /** - interpolation in z */
  else {

    class_test(pfo->ln_tau_size == 1,
               pfo->error_message,
               "You are asking for the matter power spectrum at z=%e but the code was asked to store it only at z=0. You probably forgot to pass the input parameter z_max_pk (see explanatory.ini)",z);

    /** --> get value of contormal time tau */
    class_call(background_tau_of_z(pba,
                                   z,
                                   &tau),
               pba->error_message,
               pfo->error_message);

    ln_tau = log(tau);
    last_index = pfo->ln_tau_size-1;

    /** -> check that tau is in pre-computed table */

    if (ln_tau <= pfo->ln_tau[0]) {

      /** --> if ln(tau) much too small, raise an error */
      class_test(ln_tau<pfo->ln_tau[0]-_EPSILON_,
                 pfo->error_message,
                 "requested z was not inside of tau tabulation range (Requested ln(tau_=%.10e, Min %.10e). Solution might be to increase input parameter z_max_pk (see explanatory.ini)",ln_tau,pfo->ln_tau[0]);

      /** --> if ln(tau) too small but within tolerance, round it and get right values without interpolating */
      ln_tau = pfo->ln_tau[0];

      for (index_k = 0 ; index_k < pfo->k_size; index_k++) {
        if (pk_output == pk_linear) {
          out_pk[index_k] = pfo->ln_pk_l[index_pk][index_k];
          if (do_ic == _TRUE_) {
            for (index_ic1_ic2 = 0; index_ic1_ic2 < pfo->ic_ic_size; index_ic1_ic2++) {
              out_pk_ic[index_k * pfo->ic_ic_size + index_ic1_ic2] = pfo->ln_pk_ic_l[index_pk][index_k * pfo->ic_ic_size + index_ic1_ic2];
            }
          }
        }
        else {
          out_pk[index_k] = pfo->ln_pk_nl[index_pk][index_k];
        }
      }
    }

    else if (ln_tau >= pfo->ln_tau[pfo->ln_tau_size-1]) {

      /** --> if ln(tau) much too large, raise an error */
      class_test(ln_tau>pfo->ln_tau[pfo->ln_tau_size-1]+_EPSILON_,
                 pfo->error_message,
                 "requested z was not inside of tau tabulation range (Requested ln(tau_=%.10e, Max %.10e) ",ln_tau,pfo->ln_tau[pfo->ln_tau_size-1]);

      /** --> if ln(tau) too large but within tolerance, round it and get right values without interpolating */
      ln_tau = pfo->ln_tau[pfo->ln_tau_size-1];

      for (index_k = 0 ; index_k < pfo->k_size; index_k++) {
        if (pk_output == pk_linear) {
          out_pk[index_k] = pfo->ln_pk_l[index_pk][(pfo->ln_tau_size-1) * pfo->k_size + index_k];
          if (do_ic == _TRUE_) {
            for (index_ic1_ic2 = 0; index_ic1_ic2 < pfo->ic_ic_size; index_ic1_ic2++) {
              out_pk_ic[index_k * pfo->ic_ic_size + index_ic1_ic2] = pfo->ln_pk_ic_l[index_pk][((pfo->ln_tau_size-1) * pfo->k_size + index_k) * pfo->ic_ic_size + index_ic1_ic2];
            }
          }
        }
        else {
          out_pk[index_k] = pfo->ln_pk_nl[index_pk][(pfo->ln_tau_size-1) * pfo->k_size + index_k];
        }
      }
    }

    /** -> tau is in pre-computed table: interpolate */
    else {

      if (pk_output == pk_linear) {

        /** --> interpolate P_l(k) at tau from pre-computed array */
        class_call(array_interpolate_spline(pfo->ln_tau,
                                            pfo->ln_tau_size,
                                            pfo->ln_pk_l[index_pk],
                                            pfo->ddln_pk_l[index_pk],
                                            pfo->k_size,
                                            ln_tau,
                                            &last_index,
                                            out_pk,
                                            pfo->k_size,
                                            pfo->error_message),
                   pfo->error_message,
                   pfo->error_message);

        /** --> interpolate P_ic_l(k) at tau from pre-computed array */
        if (do_ic == _TRUE_) {
          class_call(array_interpolate_spline(pfo->ln_tau,
                                              pfo->ln_tau_size,
                                              pfo->ln_pk_ic_l[index_pk],
                                              pfo->ddln_pk_ic_l[index_pk],
                                              pfo->k_size*pfo->ic_ic_size,
                                              ln_tau,
                                              &last_index,
                                              out_pk_ic,
                                              pfo->k_size*pfo->ic_ic_size,
                                              pfo->error_message),
                     pfo->error_message,
                     pfo->error_message);
        }
      }
      else {

        /** --> interpolate P_nl(k) at tau from pre-computed array */
        class_call(array_interpolate_spline(pfo->ln_tau,
                                            pfo->ln_tau_size,
                                            pfo->ln_pk_nl[index_pk],
                                            pfo->ddln_pk_nl[index_pk],
                                            pfo->k_size,
                                            ln_tau,
                                            &last_index,
                                            out_pk,
                                            pfo->k_size,
                                            pfo->error_message),
                   pfo->error_message,
                   pfo->error_message);
      }
    }
  }

  /** - so far, all output stored in logarithmic format. Eventually, convert to linear one. */

  if (mode == linear) {

    /** --> loop over k */
    for (index_k=0; index_k<pfo->k_size; index_k++) {

      /** --> convert total spectrum */
      out_pk[index_k] = exp(out_pk[index_k]);

      if (do_ic == _TRUE_) {
        /** --> convert contribution of each ic (diagonal elements) */
        for (index_ic1=0; index_ic1 < pfo->ic_size; index_ic1++) {
          index_ic1_ic1 = index_symmetric_matrix(index_ic1,index_ic1,pfo->ic_size);

          out_pk_ic[index_k * pfo->ic_ic_size + index_ic1_ic1] = exp(out_pk_ic[index_k * pfo->ic_ic_size + index_ic1_ic1]);
        }

        /** --> convert contribution of each ic (non-diagonal elements) */
        for (index_ic1=0; index_ic1 < pfo->ic_size; index_ic1++) {
          for (index_ic2=index_ic1+1; index_ic2 < pfo->ic_size; index_ic2++) {
            index_ic1_ic1 = index_symmetric_matrix(index_ic1,index_ic1,pfo->ic_size);
            index_ic2_ic2 = index_symmetric_matrix(index_ic2,index_ic2,pfo->ic_size);
            index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,pfo->ic_size);

            /* P_ic1xic2 = cos(angle) * sqrt(P_ic1 * P_ic2) */
            out_pk_ic[index_k * pfo->ic_ic_size + index_ic1_ic2]
              = out_pk_ic[index_k * pfo->ic_ic_size + index_ic1_ic2]
              *sqrt(out_pk_ic[index_k * pfo->ic_ic_size + index_ic1_ic1]
                    *out_pk_ic[index_k * pfo->ic_ic_size + index_ic2_ic2]);
          }
        }
      }
    }
  }

  return _SUCCESS_;
}

/*
 * Same as fourier_pk_at_z() (see the comments there about
 * the input/output format), excepted that we don't pass in input one
 * type of P(k) through index_pk. Instead, we get all existing types
 * in output. This function is maintained to enhance compatibility
 * with old versions, but the use of fourier_pk_at_z() should
 * be preferred.
 *
 * @param pba            Input: pointer to background structure
 * @param pfo            Input: pointer to fourier structure
 * @param mode           Input: linear or logarithmic
 * @param z              Input: redshift
 * @param out_pk_l       Output: P_m(k) returned as out_pk_l[index_k]
 * @param out_pk_ic_l    Output:  P_m_ic(k) returned as  out_pk_ic_l[index_k * pfo->ic_ic_size + index_ic1_ic2]
 * @param out_pk_cb_l    Output: P_cb(k) returned as out_pk_cb_l[index_k]
 * @param out_pk_cb_ic_l Output:  P_cb_ic(k) returned as  out_pk_cb_ic_l[index_k * pfo->ic_ic_size + index_ic1_ic2]
 * @return the error status
 */

int fourier_pks_at_z(
                     struct background * pba,
                     struct fourier * pfo,
                     enum linear_or_logarithmic mode,
                     enum pk_outputs pk_output,
                     double z,
                     double * out_pk,      // array out_pk[index_k]
                     double * out_pk_ic,   // array out_pk_ic[index_k * pfo->ic_ic_size + index_ic1_ic2]
                     double * out_pk_cb,   // array out_pk_cb[index_k]
                     double * out_pk_cb_ic // array out_pk_cb_ic[index_k * pfo->ic_ic_size + index_ic1_ic2]
                     ) {

  if (pfo->has_pk_cb == _TRUE_) {

    class_call(fourier_pk_at_z(pba,
                               pfo,
                               mode,
                               pk_output,
                               z,
                               pfo->index_pk_cb,
                               out_pk_cb,
                               out_pk_cb_ic
                               ),
               pfo->error_message,
               pfo->error_message);
  }

  if (pfo->has_pk_m == _TRUE_) {

    class_call(fourier_pk_at_z(pba,
                               pfo,
                               mode,
                               pk_output,
                               z,
                               pfo->index_pk_m,
                               out_pk,
                               out_pk_ic
                               ),
               pfo->error_message,
               pfo->error_message);
  }

  return _SUCCESS_;
}

/**
 * Return the P(k,z) for a given (k,z) and pk type (_m, _cb)
 * (linear if pk_output = pk_linear, nonlinear if pk_output = pk_nonlinear)
 *
 * In the linear case, if there are several initial conditions *and* the
 * input pointer out_pk_ic is not set to NULL, the function also
 * returns the decomposition into different IC contributions.
 *
 * Hints on input index_pk:
 *
 * a. if you want the total matter spectrum P_m(k,z), pass in input
 *    pfo->index_pk_total
 *    (this index is always defined)
 *
 * b. if you want the power spectrum relevant for galaxy or halos,
 *    given by P_cb if there is non-cold-dark-matter (e.g. massive neutrinos)
 *    and to P_m otherwise, pass in input
 *    pfo->index_pk_cluster
 *    (this index is always defined)
 *
 * c. there is another possible syntax (use it only if you know what you are doing):
 *    if pfo->has_pk_m == _TRUE_ you may pass pfo->index_pk_m to get P_m
 *    if pfo->has_pk_cb == _TRUE_ you may pass pfo->index_pk_cb to get P_cb
 *
 * Output format:
 *
 *     out_pk = P(k)
 *     out_pk_ic[diagonal] = P_ic(k)
 *     out_pk_ic[non-diagonal] = P_icxic(k)
 *
 * @param pba         Input: pointer to background structure
 * @param ppm         Input: pointer to primordial structure
 * @param pfo         Input: pointer to fourier structure
 * @param pk_output   Input: linear or nonlinear
 * @param k           Input: wavenumber in 1/Mpc
 * @param z           Input: redshift
 * @param index_pk    Input: index of pk type (_m, _cb)
 * @param out_pk      Output: pointer to P
 * @param out_pk_ic   Ouput:  P_ic returned as out_pk_ic_l[index_ic1_ic2]
 * @return the error status
 */

int fourier_pk_at_k_and_z(
                          struct background * pba,
                          struct primordial * ppm,
                          struct fourier *pfo,
                          enum pk_outputs pk_output,
                          double k,
                          double z,
                          int index_pk,
                          double * out_pk, // number *out_pk_l
                          double * out_pk_ic // array out_pk_ic_l[index_ic_ic]
                          ) {

  double * out_pk_at_z;
  double * out_pk_ic_at_z = NULL;
  double * ddout_pk_at_z;
  double * ddout_pk_ic_at_z;
  int last_index;
  int index_ic1;
  int index_ic2;
  int index_ic1_ic1;
  int index_ic2_ic2;
  int index_ic1_ic2;
  double kmin;
  double * pk_primordial_k;
  double * pk_primordial_kmin;
  short do_ic = _FALSE_;

  /** - preliminary: check whether we need the decomposition into contributions from each initial condition */

  if ((pk_output == pk_linear) && (pfo->ic_size > 1) && (out_pk_ic != NULL))
    do_ic = _TRUE_;

  /** - first step: check that k is in valid range [0:kmax]
      (the test for z will be done when calling fourier_pk_linear_at_z()) */

  class_test((k < 0.) || (k > exp(pfo->ln_k[pfo->k_size-1])),
             pfo->error_message,
             "k=%e out of bounds [%e:%e]",k,0.,exp(pfo->ln_k[pfo->k_size-1]));

  /** - deal with case k = 0 for which P(k) is set to zero
      (this non-physical result can be useful for interpolations) */

  if (k == 0.) {
    *out_pk = 0.;

    if (do_ic == _TRUE_) {
      for (index_ic1_ic2=0; index_ic1_ic2<pfo->ic_ic_size; index_ic1_ic2++) {
        out_pk_ic[index_ic1_ic2] = 0.;
      }
    }
  }

  /** - deal with 0 < k <= kmax */

  else {

    class_alloc(out_pk_at_z,
                pfo->k_size*sizeof(double),
                pfo->error_message);

    if (do_ic == _TRUE_) {
      class_alloc(out_pk_ic_at_z,
                  pfo->k_size*pfo->ic_ic_size*sizeof(double),
                  pfo->error_message);
    }

    /** - deal with standard case kmin <= k <= kmax */

    if (k > exp(pfo->ln_k[0])) {

      /** --> First, get P(k) at the right z (in logarithmic format for more accurate interpolation, and then convert to linear format) */

      class_call(fourier_pk_at_z(pba,
                                 pfo,
                                 logarithmic,
                                 pk_output,
                                 z,
                                 index_pk,
                                 out_pk_at_z,
                                 out_pk_ic_at_z
                                 ),
                 pfo->error_message,
                 pfo->error_message);

      /** --> interpolate total spectrum */

      class_alloc(ddout_pk_at_z,
                  pfo->k_size*sizeof(double),
                  pfo->error_message);

      class_call(array_spline_table_lines(pfo->ln_k,
                                          pfo->k_size,
                                          out_pk_at_z,
                                          1,
                                          ddout_pk_at_z,
                                          _SPLINE_NATURAL_,
                                          pfo->error_message),
                 pfo->error_message,
                 pfo->error_message);

      class_call(array_interpolate_spline(pfo->ln_k,
                                          pfo->k_size,
                                          out_pk_at_z,
                                          ddout_pk_at_z,
                                          1,
                                          log(k),
                                          &last_index,
                                          out_pk,
                                          1,
                                          pfo->error_message),
                 pfo->error_message,
                 pfo->error_message);

      free(ddout_pk_at_z);

      // uncomment this part if you prefer a linear interpolation

      /*
        class_call(array_interpolate_linear(pfo->ln_k,
        pfo->k_size,
        out_pk_at_z,
        1,
        log(k),
        &last_index,
        out_pk,
        1,
        pfo->error_message),
        pfo->error_message,
        pfo->error_message);
      */

      /** --> convert from logarithmic to linear format */

      *out_pk = exp(*out_pk);

      /** --> interpolate each ic component */

      if (do_ic == _TRUE_) {

        class_alloc(ddout_pk_ic_at_z,
                    pfo->k_size*pfo->ic_ic_size*sizeof(double),
                    pfo->error_message);

        class_call(array_spline_table_lines(pfo->ln_k,
                                            pfo->k_size,
                                            out_pk_ic_at_z,
                                            pfo->ic_ic_size,
                                            ddout_pk_ic_at_z,
                                            _SPLINE_NATURAL_,
                                            pfo->error_message),
                   pfo->error_message,
                   pfo->error_message);

        class_call(array_interpolate_spline(pfo->ln_k,
                                            pfo->k_size,
                                            out_pk_ic_at_z,
                                            ddout_pk_ic_at_z,
                                            pfo->ic_ic_size,
                                            log(k),
                                            &last_index,
                                            out_pk_ic,
                                            pfo->ic_ic_size,
                                            pfo->error_message),
                   pfo->error_message,
                   pfo->error_message);

        free(ddout_pk_ic_at_z);

        /** --> convert each ic component from logarithmic to linear format */

        for (index_ic1=0; index_ic1 < pfo->ic_size; index_ic1++) {
          index_ic1_ic1 = index_symmetric_matrix(index_ic1,index_ic1,pfo->ic_size);
          out_pk_ic[index_ic1_ic1] = exp(out_pk_ic[index_ic1_ic1]);
        }
        for (index_ic1=0; index_ic1 < pfo->ic_size; index_ic1++) {
          for (index_ic2=index_ic1+1; index_ic2 < pfo->ic_size; index_ic2++) {
            index_ic1_ic1 = index_symmetric_matrix(index_ic1,index_ic1,pfo->ic_size);
            index_ic2_ic2 = index_symmetric_matrix(index_ic2,index_ic2,pfo->ic_size);
            index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,pfo->ic_size);
            out_pk_ic[index_ic1_ic2]
              = out_pk_ic[index_ic1_ic2]*sqrt(out_pk_ic[index_ic1_ic1]*out_pk_ic[index_ic2_ic2]);
          }
        }
      }
    }

    /** --> deal with case 0 < k < kmin that requires extrapolation
     *    P(k) = [some number] * k  * P_primordial(k)
     *    so
     *    P(k) = P(kmin) * (k P_primordial(k)) / (kmin P_primordial(kmin))
     *    (note that the result is accurate only if kmin is such that [a0 kmin] << H0)
     *
     *    This is accurate for the synchronous gauge; TODO: write
     *    newtonian gauge case. Also, In presence of isocurvature
     *    modes, we assumes for simplicity that the mode with
     *    index_ic1_ic2=0 dominates at small k: exact treatment should be
     *    written if needed.
     */

    else {

      /** --> First, get P(k) at the right z (in linear format) */

      class_call(fourier_pk_at_z(pba,
                                 pfo,
                                 linear,
                                 pk_output,
                                 z,
                                 index_pk,
                                 out_pk_at_z,
                                 out_pk_ic_at_z
                                 ),
                 pfo->error_message,
                 pfo->error_message);

      /* get P(kmin) */
      *out_pk = out_pk_at_z[0];

      if (do_ic == _TRUE_) {
        for (index_ic1_ic2=0; index_ic1_ic2<pfo->ic_ic_size; index_ic1_ic2++) {
          out_pk_ic[index_ic1_ic2] = out_pk_ic_at_z[index_ic1_ic2];
        }
      }

      /* compute P_primordial(k) */

      class_alloc(pk_primordial_k,
                  sizeof(double)*pfo->ic_ic_size,
                  pfo->error_message);

      class_call(primordial_spectrum_at_k(ppm,
                                          pfo->index_md_scalars,
                                          linear,
                                          k,
                                          pk_primordial_k),
                 ppm->error_message,
                 pfo->error_message);

      /* compute P_primordial(kmin) */

      kmin = exp(pfo->ln_k[0]);

      class_alloc(pk_primordial_kmin,
                  sizeof(double)*pfo->ic_ic_size,
                  pfo->error_message);

      class_call(primordial_spectrum_at_k(ppm,
                                          pfo->index_md_scalars,
                                          linear,
                                          kmin,
                                          pk_primordial_kmin),
                 ppm->error_message,
                 pfo->error_message);

      /* finally, infer P(k) */

      *out_pk *= (k*pk_primordial_k[0]/kmin/pk_primordial_kmin[0]);

      if (do_ic == _TRUE_) {
        for (index_ic1_ic2=0; index_ic1_ic2<pfo->ic_ic_size; index_ic1_ic2++) {
          out_pk_ic[index_ic1_ic2] *= (k*pk_primordial_k[index_ic1_ic2]
                                       /kmin/pk_primordial_kmin[index_ic1_ic2]);
        }
      }

      free(pk_primordial_k);
      free(pk_primordial_kmin);
    }

    free(out_pk_at_z);
    if (do_ic == _TRUE_) {
      free(out_pk_ic_at_z);
    }
  }

  return _SUCCESS_;
}

/*
 * Same as fourier_pk_at_k_and_z() (see the comments there about
 * the input/output format), excepted that we don't pass in input one
 * type of P(k) through index_pk. Instead, we get all existing types
 * in output. This function is maintained to enhance compatibility
 * with old versions, but the use of fourier_pk_at_k_and_z() should
 * be preferred.
 *
 * @param pba            Input: pointer to background structure
 * @param ppm            Input: pointer to primordial structure
 * @param pfo            Input: pointer to fourier structure
 * @param k              Input: wavenumber
 * @param z              Input: redshift
 * @param out_pk_l       Output: P_m(k) returned as out_pk_l[index_k]
 * @param out_pk_ic_l    Output:  P_m_ic(k) returned as  out_pk_ic_l[index_k * pfo->ic_ic_size + index_ic1_ic2]
 * @param out_pk_cb_l    Output: P_cb(k) returned as out_pk_cb_l[index_k]
 * @param out_pk_cb_ic_l Output:  P_cb_ic(k) returned as  out_pk_cb_ic_l[index_k * pfo->ic_ic_size + index_ic1_ic2]
 * @return the error status
 */

int fourier_pks_at_k_and_z(
                           struct background * pba,
                           struct primordial * ppm,
                           struct fourier *pfo,
                           enum pk_outputs pk_output,
                           double k,
                           double z,
                           double * out_pk, // number P_m(k)
                           double * out_pk_ic, // array P_m_ic(k) of index [index_ic1_ic2]
                           double * out_pk_cb, // number P_cb(k)
                           double * out_pk_cb_ic // array P__cb_ic(k)of index [index_ic1_ic2]
                           ) {

  if (pfo->has_pk_cb == _TRUE_) {

    class_call(fourier_pk_at_k_and_z(pba,
                                     ppm,
                                     pfo,
                                     pk_output,
                                     k,
                                     z,
                                     pfo->index_pk_cb,
                                     out_pk_cb,
                                     out_pk_cb_ic
                                     ),
               pfo->error_message,
               pfo->error_message);
  }
  if (pfo->has_pk_m == _TRUE_) {

    class_call(fourier_pk_at_k_and_z(pba,
                                     ppm,
                                     pfo,
                                     pk_output,
                                     k,
                                     z,
                                     pfo->index_pk_m,
                                     out_pk,
                                     out_pk_ic
                                     ),
               pfo->error_message,
               pfo->error_message);
  }

  return _SUCCESS_;
}

/**
 * Return the P(k,z) for a grid of (k_i,z_j) passed in input,
 * for all available pk types (_m, _cb),
 * either linear or nonlinear depending on input.
 *
 * If there are several initial conditions, this function is not
 * designed to return individual contributions.
 *
 * The main goal of this routine is speed. Unlike
 * fourier_pk_at_k_and_z(), it performs no extrapolation when an
 * input k_i falls outside the pre-computed range [kmin,kmax]: in that
 * case, it just returns P(k,z)=0 for such a k_i
 *
 * @param pba            Input: pointer to background structure
 * @param pfo            Input: pointer to fourier structure
 * @param pk_output      Input: pk_linear or pk_nonlinear
 * @param kvec           Input: array of wavenumbers in ascending order (in 1/Mpc)
 * @param kvec_size      Input: size of array of wavenumbers
 * @param zvec           Input: array of redshifts in arbitrary order
 * @param zvec_size      Input: size of array of redshifts
 * @param out_pk         Output: P(k_i,z_j) for total matter (if available) in Mpc**3
 * @param out_pk_cb      Output: P_cb(k_i,z_j) for cdm+baryons (if available) in Mpc**3
 * @return the error status
 */

int fourier_pks_at_kvec_and_zvec(
                                 struct background * pba,
                                 struct fourier * pfo,
                                 enum pk_outputs pk_output,
                                 double * kvec, // kvec[index_kvec]
                                 int kvec_size,
                                 double * zvec, // zvec[index_zvec]
                                 int zvec_size,
                                 double * out_pk,   // output_pk[index_zvec*kvec_size+index_kvec],
                                                    // already allocated
                                                    //(or NULL if user knows there is no _m output)
                                 double * out_pk_cb // output_pk[index_zvec*kvec_size+index_kvec],
                                                    //already allocated
                                                    //(or NULL if user knows there is no _cb output)
                                 ) {

  /** Summary: */

  /** - define local variables */

  int index_k, index_kvec, index_zvec;
  double * ln_kvec;
  double * ln_pk_table = NULL;
  double * ddln_pk_table = NULL;
  double * ln_pk_cb_table = NULL;
  double * ddln_pk_cb_table = NULL;
  double h, a, b;

  /** - Allocate arrays */

  class_alloc(ln_kvec, sizeof(double)*kvec_size,
              pfo->error_message);

  if (pfo->has_pk_m == _TRUE_) {
    class_alloc(ln_pk_table, sizeof(double)*pfo->k_size*zvec_size,
                pfo->error_message);
    class_alloc(ddln_pk_table, sizeof(double)*pfo->k_size*zvec_size,
                pfo->error_message);
  }
  if (pfo->has_pk_cb == _TRUE_) {
    class_alloc(ln_pk_cb_table, sizeof(double)*pfo->k_size*zvec_size,
                pfo->error_message);
    class_alloc(ddln_pk_cb_table, sizeof(double)*pfo->k_size*zvec_size,
                pfo->error_message);
  }

  /** - Construct table of log(P(k_n,z_j)) for pre-computed wavenumbers but requested redshifts: */

  for (index_zvec=0; index_zvec<zvec_size; index_zvec++){

    if (pfo->has_pk_m == _TRUE_) {
      class_call(fourier_pk_at_z(pba,
                                 pfo,
                                 logarithmic,
                                 pk_output,
                                 zvec[index_zvec],
                                 pfo->index_pk_m,
                                 &(ln_pk_table[index_zvec * pfo->k_size]),
                                 NULL),
                 pfo->error_message,
                 pfo->error_message);
    }
    if (pfo->has_pk_cb == _TRUE_) {
      class_call(fourier_pk_at_z(pba,
                                 pfo,
                                 logarithmic,
                                 pk_output,
                                 zvec[index_zvec],
                                 pfo->index_pk_cb,
                                 &(ln_pk_cb_table[index_zvec * pfo->k_size]),
                                 NULL),
                 pfo->error_message,
                 pfo->error_message);
    }
  }

  /** - Spline it for interpolation along k */

  if (pfo->has_pk_m == _TRUE_) {

    class_call(array_spline_table_columns2(pfo->ln_k,
                                           pfo->k_size,
                                           ln_pk_table,
                                           zvec_size,
                                           ddln_pk_table,
                                           _SPLINE_NATURAL_,
                                           pfo->error_message),
               pfo->error_message,
               pfo->error_message);
  }
  if (pfo->has_pk_cb == _TRUE_) {

    class_call(array_spline_table_columns2(pfo->ln_k,
                                           pfo->k_size,
                                           ln_pk_cb_table,
                                           zvec_size,
                                           ddln_pk_cb_table,
                                           _SPLINE_NATURAL_,
                                           pfo->error_message),
               pfo->error_message,
               pfo->error_message);
  }

  /** - Construct ln(kvec): */

  for (index_kvec=0; index_kvec<kvec_size; index_kvec++)
    ln_kvec[index_kvec] = log(kvec[index_kvec]);

  /** - Loop over first k values. If k<kmin, fill output with zeros. If not, go to next step. */

  for (index_kvec = 0; index_kvec<kvec_size; index_kvec++){

    /* check whether one should go to next step */
    if (ln_kvec[index_kvec] >= pfo->ln_k[0])
      break;

    /* deal with k<k_min */
    for (index_zvec = 0; index_zvec < zvec_size; index_zvec++) {
      if (pfo->has_pk_m == _TRUE_)     out_pk[index_zvec*kvec_size+index_kvec] = 0.;
      if (pfo->has_pk_cb == _TRUE_) out_pk_cb[index_zvec*kvec_size+index_kvec] = 0.;
      /* (If needed, one could add instead some extrapolation here) */
    }
  }

  /** - Deal with case kmin<=k<=kmax. For better performance, do not
      loop through kvec, but through pre-computed k values. */

  for (index_k=0; index_k < (pfo->k_size-1); index_k++){

    /** --> Loop through k_i's that fall in interval [k_n,k_n+1] */

    while ((index_kvec < kvec_size) && (ln_kvec[index_kvec] <= pfo->ln_k[index_k+1])){

      /** --> for each of them, perform spine interpolation */

      h = pfo->ln_k[index_k+1]-pfo->ln_k[index_k];
      b = (ln_kvec[index_kvec] - pfo->ln_k[index_k])/h;
      a = 1.-b;

      for (index_zvec = 0; index_zvec < zvec_size; index_zvec++) {

        if (pfo->has_pk_m == _TRUE_) {

          out_pk[index_zvec*kvec_size+index_kvec] =
            exp(
                array_spline_eval(ln_pk_table,
                                  ddln_pk_table,
                                  (index_zvec * pfo->k_size + index_k),
                                  (index_zvec * pfo->k_size + index_k+1),
                                  h,a,b)
                );
        }
        if (pfo->has_pk_cb == _TRUE_) {

          out_pk_cb[index_zvec*kvec_size+index_kvec] =
            exp(
                array_spline_eval(ln_pk_cb_table,
                                  ddln_pk_cb_table,
                                  (index_zvec * pfo->k_size + index_k),
                                  (index_zvec * pfo->k_size + index_k+1),
                                  h,a,b)
                );
        }
      }
      index_kvec++;
    }
  }

  /** - Loop over possible remaining k values with k > kmax, to fill output with zeros. */

  while (index_kvec < kvec_size) {

    for (index_zvec = 0; index_zvec < zvec_size; index_zvec++) {
      if (pfo->has_pk_m == _TRUE_)     out_pk[index_zvec*kvec_size+index_kvec] = 0.;
      if (pfo->has_pk_cb == _TRUE_) out_pk_cb[index_zvec*kvec_size+index_kvec] = 0.;
      /* (If needed, one could add instead some extrapolation here) */
    }
    index_kvec++;
  }

  free(ln_kvec);
  if (pfo->has_pk_m == _TRUE_) {
    free(ln_pk_table);
    free(ddln_pk_table);
  }
  if (pfo->has_pk_cb == _TRUE_) {
    free(ln_pk_cb_table);
    free(ddln_pk_cb_table);
  }

  return _SUCCESS_;
}

/**
 * Return the logarithmic slope of P(k,z) for a given (k,z), a given pk type (_m, _cb)
 * (computed with linear P_L if pk_output = pk_linear, nonlinear P_NL if pk_output = pk_nonlinear)
 *
 * @param pba         Input: pointer to background structure
 * @param ppm         Input: pointer to primordial structure
 * @param pfo         Input: pointer to fourier structure
 * @param pk_output   Input: linear or nonlinear
 * @param k           Input: wavenumber in 1/Mpc
 * @param z           Input: redshift
 * @param index_pk    Input: index of pk type (_m, _cb)
 * @param pk_tilt     Output: logarithmic slope of P(k,z)
 * @return the error status
 */

int fourier_pk_tilt_at_k_and_z(
                               struct background * pba,
                               struct primordial * ppm,
                               struct fourier * pfo,
                               enum pk_outputs pk_output,
                               double k,
                               double z,
                               int index_pk,
                               double * pk_tilt
                               ) {

  double dlnk;
  double out_pk1,out_pk2;

  /* typical step dln(k) on which we believe that out results are not
     dominated by numerical errors and that the P(k,z) is slowly
     varying */

  dlnk = pfo->ln_k[pfo->k_size-1] - pfo->ln_k[pfo->k_size-2];

  class_call(fourier_pk_at_k_and_z(pba,
                                   ppm,
                                   pfo,
                                   pk_output,
                                   k/(1.+dlnk),
                                   z,
                                   index_pk,
                                   &out_pk1,
                                   NULL),
             pfo->error_message,
             pfo->error_message);

  class_call(fourier_pk_at_k_and_z(pba,
                                   ppm,
                                   pfo,
                                   pk_output,
                                   k*(1.+dlnk),
                                   z,
                                   index_pk,
                                   &out_pk2,
                                   NULL),
             pfo->error_message,
             pfo->error_message);

  /* logarithmic derivative: n_eff = (logPk2 - logPk1)/(logk2-logk1) */

  *pk_tilt = (log(out_pk2)-log(out_pk1))/(2.*log(1.+dlnk));

  return _SUCCESS_;

}

/**
 * This routine computes the variance of density fluctuations in a
 * sphere of radius R at redshift z, sigma(R,z), or other similar derived
 * quantitites, for one given pk type (_m, _cb).
 *
 * The integral is performed until the maximum value of k_max defined
 * in the perturbation module. Here there is not automatic checking
 * that k_max is large enough for the result to be well
 * converged. E.g. to get an accurate sigma8 at R = 8 Mpc/h, the user
 * should pass at least about P_k_max_h/Mpc = 1.
 *
 * @param ppr          Input: pointer to precision structure
 * @param pba          Input: pointer to background structure
 * @param pfo          Input: pointer to fourier structure
 * @param R            Input: radius in Mpc
 * @param z            Input: redshift
 * @param index_pk     Input: type of pk (_m, _cb)
 * @param sigma_output Input: quantity to be computed (sigma, sigma', ...)
 * @param result       Output: result
 * @return the error status
 */

int fourier_sigmas_at_z(
                        struct precision * ppr,
                        struct background * pba,
                        struct fourier * pfo,
                        double R,
                        double z,
                        int index_pk,
                        enum out_sigmas sigma_output,
                        double * result
                        ) {

  double * out_pk;
  double * ddout_pk;

  /** - allocate temporary array for P(k,z) as a function of k */

  class_alloc(out_pk, pfo->k_size*sizeof(double), pfo->error_message);
  class_alloc(ddout_pk, pfo->k_size*sizeof(double), pfo->error_message);

  /** - get P(k,z) as a function of k, for the right z */

  class_call(fourier_pk_at_z(pba,
                             pfo,
                             logarithmic,
                             pk_linear,
                             z,
                             index_pk,
                             out_pk,
                             NULL),
             pfo->error_message,
             pfo->error_message);

  /** - spline it along k */

  class_call(array_spline_table_columns(pfo->ln_k,
                                        pfo->k_size,
                                        out_pk,
                                        1,
                                        ddout_pk,
                                        _SPLINE_EST_DERIV_,
                                        pfo->error_message),
             pfo->error_message,
             pfo->error_message);

  /** - calll the function computing the sigmas */

  class_call(fourier_sigmas(pfo,
                            R,
                            out_pk,
                            ddout_pk,
                            pfo->k_size,
                            ppr->sigma_k_per_decade,
                            sigma_output,
                            result),
             pfo->error_message,
             pfo->error_message);

  /** - free allocated arrays */

  free(out_pk);
  free(ddout_pk);

  return _SUCCESS_;
}

/**
 * Return the value of the non-linearity wavenumber k_nl for a given redshift z
 *
 * @param pba     Input: pointer to background structure
 * @param pfo     Input: pointer to fourier structure
 * @param z       Input: redshift
 * @param k_nl    Output: k_nl value
 * @param k_nl_cb Ouput: k_nl value of the cdm+baryon part only, if there is ncdm
 * @return the error status
 */

int fourier_k_nl_at_z(
                      struct background *pba,
                      struct fourier * pfo,
                      double z,
                      double * k_nl,
                      double * k_nl_cb
                      ) {

  double tau;

  /** - convert input redshift into a conformal time */

  class_call(background_tau_of_z(pba,
                                 z,
                                 &tau),
             pba->error_message,
             pfo->error_message);

  /** - interpolate the precomputed k_nl array at the needed valuetime */

  if (pfo->has_pk_m == _TRUE_) {

    if (pfo->tau_size == 1) {
      *k_nl = pfo->k_nl[pfo->index_pk_m][0];
    }
    else {
      class_call(array_interpolate_two(pfo->tau,
                                       1,
                                       0,
                                       pfo->k_nl[pfo->index_pk_m],
                                       1,
                                       pfo->tau_size,
                                       tau,
                                       k_nl,
                                       1,
                                       pfo->error_message),
                 pfo->error_message,
                 pfo->error_message);
    }
  }

  /** - if needed, do the same for the baryon part only */

  if (pfo->has_pk_cb == _TRUE_) {

    if (pfo->tau_size == 1) {
      *k_nl_cb = pfo->k_nl[pfo->index_pk_cb][0];
    }
    else {
      class_call(array_interpolate_two(pfo->tau,
                                       1,
                                       0,
                                       pfo->k_nl[pfo->index_pk_cb],
                                       1,
                                       pfo->tau_size,
                                       tau,
                                       k_nl_cb,
                                       1,
                                       pfo->error_message),
                 pfo->error_message,
                 pfo->error_message);
    }
  }

  /* otherwise, return the same for k_nl_cb as for k_nl */

  else{
    *k_nl_cb = *k_nl;
  }

  return _SUCCESS_;
}

/**
 * Initialize the fourier structure, and in particular the
 * nl_corr_density and k_nl interpolation tables.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param pth Input: pointer to therodynamics structure
 * @param ppt Input: pointer to perturbation structure
 * @param ppm Input: pointer to primordial structure
 * @param pfo Input/Output: pointer to initialized fourier structure
 * @return the error status
 */

int fourier_init(
                 struct precision *ppr,
                 struct background *pba,
                 struct thermodynamics *pth,
                 struct perturbations *ppt,
                 struct primordial *ppm,
                 struct fourier *pfo
                 ) {

  int index_ncdm;
  int index_k;
  int index_tau;
  int index_tau_sources;
  int index_tau_late;
  int index_pk;

  double **pk_nl;
  double **lnpk_l;
  double **ddlnpk_l;

  short nl_corr_not_computable_at_this_k = _FALSE_;

  double * pvecback;
  int last_index;
  double a,z;

  struct fourier_workspace nw;
  struct fourier_workspace * pnw;

  /** - Do we want to compute P(k,z)? Propagate the flag has_pk_matter
      from the perturbations structure to the fourier structure */
  pfo->has_pk_matter = ppt->has_pk_matter;

  /** - preliminary tests */

  /** --> This module only makes sense for dealing with scalar
      perturbations, so it should do nothing if there are no
      scalars */
  if (ppt->has_scalars == _FALSE_) {
    pfo->method = nl_none;
    if (pfo->fourier_verbose > 0)
      printf("No scalar modes requested. Nonlinear module skipped.\n");
    return _SUCCESS_;
  }

  /** --> Nothing to be done if we don't want the matter power spectrum */

  if ((pfo->has_pk_matter == _FALSE_) && (pfo->method == nl_none)) {
    if (pfo->fourier_verbose > 0)
      printf("No Fourier spectra nor nonlinear corrections requested. Nonlinear module skipped.\n");
    return _SUCCESS_;
  }
  else {
    if (pfo->fourier_verbose > 0)
      printf("Computing linear Fourier spectra.\n");
  }

  /** --> check applicability of Halofit and HMcode */
  if (pfo->method > nl_none) {

    if (pba->has_ncdm == _TRUE_) {
      for (index_ncdm=0;index_ncdm < pba->N_ncdm; index_ncdm++){
        if (pba->m_ncdm_in_eV[index_ncdm] >  _M_EV_TOO_BIG_FOR_HALOFIT_)
          fprintf(stdout,"Warning: Halofit and HMcode are proved to work for CDM, and also with a small HDM component. But it sounds like you are running with a WDM component of mass %f eV, which makes the use of Halofit suspicious.\n",pba->m_ncdm_in_eV[index_ncdm]);
      }
    }
    if ((pba->has_idm == _TRUE_) && (ppt->perturbations_verbose > 0)){
      fprintf(stdout,"Warning: Halofit and HMcode are proved to work for CDM, and also with a small HDM component. But you have requested interacting dark matter (idm), which makes the use of Halofit or HMCode unreliable.\n");
    }
  }

  /** - define indices in fourier structure (and allocate some arrays in the structure) */

  class_call(fourier_indices(
                             ppr,
                             pba,
                             ppt,
                             ppm,
                             pfo),
             pfo->error_message,
             pfo->error_message);

  /** - get the linear power spectrum at each time */

  for (index_tau=0; index_tau<pfo->ln_tau_size;index_tau++) {

    /* If the user only wants z=0, then pfo->ln_tau_size=1 and we go
       only through index_tau=0. However we must pick up the last
       value of the source, index_tau_sources = ppt->tau_size-1. If
       the user wants several values of z, they correspond to the last
       ppt->ln_tau_size values of the sources (those that were also
       part of the array ppt->late_sources in the perturbation
       module). In all cases, the following formula gives the
       correspondance between index_tau in the current array and in
       the sources array:
    */

    index_tau_sources = ppt->tau_size-ppt->ln_tau_size+index_tau;

    /** --> loop over required pk types (_m, _cb) */

    for (index_pk=0; index_pk<pfo->pk_size; index_pk++) {

      /** --> get the linear power spectrum for this time and this type */

      class_call(fourier_pk_linear(
                                   pba,
                                   ppt,
                                   ppm,
                                   pfo,
                                   index_pk,
                                   index_tau_sources,
                                   pfo->k_size,
                                   &(pfo->ln_pk_l[index_pk][index_tau * pfo->k_size]),
                                   &(pfo->ln_pk_ic_l[index_pk][index_tau * pfo->k_size * pfo->ic_ic_size])
                                   ),
                 pfo->error_message,
                 pfo->error_message);


      /** --> if interpolation of \f$P(k,\tau)\f$ will be needed (as a
          function of tau), compute array of second derivatives in view of
          spline interpolation */

      if (pfo->ln_tau_size > 1) {

        class_call(array_spline_table_lines(pfo->ln_tau,
                                            pfo->ln_tau_size,
                                            pfo->ln_pk_l[index_pk],
                                            pfo->k_size,
                                            pfo->ddln_pk_l[index_pk],
                                            _SPLINE_EST_DERIV_,
                                            pfo->error_message),
                   pfo->error_message,
                   pfo->error_message);

        class_call(array_spline_table_lines(pfo->ln_tau,
                                            pfo->ln_tau_size,
                                            pfo->ln_pk_ic_l[index_pk],
                                            pfo->k_size*pfo->ic_ic_size,
                                            pfo->ddln_pk_ic_l[index_pk],
                                            _SPLINE_EST_DERIV_,
                                            pfo->error_message),
                   pfo->error_message,
                   pfo->error_message);
      }
    }
  }

  /** - compute and store sigma8 (variance of density fluctuations in
      spheres of radius 8/h Mpc at z=0, always computed by
      convention using the linear power spectrum) */

  for (index_pk=0; index_pk<pfo->pk_size; index_pk++) {

    class_call(fourier_sigmas_at_z(ppr,
                                   pba,
                                   pfo,
                                   8./pba->h,
                                   0.,
                                   index_pk,
                                   out_sigma,
                                   &(pfo->sigma8[index_pk])),
               pfo->error_message,
               pfo->error_message);
  }

  if (pfo->fourier_verbose>0) {

    if (pfo->has_pk_m == _TRUE_)
      fprintf(stdout," -> sigma8=%g for total matter (computed till k = %g h/Mpc)\n",
              pfo->sigma8[pfo->index_pk_m],
              pfo->k[pfo->k_size-1]/pba->h);

    if (pfo->has_pk_cb == _TRUE_)
      fprintf(stdout," -> sigma8=%g for baryons+cdm  (computed till k = %g h/Mpc)\n",
              pfo->sigma8[pfo->index_pk_cb],
              pfo->k[pfo->k_size-1]/pba->h);
  }

  /** - get the non-linear power spectrum at each time */

  /** --> First deal with the case where non non-linear corrections requested */

  if (pfo->method == nl_none) {
    if (pfo->fourier_verbose > 0)
      printf("No non-linear spectra requested. Nonlinear calculations skipped.\n");
  }

  /** --> Then go through common preliminary steps to the HALOFIT and HMcode methods */
  else if ((pfo->method == nl_halofit) || ((pfo->method == nl_HMcode))) {

    if ((pfo->fourier_verbose > 0) && (pfo->method == nl_halofit))
      printf("Computing non-linear matter power spectrum with Halofit (including update Takahashi et al. 2012 and Bird 2014)\n");

	if ((pfo->fourier_verbose > 0) && (pfo->method == nl_HMcode))
      printf("Computing non-linear matter power spectrum with HMcode \n");

    /** --> allocate temporary arrays for spectra at each given time/redshift */

    class_alloc(pk_nl,
                pfo->k_size*sizeof(double),
                pfo->error_message);

    class_alloc(lnpk_l,
                pfo->k_size*sizeof(double),
                pfo->error_message);

    class_alloc(ddlnpk_l,
                pfo->k_size*sizeof(double),
                pfo->error_message);

    for (index_pk=0; index_pk<pfo->pk_size; index_pk++){
      class_alloc(pk_nl[index_pk],pfo->k_size*sizeof(double),pfo->error_message);
      class_alloc(lnpk_l[index_pk],pfo->k_size_extra*sizeof(double),pfo->error_message);
      class_alloc(ddlnpk_l[index_pk],pfo->k_size_extra*sizeof(double),pfo->error_message);
    }

    /** --> Then go through preliminary steps specific to HMcode */

    if (pfo->method == nl_HMcode){

      pnw = &nw;

      class_call(fourier_hmcode_workspace_init(ppr,pba,pfo,pnw),
                 pfo->error_message,
                 pfo->error_message);

      class_call(fourier_hmcode_dark_energy_correction(ppr,pba,pfo,pnw),
                 pfo->error_message,
                 pfo->error_message);

      class_call(fourier_hmcode_baryonic_feedback(pfo),
                 pfo->error_message,
                 pfo->error_message);
    }

    /** --> Loop over decreasing time/growing redhsift. For each
        time/redshift, compute P_NL(k,z) using either Halofit or
        HMcode */

    /* this flag will become _TRUE_ at the minimum redshift such that
       the non-lienar corrections cannot be consistently computed */
    nl_corr_not_computable_at_this_k = _FALSE_;

    /* this index will refer to the value of time corresponding to
       that redhsift */
    pfo->index_tau_min_nl = 0;

    for (index_tau = pfo->tau_size-1; index_tau>=0; index_tau--) {

      /* loop over index_pk, defined such that it is ensured
       * that index_pk starts at index_pk_cb when neutrinos are
       * included. This is necessary for hmcode, since the sigmatable
       * needs to be filled for sigma_cb only. Thus, when HMcode
       * evalutes P_m_nl, it needs both P_m_l and P_cb_l. */

      for (index_pk=0; index_pk<pfo->pk_size; index_pk++) {

        /* get P_L(k) at this time */
        class_call(fourier_pk_linear(
                                     pba,
                                     ppt,
                                     ppm,
                                     pfo,
                                     index_pk,
                                     index_tau,
                                     pfo->k_size_extra,
                                     lnpk_l[index_pk],
                                     NULL
                                     ),
                   pfo->error_message,
                   pfo->error_message);

        /* spline P_L(k) at this time along k */
        class_call(array_spline_table_columns(
                                              pfo->ln_k,
                                              pfo->k_size_extra,
                                              lnpk_l[index_pk],
                                              1,
                                              ddlnpk_l[index_pk],
                                              _SPLINE_NATURAL_,
                                              pfo->error_message),
                   pfo->error_message,
                   pfo->error_message);

        /* if we are still in a range of time where P_NL(k) should be computable */
        if (nl_corr_not_computable_at_this_k == _FALSE_) {

          /* get P_NL(k) at this time with Halofit */
          if (pfo->method == nl_halofit) {

            class_call(fourier_halofit(
                                       ppr,
                                       pba,
                                       ppt,
                                       ppm,
                                       pfo,
                                       index_pk,
                                       pfo->tau[index_tau],
                                       pk_nl[index_pk],
                                       lnpk_l[index_pk],
                                       ddlnpk_l[index_pk],
                                       &(pfo->k_nl[index_pk][index_tau]),
                                       &nl_corr_not_computable_at_this_k),
                       pfo->error_message,
                       pfo->error_message);

          }

          /* get P_NL(k) at this time with HMcode */
          else if (pfo->method == nl_HMcode) {

            /* (preliminary step: fill table of sigma's, only for _cb if there is both _cb and _m) */
            if (index_pk == 0) {
              class_call(fourier_hmcode_fill_sigtab(ppr,
                                                    pba,
                                                    ppt,
                                                    ppm,
                                                    pfo,
                                                    index_tau,
                                                    lnpk_l[index_pk],
                                                    ddlnpk_l[index_pk],
                                                    pnw),
                         pfo->error_message, pfo->error_message);
            }

            class_call(fourier_hmcode(ppr,
                                      pba,
                                      ppt,
                                      ppm,
                                      pfo,
                                      index_pk,
                                      index_tau,
                                      pfo->tau[index_tau],
                                      pk_nl[index_pk],
                                      lnpk_l,
                                      ddlnpk_l,
                                      &(pfo->k_nl[index_pk][index_tau]),
                                      &nl_corr_not_computable_at_this_k,
                                      pnw),
                       pfo->error_message,
                       pfo->error_message);
          }

          /* infer and store R_NL=(P_NL/P_L)^1/2 */
          if (nl_corr_not_computable_at_this_k == _FALSE_) {
            for (index_k=0; index_k<pfo->k_size; index_k++) {
              pfo->nl_corr_density[index_pk][index_tau * pfo->k_size + index_k] = sqrt(pk_nl[index_pk][index_k]/exp(lnpk_l[index_pk][index_k]));
            }
          }

          /* otherwise we met the first problematic value of time */
          else {

            /* store the index of that value */
            pfo->index_tau_min_nl = MIN(pfo->tau_size-1,index_tau+1); //this MIN() ensures that index_tau_min_nl is never out of bounds

            /* store R_NL=1 for that time */
            for (index_k=0; index_k<pfo->k_size; index_k++) {
              pfo->nl_corr_density[index_pk][index_tau * pfo->k_size + index_k] = 1.;
            }

            /* send a warning to inform user about the corresponding value of redshift */
            if (pfo->fourier_verbose > 0) {
              class_alloc(pvecback,pba->bg_size*sizeof(double),pfo->error_message);
              class_call(background_at_tau(pba,pfo->tau[index_tau],short_info,inter_normal,&last_index,pvecback),
                         pba->error_message,
                         pfo->error_message);
              a = pvecback[pba->index_bg_a];
              /* redshift (remeber that a in the code stands for (a/a_0)) */
              z = 1./a-1.;
              fprintf(stdout,
                      " -> [WARNING:] Non-linear corrections could not be computed at redshift z=%5.2f and higher.\n    This is because k_max is too small for the algorithm (Halofit or HMcode) to be able to compute the scale k_NL at this redshift.\n    If non-linear corrections at such high redshift really matter for you,\n    just try to increase the precision parameter nonlinear_min_k_max (currently at %e) until k_NL can be computed at the desired z.\n",z,ppr->nonlinear_min_k_max);

              free(pvecback);
            }
          }
        }

        /* if we are still in a range of time where P_NL(k) should NOT be computable */
        else {
          /* store R_NL=1 for that time */
          for (index_k=0; index_k<pfo->k_size; index_k++) {
            pfo->nl_corr_density[index_pk][index_tau * pfo->k_size + index_k] = 1.;
          }

        }

        /** --> fill the array of nonlinear power spectra (only if we
            are at a late time where P(k) and T(k) are supposed to
            be stored, i.e., such that z(tau < z_max_pk) */

        if (index_tau >= pfo->tau_size - pfo->ln_tau_size) {

          index_tau_late = index_tau - (pfo->tau_size - pfo->ln_tau_size);

          for (index_k=0; index_k<pfo->k_size; index_k++) {
            pfo->ln_pk_nl[index_pk][index_tau_late * pfo->k_size + index_k] = pfo->ln_pk_l[index_pk][index_tau_late * pfo->k_size + index_k] + 2.*log(pfo->nl_corr_density[index_pk][index_tau * pfo->k_size + index_k]);
          }
        }

      } // end loop over index_pk
    } //end loop over index_tau


    /** --> spline the array of nonlinear power spectrum */

    if (pfo->ln_tau_size > 1) {
      for (index_pk=0; index_pk<pfo->pk_size; index_pk++) {

        class_call(array_spline_table_lines(pfo->ln_tau,
                                            pfo->ln_tau_size,
                                            pfo->ln_pk_nl[index_pk],
                                            pfo->k_size,
                                            pfo->ddln_pk_nl[index_pk],
                                            _SPLINE_EST_DERIV_,
                                            pfo->error_message),
                   pfo->error_message,
                   pfo->error_message);
      }
    }

    /* --> free temporary arrays */

    for (index_pk=0; index_pk<pfo->pk_size; index_pk++){
      free(pk_nl[index_pk]);
      free(lnpk_l[index_pk]);
      free(ddlnpk_l[index_pk]);
    }

    free(pk_nl);
    free(lnpk_l);
    free(ddlnpk_l);

    /** --> free the nonlinear workspace */

    if (pfo->method == nl_HMcode) {

      class_call(fourier_hmcode_workspace_free(pfo,pnw),
                 pfo->error_message,
                 pfo->error_message);
    }
  }

  /** - if the nl_method could not be identified */
  else {
    class_stop(pfo->error_message,
               "Your non-linear method variable is set to %d, out of the range defined in fourier.h",pfo->method);
  }

  return _SUCCESS_;
}

/**
 * Free all memory space allocated by fourier_init().
 *
 *
 * @param pfo Input: pointer to fourier structure (to be freed)
 * @return the error status
 */

int fourier_free(
                 struct fourier *pfo
                 ) {
  int index_pk;

  if ((pfo->has_pk_matter == _TRUE_) || (pfo->method > nl_none)) {

    free(pfo->k);
    free(pfo->ln_k);

    for (index_pk=0; index_pk<pfo->pk_size; index_pk++) {
      free(pfo->ln_pk_ic_l[index_pk]);
      free(pfo->ln_pk_l[index_pk]);
      if (pfo->ln_tau_size>1) {
        free(pfo->ddln_pk_ic_l[index_pk]);
        free(pfo->ddln_pk_l[index_pk]);
      }
    }
    free(pfo->ln_pk_ic_l);
    free(pfo->ln_pk_l);

    free (pfo->sigma8);

    if (pfo->ln_tau_size>1) {
      free(pfo->ddln_pk_ic_l);
      free(pfo->ddln_pk_l);
      free(pfo->ln_tau);
    }

    free(pfo->is_non_zero);
  }

  if (pfo->method > nl_none) {

    free(pfo->tau);
    for (index_pk=0;index_pk<pfo->pk_size;index_pk++){
      free(pfo->nl_corr_density[index_pk]);
      free(pfo->k_nl[index_pk]);
      free(pfo->ln_pk_nl[index_pk]);
      if (pfo->ln_tau_size > 1)
        free(pfo->ddln_pk_nl[index_pk]);
    }
    free(pfo->nl_corr_density);
    free(pfo->k_nl);
    free(pfo->ln_pk_nl);
    if (pfo->ln_tau_size > 1)
      free(pfo->ddln_pk_nl);
  }

  if (pfo->has_pk_eq == _TRUE_) {
    free(pfo->pk_eq_tau);
    free(pfo->pk_eq_w_and_Omega);
    free(pfo->pk_eq_ddw_and_ddOmega);
  }

  return _SUCCESS_;
}

/**
 * Define indices in the fourier structure, and when possible, allocate
 * arrays in this structure given the index sizes found here
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param ppt Input: pointer to perturbation structure
 * @param ppm Input: pointer to primordial structure
 * @param pfo Input/Output: pointer to fourier structure
 * @return the error status
 */

int fourier_indices(
                    struct precision *ppr,
                    struct background *pba,
                    struct perturbations * ppt,
                    struct primordial * ppm,
                    struct fourier * pfo
                    ) {

  int index_ic1_ic2;
  int index_pk;

  /** - define indices for initial conditions (and allocate related arrays) */
  pfo->index_md_scalars = ppt->index_md_scalars;
  pfo->ic_size = ppm->ic_size[pfo->index_md_scalars];
  pfo->ic_ic_size = ppm->ic_ic_size[pfo->index_md_scalars];
  class_alloc(pfo->is_non_zero,sizeof(short)*pfo->ic_ic_size,pfo->error_message);
  for (index_ic1_ic2=0; index_ic1_ic2 < pfo->ic_ic_size; index_ic1_ic2++)
    pfo->is_non_zero[index_ic1_ic2] = ppm->is_non_zero[pfo->index_md_scalars][index_ic1_ic2];

  /** - define flags indices for pk types (_m, _cb). Note: due to some
      dependencies in HMcode, when pfo->index_pk_cb exists, it must
      come first (e.g. the calculation of the non-linear P_m depends on
      sigma_cb so the cb-related quantitites must be evaluated
      first) */

  pfo->has_pk_m = _TRUE_;
  if (pba->has_ncdm == _TRUE_) {
    pfo->has_pk_cb = _TRUE_;
  }
  else {
    pfo->has_pk_cb = _FALSE_;
  }

  index_pk = 0;
  class_define_index(pfo->index_pk_cb, pfo->has_pk_cb, index_pk,1);
  class_define_index(pfo->index_pk_m, pfo->has_pk_m, index_pk,1);
  pfo->pk_size = index_pk;

  /* and two redundent but useful indices: */

  if (pfo->has_pk_cb == _TRUE_) {
    pfo->index_pk_total = pfo->index_pk_m;
    pfo->index_pk_cluster = pfo->index_pk_cb;
  }
  else {
    pfo->index_pk_total = pfo->index_pk_m;
    pfo->index_pk_cluster = pfo->index_pk_m;
  }

  /** - get list of k values */

  class_call(fourier_get_k_list(ppr,ppt,pfo),
             pfo->error_message,
             pfo->error_message);

  /** - get list of tau values */

  class_call(fourier_get_tau_list(ppt,pfo),
             pfo->error_message,
             pfo->error_message);

  /** - given previous indices, we can allocate the array of linear power spectrum values */

  class_alloc(pfo->ln_pk_ic_l,pfo->pk_size*sizeof(double*),pfo->error_message);
  class_alloc(pfo->ln_pk_l   ,pfo->pk_size*sizeof(double*),pfo->error_message);

  for (index_pk=0; index_pk<pfo->pk_size; index_pk++) {
    class_alloc(pfo->ln_pk_ic_l[index_pk],pfo->ln_tau_size*pfo->k_size*pfo->ic_ic_size*sizeof(double*),pfo->error_message);
    class_alloc(pfo->ln_pk_l[index_pk]   ,pfo->ln_tau_size*pfo->k_size*sizeof(double*),pfo->error_message);
  }

  /** - if interpolation of \f$P(k,\tau)\f$ will be needed (as a function of tau),
      compute also the array of second derivatives in view of spline interpolation */

  if (pfo->ln_tau_size > 1) {

    class_alloc(pfo->ddln_pk_ic_l,pfo->pk_size*sizeof(double*),pfo->error_message);
    class_alloc(pfo->ddln_pk_l   ,pfo->pk_size*sizeof(double*),pfo->error_message);

    for (index_pk=0; index_pk<pfo->pk_size; index_pk++) {
      class_alloc(pfo->ddln_pk_ic_l[index_pk],pfo->ln_tau_size*pfo->k_size*pfo->ic_ic_size*sizeof(double*),pfo->error_message);
      class_alloc(pfo->ddln_pk_l[index_pk]   ,pfo->ln_tau_size*pfo->k_size*sizeof(double*),pfo->error_message);
    }
  }

  /** - array of sigma8 values */

  class_alloc(pfo->sigma8,pfo->pk_size*sizeof(double*),pfo->error_message);

  /** - if non-linear computations needed, allocate array of
      non-linear correction ratio R_nl(k,z), k_nl(z) and P_nl(k,z)
      for each P(k) type */

  if (pfo->method > nl_none) {

    class_alloc(pfo->k_nl,
                pfo->pk_size*sizeof(double *),
                pfo->error_message);

    class_alloc(pfo->nl_corr_density,
                pfo->pk_size*sizeof(double *),
                pfo->error_message);

    class_alloc(pfo->ln_pk_nl,pfo->pk_size*sizeof(double*),pfo->error_message);
    if (pfo->ln_tau_size > 1)
      class_alloc(pfo->ddln_pk_nl,pfo->pk_size*sizeof(double*),pfo->error_message);

    for (index_pk=0; index_pk<pfo->pk_size; index_pk++){
      class_alloc(pfo->k_nl[index_pk],pfo->tau_size*sizeof(double),pfo->error_message);
      class_alloc(pfo->nl_corr_density[index_pk],pfo->tau_size*pfo->k_size*sizeof(double),pfo->error_message);
      class_alloc(pfo->ln_pk_nl[index_pk],pfo->ln_tau_size*pfo->k_size*sizeof(double*),pfo->error_message);
      if (pfo->ln_tau_size > 1)
        class_alloc(pfo->ddln_pk_nl[index_pk],pfo->ln_tau_size*pfo->k_size*sizeof(double*),pfo->error_message);
    }
  }

  return _SUCCESS_;
}

/**
 * Copy list of k from perturbation module, and extended it if
 * necessary to larger k for extrapolation (currently this
 * extrapolation is required only by HMcode)
 *
 * @param ppr Input: pointer to precision structure
 * @param ppt Input: pointer to perturbation structure
 * @param pfo Input/Output: pointer to fourier structure
 * @return the error status
 */

int fourier_get_k_list(
                       struct precision *ppr,
                       struct perturbations * ppt,
                       struct fourier * pfo
                       ) {

  double k=0;
  double k_max,exponent;
  int index_k;

  pfo->k_size = ppt->k_size[pfo->index_md_scalars];
  pfo->k_size_pk = ppt->k_size_pk;
  k_max = ppt->k[pfo->index_md_scalars][pfo->k_size-1];

  /** - if k extrapolation necessary, compute number of required extra values */
  if (pfo->method == nl_HMcode){
    index_k=0;
    while(k < ppr->hmcode_max_k_extra && index_k < _MAX_NUM_EXTRAPOLATION_){
      index_k++;
      k = k_max * pow(10,(double)index_k/ppr->k_per_decade_for_pk);
    }
    class_test(index_k == _MAX_NUM_EXTRAPOLATION_,
               pfo->error_message,
               "could not reach extrapolated value k = %.10e starting from k = %.10e with k_per_decade of %.10e in _MAX_NUM_INTERPOLATION_=%i steps",
               ppr->hmcode_max_k_extra,k_max,ppr->k_per_decade_for_pk,_MAX_NUM_EXTRAPOLATION_
               );
    pfo->k_size_extra = pfo->k_size+index_k;
  }
  /** - otherwise, same number of values as in perturbation module */
  else {
    pfo->k_size_extra = pfo->k_size;
  }

  /** - allocate array of k */
  class_alloc(pfo->k,   pfo->k_size_extra*sizeof(double),pfo->error_message);
  class_alloc(pfo->ln_k,pfo->k_size_extra*sizeof(double),pfo->error_message);

  /** - fill array of k (not extrapolated) */
  for (index_k=0; index_k<pfo->k_size; index_k++) {
    k = ppt->k[pfo->index_md_scalars][index_k];
    pfo->k[index_k] = k;
    pfo->ln_k[index_k] = log(k);
  }

  /** - fill additional values of k (extrapolated) */
  for (index_k=pfo->k_size; index_k<pfo->k_size_extra; index_k++) {
    exponent = (double)(index_k-(pfo->k_size-1))/ppr->k_per_decade_for_pk;
    pfo->k[index_k] = k * pow(10,exponent);
    pfo->ln_k[index_k] = log(k) + exponent*log(10.);
  }

  return _SUCCESS_;
}

/**
 * Copy list of tau from perturbation module
 *
 * @param ppt Input: pointer to perturbation structure
 * @param pfo Input/Output: pointer to fourier structure
 * @return the error status
 */

int fourier_get_tau_list(
                         struct perturbations * ppt,
                         struct fourier * pfo
                         ) {

  int index_tau;

  /** -> for linear calculations: only late times are considered, given the value z_max_pk inferred from the ionput */
  pfo->ln_tau_size = ppt->ln_tau_size;
  pfo->index_ln_tau_pk = ppt->index_ln_tau_pk;

  if (ppt->ln_tau_size > 1) {

    class_alloc(pfo->ln_tau,pfo->ln_tau_size*sizeof(double),pfo->error_message);

    for (index_tau=0; index_tau<pfo->ln_tau_size;index_tau++) {
      pfo->ln_tau[index_tau] = ppt->ln_tau[index_tau];
    }
  }

  /** -> for non-linear calculations: we wills store a correction factor for all times */
  if (pfo->method > nl_none) {

    pfo->tau_size = ppt->tau_size;

    class_alloc(pfo->tau,pfo->tau_size*sizeof(double),pfo->error_message);

    for (index_tau=0; index_tau<pfo->tau_size; index_tau++) {
      pfo->tau[index_tau] = ppt->tau_sampling[index_tau];
    }
  }
  return _SUCCESS_;
}

/**
 * Get sources for a given wavenumber (and for a given time, type, ic,
 * mode...) either directly from precomputed valkues (computed ain
 * perturbation module), or by analytic extrapolation
 *
 * @param pba             Input: pointer to background structure
 * @param ppt             Input: pointer to perturbation structure
 * @param pfo             Input: pointer to fourier structure
 * @param index_k         Input: index of required k value
 * @param index_ic        Input: index of required ic value
 * @param index_tp        Input: index of required tp value
 * @param index_tau       Input: index of required tau value
 * @param sources         Input: array containing the original sources
 * @param source          Output: desired value of source
 * @return the error status
 */

int fourier_get_source(
                       struct background * pba,
                       struct perturbations * ppt,
                       struct fourier * pfo,
                       int index_k,
                       int index_ic,
                       int index_tp,
                       int index_tau,
                       double ** sources,
                       double * source
                       ) {

  double k,k_max,k_previous;
  double source_max,source_previous;
  double scaled_factor,log_scaled_factor;

  /** - use precomputed values */
  if (index_k < pfo->k_size) {
    *source = sources[index_ic * ppt->tp_size[pfo->index_md_scalars] + index_tp][index_tau * pfo->k_size + index_k];
  }
  /** - extrapolate **/
  else {

    k = pfo->k[index_k];

    /**
     * --> Get last source and k, which are used in (almost) all methods
     */
    k_max = pfo->k[pfo->k_size-1];
    source_max = sources[index_ic * ppt->tp_size[pfo->index_md_scalars] + index_tp][index_tau * pfo->k_size + pfo->k_size - 1];

    /**
     * --> Get previous source and k, which are used in best methods
     */
    k_previous = pfo->k[pfo->k_size-2];
    source_previous = sources[index_ic * ppt->tp_size[pfo->index_md_scalars] + index_tp][index_tau * pfo->k_size + pfo->k_size - 2];

    switch(pfo->extrapolation_method){
      /**
       * --> Extrapolate by assuming the source to vanish Has terrible
       * discontinuity
       */
    case extrap_zero:
      {
        *source=0.0;
        break;
      }
      /**
       * --> Extrapolate starting from the maximum value, assuming  growth ~ ln(k)
       * Has a terrible bend in log slope, discontinuity only in derivative
       */
    case extrap_only_max:
      {
        *source = source_max*(log(k)/log(k_max));
        break;
      }
      /**
       * --> Extrapolate starting from the maximum value, assuming
       * growth ~ ln(k) Here we use k in h/Mpc instead of 1/Mpc as it
       * is done in the CAMB implementation of HMcode Has a terrible
       * bend in log slope, discontinuity only in derivative
       */
    case extrap_only_max_units:
      {
        *source = source_max*(log(k/pba->h)/log(k_max/pba->h));
        break;
      }
      /**
       * --> Extrapolate assuming source ~ ln(a*k) where a is obtained
       * from the data at k_0 Mostly continuous derivative, quite good
       */
    case extrap_max_scaled:
      {
        log_scaled_factor = (source_previous*log(k_max)-source_max*log(k_previous))/(source_max-source_previous);
        *source = source_max*((log_scaled_factor+log(k))/(log_scaled_factor+log(k_max)));
        break;
      }
      /**
       * --> Extrapolate assuming source ~ ln(e+a*k) where a is
       * estimated like is done in original HMCode
       */
    case extrap_hmcode:
      {
        scaled_factor = 1.8/(13.41*pba->a_eq*pba->H_eq);
        *source = source_max*(log(_E_+scaled_factor*k)/log(_E_+scaled_factor*k_max));
        break;
      }
      /**
       * --> If the user has a complicated model and wants to
       * interpolate differently, they can define their interpolation
       * here and switch to using it instead
       */
    case extrap_user_defined:
      {
        class_stop(pfo->error_message,"Method of source extrapolation 'user_defined' was not yet defined.");
        break;
      }
    }
  }
  return _SUCCESS_;
}

/**
 * This routine computes all the components of the matter power
 * spectrum P(k), given the source functions and the primordial
 * spectra, at a given time within the pre-computed table of sources
 * (= Fourier transfer functions) of the perturbation module, for a
 * given type (total matter _m or baryon+CDM _cb), and for the same
 * array of k values as in the pre-computed table.
 *
 * If the input array of k values pfo->ln_k contains wavemumbers
 * larger than those of the pre-computed table, the sources will be
 * extrapolated analytically.
 *
 * On the opther hand, if the primordial spectrum has sharp features
 * and needs to be sampled on a finer grid than the sources, this
 * function has to be modified to capture the features.
 *
 * There are two output arrays, because we consider:
 *
 * - the total matter (_m) or CDM+baryon (_cb) power spectrum
 *
 * - in the quantitites labelled _ic, the splitting of one of these
 * spectra in different modes for different initial conditions. If the
 * pointer ln_pk_ic is NULL in input, the function will ignore this
 * part; thus, to get the result, one should allocate the array before
 * calling the function. Then the convention is the following:
 *
 * -- the index_ic1_ic2 labels ordered pairs (index_ic1, index_ic2)
 * (since the primordial spectrum is symmetric in (index_ic1,
 * index_ic2)).
 *
 * -- for diagonal elements (index_ic1 = index_ic2) this
 * arrays contains ln[P(k)] where P(k) is positive by construction.
 *
 * -- for non-diagonal elements this arrays contains the k-dependent
 * cosine of the correlation angle, namely P(k)_(index_ic1,
 * index_ic2)/sqrt[P(k)_index_ic1 P(k)_index_ic2]. E.g. for fully
 * correlated or anti-correlated initial conditions, this non-diagonal
 * element is independent on k, and equal to +1 or -1.
 *
 * @param pba           Input: pointer to background structure
 * @param ppt           Input: pointer to perturbation structure
 * @param ppm           Input: pointer to primordial structure
 * @param pfo           Input: pointer to fourier structure
 * @param index_pk      Input: index of required P(k) type (_m, _cb)
 * @param index_tau     Input: index of time
 * @param k_size        Input: wavenumber array size
 * @param lnpk         Output: log of matter power spectrum for given type/time, for all wavenumbers
 * @param lnpk_ic      Output: log of matter power spectrum for given type/time, for all wavenumbers and initial conditions
 * @return the error status
 */

int fourier_pk_linear(
                      struct background *pba,
                      struct perturbations *ppt,
                      struct primordial *ppm,
                      struct fourier *pfo,
                      int index_pk,
                      int index_tau,
                      int k_size,
                      double * lnpk,    //lnpk[index_k]
                      double * lnpk_ic  //lnpk[index_k * pfo->ic_ic_size + index_ic1_ic2]
                      ) {

  int index_k;
  int index_tp;
  int index_ic1,index_ic2,index_ic1_ic1,index_ic1_ic2,index_ic2_ic2;
  double * primordial_pk;
  double pk;
  double * pk_ic;
  double source_ic1;
  double source_ic2;
  double cosine_correlation;

  /** - allocate temporary vector where the primordial spectrum will be stored */

  class_alloc(primordial_pk,pfo->ic_ic_size*sizeof(double),pfo->error_message);

  class_alloc(pk_ic,pfo->ic_ic_size*sizeof(double),pfo->error_message);

  if ((pfo->has_pk_m == _TRUE_) && (index_pk == pfo->index_pk_m)) {
    index_tp = ppt->index_tp_delta_m;
  }
  else if ((pfo->has_pk_cb == _TRUE_) && (index_pk == pfo->index_pk_cb)) {
    index_tp = ppt->index_tp_delta_cb;
  }
  else {
    class_stop(pfo->error_message,"P(k) is set neither to total matter nor to cold dark matter + baryons");
  }

  /** - loop over k values */

  for (index_k=0; index_k<k_size; index_k++) {

    /** --> get primordial spectrum */
    class_call(primordial_spectrum_at_k(ppm,pfo->index_md_scalars,logarithmic,pfo->ln_k[index_k],primordial_pk),
               ppm->error_message,
               pfo->error_message);

    /** --> initialize a local variable for P_m(k) and P_cb(k) to zero */
    pk = 0.;

    /** --> here we recall the relations relevant for the nomalization fo the power spectrum:
        For adiabatic modes, the curvature primordial spectrum thnat we just read was:
        P_R(k) = 1/(2pi^2) k^3 < R R >
        Thus the primordial curvature correlator is given by:
        < R R > = (2pi^2) k^-3 P_R(k)
        So the delta_m correlator reads:
        P(k) = < delta_m delta_m > = (source_m)^2 < R R > = (2pi^2) k^-3 (source_m)^2 P_R(k)

        For isocurvature or cross adiabatic-isocurvature parts,
        one would just replace one or two 'R' by 'S_i's */

    /** --> get contributions to P(k) diagonal in the initial conditions */
    for (index_ic1 = 0; index_ic1 < pfo->ic_size; index_ic1++) {

      index_ic1_ic1 = index_symmetric_matrix(index_ic1,index_ic1,pfo->ic_size);

      class_call(fourier_get_source(pba,
                                    ppt,
                                    pfo,
                                    index_k,
                                    index_ic1,
                                    index_tp,
                                    index_tau,
                                    ppt->sources[pfo->index_md_scalars],
                                    &source_ic1),
                 pfo->error_message,
                 pfo->error_message);

      pk_ic[index_ic1_ic1] = 2.*_PI_*_PI_/exp(3.*pfo->ln_k[index_k])
        *source_ic1*source_ic1
        *exp(primordial_pk[index_ic1_ic1]);

      pk += pk_ic[index_ic1_ic1];

      if (lnpk_ic != NULL) {
        lnpk_ic[index_k * pfo->ic_ic_size + index_ic1_ic1] = log(pk_ic[index_ic1_ic1]);
      }
    }

    /** --> get contributions to P(k) non-diagonal in the initial conditions */
    for (index_ic1 = 0; index_ic1 < pfo->ic_size; index_ic1++) {
      for (index_ic2 = index_ic1+1; index_ic2 < pfo->ic_size; index_ic2++) {

        index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,pfo->ic_size);
        index_ic1_ic1 = index_symmetric_matrix(index_ic1,index_ic1,pfo->ic_size);
        index_ic2_ic2 = index_symmetric_matrix(index_ic2,index_ic2,pfo->ic_size);

        if (pfo->is_non_zero[index_ic1_ic2] == _TRUE_) {

          class_call(fourier_get_source(pba,
                                        ppt,
                                        pfo,
                                        index_k,
                                        index_ic1,
                                        index_tp,
                                        index_tau,
                                        ppt->sources[pfo->index_md_scalars],
                                        &source_ic1),
                     pfo->error_message,
                     pfo->error_message);

          class_call(fourier_get_source(pba,
                                        ppt,
                                        pfo,
                                        index_k,
                                        index_ic2,
                                        index_tp,
                                        index_tau,
                                        ppt->sources[pfo->index_md_scalars],
                                        &source_ic2),
                     pfo->error_message,
                     pfo->error_message);

          cosine_correlation = primordial_pk[index_ic1_ic2]*SIGN(source_ic1)*SIGN(source_ic2);

          pk_ic[index_ic1_ic2] = cosine_correlation * sqrt(pk_ic[index_ic1_ic1]*pk_ic[index_ic2_ic2]);

          pk += 2.*pk_ic[index_ic1_ic2];

          if (lnpk_ic != NULL) {
            lnpk_ic[index_k * pfo->ic_ic_size + index_ic1_ic2] = cosine_correlation;
          }
        }
        else {
          if (lnpk_ic != NULL) {
            lnpk_ic[index_k * pfo->ic_ic_size + index_ic1_ic2] = 0.;
          }
        }
      }
    }

    lnpk[index_k] = log(pk);
  }

  free(primordial_pk);
  free(pk_ic);

  return _SUCCESS_;

}

/**
 * Calculate intermediate quantities for hmcode (sigma, sigma', ...)
 * for a given scale R and a given input P(k).
 *
 * This function has several differences w.r.t. the standard external
 * function non_linear_sigma (format of input, of output, integration
 * stepsize, management of extrapolation at large k, ...) and is
 * overall more precise for sigma(R).
 *
 * @param pfo          Input: pointer to fourier structure
 * @param R            Input: scale at which to compute sigma
 * @param lnpk_l       Input: array of ln(P(k))
 * @param ddlnpk_l     Input: its spline along k
 * @param k_size       Input: dimension of array lnpk_l, normally pfo->k_size, but inside hmcode it its increased by extrapolation to pfo->k_extra_size
 * @param k_per_decade Input: logarithmic step for the integral (recommended: pass ppr->sigma_k_per_decade)
 * @param sigma_output Input: quantity to be computed (sigma, sigma', ...)
 * @param result       Output: result
 * @return the error status
 */

int fourier_sigmas(
                   struct fourier * pfo,
                   double R,
                   double * lnpk_l,
                   double * ddlnpk_l,
                   int k_size,
                   double k_per_decade,
                   enum out_sigmas sigma_output,
                   double * result
                   ) {
  double pk, lnpk;

  double * array_for_sigma;
  int index_num;
  int index_x;
  int index_y;
  int index_ddy;
  int i=0;
  int integrand_size;
  int last_index=0;

  double k,W,W_prime,x,t;

  /** - allocate temporary array for an integral over y(x) */

  class_define_index(index_x,  _TRUE_,i,1); // index for x
  class_define_index(index_y,  _TRUE_,i,1); // index for integrand
  class_define_index(index_ddy,_TRUE_,i,1); // index for its second derivative (spline method)
  index_num=i;                              // number of columns in the array

  integrand_size=(int)(log(pfo->k[k_size-1]/pfo->k[0])/log(10.)*k_per_decade)+1;
  class_alloc(array_for_sigma,
              integrand_size*index_num*sizeof(double),
              pfo->error_message);

  /** - fill the array with values of k and of the integrand */

  for (i=0; i<integrand_size; i++) {

    k=pfo->k[0]*pow(10.,i/k_per_decade);

    if (i==0) {
      pk = exp(lnpk_l[0]);
    }
    else {
      class_call(array_interpolate_spline(
                                          pfo->ln_k,
                                          k_size,
                                          lnpk_l,
                                          ddlnpk_l,
                                          1,
                                          log(k),
                                          &last_index,
                                          &lnpk,
                                          1,
                                          pfo->error_message),
                 pfo->error_message,
                 pfo->error_message);

      pk = exp(lnpk);
    }

    t = 1./(1.+k);
    if (i == (integrand_size-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
    x=k*R;

    switch (sigma_output) {

    case out_sigma:
      if (x<0.01)
        W = 1.-x*x/10.;
      else
        W = 3./x/x/x*(sin(x)-x*cos(x));
      array_for_sigma[(integrand_size-1-i)*index_num+index_x] = t;
      array_for_sigma[(integrand_size-1-i)*index_num+index_y] = k*k*k*pk*W*W/(t*(1.-t));
      break;

    case out_sigma_prime:
      if (x<0.01) {
        W = 1.-x*x/10.;
        W_prime = -0.2*x;
      }
      else {
        W = 3./x/x/x*(sin(x)-x*cos(x));
        W_prime = 3./x/x*sin(x)-9./x/x/x/x*(sin(x)-x*cos(x));
      }
      array_for_sigma[(integrand_size-1-i)*index_num+index_x] = t;
      array_for_sigma[(integrand_size-1-i)*index_num+index_y] = k*k*k*pk*2.*k*W*W_prime/(t*(1.-t));
      break;

    case out_sigma_disp:
      if (x<0.01)
        W = 1.-x*x/10.;
      else
        W = 3./x/x/x*(sin(x)-x*cos(x));
      array_for_sigma[(integrand_size-1-i)*index_num+index_x] = k;
      array_for_sigma[(integrand_size-1-i)*index_num+index_y] = -pk*W*W;
      break;
    }
  }

  /** - spline the integrand */

  class_call(array_spline(array_for_sigma,
                          index_num,
                          integrand_size,
                          index_x,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          pfo->error_message),
             pfo->error_message,
             pfo->error_message);

  /** - integrate */

  class_call(array_integrate_all_trapzd_or_spline(array_for_sigma,
                                                  index_num,
                                                  integrand_size,
                                                  0, //integrand_size-1,
                                                  index_x,
                                                  index_y,
                                                  index_ddy,
                                                  result,
                                                  pfo->error_message),
             pfo->error_message,
             pfo->error_message);

  /** - preperly normalize the final result */

  switch (sigma_output) {

  case out_sigma:
    *result = sqrt(*result/(2.*_PI_*_PI_));
    break;

  case out_sigma_prime:
    *result = *result/(2.*_PI_*_PI_);
    break;

  case out_sigma_disp:
    *result = sqrt(*result/(2.*_PI_*_PI_*3.));
    break;
  }

  /** - free allocated array */

  free(array_for_sigma);

  return _SUCCESS_;
}

/**
 * This routine computes the variance of density fluctuations in a
 * sphere of radius R at redshift z, sigma(R,z) for one given pk type (_m, _cb).
 *
 * Try to use instead fourier_sigmas_at_z(). This function is just
 * maintained for compatibility with the deprecated function
 * harmonic_sigma()
 *
 * The integral is performed until the maximum value of k_max defined
 * in the perturbation module. Here there is not automatic checking
 * that k_max is large enough for the result to be well
 * converged. E.g. to get an accurate sigma8 at R = 8 Mpc/h, the user
 * should pass at least about P_k_max_h/Mpc = 1.
 *
 * @param pba          Input: pointer to background structure
 * @param pfo          Input: pointer to fourier structure
 * @param R            Input: radius in Mpc
 * @param z            Input: redshift
 * @param index_pk     Input: type of pk (_m, _cb)
 * @param k_per_decade Input: logarithmic step for the integral (recommended: pass ppr->sigma_k_per_decade)
 * @param result       Output: result
 * @return the error status
 */

int fourier_sigma_at_z(
                       struct background * pba,
                       struct fourier * pfo,
                       double R,
                       double z,
                       int index_pk,
                       double k_per_decade,
                       double * result
                       ) {

  double * out_pk;
  double * ddout_pk;

  /** - allocate temporary array for P(k,z) as a function of k */

  class_alloc(out_pk, pfo->k_size*sizeof(double), pfo->error_message);
  class_alloc(ddout_pk, pfo->k_size*sizeof(double), pfo->error_message);

  /** - get P(k,z) as a function of k, for the right z */

  class_call(fourier_pk_at_z(pba,
                             pfo,
                             logarithmic,
                             pk_linear,
                             z,
                             index_pk,
                             out_pk,
                             NULL),
             pfo->error_message,
             pfo->error_message);

  /** - spline it along k */

  class_call(array_spline_table_columns(pfo->ln_k,
                                        pfo->k_size,
                                        out_pk,
                                        1,
                                        ddout_pk,
                                        _SPLINE_EST_DERIV_,
                                        pfo->error_message),
             pfo->error_message,
             pfo->error_message);

  /** - calll the function computing the sigmas */

  class_call(fourier_sigmas(pfo,
                            R,
                            out_pk,
                            ddout_pk,
                            pfo->k_size,
                            k_per_decade,
                            out_sigma,
                            result),
             pfo->error_message,
             pfo->error_message);

  /** - free allocated arrays */

  free(out_pk);
  free(ddout_pk);

  return _SUCCESS_;
}

/**
 * Calculation of the nonlinear matter power spectrum with Halofit
 * (includes Takahashi 2012 + Bird 2013 revisions).
 *
 * At high redshift it is possible that the non-linear corrections are
 * so small that they can be computed only by going to very large
 * wavenumbers. Thius, for some combination of (z, k_max), the
 * calculation is not possible. In this case a _FALSE_ will be
 * returned in the flag halofit_found_k_max.
 *
 * @param ppr         Input: pointer to precision structure
 * @param pba         Input: pointer to background structure
 * @param ppt         Input: pointer to perturbation structure
 * @param ppm         Input: pointer to primordial structure
 * @param pfo         Input: pointer to fourier structure
 * @param index_pk    Input: index of component are we looking at (total matter or cdm+baryons?)
 * @param tau         Input: conformal time at which we want to do the calculation
 * @param pk_nl       Output: non linear spectrum at the relevant time
 * @param lnpk_l      Input: array of log(P(k)_linear)
 * @param ddlnpk_l    Input: array of second derivative of log(P(k)_linear) wrt k, for spline interpolation
 * @param k_nl        Output: non-linear wavenumber
 * @param nl_corr_not_computable_at_this_k Ouput: flag concerning the status of the calculation (_TRUE_ if not possible)
 * @return the error status
 */

int fourier_halofit(
                    struct precision *ppr,
                    struct background *pba,
                    struct perturbations *ppt,
                    struct primordial *ppm,
                    struct fourier *pfo,
                    int index_pk,
                    double tau,
                    double *pk_nl,
                    double *lnpk_l,
                    double *ddlnpk_l,
                    double *k_nl,
                    short * nl_corr_not_computable_at_this_k
                    ) {

  double Omega_m,Omega_v,fnu,w, dw_over_da_fld, integral_fld;

  /** Determine non linear ratios (from pk) **/

  int index_k;
  double pk_lin,pk_quasi,pk_halo,rk;
  double sigma,rknl,rneff,rncur,d1,d2;
  double diff,xlogr1,xlogr2,rmid;

  double gam,a,b,c,xmu,xnu,alpha,beta,f1,f2,f3;
  double pk_linaa;
  double y;
  double f1a,f2a,f3a,f1b,f2b,f3b,frac;

  double * pvecback;

  int last_index=0;
  int counter;
  double sum1,sum2,sum3;
  double anorm;

  double *integrand_array;
  int integrand_size;
  int index_ia_k;
  int index_ia_pk;
  int index_ia_sum;
  int index_ia_ddsum;
  /*
    int index_ia_sum2;
    int index_ia_ddsum2;
    int index_ia_sum3;
    int index_ia_ddsum3;
  */
  int ia_size;
  int index_ia;

  double k_integrand;
  double lnpk_integrand;

  double R;

  double * w_and_Omega;

  class_alloc(pvecback,pba->bg_size*sizeof(double),pfo->error_message);

  if ((pfo->has_pk_m == _TRUE_) && (index_pk == pfo->index_pk_m)) {
    fnu = pba->Omega0_ncdm_tot/pba->Omega0_m;
  }
  else if ((pfo->has_pk_cb == _TRUE_) && (index_pk == pfo->index_pk_cb)) {
    fnu = 0.;
  }
  else {
    class_stop(pfo->error_message,"P(k) is set neither to total matter nor to cold dark matter + baryons");
  }

  if (pfo->has_pk_eq == _FALSE_) {

    /* default method: compute w(tau) = w_fld(tau), Omega_m(tau) and Omega_v=Omega_DE(tau), all required by HALOFIT fitting formulas */

    class_call(background_at_tau(pba,tau,long_info,inter_normal,&last_index,pvecback),
               pba->error_message,
               pfo->error_message);

    Omega_m = pvecback[pba->index_bg_Omega_m];
    Omega_v = 1.-pvecback[pba->index_bg_Omega_m]-pvecback[pba->index_bg_Omega_r];
    /* until v2.9.3 this function was called at a_0=1 instead of a=pvecback[pba->index_bg_a] */
    class_call(background_w_fld(pba,pvecback[pba->index_bg_a],&w,&dw_over_da_fld,&integral_fld), pba->error_message, pfo->error_message);

  }
  else {

    /* alternative method called Pk_equal, described in 0810.0190 and
       1601.07230, extending the range of validity of
       HALOFIT from constant w to (w0,wa) models. In that
       case, some effective values of w(tau_i) and
       Omega_m(tau_i) have been pre-computed in the
       input module, and we just ned to interpolate
       within tabulated arrays, to get them at the
       current tau value. */

    class_alloc(w_and_Omega,pfo->pk_eq_size*sizeof(double),pfo->error_message);

    class_call(array_interpolate_spline(
                                        pfo->pk_eq_tau,
                                        pfo->pk_eq_tau_size,
                                        pfo->pk_eq_w_and_Omega,
                                        pfo->pk_eq_ddw_and_ddOmega,
                                        pfo->pk_eq_size,
                                        tau,
                                        &last_index,
                                        w_and_Omega,
                                        pfo->pk_eq_size,
                                        pfo->error_message),
               pfo->error_message,
               pfo->error_message);

    w = w_and_Omega[pfo->index_pk_eq_w];
    Omega_m = w_and_Omega[pfo->index_pk_eq_Omega_m];
    Omega_v = 1.-Omega_m;

    free(w_and_Omega);
  }

  anorm    = 1./(2*pow(_PI_,2));

  /*      Until the 17.02.2015 the values of k used for integrating sigma(R) quantities needed by Halofit where the same as in the perturbation module.
          Since then, we sample these integrals on more values, in order to get more precise integrals (thanks Matteo Zennaro for noticing the need for this).

          We create a temporary integrand_array which columns will be:
          - k in 1/Mpc
          - just linear P(k) in Mpc**3
          - 1/(2(pi**2)) P(k) k**2 exp(-(kR)**2) or 1/(2(pi**2)) P(k) k**2 2 (kR) exp(-(kR)**2) or 1/(2(pi**2)) P(k) k**2 4 (kR)(1-kR) exp(-(kR)**2)
          - second derivative of previous line with spline
  */

  index_ia=0;
  class_define_index(index_ia_k,     _TRUE_,index_ia,1);
  class_define_index(index_ia_pk,    _TRUE_,index_ia,1);
  class_define_index(index_ia_sum,   _TRUE_,index_ia,1);
  class_define_index(index_ia_ddsum, _TRUE_,index_ia,1);
  ia_size = index_ia;

  integrand_size=(int)(log(pfo->k[pfo->k_size-1]/pfo->k[0])/log(10.)*ppr->halofit_k_per_decade)+1;

  class_alloc(integrand_array,integrand_size*ia_size*sizeof(double),pfo->error_message);


  /* we fill integrand_array with values of k and P(k) using interpolation */

  last_index=0;

  for (index_k=0; index_k < integrand_size; index_k++) {

    k_integrand=pfo->k[0]*pow(10.,index_k/ppr->halofit_k_per_decade);

    if (index_k ==0 ) {
      lnpk_integrand = lnpk_l[0];
    }
    else {

      class_call(array_interpolate_spline(
                                          pfo->ln_k,
                                          pfo->k_size,
                                          lnpk_l,
                                          ddlnpk_l,
                                          1,
                                          log(k_integrand),
                                          &last_index,
                                          &lnpk_integrand,
                                          1,
                                          pfo->error_message),
                 pfo->error_message,
                 pfo->error_message);
    }

    integrand_array[index_k*ia_size + index_ia_k] = k_integrand;
    integrand_array[index_k*ia_size + index_ia_pk] = exp(lnpk_integrand);

  }

  class_call(background_at_tau(pba,tau,long_info,inter_normal,&last_index,pvecback),
             pba->error_message,
             pfo->error_message);

  Omega_m = pvecback[pba->index_bg_Omega_m];
  Omega_v = 1.-pvecback[pba->index_bg_Omega_m]-pvecback[pba->index_bg_Omega_r];

  // for debugging:
  //printf("Call Halofit at z=%e\n",1./pvecback[pba->index_bg_a]-1.);

  /* minimum value of R such that the integral giving sigma_R is
     converged.  The parameter halofit_sigma_precision should be
     understood as follows: we trust our calculation of sigma(R) as
     long as the integral reaches a value k_max such that the factor
     exp(-(Rk_max)**2) is already as low as halofit_sigma_precisio,
     shoing that the integreal is converged.  In practise this
     condition is tested only for R_max, the highest value of R in our
     bisection algorithm. Hence a smaller value of
     halofit_sigma_precision will lead to a more precise halofit
     result at the *highest* redshift at which halofit can make
     computations, at the expense of requiring a larger k_max; but
     this parameter is not relevant for the precision on P_nl(k,z) at
     other redshifts, so there is normally no need to change i
  */

  R=sqrt(-log(ppr->halofit_sigma_precision))/integrand_array[(integrand_size-1)*ia_size + index_ia_k];

  class_call(fourier_halofit_integrate(
                                       pfo,
                                       integrand_array,
                                       integrand_size,
                                       ia_size,
                                       index_ia_k,
                                       index_ia_pk,
                                       index_ia_sum,
                                       index_ia_ddsum,
                                       R,
                                       halofit_integral_one,
                                       &sum1
                                       ),
             pfo->error_message,
             pfo->error_message);

  sigma  = sqrt(sum1);

  /* the following error should not stop the code: it will arrive
     inevitably at some large redshift, and then the code should not
     stop, but just give up computing P_NL(k,z). This is why we have a
     special error handling here (using class_test_except and free()
     commands to avoid memory leaks, and calling this whole function
     not through a class_call) */

  /*
    class_test_except(sigma < 1.,
    pfo->error_message,
    free(pvecback);free(integrand_array),
    "Your k_max=%g 1/Mpc is too small for Halofit to find the non-linearity scale z_nl at z=%g. Increase input parameter P_k_max_h/Mpc or P_k_max_1/Mpc",
    pfo->k[pfo->k_size-1],
    1./pvecback[pba->index_bg_a]-1.);
  */

  if (sigma < 1.) {
    * nl_corr_not_computable_at_this_k = _TRUE_;
    free(pvecback);
    free(integrand_array);
    return _SUCCESS_;
  }
  else {
    * nl_corr_not_computable_at_this_k = _FALSE_;
  }

  xlogr1 = log(R)/log(10.);

  /* maximum value of R in the bisection algorithm leading to the
     determination of R_nl.  For this value we can make a
     conservaitive guess: 1/halofit_min_k_nonlinear, where
     halofit_min_k_nonlinear is the minimum value of k at which we ask
     halofit to give us an estimate of P_nl(k,z). By assumption we
     treat all smaller k's as linear, so we know that
     sigma(1/halofit_min_k_nonlinear) must be <<1 (and if it is not
     the test below will alert us) */

  R=1./ppr->halofit_min_k_nonlinear;

  /* corresponding value of sigma_R */
  class_call(fourier_halofit_integrate(
                                       pfo,
                                       integrand_array,
                                       integrand_size,
                                       ia_size,
                                       index_ia_k,
                                       index_ia_pk,
                                       index_ia_sum,
                                       index_ia_ddsum,
                                       R,
                                       halofit_integral_one,
                                       &sum1
                                       ),
             pfo->error_message,
             pfo->error_message);

  sigma  = sqrt(sum1);

  class_test(sigma > 1.,
             pfo->error_message,
             "Your input value for the precision parameter halofit_min_k_nonlinear=%e is too large, such that sigma(R=1/halofit_min_k_nonlinear)=% > 1. For self-consistency, it should have been <1. Decrease halofit_min_k_nonlinear",
             ppr->halofit_min_k_nonlinear,sigma);

  xlogr2 = log(R)/log(10.);

  counter = 0;
  do {
    rmid = pow(10,(xlogr2+xlogr1)/2.0);
    counter ++;

    class_call(fourier_halofit_integrate(
                                         pfo,
                                         integrand_array,
                                         integrand_size,
                                         ia_size,
                                         index_ia_k,
                                         index_ia_pk,
                                         index_ia_sum,
                                         index_ia_ddsum,
                                         rmid,
                                         halofit_integral_one,
                                         &sum1
                                         ),
               pfo->error_message,
               pfo->error_message);

    sigma  = sqrt(sum1);

    diff = sigma - 1.0;

    if (diff > ppr->halofit_tol_sigma){
      xlogr1=log10(rmid);
    }
    else if (diff < -ppr->halofit_tol_sigma) {
      xlogr2 = log10(rmid);
    }

    /* The first version of this test woukld let the code continue: */
    /*
      class_test_except(counter > _MAX_IT_,
      pfo->error_message,
      free(pvecback);free(integrand_array),
      "could not converge within maximum allowed number of iterations");
    */
    /* ... but in this situation it sounds better to make it stop and return an error! */
    class_test(counter > _MAX_IT_,
               pfo->error_message,
               "could not converge within maximum allowed number of iterations");

  } while (fabs(diff) > ppr->halofit_tol_sigma);

  /* evaluate all the other integrals at R=rmid */

  class_call(fourier_halofit_integrate(
                                       pfo,
                                       integrand_array,
                                       integrand_size,
                                       ia_size,
                                       index_ia_k,
                                       index_ia_pk,
                                       index_ia_sum,
                                       index_ia_ddsum,
                                       rmid,
                                       halofit_integral_two,
                                       &sum2
                                       ),
             pfo->error_message,
             pfo->error_message);

  class_call(fourier_halofit_integrate(
                                       pfo,
                                       integrand_array,
                                       integrand_size,
                                       ia_size,
                                       index_ia_k,
                                       index_ia_pk,
                                       index_ia_sum,
                                       index_ia_ddsum,
                                       rmid,
                                       halofit_integral_three,
                                       &sum3
                                       ),
             pfo->error_message,
             pfo->error_message);

  sigma  = sqrt(sum1);
  d1 = -sum2/sum1;
  d2 = -sum2*sum2/sum1/sum1 - sum3/sum1;

  rknl  = 1./rmid;
  rneff = -3.-d1;
  rncur = -d2;

  *k_nl = rknl;

  for (index_k = 0; index_k < pfo->k_size; index_k++){

    rk = pfo->k[index_k];

    if (rk > ppr->halofit_min_k_nonlinear) {

      pk_lin = exp(lnpk_l[index_k])*pow(pfo->k[index_k],3)*anorm;

      /* in original halofit, this is the beginning of the function halofit() */

      /*SPB11: Standard halofit underestimates the power on the smallest
       * scales by a factor of two. Add an extra correction from the
       * simulations in Bird, Viel,Haehnelt 2011 which partially accounts for
       * this.*/
      /*SPB14: This version of halofit is an updated version of the fit to the massive neutrinos
       * based on the results of Takahashi 2012, (arXiv:1208.2701).
       */
      gam=0.1971-0.0843*rneff+0.8460*rncur;
      a=1.5222+2.8553*rneff+2.3706*rneff*rneff+0.9903*rneff*rneff*rneff+ 0.2250*rneff*rneff*rneff*rneff-0.6038*rncur+0.1749*Omega_v*(1.+w);
      a=pow(10,a);
      b=pow(10, (-0.5642+0.5864*rneff+0.5716*rneff*rneff-1.5474*rncur+0.2279*Omega_v*(1.+w)));
      c=pow(10, 0.3698+2.0404*rneff+0.8161*rneff*rneff+0.5869*rncur);
      xmu=0.;
      xnu=pow(10,5.2105+3.6902*rneff);
      alpha=fabs(6.0835+1.3373*rneff-0.1959*rneff*rneff-5.5274*rncur);
      beta=2.0379-0.7354*rneff+0.3157*pow(rneff,2)+1.2490*pow(rneff,3)+0.3980*pow(rneff,4)-0.1682*rncur + fnu*(1.081 + 0.395*pow(rneff,2));

      if (fabs(1-Omega_m)>0.01) { /*then omega evolution */
        f1a=pow(Omega_m,(-0.0732));
        f2a=pow(Omega_m,(-0.1423));
        f3a=pow(Omega_m,(0.0725));
        f1b=pow(Omega_m,(-0.0307));
        f2b=pow(Omega_m,(-0.0585));
        f3b=pow(Omega_m,(0.0743));
        frac=Omega_v/(1.-Omega_m);
        f1=frac*f1b + (1-frac)*f1a;
        f2=frac*f2b + (1-frac)*f2a;
        f3=frac*f3b + (1-frac)*f3a;
      }
      else {
        f1=1.;
        f2=1.;
        f3=1.;
      }

      y=(rk/rknl);
      pk_halo = a*pow(y,f1*3.)/(1.+b*pow(y,f2)+pow(f3*c*y,3.-gam));
      pk_halo=pk_halo/(1+xmu*pow(y,-1)+xnu*pow(y,-2))*(1+fnu*0.977);

      /* until v2.9.3 pk_halo did contain an additional correction
         coming from Simeon Bird: the last factor was
         (1+fnu*(0.977-18.015*(pba->Omega0_m-0.3))). It seems that Bird
         gave it up later in his CAMB implementation and thus we also
         removed it. */
      // rk is in 1/Mpc, 47.48and 1.5 in Mpc**-2, so we need an h**2 here (Credits Antonio J. Cuesta)
      pk_linaa=pk_lin*(1+fnu*47.48*pow(rk/pba->h,2)/(1+1.5*pow(rk/pba->h,2)));
      pk_quasi=pk_lin*pow((1+pk_linaa),beta)/(1+pk_linaa*alpha)*exp(-y/4.0-pow(y,2)/8.0);

      pk_nl[index_k] = (pk_halo+pk_quasi)/pow(pfo->k[index_k],3)/anorm;

      /* in original halofit, this is the end of the function halofit() */
    }
    else {
      pk_nl[index_k] = exp(lnpk_l[index_k]);
    }
  }

  free(pvecback);
  free(integrand_array);
  return _SUCCESS_;
}

/**
 * Internal routione of Halofit. In original Halofit, this is
 * equivalent to the function wint(). It performs convolutions of the
 * linear spectrum with two window functions.
 *
 * @param pfo             Input: pointer to non linear structure
 * @param integrand_array Input: array with k, P_L(k) values
 * @param integrand_size  Input: one dimension of that array
 * @param ia_size         Input: other dimension of that array
 * @param index_ia_k      Input: index for k
 * @param index_ia_pk     Input: index for pk
 * @param index_ia_sum    Input: index for the result
 * @param index_ia_ddsum  Input: index for its spline
 * @param R               Input: radius
 * @param type            Input: which window function to use
 * @param sum             Output: result of the integral
 * @return the error status
 */

int fourier_halofit_integrate(
                              struct fourier *pfo,
                              double * integrand_array,
                              int integrand_size,
                              int ia_size,
                              int index_ia_k,
                              int index_ia_pk,
                              int index_ia_sum,
                              int index_ia_ddsum,
                              double R,
                              enum halofit_integral_type type,
                              double * sum
                              ) {

  double k,pk,x2,integrand;
  int index_k;
  double anorm = 1./(2*pow(_PI_,2));

  for (index_k=0; index_k < integrand_size; index_k++) {
    k = integrand_array[index_k*ia_size + index_ia_k];
    pk = integrand_array[index_k*ia_size + index_ia_pk];
    x2 = k*k*R*R;

    integrand = pk*k*k*anorm*exp(-x2);
    if (type == halofit_integral_two) integrand *= 2.*x2;
    if (type == halofit_integral_three) integrand *= 4.*x2*(1.-x2);

    integrand_array[index_k*ia_size + index_ia_sum] = integrand;
  }

  /* fill in second derivatives */
  class_call(array_spline(integrand_array,
                          ia_size,
                          integrand_size,
                          index_ia_k,
                          index_ia_sum,
                          index_ia_ddsum,
                          _SPLINE_NATURAL_,
                          pfo->error_message),
             pfo->error_message,
             pfo->error_message);

  /* integrate */
  class_call(array_integrate_all_spline(integrand_array,
                                        ia_size,
                                        integrand_size,
                                        index_ia_k,
                                        index_ia_sum,
                                        index_ia_ddsum,
                                        sum,
                                        pfo->error_message),
             pfo->error_message,
             pfo->error_message);

  return _SUCCESS_;
}

/**
 * Computes the nonlinear correction on the linear power spectrum via
 * the method presented in Mead et al. 1505.07833
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param ppt Input: pointer to perturbation structure
 * @param ppm Input: pointer to primordial structure
 * @param pfo Input: pointer to fourier structure
 * @param index_pk   Input: index of the pk type, either index_m or index_cb
 * @param index_tau  Input: index of tau, at which to compute the nl correction
 * @param tau        Input: tau, at which to compute the nl correction
 * @param pk_nl      Output:nonlinear power spectrum
 * @param lnpk_l     Input: logarithm of the linear power spectrum for both index_m and index_cb
 * @param ddlnpk_l   Input: spline of the logarithm of the linear power spectrum for both index_m and index_cb
 * @param nl_corr_not_computable_at_this_k Ouput: was the computation doable?
 * @param k_nl       Output: nonlinear scale for index_m and index_cb
 * @param pnw        Input/Output: pointer to nonlinear workspace
 * @return the error status
 */

int fourier_hmcode(
                   struct precision *ppr,
                   struct background *pba,
                   struct perturbations *ppt,
                   struct primordial *ppm,
                   struct fourier *pfo,
                   int index_pk,
                   int index_tau,
                   double tau,
                   double *pk_nl,
                   double **lnpk_l,
                   double **ddlnpk_l,
                   double *k_nl,
                   short * nl_corr_not_computable_at_this_k,
                   struct fourier_workspace * pnw
                   ) {

  /* integers */
  int index_mass, i, ng, nsig;
  int index_k, index_ncol;
  int last_index=0;
  int index_pk_cb;
  int counter, index_nl;

  int index_nu, index_cut;
  int index_y;
  int index_ddy;

  /* Background parameters */
  double Omega_m,fnu,Omega0_m;
  double z_at_tau;
  double rho_crit_today_in_msun_mpc3;
  double growth;
  double anorm;

  /* temporary numbers */
  double m, r, nu, sig, sigf;
  double diff, r1, r2;

  /* HMcode parameters */
  double mmin, mmax, nu_min;

  double sigma_disp, sigma_disp100, sigma8;
  double delta_c, Delta_v;
  double fraction;

  double sigma_nl, nu_nl, r_nl;
  double sigma_prime;
  double dlnsigdlnR;
  double n_eff;
  double alpha;

  double z_form, g_form;

  double eta;
  double gst, window_nfw;
  double nu_cut;
  double fac, k_star, fdamp;
  double pk_lin, pk_2h, pk_1h;

  /* data fields */
  double * pvecback;
  double * conc;
  double * mass;
  double * sigma_r;
  double * sigmaf_r;
  double * r_virial;
  double * r_real;
  double * nu_arr;

  double * p1h_integrand;


  /** include precision parameters that control the number of entries in the growth and sigma tables */
  ng = ppr->n_hmcode_tables;
  nsig = ppr->n_hmcode_tables;

  /** Compute background quantitites today */

  Omega0_m = pba->Omega0_m;
  fnu      = pba->Omega0_ncdm_tot/Omega0_m;

  /** If index_pk_cb, choose Omega0_cb as the matter density parameter.
   * If index_pk_m, choose Omega0_cbn as the matter density parameter. */
  if (index_pk==pfo->index_pk_cb){
    Omega0_m = Omega0_m - pba->Omega0_ncdm_tot;
  }

  anorm    = 1./(2*pow(_PI_,2));

  /** Call all the relevant background parameters at this tau */
  class_alloc(pvecback,pba->bg_size*sizeof(double),pfo->error_message);

  class_call(background_at_tau(pba,tau,long_info,inter_normal,&last_index,pvecback),
             pba->error_message,
             pfo->error_message);

  Omega_m = pvecback[pba->index_bg_Omega_m];//TBC (i.e. check if for P_cb here we should use Omega_cb) here the total time varying Omega_m is used for delta_c and for Delta_v according to the Mead fit of the Massara simulations.

  growth = pvecback[pba->index_bg_D];

  z_at_tau = 1./pvecback[pba->index_bg_a]-1.;

  /* The number below is the critical density today, rho_c = 3 * H0^2 / 8*pi*G, in units of M_sun over Mpc^3 */
  rho_crit_today_in_msun_mpc3 = 3.*pow(1.e5*pba->h, 2)/8./_PI_/_G_*_Mpc_over_m_/_M_SUN_;

  free(pvecback);

  /** Test whether pk_cb has to be taken into account (only if we have massive neutrinos)*/
  if (pba->has_ncdm==_TRUE_){
    index_pk_cb = pfo->index_pk_cb;
  }
  else {
    index_pk_cb = index_pk;
  }


  /** Get sigma(R=8 Mpc/h), sigma_disp(R=0), sigma_disp(R=100 Mpc/h) and write them into pfo structure */

  class_call(fourier_sigmas(pfo,
                            8./pba->h,
                            lnpk_l[index_pk],ddlnpk_l[index_pk],
                            pfo->k_size_extra,
                            ppr->sigma_k_per_decade,
                            out_sigma,
                            &sigma8),
             pfo->error_message,
             pfo->error_message);

  class_call(fourier_sigmas(pfo,
                            0.,
                            lnpk_l[index_pk],ddlnpk_l[index_pk],
                            pfo->k_size_extra,
                            ppr->sigma_k_per_decade,
                            out_sigma_disp,
                            &sigma_disp),
             pfo->error_message,
             pfo->error_message);

  class_call(fourier_sigmas(pfo,
                            100./pba->h,
                            lnpk_l[index_pk],ddlnpk_l[index_pk],
                            pfo->k_size_extra,
                            ppr->sigma_k_per_decade,
                            out_sigma_disp,
                            &sigma_disp100),
             pfo->error_message,
             pfo->error_message);

  pnw->sigma_8[index_pk][index_tau] = sigma8;
  pnw->sigma_disp[index_pk][index_tau] = sigma_disp;
  pnw->sigma_disp_100[index_pk][index_tau] = sigma_disp100;

  /** Initialisation steps for the 1-Halo Power Integral */
  mmin=ppr->mmin_for_p1h_integral/pba->h; //Minimum mass for integration; (unit conversion from  m[Msun/h] to m[Msun]  )
  mmax=ppr->mmax_for_p1h_integral/pba->h; //Maximum mass for integration;

  class_alloc(mass,ppr->nsteps_for_p1h_integral*sizeof(double),pfo->error_message);
  class_alloc(r_real,ppr->nsteps_for_p1h_integral*sizeof(double),pfo->error_message);
  class_alloc(r_virial,ppr->nsteps_for_p1h_integral*sizeof(double),pfo->error_message);
  class_alloc(sigma_r,ppr->nsteps_for_p1h_integral*sizeof(double),pfo->error_message);
  class_alloc(sigmaf_r,ppr->nsteps_for_p1h_integral*sizeof(double),pfo->error_message);
  class_alloc(nu_arr,ppr->nsteps_for_p1h_integral*sizeof(double),pfo->error_message);

  // Linear theory density perturbation threshold for spherical collapse
  delta_c = 1.59+0.0314*log(sigma8); //Mead et al. (2015; arXiv 1505.07833)
  delta_c = delta_c*(1.+0.0123*log10(Omega_m)); //Nakamura & Suto (1997) fitting formula for LCDM models (as in Mead 2016)
  delta_c = delta_c*(1.+0.262*fnu); //Mead et al. (2016; arXiv 1602.02154) neutrino addition

  // virialized overdensity
  Delta_v=418.*pow(Omega_m, -0.352); //Mead et al. (2015; arXiv 1505.07833)
  Delta_v=Delta_v*(1.+0.916*fnu); //Mead et al. (2016; arXiv 1602.02154) neutrino addition

  // mass or radius fraction respectively
  fraction = pow(0.01, 1./3.);

  /* Fill the arrays needed for the P1H Integral: mass, r_real, r_virial, nu_arr, sigma_r, sigmaf_r
   * The P1H Integral is an integral over nu=delta_c/sigma(M), where M is connected to R via R=(3M)/(4*pi*rho_m).
   * The Integrand is M*Window^2{nu(M)*k, Rv(M), c(M)}*f(nu) with the window being the fouriertransformed
   * NFW profile, Rv = R/Delta_v^(1/3) and Sheth-Thormen halo mass function f.
   * The halo concentration-mass-relation c(M) will be found later.  */

  for (index_mass=0;index_mass<ppr->nsteps_for_p1h_integral;index_mass++){

    m = exp(log(mmin)+log(mmax/mmin)*(index_mass)/(ppr->nsteps_for_p1h_integral-1));
    r = pow((3.*m/(4.*_PI_*rho_crit_today_in_msun_mpc3*Omega0_m)), (1./3.));
    mass[index_mass] = m;
    r_real[index_mass] = r;
    r_virial[index_mass] = r_real[index_mass]/pow(Delta_v, 1./3.);

    class_call(array_interpolate_spline(pnw->rtab,
                                        nsig,
                                        pnw->stab,
                                        pnw->ddstab,
                                        1,
                                        r,
                                        &last_index,
                                        &sig,
                                        1,
                                        pfo->error_message),
               pfo->error_message, pfo->error_message);

    class_call(array_interpolate_spline(pnw->rtab,
                                        nsig,
                                        pnw->stab,
                                        pnw->ddstab,
                                        1,
                                        r*fraction,
                                        &last_index,
                                        &sigf,
                                        1,
                                        pfo->error_message),
               pfo->error_message, pfo->error_message);

    nu=delta_c/sig;
    sigma_r[index_mass] = sig;
    sigmaf_r[index_mass] = sigf;
    nu_arr[index_mass] = nu;
  }

  /** find nonlinear scales k_nl and r_nl and the effective spectral index n_eff */
  nu_nl = 1.;
  nu_min = nu_arr[0];

  /* stop calculating the nonlinear correction if the nonlinear scale is not reached in the table: */
  if (nu_min > nu_nl) {
    if (pfo->fourier_verbose>0) fprintf(stdout, " -> [WARNING:] the minimum mass in the mass-table is too large to find the nonlinear scale at this redshift.\n   Decrease mmin_for_p1h_integral\n");
    * nl_corr_not_computable_at_this_k = _TRUE_;
    free(mass);
    free(r_real);
    free(r_virial);
    free(sigma_r);
    free(sigmaf_r);
    free(nu_arr);
    return _SUCCESS_;
  }

  /* make a first guess for the nonlinear scale */
  class_call(array_interpolate_two_arrays_one_column(
                                                     nu_arr,
                                                     r_real,
                                                     1,
                                                     0,
                                                     ppr->nsteps_for_p1h_integral,
                                                     nu_nl,
                                                     &r_nl,
                                                     pfo->error_message),
             pfo->error_message, pfo->error_message);

  class_call(array_search_bisect(ppr->nsteps_for_p1h_integral,nu_arr,nu_nl,&index_nl,pfo->error_message), pfo->error_message, pfo->error_message);

  r1 = r_real[index_nl-1];
  r2 = r_real[index_nl+2];

  /* // for debugging: (if it happens that r_nl is not between r1 and r2, which should never be the case)
     fprintf(stdout, "%e %e %e %e\n", r1, nu_arr[index_nl-1], r2, nu_arr[index_nl+2]);
  */

  /* do bisectional iteration between r1 and r2 to find the precise value of r_nl */
  counter = 0;
  do {
    r_nl = (r1+r2)/2.;
    counter ++;

    class_call(fourier_sigmas(pfo,
                              r_nl,
                              lnpk_l[index_pk_cb],ddlnpk_l[index_pk_cb],
                              pfo->k_size_extra,
                              ppr->sigma_k_per_decade,
                              out_sigma,
                              &sigma_nl),
               pfo->error_message, pfo->error_message);

    diff = sigma_nl - delta_c;

    if (diff > ppr->hmcode_tol_sigma){
      r1=r_nl;
    }
    else if (diff < -ppr->hmcode_tol_sigma) {
      r2 = r_nl;
    }

    class_test(counter > _MAX_IT_,
               pfo->error_message,
               "could not converge within maximum allowed number of iterations");

  } while (fabs(diff) > ppr->hmcode_tol_sigma);

  if (pfo->fourier_verbose>5){
    fprintf(stdout, "number of iterations for r_nl at z = %e: %d\n", z_at_tau, counter);
  }
  *k_nl = 1./r_nl;

  if (*k_nl > pfo->k[pfo->k_size-1]) {
    * nl_corr_not_computable_at_this_k = _TRUE_;
    free(mass);
    free(r_real);
    free(r_virial);
    free(sigma_r);
    free(sigmaf_r);
    free(nu_arr);
    return _SUCCESS_;
  }
  else {
    * nl_corr_not_computable_at_this_k = _FALSE_;
  }

  /* call sigma_prime function at r_nl to find the effective spectral index n_eff */

  class_call(fourier_sigmas(pfo,
                            r_nl,
                            lnpk_l[index_pk_cb],ddlnpk_l[index_pk_cb],
                            pfo->k_size_extra,
                            ppr->sigma_k_per_decade,
                            out_sigma_prime,
                            &sigma_prime),
             pfo->error_message,
             pfo->error_message);

  dlnsigdlnR = r_nl*pow(sigma_nl, -2)*sigma_prime;
  n_eff = -3.- dlnsigdlnR;
  alpha = 3.24*pow(1.85, n_eff);

  pnw->sigma_prime[index_pk][index_tau] = sigma_prime;

  /** Calculate halo concentration-mass relation conc(mass) (Bullock et al. 2001) */
  class_alloc(conc,ppr->nsteps_for_p1h_integral*sizeof(double),pfo->error_message);

  for (index_mass=0;index_mass<ppr->nsteps_for_p1h_integral;index_mass++){
    //find growth rate at formation
    g_form = delta_c*growth/sigmaf_r[index_mass];
    if (g_form > 1.) g_form = 1.;

    //
    class_call(array_interpolate_two_arrays_one_column(
                                                       pnw->growtable,
                                                       pnw->ztable,
                                                       1,
                                                       0,
                                                       ng,
                                                       g_form,
                                                       &z_form,
                                                       pfo->error_message),
               pfo->error_message, pfo->error_message);
    if (z_form < z_at_tau){
      conc[index_mass] = pfo->c_min;
    } else {
      conc[index_mass] = pfo->c_min*(1.+z_form)/(1.+z_at_tau)*pnw->dark_energy_correction;
    }
  }


  /** Compute the nonlinear correction */
  eta = pfo->eta_0 - 0.3*sigma8; // halo bloating parameter
  k_star=0.584/sigma_disp;   // Damping wavenumber of the 1-halo term at very large scales;
  fdamp = 0.0095*pow(sigma_disp100*pba->h, 1.37); // Damping factor for 2-halo term
  if (fdamp<1.e-3) fdamp=1.e-3;
  if (fdamp>0.99)  fdamp=0.99;

  /* the 1h integral contains the halo mass function proportional to exp(-nu^2).
   * To save time, the integration loop cuts, when nu exceeds a large value,
   * where the integrand is 0 anyhow. This cut index is found here. */
  nu_cut = 10.;
  if (nu_cut < nu_arr[ppr->nsteps_for_p1h_integral-1]){
    class_call(array_search_bisect(ppr->nsteps_for_p1h_integral,nu_arr,nu_cut,&index_cut,pfo->error_message), pfo->error_message, pfo->error_message);
  }
  else {
    index_cut = ppr->nsteps_for_p1h_integral;
  }

  i=0;
  index_nu=i;
  i++;
  index_y=i;
  i++;
  index_ddy=i;
  i++;
  index_ncol=i;

  for (index_k = 0; index_k < pfo->k_size; index_k++){

    class_alloc(p1h_integrand,index_cut*index_ncol*sizeof(double),pfo->error_message);

    pk_lin = exp(lnpk_l[index_pk][index_k])*pow(pfo->k[index_k],3)*anorm; //convert P_k to Delta_k^2

    for (index_mass=0; index_mass<index_cut; index_mass++){ //Calculates the integrand for the ph1 integral at all nu values
      //get the nu^eta-value of the window
      class_call(fourier_hmcode_window_nfw(
                                           pfo,
                                           pow(nu_arr[index_mass], eta)*pfo->k[index_k],
                                           r_virial[index_mass],
                                           conc[index_mass],
                                           &window_nfw),
                 pfo->error_message, pfo->error_message);
      //get the value of the halo mass function
      class_call(fourier_hmcode_halomassfunction(
                                                 nu_arr[index_mass],
                                                 &gst),
                 pfo->error_message, pfo->error_message);

      p1h_integrand[index_mass*index_ncol+index_nu] = nu_arr[index_mass];

      p1h_integrand[index_mass*index_ncol+index_y] = mass[index_mass]*gst*pow(window_nfw, 2.);
      //if ((tau==pba->conformal_age) && (index_k == 0)) {
      //fprintf(stdout, "%d %e %e\n", index_cut, p1h_integrand[index_mass*index_ncol+index_nu], p1h_integrand[index_mass*index_ncol+index_y]);
      //}
    }
    class_call(array_spline(p1h_integrand,
                            index_ncol,
                            index_cut,
                            index_nu,
                            index_y,
                            index_ddy,
                            _SPLINE_EST_DERIV_,
                            pfo->error_message),
               pfo->error_message,
               pfo->error_message);

    class_call(array_integrate_all_trapzd_or_spline(
                                                    p1h_integrand,
                                                    index_ncol,
                                                    index_cut,
                                                    index_cut-1, //0 or n-1
                                                    index_nu,
                                                    index_y,
                                                    index_ddy,
                                                    &pk_1h,
                                                    pfo->error_message),
               pfo->error_message,
               pfo->error_message);


    if (pow(pfo->k[index_k]/k_star, 2)>7.){
      fac = 0.;     //prevents problems if (k/k*)^2 is large
    }
    else{
      fac=exp(-pow((pfo->k[index_k]/k_star), 2.));
    }

    pk_1h = pk_1h*anorm*pow(pfo->k[index_k],3)*(1.-fac)/(rho_crit_today_in_msun_mpc3*Omega0_m);  // dimensionless power

    if (fdamp==0){
      pk_2h=pk_lin;
    }else{
      pk_2h=pk_lin*(1.-fdamp*pow(tanh(pfo->k[index_k]*sigma_disp/sqrt(fdamp)), 2.)); //dimensionless power
    }
    if (pk_2h<0.) pk_2h=0.;
    pk_nl[index_k] = pow((pow(pk_1h, alpha) + pow(pk_2h, alpha)), (1./alpha))/pow(pfo->k[index_k],3)/anorm; //converted back to P_k

    free(p1h_integrand);
  }

  // print parameter values
  if ((pfo->fourier_verbose > 1 && tau==pba->conformal_age) || pfo->fourier_verbose > 3){
    fprintf(stdout, " -> Parameters at redshift z = %e:\n", z_at_tau);
    fprintf(stdout, "    fnu:		%e\n", fnu);
    fprintf(stdout, "    sigd [Mpc/h]:	%e\n", sigma_disp*pba->h);
    fprintf(stdout, "    sigd100 [Mpc/h]:    %e\n", sigma_disp100*pba->h);
    fprintf(stdout, "    sigma8:		%e\n", sigma8);
    fprintf(stdout, "    nu min:		%e\n", nu_arr[0]);
    fprintf(stdout, "    nu max:		%e\n", nu_arr[ppr->nsteps_for_p1h_integral-1]);
    fprintf(stdout, "    r_v min [Mpc/h]:    %e\n", r_virial[0]*pba->h);
    fprintf(stdout, "    r_v max [Mpc/h]:    %e\n", r_virial[ppr->nsteps_for_p1h_integral-1]*pba->h);
    fprintf(stdout, "    r_nl [Mpc/h]:	%e\n", r_nl*pba->h);
    fprintf(stdout, "    k_nl [h/Mpc]:	%e\n", *k_nl/pba->h);
    fprintf(stdout, "    sigma_nl:		%e\n", sigma_nl/delta_c);
    fprintf(stdout, "    neff:		%e\n", n_eff);
    fprintf(stdout, "    c min:		%e\n", conc[ppr->nsteps_for_p1h_integral-1]);
    fprintf(stdout, "    c max:		%e\n", conc[0]);
    fprintf(stdout, "    Dv:			%e\n", Delta_v);
    fprintf(stdout, "    dc:			%e\n", delta_c);
    fprintf(stdout, "    eta:		%e\n", eta);
    fprintf(stdout, "    k*:			%e\n", k_star/pba->h);
    fprintf(stdout, "    Abary:		%e\n", pfo->c_min);
    fprintf(stdout, "    fdamp:		%e\n", fdamp);
    fprintf(stdout, "    alpha:		%e\n", alpha);
    fprintf(stdout, "    ksize, kmin, kmax:   %d, %e, %e\n", pfo->k_size, pfo->k[0]/pba->h, pfo->k[pfo->k_size-1]/pba->h);

  }

  free(conc);
  free(mass);
  free(r_real);
  free(r_virial);
  free(sigma_r);
  free(sigmaf_r);
  free(nu_arr);

  return _SUCCESS_;
}

/**
 * allocate and fill arrays of nonlinear workspace (currently used only by HMcode)
 *
 * @param ppr         Input: pointer to precision structure
 * @param pba         Input: pointer to background structure
 * @param pfo         Input: pointer to fourier structure
 * @param pnw         Output: pointer to nonlinear workspace
 * @return the error status
 */

int fourier_hmcode_workspace_init(
                                  struct precision *ppr,
                                  struct background *pba,
                                  struct fourier *pfo,
                                  struct fourier_workspace * pnw
                                  ){

  int ng;
  int index_pk;

  /** - allocate arrays of the nonlinear workspace */

  class_alloc(pnw->rtab,ppr->n_hmcode_tables*sizeof(double),pfo->error_message);
  class_alloc(pnw->stab,ppr->n_hmcode_tables*sizeof(double),pfo->error_message);
  class_alloc(pnw->ddstab,ppr->n_hmcode_tables*sizeof(double),pfo->error_message);

  ng = ppr->n_hmcode_tables;

  class_alloc(pnw->growtable,ng*sizeof(double),pfo->error_message);
  class_alloc(pnw->ztable,ng*sizeof(double),pfo->error_message);
  class_alloc(pnw->tautable,ng*sizeof(double),pfo->error_message);

  class_alloc(pnw->sigma_8,pfo->pk_size*sizeof(double *),pfo->error_message);
  class_alloc(pnw->sigma_disp,pfo->pk_size*sizeof(double *),pfo->error_message);
  class_alloc(pnw->sigma_disp_100,pfo->pk_size*sizeof(double *),pfo->error_message);
  class_alloc(pnw->sigma_prime,pfo->pk_size*sizeof(double *),pfo->error_message);

  for (index_pk=0; index_pk<pfo->pk_size; index_pk++){
    class_alloc(pnw->sigma_8[index_pk],pfo->tau_size*sizeof(double),pfo->error_message);
    class_alloc(pnw->sigma_disp[index_pk],pfo->tau_size*sizeof(double),pfo->error_message);
    class_alloc(pnw->sigma_disp_100[index_pk],pfo->tau_size*sizeof(double),pfo->error_message);
    class_alloc(pnw->sigma_prime[index_pk],pfo->tau_size*sizeof(double),pfo->error_message);
  }

  /** - fill table with scale independent growth factor */

  class_call(fourier_hmcode_fill_growtab(ppr,pba,pfo,pnw),
             pfo->error_message,
             pfo->error_message);

  return _SUCCESS_;
}

/**
 * deallocate arrays in the nonlinear worksapce (currently used only
 * by HMcode)
 *
 * @param pfo Input: pointer to fourier structure
 * @param pnw Input: pointer to nonlinear workspace
 * @return the error status
 */

int fourier_hmcode_workspace_free(
                                  struct fourier *pfo,
                                  struct fourier_workspace * pnw
                                  ) {
  int index_pk;

  free(pnw->rtab);
  free(pnw->stab);
  free(pnw->ddstab);

  free(pnw->growtable);
  free(pnw->ztable);
  free(pnw->tautable);

  for (index_pk=0; index_pk<pfo->pk_size; index_pk++){
    free(pnw->sigma_8[index_pk]);
    free(pnw->sigma_disp[index_pk]);
    free(pnw->sigma_disp_100[index_pk]);
    free(pnw->sigma_prime[index_pk]);
  }

  free(pnw->sigma_8);
  free(pnw->sigma_disp);
  free(pnw->sigma_disp_100);
  free(pnw->sigma_prime);

  return _SUCCESS_;
}

/**
 * set the HMcode dark energy correction (if w is not -1)
 *
 * @param ppr         Input: pointer to precision structure
 * @param pba         Input: pointer to background structure
 * @param pfo         Input: pointer to fourier structure
 * @param pnw         Output: pointer to nonlinear workspace
 * @return the error status
 */

int fourier_hmcode_dark_energy_correction(
                                          struct precision *ppr,
                                          struct background *pba,
                                          struct fourier *pfo,
                                          struct fourier_workspace * pnw
                                          ) {

  int last_index;
  double * pvecback;
  double tau_growth;
  double g_lcdm,g_wcdm;
  double w0,dw_over_da_fld,integral_fld;

  /** - if there is dynamical Dark Energy (w is not -1) modeled as a fluid */

  if (pba->has_fld==_TRUE_){

    class_alloc(pvecback,pba->bg_size*sizeof(double),pfo->error_message);

    class_call(background_tau_of_z(
                                   pba,
                                   pfo->z_infinity,
                                   &tau_growth
                                   ),
               pba->error_message,
               pfo->error_message);

    class_call(background_at_tau(pba,tau_growth,long_info,inter_normal,&last_index,pvecback),
               pba->error_message,
               pfo->error_message);

    class_call(background_w_fld(pba,1.,&w0,&dw_over_da_fld,&integral_fld),
               pba->error_message,
               pfo->error_message);

    class_call(fourier_hmcode_growint(ppr,pba,pfo,1./(1.+pfo->z_infinity),-1.,0.,&g_lcdm),
               pfo->error_message, pfo->error_message);

    class_call(fourier_hmcode_growint(ppr,pba,pfo,1./(1.+pfo->z_infinity),w0,dw_over_da_fld*(-1.),&g_wcdm),
               pfo->error_message,
               pfo->error_message);

    free(pvecback);

    pnw->dark_energy_correction = pow(g_wcdm/g_lcdm, 1.5);
  }

  /** - otherwise, we assume no dynamical Dark Energy (w is -1) */

  else {
    pnw->dark_energy_correction = 1.;
  }

  return _SUCCESS_;
}

/**
 * set the HMcode baryonic feedback parameters according to the chosen feedback model
 *
 * @param pfo   Output: pointer to fourier structure
 * @return the error status
 */

int fourier_hmcode_baryonic_feedback(
                                     struct fourier *pfo
                                     ) {

  switch (pfo->feedback) {

  case nl_emu_dmonly:
    {
      pfo->eta_0 = 0.603;
      pfo->c_min = 3.13;
      break;
    }

  case nl_owls_dmonly:
    {
      pfo->eta_0 = 0.64;
      pfo->c_min = 3.43;
      break;
    }

  case nl_owls_ref:
    {
      pfo->eta_0 = 0.68;
      pfo->c_min = 3.91;
      break;
    }

  case nl_owls_agn:
    {
      pfo->eta_0 = 0.76;
      pfo->c_min = 2.32;
      break;
    }

  case nl_owls_dblim:
    {
      pfo->eta_0 = 0.70;
      pfo->c_min = 3.01;
      break;
    }

  case nl_user_defined:
    {
      /* eta_0 and c_min already passed in input */
      break;
    }
  }
  return _SUCCESS_;
}

/**
 * Function that fills pnw->rtab, pnw->stab and pnw->ddstab with (r,
 * sigma, ddsigma) logarithmically spaced in r.  Called by
 * fourier_init at for all tau to account for scale-dependant growth
 * before fourier_hmcode is called
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param ppt Input: pointer to perturbation structure
 * @param ppm Input: pointer to primordial structure
 * @param pfo Input: pointer to fourier structure
 * @param index_tau  Input: index of tau, at which to compute the nl correction
 * @param lnpk_l   Input: logarithm of the linear power spectrum for either index_m or index_cb
 * @param ddlnpk_l Input: spline of the logarithm of the linear power spectrum for either index_m or index_cb
 * @param pnw Output: pointer to nonlinear workspace
 * @return the error status
 */

int fourier_hmcode_fill_sigtab(
                               struct precision * ppr,
                               struct background * pba,
                               struct perturbations * ppt,
                               struct primordial * ppm,
                               struct fourier * pfo,
                               int index_tau,
                               double *lnpk_l,
                               double *ddlnpk_l,
                               struct fourier_workspace * pnw
                               ) {

  double r;
  double rmin, rmax;
  double sig;
  double * sigtab;
  int i, index_r, index_sig, index_ddsig, index_n, nsig;

  rmin = ppr->rmin_for_sigtab/pba->h;
  rmax = ppr->rmax_for_sigtab/pba->h;
  nsig = ppr->n_hmcode_tables;

  i=0;
  index_r=i;
  i++;
  index_sig=i;
  i++;
  index_ddsig=i;
  i++;
  index_n=i;

  class_alloc((sigtab),(nsig*index_n*sizeof(double)),pfo->error_message);

  for (i=0;i<nsig;i++){
    r=exp(log(rmin)+log(rmax/rmin)*i/(nsig-1));

    class_call(fourier_sigmas(pfo,
                              r,
                              lnpk_l,
                              ddlnpk_l,
                              pfo->k_size_extra,
                              ppr->sigma_k_per_decade,
                              out_sigma,
                              &sig),
               pfo->error_message,
               pfo->error_message);

    sigtab[i*index_n+index_r]=r;
    sigtab[i*index_n+index_sig]=sig;
  }

  class_call(array_spline(sigtab,
						  index_n,
						  nsig,
						  index_r,
						  index_sig,
						  index_ddsig,
						  _SPLINE_EST_DERIV_,
						  pfo->error_message),
             pfo->error_message,
             pfo->error_message);
  if (index_tau == pfo->tau_size-1){
    for (i=0;i<nsig;i++){
      pnw->rtab[i] = sigtab[i*index_n+index_r];
      pnw->stab[i] = sigtab[i*index_n+index_sig];
      pnw->ddstab[i] = sigtab[i*index_n+index_ddsig];
    }
  }
  else{
    for (i=0;i<nsig;i++){
      pnw->stab[i] = sigtab[i*index_n+index_sig];
      pnw->ddstab[i] = sigtab[i*index_n+index_ddsig];
    }
  }

  free(sigtab);

  return _SUCCESS_;
}


/**
 * Function that fills pnw->tautable and pnw->growtable with (tau, D(tau))
 * linearly spaced in scalefactor a.
 * Called by fourier_init at before the loop over tau
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure (will provide the scale independent growth factor)
 * @param pfo Input/Output: pointer to fourier structure
 * @param pnw Output: pointer to nonlinear workspace
 * @return the error status
 */

int fourier_hmcode_fill_growtab(
                                struct precision * ppr,
                                struct background * pba,
                                struct fourier * pfo,
                                struct fourier_workspace * pnw
                                ){

  double z, ainit, amax, scalefactor, tau_growth;
  int index_scalefactor, last_index, ng;
  double * pvecback;

  ng = ppr->n_hmcode_tables;
  ainit = ppr->ainit_for_growtab;
  amax = ppr->amax_for_growtab;

  last_index = 0;

  class_alloc(pvecback,pba->bg_size*sizeof(double),pfo->error_message);

  for (index_scalefactor=0;index_scalefactor<ng;index_scalefactor++){
    scalefactor = ainit+(amax-ainit)*(index_scalefactor)/(ng-1);
    z = 1./scalefactor-1.;

    pnw->ztable[index_scalefactor] = z;

    class_call(background_tau_of_z(
                                   pba,
                                   z,
                                   &tau_growth
                                   ),
               pba->error_message, pfo->error_message);

    pnw->tautable[index_scalefactor] = tau_growth;

    class_call(background_at_tau(pba,tau_growth,long_info,inter_normal,&last_index,pvecback),
               pba->error_message,
               pfo->error_message);

    pnw->growtable[index_scalefactor] = pvecback[pba->index_bg_D];

  }

  free(pvecback);

  return _SUCCESS_;
}

/**
 * This function finds the scale independent growth factor by
 * integrating the approximate relation d(lnD)/d(lna) =
 * Omega_m(z)^gamma by Linder & Cahn 2007
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param pfo Input: pointer to fourier structure
 * @param a   Input: scalefactor
 * @param w0  Input: dark energy equation of state today
 * @param wa  Input: dark energy equation of state varying with a: w=w0+(1-a)wa
 * @param growth Output: scale independent growth factor at a
 * @return the error status
 */

int fourier_hmcode_growint(
                           struct precision * ppr,
                           struct background * pba,
                           struct fourier * pfo,
                           double a,
                           double w0,
                           double wa,
                           double * growth
                           ){

  double z, ainit, amax, scalefactor, gamma, X_de, Hubble2, Omega_m;
  int i, index_scalefactor, index_a, index_growth, index_ddgrowth, index_gcol, ng; // index_scalefactor is a running index while index_a is a column index
  double * pvecback;
  double * integrand;

  ng = 1024; // number of growth values (stepsize of the integral), should not be hardcoded and replaced by a precision parameter
  ainit = a;
  amax = 1.;

  i=0;
  index_a = i;
  i++;
  index_growth = i;
  i++;
  index_ddgrowth = i;
  i++;
  index_gcol = i;

  class_alloc(integrand,ng*index_gcol*sizeof(double),pfo->error_message);
  class_alloc(pvecback,pba->bg_size*sizeof(double),pfo->error_message);

  if (ainit == amax) {
    *growth = 1.;
  }
  else {

    for (index_scalefactor=0;index_scalefactor<ng;index_scalefactor++){

      scalefactor = ainit+(amax-ainit)*(index_scalefactor)/(ng-1);
      z = 1./scalefactor-1.;

      /* This will compute Omega_m(z) for the input values of w0 and wa, to let the user compare the wCDM and LCDM cases. This is why we cannot extract Omega_m(z) fromn the background module in this place. */
      X_de = pow(scalefactor, -3.*(1.+w0+wa))*exp(-3.*wa*(1.-scalefactor));
      Hubble2 = (pba->Omega0_m*pow((1.+z), 3.) + pba->Omega0_k*pow((1.+z), 2.) + pba->Omega0_de*X_de);
      Omega_m = (pba->Omega0_m*pow((1.+z), 3.))/Hubble2;
      /* Samuel brieden: TBC: check that the matching between the
         background quantity and this fitting formula improves by
         using Omega_cb (as it is done in background). Carefull:
         Hubble remains with Omega0_m */

      if (w0 == -1.){
        gamma = 0.55;
      }
      else if (w0 < -1.){
        gamma = 0.55+0.02*(1+w0);
      }
      else {
        gamma = 0.55+0.05*(1+w0);
      }
      integrand[index_scalefactor*index_gcol+index_a] = scalefactor;
      integrand[index_scalefactor*index_gcol+index_growth]= -pow(Omega_m, gamma)/scalefactor;
    }

    class_call(array_spline(integrand,
                            index_gcol,
                            ng,
                            index_a,
                            index_growth,
                            index_ddgrowth,
                            _SPLINE_EST_DERIV_,
                            pfo->error_message),
               pfo->error_message,
               pfo->error_message);

    class_call(array_integrate_all_trapzd_or_spline(integrand,
                                                    index_gcol,
                                                    ng,
                                                    0, //ng-1,
                                                    index_a,
                                                    index_growth,
                                                    index_ddgrowth,
                                                    growth,
                                                    pfo->error_message),
               pfo->error_message,
               pfo->error_message);

    *growth = exp(*growth);

  }
  //fprintf(stdout, "%e %e \n", a, *growth);
  free(pvecback);
  free(integrand);

  return _SUCCESS_;
}

/**
 * This is the fourier transform of the NFW density profile.
 *
 * @param pfo Input: pointer to fourier structure
 * @param k   Input: wave vector
 * @param rv  Input: virial radius
 * @param c   Input: concentration = rv/rs (with scale radius rs)
 * @param window_nfw Output: Window Function of the NFW profile
 * @return the error status
 */

int fourier_hmcode_window_nfw(
                              struct fourier * pfo,
                              double k,
                              double rv,
                              double c,
                              double *window_nfw
                              ){
  double si1, si2, ci1, ci2, ks;
  double p1, p2, p3;

  ks = k*rv/c;

  class_call(sine_integral(
                           ks*(1.+c),
                           &si2,
                           pfo->error_message
                           ),
             pfo->error_message, pfo->error_message);

  class_call(sine_integral(
                           ks,
                           &si1,
                           pfo->error_message
                           ),
             pfo->error_message, pfo->error_message);

  class_call(cosine_integral(
                             ks*(1.+c),
                             &ci2,
                             pfo->error_message
                             ),
             pfo->error_message, pfo->error_message);

  class_call(cosine_integral(
                             ks,
                             &ci1,
                             pfo->error_message
                             ),
             pfo->error_message, pfo->error_message);

  p1=cos(ks)*(ci2-ci1);
  p2=sin(ks)*(si2-si1);
  p3=sin(ks*c)/(ks*(1.+c));

  *window_nfw=p1+p2-p3;
  *window_nfw=*window_nfw/(log(1.+c)-c/(1.+c));

  return _SUCCESS_;
}

/**
 * This is the Sheth-Tormen halo mass function (1999, MNRAS, 308, 119)
 *
 * @param nu   Input: the \f$ \nu \f$ parameter that depends on the halo mass via \f$ \nu(M) = \delta_c/\sigma(M) \f$
 * @param hmf  Output: Value of the halo mass function at this \f$ \nu \f$
 * @return the error status
 */

int fourier_hmcode_halomassfunction(
                                    double nu,
                                    double * hmf
                                    ){

  double p, q, A;

  p=0.3;
  q=0.707;
  A=0.21616;

  *hmf=A*(1.+(pow(q*nu*nu, -p)))*exp(-q*nu*nu/2.);

  return _SUCCESS_;
}

/**
 * Compute sigma8(z)
 *
 * @param pba        Input: pointer to background structure
 * @param pfo        Input: pointer to fourier structure
 * @param z          Input: redshift
 * @param sigma_8    Output: sigma8(z)
 * @param sigma_8_cb Output: sigma8_cb(z)
 * @param pnw        Output: pointer to nonlinear workspace
 * @return the error status
 */

int fourier_hmcode_sigma8_at_z(
                               struct background *pba,
                               struct fourier * pfo,
                               double z,
                               double * sigma_8,
                               double * sigma_8_cb,
                               struct fourier_workspace * pnw
                               ) {

  double tau;

  class_call(background_tau_of_z(pba,
                                 z,
                                 &tau),
             pba->error_message,
             pfo->error_message);

  if (pfo->tau_size == 1) {
    *sigma_8 = pnw->sigma_8[pfo->index_pk_m][0];
  }
  else {
    class_call(array_interpolate_two(pfo->tau,
                                     1,
                                     0,
                                     pnw->sigma_8[pfo->index_pk_m],
                                     1,
                                     pfo->tau_size,
                                     tau,
                                     sigma_8,
                                     1,
                                     pfo->error_message),
               pfo->error_message,
               pfo->error_message);
  }


  if (pba->has_ncdm == _TRUE_){

    if (pfo->tau_size == 1) {
      *sigma_8_cb = pnw->sigma_8[pfo->index_pk_cb][0];
    }
    else {
      class_call(array_interpolate_two(pfo->tau,
                                       1,
                                       0,
                                       pnw->sigma_8[pfo->index_pk_cb],
                                       1,
                                       pfo->tau_size,
                                       tau,
                                       sigma_8_cb,
                                       1,
                                       pfo->error_message),
                 pfo->error_message,
                 pfo->error_message);
    }

  }
  else{
    *sigma_8_cb = *sigma_8;
  }



  return _SUCCESS_;
}

/**
 * Compute sigmadisp(z)
 *
 * @param pba           Input: pointer to background structure
 * @param pfo           Input: pointer to fourier structure
 * @param z             Input: redshift
 * @param sigma_disp    Output: sigmadisp(z)
 * @param sigma_disp_cb Output: sigmadisp_cb(z)
 * @param pnw           Output: pointer to nonlinear workspace
 * @return the error status
 */

int fourier_hmcode_sigmadisp_at_z(
                                  struct background *pba,
                                  struct fourier * pfo,
                                  double z,
                                  double * sigma_disp,
                                  double * sigma_disp_cb,
                                  struct fourier_workspace * pnw
                                  ) {

  double tau;

  class_call(background_tau_of_z(pba,
                                 z,
                                 &tau),
             pba->error_message,
             pfo->error_message);

  if (pfo->tau_size == 1) {
    *sigma_disp = pnw->sigma_disp[pfo->index_pk_m][0];
  }
  else {
    class_call(array_interpolate_two(pfo->tau,
                                     1,
                                     0,
                                     pnw->sigma_disp[pfo->index_pk_m],
                                     1,
                                     pfo->tau_size,
                                     tau,
                                     sigma_disp,
                                     1,
                                     pfo->error_message),
               pfo->error_message,
               pfo->error_message);
  }

  if (pba->has_ncdm == _TRUE_){

    if (pfo->tau_size == 1) {
      *sigma_disp_cb = pnw->sigma_disp[pfo->index_pk_cb][0];
    }
    else {
      class_call(array_interpolate_two(pfo->tau,
                                       1,
                                       0,
                                       pnw->sigma_disp[pfo->index_pk_cb],
                                       1,
                                       pfo->tau_size,
                                       tau,
                                       sigma_disp_cb,
                                       1,
                                       pfo->error_message),
                 pfo->error_message,
                 pfo->error_message);
    }

  }
  else{
    *sigma_disp_cb = *sigma_disp;
  }



  return _SUCCESS_;
}

/**
 * Compute sigmadisp100(z)
 *
 * @param pba               Input: pointer to background structure
 * @param pfo               Input: pointer to fourier structure
 * @param z                 Input: redshift
 * @param sigma_disp_100    Output: sigmadisp100(z)
 * @param sigma_disp_100_cb Output: sigmadisp100_cb(z)
 * @param pnw           Output: pointer to nonlinear workspace
 * @return the error status
 */

int fourier_hmcode_sigmadisp100_at_z(
                                     struct background *pba,
                                     struct fourier * pfo,
                                     double z,
                                     double * sigma_disp_100,
                                     double * sigma_disp_100_cb,
                                     struct fourier_workspace * pnw
                                     ) {

  double tau;

  class_call(background_tau_of_z(pba,
                                 z,
                                 &tau),
             pba->error_message,
             pfo->error_message);

  if (pfo->tau_size == 1) {
    *sigma_disp_100 = pnw->sigma_disp_100[pfo->index_pk_m][0];
  }
  else {
    class_call(array_interpolate_two(pfo->tau,
                                     1,
                                     0,
                                     pnw->sigma_disp_100[pfo->index_pk_m],
                                     1,
                                     pfo->tau_size,
                                     tau,
                                     sigma_disp_100,
                                     1,
                                     pfo->error_message),
               pfo->error_message,
               pfo->error_message);
  }

  if (pba->has_ncdm == _TRUE_){

    if (pfo->tau_size == 1) {
      *sigma_disp_100_cb = pnw->sigma_disp_100[pfo->index_pk_cb][0];
    }
    else {
      class_call(array_interpolate_two(pfo->tau,
                                       1,
                                       0,
                                       pnw->sigma_disp_100[pfo->index_pk_cb],
                                       1,
                                       pfo->tau_size,
                                       tau,
                                       sigma_disp_100_cb,
                                       1,
                                       pfo->error_message),
                 pfo->error_message,
                 pfo->error_message);
    }

  }
  else{
    *sigma_disp_100_cb = *sigma_disp_100;
  }


  return _SUCCESS_;
}

/**
 * Compute sigma'(z)
 *
 * @param pba            Input: pointer to background structure
 * @param pfo            Input: pointer to fourier structure
 * @param z              Input: redshift
 * @param sigma_prime    Output: sigma'(z)
 * @param sigma_prime_cb Output: sigma'_cb(z)
 * @param pnw            Output: pointer to nonlinear workspace
 * @return the error status
 */

int fourier_hmcode_sigmaprime_at_z(
                                   struct background *pba,
                                   struct fourier * pfo,
                                   double z,
                                   double * sigma_prime,
                                   double * sigma_prime_cb,
                                   struct fourier_workspace * pnw
                                   ) {

  double tau;

  class_call(background_tau_of_z(pba,
                                 z,
                                 &tau),
             pba->error_message,
             pfo->error_message);

  if (pfo->tau_size == 1) {
    *sigma_prime = pnw->sigma_prime[pfo->index_pk_m][0];
  }
  else {
    class_call(array_interpolate_two(pfo->tau,
                                     1,
                                     0,
                                     pnw->sigma_prime[pfo->index_pk_m],
                                     1,
                                     pfo->tau_size,
                                     tau,
                                     sigma_prime,
                                     1,
                                     pfo->error_message),
               pfo->error_message,
               pfo->error_message);
  }

  if (pba->has_ncdm == _TRUE_){

    if (pfo->tau_size == 1) {
      *sigma_prime_cb = pnw->sigma_prime[pfo->index_pk_cb][0];
    }
    else {
      class_call(array_interpolate_two(pfo->tau,
                                       1,
                                       0,
                                       pnw->sigma_prime[pfo->index_pk_cb],
                                       1,
                                       pfo->tau_size,
                                       tau,
                                       sigma_prime_cb,
                                       1,
                                       pfo->error_message),
                 pfo->error_message,
                 pfo->error_message);
    }

  }
  else{
    *sigma_prime_cb = *sigma_prime;
  }


  return _SUCCESS_;
}
