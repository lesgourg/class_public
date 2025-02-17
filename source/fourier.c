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
#include "halofit.h"
#include "hmcode.h"

/**
 * Return the P(k,z) for a given redshift z and pk type (_m, _cb)
 * (linear if pk_output = pk_linear, nonlinear if pk_output = pk_nonlinear,
 * nowiggle linear spectrum if pk_output = pk_numerical_nowiggle,
 * analytic approximation to nowiggle linear spectrum if pk_output = pk_analytic_nowiggle)
 *
 * In the linear case, if there are several initial conditions *and* the
 * input pointer out_pk_ic is not set to NULL, the function also
 * returns the decomposition into different IC contributions.
 *
 * In the pk_analytic_nowiggle case, the overall normalisation of the
 * spectrum is currently arbitrary and independent of redhsift.
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
 * @param pk_output   Input: linear, nonlinear, nowiggle...
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

  class_test(pk_output == pk_nonlinear && pfo->method == nl_none, pfo->error_message, "Cannot get nonlinear power spectrum when no nonlinear method is employed");

  class_test(pk_output == pk_numerical_nowiggle && pfo->has_pk_numerical_nowiggle == _FALSE_, pfo->error_message, "Cannot get nowiggle power spectrum since it was not requested in input");

  class_test(pk_output == pk_analytic_nowiggle && pfo->has_pk_analytic_nowiggle == _FALSE_, pfo->error_message, "Cannot get analytic nowiggle power spectrum since it was not requested in input");

  /** - case z=0 requiring no interpolation in z. The
        pk_analytic_nowiggle can also be computed here because, in the
        current implementation, it is only computed at z=0. */
  if ((z == 0) || (pk_output == pk_analytic_nowiggle)) {

    for (index_k=0; index_k<pfo->k_size; index_k++) {

      switch (pk_output) {

      case pk_linear:
        out_pk[index_k] = pfo->ln_pk_l[index_pk][(pfo->ln_tau_size-1)*pfo->k_size+index_k];
        if (do_ic == _TRUE_) {
          for (index_ic1_ic2 = 0; index_ic1_ic2 < pfo->ic_ic_size; index_ic1_ic2++) {
            out_pk_ic[index_k * pfo->ic_ic_size + index_ic1_ic2] =
              pfo->ln_pk_ic_l[index_pk][((pfo->ln_tau_size-1)*pfo->k_size+index_k)*pfo->ic_ic_size+index_ic1_ic2];
          }
        }
        break;

      case pk_nonlinear:
        out_pk[index_k] = pfo->ln_pk_nl[index_pk][(pfo->ln_tau_size-1)*pfo->k_size+index_k];
        break;

      case pk_numerical_nowiggle:
        out_pk[index_k] = pfo->ln_pk_l_nw_extra[(pfo->ln_tau_size-1)*pfo->k_size_extra+index_k];
        break;

      case pk_analytic_nowiggle:
        out_pk[index_k] = pfo->ln_pk_l_an_extra[index_k];
        break;

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
      class_test(ln_tau<pfo->ln_tau[0]-100.*_EPSILON_,
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
        else if (pk_output == pk_nonlinear) {
          out_pk[index_k] = pfo->ln_pk_nl[index_pk][index_k];
        }
        else {
          out_pk[index_k] = pfo->ln_pk_l_nw_extra[index_k];
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
        else if (pk_output == pk_nonlinear) {
          out_pk[index_k] = pfo->ln_pk_nl[index_pk][(pfo->ln_tau_size-1) * pfo->k_size + index_k];
        }
        else {
          out_pk[index_k] = pfo->ln_pk_l_nw_extra[(pfo->ln_tau_size-1) * pfo->k_size_extra + index_k];
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

      else if (pk_output == pk_nonlinear) {

        if (ln_tau < pfo->ln_tau[pfo->ln_tau_size-pfo->ln_tau_size_nl]) {
          /** --> we requested P_nl(k) at tau where P_l is computed but the NL correction is NOT -> Return P_l */
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

        }
        else {
          /** --> interpolate P_nl(k) at tau from pre-computed array */
          class_call(array_interpolate_spline(pfo->ln_tau+(pfo->ln_tau_size-pfo->ln_tau_size_nl),
                                              pfo->ln_tau_size_nl,
                                              pfo->ln_pk_nl[index_pk]+(pfo->ln_tau_size-pfo->ln_tau_size_nl)*pfo->k_size,
                                              pfo->ddln_pk_nl[index_pk]+(pfo->ln_tau_size-pfo->ln_tau_size_nl)*pfo->k_size,
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
      else {

        /** --> interpolate P_l_nw(k) at tau from pre-computed array */
        class_call(array_interpolate_spline(pfo->ln_tau,
                                            pfo->ln_tau_size,
                                            pfo->ln_pk_l_nw_extra,
                                            pfo->ddln_pk_l_nw_extra,
                                            pfo->k_size_extra,
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
 * (linear if pk_output = pk_linear, nonlinear if pk_output = pk_nonlinear,
 * nowiggle linear spectrum if pk_output = pk_numerical_nowiggle,
 * analytic approximation to linear nowiggle spectrum if pk_output = pk_analytic_nowiggle)
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
 * @param pk_output   Input: linear, nonlinear, nowiggle...
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
 * @param pk_output      Input: pk_linear, pk_nonlinear, nowiggle...
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
 * Return the logarithmic slope of P(k,z) for a given (k,z), a given
 * pk type (_m, _cb) (computed with linear P_L if pk_output =
 * pk_linear, nonlinear P_NL if pk_output = pk_nonlinear,
 * nowiggle linear spectrum if pk_output = pk_numerical_nowiggle,
 * analytic approximation to linear nowiggle spectrum if pk_output = pk_analytic_nowiggle)
 *
 * @param pba         Input: pointer to background structure
 * @param ppm         Input: pointer to primordial structure
 * @param pfo         Input: pointer to fourier structure
 * @param pk_output   Input: linear, nonlinear, nowiggle...
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
  int index_tau_desired_nl;

  struct hmcode_workspace hw;
  struct hmcode_workspace * phw;

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

  class_call(fourier_indices(ppr,
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

      class_call(fourier_pk_linear(pba,
                                   ppt,
                                   ppm,
                                   pfo,
                                   index_pk,
                                   index_tau_sources,
                                   pfo->k_size,
                                   &(pfo->ln_pk_l[index_pk][index_tau * pfo->k_size]),
                                   &(pfo->ln_pk_ic_l[index_pk][index_tau * pfo->k_size * pfo->ic_ic_size])),
                 pfo->error_message,
                 pfo->error_message);

      /** --> one more call to get the linear power spectrum
          extrapolated up to very large k, for this time and type,
          but ignoring the case of multiple initial
          conditions. Result stored in ln_pk_l_extra (different
          from non-extrapolated ln_pk_l) */

      class_call(fourier_pk_linear(
                                   pba,
                                   ppt,
                                   ppm,
                                   pfo,
                                   index_pk,
                                   index_tau_sources,
                                   pfo->k_size_extra,
                                   &(pfo->ln_pk_l_extra[index_pk][index_tau * pfo->k_size_extra]),
                                   NULL
                                   ),
                 pfo->error_message,
                 pfo->error_message);

    }
  }

  /** - if interpolation of \f$P(k,\tau)\f$ will be needed (as a
      function of tau), compute array of second derivatives in view of
      spline interpolation */

  if (pfo->ln_tau_size > 1) {
    for (index_pk=0; index_pk<pfo->pk_size; index_pk++) {

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

      class_call(array_spline_table_lines(pfo->ln_tau,
                                          pfo->ln_tau_size,
                                          pfo->ln_pk_l_extra[index_pk],
                                          pfo->k_size_extra,
                                          pfo->ddln_pk_l_extra[index_pk],
                                          _SPLINE_EST_DERIV_,
                                          pfo->error_message),
                 pfo->error_message,
                 pfo->error_message);
    }
  }

  /** - get the analytic_nowiggle power spectrum */

  if (pfo->has_pk_analytic_nowiggle == _TRUE_) {

    class_call(fourier_pk_analytic_nowiggle(ppr,pba,ppm,pfo),
               pfo->error_message,
               pfo->error_message);
  }

  /** - get the dewiggled power spectrum at each time in ln_tau */
  if (pfo->has_pk_numerical_nowiggle == _TRUE_) {

    if (pfo->fourier_verbose > 2)
      printf("Computing nowiggle power spectra.\n");

    class_call(fourier_wnw_split(ppr,pba,ppm,pfo),
               pfo->error_message,
               pfo->error_message);
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
                pfo->pk_size*sizeof(double*),
                pfo->error_message);

    class_alloc(lnpk_l,
                pfo->pk_size*sizeof(double*),
                pfo->error_message);

    class_alloc(ddlnpk_l,
                pfo->pk_size*sizeof(double*),
                pfo->error_message);

    for (index_pk=0; index_pk<pfo->pk_size; index_pk++){
      class_alloc(pk_nl[index_pk],pfo->k_size*sizeof(double),pfo->error_message);
      class_alloc(lnpk_l[index_pk],pfo->k_size_extra*sizeof(double),pfo->error_message);
      class_alloc(ddlnpk_l[index_pk],pfo->k_size_extra*sizeof(double),pfo->error_message);
    }

    /** --> Then go through preliminary steps specific to HMcode */

    if (pfo->method == nl_HMcode){

      phw = &hw;

      class_call(hmcode_workspace_init(ppr,pba,pfo,phw),
                 pfo->error_message,
                 pfo->error_message);

      class_call(hmcode_noradiation_growth_init(ppr,pba,pfo,phw),
                 pfo->error_message,
                 pfo->error_message);

      class_call(hmcode_dark_energy_correction(ppr,pba,pfo,phw),
                 pfo->error_message,
                 pfo->error_message);

      class_call(hmcode_baryonic_feedback(pfo),
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

    /* this size will refer to the computed Pk_NL range. We will reduce it if the corrections cannot be computed. */
    pfo->ln_tau_size_nl = pfo->ln_tau_size;

    /* if we do not want any Cl’s derived from large scale structure source functions, try to compute non-linear
	     corrections only in the range z<z_max_pk, that is, index_tau >= pfo->tau_size - pfo->ln_tau_size */
    if ((ppt->has_cl_cmb_lensing_potential == _FALSE_) &&
        (ppt->has_cl_lensing_potential == _FALSE_) &&
        (ppt->has_cl_number_count == _FALSE_)) {
		  index_tau_desired_nl = pfo->tau_size - pfo->ln_tau_size;
    }
    /* otherwise, try to compute non-linear
       corrections up the the highest possible redshift, that is, starting from the smallest possible time */
    else {
      index_tau_desired_nl = 0;
    }

    /* loop over time. Go backward, starting from today and going back to earlier times. */
    for (index_tau = pfo->tau_size-1; index_tau>=index_tau_desired_nl; index_tau--) {

      /* loop over index_pk, defined such that it is ensured
       * that index_pk starts at index_pk_cb when neutrinos are
       * included. This is necessary for hmcode, since the sigmatable
       * needs to be filled for sigma_cb only. Thus, when HMcode
       * evalutes P_m_nl, it needs both P_m_l and P_cb_l. */

      for (index_pk=0; index_pk<pfo->pk_size; index_pk++) {

        /* if we are still in a range of time where P_NL(k) should be computable */
        if (nl_corr_not_computable_at_this_k == _FALSE_) {

          /* get P_L(k) at this time in order to estimate P_NL(k) */

          /* we call fourier_pk_linear() once more. Note that we do
             not use the results stored in pfo->ln_pk_l_extra, because
             here we investigate different values of tau than when
             storing pfo->ln_pk_l, pfo->ln_pk_l_extra, ... */
          class_call(fourier_pk_linear(pba,
                                       ppt,
                                       ppm,
                                       pfo,
                                       index_pk,
                                       index_tau,
                                       pfo->k_size_extra,
                                       lnpk_l[index_pk],
                                       NULL),
                     pfo->error_message,
                     pfo->error_message);

          /* spline P_L(k) at this time along k */
          class_call(array_spline_table_columns(pfo->ln_k,
                                                pfo->k_size_extra,
                                                lnpk_l[index_pk],
                                                1,
                                                ddlnpk_l[index_pk],
                                                _SPLINE_NATURAL_,
                                                pfo->error_message),
                     pfo->error_message,
                     pfo->error_message);

          /* get P_NL(k) at this time with Halofit */
          if (pfo->method == nl_halofit) {

            class_call(halofit(ppr,
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

            /* (preliminary step: fill table of sigma(R) */
            class_call(hmcode_fill_sigtab(ppr,
                                          pba,
                                          ppt,
                                          ppm,
                                          pfo,
                                          index_tau,
                                          index_pk,
                                          lnpk_l,
                                          ddlnpk_l,
                                          phw),
                       pfo->error_message, pfo->error_message);

            class_call(hmcode(ppr,
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
                              phw),
                       pfo->error_message,
                       pfo->error_message);
          }
          else {
            class_stop(pfo->error_message,"nonlinear method not recognized.");
          }

          /* Above, we have checked the computability of NL corrections.
             If they could be computed, infer and store R_NL=(P_NL/P_L)^1/2, and also store P_NL */
           if (nl_corr_not_computable_at_this_k == _FALSE_) {

            /* Store R_NL */
             for (index_k=0; index_k<pfo->k_size; index_k++) {
               pfo->nl_corr_density[index_pk][index_tau * pfo->k_size + index_k] = sqrt(pk_nl[index_pk][index_k]/exp(lnpk_l[index_pk][index_k]));
             }

            /* Store P_NL (only if the output is requested, i.e. z < z_max_pk),
               which is equivalent to index_tau >= pfo->tau_size - pfo->ln_tau_size */
            if (index_tau >= pfo->tau_size - pfo->ln_tau_size) {

              index_tau_late = index_tau - (pfo->tau_size - pfo->ln_tau_size);

              for (index_k=0; index_k<pfo->k_size; index_k++) {
                pfo->ln_pk_nl[index_pk][index_tau_late * pfo->k_size + index_k] = log(pk_nl[index_pk][index_k]);
              }
            }
          }

          /* otherwise we met the first problematic value of time */
          else {

            /* store the index of that value */
            pfo->index_tau_min_nl = MIN(pfo->tau_size-1,index_tau+1); //this MIN() ensures that index_tau_min_nl is never out of bounds
            pfo->ln_tau_size_nl = MIN(pfo->tau_size-1-index_tau,pfo->ln_tau_size); //this will be at most ln_tau_size, and at least 0 (if NL corr cannot be computed at first step)

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

            class_test(pfo->ln_tau_size_nl >1 && pfo->ln_tau_size_nl< 3, pfo->error_message, "Not enough redshifts had a non-linear correction computed, so spline-interpolation of the non-linear correction is impossible. Increase P_k_max_h/Mpc or P_k_max_1/Mpc or fourier_min_k_max.");
          }
        }

        /* if we are still in a range of time where P_NL(k) should NOT be computable */
        else {

          /* store R_NL=1 for that time (very fast) */
          for (index_k=0; index_k<pfo->k_size; index_k++) {
            pfo->nl_corr_density[index_pk][index_tau * pfo->k_size + index_k] = 1.;
          }

        }

      } // end loop over index_pk
    } //end loop over index_tau


    /** --> spline the array of nonlinear power spectrum
            (only above the first index where it could be computed) */
    if (pfo->ln_tau_size_nl > 1) {
      for (index_pk=0; index_pk<pfo->pk_size; index_pk++) {

        class_call(array_spline_table_lines(pfo->ln_tau+(pfo->ln_tau_size-pfo->ln_tau_size_nl),
                                            pfo->ln_tau_size_nl,
                                            pfo->ln_pk_nl[index_pk]+(pfo->ln_tau_size-pfo->ln_tau_size_nl)*pfo->k_size,
                                            pfo->k_size,
                                            pfo->ddln_pk_nl[index_pk]+(pfo->ln_tau_size-pfo->ln_tau_size_nl)*pfo->k_size,
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

      class_call(hmcode_workspace_free(pfo,phw),
                 pfo->error_message,
                 pfo->error_message);

      class_call(hmcode_noradiation_growth_free(pfo,phw),
                 pfo->error_message,
                 pfo->error_message);
    }
  }

  /** - if the nl_method could not be identified */
  else {
    class_stop(pfo->error_message,
               "Your non-linear method variable is set to %d, out of the range defined in fourier.h",pfo->method);
  }

  pfo->is_allocated = _TRUE_;
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
      free(pfo->ln_pk_l_extra[index_pk]);
      if (pfo->ln_tau_size>1) {
        free(pfo->ddln_pk_ic_l[index_pk]);
        free(pfo->ddln_pk_l[index_pk]);
        free(pfo->ddln_pk_l_extra[index_pk]);
      }
    }
    free(pfo->ln_pk_ic_l);
    free(pfo->ln_pk_l);
    free(pfo->ln_pk_l_extra);

    if (pfo->has_pk_analytic_nowiggle == _TRUE_) {
      free(pfo->ln_pk_l_an_extra);
      free(pfo->ddln_pk_l_an_extra);
    }

    free (pfo->sigma8);

    if (pfo->ln_tau_size>1) {
      free(pfo->ddln_pk_ic_l);
      free(pfo->ddln_pk_l);
      free(pfo->ddln_pk_l_extra);
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

  if (pfo->has_pk_numerical_nowiggle) {
    free(pfo->ln_pk_l_nw_extra);
    if (pfo->ln_tau_size > 1)
      free(pfo->ddln_pk_l_nw_extra);
  }

  pfo->is_allocated = _FALSE_;
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

  class_call(fourier_get_k_list(ppr,ppm,ppt,pfo),
             pfo->error_message,
             pfo->error_message);

  /** - get list of tau values */

  class_call(fourier_get_tau_list(ppt,pfo),
             pfo->error_message,
             pfo->error_message);

  /** - given previous indices, we can allocate the array of linear power spectrum values */

  class_alloc(pfo->ln_pk_ic_l,pfo->pk_size*sizeof(double*),pfo->error_message);
  class_alloc(pfo->ln_pk_l   ,pfo->pk_size*sizeof(double*),pfo->error_message);
  class_alloc(pfo->ln_pk_l_extra,pfo->pk_size*sizeof(double*),pfo->error_message);

  for (index_pk=0; index_pk<pfo->pk_size; index_pk++) {
    class_alloc(pfo->ln_pk_ic_l[index_pk],pfo->ln_tau_size*pfo->k_size*pfo->ic_ic_size*sizeof(double*),pfo->error_message);
    class_alloc(pfo->ln_pk_l[index_pk]   ,pfo->ln_tau_size*pfo->k_size*sizeof(double*),pfo->error_message);
    class_alloc(pfo->ln_pk_l_extra[index_pk],pfo->ln_tau_size*pfo->k_size_extra*sizeof(double),pfo->error_message);
  }

  /** - do we we want to compute and store the analytic_nowiggle power
        spectrum? This flag was already set by the user in input.c. If
        we need HMcode, overwrite the user request and compute
        it. Will soon do the same for Oneloop. */

  if ((pfo->method == nl_HMcode) || (pfo->has_pk_numerical_nowiggle == _TRUE_)) {
    pfo->has_pk_analytic_nowiggle = _TRUE_;
  }

  if (pfo->has_pk_analytic_nowiggle == _TRUE_) {
    class_alloc(pfo->ln_pk_l_an_extra,pfo->k_size_extra*sizeof(double),pfo->error_message);
    class_alloc(pfo->ddln_pk_l_an_extra,pfo->k_size_extra*sizeof(double),pfo->error_message);
  }

  /** - if interpolation of \f$P(k,\tau)\f$ will be needed (as a function of tau),
      compute also the array of second derivatives in view of spline interpolation */

  if (pfo->ln_tau_size > 1) {

    class_alloc(pfo->ddln_pk_ic_l,pfo->pk_size*sizeof(double*),pfo->error_message);
    class_alloc(pfo->ddln_pk_l   ,pfo->pk_size*sizeof(double*),pfo->error_message);
    class_alloc(pfo->ddln_pk_l_extra,pfo->pk_size*sizeof(double*),pfo->error_message);

    for (index_pk=0; index_pk<pfo->pk_size; index_pk++) {
      class_alloc(pfo->ddln_pk_ic_l[index_pk],pfo->ln_tau_size*pfo->k_size*pfo->ic_ic_size*sizeof(double*),pfo->error_message);
      class_alloc(pfo->ddln_pk_l[index_pk]   ,pfo->ln_tau_size*pfo->k_size*sizeof(double*),pfo->error_message);
      class_alloc(pfo->ddln_pk_l_extra[index_pk],pfo->ln_tau_size*pfo->k_size_extra*sizeof(double),pfo->error_message);
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

  /** - allocate arrays for P_nw(k,z) and its splines in log(tau) */

  if (pfo->has_pk_numerical_nowiggle == _TRUE_) {
    class_alloc(pfo->ln_pk_l_nw_extra, pfo->ln_tau_size*pfo->k_size_extra*sizeof(double), pfo->error_message);

    if (pfo->ln_tau_size > 1)
      class_alloc(pfo->ddln_pk_l_nw_extra, pfo->ln_tau_size*pfo->k_size_extra*sizeof(double), pfo->error_message);
  }

  return _SUCCESS_;
}

/**
 * Copy list of k from perturbation module, and extended it if
 * necessary to larger k for extrapolation (currently this
 * extrapolation is required only by HMcode)
 *
 * @param ppr Input: pointer to precision structure
 * @param ppm Input: pointer to primordial structure
 * @param ppt Input: pointer to perturbation structure
 * @param pfo Input/Output: pointer to fourier structure
 * @return the error status
 */

int fourier_get_k_list(
                       struct precision *ppr,
                       struct primordial * ppm,
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
  if ((pfo->method == nl_HMcode) || (pfo->has_pk_numerical_nowiggle == _TRUE_)){
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

  class_test(pfo->k[pfo->k_size_extra-1]>exp(ppm->lnk[ppm->lnk_size-1]) && ppm->primordial_spec_type != analytic_Pk,
             pfo->error_message,
             "Setting the output to HMcode with a large 'hmcode_max_k_extra' and using the primordial spectrum to not analytic is incompatible. Either use the analytic power spectrum or set a smaller 'hmcode_max_k_extra' (k_max_hmcode=%.5e , k_max_primordial=%.5e)",
             pfo->k[pfo->k_size_extra-1],
             exp(ppm->lnk[ppm->lnk_size-1])
             );

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
 * Compute a smooth analytic approximation to the matter power
 * spectrum today. Store the results in the fourier structure.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param ppm Input: pointer to primordial structure
 * @param pfo Input/Output: pointer to fourier structure
 * @return the error status
 */

int fourier_pk_analytic_nowiggle(
                                 struct precision *ppr,
                                 struct background *pba,
                                 struct primordial * ppm,
                                 struct fourier *pfo
                                 ) {

  /* Eisenstein & Hu, implemented exactly like in HMcode, except
     for contribution of primodial spectrum (amplitude, tilt), now
     more accurate */
  class_call(hmcode_eisenstein_hu(ppr,
                                  pba,
                                  ppm,
                                  pfo,
                                  pfo->ln_k,
                                  pfo->k_size_extra,
                                  pfo->ln_pk_l_an_extra),
             pfo->error_message,
             pfo->error_message);

  class_call(array_spline_table_lines(pfo->ln_k,
                                      pfo->k_size_extra,
                                      pfo->ln_pk_l_an_extra,
                                      1,
                                      pfo->ddln_pk_l_an_extra,
                                      _SPLINE_NATURAL_,
                                      pfo->error_message),
             pfo->error_message,
             pfo->error_message);

  return _SUCCESS_;
}

/**
 * Compute the decomposition of the linear power spectrum into a
 * wiggly and a non-wiggly part. Store the results in the fourier
 * structure.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param ppm Input: pointer to primordial structure
 * @param pfo Input/Output: pointer to fourier structure
 * @return the error status
 */

int fourier_wnw_split(
                      struct precision *ppr,
                      struct background *pba,
                      struct primordial * ppm,
                      struct fourier *pfo
                      ) {

  class_call(hmcode_wnw_split(ppr,pba,ppm,pfo),
             pfo->error_message,
             pfo->error_message);

  if (pfo->ln_tau_size > 1) {
    /** - spline the nowiggle spectrum with respect to time */
    class_call(array_spline_table_lines(pfo->ln_tau,
                                        pfo->ln_tau_size,
                                        pfo->ln_pk_l_nw_extra,
                                        pfo->k_size_extra,
                                        pfo->ddln_pk_l_nw_extra,
                                        _SPLINE_EST_DERIV_,
                                        pfo->error_message),
               pfo->error_message,
               pfo->error_message);
  }

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
