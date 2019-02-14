/** @file spectra.c Documented spectra module
 *
 * Julien Lesgourgues, 25.08.2010
 *
 * This module computes the anisotropy and Fourier power spectra
 * \f$ C_l^{X}, P(k), ... \f$'s given the transfer and Bessel functions
 * (for anisotropy spectra), the source functions (for Fourier spectra)
 * and the primordial spectra.
 *
 * The following functions can be called from other modules:
 *
 * -# spectra_init() at the beginning (but after transfer_init())
 * -# spectra_cl_at_l() at any time for computing \f$ C_l \f$ at any l
 * -# spectra_spectrum_at_z() at any time for computing P(k) at any z
 * -# spectra_spectrum_at_k_and z() at any time for computing P at any k and z
 * -# spectra_free() at the end
 */

#include "spectra.h"



int spectra_bandpower(struct spectra * psp,
                      int l1,
                      int l2,
                      double * TT_II,
                      double * TT_RI,
                      double * TT_RR
                      ) {

  int l;
  int index_md;
  double * cl_tot;
  double ** cl_md;
  double ** cl_md_ic;

  class_alloc(cl_tot,psp->ct_size*sizeof(double),psp->error_message);
  class_alloc(cl_md,psp->md_size*sizeof(double*),psp->error_message);
  class_alloc(cl_md_ic,psp->md_size*sizeof(double*),psp->error_message);
  for (index_md=0;index_md<psp->md_size; index_md++) {
    class_alloc(cl_md[index_md],psp->ct_size*sizeof(double),psp->error_message);
    class_alloc(cl_md_ic[index_md],psp->ct_size*psp->ic_ic_size[index_md]*sizeof(double),psp->error_message);
  }

  *TT_RR=0.;
  *TT_RI=0.;
  *TT_II=0.;

  for (l=l1; l<=l2; l++) {

    class_call(spectra_cl_at_l(psp,
                               (double)l,
                               cl_tot,
                               cl_md,
                               cl_md_ic),
               psp->error_message,
               psp->error_message);

    *TT_RR += (double)(2*l+1)*cl_md_ic[psp->index_md_scalars][index_symmetric_matrix(0,0,psp->ic_size[psp->index_md_scalars])*psp->ct_size+psp->index_ct_tt];
    *TT_RI += (double)(2*l+1)*cl_md_ic[psp->index_md_scalars][index_symmetric_matrix(0,1,psp->ic_size[psp->index_md_scalars])*psp->ct_size+psp->index_ct_tt]*2.;
    *TT_II += (double)(2*l+1)*cl_md_ic[psp->index_md_scalars][index_symmetric_matrix(1,1,psp->ic_size[psp->index_md_scalars])*psp->ct_size+psp->index_ct_tt];

  }

  for (index_md=0;index_md<psp->md_size; index_md++) {
    free(cl_md[index_md]);
    free(cl_md_ic[index_md]);
  }
  free(cl_tot);
  free(cl_md);
  free(cl_md_ic);

  return _SUCCESS_;

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

int spectra_cl_at_l(
                    struct spectra * psp,
                    double l,
                    double * cl_tot,    /* array with argument cl_tot[index_ct] (must be already allocated) */
                    double * * cl_md,   /* array with argument cl_md[index_md][index_ct] (must be already allocated only if several modes) */
                    double * * cl_md_ic /* array with argument cl_md_ic[index_md][index_ic1_ic2*psp->ct_size+index_ct] (must be already allocated for a given mode only if several ic's) */
                    ) {

  /** Summary: */

  /** - define local variables */

  int last_index;
  int index_md;
  int index_ic1,index_ic2,index_ic1_ic2;
  int index_ct;

  /** - (a) treat case in which there is only one mode and one initial condition.
      Then, only cl_tot needs to be filled. */

  if ((psp->md_size == 1) && (psp->ic_size[0] == 1)) {
    index_md = 0;
    if ((int)l <= psp->l[psp->l_size[index_md]-1]) {

      /* interpolate at l */
      class_call(array_interpolate_spline(psp->l,
                                          psp->l_size[index_md],
                                          psp->cl[index_md],
                                          psp->ddcl[index_md],
                                          psp->ct_size,
                                          l,
                                          &last_index,
                                          cl_tot,
                                          psp->ct_size,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

      /* set to zero for the types such that l<l_max */
      for (index_ct=0; index_ct<psp->ct_size; index_ct++)
        if ((int)l > psp->l_max_ct[index_md][index_ct])
          cl_tot[index_ct]=0.;
    }
    else {
      for (index_ct=0; index_ct<psp->ct_size; index_ct++)
        cl_tot[index_ct]=0.;
    }
  }

  /** - (b) treat case in which there is only one mode
      with several initial condition.
      Fill cl_md_ic[index_md=0] and sum it to get cl_tot. */

  if ((psp->md_size == 1) && (psp->ic_size[0] > 1)) {
    index_md = 0;
    for (index_ct=0; index_ct<psp->ct_size; index_ct++)
      cl_tot[index_ct]=0.;
    for (index_ic1 = 0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
      for (index_ic2 = index_ic1; index_ic2 < psp->ic_size[index_md]; index_ic2++) {
        index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md]);
        if (((int)l <= psp->l[psp->l_size[index_md]-1]) &&
            (psp->is_non_zero[index_md][index_ic1_ic2] == _TRUE_)) {

          class_call(array_interpolate_spline(psp->l,
                                              psp->l_size[index_md],
                                              psp->cl[index_md],
                                              psp->ddcl[index_md],
                                              psp->ic_ic_size[index_md]*psp->ct_size,
                                              l,
                                              &last_index,
                                              cl_md_ic[index_md],
                                              psp->ic_ic_size[index_md]*psp->ct_size,
                                              psp->error_message),
                     psp->error_message,
                     psp->error_message);

          for (index_ct=0; index_ct<psp->ct_size; index_ct++)
            if ((int)l > psp->l_max_ct[index_md][index_ct])
              cl_md_ic[index_md][index_ic1_ic2*psp->ct_size+index_ct]=0.;
        }
        else {
          for (index_ct=0; index_ct<psp->ct_size; index_ct++)
            cl_md_ic[index_md][index_ic1_ic2*psp->ct_size+index_ct]=0.;
        }

        /* compute cl_tot by summing over cl_md_ic */
        for (index_ct=0; index_ct<psp->ct_size; index_ct++) {
          if (index_ic1 == index_ic2)
            cl_tot[index_ct]+=cl_md_ic[index_md][index_ic1_ic2*psp->ct_size+index_ct];
          else
            cl_tot[index_ct]+=2.*cl_md_ic[index_md][index_ic1_ic2*psp->ct_size+index_ct];
        }
      }
    }
  }

  /** - (c) loop over modes */

  if (psp->md_size > 1) {

    for (index_ct=0; index_ct<psp->ct_size; index_ct++)
      cl_tot[index_ct]=0.;

    for (index_md = 0; index_md < psp->md_size; index_md++) {

      /** - --> (c.1.) treat case in which the mode under consideration
          has only one initial condition.
          Fill cl_md[index_md]. */

      if (psp->ic_size[index_md] == 1) {
        if ((int)l <= psp->l[psp->l_size[index_md]-1]) {

          class_call(array_interpolate_spline(psp->l,
                                              psp->l_size[index_md],
                                              psp->cl[index_md],
                                              psp->ddcl[index_md],
                                              psp->ct_size,
                                              l,
                                              &last_index,
                                              cl_md[index_md],
                                              psp->ct_size,
                                              psp->error_message),
                     psp->error_message,
                     psp->error_message);

          for (index_ct=0; index_ct<psp->ct_size; index_ct++)
            if ((int)l > psp->l_max_ct[index_md][index_ct])
              cl_md[index_md][index_ct]=0.;
        }
        else {
          for (index_ct=0; index_ct<psp->ct_size; index_ct++)
            cl_md[index_md][index_ct]=0.;
        }
      }

      /** - --> (c.2.) treat case in which the mode under consideration
          has several initial conditions.
          Fill cl_md_ic[index_md] and sum it to get cl_md[index_md] */

      if (psp->ic_size[index_md] > 1) {

        if ((int)l <= psp->l[psp->l_size[index_md]-1]) {

          /* interpolate all ic and ct */
          class_call(array_interpolate_spline(psp->l,
                                              psp->l_size[index_md],
                                              psp->cl[index_md],
                                              psp->ddcl[index_md],
                                              psp->ic_ic_size[index_md]*psp->ct_size,
                                              l,
                                              &last_index,
                                              cl_md_ic[index_md],
                                              psp->ic_ic_size[index_md]*psp->ct_size,
                                              psp->error_message),
                     psp->error_message,
                     psp->error_message);

          /* set to zero some of the components */
          for (index_ic1 = 0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
            for (index_ic2 = index_ic1; index_ic2 < psp->ic_size[index_md]; index_ic2++) {
              index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md]);
              for (index_ct=0; index_ct<psp->ct_size; index_ct++) {

                if (((int)l > psp->l_max_ct[index_md][index_ct]) || (psp->is_non_zero[index_md][index_ic1_ic2] == _FALSE_))
                  cl_md_ic[index_md][index_ic1_ic2*psp->ct_size+index_ct]=0.;
              }
            }
          }
        }
        /* if l was too big, set anyway all components to zero */
        else {
          for (index_ic1 = 0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
            for (index_ic2 = index_ic1; index_ic2 < psp->ic_size[index_md]; index_ic2++) {
              index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md]);
              for (index_ct=0; index_ct<psp->ct_size; index_ct++) {
                cl_md_ic[index_md][index_ic1_ic2*psp->ct_size+index_ct]=0.;
              }
            }
          }
        }

        /* sum up all ic for each mode */

        for (index_ct=0; index_ct<psp->ct_size; index_ct++) {

          cl_md[index_md][index_ct]=0.;

          for (index_ic1 = 0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
            for (index_ic2 = index_ic1; index_ic2 < psp->ic_size[index_md]; index_ic2++) {
              index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md]);

              if (index_ic1 == index_ic2)
                cl_md[index_md][index_ct]+=cl_md_ic[index_md][index_ic1_ic2*psp->ct_size+index_ct];
              else
                cl_md[index_md][index_ct]+=2.*cl_md_ic[index_md][index_ic1_ic2*psp->ct_size+index_ct];
            }
          }
        }
      }

      /** - --> (c.3.) add contribution of cl_md[index_md] to cl_tot */

      for (index_ct=0; index_ct<psp->ct_size; index_ct++)
        cl_tot[index_ct]+=cl_md[index_md][index_ct];
    }
  }

  return _SUCCESS_;

}

/**
 * Matter power spectrum for arbitrary redshift and for all initial conditions.
 *
 * This routine evaluates the matter power spectrum at a given value of z by
 * interpolating in the pre-computed table (if several values of z have been stored)
 * or by directly reading it (if it only contains values at z=0 and we want P(k,z=0))
 *
 *
 * Can be called in two modes: linear or logarithmic.
 *
 * - linear: returns P(k) (units: \f$ Mpc^3\f$)
 *
 * - logarithmic: returns \f$\ln{P(k)}\f$
 *
 * One little subtlety: in case of several correlated initial conditions,
 * the cross-correlation spectrum can be negative. Then, in logarithmic mode,
 * the non-diagonal elements contain the cross-correlation angle \f$ P_{12}/\sqrt{P_{11} P_{22}}\f$
 * (from -1 to 1) instead of \f$\ln{P_{12}}\f$
 *
 * This function can be
 * called from whatever module at whatever time, provided that
 * spectra_init() has been called before, and spectra_free() has not
 * been called yet.
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

int spectra_pk_at_z(
                    struct background * pba,
                    struct spectra * psp,
                    enum linear_or_logarithmic mode,
                    double z,
                    double * output_tot,    /* array with argument output_tot[index_k] (must be already allocated) */
                    double * output_ic,     /* array with argument output_tot[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2] (must be already allocated only if more than one initial condition) */
                    double * output_cb_tot, /* same as output_tot for the baryon+CDM only */
                    double * output_cb_ic   /* same as output_ic  for the baryon+CDM only */
                    ) {

  /* JL 21.09.2017: TODO: now, P(k) total is already calculated and
     stored in spectra_pk(), in the array psp->pk_l. WE should use
     that here too to compute output_tot, inmstead of redoing the sum
     over ICs. */

  /** Summary: */

  /** - define local variables */

  int index_md;
  int last_index;
  int index_k;
  double tau,ln_tau;
  int index_ic1,index_ic2,index_ic1_ic2;

  index_md = psp->index_md_scalars;

  /** - first step: convert z into \f$\ln{\tau}\f$ */

  class_call(background_tau_of_z(pba,z,&tau),
             pba->error_message,
             psp->error_message);

  class_test(tau <= 0.,
             psp->error_message,
             "negative or null value of conformal time: cannot interpolate");

  ln_tau = log(tau);

  double small_deviation = 1e-10;
  class_test(ln_tau<psp->ln_tau[0]-small_deviation,
             psp->error_message,
             "requested z was not inside of tau tabulation range (Requested %.10e, Min %.10e) ",ln_tau,psp->ln_tau[0]-small_deviation);
  if(ln_tau<psp->ln_tau[0]){
    //Case of small deviation caused by rounding
    ln_tau = psp->ln_tau[0];
  }
  class_test(ln_tau>psp->ln_tau[psp->ln_tau_size-1]+small_deviation,
             psp->error_message,
             "requested z was not inside of tau tabulation range (Requested %.10e, Max %.10e) ",ln_tau,psp->ln_tau[psp->ln_tau_size-1]+small_deviation);

  if(ln_tau>psp->ln_tau[psp->ln_tau_size-1]){
    //Case of small deviation caused by rounding
    ln_tau = psp->ln_tau[psp->ln_tau_size-1];
  }
  /** - second step: for both modes (linear or logarithmic), store the spectrum in logarithmic format in the output array(s) */

  /** - --> (a) if only values at tau=tau_today are stored and we want \f$ P(k,z=0)\f$, no need to interpolate */

  if (psp->ln_tau_size == 1) {

    class_test(z != 0.,
               psp->error_message,
               "asked z=%e but only P(k,z=0) has been tabulated",z);

    for (index_k=0; index_k<psp->ln_k_size; index_k++)
      if (psp->ic_size[index_md] == 1) {
      	output_tot[index_k] = psp->ln_pk[index_k];
        if(pba->has_ncdm) output_cb_tot[index_k] = psp->ln_pk_cb[index_k];
      }
      else {
        for (index_ic1_ic2 = 0; index_ic1_ic2 < psp->ic_ic_size[index_md]; index_ic1_ic2++) {
          output_ic[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2] =
            psp->ln_pk[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2];
          if(pba->has_ncdm)
            output_cb_ic[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2] =
              psp->ln_pk_cb[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2];
        }
      }
  }

  /** - --> (b) if several values of tau have been stored, use interpolation routine to get spectra at correct redshift */

  else {

    if (psp->ic_ic_size[index_md] == 1) {

      class_call(array_interpolate_spline(psp->ln_tau,
                                          psp->ln_tau_size,
                                          psp->ln_pk,
                                          psp->ddln_pk,
                                          psp->ln_k_size,
                                          ln_tau,
                                          &last_index,
                                          output_tot,
                                          psp->ln_k_size,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

      if(pba->has_ncdm){
        class_call(array_interpolate_spline(psp->ln_tau,
                                            psp->ln_tau_size,
                                            psp->ln_pk_cb,
                                            psp->ddln_pk_cb,
                                            psp->ln_k_size,
                                            ln_tau,
                                            &last_index,
                                            output_cb_tot,
                                            psp->ln_k_size,
                                            psp->error_message),
                   psp->error_message,
                   psp->error_message);
      }

    }
    else {

      class_call(array_interpolate_spline(psp->ln_tau,
                                          psp->ln_tau_size,
                                          psp->ln_pk,
                                          psp->ddln_pk,
                                          psp->ic_ic_size[index_md]*psp->ln_k_size,
                                          ln_tau,
                                          &last_index,
                                          output_ic,
                                          psp->ic_ic_size[index_md]*psp->ln_k_size,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

      if(pba->has_ncdm){
        class_call(array_interpolate_spline(psp->ln_tau,
                                            psp->ln_tau_size,
                                            psp->ln_pk_cb,
                                            psp->ddln_pk_cb,
                                            psp->ic_ic_size[index_md]*psp->ln_k_size,
                                            ln_tau,
                                            &last_index,
                                            output_cb_ic,
                                            psp->ic_ic_size[index_md]*psp->ln_k_size,
                                            psp->error_message),
                   psp->error_message,
                   psp->error_message);
      }

    }
  }

  /** - third step: if there are several initial conditions, compute the total P(k) and set back all uncorrelated coefficients to exactly zero. Check positivity of total P(k). */

  if (psp->ic_size[index_md] > 1) {
    for (index_k=0; index_k<psp->ln_k_size; index_k++) {
      output_tot[index_k] = 0.;
      if (pba->has_ncdm) output_cb_tot[index_k] = 0.;
      for (index_ic1=0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
        for (index_ic2 = index_ic1; index_ic2 < psp->ic_size[index_md]; index_ic2++) {
          index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md]);
          if (index_ic1 == index_ic2) {
            output_tot[index_k] += exp(output_ic[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2]);
            if(pba->has_ncdm) output_cb_tot[index_k] += exp(output_cb_ic[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2]);
          }
          else {
            if (psp->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {
              output_tot[index_k] +=
                2. * output_ic[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2] *
                sqrt(exp(output_ic[index_k * psp->ic_ic_size[index_md] + index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_md])]) *
                     exp(output_ic[index_k * psp->ic_ic_size[index_md] + index_symmetric_matrix(index_ic2,index_ic2,psp->ic_size[index_md])]));
              if(pba->has_ncdm){
                output_cb_tot[index_k] +=
                  2. * output_cb_ic[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2] *
                  sqrt(exp(output_cb_ic[index_k * psp->ic_ic_size[index_md] + index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_md])]) *
                       exp(output_cb_ic[index_k * psp->ic_ic_size[index_md] + index_symmetric_matrix(index_ic2,index_ic2,psp->ic_size[index_md])]));
              }
            }
            else{
              output_ic[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2] = 0.;
              if(pba->has_ncdm) output_cb_ic[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2] = 0.;
            }
          }
        }
      }

      class_test(output_tot[index_k] <= 0.,
                 psp->error_message,
                 "for k=%e, z=%e, the matrix of initial condition amplitudes was not positive definite, hence P(k)_total=%e results negative",
                 exp(psp->ln_k[index_k]),z,output_tot[index_k]);

      if(pba->has_ncdm){
        class_test(output_cb_tot[index_k] <= 0.,
                   psp->error_message,
                   "for k=%e, z=%e, the matrix of initial condition amplitudes was not positive definite, hence P(k)_cb_total=%e results negative",
                   exp(psp->ln_k[index_k]),z,output_cb_tot[index_k]);
      }

    }
  }

  /** - fourth step: depending on requested mode (linear or logarithmic), apply necessary transformation to the output arrays */

  /** - --> (a) linear mode: if only one initial condition, convert output_tot to linear format; if several initial conditions, convert output_ic to linear format, output_tot is already in this format */

  if (mode == linear) {

    if (psp->ic_size[index_md] == 1) {
      for (index_k=0; index_k<psp->ln_k_size; index_k++) {
        output_tot[index_k] = exp(output_tot[index_k]);
        if(pba->has_ncdm) output_cb_tot[index_k] = exp(output_cb_tot[index_k]);
      }
    }

    else {
      for (index_k=0; index_k<psp->ln_k_size; index_k++) {
        for (index_ic1=0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
          index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_md]);
          output_ic[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2] = exp(output_ic[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2]);
          if(pba->has_ncdm)
            output_cb_ic[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2] = exp(output_cb_ic[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2]);
        }
        for (index_ic1=0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
          for (index_ic2 = index_ic1+1; index_ic2 < psp->ic_size[index_md]; index_ic2++) {

            output_ic[index_k * psp->ic_ic_size[index_md] + index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md])] =
              output_ic[index_k * psp->ic_ic_size[index_md] + index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md])]
              *sqrt(output_ic[index_k * psp->ic_ic_size[index_md] + index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_md])] *
                    output_ic[index_k * psp->ic_ic_size[index_md] + index_symmetric_matrix(index_ic2,index_ic2,psp->ic_size[index_md])]);

            if(pba->has_ncdm){
              output_cb_ic[index_k * psp->ic_ic_size[index_md] + index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md])] =
                output_cb_ic[index_k * psp->ic_ic_size[index_md] + index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md])]
                *sqrt(output_cb_ic[index_k * psp->ic_ic_size[index_md] + index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_md])] *
                      output_cb_ic[index_k * psp->ic_ic_size[index_md] + index_symmetric_matrix(index_ic2,index_ic2,psp->ic_size[index_md])]);
            }

          }
        }
      }
    }
  }

  /** - --> (b) logarithmic mode: if only one initial condition, nothing to be done; if several initial conditions, convert output_tot to logarithmic format, output_ic is already in this format */

  else {

    if (psp->ic_size[index_md] > 1) {
      for (index_k=0; index_k<psp->ln_k_size; index_k++) {
        /* we have already checked above that output_tot was positive */
        output_tot[index_k] = log(output_tot[index_k]);
        if(pba->has_ncdm) output_cb_tot[index_k] = log(output_cb_tot[index_k]);
      }
    }
  }

  return _SUCCESS_;

}

/**
 * Matter power spectrum for arbitrary wavenumber, redshift and initial condition.
 *
 * This routine evaluates the matter power spectrum at a given value of k and z by
 * interpolating in a table of all P(k)'s computed at this z by spectra_pk_at_z() (when kmin <= k <= kmax),
 * or eventually by using directly the primordial spectrum (when 0 <= k < kmin):
 * the latter case is an approximation, valid when kmin << comoving Hubble scale today.
 * Returns zero when k=0. Returns an error when k<0 or k > kmax.
 *
 * This function can be
 * called from whatever module at whatever time, provided that
 * spectra_init() has been called before, and spectra_free() has not
 * been called yet.
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

int spectra_pk_at_k_and_z(
                          struct background * pba,
                          struct primordial * ppm,
                          struct spectra * psp,
                          double k,
                          double z,
                          double * pk_tot,    /* pointer to a single number (must be already allocated) */
                          double * pk_ic,     /* array of argument pk_ic[index_ic1_ic2] (must be already allocated only if several initial conditions) */
                          double * pk_cb_tot, /* same as pk_tot for baryon+CDM part only */
                          double * pk_cb_ic   /* same as pk_ic  for baryon+CDM part only */
                          ) {

  /** Summary: */

  /** - define local variables */

  int index_md;
  int index_k;
  int last_index;
  int index_ic1,index_ic2,index_ic1_ic2;

  double * spectrum_at_z = NULL;
  double * spectrum_at_z_ic = NULL;
  double * spline;
  double * pk_primordial_k = NULL;
  double kmin;
  double * pk_primordial_kmin = NULL;

  double * spectrum_cb_at_z = NULL;
  double * spectrum_cb_at_z_ic = NULL;
  double * spline_cb = NULL;

  index_md = psp->index_md_scalars;

  /** - first step: check that k is in valid range [0:kmax] (the test for z will be done when calling spectra_pk_at_z()) */

  class_test((k < 0.) || (k > exp(psp->ln_k[psp->ln_k_size-1])),
             psp->error_message,
             "k=%e out of bounds [%e:%e]",k,0.,exp(psp->ln_k[psp->ln_k_size-1]));

  /** - deal with case 0 <= k < kmin */

  if (k < exp(psp->ln_k[0])) {

    /** - --> (a) subcase k=0: then P(k)=0 */

    if (k == 0.) {
      if (psp->ic_size[index_md] == 1) {
        *pk_tot=0.;
        if (pba->has_ncdm) *pk_cb_tot=0.;
      }
      else {
        for (index_ic1_ic2 = 0; index_ic1_ic2 < psp->ic_ic_size[index_md]; index_ic1_ic2++) {
          pk_ic[index_ic1_ic2] = 0.;
          if (pba->has_ncdm) pk_cb_ic[index_ic1_ic2] = 0.;
        }
      }
    }

    /** - --> (b) subcase 0<k<kmin: in this case we know that on super-Hubble scales:
     *          P(k) = [some number] * k  * P_primordial(k)
     *          so
     *          P(k) = P(kmin) * (k P_primordial(k)) / (kmin P_primordial(kmin))
     *          (note that the result is accurate only if kmin is such that [a0 kmin] << H0)
     */

    else {

      /* compute P(k,z) which contains P(kmin,z)*/
      class_alloc(spectrum_at_z,
                  psp->ln_k_size*sizeof(double),
                  psp->error_message);
      class_alloc(spectrum_cb_at_z,
                  psp->ln_k_size*sizeof(double),
                  psp->error_message);
      if (psp->ic_size[index_md] > 1) {
        class_alloc(spectrum_at_z_ic,
                    sizeof(double)*psp->ic_ic_size[index_md]*psp->ln_k_size,
                    psp->error_message);

        class_alloc(spectrum_cb_at_z_ic,
                    sizeof(double)*psp->ic_ic_size[index_md]*psp->ln_k_size,
                    psp->error_message);

      }
      class_call(spectra_pk_at_z(pba,
                                 psp,
                                 linear,
                                 z,
                                 spectrum_at_z,
                                 spectrum_at_z_ic,
                                 spectrum_cb_at_z,
                                 spectrum_cb_at_z_ic),
                 psp->error_message,
                 psp->error_message);

      /* compute P_primordial(k) */
      class_alloc(pk_primordial_k,
                  sizeof(double)*psp->ic_ic_size[index_md],
                  psp->error_message);
      class_call(primordial_spectrum_at_k(ppm,
                                          index_md,
                                          linear,
                                          k,
                                          pk_primordial_k),
                 ppm->error_message,psp->error_message);

      /* compute P_primordial(kmin) */
      kmin = exp(psp->ln_k[0]);
      class_alloc(pk_primordial_kmin,
                  sizeof(double)*psp->ic_ic_size[index_md],
                  psp->error_message);
      class_call(primordial_spectrum_at_k(ppm,
                                          index_md,
                                          linear,
                                          kmin,
                                          pk_primordial_kmin),
                 ppm->error_message,
                 psp->error_message);

      /* apply above analytic approximation for P(k) */
      index_k=0;
      if (psp->ic_size[index_md] == 1) {
        index_ic1_ic2 = 0;
        *pk_tot = spectrum_at_z[index_k]
          *k*pk_primordial_k[index_ic1_ic2]
          /kmin/pk_primordial_kmin[index_ic1_ic2];
        if (pba->has_ncdm){
          *pk_cb_tot = spectrum_cb_at_z[index_k]
            *k*pk_primordial_k[index_ic1_ic2]
            /kmin/pk_primordial_kmin[index_ic1_ic2];
        }
      }
      else {
      	for (index_ic1_ic2 = 0; index_ic1_ic2 < psp->ic_ic_size[index_md]; index_ic1_ic2++) {
          pk_ic[index_ic1_ic2] = spectrum_at_z_ic[index_ic1_ic2]
            *k*pk_primordial_k[index_ic1_ic2]
            /kmin/pk_primordial_kmin[index_ic1_ic2];
          if (pba->has_ncdm){
            pk_cb_ic[index_ic1_ic2] = spectrum_cb_at_z_ic[index_ic1_ic2]
              *k*pk_primordial_k[index_ic1_ic2]
              /kmin/pk_primordial_kmin[index_ic1_ic2];
          }
        }
      }

      free(spectrum_at_z);
      free(spectrum_cb_at_z);
      if (psp->ic_size[index_md] > 1){
        free(spectrum_at_z_ic);
        free(spectrum_cb_at_z_ic);
      }
      free(pk_primordial_k);
      free(pk_primordial_kmin);

    }
  }

  /** - deal with case kmin <= k <= kmax */

  else {

    /* compute P(k,z) (in logarithmic format for more accurate interpolation) */
    class_alloc(spectrum_at_z,
                psp->ln_k_size*sizeof(double),
                psp->error_message);
    class_alloc(spectrum_cb_at_z,
                psp->ln_k_size*sizeof(double),
                psp->error_message);
    if (psp->ic_size[index_md] > 1) {
      class_alloc(spectrum_at_z_ic,
                  sizeof(double)*psp->ic_ic_size[index_md]*psp->ln_k_size,
                  psp->error_message);
      class_alloc(spectrum_cb_at_z_ic,
                  sizeof(double)*psp->ic_ic_size[index_md]*psp->ln_k_size,
                  psp->error_message);
    }
    class_call(spectra_pk_at_z(pba,
                               psp,
                               logarithmic,
                               z,
                               spectrum_at_z,
                               spectrum_at_z_ic,
                               spectrum_cb_at_z,
                               spectrum_cb_at_z_ic),
               psp->error_message,
               psp->error_message);

    /* get its second derivatives with spline, then interpolate, then convert to linear format */

    class_alloc(spline,
                sizeof(double)*psp->ic_ic_size[index_md]*psp->ln_k_size,
                psp->error_message);

    if(pba->has_ncdm)
      class_alloc(spline_cb,
                  sizeof(double)*psp->ic_ic_size[index_md]*psp->ln_k_size,
                  psp->error_message);

    if (psp->ic_size[index_md] == 1) {

      class_call(array_spline_table_lines(psp->ln_k,
                                          psp->ln_k_size,
                                          spectrum_at_z,
                                          1,
                                          spline,
                                          _SPLINE_NATURAL_,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

      class_call(array_interpolate_spline(psp->ln_k,
                                          psp->ln_k_size,
                                          spectrum_at_z,
                                          spline,
                                          1,
                                          log(k),
                                          &last_index,
                                          pk_tot,
                                          1,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

      *pk_tot = exp(*pk_tot);

      if(pba->has_ncdm){
        class_call(array_spline_table_lines(psp->ln_k,
                                            psp->ln_k_size,
                                            spectrum_cb_at_z,
                                            1,
                                            spline_cb,
                                            _SPLINE_NATURAL_,
                                            psp->error_message),
                   psp->error_message,
                   psp->error_message);

        class_call(array_interpolate_spline(psp->ln_k,
                                            psp->ln_k_size,
                                            spectrum_cb_at_z,
                                            spline_cb,
                                            1,
                                            log(k),
                                            &last_index,
                                            pk_cb_tot,
                                            1,
                                            psp->error_message),
                   psp->error_message,
                   psp->error_message);

        *pk_cb_tot = exp(*pk_cb_tot);
      }

    }
    else {

      class_call(array_spline_table_lines(psp->ln_k,
                                          psp->ln_k_size,
                                          spectrum_at_z_ic,
                                          psp->ic_ic_size[index_md],
                                          spline,
                                          _SPLINE_NATURAL_,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

      class_call(array_interpolate_spline(psp->ln_k,
                                          psp->ln_k_size,
                                          spectrum_at_z_ic,
                                          spline,
                                          psp->ic_ic_size[index_md],
                                          log(k),
                                          &last_index,
                                          pk_ic,
                                          psp->ic_ic_size[index_md],
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

      if(pba->has_ncdm){
        class_call(array_spline_table_lines(psp->ln_k,
                                            psp->ln_k_size,
                                            spectrum_cb_at_z_ic,
                                            psp->ic_ic_size[index_md],
                                            spline_cb,
                                            _SPLINE_NATURAL_,
                                            psp->error_message),
                   psp->error_message,
                   psp->error_message);

        class_call(array_interpolate_spline(psp->ln_k,
                                            psp->ln_k_size,
                                            spectrum_cb_at_z_ic,
                                            spline_cb,
                                            psp->ic_ic_size[index_md],
                                            log(k),
                                            &last_index,
                                            pk_cb_ic,
                                            psp->ic_ic_size[index_md],
                                            psp->error_message),
                   psp->error_message,
                   psp->error_message);
      }

      for (index_ic1 = 0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
        index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_md]);
        pk_ic[index_ic1_ic2] = exp(pk_ic[index_ic1_ic2]);
        if(pba->has_ncdm) pk_cb_ic[index_ic1_ic2] = exp(pk_cb_ic[index_ic1_ic2]);
      }
      for (index_ic1 = 0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
        for (index_ic2 = index_ic1+1; index_ic2 < psp->ic_size[index_md]; index_ic2++) {
          index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md]);
          if (psp->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {
            pk_ic[index_ic1_ic2] = pk_ic[index_ic1_ic2]*
              sqrt(pk_ic[index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_md])]*
                   pk_ic[index_symmetric_matrix(index_ic2,index_ic2,psp->ic_size[index_md])]);
            if(pba->has_ncdm){
              pk_cb_ic[index_ic1_ic2] = pk_cb_ic[index_ic1_ic2]*
                sqrt(pk_cb_ic[index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_md])]*
                     pk_cb_ic[index_symmetric_matrix(index_ic2,index_ic2,psp->ic_size[index_md])]);
            }
          }
          else {
            pk_ic[index_ic1_ic2] = 0.;
            if (pba->has_ncdm) pk_cb_ic[index_ic1_ic2] = 0.;
          }
        }
      }
      free(spectrum_at_z_ic);
      free(spectrum_cb_at_z_ic);
    }

    free(spectrum_at_z);
    free(spectrum_cb_at_z);
    free(spline);
    if(pba->has_ncdm) free(spline_cb);
  }

  /** - last step: if more than one condition, sum over pk_ic to get pk_tot, and set back coefficients of non-correlated pairs to exactly zero. */

  if (psp->ic_size[index_md] > 1) {

    *pk_tot = 0.;

    if (pba->has_ncdm) *pk_cb_tot = 0.;

    for (index_ic1 = 0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
      for (index_ic2 = index_ic1; index_ic2 < psp->ic_size[index_md]; index_ic2++) {
        index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md]);

        if (psp->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {

          if (index_ic1 == index_ic2){
            *pk_tot += pk_ic[index_ic1_ic2];
            if(pba->has_ncdm) *pk_cb_tot += pk_cb_ic[index_ic1_ic2];
          }
          else{
            *pk_tot += 2.*pk_ic[index_ic1_ic2];
            if(pba->has_ncdm) *pk_cb_tot += 2.*pk_cb_ic[index_ic1_ic2];
          }
        }
        else {
          pk_ic[index_ic1_ic2] = 0.;
          if(pba->has_ncdm) pk_cb_ic[index_ic1_ic2] = 0.;
        }
      }
    }

    class_test(*pk_tot <= 0.,
               psp->error_message,
               "for k=%e, the matrix of initial condition amplitudes was not positive definite, hence P(k)_total results negative",k);
    if(pba->has_ncdm){
      class_test(*pk_cb_tot <= 0.,
                 psp->error_message,
                 "for k=%e, the matrix of initial condition amplitudes was not positive definite, hence P(k)_cb_total results negative",k);
    }
  }

  return _SUCCESS_;

}

/**
 * Non-linear total matter power spectrum for arbitrary redshift.
 *
 * This routine evaluates the non-linear matter power spectrum at a given value of z by
 * interpolating in the pre-computed table (if several values of z have been stored)
 * or by directly reading it (if it only contains values at z=0 and we want P(k,z=0))
 *
 *
 * Can be called in two modes: linear or logarithmic.
 *
 * - linear: returns P(k) (units: Mpc^3)
 *
 * - logarithmic: returns ln(P(k))
 *
 * This function can be
 * called from whatever module at whatever time, provided that
 * spectra_init() has been called before, and spectra_free() has not
 * been called yet.
 *
 * @param pba           Input: pointer to background structure (used for converting z into tau)
 * @param psp           Input: pointer to spectra structure (containing pre-computed table)
 * @param mode          Input: linear or logarithmic
 * @param z             Input: redshift
 * @param output_tot    Output: total matter power spectrum P(k) in \f$ Mpc^3\f$ (linear mode), or its logarithms (logarithmic mode)
 * @param output_cb_tot Output: b+CDM power spectrum P(k) in \f$ Mpc^3\f$ (linear mode), or its logarithms (logarithmic mode)
 * @return the error status
 */

int spectra_pk_nl_at_z(
                       struct background * pba,
                       struct spectra * psp,
                       enum linear_or_logarithmic mode,
                       double z,
                       double * output_tot,   /* array with argument output_tot[index_k] (must be already allocated) */
                       double * output_cb_tot /* same as output_tot for baryon+CDM only */
                       ) {

  /** Summary: */

  /** - define local variables */

  int last_index;
  int index_k;
  double tau,ln_tau;

  /** - first step: convert z into ln(tau) */

  class_call(background_tau_of_z(pba,z,&tau),
             pba->error_message,
             psp->error_message);

  class_test(tau <= 0.,
             psp->error_message,
             "negative or null value of conformal time: cannot interpolate");

  ln_tau = log(tau);

  /** - second step: for both modes (linear or logarithmic), store the spectrum in logarithmic format in the output array(s) */

  /** - --> (a) if only values at tau=tau_today are stored and we want P(k,z=0), no need to interpolate */

  if (psp->ln_tau_size == 1) {

    class_test(z != 0.,
               psp->error_message,
               "asked z=%e but only P(k,z=0) has been tabulated",z);

    for (index_k=0; index_k<psp->ln_k_size; index_k++) {
      output_tot[index_k] = psp->ln_pk_nl[index_k];
      if (pba->has_ncdm) output_cb_tot[index_k] = psp->ln_pk_cb_nl[index_k];
    }
  }

  /** - --> (b) if several values of tau have been stored, use interpolation routine to get spectra at correct redshift */

  else {

    class_test(ln_tau < psp->ln_tau[0],
               "This should never happen",
               psp->error_message,
               psp->error_message);

    if (ln_tau < psp->ln_tau_nl[0]) {

      class_call(array_interpolate_spline(psp->ln_tau,
                                          psp->ln_tau_size,
                                          psp->ln_pk_l,
                                          psp->ddln_pk_l,
                                          psp->ln_k_size,
                                          ln_tau,
                                          &last_index,
                                          output_tot,
                                          psp->ln_k_size,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);
    }
    else {

      class_call(array_interpolate_spline(psp->ln_tau_nl,
                                          psp->ln_tau_nl_size,
                                          psp->ln_pk_nl,
                                          psp->ddln_pk_nl,
                                          psp->ln_k_size,
                                          ln_tau,
                                          &last_index,
                                          output_tot,
                                          psp->ln_k_size,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);
    }
    if(pba->has_ncdm){
      if(ln_tau < psp->ln_tau_nl[0]){
        class_call(array_interpolate_spline(psp->ln_tau,
                                            psp->ln_tau_size,
                                            psp->ln_pk_cb_l,
                                            psp->ddln_pk_cb_l,
                                            psp->ln_k_size,
                                            ln_tau,
                                            &last_index,
                                            output_cb_tot,
                                            psp->ln_k_size,
                                            psp->error_message),
                   psp->error_message,
                   psp->error_message);

      }
      else{
        class_call(array_interpolate_spline(psp->ln_tau_nl,
                                            psp->ln_tau_nl_size,
                                            psp->ln_pk_cb_nl,
                                            psp->ddln_pk_cb_nl,
                                            psp->ln_k_size,
                                            ln_tau,
                                            &last_index,
                                            output_cb_tot,
                                            psp->ln_k_size,
                                            psp->error_message),
                   psp->error_message,
                   psp->error_message);

      }
    }

  }

  /** - fourth step: eventually convert to linear format */

  if (mode == linear) {
    for (index_k=0; index_k<psp->ln_k_size; index_k++) {
      output_tot[index_k] = exp(output_tot[index_k]);
      if(pba->has_ncdm) output_cb_tot[index_k] = exp(output_cb_tot[index_k]);
    }
  }

  return _SUCCESS_;

}

/**
 * Non-linear total matter power spectrum for arbitrary wavenumber and redshift.
 *
 * This routine evaluates the matter power spectrum at a given value of k and z by
 * interpolating in a table of all P(k)'s computed at this z by spectra_pk_nl_at_z() (when kmin <= k <= kmax),
 * or eventually by using directly the primordial spectrum (when 0 <= k < kmin):
 * the latter case is an approximation, valid when kmin << comoving Hubble scale today.
 * Returns zero when k=0. Returns an error when k<0 or k > kmax.
 *
 * This function can be
 * called from whatever module at whatever time, provided that
 * spectra_init() has been called before, and spectra_free() has not
 * been called yet.
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

int spectra_pk_nl_at_k_and_z(
                             struct background * pba,
                             struct primordial * ppm,
                             struct spectra * psp,
                             double k,
                             double z,
                             double * pk_tot,   /* pointer to a single number (must be already allocated) */
                             double * pk_cb_tot /* same as pk_tot for baryon+CDM only */
                             ) {

  /** Summary: */

  /** - define local variables */

  int index_md;
  int last_index;

  double * spectrum_at_z = NULL;
  double * spline;

  double * spectrum_cb_at_z = NULL;
  double * spline_cb = NULL;

  index_md = psp->index_md_scalars;

  /** - check that k is in valid range [0:kmax] (the test for z will be done when calling spectra_pk_at_z()) */

  class_test((k < exp(psp->ln_k[0])) || (k > exp(psp->ln_k[psp->ln_k_size-1])),
             psp->error_message,
             "k=%e out of bounds [%e:%e]",k,0.,exp(psp->ln_k[psp->ln_k_size-1]));

  /** - compute P(k,z) (in logarithmic format for more accurate interpolation) */
  class_alloc(spectrum_at_z,
              psp->ln_k_size*sizeof(double),
              psp->error_message);

  class_alloc(spectrum_cb_at_z,
              psp->ln_k_size*sizeof(double),
              psp->error_message);

  class_call(spectra_pk_nl_at_z(pba,
                                psp,
                                logarithmic,
                                z,
                                spectrum_at_z,
                                spectrum_cb_at_z),
             psp->error_message,
             psp->error_message);

  /** - get its second derivatives with spline, then interpolate, then convert to linear format */

  class_alloc(spline,
              sizeof(double)*psp->ic_ic_size[index_md]*psp->ln_k_size,
              psp->error_message);

  if(pba->has_ncdm){
    class_alloc(spline_cb,
                sizeof(double)*psp->ic_ic_size[index_md]*psp->ln_k_size,
                psp->error_message);
  }

  class_call(array_spline_table_lines(psp->ln_k,
                                      psp->ln_k_size,
                                      spectrum_at_z,
                                      1,
                                      spline,
                                      _SPLINE_NATURAL_,
                                      psp->error_message),
             psp->error_message,
             psp->error_message);

  class_call(array_interpolate_spline(psp->ln_k,
                                      psp->ln_k_size,
                                      spectrum_at_z,
                                      spline,
                                      1,
                                      log(k),
                                      &last_index,
                                      pk_tot,
                                      1,
                                      psp->error_message),
             psp->error_message,
             psp->error_message);

  *pk_tot = exp(*pk_tot);

  if(pba->has_ncdm){

    class_call(array_spline_table_lines(psp->ln_k,
                                        psp->ln_k_size,
                                        spectrum_cb_at_z,
                                        1,
                                        spline_cb,
                                        _SPLINE_NATURAL_,
                                        psp->error_message),
               psp->error_message,
               psp->error_message);

    class_call(array_interpolate_spline(psp->ln_k,
                                        psp->ln_k_size,
                                        spectrum_cb_at_z,
                                        spline_cb,
                                        1,
                                        log(k),
                                        &last_index,
                                        pk_cb_tot,
                                        1,
                                        psp->error_message),
               psp->error_message,
               psp->error_message);

    *pk_cb_tot = exp(*pk_cb_tot);

  }

  free(spectrum_at_z);
  free(spectrum_cb_at_z);
  free(spline);
  if(pba->has_ncdm) free(spline_cb);

  return _SUCCESS_;

}


/**
 * Matter transfer functions \f$ T_i(k) \f$ for arbitrary redshift and for all
 * initial conditions.
 *
 * This routine evaluates the matter transfer functions at a given value of z by
 * interpolating in the pre-computed table (if several values of z have been stored)
 * or by directly reading it (if it only contains values at z=0 and we want \f$ T_i(k,z=0)\f$)
 *
 *
 * This function can be
 * called from whatever module at whatever time, provided that
 * spectra_init() has been called before, and spectra_free() has not
 * been called yet.
 *
 * @param pba        Input: pointer to background structure (used for converting z into tau)
 * @param psp        Input: pointer to spectra structure (containing pre-computed table)
 * @param z          Input: redshift
 * @param output     Output: matter transfer functions
 * @return the error status
 */

int spectra_tk_at_z(
                    struct background * pba,
                    struct spectra * psp,
                    double z,
                    double * output /* array with argument output[(index_k*psp->ic_size[index_md]+index_ic)*psp->tr_size+index_tr] (must be already allocated) */
                    ) {

  /** Summary: */

  /** - define local variables */

  int index_md;
  int last_index;
  int index_k;
  int index_tr;
  double tau,ln_tau;
  int index_ic;

  index_md = psp->index_md_scalars;

  /** - first step: convert z into ln(tau) */

  class_call(background_tau_of_z(pba,z,&tau),
             pba->error_message,
             psp->error_message);

  class_test(tau <= 0.,
             psp->error_message,
             "negative or null value of conformal time: cannot interpolate");

  ln_tau = log(tau);

  /** - second step: store the matter transfer functions in the output array */

  /** - --> (a) if only values at tau=tau_today are stored and we want \f$ T_i(k,z=0)\f$, no need to interpolate */

  if (psp->ln_tau_size == 1) {

    class_test(z != 0.,
               psp->error_message,
               "asked z=%e but only T_i(k,z=0) has been tabulated",z);

    for (index_k=0; index_k<psp->ln_k_size; index_k++)
      for (index_tr=0; index_tr<psp->tr_size; index_tr++)
        for (index_ic = 0; index_ic < psp->ic_size[index_md]; index_ic++)
          output[(index_k*psp->ic_size[index_md]+index_ic)*psp->tr_size+index_tr]
            = psp->matter_transfer[(index_k*psp->ic_size[index_md]+index_ic)*psp->tr_size+index_tr];

  }

  /** - --> (b) if several values of tau have been stored, use interpolation routine to get spectra at correct redshift */

  else {

    class_call(array_interpolate_spline(psp->ln_tau,
                                        psp->ln_tau_size,
                                        psp->matter_transfer,
                                        psp->ddmatter_transfer,
                                        psp->ic_size[index_md]*psp->tr_size*psp->ln_k_size,
                                        ln_tau,
                                        &last_index,
                                        output,
                                        psp->ic_size[index_md]*psp->tr_size*psp->ln_k_size,
                                        psp->error_message),
               psp->error_message,
               psp->error_message);

  }

  return _SUCCESS_;

}

/**
 * Matter transfer functions \f$ T_i(k)\f$ for arbitrary wavenumber, redshift
 * and initial condition.
 *
 * This routine evaluates the matter transfer functions at a given
 * value of k and z by interpolating in a table of all \f$ T_i(k,z)\f$'s
 * computed at this z by spectra_tk_at_z() (when kmin <= k <= kmax).
 * Returns an error when k<kmin or k > kmax.
 *
 * This function can be called from whatever module at whatever time,
 * provided that spectra_init() has been called before, and
 * spectra_free() has not been called yet.
 *
 * @param pba        Input: pointer to background structure (used for converting z into tau)
 * @param psp        Input: pointer to spectra structure (containing pre-computed table)
 * @param k          Input: wavenumber in 1/Mpc
 * @param z          Input: redshift
 * @param output     Output: matter transfer functions
 * @return the error status
 */

int spectra_tk_at_k_and_z(
                          struct background * pba,
                          struct spectra * psp,
                          double k,
                          double z,
                          double * output  /* array with argument output[index_ic*psp->tr_size+index_tr] (must be already allocated) */
                          ) {

  /** Summary: */

  /** - define local variables */

  int index_md;
  int last_index;
  double * tks_at_z;
  double * ddtks_at_z;

  index_md = psp->index_md_scalars;

  /** - check that k is in valid range [0:kmax] (the test for z will be done when calling spectra_tk_at_z()) */

  class_test((k < 0.) || (k > exp(psp->ln_k[psp->ln_k_size-1])),
             psp->error_message,
             "k=%e out of bounds [%e:%e]",k,0.,exp(psp->ln_k[psp->ln_k_size-1]));

  /** - compute T_i(k,z) */

  class_alloc(tks_at_z,
              psp->ln_k_size*psp->tr_size*psp->ic_size[index_md]*sizeof(double),
              psp->error_message);

  class_call(spectra_tk_at_z(pba,
                             psp,
                             z,
                             tks_at_z),
             psp->error_message,
             psp->error_message);

  /** - get its second derivatives w.r.t. k with spline, then interpolate */

  class_alloc(ddtks_at_z,
              psp->ln_k_size*psp->tr_size*psp->ic_size[index_md]*sizeof(double),
              psp->error_message);

  class_call(array_spline_table_lines(psp->ln_k,
                                      psp->ln_k_size,
                                      tks_at_z,
                                      psp->tr_size*psp->ic_size[index_md],
                                      ddtks_at_z,
                                      _SPLINE_NATURAL_,
                                      psp->error_message),
             psp->error_message,
             psp->error_message);

  class_call(array_interpolate_spline(psp->ln_k,
                                      psp->ln_k_size,
                                      tks_at_z,
                                      ddtks_at_z,
                                      psp->tr_size*psp->ic_size[index_md],
                                      log(k),
                                      &last_index,
                                      output,
                                      psp->tr_size*psp->ic_size[index_md],
                                      psp->error_message),
             psp->error_message,
             psp->error_message);

  free(tks_at_z);
  free(ddtks_at_z);

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

int spectra_init(
                 struct precision * ppr,
                 struct background * pba,
                 struct perturbs * ppt,
                 struct primordial * ppm,
                 struct nonlinear *pnl,
                 struct transfers * ptr,
                 struct spectra * psp
                 ) {

  /** Summary: */

  double TT_II,TT_RI,TT_RR;
  int l1,l2;

  /** - check that we really want to compute at least one spectrum */

  if ((ppt->has_cls == _FALSE_) &&
      (ppt->has_pk_matter == _FALSE_) &&
      (ppt->has_density_transfers == _FALSE_) &&
      (ppt->has_velocity_transfers == _FALSE_)) {
    psp->md_size = 0;
    if (psp->spectra_verbose > 0)
      printf("No spectra requested. Spectra module skipped.\n");
    return _SUCCESS_;
  }
  else {
    if (psp->spectra_verbose > 0)
      printf("Computing unlensed linear spectra\n");
  }

  /** - initialize indices and allocate some of the arrays in the
      spectra structure */

  class_call(spectra_indices(pba,ppt,ptr,ppm,psp),
             psp->error_message,
             psp->error_message);

  /** - deal with \f$ C_l\f$'s, if any */

  if (ppt->has_cls == _TRUE_) {

    class_call(spectra_cls(pba,ppt,ptr,ppm,psp),
               psp->error_message,
               psp->error_message);

  }
  else {
    psp->ct_size=0;
  }

  /** - deal with \f$ P(k,\tau)\f$ and \f$ T_i(k,\tau)\f$ */

  if ((ppt->has_pk_matter == _TRUE_) || (ppt->has_density_transfers == _TRUE_) || (ppt->has_velocity_transfers == _TRUE_)) {

    class_call(spectra_k_and_tau(pba,ppt,pnl,psp),
               psp->error_message,
               psp->error_message);

    if (ppt->has_pk_matter == _TRUE_) {

      class_call(spectra_pk(pba,ppt,ppm,pnl,psp),
                 psp->error_message,
                 psp->error_message);

    }
    else {
      psp->ln_pk=NULL;
      psp->ln_pk_cb=NULL;
    }

    if ((ppt->has_density_transfers == _TRUE_) || (ppt->has_velocity_transfers == _TRUE_)) {

      class_call(spectra_matter_transfers(pba,ppt,psp),
                 psp->error_message,
                 psp->error_message);
    }
    else {
      psp->matter_transfer=NULL;
    }

  }
  else {
    psp->ln_k_size=0;
  }

  /* if there is one isocurvature mode, compute and store in the psp
     structure the isocurvature contribution to some bandpowers in
     different ranges of l, and the contribution to the primordial
     spectrum at different wavenumbers (used in the Planck
     analysis) */

  if ((ppt->has_scalars == _TRUE_) && (ppt->has_cls == _TRUE_) && (ppt->ic_size[ppt->index_md_scalars] == 2)) {

    l1=2;
    l2=20;

    class_call(spectra_bandpower(psp,l1,l2,&TT_II,&TT_RI,&TT_RR),
               psp->error_message,
               psp->error_message);

    class_test(TT_II+TT_RI+TT_RR==0.,
               psp->error_message,
               "should never happen");

    psp->alpha_II_2_20=TT_II/(TT_II+TT_RI+TT_RR);
    psp->alpha_RI_2_20=TT_RI/(TT_II+TT_RI+TT_RR);
    psp->alpha_RR_2_20=TT_RR/(TT_II+TT_RI+TT_RR);

    l1=21;
    l2=200;

    class_call(spectra_bandpower(psp,l1,l2,&TT_II,&TT_RI,&TT_RR),
               psp->error_message,
               psp->error_message);

    class_test(TT_II+TT_RI+TT_RR==0.,
               psp->error_message,
               "should never happen");

    psp->alpha_II_21_200=TT_II/(TT_II+TT_RI+TT_RR);
    psp->alpha_RI_21_200=TT_RI/(TT_II+TT_RI+TT_RR);
    psp->alpha_RR_21_200=TT_RR/(TT_II+TT_RI+TT_RR);

    l1=201;
    l2=2500;

    class_call(spectra_bandpower(psp,l1,l2,&TT_II,&TT_RI,&TT_RR),
               psp->error_message,
               psp->error_message);

    class_test(TT_II+TT_RI+TT_RR==0.,
               psp->error_message,
               "should never happen");

    psp->alpha_II_201_2500=TT_II/(TT_II+TT_RI+TT_RR);
    psp->alpha_RI_201_2500=TT_RI/(TT_II+TT_RI+TT_RR);
    psp->alpha_RR_201_2500=TT_RR/(TT_II+TT_RI+TT_RR);

    l1=2;
    l2=2500;

    class_call(spectra_bandpower(psp,l1,l2,&TT_II,&TT_RI,&TT_RR),
               psp->error_message,
               psp->error_message);

    class_test(TT_II+TT_RI+TT_RR==0.,
               psp->error_message,
               "should never happen");

    psp->alpha_II_2_2500=TT_II/(TT_II+TT_RI+TT_RR);
    psp->alpha_RI_2_2500=TT_RI/(TT_II+TT_RI+TT_RR);
    psp->alpha_RR_2_2500=TT_RR/(TT_II+TT_RI+TT_RR);

    if (ppt->has_cdi==_TRUE_) {

      psp->alpha_kp=ppm->f_cdi*ppm->f_cdi
        /(1.+ppm->f_cdi*ppm->f_cdi);

      psp->alpha_k1=ppm->f_cdi*ppm->f_cdi*exp((ppm->n_cdi-ppm->n_s)*log(0.002/ppm->k_pivot))
        /(1.+ppm->f_cdi*ppm->f_cdi*exp((ppm->n_cdi-ppm->n_s)*log(0.002/ppm->k_pivot)));

      psp->alpha_k2=ppm->f_cdi*ppm->f_cdi*exp((ppm->n_cdi-ppm->n_s)*log(0.1/ppm->k_pivot))
        /(1.+ppm->f_cdi*ppm->f_cdi*exp((ppm->n_cdi-ppm->n_s)*log(0.1/ppm->k_pivot)));
    }

    if (ppt->has_nid==_TRUE_) {

      psp->alpha_kp=ppm->f_nid*ppm->f_nid
        /(1.+ppm->f_nid*ppm->f_nid);

      psp->alpha_k1=ppm->f_nid*ppm->f_nid*exp((ppm->n_nid-ppm->n_s)*log(0.002/ppm->k_pivot))
        /(1.+ppm->f_nid*ppm->f_nid*exp((ppm->n_nid-ppm->n_s)*log(0.002/ppm->k_pivot)));

      psp->alpha_k2=ppm->f_nid*ppm->f_nid*exp((ppm->n_nid-ppm->n_s)*log(0.1/ppm->k_pivot))
        /(1.+ppm->f_nid*ppm->f_nid*exp((ppm->n_nid-ppm->n_s)*log(0.1/ppm->k_pivot)));
    }

    if (ppt->has_niv==_TRUE_) {

      psp->alpha_kp=ppm->f_niv*ppm->f_niv
        /(1.+ppm->f_niv*ppm->f_niv);

      psp->alpha_k1=ppm->f_niv*ppm->f_niv*exp((ppm->n_niv-ppm->n_s)*log(0.002/ppm->k_pivot))
        /(1.+ppm->f_niv*ppm->f_niv*exp((ppm->n_niv-ppm->n_s)*log(0.002/ppm->k_pivot)));

      psp->alpha_k2=ppm->f_niv*ppm->f_niv*exp((ppm->n_niv-ppm->n_s)*log(0.1/ppm->k_pivot))
        /(1.+ppm->f_niv*ppm->f_niv*exp((ppm->n_niv-ppm->n_s)*log(0.1/ppm->k_pivot)));
    }
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

int spectra_free(
                 struct spectra * psp
                 ) {

  int index_md;

  if (psp->md_size > 0) {

    if (psp->ct_size > 0) {

      for (index_md = 0; index_md < psp->md_size; index_md++) {
        free(psp->l_max_ct[index_md]);
        free(psp->cl[index_md]);
        free(psp->ddcl[index_md]);
      }
      free(psp->l);
      free(psp->l_size);
      free(psp->l_max_ct);
      free(psp->l_max);
      free(psp->cl);
      free(psp->ddcl);
    }

    if (psp->ln_k_size > 0) {

      free(psp->ln_tau);
      free(psp->ln_k);

      if (psp->ln_pk != NULL) {

        free(psp->ln_pk);

        if (psp->ln_tau_size > 1) {
          free(psp->ddln_pk);
        }

        free(psp->ln_pk_l);

        if (psp->ln_tau_size > 1) {
          free(psp->ddln_pk_l);
        }

        if (psp->ln_pk_nl != NULL) {

          free(psp->ln_tau_nl);
          free(psp->ln_pk_nl);

          if (psp->ln_tau_nl_size > 1) {
            free(psp->ddln_pk_nl);
          }
        }

      }

      if (psp->ln_pk_cb != NULL) {

        free(psp->ln_pk_cb);

        if (psp->ln_tau_size > 1) {
          free(psp->ddln_pk_cb);
        }

        free(psp->ln_pk_cb_l);

        if (psp->ln_tau_size > 1) {
          free(psp->ddln_pk_cb_l);
        }

        if (psp->ln_pk_cb_nl != NULL) {

          free(psp->ln_pk_cb_nl);

          if (psp->ln_tau_nl_size > 1) {
            free(psp->ddln_pk_cb_nl);
          }
        }

      }

      if (psp->matter_transfer != NULL) {

        free(psp->matter_transfer);
        if (psp->ln_tau_size > 1) {
          free(psp->ddmatter_transfer);
        }
      }
    }
  }

  for (index_md=0; index_md < psp->md_size; index_md++)
    free(psp->is_non_zero[index_md]);
  free(psp->is_non_zero);
  free(psp->ic_size);
  free(psp->ic_ic_size);

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

int spectra_indices(
                    struct background * pba,
                    struct perturbs * ppt,
                    struct transfers * ptr,
                    struct primordial * ppm,
                    struct spectra * psp
                    ){

  int index_ct;
  int index_md;
  int index_ic1_ic2;
  int index_tr;

  psp->md_size = ppt->md_size;
  if (ppt->has_scalars == _TRUE_)
    psp->index_md_scalars = ppt->index_md_scalars;

  class_alloc(psp->ic_size,
              sizeof(int)*psp->md_size,
              psp->error_message);

  class_alloc(psp->ic_ic_size,
              sizeof(int)*psp->md_size,
              psp->error_message);

  class_alloc(psp->is_non_zero,
              sizeof(short *)*psp->md_size,
              psp->error_message);

  for (index_md=0; index_md < psp->md_size; index_md++) {
    psp->ic_size[index_md] = ppm->ic_size[index_md];
    psp->ic_ic_size[index_md] = ppm->ic_ic_size[index_md];
    class_alloc(psp->is_non_zero[index_md],
                sizeof(short)*psp->ic_ic_size[index_md],
                psp->error_message);
    for (index_ic1_ic2=0; index_ic1_ic2 < psp->ic_ic_size[index_md]; index_ic1_ic2++)
      psp->is_non_zero[index_md][index_ic1_ic2] = ppm->is_non_zero[index_md][index_ic1_ic2];
  }

  if (ppt->has_cls == _TRUE_) {

    /* types of C_l's relevant for both scalars and tensors: TT, EE, TE */

    index_ct=0;

    if (ppt->has_cl_cmb_temperature == _TRUE_) {
      psp->has_tt = _TRUE_;
      psp->index_ct_tt=index_ct;
      index_ct++;
    }
    else {
      psp->has_tt = _FALSE_;
    }

    if (ppt->has_cl_cmb_polarization == _TRUE_) {
      psp->has_ee = _TRUE_;
      psp->index_ct_ee=index_ct;
      index_ct++;
    }
    else {
      psp->has_ee = _FALSE_;
    }

    if ((ppt->has_cl_cmb_temperature == _TRUE_) &&
        (ppt->has_cl_cmb_polarization == _TRUE_)) {
      psp->has_te = _TRUE_;
      psp->index_ct_te=index_ct;
      index_ct++;
    }
    else {
      psp->has_te = _FALSE_;
    }

    if (ppt->has_cl_cmb_polarization == _TRUE_) {
      psp->has_bb = _TRUE_;
      psp->index_ct_bb=index_ct;
      index_ct++;
    }
    else {
      psp->has_bb = _FALSE_;
    }

    /* types of C_l's relevant only for scalars: phi-phi, T-phi, E-phi, d-d, T-d */

    if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
      psp->has_pp = _TRUE_;
      psp->index_ct_pp=index_ct;
      index_ct++;
    }
    else {
      psp->has_pp = _FALSE_;
    }

    if ((ppt->has_cl_cmb_temperature == _TRUE_) && (ppt->has_cl_cmb_lensing_potential == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
      psp->has_tp = _TRUE_;
      psp->index_ct_tp=index_ct;
      index_ct++;
    }
    else {
      psp->has_tp = _FALSE_;
    }

    psp->ct_size = index_ct;

    if ((ppt->has_cl_cmb_polarization == _TRUE_) && (ppt->has_cl_cmb_lensing_potential == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
      psp->has_ep = _TRUE_;
      psp->index_ct_ep=index_ct;
      index_ct++;
    }
    else {
      psp->has_ep = _FALSE_;
    }

    if ((ppt->has_scalars == _TRUE_) &&
        ((ppt->has_cl_number_count == _TRUE_) || (ppt->has_cl_lensing_potential == _TRUE_)))
      psp->d_size=ppt->selection_num;
    else
      psp->d_size=0;

    if ((ppt->has_cl_number_count == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
      psp->has_dd = _TRUE_;
      psp->index_ct_dd=index_ct;
      index_ct+=(psp->d_size*(psp->d_size+1)-(psp->d_size-psp->non_diag)*(psp->d_size-1-psp->non_diag))/2;
    }
    else {
      psp->has_dd = _FALSE_;
    }

    /* the computation of C_l^Td would require a very good sampling of
       transfer functions over a wide range, and a huge computation
       time. In the current version, we prefer to switch it off, rather
       than either slowing down the code considerably, or producing
       very inaccurate spectra.

       if ((ppt->has_cl_cmb_temperature == _TRUE_) && (ppt->has_cl_number_count == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
       psp->has_td = _TRUE_;
       psp->index_ct_td=index_ct;
       index_ct+=psp->d_size;
       }
       else {
       psp->has_td = _FALSE_;
       }
    */
    psp->has_td = _FALSE_;

    if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (ppt->has_cl_number_count == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
      psp->has_pd = _TRUE_;
      psp->index_ct_pd=index_ct;
      index_ct+=psp->d_size;
    }
    else {
      psp->has_pd = _FALSE_;
    }

    psp->has_td = _FALSE_;

    if ((ppt->has_cl_lensing_potential == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
      psp->has_ll = _TRUE_;
      psp->index_ct_ll=index_ct;
      index_ct+=(psp->d_size*(psp->d_size+1)-(psp->d_size-psp->non_diag)*(psp->d_size-1-psp->non_diag))/2;
    }
    else {
      psp->has_ll = _FALSE_;
    }

    /* the computation of C_l^Tl would require a very good sampling of
       transfer functions over a wide range, and a huge computation
       time. In the current version, we prefer to switch it off, rather
       than either slowing down the code considerably, or producing
       very inaccurate spectra.

       if ((ppt->has_cl_cmb_temperature == _TRUE_) && (ppt->has_cl_lensing_potential == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
       psp->has_tl = _TRUE_;
       psp->index_ct_tl=index_ct;
       index_ct+=psp->d_size;
       }
       else {
       psp->has_tl = _FALSE_;
       }
    */
    psp->has_tl = _FALSE_;

    if ((ppt->has_cl_number_count == _TRUE_) && (ppt->has_cl_lensing_potential == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
      psp->has_dl = _TRUE_;
      psp->index_ct_dl=index_ct;
      index_ct += psp->d_size*psp->d_size - (psp->d_size-psp->non_diag)*(psp->d_size-1-psp->non_diag);
    }
    else {
      psp->has_dl = _FALSE_;
    }

    psp->ct_size = index_ct;

    /* infer from input quantities the l_max for each mode and type,
       l_max_ct[index_md][index_type].  Maximize it over index_ct, and
       then over index_md. */

    class_alloc(psp->l_max,sizeof(int*)*psp->md_size,psp->error_message);
    class_alloc(psp->l_max_ct,sizeof(int*)*psp->md_size,psp->error_message);
    for (index_md=0; index_md<psp->md_size; index_md++) {
      class_calloc(psp->l_max_ct[index_md],psp->ct_size,sizeof(int),psp->error_message);
    }

    if (ppt->has_scalars == _TRUE_) {

      /* spectra computed up to l_scalar_max */

      if (psp->has_tt == _TRUE_) psp->l_max_ct[ppt->index_md_scalars][psp->index_ct_tt] = ppt->l_scalar_max;
      if (psp->has_ee == _TRUE_) psp->l_max_ct[ppt->index_md_scalars][psp->index_ct_ee] = ppt->l_scalar_max;
      if (psp->has_te == _TRUE_) psp->l_max_ct[ppt->index_md_scalars][psp->index_ct_te] = ppt->l_scalar_max;
      if (psp->has_pp == _TRUE_) psp->l_max_ct[ppt->index_md_scalars][psp->index_ct_pp] = ppt->l_scalar_max;
      if (psp->has_tp == _TRUE_) psp->l_max_ct[ppt->index_md_scalars][psp->index_ct_tp] = ppt->l_scalar_max;
      if (psp->has_ep == _TRUE_) psp->l_max_ct[ppt->index_md_scalars][psp->index_ct_ep] = ppt->l_scalar_max;

      /* spectra computed up to l_lss_max */

      if (psp->has_dd == _TRUE_)
        for (index_ct=psp->index_ct_dd;
             index_ct<psp->index_ct_dd+(psp->d_size*(psp->d_size+1)-(psp->d_size-psp->non_diag)*(psp->d_size-1-psp->non_diag))/2;
             index_ct++)
          psp->l_max_ct[ppt->index_md_scalars][index_ct] = ppt->l_lss_max;

      if (psp->has_td == _TRUE_)
        for (index_ct=psp->index_ct_td;
             index_ct<psp->index_ct_td+psp->d_size;
             index_ct++)
          psp->l_max_ct[ppt->index_md_scalars][index_ct] = MIN(ppt->l_scalar_max,ppt->l_lss_max);

      if (psp->has_pd == _TRUE_)
        for (index_ct=psp->index_ct_pd;
             index_ct<psp->index_ct_pd+psp->d_size;
             index_ct++)
          psp->l_max_ct[ppt->index_md_scalars][index_ct] = MIN(ppt->l_scalar_max,ppt->l_lss_max);

      if (psp->has_ll == _TRUE_)
        for (index_ct=psp->index_ct_ll;
             index_ct<psp->index_ct_ll+(psp->d_size*(psp->d_size+1)-(psp->d_size-psp->non_diag)*(psp->d_size-1-psp->non_diag))/2;
             index_ct++)
          psp->l_max_ct[ppt->index_md_scalars][index_ct] = ppt->l_lss_max;

      if (psp->has_tl == _TRUE_)
        for (index_ct=psp->index_ct_tl;
             index_ct<psp->index_ct_tl+psp->d_size;
             index_ct++)
          psp->l_max_ct[ppt->index_md_scalars][index_ct] = MIN(ppt->l_scalar_max,ppt->l_lss_max);

      if (psp->has_dl == _TRUE_)
        for (index_ct=psp->index_ct_dl;
             index_ct < psp->index_ct_dl+(psp->d_size*psp->d_size - (psp->d_size-psp->non_diag)*(psp->d_size-1-psp->non_diag));
             index_ct++)
          psp->l_max_ct[ppt->index_md_scalars][index_ct] = ppt->l_lss_max;

    }
    if (ppt->has_tensors == _TRUE_) {

      /* spectra computed up to l_tensor_max */

      if (psp->has_tt == _TRUE_) psp->l_max_ct[ppt->index_md_tensors][psp->index_ct_tt] = ppt->l_tensor_max;
      if (psp->has_ee == _TRUE_) psp->l_max_ct[ppt->index_md_tensors][psp->index_ct_ee] = ppt->l_tensor_max;
      if (psp->has_te == _TRUE_) psp->l_max_ct[ppt->index_md_tensors][psp->index_ct_te] = ppt->l_tensor_max;
      if (psp->has_bb == _TRUE_) psp->l_max_ct[ppt->index_md_tensors][psp->index_ct_bb] = ppt->l_tensor_max;
    }

    /* maximizations */
    psp->l_max_tot = 0.;
    for (index_md=0; index_md < psp->md_size; index_md++) {
      psp->l_max[index_md] = 0.;
      for (index_ct=0.; index_ct<psp->ct_size; index_ct++)
        psp->l_max[index_md] = MAX(psp->l_max[index_md],psp->l_max_ct[index_md][index_ct]);
      psp->l_max_tot = MAX(psp->l_max_tot,psp->l_max[index_md]);
    }
  }

  /* indices for species associated with a matter transfer function in Fourier space */

  index_tr=0;
  class_define_index(psp->index_tr_delta_g,ppt->has_source_delta_g,index_tr,1);
  class_define_index(psp->index_tr_delta_b,ppt->has_source_delta_b,index_tr,1);
  class_define_index(psp->index_tr_delta_cdm,ppt->has_source_delta_cdm,index_tr,1);
  class_define_index(psp->index_tr_delta_dcdm,ppt->has_source_delta_dcdm,index_tr,1);
  class_define_index(psp->index_tr_delta_scf,ppt->has_source_delta_scf,index_tr,1);
  class_define_index(psp->index_tr_delta_fld,ppt->has_source_delta_fld,index_tr,1);
  class_define_index(psp->index_tr_delta_ur,ppt->has_source_delta_ur,index_tr,1);
  class_define_index(psp->index_tr_delta_dr,ppt->has_source_delta_dr,index_tr,1);
  class_define_index(psp->index_tr_delta_ncdm1,ppt->has_source_delta_ncdm,index_tr,pba->N_ncdm);
  class_define_index(psp->index_tr_delta_tot,ppt->has_density_transfers,index_tr,1);
  class_define_index(psp->index_tr_phi,ppt->has_source_phi,index_tr,1);
  class_define_index(psp->index_tr_psi,ppt->has_source_psi,index_tr,1);
  class_define_index(psp->index_tr_phi,ppt->has_source_phi_prime,index_tr,1);
  class_define_index(psp->index_tr_h,ppt->has_source_h,index_tr,1);
  class_define_index(psp->index_tr_h_prime,ppt->has_source_h_prime,index_tr,1);
  class_define_index(psp->index_tr_eta,ppt->has_source_eta,index_tr,1);
  class_define_index(psp->index_tr_eta_prime,ppt->has_source_eta_prime,index_tr,1);

  /* indices for species associated with a velocity transfer function in Fourier space */

  class_define_index(psp->index_tr_theta_g,ppt->has_source_theta_g,index_tr,1);
  class_define_index(psp->index_tr_theta_b,ppt->has_source_theta_b,index_tr,1);
  class_define_index(psp->index_tr_theta_cdm,ppt->has_source_theta_cdm,index_tr,1);
  class_define_index(psp->index_tr_theta_dcdm,ppt->has_source_theta_dcdm,index_tr,1);
  class_define_index(psp->index_tr_theta_scf,ppt->has_source_theta_scf,index_tr,1);
  class_define_index(psp->index_tr_theta_fld,ppt->has_source_theta_fld,index_tr,1);
  class_define_index(psp->index_tr_theta_ur,ppt->has_source_theta_ur,index_tr,1);
  class_define_index(psp->index_tr_theta_dr,ppt->has_source_theta_dr,index_tr,1);
  class_define_index(psp->index_tr_theta_ncdm1,ppt->has_source_theta_ncdm,index_tr,pba->N_ncdm);
  class_define_index(psp->index_tr_theta_tot,ppt->has_velocity_transfers,index_tr,1);

  psp->tr_size = index_tr;

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

int spectra_cls(
                struct background * pba,
                struct perturbs * ppt,
                struct transfers * ptr,
                struct primordial * ppm,
                struct spectra * psp
                ) {

  /** Summary: */

  /** - define local variables */

  int index_md;
  int index_ic1,index_ic2,index_ic1_ic2;
  int index_l;
  int index_ct;
  int cl_integrand_num_columns;

  double * cl_integrand; /* array with argument cl_integrand[index_k*cl_integrand_num_columns+1+psp->index_ct] */
  double * transfer_ic1; /* array with argument transfer_ic1[index_tt] */
  double * transfer_ic2; /* idem */
  double * primordial_pk;  /* array with argument primordial_pk[index_ic_ic]*/

  /* This code can be optionally compiled with the openmp option for parallel computation.
     Inside parallel regions, the use of the command "return" is forbidden.
     For error management, instead of "return _FAILURE_", we will set the variable below
     to "abort = _TRUE_". This will lead to a "return _FAILURE_" jus after leaving the
     parallel region. */
  int abort;

#ifdef _OPENMP
  /* instrumentation times */
  double tstart, tstop;
#endif

  /** - allocate pointers to arrays where results will be stored */

  class_alloc(psp->l_size,sizeof(int)*psp->md_size,psp->error_message);
  class_alloc(psp->cl,sizeof(double *)*psp->md_size,psp->error_message);
  class_alloc(psp->ddcl,sizeof(double *)*psp->md_size,psp->error_message);

  psp->l_size_max = ptr->l_size_max;
  class_alloc(psp->l,sizeof(double)*psp->l_size_max,psp->error_message);

  /** - store values of l */
  for (index_l=0; index_l < psp->l_size_max; index_l++) {
    psp->l[index_l] = (double)ptr->l[index_l];
  }

  /** - loop over modes (scalar, tensors, etc). For each mode: */

  for (index_md = 0; index_md < psp->md_size; index_md++) {

    /** - --> (a) store number of l values for this mode */

    psp->l_size[index_md] = ptr->l_size[index_md];

    /** - --> (b) allocate arrays where results will be stored */

    class_alloc(psp->cl[index_md],sizeof(double)*psp->l_size[index_md]*psp->ct_size*psp->ic_ic_size[index_md],psp->error_message);
    class_alloc(psp->ddcl[index_md],sizeof(double)*psp->l_size[index_md]*psp->ct_size*psp->ic_ic_size[index_md],psp->error_message);
    cl_integrand_num_columns = 1+psp->ct_size*2; /* one for k, ct_size for each type, ct_size for each second derivative of each type */

    /** - --> (c) loop over initial conditions */

    for (index_ic1 = 0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
      for (index_ic2 = index_ic1; index_ic2 < psp->ic_size[index_md]; index_ic2++) {
        index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md]);

        /* non-diagonal coefficients should be computed only if non-zero correlation */
        if (psp->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {

          /* initialize error management flag */
          abort = _FALSE_;

          /* beginning of parallel region */

#pragma omp parallel                                                    \
  shared(ptr,ppm,index_md,psp,ppt,cl_integrand_num_columns,index_ic1,index_ic2,abort) \
  private(tstart,cl_integrand,primordial_pk,transfer_ic1,transfer_ic2,index_l,tstop)

          {

#ifdef _OPENMP
            tstart = omp_get_wtime();
#endif

            class_alloc_parallel(cl_integrand,
                                 ptr->q_size*cl_integrand_num_columns*sizeof(double),
                                 psp->error_message);

            class_alloc_parallel(primordial_pk,
                                 psp->ic_ic_size[index_md]*sizeof(double),
                                 psp->error_message);

            class_alloc_parallel(transfer_ic1,
                                 ptr->tt_size[index_md]*sizeof(double),
                                 psp->error_message);

            class_alloc_parallel(transfer_ic2,
                                 ptr->tt_size[index_md]*sizeof(double),
                                 psp->error_message);

#pragma omp for schedule (dynamic)

            /** - ---> loop over l values defined in the transfer module.
                For each l, compute the \f$ C_l\f$'s for all types (TT, TE, ...)
                by convolving primordial spectra with transfer  functions.
                This elementary task is assigned to spectra_compute_cl() */

            for (index_l=0; index_l < ptr->l_size[index_md]; index_l++) {

#pragma omp flush(abort)

              class_call_parallel(spectra_compute_cl(pba,
                                                     ppt,
                                                     ptr,
                                                     ppm,
                                                     psp,
                                                     index_md,
                                                     index_ic1,
                                                     index_ic2,
                                                     index_l,
                                                     cl_integrand_num_columns,
                                                     cl_integrand,
                                                     primordial_pk,
                                                     transfer_ic1,
                                                     transfer_ic2),
                                  psp->error_message,
                                  psp->error_message);

            } /* end of loop over l */

#ifdef _OPENMP
            tstop = omp_get_wtime();
            if (psp->spectra_verbose > 1)
              printf("In %s: time spent in parallel region (loop over l's) = %e s for thread %d\n",
                     __func__,tstop-tstart,omp_get_thread_num());
#endif
            free(cl_integrand);

            free(primordial_pk);

            free(transfer_ic1);

            free(transfer_ic2);

          } /* end of parallel region */

          if (abort == _TRUE_) return _FAILURE_;

        }
        else {

          /* set non-diagonal coefficients to zero if pair of ic's uncorrelated */

          for (index_l=0; index_l < ptr->l_size[index_md]; index_l++) {
            for (index_ct=0; index_ct<psp->ct_size; index_ct++) {
              psp->cl[index_md]
                [(index_l * psp->ic_ic_size[index_md] + index_ic1_ic2) * psp->ct_size + index_ct]
                = 0.;
            }
          }
        }
      }
    }

    /** - --> (d) now that for a given mode, all possible \f$ C_l\f$'s have been computed,
        compute second derivative of the array in which they are stored,
        in view of spline interpolation. */

    class_call(array_spline_table_lines(psp->l,
                                        psp->l_size[index_md],
                                        psp->cl[index_md],
                                        psp->ic_ic_size[index_md]*psp->ct_size,
                                        psp->ddcl[index_md],
                                        _SPLINE_EST_DERIV_,
                                        psp->error_message),
               psp->error_message,
               psp->error_message);
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

int spectra_compute_cl(
                       struct background * pba,
                       struct perturbs * ppt,
                       struct transfers * ptr,
                       struct primordial * ppm,
                       struct spectra * psp,
                       int index_md,
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

  index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md]);

  if (ppt->has_cl_number_count == _TRUE_) {
    class_alloc(transfer_ic1_nc,psp->d_size*sizeof(double),psp->error_message);
    class_alloc(transfer_ic2_nc,psp->d_size*sizeof(double),psp->error_message);
  }

  for (index_q=0; index_q < ptr->q_size; index_q++) {

    //q = ptr->q[index_q];
    k = ptr->k[index_md][index_q];

    cl_integrand[index_q*cl_integrand_num_columns+0] = k;

    class_call(primordial_spectrum_at_k(ppm,index_md,linear,k,primordial_pk),
               ppm->error_message,
               psp->error_message);

    /* above routine checks that k>0: no possible division by zero below */

    for (index_tt=0; index_tt < ptr->tt_size[index_md]; index_tt++) {

      transfer_ic1[index_tt] =
        ptr->transfer[index_md]
        [((index_ic1 * ptr->tt_size[index_md] + index_tt)
          * ptr->l_size[index_md] + index_l)
         * ptr->q_size + index_q];

      if (index_ic1 == index_ic2) {
        transfer_ic2[index_tt] = transfer_ic1[index_tt];
      }
      else {
        transfer_ic2[index_tt] = ptr->transfer[index_md]
          [((index_ic2 * ptr->tt_size[index_md] + index_tt)
            * ptr->l_size[index_md] + index_l)
           * ptr->q_size + index_q];
      }
    }

    /* define combinations of transfer functions */

    if (ppt->has_cl_cmb_temperature == _TRUE_) {

      if (_scalars_) {

        transfer_ic1_temp = transfer_ic1[ptr->index_tt_t0] + transfer_ic1[ptr->index_tt_t1] + transfer_ic1[ptr->index_tt_t2];
        transfer_ic2_temp = transfer_ic2[ptr->index_tt_t0] + transfer_ic2[ptr->index_tt_t1] + transfer_ic2[ptr->index_tt_t2];

      }

      if (_vectors_) {

        transfer_ic1_temp = transfer_ic1[ptr->index_tt_t1] + transfer_ic1[ptr->index_tt_t2];
        transfer_ic2_temp = transfer_ic2[ptr->index_tt_t1] + transfer_ic2[ptr->index_tt_t2];

      }

      if (_tensors_) {

        transfer_ic1_temp = transfer_ic1[ptr->index_tt_t2];
        transfer_ic2_temp = transfer_ic2[ptr->index_tt_t2];

      }
    }

    if (ppt->has_cl_number_count == _TRUE_) {

      for (index_d1=0; index_d1<psp->d_size; index_d1++) {

        transfer_ic1_nc[index_d1] = 0.;
        transfer_ic2_nc[index_d1] = 0.;

        if (ppt->has_nc_density == _TRUE_) {
          transfer_ic1_nc[index_d1] += transfer_ic1[ptr->index_tt_density+index_d1];
          transfer_ic2_nc[index_d1] += transfer_ic2[ptr->index_tt_density+index_d1];
        }

        if (ppt->has_nc_rsd     == _TRUE_) {
          transfer_ic1_nc[index_d1]
            += transfer_ic1[ptr->index_tt_rsd+index_d1]
            + transfer_ic1[ptr->index_tt_d0+index_d1]
            + transfer_ic1[ptr->index_tt_d1+index_d1];
          transfer_ic2_nc[index_d1]
            += transfer_ic2[ptr->index_tt_rsd+index_d1]
            + transfer_ic2[ptr->index_tt_d0+index_d1]
            + transfer_ic2[ptr->index_tt_d1+index_d1];
        }

        if (ppt->has_nc_lens == _TRUE_) {
          transfer_ic1_nc[index_d1] +=
            psp->l[index_l]*(psp->l[index_l]+1.)*transfer_ic1[ptr->index_tt_nc_lens+index_d1];
          transfer_ic2_nc[index_d1] +=
            psp->l[index_l]*(psp->l[index_l]+1.)*transfer_ic2[ptr->index_tt_nc_lens+index_d1];
        }

        if (ppt->has_nc_gr == _TRUE_) {
          transfer_ic1_nc[index_d1]
            += transfer_ic1[ptr->index_tt_nc_g1+index_d1]
            + transfer_ic1[ptr->index_tt_nc_g2+index_d1]
            + transfer_ic1[ptr->index_tt_nc_g3+index_d1]
            + transfer_ic1[ptr->index_tt_nc_g4+index_d1]
            + transfer_ic1[ptr->index_tt_nc_g5+index_d1];
          transfer_ic2_nc[index_d1]
            += transfer_ic2[ptr->index_tt_nc_g1+index_d1]
            + transfer_ic2[ptr->index_tt_nc_g2+index_d1]
            + transfer_ic2[ptr->index_tt_nc_g3+index_d1]
            + transfer_ic2[ptr->index_tt_nc_g4+index_d1]
            + transfer_ic2[ptr->index_tt_nc_g5+index_d1];
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

    if (psp->has_tt == _TRUE_)
      cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_tt]=
        primordial_pk[index_ic1_ic2]
        * transfer_ic1_temp
        * transfer_ic2_temp
        * factor;

    if (psp->has_ee == _TRUE_)
      cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_ee]=
        primordial_pk[index_ic1_ic2]
        * transfer_ic1[ptr->index_tt_e]
        * transfer_ic2[ptr->index_tt_e]
        * factor;

    if (psp->has_te == _TRUE_)
      cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_te]=
        primordial_pk[index_ic1_ic2]
        * 0.5*(transfer_ic1_temp * transfer_ic2[ptr->index_tt_e] +
               transfer_ic1[ptr->index_tt_e] * transfer_ic2_temp)
        * factor;

    if (_tensors_ && (psp->has_bb == _TRUE_))
      cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_bb]=
        primordial_pk[index_ic1_ic2]
        * transfer_ic1[ptr->index_tt_b]
        * transfer_ic2[ptr->index_tt_b]
        * factor;

    if (_scalars_ && (psp->has_pp == _TRUE_))
      cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_pp]=
        primordial_pk[index_ic1_ic2]
        * transfer_ic1[ptr->index_tt_lcmb]
        * transfer_ic2[ptr->index_tt_lcmb]
        * factor;

    if (_scalars_ && (psp->has_tp == _TRUE_))
      cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_tp]=
        primordial_pk[index_ic1_ic2]
        * 0.5*(transfer_ic1_temp * transfer_ic2[ptr->index_tt_lcmb] +
               transfer_ic1[ptr->index_tt_lcmb] * transfer_ic2_temp)
        * factor;

    if (_scalars_ && (psp->has_ep == _TRUE_))
      cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_ep]=
        primordial_pk[index_ic1_ic2]
        * 0.5*(transfer_ic1[ptr->index_tt_e] * transfer_ic2[ptr->index_tt_lcmb] +
               transfer_ic1[ptr->index_tt_lcmb] * transfer_ic2[ptr->index_tt_e])
        * factor;

    if (_scalars_ && (psp->has_dd == _TRUE_)) {
      index_ct=0;
      for (index_d1=0; index_d1<psp->d_size; index_d1++) {
        for (index_d2=index_d1; index_d2<=MIN(index_d1+psp->non_diag,psp->d_size-1); index_d2++) {
          cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_dd+index_ct]=
            primordial_pk[index_ic1_ic2]
            * transfer_ic1_nc[index_d1]
            * transfer_ic2_nc[index_d2]
            * factor;
          index_ct++;
        }
      }
    }

    if (_scalars_ && (psp->has_td == _TRUE_)) {
      for (index_d1=0; index_d1<psp->d_size; index_d1++) {
        cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_td+index_d1]=
          primordial_pk[index_ic1_ic2]
          * 0.5*(transfer_ic1_temp * transfer_ic2_nc[index_d1] +
                 transfer_ic1_nc[index_d1] * transfer_ic2_temp)
          * factor;
      }
    }

    if (_scalars_ && (psp->has_pd == _TRUE_)) {
      for (index_d1=0; index_d1<psp->d_size; index_d1++) {
        cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_pd+index_d1]=
          primordial_pk[index_ic1_ic2]
          * 0.5*(transfer_ic1[ptr->index_tt_lcmb] * transfer_ic2_nc[index_d1] +
                 transfer_ic1_nc[index_d1] * transfer_ic2[ptr->index_tt_lcmb])
          * factor;
      }
    }

    if (_scalars_ && (psp->has_ll == _TRUE_)) {
      index_ct=0;
      for (index_d1=0; index_d1<psp->d_size; index_d1++) {
        for (index_d2=index_d1; index_d2<=MIN(index_d1+psp->non_diag,psp->d_size-1); index_d2++) {
          cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_ll+index_ct]=
            primordial_pk[index_ic1_ic2]
            * transfer_ic1[ptr->index_tt_lensing+index_d1]
            * transfer_ic2[ptr->index_tt_lensing+index_d2]
            * factor;
          index_ct++;
        }
      }
    }

    if (_scalars_ && (psp->has_tl == _TRUE_)) {
      for (index_d1=0; index_d1<psp->d_size; index_d1++) {
        cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_tl+index_d1]=
          primordial_pk[index_ic1_ic2]
          * 0.5*(transfer_ic1_temp * transfer_ic2[ptr->index_tt_lensing+index_d1] +
                 transfer_ic1[ptr->index_tt_lensing+index_d1] * transfer_ic2_temp)
          * factor;
      }
    }

    if (_scalars_ && (psp->has_dl == _TRUE_)) {
      index_ct=0;
      for (index_d1=0; index_d1<psp->d_size; index_d1++) {
        for (index_d2=MAX(index_d1-psp->non_diag,0); index_d2<=MIN(index_d1+psp->non_diag,psp->d_size-1); index_d2++) {
          cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_dl+index_ct]=
            primordial_pk[index_ic1_ic2]
            * transfer_ic1_nc[index_d1] * transfer_ic2[ptr->index_tt_lensing+index_d2]
            * factor;
          index_ct++;
        }
      }
    }
  }

  for (index_ct=0; index_ct<psp->ct_size; index_ct++) {

    /* treat null spectra (C_l^BB of scalars, C_l^pp of tensors, etc. */

    if ((_scalars_ && (psp->has_bb == _TRUE_) && (index_ct == psp->index_ct_bb)) ||
        (_tensors_ && (psp->has_pp == _TRUE_) && (index_ct == psp->index_ct_pp)) ||
        (_tensors_ && (psp->has_tp == _TRUE_) && (index_ct == psp->index_ct_tp)) ||
        (_tensors_ && (psp->has_ep == _TRUE_) && (index_ct == psp->index_ct_ep)) ||
        (_tensors_ && (psp->has_dd == _TRUE_) && (index_ct == psp->index_ct_dd)) ||
        (_tensors_ && (psp->has_td == _TRUE_) && (index_ct == psp->index_ct_td)) ||
        (_tensors_ && (psp->has_pd == _TRUE_) && (index_ct == psp->index_ct_pd)) ||
        (_tensors_ && (psp->has_ll == _TRUE_) && (index_ct == psp->index_ct_ll)) ||
        (_tensors_ && (psp->has_tl == _TRUE_) && (index_ct == psp->index_ct_tl)) ||
        (_tensors_ && (psp->has_dl == _TRUE_) && (index_ct == psp->index_ct_dl))
        ) {

      psp->cl[index_md]
        [(index_l * psp->ic_ic_size[index_md] + index_ic1_ic2) * psp->ct_size + index_ct] = 0.;

    }
    /* for non-zero spectra, integrate over q */
    else {

      /* spline the integrand over the whole range of k's */

      class_call(array_spline(cl_integrand,
                              cl_integrand_num_columns,
                              ptr->q_size,
                              0,
                              1+index_ct,
                              1+psp->ct_size+index_ct,
                              _SPLINE_EST_DERIV_,
                              psp->error_message),
                 psp->error_message,
                 psp->error_message);

      /* Technical point: we will now do a spline integral over the
         whole range of k's, excepted in the closed (K>0) case. In
         that case, it is a bad idea to spline over the values of k
         corresponding to nu<nu_flat_approximation. In this region, nu
         values are integer values, so the steps dq and dk have some
         discrete jumps. This makes the spline routine less accurate
         than a trapezoidal integral with finer sampling. So, in the
         closed case, we set index_q_spline to
         ptr->index_q_flat_approximation, to tell the integration
         routine that below this index, it should treat the integral
         as a trapezoidal one. For testing, one is free to set
         index_q_spline to 0, to enforce spline integration
         everywhere, or to (ptr->q_size-1), to enforce trapezoidal
         integration everywhere. */

      if (pba->sgnK == 1) {
        index_q_spline = ptr->index_q_flat_approximation;
      }

      class_call(array_integrate_all_trapzd_or_spline(cl_integrand,
                                                      cl_integrand_num_columns,
                                                      ptr->q_size,
                                                      index_q_spline,
                                                      0,
                                                      1+index_ct,
                                                      1+psp->ct_size+index_ct,
                                                      &clvalue,
                                                      psp->error_message),
                 psp->error_message,
                 psp->error_message);

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
        clvalue += cl_integrand[1+index_ct] * ptr->q[0]/ptr->k[0][0]*sqrt(pba->K)/2.;
      }

      /* we have the correct C_l now. We can store it in the transfer structure. */

      psp->cl[index_md]
        [(index_l * psp->ic_ic_size[index_md] + index_ic1_ic2) * psp->ct_size + index_ct]
        = clvalue;

    }
  }

  if (ppt->has_cl_number_count == _TRUE_) {
    free(transfer_ic1_nc);
    free(transfer_ic2_nc);
  }

  return _SUCCESS_;

}

/**
 * This routine computes the values of k and tau at which the matter
 * power spectra \f$ P(k,\tau)\f$ and the matter transfer functions \f$ T_i(k,\tau)\f$
 * will be stored.
 *
 * @param pba Input: pointer to background structure (for z to tau conversion)
 * @param ppt Input: pointer to perturbation structure (contain source functions)
 * @param pnl Input: pointer to nonlinear structure (contain nonlinear corrections)
 * @param psp Input/Output: pointer to spectra structure
 * @return the error status
 */

int spectra_k_and_tau(
                      struct background * pba,
                      struct perturbs * ppt,
                      struct nonlinear * pnl,
                      struct spectra * psp
                      ) {

  /** Summary: */

  /** - define local variables */

  int index_k;
  int index_tau;
  int index_tau_min_nl;
  double tau_min;

  /** - check the presence of scalar modes */

  class_test((ppt->has_scalars == _FALSE_),
             psp->error_message,
             "you cannot ask for matter power spectrum since you turned off scalar modes");

  /** - check the maximum redshift z_max_pk at which \f$P(k,z)\f$ and \f$ T_i(k,z)\f$ should be
      computable by interpolation. If it is equal to zero, only \f$ P(k,z=0)\f$
      needs to be computed. If it is higher, we will store in a table
      various P(k,tau) at several values of tau generously encompassing
      the range 0<z<z_max_pk */

  /* if z_max_pk<0, return error */
  class_test((psp->z_max_pk < 0),
             psp->error_message,
             "asked for negative redshift z=%e",psp->z_max_pk);

  /* if z_max_pk=0, there is just one value to store */
  if (psp->z_max_pk == 0.) {
    psp->ln_tau_size=1;
  }

  /* if z_max_pk>0, store several values (with a comfortable margin above z_max_pk) in view of interpolation */
  else{

    /* find the first relevant value of tau (last value in the table tau_ampling before tau(z_max)) and infer the number of values of tau at which P(k) must be stored */

    class_call(background_tau_of_z(pba,psp->z_max_pk,&tau_min),
               pba->error_message,
               psp->error_message);

    index_tau=0;
    class_test((tau_min <= ppt->tau_sampling[index_tau]),
               psp->error_message,
               "you asked for zmax=%e, i.e. taumin=%e, smaller than or equal to the first possible value =%e; it should be strictly bigger for a successfull interpolation",psp->z_max_pk,tau_min,ppt->tau_sampling[0]);

    while (ppt->tau_sampling[index_tau] < tau_min){
      index_tau++;
    }
    index_tau --;
    class_test(index_tau<0,
               psp->error_message,
               "by construction, this should never happen, a bug must have been introduced somewhere");

    /* whenever possible, take a few more values in to avoid boundary effects in the interpolation */
    if (index_tau>0) index_tau--;
    if (index_tau>0) index_tau--;
    if (index_tau>0) index_tau--;
    if (index_tau>0) index_tau--;
    psp->ln_tau_size=ppt->tau_size-index_tau;

  }

  /** - allocate and fill table of tau values at which \f$P(k,\tau)\f$ and \f$T_i(k,\tau)\f$ are stored */

  class_alloc(psp->ln_tau,sizeof(double)*psp->ln_tau_size,psp->error_message);

  for (index_tau=0; index_tau<psp->ln_tau_size; index_tau++) {
    psp->ln_tau[index_tau]=log(ppt->tau_sampling[index_tau-psp->ln_tau_size+ppt->tau_size]);
  }

  /** - allocate and fill table of k values at which \f$ P(k,\tau)\f$ is stored */

  psp->ln_k_size = ppt->k_size[ppt->index_md_scalars];
  class_alloc(psp->ln_k,sizeof(double)*psp->ln_k_size,psp->error_message);

  for (index_k=0; index_k<psp->ln_k_size; index_k++) {
    class_test(ppt->k[ppt->index_md_scalars][index_k] <= 0.,
               psp->error_message,
               "stop to avoid segmentation fault");
    psp->ln_k[index_k]=log(ppt->k[ppt->index_md_scalars][index_k]);
  }

  /** - if the non-linear power spectrum is requested, we should store
      it only at values of tau where non-linear corrections were
      really computed and not brutally set to one. Hence we must
      find here ln_tau_nl_size which might be smaller than
      ln_tau_size. But the same table ln_tau will be used for
      both. */

  if (pnl->method != nl_none) {

    index_tau=ppt->tau_size-psp->ln_tau_size;
    index_tau_min_nl=pnl->index_tau_min_nl;
    while (ppt->tau_sampling[index_tau] < pnl->tau[index_tau_min_nl]) {
      index_tau++;
    }
    psp->ln_tau_nl_size=ppt->tau_size-index_tau;

    class_alloc(psp->ln_tau_nl,sizeof(double)*psp->ln_tau_nl_size,psp->error_message);

    for (index_tau=0; index_tau<psp->ln_tau_nl_size; index_tau++) {
      psp->ln_tau_nl[index_tau]=log(ppt->tau_sampling[index_tau-psp->ln_tau_nl_size+ppt->tau_size]);
    }
  }

  return _SUCCESS_;
}

/**
 * This routine computes a table of values for all matter power spectra P(k),
 * given the source functions and primordial spectra.
 *
 * @param pba Input: pointer to background structure (will provide H, Omega_m at redshift of interest)
 * @param ppt Input: pointer to perturbation structure (contain source functions)
 * @param ppm Input: pointer to primordial structure
 * @param pnl Input: pointer to nonlinear structure
 * @param psp Input/Output: pointer to spectra structure
 * @return the error status
 */

int spectra_pk(
               struct background * pba,
               struct perturbs * ppt,
               struct primordial * ppm,
               struct nonlinear *pnl,
               struct spectra * psp
               ) {

  /** Summary: */

  /** - define local variables */

  int index_md;
  int index_ic1,index_ic2,index_ic1_ic1,index_ic2_ic2,index_ic1_ic2;
  int index_k;
  int index_tau;
  int delta_index_nl=0;
  int delta_index_nl_cb=0;
  double * primordial_pk; /* array with argument primordial_pk[index_ic_ic] */
  double source_ic1;
  double source_ic2;
  double pk_tot=0.,ln_pk_tot=0.;
  double source_ic1_cb;
  double source_ic2_cb;
  double pk_cb_tot=0.,ln_pk_cb_tot=0.;

  /** - check the presence of scalar modes */

  class_test((ppt->has_scalars == _FALSE_),
             psp->error_message,
             "you cannot ask for matter power spectrum since you turned off scalar modes");

  index_md = psp->index_md_scalars;

  /** - allocate temporary vectors where the primordial spectrum and the background quantities will be stored */

  class_alloc(primordial_pk,psp->ic_ic_size[index_md]*sizeof(double),psp->error_message);

  /** - allocate and fill array of \f$P(k,\tau)\f$ values */

  class_alloc(psp->ln_pk,
              sizeof(double)*psp->ln_tau_size*psp->ln_k_size*psp->ic_ic_size[index_md],
              psp->error_message);

  class_alloc(psp->ln_pk_l,
              sizeof(double)*psp->ln_tau_size*psp->ln_k_size,
              psp->error_message);

  if (pnl->method != nl_none) {

    class_alloc(psp->ln_pk_nl,
                sizeof(double)*psp->ln_tau_nl_size*psp->ln_k_size,
                psp->error_message);

    /* possible index shift between the first value of time used for
       the linear spectrum and that for the non-linear power
       spectrum (0 if no shift) */
    delta_index_nl = psp->ln_tau_size-psp->ln_tau_nl_size;
    class_test(delta_index_nl<0,
               "This should never happen",
               psp->error_message,
               psp->error_message);
  }
  else {
    psp->ln_pk_nl = NULL;
  }

  if (pba->has_ncdm) {

    class_alloc(psp->ln_pk_cb,
                sizeof(double)*psp->ln_tau_size*psp->ln_k_size*psp->ic_ic_size[index_md],
                psp->error_message);

    class_alloc(psp->ln_pk_cb_l,
                sizeof(double)*psp->ln_tau_size*psp->ln_k_size,
                psp->error_message);

    if (pnl->method != nl_none) {
      class_alloc(psp->ln_pk_cb_nl,
                  sizeof(double)*psp->ln_tau_nl_size*psp->ln_k_size,
                  psp->error_message);
      /* possible index shift between the first value of time used for
         the linear spectrum and that for the non-linear power
         spectrum (0 if no shift) */
      /* this is not really necessary, since m and cb share the same ln_tau_nl_size and ln_tau_nl */
      delta_index_nl_cb = psp->ln_tau_size-psp->ln_tau_nl_size;
      class_test(delta_index_nl_cb<0,
                 "This should never happen",
                 psp->error_message,
                 psp->error_message);
    }
    else{
      psp->ln_pk_cb_nl = NULL;
    }

  }
  else {
    psp->ln_pk_cb = NULL;
    psp->ln_pk_cb_l = NULL;
    psp->ln_pk_cb_nl = NULL;
  }

  for (index_tau=0 ; index_tau < psp->ln_tau_size; index_tau++) {
    for (index_k=0; index_k<psp->ln_k_size; index_k++) {

      class_call(primordial_spectrum_at_k(ppm,index_md,logarithmic,psp->ln_k[index_k],primordial_pk),
                 ppm->error_message,
                 psp->error_message);

      pk_tot =0;
      pk_cb_tot = 0.;
      /* curvature primordial spectrum:
         P_R(k) = 1/(2pi^2) k^3 <R R>
         so, primordial curvature correlator:
         <R R> = (2pi^2) k^-3 P_R(k)
         so, delta_m correlator:
         P(k) = <delta_m delta_m> = (2pi^2) k^-3 (source_m)^2 P_R(k)

         For isocurvature or cross adiabatic-isocurvature parts,
         replace one or two 'R' by 'S_i's */

      /* part diagonal in initial conditions */
      for (index_ic1 = 0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {

        index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_md]);

        source_ic1 = ppt->sources[index_md]
          [index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_m]
          [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

        psp->ln_pk[(index_tau * psp->ln_k_size + index_k)* psp->ic_ic_size[index_md] + index_ic1_ic2] =
          log(2.*_PI_*_PI_/exp(3.*psp->ln_k[index_k])
              *source_ic1*source_ic1
              *exp(primordial_pk[index_ic1_ic2]));

        pk_tot += exp(psp->ln_pk[(index_tau * psp->ln_k_size + index_k)* psp->ic_ic_size[index_md] + index_ic1_ic2]);

        if(pba->has_ncdm){

          source_ic1_cb = ppt->sources[index_md]
            [index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_cb]
            [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

          psp->ln_pk_cb[(index_tau * psp->ln_k_size + index_k)* psp->ic_ic_size[index_md] + index_ic1_ic2] =
            log(2.*_PI_*_PI_/exp(3.*psp->ln_k[index_k])
                *source_ic1_cb*source_ic1_cb
                *exp(primordial_pk[index_ic1_ic2]));

          pk_cb_tot += exp(psp->ln_pk_cb[(index_tau * psp->ln_k_size + index_k)* psp->ic_ic_size[index_md] + index_ic1_ic2]);

        }

      }

      /* part non-diagonal in initial conditions */
      for (index_ic1 = 0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
        for (index_ic2 = index_ic1+1; index_ic2 < psp->ic_size[index_md]; index_ic2++) {

          index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md]);
          index_ic1_ic1 = index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_md]);
          index_ic2_ic2 = index_symmetric_matrix(index_ic2,index_ic2,psp->ic_size[index_md]);

          if (psp->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {

            source_ic1 = ppt->sources[index_md]
              [index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_m]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            source_ic2 = ppt->sources[index_md]
              [index_ic2 * ppt->tp_size[index_md] + ppt->index_tp_delta_m]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->ln_pk[(index_tau * psp->ln_k_size + index_k)* psp->ic_ic_size[index_md] + index_ic1_ic2] =
              primordial_pk[index_ic1_ic2]*SIGN(source_ic1)*SIGN(source_ic2);

            pk_tot += psp->ln_pk[(index_tau * psp->ln_k_size + index_k)* psp->ic_ic_size[index_md] + index_ic1_ic2]
              * sqrt(psp->ln_pk[(index_tau * psp->ln_k_size + index_k)* psp->ic_ic_size[index_md] + index_ic1_ic1]
                     * psp->ln_pk[(index_tau * psp->ln_k_size + index_k)* psp->ic_ic_size[index_md] + index_ic2_ic2]);

            if(pba->has_ncdm){

              source_ic1_cb = ppt->sources[index_md]
                [index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_cb]
                [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

              source_ic2_cb = ppt->sources[index_md]
                [index_ic2 * ppt->tp_size[index_md] + ppt->index_tp_delta_cb]
                [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

              psp->ln_pk_cb[(index_tau * psp->ln_k_size + index_k)* psp->ic_ic_size[index_md] + index_ic1_ic2] =
                primordial_pk[index_ic1_ic2]*SIGN(source_ic1_cb)*SIGN(source_ic2_cb);

              pk_cb_tot += psp->ln_pk_cb[(index_tau * psp->ln_k_size + index_k)* psp->ic_ic_size[index_md] + index_ic1_ic2]
                * sqrt(psp->ln_pk_cb[(index_tau * psp->ln_k_size + index_k)* psp->ic_ic_size[index_md] + index_ic1_ic1]
                       * psp->ln_pk_cb[(index_tau * psp->ln_k_size + index_k)* psp->ic_ic_size[index_md] + index_ic2_ic2]);

            }

          }
          else {
            psp->ln_pk[(index_tau * psp->ln_k_size + index_k)* psp->ic_ic_size[index_md] + index_ic1_ic2] = 0.;
          }
        }
      }

      ln_pk_tot = log(pk_tot);

      if(pba->has_ncdm) ln_pk_cb_tot = log(pk_cb_tot);

      psp->ln_pk_l[index_tau * psp->ln_k_size + index_k] = ln_pk_tot;

      if(pba->has_ncdm) psp->ln_pk_cb_l[index_tau * psp->ln_k_size + index_k] = ln_pk_cb_tot;

      /* if non-linear corrections required, compute the total non-linear matter power spectrum */

      if ((pnl->method != nl_none) && (index_tau >= delta_index_nl)) {

        psp->ln_pk_nl[(index_tau-delta_index_nl) * psp->ln_k_size + index_k] =
          ln_pk_tot
          + 2.*log(pnl->nl_corr_density[pnl->index_pk_m][(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k]);

      }

      if((pba->has_ncdm) && (pnl->method != nl_none) && (index_tau >= delta_index_nl_cb)){

        psp->ln_pk_cb_nl[(index_tau-delta_index_nl_cb) * psp->ln_k_size + index_k] =
          ln_pk_cb_tot
          + 2.*log(pnl->nl_corr_density[pnl->index_pk_cb][(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k]);

      }

    }
  }

  /**- if interpolation of \f$P(k,\tau)\f$ will be needed (as a function of tau),
     compute array of second derivatives in view of spline interpolation */

  if (psp->ln_tau_size > 1) {

    class_alloc(psp->ddln_pk,sizeof(double)*psp->ln_tau_size*psp->ln_k_size*psp->ic_ic_size[index_md],psp->error_message);

    class_call(array_spline_table_lines(psp->ln_tau,
                                        psp->ln_tau_size,
                                        psp->ln_pk,
                                        psp->ic_ic_size[index_md]*psp->ln_k_size,
                                        psp->ddln_pk,
                                        _SPLINE_EST_DERIV_,
                                        psp->error_message),
               psp->error_message,
               psp->error_message);

    class_alloc(psp->ddln_pk_l,sizeof(double)*psp->ln_tau_size*psp->ln_k_size,psp->error_message);

    class_call(array_spline_table_lines(psp->ln_tau,
                                        psp->ln_tau_size,
                                        psp->ln_pk_l,
                                        psp->ln_k_size,
                                        psp->ddln_pk_l,
                                        _SPLINE_EST_DERIV_,
                                        psp->error_message),
               psp->error_message,
               psp->error_message);

    if(pba->has_ncdm){

      class_alloc(psp->ddln_pk_cb,sizeof(double)*psp->ln_tau_size*psp->ln_k_size*psp->ic_ic_size[index_md],psp->error_message);

      class_call(array_spline_table_lines(psp->ln_tau,
                                          psp->ln_tau_size,
                                          psp->ln_pk_cb,
                                          psp->ic_ic_size[index_md]*psp->ln_k_size,
                                          psp->ddln_pk_cb,
                                          _SPLINE_EST_DERIV_,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

      class_alloc(psp->ddln_pk_cb_l,sizeof(double)*psp->ln_tau_size*psp->ln_k_size,psp->error_message);

      class_call(array_spline_table_lines(psp->ln_tau,
                                          psp->ln_tau_size,
                                          psp->ln_pk_cb_l,
                                          psp->ln_k_size,
                                          psp->ddln_pk_cb_l,
                                          _SPLINE_EST_DERIV_,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);
    }

  }

  /* compute sigma8 (mean variance today in sphere of radius 8/h Mpc */

  class_call(spectra_sigma(pba,ppm,psp,8./pba->h,0.,&(psp->sigma8)),
             psp->error_message,
             psp->error_message);

  if (psp->spectra_verbose>0)
    fprintf(stdout," -> sigma8=%g (computed till k = %g h/Mpc)\n",
            psp->sigma8,
            exp(psp->ln_k[psp->ln_k_size-1])/pba->h);

  if(pba->has_ncdm){
    class_call(spectra_sigma_cb(pba,ppm,psp,8./pba->h,0.,&(psp->sigma8_cb)),
               psp->error_message,
               psp->error_message);
    if (psp->spectra_verbose>0){
      fprintf(stdout," -> sigma8 (ONLY CDM+BARYON)=%g (computed till k = %g h/Mpc)\n",
              psp->sigma8_cb,
              exp(psp->ln_k[psp->ln_k_size-1])/pba->h);
    }
  }

  /**- if interpolation of \f$ P_{NL}(k,\tau)\f$ will be needed (as a function of tau),
     compute array of second derivatives in view of spline interpolation */

  if (pnl->method != nl_none) {
    if (psp->ln_tau_nl_size > 1) {

      class_alloc(psp->ddln_pk_nl,sizeof(double)*psp->ln_tau_nl_size*psp->ln_k_size,psp->error_message);

      class_call(array_spline_table_lines(psp->ln_tau_nl,
                                          psp->ln_tau_nl_size,
                                          psp->ln_pk_nl,
                                          psp->ln_k_size,
                                          psp->ddln_pk_nl,
                                          _SPLINE_EST_DERIV_,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

      if(pba->has_ncdm){

        class_alloc(psp->ddln_pk_cb_nl,sizeof(double)*psp->ln_tau_nl_size*psp->ln_k_size,psp->error_message);

        class_call(array_spline_table_lines(psp->ln_tau_nl,
                                            psp->ln_tau_nl_size,
                                            psp->ln_pk_cb_nl,
                                            psp->ln_k_size,
                                            psp->ddln_pk_cb_nl,
                                            _SPLINE_EST_DERIV_,
                                            psp->error_message),
                   psp->error_message,
                   psp->error_message);

      }
    }
  }

  free (primordial_pk);

  return _SUCCESS_;
}

/**
 * This routine computes sigma(R) given P(k) (does not check that k_max is large
 * enough)
 *
 * @param pba   Input: pointer to background structure
 * @param ppm   Input: pointer to primordial structure
 * @param psp   Input: pointer to spectra structure
 * @param z     Input: redshift
 * @param R     Input: radius in Mpc
 * @param sigma Output: variance in a sphere of radius R (dimensionless)
 */

int spectra_sigma(
                  struct background * pba,
                  struct primordial * ppm,
                  struct spectra * psp,
                  double R,
                  double z,
                  double * sigma
                  ) {

  double pk;
  double * pk_ic = NULL;

  double pk_cb;
  double * pk_cb_ic = NULL;

  double * array_for_sigma;
  int index_num;
  int index_k;
  int index_y;
  int index_ddy;
  int i;

  double k,W,x;

  if (psp->ic_ic_size[psp->index_md_scalars]>1){
    class_alloc(pk_ic,
                psp->ic_ic_size[psp->index_md_scalars]*sizeof(double),
                psp->error_message);
    if (pba->has_ncdm)
      class_alloc(pk_cb_ic,
                  psp->ic_ic_size[psp->index_md_scalars]*sizeof(double),
                  psp->error_message);
  }

  i=0;
  index_k=i;
  i++;
  index_y=i;
  i++;
  index_ddy=i;
  i++;
  index_num=i;

  class_alloc(array_for_sigma,
              psp->ln_k_size*index_num*sizeof(double),
              psp->error_message);

  for (i=0;i<psp->ln_k_size;i++) {
    k=exp(psp->ln_k[i]);
    if (i == (psp->ln_k_size-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
    x=k*R;
    W=3./x/x/x*(sin(x)-x*cos(x));
    class_call(spectra_pk_at_k_and_z(pba,ppm,psp,k,z,&pk,pk_ic,&pk_cb,pk_cb_ic),
               psp->error_message,
               psp->error_message);
    array_for_sigma[i*index_num+index_k]=k;
    array_for_sigma[i*index_num+index_y]=k*k*pk*W*W;
  }

  class_call(array_spline(array_for_sigma,
                          index_num,
                          psp->ln_k_size,
                          index_k,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          psp->error_message),
             psp->error_message,
             psp->error_message);

  class_call(array_integrate_all_spline(array_for_sigma,
                                        index_num,
                                        psp->ln_k_size,
                                        index_k,
                                        index_y,
                                        index_ddy,
                                        sigma,
                                        psp->error_message),
             psp->error_message,
             psp->error_message);

  free(array_for_sigma);

  if (psp->ic_ic_size[psp->index_md_scalars]>1){
    free(pk_ic);
    if (pba->has_ncdm)
      free(pk_cb_ic);
  }

  *sigma = sqrt(*sigma/(2.*_PI_*_PI_));

  return _SUCCESS_;

}

//**/
int spectra_sigma_cb(
                     struct background * pba,
                     struct primordial * ppm,
                     struct spectra * psp,
                     double R,
                     double z,
                     double * sigma_cb
                     ) {

  double pk;
  double * pk_ic = NULL;

  double pk_cb;
  double * pk_cb_ic = NULL;

  double * array_for_sigma;
  int index_num;
  int index_k;
  int index_y;
  int index_ddy;
  int i;

  double k,W,x;

  if (psp->ic_ic_size[psp->index_md_scalars]>1){
    class_alloc(pk_ic,
                psp->ic_ic_size[psp->index_md_scalars]*sizeof(double),
                psp->error_message);
    if (pba->has_ncdm)
      class_alloc(pk_cb_ic,
                  psp->ic_ic_size[psp->index_md_scalars]*sizeof(double),
                  psp->error_message);
  }

  i=0;
  index_k=i;
  i++;
  index_y=i;
  i++;
  index_ddy=i;
  i++;
  index_num=i;

  class_alloc(array_for_sigma,
              psp->ln_k_size*index_num*sizeof(double),
              psp->error_message);

  for (i=0;i<psp->ln_k_size;i++) {
    k=exp(psp->ln_k[i]);
    if (i == (psp->ln_k_size-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
    x=k*R;
    W=3./x/x/x*(sin(x)-x*cos(x));
    class_call(spectra_pk_at_k_and_z(pba,ppm,psp,k,z,&pk,pk_ic,&pk_cb,pk_cb_ic),
               psp->error_message,
               psp->error_message);
    array_for_sigma[i*index_num+index_k]=k;
    array_for_sigma[i*index_num+index_y]=k*k*pk_cb*W*W;
  }

  class_call(array_spline(array_for_sigma,
                          index_num,
                          psp->ln_k_size,
                          index_k,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          psp->error_message),
             psp->error_message,
             psp->error_message);

  class_call(array_integrate_all_spline(array_for_sigma,
                                        index_num,
                                        psp->ln_k_size,
                                        index_k,
                                        index_y,
                                        index_ddy,
                                        sigma_cb,
                                        psp->error_message),
             psp->error_message,
             psp->error_message);

  free(array_for_sigma);

  if (psp->ic_ic_size[psp->index_md_scalars]>1){
    free(pk_ic);
    if (pba->has_ncdm)
      free(pk_cb_ic);
  }

  *sigma_cb = sqrt(*sigma_cb/(2.*_PI_*_PI_));

  return _SUCCESS_;

}

/**
 * This routine computes a table of values for all matter power spectra P(k),
 * given the source functions and primordial spectra.
 *
 * @param pba Input: pointer to background structure (will provide density of each species)
 * @param ppt Input: pointer to perturbation structure (contain source functions)
 * @param psp Input/Output: pointer to spectra structure
 * @return the error status
 */

int spectra_matter_transfers(
                             struct background * pba,
                             struct perturbs * ppt,
                             struct spectra * psp
                             ) {

  /** Summary: */

  /** - define local variables */

  int index_md;
  int index_ic;
  int index_k;
  int index_tau;
  int last_index_back;
  double * pvecback_sp_long; /* array with argument pvecback_sp_long[pba->index_bg] */
  double delta_i,theta_i,rho_i;
  double delta_rho_tot,rho_tot;
  double rho_plus_p_theta_tot,rho_plus_p_tot;
  int n_ncdm;
  double w_fld,dw_over_da_fld,integral_fld;

  /** - check the presence of scalar modes */

  class_test((ppt->has_scalars == _FALSE_),
             psp->error_message,
             "you cannot ask for matter power spectrum since you turned off scalar modes");

  index_md = psp->index_md_scalars;

  /** - allocate and fill array of \f$ T_i(k,\tau)\f$ values */

  class_alloc(psp->matter_transfer,sizeof(double)*psp->ln_tau_size*psp->ln_k_size*psp->ic_size[index_md]*psp->tr_size,psp->error_message);

  /** - allocate temporary vectors where the background quantities will be stored */

  class_alloc(pvecback_sp_long,pba->bg_size*sizeof(double),psp->error_message);

  for (index_tau=0 ; index_tau < psp->ln_tau_size; index_tau++) {

    class_call(background_at_tau(pba,
                                 ppt->tau_sampling[index_tau-psp->ln_tau_size+ppt->tau_size],
                                 /* for this last argument we could have passed
                                    exp(psp->ln_tau[index_tau]) but we would then loose
                                    precision in the exp(log(x)) operation) */
                                 pba->long_info,
                                 pba->inter_normal,
                                 &last_index_back,
                                 pvecback_sp_long),
               pba->error_message,
               psp->error_message);

    for (index_k=0; index_k<psp->ln_k_size; index_k++) {

      for (index_ic = 0; index_ic < psp->ic_size[index_md]; index_ic++) {

        delta_rho_tot=0.;
        rho_tot=0.;
        rho_plus_p_theta_tot=0.;
        rho_plus_p_tot=0.;

        /* T_g(k,tau) */

        rho_i = pvecback_sp_long[pba->index_bg_rho_g];

        if (ppt->has_source_delta_g == _TRUE_) {

          delta_i = ppt->sources[index_md]
            [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_g]
            [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

          psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_delta_g] = delta_i;

          delta_rho_tot += rho_i * delta_i;

          rho_tot += rho_i;

        }

        if (ppt->has_source_theta_g == _TRUE_) {

          theta_i = ppt->sources[index_md]
            [index_ic * ppt->tp_size[index_md] + ppt->index_tp_theta_g]
            [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

          psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_theta_g] = theta_i;

          rho_plus_p_theta_tot += 4./3. * rho_i * theta_i;

          rho_plus_p_tot += 4./3. * rho_i;

        }

        /* T_b(k,tau) */

        rho_i = pvecback_sp_long[pba->index_bg_rho_b];

        if (ppt->has_source_delta_b == _TRUE_) {

          delta_i = ppt->sources[index_md]
            [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_b]
            [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

          psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_delta_b] = delta_i;

          delta_rho_tot += rho_i * delta_i;

        }

        rho_tot += rho_i;

        if (ppt->has_source_theta_b == _TRUE_) {

          theta_i = ppt->sources[index_md]
            [index_ic * ppt->tp_size[index_md] + ppt->index_tp_theta_b]
            [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

          psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_theta_b] = theta_i;

          rho_plus_p_theta_tot += rho_i * theta_i;

        }

        rho_plus_p_tot += rho_i;

        /* T_cdm(k,tau) */

        if (pba->has_cdm == _TRUE_) {

          rho_i = pvecback_sp_long[pba->index_bg_rho_cdm];

          if (ppt->has_source_delta_cdm == _TRUE_) {

            delta_i = ppt->sources[index_md]
              [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_cdm]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_delta_cdm] = delta_i;

            delta_rho_tot += rho_i * delta_i;

          }

          rho_tot += rho_i;

          if (ppt->has_source_theta_cdm == _TRUE_) {

            theta_i = ppt->sources[index_md]
              [index_ic * ppt->tp_size[index_md] + ppt->index_tp_theta_cdm]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_theta_cdm] = theta_i;

            rho_plus_p_theta_tot += rho_i * theta_i;

          }

          rho_plus_p_tot += rho_i;

        }

        /* T_dcdm(k,tau) */

        if (pba->has_dcdm == _TRUE_) {

          rho_i = pvecback_sp_long[pba->index_bg_rho_dcdm];

          if (ppt->has_source_delta_dcdm == _TRUE_) {

            delta_i = ppt->sources[index_md]
              [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_dcdm]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_delta_dcdm] = delta_i;

            delta_rho_tot += rho_i * delta_i;

          }

          rho_tot += rho_i;

          if (ppt->has_source_theta_dcdm == _TRUE_) {

            theta_i = ppt->sources[index_md]
              [index_ic * ppt->tp_size[index_md] + ppt->index_tp_theta_dcdm]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_theta_dcdm] = theta_i;

            rho_plus_p_theta_tot += rho_i * theta_i;

          }

          rho_plus_p_tot += rho_i;

        }

        /* T_scf(k,tau) */

        if (pba->has_scf == _TRUE_) {

          rho_i = pvecback_sp_long[pba->index_bg_rho_scf];

          if (ppt->has_source_delta_scf == _TRUE_) {

            delta_i = ppt->sources[index_md]
              [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_scf]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_delta_scf] = delta_i;

            delta_rho_tot += rho_i * delta_i;

          }

          rho_tot += rho_i;

          if (ppt->has_source_theta_scf == _TRUE_) {

            theta_i = ppt->sources[index_md]
              [index_ic * ppt->tp_size[index_md] + ppt->index_tp_theta_scf]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_theta_scf] = theta_i;

            rho_plus_p_theta_tot += (rho_i + pvecback_sp_long[pba->index_bg_p_scf]) * theta_i;

          }

          rho_plus_p_tot += (rho_i + pvecback_sp_long[pba->index_bg_p_scf]);

        }


        /* T_fld(k,tau) */

        if (pba->has_fld == _TRUE_) {

          rho_i = pvecback_sp_long[pba->index_bg_rho_fld];

          class_call(background_w_fld(pba,0.,&w_fld,&dw_over_da_fld,&integral_fld), pba->error_message, psp->error_message);

          if (ppt->has_source_delta_fld == _TRUE_) {

            delta_i = ppt->sources[index_md]
              [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_fld]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_delta_fld] = delta_i;

            delta_rho_tot += rho_i * delta_i;

          }

          rho_tot += rho_i;

          if (ppt->has_source_theta_fld == _TRUE_) {

            theta_i = ppt->sources[index_md]
              [index_ic * ppt->tp_size[index_md] + ppt->index_tp_theta_fld]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_theta_fld] = theta_i;

            rho_plus_p_theta_tot += (1. + w_fld) * rho_i * theta_i;

          }

          rho_plus_p_tot += (1. + w_fld) * rho_i;

        }

        /* T_ur(k,tau) */

        if (pba->has_ur == _TRUE_) {

          rho_i = pvecback_sp_long[pba->index_bg_rho_ur];

          if (ppt->has_source_delta_ur == _TRUE_) {

            delta_i = ppt->sources[index_md]
              [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_ur]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_delta_ur] = delta_i;

            delta_rho_tot += rho_i * delta_i;

          }

          rho_tot += rho_i;

          if (ppt->has_source_theta_ur == _TRUE_) {

            theta_i = ppt->sources[index_md]
              [index_ic * ppt->tp_size[index_md] + ppt->index_tp_theta_ur]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_theta_ur] = theta_i;

            rho_plus_p_theta_tot += 4./3. * rho_i * theta_i;

          }

          rho_plus_p_tot += 4./3. * rho_i;

        }

        /* T_dr(k,tau) */

        if (pba->has_dr == _TRUE_) {

          rho_i = pvecback_sp_long[pba->index_bg_rho_dr];

          if (ppt->has_source_delta_dr == _TRUE_) {

            delta_i = ppt->sources[index_md]
              [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_dr]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_delta_dr] = delta_i;

            delta_rho_tot += rho_i * delta_i;

          }

          rho_tot += rho_i;

          if (ppt->has_source_theta_dr == _TRUE_) {

            theta_i = ppt->sources[index_md]
              [index_ic * ppt->tp_size[index_md] + ppt->index_tp_theta_dr]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_theta_dr] = theta_i;

            rho_plus_p_theta_tot += 4./3. * rho_i * theta_i;

          }

          rho_plus_p_tot += 4./3. * rho_i;

        }

        /* T_ncdm_i(k,tau) */

        if (pba->has_ncdm == _TRUE_) {

          for (n_ncdm=0; n_ncdm < pba->N_ncdm; n_ncdm++) {

            rho_i = pvecback_sp_long[pba->index_bg_rho_ncdm1+n_ncdm];

            if (ppt->has_source_delta_ncdm == _TRUE_) {

              delta_i = ppt->sources[index_md]
                [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_ncdm1+n_ncdm]
                [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

              psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_delta_ncdm1+n_ncdm] = delta_i;

              delta_rho_tot += rho_i * delta_i;

            }

            rho_tot += rho_i;

            if (ppt->has_source_theta_ncdm == _TRUE_) {

              theta_i = ppt->sources[index_md]
                [index_ic * ppt->tp_size[index_md] + ppt->index_tp_theta_ncdm1+n_ncdm]
                [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

              psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_theta_ncdm1+n_ncdm] = theta_i;

              rho_plus_p_theta_tot += (rho_i + pvecback_sp_long[pba->index_bg_p_ncdm1+n_ncdm]) * theta_i;

            }

            rho_plus_p_tot += (rho_i + pvecback_sp_long[pba->index_bg_p_ncdm1+n_ncdm]);

          }

        }

        if (ppt->has_source_phi == _TRUE_) {

          psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_phi] = ppt->sources[index_md]
            [index_ic * ppt->tp_size[index_md] + ppt->index_tp_phi]
            [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

        }

        if (ppt->has_source_psi == _TRUE_) {

          psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_psi] = ppt->sources[index_md]
            [index_ic * ppt->tp_size[index_md] + ppt->index_tp_psi]
            [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

        }

        if (ppt->has_source_phi_prime == _TRUE_) {

          psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_phi_prime] = ppt->sources[index_md]
            [index_ic * ppt->tp_size[index_md] + ppt->index_tp_phi_prime]
            [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

        }

        if (ppt->has_source_h == _TRUE_) {

          psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_h] = ppt->sources[index_md]
            [index_ic * ppt->tp_size[index_md] + ppt->index_tp_h]
            [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

        }

        if (ppt->has_source_h_prime == _TRUE_) {

          psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_h_prime] = ppt->sources[index_md]
            [index_ic * ppt->tp_size[index_md] + ppt->index_tp_h_prime]
            [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

        }

        if (ppt->has_source_eta == _TRUE_) {

          psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_eta] = ppt->sources[index_md]
            [index_ic * ppt->tp_size[index_md] + ppt->index_tp_eta]
            [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

        }

        if (ppt->has_source_eta_prime == _TRUE_) {

          psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_eta_prime] = ppt->sources[index_md]
            [index_ic * ppt->tp_size[index_md] + ppt->index_tp_eta_prime]
            [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

        }

        /* could include homogeneous component in rho_tot if uncommented (leave commented to match CMBFAST/CAMB definition) */

        /* 	if (pba->has_lambda == _TRUE_) { */

        /* 	  rho_i = pvecback_sp_long[pba->index_bg_rho_lambda]; */

        /* 	  rho_tot += rho_i; */
        /* 	} */

        /* T_tot(k,tau) */

        if (ppt->has_density_transfers == _TRUE_) {

          psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_delta_tot] = delta_rho_tot/rho_tot;

        }

        if (ppt->has_velocity_transfers == _TRUE_) {

          psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_theta_tot] = rho_plus_p_theta_tot/rho_plus_p_tot;

        }

      }
    }
  }

  /**- if interpolation of \f$ P(k,\tau)\f$ will be needed (as a function of tau),
     compute array of second derivatives in view of spline interpolation */

  if (psp->ln_tau_size > 1) {

    class_alloc(psp->ddmatter_transfer,sizeof(double)*psp->ln_tau_size*psp->ln_k_size*psp->ic_size[index_md]*psp->tr_size,psp->error_message);

    class_call(array_spline_table_lines(psp->ln_tau,
                                        psp->ln_tau_size,
                                        psp->matter_transfer,
                                        psp->ic_size[index_md]*psp->ln_k_size*psp->tr_size,
                                        psp->ddmatter_transfer,
                                        _SPLINE_EST_DERIV_,
                                        psp->error_message),
               psp->error_message,
               psp->error_message);

  }

  free (pvecback_sp_long);

  return _SUCCESS_;
}

int spectra_output_tk_titles(struct background *pba,
                             struct perturbs *ppt,
                             enum file_format output_format,
                             char titles[_MAXTITLESTRINGLENGTH_]
                             ){
  int n_ncdm;
  char tmp[40];

  if (output_format == class_format) {
    class_store_columntitle(titles,"k (h/Mpc)",_TRUE_);
    if (ppt->has_density_transfers == _TRUE_) {
      class_store_columntitle(titles,"d_g",_TRUE_);
      class_store_columntitle(titles,"d_b",_TRUE_);
      class_store_columntitle(titles,"d_cdm",pba->has_cdm);
      class_store_columntitle(titles,"d_fld",pba->has_fld);
      class_store_columntitle(titles,"d_ur",pba->has_ur);
      if (pba->has_ncdm == _TRUE_) {
        for (n_ncdm=0; n_ncdm < pba->N_ncdm; n_ncdm++) {
          sprintf(tmp,"d_ncdm[%d]",n_ncdm);
          class_store_columntitle(titles,tmp,_TRUE_);
        }
      }
      class_store_columntitle(titles,"d_dcdm",pba->has_dcdm);
      class_store_columntitle(titles,"d_dr",pba->has_dr);
      class_store_columntitle(titles,"d_scf",pba->has_scf);
      class_store_columntitle(titles,"d_tot",_TRUE_);
      class_store_columntitle(titles,"phi",ppt->has_source_phi);
      class_store_columntitle(titles,"psi",ppt->has_source_psi);
      class_store_columntitle(titles,"phi_prime",ppt->has_source_phi_prime);
      class_store_columntitle(titles,"h",ppt->has_source_h);
      class_store_columntitle(titles,"h_prime",ppt->has_source_h_prime);
      class_store_columntitle(titles,"eta",ppt->has_source_eta);
      class_store_columntitle(titles,"eta_prime",ppt->has_source_eta_prime);
    }
    if (ppt->has_velocity_transfers == _TRUE_) {
      class_store_columntitle(titles,"t_g",_TRUE_);
      class_store_columntitle(titles,"t_b",_TRUE_);
      class_store_columntitle(titles,"t_cdm",((pba->has_cdm == _TRUE_) && (ppt->gauge != synchronous)));
      class_store_columntitle(titles,"t_fld",pba->has_fld);
      class_store_columntitle(titles,"t_ur",pba->has_ur);
      if (pba->has_ncdm == _TRUE_) {
        for (n_ncdm=0; n_ncdm < pba->N_ncdm; n_ncdm++) {
          sprintf(tmp,"t_ncdm[%d]",n_ncdm);
          class_store_columntitle(titles,tmp,_TRUE_);
        }
      }
      class_store_columntitle(titles,"t_dcdm",pba->has_dcdm);
      class_store_columntitle(titles,"t_dr",pba->has_dr);
      class_store_columntitle(titles,"t__scf",pba->has_scf);
      class_store_columntitle(titles,"t_tot",_TRUE_);
    }
  }

  else if (output_format == camb_format) {

    class_store_columntitle(titles,"k (h/Mpc)",_TRUE_);
    class_store_columntitle(titles,"-T_cdm/k2",_TRUE_);
    class_store_columntitle(titles,"-T_b/k2",_TRUE_);
    class_store_columntitle(titles,"-T_g/k2",_TRUE_);
    class_store_columntitle(titles,"-T_ur/k2",_TRUE_);
    class_store_columntitle(titles,"-T_ncdm/k2",_TRUE_);
    class_store_columntitle(titles,"-T_tot/k2",_TRUE_);

  }

  return _SUCCESS_;

}

int spectra_output_tk_data(
                           struct background * pba,
                           struct perturbs * ppt,
                           struct spectra * psp,
                           enum file_format output_format,
                           double z,
                           int number_of_titles,
                           double *data
                           ) {

  int n_ncdm;
  double k, k_over_h, k2;
  double * tkfull=NULL;  /* array with argument
                            pk_ic[(index_k * psp->ic_size[index_md] + index_ic)*psp->tr_size+index_tr] */
  double *tk;
  double *dataptr;

  int index_md=0;
  int index_ic;
  int index_k;
  int index_tr;
  int storeidx;

  if (psp->ln_k_size*psp->ic_size[index_md]*psp->tr_size > 0){
    class_alloc(tkfull,
                psp->ln_k_size*psp->ic_size[index_md]*psp->tr_size*sizeof(double),
                psp->error_message);
  }

  /** - compute \f$T_i(k)\f$ for each k (if several ic's, compute it for each ic; if z_pk = 0, this is done by directly reading inside the pre-computed table; if not, this is done by interpolating the table at the correct value of tau. */

  /* if z_pk = 0, no interpolation needed */

  if (z == 0.) {

    for (index_k=0; index_k<psp->ln_k_size; index_k++) {
      for (index_tr=0; index_tr<psp->tr_size; index_tr++) {
        for (index_ic=0; index_ic<psp->ic_size[index_md]; index_ic++) {
          tkfull[(index_k * psp->ic_size[index_md] + index_ic) * psp->tr_size + index_tr] = psp->matter_transfer[(((psp->ln_tau_size-1)*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + index_tr];
        }
      }
    }
  }

  /* if 0 <= z_pk <= z_max_pk, interpolation needed, */
  else {

    class_call(spectra_tk_at_z(pba,
                               psp,
                               z,
                               tkfull),
               psp->error_message,
               psp->error_message);
  }

  /** - store data */

  for (index_ic = 0; index_ic < psp->ic_size[index_md]; index_ic++) {

    for (index_k=0; index_k<psp->ln_k_size; index_k++) {

      storeidx = 0;
      dataptr = data+index_ic*(psp->ln_k_size*number_of_titles)+index_k*number_of_titles;
      tk = &(tkfull[(index_k * psp->ic_size[index_md] + index_ic) * psp->tr_size]);
      k = exp(psp->ln_k[index_k]);
      k2 = k*k;
      k_over_h = k/pba->h;

      class_store_double(dataptr, k_over_h, _TRUE_,storeidx);

      /* indices for species associated with a velocity transfer function in Fourier space */

      if (output_format == class_format) {

        if (ppt->has_density_transfers == _TRUE_) {

          class_store_double(dataptr,tk[psp->index_tr_delta_g],ppt->has_source_delta_g,storeidx);
          class_store_double(dataptr,tk[psp->index_tr_delta_b],ppt->has_source_delta_b,storeidx);
          class_store_double(dataptr,tk[psp->index_tr_delta_cdm],ppt->has_source_delta_cdm,storeidx);
          class_store_double(dataptr,tk[psp->index_tr_delta_fld],ppt->has_source_delta_fld,storeidx);
          class_store_double(dataptr,tk[psp->index_tr_delta_ur],ppt->has_source_delta_ur,storeidx);
          if (pba->has_ncdm == _TRUE_){
            for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++){
              class_store_double(dataptr,tk[psp->index_tr_delta_ncdm1+n_ncdm],ppt->has_source_delta_ncdm,storeidx);
            }
          }
          class_store_double(dataptr,tk[psp->index_tr_delta_dcdm],ppt->has_source_delta_dcdm,storeidx);
          class_store_double(dataptr,tk[psp->index_tr_delta_dr],ppt->has_source_delta_dr,storeidx);
          class_store_double(dataptr,tk[psp->index_tr_delta_scf],ppt->has_source_delta_scf,storeidx);
          class_store_double(dataptr,tk[psp->index_tr_delta_tot],_TRUE_,storeidx);
          class_store_double(dataptr,tk[psp->index_tr_phi],ppt->has_source_phi,storeidx);
          class_store_double(dataptr,tk[psp->index_tr_psi],ppt->has_source_psi,storeidx);
          class_store_double(dataptr,tk[psp->index_tr_phi_prime],ppt->has_source_phi_prime,storeidx);
          class_store_double(dataptr,tk[psp->index_tr_h],ppt->has_source_h,storeidx);
          class_store_double(dataptr,tk[psp->index_tr_h_prime],ppt->has_source_h_prime,storeidx);
          class_store_double(dataptr,tk[psp->index_tr_eta],ppt->has_source_eta,storeidx);
          class_store_double(dataptr,tk[psp->index_tr_eta_prime],ppt->has_source_eta_prime,storeidx);
        }
        if (ppt->has_velocity_transfers == _TRUE_) {

          class_store_double(dataptr,tk[psp->index_tr_theta_g],ppt->has_source_theta_g,storeidx);
          class_store_double(dataptr,tk[psp->index_tr_theta_b],ppt->has_source_theta_b,storeidx);
          class_store_double(dataptr,tk[psp->index_tr_theta_cdm],ppt->has_source_theta_cdm,storeidx);
          class_store_double(dataptr,tk[psp->index_tr_theta_fld],ppt->has_source_theta_fld,storeidx);
          class_store_double(dataptr,tk[psp->index_tr_theta_ur],ppt->has_source_theta_ur,storeidx);
          if (pba->has_ncdm == _TRUE_){
            for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++){
              class_store_double(dataptr,tk[psp->index_tr_theta_ncdm1+n_ncdm],ppt->has_source_theta_ncdm,storeidx);
            }
          }
          class_store_double(dataptr,tk[psp->index_tr_theta_dcdm],ppt->has_source_theta_dcdm,storeidx);
          class_store_double(dataptr,tk[psp->index_tr_theta_dr],ppt->has_source_theta_dr,storeidx);
          class_store_double(dataptr,tk[psp->index_tr_theta_scf],ppt->has_source_theta_scf,storeidx);
          class_store_double(dataptr,tk[psp->index_tr_theta_tot],_TRUE_,storeidx);

        }

      }
      else if (output_format == camb_format) {

        /* rescale and reorder the matter transfer functions following the CMBFAST/CAMB convention */
        class_store_double_or_default(dataptr,-tk[psp->index_tr_delta_cdm]/k2,ppt->has_source_delta_cdm,storeidx,0.0);
        class_store_double_or_default(dataptr,-tk[psp->index_tr_delta_b]/k2,ppt->has_source_delta_b,storeidx,0.0);
        class_store_double_or_default(dataptr,-tk[psp->index_tr_delta_g]/k2,ppt->has_source_delta_g,storeidx,0.0);
        class_store_double_or_default(dataptr,-tk[psp->index_tr_delta_ur]/k2,ppt->has_source_delta_ur,storeidx,0.0);
        class_store_double_or_default(dataptr,-tk[psp->index_tr_delta_ncdm1]/k2,ppt->has_source_delta_ncdm,storeidx,0.0);
        class_store_double_or_default(dataptr,-tk[psp->index_tr_delta_tot]/k2,_TRUE_,storeidx,0.0);
      }
    }
  }

  //Necessary because the size could be zero (if psp->tr_size is zero)
  if (tkfull != NULL)
    free(tkfull);

  return _SUCCESS_;
}

int spectra_firstline_and_ic_suffix(struct perturbs *ppt,
                                    int index_ic,
                                    char first_line[_LINE_LENGTH_MAX_],
                                    FileName ic_suffix){

  first_line[0]='\0';
  ic_suffix[0]='\0';


  if ((ppt->has_ad == _TRUE_) && (index_ic == ppt->index_ic_ad)) {
    strcpy(ic_suffix,"ad");
    strcpy(first_line,"for adiabatic (AD) mode (normalized to initial curvature=1) ");
  }

  if ((ppt->has_bi == _TRUE_) && (index_ic == ppt->index_ic_bi)) {
    strcpy(ic_suffix,"bi");
    strcpy(first_line,"for baryon isocurvature (BI) mode (normalized to initial entropy=1)");
  }

  if ((ppt->has_cdi == _TRUE_) && (index_ic == ppt->index_ic_cdi)) {
    strcpy(ic_suffix,"cdi");
    strcpy(first_line,"for CDM isocurvature (CDI) mode (normalized to initial entropy=1)");
  }

  if ((ppt->has_nid == _TRUE_) && (index_ic == ppt->index_ic_nid)) {
    strcpy(ic_suffix,"nid");
    strcpy(first_line,"for neutrino density isocurvature (NID) mode (normalized to initial entropy=1)");
  }

  if ((ppt->has_niv == _TRUE_) && (index_ic == ppt->index_ic_niv)) {
    strcpy(ic_suffix,"niv");
    strcpy(first_line,"for neutrino velocity isocurvature (NIV) mode (normalized to initial entropy=1)");
  }
  return _SUCCESS_;
}

int spectra_fast_pk_at_kvec_and_zvec(
                                     struct background * pba,
                                     struct spectra * psp,
                                     double * kvec,
                                     int kvec_size,
                                     double * zvec,
                                     int zvec_size,
                                     double * pk_tot_out, /* (must be already allocated with kvec_size*zvec_size) */
                                     double * pk_cb_tot_out,
                                     int nonlinear
                                     ) {

  /** Summary: */

  /** - define local variables */

  int index_md;
  int index_k, index_knode, index_z;
  double *spline, *pk_at_k, *ln_kvec, *ln_pk_table;
  double h, a, b;
  double *spline_cb, *pk_cb_at_k, *ln_pk_cb_table;

  index_md = psp->index_md_scalars;
  class_test(psp->ic_size[index_md] != 1,
             psp->error_message,
             "This function has only been coded for pure adiabatic ICs, sorry.");

  /** Compute spline over ln(k) */
  class_alloc(ln_kvec, sizeof(double)*kvec_size,psp->error_message);
  class_alloc(ln_pk_table, sizeof(double)*psp->ln_k_size*zvec_size,psp->error_message);
  class_alloc(spline, sizeof(double)*psp->ln_k_size*zvec_size,psp->error_message);
  class_alloc(pk_at_k, sizeof(double)*psp->ln_tau_size,psp->error_message);

  class_alloc(ln_pk_cb_table, sizeof(double)*psp->ln_k_size*zvec_size,psp->error_message);
  class_alloc(spline_cb, sizeof(double)*psp->ln_k_size*zvec_size,psp->error_message);
  class_alloc(pk_cb_at_k, sizeof(double)*psp->ln_tau_size,psp->error_message);

  /** Construct table of log(pk) on the computed k nodes but requested redshifts: */
  for (index_z=0; index_z<zvec_size; index_z++){
    if (nonlinear==_TRUE_) {
      class_call(spectra_pk_nl_at_z(pba,psp, logarithmic,zvec[index_z],ln_pk_table+index_z*psp->ln_k_size,ln_pk_cb_table+index_z*psp->ln_k_size),
                 psp->error_message,
                 psp->error_message);
    }
    else{
      class_call(spectra_pk_at_z(pba,psp, logarithmic,zvec[index_z],ln_pk_table+index_z*psp->ln_k_size, NULL,ln_pk_cb_table+index_z*psp->ln_k_size,NULL),
                 psp->error_message,
                 psp->error_message);
    }
  }

  class_call(array_spline_table_columns2(psp->ln_k,
                                         psp->ln_k_size,
                                         ln_pk_table,
                                         zvec_size,
                                         spline,
                                         _SPLINE_NATURAL_,
                                         psp->error_message),
             psp->error_message,
             psp->error_message);

  if(pba->has_ncdm){
    class_call(array_spline_table_columns2(psp->ln_k,
                                           psp->ln_k_size,
                                           ln_pk_cb_table,
                                           zvec_size,
                                           spline_cb,
                                           _SPLINE_NATURAL_,
                                           psp->error_message),
               psp->error_message,
               psp->error_message);
  }
  /** Construct ln(kvec): */
  for (index_k=0; index_k<kvec_size; index_k++){
    ln_kvec[index_k] = log(kvec[index_k]);
  }

  /** I will assume that the k vector is sorted in ascending order.
      Case k<kmin: */
  for(index_k = 0; index_k<kvec_size; index_k++){
    if (ln_kvec[index_k] >= psp->ln_k[0])
      break;
    for (index_z = 0; index_z < zvec_size; index_z++) {
      /** If needed, add some extrapolation here */
      pk_tot_out[index_z*kvec_size+index_k] = 0.;
      if(pba->has_ncdm) pk_cb_tot_out[index_z*kvec_size+index_k] = 0.;
    }
    /** Implement some extrapolation perhaps */
  }

  /** Case kmin<=k<=kmax. Do not loop through kvec, but loop through the
      interpolation nodes. */
  for (index_knode=0; index_knode < (psp->ln_k_size-1); index_knode++){
    /** Loop through k's that fall in this interval */
    //printf("index _k is %d, do we have %g < %g <%g?\n",index_k, psp->ln_k[index_knode],ln_kvec[index_k],psp->ln_k[index_knode+1]);
    while ((index_k < kvec_size) && (ln_kvec[index_k] <= psp->ln_k[index_knode+1])){
      /** Perform interpolation */
      h = psp->ln_k[index_knode+1]-psp->ln_k[index_knode];
      b = (ln_kvec[index_k] - psp->ln_k[index_knode])/h;
      a = 1.-b;
      for (index_z = 0; index_z < zvec_size; index_z++) {
        pk_tot_out[index_z*kvec_size+index_k] =
          exp(
              a * ln_pk_table[index_z*psp->ln_k_size + index_knode]
              + b * ln_pk_table[index_z*psp->ln_k_size + index_knode+1]
              + ((a*a*a-a) * spline[index_z*psp->ln_k_size + index_knode]
                 +(b*b*b-b) * spline[index_z*psp->ln_k_size + index_knode+1])*h*h/6.0
              );
        if(pba->has_ncdm){
          pk_cb_tot_out[index_z*kvec_size+index_k] =
            exp(
                a * ln_pk_cb_table[index_z*psp->ln_k_size + index_knode]
                + b * ln_pk_cb_table[index_z*psp->ln_k_size + index_knode+1]
                + ((a*a*a-a) * spline_cb[index_z*psp->ln_k_size + index_knode]
                   +(b*b*b-b) * spline_cb[index_z*psp->ln_k_size + index_knode+1])*h*h/6.0
                );
        }
      }
      index_k++;
    }
  }

  /** case k>kmax */
  while (index_k<kvec_size){
    for (index_z = 0; index_z < zvec_size; index_z++) {
      /** If needed, add some extrapolation here */
      pk_tot_out[index_z*kvec_size+index_k] = 0.;
      if(pba->has_ncdm) pk_cb_tot_out[index_z*kvec_size+index_k] = 0.;
    }
    index_k++;
  }

  free(ln_kvec);
  free(ln_pk_table);
  free(ln_pk_cb_table);
  free(spline);
  if(pba->has_ncdm) free(spline_cb);
  free(pk_at_k);
  return _SUCCESS_;
}
