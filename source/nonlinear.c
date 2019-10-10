/** @file nonlinear.c Documented nonlinear module
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

#include "nonlinear.h"

/**
 * Return the value of the non-linearity wavenumber k_nl for a given redshift z
 * @param pba Input: pointer to background structure
 * @param pnl Input: pointer to nonlinear structure
 * @param z Input: redshift
 * @param k_nl Output: k_nl value
 * @param k_nl_cb Ouput: k_nl value of the cdm+baryon part only, if there is ncdm
 * @return the error status
 */

int nonlinear_k_nl_at_z(
                        struct background *pba,
                        struct nonlinear * pnl,
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
             pnl->error_message);

  /** - interpolate the precomputed k_nl array at the needed valuetime */

  if (pnl->has_pk_m == _TRUE_) {

    if (pnl->tau_size == 1) {
      *k_nl = pnl->k_nl[pnl->index_pk_m][0];
    }
    else {
      class_call(array_interpolate_two(pnl->tau,
                                       1,
                                       0,
                                       pnl->k_nl[pnl->index_pk_m],
                                       1,
                                       pnl->tau_size,
                                       tau,
                                       k_nl,
                                       1,
                                       pnl->error_message),
                 pnl->error_message,
                 pnl->error_message);
    }
  }

  /** - if needed, do the same for the baryon part only */

  if (pnl->has_pk_cb) {

    if (pnl->tau_size == 1) {
      *k_nl_cb = pnl->k_nl[pnl->index_pk_cb][0];
    }
    else {
      class_call(array_interpolate_two(pnl->tau,
                                       1,
                                       0,
                                       pnl->k_nl[pnl->index_pk_cb],
                                       1,
                                       pnl->tau_size,
                                       tau,
                                       k_nl_cb,
                                       1,
                                       pnl->error_message),
                 pnl->error_message,
                 pnl->error_message);
    }
  }

  /* otherwise, return the same for k_nl_cb as for k_nl */

  else{
    *k_nl_cb = *k_nl;
  }

  return _SUCCESS_;
}

/**
 * Initialize the nonlinear structure, and in particular the
 * nl_corr_density and k_nl interpolation tables.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param pth Input: pointer to therodynamics structure
 * @param ppt Input: pointer to perturbation structure
 * @param ppm Input: pointer to primordial structure
 * @param pnl Input/Output: pointer to initialized nonlinear structure
 * @return the error status
 */

int nonlinear_init(
                   struct precision *ppr,
                   struct background *pba,
                   struct thermo *pth,
                   struct perturbs *ppt,
                   struct primordial *ppm,
                   struct nonlinear *pnl
                   ) {

  int index_ncdm;
  int index_k;
  int index_tau;
  int index_pk;
  double *pk_l;
  double *pk_nl;
  double *lnk_l;
  double *lnpk_l;
  double *ddlnpk_l;
  short print_warning=_FALSE_;
  double * pvecback;
  int last_index;
  double a,z;
  short halofit_found_k_max;

  /** Define flags and indices (so few that no dedicated routine needed) */

  pnl->has_pk_m = _TRUE_;
  if (pba->has_ncdm == _TRUE_) {
    pnl->has_pk_cb = _TRUE_;
  }
  else {
    pnl->has_pk_cb = _FALSE_;
  }

  index_pk = 0;
  class_define_index(pnl->index_pk_m, pnl->has_pk_m, index_pk,1);
  class_define_index(pnl->index_pk_cb, pnl->has_pk_cb, index_pk,1);
  pnl->pk_size = index_pk;

  /** (a) First deal with the case where non non-linear corrections requested */

  if (pnl->method == nl_none) {
    if (pnl->nonlinear_verbose > 0)
      printf("No non-linear spectra requested. Nonlinear module skipped.\n");
  }

  /** (b) Compute for HALOFIT non-linear spectrum */

  else if (pnl->method == nl_halofit) {
    if (pnl->nonlinear_verbose > 0)
      printf("Computing non-linear matter power spectrum with Halofit (including update Takahashi et al. 2012 and Bird 2014)\n");

    if (pba->has_ncdm) {
      for (index_ncdm=0;index_ncdm < pba->N_ncdm; index_ncdm++){
        if (pba->m_ncdm_in_eV[index_ncdm] >  _M_EV_TOO_BIG_FOR_HALOFIT_)
          fprintf(stdout,"Warning: Halofit is proved to work for CDM, and also with a small HDM component thanks to Bird et al.'s update. But it sounds like you are running with a WDM component of mass %f eV, which makes the use of Halofit suspicious.\n",pba->m_ncdm_in_eV[index_ncdm]);
      }
    }

    /** - copy list of (k,tau) from perturbation module */

    pnl->k_size = ppt->k_size[ppt->index_md_scalars];
    class_alloc(pnl->k,pnl->k_size*sizeof(double),pnl->error_message);
    for (index_k=0; index_k<pnl->k_size; index_k++)
      pnl->k[index_k] = ppt->k[ppt->index_md_scalars][index_k];

    pnl->tau_size = ppt->tau_size;
    class_alloc(pnl->tau,pnl->tau_size*sizeof(double),pnl->error_message);
    for (index_tau=0; index_tau<pnl->tau_size; index_tau++)
      pnl->tau[index_tau] = ppt->tau_sampling[index_tau];

    class_alloc(pnl->nl_corr_density,
                pnl->pk_size*sizeof(double *),
                pnl->error_message);

    class_alloc(pnl->k_nl,
                pnl->pk_size*sizeof(double *),
                pnl->error_message);

    class_alloc(pk_l,
                pnl->k_size*sizeof(double),
                pnl->error_message);

    class_alloc(pk_nl,
                pnl->k_size*sizeof(double),
                pnl->error_message);

    class_alloc(lnk_l,
                pnl->k_size*sizeof(double),
                pnl->error_message);

    class_alloc(lnpk_l,
                pnl->k_size*sizeof(double),
                pnl->error_message);

    class_alloc(ddlnpk_l,
                pnl->k_size*sizeof(double),
                pnl->error_message);

    for (index_pk=0; index_pk<pnl->pk_size; index_pk++){
      class_alloc(pnl->nl_corr_density[index_pk],pnl->tau_size*pnl->k_size*sizeof(double),pnl->error_message);
      class_alloc(pnl->k_nl[index_pk],pnl->tau_size*sizeof(double),pnl->error_message);
    }

    print_warning=_FALSE_;

    pnl->index_tau_min_nl = 0;

    /** - loop over time */

    for (index_tau = pnl->tau_size-1; index_tau>=0; index_tau--) {

      for (index_pk=0; index_pk<pnl->pk_size; index_pk++) {

        /* get P_L(k) at this time */
        class_call(nonlinear_pk_l(pba,ppt,ppm,pnl,index_pk,index_tau,pk_l,lnk_l,lnpk_l,ddlnpk_l),
                   pnl->error_message,
                   pnl->error_message);

        /* get P_NL(k) at this time */
        if (print_warning == _FALSE_) {

          class_call(nonlinear_halofit(ppr,
                                       pba,
                                       ppt,
                                       ppm,
                                       pnl,
                                       index_pk,
                                       pnl->tau[index_tau],
                                       pk_l,
                                       pk_nl,
                                       lnk_l,
                                       lnpk_l,
                                       ddlnpk_l,
                                       &(pnl->k_nl[index_pk][index_tau]),
                                       &halofit_found_k_max),
                     pnl->error_message,
                     pnl->error_message);

          if (halofit_found_k_max == _TRUE_) {

            // for debugging:
            /*if ((index_tau == pnl->tau_size-1)){
              for (index_k=0; index_k<pnl->k_size; index_k++) {
              fprintf(stdout,"%d %e  %e  %e\n",index_pk,pnl->k[index_k],pk_l[index_pk][index_k],pk_nl[index_pk][index_k]);
              }
              fprintf(stdout,"\n\n\n");
              }*/

            for (index_k=0; index_k<pnl->k_size; index_k++) {
              pnl->nl_corr_density[index_pk][index_tau * pnl->k_size + index_k] = sqrt(pk_nl[index_k]/pk_l[index_k]);
            }
          }
          else {
            /* when Halofit found k_max too small, use 1 as the
               non-linear correction for this redshift/time, store the
               last index which worked, and print a warning. */
            print_warning = _TRUE_;
            pnl->index_tau_min_nl = index_tau+1;
            for (index_k=0; index_k<pnl->k_size; index_k++) {
              pnl->nl_corr_density[index_pk][index_tau * pnl->k_size + index_k] = 1.;
            }
            if (pnl->nonlinear_verbose > 0) {
              class_alloc(pvecback,pba->bg_size*sizeof(double),pnl->error_message);
              class_call(background_at_tau(pba,pnl->tau[index_tau],pba->short_info,pba->inter_normal,&last_index,pvecback),
                         pba->error_message,
                         pnl->error_message);
              a = pvecback[pba->index_bg_a];
              z = pba->a_today/a-1.;
              fprintf(stdout,
                      " -> [WARNING:] index_pk=%d Halofit non-linear corrections could not be computed at redshift z=%5.2f and higher.\n    This is because k_max is too small for Halofit to be able to compute the scale k_NL at this redshift.\n    If non-linear corrections at such high redshift really matter for you,\n    just try to increase one of the parameters P_k_max_h/Mpc or P_k_max_1/Mpc or halofit_min_k_max (the code will take the max of these parameters) until reaching desired z.\n",
                      index_pk,z);
              free(pvecback);
            }
          }
        }
        else {
          /* if Halofit found k_max too small at a previous
             time/redhsift, use 1 as the non-linear correction for all
             higher redshifts/earlier times. */
          for (index_k=0; index_k<pnl->k_size; index_k++) {
            pnl->nl_corr_density[index_pk][index_tau * pnl->k_size + index_k] = 1.;
          }
        }

      }//end loop over pk_type

    }//end loop over tau

    /* free allocated arrays */

    free(pk_l);
    free(pk_nl);
    free(lnk_l);
    free(lnpk_l);
    free(ddlnpk_l);

  } // end of Halofit part

  else {
    class_stop(pnl->error_message,
               "Your non-linear method variable is set to %d, out of the range defined in nonlinear.h",pnl->method);
  }

  return _SUCCESS_;
}

/**
 * Free all memory space allocated by nonlinear_init().
 *
 *
 * @param pnl Input: pointer to nonlineard structure (to be freed)
 * @return the error status
 */

int nonlinear_free(
                   struct nonlinear *pnl
                   ) {
  int index_pk;

  if (pnl->method > nl_none) {

    if (pnl->method == nl_halofit) {
      free(pnl->k);
      free(pnl->tau);
      for(index_pk=0;index_pk<pnl->pk_size;++index_pk){
        free(pnl->nl_corr_density[index_pk]);
        free(pnl->k_nl[index_pk]);
      }
      free(pnl->nl_corr_density);
      free(pnl->k_nl);
    }
  }

  if (pnl->has_pk_eq == _TRUE_) {
    free(pnl->pk_eq_tau);
    free(pnl->pk_eq_w_and_Omega);
    free(pnl->pk_eq_ddw_and_ddOmega);
  }

  return _SUCCESS_;

}

/**
 * Calculation of the linear matter power spectrum, used to get the
 * nonlinear one.  This is partially redundent with a more elaborate
 * version of this calculation performed later in the spectra
 * module. At some point the organisation will change to avoid this
 * redundency.
 *
 * @param pba       Input: pointer to background structure
 * @param ppt       Input: pointer to perturbation structure
 * @param ppm       Input: pointer to primordial structure
 * @param pnl       Input: pointer to nonlinear structure
 * @param index_pk  Input: index of component are we looking at (total matter or cdm+baryons?)
 * @param index_tau Input: index of conformal time at which we want to do the calculation
 * @param pk_l      Output: linear spectrum at the relevant time
 * @param lnk       Output: array log(wavenumber)
 * @param lnpk      Output: array of log(P(k)_linear)
 * @param ddlnpk    Output: array of second derivative of log(P(k)_linear) wrt k, for spline interpolation
 * @return the error status
 */

int nonlinear_pk_l(
                   struct background *pba,
                   struct perturbs *ppt,
                   struct primordial *ppm,
                   struct nonlinear *pnl,
                   int index_pk,
                   int index_tau,
                   double *pk_l,
                   double *lnk,
                   double *lnpk,
                   double *ddlnpk) {

  int index_md;
  int index_k;
  int index_delta;
  int index_ic1,index_ic2,index_ic1_ic2;
  double * primordial_pk;
  double source_ic1,source_ic2;

  index_md = ppt->index_md_scalars;

  if ((pnl->has_pk_m == _TRUE_) && (index_pk == pnl->index_pk_m)) {
    index_delta = ppt->index_tp_delta_m;
  }
  else if ((pnl->has_pk_cb == _TRUE_) && (index_pk == pnl->index_pk_cb)) {
    index_delta = ppt->index_tp_delta_cb;
  }
  else {
    class_stop(pnl->error_message,"P(k) is set neither to total matter nor to cold dark matter + baryons");
  }

  class_alloc(primordial_pk,ppm->ic_ic_size[index_md]*sizeof(double),pnl->error_message);

  for (index_k=0; index_k<pnl->k_size; index_k++) {

    class_call(primordial_spectrum_at_k(ppm,
                                        index_md,
                                        linear,
                                        pnl->k[index_k],
                                        primordial_pk),
               ppm->error_message,
               pnl->error_message);

    pk_l[index_k] = 0;

    /* part diagonal in initial conditions */
    for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++) {

      index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic1,ppm->ic_size[index_md]);

      source_ic1 = ppt->sources[index_md]
        [index_ic1 * ppt->tp_size[index_md] + index_delta]
        [index_tau * ppt->k_size[index_md] + index_k];

      pk_l[index_k] += 2.*_PI_*_PI_/pow(pnl->k[index_k],3)
        *source_ic1*source_ic1
        *primordial_pk[index_ic1_ic2];

    }

    /* part non-diagonal in initial conditions */
    for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++) {
      for (index_ic2 = index_ic1+1; index_ic2 < ppm->ic_size[index_md]; index_ic2++) {

        index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,ppm->ic_size[index_md]);

        if (ppm->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {

          source_ic1 = ppt->sources[index_md]
            [index_ic1 * ppt->tp_size[index_md] + index_delta]
            [index_tau * ppt->k_size[index_md] + index_k];

          source_ic2 = ppt->sources[index_md]
            [index_ic2 * ppt->tp_size[index_md] + index_delta]
            [index_tau * ppt->k_size[index_md] + index_k];

          pk_l[index_k] += 2.*2.*_PI_*_PI_/pow(pnl->k[index_k],3)
            *source_ic1*source_ic2
            *primordial_pk[index_ic1_ic2]; // extra 2 factor (to include the symmetric term ic2,ic1)

        }
      }
    }

    lnk[index_k] = log(pnl->k[index_k]);
    lnpk[index_k] = log(pk_l[index_k]);

  }

  class_call(array_spline_table_columns(lnk,
                                        pnl->k_size,
                                        lnpk,
                                        1,
                                        ddlnpk,
                                        _SPLINE_NATURAL_,
                                        pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  free(primordial_pk);

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
 * @param pnl         Input/Output: pointer to nonlinear structure
 * @param index_pk    Input: index of component are we looking at (total matter or cdm+baryons?)
 * @param tau         Input: conformal time at which we want to do the calculation
 * @param pk_l        Input: linear spectrum at the relevant time
 * @param pk_nl       Output: non linear spectrum at the relevant time
 * @param lnk_l       Input: array log(wavenumber)
 * @param lnpk_l      Input: array of log(P(k)_linear)
 * @param ddlnpk_l    Input: array of second derivative of log(P(k)_linear) wrt k, for spline interpolation
 * @param k_nl        Output: non-linear wavenumber
 * @param halofit_found_k_max Ouput: flag cocnerning the status of the calculation (_FALSE_ if not possible)
 * @return the error status
 */

int nonlinear_halofit(
                      struct precision *ppr,
                      struct background *pba,
                      struct perturbs *ppt,
                      struct primordial *ppm,
                      struct nonlinear *pnl,
                      int index_pk,
                      double tau,
                      double *pk_l,
                      double *pk_nl,
                      double *lnk_l,
                      double *lnpk_l,
                      double *ddlnpk_l,
                      double *k_nl,
                      short * halofit_found_k_max
                      ) {

  double Omega_m,Omega_v,fnu,Omega0_m, w0, dw_over_da_fld, integral_fld;

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

  class_alloc(pvecback,pba->bg_size*sizeof(double),pnl->error_message);

  Omega0_m = pba->Omega0_cdm + pba->Omega0_b + pba->Omega0_ncdm_tot + pba->Omega0_dcdm;

  if ((pnl->has_pk_m == _TRUE_) && (index_pk == pnl->index_pk_m)) {
    fnu = pba->Omega0_ncdm_tot/Omega0_m;
  }
  else if ((pnl->has_pk_cb == _TRUE_) && (index_pk == pnl->index_pk_cb)) {
    fnu = 0.;
  }
  else {
    class_stop(pnl->error_message,"P(k) is set neither to total matter nor to cold dark matter + baryons");
  }

  if (pnl->has_pk_eq == _FALSE_) {

    /* default method to compute w0 = w_fld today, Omega_m(tau) and Omega_v=Omega_DE(tau),
       all required by HALFIT fitting formulas */

    class_call(background_w_fld(pba,pba->a_today,&w0,&dw_over_da_fld,&integral_fld), pba->error_message, pnl->error_message);

    class_call(background_at_tau(pba,tau,pba->long_info,pba->inter_normal,&last_index,pvecback),
               pba->error_message,
               pnl->error_message);

    Omega_m = pvecback[pba->index_bg_Omega_m];
    Omega_v = 1.-pvecback[pba->index_bg_Omega_m]-pvecback[pba->index_bg_Omega_r];

  }
  else {

    /* alternative method called Pk_equal, described in 0810.0190 and
                      1601.07230, extending the range of validity of
                      HALOFIT from constant w to (w0,wa) models. In that
                      case, some effective values of w0(tau_i) and
                      Omega_m(tau_i) have been pre-computed in the
                      input module, and we just ned to interpolate
                      within tabulated arrays, to get them at the
                      current tau value. */

    class_alloc(w_and_Omega,pnl->pk_eq_size*sizeof(double),pnl->error_message);

    class_call(array_interpolate_spline(
                                        pnl->pk_eq_tau,
                                        pnl->pk_eq_tau_size,
                                        pnl->pk_eq_w_and_Omega,
                                        pnl->pk_eq_ddw_and_ddOmega,
                                        pnl->pk_eq_size,
                                        tau,
                                        &last_index,
                                        w_and_Omega,
                                        pnl->pk_eq_size,
                                        pnl->error_message),
               pnl->error_message,
               pnl->error_message);

    w0 = w_and_Omega[pnl->index_pk_eq_w];
    Omega_m = w_and_Omega[pnl->index_pk_eq_Omega_m];
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

  integrand_size=(int)(log(pnl->k[pnl->k_size-1]/pnl->k[0])/log(10.)*ppr->halofit_k_per_decade)+1;

  class_alloc(integrand_array,integrand_size*ia_size*sizeof(double),pnl->error_message);

  //fprintf(stderr,"Omega_m=%e,  fnu=%e\n",Omega0_m,fnu);

  /* we fill integrand_array with values of k and P(k) using interpolation */

  last_index=0;

  for (index_k=0; index_k < integrand_size; index_k++) {

    k_integrand=pnl->k[0]*pow(10.,index_k/ppr->halofit_k_per_decade);

    class_call(array_interpolate_spline(lnk_l,
                                        pnl->k_size,
                                        lnpk_l,
                                        ddlnpk_l,
                                        1,
                                        log(k_integrand),
                                        &last_index,
                                        &lnpk_integrand,
                                        1,
                                        pnl->error_message),
               pnl->error_message,
               pnl->error_message);

    integrand_array[index_k*ia_size + index_ia_k] = k_integrand;
    integrand_array[index_k*ia_size + index_ia_pk] = exp(lnpk_integrand);

  }

  class_call(background_at_tau(pba,tau,pba->long_info,pba->inter_normal,&last_index,pvecback),
             pba->error_message,
             pnl->error_message);

  Omega_m = pvecback[pba->index_bg_Omega_m];
  Omega_v = 1.-pvecback[pba->index_bg_Omega_m]-pvecback[pba->index_bg_Omega_r];

  // for debugging:
  //printf("Call Halofit at z=%e\n",pba->a_today/pvecback[pba->index_bg_a]-1.);

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

  class_call(nonlinear_halofit_integrate(
                                         pnl,
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
             pnl->error_message,
             pnl->error_message);

  sigma  = sqrt(sum1);

  /* the following error should not stop the code: it will arrive
     inevitably at some large redshift, and then the code should not
     stop, but just give up computing P_NL(k,z). This is why we have a
     special error handling here (using class_test_except and free()
     commands to avoid memory leaks, and calling this whole function
     not through a class_call) */

  /*
  class_test_except(sigma < 1.,
                    pnl->error_message,
                    free(pvecback);free(integrand_array),
                    "Your k_max=%g 1/Mpc is too small for Halofit to find the non-linearity scale z_nl at z=%g. Increase input parameter P_k_max_h/Mpc or P_k_max_1/Mpc",
                    pnl->k[pnl->k_size-1],
                    pba->a_today/pvecback[pba->index_bg_a]-1.);
  */

  if (sigma < 1.) {
    * halofit_found_k_max = _FALSE_;
    free(pvecback);
    free(integrand_array);
    return _SUCCESS_;
  }
  else {
    * halofit_found_k_max = _TRUE_;
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
  class_call(nonlinear_halofit_integrate(
                                         pnl,
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
             pnl->error_message,
             pnl->error_message);

  sigma  = sqrt(sum1);

  class_test(sigma > 1.,
             pnl->error_message,
             "Your input value for the precision parameter halofit_min_k_nonlinear=%e is too large, such that sigma(R=1/halofit_min_k_nonlinear)=% > 1. For self-consistency, it should have been <1. Decrease halofit_min_k_nonlinear",
             ppr->halofit_min_k_nonlinear,sigma);

  xlogr2 = log(R)/log(10.);

  counter = 0;
  do {
    rmid = pow(10,(xlogr2+xlogr1)/2.0);
    counter ++;

    class_call(nonlinear_halofit_integrate(
                                           pnl,
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
             pnl->error_message,
             pnl->error_message);

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
                      pnl->error_message,
                      free(pvecback);free(integrand_array),
                      "could not converge within maximum allowed number of iterations");
    */
    /* ... but in this situation it sounds better to make it stop and return an error! */
    class_test(counter > _MAX_IT_,
               pnl->error_message,
               "could not converge within maximum allowed number of iterations");

  } while (fabs(diff) > ppr->halofit_tol_sigma);

  /* evaluate all the other integrals at R=rmid */

  class_call(nonlinear_halofit_integrate(
                                         pnl,
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
             pnl->error_message,
             pnl->error_message);

  class_call(nonlinear_halofit_integrate(
                                         pnl,
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
             pnl->error_message,
             pnl->error_message);

  sigma  = sqrt(sum1);
  d1 = -sum2/sum1;
  d2 = -sum2*sum2/sum1/sum1 - sum3/sum1;

  rknl  = 1./rmid;
  rneff = -3.-d1;
  rncur = -d2;

  *k_nl = rknl;

  //fprintf(stderr,"Here\n");

  for (index_k = 0; index_k < pnl->k_size; index_k++){

    rk = pnl->k[index_k];

    if (rk > ppr->halofit_min_k_nonlinear) {

      pk_lin = pk_l[index_k]*pow(pnl->k[index_k],3)*anorm;

      /* in original halofit, this is the beginning of the function halofit() */

      /*SPB11: Standard halofit underestimates the power on the smallest
       * scales by a factor of two. Add an extra correction from the
       * simulations in Bird, Viel,Haehnelt 2011 which partially accounts for
       * this.*/
      /*SPB14: This version of halofit is an updated version of the fit to the massive neutrinos
       * based on the results of Takahashi 2012, (arXiv:1208.2701).
       */
      gam=0.1971-0.0843*rneff+0.8460*rncur;
      a=1.5222+2.8553*rneff+2.3706*rneff*rneff+0.9903*rneff*rneff*rneff+ 0.2250*rneff*rneff*rneff*rneff-0.6038*rncur+0.1749*Omega_v*(1.+w0);
      a=pow(10,a);
      b=pow(10, (-0.5642+0.5864*rneff+0.5716*rneff*rneff-1.5474*rncur+0.2279*Omega_v*(1.+w0)));
      c=pow(10, 0.3698+2.0404*rneff+0.8161*rneff*rneff+0.5869*rncur);
      xmu=0.;
      xnu=pow(10,5.2105+3.6902*rneff);
      alpha=fabs(6.0835+1.3373*rneff-0.1959*rneff*rneff-5.5274*rncur);
      beta=2.0379-0.7354*rneff+0.3157*pow(rneff,2)+1.2490*pow(rneff,3)+0.3980*pow(rneff,4)-0.1682*rncur + fnu*(1.081 + 0.395*pow(rneff,2));

      if(fabs(1-Omega_m)>0.01) { /*then omega evolution */
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
      pk_halo=pk_halo/(1+xmu*pow(y,-1)+xnu*pow(y,-2))*(1+fnu*(0.977-18.015*(Omega0_m-0.3)));
      // rk is in 1/Mpc, 47.48and 1.5 in Mpc**-2, so we need an h**2 here (Credits Antonio J. Cuesta)
      pk_linaa=pk_lin*(1+fnu*47.48*pow(rk/pba->h,2)/(1+1.5*pow(rk/pba->h,2)));
      pk_quasi=pk_lin*pow((1+pk_linaa),beta)/(1+pk_linaa*alpha)*exp(-y/4.0-pow(y,2)/8.0);

      pk_nl[index_k] = (pk_halo+pk_quasi)/pow(pnl->k[index_k],3)/anorm;

      /* in original halofit, this is the end of the function halofit() */
    }
    else {
      pk_nl[index_k] = pk_l[index_k];
    }
  }

  free(pvecback);
  free(integrand_array);
  return _SUCCESS_;
}

/**
 * Internal routione of Halofit. In original Halofit, this is
 * equivalent to the function wint()
*/

int nonlinear_halofit_integrate(
                                struct nonlinear *pnl,
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
                            pnl->error_message),
               pnl->error_message,
               pnl->error_message);

    /* integrate */
    class_call(array_integrate_all_spline(integrand_array,
                                          ia_size,
                                          integrand_size,
                                          index_ia_k,
                                          index_ia_sum,
                                          index_ia_ddsum,
                                          sum,
                                          pnl->error_message),
               pnl->error_message,
               pnl->error_message);

  return _SUCCESS_;
}
