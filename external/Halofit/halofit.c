#include "halofit.h"
#include "fourier.h"

/**
 * Calculation of the nonlinear matter power spectrum with Halofit
 * (includes Takahashi 2012 + Bird 2013 revisions).
 *
 * At high redshift it is possible that the non-linear corrections are
 * so small that they can be computed only by going to very large
 * wavenumbers. Thus, for some combination of (z, k_max), the
 * calculation is not possible. In this case a _TRUE_ will be
 * returned in the flag nl_corr_not_computable_at_this_k.
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

int halofit(
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
            short *nl_corr_not_computable_at_this_k
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

  double *w_and_Omega;

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

  class_call(halofit_integrate(
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
  class_call(halofit_integrate(
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

    class_call(halofit_integrate(
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

  class_call(halofit_integrate(
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

  class_call(halofit_integrate(
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

int halofit_integrate(
                      struct fourier *pfo,
                      double *integrand_array,
                      int integrand_size,
                      int ia_size,
                      int index_ia_k,
                      int index_ia_pk,
                      int index_ia_sum,
                      int index_ia_ddsum,
                      double R,
                      enum halofit_integral_type type,
                      double *sum
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
