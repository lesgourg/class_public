/** @file cl.c Documented nonlinear module
 *
 * Julien Lesgourgues, 6.03.2014
 *
 * New module replacing an older one present up to version 2.0 The new
 * module is located in a better place in the main, allowing it to
 * compute non-linear correction to Cl's and not just P(k).  It will
 * also be easier to generalise to new methods.  The old implemenation
 * of one-loop calculations and TRG calculations has been dropped from
 * this version, they can still be found in older versions.
 *
 */

#include "nonlinear.h"

int nonlinear_k_nl_at_z(
                        struct background *pba,
                        struct nonlinear * pnl,
                        double z,
                        double * k_nl
                        ) {

  double tau;

  class_call(background_tau_of_z(pba,
                                 z,
                                 &tau),
             pba->error_message,
             pnl->error_message);

  if (pnl->tau_size == 1) {
    *k_nl = pnl->k_nl[0];
  }
  else {
    class_call(array_interpolate_two(pnl->tau,
                                     1,
                                     0,
                                     pnl->k_nl,
                                     1,
                                     pnl->tau_size,
                                     tau,
                                     k_nl,
                                     1,
                                     pnl->error_message),
               pnl->error_message,
               pnl->error_message);
  }

  return _SUCCESS_;
}

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
  double *pk_l;
  double *pk_nl;

  /** (a) if non non-linear corrections requested */

  if (pnl->method == nl_none) {
    if (pnl->nonlinear_verbose > 0)
      printf("No non-linear spectra requested. Nonlinear module skipped.\n");
  }

  /** (b) for HALOFIT non-linear spectrum */

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

    class_alloc(pnl->nl_corr_density,pnl->tau_size*pnl->k_size*sizeof(double),pnl->error_message);
    class_alloc(pnl->k_nl,pnl->tau_size*sizeof(double),pnl->error_message);

    class_alloc(pk_l,pnl->k_size*sizeof(double),pnl->error_message);
    class_alloc(pk_nl,pnl->k_size*sizeof(double),pnl->error_message);

    /** - loop over time */

    for (index_tau = pnl->tau_size-1; index_tau>=0; index_tau--) {

      /* get P_L(k) at this time */
      class_call(nonlinear_pk_l(ppt,ppm,pnl,index_tau,pk_l),
                 pnl->error_message,
                 pnl->error_message);

      /*
      for (index_k=0; index_k<pnl->k_size; index_k++) {
        fprintf(stdout,"%e  %e\n",pnl->k[index_k],pk_l[index_k]);
      }
      */

      //class_stop(pnl->error_message,"stop here");

      /* get P_NL(k) at this time */
      if (nonlinear_halofit(ppr,
                            pba,
                            ppm,
                            pnl,
                            pnl->tau[index_tau],
                            pk_l,
                            pk_nl,
                            &(pnl->k_nl[index_tau])) == _SUCCESS_) {

        /*
          for (index_k=0; index_k<pnl->k_size; index_k++) {
          fprintf(stdout,"%e  %e  %e\n",pnl->k[index_k],pk_l[index_k],pk_nl[index_k]);
          }
          fprintf(stdout,"\n\n");
        */

        for (index_k=0; index_k<pnl->k_size; index_k++) {
          pnl->nl_corr_density[index_tau * pnl->k_size + index_k] = sqrt(pk_nl[index_k]/pk_l[index_k]);
        }
      }
      else {
        for (index_k=0; index_k<pnl->k_size; index_k++) {
          pnl->nl_corr_density[index_tau * pnl->k_size + index_k] = 1.;
        }
      }
    }

    /*
    for (index_tau = pnl->tau_size-1; index_tau>=0; index_tau--) {
      for (index_k=0; index_k<pnl->k_size; index_k++) {
        fprintf(stdout,"%e  %e\n",pnl->k[index_k],pnl->nl_corr_density[index_tau * pnl->k_size + index_k]);
      }
      fprintf(stdout,"\n\n");
    }
    */

    free(pk_l);
    free(pk_nl);
  }

  else {
    class_stop(pnl->error_message,
               "Your non-linear method variable is set to %d, out of the range defined in nonlinear.h",pnl->method);
  }

  return _SUCCESS_;
}

int nonlinear_free(
                   struct nonlinear *pnl
                   ) {

  if (pnl->method > nl_none) {

    if (pnl->method == nl_halofit) {
      /* free here */
      free(pnl->k);
      free(pnl->tau);
      free(pnl->nl_corr_density);
      free(pnl->k_nl);
    }
  }

  return _SUCCESS_;

}

int nonlinear_pk_l(
                   struct perturbs *ppt,
                   struct primordial *ppm,
                   struct nonlinear *pnl,
                   int index_tau,
                   double *pk_l
                   ) {

  int index_md;
  int index_k;
  int index_ic1,index_ic2,index_ic1_ic2;
  double * primordial_pk;
  double source_ic1,source_ic2;

  index_md = ppt->index_md_scalars;

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
        [index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_m]
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
            [index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_m]
            [index_tau * ppt->k_size[index_md] + index_k];

          source_ic2 = ppt->sources[index_md]
            [index_ic2 * ppt->tp_size[index_md] + ppt->index_tp_delta_m]
            [index_tau * ppt->k_size[index_md] + index_k];

          pk_l[index_k] += 2.*2.*_PI_*_PI_/pow(pnl->k[index_k],3)
            *source_ic1*source_ic2
            *primordial_pk[index_ic1_ic2]; // extra 2 factor (to include the symetric term ic2,ic1)

        }
      }
    }

  }

  free(primordial_pk);

  return _SUCCESS_;

}

int nonlinear_halofit(
                      struct precision *ppr,
                      struct background *pba,
                      struct primordial *ppm,
                      struct nonlinear *pnl,
                      double tau,
                      double *pk_l,
                      double *pk_nl,
                      double *k_nl
                      ) {

  double Omega_m,Omega_v,fnu,Omega0_m, w0;

  /** determine non linear ratios (from pk) **/

  int index_k;
  double pk_lin,pk_quasi,pk_halo,rk;
  double sigma,rknl,rneff,rncur,d1,d2;
  double diff,xlogr1,xlogr2,rmid;

  double gam,a,b,c,xmu,xnu,alpha,beta,f1,f2,f3;
  double pk_linaa;
  double y;
  double f1a,f2a,f3a,f1b,f2b,f3b,frac;

  double * pvecback;

  int last_index;
  int counter;
  double sum1,sum2,sum3;
  double anorm;

  double *integrand_array;

  double x2,R;

  class_alloc(pvecback,pba->bg_size*sizeof(double),pnl->error_message);

  Omega0_m = (pba->Omega0_cdm + pba->Omega0_b + pba->Omega0_ncdm_tot);
  w0 = pba->w0_fld;
  fnu      = pba->Omega0_ncdm_tot/(pba->Omega0_b + pba->Omega0_cdm+pba->Omega0_ncdm_tot);
  anorm    = 1./(2*pow(_PI_,2));

  class_alloc(integrand_array,pnl->k_size*7*sizeof(double),pnl->error_message);

  for (index_k=0; index_k < pnl->k_size; index_k++) {
    integrand_array[index_k*7 + 0] = pnl->k[index_k];
  }

  class_call(background_at_tau(pba,tau,pba->long_info,pba->inter_normal,&last_index,pvecback),
             pba->error_message,
             pnl->error_message);

  Omega_m = pvecback[pba->index_bg_Omega_m];
  Omega_v = 1.-pvecback[pba->index_bg_Omega_m]-pvecback[pba->index_bg_Omega_r];

  /* minimum value of R such that the integral giving sigma_R is converged */
  R=sqrt(-log(ppr->halofit_sigma_precision))/pnl->k[pnl->k_size-1];

  /* corresponding value of sigma_R */
  sum1=0.;
  for (index_k=0; index_k < pnl->k_size; index_k++) {
    x2 = pow(pnl->k[index_k]*R,2);
    integrand_array[index_k*7 + 1] = pk_l[index_k]*pow(pnl->k[index_k],2)*anorm*exp(-x2);
  }
  /* fill in second derivatives */
  class_call(array_spline(integrand_array,7,pnl->k_size,0,1,2,_SPLINE_EST_DERIV_,pnl->error_message),
             pnl->error_message,
             pnl->error_message);
  /* integrate */
  class_call(array_integrate_all_spline(integrand_array,7,pnl->k_size,0,1,2,&sum1,pnl->error_message),
             pnl->error_message,
             pnl->error_message);
  sigma  = sqrt(sum1);

  class_test_except(sigma < 1.,
                    pnl->error_message,
                    free(pvecback);free(integrand_array),
                    "Your k_max=%g 1/Mpc is too small for Halofit to find the non-linearity scale z_nl at z=%g. Increase input parameter P_k_max_h/Mpc or P_k_max_1/Mpc",
                    pnl->k[pnl->k_size-1],
                    pba->a_today/pvecback[pba->index_bg_a]-1.);

  xlogr1 = log(R)/log(10.);

  /* maximum value of R in the bisection algorithm leading to the determination of R_nl */
  R=1./ppr->halofit_min_k_nonlinear;

  /* corresponding value of sigma_R */
  sum1=0.;
  for (index_k=0; index_k < pnl->k_size; index_k++) {
    x2 = pow(pnl->k[index_k]*R,2);
    integrand_array[index_k*7 + 1] = pk_l[index_k]*pow(pnl->k[index_k],2)*anorm*exp(-x2);
  }
  /* fill in second derivatives */
  class_call(array_spline(integrand_array,7,pnl->k_size,0,1,2,_SPLINE_EST_DERIV_,pnl->error_message),
             pnl->error_message,
             pnl->error_message);
  /* integrate */
  class_call(array_integrate_all_spline(integrand_array,7,pnl->k_size,0,1,2,&sum1,pnl->error_message),
             pnl->error_message,
             pnl->error_message);
  sigma  = sqrt(sum1);

  class_test_except(sigma > 1.,
                    pnl->error_message,
                    free(pvecback);free(integrand_array),
                    "Your input value for the precision parameter halofit_min_k_nonlinear=%e is too large, the non-linear wavenumber k_nl must be smaller than that",
                    ppr->halofit_min_k_nonlinear);

  xlogr2 = log(R)/log(10.);

  counter = 0;
  do {
    rmid = pow(10,(xlogr2+xlogr1)/2.0);
    counter ++;

    /* in original halofit, this is the function wint() */
    sum1=0.;
    sum2=0.;
    sum3=0.;

    for (index_k=0; index_k < pnl->k_size; index_k++) {
      x2 = pnl->k[index_k]*pnl->k[index_k]*rmid*rmid;
      integrand_array[index_k*7 + 1] = pk_l[index_k]*pow(pnl->k[index_k],2)*anorm*exp(-x2);
      integrand_array[index_k*7 + 3] = pk_l[index_k]*pow(pnl->k[index_k],2)*anorm*2.*x2*exp(-x2);
      integrand_array[index_k*7 + 5] = pk_l[index_k]*pow(pnl->k[index_k],2)*anorm*4.*x2*(1.-x2)*exp(-x2);
    }
    /* fill in second derivatives */
    class_call(array_spline(integrand_array,7,pnl->k_size,0,1,2,_SPLINE_EST_DERIV_,pnl->error_message),
               pnl->error_message,
               pnl->error_message);
    class_call(array_spline(integrand_array,7,pnl->k_size,0,3,4,_SPLINE_EST_DERIV_,pnl->error_message),
               pnl->error_message,
               pnl->error_message);
    class_call(array_spline(integrand_array,7,pnl->k_size,0,5,6,_SPLINE_EST_DERIV_,pnl->error_message),
               pnl->error_message,
               pnl->error_message);

    /* integrate */
    class_call(array_integrate_all_spline(integrand_array,7,pnl->k_size,0,1,2,&sum1,pnl->error_message),
               pnl->error_message,
               pnl->error_message);
    class_call(array_integrate_all_spline(integrand_array,7,pnl->k_size,0,3,4,&sum2,pnl->error_message),
               pnl->error_message,
               pnl->error_message);
    class_call(array_integrate_all_spline(integrand_array,7,pnl->k_size,0,5,6,&sum3,pnl->error_message),
               pnl->error_message,
               pnl->error_message);

    sigma  = sqrt(sum1);
    d1 = -sum2/sum1;
    d2 = -sum2*sum2/sum1/sum1 - sum3/sum1;
    /* in original halofit, this is the end of the function wint() */

    diff = sigma - 1.0;
    /*fprintf(stderr,"xlogr1 = %g, xlogr2 = %g, rmid = %g, diff: =%g, abs(diff) = %g\n",xlogr1,xlogr2,log10(rmid),diff,fabs(diff));*/
    if (diff>0.001){
      xlogr1=log10(rmid);
      /*fprintf(stderr,"going up  , new xlogr1=%g\n",xlogr1);*/
    }
    else if (diff < -0.001) {
      xlogr2 = log10(rmid);
      /*fprintf(stderr,"going down, new xlogr2=%g\n",xlogr2);*/
    }

    class_test_except(counter > _MAX_IT_,
                      pnl->error_message,
                      free(pvecback);free(integrand_array),
                      "could not converge within maximum allowed number of iterations");

  } while (fabs(diff) > 0.001);

  rknl  = 1./rmid;
  rneff = -3.-d1;
  rncur = -d2;

  *k_nl = rknl;

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
      alpha=abs(6.0835+1.3373*rneff-0.1959*rneff*rneff-5.5274*rncur);
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
      pk_linaa=pk_lin*(1+fnu*47.48*pow(rk,2)/(1+1.5*pow(rk,2)));
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
