/** @file lensing.c Documented lensing module
 *
 * Simon Prunet and Julien Lesgourgues, 6.12.2010
 *
 * This module computes the lensed temperature and polarization
 * anisotropy power spectra \f$ C_l^{X}, P(k), ... \f$'s given the
 * unlensed temperature, polarization and lensing potential spectra.
 *
 * Follows Challinor and Lewis full-sky method, astro-ph/0502425
 *
 * The following functions can be called from other modules:
 *
 * -# lensing_init() at the beginning (but after harmonic_init())
 * -# lensing_cl_at_l() at any time for computing Cl_lensed at any l
 * -# lensing_free() at the end
 */

#include "lensing.h"
#include <time.h>

/**
 * Anisotropy power spectra \f$ C_l\f$'s for all types, modes and initial conditions.
 * SO FAR: ONLY SCALAR
 *
 * This routine evaluates all the lensed \f$ C_l\f$'s at a given value of l by
 * picking it in the pre-computed table. When relevant, it also
 * sums over all initial conditions for each mode, and over all modes.
 *
 * This function can be called from whatever module at whatever time,
 * provided that lensing_init() has been called before, and
 * lensing_free() has not been called yet.
 *
 * @param ple        Input: pointer to lensing structure
 * @param l          Input: multipole number
 * @param cl_lensed  Output: lensed \f$ C_l\f$'s for all types (TT, TE, EE, etc..)
 * @return the error status
 */

int lensing_cl_at_l(
                    struct lensing * ple,
                    int l,
                    double * cl_lensed    /* array with argument cl_lensed[index_ct] (must be already allocated) */
                    ) {
  int last_index;
  int index_lt;

  class_test(l > ple->l_lensed_max,
             ple->error_message,
             "you asked for lensed Cls at l=%d, they were computed only up to l=%d, you should increase l_max_scalars or decrease the precision parameter delta_l_max",l,ple->l_lensed_max);

  class_call(array_interpolate_spline(ple->l,
                                      ple->l_size,
                                      ple->cl_lens,
                                      ple->ddcl_lens,
                                      ple->lt_size,
                                      l,
                                      &last_index,
                                      cl_lensed,
                                      ple->lt_size,
                                      ple->error_message),
             ple->error_message,
             ple->error_message);

  /* set to zero for the types such that l<l_max */
  for (index_lt=0; index_lt<ple->lt_size; index_lt++)
    if ((int)l > ple->l_max_lt[index_lt])
      cl_lensed[index_lt]=0.;

  return _SUCCESS_;
}

/**
 * This routine initializes the lensing structure (in particular,
 * computes table of lensed anisotropy spectra \f$ C_l^{X} \f$)
 *
 * @param ppr Input: pointer to precision structure
 * @param ppt Input: pointer to perturbation structure (just in case, not used in current version...)
 * @param phr Input: pointer to harmonic structure
 * @param pfo Input: pointer to fourier structure
 * @param ple Output: pointer to initialized lensing structure
 * @return the error status
 */

int lensing_init(
                 struct precision * ppr,
                 struct perturbations * ppt,
                 struct harmonic * phr,
                 struct fourier * pfo,
                 struct lensing * ple
                 ) {

  /** Summary: */
  /** - Define local variables */

  double * mu; /* mu[index_mu]: discretized values of mu
                  between -1 and 1, roots of Legendre polynomial */
  double * w8; /* Corresponding Gauss-Legendre quadrature weights */
  double theta,delta_theta;

  double ** d00;  /* dmn[index_mu][index_l] */
  double ** d11;
  double ** d2m2;
  double ** d22 = NULL;
  double ** d20 = NULL;
  double ** d1m1;
  double ** d31 = NULL;
  double ** d40 = NULL;
  double ** d3m1 = NULL;
  double ** d3m3 = NULL;
  double ** d4m2 = NULL;
  double ** d4m4 = NULL;
  double * buf_dxx; /* buffer */

  double * Cgl;   /* Cgl[index_mu] */
  double * Cgl2;  /* Cgl2[index_mu] */
  double * sigma2; /* sigma[index_mu] */

  double * ksi = NULL;  /* ksi[index_mu] */
  double * ksiX = NULL;  /* ksiX[index_mu] */
  double * ksip = NULL;  /* ksip[index_mu] */
  double * ksim = NULL;  /* ksim[index_mu] */

  double fac,fac1;
  double X_000;
  double X_p000;
  double X_220;
  double X_022;
  double X_p022;
  double X_121;
  double X_132;
  double X_242;

  int num_mu,index_mu,icount;
  int l;
  double ll;
  double * cl_unlensed;  /* cl_unlensed[index_ct] */
  double * cl_tt; /* unlensed  cl, to be filled to avoid repeated calls to harmonic_cl_at_l */
  double * cl_te = NULL; /* unlensed  cl, to be filled to avoid repeated calls to harmonic_cl_at_l */
  double * cl_ee = NULL; /* unlensed  cl, to be filled to avoid repeated calls to harmonic_cl_at_l */
  double * cl_bb = NULL; /* unlensed  cl, to be filled to avoid repeated calls to harmonic_cl_at_l */
  double * cl_pp; /* potential cl, to be filled to avoid repeated calls to harmonic_cl_at_l */

  double res,resX,lens;
  double resp, resm, lensp, lensm;

  double * sqrt1;
  double * sqrt2;
  double * sqrt3;
  double * sqrt4;
  double * sqrt5;

  double ** cl_md_ic; /* array with argument
                         cl_md_ic[index_md][index_ic1_ic2*phr->ct_size+index_ct] */

  double ** cl_md;    /* array with argument
                         cl_md[index_md][index_ct] */

  int index_md;

  /* Timing */
  //double debut, fin;
  //double cpu_time;

  /** - check that we really want to compute at least one spectrum */

  if (ple->has_lensed_cls == _FALSE_) {
    if (ple->lensing_verbose > 0)
      printf("No lensing requested. Lensing module skipped.\n");
    return _SUCCESS_;
  }
  else {
    if (ple->lensing_verbose > 0) {
      printf("Computing lensed spectra ");
      if (ppr->accurate_lensing==_TRUE_)
        printf("(accurate mode)\n");
      else
        printf("(fast mode)\n");
    }
  }

  /** - initialize indices and allocate some of the arrays in the
      lensing structure */

  class_call(lensing_indices(ppr,phr,ple),
             ple->error_message,
             ple->error_message);

  /** - put all precision variables hare; will be stored later in precision structure */
  /** - Last element in \f$ \mu \f$ will be for \f$ \mu=1 \f$, needed for sigma2.
      The rest will be chosen as roots of a Gauss-Legendre quadrature **/

  if (ppr->accurate_lensing == _TRUE_) {
    num_mu=(ple->l_unlensed_max+ppr->num_mu_minus_lmax); /* Must be even ?? CHECK */
    num_mu += num_mu%2; /* Force it to be even */
  } else {
    /* Integrate correlation function difference on [0,pi/16] */
    num_mu = (ple->l_unlensed_max * 2 )/16;
  }
  /** - allocate array of \f$ \mu \f$ values, as well as quadrature weights */

  class_alloc(mu,
              num_mu*sizeof(double),
              ple->error_message);
  /* Reserve last element of mu for mu=1, needed for sigma2 */
  mu[num_mu-1] = 1.0;

  class_alloc(w8,
              (num_mu-1)*sizeof(double),
              ple->error_message);

  if (ppr->accurate_lensing == _TRUE_) {

    //debut = omp_get_wtime();
    class_call(quadrature_gauss_legendre(mu,
                                         w8,
                                         num_mu-1,
                                         ppr->tol_gauss_legendre,
                                         ple->error_message),
               ple->error_message,
               ple->error_message);
    //fin = omp_get_wtime();
    //cpu_time = (fin-debut);
    //printf("time in quadrature_gauss_legendre=%4.3f s\n",cpu_time);

  } else { /* Crude integration on [0,pi/16]: Riemann sum on theta */

    delta_theta = _PI_/16. / (double)(num_mu-1);
    for (index_mu=0;index_mu<num_mu-1;index_mu++) {
      theta = (index_mu+1)*delta_theta;
      mu[index_mu] = cos(theta);
      w8[index_mu] = sin(theta)*delta_theta; /* We integrate on mu */
    }
  }

  /** - Compute \f$ d^l_{mm'} (\mu) \f$*/

  icount = 0;
  class_alloc(d00,
              num_mu*sizeof(double*),
              ple->error_message);

  class_alloc(d11,
              num_mu*sizeof(double*),
              ple->error_message);

  class_alloc(d1m1,
              num_mu*sizeof(double*),
              ple->error_message);

  class_alloc(d2m2,
              num_mu*sizeof(double*),
              ple->error_message);
  icount += 4*num_mu*(ple->l_unlensed_max+1);

  if(ple->has_te==_TRUE_) {

    class_alloc(d20,
                num_mu*sizeof(double*),
                ple->error_message);

    class_alloc(d3m1,
                num_mu*sizeof(double*),
                ple->error_message);

    class_alloc(d4m2,
                num_mu*sizeof(double*),
                ple->error_message);
    icount += 3*num_mu*(ple->l_unlensed_max+1);
  }

  if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {

    class_alloc(d22,
                num_mu*sizeof(double*),
                ple->error_message);

    class_alloc(d31,
                num_mu*sizeof(double*),
                ple->error_message);

    class_alloc(d3m3,
                num_mu*sizeof(double*),
                ple->error_message);

    class_alloc(d40,
                num_mu*sizeof(double*),
                ple->error_message);

    class_alloc(d4m4,
                num_mu*sizeof(double*),
                ple->error_message);
    icount += 5*num_mu*(ple->l_unlensed_max+1);
  }

  icount += 5*(ple->l_unlensed_max+1); /* for arrays sqrt1[l] to sqrt5[l] */

  /** - Allocate main contiguous buffer **/
  class_alloc(buf_dxx,
              icount * sizeof(double),
              ple->error_message);

  icount = 0;
  for (index_mu=0; index_mu<num_mu; index_mu++) {

    d00[index_mu] = &(buf_dxx[icount+index_mu            * (ple->l_unlensed_max+1)]);
    d11[index_mu] = &(buf_dxx[icount+(index_mu+num_mu)   * (ple->l_unlensed_max+1)]);
    d1m1[index_mu]= &(buf_dxx[icount+(index_mu+2*num_mu) * (ple->l_unlensed_max+1)]);
    d2m2[index_mu]= &(buf_dxx[icount+(index_mu+3*num_mu) * (ple->l_unlensed_max+1)]);
  }
  icount += 4*num_mu*(ple->l_unlensed_max+1);

  if (ple->has_te==_TRUE_) {
    for (index_mu=0; index_mu<num_mu; index_mu++) {
      d20[index_mu] = &(buf_dxx[icount+index_mu            * (ple->l_unlensed_max+1)]);
      d3m1[index_mu]= &(buf_dxx[icount+(index_mu+num_mu)   * (ple->l_unlensed_max+1)]);
      d4m2[index_mu]= &(buf_dxx[icount+(index_mu+2*num_mu) * (ple->l_unlensed_max+1)]);
    }
    icount += 3*num_mu*(ple->l_unlensed_max+1);
  }

  if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {

    for (index_mu=0; index_mu<num_mu; index_mu++) {
      d22[index_mu] = &(buf_dxx[icount+index_mu            * (ple->l_unlensed_max+1)]);
      d31[index_mu] = &(buf_dxx[icount+(index_mu+num_mu)   * (ple->l_unlensed_max+1)]);
      d3m3[index_mu]= &(buf_dxx[icount+(index_mu+2*num_mu) * (ple->l_unlensed_max+1)]);
      d40[index_mu] = &(buf_dxx[icount+(index_mu+3*num_mu) * (ple->l_unlensed_max+1)]);
      d4m4[index_mu]= &(buf_dxx[icount+(index_mu+4*num_mu) * (ple->l_unlensed_max+1)]);
    }
    icount += 5*num_mu*(ple->l_unlensed_max+1);
  }

  sqrt1 = &(buf_dxx[icount]);
  icount += ple->l_unlensed_max+1;
  sqrt2 = &(buf_dxx[icount]);
  icount += ple->l_unlensed_max+1;
  sqrt3 = &(buf_dxx[icount]);
  icount += ple->l_unlensed_max+1;
  sqrt4 = &(buf_dxx[icount]);
  icount += ple->l_unlensed_max+1;
  sqrt5 = &(buf_dxx[icount]);
  icount += ple->l_unlensed_max+1;

  //debut = omp_get_wtime();
  class_call(lensing_d00(mu,num_mu,ple->l_unlensed_max,d00),
             ple->error_message,
             ple->error_message);

  class_call(lensing_d11(mu,num_mu,ple->l_unlensed_max,d11),
             ple->error_message,
             ple->error_message);

  class_call(lensing_d1m1(mu,num_mu,ple->l_unlensed_max,d1m1),
             ple->error_message,
             ple->error_message);

  class_call(lensing_d2m2(mu,num_mu,ple->l_unlensed_max,d2m2),
             ple->error_message,
             ple->error_message);
  //fin = omp_get_wtime();
  //cpu_time = (fin-debut);
  //printf("time in lensing_dxx=%4.3f s\n",cpu_time);


  if (ple->has_te==_TRUE_) {

    class_call(lensing_d20(mu,num_mu,ple->l_unlensed_max,d20),
               ple->error_message,
               ple->error_message);

    class_call(lensing_d3m1(mu,num_mu,ple->l_unlensed_max,d3m1),
               ple->error_message,
               ple->error_message);

    class_call(lensing_d4m2(mu,num_mu,ple->l_unlensed_max,d4m2),
               ple->error_message,
               ple->error_message);

  }

  if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {

    class_call(lensing_d22(mu,num_mu,ple->l_unlensed_max,d22),
               ple->error_message,
               ple->error_message);

    class_call(lensing_d31(mu,num_mu,ple->l_unlensed_max,d31),
               ple->error_message,
               ple->error_message);

    class_call(lensing_d3m3(mu,num_mu,ple->l_unlensed_max,d3m3),
               ple->error_message,
               ple->error_message);

    class_call(lensing_d40(mu,num_mu,ple->l_unlensed_max,d40),
               ple->error_message,
               ple->error_message);

    class_call(lensing_d4m4(mu,num_mu,ple->l_unlensed_max,d4m4),
               ple->error_message,
               ple->error_message);
  }

  /** - compute \f$ Cgl(\mu)\f$, \f$ Cgl2(\mu) \f$ and sigma2(\f$\mu\f$) */

  class_alloc(Cgl,
              num_mu*sizeof(double),
              ple->error_message);

  class_alloc(Cgl2,
              num_mu*sizeof(double),
              ple->error_message);

  class_alloc(sigma2,
              (num_mu-1)*sizeof(double), /* Zero separation is omitted */
              ple->error_message);

  class_alloc(cl_unlensed,
              phr->ct_size*sizeof(double),
              ple->error_message);


  /** - Locally store unlensed temperature \f$ cl_{tt}\f$ and potential \f$ cl_{pp}\f$ spectra **/
  class_alloc(cl_tt,
              (ple->l_unlensed_max+1)*sizeof(double),
              ple->error_message);
  if (ple->has_te==_TRUE_) {
    class_alloc(cl_te,
                (ple->l_unlensed_max+1)*sizeof(double),
                ple->error_message);
  }
  if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
    class_alloc(cl_ee,
                (ple->l_unlensed_max+1)*sizeof(double),
                ple->error_message);

    class_alloc(cl_bb,
                (ple->l_unlensed_max+1)*sizeof(double),
                ple->error_message);
  }
  class_alloc(cl_pp,
              (ple->l_unlensed_max+1)*sizeof(double),
              ple->error_message);

  class_alloc(cl_md_ic,
              phr->md_size*sizeof(double *),
              ple->error_message);

  class_alloc(cl_md,
              phr->md_size*sizeof(double *),
              ple->error_message);

  for (index_md = 0; index_md < phr->md_size; index_md++) {

    if (phr->md_size > 1)

      class_alloc(cl_md[index_md],
                  phr->ct_size*sizeof(double),
                  ple->error_message);

    if (phr->ic_size[index_md] > 1)

      class_alloc(cl_md_ic[index_md],
                  phr->ic_ic_size[index_md]*phr->ct_size*sizeof(double),
                  ple->error_message);
  }

  for (l=2; l<=ple->l_unlensed_max; l++) {
    class_call(harmonic_cl_at_l(phr,l,cl_unlensed,cl_md,cl_md_ic),
               phr->error_message,
               ple->error_message);
    cl_tt[l] = cl_unlensed[ple->index_lt_tt];
    cl_pp[l] = cl_unlensed[ple->index_lt_pp];
    if (ple->has_te==_TRUE_) {
      cl_te[l] = cl_unlensed[ple->index_lt_te];
    }
    if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
      cl_ee[l] = cl_unlensed[ple->index_lt_ee];
      cl_bb[l] = cl_unlensed[ple->index_lt_bb];
    }
  }

  for (index_md = 0; index_md < phr->md_size; index_md++) {

    if (phr->md_size > 1)
      free(cl_md[index_md]);

    if (phr->ic_size[index_md] > 1)
      free(cl_md_ic[index_md]);

  }

  free(cl_md_ic);
  free(cl_md);

  /** - Compute sigma2\f$(\mu)\f$ and Cgl2(\f$\mu\f$) **/

  //debut = omp_get_wtime();
#pragma omp parallel for                        \
  private (index_mu,l)                          \
  schedule (static)
  for (index_mu=0; index_mu<num_mu; index_mu++) {

    Cgl[index_mu]=0;
    Cgl2[index_mu]=0;

    for (l=2; l<=ple->l_unlensed_max; l++) {

      Cgl[index_mu] += (2.*l+1.)*l*(l+1.)*
        cl_pp[l]*d11[index_mu][l];

      Cgl2[index_mu] += (2.*l+1.)*l*(l+1.)*
        cl_pp[l]*d1m1[index_mu][l];

    }

    Cgl[index_mu] /= 4.*_PI_;
    Cgl2[index_mu] /= 4.*_PI_;

  }

  for (index_mu=0; index_mu<num_mu-1; index_mu++) {
    /* Cgl(1.0) - Cgl(mu) */
    sigma2[index_mu] = Cgl[num_mu-1] - Cgl[index_mu];
  }
  //fin = omp_get_wtime();
  //cpu_time = (fin-debut);
  //printf("time in Cgl,Cgl2,sigma2=%4.3f s\n",cpu_time);


  /** - compute ksi, ksi+, ksi-, ksiX */

  /** - --> ksi is for TT **/
  if (ple->has_tt==_TRUE_) {

    class_calloc(ksi,
                 (num_mu-1),
                 sizeof(double),
                 ple->error_message);
  }

  /** - --> ksiX is for TE **/
  if (ple->has_te==_TRUE_) {

    class_calloc(ksiX,
                 (num_mu-1),
                 sizeof(double),
                 ple->error_message);
  }

  /** - --> ksip, ksim for EE, BB **/
  if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {

    class_calloc(ksip,
                 (num_mu-1),
                 sizeof(double),
                 ple->error_message);

    class_calloc(ksim,
                 (num_mu-1),
                 sizeof(double),
                 ple->error_message);
  }

  for (l=2;l<=ple->l_unlensed_max;l++) {

    ll = (double)l;
    sqrt1[l]=sqrt((ll+2)*(ll+1)*ll*(ll-1));
    sqrt2[l]=sqrt((ll+2)*(ll-1));
    sqrt3[l]=sqrt((ll+3)*(ll-2));
    sqrt4[l]=sqrt((ll+4)*(ll+3)*(ll-2.)*(ll-3));
    sqrt5[l]=sqrt(ll*(ll+1));
  }


  //debut = omp_get_wtime();
#pragma omp parallel for                                                \
  private (index_mu,l,ll,res,resX,resp,resm,lens,lensp,lensm,           \
           fac,fac1,X_000,X_p000,X_220,X_022,X_p022,X_121,X_132,X_242)	\
  schedule (static)

  for (index_mu=0;index_mu<num_mu-1;index_mu++) {

    for (l=2;l<=ple->l_unlensed_max;l++) {

      ll = (double)l;

      fac = ll*(ll+1)/4.;
      fac1 = (2*ll+1)/(4.*_PI_);

      /* In the following we will keep terms of the form (sigma2)^k*(Cgl2)^m
         with k+m <= 2 */

      X_000 = exp(-fac*sigma2[index_mu]);
      X_p000 = -fac*X_000;
      /* X_220 = 0.25*sqrt1[l] * exp(-(fac-0.5)*sigma2[index_mu]); */
      X_220 = 0.25*sqrt1[l] * X_000; /* Order 0 */
      /* next 5 lines useless, but avoid compiler warning 'may be used uninitialized' */
      X_242=0.;
      X_132=0.;
      X_121=0.;
      X_p022=0.;
      X_022=0.;

      if (ple->has_te==_TRUE_ || ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
        /* X_022 = exp(-(fac-1.)*sigma2[index_mu]); */
        X_022 = X_000 * (1+sigma2[index_mu]*(1+0.5*sigma2[index_mu])); /* Order 2 */
        X_p022 = -(fac-1.)*X_022; /* Old versions were missing the
        minus sign in this line, which introduced a very small error
        on the high-l C_l^TE lensed spectrum [credits for bug fix:
        Selim Hotinli] */

        /* X_242 = 0.25*sqrt4[l] * exp(-(fac-5./2.)*sigma2[index_mu]); */
        X_242 = 0.25*sqrt4[l] * X_000; /* Order 0 */
        if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {

          /* X_121 = - 0.5*sqrt2[l] * exp(-(fac-2./3.)*sigma2[index_mu]);
             X_132 = - 0.5*sqrt3[l] * exp(-(fac-5./3.)*sigma2[index_mu]); */
          X_121 = -0.5*sqrt2[l] * X_000 * (1+2./3.*sigma2[index_mu]); /* Order 1 */
          X_132 = -0.5*sqrt3[l] * X_000 * (1+5./3.*sigma2[index_mu]); /* Order 1 */
        }
      }


      if (ple->has_tt==_TRUE_) {

        res = fac1*cl_tt[l];

        lens = (X_000*X_000*d00[index_mu][l] +
                X_p000*X_p000*d1m1[index_mu][l]
                *Cgl2[index_mu]*8./(ll*(ll+1)) +
                (X_p000*X_p000*d00[index_mu][l] +
                 X_220*X_220*d2m2[index_mu][l])
                *Cgl2[index_mu]*Cgl2[index_mu]);
        if (ppr->accurate_lensing == _FALSE_) {
          /* Remove unlensed correlation function */
          lens -= d00[index_mu][l];
        }
        res *= lens;
        ksi[index_mu] += res;
      }

      if (ple->has_te==_TRUE_) {

        resX = fac1*cl_te[l];


        lens = ( X_022*X_000*d20[index_mu][l] +
                 Cgl2[index_mu]*2.*X_p000/sqrt5[l] *
                 (X_121*d11[index_mu][l] + X_132*d3m1[index_mu][l]) +
                 0.5 * Cgl2[index_mu] * Cgl2[index_mu] *
                 ( ( 2.*X_p022*X_p000+X_220*X_220 ) *
                   d20[index_mu][l] + X_220*X_242*d4m2[index_mu][l] ) );
        if (ppr->accurate_lensing == _FALSE_) {
          lens -= d20[index_mu][l];
        }
        resX *= lens;
        ksiX[index_mu] += resX;
      }

      if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {

        resp = fac1*(cl_ee[l]+cl_bb[l]);
        resm = fac1*(cl_ee[l]-cl_bb[l]);

        lensp = ( X_022*X_022*d22[index_mu][l] +
                  2.*Cgl2[index_mu]*X_132*X_121*d31[index_mu][l] +
                  Cgl2[index_mu]*Cgl2[index_mu] *
                  ( X_p022*X_p022*d22[index_mu][l] +
                    X_242*X_220*d40[index_mu][l] ) );

        lensm = ( X_022*X_022*d2m2[index_mu][l] +
                  Cgl2[index_mu] *
                  ( X_121*X_121*d1m1[index_mu][l] +
                    X_132*X_132*d3m3[index_mu][l] ) +
                  0.5 * Cgl2[index_mu] * Cgl2[index_mu] *
                  ( 2.*X_p022*X_p022*d2m2[index_mu][l] +
                    X_220*X_220*d00[index_mu][l] +
                    X_242*X_242*d4m4[index_mu][l] ) );
        if (ppr->accurate_lensing == _FALSE_) {
          lensp -= d22[index_mu][l];
          lensm -= d2m2[index_mu][l];
        }
        resp *= lensp;
        resm *= lensm;
        ksip[index_mu] += resp;
        ksim[index_mu] += resm;
      }
    }
  }
  //fin = omp_get_wtime();
  //cpu_time = (fin-debut);
  //printf("time in ksi=%4.3f s\n",cpu_time);


  /** - compute lensed \f$ C_l\f$'s by integration */
  //debut = omp_get_wtime();
  if (ple->has_tt==_TRUE_) {
    class_call(lensing_lensed_cl_tt(ksi,d00,w8,num_mu-1,ple),
               ple->error_message,
               ple->error_message);
    if (ppr->accurate_lensing == _FALSE_) {
      class_call(lensing_addback_cl_tt(ple,cl_tt),
                 ple->error_message,
                 ple->error_message);
    }
  }

  if (ple->has_te==_TRUE_) {
    class_call(lensing_lensed_cl_te(ksiX,d20,w8,num_mu-1,ple),
               ple->error_message,
               ple->error_message);
    if (ppr->accurate_lensing == _FALSE_) {
      class_call(lensing_addback_cl_te(ple,cl_te),
                 ple->error_message,
                 ple->error_message);
    }
  }

  if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {

    class_call(lensing_lensed_cl_ee_bb(ksip,ksim,d22,d2m2,w8,num_mu-1,ple),
               ple->error_message,
               ple->error_message);
    if (ppr->accurate_lensing == _FALSE_) {
      class_call(lensing_addback_cl_ee_bb(ple,cl_ee,cl_bb),
                 ple->error_message,
                 ple->error_message);
    }
  }
  //fin=omp_get_wtime();
  //cpu_time = (fin-debut);
  //printf("time in final lensing computation=%4.3f s\n",cpu_time);

  /** - spline computed \f$ C_l\f$'s in view of interpolation */

  class_call(array_spline_table_lines(ple->l,
                                      ple->l_size,
                                      ple->cl_lens,
                                      ple->lt_size,
                                      ple->ddcl_lens,
                                      _SPLINE_EST_DERIV_,
                                      ple->error_message),
             ple->error_message,
             ple->error_message);

  /** - Free lots of stuff **/
  free(buf_dxx);

  free(d00);
  free(d11);
  free(d1m1);
  free(d2m2);
  if (ple->has_te==_TRUE_) {
    free(d20);
    free(d3m1);
    free(d4m2);
  }
  if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
    free(d22);
    free(d31);
    free(d3m3);
    free(d40);
    free(d4m4);
  }

  if (ple->has_tt==_TRUE_)
    free(ksi);
  if (ple->has_te==_TRUE_)
    free(ksiX);
  if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
    free(ksip);
    free(ksim);
  }
  free(Cgl);
  free(Cgl2);
  free(sigma2);

  free(mu);
  free(w8);

  free(cl_unlensed);
  free(cl_tt);
  if (ple->has_te==_TRUE_)
    free(cl_te);
  if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
    free(cl_ee);
    free(cl_bb);
  }
  free(cl_pp);
  /** - Exit **/

  return _SUCCESS_;

}

/**
 * This routine frees all the memory space allocated by lensing_init().
 *
 * To be called at the end of each run, only when no further calls to
 * lensing_cl_at_l() are needed.
 *
 * @param ple Input: pointer to lensing structure (which fields must be freed)
 * @return the error status
 */

int lensing_free(
                 struct lensing * ple
                 ) {

  if (ple->has_lensed_cls == _TRUE_) {

    free(ple->l);
    free(ple->cl_lens);
    free(ple->ddcl_lens);
    free(ple->l_max_lt);

  }

  return _SUCCESS_;

}

/**
 * This routine defines indices and allocates tables in the lensing structure
 *
 * @param ppr  Input: pointer to precision structure
 * @param phr  Input: pointer to harmonic structure
 * @param ple  Input/output: pointer to lensing structure
 * @return the error status
 */

int lensing_indices(
                    struct precision * ppr,
                    struct harmonic * phr,
                    struct lensing * ple
                    ){

  int index_l;

  double ** cl_md_ic; /* array with argument
                         cl_md_ic[index_md][index_ic1_ic2*phr->ct_size+index_ct] */

  double ** cl_md;    /* array with argument
                         cl_md[index_md][index_ct] */

  int index_md;
  int index_lt;

  /* indices of all Cl types (lensed and unlensed) */

  if (phr->has_tt == _TRUE_) {
    ple->has_tt = _TRUE_;
    ple->index_lt_tt=phr->index_ct_tt;
  }
  else {
    ple->has_tt = _FALSE_;
  }

  if (phr->has_ee == _TRUE_) {
    ple->has_ee = _TRUE_;
    ple->index_lt_ee=phr->index_ct_ee;
  }
  else {
    ple->has_ee = _FALSE_;
  }

  if (phr->has_te == _TRUE_) {
    ple->has_te = _TRUE_;
    ple->index_lt_te=phr->index_ct_te;
  }
  else {
    ple->has_te = _FALSE_;
  }

  if (phr->has_bb == _TRUE_) {
    ple->has_bb = _TRUE_;
    ple->index_lt_bb=phr->index_ct_bb;
  }
  else {
    ple->has_bb = _FALSE_;
  }

  if (phr->has_pp == _TRUE_) {
    ple->has_pp = _TRUE_;
    ple->index_lt_pp=phr->index_ct_pp;
  }
  else {
    ple->has_pp = _FALSE_;
  }

  if (phr->has_tp == _TRUE_) {
    ple->has_tp = _TRUE_;
    ple->index_lt_tp=phr->index_ct_tp;
  }
  else {
    ple->has_tp = _FALSE_;
  }

  if (phr->has_dd == _TRUE_) {
    ple->has_dd = _TRUE_;
    ple->index_lt_dd=phr->index_ct_dd;
  }
  else {
    ple->has_dd = _FALSE_;
  }

  if (phr->has_td == _TRUE_) {
    ple->has_td = _TRUE_;
    ple->index_lt_td=phr->index_ct_td;
  }
  else {
    ple->has_td = _FALSE_;
  }

  if (phr->has_ll == _TRUE_) {
    ple->has_ll = _TRUE_;
    ple->index_lt_ll=phr->index_ct_ll;
  }
  else {
    ple->has_ll = _FALSE_;
  }

  if (phr->has_tl == _TRUE_) {
    ple->has_tl = _TRUE_;
    ple->index_lt_tl=phr->index_ct_tl;
  }
  else {
    ple->has_tl = _FALSE_;
  }

  ple->lt_size = phr->ct_size;

  /* number of multipoles */

  ple->l_unlensed_max = phr->l_max_tot;

  ple->l_lensed_max = ple->l_unlensed_max - ppr->delta_l_max;

  for (index_l=0; (index_l < phr->l_size_max) && (phr->l[index_l] <= ple->l_lensed_max); index_l++);

  if (index_l < phr->l_size_max) index_l++; /* one more point in order to be able to interpolate till ple->l_lensed_max */

  ple->l_size = index_l+1;

  class_alloc(ple->l,ple->l_size*sizeof(double),ple->error_message);

  for (index_l=0; index_l < ple->l_size; index_l++) {

    ple->l[index_l] = phr->l[index_l];

  }

  /* allocate table where results will be stored */

  class_alloc(ple->cl_lens,
              ple->l_size*ple->lt_size*sizeof(double),
              ple->error_message);

  class_alloc(ple->ddcl_lens,
              ple->l_size*ple->lt_size*sizeof(double),
              ple->error_message);

  /* fill with unlensed cls */

  class_alloc(cl_md_ic,
              phr->md_size*sizeof(double *),
              ple->error_message);

  class_alloc(cl_md,
              phr->md_size*sizeof(double *),
              ple->error_message);

  for (index_md = 0; index_md < phr->md_size; index_md++) {

    if (phr->md_size > 1)

      class_alloc(cl_md[index_md],
                  phr->ct_size*sizeof(double),
                  ple->error_message);

    if (phr->ic_size[index_md] > 1)

      class_alloc(cl_md_ic[index_md],
                  phr->ic_ic_size[index_md]*phr->ct_size*sizeof(double),
                  ple->error_message);
  }

  for (index_l=0; index_l<ple->l_size; index_l++) {

    class_call(harmonic_cl_at_l(phr,ple->l[index_l],&(ple->cl_lens[index_l*ple->lt_size]),cl_md,cl_md_ic),
               phr->error_message,
               ple->error_message);

  }

  for (index_md = 0; index_md < phr->md_size; index_md++) {

    if (phr->md_size > 1)
      free(cl_md[index_md]);

    if (phr->ic_size[index_md] > 1)
      free(cl_md_ic[index_md]);

  }

  free(cl_md_ic);
  free(cl_md);

  /* we want to output Cl_lensed up to the same l_max as Cl_unlensed
     (even if a number delta_l_max of extra values of l have been used
     internally for more accurate results). Notable exception to the
     above rule: ClBB_lensed(scalars) must be outputed at least up to the same l_max as
     ClEE_unlensed(scalars) (since ClBB_unlensed is null for scalars)
  */

  class_alloc(ple->l_max_lt,ple->lt_size*sizeof(double),ple->error_message);
  for (index_lt = 0; index_lt < ple->lt_size; index_lt++) {
    ple->l_max_lt[index_lt]=0.;
    for (index_md = 0; index_md < phr->md_size; index_md++) {
      ple->l_max_lt[index_lt]=MAX(ple->l_max_lt[index_lt],phr->l_max_ct[index_md][index_lt]);

      if ((ple->has_bb == _TRUE_) && (ple->has_ee == _TRUE_) && (index_lt == ple->index_lt_bb)) {
        ple->l_max_lt[index_lt]=MAX(ple->l_max_lt[index_lt],phr->l_max_ct[index_md][ple->index_lt_ee]);
      }

    }
  }

  return _SUCCESS_;

}

/**
 * This routine computes the lensed power spectra by Gaussian quadrature
 *
 * @param ksi  Input: Lensed correlation function (ksi[index_mu])
 * @param d00  Input: Legendre polynomials (\f$ d^l_{00}\f$[l][index_mu])
 * @param w8   Input: Legendre quadrature weights (w8[index_mu])
 * @param nmu  Input: Number of quadrature points (0<=index_mu<=nmu)
 * @param ple  Input/output: Pointer to the lensing structure
 * @return the error status
 */


int lensing_lensed_cl_tt(
                         double *ksi,
                         double **d00,
                         double *w8,
                         int nmu,
                         struct lensing * ple
                         ) {

  double cle;
  int imu;
  int index_l;

  /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
  private (imu,index_l,cle)                     \
  schedule (static)

  for(index_l=0; index_l<ple->l_size; index_l++){
    cle=0;
    for (imu=0;imu<nmu;imu++) {
      cle += ksi[imu]*d00[imu][(int)ple->l[index_l]]*w8[imu]; /* loop could be optimized */
    }
    ple->cl_lens[index_l*ple->lt_size+ple->index_lt_tt]=cle*2.0*_PI_;
  }

  return _SUCCESS_;
}

/**
 * This routine adds back the unlensed \f$ cl_{tt}\f$ power spectrum
 * Used in case of fast (and BB inaccurate) integration of
 * correlation functions.
 *
 * @param ple   Input/output: Pointer to the lensing structure
 * @param cl_tt Input: Array of unlensed power spectrum
 * @return the error status
 */

int lensing_addback_cl_tt(
                          struct lensing * ple,
                          double *cl_tt) {
  int index_l, l;

  for (index_l=0; index_l<ple->l_size; index_l++) {
    l = (int)ple->l[index_l];
    ple->cl_lens[index_l*ple->lt_size+ple->index_lt_tt] += cl_tt[l];
  }
  return _SUCCESS_;

}

/**
 * This routine computes the lensed power spectra by Gaussian quadrature
 *
 * @param ksiX Input: Lensed correlation function (ksiX[index_mu])
 * @param d20  Input: Wigner d-function (\f$ d^l_{20}\f$[l][index_mu])
 * @param w8   Input: Legendre quadrature weights (w8[index_mu])
 * @param nmu  Input: Number of quadrature points (0<=index_mu<=nmu)
 * @param ple  Input/output: Pointer to the lensing structure
 * @return the error status
 */


int lensing_lensed_cl_te(
                         double *ksiX,
                         double **d20,
                         double *w8,
                         int nmu,
                         struct lensing * ple
                         ) {

  double clte;
  int imu;
  int index_l;

  /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
  private (imu,index_l,clte)                    \
  schedule (static)

  for(index_l=0; index_l < ple->l_size; index_l++){
    clte=0;
    for (imu=0;imu<nmu;imu++) {
      clte += ksiX[imu]*d20[imu][(int)ple->l[index_l]]*w8[imu]; /* loop could be optimized */
    }
    ple->cl_lens[index_l*ple->lt_size+ple->index_lt_te]=clte*2.0*_PI_;
  }

  return _SUCCESS_;
}

/**
 * This routine adds back the unlensed \f$ cl_{te}\f$ power spectrum
 * Used in case of fast (and BB inaccurate) integration of
 * correlation functions.
 *
 * @param ple   Input/output: Pointer to the lensing structure
 * @param cl_te Input: Array of unlensed power spectrum
 * @return the error status
 */

int lensing_addback_cl_te(
                          struct lensing * ple,
                          double *cl_te) {
  int index_l, l;

  for (index_l=0; index_l<ple->l_size; index_l++) {
    l = (int)ple->l[index_l];
    ple->cl_lens[index_l*ple->lt_size+ple->index_lt_te] += cl_te[l];
  }
  return _SUCCESS_;

}

/**
 * This routine computes the lensed power spectra by Gaussian quadrature
 *
 * @param ksip Input: Lensed correlation function (ksi+[index_mu])
 * @param ksim Input: Lensed correlation function (ksi-[index_mu])
 * @param d22  Input: Wigner d-function (\f$ d^l_{22}\f$[l][index_mu])
 * @param d2m2 Input: Wigner d-function (\f$ d^l_{2-2}\f$[l][index_mu])
 * @param w8   Input: Legendre quadrature weights (w8[index_mu])
 * @param nmu  Input: Number of quadrature points (0<=index_mu<=nmu)
 * @param ple  Input/output: Pointer to the lensing structure
 * @return the error status
 */


int lensing_lensed_cl_ee_bb(
                            double *ksip,
                            double *ksim,
                            double **d22,
                            double **d2m2,
                            double *w8,
                            int nmu,
                            struct lensing * ple
                            ) {

  double clp, clm;
  int imu;
  int index_l;

  /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
  private (imu,index_l,clp,clm)                 \
  schedule (static)

  for(index_l=0; index_l < ple->l_size; index_l++){
    clp=0; clm=0;
    for (imu=0;imu<nmu;imu++) {
      clp += ksip[imu]*d22[imu][(int)ple->l[index_l]]*w8[imu]; /* loop could be optimized */
      clm += ksim[imu]*d2m2[imu][(int)ple->l[index_l]]*w8[imu]; /* loop could be optimized */
    }
    ple->cl_lens[index_l*ple->lt_size+ple->index_lt_ee]=(clp+clm)*_PI_;
    ple->cl_lens[index_l*ple->lt_size+ple->index_lt_bb]=(clp-clm)*_PI_;
  }

  return _SUCCESS_;
}

/**
 * This routine adds back the unlensed \f$ cl_{ee}\f$, \f$ cl_{bb}\f$ power spectra
 * Used in case of fast (and BB inaccurate) integration of
 * correlation functions.
 *
 * @param ple   Input/output: Pointer to the lensing structure
 * @param cl_ee Input: Array of unlensed power spectrum
 * @param cl_bb Input: Array of unlensed power spectrum
 * @return the error status
 */

int lensing_addback_cl_ee_bb(
                             struct lensing * ple,
                             double * cl_ee,
                             double * cl_bb) {

  int index_l, l;

  for (index_l=0; index_l<ple->l_size; index_l++) {
    l = (int)ple->l[index_l];
    ple->cl_lens[index_l*ple->lt_size+ple->index_lt_ee] += cl_ee[l];
    ple->cl_lens[index_l*ple->lt_size+ple->index_lt_bb] += cl_bb[l];
  }
  return _SUCCESS_;

}

/**
 * This routine computes the d00 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d00    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d00(
                double * mu,
                int num_mu,
                int lmax,
                double ** d00
                ) {
  double ll, dlm1, dl, dlp1;
  int index_mu, l;
  double *fac1, *fac2, *fac3;
  ErrorMsg erreur;

  class_alloc(fac1,lmax*sizeof(double),erreur);
  class_alloc(fac2,lmax*sizeof(double),erreur);
  class_alloc(fac3,lmax*sizeof(double),erreur);
  for (l=1; l<lmax; l++) {
    ll = (double) l;
    fac1[l] = sqrt((2*ll+3)/(2*ll+1))*(2*ll+1)/(ll+1);
    fac2[l] = sqrt((2*ll+3)/(2*ll-1))*ll/(ll+1);
    fac3[l] = sqrt(2./(2*ll+3));
  }

#pragma omp parallel for                        \
  private (index_mu,dlm1,dl,dlp1,l,ll)          \
  schedule (static)

  for (index_mu=0;index_mu<num_mu;index_mu++) {
    dlm1=1.0/sqrt(2.); /* l=0 */
    d00[index_mu][0]=dlm1*sqrt(2.);
    dl=mu[index_mu] * sqrt(3./2.); /*l=1*/
    d00[index_mu][1]=dl*sqrt(2./3.);
    for(l=1;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d00 recurrence, supposed to be more stable */
      dlp1 = fac1[l]*mu[index_mu]*dl - fac2[l]*dlm1;
      d00[index_mu][l+1] = dlp1 * fac3[l];
      dlm1 = dl;
      dl = dlp1;
    }
  }
  free(fac1); free(fac2); free(fac3);
  return _SUCCESS_;
}


/**
 * This routine computes the d11 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d11    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d11(
                double * mu,
                int num_mu,
                int lmax,
                double ** d11
                ) {
  double ll, dlm1, dl, dlp1;
  int index_mu, l;
  double *fac1, *fac2, *fac3, *fac4;
  ErrorMsg erreur;
  class_alloc(fac1,lmax*sizeof(double),erreur);
  class_alloc(fac2,lmax*sizeof(double),erreur);
  class_alloc(fac3,lmax*sizeof(double),erreur);
  class_alloc(fac4,lmax*sizeof(double),erreur);
  for (l=2;l<lmax;l++) {
    ll = (double) l;
    fac1[l] = sqrt((2*ll+3)/(2*ll+1))*(ll+1)*(2*ll+1)/(ll*(ll+2));
    fac2[l] = 1.0/(ll*(ll+1.));
    fac3[l] = sqrt((2*ll+3)/(2*ll-1))*(ll-1)*(ll+1)/(ll*(ll+2))*(ll+1)/ll;
    fac4[l] = sqrt(2./(2*ll+3));
  }
#pragma omp parallel for                        \
  private (index_mu,dlm1,dl,dlp1,l,ll)          \
  schedule (static)

  for (index_mu=0;index_mu<num_mu;index_mu++) {
    d11[index_mu][0]=0;
    dlm1=(1.0+mu[index_mu])/2. * sqrt(3./2.); /*l=1*/
    d11[index_mu][1]=dlm1 * sqrt(2./3.);
    dl=(1.0+mu[index_mu])/2.*(2.0*mu[index_mu]-1.0) * sqrt(5./2.); /*l=2*/
    d11[index_mu][2] = dl * sqrt(2./5.);
    for(l=2;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d11 recurrence, supposed to be more stable */
      dlp1 = fac1[l]*(mu[index_mu]-fac2[l])*dl - fac3[l]*dlm1;
      d11[index_mu][l+1] = dlp1 * fac4[l];
      dlm1 = dl;
      dl = dlp1;
    }
  }
  free(fac1); free(fac2); free(fac3); free(fac4);
  return _SUCCESS_;
}

/**
 * This routine computes the d1m1 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d1m1    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d1m1(
                 double * mu,
                 int num_mu,
                 int lmax,
                 double ** d1m1
                 ) {
  double ll, dlm1, dl, dlp1;
  int index_mu, l;
  double *fac1, *fac2, *fac3, *fac4;
  ErrorMsg erreur;
  class_alloc(fac1,lmax*sizeof(double),erreur);
  class_alloc(fac2,lmax*sizeof(double),erreur);
  class_alloc(fac3,lmax*sizeof(double),erreur);
  class_alloc(fac4,lmax*sizeof(double),erreur);
  for (l=2;l<lmax;l++) {
    ll = (double) l;
    fac1[l] = sqrt((2*ll+3)/(2*ll+1))*(ll+1)*(2*ll+1)/(ll*(ll+2));
    fac2[l] = 1.0/(ll*(ll+1.));
    fac3[l] = sqrt((2*ll+3)/(2*ll-1))*(ll-1)*(ll+1)/(ll*(ll+2))*(ll+1)/ll;
    fac4[l] = sqrt(2./(2*ll+3));
  }
#pragma omp parallel for                        \
  private (index_mu,dlm1,dl,dlp1,l,ll)          \
  schedule (static)

  for (index_mu=0;index_mu<num_mu;index_mu++) {
    d1m1[index_mu][0]=0;
    dlm1=(1.0-mu[index_mu])/2. * sqrt(3./2.); /*l=1*/
    d1m1[index_mu][1]=dlm1 * sqrt(2./3.);
    dl=(1.0-mu[index_mu])/2.*(2.0*mu[index_mu]+1.0) * sqrt(5./2.); /*l=2*/
    d1m1[index_mu][2] = dl * sqrt(2./5.);
    for(l=2;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d1m1 recurrence, supposed to be more stable */
      dlp1 = fac1[l]*(mu[index_mu]+fac2[l])*dl - fac3[l]*dlm1;
      d1m1[index_mu][l+1] = dlp1 * fac4[l];
      dlm1 = dl;
      dl = dlp1;
    }
  }
  free(fac1); free(fac2); free(fac3); free(fac4);
  return _SUCCESS_;
}

/**
 * This routine computes the d2m2 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d2m2   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d2m2(
                 double * mu,
                 int num_mu,
                 int lmax,
                 double ** d2m2
                 ) {
  double ll, dlm1, dl, dlp1;
  int index_mu, l;
  double *fac1, *fac2, *fac3, *fac4;
  ErrorMsg erreur;
  class_alloc(fac1,lmax*sizeof(double),erreur);
  class_alloc(fac2,lmax*sizeof(double),erreur);
  class_alloc(fac3,lmax*sizeof(double),erreur);
  class_alloc(fac4,lmax*sizeof(double),erreur);
  for (l=2;l<lmax;l++) {
    ll = (double) l;
    fac1[l] = sqrt((2*ll+3)/(2*ll+1))*(ll+1)*(2*ll+1)/((ll-1)*(ll+3));
    fac2[l] = 4.0/(ll*(ll+1));
    fac3[l] = sqrt((2*ll+3)/(2*ll-1))*(ll-2)*(ll+2)/((ll-1)*(ll+3))*(ll+1)/ll;
    fac4[l] = sqrt(2./(2*ll+3));
  }
#pragma omp parallel for                        \
  private (index_mu,dlm1,dl,dlp1,l,ll)          \
  schedule (static)

  for (index_mu=0;index_mu<num_mu;index_mu++) {
    d2m2[index_mu][0]=0;
    dlm1=0.; /*l=1*/
    d2m2[index_mu][1]=0;
    dl=(1.0-mu[index_mu])*(1.0-mu[index_mu])/4. * sqrt(5./2.); /*l=2*/
    d2m2[index_mu][2] = dl * sqrt(2./5.);
    for(l=2;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d2m2 recurrence, supposed to be more stable */
      dlp1 = fac1[l]*(mu[index_mu]+fac2[l])*dl - fac3[l]*dlm1;
      d2m2[index_mu][l+1] = dlp1 * fac4[l];
      dlm1 = dl;
      dl = dlp1;
    }
  }
  free(fac1); free(fac2); free(fac3); free(fac4);
  return _SUCCESS_;
}

/**
 * This routine computes the d22 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d22    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d22(
                double * mu,
                int num_mu,
                int lmax,
                double ** d22
                ) {
  double ll, dlm1, dl, dlp1;
  int index_mu, l;
  double *fac1, *fac2, *fac3, *fac4;
  ErrorMsg erreur;
  class_alloc(fac1,lmax*sizeof(double),erreur);
  class_alloc(fac2,lmax*sizeof(double),erreur);
  class_alloc(fac3,lmax*sizeof(double),erreur);
  class_alloc(fac4,lmax*sizeof(double),erreur);
  for (l=2;l<lmax;l++) {
    ll = (double) l;
    fac1[l] = sqrt((2*ll+3)/(2*ll+1))*(ll+1)*(2*ll+1)/((ll-1)*(ll+3));
    fac2[l] = 4.0/(ll*(ll+1));
    fac3[l] = sqrt((2*ll+3)/(2*ll-1))*(ll-2)*(ll+2)/((ll-1)*(ll+3))*(ll+1)/ll;
    fac4[l] = sqrt(2./(2*ll+3));
  }
#pragma omp parallel for                        \
  private (index_mu,dlm1,dl,dlp1,l,ll)          \
  schedule (static)

  for (index_mu=0;index_mu<num_mu;index_mu++) {
    d22[index_mu][0]=0;
    dlm1=0.; /*l=1*/
    d22[index_mu][1]=0;
    dl=(1.0+mu[index_mu])*(1.0+mu[index_mu])/4. * sqrt(5./2.); /*l=2*/
    d22[index_mu][2] = dl * sqrt(2./5.);
    for(l=2;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
      dlp1 = fac1[l]*(mu[index_mu]-fac2[l])*dl - fac3[l]*dlm1;
      d22[index_mu][l+1] = dlp1 * fac4[l];
      dlm1 = dl;
      dl = dlp1;
    }
  }
  free(fac1); free(fac2); free(fac3); free(fac4);
  return _SUCCESS_;
}

/**
 * This routine computes the d20 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d20    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d20(
                double * mu,
                int num_mu,
                int lmax,
                double ** d20
                ) {
  double ll, dlm1, dl, dlp1;
  int index_mu, l;
  double *fac1, *fac3, *fac4;
  ErrorMsg erreur;
  class_alloc(fac1,lmax*sizeof(double),erreur);
  class_alloc(fac3,lmax*sizeof(double),erreur);
  class_alloc(fac4,lmax*sizeof(double),erreur);
  for (l=2;l<lmax;l++) {
    ll = (double) l;
    fac1[l] = sqrt((2*ll+3)*(2*ll+1)/((ll-1)*(ll+3)));
    fac3[l] = sqrt((2*ll+3)*(ll-2)*(ll+2)/((2*ll-1)*(ll-1)*(ll+3)));
    fac4[l] = sqrt(2./(2*ll+3));
  }
#pragma omp parallel for                        \
  private (index_mu,dlm1,dl,dlp1,l,ll)          \
  schedule (static)

  for (index_mu=0;index_mu<num_mu;index_mu++) {
    d20[index_mu][0]=0;
    dlm1=0.; /*l=1*/
    d20[index_mu][1]=0;
    dl=sqrt(15.)/4.*(1-mu[index_mu]*mu[index_mu]); /*l=2*/
    d20[index_mu][2] = dl * sqrt(2./5.);
    for(l=2;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
      dlp1 = fac1[l]*mu[index_mu]*dl - fac3[l]*dlm1;
      d20[index_mu][l+1] = dlp1 * fac4[l];
      dlm1 = dl;
      dl = dlp1;
    }
  }
  free(fac1); free(fac3); free(fac4);
  return _SUCCESS_;
}

/**
 * This routine computes the d31 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d31    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d31(
                double * mu,
                int num_mu,
                int lmax,
                double ** d31
                ) {
  double ll, dlm1, dl, dlp1;
  int index_mu, l;
  double *fac1, *fac2, *fac3, *fac4;
  ErrorMsg erreur;
  class_alloc(fac1,lmax*sizeof(double),erreur);
  class_alloc(fac2,lmax*sizeof(double),erreur);
  class_alloc(fac3,lmax*sizeof(double),erreur);
  class_alloc(fac4,lmax*sizeof(double),erreur);
  for (l=3;l<lmax;l++) {
    ll = (double) l;
    fac1[l] = sqrt((2*ll+3)*(2*ll+1)/((ll-2)*(ll+4)*ll*(ll+2))) * (ll+1);
    fac2[l] = 3.0/(ll*(ll+1));
    fac3[l] = sqrt((2*ll+3)/(2*ll-1)*(ll-3)*(ll+3)*(ll-1)*(ll+1)/((ll-2)*(ll+4)*ll*(ll+2)))*(ll+1)/ll;
    fac4[l] = sqrt(2./(2*ll+3));
  }
#pragma omp parallel for                        \
  private (index_mu,dlm1,dl,dlp1,l,ll)          \
  schedule (static)

  for (index_mu=0;index_mu<num_mu;index_mu++) {
    d31[index_mu][0]=0;
    d31[index_mu][1]=0;
    dlm1=0.; /*l=2*/
    d31[index_mu][2]=0;
    dl=sqrt(105./2.)*(1+mu[index_mu])*(1+mu[index_mu])*(1-mu[index_mu])/8.; /*l=3*/
    d31[index_mu][3] = dl * sqrt(2./7.);
    for(l=3;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
      dlp1 = fac1[l]*(mu[index_mu]-fac2[l])*dl - fac3[l]*dlm1;
      d31[index_mu][l+1] = dlp1 * fac4[l];
      dlm1 = dl;
      dl = dlp1;
    }
  }
  free(fac1); free(fac2); free(fac3); free(fac4);
  return _SUCCESS_;
}

/**
 * This routine computes the d3m1 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d3m1   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d3m1(
                 double * mu,
                 int num_mu,
                 int lmax,
                 double ** d3m1
                 ) {
  double ll, dlm1, dl, dlp1;
  int index_mu, l;
  double *fac1, *fac2, *fac3, *fac4;
  ErrorMsg erreur;
  class_alloc(fac1,lmax*sizeof(double),erreur);
  class_alloc(fac2,lmax*sizeof(double),erreur);
  class_alloc(fac3,lmax*sizeof(double),erreur);
  class_alloc(fac4,lmax*sizeof(double),erreur);
  for (l=3;l<lmax;l++) {
    ll = (double) l;
    fac1[l] = sqrt((2*ll+3)*(2*ll+1)/((ll-2)*(ll+4)*ll*(ll+2))) * (ll+1);
    fac2[l] = 3.0/(ll*(ll+1));
    fac3[l] = sqrt((2*ll+3)/(2*ll-1)*(ll-3)*(ll+3)*(ll-1)*(ll+1)/((ll-2)*(ll+4)*ll*(ll+2)))*(ll+1)/ll;
    fac4[l] = sqrt(2./(2*ll+3));
  }
#pragma omp parallel for                        \
  private (index_mu,dlm1,dl,dlp1,l,ll)          \
  schedule (static)

  for (index_mu=0;index_mu<num_mu;index_mu++) {
    d3m1[index_mu][0]=0;
    d3m1[index_mu][1]=0;
    dlm1=0.; /*l=2*/
    d3m1[index_mu][2]=0;
    dl=sqrt(105./2.)*(1+mu[index_mu])*(1-mu[index_mu])*(1-mu[index_mu])/8.; /*l=3*/
    d3m1[index_mu][3] = dl * sqrt(2./7.);
    for(l=3;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
      dlp1 = fac1[l]*(mu[index_mu]+fac2[l])*dl - fac3[l]*dlm1;
      d3m1[index_mu][l+1] = dlp1 * fac4[l];
      dlm1 = dl;
      dl = dlp1;
    }
  }
  free(fac1); free(fac2); free(fac3); free(fac4);
  return _SUCCESS_;
}

/**
 * This routine computes the d3m3 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d3m3   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d3m3(
                 double * mu,
                 int num_mu,
                 int lmax,
                 double ** d3m3
                 ) {
  double ll, dlm1, dl, dlp1;
  int index_mu, l;
  double *fac1, *fac2, *fac3, *fac4;
  ErrorMsg erreur;
  class_alloc(fac1,lmax*sizeof(double),erreur);
  class_alloc(fac2,lmax*sizeof(double),erreur);
  class_alloc(fac3,lmax*sizeof(double),erreur);
  class_alloc(fac4,lmax*sizeof(double),erreur);
  for (l=3;l<lmax;l++) {
    ll = (double) l;
    fac1[l] = sqrt((2*ll+3)*(2*ll+1))*(ll+1)/((ll-2)*(ll+4));
    fac2[l] = 9.0/(ll*(ll+1));
    fac3[l] = sqrt((2*ll+3)/(2*ll-1))*(ll-3)*(ll+3)*(l+1)/((ll-2)*(ll+4)*ll);
    fac4[l] = sqrt(2./(2*ll+3));
  }
#pragma omp parallel for                        \
  private (index_mu,dlm1,dl,dlp1,l,ll)          \
  schedule (static)

  for (index_mu=0;index_mu<num_mu;index_mu++) {
    d3m3[index_mu][0]=0;
    d3m3[index_mu][1]=0;
    dlm1=0.; /*l=2*/
    d3m3[index_mu][2]=0;
    dl=sqrt(7./2.)*(1-mu[index_mu])*(1-mu[index_mu])*(1-mu[index_mu])/8.; /*l=3*/
    d3m3[index_mu][3] = dl * sqrt(2./7.);
    for(l=3;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
      dlp1 = fac1[l]*(mu[index_mu]+fac2[l])*dl - fac3[l]*dlm1;
      d3m3[index_mu][l+1] = dlp1 * fac4[l];
      dlm1 = dl;
      dl = dlp1;
    }
  }
  free(fac1); free(fac2); free(fac3); free(fac4);
  return _SUCCESS_;
}

/**
 * This routine computes the d40 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d40    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d40(
                double * mu,
                int num_mu,
                int lmax,
                double ** d40
                ) {
  double ll, dlm1, dl, dlp1;
  int index_mu, l;
  double *fac1, *fac3, *fac4;
  ErrorMsg erreur;
  class_alloc(fac1,lmax*sizeof(double),erreur);
  class_alloc(fac3,lmax*sizeof(double),erreur);
  class_alloc(fac4,lmax*sizeof(double),erreur);
  for (l=4;l<lmax;l++) {
    ll = (double) l;
    fac1[l] = sqrt((2*ll+3)*(2*ll+1)/((ll-3)*(ll+5)));
    fac3[l] = sqrt((2*ll+3)*(ll-4)*(ll+4)/((2*ll-1)*(ll-3)*(ll+5)));
    fac4[l] = sqrt(2./(2*ll+3));
  }
#pragma omp parallel for                        \
  private (index_mu,dlm1,dl,dlp1,l,ll)          \
  schedule (static)

  for (index_mu=0;index_mu<num_mu;index_mu++) {
    d40[index_mu][0]=0;
    d40[index_mu][1]=0;
    d40[index_mu][2]=0;
    dlm1=0.; /*l=3*/
    d40[index_mu][3]=0;
    dl=sqrt(315.)*(1+mu[index_mu])*(1+mu[index_mu])*(1-mu[index_mu])*(1-mu[index_mu])/16.; /*l=4*/
    d40[index_mu][4] = dl * sqrt(2./9.);
    for(l=4;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
      dlp1 = fac1[l]*mu[index_mu]*dl - fac3[l]*dlm1;
      d40[index_mu][l+1] = dlp1 * fac4[l];
      dlm1 = dl;
      dl = dlp1;
    }
  }
  free(fac1); free(fac3); free(fac4);
  return _SUCCESS_;
}

/**
 * This routine computes the d4m2 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d4m2   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d4m2(
                 double * mu,
                 int num_mu,
                 int lmax,
                 double ** d4m2
                 ) {
  double ll, dlm1, dl, dlp1;
  int index_mu, l;
  double *fac1, *fac2, *fac3, *fac4;
  ErrorMsg erreur;
  class_alloc(fac1,lmax*sizeof(double),erreur);
  class_alloc(fac2,lmax*sizeof(double),erreur);
  class_alloc(fac3,lmax*sizeof(double),erreur);
  class_alloc(fac4,lmax*sizeof(double),erreur);
  for (l=4;l<lmax;l++) {
    ll = (double) l;
    fac1[l] = sqrt((2*ll+3)*(2*ll+1)/((ll-3)*(ll+5)*(ll-1)*(ll+3))) * (ll+1.);
    fac2[l] = 8./(ll*(ll+1));
    fac3[l] = sqrt((2*ll+3)*(ll-4)*(ll+4)*(ll-2)*(ll+2)/((2*ll-1)*(ll-3)*(ll+5)*(ll-1)*(ll+3)))*(ll+1)/ll;
    fac4[l] = sqrt(2./(2*ll+3));
  }
#pragma omp parallel for                        \
  private (index_mu,dlm1,dl,dlp1,l,ll)          \
  schedule (static)

  for (index_mu=0;index_mu<num_mu;index_mu++) {
    d4m2[index_mu][0]=0;
    d4m2[index_mu][1]=0;
    d4m2[index_mu][2]=0;
    dlm1=0.; /*l=3*/
    d4m2[index_mu][3]=0;
    dl=sqrt(126.)*(1+mu[index_mu])*(1-mu[index_mu])*(1-mu[index_mu])*(1-mu[index_mu])/16.; /*l=4*/
    d4m2[index_mu][4] = dl * sqrt(2./9.);
    for(l=4;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
      dlp1 = fac1[l]*(mu[index_mu]+fac2[l])*dl - fac3[l]*dlm1;
      d4m2[index_mu][l+1] = dlp1 * fac4[l];
      dlm1 = dl;
      dl = dlp1;
    }
  }
  free(fac1); free(fac2); free(fac3); free(fac4);
  return _SUCCESS_;
}

/**
 * This routine computes the d4m4 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d4m4   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d4m4(
                 double * mu,
                 int num_mu,
                 int lmax,
                 double ** d4m4
                 ) {
  double ll, dlm1, dl, dlp1;
  int index_mu, l;
  double *fac1, *fac2, *fac3, *fac4;
  ErrorMsg erreur;
  class_alloc(fac1,lmax*sizeof(double),erreur);
  class_alloc(fac2,lmax*sizeof(double),erreur);
  class_alloc(fac3,lmax*sizeof(double),erreur);
  class_alloc(fac4,lmax*sizeof(double),erreur);
  for (l=4;l<lmax;l++) {
    ll = (double) l;
    fac1[l] = sqrt((2*ll+3)*(2*ll+1))*(ll+1)/((ll-3)*(ll+5));
    fac2[l] = 16./(ll*(ll+1));
    fac3[l] = sqrt((2*ll+3)/(2*ll-1))*(ll-4)*(ll+4)*(ll+1)/((ll-3)*(ll+5)*ll);
    fac4[l] = sqrt(2./(2*ll+3));
  }
#pragma omp parallel for                        \
  private (index_mu,dlm1,dl,dlp1,l,ll)          \
  schedule (static)

  for (index_mu=0;index_mu<num_mu;index_mu++) {
    d4m4[index_mu][0]=0;
    d4m4[index_mu][1]=0;
    d4m4[index_mu][2]=0;
    dlm1=0.; /*l=3*/
    d4m4[index_mu][3]=0;
    dl=sqrt(9./2.)*(1-mu[index_mu])*(1-mu[index_mu])*(1-mu[index_mu])*(1-mu[index_mu])/16.; /*l=4*/
    d4m4[index_mu][4] = dl * sqrt(2./9.);
    for(l=4;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
      dlp1 = fac1[l]*(mu[index_mu]+fac2[l])*dl - fac3[l]*dlm1;
      d4m4[index_mu][l+1] = dlp1 * fac4[l];
      dlm1 = dl;
      dl = dlp1;
    }
  }
  free(fac1); free(fac2); free(fac3); free(fac4);
  return _SUCCESS_;
}
