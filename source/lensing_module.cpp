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
 * -# lensing_init() at the beginning (but after spectra_init())
 * -# lensing_cl_at_l() at any time for computing Cl_lensed at any l
 * -# lensing_free() at the end
 */

#include "lensing_module.h"
#include "exceptions.h"
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

LensingModule::LensingModule(const Input& input, const SpectraModule& spectra_module)
: BaseModule(input)
, spectra_module_(spectra_module) {
  ThrowInvalidArgumentIf(lensing_init() != _SUCCESS_, error_message_);
}

LensingModule::~LensingModule() {
  lensing_free();
}

int LensingModule::lensing_cl_at_l(int l, double * cl_lensed) const {
  int last_index;
  int index_lt;

  class_test(l > l_lensed_max_,
             error_message_,
             "you asked for lensed Cls at l=%d, they were computed only up to l=%d, you should increase l_max_scalars or decrease the precision parameter delta_l_max",
             l,
             l_lensed_max_);

  class_call(array_interpolate_spline(l_,
                                      l_size_,
                                      cl_lens_,
                                      ddcl_lens_,
                                      lt_size_,
                                      l,
                                      &last_index,
                                      cl_lensed,
                                      lt_size_,
                                      error_message_),
             error_message_,
             error_message_);

  /* set to zero for the types such that l<l_max */
  for (index_lt=0; index_lt<lt_size_; index_lt++)
    if ((int)l > l_max_lt_[index_lt])
      cl_lensed[index_lt]=0.;

  return _SUCCESS_;
}

/**
 * This routine initializes the lensing structure (in particular,
 * computes table of lensed anisotropy spectra \f$ C_l^{X} \f$)
 *
 * @param ppr Input: pointer to precision structure
 * @param ppt Input: pointer to perturbation structure (just in case, not used in current version...)
 * @param psp Input: pointer to spectra structure
 * @param pnl Input: pointer to nonlinear structure
 * @param ple Output: pointer to initialized lensing structure
 * @return the error status
 */

int LensingModule::lensing_init() {

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
  double * cl_tt; /* unlensed  cl, to be filled to avoid repeated calls to spectra_cl_at_l */
  double * cl_te = NULL; /* unlensed  cl, to be filled to avoid repeated calls to spectra_cl_at_l */
  double * cl_ee = NULL; /* unlensed  cl, to be filled to avoid repeated calls to spectra_cl_at_l */
  double * cl_bb = NULL; /* unlensed  cl, to be filled to avoid repeated calls to spectra_cl_at_l */
  double * cl_pp; /* potential cl, to be filled to avoid repeated calls to spectra_cl_at_l */

  double res,resX,lens;
  double resp, resm, lensp, lensm;

  double * sqrt1;
  double * sqrt2;
  double * sqrt3;
  double * sqrt4;
  double * sqrt5;

  double ** cl_md_ic; /* array with argument
                         cl_md_ic[index_md][index_ic1_ic2*spectra_module_.ct_size_+index_ct] */

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

  class_call(lensing_indices(),
             error_message_,
             error_message_);

  /** - put all precision variables hare; will be stored later in precision structure */
  /** - Last element in \f$ \mu \f$ will be for \f$ \mu=1 \f$, needed for sigma2.
      The rest will be chosen as roots of a Gauss-Legendre quadrature **/

  if (ppr->accurate_lensing == _TRUE_) {
    num_mu = (l_unlensed_max_ + ppr->num_mu_minus_lmax); /* Must be even ?? CHECK */
    num_mu += num_mu%2; /* Force it to be even */
  } else {
    /* Integrate correlation function difference on [0,pi/16] */
    num_mu = (l_unlensed_max_ * 2 )/16;
  }
  /** - allocate array of \f$ \mu \f$ values, as well as quadrature weights */

  class_alloc(mu,
              num_mu*sizeof(double),
              error_message_);
  /* Reserve last element of mu for mu=1, needed for sigma2 */
  mu[num_mu-1] = 1.0;

  class_alloc(w8,
              (num_mu-1)*sizeof(double),
              error_message_);

  if (ppr->accurate_lensing == _TRUE_) {

    //debut = omp_get_wtime();
    class_call(quadrature_gauss_legendre(mu,
                                         w8,
                                         num_mu-1,
                                         ppr->tol_gauss_legendre,
                                         error_message_),
               error_message_,
               error_message_);
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
              error_message_);

  class_alloc(d11,
              num_mu*sizeof(double*),
              error_message_);

  class_alloc(d1m1,
              num_mu*sizeof(double*),
              error_message_);

  class_alloc(d2m2,
              num_mu*sizeof(double*),
              error_message_);
  icount += 4*num_mu*(l_unlensed_max_ + 1);

  if (has_te_ == _TRUE_) {

    class_alloc(d20,
                num_mu*sizeof(double*),
                error_message_);

    class_alloc(d3m1,
                num_mu*sizeof(double*),
                error_message_);

    class_alloc(d4m2,
                num_mu*sizeof(double*),
                error_message_);
    icount += 3*num_mu*(l_unlensed_max_ + 1);
  }

  if (has_ee_ == _TRUE_ || has_bb_ == _TRUE_) {

    class_alloc(d22,
                num_mu*sizeof(double*),
                error_message_);

    class_alloc(d31,
                num_mu*sizeof(double*),
                error_message_);

    class_alloc(d3m3,
                num_mu*sizeof(double*),
                error_message_);

    class_alloc(d40,
                num_mu*sizeof(double*),
                error_message_);

    class_alloc(d4m4,
                num_mu*sizeof(double*),
                error_message_);
    icount += 5*num_mu*(l_unlensed_max_ + 1);
  }

  icount += 5*(l_unlensed_max_ + 1); /* for arrays sqrt1[l] to sqrt5[l] */

  /** - Allocate main contiguous buffer **/
  class_alloc(buf_dxx,
              icount * sizeof(double),
              error_message_);

  icount = 0;
  for (index_mu=0; index_mu<num_mu; index_mu++) {

    d00 [index_mu] = &(buf_dxx[icount + (index_mu + 0*num_mu)*(l_unlensed_max_ + 1)]);
    d11 [index_mu] = &(buf_dxx[icount + (index_mu + 1*num_mu)*(l_unlensed_max_ + 1)]);
    d1m1[index_mu] = &(buf_dxx[icount + (index_mu + 2*num_mu)*(l_unlensed_max_ + 1)]);
    d2m2[index_mu] = &(buf_dxx[icount + (index_mu + 3*num_mu)*(l_unlensed_max_ + 1)]);
  }
  icount += 4*num_mu*(l_unlensed_max_ + 1);

  if (has_te_ == _TRUE_) {
    for (index_mu=0; index_mu<num_mu; index_mu++) {
      d20 [index_mu] = &(buf_dxx[icount + (index_mu + 0*num_mu)*(l_unlensed_max_ + 1)]);
      d3m1[index_mu] = &(buf_dxx[icount + (index_mu + 1*num_mu)*(l_unlensed_max_ + 1)]);
      d4m2[index_mu] = &(buf_dxx[icount + (index_mu + 2*num_mu)*(l_unlensed_max_ + 1)]);
    }
    icount += 3*num_mu*(l_unlensed_max_ + 1);
  }

  if (has_ee_ == _TRUE_ || has_bb_ == _TRUE_) {

    for (index_mu=0; index_mu<num_mu; index_mu++) {
      d22 [index_mu] = &(buf_dxx[icount + (index_mu + 0*num_mu)*(l_unlensed_max_ + 1)]);
      d31 [index_mu] = &(buf_dxx[icount + (index_mu + 1*num_mu)*(l_unlensed_max_ + 1)]);
      d3m3[index_mu] = &(buf_dxx[icount + (index_mu + 2*num_mu)*(l_unlensed_max_ + 1)]);
      d40 [index_mu] = &(buf_dxx[icount + (index_mu + 3*num_mu)*(l_unlensed_max_ + 1)]);
      d4m4[index_mu] = &(buf_dxx[icount + (index_mu + 4*num_mu)*(l_unlensed_max_ + 1)]);
    }
    icount += 5*num_mu*(l_unlensed_max_ + 1);
  }

  sqrt1 = &(buf_dxx[icount]);
  icount += l_unlensed_max_ + 1;
  sqrt2 = &(buf_dxx[icount]);
  icount += l_unlensed_max_ + 1;
  sqrt3 = &(buf_dxx[icount]);
  icount += l_unlensed_max_ + 1;
  sqrt4 = &(buf_dxx[icount]);
  icount += l_unlensed_max_ + 1;
  sqrt5 = &(buf_dxx[icount]);
  icount += l_unlensed_max_ + 1;

  //debut = omp_get_wtime();
  class_call(lensing_d00(mu, num_mu, l_unlensed_max_, d00),
             error_message_,
             error_message_);

  class_call(lensing_d11(mu, num_mu, l_unlensed_max_, d11),
             error_message_,
             error_message_);

  class_call(lensing_d1m1(mu, num_mu, l_unlensed_max_, d1m1),
             error_message_,
             error_message_);

  class_call(lensing_d2m2(mu, num_mu, l_unlensed_max_, d2m2),
             error_message_,
             error_message_);
  //fin = omp_get_wtime();
  //cpu_time = (fin-debut);
  //printf("time in lensing_dxx=%4.3f s\n",cpu_time);


  if (has_te_ == _TRUE_) {

    class_call(lensing_d20(mu, num_mu, l_unlensed_max_, d20),
               error_message_,
               error_message_);

    class_call(lensing_d3m1(mu, num_mu, l_unlensed_max_, d3m1),
               error_message_,
               error_message_);

    class_call(lensing_d4m2(mu, num_mu, l_unlensed_max_, d4m2),
               error_message_,
               error_message_);

  }

  if (has_ee_ == _TRUE_ || has_bb_ == _TRUE_) {

    class_call(lensing_d22(mu, num_mu, l_unlensed_max_, d22),
               error_message_,
               error_message_);

    class_call(lensing_d31(mu, num_mu, l_unlensed_max_, d31),
               error_message_,
               error_message_);

    class_call(lensing_d3m3(mu, num_mu, l_unlensed_max_, d3m3),
               error_message_,
               error_message_);

    class_call(lensing_d40(mu, num_mu, l_unlensed_max_, d40),
               error_message_,
               error_message_);

    class_call(lensing_d4m4(mu, num_mu, l_unlensed_max_, d4m4),
               error_message_,
               error_message_);
  }

  /** - compute \f$ Cgl(\mu)\f$, \f$ Cgl2(\mu) \f$ and sigma2(\f$\mu\f$) */

  class_alloc(Cgl,
              num_mu*sizeof(double),
              error_message_);

  class_alloc(Cgl2,
              num_mu*sizeof(double),
              error_message_);

  class_alloc(sigma2,
              (num_mu-1)*sizeof(double), /* Zero separation is omitted */
              error_message_);

  class_alloc(cl_unlensed,
              spectra_module_.ct_size_*sizeof(double),
              error_message_);


  /** - Locally store unlensed temperature \f$ cl_{tt}\f$ and potential \f$ cl_{pp}\f$ spectra **/
  class_alloc(cl_tt,
              (l_unlensed_max_ + 1)*sizeof(double),
              error_message_);
  if (has_te_ == _TRUE_) {
    class_alloc(cl_te,
                (l_unlensed_max_ + 1)*sizeof(double),
                error_message_);
  }
  if (has_ee_ == _TRUE_ || has_bb_ == _TRUE_) {
    class_alloc(cl_ee,
                (l_unlensed_max_ + 1)*sizeof(double),
                error_message_);

    class_alloc(cl_bb,
                (l_unlensed_max_ + 1)*sizeof(double),
                error_message_);
  }
  class_alloc(cl_pp,
              (l_unlensed_max_ + 1)*sizeof(double),
              error_message_);

  class_alloc(cl_md_ic,
              spectra_module_.md_size_*sizeof(double *),
              error_message_);

  class_alloc(cl_md,
              spectra_module_.md_size_*sizeof(double *),
              error_message_);

  for (index_md = 0; index_md < spectra_module_.md_size_; index_md++) {

    if (spectra_module_.md_size_ > 1)

      class_alloc(cl_md[index_md],
                  spectra_module_.ct_size_*sizeof(double),
                  error_message_);

    if (spectra_module_.ic_size_[index_md] > 1)

      class_alloc(cl_md_ic[index_md],
                  spectra_module_.ic_ic_size_[index_md]*spectra_module_.ct_size_*sizeof(double),
                  error_message_);
  }

  for (l=2; l<=l_unlensed_max_; l++) {
    class_call(spectra_module_.spectra_cl_at_l(l, cl_unlensed, cl_md, cl_md_ic),
               psp->error_message,
               error_message_);
    cl_tt[l] = cl_unlensed[index_lt_tt_];
    cl_pp[l] = cl_unlensed[index_lt_pp_];
    if (has_te_ == _TRUE_) {
      cl_te[l] = cl_unlensed[index_lt_te_];
    }
    if (has_ee_ == _TRUE_ || has_bb_ == _TRUE_) {
      cl_ee[l] = cl_unlensed[index_lt_ee_];
      cl_bb[l] = cl_unlensed[index_lt_bb_];
    }
  }

  for (index_md = 0; index_md < spectra_module_.md_size_; index_md++) {

    if (spectra_module_.md_size_ > 1)
      free(cl_md[index_md]);

    if (spectra_module_.ic_size_[index_md] > 1)
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

    for (l=2; l<=l_unlensed_max_; l++) {

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
  if (has_tt_ == _TRUE_) {

    class_calloc(ksi,
                 (num_mu-1),
                 sizeof(double),
                 error_message_);
  }

  /** - --> ksiX is for TE **/
  if (has_te_ == _TRUE_) {

    class_calloc(ksiX,
                 (num_mu-1),
                 sizeof(double),
                 error_message_);
  }

  /** - --> ksip, ksim for EE, BB **/
  if (has_ee_ == _TRUE_ || has_bb_ == _TRUE_) {

    class_calloc(ksip,
                 (num_mu-1),
                 sizeof(double),
                 error_message_);

    class_calloc(ksim,
                 (num_mu-1),
                 sizeof(double),
                 error_message_);
  }

  for (l=2; l<=l_unlensed_max_; l++) {

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

    for (l=2; l<=l_unlensed_max_; l++) {

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

      if (has_te_ == _TRUE_ || has_ee_ == _TRUE_ || has_bb_ == _TRUE_) {
        /* X_022 = exp(-(fac-1.)*sigma2[index_mu]); */
        X_022 = X_000 * (1+sigma2[index_mu]*(1+0.5*sigma2[index_mu])); /* Order 2 */
        X_p022 = -(fac-1.)*X_022; /* Old versions were missing the
        minus sign in this line, which introduced a very small error
        on the high-l C_l^TE lensed spectrum [credits for bug fix:
        Selim Hotinli] */

        /* X_242 = 0.25*sqrt4[l] * exp(-(fac-5./2.)*sigma2[index_mu]); */
        X_242 = 0.25*sqrt4[l] * X_000; /* Order 0 */
        if (has_ee_ == _TRUE_ || has_bb_ == _TRUE_) {

          /* X_121 = - 0.5*sqrt2[l] * exp(-(fac-2./3.)*sigma2[index_mu]);
             X_132 = - 0.5*sqrt3[l] * exp(-(fac-5./3.)*sigma2[index_mu]); */
          X_121 = -0.5*sqrt2[l] * X_000 * (1+2./3.*sigma2[index_mu]); /* Order 1 */
          X_132 = -0.5*sqrt3[l] * X_000 * (1+5./3.*sigma2[index_mu]); /* Order 1 */
        }
      }


      if (has_tt_ == _TRUE_) {

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

      if (has_te_ == _TRUE_) {

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

      if (has_ee_ == _TRUE_ || has_bb_ == _TRUE_) {

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
  if (has_tt_ == _TRUE_) {
    class_call(lensing_lensed_cl_tt(ksi, d00, w8, num_mu - 1),
               error_message_,
               error_message_);
    if (ppr->accurate_lensing == _FALSE_) {
      class_call(lensing_addback_cl_tt(cl_tt),
                 error_message_,
                 error_message_);
    }
  }

  if (has_te_ == _TRUE_) {
    class_call(lensing_lensed_cl_te(ksiX, d20, w8, num_mu - 1),
               error_message_,
               error_message_);
    if (ppr->accurate_lensing == _FALSE_) {
      class_call(lensing_addback_cl_te(cl_te),
                 error_message_,
                 error_message_);
    }
  }

  if (has_ee_ == _TRUE_ || has_bb_ == _TRUE_) {

    class_call(lensing_lensed_cl_ee_bb(ksip, ksim, d22, d2m2, w8, num_mu - 1),
               error_message_,
               error_message_);
    if (ppr->accurate_lensing == _FALSE_) {
      class_call(lensing_addback_cl_ee_bb(cl_ee,cl_bb),
                 error_message_,
                 error_message_);
    }
  }
  //fin=omp_get_wtime();
  //cpu_time = (fin-debut);
  //printf("time in final lensing computation=%4.3f s\n",cpu_time);

  /** - spline computed \f$ C_l\f$'s in view of interpolation */

  class_call(array_spline_table_lines(l_,
                                      l_size_,
                                      cl_lens_,
                                      lt_size_,
                                      ddcl_lens_,
                                      _SPLINE_EST_DERIV_,
                                      error_message_),
             error_message_,
             error_message_);

  /** - Free lots of stuff **/
  free(buf_dxx);

  free(d00);
  free(d11);
  free(d1m1);
  free(d2m2);
  if (has_te_ == _TRUE_) {
    free(d20);
    free(d3m1);
    free(d4m2);
  }
  if (has_ee_ == _TRUE_ || has_bb_ == _TRUE_) {
    free(d22);
    free(d31);
    free(d3m3);
    free(d40);
    free(d4m4);
  }

  if (has_tt_ == _TRUE_)
    free(ksi);
  if (has_te_ == _TRUE_)
    free(ksiX);
  if (has_ee_ == _TRUE_ || has_bb_ == _TRUE_) {
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
  if (has_te_ == _TRUE_)
    free(cl_te);
  if (has_ee_ == _TRUE_ || has_bb_ == _TRUE_) {
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

int LensingModule::lensing_free() {

  if (ple->has_lensed_cls == _TRUE_) {

    free(l_);
    free(cl_lens_);
    free(ddcl_lens_);
    free(l_max_lt_);

  }

  return _SUCCESS_;

}

/**
 * This routine defines indices and allocates tables in the lensing structure
 *
 * @param ppr  Input: pointer to precision structure
 * @param psp  Input: pointer to spectra structure
 * @param ple  Input/output: pointer to lensing structure
 * @return the error status
 */

int LensingModule::lensing_indices(){

  int index_l;

  double ** cl_md_ic; /* array with argument
                         cl_md_ic[index_md][index_ic1_ic2*spectra_module_.ct_size_+index_ct] */

  double ** cl_md;    /* array with argument
                         cl_md[index_md][index_ct] */

  int index_md;
  int index_lt;

  /* indices of all Cl types (lensed and unlensed) */

  if (spectra_module_.has_tt_ == _TRUE_) {
    has_tt_ = _TRUE_;
    index_lt_tt_ = spectra_module_.index_ct_tt_;
  }
  else {
    has_tt_ = _FALSE_;
  }

  if (spectra_module_.has_ee_ == _TRUE_) {
    has_ee_ = _TRUE_;
    index_lt_ee_ = spectra_module_.index_ct_ee_;
  }
  else {
    has_ee_ = _FALSE_;
  }

  if (spectra_module_.has_te_ == _TRUE_) {
    has_te_ = _TRUE_;
    index_lt_te_ = spectra_module_.index_ct_te_;
  }
  else {
    has_te_ = _FALSE_;
  }

  if (spectra_module_.has_bb_ == _TRUE_) {
    has_bb_ = _TRUE_;
    index_lt_bb_ = spectra_module_.index_ct_bb_;
  }
  else {
    has_bb_ = _FALSE_;
  }

  if (spectra_module_.has_pp_ == _TRUE_) {
    has_pp_ = _TRUE_;
    index_lt_pp_ = spectra_module_.index_ct_pp_;
  }
  else {
    has_pp_ = _FALSE_;
  }

  if (spectra_module_.has_tp_ == _TRUE_) {
    has_tp_ = _TRUE_;
    index_lt_tp_ = spectra_module_.index_ct_tp_;
  }
  else {
    has_tp_ = _FALSE_;
  }

  if (spectra_module_.has_dd_ == _TRUE_) {
    has_dd_ = _TRUE_;
    index_lt_dd_ = spectra_module_.index_ct_dd_;
  }
  else {
    has_dd_ = _FALSE_;
  }

  if (spectra_module_.has_td_ == _TRUE_) {
    has_td_ = _TRUE_;
    index_lt_td_ = spectra_module_.index_ct_td_;
  }
  else {
    has_td_ = _FALSE_;
  }

  if (spectra_module_.has_ll_ == _TRUE_) {
    has_ll_ = _TRUE_;
    index_lt_ll_ = spectra_module_.index_ct_ll_;
  }
  else {
    has_ll_ = _FALSE_;
  }

  if (spectra_module_.has_tl_ == _TRUE_) {
    has_tl_ = _TRUE_;
    index_lt_tl_ = spectra_module_.index_ct_tl_;
  }
  else {
    has_tl_ = _FALSE_;
  }

  lt_size_ = spectra_module_.ct_size_;

  /* number of multipoles */

  l_unlensed_max_ = spectra_module_.l_max_tot_;

  l_lensed_max_ = l_unlensed_max_ - ppr->delta_l_max;

  for (index_l = 0; (index_l < spectra_module_.l_size_max_) && (spectra_module_.l_[index_l] <= l_lensed_max_); index_l++);

  if (index_l < spectra_module_.l_size_max_) index_l++; /* one more point in order to be able to interpolate till l_lensed_max_ */

  l_size_ = index_l+1;

  class_alloc(l_, l_size_*sizeof(double), error_message_);

  for (index_l=0; index_l < l_size_; index_l++) {

    l_[index_l] = spectra_module_.l_[index_l];

  }

  /* allocate table where results will be stored */

  class_alloc(cl_lens_,
              l_size_*lt_size_*sizeof(double),
              error_message_);

  class_alloc(ddcl_lens_,
              l_size_*lt_size_*sizeof(double),
              error_message_);

  /* fill with unlensed cls */

  class_alloc(cl_md_ic,
              spectra_module_.md_size_*sizeof(double *),
              error_message_);

  class_alloc(cl_md,
              spectra_module_.md_size_*sizeof(double *),
              error_message_);

  for (index_md = 0; index_md < spectra_module_.md_size_; index_md++) {

    if (spectra_module_.md_size_ > 1)

      class_alloc(cl_md[index_md],
                  spectra_module_.ct_size_*sizeof(double),
                  error_message_);

    if (spectra_module_.ic_size_[index_md] > 1)

      class_alloc(cl_md_ic[index_md],
                  spectra_module_.ic_ic_size_[index_md]*spectra_module_.ct_size_*sizeof(double),
                  error_message_);
  }

  for (index_l=0; index_l<l_size_; index_l++) {

    class_call(spectra_module_.spectra_cl_at_l(l_[index_l], &(cl_lens_[index_l*lt_size_]), cl_md, cl_md_ic),
               psp->error_message,
               error_message_);

  }

  for (index_md = 0; index_md < spectra_module_.md_size_; index_md++) {

    if (spectra_module_.md_size_ > 1)
      free(cl_md[index_md]);

    if (spectra_module_.ic_size_[index_md] > 1)
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

  class_alloc(l_max_lt_, lt_size_*sizeof(double), error_message_);
  for (index_lt = 0; index_lt < lt_size_; index_lt++) {
    l_max_lt_[index_lt] = 0.;
    for (index_md = 0; index_md < spectra_module_.md_size_; index_md++) {
      l_max_lt_[index_lt] = MAX(l_max_lt_[index_lt], spectra_module_.l_max_ct_[index_md][index_lt]);

      if ((has_bb_ == _TRUE_) && (has_ee_ == _TRUE_) && (index_lt == index_lt_bb_)) {
        l_max_lt_[index_lt] = MAX(l_max_lt_[index_lt], spectra_module_.l_max_ct_[index_md][index_lt_ee_]);
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


int LensingModule::lensing_lensed_cl_tt(double *ksi, double **d00, double *w8, int nmu) {

  double cle;
  int imu;
  int index_l;

  /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
  private (imu,index_l,cle)                     \
  schedule (static)

  for(index_l=0; index_l<l_size_; index_l++){
    cle=0;
    for (imu=0;imu<nmu;imu++) {
      cle += ksi[imu]*d00[imu][(int)l_[index_l]]*w8[imu]; /* loop could be optimized */
    }
    cl_lens_[index_l*lt_size_ + index_lt_tt_] = cle*2.0*_PI_;
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

int LensingModule::lensing_addback_cl_tt(double *cl_tt) {
  int index_l, l;

  for (index_l=0; index_l<l_size_; index_l++) {
    l = (int)l_[index_l];
    cl_lens_[index_l*lt_size_ + index_lt_tt_] += cl_tt[l];
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


int LensingModule::lensing_lensed_cl_te(double *ksiX, double **d20, double *w8, int nmu) {

  double clte;
  int imu;
  int index_l;

  /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
  private (imu,index_l,clte)                    \
  schedule (static)

  for(index_l=0; index_l < l_size_; index_l++){
    clte=0;
    for (imu=0;imu<nmu;imu++) {
      clte += ksiX[imu]*d20[imu][(int)l_[index_l]]*w8[imu]; /* loop could be optimized */
    }
    cl_lens_[index_l*lt_size_ + index_lt_te_] = clte*2.0*_PI_;
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

int LensingModule::lensing_addback_cl_te(double *cl_te) {
  int index_l, l;

  for (index_l=0; index_l<l_size_; index_l++) {
    l = (int)l_[index_l];
    cl_lens_[index_l*lt_size_ + index_lt_te_] += cl_te[l];
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


int LensingModule::lensing_lensed_cl_ee_bb(double *ksip, double *ksim, double **d22, double **d2m2, double *w8, int nmu) {

  double clp, clm;
  int imu;
  int index_l;

  /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
  private (imu,index_l,clp,clm)                 \
  schedule (static)

  for(index_l=0; index_l < l_size_; index_l++){
    clp=0; clm=0;
    for (imu=0;imu<nmu;imu++) {
      clp += ksip[imu]*d22[imu][(int)l_[index_l]]*w8[imu]; /* loop could be optimized */
      clm += ksim[imu]*d2m2[imu][(int)l_[index_l]]*w8[imu]; /* loop could be optimized */
    }
    cl_lens_[index_l*lt_size_ + index_lt_ee_] = (clp + clm)*_PI_;
    cl_lens_[index_l*lt_size_ + index_lt_bb_] = (clp - clm)*_PI_;
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

int LensingModule::lensing_addback_cl_ee_bb(double * cl_ee, double * cl_bb) {

  int index_l, l;

  for (index_l=0; index_l<l_size_; index_l++) {
    l = (int)l_[index_l];
    cl_lens_[index_l*lt_size_ + index_lt_ee_] += cl_ee[l];
    cl_lens_[index_l*lt_size_ + index_lt_bb_] += cl_bb[l];
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

int LensingModule::lensing_d00(
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

int LensingModule::lensing_d11(
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

int LensingModule::lensing_d1m1(
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

int LensingModule::lensing_d2m2(
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

int LensingModule::lensing_d22(
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

int LensingModule::lensing_d20(
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

int LensingModule::lensing_d31(
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

int LensingModule::lensing_d3m1(
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

int LensingModule::lensing_d3m3(
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

int LensingModule::lensing_d40(
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

int LensingModule::lensing_d4m2(
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

int LensingModule::lensing_d4m4(
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
