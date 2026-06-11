/** @file lensing.c Documented lensing module
 *
 * Simon Prunet and Julien Lesgourgues, 6.12.2010
 * Improved by Cyril Pitrou, 6.2026
 *
 * This module computes the lensed temperature and polarization
 * anisotropy power spectra \f$ C_l^{X}, P(k), ... \f$'s given the
 * unlensed temperature, polarization and lensing potential spectra.
 *
 * Follows Challinor & Lewis full-sky method, astro-ph/0502425.
 * Follows additionally Lewis & Challinor astro-ph/0601594 for
 * higher order corrections.
 *
 * The following functions can be called from other modules:
 *
 * -# lensing_init() at the beginning (but after harmonic_init())
 * -# lensing_cl_at_l() at any time for computing Cl_lensed at any l
 * -# lensing_free() at the end
 */

#include "lensing.h"
#include "parallel.h"

/**
 * Anisotropy power spectra \f$ C_l\f$'s for all types, modes and initial conditions.
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
  /* The following Wigner dm1m2 functions are needed when using higher order lensing according to 9.12 and 9.16-9.18 of astro-ph/0601594 */
  double ** d5m1 = NULL;
  double ** d5m3 = NULL;
  double ** d6m2 = NULL;
  double ** d33 = NULL;
  double * buf_dxx; /* buffer */

  double * Cgl;   /* Cgl[index_mu] */
  double * Cgl2;  /* Cgl2[index_mu] */
  double * sigma2; /* sigma[index_mu] */

  double * ksi = NULL;  /* ksi[index_mu] */
  double * ksiX = NULL;  /* ksiX[index_mu] */
  double * ksip = NULL;  /* ksip[index_mu] */
  double * ksim = NULL;  /* ksim[index_mu] */

  int num_mu,index_mu;
  long long icount;
  int l;
  double ll;
  double * cl_unlensed;  /* cl_unlensed[index_ct] */
  double * cl_tt; /* unlensed  cl, to be filled to avoid repeated calls to harmonic_cl_at_l */
  double * cl_te = NULL; /* unlensed  cl, to be filled to avoid repeated calls to harmonic_cl_at_l */
  double * cl_ee = NULL; /* unlensed  cl, to be filled to avoid repeated calls to harmonic_cl_at_l */
  double * cl_bb = NULL; /* unlensed  cl, to be filled to avoid repeated calls to harmonic_cl_at_l */
  double * cl_pp; /* potential cl, to be filled to avoid repeated calls to harmonic_cl_at_l */

  double * sqrt1;
  double * sqrt2;
  double * sqrt3;
  double * sqrt4;
  double * sqrt5; /* We now invert the definition */

  double ** cl_md_ic; /* array with argument
                         cl_md_ic[index_md][index_ic1_ic2*phr->ct_size+index_ct] */

  double ** cl_md;    /* array with argument
                         cl_md[index_md][index_ct] */

  int index_md;

  /** - check that we really want to compute at least one spectrum */

  if (ple->has_lensed_cls == _FALSE_) {
    if (ple->lensing_verbose > 0)
      printf("No lensing requested. Lensing module skipped.\n");
    return _SUCCESS_;
  }
  else {
    if (ple->lensing_verbose > 0) {
      printf("Computing lensed spectra at order %d in Cgl and %d in Cgl2 ",ppr->lensing_C0_order,ppr->lensing_C2_order);
      if (ppr->accurate_lensing==_TRUE_)
        printf("(accurate mode)\n");
      else
        printf("(fast mode)\n");
    }
  }

  /** - Consistency checks */

  class_test(((ppr->lensing_C0_order > 1)||(ppr->lensing_C0_order < 0)),
	     ple->error_message,
	     "The order of the monopole of lensing correlation C0 should be either 0 or 1 but you passed  C0 = %d.", ppr->lensing_C0_order);

  class_test(((ppr->lensing_C2_order > 4)||(ppr->lensing_C2_order < 0)),
	     ple->error_message,
	     "The order of the quadrupole of lensing correlation C2 should be either 0, 1, 2, 3 or 4 but you passed  C2 = %d.", ppr->lensing_C2_order);

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
    /* Integrate correlation function difference on [0,pi/ppr->non_accurate_lensing_boundary]  (with default setting ppr->non_accurate_lensing_boundary = 16) */
    num_mu = (ple->l_unlensed_max * 2 )/ppr->non_accurate_lensing_boundary;
    /* When l_unlensed_max is large (typically above 2000), this
	method can be very inaccurate. As explained in II.B.3 of
	astro-ph/0502425, this reduction of the upper boundary introduces
	ringing on small scales. It is advised to use the accurate method
	(Gauss-Legendre quadrature) when l_max > 2000. Also note that this
	fast method might not be very accurate around l_unlensed_max since
	the density of points in the Riemann integration is close to the
	Shannon criterium. One could increase the boundary and the number
	of points with ppr->non_accurate_lensing_boundary=8, 4, 2, 1 to
	avoid this issue, but again at the price of slowing the numerical
	evaluation of the Riemann integral. If one wants a precise result
	one should switch to the accurate lensing method (the
	Gauss-Legendre quadrature)
     */
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

    class_call(quadrature_gauss_legendre(mu,
                                         w8,
                                         num_mu-1,
                                         ppr->tol_gauss_legendre,
                                         ple->error_message),
               ple->error_message,
               ple->error_message);

  } else { /* Crude integration on [0,pi/16]: Riemann sum on theta */

    delta_theta = _PI_/ppr->non_accurate_lensing_boundary / (double)(num_mu-1);

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

  if (ple->has_te==_TRUE_) {

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

  /* The following Wigner dlm1m2 functions are only needed when using
     higher order lensing (order 3 or 4 in the expansion in
     Cgl_2*l(l+1) see 9.12 - 9.16 of astro-ph/0601594) We do not
     distinguish the case where the user only wants TT correlations
     and no polarization and we compute all d-Wigner needed for TT TE
     EE and BB for these higher orders.
   */
  if ((ppr->lensing_C2_order > 2) || (ppr->lensing_C0_order > 0)) {
    class_alloc(d5m1,
                num_mu*sizeof(double*),
                ple->error_message);

    class_alloc(d5m3,
                num_mu*sizeof(double*),
                ple->error_message);

    class_alloc(d6m2,
                num_mu*sizeof(double*),
                ple->error_message);

    class_alloc(d33,
                num_mu*sizeof(double*),
                ple->error_message);

    icount += 4*num_mu*(ple->l_unlensed_max+1);
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

  /* These Wigner d functions are only needed when considering higher order terms in the expansion in Cgl_2*l*(l+1) (third and fourth order). */
  if ((ppr->lensing_C2_order > 2) || (ppr->lensing_C0_order > 0)) {

    for (index_mu=0; index_mu<num_mu; index_mu++) {
      d5m1[index_mu] = &(buf_dxx[icount+index_mu           * (ple->l_unlensed_max+1)]);
      d5m3[index_mu] = &(buf_dxx[icount+(index_mu+num_mu)  * (ple->l_unlensed_max+1)]);
      d6m2[index_mu]= &(buf_dxx[icount+(index_mu+2*num_mu) * (ple->l_unlensed_max+1)]);
      d33[index_mu]= &(buf_dxx[icount+(index_mu+3*num_mu)  * (ple->l_unlensed_max+1)]);
    }
    icount += 4*num_mu*(ple->l_unlensed_max+1);
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

  class_call(lensing_dm1m2(ple,mu,num_mu,ple->l_unlensed_max,0,0,d00),
             ple->error_message,
             ple->error_message);

  class_call(lensing_dm1m2(ple,mu,num_mu,ple->l_unlensed_max,1,1,d11),
             ple->error_message,
             ple->error_message);

  class_call(lensing_dm1m2(ple,mu,num_mu,ple->l_unlensed_max,1,-1,d1m1),
             ple->error_message,
             ple->error_message);

  class_call(lensing_dm1m2(ple,mu,num_mu,ple->l_unlensed_max,2,-2,d2m2),
             ple->error_message,
             ple->error_message);

  if (ple->has_te==_TRUE_) {

    class_call(lensing_dm1m2(ple,mu,num_mu,ple->l_unlensed_max,2,0,d20),
               ple->error_message,
               ple->error_message);

    class_call(lensing_dm1m2(ple,mu,num_mu,ple->l_unlensed_max,3,-1,d3m1),
               ple->error_message,
               ple->error_message);

    class_call(lensing_dm1m2(ple,mu,num_mu,ple->l_unlensed_max,4,-2,d4m2),
               ple->error_message,
               ple->error_message);

  }

  if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {

    class_call(lensing_dm1m2(ple,mu,num_mu,ple->l_unlensed_max,2,2,d22),
               ple->error_message,
               ple->error_message);

    class_call(lensing_dm1m2(ple,mu,num_mu,ple->l_unlensed_max,3,1,d31),
               ple->error_message,
               ple->error_message);

    class_call(lensing_dm1m2(ple,mu,num_mu,ple->l_unlensed_max,3,-3,d3m3),
               ple->error_message,
               ple->error_message);

    class_call(lensing_dm1m2(ple,mu,num_mu,ple->l_unlensed_max,4,0,d40),
               ple->error_message,
               ple->error_message);

    class_call(lensing_dm1m2(ple,mu,num_mu,ple->l_unlensed_max,4,-4,d4m4),
               ple->error_message,
               ple->error_message);

  }

  /* We also add the Wigner d functions needed for higher order lensing */
  if ((ppr->lensing_C2_order > 2) || (ppr->lensing_C0_order > 0)) {

    class_call(lensing_dm1m2(ple,mu,num_mu,ple->l_unlensed_max,5,-1,d5m1),
               ple->error_message,
               ple->error_message);

    class_call(lensing_dm1m2(ple,mu,num_mu,ple->l_unlensed_max,5,-3,d5m3),
               ple->error_message,
               ple->error_message);

    class_call(lensing_dm1m2(ple,mu,num_mu,ple->l_unlensed_max,6,-2,d6m2),
               ple->error_message,
               ple->error_message);

    class_call(lensing_dm1m2(ple,mu,num_mu,ple->l_unlensed_max,3,3,d33),
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

  class_setup_parallel();

  for (index_mu=0; index_mu<num_mu; index_mu++) {

    int l_unlensed_max;
    l_unlensed_max = ple->l_unlensed_max;
    class_run_parallel(with_arguments(index_mu,l_unlensed_max,Cgl,Cgl2,cl_pp,d11,d1m1),
    int l;

    Cgl[index_mu]=0;
    Cgl2[index_mu]=0;

    for (l=2; l<=l_unlensed_max; l++) {

      /* Eqs. 35 of astro-ph/0502425. These are the monopole and
         quadrupole part of the correlation of lensing
         displacement. */
      Cgl[index_mu] += (2.*l+1.)*l*(l+1.)*
        cl_pp[l]*d11[index_mu][l];

      Cgl2[index_mu] += (2.*l+1.)*l*(l+1.)*
        cl_pp[l]*d1m1[index_mu][l];

    }

    Cgl[index_mu] /= 4.*_PI_;
    Cgl2[index_mu] /= 4.*_PI_;
    return _SUCCESS_;
    );

  }

  class_finish_parallel();

  for (index_mu=0; index_mu<num_mu-1; index_mu++) {
    /* Cgl(1.0) - Cgl(mu) */
    sigma2[index_mu] = Cgl[num_mu-1] - Cgl[index_mu];
  }

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
    sqrt5[l]=1./sqrt(ll*(ll+1)); /* The definition is now inverted. */
  }

  /* We now compute the lensed correlation functions. It is based on
      an expansion in Clg (monopole of correlation of lensing) and
      Clg2 (quadrupole of correlation lensing). The strategy we adopt
      is to use astro-ph/0502425 for the order Cgl2^0, Cgl2^1 and
      Cgl2^2 corrections, but to use 9.12 and 9.16-018 of
      astro-ph/0601594 for order Cgl2^3, Cgl2^4 and Cgl^1. The method
      of astro-ph/0502425 includes corrections which are of order 1/l
      or 1/l^2, and which are not in astro-ph/0601594 so we must make
      sure not to ignore these 1/l or 1/l^2 corrections.
   */

  for (index_mu=0;index_mu<num_mu-1;index_mu++) {

    // = means that all dependencies are captured.
    class_run_parallel(=,

    int l;
    double declare_list_of_variables_inside_parallel_region(ll,fac, fac1, fac2, X_000, X_p000, X_220,X_022,X_p022,X_121,X_132,X_242,x,X_000_square);
    double declare_list_of_variables_inside_parallel_region(res,resX,resp,resm,lens,lensp,lensm);
    for (l=2;l<=ple->l_unlensed_max;l++) {

      ll = (double)l;

      fac = ll*(ll+1)/4.;
      fac1 = (2*ll+1)/(4.*_PI_);
      fac2 = 1/(ll*(ll+1));

      /* Up to version 3.3.4 the strategy was to keep terms of the
	  form (sigma2)^k*(Cgl2)^m with k+m <= 2. Now that we allow the
	  possibility to consider third and fourth order in Cgl2*l*(l+1),
	  we do not perform such crude truncation. However we also do not
	  want to compute too many terms which involve an exponential
	  since this is slow. The dominant contribution from sigma2 is the
	  global prefactor exp(-l(l+1)sigma^2/2) which is the square of
	  X_000. But we also encounter other factors of the type
	  exp(-number * sigma2), where number does not grow with l, and we
	  replace these by (1- number * sigma2).  This should always be
	  sufficient since sigma is never larger than 10^-7 for standard
	  cosmology (see e.g. Fig. 2 in astro-ph/0502425)
       */

      /* These X_{ijk} functions are defined in eqs 39, 40, 57-60 of astro-ph/0502425 */
      X_000 = exp(-fac*sigma2[index_mu]); /* Eq. 39 in astro-ph/0502425. The numerical evaluation of this exponential takes quite some time */
      X_p000 = -fac*X_000;
      X_000_square = X_000 * X_000;

      /* Up to version 3.3.4, the expression X_220 = 0.25*sqrt1[l] *
         X_000 was used, but it is an approximation which consists in
         neglecting a 1/l^2 correction, and we do not want that. The
         correct expression would be X_220 = 0.25 * sqrt1[l] * X_000 *
         exp(0.5*sigma2[index_mu]) (Eq. 40 in astro-ph/0502425) and
         according to the method detailed above we use: */
      X_220 = 0.25*sqrt1[l] * X_000 * (1 +  0.5*sigma2[index_mu]);


      /* next 5 lines useless, but avoid compiler warning 'may be used uninitialized' */
      X_242=0.;
      X_132=0.;
      X_121=0.;
      X_p022=0.;
      X_022=0.;

      if ((ppr->lensing_C2_order > 2) || (ppr->lensing_C0_order > 0))
        x = 2*fac*Cgl2[index_mu]; /* l(l+1)/2 C_2 */

      if (ple->has_te==_TRUE_ || ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
        /* X_022 = exp(-(fac-1.)*sigma2[index_mu]); */
        X_022 = X_000 * (1+sigma2[index_mu]*(1+0.5*sigma2[index_mu])); /* Order 2 of (57) in astro-ph/0502425 */

        X_p022 = -(fac-1.)*X_022; /* Old versions were missing the
                                     minus sign in this line, which introduced a very small error
                                     on the high-l C_l^TE lensed spectrum [credits for bug fix:
                                     Selim Hotinli] */

	    /* Up to version 3.3.4 we used X_242 = 0.25*sqrt4[l] * X_000,
	    but the full expression is:
        X_242 = 0.25*sqrt4[l] * exp(-(fac-5./2.)*sigma2[index_mu])
        which we approximate by: */
        X_242 = 0.25*sqrt4[l] * X_000 * (1+ 2.5* sigma2[index_mu]);

        if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {

          /* X_121 = - 0.5*sqrt2[l] * exp(-(fac-2./3.)*sigma2[index_mu]);
             X_132 = - 0.5*sqrt3[l] * exp(-(fac-5./3.)*sigma2[index_mu]); */
          X_121 = -0.5*sqrt2[l] * X_000 * (1+2./3.*sigma2[index_mu]); /* Order 1 */
          X_132 = -0.5*sqrt3[l] * X_000 * (1+5./3.*sigma2[index_mu]); /* Order 1 */
        }
      }

      if (ple->has_tt==_TRUE_) {

        res = fac1*cl_tt[l];

        lens=X_000_square *d00[index_mu][l]; /* Order Cgl2^0, meaning we only have the effect of exp(-l*(l+1)sigma2/2) */

        /* All sky method of astro-ph/0502425 */
        if (ppr->lensing_C2_order >= 1) {
          lens += X_p000*X_p000*d1m1[index_mu][l]*Cgl2[index_mu]*8.*fac2; /* First order in Cgl2 in Eq. 38 of 0502425 */
          if (ppr->lensing_C2_order >= 2) {
            lens += (X_p000*X_p000*d00[index_mu][l] + X_220*X_220*d2m2[index_mu][l])*Cgl2[index_mu]*Cgl2[index_mu]; /* second order in Cgl2 in Eq. 38 of 0502425 */
            if (ppr->lensing_C2_order >= 3) {
              lens += X_000_square * x*x*x * (1/8.*d1m1[index_mu][l] + 1/24.*d3m3[index_mu][l]); /* Third order in Cgl2 deduced from expansion of 9.12 of 0601594 */
              if (ppr->lensing_C2_order >= 4) {
                lens += X_000_square * x*x*x*x * (1/64.*d00[index_mu][l] + 1/48.*d2m2[index_mu][l] + 1/192.*d4m4[index_mu][l]); /* Fourth order in Cgl2 deduced from expansion of 9.12 of 0601594 */
              }
            }
          }
        }
        if (ppr->lensing_C0_order == 1) {
          lens += (2*X_000*X_p000*d00[index_mu][l]  + 8.*fac2*X_p000*X_p000*d11[index_mu][l])*Cgl[index_mu];
          /* first order in Cgl (no coupling to Cgl2 since very small) from C1 of 0502425. */
        }

        /* Old implementation: */
        /*lens = (X_000*X_000*d00[index_mu][l] +
                X_p000*X_p000*d1m1[index_mu][l]
                *Cgl2[index_mu]*8./(ll*(ll+1)) +
                (X_p000*X_p000*d00[index_mu][l] +
                 X_220*X_220*d2m2[index_mu][l])
                 *Cgl2[index_mu]*Cgl2[index_mu]);*/
        if (ppr->accurate_lensing == _FALSE_) {
          /* Remove unlensed correlation function */
          lens -= d00[index_mu][l];
        }
        res *= lens;
        ksi[index_mu] += res;
      }

      if (ple->has_te==_TRUE_) {

        resX = fac1*cl_te[l];

        lens = X_022*X_000*d20[index_mu][l]; /* Order Cgl2^0, meaning we essentially only have the effect of exp(-l*(l+1)sigma2/2) */

        if (ppr->lensing_C2_order >= 1 ) {
          lens += Cgl2[index_mu]*2.*X_p000*sqrt5[l] * (X_121*d11[index_mu][l] + X_132*d3m1[index_mu][l]); /* First order in Cgl2 from 56 of 0502425 but with typo corrected (X_112 replaced by X_121) */
          if (ppr->lensing_C2_order >= 2 ) {
            lens += 0.5 * Cgl2[index_mu] * Cgl2[index_mu] *  ( ( 2.*X_p022*X_p000+X_220*X_220 ) * d20[index_mu][l] + X_220*X_242*d4m2[index_mu][l] ); /* Second order in Cgl2 from 56 of 0502425 */
            if (ppr->lensing_C2_order >= 3 ) {
              lens+= X_000_square * x*x*x * (1/16.*d11[index_mu][l] + 1/12.*d3m1[index_mu][l] + 1/48.*d5m3[index_mu][l]); /* Third order in Cgl2 from 9.18 of 0601594 */
              if (ppr->lensing_C2_order >= 4 ) {
                lens += X_000_square * x*x*x*x * (10/384.*d20[index_mu][l]  +5/384.*d4m2[index_mu][l] ); /* Fourth order in Cgl2 from 9.18 of 0601594 */
              }
            }
          }
        }
        if (ppr->lensing_C0_order == 1) {
          lens += ( (X_022*X_p000 + X_p022*X_000)*d20[index_mu][l]  + 2.*sqrt5[l]*(X_p000*X_132*d31[index_mu][l] + X_121*X_p000*d1m1[index_mu][l]) )*Cgl[index_mu];
          /* first order in Cgl (no coupling to Cgl2 since it is very small) from C4 of 0502425 */
        }

        /* Old implementation: */
        /* lens = ( X_022*X_000*d20[index_mu][l] +
                 Cgl2[index_mu]*2.*X_p000/sqrt5[l] * //If we want to revert to this old implementation we must also invert ssqrt5 since we have inverted definition.
                 (X_121*d11[index_mu][l] + X_132*d3m1[index_mu][l]) +
                 0.5 * Cgl2[index_mu] * Cgl2[index_mu] *
                 ( ( 2.*X_p022*X_p000+X_220*X_220 ) *
		 d20[index_mu][l] + X_220*X_242*d4m2[index_mu][l] ) );*/
        if (ppr->accurate_lensing == _FALSE_) {
          lens -= d20[index_mu][l];
        }
        resX *= lens;
        ksiX[index_mu] += resX;
      }

      if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {

        resp = fac1*(cl_ee[l]+cl_bb[l]);
        resm = fac1*(cl_ee[l]-cl_bb[l]);

        /* Order Cgl2 meaning we essentially only have the effect of exp(-l*(l+1)sigma2/2) */
        lensp = X_022*X_022*d22[index_mu][l];
        lensm = X_022*X_022*d2m2[index_mu][l];

        /* First order in Cgl2 from 54-55 of 0502425 */
        if (ppr->lensing_C2_order >= 1) {
          lensp += 2.*Cgl2[index_mu]*X_132*X_121*d31[index_mu][l];
          lensm += Cgl2[index_mu] * ( X_121*X_121*d1m1[index_mu][l] + X_132*X_132*d3m3[index_mu][l] );
          /* Second order in Cgl2 (same equations) */
          if (ppr->lensing_C2_order >= 2 ) {
            lensp += Cgl2[index_mu]*Cgl2[index_mu] * ( X_p022*X_p022*d22[index_mu][l] + X_242*X_220*d40[index_mu][l] );
            lensm += 0.5 * Cgl2[index_mu] * Cgl2[index_mu] * ( 2.*X_p022*X_p022*d2m2[index_mu][l] + X_220*X_220*d00[index_mu][l] + X_242*X_242*d4m4[index_mu][l] );
            /* Third and fourth order from 9.16-9.17 of 0601594 */
            if (ppr->lensing_C2_order >= 3 ) {
              lensp += X_000_square * x*x*x * (1/8.*d31[index_mu][l] + 1/24.*d5m1[index_mu][l]) ;
              lensm += X_000_square * x*x*x * (1/12.*d1m1[index_mu][l]+ 1/16.*d3m3[index_mu][l])  ;
              if (ppr->lensing_C2_order >= 4 ) {
                lensp += X_000_square * x*x*x*x * (1/64.*d22[index_mu][l] + 1/48.*d40[index_mu][l] + 1/192.*d6m2[index_mu][l]) ;
                lensm += X_000_square * x*x*x*x * (1/96.*d00[index_mu][l]+ 7/384.*d2m2[index_mu][l]+ 1/96.*d4m4[index_mu][l]);
              }
            }
          }
        }
        /* First order in Cgl from C2-C3 of 0502425. No coupling to Cgl2 because it is very very small. */
        if (ppr->lensing_C0_order == 1) {
          lensp += (2*X_022*X_p022*d22[index_mu][l]  + X_132*X_132*d33[index_mu][l] + X_121*X_121*d11[index_mu][l]) *Cgl[index_mu];
          lensm += (2*X_022*X_p022*d2m2[index_mu][l]  + 2*X_121*X_132*d3m1[index_mu][l]) *Cgl[index_mu];
        }

        /* Old implementation: */
	    /*
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
                    X_242*X_242*d4m4[index_mu][l] ) );*/
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
    return _SUCCESS_;
    );
  }

  class_finish_parallel();

  /** - compute lensed \f$ C_l\f$'s by integration */

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

  if ((ppr->lensing_C2_order > 2) || (ppr->lensing_C0_order > 0)) {
    free(d5m1);
    free(d5m3);
    free(d6m2);
    free(d33);
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

  ple->is_allocated = _TRUE_;

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

  ple->is_allocated = _FALSE_;

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

  int index_l;

  /** Integration by Gauss-Legendre quadrature. **/
  class_setup_parallel();

  for (index_l=0; index_l<ple->l_size; index_l++){
    class_run_parallel(=,
    double cle;
    int imu;
    cle=0;
    for (imu=0;imu<nmu;imu++) {
      cle += ksi[imu]*d00[imu][(int)ple->l[index_l]]*w8[imu]; /* loop could be optimized */
    }
    ple->cl_lens[index_l*ple->lt_size+ple->index_lt_tt]=cle*2.0*_PI_;
    return _SUCCESS_;
    );
  }

  class_finish_parallel();

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

  int index_l;

  /** Integration by Gauss-Legendre quadrature. **/
  class_setup_parallel();

  for (index_l=0; index_l < ple->l_size; index_l++){
    class_run_parallel(=,
    double clte;
    int imu;
    clte=0;
    for (imu=0;imu<nmu;imu++) {
      clte += ksiX[imu]*d20[imu][(int)ple->l[index_l]]*w8[imu]; /* loop could be optimized */
    }
    ple->cl_lens[index_l*ple->lt_size+ple->index_lt_te]=clte*2.0*_PI_;
    return _SUCCESS_;
    );
  }

  class_finish_parallel();
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

  int index_l;

  class_setup_parallel();
  /** Integration by Gauss-Legendre quadrature. **/
  for (index_l=0; index_l < ple->l_size; index_l++){
    class_run_parallel(=,
    double clp;
    double clm;
    int imu;
    clp=0; clm=0;
    for (imu=0;imu<nmu;imu++) {
      clp += ksip[imu]*d22[imu][(int)ple->l[index_l]]*w8[imu]; /* loop could be optimized */
      clm += ksim[imu]*d2m2[imu][(int)ple->l[index_l]]*w8[imu]; /* loop could be optimized */
    }
    ple->cl_lens[index_l*ple->lt_size+ple->index_lt_ee]=(clp+clm)*_PI_;
    ple->cl_lens[index_l*ple->lt_size+ple->index_lt_bb]=(clp-clm)*_PI_;
    return _SUCCESS_;
    );
  }
  class_finish_parallel();

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
 * This routine computes the d^l_{m1 m2}(cos beta) when the difference between m1 and m2 is an even number and when m1 >=0 and abs(m2) < m1.
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param dm1m2          Input/output: Result is stored here
 * @param error_message  Output: error message
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_dm1m2(
                  struct lensing * ple,
                  double * mu,
                  int num_mu,
                  int lmax,
                  int m1,
                  int m2,
                  double ** dm1m2
                  ) {
  double ll;
  int index_mu, l, lmin, i ;
  double *fac1, *fac2, *fac3, *fac4;
  double argsqrt=1.,pref;

  class_alloc(fac1,lmax*sizeof(double),ple->error_message);
  class_alloc(fac2,lmax*sizeof(double),ple->error_message);
  class_alloc(fac3,lmax*sizeof(double),ple->error_message);
  class_alloc(fac4,lmax*sizeof(double),ple->error_message);

  class_test(m1<0,ple->error_message,"You should not use this function with m1<0 and you used m1 = %d \n",m1);
  class_test(abs(m2)>m1,ple->error_message,"You should not use this function with abs(m2)>m1 and you used m1 = %d m2 = %d \n",m1,m2);

  lmin = m1;

  /** d^j_{jm}(beta), when (j-m) is even, is equal to sqrt((2j)!/(j-m)!/(j+m)!) * [(1+ cos(beta))/2]^((j+m)/2) *  [(1-cos(beta))/2]^((j-m)/2)
      We compute the prefactor which is common to all beta, and then we compute the powers involving the cos later when initializing the recurrence.
   */
  if (-m2==lmin) {
    pref=1.;
  }
  else {
    argsqrt=1.;
    for (i=1;i<=lmin+m2;i++)
      argsqrt /= i; /* = (2j)! */
    for (i=1;i<=lmin-m2;i++)
      argsqrt /= i; /* = 1/(j-m2)! */
    for (i=1;i<=2*lmin;i++)
      argsqrt *= i; /* = 1.(j+m2)! */
    pref = sqrt(argsqrt); /* Final sqrt. */
  }

  /* We must separate the case for d^l_00 for which lmin = 0, since there are some l^2/l which must be understood as being 0. */
  for (l=MAX(1,lmin);l<lmax;l++) {
    ll = (double) l;
    fac1[l] = sqrt((2*ll+3)*(2*ll+1)/((ll+1+m1)*(ll+1-m1)*(ll+1+m2)*(ll+1-m2)))*(ll+1);
    fac2[l] = -m1*m2/(ll*(ll+1));
    fac3[l] = sqrt((2*ll+3)/(2*ll-1)*(ll+m1)*(ll-m1)*(ll+m2)*(ll-m2)/((ll+1+m1)*(ll+1-m1)*(ll+1+m2)*(ll+1-m2)))*(ll+1)/ll;
    fac4[l] = sqrt(2./(2*ll+3));
  }
  if (lmin == 0) {
    fac1[0] = sqrt(3.);
    fac2[0] = 0.;
    fac3[0] = 0.;
    fac4[0] = sqrt(2./3.);
  }

  class_setup_parallel();
  for (index_mu=0;index_mu<num_mu;index_mu++) {
    class_run_parallel(=,
    int i;
    int l;
    double declare_list_of_variables_inside_parallel_region(mmu,m1,dl,dlp1,dlm1);
    //double mmu,m1,dl,dlp1,dlm1;

    mmu = mu[index_mu];
    for (l=0;l<lmin;l++){
      dm1m2[index_mu][l]=0;
    }
    m1=0.; /*l=lmin-1*/
    dl=sqrt((2.*lmin+1)/2.) * pref; /*l=lmin*/
    /* We now multiply by the [(1+ cos(beta))/2]^((j+m)/2) *  [(1-cos(beta))/2]^((j-m)/2) factor */
    for (i=1;i<=(lmin+m2)/2;i++)
      dl *= (1.+ mmu)/2.;
    for (i=1;i<=(lmin-m2)/2;i++)
      dl *= (1.- mmu)/2.;
    dm1m2[index_mu][lmin] = dl * sqrt(2./(2.*lmin+1));
    for (l=lmin;l<lmax;l++){
      /* sqrt((2l+1)/2)*dm1m2 recurrence, supposed to be more stable */
      dlp1 = fac1[l]*(mmu+fac2[l])*dl - fac3[l]*dlm1;
      dm1m2[index_mu][l+1] = dlp1 * fac4[l];
      dlm1 = dl;
      dl = dlp1;
    }
    return _SUCCESS_;
    );
  }
  class_finish_parallel();

  free(fac1); free(fac2); free(fac3); free(fac4);
  return _SUCCESS_;
}
