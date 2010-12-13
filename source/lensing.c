/** @file lensing.c Documented lensing module
 *
 * Julien Lesgourgues and Simon Prunet, 6.12.2010    
 *
 * This module computes the lensed temperature and polarization
 * anisotropy power spectra \f$ C_l^{X}, P(k), ... \f$'s given the
 * unlensed temperature, polarization and lensing potential spectra.
 *
 * The following functions can be called from other modules:
 *
 * -# lensing_init() at the beginning (but after spectra_init())
 * -# lensing_cl_at_l() at any time for computing Cl_lensed at any l
 * -# lensing_free() at the end
 */

#include "lensing.h"

/** 
 * Anisotropy power spectra C_l's for all types, modes and initial conditions.
 * SO FAR: ONLY SCALAR ADIABATIC 
 *
 * This routine evaluates all the lensed C_l's at a given value of l by
 * picking it in the pre-computed table. When relevant, it also
 * sums over all initial conditions for each mode, and over all modes.
 * 
 * This function can be called from whatever module at whatever time,
 * provided that lensing_init() has been called before, and
 * lensing_free() has not been called yet.
 *
 * @param ple        Input : pointer to lensing structure
 * @param l          Input : multipole number 
 * @param cl_lensed  Output: lensed C_l's for all types (TT, TE, EE, etc..)
 * @return the error status
 */

int lensing_cl_at_l(
		    struct lensing * ple,
		    int l,
		    double * cl_lensed    /* array with argument cl_lensed[index_ct] (must be already allocated) */
		    ) {
  int index_lt;

  class_test(l > ple->l_lensed_max,
	     "you asked for lensed Cls at l=%d, they were computed only up to l=%d, you should increase l_max_scalars or decrease the precision parameter delta_l_max",
	     ple->error_message);

  for (index_lt=0; index_lt<ple->lt_size; index_lt++) {
    cl_lensed[index_lt]=ple->cl_lensed[l*ple->lt_size+index_lt];
  }
  
  return _SUCCESS_;
}

/**
 * This routine initializes the lensing structure (in particular, 
 * computes table of lensed anisotropy spectra \f$ C_l^{X} \f$)
 * 
 * @param ppt Input : pointer to perturbation structure
 * @param psp Input : pointer to spectra structure
 * @param ple Output: pointer to initialized lensing structure
 * @return the error status
 */

int lensing_init(
		 struct precision * ppr,
		 struct perturbs * ppt,
		 struct spectra * psp,
		 struct lensing * ple
		 ) {

  /** local variables */

  double * mu; /* mu[index_mu]: discretized values of mu
		    between -1 and 1, roots of Legendre polynomial */
  double * w8; /* Corresponding Gauss-Legendre quadrature weights */
  
  double ** d00;  /* dmn[index_mu][index_l] */
  double ** d11;
  double ** d2m2; 
  double ** d22;  
  double ** d20;  
  double ** d1m1; 
  double ** d31;  
  double ** d40; 
  double ** d3m1;
  double ** d3m3;
  double ** d4m2;
  double ** d4m4; 

  double * Cgl;   /* Cgl[index_mu] */
  double * Cgl2;  /* Cgl2[index_mu] */
  double * sigma2; /* sigma[index_mu] */
  
  double * ksi;  /* ksi[index_mu] */
  double * ksiX;  /* ksiX[index_mu] */
  double * ksip;  /* ksip[index_mu] */
  double * ksim;  /* ksim[index_mu] */

  double ** X000;  /* Ximn[index_mu][index_l] */ 
  double ** Xp000;
  double ** X220; 
  double ** X022; 
  double ** X132; 
  double ** X121; 
  double ** Xp022; 
  double ** X242; 

  int num_mu,index_mu;
  int l;
  double ll;
  double * cl_unlensed;  /* cl_unlensed[index_ct] */
  double * cl_tt; /* unlensed  cl, to be filled to avoid repeated calls to spectra_cl_at_l */
  double * cl_te; /* unlensed  cl, to be filled to avoid repeated calls to spectra_cl_at_l */
  double * cl_ee; /* unlensed  cl, to be filled to avoid repeated calls to spectra_cl_at_l */
  double * cl_bb; /* unlensed  cl, to be filled to avoid repeated calls to spectra_cl_at_l */
  double * cl_pp; /* potential cl, to be filled to avoid repeated calls to spectra_cl_at_l */  
  double ** junk1=NULL, ** junk2=NULL;


  /** Summary: */

  /** - check that we really want to compute at least one spectrum */

  if (ple->has_lensed_cls == _FALSE_) {
    if (ple->lensing_verbose > 0)
      printf("No lensing requested. Lensing module skipped.\n");
    return _SUCCESS_;
  }
  else {
    if (ple->lensing_verbose > 0)
      printf("Computing lensed spectra\n");
  }
  
  /** - initialize indices and allocate some of the arrays in the 
      spectra structure */
  
  class_call(lensing_indices(ppr,ppt,psp,ple),
	     ple->error_message,
	     ple->error_message);

  /** put here all precision variables; will be stored later in precision structure */
  /** Last element in mu will be for mu=1, needed for sigma2 
   The rest will be chosen as roots of a Gauss-Legendre quadrature **/
  
  num_mu=(ple->l_unlensed_max+ppr->num_mu_minus_lmax); /* Must be even ?? CHECK */
  num_mu += num_mu%2; /* Force it to be even */ 
  /** - allocate array of mu values, as well as quadrature weights */

  class_alloc(mu,
              num_mu*sizeof(double),
              ple->error_message);
  /* Reserve last element of mu for mu=1, needed for sigma2 */
  mu[num_mu-1] = 1.0;
  
  class_alloc(w8,
              (num_mu-1)*sizeof(double),
              ple->error_message);

  class_call(lensing_gauss_legendre(mu,w8,num_mu-1),
             ple->error_message,
             ple->error_message); 

  /** - compute d^l_mm'(mu) */

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
  }
	
  for (index_mu=0; index_mu<num_mu; index_mu++) {

    class_alloc(d00[index_mu],
		(ple->l_unlensed_max+1)*sizeof(double),
		ple->error_message);  

    class_alloc(d11[index_mu],
		(ple->l_unlensed_max+1)*sizeof(double),
		ple->error_message);  

    class_alloc(d1m1[index_mu],
		(ple->l_unlensed_max+1)*sizeof(double),
		ple->error_message);  

    class_alloc(d2m2[index_mu],
		(ple->l_unlensed_max+1)*sizeof(double),
		ple->error_message);  
  }
  
  if (ple->has_te==_TRUE_) {
    for (index_mu=0; index_mu<num_mu; index_mu++) {
      class_alloc(d20[index_mu],
                  (ple->l_unlensed_max+1)*sizeof(double),
                  ple->error_message);

      class_alloc(d3m1[index_mu],
                  (ple->l_unlensed_max+1)*sizeof(double),
                  ple->error_message);
      
      class_alloc(d4m2[index_mu],
                  (ple->l_unlensed_max+1)*sizeof(double),
                  ple->error_message);
    }
  }

  if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {

    for (index_mu=0; index_mu<num_mu; index_mu++) {
      class_alloc(d22[index_mu],
                  (ple->l_unlensed_max+1)*sizeof(double),
                  ple->error_message);
      
      class_alloc(d31[index_mu],
                  (ple->l_unlensed_max+1)*sizeof(double),
                  ple->error_message);

      class_alloc(d3m3[index_mu],
                  (ple->l_unlensed_max+1)*sizeof(double),
                  ple->error_message);

      class_alloc(d40[index_mu],
                  (ple->l_unlensed_max+1)*sizeof(double),
                  ple->error_message);

      class_alloc(d4m4[index_mu],
                  (ple->l_unlensed_max+1)*sizeof(double),
                  ple->error_message);
    }    
  } 

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
  
  /** - compute Cgl(mu), Cgl2(mu) and sigma2(mu) */

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
	      psp->ct_size*sizeof(double),
	      ple->error_message);


  /** Locally store unlensed temperature cl_tt and potential cl_pp spectra **/
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

  for (l=2; l<=ple->l_unlensed_max; l++) {
    class_call(spectra_cl_at_l(psp,l,cl_unlensed,junk1,junk2),
               psp->error_message,
               ple->error_message);
    cl_tt[l] = cl_unlensed[psp->index_ct_tt];
    cl_pp[l] = cl_unlensed[psp->index_ct_pp];
    if (ple->has_te==_TRUE_) {
      cl_te[l] = cl_unlensed[psp->index_ct_te];
    }
    if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
      cl_ee[l] = cl_unlensed[psp->index_ct_ee];    
      cl_bb[l] = cl_unlensed[psp->index_ct_bb];
    }
  }

  /** Compute sigma2(mu) and Cgl2(mu) **/

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
  
  /** - compute X000(mu), X'000(mu), X220 and other Ximn */
  /** Zero separation is not used from now on, hence num_mu-1 **/
  class_alloc(X000,
              (num_mu-1)*sizeof(double*),
              ple->error_message);
  class_alloc(Xp000,
              (num_mu-1)*sizeof(double*),
              ple->error_message);
  class_alloc(X220,
              (num_mu-1)*sizeof(double*),
              ple->error_message);
  
  if (ple->has_te==_TRUE_ || ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {

    class_alloc(X022,
                (num_mu-1)*sizeof(double*),
                ple->error_message);
    class_alloc(Xp022,
                (num_mu-1)*sizeof(double*),
                ple->error_message);
    class_alloc(X121,
                (num_mu-1)*sizeof(double*),
                ple->error_message);
    class_alloc(X132,
                (num_mu-1)*sizeof(double*),
                ple->error_message);
    class_alloc(X242,
                (num_mu-1)*sizeof(double*),
                ple->error_message);
  }
  
  for (index_mu=0; index_mu<num_mu-1; index_mu++) {
    
    class_alloc(X000[index_mu],
                (ple->l_unlensed_max+1)*sizeof(double),
                ple->error_message);  
    class_alloc(Xp000[index_mu],
                (ple->l_unlensed_max+1)*sizeof(double),
                ple->error_message);  
    class_alloc(X220[index_mu],
                (ple->l_unlensed_max+1)*sizeof(double),
                ple->error_message);  
  }
  
  if (ple->has_te==_TRUE_ || ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {

    for (index_mu=0; index_mu<num_mu-1; index_mu++) {
      class_alloc(X022[index_mu],
                  (ple->l_unlensed_max+1)*sizeof(double),
                  ple->error_message);  
      class_alloc(Xp022[index_mu],
                  (ple->l_unlensed_max+1)*sizeof(double),
                  ple->error_message);  
      class_alloc(X121[index_mu],
                  (ple->l_unlensed_max+1)*sizeof(double),
                  ple->error_message);  
      class_alloc(X132[index_mu],
                  (ple->l_unlensed_max+1)*sizeof(double),
                  ple->error_message);  
      class_alloc(X242[index_mu],
                  (ple->l_unlensed_max+1)*sizeof(double),
                  ple->error_message);
    }
  }
  
  class_call(lensing_X000(mu,num_mu-1,ple->l_unlensed_max,sigma2,X000),
             ple->error_message,
             ple->error_message);
  class_call(lensing_Xp000(mu,num_mu-1,ple->l_unlensed_max,sigma2,Xp000),
             ple->error_message,
             ple->error_message);
  class_call(lensing_X220(mu,num_mu-1,ple->l_unlensed_max,sigma2,X220),
             ple->error_message,
             ple->error_message);
  
  if (ple->has_te==_TRUE_ || ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {

    class_call(lensing_X022(mu,num_mu-1,ple->l_unlensed_max,sigma2,X022),
               ple->error_message,
               ple->error_message);

    class_call(lensing_Xp022(mu,num_mu-1,ple->l_unlensed_max,sigma2,Xp022),
               ple->error_message,
               ple->error_message);

    class_call(lensing_X121(mu,num_mu-1,ple->l_unlensed_max,sigma2,X121),
               ple->error_message,
               ple->error_message);

    class_call(lensing_X132(mu,num_mu-1,ple->l_unlensed_max,sigma2,X132),
               ple->error_message,
               ple->error_message);

    class_call(lensing_X242(mu,num_mu-1,ple->l_unlensed_max,sigma2,X242),
               ple->error_message,
               ple->error_message);
  }
  

  /** - compute ksi, ksi+, ksi-, ksiX */
  
  /** ksi is for TT **/
  class_alloc(ksi,
              (num_mu-1)*sizeof(double),
              ple->error_message);
  {
    double res;
    for (index_mu=0;index_mu<num_mu-1;index_mu++) {
      ksi[index_mu]=0;
      for (l=2;l<=ple->l_unlensed_max;l++) {
        ll = (double) l;
        res = (2*ll+1)/(4.*_PI_)*cl_tt[l];
        /* res *= d00[l][index_mu]; */ /*DEBUG: TO BE REMOVED */
        
        res *= (X000[index_mu][l]*X000[index_mu][l]*d00[index_mu][l] +
                Xp000[index_mu][l]*Xp000[index_mu][l]*d1m1[index_mu][l]
                *Cgl2[index_mu]*8./(ll*(ll+1)) +
                (Xp000[index_mu][l]*Xp000[index_mu][l]*d00[index_mu][l] +
                 X220[index_mu][l]*X220[index_mu][l]*d2m2[index_mu][l])
                *Cgl2[index_mu]*Cgl2[index_mu]);
        
        ksi[index_mu] += res;
      }
    }
  }
  
  /** ksiX is for TE **/
  if (ple->has_te==_TRUE_) {
    class_alloc(ksiX,
                (num_mu-1)*sizeof(double),
                ple->error_message);
    {
      double res;
      double *sqllp1;
      class_alloc(sqllp1,
                  (ple->l_unlensed_max+1)*sizeof(double),
                  ple->error_message);
      for (l=2;l<=ple->l_unlensed_max;l++) {
        ll=(double)l;
        sqllp1[l] = sqrt(ll*(ll+1));
      }
      
      for (index_mu=0;index_mu<num_mu-1;index_mu++) {
        ksiX[index_mu]=0;
        for (l=2;l<=ple->l_unlensed_max;l++) {
          ll = (double)l;
          res = (2*ll+1)/(4.*_PI_)*cl_te[l];
          res *= ( X022[index_mu][l]*X000[index_mu][l]*d20[index_mu][l] +
                   Cgl2[index_mu]*2.*Xp000[index_mu][l]/sqllp1[l] *
                   (X121[index_mu][l]*d11[index_mu][l] + X132[index_mu][l]*d3m1[index_mu][l]) +
                   0.5 * Cgl2[index_mu] * Cgl2[index_mu] *
                   ( ( 2.*Xp022[index_mu][l]*Xp000[index_mu][l]+X220[index_mu][l]*X220[index_mu][l] ) *
                    d20[index_mu][l] + X220[index_mu][l]*X242[index_mu][l]*d4m2[index_mu][l] ) );
          ksiX[index_mu] += res;
        }
      }
      free(sqllp1);
    }
  }

  /** ksip, ksim for EE, BB **/
  if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
    class_alloc(ksip,
		(num_mu-1)*sizeof(double),
		ple->error_message);
    class_alloc(ksim,
		(num_mu-1)*sizeof(double),
		ple->error_message);
    {
      double resp, resm;
      for (index_mu=0;index_mu<num_mu-1;index_mu++) {
	ksip[index_mu]=0;
	ksim[index_mu]=0;
	for (l=2;l<=ple->l_unlensed_max;l++) {
	  ll = (double) l;
	  resp = (2*ll+1)/(4.*_PI_)*(cl_ee[l]+cl_bb[l]);
	  resm = (2*ll+1)/(4.*_PI_)*(cl_ee[l]-cl_bb[l]);
	  
	  resp *= ( X022[index_mu][l]*X022[index_mu][l]*d22[index_mu][l] +
		    2.*Cgl2[index_mu]*X132[index_mu][l]*X121[index_mu][l]*d31[index_mu][l] +
		    Cgl2[index_mu]*Cgl2[index_mu] *
		    ( Xp022[index_mu][l]*Xp022[index_mu][l]*d22[index_mu][l] +
		      X242[index_mu][l]*X220[index_mu][l]*d40[index_mu][l] ) );

	  resm *= ( X022[index_mu][l]*X022[index_mu][l]*d2m2[index_mu][l] +
		    Cgl2[index_mu] *
		    ( X121[index_mu][l]*X121[index_mu][l]*d1m1[index_mu][l] +
		      X132[index_mu][l]*X132[index_mu][l]*d3m3[index_mu][l] ) +
		    0.5 * Cgl2[index_mu] * Cgl2[index_mu] *
		    ( 2.*Xp022[index_mu][l]*Xp022[index_mu][l]*d2m2[index_mu][l] +
		      X220[index_mu][l]*X220[index_mu][l]*d00[index_mu][l] +
		      X242[index_mu][l]*X242[index_mu][l]*d4m4[index_mu][l] ) );

	  ksip[index_mu] += resp;
	  ksim[index_mu] += resm;
	}
      }
    }
  }

  /** - compute lensed Cls by integration */
  class_call(lensing_lensed_cl_tt(ksi,d00,w8,num_mu-1,ple),
             ple->error_message,
             ple->error_message);

  if (ple->has_te==_TRUE_) {
    class_call(lensing_lensed_cl_te(ksiX,d20,w8,num_mu-1,ple),
               ple->error_message,
               ple->error_message);
  }
  
  if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {

    class_call(lensing_lensed_cl_ee_bb(ksip,ksim,d22,d2m2,w8,num_mu-1,ple),
	       ple->error_message,
	       ple->error_message);
  } 

  /** DEBUG **/
  /*
  FILE *fpp;
  fpp=fopen("toto.txt","w");
  for (l=2; l<=ple->l_lensed_max; l++) {
    fprintf(fpp,"%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",l,
            cl_tt[l],
            ple->cl_lensed[l*ple->lt_size+ple->index_lt_tt],
            cl_te[l],
            ple->cl_lensed[l*ple->lt_size+ple->index_lt_te],
	    cl_ee[l],
	    ple->cl_lensed[l*ple->lt_size+ple->index_lt_ee],
	    cl_bb[l],
	    ple->cl_lensed[l*ple->lt_size+ple->index_lt_bb]
	    );
  }
  fclose(fpp);
  */
  /** Free lots of stuff **/
  for (index_mu=0; index_mu<num_mu-1; index_mu++) {
    free(d00[index_mu]);
    free(d11[index_mu]);
    free(d1m1[index_mu]);
    free(d2m2[index_mu]);
    free(X000[index_mu]);
    free(Xp000[index_mu]);
    free(X220[index_mu]);
  }
  /** Free zero-separation vectors **/
  free(d00[num_mu-1]);
  free(d11[num_mu-1]);
  free(d1m1[num_mu-1]);
  free(d2m2[num_mu-1]);

  /* Free the remaining vectors */
  free(d00);
  free(d11);
  free(d1m1);
  free(d2m2);

  free(X000);
  free(Xp000);
  free(X220);

  free(ksi);
  free(Cgl);
  free(Cgl2);

  free(mu);
  free(w8);

  free(cl_unlensed);
  free(cl_tt);
  free(cl_pp);
  /** Exits **/
  
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

    free(ple->cl_lensed);
  }

  return _SUCCESS_;
 
}

/**
 * This routine defines indices and allocates tables in the lensing structure 
 *
 * @param ppt  Input       : pointer to perturbation structure
 * @param psp  Input       : pointer to spectra structure 
 * @param ple  Input/output: pointer to lensing structure 
 * @return the error status
 */

int lensing_indices(
		    struct precision * ppr,
		    struct perturbs * ppt,
		    struct spectra * psp,
		    struct lensing * ple
		    ){
  
  int l;
  double ** junk1=NULL;
  double ** junk2=NULL;
  
  /* indices of all Cl types (lensed and unlensed) */
  
  if (psp->has_tt == _TRUE_) {
    ple->has_tt = _TRUE_;
    ple->index_lt_tt=psp->index_ct_tt;      
  }
  else {
    ple->has_tt = _FALSE_;
  }

  if (psp->has_ee == _TRUE_) {
    ple->has_ee = _TRUE_;
    ple->index_lt_ee=psp->index_ct_ee;      
  }
  else {
    ple->has_ee = _FALSE_;
  }

  if (psp->has_te == _TRUE_) {
    ple->has_te = _TRUE_;
    ple->index_lt_te=psp->index_ct_te;      
  }
  else {
    ple->has_te = _FALSE_;
  }

  if (psp->has_bb == _TRUE_) {
    ple->has_bb = _TRUE_;
    ple->index_lt_bb=psp->index_ct_bb;      
  }
  else {
    ple->has_bb = _FALSE_;
  }

  if (psp->has_pp == _TRUE_) {
    ple->has_pp = _TRUE_;
    ple->index_lt_te=psp->index_ct_pp;      
  }
  else {
    ple->has_pp = _FALSE_;
  }

  if (psp->has_tp == _TRUE_) {
    ple->has_tp = _TRUE_;
    ple->index_lt_tp=psp->index_ct_tp;      
  }
  else {
    ple->has_tp = _FALSE_;
  }

  ple->lt_size = psp->ct_size;

  /* number of multipoles */

  ple->l_unlensed_max = psp->l_max[ppt->index_md_scalars];

  ple->l_lensed_max = ple->l_unlensed_max - ppr->delta_l_max;
  
  /* allocate table where results will be stored */
  
  class_alloc(ple->cl_lensed,
	      (ple->l_lensed_max+1)*ple->lt_size*sizeof(double),
	      ple->error_message);
  
  /* fill with unlensed cls */
  /* Should be removed when polarized lensed cls are in place */
  for (l=2; l<=ple->l_lensed_max; l++) { 
    
    class_call(spectra_cl_at_l(psp,l,&(ple->cl_lensed[l*ple->lt_size]),junk1,junk2),
	       psp->error_message,
	       ple->error_message);
    
  }

  return _SUCCESS_;
  
}
  
/**
 * This routine computes the lensed power spectra by Gaussian quadrature 
 *
 * @param ksi  Input       : Lensed correlation function (ksi[index_mu])
 * @param d00  Input       : Legendre polynomials (d^l_{00}[l][index_mu]) 
 * @param w8   Input       : Legendre quadrature weights (w8[index_mu])
 * @param nmu  Input       : Number of quadrature points (0<=index_mu<=nmu)
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
  int l, imu;
  /** Integration by Gauss-Legendre quadrature **/
  for(l=2;l<=ple->l_lensed_max;l++){
    cle=0;
    for (imu=0;imu<nmu;imu++) {
      cle += ksi[imu]*d00[imu][l]*w8[imu]; /* loop could be optimized */
    }
    ple->cl_lensed[l*ple->lt_size+ple->index_lt_tt]=cle*2.0*_PI_;
  }
  return _SUCCESS_;
}

/**
 * This routine computes the lensed power spectra by Gaussian quadrature 
 *
 * @param ksiX Input       : Lensed correlation function (ksiX[index_mu])
 * @param d20  Input       : Wigner d-function (d^l_{20}[l][index_mu]) 
 * @param w8   Input       : Legendre quadrature weights (w8[index_mu])
 * @param nmu  Input       : Number of quadrature points (0<=index_mu<=nmu)
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
  int l, imu;
  /** Integration by Gauss-Legendre quadrature **/
  for(l=2;l<=ple->l_lensed_max;l++){
    clte=0;
    for (imu=0;imu<nmu;imu++) {
      clte += ksiX[imu]*d20[imu][l]*w8[imu]; /* loop could be optimized */
    }
    ple->cl_lensed[l*ple->lt_size+ple->index_lt_te]=clte*2.0*_PI_;
  }
  return _SUCCESS_;
}

/**
 * This routine computes the lensed power spectra by Gaussian quadrature 
 *
 * @param ksip Input       : Lensed correlation function (ksi+[index_mu])
 * @param ksim Input       : Lensed correlation function (ksi-[index_mu])
 * @param d22  Input       : Wigner d-function (d^l_{22}[l][index_mu]) 
 * @param d2m2 Input       : Wigner d-function (d^l_{2-2}[l][index_mu]) 
 * @param w8   Input       : Legendre quadrature weights (w8[index_mu])
 * @param nmu  Input       : Number of quadrature points (0<=index_mu<=nmu)
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
  int l, imu;
  /** Integration by Gauss-Legendre quadrature **/
  for(l=2;l<=ple->l_lensed_max;l++){
    clp=0; clm=0; 
    for (imu=0;imu<nmu;imu++) {
      clp += ksip[imu]*d22[imu][l]*w8[imu]; /* loop could be optimized */
      clm += ksim[imu]*d2m2[imu][l]*w8[imu]; /* loop could be optimized */
    }
    ple->cl_lensed[l*ple->lt_size+ple->index_lt_ee]=(clp+clm)*_PI_;
    ple->cl_lensed[l*ple->lt_size+ple->index_lt_bb]=(clp-clm)*_PI_;
  }
  return _SUCCESS_;
}

/**
 * This routine computes the X000 term
 *
 * @param mu     Input       : Vector of cos(beta) values
 * @param num_mu Input       : Number of cos(beta) values
 * @param lmax   Input       : maximum multipole
 * @param sigma2 Input       : Vector of sigma2(mu) values
 * @param X000   Input/output: Result is stored here
 
 **/

int lensing_X000(
        double * mu,
        int num_mu,
        int lmax,
        double * sigma2,
        double ** X000
        ) {
  int index_mu, l;
  double ll;
  for (index_mu=0;index_mu<num_mu;index_mu++) {
    for (l=0;l<=lmax;l++) {
      ll = (double) l;
      X000[index_mu][l]=exp(-ll*(ll+1)*sigma2[index_mu]/4.);
    }
  }
  return _SUCCESS_;
}

/**
 * This routine computes the Xp000 term
 *
 * @param mu     Input       : Vector of cos(beta) values
 * @param num_mu Input       : Number of cos(beta) values
 * @param lmax   Input       : maximum multipole
 * @param sigma2 Input       : Vector of sigma2(mu) values
 * @param Xp000  Input/output: Result is stored here
 
 **/

int lensing_Xp000(
                 double * mu,
                 int num_mu,
                 int lmax,
                 double * sigma2,
                 double ** Xp000
                 ) {
  int index_mu, l;
  double ll;
  double *fac;
  ErrorMsg erreur;
  class_alloc(fac,(lmax+1)*sizeof(double),erreur);
  for (l=0;l<=lmax;l++) {
    ll = (double) l;
    fac[l]=ll*(ll+1)/4.;
  }
  for (index_mu=0;index_mu<num_mu;index_mu++) {
    for (l=0;l<=lmax;l++) {
      ll = (double) l;
      Xp000[index_mu][l]=-fac[l]*exp(-fac[l]*sigma2[index_mu]);
    }
  }
  free(fac);
  return _SUCCESS_;
}

/**
 * This routine computes the X220 term
 *
 * @param mu     Input       : Vector of cos(beta) values
 * @param num_mu Input       : Number of cos(beta) values
 * @param lmax   Input       : maximum multipole
 * @param sigma2 Input       : Vector of sigma2(mu) values
 * @param X220   Input/output: Result is stored here
 
 **/

int lensing_X220(
                 double * mu,
                 int num_mu,
                 int lmax,
                 double * sigma2,
                 double ** X220
                 ) {
  int index_mu, l;
  double ll;
  double *fac1, *fac2;
  ErrorMsg erreur;
  class_alloc(fac1,(lmax+1)*sizeof(double),erreur);
  class_alloc(fac2,(lmax+1)*sizeof(double),erreur);
  for (l=2; l<=lmax; l++) {
    ll = (double) l;
    fac1[l] = 0.25*sqrt((ll+2)*(ll+1)*ll*(ll-1));
    fac2[l] = ll*(ll+1)/4.;
  }
  for (index_mu=0;index_mu<num_mu;index_mu++) {
    X220[index_mu][0]=0;
    X220[index_mu][1]=0;
    for (l=2;l<=lmax;l++) {
      ll = (double) l;
      X220[index_mu][l]=fac1[l] * exp(-fac2[l]*sigma2[index_mu]);
    }
  }
  free(fac1); free(fac2);
  return _SUCCESS_;
}

/**
 * This routine computes the X022 term
 *
 * @param mu     Input       : Vector of cos(beta) values
 * @param num_mu Input       : Number of cos(beta) values
 * @param lmax   Input       : maximum multipole
 * @param sigma2 Input       : Vector of sigma2(mu) values
 * @param X022   Input/output: Result is stored here
 
 **/

int lensing_X022(
                 double * mu,
                 int num_mu,
                 int lmax,
                 double * sigma2,
                 double ** X022
                 ) {
  int index_mu, l;
  double ll;
  double *fac;
  ErrorMsg erreur;
  class_alloc(fac,(lmax+1)*sizeof(double),erreur);
  for (l=2; l<=lmax; l++) {
    ll = (double) l;
    fac[l] = (ll*(ll+1)-4.)/4.;
  }
  for (index_mu=0;index_mu<num_mu;index_mu++) {
    X022[index_mu][0]=0;
    X022[index_mu][1]=0;
    for (l=2;l<=lmax;l++) {
      ll = (double) l;
      X022[index_mu][l]=exp(-fac[l]*sigma2[index_mu]);
    }
  }
  free(fac);
  return _SUCCESS_;
}

/**
 * This routine computes the Xp022 term
 *
 * @param mu     Input       : Vector of cos(beta) values
 * @param num_mu Input       : Number of cos(beta) values
 * @param lmax   Input       : maximum multipole
 * @param sigma2 Input       : Vector of sigma2(mu) values
 * @param Xp022  Input/output: Result is stored here
 
 **/

int lensing_Xp022(
                 double * mu,
                 int num_mu,
                 int lmax,
                 double * sigma2,
                 double ** Xp022
                 ) {
  int index_mu, l;
  double ll;
  double *fac;
  ErrorMsg erreur;
  class_alloc(fac,(lmax+1)*sizeof(double),erreur);
  for (l=2; l<=lmax; l++) {
    ll = (double) l;
    fac[l] = (ll*(ll+1)-4.)/4.;
  }
  for (index_mu=0;index_mu<num_mu;index_mu++) {
    Xp022[index_mu][0]=0;
    Xp022[index_mu][1]=0;
    for (l=2;l<=lmax;l++) {
      ll = (double) l;
      Xp022[index_mu][l]= -fac[l]*exp(-fac[l]*sigma2[index_mu]);
    }
  }
  free(fac);
  return _SUCCESS_;
}

/**
 * This routine computes the X121 term
 *
 * @param mu     Input       : Vector of cos(beta) values
 * @param num_mu Input       : Number of cos(beta) values
 * @param lmax   Input       : maximum multipole
 * @param sigma2 Input       : Vector of sigma2(mu) values
 * @param X121   Input/output: Result is stored here
 
 **/

int lensing_X121(
                 double * mu,
                 int num_mu,
                 int lmax,
                 double * sigma2,
                 double ** X121
                 ) {
  int index_mu, l;
  double ll;
  double *fac1, *fac2;
  ErrorMsg erreur;
  class_alloc(fac1,(lmax+1)*sizeof(double),erreur);
  class_alloc(fac2,(lmax+1)*sizeof(double),erreur);
  for (l=2; l<=lmax; l++) {
    ll = (double) l;
    fac1[l] = 0.5*sqrt((ll+2)*(ll-1));
    fac2[l] = (ll*(ll+1)-8./3.)/4.;
  }
  for (index_mu=0;index_mu<num_mu;index_mu++) {
    X121[index_mu][0]=0;
    X121[index_mu][1]=0;
    for (l=2;l<=lmax;l++) {
      ll = (double) l;
      X121[index_mu][l]= -fac1[l] * exp(-fac2[l]*sigma2[index_mu]);
    }
  }
  free(fac1); free(fac2);
  return _SUCCESS_;
}

/**
 * This routine computes the X132 term
 *
 * @param mu     Input       : Vector of cos(beta) values
 * @param num_mu Input       : Number of cos(beta) values
 * @param lmax   Input       : maximum multipole
 * @param sigma2 Input       : Vector of sigma2(mu) values
 * @param X132   Input/output: Result is stored here
 
 **/

int lensing_X132(
                 double * mu,
                 int num_mu,
                 int lmax,
                 double * sigma2,
                 double ** X132
                 ) {
  int index_mu, l;
  double ll;
  double *fac1, *fac2;
  ErrorMsg erreur;
  class_alloc(fac1,(lmax+1)*sizeof(double),erreur);
  class_alloc(fac2,(lmax+1)*sizeof(double),erreur);
  for (l=3; l<=lmax; l++) {
    ll = (double) l;
    fac1[l] = 0.5*sqrt((ll+3)*(ll-2));
    fac2[l] = (ll*(ll+1)-20./3.)/4.;
  }
  for (index_mu=0;index_mu<num_mu;index_mu++) {
    X132[index_mu][0]=0;
    X132[index_mu][1]=0;
    X132[index_mu][2]=0;
    for (l=3;l<=lmax;l++) {
      ll = (double) l;
      X132[index_mu][l]= -fac1[l] * exp(-fac2[l]*sigma2[index_mu]);
    }
  }
  free(fac1); free(fac2);
  return _SUCCESS_;
}

/**
 * This routine computes the X242 term
 *
 * @param mu     Input       : Vector of cos(beta) values
 * @param num_mu Input       : Number of cos(beta) values
 * @param lmax   Input       : maximum multipole
 * @param sigma2 Input       : Vector of sigma2(mu) values
 * @param X242   Input/output: Result is stored here
 
 **/

int lensing_X242(
                 double * mu,
                 int num_mu,
                 int lmax,
                 double * sigma2,
                 double ** X242
                 ) {
  int index_mu, l;
  double ll;
  double *fac1, *fac2;
  ErrorMsg erreur;
  class_alloc(fac1,(lmax+1)*sizeof(double),erreur);
  class_alloc(fac2,(lmax+1)*sizeof(double),erreur);
  for (l=4; l<=lmax; l++) {
    ll = (double) l;
    fac1[l] = 0.25*sqrt((ll+4)*(ll+3)*(ll-2.)*(ll-3));
    fac2[l] = (ll*(ll+1)-10.)/4.;
  }
  for (index_mu=0;index_mu<num_mu;index_mu++) {
    X242[index_mu][0]=0;
    X242[index_mu][1]=0;
    X242[index_mu][2]=0;
    X242[index_mu][3]=0;
    for (l=4;l<=lmax;l++) {
      ll = (double) l;
      X242[index_mu][l]=fac1[l] * exp(-fac2[l]*sigma2[index_mu]);
    }
  }
  free(fac1); free(fac2);
  return _SUCCESS_;
}

/**
 * This routine computes the d00 term
 *
 * @param mu     Input       : Vector of cos(beta) values
 * @param num_mu Input       : Number of cos(beta) values
 * @param lmax   Input       : maximum multipole
 * @param d00    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on sqrt((2l+1)/2) d^l_{mm'} for stability
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
 * @param mu     Input       : Vector of cos(beta) values
 * @param num_mu Input       : Number of cos(beta) values
 * @param lmax   Input       : maximum multipole
 * @param d11    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on sqrt((2l+1)/2) d^l_{mm'} for stability
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
 * @param mu     Input       : Vector of cos(beta) values
 * @param num_mu Input       : Number of cos(beta) values
 * @param lmax   Input       : maximum multipole
 * @param d1m1    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on sqrt((2l+1)/2) d^l_{mm'} for stability
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
 * @param mu     Input       : Vector of cos(beta) values
 * @param num_mu Input       : Number of cos(beta) values
 * @param lmax   Input       : maximum multipole
 * @param d2m2    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on sqrt((2l+1)/2) d^l_{mm'} for stability
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
 * @param mu     Input       : Vector of cos(beta) values
 * @param num_mu Input       : Number of cos(beta) values
 * @param lmax   Input       : maximum multipole
 * @param d22    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on sqrt((2l+1)/2) d^l_{mm'} for stability
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
 * @param mu     Input       : Vector of cos(beta) values
 * @param num_mu Input       : Number of cos(beta) values
 * @param lmax   Input       : maximum multipole
 * @param d20    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on sqrt((2l+1)/2) d^l_{mm'} for stability
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
 * @param mu     Input       : Vector of cos(beta) values
 * @param num_mu Input       : Number of cos(beta) values
 * @param lmax   Input       : maximum multipole
 * @param d31    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on sqrt((2l+1)/2) d^l_{mm'} for stability
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
 * @param mu     Input       : Vector of cos(beta) values
 * @param num_mu Input       : Number of cos(beta) values
 * @param lmax   Input       : maximum multipole
 * @param d3m1   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on sqrt((2l+1)/2) d^l_{mm'} for stability
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
 * @param mu     Input       : Vector of cos(beta) values
 * @param num_mu Input       : Number of cos(beta) values
 * @param lmax   Input       : maximum multipole
 * @param d3m3   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on sqrt((2l+1)/2) d^l_{mm'} for stability
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
 * @param mu     Input       : Vector of cos(beta) values
 * @param num_mu Input       : Number of cos(beta) values
 * @param lmax   Input       : maximum multipole
 * @param d40    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on sqrt((2l+1)/2) d^l_{mm'} for stability
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
 * @param mu     Input       : Vector of cos(beta) values
 * @param num_mu Input       : Number of cos(beta) values
 * @param lmax   Input       : maximum multipole
 * @param d4m2   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on sqrt((2l+1)/2) d^l_{mm'} for stability
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
 * @param mu     Input       : Vector of cos(beta) values
 * @param num_mu Input       : Number of cos(beta) values
 * @param lmax   Input       : maximum multipole
 * @param d4m4   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on sqrt((2l+1)/2) d^l_{mm'} for stability
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

/**
 * This routine computes the weights and abscissas of a Gauss-Legendre quadrature
 *
 * @param mu     Input/output: Vector of cos(beta) values
 * @param w8     Input/output: Vector of quadrature weights
 * @param n      Input       : Number of quadrature points
 *
 * From Numerical recipes 
 **/

int lensing_gauss_legendre(
                           double *mu,
                           double *w8,
                           int n) {
  
  int m,j,i;
  double z1,z,xm,xl,pp,p3,p2,p1;
  double x1,x2,EPS;
  x1 = -1.0;
  x2 = 1.0;
  EPS = 1.0e-11;
  
  m=(n+1)/2;
  xm=0.5*(x2+x1);
  xl=0.5*(x2-x1);
  for (i=1;i<=m;i++) {
    z=cos(_PI_*(i-0.25)/(n+0.5));
    do {
      p1=1.0;
      p2=0.0;
      for (j=1;j<=n;j++) {
        p3=p2;
        p2=p1;
        p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp=n*(z*p1-p2)/(z*z-1.0);
      z1=z;
      z=z1-p1/pp;
    } while (fabs(z-z1) > EPS);
    mu[i-1]=xm-xl*z;
    mu[n-i]=xm+xl*z;
    w8[i-1]=2.0*xl/((1.0-z*z)*pp*pp);
    w8[n-i]=w8[i-1];
    
  }
  return _SUCCESS_;
}

