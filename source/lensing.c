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
		 struct perturbs * ppt,
		 struct spectra * psp,
		 struct lensing * ple
		 ) {

  /** local variables */

  double * mu; /* mu[index_mu]: discretized values of mu
		    between 0 and pi */

  double ** d00;  /* dmn[index_l][index_mu] */
  double ** d11;
  double ** d2m2; 
  double ** d22;  
  double ** d20;  
  double ** d1m1; 
  double ** d31;  
  double ** d40;  
  double ** d3m3; 
  double ** d4m4; 

  double * Cgl;   /* Cgl[index_mu] */
  double * Cgl2;  /* Cgl2[index_mu] */
  double * sigma2; /* sigma[index_mu] */

  double ** X000;  /* Ximn[index_l][index_mu] */ 
  double ** Xp000;
  double ** X220; 
  double ** X022; 
  double ** X132; 
  double ** X121; 
  double ** Xp022; 
  double ** X242; 
  double ** X112; 

  int num_mu,index_mu;
  int l;
  double ll;
  double * cl_unlensed;  /* cl_unlensed[index_ct] */
  double ** junk1, ** junk2;

  /** put here all precision variables; will be stored later in precision structure */

  /** Last element in mu will be for mu=1, needed for sigma2 
      The rest will be chosen as roots of a Gauss-Legendre quadrature **/
  
  int num_mu=2001;
   

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
  
  class_call(lensing_indices(ppt,psp,ple),
	     ple->error_message,
	     ple->error_message);

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
	      (ple->l_max+1)*sizeof(double*),
	      ple->error_message);

  class_alloc(d11,
	      (ple->l_max+1)*sizeof(double*),
	      ple->error_message);
	
  class_alloc(d1m1,
	      (ple->l_max+1)*sizeof(double*),
	      ple->error_message);
	
  class_alloc(d2m2,
	      (ple->l_max+1)*sizeof(double*),
	      ple->error_message);
	
  for (l=2; l<=ple->l_max; l++) {

    class_alloc(d00[l],
		num_mu*sizeof(double),
		ple->error_message);  

    class_alloc(d11[l],
		num_mu*sizeof(double),
		ple->error_message);  

    class_alloc(d1m1[l],
		num_mu*sizeof(double),
		ple->error_message);  

    class_alloc(d2m2[l],
		num_mu*sizeof(double),
		ple->error_message);  
  }

  class_call(lensing_d00(mu,num_mu,d00),
	     ple->error_message,
	     ple->error_message);

  class_call(lensing_d11(mu,num_mu,d11),
	     ple->error_message,
	     ple->error_message);

  class_call(lensing_d1m1(mu,num_mu,d1m1),
	     ple->error_message,
	     ple->error_message);

  class_call(lensing_d2m2(mu,num_mu,d2m2),
	     ple->error_message,
	     ple->error_message);

  /** - compute Cgl(mu), Cgl2(mu) and sigma2(mu) */

  class_alloc(Cgl,
	      num_mu*sizeof(double),
	      ple->error_message);

  class_alloc(Cgl2,
	      num_mu*sizeof(double),
	      ple->error_message);

  class_alloc(sigma2,
	      num_mu*sizeof(double),
	      ple->error_message);

  class_alloc(cl_unlensed,
	      psp->ct_size*sizeof(double),
	      ple->error_message);

  for (index_mu=0; index_mu<num_mu; index_mu++) {

    Cgl[index_mu]=0;
    Cgl2[index_mu]=0;

    for (l=2; l<=ple->l_max; l++) {

      class_call(spectra_cl_at_l(psp,l,cl_unlensed,junk1,junk2),
		 psp->error_message,
		 ple->error_message);

      Cgl[index_mu] += (2.*l+1.)*l*(l+1.)*
	cl_unlensed[psp->index_ct_pp]*d11[l][index_mu];

      Cgl2[index_mu] += (2.*l+1.)*l*(l+1.)*
	cl_unlensed[psp->index_ct_pp]*dm11[l][index_mu];

    }

    Cgl[index_mu] /= 4.*_PI_;
    Cgl2[index_mu] /= 4.*_PI_;

  }

  for (index_mu=0; index_mu<num_mu; index_mu++) {

    sigma2[index_mu] = Cgl[0] - Cgl[index_mu];

  }

  /** - compute X000(mu), X'000(mu), X220 and other Ximn */
  class_alloc(X000,
              (ple->l_max+1)*sizeof(double*),
              ple->error_message);
  class_alloc(Xp000,
              (ple->l_max+1)*sizeof(double*),
              ple->error_message);
  class_alloc(X220,
              (ple->l_max+1)*sizeof(double*),
              ple->error_message);
  
  for (l=2; l<=ple->l_max; l++) {
    
    class_alloc(X000[l],
                num_mu*sizeof(double),
                ple->error_message);  
    class_alloc(Xp000[l],
                num_mu*sizeof(double),
                ple->error_message);  
    class_alloc(X220[l],
                num_mu*sizeof(double),
                ple->error_message);  
  }
  
  class_call(lensing_X000(mu,num_mu,sigma2,X000),
             ple->error_message,
             ple->error_message);
  class_call(lensing_Xp000(mu,num_mu,sigma2,Xp000),
             ple->error_message,
             ple->error_message);
  class_call(lensing_X220(mu,num_mu,sigma2,X220),
             ple->error_message,
             ple->error_message);
  
  

  /** - compute ksi, ksi+, ksi-, ksiX */
  class_alloc(ksi,
              num_mu*sizeof(double),
              ple->error_message);
  
  for (index_mu=0;index_mu<num_mu;index_mu++) {
    ksi[index_mu]=0;
    for (l=2;l<=ple->l_max;l++) {
      ll = (double) l;
      class_call(spectra_cl_at_l(psp,l,cl_unlensed,junk1,junk2),
                 psp->error_message,
                 ple->error_message);
      res = (2*ll+1)/(4.*_PI_);
      res *= (X000[l][index_mu]*X000[l][index_mu]*d00[l][index_mu] +
              Xp000[l][index_mu]*Xp000[l][index_mu]*d1m1[l][index_mu]
              *Cgl2[index_mu]*8./(ll*(ll+1)) +
              (Xp000[l][index_mu]*Xp000[l][index_mu]*d00[l][index_mu] +
               X220[l][index_mu]*X220[l][index_mu]*d2m2[l][index_mu])
              *Cgl2[index_mu]*Cgl2[index_mu])
      ksi[index_mu] += res;
    }
  }
  /** - compute lensed Cls by integration */
  class_call(lensing_lensed_cl(ple),
             ple->error_message,
             ple->error_message);
  

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
		    struct perturbs * ppt,
		    struct spectra * psp,
		    struct lensing * ple
		    ){
  
  int index_lt;
  
  /* types of C_l's relevant for lensing: TT, EE, TE */

  index_lt=0;

  if (psp->has_tt == _TRUE_) {
    ple->has_tt = _TRUE_;
    ple->index_lt_tt=index_lt;      
    index_lt++;
  }
  else {
    ple->has_tt = _FALSE_;
  }

  if (psp->has_ee == _TRUE_) {
    ple->has_ee = _TRUE_;
    ple->index_lt_ee=index_lt;      
    index_lt++;
  }
  else {
    ple->has_ee = _FALSE_;
  }

  if (psp->has_te == _TRUE_) {
    ple->has_te = _TRUE_;
    ple->index_lt_te=index_lt;      
    index_lt++;
  }
  else {
    ple->has_te = _FALSE_;
  }

  ple->lt_size = index_lt;

  /* number of multipoles */

  ple->l_max = psp->l_max[ppt->index_md_scalars];

  /* allocate table where results will be stored */

  class_alloc(ple->cl_lensed,
	      (ple->l_max+1)*ple->lt_size*sizeof(double),
	      ple->error_message);
  
  return _SUCCESS_;
  
}

int lensing_lensed_cl(
        double *ksi, 
        double **d00, 
        lensing * ple
        ) {
  
  double cle;
  cle =0;
  /* Integration to be replaced by a quadrature rule ... */
  
}


int lensing_X000(
        double * mu,
        int num_mu,
        double * sigma2,
        double * X000
        ) {
  int index_mu, l;
  double ll;
  for (l=2;l<=ple->l_max;l++) {
    ll = (double) l;
    for (index_mu=0;index_mu<num_mu;index_mu++) {
      X000[l][index_mu]=exp(-ll*(ll+1)*sigma2[index_mu]/4.);
    }
  }
  return _SUCCESS_;
}

int lensing_Xp000(
                 double * mu,
                 int num_mu,
                 double * sigma2,
                 double * Xp000
                 ) {
  int index_mu, l;
  double ll;
  for (l=2;l<=ple->l_max;l++) {
    ll = (double) l;
    for (index_mu=0;index_mu<num_mu;index_mu++) {
      Xp000[l][index_mu]=-ll*(ll+1)/4.*exp(-ll*(ll+1)*sigma2[index_mu]/4.);
    }
  }
  return _SUCCESS_;
}

int lensing_X220(
                 double * mu,
                 int num_mu,
                 double * sigma2,
                 double * X220
                 ) {
  int index_mu, l;
  double ll;
  for (l=2;l<=ple->l_max;l++) {
    ll = (double) l;
    for (index_mu=0;index_mu<num_mu;index_mu++) {
      X220[l][index_mu]=0.25*sqrt((ll+2)*(ll+1)*ll*(ll-1)) 
      * exp(-ll*(ll+1)*sigma2[index_mu]/4.);
    }
  }
  return _SUCCESS_;
}

/* Wigner d-functions, computed by recurrence */
/* actual recurrence on sqrt((2l+1)/2) d^l_{mm'} for stability */
/* Formulae from Kostelec & Rockmore 2003 */

int lensing_d00(
		double * mu,
		int num_mu,
		double ** d00
		) {
  double ll, dlm1, dl, dlp1;
  int index_mu, l;
  for (index_mu=0;index_mu<num_mu;index_mu++) {
    dlm1=1.0/sqrt(2.); /* l=0 */
    dl=mu[index_mu] * sqrt(3./2.); /*l=1*/
    for(l=1;l<ple->l_max;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d00 recurrence, supposed to be more stable */ 
      dlp1 = sqrt((2*ll+3)/(2*ll+1))*(2*ll+1)/(ll+1)*mu[index+mu]*dl
      - sqrt((2*ll+3)/(2*ll-1))*ll/(ll+1)*dlm1;
      d00[l+1][index_mu] = dlp1 * sqrt(2./(2*ll+3));
      dlm1 = dl;
      dl = dlp1;
    }
  }
  return _SUCCESS_;
}

int lensing_d11(
                double * mu,
                int num_mu,
                double ** d11
                ) {
  double ll, dlm1, dl, dlp1;
  int index_mu, l;
  for (index_mu=0;index_mu<num_mu;index_mu++) {
    dlm1=(1.0+mu[index_mu])/2. * sqrt(3./2.); /*l=1*/
    dl=(1.0+mu[index_mu])/2.*(2.0*mu[index_mu]-1.0) * sqrt(5./2.); /*l=2*/
    d11[2][index_mu] = dl / sqrt(2./5.);
    for(l=2;l<ple->l_max;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d11 recurrence, supposed to be more stable */
      dlp1 = sqrt((2*ll+3)/(2*ll+1))*(ll+1)*(2*ll+1)/(ll*(ll+2))*(mu[index_mu]-1.0/(ll*(ll+1.)))*dlp1
      - sqrt((2*ll+3)/(2*ll-1))*(ll-1)*(ll+1)/(ll*(ll+2))*(ll+1)/ll*dlm1;
      d11[l+1][index_mu] = dlp1 * sqrt(2./(2*ll+3));
      dlm1 = dl;
      dl = dlp1;
    }
  }
  return _SUCCESS_;
}

int lensing_d1m1(
                double * mu,
                int num_mu,
                double ** d1m1
                ) {
  double ll, dlm1, dl, dlp1;
  int index_mu, l;
  for (index_mu=0;index_mu<num_mu;index_mu++) {
    dlm1=(1.0-mu[index_mu])/2. * sqrt(3./2.); /*l=1*/
    dl=(1.0-mu[index_mu])/2.*(2.0*mu[index_mu]+1.0) * sqrt(5./2.); /*l=2*/
    d1m1[2][index_mu] = dl / sqrt(2./5.);
    for(l=2;l<ple->l_max;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d1m1 recurrence, supposed to be more stable */
      dlp1 = sqrt((2*ll+3)/(2*ll+1))*(ll+1)*(2*ll+1)/(ll*(ll+2))*(mu[index_mu]+1.0/(ll*(ll+1.)))*dl
      - sqrt((2*ll+3)/(2*ll-1))*(ll-1)*(ll+1)/(ll*(ll+2))*(ll+1)/ll*dlm1;
      d1m1[l+1][index_mu] = dlp1 * sqrt(2./(2*ll+3));
      dlm1 = dl;
      dl = dlp1;
    }
  }
  return _SUCCESS_;
}

int lensing_d2m2(
                 double * mu,
                 int num_mu,
                 double ** d2m2
                 ) {
  double ll, dlm1, dl, dlp1;
  int index_mu, l;
  for (index_mu=0;index_mu<num_mu;index_mu++) {
    dlm1=0.; /*l=1*/
    dl=(1.0-mu[index_mu])*(1.0-mu[index_mu])/4. * sqrt(5./2.); /*l=2*/
    d2m2[2][index_mu] = dl / sqrt(2./5.);
    for(l=2;l<ple->l_max;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d2m2 recurrence, supposed to be more stable */
      dlp1 = sqrt((2*ll+3)/(2*ll+1))*(ll+1)*(2*ll+1)/((ll-1)*(ll+3))*(mu[index_mu]+4.0/(ll*(ll+1)))*dl
      - sqrt((2*ll+3)/(2*ll-1))*(ll-2)*(ll+2)/((ll-1)*(ll+3))*(ll+1)/ll*dlm1;
      d2m2[l+1][index_mu] = dlp1 * sqrt(2./(2*ll+3));
      dlm1 = dl;
      dl = dlp1;
    }
  }
  return _SUCCESS_;
}


/** Gauss-Legendre quadrature formulae,
    Strongly inspired from Numerical recipes **/

int lensing_gauss_legendre(
                           double *mu,
                           double *w8,
                           int nmu) {
  
  
}

