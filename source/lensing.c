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
		    between -1 and 1, roots of Legendre polynomial */
  double * w8; /* Corresponding Gauss-Legendre quadrature weights */
  
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
  
  double * ksi;  /* ksi[index_mu] */

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

  /** put here all precision variables; will be stored later in precision structure */
  /** Last element in mu will be for mu=1, needed for sigma2 
   The rest will be chosen as roots of a Gauss-Legendre quadrature **/
  
  num_mu=(ple->l_max+1000); /* Must be even ?? CHECK */
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

  class_call(lensing_d00(mu,num_mu,ple->l_max,d00),
	     ple->error_message,
	     ple->error_message);

  class_call(lensing_d11(mu,num_mu,ple->l_max,d11),
	     ple->error_message,
	     ple->error_message);

  class_call(lensing_d1m1(mu,num_mu,ple->l_max,d1m1),
	     ple->error_message,
	     ple->error_message);

  class_call(lensing_d2m2(mu,num_mu,ple->l_max,d2m2),
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
	cl_unlensed[psp->index_ct_pp]*d1m1[l][index_mu];

    }

    Cgl[index_mu] /= 4.*_PI_;
    Cgl2[index_mu] /= 4.*_PI_;

  }

  for (index_mu=0; index_mu<num_mu; index_mu++) {
    /* Cgl(1.0) - Cgl(mu) */
    sigma2[index_mu] = Cgl[num_mu-1] - Cgl[index_mu];
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
  
  class_call(lensing_X000(mu,num_mu-1,ple->l_max,sigma2,X000),
             ple->error_message,
             ple->error_message);
  class_call(lensing_Xp000(mu,num_mu-1,ple->l_max,sigma2,Xp000),
             ple->error_message,
             ple->error_message);
  class_call(lensing_X220(mu,num_mu-1,ple->l_max,sigma2,X220),
             ple->error_message,
             ple->error_message);
  
  

  /** - compute ksi, ksi+, ksi-, ksiX */
  class_alloc(ksi,
              (num_mu-1)*sizeof(double),
              ple->error_message);
  {
    double res;
    for (index_mu=0;index_mu<num_mu-1;index_mu++) {
      ksi[index_mu]=0;
      for (l=2;l<=ple->l_max;l++) {
        ll = (double) l;
        class_call(spectra_cl_at_l(psp,l,cl_unlensed,junk1,junk2),
                   psp->error_message,
                   ple->error_message);
        res = (2*ll+1)/(4.*_PI_)*cl_unlensed[psp->index_ct_tt];
        /* res *= d00[l][index_mu]; */ /*DEBUG: TO BE REMOVED */
        
        res *= (X000[l][index_mu]*X000[l][index_mu]*d00[l][index_mu] +
                Xp000[l][index_mu]*Xp000[l][index_mu]*d1m1[l][index_mu]
                *Cgl2[index_mu]*8./(ll*(ll+1)) +
                (Xp000[l][index_mu]*Xp000[l][index_mu]*d00[l][index_mu] +
                 X220[l][index_mu]*X220[l][index_mu]*d2m2[l][index_mu])
                *Cgl2[index_mu]*Cgl2[index_mu]);
        
        ksi[index_mu] += res;
      }
    }
  }
  /** - compute lensed Cls by integration */
  class_call(lensing_lensed_cl(ksi,d00,w8,num_mu-1,ple),
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


int lensing_lensed_cl(
        double *ksi, 
        double **d00,
        double *w8,
        int nmu,
        struct lensing * ple
        ) {
  
  double cle;
  int l, imu;
  /** Integration by Gauss-Legendre quadrature **/
  for(l=2;l<=ple->l_max;l++){
    cle=0;
    for (imu=0;imu<nmu;imu++) {
      cle += ksi[imu]*d00[l][imu]*w8[imu];
    }
    ple->cl_lensed[l*ple->lt_size+ple->index_lt_tt]=cle*2.0*_PI_;
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
  for (l=2;l<=lmax;l++) {
    ll = (double) l;
    for (index_mu=0;index_mu<num_mu;index_mu++) {
      X000[l][index_mu]=exp(-ll*(ll+1)*sigma2[index_mu]/4.);
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
  for (l=2;l<=lmax;l++) {
    ll = (double) l;
    for (index_mu=0;index_mu<num_mu;index_mu++) {
      Xp000[l][index_mu]=-ll*(ll+1)/4.*exp(-ll*(ll+1)*sigma2[index_mu]/4.);
    }
  }
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
  for (l=2;l<=lmax;l++) {
    ll = (double) l;
    for (index_mu=0;index_mu<num_mu;index_mu++) {
      X220[l][index_mu]=0.25*sqrt((ll+2)*(ll+1)*ll*(ll-1)) 
      * exp(-ll*(ll+1)*sigma2[index_mu]/4.);
    }
  }
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
  for (index_mu=0;index_mu<num_mu;index_mu++) {
    dlm1=1.0/sqrt(2.); /* l=0 */
    dl=mu[index_mu] * sqrt(3./2.); /*l=1*/
    for(l=1;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d00 recurrence, supposed to be more stable */ 
      dlp1 = sqrt((2*ll+3)/(2*ll+1))*(2*ll+1)/(ll+1)*mu[index_mu]*dl
      - sqrt((2*ll+3)/(2*ll-1))*ll/(ll+1)*dlm1;
      d00[l+1][index_mu] = dlp1 * sqrt(2./(2*ll+3));
      dlm1 = dl;
      dl = dlp1;
    }
  }
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
  for (index_mu=0;index_mu<num_mu;index_mu++) {
    dlm1=(1.0+mu[index_mu])/2. * sqrt(3./2.); /*l=1*/
    dl=(1.0+mu[index_mu])/2.*(2.0*mu[index_mu]-1.0) * sqrt(5./2.); /*l=2*/
    d11[2][index_mu] = dl / sqrt(2./5.);
    for(l=2;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d11 recurrence, supposed to be more stable */
      dlp1 = sqrt((2*ll+3)/(2*ll+1))*(ll+1)*(2*ll+1)/(ll*(ll+2))*(mu[index_mu]-1.0/(ll*(ll+1.)))*dl
      - sqrt((2*ll+3)/(2*ll-1))*(ll-1)*(ll+1)/(ll*(ll+2))*(ll+1)/ll*dlm1;
      d11[l+1][index_mu] = dlp1 * sqrt(2./(2*ll+3));
      dlm1 = dl;
      dl = dlp1;
    }
  }
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
  for (index_mu=0;index_mu<num_mu;index_mu++) {
    dlm1=(1.0-mu[index_mu])/2. * sqrt(3./2.); /*l=1*/
    dl=(1.0-mu[index_mu])/2.*(2.0*mu[index_mu]+1.0) * sqrt(5./2.); /*l=2*/
    d1m1[2][index_mu] = dl / sqrt(2./5.);
    for(l=2;l<lmax;l++){
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
  for (index_mu=0;index_mu<num_mu;index_mu++) {
    dlm1=0.; /*l=1*/
    dl=(1.0-mu[index_mu])*(1.0-mu[index_mu])/4. * sqrt(5./2.); /*l=2*/
    d2m2[2][index_mu] = dl / sqrt(2./5.);
    for(l=2;l<lmax;l++){
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

