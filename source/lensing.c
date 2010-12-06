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

  double * beta; /* beta[index_beta]: discretized values of beta
		    between -pi and pi */

  double ** d00;  /* dmn[index_l][index_beta] */
  double ** d11;
  double ** dm11;
  double ** d2m2; 
  double ** d22;  
  double ** d20;  
  double ** d1m1; 
  double ** d31;  
  double ** d40;  
  double ** d3m3; 
  double ** d4m4; 

  double * Cgl;   /* Cgl[index_beta] */
  double * Cgl2;  /* Cgl2[index_beta] */
  double * sigma2; /* sigma[index_beta] */

  double ** X000;  /* Ximn[index_l][index_beta] */ 
  double ** Xp000;
  double ** X220; 
  double ** X022; 
  double ** X132; 
  double ** X121; 
  double ** Xp022; 
  double ** X242; 
  double ** X112; 

  int num_beta,index_beta;
  int l;
  double * cl_unlensed;  /* cl_unlensed[index_ct] */
  double ** junk1, ** junk2;

  /** put here all precision variables; will be stored later in precision structure */

  int beta_sampling=20;
   

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

  /** - allocate array of beta values */

  num_beta = 2*beta_sampling+1;

  class_alloc(beta,
	      num_beta*sizeof(double),
	      ple->error_message);

  for (index_beta=0; index_beta<num_beta; index_beta++) {

    beta[index_beta] = _PI_*(index_beta-beta_sampling)/beta_sampling;

  }

  /** - compute d^l_mm'(beta) */

  class_alloc(d00,
	      (ple->l_max+1)*sizeof(double*),
	      ple->error_message);

  for (l=2; l<=ple->l_max; l++) {

    class_alloc(d00[l],
		num_beta*sizeof(double),
		ple->error_message);  
  }

  class_call(lensing_d00(beta,num_beta,d00),
	     ple->error_message,
	     ple->error_message);

  /** - compute Cgl(beta), Cgl2(beta) and sigma(beta) */

  class_alloc(Cgl,
	      num_beta*sizeof(double),
	      ple->error_message);

  class_alloc(Cgl2,
	      num_beta*sizeof(double),
	      ple->error_message);

  class_alloc(sigma2,
	      num_beta*sizeof(double),
	      ple->error_message);

  class_alloc(cl_unlensed,
	      psp->ct_size*sizeof(double),
	      ple->error_message);

  for (index_beta=0; index_beta<num_beta; index_beta++) {

    Cgl[index_beta]=0;
    Cgl2[index_beta]=0;

    for (l=2; l<=ple->l_max; l++) {

      class_call(spectra_cl_at_l(psp,l,cl_unlensed,junk1,junk2),
		 psp->error_message,
		 ple->error_message);

      Cgl[index_beta] += (2.*l+1.)*l*(l+1.)*
	cl_unlensed[psp->index_ct_pp]*d11[l][index_beta];

      Cgl2[index_beta] += (2.*l+1.)*l*(l+1.)*
	cl_unlensed[psp->index_ct_pp]*dm11[l][index_beta];

    }

    Cgl[index_beta] /= 4.*_PI_;
    Cgl[index_beta] /= 4.*_PI_;

  }

  for (index_beta=0; index_beta<num_beta; index_beta++) {

    sigma2[index_beta] = Cgl[beta_sampling] - Cgl[index_beta];

  }

  /** - compute X000(beta), X'000(beta), X220 and other Ximn */

  /** - compute ksi, ksi+, ksi-, ksiX */

  /** - compute lensed Cls */

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

int lensing_d00(
		double * beta,
		int num_beta,
		double ** d00
		) {
  
  return _SUCCESS_;

}
