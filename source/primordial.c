/** @file primordial.c primordial module
 * Julien Lesgourgues, 6.05.2010    
 *
 * This module computes the primordial spectra
 * (from simple parametric form, or by evolving inflaton perturbations).
 *
 * Hence the following functions can be called from other modules:
 *
 * -# primordial_init() at the beginning (but after precision_init(), input_init())
 * -# primordial_spectrum_at_k() at any time for computing a value P(k) at any k by interpolation
 * -# primordial_free() at the end
 */

#include "primordial.h"

/** 
 * Primordial spectra for arbitrary argument k and for all initial conditions. 
 *
 * Evaluates the primordial spectrum at a given value of k by
 * interpolating in the pre-computed table of (lnk, lnPk).  This function can be
 * called from whatever module at whatever time, provided that
 * primordial_init() has been called before, and primordial_free() has not
 * been called yet.
 *
 * @param ppm        Input: pointer to primordial structure containing tabulated primordial spectrum 
 * @param index_mode Input: index of mode (scalar, tensor, ...) 
 * @param k          Input: wavenumber in 1/Mpc
 * @param pk         Ouput: primordial spectrum P(k) in Mpc**3 for each initial condition, pk[index_ic]
 * @return the error status
 */
int primordial_spectrum_at_k(
			     struct primordial * ppm,
			     int index_mode,
			     double k,
			     double * pk
		) {

  /** Summary: */

  /** - define local variables */

  int index_ic;
  double lnk,lnpk;
  int last_index;

  /** if k negative or null return an error */

  class_test(k<=0.,
	     ppm->error_message,
	     "k = %e",k);

  lnk=log(k);

  /** - if ln(k) is too large to be in the interpolation table, return  an error unless we are in the coase of a analytic spectrum, for which extrapolation is possible */

  if ((lnk > ppm->lnk[ppm->lnk_size-1]) || (lnk < ppm->lnk[0])) {

    class_test(ppm->primordial_spec_type != analytic_Pk,
	       ppm->error_message,
	       "k=%e out of range [%e : %e]",k,exp(ppm->lnk[0]),exp(ppm->lnk[ppm->lnk_size-1]));

    /* extrapolate */
    for (index_ic=0; index_ic < ppm->ic_size[index_mode]; index_ic++) {
      class_call(primordial_analytic_spectrum(ppm,
					      index_mode,
					      index_ic,
					      k,
					      &lnpk),
		 ppm->error_message,
		 ppm->error_message);

      pk[index_ic]=exp(lnpk);
    }
    
    return _SUCCESS_;

  } 

  /** - otherwise, interpolate: */

  class_call(array_interpolate_spline(
				      ppm->lnk,
				      ppm->lnk_size,
				      ppm->lnpk[index_mode],
				      ppm->ddlnpk[index_mode],
				      ppm->ic_size[index_mode],
				      lnk,
				      &last_index,
				      pk,
				      ppm->ic_size[index_mode],
				      ppm->error_message),
	     ppm->error_message,
	     ppm->error_message);

  /** - output P(k), not ln(P(k)) */

  for (index_ic=0; index_ic < ppm->ic_size[index_mode]; index_ic++)
    pk[index_ic]=exp(pk[index_ic]);

  return _SUCCESS_;

}

/**
 * Initialize primordial structure (in particular, compute table of primordial spectrum values)
 * 
 * @param ppr Input : pointer to precision structure (defines method and precision for all computations)
 * @param ppt Input : pointer to perturbation structure (useful for knowing k_min, k_max, etc.)
 * @param ppm Output: pointer to initialized primordial structure 
 * @return the error status
 */
int primordial_init(
		    struct precision  * ppr,
		    struct perturbs   * ppt,
		    struct primordial * ppm
		    ) {

  /** Summary: */

  /** - define local variables */

  double k,k_min,k_max;
  int index_mode,index_ic,index_k;
  double lnpk;

  ppm->lnk_size=0;

  if (ppt->tp_size == NULL) {
    if (ppm->primordial_verbose > 0)
      printf("No perturbations requested. Primordial module skipped.\n");
    return _SUCCESS_;
  }
  else {
    if (ppm->primordial_verbose > 0)
      printf("Computing primordial spectra\n");
  }

  /** - get kmin and kmax from perturbation structure */

  k_min=_HUGE_; /* huge initial value before scanning all modes */
  k_max=0.;    /* zero initial value before scanning all modes */

  for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {

    k_min = min(k_min,ppt->k[index_mode][0]); 
    /* first value, inferred from perturbations structure */

    k_max = max(k_max,ppt->k[index_mode][ppt->k_size[index_mode]-1]); 
    /* last value, inferred from perturbations structure */

  }
  
  class_test(k_min <= 0.,
	     ppm->error_message,
	     "k_min negative or null: stop to avoid segmentation fault");

  class_test(k_max <= 0.,
	     ppm->error_message,
	     "k_max negative or null: stop to avoid segmentation fault");

  class_test(ppm->k_pivot <= 0.,
	     ppm->error_message,
	     "k_pivot negative or null: stop to avoid segmentation fault");

  class_test(ppr->k_per_decade_primordial <= 0.,
	     ppm->error_message,
	     "k_per_decade_primordial negative or null: stop to avoid segmentation fault");
    
  /** - allocate and fill values of lnk's */

  class_call(primordial_get_lnk_list(ppm,
				     k_min,
				     k_max,
				     ppr->k_per_decade_primordial
				     ),
	     ppm->error_message,
	     ppm->error_message);

  class_call(primordial_indices(ppt,
				ppm),
	     ppm->error_message,
	     ppm->error_message);
				
  /** - deal with case of analytic primordial spectra (with titl, running etc.) */

  if (ppm->primordial_spec_type == analytic_Pk) {

    if ((ppm->has_scalars == _TRUE_) && (ppm->has_ad == _TRUE_)) {
      class_test(ppm->A_s_ad <= 0.,
		 ppm->error_message,
		 "stop to avoid segmentation fault");
    }

    if (ppm->has_tensors == _TRUE_) {
      class_test(ppm->r <= 0.,
		 ppm->error_message,
		 "stop to avoid segmentation fault");
    }

    for (index_k = 0; index_k < ppm->lnk_size; index_k++) {

      k=exp(ppm->lnk[index_k]);

      for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {
	
	for (index_ic = 0; index_ic < ppm->ic_size[index_mode]; index_ic++) {

	  class_call(primordial_analytic_spectrum(ppm,
						  index_mode,
						  index_ic,
						  k,
						  &lnpk),
		     ppm->error_message,
		     ppm->error_message);

	  ppm->lnpk[index_mode][index_k*ppm->ic_size[index_mode]+index_ic] = lnpk;

	}

      }

    }

  }

  else {

    class_test(0==0,
	       ppm->error_message,
	       "only analytic primordial spectrum coded yet");

  }     

  /** - compute second derivative of each P(k) with spline, in view of interpolation */

  for (index_mode = 0; index_mode < ppm->md_size; index_mode++) {

    class_call(array_spline_table_lines(ppm->lnk,
					ppm->lnk_size,
					ppm->lnpk[index_mode],
					ppm->ic_size[index_mode],
					ppm->ddlnpk[index_mode],
					_SPLINE_EST_DERIV_,
					ppm->error_message),
	       ppm->error_message,
	       ppm->error_message);

  }

  return _SUCCESS_;

}

/**
 * Free all memory space allocated by primordial_init().
 *
 * To be called at the end of each run.
 *
 * @param ppm Input: pointer to primordial structure (which fields must be freed)
 * @return the error status
 */
int primordial_free(
		    struct primordial * ppm
		    ) {

  int index_mode;

  if (ppm->lnk_size > 0) {

    for (index_mode = 0; index_mode < ppm->md_size; index_mode++) {
      free(ppm->lnpk[index_mode]);
      free(ppm->ddlnpk[index_mode]);
    }
    free(ppm->ic_size);
    free(ppm->lnpk);
    free(ppm->ddlnpk);
    free(ppm->lnk);    
  }

  return _SUCCESS_; 
}

int primordial_indices(
		       struct perturbs   * ppt,
		       struct primordial * ppm
		       ) {

  int index_mode;

  ppm->md_size = ppt->md_size;

  ppm->has_scalars = ppt->has_scalars;

  if (ppm->has_scalars == _TRUE_) {

    ppm->index_md_scalars = ppt->index_md_scalars;

    ppm->has_ad = ppt->has_ad;
    if (ppm->has_ad == _TRUE_) ppm->index_ic_ad = ppt->index_ic_ad;

    ppm->has_bi = ppt->has_bi;
    if (ppm->has_bi == _TRUE_) ppm->index_ic_bi = ppt->index_ic_bi;

    ppm->has_cdi = ppt->has_cdi;
    if (ppm->has_cdi == _TRUE_) ppm->index_ic_cdi = ppt->index_ic_cdi;

    ppm->has_nid = ppt->has_nid;
    if (ppm->has_nid == _TRUE_) ppm->index_ic_nid = ppt->index_ic_nid;

    ppm->has_niv = ppt->has_niv;
    if (ppm->has_niv == _TRUE_) ppm->index_ic_niv = ppt->index_ic_niv;

  }

  ppm->has_tensors = ppt->has_tensors;

  if (ppm->has_tensors == _TRUE_) {

    ppm->index_md_tensors = ppt->index_md_tensors;

    ppm->index_ic_ten = ppt->index_ic_ten;

  }

  class_alloc(ppm->lnpk,ppt->md_size*sizeof(double*),ppm->error_message);

  class_alloc(ppm->ddlnpk,ppt->md_size*sizeof(double*),ppm->error_message);

  class_alloc(ppm->ic_size,ppt->md_size*sizeof(int*),ppm->error_message);

  for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {		     

    class_alloc(ppm->lnpk[index_mode],ppm->lnk_size*ppt->ic_size[index_mode]*sizeof(double),ppm->error_message);

    class_alloc(ppm->ddlnpk[index_mode],ppm->lnk_size*ppt->ic_size[index_mode]*sizeof(double),ppm->error_message);

    ppm->ic_size[index_mode] = ppt->ic_size[index_mode];

  }

  return _SUCCESS_;

}

int primordial_analytic_spectrum(
				 struct primordial * ppm,
				 int index_mode,
				 int index_ic,
				 double k,
				 double * lnpk
				 ) {  

  if ((ppm->has_scalars == _TRUE_) && (index_mode == ppm->index_md_scalars)) {

    if ((ppm->has_ad == _TRUE_) && (index_ic == ppm->index_ic_ad)) {

      /** (a) scalar adiabatic primordial spectrum */
      *lnpk = log(ppm->A_s_ad) 
	+ (ppm->n_s_ad-1.)*log(k/ppm->k_pivot)
	+ 0.5 * ppm->alpha_s_ad * pow(log(k/ppm->k_pivot), 2.);
	      
      return _SUCCESS_;
 
    }

  }

  if ((ppm->has_tensors == _TRUE_) && (index_mode == ppm->index_md_tensors)) {

    /** (b) tensor primordial spectrum */
    *lnpk = log(ppm->A_s_ad*ppm->r/16.) 
      + ppm->n_t * log(k/ppm->k_pivot)
      + 0.5 * ppm->alpha_t * pow(log(k/ppm->k_pivot), 2.);

    return _SUCCESS_;

  }

  class_test(0 == 0,
	     ppm->error_message,
	     "could not recognize which primordial spectrum you want; maybe yet uncoded isocurvature? Or vectors? \n");

}

/**
 * Allocate and fill list of wavenumbers k
 *
 *
 * @param ppm  Input/output: pointer to primordial structure 
 * @param kmin Input : first value
 * @param kmax Input : last value that we should encompass
 * @param k_per_decade Input : number of k per decade
 * @return the error status
 */
int primordial_get_lnk_list(
			    struct primordial * ppm,
			    double kmin,
			    double kmax,
			    double k_per_decade
			    ) {

  int i;

  class_test((kmin <= 0.) || (kmax <= kmin),
	     ppm->error_message,
	     "inconsistent values of kmin=%e, kmax=%e",kmin,kmax);

  ppm->lnk_size = (int)(log(kmax/kmin)/log(10.)*k_per_decade) + 2;

  class_alloc(ppm->lnk,ppm->lnk_size*sizeof(double),ppm->error_message);

  for (i=0; i<ppm->lnk_size; i++)
    ppm->lnk[i]=log(kmin)+i*log(10.)/k_per_decade;
        
  return _SUCCESS_;
  
}
