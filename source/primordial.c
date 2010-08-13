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

  int index_ic1,index_ic2,index_ic1_ic2;
  double lnk;
  int last_index;

  /** if k negative or null return an error */

  class_test(k<=0.,
	     ppm->error_message,
	     "k = %e",k);

  lnk=log(k);

  /** - if ln(k) is too large to be in the interpolation table, return  an error 
      unless we are in the case of a analytic spectrum, for which an extrapolation is possible */
  
  if ((lnk > ppm->lnk[ppm->lnk_size-1]) || (lnk < ppm->lnk[0])) {

    class_test(ppm->primordial_spec_type != analytic_Pk,
	       ppm->error_message,
	       "k=%e out of range [%e : %e]",k,exp(ppm->lnk[0]),exp(ppm->lnk[ppm->lnk_size-1]));

    /* extrapolate */

    for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_mode]; index_ic1++) {
      for (index_ic2 = index_ic1; index_ic2 < ppm->ic_size[index_mode]; index_ic2++) {
	
/* index value for the coefficients of the symmetric index_ic1*index_ic2 matrix; 
	   takes values between 0 and N(N+1)/2-1 with N=ppm->ic_size[index_mode] */
	index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,ppm->ic_size[index_mode]);

	if (ppm->is_non_zero[index_mode][index_ic1_ic2] == _TRUE_) {

	  class_call(primordial_analytic_spectrum(ppm,
						  index_mode,
						  index_ic1_ic2,
						  k,
						  &(pk[index_ic1_ic2])),
		     ppm->error_message,
		     ppm->error_message);
	}
	else {
	  pk[index_ic1_ic2] = 0.;
	}
      }
    }
    
    return _SUCCESS_;

  } 

  /** - otherwise, interpolate: */

  class_call(array_interpolate_spline(
				      ppm->lnk,
				      ppm->lnk_size,
				      ppm->lnpk[index_mode],
				      ppm->ddlnpk[index_mode],
				      ppm->ic_ic_size[index_mode],
				      lnk,
				      &last_index,
				      pk,
				      ppm->ic_ic_size[index_mode],
				      ppm->error_message),
	     ppm->error_message,
	     ppm->error_message);

  /** - diagonal coefficients in the initial conditions: output P(k), not ln(P(k)) */

  for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_mode]; index_ic1++) {

    index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic1,ppm->ic_size[index_mode]);
    pk[index_ic1_ic2]=exp(pk[index_ic1_ic2]);

  }

  
  /** - non-diagonal coefficients in the initial conditions: output P(k), not cosDelta(k) (Delta = cross-correlation angle) */
  
  for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_mode]; index_ic1++) {
    for (index_ic2 = index_ic1+1; index_ic2 < ppm->ic_size[index_mode]; index_ic2++) {
      
      index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,ppm->ic_size[index_mode]);
      
      if (ppm->is_non_zero[index_mode][index_ic1_ic2] == _TRUE_) {

	pk[index_ic1_ic2] *= sqrt(pk[index_symmetric_matrix(index_ic1,index_ic1,ppm->ic_size[index_mode])]*
				  pk[index_symmetric_matrix(index_ic2,index_ic2,ppm->ic_size[index_mode])]);
      }
      else {
	pk[index_ic1_ic2] = 0.;
      }
    }
  }
  
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
  int index_mode,index_ic1,index_ic2,index_ic1_ic2,index_k;
  double lnpk,pk,pk1,pk2;

  ppm->lnk_size=0;

  if (ppt->has_perturbations == _FALSE_) {
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
		
  /** - deal with case of analytic primordial spectra (with tilt, running etc.) */

  if (ppm->primordial_spec_type == analytic_Pk) {

    class_call(primordial_analytic_spectrum_init(ppt,
						 ppm),
	       ppm->error_message,
	       ppm->error_message);
    
    for (index_k = 0; index_k < ppm->lnk_size; index_k++) {

      k=exp(ppm->lnk[index_k]);

      for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {
	for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_mode]; index_ic1++) {
	  for (index_ic2 = index_ic1; index_ic2 < ppm->ic_size[index_mode]; index_ic2++) {

	    index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,ppm->ic_size[index_mode]);

	    if (ppm->is_non_zero[index_mode][index_ic1_ic2] == _TRUE_) {

	      class_call(primordial_analytic_spectrum(ppm,
						      index_mode,
						      index_ic1_ic2,
						      k,
						      &pk),
			 ppm->error_message,
			 ppm->error_message);
	      
	      if (index_ic1 == index_ic2) {

		/* diagonal coefficients: ln[P(k)] */

		ppm->lnpk[index_mode][index_k*ppm->ic_ic_size[index_mode]+index_ic1_ic2] = log(pk);
	      }
	      else {

		/* non-diagonal coefficients: cosDelta(k) = P(k)_12/sqrt[P(k)_1 P(k)_2] */

		class_call(primordial_analytic_spectrum(ppm,
							index_mode,
							index_symmetric_matrix(index_ic1,index_ic1,ppm->ic_size[index_mode]),
							k,
							&pk1),
			   ppm->error_message,
			   ppm->error_message);			     

		class_call(primordial_analytic_spectrum(ppm,
							index_mode,
							index_symmetric_matrix(index_ic2,index_ic2,ppm->ic_size[index_mode]),
							k,
							&pk2),
			   ppm->error_message,
			   ppm->error_message);	

		ppm->lnpk[index_mode][index_k*ppm->ic_ic_size[index_mode]+index_ic1_ic2] = pk/sqrt(pk1*pk2);
	      }
	    }
	    else {
	      
	      /* non-diagonal coefficients when ic's are uncorrelated */

	      ppm->lnpk[index_mode][index_k*ppm->ic_ic_size[index_mode]+index_ic1_ic2] = 0.;
	    }
	  }
	}
      }
    }
  }

  else {

    class_test(0==0,
	       ppm->error_message,
	       "only analytic primordial spectrum coded yet");

  }     

  /** - compute second derivative of each lnpk versus lnk  with spline, in view of interpolation */

  for (index_mode = 0; index_mode < ppm->md_size; index_mode++) {

    class_call(array_spline_table_lines(ppm->lnk,
					ppm->lnk_size,
					ppm->lnpk[index_mode],
					ppm->ic_ic_size[index_mode],
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

    if (ppm->primordial_spec_type == analytic_Pk) {
      for (index_mode = 0; index_mode < ppm->md_size; index_mode++) {
	free(ppm->is_non_zero[index_mode]);
	free(ppm->amplitude[index_mode]);
	free(ppm->tilt[index_mode]);
	free(ppm->running[index_mode]);
      }
      free(ppm->is_non_zero);
      free(ppm->amplitude);
      free(ppm->tilt);
      free(ppm->running);
    }

    for (index_mode = 0; index_mode < ppm->md_size; index_mode++) {
      free(ppm->lnpk[index_mode]);
      free(ppm->ddlnpk[index_mode]);
    }

    free(ppm->lnpk);
    free(ppm->ddlnpk);
    free(ppm->ic_size);
    free(ppm->ic_ic_size);

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

  class_alloc(ppm->lnpk,ppt->md_size*sizeof(double*),ppm->error_message);

  class_alloc(ppm->ddlnpk,ppt->md_size*sizeof(double*),ppm->error_message);

  class_alloc(ppm->ic_size,ppt->md_size*sizeof(int*),ppm->error_message);

  class_alloc(ppm->ic_ic_size,ppt->md_size*sizeof(int*),ppm->error_message);

  for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {		     

    ppm->ic_size[index_mode] = ppt->ic_size[index_mode];

    ppm->ic_ic_size[index_mode] = (ppm->ic_size[index_mode]*(ppm->ic_size[index_mode]+1))/2;

    class_alloc(ppm->lnpk[index_mode],ppm->lnk_size*ppm->ic_ic_size[index_mode]*sizeof(double),ppm->error_message);

    class_alloc(ppm->ddlnpk[index_mode],ppm->lnk_size*ppm->ic_ic_size[index_mode]*sizeof(double),ppm->error_message);

  }

  return _SUCCESS_;

}

int primordial_analytic_spectrum_init(
				      struct perturbs   * ppt,
				      struct primordial * ppm
				      ) {

  int index_mode,index_ic1,index_ic2;
  int index_ic1_ic2,index_ic1_ic1,index_ic2_ic2;
  double one_amplitude,one_tilt,one_running,one_correlation;

  class_alloc(ppm->is_non_zero,
	      ppm->md_size*sizeof(short *),
	      ppm->error_message);

  class_alloc(ppm->amplitude,
	      ppm->md_size*sizeof(double *),
	      ppm->error_message);

  class_alloc(ppm->tilt,
	      ppm->md_size*sizeof(double *),
	      ppm->error_message);

  class_alloc(ppm->running,
	      ppm->md_size*sizeof(double *),
	      ppm->error_message);

  for (index_mode = 0; index_mode < ppm->md_size; index_mode++) {

    class_alloc(ppm->is_non_zero[index_mode],
		ppm->ic_ic_size[index_mode]*sizeof(short),
		ppm->error_message);

    class_alloc(ppm->amplitude[index_mode],
		ppm->ic_ic_size[index_mode]*sizeof(double),
		ppm->error_message);

    class_alloc(ppm->tilt[index_mode],
		ppm->ic_ic_size[index_mode]*sizeof(double),
		ppm->error_message);

    class_alloc(ppm->running[index_mode],
		ppm->ic_ic_size[index_mode]*sizeof(double),
		ppm->error_message);

  }
      
  for (index_mode = 0; index_mode < ppm->md_size; index_mode++) {

    /* digonal coefficients */
    
    for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_mode]; index_ic1++) {

      if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {

	if ((ppt->has_ad == _TRUE_) && (index_ic1 == ppt->index_ic_ad)) {
	  one_amplitude = ppm->A_s;
	  one_tilt = ppm->n_s;
	  one_running = ppm->alpha_s;
	}

	if ((ppt->has_bi == _TRUE_) && (index_ic1 == ppt->index_ic_bi)) {
	  one_amplitude = ppm->A_s*ppm->f_bi*ppm->f_bi;
	  one_tilt = ppm->n_bi;
	  one_running = ppm->alpha_bi;
	}

	if ((ppt->has_cdi == _TRUE_) && (index_ic1 == ppt->index_ic_cdi)) {
	  one_amplitude = ppm->A_s*ppm->f_cdi*ppm->f_cdi;
	  one_tilt = ppm->n_cdi;
	  one_running = ppm->alpha_cdi;
	}

	if ((ppt->has_nid == _TRUE_) && (index_ic1 == ppt->index_ic_nid)) {
	  one_amplitude = ppm->A_s*ppm->f_nid*ppm->f_nid;
	  one_tilt = ppm->n_nid;
	  one_running = ppm->alpha_nid;
	}
    
	if ((ppt->has_niv == _TRUE_) && (index_ic1 == ppt->index_ic_niv)) {
	  one_amplitude = ppm->A_s*ppm->f_niv*ppm->f_niv;
	  one_tilt = ppm->n_niv;
	  one_running = ppm->alpha_niv;
	}
      }

      if ((ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors)) {

	if (index_ic1 == ppt->index_ic_ten) {
	  one_amplitude = ppm->A_s*ppm->r;
	  one_tilt = ppm->n_t+1.; /* +1 to match usual definition of n_t (equivalent to n_s-1) */
	  one_running = ppm->alpha_t;
	}
      }

      class_test(one_amplitude <= 0.,
		 ppm->error_message,
		 "inconsistent input for primordial amplitude: %g for index_mode=%d, index_ic=%d\n",
		 one_amplitude,index_mode,index_ic1);

      index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic1,ppm->ic_size[index_mode]);
      
      ppm->is_non_zero[index_mode][index_ic1_ic2] = _TRUE_;
      ppm->amplitude[index_mode][index_ic1_ic2] = one_amplitude;
      ppm->tilt[index_mode][index_ic1_ic2] = one_tilt;
      ppm->running[index_mode][index_ic1_ic2] = one_running;
    }

    /* non-diagonal coefficients */

    for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_mode]; index_ic1++) {
      for (index_ic2 = index_ic1+1; index_ic2 < ppm->ic_size[index_mode]; index_ic2++) {
     
	if ((ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) {
 
	  if ((ppt->has_ad == _TRUE_) && (ppt->has_bi == _TRUE_) && 
	      (((index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_bi)) ||
	       ((index_ic1 == ppt->index_ic_ad) && (index_ic1 == ppt->index_ic_bi)))) {
	    one_correlation = ppm->c_ad_bi;
	    one_tilt = ppm->n_ad_bi;
	    one_running = ppm->alpha_ad_bi;
	  }

	  if ((ppt->has_ad == _TRUE_) && (ppt->has_cdi == _TRUE_) && 
	      (((index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_cdi)) ||
	       ((index_ic2 == ppt->index_ic_ad) && (index_ic1 == ppt->index_ic_cdi)))) {
	    one_correlation = ppm->c_ad_cdi;
	    one_tilt = ppm->n_ad_cdi;
	    one_running = ppm->alpha_ad_cdi;
	  }

	  if ((ppt->has_ad == _TRUE_) && (ppt->has_nid == _TRUE_) && 
	      (((index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_nid)) ||
	       ((index_ic2 == ppt->index_ic_ad) && (index_ic1 == ppt->index_ic_nid)))) {
	    one_correlation = ppm->c_ad_nid;
	    one_tilt = ppm->n_ad_nid;
	    one_running = ppm->alpha_ad_nid;
	  }

	  if ((ppt->has_ad == _TRUE_) && (ppt->has_niv == _TRUE_) && 
	      (((index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_niv)) ||
	       ((index_ic2 == ppt->index_ic_ad) && (index_ic1 == ppt->index_ic_niv)))) {
	    one_correlation = ppm->c_ad_niv;
	    one_tilt = ppm->n_ad_niv;
	    one_running = ppm->alpha_ad_niv;
	  }

	  if ((ppt->has_bi == _TRUE_) && (ppt->has_cdi == _TRUE_) && 
	      (((index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_cdi)) ||
	       ((index_ic2 == ppt->index_ic_bi) && (index_ic1 == ppt->index_ic_cdi)))) {
	    one_correlation = ppm->c_bi_cdi;
	    one_tilt = ppm->n_bi_cdi;
	    one_running = ppm->alpha_bi_cdi;
	  }

	  if ((ppt->has_bi == _TRUE_) && (ppt->has_nid == _TRUE_) && 
	      (((index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_nid)) ||
	       ((index_ic2 == ppt->index_ic_bi) && (index_ic1 == ppt->index_ic_nid)))) {
	    one_correlation = ppm->c_bi_nid;
	    one_tilt = ppm->n_bi_nid;
	    one_running = ppm->alpha_bi_nid;
	  }

	  if ((ppt->has_bi == _TRUE_) && (ppt->has_niv == _TRUE_) && 
	      (((index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_niv)) ||
	       ((index_ic2 == ppt->index_ic_bi) && (index_ic1 == ppt->index_ic_niv)))) {
	    one_correlation = ppm->c_bi_niv;
	    one_tilt = ppm->n_bi_niv;
	    one_running = ppm->alpha_bi_niv;
	  }

	  if ((ppt->has_cdi == _TRUE_) && (ppt->has_nid == _TRUE_) && 
	      (((index_ic1 == ppt->index_ic_cdi) && (index_ic2 == ppt->index_ic_nid)) ||
	       ((index_ic2 == ppt->index_ic_cdi) && (index_ic1 == ppt->index_ic_nid)))) {
	    one_correlation = ppm->c_cdi_nid;
	    one_tilt = ppm->n_cdi_nid;
	    one_running = ppm->alpha_cdi_nid;
	  }

	  if ((ppt->has_cdi == _TRUE_) && (ppt->has_niv == _TRUE_) && 
	      (((index_ic1 == ppt->index_ic_cdi) && (index_ic2 == ppt->index_ic_niv)) ||
	       ((index_ic2 == ppt->index_ic_cdi) && (index_ic1 == ppt->index_ic_niv)))) {
	    one_correlation = ppm->c_cdi_niv;
	    one_tilt = ppm->n_cdi_niv;
	    one_running = ppm->alpha_cdi_niv;
	  }

	  if ((ppt->has_nid == _TRUE_) && (ppt->has_niv == _TRUE_) && 
	      (((index_ic1 == ppt->index_ic_nid) && (index_ic2 == ppt->index_ic_niv)) ||
	       ((index_ic2 == ppt->index_ic_nid) && (index_ic1 == ppt->index_ic_niv)))) {
	    one_correlation = ppm->c_nid_niv;
	    one_tilt = ppm->n_nid_niv;
	    one_running = ppm->alpha_nid_niv;
	  }

	}

	class_test((one_correlation < -1) || (one_correlation > 1),
		   ppm->error_message,
		   "inconsistent input for isocurvature cross-correlation\n");

	index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,ppm->ic_size[index_mode]);
	index_ic1_ic1 = index_symmetric_matrix(index_ic1,index_ic1,ppm->ic_size[index_mode]);
	index_ic2_ic2 = index_symmetric_matrix(index_ic2,index_ic2,ppm->ic_size[index_mode]);

	if (one_correlation == 0.) {
	  ppm->is_non_zero[index_mode][index_ic1_ic2] = _FALSE_;
	  ppm->amplitude[index_mode][index_ic1_ic2] = 0.;
	  ppm->tilt[index_mode][index_ic1_ic2] = 0.;
	  ppm->running[index_mode][index_ic1_ic2] = 0.;
	}
	else {
	  ppm->is_non_zero[index_mode][index_ic1_ic2] = _TRUE_;
	  ppm->amplitude[index_mode][index_ic1_ic2] = 
	    sqrt(ppm->amplitude[index_mode][index_ic1_ic1]*
		 ppm->amplitude[index_mode][index_ic2_ic2])*
	    one_correlation;
	  ppm->tilt[index_mode][index_ic1_ic2] = 
	    0.5*(ppm->tilt[index_mode][index_ic1_ic1]
		 +ppm->tilt[index_mode][index_ic2_ic2])
	    + one_tilt;
	  ppm->running[index_mode][index_ic1_ic2] = 
	    0.5*(ppm->running[index_mode][index_ic1_ic1]
		 +ppm->running[index_mode][index_ic2_ic2])
	    + one_running;
	}
      }
    }
  }
  
  return _SUCCESS_;

}

int primordial_analytic_spectrum(
				 struct primordial * ppm,
				 int index_mode,
				 int index_ic1_ic2,
				 double k,
				 double * pk
				 ) {  
  
  if (ppm->is_non_zero[index_mode][index_ic1_ic2] == _TRUE_) {

    *pk = ppm->amplitude[index_mode][index_ic1_ic2]
      *exp((ppm->tilt[index_mode][index_ic1_ic2]-1.)*log(k/ppm->k_pivot)
	   + 0.5 * ppm->running[index_mode][index_ic1_ic2] * pow(log(k/ppm->k_pivot), 2.));
    
  }

  else {

    class_test(0==0,
	       ppm->error_message,
	       "you are calling this routine for an uncorrelated pair of ic's, this should never happen\n");
    
  }

  return _SUCCESS_;
  
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
