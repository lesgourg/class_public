/** @file primordial.c Documented Primordial module
 * Julien Lesgourgues, 18.04.2010    
 *
 * This module computes primordial spectra
 * (from simple parametric form, or by evolving inflaton perturbations).
 *
 * Hence the following functions can be called from other modules:
 *
 * -# primordial_init() at the beginning (but after precision_init(), input_init())
 * -# primordial_at_k() at any time for computing a value P(k) at any k by interpolation
 * -# primordial_free() at the end
 */

#include "primordial.h"

/** @name - structures used within the transfer module: */

//@{

struct primordial * ppm; /**< a primordial structure pointer for internal use in the perturbation module */

//@}

/** @name - miscellaneous: */

//@{

ErrorMsg Transmit_Error_Message; /**< contains error message */

//@}

/** 
 * Primordial spectra for arbitrary argument k and for all initial conditions. 
 *
 * Evaluates the primordial spectrum at a given value of k by
 * interpolating in the pre-computed table of (lnk, lnPk).  This function can be
 * called from whatever module at whatever time, provided that
 * primordial_init() has been called before, and primordial_free() has not
 * been called yet.
 *
 * @param index_mode Input: index of mode (scalar, tensor, ...) 
 * @param k          Input: wavenumber in 1/Mpc
 * @param pk         Ouput: primordial spectrum P(k) in Mpc**3 for each initial condition
 * @return the error status
 */
int primordial_at_k(
		    int index_mode,
		    double k,
		    double * pk  /* pk[index_ic] */
		) {

  /** Summary: */

  /** - define local variables */

  int index_ic;
  double lnk;
  int last_index;

  /** if k negative or null return an error */

  if (k<=0.) {
    sprintf(ppm->error_message,"%s(L:%d) : k negative or null",__func__,__LINE__);
    return _FAILURE_;
  } 

  lnk=log(k);

  /** - if ln(k) is too large to be in the interpolation table, return  an error */

  if (lnk > ppm->lnk[ppm->lnk_size-1]) {
    sprintf(ppm->error_message,"%s(L:%d) : k=%e > k_max=%e",__func__,__LINE__,k,exp(ppm->lnk[ppm->lnk_size-1]));
    return _FAILURE_;
  } 

  /** - if ln(k) is too small to be in the interpolation table, return an error */

  if (lnk < ppm->lnk[0]) {
    sprintf(ppm->error_message,"%s(L:%d) : k=%e < k_min=%e",__func__,__LINE__,k,exp(ppm->lnk[0]));
    return _FAILURE_;
  } 

  /** - otherwise, interpolate: */

  if (array_interpolate_spline(
			       ppm->lnk,
			       ppm->lnk_size,
			       ppm->lnpk[index_mode],
			       ppm->ddlnpk[index_mode],
			       ppm->ic_size[index_mode],
			       lnk,
			       &last_index,
			       pk,
			       ppm->ic_size[index_mode],
			       Transmit_Error_Message) == _FAILURE_) {
    sprintf(ppm->error_message,"%s(L:%d) : error in array_interpolate_spline() \n=>%s",__func__,__LINE__,Transmit_Error_Message);
    return _FAILURE_;
  }

  /** - output P(k), not ln(P(k)) */

  for (index_ic=0; index_ic < ppm->ic_size[index_mode]; index_ic++)
    pk[index_ic]=exp(pk[index_ic]);

  return _SUCCESS_;

}

/**
 *
 * @param ppt_input Input : Initialized perturbation structure (useful for knowing k_min, k_max, etc.)
 * @param ppr_input Input : Parameters describing how the computation is to be performed
 * @param ppm_output Output : Initialized primordial structure 
 * @return the error status
 */
int primordial_init(
		struct perturbs * ppt_input,
		struct precision * ppr_input,
		struct primordial * ppm_output
		) {

  /** Summary: */

  /** - define local variables */

  double k,k_min,k_max;
  int index_mode,index_ic,index_k;

  /** - identify the primordial structure ppm (used throughout primordial.c as global variables) to the output structure ppm_output of this function */
  ppm = ppm_output; 
  if (ppm->primordial_verbose > 0)
    printf("Computing primordial spectra\n");

  /** - get kmin and kmax from perturbation structure */

  k_min=1.e10; /* huge initial value before scanning all modes */
  k_max=0.;    /* zero initial value before scanning all modes */

  for (index_mode = 0; index_mode < ppt_input->md_size; index_mode++) {

    k_min = min(k_min,ppt_input->k[index_mode][0]); 
    /* first value, inferred from perturbations structure */

    k_max = max(k_max,ppt_input->k[index_mode][ppt_input->k_size[index_mode]-1]); 
    /* last value, inferred from perturbations structure */

  }
  
  if (k_min <= 0.) {
    sprintf(ppm->error_message,"%s(L:%d) : k_min negative or null: stop to avoid segmentation fault",__func__,__LINE__);
    return _FAILURE_;
  }

  if (k_max <= 0.) {
    sprintf(ppm->error_message,"%s(L:%d) : k_max negative or null: stop to avoid segmentation fault",__func__,__LINE__);
    return _FAILURE_;
  }

   if (ppm->k_pivot <= 0.) {
    sprintf(ppm->error_message,"%s(L:%d) : k_pivot negative or null: stop to avoid segmentation fault",__func__,__LINE__);
    return _FAILURE_;
  }

   if (ppr_input->k_per_decade_primordial <= 0.) {
     sprintf(ppm->error_message,"%s(L:%d) : k_per_decade_primordial negative or null: stop to avoid segmentation fault",__func__,__LINE__);
     return _FAILURE_;
   }

   /** - get number of lnk's */

  if (primordial_get_lnk_list_size(k_min,
				   k_max,
				   ppr_input->k_per_decade_primordial,
				   &(ppm->lnk_size)) == _FAILURE_) {
    sprintf(Transmit_Error_Message,"%s(L:%d) : error in primordial_get_lnk_list_size()\n=>%s",__func__,__LINE__,ppm->error_message);
    sprintf(ppm->error_message,"%s",Transmit_Error_Message);
    return _FAILURE_;
  }
    
  /** - get values of lnk's */

  ppm->lnk =  malloc(ppm->lnk_size*sizeof(double));
  if (ppm->lnk == NULL) {
    sprintf(ppm->error_message,"%s(L:%d) : Could not allocate lnk",__func__,__LINE__);
    return _FAILURE_;
  }

  if (primordial_get_lnk_list(k_min,
			      ppr_input->k_per_decade_primordial,
			      ppm->lnk_size,
			      ppm->lnk) == _FAILURE_) {
    sprintf(Transmit_Error_Message,"%s(L:%d) : error in primordial_get_lnk_list()\n=>%s",__func__,__LINE__,ppm->error_message);
    sprintf(ppm->error_message,"%s",Transmit_Error_Message);
    return _FAILURE_;
  }

  /** - allocate array for lnPk and its second derivative */

  ppm->md_size = ppt_input->md_size;

  ppm->lnpk = malloc(ppt_input->md_size*sizeof(double*));
  if (ppm->lnpk == NULL) {
    sprintf(ppm->error_message,"%s(L:%d) : Could not allocate lnpk",__func__,__LINE__);
    return _FAILURE_;
  }

  ppm->ddlnpk = malloc(ppt_input->md_size*sizeof(double*));
  if (ppm->ddlnpk == NULL) {
    sprintf(ppm->error_message,"%s(L:%d) : Could not allocate ddlnpk",__func__,__LINE__);
    return _FAILURE_;
  }

  ppm->ic_size = malloc(ppt_input->md_size*sizeof(int*));
  if (ppm->ic_size == NULL) {
    sprintf(ppm->error_message,"%s(L:%d) : Could not allocate ic_size",__func__,__LINE__);
    return _FAILURE_;
  }


  for (index_mode = 0; index_mode < ppt_input->md_size; index_mode++) {		     

    ppm->lnpk[index_mode] = malloc(ppm->lnk_size*ppt_input->ic_size[index_mode]*sizeof(double));
    if (ppm->lnpk[index_mode] == NULL) {
      sprintf(ppm->error_message,"%s(L:%d) : Could not allocate lnpk[index_mode]",__func__,__LINE__);
      return _FAILURE_;
    }

    ppm->ddlnpk[index_mode] = malloc(ppm->lnk_size*ppt_input->ic_size[index_mode]*sizeof(double));
    if (ppm->ddlnpk[index_mode] == NULL) {
      sprintf(ppm->error_message,"%s(L:%d) : Could not allocate ddlnpk[index_mode]",__func__,__LINE__);
      return _FAILURE_;
    }

    ppm->ic_size[index_mode] = ppt_input->ic_size[index_mode];

  }
		     
  /** - deal with case of smooth primordial spectra (with titl, running etc.) */

  if (ppm->primordial_spec_type == smooth_Pk) {

    if ((ppt_input->has_ad) && (ppt_input->has_ad)) {
      if (ppm->A_s_ad <= 0.) {
	sprintf(ppm->error_message,"%s(L:%d) : A_s_ad negative or null: stop to avoid segmentation fault",__func__,__LINE__);
	return _FAILURE_;
      }
    }

    for (index_k = 0; index_k < ppm->lnk_size; index_k++) {

      k=exp(ppm->lnk[index_k]);

      for (index_mode = 0; index_mode < ppt_input->md_size; index_mode++) {
	
	for (index_ic = 0; index_ic < ppt_input->ic_size[index_mode]; index_ic++) {

	  if ((ppt_input->has_ad) && (index_ic == ppt_input->index_ic_ad)) {

	    if ((ppt_input->has_scalars == _TRUE_) && (index_mode == ppt_input->index_md_scalars)) {

	      /** (a) scalar adiabatic primordial spectrum */
	      ppm->lnpk[index_mode][index_k*ppm->ic_size[index_mode]+index_ic] = 
		log(ppm->A_s_ad) 
		+ (ppm->n_s_ad-1.) * log(k/ppm->k_pivot)
		+ 0.5 * ppm->alpha_s_ad * pow(log(k/ppm->k_pivot), 2.); 

	    }
	    
	    else {
	      sprintf(ppm->error_message,"%s(L:%d) : isocurvature primordial spectrum not coded yet (although trivial)",__func__,__LINE__);
	      return _FAILURE_;
	    }

	  }

	  else {
	    sprintf(ppm->error_message,"%s(L:%d) : tensor primordial spectrum not coded yet (although trivial)",__func__,__LINE__);
	    return _FAILURE_;
	  }

	}

      }

    }

  }

  else {
    sprintf(ppm->error_message,"%s(L:%d) : only smooth primordial spectrum coded yet",__func__,__LINE__);
    return _FAILURE_;
  }     

  /** - compute second derivative of each P(k) with spline, in view of interpolation */

  for (index_mode = 0; index_mode < ppt_input->md_size; index_mode++) {

    if (array_spline_table_lines(ppm->lnk,
				 ppm->lnk_size,
				 ppm->lnpk[index_mode],
				 ppm->ic_size[index_mode],
				 ppm->ddlnpk[index_mode],
				 _SPLINE_EST_DERIV_,
				 Transmit_Error_Message) == _FAILURE_) {
      sprintf(ppm->error_message,"%s(L:%d) : error in array_spline_table_lines \n=>%s",__func__,__LINE__,Transmit_Error_Message);
      return _FAILURE_;
    }

  }

  return _SUCCESS_;

}

/**
 * Free all memory space allocated by primordial_init().
 *
 * To be called at the end of each run.
 *
 * @return the error status
 */
int primordial_free() {

  int index_mode;

  for (index_mode = 0; index_mode < ppm->md_size; index_mode++) {
    free(ppm->lnpk[index_mode]);
    free(ppm->ddlnpk[index_mode]);
  }
  free(ppm->ic_size);
  free(ppm->lnpk);
  free(ppm->ddlnpk);
  free(ppm->lnk);

  return _SUCCESS_; 
}

/**
 * Define number of wavenumbers k.
 *
 *
 * @param kmin Input : first value
 * @param kmax Input : last value that we should encompass
 * @param k_per_decade Input : number of k per decade
 * @param lnk_list_size Output: number of wavenumbers 
 * @return the error status
 */
int primordial_get_lnk_list_size(
			     double kmin,
			     double kmax,
			     double k_per_decade,
			     int * lnk_list_size
			     ) {

  double number;

  number = log(kmax/kmin)/log(10.)*k_per_decade;

  * lnk_list_size = (int)number + 2;

}

/**
 * Define list of wavenumbers k
 *
 *
 * @param kmin Input : first value
 * @param k_per_decade Input : number of k per decade
 * @param lnk_list_size Input: number of wavenumbers 
 * @param lnpk_list Output: list of wavenumbers 
 * @return the error status
 */
int primordial_get_lnk_list(
			double kmin,
			double k_per_decade,
			int lnk_list_size,
			double * lnk_list
			) {

  int i;

  for (i=0; i<lnk_list_size; i++) {

    lnk_list[i]=log(kmin)+i*log(10.)/k_per_decade;

  }
        
  return _SUCCESS_;
  
}
