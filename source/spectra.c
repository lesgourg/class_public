/** @file cl.c Documented \f$ C_l^{X}, P(k), ... \f$ module
 * Julien Lesgourgues, 18.04.2010    
 *
 * This module computes the power spectra \f$ C_l^{X}, P(k), ... \f$'s given the transfer functions and
 * primordial spectra.
 */

#include "spectra.h"

/** @name - structures used within the transfer module: */

//@{

struct perturbs * ppt; /**< a perturbs structure pointer for internal use in the perturbation module */
struct transfers * ptr; /**< a transfers structure pointer for internal use in the perturbation module */
struct spectra * psp; /**< a spectra structure pointer for internal use in the perturbation module */

//@}
/** @name - miscellaneous: */

//@{

char * errmsg; /**< error management pointer */
char Transmit_Error_Message[2048]; /**< contains error message */

//@}

int spectra_cl_at_l(
		    double l,
		    int index_mode,
		    double * cl
		    ) {

  int * last_index;

  if ((index_mode < 0) && (index_mode >= ppt->md_size)) {
    sprintf(psp->error_message,"%s(L:%d) : index_mode=%d out of bounds",__func__,__LINE__,index_mode);
    return _FAILURE_;
  }

  if ((l < ptr->l[index_mode][0]) && (l >= ptr->l[index_mode][ptr->l_size[index_mode]])) {
    sprintf(psp->error_message,"%s(L:%d) : l=%d out of range [%d:%d]",
	    __func__,__LINE__,l,ptr->l[index_mode][0],ptr->l[index_mode][ptr->l_size[index_mode]]);
    return _FAILURE_;
  }

  if (array_interpolate_spline(psp->l[index_mode],
			       psp->l_size[index_mode],
			       psp->cl[index_mode],
			       psp->ddcl[index_mode],
			       ppt->ic_size[index_mode]*psp->ct_size,
			       l,
			       last_index,
			       cl,
			       ppt->ic_size[index_mode]*psp->ct_size,
			       errmsg) == _FAILURE_) {
    sprintf(psp->error_message,"%s(L:%d) : error in array_interpolate_spline()\n=>%s",__func__,__LINE__,errmsg);
    return _FAILURE_;
  }

  return _SUCCESS_;

}

/**
 * Computes the \f$ C_l^{X}, P(k), ... \f$'s
 *
 * @param ppt_intput Intput : Initialized perturbation structure
 * @param ptr_input Input : Initialized transfers structure
 * @param ppm_input Input : Initialized primordial structure
 * @param pcl_output Output : Initialized cls structure
 * @return the error status
 */
int spectra_init(
		 struct perturbs * ppt_input,
		 struct transfers * ptr_input,
		 struct primordial * ppm_input,
		 struct spectra * psp_output
		 ) {

  /** - define local variables */
  int index_mode; /* index running over modes (scalar, tensor, ...) */
  int index_ic; /* index running over initial conditions */
  int index_type; /* index running over types (temperature, polarisation, ...) */
  int index_k; /* index running over wavenumber */
  int index_l;  /* multipoles */
  int index_ct;
  double k; /* wavenumber */
  double clvalue;
  int cl_integrand_num_columns;
  double * cl_integrand;
  double * trsf;
  double * primordial_pk;  /*pk[index_ic]*/

  /** - identify the spectra structure (used throughout transfer.c as global variable) to the input/output structure of this function */
  ppt = ppt_input; 
  ptr = ptr_input; 
  psp = psp_output; 

  if (psp->spectra_verbose > 0)
    printf("Computing output spectra\n");

  if (spectra_indices() == _FAILURE_) {
    sprintf(Transmit_Error_Message,"%s(L:%d) : error in spectra_indices()\n=>%s",__func__,__LINE__,psp->error_message);
    sprintf(psp->error_message,Transmit_Error_Message);
    return _FAILURE_;
  }

  /** - deal with cl's, if any */

  if (ptr_input->has_cls == _TRUE_) {

    psp->l_size = malloc (sizeof(int)*ppt->md_size);
    if (psp->l_size == NULL) {
      sprintf(psp->error_message,"%s(L:%d) : Could not allocate l_size",__func__,__LINE__);
      return _FAILURE_;
    }

    psp->l = malloc (sizeof(double *)*ppt->md_size);
    if (psp->l == NULL) {
      sprintf(psp->error_message,"%s(L:%d) : Could not allocate l",__func__,__LINE__);
      return _FAILURE_;
    }

    psp->cl = malloc (sizeof(double *)*ppt_input->md_size);
    if (psp->cl == NULL) {
      sprintf(psp->error_message,"%s(L:%d) : Could not allocate cl",__func__,__LINE__);
      return _FAILURE_;
    }

    psp->ddcl = malloc (sizeof(double *)*ppt_input->md_size);
    if (psp->ddcl == NULL) {
      sprintf(psp->error_message,"%s(L:%d) : Could not allocate ddcl",__func__,__LINE__);
      return _FAILURE_;
    }

    /** - loop over modes (scalar, tensors, etc) */
    for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {

      psp->l_size[index_mode] = ptr->l_size[index_mode];

      psp->l[index_mode] = malloc(sizeof(double)*psp->l_size[index_mode]);
      if (psp->l[index_mode] == NULL) {
	sprintf(psp->error_message,"%s(L:%d) : Could not allocate l[index_mode]",__func__,__LINE__);
	return _FAILURE_;
      }


      for (index_l=0; index_l < ptr->l_size[index_mode]; index_l++) {
	psp->l[index_mode][index_l] = (double)ptr->l[index_mode][index_l];
      }

      psp->cl[index_mode] = malloc(sizeof(double)*psp->ct_size*ppt_input->ic_size[index_mode]*ptr->l_size[index_mode]);
      if (psp->cl[index_mode] == NULL) {
	sprintf(psp->error_message,"%s(L:%d) : Could not allocate cl[index_mode]",__func__,__LINE__);
	return _FAILURE_;
      }

      psp->ddcl[index_mode] = malloc(sizeof(double)*psp->ct_size*ppt_input->ic_size[index_mode]*ptr->l_size[index_mode]);
      if (psp->ddcl[index_mode] == NULL) {
	sprintf(psp->error_message,"%s(L:%d) : Could not allocate ddcl[index_mode]",__func__,__LINE__);
	return _FAILURE_;
      }

      cl_integrand_num_columns = 1+psp->ct_size*2; /* one for k, ct_size for each type, ct_size for each second derivative of each type */

      cl_integrand = malloc(ptr->k_size[index_mode]*cl_integrand_num_columns*sizeof(double));
      if (cl_integrand == NULL) {
	sprintf(psp->error_message,"%s(L:%d) : Could not allocate cl_integrand",__func__,__LINE__);
	return _FAILURE_;
      }

      primordial_pk=malloc(ppt_input->ic_size[index_mode]);
      if (primordial_pk == NULL) {
	sprintf(psp->error_message,"%s(L:%d) : Could not allocate primordial_pk",__func__,__LINE__);
	return _FAILURE_;
      }

      trsf=malloc(ppt_input->tp_size);
      if (trsf == NULL) {
	sprintf(psp->error_message,"%s(L:%d) : Could not allocate trsf",__func__,__LINE__);
	return _FAILURE_;
      }

      /** - loop over initial conditions */
      for (index_ic = 0; index_ic < ppt->ic_size[index_mode]; index_ic++) {

	/** - loop over l values defined in the transfer module. For each l: */
	for (index_l=0; index_l < ptr->l_size[index_mode]; index_l++) {

	  /** - compute integrand \f$ \Delta_l(k)^2 / k \f$ for each k in the transfer function's table (assumes flat primordial spectrum) */

	  for (index_k=0; index_k < ptr->k_size[index_mode]; index_k++) {

	    k = ptr->k[index_mode][index_k];

	    cl_integrand[index_k*cl_integrand_num_columns+0] = k;
	    
	    if (primordial_at_k(index_mode,k,primordial_pk) == _FAILURE_) {
	      sprintf(psp->error_message,"%s(L:%d) : error in primordial_at_k()\n=>%s",__func__,__LINE__,ppm_input->error_message);
	      return _FAILURE_;
	    }

	    /* above routine checks that k>0: no possible division by zero below */

	    for (index_type=0; index_type < ppt->tp_size; index_type++) {

	      trsf[index_type] = 
		ptr->transfer[index_mode]
		[((index_ic * ppt->tp_size + index_type)
		  * ptr->l_size[index_mode] + index_l)
		 * ptr->k_size[index_mode] + index_k];

	    }

	    if (ppt->has_source_t == _TRUE_) {
	      cl_integrand[index_k*cl_integrand_num_columns+1+psp->index_ct_tt]=primordial_pk[index_ic]
		* trsf[ppt->index_tp_t]
		* trsf[ppt->index_tp_t]
		/ k;
	    }
	      
	    if (ppt->has_source_p == _TRUE_) {
	      cl_integrand[index_k*cl_integrand_num_columns+1+psp->index_ct_ee]=primordial_pk[index_ic]
		* trsf[ppt->index_tp_p]
		* trsf[ppt->index_tp_p]
		/ k;
	    }

	    if ((ppt->has_source_t == _TRUE_) && (ppt->has_source_p == _TRUE_)) {
	      cl_integrand[index_k*cl_integrand_num_columns+1+psp->index_ct_te]=primordial_pk[index_ic]
		* trsf[ppt->index_tp_t]
		* trsf[ppt->index_tp_p]
		/ k;
	    }

	    if ((ppt->has_source_p == _TRUE_) && (ppt->has_tensors == _TRUE_)) {

	      if (index_mode == ppt->index_md_scalars) {
		cl_integrand[index_k*cl_integrand_num_columns+1+psp->index_ct_bb]=0.;
	      }
	      if (index_mode == ppt->index_md_tensors) {
		sprintf(psp->error_message,"%s(L:%d) : tensors not coded yet",__func__,__LINE__);
		return _FAILURE_;
	      }
	      

	    }
	  }

	  for (index_ct=0; index_ct<psp->ct_size; index_ct++) {

	    if (array_spline(cl_integrand,
			     cl_integrand_num_columns,
			     ptr_input->k_size[index_mode],
			     0,
			     1+index_ct,
			     1+psp->ct_size+index_ct,
			     _SPLINE_EST_DERIV_,
			     errmsg) == _FAILURE_) {
	      sprintf(psp->error_message,"%s(L:%d) : error in array_spline_table_lines \n=>%s",__func__,__LINE__,errmsg);
	      return _FAILURE_;
	    }
	    
	    if (array_integrate_all_spline(cl_integrand,
					   cl_integrand_num_columns,
					   ptr_input->k_size[index_mode],
					   0,
					   1+index_ct,
					   1+psp->ct_size+index_ct,
					   &clvalue,
					   errmsg) == _FAILURE_) {
	      sprintf(psp->error_message,"%s(L:%d) : error in array_spline_table_lines \n=>%s",__func__,__LINE__,errmsg);
	      return _FAILURE_;
	    }

	    psp->cl[index_mode]
	      [(index_l * ppt->ic_size[index_mode] + index_ic) * psp->ct_size + index_ct]
	      = clvalue;

	  }

	}

      }

      free(cl_integrand);
      free(primordial_pk);
      free(trsf);

      if (array_spline_table_lines(psp->l[index_mode],
				   psp->l_size[index_mode],
				   psp->cl[index_mode],
				   ppt->ic_size[index_mode]*psp->ct_size,
				   psp->ddcl[index_mode],
				   _SPLINE_EST_DERIV_,
				   errmsg) == _FAILURE_) {
	sprintf(psp->error_message,"%s(L:%d) : error in array_spline_table_lines()\n=>%s",__func__,__LINE__,errmsg);
	return _FAILURE_;
      }

    }

  }

  return _SUCCESS_;
}

int spectra_indices(){

  int index_ct;

  if (ptr->has_cls == _TRUE_) {

    index_ct=0;
    if (ppt->has_source_t == _TRUE_) {
      psp->index_ct_tt=index_ct;
      index_ct++;
    }
    if (ppt->has_source_p == _TRUE_) {
      psp->index_ct_ee=index_ct;
      index_ct++;
    }
    if ((ppt->has_source_t == _TRUE_) && (ppt->has_source_p == _TRUE_)) {
      psp->index_ct_te=index_ct;
      index_ct++;
    }
    if ((ppt->has_source_p == _TRUE_) && (ppt->has_tensors == _TRUE_)) {
      psp->index_ct_bb=index_ct;
      index_ct++;
    }
    psp->ct_size = index_ct;

  }

  return _SUCCESS_;

}

int spectra_free() {

  int index_mode;

  for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {
    free(psp->l[index_mode]);
    free(psp->cl[index_mode]);
    free(psp->ddcl[index_mode]);
  }
  free(psp->l);
  free(psp->l_size);
  free(psp->cl);
  free(psp->ddcl);

  return _SUCCESS_;
 
}
