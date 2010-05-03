/** @file cl.c Documented \f$ C_l^{X}, P(k), ... \f$ module
 * Julien Lesgourgues, 18.04.2010    
 *
 * This module computes the power spectra \f$ C_l^{X}, P(k), ... \f$'s given the transfer functions and
 * primordial spectra.
 */

#include "spectra.h"

/** @name - structures used within the transfer module: */

//@{

struct precision * ppr; /**< a precision_params structure pointer for internal use in the perturbation module */
struct background * pba; /**< a cosmo structure pointer for internal use in the thermodynamics module */
struct perturbs * ppt; /**< a perturbs structure pointer for internal use in the perturbation module */
struct transfers * ptr; /**< a transfers structure pointer for internal use in the perturbation module */
struct primordial * ppm; /**< a primordial structure pointer for internal use in the perturbation module */
struct spectra * psp; /**< a spectra structure pointer for internal use in the perturbation module */

//@}
/** @name - miscellaneous: */

//@{

ErrorMsg Transmit_Error_Message; /**< contains error message */

//@}

int spectra_cl_at_l(
		    double l,
		    int index_mode,
		    double * cl
		    ) {

  int last_index;

  if ((index_mode < 0) && (index_mode >= ppt->md_size)) {
    sprintf(psp->error_message,"%s(L:%d) : index_mode=%d out of bounds",__func__,__LINE__,index_mode);
    return _FAILURE_;
  }

  if (((int)l < ptr->l[index_mode][0]) && ((int)l >= ptr->l[index_mode][ptr->l_size[index_mode]])) {
    sprintf(psp->error_message,"%s(L:%d) : l=%d out of range [%d:%d]",
	    __func__,__LINE__,(int)l,ptr->l[index_mode][0],ptr->l[index_mode][ptr->l_size[index_mode]]);
    return _FAILURE_;
  }

  if (array_interpolate_spline(psp->l[index_mode],
			       psp->l_size[index_mode],
			       psp->cl[index_mode],
			       psp->ddcl[index_mode],
			       ppt->ic_size[index_mode]*psp->ct_size,
			       l,
			       &last_index,
			       cl,
			       ppt->ic_size[index_mode]*psp->ct_size,
			       Transmit_Error_Message) == _FAILURE_) {
    sprintf(psp->error_message,"%s(L:%d) : error in array_interpolate_spline()\n=>%s",__func__,__LINE__,Transmit_Error_Message);
    return _FAILURE_;
  }

  return _SUCCESS_;

}

int spectra_pk_at_z(
		    double z,
		    double * pk
		    ) {

  int last_index,index_mode;
  double eta_requested;

  if (background_eta_of_z(z,&eta_requested) == _FAILURE_) {
    sprintf(psp->error_message,"%s(L:%d) : error in background_at_eta()\n=>%s",__func__,__LINE__,pba->error_message);
    return _FAILURE_;
  } 

  if ((eta_requested < psp->eta[0]) || (eta_requested > psp->eta[psp->eta_size-1])) {
    sprintf(psp->error_message,"%s(L:%d) : eta(z)=%e out of bounds",__func__,__LINE__,eta_requested);
    return _FAILURE_;
  }

  index_mode=0;

  if (array_interpolate_spline(psp->eta,
			       psp->eta_size,
			       psp->pk,
			       psp->ddpk,
			       ppt->ic_size[index_mode]*psp->k_size,
			       eta_requested,
			       &last_index,
			       pk,
			       ppt->ic_size[index_mode]*psp->k_size,
			       Transmit_Error_Message) == _FAILURE_) {
    sprintf(psp->error_message,"%s(L:%d) : error in array_interpolate_spline()\n=>%s",__func__,__LINE__,Transmit_Error_Message);
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
		 struct precision * ppr_input,
		 struct background * pba_input,
		 struct perturbs * ppt_input,
		 struct transfers * ptr_input,
		 struct primordial * ppm_input,
		 struct spectra * psp_output
		 ) {

  /** - define local variables */
  int index_mode; /* index running over modes (scalar, tensor, ...) */
  int index_ic; /* index running over initial conditions */
  int index_tt; /* index running over transfer type (temperature, polarisation, ...) */
  int index_k; /* index running over wavenumber */
  int index_l;  /* multipoles */
  int index_eta;
  int index_ct;
  double k; /* wavenumber */
  double clvalue;
  int cl_integrand_num_columns;
  double * cl_integrand;
  double * transfer;
  double * primordial_pk;  /*pk[index_ic]*/
  int last_index_back;
  double * pvecback_sp_long;
  double Omega_m;

  /** - identify the spectra structure (used throughout transfer.c as global variable) to the input/output structure of this function */
  ppr = ppr_input;
  pba = pba_input; 
  ppt = ppt_input; 
  ptr = ptr_input; 
  ppm = ppm_input; 
  psp = psp_output; 

  if (psp->spectra_verbose > 0)
    printf("Computing output spectra\n");

  if (spectra_indices() == _FAILURE_) {
    sprintf(Transmit_Error_Message,"%s(L:%d) : error in spectra_indices()\n=>%s",__func__,__LINE__,psp->error_message);
    sprintf(psp->error_message,"%s",Transmit_Error_Message);
    return _FAILURE_;
  }

  /** - deal with cl's, if any */

  if (ptr->tt_size > 0) {

    if (spectra_cl() == _FAILURE_) {
      sprintf(Transmit_Error_Message,"%s(L:%d) : error in spectra_cl()\n=>%s",__func__,__LINE__,psp->error_message);
      sprintf(psp->error_message,"%s",Transmit_Error_Message);
      return _FAILURE_;
    }

  }

  /** - deal with pk's, if any */

  if (ppt->has_pk_matter == _TRUE_) {

    if (spectra_pk() == _FAILURE_) {
      sprintf(Transmit_Error_Message,"%s(L:%d) : error in spectra_pk()\n=>%s",__func__,__LINE__,psp->error_message);
      sprintf(psp->error_message,"%s",Transmit_Error_Message);
      return _FAILURE_;
    }

  }

  return _SUCCESS_;
}

int spectra_cl() {

  /** - define local variables */
  int index_mode; /* index running over modes (scalar, tensor, ...) */
  int index_ic; /* index running over initial conditions */
  int index_tt; /* index running over transfer type (temperature, polarisation, ...) */
  int index_k; /* index running over wavenumber */
  int index_l;  /* multipoles */
  int index_ct;
  double k; /* wavenumber */
  double clvalue;
  int cl_integrand_num_columns;
  double * cl_integrand;
  double * transfer;
  double * primordial_pk;  /*pk[index_ic]*/

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

  psp->cl = malloc (sizeof(double *)*ppt->md_size);
  if (psp->cl == NULL) {
    sprintf(psp->error_message,"%s(L:%d) : Could not allocate cl",__func__,__LINE__);
    return _FAILURE_;
  }

  psp->ddcl = malloc (sizeof(double *)*ppt->md_size);
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

    psp->cl[index_mode] = malloc(sizeof(double)*psp->ct_size*ppt->ic_size[index_mode]*ptr->l_size[index_mode]);
    if (psp->cl[index_mode] == NULL) {
      sprintf(psp->error_message,"%s(L:%d) : Could not allocate cl[index_mode]",__func__,__LINE__);
      return _FAILURE_;
    }

    psp->ddcl[index_mode] = malloc(sizeof(double)*psp->ct_size*ppt->ic_size[index_mode]*ptr->l_size[index_mode]);
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

    primordial_pk=malloc(ppt->ic_size[index_mode]*sizeof(double));
    if (primordial_pk == NULL) {
      sprintf(psp->error_message,"%s(L:%d) : Could not allocate primordial_pk",__func__,__LINE__);
      return _FAILURE_;
    }

    transfer=malloc(ppt->tp_size*sizeof(double));
    if (transfer == NULL) {
      sprintf(psp->error_message,"%s(L:%d) : Could not allocate transfer",__func__,__LINE__);
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
	    sprintf(psp->error_message,"%s(L:%d) : error in primordial_at_k()\n=>%s",__func__,__LINE__,ppm->error_message);
	    return _FAILURE_;
	  }

	  /* above routine checks that k>0: no possible division by zero below */

	  for (index_tt=0; index_tt < ptr->tt_size; index_tt++) {

	    transfer[index_tt] = 
	      ptr->transfer[index_mode]
	      [((index_ic * ptr->tt_size + index_tt)
		* ptr->l_size[index_mode] + index_l)
	       * ptr->k_size[index_mode] + index_k];

	  }

	  if (ppt->has_cl_cmb_temperature == _TRUE_) {
	    cl_integrand[index_k*cl_integrand_num_columns+1+psp->index_ct_tt]=primordial_pk[index_ic]
	      * transfer[ptr->index_tt_t]
	      * transfer[ptr->index_tt_t]
	      * 4. * _PI_ / k;
	  }
	      
	  if (ppt->has_cl_cmb_polarization == _TRUE_) {
	    cl_integrand[index_k*cl_integrand_num_columns+1+psp->index_ct_ee]=primordial_pk[index_ic]
	      * transfer[ptr->index_tt_p]
	      * transfer[ptr->index_tt_p]
	      * 4. * _PI_ / k;
	  }

	  if ((ppt->has_cl_cmb_temperature == _TRUE_) && (ppt->has_cl_cmb_polarization == _TRUE_)) {
	    cl_integrand[index_k*cl_integrand_num_columns+1+psp->index_ct_te]=primordial_pk[index_ic]
	      * transfer[ptr->index_tt_t]
	      * transfer[ptr->index_tt_p]
	      * 4. * _PI_ / k;
	  }

	  if ((ppt->has_cl_cmb_polarization == _TRUE_) && (ppt->has_tensors == _TRUE_)) {

	    if (index_mode == ppt->index_md_scalars) {
	      cl_integrand[index_k*cl_integrand_num_columns+1+psp->index_ct_bb]=0.;
	    }
	    if (index_mode == ppt->index_md_tensors) {
	      sprintf(psp->error_message,"%s(L:%d) : tensors not coded yet",__func__,__LINE__);
	      return _FAILURE_;
	    }
	      
	  }
	      
	  if (ppt->has_cl_cmb_lensing_potential == _TRUE_) {
	    cl_integrand[index_k*cl_integrand_num_columns+1+psp->index_ct_pp]=primordial_pk[index_ic]
	      * transfer[ptr->index_tt_lcmb]
	      * transfer[ptr->index_tt_lcmb]
	      * 4. * _PI_ / k;
	  }

	  if ((ppt->has_cl_cmb_temperature == _TRUE_) && (ppt->has_cl_cmb_lensing_potential == _TRUE_)) {
	    cl_integrand[index_k*cl_integrand_num_columns+1+psp->index_ct_tp]=primordial_pk[index_ic]
	      * transfer[ptr->index_tt_t]
	      * transfer[ptr->index_tt_lcmb]
	      * 4. * _PI_ / k;
	  }

	}

	for (index_ct=0; index_ct<psp->ct_size; index_ct++) {

	  if (array_spline(cl_integrand,
			   cl_integrand_num_columns,
			   ptr->k_size[index_mode],
			   0,
			   1+index_ct,
			   1+psp->ct_size+index_ct,
			   _SPLINE_EST_DERIV_,
			   Transmit_Error_Message) == _FAILURE_) {
	    sprintf(psp->error_message,"%s(L:%d) : error in array_spline_table_lines \n=>%s",__func__,__LINE__,Transmit_Error_Message);
	    return _FAILURE_;
	  }
	    
	  if (array_integrate_all_spline(cl_integrand,
					 cl_integrand_num_columns,
					 ptr->k_size[index_mode],
					 0,
					 1+index_ct,
					 1+psp->ct_size+index_ct,
					 &clvalue,
					 Transmit_Error_Message) == _FAILURE_) {
	    sprintf(psp->error_message,"%s(L:%d) : error in array_spline_table_lines \n=>%s",__func__,__LINE__,Transmit_Error_Message);
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

    free(transfer);

    if (array_spline_table_lines(psp->l[index_mode],
				 psp->l_size[index_mode],
				 psp->cl[index_mode],
				 ppt->ic_size[index_mode]*psp->ct_size,
				 psp->ddcl[index_mode],
				 _SPLINE_EST_DERIV_,
				 Transmit_Error_Message) == _FAILURE_) {
      sprintf(psp->error_message,"%s(L:%d) : error in array_spline_table_lines()\n=>%s",__func__,__LINE__,Transmit_Error_Message);
      return _FAILURE_;
    }

  }

  return _SUCCESS_;
}

int spectra_pk() {

  int index_mode; /* index running over modes (scalar, tensor, ...) */
  int index_ic; /* index running over initial conditions */
  int index_k; /* index running over wavenumber */
  int index_eta; /* index running over conformal time */
  double * primordial_pk;  /*pk[index_ic]*/
  int last_index_back;
  double * pvecback_sp_long;
  double Omega_m;
  double eta_min;

  if (ppt->has_scalars == _FALSE_) {
    sprintf(psp->error_message,"%s(L:%d) : you cannot ask for matter power spectrum since you turned off scalar modes",__func__,__LINE__);
    return _FAILURE_;
  }

  index_mode = ppt->index_md_scalars;

  /* if z_max_pk<0, return error */
  if (ppr->z_max_pk < 0) {
    sprintf(psp->error_message,"%s(L:%d) : aksed for z=%e, cannot compute P(k) in the future",__func__,__LINE__,ppr->z_max_pk);
    return _FAILURE_;
  }

  /* if z_max_pk=0, there is just one value to store */
  if (ppr->z_max_pk == 0.) {
    psp->eta_size=1;
  }

  /* if z_max_pk>0, store several values (with a confortable margin above z_max_pk) in view of interpolation */
  else{

    /* find the first relevant value of eta (last value in the table eta_ampling before eta(z_max)) and infer the number of vlaues of eta at which P(k) must be stored */

    if (background_eta_of_z(ppr->z_max_pk,&eta_min) == _FAILURE_) {
      sprintf(psp->error_message,"%s(L:%d) : error in background_at_eta()\n=>%s",__func__,__LINE__,pba->error_message);
      return _FAILURE_;
    }   

    index_eta=0;
    if (eta_min < ppt->eta_sampling[index_eta]) {
      sprintf(psp->error_message,"%s(L:%d) : you asked for zmax=%e, i.e. etamin=%e, smaller than first possible value =%e",__func__,__LINE__,ppr->z_max_pk,eta_min,ppt->eta_sampling[0]);
      return _FAILURE_;
    }
    while (ppt->eta_sampling[index_eta] < eta_min){
      index_eta++;
    }
    index_eta --;
    /* whenever possible, take a few more values in to avoid boundary effects in the interpolation */
    if (index_eta>0) index_eta--; 
    if (index_eta>0) index_eta--; 
    if (index_eta>0) index_eta--; 
    if (index_eta>0) index_eta--; 
    psp->eta_size=ppt->eta_size-index_eta;

  }

  /** - allocate and fill table of eta values at which P(k,eta) is stored */
  psp->eta=malloc(sizeof(double)*psp->eta_size);
  if (psp->eta == NULL) {
    sprintf(psp->error_message,"%s(L:%d) : Could not allocate eta",__func__,__LINE__);
    return _FAILURE_;
  }
  for (index_eta=0; index_eta<psp->eta_size; index_eta++) {
    psp->eta[index_eta]=ppt->eta_sampling[index_eta-psp->eta_size+ppt->eta_size];
  }

  /** - allocate and fill table of k values at which P(k,eta) is stored */
  psp->k_size = ppt->k_size[index_mode];
  psp->k = malloc(sizeof(double)*psp->k_size);
  if (psp->k == NULL) {
    sprintf(psp->error_message,"%s(L:%d) : Could not allocate k",__func__,__LINE__);
    return _FAILURE_;
  }
  for (index_k=0; index_k<psp->k_size; index_k++) {
    psp->k[index_k]=ppt->k[index_mode][index_k];
  }

  /** - allocate temporary vectors where the primordial spectrum and the background quantitites will be stored */
  primordial_pk=malloc(ppt->ic_size[index_mode]*sizeof(double));
  if (primordial_pk == NULL) {
    sprintf(psp->error_message,"%s(L:%d) : Could not allocate primordial_pk",__func__,__LINE__);
    return _FAILURE_;
  }
  pvecback_sp_long = malloc(pba->bg_size*sizeof(double));
  if (pvecback_sp_long==NULL) {
    sprintf(psp->error_message,"%s(L:%d): Cannot allocate pvecback_sp_long \n",__func__,__LINE__);
    return _FAILURE_;
  }

  /** - allocate and fill array of P(k,eta) values */
  psp->pk = malloc(sizeof(double)*psp->eta_size*psp->k_size*ppt->ic_size[index_mode]);
  if (psp->pk == NULL) {
    sprintf(psp->error_message,"%s(L:%d) : Could not allocate pk",__func__,__LINE__);
    return _FAILURE_;
  }

  for (index_eta=0 ; index_eta < psp->eta_size; index_eta++) {

    if (background_at_eta(ppt->eta_sampling[index_eta-psp->eta_size+ppt->eta_size], long_info, normal, &last_index_back, pvecback_sp_long) == _FAILURE_) {
      sprintf(psp->error_message,"%s(L:%d) : error in background_at_eta()\n=>%s",__func__,__LINE__,pba->error_message);
      return _FAILURE_;
    }  

    Omega_m = pvecback_sp_long[pba->index_bg_Omega_b];
    if (pba->has_cdm == _TRUE_) {
      Omega_m += pvecback_sp_long[pba->index_bg_Omega_cdm];
    }

    for (index_k=0; index_k<psp->k_size; index_k++) {

      if (primordial_at_k(index_mode,psp->k[index_k],primordial_pk) == _FAILURE_) {
	sprintf(psp->error_message,"%s(L:%d) : error in primordial_at_k()\n=>%s",__func__,__LINE__,ppm->error_message);
	return _FAILURE_;
      }

      /* loop over initial conditions */
      for (index_ic = 0; index_ic < ppt->ic_size[index_mode]; index_ic++) {
	
	/* primordial spectrum: 
	   P_R(k) = 1/(2pi^2) k^3 <R R>
	   so, primordial curvature correlator: 
	   <R R> = (2pi^2) k^-3 P_R(k) 
	   so, gravitational potential correlator:
	   <phi phi> = (2pi^2) k^-3 (source_phi)^2 P_R(k) 
	   so, matter power spectrum (using Poisson):
	   P(k) = <delta_m delta_m>
	        = 4/9 H^-4 Omega_m^-2 (k/a)^4 <phi phi>
	        = 4/9 H^-4 Omega_m^-2 (k/a)^4 (source_phi)^2 <R R> 
		= 8pi^2/9 H^-4 Omega_m^-2 k/a^4 (source_phi)^2 <R R> */

	psp->pk[(index_eta * ppt->ic_size[index_mode] + index_ic) * psp->k_size + index_k] =
	  8.*_PI_*_PI_/9./pow(pvecback_sp_long[pba->index_bg_H],4)/pow(Omega_m,2)*psp->k[index_k]/pow(pvecback_sp_long[pba->index_bg_a],4)
	  *primordial_pk[index_ic]
	  *pow(ppt->sources[index_mode]
	       [index_ic * ppt->tp_size + ppt->index_tp_g]
	       [(index_eta-psp->eta_size+ppt->eta_size) * ppt->k_size[index_mode] + index_k],2);
	
      }
      
    }

  }

  /* if interpolation of P(k,eta) needed (as a function of eta), spline
     the table */  
  if (psp->eta_size > 1) {

    psp->ddpk = malloc(sizeof(double)*psp->eta_size*psp->k_size*ppt->ic_size[index_mode]);
    if (psp->ddpk == NULL) {
      sprintf(psp->error_message,"%s(L:%d) : Could not allocate pk",__func__,__LINE__);
      return _FAILURE_;
    }

    if (array_spline_table_lines(psp->eta,
				 psp->eta_size,
				 psp->pk,
				 ppt->ic_size[index_mode]*psp->k_size,
				 psp->ddpk,
				 _SPLINE_EST_DERIV_,
				 Transmit_Error_Message) == _FAILURE_) {
      sprintf(psp->error_message,"%s(L:%d) : error in array_spline_table_lines()\n=>%s",__func__,__LINE__,Transmit_Error_Message);
      return _FAILURE_;
    }

  }

  free (primordial_pk);
  free (pvecback_sp_long);

  return _SUCCESS_;
}

int spectra_indices(){

  int index_ct;

  if (ptr->tt_size > 0) {

    index_ct=0;
    if (ppt->has_cl_cmb_temperature == _TRUE_) {
      psp->index_ct_tt=index_ct;
      index_ct++;
    }
    if (ppt->has_cl_cmb_polarization == _TRUE_) {
      psp->index_ct_ee=index_ct;
      index_ct++;
    }
    if ((ppt->has_cl_cmb_temperature == _TRUE_) && 
	(ppt->has_cl_cmb_polarization == _TRUE_)) {
      psp->index_ct_te=index_ct;
      index_ct++;
    }
    if ((ppt->has_cl_cmb_polarization == _TRUE_) && 
	(ppt->has_tensors == _TRUE_)) {
      psp->index_ct_bb=index_ct;
      index_ct++;
    }
    if (ppt->has_cl_cmb_lensing_potential == _TRUE_) {
      psp->index_ct_pp=index_ct;
      index_ct++;
    }
    if ((ppt->has_cl_cmb_temperature == _TRUE_) && 
	(ppt->has_cl_cmb_lensing_potential == _TRUE_)) {
      psp->index_ct_tp=index_ct;
      index_ct++;
    }
    psp->ct_size = index_ct;

  }

  return _SUCCESS_;

}

int spectra_free() {

  int index_mode;

  if (ptr->tt_size > 0) {

    for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {
      free(psp->l[index_mode]);
      free(psp->cl[index_mode]);
      free(psp->ddcl[index_mode]);
    }
    free(psp->l);
    free(psp->l_size);
    free(psp->cl);
    free(psp->ddcl);
  }

  if (ppt->has_pk_matter == _TRUE_) {

    free(psp->eta);
    free(psp->k);
    free(psp->pk);
    if (psp->eta_size > 0) {
      free(psp->ddpk);
    }
  }    

  return _SUCCESS_;
 
}
