/** @file cl.c Documented \f$ C_l^{X}, P(k), ... \f$ module
 * Julien Lesgourgues, 18.04.2010    
 *
 * This module computes the power spectra \f$ C_l^{X}, P(k), ... \f$'s given the transfer functions and
 * primordial spectra.
 */

#include "spectra.h"
#ifdef _OPENMP
#include "omp.h"
#endif

int spectra_cl_at_l(
		    struct spectra * psp,
		    double l,
		    double * cl,      /* cl[index_ct] (already allocated) */
		    double * * cl_md,   /* cl_md[index_mode][index_ct] (already allocated only if several modes) */
		    double * * cl_md_ic /* cl_md_ic[index_mode][index_ic1_ic2*psp->ct_size+index_ct] (already allocated for a given mode only if several ic's) */
		    ) {

  int last_index;
  int index_mode, index_ic1,index_ic2,index_ic1_ic2,index_ct;

    if ((psp->md_size == 1) && (psp->ic_size[0] == 1)) {
      index_mode = 0;
      if (l <= psp->l_max_tot) {
	class_call(array_interpolate_spline(psp->l[index_mode],
					    psp->l_size[index_mode],
					    psp->cl[index_mode],
					    psp->ddcl[index_mode],
					    psp->ic_ic_size[index_mode]*psp->ct_size,
					    l,
					    &last_index,
					    cl,
					    psp->ic_size[index_mode]*psp->ct_size,
					    psp->error_message),
		   psp->error_message,
		   psp->error_message);
      }
      else {
	for (index_ct=0; index_ct<psp->ct_size; index_ct++) 
	  cl[index_ct]=0.;
      }
    return _SUCCESS_;
  }
  if (psp->md_size > 1) {
    for (index_ct=0; index_ct<psp->ct_size; index_ct++) 
      cl[index_ct]=0.;
    for (index_mode = 0; index_mode < psp->md_size; index_mode++) {
      if (psp->ic_size[index_mode] == 1) {
	if (l <= psp->l_max[index_mode]) {
	  class_call(array_interpolate_spline(psp->l[index_mode],
					      psp->l_size[index_mode],
					      psp->cl[index_mode],
					      psp->ddcl[index_mode],
					      psp->ic_ic_size[index_mode]*psp->ct_size,
					      l,
					      &last_index,
					      cl_md[index_mode],
					      psp->ic_size[index_mode]*psp->ct_size,
					      psp->error_message),
		     psp->error_message,
		     psp->error_message);
	}
	else {
	  for (index_ct=0; index_ct<psp->ct_size; index_ct++) 
	    cl_md[index_mode][index_ct]=0.;
	}
      }
      else {
	for (index_ct=0; index_ct<psp->ct_size; index_ct++) 
	  cl_md[index_mode][index_ct]=0.;
	for (index_ic1 = 0; index_ic1 < psp->ic_size[index_mode]; index_ic1++) {
	  for (index_ic2 = index_ic1; index_ic2 < psp->ic_size[index_mode]; index_ic2++) {
	    /* index value for the coefficients of the symmetric index_ic1*index_ic2 matrix; 
	       takes values between 0 and N(N+1)/2-1 with N=ppt->ic_size[index_mode] */
	    index_ic1_ic2 = index_ic1 + psp->ic_size[index_mode]*index_ic2 - (index_ic2*(index_ic2+1))/2;
	    if ((l <= psp->l_max[index_mode]) && 
		((index_ic1 == index_ic2) || (psp->is_non_zero[index_mode][index_ic1_ic2] == _TRUE_))) {
	      class_call(array_interpolate_spline(psp->l[index_mode],
						  psp->l_size[index_mode],
						  psp->cl[index_mode],
						  psp->ddcl[index_mode],
						  psp->ic_ic_size[index_mode]*psp->ct_size,
						  l,
						  &last_index,
						  cl_md_ic[index_mode],
						  psp->ic_size[index_mode]*psp->ct_size,
						  psp->error_message),
			 psp->error_message,
			 psp->error_message);
	    }
	    else {
	      for (index_ct=0; index_ct<psp->ct_size; index_ct++)
		cl_md_ic[index_mode][index_ic1_ic2*psp->ct_size+index_ct]=0.;
	    }
	    for (index_ct=0; index_ct<psp->ct_size; index_ct++) {
	      if (index_ic1 == index_ic2)
		cl_md[index_mode][index_ct]+=cl_md_ic[index_mode][index_ic1_ic2*psp->ct_size+index_ct];
	      else
		cl_md[index_mode][index_ct]+=2.*cl_md_ic[index_mode][index_ic1_ic2*psp->ct_size+index_ct];
	    }
	  }
	}
      }
      for (index_ct=0; index_ct<psp->ct_size; index_ct++) 
	cl[index_ct]+=cl_md[index_mode][index_ct];
    }
  }

  return _SUCCESS_;

}

/**
 * Interpolate the spectrum at an arbitrary redhsifts, for all tabulated values of k between kmin and kmax.
 *
 * The output vector pk[] is in the format:
 *
 *     pk[index_ic * psp->k_size + index_k]
 *
 * with index_ic=0, ..., ppt->ic_size[index_mode]
 *      index_k =0, ..., psp->k_size
 *
 * it corresponds physically to the values of the power spectrum P(psp->k[index_k]) 
 * for the initial condtions (adiabatic, cdi, etc...) of index_ic
 *
 * k[] is in units of 1/Mpc, pk[] in Mpc^3
 *
 */
int spectra_pk_at_z(
		    struct background * pba,
		    struct spectra * psp,
		    int index_mode,
		    double z,
		    double * pk
		    ) {

  int last_index;
  double eta_requested;
  int index_ic;

  class_call(background_eta_of_z(pba,z,&eta_requested),
	     pba->error_message,
	     psp->error_message);

  /* interpolation makes sense only if there are at least two values of eta. Deal with case of one value. */
  if (psp->eta_size == 1) {

    class_test(z!=0.,
	       psp->error_message,
	       "asked z=%e but only P(k,z=0) has been tabulated",z);

    /* case z=0 */
    for (index_ic=0; index_ic < psp->ic_size[index_mode]*psp->k_size; index_ic++) {
      pk[index_ic] = psp->pk[index_ic];
    }
    return _SUCCESS_;
    
  }

  class_test((eta_requested < psp->eta[0]) || (eta_requested > psp->eta[psp->eta_size-1]),
	     psp->error_message,
	     "eta(z)=%e out of bounds [%e:%e]",
	     eta_requested,psp->eta[0],psp->eta[psp->eta_size-1]);


  class_call(array_interpolate_logspline(psp->eta,
				      psp->eta_size,
				      psp->pk,
				      psp->ddlnpk,
				      psp->ic_size[index_mode]*psp->k_size,
				      eta_requested,
				      &last_index,
				      pk,
				      psp->ic_size[index_mode]*psp->k_size,
				      psp->error_message),
	     psp->error_message,
	     psp->error_message);

  return _SUCCESS_;

}

/**
 * Interpolate the spectrum at an arbitrary redhsifts and wavenumber between kmin and kmax
 * (in the case of an analytic primordial spectrum, k may be between 0 and kmax: extrapolation 
 *  is performed for suoper-hubble scales down to k=0, given the values of the tilt, running, etc.)
 *
 *  @param k Input : wavenumber k in units of 1/Mpc 
 *  @param z Input : redhsift z
 *  @param index_ic : index of initial condition (ppt->index_ic_ad for adiabatic, etc...)
 *  @param pk Ouput : "P(k,z) in units of Mpc^3  
 */
int spectra_pk_at_k_and_z(
			  struct background * pba,
			  struct primordial * ppm,
			  struct spectra * psp,
			  int index_mode,
			  int index_ic,
			  double k,
			  double z,
			  double * pk
			  ) {

  double * temporary_pk;
  double * spline_pk;
  double * spline_ddlnpk;
  int index_k;
  int last_index;
  double * pkini_k;
  double * pkini_kmin;

  /* check input parameters (z will be checked in spectra_pk_at_z) */

  class_test((index_ic < 0) || (index_ic >= psp->ic_size[index_mode]),
	     psp->error_message,
	     "index_ic=%d out of bounds [%d:%d]",
	     index_ic,0,psp->ic_size[index_mode]);

  class_test((k < 0) || (k > psp->k[psp->k_size-1]),
	     psp->error_message,
	     "k=%e out of bounds [%e:%e]",k,0.,psp->k[psp->k_size-1]);

  /* get P(k) at the right value of z */

  class_alloc(temporary_pk,
	      sizeof(double)*psp->ic_size[index_mode]*psp->k_size,
	      psp->error_message);

  class_call(spectra_pk_at_z(pba,psp,index_mode,z,temporary_pk),
	     psp->error_message,
	     psp->error_message);

  /* deal with case k < kmin */

  if (k < psp->k[0]) {

    /* extrapolate down to zero if analytical primordial spectrum: in this case we know that on super-Hubble scales:        
     * P(k) = A k P_ini(k) 
     * so     
     * P(k) = P(kmin) * (k P_ini(k)) / (kmin P_ini(kmin)) 
     */

    class_test((ppm->primordial_spec_type != analytic_Pk),
	       psp->error_message,
	       "in this case, exptrapolation for k=%e below k_min=%e is impossible",k,psp->k[0]);
      
    if (k == 0.) {
      *pk=0.;
      return _SUCCESS_;
    }
	
    class_alloc(pkini_k,sizeof(double)*psp->ic_size[index_mode],psp->error_message);
    class_alloc(pkini_kmin,sizeof(double)*psp->ic_size[index_mode],psp->error_message);

    class_call(primordial_spectrum_at_k(ppm,index_mode,k,pkini_k),ppm->error_message,psp->error_message);
    class_call(primordial_spectrum_at_k(ppm,index_mode,psp->k[0],pkini_kmin),ppm->error_message,psp->error_message);
    
    *pk=temporary_pk[index_ic*psp->k_size]*k*pkini_k[index_ic]/psp->k[0]/pkini_kmin[index_ic];
    
    free(temporary_pk);
    free(pkini_k);
    free(pkini_kmin);
    
    return _SUCCESS_;

  }

  /* now, spline and interpolate */

  class_alloc(spline_pk,sizeof(double)*psp->k_size,psp->error_message);

  for (index_k=0; index_k < psp->k_size; index_k++) {
    spline_pk[index_k] = temporary_pk[index_ic*psp->k_size+index_k];
  }

  free(temporary_pk);
  
  class_alloc(spline_ddlnpk,sizeof(double)*psp->k_size,psp->error_message);

  class_call(array_logspline_table_lines(psp->k,
				      psp->k_size,
				      spline_pk,
				      1,
				      spline_ddlnpk,
				      _SPLINE_NATURAL_,
				      psp->error_message),
	     psp->error_message,
	     psp->error_message);

  class_call(array_interpolate_logspline(psp->k,
				      psp->k_size,
				      spline_pk,
				      spline_ddlnpk,
				      1,
				      k,
				      &last_index,
				      pk,
				      1,
				      psp->error_message),
	     psp->error_message,
	     psp->error_message);

  free(spline_pk);
  free(spline_ddlnpk);

  return _SUCCESS_;

}


/**
 * Computes the \f$ C_l^{X}, P(k), ... \f$'s
 *
 * @param ppt Input : Initialized perturbation structure
 * @param ptr Input : Initialized transfers structure
 * @param ppm Input : Initialized primordial structure
 * @param pcl Output : Initialized cls structure
 * @return the error status
 */
int spectra_init(
		 struct background * pba,
		 struct perturbs * ppt,
		 struct transfers * ptr,
		 struct primordial * ppm,
		 struct spectra * psp
		 ) {

  if (ppt->has_perturbations == _FALSE_) {
    psp->md_size = 0;
    if (psp->spectra_verbose > 0)
      printf("No spectra requested. Spectra module skipped.\n");
    return _SUCCESS_;
  }
  else {
    if (psp->spectra_verbose > 0)
      printf("Computing output spectra\n");
  }

  class_call(spectra_indices(ppt,ptr,psp),
	     psp->error_message,
	     psp->error_message);

  /** - deal with cl's, if any */

  if (ppt->has_cls == _TRUE_) {

    class_call(spectra_cls(ppt,ptr,ppm,psp),
	       psp->error_message,
	       psp->error_message);

  }
  else {
    psp->ct_size=0;
  }

  /** - deal with pk's, if any */

  if (ppt->has_pk_matter == _TRUE_) {

    class_call(spectra_pk(pba,ppt,ppm,psp),
	       psp->error_message,
	       psp->error_message);
  }
  else {
    psp->k_size=0;
  }

  return _SUCCESS_;
}

int spectra_free(
		 struct spectra * psp
		 ) {

  int index_mode;

  if (psp->md_size > 0) {

    if (psp->ct_size > 0) {
    
      for (index_mode = 0; index_mode < psp->md_size; index_mode++) {
	free(psp->l[index_mode]);
	free(psp->cl[index_mode]);
	free(psp->ddcl[index_mode]);
      }
      free(psp->l);
      free(psp->l_size);
      free(psp->l_max);
      free(psp->cl);
      free(psp->ddcl);
    }

    if (psp->k_size > 0) {

      free(psp->eta);
      free(psp->k);
      free(psp->pk);
      if (psp->eta_size > 1) {
	free(psp->ddlnpk);
      }
    }    

  }

  for (index_mode=0; index_mode < psp->md_size; index_mode++)
    free(psp->is_non_zero[index_mode]);
  free(psp->is_non_zero);
  free(psp->ic_size);
  free(psp->ic_ic_size);

  return _SUCCESS_;
 
}

int spectra_indices(
		    struct perturbs * ppt,
		    struct transfers * ptr,
		    struct spectra * psp
		    ){

  int index_ct;
  int index_mode;

  psp->md_size = ppt->md_size;

  class_alloc(psp->ic_size,
	      sizeof(int)*psp->md_size,
	      psp->error_message);

  class_alloc(psp->ic_ic_size,
	      sizeof(int)*psp->md_size,
	      psp->error_message);

  class_alloc(psp->is_non_zero,
	      sizeof(short *)*psp->md_size,
	      psp->error_message);

  for (index_mode=0; index_mode < psp->md_size; index_mode++) {
    psp->ic_size[index_mode] = ppt->ic_size[index_mode];
    psp->ic_ic_size[index_mode] = (psp->ic_size[index_mode]*(psp->ic_size[index_mode]+1))/2;
    class_alloc(psp->is_non_zero[index_mode],
		sizeof(short)*psp->ic_ic_size[index_mode],
		psp->error_message);
  }

  if (ppt->has_cls == _TRUE_) {

    /* types of C_l's relevant for both scalars and tensors: TT, EE, TE */

    index_ct=0;

    if (ppt->has_cl_cmb_temperature == _TRUE_) {
      psp->has_tt = _TRUE_;
      psp->index_ct_tt=index_ct;      
      index_ct++;
    }
    else {
      psp->has_tt = _FALSE_;
    }

    if (ppt->has_cl_cmb_polarization == _TRUE_) {
      psp->has_ee = _TRUE_;
      psp->index_ct_ee=index_ct;
      index_ct++;
    }
    else {
      psp->has_ee = _FALSE_;
    }

    if ((ppt->has_cl_cmb_temperature == _TRUE_) && 
	(ppt->has_cl_cmb_polarization == _TRUE_)) {
      psp->has_te = _TRUE_;
      psp->index_ct_te=index_ct;      
      index_ct++;
    }
    else {
      psp->has_te = _FALSE_;
    }

    /* types of C_l's relevant only for tensors: BB */

    if ((ppt->has_cl_cmb_polarization == _TRUE_) && (ppt->has_tensors == _TRUE_)) {
      psp->has_bb = _TRUE_;
      psp->index_ct_bb=index_ct;
      index_ct++;
    }
    else {
      psp->has_bb = _FALSE_;
    }

    /* types of C_l's relevant only for scalars: phi-phi, T-phi */
    
    if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
      psp->has_pp = _TRUE_;
      psp->index_ct_pp=index_ct;
      index_ct++;
    }
    else {
      psp->has_pp = _FALSE_;
    }

    if ((ppt->has_cl_cmb_temperature == _TRUE_) && (ppt->has_cl_cmb_lensing_potential == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
      psp->has_tp = _TRUE_;
      psp->index_ct_tp=index_ct;
      index_ct++;
    }
    else {
      psp->has_tp = _FALSE_;
    }

    psp->ct_size = index_ct;

  }

  return _SUCCESS_;

}

int spectra_cls(
		struct perturbs * ppt,
		struct transfers * ptr,
		struct primordial * ppm,
		struct spectra * psp
		) {

  /** - define local variables */
  int index_mode; /* index running over modes (scalar, tensor, ...) */
  int index_ic1,index_ic2,index_ic1_ic2; /* index running over initial conditions */
  int index_l;  /* multipoles */
  int index_ct;
  int cl_integrand_num_columns;
  double * cl_integrand;
  double * transfer_ic1;
  double * transfer_ic2;
  double * primordial_pk;  /*pk[index_ic]*/

  /* This code can be optionally compiled with the openmp option for parallel computation.
     Inside parallel regions, the use of the command "return" is forbidden.
     For error management, instead of "return _FAILURE_", we will set the variable below
     to "abort = _TRUE_". This will lead to a "return _FAILURE_" jus after leaving the 
     parallel region. */
  int abort;

#ifdef _OPENMP
  /* instrumentation times */
  double tstart, tstop;
#endif

  class_alloc(psp->l_size,sizeof(int)*psp->md_size,psp->error_message);
  class_alloc(psp->l,sizeof(double *)*psp->md_size,psp->error_message);
  class_alloc(psp->l_max,sizeof(int)*psp->md_size,psp->error_message);
  class_alloc(psp->cl,sizeof(double *)*psp->md_size,psp->error_message);
  class_alloc(psp->ddcl,sizeof(double *)*psp->md_size,psp->error_message);

  /** - loop over modes (scalar, tensors, etc) */
  for (index_mode = 0; index_mode < psp->md_size; index_mode++) {

    psp->l_size[index_mode] = ptr->l_size[index_mode];

    class_alloc(psp->l[index_mode],sizeof(double)*psp->l_size[index_mode],psp->error_message);

    for (index_l=0; index_l < psp->l_size[index_mode]; index_l++) {
      psp->l[index_mode][index_l] = (double)ptr->l[index_mode][index_l];
    }

    class_alloc(psp->cl[index_mode],sizeof(double)*psp->l_size[index_mode]*psp->ct_size*psp->ic_ic_size[index_mode],psp->error_message);
    class_alloc(psp->ddcl[index_mode],sizeof(double)*psp->l_size[index_mode]*psp->ct_size*psp->ic_ic_size[index_mode],psp->error_message);
    cl_integrand_num_columns = 1+psp->ct_size*2; /* one for k, ct_size for each type, ct_size for each second derivative of each type */

    /* last multipole for each mode (given as an input) at which we trust our C_ls;
       (l[index_mode][l_size[index_mode]-1] can be larger than l_max[index_mode], 
       in order to ensure better interpolation with no boundary effects).
       Compute also the max over all modes */
    psp->l_max_tot=0;
    if ((ppt->has_scalars) && (index_mode == ppt->index_md_scalars)) {
      psp->l_max[index_mode] = ptr->l_scalar_max;
      if (psp->l_max[index_mode] > psp->l_max_tot) psp->l_max_tot=psp->l_max[index_mode];
    }
    
    if ((ppt->has_tensors) && (index_mode == ppt->index_md_tensors)) {
      psp->l_max[index_mode] = ptr->l_tensor_max;
      if (psp->l_max[index_mode] > psp->l_max_tot) psp->l_max_tot=psp->l_max[index_mode];
    }

    /** - loop over initial conditions */
    for (index_ic1 = 0; index_ic1 < psp->ic_size[index_mode]; index_ic1++) {
      for (index_ic2 = index_ic1; index_ic2 < psp->ic_size[index_mode]; index_ic2++) {

	/* non-diagonal coefficients should be computed only if non-zero correlation */
	if ((index_ic1 == index_ic2) || (ppm->has_correlation[index_mode][index_ic1*ppm->ic_size[index_mode]+index_ic2] == _TRUE_)) {

	  /* initialize error management flag */
	  abort = _FALSE_;

	  /*** beginning of parallel region ***/

#pragma omp parallel							\
  shared(ptr,ppm,index_mode,psp,ppt,cl_integrand_num_columns,index_ic1,index_ic2,abort) \
  private(tstart,cl_integrand,primordial_pk,transfer_ic1,transfer_ic2,index_l,tstop)
      
	  {

#ifdef _OPENMP
	    tstart = omp_get_wtime();
#endif

	    class_alloc_parallel(cl_integrand,
				 ptr->k_size[index_mode]*cl_integrand_num_columns*sizeof(double),
				 psp->error_message);
	    
	    class_alloc_parallel(primordial_pk,
				 psp->ic_ic_size[index_mode]*sizeof(double),
				 psp->error_message);
	    
	    class_alloc_parallel(transfer_ic1,
				 ptr->tt_size[index_mode]*sizeof(double),
				 psp->error_message);
	    
	    class_alloc_parallel(transfer_ic2,
				 ptr->tt_size[index_mode]*sizeof(double),
				 psp->error_message);

#pragma omp for schedule (dynamic)

	    /** - loop over l values defined in the transfer module. For each l: */
	    for (index_l=0; index_l < ptr->l_size[index_mode]; index_l++) {

#pragma omp flush(abort)

	      class_call_parallel(spectra_compute_cl(ppt,
						     ptr,
						     ppm,
						     psp,
						     index_mode,
						     index_ic1,
						     index_ic2,
						     index_l,
						     cl_integrand_num_columns,
						     cl_integrand,
						     primordial_pk,
						     transfer_ic1,
						     transfer_ic2),
				  psp->error_message,
				  psp->error_message);

	    } /* end of loop over l */
	
#ifdef _OPENMP
	    tstop = omp_get_wtime();
	    if (psp->spectra_verbose > 1)
	      printf("In %s: time spent in parallel region (loop over l's) = %e s for thread %d\n",
		     __func__,tstop-tstart,omp_get_thread_num());
#endif
	    free(cl_integrand);
    
	    free(primordial_pk);
	
	    free(transfer_ic1);
	  
	    free(transfer_ic2);

	  } /* end of parallel region */

	  if (abort == _TRUE_) return _FAILURE_;

	}
	else {
          /* set non-diagonal coefficients to zero if pair of ic's uncorrelated */
	  index_ic1_ic2 = index_ic1 + ppt->ic_size[index_mode]*index_ic2 - (index_ic2*(index_ic2+1))/2;
	  for (index_l=0; index_l < ptr->l_size[index_mode]; index_l++) {
	    for (index_ct=0; index_ct<psp->ct_size; index_ct++) {
	      psp->cl[index_mode]
		[(index_l * psp->ic_ic_size[index_mode] + index_ic1_ic2) * psp->ct_size + index_ct]
		= 0.;
	      psp->is_non_zero[index_mode][index_ic1_ic2] = _FALSE_;

	    }
	  }
	}
      }
    }

    class_call(array_spline_table_lines(psp->l[index_mode],
					psp->l_size[index_mode],
					psp->cl[index_mode],
					psp->ic_ic_size[index_mode]*psp->ct_size,
					psp->ddcl[index_mode],
					_SPLINE_EST_DERIV_,
					psp->error_message),
	       psp->error_message,
	       psp->error_message);
  }

  return _SUCCESS_;

}

int spectra_compute_cl(
		       struct perturbs * ppt,
		       struct transfers * ptr,
		       struct primordial * ppm,
		       struct spectra * psp,
		       int index_mode,
		       int index_ic1,
		       int index_ic2,
		       int index_l,
		       int cl_integrand_num_columns,
		       double * cl_integrand,
		       double * primordial_pk,
		       double * transfer_ic1,
		       double * transfer_ic2
		       ) {

  int index_k;
  int index_tt;
  int index_ct;
  double k;
  double clvalue;
  int nonzero;
  int index_ic1_ic2;

  /* index value for the coefficients of the symmetric index_ic1*index_ic2 matrix; 
     takes values between 0 and N(N+1)/2-1 with N=ppt->ic_size[index_mode] */
  index_ic1_ic2 = index_ic1 + ppt->ic_size[index_mode]*index_ic2 - (index_ic2*(index_ic2+1))/2; 

  for (index_k=0; index_k < ptr->k_size[index_mode]; index_k++) {

    k = ptr->k[index_mode][index_k];

    cl_integrand[index_k*cl_integrand_num_columns+0] = k;

    class_call(primordial_spectrum_at_k(ppm,index_mode,k,primordial_pk),
	       ppm->error_message,
	       psp->error_message);

    /* above routine checks that k>0: no possible division by zero below */

    for (index_tt=0; index_tt < ptr->tt_size[index_mode]; index_tt++) {
		  
      transfer_ic1[index_tt] = 
	ptr->transfer[index_mode]
	[((index_ic1 * ptr->tt_size[index_mode] + index_tt)
	  * ptr->l_size[index_mode] + index_l)
	 * ptr->k_size[index_mode] + index_k];
      
      if (index_ic1 == index_ic2) {
	transfer_ic2[index_tt] = transfer_ic1[index_tt];
      }
      else {
	transfer_ic2[index_tt] = ptr->transfer[index_mode]
	  [((index_ic2 * ptr->tt_size[index_mode] + index_tt)
	    * ptr->l_size[index_mode] + index_l)
	   * ptr->k_size[index_mode] + index_k];
      }
    }

    if (psp->has_tt == _TRUE_)
      cl_integrand[index_k*cl_integrand_num_columns+1+psp->index_ct_tt]=
	primordial_pk[index_ic1_ic2]
	* transfer_ic1[ptr->index_tt_t]
	* transfer_ic2[ptr->index_tt_t]
	* 4. * _PI_ / k;
		  
    if (psp->has_ee == _TRUE_)
      cl_integrand[index_k*cl_integrand_num_columns+1+psp->index_ct_ee]=
	primordial_pk[index_ic1_ic2]
	* transfer_ic1[ptr->index_tt_e]
	* transfer_ic2[ptr->index_tt_e]
	* 4. * _PI_ / k;

    if (psp->has_te == _TRUE_)
      cl_integrand[index_k*cl_integrand_num_columns+1+psp->index_ct_te]=
	primordial_pk[index_ic1_ic2]
	* 0.5*(transfer_ic1[ptr->index_tt_t] * transfer_ic2[ptr->index_tt_e] +
	       transfer_ic1[ptr->index_tt_e] * transfer_ic2[ptr->index_tt_t])
	* 4. * _PI_ / k;
    
    if ((psp->has_bb == _TRUE_) && (ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors))
	cl_integrand[index_k*cl_integrand_num_columns+1+psp->index_ct_bb]=
	  primordial_pk[index_ic1_ic2]
	* transfer_ic1[ptr->index_tt_b]
	* transfer_ic2[ptr->index_tt_b]
	* 4. * _PI_ / k;

    if ((psp->has_pp == _TRUE_) && (ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars))
      cl_integrand[index_k*cl_integrand_num_columns+1+psp->index_ct_pp]=
	primordial_pk[index_ic1_ic2]
	* transfer_ic1[ptr->index_tt_lcmb]
	* transfer_ic2[ptr->index_tt_lcmb]
	* 4. * _PI_ / k;
    
    if ((psp->has_tp == _TRUE_) && (ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars))
      cl_integrand[index_k*cl_integrand_num_columns+1+psp->index_ct_tp]=
	primordial_pk[index_ic1_ic2]
	* 0.5*(transfer_ic1[ptr->index_tt_t] * transfer_ic2[ptr->index_tt_lcmb] +
	       transfer_ic1[ptr->index_tt_lcmb] * transfer_ic2[ptr->index_tt_t])
	* 4. * _PI_ / k;

  }
  
  for (index_ct=0; index_ct<psp->ct_size; index_ct++) {

    /* treat null spectra (C_l^BB of scalars, C_l^pp of tensors, etc. */

    if (((index_ct == psp->index_ct_bb) && (psp->has_bb == _TRUE_) && (ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) ||
	((index_ct == psp->index_ct_pp) && (psp->has_pp == _TRUE_) && (ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors)) ||
	((index_ct == psp->index_ct_tp) && (psp->has_tp == _TRUE_) && (ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors))) {

      index_ic1_ic2 = index_ic1 + psp->ic_size[index_mode]*index_ic2 - (index_ic2*(index_ic2-1))/2;

      psp->cl[index_mode]
	[(index_l * psp->ic_ic_size[index_mode] + index_ic1_ic2) * psp->ct_size + index_ct] = 0.;

    }
    /* for non-zero spectra, integrate over k */
    else {

      class_call(array_spline(cl_integrand,
			      cl_integrand_num_columns,
			      ptr->k_size[index_mode],
			      0,
			      1+index_ct,
			      1+psp->ct_size+index_ct,
			      _SPLINE_EST_DERIV_,
			      psp->error_message),
		 psp->error_message,
		 psp->error_message);
      
      class_call(array_integrate_all_spline(cl_integrand,
					    cl_integrand_num_columns,
					    ptr->k_size[index_mode],
					    0,
					    1+index_ct,
					    1+psp->ct_size+index_ct,
					    &clvalue,
					    psp->error_message),
		 psp->error_message,
		 psp->error_message);
      
      psp->cl[index_mode]
	[(index_l * psp->ic_ic_size[index_mode] + index_ic1_ic2) * psp->ct_size + index_ct]
	= clvalue;
      psp->is_non_zero[index_mode][index_ic1_ic2] = _TRUE_;

      fprintf(stderr,"%d %d %d %d %d %g %g\n",index_mode,index_ic1_ic2,index_l,index_ct,
	      (index_l * psp->ic_ic_size[index_mode] + index_ic1_ic2) * psp->ct_size + index_ct,
	      psp->l[index_mode][index_l],
	      psp->cl[index_mode][(index_l * psp->ic_ic_size[index_mode] + index_ic1_ic2) * psp->ct_size + index_ct]);

    }

  }

  return _SUCCESS_;

}

int spectra_pk(
	       struct background * pba,
	       struct perturbs * ppt,
	       struct primordial * ppm,
	       struct spectra * psp
	       ) {

  int index_mode; /* index running over modes (scalar, tensor, ...) */
  int index_ic; /* index running over initial conditions */
  int index_k; /* index running over wavenumber */
  int index_eta; /* index running over conformal time */
  double * primordial_pk;  /*pk[index_ic]*/
  int last_index_back;
  double * pvecback_sp_long;
  double Omega_m;
  double eta_min;

  class_test((ppt->has_scalars == _FALSE_),
	     psp->error_message,
	     "you cannot ask for matter power spectrum since you turned off scalar modes");

  index_mode = ppt->index_md_scalars;

  /* if z_max_pk<0, return error */
  class_test((psp->z_max_pk < 0),
	     psp->error_message,
	     "asked for z=%e. Code not designed to compute anything in the future...",psp->z_max_pk);

  /* if z_max_pk=0, there is just one value to store */
  if (psp->z_max_pk == 0.) {
    psp->eta_size=1;
  }

  /* if z_max_pk>0, store several values (with a confortable margin above z_max_pk) in view of interpolation */
  else{

    /* find the first relevant value of eta (last value in the table eta_ampling before eta(z_max)) and infer the number of values of eta at which P(k) must be stored */

    class_call(background_eta_of_z(pba,psp->z_max_pk,&eta_min),
	       pba->error_message,
	       psp->error_message);

    index_eta=0;
    class_test((eta_min < ppt->eta_sampling[index_eta]),
	       psp->error_message,
	       "you asked for zmax=%e, i.e. etamin=%e, smaller than first possible value =%e",psp->z_max_pk,eta_min,ppt->eta_sampling[0]);
  
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
  class_alloc(psp->eta,sizeof(double)*psp->eta_size,psp->error_message);

  for (index_eta=0; index_eta<psp->eta_size; index_eta++) {
    psp->eta[index_eta]=ppt->eta_sampling[index_eta-psp->eta_size+ppt->eta_size];
  }

  /** - allocate and fill table of k values at which P(k,eta) is stored */
  psp->k_size = ppt->k_size[index_mode];
  class_alloc(psp->k,sizeof(double)*psp->k_size,psp->error_message);

  for (index_k=0; index_k<psp->k_size; index_k++) {
    psp->k[index_k]=ppt->k[index_mode][index_k];
  }

  /** - allocate temporary vectors where the primordial spectrum and the background quantitites will be stored */
  class_alloc(primordial_pk,psp->ic_size[index_mode]*sizeof(double),psp->error_message);
  class_alloc(pvecback_sp_long,pba->bg_size*sizeof(double),psp->error_message);

  /** - allocate and fill array of P(k,eta) values */
  class_alloc(psp->pk,sizeof(double)*psp->eta_size*psp->k_size*psp->ic_size[index_mode],psp->error_message);

  for (index_eta=0 ; index_eta < psp->eta_size; index_eta++) {

    class_call(background_at_eta(pba,ppt->eta_sampling[index_eta-psp->eta_size+ppt->eta_size], 
				 long_info, 
				 normal, 
				 &last_index_back, 
				 pvecback_sp_long),
	       pba->error_message,
	       psp->error_message);

    Omega_m = pvecback_sp_long[pba->index_bg_Omega_b];
    if (pba->has_cdm == _TRUE_) {
      Omega_m += pvecback_sp_long[pba->index_bg_Omega_cdm];
    }

    for (index_k=0; index_k<psp->k_size; index_k++) {

      class_call(primordial_spectrum_at_k(ppm,index_mode,psp->k[index_k],primordial_pk),
		 ppm->error_message,
		 psp->error_message);

      /* loop over initial conditions */
      for (index_ic = 0; index_ic < psp->ic_size[index_mode]; index_ic++) {
	
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

	psp->pk[(index_eta * psp->ic_size[index_mode] + index_ic) * psp->k_size + index_k] =
	  8.*_PI_*_PI_/9./pow(pvecback_sp_long[pba->index_bg_H],4)/pow(Omega_m,2)*psp->k[index_k]/pow(pvecback_sp_long[pba->index_bg_a],4)
	  *primordial_pk[index_ic]
	  *pow(ppt->sources[index_mode]
	       [index_ic * ppt->tp_size[index_mode] + ppt->index_tp_g]
	       [(index_eta-psp->eta_size+ppt->eta_size) * ppt->k_size[index_mode] + index_k],2);
	
      }
      
    }

  }

  /* if interpolation of P(k,eta) needed (as a function of eta), spline
     the table */  
  if (psp->eta_size > 1) {

    class_alloc(psp->ddlnpk,sizeof(double)*psp->eta_size*psp->k_size*psp->ic_size[index_mode],psp->error_message);

    class_call(array_logspline_table_lines(psp->eta,
					psp->eta_size,
					psp->pk,
					psp->ic_size[index_mode]*psp->k_size,
					psp->ddlnpk,
					_SPLINE_EST_DERIV_,
					psp->error_message),
	       psp->error_message,
	       psp->error_message);

  }

  free (primordial_pk);
  free (pvecback_sp_long);

  return _SUCCESS_;
}
