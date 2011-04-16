/** @file cl.c Documented spectra module
 *
 * Julien Lesgourgues, 25.08.2010    
 *
 * This module computes the anisotropy and Fourier power spectra 
 * \f$ C_l^{X}, P(k), ... \f$'s given the transfer and Bessel functions 
 * (for anisotropy spectra), the source functions (for Fourier spectra) 
 * and the primordial spectra.
 *
 * The following functions can be called from other modules:
 *
 * -# spectra_init() at the beginning (but after transfer_init())
 * -# spectra_cl_at_l() at any time for computing C at any l
 * -# spectra_spectrum_at_z() at any time for computing P(k) at any z
 * -# spectra_spectrum_at_k_and z() at any time for computing P at any k and z
 * -# spectra_free() at the end
 */

#include "spectra.h"

/** 
 * Anisotropy power spectra C_l's for all types, modes and initial conditions. 
 *
 * This routine evaluates all the C_l's at a given value of l by
 * interpolating in the pre-computed table. When relevant, it also
 * sums over all initial conditions for each mode, and over all modes.
 * 
 * This function can be
 * called from whatever module at whatever time, provided that
 * spectra_init() has been called before, and spectra_free() has not
 * been called yet.
 *
 * @param psp        Input: pointer to spectra structure (containing pre-computed table)
 * @param l          Input: multipole number 
 * @param cl_tot     Ouput: total C_l's for all types (TT, TE, EE, etc..)
 * @param cl_md      Ouput: C_l's for all types (TT, TE, EE, etc..) decomposed mode by mode (scalar, tensor, ...) when relevant
 * @param cl_md_ic   Ouput: C_l's for all types (TT, TE, EE, etc..) decomposed by pairs of initial conditions (adiabatic, isocurvatures) for each mode (usually, only for the scalar mode) when relevant
 * @return the error status
 */

int spectra_cl_at_l(
		    struct spectra * psp,
		    double l,
		    double * cl_tot,    /* array with argument cl_tot[index_ct] (must be already allocated) */
		    double * * cl_md,   /* array with argument cl_md[index_mode][index_ct] (must be already allocated only if several modes) */
		    double * * cl_md_ic /* array with argument cl_md_ic[index_mode][index_ic1_ic2*psp->ct_size+index_ct] (must be already allocated for a given mode only if several ic's) */
		    ) {

  /** Summary: */

  /** - define local variables */

  int last_index;
  int index_mode;
  int index_ic1,index_ic2,index_ic1_ic2;
  int index_ct;

  /** A) treat case in which there is only one mode and one initial condition. 
         Then, only cl_tot needs to be filled. */

  if ((psp->md_size == 1) && (psp->ic_size[0] == 1)) {
    index_mode = 0;
    if ((int)l <= psp->l_max_tot) {

      class_call(array_interpolate_spline(psp->l,
					  psp->l_size[index_mode],
					  psp->cl[index_mode],
					  psp->ddcl[index_mode],
					  psp->ct_size,
					  l,
					  &last_index,
					  cl_tot,
					  psp->ct_size,
					  psp->error_message),
		 psp->error_message,
		 psp->error_message);
    }
    else {
      for (index_ct=0; index_ct<psp->ct_size; index_ct++) 
	cl_tot[index_ct]=0.;
    }
  }
    
  /** B) treat case in which there is only one mode 
         with several initial condition. 
         Fill cl_md_ic[index_mode=0] and sum it to get cl_tot. */

  if ((psp->md_size == 1) && (psp->ic_size[0] > 1)) {
    index_mode = 0;
    for (index_ct=0; index_ct<psp->ct_size; index_ct++) 
      cl_tot[index_ct]=0.;
    for (index_ic1 = 0; index_ic1 < psp->ic_size[index_mode]; index_ic1++) {
      for (index_ic2 = index_ic1; index_ic2 < psp->ic_size[index_mode]; index_ic2++) {
	index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_mode]);
	if (((int)l <= psp->l_max[index_mode]) && 
	    (psp->is_non_zero[index_mode][index_ic1_ic2] == _TRUE_)) {

	  class_call(array_interpolate_spline(psp->l,
					      psp->l_size[index_mode],
					      psp->cl[index_mode],
					      psp->ddcl[index_mode],
					      psp->ic_ic_size[index_mode]*psp->ct_size,
					      l,
					      &last_index,
					      cl_md_ic[index_mode],
					      psp->ic_ic_size[index_mode]*psp->ct_size,
					      psp->error_message),
		     psp->error_message,
		     psp->error_message);
	}
	else {
	  for (index_ct=0; index_ct<psp->ct_size; index_ct++)
	    cl_md_ic[index_mode][index_ic1_ic2*psp->ct_size+index_ct]=0.;
	}

        /* compute cl_tot by summing over cl_md_ic */
	for (index_ct=0; index_ct<psp->ct_size; index_ct++) {
	  if (index_ic1 == index_ic2)
	    cl_tot[index_ct]+=cl_md_ic[index_mode][index_ic1_ic2*psp->ct_size+index_ct];
	  else
	    cl_tot[index_ct]+=2.*cl_md_ic[index_mode][index_ic1_ic2*psp->ct_size+index_ct];
	}
      }
    }
  }

  /** C) loop over modes */

  if (psp->md_size > 1) {

    for (index_ct=0; index_ct<psp->ct_size; index_ct++) 
      cl_tot[index_ct]=0.;

    for (index_mode = 0; index_mode < psp->md_size; index_mode++) {

  /** C.1) treat case in which the mode under consideration 
           has only one initial condition. 
	   Fill cl_md[index_mode]. */

      if (psp->ic_size[index_mode] == 1) {
	if ((int)l <= psp->l_max[index_mode]) {

	  class_call(array_interpolate_spline(psp->l,
					      psp->l_size[index_mode],
					      psp->cl[index_mode],
					      psp->ddcl[index_mode],
					      psp->ct_size,
					      l,
					      &last_index,
					      cl_md[index_mode],
					      psp->ct_size,
					      psp->error_message),
		     psp->error_message,
		     psp->error_message);
	}
	else {
	  for (index_ct=0; index_ct<psp->ct_size; index_ct++) 
	    cl_md[index_mode][index_ct]=0.;
	}
      }

  /** C.2) treat case in which the mode under consideration 
           has several initial conditions. 
	   Fill cl_md_ic[index_mode] and sum it to get cl_md[index_mode] */

      if (psp->ic_size[index_mode] > 1) {
	for (index_ct=0; index_ct<psp->ct_size; index_ct++) 
	  cl_md[index_mode][index_ct]=0.;
	for (index_ic1 = 0; index_ic1 < psp->ic_size[index_mode]; index_ic1++) {
	  for (index_ic2 = index_ic1; index_ic2 < psp->ic_size[index_mode]; index_ic2++) {
	    index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_mode]);
	    if (((int)l <= psp->l_max[index_mode]) && 
		(psp->is_non_zero[index_mode][index_ic1_ic2] == _TRUE_)) {

	      class_call(array_interpolate_spline(psp->l,
						  psp->l_size[index_mode],
						  psp->cl[index_mode],
						  psp->ddcl[index_mode],
						  psp->ic_ic_size[index_mode]*psp->ct_size,
						  l,
						  &last_index,
						  cl_md_ic[index_mode],
						  psp->ic_ic_size[index_mode]*psp->ct_size,
						  psp->error_message),
			 psp->error_message,
			 psp->error_message);
	    }
	    else {
	      for (index_ct=0; index_ct<psp->ct_size; index_ct++)
		cl_md_ic[index_mode][index_ic1_ic2*psp->ct_size+index_ct]=0.;
	    }

	    /* compute cl_md[index_mode] by summing over cl_md_ic[index_mode] */
	    for (index_ct=0; index_ct<psp->ct_size; index_ct++) {
	      if (index_ic1 == index_ic2)
		cl_md[index_mode][index_ct]+=cl_md_ic[index_mode][index_ic1_ic2*psp->ct_size+index_ct];
	      else
		cl_md[index_mode][index_ct]+=2.*cl_md_ic[index_mode][index_ic1_ic2*psp->ct_size+index_ct];
	    }
	  }
	}
      }

  /** C.3) add contribution of cl_md[index_mode] to cl_tot */

      for (index_ct=0; index_ct<psp->ct_size; index_ct++) 
	cl_tot[index_ct]+=cl_md[index_mode][index_ct];
    }
  }

  return _SUCCESS_;

}

/** 
 * Matter power spectrum for arbitrary redshift and for all initial conditions. 
 *
 * This routine evaluates the matter power spectrum at a given value of z by
 * interpolating in the pre-computed table (if several values of z have been stored)
 * or by directly reading it (if it only contains values at z=0 and we want P(k,z=0)) 
 * 
 * 
 * Can be called in two modes: linear or logarithmic.
 * 
 * - linear: returns P(k) (units: Mpc^3)
 * 
 * - logarithmic: returns ln(P(k))
 *
 * One little subtlety: in case of several correlated initial conditions,
 * the cross-correlation spectrum can be negative. Then, in logarithmic mode, 
 * the non-diagonal elements contain the cross-correlation angle P_12/sqrt(P_11 P_22) 
 * (from -1 to 1) instead of ln(P_12)
 *
 * This function can be
 * called from whatever module at whatever time, provided that
 * spectra_init() has been called before, and spectra_free() has not
 * been called yet.
 *
 * @param pba        Input: pointer to background structure (used for converting z into tau)
 * @param psp        Input: pointer to spectra structure (containing pre-computed table)
 * @param mode       Input: linear or logarithmic
 * @param z          Input: redshift
 * @param output_tot Ouput: total matter power spectrum P(k) in Mpc**3 (linear mode), or its logarithms (logarithmic mode)
 * @param output_ic  Ouput: for each pair of initial conditions, matter power spectra P(k) in Mpc**3 (linear mode), or their logarithms and cross-correlation angles (logarithmic mode)
 * @return the error status
 */

int spectra_pk_at_z(
		    struct background * pba,
		    struct spectra * psp,
		    enum linear_or_logarithmic mode, 
		    double z,
		    double * output_tot, /* array with argument output_tot[index_k] (must be already allocated) */
		    double * output_ic   /* array with argument output_tot[index_k * psp->ic_ic_size[index_mode] + index_ic1_ic2] (must be already allocated only if more than one initial condition) */
		    ) {

  /** Summary: */

  /** - define local variables */

  int index_mode;
  int last_index;
  int index_k;
  double tau,ln_tau;
  int index_ic1,index_ic2,index_ic1_ic2;

  index_mode = psp->index_md_scalars;

  /** - first step: convert z into ln(tau) */

  class_call(background_tau_of_z(pba,z,&tau),
	     pba->error_message,
	     psp->error_message);

  class_test(tau <= 0.,
	     psp->error_message,
	     "negative or null value of conformal time: cannot interpolate");

  ln_tau = log(tau);

  /** - second step: for both modes (linear or logarithmic), store the spectrum in logarithmic format in the output array(s) */

  /**   (a.) if only values at tau=tau_today are stored and we want P(k,z=0), no need to interpolate */

  if (psp->ln_tau_size == 1) {

    class_test(z != 0.,
	       psp->error_message,
	       "asked z=%e but only P(k,z=0) has been tabulated",z);

    for (index_k=0; index_k<psp->ln_k_size; index_k++)
      if (psp->ic_size[index_mode] == 1) {
      	output_tot[index_k] = psp->ln_pk[index_k];
      }
      else {
	for (index_ic1_ic2 = 0; index_ic1_ic2 < psp->ic_ic_size[index_mode]; index_ic1_ic2++) {
	  output_ic[index_k * psp->ic_ic_size[index_mode] + index_ic1_ic2] = 
	    psp->ln_pk[index_k * psp->ic_ic_size[index_mode] + index_ic1_ic2];
	}
      }
  }

  /**   (b.) if several values of tau have been stored, use interpolation routine to get spectra at correct redshift */

  else {

    if (psp->ic_ic_size[index_mode] == 1) {
    
      class_call(array_interpolate_spline(psp->ln_tau,
					  psp->ln_tau_size,
					  psp->ln_pk,
					  psp->ddln_pk,
					  psp->ln_k_size,
					  ln_tau,
					  &last_index,
					  output_tot,
					  psp->ln_k_size,
					  psp->error_message),
		 psp->error_message,
		 psp->error_message);

    }
    else {
      
      class_call(array_interpolate_spline(psp->ln_tau,
					  psp->ln_tau_size,
					  psp->ln_pk,
					  psp->ddln_pk,
					  psp->ic_ic_size[index_mode]*psp->ln_k_size,
					  ln_tau,
					  &last_index,
					  output_ic,
					  psp->ic_ic_size[index_mode]*psp->ln_k_size,
					  psp->error_message),
		 psp->error_message,
		 psp->error_message);
    }
  }

  /** - third step: if there are several initial conditions, compute the total P(k) and set back all uncorrelated coefficients to exactly zero. Check positivity of total P(k). */

  if (psp->ic_size[index_mode] > 1) {
    for (index_k=0; index_k<psp->ln_k_size; index_k++) {
      output_tot[index_k] = 0.;
      for (index_ic1=0; index_ic1 < psp->ic_size[index_mode]; index_ic1++) {
	for (index_ic2 = index_ic1; index_ic2 < psp->ic_size[index_mode]; index_ic2++) {
	  index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_mode]);
	  if (index_ic1 == index_ic2) {
	    output_tot[index_k] += exp(output_ic[index_k * psp->ic_ic_size[index_mode] + index_ic1_ic2]);
	  }
	  else {
	    if (psp->is_non_zero[index_mode][index_ic1_ic2] == _TRUE_) {
	      output_tot[index_k] += 
		2. * output_ic[index_k * psp->ic_ic_size[index_mode] + index_ic1_ic2] *
		sqrt(exp(output_ic[index_k * psp->ic_ic_size[index_mode] + index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_mode])]) *
		     exp(output_ic[index_k * psp->ic_ic_size[index_mode] + index_symmetric_matrix(index_ic2,index_ic2,psp->ic_size[index_mode])]));
	}
	    else
	      output_ic[index_k * psp->ic_ic_size[index_mode] + index_ic1_ic2] = 0.;
	  }
	}
      }

      class_test(output_tot[index_k] <= 0.,
		 psp->error_message,
		 "for k=%e, z=%e, the matrix of initial condition amplitudes was not positive definite, hence P(k)_total=%e results negative",
		 exp(psp->ln_k[index_k]),z,output_tot[index_k]);

    }
  }

  /** - fourth step: depending on requested mode (linear or logarithmic), apply necessary transformation to the output arrays */

  /**   (a.) linear mode: if only one initial condition, convert output_pk to linear format; if several initial conditions, convert output_ic to linear format, output_tot is already in this format */

  if (mode == linear) {

    if (psp->ic_size[index_mode] == 1) {
      for (index_k=0; index_k<psp->ln_k_size; index_k++) {
	output_tot[index_k] = exp(output_tot[index_k]); 
      }
    }

    else {
      for (index_k=0; index_k<psp->ln_k_size; index_k++) {
	for (index_ic1=0; index_ic1 < psp->ic_size[index_mode]; index_ic1++) {
	  index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_mode]);
	  output_ic[index_k * psp->ic_ic_size[index_mode] + index_ic1_ic2] = exp(output_ic[index_k * psp->ic_ic_size[index_mode] + index_ic1_ic2]);
	}
	for (index_ic1=0; index_ic1 < psp->ic_size[index_mode]; index_ic1++) {
	  for (index_ic2 = index_ic1+1; index_ic2 < psp->ic_size[index_mode]; index_ic2++) {

	    output_ic[index_k * psp->ic_ic_size[index_mode] + index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_mode])] =
	      output_ic[index_k * psp->ic_ic_size[index_mode] + index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_mode])]
	      *sqrt(output_ic[index_k * psp->ic_ic_size[index_mode] + index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_mode])] *
		    output_ic[index_k * psp->ic_ic_size[index_mode] + index_symmetric_matrix(index_ic2,index_ic2,psp->ic_size[index_mode])]);
	  }
	}
      }
    }
  }

  /**   (b.) logarithmic mode: if only one initial condition, nothing to be done; if several initial conditions, convert output_tot to logarithmic format, output_ic is already in this format */

  else {

    if (psp->ic_size[index_mode] > 1) {
      for (index_k=0; index_k<psp->ln_k_size; index_k++) {
	/* we have already checked above that output_tot was positive */
	output_tot[index_k] = log(output_tot[index_k]); 
      }
    }
  }
  
  return _SUCCESS_;

}

/** 
 * Matter power spectrum for arbitrary wavenumber, redshift and initial condition. 
 *
 * This routine evaluates the matter power spectrum at a given value of k and z by
 * interpolating in a table of all P(k)'s computed at this z by spectra_pk_at_z() (when kmin <= k <= kmax), 
 * or eventually by using directly the primordial spectrum (when 0 <= k < kmin): 
 * the latter case is an approximation, valid when kmin << comoving Hubble scale today.
 * Returns zero when k=0. Returns an error when k<0 or k > kmax.
 * 
 * This function can be
 * called from whatever module at whatever time, provided that
 * spectra_init() has been called before, and spectra_free() has not
 * been called yet.
 *
 * @param pba        Input: pointer to background structure (used for converting z into tau)
 * @param ppm        Input: pointer to primordial structure (used only in the case 0 < k < kmin)
 * @param psp        Input: pointer to spectra structure (containing pre-computed table)
 * @param k          Input: wavenumber in 1/Mpc
 * @param z          Input: redshift
 * @param pk_tot     Ouput: total matter power spectrum P(k) in Mpc**3
 * @param pk_ic      Ouput: for each pair of initial conditions, matter power spectra P(k) in Mpc**3
 * @return the error status
 */

int spectra_pk_at_k_and_z(
			  struct background * pba,
			  struct primordial * ppm,
			  struct spectra * psp,
			  double k,
			  double z,
			  double * pk_tot, /* pointer to a single number (must be already allocated) */
			  double * pk_ic   /* array of argument pk_ic[index_ic1_ic2] (must be already allocated only if several initial conditions) */
			  ) {

  /** Summary: */

  /** - define local variables */

  int index_mode;
  int index_k;
  int last_index;
  int index_ic1,index_ic2,index_ic1_ic2;

  double * spectrum_at_z = NULL;
  double * spectrum_at_z_ic = NULL;
  double * spline;
  double * pk_primordial_k = NULL;
  double kmin;
  double * pk_primordial_kmin = NULL;

  index_mode = psp->index_md_scalars;

  /** - first step: check that k is in valid range [0:kmax] (the test for z will be done when calling spectra_pk_at_z()) */
 
  class_test((k < 0.) || (k > exp(psp->ln_k[psp->ln_k_size-1])),
	     psp->error_message,
	     "k=%e out of bounds [%e:%e]",k,0.,exp(psp->ln_k[psp->ln_k_size-1]));
  
  /** - deal with case 0 <= k < kmin */

  if (k < exp(psp->ln_k[0])) {
    
    /**   (a.) subcase k=0: then P(k)=0 */
  
    if (k == 0.) {
      if (psp->ic_size[index_mode] == 1) {
	*pk_tot=0.;
      }
      else {
	for (index_ic1_ic2 = 0; index_ic1_ic2 < psp->ic_ic_size[index_mode]; index_ic1_ic2++) {
	  pk_ic[index_ic1_ic2] = 0.;
	}
      }
    }

    /**    (b.) subcase 0<k<kmin: in this case we know that on super-Hubble scales:        
     *          P(k) = [some number] * k  * P_primordial(k) 
     *          so     
     *          P(k) = P(kmin) * (k P_primordial(k)) / (kmin P_primordial(kmin)) 
     *          (note that the result is accurate only if kmin is such that [a0 kmin] << H0) 
     */

    else {

      /* compute P(k,z) which contains P(kmin,z)*/
      class_alloc(spectrum_at_z,
		  psp->ln_k_size*sizeof(double),
		  psp->error_message);
      if (psp->ic_size[index_mode] > 1) {
	class_alloc(spectrum_at_z_ic,
		    sizeof(double)*psp->ic_ic_size[index_mode]*psp->ln_k_size,
		    psp->error_message); 
      }
      class_call(spectra_pk_at_z(pba,
				 psp,
				 linear,
				 z,
				 spectrum_at_z,
				 spectrum_at_z_ic),
		 psp->error_message,
		 psp->error_message);
      
      /* compute P_primordial(k) */
      class_alloc(pk_primordial_k,
		  sizeof(double)*psp->ic_ic_size[index_mode],
		  psp->error_message);
      class_call(primordial_spectrum_at_k(ppm,
					  index_mode,
					  linear,
					  k,
					  pk_primordial_k),
		 ppm->error_message,psp->error_message);

      /* compute P_primordial(kmin) */
      kmin = exp(psp->ln_k[0]);
      class_alloc(pk_primordial_kmin,
		  sizeof(double)*psp->ic_ic_size[index_mode],
		  psp->error_message);
      class_call(primordial_spectrum_at_k(ppm,
					  index_mode,
					  linear,
					  kmin,
					  pk_primordial_kmin),
		 ppm->error_message,
		 psp->error_message);
    
      /* apply above analytic approximation for P(k) */
      index_k=0;
      if (psp->ic_size[index_mode] == 1) {
	index_ic1_ic2 = 0;
	*pk_tot = spectrum_at_z[index_k]
	  *k*pk_primordial_k[index_ic1_ic2]
	  /kmin/pk_primordial_kmin[index_ic1_ic2];
      }
      else {
      	for (index_ic1_ic2 = 0; index_ic1_ic2 < psp->ic_ic_size[index_mode]; index_ic1_ic2++) {
	  pk_ic[index_ic1_ic2] = spectrum_at_z_ic[index_ic1_ic2]
	    *k*pk_primordial_k[index_ic1_ic2]
	    /kmin/pk_primordial_kmin[index_ic1_ic2];
	}
      }
    }

    free(spectrum_at_z);
    if (psp->ic_size[index_mode] > 1)
      free(spectrum_at_z_ic);
    free(pk_primordial_k);
    free(pk_primordial_kmin);

  }

  /** - deal with case kmin <= k <= kmax */

  else {

    /* compute P(k,z) (in logarithmic format for more accurate interpolation) */
    class_alloc(spectrum_at_z,
		psp->ln_k_size*sizeof(double),
		psp->error_message);
    if (psp->ic_size[index_mode] > 1) {
      class_alloc(spectrum_at_z_ic,
		  sizeof(double)*psp->ic_ic_size[index_mode]*psp->ln_k_size,
		  psp->error_message); 
    }
    class_call(spectra_pk_at_z(pba,
			       psp,
			       logarithmic,
			       z,
			       spectrum_at_z,
			       spectrum_at_z_ic),
	       psp->error_message,
	       psp->error_message);

    /* get its second derivatives with spline, then interpolate, then convert to linear format */

    class_alloc(spline,
		sizeof(double)*psp->ic_ic_size[index_mode]*psp->ln_k_size,
		psp->error_message);

    if (psp->ic_size[index_mode] == 1) {

      class_call(array_spline_table_lines(psp->ln_k,
					  psp->ln_k_size,
					  spectrum_at_z,
					  1,
					  spline,
					  _SPLINE_NATURAL_,
					  psp->error_message),
		 psp->error_message,
		 psp->error_message);
      
      class_call(array_interpolate_spline(psp->ln_k,
					  psp->ln_k_size,
					  spectrum_at_z,
					  spline,
					  1,
					  log(k),
					  &last_index,
					  pk_tot,
					  1,
					  psp->error_message),
		 psp->error_message,
		 psp->error_message);

      *pk_tot = exp(*pk_tot);

    }
    else {
      
      class_call(array_spline_table_lines(psp->ln_k,
					  psp->ln_k_size,
					  spectrum_at_z_ic,
					  psp->ic_ic_size[index_mode],
					  spline,
					  _SPLINE_NATURAL_,
					  psp->error_message),
		 psp->error_message,
		 psp->error_message);

      class_call(array_interpolate_spline(psp->ln_k,
					  psp->ln_k_size,
					  spectrum_at_z_ic,
					  spline,
					  psp->ic_ic_size[index_mode],
					  log(k),
					  &last_index,
					  pk_ic,
					  psp->ic_ic_size[index_mode],
					  psp->error_message),
		 psp->error_message,
		 psp->error_message);

      for (index_ic1 = 0; index_ic1 < psp->ic_size[index_mode]; index_ic1++) {
	index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_mode]);
	pk_ic[index_ic1_ic2] = exp(pk_ic[index_ic1_ic2]);
      }
      for (index_ic1 = 0; index_ic1 < psp->ic_size[index_mode]; index_ic1++) {
	for (index_ic2 = index_ic1+1; index_ic2 < psp->ic_size[index_mode]; index_ic2++) {
	  index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_mode]);
	  if (psp->is_non_zero[index_mode][index_ic1_ic2] == _TRUE_) {
	    pk_ic[index_ic1_ic2] = pk_ic[index_ic1_ic2]*
	      sqrt(pk_ic[index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_mode])]*
		   pk_ic[index_symmetric_matrix(index_ic2,index_ic2,psp->ic_size[index_mode])]);
	  }
	  else {
	    pk_ic[index_ic1_ic2] = 0.;
	  }
	}
      }
      free(spectrum_at_z_ic);
    }

    free(spectrum_at_z);
    free(spline);
  }

  /** - last step: if more than one condition, sum over pk_ic to get pk_tot, and set back coefficients of non-correlated pairs to exactly zero. */

  if (psp->ic_size[index_mode] > 1) {

    *pk_tot = 0.;
      
    for (index_ic1 = 0; index_ic1 < psp->ic_size[index_mode]; index_ic1++) {
      for (index_ic2 = index_ic1; index_ic2 < psp->ic_size[index_mode]; index_ic2++) {
	index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_mode]);
	
	if (psp->is_non_zero[index_mode][index_ic1_ic2] == _TRUE_) {
	  
	  if (index_ic1 == index_ic2)
	    *pk_tot += pk_ic[index_ic1_ic2];
	  else
	    *pk_tot += 2.*pk_ic[index_ic1_ic2];
	}
	else {
	  pk_ic[index_ic1_ic2] = 0.;
	}
      }
    }

    class_test(*pk_tot <= 0.,
	       psp->error_message,
	       "for k=%e, the matrix of initial condition amplitudes was not positive definite, hence P(k)_total results negative",k);
    
  }

  return _SUCCESS_;

}

/** 
 * Matter transfer functions T_i(k) for arbitrary redshift and for all
 * initial conditions.
 *
 * This routine evaluates the matter transfer functions at a given value of z by
 * interpolating in the pre-computed table (if several values of z have been stored)
 * or by directly reading it (if it only contains values at z=0 and we want T_i(k,z=0)) 
 * 
 * 
 * This function can be
 * called from whatever module at whatever time, provided that
 * spectra_init() has been called before, and spectra_free() has not
 * been called yet.
 *
 * @param pba        Input: pointer to background structure (used for converting z into tau)
 * @param psp        Input: pointer to spectra structure (containing pre-computed table)
 * @param z          Input: redshift
 * @param output     Ouput: matter transfer functions
 * @return the error status
 */

int spectra_tk_at_z(
			   struct background * pba,
			   struct spectra * psp,
			   double z,
			   double * output /* array with argument output[(index_k*psp->ic_size[index_mode]+index_ic)*psp->tr_size+index_tr] (must be already allocated) */
			   ) {

  /** Summary: */

  /** - define local variables */

  int index_mode;
  int last_index;
  int index_k;
  int index_tr;
  double tau,ln_tau;
  int index_ic;

  index_mode = psp->index_md_scalars;

  /** - first step: convert z into ln(tau) */

  class_call(background_tau_of_z(pba,z,&tau),
	     pba->error_message,
	     psp->error_message);

  class_test(tau <= 0.,
	     psp->error_message,
	     "negative or null value of conformal time: cannot interpolate");

  ln_tau = log(tau);

  /** - second step: store the matter transfer functions in the output array */

  /**   (a.) if only values at tau=tau_today are stored and we want T_i(k,z=0), no need to interpolate */

  if (psp->ln_tau_size == 1) {
    
    class_test(z != 0.,
	       psp->error_message,
	       "asked z=%e but only T_i(k,z=0) has been tabulated",z);
    
    for (index_k=0; index_k<psp->ln_k_size; index_k++)
      for (index_tr=0; index_tr<psp->tr_size; index_tr++)
	for (index_ic = 0; index_ic < psp->ic_size[index_mode]; index_ic++)
	  output[(index_k*psp->ic_size[index_mode]+index_ic)*psp->tr_size+index_tr]
	    = psp->matter_transfer[(index_k*psp->ic_size[index_mode]+index_ic)*psp->tr_size+index_tr];
    
  }

  /**   (b.) if several values of tau have been stored, use interpolation routine to get spectra at correct redshift */

  else {
      
    class_call(array_interpolate_spline(psp->ln_tau,
					psp->ln_tau_size,
					psp->matter_transfer,
					psp->ddmatter_transfer,
					psp->ic_size[index_mode]*psp->tr_size*psp->ln_k_size,
					ln_tau,
					&last_index,
					output,
					psp->ic_size[index_mode]*psp->tr_size*psp->ln_k_size,
					psp->error_message),
	       psp->error_message,
	       psp->error_message);
    
  }

  return _SUCCESS_;

}

/** 
 * Matter transfer functions T_i(k) for arbitrary wavenumber, redshift
 * and initial condition.
 *
 * This routine evaluates the matter transfer functions at a given
 * value of k and z by interpolating in a table of all T_i(k,z)'s
 * computed at this z by spectra_tk_at_z() (when kmin <= k <= kmax).
 * Returns an error when k<kmin or k > kmax.
 * 
 * This function can be called from whatever module at whatever time,
 * provided that spectra_init() has been called before, and
 * spectra_free() has not been called yet.
 *
 * @param pba        Input: pointer to background structure (used for converting z into tau)
 * @param psp        Input: pointer to spectra structure (containing pre-computed table)
 * @param k          Input: wavenumber in 1/Mpc
 * @param z          Input: redshift
 * @param output     Ouput: matter transfer functions
 * @return the error status
 */

int spectra_tk_at_k_and_z(
			  struct background * pba,
			  struct spectra * psp,
			  double k,
			  double z,
			  double * output  /* array with argument output[index_ic*psp->tr_size+index_tr] (must be already allocated) */
			  ) {

  /** Summary: */

  /** - define local variables */

  int index_mode;
  int last_index;
  double * tks_at_z;
  double * ddtks_at_z;

  index_mode = psp->index_md_scalars;

  /** - first step: check that k is in valid range [0:kmax] (the test for z will be done when calling spectra_tk_at_z()) */
 
  class_test((k < 0.) || (k > exp(psp->ln_k[psp->ln_k_size-1])),
	     psp->error_message,
	     "k=%e out of bounds [%e:%e]",k,0.,exp(psp->ln_k[psp->ln_k_size-1]));
  
  /* compute T_i(k,z) */

  class_alloc(tks_at_z,
	      psp->ln_k_size*psp->tr_size*psp->ic_size[index_mode]*sizeof(double),
	      psp->error_message);

  class_call(spectra_tk_at_z(pba,
			     psp,
			     z,
			     tks_at_z),
	     psp->error_message,
	     psp->error_message);

  /* get its second derivatives w.r.t. k with spline, then interpolate */

  class_alloc(ddtks_at_z,
	      psp->ln_k_size*psp->tr_size*psp->ic_size[index_mode]*sizeof(double),
	      psp->error_message);

  class_call(array_spline_table_lines(psp->ln_k,
				      psp->ln_k_size,
				      tks_at_z,
				      psp->tr_size*psp->ic_size[index_mode],
				      ddtks_at_z,
				      _SPLINE_NATURAL_,
				      psp->error_message),
	     psp->error_message,
	     psp->error_message);
      
  class_call(array_interpolate_spline(psp->ln_k,
				      psp->ln_k_size,
				      tks_at_z,
				      ddtks_at_z,
				      psp->tr_size*psp->ic_size[index_mode],
				      log(k),
				      &last_index,
				      output,
				      psp->tr_size*psp->ic_size[index_mode],
				      psp->error_message),
	     psp->error_message,
	     psp->error_message);

  free(tks_at_z);
  free(ddtks_at_z);

  return _SUCCESS_;

}

/**
 * This routine initializes the spectra structure (in particular, 
 * computes table of anisotropy and Fourier spectra \f$ C_l^{X}, P(k), ... \f$)
 * 
 * @param pba Input : pointer to background structure (will provide H, Omega_m at redshift of interest)
 * @param ppt Input : pointer to perturbation structure
 * @param ptr Input : pointer to transfer structure
 * @param ppm Input : pointer to primordial structure
 * @param psp Output: pointer to initialized spectra structure
 * @return the error status
 */

int spectra_init(
		 struct background * pba,
		 struct perturbs * ppt,
		 struct transfers * ptr,
		 struct primordial * ppm,
		 struct spectra * psp
		 ) {

  /** Summary: */

  /** - check that we really want to compute at least one spectrum */

  if ((ppt->has_cls == _FALSE_) && (ppt->has_pk_matter == _FALSE_)) {
    psp->md_size = 0;
    if (psp->spectra_verbose > 0)
      printf("No spectra requested. Spectra module skipped.\n");
    return _SUCCESS_;
  }
  else {
    if (psp->spectra_verbose > 0)
      printf("Computing unlensed linear spectra\n");
  }

  /** - initialize indices and allocate some of the arrays in the 
        spectra structure */

  class_call(spectra_indices(pba,ppt,ptr,ppm,psp),
	     psp->error_message,
	     psp->error_message);

  /** - deal with C_l's, if any */

  if (ppt->has_cls == _TRUE_) {

    class_call(spectra_cls(ppt,ptr,ppm,psp),
	       psp->error_message,
	       psp->error_message);

  }
  else {
    psp->ct_size=0;
  }

  /** - deal with P(k,tau) and T_i(k,tau) */

  if ((ppt->has_pk_matter == _TRUE_) || (ppt->has_pk_matter == _TRUE_)) {

    class_call(spectra_k_and_tau(pba,ppt,psp),
	       psp->error_message,
	       psp->error_message);

    if (ppt->has_pk_matter == _TRUE_) {

      class_call(spectra_pk(pba,ppt,ppm,psp),
		 psp->error_message,
		 psp->error_message);
    }
    else {
      psp->ln_pk=NULL;
    }

    if (ppt->has_matter_transfers == _TRUE_) {

      class_call(spectra_matter_transfers(pba,ppt,psp),
		 psp->error_message,
		 psp->error_message);
    }
    else {
      psp->matter_transfer=NULL;
    }

  }
  else {
    psp->ln_k_size=0;
  }

  return _SUCCESS_;
}

/**
 * This routine frees all the memory space allocated by spectra_init().
 *
 * To be called at the end of each run, only when no further calls to
 * spectra_cls_at_l(), spectra_pk_at_z(), spectra_pk_at_k_and_z() are needed.
 *
 * @param psp Input: pointer to spectra structure (which fields must be freed)
 * @return the error status
 */

int spectra_free(
		 struct spectra * psp
		 ) {

  int index_mode;

  if (psp->md_size > 0) {

    if (psp->ct_size > 0) {
    
      for (index_mode = 0; index_mode < psp->md_size; index_mode++) {
	free(psp->cl[index_mode]);
	free(psp->ddcl[index_mode]);
      }
      free(psp->l);
      free(psp->l_size);
      free(psp->l_max);
      free(psp->cl);
      free(psp->ddcl);
    }

    if (psp->ln_k_size > 0) {

      free(psp->ln_tau);
      free(psp->ln_k);

      if (psp->ln_pk != NULL) { 

	free(psp->ln_pk);

	if (psp->ln_tau_size > 1) {
	  free(psp->ddln_pk);
	}
      }

      if (psp->matter_transfer != NULL) {
	
	free(psp->matter_transfer);
	if (psp->ln_tau_size > 1) {
	  free(psp->ddmatter_transfer);
	}
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

/**
 * This routine defines indices and allocates tables in the spectra structure 
 *
 * @param ppt  Input : pointer to perturbation structure
 * @param ptr  Input : pointer to transfers structure
 * @param ppm  Input : pointer to primordial structure
 * @param psp  Input/output: pointer to spectra structure 
 * @return the error status
 */

int spectra_indices(
		    struct background * pba,
		    struct perturbs * ppt,
		    struct transfers * ptr,
		    struct primordial * ppm,
		    struct spectra * psp
		    ){

  int index_ct;
  int index_mode;
  int index_ic1_ic2;
  int index_tr;

  psp->md_size = ppm->md_size;

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
    psp->ic_size[index_mode] = ppm->ic_size[index_mode];
    psp->ic_ic_size[index_mode] = ppm->ic_ic_size[index_mode];
    class_alloc(psp->is_non_zero[index_mode],
		sizeof(short)*psp->ic_ic_size[index_mode],
		psp->error_message);
    for (index_ic1_ic2=0; index_ic1_ic2 < psp->ic_ic_size[index_mode]; index_ic1_ic2++)
      psp->is_non_zero[index_mode][index_ic1_ic2] = ppm->is_non_zero[index_mode][index_ic1_ic2];
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

    if (ppt->has_cl_cmb_polarization == _TRUE_) {
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

  if (ppt->has_matter_transfers == _TRUE_) {

    /* indices for species associated with a matter transfer function in Fourier space */

    index_tr=0;
    psp->index_tr_g = index_tr;
    index_tr++;
    psp->index_tr_b = index_tr;
    index_tr++;
    if (pba->has_cdm == _TRUE_) {
      psp->index_tr_cdm = index_tr;
      index_tr++;
    }
    if (pba->has_fld == _TRUE_) {
      psp->index_tr_de = index_tr;
      index_tr++;
    }
    if (pba->has_ur == _TRUE_) {
      psp->index_tr_ur = index_tr;
      index_tr++;
    }
    if (pba->has_ncdm == _TRUE_) {
      psp->index_tr_ncdm1 = index_tr;
      index_tr+=pba->N_ncdm;
    }
    psp->index_tr_tot = index_tr;
    index_tr++;
    psp->tr_size = index_tr;

  }

  return _SUCCESS_;

}

/**
 * This routine computes a table of values for all harmonic spectra C_l's,
 * given the transfer functions and primordial spectra.
 * 
 * @param ppt Input : pointer to perturbation structure
 * @param ptr Input : pointer to transfers structure
 * @param ppm Input : pointer to primordial structure
 * @param psp Input/Output: pointer to spectra structure 
 * @return the error status
 */

int spectra_cls(
		struct perturbs * ppt,
		struct transfers * ptr,
		struct primordial * ppm,
		struct spectra * psp
		) {

  /** Summary: */

  /** - define local variables */

  int index_mode;
  int index_ic1,index_ic2,index_ic1_ic2;
  int index_l;
  int index_ct;
  int cl_integrand_num_columns;

  double * cl_integrand; /* array with argument cl_integrand[index_k*cl_integrand_num_columns+1+psp->index_ct] */
  double * transfer_ic1; /* array with argument transfer_ic1[index_tt] */
  double * transfer_ic2; /* idem */
  double * primordial_pk;  /* array with argument primordial_pk[index_ic_ic]*/

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

  /** - allocate pointers to arrays where results will be stored */

  class_alloc(psp->l_size,sizeof(int)*psp->md_size,psp->error_message);
  class_alloc(psp->l_max,sizeof(int)*psp->md_size,psp->error_message);
  class_alloc(psp->cl,sizeof(double *)*psp->md_size,psp->error_message);
  class_alloc(psp->ddcl,sizeof(double *)*psp->md_size,psp->error_message);

  psp->l_size_max = ptr->l_size_max;
  class_alloc(psp->l,sizeof(double)*psp->l_size_max,psp->error_message);

  /** - store values of l */
  for (index_l=0; index_l < psp->l_size_max; index_l++) {
    psp->l[index_l] = (double)ptr->l[index_l];
  }

  psp->l_max_tot = 0;

  /** - loop over modes (scalar, tensors, etc). For each mode: */

  for (index_mode = 0; index_mode < psp->md_size; index_mode++) {

    /** - a) store number of l values for this mode */

    psp->l_size[index_mode] = ptr->l_size[index_mode];

    /** - b) allocate arrays where results will be stored */

    class_alloc(psp->cl[index_mode],sizeof(double)*psp->l_size[index_mode]*psp->ct_size*psp->ic_ic_size[index_mode],psp->error_message);
    class_alloc(psp->ddcl[index_mode],sizeof(double)*psp->l_size[index_mode]*psp->ct_size*psp->ic_ic_size[index_mode],psp->error_message);
    cl_integrand_num_columns = 1+psp->ct_size*2; /* one for k, ct_size for each type, ct_size for each second derivative of each type */

    /** c) get from input the last multipole for each mode (given as an input) 
	at which we trust our C_ls; this number 
	l[l_size[index_mode]-1] can be larger than 
	l_max[index_mode], 
	in order to ensure better interpolation with no boundary effects.
	Compute also the max over all modes */

    if ((ppt->has_scalars) && (index_mode == ppt->index_md_scalars)) {
      psp->l_max[index_mode] = ppt->l_scalar_max;
    }
    
    if ((ppt->has_tensors) && (index_mode == ppt->index_md_tensors)) {
      psp->l_max[index_mode] = ppt->l_tensor_max;
    }

    psp->l_max_tot=max(psp->l_max_tot,psp->l_max[index_mode]);

    /** d) loop over initial conditions */

    for (index_ic1 = 0; index_ic1 < psp->ic_size[index_mode]; index_ic1++) {
      for (index_ic2 = index_ic1; index_ic2 < psp->ic_size[index_mode]; index_ic2++) {
	index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_mode]);

	/* non-diagonal coefficients should be computed only if non-zero correlation */
	if (psp->is_non_zero[index_mode][index_ic1_ic2] == _TRUE_) {

	  /* initialize error management flag */
	  abort = _FALSE_;

	  /* beginning of parallel region */

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

	    /** - loop over l values defined in the transfer module. 
		For each l, compute the C_l's for all types (TT, TE, ...) 
		by convolving primordial spectra with transfer  functions. 
		This elementary task is assigned to spectra_compute_cl() */

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
	
	  for (index_l=0; index_l < ptr->l_size[index_mode]; index_l++) {
	    for (index_ct=0; index_ct<psp->ct_size; index_ct++) {
	      psp->cl[index_mode]
		[(index_l * psp->ic_ic_size[index_mode] + index_ic1_ic2) * psp->ct_size + index_ct]
		= 0.;
	    }
	  }
	}
      }
    }

    /** - e) now that for a given mode, all possible C_l's have been computed, 
	compute second derivative of the array in which they are stored, 
	in view of spline interpolation. */

    class_call(array_spline_table_lines(psp->l,
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

/**
 * This routine computes the C_l's for a given mode, pair of initial conditions
 * and multipole, but for all types (TT, TE...), by convolving the
 * transfer functions with the primordial spectra.
 * 
 * @param ppt           Input : pointer to perturbation structure
 * @param ptr           Input : pointer to transfers structure
 * @param ppm           Input : pointer to primordial structure
 * @param psp           Input/Output: pointer to spectra structure (result stored here)
 * @param index_mode    Input : index of mode under consideration
 * @param index_ic1     Input : index of first initial condition in the correlator
 * @param index_ic2     Input : index of second initial condition in the correlato
 * @param index_l       Input : index of multipole under consideration
 * @param cl_integrand_num_column Input : number of columns in cl_integrand 
 * @param cl_integrand  Input : an allocated workspace
 * @param primordial_pk Input : table of primordial spectrum values
 * @param transfer_ic1  Input : table of transfer function values for first initial condition
 * @param transfer_ic2  Input : table of transfer function values for second initial condition
 * @return the error status
 */

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
  int index_ic1_ic2;

  index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_mode]);

  for (index_k=0; index_k < ptr->k_size[index_mode]; index_k++) {

    k = ptr->k[index_mode][index_k];

    cl_integrand[index_k*cl_integrand_num_columns+0] = k;

    class_call(primordial_spectrum_at_k(ppm,index_mode,linear,k,primordial_pk),
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

    if (((psp->has_bb == _TRUE_) && (index_ct == psp->index_ct_bb) && (ppt->has_scalars == _TRUE_) && (index_mode == ppt->index_md_scalars)) ||
	((psp->has_pp == _TRUE_) && (index_ct == psp->index_ct_pp) && (ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors)) ||
	((psp->has_tp == _TRUE_) && (index_ct == psp->index_ct_tp) && (ppt->has_tensors == _TRUE_) && (index_mode == ppt->index_md_tensors))) {

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

    }
  }

  return _SUCCESS_;

}

/**
 * This routine computes the values of k and tau at which the matter
 * power spectra P(k,tau) and the matter transfer functions T_i(k,tau)
 * will be stored.
 * 
 * @param pba Input : pointer to background structure (for z to tau conversion)
 * @param ppt Input : pointer to perturbation structure (contain source functions)
 * @param psp Input/Output: pointer to spectra structure 
 * @return the error status
 */

int spectra_k_and_tau(
		      struct background * pba,
		      struct perturbs * ppt,
		      struct spectra * psp
		      ) {

  /** Summary: */

  /** - define local variables */

  int index_mode;
  int index_k;
  int index_tau;
  double tau_min;

  /** - check the presence of scalar modes */

  class_test((ppt->has_scalars == _FALSE_),
	     psp->error_message,
	     "you cannot ask for matter power spectrum since you turned off scalar modes");
  
  psp->index_md_scalars = ppt->index_md_scalars;
  index_mode = psp->index_md_scalars;

  /** - check the maximum redshift z_max_pk at which P(k,z) and T_i(k,z) should be 
      computable by interpolation. If it is equal to zero, only P(k,z=0) 
      needs to be computed. If it is higher, we will store in a table 
      various P(k,tau) at several values of tau generously encompassing 
      the range 0<z<z_max_pk */

  /* if z_max_pk<0, return error */
  class_test((psp->z_max_pk < 0),
	     psp->error_message,
	     "asked for negative redshift z=%e",psp->z_max_pk);

  /* if z_max_pk=0, there is just one value to store */
  if (psp->z_max_pk == 0.) {
    psp->ln_tau_size=1;
  }

  /* if z_max_pk>0, store several values (with a confortable margin above z_max_pk) in view of interpolation */
  else{

    /* find the first relevant value of tau (last value in the table tau_ampling before tau(z_max)) and infer the number of values of tau at which P(k) must be stored */

    class_call(background_tau_of_z(pba,psp->z_max_pk,&tau_min),
	       pba->error_message,
	       psp->error_message);

    index_tau=0;
    class_test((tau_min < ppt->tau_sampling[index_tau]),
	       psp->error_message,
	       "you asked for zmax=%e, i.e. taumin=%e, smaller than first possible value =%e",psp->z_max_pk,tau_min,ppt->tau_sampling[0]);
  
    while (ppt->tau_sampling[index_tau] < tau_min){
      index_tau++;
    }
    index_tau --;
    /* whenever possible, take a few more values in to avoid boundary effects in the interpolation */
    if (index_tau>0) index_tau--; 
    if (index_tau>0) index_tau--; 
    if (index_tau>0) index_tau--; 
    if (index_tau>0) index_tau--; 
    psp->ln_tau_size=ppt->tau_size-index_tau;

  }

  /** - allocate and fill table of tau values at which P(k,tau) and T_i(k,tau) are stored */

  class_alloc(psp->ln_tau,sizeof(double)*psp->ln_tau_size,psp->error_message);

  for (index_tau=0; index_tau<psp->ln_tau_size; index_tau++) {
    psp->ln_tau[index_tau]=log(ppt->tau_sampling[index_tau-psp->ln_tau_size+ppt->tau_size]);
  }

  /** - allocate and fill table of k values at which P(k,tau) is stored */

  psp->ln_k_size = ppt->k_size[index_mode];
  class_alloc(psp->ln_k,sizeof(double)*psp->ln_k_size,psp->error_message);

  for (index_k=0; index_k<psp->ln_k_size; index_k++) {
    class_test(ppt->k[index_mode][index_k] <= 0.,
	       psp->error_message,
	       "stop to avoid segmentation fault");
    psp->ln_k[index_k]=log(ppt->k[index_mode][index_k]);
  }

  return _SUCCESS_;
}

/**
 * This routine computes a table of values for all matter power spectra P(k),
 * given the source functions and primordial spectra.
 * 
 * @param pba Input : pointer to background structure (will provide H, Omega_m at redshift of interest)
 * @param ppt Input : pointer to perturbation structure (contain source functions)
 * @param ppm Input : pointer to primordial structure
 * @param psp Input/Output: pointer to spectra structure 
 * @return the error status
 */

int spectra_pk(
	       struct background * pba,
	       struct perturbs * ppt,
	       struct primordial * ppm,
	       struct spectra * psp
	       ) {

  /** Summary: */

  /** - define local variables */

  int index_mode;
  int index_ic1,index_ic2,index_ic1_ic2;
  int index_k;
  int index_tau;
  double * primordial_pk; /* array with argument primordial_pk[index_ic_ic] */
  int last_index_back;
  double * pvecback_sp_long; /* array with argument pvecback_sp_long[pba->index_bg] */
  double source_ic1;
  double source_ic2;

  /** - check the presence of scalar modes */

  class_test((ppt->has_scalars == _FALSE_),
	     psp->error_message,
	     "you cannot ask for matter power spectrum since you turned off scalar modes");
  
  psp->index_md_scalars = ppt->index_md_scalars;
  index_mode = psp->index_md_scalars;

  /** - allocate temporary vectors where the primordial spectrum and the background quantitites will be stored */

  class_alloc(primordial_pk,psp->ic_ic_size[index_mode]*sizeof(double),psp->error_message);
  class_alloc(pvecback_sp_long,pba->bg_size*sizeof(double),psp->error_message);

  /** - allocate and fill array of P(k,tau) values */

  class_alloc(psp->ln_pk,sizeof(double)*psp->ln_tau_size*psp->ln_k_size*psp->ic_ic_size[index_mode],psp->error_message);

  for (index_tau=0 ; index_tau < psp->ln_tau_size; index_tau++) {

    class_call(background_at_tau(pba,
				 ppt->tau_sampling[index_tau-psp->ln_tau_size+ppt->tau_size], 
    /* for this last argument we could have passed 
       exp(psp->ln_tau[index_tau]) but we would then loose 
       precision in the exp(log(x)) operation) */
				 long_info, 
				 normal, 
				 &last_index_back, 
				 pvecback_sp_long),
	       pba->error_message,
	       psp->error_message);

    for (index_k=0; index_k<psp->ln_k_size; index_k++) {

      class_call(primordial_spectrum_at_k(ppm,index_mode,logarithmic,psp->ln_k[index_k],primordial_pk),
		 ppm->error_message,
		 psp->error_message);

      /* curvature primordial spectrum: 
	 P_R(k) = 1/(2pi^2) k^3 <R R>
	 so, primordial curvature correlator: 
	 <R R> = (2pi^2) k^-3 P_R(k) 
	 so, gravitational potential correlator:
	 <phi phi> = (2pi^2) k^-3 (source_phi)^2 P_R(k)

	 Matter power spectrum can be computed in two ways:

	 1) from source function = delta_pk (default)
	 P(k) = <delta_pk delta_pk>
	 = (source_delta_pk)^2 <R R>
	 = (2pi^2) k^-3 (source_delta_pk)^2 P_R(k)
	 
	 2) from source function = gravitational potential (using Poisson):
	 P(k) = <delta_m delta_m>
	 = 4/9 H^-4 Omega_m^-2 (k/a)^4 <phi phi>
	 = 4/9 H^-4 Omega_m^-2 (k/a)^4 (source_phi)^2 <R R> 
	 = 8pi^2/9 H^-4 Omega_m^-2 k/a^4 (source_phi)^2 <R R> 

	 CLASS now uses the first way. You could easily use the second
	 one (by uncommenting the lines below) if you prefer the
	 second way, and if you want to save a bit of computing time
	 (since in the second case we need to compute a single source
	 function for both P(k) and lensing). The result will
	 differ only on scales close to Hubble scale.
 
	 For isocurvature or cross adiabatic-isocurvature parts, 
	 replace one or two 'R' by 'S_i's */

      /* part diagonal in initial conditions */
      for (index_ic1 = 0; index_ic1 < psp->ic_size[index_mode]; index_ic1++) {
	
	index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_mode]);

	if (ppt->has_source_delta_pk == _TRUE_) {

	  source_ic1 = ppt->sources[index_mode]
	    [index_ic1 * ppt->tp_size[index_mode] + ppt->index_tp_delta_pk]
	    [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_mode] + index_k];
	  
	  psp->ln_pk[(index_tau * psp->ln_k_size + index_k)* psp->ic_ic_size[index_mode] + index_ic1_ic2] =
	  log(2.*_PI_*_PI_/exp(3.*psp->ln_k[index_k])
	      *source_ic1*source_ic1
	      *exp(primordial_pk[index_ic1_ic2]));
	}
	else {

	  source_ic1 = ppt->sources[index_mode]
	    [index_ic1 * ppt->tp_size[index_mode] + ppt->index_tp_g]
	    [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_mode] + index_k];
	
	  psp->ln_pk[(index_tau * psp->ln_k_size + index_k)* psp->ic_ic_size[index_mode] + index_ic1_ic2] =
	    log(8.*_PI_*_PI_/9./pow(pvecback_sp_long[pba->index_bg_H],4)
		/pow(pvecback_sp_long[pba->index_bg_Omega_m],2)
		*exp(psp->ln_k[index_k])
	      /pow(pvecback_sp_long[pba->index_bg_a],4)
		*exp(primordial_pk[index_ic1_ic2])*source_ic1*source_ic1);
	}
      }

      /* part non-diagonal in initial conditions */
      for (index_ic1 = 0; index_ic1 < psp->ic_size[index_mode]; index_ic1++) {
	for (index_ic2 = index_ic1+1; index_ic2 < psp->ic_size[index_mode]; index_ic2++) {
	  
	  index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_mode]);

	  if (psp->is_non_zero[index_mode][index_ic1_ic2] == _TRUE_) {

	    if (ppt->has_source_delta_pk == _TRUE_) {

	      source_ic1 = ppt->sources[index_mode]
		[index_ic1 * ppt->tp_size[index_mode] + ppt->index_tp_delta_pk]
		[(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_mode] + index_k];
	      
	      source_ic2 = ppt->sources[index_mode]
		[index_ic2 * ppt->tp_size[index_mode] + ppt->index_tp_delta_pk]
		[(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_mode] + index_k];
	      	      
	    }
	    else {
	      
	      source_ic1 = ppt->sources[index_mode]
		[index_ic1 * ppt->tp_size[index_mode] + ppt->index_tp_g]
		[(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_mode] + index_k];
	      
	      source_ic2 = ppt->sources[index_mode]
		[index_ic2 * ppt->tp_size[index_mode] + ppt->index_tp_g]
		[(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_mode] + index_k]; 
	    
	    }

	    psp->ln_pk[(index_tau * psp->ln_k_size + index_k)* psp->ic_ic_size[index_mode] + index_ic1_ic2] = 
	      primordial_pk[index_ic1_ic2]*sign(source_ic1)*sign(source_ic2);
	    
	  }
	  else {
	    psp->ln_pk[(index_tau * psp->ln_k_size + index_k)* psp->ic_ic_size[index_mode] + index_ic1_ic2] = 0.;
	  }
	}
      }
    }
  }

  /**- if interpolation of P(k,tau) will be needed (as a function of tau), 
     compute array of second derivatives in view of spline interpolation */  
 
  if (psp->ln_tau_size > 1) {

    class_alloc(psp->ddln_pk,sizeof(double)*psp->ln_tau_size*psp->ln_k_size*psp->ic_ic_size[index_mode],psp->error_message);

    class_call(array_spline_table_lines(psp->ln_tau,
					psp->ln_tau_size,
					psp->ln_pk,
					psp->ic_ic_size[index_mode]*psp->ln_k_size,
					psp->ddln_pk,
					_SPLINE_EST_DERIV_,
					psp->error_message),
	       psp->error_message,
	       psp->error_message);
    
  }
  
  free (primordial_pk);
  free (pvecback_sp_long);
  
  return _SUCCESS_;
}

/**
 * This routine computes a table of values for all matter power spectra P(k),
 * given the source functions and primordial spectra.
 * 
 * @param pba Input : pointer to background structure (will provide density of each species)
 * @param ppt Input : pointer to perturbation structure (contain source functions)
 * @param psp Input/Output: pointer to spectra structure 
 * @return the error status
 */

int spectra_matter_transfers(
	       struct background * pba,
	       struct perturbs * ppt,
	       struct spectra * psp
	       ) {

  /** Summary: */

  /** - define local variables */

  int index_mode;
  int index_ic;
  int index_k;
  int index_tau;
  int last_index_back;
  double * pvecback_sp_long; /* array with argument pvecback_sp_long[pba->index_bg] */
  double delta_i,rho_i;
  double delta_rho_tot,rho_tot;
  int n_ncdm;

  /** - check the presence of scalar modes */

  class_test((ppt->has_scalars == _FALSE_),
	     psp->error_message,
	     "you cannot ask for matter power spectrum since you turned off scalar modes");
  
  psp->index_md_scalars = ppt->index_md_scalars;
  index_mode = psp->index_md_scalars;

  /** - allocate and fill array of T_i(k,tau) values */

  class_alloc(psp->matter_transfer,sizeof(double)*psp->ln_tau_size*psp->ln_k_size*psp->ic_size[index_mode]*psp->tr_size,psp->error_message);

  /** - allocate temporary vectors where the background quantitites will be stored */

  class_alloc(pvecback_sp_long,pba->bg_size*sizeof(double),psp->error_message);

  for (index_tau=0 ; index_tau < psp->ln_tau_size; index_tau++) {

    class_call(background_at_tau(pba,
				 ppt->tau_sampling[index_tau-psp->ln_tau_size+ppt->tau_size], 
    /* for this last argument we could have passed 
       exp(psp->ln_tau[index_tau]) but we would then loose 
       precision in the exp(log(x)) operation) */
				 long_info, 
				 normal, 
				 &last_index_back, 
				 pvecback_sp_long),
	       pba->error_message,
	       psp->error_message);

    for (index_k=0; index_k<psp->ln_k_size; index_k++) {

      for (index_ic = 0; index_ic < psp->ic_size[index_mode]; index_ic++) {

	delta_rho_tot=0.;
	rho_tot=0.;

	/* T_g(k,tau) */

	delta_i = ppt->sources[index_mode]
	  [index_ic * ppt->tp_size[index_mode] + ppt->index_tp_delta_g]
	  [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_mode] + index_k];
	
	rho_i = pvecback_sp_long[pba->index_bg_rho_g];

	psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_mode] + index_ic) * psp->tr_size + psp->index_tr_g] = delta_i;

	delta_rho_tot += rho_i * delta_i;

	rho_tot += rho_i;

	/* T_b(k,tau) */

	delta_i = ppt->sources[index_mode]
	  [index_ic * ppt->tp_size[index_mode] + ppt->index_tp_delta_b]
	  [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_mode] + index_k];
	
	rho_i = pvecback_sp_long[pba->index_bg_rho_b];

	psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_mode] + index_ic) * psp->tr_size + psp->index_tr_b] = delta_i;

	delta_rho_tot += rho_i * delta_i;

	rho_tot += rho_i;
	
	/* T_cdm(k,tau) */
	
	if (pba->has_cdm == _TRUE_) {
	  
	  delta_i = ppt->sources[index_mode]
	    [index_ic * ppt->tp_size[index_mode] + ppt->index_tp_delta_cdm]
	    [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_mode] + index_k];
	  
	  rho_i = pvecback_sp_long[pba->index_bg_rho_cdm];

	psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_mode] + index_ic) * psp->tr_size + psp->index_tr_cdm] = delta_i;
	  
	  delta_rho_tot += rho_i * delta_i;
	  
	  rho_tot += rho_i;
	  
	}
	
	/* T_fld(k,tau) */

	if (pba->has_fld == _TRUE_) {

	  delta_i = ppt->sources[index_mode]
	    [index_ic * ppt->tp_size[index_mode] + ppt->index_tp_delta_fld]
	    [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_mode] + index_k];

	  rho_i = pvecback_sp_long[pba->index_bg_rho_fld];
	  
	  psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_mode] + index_ic) * psp->tr_size + psp->index_tr_de] = delta_i;
	  
	  delta_rho_tot += rho_i * delta_i;
	  
	  rho_tot += rho_i;

	}

	/* T_ur(k,tau) */

	if (pba->has_ur == _TRUE_) {
	
	  delta_i = ppt->sources[index_mode]
	    [index_ic * ppt->tp_size[index_mode] + ppt->index_tp_delta_ur]
	    [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_mode] + index_k];
	  
	  rho_i = pvecback_sp_long[pba->index_bg_rho_ur];
	  
	  psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_mode] + index_ic) * psp->tr_size + psp->index_tr_ur] = delta_i;
	  
	  delta_rho_tot += rho_i * delta_i;
	  
	  rho_tot += rho_i;

	}

	/* T_ncdm_i(k,tau) */

	if (pba->has_ncdm == _TRUE_) {

	  for (n_ncdm=0; n_ncdm < pba->N_ncdm; n_ncdm++) {

	    delta_i = ppt->sources[index_mode]
	      [index_ic * ppt->tp_size[index_mode] + ppt->index_tp_delta_ncdm1+n_ncdm]
	      [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_mode] + index_k];
	  
	    rho_i = pvecback_sp_long[pba->index_bg_rho_ncdm1+n_ncdm];
	  
	    psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_mode] + index_ic) * psp->tr_size + psp->index_tr_ncdm1+n_ncdm] = delta_i;
	  
	    delta_rho_tot += rho_i * delta_i;
	    
	    rho_tot += rho_i;
	  }
	}

	/* include homogeneous component in rho_tot */

	if (pba->has_lambda == _TRUE_) {

	  rho_i = pvecback_sp_long[pba->index_bg_rho_lambda];

	  rho_tot += rho_i;
	}

	/* T_tot(k,tau) */

	psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_mode] + index_ic) * psp->tr_size + psp->index_tr_tot] = delta_rho_tot/rho_tot;
	
      }
    }
  }

  /**- if interpolation of P(k,tau) will be needed (as a function of tau), 
     compute array of second derivatives in view of spline interpolation */  
 
  if (psp->ln_tau_size > 1) {

    class_alloc(psp->ddmatter_transfer,sizeof(double)*psp->ln_tau_size*psp->ln_k_size*psp->ic_size[index_mode]*psp->tr_size,psp->error_message);

    class_call(array_spline_table_lines(psp->ln_tau,
					psp->ln_tau_size,
					psp->matter_transfer,
					psp->ic_size[index_mode]*psp->ln_k_size*psp->tr_size,
					psp->ddmatter_transfer,
					_SPLINE_EST_DERIV_,
					psp->error_message),
	       psp->error_message,
	       psp->error_message);
    
  }

  free (pvecback_sp_long);

  return _SUCCESS_;
}
