/** @file background.c Documented background module
 *
 * Julien Lesgourgues, 27.08.2010    
 *
 * Deals with the cosmological background evolution.
 * This module has two purposes: 
 *
 * - at the beginning, to initialize the background, i.e. to integrate
 *    the background equations, and store all background quantities
 *    as a function of conformal time inside an interpolation table.
 *
 * - to provide routines which allow other modules to evaluate any
 *    background quantity for a given value of the conformal time (by
 *    interpolating within the interpolation table), or to find the
 *    correspondance between redhsift and conformal time.
 *
 *
 * The overall logic in this module is the following: 
 *
 * 1. most background parameters that we will call {A}
 * (e.g. rho_gamma, ..) can be expressed as simple analytical
 * functions of a few variables that we will call {B} (in simplest
 * models, of the scale factor 'a'; in extended cosmologies, of 'a'
 * plus e.g. (phi, phidot) for quintessence, or some temperature for
 * exotic particles, etc...).
 *
 * 2. in turn, quantitites {B} can be found as a function of conformal
 * time by integrating the background equations.
 * 
 * 3. some other quantitites that we will call {C} (like e.g. the
 * sound horizon or proper time) also require an integration with
 * respect to time, that cannot be infered analytically from
 * parameters {B}.
 *
 * So, we define the following routines:
 *
 * - background_functions() returns all background
 *    quantitites {A} as a function of quantitites {B}.
 * 
 * - background_solve() integrates the quantities {B} and {C} with
 *    respect to conformal time; this integration requires many calls
 *    to background_functions().
 *
 * - the result is stored in the form of a big table in the background
 *    structure. There is one column for conformal time 'eta'; one or
 *    more for quantitites {B}; then several columns for quantities {A}
 *    and {C}.
 *
 * Later in the code, if we know the variables {B} and need some
 * quantity {A}, the quickest and most procise way is to call directly
 * background_functions() (for instance, in simple models, if we want
 * H at a given value of the scale factor). If we know 'eta' and want
 * any other qunatity, we can call background_at_eta(), which
 * interpolates in the table and returns all values. Finally it can be
 * useful to get 'eta' for a given redshift 'z': this can be done with
 * background_eta_of_z(). So, in principle, if we know z, calling
 * background_eta_of_z() and then background_at_eta() will provide
 * roughly the same information as calling directly
 * background_functions() with a=a_0/(1+z).
 *
 *
 * In order to save time, background_at_eta() can be called in three
 * modes: short_info, normal_info, long_info (returning only essential
 * quantities, or useful quantitites, or rarely useful
 * quantities). Each line in the interpolation table is a vector which
 * first few elements correspond to the short_info format; a larger
 * fraction contribute to the normal format; and the full vector
 * corresponds to the long format. The guideline is that short_info
 * returns only geometric quantitites like a, H, H'; normal format
 * returns quantities strictly needed at each step in the integration
 * of perturbations; long_info returns quantitites needed only
 * occasionally.
 *
 * In summary, the following functions can be called from other modules:
 *
 * -# background_init() at the beginning  
 * -# background_at_eta(), background_functions(), background_eta_of_z() at any later time
 * -# background_free() at the end, when no more calls to the previous functions are needed
 */

#include "background.h"

/**
 * Background quantities at given conformal time eta.
 *
 * Evaluates all background quantities at a given value of
 * conformal time by reading the pre-computed table ant interpolating.
 *
 * @param pba           Input: pointer to background structure (containing pre-computed table)
 * @param eta           Input: value of conformal time
 * @param return_format Input: format of output vector (short, normal, long)
 * @param intermode     Input: interpolation mode (normal or closeby)
 * @param last_index    Input/Ouput: index of the previous/current point in the interpolation array (input only for closeby mode, output for both) 
 * @param pvecback      Output: vector (assumed to be already allocated)
 * @return the error status
 */

int background_at_eta(
		      struct background *pba,
		      double eta,
		      enum format_info return_format,
		      enum interpolation_mode intermode,
		      int * last_index,
		      double * pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size comptible with return_format) */
		      ) {

  /** Summary: */

  /** - define local variables */

  /* size of output vector, controlled by input parameter return_format */
  int pvecback_size;

  /** - check that eta is in the pre-computed range */

  class_test(eta < pba->eta_table[0],
	     pba->error_message,
	     "out of range: eta=%e < eta_min=%e, you should decrease the precision parameter a_ini_over_a_today_default\n",eta,pba->eta_table[0]);

  class_test(eta > pba->eta_table[pba->bt_size-1],
	     pba->error_message,
	     "out of range: eta=%e > eta_max=%e\n",eta,pba->eta_table[pba->bt_size-1]);

  /** - deduce length of returned vector from format mode */ 

  if (return_format == normal_info) {
    pvecback_size=pba->bg_size_normal;
  }
  else {
    if (return_format == short_info) {
      pvecback_size=pba->bg_size_short;
    }
    else { 
      pvecback_size=pba->bg_size;
    }
  }

  /** - interpolate from pre-computed table with array_interpolate()
      or array_interpolate_growing_closeby() (depending on
      interpolation mode) */

  if (intermode == normal) {
    class_call(array_interpolate_spline(
					pba->eta_table,
					pba->bt_size,
					pba->background_table,
					pba->d2background_deta2_table,
					pba->bg_size,
					eta,
					last_index,
					pvecback,
					pvecback_size,
					pba->error_message),
	       pba->error_message,
	       pba->error_message);
  }
  if (intermode == closeby) {
    class_call(array_interpolate_spline_growing_closeby(
							pba->eta_table,
							pba->bt_size,
							pba->background_table,
							pba->d2background_deta2_table,
							pba->bg_size,
							eta,
							last_index,
							pvecback,
							pvecback_size,
							pba->error_message),
	       pba->error_message,
	       pba->error_message);
  }

  return _SUCCESS_;
}

/** 
 * Conformal time at given redhsift.
 *
 * Returns eta(z) by interpolation from pre-computed table.
 *
 * @param pba Input: pointer to background structure
 * @param z   Input: redshift
 * @param eta Output: conformal time
 * @return the error status
 */

int background_eta_of_z(
			struct background *pba,
			double z,
			double * eta
			) {

  /** Summary: */

  /** - define local variables */

  /* necessary for calling array_interpolate(), but never used */
  int last_index; 

  /** - check that \f$ z \f$ is in the pre-computed range */
  class_test(z < pba->z_table[pba->bt_size-1],
	     pba->error_message,
	     "out of range: z=%e < z_min=%e\n",z,pba->z_table[pba->bt_size-1]);

  class_test(z > pba->z_table[0],
	     pba->error_message,
	     "out of range: a=%e > a_max=%e\n",z,pba->z_table[0]);

  /** - interpolate from pre-computed table with array_interpolate() */
  class_call(array_interpolate_spline(
				      pba->z_table,
				      pba->bt_size, 
				      pba->eta_table,
				      pba->d2eta_dz2_table,
				      1,
				      z,
				      &last_index,
				      eta,
				      1,
				      pba->error_message),
	     pba->error_message,
	     pba->error_message);

  return _SUCCESS_;
}

/**
 * Background quantities at given a.
 *
 * Function evaluating all background quantities which can be computed
 * analytically as a function of parameters like the scale factor 'a'
 * (see discussion at the beginnign of this file). In extended
 * comsological models, we would include other input parameters than
 * just 'a', e.g. (phi, phidot) for quintessence, some temperature of
 * exotic relics, etc...
 *
 * @param pba           Input: pointer to background structure
 * @param a             Input: value of scale factor 
 * @param return_format Input: format of output vector 
 * @param pvecback      Output: vector of background quantities (assmued to be already allocated) 
 * @return the error status
 */

int background_functions(
			 struct background *pba,
			 double a, /* in extended models there could be more than one argument: phi, phidot of quintessence; temperature of some particles; etc. */
			 enum format_info return_format,
			 double * pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size comptible with return_format) */
			 ) {
  
  /** Summary: */

  /** - define local variables */

  /* total density */
  double rho_tot;
  /* total pressure */
  double p_tot;
  /* total relativistic density */
  double rho_r;
  /* total non-relativistic density */
  double rho_m;
  /* scale factor relative to scale factor today */
  double a_rel;

  double rho_ncdm,p_ncdm,pseudo_p_ncdm;
  int n_ncdm;

  /** - initialize local variables */
  rho_tot = 0.;
  p_tot = 0.;
  rho_r=0.;
  rho_m=0.;
  a_rel = a / pba->a_today;

  class_test(a_rel <= 0.,
	     pba->error_message,
	     "a = %e instead of strictly positive",a_rel);

  /** - pass value of a to output */
  pvecback[pba->index_bg_a] = a;

  /** - compute each component's density and pressure */

  /* photons */
  pvecback[pba->index_bg_rho_g] = pba->Omega0_g * pow(pba->H0,2) / pow(a_rel,4);
  rho_tot += pvecback[pba->index_bg_rho_g];
  p_tot += (1./3.) * pvecback[pba->index_bg_rho_g];
  rho_r += pvecback[pba->index_bg_rho_g];

  /* baryons */
  pvecback[pba->index_bg_rho_b] = pba->Omega0_b * pow(pba->H0,2) / pow(a_rel,3);
  rho_tot += pvecback[pba->index_bg_rho_b];
  p_tot += 0;
  rho_m += pvecback[pba->index_bg_rho_b];

  /* cdm */
  if (pba->has_cdm == _TRUE_) {
    pvecback[pba->index_bg_rho_cdm] = pba->Omega0_cdm * pow(pba->H0,2) / pow(a_rel,3);
    rho_tot += pvecback[pba->index_bg_rho_cdm];
    p_tot += 0.;
    rho_m += pvecback[pba->index_bg_rho_cdm];
  }

  /* ncdm */
  if (pba->has_ncdm == _TRUE_) {
    for(n_ncdm=0; n_ncdm<pba->N_ncdm; n_ncdm++){
      /* Loop over species: */
      class_call(background_ncdm_momenta(
					 /* Only calculate for non-NULL pointers: */
					 pba->q_ncdm_bg[n_ncdm],
					 pba->w_ncdm_bg[n_ncdm],
					 pba->q_size_ncdm_bg[n_ncdm],
					 pba->M_ncdm[n_ncdm],
					 pba->factor_ncdm[n_ncdm],
					 1./a_rel-1.,
					 NULL,
					 &rho_ncdm,
					 &p_ncdm,
					 NULL,
					 &pseudo_p_ncdm),
		 pba->error_message,
		 pba->error_message);

      pvecback[pba->index_bg_rho_ncdm1+n_ncdm] = rho_ncdm;
      rho_tot += rho_ncdm;
      pvecback[pba->index_bg_p_ncdm1+n_ncdm] = p_ncdm;
      p_tot += p_ncdm;
      pvecback[pba->index_bg_pseudo_p_ncdm1+n_ncdm] = pseudo_p_ncdm;
      /* (3 p_ncdm1) is the "relativistic" contrinution to rho_ncdm1 */
      rho_r += 3.* p_ncdm;
      /* (rho_ncdm1 - 3 p_ncdm1) is the "non-relativistic" contribution to rho_ncdm1 */
      rho_m += rho_ncdm - 3.* p_ncdm;
    }
  }

  /* Lambda */
  if (pba->has_lambda == _TRUE_) {
    pvecback[pba->index_bg_rho_lambda] = pba->Omega0_lambda * pow(pba->H0,2);
    rho_tot += pvecback[pba->index_bg_rho_lambda];
    p_tot -= pvecback[pba->index_bg_rho_lambda];
  }

  /* dark energy fluid with constant w */
  if (pba->has_dark_energy_fluid == _TRUE_) {
    pvecback[pba->index_bg_rho_de] = pba->Omega0_de * pow(pba->H0,2) / pow(a_rel,3.*(1.+pba->w_de));
    rho_tot += pvecback[pba->index_bg_rho_de];
    p_tot += pba->w_de * pvecback[pba->index_bg_rho_de];
  }

  /* relativistic neutrinos (and all relativistic relics) */
  if (pba->has_nur == _TRUE_) {
    pvecback[pba->index_bg_rho_nur] = pba->Omega0_nur * pow(pba->H0,2) / pow(a_rel,4);
    rho_tot += pvecback[pba->index_bg_rho_nur];
    p_tot += (1./3.) * pvecback[pba->index_bg_rho_nur];
    rho_r += pvecback[pba->index_bg_rho_nur];
  }

  /** - compute expansion rate H from Friedmann equation */
  pvecback[pba->index_bg_H] = sqrt(rho_tot);

  /** - compute derivative of H with respect to conformal time */
  pvecback[pba->index_bg_H_prime] = - (3./2.) * (rho_tot + p_tot) * a;

  /** - compute relativistic density to total density ratio */
  pvecback[pba->index_bg_Omega_r] = rho_r / rho_tot;

  /** - compute other quantities in the exhaustive, redundent format: */
  if (return_format == long_info) {
    
    /** - compute critical density */
    pvecback[pba->index_bg_rho_crit] = rho_tot;
    class_test(pvecback[pba->index_bg_rho_crit] <= 0.,
	       pba->error_message,
	       "rho_tot = %e instead of strictly positive",pvecback[pba->index_bg_rho_crit]);

    /** - compute Omega's */

    /** - compute relativistic density to total density ratio */
    pvecback[pba->index_bg_Omega_m] = rho_m / rho_tot;
    
    /* one can put other variables here */
    /*  */
    /*  */

  }

  return _SUCCESS_;
  
}

/** 
 * Initialize the background structure, and in particular the
 * background interpolation table.
 * 
 * @param ppr Input : pointer to precision structure
 * @param pba Input/Output : pointer to initialized background structure
 * @return the error status
 */

int background_init(
		    struct precision * ppr,
		    struct background * pba
		    ) {

  /** Summary: */

  /** - local variables : */
  int n_ncdm;	
  double Omega0_tot;
  int filenum=0;

  if (pba->background_verbose > 0) {
    printf("Running CLASS version %s\n",_VERSION_);
    printf("Computing background\n");

    if (pba->N_ncdm > 0) {
      for (n_ncdm=0;n_ncdm<pba->N_ncdm; n_ncdm++) {

	if (pba->got_files[n_ncdm] == _TRUE_) {
	  printf(" -> ncdm species i=%d read from file %s\n",n_ncdm+1,pba->ncdm_psd_files+filenum*_ARGUMENT_LENGTH_MAX_);
	  filenum++;
	}

	printf(" -> ncdm species i=%d sampled with %d (resp. %d) points for purpose of background (resp. perturbation) integration\n",
	       n_ncdm+1,
	       pba->q_size_ncdm_bg[n_ncdm],
	       pba->q_size_ncdm[n_ncdm]);
      }
      if (ppr->evolver == rk)
	printf(" -> WARNING: you are using the ndf15 integrator, with ncdm species it is recommended to use the Runge-Kutta one (write evolver=0 in one of your input files)\n");
    }
  }

  /** - assign values to all indices in vectors of background quantities with background_indices()*/
  class_call(background_indices(pba),
	     pba->error_message,
	     pba->error_message);
  
  /** - control that cosmological parameter values make sense */

  /* H0 in Mpc^{-1} */
  class_test((pba->H0 < _H0_SMALL_)||(pba->H0 > _H0_BIG_),
	     pba->error_message,
	     "H0=%g out of bounds (%g<H0<%g) \n",pba->H0,_H0_SMALL_,_H0_BIG_);

  class_test(fabs(pba->h * 1.e5 / _c_  / pba->H0 -1.)>ppr->smallest_allowed_variation,
	     pba->error_message,
	     "inconsistency between Hubble and reduced Hubble parameters: you have H0=%f/Mpc=%fkm/s/Mpc, but h=%f",pba->H0,pba->H0/1.e5* _c_,pba->h);

  /* Tcmb in K */
  class_test((pba->Tcmb < _TCMB_SMALL_)||(pba->Tcmb > _TCMB_BIG_),
	     pba->error_message,
	     "Tcmb=%g out of bounds (%g<Tcmb<%g)",pba->Tcmb,_TCMB_SMALL_,_TCMB_BIG_);

  if ((pba->background_verbose > 0) && (pba->has_ncdm == _TRUE_)) {
    for (n_ncdm=0; n_ncdm < pba->N_ncdm; n_ncdm++) {
      printf(" -> non-cold dark matter species with i=%d has m_i = %e eV (so m_i / omega_i =%e eV)\n",
	     n_ncdm+1,
	     pba->m_ncdm_in_eV[n_ncdm],
	     pba->m_ncdm_in_eV[n_ncdm]/pba->Omega0_ncdm[n_ncdm]/pba->h/pba->h);
    }
  }

  /* curvature */
  Omega0_tot = pba->Omega0_g + pba->Omega0_b;
  if (pba->has_cdm == _TRUE_) {
    Omega0_tot += pba->Omega0_cdm;
  }
  if (pba->has_ncdm == _TRUE_) {
    Omega0_tot += pba->Omega0_ncdm_tot;
  }
  if (pba->has_lambda == _TRUE_) {
    Omega0_tot += pba->Omega0_lambda;
  }
  if (pba->has_dark_energy_fluid == _TRUE_) {
    Omega0_tot += pba->Omega0_de;
  }
  if (pba->has_nur == _TRUE_) {
    Omega0_tot += pba->Omega0_nur;
  }
  class_test(fabs(Omega0_tot-1.) > _TOLERANCE_ON_CURVATURE_,
	     pba->error_message,
	     "non zero spatial curvature not available (Omega0 = %g)",Omega0_tot);

  /* other quantities which would lead to segmentation fault if zero */
  class_test(pba->a_today <= 0,
	     pba->error_message,
	     "input a_today = %e instead of strictly positive",pba->a_today);

  class_test(_Gyr_over_Mpc_ <= 0,
	     pba->error_message,
	     "_Gyr_over_Mpc = %e instead of strictly positive",_Gyr_over_Mpc_);

  /** - integrate background, allocate and fill the background table with background_solve()*/
  class_call(background_solve(ppr,pba),
	     pba->error_message,
	     pba->error_message);

  return _SUCCESS_;

}

/**
 * Free all memory space allocated by background_init().
 * 
 *
 * @param pba Input : pointer to background structure (to be freed)
 * @return the error status
 */

int background_free(
		    struct background *pba
		    ) {
  int k;
  
  free(pba->eta_table);
  free(pba->z_table);
  free(pba->d2eta_dz2_table);
  free(pba->background_table);
  free(pba->d2background_deta2_table);

  if (pba->has_ncdm == _TRUE_) {
    for(k=0; k<pba->N_ncdm; k++){
      free(pba->q_ncdm[k]);
      free(pba->w_ncdm[k]);
      free(pba->q_ncdm_bg[k]);
      free(pba->w_ncdm_bg[k]);
      free(pba->dlnf0_dlnq_ncdm[k]);
    }
    free(pba->q_ncdm);
    free(pba->w_ncdm);
    free(pba->q_ncdm_bg);
    free(pba->w_ncdm_bg);
    free(pba->dlnf0_dlnq_ncdm);
    free(pba->q_size_ncdm);
    free(pba->q_size_ncdm_bg);
    free(pba->M_ncdm);
    free(pba->T_ncdm);
    free(pba->ksi_ncdm);
    free(pba->deg_ncdm);
    free(pba->Omega0_ncdm);
    free(pba->m_ncdm_in_eV);
    free(pba->factor_ncdm);
    if(pba->got_files!=NULL) free(pba->got_files);
    if(pba->ncdm_psd_files!=NULL)  free(pba->ncdm_psd_files);
  }

  return _SUCCESS_;
}

/**
 * Assign value to each relevant index in vectors of background quantities.
 *
 * @param pba Input : pointer to background structure
 * @return the error status
 */

int background_indices(
		       struct background *pba
		       ) {

  /** Summary: */

  /** - define local variables */

  /* a running index for the vector of background quantities */
  int index_bg;
  /* a running index for the vector of background quantities to be integrated */
  int index_bi;

  /** - intialization of all indices and flags */
  index_bg=0;
  index_bi=0;

  /* index for scale factor */
  pba->index_bg_a = index_bg; 
  index_bg++;

  /* - indices for H and its conformal-time-derivative */
  pba->index_bg_H = index_bg; 
  index_bg++;
  pba->index_bg_H_prime = index_bg; 
  index_bg++;

  /* - end of indices in the short vector of background values */
  pba->bg_size_short = index_bg;

  /* - index for rho_g (photon density) */
  pba->index_bg_rho_g = index_bg; 
  index_bg++;

  /* - index for rho_b (baryon density) */
  pba->index_bg_rho_b = index_bg; 
  index_bg++;

  /* - cdm */
  if (pba->Omega0_cdm != 0.) {
    pba->has_cdm = _TRUE_;
    /* -> index for rho_cdm (cdm density) */
    pba->index_bg_rho_cdm = index_bg; 
    index_bg++;
  }
  else { 
    pba->has_cdm = _FALSE_;
  }

  /* - ncdm1 */
  if (pba->Omega0_ncdm_tot != 0.) {
    pba->has_ncdm = _TRUE_;
    /* -> index for rho_ncdm1 (ncdm1 density) */
    pba->index_bg_rho_ncdm1 = index_bg; 
    index_bg +=pba->N_ncdm;
    pba->index_bg_p_ncdm1 = index_bg; 
    index_bg +=pba->N_ncdm;
    pba->index_bg_pseudo_p_ncdm1 = index_bg;
    index_bg +=pba->N_ncdm;
  }
  else { 
    pba->has_ncdm = _FALSE_;
  }
  
  /* - Lambda */
  if (pba->Omega0_lambda != 0.) {
    pba->has_lambda = _TRUE_;
    /* -> index for rho_Lambda (Lambda density) */
    pba->index_bg_rho_lambda = index_bg; 
    index_bg++;
  }
  else {
    pba->has_lambda = _FALSE_;
  }
  
  /* - Dark energy fluid */
  if (pba->Omega0_de != 0.) {
    pba->has_dark_energy_fluid = _TRUE_;
    /* -> index for rho_de (dark energy density) */
    pba->index_bg_rho_de = index_bg; 
    index_bg++;
  }
  else {
    pba->has_dark_energy_fluid = _FALSE_;
  }
  
  /* - relativistic neutrinos (and any relativistic relics) */
  if (pba->Omega0_nur != 0.) {
    pba->has_nur = _TRUE_;
    /* -> index for rho_nur (relativistic neutrino/relics density) */
    pba->index_bg_rho_nur = index_bg; 
    index_bg++;
  }
  else {
    pba->has_nur = _FALSE_;
  }

  /* - index for Omega_r (relativistic density fraction) */
  pba->index_bg_Omega_r = index_bg; 
  index_bg++;

  /* - put here additional ingredients that you want to appear in the
     normal vector */
  /*    */
  /*    */

  /* - end of indices in the normal vector of background values */
  pba->bg_size_normal = index_bg;

  /* - indices in the long version : */

  /* -> critical density */
  pba->index_bg_rho_crit = index_bg; 
  index_bg++;
  
  /* - index for Omega_ (non-relativistic density fraction) */
  pba->index_bg_Omega_m = index_bg; 
  index_bg++;

  /* -> conformal distance */
  pba->index_bg_conf_distance = index_bg; 
  index_bg++;

  /* -> proper time (for age of the Universe) */
  pba->index_bg_time = index_bg; 
  index_bg++;

  /* -> conformal sound horizon */
  pba->index_bg_rs = index_bg; 
  index_bg++;

  /* -> put here additional quantities describing background */
  /*    */
  /*    */

  /* -> end of indices in the long vector of background values */
  pba->bg_size = index_bg;

  /* - now, indices in vector of variables to integrate */

  /* -> scale factor */
  pba->index_bi_a = index_bi; 
  index_bi++;

  /* -> proper time (for age of the Universe) */
  pba->index_bi_time = index_bi; 
  index_bi++;

  /* -> sound horizon */
  pba->index_bi_rs = index_bi; 
  index_bi++;

  /* -> index for conformal time in vector of variables to integrate */
  pba->index_bi_eta = index_bi; 
  index_bi++;

  /* -> end of indices in the vector of variables to integrate */
  pba->bi_size = index_bi;

  /* index_bi_eta must be the last index, because eta is part of this vector for the purpose of being stored, */
  /* but it is not a quantity to be integrated (since integration is over eta itself) */      
  class_test(pba->index_bi_eta != index_bi-1,
	     pba->error_message,
	     "background integration requires index_bi_eta to be the last of all index_bi's");

  return _SUCCESS_;

}

int background_ncdm_distribution(
				 void * pbadist,
				 double q,
				 double * f0
				 ) {
  struct background * pba;
  struct background_parameters_for_distributions * pbadist_local;
  int n_ncdm,lastidx;
  double ksi,g;
  double qlast,dqlast,f0last,df0last;
  
  pbadist_local = pbadist;
  pba = pbadist_local->pba;

  n_ncdm = pbadist_local->n_ncdm;
  ksi = pba->ksi_ncdm[n_ncdm];
  g = pba->deg_ncdm[n_ncdm];
  
  /* Do we interpolate or use analytical formula below?  */
  if ((pba->got_files!=NULL)&&(pba->got_files[n_ncdm]==_TRUE_)){
    lastidx = pbadist_local->tablesize-1;
    if(q<pbadist_local->q[0]){
      //Handle q->0 case:
      *f0 = pbadist_local->f0[0];
    }
    else if(q>pbadist_local->q[lastidx]){
      //Handle q>qmax case (ensure continuous and derivable function with Boltzmann tail):
      qlast=pbadist_local->q[lastidx];
      f0last=pbadist_local->f0[lastidx];
      dqlast=qlast - pbadist_local->q[lastidx-1];
      df0last=f0last - pbadist_local->f0[lastidx-1];
       
      *f0 = f0last*exp(-(qlast-q)*df0last/f0last/dqlast);
    }
    else{
      //Do interpolation:
      class_call(array_interpolate_spline(
					  pbadist_local->q,
					  pbadist_local->tablesize,
					  pbadist_local->f0,
					  pbadist_local->d2f0,
					  1,
					  q,
					  &pbadist_local->last_index,
					  f0,
					  1,
					  pba->error_message),
		 pba->error_message,     pba->error_message);
    }
    *f0 = *f0 *g;
  }
  else{
    //FD distribution:
    *f0 = g/pow(2*_PI_,3)*(1./(exp(q-ksi)+1.) +1./(exp(q+ksi)+1.));
  }

  return _SUCCESS_;
}

int background_ncdm_test_function(
				  void * pbadist,
				  double q,
				  double * test
				  ) {
  double zeta3=1.2020569031595942853997381615114499907649862923404988817922;
  double zeta5=1.0369277551433699263313654864570341680570809195019128119741;
  double a=1.0/log(2), b=12.0/(_PI_*_PI_), c=2.0/(3.0*zeta3);
  double d=120.0/(7.0*pow(_PI_,4)), e=2.0/(45.0*zeta5);




  *test = pow(2.0*_PI_,3)/10.0*(a+b*q+c*q*q+d*q*q*q+e*q*q*q*q);

  return _SUCCESS_;
}

int background_ncdm_init(
			 struct precision *ppr,
			 struct background *pba
			 ) {
  
  int index_q, k,tolexp,row,status,filenum;
  double f0m2,f0m1,f0,f0p1,f0p2,dq,q,df0dq,tmp1,tmp2;
  struct background_parameters_for_distributions pbadist;
  FILE *psdfile;

  pbadist.pba = pba;

  /* Allocate pointer arrays: */
  class_alloc(pba->q_ncdm, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->w_ncdm, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->q_ncdm_bg, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->w_ncdm_bg, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->dlnf0_dlnq_ncdm, sizeof(double*)*pba->N_ncdm,pba->error_message);

  /* Allocate pointers: */
  class_alloc(pba->q_size_ncdm,sizeof(int)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->q_size_ncdm_bg,sizeof(int)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->factor_ncdm,sizeof(double)*pba->N_ncdm,pba->error_message);

  for(k=0, filenum=0; k<pba->N_ncdm; k++){
    pbadist.n_ncdm = k;
    pbadist.q = NULL;
    pbadist.tablesize = 0;
    /*Do we need to read in a file to interpolate the distribution function? */
    if ((pba->got_files!=NULL)&&(pba->got_files[k]==_TRUE_)){
      psdfile = fopen(pba->ncdm_psd_files+filenum*_ARGUMENT_LENGTH_MAX_,"r");
      class_test(psdfile == NULL,pba->error_message,
		 "Could not open file %s!",pba->ncdm_psd_files+filenum*_ARGUMENT_LENGTH_MAX_);
      // Find size of table:
      for (row=0,status=2; status==2; row++){
	status = fscanf(psdfile,"%lf %lf",&tmp1,&tmp2);
      }
      rewind(psdfile);
      pbadist.tablesize = row-1;

      /*Allocate room for interpolation table: */
      class_alloc(pbadist.q,sizeof(double)*pbadist.tablesize,pba->error_message);
      class_alloc(pbadist.f0,sizeof(double)*pbadist.tablesize,pba->error_message);
      class_alloc(pbadist.d2f0,sizeof(double)*pbadist.tablesize,pba->error_message);
      for (row=0; row<pbadist.tablesize; row++){
	status = fscanf(psdfile,"%lf %lf",
			&pbadist.q[row],&pbadist.f0[row]);
	//		printf("(q,f0) = (%g,%g)\n",pbadist.q[row],pbadist.f0[row]);
      }
      fclose(psdfile);
      /* Call spline interpolation: */
      class_call(array_spline_table_lines(pbadist.q,
					  pbadist.tablesize,
					  pbadist.f0,
					  1,
					  pbadist.d2f0,
					  _SPLINE_EST_DERIV_,
					  pba->error_message), 
		 pba->error_message,
		 pba->error_message);
      filenum++;
    }

    /* Handle perturbation qsampling: */
    class_alloc(pba->q_ncdm[k],_QUADRATURE_MAX_*sizeof(double),pba->error_message);
    class_alloc(pba->w_ncdm[k],_QUADRATURE_MAX_*sizeof(double),pba->error_message); 	

    class_call(get_qsampling(pba->q_ncdm[k],
			     pba->w_ncdm[k],
			     &(pba->q_size_ncdm[k]),
			     _QUADRATURE_MAX_,
			     ppr->tol_ncdm,
			     pbadist.q,
			     pbadist.tablesize,
			     background_ncdm_test_function,
			     background_ncdm_distribution,
			     &pbadist,
			     pba->error_message),
	       pba->error_message,
	       pba->error_message);
    pba->q_ncdm[k]=realloc(pba->q_ncdm[k],pba->q_size_ncdm[k]*sizeof(double));
    pba->w_ncdm[k]=realloc(pba->w_ncdm[k],pba->q_size_ncdm[k]*sizeof(double));
  

    if (pba->background_verbose > 0)
      printf("ncdm species i=%d sampled with %d points for purpose of perturbation integration\n",
	     k+1,
	     pba->q_size_ncdm[k]);

    /* Handle background q_sampling: */
    class_alloc(pba->q_ncdm_bg[k],_QUADRATURE_MAX_BG_*sizeof(double),pba->error_message);
    class_alloc(pba->w_ncdm_bg[k],_QUADRATURE_MAX_BG_*sizeof(double),pba->error_message);

    class_call(get_qsampling(pba->q_ncdm_bg[k],
			     pba->w_ncdm_bg[k],
			     &(pba->q_size_ncdm_bg[k]),
			     _QUADRATURE_MAX_BG_,
			     ppr->tol_ncdm_bg,
			     pbadist.q,
			     pbadist.tablesize,
			     background_ncdm_test_function,
			     background_ncdm_distribution,
			     &pbadist,
			     pba->error_message),
	       pba->error_message,
	       pba->error_message);

    
    pba->q_ncdm_bg[k]=realloc(pba->q_ncdm_bg[k],pba->q_size_ncdm_bg[k]*sizeof(double));
    pba->w_ncdm_bg[k]=realloc(pba->w_ncdm_bg[k],pba->q_size_ncdm_bg[k]*sizeof(double));

    if (pba->background_verbose > 0)
      printf("ncdm species i=%d sampled with %d points for purpose of background integration\n",
	     k+1,
	     pba->q_size_ncdm_bg[k]);

    class_alloc(pba->dlnf0_dlnq_ncdm[k],
		pba->q_size_ncdm[k]*sizeof(double),
		pba->error_message);


    for (index_q=0; index_q<pba->q_size_ncdm[k]; index_q++) {
      q = pba->q_ncdm[k][index_q];
      class_call(background_ncdm_distribution(&pbadist,q,&f0),
		 pba->error_message,pba->error_message);
  
      //Loop to find appropriate dq:
      for(tolexp=_PSD_DERIVATIVE_EXP_MIN_; tolexp<_PSD_DERIVATIVE_EXP_MAX_; tolexp++){
   
        if (index_q == 0){
          dq = min((0.5-ppr->smallest_allowed_variation)*q,2*exp(tolexp)*(pba->q_ncdm[k][index_q+1]-q));
        }       
        else if (index_q == pba->q_size_ncdm[k]-1){
          dq = exp(tolexp)*2.0*(pba->q_ncdm[k][index_q]-pba->q_ncdm[k][index_q-1]);
        }       
        else{
          dq = exp(tolexp)*(pba->q_ncdm[k][index_q+1]-pba->q_ncdm[k][index_q-1]);
        }

	class_call(background_ncdm_distribution(&pbadist,q-2*dq,&f0m2),
		   pba->error_message,pba->error_message);
	class_call(background_ncdm_distribution(&pbadist,q+2*dq,&f0p2),
		   pba->error_message,pba->error_message);

	if (fabs((f0p2-f0m2)/f0)>sqrt(ppr->smallest_allowed_variation)) break;
      }

      class_call(background_ncdm_distribution(&pbadist,q-dq,&f0m1),
		 pba->error_message,pba->error_message);
      class_call(background_ncdm_distribution(&pbadist,q+dq,&f0p1),
		 pba->error_message,pba->error_message);
      //5 point estimate of the derivative:    
      df0dq = (+f0m2-8*f0m1+8*f0p1-f0p2)/12.0/dq;
      //printf("df0dq[%g] = %g. dlf=%g ?= %g. f0 =%g.\n",q,df0dq,q/f0*df0dq,
      //-q/(1.0+exp(-q)),f0);
      //Avoid underflow in extreme tail:
      if (fabs(f0)==0.)
	pba->dlnf0_dlnq_ncdm[k][index_q] = -q; /* valid for whatever f0 with exponential tail in exp(-q) */
      else
   	pba->dlnf0_dlnq_ncdm[k][index_q] = q/f0*df0dq;
    }

    pba->factor_ncdm[k]=4*_PI_*pow(pba->Tcmb*pba->T_ncdm[k]*_k_B_,4)*8*_PI_*_G_
      /3./pow(_h_P_/2./_PI_,3)/pow(_c_,7)*_Mpc_over_m_*_Mpc_over_m_;

    /* If allocated, deallocate interpolation table:  */
    if ((pba->got_files!=NULL)&&(pba->got_files[k]==_TRUE_)){
      free(pbadist.q);
      free(pbadist.f0);
      free(pbadist.d2f0);
    }
  }


  return _SUCCESS_;
}

int background_ncdm_momenta(
			    /* Only calculate for non-NULL pointers: */
			    double * qvec,
			    double * wvec,
			    int qsize,
			    double M,
			    double factor,
			    double z,
			    double * n,
			    double * rho,
			    double * p,
			    double * drho_dM,
			    double * pseudo_p
			    ) {

  int index_q;
  double epsilon;
  double q2;
  double factor2;

  factor2 = factor*pow(1+z,4);

  if (n!=NULL) *n = 0.;
  if (rho!=NULL) *rho = 0.;
  if (p!=NULL) *p = 0.;
  if (drho_dM!=NULL) *drho_dM = 0.;
  if (pseudo_p!=NULL) *pseudo_p = 0.;

  for (index_q=0; index_q<qsize; index_q++) {

    q2 = qvec[index_q]*qvec[index_q];
    epsilon = sqrt(q2+M*M/(1.+z)/(1.+z));

    if (n!=NULL) *n += q2*wvec[index_q];
    if (rho!=NULL) *rho += q2*epsilon*wvec[index_q];
    if (p!=NULL) *p += q2*q2/3./epsilon*wvec[index_q];
    if (drho_dM!=NULL) *drho_dM += q2*M/(1.+z)/(1.+z)/epsilon*wvec[index_q];
    if (pseudo_p!=NULL) *pseudo_p += pow(q2/epsilon,3)/3.0*wvec[index_q]; 
  }
  if (n!=NULL) *n *= factor2;
  if (rho!=NULL) *rho *= factor2;
  if (p!=NULL) *p *= factor2;
  if (drho_dM!=NULL) *drho_dM *= factor2; 
  if (pseudo_p!=NULL) *pseudo_p *=factor2;

  return _SUCCESS_;
}

int background_ncdm_M_from_Omega(
				 struct precision *ppr,
				 struct background *pba,
				 int k
				 ) {
  double rho0,rho,n,M,deltaM,drhodM;
  int iter,maxiter=50;
	
  rho0 = pba->H0*pba->H0*pba->Omega0_ncdm[k]; /*Remember that rho is defined such that H^2=sum(rho_i) */
  M = 0.0;

  background_ncdm_momenta(pba->q_ncdm_bg[k],
			  pba->w_ncdm_bg[k],
			  pba->q_size_ncdm_bg[k],
			  M,
			  pba->factor_ncdm[k],
              		  0.,
			  &n,
			  &rho,
			  NULL,
			  NULL,
			  NULL);

  /* Is the value of Omega less than a massless species?*/
  class_test(rho0<rho,pba->error_message,
	     "The value of Omega for the %dth species, %g, is less than for a massless species! It should be atleast %g. Check your input.",
	     k,pba->Omega0_ncdm[k],pba->Omega0_ncdm[k]*rho/rho0);

  /* In the strict NR limit we have rho = n*(M) today, giving a zero'th order guess: */
  M = rho0/n; /* This is our guess for M. */
  for (iter=1; iter<=maxiter; iter++){
    //    printf("iter=%d, M=%g, drhodM=%g, err:%g.\n",iter,M,drhodM,deltaM/M);
    /* Newton iteration. First get relevant quantities at M: */
    background_ncdm_momenta(pba->q_ncdm_bg[k],
			    pba->w_ncdm_bg[k],
			    pba->q_size_ncdm_bg[k],
			    M,
			    pba->factor_ncdm[k],
			    0.,
			    NULL,
			    &rho,
			    NULL,
			    &drhodM,
			    NULL);

    deltaM = (rho0-rho)/drhodM; /* By definition of the derivative */
    if ((M+deltaM)<0.0) deltaM = -M/2.0; /* Avoid overshooting to negative M value. */
    M += deltaM; /* Update value of M.. */
    if (fabs(deltaM/M)<ppr->tol_M_ncdm){
      /* Accuracy reached.. */
      pba->M_ncdm[k] = M;
      break;
    }
  }
  class_test(iter>=maxiter,pba->error_message,
	     "Newton iteration could not converge on a mass for some reason.");
  return _SUCCESS_;
}

int background_solve(
		     struct precision *ppr,
		     struct background *pba
		     ) {

  /** Summary: */

  /** - define local variables */

  /* contains all quantities relevant for the integration algorithm */
  struct generic_integrator_workspace gi;
  /* parameters and workspace for the background_derivs function */
  struct background_parameters_and_workspace bpaw;
  /* a growing table (since the number of time steps is not known a priori) */
  growTable gTable;
  /* needed for growing table */
  double * pData;
  /* needed for growing table */
  void * memcopy_result;
  /* initial conformal time */
  double eta_start;
  /* final conformal time */
  double eta_end;
  /* an index running over bi indices */
  int i;
  /* vector of quantities to be integrated */
  double * pvecback_integration;
  
  double * pvecback;

  int last_index=0; /* necessary for calling array_interpolate(), but never used */

  bpaw.pba = pba;
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);
  bpaw.pvecback = pvecback;  

  /** - allocate vector of quantities to be integrated */
  class_alloc(pvecback_integration,pba->bi_size*sizeof(double),pba->error_message);

  /** - initialize generic integrator with initialize_generic_integrator() */ 

  /* Size of vector to integrate is (pba->bi_size-1) rather than
   * (pba->bi_size), since eta is not integrated.
   */
  class_call(initialize_generic_integrator((pba->bi_size-1),&gi),
	     gi.error_message,
	     pba->error_message);

  /** - impose initial conditions with background_initial_conditions() */
  class_call(background_initial_conditions(ppr,pba,pvecback,pvecback_integration),
	     pba->error_message,
	     pba->error_message);

  /* here eta_end is in fact the initial time (in the next loop
     eta_start = eta_end) */
  eta_end=pvecback_integration[pba->index_bi_eta];

  /** - create a growTable with gt_init() */
  class_call(gt_init(&gTable),
	     gTable.error_message,
	     pba->error_message);
  
  /* initialize the counter for the number of steps */
  pba->bt_size=0;

  /** - loop over integration steps : call background_functions(), find step size, save data in growTable with gt_add(), perform one step with generic_integrator(), store new value of eta */

  while (pvecback_integration[pba->index_bi_a] < pba->a_today) {

    eta_start = eta_end;

    /* -> find step size (trying to adjust the last step as close as possible to the one needed to reach a=a_today; need not be exact, difference corrected later) */
    class_call(background_functions(pba,pvecback_integration[pba->index_bi_a], short_info, pvecback),
	       pba->error_message,
	       pba->error_message);

    if ((pvecback_integration[pba->index_bi_a]*(1.+ppr->back_integration_stepsize)) < pba->a_today) {
      eta_end = eta_start + ppr->back_integration_stepsize / (pvecback_integration[pba->index_bi_a]*pvecback[pba->index_bg_H]); 
      /* no possible segmentation fault here: non-zeroness of "a" has been checked in background_functions() */
    }
    else {
      eta_end = eta_start + (1./pvecback_integration[pba->index_bi_a]-1.) / (pvecback_integration[pba->index_bi_a]*pvecback[pba->index_bg_H]);  
      /* no possible segmentation fault here: non-zeroness of "a" has been checked in background_functions() */
    }

    class_test((eta_end-eta_start)/eta_start < ppr->smallest_allowed_variation,
	       pba->error_message,
	       "integration step: relative change in time =%e < machine precision : leads either to numerical error or infinite loop",(eta_end-eta_start)/eta_start);

    /* -> save data in growTable */
    class_call(gt_add(&gTable,_GT_END_,(void *) pvecback_integration,sizeof(double)*pba->bi_size),
	       gTable.error_message,
	       pba->error_message);
    pba->bt_size++;

    /* -> perform one step */
    class_call(generic_integrator(background_derivs,
				  eta_start,
				  eta_end,
				  pvecback_integration,
				  &bpaw,
				  ppr->tol_background_integration,
				  ppr->smallest_allowed_variation,
				  &gi),
	       gi.error_message,
	       pba->error_message);
    
    /* -> store value of eta */
    pvecback_integration[pba->index_bi_eta]=eta_end;

  }

  /** - save last data in growTable with gt_add() */
  class_call(gt_add(&gTable,_GT_END_,(void *) pvecback_integration,sizeof(double)*pba->bi_size),
	     gTable.error_message,
	     pba->error_message);
  pba->bt_size++;


  /* integration finished */

  /** - clean up generic integrator with cleanup_generic_integrator() */
  class_call(cleanup_generic_integrator(&gi),
	     gi.error_message,
	     pba->error_message);

  /** - retrieve data stored in the growTable with gt_getPtr() */
  class_call(gt_getPtr(&gTable,(void**)&pData),
	     gTable.error_message,
	     pba->error_message);

  /** - interpolate to get quantities precisely today with array_interpolate() */
  class_call(array_interpolate(
			       pData,
			       pba->bi_size,
			       pba->bt_size,
			       pba->index_bi_a, 
			       pba->a_today,
			       &last_index,
			       pvecback_integration,
			       pba->bi_size,
			       pba->error_message),
	     pba->error_message,
	     pba->error_message);

  /* substitute last line with quantities today */
  for (i=0; i<pba->bi_size; i++) 
    pData[(pba->bt_size-1)*pba->bi_size+i]=pvecback_integration[i];

  /** - deduce age of the Universe */
  /* -> age in Gyears */
  pba->age = pvecback_integration[pba->index_bi_time]/_Gyr_over_Mpc_;
  /* -> conformal age in Mpc */
  pba->conformal_age = pvecback_integration[pba->index_bi_eta];

  /** - allocate background tables */
  class_alloc(pba->eta_table,pba->bt_size * sizeof(double),pba->error_message);

  class_alloc(pba->z_table,pba->bt_size * sizeof(double),pba->error_message);

  class_alloc(pba->d2eta_dz2_table,pba->bt_size * sizeof(double),pba->error_message);

  class_alloc(pba->background_table,pba->bt_size * pba->bg_size * sizeof(double),pba->error_message);

  class_alloc(pba->d2background_deta2_table,pba->bt_size * pba->bg_size * sizeof(double),pba->error_message);
  
  /** - In a loop over lines, fill background table using the result of the integration plus background_functions() */
  for (i=0; i < pba->bt_size; i++) {
    
    /* -> establish correspondance between the integrated variable and the bg variables */

    pba->eta_table[i] = pData[i*pba->bi_size+pba->index_bi_eta];

    class_test(pData[i*pba->bi_size+pba->index_bi_a] <= 0.,
	       pba->error_message,
	       "a = %e instead of strictly positiv",pData[i*pba->bi_size+pba->index_bi_a]);

    pba->z_table[i] = pba->a_today/pData[i*pba->bi_size+pba->index_bi_a]-1.;

    pvecback[pba->index_bg_time] = pData[i*pba->bi_size+pba->index_bi_time];
    pvecback[pba->index_bg_conf_distance] = pba->conformal_age - pData[i*pba->bi_size+pba->index_bi_eta];
    pvecback[pba->index_bg_rs] = pData[i*pba->bi_size+pba->index_bi_rs];

    /* -> compute all other quantities */
    class_call(background_functions(pba,pData[i*pba->bi_size+pba->index_bi_a], long_info, pvecback),
	       pba->error_message,
	       pba->error_message);
    
    /* -> write in the table */
    memcopy_result = memcpy(pba->background_table + i*pba->bg_size,pvecback,pba->bg_size*sizeof(double));

    class_test(memcopy_result != pba->background_table + i*pba->bg_size,
	       pba->error_message,
	       "cannot copy data back to pba->background_table");
  }

  /** - free the growTable with gt_free() */

  class_call(gt_free(&gTable),
	     gTable.error_message,
	     pba->error_message);
  
  /** - fill tables of second derivatives (in view of spline interpolation) */
  class_call(array_spline_table_lines(pba->z_table,
				      pba->bt_size,
				      pba->eta_table,
				      1,
				      pba->d2eta_dz2_table,
				      _SPLINE_EST_DERIV_,
				      pba->error_message),
	     pba->error_message,
	     pba->error_message);

  class_call(array_spline_table_lines(pba->eta_table,
				      pba->bt_size,
				      pba->background_table,
				      pba->bg_size,
				      pba->d2background_deta2_table,
				      _SPLINE_EST_DERIV_,
				      pba->error_message),
	     pba->error_message,
	     pba->error_message);

  if (pba->background_verbose > 0) {
    printf(" -> age = %f Gyr\n",pba->age);
    printf(" -> conformal age = %f Mpc\n",pba->conformal_age);    
  }

  free(pvecback);
  free(pvecback_integration);

  return _SUCCESS_;
  
}

/**
 * Assign initial values to background integrated variables.
 *
 * @param ppr                  Input : pointer to precision structure
 * @param pba                  Input : pointer to background structure
 * @param pvecback             Input : vector of background quantitites used as workspace
 * @param pvecback_integration Output : vector of background quantitites to be integrated, returned with proper initial values
 * @return the error status
 */

int background_initial_conditions(
				  struct precision *ppr,
				  struct background *pba,
				  double * pvecback, /* vector with argument pvecback[index_bg] (must be already allocated, normal format is sufficient) */
				  double * pvecback_integration /* vector with argument pvecback_integration[index_bi] (must be already allocated with size pba->bi_size) */
				  ) {

  /** Summary: */

  /** - define local variables */

  /* scale factor */
  double a;

  int counter,is_early_enough,n_ncdm;

  /** - fix initial value of \f$ a \f$ */
  a = ppr->a_ini_over_a_today_default * pba->a_today / _SCALE_BACK_;

  /* test the validity of this choice (e.g.: are massive neutrinos
     relativistic? etc.) */

  for (counter=0, is_early_enough = _FALSE_; 
       (counter < _MAX_IT_) && (is_early_enough == _FALSE_ ); 
       counter++) {

    a *= _SCALE_BACK_;

    /** - compute initial H with background_functions() */
    class_call(background_functions(pba,a, normal_info, pvecback),
	       pba->error_message,
	       pba->error_message);

    is_early_enough = _TRUE_;

    /* check that we are deep inside radiation domination, in order to
       use the approximation \f$ t=1/(2H) \f$ below */
    if (fabs(pvecback[pba->index_bg_Omega_r]-1.) > ppr->tol_initial_Omega_r)
      is_early_enough = _FALSE_;

    /* check that non-cold relics are relativistic, in view of
       integrating their perturbations starting from correct initial
       conditions */
    if (pba->has_ncdm == _TRUE_) {
      for (n_ncdm=0; n_ncdm<pba->N_ncdm; n_ncdm++) {
	if (fabs(pvecback[pba->index_bg_p_ncdm1+n_ncdm]/pvecback[pba->index_bg_rho_ncdm1+n_ncdm]-1./3.) > ppr->tol_ncdm_initial_w)
	  is_early_enough = _FALSE_;					     
      }
    }
  }

  pvecback_integration[pba->index_bi_a] = a;

  /** - compute initial proper time, assuming radiation-dominated
      universe since Big Bang and therefore \f$ t=1/(2H) \f$ (good
      approximation for most purposes) */

  class_test(pvecback[pba->index_bg_H] <= 0.,
	     pba->error_message,
	     "H = %e instead of strictly positive",pvecback[pba->index_bg_H]);

  pvecback_integration[pba->index_bi_time] = 1./(2.* pvecback[pba->index_bg_H]);

  /** - compute initial conformal time, assuming radiation-dominated
      universe since Big Bang and therefore \f$ \eta=1/(aH) \f$
      (good approximation for most purposes) */
  pvecback_integration[pba->index_bi_eta] = 1./(a * pvecback[pba->index_bg_H]);



  /** - compute initial sound horizon, assuming c_s=1/sqrt(3) initially */
  pvecback_integration[pba->index_bi_rs] = pvecback_integration[pba->index_bi_eta]/sqrt(3.);

  return _SUCCESS_;

}

/** 
 * Subroutine evaluating the derivative with respect to conformal time
 * of quantities which are integrated (a, t, etc).
 *
 * This is one of the few functions in the code which are passed to
 * the generic_integrator() routine.  Since generic_integrator()
 * should work with functions passed from various modules, the format
 * of the arguments is a bit special:
 *
 * - fixed input parameters and wokspaces are passed through a generic
 * pointer. Here, this is just a pointer to the background structure
 * and to a background vector, but generic_integrator() doesn't know
 * its fine structure.  
 *
 * - the error management is a bit special: errors are not written as
 * usual to pba->error_message, but to a generic error_message passed
 * in the list of arguments.
 *
 * @param eta                      Input : conformal time
 * @param y                        Input : vector of variable
 * @param dy                       Output : its derivative (already allocated)
 * @param parameters_and_workspace Input: pointer to fixed parameters (e.g. indices)
 * @param error_message            Output : error message
 */
int background_derivs(
		      double eta,
		      double* y, /* vector with argument y[index_bi] (must be already allocated with size pba->bi_size) */
		      double* dy, /* vector with argument dy[index_bi]
				     (must be already allocated with
				     size pba->bi_size) */
		      void * parameters_and_workspace,
		      ErrorMsg error_message
		      ) {

  /** Summary: */

  /** - define local variables */

  struct background_parameters_and_workspace * pbpaw;
  struct background * pba;
  double * pvecback;

  pbpaw = parameters_and_workspace;
  pba =  pbpaw->pba;
  pvecback = pbpaw->pvecback;

  /** - Calculates functions of /f$ a /f$ with background_functions() */
  class_call(background_functions((struct background *)pba,y[pba->index_bi_a], normal_info, pvecback),
	     pba->error_message,
	     error_message);

  /** - calculate /f$ a'=a^2 H /f$ */
  dy[pba->index_bi_a] = y[pba->index_bi_a] * y[pba->index_bi_a] * pvecback[pba->index_bg_H];

  /** - calculate /f$ t' = a /f$ */
  dy[pba->index_bi_time] = y[pba->index_bi_a];

  class_test(pvecback[pba->index_bg_rho_g] <= 0.,
	     error_message,
	     "rho_g = %e instead of strictly positive",pvecback[pba->index_bg_rho_g]);

  dy[pba->index_bi_rs] = 1./sqrt(3.*(1.+3.*pvecback[pba->index_bg_rho_b]/4./pvecback[pba->index_bg_rho_g]));

  return _SUCCESS_;

}
