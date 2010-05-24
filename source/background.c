/** @file background.c Documented background module
 * Julien Lesgourgues, 18.04.2010    
 *
 * Deals with the cosmological background evolution issues.
 * This module has two purposes: 
 *
 * - at the beginning, to initialize the background, i.e. to integrate
     the background equations, and store all background quantities
     as a function of time and scale factor, inside an interpolation
     table.
 *
 * - at any time in the code, to evaluate any background quantity for
     a given value of the scale factor, redshift or conformal time (by
     interpolating within the interpolation table).
 *
 * Hence the following functions can be called from other modules:
 *
 * -# background_init() at the beginning  
 * -# background_at_eta(), background_functions_of_a(), background_eta_of_z() at any later time
 * -# background_free() at the end, when no more calls to the previous functions are needed
 */

#include "background.h"


double * pvecback; /**< vector of background quantities, used
		      throughout the background module. Use a global
		      variable in order to avoid reallocating it many
		      times. */

//@}

/** @name - miscellaneous: */

//@{

ErrorMsg Transmit_Error_Message; /**< contains error message */

//@}

  /**
   * Background quantities at given conformal time eta.
   *
  * Evaluates all background quantities at a given value of
  * conformal time by reading the pre-computed table ant interpolating.
  * This function can be called from whatever module at whatever time,
  * provided that background_init() has been called before, and 
  * background_free() has not been called yet.
  *
  * @param eta Input: value of conformal time
  * @param return_format Input: format of output vector
  * @param intermode Input: interpolation mode (normal or growing_closeby)
  * @param last_index Input/Ouput: index of the previous/current point in the interpolation array (input only for closeby mode, output for both) 
  * @param pvecback_local Output: vector (assumed to be already allocated)
  * @return the error status
  */
int background_at_eta(
		      struct background *pba,
		      double eta,
		      enum format_info return_format,
		      enum interpolation_mode intermode,
		      int * last_index,
		      double * pvecback_local
		      ) {

  /** Summary: */

  /** - define local variables */

  /* size of output vector, controlled by input parameter return_format */
  int pvecback_size;

  /** - check that eta is in the pre-computed range */

  if (eta < pba->eta_table[0]) {
    sprintf(pba->error_message,"%s(L:%d) : eta=%e < eta_min=%e, you should decrease the precision parameter a_ini_over_a_today_default\n",__func__,__LINE__,eta,pba->eta_table[0]);
    return _FAILURE_;
  }
  if (eta > pba->eta_table[pba->bt_size-1]) {
    sprintf(pba->error_message,"%s(L:%d) : eta=%e > eta_max=%e\n",__func__,__LINE__,eta,pba->eta_table[pba->bt_size-1]);
    return _FAILURE_;
  }

  /** - deduce length of returned vector from format mode */ 

  if (return_format == long_info) {
    pvecback_size=pba->bg_size;
  }
  else { 
    pvecback_size=pba->bg_size_short;
  }

  /** - interpolate from pre-computed table with array_interpolate() or array_interpolate_growing_closeby() (depending on interpolation mode) */

  if (intermode == normal) {
    if (array_interpolate_spline(
				 pba->eta_table,
				 pba->bt_size,
				 pba->background_table,
				 pba->d2background_deta2_table,
				 pba->bg_size,
				 eta,
				 last_index,
				 pvecback_local,
				 pvecback_size,
				 Transmit_Error_Message) == _FAILURE_) {
      sprintf(pba->error_message,"%s(L:%d) : error in array_interpolate_spline() \n=>%s",__func__,__LINE__,Transmit_Error_Message);
      return _FAILURE_;
    }
  }
  if (intermode == closeby) {
    if (array_interpolate_spline_growing_closeby(
						 pba->eta_table,
						 pba->bt_size,
						 pba->background_table,
						 pba->d2background_deta2_table,
						 pba->bg_size,
						 eta,
						 last_index,
						 pvecback_local,
						 pvecback_size,
						 Transmit_Error_Message) == _FAILURE_) {
      sprintf(pba->error_message,"%s(L:%d) : error in array_interpolate_spline_growing_closeby() \n=>%s",__func__,__LINE__,Transmit_Error_Message);
      return _FAILURE_;  
    }
  }

  return _SUCCESS_;
}

/** 
  * Conformal time at given redhsift.
  *
  * Returns eta(z) by interpolation from pre-computed table.
  * This function can be called from whatever module at whatever time,
  * provided that background_init() has been called before, and 
  * background_free() has not been called yet.
  *
  * @param z Input: redshift
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

  /* scale factor */
  double a; 

  /* vector of background quantities */
  double * pvecback_local;

  int last_index; /* necessary for calling array_interpolate(), but never used */

  /** - check that \f$ z \f$ is in the pre-computed range */
  if (z < pba->z_table[pba->bt_size-1]) {
    sprintf(pba->error_message,"%s(L%d) : z=%e < z_min=%e\n",__func__,__LINE__,z,pba->z_table[pba->bt_size-1]);
    return _FAILURE_;
  }
  if (z > pba->z_table[0]) {
    sprintf(pba->error_message,"%s(L%d) : a=%e > a_max=%e\n",__func__,__LINE__,z,pba->z_table[0]);
    return _FAILURE_;
  }

  /** - interpolate from pre-computed table with array_interpolate() */
  if (array_interpolate_spline(
			       pba->z_table,
			       pba->bt_size, 
			       pba->eta_table,
			       pba->d2eta_dz2_table,
			       1,
			       z,
			       &last_index,
			       eta,
			       1,
			       Transmit_Error_Message) == _FAILURE_) {
    sprintf(pba->error_message,"%s(L:%d) : error in array_interpolate_spline() \n=>%s",__func__,__LINE__,Transmit_Error_Message);
    return _FAILURE_;
  }

  return _SUCCESS_;
}

  /**
   * Background quantities at given a.
   *
  * Function evaluating all background quantities which can be
  * computed analytically as a function of a. Hence, it does not
  * interpolate within the table, and it does not return the variables
  * which require a time integration. In particular, it cannot provide
  * the time as a function of a: for this purpose one should call
  * background_eta_of_z().  This function can be called at any time,
  * even before background_init().  It is called at each step of
  * integration by background_solve() and background_derivs().
  *
  * @param a Input: value of scale factor 
  * @param return_format Input: format of output vector 
  * @param pvecback_local Output: vector of background quantities (assmued to be already allocated) 
  * @return the error status
  */
int background_functions_of_a(
			      struct background *pba,
			      double a,
			      enum format_info return_format,
			      double * pvecback_local
			      ) {

  /** Summary: */

  /** - define local variables */

  /* total density */
  double rho_tot;
  /* total pressure */
  double p_tot;
  /* total relativistic density */
  double rho_r;
  /* scale factor relative to scale factor today */
  double a_rel;

  /** - initialize local variables */
  rho_tot = 0.;
  p_tot = 0.;
  rho_r=0.;
  a_rel = a / pba->a_today;

  if (a_rel <= 0.) { /* mainly to avoid segmentation faults */
    sprintf(pba->error_message,"%s(L:%d): a = %e instead of strictly positive \n",__func__,__LINE__,a_rel);
    return _FAILURE_;
  }

  /** - pass value of a to output */
  pvecback_local[pba->index_bg_a] = a;

  /** - compute each component's density and pressure */

  /* photons */
  pvecback_local[pba->index_bg_rho_g] = pba->Omega0_g * pow(pba->H0,2) / pow(a_rel,4);
  rho_tot += pvecback_local[pba->index_bg_rho_g];
  p_tot += (1./3.) * pvecback_local[pba->index_bg_rho_g];
  rho_r += pvecback_local[pba->index_bg_rho_g];

  /* baryons */
  pvecback_local[pba->index_bg_rho_b] = pba->Omega0_b * pow(pba->H0,2) / pow(a_rel,3);
  rho_tot += pvecback_local[pba->index_bg_rho_b];
  p_tot += 0;

  /* cdm */
  if (pba->has_cdm == _TRUE_) {
    pvecback_local[pba->index_bg_rho_cdm] = pba->Omega0_cdm * pow(pba->H0,2) / pow(a_rel,3);
    rho_tot += pvecback_local[pba->index_bg_rho_cdm];
    p_tot += 0.;
  }

  /* Lambda */
  if (pba->has_lambda == _TRUE_) {
    pvecback_local[pba->index_bg_rho_lambda] = pba->Omega0_lambda * pow(pba->H0,2);
    rho_tot += pvecback_local[pba->index_bg_rho_lambda];
    p_tot -= pvecback_local[pba->index_bg_rho_lambda];
  }

  /* dark energy fluid with constant w */
  if (pba->has_dark_energy_fluid == _TRUE_) {
    pvecback_local[pba->index_bg_rho_de] = pba->Omega0_de * pow(pba->H0,2) / pow(a_rel,3.*(1.+pba->w_de));
    rho_tot += pvecback_local[pba->index_bg_rho_de];
    p_tot += pba->w_de * pvecback_local[pba->index_bg_rho_de];
  }

  /* relativistic neutrinos (and all relativistic relics) */
  if (pba->has_nur == _TRUE_) {
    pvecback_local[pba->index_bg_rho_nur] = pba->Omega0_nur * pow(pba->H0,2) / pow(a_rel,4);
    rho_tot += pvecback_local[pba->index_bg_rho_nur];
    p_tot += (1./3.) * pvecback_local[pba->index_bg_rho_nur];
    rho_r += pvecback_local[pba->index_bg_rho_nur];
  }

  /** - compute relativistic density to total density ratio */
  pvecback_local[pba->index_bg_Omega_r] = rho_r / rho_tot;

  /** - compute expansion rate H from Friedmann equation */
  pvecback_local[pba->index_bg_H] = sqrt(rho_tot);

  /** - compute derivative of H with respect to conformal time */
  pvecback_local[pba->index_bg_H_prime] = - (3./2.) * (rho_tot + p_tot) * a;

  /** - compute other quantities in the exhaustive, redundent format: */
  if (return_format == long_info) {
    
    /** - compute critical density */
    pvecback_local[pba->index_bg_rho_crit] = rho_tot;
    if (pvecback_local[pba->index_bg_rho_crit] <= 0.) { /* to avoid segmentation fault in next lines */
      sprintf(pba->error_message,"%s(L:%d): rho_tot = %e instead of strictly positive \n",__func__,__LINE__,pvecback_local[pba->index_bg_rho_crit]);
      return _FAILURE_;
    }

    /** - compute Omega's */

    /* photons */
    pvecback_local[pba->index_bg_Omega_g] = pvecback_local[pba->index_bg_rho_g] / pvecback_local[pba->index_bg_rho_crit];
    /* baryons */
    pvecback_local[pba->index_bg_Omega_b] = pvecback_local[pba->index_bg_rho_b] / pvecback_local[pba->index_bg_rho_crit];
    /* cdm */
    if (pba->has_cdm == _TRUE_)
      pvecback_local[pba->index_bg_Omega_cdm] = pvecback_local[pba->index_bg_rho_cdm] / pvecback_local[pba->index_bg_rho_crit];
    /* Lambda */
    if (pba->has_lambda == _TRUE_)
      pvecback_local[pba->index_bg_Omega_lambda] = pvecback_local[pba->index_bg_rho_lambda] / pvecback_local[pba->index_bg_rho_crit];
    /* dark energy fluid with constant w */
    if (pba->has_dark_energy_fluid == _TRUE_)
      pvecback_local[pba->index_bg_Omega_de] = pvecback_local[pba->index_bg_rho_de] / pvecback_local[pba->index_bg_rho_crit];
    /* relativistic neutrinos */
    if (pba->has_nur == _TRUE_)
      pvecback_local[pba->index_bg_Omega_nur] = pvecback_local[pba->index_bg_rho_nur] / pvecback_local[pba->index_bg_rho_crit];
    
    /* one can put other variables here */
    /*  */
    /*  */

  }

  return _SUCCESS_;
  
}

/** 
  * Initialize the cosmo structure, including interpolation table.
  * 
  * Reads the cosmological parameters and initialize all fields in the
  * structure cosmo, in particular:
  *
  * - initializes all indices given the list of cosmological ingredients
  *
  * - integrates the background equations (Friedmann, ...) and initializes the background interpolation table.
  *
  * This function shall be called at the beginning of each run. It allocates memory spaces which should be freed later with background_free().
  *
  * @param ppr Input : Parameters describing how the computation is to be performed
  * @param pba Output : Initialized background structure
  * @return the error status
  */
int background_init(
  struct precision * ppr,
  struct background * pba
  ) {

  /** Summary: */

  /** - local variables : */

  double Omega0_tot;

  if (pba->background_verbose > 0)
    printf("Computing background\n");

  /** - assign values to all indices in vectors of background quantities with background_indices()*/
  if (background_indices(pba) == _FAILURE_) {
    sprintf(Transmit_Error_Message,"%s",pba->error_message);
    sprintf(pba->error_message,"%s(L:%d) : error in background_indices() \n=>%s",__func__,__LINE__,Transmit_Error_Message);
    return _FAILURE_;
  }

  /** - allocate memory for the two useful background vectors pvecback and pvecback_integration (global variables in background.c, used as long as background_init() is not finished) */
  pvecback = malloc(pba->bg_size*sizeof(double));
  if (pvecback==NULL) {
    sprintf(pba->error_message,"%s(L:%d): error : cannot allocate pvecback \n",__func__,__LINE__);
    return _FAILURE_;
  }

  /** - control that cosmological parameter values make sense */

  /* H0 in Mpc^{-1} */
  if ((pba->H0 < _H0_SMALL_)||(pba->H0 > _H0_BIG_)) {
    sprintf(pba->error_message,"%s(L:%d): H0 out of bounds (got %g, expected %g<H0<%g) \n",__func__,__LINE__,pba->H0,_H0_SMALL_,_H0_BIG_);
    return _FAILURE_;
  }

  pba->h = pba->H0 * _c_ /100.;

  /* curvature */
  Omega0_tot = pba->Omega0_g + pba->Omega0_b;
  if (pba->has_cdm == _TRUE_) {
    Omega0_tot += pba->Omega0_cdm;
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
  if (fabs(Omega0_tot-1.) > _TOLERANCE_ON_CURVATURE_) {
    sprintf(pba->error_message,"%s(L:%d): Non zero spatial curvature not available (Omega0 = %g) \n",__func__,__LINE__,Omega0_tot);
    return _FAILURE_;
  }

  /* other quantities which would lead to segmentation fault if zero */
  if (pba->a_today <= 0) {
    sprintf(pba->error_message,"%s(L:%d): input a_today = %e instead of strictly positive \n",__func__,__LINE__,pba->a_today);
    return _FAILURE_;
  }
  if (_Gyr_over_Mpc_ <= 0) {
    sprintf(pba->error_message,"%s(L:%d): _Gyr_over_Mpc = %e instead of strictly positive \n",__func__,__LINE__,_Gyr_over_Mpc_);
    return _FAILURE_;
  }

  /** - integrate background, allocate and fill the background table with background_solve()*/
  if (background_solve(ppr,pba) == _FAILURE_) {
    sprintf(Transmit_Error_Message,"%s",pba->error_message);
    sprintf(pba->error_message,"%s(L:%d) : error in background_solve() \n=>%s",__func__,__LINE__,Transmit_Error_Message);
    return _FAILURE_;
  }

  /** - free memory allocated for pvecback and pvecback_integration */
 
  free(pvecback);

  return _SUCCESS_;

}

 /**
 * Free all memory space allocated by background_init() for the
 * background table.
 * 
 * To be called at the end of each run, only when no further calls to
 * background_at_eta() or  background_eta_of_z() are needed. 
 *
 * @return the error status
 */
int background_free(
		    struct background *pba
		    ) {

  free(pba->eta_table);
  free(pba->z_table);
  free(pba->d2eta_dz2_table);
  free(pba->background_table);
  free(pba->d2background_deta2_table);

  return _SUCCESS_;
}

/**
 * Assign value to each relevant index in vectors of background quantities.
 *
 * Called once by background_init().
 *
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
  
  /* - put here additional ingredients that you want to appear in the short array */
  /*    */
  /*    */

  /* - end of indices in the short vector of background values */
  pba->bg_size_short = index_bg;

  /* - indices in the long version : */

  /* -> critical density */
  pba->index_bg_rho_crit = index_bg; 
  index_bg++;

  /* -> all Omega's : */

  /* photons */
  pba->index_bg_Omega_g = index_bg; 
  index_bg++;

  /* baryons */
  pba->index_bg_Omega_b = index_bg; 
  index_bg++;
  
  /* cdm */
  if (pba->has_cdm==_TRUE_) {
    pba->index_bg_Omega_cdm = index_bg; 
    index_bg++;
  }

  /* lambda */
  if (pba->has_lambda==_TRUE_) {
    pba->index_bg_Omega_lambda = index_bg; 
    index_bg++;
  }

  /* dark energy fluid */
  if (pba->has_dark_energy_fluid==_TRUE_) {
    pba->index_bg_Omega_de = index_bg; 
    index_bg++;
  }

  /* relativistic neutrinos and relics */
  if (pba->has_nur==_TRUE_) {
    pba->index_bg_Omega_nur = index_bg; 
    index_bg++;
  }
  
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
  if (pba->index_bi_eta != index_bi-1) {
    sprintf(pba->error_message,"%s(L%d) : background integration requires index_bi_eta to be the last of all index_bi's. \n",__func__,__LINE__);
    return _FAILURE_;
  }

  return _SUCCESS_;

}

/** 
 * Integrate background, allocate and fill the background interpolation table
 *
 * Called once by background_init().
 *
 * @return the error status
 */
int background_solve(
		     struct precision *ppr,
		     struct background *pba
		     ) {

  /** Summary: */

  /** - define local variables */

  /* contains all quantities relevant for the integration algorithm */
  struct generic_integrator_workspace gi;
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
  
  int last_index=0; /* necessary for calling array_interpolate(), but never used */

  /** - allocate vector of quantities to be integrated */
  pvecback_integration = malloc(pba->bi_size*sizeof(double));
  if (pvecback_integration==NULL) {
    sprintf(pba->error_message,"%s(L:%d): error : cannot allocate pvecback_integration \n",__func__,__LINE__);
    return _FAILURE_;
  }

  /** - initialize generic integrator with initialize_generic_integrator() */ 

  /* Size of vector to integrate is (pba->bi_size-1) rather than
   * (pba->bi_size), since eta is not integrated.
   */
  class_call(initialize_generic_integrator((pba->bi_size-1),&gi),
	     gi.error_message,
	     pba->error_message);

  /** - impose initial conditions with background_initial_conditions() */
  if (background_initial_conditions(ppr,pba,pvecback_integration) == _FAILURE_) {
    sprintf(Transmit_Error_Message,"%s",pba->error_message);
    sprintf(pba->error_message,"%s(L:%d) : probelm in background_initial_conditions() \n=>%s",__func__,__LINE__,Transmit_Error_Message);
    return _FAILURE_;
  }

  /* here eta_end is in fact the initial time (in the next loop
     eta_start = eta_end) */
  eta_end=pvecback_integration[pba->index_bi_eta];

  /** - create a growTable with gt_init() */
  if (gt_init(&gTable) == _FAILURE_) {
    sprintf(pba->error_message,"%s(L:%d) : error in gt_init() \n=>%s",__func__,__LINE__,gTable.error_message);
  }
  
  /* initialize the counter for the number of steps */
  pba->bt_size=0;

  /** - loop over integration steps : call background_functions_of_a(), find step size, save data in growTable with gt_add(), perform one step with generic_integrator(), store new value of eta */

  while (pvecback_integration[pba->index_bi_a] < pba->a_today) {

    eta_start = eta_end;

    /* -> find step size (trying to adjust the last step as close as possible to the one needed to reach a=a_today; need not be exact, difference corrected later) */
    if (background_functions_of_a(pba,pvecback_integration[pba->index_bi_a], short_info, pvecback) == _FAILURE_) {
      sprintf(Transmit_Error_Message,"%s",pba->error_message);
      sprintf(pba->error_message,"%s(L:%d) : error in background_functions_of_a() \n=>%s",__func__,__LINE__,Transmit_Error_Message);
      return _FAILURE_;
    }

    if ((pvecback_integration[pba->index_bi_a]*(1.+ppr->back_integration_stepsize)) < pba->a_today) {
      eta_end = eta_start + ppr->back_integration_stepsize / (pvecback_integration[pba->index_bi_a]*pvecback[pba->index_bg_H]); 
      /* no possible segmentation fault here: non-zeroness of "a" has been checked in background_functions_of_a() */
    }
    else {
      eta_end = eta_start + (1./pvecback_integration[pba->index_bi_a]-1.) / (pvecback_integration[pba->index_bi_a]*pvecback[pba->index_bg_H]);  
      /* no possible segmentation fault here: non-zeroness of "a" has been checked in background_functions_of_a() */
    }

    if ((eta_end-eta_start) < ppr->smallest_allowed_variation) {
      sprintf(pba->error_message,"%s(L:%d) : error : integration step =%e < machine precision : leads either to numerical error or infinite loop \n",__func__,__LINE__,eta_end-eta_start,gTable.error_message);
      return _FAILURE_;
    }

    /* -> save data in growTable */
    if (gt_add(&gTable,_GT_END_,(void *) pvecback_integration,sizeof(double)*pba->bi_size) == _FAILURE_) { 
      sprintf(pba->error_message,"%s(L:%d) : error in gt_add() \n=>%s",__func__,__LINE__,gTable.error_message);
      return _FAILURE_;
    }
    pba->bt_size++;

    /* -> perform one step */
    class_call(generic_integrator(background_derivs,
				  eta_start,
				  eta_end,
				  pvecback_integration,
				  pba,
				  ppr->tol_background_integration,
				  ppr->smallest_allowed_variation,
				  &gi),
	       gi.error_message,
	       pba->error_message);
    
    /* -> store value of eta */
    pvecback_integration[pba->index_bi_eta]=eta_end;

  }

  /** - save last data in growTable with gt_add() */
  if (gt_add(&gTable,_GT_END_,(void *) pvecback_integration,sizeof(double)*pba->bi_size) == _FAILURE_) {
    sprintf(pba->error_message,"%s(L:%d) : error in gt_add() \n=>%s",__func__,__LINE__,gTable.error_message);
    return _FAILURE_;
  }
  pba->bt_size++;


  /* integration finished */

  /** - clean up generic integrator with cleanup_generic_integrator() */
  class_call(cleanup_generic_integrator(&gi),
	     gi.error_message,
	     pba->error_message);

  /** - retrieve data stored in the growTable with gt_getPtr() */
  if (gt_getPtr(&gTable,(void**)&pData) == _FAILURE_) {
    sprintf(pba->error_message,"%s(L:%d) : error in gt_getPtr() \n=>%s",__func__,__LINE__,gTable.error_message);
    return _FAILURE_;
  }

  /** - interpolate to get quantities precisely today with array_interpolate() */
  if (array_interpolate(
		   pData,
		   pba->bi_size,
		   pba->bt_size,
		   pba->index_bi_a, 
		   pba->a_today,
		   &last_index,
		   pvecback_integration,
		   pba->bi_size,
		   Transmit_Error_Message) == _FAILURE_) {
    sprintf(pba->error_message,"%s(L:%d) : error in array_interpolate() \n=>%s",__func__,__LINE__,Transmit_Error_Message);
    return _FAILURE_;
  }

  /* substitute last line with quantities today */
  for (i=0; i<pba->bi_size; i++) 
    pData[(pba->bt_size-1)*pba->bi_size+i]=pvecback_integration[i];

  /** - deduce age of the Universe */
  /* -> age in Gyears */
  pba->age = pvecback_integration[pba->index_bi_time]/_Gyr_over_Mpc_;
  /* -> conformal age in Mpc */
  pba->conformal_age = pvecback_integration[pba->index_bi_eta];

  /** - allocate background tables */
  pba->eta_table = malloc(pba->bt_size * sizeof(double));
  if (pba->eta_table==NULL) {
    sprintf(pba->error_message,"%s(L:%d): error : cannot allocate pba->eta_table \n",__func__,__LINE__);
    return _FAILURE_;
  }
  pba->z_table = malloc(pba->bt_size * sizeof(double));
  if (pba->z_table==NULL) {
    sprintf(pba->error_message,"%s(L:%d): error : cannot allocate pba->z_table \n",__func__,__LINE__);
    return _FAILURE_;
  }
  pba->d2eta_dz2_table = malloc(pba->bt_size * sizeof(double));
  if (pba->d2eta_dz2_table==NULL) {
    sprintf(pba->error_message,"%s(L:%d): error : cannot allocate pba->d2eta_dz2_table \n",__func__,__LINE__);
    return _FAILURE_;
  }
  pba->background_table = malloc(pba->bt_size * pba->bg_size * sizeof(double));
  if (pba->background_table==NULL) {
    sprintf(pba->error_message,"%s(L:%d): error : cannot allocate pba->background_table \n",__func__,__LINE__);
    return _FAILURE_;
  }
  pba->d2background_deta2_table = malloc(pba->bt_size * pba->bg_size * sizeof(double));
  if (pba->d2background_deta2_table==NULL) {
    sprintf(pba->error_message,"%s(L:%d): error : cannot allocate pba->d2background_deta2_table \n",__func__,__LINE__);
    return _FAILURE_;
  }
  
  /** - In a loop over lines, fill background table using the result of the integration plus background_functions_of_a() */
  for (i=0; i < pba->bt_size; i++) {
    
    /* -> establish correspondance between the integrated variable and the bg variables */

    pba->eta_table[i] = pData[i*pba->bi_size+pba->index_bi_eta];

    if (pData[i*pba->bi_size+pba->index_bi_a] > 0) /* mainly to avoid segmentation fault */
      pba->z_table[i] = pba->a_today/pData[i*pba->bi_size+pba->index_bi_a]-1.;
    else {
      sprintf(pba->error_message,"%s(L:%d): a = %e instead of strictly positive \n",__func__,__LINE__,pData[i*pba->bi_size+pba->index_bi_a]);
      return _FAILURE_;
    }

    pvecback[pba->index_bg_a] = pData[i*pba->bi_size+pba->index_bi_a];
    pvecback[pba->index_bg_time] = pData[i*pba->bi_size+pba->index_bi_time];
    pvecback[pba->index_bg_conf_distance] = pba->conformal_age - pData[i*pba->bi_size+pba->index_bi_eta];
    pvecback[pba->index_bg_rs] = pData[i*pba->bi_size+pba->index_bi_rs];

    /* -> compute all other quantities */
    if (background_functions_of_a(pba,pvecback[pba->index_bg_a], long_info, pvecback) == _FAILURE_) {
      sprintf(Transmit_Error_Message,"%s",pba->error_message);
      sprintf(pba->error_message,"%s(L:%d) : error in background_functions_of_a() \n=>%s",__func__,__LINE__,Transmit_Error_Message);
      return _FAILURE_;
    }
    
    /* -> write in the table */
    memcopy_result = memcpy(pba->background_table + i*pba->bg_size,pvecback,pba->bg_size*sizeof(double));
    if (memcopy_result != pba->background_table + i*pba->bg_size) {
      sprintf(pba->error_message,"%s(L:%d) : error : cannot copy data back to pba->background_table \n",__func__,__LINE__);
      return _FAILURE_;
                                                                                                                                                                                                                   }
  }

  /** - free the growTable with gt_free() */

  gt_free(&gTable);
  
  /** - fill tables of second derivatives (in view of spline interpolation) */
  if (array_spline_table_lines(pba->z_table,
			       pba->bt_size,
			       pba->eta_table,
			       1,
			       pba->d2eta_dz2_table,
			       _SPLINE_EST_DERIV_,
			       Transmit_Error_Message) == _FAILURE_) {
    sprintf(pba->error_message,"%s(L:%d) : error in array_spline_table_lines() \n=>%s",__func__,__LINE__,Transmit_Error_Message);
    return _FAILURE_;
  }

  if (array_spline_table_lines(pba->eta_table,
			       pba->bt_size,
			       pba->background_table,
			       pba->bg_size,
			       pba->d2background_deta2_table,
			       _SPLINE_EST_DERIV_,
			       Transmit_Error_Message) == _FAILURE_) {
    sprintf(pba->error_message,"%s(L:%d) : error in array_spline_table_lines \n=>%s",__func__,__LINE__,Transmit_Error_Message);
    return _FAILURE_;
  }

  if (pba->background_verbose > 0) {
    printf(" -> age = %f Gyr\n",pba->age);
    printf(" -> conformal age = %f Mpc\n",pba->conformal_age);    
  }

  free(pvecback_integration);

  return _SUCCESS_;
  
}

/**
 * Assign initial values to background integrated variables (initial a, eta, t, etc).
 *
 * Called once by background_solve().
 *
 * @return the error status
 */
int background_initial_conditions(
				  struct precision *ppr,
				  struct background *pba,
				  double * pvecback_integration
				  ) {

  /** Summary: */

  /** - define local variables */

  /* scale factor */
  double a;

  /** - fix initial value of \f$ a \f$ */
  a = ppr->a_ini_over_a_today_default * pba->a_today;
  pvecback_integration[pba->index_bi_a] = a;

  /* for some models, we will need to add here some tests on the
     validity of this choice (e.g.: are massive neutrinos
     relativistic? etc.) If the test is OK, fix other initial values: */

  /** - compute initial H with background_functions_of_a() */
  if (background_functions_of_a(pba,a, short_info, pvecback) == _FAILURE_) {
    sprintf(Transmit_Error_Message,"%s",pba->error_message);
    sprintf(pba->error_message,"%s(L:%d) : error in background_functions_of_a() \n=>%s",__func__,__LINE__,Transmit_Error_Message);
    return _FAILURE_;
  }

  /** - compute initial proper time, assuming radiation-dominated
        universe since Big Bang and therefore \f$ t=1/(2H) \f$ (good
        approximation for most purposes) */
  if (pvecback[pba->index_bg_H] > 0.) /* mainly to avoid segmentation fault */
    pvecback_integration[pba->index_bi_time] = 1./(2.*  pvecback[pba->index_bg_H]);
  else { 
    sprintf(pba->error_message,"%s(L:%d): H = %e instead of strictly positive \n",__func__,__LINE__,pvecback[pba->index_bg_H]);
    return _FAILURE_;
  }

  /** - compute initial conformal time, assuming radiation-dominated
        universe since Big Bang and therefore \f$ \eta=1/(aH) \f$
        (good approximation for most purposes) */
  pvecback_integration[pba->index_bi_eta] = 1./(a * pvecback[pba->index_bg_H]);

  /** - compute initial sound horizon, assuming c_s=1/sqrt(3) initially */
  pvecback_integration[pba->index_bi_rs] = pvecback_integration[pba->index_bi_eta]/sqrt(3.);

  return _SUCCESS_;

}

/** 
  * Subroutine evaluating the derivative with respect to conformal time of quantities which are integrated (a, t, etc). 
  *
  * Passed as an argument to the generic_integrator() function.
  *
  * @param eta Input : conformal time
  * @param y Input : vector of variable
  * @param dy Output : its derivative (already allocated)
  * @param fixed_parameters Input: pointer to fixed parameters (e.g. indices); here, this is just a pointer to the background structure; passed as a generic pointer in order to match declarations in generic integrator module.
  */
int background_derivs(
		      double eta,
		      double* y, 
		      double* dy,
		      void * fixed_parameters
		      ) {

  struct background * pba;

  pba =  fixed_parameters;

  /** - Calculates functions of /f$ a /f$ with background_functions_of_a() */
  if (background_functions_of_a((struct background *)pba,y[pba->index_bi_a], short_info, pvecback) == _FAILURE_) {
    sprintf(Transmit_Error_Message,"%s",pba->error_message);
    sprintf(pba->error_message,"%s(L:%d) : error in calling background_derivs() \n=>%s",__func__,__LINE__,Transmit_Error_Message);
      return;
    }

    /** - calculate /f$ a'=a^2 H /f$ */
    dy[pba->index_bi_a] = y[pba->index_bi_a] * y[pba->index_bi_a] * pvecback[pba->index_bg_H];

    /** - calculate /f$ t' = a /f$ */
    dy[pba->index_bi_time] = y[pba->index_bi_a];

    /** - calculate /f$ r_s' = c_s /f$ */
    if (pvecback[pba->index_bg_rho_g] > 0.) /* mainly to avoid segmentation fault */
      dy[pba->index_bi_rs] = 1./sqrt(3.*(1.+3.*pvecback[pba->index_bg_rho_b]/4./pvecback[pba->index_bg_rho_g]));
    else {
      sprintf(pba->error_message,"%s(L:%d): rho_g = %e instead of strictly positive \n",__func__,__LINE__,pvecback[pba->index_bg_rho_g]);
      return;
    }

    return;

}
