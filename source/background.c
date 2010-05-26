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
  * @param pvecback Output: vector (assumed to be already allocated)
  * @return the error status
  */
int background_at_eta(
		      struct background *pba,
		      double eta,
		      enum format_info return_format,
		      enum interpolation_mode intermode,
		      int * last_index,
		      double * pvecback
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

  if (return_format == long_info) {
    pvecback_size=pba->bg_size;
  }
  else { 
    pvecback_size=pba->bg_size_short;
  }

  /** - interpolate from pre-computed table with array_interpolate() or array_interpolate_growing_closeby() (depending on interpolation mode) */

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

  int last_index; /* necessary for calling array_interpolate(), but never used */

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
  * @param pvecback Output: vector of background quantities (assmued to be already allocated) 
  * @return the error status
  */
int background_functions_of_a(
			      struct background *pba,
			      double a,
			      enum format_info return_format,
			      double * pvecback
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

  /* cdm */
  if (pba->has_cdm == _TRUE_) {
    pvecback[pba->index_bg_rho_cdm] = pba->Omega0_cdm * pow(pba->H0,2) / pow(a_rel,3);
    rho_tot += pvecback[pba->index_bg_rho_cdm];
    p_tot += 0.;
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

  /** - compute relativistic density to total density ratio */
  pvecback[pba->index_bg_Omega_r] = rho_r / rho_tot;

  /** - compute expansion rate H from Friedmann equation */
  pvecback[pba->index_bg_H] = sqrt(rho_tot);

  /** - compute derivative of H with respect to conformal time */
  pvecback[pba->index_bg_H_prime] = - (3./2.) * (rho_tot + p_tot) * a;

  /** - compute other quantities in the exhaustive, redundent format: */
  if (return_format == long_info) {
    
    /** - compute critical density */
    pvecback[pba->index_bg_rho_crit] = rho_tot;
    class_test(pvecback[pba->index_bg_rho_crit] <= 0.,
	       pba->error_message,
	       "rho_tot = %e instead of strictly positive",pvecback[pba->index_bg_rho_crit]);

    /** - compute Omega's */

    /* photons */
    pvecback[pba->index_bg_Omega_g] = pvecback[pba->index_bg_rho_g] / pvecback[pba->index_bg_rho_crit];
    /* baryons */
    pvecback[pba->index_bg_Omega_b] = pvecback[pba->index_bg_rho_b] / pvecback[pba->index_bg_rho_crit];
    /* cdm */
    if (pba->has_cdm == _TRUE_)
      pvecback[pba->index_bg_Omega_cdm] = pvecback[pba->index_bg_rho_cdm] / pvecback[pba->index_bg_rho_crit];
    /* Lambda */
    if (pba->has_lambda == _TRUE_)
      pvecback[pba->index_bg_Omega_lambda] = pvecback[pba->index_bg_rho_lambda] / pvecback[pba->index_bg_rho_crit];
    /* dark energy fluid with constant w */
    if (pba->has_dark_energy_fluid == _TRUE_)
      pvecback[pba->index_bg_Omega_de] = pvecback[pba->index_bg_rho_de] / pvecback[pba->index_bg_rho_crit];
    /* relativistic neutrinos */
    if (pba->has_nur == _TRUE_)
      pvecback[pba->index_bg_Omega_nur] = pvecback[pba->index_bg_rho_nur] / pvecback[pba->index_bg_rho_crit];
    
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
  class_call(background_indices(pba),
	     pba->error_message,
	     pba->error_message);
  
  /** - control that cosmological parameter values make sense */

  /* H0 in Mpc^{-1} */
  class_test((pba->H0 < _H0_SMALL_)||(pba->H0 > _H0_BIG_),
    pba->error_message,
    "H0=%g out of bounds (%g<H0<%g) \n",pba->H0,_H0_SMALL_,_H0_BIG_);

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
  class_test(pba->index_bi_eta != index_bi-1,
	     pba->error_message,
	     "background integration requires index_bi_eta to be the last of all index_bi's");

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

  /** - loop over integration steps : call background_functions_of_a(), find step size, save data in growTable with gt_add(), perform one step with generic_integrator(), store new value of eta */

  while (pvecback_integration[pba->index_bi_a] < pba->a_today) {

    eta_start = eta_end;

    /* -> find step size (trying to adjust the last step as close as possible to the one needed to reach a=a_today; need not be exact, difference corrected later) */
    class_call(background_functions_of_a(pba,pvecback_integration[pba->index_bi_a], short_info, pvecback),
	       pba->error_message,
	       pba->error_message);

    if ((pvecback_integration[pba->index_bi_a]*(1.+ppr->back_integration_stepsize)) < pba->a_today) {
      eta_end = eta_start + ppr->back_integration_stepsize / (pvecback_integration[pba->index_bi_a]*pvecback[pba->index_bg_H]); 
      /* no possible segmentation fault here: non-zeroness of "a" has been checked in background_functions_of_a() */
    }
    else {
      eta_end = eta_start + (1./pvecback_integration[pba->index_bi_a]-1.) / (pvecback_integration[pba->index_bi_a]*pvecback[pba->index_bg_H]);  
      /* no possible segmentation fault here: non-zeroness of "a" has been checked in background_functions_of_a() */
    }

    class_test((eta_end-eta_start) < ppr->smallest_allowed_variation,
	       pba->error_message,
	       "integration step =%e < machine precision : leads either to numerical error or infinite loop",eta_end-eta_start);

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
  
  /** - In a loop over lines, fill background table using the result of the integration plus background_functions_of_a() */
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
    class_call(background_functions_of_a(pba,pData[i*pba->bi_size+pba->index_bi_a], long_info, pvecback),
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
 * Assign initial values to background integrated variables (initial a, eta, t, etc).
 *
 * Called once by background_solve().
 *
 * @return the error status
 */
int background_initial_conditions(
				  struct precision *ppr,
				  struct background *pba,
				  double * pvecback,
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
  class_call(background_functions_of_a(pba,a, short_info, pvecback),
	     pba->error_message,
	     pba->error_message);

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
  * Subroutine evaluating the derivative with respect to conformal time of quantities which are integrated (a, t, etc). 
  *
  * This is one of the few functions in the code which are passed to the generic_integrator() routine. 
  * Since generic_integrator() should work with functions passed from various modules, the format of the arguments
  * is a bit special:
  * - fixed input parameters and wokspaces are passed through a generic pointer. Here, this is just a pointer to the 
  *   background structure and to a background vector, but generic_integrator() doesn't know its fine structure.
  * - the error management is a bit special: errors are not written as usual to pba->error_message, but to a generic 
  *   error_message passed in the list of arguments.
  *
  * @param eta Input : conformal time
  * @param y Input : vector of variable
  * @param dy Output : its derivative (already allocated)
  * @param parameters_and_workspace Input: pointer to fixed parameters (e.g. indices)
  * @param error_message Output : error message
  */
int background_derivs(
		      double eta,
		      double* y, 
		      double* dy,
		      void * parameters_and_workspace,
		      ErrorMsg error_message
		      ) {

  struct background_parameters_and_workspace * pbpaw;
  struct background * pba;
  double * pvecback;

  pbpaw = parameters_and_workspace;
  pba =  pbpaw->pba;
  pvecback = pbpaw->pvecback;

  /** - Calculates functions of /f$ a /f$ with background_functions_of_a() */
  class_call(background_functions_of_a((struct background *)pba,y[pba->index_bi_a], short_info, pvecback),
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
