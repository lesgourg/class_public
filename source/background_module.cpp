/** @file background.c Documented background module
 *
 * * Julien Lesgourgues, 17.04.2011
 * * routines related to ncdm written by T. Tram in 2011
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
 *    correspondence between redshift and conformal time.
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
 * 2. in turn, quantities {B} can be found as a function of conformal
 * time by integrating the background equations.
 *
 * 3. some other quantities that we will call {C} (like e.g. the
 * sound horizon or proper time) also require an integration with
 * respect to time, that cannot be inferred analytically from
 * parameters {B}.
 *
 * So, we define the following routines:
 *
 * - background_functions() returns all background
 *    quantities {A} as a function of quantities {B}.
 *
 * - background_solve() integrates the quantities {B} and {C} with
 *    respect to conformal time; this integration requires many calls
 *    to background_functions().
 *
 * - the result is stored in the form of a big table in the background
 *    structure. There is one column for conformal time 'tau'; one or
 *    more for quantities {B}; then several columns for quantities {A}
 *    and {C}.
 *
 * Later in the code, if we know the variables {B} and need some
 * quantity {A}, the quickest and most precise way is to call directly
 * background_functions() (for instance, in simple models, if we want
 * H at a given value of the scale factor). If we know 'tau' and want
 * any other quantity, we can call background_at_tau(), which
 * interpolates in the table and returns all values. Finally it can be
 * useful to get 'tau' for a given redshift 'z': this can be done with
 * background_tau_of_z(). So if we are somewhere in the code, knowing
 * z and willing to get background quantities, we should call first
 * background_tau_of_z() and then background_at_tau().
 *
 *
 * In order to save time, background_at_tau() can be called in three
 * modes: short_info, normal_info, long_info (returning only essential
 * quantities, or useful quantities, or rarely useful
 * quantities). Each line in the interpolation table is a vector whose
 * first few elements correspond to the short_info format; a larger
 * fraction contribute to the normal format; and the full vector
 * corresponds to the long format. The guideline is that short_info
 * returns only geometric quantities like a, H, H'; normal format
 * returns quantities strictly needed at each step in the integration
 * of perturbations; long_info returns quantities needed only
 * occasionally.
 *
 * In summary, the following functions can be called from other modules:
 *
 * -# background_init() at the beginning
 * -# background_at_tau(), background_tau_of_z() at any later time
 * -# background_free() at the end, when no more calls to the previous functions are needed
 */

#include "background_module.h"
#include "non_cold_dark_matter.h"

BackgroundModule::BackgroundModule(InputModulePtr input_module)
: BaseModule(input_module) {
  if (background_init() != _SUCCESS_) {
    throw std::runtime_error(error_message_);
  }
}

BackgroundModule::~BackgroundModule() {
  background_free();
}

// Wrapper functions to pass non-static member functions
int BackgroundModule::background_derivs(double z, double* y, double* dy, void* parameters_and_workspace, ErrorMsg error_message) {
  auto pbpaw = static_cast<background_parameters_and_workspace*>(parameters_and_workspace);
  return pbpaw->background_module->background_derivs_member(z, y, dy, parameters_and_workspace, error_message);
}
int BackgroundModule::background_derivs_loga(double loga, double* y, double* dy, void* parameters_and_workspace, ErrorMsg error_message) {
  auto pbpaw = static_cast<background_parameters_and_workspace*>(parameters_and_workspace);
  return pbpaw->background_module->background_derivs_loga_member(loga, y, dy, parameters_and_workspace, error_message);
}
int BackgroundModule::background_add_line_to_bg_table(double loga, double* y, double* dy, int index_loga, void* parameters_and_workspace, ErrorMsg error_message) {
  auto pbpaw = static_cast<background_parameters_and_workspace*>(parameters_and_workspace);
  return pbpaw->background_module->background_add_line_to_bg_table_member(loga, y, dy, index_loga, parameters_and_workspace, error_message);
}


/**
 * Background quantities at given conformal time tau.
 *
 * Evaluates all background quantities at a given value of
 * conformal time by reading the pre-computed table and interpolating.
 *
 * @param pba           Input: pointer to background structure (containing pre-computed table)
 * @param tau           Input: value of conformal time
 * @param return_format Input: format of output vector (short, normal, long)
 * @param intermode     Input: interpolation mode (normal or closeby)
 * @param last_index    Input/Output: index of the previous/current point in the interpolation array (input only for closeby mode, output for both)
 * @param pvecback      Output: vector (assumed to be already allocated)
 * @return the error status
 */

int BackgroundModule::background_at_tau(double tau,
                                        short return_format,
                                        short intermode,
                                        int* last_index,
                                        double* pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size compatible with return_format) */
                                        ) const {

  /** Summary: */

  /** - define local variables */

  /* size of output vector, controlled by input parameter return_format */
  int pvecback_size;

  /** - check that tau is in the pre-computed range */

  class_test(tau < tau_table_[0],
             error_message_,
             "out of range: tau=%e < tau_min=%e, you should decrease the precision parameter a_ini_over_a_today_default\n",tau,tau_table_[0]);

  class_test(tau > tau_table_[bt_size_-1],
             error_message_,
             "out of range: tau=%e > tau_max=%e\n",tau,tau_table_[bt_size_-1]);

  /** - deduce length of returned vector from format mode */

  if (return_format == pba->normal_info) {
    pvecback_size=bg_size_normal_;
  }
  else {
    if (return_format == pba->short_info) {
      pvecback_size=bg_size_short_;
    }
    else {
      pvecback_size=bg_size_;
    }
  }

  /** - interpolate from pre-computed table with array_interpolate()
      or array_interpolate_growing_closeby() (depending on
      interpolation mode) */

  if (intermode == pba->inter_normal) {
    class_call(array_interpolate_spline(
                                        tau_table_,
                                        bt_size_,
                                        background_table_,
                                        d2background_dtau2_table_,
                                        bg_size_,
                                        tau,
                                        last_index,
                                        pvecback,
                                        pvecback_size,
                                        error_message_),
               error_message_,
               error_message_);
  }
  if (intermode == pba->inter_closeby) {
    class_call(array_interpolate_spline_growing_closeby(
                                                        tau_table_,
                                                        bt_size_,
                                                        background_table_,
                                                        d2background_dtau2_table_,
                                                        bg_size_,
                                                        tau,
                                                        last_index,
                                                        pvecback,
                                                        pvecback_size,
                                                        error_message_),
               error_message_,
               error_message_);
  }

  return _SUCCESS_;
}

/**
 * Conformal time at given redshift.
 *
 * Returns tau(z) by interpolation from pre-computed table.
 *
 * @param pba Input: pointer to background structure
 * @param z   Input: redshift
 * @param tau Output: conformal time
 * @return the error status
 */

int BackgroundModule::background_tau_of_z(double z, double* tau) const {

  /** Summary: */

  /** - define local variables */

  /* necessary for calling array_interpolate(), but never used */
  int last_index;

  /** - check that \f$ z \f$ is in the pre-computed range */
  class_test(z < z_table_[bt_size_ - 1],
             error_message_,
             "out of range: z=%e < z_min=%e\n", z, z_table_[bt_size_ - 1]);

  class_test(z > z_table_[0],
             error_message_,
             "out of range: a=%e > a_max=%e\n", z, z_table_[0]);

  /** - interpolate from pre-computed table with array_interpolate() */
  class_call(array_interpolate_spline(
                                      z_table_,
                                      bt_size_,
                                      tau_table_,
                                      d2tau_dz2_table_,
                                      1,
                                      z,
                                      &last_index,
                                      tau,
                                      1,
                                      error_message_),
             error_message_,
             error_message_);

  return _SUCCESS_;
}

/**
 * Background quantities at given \f$ a \f$.
 *
 * Function evaluating all background quantities which can be computed
 * analytically as a function of {B} parameters such as the scale factor 'a'
 * (see discussion at the beginning of this file). In extended
 * cosmological models, the pvecback_B vector contains other input parameters than
 * just 'a', e.g. (phi, phidot) for quintessence, some temperature of
 * exotic relics, etc...
 *
 * @param pba           Input: pointer to background structure
 * @param pvecback_B    Input: vector containing all {B} type quantities (scale factor, ...)
 * @param return_format Input: format of output vector
 * @param pvecback      Output: vector of background quantities (assumed to be already allocated)
 * @return the error status
 */

int BackgroundModule::background_functions(double* pvecback_B, /* Vector containing all {B} quantities. */
                                           short return_format,
                                           double* pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size compatible with return_format) */
                                           ) {

  /** Summary: */

  /** - define local variables */

  /* total density */
  double rho_tot;
  /* critical density */
  double rho_crit;
  /* total pressure */
  double p_tot;
  /* total relativistic density */
  double rho_r;
  /* total non-relativistic density */
  double rho_m;
  /* scale factor relative to scale factor today */
  double a_rel;
  /* background ncdm quantities */
  double rho_ncdm,p_ncdm,pseudo_p_ncdm;
  /* index for n_ncdm species */
  int n_ncdm;
  /* fluid's time-dependent equation of state parameter */
  double w_fld, dw_over_da, integral_fld;
  /* scale factor */
  double a;
  /* scalar field quantities */
  double phi, phi_prime;
  /* Since we only know a_prime_over_a after we have rho_tot,
     it is not possible to simply sum up p_tot_prime directly.
     Instead we sum up dp_dloga = p_prime/a_prime_over_a. The formula is
     p_prime = a_prime_over_a * dp_dloga = a_prime_over_a * Sum [ (w_prime/a_prime_over_a -3(1+w)w)rho].
     Note: The scalar field contribution must be added in the end, as an exception!*/
  double dp_dloga;

  /** - initialize local variables */
  a = pvecback_B[index_bi_a_];
  rho_tot = 0.;
  p_tot = 0.;
  dp_dloga = 0.;
  rho_r=0.;
  rho_m=0.;
  a_rel = a / pba->a_today;

  class_test(a_rel <= 0.,
             error_message_,
             "a = %e instead of strictly positive",a_rel);

  /** - pass value of \f$ a\f$ to output */
  pvecback[index_bg_a_] = a;

  /** - compute each component's density and pressure */

  /* photons */
  pvecback[index_bg_rho_g_] = pba->Omega0_g*pow(pba->H0, 2)/pow(a_rel, 4);
  rho_tot += pvecback[index_bg_rho_g_];
  p_tot += 1./3.*pvecback[index_bg_rho_g_];
  dp_dloga += -4./3.*pvecback[index_bg_rho_g_];
  rho_r += pvecback[index_bg_rho_g_];

  /* baryons */
  pvecback[index_bg_rho_b_] = pba->Omega0_b*pow(pba->H0, 2)/pow(a_rel, 3);
  rho_tot += pvecback[index_bg_rho_b_];
  p_tot += 0;
  rho_m += pvecback[index_bg_rho_b_];

  /* cdm */
  if (pba->has_cdm == _TRUE_) {
    pvecback[index_bg_rho_cdm_] = pba->Omega0_cdm*pow(pba->H0, 2)/pow(a_rel, 3);
    rho_tot += pvecback[index_bg_rho_cdm_];
    p_tot += 0.;
    rho_m += pvecback[index_bg_rho_cdm_];
  }

  /* dcdm */
  if (pba->has_dcdm == _TRUE_) {
    /* Pass value of rho_dcdm to output */
    pvecback[index_bg_rho_dcdm_] = pvecback_B[index_bi_rho_dcdm_];
    rho_tot += pvecback[index_bg_rho_dcdm_];
    p_tot += 0.;
    rho_m += pvecback[index_bg_rho_dcdm_];
  }

  /* dr */
  if (pba->has_dr == _TRUE_) {
    /* Pass value of rho_dr to output */
    pvecback[index_bg_rho_dr_] = pvecback_B[index_bi_rho_dr_];
    rho_tot += pvecback[index_bg_rho_dr_];
    p_tot += 1./3.*pvecback[index_bg_rho_dr_];
    dp_dloga += -4./3.*pvecback[index_bg_rho_dr_];
    rho_r += pvecback[index_bg_rho_dr_];
  }

  /* Scalar field */
  if (pba->has_scf == _TRUE_) {
    phi = pvecback_B[index_bi_phi_scf_];
    phi_prime = pvecback_B[index_bi_phi_prime_scf_];
    pvecback[index_bg_phi_scf_] = phi; // value of the scalar field phi
    pvecback[index_bg_phi_prime_scf_] = phi_prime; // value of the scalar field phi derivative wrt conformal time
    pvecback[index_bg_V_scf_] = V_scf(phi); //V_scf(phi); //write here potential as function of phi
    pvecback[index_bg_dV_scf_] = dV_scf(phi); // dV_scf(phi); //potential' as function of phi
    pvecback[index_bg_ddV_scf_] = ddV_scf(phi); // ddV_scf(phi); //potential'' as function of phi
    pvecback[index_bg_rho_scf_] = (phi_prime*phi_prime/(2*a*a) + V_scf(phi))/3.; // energy of the scalar field. The field units are set automatically by setting the initial conditions
    pvecback[index_bg_p_scf_] =(phi_prime*phi_prime/(2*a*a) - V_scf(phi))/3.; // pressure of the scalar field
    rho_tot += pvecback[index_bg_rho_scf_];
    p_tot += pvecback[index_bg_p_scf_];
    dp_dloga += 0.0; /** <-- This depends on a_prime_over_a, so we cannot add it now! */
    //divide relativistic & nonrelativistic (not very meaningful for oscillatory models)
    rho_r += 3.*pvecback[index_bg_p_scf_]; //field pressure contributes radiation
    rho_m += pvecback[index_bg_rho_scf_] - 3.*pvecback[index_bg_p_scf_]; //the rest contributes matter
    //printf(" a= %e, Omega_scf = %f, \n ",a_rel, pvecback[index_bg_rho_scf_]/rho_tot );
  }

  /* ncdm */
  if (pba->has_ncdm == _TRUE_) {

    /* Loop over species: */
    for(n_ncdm=0; n_ncdm<pba->N_ncdm; n_ncdm++){

      /* function returning background ncdm[n_ncdm] quantities (only
         those for which non-NULL pointers are passed) */
      class_call(pba->ncdm->background_ncdm_momenta(n_ncdm,
                                                    1./a_rel - 1.,
                                                    NULL,
                                                    &rho_ncdm,
                                                    &p_ncdm,
                                                    NULL,
                                                    &pseudo_p_ncdm),
                 error_message_,
                 error_message_);

      pvecback[index_bg_rho_ncdm1_ + n_ncdm] = rho_ncdm;
      rho_tot += rho_ncdm;
      pvecback[index_bg_p_ncdm1_ + n_ncdm] = p_ncdm;
      p_tot += p_ncdm;
      pvecback[index_bg_pseudo_p_ncdm1_ + n_ncdm] = pseudo_p_ncdm;
      /** See e.g. Eq. A6 in 1811.00904. */
      dp_dloga += (pseudo_p_ncdm - 5*p_ncdm);

      /* (3 p_ncdm1) is the "relativistic" contribution to rho_ncdm1 */
      rho_r += 3.* p_ncdm;

      /* (rho_ncdm1 - 3 p_ncdm1) is the "non-relativistic" contribution
         to rho_ncdm1 */
      rho_m += rho_ncdm - 3.* p_ncdm;
    }
  }

  /* Lambda */
  if (pba->has_lambda == _TRUE_) {
    pvecback[index_bg_rho_lambda_] = pba->Omega0_lambda*pow(pba->H0, 2);
    rho_tot += pvecback[index_bg_rho_lambda_];
    p_tot -= pvecback[index_bg_rho_lambda_];
  }

  /* fluid with w(a) and constant cs2 */
  if (pba->has_fld == _TRUE_) {

    /* get rho_fld from vector of integrated variables */
    pvecback[index_bg_rho_fld_] = pvecback_B[index_bi_rho_fld_];

    /* get w_fld from dedicated function */
    class_call(background_w_fld(a, &w_fld, &dw_over_da, &integral_fld), error_message_, error_message_);
    pvecback[index_bg_w_fld_] = w_fld;

    // Obsolete: at the beginning, we had here the analytic integral solution corresponding to the case w=w0+w1(1-a/a0):
    // pvecback[index_bg_rho_fld_] = pba->Omega0_fld * pow(pba->H0,2) / pow(a_rel,3.*(1.+pba->w0_fld+pba->wa_fld)) * exp(3.*pba->wa_fld*(a_rel-1.));
    // But now everthing is integrated numerically for a given w_fld(a) defined in the function background_w_fld.

    rho_tot += pvecback[index_bg_rho_fld_];
    p_tot += w_fld*pvecback[index_bg_rho_fld_];
    dp_dloga += (a*dw_over_da - 3*(1 + w_fld)*w_fld)*pvecback[index_bg_rho_fld_];
  }

  /* relativistic neutrinos (and all relativistic relics) */
  if (pba->has_ur == _TRUE_) {
    pvecback[index_bg_rho_ur_] = pba->Omega0_ur*pow(pba->H0, 2)/pow(a_rel, 4);
    rho_tot += pvecback[index_bg_rho_ur_];
    p_tot += 1./3.*pvecback[index_bg_rho_ur_];
    dp_dloga += -4./3.*pvecback[index_bg_rho_ur_];
    rho_r += pvecback[index_bg_rho_ur_];
  }

  /* interacting dark matter */
  if (pba->has_idm_dr == _TRUE_) {
    pvecback[index_bg_rho_idm_dr_] = pba->Omega0_idm_dr*pow(pba->H0, 2)/pow(a_rel, 3);
    rho_tot += pvecback[index_bg_rho_idm_dr_];
    p_tot += 0.;
    rho_m += pvecback[index_bg_rho_idm_dr_];
  }

  /* interacting dark radiation */
  if (pba->has_idr == _TRUE_) {
    pvecback[index_bg_rho_idr_] = pba->Omega0_idr*pow(pba->H0, 2)/pow(a_rel, 4);
    rho_tot += pvecback[index_bg_rho_idr_];
    p_tot += 1./3.*pvecback[index_bg_rho_idr_];
    rho_r += pvecback[index_bg_rho_idr_];
  }

  /** - compute expansion rate H from Friedmann equation: this is the
      only place where the Friedmann equation is assumed. Remember
      that densities are all expressed in units of \f$ [3c^2/8\pi G] \f$, ie
      \f$ \rho_{class} = [8 \pi G \rho_{physical} / 3 c^2]\f$ */
  pvecback[index_bg_H_] = sqrt(rho_tot-pba->K/a/a);

  /** - compute derivative of H with respect to conformal time */
  pvecback[index_bg_H_prime_] = -3./2.*(rho_tot + p_tot)*a + pba->K/a;

  /* Total energy density*/
  pvecback[index_bg_rho_tot_] = rho_tot;

  /* Total pressure */
  pvecback[index_bg_p_tot_] = p_tot;

  /* Derivative of total pressure w.r.t. conformal time */
  pvecback[index_bg_p_tot_prime_] = a*pvecback[index_bg_H_]*dp_dloga;
  if (pba->has_scf == _TRUE_){
    /** The contribution of scf was not added to dp_dloga, add p_scf_prime here: */
    pvecback[index_bg_p_prime_scf_] = pvecback[index_bg_phi_prime_scf_]*
      (-pvecback[index_bg_phi_prime_scf_]*pvecback[index_bg_H_]/a - 2./3.*pvecback[index_bg_dV_scf_]);
    pvecback[index_bg_p_tot_prime_] += pvecback[index_bg_p_prime_scf_];
  }

  /** - compute critical density */
  rho_crit = rho_tot-pba->K/a/a;
  class_test(rho_crit <= 0.,
             error_message_,
             "rho_crit = %e instead of strictly positive",rho_crit);

  /** - compute relativistic density to total density ratio */
  pvecback[index_bg_Omega_r_] = rho_r/rho_crit;

  /** - compute other quantities in the exhaustive, redundant format */
  if (return_format == pba->long_info) {

    /** - store critical density */
    pvecback[index_bg_rho_crit_] = rho_crit;

    /** - compute Omega_m */
    pvecback[index_bg_Omega_m_] = rho_m/rho_crit;

    /* one can put other variables here */
    /*  */
    /*  */

  }

  return _SUCCESS_;

}

/**
 * Single place where the fluid equation of state is
 * defined. Parameters of the function are passed through the
 * background structure. Generalisation to arbitrary functions should
 * be simple.
 *
 * @param pba            Input: pointer to background structure
 * @param a              Input: current value of scale factor
 * @param w_fld          Output: equation of state parameter w_fld(a)
 * @param dw_over_da_fld Output: function dw_fld/da
 * @param integral_fld   Output: function \f$ \int_{a}^{a_0} da 3(1+w_{fld})/a \f$
 * @return the error status
 */

int BackgroundModule::background_w_fld(double a, double* w_fld, double* dw_over_da_fld, double* integral_fld) const {

  double Omega_ede = 0.;
  double dOmega_ede_over_da = 0.;
  double d2Omega_ede_over_da2 = 0.;
  double a_eq, Omega_r, Omega_m;

  /** - first, define the function w(a) */
  switch (pba->fluid_equation_of_state) {
  case CLP:
    *w_fld = pba->w0_fld + pba->wa_fld * (1. - a / pba->a_today);
    break;
  case EDE:
    // Omega_ede(a) taken from eq. (10) in 1706.00730
    Omega_ede = (pba->Omega0_fld - pba->Omega_EDE*(1.-pow(a,-3.*pba->w0_fld)))
      /(pba->Omega0_fld+(1.-pba->Omega0_fld)*pow(a,3.*pba->w0_fld))
      + pba->Omega_EDE*(1.-pow(a,-3.*pba->w0_fld));

    // d Omega_ede / d a taken analytically from the above
    dOmega_ede_over_da = - pba->Omega_EDE* 3.*pba->w0_fld*pow(a,-3.*pba->w0_fld-1.)/(pba->Omega0_fld+(1.-pba->Omega0_fld)*pow(a,3.*pba->w0_fld))
      - (pba->Omega0_fld - pba->Omega_EDE*(1.-pow(a,-3.*pba->w0_fld)))*(1.-pba->Omega0_fld)*3.*pba->w0_fld*pow(a,3.*pba->w0_fld-1.)/pow(pba->Omega0_fld+(1.-pba->Omega0_fld)*pow(a,3.*pba->w0_fld),2)
      + pba->Omega_EDE*3.*pba->w0_fld*pow(a,-3.*pba->w0_fld-1.);

    // find a_equality (needed because EDE tracks first radiation, then matter)
    Omega_r = pba->Omega0_g * (1. + 3.046 * 7./8.*pow(4./11.,4./3.)); // assumes LambdaCDM + eventually massive neutrinos so light that they are relativistic at equality; needs to be generalised later on.
    Omega_m = pba->Omega0_b;
    if (pba->has_cdm == _TRUE_) Omega_m += pba->Omega0_cdm;
    if (pba->has_idm_dr == _TRUE_) Omega_m += pba->Omega0_idm_dr;
    if (pba->has_dcdm == _TRUE_)
      class_stop(error_message_, "Early Dark Energy not compatible with decaying Dark Matter because we omitted to code the calculation of a_eq in that case, but it would not be difficult to add it if necessary, should be a matter of 5 minutes");
    a_eq = Omega_r/Omega_m; // assumes a flat universe with a=1 today

    // w_ede(a) taken from eq. (11) in 1706.00730
    *w_fld = - dOmega_ede_over_da*a/Omega_ede/3./(1.-Omega_ede)+a_eq/3./(a+a_eq);
    break;
  }


  /** - then, give the corresponding analytic derivative dw/da (used
      by perturbation equations; we could compute it numerically,
      but with a loss of precision; as long as there is a simple
      analytic expression of the derivative of the previous
      function, let's use it! */
  switch (pba->fluid_equation_of_state) {
  case CLP:
    *dw_over_da_fld = - pba->wa_fld / pba->a_today;
    break;
  case EDE:
    d2Omega_ede_over_da2 = 0.;
    *dw_over_da_fld = - d2Omega_ede_over_da2*a/3./(1.-Omega_ede)/Omega_ede
      - dOmega_ede_over_da/3./(1.-Omega_ede)/Omega_ede
      + dOmega_ede_over_da*dOmega_ede_over_da*a/3./(1.-Omega_ede)/(1.-Omega_ede)/Omega_ede
      + a_eq/3./(a+a_eq)/(a+a_eq);
    break;
  }

  /** - finally, give the analytic solution of the following integral:
        \f$ \int_{a}^{a0} da 3(1+w_{fld})/a \f$. This is used in only
        one place, in the initial conditions for the background, and
        with a=a_ini. If your w(a) does not lead to a simple analytic
        solution of this integral, no worry: instead of writing
        something here, the best would then be to leave it equal to
        zero, and then in background_initial_conditions() you should
        implement a numerical calculation of this integral only for
        a=a_ini, using for instance Romberg integration. It should be
        fast, simple, and accurate enough. */
  switch (pba->fluid_equation_of_state) {
  case CLP:
    *integral_fld = 3.*((1.+pba->w0_fld+pba->wa_fld)*log(pba->a_today/a) + pba->wa_fld*(a/pba->a_today-1.));
    break;
  case EDE:
    class_stop(error_message_, "EDE implementation not finished: to finish it, read the comments in background.c just before this line\n");
    break;
  }

  /** note: of course you can generalise these formulas to anything,
      defining new parameters pba->w..._fld. Just remember that so
      far, HyRec explicitely assumes that w(a)= w0 + wa (1-a/a0); but
      Recfast does not assume anything */

  return _SUCCESS_;
}

/**
 * Initialize the background structure, and in particular the
 * background interpolation table.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input/Output: pointer to initialized background structure
 * @return the error status
 */

int BackgroundModule::background_init() {

  /** Summary: */

  /** - define local variables */
  double rho_nu_rel;
  double Neff, N_dark;
  double w_fld, dw_over_da, integral_fld;

  /** - in verbose mode, provide some information */
  if (pba->background_verbose > 0) {
    printf("Running CLASS version %s\n",_VERSION_);
    printf("Computing background\n");

    /* below we want to inform the user about ncdm species and/or the total N_eff */
    if ((pba->N_ncdm > 0) || (pba->Omega0_idr != 0.))  {

      /* contribution of ultra-relativistic species _ur to N_eff */
      Neff = pba->Omega0_ur/7.*8./pow(4./11.,4./3.)/pba->Omega0_g;

      /* contribution of ncdm species to N_eff*/
      if (pba->N_ncdm > 0){
        Neff += pba->ncdm->GetNeff(0.);
        pba->ncdm->PrintNeffInfo();
      }

      /* contribution of interacting dark radiation _idr to N_eff */
      if (pba->Omega0_idr != 0.) {
        N_dark = pba->Omega0_idr/7.*8./pow(4./11.,4./3.)/pba->Omega0_g;
        Neff += N_dark;
        printf(" -> dark radiation Delta Neff %e\n",N_dark);
      }

      printf(" -> total N_eff = %g (sumed over ultra-relativistic species, ncdm and dark radiation)\n",Neff);

    }
  }

  /** - assign values to all indices in vectors of background quantities with background_indices()*/
  class_call(background_indices(),
             error_message_,
             error_message_);

  /* fluid equation of state */
  if (pba->has_fld == _TRUE_) {

    class_call(background_w_fld(0., &w_fld, &dw_over_da, &integral_fld), error_message_, error_message_);

    class_test(w_fld >= 1./3.,
               error_message_,
               "Your choice for w(a--->0)=%g is suspicious, since it is bigger than -1/3 there cannot be radiation domination at early times\n",
               w_fld);
  }

  /* in verbose mode, inform the user about the value of the ncdm
     masses in eV and about the ratio [m/omega_ncdm] in eV (the usual
     93 point something)*/
  if ((pba->background_verbose > 0) && (pba->has_ncdm == _TRUE_)) {
    pba->ncdm->PrintMassInfo();
  }

  /* check other quantities which would lead to segmentation fault if zero */
  class_test(pba->a_today <= 0,
             error_message_,
             "input a_today = %e instead of strictly positive",pba->a_today);

  class_test(_Gyr_over_Mpc_ <= 0,
             error_message_,
             "_Gyr_over_Mpc = %e instead of strictly positive",_Gyr_over_Mpc_);

  /** - this function integrates the background over time, allocates
      and fills the background table */
  switch (pba->background_method) {
  case (bgevo_rk):
      class_call(background_solve(), error_message_, error_message_);
      break;
  case (bgevo_evolver):
      class_call(background_solve_evolver(), error_message_, error_message_);
      break;
  default:
    printf("No background evolution method selected!\n");
  }

  /** - this function finds and stores a few derived parameters at radiation-matter equality */
  class_call(background_find_equality(), error_message_, error_message_);

  class_call(background_output_budget(), error_message_, error_message_);

  return _SUCCESS_;

}

/**
 * Free all memory space allocated by background_init().
 *
 *
 * @param pba Input: pointer to background structure (to be freed)
 * @return the error status
 */

int BackgroundModule::background_free() {

  class_call(background_free_noinput(),
              error_message_,
              error_message_);

  return _SUCCESS_;
}

/**
 * Free only the memory space NOT allocated through input_read_parameters()
 *
 * @param pba Input: pointer to background structure (to be freed)
 * @return the error status
 */

int BackgroundModule::background_free_noinput() const {

  free(tau_table_);
  free(z_table_);
  free(d2tau_dz2_table_);
  free(background_table_);
  free(d2background_dtau2_table_);

  return _SUCCESS_;
}

/**
 * Assign value to each relevant index in vectors of background quantities.
 *
 * @param pba Input: pointer to background structure
 * @return the error status
 */

int BackgroundModule::background_indices() {

  /** Summary: */

  /** - define local variables */

  /* a running index for the vector of background quantities */
  int index_bg;
  /* a running index for the vector of background quantities to be integrated */
  int index_bi;

  /** - initialize all indices */

  index_bg=0;

  /* index for scale factor */
  class_define_index(index_bg_a_, _TRUE_, index_bg, 1);

  /* - indices for H and its conformal-time-derivative */
  class_define_index(index_bg_H_, _TRUE_, index_bg, 1);
  class_define_index(index_bg_H_prime_, _TRUE_, index_bg, 1);

  /* - end of indices in the short vector of background values */
  bg_size_short_ = index_bg;

  /* - index for rho_g (photon density) */
  class_define_index(index_bg_rho_g_, _TRUE_, index_bg, 1);

  /* - index for rho_b (baryon density) */
  class_define_index(index_bg_rho_b_, _TRUE_, index_bg, 1);

  /* - index for rho_cdm */
  class_define_index(index_bg_rho_cdm_, pba->has_cdm, index_bg, 1);

  /* - indices for ncdm. We only define the indices for ncdm1
     (density, pressure, pseudo-pressure), the other ncdm indices
     are contiguous */
  class_define_index(index_bg_rho_ncdm1_ ,pba->has_ncdm, index_bg, pba->N_ncdm);
  class_define_index(index_bg_p_ncdm1_, pba->has_ncdm, index_bg, pba->N_ncdm);
  class_define_index(index_bg_pseudo_p_ncdm1_, pba->has_ncdm, index_bg, pba->N_ncdm);

  /* - index for dcdm */
  class_define_index(index_bg_rho_dcdm_, pba->has_dcdm, index_bg, 1);

  /* - index for dr */
  class_define_index(index_bg_rho_dr_, pba->has_dr, index_bg, 1);

  /* - indices for scalar field */
  class_define_index(index_bg_phi_scf_, pba->has_scf, index_bg, 1);
  class_define_index(index_bg_phi_prime_scf_, pba->has_scf, index_bg, 1);
  class_define_index(index_bg_V_scf_, pba->has_scf, index_bg, 1);
  class_define_index(index_bg_dV_scf_, pba->has_scf,index_bg, 1);
  class_define_index(index_bg_ddV_scf_, pba->has_scf,index_bg, 1);
  class_define_index(index_bg_rho_scf_, pba->has_scf,index_bg, 1);
  class_define_index(index_bg_p_scf_, pba->has_scf, index_bg, 1);
  class_define_index(index_bg_p_prime_scf_, pba->has_scf, index_bg, 1);

  /* - index for Lambda */
  class_define_index(index_bg_rho_lambda_, pba->has_lambda,index_bg, 1);

  /* - index for fluid */
  class_define_index(index_bg_rho_fld_, pba->has_fld, index_bg, 1);
  class_define_index(index_bg_w_fld_, pba->has_fld, index_bg, 1);

  /* - index for ultra-relativistic neutrinos/species */
  class_define_index(index_bg_rho_ur_, pba->has_ur, index_bg, 1);

  /* - index for total density */
  class_define_index(index_bg_rho_tot_, _TRUE_, index_bg, 1);

  /* - index for total pressure */
  class_define_index(index_bg_p_tot_, _TRUE_, index_bg, 1);

  /* - index for derivative of total pressure */
  class_define_index(index_bg_p_tot_prime_, _TRUE_, index_bg, 1);

  /* - index for Omega_r (relativistic density fraction) */
  class_define_index(index_bg_Omega_r_, _TRUE_, index_bg, 1);

  /* - index interacting for dark radiation */
  class_define_index(index_bg_rho_idr_, pba->has_idr, index_bg, 1);

  /* - index for interacting dark matter */
  class_define_index(index_bg_rho_idm_dr_, pba->has_idm_dr, index_bg, 1);

  /* - put here additional ingredients that you want to appear in the
     normal vector */
  /*    */
  /*    */

  /* - end of indices in the normal vector of background values */
  bg_size_normal_ = index_bg;

  /* - indices in the long version : */

  /* -> critical density */
  class_define_index(index_bg_rho_crit_, _TRUE_, index_bg, 1);

  /* - index for Omega_m (non-relativistic density fraction) */
  class_define_index(index_bg_Omega_m_, _TRUE_, index_bg, 1);

  /* -> conformal distance */
  class_define_index(index_bg_conf_distance_, _TRUE_, index_bg, 1);

  /* -> angular diameter distance */
  class_define_index(index_bg_ang_distance_, _TRUE_, index_bg, 1);

  /* -> luminosity distance */
  class_define_index(index_bg_lum_distance_, _TRUE_, index_bg, 1);

  /* -> proper time (for age of the Universe) */
  class_define_index(index_bg_time_, _TRUE_, index_bg, 1);

  /* -> conformal sound horizon */
  class_define_index(index_bg_rs_, _TRUE_, index_bg, 1);

  /* -> density growth factor in dust universe */
  class_define_index(index_bg_D_, _TRUE_, index_bg, 1);

  /* -> velocity growth factor in dust universe */
  class_define_index(index_bg_f_, _TRUE_, index_bg, 1);

  /* -> put here additional quantities describing background */
  /*    */
  /*    */

  /* -> end of indices in the long vector of background values */
  bg_size_ = index_bg;

  /* - now, indices in vector of variables to integrate.
     First {B} variables, then {C} variables. */

  index_bi=0;

  /* -> scale factor */
  class_define_index(index_bi_a_, _TRUE_, index_bi, 1);

  /* -> energy density in DCDM */
  class_define_index(index_bi_rho_dcdm_, pba->has_dcdm, index_bi, 1);

  /* -> energy density in DR */
  class_define_index(index_bi_rho_dr_, pba->has_dr, index_bi, 1);

  /* -> energy density in fluid */
  class_define_index(index_bi_rho_fld_, pba->has_fld, index_bi, 1);

  /* -> scalar field and its derivative wrt conformal time (Zuma) */
  class_define_index(index_bi_phi_scf_, pba->has_scf, index_bi, 1);
  class_define_index(index_bi_phi_prime_scf_, pba->has_scf, index_bi, 1);

  /* End of {B} variables, now continue with {C} variables */
  bi_B_size_ = index_bi;

  /* -> proper time (for age of the Universe) */
  class_define_index(index_bi_time_, _TRUE_, index_bi, 1);

  /* -> sound horizon */
  class_define_index(index_bi_rs_, _TRUE_, index_bi, 1);

  /* -> Second order equation for growth factor */
  class_define_index(index_bi_D_, _TRUE_, index_bi, 1);
  class_define_index(index_bi_D_prime_, _TRUE_, index_bi, 1);

  /* -> index for conformal time in vector of variables to integrate */
  class_define_index(index_bi_tau_, _TRUE_, index_bi, 1);

  /* -> end of indices in the vector of variables to integrate */
  bi_size_ = index_bi;

  /* index_bi_tau must be the last index, because tau is part of this vector for the purpose of being stored, */
  /* but it is not a quantity to be integrated (since integration is over tau itself) */
  class_test(index_bi_tau_ != index_bi - 1,
             error_message_,
             "background integration requires index_bi_tau to be the last of all index_bi's");

  return _SUCCESS_;

}

/**
 *  This function integrates the background over time, allocates and
 *  fills the background table
 *
 * @param ppr Input: precision structure
 * @param pba Input/Output: background structure
 */

int BackgroundModule::background_solve() {

  /** Summary: */

  /** - define local variables */

  /* contains all quantities relevant for the integration algorithm */
  struct generic_integrator_workspace gi;
  /* parameters and workspace for the background_derivs function */
  background_parameters_and_workspace bpaw{this};
  /* a growing table (since the number of time steps is not known a priori) */
  growTable gTable;
  /* needed for growing table */
  double * pData;
  /* needed for growing table */
  void * memcopy_result;
  /* initial conformal time */
  double tau_start;
  /* final conformal time */
  double tau_end;
  /* an index running over bi indices */
  int i;
  /* vector of quantities to be integrated */
  double * pvecback_integration;
  /* vector of all background quantities */
  double * pvecback;
  /* necessary for calling array_interpolate(), but never used */
  int last_index=0;
  /* comoving radius coordinate in Mpc (equal to conformal distance in flat case) */
  double comoving_radius=0.;

  class_alloc(pvecback,bg_size_*sizeof(double), error_message_);
  bpaw.pvecback = pvecback;

  /** - allocate vector of quantities to be integrated */
  class_alloc(pvecback_integration, bi_size_*sizeof(double), error_message_);

  /** - initialize generic integrator with initialize_generic_integrator() */

  /* Size of vector to integrate is (bi_size_ - 1) rather than
   * (bi_size_), since tau is not integrated.
   */
  class_call(initialize_generic_integrator(bi_size_ - 1, &gi),
             gi.error_message,
             error_message_);

  /** - impose initial conditions with background_initial_conditions() */
  class_call(background_initial_conditions(pvecback, pvecback_integration),
             error_message_,
             error_message_);

  /* here tau_end is in fact the initial time (in the next loop
     tau_start = tau_end) */
  tau_end=pvecback_integration[index_bi_tau_];

  /** - create a growTable with gt_init() */
  class_call(gt_init(&gTable),
             gTable.error_message,
             error_message_);

  /* initialize the counter for the number of steps */
  bt_size_=0;

  /** - loop over integration steps: call background_functions(), find step size, save data in growTable with gt_add(), perform one step with generic_integrator(), store new value of tau */

  while (pvecback_integration[index_bi_a_] < pba->a_today) {

    tau_start = tau_end;

    /* -> find step size (trying to adjust the last step as close as possible to the one needed to reach a=a_today; need not be exact, difference corrected later) */
    class_call(background_functions(pvecback_integration, pba->short_info, pvecback),
               error_message_,
               error_message_);

    if ((pvecback_integration[index_bi_a_]*(1. + ppr->back_integration_stepsize)) < pba->a_today) {
      tau_end = tau_start + ppr->back_integration_stepsize/(pvecback_integration[index_bi_a_]*pvecback[index_bg_H_]);
      /* no possible segmentation fault here: non-zeroness of "a" has been checked in background_functions() */
    }
    else {
      tau_end = tau_start + (pba->a_today/pvecback_integration[index_bi_a_] - 1.)/(pvecback_integration[index_bi_a_]*pvecback[index_bg_H_]);
      /* no possible segmentation fault here: non-zeroness of "a" has been checked in background_functions() */
    }

    class_test((tau_end-tau_start)/tau_start < ppr->smallest_allowed_variation,
               error_message_,
               "integration step: relative change in time =%e < machine precision : leads either to numerical error or infinite loop",(tau_end-tau_start)/tau_start);

    /* -> save data in growTable */
    class_call(gt_add(&gTable,_GT_END_, (void*) pvecback_integration, sizeof(double)*bi_size_),
               gTable.error_message,
               error_message_);
    bt_size_++;

    /* -> perform one step */
    class_call(generic_integrator(background_derivs,
                                  tau_start,
                                  tau_end,
                                  pvecback_integration,
                                  &bpaw,
                                  ppr->tol_background_integration,
                                  ppr->smallest_allowed_variation,
                                  &gi),
               gi.error_message,
               error_message_);

    /* -> store value of tau */
    pvecback_integration[index_bi_tau_] = tau_end;

  }

  /** - save last data in growTable with gt_add() */
  class_call(gt_add(&gTable, _GT_END_, (void*) pvecback_integration, sizeof(double)*bi_size_),
             gTable.error_message,
             error_message_);
  bt_size_++;


  /* integration finished */

  /** - clean up generic integrator with cleanup_generic_integrator() */
  class_call(cleanup_generic_integrator(&gi),
             gi.error_message,
             error_message_);

  /** - retrieve data stored in the growTable with gt_getPtr() */
  class_call(gt_getPtr(&gTable,(void**)&pData),
             gTable.error_message,
             error_message_);

  /** - interpolate to get quantities precisely today with array_interpolate() */
  class_call(array_interpolate(
                               pData,
                               bi_size_,
                               bt_size_,
                               index_bi_a_,
                               pba->a_today,
                               &last_index,
                               pvecback_integration,
                               bi_size_,
                               error_message_),
             error_message_,
             error_message_);

  /* substitute last line with quantities today */
  for (i = 0; i < bi_size_; i++)
    pData[(bt_size_ - 1)*bi_size_ + i] = pvecback_integration[i];

  /** - deduce age of the Universe */
  /* -> age in Gyears */
  age_ = pvecback_integration[index_bi_time_]/_Gyr_over_Mpc_;
  /* -> conformal age in Mpc */
  conformal_age_ = pvecback_integration[index_bi_tau_];
  /* -> contribution of decaying dark matter and dark radiation to the critical density today: */
  if (pba->has_dcdm == _TRUE_){
    Omega0_dcdm_ = pvecback_integration[index_bi_rho_dcdm_]/pba->H0/pba->H0;
  }
  if (pba->has_dr == _TRUE_){
    Omega0_dr_ = pvecback_integration[index_bi_rho_dr_]/pba->H0/pba->H0;
  }

  /** - allocate background tables */
  class_alloc(tau_table_, bt_size_*sizeof(double), error_message_);

  class_alloc(z_table_, bt_size_*sizeof(double), error_message_);

  class_alloc(d2tau_dz2_table_, bt_size_*sizeof(double), error_message_);

  class_alloc(background_table_, bt_size_*bg_size_*sizeof(double), error_message_);

  class_alloc(d2background_dtau2_table_, bt_size_*bg_size_*sizeof(double), error_message_);

  /** - In a loop over lines, fill background table using the result of the integration plus background_functions() */
  for (i = 0; i < bt_size_; i++) {

    /* -> establish correspondence between the integrated variable and the bg variables */

    tau_table_[i] = pData[i*bi_size_ + index_bi_tau_];

    class_test(pData[i*bi_size_ + index_bi_a_] <= 0.,
               error_message_,
               "a = %e instead of strictly positiv",
               pData[i*bi_size_ + index_bi_a_]);

    z_table_[i] = pba->a_today/pData[i*bi_size_ + index_bi_a_] - 1.;

    pvecback[index_bg_time_] = pData[i*bi_size_ + index_bi_time_];
    pvecback[index_bg_conf_distance_] = conformal_age_ - pData[i*bi_size_ + index_bi_tau_];

    if (pba->sgnK == 0) comoving_radius = pvecback[index_bg_conf_distance_];
    else if (pba->sgnK == 1) comoving_radius = sin(sqrt(pba->K)*pvecback[index_bg_conf_distance_])/sqrt(pba->K);
    else if (pba->sgnK == -1) comoving_radius = sinh(sqrt(-pba->K)*pvecback[index_bg_conf_distance_])/sqrt(-pba->K);

    pvecback[index_bg_ang_distance_] = pba->a_today*comoving_radius/(1. + z_table_[i]);
    pvecback[index_bg_lum_distance_] = pba->a_today*comoving_radius*(1. + z_table_[i]);
    pvecback[index_bg_rs_] = pData[i*bi_size_ + index_bi_rs_];

    /* -> compute all other quantities depending only on {B} variables.
       The value of {B} variables in pData are also copied to pvecback.*/
    class_call(background_functions(pData + i*bi_size_, pba->long_info, pvecback),
               error_message_,
               error_message_);

    /* -> compute growth functions (valid in dust universe) */

    /* Normalise D(z=0)=1 and construct f = D_prime/(aHD) */
    pvecback[index_bg_D_] = pData[i*bi_size_ + index_bi_D_]/pData[(bt_size_ - 1)*bi_size_ + index_bi_D_];
    pvecback[index_bg_f_] = pData[i*bi_size_ + index_bi_D_prime_]/
      (pData[i*bi_size_ + index_bi_D_]*pvecback[index_bg_a_]*pvecback[index_bg_H_]);

    /* -> write in the table */
    memcopy_result = memcpy(background_table_ + i*bg_size_, pvecback, bg_size_*sizeof(double));

    class_test(memcopy_result != background_table_ + i*bg_size_,
               error_message_,
               "cannot copy data back to background_table_");
  }

  /** - free the growTable with gt_free() */

  class_call(gt_free(&gTable),
             gTable.error_message,
             error_message_);

  /** - fill tables of second derivatives (in view of spline interpolation) */
  class_call(array_spline_table_lines(z_table_,
                                      bt_size_,
                                      tau_table_,
                                      1,
                                      d2tau_dz2_table_,
                                      _SPLINE_EST_DERIV_,
                                      error_message_),
             error_message_,
             error_message_);

  class_call(array_spline_table_lines(tau_table_,
                                      bt_size_,
                                      background_table_,
                                      bg_size_,
                                      d2background_dtau2_table_,
                                      _SPLINE_EST_DERIV_,
                                      error_message_),
             error_message_,
             error_message_);

  /** - compute remaining "related parameters" */

  /**  - so-called "effective neutrino number", computed at earliest
      time in interpolation table. This should be seen as a
      definition: Neff is the equivalent number of
      instantaneously-decoupled neutrinos accounting for the
      radiation density, beyond photons */

  Neff_ = (background_table_[index_bg_Omega_r_]*background_table_[index_bg_rho_crit_] - background_table_[index_bg_rho_g_])
    /(7./8.*pow(4./11.,4./3.)*background_table_[index_bg_rho_g_]);

  /** - done */
  if (pba->background_verbose > 0) {
    printf(" -> age = %f Gyr\n", age_);
    printf(" -> conformal age = %f Mpc\n", conformal_age_);
  }

  if (pba->background_verbose > 2) {
    printf(" -> Neff_ = %f\n", Neff_);
    if ((pba->has_dcdm == _TRUE_)&&(pba->has_dr == _TRUE_)){
      printf("    Decaying Cold Dark Matter details: (DCDM --> DR)\n");
      printf("     -> Omega0_dcdm = %f\n", Omega0_dcdm_);
      printf("     -> Omega0_dr = %f\n", Omega0_dr_);
      printf("     -> Omega0_dr+Omega0_dcdm = %f, input value = %f\n", Omega0_dr_ + Omega0_dcdm_, pba->Omega0_dcdmdr);
      printf("     -> Omega_ini_dcdm/Omega_b = %f\n",pba->Omega_ini_dcdm/pba->Omega0_b);
    }
    if (pba->has_scf == _TRUE_){
      printf("    Scalar field details:\n");
      printf("     -> Omega_scf = %g, wished %g\n", pvecback[index_bg_rho_scf_]/pvecback[index_bg_rho_crit_], pba->Omega0_scf);
      if(pba->has_lambda == _TRUE_)
        printf("     -> Omega_Lambda = %g, wished %g\n", pvecback[index_bg_rho_lambda_]/pvecback[index_bg_rho_crit_], pba->Omega0_lambda);
      printf("     -> parameters: [lambda, alpha, A, B] = \n");
      printf("                    [");
      for (i = 0; i < pba->scf_parameters.size() - 1; i++){
        printf("%.3f, ",pba->scf_parameters[i]);
      }
      printf("%.3f]\n", pba->scf_parameters[pba->scf_parameters.size() - 1]);
    }
  }

  /**  - total matter, radiation, dark energy today */
  Omega0_m_ = background_table_[(bt_size_ - 1)*bg_size_ + index_bg_Omega_m_];
  Omega0_r_ = background_table_[(bt_size_ - 1)*bg_size_ + index_bg_Omega_r_];
  Omega0_de_ = 1. - (Omega0_m_ + Omega0_r_ + pba->Omega0_k);

  free(pvecback);
  free(pvecback_integration);

  return _SUCCESS_;

}

int BackgroundModule::background_solve_evolver() {

  /** Summary: */

  /** - define local variables */

  /* contains all quantities relevant for the integration algorithm */
  struct generic_integrator_workspace gi;
  /* parameters and workspace for the background_derivs function */
  struct background_parameters_and_workspace bpaw{this};
  /* vector of quantities to be integrated */
  double* pvecback_integration;
  /* vector of all background quantities */
  double* pvecback;
  /* necessary for calling array_interpolate(), but never used */
  int last_index = 0;
  /* comoving radius coordinate in Mpc (equal to conformal distance in flat case) */
  double comoving_radius = 0.;

  class_alloc(pvecback, bg_size_*sizeof(double), error_message_);
  bpaw.pvecback = pvecback;
  /** - allocate vector of quantities to be integrated */
  class_alloc(pvecback_integration, bi_size_*sizeof(double), error_message_);

  /* Size of vector to integrate is (bi_size_-1) rather than
   * (bi_size_), since a is not integrated.
   */

  /** - impose initial conditions with background_initial_conditions() */
  class_call(background_initial_conditions(pvecback, pvecback_integration),
             error_message_,
             error_message_);

  /** - Determine output vector */
  double loga_ini = log(pvecback_integration[index_bi_a_]);
  double loga_final = log(pba->a_today);
  bt_size_ = (loga_final - loga_ini)/ppr->back_integration_stepsize;
  double* loga;
  int* used_in_output;
  class_alloc(loga, bt_size_*sizeof(double), error_message_);
  class_alloc(used_in_output, bt_size_*sizeof(int), error_message_);
  for (int index_loga = 0; index_loga < bt_size_; index_loga++){
    loga[index_loga] = loga_ini + index_loga*(loga_final - loga_ini)/(bt_size_ - 1);
    used_in_output[index_loga] = 1;
  }

  /** - Remember that we evolve tau at index_bi_a: */
  pvecback_integration[index_bi_a_] =  pvecback_integration[index_bi_tau_];

  /** - allocate background tables */
  class_alloc(tau_table_, bt_size_*sizeof(double), error_message_);
  class_alloc(z_table_, bt_size_*sizeof(double), error_message_);
  class_alloc(d2tau_dz2_table_, bt_size_*sizeof(double), error_message_);
  class_alloc(background_table_, bt_size_*bg_size_*sizeof(double), error_message_);
  class_alloc(d2background_dtau2_table_, bt_size_*bg_size_*sizeof(double), error_message_);

  auto generic_evolver = &evolver_ndf15;
  if (ppr->evolver == rk){
      generic_evolver = &evolver_rk;
  }

  class_call(generic_evolver(background_derivs_loga,
           loga_ini,
           loga_final,
           pvecback_integration,
           used_in_output,
           bi_size_ - 1,
           &bpaw,
           1e-6,
           ppr->smallest_allowed_variation,
           NULL,
           ppr->perturb_integration_stepsize,
           loga,
           bt_size_,
           background_add_line_to_bg_table,
           NULL,
           error_message_),
       error_message_,
       error_message_);


  /** - deduce age of the Universe */
  /* -> age in Gyears */
  age_ = pvecback_integration[index_bi_time_]/_Gyr_over_Mpc_;
  /* -> conformal age in Mpc. Remember that tau is stored at index_bi_a now */
  conformal_age_ = pvecback_integration[index_bi_a_];
  /* -> contribution of decaying dark matter and dark radiation to the critical density today: */
  if (pba->has_dcdm == _TRUE_){
    Omega0_dcdm_ = pvecback_integration[index_bi_rho_dcdm_]/pba->H0/pba->H0;
  }
  if (pba->has_dr == _TRUE_){
    Omega0_dr_ = pvecback_integration[index_bi_rho_dr_]/pba->H0/pba->H0;
  }

  /** - In a loop over lines, fill rest of background table using the result of the integration plus background_functions() */
  /** Recover some quantities today */
  double D_today = pvecback_integration[index_bi_D_];
  for (int i = 0; i < bt_size_; i++) {
    double* bg_table_row = background_table_ + i*bg_size_;
    /** Recover integration vector at this time */
    for (int index_bi = 0; index_bi < bi_size_; index_bi++) {
      pvecback_integration[index_bi] = bg_table_row[index_bi];
    }
    /** bg_table_row is the row in the background table, pvecback_integration holds the contents of the integrated vector
  at the current time */

    /* -> establish correspondence between the integrated variable and the bg variables */
    bg_table_row[index_bg_time_] = pvecback_integration[index_bi_time_];
    bg_table_row[index_bg_conf_distance_] = conformal_age_ - pvecback_integration[index_bi_tau_];

    if (pba->sgnK == 0) comoving_radius = bg_table_row[index_bg_conf_distance_];
    else if (pba->sgnK == 1) comoving_radius = sin(sqrt(pba->K)*bg_table_row[index_bg_conf_distance_])/sqrt(pba->K);
    else if (pba->sgnK == -1) comoving_radius = sinh(sqrt(-pba->K)*bg_table_row[index_bg_conf_distance_])/sqrt(-pba->K);

    bg_table_row[index_bg_ang_distance_] = pba->a_today*comoving_radius/(1. + z_table_[i]);
    bg_table_row[index_bg_lum_distance_] = pba->a_today*comoving_radius*(1. + z_table_[i]);
    bg_table_row[index_bg_rs_] = pvecback_integration[index_bi_rs_];

    /* -> compute all other quantities depending only on {B} variables.
       The value of {B} variables in pData are also copied to bg_table_row.*/
    class_call(background_functions(pvecback_integration, pba->long_info, bg_table_row),
               error_message_,
               error_message_);

    /* Normalise D(z=0)=1 and construct f = D_prime/(aHD) */
    bg_table_row[index_bg_D_] = pvecback_integration[index_bi_D_]/D_today;
    bg_table_row[index_bg_f_] = pvecback_integration[index_bi_D_prime_]/
      (pvecback_integration[index_bi_D_]*bg_table_row[index_bg_a_]*bg_table_row[index_bg_H_]);
  }

  /** - fill tables of second derivatives (in view of spline interpolation) */
  class_call(array_spline_table_lines(z_table_,
                                      bt_size_,
                                      tau_table_,
                                      1,
                                      d2tau_dz2_table_,
                                      _SPLINE_EST_DERIV_,
                                      error_message_),
             error_message_,
             error_message_);

  class_call(array_spline_table_lines(tau_table_,
                                      bt_size_,
                                      background_table_,
                                      bg_size_,
                                      d2background_dtau2_table_,
                                      _SPLINE_EST_DERIV_,
                                      error_message_),
             error_message_,
             error_message_);

  /** - compute remaining "related parameters"
   *     - so-called "effective neutrino number", computed at earliest
      time in interpolation table. This should be seen as a
      definition: Neff is the equivalent number of
      instantaneously-decoupled neutrinos accounting for the
      radiation density, beyond photons */
  Neff_ = (background_table_[index_bg_Omega_r_]
               *background_table_[index_bg_rho_crit_]
               -background_table_[index_bg_rho_g_])
    /(7./8.*pow(4./11.,4./3.)*background_table_[index_bg_rho_g_]);

  /** - done */
  if (pba->background_verbose > 0) {
    printf(" -> age = %f Gyr\n", age_);
    printf(" -> conformal age = %f Mpc\n", conformal_age_);
  }

  if (pba->background_verbose > 2) {
    if ((pba->has_dcdm == _TRUE_)&&(pba->has_dr == _TRUE_)){
      printf("    Decaying Cold Dark Matter details: (DCDM --> DR)\n");
      printf("     -> Omega0_dcdm = %f\n", Omega0_dcdm_);
      printf("     -> Omega0_dr = %f\n", Omega0_dr_);
      printf("     -> Omega0_dr+Omega0_dcdm = %f, input value = %f\n",
              Omega0_dr_ + Omega0_dcdm_, pba->Omega0_dcdmdr);
      printf("     -> Omega_ini_dcdm/Omega_b = %f\n", pba->Omega_ini_dcdm/pba->Omega0_b);
    }
    if (pba->has_scf == _TRUE_){
      printf("    Scalar field details:\n");
      printf("     -> Omega_scf = %g, wished %g\n", pvecback[index_bg_rho_scf_]/pvecback[index_bg_rho_crit_], pba->Omega0_scf);
      if(pba->has_lambda == _TRUE_)
        printf("     -> Omega_Lambda = %g, wished %g\n", pvecback[index_bg_rho_lambda_]/pvecback[index_bg_rho_crit_], pba->Omega0_lambda);
      printf("     -> parameters: [lambda, alpha, A, B] = \n");
      printf("                    [");
      for (int i = 0; i < pba->scf_parameters.size() - 1; i++){
        printf("%.3f, ",pba->scf_parameters[i]);
      }
      printf("%.3f]\n", pba->scf_parameters[pba->scf_parameters.size() - 1]);
    }
  }

  /**  - total matter, radiation, dark energy today */
  Omega0_m_ = background_table_[(bt_size_ - 1)*bg_size_ + index_bg_Omega_m_];
  Omega0_r_ = background_table_[(bt_size_ - 1)*bg_size_ + index_bg_Omega_r_];
  Omega0_de_ = 1. - (Omega0_m_ + Omega0_r_ + pba->Omega0_k);

  free(loga);
  free(used_in_output);
  free(pvecback);
  free(pvecback_integration);

  return _SUCCESS_;

}

/**
 * Assign initial values to background integrated variables.
 *
 * @param ppr                  Input: pointer to precision structure
 * @param pba                  Input: pointer to background structure
 * @param pvecback             Input: vector of background quantities used as workspace
 * @param pvecback_integration Output: vector of background quantities to be integrated, returned with proper initial values
 * @return the error status
 */

int BackgroundModule::background_initial_conditions(double* pvecback, /* vector with argument pvecback[index_bg] (must be already allocated, normal format is sufficient) */
                                                    double* pvecback_integration /* vector with argument pvecback_integration[index_bi] (must be already allocated with size bi_size_) */
                                                    ) {

  /** Summary: */

  /** - define local variables */

  /* scale factor */
  double a;

  double rho_ncdm, p_ncdm, rho_ncdm_rel_tot=0.;
  double f,Omega_rad, rho_rad;
  int counter,is_early_enough,n_ncdm;
  double scf_lambda;
  double rho_fld_today;
  double w_fld,dw_over_da_fld,integral_fld;

  /** - fix initial value of \f$ a \f$ */
  a = ppr->a_ini_over_a_today_default * pba->a_today;

  /**  If we have ncdm species, perhaps we need to start earlier
       than the standard value for the species to be relativistic.
       This could happen for some WDM models.
  */

  if (pba->has_ncdm == _TRUE_) {

    for (counter=0; counter < _MAX_IT_; counter++) {

      is_early_enough = _TRUE_;
      rho_ncdm_rel_tot = 0.;

      for (n_ncdm=0; n_ncdm<pba->N_ncdm; n_ncdm++) {

        class_call(pba->ncdm->background_ncdm_momenta(n_ncdm, pba->a_today/a - 1.0, NULL, &rho_ncdm, &p_ncdm, NULL, NULL),
                   error_message_,
                   error_message_);
        rho_ncdm_rel_tot += 3.*p_ncdm;
        if (fabs(p_ncdm/rho_ncdm-1./3.)>ppr->tol_ncdm_initial_w)
          is_early_enough = _FALSE_;
      }
      if (is_early_enough == _TRUE_)
        break;
      else
        a *= _SCALE_BACK_;
    }
    class_test(counter == _MAX_IT_,
               error_message_,
               "Search for initial scale factor a such that all ncdm species are relativistic failed.");
  }

  pvecback_integration[index_bi_a_] = a;

  /* Set initial values of {B} variables: */
  Omega_rad = pba->Omega0_g;
  if (pba->has_ur == _TRUE_)
    Omega_rad += pba->Omega0_ur;
  if (pba->has_idr == _TRUE_)
    Omega_rad += pba->Omega0_idr;
  rho_rad = Omega_rad*pow(pba->H0,2)/pow(a/pba->a_today,4);
  if (pba->has_ncdm == _TRUE_){
    /** - We must add the relativistic contribution from NCDM species */
    rho_rad += rho_ncdm_rel_tot;
  }
  if (pba->has_dcdm == _TRUE_){
    /* Remember that the critical density today in CLASS conventions is H0^2 */
    pvecback_integration[index_bi_rho_dcdm_] =
      pba->Omega_ini_dcdm*pba->H0*pba->H0*pow(pba->a_today/a,3);
    if (pba->background_verbose > 3)
      printf("Density is %g. a_today=%g. Omega_ini=%g\n", pvecback_integration[index_bi_rho_dcdm_], pba->a_today, pba->Omega_ini_dcdm);
  }

  if (pba->has_dr == _TRUE_){
    if (pba->has_dcdm == _TRUE_){
      /**  - f is the critical density fraction of DR. The exact solution is:
       *
       * `f = -Omega_rad+pow(pow(Omega_rad,3./2.)+0.5*pow(a/pba->a_today,6)*pvecback_integration[index_bi_rho_dcdm_]*pba->Gamma_dcdm/pow(pba->H0,3),2./3.);`
       *
       * but it is not numerically stable for very small f which is always the case.
       * Instead we use the Taylor expansion of this equation, which is equivalent to
       * ignoring f(a) in the Hubble rate.
       */
      f = 1./3.*pow(a/pba->a_today,6)*pvecback_integration[index_bi_rho_dcdm_]*pba->Gamma_dcdm/pow(pba->H0,3)/sqrt(Omega_rad);
      pvecback_integration[index_bi_rho_dr_] = f*pba->H0*pba->H0/pow(a/pba->a_today,4);
    }
    else{
      /** There is also a space reserved for a future case where dr is not sourced by dcdm */
      pvecback_integration[index_bi_rho_dr_] = 0.0;
    }
  }

  if (pba->has_fld == _TRUE_){

    /* rho_fld today */
    rho_fld_today = pba->Omega0_fld * pow(pba->H0,2);

    /* integrate rho_fld(a) from a_ini to a_0, to get rho_fld(a_ini) given rho_fld(a0) */
    class_call(background_w_fld(a, &w_fld, &dw_over_da_fld, &integral_fld), error_message_, error_message_);

    /* Note: for complicated w_fld(a) functions with no simple
       analytic integral, this is the place were you should compute
       numerically the simple 1d integral [int_{a_ini}^{a_0} 3
       [(1+w_fld)/a] da] (e.g. with the Romberg method?) instead of
       calling background_w_fld */

    /* rho_fld at initial time */
    pvecback_integration[index_bi_rho_fld_] = rho_fld_today*exp(integral_fld);

  }

  /** - Fix initial value of \f$ \phi, \phi' \f$
   * set directly in the radiation attractor => fixes the units in terms of rho_ur
   *
   * TODO:
   * - There seems to be some small oscillation when it starts.
   * - Check equations and signs. Sign of phi_prime?
   * - is rho_ur all there is early on?
   */
  if(pba->has_scf == _TRUE_){
    scf_lambda = pba->scf_parameters[0];
    if(pba->attractor_ic_scf == _TRUE_){
      pvecback_integration[index_bi_phi_scf_] = -1./scf_lambda*
        log(rho_rad*4./(3*pow(scf_lambda,2)-12))*pba->phi_ini_scf;
      if (3.*pow(scf_lambda,2)-12. < 0){
        /** - --> If there is no attractor solution for scf_lambda, assign some value. Otherwise would give a nan.*/
    	pvecback_integration[index_bi_phi_scf_] = 1./scf_lambda;//seems to the work
        if (pba->background_verbose > 0)
          printf(" No attractor IC for lambda = %.3e ! \n ",scf_lambda);
      }
      pvecback_integration[index_bi_phi_prime_scf_] = 2*pvecback_integration[index_bi_a_]*
        sqrt(V_scf(pvecback_integration[index_bi_phi_scf_]))*pba->phi_prime_ini_scf;
    }
    else{
      printf("Not using attractor initial conditions\n");
      /** - --> If no attractor initial conditions are assigned, gets the provided ones. */
      pvecback_integration[index_bi_phi_scf_] = pba->phi_ini_scf;
      pvecback_integration[index_bi_phi_prime_scf_] = pba->phi_prime_ini_scf;
    }
    class_test(!isfinite(pvecback_integration[index_bi_phi_scf_]) ||
               !isfinite(pvecback_integration[index_bi_phi_scf_]),
               error_message_,
               "initial phi = %e phi_prime = %e -> check initial conditions",
               pvecback_integration[index_bi_phi_scf_],
               pvecback_integration[index_bi_phi_scf_]);
  }

  /* Infer pvecback from pvecback_integration */
  class_call(background_functions(pvecback_integration, pba->normal_info, pvecback),
             error_message_,
             error_message_);

  /* Just checking that our initial time indeed is deep enough in the radiation
     dominated regime */
  class_test(fabs(pvecback[index_bg_Omega_r_] - 1.) > ppr->tol_initial_Omega_r,
             error_message_,
             "Omega_r = %e, not close enough to 1. Decrease a_ini_over_a_today_default in order to start from radiation domination.",
             pvecback[index_bg_Omega_r_]);

  /** - compute initial proper time, assuming radiation-dominated
      universe since Big Bang and therefore \f$ t=1/(2H) \f$ (good
      approximation for most purposes) */

  class_test(pvecback[index_bg_H_] <= 0.,
             error_message_,
             "H = %e instead of strictly positive",pvecback[index_bg_H_]);

  pvecback_integration[index_bi_time_] = 1./(2.* pvecback[index_bg_H_]);

  /** - compute initial conformal time, assuming radiation-dominated
      universe since Big Bang and therefore \f$ \tau=1/(aH) \f$
      (good approximation for most purposes) */
  pvecback_integration[index_bi_tau_] = 1./(a*pvecback[index_bg_H_]);

  /** - compute initial sound horizon, assuming \f$ c_s=1/\sqrt{3} \f$ initially */
  pvecback_integration[index_bi_rs_] = pvecback_integration[index_bi_tau_]/sqrt(3.);

  /** - set initial value of D and D' in RD. D will be renormalised later, but D' must be correct. */
  pvecback_integration[index_bi_D_] = a;
  pvecback_integration[index_bi_D_prime_] = 2*pvecback_integration[index_bi_D_]*pvecback[index_bg_H_];

  return _SUCCESS_;

}

/**
 * Find the time of radiation/matter equality and store characteristic
 * quantitites at that time in the background structure..
 *
 * @param ppr                  Input: pointer to precision structure
 * @param pba                  Input/Output: pointer to background structure
 * @return the error status
 */

int BackgroundModule::background_find_equality() {

  double Omega_m_over_Omega_r=0.;
  int index_tau_minus = 0;
  int index_tau_plus = bt_size_ - 1;
  int index_tau_mid = 0;
  double tau_minus,tau_plus,tau_mid=0.;
  double * pvecback;

  /* first bracket the right tau value between two consecutive indices in the table */

  while ((index_tau_plus - index_tau_minus) > 1) {

    index_tau_mid = (int)(0.5*(index_tau_plus+index_tau_minus));

    Omega_m_over_Omega_r = background_table_[index_tau_mid*bg_size_ + index_bg_Omega_m_]/background_table_[index_tau_mid*bg_size_ + index_bg_Omega_r_];

    if (Omega_m_over_Omega_r > 1)
      index_tau_plus = index_tau_mid;
    else
      index_tau_minus = index_tau_mid;

  }

  /* then get a better estimate within this range */

  tau_minus = tau_table_[index_tau_minus];
  tau_plus =  tau_table_[index_tau_plus];

  class_alloc(pvecback, bg_size_*sizeof(double), error_message_);

  while ((tau_plus - tau_minus) > ppr->tol_tau_eq) {

    tau_mid = 0.5*(tau_plus+tau_minus);

    class_call(background_at_tau(tau_mid, pba->long_info, pba->inter_closeby, &index_tau_minus, pvecback),
               error_message_,
               error_message_);

    Omega_m_over_Omega_r = pvecback[index_bg_Omega_m_]/pvecback[index_bg_Omega_r_];

    if (Omega_m_over_Omega_r > 1)
      tau_plus = tau_mid;
    else
      tau_minus = tau_mid;

  }

  a_eq_ = pvecback[index_bg_a_];
  H_eq_ = pvecback[index_bg_H_];
  z_eq_ = pba->a_today/a_eq_ - 1.;
  tau_eq_ = tau_mid;

  if (pba->background_verbose > 0) {
    printf(" -> radiation/matter equality at z = %f\n", z_eq_);
    printf("    corresponding to conformal time = %f Mpc\n", tau_eq_);
  }

  free(pvecback);

  return _SUCCESS_;

}


/**
 * Subroutine for formatting background output
 *
 */

int BackgroundModule::background_output_titles(char titles[_MAXTITLESTRINGLENGTH_]) const {

  /** - Length of the column title should be less than _OUTPUTPRECISION_+6
      to be indented correctly, but it can be as long as . */
  int n;
  char tmp[40];

  class_store_columntitle(titles,"z",_TRUE_);
  class_store_columntitle(titles,"proper time [Gyr]",_TRUE_);
  class_store_columntitle(titles,"conf. time [Mpc]",_TRUE_);
  class_store_columntitle(titles,"H [1/Mpc]",_TRUE_);
  class_store_columntitle(titles,"comov. dist.",_TRUE_);
  class_store_columntitle(titles,"ang.diam.dist.",_TRUE_);
  class_store_columntitle(titles,"lum. dist.",_TRUE_);
  class_store_columntitle(titles,"comov.snd.hrz.",_TRUE_);
  class_store_columntitle(titles,"(.)rho_g",_TRUE_);
  class_store_columntitle(titles,"(.)rho_b",_TRUE_);
  class_store_columntitle(titles,"(.)rho_cdm",pba->has_cdm);
  if (pba->has_ncdm == _TRUE_){
    for (n=0; n<pba->N_ncdm; n++){
      sprintf(tmp,"(.)rho_ncdm[%d]",n);
      class_store_columntitle(titles,tmp,_TRUE_);
      sprintf(tmp,"(.)p_ncdm[%d]",n);
      class_store_columntitle(titles,tmp,_TRUE_);
    }
  }
  class_store_columntitle(titles,"(.)rho_lambda",pba->has_lambda);
  class_store_columntitle(titles,"(.)rho_fld",pba->has_fld);
  class_store_columntitle(titles,"(.)w_fld",pba->has_fld);
  class_store_columntitle(titles,"(.)rho_ur",pba->has_ur);
  class_store_columntitle(titles,"(.)rho_idr",pba->has_idr);
  class_store_columntitle(titles,"(.)rho_idm_dr",pba->has_idm_dr);
  class_store_columntitle(titles,"(.)rho_crit",_TRUE_);
  class_store_columntitle(titles,"(.)rho_dcdm",pba->has_dcdm);
  class_store_columntitle(titles,"(.)rho_dr",pba->has_dr);

  class_store_columntitle(titles,"(.)rho_scf",pba->has_scf);
  class_store_columntitle(titles,"(.)p_scf",pba->has_scf);
  class_store_columntitle(titles,"(.)p_prime_scf",pba->has_scf);
  class_store_columntitle(titles,"phi_scf",pba->has_scf);
  class_store_columntitle(titles,"phi'_scf",pba->has_scf);
  class_store_columntitle(titles,"V_scf",pba->has_scf);
  class_store_columntitle(titles,"V'_scf",pba->has_scf);
  class_store_columntitle(titles,"V''_scf",pba->has_scf);

  class_store_columntitle(titles,"(.)rho_tot",_TRUE_);
  class_store_columntitle(titles,"(.)p_tot",_TRUE_);
  class_store_columntitle(titles,"(.)p_tot_prime",_TRUE_);

  class_store_columntitle(titles,"gr.fac. D",_TRUE_);
  class_store_columntitle(titles,"gr.fac. f",_TRUE_);

  return _SUCCESS_;
}

int BackgroundModule::background_output_data(int number_of_titles, double* data) const {
  int index_tau, storeidx, n;
  double *dataptr, *pvecback;

  /** Stores quantities */
  for (index_tau = 0; index_tau < bt_size_; index_tau++){
    dataptr = data + index_tau*number_of_titles;
    pvecback = background_table_ + index_tau*bg_size_;
    storeidx = 0;

    class_store_double(dataptr, pba->a_today/pvecback[index_bg_a_]-1., _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_time_]/_Gyr_over_Mpc_, _TRUE_, storeidx);
    class_store_double(dataptr, conformal_age_-pvecback[index_bg_conf_distance_], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_H_], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_conf_distance_], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_ang_distance_], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_lum_distance_], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_rs_], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_rho_g_], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_rho_b_], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_rho_cdm_], pba->has_cdm, storeidx);
    if (pba->has_ncdm == _TRUE_){
      for (n=0; n<pba->N_ncdm; n++){
        class_store_double(dataptr, pvecback[index_bg_rho_ncdm1_+n], _TRUE_, storeidx);
        class_store_double(dataptr, pvecback[index_bg_p_ncdm1_+n], _TRUE_, storeidx);
      }
    }
    class_store_double(dataptr, pvecback[index_bg_rho_lambda_], pba->has_lambda, storeidx);
    class_store_double(dataptr, pvecback[index_bg_rho_fld_], pba->has_fld, storeidx);
    class_store_double(dataptr, pvecback[index_bg_w_fld_], pba->has_fld, storeidx);
    class_store_double(dataptr, pvecback[index_bg_rho_ur_], pba->has_ur, storeidx);
    class_store_double(dataptr, pvecback[index_bg_rho_idr_], pba->has_idr, storeidx);
    class_store_double(dataptr, pvecback[index_bg_rho_idm_dr_], pba->has_idm_dr, storeidx);
    class_store_double(dataptr, pvecback[index_bg_rho_crit_], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_rho_dcdm_], pba->has_dcdm, storeidx);
    class_store_double(dataptr, pvecback[index_bg_rho_dr_], pba->has_dr, storeidx);

    class_store_double(dataptr, pvecback[index_bg_rho_scf_], pba->has_scf, storeidx);
    class_store_double(dataptr, pvecback[index_bg_p_scf_], pba->has_scf, storeidx);
    class_store_double(dataptr, pvecback[index_bg_p_prime_scf_], pba->has_scf, storeidx);
    class_store_double(dataptr, pvecback[index_bg_phi_scf_], pba->has_scf, storeidx);
    class_store_double(dataptr, pvecback[index_bg_phi_prime_scf_], pba->has_scf, storeidx);
    class_store_double(dataptr, pvecback[index_bg_V_scf_], pba->has_scf, storeidx);
    class_store_double(dataptr, pvecback[index_bg_dV_scf_], pba->has_scf, storeidx);
    class_store_double(dataptr, pvecback[index_bg_ddV_scf_], pba->has_scf, storeidx);

    class_store_double(dataptr, pvecback[index_bg_rho_tot_], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_p_tot_], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_p_tot_prime_], _TRUE_, storeidx);

    class_store_double(dataptr, pvecback[index_bg_D_], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[index_bg_f_], _TRUE_, storeidx);
  }

  return _SUCCESS_;
}


/**
 * Subroutine evaluating the derivative with respect to conformal time
 * of quantities which are integrated (a, t, etc).
 *
 * This is one of the few functions in the code which is passed to
 * the generic_integrator() routine.  Since generic_integrator()
 * should work with functions passed from various modules, the format
 * of the arguments is a bit special:
 *
 * - fixed input parameters and workspaces are passed through a generic
 * pointer. Here, this is just a pointer to the background structure
 * and to a background vector, but generic_integrator() doesn't know
 * its fine structure.
 *
 * - the error management is a bit special: errors are not written as
 * usual to error_message_, but to a generic error_message passed
 * in the list of arguments.
 *
 * @param tau                      Input: conformal time
 * @param y                        Input: vector of variable
 * @param dy                       Output: its derivative (already allocated)
 * @param parameters_and_workspace Input: pointer to fixed parameters (e.g. indices)
 * @param error_message            Output: error message
 */
int BackgroundModule::background_derivs_member(
                                               double tau,
                                               double* y, /* vector with argument y[index_bi] (must be already allocated with size bi_size_) */
                                               double* dy, /* vector with argument dy[index_bi]
                                                            (must be already allocated with
                                                            size bi_size_) */
                                               void* parameters_and_workspace,
                                               ErrorMsg error_message
                                               ) {

  /** Summary: */

  /** - define local variables */

  double * pvecback, a, H, rho_M;

  background_parameters_and_workspace* pbpaw = static_cast<background_parameters_and_workspace*>(parameters_and_workspace);
  pvecback = pbpaw->pvecback;

  /** - calculate functions of \f$ a \f$ with background_functions() */
  class_call(background_functions(y, pba->normal_info, pvecback),
             error_message_,
             error_message);

  /** - Short hand notation */
  a = y[index_bi_a_];
  H = pvecback[index_bg_H_];

  /** - calculate \f$ a'=a^2 H \f$ */
  dy[index_bi_a_] = y[index_bi_a_]*y[index_bi_a_]*pvecback[index_bg_H_];

  /** - calculate \f$ t' = a \f$ */
  dy[index_bi_time_] = y[index_bi_a_];

  class_test(pvecback[index_bg_rho_g_] <= 0.,
             error_message,
             "rho_g = %e instead of strictly positive", pvecback[index_bg_rho_g_]);

  /** - calculate \f$ rs' = c_s \f$*/
  dy[index_bi_rs_] = 1./sqrt(3.*(1. + 3.*pvecback[index_bg_rho_b_]/4./pvecback[index_bg_rho_g_]))*sqrt(1. - pba->K*y[index_bi_rs_]*y[index_bi_rs_]); // TBC: curvature correction

  /** - solve second order growth equation  \f$ [D''(\tau)=-aHD'(\tau)+3/2 a^2 \rho_M D(\tau) \f$ */
  rho_M = pvecback[index_bg_rho_b_];
  if (pba->has_cdm)
    rho_M += pvecback[index_bg_rho_cdm_];
  if (pba->has_idm_dr)
    rho_M += pvecback[index_bg_rho_idm_dr_];

  dy[index_bi_D_] = y[index_bi_D_prime_];
  dy[index_bi_D_prime_] = -a*H*y[index_bi_D_prime_] + 1.5*a*a*rho_M*y[index_bi_D_];

  if (pba->has_dcdm == _TRUE_){
    /** - compute dcdm density \f$ \rho' = -3aH \rho - a \Gamma \rho \f$*/
    dy[index_bi_rho_dcdm_] = -3.*y[index_bi_a_]*pvecback[index_bg_H_]*y[index_bi_rho_dcdm_]-
      y[index_bi_a_]*pba->Gamma_dcdm*y[index_bi_rho_dcdm_];
  }

  if ((pba->has_dcdm == _TRUE_) && (pba->has_dr == _TRUE_)){
    /** - Compute dr density \f$ \rho' = -4aH \rho - a \Gamma \rho \f$ */
    dy[index_bi_rho_dr_] = -4.*y[index_bi_a_]*pvecback[index_bg_H_]*y[index_bi_rho_dr_]+
      y[index_bi_a_]*pba->Gamma_dcdm*y[index_bi_rho_dcdm_];
  }

  if (pba->has_fld == _TRUE_) {
    /** - Compute fld density \f$ \rho' = -3aH (1+w_{fld}(a)) \rho \f$ */
    dy[index_bi_rho_fld_] = -3.*y[index_bi_a_]*pvecback[index_bg_H_]*(1. + pvecback[index_bg_w_fld_])*y[index_bi_rho_fld_];
  }

  if (pba->has_scf == _TRUE_){
    /** - Scalar field equation: \f$ \phi'' + 2 a H \phi' + a^2 dV = 0 \f$  (note H is wrt cosmic time) */
    dy[index_bi_phi_scf_] = y[index_bi_phi_prime_scf_];
    dy[index_bi_phi_prime_scf_] = -y[index_bi_a_]*(2*pvecback[index_bg_H_]*y[index_bi_phi_prime_scf_] + y[index_bi_a_]*dV_scf(y[index_bi_phi_scf_]));
  }

  return _SUCCESS_;

}

/**
 * Scalar field potential and its derivatives with respect to the field _scf
 * For Albrecht & Skordis model: 9908085
 * - \f$ V = V_{p_{scf}}*V_{e_{scf}} \f$
 * - \f$ V_e =  \exp(-\lambda \phi) \f$ (exponential)
 * - \f$ V_p = (\phi - B)^\alpha + A \f$ (polynomial bump)
 *
 * TODO:
 * - Add some functionality to include different models/potentials (tuning would be difficult, though)
 * - Generalize to Kessence/Horndeski/PPF and/or couplings
 * - A default module to numerically compute the derivatives when no analytic functions are given should be added.
 * - Numerical derivatives may further serve as a consistency check.
 *
 */

/**
 *
 * The units of phi, tau in the derivatives and the potential V are the following:
 * - phi is given in units of the reduced Planck mass \f$ m_{pl} = (8 \pi G)^{(-1/2)}\f$
 * - tau in the derivative is given in units of Mpc.
 * - the potential \f$ V(\phi) \f$ is given in units of \f$ m_{pl}^2/Mpc^2 \f$.
 * With this convention, we have
 * \f$ \rho^{class} = (8 \pi G)/3 \rho^{physical} = 1/(3 m_{pl}^2) \rho^{physical} = 1/3 * [ 1/(2a^2) (\phi')^2 + V(\phi) ] \f$
 and \f$ \rho^{class} \f$ has the proper dimension \f$ Mpc^-2 \f$.
*/

double BackgroundModule::V_e_scf(double phi) const {
  double scf_lambda = pba->scf_parameters[0];
  //  double scf_alpha  = pba->scf_parameters[1];
  //  double scf_A      = pba->scf_parameters[2];
  //  double scf_B      = pba->scf_parameters[3];

  return  exp(-scf_lambda*phi);
}

double BackgroundModule::dV_e_scf(double phi) const {
  double scf_lambda = pba->scf_parameters[0];
  //  double scf_alpha  = pba->scf_parameters[1];
  //  double scf_A      = pba->scf_parameters[2];
  //  double scf_B      = pba->scf_parameters[3];

  return -scf_lambda*V_scf(phi);
}

double BackgroundModule::ddV_e_scf(double phi) const {
  double scf_lambda = pba->scf_parameters[0];
  //  double scf_alpha  = pba->scf_parameters[1];
  //  double scf_A      = pba->scf_parameters[2];
  //  double scf_B      = pba->scf_parameters[3];

  return pow(-scf_lambda,2)*V_scf(phi);
}


/** parameters and functions for the polynomial coefficient
 * \f$ V_p = (\phi - B)^\alpha + A \f$(polynomial bump)
 *
 * double scf_alpha = 2;
 *
 * double scf_B = 34.8;
 *
 * double scf_A = 0.01; (values for their Figure 2)
 */

double BackgroundModule::V_p_scf(double phi) const {
  //  double scf_lambda = pba->scf_parameters[0];
  double scf_alpha  = pba->scf_parameters[1];
  double scf_A      = pba->scf_parameters[2];
  double scf_B      = pba->scf_parameters[3];

  return  pow(phi - scf_B,  scf_alpha) +  scf_A;
}

double BackgroundModule::dV_p_scf(double phi) const {

  //  double scf_lambda = pba->scf_parameters[0];
  double scf_alpha  = pba->scf_parameters[1];
  //  double scf_A      = pba->scf_parameters[2];
  double scf_B      = pba->scf_parameters[3];

  return   scf_alpha*pow(phi -  scf_B,  scf_alpha - 1);
}

double BackgroundModule::ddV_p_scf(double phi) const {
  //  double scf_lambda = pba->scf_parameters[0];
  double scf_alpha  = pba->scf_parameters[1];
  //  double scf_A      = pba->scf_parameters[2];
  double scf_B      = pba->scf_parameters[3];

  return  scf_alpha*(scf_alpha - 1.)*pow(phi -  scf_B,  scf_alpha - 2);
}

/** Fianlly we can obtain the overall potential \f$ V = V_p*V_e \f$
 */

double BackgroundModule::V_scf(double phi) const {
  return  V_e_scf(phi)*V_p_scf(phi);
}

double BackgroundModule::dV_scf(double phi) const {
  return dV_e_scf(phi)*V_p_scf(phi) + V_e_scf(phi)*dV_p_scf(phi);
}

double BackgroundModule::ddV_scf(double phi) const {
  return ddV_e_scf(phi)*V_p_scf(phi) + 2*dV_e_scf(phi)*dV_p_scf(phi) + V_e_scf(phi)*ddV_p_scf(phi);
}

/**
 * Function outputting the fractions Omega of the total critical density
 * today, and also the reduced fractions omega=Omega*h*h
 *
 * It also prints the total budgets of non-relativistic, relativistic,
 * and other contents, and of the total
 *
 * @param pba                      Input: Pointer to background structure
 * @return the error status
 */

int BackgroundModule::background_output_budget() {

  double budget_matter, budget_radiation, budget_other,budget_neutrino;

  budget_matter = 0;
  budget_radiation = 0;
  budget_other = 0;
  budget_neutrino = 0;

  //The name for the _class_print_species_ macro can be at most 30 characters total
  if(pba->background_verbose > 1){

    printf(" ---------------------------- Budget equation ----------------------- \n");

    printf(" ---> Nonrelativistic Species \n");
    _class_print_species_("Bayrons",b);
    budget_matter+=pba->Omega0_b;
    if(pba->has_cdm){
      _class_print_species_("Cold Dark Matter",cdm);
      budget_matter+=pba->Omega0_cdm;
    }
    if(pba->has_idm_dr){
      _class_print_species_("Interacting Dark Matter - DR ",idm_dr);
      budget_matter+=pba->Omega0_idm_dr;
    }
    if(pba->has_dcdm){
      printf("-> %-30s Omega = %-15g , omega = %-15g\n", "Decaying Cold Dark Matter", Omega0_dcdm_, Omega0_dcdm_*pba->h*pba->h);
      budget_matter += Omega0_dcdm_;
    }


    printf(" ---> Relativistic Species \n");
    _class_print_species_("Photons",g);
    budget_radiation+=pba->Omega0_g;
    if(pba->has_ur){
      _class_print_species_("Ultra-relativistic relics",ur);
      budget_radiation+=pba->Omega0_ur;
    }
    if(pba->has_dr){
      printf("-> %-30s Omega = %-15g , omega = %-15g\n", "Dark Radiation (from decay)", Omega0_dr_, Omega0_dr_*pba->h*pba->h);
      budget_radiation += Omega0_dr_;
    }
    if(pba->has_idr){
      _class_print_species_("Interacting Dark Radiation",idr);
      budget_radiation+=pba->Omega0_idr;
    }

    if(pba->N_ncdm > 0){
      printf(" ---> Massive Neutrino Species \n");
    }
    if(pba->N_ncdm > 0){
      pba->ncdm->PrintOmegaInfo();
      budget_neutrino += pba->ncdm->GetOmega0();
    }

    if(pba->has_lambda || pba->has_fld || pba->has_scf || pba->has_curvature){
      printf(" ---> Other Content \n");
    }
    if(pba->has_lambda){
      _class_print_species_("Cosmological Constant",lambda);
      budget_other+=pba->Omega0_lambda;
    }
    if(pba->has_fld){
      _class_print_species_("Dark Energy Fluid",fld);
      budget_other+=pba->Omega0_fld;
    }
    if(pba->has_scf){
      _class_print_species_("Scalar Field",scf);
      budget_other+=pba->Omega0_scf;
    }
    if(pba->has_curvature){
      _class_print_species_("Spatial Curvature",k);
      budget_other+=pba->Omega0_k;
    }

    printf(" ---> Total budgets \n");
    printf(" Radiation                        Omega = %-15g , omega = %-15g \n",budget_radiation,budget_radiation*pba->h*pba->h);
    printf(" Non-relativistic                 Omega = %-15g , omega = %-15g \n",budget_matter,budget_matter*pba->h*pba->h);
    if(pba->N_ncdm > 0){
      printf(" Neutrinos                        Omega = %-15g , omega = %-15g \n",budget_neutrino,budget_neutrino*pba->h*pba->h);
    }
    if(pba->has_lambda || pba->has_fld || pba->has_scf || pba->has_curvature){
      printf(" Other Content                    Omega = %-15g , omega = %-15g \n",budget_other,budget_other*pba->h*pba->h);
    }
    printf(" TOTAL                            Omega = %-15g , omega = %-15g \n",budget_radiation+budget_matter+budget_neutrino+budget_other,(budget_radiation+budget_matter+budget_neutrino+budget_other)*pba->h*pba->h);

    printf(" -------------------------------------------------------------------- \n");
  }

  return _SUCCESS_;
}

/**
 * Subroutine evaluating the derivative with respect to log(a)
 * of quantities which are integrated (tau, t, etc).
 *
 * This is one of the few functions in the code which is passed to
 * the generic_integrator() routine.  Since generic_integrator()
 * should work with functions passed from various modules, the format
 * of the arguments is a bit special:
 *
 * - fixed input parameters and workspaces are passed through a generic
 * pointer. Here, this is just a pointer to the background structure
 * and to a background vector, but generic_integrator() doesn't know
 * its fine structure.
 *
 * - the error management is a bit special: errors are not written as
 * usual to pba->error_message, but to a generic error_message passed
 * in the list of arguments.
 *
 * @param loga                        Input: scale factor
 * @param y                        Input: vector of variable
 * @param dy                       Output: its derivative (already allocated)
 * @param parameters_and_workspace Input: pointer to fixed parameters (e.g. indices)
 * @param error_message            Output: error message
 */
int BackgroundModule::background_derivs_loga_member(double loga,
                                                    double* y, /* vector with argument y[index_bi] (must be already allocated with size bi_size_) */
                                                    double* dy, /* vector with argument dy[index_bi]
                                                                 (must be already allocated with
                                                                 size bi_size_) */
                                                    void* parameters_and_workspace,
                                                    ErrorMsg error_message
                                                    ) {

  double a, tau;
  double H;
  int index_bi;

  background_parameters_and_workspace* pbpaw = static_cast<background_parameters_and_workspace*>(parameters_and_workspace);
  double* pvecback = pbpaw->pvecback;

  /** Note that we want to reuse as much as possible of the usual tau
      integration, so inside y, index_bi_a is really the index of tau. */
  a = exp(loga);
  tau = y[index_bi_a_];
  y[index_bi_a_] = a;

  /** Get derivatives w.r.t. conformal time */
  class_call(background_derivs_member(tau, y, dy, parameters_and_workspace, error_message),
       error_message_,
       error_message_);

  /** Swap a and tau again */
  y[index_bi_a_] = tau;
  /** Set dtau/dloga: */
  dy[index_bi_a_] = 1.0;

  /** Convert to scale factor derivatives */
  H = pvecback[index_bg_H_];
  for (index_bi = 0; index_bi < bi_size_ - 1; index_bi++) {
    dy[index_bi] *= 1./(a*H);
  }
  return _SUCCESS_;
}

int BackgroundModule::background_add_line_to_bg_table_member(
                                                             double loga,
                                                             double* y,
                                                             double* dy,
                                                             int index_loga,
                                                             void* parameters_and_workspace,
                                                             ErrorMsg error_message) {

  z_table_[index_loga] = MAX(0., pba->a_today/exp(loga) - 1.);
  tau_table_[index_loga] = y[index_bi_a_]; /* tau at a's spot..*/

  /* -> write in the table. We are only storing the integration vector at this point. */
  for (int index_bi = 0; index_bi < bi_size_ - 1; index_bi++) {
    background_table_[index_loga*bg_size_ + index_bi] = y[index_bi];
  }

  /* tau is evolved at a's spot...*/
  background_table_[index_loga*bg_size_ + index_bi_tau_] = y[index_bi_a_];
  background_table_[index_loga*bg_size_ + index_bi_a_] = pba->a_today*1./(1. + z_table_[index_loga]);
  return _SUCCESS_;
}
