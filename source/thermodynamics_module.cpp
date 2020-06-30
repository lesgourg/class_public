/** @file thermodynamics.c Documented thermodynamics module
 *
 * Julien Lesgourgues, 6.09.2010
 *
 * Deals with the thermodynamical evolution.
 * This module has two purposes:
 *
 * - at the beginning, to initialize the thermodynamics, i.e. to
 *   integrate the thermodynamical equations, and store all
 *   thermodynamical quantities as a function of redshift inside an
 *   interpolation table. The current version of recombination is
 *   based on RECFAST v1.5. The current version of reionization is
 *   based on exactly the same reionization function as in CAMB, in
 *   order to make allow for comparison. It should be easy to
 *   generalize the module to more complicated reionization histories.
 *
 * - to provide a routine which allow other modules to evaluate any
 *   thermodynamical quantities at a given redshift value (by
 *   interpolating within the interpolation table).
 *
 *
 * The logic is the following:
 *
 * - in a first step, the code assumes that there is no reionization,
 *   and computes the ionization fraction, Thomson scattering rate,
 *   baryon temperature, etc., using RECFAST. The result is stored in
 *   a temporary table 'recombination_table' (within a temporary
 *   structure of type 'recombination') for each redshift in a range 0
 *   < z < z_initial.  The sampling in z space is done with a simple
 *   linear step size.
 * - in a second step, the code adds the reionization history,
 *   starting from a redshift z_reio_start. The ionization fraction at
 *   this redshift is read in the previous recombination table in
 *   order to ensure a perfect matching. The code computes the
 *   ionization fraction, Thomson scattering rate, baryon temperature,
 *   etc., using a given parametrization of the reionization
 *   history. The result is stored in a temporary table
 *   'reionization_table' (within a temporary structure of type
 *   'reionization') for each redshift in the range 0 < z <
 *   z_reio_start. The sampling in z space is found automatically,
 *   given the precision parameter 'reionization_sampling'.
 *
 * - in a third step, the code merges the two tables
 *   'recombination_table' and 'reionization_table' inside the table
 *   'thermodynamics_table', and the temporary structures
 *   'recombination' and 'reionization' are freed. In
 *   'thermodynamics_table', the sampling in z space is the one
 *   defined in the recombination algorithm for z_reio_start < z <
 *   z_initial, and the one defined in the reionization algorithm for
 *   0 < z < z_reio_start.
 *
 * - at this stage, only a few columns in the table
 *   'thermodynamics_table' have been filled. In a fourth step, the
 *   remaining columns are filled, using some numerical
 *   integration/derivation routines from the 'array.c' tools module.
 *
 * - small detail: one of the columns contains the maximum variation
 *   rate of a few relevant thermodynamical quantities. This rate
 *   will be used for defining automatically the sampling step size in
 *   the perturbation module. Hence, the exact value of this rate is
 *   unimportant, but its order of magnitude at a given z defines the
 *   sampling precision of the perturbation module. Hence, it is
 *   harmless to use a smoothing routine in order to make this rate
 *   look nicer, although this will not affect the final result
 *   significantly. The last step in the thermodynamics_init module is
 *   to perform this smoothing.
 *
 * In summary, the following functions can be called from other modules:
 *
 * -# thermodynamics_init() at the beginning (but after background_init())
 * -# thermodynamics_at_z() at any later time
 * -# thermodynamics_free() at the end, when no more calls to thermodynamics_at_z() are needed
 */

#include "background_module.h"
#include "thermodynamics_module.h"

#ifdef HYREC
#include "hyrec.h"
#endif

ThermodynamicsModule::ThermodynamicsModule(InputModulePtr input_module, BackgroundModulePtr background_module)
: BaseModule(std::move(input_module))
, background_module_(std::move(background_module)) {
  if (thermodynamics_init() != _SUCCESS_) {
    throw std::runtime_error(error_message_);
  }
}

ThermodynamicsModule::~ThermodynamicsModule() {
  thermodynamics_free();
}

// Wrapper functions to pass non-static member functions
int ThermodynamicsModule::thermodynamics_derivs_with_recfast(double z, double* y, double* dy, void* fixed_parameters, ErrorMsg error_message) {
  auto tppaw = static_cast<thermodynamics_parameters_and_workspace*>(fixed_parameters);
  return tppaw->thermodynamics_module->thermodynamics_derivs_with_recfast_member(z, y, dy, fixed_parameters, error_message);
}

/**
 * Thermodynamics quantities at given redshift z.
 *
 * Evaluates all thermodynamics quantities at a given value of
 * the redshift by reading the pre-computed table and interpolating.
 *
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to the thermodynamics structure (containing pre-computed table)
 * @param z          Input: redshift
 * @param inter_mode Input: interpolation mode (normal or growing_closeby)
 * @param last_index Input/Output: index of the previous/current point in the interpolation array (input only for closeby mode, output for both)
 * @param pvecback   Input: vector of background quantities (used only in case z>z_initial for getting ddkappa and dddkappa; in that case, should be already allocated and filled, with format short_info or larger; in other cases, will be ignored)
 * @param pvecthermo Output: vector of thermodynamics quantities (assumed to be already allocated)
 * @return the error status
 */

int ThermodynamicsModule::thermodynamics_at_z(double z, short inter_mode, int* last_index, double* pvecback, double* pvecthermo) const {

  /** Summary: */

  /** - define local variables */

  double x0;

  /* - the fact that z is in the pre-computed range 0 <= z <= z_initial
     will be checked in the interpolation routines below. Before
     trying to interpolate, allow the routine to deal with the case z
     > z_intial: then, all relevant quantities can be extrapolated
     using simple analytic approximations */

  if (z >= z_table_[tt_size_ - 1]) {

    /* ionization fraction assumed to remain constant at large z */
    x0 = thermodynamics_table_[(tt_size_ - 1)*th_size_ + index_th_xe_];
    pvecthermo[index_th_xe_] = x0;

    /* Calculate dkappa/dtau (dkappa/dtau = a n_e x_e sigma_T = a^{-2} n_e(today) x_e sigma_T in units of 1/Mpc) */
    pvecthermo[index_th_dkappa_] = (1. + z)*(1. + z)*n_e_*x0*_sigma_*_Mpc_over_m_;

    /* tau_d scales like (1+z)**2 */
    pvecthermo[index_th_tau_d_] = thermodynamics_table_[(tt_size_ - 1)*th_size_ + index_th_tau_d_]*pow((1 + z)/(1. + z_table_[tt_size_ - 1]), 2);

    if (pth->compute_damping_scale == _TRUE_) {

      /* r_d scales like (1+z)**-3/2 */
      pvecthermo[index_th_r_d_] = thermodynamics_table_[(tt_size_ - 1)*th_size_ + index_th_r_d_]*pow((1 + z)/(1. + z_table_[tt_size_ - 1]), -1.5);

    }

    /* Calculate d2kappa/dtau2 = dz/dtau d/dz[dkappa/dtau] given that [dkappa/dtau] proportional to (1+z)^2 and dz/dtau = -H */
    pvecthermo[index_th_ddkappa_] = -pvecback[background_module_->index_bg_H_]*2./(1. + z)*pvecthermo[index_th_dkappa_];

    /* Calculate d3kappa/dtau3 given that [dkappa/dtau] proportional to (1+z)^2 */
    pvecthermo[index_th_dddkappa_] = (pvecback[background_module_->index_bg_H_]*pvecback[background_module_->index_bg_H_]/(1. + z) - pvecback[background_module_->index_bg_H_prime_])*2./(1. + z)*pvecthermo[index_th_dkappa_];

    /* \f$ exp^{-\kappa}, g, g', g'' \f$ can be set to zero: they are
       used only for computing the source functions in the
       perturbation module; but source functions only need to be
       sampled below z_initial (the condition that
       z_start_sources<z_initial is checked in the perturbation
       module) */
    pvecthermo[index_th_exp_m_kappa_] = 0.;
    pvecthermo[index_th_g_] = 0.;
    pvecthermo[index_th_dg_] = 0.;
    pvecthermo[index_th_ddg_] = 0.;

    /* Calculate Tb */
    pvecthermo[index_th_Tb_] = pba->T_cmb*(1. + z);

    /* Calculate baryon equation of state parameter wb = (k_B/mu) Tb */
    /* note that m_H / mu = 1 + (m_H/m_He-1) Y_p + x_e (1-Y_p) */
    pvecthermo[index_th_wb_] = _k_B_/(_c_*_c_*_m_H_)*(1. + (1./_not4_ - 1.)*YHe_ + x0*(1. - YHe_))*pba->T_cmb*(1. + z);

    /* Calculate baryon adiabatic sound speed cb2 = (k_B/mu) Tb (1-1/3 dlnTb/dlna) = (k_B/mu) Tb (1+1/3 (1+z) dlnTb/dz) */
    /* note that m_H / mu = 1 + (m_H/m_He-1) Y_p + x_e (1-Y_p) */
    pvecthermo[index_th_cb2_] = pvecthermo[index_th_wb_]*4./3.;

    /* derivatives of baryon sound speed (only computed if some non-minimal tight-coupling schemes is requested) */
    if (pth->compute_cb2_derivatives == _TRUE_) {

      /* since cb2 proportional to (1+z) or 1/a, its derivative wrt conformal time is given by dcb2 = - a H cb2 */
      pvecthermo[index_th_dcb2_] = -pvecback[background_module_->index_bg_H_]*pvecback[background_module_->index_bg_a_]*pvecthermo[index_th_cb2_];

      /* then its second derivative is given by ddcb2 = - a H' cb2 */
      pvecthermo[index_th_ddcb2_] = -pvecback[background_module_->index_bg_H_prime_]*pvecback[background_module_->index_bg_a_]*pvecthermo[index_th_cb2_];
    }

    /* in this regime, variation rate = dkappa/dtau */
    pvecthermo[index_th_rate_] = pvecthermo[index_th_dkappa_];

    /* quantities related to DM interacting with DR */
    if(pba->has_idm_dr == _TRUE_){

      /* calculate dmu_idm_dr and approximate its derivatives as zero */
      pvecthermo[index_th_dmu_idm_dr_] = pth->a_idm_dr*pow((1. + z)/1.e7, pth->nindex_idm_dr)*pba->Omega0_idm_dr*pow(pba->h, 2);
      pvecthermo[index_th_ddmu_idm_dr_] =  -pvecback[background_module_->index_bg_H_]*pth->nindex_idm_dr/(1 + z)*pvecthermo[index_th_dmu_idm_dr_];
      pvecthermo[index_th_dddmu_idm_dr_] = (pvecback[background_module_->index_bg_H_]*pvecback[background_module_->index_bg_H_]/(1. + z) -
                                            pvecback[background_module_->index_bg_H_prime_])*pth->nindex_idm_dr/(1. + z)*pvecthermo[index_th_dmu_idm_dr_];

      /* calculate dmu_idr (self interaction) */
      pvecthermo[index_th_dmu_idr_] = pth->b_idr*pow((1. + z)/1.e7, pth->nindex_idm_dr)*pba->Omega0_idr*pow(pba->h, 2);

      /* extrapolate optical depth of idm_dr and idr */
      pvecthermo[index_th_tau_idm_dr_] = thermodynamics_table_[(tt_size_ - 1)*th_size_+index_th_tau_idm_dr_] +
        (thermodynamics_table_[(tt_size_ - 1)*th_size_ + index_th_tau_idm_dr_] - thermodynamics_table_[(tt_size_ - 2)*th_size_ + index_th_tau_idm_dr_])
        *(z - z_table_[tt_size_ - 1])/(z_table_[tt_size_ - 1] - z_table_[tt_size_ - 2]);

      pvecthermo[index_th_tau_idr_] = thermodynamics_table_[(tt_size_ - 1)*th_size_ + index_th_tau_idr_] +
        (thermodynamics_table_[(tt_size_ - 1)*th_size_ + index_th_tau_idr_] - thermodynamics_table_[(tt_size_ - 2)*th_size_ + index_th_tau_idr_])
        *(z - z_table_[tt_size_ - 1])/(z_table_[tt_size_ - 1] - z_table_[tt_size_ - 2]);

      /* extrapolate idm_dr visibility function */
      pvecthermo[index_th_g_idm_dr_] = thermodynamics_table_[(tt_size_ - 1)*th_size_ + index_th_g_idm_dr_];

      /* calculate interacting dark matter sound speed */
      pvecthermo[index_th_cidm_dr2_] = 4*_k_B_*pba->T_idr*(1. + z)/_eV_/3./pth->m_idm;

      /* calculate interacting dark matter temperature (equal to idr temperature at this redhsift) */
      pvecthermo[index_th_Tidm_dr_] = pba->T_idr*(1. + z);
    }

  }

  /** - interpolate in table with array_interpolate_spline() (normal
      mode) or array_interpolate_spline_growing_closeby() (closeby
      mode) */

  else {

    /* some very specific cases require linear interpolation because of a break in the derivative of the functions */
    if (((pth->reio_parametrization == reio_half_tanh) && (z < 2*z_reionization_))
        || ((pth->reio_parametrization == reio_inter) && (z < 50.))) {

      class_call(array_interpolate_linear(
                                          z_table_,
                                          tt_size_,
                                          thermodynamics_table_,
                                          th_size_,
                                          z,
                                          last_index,
                                          pvecthermo,
                                          th_size_,
                                          error_message_),
                 error_message_,
                 error_message_);
    }

    /* in the "normal" case, use spline interpolation */
    else {

      if (inter_mode == inter_normal_) {

        class_call(array_interpolate_spline(
                                            z_table_,
                                            tt_size_,
                                            thermodynamics_table_,
                                            d2thermodynamics_dz2_table_,
                                            th_size_,
                                            z,
                                            last_index,
                                            pvecthermo,
                                            th_size_,
                                            error_message_),
                   error_message_,
                   error_message_);
      }

      if (inter_mode == inter_closeby_) {

        class_call(array_interpolate_spline_growing_closeby(
                                                            z_table_,
                                                            tt_size_,
                                                            thermodynamics_table_,
                                                            d2thermodynamics_dz2_table_,
                                                            th_size_,
                                                            z,
                                                            last_index,
                                                            pvecthermo,
                                                            th_size_,
                                                            error_message_),
                   error_message_,
                   error_message_);

      }
    }
  }
  return _SUCCESS_;
}

/**
 * Initialize the thermo structure, and in particular the
 * thermodynamics interpolation table.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param pth Input/Output: pointer to initialized thermo structure
 * @return the error status
 */
int ThermodynamicsModule::thermodynamics_init() {

  /** Summary: */

  /** - define local variables */

  /* index running over time*/
  int index_tau;
  /* temporary variables related to visibility function */
  double g;
  /* vector of background values for calling background_at_tau() */
  double * pvecback;
  /* index for calling background_at_tau() */
  int last_index_back;
  /* temporary table of values of tau associated with z values in z_table_ */
  double * tau_table;
  /* same ordered in growing time rather than growing redshift */
  double * tau_table_growing;
  /* conformal time of reionization */
  double tau_reio;
  /* R = (3./4.)*(rho_b/rho_g) */
  double R;
  /* structures for storing temporarily information on recombination
     and reionization */
  struct recombination reco;
  struct reionization reio;
  struct recombination * preco;
  struct reionization * preio;

  double tau,tau_ini;
  double g_max;
  int index_tau_max;
  double dkappa_ini;

  double z_idm_dr, z_idr, tau_idm_dr, tau_idr, Gamma_heat_idm_dr, dTdz_idm_dr, T_idm_dr, z, T_idr, dz, T_adia, z_adia;

  double tau_idm_dr_fs=0.;
  int index_tau_fs;
  int n, N_sub_steps;
  double dz_sub_step;

  tau_reionization_ = pth->tau_reio;
  z_reionization_ = pth->z_reio;
  YHe_ = pth->YHe;

  if (pth->thermodynamics_verbose > 0)
    printf("Computing thermodynamics");

  /** - compute and check primordial Helium fraction  */

  /* Y_He */
  if (YHe_ == _BBN_) {
    class_call(thermodynamics_helium_from_bbn(), error_message_, error_message_);
    if (pth->thermodynamics_verbose > 0)
      printf(" with Y_He=%.4f\n", YHe_);
  }
  else {
    if (pth->thermodynamics_verbose > 0)
      printf("\n");
  }

  class_test((YHe_ < _YHE_SMALL_) || (YHe_ > _YHE_BIG_),
             error_message_,
             "Y_He=%g out of bounds (%g<Y_He<%g)", YHe_, _YHE_SMALL_, _YHE_BIG_);

  /** - check energy injection parameters */

  class_test((pth->annihilation<0),
             error_message_,
             "annihilation parameter cannot be negative");

  class_test((pth->annihilation>1.e-4),
             error_message_,
             "annihilation parameter suspiciously large (%e, while typical bounds are in the range of 1e-7 to 1e-6)",
             pth->annihilation);

  class_test((pth->annihilation_variation>0),
             error_message_,
             "annihilation variation parameter must be negative (decreasing annihilation rate)");

  class_test((pth->annihilation_z<0),
             error_message_,
             "characteristic annihilation redshift cannot be negative");

  class_test((pth->annihilation_zmin<0),
             error_message_,
             "characteristic annihilation redshift cannot be negative");

  class_test((pth->annihilation_zmax<0),
             error_message_,
             "characteristic annihilation redshift cannot be negative");

  class_test((pth->annihilation>0) && ((pba->has_cdm==_FALSE_)&&(pba->has_idm_dr==_FALSE_)),
             error_message_,
             "CDM annihilation effects require the presence of CDM or IDM!");

  class_test((pth->annihilation_f_halo>0) && (pth->recombination==recfast),
             error_message_,
             "Switching on DM annihilation in halos requires using HyRec instead of RECFAST. Otherwise some values go beyond their range of validity in the RECFAST fits, and the thermodynamics module fails. Two solutions: add 'recombination = HyRec' to your input, or set 'annihilation_f_halo = 0.' (default).");

  class_test((pth->annihilation_f_halo<0),
             error_message_,
             "Parameter for DM annihilation in halos cannot be negative");

  class_test((pth->annihilation_z_halo<0),
             error_message_,
             "Parameter for DM annihilation in halos cannot be negative");

  if (pth->thermodynamics_verbose > 0)
    if ((pth->annihilation >0) && (pth->reio_parametrization == reio_none) && (ppr->recfast_Heswitch >= 3) && (pth->recombination==recfast))
      printf("Warning: if you have DM annihilation and you use recfast with option recfast_Heswitch >= 3, then the expression for CfHe_t and dy[1] becomes undefined at late times, producing nan's. This is however masked by reionization if you are not in reio_none mode.");

  class_test((pth->decay<0),
             error_message_,
             "decay parameter cannot be negative");

  class_test((pth->decay>0)&&((pba->has_cdm==_FALSE_)&&(pba->has_idm_dr==_FALSE_)),
             error_message_,
             "CDM decay effects require the presence of CDM or IDM!");

  /* tests in order to prevent segmentation fault in the following */
  class_test(_not4_ == 0.,
             error_message_,
             "stop to avoid division by zero");
  class_test(YHe_ == 1.,
             error_message_,
             "stop to avoid division by zero");

  /** - initialize pointers */

  preco=&reco;
  preio=&reio;

  /** - assign values to all indices in the structures with thermodynamics_indices()*/

  class_call(thermodynamics_indices(preco, preio),
             error_message_,
             error_message_);

  /** - allocate background vector */

  class_alloc(pvecback, background_module_->bg_size_*sizeof(double), error_message_);

  /** - solve recombination and store values of \f$ z, x_e, d \kappa / d \tau, T_b, c_b^2 \f$ with thermodynamics_recombination() */

  class_call_except(thermodynamics_recombination(preco, pvecback),
                    error_message_,
                    error_message_,
                    free(pvecback));

  /** - if there is reionization, solve reionization and store values of \f$ z, x_e, d \kappa / d \tau, T_b, c_b^2 \f$ with thermodynamics_reionization()*/

  if (pth->reio_parametrization != reio_none) {
    class_call_except(thermodynamics_reionization(preco, preio, pvecback),
                      error_message_,
                      error_message_,
                      free(preco->recombination_table);free(pvecback));
  }
  else {
    preio->rt_size=0;
    preio->index_reco_when_reio_start=-1;
  }

  /** - merge tables in recombination and reionization structures into
      a single table in thermo structure */

  class_call(thermodynamics_merge_reco_and_reio(preco, preio),
             error_message_,
             error_message_);

  /** - compute table of corresponding conformal times */

  class_alloc(tau_table, tt_size_*sizeof(double), error_message_);

  for (index_tau = 0; index_tau < tt_size_; index_tau++) {
    class_call(background_module_->background_tau_of_z(z_table_[index_tau], tau_table + index_tau),
               background_module_->error_message_,
               error_message_);
  }

  /** - store initial value of conformal time in the structure */

  tau_ini_ = tau_table[tt_size_ - 1];

  /** - fill missing columns (quantities not computed previously but related) */

  /** - --> minus the baryon drag interaction rate, -dkappa_d/dtau = -[1/R * kappa'], with R = 3 rho_b / 4 rho_gamma, stored temporarily in column ddkappa */

  last_index_back = background_module_->bg_size_ - 1;

  for (index_tau = 0; index_tau < tt_size_; index_tau++) {

    class_call(background_module_->background_at_tau(
                                 tau_table[index_tau],
                                 pba->normal_info,
                                 pba->inter_closeby,
                                 &last_index_back,
                                 pvecback),
               background_module_->error_message_,
               error_message_);

    R = 3./4.*pvecback[background_module_->index_bg_rho_b_]/pvecback[background_module_->index_bg_rho_g_];

    thermodynamics_table_[index_tau*th_size_+index_th_ddkappa_] = -1./R*thermodynamics_table_[index_tau*th_size_ + index_th_dkappa_];

    if(pba->has_idm_dr == _TRUE_) {

      /* - --> idr interaction rate with idm_dr (i.e. idr opacity to idm_dr scattering) */
      thermodynamics_table_[index_tau*th_size_ + index_th_dmu_idm_dr_] =
        pth->a_idm_dr*pow((1. + z_table_[index_tau])/1.e7, pth->nindex_idm_dr)*pba->Omega0_idm_dr*pow(pba->h, 2);

      /* - --> idm_dr interaction rate with idr (i.e. idm_dr opacity
               to idr scattering), [Sinv*dmu_idm_dr] with Sinv = (4
               rho_idr) / (3 rho_idm_dr), stored temporarily in
               ddmu_idm_dr */
      thermodynamics_table_[index_tau*th_size_ + index_th_ddmu_idm_dr_] =
        4./3.*pvecback[background_module_->index_bg_rho_idr_]/pvecback[background_module_->index_bg_rho_idm_dr_]
        *thermodynamics_table_[index_tau*th_size_ + index_th_dmu_idm_dr_];

      /* - --> idr self-interaction rate */
      thermodynamics_table_[index_tau*th_size_ + index_th_dmu_idr_] =
        pth->b_idr*pow((1. + z_table_[index_tau])/1.e7, pth->nindex_idm_dr)*pba->Omega0_idr*pow(pba->h, 2);
    }
  }

  /** - --> second derivative of this rate, -[1/R * kappa']'', stored temporarily in column dddkappa */
  class_call(array_spline_table_line_to_line(tau_table,
                                             tt_size_,
                                             thermodynamics_table_,
                                             th_size_,
                                             index_th_ddkappa_,
                                             index_th_dddkappa_,
                                             _SPLINE_EST_DERIV_,
                                             error_message_),
             error_message_,
             error_message_);

  /** - --> compute tau_d = [int_{tau_today}^{tau} dtau -dkappa_d/dtau] */
  class_call(array_integrate_spline_table_line_to_line(tau_table,
                                                       tt_size_,
                                                       thermodynamics_table_,
                                                       th_size_,
                                                       index_th_ddkappa_,
                                                       index_th_dddkappa_,
                                                       index_th_tau_d_,
                                                       error_message_),
             error_message_,
             error_message_);

  /* the temporary quantities stored in columns ddkappa and dddkappa
     will not be used anymore, so they can be overwritten by other
     intermediate steps of other computations */

  if(pba->has_idm_dr == _TRUE_){

    /** --> second derivative of idm_dr interaction rate (with idr), [Sinv*dmu_idm_dr]'', stored temporarily in column dddmu */
    class_call(array_spline_table_line_to_line(tau_table,
                                               tt_size_,
                                               thermodynamics_table_,
                                               th_size_,
                                               index_th_ddmu_idm_dr_,
                                               index_th_dddmu_idm_dr_,
                                               _SPLINE_EST_DERIV_,
                                               error_message_),
               error_message_,
               error_message_);

    /** - --> compute optical depth of idm, tau_idm_dr = [int_{tau_today}^{tau} dtau [Sinv*dmu_idm_dr] ].
              This step gives -tau_idm_dr. The resulty is mutiplied by -1 later on. */
    class_call(array_integrate_spline_table_line_to_line(tau_table,
                                                         tt_size_,
                                                         thermodynamics_table_,
                                                         th_size_,
                                                         index_th_ddmu_idm_dr_,
                                                         index_th_dddmu_idm_dr_,
                                                         index_th_tau_idm_dr_,
                                                         error_message_),
               error_message_,
               error_message_);


    /** - --> second derivative of idr interaction rate (with idm_dr), [dmu_idm_idr]'', stored temporarily in column dddmu */
    class_call(array_spline_table_line_to_line(tau_table,
                                               tt_size_,
                                               thermodynamics_table_,
                                               th_size_,
                                               index_th_dmu_idm_dr_,
                                               index_th_dddmu_idm_dr_,
                                               _SPLINE_EST_DERIV_,
                                               error_message_),
               error_message_,
               error_message_);

    /** - --> compute optical depth of idr, tau_idr = [int_{tau_today}^{tau} dtau [dmu_idm_idr] ].
              This step gives -tau_idr. The resulty is mutiplied by -1 later on. */
    class_call(array_integrate_spline_table_line_to_line(tau_table,
                                                         tt_size_,
                                                         thermodynamics_table_,
                                                         th_size_,
                                                         index_th_dmu_idm_dr_,
                                                         index_th_dddmu_idm_dr_,
                                                         index_th_tau_idr_,
                                                         error_message_),
               error_message_,
               error_message_);
  }

  /** - --> compute damping scale:

      r_d = 2pi/k_d = 2pi * [int_{tau_ini}^{tau} dtau (1/kappa') 1/6 (R^2+16/15(1+R))/(1+R)^2]^1/2
      = 2pi * [int_{tau_ini}^{tau} dtau (1/kappa') 1/6 (R^2/(1+R)+16/15)/(1+R)]^1/2

      which is like in CosmoTherm (CT), but slightly
      different from Wayne Hu (WH)'s thesis eq. (5.59):
      the factor 16/15 in CT is 4/5 in WH */

  if (pth->compute_damping_scale == _TRUE_) {

    class_alloc(tau_table_growing, tt_size_*sizeof(double), error_message_);

    /* compute integrand and store temporarily in column "ddkappa" */
    for (index_tau = 0; index_tau < tt_size_; index_tau++) {

      tau_table_growing[index_tau] = tau_table[tt_size_ - 1 - index_tau];

      class_call(background_module_->background_at_tau(
                                   tau_table_growing[index_tau],
                                   pba->normal_info,
                                   pba->inter_closeby,
                                   &last_index_back,
                                   pvecback),
                 background_module_->error_message_,
                 error_message_);

      R = 3./4.*pvecback[background_module_->index_bg_rho_b_]/pvecback[background_module_->index_bg_rho_g_];

      thermodynamics_table_[index_tau*th_size_ + index_th_ddkappa_] =
        1./6./thermodynamics_table_[(tt_size_ - 1 - index_tau)*th_size_+index_th_dkappa_]
        *(R*R/(1 + R)+16./15.)/(1. + R);

    }

    /* compute second derivative of integrand, and store temporarily in column "dddkappa" */
    class_call(array_spline_table_line_to_line(tau_table_growing,
                                               tt_size_,
                                               thermodynamics_table_,
                                               th_size_,
                                               index_th_ddkappa_,
                                               index_th_dddkappa_,
                                               _SPLINE_EST_DERIV_,
                                               error_message_),
               error_message_,
               error_message_);


    /* compute integral and store temporarily in column "g" */
    class_call(array_integrate_spline_table_line_to_line(tau_table_growing,
                                                         tt_size_,
                                                         thermodynamics_table_,
                                                         th_size_,
                                                         index_th_ddkappa_,
                                                         index_th_dddkappa_,
                                                         index_th_g_,
                                                         error_message_),
               error_message_,
               error_message_);

    free(tau_table_growing);

    /* we could now write the result as r_d = 2pi * sqrt(integral),
       but we will first better acount for the contribution frokm the tau_ini boundary.
       Close to this boundary, R=0 and the integrand is just 16/(15*6)/kappa'
       Using kappa' propto 1/a^2 and tau propro a during RD, we get the analytic result:
       int_0^{tau_ini} dtau / kappa' = tau_ini / 3 / kappa'_ini
       Thus r_d = 2pi * sqrt( 16/(15*6*3) * (tau_ini/ kappa'_ini) * integral) */

    tau_ini = tau_table[tt_size_ - 1];
    dkappa_ini = thermodynamics_table_[(tt_size_ - 1)*th_size_ + index_th_dkappa_];

    for (index_tau = 0; index_tau < tt_size_; index_tau++) {

      thermodynamics_table_[index_tau*th_size_ + index_th_r_d_] =
        2.*_PI_*sqrt(16./(15.*6.*3.)*tau_ini/dkappa_ini
                     +thermodynamics_table_[(tt_size_ - 1 - index_tau)*th_size_+index_th_g_]);
    }

  } // end of damping scale calculation

  /** - --> second derivative with respect to tau of dkappa (in view of spline interpolation) */
  class_call(array_spline_table_line_to_line(tau_table,
                                             tt_size_,
                                             thermodynamics_table_,
                                             th_size_,
                                             index_th_dkappa_,
                                             index_th_dddkappa_,
                                             _SPLINE_EST_DERIV_,
                                             error_message_),
             error_message_,
             error_message_);

  /** - --> first derivative with respect to tau of dkappa (using spline interpolation) */
  class_call(array_derive_spline_table_line_to_line(tau_table,
                                                    tt_size_,
                                                    thermodynamics_table_,
                                                    th_size_,
                                                    index_th_dkappa_,
                                                    index_th_dddkappa_,
                                                    index_th_ddkappa_,
                                                    error_message_),
             error_message_,
             error_message_);

  /** - --> compute -kappa = [int_{tau_today}^{tau} dtau dkappa/dtau], store temporarily in column "g" */
  class_call(array_integrate_spline_table_line_to_line(tau_table,
                                                       tt_size_,
                                                       thermodynamics_table_,
                                                       th_size_,
                                                       index_th_dkappa_,
                                                       index_th_dddkappa_,
                                                       index_th_g_,
                                                       error_message_),
             error_message_,
             error_message_);

  /** - --> derivatives of baryon sound speed (only computed if some non-minimal tight-coupling schemes is requested) */
  if (pth->compute_cb2_derivatives == _TRUE_) {

    /** - ---> second derivative with respect to tau of cb2 */
    class_call(array_spline_table_line_to_line(tau_table,
                                               tt_size_,
                                               thermodynamics_table_,
                                               th_size_,
                                               index_th_cb2_,
                                               index_th_ddcb2_,
                                               _SPLINE_EST_DERIV_,
                                               error_message_),
               error_message_,
               error_message_);


    /** - ---> first derivative with respect to tau of cb2 (using spline interpolation) */
    class_call(array_derive_spline_table_line_to_line(tau_table,
                                                      tt_size_,
                                                      thermodynamics_table_,
                                                      th_size_,
                                                      index_th_cb2_,
                                                      index_th_ddcb2_,
                                                      index_th_dcb2_,
                                                      error_message_),
               error_message_,
               error_message_);
  }

  /** - --> compute visibility: \f$ g= (d \kappa/d \tau) e^{- \kappa} \f$ */

  /* loop on z (decreasing z, increasing time) */
  for (index_tau = tt_size_ - 1; index_tau >= 0; index_tau--) {

    /** - ---> compute g */
    g = thermodynamics_table_[index_tau*th_size_ + index_th_dkappa_]*
      exp(thermodynamics_table_[index_tau*th_size_ + index_th_g_]);

    /** - ---> compute exp(-kappa) */
    thermodynamics_table_[index_tau*th_size_+index_th_exp_m_kappa_] =
      exp(thermodynamics_table_[index_tau*th_size_ + index_th_g_]);

    /** - ---> compute g' (the plus sign of the second term is correct, see def of -kappa in thermodynamics module!) */
    thermodynamics_table_[index_tau*th_size_ + index_th_dg_] =
      (thermodynamics_table_[index_tau*th_size_ + index_th_ddkappa_] +
       thermodynamics_table_[index_tau*th_size_ + index_th_dkappa_] *
       thermodynamics_table_[index_tau*th_size_ + index_th_dkappa_]) *
      exp(thermodynamics_table_[index_tau*th_size_ + index_th_g_]);

    /** - ---> compute g''  */
    thermodynamics_table_[index_tau*th_size_ + index_th_ddg_] =
      (thermodynamics_table_[index_tau*th_size_ + index_th_dddkappa_] +
       thermodynamics_table_[index_tau*th_size_ + index_th_dkappa_] *
       thermodynamics_table_[index_tau*th_size_ + index_th_ddkappa_]*3. +
       thermodynamics_table_[index_tau*th_size_ + index_th_dkappa_] *
       thermodynamics_table_[index_tau*th_size_ + index_th_dkappa_] *
       thermodynamics_table_[index_tau*th_size_ + index_th_dkappa_]) *
      exp(thermodynamics_table_[index_tau*th_size_ + index_th_g_]);

    /** - ---> store g */
    thermodynamics_table_[index_tau*th_size_ + index_th_g_] = g;

    /** - ---> compute variation rate */
    class_test(thermodynamics_table_[index_tau*th_size_ + index_th_dkappa_] == 0.,
               error_message_,
               "variation rate diverges");

    thermodynamics_table_[index_tau*th_size_ + index_th_rate_] =
      sqrt(pow(thermodynamics_table_[index_tau*th_size_ + index_th_dkappa_], 2)
           + pow(thermodynamics_table_[index_tau*th_size_ + index_th_ddkappa_]/
                thermodynamics_table_[index_tau*th_size_ + index_th_dkappa_], 2)
           + fabs(thermodynamics_table_[index_tau*th_size_ + index_th_dddkappa_]/
                 thermodynamics_table_[index_tau*th_size_ + index_th_dkappa_]));

    /* - ---> restore correct sign for idm_dr and idr optical depth, and calculate idm_dr visibility function */
    if(pba->has_idm_dr == _TRUE_){

      /* restore the correct sign for tau_idm_dr */
      thermodynamics_table_[index_tau*th_size_ + index_th_tau_idm_dr_] *= -1.;

      /* restore the correct sign for tau_idr */
      thermodynamics_table_[index_tau*th_size_ + index_th_tau_idr_] *= -1.;

      /* visibility function for idm_dr : g_idm_dr = [Sinv*dmu_idm_dr] * exp(-tau_idm_dr) */
      thermodynamics_table_[index_tau*th_size_ + index_th_g_idm_dr_] =
        thermodynamics_table_[index_tau*th_size_ + index_th_ddmu_idm_dr_]
        * exp(-thermodynamics_table_[index_tau*th_size_ + index_th_tau_idm_dr_]);
    }
  }

  /** - smooth the rate (details of smoothing unimportant: only the
      order of magnitude of the rate matters) */
  class_call(array_smooth(thermodynamics_table_,
                          th_size_,
                          tt_size_,
                          index_th_rate_,
                          ppr->thermo_rate_smoothing_radius,
                          error_message_),
             error_message_,
             error_message_);

  /* - ---> fill columns for ddmu_idm_dr and dddmu_idm_dr with true values, and compute idm_dr temperature and sound speed */
  if(pba->has_idm_dr == _TRUE_){

    /** - --> second derivative with respect to tau of dmu_idm_dr (in view of spline interpolation) */
    class_call(array_spline_table_line_to_line(tau_table,
                                               tt_size_,
                                               thermodynamics_table_,
                                               th_size_,
                                               index_th_dmu_idm_dr_,
                                               index_th_dddmu_idm_dr_,
                                               _SPLINE_EST_DERIV_,
                                               error_message_),
               error_message_,
               error_message_);

    /** - --> first derivative with respect to tau of dmu_idm_dr (using spline interpolation) */
    class_call(array_derive_spline_table_line_to_line(tau_table,
                                                      tt_size_,
                                                      thermodynamics_table_,
                                                      th_size_,
                                                      index_th_dmu_idm_dr_,
                                                      index_th_dddmu_idm_dr_,
                                                      index_th_ddmu_idm_dr_,
                                                      error_message_),
               error_message_,
               error_message_);

    /** - --> now compute idm_dr temperature and sound speed in various regimes */

    /* (A) - initial value of T_idm_dr at the maximum z (minimum tau) */

    z = z_table_[tt_size_ - 1];

    class_call(background_module_->background_tau_of_z(z, &(tau)),
               background_module_->error_message_,
               error_message_);

    class_call(background_module_->background_at_tau(tau, pba->short_info, pba->inter_normal, &last_index_back, pvecback),
               background_module_->error_message_,
               error_message_);

    Gamma_heat_idm_dr = 2.*pba->Omega0_idr*pow(pba->h,2)*pth->a_idm_dr*pow((1.+z),(pth->nindex_idm_dr+1.))/pow(1.e7,pth->nindex_idm_dr);

    /* (A1) --> if Gamma is not much smaller than H, set T_idm_dr to T_idm_dr = T_idr = xi*T_gamma (tight coupling solution) */
    if (Gamma_heat_idm_dr > 1.e-3*pvecback[background_module_->index_bg_a_]*pvecback[background_module_->index_bg_H_]) {
      T_idm_dr = pba->T_idr*(1.+z);
      dTdz_idm_dr = pba->T_idr;
    }

    /* (A2) --> otherwise, if Gamma << H, set initial T_idm_dr to the
       approximate analytic solution (Gamma/aH)/(1+(Gamma/aH)*T_idr)
       (eq. (A62) in ETHOS I ) */
    else {
      T_idr = pba->T_idr*(1.+z);
      T_idm_dr = Gamma_heat_idm_dr/(pvecback[background_module_->index_bg_a_]*pvecback[background_module_->index_bg_H_])
        /(1. + Gamma_heat_idm_dr/(pvecback[background_module_->index_bg_a_]*pvecback[background_module_->index_bg_H_]))*T_idr;
      dTdz_idm_dr = 2.*T_idm_dr - Gamma_heat_idm_dr/pvecback[background_module_->index_bg_H_]*(T_idr - T_idm_dr);
    }

    thermodynamics_table_[(tt_size_ - 1)*th_size_ + index_th_Tidm_dr_] = T_idm_dr;
    thermodynamics_table_[(tt_size_ - 1)*th_size_ + index_th_cidm_dr2_] = _k_B_*T_idm_dr/_eV_/pth->m_idm*(1. + dTdz_idm_dr/3./T_idm_dr);

    /* T_adia and z_adia will be used later. They are defined as "the
       last T_idm_dr(z) at which the temperature was evaluated
       explicitely, rather than scaled like a^{-2} (decoupled DM
       regime)". Here we just initialize them. They will be updated
       each time that we recompte T_idm_dr explicitely. */
    T_adia = T_idm_dr;
    z_adia = z;

    /* (B) - iterate over growing tau / decreasing z to find other
       values. At each new z we need to compute the following
       quantities: T_idr, T_idm_dr, Gamma_heat_idm_dr, a, H, dT_idm_dr,/dz,
       c_s_idm_dr^2. They all needed to be known from step to step, even
       if the final goal is only to store T_idm_dr, c_s_idm^2 */
    for (index_tau = tt_size_ - 2; index_tau >= 0; index_tau--) {

      /* (B1) --> tight-coupling solution: Gamma >> H implies T_idm_dr=T_idr=xi*T_gamma */
      if (Gamma_heat_idm_dr > 1.e3*pvecback[background_module_->index_bg_a_]*pvecback[background_module_->index_bg_H_]) {
        z = z_table_[index_tau];
        T_idr = pba->T_idr*(1.+z);
        T_idm_dr = T_idr;
        Gamma_heat_idm_dr = 2.*pba->Omega0_idr*pow(pba->h,2)*pth->a_idm_dr*pow((1.+z),(pth->nindex_idm_dr+1.))/pow(1.e7,pth->nindex_idm_dr);
        class_call(background_module_->background_tau_of_z(z, &(tau)),
                   background_module_->error_message_,
                   error_message_);
        class_call(background_module_->background_at_tau(tau, pba->short_info, pba->inter_normal, &last_index_back, pvecback),
                   background_module_->error_message_,
                   error_message_);
        dTdz_idm_dr =pba->T_idr;
      }

      /* (B2) --> intermediate solution: integrate differential equation equation dT_idm_dr/dz = 2 a T_DM - Gamma/H (T_idr - T_idm_dr) */
      else if (Gamma_heat_idm_dr > 1.e-3*pvecback[background_module_->index_bg_a_]*pvecback[background_module_->index_bg_H_]) {

        dz = z_table_[index_tau + 1] - z_table_[index_tau];

        /* (B2a) ----> if dz << H/Gamma the equation is not too stiff and the traditional forward Euler method converges */
        if (dz < pvecback[background_module_->index_bg_H_]/Gamma_heat_idm_dr/10.) {
          z = z_table_[index_tau];
          T_idr = pba->T_idr*(1.+z);
          T_idm_dr -= dTdz_idm_dr*dz;
          Gamma_heat_idm_dr = 2.*pba->Omega0_idr*pow(pba->h,2)*pth->a_idm_dr*pow((1.+z),(pth->nindex_idm_dr+1.))/pow(1.e7,pth->nindex_idm_dr);
          class_call(background_module_->background_tau_of_z(z, &(tau)),
                     background_module_->error_message_,
                     error_message_);
          class_call(background_module_->background_at_tau(tau, pba->short_info, pba->inter_normal, &last_index_back, pvecback),
                     background_module_->error_message_,
                     error_message_);
          dTdz_idm_dr = 2.*pvecback[background_module_->index_bg_a_]*T_idm_dr - Gamma_heat_idm_dr/(pvecback[background_module_->index_bg_H_])*(T_idr - T_idm_dr);
        }

        /* (B2b) ----> otherwise, the equation is too stiff and the
           traditional forward Euler method diverges with this
           stepsize. But we can just decreasee dz to bring it back
           well within the convergence radius H/Gamma of the
           equation. */
        else {
          N_sub_steps = (int)(dz/ (pvecback[background_module_->index_bg_H_]/Gamma_heat_idm_dr/10.)) + 1;
          dz_sub_step = dz/N_sub_steps;

          /* loop over sub-steps */
          for (n=0; n<N_sub_steps; n++) {

            /* evolve quantities over  sub-step wioth forward Euler method */

            z -= dz_sub_step;
            /* final redshift last sub-step overwritten to avoid small rounding error */
            if (n == (N_sub_steps - 1)) z = z_table_[index_tau];

            T_idr = pba->T_idr*(1.+z);
            T_idm_dr -= dTdz_idm_dr*dz_sub_step;
            Gamma_heat_idm_dr = 2.*pba->Omega0_idr*pow(pba->h,2)*pth->a_idm_dr*pow((1.+z),(pth->nindex_idm_dr+1.))/pow(1.e7,pth->nindex_idm_dr);
            class_call(background_module_->background_tau_of_z(z, &(tau)),
                       background_module_->error_message_,
                       error_message_);
            class_call(background_module_->background_at_tau( tau, pba->short_info, pba->inter_normal, &last_index_back, pvecback),
                       background_module_->error_message_,
                       error_message_);
            dTdz_idm_dr = 2.*pvecback[background_module_->index_bg_a_]*T_idm_dr - Gamma_heat_idm_dr/(pvecback[background_module_->index_bg_H_])*(T_idr - T_idm_dr);
          }
        }

        /* update T_adia, z_adia */
        T_adia = T_idm_dr;
        z_adia = z;
      }

      /* (B3) --> decoupled solution: T_idm_dr scales like a^-2 */
      else {
        z = z_table_[index_tau];
        T_idr = pba->T_idr*(1.+z);
        T_idm_dr = T_adia * pow((1.+z)/(1.+z_adia),2);
        Gamma_heat_idm_dr = 2.*pba->Omega0_idr*pow(pba->h,2)*pth->a_idm_dr*pow((1.+z),(pth->nindex_idm_dr+1.))/pow(1.e7,pth->nindex_idm_dr);
        class_call(background_module_->background_tau_of_z(z, &(tau)),
                   background_module_->error_message_,
                   error_message_);
        class_call(background_module_->background_at_tau(tau, pba->short_info, pba->inter_normal, &last_index_back, pvecback),
                   background_module_->error_message_,
                   error_message_);
        dTdz_idm_dr = 2./(1+z)*T_idm_dr;
      }

      thermodynamics_table_[index_tau*th_size_ + index_th_Tidm_dr_] = T_idm_dr;
      thermodynamics_table_[index_tau*th_size_ + index_th_cidm_dr2_] = _k_B_*T_idm_dr/_eV_/pth->m_idm*(1. + dTdz_idm_dr/3./T_idm_dr);
    }
  }

  free(tau_table);

  /** - fill tables of second derivatives with respect to z (in view of spline interpolation) */

  class_call(array_spline_table_lines(z_table_,
                                      tt_size_,
                                      thermodynamics_table_,
                                      th_size_,
                                      d2thermodynamics_dz2_table_,
                                      _SPLINE_EST_DERIV_,
                                      error_message_),
             error_message_,
             error_message_);

  /** - find maximum of g */

  index_tau = tt_size_ - 1;
  while (z_table_[index_tau] > _Z_REC_MAX_) {
    index_tau--;
  }

  class_test(thermodynamics_table_[(index_tau + 1)*th_size_ + index_th_g_] >
             thermodynamics_table_[index_tau*th_size_ + index_th_g_],
             error_message_,
             "found a recombination redshift greater or equal to the maximum value imposed in thermodynamics.h, z_rec_max=%g",_Z_REC_MAX_);

  while (thermodynamics_table_[(index_tau + 1)*th_size_ + index_th_g_] <
         thermodynamics_table_[index_tau*th_size_ + index_th_g_]) {
    index_tau--;
  }

  g_max = thermodynamics_table_[index_tau*th_size_ + index_th_g_];
  index_tau_max = index_tau;

  /* approximation for maximum of g, using cubic interpolation, assuming equally spaced z's */
  z_rec_ = z_table_[index_tau + 1] + 0.5*(z_table_[index_tau + 1] - z_table_[index_tau])
    *(thermodynamics_table_[(index_tau)*th_size_ + index_th_g_] - 1.*thermodynamics_table_[(index_tau + 2)*th_size_ + index_th_g_])
    /(thermodynamics_table_[(index_tau)*th_size_ + index_th_g_] - 2.*thermodynamics_table_[(index_tau + 1)*th_size_ + index_th_g_]
      + thermodynamics_table_[(index_tau + 2)*th_size_ + index_th_g_]);

  class_test(z_rec_ + ppr->smallest_allowed_variation >= _Z_REC_MAX_,
             error_message_,
             "found a recombination redshift greater or equal to the maximum value imposed in thermodynamics.h, z_rec_max=%g",_Z_REC_MAX_);

  class_test(z_rec_ - ppr->smallest_allowed_variation <= _Z_REC_MIN_,
             error_message_,
             "found a recombination redshift smaller or equal to the maximum value imposed in thermodynamics.h, z_rec_min=%g",_Z_REC_MIN_);

  /** - find conformal recombination time using background_tau_of_z() **/

  class_call(background_module_->background_tau_of_z(z_rec_, &(tau_rec_)),
             background_module_->error_message_,
             error_message_);

  class_call(background_module_->background_at_tau(tau_rec_, pba->long_info, pba->inter_normal, &last_index_back, pvecback),
             background_module_->error_message_,
             error_message_);

  rs_rec_ = pvecback[background_module_->index_bg_rs_];
  ds_rec_ = rs_rec_*pba->a_today/(1. + z_rec_);
  da_rec_ = pvecback[background_module_->index_bg_ang_distance_];
  ra_rec_ = da_rec_*(1. + z_rec_)/pba->a_today;
  angular_rescaling_ = ra_rec_/(background_module_->conformal_age_ - tau_rec_);

  /** - find damping scale at recombination (using linear interpolation) */

  if (pth->compute_damping_scale == _TRUE_) {

    rd_rec_ = (z_table_[index_tau + 1] - z_rec_)/(z_table_[index_tau + 1] - z_table_[index_tau])*thermodynamics_table_[(index_tau)*th_size_ + index_th_r_d_]
      + (z_rec_ - z_table_[index_tau])/(z_table_[index_tau + 1] - z_table_[index_tau])*thermodynamics_table_[(index_tau + 1)*th_size_ + index_th_r_d_];

  }

  /** - find time (always after recombination) at which tau_c/tau
      falls below some threshold, defining tau_free_streaming */

  class_call(background_module_->background_tau_of_z(z_table_[index_tau], &tau),
             background_module_->error_message_,
             error_message_);

  while ((1./thermodynamics_table_[(index_tau)*th_size_ + index_th_dkappa_]/tau < ppr->radiation_streaming_trigger_tau_c_over_tau)
         && (index_tau>0)) {

    index_tau--;

    class_call(background_module_->background_tau_of_z(z_table_[index_tau], &tau),
               background_module_->error_message_,
               error_message_);

  }

  tau_free_streaming_ = tau;

  /** - Find interacting dark radiation free-streaming time */
  index_tau_fs = index_tau;

  if(pba->has_idr == _TRUE_) {

    if(pba->has_idm_dr == _TRUE_) {

      if(pth->nindex_idm_dr>=2){
        index_tau=index_tau_fs-1;
        /* comment: using index_tau_max (index_tau_fs) instead of tt_size_ - 1 ensures that the switch is always after recombination (free streaming) */
      }
      else{
        index_tau=0;
      }

      class_call(background_module_->background_tau_of_z(z_table_[index_tau], &tau),
                 background_module_->error_message_,
                 error_message_);

      while ((1./thermodynamics_table_[(index_tau)*th_size_ + index_th_dmu_idm_dr_]/tau
              < ppr->idr_streaming_trigger_tau_c_over_tau) &&
             ((pth->nindex_idm_dr >= 2 && index_tau > 0) ||
              (pth->nindex_idm_dr < 2 && index_tau < tt_size_ - 1))) {

        if(pth->nindex_idm_dr>=2){
          index_tau--;
        }
        else{
          index_tau++;
        }

        class_call(background_module_->background_tau_of_z(z_table_[index_tau], &tau),
                   background_module_->error_message_,
                   error_message_);

      }

      tau_idm_dr_fs = tau;
      tau_idr_free_streaming_ = tau;
    }

    /* case of idr alone without idm_dr */
    else {
      index_tau=index_tau_fs-1;
      class_call(background_module_->background_tau_of_z(z_table_[index_tau], &tau),
                 background_module_->error_message_,
                 error_message_)
        tau_idm_dr_fs = tau;
      tau_idr_free_streaming_ = tau;
    }
  }

  /** - find z_star (when optical depth kappa crosses one, using linear
      interpolation) and sound horizon at that time */

  index_tau=0;
  while ((thermodynamics_table_[(index_tau)*th_size_ + index_th_exp_m_kappa_] > 1./_E_) && (index_tau < tt_size_))
    index_tau++;

  z_star_ = z_table_[index_tau - 1]+
    (1./_E_-thermodynamics_table_[(index_tau - 1)*th_size_ + index_th_exp_m_kappa_])
    /(thermodynamics_table_[(index_tau)*th_size_ + index_th_exp_m_kappa_] - thermodynamics_table_[(index_tau - 1)*th_size_ + index_th_exp_m_kappa_])
    *(z_table_[index_tau] - z_table_[index_tau - 1]);

  class_call(background_module_->background_tau_of_z(z_star_, &(tau_star_)),
             background_module_->error_message_,
             error_message_);

  class_call(background_module_->background_at_tau(tau_star_, pba->long_info, pba->inter_normal, &last_index_back, pvecback),
             background_module_->error_message_,
             error_message_);

  rs_star_ = pvecback[background_module_->index_bg_rs_];
  ds_star_ = rs_star_*pba->a_today/(1. + z_star_);
  da_star_ = pvecback[background_module_->index_bg_ang_distance_];
  ra_star_ = da_star_*(1. + z_star_)/pba->a_today;

  if (pth->compute_damping_scale == _TRUE_) {

    rd_star_ = (z_table_[index_tau + 1] - z_star_)/(z_table_[index_tau + 1] - z_table_[index_tau])*thermodynamics_table_[(index_tau)*th_size_ + index_th_r_d_]
      + (z_star_ - z_table_[index_tau])/(z_table_[index_tau + 1] - z_table_[index_tau])*thermodynamics_table_[(index_tau + 1)*th_size_ + index_th_r_d_];

  }

  /** - find baryon drag time (when tau_d crosses one, using linear
      interpolation) and sound horizon at that time */

  index_tau=0;
  while ((thermodynamics_table_[(index_tau)*th_size_ + index_th_tau_d_] < 1.) && (index_tau < tt_size_))
    index_tau++;

  z_d_ = z_table_[index_tau - 1] +
    (1. - thermodynamics_table_[(index_tau - 1)*th_size_ + index_th_tau_d_])
    /(thermodynamics_table_[(index_tau)*th_size_ + index_th_tau_d_] - thermodynamics_table_[(index_tau - 1)*th_size_ + index_th_tau_d_])
    *(z_table_[index_tau] - z_table_[index_tau - 1]);

  class_call(background_module_->background_tau_of_z(z_d_, &(tau_d_)),
             background_module_->error_message_,
             error_message_);

  class_call(background_module_->background_at_tau(tau_d_, pba->long_info, pba->inter_normal, &last_index_back, pvecback),
             background_module_->error_message_,
             error_message_);

  rs_d_ = pvecback[background_module_->index_bg_rs_];
  ds_d_ = rs_d_*pba->a_today/(1. + z_d_);

  /** - find idm_dr and idr drag times */
  if(pba->has_idm_dr == _TRUE_){

    if((thermodynamics_table_[(tt_size_ - 1)*th_size_ + index_th_tau_idm_dr_] > 1.) &&
       (thermodynamics_table_[(tt_size_ - 1)*th_size_ + index_th_tau_idr_] > 1.)) {
      index_tau=0;

      while ((thermodynamics_table_[(index_tau)*th_size_ + index_th_tau_idm_dr_] < 1.) && (index_tau < tt_size_ - 1))
        index_tau++;

      z_idm_dr = z_table_[index_tau - 1] + (1. - thermodynamics_table_[(index_tau - 1)*th_size_ + index_th_tau_idm_dr_])
        /(thermodynamics_table_[(index_tau)*th_size_ + index_th_tau_idm_dr_] - thermodynamics_table_[(index_tau - 1)*th_size_ + index_th_tau_idm_dr_])
        *(z_table_[index_tau] - z_table_[index_tau - 1]);

      class_call(background_module_->background_tau_of_z(z_idm_dr, &(tau_idm_dr)),
                 background_module_->error_message_,
                 error_message_);

      index_tau = 0;

      while ((thermodynamics_table_[(index_tau)*th_size_ + index_th_tau_idr_] < 1.) && (index_tau < tt_size_ - 1))
        index_tau++;

      z_idr = z_table_[index_tau - 1] +
        (1. - thermodynamics_table_[(index_tau - 1)*th_size_ + index_th_tau_idr_])
        /(thermodynamics_table_[(index_tau)*th_size_ + index_th_tau_idr_] - thermodynamics_table_[(index_tau - 1)*th_size_ + index_th_tau_idr_])
        *(z_table_[index_tau] - z_table_[index_tau - 1]);

      class_call(background_module_->background_tau_of_z(z_idr, &(tau_idr)),
                 background_module_->error_message_,
                 error_message_);
    }
  }

  /** - find time above which visibility falls below a given fraction of its maximum */

  index_tau=index_tau_max;
  while ((thermodynamics_table_[(index_tau)*th_size_ + index_th_g_] > g_max*ppr->neglect_CMB_sources_below_visibility)
         && (index_tau > 0))
    index_tau--;

  class_call(background_module_->background_tau_of_z(z_table_[index_tau], &(tau_cut_)),
             background_module_->error_message_,
             error_message_);

  /** - if verbose flag set to next-to-minimum value, print the main results */

  if (pth->thermodynamics_verbose > 0) {

    if(pba->has_idm_dr == _TRUE_) {
      if((thermodynamics_table_[(tt_size_ - 1)*th_size_ + index_th_tau_idm_dr_] > 1.) &&
         (thermodynamics_table_[(tt_size_ - 1)*th_size_ + index_th_tau_idr_] > 1.)) {
        printf(" -> idr decouples at tau_idr = %e Mpc\n",tau_idr);
        printf(" -> idm_dr decouples at tau_idm_dr = %e Mpc\n",tau_idm_dr);
      }
      else{
        printf(" -> computation of decoupling time of idm_dr and idr skipped, because z would not be in z_table\n");
      }
    }
    printf(" -> recombination at z = %f (max of visibility function)\n", z_rec_);
    printf("    corresponding to conformal time = %f Mpc\n", tau_rec_);
    printf("    with comoving sound horizon = %f Mpc\n", rs_rec_);
    printf("    angular diameter distance = %f Mpc\n", da_rec_);
    printf("    and sound horizon angle 100*theta_s = %f\n", 100.*rs_rec_/ra_rec_);
    if (pth->compute_damping_scale == _TRUE_) {
      printf("    and with comoving photon damping scale = %f Mpc\n", rd_rec_);
      printf("    or comoving damping wavenumber k_d = %f 1/Mpc\n", 2.*_PI_/rd_rec_);
    }
    printf("    Thomson optical depth crosses one at z_* = %f\n", z_star_);
    printf("    giving an angle 100*theta_* = %f\n", 100.*rs_star_/ra_star_);
    printf(" -> baryon drag stops at z = %f\n", z_d_);
    printf("    corresponding to conformal time = %f Mpc\n" ,tau_d_);
    printf("    with comoving sound horizon rs = %f Mpc\n" ,rs_d_);
    if ((pth->reio_parametrization == reio_camb) || (pth->reio_parametrization == reio_half_tanh)) {
      if (pth->reio_z_or_tau==reio_tau)
        printf(" -> reionization  at z = %f\n", z_reionization_);
      if (pth->reio_z_or_tau==reio_z)
        printf(" -> reionization with optical depth = %f\n", tau_reionization_);
      class_call(background_module_->background_tau_of_z(z_reionization_, &tau_reio),
                 background_module_->error_message_,
                 error_message_);
      printf("    corresponding to conformal time = %f Mpc\n",tau_reio);
    }
    if (pth->reio_parametrization == reio_bins_tanh) {
      printf(" -> binned reionization gives optical depth = %f\n", tau_reionization_);
    }
    if (pth->reio_parametrization == reio_many_tanh) {
      printf(" -> many-step reionization gives optical depth = %f\n", tau_reionization_);
    }
    if (pth->reio_parametrization == reio_inter) {
      printf(" -> interpolated reionization history gives optical depth = %f\n", tau_reionization_);
    }
    if (pth->thermodynamics_verbose > 1) {
      printf(" -> free-streaming approximation can be turned on as soon as tau=%g Mpc\n", tau_free_streaming_);
    }
    if ((pba->has_idr)&&(pth->thermodynamics_verbose > 1)) {
      printf(" -> dark free-streaming approximation can be turned on as soon as tau=%g Mpc\n",
             tau_idm_dr_fs);
    }
  }

  free(pvecback);

  return _SUCCESS_;
}

/**
 * Free all memory space allocated by thermodynamics_init().
 *
 *
 * @param pth Input/Output: pointer to thermo structure (to be freed)
 * @return the error status
 */

int ThermodynamicsModule::thermodynamics_free() {

  free(z_table_);
  free(thermodynamics_table_);
  free(d2thermodynamics_dz2_table_);

  return _SUCCESS_;
}

/**
 * Assign value to each relevant index in vectors of thermodynamical quantities,
 * as well as in vector containing reionization parameters.
 *
 *
 * @param pba   Input: pointer to background structure
 * @param pth   Input/Output: pointer to thermo structure
 * @param preco Input/Output: pointer to recombination structure
 * @param preio Input/Output: pointer to reionization structure
 * @return the error status
 */

int ThermodynamicsModule::thermodynamics_indices(recombination* preco, reionization* preio) {

  /** Summary: */

  /** - define local variables */

  /* a running index for the vector of thermodynamics quantities */
  int index;

  /** - initialization of all indices and flags in thermo structure */
  index = 0;

  index_th_xe_ = index;
  index++;
  index_th_dkappa_ = index;
  index++;
  index_th_tau_d_ = index;
  index++;
  index_th_ddkappa_ = index;
  index++;
  index_th_dddkappa_ = index;
  index++;
  index_th_exp_m_kappa_ = index;
  index++;
  index_th_g_ = index;
  index++;
  index_th_dg_ = index;
  index++;
  index_th_ddg_ = index;
  index++;
  index_th_Tb_ = index;
  index++;
  index_th_wb_ = index;
  index++;
  index_th_cb2_ = index;
  index++;

  if(pba->has_idm_dr == _TRUE_){
    index_th_dmu_idm_dr_ = index;
    index++;
    index_th_ddmu_idm_dr_ = index;
    index++;
    index_th_dddmu_idm_dr_ = index;
    index++;
    index_th_tau_idm_dr_ = index;
    index++;
    index_th_tau_idr_ = index;
    index++;
    index_th_g_idm_dr_ = index;
    index++;
    index_th_cidm_dr2_ = index;
    index++;
    index_th_Tidm_dr_ = index;
    index++;
    index_th_dmu_idr_ = index;
    index++;
  }

  /* derivatives of baryon sound speed (only computed if some non-minimal tight-coupling schemes is requested) */
  if (pth->compute_cb2_derivatives == _TRUE_) {
    index_th_dcb2_ = index;
    index++;
    index_th_ddcb2_ = index;
    index++;
  }

  index_th_rate_ = index;
  index++;

  if (pth->compute_damping_scale == _TRUE_) {
    index_th_r_d_ = index;
    index++;
  }

  /* end of indices */
  th_size_ = index;

  /** - initialization of all indices and flags in recombination structure */
  index = 0;

  preco->index_re_z = index;
  index++;
  preco->index_re_xe = index;
  index++;
  preco->index_re_dkappadtau = index;
  index++;
  preco->index_re_Tb = index;
  index++;
  preco->index_re_wb = index;
  index++;
  preco->index_re_cb2 = index;
  index++;

  /* end of indices */
  preco->re_size = index;

  /** - initialization of all indices and flags in reionization structure */
  index = 0;

  preio->index_re_z = index;
  index++;
  preio->index_re_xe = index;
  index++;
  preio->index_re_Tb = index;
  index++;
  preio->index_re_wb = index;
  index++;
  preio->index_re_cb2 = index;
  index++;
  preio->index_re_dkappadtau = index;
  index++;
  preio->index_re_dkappadz = index;
  index++;
  preio->index_re_d3kappadz3 = index;
  index++;

  /* end of indices */
  preio->re_size = index;

  /** - same with parameters of the function \f$ X_e(z)\f$ */

  index=0;

  preio->index_reio_start = index;
  index++;

  /* case where x_e(z) taken like in CAMB (other cases can be added) */
  if ((pth->reio_parametrization == reio_camb) || (pth->reio_parametrization == reio_half_tanh)) {

    preio->index_reio_redshift = index;
    index++;
    preio->index_reio_exponent = index;
    index++;
    preio->index_reio_width = index;
    index++;
    preio->index_reio_xe_before = index;
    index++;
    preio->index_reio_xe_after = index;
    index++;
    preio->index_helium_fullreio_fraction = index;
    index++;
    preio->index_helium_fullreio_redshift = index;
    index++;
    preio->index_helium_fullreio_width = index;
    index++;

  }

  /* case where x_e(z) is binned */
  if (pth->reio_parametrization == reio_bins_tanh) {

    /* the code will not only copy here the "bin centers" passed in
       input. It will add an initial and final value for (z,xe). So
       this array has a dimension bigger than the bin center array */

    preio->reio_num_z=pth->binned_reio_num+2; /* add two values: beginning and end of reio */

    preio->index_reio_first_z = index;
    index+= preio->reio_num_z;
    preio->index_reio_first_xe = index;
    index+= preio->reio_num_z;
    preio->index_reio_step_sharpness = index;
    index++;

  }

  /* case where x_e(z) has many tanh jumps */
  if (pth->reio_parametrization == reio_many_tanh) {

    /* the code will not only copy here the "jump centers" passed in
       input. It will add an initial and final value for (z,xe). So
       this array has a dimension bigger than the jump center array */

    preio->reio_num_z=pth->many_tanh_num+2; /* add two values: beginning and end of reio */

    preio->index_reio_first_z = index;
    index+= preio->reio_num_z;
    preio->index_reio_first_xe = index;
    index+= preio->reio_num_z;
    preio->index_reio_step_sharpness = index;
    index++;

  }

  /* case where x_e(z) must be interpolated */
  if (pth->reio_parametrization == reio_inter) {

    preio->reio_num_z=pth->reio_inter_num;

    preio->index_reio_first_z = index;
    index+= preio->reio_num_z;
    preio->index_reio_first_xe = index;
    index+= preio->reio_num_z;

  }

  preio->reio_num_params = index;

  /* flags for calling the interpolation routine */

  inter_normal_ = 0;
  inter_closeby_ = 1;

  return _SUCCESS_;
}

/**
 * Infer the primordial helium fraction from standard BBN, as a
 * function of the baryon density and expansion rate during BBN.
 *
 * This module is simpler then the one used in arXiv:0712.2826 because
 * it neglects the impact of a possible significant chemical
 * potentials for electron neutrinos. The full code with xi_nu_e could
 * be introduced here later.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param pth Input/Output: pointer to initialized thermo structure
 * @return the error status
 */
int ThermodynamicsModule::thermodynamics_helium_from_bbn() {

  FILE * fA;
  char line[_LINE_LENGTH_MAX_];
  char * left;

  int num_omegab=0;
  int num_deltaN=0;

  double * omegab=NULL;
  double * deltaN=NULL;
  double * YHe=NULL;
  double * ddYHe=NULL;
  double * YHe_at_deltaN=NULL;
  double * ddYHe_at_deltaN=NULL;

  int array_line=0;
  double DeltaNeff;
  double omega_b;
  int last_index;
  double Neff_bbn, z_bbn, tau_bbn, *pvecback;

  /**Summary: */
  /** - Infer effective number of neutrinos at the time of BBN */
  class_alloc(pvecback, background_module_->bg_size_*sizeof(double), error_message_);

  /** - 8.6173e-11 converts from Kelvin to MeV. We randomly choose 0.1 MeV to be the temperature of BBN */
  z_bbn = 0.1/(8.6173e-11*pba->T_cmb)-1.0;

  class_call(background_module_->background_tau_of_z(z_bbn, &tau_bbn),
             background_module_->error_message_,
             error_message_);

  class_call(background_module_->background_at_tau(
                                                  tau_bbn,
                                                  pba->long_info,
                                                  pba->inter_normal,
                                                  &last_index,
                                                  pvecback),
             background_module_->error_message_,
             error_message_);

  Neff_bbn = (pvecback[background_module_->index_bg_Omega_r_]*pvecback[background_module_->index_bg_rho_crit_]
              - pvecback[background_module_->index_bg_rho_g_])
    /(7./8.*pow(4./11.,4./3.)*pvecback[background_module_->index_bg_rho_g_]);

  free(pvecback);

  //  printf("Neff early = %g, Neff at bbn: %g\n", background_module_->Neff_, Neff_bbn);

  /** - compute Delta N_eff as defined in bbn file, i.e. \f$ \Delta N_{eff}=0\f$ means \f$ N_{eff}=3.046\f$ */
  DeltaNeff = Neff_bbn - 3.046;

  /* the following file is assumed to contain (apart from comments and blank lines):
     - the two numbers (num_omegab, num_deltaN) = number of values of BBN free parameters
     - three columns (omegab, deltaN, YHe) where omegab = Omega0_b h^2 and deltaN = Neff-3.046 by definition
     - omegab and deltaN are assumed to be arranged as:
     omegab1 deltaN1 YHe
     omegab2 deltaN1 YHe
     .....
     omegab1 delatN2 YHe
     omegab2 deltaN2 YHe
     .....
  */

  class_open(fA, ppr->sBBN_file, "r", error_message_);

  /* go through each line */
  while (fgets(line,_LINE_LENGTH_MAX_-1,fA) != NULL) {

    /* eliminate blank spaces at beginning of line */
    left=line;
    while (left[0]==' ') {
      left++;
    }

    /* check that the line is neither blank neither a comment. In
       ASCII, left[0]>39 means that first non-blank character might
       be the beginning of some data (it is not a newline, a #, a %,
       etc.) */
    if (left[0] > 39) {

      /* if the line contains data, we must interpret it. If
         (num_omegab, num_deltaN)=(0,0), the current line must contain
         their values. Otherwise, it must contain (omegab, delatN,
         YHe). */
      if ((num_omegab==0) && (num_deltaN==0)) {

        /* read (num_omegab, num_deltaN), infer size of arrays and allocate them */
        class_test(sscanf(line,"%d %d",&num_omegab,&num_deltaN) != 2,
                   error_message_,
                   "could not read value of parameters (num_omegab,num_deltaN) in file %s\n",ppr->sBBN_file);

        class_alloc(omegab, num_omegab*sizeof(double), error_message_);
        class_alloc(deltaN, num_deltaN*sizeof(double), error_message_);
        class_alloc(YHe, num_omegab*num_deltaN*sizeof(double), error_message_);
        class_alloc(ddYHe, num_omegab*num_deltaN*sizeof(double), error_message_);
        class_alloc(YHe_at_deltaN, num_omegab*sizeof(double), error_message_);
        class_alloc(ddYHe_at_deltaN, num_omegab*sizeof(double), error_message_);
        array_line=0;

      }
      else {

        /* read (omegab, deltaN, YHe) */
        class_test(sscanf(line,"%lg %lg %lg",
                          &(omegab[array_line%num_omegab]),
                          &(deltaN[array_line/num_omegab]),
                          &(YHe[array_line])
                          ) != 3,
                   error_message_,
                   "could not read value of parameters (omegab,deltaN,YHe) in file %s\n",ppr->sBBN_file);
        array_line ++;
      }
    }
  }

  fclose(fA);

  /** - spline in one dimension (along deltaN) */
  class_call(array_spline_table_lines(deltaN,
                                      num_deltaN,
                                      YHe,
                                      num_omegab,
                                      ddYHe,
                                      _SPLINE_NATURAL_,
                                      error_message_),
             error_message_,
             error_message_);

  omega_b=pba->Omega0_b*pba->h*pba->h;

  class_test_except(omega_b < omegab[0],
                    error_message_,
                    free(omegab);free(deltaN);free(YHe);free(ddYHe);free(YHe_at_deltaN);free(ddYHe_at_deltaN),
                    "You have asked for an unrealistic small value omega_b = %e. The corresponding value of the primordial helium fraction cannot be found in the interpolation table. If you really want this value, you should fix YHe to a given value rather than to BBN",
                    omega_b);

  class_test_except(omega_b > omegab[num_omegab-1],
                    error_message_,
                    free(omegab);free(deltaN);free(YHe);free(ddYHe);free(YHe_at_deltaN);free(ddYHe_at_deltaN),
                    "You have asked for an unrealistic high value omega_b = %e. The corresponding value of the primordial helium fraction cannot be found in the interpolation table. If you really want this value, you should fix YHe to a given value rather than to BBN",
                    omega_b);

  class_test_except(DeltaNeff < deltaN[0],
                    error_message_,
                    free(omegab);free(deltaN);free(YHe);free(ddYHe);free(YHe_at_deltaN);free(ddYHe_at_deltaN),
                    "You have asked for an unrealistic small value of Delta N_eff = %e. The corresponding value of the primordial helium fraction cannot be found in the interpolation table. If you really want this value, you should fix YHe to a given value rather than to BBN",
                    DeltaNeff);

  class_test_except(DeltaNeff > deltaN[num_deltaN-1],
                    error_message_,
                    free(omegab);free(deltaN);free(YHe);free(ddYHe);free(YHe_at_deltaN);free(ddYHe_at_deltaN),
                    "You have asked for an unrealistic high value of Delta N_eff = %e. The corresponding value of the primordial helium fraction cannot be found in the interpolation table. If you really want this value, you should fix YHe to a given value rather than to BBN",
                    DeltaNeff);

  /** - interpolate in one dimension (along deltaN) */
  class_call(array_interpolate_spline(deltaN,
                                      num_deltaN,
                                      YHe,
                                      ddYHe,
                                      num_omegab,
                                      DeltaNeff,
                                      &last_index,
                                      YHe_at_deltaN,
                                      num_omegab,
                                      error_message_),
             error_message_,
             error_message_);

  /** - spline in remaining dimension (along omegab) */
  class_call(array_spline_table_lines(omegab,
                                      num_omegab,
                                      YHe_at_deltaN,
                                      1,
                                      ddYHe_at_deltaN,
                                      _SPLINE_NATURAL_,
                                      error_message_),
             error_message_,
             error_message_);

  /** - interpolate in remaining dimension (along omegab) */
  class_call(array_interpolate_spline(omegab,
                                      num_omegab,
                                      YHe_at_deltaN,
                                      ddYHe_at_deltaN,
                                      1,
                                      omega_b,
                                      &last_index,
                                      &(YHe_),
                                      1,
                                      error_message_),
             error_message_,
             error_message_);

  /** - deallocate arrays */
  free(omegab);
  free(deltaN);
  free(YHe);
  free(ddYHe);
  free(YHe_at_deltaN);
  free(ddYHe_at_deltaN);

  return _SUCCESS_;

}

/**
 * In case of non-minimal cosmology, this function determines the
 * energy rate injected in the IGM at a given redshift z (= on-the-spot
 * annihilation). This energy injection may come e.g. from dark matter
 * annihilation or decay.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param preco Input: pointer to recombination structure
 * @param z Input: redshift
 * @param energy_rate Output: energy density injection rate
 * @param error_message Output: error message
 * @return the error status
 */

int ThermodynamicsModule::thermodynamics_onthespot_energy_injection(recombination* preco, double z, double* energy_rate, ErrorMsg error_message) {

  double annihilation_at_z;
  double rho_cdm_today;
  double u_min;
  double erfc;

  /*redshift-dependent annihilation parameter*/

  if (z>preco->annihilation_zmax) {

    annihilation_at_z = preco->annihilation*
      exp(-preco->annihilation_variation*pow(log((preco->annihilation_z+1.)/(preco->annihilation_zmax+1.)),2));
  }
  else if (z>preco->annihilation_zmin) {

    annihilation_at_z = preco->annihilation*
      exp(preco->annihilation_variation*(-pow(log((preco->annihilation_z+1.)/(preco->annihilation_zmax+1.)),2)
                                         +pow(log((z+1.)/(preco->annihilation_zmax+1.)),2)));
  }
  else {

    annihilation_at_z = preco->annihilation*
      exp(preco->annihilation_variation*(-pow(log((preco->annihilation_z+1.)/(preco->annihilation_zmax+1.)),2)
                                         +pow(log((preco->annihilation_zmin+1.)/(preco->annihilation_zmax+1.)),2)));
  }

  rho_cdm_today = pow(pba->H0*_c_/_Mpc_over_m_,2)*3/8./_PI_/_G_*(pba->Omega0_idm_dr+pba->Omega0_cdm)*_c_*_c_; /* energy density in J/m^3 */

  u_min = (1+z)/(1+preco->annihilation_z_halo);

  erfc = pow(1.+0.278393*u_min+0.230389*u_min*u_min+0.000972*u_min*u_min*u_min+0.078108*u_min*u_min*u_min*u_min,-4);

  *energy_rate = pow(rho_cdm_today,2)/_c_/_c_*pow((1+z),3)*
    (pow((1.+z),3)*annihilation_at_z+preco->annihilation_f_halo*erfc)
    +rho_cdm_today*pow((1+z),3)*preco->decay;
  /* energy density rate in J/m^3/s (remember that annihilation_at_z is in m^3/s/Kg and decay in s^-1) */

  return _SUCCESS_;

}

/**
 * In case of non-minimal cosmology, this function determines the
 * effective energy rate absorbed by the IGM at a given redshift
 * (beyond the on-the-spot annihilation). This energy injection may
 * come e.g. from dark matter annihilation or decay.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param preco Input: pointer to recombination structure
 * @param z Input: redshift
 * @param energy_rate Output: energy density injection rate
 * @param error_message Output: error message
 * @return the error status
 */

int ThermodynamicsModule::thermodynamics_energy_injection(recombination* preco, double z, double* energy_rate, ErrorMsg error_message) {

  double zp,dz;
  double integrand,first_integrand;
  double factor,result;
  double nH0;
  double onthespot;

  if (preco->annihilation > 0) {

    if (preco->has_on_the_spot == _FALSE_) {

      /* number of hydrogen nuclei today in m**-3 */
      nH0 = 3.*preco->H0*preco->H0*pba->Omega0_b/(8.*_PI_*_G_*_m_H_)*(1.-preco->YHe);

      /* factor = c sigma_T n_H(0) / (H(0) \sqrt(Omega_m)) (dimensionless) */
      factor = _sigma_ * nH0 / pba->H0 * _Mpc_over_m_ / sqrt(pba->Omega0_b+pba->Omega0_cdm+pba->Omega0_idm_dr);

      /* integral over z'(=zp) with step dz */
      dz=1.;

      /* first point in trapezoidal integral */
      zp = z;
      class_call(thermodynamics_onthespot_energy_injection(preco, zp, &onthespot, error_message),
                 error_message,
                 error_message);
      first_integrand = factor*pow(1+z,8)/pow(1+zp,7.5)*exp(2./3.*factor*(pow(1+z,1.5)-pow(1+zp,1.5)))*onthespot; // beware: versions before 2.4.3, there were wrong exponents: 6 and 5.5 instead of 8 and 7.5
      result = 0.5*dz*first_integrand;

      /* other points in trapezoidal integral */
      do {

        zp += dz;
        class_call(thermodynamics_onthespot_energy_injection(preco, zp, &onthespot, error_message),
                   error_message,
                   error_message);
        integrand = factor*pow(1+z,8)/pow(1+zp,7.5)*exp(2./3.*factor*(pow(1+z,1.5)-pow(1+zp,1.5)))*onthespot; // beware: versions before 2.4.3, there were wrong exponents: 6 and 5.5 instead of 8 and 7.5
        result += dz*integrand;

      } while (integrand/first_integrand > 0.02);

      /* uncomment these lines if you also want to compute the on-the-spot for comparison */
      class_call(thermodynamics_onthespot_energy_injection(preco, z, &onthespot, error_message),
                 error_message,
                 error_message);

    }
    else {
      class_call(thermodynamics_onthespot_energy_injection(preco, z, &result, error_message),
                 error_message,
                 error_message);
    }

    /* these test lines print the energy rate rescaled by (1+z)^6 in J/m^3/s, with or without the on-the-spot approximation */
    /*
      fprintf(stdout,"%e  %e  %e \n",
      1.+z,
      result/pow(1.+z,6),
      onthespot/pow(1.+z,6));
    */

    /* effective energy density rate in J/m^3/s  */
    *energy_rate = result;

  }
  else {
    *energy_rate = 0.;
  }

  return _SUCCESS_;

}

/**
 * This subroutine contains the reionization function \f$ X_e(z) \f$
 * (one for each scheme; so far, only the function corresponding to
 * the reio_camb scheme is coded)
 *
 * @param z     Input: redshift
 * @param pth   Input: pointer to thermo structure, to know which scheme is used
 * @param preio Input: pointer to reionization structure, containing the parameters of the function \f$ X_e(z) \f$
 * @param xe    Output: \f$ X_e(z) \f$
 */

int ThermodynamicsModule::thermodynamics_reionization_function(double z, reionization* preio, double* xe) {

  /** Summary: */

  /** - define local variables */
  double argument;
  int i;
  double z_jump;

  int jump;
  double center,before, after,width,one_jump;

  /** - implementation of ionization function similar to the one in CAMB */

  if ((pth->reio_parametrization == reio_camb) || (pth->reio_parametrization == reio_half_tanh)) {

    /** - --> case z > z_reio_start */

    if (z > preio->reionization_parameters[preio->index_reio_start]) {

      *xe = preio->reionization_parameters[preio->index_reio_xe_before];

    }

    else {

      /** - --> case z < z_reio_start: hydrogen contribution (tanh of complicated argument) */

      argument = (pow((1.+preio->reionization_parameters[preio->index_reio_redshift]),
                      preio->reionization_parameters[preio->index_reio_exponent])
                  - pow((1.+z),preio->reionization_parameters[preio->index_reio_exponent]))
        /(preio->reionization_parameters[preio->index_reio_exponent]
          /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
          *pow((1.+preio->reionization_parameters[preio->index_reio_redshift]),
               (preio->reionization_parameters[preio->index_reio_exponent]-1.)))
        /preio->reionization_parameters[preio->index_reio_width];
      /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */

      if (pth->reio_parametrization == reio_camb) {
        *xe = (preio->reionization_parameters[preio->index_reio_xe_after]
               -preio->reionization_parameters[preio->index_reio_xe_before])
          *(tanh(argument)+1.)/2.
          +preio->reionization_parameters[preio->index_reio_xe_before];
      }
      else {
        *xe = (preio->reionization_parameters[preio->index_reio_xe_after]
               -preio->reionization_parameters[preio->index_reio_xe_before])
          *tanh(argument)
          +preio->reionization_parameters[preio->index_reio_xe_before];
      }

      /** - --> case z < z_reio_start: helium contribution (tanh of simpler argument) */

      if (pth->reio_parametrization == reio_camb) {
        argument = (preio->reionization_parameters[preio->index_helium_fullreio_redshift] - z)
          /preio->reionization_parameters[preio->index_helium_fullreio_width];
        /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
        *xe += preio->reionization_parameters[preio->index_helium_fullreio_fraction]
          * (tanh(argument)+1.)/2.;
      }
    }

    return _SUCCESS_;

  }

  /** - implementation of binned ionization function similar to astro-ph/0606552 */

  if (pth->reio_parametrization == reio_bins_tanh) {

    /** - --> case z > z_reio_start */

    if (z > preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1]) {
      *xe = preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-1];
    }

    else if (z < preio->reionization_parameters[preio->index_reio_first_z]) {
      *xe = preio->reionization_parameters[preio->index_reio_first_xe];
    }

    else {

      i = 0;
      while (preio->reionization_parameters[preio->index_reio_first_z+i+1]<z) i++;

      /* This is the expression of the tanh-like jumps of the
         reio_bins_tanh scheme until the 10.06.2015. It appeared to be
         not robust enough. It could lead to a kink in xe(z) near the
         maximum value of z at which reionisation is sampled. It has
         been replaced by the simpler and more robust expression
         below.

         *xe = preio->reionization_parameters[preio->index_reio_first_xe+i]
         +0.5*(tanh((2.*(z-preio->reionization_parameters[preio->index_reio_first_z+i])
         /(preio->reionization_parameters[preio->index_reio_first_z+i+1]
         -preio->reionization_parameters[preio->index_reio_first_z+i])-1.)
         /preio->reionization_parameters[preio->index_reio_step_sharpness])
         /tanh(1./preio->reionization_parameters[preio->index_reio_step_sharpness])+1.)
         *(preio->reionization_parameters[preio->index_reio_first_xe+i+1]
         -preio->reionization_parameters[preio->index_reio_first_xe+i]);
      */

      /* compute the central redshift value of the tanh jump */

      if (i == preio->reio_num_z-2) {
        z_jump = preio->reionization_parameters[preio->index_reio_first_z+i]
          + 0.5*(preio->reionization_parameters[preio->index_reio_first_z+i]
                 -preio->reionization_parameters[preio->index_reio_first_z+i-1]);
      }
      else  {
        z_jump =  0.5*(preio->reionization_parameters[preio->index_reio_first_z+i+1]
                       + preio->reionization_parameters[preio->index_reio_first_z+i]);
      }

      /* implementation of the tanh jump */

      *xe = preio->reionization_parameters[preio->index_reio_first_xe+i]
        +0.5*(tanh((z-z_jump)
                   /preio->reionization_parameters[preio->index_reio_step_sharpness])+1.)
        *(preio->reionization_parameters[preio->index_reio_first_xe+i+1]
          -preio->reionization_parameters[preio->index_reio_first_xe+i]);

    }

    return _SUCCESS_;

  }

  /** - implementation of many tanh jumps */

  if (pth->reio_parametrization == reio_many_tanh) {

    /** - --> case z > z_reio_start */

    if (z > preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1]) {
      *xe = preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-1];
    }

    else if (z > preio->reionization_parameters[preio->index_reio_first_z]) {

      *xe = preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-1];

      for (jump=1; jump<preio->reio_num_z-1; jump++){

        center = preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1-jump];
        // before and after are meant with respect to growing z, not growing time
        before = preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-1-jump]
          -preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-jump];
        after = 0.;
        width = preio->reionization_parameters[preio->index_reio_step_sharpness];

        class_call(thermodynamics_tanh(z,center,before,after,width,&one_jump),
                   error_message_,
                   error_message_);

        *xe += one_jump;

      }

    }

    else {
      *xe = preio->reionization_parameters[preio->index_reio_first_xe];
    }

    return _SUCCESS_;

  }

  /** - implementation of reio_inter */

  if (pth->reio_parametrization == reio_inter) {

    /** - --> case z > z_reio_start */

    if (z > preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1]) {
      *xe = preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-1];
      class_stop(error_message_, "Check: is it normal that we are here?");
    }

    else {

      i=0;
      while (preio->reionization_parameters[preio->index_reio_first_z+i+1] < z) i++;

      double z_min = preio->reionization_parameters[preio->index_reio_first_z+i];
      double z_max = preio->reionization_parameters[preio->index_reio_first_z+i+1];

      class_test(z<z_min,
                 error_message_,
                 "");

      class_test(z>z_max,
                 error_message_,
                 "");

      double x=(z-preio->reionization_parameters[preio->index_reio_first_z+i])
        /(preio->reionization_parameters[preio->index_reio_first_z+i+1]
          -preio->reionization_parameters[preio->index_reio_first_z+i]);

      *xe = preio->reionization_parameters[preio->index_reio_first_xe+i]
        + x*(preio->reionization_parameters[preio->index_reio_first_xe+i+1]
             -preio->reionization_parameters[preio->index_reio_first_xe+i]);

      class_test(*xe<0.,
                 error_message_,
                 "%e %e %e\n",
                 x,
                 preio->reionization_parameters[preio->index_reio_first_xe+i],
                 preio->reionization_parameters[preio->index_reio_first_xe+i+1]);

    }

    return _SUCCESS_;

  }

  class_test(0 == 0,
             error_message_,
             "value of reio_parametrization=%d unclear",pth->reio_parametrization);
}

/**
 * This subroutine reads \f$ X_e(z) \f$ in the recombination table at
 * the time at which reionization starts. Hence it provides correct
 * initial conditions for the reionization function.
 *
 * @param ppr   Input: pointer to precision structure
 * @param pth   Input: pointer to thermo structure
 * @param preco Input: pointer to recombination structure
 * @param z     Input: redshift z_reio_start
 * @param xe    Output: \f$ X_e(z) \f$ at z
 */

int ThermodynamicsModule::thermodynamics_get_xe_before_reionization(recombination* preco, double z, double* xe) {

  int last_index=0;

  class_call(array_interpolate_one_growing_closeby(preco->recombination_table,
                                                   preco->re_size,
                                                   preco->rt_size,
                                                   preco->index_re_z,
                                                   z,
                                                   &last_index,
                                                   preco->index_re_xe,
                                                   xe,
                                                   error_message_),
             error_message_,
             error_message_);

  return _SUCCESS_;

}


/**
 * This routine computes the reionization history. In the reio_camb
 * scheme, this is straightforward if the input parameter is the
 * reionization redshift. If the input is the optical depth, need to
 * find z_reio by dichotomy (trying several z_reio until the correct
 * tau_reio is approached).
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param pth Input: pointer to thermo structure
 * @param preco Input: pointer to filled recombination structure
 * @param preio Input/Output: pointer to reionization structure (to be filled)
 * @param pvecback   Input: vector of background quantities (used as workspace: must be already allocated, with format short_info or larger, but does not need to be filled)
 * @return the error status
 */

int ThermodynamicsModule::thermodynamics_reionization(recombination* preco, reionization* preio, double* pvecback) {

  /** Summary: */

  /** - define local variables */

  int counter;
  double z_sup,z_mid,z_inf;
  double tau_sup,tau_mid,tau_inf;
  int bin;
  int point;
  double xe_input,xe_actual;

  /** - allocate the vector of parameters defining the function \f$ X_e(z) \f$ */

  class_alloc(preio->reionization_parameters, preio->reio_num_params*sizeof(double), error_message_);

  /** - (a) if reionization implemented like in CAMB */

  if ((pth->reio_parametrization == reio_camb) || (pth->reio_parametrization == reio_half_tanh)) {

    /** - --> set values of these parameters, excepted those depending on the reionization redshift */

    if (pth->reio_parametrization == reio_camb) {
      preio->reionization_parameters[preio->index_reio_xe_after] = 1. + YHe_/(_not4_*(1. - YHe_));    /* xe_after_reio: H + singly ionized He (note: segmentation fault impossible, checked before that denominator is non-zero) */
    }
    if (pth->reio_parametrization == reio_half_tanh) {
      preio->reionization_parameters[preio->index_reio_xe_after] = 1.; /* xe_after_reio: neglect He ionization */
      //+ 2*YHe_/(_not4_*(1.-YHe_));    /* xe_after_reio: H + fully ionized He */
    }
    preio->reionization_parameters[preio->index_reio_exponent] = pth->reionization_exponent; /* reio_exponent */
    preio->reionization_parameters[preio->index_reio_width] = pth->reionization_width;    /* reio_width */
    preio->reionization_parameters[preio->index_helium_fullreio_fraction] = YHe_/(_not4_*(1. - YHe_)); /* helium_fullreio_fraction (note: segmentation fault impossible, checked before that denominator is non-zero) */
    preio->reionization_parameters[preio->index_helium_fullreio_redshift] = pth->helium_fullreio_redshift; /* helium_fullreio_redshift */
    preio->reionization_parameters[preio->index_helium_fullreio_width] = pth->helium_fullreio_width;    /* helium_fullreio_width */

    class_test(preio->reionization_parameters[preio->index_reio_exponent]==0,
               error_message_,
               "stop to avoid division by zero");

    class_test(preio->reionization_parameters[preio->index_reio_width]==0,
               error_message_,
               "stop to avoid division by zero");

    class_test(preio->reionization_parameters[preio->index_helium_fullreio_width]==0,
               error_message_,
               "stop to avoid division by zero");

    /** - --> if reionization redshift given as an input, initialize the remaining values and fill reionization table*/

    if (pth->reio_z_or_tau == reio_z) {

      /* reionization redshift */
      preio->reionization_parameters[preio->index_reio_redshift] = z_reionization_;

      /* infer starting redshift for hydrogen */

      if (pth->reio_parametrization == reio_camb) {

        preio->reionization_parameters[preio->index_reio_start] = preio->reionization_parameters[preio->index_reio_redshift]+ppr->reionization_start_factor*pth->reionization_width;

        /* if starting redshift for helium is larger, take that one
           (does not happen in realistic models) */
        if (preio->reionization_parameters[preio->index_reio_start] <
            pth->helium_fullreio_redshift+ppr->reionization_start_factor*pth->helium_fullreio_width)

          preio->reionization_parameters[preio->index_reio_start] =
            pth->helium_fullreio_redshift+ppr->reionization_start_factor*pth->helium_fullreio_width;

      }
      else {

        preio->reionization_parameters[preio->index_reio_start] = z_reionization_;
      }

      class_test(preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max,
                 error_message_,
                 "starting redshift for reionization > reionization_z_start_max = %e\n",ppr->reionization_z_start_max);

      /* infer xe_before_reio */
      class_call(thermodynamics_get_xe_before_reionization(preco,
                                                           preio->reionization_parameters[preio->index_reio_start],
                                                           &(preio->reionization_parameters[preio->index_reio_xe_before])),
                 error_message_,
                 error_message_);

      /* fill reionization table */
      class_call(thermodynamics_reionization_sample(preco, preio, pvecback),
                 error_message_,
                 error_message_);

      tau_reionization_ = preio->reionization_optical_depth;

    }

    /** - --> if reionization optical depth given as an input, find reionization redshift by dichotomy and initialize the remaining values */

    if (pth->reio_z_or_tau == reio_tau) {

      /* upper value */

      z_sup = ppr->reionization_z_start_max-ppr->reionization_start_factor*pth->reionization_width;
      class_test(z_sup < 0.,
                 error_message_,
                 "parameters are such that reionization cannot take place before today while starting after z_start_max; need to increase z_start_max");

      /* maximum possible reionization redshift */
      preio->reionization_parameters[preio->index_reio_redshift] = z_sup;
      /* maximum possible starting redshift */
      preio->reionization_parameters[preio->index_reio_start] = ppr->reionization_z_start_max;
      /* infer xe_before_reio */
      class_call(thermodynamics_get_xe_before_reionization(preco,
                                                           preio->reionization_parameters[preio->index_reio_start],
                                                           &(preio->reionization_parameters[preio->index_reio_xe_before])),
                 error_message_,
                 error_message_);

      /* fill reionization table */
      class_call(thermodynamics_reionization_sample(preco, preio, pvecback),
                 error_message_,
                 error_message_);

      tau_sup=preio->reionization_optical_depth;

      class_test_except(tau_sup < tau_reionization_,
                        error_message_,
                        free(preio->reionization_parameters);free(preio->reionization_table),
                        "parameters are such that reionization cannot start after z_start_max");

      /* lower value */

      z_inf = 0.;
      tau_inf = 0.;

      /* try intermediate values */

      counter=0;
      while ((tau_sup - tau_inf) > tau_reionization_*ppr->reionization_optical_depth_tol) {
        z_mid=0.5*(z_sup+z_inf);

        /* reionization redshift */
        preio->reionization_parameters[preio->index_reio_redshift] = z_mid;
        /* infer starting redshift for hygrogen */
        preio->reionization_parameters[preio->index_reio_start] = preio->reionization_parameters[preio->index_reio_redshift]+ppr->reionization_start_factor*pth->reionization_width;
        /* if starting redshift for helium is larger, take that one
           (does not happen in realistic models) */
        if (preio->reionization_parameters[preio->index_reio_start] <
            pth->helium_fullreio_redshift+ppr->reionization_start_factor*pth->helium_fullreio_width)

          preio->reionization_parameters[preio->index_reio_start] =
            pth->helium_fullreio_redshift+ppr->reionization_start_factor*pth->helium_fullreio_width;

        class_test(preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max,
                   error_message_,
                   "starting redshift for reionization > reionization_z_start_max = %e",ppr->reionization_z_start_max);

        /* infer xe_before_reio */
        class_call(thermodynamics_get_xe_before_reionization(preco,
                                                             preio->reionization_parameters[preio->index_reio_start],
                                                             &(preio->reionization_parameters[preio->index_reio_xe_before])),
                   error_message_,
                   error_message_);

        /* clean and fill reionization table */
        free(preio->reionization_table);
        class_call(thermodynamics_reionization_sample(preco, preio, pvecback),
                   error_message_,
                   error_message_);

        tau_mid=preio->reionization_optical_depth;

        /* trial */

        if (tau_mid > tau_reionization_) {
          z_sup=z_mid;
          tau_sup=tau_mid;
        }
        else {
          z_inf=z_mid;
          tau_inf=tau_mid;
        }

        counter++;
        class_test(counter > _MAX_IT_,
                   error_message_,
                   "while searching for reionization_optical_depth, maximum number of iterations exceeded");
      }

      /* store z_reionization in thermodynamics structure */
      z_reionization_ = preio->reionization_parameters[preio->index_reio_redshift];

    }

    free(preio->reionization_parameters);

    return _SUCCESS_;

  }

  /** - (b) if reionization implemented with reio_bins_tanh scheme */

  if (pth->reio_parametrization == reio_bins_tanh) {

    /* this algorithm requires at least two bin centers (i.e. at least
       4 values in the (z,xe) array, counting the edges). */
    class_test(pth->binned_reio_num<2,
               error_message_,
               "current implementation of binned reio requires at least two bin centers");

    /* check that this input can be interpreted by the code */
    for (bin=1; bin<pth->binned_reio_num; bin++) {
      class_test(pth->binned_reio_z[bin-1]>=pth->binned_reio_z[bin],
                 error_message_,
                 "value of reionization bin centers z_i expected to be passed in growing order: %e, %e",
                 pth->binned_reio_z[bin-1],
                 pth->binned_reio_z[bin]);
    }

    /* the code will not only copy here the "bin centers" passed in
       input. It will add an initial and final value for (z,xe).
       First, fill all entries except the first and the last */

    for (bin=1; bin<preio->reio_num_z-1; bin++) {
      preio->reionization_parameters[preio->index_reio_first_z+bin] = pth->binned_reio_z[bin-1];
      preio->reionization_parameters[preio->index_reio_first_xe+bin] = pth->binned_reio_xe[bin-1];
    }


    /* find largest value of z in the array. We choose to define it as
       z_(i_max) + 2*(the distance between z_(i_max) and z_(i_max-1)). E.g. if
       the bins are in 10,12,14, the largest z will be 18. */
    preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1] =
      preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-2]
      +2.*(preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-2]
           -preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-3]);

    /* copy this value in reio_start */
    preio->reionization_parameters[preio->index_reio_start] = preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1];

    /* check it's not too big */
    class_test(preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max,
               error_message_,
               "starting redshift for reionization = %e, reionization_z_start_max = %e, you must change the binning or increase reionization_z_start_max",
               preio->reionization_parameters[preio->index_reio_start],
               ppr->reionization_z_start_max);

    /* find smallest value of z in the array. We choose
       to define it as z_0 - (the distance between z_1 and z_0). E.g. if
       the bins are in 10,12,14, the stop redshift will be 8. */

    preio->reionization_parameters[preio->index_reio_first_z] =
      2.*preio->reionization_parameters[preio->index_reio_first_z+1]
      -preio->reionization_parameters[preio->index_reio_first_z+2];

    /* check it's not too small */
    /* 6.06.2015: changed this test to simply imposing that the first z is at least zero */
    /*
      class_test(preio->reionization_parameters[preio->index_reio_first_z] < 0,
      error_message_,
      "final redshift for reionization = %e, you must change the binning or redefine the way in which the code extrapolates below the first value of z_i",preio->reionization_parameters[preio->index_reio_first_z]);
    */
    if (preio->reionization_parameters[preio->index_reio_first_z] < 0) {
      preio->reionization_parameters[preio->index_reio_first_z] = 0.;
    }

    /* infer xe before reio */
    class_call(thermodynamics_get_xe_before_reionization(preco,
                                                         preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1],
                                                         &(preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-1])),
               error_message_,
               error_message_);

    /* infer xe after reio */
    preio->reionization_parameters[preio->index_reio_first_xe] = 1. + YHe_/(_not4_*(1. - YHe_));    /* xe_after_reio: H + singly ionized He (note: segmentation fault impossible, checked before that denominator is non-zero) */

    /* pass step sharpness parameter */
    preio->reionization_parameters[preio->index_reio_step_sharpness] = pth->binned_reio_step_sharpness;

    /* fill reionization table */
    class_call(thermodynamics_reionization_sample(preco, preio, pvecback),
               error_message_,
               error_message_);

    tau_reionization_ = preio->reionization_optical_depth;

    return _SUCCESS_;

  }

  /** - (c) if reionization implemented with reio_many_tanh scheme */

  if (pth->reio_parametrization == reio_many_tanh) {

    /* this algorithm requires at least one jump centers */
    class_test(pth->many_tanh_num<1,
               error_message_,
               "current implementation of reio_many_tanh requires at least one jump center");

    /* check that z input can be interpreted by the code */
    for (bin=1; bin<pth->many_tanh_num; bin++) {
      class_test(pth->many_tanh_z[bin-1]>=pth->many_tanh_z[bin],
                 error_message_,
                 "value of reionization bin centers z_i expected to be passed in growing order: %e, %e",
                 pth->many_tanh_z[bin-1],
                 pth->many_tanh_z[bin]);
    }

    /* the code will not only copy here the "jump centers" passed in
       input. It will add an initial and final value for (z,xe).
       First, fill all entries except the first and the last */

    for (bin=1; bin<preio->reio_num_z-1; bin++) {

      preio->reionization_parameters[preio->index_reio_first_z+bin] = pth->many_tanh_z[bin-1];

      /* check that xe input can be interpreted by the code */
      xe_input = pth->many_tanh_xe[bin-1];
      if (xe_input >= 0.) {
        xe_actual = xe_input;
      }
      //-1 means "after hydrogen + first helium recombination"
      else if ((xe_input<-0.9) && (xe_input>-1.1)) {
        xe_actual = 1. + YHe_/(_not4_*(1. - YHe_));
      }
      //-2 means "after hydrogen + second helium recombination"
      else if ((xe_input<-1.9) && (xe_input>-2.1)) {
        xe_actual = 1. + 2.*YHe_/(_not4_*(1. - YHe_));
      }
      //other negative number is nonsense
      else {
        class_stop(error_message_,
                   "Your entry for many_tanh_xe[%d] is %e, this makes no sense (either positive or 0,-1,-2)",
                   bin-1,pth->many_tanh_xe[bin-1]);
      }

      preio->reionization_parameters[preio->index_reio_first_xe+bin] = xe_actual;
    }

    /* find largest value of z in the array. We choose to define it as
       z_(i_max) + ppr->reionization_start_factor*step_sharpness. */
    preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1] =
      preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-2]
      +ppr->reionization_start_factor*pth->many_tanh_width;

    /* copy this value in reio_start */
    preio->reionization_parameters[preio->index_reio_start] = preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1];

    /* check it's not too big */
    class_test(preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max,
               error_message_,
               "starting redshift for reionization = %e, reionization_z_start_max = %e, you must change the binning or increase reionization_z_start_max",
               preio->reionization_parameters[preio->index_reio_start],
               ppr->reionization_z_start_max);

    /* find smallest value of z in the array. We choose
       to define it as z_0 - ppr->reionization_start_factor*step_sharpness, but at least zero. */

    preio->reionization_parameters[preio->index_reio_first_z] =
      preio->reionization_parameters[preio->index_reio_first_z+1]
      -ppr->reionization_start_factor*pth->many_tanh_width;

    if (preio->reionization_parameters[preio->index_reio_first_z] < 0) {
      preio->reionization_parameters[preio->index_reio_first_z] = 0.;
    }

    /* infer xe before reio */
    class_call(thermodynamics_get_xe_before_reionization(preco,
                                                         preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1],
                                                         &(preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-1])),
               error_message_,
               error_message_);

    /* infer xe after reio */

    preio->reionization_parameters[preio->index_reio_first_xe] = preio->reionization_parameters[preio->index_reio_first_xe+1];

    /* if we want to model only hydrogen reionization and neglect both helium reionization */
    //preio->reionization_parameters[preio->index_reio_first_xe] = 1.;

    /* if we want to model only hydrogen + first helium reionization and neglect second helium reionization */
    //preio->reionization_parameters[preio->index_reio_first_xe] = 1. + YHe_/(_not4_*(1.-YHe_));

    /* if we want to model hydrogen + two helium reionization */
    //preio->reionization_parameters[preio->index_reio_first_xe] = 1. + 2.*YHe_/(_not4_*(1.-YHe_));

    /* pass step sharpness parameter */
    class_test(pth->many_tanh_width<=0,
               error_message_,
               "many_tanh_width must be strictly positive, you passed %e",
               pth->many_tanh_width);

    preio->reionization_parameters[preio->index_reio_step_sharpness] = pth->many_tanh_width;

    /* fill reionization table */
    class_call(thermodynamics_reionization_sample(preco, preio, pvecback),
               error_message_,
               error_message_);

    tau_reionization_ = preio->reionization_optical_depth;

    return _SUCCESS_;

  }

  /** - (d) if reionization implemented with reio_inter scheme */

  if (pth->reio_parametrization == reio_inter) {

    /* this parametrization requires at least one point (z,xe) */
    class_test(pth->reio_inter_num<1,
               error_message_,
               "current implementation of reio_inter requires at least one point (z,xe)");

    /* this parametrization requires that the first z value is zero */
    class_test(pth->reio_inter_z[0] != 0.,
               error_message_,
               "For reio_inter scheme, the first value of reio_inter_z[...]  should always be zero, you passed %e",
               pth->reio_inter_z[0]);

    /* check that z input can be interpreted by the code */
    for (point=1; point<pth->reio_inter_num; point++) {
      class_test(pth->reio_inter_z[point-1]>=pth->reio_inter_z[point],
                 error_message_,
                 "value of reionization bin centers z_i expected to be passed in growing order, unlike: %e, %e",
                 pth->reio_inter_z[point-1],
                 pth->reio_inter_z[point]);
    }

    /* this parametrization requires that the last x_i value is zero
       (the code will substitute it with the value that one would get in
       absence of reionization, as compute by the recombination code) */
    class_test(pth->reio_inter_xe[pth->reio_inter_num-1] != 0.,
               error_message_,
               "For reio_inter scheme, the last value of reio_inter_xe[...]  should always be zero, you passed %e",
               pth->reio_inter_xe[pth->reio_inter_num-1]);

    /* copy here the (z,xe) values passed in input. */

    for (point=0; point<preio->reio_num_z; point++) {

      preio->reionization_parameters[preio->index_reio_first_z + point] = pth->reio_inter_z[point];

      /* check that xe input can be interpreted by the code */
      xe_input = pth->reio_inter_xe[point];
      if (xe_input >= 0.) {
        xe_actual = xe_input;
      }
      //-1 means "after hydrogen + first helium recombination"
      else if ((xe_input<-0.9) && (xe_input>-1.1)) {
        xe_actual = 1. + YHe_/(_not4_*(1. - YHe_));
      }
      //-2 means "after hydrogen + second helium recombination"
      else if ((xe_input<-1.9) && (xe_input>-2.1)) {
        xe_actual = 1. + 2.*YHe_/(_not4_*(1. - YHe_));
      }
      //other negative number is nonsense
      else {
        class_stop(error_message_,
                   "Your entry for reio_inter_xe[%d] is %e, this makes no sense (either positive or 0,-1,-2)",
                   point,pth->reio_inter_xe[point]);
      }

      preio->reionization_parameters[preio->index_reio_first_xe+point] = xe_actual;
    }

    /* copy highest redshift in reio_start */
    preio->reionization_parameters[preio->index_reio_start] = preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1];

    /* check it's not too big */
    class_test(preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max,
               error_message_,
               "starting redshift for reionization = %e, reionization_z_start_max = %e, you must change the binning or increase reionization_z_start_max",
               preio->reionization_parameters[preio->index_reio_start],
               ppr->reionization_z_start_max);

    /* infer xe before reio */
    class_call(thermodynamics_get_xe_before_reionization(preco,
                                                         preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1],
                                                         &(preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-1])),
               error_message_,
               error_message_);

    /* fill reionization table */
    class_call(thermodynamics_reionization_sample(preco, preio, pvecback),
               error_message_,
               error_message_);

    tau_reionization_ = preio->reionization_optical_depth;

    return _SUCCESS_;

  }

  class_test(0 == 0,
             error_message_,
             "value of reio_z_or_tau=%d unclear",pth->reio_z_or_tau);

}

/**
 * For fixed input reionization parameters, this routine computes the
 * reionization history and fills the reionization table.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param pth Input: pointer to thermo structure
 * @param preco Input: pointer to filled recombination structure
 * @param preio Input/Output: pointer to reionization structure (to be filled)
 * @param pvecback   Input: vector of background quantities (used as workspace: must be already allocated, with format short_info or larger, but does not need to be filled)
 * @return the error status
 */

int ThermodynamicsModule::thermodynamics_reionization_sample(recombination* preco, reionization* preio, double* pvecback) {

  /** Summary: */

  /** - define local variables */

  /* a growing table (since the number of redshift steps is not known a priori) */
  growTable gTable;
  /* needed for growing table */
  double * pData;
  /* needed for growing table */
  void * memcopy_result;
  /* current vector of values related to reionization */
  double * reio_vector;
  /* running index inside thermodynamics table */
  int i;
  int number_of_redshifts;
  /* values of z, dz, X_e */
  double dz,dz_max;
  double z,z_next;
  double xe,xe_next;
  double dkappadz,dkappadz_next;
  double Tb,Yp,dTdz,opacity,mu;
  double dkappadtau,dkappadtau_next;
  double energy_rate;
  double tau;
  double chi_heat;
  int last_index_back;
  double relative_variation;

  Yp = YHe_;

  /** - (a) allocate vector of values related to reionization */
  class_alloc(reio_vector, preio->re_size*sizeof(double), error_message_);

  /** - (b) create a growTable with gt_init() */
  class_call(gt_init(&gTable),
             gTable.error_message,
             error_message_);

  /** - (c) first line is taken from thermodynamics table, just before reionization starts */

  /** - --> look where to start in current thermodynamics table */
  i=0;
  while (preco->recombination_table[i*preco->re_size+preco->index_re_z] < preio->reionization_parameters[preio->index_reio_start]) {
    i++;
    class_test(i == ppr->recfast_Nz0,
               error_message_,
               "reionization_z_start_max = %e > largest redshift in thermodynamics table",ppr->reionization_z_start_max);
  }

  /** - --> get redshift */
  z=preco->recombination_table[i*preco->re_size+preco->index_re_z];
  reio_vector[preio->index_re_z]=z;
  preio->index_reco_when_reio_start=i;

  /** - --> get \f$ X_e \f$ */
  class_call(thermodynamics_reionization_function(z, preio, &xe),
             error_message_,
             error_message_);

  reio_vector[preio->index_re_xe] = xe;

  /** -  --> get \f$ d \kappa / d z = (d \kappa / d \tau) * (d \tau / d z) = - (d \kappa / d \tau) / H \f$ */

  class_call(background_module_->background_tau_of_z(z, &tau),
             background_module_->error_message_,
             error_message_);

  class_call(background_module_->background_at_tau(
                                                  tau,
                                                  pba->short_info,
                                                  pba->inter_normal,
                                                  &last_index_back,
                                                  pvecback),
             background_module_->error_message_,
             error_message_);

  reio_vector[preio->index_re_dkappadtau] = (1. + z)*(1. + z)*n_e_*xe*_sigma_*_Mpc_over_m_;

  class_test(pvecback[background_module_->index_bg_H_] == 0.,
             error_message_,
             "stop to avoid division by zero");

  reio_vector[preio->index_re_dkappadz] = reio_vector[preio->index_re_dkappadtau]/pvecback[background_module_->index_bg_H_];

  dkappadz = reio_vector[preio->index_re_dkappadz];
  dkappadtau = reio_vector[preio->index_re_dkappadtau];

  /** - --> get baryon temperature **/
  Tb = preco->recombination_table[i*preco->re_size+preco->index_re_Tb];
  reio_vector[preio->index_re_Tb] = Tb;

  /** - --> after recombination, Tb scales like (1+z)**2. Compute constant factor Tb/(1+z)**2. */
  //Tba2 = Tb/(1+z)/(1+z);

  /** - --> get baryon equation of state */
  reio_vector[preio->index_re_wb] = _k_B_ / ( _c_ * _c_ * _m_H_ ) * (1. + (1./_not4_ - 1.) * Yp + xe * (1.-Yp)) * Tb;

  /** - --> get baryon adiabatic sound speed */
  reio_vector[preio->index_re_cb2] = 5./3. * reio_vector[preio->index_re_wb];

  /** - --> store these values in growing table */
  class_call(gt_add(&gTable,_GT_END_,(void *) reio_vector,sizeof(double)*(preio->re_size)),
             gTable.error_message,
             error_message_);

  number_of_redshifts=1;

  /** - (d) set the maximum step value (equal to the step in thermodynamics table) */
  dz_max=preco->recombination_table[i*preco->re_size+preco->index_re_z]
    -preco->recombination_table[(i-1)*preco->re_size+preco->index_re_z];

  /** - (e) loop over redshift values in order to find values of z, x_e, kappa' (Tb and cb2 found later by integration). The sampling in z space is found here. */

  /* initial step */
  dz = dz_max;

  while (z > 0.) {

    class_test(dz < ppr->smallest_allowed_variation,
               error_message_,
               "stuck in the loop for reionization sampling, as if you were trying to impose a discontinuous evolution for xe(z)");

    /* - try next step */
    z_next=z-dz;
    if (z_next < 0.) z_next=0.;

    class_call(thermodynamics_reionization_function(z_next, preio, &xe_next),
               error_message_,
               error_message_);

    class_call(background_module_->background_tau_of_z(z_next, &tau),
               background_module_->error_message_,
               error_message_);

    class_call(background_module_->background_at_tau(
                                                    tau,
                                                    pba->short_info,
                                                    pba->inter_normal,
                                                    &last_index_back,
                                                    pvecback),
               background_module_->error_message_,
               error_message_);

    class_test(pvecback[background_module_->index_bg_H_] == 0.,
               error_message_,
               "stop to avoid division by zero");

    dkappadz_next = (1. + z_next)*(1. + z_next)*n_e_*xe_next*_sigma_*_Mpc_over_m_/pvecback[background_module_->index_bg_H_];

    dkappadtau_next = (1. + z_next)*(1. + z_next)*n_e_*xe_next*_sigma_*_Mpc_over_m_;

    class_test((dkappadz == 0.) || (dkappadtau == 0.),
               error_message_,
               "stop to avoid division by zero");

    relative_variation = fabs((dkappadz_next-dkappadz)/dkappadz) +
      fabs((dkappadtau_next-dkappadtau)/dkappadtau);

    if (relative_variation < ppr->reionization_sampling) {
      /* accept the step: get \f$ z, X_e, d kappa / d z \f$ and store in growing table */

      z=z_next;
      xe=xe_next;
      dkappadz=dkappadz_next;
      dkappadtau= dkappadtau_next;

      class_test((dkappadz == 0.) || (dkappadtau == 0.),
                 error_message_,
                 "dkappadz=%e, dkappadtau=%e, stop to avoid division by zero",dkappadz,dkappadtau);

      reio_vector[preio->index_re_z] = z;
      reio_vector[preio->index_re_xe] = xe;
      reio_vector[preio->index_re_dkappadz] = dkappadz;
      reio_vector[preio->index_re_dkappadtau] = dkappadz*pvecback[background_module_->index_bg_H_];

      class_call(gt_add(&gTable,_GT_END_,(void *) reio_vector,sizeof(double)*(preio->re_size)),
                 gTable.error_message,
                 error_message_);

      number_of_redshifts++;

      dz = MIN(0.9*(ppr->reionization_sampling/relative_variation),5.)*dz;
      dz = MIN(dz,dz_max);
    }
    else {
      /* do not accept the step and update dz */
      dz = 0.9*(ppr->reionization_sampling/relative_variation)*dz;
    }
  }

  /** - (f) allocate reionization_table with correct size */
  class_alloc(preio->reionization_table, preio->re_size*number_of_redshifts*sizeof(double), error_message_);

  preio->rt_size=number_of_redshifts;

  /** - (g) retrieve data stored in the growTable with gt_getPtr() */
  class_call(gt_getPtr(&gTable,(void**)&pData),
             gTable.error_message,
             error_message_);

  /** - (h) copy growTable to reionization_temporary_table (invert order of lines, so that redshift is growing, like in recombination table) */
  for (i=0; i < preio->rt_size; i++) {
    memcopy_result = memcpy(preio->reionization_table+i*preio->re_size,pData+(preio->rt_size-i-1)*preio->re_size,preio->re_size*sizeof(double));
    class_test(memcopy_result != preio->reionization_table+i*preio->re_size,
               error_message_,
               "cannot copy data back to reionization_temporary_table");

  }

  /** - (i) free the growTable with gt_free() , free vector of reionization variables */
  class_call(gt_free(&gTable),
             gTable.error_message,
             error_message_);

  free(reio_vector);

  /** - (j) another loop on z, to integrate equation for Tb and to compute cb2 */
  for (i=preio->rt_size-1; i >0 ; i--) {

    z = preio->reionization_table[i*preio->re_size+preio->index_re_z];

    class_call(background_module_->background_tau_of_z(z, &tau),
               background_module_->error_message_,
               error_message_);

    class_call(background_module_->background_at_tau(
                                                    tau,
                                                    pba->normal_info,
                                                    pba->inter_normal,
                                                    &last_index_back,
                                                    pvecback),
               background_module_->error_message_,
               error_message_);

    dz = (preio->reionization_table[i*preio->re_size+preio->index_re_z]-preio->reionization_table[(i-1)*preio->re_size+preio->index_re_z]);

    opacity = (1. + z)*(1. + z)*n_e_*preio->reionization_table[i*preio->re_size + preio->index_re_xe]*_sigma_*_Mpc_over_m_;

    mu = _m_H_/(1. + (1./_not4_ - 1.)*YHe_ + preio->reionization_table[i*preio->re_size + preio->index_re_xe]*(1. - YHe_));

    /** - --> derivative of baryon temperature */

    dTdz=2./(1+z)*preio->reionization_table[i*preio->re_size+preio->index_re_Tb]
      - 2.*mu/_m_e_*4.*pvecback[background_module_->index_bg_rho_g_]/3./pvecback[background_module_->index_bg_rho_b_]*opacity*
      (pba->T_cmb*(1. + z) - preio->reionization_table[i*preio->re_size + preio->index_re_Tb])/pvecback[background_module_->index_bg_H_];

    if (preco->annihilation > 0) {

      class_call(thermodynamics_energy_injection(preco, z, &energy_rate, error_message_),
                 error_message_,
                 error_message_);

      // old approximation from Chen and Kamionkowski:
      // chi_heat = (1.+2.*preio->reionization_table[i*preio->re_size+preio->index_re_xe])/3.;

      // coefficient as revised by Slatyer et al. 2013
      // (in fact it is a fit by Vivian Poulin of columns 1 and 2 in Table V
      // of Slatyer et al. 2013):
      xe = preio->reionization_table[i*preio->re_size+preio->index_re_xe];
      if (xe < 1.)
        chi_heat = MIN(0.996857*(1.-pow(1.-pow(xe,0.300134),1.51035)),1);
      else
        chi_heat = 1.;

      dTdz+= -2./(3.*_k_B_)*energy_rate*chi_heat
        /(preco->Nnow*pow(1.+z,3))/(1.+preco->fHe+preio->reionization_table[i*preio->re_size+preio->index_re_xe])
        /(pvecback[background_module_->index_bg_H_]*_c_/_Mpc_over_m_*(1. + z)); /* energy injection */

    }

    /** - --> increment baryon temperature */

    preio->reionization_table[(i-1)*preio->re_size+preio->index_re_Tb] =
      preio->reionization_table[i*preio->re_size+preio->index_re_Tb]-dTdz*dz;

    /** - --> get baryon equation of state */

    preio->reionization_table[(i-1)*preio->re_size+preio->index_re_wb] =
      _k_B_/ ( _c_ * _c_ * mu)
      * preio->reionization_table[(i-1)*preio->re_size+preio->index_re_Tb];

    /** - --> get baryon adiabatic sound speed */

    preio->reionization_table[(i-1)*preio->re_size+preio->index_re_cb2] =
      preio->reionization_table[(i-1)*preio->re_size+preio->index_re_wb]
      *(1.+(1+z)/3.*dTdz/preio->reionization_table[(i-1)*preio->re_size+preio->index_re_Tb]);
  }

  /** - --> spline \f$ d \tau / dz \f$ with respect to z in view of integrating for optical depth */
  class_call(array_spline(preio->reionization_table,
                          preio->re_size,
                          preio->rt_size,
                          preio->index_re_z,
                          preio->index_re_dkappadz,
                          preio->index_re_d3kappadz3,
                          _SPLINE_EST_DERIV_,
                          error_message_),
             error_message_,
             error_message_);

  /** - --> integrate for optical depth */
  class_call(array_integrate_all_spline(preio->reionization_table,
                                        preio->re_size,
                                        preio->rt_size,
                                        preio->index_re_z,
                                        preio->index_re_dkappadz,
                                        preio->index_re_d3kappadz3,
                                        &(preio->reionization_optical_depth),
                                        error_message_),
             error_message_,
             error_message_);

  return _SUCCESS_;

}

/**
 * Integrate thermodynamics with your favorite recombination code.
 *
 */

int ThermodynamicsModule::thermodynamics_recombination(recombination* preco, double* pvecback) {

  if (pth->recombination==hyrec) {

    class_call(thermodynamics_recombination_with_hyrec(preco, pvecback),
               error_message_,
               error_message_);

  }

  if (pth->recombination==recfast) {

    class_call(thermodynamics_recombination_with_recfast(preco, pvecback),
               error_message_,
               error_message_);

  }

  return _SUCCESS_;

}

/**
 * Integrate thermodynamics with HyRec.
 *
 * Integrate thermodynamics with HyRec, allocate and fill the part
 * of the thermodynamics interpolation table (the rest is filled in
 * thermodynamics_init()). Called once by
 * thermodynamics_recombination(), from thermodynamics_init().
 *
 *************************************************************************************************
 *                 HYREC: Hydrogen and Helium Recombination Code
 *         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)
 *************************************************************************************************
 *
 *
 * @param ppr      Input: pointer to precision structure
 * @param pba      Input: pointer to background structure
 * @param pth      Input: pointer to thermodynamics structure
 * @param preco    Output: pointer to recombination structure
 * @param pvecback Input: pointer to an allocated (but empty) vector of background variables
 */

int ThermodynamicsModule::thermodynamics_recombination_with_hyrec(recombination* preco, double* pvecback) {
  /** Summary: */
#ifdef HYREC

  REC_COSMOPARAMS param;
  HRATEEFF rate_table;
  TWO_PHOTON_PARAMS twog_params;
  double *xe_output, *Tm_output;
  int i,j,l,Nz,b;
  double z, xe, Tm, Hz;
  FILE *fA;
  FILE *fR;
  double L2s1s_current;
  void * buffer;
  int buf_size;
  double tau;
  int last_index_back;
  double w_fld,dw_over_da_fld,integral_fld;

  /** - Fill hyrec parameter structure */

  param.T0 = pba->T_cmb;
  param.obh2 = pba->Omega0_b*pba->h*pba->h;
  param.omh2 = (pba->Omega0_b+pba->Omega0_cdm+pba->Omega0_idm_dr+pba->Omega0_ncdm_tot)*pba->h*pba->h;
  param.okh2 = pba->Omega0_k*pba->h*pba->h;
  param.odeh2 = (pba->Omega0_lambda+pba->Omega0_fld)*pba->h*pba->h;
  class_call(background_module_->background_w_fld(pba->a_today, &w_fld, &dw_over_da_fld, &integral_fld), background_module_->error_message_, error_message_);
  param.w0 = w_fld;
  param.wa = -dw_over_da_fld*pba->a_today;
  param.Y = YHe_;
  param.Nnueff = background_module_->Neff_;
  param.nH0 = 11.223846333047*param.obh2*(1.-param.Y);  /* number density of hydrogen today in m-3 */
  param.fHe = param.Y/(1-param.Y)/3.97153;              /* abundance of helium by number */
  param.zstart = ppr->recfast_z_initial; /* Redshift range */
  param.zend = 0.;
  param.dlna = 8.49e-5;
  param.nz = (long) floor(2+log((1.+param.zstart)/(1.+param.zend))/param.dlna);
  param.annihilation = pth->annihilation;
  param.has_on_the_spot = pth->has_on_the_spot;
  param.decay = pth->decay;
  param.annihilation_variation = pth->annihilation_variation;
  param.annihilation_z = pth->annihilation_z;
  param.annihilation_zmax = pth->annihilation_zmax;
  param.annihilation_zmin = pth->annihilation_zmin;
  param.annihilation_f_halo = pth->annihilation_f_halo;
  param.annihilation_z_halo = pth->annihilation_z_halo;

  /** - Build effective rate tables */

  /* allocate contiguous memory zone */

  buf_size = (2*NTR+NTM+2*NTR*NTM+2*param.nz)*sizeof(double) + 2*NTM*sizeof(double*);

  class_alloc(buffer,
              buf_size,
              error_message_);

  /** - distribute addresses for each table */

  rate_table.logTR_tab = (double*)buffer;
  rate_table.TM_TR_tab = (double*)(rate_table.logTR_tab + NTR);
  rate_table.logAlpha_tab[0] = (double**)(rate_table.TM_TR_tab+NTM);
  rate_table.logAlpha_tab[1] = (double**)(rate_table.logAlpha_tab[0]+NTM);
  rate_table.logAlpha_tab[0][0] = (double*)(rate_table.logAlpha_tab[1]+NTM);
  for (j=1;j<NTM;j++) {
    rate_table.logAlpha_tab[0][j] = (double*)(rate_table.logAlpha_tab[0][j-1]+NTR);
  }
  rate_table.logAlpha_tab[1][0] = (double*)(rate_table.logAlpha_tab[0][NTM-1]+NTR);
  for (j=1;j<NTM;j++) {
    rate_table.logAlpha_tab[1][j] = (double*)(rate_table.logAlpha_tab[1][j-1]+NTR);
  }
  rate_table.logR2p2s_tab = (double*)(rate_table.logAlpha_tab[1][NTM-1]+NTR);

  xe_output = (double*)(rate_table.logR2p2s_tab+NTR);
  Tm_output = (double*)(xe_output+param.nz);

  /* store sampled values of temperatures */

  for (i = 0; i < NTR; i++)
    rate_table.logTR_tab[i] = log(TR_MIN) + i * (log(TR_MAX)-log(TR_MIN))/(NTR-1.);
  for (i = 0; i < NTM; i++)
    rate_table.TM_TR_tab[i] = TM_TR_MIN + i * (TM_TR_MAX-TM_TR_MIN)/(NTM-1.);

  rate_table.DlogTR = rate_table.logTR_tab[1] - rate_table.logTR_tab[0];
  rate_table.DTM_TR = rate_table.TM_TR_tab[1] - rate_table.TM_TR_tab[0];

  /* read in file */

  class_open(fA, ppr->hyrec_Alpha_inf_file, "r", error_message_);
  class_open(fR, ppr->hyrec_R_inf_file, "r", error_message_);

  for (i = 0; i < NTR; i++) {
    for (j = 0; j < NTM; j++) {
      for (l = 0; l <= 1; l++) {
        if (fscanf(fA, "%le", &(rate_table.logAlpha_tab[l][j][i])) != 1)
          class_stop(error_message_, "Error reading hyrec data file %s", ppr->hyrec_Alpha_inf_file);
        rate_table.logAlpha_tab[l][j][i] = log(rate_table.logAlpha_tab[l][j][i]);
      }
    }

    if (fscanf(fR, "%le", &(rate_table.logR2p2s_tab[i])) !=1)
      class_stop(error_message_, "Error reading hyrec data file %s", ppr->hyrec_R_inf_file);
    rate_table.logR2p2s_tab[i] = log(rate_table.logR2p2s_tab[i]);

  }
  fclose(fA);
  fclose(fR);

  /* Read two-photon rate tables */

  class_open(fA, ppr->hyrec_two_photon_tables_file, "r", error_message_);

  for (b = 0; b < NVIRT; b++) {
    if ((fscanf(fA, "%le", &(twog_params.Eb_tab[b])) != 1) ||
        (fscanf(fA, "%le", &(twog_params.A1s_tab[b])) != 1) ||
        (fscanf(fA, "%le", &(twog_params.A2s_tab[b])) != 1) ||
        (fscanf(fA, "%le", &(twog_params.A3s3d_tab[b])) != 1) ||
        (fscanf(fA, "%le", &(twog_params.A4s4d_tab[b])) != 1))
      class_stop(error_message_, "Error reading hyrec data file %s", ppr->hyrec_two_photon_tables_file);
  }

  fclose(fA);

  /** - Normalize 2s--1s differential decay rate to L2s1s (can be set by user in hydrogen.h) */
  L2s1s_current = 0.;
  for (b = 0; b < NSUBLYA; b++) L2s1s_current += twog_params.A2s_tab[b];
  for (b = 0; b < NSUBLYA; b++) twog_params.A2s_tab[b] *= L2s1s/L2s1s_current;

  /*  In CLASS, we have neutralized the switches for the various
      effects considered in Hirata (2008), keeping the full
      calculation as a default; but you could restore their
      functionality by copying a few lines from hyrec/hyrec.c to
      here */

  /** - Compute the recombination history by calling a function in hyrec (no CLASS-like error management here) */

  if (pth->thermodynamics_verbose > 0)
    printf(" -> calling HyRec version %s,\n",HYREC_VERSION);

  rec_build_history(&param, &rate_table, &twog_params, xe_output, Tm_output);

  if (pth->thermodynamics_verbose > 0)
    printf("    by Y. Ali-Haïmoud & C. Hirata\n");

  /** - fill a few parameters in preco and pth */

  Nz=ppr->recfast_Nz0;

  preco->rt_size = Nz;
  preco->H0 = pba->H0 * _c_ / _Mpc_over_m_;
  /* preco->H0 in inverse seconds (while pba->H0 is [H0/c] in inverse Mpcs) */
  preco->YHe = YHe_;
  preco->Nnow = 3.*preco->H0*preco->H0*pba->Omega0_b*(1.-preco->YHe)/(8.*_PI_*_G_*_m_H_);
  /* energy injection parameters */
  preco->annihilation = pth->annihilation;
  preco->has_on_the_spot = pth->has_on_the_spot;
  preco->annihilation_variation = pth->annihilation_variation;
  preco->annihilation_z = pth->annihilation_z;
  preco->annihilation_zmax = pth->annihilation_zmax;
  preco->annihilation_zmin = pth->annihilation_zmin;
  preco->decay = pth->decay;
  preco->annihilation_f_halo = pth->annihilation_f_halo;
  preco->annihilation_z_halo = pth->annihilation_z_halo;
  n_e_ = preco->Nnow;

  /** - allocate memory for thermodynamics interpolation tables (size known in advance) and fill it */

  class_alloc(preco->recombination_table, preco->re_size*preco->rt_size*sizeof(double), error_message_);

  for(i=0; i <Nz; i++) {

    /** - --> get redshift, corresponding results from hyrec, and background quantities */

    z = param.zstart * (1. - (double)(i+1) / (double)Nz);

    /* get (xe,Tm) by interpolating in pre-computed tables */


    class_call(array_interpolate_cubic_equal(-log(1.+param.zstart),
                                             param.dlna,
                                             xe_output,
                                             param.nz,
                                             -log(1.+z),
                                             &xe,
                                             error_message_),
               error_message_,
               error_message_);

    class_call(array_interpolate_cubic_equal(-log(1.+param.zstart),
                                             param.dlna,
                                             Tm_output,
                                             param.nz,
                                             -log(1.+z),
                                             &Tm,
                                             error_message_),
               error_message_,
               error_message_);

    class_call(background_module_->background_tau_of_z(z, &tau),
               background_module_->error_message_,
               error_message_);

    class_call(background_module_->background_at_tau(
                                                    tau,
                                                    pba->short_info,
                                                    pba->inter_normal,
                                                    &last_index_back,
                                                    pvecback),
               background_module_->error_message_,
               error_message_);

    /*   class_call(thermodynamics_energy_injection(ppr,pba,preco,z,&energy_rate,error_message_),
         error_message_,
         error_message_);
    */

    /* Hz is H in inverse seconds (while pvecback returns [H0/c] in inverse Mpcs) */
    Hz = pvecback[background_module_->index_bg_H_]*_c_/_Mpc_over_m_;

    /** - --> store the results in the table */

    /* results are obtained in order of decreasing z, and stored in order of growing z */

    /* redshift */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_z)=z;

    /* ionization fraction */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_xe)=xe;

    /* Tb */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_Tb)=Tm;

    /* wb = (k_B/mu) Tb */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_wb)
      = _k_B_/(_c_*_c_*_m_H_)*(1. + (1./_not4_ - 1.)*YHe_ + xe*(1. - YHe_))*Tm;

    /* cb2 = (k_B/mu) Tb (1-1/3 dlnTb/dlna) = (k_B/mu) Tb (1+1/3 (1+z) dlnTb/dz)
       with (1+z)dlnTb/dz= - [dlnTb/dlna] */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_cb2)
      = *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_wb)
      * (1. - rec_dTmdlna(xe, Tm, pba->T_cmb*(1.+z), Hz, param.fHe, param.nH0*pow((1+z),3)*1e-6, energy_injection_rate(&param,z)) / Tm / 3.);

    /* dkappa/dtau = a n_e x_e sigma_T = a^{-2} n_e(today) x_e sigma_T (in units of 1/Mpc) */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_dkappadtau)
      = (1.+z) * (1.+z) * preco->Nnow * xe * _sigma_ * _Mpc_over_m_;

  }

  /* Cleanup */

  free(buffer);

#else

  class_stop(error_message_,
             "you compiled without including the HyRec code, and now wish to use it. Either set the input parameter 'recombination' to something else than 'HyRec', or recompile after setting in the Makefile the appropriate path HYREC=... ");

#endif

  return _SUCCESS_;
}

/**
 * Integrate thermodynamics with RECFAST.
 *
 * Integrate thermodynamics with RECFAST, allocate and fill the part
 * of the thermodynamics interpolation table (the rest is filled in
 * thermodynamics_init()). Called once by
 * thermodynamics_recombination, from thermodynamics_init().
 *
 *
 *******************************************************************************
 * RECFAST is an integrator for Cosmic Recombination of Hydrogen and Helium,
 * developed by Douglas Scott (dscott@astro.ubc.ca)
 * based on calculations in the paper Seager, Sasselov & Scott
 * (ApJ, 523, L1, 1999).
 * and "fudge" updates in Wong, Moss & Scott (2008).
 *
 * Permission to use, copy, modify and distribute without fee or royalty at
 * any tier, this software and its documentation, for any purpose and without
 * fee or royalty is hereby granted, provided that you agree to comply with
 * the following copyright notice and statements, including the disclaimer,
 * and that the same appear on ALL copies of the software and documentation,
 * including modifications that you make for internal use or for distribution:
 *
 * Copyright 1999-2010 by University of British Columbia.  All rights reserved.
 *
 * THIS SOFTWARE IS PROVIDED "AS IS", AND U.B.C. MAKES NO
 * REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
 * BY WAY OF EXAMPLE, BUT NOT LIMITATION,
 * U.B.C. MAKES NO REPRESENTATIONS OR WARRANTIES OF
 * MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT
 * THE USE OF THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE
 * ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.
 *******************************************************************************
 *
 * Version 1.5: includes extra fitting function from
 *              Rubino-Martin et al. arXiv:0910.4383v1 [astro-ph.CO]
 *
 * @param ppr      Input: pointer to precision structure
 * @param pba      Input: pointer to background structure
 * @param pth      Input: pointer to thermodynamics structure
 * @param preco    Output: pointer to recombination structure
 * @param pvecback Input: pointer to an allocated (but empty) vector of background variables
 * @return the error status
 */

int ThermodynamicsModule::thermodynamics_recombination_with_recfast(recombination* preco, double* pvecback) {

  /** Summary: */

  /** - define local variables */

  /* vector of variables to be integrated: x_H, x_He, Tmat */
  double y[3],dy[3];

  /* other recfast variables */
  double OmegaB,zinitial,x_He0,x0;
  double x_H0=0.;
  double z,mu_H,Lalpha,Lalpha_He,DeltaB,DeltaB_He;
  double zstart,zend,rhs;
  int i,Nz;

  /* introduced by JL for smoothing the various steps */
  double x0_previous,x0_new,s,weight;

  /* contains all quantities relevant for the integration algorithm */
  struct generic_integrator_workspace gi;

  /* contains all fixed parameters which should be passed to thermodynamics_derivs_with_recfast */
  thermodynamics_parameters_and_workspace tpaw{this};

  /** - allocate memory for thermodynamics interpolation tables (size known in advance) */
  preco->rt_size = ppr->recfast_Nz0;
  class_alloc(preco->recombination_table, preco->re_size*preco->rt_size*sizeof(double), error_message_);

  /** - initialize generic integrator with initialize_generic_integrator() */
  class_call(initialize_generic_integrator(_RECFAST_INTEG_SIZE_, &gi),
             gi.error_message,
             error_message_);

  /** - read a few precision/cosmological parameters */

  /* Nz */
  Nz=ppr->recfast_Nz0;

  /* preco->H0 is H0 in inverse seconds (while pba->H0 is [H0/c] in inverse Mpcs) */
  preco->H0 = pba->H0 * _c_ / _Mpc_over_m_;

  /* Omega_b */
  OmegaB = pba->Omega0_b;

  /* Yp */
  preco->YHe = YHe_;

  /* Tnow */
  preco->Tnow = pba->T_cmb;

  /* z_initial */
  zinitial=ppr->recfast_z_initial;

  /* H_frac */
  preco->H_frac = ppr->recfast_H_frac;

  /* H fudging */
  class_test((ppr->recfast_Hswitch != _TRUE_) && (ppr->recfast_Hswitch != _FALSE_),
             error_message_,
             "RECFAST error: unknown H fudging scheme");
  preco->fu = ppr->recfast_fudge_H;
  if (ppr->recfast_Hswitch == _TRUE_)
    preco->fu += ppr->recfast_delta_fudge_H;

  /* He fudging */
  class_test((ppr->recfast_Heswitch < 0) || (ppr->recfast_Heswitch > 6),
             error_message_,
             "RECFAST error: unknown He fudging scheme");

  /* related quantities */
  z=zinitial;
  mu_H = 1./(1.-preco->YHe);
  //mu_T = _not4_ /(_not4_ - (_not4_-1.)*preco->YHe); /* recfast 1.4*/
  preco->fHe = preco->YHe/(_not4_ *(1.-preco->YHe)); /* recfast 1.4 */
  preco->Nnow = 3.*preco->H0*preco->H0*OmegaB/(8.*_PI_*_G_*mu_H*_m_H_);
  n_e_ = preco->Nnow;

  /* energy injection parameters */
  preco->annihilation = pth->annihilation;
  preco->has_on_the_spot = pth->has_on_the_spot;
  preco->annihilation_variation = pth->annihilation_variation;
  preco->annihilation_z = pth->annihilation_z;
  preco->annihilation_zmax = pth->annihilation_zmax;
  preco->annihilation_zmin = pth->annihilation_zmin;
  preco->decay = pth->decay;
  preco->annihilation_f_halo = pth->annihilation_f_halo;
  preco->annihilation_z_halo = pth->annihilation_z_halo;

  /* quantities related to constants defined in thermodynamics.h */
  //n = preco->Nnow * pow((1.+z),3);
  Lalpha = 1./_L_H_alpha_;
  Lalpha_He = 1./_L_He_2p_;
  DeltaB = _h_P_*_c_*(_L_H_ion_-_L_H_alpha_);
  preco->CDB = DeltaB/_k_B_;
  DeltaB_He = _h_P_*_c_*(_L_He1_ion_-_L_He_2s_);
  preco->CDB_He = DeltaB_He/_k_B_;
  preco->CB1 = _h_P_*_c_*_L_H_ion_/_k_B_;
  preco->CB1_He1 = _h_P_*_c_*_L_He1_ion_/_k_B_;
  preco->CB1_He2 = _h_P_*_c_*_L_He2_ion_/_k_B_;
  preco->CR = 2.*_PI_*(_m_e_/_h_P_)*(_k_B_/_h_P_);
  preco->CK = pow(Lalpha,3)/(8.*_PI_);
  preco->CK_He = pow(Lalpha_He,3)/(8.*_PI_);
  preco->CL = _c_*_h_P_/(_k_B_*Lalpha);
  preco->CL_He = _c_*_h_P_/(_k_B_/_L_He_2s_);
  preco->CT = (8./3.) * (_sigma_/(_m_e_*_c_)) *
    (8.*pow(_PI_,5)*pow(_k_B_,4)/ 15./ pow(_h_P_,3)/pow(_c_,3));

  preco->Bfact = _h_P_*_c_*(_L_He_2p_-_L_He_2s_)/_k_B_;

  /** - define the fields of the 'thermodynamics parameter and workspace' structure */
  tpaw.preco = preco;
  tpaw.pvecback = pvecback;

  /** - impose initial conditions at early times */

  class_test(zinitial < ppr->recfast_z_He_3,
             error_message_,
             "increase zinitial, otherwise should get initial conditions from recfast's get_init routine (less precise anyway)");

  y[0] = 1.;
  y[1] = 1.;
  x0 = 1.+2.*preco->fHe;
  y[2] = preco->Tnow*(1.+z);

  /** - loop over redshift steps Nz; integrate over each step with
      generic_integrator(), store the results in the table using
      thermodynamics_derivs_with_recfast()*/

  for(i=0; i <Nz; i++) {

    zstart = zinitial * (double)(Nz-i) / (double)Nz;
    zend   = zinitial * (double)(Nz-i-1) / (double)Nz;

    z = zend;

    /** - --> first approximation: H and Helium fully ionized */

    if (z > ppr->recfast_z_He_1+ppr->recfast_delta_z_He_1) {
      x_H0 = 1.;
      x_He0 = 1.;
      x0 = 1.+2.*preco->fHe;
      y[0] = x_H0;
      y[1] = x_He0;
      y[2] = preco->Tnow*(1.+z);
    }

    /** - --> second approximation: first Helium recombination (analytic approximation) */

    else if (z > ppr->recfast_z_He_2+ppr->recfast_delta_z_He_2) {
      x_H0 = 1.;
      x_He0 = 1.;

      rhs = exp( 1.5*log(preco->CR*preco->Tnow/(1.+z)) - preco->CB1_He2/(preco->Tnow*(1.+z)) ) / preco->Nnow;

      /* smoothed transition */
      if (z > ppr->recfast_z_He_1-ppr->recfast_delta_z_He_1) {
        x0_previous = 1.+2.*preco->fHe;
        x0_new = 0.5*(sqrt(pow((rhs-1.-preco->fHe),2) + 4.*(1.+2.*preco->fHe)*rhs) - (rhs-1.-preco->fHe));

        /* get s from -1 to 1 */
        s = (ppr->recfast_z_He_1-z)/ppr->recfast_delta_z_He_1;
        /* infer f1(s) = smooth function interpolating from 0 to 1 */
        weight = f1(s);

        x0 = weight*x0_new+(1.-weight)*x0_previous;
      }
      /* transition finished */
      else {
        x0 = 0.5*(sqrt(pow((rhs-1.-preco->fHe),2) + 4.*(1.+2.*preco->fHe)*rhs) - (rhs-1.-preco->fHe));
      }

      y[0] = x_H0;
      y[1] = x_He0;
      y[2] = preco->Tnow*(1.+z);
    }

    /** - --> third approximation: first Helium recombination completed */

    else if (z > ppr->recfast_z_He_3+ppr->recfast_delta_z_He_3) {
      x_H0 = 1.;
      x_He0 = 1.;

      /* smoothed transition */
      if (z > ppr->recfast_z_He_2-ppr->recfast_delta_z_He_2) {
        rhs = exp( 1.5*log(preco->CR*preco->Tnow/(1.+z)) - preco->CB1_He2/(preco->Tnow*(1.+z)) ) / preco->Nnow;
        x0_previous = 0.5*(sqrt(pow((rhs-1.-preco->fHe),2) + 4.*(1.+2.*preco->fHe)*rhs) - (rhs-1.-preco->fHe));
        x0_new = 1. + preco->fHe;
        /* get s from -1 to 1 */
        s = (ppr->recfast_z_He_2-z)/ppr->recfast_delta_z_He_2;
        /* infer f1(s) = smooth function interpolating from 0 to 1 */
        weight = f1(s);

        x0 = weight*x0_new+(1.-weight)*x0_previous;

      }
      /* transition finished */
      else {
        x0 = 1.+preco->fHe;
      }

      y[0] = x_H0;
      y[1] = x_He0;
      y[2] = preco->Tnow*(1.+z);
    }

    /** - --> fourth approximation: second Helium recombination starts (analytic approximation) */

    else if (y[1] > ppr->recfast_x_He0_trigger) {
      x_H0 = 1.;

      rhs = 4.*exp(1.5*log(preco->CR*preco->Tnow/(1.+z)) - preco->CB1_He1/(preco->Tnow*(1.+z)))/preco->Nnow;
      x_He0 = 0.5*(sqrt(pow((rhs-1.),2) + 4.*(1.+preco->fHe)*rhs )- (rhs-1.));

      /* smoothed transition */
      if (z > ppr->recfast_z_He_3-ppr->recfast_delta_z_He_3) {
        x0_previous = 1. + preco->fHe;
        x0_new = x_He0;
        /* get s from -1 to 1 */
        s = (ppr->recfast_z_He_3-z)/ppr->recfast_delta_z_He_3;
        /* infer f1(x) = smooth function interpolating from 0 to 1 */
        weight = f1(s);

        x0 = weight*x0_new+(1.-weight)*x0_previous;
      }
      /* transition finished */
      else {
        x0 = x_He0;
      }

      x_He0 = (x0-1.)/preco->fHe;
      y[0] = x_H0;
      y[1] = x_He0;
      y[2] = preco->Tnow*(1.+z);
    }

    /** - --> fifth approximation: second Helium recombination (full
        evolution for Helium), H recombination starts (analytic
        approximation) */

    else if (y[0] > ppr->recfast_x_H0_trigger) {

      rhs = exp(1.5*log(preco->CR*preco->Tnow/(1.+z)) - preco->CB1/(preco->Tnow*(1.+z)))/preco->Nnow;
      x_H0 = 0.5*(sqrt(pow(rhs,2)+4.*rhs) - rhs);

      class_call(generic_integrator(thermodynamics_derivs_with_recfast,
                                    zstart,
                                    zend,
                                    y,
                                    &tpaw,
                                    ppr->tol_thermo_integration,
                                    ppr->smallest_allowed_variation,
                                    &gi),
                 gi.error_message,
                 error_message_);

      y[0] = x_H0;

      /* smoothed transition */
      if (ppr->recfast_x_He0_trigger - y[1] < ppr->recfast_x_He0_trigger_delta) {
        rhs = 4.*exp(1.5*log(preco->CR*preco->Tnow/(1.+z)) - preco->CB1_He1/(preco->Tnow*(1.+z)))/preco->Nnow;
        x0_previous = 0.5*(sqrt(pow((rhs-1.),2) + 4.*(1.+preco->fHe)*rhs )- (rhs-1.));
        x0_new = y[0] + preco->fHe*y[1];
        /* get s from 0 to 1 */
        s = (ppr->recfast_x_He0_trigger - y[1])/ppr->recfast_x_He0_trigger_delta;
        /* infer f2(x) = smooth function interpolating from 0 to 1 */
        weight = f2(s);

        x0 = weight*x0_new+(1.-weight)*x0_previous;
      }
      /* transition finished */
      else {
        x0 = y[0] + preco->fHe*y[1];
      }

    }

    /** - --> last case: full evolution for H and Helium */

    else {

      /* quantities used for smoothed transition */
      if (ppr->recfast_x_H0_trigger - y[0] < ppr->recfast_x_H0_trigger_delta) {
        rhs = exp(1.5*log(preco->CR*preco->Tnow/(1.+z)) - preco->CB1/(preco->Tnow*(1.+z)))/preco->Nnow;
        x_H0 = 0.5*(sqrt(pow(rhs,2)+4.*rhs) - rhs);
      }

      class_call(generic_integrator(thermodynamics_derivs_with_recfast,
                                    zstart,
                                    zend,
                                    y,
                                    &tpaw,
                                    ppr->tol_thermo_integration,
                                    ppr->smallest_allowed_variation,
                                    &gi),
                 gi.error_message,
                 error_message_);

      /* smoothed transition */
      if (ppr->recfast_x_H0_trigger - y[0] < ppr->recfast_x_H0_trigger_delta) {
        /* get s from 0 to 1 */
        s = (ppr->recfast_x_H0_trigger - y[0])/ppr->recfast_x_H0_trigger_delta;
        /* infer f2(s) = smooth function interpolating from 0 to 1 */
        weight = f2(s);

        x0 = weight*y[0]+(1.-weight)*x_H0 + preco->fHe*y[1];

      }
      /* transition finished */
      else {
        x0 = y[0] + preco->fHe*y[1];
      }
    }

    /** - --> store the results in the table */
    /* results are obtained in order of decreasing z, and stored in order of growing z */

    /* redshift */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_z)=zend;

    /* ionization fraction */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_xe)=x0;

    /* Tb */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_Tb)=y[2];

    /* get dTb/dz=dy[2] */
    class_call(thermodynamics_derivs_with_recfast(zend, y, dy, &tpaw, error_message_),
               error_message_,
               error_message_);

    /* wb = (k_B/mu) Tb  = (k_B/mu) Tb */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_wb)
      = _k_B_ / ( _c_ * _c_ * _m_H_ ) * (1. + (1./_not4_ - 1.) * preco->YHe + x0 * (1.-preco->YHe)) * y[2];

    /* cb2 = (k_B/mu) Tb (1-1/3 dlnTb/dlna) = (k_B/mu) Tb (1+1/3 (1+z) dlnTb/dz) */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_cb2)
      = *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_wb)
      * (1. + (1.+zend) * dy[2] / y[2] / 3.);

    /* dkappa/dtau = a n_e x_e sigma_T = a^{-2} n_e(today) x_e sigma_T (in units of 1/Mpc) */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_dkappadtau)
      = (1.+zend) * (1.+zend) * preco->Nnow * x0 * _sigma_ * _Mpc_over_m_;

    /* fprintf(stdout,"%e %e %e %e %e %e\n", */
    /* 	    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_z), */
    /* 	    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_xe), */
    /* 	    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_Tb), */
    /* 	    (1.+zend) * dy[2], */
    /* 	    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_cb2), */
    /* 	    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_dkappadtau) */
    /* 	    ); */

  }

  /** - cleanup generic integrator with cleanup_generic_integrator() */

  class_call(cleanup_generic_integrator(&gi),
             gi.error_message,
             error_message_);

  return _SUCCESS_;
}

/**
 * Subroutine evaluating the derivative with respect to redshift of
 * thermodynamical quantities (from RECFAST version 1.4).
 *
 * Computes derivatives of the three variables to integrate: \f$ d x_H
 * / dz, d x_{He} / dz, d T_{mat} / dz \f$.
 *
 * This is one of the few functions in the code which are passed to
 * the generic_integrator() routine.  Since generic_integrator()
 * should work with functions passed from various modules, the format
 * of the arguments is a bit special:
 *
 * - fixed parameters and workspaces are passed through a generic
 *   pointer. Here, this pointer contains the precision, background
 *   and recombination structures, plus a background vector, but
 *   generic_integrator() doesn't know its fine structure.
 *
 * - the error management is a bit special: errors are not written as
 *   usual to error_message_, but to a generic error_message
 *   passed in the list of arguments.
 *
 * @param z                        Input: redshift
 * @param y                        Input: vector of variable to integrate
 * @param dy                       Output: its derivative (already allocated)
 * @param parameters_and_workspace Input: pointer to fixed parameters (e.g. indices) and workspace (already allocated)
 * @param error_message            Output: error message
 */

int ThermodynamicsModule::thermodynamics_derivs_with_recfast_member(double z, double* y, double* dy, void* parameters_and_workspace, ErrorMsg error_message) {


  /* define local variables */

  double x,n,n_He,Trad,Tmat,x_H,x_He,Hz,dHdz,epsilon;
  double Rup,Rdown,K,K_He,Rup_He,Rdown_He,He_Boltz;
  double timeTh,timeH;
  double sq_0,sq_1;

  /* new in recfast 1.4: */
  double Rdown_trip,Rup_trip,tauHe_s,pHe_s,Doppler,gamma_2Ps,pb,qb,AHcon;
  double tauHe_t,pHe_t,CL_PSt,gamma_2Pt;
  double CfHe_t=0.;
  int Heflag;

  struct thermodynamics_parameters_and_workspace * ptpaw;
  struct recombination * preco;
  double * pvecback;

  /* used for energy injection from dark matter */
  double C;
  //double C_He;
  double energy_rate;

  double tau;
  double chi_heat;
  double chi_ion_H;
  int last_index_back;

  ptpaw = (struct thermodynamics_parameters_and_workspace*)parameters_and_workspace;
  preco = ptpaw->preco;
  pvecback = ptpaw->pvecback;

  x_H = y[0];
  x_He = y[1];
  x = x_H + preco->fHe * x_He;
  Tmat = y[2];

  n = preco->Nnow * (1.+z) * (1.+z) * (1.+z);
  n_He = preco->fHe * n;
  Trad = preco->Tnow * (1.+z);

  class_call(background_module_->background_tau_of_z(z, &tau),
             background_module_->error_message_,
             error_message);

  class_call(background_module_->background_at_tau(
                               tau,
                               pba->short_info,
                               pba->inter_normal,
                               &last_index_back,
                               pvecback),
             background_module_->error_message_,
             error_message);

  class_call(thermodynamics_energy_injection(preco, z, &energy_rate, error_message),
             error_message,
             error_message);

  /* Hz is H in inverse seconds (while pvecback returns [H0/c] in inverse Mpcs) */
  Hz = pvecback[background_module_->index_bg_H_]*_c_/_Mpc_over_m_;

  Rdown=1.e-19*_a_PPB_*pow((Tmat/1.e4),_b_PPB_)/(1.+_c_PPB_*pow((Tmat/1.e4),_d_PPB_));
  Rup = Rdown * pow((preco->CR*Tmat),1.5)*exp(-preco->CDB/Tmat);

  sq_0 = sqrt(Tmat/_T_0_);
  sq_1 = sqrt(Tmat/_T_1_);
  Rdown_He = _a_VF_/(sq_0 * pow((1.+sq_0),(1.-_b_VF_)) * pow((1. + sq_1),(1. + _b_VF_)));
  Rup_He = 4.*Rdown_He*pow((preco->CR*Tmat),1.5)*exp(-preco->CDB_He/Tmat);
  K = preco->CK/Hz;

  /* following is from recfast 1.5 */

  if (ppr->recfast_Hswitch == _TRUE_ )
    K *= 1.
      + ppr->recfast_AGauss1*exp(-pow((log(1.+z)-ppr->recfast_zGauss1)/ppr->recfast_wGauss1,2))
      + ppr->recfast_AGauss2*exp(-pow((log(1.+z)-ppr->recfast_zGauss2)/ppr->recfast_wGauss2,2));

  /* end of new recfast 1.5 piece */

  /* following is from recfast 1.4 */

  Rdown_trip = _a_trip_/(sq_0*pow((1.+sq_0),(1.-_b_trip_)) * pow((1.+sq_1),(1.+_b_trip_)));
  Rup_trip = Rdown_trip*exp(-_h_P_*_c_*_L_He2St_ion_/(_k_B_*Tmat))*pow(preco->CR*Tmat,1.5)*4./3.;

  if ((x_He < 5.e-9) || (x_He > ppr->recfast_x_He0_trigger2))
    Heflag = 0;
  else
    Heflag = ppr->recfast_Heswitch;

  if (Heflag == 0)
    K_He = preco->CK_He/Hz;
  else {
    tauHe_s = _A2P_s_*preco->CK_He*3.*n_He*(1.-x_He)/Hz;
    pHe_s = (1.-exp(-tauHe_s))/tauHe_s;
    K_He = 1./(_A2P_s_*pHe_s*3.*n_He*(1.-x_He));

    /*    if (((Heflag == 2) || (Heflag >= 5)) && (x_H < 0.99999)) { */
    if (((Heflag == 2) || (Heflag >= 5)) && (x_H < 0.9999999)) { /* threshold changed by Antony Lewis in 2008 to get smoother Helium */

      Doppler = 2.*_k_B_*Tmat/(_m_H_*_not4_*_c_*_c_);
      Doppler = _c_*_L_He_2p_*sqrt(Doppler);
      gamma_2Ps = 3.*_A2P_s_*preco->fHe*(1.-x_He)*_c_*_c_
        /(sqrt(_PI_)*_sigma_He_2Ps_*8.*_PI_*Doppler*(1.-x_H))
        /pow(_c_*_L_He_2p_,2);
      pb = 0.36;
      qb = ppr->recfast_fudge_He;
      AHcon = _A2P_s_/(1.+pb*pow(gamma_2Ps,qb));
      K_He=1./((_A2P_s_*pHe_s+AHcon)*3.*n_He*(1.-x_He));
    }

    if (Heflag >= 3) {
      tauHe_t = _A2P_t_*n_He*(1.-x_He)*3./(8.*_PI_*Hz*pow(_L_He_2Pt_,3));
      pHe_t = (1. - exp(-tauHe_t))/tauHe_t;
      CL_PSt = _h_P_*_c_*(_L_He_2Pt_ - _L_He_2St_)/_k_B_;
      if ((Heflag == 3) || (Heflag == 5) || (x_H >= 0.99999)) {
        CfHe_t = _A2P_t_*pHe_t*exp(-CL_PSt/Tmat);
        CfHe_t = CfHe_t/(Rup_trip+CfHe_t);
      }
      else {
        Doppler = 2.*_k_B_*Tmat/(_m_H_*_not4_*_c_*_c_);
        Doppler = _c_*_L_He_2Pt_*sqrt(Doppler);
        gamma_2Pt = 3.*_A2P_t_*preco->fHe*(1.-x_He)*_c_*_c_
          /(sqrt(_PI_)*_sigma_He_2Pt_*8.*_PI_*Doppler*(1.-x_H))
          /pow(_c_*_L_He_2Pt_,2);
        pb = 0.66;
        qb = 0.9;
        AHcon = _A2P_t_/(1.+pb*pow(gamma_2Pt,qb))/3.;
        CfHe_t = (_A2P_t_*pHe_t+AHcon)*exp(-CL_PSt/Tmat);
        CfHe_t = CfHe_t/(Rup_trip+CfHe_t);
      }
    }
  }

  /* end of new recfast 1.4 piece */

  timeTh=(1./(preco->CT*pow(Trad,4)))*(1.+x+preco->fHe)/x;
  timeH=2./(3.*preco->H0*pow(1.+z,1.5));

  /************/
  /* hydrogen */
  /************/

  if (x_H > ppr->recfast_x_H0_trigger)
    dy[0] = 0.;
  else {

    /* Peebles' coefficient (approximated as one when the Hydrogen
       ionization fraction is very close to one) */
    if (x_H < ppr->recfast_x_H0_trigger2) {
      C = (1. + K*_Lambda_*n*(1.-x_H))/(1./preco->fu+K*_Lambda_*n*(1.-x_H)/preco->fu +K*Rup*n*(1.-x_H));
    }
    else {
      C = 1.;
    }

    /* For DM annihilation: fraction of injected energy going into
       ionization and Lya excitation */

    /* - old approximation from Chen and Kamionkowski: */

    //chi_ion_H = (1.-x)/3.;

    /* coefficient as revised by Slatyer et al. 2013 (in fact it is a fit by Vivian Poulin of columns 1 and 2 in Table V of Slatyer et al. 2013): */

    if (x < 1.)
      chi_ion_H = 0.369202*pow(1.-pow(x,0.463929),1.70237);
    else
      chi_ion_H = 0.;

    /* evolution of hydrogen ionisation fraction: */

    // JL: test for debugginf reio_inter
    //fprintf(stdout,"%e  %e  %e  %e\n",z,Tmat,K*_Lambda_*n,K*Rup*n);

    dy[0] = (x*x_H*n*Rdown - Rup*(1.-x_H)*exp(-preco->CL/Tmat)) * C / (Hz*(1.+z))       /* Peeble's equation with fudged factors */
      -energy_rate*chi_ion_H/n*(1./_L_H_ion_+(1.-C)/_L_H_alpha_)/(_h_P_*_c_*Hz*(1.+z)); /* energy injection (neglect fraction going to helium) */

  }

  /************/
  /* helium   */
  /************/

  if (x_He < 1.e-15)
    dy[1]=0.;
  else {

    if (preco->Bfact/Tmat < 680.)
      He_Boltz=exp(preco->Bfact/Tmat);
    else
      He_Boltz=exp(680.);

    /* equations modified to take into account energy injection from dark matter */
    //C_He=(1. + K_He*_Lambda_He_*n_He*(1.-x_He)*He_Boltz)/(1. + K_He*(_Lambda_He_+Rup_He)*n_He*(1.-x_He)*He_Boltz);

    dy[1] = ((x*x_He*n*Rdown_He - Rup_He*(1.-x_He)*exp(-preco->CL_He/Tmat))
             *(1. + K_He*_Lambda_He_*n_He*(1.-x_He)*He_Boltz))
      /(Hz*(1+z)* (1. + K_He*(_Lambda_He_+Rup_He)*n_He*(1.-x_He)*He_Boltz)); /* in case of energy injection due to DM, we neglect the contribution to helium ionization */

    /* following is from recfast 1.4 */
    /* this correction is not self-consistent when there is energy injection  from dark matter, and leads to nan's  at small redshift (unimportant when reionization takes over before that redshift) */

    if (Heflag >= 3)
      dy[1] = dy[1] +
        (x*x_He*n*Rdown_trip
         - (1.-x_He)*3.*Rup_trip*exp(-_h_P_*_c_*_L_He_2St_/(_k_B_*Tmat)))
        *CfHe_t/(Hz*(1.+z));

    /* end of new recfast 1.4 piece */

  }

  if (timeTh < preco->H_frac*timeH) {
    /*   dy[2]=Tmat/(1.+z); */
    /* v 1.5: like in camb, add here a smoothing term as suggested by Adam Moss */
    dHdz = -pvecback[background_module_->index_bg_H_prime_]/pvecback[background_module_->index_bg_H_]/pba->a_today*_c_/_Mpc_over_m_;
    epsilon = Hz * (1.+x+preco->fHe) / (preco->CT*pow(Trad,3)*x);
    dy[2] = preco->Tnow + epsilon*((1.+preco->fHe)/(1.+preco->fHe+x))*((dy[0]+preco->fHe*dy[1])/x)
      - epsilon* dHdz/Hz + 3.*epsilon/(1.+z);
  }
  else {
    /* equations modified to take into account energy injection from dark matter */

    //chi_heat = (1.+2.*preio->reionization_table[i*preio->re_size+preio->index_re_xe])/3.; // old approximation from Chen and Kamionkowski

    // coefficient as revised by Slatyer et al. 2013 (in fact it is a fit by Vivian Poulin of columns 1 and 2 in Table V of Slatyer et al. 2013)
    if (x < 1.)
      chi_heat = MIN(0.996857*(1.-pow(1.-pow(x,0.300134),1.51035)),1);
    else
      chi_heat = 1.;

    dy[2]= preco->CT * pow(Trad,4) * x / (1.+x+preco->fHe) * (Tmat-Trad) / (Hz*(1.+z)) + 2.*Tmat/(1.+z)
      -2./(3.*_k_B_)*energy_rate*chi_heat/n/(1.+preco->fHe+x)/(Hz*(1.+z)); /* energy injection */
  }

  return _SUCCESS_;
}

/**
 * This routine merges the two tables 'recombination_table' and
 * 'reionization_table' inside the table 'thermodynamics_table', and
 * frees the temporary structures 'recombination' and 'reionization'.
 *
 * @param ppr   Input: pointer to precision structure
 * @param pth   Input/Output: pointer to thermo structure
 * @param preco Input: pointer to filled recombination structure
 * @param preio Input: pointer to reionization structure
 * @return the error status
 */

int ThermodynamicsModule::thermodynamics_merge_reco_and_reio(recombination* preco, reionization* preio) {
  /** Summary: */

  /** - define local variables */

  int i,index_th,index_re;
  double x0;

  /** - first, a little check that the two tables match each other and can be merged */

  if (pth->reio_parametrization != reio_none) {
    class_test(preco->recombination_table[preio->index_reco_when_reio_start*preco->re_size+preco->index_re_z] !=
               preio->reionization_table[(preio->rt_size -1)*preio->re_size+preio->index_re_z],
               error_message_,
               "mismatch which should never happen");
  }

  /** - find number of redshift in full table = number in reco + number in reio - overlap */

  tt_size_ = ppr->recfast_Nz0 + preio->rt_size - preio->index_reco_when_reio_start - 1;

  /** - add  more points to start earlier in presence of interacting DM */

  if(pba->has_idm_dr == _TRUE_) tt_size_ += ppr->thermo_Nz1_idm_dr + ppr->thermo_Nz2_idm_dr - 1;

  /** - allocate arrays in thermo structure */

  class_alloc(z_table_, tt_size_*sizeof(double), error_message_);
  class_alloc(thermodynamics_table_, th_size_*tt_size_*sizeof(double), error_message_);
  class_alloc(d2thermodynamics_dz2_table_, th_size_*tt_size_*sizeof(double), error_message_);

  /** - fill these arrays */

  for (i=0; i < preio->rt_size; i++) {
    z_table_[i] =
      preio->reionization_table[i*preio->re_size+preio->index_re_z];
    thermodynamics_table_[i*th_size_ + index_th_xe_] =
      preio->reionization_table[i*preio->re_size+preio->index_re_xe];
    thermodynamics_table_[i*th_size_ + index_th_dkappa_] =
      preio->reionization_table[i*preio->re_size+preio->index_re_dkappadtau];
    thermodynamics_table_[i*th_size_ + index_th_Tb_] =
      preio->reionization_table[i*preio->re_size+preio->index_re_Tb];
    thermodynamics_table_[i*th_size_ + index_th_wb_] =
      preio->reionization_table[i*preio->re_size+preio->index_re_wb];
    thermodynamics_table_[i*th_size_ + index_th_cb2_] =
      preio->reionization_table[i*preio->re_size+preio->index_re_cb2];
  }
  for (i=0; i < ppr->recfast_Nz0 - preio->index_reco_when_reio_start - 1; i++) {
    index_th=i+preio->rt_size;
    index_re=i+preio->index_reco_when_reio_start+1;
    z_table_[index_th] =
      preco->recombination_table[index_re*preco->re_size + preco->index_re_z];
    thermodynamics_table_[index_th*th_size_ + index_th_xe_] =
      preco->recombination_table[index_re*preco->re_size + preco->index_re_xe];
    thermodynamics_table_[index_th*th_size_ + index_th_dkappa_] =
      preco->recombination_table[index_re*preco->re_size + preco->index_re_dkappadtau];
    thermodynamics_table_[index_th*th_size_ + index_th_Tb_] =
      preco->recombination_table[index_re*preco->re_size + preco->index_re_Tb];
    thermodynamics_table_[index_th*th_size_ + index_th_wb_] =
      preco->recombination_table[index_re*preco->re_size + preco->index_re_wb];
    thermodynamics_table_[index_th*th_size_ + index_th_cb2_] =
      preco->recombination_table[index_re*preco->re_size + preco->index_re_cb2];
  }

  /** - add more points at larger redshift in presence of interacting
        DM. This is necessary because the value of integrated
        quantitites like tau_idm_dr or tau_idr will then be computed
        exactly up to high redshift. With extrapolations in
        thermodynamics_at_z() we could not obtain this. */

    if(pba->has_idm_dr == _TRUE_){

      for (i = 0; i < ppr->thermo_Nz2_idm_dr + ppr->thermo_Nz1_idm_dr - 1; i++){

        /* with an intermediate step Delta z = (thermo_z_initial_idm_dr-recfast_z_initial)/thermo_Nz1_idm_dr/thermo_Nz1_idm_dr */
        if (i < ppr->thermo_Nz2_idm_dr - 1) {
          index_th = i + preio->rt_size + ppr->recfast_Nz0 - preio->index_reco_when_reio_start - 1;
          z_table_[index_th] = ppr->recfast_z_initial + ((double)i + 1.)*(ppr->thermo_z_initial_idm_dr - ppr->recfast_z_initial)/(double)ppr->thermo_Nz1_idm_dr/(double)ppr->thermo_Nz2_idm_dr;
        }
        /* with a large step Delta z  = (thermo_z_initial_idm_dr-recfast_z_initial)/thermo_Nz1_idm_dr */
        else {
          index_th = (i - ppr->thermo_Nz2_idm_dr + 1) + preio->rt_size + ppr->recfast_Nz0 - preio->index_reco_when_reio_start - 1 + ppr->thermo_Nz2_idm_dr - 1;
          z_table_[index_th] = ppr->recfast_z_initial + ((double)(i - ppr->thermo_Nz2_idm_dr + 1) + 1.)*
          (ppr->thermo_z_initial_idm_dr - ppr->recfast_z_initial)/(double)ppr->thermo_Nz1_idm_dr;
        }
        /* same extrapolation formulas as in thermodynamics_at_z() */
        x0 = thermodynamics_table_[(preio->rt_size + ppr->recfast_Nz0 - preio->index_reco_when_reio_start - 2)*th_size_ + index_th_xe_];
        thermodynamics_table_[index_th*th_size_ + index_th_xe_] = x0;
        thermodynamics_table_[index_th*th_size_ + index_th_dkappa_] = (1. + z_table_[index_th])*(1. + z_table_[index_th])*n_e_*x0*_sigma_*_Mpc_over_m_;
        thermodynamics_table_[index_th*th_size_ + index_th_Tb_] = pba->T_cmb*(1. + z_table_[index_th]);
        thermodynamics_table_[index_th*th_size_ + index_th_wb_] = _k_B_/(_c_*_c_*_m_H_)*(1. + (1./_not4_ - 1.)*YHe_ + x0*(1. - YHe_))*pba->T_cmb*(1. + z_table_[index_th]);
        thermodynamics_table_[index_th*th_size_ + index_th_cb2_] = thermodynamics_table_[index_th*th_size_+index_th_wb_]*4./3.;
      }

  }

  /** - free the temporary structures */

  free(preco->recombination_table);

  if (pth->reio_parametrization != reio_none)
    free(preio->reionization_table);

  return _SUCCESS_;
}

/**
 * Subroutine for formatting thermodynamics output
 */

int ThermodynamicsModule::thermodynamics_output_titles(char titles[_MAXTITLESTRINGLENGTH_]) const {

  class_store_columntitle(titles,"z",_TRUE_);
  class_store_columntitle(titles,"conf. time [Mpc]",_TRUE_);
  class_store_columntitle(titles,"x_e",_TRUE_);
  class_store_columntitle(titles,"kappa' [Mpc^-1]",_TRUE_);
  //class_store_columntitle(titles,"kappa''",_TRUE_);
  //class_store_columntitle(titles,"kappa'''",_TRUE_);
  class_store_columntitle(titles,"exp(-kappa)",_TRUE_);
  class_store_columntitle(titles,"g [Mpc^-1]",_TRUE_);
  //class_store_columntitle(titles,"g'",_TRUE_);
  //class_store_columntitle(titles,"g''",_TRUE_);
  class_store_columntitle(titles,"Tb [K]",_TRUE_);
  class_store_columntitle(titles,"w_b",_TRUE_);
  class_store_columntitle(titles,"c_b^2",_TRUE_);
  class_store_columntitle(titles,"tau_d",_TRUE_);
  //class_store_columntitle(titles,"max. rate",_TRUE_);
  class_store_columntitle(titles,"r_d",pth->compute_damping_scale);

  if(pba->has_idm_dr == _TRUE_){
    class_store_columntitle(titles,"dmu_idm_dr",_TRUE_);
    //class_store_columntitle(titles,"ddmu_idm_dr",_TRUE_);
    //class_store_columntitle(titles,"dddmu_idm_dr",_TRUE_);
    class_store_columntitle(titles,"tau_idm_dr",_TRUE_);
    class_store_columntitle(titles,"tau_idr",_TRUE_);
    class_store_columntitle(titles,"g_idm_dr [Mpc^-1]",_TRUE_);
    class_store_columntitle(titles,"c_idm_dr^2",_TRUE_);
    class_store_columntitle(titles,"T_idm_dr",_TRUE_);
    class_store_columntitle(titles,"dmu_idr",_TRUE_);
  }

  return _SUCCESS_;
}

int ThermodynamicsModule::thermodynamics_output_data(int number_of_titles, double *data) const {

  int index_z, storeidx;
  double *dataptr, *pvecthermo;
  double z,tau;

  //  pth->number_of_thermodynamics_titles = get_number_of_titles(pth->thermodynamics_titles);
  //pth->size_thermodynamics_data = pth->number_of_thermodynamics_titles*tt_size_;


  /* Store quantities: */
  for (index_z = 0; index_z < tt_size_; index_z++){
    dataptr = data + index_z*number_of_titles;
    pvecthermo = thermodynamics_table_ + index_z*th_size_;
    z = z_table_[index_z];
    storeidx=0;

    class_call(background_module_->background_tau_of_z(z, &tau),
               background_module_->error_message_,
               error_message_);

    class_store_double(dataptr, z,_TRUE_, storeidx);
    class_store_double(dataptr, tau,_TRUE_, storeidx);
    class_store_double(dataptr, pvecthermo[index_th_xe_],_TRUE_, storeidx);
    class_store_double(dataptr, pvecthermo[index_th_dkappa_],_TRUE_, storeidx);
    //class_store_double(dataptr, pvecthermo[index_th_ddkappa_],_TRUE_, storeidx);
    //class_store_double(dataptr, pvecthermo[index_th_dddkappa_],_TRUE_, storeidx);
    class_store_double(dataptr, pvecthermo[index_th_exp_m_kappa_],_TRUE_, storeidx);
    class_store_double(dataptr, pvecthermo[index_th_g_],_TRUE_, storeidx);
    //class_store_double(dataptr, pvecthermo[index_th_dg_],_TRUE_, storeidx);
    //class_store_double(dataptr, pvecthermo[index_th_ddg_],_TRUE_, storeidx);
    class_store_double(dataptr, pvecthermo[index_th_Tb_],_TRUE_, storeidx);
    class_store_double(dataptr, pvecthermo[index_th_wb_],_TRUE_, storeidx);
    class_store_double(dataptr, pvecthermo[index_th_cb2_],_TRUE_, storeidx);
    class_store_double(dataptr, pvecthermo[index_th_tau_d_],_TRUE_, storeidx);
    //class_store_double(dataptr, pvecthermo[index_th_rate_],_TRUE_, storeidx);
    class_store_double(dataptr, pvecthermo[index_th_r_d_],pth->compute_damping_scale, storeidx);

    if(pba->has_idm_dr == _TRUE_){
      class_store_double(dataptr, pvecthermo[index_th_dmu_idm_dr_],_TRUE_, storeidx);
      //class_store_double(dataptr, pvecthermo[index_th_ddmu_idm_dr_],_TRUE_, storeidx);
      //class_store_double(dataptr, pvecthermo[index_th_dddmu_idm_dr_],_TRUE_, storeidx);
      class_store_double(dataptr, pvecthermo[index_th_tau_idm_dr_],_TRUE_, storeidx);
      class_store_double(dataptr, pvecthermo[index_th_tau_idr_],_TRUE_, storeidx);
      class_store_double(dataptr, pvecthermo[index_th_g_idm_dr_],_TRUE_, storeidx);
      class_store_double(dataptr, pvecthermo[index_th_cidm_dr2_],_TRUE_, storeidx);
      class_store_double(dataptr, pvecthermo[index_th_Tidm_dr_],_TRUE_, storeidx);
      class_store_double(dataptr, pvecthermo[index_th_dmu_idr_],_TRUE_, storeidx);
    }
  }

  return _SUCCESS_;
}

int ThermodynamicsModule::thermodynamics_tanh(double x, double center, double before, double after, double width, double* result) {

  *result = before + (after-before)*(tanh((x-center)/width)+1.)/2.;

  return _SUCCESS_;
}
