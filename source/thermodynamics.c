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

#include "thermodynamics.h"

#ifdef HYREC
#include "hyrec.h"
#endif

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

int thermodynamics_at_z(
                        struct background * pba,
                        struct thermo * pth,
                        double z,
                        short inter_mode,
                        int * last_index,
                        double * pvecback,
                        double * pvecthermo
                        ) {

  /** Summary: */

  /** - define local variables */

  double x0;

  /* - the fact that z is in the pre-computed range 0 <= z <= z_initial
     will be checked in the interpolation routines below. Before
     trying to interpolate, allow the routine to deal with the case z
     > z_intial: then, all relevant quantities can be extrapolated
     using simple analytic approximations */

  if (z >= pth->z_table[pth->tt_size-1]) {

    /* ionization fraction assumed to remain constant at large z */
    x0= pth->thermodynamics_table[(pth->tt_size-1)*pth->th_size+pth->index_th_xe];
    pvecthermo[pth->index_th_xe] = x0;

    /* Calculate dkappa/dtau (dkappa/dtau = a n_e x_e sigma_T = a^{-2} n_e(today) x_e sigma_T in units of 1/Mpc) */
    pvecthermo[pth->index_th_dkappa] = (1.+z) * (1.+z) * pth->n_e * x0 * _sigma_ * _Mpc_over_m_;

    /* tau_d scales like (1+z)**2 */
    pvecthermo[pth->index_th_tau_d] = pth->thermodynamics_table[(pth->tt_size-1)*pth->th_size+pth->index_th_tau_d]*pow((1+z)/(1.+pth->z_table[pth->tt_size-1]),2);

    if (pth->compute_damping_scale == _TRUE_) {

      /* r_d scales like (1+z)**-3/2 */
      pvecthermo[pth->index_th_r_d] = pth->thermodynamics_table[(pth->tt_size-1)*pth->th_size+pth->index_th_r_d]*pow((1+z)/(1.+pth->z_table[pth->tt_size-1]),-1.5);

    }

    /* Calculate d2kappa/dtau2 = dz/dtau d/dz[dkappa/dtau] given that [dkappa/dtau] proportional to (1+z)^2 and dz/dtau = -H */
    pvecthermo[pth->index_th_ddkappa] = -pvecback[pba->index_bg_H] * 2. / (1.+z) * pvecthermo[pth->index_th_dkappa];

    /* Calculate d3kappa/dtau3 given that [dkappa/dtau] proportional to (1+z)^2 */
    pvecthermo[pth->index_th_dddkappa] = (pvecback[pba->index_bg_H]*pvecback[pba->index_bg_H]/ (1.+z) - pvecback[pba->index_bg_H_prime]) * 2. / (1.+z) * pvecthermo[pth->index_th_dkappa];

    /* \f$ exp^{-\kappa}, g, g', g'' \f$ can be set to zero: they are
       used only for computing the source functions in the
       perturbation module; but source functions only need to be
       sampled below z_initial (the condition that
       z_start_sources<z_initial is checked in the perturbation
       module) */
    pvecthermo[pth->index_th_exp_m_kappa] = 0.;
    pvecthermo[pth->index_th_g]=0.;
    pvecthermo[pth->index_th_dg]=0.;
    pvecthermo[pth->index_th_ddg]=0.;

    /* Calculate Tb */
    pvecthermo[pth->index_th_Tb] = pba->T_cmb*(1.+z);

    /* Calculate cb2 (cb2 = (k_B/mu) Tb (1-1/3 dlnTb/dlna) = (k_B/mu) Tb (1+1/3 (1+z) dlnTb/dz)) */
    /* note that m_H / mu = 1 + (m_H/m_He-1) Y_p + x_e (1-Y_p) */
    pvecthermo[pth->index_th_cb2] = _k_B_ / ( _c_ * _c_ * _m_H_ ) * (1. + (1./_not4_ - 1.) * pth->YHe + x0 * (1.-pth->YHe)) * pba->T_cmb * (1.+z) * 4. / 3.;

    /* derivatives of baryon sound speed (only computed if some non-minimal tight-coupling schemes is requested) */
    if (pth->compute_cb2_derivatives == _TRUE_) {

      /* since cb2 proportional to (1+z) or 1/a, its derivative wrt conformal time is given by dcb2 = - a H cb2 */
      pvecthermo[pth->index_th_dcb2] = - pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a] * pvecthermo[pth->index_th_cb2];

      /* then its second derivative is given by ddcb2 = - a H' cb2 */
      pvecthermo[pth->index_th_ddcb2] = - pvecback[pba->index_bg_H_prime] * pvecback[pba->index_bg_a] * pvecthermo[pth->index_th_cb2];
    }

    /* in this regime, variation rate = dkappa/dtau */
    pvecthermo[pth->index_th_rate] = pvecthermo[pth->index_th_dkappa];

  }

  /** - interpolate in table with array_interpolate_spline() (normal
      mode) or array_interpolate_spline_growing_closeby() (closeby
      mode) */

  else {

    /* some very specific cases require linear interpolation because of a break in the derivative of the functions */
    if (((pth->reio_parametrization == reio_half_tanh) && (z < 2*pth->z_reio))
        || ((pth->reio_parametrization == reio_inter) && (z < 50.))) {

      class_call(array_interpolate_linear(
                                          pth->z_table,
                                          pth->tt_size,
                                          pth->thermodynamics_table,
                                          pth->th_size,
                                          z,
                                          last_index,
                                          pvecthermo,
                                          pth->th_size,
                                          pth->error_message),
                 pth->error_message,
                 pth->error_message);
    }

    /* in the "normal" case, use spline interpolation */
    else {

      if (inter_mode == pth->inter_normal) {

        class_call(array_interpolate_spline(
                                            pth->z_table,
                                            pth->tt_size,
                                            pth->thermodynamics_table,
                                            pth->d2thermodynamics_dz2_table,
                                            pth->th_size,
                                            z,
                                            last_index,
                                            pvecthermo,
                                            pth->th_size,
                                            pth->error_message),
                   pth->error_message,
                   pth->error_message);
      }

      if (inter_mode == pth->inter_closeby) {

        class_call(array_interpolate_spline_growing_closeby(
                                                            pth->z_table,
                                                            pth->tt_size,
                                                            pth->thermodynamics_table,
                                                            pth->d2thermodynamics_dz2_table,
                                                            pth->th_size,
                                                            z,
                                                            last_index,
                                                            pvecthermo,
                                                            pth->th_size,
                                                            pth->error_message),
                   pth->error_message,
                   pth->error_message);

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
int thermodynamics_init(
                        struct precision * ppr,
                        struct background * pba,
                        struct thermo * pth
                        ) {

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
  /* temporary table of values of tau associated with z values in pth->z_table */
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

  double tau;
  double g_max;
  int index_tau_max;

  /** - initialize pointers, allocate background vector */

  preco=&reco;
  preio=&reio;
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);

  if (pth->thermodynamics_verbose > 0)
    printf("Computing thermodynamics");

  /** - compute and check primordial Helium fraction  */

  /* Y_He */
  if (pth->YHe == _BBN_) {
    class_call(thermodynamics_helium_from_bbn(ppr,pba,pth),
               pth->error_message,
               pth->error_message);
    if (pth->thermodynamics_verbose > 0)
      printf(" with Y_He=%.4f\n",pth->YHe);
      
  }
  else {
    if (pth->thermodynamics_verbose > 0)
      printf("\n");
  }

  class_test((pth->YHe < _YHE_SMALL_)||(pth->YHe > _YHE_BIG_),
             pth->error_message,
             "Y_He=%g out of bounds (%g<Y_He<%g)",pth->YHe,_YHE_SMALL_,_YHE_BIG_);

  /** - check energy injection parameters */

  class_test((pth->annihilation<0),
             pth->error_message,
             "annihilation parameter cannot be negative");

  class_test((pth->annihilation>1.e-4),
             pth->error_message,
             "annihilation parameter suspiciously large (%e, while typical bounds are in the range of 1e-7 to 1e-6)",
             pth->annihilation);

  class_test((pth->annihilation_variation>0),
             pth->error_message,
             "annihilation variation parameter must be negative (decreasing annihilation rate)");

  class_test((pth->annihilation_z<0),
             pth->error_message,
             "characteristic annihilation redshift cannot be negative");

  class_test((pth->annihilation_zmin<0),
             pth->error_message,
             "characteristic annihilation redshift cannot be negative");

  class_test((pth->annihilation_zmax<0),
             pth->error_message,
             "characteristic annihilation redshift cannot be negative");

  class_test((pth->annihilation>0)&&(pba->has_cdm==_FALSE_),
             pth->error_message,
             "CDM annihilation effects require the presence of CDM!");

  class_test((pth->annihilation_f_halo>0) && (pth->recombination==recfast),
             pth->error_message,
             "Switching on DM annihilation in halos requires using HyRec instead of RECFAST. Otherwise some values go beyond their range of validity in the RECFAST fits, and the thermodynamics module fails. Two solutions: add 'recombination = HyRec' to your input, or set 'annihilation_f_halo = 0.' (default).");

  class_test((pth->annihilation_f_halo<0),
             pth->error_message,
             "Parameter for DM annihilation in halos cannot be negative");

  class_test((pth->annihilation_z_halo<0),
             pth->error_message,
             "Parameter for DM annihilation in halos cannot be negative");

  if (pth->thermodynamics_verbose > 0)
    if ((pth->annihilation >0) && (pth->reio_parametrization == reio_none) && (ppr->recfast_Heswitch >= 3) && (pth->recombination==recfast))
      printf("Warning: if you have DM annihilation and you use recfast with option recfast_Heswitch >= 3, then the expression for CfHe_t and dy[1] becomes undefined at late times, producing nan's. This is however masked by reionization if you are not in reio_none mode.");

  class_test((pth->decay<0),
             pth->error_message,
             "decay parameter cannot be negative");

  class_test((pth->decay>0)&&(pba->has_cdm==_FALSE_),
             pth->error_message,
             "CDM decay effects require the presence of CDM!");

  /* tests in order to prevent segmentation fault in the following */
  class_test(_not4_ == 0.,
             pth->error_message,
             "stop to avoid division by zero");
  class_test(pth->YHe == 1.,
             pth->error_message,
             "stop to avoid division by zero");

  /** - assign values to all indices in the structures with thermodynamics_indices()*/

  class_call(thermodynamics_indices(pth,preco,preio),
             pth->error_message,
             pth->error_message);

  /** - solve recombination and reionization (if needed) and store values of \f$ z, x_e, d \kappa / d \tau, T_b, c_b^2 \f$ with thermodynamics_solve() */

  class_call(thermodynamics_solve(ppr,pba,pth,preco,preio,pvecback),
             pth->error_message,
             pth->error_message);


  /** - merge tables in recombination and reionization structures into
      a single table in thermo structure */

  class_call(thermodynamics_merge_reco_and_reio(ppr,pth,preco,preio),
             pth->error_message,
             pth->error_message);

  /** - compute table of corresponding conformal times */

  class_alloc(tau_table,pth->tt_size*sizeof(double),pth->error_message);

  for (index_tau=0; index_tau < pth->tt_size; index_tau++) {
    class_call(background_tau_of_z(pba,
                                   pth->z_table[index_tau],
                                   tau_table+index_tau),
               pba->error_message,
               pth->error_message);
  }

  /** - store initial value of conformal time in the structure */

  pth->tau_ini = tau_table[pth->tt_size-1];

  /** - fill missing columns (quantities not computed previously but related) */

  /** - --> baryon drag interaction rate time minus one, -[1/R * kappa'], with R = 3 rho_b / 4 rho_gamma, stored temporarily in column ddkappa */

  last_index_back = pba->bg_size-1;

  for (index_tau=0; index_tau < pth->tt_size; index_tau++) {

    class_call(background_at_tau(pba,
                                 tau_table[index_tau],
                                 pba->normal_info,
                                 pba->inter_closeby,
                                 &last_index_back,
                                 pvecback),
               pba->error_message,
               pth->error_message);

    R = 3./4.*pvecback[pba->index_bg_rho_b]/pvecback[pba->index_bg_rho_g];

    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_ddkappa] =
      -1./R*pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa];

  }

  /** - --> second derivative of this rate, -[1/R * kappa']'', stored temporarily in column dddkappa */
  class_call(array_spline_table_line_to_line(tau_table,
                                             pth->tt_size,
                                             pth->thermodynamics_table,
                                             pth->th_size,
                                             pth->index_th_ddkappa,
                                             pth->index_th_dddkappa,
                                             _SPLINE_EST_DERIV_,
                                             pth->error_message),
             pth->error_message,
             pth->error_message);

  /** - --> compute tau_d = [int_{tau_today}^{tau} dtau -dkappa_d/dtau] */
  class_call(array_integrate_spline_table_line_to_line(tau_table,
                                                       pth->tt_size,
                                                       pth->thermodynamics_table,
                                                       pth->th_size,
                                                       pth->index_th_ddkappa,
                                                       pth->index_th_dddkappa,
                                                       pth->index_th_tau_d,
                                                       pth->error_message),
             pth->error_message,
             pth->error_message);

  /* the temporary quantities stored in columns ddkappa and dddkappa
     will not be used anymore, they will be overwritten */

  /** - --> compute r_d = 2pi/k_d = 2pi * [int_{tau_ini}^{tau} dtau (1/kappa') (R^2+4/5(1+R))/(1+R^2)/6 ]^1/2 (see e.g. Wayne Hu's thesis eq. (5.59) */

  if (pth->compute_damping_scale == _TRUE_) {

    class_alloc(tau_table_growing,pth->tt_size*sizeof(double),pth->error_message);

    /* compute integrand 1/kappa' (R^2+4/5(1+R))/(1+R^2)/6 and store temporarily in column "ddkappa" */
    for (index_tau=0; index_tau < pth->tt_size; index_tau++) {

      tau_table_growing[index_tau]=tau_table[pth->tt_size-1-index_tau];

      class_call(background_at_tau(pba,
                                 tau_table_growing[index_tau],
                                 pba->normal_info,
                                 pba->inter_closeby,
                                 &last_index_back,
                                 pvecback),
               pba->error_message,
               pth->error_message);

      R = 3./4.*pvecback[pba->index_bg_rho_b]/pvecback[pba->index_bg_rho_g];

      pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_ddkappa] =
                 1./pth->thermodynamics_table[(pth->tt_size-1-index_tau)*pth->th_size+pth->index_th_dkappa]
                 *(R*R+4./5.*(1.+R))/(1.+R*R)/6.;

    }

    /* compute second derivative of integrand 1/kappa' and store temporarily in column "dddkappa" */
    class_call(array_spline_table_line_to_line(tau_table_growing,
                                               pth->tt_size,
                                               pth->thermodynamics_table,
                                               pth->th_size,
                                               pth->index_th_ddkappa,
                                               pth->index_th_dddkappa,
                                               _SPLINE_EST_DERIV_,
                                               pth->error_message),
               pth->error_message,
               pth->error_message);

    /* compute integrated quantity r_d^2 and store temporarily in column "g" */
     class_call(array_integrate_spline_table_line_to_line(tau_table_growing,
                                                          pth->tt_size,
                                                          pth->thermodynamics_table,
                                                          pth->th_size,
                                                          pth->index_th_ddkappa,
                                                          pth->index_th_dddkappa,
                                                          pth->index_th_g,
                                                          pth->error_message),
                pth->error_message,
                pth->error_message);

     free(tau_table_growing);

     /* an analytic calculation shows that in the early
        radiation-dominated and ionized universe, when kappa' is
        proportional to (1+z)^2 and tau is proportional to the scale
        factor, the integral is equal to eta/(3 kappa')*2./15. So
        [tau_ini/3/kappa'_ini*2./15.] should be added to the integral in
        order to account for the integration between 0 and tau_ini */

     /* compute r_d */
     for (index_tau=0; index_tau < pth->tt_size; index_tau++) {

       pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_r_d] =
         2.*_PI_*sqrt(tau_table[pth->tt_size-1]/3.
                      /pth->thermodynamics_table[(pth->tt_size-1)*pth->th_size+pth->index_th_dkappa]*2./15.
                      +pth->thermodynamics_table[(pth->tt_size-1-index_tau)*pth->th_size+pth->index_th_g]);

     }

  }

  /** - --> second derivative with respect to tau of dkappa (in view of spline interpolation) */
  class_call(array_spline_table_line_to_line(tau_table,
                                             pth->tt_size,
                                             pth->thermodynamics_table,
                                             pth->th_size,
                                             pth->index_th_dkappa,
                                             pth->index_th_dddkappa,
                                             _SPLINE_EST_DERIV_,
                                             pth->error_message),
             pth->error_message,
             pth->error_message);

  /** - --> first derivative with respect to tau of dkappa (using spline interpolation) */
  class_call(array_derive_spline_table_line_to_line(tau_table,
                                                    pth->tt_size,
                                                    pth->thermodynamics_table,
                                                    pth->th_size,
                                                    pth->index_th_dkappa,
                                                    pth->index_th_dddkappa,
                                                    pth->index_th_ddkappa,
                                                    pth->error_message),
             pth->error_message,
             pth->error_message);

  /** - --> compute -kappa = [int_{tau_today}^{tau} dtau dkappa/dtau], store temporarily in column "g" */
  class_call(array_integrate_spline_table_line_to_line(tau_table,
                                                       pth->tt_size,
                                                       pth->thermodynamics_table,
                                                       pth->th_size,
                                                       pth->index_th_dkappa,
                                                       pth->index_th_dddkappa,
                                                       pth->index_th_g,
                                                       pth->error_message),
             pth->error_message,
             pth->error_message);

  /** - --> derivatives of baryon sound speed (only computed if some non-minimal tight-coupling schemes is requested) */
  if (pth->compute_cb2_derivatives == _TRUE_) {

    /** - ---> second derivative with respect to tau of cb2 */
    class_call(array_spline_table_line_to_line(tau_table,
                                               pth->tt_size,
                                               pth->thermodynamics_table,
                                               pth->th_size,
                                               pth->index_th_cb2,
                                               pth->index_th_ddcb2,
                                               _SPLINE_EST_DERIV_,
                                               pth->error_message),
               pth->error_message,
               pth->error_message);


    /** - ---> first derivative with respect to tau of cb2 (using spline interpolation) */
    class_call(array_derive_spline_table_line_to_line(tau_table,
                                                      pth->tt_size,
                                                      pth->thermodynamics_table,
                                                      pth->th_size,
                                                      pth->index_th_cb2,
                                                      pth->index_th_ddcb2,
                                                      pth->index_th_dcb2,
                                                      pth->error_message),
               pth->error_message,
               pth->error_message);
  }

  free(tau_table);

  /** - --> compute visibility: \f$ g= (d \kappa/d \tau) e^{- \kappa} \f$ */

  /* loop on z (decreasing z, increasing time) */
  for (index_tau=pth->tt_size-1; index_tau>=0; index_tau--) {

    /** - ---> compute g */
    g = pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa] *
      exp(pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g]);

    /** - ---> compute exp(-kappa) */
    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_exp_m_kappa] =
      exp(pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g]);

    /** - ---> compute g' (the plus sign of the second term is correct, see def of -kappa in thermodynamics module!) */
    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dg] =
      (pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_ddkappa] +
       pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa] *
       pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa]) *
      exp(pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g]);

    /** - ---> compute g''  */
    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_ddg] =
      (pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dddkappa] +
       pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa] *
       pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_ddkappa] * 3. +
       pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa] *
       pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa] *
       pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa]) *
      exp(pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g]);

    /** - ---> store g */
    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g] = g;

    /** - ---> compute variation rate */
    class_test(pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa] == 0.,
               pth->error_message,
               "variation rate diverges");

    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_rate] =
      sqrt(pow(pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa],2)
           +pow(pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_ddkappa]/
                pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa],2)
           +fabs(pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dddkappa]/
                 pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa]));

  }

  /** - smooth the rate (details of smoothing unimportant: only the
      order of magnitude of the rate matters) */
  class_call(array_smooth(pth->thermodynamics_table,
                          pth->th_size,
                          pth->tt_size,
                          pth->index_th_rate,
                          ppr->thermo_rate_smoothing_radius,
                          pth->error_message),
             pth->error_message,
             pth->error_message);

  /** - fill tables of second derivatives with respect to z (in view of spline interpolation) */

  class_call(array_spline_table_lines(pth->z_table,
                                      pth->tt_size,
                                      pth->thermodynamics_table,
                                      pth->th_size,
                                      pth->d2thermodynamics_dz2_table,
                                      _SPLINE_EST_DERIV_,
                                      pth->error_message),
             pth->error_message,
             pth->error_message);

  /** - find maximum of g */

  index_tau=pth->tt_size-1;
  while (pth->z_table[index_tau]>_Z_REC_MAX_) {
    index_tau--;
  }

  class_test(pth->thermodynamics_table[(index_tau+1)*pth->th_size+pth->index_th_g] >
             pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g],
             pth->error_message,
             "found a recombination redshift greater or equal to the maximum value imposed in thermodynamics.h, z_rec_max=%g",_Z_REC_MAX_);

  while (pth->thermodynamics_table[(index_tau+1)*pth->th_size+pth->index_th_g] <
         pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g]) {
    index_tau--;
  }

  g_max = pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g];
  index_tau_max = index_tau;

  /* approximation for maximum of g, using cubic interpolation, assuming equally spaced z's */
  pth->z_rec=pth->z_table[index_tau+1]+0.5*(pth->z_table[index_tau+1]-pth->z_table[index_tau])*(pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_g]-pth->thermodynamics_table[(index_tau+2)*pth->th_size+pth->index_th_g])/(pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_g]-2.*pth->thermodynamics_table[(index_tau+1)*pth->th_size+pth->index_th_g]+pth->thermodynamics_table[(index_tau+2)*pth->th_size+pth->index_th_g]);

  class_test(pth->z_rec+ppr->smallest_allowed_variation >= _Z_REC_MAX_,
             pth->error_message,
             "found a recombination redshift greater or equal to the maximum value imposed in thermodynamics.h, z_rec_max=%g",_Z_REC_MAX_);

  class_test(pth->z_rec-ppr->smallest_allowed_variation <= _Z_REC_MIN_,
             pth->error_message,
             "found a recombination redshift smaller or equal to the maximum value imposed in thermodynamics.h, z_rec_min=%g",_Z_REC_MIN_);

  /** - find conformal recombination time using background_tau_of_z() **/

  class_call(background_tau_of_z(pba,pth->z_rec,&(pth->tau_rec)),
             pba->error_message,
             pth->error_message);

  class_call(background_at_tau(pba,pth->tau_rec, pba->long_info, pba->inter_normal, &last_index_back, pvecback),
             pba->error_message,
             pth->error_message);

  pth->rs_rec=pvecback[pba->index_bg_rs];
  pth->ds_rec=pth->rs_rec*pba->a_today/(1.+pth->z_rec);
  pth->da_rec=pvecback[pba->index_bg_ang_distance];
  pth->ra_rec=pth->da_rec*(1.+pth->z_rec)/pba->a_today;
  pth->angular_rescaling=pth->ra_rec/(pba->conformal_age-pth->tau_rec);

  /** - find damping scale at recombination (using linear interpolation) */

  if (pth->compute_damping_scale == _TRUE_) {

    pth->rd_rec = (pth->z_table[index_tau+1]-pth->z_rec)/(pth->z_table[index_tau+1]-pth->z_table[index_tau])*pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_r_d]
      +(pth->z_rec-pth->z_table[index_tau])/(pth->z_table[index_tau+1]-pth->z_table[index_tau])*pth->thermodynamics_table[(index_tau+1)*pth->th_size+pth->index_th_r_d];

  }

  /** - find time (always after recombination) at which tau_c/tau
      falls below some threshold, defining tau_free_streaming */

  class_call(background_tau_of_z(pba,pth->z_table[index_tau],&tau),
             pba->error_message,
             pth->error_message);

  while (1./pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_dkappa]/tau
         < ppr->radiation_streaming_trigger_tau_c_over_tau) {

    index_tau--;

    class_call(background_tau_of_z(pba,pth->z_table[index_tau],&tau),
               pba->error_message,
               pth->error_message);

  }

  pth->tau_free_streaming = tau;

  /** - find baryon drag time (when tau_d crosses one, using linear
      interpolation) and sound horizon at that time */

  index_tau=0;
  while ((pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_tau_d] < 1.) && (index_tau < pth->tt_size))
    index_tau++;

  pth->z_d = pth->z_table[index_tau-1]+
    (1.-pth->thermodynamics_table[(index_tau-1)*pth->th_size+pth->index_th_tau_d])
    /(pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_tau_d]-pth->thermodynamics_table[(index_tau-1)*pth->th_size+pth->index_th_tau_d])
    *(pth->z_table[index_tau]-pth->z_table[index_tau-1]);

  class_call(background_tau_of_z(pba,pth->z_d,&(pth->tau_d)),
             pba->error_message,
             pth->error_message);

  class_call(background_at_tau(pba,pth->tau_d, pba->long_info, pba->inter_normal, &last_index_back, pvecback),
             pba->error_message,
             pth->error_message);

  pth->rs_d=pvecback[pba->index_bg_rs];
  pth->ds_d=pth->rs_d*pba->a_today/(1.+pth->z_d);

  /** - find time above which visibility falls below a given fraction of its maximum */

  index_tau=index_tau_max;
  while ((pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_g] >
          g_max * ppr->neglect_CMB_sources_below_visibility)
         && (index_tau > 0))
    index_tau--;

  class_call(background_tau_of_z(pba,pth->z_table[index_tau],&(pth->tau_cut)),
             pba->error_message,
             pth->error_message);

  /** - if verbose flag set to next-to-minimum value, print the main results */

  if (pth->thermodynamics_verbose > 0) {
    printf(" -> recombination at z = %f\n",pth->z_rec);
    printf("    corresponding to conformal time = %f Mpc\n",pth->tau_rec);
    printf("    with comoving sound horizon = %f Mpc\n",pth->rs_rec);
    printf("    angular diameter distance = %f Mpc\n",pth->da_rec);
    printf("    and sound horizon angle 100*theta_s = %f\n",100.*pth->rs_rec/pth->ra_rec);
    if (pth->compute_damping_scale == _TRUE_) {
      printf("    and with comoving photon damping scale = %f Mpc\n",pth->rd_rec);
      printf("    or comoving damping wavenumber k_d = %f 1/Mpc\n",2.*_PI_/pth->rd_rec);
    }
    printf(" -> baryon drag stops at z = %f\n",pth->z_d);
    printf("    corresponding to conformal time = %f Mpc\n",pth->tau_d);
    printf("    with comoving sound horizon rs = %f Mpc\n",pth->rs_d);
    if ((pth->reio_parametrization == reio_camb) || (pth->reio_parametrization == reio_half_tanh)) {
      if (pth->reio_z_or_tau==reio_tau)
        printf(" -> reionization  at z = %f\n",pth->z_reio);
      if (pth->reio_z_or_tau==reio_z)
        printf(" -> reionization with optical depth = %f\n",pth->tau_reio);
      class_call(background_tau_of_z(pba,pth->z_reio,&tau_reio),
                 pba->error_message,
                 pth->error_message);
      printf("    corresponding to conformal time = %f Mpc\n",tau_reio);
    }
    if (pth->reio_parametrization == reio_bins_tanh) {
      printf(" -> binned reionization gives optical depth = %f\n",pth->tau_reio);
    }
    if (pth->reio_parametrization == reio_many_tanh) {
      printf(" -> many-step reionization gives optical depth = %f\n",pth->tau_reio);
    }
    if (pth->reio_parametrization == reio_inter) {
      printf(" -> interpolated reionization history gives optical depth = %f\n",pth->tau_reio);
    }
    if (pth->thermodynamics_verbose > 1) {
      printf(" -> free-streaming approximation can be turned on as soon as tau=%g Mpc\n",
             pth->tau_free_streaming);
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

int thermodynamics_free(
                        struct thermo * pth
                        ) {

  free(pth->z_table);
  free(pth->thermodynamics_table);
  free(pth->d2thermodynamics_dz2_table);

  return _SUCCESS_;
}

/**
 * Assign value to each relevant index in vectors of thermodynamical quantities,
 * as well as in vector containing reionization parameters.
 *
 *
 * @param pth   Input/Output: pointer to thermo structure
 * @param preco Input/Output: pointer to recombination structure
 * @param preio Input/Output: pointer to reionization structure
 * @return the error status
 */

int thermodynamics_indices(
                           struct thermo * pth,
                           struct recombination * preco,
                           struct reionization * preio
                           ) {

  /** Summary: */

  /** - define local variables */

  /* a running index for the vector of thermodynamics quantities */
  int index;

  /** - initialization of all indices and flags in thermo structure */
  index = 0;

  pth->index_th_xe = index;
  index++;
  pth->index_th_dkappa = index;
  index++;
  pth->index_th_tau_d = index;
  index++;
  pth->index_th_ddkappa = index;
  index++;
  pth->index_th_dddkappa = index;
  index++;
  pth->index_th_exp_m_kappa = index;
  index++;
  pth->index_th_g = index;
  index++;
  pth->index_th_dg = index;
  index++;
  pth->index_th_ddg = index;
  index++;
  pth->index_th_Tb = index;
  index++;
  pth->index_th_cb2 = index;
  index++;

  /* derivatives of baryon sound speed (only computed if some non-minimal tight-coupling schemes is requested) */
  if (pth->compute_cb2_derivatives == _TRUE_) {
    pth->index_th_dcb2 = index;
    index++;
    pth->index_th_ddcb2 = index;
    index++;
  }

  pth->index_th_rate = index;
  index++;

  if (pth->compute_damping_scale == _TRUE_) {
    pth->index_th_r_d = index;
    index++;
  }

  /* end of indices */
  pth->th_size = index;

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
    preio->index_reio_xe_before = index;
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
    preio->index_reio_xe_before = index;
    index++;
    
  }

    /* case where x_e(z) must be interpolated */
  if (pth->reio_parametrization == reio_inter) {

    preio->reio_num_z=pth->reio_inter_num;

    preio->index_reio_first_z = index;
    index+= preio->reio_num_z;
    preio->index_reio_first_xe = index;
    index+= preio->reio_num_z;
    preio->index_reio_xe_before = index;
    index++;
    
  }

  preio->reio_num_params = index;

  /* flags for calling the interpolation routine */

  pth->inter_normal=0;
  pth->inter_closeby=1;

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
int thermodynamics_helium_from_bbn(
                                   struct precision * ppr,
                                   struct background * pba,
                                   struct thermo * pth
                                   ) {

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
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);

  /** - 8.6173e-11 converts from Kelvin to MeV. We randomly choose 0.1 MeV to be the temperature of BBN */
  z_bbn = 0.1/(8.6173e-11*pba->T_cmb)-1.0;

  class_call(background_tau_of_z(pba,
                                 z_bbn,
                                 &tau_bbn),
             pba->error_message,
             pth->error_message);

  class_call(background_at_tau(pba,
                               tau_bbn,
                               pba->long_info,
                               pba->inter_normal,
                               &last_index,
                               pvecback),
             pba->error_message,
             pth->error_message);

  Neff_bbn = (pvecback[pba->index_bg_Omega_r]
	      *pvecback[pba->index_bg_rho_crit]
	      -pvecback[pba->index_bg_rho_g])
    /(7./8.*pow(4./11.,4./3.)*pvecback[pba->index_bg_rho_g]);

  free(pvecback);

  //  printf("Neff early = %g, Neff at bbn: %g\n",pba->Neff,Neff_bbn);

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

  class_open(fA,ppr->sBBN_file, "r",pth->error_message);

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
                   pth->error_message,
                   "could not read value of parameters (num_omegab,num_deltaN) in file %s\n",ppr->sBBN_file);

        class_alloc(omegab,num_omegab*sizeof(double),pth->error_message);
        class_alloc(deltaN,num_deltaN*sizeof(double),pth->error_message);
        class_alloc(YHe,num_omegab*num_deltaN*sizeof(double),pth->error_message);
        class_alloc(ddYHe,num_omegab*num_deltaN*sizeof(double),pth->error_message);
        class_alloc(YHe_at_deltaN,num_omegab*sizeof(double),pth->error_message);
        class_alloc(ddYHe_at_deltaN,num_omegab*sizeof(double),pth->error_message);
        array_line=0;

      }
      else {

        /* read (omegab, deltaN, YHe) */
        class_test(sscanf(line,"%lg %lg %lg",
                          &(omegab[array_line%num_omegab]),
                          &(deltaN[array_line/num_omegab]),
                          &(YHe[array_line])
                          ) != 3,
                   pth->error_message,
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
                                      pth->error_message),
             pth->error_message,
             pth->error_message);

  omega_b=pba->Omega0_b*pba->h*pba->h;

  class_test(omega_b < omegab[0],
             pth->error_message,
             "You have asked for an unrealistic small value omega_b = %e. The corresponding value of the primordial helium fraction cannot be found in the interpolation table. If you really want this value, you should fix YHe to a given value rather than to BBN",
             omega_b);

  class_test(omega_b > omegab[num_omegab-1],
             pth->error_message,
             "You have asked for an unrealistic high value omega_b = %e. The corresponding value of the primordial helium fraction cannot be found in the interpolation table. If you really want this value, you should fix YHe to a given value rather than to BBN",
             omega_b);

  class_test(DeltaNeff < deltaN[0],
             pth->error_message,
             "You have asked for an unrealistic small value of Delta N_eff = %e. The corresponding value of the primordial helium fraction cannot be found in the interpolation table. If you really want this value, you should fix YHe to a given value rather than to BBN",
             DeltaNeff);

  class_test(DeltaNeff > deltaN[num_deltaN-1],
             pth->error_message,
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
                                      pth->error_message),
             pth->error_message,
             pth->error_message);

  /** - spline in remaining dimension (along omegab) */
  class_call(array_spline_table_lines(omegab,
                                      num_omegab,
                                      YHe_at_deltaN,
                                      1,
                                      ddYHe_at_deltaN,
                                      _SPLINE_NATURAL_,
                                      pth->error_message),
             pth->error_message,
             pth->error_message);

  /** - interpolate in remaining dimension (along omegab) */
  class_call(array_interpolate_spline(omegab,
                                      num_omegab,
                                      YHe_at_deltaN,
                                      ddYHe_at_deltaN,
                                      1,
                                      omega_b,
                                      &last_index,
                                      &(pth->YHe),
                                      1,
                                      pth->error_message),
             pth->error_message,
             pth->error_message);

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

int thermodynamics_onthespot_energy_injection(
                                              struct precision * ppr,
                                              struct background * pba,
                                              struct recombination * preco,
                                              double z,
                                              double * energy_rate,
                                              ErrorMsg error_message
                                              ) {

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

  rho_cdm_today = pow(pba->H0*_c_/_Mpc_over_m_,2)*3/8./_PI_/_G_*pba->Omega0_cdm*_c_*_c_; /* energy density in J/m^3 */

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

int thermodynamics_energy_injection(
                                    struct precision * ppr,
                                    struct background * pba,
                                    struct recombination * preco,
                                    double z,
                                    double * energy_rate,
                                    ErrorMsg error_message
                                    ) {

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
      factor = _sigma_ * nH0 / pba->H0 * _Mpc_over_m_ / sqrt(pba->Omega0_b+pba->Omega0_cdm);

      /* integral over z'(=zp) with step dz */
      dz=1.;

      /* first point in trapezoidal integral */
      zp = z;
      class_call(thermodynamics_onthespot_energy_injection(ppr,pba,preco,zp,&onthespot,error_message),
                 error_message,
                 error_message);
      first_integrand = factor*pow(1+z,8)/pow(1+zp,7.5)*exp(2./3.*factor*(pow(1+z,1.5)-pow(1+zp,1.5)))*onthespot; // beware: versions before 2.4.3, there were wrong exponents: 6 and 5.5 instead of 8 and 7.5
      result = 0.5*dz*first_integrand;

      /* other points in trapezoidal integral */
      do {

        zp += dz;
        class_call(thermodynamics_onthespot_energy_injection(ppr,pba,preco,zp,&onthespot,error_message),
                   error_message,
                   error_message);
        integrand = factor*pow(1+z,8)/pow(1+zp,7.5)*exp(2./3.*factor*(pow(1+z,1.5)-pow(1+zp,1.5)))*onthespot; // beware: versions before 2.4.3, there were wrong exponents: 6 and 5.5 instead of 8 and 7.5
        result += dz*integrand;

      } while (integrand/first_integrand > 0.02);

      /* uncomment these lines if you also want to compute the on-the-spot for comparison */
      class_call(thermodynamics_onthespot_energy_injection(ppr,pba,preco,z,&onthespot,error_message),
                 error_message,
                 error_message);

    }
    else {
      class_call(thermodynamics_onthespot_energy_injection(ppr,pba,preco,z,&result,error_message),
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
 * @param x    Output: \f$ X_e(z) \f$
 * @param dx    Output: \f$ X_e(z) \f$
 */

int thermodynamics_reionization_function(
                                         double z,
                                         struct thermo * pth,
                                         struct reionization * preio,
                                         double * x,
                                         double * dx
                                         ) {

  /** Summary: */

  /** - define local variables */
  double argument,dargument;
  int i;
  double z_jump;

  int jump;
  double center,before, after,width,one_jump;
  double z_min, z_max;
  

  /** - implementation of ionization function similar to the one in CAMB */

  if ((pth->reio_parametrization == reio_camb) || (pth->reio_parametrization == reio_half_tanh)) {

    /** - --> case z > z_reio_start */

    if (z > preio->reionization_parameters[preio->index_reio_start]) {

      *x = preio->reionization_parameters[preio->index_reio_xe_before];
      *dx = 0.0;

    }

    else {

      /** - --> case z < z_reio_start: hydrogen contribution (tanh of complicated argument) */

      argument = (pow((1.+preio->reionization_parameters[preio->index_reio_redshift]),
                    preio->reionization_parameters[preio->index_reio_exponent])
                   -pow((1.+z),preio->reionization_parameters[preio->index_reio_exponent]))
                 /(preio->reionization_parameters[preio->index_reio_exponent]
          /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
                 *pow((1.+preio->reionization_parameters[preio->index_reio_redshift]),
                    (preio->reionization_parameters[preio->index_reio_exponent]-1.)))
        /preio->reionization_parameters[preio->index_reio_width];
      /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */

      dargument = -pow((1.+z),(preio->reionization_parameters[preio->index_reio_exponent]-1.))
          /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
                 /pow((1.+preio->reionization_parameters[preio->index_reio_redshift]),
                   (preio->reionization_parameters[preio->index_reio_exponent]-1.))
                 /preio->reionization_parameters[preio->index_reio_width];
      
      
      if (pth->reio_parametrization == reio_camb) {
        *x = (preio->reionization_parameters[preio->index_reio_xe_after]
               -preio->reionization_parameters[preio->index_reio_xe_before])
             *(tanh(argument)+1.)/2.
             +preio->reionization_parameters[preio->index_reio_xe_before];
          
        *dx = (preio->reionization_parameters[preio->index_reio_xe_after]
               -preio->reionization_parameters[preio->index_reio_xe_before])
              *(1-tanh(argument)*tanh(argument))/2.*dargument;  
      }
      else {
        *x = (preio->reionization_parameters[preio->index_reio_xe_after]
               -preio->reionization_parameters[preio->index_reio_xe_before])
             *tanh(argument)
             +preio->reionization_parameters[preio->index_reio_xe_before];
          
        *dx = (preio->reionization_parameters[preio->index_reio_xe_after]
                -preio->reionization_parameters[preio->index_reio_xe_before])
              *(1-tanh(argument)*tanh(argument))*dargument;
      }

      /** - --> case z < z_reio_start: helium contribution (tanh of simpler argument) */

      if (pth->reio_parametrization == reio_camb) {
        argument = (preio->reionization_parameters[preio->index_helium_fullreio_redshift] - z)
          /preio->reionization_parameters[preio->index_helium_fullreio_width];
          
        dargument = -1./preio->reionization_parameters[preio->index_helium_fullreio_width];
        /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
        *x += preio->reionization_parameters[preio->index_helium_fullreio_fraction]
              *(tanh(argument)+1.)/2.;
          
        *dx += preio->reionization_parameters[preio->index_helium_fullreio_fraction]
               *(1-tanh(argument)*tanh(argument))/2.*dargument;
          
      }
    }

    return _SUCCESS_;

  }

  /** - implementation of binned ionization function similar to astro-ph/0606552 */

  if (pth->reio_parametrization == reio_bins_tanh) {

    /** - --> case z > z_reio_start */
 
    if (z > preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1]) {
      *x = preio->reionization_parameters[preio->index_reio_xe_before];
      *dx = 0.0;
      
    }

    else if (z < preio->reionization_parameters[preio->index_reio_first_z]) {
      *x = preio->reionization_parameters[preio->index_reio_first_xe];
      *dx = 0.0;
    }

    else {

      i = 0;
      while (preio->reionization_parameters[preio->index_reio_first_z+i+1]<z) i++;

      
      /* fix the final xe to xe_before*/
      preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-1] = preio->reionization_parameters[preio->index_reio_xe_before];

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

      *x = preio->reionization_parameters[preio->index_reio_first_xe+i]
        +0.5*(tanh((z-z_jump)
                   /preio->reionization_parameters[preio->index_reio_step_sharpness])+1.)
        *(preio->reionization_parameters[preio->index_reio_first_xe+i+1]
          -preio->reionization_parameters[preio->index_reio_first_xe+i]);
     
      *dx = 0.5*(1-tanh((z-z_jump)/preio->reionization_parameters[preio->index_reio_step_sharpness])
                  *tanh((z-z_jump)/preio->reionization_parameters[preio->index_reio_step_sharpness]))
            *(preio->reionization_parameters[preio->index_reio_first_xe+i+1]
              -preio->reionization_parameters[preio->index_reio_first_xe+i]) 
            /preio->reionization_parameters[preio->index_reio_step_sharpness];
    }

    return _SUCCESS_;

  }

  /** - implementation of many tanh jumps */

  if (pth->reio_parametrization == reio_many_tanh) {

    /** - --> case z > z_reio_start */

    if (z > preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1]) {
      *x = preio->reionization_parameters[preio->index_reio_xe_before];
      *dx = 0.0;
      
    }

    else if (z > preio->reionization_parameters[preio->index_reio_first_z]) {

      *x = preio->reionization_parameters[preio->index_reio_xe_before];
      *dx = 0.0;
      
      /* fix the final xe to xe_before*/
      preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-1] = preio->reionization_parameters[preio->index_reio_xe_before];

      
      for (jump=1; jump<preio->reio_num_z-1; jump++){

        center = preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1-jump];
        // before and after are meant with respect to growing z, not growing time
        before = preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-1-jump]
          -preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-jump];
        after = 0.;
        width = preio->reionization_parameters[preio->index_reio_step_sharpness];

        one_jump = before + (after-before)*(tanh((z-center)/width)+1.)/2.;

        *x += one_jump;
        *dx += (after-before)*(1-tanh((z-center)/width)*tanh((z-center)/width))/2./width;

      }

    }

    else {
      *x = preio->reionization_parameters[preio->index_reio_first_xe];
      *dx = 0.0;
    }

    return _SUCCESS_;

  }

    /** - implementation of reio_inter */

  if (pth->reio_parametrization == reio_inter) {

    /** - --> case z > z_reio_start */

    if (z > preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1]){
      *x = preio->reionization_parameters[preio->index_reio_xe_before];
      *dx = 0.0;
    }
    else{

      i=0;
      while (preio->reionization_parameters[preio->index_reio_first_z+i+1] < z) i++;

      z_min = preio->reionization_parameters[preio->index_reio_first_z+i];
      z_max = preio->reionization_parameters[preio->index_reio_first_z+i+1];

      /* fix the final xe to xe_before*/
      preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-1] = preio->reionization_parameters[preio->index_reio_xe_before];

      
      class_test(z<z_min,
                 pth->error_message,
                 "z out of range for reionization interpolation");

      class_test(z>z_max,
                 pth->error_message,
                 "z out of range for reionization interpolation");

      argument =(z-preio->reionization_parameters[preio->index_reio_first_z+i])
        /(preio->reionization_parameters[preio->index_reio_first_z+i+1]
          -preio->reionization_parameters[preio->index_reio_first_z+i]);

      dargument = 1./(preio->reionization_parameters[preio->index_reio_first_z+i+1]
          -preio->reionization_parameters[preio->index_reio_first_z+i]);
        
      *x = preio->reionization_parameters[preio->index_reio_first_xe+i]
        + argument*(preio->reionization_parameters[preio->index_reio_first_xe+i+1]
             -preio->reionization_parameters[preio->index_reio_first_xe+i]);

      *dx = dargument*(preio->reionization_parameters[preio->index_reio_first_xe+i+1]
             -preio->reionization_parameters[preio->index_reio_first_xe+i]);  
        
      class_test(*x<0.,
                 pth->error_message,
                 "Interpolation gives negative ionization fraction\n",
                 argument,
                 preio->reionization_parameters[preio->index_reio_first_xe+i],
                 preio->reionization_parameters[preio->index_reio_first_xe+i+1]);

    }

    return _SUCCESS_;

  }

  class_test(0 == 0,
             pth->error_message,
             "value of reio_parametrization=%d unclear",pth->reio_parametrization);
}

/**
 * This subroutine reads \f$ X_e(z) \f$ and \f$ Tmat(z) \f$ in the recombination table at
 * the time at which reionization starts. Hence it provides correct
 * initial conditions for the reionization function.
 *
 * @param ppr   Input: pointer to precision structure
 * @param pth   Input: pointer to thermo structure
 * @param preco Input: pointer to recombination structure
 * @param z     Input: redshift z_reio_start
 * @param x    Output: \f$ X_e(z) \f$ at z
 * @param Tmat  Output: \f$ Tmat(z) \f$ at z
 */

int thermodynamics_interpolate_recombination_table(
                                              struct precision * ppr,
                                              struct thermo * pth,
                                              struct recombination * preco,
                                              double z,
                                              double * x,
                                              double * Tmat
                                              ) {

  int last_index=0;

  class_call(array_interpolate_one_growing_closeby(preco->recombination_table,
                                                   preco->re_size,
                                                   preco->rt_size,
                                                   preco->index_re_z,
                                                   z,
                                                   &last_index,
                                                   preco->index_re_xe,
                                                   x,
                                                   pth->error_message),
             pth->error_message,
             pth->error_message);
  
  last_index=0;
  
  class_call(array_interpolate_one_growing_closeby(preco->recombination_table,
                                                   preco->re_size,
                                                   preco->rt_size,
                                                   preco->index_re_z,
                                                   z,
                                                   &last_index,
                                                   preco->index_re_Tb,
                                                   Tmat,
                                                   pth->error_message),
             pth->error_message,
             pth->error_message);

  return _SUCCESS_;

}

/**
 * Integrate thermodynamics with your favorite recombination code.
 *
 */

int thermodynamics_solve(
                         struct precision * ppr,
                         struct background * pba,
                         struct thermo * pth,
                         struct recombination * preco,
                         struct reionization * preio,
                         double * pvecback
                         ) {

  struct thermo_workspace * ptw;
  
  /** - allocate and initialize the 'thermo workspace' structure */
  class_alloc(ptw,sizeof(struct thermo_workspace),pth->error_message);
  
  class_call(thermo_workspace_init(ppr,
                                   pba,
                                   pth,
                                   ptw),
             pth->error_message,
             pth->error_message);
  
  if (pth->recombination == hyrec) {

    /* Compute recombination */
    class_call(thermodynamics_recombination_with_hyrec(ppr,pba,pth,preco,pvecback),
               pth->error_message,
               pth->error_message);

    /* Compute reionization if needed*/
    if(pth->reio_parametrization != reio_none){    
    class_call(thermodynamics_solve_with_recfast(ppr,pba,pth,preco,preio,ptw,pvecback),
               pth->error_message,
               pth->error_message);
    }
    
  }
  else if (pth->recombination == recfast) {
    
    /* Compute thermodynamics */
    class_call(thermodynamics_solve_with_recfast(ppr,pba,pth,preco,preio,ptw,pvecback),
               pth->error_message,
               pth->error_message);
  }
  
  class_call(thermo_workspace_free(ptw),
             pth->error_message,
             pth->error_message);    
  
 return _SUCCESS_;

}

/**
 * Integrate thermodynamics with HyRec.
 *
 * Integrate thermodynamics with HyRec, allocate and fill the part
 * of the thermodynamics interpolation table (the rest is filled in
 * thermodynamics_init()). Called once by
 * thermodynamics_solve(), from thermodynamics_init().
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

int thermodynamics_recombination_with_hyrec(
                                            struct precision * ppr,
                                            struct background * pba,
                                            struct thermo * pth,
                                            struct recombination * preco,
                                            double * pvecback
                                            ) {
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
  param.omh2 = (pba->Omega0_b+pba->Omega0_cdm+pba->Omega0_ncdm_tot)*pba->h*pba->h;
  param.okh2 = pba->Omega0_k*pba->h*pba->h;
  param.odeh2 = (pba->Omega0_lambda+pba->Omega0_fld)*pba->h*pba->h;
  class_call(background_w_fld(pba,pba->a_today,&w_fld,&dw_over_da_fld,&integral_fld), pba->error_message, pth->error_message);
  param.w0 = w_fld;
  param.wa = -dw_over_da_fld*pba->a_today;
  param.Y = pth->YHe;
  param.Nnueff = pba->Neff;
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
              pth->error_message);

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

  class_open(fA,ppr->hyrec_Alpha_inf_file, "r",pth->error_message);
  class_open(fR,ppr->hyrec_R_inf_file, "r",pth->error_message);

  for (i = 0; i < NTR; i++) {
    for (j = 0; j < NTM; j++) {
      for (l = 0; l <= 1; l++) {
        if (fscanf(fA, "%le", &(rate_table.logAlpha_tab[l][j][i])) != 1)
          class_stop(pth->error_message,"Error reading hyrec data file %s",ppr->hyrec_Alpha_inf_file);
        rate_table.logAlpha_tab[l][j][i] = log(rate_table.logAlpha_tab[l][j][i]);
      }
    }

    if (fscanf(fR, "%le", &(rate_table.logR2p2s_tab[i])) !=1)
      class_stop(pth->error_message,"Error reading hyrec data file %s",ppr->hyrec_R_inf_file);
    rate_table.logR2p2s_tab[i] = log(rate_table.logR2p2s_tab[i]);

  }
  fclose(fA);
  fclose(fR);

  /* Read two-photon rate tables */

  class_open(fA,ppr->hyrec_two_photon_tables_file, "r",pth->error_message);

  for (b = 0; b < NVIRT; b++) {
    if ((fscanf(fA, "%le", &(twog_params.Eb_tab[b])) != 1) ||
        (fscanf(fA, "%le", &(twog_params.A1s_tab[b])) != 1) ||
        (fscanf(fA, "%le", &(twog_params.A2s_tab[b])) != 1) ||
        (fscanf(fA, "%le", &(twog_params.A3s3d_tab[b])) != 1) ||
        (fscanf(fA, "%le", &(twog_params.A4s4d_tab[b])) != 1))
      class_stop(pth->error_message,"Error reading hyrec data file %s",ppr->hyrec_two_photon_tables_file);
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
  preco->YHe = pth->YHe;
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
  pth->n_e=preco->Nnow;

  /** - allocate memory for thermodynamics interpolation tables (size known in advance) and fill it */

  class_alloc(preco->recombination_table,preco->re_size*preco->rt_size*sizeof(double),pth->error_message);

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
                                             pth->error_message),
               pth->error_message,
               pth->error_message);

    class_call(array_interpolate_cubic_equal(-log(1.+param.zstart),
                                             param.dlna,
                                             Tm_output,
                                             param.nz,
                                             -log(1.+z),
                                             &Tm,
                                             pth->error_message),
               pth->error_message,
               pth->error_message);

    class_call(background_tau_of_z(pba,
                                   z,
                                   &tau),
               pba->error_message,
               pth->error_message);

    class_call(background_at_tau(pba,
                                 tau,
                                 pba->short_info,
                                 pba->inter_normal,
                                 &last_index_back,
                                 pvecback),
               pba->error_message,
               pth->error_message);

    /*   class_call(thermodynamics_energy_injection(ppr,pba,preco,z,&energy_rate,pth->error_message),
         pth->error_message,
         pth->error_message);
    */

    /* Hz is H in inverse seconds (while pvecback returns [H0/c] in inverse Mpcs) */
    Hz=pvecback[pba->index_bg_H] * _c_ / _Mpc_over_m_;

    /** - --> store the results in the table */

    /* results are obtained in order of decreasing z, and stored in order of growing z */

    /* redshift */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_z)=z;

    /* ionization fraction */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_xe)=xe;

    /* Tb */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_Tb)=Tm;

    /* cb2 = (k_B/mu) Tb (1-1/3 dlnTb/dlna) = (k_B/mu) Tb (1+1/3 (1+z) dlnTb/dz)
       with (1+z)dlnTb/dz= - [dlnTb/dlna] */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_cb2)
      = _k_B_ / ( _c_ * _c_ * _m_H_ ) * (1. + (1./_not4_ - 1.) * pth->YHe + xe * (1.-pth->YHe)) * Tm * (1. - rec_dTmdlna(xe, Tm, pba->T_cmb*(1.+z), Hz, param.fHe, param.nH0*pow((1+z),3)*1e-6, energy_injection_rate(&param,z)) / Tm / 3.);

    /* dkappa/dtau = a n_e x_e sigma_T = a^{-2} n_e(today) x_e sigma_T (in units of 1/Mpc) */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_dkappadtau)
      = (1.+z) * (1.+z) * preco->Nnow * xe * _sigma_ * _Mpc_over_m_;

  }

  /* Cleanup */

  free(buffer);

#else

  class_stop(pth->error_message,
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

int thermodynamics_solve_with_recfast(
                                      struct precision * ppr,
                                      struct background * pba,
                                      struct thermo * pth,
                                      struct recombination * preco,
                                      struct reionization * preio,
                                      struct thermo_workspace * ptw,
                                      double * pvecback
                                     ) {

  /** Summary: */

  /** - define local variables */

  /* number of time intervals of one approximation scheme*/
  int interval_number;

  /* index running over such time intervals */
  int index_interval;

  /* edge of intervals where approximation scheme is uniform: z_ini, z_switch_1, ..., z_end */
  double * interval_limit;

  /* other recfast variables */
  double zinitial;
  int i,Nz;

  /* contains all fixed parameters which should be passed to thermodynamics_derivs_with_recfast */
  struct thermodynamics_parameters_and_workspace tpaw;

  /* function pointer to ODE evolver and names of possible evolvers */
  extern int evolver_rk();
  extern int evolver_ndf15();
  int (*generic_evolver)();
  double * mz_output;


  /** - read a few precision/cosmological parameters */

  class_call(thermodynamics_recombination_set_parameters(ppr,
                                                         pba,
                                                         pth,
                                                         preco),
             pth->error_message,
             pth->error_message);
  
  class_call(thermodynamics_reionization_set_parameters(ppr,
                                                        pba,
                                                        pth,
                                                        preco,
                                                        preio),
             pth->error_message,
             pth->error_message);
  
  
  
  /* z_initial */
  zinitial=ppr->recfast_z_initial;

  /** - allocate memory for thermodynamics interpolation tables (size known in advance) */
  Nz = preco->Nz_reco + preio->Nz_reio;
  
  preio->rt_size = Nz;
  class_alloc(preio->reionization_table,preio->re_size*preio->rt_size*sizeof(double),pth->error_message);
  
  /** - define the fields of the 'thermodynamics parameter and workspace' structure */
  tpaw.pba = pba;
  tpaw.ppr = ppr;
  tpaw.pth = pth;
  tpaw.preco = preco;
  tpaw.preio = preio;
  tpaw.pvecback = pvecback;
  tpaw.ptw = ptw;

 
  /* create mz_output array (inverted array of decreasing negative redshift) and interval_limit array*/
  
  class_alloc(mz_output,Nz*sizeof(double), pth->error_message);
  class_alloc(interval_limit,(ptw->ap_size+1)*sizeof(double),pth->error_message);
  
  
  for(i=0; i <preco->Nz_reco; i++) {
    mz_output[i] = -(zinitial-ppr->reionization_z_start_max) * (double)(preco->Nz_reco-1-i) / (double)(preco->Nz_reco-1) - ppr->reionization_z_start_max;
  }
  for(i=0; i <preio->Nz_reio; i++) {
    mz_output[i+preco->Nz_reco] = -ppr->reionization_z_start_max * (double)(preio->Nz_reio-1-i) / (double)(preio->Nz_reio);
  }
  
  /*set the switching z's for the approximations */

  class_call(thermodynamics_set_approximation_limits(ppr,
                                                     pba,
                                                     pth,
                                                     preco,
                                                     ptw,
                                                     mz_output[0],
                                                     mz_output[Nz-1],
                                                     &interval_number,
                                                     interval_limit),
             pth->error_message,
             pth->error_message);

 
  /** - loop over intervals over which approximation scheme is uniform. For each interval: */
  for (index_interval=0; index_interval<interval_number; index_interval++) {

    /** - --> (a) fix current approximation scheme. */

    ptw->ap_current = index_interval;
    //printf("\n index: %i, interval number: %i \n",index_interval,interval_number);
    //printf(" --- limit: %e  ---- \n", interval_limit[index_interval+1]);
    /** - --> (b) define the vector of quantities to be integrated
        over. If the current interval starts from the initial time
        zinitial, fill the vector with initial conditions for. If 
        it starts from an approximation switching point,
        redistribute correctly the values from the previous to
        the new vector. */

    class_call(thermo_vector_init(ppr,
                                  pba,
                                  pth,
                                  preco,
                                  preio,
                                  interval_limit[index_interval],
                                  ptw),
               pth->error_message,
               pth->error_message);
    
    /** - --> (c) integrate the quantities over the current interval. */
    if(ppr->evolver == rk){
      generic_evolver = evolver_rk;
    }
    else{
      generic_evolver = evolver_ndf15;
    }
    
    /* If we have tau_reio as input the last evolver step has to be done separately in a loop to approximate tau by dichotomy */
    
    if(pth->reio_z_or_tau == reio_tau && index_interval == ptw->index_ap_reio){
        class_call(thermodynamics_reionization_evolve_with_tau(&tpaw,
                                                               interval_limit[index_interval],
                                                               interval_limit[index_interval+1],
                                                               mz_output,
                                                               Nz),
                   pth->error_message,
                   pth->error_message);                   
        
    }
    else{
      class_call(generic_evolver(thermodynamics_derivs_with_recfast,
                                 interval_limit[index_interval],
                                 interval_limit[index_interval+1],
                                 ptw->tv->y,
                                 ptw->tv->used_in_output,
                                 ptw->tv->tv_size,
                                 &tpaw,
                                 ppr->tol_thermo_integration,
                                 ppr->smallest_allowed_variation,
                                 thermodynamics_timescale_with_recfast,  // timescale
                                 1., // stepsize
                                 mz_output, // values of z for output
                                 Nz, // size of previous array
                                 thermodynamics_sources_with_recfast, // function for output
                                 NULL, // print variables
                                 pth->error_message),
                   pth->error_message,
                   pth->error_message);
    }

  }

  /* Compute reionization optical depth */
  if (pth->reio_parametrization != reio_none && pth->reio_z_or_tau == reio_z) {

    class_call(thermodynamics_reionization_get_tau(ppr,
                                                   pba,
                                                   pth,
                                                   preco,
                                                   preio),
               pth->error_message,
               pth->error_message);
    
    pth->tau_reio=preio->reionization_optical_depth;
      
  }
  

  /** - free quantities allocated at the beginning of the routine */
  if(ptw->ap_size != 0){
    class_call(thermo_vector_free(ptw->tv),
               pth->error_message,
               pth->error_message);
  }
  
  free(interval_limit);
  free(preio->reionization_parameters);
  free(mz_output);


  return _SUCCESS_;
}

/**
 * Subroutine evaluating the derivative with respect to redshift of
 * thermodynamical quantities (from RECFAST version 1.4).
 *
 * Automatically recognizes the current approximation interval and computes
 * the needed derivatives for this interval: \f$ d x_H / dz, d x_{He} / dz, \
 * d T_{mat} / dz \f$.
 *
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
 *   usual to pth->error_message, but to a generic error_message
 *   passed in the list of arguments.
 *
 * @param z                        Input: redshift
 * @param y                        Input: vector of variable to integrate
 * @param dy                       Output: its derivative (already allocated)
 * @param parameters_and_workspace Input: pointer to fixed parameters (e.g. indices) and workspace (already allocated)
 * @param error_message            Output: error message
 */

int thermodynamics_derivs_with_recfast(
                                       double mz,
                                       double * y,
                                       double * dy,
                                       void * parameters_and_workspace,
                                       ErrorMsg error_message
                                       ) {


  /* define local variables */

  double z;
  double x,n,n_He,Trad,Tmat,x_H,x_He,dx_H,dx_He,dx,Hz,dHdz,epsilon;
  double Rup,Rdown,K,K_He,Rup_He,Rdown_He,He_Boltz;
  double timeTh,timeH;
  double sq_0,sq_1;
  int index_y;



  /* new in recfast 1.4: */
  double Rdown_trip,Rup_trip,tauHe_s,pHe_s,Doppler,gamma_2Ps,pb,qb,AHcon;
  double tauHe_t,pHe_t,CL_PSt,gamma_2Pt;
  double CfHe_t=0.;
  int Heflag;

  struct thermodynamics_parameters_and_workspace * ptpaw;
  struct precision * ppr;
  struct background * pba;
  struct thermo * pth;
  struct recombination * preco;
  struct reionization * preio;
  double * pvecback;
  struct thermo_workspace * ptw;
  int ap_current;
 
  /* used for energy injection from dark matter */
  double C;
  //double C_He;
  double energy_rate;

  double tau;
  double chi_heat;
  double chi_ion_H;
  int last_index_back;

  z = -mz;

  /** - rename structure fields (just to avoid heavy notations) */
  ptpaw = parameters_and_workspace;
  ppr = ptpaw->ppr;
  pba = ptpaw->pba;
  pth = ptpaw->pth;
  preco = ptpaw->preco;
  preio = ptpaw->preio;
  pvecback = ptpaw->pvecback;
  ptw = ptpaw->ptw;
  ap_current = ptw->ap_current;

  /** - get background/thermo quantities in this point */
  class_call(background_tau_of_z(pba,
                                 z,
                                 &tau),
             pba->error_message,
             error_message);

  class_call(background_at_tau(pba,
                               tau,
                               pba->short_info,
                               pba->inter_normal,
                               &last_index_back,
                               pvecback),
             pba->error_message,
             error_message);

  class_call(thermodynamics_energy_injection(ppr,pba,preco,z,&energy_rate,error_message),
             error_message,
             error_message);


  /* Hz is H in inverse seconds (while pvecback returns [H0/c] in inverse Mpcs) */
  Hz=pvecback[pba->index_bg_H]* _c_ / _Mpc_over_m_;
  
  n = preco->Nnow * (1.+z) * (1.+z) * (1.+z);
  n_He = preco->fHe * n;
  Trad = preco->Tnow * (1.+z);

  /* Hydrogen and helium */
  /* First check for the current approximation scheme. As long as there is 
     no full recombination, x_H, x_He and x are evolved with analytic functions.
  */
  Tmat = ptw->tv->y[ptw->tv->index_Tmat];

  ptw->Tmat = Tmat;
  ptw->dTmat = -ptw->tv->dy[ptw->tv->index_Tmat];

  class_call(thermodynamics_x_analytic(z,
                                       ppr,
                                       pth,
                                       preco,
                                       preio,
                                       ptw,
                                       ap_current),  
             error_message,
             error_message);
  /* 
   * Obtain values for x_H and x_He and x
   * (H can never be required without He being also required)
   * 
   * */
  
  if(ptw->require_He){
    x_He = y[ptw->tv->index_x_He];
    dx_He = -dy[ptw->tv->index_x_He];
    x = ptw->x_H + preco->fHe * x_He;
  }
  else{
    x_He = ptw->x_He;
    dx_He = ptw->dx_He;
    x = ptw->x;
  }
  
  if(ptw->require_H){
    x_H = y[ptw->tv->index_x_H];
    dx_H = -dy[ptw->tv->index_x_H];
    x = x_H + preco->fHe * x_He;
  }
  else{
    x_H = ptw->x_H;
    dx_H = ptw->dx_H;
  }

  /* During reionization, recalculate x */
  if(ap_current == ptw->index_ap_reio){
    ptw->x = x;
    class_call(thermodynamics_x_analytic(z,
                                         ppr,
                                         pth,
                                         preco,
                                         preio,
                                         ptw,
                                         ap_current),  
               error_message,
               error_message);
  
    x = ptw->x;
    dx = ptw->dx + dx_H + preco->fHe * dx_He;
    
  }
  else if(ap_current == ptw->index_ap_reio_hyrec){
    
    /* get x from the recombination */
    class_call(thermodynamics_interpolate_recombination_table(ppr,
                                                              pth,
                                                              preco,
                                                              z,
                                                              &(ptw->x),
                                                              &(ptw->Tmat)),
               pth->error_message,
               pth->error_message);
  
    class_call(thermodynamics_x_analytic(z,
                                         ppr,
                                         pth,
                                         preco,
                                         preio,
                                         ptw,
                                         ptw->index_ap_reio),  
               error_message,
               error_message);
    
    x = ptw->x;
    dx = ptw->dx;
    
  }
  
  
  if(ptw->require_H){
    
    Rdown=1.e-19*_a_PPB_*pow((Tmat/1.e4),_b_PPB_)/(1.+_c_PPB_*pow((Tmat/1.e4),_d_PPB_));
    Rup = Rdown * pow((preco->CR*Tmat),1.5)*exp(-preco->CDB/Tmat);
    K = preco->CK/Hz;
    
    /* following is from recfast 1.5 */
  
    if (ppr->recfast_Hswitch == _TRUE_ ){
      K *= 1.
        + ppr->recfast_AGauss1*exp(-pow((log(1.+z)-ppr->recfast_zGauss1)/ppr->recfast_wGauss1,2))
        + ppr->recfast_AGauss2*exp(-pow((log(1.+z)-ppr->recfast_zGauss2)/ppr->recfast_wGauss2,2));
    }
    /* end of new recfast 1.5 piece */
    
    /************/
    /* hydrogen */
    /************/
    /* Peebles' coefficient (approximated as one when the Hydrogen
     * ionization fraction is very close to one) */
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

    if (x < 1.){
      chi_ion_H = 0.369202*pow(1.-pow(x,0.463929),1.70237);
    }
    else{
      chi_ion_H = 0.;
    }
    
    /* evolution of hydrogen ionisation fraction: */

    // JL: test for debugginf reio_inter
    //fprintf(stdout,"%e  %e  %e  %e\n",z,Tmat,K*_Lambda_*n,K*Rup*n);
    dy[ptw->tv->index_x_H] = (x*x_H*n*Rdown - Rup*(1.-x_H)*exp(-preco->CL/Tmat)) * C / (Hz*(1.+z));      /* Peeble's equation with fudged factors */

    dy[ptw->tv->index_x_H]+= -energy_rate*chi_ion_H/n*(1./_L_H_ion_+(1.-C)/_L_H_alpha_)/(_h_P_*_c_*Hz*(1.+z));  /* energy injection (neglect fraction going to helium) */
    
    /* Store derivatives for temperature integration */
    dx_H = dy[ptw->tv->index_x_H];
    
  }
  if(ptw->require_He){
    
    sq_0 = sqrt(Tmat/_T_0_);
    sq_1 = sqrt(Tmat/_T_1_);
    Rdown_He = _a_VF_/(sq_0 * pow((1.+sq_0),(1.-_b_VF_)) * pow((1. + sq_1),(1. + _b_VF_)));
    Rup_He = 4.*Rdown_He*pow((preco->CR*Tmat),1.5)*exp(-preco->CDB_He/Tmat);
    
    if ((x_He < 5.e-9) || (x_He > ppr->recfast_x_He0_trigger2)){
      Heflag = 0;
    }
    else{
      Heflag = ppr->recfast_Heswitch;
    }
    if (Heflag == 0){
      K_He = preco->CK_He/Hz;
    }
    else {
      tauHe_s = _A2P_s_*preco->CK_He*3.*n_He*(1.-x_He)/Hz;
      pHe_s = (1.-exp(-tauHe_s))/tauHe_s;
      K_He = 1./(_A2P_s_*pHe_s*3.*n_He*(1.-x_He));

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
        Rdown_trip = _a_trip_/(sq_0*pow((1.+sq_0),(1.-_b_trip_)) * pow((1.+sq_1),(1.+_b_trip_)));
        Rup_trip = Rdown_trip*exp(-_h_P_*_c_*_L_He2St_ion_/(_k_B_*Tmat))*pow(preco->CR*Tmat,1.5)*4./3.;

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
    /************/
    /* helium   */
    /************/

    if (x_He < 1.e-15){
      dy[ptw->tv->index_x_He]=0.;
    }
    else {

      if (preco->Bfact/Tmat < 680.){
        He_Boltz=exp(preco->Bfact/Tmat);
      }
      else{
        He_Boltz=exp(680.);
      }

      /* equations modified to take into account energy injection from dark matter */
      //C_He=(1. + K_He*_Lambda_He_*n_He*(1.-x_He)*He_Boltz)/(1. + K_He*(_Lambda_He_+Rup_He)*n_He*(1.-x_He)*He_Boltz);

      dy[ptw->tv->index_x_He] = ((x*x_He*n*Rdown_He - Rup_He*(1.-x_He)*exp(-preco->CL_He/Tmat))
               *(1. + K_He*_Lambda_He_*n_He*(1.-x_He)*He_Boltz))
        /(Hz*(1+z)* (1. + K_He*(_Lambda_He_+Rup_He)*n_He*(1.-x_He)*He_Boltz)); /* in case of energy injection due to DM, we neglect the contribution to helium ionization */

      /* following is from recfast 1.4 (now reordered) */
      /* this correction is not self-consistent when there is energy injection  from dark matter, and leads to nan's  at small redshift (unimportant when reionization takes   over before that redshift) */
     
      if (Heflag >= 3){
        dy[ptw->tv->index_x_He] = dy[ptw->tv->index_x_He] +
            (x*x_He*n*Rdown_trip
           - (1.-x_He)*3.*Rup_trip*exp(-_h_P_*_c_*_L_He_2St_/(_k_B_*Tmat)))
          *CfHe_t/(Hz*(1.+z));
      }
      /* end of new recfast 1.4 piece */
    }
    /* Store derivatives for temperature integration */
    dx_He = dy[ptw->tv->index_x_He];
  }

  /* Calculate dx depending on approximation scheme */
  if(ap_current == ptw->index_ap_H){
    dx = ptw->dx_H + preco->fHe * dx_He;
  }
  else if(ap_current == ptw->index_ap_frec){
    dx = dx_H + preco->fHe * dx_He;
  }
  else if(ap_current == ptw->index_ap_reio){
    dx = dx_H + preco->fHe * dx_He;// + ptw->dx;
  }
  else{
    dx = ptw->dx;  
  }
  
  /* Matter temperature */
  /* Tmat is integrated always */
  
  timeTh=(1./(preco->CT*pow(Trad,4)))*(1.+x+preco->fHe)/x;
  timeH=2./(3.*preco->H0*pow(1.+z,1.5));
  
  if (timeTh < preco->H_frac*timeH) {
    /*   dy[2]=Tmat/(1.+z); */
    /* v 1.5: like in camb, add here a smoothing term as suggested by Adam Moss */
    dHdz=-pvecback[pba->index_bg_H_prime]/pvecback[pba->index_bg_H]/pba->a_today* _c_ / _Mpc_over_m_;
    epsilon = Hz * (1.+x+preco->fHe) / (preco->CT*pow(Trad,3)*x);
    dy[ptw->tv->index_Tmat] = preco->Tnow + epsilon*((1.+preco->fHe)/(1.+preco->fHe+x))*(dx/x)
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
    
    dy[ptw->tv->index_Tmat]= preco->CT * pow(Trad,4) * x / (1.+x+preco->fHe) * (Tmat-Trad) / (Hz*(1.+z)) + 2.*Tmat/(1.+z)
      -2./(3.*_k_B_)*energy_rate*chi_heat/n/(1.+preco->fHe+x)/(Hz*(1.+z)); /* energy injection */

  }

  /* Store ionization fraction to workspace without smoothing for x*/
  ptw->x = x;
  ptw->dx = dx;
  ptw->x_H = x_H;
  ptw->x_He = x_He;
  ptw->dx_H = dx_H;
  ptw->dx_He = dx_He;
  
  
  /* time-invert derivatives */
  for(index_y=0;index_y<ptw->tv->tv_size;index_y++){
    dy[index_y]=-dy[index_y];
  }

  return _SUCCESS_;
}

/**
 * This routine computes analytic values of x_H, x_He, x and redshift 
 * derivatives dx_H, dx_He and dx (in positive redshift direction) 
 * depending on the input approximation scheme current_ap. Because
 * some functions depend on the current Tmat, directly before calling
 * this routine the thermo_workspace should be updated with the current
 * Tmat at this z.
 *
 * @param z   Input: redshift
 * @param preco Input: pointer to recombination structure
 * @param ptw Input/Output: pointer to thermo workspace
 * @param current_ap Input: index of the wished approximation scheme
 * @return the error status
 */


int thermodynamics_x_analytic(
                              double z,
                              struct precision * ppr,
                              struct thermo * pth,
                              struct recombination * preco,
                              struct reionization * preio,
                              struct thermo_workspace * ptw,
                              int current_ap                          
                              ) {
  
  double x_H, x_He, x, rhs,dx_H,dx_He,dx,sqrt_val,drhs,argument;
  
  /** calculate Hydrogen and Helium fraction with analytical approximations specific for each current approximation interval.
  /** - --> first approximation: H and Helium fully ionized */
  if(current_ap == ptw->index_ap_brec){
    
    x_H = 1.;
    x_He = 1.;
    x = 1. + 2.*preco->fHe;
    dx = 0.;
    dx_H=0.;
    dx_He=0.;
    
  }
  /** - --> second approximation: first Helium recombination (analytic approximation) */
  else if (current_ap == ptw->index_ap_He1) {
    
    /* analytic approximations */
    x_H = 1.;
    x_He = 1.;
    rhs = exp( 1.5*log(preco->CR*ptw->Tmat/(1.+z)/(1.+z)) - preco->CB1_He2/ptw->Tmat ) / preco->Nnow;
    sqrt_val = sqrt(pow((rhs-1.-preco->fHe),2) + 4.*(1.+2.*preco->fHe)*rhs);
    drhs = rhs*((preco->CB1_He2*ptw->dTmat/ptw->Tmat/ptw->Tmat)+1.5*(ptw->dTmat/ptw->Tmat-2./(1.+z)) );
    x = 0.5*(sqrt_val - (rhs-1.-preco->fHe));
    dx_H=0.;
    dx_He=0.;
    dx = 0.5*(  ((rhs-1.-preco->fHe) + 2.*(1.+2.*preco->fHe))/sqrt_val   -   1.  )*drhs;
  }
  /** - --> third approximation: first Helium recombination finished */
  else if (current_ap == ptw->index_ap_He1f) {
    
    /* analytic approximations */
    x_H = 1.;
    x_He = 1.;
    x = 1.+preco->fHe;
    dx_H=0.;
    dx_He=0.;
    dx = 0.;
  }
  /** - --> fourth approximation: second Helium recombination starts */
  else if (current_ap == ptw->index_ap_He2) {
    
    /* analytic approximations */
    x_H = 1.;
    
    rhs = 4.*exp(1.5*log(preco->CR*ptw->Tmat/(1.+z)/(1.+z)) - preco->CB1_He1/ptw->Tmat )/preco->Nnow;
    sqrt_val = sqrt(pow((rhs-1.),2) + 4.*(1.+preco->fHe)*rhs );
    drhs = rhs*((preco->CB1_He1*ptw->dTmat/ptw->Tmat/ptw->Tmat)+1.5*(ptw->dTmat/ptw->Tmat-2./(1.+z)) );
    
    x = 0.5*(sqrt_val - (rhs-1.));
    x_He = (x-1.)/preco->fHe;
    dx_H=0.;
    dx = 0.5*(  ((rhs-1.) + 2.*(1.+ preco->fHe))/sqrt_val   -   1.  )*drhs;
    dx_He = (dx/preco->fHe);
  }
  /** - --> fifth approximation: Hydrogen recombination starts */
  else if (current_ap == ptw->index_ap_H) {
    
    /* analytic approximations */
    rhs = exp(1.5*log(preco->CR*ptw->Tmat/(1.+z)/(1.+z)) - preco->CB1/ptw->Tmat)/preco->Nnow;
    
    sqrt_val = sqrt(pow(rhs,2)+4.*rhs);
    drhs = rhs*((preco->CB1*ptw->dTmat/ptw->Tmat/ptw->Tmat)+1.5*(ptw->dTmat/ptw->Tmat-2./(1.+z)) );
    x_H = 0.5*(sqrt_val - rhs);
    dx_H = 0.5*(  (rhs + 2.)/sqrt_val   -   1.  )*drhs;
    
  }
  /** - --> sixth approximation: reionization */
  else if (current_ap == ptw->index_ap_reio) {
    /* analytic approximations */

    preio->reionization_parameters[preio->index_reio_xe_before] = ptw->x;
      
    class_call(thermodynamics_reionization_function(z,pth,preio,&x,&dx),
             pth->error_message,
             pth->error_message);
    
    
  }
  
  
  /** Save x_H, x_He and x and their derivatives into workspace */
  ptw->x_H = x_H;
  ptw->x_He = x_He;
  ptw->x = x;
  ptw->dx_H = dx_H;
  ptw->dx_He = dx_He;
  ptw->dx = dx;

  return _SUCCESS_;
}    

/**
 * Initialize the field '-->tv' of a thermo_workspace structure, which
 * is a thermo_vector structure. This structure contains indices and
 * values of all quantities which need to be integrated with respect
 * to time (and only them: quantities fixed analytically or obeying
 * constraint equations are NOT included in this vector). 
 * 
 * The routine sets and allocates the vector y, dy and used_in_output 
 * with the right size depending on the current approximation scheme 
 * stored in the workspace. Moreover the initial conditions for each
 * approximation scheme are calculated and set correctly.
 * 
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param preco      Input: pointer to the recombination structure
 * @param mz         Input: negative redshift
 * @param ptw        Input/Output: workspace containing in input the approximation scheme, the background/thermodynamics/metric quantities, and possibly the previous vector y; and in output the new vector y.
 * @return the error status
 */

int thermo_vector_init(
                       struct precision * ppr,
                       struct background * pba,
                       struct thermo * pth,
                       struct recombination * preco,
                       struct reionization * preio,
                       double mz,
                       struct thermo_workspace * ptw /* ptw->tv unallocated if ap_current == index_ap_brec, allocated and filled otherwise */
                       ) {
  /* Initialize variables */
  int index_tv;
  struct thermo_vector * ptv;
  
  double z,x,Tmat;
  
  class_alloc(ptv,sizeof(struct thermo_vector),pth->error_message);
  
  /* mz = Minus z is inverted*/
  z = -mz;
  
  /* Start from no component */
  index_tv = 0;
  
  /* Add common indices (Have to be added before) */
  class_define_index(ptv->index_Tmat,_TRUE_,index_tv,1);
  
  /* Add all components that should be evolved */
  if(ptw->ap_current == ptw->index_ap_brec){
    /* Nothing else to add */
  }
  else if(ptw->ap_current == ptw->index_ap_He1){ 
    /* Nothing else to add */
  }
  else if(ptw->ap_current == ptw->index_ap_He1f){ 
    /* Nothing else to add */
  }
  else if(ptw->ap_current == ptw->index_ap_He2){ 
    /* Nothing else to add */
  }
  else if(ptw->ap_current == ptw->index_ap_H){ 
    class_define_index(ptv->index_x_He,_TRUE_,index_tv,1);
  }
  else if(ptw->ap_current == ptw->index_ap_frec){ 
    class_define_index(ptv->index_x_He,_TRUE_,index_tv,1);
    class_define_index(ptv->index_x_H,_TRUE_,index_tv,1); 
  }
  else if(ptw->ap_current == ptw->index_ap_reio){ 
    class_define_index(ptv->index_x_He,_TRUE_,index_tv,1);
    class_define_index(ptv->index_x_H,_TRUE_,index_tv,1); 
  }
  else if(ptw->ap_current == ptw->index_ap_reio){ 
    /* Nothing else to add */
  }
  
  /* We have now obtained the full size */
  ptv->tv_size = index_tv;
  
  /* Allocate all arrays used during the evolution */
  class_calloc(ptv->y,ptv->tv_size,sizeof(double),pth->error_message);
  class_alloc(ptv->dy,ptv->tv_size*sizeof(double),pth->error_message);
  class_alloc(ptv->used_in_output,ptv->tv_size*sizeof(int),pth->error_message);
  
  for (index_tv=0; index_tv<ptv->tv_size; index_tv++){
    ptv->used_in_output[index_tv] = _TRUE_;
  }
  
  /* setting intial conditions for each approximation */

  /* - in the first scheme (brec = before recombination), we need initial condition for the matter temperature given by the photon temperature */
  if(ptw->ap_current == ptw->index_ap_brec){
    /* Store Tmat in workspace for later use */
    ptw->Tmat = preco->Tnow*(1.+z);
    ptw->dTmat = preco->Tnow;
    
    /* Set the new vector and its indices */ 
    ptw->tv = ptv;

    ptw->tv->y[ptw->tv->index_Tmat] = preco->Tnow*(1.+z);
    
    ptw->tv->dy[ptw->tv->index_Tmat] = -preco->Tnow;
    
    ptw->require_H = _FALSE_;
    ptw->require_He = _FALSE_;
  }
  /* - in this scheme we start to evolve Helium and thus need to set its initial condition via the analytic function */
  else if(ptw->ap_current == ptw->index_ap_H){
    /* Store Tmat in workspace for later use */
    ptw->Tmat = ptw->tv->y[ptw->tv->index_Tmat];
    ptw->dTmat = -ptw->tv->dy[ptw->tv->index_Tmat];
    
    /* Obtain initial contents of new vector analytically, especially x_He */
    class_call(thermodynamics_x_analytic(z,
                                         ppr,
                                         pth,
                                         preco,
                                         preio,
                                         ptw,
                                         ptw->ap_current-1),  
               pth->error_message,
               pth->error_message);
    
    /* Set the new vector and its indices */ 
    ptv->y[ptv->index_Tmat] = ptw->tv->y[ptw->tv->index_Tmat];  
    ptv->dy[ptv->index_Tmat] = ptw->tv->dy[ptw->tv->index_Tmat];
    ptv->y[ptv->index_x_He] = ptw->x_He;
    ptv->dy[ptv->index_x_He] = -ptw->dx_He;
    
    /* Free the old vector and its indices */ 
    class_call(thermo_vector_free(ptw->tv),
               pth->error_message,
               pth->error_message);
               
    /* Copy the new vector into the position of the old one*/
    ptw->tv = ptv;
    
    ptw->require_H = _FALSE_;
    ptw->require_He = _TRUE_;
  }
  /* - in the scheme of full recombination (=frec) we evolve all quantities and thus need to set their initial conditions. Tmat and x_He are solely taken from the previous scheme, x_H is set via the analytic function */
  else if(ptw->ap_current == ptw->index_ap_frec){
    /* Store Tmat in workspace for later use */
    ptw->Tmat = ptw->tv->y[ptw->tv->index_Tmat];
    ptw->dTmat = -ptw->tv->dy[ptw->tv->index_Tmat];
    
    /* Obtain initial contents of new vector analytically, especially x_H */
    class_call(thermodynamics_x_analytic(z,
                                         ppr,
                                         pth,
                                         preco,
                                         preio,
                                         ptw,
                                         ptw->ap_current-1),  
               pth->error_message,
               pth->error_message);
    
    /* Set the new vector and its indices */ 
    ptv->y[ptv->index_Tmat] = ptw->tv->y[ptw->tv->index_Tmat];  
    ptv->dy[ptv->index_Tmat] = ptw->tv->dy[ptw->tv->index_Tmat];
    ptv->y[ptv->index_x_H] = ptw->x_H;
    ptv->dy[ptv->index_x_H] = -ptw->dx_H;
    ptv->y[ptv->index_x_He] = ptw->tv->y[ptw->tv->index_x_He];
    ptv->dy[ptv->index_x_He] = ptw->tv->dy[ptw->tv->index_x_He];
    
    /* Free the old vector and its indices */
    class_call(thermo_vector_free(ptw->tv),
               pth->error_message,
               pth->error_message);
    
    /* Copy the new vector into the position of the old one*/
    
    ptw->tv = ptv;
    
    ptw->require_H = _TRUE_;
    ptw->require_He = _TRUE_;
  }
  /* - during reionization we continue to evolve all quantities. Now all three intial conditions are just taken from the previous scheme */
  else if(ptw->ap_current == ptw->index_ap_reio){

    /* Set the new vector and its indices */ 
    ptv->y[ptv->index_Tmat] = ptw->tv->y[ptw->tv->index_Tmat];  
    ptv->dy[ptv->index_Tmat] = ptw->tv->dy[ptw->tv->index_Tmat];
    ptv->y[ptv->index_x_H] = ptw->tv->y[ptw->tv->index_x_H];
    ptv->dy[ptv->index_x_H] = ptw->tv->dy[ptw->tv->index_x_H];
    ptv->y[ptv->index_x_He] = ptw->tv->y[ptw->tv->index_x_He];
    ptv->dy[ptv->index_x_He] = ptw->tv->dy[ptw->tv->index_x_He];
    
    /* Free the old vector and its indices */
    class_call(thermo_vector_free(ptw->tv),
               pth->error_message,
               pth->error_message);
    
    /* Copy the new vector into the position of the old one*/
    
    ptw->tv = ptv;
    
    ptw->require_H = _TRUE_;
    ptw->require_He = _TRUE_;
  }
  else if(ptw->ap_current == ptw->index_ap_reio_hyrec){
    
    class_call(thermodynamics_interpolate_recombination_table(ppr,
                                                              pth,
                                                              preco,
                                                              z,
                                                              &x,
                                                              &Tmat),
               pth->error_message,
               pth->error_message);
    
    /* Store Tmat in workspace for later use */
    ptw->Tmat = Tmat;
        
    /* Set the new vector and its indices */ 
    ptw->tv = ptv;

    ptw->tv->y[ptw->tv->index_Tmat] = Tmat;
    
    /* Approximated dTmat by adiabatic cooling */
    ptw->tv->dy[ptw->tv->index_Tmat] = -preco->Tnow*(1.+z)*2.;
    
    ptw->require_H = _FALSE_;
    ptw->require_He = _FALSE_;
  }
  /* - in all other approximations we only evolve Tmat and set its initial conditions from the previous scheme */
  else{
    /* Store Tmat in workspace for later use */
    ptw->Tmat = ptw->tv->y[ptw->tv->index_Tmat];
    ptw->dTmat = -ptw->tv->dy[ptw->tv->index_Tmat];
      
    /* Set the new vector and its indices */ 
    ptv->y[ptv->index_Tmat] = ptw->tv->y[ptw->tv->index_Tmat];
    ptv->dy[ptv->index_Tmat] = ptw->tv->dy[ptw->tv->index_Tmat];
    
    /* Free the old vector and its indices */ 
    class_call(thermo_vector_free(ptw->tv),
               pth->error_message,
               pth->error_message);
               
    /* Copy the new vector into the position of the old one*/
    ptw->tv = ptv;
    
    ptw->require_H = _FALSE_;
    ptw->require_He = _FALSE_;
  }
  
  return _SUCCESS_;
}

/**
 * Free the thermo_vector structure.
 *
 * @param tv        Input: pointer to thermo_vector structure to be freed
 * @return the error status
 */

int thermo_vector_free(
                        struct thermo_vector * tv
                        ) {

  free(tv->y);
  free(tv->dy);
  free(tv->used_in_output);
  free(tv);

  return _SUCCESS_;
}

/**
 * Initialize a thermo_workspace structure. All fields are allocated
 * here, with the exception of the thermo_vector '-->tv' field, which
 * is allocated separately in thermo_vector_init. We allocate the
 * thermo_workspace structure before calling thermodynamics_solve_with_recfast.
 * 
 * Here the approximation schemes are set with their respective ending
 * redshifts and their smoothing parameters. 
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param ptw        Input/Output: pointer to thermo_workspace structure which fields are allocated or filled here
 * @return the error status
 */

int thermo_workspace_init(
                           struct precision * ppr,
                           struct background * pba,
                           struct thermo * pth,
                           struct thermo_workspace * ptw
                           ) {

  int index_ap;

  /** - count number of approximations, initialize their indices */
  index_ap=0;

  /* With recombination computed by HyRec, the evolver only has to integrate reionization */
  if(pth->recombination == hyrec){
    
    /** approximations have to appear in chronological order here*/
    if (pth->reio_parametrization != reio_none) {
      class_define_index(ptw->index_ap_reio_hyrec,_TRUE_,index_ap,1);
      ptw->ap_size=index_ap;
    }
    else{
      ptw->ap_size=index_ap;
      
      /*Initialize approximation schemes not included in the evolver loop */
      class_define_index(ptw->index_ap_reio_hyrec,_TRUE_,index_ap,1);
    }
    
    class_define_index(ptw->index_ap_brec,_TRUE_,index_ap,1);
    class_define_index(ptw->index_ap_He1,_TRUE_,index_ap,1);
    class_define_index(ptw->index_ap_He1f,_TRUE_,index_ap,1);
    class_define_index(ptw->index_ap_He2,_TRUE_,index_ap,1);
    class_define_index(ptw->index_ap_H,_TRUE_,index_ap,1);
    class_define_index(ptw->index_ap_frec,_TRUE_,index_ap,1);
    class_define_index(ptw->index_ap_reio,_TRUE_,index_ap,1);

    
    ptw->ap_size_loaded=index_ap;

  }
  else{
    
    /** approximations have to appear in chronological order here*/
    class_define_index(ptw->index_ap_brec,_TRUE_,index_ap,1);
    class_define_index(ptw->index_ap_He1,_TRUE_,index_ap,1);
    class_define_index(ptw->index_ap_He1f,_TRUE_,index_ap,1);
    class_define_index(ptw->index_ap_He2,_TRUE_,index_ap,1);
    class_define_index(ptw->index_ap_H,_TRUE_,index_ap,1);
    class_define_index(ptw->index_ap_frec,_TRUE_,index_ap,1);
    
    if (pth->reio_parametrization != reio_none) {
      class_define_index(ptw->index_ap_reio,_TRUE_,index_ap,1);
      ptw->ap_size=index_ap;
    }
    else{
      ptw->ap_size=index_ap;
      
      /*Initialize approximation schemes not included in the evolver loop */
      class_define_index(ptw->index_ap_reio,_TRUE_,index_ap,1);
    }
    
    class_define_index(ptw->index_ap_reio_hyrec,_TRUE_,index_ap,1);

    ptw->ap_size_loaded=index_ap;

  }  
    

    /** store all ending redshifts for each approximation */
    class_alloc(ptw->ap_z_limits,ptw->ap_size_loaded*sizeof(double),pth->error_message);

    ptw->ap_z_limits[ptw->index_ap_brec] = ppr->recfast_z_He_1+ppr->recfast_delta_z_He_1;
    ptw->ap_z_limits[ptw->index_ap_He1] = ppr->recfast_z_He_2+ppr->recfast_delta_z_He_2;
    ptw->ap_z_limits[ptw->index_ap_He1f] = ppr->recfast_z_He_3+ppr->recfast_delta_z_He_3;
    ptw->ap_z_limits[ptw->index_ap_He2] = 2870.;// TODO :: set correctly
    ptw->ap_z_limits[ptw->index_ap_H] = 1600.;// TODO :: set correctly
    ptw->ap_z_limits[ptw->index_ap_frec] = ppr->reionization_z_start_max;
    ptw->ap_z_limits[ptw->index_ap_reio] = 0.0;
    ptw->ap_z_limits[ptw->index_ap_reio_hyrec] = 0.0;

    
    /** store smoothing deltas for transitions at the beginning of each aproximation */
    class_alloc(ptw->ap_z_limits_delta,ptw->ap_size_loaded*sizeof(double),pth->error_message);
    
    ptw->ap_z_limits_delta[ptw->index_ap_brec] = 0.;
    ptw->ap_z_limits_delta[ptw->index_ap_He1] = ppr->recfast_delta_z_He_1;
    ptw->ap_z_limits_delta[ptw->index_ap_He1f] = ppr->recfast_delta_z_He_2;
    ptw->ap_z_limits_delta[ptw->index_ap_He2] = ppr->recfast_delta_z_He_3;
    ptw->ap_z_limits_delta[ptw->index_ap_H] = 50.; // TODO :: set correctly
    ptw->ap_z_limits_delta[ptw->index_ap_frec] = 50.;// TODO :: set correctly
    ptw->ap_z_limits_delta[ptw->index_ap_reio] = 2.0;// TODO :: set correctly
    ptw->ap_z_limits_delta[ptw->index_ap_reio_hyrec] = 2.0;// TODO :: set correctly


  /*fix current approximation scheme abritrarily */
  ptw->ap_current = ptw->index_ap_brec;
  
  return _SUCCESS_;
}

/**
 * Free the thermo_workspace structure (with the exception of the
 * thermo_vector '-->tv' field, which is freed separately in
 * thermo_vector_free).
 *
 * @param ptw        Input: pointer to perturb_workspace structure to be freed
 * @return the error status
 */

int thermo_workspace_free (
                            struct thermo_workspace * ptw
                            ) {

  //free(ptw->pvecback);
  free(ptw->ap_z_limits);
  free(ptw->ap_z_limits_delta);

  free(ptw);

  return _SUCCESS_;
}

int thermodynamics_recombination_set_parameters(
                       struct precision * ppr,
                       struct background * pba,
                       struct thermo * pth,
                       struct recombination * preco
                       ) {

  /* intern recfast variables */

  double OmegaB,mu_H,Lalpha,Lalpha_He,DeltaB,DeltaB_He;

  /** - read a few precision/cosmological parameters */
  
  if(pth->recombination == hyrec){
    preco->Nz_reco = 0;
  }
  else{
    preco->Nz_reco = ppr->recfast_Nz0;
  }
    
  /* preco->H0 is H0 in inverse seconds (while pba->H0 is [H0/c] in inverse Mpcs) */
  preco->H0 = pba->H0 * _c_ / _Mpc_over_m_;

  /* Omega_b */
  OmegaB = pba->Omega0_b;

  /* Yp */
  preco->YHe = pth->YHe;

  /* Tnow */
  preco->Tnow = pba->T_cmb;

  /* H_frac */
  preco->H_frac = ppr->recfast_H_frac;

  /* H fudging */
  class_test((ppr->recfast_Hswitch != _TRUE_) && (ppr->recfast_Hswitch != _FALSE_),
             pth->error_message,
             "RECFAST error: unknown H fudging scheme");
  preco->fu = ppr->recfast_fudge_H;
  if (ppr->recfast_Hswitch == _TRUE_)
    preco->fu += ppr->recfast_delta_fudge_H;

  /* He fudging */
  class_test((ppr->recfast_Heswitch < 0) || (ppr->recfast_Heswitch > 6),
             pth->error_message,
             "RECFAST error: unknown He fudging scheme");

  /* related quantities */
  mu_H = 1./(1.-preco->YHe);
  //mu_T = _not4_ /(_not4_ - (_not4_-1.)*preco->YHe); /* recfast 1.4*/
  preco->fHe = preco->YHe/(_not4_ *(1.-preco->YHe)); /* recfast 1.4 */
  preco->Nnow = 3.*preco->H0*preco->H0*OmegaB/(8.*_PI_*_G_*mu_H*_m_H_);
  pth->n_e = preco->Nnow;

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

 return _SUCCESS_;

}

int thermodynamics_reionization_set_parameters(
                                struct precision * ppr,
                                struct background * pba,
                                struct thermo * pth,
                                struct recombination * preco,
                                struct reionization * preio
                                ) {
    
  int bin;
  int point;
  double xe_input,xe_actual,z_sup;
  
  /** - allocate the vector of parameters defining the function \f$ X_e(z) \f$ */

  class_alloc(preio->reionization_parameters,preio->reio_num_params*sizeof(double),pth->error_message);

  
  class_test(ppr->reionization_sampling <= 0.0,
             pth->error_message,
             "stop to avoid division by zero. Reionization stepsize has to be larger than zero");
  
  preio->Nz_reio = ppr->reionization_z_start_max / ppr->reionization_sampling;
  /** - (a) if reionization implemented like in CAMB */

  if ((pth->reio_parametrization == reio_camb) || (pth->reio_parametrization == reio_half_tanh)) {

    /** - --> set values of these parameters, excepted those depending on the reionization redshift */

    if (pth->reio_parametrization == reio_camb) {
      preio->reionization_parameters[preio->index_reio_xe_after] = 1. + pth->YHe/(_not4_*(1.-pth->YHe));    /* xe_after_reio: H + singly ionized He (note: segmentation fault impossible, checked before that denominator is non-zero) */
    }
    if (pth->reio_parametrization == reio_half_tanh) {
      preio->reionization_parameters[preio->index_reio_xe_after] = 1.; /* xe_after_reio: neglect He ionization */
      //+ 2*pth->YHe/(_not4_*(1.-pth->YHe));    /* xe_after_reio: H + fully ionized He */
    }
    preio->reionization_parameters[preio->index_reio_exponent] = pth->reionization_exponent; /* reio_exponent */
    preio->reionization_parameters[preio->index_reio_width] = pth->reionization_width;    /* reio_width */
    preio->reionization_parameters[preio->index_helium_fullreio_fraction] = pth->YHe/(_not4_*(1.-pth->YHe)); /* helium_fullreio_fraction (note: segmentation fault impossible, checked before that denominator is non-zero) */
    preio->reionization_parameters[preio->index_helium_fullreio_redshift] = pth->helium_fullreio_redshift; /* helium_fullreio_redshift */
    preio->reionization_parameters[preio->index_helium_fullreio_width] = pth->helium_fullreio_width;    /* helium_fullreio_width */

    class_test(preio->reionization_parameters[preio->index_reio_exponent]==0,
               pth->error_message,
               "stop to avoid division by zero");

    class_test(preio->reionization_parameters[preio->index_reio_width]==0,
               pth->error_message,
               "stop to avoid division by zero");

    class_test(preio->reionization_parameters[preio->index_helium_fullreio_width]==0,
               pth->error_message,
               "stop to avoid division by zero");

    /** - --> if reionization redshift given as an input, initialize the remaining values*/

    if (pth->reio_z_or_tau == reio_z) {

      /* reionization redshift */
      preio->reionization_parameters[preio->index_reio_redshift] = pth->z_reio;

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

        preio->reionization_parameters[preio->index_reio_start] = pth->z_reio;
      }

      class_test(preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max,
                 pth->error_message,
                 "starting redshift for reionization > reionization_z_start_max = %e\n",ppr->reionization_z_start_max);

    }

    /** - --> if reionization optical depth given as an input, find reionization redshift by dichotomy and initialize the remaining values */

    if (pth->reio_z_or_tau == reio_tau) {
           z_sup = ppr->reionization_z_start_max-ppr->reionization_start_factor*pth->reionization_width;
      class_test(z_sup < 0.,
                 pth->error_message,
                 "parameters are such that reionization cannot take place before today while starting after z_start_max; need to increase z_start_max");

      /* maximum possible reionization redshift */
      preio->reionization_parameters[preio->index_reio_redshift] = z_sup;
      /* maximum possible starting redshift */
      preio->reionization_parameters[preio->index_reio_start] = ppr->reionization_z_start_max;
      /* infer xe_before_reio */   
      
    }


  }
  /** - (b) if reionization implemented with reio_bins_tanh scheme */

  else if (pth->reio_parametrization == reio_bins_tanh) {

    /* this algorithm requires at least two bin centers (i.e. at least
       4 values in the (z,xe) array, counting the edges). */
    class_test(pth->binned_reio_num<2,
               pth->error_message,
               "current implementation of binned reio requires at least two bin centers");

    /* check that this input can be interpreted by the code */
    for (bin=1; bin<pth->binned_reio_num; bin++) {
      class_test(pth->binned_reio_z[bin-1]>=pth->binned_reio_z[bin],
                 pth->error_message,
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
               pth->error_message,
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
               pth->error_message,
               "final redshift for reionization = %e, you must change the binning or redefine the way in which the code extrapolates below the first value of z_i",preio->reionization_parameters[preio->index_reio_first_z]);
    */
    if (preio->reionization_parameters[preio->index_reio_first_z] < 0) {
      preio->reionization_parameters[preio->index_reio_first_z] = 0.;
    }

    /* infer xe after reio */
    preio->reionization_parameters[preio->index_reio_first_xe] = 1. + pth->YHe/(_not4_*(1.-pth->YHe));    /* xe_after_reio: H + singly ionized He (note: segmentation fault impossible, checked before that denominator is non-zero) */

    /* pass step sharpness parameter */
    preio->reionization_parameters[preio->index_reio_step_sharpness] = pth->binned_reio_step_sharpness;

  }

  /** - (c) if reionization implemented with reio_many_tanh scheme */

  else if (pth->reio_parametrization == reio_many_tanh) {

    /* this algorithm requires at least one jump centers */
    class_test(pth->many_tanh_num<1,
               pth->error_message,
               "current implementation of reio_many_tanh requires at least one jump center");

    /* check that z input can be interpreted by the code */
    for (bin=1; bin<pth->many_tanh_num; bin++) {
      class_test(pth->many_tanh_z[bin-1]>=pth->many_tanh_z[bin],
                 pth->error_message,
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
        xe_actual = 1. + pth->YHe/(_not4_*(1.-pth->YHe));
      }
      //-2 means "after hydrogen + second helium recombination"
      else if ((xe_input<-1.9) && (xe_input>-2.1)) {
        xe_actual = 1. + 2.*pth->YHe/(_not4_*(1.-pth->YHe));
      }
      //other negative number is nonsense
      else {
        class_stop(pth->error_message,
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
               pth->error_message,
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

    /* infer xe after reio */

    preio->reionization_parameters[preio->index_reio_first_xe] = preio->reionization_parameters[preio->index_reio_first_xe+1];

    /* if we want to model only hydrogen reionization and neglect both helium reionization */
    //preio->reionization_parameters[preio->index_reio_first_xe] = 1.;

    /* if we want to model only hydrogen + first helium reionization and neglect second helium reionization */
    //preio->reionization_parameters[preio->index_reio_first_xe] = 1. + pth->YHe/(_not4_*(1.-pth->YHe));

    /* if we want to model hydrogen + two helium reionization */
    //preio->reionization_parameters[preio->index_reio_first_xe] = 1. + 2.*pth->YHe/(_not4_*(1.-pth->YHe));

    /* pass step sharpness parameter */
    class_test(pth->many_tanh_width<=0,
               pth->error_message,
               "many_tanh_width must be strictly positive, you passed %e",
               pth->many_tanh_width);

    preio->reionization_parameters[preio->index_reio_step_sharpness] = pth->many_tanh_width;

  }

  /** - (d) if reionization implemented with reio_inter scheme */

  else if (pth->reio_parametrization == reio_inter) {

    /* this parametrization requires at least one point (z,xe) */
    class_test(pth->reio_inter_num<1,
               pth->error_message,
               "current implementation of reio_inter requires at least one point (z,xe)");

    /* this parametrization requires that the first z value is zero */
    class_test(pth->reio_inter_z[0] != 0.,
               pth->error_message,
               "For reio_inter scheme, the first value of reio_inter_z[...]  should always be zero, you passed %e",
               pth->reio_inter_z[0]);

    /* check that z input can be interpreted by the code */
    for (point=1; point<pth->reio_inter_num; point++) {
      class_test(pth->reio_inter_z[point-1]>=pth->reio_inter_z[point],
                 pth->error_message,
                 "value of reionization bin centers z_i expected to be passed in growing order, unlike: %e, %e",
                 pth->reio_inter_z[point-1],
                 pth->reio_inter_z[point]);
    }

    /* this parametrization requires that the last x_i value is zero
       (the code will substitute it with the value that one would get in
       absence of reionization, as compute by the recombination code) */
    class_test(pth->reio_inter_xe[pth->reio_inter_num-1] != 0.,
               pth->error_message,
               "For reio_inter scheme, the last value of reio_inter_xe[...]  should always be zero, you passed %e",
               pth->reio_inter_xe[pth->reio_inter_num-1]);

    /* copy here the (z,xe) values passed in input. */

    for (point=0; point<preio->reio_num_z; point++) {

      preio->reionization_parameters[preio->index_reio_first_z+point] = pth->reio_inter_z[point];

      /* check that xe input can be interpreted by the code */
      xe_input = pth->reio_inter_xe[point];
      if (xe_input >= 0.) {
        xe_actual = xe_input;
      }
      //-1 means "after hydrogen + first helium recombination"
      else if ((xe_input<-0.9) && (xe_input>-1.1)) {
        xe_actual = 1. + pth->YHe/(_not4_*(1.-pth->YHe));
      }
      //-2 means "after hydrogen + second helium recombination"
      else if ((xe_input<-1.9) && (xe_input>-2.1)) {
        xe_actual = 1. + 2.*pth->YHe/(_not4_*(1.-pth->YHe));
      }
      //other negative number is nonsense
      else {
        class_stop(pth->error_message,
                   "Your entry for reio_inter_xe[%d] is %e, this makes no sense (either positive or 0,-1,-2)",
                   point,pth->reio_inter_xe[point]);
      }

      preio->reionization_parameters[preio->index_reio_first_xe+point] = xe_actual;
    }

    /* copy highest redshift in reio_start */
    preio->reionization_parameters[preio->index_reio_start] = preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1];

    /* check it's not too big */
    class_test(preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max,
               pth->error_message,
               "starting redshift for reionization = %e, reionization_z_start_max = %e, you must change the binning or increase reionization_z_start_max",
               preio->reionization_parameters[preio->index_reio_start],
               ppr->reionization_z_start_max);

  }
  else if (pth->reio_parametrization == reio_none) {
  printf("Computing thermodynamics without reionization \n");
  }
  else{
    class_stop(pth->error_message,
               "Other reionization schemes not supported yet.");
  }
       
 return _SUCCESS_;

}

int thermodynamics_reionization_evolve_with_tau(
                                struct thermodynamics_parameters_and_workspace * ptpaw,
                                double mz_ini,
                                double mz_end,
                                double * mz_output,
                                int Nz
                                ) {
    
  /** - define local variables */

  int counter;
  double z_sup,z_mid,z_inf;
  double tau_sup,tau_mid,tau_inf;
  
  int index_tv;
  
  struct precision * ppr;
  struct background * pba;
  struct thermo * pth;
  struct recombination * preco;
  struct reionization * preio;
  struct thermo_workspace * ptw;

  /* function pointer to ODE evolver and names of possible evolvers */
  extern int evolver_rk();
  extern int evolver_ndf15();
  int (*generic_evolver)();
    

  /* Remame fields to avoid heavy notations */
  ppr = ptpaw->ppr;
  pba = ptpaw->pba;
  pth = ptpaw->pth;
  preco = ptpaw->preco;
  preio = ptpaw->preio;
  ptw = ptpaw->ptw;
   
  
  /* Initialiaze two thermo vectors to store initial conditions and one temporal vector for the calculations in the bisection*/
  struct thermo_vector * ptv; // Temporal vector as workspace
  struct thermo_vector * ptv_store; // Vector for storing the initial conditions
  
  ptv_store = ptw->tv;
  
  class_alloc(ptv,sizeof(struct thermo_vector),pth->error_message);
  
  class_define_index(ptv->index_Tmat,_TRUE_,index_tv,1);
  class_define_index(ptv->index_x_He,_TRUE_,index_tv,1);
  class_define_index(ptv->index_x_H,_TRUE_,index_tv,1);
  
  /* We have now obtained the full size */
  ptv->tv_size = index_tv;
    
  /* Allocate all arrays used during the evolution */
  class_calloc(ptv->y,ptv->tv_size,sizeof(double),pth->error_message);
  class_alloc(ptv->dy,ptv->tv_size*sizeof(double),pth->error_message);
  class_alloc(ptv->used_in_output,ptv->tv_size*sizeof(int),pth->error_message);
  
  for (index_tv=0; index_tv<ptv->tv_size; index_tv++){
    ptv->used_in_output[index_tv] = _TRUE_;
  }
  /* Set the new vector */ 
  ptv->y[ptv->index_Tmat] = ptw->tv->y[ptw->tv->index_Tmat];  
  ptv->dy[ptv->index_Tmat] = ptw->tv->dy[ptw->tv->index_Tmat];
  ptv->y[ptv->index_x_H] = ptw->tv->y[ptw->tv->index_x_H];
  ptv->dy[ptv->index_x_H] = ptw->tv->dy[ptw->tv->index_x_H];
  ptv->y[ptv->index_x_He] = ptw->tv->y[ptw->tv->index_x_He];
  ptv->dy[ptv->index_x_He] = ptw->tv->dy[ptw->tv->index_x_He];
  
  /*
  for (index_tv=0; index_tv<ptv->tv_size; index_tv++){
    printf("i=%i , y=%e dy=%e used %d \n", index_tv, ptw->tv->y[index_tv],ptw->tv->dy[index_tv],ptw->tv->used_in_output[index_tv]);
  }*/
  
  ptw->tv = ptv;
  
  /* upper value */
  
  z_sup = ppr->reionization_z_start_max-ppr->reionization_start_factor*pth->reionization_width;
  class_test(z_sup < 0.,
             pth->error_message,
             "parameters are such that reionization cannot take place before today while starting after z_start_max; need to increase z_start_max");
  
  /* maximum possible reionization redshift */
  preio->reionization_parameters[preio->index_reio_redshift] = z_sup;
  /* maximum possible starting redshift */
  preio->reionization_parameters[preio->index_reio_start] = ppr->reionization_z_start_max;
  /* infer xe_before_reio */
  
  
  if(ppr->evolver == rk){
    generic_evolver = evolver_rk;
  }
  else{
    generic_evolver = evolver_ndf15;
  }
  
  /* Calculate a first ionization history */ 
  class_call(generic_evolver(thermodynamics_derivs_with_recfast,
                             mz_ini,
                             mz_end,
                             ptv->y,
                             ptv->used_in_output,
                             ptv->tv_size,
                             ptpaw,
                             ppr->tol_thermo_integration,
                             ppr->smallest_allowed_variation,
                             thermodynamics_timescale_with_recfast,  // timescale
                             1., // stepsize
                             mz_output, // values of z for output
                             Nz, // size of previous array
                             thermodynamics_sources_with_recfast, // function for output
                             NULL, // print variables
                             pth->error_message),
             pth->error_message,
             pth->error_message);
  
  class_call(thermodynamics_reionization_get_tau(ppr,
                                                 pba,
                                                 pth,
                                                 preco,
                                                 preio),
             pth->error_message,
             pth->error_message);
    
  tau_sup=preio->reionization_optical_depth;

  class_test(tau_sup < pth->tau_reio,
             pth->error_message,
             "parameters are such that reionization cannot start after z_start_max");

  /* lower value */

  z_inf = 0.;
  tau_inf = 0.;
  
  /* Restore initial conditions */ 
  ptv->y[ptv->index_Tmat] = ptv_store->y[ptv_store->index_Tmat];  
  ptv->dy[ptv->index_Tmat] = ptv_store->dy[ptv_store->index_Tmat];
  ptv->y[ptv->index_x_H] = ptv_store->y[ptv_store->index_x_H];
  ptv->dy[ptv->index_x_H] = ptv_store->dy[ptv_store->index_x_H];
  ptv->y[ptv->index_x_He] = ptv_store->y[ptv_store->index_x_He];
  ptv->dy[ptv->index_x_He] = ptv_store->dy[ptv_store->index_x_He];

  /* try intermediate values */

  counter=0;
  while ((tau_sup-tau_inf) > pth->tau_reio * ppr->reionization_optical_depth_tol) {
      z_mid=0.5*(z_sup+z_inf);
      
      /* reionization redshift */
      preio->reionization_parameters[preio->index_reio_redshift] = z_mid;
      /* infer starting redshift for hygrogen */
      preio->reionization_parameters[preio->index_reio_start] = preio->reionization_parameters[preio->index_reio_redshift]+ppr->reionization_start_factor*pth->reionization_width;
      /* if starting redshift for helium is larger, take that one
       *    (does not happen in realistic models) */
      if(preio->reionization_parameters[preio->index_reio_start] < pth->helium_fullreio_redshift+ppr->reionization_start_factor*pth->helium_fullreio_width){      
         
          preio->reionization_parameters[preio->index_reio_start] = pth->helium_fullreio_redshift+ppr->reionization_start_factor*pth->helium_fullreio_width;
      }
      
      class_test(preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max,
                 pth->error_message,
                 "starting redshift for reionization > reionization_z_start_max = %e",ppr->reionization_z_start_max);
 
     
      /* Compute a new ionization history */
      
      class_call(generic_evolver(thermodynamics_derivs_with_recfast,
                                 mz_ini,
                                 mz_end,
                                 ptv->y,
                                 ptv->used_in_output,
                                 ptv->tv_size,
                                 ptpaw,
                                 ppr->tol_thermo_integration,
                                 ppr->smallest_allowed_variation,
                                 thermodynamics_timescale_with_recfast,  // timescale
                                 1., // stepsize
                                 mz_output, // values of z for output
                                 Nz, // size of previous array
                                 thermodynamics_sources_with_recfast, // function for output
                                 NULL, // print variables
                                 pth->error_message),
                 pth->error_message,
                 pth->error_message);
      
      /* Restore initial conditions */ 
      ptv->y[ptv->index_Tmat] = ptv_store->y[ptv_store->index_Tmat];  
      ptv->dy[ptv->index_Tmat] = ptv_store->dy[ptv_store->index_Tmat];
      ptv->y[ptv->index_x_H] = ptv_store->y[ptv_store->index_x_H];
      ptv->dy[ptv->index_x_H] = ptv_store->dy[ptv_store->index_x_H];
      ptv->y[ptv->index_x_He] = ptv_store->y[ptv_store->index_x_He];
      ptv->dy[ptv->index_x_He] = ptv_store->dy[ptv_store->index_x_He];
      
      
      class_call(thermodynamics_reionization_get_tau(ppr,
                                                     pba,
                                                     pth,
                                                     preco,
                                                     preio),
                 pth->error_message,
                 pth->error_message);
      
      tau_mid=preio->reionization_optical_depth;
      
      //printf("tau %d: %e, z_mid: %e tausup %e\n", counter+2,tau_mid,z_mid,tau_sup);

      /* trial */
      
      if (tau_mid > pth->tau_reio) {
          z_sup=z_mid;
          tau_sup=tau_mid;
      }
      else {
          z_inf=z_mid;
          tau_inf=tau_mid;
      }

      
         
      counter++;
      class_test(counter > _MAX_IT_,
                 pth->error_message,
                 "while searching for reionization_optical_depth, maximum number of iterations exceeded");
  }

  /* store z_reionization in thermodynamics structure */
  pth->z_reio=preio->reionization_parameters[preio->index_reio_redshift];
  
  
  class_call(thermo_vector_free(ptv),
             pth->error_message,
             pth->error_message);
  
  ptw->tv = ptv_store;
  
  
 return _SUCCESS_;

}

/*Routine to get the optical depth of reionization */

int thermodynamics_reionization_get_tau(
                                struct precision * ppr,
                                struct background * pba,
                                struct thermo * pth,
                                struct recombination * preco,
                                struct reionization * preio
                                ) {
 
  /* running index inside thermodynamics table */
  int i,integration_index;

  /** - --> search of index of reionization start in current table */
  i=0;
  while (preio->reionization_table[i*preio->re_size+preio->index_re_z] < preio->reionization_parameters[preio->index_reio_start]){
    i++;
    class_test(i == preio->rt_size,
               pth->error_message,
               "reionization_z_start_max = %e > largest redshift in thermodynamics table",ppr->reionization_z_start_max);
  }
 
  integration_index=i;
  
  /** - --> spline \f$ d \tau / dz \f$ with respect to z in view of integrating for optical depth between 0 and the just found starting index */
  class_call(array_spline(preio->reionization_table,
                          preio->re_size,
                          integration_index,
                          preio->index_re_z,
                          preio->index_re_dkappadz,
                          preio->index_re_d3kappadz3,
                          _SPLINE_EST_DERIV_,
                          pth->error_message),
             pth->error_message,
             pth->error_message);

  /** - --> integrate for optical depth */
  class_call(array_integrate_all_spline(preio->reionization_table,
                                        preio->re_size,
                                        integration_index,
                                        preio->index_re_z,
                                        preio->index_re_dkappadz,
                                        preio->index_re_d3kappadz3,
                                        &(preio->reionization_optical_depth),
                                        pth->error_message),
             pth->error_message,
             pth->error_message); 
  
 return _SUCCESS_;

}

/**
 * This routine is called in thermodynamics_solve_with_recfast
 * before the loop over the approximation schemes to set their boundary
 * redshifts in the vector interval_limit. Note that the the redshifts
 * are stored in negative values as we integrate in negative direction.
 *
 * Moreover the number of intervals is fixed to the number of approximation
 * schemes used.
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param preco      Input: pointer to the recombination structure
 * @param ptw        Input: pointer to perturb_workspace structure which fields are allocated or filled here
 * @param mz_ini     Input: negative redshift of the starting point of integration
 * @param mz_end     Input: negative redshift of the ending point of integration
 * @param interval_number Output: Number of intervals over which the evolver loop goes
 * @param interval_limit  Output: Vector which stores the negative redshifts which are the boundaries of each approximation scheme
 * @return the error status
 */

int thermodynamics_set_approximation_limits(
                                      struct precision * ppr,
                                      struct background * pba,
                                      struct thermo * pth,
                                      struct recombination * preco,
                                      struct thermo_workspace * ptw,
                                      double mz_ini,
                                      double mz_end,
                                      int* interval_number,
                                      double * interval_limit
                                      ){

  int index_ap; 
  
  /*fix interval number to number of approximations*/

  *interval_number = ptw->ap_size;

  class_test(pth->recombination == recfast && -mz_ini < ppr->recfast_z_He_3,
             pth->error_message,
             "increase zinitial, otherwise should get initial conditions from recfast's get_init routine (less precise anyway)");


  /*set limits for the intervals. Redshift is set negative for the evolver.*/
  
  /*integration starts at z_ini and ends at z_end*/	
  interval_limit[0]= mz_ini;
  interval_limit[ptw->ap_size] = mz_end;	
    
  /*each interval ends with the proper ending redshift of its approximation */
  for(index_ap=0; index_ap < ptw->ap_size-1; index_ap++){
    interval_limit[index_ap+1] = -ptw->ap_z_limits[index_ap];
  }
  
  class_test(interval_limit[ptw->ap_size-1] > mz_end,
             pth->error_message,
             "you want the second to last approximation scheme to end after z_end, that should not happen");
              
  
  for(index_ap=0; index_ap < ptw->ap_size; index_ap++){
    printf("Interval %i ending at %e \n",index_ap+1,-interval_limit[index_ap+1]);
  }
  
	

 return _SUCCESS_;

}

/**
 * This function is passed to the generic evolver and is called whenever
 * we want to store values for a given z that is passed in mz_output.
 * Depending on the current approximation scheme the ionization fraction
 * is either computed analytically, semi-analytically and from the 
 * (interpolated) output values of y. Moreover there is an automatic
 * smoothing enabled which smoothes out the the ionization_fraction 
 * after each approximation switch.
 *
 * This is one of the few functions in the code which is passed to
 * the generic_evolver() routine. Since generic_evolver()
 * should work with functions passed from various modules, the format
 * of the arguments is a bit special:
 *
 * - fixed parameters and workspaces are passed through a generic
 * pointer.  generic_evolver() doesn't know the content of this
 * pointer.
 *
 * - the error management is a bit special: errors are not written as
 * usual to pth->error_message, but to a generic error_message passed
 * in the list of arguments.
 *
 * @param mz                       Input: negative redshift
 * @param y                        Input: vector of thermodynamical quantities
 * @param dy                       Input: vector of redshift derivatives of theses quantities
 * @param index_z                  Input: index in the array mz_output
 * @param parameters_and_workspace Input/Output: in input, all parameters needed by thermodynamics_derivs_with_recfast; in output, recombination table
 * @param error_message            Output: error message
 * @return the error status
 */

int thermodynamics_sources_with_recfast(
                                        double mz,
                                        double * y,
                                        double * dy,
                                        int index_z,
                                        void * thermodynamics_parameters_and_workspace,
                                        ErrorMsg error_message
                    ) {
  int Nz;
  double z;
  
  double x,x_previous, weight,s,tau;
  int last_index_back;
  
  struct thermodynamics_parameters_and_workspace * ptpaw;
  struct precision * ppr;
  struct background * pba;
  struct thermo * pth;
  struct recombination * preco;
  struct reionization * preio;
  double * pvecback;
  struct thermo_workspace * ptw;
  int ap_current;

  /** - rename structure fields (just to avoid heavy notations) */
  ptpaw = thermodynamics_parameters_and_workspace;
  ppr = ptpaw->ppr;
  pba = ptpaw->pba;
  pth = ptpaw->pth;
  preco = ptpaw->preco;
  preio = ptpaw->preio;
  pvecback = ptpaw->pvecback;
  ptw = ptpaw->ptw;
  ap_current = ptw->ap_current;
  
  Nz = preio->rt_size;  
  z = -mz;
  
  
  class_call(background_tau_of_z(pba,
                                 z,
                                 &tau),
             pba->error_message,
             pth->error_message);

  class_call(background_at_tau(pba,
                               tau,
                               pba->short_info,
                               pba->inter_normal,
                               &last_index_back,
                               pvecback),
             pba->error_message,
             pth->error_message);

  class_test(pvecback[pba->index_bg_H] == 0.,
             pth->error_message,
             "stop to avoid division by zero");

    
  ptw->Tmat = y[ptw->tv->index_Tmat];
  ptw->dTmat = -dy[ptw->tv->index_Tmat];
  
  /* Depending on the current approximation scheme x has to be calculated for this specific z */
  if(ap_current == ptw->index_ap_H){
    class_call(thermodynamics_x_analytic(z,
                                         ppr,
                                         pth,
                                         preco,
                                         preio,
                                         ptw,
                                         ap_current),  
               error_message,
               error_message);
    
    x = ptw->x_H+preco->fHe*y[ptw->tv->index_x_He];
  }
  else if(ap_current == ptw->index_ap_frec){
    x = y[ptw->tv->index_x_H]+preco->fHe*y[ptw->tv->index_x_He];
  }
  else if(ap_current == ptw->index_ap_reio){
    ptw->x = y[ptw->tv->index_x_H]+preco->fHe*y[ptw->tv->index_x_He];
      
    class_call(thermodynamics_x_analytic(z,
                                         ppr,
                                         pth,
                                         preco,
                                         preio,
                                         ptw,
                                         ap_current),  
               error_message,
               error_message);
    
    x = ptw->x;
  }
  else if(ap_current == ptw->index_ap_reio_hyrec){
    /* get x from the recombination */
    class_call(thermodynamics_interpolate_recombination_table(ppr,
                                                              pth,
                                                              preco,
                                                              z,
                                                              &(ptw->x),
                                                              &(ptw->Tmat)),
               pth->error_message,
               pth->error_message);
  
    class_call(thermodynamics_x_analytic(z,
                                         ppr,
                                         pth,
                                         preco,
                                         preio,
                                         ptw,
                                         ptw->index_ap_reio),  
               error_message,
               error_message);
    
    x = ptw->x;
  }
  else{
    class_call(thermodynamics_x_analytic(z,
                                         ppr,
                                         pth,
                                         preco,
                                         preio,
                                         ptw,
                                         ap_current),  
               error_message,
               error_message);
  
    x = ptw->x;
  }
  /* Smoothing if we are shortly after an approximation switch, i.e. if z is within 2 delta after the switch*/
  if(ap_current != 0 && z > ptw->ap_z_limits[ap_current-1]-2*ptw->ap_z_limits_delta[ap_current]){
    
    class_call(thermodynamics_x_analytic(z,
                                         ppr,
                                         pth,
                                         preco,
                                         preio,
                                         ptw,
                                         ap_current-1),  
               error_message,
               error_message);
    if(ap_current-1 == ptw->index_ap_H){
      x_previous = ptw->x_H+preco->fHe*y[ptw->tv->index_x_He];
    }
    else if(ap_current-1 == ptw->index_ap_frec){
      x_previous = y[ptw->tv->index_x_H]+preco->fHe*y[ptw->tv->index_x_He];
    }
    else{
      x_previous = ptw->x;
    }
    /* get s from 0 to 1 */
    s = (ptw->ap_z_limits[ap_current-1]-z)/(2*ptw->ap_z_limits_delta[ap_current]);
    /* infer f2(x) = smooth function interpolating from 0 to 1 */
    weight = f2(s);


    x = weight*x+(1.-weight)*x_previous;
    //printf("sources z %e, x %e, ap %d, weight %e, x_previous %e \n", z,x,ap_current,weight,x_previous);
  }
  
  /** - --> store the results in the table */
  /* results are obtained in order of decreasing z, and stored in order of growing z */

  /* redshift */
  *(preio->reionization_table+(Nz-index_z-1)*preio->re_size+preio->index_re_z)=z;

  /* ionization fraction */
  *(preio->reionization_table+(Nz-index_z-1)*preio->re_size+preio->index_re_xe)=x;

  /* Tb */
  *(preio->reionization_table+(Nz-index_z-1)*preio->re_size+preio->index_re_Tb)=y[ptw->tv->index_Tmat];

  /* cb2 = (k_B/mu) Tb (1-1/3 dlnTb/dlna) = (k_B/mu) Tb (1+1/3 (1+z) dlnTb/dz) */
  *(preio->reionization_table+(Nz-index_z-1)*preio->re_size+preio->index_re_cb2)
    = _k_B_ / ( _c_ * _c_ * _m_H_ ) * (1. + (1./_not4_ - 1.) * preco->YHe + x * (1.-preco->YHe)) * y[ptw->tv->index_Tmat] * (1. - (1.+z) * dy[ptw->tv->index_Tmat] / y[ptw->tv->index_Tmat] / 3.);

  /* dkappa/dtau = a n_e x_e sigma_T = a^{-2} n_e(today) x_e sigma_T (in units of 1/Mpc) */
  *(preio->reionization_table+(Nz-index_z-1)*preio->re_size+preio->index_re_dkappadtau)
    = (1.+z) * (1.+z) * preco->Nnow * x * _sigma_ * _Mpc_over_m_;
  
  /* dkappa/dz = (dkappa/dtau) * (dtau/dz) = - (dkappa/dtau) / H  */
  *(preio->reionization_table+(Nz-index_z-1)*preio->re_size+preio->index_re_dkappadz)
    = (1.+z) * (1.+z) * preco->Nnow * x * _sigma_ * _Mpc_over_m_ / pvecback[pba->index_bg_H];
    
  //  printf("%e ,%e \n",x, *(preco->recombination_table+(Nz-index_z-1)*preco->re_size+preco->index_re_dkappadtau));
  
  return _SUCCESS_;

}

int thermodynamics_timescale_with_recfast(
                                        double z,
                                        void * thermodynamics_parameters_and_workspace,
                                        double * timescale,
                                        ErrorMsg error_message
                    ) {
  *timescale = 1.;
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

int thermodynamics_merge_reco_and_reio(
                                       struct precision * ppr,
                                       struct thermo * pth,
                                       struct recombination * preco,
                                       struct reionization * preio
                                       ) {
  /** Summary: */

  /** - define local variables */

  int i,index_th,index_re;

  /** - first, a little check that the two tables match each other and can be merged */
  
  if (pth->recombination == hyrec){ 
    if(pth->reio_parametrization != reio_none) {
      /** - --> look where to start in current recombination table */
      i=0;
      while (preco->recombination_table[i*preco->re_size+preco->index_re_z] < preio->reionization_table[(preio->rt_size -1)*preio->re_size+preio->index_re_z]) {
        i++;
        class_test(i == preco->rt_size,
                   pth->error_message,
                   "reionization_z_start_max = %e > largest redshift in thermodynamics table",ppr->reionization_z_start_max);
      }
      
      preio->index_reco_when_reio_start=i;   
      /*
       *  class_test(preco->recombination_table[preio->index_reco_when_reio_start*preco->re_size+preco->index_re_z] !=
       *             preio->reionization_table[(preio->rt_size -1)*preio->re_size+preio->index_re_z],
       *             pth->error_message,
       *             "mismatch which should never happen");*/
      
    }
    else{
      preio->rt_size = 0;
      preio->index_reco_when_reio_start = -1;
    }
  }
  else{
    /* When not in hyrec mode there is no recombination table */      
    preco->rt_size = 0;
    preio->index_reco_when_reio_start = -1;
  }

  /** - find number of redshift in full table = number in reco + number in reio - overlap */

  pth->tt_size = preco->rt_size + preio->rt_size - preio->index_reco_when_reio_start - 1;

  //printf("tt size %d \n", pth->tt_size);
  /** - allocate arrays in thermo structure */

  class_alloc(pth->z_table,pth->tt_size*sizeof(double),pth->error_message);
  class_alloc(pth->thermodynamics_table,pth->th_size*pth->tt_size*sizeof(double),pth->error_message);
  class_alloc(pth->d2thermodynamics_dz2_table,pth->th_size*pth->tt_size*sizeof(double),pth->error_message);

  /** - fill these arrays */

  for (i=0; i < preio->rt_size; i++) {
    pth->z_table[i]=
      preio->reionization_table[i*preio->re_size+preio->index_re_z];
    pth->thermodynamics_table[i*pth->th_size+pth->index_th_xe]=
      preio->reionization_table[i*preio->re_size+preio->index_re_xe];
    pth->thermodynamics_table[i*pth->th_size+pth->index_th_dkappa]=
      preio->reionization_table[i*preio->re_size+preio->index_re_dkappadtau];
    pth->thermodynamics_table[i*pth->th_size+pth->index_th_Tb]=
      preio->reionization_table[i*preio->re_size+preio->index_re_Tb];
    pth->thermodynamics_table[i*pth->th_size+pth->index_th_cb2]=
      preio->reionization_table[i*preio->re_size+preio->index_re_cb2];
  }
  for (i=0; i < preco->rt_size - preio->index_reco_when_reio_start - 1; i++) {
    index_th=i+preio->rt_size;
    index_re=i+preio->index_reco_when_reio_start+1;
    pth->z_table[index_th]=
      preco->recombination_table[index_re*preco->re_size+preco->index_re_z];
    pth->thermodynamics_table[index_th*pth->th_size+pth->index_th_xe]=
      preco->recombination_table[index_re*preco->re_size+preco->index_re_xe];
    pth->thermodynamics_table[index_th*pth->th_size+pth->index_th_dkappa]=
      preco->recombination_table[index_re*preco->re_size+preco->index_re_dkappadtau];
    pth->thermodynamics_table[index_th*pth->th_size+pth->index_th_Tb]=
      preco->recombination_table[index_re*preco->re_size+preco->index_re_Tb];
    pth->thermodynamics_table[index_th*pth->th_size+pth->index_th_cb2]=
      preco->recombination_table[index_re*preco->re_size+preco->index_re_cb2];
  } 

  /** - free the temporary structures */
  if(pth->recombination == hyrec){
    free(preco->recombination_table);
    
    if (pth->reio_parametrization != reio_none){
      free(preio->reionization_table);
    }
  }
  else{
    free(preio->reionization_table);    
  }
  
  
  
  return _SUCCESS_;
}

/**
 * Subroutine for formatting thermodynamics output
 */

int thermodynamics_output_titles(struct background * pba,
                                 struct thermo *pth,
                                 char titles[_MAXTITLESTRINGLENGTH_]
                                 ){

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
  class_store_columntitle(titles,"c_b^2",_TRUE_);
  class_store_columntitle(titles,"tau_d",_TRUE_);
  //class_store_columntitle(titles,"max. rate",_TRUE_);
  class_store_columntitle(titles,"r_d",pth->compute_damping_scale);

  return _SUCCESS_;
}

int thermodynamics_output_data(struct background * pba,
                               struct thermo *pth,
                               int number_of_titles,
                               double *data
                               ){

  int index_z, storeidx;
  double *dataptr, *pvecthermo;
  double z,tau;

  //  pth->number_of_thermodynamics_titles = get_number_of_titles(pth->thermodynamics_titles);
  //pth->size_thermodynamics_data = pth->number_of_thermodynamics_titles*pth->tt_size;


  /* Store quantities: */
  for (index_z=0; index_z<pth->tt_size; index_z++){
    dataptr = data + index_z*number_of_titles;
    pvecthermo = pth->thermodynamics_table+index_z*pth->th_size;
    z = pth->z_table[index_z];
    storeidx=0;

    class_call(background_tau_of_z(
                                   pba,
                                   z,
                                   &tau
                                   ),
               pba->error_message,
               pth->error_message);

    class_store_double(dataptr,z,_TRUE_,storeidx);
    class_store_double(dataptr,tau,_TRUE_,storeidx);
    class_store_double(dataptr,pvecthermo[pth->index_th_xe],_TRUE_,storeidx);
    class_store_double(dataptr,pvecthermo[pth->index_th_dkappa],_TRUE_,storeidx);
    //class_store_double(dataptr,pvecthermo[pth->index_th_ddkappa],_TRUE_,storeidx);
    //class_store_double(dataptr,pvecthermo[pth->index_th_dddkappa],_TRUE_,storeidx);
    class_store_double(dataptr,pvecthermo[pth->index_th_exp_m_kappa],_TRUE_,storeidx);
    class_store_double(dataptr,pvecthermo[pth->index_th_g],_TRUE_,storeidx);
    //class_store_double(dataptr,pvecthermo[pth->index_th_dg],_TRUE_,storeidx);
    //class_store_double(dataptr,pvecthermo[pth->index_th_ddg],_TRUE_,storeidx);
    class_store_double(dataptr,pvecthermo[pth->index_th_Tb],_TRUE_,storeidx);
    class_store_double(dataptr,pvecthermo[pth->index_th_cb2],_TRUE_,storeidx);
    class_store_double(dataptr,pvecthermo[pth->index_th_tau_d],_TRUE_,storeidx);
    //class_store_double(dataptr,pvecthermo[pth->index_th_rate],_TRUE_,storeidx);
    class_store_double(dataptr,pvecthermo[pth->index_th_r_d],pth->compute_damping_scale,storeidx);

  }

  return _SUCCESS_;
}

