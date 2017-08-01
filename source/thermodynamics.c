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

#ifdef COSMOREC
#include "CosmoRec.h"
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

    if ((((pth->reio_parametrization == reio_half_tanh) || (pth->reio_stars_and_dark_matter == _TRUE_))&& (z < 2*pth->z_reio))
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
  /** Initialize annihilation coefficient */



  if(pth->energy_repart_functions==Galli_et_al_interpolation || pth->energy_repart_functions==no_factorization){
    class_call(thermodynamics_annihilation_coefficients_init(ppr,pba,pth),
               pth->error_message,
               pth->error_message);
  }

  if(pth->has_on_the_spot==_FALSE_ && pth->energy_repart_functions!=no_factorization){
    // fprintf(stdout, "here\n" );
    class_call(thermodynamics_annihilation_f_eff_init(ppr,pba,preco),
               preco->error_message,
               preco->error_message);

  }
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

  // class_test((pth->annihilation_f_halo>0) && (pth->recombination==recfast),
  //            pth->error_message,
  //            "Switching on DM annihilation in halos requires using HyRec instead of RECFAST. Otherwise some values go beyond their range of validity in the RECFAST fits, and the thermodynamics module fails. Two solutions: add 'recombination = HyRec' to your input, or set 'annihilation_f_halo = 0.' (default).");

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

  // class_test((pth->annihilation_f_halo>0) && (pth->recombination==recfast),
  //            pth->error_message,
  //            "Switching on DM annihilation in halos requires using HyRec instead of RECFAST. Otherwise some values go beyond their range of validity in the RECFAST fits, and the thermodynamics module fails. Two solutions: add 'recombination = HyRec' to your input, or set 'annihilation_f_halo = 0.' (default).");

  class_test((pth->annihilation_f_halo<0),
             pth->error_message,
             "Parameter for DM annihilation in halos cannot be negative");

  class_test((pth->annihilation_z_halo<0),
             pth->error_message,
             "Parameter for DM annihilation in halos cannot be negative");

  if (pth->thermodynamics_verbose > 0)
    if ((pth->annihilation >0) && (pth->reio_parametrization == reio_none) && (ppr->recfast_Heswitch >= 3) && (pth->recombination==recfast))
      printf("Warning: if you have DM annihilation and you use recfast with option recfast_Heswitch >= 3, then the expression for CfHe_t and dy[1] becomes undefined at late times, producing nan's. This is however masked by reionization if you are not in reio_none mode.");

  class_test((pth->decay_fraction<0),
             pth->error_message,
             "decay parameter cannot be negative");

  class_test((pth->decay_fraction>0)&&(pba->has_cdm==_FALSE_),
             pth->error_message,
             "CDM decay effects require the presence of CDM!");

  /* tests in order to prevent segmentation fault in the following */
  class_test(_not4_ == 0.,
             pth->error_message,
             "stop to avoid division by zero");
  class_test(pth->YHe == 1.,
             pth->error_message,
             "stop to avoid division by zero");
 class_test(pth->alpha_asymmetric_planck_16<1.5 || pth->alpha_asymmetric_planck_16>50 ,pth->error_message,
   "alpha_asymmetric_planck_16 out of range [1.5,50]: rejected to avoid memory leakage.");
  /** - assign values to all indices in the structures with thermodynamics_indices()*/

  class_call(thermodynamics_indices(pth,preco,preio),
             pth->error_message,
             pth->error_message);

  /** - solve recombination and store values of \f$ z, x_e, d \kappa / d \tau, T_b, c_b^2 \f$ with thermodynamics_recombination() */

  class_call(thermodynamics_recombination(ppr,pba,pth,preco,pvecback),
             pth->error_message,
             pth->error_message);

  /** - if there is reionization, solve reionization and store values of \f$ z, x_e, d \kappa / d \tau, T_b, c_b^2 \f$ with thermodynamics_reionization()*/

  if ((pth->reio_parametrization != reio_none)) {
    class_call(thermodynamics_reionization(ppr,pba,pth,preco,preio,pvecback),
               pth->error_message,
               pth->error_message);
  }
  else {
    preio->rt_size=0;
    preio->index_reco_when_reio_start=-1;
  }

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

  /** - --> baryon drag interaction rate time minus one, -[R * kappa'], stored temporarily in column ddkappa */

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

    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_ddkappa] =
      -4./3.*pvecback[pba->index_bg_rho_g]/pvecback[pba->index_bg_rho_b]
      *pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa];

  }

  /** - --> second derivative of this rate, -[R * kappa']'', stored temporarily in column dddkappa */
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

  /** - --> compute r_d = [int_{tau_ini}^{tau} dtau [1/kappa'] */
  if (pth->compute_damping_scale == _TRUE_) {

    class_alloc(tau_table_growing,pth->tt_size*sizeof(double),pth->error_message);

    /* compute integrand 1/kappa' and store temporarily in column "ddkappa" */
    for (index_tau=0; index_tau < pth->tt_size; index_tau++) {

      tau_table_growing[index_tau]=tau_table[pth->tt_size-1-index_tau];

      pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_ddkappa] =
        1./pth->thermodynamics_table[(pth->tt_size-1-index_tau)*pth->th_size+pth->index_th_dkappa];

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
        factor, r_d^2 is equal to eta/(3 kappa'). So [r_d,ini^2] =
        [tau_ini/3/kappa'_ini] should be added to the integral in
        order to account for the integration between 0 and tau_ini */

     /* compute r_d */
     for (index_tau=0; index_tau < pth->tt_size; index_tau++) {

       pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_r_d] =
         sqrt(tau_table[pth->tt_size-1]/3./pth->thermodynamics_table[(pth->tt_size-1)*pth->th_size+pth->index_th_dkappa]
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
  // fprintf(stdout, "z_rec %e %e %e %e %e %e\n",pth->z_table[index_tau+1],pth->z_table[index_tau],pth->thermodynamics_table[(index_tau+2)*pth->th_size+pth->index_th_g],2.*pth->thermodynamics_table[(index_tau+1)*pth->th_size+pth->index_th_g],pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_g]);
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
    if((pth->reio_parametrization == reio_camb) || (pth->reio_parametrization == reio_half_tanh)) {
      if (pth->reio_z_or_tau==reio_tau)
        printf(" -> reionization  at z = %f\n",pth->z_reio);
      if (pth->reio_z_or_tau==reio_z)
        printf(" -> reionization with optical depth = %f\n",pth->tau_reio);
      class_call(background_tau_of_z(pba,pth->z_reio,&tau_reio),
                 pba->error_message,
                 pth->error_message);
      printf("    corresponding to conformal time = %f Mpc\n",tau_reio);
      printf("duration of reionization = %e, with z_beg = %e, z_mid = %e, z_end = %e\n",pth->duration_of_reionization,pth->z_10_percent,pth->z_50_percent, pth->z_99_percent);
    }
    if((pth->reio_parametrization == reio_douspis_et_al) || (pth->reio_parametrization == reio_asymmetric_planck_16)){
      printf(" -> reionization with optical depth = %f\n",pth->tau_reio);
      class_call(background_tau_of_z(pba,pth->z_reio,&tau_reio),
               pba->error_message,
               pth->error_message);
      printf("    corresponding to conformal time = %f Mpc\n",tau_reio);
      printf("duration of reionization = %e, with z_beg = %e, z_mid = %e, z_end = %e\n",pth->duration_of_reionization,pth->z_10_percent,pth->z_50_percent, pth->z_99_percent);
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

    preio->reio_num_params = index;

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

  if (pth->reio_parametrization == reio_douspis_et_al){

    preio->index_reio_xe_before = index;
    index++;
    preio->index_reio_xe_after = index;
    index++;
    preio->index_lambda_douspis_et_al = index;
    index++;
    preio->index_zp_douspis_et_al = index;
    index++;
    preio->index_Qp_douspis_et_al = index;
    index++;
    preio->index_helium_fullreio_fraction = index;
    index++;
    preio->index_helium_fullreio_redshift = index;
    index++;
    preio->index_helium_fullreio_width = index;
    index++;


    preio->reio_num_params = index;

  }

  if (pth->reio_parametrization == reio_asymmetric_planck_16){

    preio->index_reio_xe_before = index;
    index++;
    preio->index_reio_xe_after = index;
    index++;
    preio->index_alpha_asymmetric_planck_16 = index;
    index++;
    preio->index_z_end_asymmetric_planck_16= index;
    index++;
    preio->index_z_start_asymmetric_planck_16= index;
    index++;
    preio->index_helium_fullreio_fraction = index;
    index++;
    preio->index_helium_fullreio_redshift = index;
    index++;
    preio->index_helium_fullreio_width = index;
    index++;


    preio->reio_num_params = index;

  }

  if (pth->reio_parametrization == reio_stars_sfr_source_term){

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


    preio->reio_num_params = index;

  }

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

  /** - Since class v2.7 we make use of a fitting formula from PArthENoPE, Iocco et al. (2009), updated with the latest observational data
  on nuclear rates and neutron lifetime. It corresponds to the fit used by the Planck team, see e.g. Planck 15 cosmological papers. */

  double a[9]={0.2311,0.9502,-11.27,0.01356,0.008581,-0.1810,-0.0009795,-0.001370,0.01746}; //Coefficient of the fit.
  double b = 0.728; //Exponent of the lifetime in the fit.
  double tmp_YHe = 0.;
  int ii = 0, m = 0, n = 0;

  omega_b=pba->Omega0_b*pba->h*pba->h;

  for(ii=0;ii<9;ii++){
    if(ii == 0){
      m=0;
      n=0;
    }
    else if(ii == 3){
      m=1;
      n=0;
    }
    else if(ii==6){
      m=2;
      n=0;
    }
    tmp_YHe += pow(_NEUTRON_LIFETIME_/880.3,b)*a[ii]*pow(omega_b,n)*pow(DeltaNeff,m);
    n++;
  }

  pth->YHe = tmp_YHe;

  return _SUCCESS_;

}

/*****MODIF Vivian Poulin : Add new functions for energy repartition from DM annihilations or decays

Modification: Patrick Stöcker (20.02.17): Adding call to external script to calculate the annihilation coefficients on the fly.

*****/
int thermodynamics_annihilation_coefficients_init(
                                                  struct precision * ppr,
                                                  struct background * pba,
                                                  struct thermo * pth
                                                  ) {

  FILE * fA;
  char line[_LINE_LENGTH_MAX_];
  char * left;

  /* BEGIN: New variables related to the use of an external code to calculate the annihilation coefficients */
  char arguments[_ARGUMENT_LENGTH_MAX_];
  char command_with_arguments[2*_ARGUMENT_LENGTH_MAX_];
  int status;
  /* END */

  int num_lines=0;
  int array_line=0;

  /*

      the following file is assumed to contain (apart from comments and blank lines):
     - One number (num_lines) = number of lines of the file
     - six columns (xe , chi_heat, chi_Lya, chi_H, chi_He, chi_lowE) where each chi represents the fraction of energy going respectively into
     heat, excitation of lyman-alpha level, Hydrogen ionisation, Helium ionisation, photons below 10.2 eV unseeable by the IGM.

  */

  /* BEGIN: Add switch (1) */
  if (!(ppr->fz_is_extern)) {
    class_open(fA,ppr->energy_injec_coeff_file, "r",pth->error_message);
  } else {
    /* Prepare the command */
    /* Pass the list of arguments */
    sprintf(arguments, "%g %g %g %g %g", ppr->param_fz_1, ppr->param_fz_2, ppr->param_fz_3, ppr->param_fz_4, ppr->param_fz_5);
    /* Write the actual command */
    sprintf(command_with_arguments, "%s %s", ppr->command_fz, arguments);
    free(ppr->command_fz);
    if (pth->thermodynamics_verbose > 0)
      printf(" -> running: %s\n", command_with_arguments);
    /* Launch the process and retrieve the output */

    fA = popen(command_with_arguments, "r");
    class_test(fA == NULL, pth->error_message, "The program failed to set the environment for the external command.");
  }

  /* END */

  /* go through each line */
  while (fgets(line,_LINE_LENGTH_MAX_-1,fA) != NULL) {
    /* eliminate blank spaces at beginning of line */
    left=line;
    while (left[0]==' ') {
      left++;
    }

    /* check that the line is neither blank nor a comment. In
       ASCII, left[0]>39 means that first non-blank charachter might
       be the beginning of some data (it is not a newline, a #, a %,
       etc.) */
    if (left[0] > 39) {

      /* if the line contains data, we must interprete it. If
         num_lines == 0 , the current line must contain
         its value. Otherwise, it must contain (xe , chi_heat, chi_Lya, chi_H, chi_He, chi_lowE). */
      if (num_lines == 0) {

        /* read num_lines, infer size of arrays and allocate them */
        class_test(sscanf(line,"%d",&num_lines) != 1,
                   pth->error_message,
                   "could not read value of parameters num_lines in file %s\n",ppr->energy_injec_coeff_file);
        class_alloc(pth->annihil_coef_xe,num_lines*sizeof(double),pth->error_message);
        class_alloc(pth->annihil_coef_heat,num_lines*sizeof(double),pth->error_message);
        class_alloc(pth->annihil_coef_lya,num_lines*sizeof(double),pth->error_message);
        class_alloc(pth->annihil_coef_ionH,num_lines*sizeof(double),pth->error_message);
        class_alloc(pth->annihil_coef_ionHe,num_lines*sizeof(double),pth->error_message);
        class_alloc(pth->annihil_coef_lowE,num_lines*sizeof(double),pth->error_message);

        class_alloc(pth->annihil_coef_dd_heat,num_lines*sizeof(double),pth->error_message);
        class_alloc(pth->annihil_coef_dd_lya,num_lines*sizeof(double),pth->error_message);
        class_alloc(pth->annihil_coef_dd_ionH,num_lines*sizeof(double),pth->error_message);
        class_alloc(pth->annihil_coef_dd_ionHe,num_lines*sizeof(double),pth->error_message);
        class_alloc(pth->annihil_coef_dd_lowE,num_lines*sizeof(double),pth->error_message);
        pth->annihil_coef_num_lines = num_lines;


        array_line=0;

      }
      else {

        /* read coefficients */
        class_test(sscanf(line,"%lg %lg %lg %lg %lg %lg",
                          &(pth->annihil_coef_xe[array_line]),
                          &(pth->annihil_coef_heat[array_line]),
                          &(pth->annihil_coef_lya[array_line]),
                          &(pth->annihil_coef_ionH[array_line]),
                          &(pth->annihil_coef_ionHe[array_line]),
                          &(pth->annihil_coef_lowE[array_line])) != 6,
                   pth->error_message,
                   "could not read value of parameters coeeficients in file %s\n\nThe line I do not understand is: >> %s <<\n",ppr->energy_injec_coeff_file,line);
        array_line ++;
      }
    }
  }
  /* BEGIN: Add switch (2) */
  if (!(ppr->fz_is_extern)) {
    fclose(fA);
  } else {
    status = pclose(fA);
	class_test(status != 0., pth->error_message, "The attempt to launch the external command was not successful. Maybe the output of the external command is not in the right format.");
  }
  /* END */

  /* spline in one dimension */
  class_call(array_spline_table_lines(pth->annihil_coef_xe,
                                      num_lines,
                                      pth->annihil_coef_heat,
                                      1,
                                      pth->annihil_coef_dd_heat,
                                      _SPLINE_NATURAL_,
                                      pth->error_message),
             pth->error_message,
             pth->error_message);

  class_call(array_spline_table_lines(pth->annihil_coef_xe,
                                      num_lines,
                                      pth->annihil_coef_lya,
                                      1,
                                      pth->annihil_coef_dd_lya,
                                      _SPLINE_NATURAL_,
                                      pth->error_message),
             pth->error_message,
             pth->error_message);

  class_call(array_spline_table_lines(pth->annihil_coef_xe,
                                      num_lines,
                                      pth->annihil_coef_ionH,
                                      1,
                                      pth->annihil_coef_dd_ionH,
                                      _SPLINE_NATURAL_,
                                      pth->error_message),
             pth->error_message,
             pth->error_message);
  class_call(array_spline_table_lines(pth->annihil_coef_xe,
                                     num_lines,
                                     pth->annihil_coef_ionHe,
                                     1,
                                     pth->annihil_coef_dd_ionHe,
                                     _SPLINE_NATURAL_,
                                     pth->error_message),
              pth->error_message,
              pth->error_message);
  class_call(array_spline_table_lines(pth->annihil_coef_xe,
                                      num_lines,
                                      pth->annihil_coef_lowE,
                                      1,
                                      pth->annihil_coef_dd_lowE,
                                      _SPLINE_NATURAL_,
                                      pth->error_message),
               pth->error_message,
               pth->error_message);

  return _SUCCESS_;

}
/**
 * This function is used by the energy injection module for two different interpolations:
 * Either to directly interpolate the f(z) functions per channels or to interpolate
 * the chi(x_e) functions  per channels when the factorisation approximation is assumed.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param pth Input/Output: pointer to initialized thermo structure
 * @return the error status
 */

int thermodynamics_annihilation_coefficients_interpolate(
                                                         struct precision * ppr,
                                                         struct background * pba,
                                                         struct thermo * pth,
                                                         double xe_or_z
                                                         ) {

  int last_index;

  class_call(array_interpolate_spline(pth->annihil_coef_xe,
                                      pth->annihil_coef_num_lines,
                                      pth->annihil_coef_heat,
                                      pth->annihil_coef_dd_heat,
                                      1,
                                      xe_or_z,
                                      &last_index,
                                      &(pth->chi_heat),
                                      1,
                                      pth->error_message),
             pth->error_message,
             pth->error_message);

  class_call(array_interpolate_spline(pth->annihil_coef_xe,
                                      pth->annihil_coef_num_lines,
                                      pth->annihil_coef_lya,
                                      pth->annihil_coef_dd_lya,
                                      1,
                                      xe_or_z,
                                      &last_index,
                                      &(pth->chi_lya),
                                      1,
                                      pth->error_message),
             pth->error_message,
             pth->error_message);

  class_call(array_interpolate_spline(pth->annihil_coef_xe,
                                      pth->annihil_coef_num_lines,
                                      pth->annihil_coef_ionH,
                                      pth->annihil_coef_dd_ionH,
                                      1,
                                      xe_or_z,
                                      &last_index,
                                      &(pth->chi_ionH),
                                      1,
                                      pth->error_message),
             pth->error_message,
             pth->error_message);
    class_call(array_interpolate_spline(pth->annihil_coef_xe,
                                       pth->annihil_coef_num_lines,
                                       pth->annihil_coef_ionHe,
                                       pth->annihil_coef_dd_ionHe,
                                       1,
                                       xe_or_z,
                                       &last_index,
                                       &(pth->chi_ionHe),
                                       1,
                                       pth->error_message),
              pth->error_message,
              pth->error_message);
    class_call(array_interpolate_spline(pth->annihil_coef_xe,
                                       pth->annihil_coef_num_lines,
                                       pth->annihil_coef_lowE,
                                       pth->annihil_coef_dd_lowE,
                                       1,
                                       xe_or_z,
                                       &last_index,
                                       &(pth->chi_lowE),
                                       1,
                                       pth->error_message),
              pth->error_message,
              pth->error_message);

  return _SUCCESS_;

}


int thermodynamics_annihilation_coefficients_free(
                                                  struct thermo * pth
                                                  ) {

  free(pth->annihil_coef_xe);
  free(pth->annihil_coef_heat);
  free(pth->annihil_coef_lya);
  free(pth->annihil_coef_ionH);
  free(pth->annihil_coef_ionHe);
  free(pth->annihil_coef_lowE);

  free(pth->annihil_coef_dd_heat);
  free(pth->annihil_coef_dd_lya);
  free(pth->annihil_coef_dd_ionH);
  free(pth->annihil_coef_dd_ionHe);
  free(pth->annihil_coef_dd_lowE);

  return _SUCCESS_;

}
// /****************MODIF Vivian Poulin 2 : Add f(z) functions in halos****************/
int thermodynamics_annihilation_f_eff_init(
                                                  struct precision * ppr,
                                                  struct background * pba,
                                                  struct recombination * preco
                                                  ) {

  FILE * fA;
  char line[_LINE_LENGTH_MAX_];
  char * left;

  int num_lines=0;
  int array_line=0;

  /*

      the following file is assumed to contain (apart from comments and blank lines):
     - One number (num_lines) = number of lines of the file
     - One column (z , f(z)) where f(z) represents the "effective" fraction of energy deposited into the medium at redshift z, in presence of halo formation.

  */
  class_open(fA,ppr->energy_injec_f_eff_file, "r",preco->error_message);

  /* go through each line */
  while (fgets(line,_LINE_LENGTH_MAX_-1,fA) != NULL) {

    /* eliminate blank spaces at beginning of line */
    left=line;
    while (left[0]==' ') {
      left++;
    }

    /* check that the line is neither blank nor a comment. In
       ASCII, left[0]>39 means that first non-blank charachter might
       be the beginning of some data (it is not a newline, a #, a %,
       etc.) */
    if (left[0] > 39) {

      /* if the line contains data, we must interprete it. If
         num_lines == 0 , the current line must contain
         its value. Otherwise, it must contain (xe , chi_heat, chi_Lya, chi_H, chi_He, chi_lowE). */
      if (num_lines == 0) {

        /* read num_lines, infer size of arrays and allocate them */
        class_test(sscanf(line,"%d",&num_lines) != 1,
                   preco->error_message,
                   "could not read value of parameters num_lines in file %s\n",ppr->energy_injec_f_eff_file);
        class_alloc(preco->annihil_z,num_lines*sizeof(double),preco->error_message);
        class_alloc(preco->annihil_f_eff,num_lines*sizeof(double),preco->error_message);

        class_alloc(preco->annihil_dd_f_eff,num_lines*sizeof(double),preco->error_message);

        preco->annihil_f_eff_num_lines = num_lines;


        array_line=0;

      }
      else {

        /* read coefficients */
        class_test(sscanf(line,"%lg %lg",
                          &(preco->annihil_z[array_line]),
                          &(preco->annihil_f_eff[array_line]))!= 2,
                   preco->error_message,
                   "could not read value of parameters coefficients in file %s\n",ppr->energy_injec_f_eff_file);
        array_line ++;
      }
    }
  }

  fclose(fA);

  /* spline in one dimension */
  class_call(array_spline_table_lines(preco->annihil_z,
                                      num_lines,
                                      preco->annihil_f_eff,
                                      1,
                                      preco->annihil_dd_f_eff,
                                      _SPLINE_NATURAL_,
                                      preco->error_message),
             preco->error_message,
             preco->error_message);


  return _SUCCESS_;

}


int thermodynamics_annihilation_f_eff_interpolate(
                                                  struct precision * ppr,
                                                  struct background * pba,
                                                  struct recombination * preco,
                                                  double z
                                                ) {

  int last_index;
  class_call(array_interpolate_spline(preco->annihil_z,
                                      preco->annihil_f_eff_num_lines,
                                      preco->annihil_f_eff,
                                      preco->annihil_dd_f_eff,
                                      1,
                                      z,
                                      &last_index,
                                      &(preco->f_eff),
                                      1,
                                      preco->error_message),
             preco->error_message,
             preco->error_message);


  return _SUCCESS_;

}


int thermodynamics_annihilation_f_eff_free(
                                                  struct recombination * preco
                                                  ) {

  free(preco->annihil_z);
  free(preco->annihil_f_eff);
  free(preco->annihil_dd_f_eff);


  return _SUCCESS_;

}
/********** END OF MODIFICATION By Vivian Poulin **************/
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

/**********************************************************************************************/
/******************************Energy Injection DM annihilation**********************************/
int thermodynamics_DM_annihilation_energy_injection(
                                                  struct precision * ppr,
                                                  struct background * pba,
                                                  struct recombination * preco,
                                                  double z,
                                                  double * energy_rate,
                                                  ErrorMsg error_message
                                                ){

  double rho_cdm_today;
  double Boost_factor;

  rho_cdm_today = pow(pba->H0*_c_/_Mpc_over_m_,2)*3/8./_PI_/_G_*(pba->Omega0_cdm)*_c_*_c_; /* energy density in J/m^3 */

  if(preco->annihilation_z_halo>0.){
  Boost_factor = preco->annihilation_f_halo*erfc((1+z)/(1+preco->annihilation_z_halo))/pow(1+z,3);
  }
  else Boost_factor = 0;

  *energy_rate = pow(rho_cdm_today,2)/_c_/_c_*(pow((1.+z),6)*preco->annihilation)*(1+Boost_factor);
  /* energy density rate in J/m^3/s (remember that sigma_thermal/(preco->annihilation_m_DM*conversion) is in m^3/s/Kg) */

}
/******************************Energy Injection DM decay**********************************/
int thermodynamics_DM_decay_energy_injection(
                                                  struct precision * ppr,
                                                  struct background * pba,
                                                  struct recombination * preco,
                                                  double z,
                                                  double * energy_rate,
                                                  ErrorMsg error_message
                                                ){
  double rho_cdm_today, rho_dcdm,decay_factor;
  double tau;
  int last_index_back;
  double * pvecback;
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);

  if(pba->Omega_ini_dcdm!=0 || pba->Omega0_dcdmdr !=0){
    class_call(background_tau_of_z(pba,
                                   z,
                                   &tau),
               pba->error_message,
               ppr->error_message);

    class_call(background_at_tau(pba,
                                 tau,
                                 pba->long_info,
                                 pba->inter_normal,
                                 &last_index_back,
                                 pvecback),
               pba->error_message,
               ppr->error_message);
     rho_dcdm = pvecback[pba->index_bg_rho_dcdm]*pow(_Mpc_over_m_,2)*3/8./_PI_/_G_*_c_*_c_; /* energy density in J/m^3 */
     /* If uncommented, these lines allow to check approximation when computing the dcdm density with analytical results. Works very well until Omega_lambda dominates, then ~10% difference. */
    //  result_integrale = exp(-pba->Gamma_dcdm*2*((pba->Omega0_b+pba->Omega0_cdm)*pow(pba->Omega0_g+(pba->Omega0_b+pba->Omega0_cdm)/(1+z),0.5)
    //  +2*pow(pba->Omega0_g,1.5)*(1+z)-2*pba->Omega0_g*pow((1+z)*(pba->Omega0_g*(1+z)+(pba->Omega0_b+pba->Omega0_cdm)),0.5))/(3*pow((pba->Omega0_b+pba->Omega0_cdm),2)*(1+z)*pba->H0));
    //   rho_dcdm_approchee = pow(pba->H0*_c_/_Mpc_over_m_,2)*3/8./_PI_/_G_*(pba->Omega_ini_dcdm)*_c_*_c_*result_integrale*pow(1+z,3);
    //   fprintf(stdout, "z = %e vrai = %e  approchee = %e relativ diff = %e\n",z,rho_dcdm, rho_dcdm_approchee,(rho_dcdm-rho_dcdm_approchee)/rho_dcdm_approchee);
  }
  else{
    if(preco->has_on_the_spot == _FALSE_)decay_factor=1; //The effect of the exponential decay is already incoporated within the f_z functions.
    else {
    class_call(background_tau_of_z(pba,
                                   z,
                                   &tau),
               pba->error_message,
               ppr->error_message);

    class_call(background_at_tau(pba,
                                 tau,
                                 pba->long_info,
                                 pba->inter_normal,
                                 &last_index_back,
                                 pvecback),
               pba->error_message,
               ppr->error_message);
    decay_factor = exp(-pba->Gamma_dcdm*pvecback[pba->index_bg_time]);
    }
    rho_cdm_today = pow(pba->H0*_c_/_Mpc_over_m_,2)*3/8./_PI_/_G_*(pba->Omega0_cdm)*_c_*_c_; /* energy density in J/m^3 */
    rho_dcdm = rho_cdm_today*pow((1+z),3)*decay_factor; // This trick avoid mixing gravitational and electromagnetic impacts of the decay on the CMB power spectra.
  }


  *energy_rate = rho_dcdm*preco->decay_fraction*(pba->Gamma_dcdm*_c_/_Mpc_over_m_);
  // fprintf(stdout, "*energy_rate %e\n",*energy_rate );
  free(pvecback);


}
/******************************Energy Injection low mass PBH (evaporation)**********************************/
int thermodynamics_low_mass_pbh_energy_injection(
                                                  struct precision * ppr,
                                                  struct background * pba,
                                                  struct recombination * preco,
                                                  double z,
                                                  double * energy_rate,
                                                  ErrorMsg error_message
                                                ){

  double rho_cdm_today;
  double tau;
  int last_index_back;
  double * pvecback;
  //Parameters related to PBH

  double f,f_neutrinos, em_branching,tau_pbh,pbh_mass;
  double dMdt;
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);
  rho_cdm_today = pow(pba->H0*_c_/_Mpc_over_m_,2)*3/8./_PI_/_G_*(pba->Omega0_cdm)*_c_*_c_; /* energy density in J/m^3 */

    f_neutrinos = 6*0.147;
    if(preco->PBH_low_mass<1e17) f = (2*0.06+6*0.147+2*0.007+2*2*0.142);//electrons
    if(preco->PBH_low_mass<1e15) f = (2*0.06+6*0.147+2*0.007+2*2*0.142+2*0.007+2*2*0.142);//muons
    if(preco->PBH_low_mass<1e14) f = (2*0.06+6*0.147+2*0.007+2*2*0.142+2*0.007+2*2*0.142+2*2*0.142);//tau
    if(preco->PBH_low_mass<5e13) f = (2*0.06+6*0.147+2*0.007+2*2*0.142+2*0.007+2*2*0.142+2*2*0.142+3*12*0.142 + 16*0.06);//u d s and gluons      zinitial = preco->z_tmp;
    em_branching = (f-f_neutrinos)/f;


    if(preco->PBH_low_mass<1e15){
      tau_pbh = 407*pow(f/15.35,-1)*pow(preco->PBH_low_mass/(1e10),3); //not used, for comparison only
      //fprintf(stdout, "T %e fem %e br_ee %e br_gg %e tau_pbh %e exp %e \n",1.06*(pow(10,10)/preco->PBH_low_mass)*pow(10,12),f-f_neutrinos,2*2*0.142/(f-f_neutrinos),2*0.06/(f-f_neutrinos),tau_pbh,exp(1e14/tau_pbh));
      class_call(background_tau_of_z(pba,
                                     z,
                                     &tau),
                 pba->error_message,
                 ppr->error_message);

      class_call(background_at_tau(pba,
                                   tau,
                                   pba->long_info,
                                   pba->inter_normal,
                                   &last_index_back,
                                   pvecback),
                 pba->error_message,
                 ppr->error_message);
      pbh_mass = pow(pow(preco->PBH_low_mass,3)-3*5.34e-5*f*1e30*pvecback[pba->index_bg_time]/(_c_ / _Mpc_over_m_),1./3);
      if(pbh_mass <= 0){
        pbh_mass = 0;
        dMdt = 0;
      }
      else dMdt=5.34e-5*f*pow(pbh_mass/1e10,-2)*1e10;
      *energy_rate = rho_cdm_today*pow((1+z),3)*preco->PBH_fraction/preco->PBH_low_mass*em_branching*(dMdt);
      // if(tau_pbh!=0)*energy_rate = rho_cdm_today*pow((1+z),3)*preco->PBH_fraction/tau_pbh*result_integrale_2;
      if(isnan(*energy_rate)==1)*energy_rate=0.;
      // fprintf(stdout, "energy_rate %e A %e B %e z %e\n",*energy_rate,rho_cdm_today*pow((1+z),3)*preco->PBH_fraction/preco->PBH_low_mass*(dMdt),z);
    }
    else {
      // fprintf(stdout, "em_branching %e \n",em_branching );
      tau_pbh = 407*pow(f/15.35,-1)*pow(preco->PBH_low_mass/(1e10),3); // if preco->PBH_low_mass>1e15, we assume that the BH (effectively) does not lose mass since its lifetime is bigger than the age of the Universe.
      *energy_rate = rho_cdm_today*pow((1+z),3)*preco->PBH_fraction/tau_pbh*em_branching;
      if(isnan(*energy_rate)==1)*energy_rate=0.;

    }
    free(pvecback);

}
/******************************Energy Injection high mass PBH (accretion)**********************************/
int thermodynamics_high_mass_pbh_energy_injection(
                                                  struct precision * ppr,
                                                  struct background * pba,
                                                  struct recombination * preco,
                                                  double z,
                                                  double * energy_rate,
                                                  ErrorMsg error_message
                                                ){

  double rho_cdm_today;
  double tau;
  int last_index_back;
  double * pvecback;
  //Parameters related to PBH
  double c_s, v_eff,v_eff_2,v_l, r_B,x_e,beta,beta_eff,beta_hat,x_cr,lambda,n_gas,M_b_dot,M_sun,M_ed_dot,epsilon,L_acc,Integrale,Normalization;
  double m_H, m_dot, m_dot_2, L_acc_2,L_ed,l,l2,M_crit;
  double rho, m_p = 938, m_e = 0.511, T_infinity = 0, rho_infinity = 0, x_e_infinity = 0, P_infinity = 0, rho_cmb = 0, t_B = 0, v_B = 0;
  double lambda_1,lambda_2,lambda_ad,lambda_iso,gamma_cooling,beta_compton_drag, T_s, T_ion, Y_s, J,tau_cooling;

  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);
  rho_cdm_today = pow(pba->H0*_c_/_Mpc_over_m_,2)*3/8./_PI_/_G_*(pba->Omega0_cdm)*_c_*_c_; /* energy density in J/m^3 */
        class_call(background_tau_of_z(pba,
                                       z,
                                       &tau),
                   pba->error_message,
                   preco->error_message);

        class_call(background_at_tau(pba,
                                     tau,
                                     pba->long_info,
                                     pba->inter_normal,
                                     &last_index_back,
                                     pvecback),
                   pba->error_message,
                   preco->error_message);


        c_s = 5.7e3*pow(preco->Tm_tmp/2730,0.5);//conversion km en m
        M_sun = 2e30; // in Kg
        n_gas = 200*1e6*pow((1+z)/1000,3); // 1e6 = conversion cm^-3 en m^-3;
        m_H= 1.67e-27; // Hydrogen mass in kg

        // x_e = 1;
        x_e = preco->xe_tmp;
        T_infinity = preco->Tm_tmp*_eV_over_Kelvin_*1e-6; //Temperature in MeV
        M_ed_dot = 1.44e17*(preco->PBH_high_mass)*1e-3; // 1e-3 = conversion g en Kg;
        L_ed = M_ed_dot*_c_*_c_; // J s^-1;


        if(preco->PBH_accretion_recipe == Ricotti_et_al || preco->PBH_accretion_recipe == Horowitz){

            if(preco->PBH_accretion_recipe == Ricotti_et_al){
              if(z<180)v_eff = pow(1+z,0.63403)*185.806;// Result of a fit on fig. 2 of Ricotti et al. 0709.0524
              if(z>=180)v_eff = pow(1+z,0.221672)*1581.39;// Result of a fit on fig. 2 of Ricotti et al. 0709.0524
            }
            else if(preco->PBH_accretion_recipe == Horowitz){


              // v_B = sqrt((1+x_e)*T_infinity/m_p)*_c_;
              v_B = 5.7e3*pow(preco->Tm_tmp/2730,0.5);//conversion km en m
              v_l = 31*MIN(1,z/1000)*1e3;
              if(v_B < v_l){
                 v_eff = sqrt(v_B*v_l);
                 v_eff *= sqrt(1+(v_l/v_B,2));

               }
              else {
                v_eff = v_B;
                v_eff *= pow(16/sqrt(2*_PI_)*pow(v_l/v_B,3),1./6);
                // fprintf(stdout, "v_B %e v_l %e z %e\n",v_B,v_l,z);
              }
              // v_eff = sqrt(v_B*v_l);
              // if(z<1500)v_eff = pow(1+z,0.793943)*0.0541752*1e3;// Result of a fit on fig. 7 of Ali-Haimoud et al. 1612.05644 // in m
              // if(z>=1500)v_eff = pow(1+z,0.21987)*3.64188*1e3;// Result of a fit on fig. 7 of Ali-Haimoud et al. 1612.05644  // in m
              // *sqrt(1+(v_l/v_B,2));
              // v_eff = 31000;
              // if(v_eff/c_s>1)v_eff = c_s*pow(16/sqrt(2*_PI_)*pow(v_eff/c_s,3),1./6); // eq 6 from Horowitz 1612.07264
              // else v_eff = c_s*pow(1+pow(v_eff/c_s,2),1./2);
              // fprintf(stdout, "v_eff %e z %e T_infinity %e x_e %e \n",v_eff,z,T_infinity,x_e);
              // if(z>=1000)v_eff = 31000;
            }

            // //First way of computing m_dot and L_acc (from eq. 22 of Ricotti et al.)
            // r_B = _G_*preco->PBH_high_mass*M_sun*pow(v_eff,-2);
            // beta = 2.06e-23*x_e*pow(1+z,4);
            // beta_eff = beta+pvecback[pba->index_bg_H]*_c_/_Mpc_over_m_;
            // beta_hat = beta_eff*r_B/c_s;
            // x_cr = (-1+pow(1+beta_hat,0.5))/beta_hat;
            // lambda = exp(4.5/(3+pow(beta_hat,0.75)))*x_cr*x_cr;
            // m_dot = (2.8e-3*lambda)*pow((1+z)/1000,3)*(preco->PBH_high_mass)*pow(v_eff/5.74e3,-3); //1.8 -> 2.8 mistake in Ricotti et al.
            // // if(isnan(m_dot)==1)m_dot = 0;
            // // L_ed = 4*_PI_*_G_*preco->PBH_high_mass*M_sun*m_p*1e6/_eV_over_joules_/(_sigma_*_c_);
            // M_ed_dot = 1.44e17*(preco->PBH_high_mass)*1e-3; // 1e-3 = conversion g en Kg;
            // L_ed = M_ed_dot*_c_*_c_; // J s^-1;
            // if(preco->PBH_high_mass>100) l = preco->PBH_fraction*MIN(0.1*m_dot,1);
            // else l = 0.011*m_dot*m_dot;
            // L_acc_2 = l*L_ed;

            //Second way of computing m_dot and L_acc (from eq. 1 of Ricotti et al.)
            r_B = _G_*preco->PBH_high_mass*M_sun*pow(v_eff,-2);
            beta = 2.06e-23*x_e*pow(1+z,4);
            beta_eff = beta+pvecback[pba->index_bg_H]*_c_/_Mpc_over_m_;
            // fprintf(stdout, "z %e 1/beta %e 1/H %e 1/beta_eff %e t_B %e \n",z,1/beta,1/(pvecback[pba->index_bg_H]*_c_/_Mpc_over_m_),1/beta_eff,r_B/v_eff);
            beta_hat = beta_eff*r_B/c_s;
            x_cr = (-1+pow(1+beta_hat,0.5))/beta_hat;
            lambda = exp(4.5/(3+pow(beta_hat,0.75)))*x_cr*x_cr;
            M_b_dot = 4*_PI_*lambda*m_H*n_gas*v_eff*r_B*r_B;
            M_ed_dot = 1.44e17*(preco->PBH_high_mass)*1e-3; // 1e-3 = conversion g en Kg;
            m_dot_2 = M_b_dot/M_ed_dot;
            // m_dot_old = M_b_dot_old/M_ed_dot;
            if(preco->PBH_high_mass>100) l2 = preco->PBH_fraction*MIN(0.1*m_dot_2,1);
            else l2 = 0.011*m_dot_2*m_dot_2;
            L_acc_2 = l2*L_ed;

            // fprintf(stdout, "z %e m_dot %e L_acc_2 %e   \n",z,m_dot,L_acc_2);
          }

        //Third way of computing m_dot and L_acc from Gaggero et al. arXiv:1612.00457
        else if(preco->PBH_accretion_recipe == Gaggero_et_al){
            //v_eff not clearly given, I use the alternative v_eff by Ali-Haimoud et al.
            v_B = sqrt((1+x_e)*T_infinity/m_p)*_c_;
            v_l = 30*MIN(1,z/1000)*1e3;
            if(v_B < v_l) v_eff = sqrt(v_B*v_l);
            else v_eff = v_B;
            lambda = 0.01;
            rho = pvecback[pba->index_bg_rho_b]/pow(_Mpc_over_m_,2)*3/8./_PI_/_G_*_c_*_c_; /* energy density in kg/m^3 */
            M_b_dot = 4*_PI_*lambda*pow(_G_*preco->PBH_high_mass*M_sun,2)*rho*pow(v_eff,-3.);

            M_crit = 0.01*4*_PI_*_G_*preco->PBH_high_mass*M_sun*m_p*1e6/_eV_over_joules_/(_sigma_*_c_)/(_c_*_c_); //1% of the eddington accretion rate.
            L_acc_2 = 0.3*0.1*M_b_dot*M_b_dot*_c_*_c_/M_crit; // 1.1e15 = Eddington accretion rate for a 7 solar mass black hole in kg s^-1
            L_ed = 4*_PI_*_G_*preco->PBH_high_mass*M_sun*m_p*1e6/_eV_over_joules_/(_sigma_*_c_);

            // fprintf(stdout, "z %e M_crit %e M_b_dot/Medd %e L_acc_2/Ledd %e   \n",z,M_crit,M_b_dot/(100*M_crit),L_acc_2/(0.3*L_ed));

          }
        //Fourth way of computing m_dot and L_acc from Ali-Haimoud et al. 1612.05644
        else if(preco->PBH_accretion_recipe == Ali_Haimoud){
          rho_cmb = pvecback[pba->index_bg_rho_g]/pow(_Mpc_over_m_,2)*3/8./_PI_/_G_*_c_*_c_*_c_*_c_* 6.241509e12; /* energy density in MeV/m^3 */
          // x_e_infinity = 1; // change to 1 for the strong-feedback case
          x_e_infinity = x_e; // change to x_e for the no-feedback case
          v_B = sqrt((1+x_e_infinity)*T_infinity/m_p)*_c_;
          v_l = 30*MIN(1,z/1000)*1e3;
          if(v_B < v_l) v_eff = sqrt(v_B*v_l);
          else v_eff = v_B;
          // if(v_eff == 0) v_eff = pow(1+z,0.21987)*3.64188*1e3;
          // fprintf(stdout, " z %e x_e %e T_infinity %e v_B %e v_l %e v_eff %e\n",z,x_e,T_infinity,v_B,v_l,v_eff);
          // fprintf(stdout, "z %e v_l %e sqrt(5/3)*v_B %e c_s %e \n", z,v_l,sqrt(5./3)*v_B,5.7e3*pow(preco->Tm_tmp/2730,0.5));
          // if(z<1500)v_eff_2 = pow(1+z,0.793943)*0.0541752*1e3;// Result of a fit on fig. 7 of Ali-Haimoud et al. 1612.05644 // in m
          // if(z>=1500)v_eff_2 = pow(1+z,0.21987)*3.64188*1e3;// Result of a fit on fig. 7 of Ali-Haimoud et al. 1612.05644  // in m
          r_B = _G_*preco->PBH_high_mass*M_sun*pow(v_eff,-2); // in m
          t_B = _G_*preco->PBH_high_mass*M_sun/pow(v_eff,3); // in s
          // fprintf(stdout, "z %e T_infinity %e rho_infinity %e P_infinity %e v_B %e t_B %e x_e_infinity %e \n",z,T_infinity,rho_infinity,P_infinity,v_B,t_B,x_e_infinity );
          beta_compton_drag = 4./3*x_e_infinity*_sigma_*rho_cmb*t_B/(m_p)*_c_; //Misprint in Ali-Haimoud et al., c should be in numerator otherwise not dimensionless
          gamma_cooling = 2*m_p/(m_e*(1+x_e_infinity))*beta_compton_drag;
          // fprintf(stdout, "z %e 1/beta %e 1/H %e 1/gamma_cooling %e t_B %e \n",z,1/beta_compton_drag,1/(pvecback[pba->index_bg_H]*_c_/_Mpc_over_m_),1/gamma_cooling,t_B);
          lambda_iso = 0.25*exp(1.5);
          lambda_ad = 0.25*pow(3./5,1.5);
          lambda_1 = lambda_ad+(lambda_iso-lambda_ad)*pow(gamma_cooling*gamma_cooling/(88+gamma_cooling*gamma_cooling),0.22);
          lambda_2 = exp(4.5/(3+pow(beta_compton_drag,0.75)))*1/(pow(pow(1+beta_compton_drag,0.5)+1,2));
          // lambda_2 = exp(4.5/(3+pow(beta_compton_drag,0.75)))*(pow(1+beta_compton_drag,0.5)-1)/beta_compton_drag;
          // fprintf(stdout, "z %e beta_compton_drag %e gamma_cooling %e lambda_1 %e lambda_2 %e lambda_3 %e \n",z,beta_compton_drag,gamma_cooling,lambda_1,lambda_2,exp(4.5/(3+pow(beta_compton_drag,0.75)))*(pow(1+beta_compton_drag,0.5)-1)/beta_compton_drag);
          lambda = lambda_1*lambda_2/lambda_iso;
          rho = pvecback[pba->index_bg_rho_b]/pow(_Mpc_over_m_,2)*3/8./_PI_/_G_*_c_*_c_; /* energy density in kg/m^3 */
          // fprintf(stdout, "rho/m_H %e n_gas %e \n", 0.75*rho/m_H,n_gas);
          M_b_dot = 4*_PI_*lambda*rho*r_B*r_B*v_eff; //in kg s^-1
          T_ion = 1.5e4*_eV_over_Kelvin_;
          tau_cooling = 1.5/(5+pow(gamma_cooling,2./3));
          Y_s = pow((1+x_e_infinity)/2,2./3*13.6/T_ion)*tau_cooling/4*pow(1-5./2*tau_cooling,1./3)*m_p/m_e;
          T_s = m_e * Y_s*pow(1+Y_s/0.27,-1./3); // in MeV
          if(T_s/m_e > 1)  J = 27/(2*_PI_)*(log(2*T_s/(m_e)*exp(-0.577)+0.08)+4./3);
          else J = 4/_PI_*sqrt(2/_PI_)*pow(T_s/m_e,-0.5)*(1+5.5*pow(T_s/m_e,1.25));
          // fprintf(stdout, "z %e J %e T_s %e Y_s %e  tau_cooling %e \n", z,J,T_s,Y_s,tau_cooling);
          L_ed = 4*_PI_*_G_*preco->PBH_high_mass*M_sun*m_p*1e6/_eV_over_joules_/(_sigma_*_c_);
          L_acc_2 = 1./137*T_s/(m_p)*J*pow(M_b_dot*_c_*_c_,2)/L_ed;
          // fprintf(stdout, "z %e M_crit %e M_b_dot %e L_acc_2 %e   \n",z,M_crit,M_b_dot,L_acc_2);
        }
        else if (preco->PBH_accretion_recipe == Hybrid){
          //Hybrid accretion model starting from a spherical accretion followed by a transition to disk accretion. We smooth the transition thanks to a tanh function.
          //In that case we follow spherical accretion recipe from AliHaimoud&Kamionkowski.
          rho_cmb = pvecback[pba->index_bg_rho_g]/pow(_Mpc_over_m_,2)*3/8./_PI_/_G_*_c_*_c_*_c_*_c_* 6.241509e12; /* energy density in MeV/m^3 */
          x_e_infinity = x_e; // change to x_e for the no-feedback case
          v_B = sqrt((1+x_e_infinity)*T_infinity/m_p)*_c_;
          v_l = 30*MIN(1,z/1000)*1e3;
          if(v_B < v_l) v_eff = sqrt(v_B*v_l);
          else v_eff = v_B;
          r_B = _G_*preco->PBH_high_mass*M_sun*pow(v_eff,-2); // in m
          t_B = _G_*preco->PBH_high_mass*M_sun/pow(v_eff,3); // in s
          beta_compton_drag = 4./3*x_e_infinity*_sigma_*rho_cmb*t_B/(m_p)*_c_; //Misprint in Ali-Haimoud et al., c should be in numerator otherwise not dimensionless
          gamma_cooling = 2*m_p/(m_e*(1+x_e_infinity))*beta_compton_drag;
          lambda_iso = 0.25*exp(1.5);
          lambda_ad = 0.25*pow(3./5,1.5);
          lambda_1 = lambda_ad+(lambda_iso-lambda_ad)*pow(gamma_cooling*gamma_cooling/(88+gamma_cooling*gamma_cooling),0.22);
          lambda_2 = exp(4.5/(3+pow(beta_compton_drag,0.75)))*1/(pow(pow(1+beta_compton_drag,0.5)+1,2)); // as defined in AliHaimoud&Kamionkowski
          // lambda_2 = exp(4.5/(3+pow(beta_compton_drag,0.75)))*(pow(1+beta_compton_drag,0.5)-1)/beta_compton_drag; // as defined in Ricotti et al.
          lambda = lambda_1*lambda_2/lambda_iso;
          rho = pvecback[pba->index_bg_rho_b]/pow(_Mpc_over_m_,2)*3/8./_PI_/_G_*_c_*_c_; /* energy density in kg/m^3 */
          M_b_dot = 4*_PI_*lambda*rho*r_B*r_B*v_eff; //in kg s^-1
          T_ion = 1.5e4*_eV_over_Kelvin_;
          tau_cooling = 1.5/(5+pow(gamma_cooling,2./3));
          Y_s = pow((1+x_e_infinity)/2,2./3*13.6/T_ion)*tau_cooling/4*pow(1-5./2*tau_cooling,1./3)*m_p/m_e;
          T_s = m_e * Y_s*pow(1+Y_s/0.27,-1./3); // in MeV
          if(T_s/m_e > 1)  J = 27/(2*_PI_)*(log(2*T_s/(m_e)*exp(-0.577)+0.08)+4./3);
          else J = 4/_PI_*sqrt(2/_PI_)*pow(T_s/m_e,-0.5)*(1+5.5*pow(T_s/m_e,1.25));
          L_ed = 4*_PI_*_G_*preco->PBH_high_mass*M_sun*m_p*1e6/_eV_over_joules_/(_sigma_*_c_);
          L_acc = 1./137*T_s/(m_p)*J*pow(M_b_dot*_c_*_c_,2)/L_ed;
          L_acc = L_acc*MAX(tanh((z-preco->PBH_disk_formation_redshift)/5),0);
          M_crit = 0.01*4*_PI_*_G_*preco->PBH_high_mass*M_sun*m_p*1e6/_eV_over_joules_/(_sigma_*_c_)/(_c_*_c_); //1% of the eddington accretion rate.

          //Here we follow Gaggero et al. recipe for disk accretion, based on observation and very conservative.
          lambda = 0.01;
          rho = pvecback[pba->index_bg_rho_b]/pow(_Mpc_over_m_,2)*3/8./_PI_/_G_*_c_*_c_; /* energy density in kg/m^3 */
          M_b_dot = 4*_PI_*lambda*pow(_G_*preco->PBH_high_mass*M_sun,2)*rho*pow(v_eff,-3.);
          M_crit = 0.01*4*_PI_*_G_*preco->PBH_high_mass*M_sun*m_p*1e6/_eV_over_joules_/(_sigma_*_c_)/(_c_*_c_); //1% of the eddington accretion rate.
          L_acc_2 = 0.3*0.1*M_b_dot*M_b_dot*_c_*_c_/M_crit; // 1.1e15 = Eddington accretion rate for a 7 solar mass black hole in kg s^-1
          L_acc_2 = L_acc_2*MIN((1-tanh((z-preco->PBH_disk_formation_redshift)/5)),1);
          L_acc_2 = L_acc+L_acc_2;

        }


        *energy_rate =  (rho_cdm_today/(preco->PBH_high_mass*M_sun*_c_*_c_))*pow(1+z,3)*L_acc_2*preco->PBH_fraction;
        // fprintf(stdout, "%e %e %e %e %e %e %e %e %e %e %e %e\n",z, beta_compton_drag, gamma_cooling,lambda,M_b_dot*_c_*_c_/L_ed,T_s*1e6/_eV_over_Kelvin_,T_s*J/m_p/137,v_eff,v_B,v_l,L_acc_2/L_ed,*energy_rate);
        // fprintf(stdout, "%e %e %e %e %e %e %e %e \n",x_e, M_b_dot,lambda,m_dot_2,l2,L_acc_2,*energy_rate,z);
        // fprintf(stdout, "%e %e %e \n", z,m_dot,L_acc_2);
        free(pvecback);

}



/*************************New version, corrected by Vivian Poulin******************************/
int thermodynamics_onthespot_energy_injection(
                                              struct precision * ppr,
                                              struct background * pba,
                                              struct recombination * preco,
                                              double z,
                                              double * energy_rate,
                                              ErrorMsg error_message
                                              ) {

  if(preco->annihilation > 0){
    thermodynamics_DM_annihilation_energy_injection(ppr,pba,preco,z,energy_rate,error_message);
  }
  if(preco->decay_fraction > 0.){
    thermodynamics_DM_decay_energy_injection(ppr,pba,preco,z,energy_rate,error_message);
  }
    // fprintf(stdout, "z = %e energy_rate = %e\n", z, *energy_rate);
  if(preco->PBH_high_mass > 0.){
    thermodynamics_high_mass_pbh_energy_injection(ppr,pba,preco,z,energy_rate,error_message);
  }
  if(preco->PBH_low_mass > 0.){
    thermodynamics_low_mass_pbh_energy_injection(ppr,pba,preco,z,energy_rate,error_message);
  }
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
  double exponent_z,exponent_zp;
  if (preco->annihilation > 0 || preco->decay_fraction > 0 || preco->PBH_high_mass > 0 || preco->PBH_low_mass > 0 ) {
    if (preco->has_on_the_spot == _FALSE_) {


      if(preco->energy_deposition_treatment == Analytical_approximation){
        // // /* number of hydrogen nuclei today in m**-3 */
        nH0 = 3.*preco->H0*preco->H0*pba->Omega0_b/(8.*_PI_*_G_*_m_H_)*(1.-preco->YHe);

        /*Value from Poulin 1508.01370*/
        /* factor = c sigma_T n_H(0) / (H(0) \sqrt(Omega_m)) (dimensionless) */
        factor = _sigma_ * nH0 / pba->H0 * _Mpc_over_m_ / sqrt(pba->Omega0_b+pba->Omega0_cdm);
        exponent_z = 8;
        exponent_zp = 7.5;

        /*Value from Ali-Haimoud & Kamionkowski 1612.05644*/
        // factor = 0.1*_sigma_ * nH0 / pba->H0 * _Mpc_over_m_ / sqrt(pba->Omega0_b+pba->Omega0_cdm);
        // exponent_z = 7;
        // exponent_zp = 6.5;


        /* integral over z'(=zp) with step dz */
        dz=1.;

        /* first point in trapezoidal integral */
        zp = z;
        class_call(thermodynamics_onthespot_energy_injection(ppr,pba,preco,zp,&onthespot,error_message),
                   error_message,
                   error_message);
        first_integrand = factor*pow(1+z,exponent_z)/pow(1+zp,exponent_zp)*exp(2./3.*factor*(pow(1+z,1.5)-pow(1+zp,1.5)))*onthespot; // beware: versions before 2.4.3, there were rwrong exponents: 6 and 5.5 instead of 7 and 7.5
        result = 0.5*dz*first_integrand;

        /* other points in trapezoidal integral */
        do {

          zp += dz;
          class_call(thermodynamics_onthespot_energy_injection(ppr,pba,preco,zp,&onthespot,error_message),
                     error_message,
                     error_message);
          integrand = factor*pow(1+z,exponent_z)/pow(1+zp,exponent_zp)*exp(2./3.*factor*(pow(1+z,1.5)-pow(1+zp,1.5)))*onthespot; // beware: versions before 2.4.3, there were rwrong exponents: 6 and 5.5 instead of 7 and 7.5
          result += dz*integrand;

        } while (integrand/first_integrand > 0.02);
      }

      // // /***********************************************************************************************************************/
      else if(preco->energy_deposition_treatment == Slatyer){

            if(preco->energy_repart_functions!=no_factorization){
              class_call(thermodynamics_annihilation_f_eff_interpolate(ppr,pba,preco,z),
                        preco->error_message,
                        preco->error_message);
              preco->f_eff=MAX(preco->f_eff,0.);
            }
            else preco->f_eff=1.;

            class_call(thermodynamics_onthespot_energy_injection(ppr,pba,preco,z,&result,error_message),
                      error_message,
                      error_message);
            *energy_rate =  result*preco->f_eff;
            // fprintf(stdout, "energy_rate %e\n", *energy_rate);
      }
      /* uncomment these lines if you also want to compute the on-the-spot for comparison */
      // class_call(thermodynamics_onthespot_energy_injection(ppr,pba,preco,z,&onthespot,error_message),
      //            error_message,
      //            error_message);
      // /* these test lines print the energy rate rescaled by (1+z)^6 in J/m^3/s, with or without the on-the-spot approximation */
      // fprintf(stdout,"%e  %e  %e  %e\n",
      // 1.+z,
      // result/pow(1.+z,6),
      // onthespot/pow(1.+z,6),result/onthespot);
    }
    else {
      class_call(thermodynamics_onthespot_energy_injection(ppr,pba,preco,z,&result,error_message),
                 error_message,
                 error_message);
    }




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

int thermodynamics_reionization_function(
                                         double z,
                                         struct thermo * pth,
                                         struct reionization * preio,
                                         struct recombination * preco,
                                         double * xe
                                         ) {
  /** Summary: */
  /** - define local variables */
  double argument, A, B, factor;
  int i;
  double z_jump;
  int jump;
  double center,before, after,width,one_jump;
  double x_tmp,z_tmp;
  /** - implementation of ionization function similar to the one in CAMB */
  if ((pth->reio_parametrization == reio_camb) || (pth->reio_parametrization == reio_half_tanh)) {
    /** -> case z > z_reio_start */

    if(z > preio->reionization_parameters[preio->index_reio_start]) {
      *xe = preio->reionization_parameters[preio->index_reio_xe_before];

    }

    else {

      /** - --> case z < z_reio_start: hydrogen contribution (tanh of complicated argument) */
      /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */

      A = (pow((1.+preio->reionization_parameters[preio->index_reio_redshift]),
                    preio->reionization_parameters[preio->index_reio_exponent])
                - pow((1.+z),preio->reionization_parameters[preio->index_reio_exponent]));
      // A = MAX (A,0.);
      B = (preio->reionization_parameters[preio->index_reio_exponent]
        *pow((1.+preio->reionization_parameters[preio->index_reio_redshift]),
             (preio->reionization_parameters[preio->index_reio_exponent]-1.)))*preio->reionization_parameters[preio->index_reio_width];
             if(B == 0)argument = 0;
             else argument = A / B;



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

        /** -> case z < z_reio_start: helium contribution (tanh of simpler argument) */
        /********Modified by Vivian Poulin to take into account helium reionization in both cases***********/
          argument = (preio->reionization_parameters[preio->index_helium_fullreio_redshift] - z)
            /preio->reionization_parameters[preio->index_helium_fullreio_width];
          /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
          *xe += preio->reionization_parameters[preio->index_helium_fullreio_fraction]
            * (tanh(argument)+1.)/2.;
      // if (pth->reio_parametrization == reio_camb) {
      //   argument = (preio->reionization_parameters[preio->index_helium_fullreio_redshift] - z)
      //     /preio->reionization_parameters[preio->index_helium_fullreio_width];
      //   /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
      //   *xe += preio->reionization_parameters[preio->index_helium_fullreio_fraction]
      //     * (tanh(argument)+1.)/2.;
      // }

      if(*xe > 0.1*(preio->reionization_parameters[preio->index_reio_xe_after]-preio->reionization_parameters[preio->index_reio_xe_before]) && pth->z_10_percent == 0){
        pth->z_10_percent = z;
        // fprintf(stdout, "pth->z_10_percent %e\n", pth->z_10_percent);
      }
      if(*xe > 0.50*(preio->reionization_parameters[preio->index_reio_xe_after]-preio->reionization_parameters[preio->index_reio_xe_before]) && pth->z_50_percent == 0 ){
        pth->z_50_percent = z;
        // fprintf(stdout, "pth->z_50_percent %e\n", pth->z_50_percent);

      }
      if(*xe > 0.99*(preio->reionization_parameters[preio->index_reio_xe_after]-preio->reionization_parameters[preio->index_reio_xe_before]) && pth->z_99_percent == 0){
      // if(*xe > 0.99*(preio->reionization_parameters[preio->index_reio_xe_after]-preio->reionization_parameters[preio->index_reio_xe_before]) && pth->z_99_percent > 200){
        pth->z_99_percent = z;
        // class_test(pth->z_99_percent<=6,
        //            pth->error_message,
        //            "z_99_percent < 6, we reject the point");

      }
      if(pth->z_99_percent != 0 && pth->z_10_percent != 0 && pth->duration_of_reionization == 0){
        pth->duration_of_reionization = pth->z_10_percent  - pth->z_99_percent;
        if(pth->duration_of_reionization < 0) pth->duration_of_reionization ==0;
        // class_test(pth->duration_of_reionization<1,
        //            pth->error_message,
        //            "duration_of_reionization < 1, we reject the point");
        // fprintf(stdout, "pth->duration_of_reionization %e v2 %e \n", pth->duration_of_reionization,pth->z_10_percent  -preio->reionization_parameters[preio->index_z_end_asymmetric_planck_16]);

      }


      x_tmp = *xe;
      z_tmp = z;
    }




    // fprintf(stdout, "xe %e z%e \n ", *xe, z);
    // fprintf(stdout, "xe %e z%e argument %e A %e B %e A / B %e \n", *xe, z, argument,A,B, A / B);
    /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
    return _SUCCESS_;

  }

  if(pth->reio_parametrization == reio_douspis_et_al){
    if(z > preio->reionization_parameters[preio->index_reio_start]) {
      *xe = preio->reionization_parameters[preio->index_reio_xe_before];

    }

    else {

    if(z < preio->reionization_parameters[preio->index_zp_douspis_et_al]){
      factor = (1-preio->reionization_parameters[preio->index_Qp_douspis_et_al])/(pow(1+preio->reionization_parameters[preio->index_zp_douspis_et_al],3)-1)
      *(pow(1+preio->reionization_parameters[preio->index_zp_douspis_et_al],3)-pow(1+z,3))
      +preio->reionization_parameters[preio->index_Qp_douspis_et_al];

    *xe =  (preio->reionization_parameters[preio->index_reio_xe_after]
             -preio->reionization_parameters[preio->index_reio_xe_before])*
             factor;

    }
    else{
      factor = preio->reionization_parameters[preio->index_Qp_douspis_et_al]*exp(-preio->reionization_parameters[preio->index_lambda_douspis_et_al]
        *pow(z-preio->reionization_parameters[preio->index_zp_douspis_et_al],3)/(pow(z-preio->reionization_parameters[preio->index_zp_douspis_et_al],2)+0.2));
      *xe = (preio->reionization_parameters[preio->index_reio_xe_after]
               -preio->reionization_parameters[preio->index_reio_xe_before])*
               factor;
      *xe = MAX(preio->reionization_parameters[preio->index_reio_xe_before],*xe);

    }

    /** -> case z < z_reio_start: helium contribution (tanh of simpler argument) */
    /********Helium reionization is taken into account with a camb-like parametrization***********/
      argument = (preio->reionization_parameters[preio->index_helium_fullreio_redshift] - z)
        /preio->reionization_parameters[preio->index_helium_fullreio_width];
      /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
      *xe += preio->reionization_parameters[preio->index_helium_fullreio_fraction]
        * (tanh(argument)+1.)/2.;
    // fprintf(stdout, "z %e x_e %e xe_before %e factor %elambda %e zp %e qp %e\n", z,*xe,preio->reionization_parameters[preio->index_reio_xe_before],factor,preio->reionization_parameters[preio->index_lambda_douspis_et_al],preio->reionization_parameters[preio->index_zp_douspis_et_al],preio->reionization_parameters[preio->index_Qp_douspis_et_al]);

    }
    if(*xe > 0.1*(preio->reionization_parameters[preio->index_reio_xe_after]-preio->reionization_parameters[preio->index_reio_xe_before]) && pth->z_10_percent == 0){
      pth->z_10_percent = z;
      // fprintf(stdout, "pth->z_10_percent %e\n", pth->z_10_percent);
    }
    if(*xe > 0.50*(preio->reionization_parameters[preio->index_reio_xe_after]-preio->reionization_parameters[preio->index_reio_xe_before]) && pth->z_50_percent == 0 ){
      pth->z_50_percent = z;
      // fprintf(stdout, "pth->z_50_percent %e\n", pth->z_50_percent);

    }
    if(*xe > 0.99*(preio->reionization_parameters[preio->index_reio_xe_after]-preio->reionization_parameters[preio->index_reio_xe_before]) && pth->z_99_percent == 0){
    // if(*xe > 0.99*(preio->reionization_parameters[preio->index_reio_xe_after]-preio->reionization_parameters[preio->index_reio_xe_before]) && pth->z_99_percent > 200){
      pth->z_99_percent = z;
    //   class_test(pth->z_99_percent<=6,
    //              pth->error_message,
    //              "z_99_percent < 6, we reject the point");
    //  fprintf(stdout, "pth->z_99_percent %e  \n", pth->z_99_percent);


    }
    if(pth->z_99_percent != 0 && pth->z_10_percent != 0 && pth->duration_of_reionization == 0){
      pth->duration_of_reionization = pth->z_10_percent  - pth->z_99_percent;
      if(pth->duration_of_reionization < 0) pth->duration_of_reionization ==0;
      // class_test(pth->duration_of_reionization<1,
      //            pth->error_message,
      //            "duration_of_reionization < 1, we reject the point");
      // fprintf(stdout, "pth->duration_of_reionization %e v2 %e \n", pth->duration_of_reionization,pth->z_10_percent  -preio->reionization_parameters[preio->index_z_end_asymmetric_planck_16]);

    }

    return _SUCCESS_;

  }
  if(pth->reio_parametrization == reio_asymmetric_planck_16){
    if(z > preio->reionization_parameters[preio->index_reio_start]) {
      *xe = preio->reionization_parameters[preio->index_reio_xe_before];

    }

    else {

    if(z < preio->reionization_parameters[preio->index_z_end_asymmetric_planck_16]){


    *xe =  (preio->reionization_parameters[preio->index_reio_xe_after]
             -preio->reionization_parameters[preio->index_reio_xe_before]);

    }
    else{
      factor = pow((preio->reionization_parameters[preio->index_reio_start]-z)/(preio->reionization_parameters[preio->index_reio_start]-preio->reionization_parameters[preio->index_z_end_asymmetric_planck_16]),preio->reionization_parameters[preio->index_alpha_asymmetric_planck_16]);
      *xe = (preio->reionization_parameters[preio->index_reio_xe_after]
               -preio->reionization_parameters[preio->index_reio_xe_before])*
               factor;
      *xe = MAX(preio->reionization_parameters[preio->index_reio_xe_before],*xe);

    }

    if(*xe > 0.1*(preio->reionization_parameters[preio->index_reio_xe_after]-preio->reionization_parameters[preio->index_reio_xe_before]) && pth->z_10_percent == 0){
      pth->z_10_percent = z;
      // fprintf(stdout, "pth->z_10_percent %e\n", pth->z_10_percent);
    }
    if(*xe > 0.50*(preio->reionization_parameters[preio->index_reio_xe_after]-preio->reionization_parameters[preio->index_reio_xe_before]) && pth->z_50_percent == 0 ){
      pth->z_50_percent = z;
      // fprintf(stdout, "pth->z_50_percent %e\n", pth->z_50_percent);

    }
    if(*xe > 0.99*(preio->reionization_parameters[preio->index_reio_xe_after]-preio->reionization_parameters[preio->index_reio_xe_before]) && pth->z_99_percent == 0){
    // if(*xe > 0.99*(preio->reionization_parameters[preio->index_reio_xe_after]-preio->reionization_parameters[preio->index_reio_xe_before]) && pth->z_99_percent > 200){
      pth->z_99_percent = z;
      // class_test(pth->z_99_percent<=6,
      //            pth->error_message,
      //            "z_99_percent < 6, we reject the point");


    }
    if(pth->z_99_percent != 0 && pth->z_10_percent != 0 && pth->duration_of_reionization == 0){
      pth->duration_of_reionization = pth->z_10_percent  - pth->z_99_percent;
      if(pth->duration_of_reionization < 0) pth->duration_of_reionization ==0;
      // class_test(pth->duration_of_reionization<1,
      //            pth->error_message,
      //            "duration_of_reionization < 1, we reject the point");
      // fprintf(stdout, "pth->duration_of_reionization %e v2 %e \n", pth->duration_of_reionization,pth->z_10_percent  -preio->reionization_parameters[preio->index_z_end_asymmetric_planck_16]);

    }
    // fprintf(stdout, "z %e xe %e factor %e\n",z,*xe,factor);

    /** -> case z < z_reio_start: helium contribution (tanh of simpler argument) */
    /********Helium reionization is taken into account with a camb-like parametrization***********/
      argument = (preio->reionization_parameters[preio->index_helium_fullreio_redshift] - z)
        /preio->reionization_parameters[preio->index_helium_fullreio_width];
      /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
      *xe += preio->reionization_parameters[preio->index_helium_fullreio_fraction]
        * (tanh(argument)+1.)/2.;
    // fprintf(stdout, "z %e x_e %e xe_before %e factor %elambda %e zp %e qp %e\n", z,*xe,preio->reionization_parameters[preio->index_reio_xe_before],factor,preio->reionization_parameters[preio->index_lambda_douspis_et_al],preio->reionization_parameters[preio->index_zp_douspis_et_al],preio->reionization_parameters[preio->index_Qp_douspis_et_al]);



    }


    return _SUCCESS_;

  }
  if(pth->reio_parametrization == reio_stars_sfr_source_term){
    if(z > preio->reionization_parameters[preio->index_reio_start]) {
      *xe = preio->reionization_parameters[preio->index_reio_xe_before];
    }

    else {

    /** -> case z < z_reio_start: helium contribution (tanh of simpler argument) */
    /********Helium reionization is taken into account with a camb-like parametrization***********/
      argument = (preio->reionization_parameters[preio->index_helium_fullreio_redshift] - z)
        /preio->reionization_parameters[preio->index_helium_fullreio_width];
      /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
      *xe =  (preio->reionization_parameters[preio->index_reio_xe_after]
             -preio->reionization_parameters[preio->index_reio_xe_before])
        *(tanh(argument)+1.)/2.
        +preio->reionization_parameters[preio->index_reio_xe_before]; //both He reionization are done simulatenouesly
      // *xe= MIN(*xe,preio->reionization_parameters[preio->index_reio_xe_after]);
      // fprintf(stderr, "z %e xe %e  xe after %e xe before %e Yhe %e argument %e \n",z, *xe,preio->reionization_parameters[preio->index_reio_xe_after],preio->reionization_parameters[preio->index_reio_xe_before],preio->reionization_parameters[preio->index_helium_fullreio_fraction],(tanh(argument)+1.)/2);
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
      @@ -1412,27 +1907,7 @@ int thermodynamics_reionization_function(
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
                   pth->error_message,
                   pth->error_message);

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
      class_stop(pth->error_message,"Check: is it normal that we are here?");
    }

    else {

      i=0;
      while (preio->reionization_parameters[preio->index_reio_first_z+i+1] < z) i++;

      double z_min = preio->reionization_parameters[preio->index_reio_first_z+i];
      double z_max = preio->reionization_parameters[preio->index_reio_first_z+i+1];

      class_test(z<z_min,
                 pth->error_message,
                 "");

      class_test(z>z_max,
                 pth->error_message,
                 "");

      double x=(z-preio->reionization_parameters[preio->index_reio_first_z+i])
        /(preio->reionization_parameters[preio->index_reio_first_z+i+1]
          -preio->reionization_parameters[preio->index_reio_first_z+i]);

      *xe = preio->reionization_parameters[preio->index_reio_first_xe+i]
        + x*(preio->reionization_parameters[preio->index_reio_first_xe+i+1]
             -preio->reionization_parameters[preio->index_reio_first_xe+i]);

      class_test(*xe<0.,
                 pth->error_message,
                 "%e %e %e\n",
                 x,
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

int thermodynamics_get_xe_before_reionization(
                                              struct precision * ppr,
                                              struct thermo * pth,
                                              struct recombination * preco,
                                              double z,
                                              double * xe
                                              ) {

  int last_index=0;

  class_call(array_interpolate_one_growing_closeby(preco->recombination_table,
                                                   preco->re_size,
                                                   preco->rt_size,
                                                   preco->index_re_z,
                                                   z,
                                                   &last_index,
                                                   preco->index_re_xe,
                                                   xe,
                                                   pth->error_message),
             pth->error_message,
             pth->error_message);

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

int thermodynamics_reionization(
                                struct precision * ppr,
                                struct background * pba,
                                struct thermo * pth,
                                struct recombination * preco,
                                struct reionization * preio,
                                double * pvecback
                                ) {

  /** Summary: */

  /** - define local variables */

  int counter;
  double z_sup,z_mid,z_inf;
  double tau_sup,tau_mid,tau_inf;
  int bin;
  int point;
  double xe_input,xe_actual;

  /** - allocate the vector of parameters defining the function \f$ X_e(z) \f$ */

  class_alloc(preio->reionization_parameters,preio->reio_num_params*sizeof(double),pth->error_message);

  /** (a) if reionization implemented like in CAMB */
  if ((pth->reio_parametrization == reio_camb) || (pth->reio_parametrization == reio_half_tanh) ) {

    /** - --> set values of these parameters, excepted those depending on the reionization redshift */

    if (pth->reio_parametrization == reio_camb ) {
      preio->reionization_parameters[preio->index_reio_xe_after] = 1. + pth->YHe/(_not4_*(1.-pth->YHe));    /* xe_after_reio: H + singly ionized He (note: segmentation fault impossible, checked before that denominator is non-zero) */
    }
    if (pth->reio_parametrization == reio_half_tanh ) {
      /********Modified by Vivian Poulin to take into account helium reionization in both cases***********/
      preio->reionization_parameters[preio->index_reio_xe_after] = 1.;
      // ; /* xe_after_reio: neglect He ionization */
      // + pth->YHe/(_not4_*(1.-pth->YHe));    /* xe_after_reio: H + singly ionized H */
      // + 2*pth->YHe/(_not4_*(1.-pth->YHe));    /* xe_after_reio: H + fully ionized He */
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

    /** - --> if reionization redshift given as an input, initialize the remaining values and fill reionization table*/

    if (pth->reio_z_or_tau == reio_z) {

      /* reionization redshift */
      preio->reionization_parameters[preio->index_reio_redshift] = pth->z_reio;

      /* infer starting redshift for hydrogen */

      if (pth->reio_parametrization == reio_camb || pth->reio_parametrization == reio_half_tanh) {

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

      /* infer xe_before_reio */
      class_call(thermodynamics_get_xe_before_reionization(ppr,
                                                           pth,
                                                           preco,
                                                           preio->reionization_parameters[preio->index_reio_start],
                                                           &(preio->reionization_parameters[preio->index_reio_xe_before])),
                 pth->error_message,
                 pth->error_message);

      /* fill reionization table */
      class_call(thermodynamics_reionization_sample(ppr,pba,pth,preco,preio,pvecback),
                 pth->error_message,
                 pth->error_message);

      pth->tau_reio=preio->reionization_optical_depth;

    }

    /** - --> if reionization optical depth given as an input, find reionization redshift by dichotomy and initialize the remaining values */

    if (pth->reio_z_or_tau == reio_tau) {

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
      class_call(thermodynamics_get_xe_before_reionization(ppr,
                                                           pth,
                                                           preco,
                                                           preio->reionization_parameters[preio->index_reio_start],
                                                           &(preio->reionization_parameters[preio->index_reio_xe_before])),
                 pth->error_message,
                 pth->error_message);

      /* fill reionization table */
      class_call(thermodynamics_reionization_sample(ppr,pba,pth,preco,preio,pvecback),
                 pth->error_message,
                 pth->error_message);

      tau_sup=preio->reionization_optical_depth;

      class_test(tau_sup < pth->tau_reio,
                 pth->error_message,
                 "parameters are such that reionization cannot start after z_start_max");

      /* lower value */

      z_inf = 0.;
      tau_inf = 0.;

      /* try intermediate values */

      counter=0;
      while ((tau_sup-tau_inf) > pth->tau_reio * ppr->reionization_optical_depth_tol) {
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
                   pth->error_message,
                   "starting redshift for reionization > reionization_z_start_max = %e",ppr->reionization_z_start_max);

        /* infer xe_before_reio */
        class_call(thermodynamics_get_xe_before_reionization(ppr,
                                                             pth,
                                                             preco,
                                                             preio->reionization_parameters[preio->index_reio_start],
                                                             &(preio->reionization_parameters[preio->index_reio_xe_before])),
                   pth->error_message,
                   pth->error_message);

        /* clean and fill reionization table */
        free(preio->reionization_table);
        class_call(thermodynamics_reionization_sample(ppr,pba,pth,preco,preio,pvecback),
                   pth->error_message,
                   pth->error_message);

        tau_mid=preio->reionization_optical_depth;

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

    }

    free(preio->reionization_parameters);

    return _SUCCESS_;

  }

  /** - (b) if reionization implemented with reio_bins_tanh scheme */

  if (pth->reio_parametrization == reio_bins_tanh) {

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

    /* infer xe before reio */
    class_call(thermodynamics_get_xe_before_reionization(ppr,
                                                         pth,
                                                         preco,
                                                         preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1],
                                                         &(preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-1])),
               pth->error_message,
               pth->error_message);

    /* infer xe after reio */
    preio->reionization_parameters[preio->index_reio_first_xe] = 1. + pth->YHe/(_not4_*(1.-pth->YHe));    /* xe_after_reio: H + singly ionized He (note: segmentation fault impossible, checked before that denominator is non-zero) */

    /* pass step sharpness parameter */
    preio->reionization_parameters[preio->index_reio_step_sharpness] = pth->binned_reio_step_sharpness;

    /* fill reionization table */
    class_call(thermodynamics_reionization_sample(ppr,pba,pth,preco,preio,pvecback),
               pth->error_message,
               pth->error_message);

    pth->tau_reio=preio->reionization_optical_depth;

    return _SUCCESS_;

  }

  /** - (c) if reionization implemented with reio_many_tanh scheme */

  if (pth->reio_parametrization == reio_many_tanh) {

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

    /* infer xe before reio */
    class_call(thermodynamics_get_xe_before_reionization(ppr,
                                                         pth,
                                                         preco,
                                                         preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1],
                                                         &(preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-1])),
               pth->error_message,
               pth->error_message);

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

    /* fill reionization table */
    class_call(thermodynamics_reionization_sample(ppr,pba,pth,preco,preio,pvecback),
               pth->error_message,
               pth->error_message);

    pth->tau_reio=preio->reionization_optical_depth;

    return _SUCCESS_;

  }
  if((pth->reio_parametrization == reio_douspis_et_al)){
    preio->reionization_parameters[preio->index_reio_xe_after] = 1. + pth->YHe/(_not4_*(1.-pth->YHe));
    preio->reionization_parameters[preio->index_helium_fullreio_fraction] = pth->YHe/(_not4_*(1.-pth->YHe)); /* helium_fullreio_fraction (note: segmentation fault impossible, checked before that denominator is non-zero) */
    preio->reionization_parameters[preio->index_helium_fullreio_redshift] = pth->helium_fullreio_redshift; /* helium_fullreio_redshift */
    preio->reionization_parameters[preio->index_helium_fullreio_width] = pth->helium_fullreio_width;    /* helium_fullreio_width */
    preio->reionization_parameters[preio->index_lambda_douspis_et_al] = pth->lambda_douspis_et_al;
    preio->reionization_parameters[preio->index_zp_douspis_et_al] = pth->zp_douspis_et_al;
    preio->reionization_parameters[preio->index_Qp_douspis_et_al] = pth->Qp_douspis_et_al;

    class_test(preio->reionization_parameters[preio->index_helium_fullreio_width]==0,
               pth->error_message,
               "stop to avoid division by zero");
     preio->reionization_parameters[preio->index_reio_start] = ppr->reionization_z_start_max;

     class_test(preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max,
                pth->error_message,
                "starting redshift for reionization > reionization_z_start_max = %e\n",ppr->reionization_z_start_max);

     /* infer xe_before_reio */
     class_call(thermodynamics_get_xe_before_reionization(ppr,
                                                          pth,
                                                          preco,
                                                          preio->reionization_parameters[preio->index_reio_start],
                                                          &(preio->reionization_parameters[preio->index_reio_xe_before])),
                pth->error_message,
                pth->error_message);
     /* fill reionization table */
     class_call(thermodynamics_reionization_sample(ppr,pba,pth,preco,preio,pvecback),
                pth->error_message,
                pth->error_message);
     pth->tau_reio=preio->reionization_optical_depth;

     return _SUCCESS_;

  }

  if((pth->reio_parametrization == reio_asymmetric_planck_16)){
    preio->reionization_parameters[preio->index_reio_xe_after] = 1. + pth->YHe/(_not4_*(1.-pth->YHe));
    preio->reionization_parameters[preio->index_helium_fullreio_fraction] = pth->YHe/(_not4_*(1.-pth->YHe)); /* helium_fullreio_fraction (note: segmentation fault impossible, checked before that denominator is non-zero) */
    preio->reionization_parameters[preio->index_helium_fullreio_redshift] = pth->helium_fullreio_redshift; /* helium_fullreio_redshift */
    preio->reionization_parameters[preio->index_helium_fullreio_width] = pth->helium_fullreio_width;    /* helium_fullreio_width */
    preio->reionization_parameters[preio->index_alpha_asymmetric_planck_16] = pth->alpha_asymmetric_planck_16;
    preio->reionization_parameters[preio->index_z_end_asymmetric_planck_16] = pth->z_end_asymmetric_planck_16;
    preio->reionization_parameters[preio->index_reio_start] = pth->z_start_asymmetric_planck_16;
    pth->z_10_percent = 0;
    pth->z_50_percent = 0;
    pth->z_99_percent = 0;
    class_test(preio->reionization_parameters[preio->index_helium_fullreio_width]==0,
               pth->error_message,
               "stop to avoid division by zero");
    //  preio->reionization_parameters[preio->index_reio_start] = ppr->reionization_z_start_max;

     class_test(preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max,
                pth->error_message,
                "starting redshift for reionization > reionization_z_start_max = %e\n",ppr->reionization_z_start_max);

     /* infer xe_before_reio */
     class_call(thermodynamics_get_xe_before_reionization(ppr,
                                                          pth,
                                                          preco,
                                                          preio->reionization_parameters[preio->index_reio_start],
                                                          &(preio->reionization_parameters[preio->index_reio_xe_before])),
                pth->error_message,
                pth->error_message);
     /* fill reionization table */
     class_call(thermodynamics_reionization_sample(ppr,pba,pth,preco,preio,pvecback),
                pth->error_message,
                pth->error_message);
     pth->tau_reio=preio->reionization_optical_depth;

     return _SUCCESS_;

  }

  if((pth->reio_parametrization == reio_stars_sfr_source_term)){//Helium reionization is still a tanh, to be improved
    preio->reionization_parameters[preio->index_reio_xe_after] = 1. + 2*pth->YHe/(_not4_*(1.-pth->YHe));
    preio->reionization_parameters[preio->index_helium_fullreio_fraction] = pth->YHe/(_not4_*(1.-pth->YHe)); /* helium_fullreio_fraction (note: segmentation fault impossible, checked before that denominator is non-zero) */
    preio->reionization_parameters[preio->index_helium_fullreio_redshift] = pth->helium_fullreio_redshift; /* helium_fullreio_redshift */
    preio->reionization_parameters[preio->index_helium_fullreio_width] = pth->helium_fullreio_width;    /* helium_fullreio_width */
    preio->reionization_parameters[preio->index_reio_start] = pth->helium_fullreio_redshift+pth->helium_fullreio_width;
    class_test(preio->reionization_parameters[preio->index_helium_fullreio_width]==0,
               pth->error_message,
               "stop to avoid division by zero");
    //  preio->reionization_parameters[preio->index_reio_start] = ppr->reionization_z_start_max;
     class_test(preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max,
                pth->error_message,
                "starting redshift for reionization > reionization_z_start_max = %e\n",ppr->reionization_z_start_max);

     /* infer xe_before_reio */
     class_call(thermodynamics_get_xe_before_reionization(ppr,
                                                          pth,
                                                          preco,
                                                          preio->reionization_parameters[preio->index_reio_start],
                                                          &(preio->reionization_parameters[preio->index_reio_xe_before])),
                pth->error_message,
                pth->error_message);
     /* fill reionization table */
     class_call(thermodynamics_reionization_sample(ppr,pba,pth,preco,preio,pvecback),
                pth->error_message,
                pth->error_message);
     pth->tau_reio=preio->reionization_optical_depth;

     return _SUCCESS_;

  }
  /** - (d) if reionization implemented with reio_inter scheme */

  if (pth->reio_parametrization == reio_inter) {

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

    /* infer xe before reio */
    class_call(thermodynamics_get_xe_before_reionization(ppr,
                                                         pth,
                                                         preco,
                                                         preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1],
                                                         &(preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-1])),
               pth->error_message,
               pth->error_message);

    /* fill reionization table */
    class_call(thermodynamics_reionization_sample(ppr,pba,pth,preco,preio,pvecback),
               pth->error_message,
               pth->error_message);

    pth->tau_reio=preio->reionization_optical_depth;

    return _SUCCESS_;

  }
  class_test(0 == 0,
             pth->error_message,
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

int thermodynamics_reionization_sample(
                                       struct precision * ppr,
                                       struct background * pba,
                                       struct thermo * pth,
                                       struct recombination * preco,
                                       struct reionization * preio,
                                       double * pvecback
                                       ) {

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
  int i,j;
  int number_of_redshifts;
  /* values of z, dz, X_e */
  double dz,dz_max;
  double z,z_next;
  double xe,xe_next,x_tmp;
  double dkappadz,dkappadz_next;
  double delta_z_old, delta_z_new;
  double Tb,Yp,dTdz,dTdz_adia,dTdz_CMB,dTdz_DM,dTdz_stars,opacity,mu;
  double dkappadtau,dkappadtau_next;
  double energy_rate;
  double tau;
  double chi_heat, chi_heat_x_ray;
  double chi_lya;
  double chi_ionH;
  double chi_ionHe;
  double chi_lowE;
  double argument;
  int last_index_back;
  double relative_variation;
  double L_x, rho_sfr;
  Yp = pth->YHe;

  /** - (a) allocate vector of values related to reionization */
  class_alloc(reio_vector,preio->re_size*sizeof(double),pth->error_message);

  /** - (b) create a growTable with gt_init() */
  class_call(gt_init(&gTable),
             gTable.error_message,
             pth->error_message);

  /** - (c) first line is taken from thermodynamics table, just before reionization starts */

  /** - --> look where to start in current thermodynamics table */
  i=0;
  while (preco->recombination_table[i*preco->re_size+preco->index_re_z] < preio->reionization_parameters[preio->index_reio_start]) {
    i++;
    class_test(i == ppr->recfast_Nz0,
               pth->error_message,
               "reionization_z_start_max = %e > largest redshift in thermodynamics table",ppr->reionization_z_start_max);
  }
  if(preco->recombination_table[i*preco->re_size+preco->index_re_z] >  preio->reionization_parameters[preio->index_reio_start])i--;
  j=i;
  /** - --> get redshift */
  z=preco->recombination_table[i*preco->re_size+preco->index_re_z];

  reio_vector[preio->index_re_z]=z;
  preio->index_reco_when_reio_start=i;

  /** - --> get \f$ X_e \f$ */
  class_call(thermodynamics_reionization_function(z,pth,preio,preco,&xe),
             pth->error_message,
             pth->error_message);

    if(pth->reio_stars_and_dark_matter == _TRUE_){
      xe=preco->recombination_table[i*preco->re_size+preco->index_re_xe];
      // xe=MAX(xe,x_tmp);
    }
  reio_vector[preio->index_re_xe] = xe;

  /** -  --> get \f$ d \kappa / d z = (d \kappa / d \tau) * (d \tau / d z) = - (d \kappa / d \tau) / H \f$ */

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

  reio_vector[preio->index_re_dkappadtau] = (1.+z) * (1.+z) * pth->n_e * xe * _sigma_ * _Mpc_over_m_;

  class_test(pvecback[pba->index_bg_H] == 0.,
             pth->error_message,
             "stop to avoid division by zero");

  reio_vector[preio->index_re_dkappadz] = reio_vector[preio->index_re_dkappadtau] / pvecback[pba->index_bg_H];

  dkappadz = reio_vector[preio->index_re_dkappadz];
  dkappadtau = reio_vector[preio->index_re_dkappadtau];

  /** - --> get baryon temperature **/
  Tb = preco->recombination_table[i*preco->re_size+preco->index_re_Tb];
  reio_vector[preio->index_re_Tb] = Tb;

  /** - --> after recombination, Tb scales like (1+z)**2. Compute constant factor Tb/(1+z)**2. */
  //Tba2 = Tb/(1+z)/(1+z);

  /** - --> get baryon sound speed */
  reio_vector[preio->index_re_cb2] = 5./3. * _k_B_ / ( _c_ * _c_ * _m_H_ ) * (1. + (1./_not4_ - 1.) * Yp + xe * (1.-Yp)) * Tb;

  /** - --> store these values in growing table */
  class_call(gt_add(&gTable,_GT_END_,(void *) reio_vector,sizeof(double)*(preio->re_size)),
             gTable.error_message,
             pth->error_message);

  number_of_redshifts=1;

  /** - (d) set the maximum step value (equal to the step in thermodynamics table) */
  dz_max=preco->recombination_table[i*preco->re_size+preco->index_re_z]
    -preco->recombination_table[(i-1)*preco->re_size+preco->index_re_z];
    // fprintf(stderr, "dz %e\n", dz_max);
  /** - (e) loop over redshift values in order to find values of z, x_e, kappa' (Tb and cb2 found later by integration). The sampling in z space is found here. */

  /* initial step */
  dz = dz_max;

  while (z > 0.) {
    if (j<0)j=0;
    // fprintf(stdout, "j %d \n",j);
    // dz = MAX(ppr->smallest_allowed_variation,dz);

    class_test(dz < ppr->smallest_allowed_variation,
               pth->error_message,
               "stuck in the loop for reionization sampling, as if you were trying to impose a discontinuous evolution for xe(z)");

    /* - try next step */
    z_next=z-dz;
    delta_z_old = z_next-preco->recombination_table[(j-1)*preco->re_size+preco->index_re_z];
    delta_z_new = z_next-preco->recombination_table[(j-2)*preco->re_size+preco->index_re_z];
    if(fabs(delta_z_old)<fabs(delta_z_new))j++;
    while(z_next > preco->recombination_table[(j-2)*preco->re_size+preco->index_re_z])j++;
    // fprintf(stdout, "z = %e z_next = %e\n",z,z_next);

    if (z_next < 0.) z_next=0.;
    class_call(thermodynamics_reionization_function(z_next,pth,preio,preco,&xe_next),
               pth->error_message,
               pth->error_message);

    if(pth->reio_stars_and_dark_matter == _TRUE_){
      /**
       * This small routine compares the reionization table to the recombination one and choose the highest x_e between the two.
       * This way allows to enables to avoid unphysical discontiniuty in the ionization fraction at low x_e.
       * First, we interpolate the value of x_e at the evaluated redshift from the recombination table. The linear interpolation has been checked
       * to work well but it could be improve for security. Then we perform the comparison.
       */

      x_tmp= (preco->recombination_table[(j-2)*preco->re_size+preco->index_re_xe]-preco->recombination_table[(j-1)*preco->re_size+preco->index_re_xe])/(preco->recombination_table[(j-2)*preco->re_size+preco->index_re_z]
        -preco->recombination_table[(j-1)*preco->re_size+preco->index_re_z])*(z_next-preco->recombination_table[(j-1)*preco->re_size+preco->index_re_z])+
        preco->recombination_table[(j-1)*preco->re_size+preco->index_re_xe]  ;
      x_tmp = MAX(0.,x_tmp); // Small check to avoid negative values of x_e.

      if(x_tmp <1. + 2.*pth->YHe/(_not4_*(1.-pth->YHe))) xe_next=MAX(xe_next,x_tmp); // Here the comparison is made.
      else x_tmp = 1. + 2.*pth->YHe/(_not4_*(1.-pth->YHe)); // the maximal value that x_e can reach.

    }

    class_call(background_tau_of_z(pba,
                                   z_next,
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
    // if(xe_next > 1.17) fprintf(stdout, "error xe next %e\n", xe_next);
    dkappadz_next= (1.+z_next) * (1.+z_next) * pth->n_e * xe_next * _sigma_ * _Mpc_over_m_ / pvecback[pba->index_bg_H];
    dkappadtau_next= (1.+z_next) * (1.+z_next) * pth->n_e * xe_next * _sigma_ * _Mpc_over_m_;


    class_test((dkappadz == 0.) || (dkappadtau == 0.),
               pth->error_message,
               "stop to avoid division by zero");

    relative_variation = fabs((dkappadz_next-dkappadz)/dkappadz) +
      fabs((dkappadtau_next-dkappadtau)/dkappadtau);

    if (relative_variation < ppr->reionization_sampling || pth->reio_stars_and_dark_matter == _TRUE_) {
      /* accept the step: get \f$ z, X_e, d kappa / d z \f$ and store in growing table */

      z=z_next;
      xe=xe_next;
      dkappadz=dkappadz_next;
      dkappadtau= dkappadtau_next;
      class_test((dkappadz == 0.) || (dkappadtau == 0.),
                 pth->error_message,
                 "dkappadz=%e, dkappadtau=%e, stop to avoid division by zero",dkappadz,dkappadtau);

      reio_vector[preio->index_re_z] = z;
      reio_vector[preio->index_re_xe] = xe;
      reio_vector[preio->index_re_dkappadz] = dkappadz;
      reio_vector[preio->index_re_dkappadtau] = dkappadz * pvecback[pba->index_bg_H];

      class_call(gt_add(&gTable,_GT_END_,(void *) reio_vector,sizeof(double)*(preio->re_size)),
                 gTable.error_message,
                 pth->error_message);

      number_of_redshifts++;
      // delta_z_old = z_next-preco->recombination_table[(j-1)*preco->re_size+preco->index_re_z];
      // delta_z_new = z_next-preco->recombination_table[(j-2)*preco->re_size+preco->index_re_z];
      // if(fabs(delta_z_old)>fabs(delta_z_new))j--;
      j--;
      dz = MIN(0.9*(ppr->reionization_sampling/relative_variation),5.)*dz;
      // dz = MIN(dz,dz_max);
      // dz = MAX(ppr->smallest_allowed_variation,dz);
    }
    else {
      /* do not accept the step and update dz */
      // delta_z_old = z_next-preco->recombination_table[(j-1)*preco->re_size+preco->index_re_z];
      // delta_z_new = z_next-preco->recombination_table[(j-2)*preco->re_size+preco->index_re_z];
      dz = 0.9*(ppr->reionization_sampling/relative_variation)*dz;
      // dz = MIN(dz,dz_max);
      // dz = MAX(ppr->smallest_allowed_variation,dz);
      // j--;
      // if(fabs(delta_z_old)>fabs(delta_z_new))j--;

    }
  }

  /** - (f) allocate reionization_table with correct size */
  class_alloc(preio->reionization_table,preio->re_size*number_of_redshifts*sizeof(double),pth->error_message);

  preio->rt_size=number_of_redshifts;

  /** - (g) retrieve data stored in the growTable with gt_getPtr() */
  class_call(gt_getPtr(&gTable,(void**)&pData),
             gTable.error_message,
             pth->error_message);

  /** - (h) copy growTable to reionization_temporary_table (invert order of lines, so that redshift is growing, like in recombination table) */
  for (i=0; i < preio->rt_size; i++) {
    memcopy_result = memcpy(preio->reionization_table+i*preio->re_size,pData+(preio->rt_size-i-1)*preio->re_size,preio->re_size*sizeof(double));
    class_test(memcopy_result != preio->reionization_table+i*preio->re_size,
               pth->error_message,
               "cannot copy data back to reionization_temporary_table");

  }

  /** - (i) free the growTable with gt_free() , free vector of reionization variables */
  class_call(gt_free(&gTable),
             gTable.error_message,
             pth->error_message);

  free(reio_vector);

  /** - (j) another loop on z, to integrate equation for Tb and to compute cb2 */
  for (i=preio->rt_size-1; i >0 ; i--) {

    z = preio->reionization_table[i*preio->re_size+preio->index_re_z];

    class_call(background_tau_of_z(pba,
                                   z,
                                   &tau),
               pba->error_message,
               pth->error_message);

    class_call(background_at_tau(pba,
                                 tau,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &last_index_back,
                                 pvecback),
               pba->error_message,
               pth->error_message);

    dz = (preio->reionization_table[i*preio->re_size+preio->index_re_z]-preio->reionization_table[(i-1)*preio->re_size+preio->index_re_z]);

    opacity = (1.+z) * (1.+z) * pth->n_e
      * preio->reionization_table[i*preio->re_size+preio->index_re_xe] * _sigma_ * _Mpc_over_m_;

    mu = _m_H_/(1. + (1./_not4_ - 1.) * pth->YHe + preio->reionization_table[i*preio->re_size+preio->index_re_xe] * (1.-pth->YHe));


    /** - derivative of baryon temperature */

      /** - First possibility: Add a tanh term in the temperature and bypass evolution equation*/
      if(pth->star_heating_parametrization== heating_reiolike_tanh){

        argument = (pow((1.+preio->reionization_parameters[preio->index_reio_redshift]),
                        preio->reionization_parameters[preio->index_reio_exponent])
                    - pow((1.+z),preio->reionization_parameters[preio->index_reio_exponent]))
          /(preio->reionization_parameters[preio->index_reio_exponent]
            /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
            *pow((1.+preio->reionization_parameters[preio->index_reio_redshift]),
                 (preio->reionization_parameters[preio->index_reio_exponent]-1.)))
          /preio->reionization_parameters[preio->index_reio_width];
        /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
          if(z< preio->reionization_parameters[preio->index_reio_start])
          preio->reionization_table[(i-1)*preio->re_size+preio->index_re_Tb] = (pth->final_IGM_temperature
                 -preio->reionization_table[i*preio->re_size+preio->index_re_Tb])
        *(tanh(argument)+1.)/2 + preio->reionization_table[i*preio->re_size+preio->index_re_Tb];

      }
      /** - Second possibility: Compute temperature evolution from each sources*/

      else {


    dTdz_adia=2./(1+z)*preio->reionization_table[i*preio->re_size+preio->index_re_Tb];

    dTdz_CMB = - 2.*mu/_m_e_*4.*pvecback[pba->index_bg_rho_g]/3./pvecback[pba->index_bg_rho_b]*opacity*
      (pba->T_cmb * (1.+z)-preio->reionization_table[i*preio->re_size+preio->index_re_Tb])/pvecback[pba->index_bg_H];

      /** - Parameters related to exotic energy injection */
      if((pth->annihilation != 0 || pth->decay_fraction != 0 || pth->PBH_high_mass != 0 || pth->PBH_low_mass != 0)){

            /** - --> derivative of baryon temperature */
              preco->xe_tmp=preio->reionization_table[i*preio->re_size+preio->index_re_xe];
              preco->Tm_tmp=preio->reionization_table[i*preio->re_size+preio->index_re_Tb];

              class_call(thermodynamics_energy_injection(ppr,pba,preco,z,&energy_rate,pth->error_message),
                         pth->error_message,
                         pth->error_message);

              preco->z_tmp=z;
               /* coefficient as revised by Slatyer et al. 2013 (in fact it is an interpolation by Vivian Poulin of columns 1 and 2 in Table V of Slatyer et al. 2013) */
               if(pth->energy_repart_functions==Galli_et_al_interpolation){
                 class_call(thermodynamics_annihilation_coefficients_interpolate(ppr,pba,pth,preio->reionization_table[i*preio->re_size+preio->index_re_xe]),
                          pth->error_message,
                          pth->error_message);
                 chi_heat = pth->chi_heat;
              }
               if(pth->energy_repart_functions==no_factorization){
                 class_call(thermodynamics_annihilation_coefficients_interpolate(ppr,pba,pth,z),
                          pth->error_message,
                          pth->error_message);
                 chi_heat = pth->chi_heat;
              }
              /* coefficient as revised by Slatyer et al. 2013 (in fact it is an fit by Vivian Poulin of columns 1 and 2 in Table V of Slatyer et al. 2013) */
               if(pth->energy_repart_functions==Galli_et_al_fit){
                 chi_heat = 0.996857*(1.-pow(1.-pow(preio->reionization_table[i*preio->re_size+preio->index_re_xe],0.300134),1.51035));
               }
               /* old approximation from Chen and Kamionkowski */
               if(pth->energy_repart_functions==SSCK){
                 chi_heat = (1.+2.*preio->reionization_table[i*preio->re_size+preio->index_re_xe])/3.;
               }

              chi_heat= MIN(chi_heat,1.);
              chi_heat = MAX(chi_heat,0.);


              dTdz_DM = - 2./(3.*_k_B_)*energy_rate*chi_heat
              /(preco->Nnow*pow(1.+z,3))/(1.+preco->fHe+preio->reionization_table[i*preio->re_size+preio->index_re_xe])
              /(pvecback[pba->index_bg_H]*_c_/_Mpc_over_m_*(1.+z)); /* energy injection */

              if(pth->thermodynamics_verbose>10){
                fprintf(stdout, "z %e dTdz_CMB %edTdz_DM %e xe %e energy_rate %e chi_heat%e\n",z,dTdz_CMB,dTdz_DM,preio->reionization_table[i*preio->re_size+preio->index_re_xe],energy_rate, chi_heat);
              }

      }
      else dTdz_DM = 0.;


      /** parametrization of reheating by stars */
      if(pth->star_heating_parametrization == heating_none){ //Standard assumption, no reheating by stars. Enough for accurate computations of CMB power spectra.
        dTdz_stars = 0;
      }

      else if(pth->star_heating_parametrization== heating_stars_sfr_source_term ){ //Reheating term based on the SFR rates. See Poulin et al. 1508.01370.
          rho_sfr = pth->ap*pow(1+z,pth->bp)/(1+pow((1+z)/pth->cp,pth->dp))/pow(_Mpc_over_m_,3)*pow(1+z,3)*(1+tanh((pth->z_start_reio_stars-z)))/2; //add a (sharp) smoothing function.
          L_x = pth->Ex* pth->fx *rho_sfr*2./(3.*_k_B_*preco->Nnow*pow(1.+z,3)*(1.+preco->fHe+preio->reionization_table[i*preio->re_size+preio->index_re_xe]))
          /(pvecback[pba->index_bg_H]*_c_/_Mpc_over_m_*(1.+z));

          dTdz_stars = -L_x*(1+2*preio->reionization_table[i*preio->re_size+preio->index_re_xe])/3.;
        }

      // else if(){
      //   // ready for other parametrization
      // }

      dTdz = dTdz_adia+dTdz_CMB+dTdz_DM+dTdz_stars;
      if(pth->thermodynamics_verbose>10){
      fprintf(stdout, "z %e dT %e Tmat %e dTdz_adia %e dTdz_CMB %e dTdz_DM %e dTdz_stars %e \n", z,dTdz, preio->reionization_table[i*preio->re_size+preio->index_re_Tb],dTdz_adia, dTdz_CMB ,dTdz_DM,dTdz_stars);
      }
      /** - --> increment baryon temperature  */

        preio->reionization_table[(i-1)*preio->re_size+preio->index_re_Tb] =
          preio->reionization_table[i*preio->re_size+preio->index_re_Tb]-dTdz*dz;
      }


    /** - get baryon sound speed */

    preio->reionization_table[(i-1)*preio->re_size+preio->index_re_cb2] = _k_B_/ ( _c_ * _c_ * mu)
      * preio->reionization_table[(i-1)*preio->re_size+preio->index_re_Tb]
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
                          pth->error_message),
             pth->error_message,
             pth->error_message);

  /** - --> integrate for optical depth */
  class_call(array_integrate_all_spline(preio->reionization_table,
                                        preio->re_size,
                                        preio->rt_size,
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
 * Integrate thermodynamics with your favorite recombination code.
 *
 */

int thermodynamics_recombination(
                                 struct precision * ppr,
                                 struct background * pba,
                                 struct thermo * pth,
                                 struct recombination * preco,
                                 double * pvecback
                                 ) {

  if (pth->recombination==hyrec) {

    class_call(thermodynamics_recombination_with_hyrec(ppr,pba,pth,preco,pvecback),
               pth->error_message,
               pth->error_message);

  }

  if (pth->recombination==recfast) {

    class_call(thermodynamics_recombination_with_recfast(ppr,pba,pth,preco,pvecback),
               pth->error_message,
               pth->error_message);

  }

  if (pth->recombination==cosmorec) {
        class_call(thermodynamics_recombination_with_cosmorec(ppr,pba,pth,preco,pvecback),
               pth->error_message,
               pth->error_message);
  }

  return _SUCCESS_;

}
/**
 * Integrate thermodynamics with CosmoRec.
 *
 * Integrate thermodynamics with CosmoRec, allocate and fill the part
 * of the thermodynamics interpolation table (the rest is filled in
 * thermodynamics_init()). Called once by
 * thermodynamics_recombination(), from thermodynamics_init().
 *
 *************************************************************************************************
 *                 CosmoRec: Cosmological Recombination Project
 *                Written by Jens Chluba (University of Manchester)
 *************************************************************************************************
 * @param ppr      Input: pointer to precision structure
 * @param pba      Input: pointer to background structure
 * @param pth      Input: pointer to thermodynamics structure
 * @param preco    Output: pointer to recombination structure
 * @param pvecback Input: pointer to an allocated (but empty) vector of background variables
 */
int thermodynamics_recombination_with_cosmorec(
                                            struct precision * ppr,
                                            struct background * pba,
                                            struct thermo * pth,
                                            struct recombination * preco,
                                            double * pvecback
                                            ) {
#ifdef COSMOREC
  int i;
  double nH0 = 11.223846333047*pba->Omega0_b*pba->h*pba->h*(1.-pth->YHe);  /* number density of hydrogen today in m-3 */
  double rho_cdm_today = pow(pba->H0*_c_/_Mpc_over_m_,2)*3/8./_PI_/_G_*pba->Omega0_cdm*_c_*_c_; /* energy density in J/m^3 */
  double DM_annihilation =  pth->annihilation*1e-6/_c_/_c_*pow(rho_cdm_today,2)/nH0*1e6/_eV_; /*conversion in cosmorec unit as described in Chluba 2010 0910.3663 (without factor 2, to respect class convention of majorana particles)*/
  double runpars[4] = {
    DM_annihilation, /* defines the dark matter annihilation efficiency in eV/s. */
    pth->cosmorec_accuracy, /* setting for cosmorec accuracy (default = default cosmorec setting) */
    pth->cosmorec_verbose, /* setting for cosmorec verbose (default = no output produced) */
    pth->Lambda_over_theoritical_Lambda *_Lambda_ /* theoritical value by Labzowsky et al 2005 for H1_A2s_1s is rescaled, by default Lambda_over_theoritical_Lambda = 1. In agreement with standard cosmorec.*/
  };

  double H0 = pba->H0 / 1e3 * _c_;
  int nz = ppr->recfast_Nz0;
  double * z_arr;
  double * Hz_arr;
  double z, xe, Tm, Hz;
  double z_start=ppr->recfast_z_initial;
  double z_end=0;
  double step;

  double tau_at_z;
  int last_index;

  double * xe_out;
  double * tb_out;

  double drho_dt = 0, Tg;
  double dlnTb_dz;
  int label=0; /* iterator for cosmorec output file name, not used in class version of cosmorec */

  /* Initialize Hubble rate for CosmoRec */
  class_alloc(z_arr, sizeof(double) * nz, pth->error_message);
  class_alloc(Hz_arr, sizeof(double) * nz, pth->error_message);

  step = (z_start - z_end) / (nz);
  for(i=0; i < nz; i++) {
    z_arr[i] = z_end + i * step;

      class_call(
        background_tau_of_z(
          pba,
          z_arr[i],
          &tau_at_z
        ),
        pba->error_message,
        pth->error_message
      );

      class_call(
        background_at_tau(
          pba,
          tau_at_z,
          pba->short_info,
          pba->inter_normal,
          &last_index,
          pvecback
        ),
        pba->error_message,
        pth->error_message
      );

      Hz_arr[i]=pvecback[pba->index_bg_H] * _c_ / _Mpc_over_m_;
  }
  /* Initialize x_e and tb output tables */

  class_alloc(xe_out, sizeof(double) * nz, pth->error_message);
  class_alloc(tb_out, sizeof(double) * nz, pth->error_message);

  /* call cosmorec */
  /* Currently we give parameters separetely, eventually to be changed for a structure, easier to modify.*/

  cosmorec_calc_h_cpp_(
    &(pth->cosmorec_runmode), runpars,
    &(pba->Omega0_cdm), &(pba->Omega0_b), &(pba->Omega0_k),
    &(pba->Neff), &H0,
    &(pba->T_cmb), &(pth->YHe),
    z_arr, Hz_arr, &nz,
    z_arr, xe_out, tb_out,
    &nz,
    &label
  );

  /** - fill a few parameters in preco and pth */



  preco->rt_size = nz;
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
  preco->decay_fraction = pth->decay_fraction;
  preco->annihilation_f_halo = pth->annihilation_f_halo;
  preco->annihilation_z_halo = pth->annihilation_z_halo;
  pth->n_e=preco->Nnow;

  /** - allocate memory for thermodynamics interpolation tables (size known in advance) and fill it */

  class_alloc(preco->recombination_table,preco->re_size*preco->rt_size*sizeof(double),pth->error_message);

  for(i=nz-1; i >= 0; i--) {

    /** - --> get redshift, corresponding results from cosmorec, and background quantities */

    z = z_arr[i];
    xe = xe_out[i];
    Tm = tb_out[i];
    Hz = Hz_arr[i];
    /** - --> store the results in the table */

    /* results are obtained in order of decreasing z, and stored in order of growing z */

    /* redshift */
    *(preco->recombination_table+(i)*preco->re_size+preco->index_re_z)=z;

    /* ionization fraction */
    *(preco->recombination_table+(i)*preco->re_size+preco->index_re_xe)=xe;

    /* Tb */
    *(preco->recombination_table+(i)*preco->re_size+preco->index_re_Tb)=Tm;

    /* cb2 = (k_B/mu) Tb (1-1/3 dlnTb/dlna) = (k_B/mu) Tb (1+1/3 (1+z) dlnTb/dz)
       with (1+z)dlnTb/dz= - [dlnTb/dlna] */
    /* note that m_H / mu = 1 + (m_H/m_He-1) Y_p + x_e (1-Y_p) */

    Tg = pba->T_cmb * (1+z);
    dlnTb_dz = - Tg/Tm*drho_dt/(1+z)/Hz+1/(1+z);

   evaluate_TM(z, xe,preco->fHe, Tm/Tg, Tg, Hz, &drho_dt);
    *(preco->recombination_table+(i)*preco->re_size+preco->index_re_cb2)
      = _k_B_ / ( _c_ * _c_ * _m_H_ ) * (1. + (1./_not4_ - 1.) * pth->YHe + xe * (1.-pth->YHe)) * Tm * (1. + (1+z)*dlnTb_dz / 3.);
    /* dkappa/dtau = a n_e x_e sigma_T = a^{-2} n_e(today) x_e sigma_T (in units of 1/Mpc) */
    *(preco->recombination_table+(i)*preco->re_size+preco->index_re_dkappadtau)
      = (1.+z) * (1.+z) * preco->Nnow * xe * _sigma_ * _Mpc_over_m_;
      //  fprintf(stdout,"xe %e Tm %e cb2 %e z %e dlnTb_dz %e *dkappa_dtau %e\n",xe,Tm,*(preco->recombination_table+(i)*preco->re_size+preco->index_re_cb2),z,dlnTb_dz,*(preco->recombination_table+(i)*preco->re_size+preco->index_re_dkappadtau));

  }

  /* clean up */

  free(xe_out);
  free(tb_out);

  free(z_arr);
  free(Hz_arr);

#else

class_stop(pth->error_message,
           "you compiled without including the CosmoRec code, and now wish to use it. Either set the input parameter 'recombination' to something else than 'CosmoRec', or recompile after setting in the Makefile the appropriate path COSMOREC=... ");


#endif /* COSMOREC */

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
 int thermodynamics_recombination_with_hyrec(
                                             struct precision * ppr,
                                             struct background * pba,
                                             struct thermo * pth,
                                             struct recombination * preco,
                                             double * pvecback
                                             ) {
   /** Summary: */
 #ifdef HYREC

   HYREC_DATA hyrec_data;
   hyrec_allocate(&hyrec_data, ppr->recfast_z_initial, 0.);

   double Omega_m = pba->Omega0_b + pba->Omega0_cdm + pba->Omega0_ncdm_tot;

   double alpha_ratio = 1.;    /* Ratio of fine-structure constant to standard value */
   double me_ratio    = 1.;    /* Ratio of electron mass to standard value */

   double pann        = 1.78266e-21 *pth->annihilation;  /* Converting from m^3/s/kg to cm^3/s/GeV */
   double pann_halo   = 1.78266e-21 *pth->annihilation_f_halo;

   int i,j,Nz;
   double z, xe, Tm, Hz;
   void * buffer;
   double tau;
   int last_index_back;
   int on_the_spot = 1;

   if(pth->has_on_the_spot == _FALSE_){
     on_the_spot = 0;
   }
   /** - Compute the recombination history by calling hyrec_compute.
         No CLASS-like error management here, but YAH working on it :) **/

   if (pth->thermodynamics_verbose > 0)
     printf(" -> calling HyRec version %s,\n",HYREC_VERSION);

   hyrec_compute(&hyrec_data, FULL,
 		pba->h, pba->T_cmb, pba->Omega0_b, Omega_m, pba->Omega0_k, pth->YHe, pba->Neff,
 		alpha_ratio, me_ratio, pann, pann_halo, pth->annihilation_z, pth->annihilation_zmax,
 		pth->annihilation_zmin, pth->annihilation_variation, pth->annihilation_z_halo,
 		pth->PBH_high_mass, pth->PBH_fraction, pth->coll_ion_pbh,on_the_spot);

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
   preco->decay_fraction = pth->decay_fraction;
   preco->annihilation_f_halo = pth->annihilation_f_halo;
   preco->annihilation_z_halo = pth->annihilation_z_halo;
   pth->n_e=preco->Nnow;

   /** - allocate memory for thermodynamics interpolation tables (size known in advance) and fill it */

   class_alloc(preco->recombination_table,preco->re_size*preco->rt_size*sizeof(double),pth->error_message);

   for(i = 0; i < Nz; i++) {

     /** - --> get redshift, corresponding results from hyrec, and background quantities */

     z = ppr->recfast_z_initial * (1. - (double)(i+1) / (double)Nz);

     xe = hyrec_xe(z, &hyrec_data);
     Tm = hyrec_Tm(z, &hyrec_data);

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
       = _k_B_ / ( _c_ * _c_ * _m_H_ ) * (1. + (1./_not4_ - 1.) * pth->YHe + xe * (1.-pth->YHe)) * Tm *(1. - hyrec_dTmdlna(z, &hyrec_data) / Tm / 3.);

     /* dkappa/dtau = a n_e x_e sigma_T = a^{-2} n_e(today) x_e sigma_T (in units of 1/Mpc) */
     *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_dkappadtau)
       = (1.+z) * (1.+z) * preco->Nnow * xe * _sigma_ * _Mpc_over_m_;

   }

   /* Cleanup */

   free(buffer);
   hyrec_free(&hyrec_data);

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

int thermodynamics_recombination_with_recfast(
                                              struct precision * ppr,
                                              struct background * pba,
                                              struct thermo * pth,
                                              struct recombination * preco,
                                              double * pvecback
                                              ) {

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
  struct thermodynamics_parameters_and_workspace tpaw;

  /** - allocate memory for thermodynamics interpolation tables (size known in advance) */
  preco->rt_size = ppr->recfast_Nz0;
  class_alloc(preco->recombination_table,preco->re_size*preco->rt_size*sizeof(double),pth->error_message);

  /** - initialize generic integrator with initialize_generic_integrator() */
  class_call(initialize_generic_integrator(_RECFAST_INTEG_SIZE_, &gi),
             gi.error_message,
             pth->error_message);

  /** - read a few precision/cosmological parameters */

  /* Nz */
  Nz=ppr->recfast_Nz0;

  /* preco->H0 is H0 in inverse seconds (while pba->H0 is [H0/c] in inverse Mpcs) */
  preco->H0 = pba->H0 * _c_ / _Mpc_over_m_;

  /* Omega_b */
  OmegaB = pba->Omega0_b;

  /* Yp */
  preco->YHe = pth->YHe;

  /* Tnow */
  preco->Tnow = pba->T_cmb;

  /* z_initial */
  zinitial=ppr->recfast_z_initial;

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
  z=zinitial;
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
  preco->decay_fraction = pth->decay_fraction;
  preco->PBH_high_mass = pth->PBH_high_mass;
  preco->PBH_disk_formation_redshift = pth->PBH_disk_formation_redshift;
  preco->PBH_accretion_recipe = pth->PBH_accretion_recipe;
  preco->energy_deposition_treatment = pth->energy_deposition_treatment;
  preco->PBH_low_mass = pth->PBH_low_mass;
  preco->PBH_fraction = pth->PBH_fraction;
  preco->energy_repart_functions = pth->energy_repart_functions;
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
  tpaw.pba = pba;
  tpaw.ppr = ppr;
  tpaw.preco = preco;
  tpaw.pth = pth;
  tpaw.pvecback = pvecback;

  /** - impose initial conditions at early times */

  class_test(zinitial < ppr->recfast_z_He_3,
             pth->error_message,
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

    else if (y[1] > ppr->recfast_x_He0_trigger ) {
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

    else if (y[0] > ppr->recfast_x_H0_trigger && z > 200) {

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
                 pth->error_message);

      y[0] = MIN(x_H0,1);   //Vivian

      if (pth->thermodynamics_verbose > 1) {
        fprintf(stdout, "in function thermodynamics_recombination_with_recfast, fifth approximation : zend %e y[0] %e\n",zend, y[0]);
      }
      /* smoothed transition */
      if (ppr->recfast_x_He0_trigger - y[1] < ppr->recfast_x_He0_trigger_delta && z > 200) {
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
      // x0 = y[0] + preco->fHe*y[1];

    }

    /** - --> last case: full evolution for H and Helium */

    else {

      /* quantities used for smoothed transition */
      if (ppr->recfast_x_H0_trigger - y[0] < ppr->recfast_x_H0_trigger_delta  && z > 200 ) {
        rhs = exp(1.5*log(preco->CR*preco->Tnow/(1.+z)) - preco->CB1/(preco->Tnow*(1.+z)))/preco->Nnow;
        x_H0 = 0.5*(sqrt(pow(rhs,2)+4.*rhs) - rhs);
      }
      else x_H0 = y[0];

      class_call(generic_integrator(thermodynamics_derivs_with_recfast,
                                    zstart,
                                    zend,
                                    y,
                                    &tpaw,
                                    ppr->tol_thermo_integration,
                                    ppr->smallest_allowed_variation,
                                    &gi),
                 gi.error_message,
                 pth->error_message);
       y[0] = MIN(y[0],1);//Vivian


      /* smoothed transition */
      if (ppr->recfast_x_H0_trigger - y[0] < ppr->recfast_x_H0_trigger_delta && z > 200) {
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

        // x0 = y[0] + preco->fHe*y[1];
        if(pth->thermodynamics_verbose>1){
          fprintf(stdout, "in function thermodynamics_recombination_with_recfast, full calculation zend %e x0 %e y[0] %e\n",zend, x0, y[0]);
        }

    }
  /*  double argument;
    argument = (pth->helium_fullreio_redshift - z)
        /pth->helium_fullreio_width;
      x0 += pth->YHe/(_not4_*(1.-pth->YHe))* (tanh(argument)+1.);*/

    /** - --> store the results in the table */
    /* results are obtained in order of decreasing z, and stored in order of growing z */

    /* redshift */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_z)=zend;

    /* ionization fraction */

    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_xe)=x0;

    /* Tb */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_Tb)=y[2];

    /* get dTb/dz=dy[2] */
    class_call(thermodynamics_derivs_with_recfast(zend, y, dy, &tpaw,pth->error_message),
               pth->error_message,
               pth->error_message);

    /* cb2 = (k_B/mu) Tb (1-1/3 dlnTb/dlna) = (k_B/mu) Tb (1+1/3 (1+z) dlnTb/dz) */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_cb2)
      = _k_B_ / ( _c_ * _c_ * _m_H_ ) * (1. + (1./_not4_ - 1.) * preco->YHe + x0 * (1.-preco->YHe)) * y[2] * (1. + (1.+zend) * dy[2] / y[2] / 3.);

    /* dkappa/dtau = a n_e x_e sigma_T = a^{-2} n_e(today) x_e sigma_T (in units of 1/Mpc) */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_dkappadtau)
      = (1.+zend) * (1.+zend) * preco->Nnow * x0 * _sigma_ * _Mpc_over_m_;
      if(pth->thermodynamics_verbose>1){
        fprintf(stdout,"%e %e %e %e %e %e\n",
             *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_z),
             *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_xe),
             *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_Tb),
             (1.+zend) * dy[2],
             *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_cb2),
             *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_dkappadtau)
             );
      }


  }

  /** - cleanup generic integrator with cleanup_generic_integrator() */

  class_call(cleanup_generic_integrator(&gi),
             gi.error_message,
             pth->error_message);

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
                                       double z,
                                       double * y,
                                       double * dy,
                                       void * parameters_and_workspace,
                                       ErrorMsg error_message
                                       ) {


  /* define local variables */
  double x,n,n_He,Trad,Tmat,x_H,x_He,Hz,dHdz,epsilon;
  double Rup,Rup_2,Rdown,K,K_He,Rup_He,Rup_He_2,Rdown_He,He_Boltz;
  double timeTh,timeH;
  double sq_0,sq_1;

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
  double * pvecback;
  int last_index_back;

  /* used for energy injection from dark matter */
  double C;
  //double C_He;
  double energy_rate;

  double tau;
  double chi_heat;
  double chi_lya;
  double chi_ionH;
  double chi_ionHe;
  double chi_lowE;
  double dTdz_DM, dTdz_CMB, dTdz_adia, dTdz_stars;
   /*used for reionization from realistic star model*/
  double rho_sfr,stars_xe,dNion_over_dt,L_x;

  ptpaw = parameters_and_workspace;
  ppr = ptpaw->ppr;
  pba = ptpaw->pba;
  pth = ptpaw->pth;
  preco = ptpaw->preco;
  pvecback = ptpaw->pvecback;
  rho_sfr = pth->ap*pow(1+z,pth->bp)/(1+pow((1+z)/pth->cp,pth->dp))/pow(_Mpc_over_m_,3)*pow(1+z,3)*(1+tanh((pth->z_start_reio_stars-z)))/2;//add a (sharp) smoothing function.

  /* security added by Vivian Poulin to avoid bug when dealing with energy injection modifying ionisation history */
  x_H = MIN(y[0],1.);
  x_H = MAX(y[0],0.);
  x_He = MIN(y[1],1.);
  x_He = MAX(y[1],0.);
  x = MIN(x_H + preco->fHe * x_He,1+preco->fHe);
  x = MAX(x_H + preco->fHe * x_He,0.);

  // x_H = y[0];
  // x_He = y[1];
  // x = x_H + preco->fHe * x_He;

  Tmat = MAX(y[2],0.);

  // fprintf(stderr, "input xH %e xHe %e x %e Tmat %e\n",x_H, x_He, x, Tmat );
  n = preco->Nnow * (1.+z) * (1.+z) * (1.+z);
  n_He = preco->fHe * n;
  Trad = preco->Tnow * (1.+z);

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

   if((pth->annihilation!=0 || pth->decay_fraction!=0 || pth->PBH_high_mass!=0 || pth->PBH_low_mass != 0)){
     preco->xe_tmp=x;
     preco->Tm_tmp=Tmat;
  // if( z > 2){//sometimes problem with interpolation
    // fprintf(stdout, "z %e,Tmat %e, x %e\n",z,Tmat,x);
  class_call(thermodynamics_energy_injection(ppr,pba,preco,z,&energy_rate,error_message),
             error_message,
             error_message);
  // fprintf(stdout, "energy_rate %e\n",energy_rate);
  // }
  // else energy_rate = 0;
     preco->z_tmp=z;
}
else energy_rate=0;
 // fprintf(stdout,"%e      %e     %e      %e      %e    \n", x,pth->chi_heat,pth->chi_lya, pth->chi_ionH,pth->chi_ionHe,pth->chi_lowE);
  /* Hz is H in inverse seconds (while pvecback returns [H0/c] in inverse Mpcs) */
  Hz=pvecback[pba->index_bg_H]* _c_ / _Mpc_over_m_;

  Rdown=1.e-19*_a_PPB_*pow((Tmat/1.e4),_b_PPB_)/(1.+_c_PPB_*pow((Tmat/1.e4),_d_PPB_));
  Rup_2 = 1.e-19*_a_PPB_*pow((Trad/1.e4),_b_PPB_)/(1.+_c_PPB_*pow((Trad/1.e4),_d_PPB_)) * pow((preco->CR*Trad),1.5)*exp(-preco->CDB/Trad);
  Rup = Rdown * pow((preco->CR*Tmat),1.5)*exp(-preco->CDB/Tmat);

  sq_0 = sqrt(Tmat/_T_0_);
  sq_1 = sqrt(Tmat/_T_1_);
  Rdown_He = _a_VF_/(sq_0 * pow((1.+sq_0),(1.-_b_VF_)) * pow((1. + sq_1),(1. + _b_VF_)));
  Rup_He_2 = 4.*Rdown_He*pow((preco->CR*Trad),1.5)*exp(-preco->CDB_He/Trad);
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
    /* equations modified to take into account energy injection from dark matter */
      chi_ionH = 0.;
      chi_ionHe = 0.;
      chi_lya = 0.;

    if(preco->annihilation > 0 || preco->decay_fraction > 0 || preco->PBH_high_mass > 0 || preco->PBH_low_mass > 0){
      if (x < 1.){
        /* coefficient as revised by Galli et al. 2013 (in fact it is an interpolation by Vivian Poulin of Table V of Galli et al. 2013) */
        if(pth->energy_repart_functions==Galli_et_al_interpolation){
          class_call(thermodynamics_annihilation_coefficients_interpolate(ppr,pba,pth,x),
                         error_message,
                         error_message);
          chi_ionH = pth->chi_ionH;
          chi_ionHe = pth->chi_ionHe;
          chi_lya = pth->chi_lya;

        }
        if(pth->energy_repart_functions==no_factorization){
          class_call(thermodynamics_annihilation_coefficients_interpolate(ppr,pba,pth,z),
                         error_message,
                         error_message);
          chi_ionH = pth->chi_ionH;
          chi_ionHe = pth->chi_ionHe;
          chi_lya = pth->chi_lya;
        }
        /* old approximation from Chen and Kamionkowski */
        if(pth->energy_repart_functions==SSCK){
          chi_ionH = (1.-x)/3.;
          chi_lya = chi_ionH;
          chi_ionHe=0;
        }
        /* coefficient as revised by Slatyer et al. 2013 (in fact it is a fit by Vivian Poulin of columns 1 and 2 in Table V of Slatyer et al. 2013): */
        if(pth->energy_repart_functions==Galli_et_al_fit){
          chi_ionH = 0.369202*pow(1.-pow(x,0.463929),1.70237);
          chi_ionHe =0.0312604*pow(1.-pow(x,0.200634),0.82247);
          chi_lya = 0.335597*pow(1.-pow(x,0.375314),1.80722);
        }

        chi_ionH = MIN(chi_ionH,1.);
        chi_ionHe = MIN(chi_ionHe,1.);
        chi_lya = MIN(chi_lya,1.);
        chi_ionH = MAX(chi_ionH,0.);
        chi_ionHe = MAX(chi_ionHe,0.);
        chi_lya = MAX(chi_lya,0.);

      }
      else {
        chi_ionH = 0.;
        chi_ionHe = 0.;
        chi_lya = 0.;
      }
      if(chi_ionH < 0 )fprintf(stdout, "chi_ionH %e \n",chi_ionH);
      if(pth->thermodynamics_verbose>10){
        fprintf(stdout, "chi_ionH %e chi_ionHe %e chi_lya %e z% e\n", chi_ionH , chi_ionHe, chi_lya, z);
      }
    }

    /* Peebles' coefficient (approximated as one when the Hydrogen
           ionization fraction is very close to one) */
    if (x_H < ppr->recfast_x_H0_trigger2) {
      C = (1. + K*pth->Lambda_over_theoritical_Lambda*_Lambda_*n*(1.-x_H))/(1./preco->fu+K*pth->Lambda_over_theoritical_Lambda*_Lambda_*n*(1.-x_H)/preco->fu +K*Rup_2*n*(1.-x_H));  /* 2 modifications : 1) Rup -> Rup_2 evaluating the coefficient using Trad instead of Tmat; 2) add pth->Lambda_over_theoritical_Lambda, 1 in the standard case, allow to constraint A2s1s otherwise*/
      // C = (1. + K*_Lambda_*n*(1.-x_H))/(1./preco->fu+K*_Lambda_*n*(1.-x_H)/preco->fu +K*Rup*n*(1.-x_H));
      // fprintf(stdout, "A2s1s %e\n", pth->Lambda_over_theoritical_Lambda*_Lambda_);
    }
    else {
      C = 1.;

      }
      /* evolution of hydrogen ionisation fraction: */

      dy[0] = (x*x_H*n*Rdown - Rup_2*(1.-x_H)*exp(-preco->CL/Tmat)) * C / (Hz*(1.+z))       /* Peeble's equation with fudged factors */
            -energy_rate/n*((chi_ionH+chi_ionHe)/_L_H_ion_+chi_lya*(1.-C)/_L_H_alpha_)/(_h_P_*_c_*Hz*(1.+z)); /* energy injection (neglect fraction going to helium) */

      // dy[0] = -5.89e-5*sqrt(Tmat/1e4)*exp(-1.58e5/Tmat)*(1-x_H)*x/ (Hz*(1.+z)) * 1e-6;     // Collisional ionisation, taken from 1503.04827, last factor is for conversion cm^3->m^3


      // fprintf(stdout, "z %e Tmat %e collision %e DM  %e standard %e\n",z, Tmat, -5.89e-5*sqrt(Tmat/1e4)*exp(-1.58e5/Tmat)*(1-x_H)*x/ (Hz*(1.+z)) * 1e-6,-energy_rate/n*((chi_ionH+chi_ionHe)/_L_H_ion_+chi_lya*(1.-C)/_L_H_alpha_)/(_h_P_*_c_*Hz*(1.+z)),(x*x_H*n*Rdown - Rup_2*(1.-x_H)*exp(-preco->CL/Tmat)) * C / (Hz*(1.+z)));
      if(pth->reio_parametrization == reio_stars_sfr_source_term){
        dNion_over_dt=pth->f_esc*pth->Zeta_ion*rho_sfr;
        stars_xe=dNion_over_dt/(Hz*(1.+z)*n);
        dy[0] -= stars_xe*(1-x)/3;
        // fprintf(stdout, " %e  %e  %e %e %e %e  %e\n",rho_sfr,stars_xe, dNion_over_dt,Hz,n,(1-x)/3,z );
      }
    // JL: test for debugginf reio_inter
    //fprintf(stdout,"%e  %e  %e  %e\n",z,Tmat,K*_Lambda_*n,K*Rup*n);

      if(pth->thermodynamics_verbose>10){
      fprintf(stdout, "z %e Tmat %e  DM  %e standard %e stars %e \n",z, Tmat,-energy_rate/n*((chi_ionH+chi_ionHe)/_L_H_ion_+chi_lya*(1.-C)/_L_H_alpha_)/(_h_P_*_c_*Hz*(1.+z)),(x*x_H*n*Rdown - Rup_2*(1.-x_H)*exp(-preco->CL/Tmat)) * C / (Hz*(1.+z)),stars_xe);
      }
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

    dy[1] = ((x*x_He*n*Rdown_He - Rup_He_2*(1.-x_He)*exp(-preco->CL_He/Tmat))
             *(1. + K_He*_Lambda_He_*n_He*(1.-x_He)*He_Boltz))
      /(Hz*(1+z)* (1. + K_He*(_Lambda_He_+Rup_He_2)*n_He*(1.-x_He)*He_Boltz)); /* in case of energy injection due to DM, we neglect the contribution to helium ionization */

    // dy[1] += -2.02e-9*sqrt(Tmat/1e4)*exp(-2.85e5/Tmat)*(1-x_He)*x / (Hz*(1.+z)) * 1e-6;     // Collisional ionisation, taken from 1503.04827, last factor is for conversion cm^3->m^3
      //
      // /*******************Helium**********************/
      // dxedlna+=stars_xe*param->fHe*(1+tanh((6-z1)/0.5));
      // if(z1<6)dxedlna+=stars_xe*param->fHe*(1+tanh((3.5-z1)/0.5));
      /***********************************************/
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
    dHdz=-pvecback[pba->index_bg_H_prime]/pvecback[pba->index_bg_H]/pba->a_today* _c_ / _Mpc_over_m_;
    epsilon = Hz * (1.+x+preco->fHe) / (preco->CT*pow(Trad,3)*x);
    dy[2] = preco->Tnow + epsilon*((1.+preco->fHe)/(1.+preco->fHe+x))*((dy[0]+preco->fHe*dy[1])/x)
      - epsilon* dHdz/Hz + 3.*epsilon/(1.+z) ;
  }
  else {
    /* equations modified to take into account energy injection from dark matter */

    if(pth->annihilation >0 || pth->decay_fraction > 0 || pth->PBH_high_mass > 0 || pth->PBH_low_mass > 0){
      if (x < 1.){
        /* coefficient as revised by Galli et al. 2013 (in fact it is an interpolation by Vivian Poulin of columns 1 and 2 in Table V of Galli et al. 2013) */
        if(pth->energy_repart_functions==Galli_et_al_interpolation){
          class_call(thermodynamics_annihilation_coefficients_interpolate(ppr,pba,pth,x),
                        error_message,
                        error_message);
            chi_heat = pth->chi_heat;
        }
        if(pth->energy_repart_functions==no_factorization){
          class_call(thermodynamics_annihilation_coefficients_interpolate(ppr,pba,pth,z),
                        error_message,
                        error_message);
            chi_heat = pth->chi_heat;
        }
        /* coefficient as revised by Slatyer et al. 2013 (in fact it is a fit by Vivian Poulin of columns 1 and 2 in Table V of Slatyer et al. 2013) */
        if(pth->energy_repart_functions==Galli_et_al_fit){
            chi_heat = 0.996857*(1.-pow(1.-pow(x,0.300134),1.51035));
        }
        /* old approximation from Chen and Kamionkowski */
        if(pth->energy_repart_functions==SSCK){
            chi_heat = (1.+2.*x)/3.;
        }

      }
      else
        chi_heat = 1.;

        chi_heat= MIN(chi_heat,1.);
        chi_heat = MAX(chi_heat,0.);
      // fprintf(stdout, "z %e chi_heat %e chi_heat_old %e xe %e\n",z ,chi_heat,(1.+2.*x)/3., x);
    }
    else chi_heat = 0.;
    dTdz_adia=2.*Tmat/(1.+z);
    dTdz_CMB = preco->CT * pow(Trad,4) * x / (1.+x+preco->fHe) * (Tmat-Trad) / (Hz*(1.+z));
    dTdz_DM = -2./(3.*_k_B_)*energy_rate*chi_heat/n/(1.+preco->fHe+x)/(Hz*(1.+z));
    if(pth->star_heating_parametrization == heating_stars_sfr_source_term){
    L_x = 2*pth->Ex * pth->fx * rho_sfr/(3*_k_B_*n*Hz*(1.+z)*(1.+x+preco->fHe));
    dTdz_stars = -L_x*(1+2*x)/3.;
    }
    else dTdz_stars = 0;
    dy[2]= dTdz_CMB + dTdz_adia + dTdz_DM + dTdz_stars; /* dTdz_DM = energy injection */

    if(pth->thermodynamics_verbose>10){
      fprintf(stdout, "chi_heat %e z% e\n", chi_heat, z);
      fprintf(stdout, "z %e dT %e Tmat %e Trad %e dTdz_adia %e dTdz_CMB %e dTdz_DM %e dTdz_stars %e x %e\n", z, dy[2], Tmat, Trad,dTdz_adia, dTdz_CMB ,dTdz_DM,dTdz_stars,x);
      }

    // dy[2] += 5.89e-5*sqrt(Tmat/1e4)*exp(-1.58e5/Tmat)*(1-x_H)*x/ (Hz*(1.+z)) * 1e-6;
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

  if ((pth->reio_parametrization != reio_none)) {
    class_test(preco->recombination_table[preio->index_reco_when_reio_start*preco->re_size+preco->index_re_z] !=
               preio->reionization_table[(preio->rt_size -1)*preio->re_size+preio->index_re_z],
               pth->error_message,
               "mismatch which should never happen");
  }

  /** - find number of redshift in full table = number in reco + number in reio - overlap */

  pth->tt_size = ppr->recfast_Nz0 + preio->rt_size - preio->index_reco_when_reio_start - 1;


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
  for (i=0; i < ppr->recfast_Nz0 - preio->index_reco_when_reio_start - 1; i++) {
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

  free(preco->recombination_table);
  if(pth->has_on_the_spot == _FALSE_ && pth->energy_repart_functions!=no_factorization){
    thermodynamics_annihilation_f_eff_free(preco);
  }
  if(pth->energy_repart_functions == Galli_et_al_interpolation || pth->energy_repart_functions==no_factorization){
    thermodynamics_annihilation_coefficients_free(pth);
  }
  if ((pth->reio_parametrization != reio_none))
    free(preio->reionization_table);

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
  //class_store_columntitle(titles,"max. rate",_TRUE_,colnum);
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

int thermodynamics_tanh(double x,
                        double center,
                        double before,
                        double after,
                        double width,
                        double * result) {

  *result = before + (after-before)*(tanh((x-center)/width)+1.)/2.;

  return _SUCCESS_;
}
