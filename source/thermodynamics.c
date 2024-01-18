/** @file thermodynamics.c Documented thermodynamics module
 *
 * * Julien Lesgourgues, 6.09.2010
 * * Restructured by Nils Schoeneberg and Matteo Lucca, 27.02.2019
 * * Evolver implementation by Daniel Meinert, spring 2019
 *
 * Deals with the thermodynamical evolution.
 * This module has two purposes:
 *
 * - at the beginning, to initialize the thermodynamics, i.e. to
 *   integrate the thermodynamical equations, and store all
 *   thermodynamical quantities as a function of redshift inside an
 *   interpolation table.
 *
 * - to provide a routine which allow other modules to evaluate any
 *   thermodynamical quantities at a given redshift value (by
 *   interpolating within the interpolation table).
 *
 * The most important differential equations to compute the free
 * electron fraction x at each step are provided either by the HyRec
 * 2020 or RecFastCLASS code, located in the external/ directory. The
 * thermodynamics module integrates these equations using the generic
 * integrator (which can be set to ndf15, rkck4, etc.) The HyRec and
 * RecFastCLASS algorithms are used and called in the same way by this
 * module.
 *
 * In summary, the following functions can be called from other modules:
 *
 * -# thermodynamics_init at the beginning (but after background_init)
 * -# thermodynamics_at_z at any later time
 * -# thermodynamics_free at the end, when no more calls to thermodynamics_at_z are needed
 */

#include "thermodynamics.h"

#include "history.h"
#include "hyrectools.h"
#include "helium.h"
#include "wrap_hyrec.h"


/**
 * Thermodynamics quantities at given redshift z.
 * Evaluates all thermodynamics quantities at a given value of the redshift by reading the pre-computed table and interpolating.
 *
 * @param pba          Input: pointer to background structure
 * @param pth          Input: pointer to the thermodynamics structure (containing pre-computed table)
 * @param z            Input: redshift
 * @param inter_mode   Input: interpolation mode (normal or growing_closeby)
 * @param last_index   Input/Output: index of the previous/current point in the interpolation array (input only for closeby mode, output for both)
 * @param pvecback     Input: vector of background quantities (used only in case z>z_initial for getting ddkappa and dddkappa; in that case,
                            should be already allocated and filled, with format short_info or larger; in other cases, will be ignored)
 * @param pvecthermo Output: vector of thermodynamics quantities (assumed to be already allocated)
 * @return the error status
 */

int thermodynamics_at_z(
                        struct background * pba,
                        struct thermodynamics * pth,
                        double z,
                        enum interpolation_method inter_mode,
                        int * last_index,
                        double * pvecback,
                        double * pvecthermo
                        ) {

  /** Summary: */

  /** - define local variables */
  double x0;

  /* The fact that z is in the pre-computed range 0 <= z <= z_initial will be checked in the interpolation routines below. Before
     trying to interpolate, allow the routine to deal with the case z > z_initial: then, all relevant quantities can be extrapolated
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

    /* \f$ exp^{-\kappa}, g, g', g'' \f$ can be set to zero: they are used only for computing the source functions in the
       perturbation module; but source functions only need to be sampled below z_initial (the condition that
       z_start_sources<z_initial is checked in the perturbation module) */
    pvecthermo[pth->index_th_exp_m_kappa] = 0.;
    pvecthermo[pth->index_th_g]=0.;
    pvecthermo[pth->index_th_dg]=0.;
    pvecthermo[pth->index_th_ddg]=0.;

    /* Calculate Tb */
    pvecthermo[pth->index_th_Tb] = pba->T_cmb*(1.+z);

    /* Tb derivative */
    pvecthermo[pth->index_th_dTb] = pba->T_cmb;

    /* Calculate baryon equation of state parameter wb = (k_B/mu) Tb */
    /* note that m_H / mu = 1 + (m_H/m_He-1) Y_p + x_e (1-Y_p) */
    pvecthermo[pth->index_th_wb] = _k_B_ / ( _c_ * _c_ * _m_H_ ) * (1. + (1./_not4_ - 1.) * pth->YHe + x0 * (1.-pth->YHe)) * pba->T_cmb * (1.+z);

    /* Calculate cb2 (cb2 = (k_B/mu) Tb (1-1/3 dlnTb/dlna) = (k_B/mu) Tb (1+1/3 (1+z) dlnTb/dz)) */
    /* note that m_H / mu = 1 + (m_H/m_He-1) Y_p + x_e (1-Y_p) */
    pvecthermo[pth->index_th_cb2] = pvecthermo[pth->index_th_wb] * 4. / 3.;

    /* derivatives of baryon sound speed (only computed if some non-minimal tight-coupling schemes is requested) */
    if (pth->compute_cb2_derivatives == _TRUE_) {

      /* since cb2 proportional to (1+z) or 1/a, its derivative wrt conformal time is given by dcb2 = - a H cb2 */
      pvecthermo[pth->index_th_dcb2] = - pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a] * pvecthermo[pth->index_th_cb2];

      /* then its second derivative is given by ddcb2 = - a H' cb2 */
      pvecthermo[pth->index_th_ddcb2] = - pvecback[pba->index_bg_H_prime] * pvecback[pba->index_bg_a] * pvecthermo[pth->index_th_cb2];
    }

    /* in this regime, variation rate = dkappa/dtau */
    pvecthermo[pth->index_th_rate] = pvecthermo[pth->index_th_dkappa];

     /* quantities related to DM interacting with DR */
    if (pba->has_idm_dr == _TRUE_) {

      /* calculate dmu_idm_dr and approximate its derivatives as zero */
      pvecthermo[pth->index_th_dmu_idm_dr] = pth->a_idm_dr*pow((1.+z)/1.e7,pth->nindex_idm_dr)*pba->Omega0_idm_dr*pow(pba->h,2);
      pvecthermo[pth->index_th_ddmu_idm_dr] =  -pvecback[pba->index_bg_H] * pth->nindex_idm_dr / (1+z) * pvecthermo[pth->index_th_dmu_idm_dr];
      pvecthermo[pth->index_th_dddmu_idm_dr] = (pvecback[pba->index_bg_H]*pvecback[pba->index_bg_H]/ (1.+z) - pvecback[pba->index_bg_H_prime])
        *  pth->nindex_idm_dr / (1.+z) * pvecthermo[pth->index_th_dmu_idm_dr];

      /* calculate dmu_idr (self interaction) */
      pvecthermo[pth->index_th_dmu_idr] = pth->b_idr*pow((1.+z)/1.e7,pth->nindex_idm_dr)*pba->Omega0_idr*pow(pba->h,2);

      /* extrapolate optical depth of idm_dr and idr */
      pvecthermo[pth->index_th_tau_idm_dr] = pth->thermodynamics_table[(pth->tt_size-1)*pth->th_size+pth->index_th_tau_idm_dr]+
        (pth->thermodynamics_table[(pth->tt_size-1)*pth->th_size+pth->index_th_tau_idm_dr]-pth->thermodynamics_table[(pth->tt_size-2)*pth->th_size+pth->index_th_tau_idm_dr])
        *(z-pth->z_table[pth->tt_size-1])/(pth->z_table[pth->tt_size-1]-pth->z_table[pth->tt_size-2]);

      pvecthermo[pth->index_th_tau_idr] = pth->thermodynamics_table[(pth->tt_size-1)*pth->th_size+pth->index_th_tau_idr]+
        (pth->thermodynamics_table[(pth->tt_size-1)*pth->th_size+pth->index_th_tau_idr]-pth->thermodynamics_table[(pth->tt_size-2)*pth->th_size+pth->index_th_tau_idr])
        *(z-pth->z_table[pth->tt_size-1])/(pth->z_table[pth->tt_size-1]-pth->z_table[pth->tt_size-2]);

      /* extrapolate idm_dr visibility function */
      pvecthermo[pth->index_th_g_idm_dr] = pth->thermodynamics_table[(pth->tt_size-1)*pth->th_size+pth->index_th_g_idm_dr];

      /* calculate interacting dark matter sound speed */
      pvecthermo[pth->index_th_cidm_dr2] = 4*_k_B_*pba->T_idr*(1.+z)/_eV_/3./pth->m_idm_dr;

      /* calculate interacting dark matter temperature (equal to idr temperature at this redhsift) */
      pvecthermo[pth->index_th_Tidm_dr] = pba->T_idr*(1.+z);
    }

  }

  /** - interpolate in table with array_interpolate_spline (normal mode) or array_interpolate_spline_growing_closeby (closeby mode) */

  else {

    /* some very specific cases require linear interpolation because of a break in the derivative of the functions */
    if (((pth->reio_parametrization == reio_half_tanh) && (z < 2*pth->z_reio))
        || ((pth->reio_parametrization == reio_inter) && (z < 50.))) {

      class_call(array_interpolate_linear(pth->z_table,
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

      if (inter_mode == inter_normal) {

        class_call(array_interpolate_spline(pth->z_table,
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

      if (inter_mode == inter_closeby) {

        class_call(array_interpolate_spline_growing_closeby(pth->z_table,
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
 * Initialize the thermodynamics structure, and in particular the
 * thermodynamics interpolation table.
 *
 * @param ppr   Input: pointer to precision structure
 * @param pba   Input: pointer to background structure
 * @param pth   Input/Output: pointer to initialized thermodynamics structure
 * @return the error status
 */

int thermodynamics_init(
                        struct precision * ppr,
                        struct background * pba,
                        struct thermodynamics * pth
                        ) {

  /** Summary: */

  /** - define local variables */

  /* vector of background values for calling background_at_tau */
  double * pvecback;

  /* structures for storing temporarily information on recombination and reionization */
  struct thermo_workspace * ptw;

  if (pth->thermodynamics_verbose > 0) {
    switch (pth->recombination) {

    case recfast:
      printf("Computing thermodynamics using RecFastCLASS (based on v1.5)\n");
      break;

    case hyrec:
      printf("Computing thermodynamics using HyRec 2020\n");
      break;

    default:
      class_stop(pth->error_message,"pth->recombination=%d different from all known cases",pth->recombination);
      break;
    }
  }

  /** - compute and check primordial Helium mass fraction rho_He/(rho_H+rho_He) */

  if (pth->YHe == _YHE_BBN_) {
    class_call(thermodynamics_helium_from_bbn(ppr,pba,pth),
               pth->error_message,
               pth->error_message);
  }
  if (pth->thermodynamics_verbose > 0) {
    printf(" -> with primordial helium mass fraction Y_He = %.4f\n",pth->YHe);
  }

  /** - infer primordial helium-to-hydrogen nucleon ratio n_He/n_H
   * It is calculated via n_He/n_H = rho_He/(m_He/m_H * rho_H) = YHe * rho_b / (m_He/m_H * (1-YHe) rho_b) = YHe / (m_He/m_H * (1-YHe))*/
  pth->fHe = pth->YHe/(_not4_ *(1.-pth->YHe));

  /** - infer number of hydrogen nuclei today in m**-3 */
  pth->n_e = 3.*pow(pba->H0 * _c_ / _Mpc_over_m_,2)*pba->Omega0_b/(8.*_PI_*_G_*_m_H_)*(1.-pth->YHe);

  /** -  If there is idm-dr, we want the thermodynamics table to
   * start at a much larger z, in order to capture the possible
   * non-trivial behavior of the dark matter interaction rate at early times*/
  if (pba->has_idm_dr == _TRUE_) {
    ppr->thermo_z_initial = ppr->thermo_z_initial_if_idm_dr;
    ppr->thermo_Nz_log = ppr->thermo_Nz_log_if_idm_dr;
  }

  /** - test whether all parameters are in the correct regime */
  class_call(thermodynamics_checks(ppr,pba,pth),
             pth->error_message,
             pth->error_message);

  /** - allocate and assign all temporary structures and indices */
  class_alloc(ptw, sizeof(struct thermo_workspace), pth->error_message);
  class_call(thermodynamics_workspace_init(ppr,pba,pth,ptw),
             pth->error_message,
             pth->error_message);

  class_call(thermodynamics_indices(pba,pth,ptw),
             pth->error_message,
             pth->error_message);

  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);

  class_call(thermodynamics_lists(ppr,pba,pth,ptw),
             pth->error_message,
             pth->error_message);

  /** - initialize injection struct (not temporary) */
  if (pth->has_exotic_injection == _TRUE_) {
    class_call(injection_init(ppr,
                              pba,
                              pth),
               (pth->in).error_message,
               pth->error_message);
  }

  /** - assign reionisation parameters */
  class_call(thermodynamics_set_parameters_reionization(ppr,
                                                        pba,
                                                        pth,
                                                        ptw->ptrp),
             pth->error_message,
             pth->error_message);

  /** - solve recombination and reionization and store values of \f$ z, x_e, d \kappa / d \tau, T_b, c_b^2 \f$  */
  class_call(thermodynamics_solve(ppr,pba,pth,ptw,pvecback),
             pth->error_message,
             pth->error_message);

  /** - the differential equation system is now completely solved  */

  /** - fill missing columns (quantities not computed during the differential evolution but related) */
  class_call(thermodynamics_calculate_remaining_quantities(ppr,pba,pth,pvecback),
             pth->error_message,
             pth->error_message);

  /** - write information on thermal history in standard output */
  if (pth->thermodynamics_verbose > 0) {
    class_call(thermodynamics_output_summary(pba,pth),
               pth->error_message,
               pth->error_message);
  }

  /** - free workspace and local variables */
  class_call(thermodynamics_workspace_free(pth,ptw),
             pth->error_message,
             pth->error_message);

  free(pvecback);

  return _SUCCESS_;
}


/**
 * Free all memory space allocated by thermodynamics_init.
 *
 * @param pth Input/Output: pointer to thermodynamics structure (to be freed)
 * @return the error status
 */
int thermodynamics_free(
                        struct thermodynamics * pth
                        ) {

  if (pth->has_exotic_injection == _TRUE_) {
    /* Free all injection-related functions */
    class_call(injection_free(pth),
               (pth->in).error_message,
               pth->error_message);
  }

  /* Free thermodynamics-related functions */
  free(pth->z_table);
  free(pth->tau_table);
  free(pth->thermodynamics_table);
  free(pth->d2thermodynamics_dz2_table);

  return _SUCCESS_;
}

/**
 * Infer the primordial helium mass fraction from standard BBN
 * calculations, as a function of the baryon density and expansion
 * rate during BBN.
 *
 * This module is simpler then the one used in arXiv:0712.2826 because
 * it neglects the impact of a possible significant chemical
 * potentials for electron neutrinos. The full code with xi_nu_e could
 * be introduced here later.
 *
 * @param ppr   Input: pointer to precision structure
 * @param pba   Input: pointer to background structure
 * @param pth   Input/Output: pointer to initialized thermodynamics structure
 * @return the error status
 */
int thermodynamics_helium_from_bbn(
                                   struct precision * ppr,
                                   struct background * pba,
                                   struct thermodynamics * pth
                                   ) {

  /** Summary: */

  /** Define local variables */
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
  double Neff_bbn, z_bbn, * pvecback;

  /** - Infer effective number of neutrinos at the time of BBN */
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);

  /** - We randomly choose 0.1 MeV to be the temperature of BBN */
  z_bbn = 0.1*1e6/(_eV_over_Kelvin_*pba->T_cmb)-1.0;

  class_call(background_at_z(pba,
                             z_bbn,
                             long_info,
                             inter_normal,
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

  /** - compute Delta N_eff as defined in bbn file, i.e. \f$ \Delta N_{eff}=0\f$ means \f$ N_{eff}=3.046\f$. Note that even if 3.044 is a better default value, we must keep 3.046 here as long as the BBN file we are using has been computed assuming 3.046. */
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

    /* check that the line is neither blank neither a comment. In ASCII, left[0]>39 means that first non-blank character might
       be the beginning of some data (it is not a newline, a #, a %, etc.) */
    if (left[0] > 39) {

      /* if the line contains data, we must interpret it. If (num_omegab, num_deltaN)=(0,0), the current line must contain
         their values. Otherwise, it must contain (omegab, delatN, YHe). */
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
      else{

        /* read (omegab, deltaN, YHe) */
        class_test(sscanf(line,"%lg %lg %lg",&(omegab[array_line%num_omegab]),
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
 * Check the thermodynamics structure parameters for bounds and critical values.
 *
 * @param ppr   Input: pointer to precision structure
 * @param pba   Input: pointer to background structure
 * @param pth   Input: pointer to initialized thermodynamics structure
 * @return the error status
 */

int thermodynamics_checks(
                          struct precision * ppr,
                          struct background* pba,
                          struct thermodynamics * pth
                          ) {

  /** Summary: */

  /** - check BBN Y_He fracion */
  class_test((pth->YHe < _YHE_SMALL_)||(pth->YHe > _YHE_BIG_),
             pth->error_message,
             "Y_He=%g out of bounds (%g<Y_He<%g)",pth->YHe,_YHE_SMALL_,_YHE_BIG_);

  /** - tests in order to prevent divisions by zero */
  class_test(pth->YHe == 1.,
             pth->error_message,
             "stop to avoid division by zero");

  /** - test initial condition for recombination */
  class_test(ppr->thermo_z_initial < ppr->recfast_z_He_3,
             pth->error_message,
             "increase z_initial in order to start before HeliumIII recombination");

  return _SUCCESS_;
}

/**
 * Initialize the thermodynamics workspace.
 *
 * The workspace contains the arrays used for solving differential
 * equations (dubbed thermo_diffeq_workspace), and storing all
 * approximations, reionization parameters, heating parameters.
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param ptw        Input/Output: pointer to thermodynamics workspace
 * @return the error status
 */

int thermodynamics_workspace_init(
                                  struct precision * ppr,
                                  struct background * pba,
                                  struct thermodynamics * pth,
                                  struct thermo_workspace * ptw
                                  ) {

  /** Summary: */

  /** Define local variables */
  int index_ap;

  /** - number of z values */
  ptw->Nz_reco_lin = ppr->thermo_Nz_lin;
  ptw->Nz_reco_log = ppr->thermo_Nz_log;
  ptw->Nz_reco = ptw->Nz_reco_lin + ptw->Nz_reco_log;
  ptw->Nz_reio = ppr->reionization_z_start_max / ppr->reionization_sampling;
  ptw->Nz_tot = ptw->Nz_reio + ptw->Nz_reco;

  /** - relevant cosmological parameters */

  /* primordial helium mass fraction */
  ptw->YHe = pth->YHe;
  /* primordial helium-to-hydrogen nucleon ratio */
  ptw->fHe = pth->fHe;
  /* Hubble parameter today in SI units */
  ptw->SIunit_H0 = pba->H0 * _c_ / _Mpc_over_m_;
  /* H number density today in SI units*/
  ptw->SIunit_nH0 = 3.*ptw->SIunit_H0*ptw->SIunit_H0*pba->Omega0_b/(8.*_PI_*_G_*_m_H_)*(1.-ptw->YHe);
  /* CMB temperature today in Kelvin */
  ptw->Tcmb = pba->T_cmb;

  /** - relevant constants */

  /* Prefactor in non-relativistic number density for temperature -- (2*pi*m_e) and unit conversion */
  ptw->const_NR_numberdens = 2.*_PI_*(_m_e_/_h_P_)*(_k_B_/_h_P_);
  /* Ionization energy for HI -- temperature equivalent in Kelvin */
  ptw->const_Tion_H = _h_P_*_c_*_L_H_ion_/_k_B_;
  /* Ionization energy for HeI -- temperature equivalent in Kelvin */
  ptw->const_Tion_HeI = _h_P_*_c_*_L_He1_ion_/_k_B_;
  /* Ionization energy for HeII -- temperature equivalent in Kelvin */
  ptw->const_Tion_HeII = _h_P_*_c_*_L_He2_ion_/_k_B_;

  /* the field reionization_optical_depth is computed and filled later */

  /** - Allocate and initialize differential equation workspace */
  class_alloc(ptw->ptdw,
              sizeof(struct thermo_diffeq_workspace),
              pth->error_message);

  // Initialize ionisation fraction.
  ptw->ptdw->x_reio = 1.+2.*ptw->fHe;
  ptw->ptdw->x_noreio = 1.+2.*ptw->fHe;

  /** - define approximations */
  index_ap=0;
  /* Approximations have to appear in chronological order here! */
  class_define_index(ptw->ptdw->index_ap_brec,_TRUE_,index_ap,1); // before H- and He-recombination
  class_define_index(ptw->ptdw->index_ap_He1,_TRUE_,index_ap,1);  // during 1st He-recombination (HeIII)
  class_define_index(ptw->ptdw->index_ap_He1f,_TRUE_,index_ap,1); // in between 1st and 2nd He recombination
  class_define_index(ptw->ptdw->index_ap_He2,_TRUE_,index_ap,1);  // beginning of 2nd He-recombination (HeII)
  class_define_index(ptw->ptdw->index_ap_H,_TRUE_,index_ap,1);    // beginning of H-recombination (HI)
  class_define_index(ptw->ptdw->index_ap_frec,_TRUE_,index_ap,1); // during and after full H- and HeII-recombination
  class_define_index(ptw->ptdw->index_ap_reio,_TRUE_,index_ap,1); // during reionization
  ptw->ptdw->ap_size=index_ap;

  /*fix current approximation scheme to before recombination */
  ptw->ptdw->ap_current = ptw->ptdw->index_ap_brec;

  /** - store all ending redshifts for each approximation */
  class_alloc(ptw->ptdw->ap_z_limits,ptw->ptdw->ap_size*sizeof(double),pth->error_message);

  ptw->ptdw->ap_z_limits[ptw->ptdw->index_ap_brec] =
    ppr->recfast_z_He_1+ppr->recfast_delta_z_He_1; // beginning 1st He-recombination
  ptw->ptdw->ap_z_limits[ptw->ptdw->index_ap_He1] =
    ppr->recfast_z_He_2+ppr->recfast_delta_z_He_2; // end 1st He-recombination
  ptw->ptdw->ap_z_limits[ptw->ptdw->index_ap_He1f] =
    ppr->recfast_z_He_3+ppr->recfast_delta_z_He_3; // beginning 2nd He-recombination
  ptw->ptdw->ap_z_limits[ptw->ptdw->index_ap_He2] =
    ppr->recfast_z_early_H_recombination;          // beginning early H-recombination
  ptw->ptdw->ap_z_limits[ptw->ptdw->index_ap_H] =
    ppr->recfast_z_full_H_recombination;           // beginning full recombination equations
  ptw->ptdw->ap_z_limits[ptw->ptdw->index_ap_frec] =
    ppr->reionization_z_start_max;                 // beginning reionization
  ptw->ptdw->ap_z_limits[ptw->ptdw->index_ap_reio] = 0.0; // today

  /** - store smoothing deltas for transitions at the beginning of each aproximation */
  class_alloc(ptw->ptdw->ap_z_limits_delta,ptw->ptdw->ap_size*sizeof(double),pth->error_message);

  ptw->ptdw->ap_z_limits_delta[ptw->ptdw->index_ap_brec] = 0.;
  ptw->ptdw->ap_z_limits_delta[ptw->ptdw->index_ap_He1] = ppr->recfast_delta_z_He_1;
  ptw->ptdw->ap_z_limits_delta[ptw->ptdw->index_ap_He1f] = ppr->recfast_delta_z_He_2;
  ptw->ptdw->ap_z_limits_delta[ptw->ptdw->index_ap_He2] = ppr->recfast_delta_z_He_3;
  ptw->ptdw->ap_z_limits_delta[ptw->ptdw->index_ap_H] = ppr->recfast_delta_z_early_H_recombination;
  ptw->ptdw->ap_z_limits_delta[ptw->ptdw->index_ap_frec] = ppr->recfast_delta_z_full_H_recombination;
  ptw->ptdw->ap_z_limits_delta[ptw->ptdw->index_ap_reio] = ppr->recfast_delta_z_reio;

  /* With recombination computed by HyRec or Recfast, we need to allocate and initialize the wrappers */

  switch (pth->recombination) {

  case hyrec:
    class_alloc(ptw->ptdw->phyrec,
                sizeof(struct thermohyrec),
                pth->error_message);

    ptw->ptdw->phyrec->thermohyrec_verbose = pth->hyrec_verbose;
    class_call(thermodynamics_hyrec_init(ppr,pba,pth,ptw->SIunit_nH0,pba->T_cmb,ptw->fHe, ptw->ptdw->ap_z_limits[ptw->ptdw->index_ap_brec],ptw->ptdw->phyrec),
               ptw->ptdw->phyrec->error_message,
               pth->error_message);
    break;

  case recfast:
    class_alloc(ptw->ptdw->precfast,
                sizeof(struct thermorecfast),
                pth->error_message);

    class_call(recfast_init(ppr,pba,pth,ptw->ptdw->precfast,pth->recfast_photoion_mode,ptw->fHe),
               ptw->ptdw->precfast->error_message,
               pth->error_message);

    break;
  }

  /** - Allocate reionisation parameter workspace */
  class_alloc(ptw->ptrp,
              sizeof(struct thermo_reionization_parameters),
              pth->error_message);

  return _SUCCESS_;

}

/**
 * Assign value to each relevant index in vectors of thermodynamical
 * quantities, and the reionization parameters
 *
 * @param pba   Input: pointer to background structure
 * @param pth   Input/Output: pointer to thermodynamics structure
 * @param ptw   Input/Output: pointer to thermo workspace
 * @return the error status
 */

int thermodynamics_indices(
                           struct background * pba,
                           struct thermodynamics * pth,
                           struct thermo_workspace * ptw
                           ) {

  /** Summary: */

  /** - define local variables */
  struct thermo_reionization_parameters* ptrp = ptw->ptrp;
  /* a running index for the vector of thermodynamics quantities */
  int index_th;
  /* a running index for the vector of reionization parameters */
  int index_re;

  /** - initialization of all indices and flags in thermodynamics structure */
  index_th = 0;

  /* Free electron fraction */
  class_define_index(pth->index_th_xe,_TRUE_,index_th,1);
  /* Optical depth and related quantities */
  class_define_index(pth->index_th_dkappa,_TRUE_,index_th,1);
  class_define_index(pth->index_th_ddkappa,_TRUE_,index_th,1);
  class_define_index(pth->index_th_dddkappa,_TRUE_,index_th,1);
  class_define_index(pth->index_th_exp_m_kappa,_TRUE_,index_th,1);
  /* Visibility function + derivatives */
  class_define_index(pth->index_th_g,_TRUE_,index_th,1);
  class_define_index(pth->index_th_dg,_TRUE_,index_th,1);
  class_define_index(pth->index_th_ddg,_TRUE_,index_th,1);
  /* Baryon quantities, Temperature, Sound Speed, Drag time end */
  class_define_index(pth->index_th_Tb,_TRUE_,index_th,1);
  class_define_index(pth->index_th_dTb,_TRUE_,index_th,1);
  class_define_index(pth->index_th_wb,_TRUE_,index_th,1);
  class_define_index(pth->index_th_cb2,_TRUE_,index_th,1);
  class_define_index(pth->index_th_tau_d,_TRUE_,index_th,1);
  /* Derivatives of baryon sound speed (only computed if some non-minimal tight-coupling schemes is requested) */
  class_define_index(pth->index_th_dcb2,pth->compute_cb2_derivatives,index_th,1);
  class_define_index(pth->index_th_ddcb2,pth->compute_cb2_derivatives,index_th,1);
  /* IDM-DR quantities */
  class_define_index(pth->index_th_dmu_idm_dr,pba->has_idm_dr,index_th,1);
  class_define_index(pth->index_th_ddmu_idm_dr,pba->has_idm_dr,index_th,1);
  class_define_index(pth->index_th_dddmu_idm_dr,pba->has_idm_dr,index_th,1);
  class_define_index(pth->index_th_tau_idm_dr,pba->has_idm_dr,index_th,1);
  class_define_index(pth->index_th_tau_idr,pba->has_idm_dr,index_th,1);
  class_define_index(pth->index_th_g_idm_dr,pba->has_idm_dr,index_th,1);
  class_define_index(pth->index_th_cidm_dr2,pba->has_idm_dr,index_th,1);
  class_define_index(pth->index_th_Tidm_dr,pba->has_idm_dr,index_th,1);
  class_define_index(pth->index_th_dmu_idr,pba->has_idm_dr,index_th,1);
  /* Quantity defining the stepsize in perturbations.c */
  class_define_index(pth->index_th_rate,_TRUE_,index_th,1);
  /* Damping scale */
  class_define_index(pth->index_th_r_d,pth->compute_damping_scale,index_th,1);

  /* end of thermodynamics indices */

  pth->th_size = index_th;

  /** - initialization of all indices of parameters of reionization function */

  index_re=0;

  class_define_index(ptrp->index_re_reio_start,_TRUE_,index_re,1);

  switch (pth->reio_parametrization) {

    /* case with no reionization requested */
  case reio_none:
    class_define_index(ptrp->index_re_xe_before,_TRUE_,index_re,1);
    break;

    /* case where x_e(z) taken like in CAMB (other cases can be added) */
  case reio_camb:
  case reio_half_tanh:
    class_define_index(ptrp->index_re_reio_redshift,_TRUE_,index_re,1);
    class_define_index(ptrp->index_re_reio_exponent,_TRUE_,index_re,1);
    class_define_index(ptrp->index_re_reio_width,_TRUE_,index_re,1);
    class_define_index(ptrp->index_re_xe_before,_TRUE_,index_re,1);
    class_define_index(ptrp->index_re_xe_after,_TRUE_,index_re,1);
    class_define_index(ptrp->index_re_helium_fullreio_fraction,_TRUE_,index_re,1);
    class_define_index(ptrp->index_re_helium_fullreio_redshift,_TRUE_,index_re,1);
    class_define_index(ptrp->index_re_helium_fullreio_width,_TRUE_,index_re,1);
    break;

    /* case where x_e(z) is interpolated by tanh between plateaux */
  case reio_bins_tanh:

    /* the code will not only copy here the "bin centers" passed in input. It will add an initial and final value for (z,xe). So
       this array has a dimension bigger than the bin center array */

    ptrp->re_z_size=pth->binned_reio_num+2; /* add two values: beginning and end of reio */

    class_define_index(ptrp->index_re_first_z,_TRUE_,index_re,ptrp->re_z_size);
    class_define_index(ptrp->index_re_first_xe,_TRUE_,index_re,ptrp->re_z_size);
    class_define_index(ptrp->index_re_step_sharpness,_TRUE_,index_re,1);
    class_define_index(ptrp->index_re_xe_before,_TRUE_,index_re,1);
    break;

    /* case where x_e(z) is interpolated by tanh between knots */
  case reio_many_tanh:

    /* the code will not only copy here the "jump centers" passed in input. It will add an initial and final value for (z,xe). So
       this array has a dimension bigger than the jump center array */

    ptrp->re_z_size=pth->many_tanh_num+2; /* add two values: beginning and end of reio */

    class_define_index(ptrp->index_re_first_z,_TRUE_,index_re,ptrp->re_z_size);
    class_define_index(ptrp->index_re_first_xe,_TRUE_,index_re,ptrp->re_z_size);
    class_define_index(ptrp->index_re_step_sharpness,_TRUE_,index_re,1);
    class_define_index(ptrp->index_re_xe_before,_TRUE_,index_re,1);
    break;

    /* case where x_e(z) is linearly interpolated between knots */
  case reio_inter:

    ptrp->re_z_size=pth->reio_inter_num;

    class_define_index(ptrp->index_re_first_z,_TRUE_,index_re,ptrp->re_z_size);
    class_define_index(ptrp->index_re_first_xe,_TRUE_,index_re,ptrp->re_z_size);
    class_define_index(ptrp->index_re_xe_before,_TRUE_,index_re,1);
    break;

  default:
    class_stop(pth->error_message,
               "value of reio_parametrization=%d unclear",pth->reio_parametrization);
    break;
  }
  ptrp->re_size = index_re;

  return _SUCCESS_;

}

/**
 * Initialize the lists (of redshift, tau, etc.) of the thermodynamics struct
 *
 * @param ppr   Input: pointer to precision structure
 * @param pba   Input: pointer to background structure
 * @param pth   Input/Output: pointer to thermodynamics structure
 * @param ptw   Input: pointer to thermo workspace
 * @return the error status
 */

int thermodynamics_lists(
                         struct precision * ppr,
                         struct background* pba,
                         struct thermodynamics* pth,
                         struct thermo_workspace* ptw
                         ) {

  /** Summary: */

  /** Define local variables */
  int index_tau, index_z;
  double zinitial,zlinear;

  pth->tt_size = ptw->Nz_tot;

  /** - allocate tables*/
  class_alloc(pth->tau_table,pth->tt_size*sizeof(double),pth->error_message);
  class_alloc(pth->z_table,pth->tt_size*sizeof(double),pth->error_message);
  class_alloc(pth->thermodynamics_table,pth->th_size*pth->tt_size*sizeof(double),pth->error_message);
  class_alloc(pth->d2thermodynamics_dz2_table,pth->th_size*pth->tt_size*sizeof(double),pth->error_message);

  /** - define time sampling */

  /* Initial z, and the z at which we switch to linear sampling */
  zinitial = ppr->thermo_z_initial;
  zlinear  = ppr->thermo_z_linear;

  /* -> Between z_initial and z_linear, we use the spacing of recombination sampling */
  for(index_z=0; index_z <ptw->Nz_reco_log; index_z++) {
    pth->z_table[(pth->tt_size-1) - index_z] = -(-exp((log(zinitial)-log(zlinear))*(double)(ptw->Nz_reco_log-1-index_z) / (double)(ptw->Nz_reco_log-1)+log(zlinear)));
  }

  /* -> Between z_linear and reionization_z_start_max, we use the spacing of recombination sampling */
  for(index_z=0; index_z <ptw->Nz_reco_lin; index_z++) {
    pth->z_table[(pth->tt_size-1)-(index_z+ptw->Nz_reco_log)] = -(-(zlinear-ppr->reionization_z_start_max) * (double)(ptw->Nz_reco_lin-1-index_z) / (double)(ptw->Nz_reco_lin) - ppr->reionization_z_start_max);
  }

  /* -> Between reionization_z_start_max and 0, we use the spacing of reionization sampling, leaving out the first point to not double-count it */
  for(index_z=0; index_z <ptw->Nz_reio; index_z++) {
    pth->z_table[(pth->tt_size-1)-(index_z+ptw->Nz_reco)] = -(-ppr->reionization_z_start_max * (double)(ptw->Nz_reio-1-index_z) / (double)(ptw->Nz_reio));
  }

  for (index_tau=0; index_tau < pth->tt_size; index_tau++) {
    class_call(background_tau_of_z(pba,
                                   pth->z_table[index_tau],
                                   pth->tau_table+index_tau),
               pba->error_message,
               pth->error_message);
  }

  /** - store initial value of conformal time in the structure */
  pth->tau_ini = pth->tau_table[pth->tt_size-1];

  return _SUCCESS_;
}

/**
 * This routine initializes reionization_parameters for the chosen scheme of reionization function.
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param preio      Input/Output: pointer to the reionization parameters structure
 * @return the error status
 */

int thermodynamics_set_parameters_reionization(
                                               struct precision * ppr,
                                               struct background * pba,
                                               struct thermodynamics * pth,
                                               struct thermo_reionization_parameters * preio
                                               ) {

  /** Summary: */

  /** Define local variables */
  int bin;
  int point;
  double xe_input,xe_actual,z_sup;

  /** - allocate the vector of parameters defining the function \f$ X_e(z) \f$ */
  class_alloc(preio->reionization_parameters,preio->re_size*sizeof(double),pth->error_message);

  class_test(ppr->reionization_sampling <= 0.0,
             pth->error_message,
             "stop to avoid division by zero. Reionization stepsize has to be larger than zero");

  switch (pth->reio_parametrization) {

    /** - (a) no reionization */
  case reio_none:
    preio->reionization_parameters[preio->index_re_reio_start] = 0.;
    break;

    /** - (b) if reionization implemented like in CAMB, or half tanh like in  1209.0247 */
  case reio_camb:
  case reio_half_tanh:

    /** - --> set values of these parameters, excepted those depending on the reionization redshift */

    if (pth->reio_parametrization == reio_camb) {
      /* xe_after_reio: H + singly ionized He (checked before that denominator is non-zero) */
      preio->reionization_parameters[preio->index_re_xe_after] = 1. + pth->YHe/(_not4_*(1.-pth->YHe));
    }
    if (pth->reio_parametrization == reio_half_tanh) {
      /* xe_after_reio: neglect He ionization */
      preio->reionization_parameters[preio->index_re_xe_after] = 1.;
      //+ 2*pth->YHe/(_not4_*(1.-pth->YHe));    /* xe_after_reio: H + fully ionized He */
    }

    preio->reionization_parameters[preio->index_re_reio_exponent] = pth->reionization_exponent; /* reio_exponent */
    preio->reionization_parameters[preio->index_re_reio_width] = pth->reionization_width;    /* reio_width */
    preio->reionization_parameters[preio->index_re_helium_fullreio_fraction] = pth->YHe/(_not4_*(1.-pth->YHe)); /* helium_fullreio_fraction (checked before that denominator is non-zero) */
    preio->reionization_parameters[preio->index_re_helium_fullreio_redshift] = pth->helium_fullreio_redshift; /* helium_fullreio_redshift */
    preio->reionization_parameters[preio->index_re_helium_fullreio_width] = pth->helium_fullreio_width;    /* helium_fullreio_width */

    class_test(preio->reionization_parameters[preio->index_re_reio_exponent]==0,
               pth->error_message,
               "stop to avoid division by zero");

    class_test(preio->reionization_parameters[preio->index_re_reio_width]==0,
               pth->error_message,
               "stop to avoid division by zero");

    class_test(preio->reionization_parameters[preio->index_re_helium_fullreio_width]==0,
               pth->error_message,
               "stop to avoid division by zero");

    /** - --> if reionization redshift given as an input, initialize the remaining values*/

    if (pth->reio_z_or_tau == reio_z) {

      /* reionization redshift */
      preio->reionization_parameters[preio->index_re_reio_redshift] = pth->z_reio;

      /* infer starting redshift for hydrogen */

      if (pth->reio_parametrization == reio_camb) {

        preio->reionization_parameters[preio->index_re_reio_start] = preio->reionization_parameters[preio->index_re_reio_redshift]+
                                                                  ppr->reionization_start_factor*pth->reionization_width;

        /* if starting redshift for helium is larger, take that one (does not happen in realistic models) */
        if (preio->reionization_parameters[preio->index_re_reio_start] <
            pth->helium_fullreio_redshift+ppr->reionization_start_factor*pth->helium_fullreio_width)

          preio->reionization_parameters[preio->index_re_reio_start] =
            pth->helium_fullreio_redshift+ppr->reionization_start_factor*pth->helium_fullreio_width;

      }
      else {

        preio->reionization_parameters[preio->index_re_reio_start] = pth->z_reio;
      }

      class_test(preio->reionization_parameters[preio->index_re_reio_start] > ppr->reionization_z_start_max,
                 pth->error_message,
                 "starting redshift for reionization > reionization_z_start_max = %e\n",ppr->reionization_z_start_max);

    }

    /** - --> if reionization optical depth given as an input, find reionization redshift by bisection and initialize the remaining values */
    if (pth->reio_z_or_tau == reio_tau) {
      z_sup = ppr->reionization_z_start_max-ppr->reionization_start_factor*pth->reionization_width;
      class_test(z_sup < 0.,
                 pth->error_message,
                 "parameters are such that reionization cannot take place before today while starting after z_start_max; need to increase z_start_max");

      /* maximum possible reionization redshift */
      preio->reionization_parameters[preio->index_re_reio_redshift] = z_sup;
      /* maximum possible starting redshift */
      preio->reionization_parameters[preio->index_re_reio_start] = ppr->reionization_z_start_max;
    }
    break;

    /** - (c) if reionization implemented with reio_bins_tanh scheme */
  case reio_bins_tanh:

    /* this algorithm requires at least two bin centers (i.e. at least 4 values in the (z,xe) array, counting the edges). */
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

    /* the code will not only copy here the "bin centers" passed in input. It will add an initial and final value for (z,xe).
       First, fill all entries except the first and the last */
    for (bin=1; bin<preio->re_z_size-1; bin++) {
      preio->reionization_parameters[preio->index_re_first_z+bin] = pth->binned_reio_z[bin-1];
      preio->reionization_parameters[preio->index_re_first_xe+bin] = pth->binned_reio_xe[bin-1];
    }

    /* find largest value of z in the array. We choose to define it as z_(i_max) + 2*(the distance between z_(i_max) and z_(i_max-1)). E.g. if
       the bins are in 10,12,14, the largest z will be 18. */
    preio->reionization_parameters[preio->index_re_first_z+preio->re_z_size-1] =
      preio->reionization_parameters[preio->index_re_first_z+preio->re_z_size-2]
      +2.*(preio->reionization_parameters[preio->index_re_first_z+preio->re_z_size-2]
        -preio->reionization_parameters[preio->index_re_first_z+preio->re_z_size-3]);

    /* copy this value in reio_start */
    preio->reionization_parameters[preio->index_re_reio_start] = preio->reionization_parameters[preio->index_re_first_z+preio->re_z_size-1];

    /* check it's not too big */
    class_test(preio->reionization_parameters[preio->index_re_reio_start] > ppr->reionization_z_start_max,
               pth->error_message,
               "starting redshift for reionization = %e, reionization_z_start_max = %e, you must change the binning or increase reionization_z_start_max",
               preio->reionization_parameters[preio->index_re_reio_start],
               ppr->reionization_z_start_max);

    /* find smallest value of z in the array. We choose to define it as z_0 - (the distance between z_1 and z_0). E.g. if
       the bins are in 10,12,14, the stop redshift will be 8. */
    preio->reionization_parameters[preio->index_re_first_z] =
      2.*preio->reionization_parameters[preio->index_re_first_z+1]
      -preio->reionization_parameters[preio->index_re_first_z+2];

    /* check it's not too small */
    /* 6.06.2015: changed this test to simply imposing that the first z is at least zero */
    /*
    class_test(preio->reionization_parameters[preio->index_re_first_z] < 0,
               pth->error_message,
               "final redshift for reionization = %e, you must change the binning or redefine the way in which the code extrapolates below the first value of z_i",preio->reionization_parameters[preio->index_re_first_z]);
    */
    if (preio->reionization_parameters[preio->index_re_first_z] < 0) {
      preio->reionization_parameters[preio->index_re_first_z] = 0.;
    }

    /* infer xe after reio */
    preio->reionization_parameters[preio->index_re_first_xe] = 1. + pth->YHe/(_not4_*(1.-pth->YHe)); /* xe_after_reio: H + singly ionized He (note: segmentation fault impossible,
                                                                                                          checked before that denominator is non-zero) */

    /* pass step sharpness parameter */
    preio->reionization_parameters[preio->index_re_step_sharpness] = pth->binned_reio_step_sharpness;
    break;

    /** - (d) if reionization implemented with reio_many_tanh scheme */
  case reio_many_tanh:

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

    /* the code will not only copy here the "jump centers" passed in input. It will add an initial and final value for (z,xe).
       First, fill all entries except the first and the last */
    for (bin=1; bin<preio->re_z_size-1; bin++) {

      preio->reionization_parameters[preio->index_re_first_z+bin] = pth->many_tanh_z[bin-1];

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

      preio->reionization_parameters[preio->index_re_first_xe+bin] = xe_actual;
    }

    /* find largest value of z in the array. We choose to define it as z_(i_max) + ppr->reionization_start_factor*step_sharpness. */
    preio->reionization_parameters[preio->index_re_first_z+preio->re_z_size-1] =
      preio->reionization_parameters[preio->index_re_first_z+preio->re_z_size-2]
      +ppr->reionization_start_factor*pth->many_tanh_width;

    /* copy this value in reio_start */
    preio->reionization_parameters[preio->index_re_reio_start] = preio->reionization_parameters[preio->index_re_first_z+preio->re_z_size-1];

    /* check it's not too big */
    class_test(preio->reionization_parameters[preio->index_re_reio_start] > ppr->reionization_z_start_max,
               pth->error_message,
               "starting redshift for reionization = %e, reionization_z_start_max = %e, you must change the binning or increase reionization_z_start_max",
               preio->reionization_parameters[preio->index_re_reio_start],
               ppr->reionization_z_start_max);

    /* find smallest value of z in the array. We choose to define it as z_0 - ppr->reionization_start_factor*step_sharpness, but at least zero. */
    preio->reionization_parameters[preio->index_re_first_z] =
      preio->reionization_parameters[preio->index_re_first_z+1]
      -ppr->reionization_start_factor*pth->many_tanh_width;

    if (preio->reionization_parameters[preio->index_re_first_z] < 0) {
      preio->reionization_parameters[preio->index_re_first_z] = 0.;
    }

    /* infer xe after reio */
    preio->reionization_parameters[preio->index_re_first_xe] = preio->reionization_parameters[preio->index_re_first_xe+1];

    /* if we want to model only hydrogen reionization and neglect both helium reionization */
    //preio->reionization_parameters[preio->index_re_first_xe] = 1.;

    /* if we want to model only hydrogen + first helium reionization and neglect second helium reionization */
    //preio->reionization_parameters[preio->index_re_first_xe] = 1. + pth->YHe/(_not4_*(1.-pth->YHe));

    /* if we want to model hydrogen + two helium reionization */
    //preio->reionization_parameters[preio->index_re_first_xe] = 1. + 2.*pth->YHe/(_not4_*(1.-pth->YHe));

    /* pass step sharpness parameter */
    class_test(pth->many_tanh_width<=0,
               pth->error_message,
               "many_tanh_width must be strictly positive, you passed %e",
               pth->many_tanh_width);

    preio->reionization_parameters[preio->index_re_step_sharpness] = pth->many_tanh_width;
    break;

    /** - (e) if reionization implemented with reio_inter scheme */
  case reio_inter:

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

    /* this parametrization requires that the last x_i value is zero (the code will substitute it with the value that one would get in
       absence of reionization, as compute by the recombination code) */
    class_test(pth->reio_inter_xe[pth->reio_inter_num-1] != 0.,
               pth->error_message,
               "For reio_inter scheme, the last value of reio_inter_xe[...]  should always be zero, you passed %e",
               pth->reio_inter_xe[pth->reio_inter_num-1]);

    /* copy here the (z,xe) values passed in input. */
    for (point=0; point<preio->re_z_size; point++) {

      preio->reionization_parameters[preio->index_re_first_z+point] = pth->reio_inter_z[point];

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

      preio->reionization_parameters[preio->index_re_first_xe+point] = xe_actual;
    }

    /* copy highest redshift in reio_start */
    preio->reionization_parameters[preio->index_re_reio_start] = preio->reionization_parameters[preio->index_re_first_z+preio->re_z_size-1];

    /* check it's not too big */
    class_test(preio->reionization_parameters[preio->index_re_reio_start] > ppr->reionization_z_start_max,
               pth->error_message,
               "starting redshift for reionization = %e, reionization_z_start_max = %e, you must change the binning or increase reionization_z_start_max",
               preio->reionization_parameters[preio->index_re_reio_start],
               ppr->reionization_z_start_max);
    break;

  default:
    class_stop(pth->error_message,
               "value of reio_parametrization=%d unclear",pth->reio_parametrization);
    break;
  }

  return _SUCCESS_;

}

/**
 * Integrate thermodynamics with your favorite recombination code. The default options are HyRec and RecFastCLASS.
 *
 * Integrate thermodynamics with HyRec or Recfast, allocate and fill part of the thermodynamics interpolation table (the rest is filled in
 * thermodynamics_calculate_remaining_quantitie).
 *
 * Version modified by Daniel Meinert and Nils Schoeneberg to use the ndf15 evolver or any other evolver inherent to CLASS,
 * modified again by Nils Schoeneberg to use wrappers.
 *
 * @param ppr      Input: pointer to precision structure
 * @param pba      Input: pointer to background structure
 * @param pth      Input/Output: pointer to thermodynamics structure where results are stored
 * @param ptw      Input: pointer to thermo_workspace structure used to communicate with generic evolver
 * @param pvecback Input: pointer to an allocated (but empty) vector of background variables
 * @return the error status
 *
 * Integrate thermodynamics with your favorite recombination code. The default options are HyRec and Recfast.
 */

int thermodynamics_solve(
                         struct precision * ppr,
                         struct background * pba,
                         struct thermodynamics * pth,
                         struct thermo_workspace * ptw,
                         double * pvecback
                         ) {
  /** Summary: */

  /** - define local variables */
  /* Index of current approximation scheme */
  int index_ap;
  /* number of time intervals of one approximation scheme */
  int interval_number;
  /* index running over such time intervals */
  int index_interval;
  /* edge of intervals where approximation scheme is uniform: z_ini, z_switch_1, ..., z_end */
  double * interval_limit;
  /* other z sampling variables */
  int i;
  double * mz_output;

  /* contains all fixed parameters which should be passed to thermodynamics_derivs */
  struct thermodynamics_parameters_and_workspace tpaw;

  /* function pointer to ODE evolver and names of possible evolvers. */
  extern int evolver_rk();
  extern int evolver_ndf15();
  int (*generic_evolver)() = evolver_ndf15;

  /** - choose evolver */
  switch (ppr->thermo_evolver) {
  case rk:
    generic_evolver = evolver_rk;
    break;
  case ndf15:
    generic_evolver = evolver_ndf15;
    break;
  }

  /** - define the fields of the 'thermodynamics parameter and workspace' structure */
  tpaw.pba = pba;
  tpaw.ppr = ppr;
  tpaw.pth = pth;
  tpaw.pvecback = pvecback;
  tpaw.ptw = ptw;

  /** - define time sampling: create a local array of minus z values
        called mz (from mz=-zinitial growing towards mz=0) */

  class_alloc(mz_output,pth->tt_size*sizeof(double), pth->error_message);
  for(i=0; i < pth->tt_size; ++i) {
    mz_output[i] = -pth->z_table[pth->tt_size-1-i];
  }

  /** - define intervals for each approximation scheme */

  /* create the array of interval limits */
  class_alloc(interval_limit,(ptw->ptdw->ap_size+1)*sizeof(double),pth->error_message);
  /* fix interval number to number of approximations */
  interval_number = ptw->ptdw->ap_size;
  /* integration starts at z_ini and ends at z_end */
  interval_limit[0]= mz_output[0];
  interval_limit[ptw->ptdw->ap_size] = mz_output[ptw->Nz_tot-1];
  /* each interval ends with the proper ending redshift of its approximation */
  for(index_ap=0; index_ap < ptw->ptdw->ap_size-1; index_ap++) {
    interval_limit[index_ap+1] = -ptw->ptdw->ap_z_limits[index_ap];
  }

  /** - loop over intervals over which approximation scheme is
        uniform. For each interval: */

  for (index_interval=0; index_interval<interval_number; index_interval++) {

    /** - --> (a) fix current approximation scheme. */

    ptw->ptdw->ap_current = index_interval;

    /** - --> (b) define the vector of quantities to be integrated
                  over. If the current interval starts from the
                  initial time zinitial, fill the vector with initial
                  conditions. If it starts from an approximation
                  switching point, redistribute correctly the values
                  from the previous to the new vector. For both
                  RECFAST and HYREC, the vector consists of Tmat, x_H,
                  x_He, + others for exotic models */

    class_call(thermodynamics_vector_init(ppr,
                                          pba,
                                          pth,
                                          interval_limit[index_interval],
                                          ptw),
               pth->error_message,
               pth->error_message);

    /* find the value of last_index_back at z = - interval_limit[index_interval], in order to speed up
       subsequent interpolations in thermodynamics_derivs */
    class_call(background_at_z(pba,
                               -interval_limit[index_interval],
                               normal_info,
                               inter_normal,
                               &(ptw->last_index_back),
                               pvecback),
               pba->error_message,
               pth->error_message);

    /** - --> (c1) If we have the optical depth tau_reio as input the
       last evolver step (reionization approximation) is done
       separately in a loop, to find the approximate redshift of
       reionization given the input of tau_reio, using a bisection
       method. This is similar to the general CLASS shooting method,
       but doing this step here is more davantageous since we only
       need to do repeatedly the last approximation step of the
       integration, instead of the full background and thermodynamics
       module */
    if ((pth->reio_z_or_tau == reio_tau) && (index_interval == ptw->ptdw->index_ap_reio)) {

      class_call(thermodynamics_reionization_evolve_with_tau(&tpaw,
                                                             interval_limit[index_interval],
                                                             interval_limit[index_interval+1],
                                                             mz_output,
                                                             pth->tt_size),
                 pth->error_message,
                 pth->error_message);
    }

    /** --> (c2) otherwise, just integrate quantities over the current interval. */
    else{

      class_call(generic_evolver(thermodynamics_derivs,
                                 interval_limit[index_interval],
                                 interval_limit[index_interval+1],
                                 ptw->ptdw->ptv->y,
                                 ptw->ptdw->ptv->used_in_output,
                                 ptw->ptdw->ptv->ti_size,
                                 &tpaw,
                                 ppr->tol_thermo_integration,
                                 ppr->smallest_allowed_variation,
                                 thermodynamics_timescale,  // timescale
                                 ppr->thermo_integration_stepsize, // stepsize = this number * timescale
                                 mz_output, // values of z for output
                                 pth->tt_size, // size of previous array
                                 thermodynamics_sources, // function for output
                                 NULL, // print variables
                                 pth->error_message),
                   pth->error_message,
                   pth->error_message);
    }

  }

  /** - Compute reionization optical depth, if not supplied as input parameter */
  if (pth->reio_z_or_tau == reio_z) {

    class_call(thermodynamics_reionization_get_tau(ppr,
                                                   pba,
                                                   pth,
                                                   ptw),
               pth->error_message,
               pth->error_message);

    pth->tau_reio=ptw->reionization_optical_depth;

  }

  /** - free quantities allocated at the beginning of the routine */
  if (ptw->ptdw->ap_size != 0) {
    class_call(thermodynamics_vector_free(ptw->ptdw->ptv),
               pth->error_message,
               pth->error_message);
  }

  free(interval_limit);
  free(mz_output);

  return _SUCCESS_;

}

/**
 * Calculate those thermodynamics quantities which are not inside of
 * the thermodynamics table already.
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param pth        Input/Output: pointer to initialized thermodynamics structure
 * @param pvecback   Input: pointer to some allocated pvecback
 * @return the error status
 */

int thermodynamics_calculate_remaining_quantities(
                                                  struct precision * ppr,
                                                  struct background * pba,
                                                  struct thermodynamics* pth,
                                                  double* pvecback
                                                  ) {

  /** Summary: */

  /* The temporary quantities stored in columns ddkappa and dddkappa will not be used anymore, so they can and WILL be overwritten by other
     intermediate steps of other computations */

  class_call(thermodynamics_calculate_conformal_drag_time(pba,pth,pvecback),
             pth->error_message,
             pth->error_message);

  if (pth->compute_damping_scale == _TRUE_) {
    class_call(thermodynamics_calculate_damping_scale(pba,pth,pvecback),
               pth->error_message,
               pth->error_message);
  }

  class_call(thermodynamics_calculate_opticals(ppr,pth),
             pth->error_message,
             pth->error_message);

  if (pba->has_idr == _TRUE_) {
    class_call(thermodynamics_calculate_idm_dr_quantities(ppr,pba,pth,pvecback),
               pth->error_message,
               pth->error_message);
  }

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

  class_call(thermodynamics_calculate_recombination_quantities(ppr,pba,pth,pvecback),
             pth->error_message,
             pth->error_message);

  class_call(thermodynamics_calculate_drag_quantities(ppr,pba,pth,pvecback),
             pth->error_message,
             pth->error_message);

  return _SUCCESS_;
}

/**
 * In verbose mode, print basic information on the thermal history
 *
 * @param pba   Input: pointer to background structure
 * @param pth   Input/Output: pointer to initialized thermodynamics structure
 * @return the error status
 */

int thermodynamics_output_summary(
                                  struct background* pba,
                                  struct thermodynamics* pth
                                  ) {

  /** Summary: */

  /** Define local variables */
  double tau_reio;

  /** - print the main results */
  if (pba->has_idm_dr == _TRUE_) {
    if ((pth->thermodynamics_table[(pth->tt_size-1)*pth->th_size+pth->index_th_tau_idm_dr]>1.) && (pth->thermodynamics_table[(pth->tt_size-1)*pth->th_size+pth->index_th_tau_idr]>1.)) {
      printf(" -> idr decouples at tau_idr = %e Mpc\n",pth->tau_idr);
      printf(" -> idm_dr decouples at tau_idm_dr = %e Mpc\n",pth->tau_idm_dr);
    }
    else{
      printf(" -> computation of decoupling time of idm_dr and idr skipped, because z would not be in z_table\n");
    }
  }
  printf(" -> recombination (maximum of visibility function) at z = %f\n",pth->z_rec);
  printf("    corresponding to conformal time = %f Mpc\n",pth->tau_rec);
  printf("    with comoving sound horizon = %f Mpc\n",pth->rs_rec);
  printf("    angular diameter distance = %f Mpc\n",pth->da_rec);
  printf("    sound horizon angle 100*theta_s = %f\n",100.*pth->rs_rec/pth->ra_rec);
  if (pth->compute_damping_scale == _TRUE_) {
    printf("    comoving photon damping scale = %f Mpc\n",pth->rd_rec);
    printf("    comoving damping wavenumber k_d = %f 1/Mpc\n",2.*_PI_/pth->rd_rec);
  }
  printf("    Thomson optical depth crosses one at z_* = %f\n",pth->z_star);
  printf("    giving an angle 100*theta_* = %f\n",100.*pth->rs_star/pth->ra_star);
  printf(" -> baryon drag stops at z = %f\n",pth->z_d);
  printf("    corresponding to conformal time = %f Mpc\n",pth->tau_d);
  printf("    with comoving sound horizon rs = %f Mpc\n",pth->rs_d);

  switch (pth->reio_parametrization) {

  case reio_none:
    /* the information returned by this line could be interesting when
       using reio_none + exotic energy ionjection, since there could
       still be a global minimum of x_e(z) */
    printf(" -> no reionization requested, optical depth = %f\n",pth->tau_reio);
    break;

  case reio_camb:
  case reio_half_tanh:
    switch (pth->reio_z_or_tau) {
    case reio_tau:
      printf(" -> reionization at z = %f\n",pth->z_reio);
      break;
    case reio_z:
      printf(" -> reionization with optical depth = %f\n",pth->tau_reio);
      break;
    }
    class_call(background_tau_of_z(pba,pth->z_reio,&tau_reio),
               pba->error_message,
               pth->error_message);
    printf("    corresponding to conformal time = %f Mpc\n",tau_reio);
    break;

  case reio_bins_tanh:
    printf(" -> binned reionization gives optical depth = %f\n",pth->tau_reio);
    break;

  case reio_many_tanh:
    printf(" -> many-step reionization gives optical depth = %f\n",pth->tau_reio);
    break;

  case reio_inter:
    printf(" -> interpolated reionization history gives optical depth = %f\n",pth->tau_reio);
    break;

  default:
    class_stop(pth->error_message,
               "value of reio_parametrization=%d unclear",pth->reio_parametrization);
    break;
  }

  if (pth->thermodynamics_verbose > 1) {
    printf(" -> free-streaming approximation can be turned on as soon as tau=%g Mpc\n",pth->tau_free_streaming);
  }
  if ((pba->has_idr)&&(pth->thermodynamics_verbose > 1)) {
    printf(" -> dark free-streaming approximation can be turned on as soon as tau=%g Mpc\n",
           pth->tau_idr_free_streaming);
  }

  return _SUCCESS_;
}

/**
 * Free the thermo_workspace structure (with the exception of the thermo_vector '->ptv' field, which is freed separately in
 * thermo_vector_free).
 *
 * @param pth        Input: pointer to initialized thermodynamics structure
 * @param ptw        Input: pointer to perturbations_workspace structure to be freed
 * @return the error status
 */
int thermodynamics_workspace_free(
                                  struct thermodynamics* pth,
                                  struct thermo_workspace * ptw
                                  ) {

  free(ptw->ptdw->ap_z_limits);
  free(ptw->ptdw->ap_z_limits_delta);

  switch (pth->recombination) {

  case hyrec:
    class_call(thermodynamics_hyrec_free(ptw->ptdw->phyrec),
               ptw->ptdw->phyrec->error_message,
               pth->error_message);
    free(ptw->ptdw->phyrec);
    break;

  case recfast:
    free(ptw->ptdw->precfast);
    break;
  }

  free(ptw->ptrp->reionization_parameters);
  free(ptw->ptdw);
  free(ptw->ptrp);

  free(ptw);

  return _SUCCESS_;
}

/**
 * Initialize the field '->ptv' of a thermo_diffeq_workspace structure, which is a thermo_vector structure. This structure contains indices
 * and values of all quantities which need to be integrated with respect to time (and only them: quantities fixed analytically or obeying
 * constraint equations are NOT included in this vector).
 *
 * The routine sets and allocates the vector y, dy and used_in_output with the right size depending on the current approximation scheme
 * stored in the workspace. Moreover the initial conditions for each approximation scheme are calculated and set correctly.
 *
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param mz         Input: negative redshift
 * @param ptw        Input/Output: pointer to thermodynamics workspace
 *
 * @return the error status
 */

int thermodynamics_vector_init(
                               struct precision * ppr,
                               struct background * pba,
                               struct thermodynamics * pth,
                               double mz,
                               struct thermo_workspace * ptw
                               ) {

  /** Summary: */

  /** Define local variables */
  int index_ti;
  /* ptdw->ptv unallocated if ap_current == index_ap_brec, allocated and filled otherwise */
  struct thermo_vector * ptv;
  struct thermo_diffeq_workspace * ptdw = ptw->ptdw;

  double z;

  class_alloc(ptv,sizeof(struct thermo_vector),pth->error_message);

  /* mz = Minus z is inverted*/
  z = -mz;

  /* Start from no component */
  index_ti = 0;

  /* Add common indices (Have to be added before) */
  class_define_index(ptv->index_ti_D_Tmat,_TRUE_,index_ti,1);

  /* Add all components that should be evolved */
  if (ptdw->ap_current == ptdw->index_ap_brec) {
    /* Nothing else to add */
  }
  else if (ptdw->ap_current == ptdw->index_ap_He1) {
    /* Nothing else to add */
  }
  else if (ptdw->ap_current == ptdw->index_ap_He1f) {
    /* Nothing else to add */
  }
  else if (ptdw->ap_current == ptdw->index_ap_He2) {
    /* Nothing else to add */
  }
  else if (ptdw->ap_current == ptdw->index_ap_H) {
    class_define_index(ptv->index_ti_x_He,_TRUE_,index_ti,1);
  }
  else if (ptdw->ap_current == ptdw->index_ap_frec) {
    class_define_index(ptv->index_ti_x_He,_TRUE_,index_ti,1);
    class_define_index(ptv->index_ti_x_H,_TRUE_,index_ti,1);
  }
  else if (ptdw->ap_current == ptdw->index_ap_reio) {
    class_define_index(ptv->index_ti_x_He,_TRUE_,index_ti,1);
    class_define_index(ptv->index_ti_x_H,_TRUE_,index_ti,1);
  }

  /* We have now obtained the full size */
  ptv->ti_size = index_ti;

  /* Allocate all arrays used during the evolution */
  class_calloc(ptv->y,ptv->ti_size,sizeof(double),pth->error_message);
  class_alloc(ptv->dy,ptv->ti_size*sizeof(double),pth->error_message);
  class_alloc(ptv->used_in_output,ptv->ti_size*sizeof(int),pth->error_message);

  for (index_ti=0; index_ti<ptv->ti_size; index_ti++) {
    ptv->used_in_output[index_ti] = _TRUE_;
  }

  /* setting intial conditions for each approximation: */
  /* - in the first scheme (brec = before recombination), we need initial condition for the matter temperature given by the photon temperature */
  if (ptdw->ap_current == ptdw->index_ap_brec) {
    /* Store Tmat in workspace for later use */
    ptdw->Tmat = ptw->Tcmb*(1.+z);

    /* Set the new vector and its indices */
    ptdw->ptv = ptv;

    ptdw->ptv->y[ptdw->ptv->index_ti_D_Tmat] = 0.;

    ptdw->require_H = _FALSE_;
    ptdw->require_He = _FALSE_;
  }
  /* - in this scheme we start to evolve Helium and thus need to set its initial condition via the analytic function */
  else if (ptdw->ap_current == ptdw->index_ap_H) {
    /* Store Tmat in workspace for later use */
    ptdw->Tmat = ptdw->ptv->y[ptdw->ptv->index_ti_D_Tmat] + ptw->Tcmb*(1.+z);

    /* Obtain initial contents of new vector analytically, especially x_He */
    class_call(thermodynamics_ionization_fractions(z,ptdw->ptv->y,pth,ptw,ptdw->ap_current-1),
               pth->error_message,
               pth->error_message);

    /* Set the new vector and its indices */
    ptv->y[ptv->index_ti_D_Tmat] = ptdw->ptv->y[ptdw->ptv->index_ti_D_Tmat];
    ptv->y[ptv->index_ti_x_He] = ptdw->x_He;

    /* Free the old vector and its indices */
    class_call(thermodynamics_vector_free(ptdw->ptv),
               pth->error_message,
               pth->error_message);

    /* Copy the new vector into the position of the old one*/
    ptdw->ptv = ptv;

    ptdw->require_H = _FALSE_;
    ptdw->require_He = _TRUE_;
  }
  /* - in the scheme of full recombination (=frec) we evolve all quantities and thus need to set their initial conditions.
       Tmat and x_He are solely taken from the previous scheme, x_H is set via the analytic function */
  else if (ptdw->ap_current == ptdw->index_ap_frec) {
    /* Store Tmat in workspace for later use */
    ptdw->Tmat = ptdw->ptv->y[ptdw->ptv->index_ti_D_Tmat] + ptw->Tcmb*(1.+z);

    /* Obtain initial contents of new vector analytically, especially x_H */
    class_call(thermodynamics_ionization_fractions(z,ptdw->ptv->y,pth,ptw,ptdw->ap_current-1),
               pth->error_message,
               pth->error_message);

    /* Set the new vector and its indices */
    ptv->y[ptv->index_ti_D_Tmat] = ptdw->ptv->y[ptdw->ptv->index_ti_D_Tmat];
    ptv->y[ptv->index_ti_x_H] = ptdw->x_H;
    ptv->y[ptv->index_ti_x_He] = ptdw->ptv->y[ptdw->ptv->index_ti_x_He];

    /* Free the old vector and its indices */
    class_call(thermodynamics_vector_free(ptdw->ptv),
               pth->error_message,
               pth->error_message);

    /* Copy the new vector into the position of the old one*/

    ptdw->ptv = ptv;

    ptdw->require_H = _TRUE_;
    ptdw->require_He = _TRUE_;
  }
  /* - during reionization we continue to evolve all quantities. Now all three intial conditions are just taken from the previous scheme */
  else if (ptdw->ap_current == ptdw->index_ap_reio) {

    /* Set the new vector and its indices */
    ptv->y[ptv->index_ti_D_Tmat] = ptdw->ptv->y[ptdw->ptv->index_ti_D_Tmat];
    ptv->y[ptv->index_ti_x_H] = ptdw->ptv->y[ptdw->ptv->index_ti_x_H];
    ptv->y[ptv->index_ti_x_He] = ptdw->ptv->y[ptdw->ptv->index_ti_x_He];

    /* Free the old vector and its indices */
    class_call(thermodynamics_vector_free(ptdw->ptv),
               pth->error_message,
               pth->error_message);

    /* Copy the new vector into the position of the old one*/

    ptdw->ptv = ptv;

    ptdw->require_H = _TRUE_;
    ptdw->require_He = _TRUE_;
  }
  /* - in all other approximations we only evolve Tmat and set its initial conditions from the previous scheme */
  else{
    /* Store Tmat in workspace for later use */
    ptdw->Tmat = ptdw->ptv->y[ptdw->ptv->index_ti_D_Tmat] + ptw->Tcmb*(1.+z);

    /* Set the new vector and its indices */
    ptv->y[ptv->index_ti_D_Tmat] = ptdw->ptv->y[ptdw->ptv->index_ti_D_Tmat];

    /* Free the old vector and its indices */
    class_call(thermodynamics_vector_free(ptdw->ptv),
               pth->error_message,
               pth->error_message);

    /* Copy the new vector into the position of the old one*/
    ptdw->ptv = ptv;

    ptdw->require_H = _FALSE_;
    ptdw->require_He = _FALSE_;
  }

  return _SUCCESS_;
}

/**
 * If the input for reionization is tau_reio, thermodynamics_solve()
 * calls this function instead of the evolver for dealing with the
 * last era (the reionization era).
 *
 * Instead of computing the evolution of quantities during
 * reionization for a fixed z_reio, as the evolver would do, this
 * function finds z_reio by bisection. First we make an initial guess
 * for z_reio with reionization_z_start_max and then find a z_reio
 * which leads to the given tau_reio (in the range of tolerance
 * reionization_optical_depth_tol).
 *
 * @param ptpaw      Input: pointer to parameters and workspace
 * @param mz_ini     Input: initial redshift
 * @param mz_end     Input: ending redshift
 * @param mz_output  Input: pointer to redshift array at which output should be written
 * @param mz_size    Input: number of redshift values in this array
 * @return the error status
 */

int thermodynamics_reionization_evolve_with_tau(
                                                struct thermodynamics_parameters_and_workspace * ptpaw,
                                                double mz_ini,
                                                double mz_end,
                                                double * mz_output,
                                                int mz_size
                                                ) {

  /** Summary: */

  /** Define local variables */
  int counter;
  double z_sup,z_mid,z_inf;
  double tau_sup,tau_mid,tau_inf;

  int index_ti;
  int last_index_back_mz_ini;

  struct precision * ppr;
  struct background * pba;
  struct thermodynamics * pth;
  struct thermo_workspace * ptw;

  /* function pointer to ODE evolver and names of possible evolvers */
  extern int evolver_rk();
  extern int evolver_ndf15();
  int (*generic_evolver)() = evolver_ndf15;

  /* pointers towards two thermo vector stuctures (see below) */

  struct thermo_vector * ptvs; // Vector for storing the initial conditions
  struct thermo_vector * ptv; // Temporary vector as workspace

  /** - Remame fields to avoid heavy notations */

  ppr = ptpaw->ppr;
  pba = ptpaw->pba;
  pth = ptpaw->pth;
  ptw = ptpaw->ptw;

  /** - Choose evolver */

  switch (ppr->thermo_evolver) {
  case rk:
    generic_evolver = evolver_rk;
    break;
  case ndf15:
    generic_evolver = evolver_ndf15;
    break;
  }

  /** - ptvs will be a pointer towards the same thermo vector that was
        used in the previous approximation schemes; it contains values
        that will serve here to set initial conditions. */

  ptvs = ptw->ptdw->ptv;

  /** - ptv is a pointer towards a whole new thermo vector used for
        the calculations in the bisection, that we must allocate and
        initialize */

  class_alloc(ptv,sizeof(struct thermo_vector),pth->error_message);

  /* allocate vector indices dynamically */
  index_ti = 0;
  class_define_index(ptv->index_ti_D_Tmat,_TRUE_,index_ti,1);
  class_define_index(ptv->index_ti_x_He,_TRUE_,index_ti,1);
  class_define_index(ptv->index_ti_x_H,_TRUE_,index_ti,1);
  ptv->ti_size = index_ti;

  /* Allocate all arrays used during the evolution */
  class_calloc(ptv->y,ptv->ti_size,sizeof(double),pth->error_message);
  class_alloc(ptv->dy,ptv->ti_size*sizeof(double),pth->error_message);
  class_alloc(ptv->used_in_output,ptv->ti_size*sizeof(int),pth->error_message);

  for (index_ti=0; index_ti<ptv->ti_size; index_ti++) {
    ptv->used_in_output[index_ti] = _TRUE_;
  }

  /** - Initialize the values of the temporary vector */

  ptv->y[ptv->index_ti_D_Tmat] = ptvs->y[ptvs->index_ti_D_Tmat];
  ptv->y[ptv->index_ti_x_H] = ptvs->y[ptvs->index_ti_x_H];
  ptv->y[ptv->index_ti_x_He] = ptvs->y[ptvs->index_ti_x_He];

  ptw->ptdw->ptv = ptv;

  /** - Evolve quantities through reionization assuming upper value of z_reio */

  z_sup = ppr->reionization_z_start_max-ppr->reionization_start_factor*pth->reionization_width;
  class_test(z_sup < 0.,
             pth->error_message,
             "parameters are such that reionization cannot take place before today while starting after z_start_max; need to increase z_start_max");

  /* maximum possible reionization redshift */
  ptw->ptrp->reionization_parameters[ptw->ptrp->index_re_reio_redshift] = z_sup;
  /* maximum possible starting redshift */
  switch (pth->reio_parametrization) {
  case reio_camb:
    ptw->ptrp->reionization_parameters[ptw->ptrp->index_re_reio_start] = ppr->reionization_z_start_max;
    break;
  case reio_half_tanh:
    ptw->ptrp->reionization_parameters[ptw->ptrp->index_re_reio_start] = z_sup;
    break;
  default:
    class_stop(pth->error_message,"Should not be there: tau_reio acan be an input only for reio_camb and reio_half_tanh");
    break;
  }

  /* ptaw->ptw->last_index_back has been properly set according to the
     redshift z = -mz_inbi, we should keep memory of it */
  last_index_back_mz_ini = ptpaw->ptw->last_index_back;

  /* Calculate a first ionization history at upper limit */
  class_call(generic_evolver(thermodynamics_derivs,
                             mz_ini,
                             mz_end,
                             ptv->y,
                             ptv->used_in_output,
                             ptv->ti_size,
                             ptpaw,
                             ppr->tol_thermo_integration,
                             ppr->smallest_allowed_variation,
                             thermodynamics_timescale,  // timescale
                             ppr->thermo_integration_stepsize, // stepsize
                             mz_output, // values of z for output
                             mz_size, // size of previous array
                             thermodynamics_sources, // function for output
                             NULL, // print variables
                             pth->error_message),
             pth->error_message,
             pth->error_message);

  /* infer corresponding tau_reio */
  class_call(thermodynamics_reionization_get_tau(ppr,
                                                 pba,
                                                 pth,
                                                 ptw),
             pth->error_message,
             pth->error_message);

  tau_sup=ptw->reionization_optical_depth;

  class_test(tau_sup < pth->tau_reio,
             pth->error_message,
             "parameters are such that reionization cannot start after z_start_max");

  /** - Restore initial conditions */

  ptv->y[ptv->index_ti_D_Tmat] = ptvs->y[ptvs->index_ti_D_Tmat];
  ptv->y[ptv->index_ti_x_H] = ptvs->y[ptvs->index_ti_x_H];
  ptv->y[ptv->index_ti_x_He] = ptvs->y[ptvs->index_ti_x_He];

  /** - Evolve quantities through reionization assuming lower value of z_reio */

  z_inf = 0.;

  /* minimum possible reionization redshift */
  ptw->ptrp->reionization_parameters[ptw->ptrp->index_re_reio_redshift] = z_inf;
  /* minimum possible starting redshift */
  switch (pth->reio_parametrization) {
  case reio_camb:
    ptw->ptrp->reionization_parameters[ptw->ptrp->index_re_reio_start] = ppr->reionization_start_factor*pth->reionization_width;
    if (ptw->ptrp->reionization_parameters[ptw->ptrp->index_re_reio_start] < pth->helium_fullreio_redshift+ppr->reionization_start_factor*pth->helium_fullreio_width) {
      ptw->ptrp->reionization_parameters[ptw->ptrp->index_re_reio_start] = pth->helium_fullreio_redshift+ppr->reionization_start_factor*pth->helium_fullreio_width;
    }
  case reio_half_tanh:
    ptw->ptrp->reionization_parameters[ptw->ptrp->index_re_reio_start] = z_inf;
    break;
  default:
    class_stop(pth->error_message,"Should not be there: tau_reio acan be an input only for reio_camb and reio_half_tanh");
    break;
  }

  /* reset ptaw->ptw->last_index_back to match the redshift z = -mz_inbi */
  ptpaw->ptw->last_index_back = last_index_back_mz_ini;

  /* Calculate a second ionization history at lower limit */
  class_call(generic_evolver(thermodynamics_derivs,
                             mz_ini,
                             mz_end,
                             ptv->y,
                             ptv->used_in_output,
                             ptv->ti_size,
                             ptpaw,
                             ppr->tol_thermo_integration,
                             ppr->smallest_allowed_variation,
                             thermodynamics_timescale,  // timescale
                             ppr->thermo_integration_stepsize, // stepsize
                             mz_output, // values of z for output
                             mz_size, // size of previous array
                             thermodynamics_sources, // function for output
                             NULL, // print variables
                             pth->error_message),
             pth->error_message,
             pth->error_message);

  /* infer corresponding tau_reio */
  class_call(thermodynamics_reionization_get_tau(ppr,
                                                 pba,
                                                 pth,
                                                 ptw),
             pth->error_message,
             pth->error_message);

  tau_inf=ptw->reionization_optical_depth;

  class_test(tau_inf > pth->tau_reio,
             pth->error_message,
             "CLASS cannot reach the low value of tau_reio that was selected, even when setting z_reio as low as 0.\nThis means that some additional physical component is requiring some minimal tau_reio_min = %.10e.\nThis is usually caused by strong energy injections or other modifications of the x_e(z) behaviour.",tau_inf);

  /** - Restore initial conditions */

  ptv->y[ptv->index_ti_D_Tmat] = ptvs->y[ptvs->index_ti_D_Tmat];
  ptv->y[ptv->index_ti_x_H] = ptvs->y[ptvs->index_ti_x_H];
  ptv->y[ptv->index_ti_x_He] = ptvs->y[ptvs->index_ti_x_He];

  /** - Evolve quantities through reionization, trying intermediate values of z_reio by bisection */

  counter=0;
  while ((tau_sup-tau_inf) > pth->tau_reio * ppr->reionization_optical_depth_tol) {
    z_mid=0.5*(z_sup+z_inf);

    /* reionization redshift */
    ptw->ptrp->reionization_parameters[ptw->ptrp->index_re_reio_redshift] = z_mid;

    /* infer starting redshift for hydrogen (Note, that this is only the start of the ADDITIONAL tanh re-ionization function)*/
    switch (pth->reio_parametrization) {
    case reio_camb:
      ptw->ptrp->reionization_parameters[ptw->ptrp->index_re_reio_start] = ptw->ptrp->reionization_parameters[ptw->ptrp->index_re_reio_redshift]+ppr->reionization_start_factor*pth->reionization_width;
      /* if starting redshift for helium is larger, take that one
       *    (does not happen in realistic models) */
      if (ptw->ptrp->reionization_parameters[ptw->ptrp->index_re_reio_start] < pth->helium_fullreio_redshift+ppr->reionization_start_factor*pth->helium_fullreio_width) {
        ptw->ptrp->reionization_parameters[ptw->ptrp->index_re_reio_start] = pth->helium_fullreio_redshift+ppr->reionization_start_factor*pth->helium_fullreio_width;
      }
      break;
    case reio_half_tanh:
      ptw->ptrp->reionization_parameters[ptw->ptrp->index_re_reio_start] = z_mid;
      break;
    default:
      class_stop(pth->error_message,"Should not be there: tau_reio acan be an input only for reio_camb and reio_half_tanh");
      break;
    }

    class_test(ptw->ptrp->reionization_parameters[ptw->ptrp->index_re_reio_start] > ppr->reionization_z_start_max,
               pth->error_message,
               "starting redshift for reionization > reionization_z_start_max = %e",ppr->reionization_z_start_max);

    /* reset ptaw->ptw->last_index_back to match the redshift z = -mz_inbi */
    ptpaw->ptw->last_index_back = last_index_back_mz_ini;

    /* Compute a new ionization history */
    class_call(generic_evolver(thermodynamics_derivs,
                               mz_ini,
                               mz_end,
                               ptv->y,
                               ptv->used_in_output,
                               ptv->ti_size,
                               ptpaw,
                               ppr->tol_thermo_integration,
                               ppr->smallest_allowed_variation,
                               thermodynamics_timescale,  // timescale
                               ppr->thermo_integration_stepsize, // stepsize
                               mz_output, // values of z for output
                               mz_size, // size of previous array
                               thermodynamics_sources, // function for output
                               NULL, // print variables
                               pth->error_message),
               pth->error_message,
               pth->error_message);

    /* infer corresponding tau_reio */
    class_call(thermodynamics_reionization_get_tau(ppr,
                                                   pba,
                                                   pth,
                                                   ptw),
               pth->error_message,
               pth->error_message);

    tau_mid=ptw->reionization_optical_depth;

    /* trial */
    if (tau_mid > pth->tau_reio) {
        z_sup=z_mid;
        tau_sup=tau_mid;
    }
    else {
        z_inf=z_mid;
        tau_inf=tau_mid;
    }

    /* Restore initial conditions */
    ptv->y[ptv->index_ti_D_Tmat] = ptvs->y[ptvs->index_ti_D_Tmat];
    ptv->y[ptv->index_ti_x_H] = ptvs->y[ptvs->index_ti_x_H];
    ptv->y[ptv->index_ti_x_He] = ptvs->y[ptvs->index_ti_x_He];

    /* counter to avoid infinite loop */
    counter++;
    class_test(counter > _MAX_IT_,
               pth->error_message,
               "while searching for reionization_optical_depth, maximum number of iterations exceeded");
  }

  /** - Store the ionization redshift in the thermodynamics structure */
  pth->z_reio = ptw->ptrp->reionization_parameters[ptw->ptrp->index_re_reio_redshift];

  /** - Free tempeoraty thermo vector */
  class_call(thermodynamics_vector_free(ptv),
             pth->error_message,
             pth->error_message);

  /* point ptw->ptdw->ptv to same location as when entering the function */
  ptw->ptdw->ptv = ptvs;

 return _SUCCESS_;
}

/**
 * Subroutine evaluating the derivative of thermodynamical quantities
 * with respect to negative redshift mz=-z.
 *
 * Automatically recognizes the current approximation interval and
 * computes the derivatives for this interval of the vector y, which
 * contains (Tmat, x_H, x_He) + others for exotic models.
 *
 * Derivatives are obtained either by calling either HyRec 2020 (Lee and Ali-Haimoud 2020, 2007.14114)
 * or RecFastCLASS (that is, RecFast version 1.5, modified by Daniel Meinert and Nils Schoeneberg for
 * better precision and smoothness at early times). See credits and licences in the wrappers (in external/...)
 *
 * This is one of the few functions in the code which are passed to the generic_evolver routine.  Since generic_evolver
 * should work with functions passed from various modules, the format of the arguments is a bit special:
 *
 * - fixed parameters and workspaces are passed through a generic pointer. Here, this pointer contains the precision, background
 *   and thermodynamics structures, plus a background vector, but generic_evolver doesn't know its precise structure.
 *
 * - the error management is a bit special: errors are not written as usual to pth->error_message, but to a generic error_message
 *   passed in the list of arguments.
 *
 * @param mz                       Input: negative redshift mz = -z
 * @param y                        Input: vector of variable to integrate
 * @param dy                       Output: its derivative (already allocated)
 * @param parameters_and_workspace Input: pointer to fixed parameters (e.g. indices) and workspace (already allocated)
 * @param error_message            Output: error message
 */

int thermodynamics_derivs(
                          double mz,
                          double * y,
                          double * dy,
                          void * parameters_and_workspace,
                          ErrorMsg error_message
                          ) {
  /** Summary: */

  /** Define local variables */
  /* Index for iterating over derivatives */
  int index_ti;

  /* Shorthand notations */
  double z,x,nH,Trad,Tmat,x_H,x_He,Hz,eps,depsdlna,dHdlna,heat_capacity,rate_gamma_b;

  /* Shorthand notations for all of the structs */
  struct thermodynamics_parameters_and_workspace * ptpaw;
  struct precision * ppr;
  struct background * pba;
  struct thermodynamics * pth;
  double * pvecback;
  struct thermo_workspace * ptw;
  struct thermo_diffeq_workspace * ptdw;
  struct thermo_vector * ptv;
  struct thermorecfast * precfast;
  struct injection * pin;
  int ap_current;

  /* Redshift */
  z = -mz;

  /** - Rename structure fields (just to avoid heavy notations) */

  /* structures */
  ptpaw = parameters_and_workspace;
  ppr = ptpaw->ppr;
  pba = ptpaw->pba;
  pth = ptpaw->pth;
  pin = &(pth->in);
  /* vector of background quantities */
  pvecback = ptpaw->pvecback;
  /* thermodynamics workspace & vector */
  ptw = ptpaw->ptw;
  ptdw = ptw->ptdw;
  ptv = ptdw->ptv;
  /* pointer to Recfast/HyRec wrappers */
  precfast = ptdw->precfast;
  /* Approximation flag */
  ap_current = ptdw->ap_current;

  /** - Get background/thermo quantities in this point */

  class_call(background_at_z(pba,
                             z,
                             long_info,
                             inter_closeby,
                             &(ptw->last_index_back),
                             pvecback),
             pba->error_message,
             error_message);

  /* Hz is H in inverse seconds (while pvecback returns [H0/c] in inverse Mpcs) */
  Hz = pvecback[pba->index_bg_H] * _c_ / _Mpc_over_m_;

  /* Total number density of Hydrogen nuclei in SI units */
  nH = ptw->SIunit_nH0 * (1.+z) * (1.+z) * (1.+z);

  /* Photon temperature in Kelvins. Modify this for some non-trivial photon temperature changes */
  Trad = ptw->Tcmb * (1.+z);

  /** Set Tmat from the evolver (it is always evolved) and store it in the workspace. */
  Tmat = y[ptv->index_ti_D_Tmat] + ptw->Tcmb*(1.+z);

  /** - The input vector y contains thermodynamic variables like
     (Tmat, x_H,x_He).  The goal of this function is: 1) Depending on
     the chosen code and current approximation, to use either
     analytical approximations or the vector y to calculate x_e; 2) To
     compute re-ionization effects on x_e; The output of this function
     is stored in the workspace ptdw */

  class_call(thermodynamics_ionization_fractions(z,y,pth,ptw,ap_current),
             pth->error_message,
             error_message);

  /* Save the output in local variables */
  x = ptdw->x_reio;

  /** - If needed, calculate heating effects (i.e. any possible energy deposition rates
        affecting the evolution equations for x and Tmat) */

  if (pth->has_exotic_injection == _TRUE_) {
    /* In case of energy injection, we currently neglect the contribution to helium ionization for RecFast ! */
    /* Note that we calculate here the energy injection INCLUDING reionization ! */
    class_call(injection_calculate_at_z(pba,pth,x,z,Tmat,pvecback),
               pin->error_message,
               error_message);
  }

  /** - Derivative of the ionization fractions */
  x_H = ptdw->x_H;
  x_He = ptdw->x_He;
  x = ptdw->x_noreio;
  switch (pth->recombination) {

    /** --> use Recfast or HyRec to get the derivatives d(x_H)/dz and
        d(x_He)/dz, and store the result directly in the vector
        dy. This gives the derivative of the ionization fractions from
        recombination only (not from reionization). Of course, the
        full treatment would involve the actual evolution equations
        for x_H and x_He during reionization, but these are not yet
        fully implemented. */
  case recfast:
    /* Hydrogen equations */
    if (ptdw->require_H == _TRUE_) {
      class_call(recfast_dx_H_dz(pth,precfast,x_H,x,nH,z,Hz,Tmat,Trad,&(dy[ptv->index_ti_x_H])),
                 precfast->error_message,
                 error_message);
    }

    /* Helium equations */
    if (ptdw->require_He == _TRUE_) {
      class_call(recfast_dx_He_dz(pth,precfast,x_He,x,x_H,nH,z,Hz,Tmat,Trad,&(dy[ptv->index_ti_x_He])),
                 precfast->error_message,
                 error_message);
    }

    break;
  case hyrec:
    /* Hydrogen equations */
    if (ptdw->require_H == _TRUE_) {
      class_call(hyrec_dx_H_dz(pth,ptw->ptdw->phyrec,x_H,x_He,x,nH,z,Hz,Tmat,Trad,&(dy[ptv->index_ti_x_H])),
                 ptw->ptdw->phyrec->error_message,
                 error_message);
    }

    /* Helium equations */
    if (ptdw->require_He == _TRUE_) {
      class_call(hyrec_dx_He_dz(pth,ptw->ptdw->phyrec,x_H,x_He,x,nH,z,Hz,Tmat,Trad,&(dy[ptv->index_ti_x_He])),
                 ptw->ptdw->phyrec->error_message,
                 error_message);
    }
    break;
  }

  /** - Derivative of the matter temperature (relevant for both Recfast and HyRec cases) */

  /* Restore the real x for the temperature equations. */
  x = ptdw->x_reio;

  /* Using the following definitions and equations, we derive a few important quantities
     Using n_e = x * n_H, n_He = f * n_H, rho_He ~ YHe * rho_b, rho_H ~ (1-YHe)*rho_b)
     - Heat capacity of the IGM
        heat_capacity = (3/2)*k_B*(n_H+n_He+n_e) = (3/2)*k_B*(1+f+x)*n_H
     - Mean baryonic molecular mass
        mu_bar = (rho_H + rho_He + rho_e)/(n_H + n_He + n_e) ~ ( (1. + YHe/(1+YHe) + 0.) * rho_H )/( (1. + f + x) * n_H) = m_H / (1+x+f) /(1-YHe)
     - Photon-baryon momentum transfer rate
        R_g = (4/3 * rho_g/rho_b) * (sigma_T * n_e) = (4/3 * rho_g * sigma_T ) * x * (1-YHe)/m_H
     - Photon-Baryon interaction rate:
        rate_gamma_b = 2 * mu_bar / m_e * R_g =  2 * (sigma_T/m_e) * (4/3 * rho_g) * x / (1. + f + x)
  */
  /* Photon-Baryon temperature change rate  */
  rate_gamma_b = ( 2. * _sigma_/_m_e_/_c_ ) * ( 4./3. * pvecback[pba->index_bg_rho_g] * _Jm3_over_Mpc2_ ) * x / (1.+x+ptw->fHe);
  /* Heat capacity of the IGM */
  heat_capacity = (3./2.)*_k_B_*nH*(1.+ptw->fHe+x);

 /*
  * A note on the temperature definition:
  *
  * All equations work with D_Tmat = Tmat - Trad
  *
  * Thus all derivatives are calculated as dD_Tmat/dz = dTmat/dz - Tcmb
  **/

 /*
  * A note on the 'early' time steady state expansion (activated here before HeIII recombination):
  *
  * Note: dTr/dz = Tr/(1+z) = Tcmb
  *
  * The early system of CMB and matter is very tightly coupled anyway, so we can expand in the following way:
  * The full equation is dTm/dz = (Tm-Tr)/e /(1+z) + 2 Tm/(1+z). Here e = H*(1+x+f)/(cT*Tr^4*x) << 1 at early times
  *
  * Find the first order solution in e, by multiplying in (1+z)*e, and approximate
  *  e*(dTm/dz)*(1+z) ~ e*(dTr/dz)*(1+z) + O(e^2) ~ e * Tr
  *
  * You find e*Tr = (Tm-Tr) + 2 Tm * e
  * Thus Tm = (1+e)/(1+2*e) * Tr = Tr * (1-e) + O(e^2)
  *
  * This is the steady state solution, which is the SAME as e.g. in HyRec
  * In our notation, eps = e*Tr, so we get Tm = Tr - eps
  *
  * So, taking the derivative of the right hand side, we obtain dTm/dz = Tcmb - eps*(dln(eps)/dz)
  *
  * Now use the form of eps = Tr*e = H*(1+x+f)/(cT*Tr^3*x) to derive the remaining terms in the below formula
  * => dln(eps)/dlna = dln(H)/dlna  - (1+f)/(1+x+f)*dln(x)/dlna + 3*dln(Tr)/dlna
  *
  * We also approximate dln(x)/dlna << 1, since we are before HeIII recombination, thus finding
  * => dln(eps)/dlna ~ dln(H)/dlna + 3
  *
  * dD_Tmat/dz = d(-eps)/dz = - eps * dln(eps)/dz = eps *dln(eps)/dlna /(1.+z)
  **/

  if ( ap_current == ptdw->index_ap_brec) {
    /* Early time steady state equation */
    dHdlna = (1.+z)*pvecback[pba->index_bg_H_prime]/pvecback[pba->index_bg_H] * _c_ / _Mpc_over_m_;
    eps =  Trad * Hz / rate_gamma_b;
    depsdlna = dHdlna/Hz + 3.;
    /* Recfast v 1.5: add here a smoothing term as suggested by Adam Moss */
    dy[ptdw->ptv->index_ti_D_Tmat] = eps * depsdlna / (1.+z);
  }

  else {
    /* Full equations at later times */
    dy[ptv->index_ti_D_Tmat] =
      + 2.*Tmat/(1.+z)                                                          /* Adiabatic expansion */
      + rate_gamma_b * (Tmat-Trad) / (Hz*(1.+z))                                /* Coupling to photons*/
      - ptw->Tcmb;                                                              /* dTrad/dz */

    if (pth->has_exotic_injection == _TRUE_) {
      dy[ptv->index_ti_D_Tmat] -= pin->pvecdeposition[pin->index_dep_heat] / heat_capacity / (Hz*(1.+z));  /* Heating from energy injection */
    }
  }

  /** - If we have extreme heatings, recombination does not fully happen
   * and/or re-ionization happens before a redshift of
   * reionization_z_start_max (default = 50).  We want to catch this
   * unphysical regime, because it would lead to further errors
   * (and/or unphysical calculations) within our recombination codes
   */

  class_test((x>1.0) && (z < ppr->z_end_reco_test) && (z > ppr->reionization_z_start_max),
             error_message,
             "At redshift %.5g : Recombination did not complete by redshift %.5g, or re-ionization happened before %.5g.\nIf this is a desired behavior, please adjust z_end_reco_test and/or reionization_z_start_max.",
             z,ppr->z_end_reco_test,ppr->reionization_z_start_max);

  /** - invert all derivatives (because the evolver evolves with -z, not with +z) */

  for(index_ti=0;index_ti<ptdw->ptv->ti_size;index_ti++) {
    dy[index_ti]=-dy[index_ti];
  }

  return _SUCCESS_;
}

/**
 * This function is relevant for the rk evolver, not ndf15. It
 * estimates a timescale 'delta z' over which quantitites vary. The rk
 * evolver divides large intervals in steps given by this timescale
 * multiplied by ppr->thermo_integration_stepsize.
 *
 * This is one of the few functions in the code which is passed to the generic_evolver routine. Since generic_evolver
 * should work with functions passed from various modules, the format of the arguments is a bit special:
 *
 * - fixed parameters and workspaces are passed through a generic pointer. generic_evolver doesn't know the content of this
 *   pointer.
 *
 * - the error management is a bit special: errors are not written as usual to pth->error_message, but to a generic error_message passed
 *   in the list of arguments.
 *
 * @param mz                              Input: minus the redshift
 * @param thermo_parameters_and_workspace Input: pointer to parameters and workspace
 * @param timescale                       Output: pointer to the timescale
 * @param error_message                   Output: possible errors are written here
 * @return the error status
 */

int thermodynamics_timescale(
                             double mz,
                             void * thermo_parameters_and_workspace,
                             double * timescale,
                             ErrorMsg error_message
                             ) {

  int index_z;
  struct thermodynamics_parameters_and_workspace * ptpaw;

  ptpaw = thermo_parameters_and_workspace;

  /* We could evaluate the timescale automatically, e.g, as [x / (dx/dz)]. */

  /* But for simplicity, we assume that the array of values of z to
     sample (pth->z_table) has been chosen in such way that the
     quantities vary only by a small amount over each step. Thus we
     define our step as delta(-z) = ((-z_i+1) - (-z_i)) = (z_i -
     z_i+1) */

  /* find index_z such that pth->z_table[index_z] > z > pth->z_table[index_z+1] */
  class_call(array_hunt_ascending(ptpaw->pth->z_table,
                                   ptpaw->pth->tt_size,
                                   -mz,
                                   &index_z,
                                   error_message),
             error_message,
             error_message);

  *timescale = ptpaw->pth->z_table[index_z+1] - ptpaw->pth->z_table[index_z];

  return _SUCCESS_;
}

/**
 * This function is passed to the generic evolver and is called whenever we want to store values for a given mz.
 *
 * The ionization fraction is either computed within a call to
 * thermodynamics_derivs(). Moreover there is an automatic smoothing
 * enabled which smoothes out the the ionization_fraction after each
 * approximation switch. This is also the place where HyRec is asked
 * to evolve x(z) using its internal system of differential equations
 * over the next range [z_i, z_i+1], and to store the result in a
 * temporary table.
 *
 * This is one of the few functions in the code which is passed to the generic_evolver routine. Since generic_evolver
 * should work with functions passed from various modules, the format of the arguments is a bit special:
 *
 * - fixed parameters and workspaces are passed through a generic pointer. generic_evolver doesn't know the content of this
 *   pointer.
 *
 * - the error management is a bit special: errors are not written as usual to pth->error_message, but to a generic error_message passed
 *   in the list of arguments.
 *
 * All quantities are computed by a simple call to thermodynamics_derivs, which computes all necessary quantities
 * and stores them in the ptdw thermo_diffeq_workspace structure
 *
 * @param mz                              Input: negative redshift, belonging to array mz_output
 * @param y                               Input: vector of evolved thermodynamical quantities
 * @param dy                              Input: derivatives of this vector w.r.t redshift
 * @param index_z                         Input: index in the array mz_output
 * @param thermo_parameters_and_workspace Input/Output: in input, all parameters needed by thermodynamics_derivs; in output, recombination table
 * @param error_message            Output: error message
 * @return the error status
 */

int thermodynamics_sources(
                           double mz,
                           double * y,
                           double * dy,
                           int index_z,
                           void * thermo_parameters_and_workspace,
                           ErrorMsg error_message
                           ) {

  /** Summary: */

  /** Define local variables */
  /* Shorthand notations */
  double z,x=0.,Tmat,Trad,dTmat;
  /* Recfast smoothing */
  double x_previous, weight,s;
  /* Structures as shorthand_notation */
  struct thermodynamics_parameters_and_workspace * ptpaw;
  struct thermodynamics * pth;
  struct thermo_workspace * ptw;
  struct thermo_diffeq_workspace * ptdw;
  struct thermo_vector * ptv;
  int ap_current;

  /* Redshift */
  z = -mz;

  /** - Rename structure fields (just to avoid heavy notations) */

  /* Structs */
  ptpaw = thermo_parameters_and_workspace;
  pth = ptpaw->pth;
  /* Thermo workspace & vector */
  ptw = ptpaw->ptw;
  ptdw = ptw->ptdw;
  ptv = ptdw->ptv;
  /* Approximation flag */
  ap_current = ptdw->ap_current;

  if (pth->has_exotic_injection == _TRUE_) {
    /* Tell heating module that it should store the heating at this z in its internal table */
    (pth->in).to_store = _TRUE_;
  }

  /** - Recalculate all quantities at this current redshift: we need
        at least pvecback, ptdw->x_reio, dy[ptv->index_ti_D_Tmat] */

  class_call(thermodynamics_derivs(mz,y,dy,thermo_parameters_and_workspace,error_message),
             error_message,
             error_message);

  /* Assign local variables (note that pvecback is filled through derivs) */
  Trad = ptw->Tcmb*(1.+z);
  Tmat = y[ptv->index_ti_D_Tmat] + Trad;
  /* Note that dy[index_ti_Q] represents dQ/d(-z), thus we need -dy here */
  dTmat = -dy[ptv->index_ti_D_Tmat] + Trad/(1.+z);

  /* get x */
  x = ptdw->x_reio;

  /** - In the recfast case, we manually smooth the results a bit */

  /* Smoothing if we are shortly after an approximation switch, i.e. if z is within 2 delta after the switch*/
  if ((ap_current != 0) && (z > ptdw->ap_z_limits[ap_current-1]-2*ptdw->ap_z_limits_delta[ap_current])) {

    class_call(thermodynamics_ionization_fractions(z,y,pth,ptw,ap_current-1),
               pth->error_message,
               error_message);

    x_previous = ptdw->x_reio;
    // get s from 0 to 1
    s = (ptdw->ap_z_limits[ap_current-1]-z)/(2*ptdw->ap_z_limits_delta[ap_current]);
    // infer f2(x) = smooth function interpolating from 0 to 1
    weight = f2(s);

    /* get smoothed x */
    x = weight*x+(1.-weight)*x_previous;
  }

  /** - Store the results in the table. Results are obtained in order of decreasing z, and stored in order of growing z */

  /* ionization fraction */
  pth->thermodynamics_table[(pth->tt_size-index_z-1)*pth->th_size+pth->index_th_xe] = x;

  /* Tb */
  pth->thermodynamics_table[(pth->tt_size-index_z-1)*pth->th_size+pth->index_th_Tb] = Tmat;

  /* Baryon temperature derivative */
  pth->thermodynamics_table[(pth->tt_size-index_z-1)*pth->th_size+pth->index_th_dTb] = dTmat;

  /* wb = (k_B/mu) Tb */
  pth->thermodynamics_table[(pth->tt_size-index_z-1)*pth->th_size+pth->index_th_wb]
    = _k_B_ / ( _c_ * _c_ * _m_H_ ) * (1. + (1./_not4_ - 1.) * ptw->YHe + x * (1.-ptw->YHe)) * Tmat;

  /* cb2 = (k_B/mu) Tb (1-1/3 dlnTb/dlna) = (k_B/mu) Tb (1 + 1/3 (1+z) dlnTb/dz) */
  pth->thermodynamics_table[(pth->tt_size-index_z-1)*pth->th_size+pth->index_th_cb2]
    = _k_B_ / ( _c_ * _c_ * _m_H_ ) * (1. + (1./_not4_ - 1.) * ptw->YHe + x * (1.-ptw->YHe)) * Tmat * (1. + (1.+z) * dTmat / Tmat / 3.);

  /* dkappa/dtau = a n_e x_e sigma_T = a^{-2} n_e(today) x_e sigma_T (in units of 1/Mpc) */
  pth->thermodynamics_table[(pth->tt_size-index_z-1)*pth->th_size+pth->index_th_dkappa]
    = (1.+z) * (1.+z) * ptw->SIunit_nH0 * x * _sigma_ * _Mpc_over_m_;

  return _SUCCESS_;
}

/**
 * Get the optical depth of reionization tau_reio for a given thermodynamical history.
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param ptw        Input: pointer to thermodynamics workspace
 * @return the error status
 */

int thermodynamics_reionization_get_tau(
                                        struct precision * ppr,
                                        struct background * pba,
                                        struct thermodynamics * pth,
                                        struct thermo_workspace * ptw
                                        ) {

  /** Summary: */

  /** Define local variables */
  /* running index inside thermodynamics table */
  int index_z,index_reio_start=0;
  double x_e_min;

  x_e_min = _HUGE_;

  /**
   * We are searching now for the start of reionization. This will be
   * the time at which the optical depth tau_reio will be computed.
   *
   * Note that the value reionization_parameters[index_reio_start] is
   * only the start of the reionization function added manually, but
   * not necessarily the total start of reionization. Reionization
   * could be longer/shifted by energy injection.
   *
   * The actual the definition of tau_reio is not unique and
   * unambiguous. We defined it here to be the optical depth up to the
   * time at which there is a global minimum in the free electron
   * fraction. We search for this time by iterating over the
   * thermodynamics table, in order to find the corresponding
   * index_reio_start.
   */

  for (index_z=0; index_z<pth->tt_size-1; index_z++) {
    if (pth->thermodynamics_table[index_z*pth->th_size+pth->index_th_xe] < x_e_min) {
      x_e_min = pth->thermodynamics_table[index_z*pth->th_size+pth->index_th_xe];
      index_reio_start = index_z;
    }
  }

  class_test(index_reio_start == pth->tt_size,
             pth->error_message,
             "reionization start = %e > largest redshift in thermodynamics table",pth->z_table[index_reio_start]);

  if (index_reio_start == 0) {
    /* the global minimum of xe(z) is at z=0. This is possible in
       models no reionization and no exotic energy
       injection. According to our definition of the reionization
       optical depth, tau_reio should then be zero. */
    ptw->reionization_optical_depth = 0;
    return _SUCCESS_;
  }
  if (index_reio_start < 3) {
    /* we need a minimum of three values to do a spline integration in the next step */
    index_reio_start =3;
  }

  /** - --> spline \f$ d \tau / dz \f$ with respect to z in view of
            integrating for optical depth between 0 and the just found
            starting index */
  class_call(array_spline_table_line_to_line(pth->tau_table,
                                             index_reio_start,
                                             pth->thermodynamics_table,
                                             pth->th_size,
                                             pth->index_th_dkappa,
                                             pth->index_th_dddkappa,
                                             _SPLINE_EST_DERIV_,
                                             pth->error_message),
             pth->error_message,
             pth->error_message);

  /** - --> integrate for optical depth */
  class_call(array_integrate_all_spline_table_line_to_line(pth->tau_table,
                                                           index_reio_start,
                                                           pth->thermodynamics_table,
                                                           pth->th_size,
                                                           pth->index_th_dkappa,
                                                           pth->index_th_dddkappa,
                                                           &(ptw->reionization_optical_depth),
                                                           pth->error_message),
             pth->error_message,
             pth->error_message);

  ptw->reionization_optical_depth *= -1; // tau and z go in reverse order, so we must flip the sign

  return _SUCCESS_;
}

/**
 * Free the thermo_vector structure, which is the '->ptv' field of the thermodynamics_differential_workspace ptdw structure
 *
 * @param tv        Input: pointer to thermo_vector structure to be freed
 * @return the error status
 */

int thermodynamics_vector_free(
                               struct thermo_vector * tv
                               ) {

  free(tv->y);
  free(tv->dy);
  free(tv->used_in_output);
  free(tv);

  return _SUCCESS_;
}

/**
 * Compute the baryon drag conformal time tau_d = [int_{tau_today}^{tau} dtau -dkappa_d/dtau]
 *
 * @param pba                Input: pointer to background structure
 * @param pth                Input/Output: pointer to initialized thermodynamics structure
 * @param pvecback           Input: Initialized vector of background quantities
 * @return the error status
 */

int thermodynamics_calculate_conformal_drag_time(
                                                 struct background* pba,
                                                 struct thermodynamics* pth,
                                                 double* pvecback
                                                 ) {

  /** Summary: */

  /** Define local variables */
  double R;
  int index_tau;
  int last_index_back;

  /** - compute minus the baryon drag interaction rate time, -dkappa_d/dtau = -[1/R * kappa'], with R = 3 rho_b / 4 rho_gamma,
        stored temporarily in column ddkappa */

  /* find the value of last_index_back for tau_table[0] in order to speed up subsequent interpolations in the loop */
  class_call(background_at_tau(pba,
                               pth->tau_table[0],
                               normal_info,
                               inter_normal,
                               &last_index_back,
                               pvecback),
             pba->error_message,
             pth->error_message);

  for (index_tau=0; index_tau < pth->tt_size; index_tau++) {

    class_call(background_at_tau(pba,
                                 pth->tau_table[index_tau],
                                 normal_info,
                                 inter_closeby,
                                 &last_index_back,
                                 pvecback),
               pba->error_message,
               pth->error_message);

    R = 3./4.*pvecback[pba->index_bg_rho_b]/pvecback[pba->index_bg_rho_g];

    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_ddkappa] =
      -1./R*pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa];

  }

  /** - compute second derivative of this rate, -[1/R * kappa']'', stored temporarily in column dddkappa */
  class_call(array_spline_table_line_to_line(pth->tau_table,
                                             pth->tt_size,
                                             pth->thermodynamics_table,
                                             pth->th_size,
                                             pth->index_th_ddkappa,
                                             pth->index_th_dddkappa,
                                             _SPLINE_EST_DERIV_,
                                             pth->error_message),
             pth->error_message,
             pth->error_message);

  /** - compute tau_d = [int_{tau_today}^{tau} dtau -dkappa_d/dtau] */
  class_call(array_integrate_spline_table_line_to_line(pth->tau_table,
                                                       pth->tt_size,
                                                       pth->thermodynamics_table,
                                                       pth->th_size,
                                                       pth->index_th_ddkappa,
                                                       pth->index_th_dddkappa,
                                                       pth->index_th_tau_d,
                                                       pth->error_message),
             pth->error_message,
             pth->error_message);

  return _SUCCESS_;
}

/**
 * Compute the damping scale
 *   r_d = 2pi/k_d = 2pi * [int_{tau_ini}^{tau} dtau (1/kappa') 1/6 (R^2+16/15(1+R))/(1+R)^2]^1/2
 *                 =  2pi * [int_{tau_ini}^{tau} dtau (1/kappa') 1/6 (R^2/(1+R)+16/15)/(1+R)]^1/2
 *
 * which is like in CosmoTherm (CT), but slightly different from Wayne Hu (WH)'s thesis eq. (5.59):
 * The factor 16/15 in CT is 4/5 in WH, but 16/15 is taking more effects into account
 *
 * @param pba                Input: pointer to background structure
 * @param pth                Input/Output: pointer to initialized thermodynamics structure
 * @param pvecback           Input: Initialized vector of background quantities
 * @return the error status
 */

int thermodynamics_calculate_damping_scale(
                                           struct background* pba,
                                           struct thermodynamics* pth,
                                           double* pvecback
                                           ) {

  /** Summary: */

  /** Define local variables */
  double R;
  /* Initial time and dkappa/dtau */
  double tau_ini,dkappa_ini;
  double* tau_table_growing;
  int index_tau;
  int last_index_back;

  class_alloc(tau_table_growing,pth->tt_size*sizeof(double),pth->error_message);

  /* find the value of last_index_back for
     tau_table[pth->tt_size-1] = tau_table_growing[0] in order to
     speed up subsequent interpolations in the loop */
  class_call(background_at_tau(pba,
                               pth->tau_table[pth->tt_size-1],
                               normal_info,
                               inter_normal,
                               &last_index_back,
                               pvecback),
             pba->error_message,
             pth->error_message);

  /* compute integrand and store temporarily in column "ddkappa" */
  for (index_tau=0; index_tau < pth->tt_size; index_tau++) {

    tau_table_growing[index_tau] = pth->tau_table[pth->tt_size-1-index_tau];

    class_call(background_at_tau(pba,
                                 tau_table_growing[index_tau],
                                 normal_info,
                                 inter_closeby,
                                 &last_index_back,
                                   pvecback),
               pba->error_message,
               pth->error_message);

    R = 3./4.*pvecback[pba->index_bg_rho_b]/pvecback[pba->index_bg_rho_g];

    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_ddkappa] =
      1./6./pth->thermodynamics_table[(pth->tt_size-1-index_tau)*pth->th_size+pth->index_th_dkappa]
      *(R*R/(1+R)+16./15.)/(1.+R);

  }

  /* compute second derivative of integrand, and store temporarily in column "dddkappa" */
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

  /* compute integratal and store temporarily in column "g" */
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

  /* we could now write the result as r_d = 2pi * sqrt(integral),
     but we will first better acount for the contribution frokm the tau_ini boundary.
     Close to this boundary, R=0 and the integrand is just 16/(15*6)/kappa'
     Using kappa' propto 1/a^2 and tau propro a during RD, we get the analytic result:
     int_0^{tau_ini} dtau / kappa' = tau_ini / 3 / kappa'_ini
     Thus r_d = 2pi * sqrt( 16/(15*6*3) * (tau_ini/ kappa'_ini) * integral) */

  tau_ini = pth->tau_table[pth->tt_size-1];
  dkappa_ini = pth->thermodynamics_table[(pth->tt_size-1)*pth->th_size+pth->index_th_dkappa];

  for (index_tau=0; index_tau < pth->tt_size; index_tau++) {

    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_r_d] =
      2.*_PI_*sqrt(16./(15.*6.*3.)*tau_ini/dkappa_ini
                   +pth->thermodynamics_table[(pth->tt_size-1-index_tau)*pth->th_size+pth->index_th_g]);

  }

  return _SUCCESS_;
}

/**
 * Calculate quantities relating to optical phenomena like kappa' and exp(-kappa) and the visibility function,
 * optical depth, etc.
 *
 * @param ppr   Input: pointer to precision structure
 * @param pth   Input/Output: pointer to thermodynamics structure
 * @return the error status
 */

int thermodynamics_calculate_opticals(
                                      struct precision* ppr,
                                      struct thermodynamics* pth
                                      ) {

  /** Summary: */

  /** Define local quantities */
  /* Visibility function value */
  double g;
  /* kappa derivative values*/
  double dkappa,ddkappa,dddkappa,expmkappa;
  int index_tau;

  /** - --> second derivative with respect to tau of dkappa (in view of spline interpolation) */
  class_call(array_spline_table_line_to_line(pth->tau_table,
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
  class_call(array_derive_spline_table_line_to_line(pth->tau_table,
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
  class_call(array_integrate_spline_table_line_to_line(pth->tau_table,
                                                       pth->tt_size,
                                                       pth->thermodynamics_table,
                                                       pth->th_size,
                                                       pth->index_th_dkappa,
                                                       pth->index_th_dddkappa,
                                                       pth->index_th_g,
                                                       pth->error_message),
             pth->error_message,
             pth->error_message);

  /** - --> compute visibility: \f$ g= (d \kappa/d \tau) e^{- \kappa} \f$ */

  /* loop on z (decreasing z, increasing time) */
  for (index_tau=pth->tt_size-1; index_tau>=0; index_tau--) {

    dkappa = pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa];
    ddkappa = pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_ddkappa];
    dddkappa = pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dddkappa];
    expmkappa = exp(pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g]);

    /** - ---> compute g */
    g = dkappa * expmkappa;

    /* for some very extreme models, in the last line, the exponential of a large negative number could go beyond the range covered by the "double" representation numbers, and be set to zero. To avoid a division by zero in the next steps, it is then better to set it to the minimum non-zero double (this has no impact on observables). */
    if (g==0.) g=DBL_MIN;

    /** - ---> compute exp(-kappa) */
    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_exp_m_kappa] = expmkappa;

    /** - ---> compute g' (the plus sign of the second term is correct, see def of -kappa in thermodynamics module!) */
    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dg] =
      (ddkappa + dkappa * dkappa) * expmkappa;

    /** - ---> compute g''  */
    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_ddg] =
      (dddkappa + dkappa * ddkappa * 3. + dkappa * dkappa * dkappa ) * expmkappa;

    /** - ---> store g */
    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g] = g;

    /** - ---> compute variation rate */
    class_test(pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa] == 0.,
               pth->error_message,
               "variation rate diverges");

    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_rate] =
      sqrt(pow(dkappa,2)+pow(ddkappa/dkappa,2)+fabs(dddkappa/dkappa));

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

  /** - --> derivatives of baryon sound speed (only computed if some non-minimal tight-coupling schemes is requested) */
  if (pth->compute_cb2_derivatives == _TRUE_) {

    /** - ---> second derivative with respect to tau of cb2 */
    class_call(array_spline_table_line_to_line(pth->tau_table,
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
    class_call(array_derive_spline_table_line_to_line(pth->tau_table,
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

  return _SUCCESS_;
}

/**
 * Compute the idm_dr quantities: idm-dr opacities, idr self-interaction rate,
 * dark optical depths, etc.
 *
 * @param ppr                Input: pointer to precision structure
 * @param pba                Input: pointer to background structure
 * @param pth                Input/Output: pointer to initialized thermodynamics structure
 * @param pvecback           Input: Initialized vector of background quantities
 * @return the error status
 */

int thermodynamics_calculate_idm_dr_quantities(
                                               struct precision* ppr,
                                               struct background * pba,
                                               struct thermodynamics* pth,
                                               double* pvecback
                                               ) {

  /** Summary: */

  /** - Define local variables */
  double z_idm_dr, z_idr, tau_idm_dr, tau_idr, Gamma_heat_idm_dr, dTdz_idm_dr, T_idm_dr, z, T_idr, dz, T_adia, z_adia;
  int index_tau_fs;
  int n, N_sub_steps;
  double dz_sub_step;
  int index_tau;
  double tau;
  int last_index_back;

  /* find the value of last_index_back for tau_table[0] in order to speed up subsequent interpolations in the loop */
  class_call(background_at_tau(pba,
                               pth->tau_table[0],
                               normal_info,
                               inter_normal,
                               &last_index_back,
                               pvecback),
             pba->error_message,
             pth->error_message);

  for (index_tau=0; index_tau < pth->tt_size; index_tau++) {

    class_call(background_at_tau(pba,
                                 pth->tau_table[index_tau],
                                 normal_info,
                                 inter_closeby,
                                 &last_index_back,
                                 pvecback),
               pba->error_message,
               pth->error_message);

    /* - idr interaction rate with idm_dr (i.e. idr opacity to idm_dr scattering) */
    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dmu_idm_dr] =
      pth->a_idm_dr*pow((1.+pth->z_table[index_tau])/1.e7,pth->nindex_idm_dr)*pba->Omega0_idm_dr*pow(pba->h,2);

    /* - idm_dr interaction rate with idr (i.e. idm_dr opacity
       to idr scattering), [Sinv*dmu_idm_dr] with Sinv = (4
       rho_idr) / (3 rho_idm_dr), stored temporarily in
       ddmu_idm_dr */
    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_ddmu_idm_dr] =
      4./3.*pvecback[pba->index_bg_rho_idr]/pvecback[pba->index_bg_rho_idm_dr]
      *pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dmu_idm_dr];

    /* - idr self-interaction rate */
    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dmu_idr] =
      pth->b_idr*pow((1.+pth->z_table[index_tau])/1.e7,pth->nindex_idm_dr)*pba->Omega0_idr*pow(pba->h,2);
  }

  /** - second derivative of idm_dr interaction rate (with idr), [Sinv*dmu_idm_dr]'', stored temporarily in column dddmu */
  class_call(array_spline_table_line_to_line(pth->tau_table,
                                             pth->tt_size,
                                             pth->thermodynamics_table,
                                             pth->th_size,
                                             pth->index_th_ddmu_idm_dr,
                                             pth->index_th_dddmu_idm_dr,
                                             _SPLINE_EST_DERIV_,
                                             pth->error_message),
             pth->error_message,
             pth->error_message);

  /** - compute optical depth of idm, tau_idm_dr = [int_{tau_today}^{tau} dtau [Sinv*dmu_idm_dr] ].
      This step gives -tau_idm_dr. The resulty is mutiplied by -1 later on. */
  class_call(array_integrate_spline_table_line_to_line(pth->tau_table,
                                                       pth->tt_size,
                                                       pth->thermodynamics_table,
                                                       pth->th_size,
                                                       pth->index_th_ddmu_idm_dr,
                                                       pth->index_th_dddmu_idm_dr,
                                                       pth->index_th_tau_idm_dr,
                                                       pth->error_message),
             pth->error_message,
             pth->error_message);


  /** - second derivative of idr interaction rate (with idm_dr), [dmu_idm_idr]'', stored temporarily in column dddmu */
  class_call(array_spline_table_line_to_line(pth->tau_table,
                                             pth->tt_size,
                                             pth->thermodynamics_table,
                                             pth->th_size,
                                             pth->index_th_dmu_idm_dr,
                                             pth->index_th_dddmu_idm_dr,
                                             _SPLINE_EST_DERIV_,
                                             pth->error_message),
             pth->error_message,
             pth->error_message);

  /** - compute optical depth of idr, tau_idr = [int_{tau_today}^{tau} dtau [dmu_idm_idr] ].
      This step gives -tau_idr. The resulty is mutiplied by -1 later on. */
  class_call(array_integrate_spline_table_line_to_line(pth->tau_table,
                                                       pth->tt_size,
                                                       pth->thermodynamics_table,
                                                       pth->th_size,
                                                       pth->index_th_dmu_idm_dr,
                                                       pth->index_th_dddmu_idm_dr,
                                                       pth->index_th_tau_idr,
                                                       pth->error_message),
             pth->error_message,
             pth->error_message);


  /* - restore correct sign for idm_dr and idr optical depth, and calculate idm_dr visibility function */
  /* loop on z (decreasing z, increasing time) */
  for (index_tau=pth->tt_size-1; index_tau>=0; index_tau--) {
    /* - --> restore the correct sign for tau_idm_dr */
    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_tau_idm_dr] *= -1.;

    /* - --> restore the correct sign for tau_idr */
    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_tau_idr] *= -1.;

    /* - --> visibility function for idm_dr : g_idm_dr = [Sinv*dmu_idm_dr] * exp(-tau_idm_dr) */
    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g_idm_dr] =
      pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_ddmu_idm_dr]
      * exp(-pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_tau_idm_dr]);
  }

  /* - fill columns for ddmu_idm_dr and dddmu_idm_dr with true values, and compute idm_dr temperature and sound speed */

  /** - --> second derivative with respect to tau of dmu_idm_dr (in view of spline interpolation) */
  class_call(array_spline_table_line_to_line(pth->tau_table,
                                             pth->tt_size,
                                             pth->thermodynamics_table,
                                             pth->th_size,
                                             pth->index_th_dmu_idm_dr,
                                             pth->index_th_dddmu_idm_dr,
                                             _SPLINE_EST_DERIV_,
                                             pth->error_message),
             pth->error_message,
             pth->error_message);

  /** - --> first derivative with respect to tau of dmu_idm_dr (using spline interpolation) */
  class_call(array_derive_spline_table_line_to_line(pth->tau_table,
                                                    pth->tt_size,
                                                    pth->thermodynamics_table,
                                                    pth->th_size,
                                                    pth->index_th_dmu_idm_dr,
                                                    pth->index_th_dddmu_idm_dr,
                                                    pth->index_th_ddmu_idm_dr,
                                                    pth->error_message),
             pth->error_message,
             pth->error_message);

  /** - --> now compute idm_dr temperature and sound speed in various regimes */

  /* (A) - initial value of T_idm_dr at the maximum z (minimum tau) */

  z = pth->z_table[pth->tt_size-1];

  class_call(background_at_z(pba, z, short_info, inter_normal, &last_index_back, pvecback),
             pba->error_message,
             pth->error_message);

  Gamma_heat_idm_dr = 2.*pba->Omega0_idr*pow(pba->h,2)*pth->a_idm_dr*pow((1.+z),(pth->nindex_idm_dr+1.))/pow(1.e7,pth->nindex_idm_dr);

  /* (A1) --> if Gamma is not much smaller than H, set T_idm_dr to T_idm_dr = T_idr = xi*T_gamma (tight coupling solution) */
  if (Gamma_heat_idm_dr > 1.e-3 * pvecback[pba->index_bg_a]*pvecback[pba->index_bg_H]) {
    T_idm_dr = pba->T_idr*(1.+z);
    dTdz_idm_dr = pba->T_idr;
  }

  /* (A2) --> otherwise, if Gamma << H, set initial T_idm_dr to the
     approximate analytic solution (Gamma/aH)/(1+(Gamma/aH)*T_idr)
     (eq. (A62) in ETHOS I ) */
  else {
    T_idr = pba->T_idr*(1.+z);
    T_idm_dr = Gamma_heat_idm_dr/(pvecback[pba->index_bg_a]*pvecback[pba->index_bg_H])
      /(1. + Gamma_heat_idm_dr/(pvecback[pba->index_bg_a]*pvecback[pba->index_bg_H]))*T_idr;
    dTdz_idm_dr = 2.*T_idm_dr - Gamma_heat_idm_dr/pvecback[pba->index_bg_H] * (T_idr - T_idm_dr);
  }

  pth->thermodynamics_table[(pth->tt_size-1)*pth->th_size+pth->index_th_Tidm_dr] = T_idm_dr;
  pth->thermodynamics_table[(pth->tt_size-1)*pth->th_size+pth->index_th_cidm_dr2] = _k_B_*T_idm_dr/_eV_/pth->m_idm_dr*(1.+dTdz_idm_dr/3./T_idm_dr);

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
  for (index_tau=pth->tt_size-2; index_tau>=0; index_tau--) {

    /* (B1) --> tight-coupling solution: Gamma >> H implies T_idm_dr=T_idr=xi*T_gamma */
    if (Gamma_heat_idm_dr > 1.e3 * pvecback[pba->index_bg_a]*pvecback[pba->index_bg_H]) {
      z = pth->z_table[index_tau];
      T_idr = pba->T_idr*(1.+z);
      T_idm_dr = T_idr;
      Gamma_heat_idm_dr = 2.*pba->Omega0_idr*pow(pba->h,2)*pth->a_idm_dr*pow((1.+z),(pth->nindex_idm_dr+1.))/pow(1.e7,pth->nindex_idm_dr);
      class_call(background_at_z(pba, z, short_info, inter_normal, &last_index_back, pvecback),
                 pba->error_message,
                 pth->error_message);
      dTdz_idm_dr =pba->T_idr;
    }

    /* (B2) --> intermediate solution: integrate differential equation equation dT_idm_dr/dz = 2 a T_DM - Gamma/H (T_idr - T_idm_dr) */
    else if (Gamma_heat_idm_dr > 1.e-3 * pvecback[pba->index_bg_a]*pvecback[pba->index_bg_H]) {

      dz = pth->z_table[index_tau+1] - pth->z_table[index_tau];

      /* (B2a) ----> if dz << H/Gamma the equation is not too stiff and the traditional forward Euler method converges */
      if (dz < pvecback[pba->index_bg_H]/Gamma_heat_idm_dr/10.) {
        z = pth->z_table[index_tau];
        T_idr = pba->T_idr*(1.+z);
        T_idm_dr -= dTdz_idm_dr*dz;
        Gamma_heat_idm_dr = 2.*pba->Omega0_idr*pow(pba->h,2)*pth->a_idm_dr*pow((1.+z),(pth->nindex_idm_dr+1.))/pow(1.e7,pth->nindex_idm_dr);
        class_call(background_at_z(pba, z, short_info, inter_normal, &last_index_back, pvecback),
                   pba->error_message,
                   pth->error_message);
        dTdz_idm_dr = 2.*pvecback[pba->index_bg_a]*T_idm_dr-Gamma_heat_idm_dr/(pvecback[pba->index_bg_H])*(T_idr-T_idm_dr);
      }

      /* (B2b) ----> otherwise, the equation is too stiff and the
         traditional forward Euler method diverges with this
         stepsize. But we can just decreasee dz to bring it back
         well within the convergence radius H/Gamma of the
         equation. */
      else {
        N_sub_steps = (int)(dz/ (pvecback[pba->index_bg_H]/Gamma_heat_idm_dr/10.))+1;
        dz_sub_step = dz/N_sub_steps;

        /* loop over sub-steps */
        for (n=0; n<N_sub_steps; n++) {

          /* evolve quantities over  sub-step with forward Euler method */

          z -= dz_sub_step;
          /* final redshift last sub-step overwritten to avoid small rounding error */
          if (n==(N_sub_steps-1)) z=pth->z_table[index_tau];

          T_idr = pba->T_idr*(1.+z);
          T_idm_dr -= dTdz_idm_dr*dz_sub_step;
          Gamma_heat_idm_dr = 2.*pba->Omega0_idr*pow(pba->h,2)*pth->a_idm_dr*pow((1.+z),(pth->nindex_idm_dr+1.))/pow(1.e7,pth->nindex_idm_dr);
          class_call(background_at_z(pba, z, short_info, inter_normal, &last_index_back, pvecback),
                     pba->error_message,
                     pth->error_message);
          dTdz_idm_dr = 2.*pvecback[pba->index_bg_a]*T_idm_dr-Gamma_heat_idm_dr/(pvecback[pba->index_bg_H])*(T_idr-T_idm_dr);
        }
      }

      /* update T_adia, z_adia */
      T_adia = T_idm_dr;
      z_adia = z;
    }

    /* (B3) --> decoupled solution: T_idm_dr scales like a^-2 */
    else {
      z = pth->z_table[index_tau];
      T_idr = pba->T_idr*(1.+z);
      T_idm_dr = T_adia * pow((1.+z)/(1.+z_adia),2);
      Gamma_heat_idm_dr = 2.*pba->Omega0_idr*pow(pba->h,2)*pth->a_idm_dr*pow((1.+z),(pth->nindex_idm_dr+1.))/pow(1.e7,pth->nindex_idm_dr);
      class_call(background_at_z(pba, z, short_info, inter_normal, &last_index_back, pvecback),
                 pba->error_message,
                 pth->error_message);
      dTdz_idm_dr = 2./(1+z)*T_idm_dr;
    }

    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_Tidm_dr] = T_idm_dr;
    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_cidm_dr2] = _k_B_*T_idm_dr/_eV_/pth->m_idm_dr*(1.+dTdz_idm_dr/3./T_idm_dr);
  }

  /** - Find interacting dark radiation free-streaming time */
  index_tau_fs = index_tau;

  if (pba->has_idm_dr == _TRUE_) {

    if (pth->nindex_idm_dr>=2) {
      index_tau=index_tau_fs-1;
      /* comment: using index_tau_max (index_tau_fs) instead of pth->tt_size-1 ensures that the switch is always after recombination (free streaming) */
    }
    else{
      index_tau=0;
    }

    class_call(background_tau_of_z(pba,pth->z_table[index_tau],&tau),
               pba->error_message,
               pth->error_message);

    while ((1./pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_dmu_idm_dr]/tau
            < ppr->idr_streaming_trigger_tau_c_over_tau) &&
           ((pth->nindex_idm_dr >= 2 && index_tau > 0) ||
            (pth->nindex_idm_dr < 2 && index_tau < pth->tt_size-1))) {

      if (pth->nindex_idm_dr>=2) {
        index_tau--;
      }
      else{
        index_tau++;
      }

      class_call(background_tau_of_z(pba,pth->z_table[index_tau],&tau),
                 pba->error_message,
                 pth->error_message);

      }

    pth->tau_idr_free_streaming = tau;
  }

  /* case of idr alone without idm_dr */
  else {
    index_tau=index_tau_fs-1;
    class_call(background_tau_of_z(pba,pth->z_table[index_tau],&tau),
               pba->error_message,
               pth->error_message)

    pth->tau_idr_free_streaming = tau;
  }

  /** - find idm_dr and idr drag times */
  if ((pth->thermodynamics_table[(pth->tt_size-1)*pth->th_size+pth->index_th_tau_idm_dr]>1.) && (pth->thermodynamics_table[(pth->tt_size-1)*pth->th_size+pth->index_th_tau_idr]>1.)) {
    index_tau=0;

    while ((pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_tau_idm_dr] < 1.) && (index_tau < pth->tt_size-1))
      index_tau++;

    z_idm_dr = pth->z_table[index_tau-1]+(1.-pth->thermodynamics_table[(index_tau-1)*pth->th_size+pth->index_th_tau_idm_dr])
      /(pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_tau_idm_dr]-pth->thermodynamics_table[(index_tau-1)*pth->th_size+pth->index_th_tau_idm_dr])
      *(pth->z_table[index_tau]-pth->z_table[index_tau-1]);

    class_call(background_tau_of_z(pba,z_idm_dr,&(tau_idm_dr)),
               pba->error_message,
               pth->error_message);

    index_tau=0;

    while ((pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_tau_idr] < 1.) && (index_tau < pth->tt_size-1))
      index_tau++;

    z_idr = pth->z_table[index_tau-1]+
      (1.-pth->thermodynamics_table[(index_tau-1)*pth->th_size+pth->index_th_tau_idr])
      /(pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_tau_idr]-pth->thermodynamics_table[(index_tau-1)*pth->th_size+pth->index_th_tau_idr])
      *(pth->z_table[index_tau]-pth->z_table[index_tau-1]);

    class_call(background_tau_of_z(pba,z_idr,&(tau_idr)),
               pba->error_message,
               pth->error_message);
  }
  /* - store final quantities */

  pth->tau_idr = tau_idr;
  pth->tau_idm_dr = tau_idm_dr;

  return _SUCCESS_;
}

/**
 * Calculate various quantities at the time of recombination, as well
 * as the time tau_cut at which visibility gets negligible and one can
 * assume pure free-streaming.
 *
 * @param ppr                Input: pointer to precision structure
 * @param pba                Input: pointer to background structure
 * @param pth                Input/Output: pointer to initialized thermodynamics structure
 * @param pvecback           Input: pointer to some allocated pvecback
 * @return the error status
 */

int thermodynamics_calculate_recombination_quantities(
                                                      struct precision* ppr,
                                                      struct background * pba,
                                                      struct thermodynamics* pth,
                                                      double* pvecback
                                                      ) {

  /** Summary: */

  /** Define local variables */
  double g_max;

  int index_tau;
  int index_tau_max;
  int last_index_back=0;

  double tau;

  /** - find maximum of g */
  index_tau=pth->tt_size-1;
  while (pth->z_table[index_tau]>_Z_REC_MAX_) {
    index_tau--;
  }

  class_test(pth->thermodynamics_table[(index_tau+1)*pth->th_size+pth->index_th_g] >
             pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g],
             pth->error_message,
             "The visibility function is not increasing at redshift _Z_REC_MAX_=%g, which is the value imposed in thermodynamics.h\n This implies that recombination must have already happened at a too early time.",_Z_REC_MAX_);

  while (pth->thermodynamics_table[(index_tau+1)*pth->th_size+pth->index_th_g] <=
         pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g]) {
    index_tau--;
  }

  g_max = pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g];
  index_tau_max = index_tau;

  /* approximation for maximum of g, using cubic interpolation, assuming equally spaced z's */
  pth->z_rec=pth->z_table[index_tau+1]+0.5*(pth->z_table[index_tau+1]-pth->z_table[index_tau])
                                          *(pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_g]
                                                   -pth->thermodynamics_table[(index_tau+2)*pth->th_size+pth->index_th_g])
                                           /(pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_g]
                                                   -2.*pth->thermodynamics_table[(index_tau+1)*pth->th_size+pth->index_th_g]
                                                   +pth->thermodynamics_table[(index_tau+2)*pth->th_size+pth->index_th_g]
                                            );

  class_test(pth->z_rec+ppr->smallest_allowed_variation >= _Z_REC_MAX_,
             pth->error_message,
             "recombination (at z=%g) happens before _Z_REC_MAX_=%g, which is the maximum value imposed in thermodynamics.h",pth->z_rec+ppr->smallest_allowed_variation,_Z_REC_MAX_);

  class_test(pth->z_rec-ppr->smallest_allowed_variation <= _Z_REC_MIN_,
             pth->error_message,
             "recombination (at z=%g) happens after _Z_REC_MIN_=%g, which is the minimum value imposed in thermodynamics.h",pth->z_rec-ppr->smallest_allowed_variation,_Z_REC_MIN_);

  /** - find conformal recombination time using background_tau_of_z **/

  class_call(background_tau_of_z(pba,pth->z_rec,&(pth->tau_rec)),
             pba->error_message,
             pth->error_message);

  class_call(background_at_z(pba,pth->z_rec, long_info, inter_normal, &last_index_back, pvecback),
             pba->error_message,
             pth->error_message);

  pth->rs_rec=pvecback[pba->index_bg_rs];
  pth->ds_rec=pth->rs_rec/(1.+pth->z_rec);
  pth->da_rec=pvecback[pba->index_bg_ang_distance];
  pth->ra_rec=pth->da_rec*(1.+pth->z_rec);
  pth->angular_rescaling=pth->ra_rec/(pba->conformal_age-pth->tau_rec);

  /** - find damping scale at recombination (using linear interpolation) */
  if (pth->compute_damping_scale == _TRUE_) {

    pth->rd_rec = (pth->z_table[index_tau+1]-pth->z_rec)/(pth->z_table[index_tau+1]-pth->z_table[index_tau])*pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_r_d]
      +(pth->z_rec-pth->z_table[index_tau])/(pth->z_table[index_tau+1]-pth->z_table[index_tau])*pth->thermodynamics_table[(index_tau+1)*pth->th_size+pth->index_th_r_d];

  }

  /** - find time (always after recombination) at which tau_c/tau falls below some threshold, defining tau_free_streaming */
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

  /** - find time above which visibility falls below a given fraction of its maximum */
  index_tau=index_tau_max;
  while ((pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_g] >
          g_max * ppr->neglect_CMB_sources_below_visibility)
         && (index_tau > 0))
    index_tau--;

  class_call(background_tau_of_z(pba,pth->z_table[index_tau],&(pth->tau_cut)),
             pba->error_message,
             pth->error_message);

  /** - find z_star (when optical depth kappa crosses one, using linear
      interpolation) and sound horizon at that time */

  index_tau=0;
  while ((pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_exp_m_kappa] > 1./_E_) && (index_tau < pth->tt_size))
    index_tau++;

  pth->z_star = pth->z_table[index_tau-1]+
    (1./_E_-pth->thermodynamics_table[(index_tau-1)*pth->th_size+pth->index_th_exp_m_kappa])
    /(pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_exp_m_kappa]-pth->thermodynamics_table[(index_tau-1)*pth->th_size+pth->index_th_exp_m_kappa])
    *(pth->z_table[index_tau]-pth->z_table[index_tau-1]);

  class_call(background_tau_of_z(pba,pth->z_star,&(pth->tau_star)),
             pba->error_message,
             pth->error_message);

  class_call(background_at_tau(pba,pth->tau_star, long_info, inter_normal, &last_index_back, pvecback),
             pba->error_message,
             pth->error_message);

  pth->rs_star=pvecback[pba->index_bg_rs];
  pth->ds_star=pth->rs_star/(1.+pth->z_star);
  pth->da_star=pvecback[pba->index_bg_ang_distance];
  pth->ra_star=pth->da_star*(1.+pth->z_star);

  if (pth->compute_damping_scale == _TRUE_) {

    pth->rd_star = (pth->z_table[index_tau+1]-pth->z_star)/(pth->z_table[index_tau+1]-pth->z_table[index_tau])*pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_r_d]
      +(pth->z_star-pth->z_table[index_tau])/(pth->z_table[index_tau+1]-pth->z_table[index_tau])*pth->thermodynamics_table[(index_tau+1)*pth->th_size+pth->index_th_r_d];

  }

  return _SUCCESS_;
}

/**
 * Calculate various quantities at the time of ending of baryon drag (It is precisely where tau_d crosses one)
 *
 * @param ppr                Input: pointer to precision structure
 * @param pba                Input: pointer to background structure
 * @param pth                Input/Output: pointer to initialized thermodynamics structure
 * @param pvecback           Input: pointer to some allocated pvecback
 * @return the error status
 */

int thermodynamics_calculate_drag_quantities(
                                             struct precision* ppr,
                                             struct background * pba,
                                             struct thermodynamics* pth,
                                             double* pvecback
                                             ) {

  /** Summary: */

  /** Define local variables */
  int index_tau;
  int last_index_back;

  /** - find baryon drag time (when tau_d crosses one, using linear interpolation) and sound horizon at that time */
  index_tau=0;
  while ((pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_tau_d] < 1.) && (index_tau < pth->tt_size)) {
    index_tau++;
  }

  pth->z_d = pth->z_table[index_tau-1]+
    (1.-pth->thermodynamics_table[(index_tau-1)*pth->th_size+pth->index_th_tau_d])
    /(pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_tau_d]-pth->thermodynamics_table[(index_tau-1)*pth->th_size+pth->index_th_tau_d])
    *(pth->z_table[index_tau]-pth->z_table[index_tau-1]);

  class_call(background_tau_of_z(pba,pth->z_d,&(pth->tau_d)),
             pba->error_message,
             pth->error_message);

  class_call(background_at_z(pba,pth->z_d, long_info, inter_normal, &last_index_back, pvecback),
             pba->error_message,
             pth->error_message);

  pth->rs_d=pvecback[pba->index_bg_rs];
  pth->ds_d=pth->rs_d/(1.+pth->z_d);

  return _SUCCESS_;
}

/**
 * Compute ionization fractions with the RecFast of HyRec algorithm.
 *
 * Compute ionization fractions using either the vector y or, in some
 * approximation schemes, some analytic approximations. The output of
 * this function is located in the ptw->ptdw workspace. We need to
 * assign:
 *
 * - in the RecFast scheme only:
 *   -- ptdw->x_H, ptdw-> x_He (all neglecting reionisation, which is accounted for later on);
 *
 * - in both schemes:
 *   -- ptdw->x_noreio (neglecting reionisation);
 *   -- ptdw->x_reio (if reionisation is going on; obtained by calling
 *                    thermodynamics_reionization_function() and adding something to ptdw->x_noreio)
 *
 * @param z            Input: redshift
 * @param y            Input: vector of quantities to integrate with evolver
 * @param pth          Input: pointer to thermodynamics structure
 * @param ptw          Input/Output: pointer to thermo workspace. Contains output for x, ...
 * @param current_ap   Input: index of current approximation scheme
 * @return the error status
 */

int thermodynamics_ionization_fractions(
                                        double z,
                                        double * y,
                                        struct thermodynamics * pth,
                                        struct thermo_workspace * ptw,
                                        int current_ap
                                        ) {

  /** Summary: */

  /** Define local variables */
  struct thermo_diffeq_workspace * ptdw = ptw->ptdw;
  struct thermo_vector * ptv = ptdw->ptv;

  /* Thermo quantities */
  double x_H, x_He, xHeII, x, Tmat;
  /* Analytical quantities */
  double rhs, sqrt_val;

  /** - Calculate x_noreio from Recfast/Hyrec in each approximation
      regime. Store the results in the workspace. */
  /** - --> For credits, see external/wrap_recfast.c */

  /* Set Tmat from the y vector (it is always evolved). */
  Tmat = y[ptv->index_ti_D_Tmat] + ptw->Tcmb*(1.+z);

  /** - --> first regime: H and Helium fully ionized */
  if (current_ap == ptdw->index_ap_brec) {

    /* This is equivalent to the formula for HeIII --> HeII in Saha, just using rhs' = 1/rhs */
    rhs = ptw->SIunit_nH0/exp( 1.5*log(ptw->const_NR_numberdens*Tmat/(1.+z)/(1.+z)) - ptw->const_Tion_HeII/Tmat );
    sqrt_val = sqrt(pow(1.-rhs*(1.+ptw->fHe),2) + 4.*rhs*(1.+2*ptw->fHe));

    x = 2.*(1+2.*ptw->fHe)/(1.-rhs*(1.+ptw->fHe) + sqrt_val);

    ptdw->x_H = 1.;
    ptdw->x_He = 1.;

  }
  /** - --> second regime: first Helium recombination (analytic approximation) */
  else if (current_ap == ptdw->index_ap_He1) {

    /* Assuming Saha equilibrium for HeIII --> HeII */
    rhs = exp( 1.5*log(ptw->const_NR_numberdens*Tmat/(1.+z)/(1.+z)) - ptw->const_Tion_HeII/Tmat ) / ptw->SIunit_nH0;
    sqrt_val = sqrt(pow((rhs-1.-ptw->fHe),2) + 4.*(1.+2.*ptw->fHe)*rhs);

    x = 0.5*(sqrt_val - (rhs-1.-ptw->fHe));

    ptdw->x_H = 1.;
    ptdw->x_He = 1.;

  }
  /** - --> third regime: first Helium recombination finished, H and Helium fully ionized */
  else if (current_ap == ptdw->index_ap_He1f) {

    /* Assuming Saha equilibrium for HeII --> HeI with HII fully ionized, again expanding in rhs' = 1/rhs compared to below */
    rhs = 0.25*ptw->SIunit_nH0/exp( 1.5*log(ptw->const_NR_numberdens*Tmat/(1.+z)/(1.+z)) - ptw->const_Tion_HeI/Tmat );
    sqrt_val = sqrt(pow(1.-rhs,2) + 4.*rhs*(1.+ptw->fHe));

    x = 2.*(1+ptw->fHe)/(1.-rhs + sqrt_val);

    ptdw->x_H = 1.;
    ptdw->x_He = 1.;

  }
  /** - --> fourth regime: second Helium recombination starts (analytic approximation) */
  else if (current_ap == ptdw->index_ap_He2) {

    /* Assuming Saha equilibrium for HeII --> HeI with HII fully ionized */
    rhs = 4.*exp(1.5*log(ptw->const_NR_numberdens*Tmat/(1.+z)/(1.+z)) - ptw->const_Tion_HeI/Tmat ) / ptw->SIunit_nH0;
    sqrt_val = sqrt(pow((rhs-1.),2) + 4.*(1.+ptw->fHe)*rhs );

    x = 0.5*(sqrt_val - (rhs-1.));

    ptdw->x_H = 1.;
    ptdw->x_He = (x-1.)/ptw->fHe;

  }
  /** - --> fifth regime: Hydrogen recombination starts (analytic approximation)
      while Helium recombination continues (full equation) */
  else if (current_ap == ptdw->index_ap_H) {

    rhs = exp(1.5*log(ptw->const_NR_numberdens*Tmat/(1.+z)/(1.+z)) - ptw->const_Tion_H/Tmat)/ptw->SIunit_nH0;

    /* Assuming Saha equilibrium for HII->HI. Includes xHeII corrections from incomplete recombination of HeII --> HeI (non-zero x_HeII) */
    xHeII = y[ptv->index_ti_x_He]*ptw->fHe;
    x_H = 2./(1.+xHeII/rhs + sqrt((1.+xHeII/rhs)*(1.+xHeII/rhs)+4./rhs));

    x_He = y[ptv->index_ti_x_He];
    x = x_H + ptw->fHe * x_He;

    ptdw->x_H = x_H;
    ptdw->x_He = x_He;

  }
  /** - --> sixth regime: full Hydrogen and Helium equations */
  else if (current_ap == ptdw->index_ap_frec) {
    x_H = y[ptv->index_ti_x_H];
    x_He = y[ptv->index_ti_x_He];
    x = x_H + ptw->fHe * x_He;

    ptdw->x_H = x_H;
    ptdw->x_He = x_He;

  }
  /** - --> seventh regime: calculate x_noreio during reionization
      (i.e. x before taking reionisation into account) */
  else if (current_ap == ptdw->index_ap_reio) {

    x_H = y[ptv->index_ti_x_H];
    x_He = y[ptv->index_ti_x_He];
    x = x_H + ptw->fHe * x_He;

    ptdw->x_H = x_H;
    ptdw->x_He = x_He;

  }

  ptdw->x_noreio = x;

  /** - If z is during reionization, also calculate the reionized x */
  if (current_ap == ptdw->index_ap_reio) {

    /* set x from the evolver (which is very low ~10^-4) as 'xe_before' */
    ptw->ptrp->reionization_parameters[ptw->ptrp->index_re_xe_before] = x;

    /* compute x */
    class_call(thermodynamics_reionization_function(z,pth,ptw->ptrp,&x),
               pth->error_message,
               pth->error_message);
  }

  ptdw->x_reio = x;

  return _SUCCESS_;
}

/**
 * This subroutine contains the reionization function \f$ X_e(z) \f$ (one for each scheme) and gives x for a given z.
 *
 * @param z     Input: redshift
 * @param pth   Input: pointer to thermodynamics structure, to know which scheme is used
 * @param preio Input: pointer to reionization parameters of the function \f$ X_e(z) \f$
 * @param x     Output: \f$ X_e(z) \f$
 */

int thermodynamics_reionization_function(
                                         double z,
                                         struct thermodynamics * pth,
                                         struct thermo_reionization_parameters * preio,
                                         double * x
                                         ) {

  /** Summary: */

  /** - define local variables */
  double argument;
  int i;
  double z_jump;

  int jump;
  double center,before, after,width,one_jump;
  double z_min, z_max;

  switch (pth->reio_parametrization) {

    /** - no reionization means nothing to be added to xe_before */
  case reio_none:
    *x = preio->reionization_parameters[preio->index_re_xe_before];
    break;

    /** - implementation of ionization function similar to the one in CAMB */
  case reio_camb:

    /** - --> case z > z_reio_start */
    if (z > preio->reionization_parameters[preio->index_re_reio_start]) {
      *x = preio->reionization_parameters[preio->index_re_xe_before];
    }
    else {
      /** - --> case z < z_reio_start: hydrogen contribution (tanh of complicated argument) */
      argument = (pow((1.+preio->reionization_parameters[preio->index_re_reio_redshift]),
                    preio->reionization_parameters[preio->index_re_reio_exponent])
                   -pow((1.+z),preio->reionization_parameters[preio->index_re_reio_exponent]))
                 /(preio->reionization_parameters[preio->index_re_reio_exponent]
                 *pow((1.+preio->reionization_parameters[preio->index_re_reio_redshift]),
                    (preio->reionization_parameters[preio->index_re_reio_exponent]-1.)))
        /preio->reionization_parameters[preio->index_re_reio_width];

      *x = (preio->reionization_parameters[preio->index_re_xe_after]
            -preio->reionization_parameters[preio->index_re_xe_before])
        *(tanh(argument)+1.)/2.
        +preio->reionization_parameters[preio->index_re_xe_before];

      /** - --> case z < z_reio_start: helium contribution (tanh of simpler argument) */
      argument = (preio->reionization_parameters[preio->index_re_helium_fullreio_redshift] - z)
        /preio->reionization_parameters[preio->index_re_helium_fullreio_width];

      *x += preio->reionization_parameters[preio->index_re_helium_fullreio_fraction]
        *(tanh(argument)+1.)/2.;
    }
    break;

    /** - implementation of half-tangent like in 1209.0247 */
  case reio_half_tanh:

    /** - --> case z > z_reio_start */
    if (z > preio->reionization_parameters[preio->index_re_reio_start]) {
      *x = preio->reionization_parameters[preio->index_re_xe_before];
    }
    else {
      /** - --> case z < z_reio_start: hydrogen contribution (tanh of complicated argument) */
      argument = (pow((1.+preio->reionization_parameters[preio->index_re_reio_redshift]),
                    preio->reionization_parameters[preio->index_re_reio_exponent])
                   -pow((1.+z),preio->reionization_parameters[preio->index_re_reio_exponent]))
                 /(preio->reionization_parameters[preio->index_re_reio_exponent]
                 *pow((1.+preio->reionization_parameters[preio->index_re_reio_redshift]),
                    (preio->reionization_parameters[preio->index_re_reio_exponent]-1.)))
        /preio->reionization_parameters[preio->index_re_reio_width];

      /* argument goes from 0 to infty, not from -infty to infty like
         in reio_camb case. Thus tanh(argument) goes from 0 to 1, not
         from -1 to 1.  */

      *x = (preio->reionization_parameters[preio->index_re_xe_after]
            -preio->reionization_parameters[preio->index_re_xe_before])
        *tanh(argument)
        +preio->reionization_parameters[preio->index_re_xe_before];
    }
    break;

    /** - implementation of binned ionization function similar to astro-ph/0606552 */
  case reio_bins_tanh:

    /** - --> case z > z_reio_start */
    if (z > preio->reionization_parameters[preio->index_re_first_z+preio->re_z_size-1]) {
      *x = preio->reionization_parameters[preio->index_re_xe_before];
    }
    else if (z < preio->reionization_parameters[preio->index_re_first_z]) {
      *x = preio->reionization_parameters[preio->index_re_first_xe];
    }
    else {
      i = 0;
      while (preio->reionization_parameters[preio->index_re_first_z+i+1]<z) i++;

      /* fix the final xe to xe_before*/
      preio->reionization_parameters[preio->index_re_first_xe+preio->re_z_size-1] = preio->reionization_parameters[preio->index_re_xe_before];

      /* This is the expression of the tanh-like jumps of the reio_bins_tanh scheme until the 10.06.2015. It appeared to be
         not robust enough. It could lead to a kink in xe(z) near the maximum value of z at which reionisation is sampled. It has
         been replaced by the simpler and more robust expression below.

      *xe = preio->reionization_parameters[preio->index_re_first_xe+i]
        +0.5*(tanh((2.*(z-preio->reionization_parameters[preio->index_re_first_z+i])
                    /(preio->reionization_parameters[preio->index_re_first_z+i+1]
                      -preio->reionization_parameters[preio->index_re_first_z+i])-1.)
                   /preio->reionization_parameters[preio->index_re_step_sharpness])
              /tanh(1./preio->reionization_parameters[preio->index_re_step_sharpness])+1.)
        *(preio->reionization_parameters[preio->index_re_first_xe+i+1]
          -preio->reionization_parameters[preio->index_re_first_xe+i]);
      */

      /* compute the central redshift value of the tanh jump */
      if (i == preio->re_z_size-2) {
        z_jump = preio->reionization_parameters[preio->index_re_first_z+i]
          + 0.5*(preio->reionization_parameters[preio->index_re_first_z+i]
                 -preio->reionization_parameters[preio->index_re_first_z+i-1]);
      }
      else{
        z_jump =  0.5*(preio->reionization_parameters[preio->index_re_first_z+i+1]
                       + preio->reionization_parameters[preio->index_re_first_z+i]);
      }

      /* implementation of the tanh jump */
      *x = preio->reionization_parameters[preio->index_re_first_xe+i]
        +0.5*(tanh((z-z_jump)
                   /preio->reionization_parameters[preio->index_re_step_sharpness])+1.)
        *(preio->reionization_parameters[preio->index_re_first_xe+i+1]
          -preio->reionization_parameters[preio->index_re_first_xe+i]);
    }
    break;

    /** - implementation of many tanh jumps */
  case reio_many_tanh:

    /** - --> case z > z_reio_start */
    if (z > preio->reionization_parameters[preio->index_re_first_z+preio->re_z_size-1]) {
      *x = preio->reionization_parameters[preio->index_re_xe_before];
    }
    else if (z > preio->reionization_parameters[preio->index_re_first_z]) {

      *x = preio->reionization_parameters[preio->index_re_xe_before];

      /* fix the final xe to xe_before*/
      preio->reionization_parameters[preio->index_re_first_xe+preio->re_z_size-1] = preio->reionization_parameters[preio->index_re_xe_before];

      for (jump=1; jump<preio->re_z_size-1; jump++) {

        center = preio->reionization_parameters[preio->index_re_first_z+preio->re_z_size-1-jump];

        /* before and after are meant with respect to growing z, not growing time */
        before = preio->reionization_parameters[preio->index_re_first_xe+preio->re_z_size-1-jump]
          -preio->reionization_parameters[preio->index_re_first_xe+preio->re_z_size-jump];
        after = 0.;
        width = preio->reionization_parameters[preio->index_re_step_sharpness];

        one_jump = before + (after-before)*(tanh((z-center)/width)+1.)/2.;

        *x += one_jump;
      }

    }
    else{
      *x = preio->reionization_parameters[preio->index_re_first_xe];
    }
    break;

    /** - implementation of reio_inter */
  case reio_inter:

    /** - --> case z > z_reio_start */
    if (z > preio->reionization_parameters[preio->index_re_first_z+preio->re_z_size-1]) {
      *x = preio->reionization_parameters[preio->index_re_xe_before];
    }
    else{
      i=0;
      while (preio->reionization_parameters[preio->index_re_first_z+i+1] < z) i++;

      z_min = preio->reionization_parameters[preio->index_re_first_z+i];
      z_max = preio->reionization_parameters[preio->index_re_first_z+i+1];

      /* fix the final xe to xe_before*/
      preio->reionization_parameters[preio->index_re_first_xe+preio->re_z_size-1] = preio->reionization_parameters[preio->index_re_xe_before];

      class_test(z<z_min,
                 pth->error_message,
                 "z out of range for reionization interpolation");

      class_test(z>z_max,
                 pth->error_message,
                 "z out of range for reionization interpolation");

      argument =(z-preio->reionization_parameters[preio->index_re_first_z+i])
        /(preio->reionization_parameters[preio->index_re_first_z+i+1]
          -preio->reionization_parameters[preio->index_re_first_z+i]);

      *x = preio->reionization_parameters[preio->index_re_first_xe+i]
        + argument*(preio->reionization_parameters[preio->index_re_first_xe+i+1]
             -preio->reionization_parameters[preio->index_re_first_xe+i]);

      class_test(*x<0.,
                 pth->error_message,
                 "Interpolation gives negative ionization fraction\n",
                 argument,
                 preio->reionization_parameters[preio->index_re_first_xe+i],
                 preio->reionization_parameters[preio->index_re_first_xe+i+1]);
    }
    break;

  default:
    class_stop(pth->error_message,
               "value of reio_parametrization=%d unclear",pth->reio_parametrization);
    break;
  }
  return _SUCCESS_;
}

/**
 * Function for formatting the titles to be output
 *
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param titles     Input: titles string containing all titles
 * @return the error status
 */

int thermodynamics_output_titles(
                                 struct background * pba,
                                 struct thermodynamics *pth,
                                 char titles[_MAXTITLESTRINGLENGTH_]
                                 ) {

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
  class_store_columntitle(titles,"dTb [K]",_TRUE_);
  class_store_columntitle(titles,"w_b",_TRUE_);
  class_store_columntitle(titles,"c_b^2",_TRUE_);
  if (pba->has_idm_dr == _TRUE_) {
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
  class_store_columntitle(titles,"tau_d",_TRUE_);
  //class_store_columntitle(titles,"max. rate",_TRUE_);
  class_store_columntitle(titles,"r_d",pth->compute_damping_scale);

  return _SUCCESS_;
}

/**
 * Output the data for the output into files
 *
 * @param pba                 Input: pointer to background structure
 * @param pth                 Input: pointer to the thermodynamics structure
 * @param number_of_titles    Input: number of titles
 * @param data                Input: pointer to data file
 * @return the error status
 */

int thermodynamics_output_data(
                               struct background * pba,
                               struct thermodynamics *pth,
                               int number_of_titles,
                               double *data
                               ) {

  int index_z, storeidx;
  double *dataptr, *pvecthermo;
  double z,tau;

  // pth->number_of_thermodynamics_titles = get_number_of_titles(pth->thermodynamics_titles);
  // pth->size_thermodynamics_data = pth->number_of_thermodynamics_titles*pth->tt_size;

  /* Store quantities: */
  for (index_z=0; index_z<pth->tt_size; index_z++) {
    dataptr = data + index_z*number_of_titles;
    pvecthermo = pth->thermodynamics_table+index_z*pth->th_size;
    z = pth->z_table[index_z];
    storeidx=0;

    class_call(background_tau_of_z(pba, z, &tau),
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
    class_store_double(dataptr,pvecthermo[pth->index_th_dTb],_TRUE_,storeidx);
    class_store_double(dataptr,pvecthermo[pth->index_th_wb],_TRUE_,storeidx);
    class_store_double(dataptr,pvecthermo[pth->index_th_cb2],_TRUE_,storeidx);
    if (pba->has_idm_dr == _TRUE_) {
      class_store_double(dataptr,pvecthermo[pth->index_th_dmu_idm_dr],_TRUE_,storeidx);
      //class_store_double(dataptr,pvecthermo[pth->index_th_ddmu_idm_dr],_TRUE_,storeidx);
      //class_store_double(dataptr,pvecthermo[pth->index_th_dddmu_idm_dr],_TRUE_,storeidx);
      class_store_double(dataptr,pvecthermo[pth->index_th_tau_idm_dr],_TRUE_,storeidx);
      class_store_double(dataptr,pvecthermo[pth->index_th_tau_idr],_TRUE_,storeidx);
      class_store_double(dataptr,pvecthermo[pth->index_th_g_idm_dr],_TRUE_,storeidx);
      class_store_double(dataptr,pvecthermo[pth->index_th_cidm_dr2],_TRUE_,storeidx);
      class_store_double(dataptr,pvecthermo[pth->index_th_Tidm_dr],_TRUE_,storeidx);
      class_store_double(dataptr,pvecthermo[pth->index_th_dmu_idr],_TRUE_,storeidx);
    }
    class_store_double(dataptr,pvecthermo[pth->index_th_tau_d],_TRUE_,storeidx);
    //class_store_double(dataptr,pvecthermo[pth->index_th_rate],_TRUE_,storeidx);
    class_store_double(dataptr,pvecthermo[pth->index_th_r_d],pth->compute_damping_scale,storeidx);

  }

  return _SUCCESS_;
}
