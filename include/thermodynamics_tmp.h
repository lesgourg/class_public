/** @file thermodynamics.h Documented includes for thermodynamics module */

#ifndef __THERMODYNAMICS__
#define __THERMODYNAMICS__

#include "background.h"
//#include "arrays.h"
//#include "helium.h"
//#include "hydrogen.h"

/**
 * List of possible recombination algorithms.
 */

enum recombination_algorithm {
  recfast,
  hyrec,
  cosmorec
};

/**
 * List of possible reionization schemes.
 */

enum reionization_parametrization {
  reio_none, /**< no reionization */
  reio_camb,  /**< reionization parameterized like in CAMB */
  reio_bins_tanh,  /**< binned reionization history with tanh inteprolation between bins */
  reio_half_tanh,  /**< half a tanh, intead of the full tanh */
  reio_many_tanh,  /**< similar to reio_camb but with more than one tanh */
  reio_stars_sfr_source_term, /**< Reionization parameterization based on the star formation rate, see Poulin et al. arXiv:1508.01370 and references therein */
  reio_duspis_et_al,/**< Redshift asymetric reionisation parametrization as introduced by Duspis et al. 1509.02785 and improved by 1605.03928 */
  reio_asymmetric_planck_16 /**< Redshift asymetric reionisation parametrization as introduced by the Planck collaboration in 2016 data release 1605.03507  */
};
/**
 * List of possible energy repartition functions.
 */
enum energy_repartition_functions {
  SSCK, /**< Shull Van Stanbeerg Chen Kamionkowski parameterization */
  Galli_et_al_fit,  /**< Fit of Galli et al 2013 functions */
  Galli_et_al_interpolation,  /**< Interpolation of Galli et al 2013 functions  */
  no_factorization
};
/**
 * List of possible modelisation of star reheating.
 */
enum heating_by_stars_parametrization {
  heating_none, /**< No reheating by stars */
  heating_reiolike_tanh, /**< reheating parameterized like reionization with normalization adjust to fit data */
  heating_stars_sfr_source_term  /**< reheating parameterization based on the star formation rate (SFR). See Poulin et al. 1508.01370.  */
};
enum modelisation_of_SFR {
  model_SFR_bestfit, /**< No reheating by stars */
  model_SFR_p1sig, /**< reheating parameterized like reionization with normalization adjust to fit data */
  model_SFR_m1sig,  /**< reheating parameterization based on the star formation rate (SFR). See Poulin et al. 1508.01370.  */
  model_SFR_free
};
enum PBH_accretion_recipe {
  Ali_Haimoud, /**< Accretion recipe from Ali_Haimoud & Kamionkowski, arXiv:1612.05644 */
  Ricotti_et_al,  /**< Accretion recipe from Ricotti et al., arXiv:0709.0524 */
  Gaggero_et_al,  /**< Accretion recipe from Gaggero et al., arXiv:1612.00457 */
  Horowitz, /**< Accretion recipe from Horowitz, arXiv:1612.07264 */
  Hybrid /**<A more realistic accretion recipe, with a transition from spherical to disk accretion at a redshift "PBH_disk_formation_redshift" */
};
enum energy_deposition_treatment {
  Analytical_approximation, /**< Analytical energy deposition treatment, introduced in 1209.0247 and corrected in 1612.05644 */
  Slatyer  /**< f(z) functions from Slatyer, introduced in 1211.0283 and updated in 1506.03812 */
};

/**
 * Is the input parameter the reionization redshift or optical depth?
 */

enum reionization_z_or_tau {
  reio_z,  /**< input = redshift */
  reio_tau /**< input = tau */
};

/**
 * Two useful smooth step functions, for smoothing transitions in recfast.
 */

#define f1(x) (-0.75*x*(x*x/3.-1.)+0.5)  /**< goes from 0 to 1 when x goes from -1 to 1 */
#define f2(x) (x*x*(0.5-x/3.)*6.)        /**< goes from 0 to 1 when x goes from  0 to 1 */

/**
 * All thermodynamics parameters and evolution that other modules need to know.
 *
 * Once initialized by thermodynamics_init(), contains all the
 * necessary information on the thermodynamics, and in particular, a
 * table of thermodynamical quantities as a function of the redshift,
 * used for interpolation in other modules.
 */

struct thermo
{
  /** @name - input parameters initialized by user in input module
   *  (all other quantities are computed in this module, given these parameters
   *   and the content of the 'precision' and 'background' structures) */

  //@{

  double YHe;  /**< \f$ Y_{He} \f$: primordial helium fraction */
  double Lambda_over_theoritical_Lambda; /**< ratio of A2s1s transition with respect to theoritical value 8.2206 */

  enum recombination_algorithm recombination; /**< recombination code */

  enum reionization_parametrization reio_parametrization; /**< reionization scheme */

  enum reionization_z_or_tau reio_z_or_tau; /**< is the input parameter the reionization redshift or optical depth? */

  enum heating_by_stars_parametrization star_heating_parametrization; /**< star heating parametrization */

  enum modelisation_of_SFR model_SFR; /**< choice of SFR modelling. Currently Roberston et al 15: bestfit, m1sig, p1sig. */

  double tau_reio; /**< if above set to tau, input value of reionization optical depth */

  double z_reio;   /**< if above set to z,   input value of reionization redshift */

  double z_opt_depth;   /**< if above set to a positive value,  compute the optical depth at the required redshift. if set to 0, compute the optical depth at the minimum of the visibility function (default). */

  short compute_cb2_derivatives; /**< do we want to include in computation derivatives of baryon sound speed? */

  short compute_damping_scale; /**< do we want to compute the simplest analytic approximation to the photon damping (or diffusion) scale? */

  /** parameters for reio_camb */

  double reionization_width; /**< width of H reionization */

  double reionization_exponent; /**< shape of H reionization */

  double helium_fullreio_redshift; /**< redshift for of helium reionization */

  double helium_fullreio_width; /**< width of helium reionization */

  /** parameters for reio_bins_tanh */

  int binned_reio_num; /**< with how many bins do we want to describe reionization? */

  double * binned_reio_z; /**< central z value for each bin */

  double * binned_reio_xe; /**< imposed \f$ X_e(z)\f$ value at center of each bin */

  double binned_reio_step_sharpness; /**< sharpness of tanh() step interpolating between binned values */

    /** parameters for reio_many_tanh */

  int many_tanh_num; /**< with how many jumps do we want to describe reionization? */

  double * many_tanh_z; /**< central z value for each tanh jump */

  double * many_tanh_xe; /**< imposed \f$ X_e(z)\f$ value at the end of each jump (ie at later times)*/

  double many_tanh_width; /**< sharpness of tanh() steps */

  /** parameters used by duspis et al. parametrization */

  double Qp_duspis_et_al;
  double zp_duspis_et_al;
  double lambda_duspis_et_al;

  /** parameters used by planck 16 asymmetric parametrization */

  double z_end_asymmetric_planck_16;
  double z_start_asymmetric_planck_16;
  double alpha_asymmetric_planck_16;

  /** parameters used by reio_stars_sfr_source_term and heating_by_stars_parametrization*/


  double f_esc; /**< fraction of photons produced by stellar populations that escape to ionize the IGM */
  double Zeta_ion; /**< Lyman continuum photon production efficiency of the stellar population */
  double Log10_Zeta_ion;/**< The log10 of former parameter. */
  double fx; /**< X-ray efficiency fudge factor of photons responsible for heating the medium. */
  double Ex; /**< Associated normalization from Pober et al. 1503.00045. */
  double ap;   /**<  a few parameters entering the fit of the star formation rate (SFR), introduced in Madau & Dickinson, Ann.Rev.Astron.Astrophys. 52 (2014) 415-486, updated in Robertson & al. 1502.02024.*/
  double bp;
  double cp;
  double dp;
  double z_start_reio_stars; /**< Controls the beginning of star reionisation, the SFR experiences is put to 0 above this value. */

  /** parameter used by the tanh reheating */
  double final_IGM_temperature; /**< Controls the final temperature of the IGM if a tanh reheating is required. Other parameters (duration, starting point...) are the same as the reionisation tanh. */

  /** some derived parameters useful to compare reionization models */

  double duration_of_reionization;   /**< it measures the duration of reionation, it is defined as z_10_percent - z_99_percent*/
  double z_10_percent;  /** <redshift at which x_e = 0.1*(1+f_He) */
  double z_50_percent;  /** <redshift at which x_e = 0.5*(1+f_He) */
  double z_99_percent; /** <redshift at which x_e = 0.99*(1+f_He) */

  /** parameters for energy injection */

  double annihilation; /** parameter describing CDM annihilation (f <sigma*v> / m_cdm, see e.g. 0905.0003) */
  double annihilation_boost_factor;
  double annihilation_m_DM;
  double increase_T_from_stars;
  short has_on_the_spot; /** flag to specify if we want to use the on-the-spot approximation **/
  short reio_stars_and_dark_matter;  /* switch that indicates if DM decay or halos are switched on to better combine star reionisation and DM */
  enum energy_repartition_functions energy_repart_functions; /**< energy repartition functions */
  double decay; /** parameter describing CDM decay (f/tau, see e.g. 1109.6322)*/
  double PBH_mass; /**< mass from the PBH, in case of Dark Matter being high masses PBH */
  enum PBH_accretion_recipe PBH_accretion_recipe; /**< recipe to compute accretion from PBH */
  double PBH_disk_formation_redshift; /**< Disk formation redshift, in case of Dark Matter being high masses PBH and realistic accretion model*/
  enum energy_deposition_treatment energy_deposition_treatment; /**< Treatment of energy deposition in the medium following DM annihilation, decay, PBH evaporation etc. */
  double PBH_low_mass; /**< mass from the PBH, in case of Dark Matter being low mass PBH */
  double PBH_fraction; /**< fraction of Dark Matter being PBH */

  double annihilation_variation; /** if this parameter is non-zero,
				     the function F(z)=(f <sigma*v> /
				     m_cdm)(z) will be a parabola in
				     log-log scale between zmin and
				     zmax, with a curvature given by
				     annihlation_variation (must be
				     negative), and with a maximum in
				     zmax; it will be constant outside
				     this range */

  double annihilation_z; /** if annihilation_variation is non-zero,
			     this is the value of z at which the
			     parameter annihilation is defined, i.e.
			     F(annihilation_z)=annihilation */

  double annihilation_zmax; /** if annihilation_variation is non-zero,
				redshift above which annihilation rate
				is maximal */

  double annihilation_zmin; /** if annihilation_variation is non-zero,
				redshift below which annihilation rate
				is constant */

  double annihilation_f_halo; /** takes the contribution of DM annihilation in halos into account*/
  double annihilation_z_halo; /** characteristic redshift for DM annihilation in halos*/

  double * annihil_coef_xe;
  double * annihil_coef_heat;
  double * annihil_coef_lya;
  double * annihil_coef_ionH;
  double * annihil_coef_ionHe;
  double * annihil_coef_lowE;
  double * annihil_coef_dd_heat;
  double * annihil_coef_dd_lya;
  double * annihil_coef_dd_ionH;
  double * annihil_coef_dd_ionHe;
  double * annihil_coef_dd_lowE;

  double chi_heat;
  double chi_lya;
  double chi_ionH;
  double chi_ionHe;
  double chi_lowE;
  int annihil_coef_num_lines;

  //@}

  /** @name - all indices for the vector of thermodynamical (=th) quantities stored in table */

  //@{

  int index_th_xe;            /**< ionization fraction \f$ x_e \f$ */
  int index_th_dkappa;        /**< Thomson scattering rate \f$ d \kappa / d \tau\f$ (units 1/Mpc) */
  int index_th_tau_d;         /**< Baryon drag optical depth */
  int index_th_ddkappa;       /**< scattering rate derivative \f$ d^2 \kappa / d \tau^2 \f$ */
  int index_th_dddkappa;      /**< scattering rate second derivative \f$ d^3 \kappa / d \tau^3 \f$ */
  int index_th_exp_m_kappa;  /**< \f$ exp^{-\kappa} \f$ */
  int index_th_g;             /**< visibility function \f$ g = (d \kappa / d \tau) * exp^{-\kappa} \f$ */
  int index_th_dg;            /**< visibility function derivative \f$ (d g / d \tau) \f$ */
  int index_th_ddg;           /**< visibility function second derivative \f$ (d^2 g / d \tau^2) \f$ */
  int index_th_Tb;            /**< baryon temperature \f$ T_b \f$ */
  int index_th_cb2;           /**< squared baryon sound speed \f$ c_b^2 \f$ */
  int index_th_dcb2;          /**< derivative wrt conformal time of squared baryon sound speed \f$ d [c_b^2] / d \tau \f$ (only computed if some non-minimal tight-coupling schemes is requested) */
  int index_th_ddcb2;         /**< second derivative wrt conformal time of squared baryon sound speed  \f$ d^2 [c_b^2] / d \tau^2 \f$ (only computed if some non0-minimal tight-coupling schemes is requested) */
  int index_th_rate;          /**< maximum variation rate of \f$ exp^{-\kappa}\f$, g and \f$ (d g / d \tau) \f$, used for computing integration step in perturbation module */
  int index_th_r_d;           /**< simple analytic approximation to the photon comoving damping scale */
  int th_size;                /**< size of thermodynamics vector */

  //@}

  /** @name - thermodynamics interpolation tables */

  //@{

  int tt_size; /**< number of lines (redshift steps) in the tables */
  double * z_table; /**< vector z_table[index_z] with values of redshift (vector of size tt_size) */
  double * thermodynamics_table; /**< table thermodynamics_table[index_z*pth->tt_size+pba->index_th] with all other quantities (array of size th_size*tt_size) */

  //@}

  /** @name - table of their second derivatives, used for spline interpolation */

  //@{

  double * d2thermodynamics_dz2_table; /**< table d2thermodynamics_dz2_table[index_z*pth->tt_size+pba->index_th] with values of \f$ d^2 t_i / dz^2 \f$ (array of size th_size*tt_size) */

  //@}


  /** @name - redshift, conformal time and sound horizon at recombination */

  //@{

  double z_rec;   /**< z at which the visibility reaches its maximum (= recombination redshift) */
  double tau_rec; /**< conformal time at which the visibility reaches its maximum (= recombination time) */
  double rs_rec;  /**< comoving sound horizon at recombination */
  double ds_rec;  /**< physical sound horizon at recombination */
  double ra_rec;  /**< conformal angular diameter distance to recombination */
  double da_rec;  /**< physical angular diameter distance to recombination */
  double rd_rec;  /**< comoving photon damping scale at recombination */
  double z_d;     /**< baryon drag redshift */
  double tau_d;   /**< baryon drag time */
  double ds_d;    /**< physical sound horizon at baryon drag */
  double rs_d;    /**< comoving sound horizon at baryon drag */
  double tau_cut; /**< at at which the visibility goes below a fixed fraction of the maximum visibility, used for an approximation in perturbation module */
  double angular_rescaling; /**< [ratio ra_rec / (tau0-tau_rec)]: gives CMB rescaling in angular space relative to flat model (=1 for curvature K=0) */

  //@}

  /** @name - redshift, conformal time and sound horizon at recombination */

  //@{

  double tau_free_streaming;   /**< minimum value of tau at which sfree-streaming approximation can be switched on */

  //@}

  /** @name - initial conformal time at which thermodynamical variables have been be integrated */

  //@{

  double tau_ini; /**< initial conformal time at which thermodynamical variables have been be integrated */

  //@}

/** @name - total number density of electrons today (free or not) */

  //@{

  double n_e; /**< total number density of electrons today (free or not) */

  //@}

  /**
   *@name - some flags needed for thermodynamics functions
   */

  //@{

  short inter_normal;  /**< flag for calling thermodynamics_at_z and find position in interpolation table normally */
  short inter_closeby; /**< flag for calling thermodynamics_at_z and find position in interpolation table starting from previous position in previous call */

  //@}

  /** @name - technical parameters */

  //@{

  short thermodynamics_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */

  //@}

};


/**
 * Temporary structure where all the recombination history is defined and stored.
 *
 * This structure is used internally by the thermodynamics module,
 * but never passed to other modules.
 */

struct recombination {

  /** @name - indices of vector of thermodynamics variables related to recombination */

  //@{

  int index_re_z;          /**< redshift \f$ z \f$ */
  int index_re_xe;         /**< ionization fraction \f$ x_e \f$ */
  int index_re_Tb;         /**< baryon temperature \f$ T_b \f$ */
  int index_re_cb2;        /**< squared baryon sound speed \f$ c_b^2 \f$ */
  int index_re_dkappadtau; /**< Thomson scattering rate \f$ d \kappa / d \tau \f$ (units 1/Mpc) */
  int re_size;             /**< size of this vector */

  //@}

  /** @name - table of the above variables at each redshift, and number of redshifts */

  //@{

  int rt_size; /**< number of lines (redshift steps) in the table */
  double * recombination_table; /**< table recombination_table[index_z*preco->re_size+index_re] with all other quantities (array of size preco->rt_size*preco->re_size) */

  //@}

  /** @name - recfast parameters needing to be passed to
      thermodynamics_derivs_with_recfast() routine */

  //@{


  /** parameters for energy injection */

  double CDB; /**< defined as in RECFAST */
  double CR;  /**< defined as in RECFAST */
  double CK;  /**< defined as in RECFAST */
  double CL;  /**< defined as in RECFAST */
  double CT;  /**< defined as in RECFAST */
  double fHe; /**< defined as in RECFAST */
  double CDB_He; /**< defined as in RECFAST */
  double CK_He;  /**< defined as in RECFAST */
  double CL_He;  /**< defined as in RECFAST */
  double fu; /**< defined as in RECFAST */
  double H_frac; /**< defined as in RECFAST */
  double Tnow;   /**< defined as in RECFAST */
  double Nnow;   /**< defined as in RECFAST */
  double Bfact;  /**< defined as in RECFAST */
  double CB1;    /**< defined as in RECFAST */
  double CB1_He1; /**< defined as in RECFAST */
  double CB1_He2; /**< defined as in RECFAST */
  double H0;  /**< defined as in RECFAST */
  double YHe; /**< defined as in RECFAST */

  /* parameters for energy injection */

  double annihilation; /**< parameter describing CDM annihilation (f <sigma*v> / m_cdm, see e.g. 0905.0003) */
  double annihilation_boost_factor;/**< alternative parameterization to annihilation parameter, describes the boost factor to annihilation cross section */
  double annihilation_m_DM; /**< in case of alternative parameterization to annihilation parameter, describes the mass of the dark matter */

  short has_on_the_spot; /**< flag to specify if we want to use the on-the-spot approximation **/

  double decay; /**< parameter describing CDM decay (f/tau, see e.g. 1109.6322)*/
  double PBH_mass; /**< mass from the PBH, in case of Dark Matter being PBH */
  enum PBH_accretion_recipe PBH_accretion_recipe; /**< recipe to compute accretion from PBH */
  double PBH_disk_formation_redshift; /**< Disk formation redshift, in case of Dark Matter being high masses PBH and realistic accretion model*/
  enum energy_deposition_treatment energy_deposition_treatment; /**< Treatment of energy deposition in the medium following DM annihilation, decay, PBH evaporation etc. */
  double PBH_low_mass; /**< mass from the PBH, in case of Dark Matter being low mass PBH */
  double PBH_fraction; /**< fraction of Dark Matter being PBH */
  double Tm_tmp; /**< To temporarily store the value of the matter temperature*/
  double xe_tmp; /**< To temporarily store the value of the free electron fraction */
  double PBH_low_mass_tmp; /**< To temporarily store the value of the (low) mass pbh*/
  double exponent_pbh_integral; /**< To temporarily store the value of the exponent in the integral of the density of low mass pbh*/
  double z_tmp; /**< To temporarily store the value of the redshift*/
  double annihilation_variation; /**< if this parameter is non-zero,
				     the function F(z)=(f <sigma*v> /
				     m_cdm)(z) will be a parabola in
				     log-log scale between zmin and
				     zmax, with a curvature given by
				     annihlation_variation (must be
				     negative), and with a maximum in
				     zmax; it will be constant outside
				     this range */

  double annihilation_z; /**< if annihilation_variation is non-zero,
			     this is the value of z at which the
			     parameter annihilation is defined, i.e.
			     F(annihilation_z)=annihilation */

  double annihilation_zmax; /**< if annihilation_variation is non-zero,
				redshift above which annihilation rate
				is maximal */

  double annihilation_zmin; /**< if annihilation_variation is non-zero,
				redshift below which annihilation rate
				is constant */


  double annihilation_f_halo; /**< takes the contribution of DM annihilation in halos into account*/
  double annihilation_z_halo; /**< characteristic redshift for DM annihilation in halos*/

  //@}
  /** A few parameters useful if realistic energy deposition is required in case of annihilations in halos or energy injection due to decay of short lived DM */
  double * annihil_z;
  double * annihil_f_eff;
  double * annihil_dd_f_eff;

  double f_eff;
  int annihil_f_eff_num_lines;
  double increase_T_from_stars; /**< To compute the increase of temperature due to star formation in a given model */

  enum energy_repartition_functions energy_repart_functions; /**< energy repartition functions */


  ErrorMsg error_message;
};

/**
 * Temporary structure where all the reionization history is defined and stored.
 *
 * This structure is used internally by the thermodynamics module,
 * but never passed to other modules.
 */

struct reionization {

  /** @name - indices of vector of thermodynamics variables related to reionization */

  //@{

  int index_re_z;          /**< redshift \f$ z \f$ */
  int index_re_xe;         /**< ionization fraction \f$ x_e \f$ */
  int index_re_Tb;         /**< baryon temperature \f$ T_b \f$ */
  int index_re_cb2;        /**< squared baryon sound speed \f$ c_b^2 \f$ */
  int index_re_dkappadtau; /**< Thomson scattering rate \f$ d \kappa / d \tau\f$ (units 1/Mpc) */
  int index_re_dkappadz;   /**< Thomson scattering rate with respect to redshift \f$ d \kappa / d z\f$ (units 1/Mpc) */
  int index_re_d3kappadz3; /**< second derivative of previous quantity with respect to redshift */
  int re_size;             /**< size of this vector */

  //@}

  /** @name - table of the above variables at each redshift, and number of redshifts */

  //@{

  int rt_size;                 /**< number of lines (redshift steps) in the table */
  double * reionization_table; /**< table reionization_table[index_z*preio->re_size+index_re] with all other quantities (array of size preio->rt_size*preio->re_size) */

  //@}

  /** @name - reionization optical depth inferred from reionization history */

  //@{

  double reionization_optical_depth; /**< reionization optical depth inferred from reionization history */

  //@}

  /** @name - indices describing input parameters used in the definition of the various possible functions x_e(z) */

  //@{

  /* parameters used by reio_camb */

  int index_reio_redshift;  /**< hydrogen reionization redshift */
  int index_reio_exponent;  /**< an exponent used in the function x_e(z) in the reio_camb scheme */
  int index_reio_width;     /**< a width defining the duration of hydrogen reionization in the reio_camb scheme */
  int index_reio_xe_before; /**< ionization fraction at redshift 'reio_start' */
  int index_reio_xe_after;  /**< ionization fraction after full reionization */
  int index_helium_fullreio_fraction; /**< helium full reionization fraction inferred from primordial helium fraction */
  int index_helium_fullreio_redshift; /**< helium full reionization redshift */
  int index_helium_fullreio_width;    /**< a width defining the duration of helium full reionization in the reio_camb scheme */

  /* parameters used by reio_bins_tanh and reio_many_tanh */

  int reio_num_z; /**< number of reionization jumps */
  int index_reio_first_z; /**< redshift at which we start to impose reionization function */
  int index_reio_first_xe; /**< ionization fraction at redshift first_z (inferred from recombination code) */
  int index_reio_step_sharpness; /**< sharpness of tanh jump */

  /* parameters used by duspis et al. parametrization */

  int index_Qp_duspis_et_al;
  int index_zp_duspis_et_al;
  int index_lambda_duspis_et_al;

  /* parameters used by planck 16 asymmetric parametrization */

  int index_z_end_asymmetric_planck_16;
  int index_z_start_asymmetric_planck_16;
  int index_alpha_asymmetric_planck_16;


  /* parameters used by all schemes */

  int index_reio_start;     /**< redshift above which hydrogen reionization neglected */

  //@}

  /** @name - vector of such parameters, and its size */

  double * reionization_parameters; /**< vector containing all reionization parameters necessary to compute xe(z) */
  int reio_num_params; /**< length of vector reionization_parameters */

  //@}

  /** @name - index of line in recombination table corresponding to first line of reionization table */

  //@{

  int index_reco_when_reio_start; /**< index of line in recombination table corresponding to first line of reionization table*/

  //@}

};

/**
 * temporary  parameters and workspace passed to the thermodynamics_derivs function
 */

struct thermodynamics_parameters_and_workspace {

  /* structures containing fixed input parameters (indices, ...) */
  struct background * pba;
  struct precision * ppr;
  struct recombination * preco;
  struct thermo * pth;
  /* workspace */
  double * pvecback;

};

/**************************************************************/
/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int thermodynamics_at_z(
			  struct background * pba,
			  struct thermo * pth,
			  double z,
			  short inter_mode,
			  int * last_index,
			  double * pvecback,
			  double * pvecthermo
			  );

  int thermodynamics_init(
			  struct precision * ppr,
			  struct background * pba,
			  struct thermo * pth
			  );

  int thermodynamics_free(
			  struct thermo * pthermo
			  );

  int thermodynamics_indices(
			     struct thermo * pthermo,
			     struct recombination * preco,
			     struct reionization * preio
			     );

  int thermodynamics_helium_from_bbn(
				     struct precision * ppr,
				     struct background * pba,
				     struct thermo * pth
				     );
  int thermodynamics_annihilation_coefficients_init(
                                                    struct precision * ppr,
                                                    struct background * pba,
                                                    struct thermo * pth
                                                  );
  int thermodynamics_annihilation_coefficients_interpolate(
                                                     struct precision * ppr,
                                                     struct background * pba,
                                                     struct thermo * pth,
                                                     double xe,
                                                     double * chi_heat,
                                                     double * chi_lya,
                                                     double * chi_ionH,
                                                     double * chi_ionHe,
                                                     double * chi_lowE
                                                   );
  int thermodynamics_annihilation_coefficients_free(
                                                   struct thermo * pth
                                                 );
  int thermodynamics_annihilation_f_eff_init(
                                                   struct precision * ppr,
                                                   struct background * pba,
                                                   struct recombination * preco
                                                 );
  int thermodynamics_annihilation_f_eff_interpolate(
                                                    struct precision * ppr,
                                                    struct background * pba,
                                                    struct recombination * preco,
                                                    double z,
                                                    double * f_eff
                                                  );
  int thermodynamics_annihilation_f_eff_free(
                                                  struct recombination * preco
                                                );
  int thermodynamics_onthespot_energy_injection(
				      struct precision * ppr,
				      struct background * pba,
				      struct recombination * preco,
				      double z,
				      double * energy_rate,
				      ErrorMsg error_message
				      );
  int thermodynamics_beyond_onthespot_energy_injection(
                                                struct precision * ppr,
                                                struct background * pba,
                                                struct recombination * preco,
                                                double z,
                                                double *energy_rate,
                                                ErrorMsg error_message
                                              );
  int thermodynamics_energy_injection(
				      struct precision * ppr,
				      struct background * pba,
				      struct recombination * preco,
 				      double z,
				      double * energy_rate,
				      ErrorMsg error_message
				      );

  int thermodynamics_reionization_function(
					   double z,
					   struct thermo * pth,
					   struct reionization * preio,
             struct recombination * preco,
					   double * xe
					   );

  int thermodynamics_reionization(
				  struct precision * ppr,
				  struct background * pba,
				  struct thermo * pth,
				  struct recombination * preco,
				  struct reionization * preio,
				  double * pvecback
				  );

  int thermodynamics_reionization_sample(
					 struct precision * ppr,
					 struct background * pba,
					 struct thermo * pth,
					 struct recombination * preco,
					 struct reionization * preio,
					 double * pvecback
					 );

  int thermodynamics_get_xe_before_reionization(
                                                struct precision * ppr,
                                                struct thermo * pth,
                                                struct recombination * preco,
                                                double z,
                                                double * xe);

  int thermodynamics_recombination(
				   struct precision * ppr,
				   struct background * pba,
				   struct thermo * pth,
				   struct recombination * prec,
				   double * pvecback
				   );

  int thermodynamics_recombination_with_hyrec(
						struct precision * ppr,
						struct background * pba,
						struct thermo * pth,
						struct recombination * prec,
						double * pvecback
						);
            
  int thermodynamics_recombination_with_cosmorec(
						struct precision * ppr,
						struct background * pba,
						struct thermo * pth,
						struct recombination * prec,
						double * pvecback
						);

  int thermodynamics_recombination_with_recfast(
						struct precision * ppr,
						struct background * pba,
						struct thermo * pth,
						struct recombination * prec,
						double * pvecback
						);

  int thermodynamics_derivs_with_recfast(
					 double z,
					 double * y,
					 double * dy,
					 void * fixed_parameters,
					 ErrorMsg error_message
					 );

  int thermodynamics_merge_reco_and_reio(
					 struct precision * ppr,
					 struct thermo * pth,
					 struct recombination * preco,
					 struct reionization * preio
					 );

  int thermodynamics_output_titles(struct background * pba,
                                   struct thermo *pth,
                                   char titles[_MAXTITLESTRINGLENGTH_]
                                   );

  int thermodynamics_output_data(struct background * pba,
                                 struct thermo *pth,
                                 int number_of_titles,
                                 double *data
                                 );
  int get_optical_depth_at_z_or_min_g(struct thermo *pth,
                               double *z_opt_depth,
                               double *opt_depth
                             );
  int thermodynamics_tanh(double x,
                          double center,
                          double before,
                          double after,
                          double width,
                          double * result);

#ifdef __cplusplus
}
#endif

/**************************************************************/

/**
 * @name some flags
 */

//@{

#define _BBN_ -1

//@}

/**
 * @name Some basic constants needed by RECFAST:
 */

//@{

#define _m_e_ 9.10938215e-31  /**< electron mass in Kg */
#define _m_p_ 1.672621637e-27 /**< proton mass in Kg */
#define _m_H_ 1.673575e-27    /**< Hydrogen mass in Kg */
#define _not4_ 3.9715         /**< Helium to Hydrogen mass ratio */
#define _sigma_ 6.6524616e-29 /**< Thomson cross-section in m^2 */

//@}

/**
 * @name Some specific constants needed by RECFAST:
 */

//@{

#define _RECFAST_INTEG_SIZE_ 3

#define _Lambda_ 8.2206 /*Updated value from (Labzowsky et al 2005)*/
// #define _Lambda_ 8.2245809 /*Old value from recfast original */
#define _Lambda_He_ 51.3
#define _L_H_ion_ 1.096787737e7
#define _L_H_alpha_ 8.225916453e6
#define _L_He1_ion_ 1.98310772e7
#define _L_He2_ion_ 4.389088863e7
#define _L_He_2s_ 1.66277434e7
#define _L_He_2p_ 1.71134891e7
#define	_A2P_s_		1.798287e9     /*updated like in recfast 1.4*/
#define	_A2P_t_		177.58e0       /*updated like in recfast 1.4*/
#define	_L_He_2Pt_	1.690871466e7  /*updated like in recfast 1.4*/
#define	_L_He_2St_	1.5985597526e7 /*updated like in recfast 1.4*/
#define	_L_He2St_ion_	3.8454693845e6 /*updated like in recfast 1.4*/
#define	_sigma_He_2Ps_	1.436289e-22   /*updated like in recfast 1.4*/
#define	_sigma_He_2Pt_	1.484872e-22   /*updated like in recfast 1.4*/

//@}

/**
 * @name Some specific constants needed by recfast_derivs:
 */

//@{

#define _a_PPB_ 4.309
#define _b_PPB_ -0.6166
#define _c_PPB_ 0.6703
#define _d_PPB_ 0.5300
#define _T_0_ pow(10.,0.477121)   /* from recfast 1.4 */
#define _a_VF_ pow(10.,-16.744)
#define _b_VF_ 0.711
#define _T_1_ pow(10.,5.114)
#define	_a_trip_ pow(10.,-16.306) /* from recfast 1.4 */
#define	_b_trip_ 0.761            /* from recfast 1.4 */

//@}

/**
 * @name Some limits imposed on cosmological parameter values:
 */
/* @endcond */
//@{

#define _YHE_BIG_ 0.5      /**< maximal \f$ Y_{He} \f$ */
#define _YHE_SMALL_ 0.01   /**< minimal \f$ Y_{He} \f$ */
#define _Z_REC_MAX_ 2000.
#define _Z_REC_MIN_ 500.

//@}

#endif
