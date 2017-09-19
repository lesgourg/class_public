/* Structure with all energy injection parameters */
/* If adding a new energy injection process
   make sure to add relevant parameters here */

/* The following four lines constitue a little patch for an easier class compatibility when dealing with exotic electromagnetic injection */
#include "../include/arrays.h"
#include "../include/common.h"
#define _SUCCESS_ 0 /**< integer returned after successful call of a function */
#define _FAILURE_ 1 /**< integer returned after non-successful call of a function */

typedef struct {

  double odmh2;                 /* Omega_dm h^2 */

  double pann, pann_halo;       /* DM annihilation parameter in the smooth background and in haloes */
                                /* Units of pann and pann_halo are cm^3/s/GeV */

  double ann_z, ann_zmax, ann_zmin, ann_var; /* Parameters for the variation of pann(z) */
  double ann_z_halo;                         /* Characteristic redshift for annihilation in haloes */
  double ann_f_halo;                         /* takes the contribution of DM annihilation in halos into account*/

  double Mpbh, fpbh;           /* Mass and fraction of DM made of primordial black holes */
  int coll_ion;                /* If 1: assume gas gest collisionally ionized.
                                  If 0: assume gas gets photoionized by PBH radiation */

  int on_the_spot;            /* if set to 1 assume energy deposition rate = injection rate */
                              /* Otherwise solves for deposition given injection with simple recipe */


   double zstart, zend, dlna;   /* initial and final redshift and step size in log a */
   long nz;                     /* total number of redshift steps */

   /** parameters for energy injection */

   double annihilation; /** parameter describing CDM annihilation (f <sigma*v> / m_cdm, see e.g. 0905.0003) */
   double annihilation_boost_factor;
   double annihilation_m_DM;
   short has_on_the_spot; /** do we want to use the on-the-spot approximation? */

   double decay_fraction; /**  fraction of decaying DM */
   double Gamma_dcdm; /** Inverse lifetime of the decaying DM */

   double PBH_ADAF_delta;
   double PBH_accretion_eigenvalue;
   int PBH_accretion_recipe;

   int annihil_coef_num_lines;
   double *annihil_coef_heat;
   double *annihil_coef_ionH;
   double *annihil_coef_ionHe;
   double *annihil_coef_lya;
   double *annihil_coef_lowE;
   double *annihil_coef_xe;
   double *annihil_coef_dd_heat;
   double *annihil_coef_dd_ionH;
   double *annihil_coef_dd_ionHe;
   double *annihil_coef_dd_lya;
   double *annihil_coef_dd_lowE;
   double chi_heat;
   double chi_lya;
   double chi_ionH;
   double chi_ionHe;
   double chi_lowE;


   double annihilation_z_halo; /*characteristic redshift for DM annihilation in halos*/
   double * annihil_z;
   double * annihil_f_eff;
   double * annihil_dd_f_eff;
   double f_eff;
   short  energy_repart_functions; /**< energy repartition functions */
   int energy_deposition_treatment;
   int annihil_f_eff_num_lines;

   /** for PBH evaporation */
   double PBH_low_mass; /**< mass from the evaporating PBH in g */

   short PBH_table_is_initialized; /**< Flag to specify if the PBH-mass evolution was calculated */
   double PBH_z_evaporation; /**< Double to store the evaporation redshift. Useful to avoid bad extrapolation at low z. */
   int PBH_table_size; /**< Length of the PBH-mass evolution table */
   double * PBH_table_z; /**< Array of redshift for the evolution of the PBH-mass (used for evaporation) */
   double * PBH_table_mass; /**< Array of the PBH-mass given the redshift in 'PBH_table_z' */
   double * PBH_table_mass_dd; /**< Array of the second derivative of PBH-mass w.r.t. the redshift */
   double * PBH_table_F; /**< Array of F(z)  given the redshift in 'PBH_table_z' */
   double * PBH_table_F_dd; /**< Array of the second derivative of F(z) w.r.t. the redshift */






   int reio_parametrization; /*Do we want the reio by stars based on SFR modeling ? 0 = no, 1 = yes*/
   int star_heating_parametrization; /*Do we want heating by stars based on SFR modeling ? 0 = no, 1 = yes*/
   /** A few parameters if the scheme reio_stars_sfr_source_term is chosen */

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


   double Omega0_r;
   double Omega0_b;
   double Omega0_cdm;
  //  double Omega0_dcdm;
  //  double Omega0_lambda;
   double H0;


} INJ_PARAMS;

double beta_pbh(double Mpbh, double z, double xe, double Teff);
double gamma_pbh(double Mpbh, double z, double xe, double Teff);
double lambda_pbh(double Mpbh, double z, double xe, double Teff);
double Mdot_pbh(double Mpbh, double z, double xe, double Teff);
double TS_over_me_pbh(double Mpbh, double z, double xe, double Teff, int coll_ion);
double eps_over_mdot_pbh(double Mpbh, double z, double xe, double Teff, int coll_ion);
double L_pbh(double Mpbh, double z, double xe, double Teff, int coll_ion);
double vbc_rms_func(double z);
double L_pbh_av(double Mpbh, double z, double xe, double Tgas, int coll_ion);
double L_pbh_ADAF(double z, double xe, double Tgas, INJ_PARAMS *params);
double dEdtdV_accreting_PBH(double z, double xe, double Tgas, INJ_PARAMS *params);
double dEdtdV_inj(double z, double xe, double Tgas, INJ_PARAMS *params);
void update_dEdtdV_dep(double z_out, double dlna, double xe, double Tgas,
		       double nH, double H, INJ_PARAMS *params, double *dEdtdV_dep);
int evaluate_chi_heat(INJ_PARAMS *param,double z, double xe);
int evaluate_chi_ionisation(INJ_PARAMS *param,double z, double xe);
int hyrec_annihilation_coefficients_interpolate(INJ_PARAMS *inj_params, double xe_or_z);
double dEdVdt_evaporating_PBH(double z, INJ_PARAMS *params);
double dEdtdV_DM_decay(double z, INJ_PARAMS *params);
double vbc_av(double z, double xe, double T);
