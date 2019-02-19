/** @file distortions.h Documented module on spectral distortions
 * Matteo Lucca, 31.10.2018
 */

#ifndef __DISTORTIONS__
#define __DISTORTIONS__

#include "arrays.h"
#include "background.h"
#include "thermodynamics.h"
#include "perturbations.h"
#include "primordial.h"


/**
 * All deistortions parameters and evolution that other modules need to know.
 */


struct distortions
{
  /** @name - input parameters initialized by user in input module
   *  (all other quantities are computed in this module, given these
   *  parameters and the content of the 'precision', 'background',
   *  'thermodynamics' and 'primordial' structures) */

  //@{

  int has_distortions;

  int branching_approx;                      /* Which approximation to use for the branching ratios? */
  int dQrho_dz_diss_approx;                  /* Use full version of dQrho_dz_diss or its approximation? */

  //@}


  /** @name - all indices for the vector of heating functions and distortions (=sd) stored in table */

  //@{

  /* Precision parameters */
  double z_min;                              /* Minimum redshift */
  double z_max;                              /* Maximum redshift */
  double z_delta;                            /* Redshift intervals */
  int z_size;                                /* Lenght of redshift array */

  double x_min;                              /* Minimum dimentionless frequency */
  double x_max;                              /* Maximum dimentionless frequency */
  double x_delta;                            /* dimentionless frequency intervals */
  int x_size;                                /* Lenght of dimentionless frequency array */

  /* z-dependent parameters */
  int index_ht_dQrho_dz_cool;                /* Heating function from cooling of electron and baryions */
  int index_ht_dQrho_dz_diss;                /* Heating function from Silk damping */
  int index_ht_dQrho_dz_ann;                 /* Heating function from particle annihilation */
  int index_ht_dQrho_dz_dec;                 /* Heating function from particle decay */
  int index_ht_dQrho_dz_eva_PBH;             /* Heating function from evaporation of primordial black holes */
  int index_ht_dQrho_dz_acc_PBH;             /* Heating function from accretion of matter into primordial black holes */
  int index_ht_dQrho_dz_tot;                 /* Total heating function */
  int ht_size;                               /* Size of the allocated space for heating quantities */

  double* bb_visibility_function;            /* bb = blackbody [NS] */
  double ** branching_ratios;                /* [index_br][index_z] [NS] */
  double * sd_parameter;
  int index_br_f_g;                          /* Branching ratios */
  int index_br_f_mu;
  int index_br_f_y;
  int index_br_f_r;
  int br_size;

  /* x-dependent parameters */
  int index_sd_Y;                            /* Shape of y distortions */
  int index_sd_M;                            /* Shape of mu distortions */
  int index_sd_G;                            /* Shape of shifted power specrum */
  int index_sd_DI;                           /* Shape of final distortions */
  int sd_size;                               /* Size of the allocated space for distortions quantities */

  /* Output parameters */
  double * z;                                /* z[index_z] = list of values */
  double * dQrho_dz_tot;                     /* dQrho_dz_tot[index_z] = list of values */

  double g;                                  /* g-parameter */
  double mu;                                 /* mu-parameter */
  double y;                                  /* y-parameter */
  double r;                                  /* r-parameter */
  double Drho_over_rho;                      /* Total emitted/injected heat */

  double * x;                                /* x[index_x] = list of values */
  double * DI;                               /* DI[index_x] = list of values */

  //@}


  /** @name - Internal parameter used only in distortions.c */

  //@{

  double * z_weights;

  /* Variable from external file Greens_data.dat */
  int Greens_Nz;
  double * Greens_z;
  int Greens_Nx;
  double * Greens_x;
  double * Greens_T_ini;
  double * Greens_T_last;
  double * Greens_rho;
  double * Greens_blackbody;
  double * Greens_function;

  /* Variable from external file branching_ratios_exact.dat */
  int br_exact_Nz;
  int br_exact_N_columns;
  double * br_exact_z;
  double * f_g_exact;
  double * f_y_exact;
  double * f_mu_exact;
  int index_e;
  double * E_vec;                /* E_vec[index_e][index_z] with index_e=1-8 */
  double * br_exact_table;         

  /* Variable from external file PCA_distortions_schape.dat */
  int PCA_Nx;
  int PCA_N_columns;
  double * PCA_x;
  double * PCA_J_T;
  double * PCA_J_y;
  double * PCA_J_mu;
  int index_s;
  double * S_vec;                /* S_vec[index_s][index_x] with index_e=1-8 */
  double * PCA_table;         


  //@}


  /** @name - technical parameters */

  //@{

  short distortions_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */

  //@}

};

/*************************************************************************************************************/
/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int distortions_init(struct precision * ppr,
                       struct background * pba,
                       struct perturbs * ppt,
                       struct thermo * pth,
                       struct primordial * ppm,
                       struct distortions * psd);

  int distortions_free(struct distortions * psd);

  int distortions_indices(struct distortions * psd);

  int distortions_get_xz_lists(struct distortions* psd);

  int distortions_branching_ratios(struct precision * ppr,
                                   struct distortions* psd);

  int distortions_visibility_function(struct background* pba,
                                      struct thermo * pth,
                                      struct distortions* psd);

  int heating_at_z(struct precision * ppr,
                   struct background* pba,
                   struct perturbs * ppt,
                   struct thermo * pth,
                   struct primordial * ppm,
                   struct distortions * psd,
                   double z,
                   double * pvecheat);

  int distortions_at_x(struct background* pba,
                       struct distortions * psd,
                       double x,
                       double * pvecdist);

  int distortions_read_Greens_data(struct precision * ppr,
                                   struct distortions * psd);
  int distortions_free_Greens_data(struct distortions * psd);

  int distortions_read_BR_exact_data(struct precision * ppr,
                                     struct distortions * psd);
  int distortions_free_BR_exact_data(struct distortions * psd);

  int distortions_read_PCA_dist_shapes_data(struct precision * ppr,
                                            struct distortions * psd);
  int distortions_free_PCA_dist_shapes_data(struct distortions * psd);

  int heating_output_titles(char titles[_MAXTITLESTRINGLENGTH_]);

  int heating_output_data(struct distortions * psd,
                          int number_of_titles,
                          double * data);

  int distortions_output_titles(char titles[_MAXTITLESTRINGLENGTH_]);

  int distortions_output_data(struct distortions * psd,
                              int number_of_titles,
                              double * data);

#ifdef __cplusplus
}
#endif

/**************************************************************/

/**
 * @name Some conversion factors needed by distortions module:
 */

//@{

#define _s_over_Mpc_ 9.7157e-15  /**< conversion factor from s to megaparsecs (1 s= const*Mpc) */
#define _Mpc_over_GeV_ 1.5637e38  /**< conversion factor from GeV to megaparsecs (1 GeV= const/Mpc) */
#define _GeV_over_kg_ 1.7827e-27  /**< conversion factor from GeV to kg  (1 GeV= const*kg) */
#define _GeVcm3_over_Mpc4_ 0.01056  /**< conversion factor from GeV/cm^3 to 1/Mpc^4 (GeV/cm^3=const/Mpc^4) */

//@}


#endif
/* @endcond */
