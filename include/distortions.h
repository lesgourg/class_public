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
enum br_approx {bra_sharp_sharp,bra_sharp_soft,bra_soft_soft,bra_soft_soft_cons,bra_exact};


struct distortions
{
  /** @name - input parameters initialized by user in input module (all other quantities are computed in this module,
   *   given these parameters and the content of the 'precision', 'background', 'thermodynamics' and 
   *  'primordial' structures) */

  //@{

  int has_distortions;

  int branching_approx;                      /* Which approximation to use for the branching ratios? */

  int dQrho_dz_diss_approx;                  /* Use full version of dQrho_dz_diss or its approximation? */

  int N_PCA;
  char * detector;                           /* Name of detector */
  double nu_min_detector;                    /* Minimum frequency of chosen detector */
  double nu_max_detector;                    /* Maximum frequency of chosen detector */
  int nu_delta_detector;                     /* Bin size of chosen detector */


  //@}


  /** @name - Public tables and parameters */

  //@{

  char ** distortion_names;                  /* Names of the distortions */

  /* Precision parameters */
  double z_muy;
  double z_th;

  double z_min;                              /* Minimum redshift */
  double z_max;                              /* Maximum redshift */
  int z_size;                                /* Lenght of redshift array */
  double z_delta;                            /* Redshift intervals */
  double * z;                                /* z[index_z] = list of values */

  double * z_weights;

  double x_min;                              /* Minimum dimentionless frequency */
  double x_max;                              /* Maximum dimentionless frequency */
  int x_size;                                /* Lenght of dimentionless frequency array */
  double x_delta;                            /* dimentionless frequency intervals */
  double * x;                                /* x[index_x] = list of values */

  double x_to_nu;                            /* Conversion factor nu[GHz] = x_to_nu * x */
  double DI_units;                           /* Conversion from unitless DI to DI[10^26 W m^-2 Hz^-1 sr^-1] */

  /* Table storing heating rates */
  double ** heating_table;
  int index_ht_dQrho_dz_cool;                /* Heating function from cooling of electron and baryions */
  int index_ht_dQrho_dz_diss;                /* Heating function from Silk damping */
  int index_ht_dQrho_dz_CRR;                 /* Heating function from cosmological recombination radiation */
  int index_ht_dQrho_dz_ann;                 /* Heating function from particle annihilation */
  int index_ht_dQrho_dz_dec;                 /* Heating function from particle decay */
  int index_ht_dQrho_dz_eva_PBH;             /* Heating function from evaporation of primordial black holes */
  int index_ht_dQrho_dz_acc_PBH;             /* Heating function from accretion of matter into primordial black holes */
  int index_ht_dQrho_dz_tot;                 /* Total heating function */
  int index_ht_dQrho_dz_tot_screened;        /* Total heating function times blackbody visibility function */
  int ht_size;                               /* Size of the allocated space for heating quantities */

  /* Tables storing branching ratios, distortions amplitudes and spectral distoritons for all types of distortios */
  double ** br_table; 
  double * sd_parameter_table;
  double ** sd_shape_table;
  double ** sd_table;

  int index_type_g; 
  int index_type_mu;
  int index_type_y;
  int index_type_PCA;
  int type_size;

  /* TODO ?? */
  double Drho_over_rho;
  double * DI;                               /* DI[index_x] = list of values */

  /* Variables to read and allocate external file Greens_data.dat */
  int Greens_Nz;
  double * Greens_z;

  double * Greens_T_ini;
  double * ddGreens_T_ini;
  double * Greens_T_last;
  double * ddGreens_T_last;
  double * Greens_rho;
  double * ddGreens_rho;

  int Greens_Nx;
  double * Greens_x;

  double * Greens_function;
  double * ddGreens_function;

  double * Greens_blackbody;
  double * ddGreens_blackbody;

  /* Variables to read, allocate and interpolate external file branching_ratios_exact.dat */
  double * br_exact_z;
  int br_exact_Nz;

  double * f_g_exact;
  double * ddf_g_exact;
  double * f_y_exact;
  double * ddf_y_exact;
  double * f_mu_exact;
  double * ddf_mu_exact;

  double * E_vec;                /* E_vec[index_e][index_z] with index_e=1-8 */
  double * ddE_vec;
  int E_vec_size;

  /* Variable to read, allocate and interpolate external file PCA_distortions_schape.dat */
  double * PCA_nu;
  int PCA_Nnu;

  double * PCA_G_T;
  double * ddPCA_G_T;
  double * PCA_Y_SZ;
  double * ddPCA_Y_SZ;
  double * PCA_M_mu;
  double * ddPCA_M_mu;

  double * S_vec;                /* S_vec[index_s][index_x] with index_e=1-8 */
  double * ddS_vec;
  int S_vec_size;

  //@}


  /** @name - technical parameters */

  //@{

  short distortions_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message;    /**< zone for writing error messages */

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

  /* Main functions */
  int distortions_init(struct precision * ppr,
                       struct background * pba,
                       struct perturbs * ppt,
                       struct thermo * pth,
                       struct primordial * ppm,
                       struct distortions * psd);

  int distortions_free(struct distortions * psd);

  /* Indices and lists */
  int distortions_indices(struct distortions * psd);

  int distortions_get_xz_lists(struct precision * ppr,
                               struct background* pba, 
                               struct thermo* pth, 
                               struct distortions* psd);

  /* The main computation methods */
  int distortions_compute_branching_ratios(struct precision * ppr,
                                           struct distortions* psd);

  int distortions_compute_heating_rate(struct precision * ppr,
                                       struct background* pba,
                                       struct perturbs * ppt,
                                       struct thermo * pth,
                                       struct primordial * ppm,
                                       struct distortions * psd);

  int distortions_compute_spectral_amplitudes(struct distortions * psd);

  int distortions_compute_spectral_shapes(struct precision * ppr,
                                          struct background * pba,
                                          struct distortions * psd);

  /* PCA decomposition (branching ratios and spectral shapes) for unknown detector */
  int distortions_read_Greens_data(struct precision * ppr,
                                   struct distortions * psd);
  int distortions_spline_Greens_data(struct distortions* psd);
  int distortions_interpolate_Greens_data(struct distortions* psd,
                                          double z,
                                          double x,
                                          double * T_ini,
                                          double * T_last,
                                          double * rho,
                                          double * Greens_function,
                                          double * Greens_blackbody,
                                          int * last_index_z,
                                          int * last_index_x);
  int distortions_free_Greens_data(struct distortions * psd);

  /* PCA decomposition (branching ratios and spectral shapes) for PIXIE */
  int distortions_read_PIXIE_br_data(struct precision * ppr,
                                     struct distortions * psd);
  int distortions_spline_PIXIE_br_data(struct distortions* psd);
  int distortions_interpolate_PIXIE_br_data(struct distortions* psd,
                                            double z,
                                            double* f_g,
                                            double* f_y,
                                            double* f_mu,
                                            double* E,
                                            int * last_index);
  int distortions_free_PIXIE_br_data(struct distortions * psd);

  int distortions_read_PIXIE_sd_data(struct precision * ppr,
                                     struct distortions * psd);
  int distortions_spline_PIXIE_sd_data(struct distortions* psd);
  int distortions_interpolate_PIXIE_sd_data(struct distortions* psd,
                                            double nu,
                                            double * G_T,
                                            double * Y_SZ,
                                            double * M_mu,
                                            double * S,
                                            int * index);
  int distortions_free_PIXIE_sd_data(struct distortions * psd);

  /* Output */
  int heating_output_titles(struct distortions * psd, char titles[_MAXTITLESTRINGLENGTH_]);
  int heating_output_data(struct distortions * psd,
                          int number_of_titles,
                          double * data);

  int distortions_output_titles(struct distortions * psd, char titles[_MAXTITLESTRINGLENGTH_]);
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
