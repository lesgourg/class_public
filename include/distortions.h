/** @file distortions.h Documented module on spectral distortions
 * Matteo Lucca, 31.10.2018
 * Nils Schoeneberg, 18.02.2019
 */

#ifndef __DISTORTIONS__
#define __DISTORTIONS__

#include "arrays.h"
#include "background.h"
#include "thermodynamics.h"
#include "perturbations.h"
#include "primordial.h"

#define _MAX_DETECTOR_NAME_LENGTH_ 100
typedef char DetectorName[_MAX_DETECTOR_NAME_LENGTH_];
typedef char DetectorFileName[_FILENAMESIZE_+_MAX_DETECTOR_NAME_LENGTH_+256];

/**
 * All deistortions parameters and evolution that other modules need to know.
 */
enum br_approx {bra_sharp_sharp,bra_sharp_soft,bra_soft_soft,bra_soft_soft_cons,bra_exact};
enum reio_approx {sd_reio_Nozawa, sd_reio_Chluba};


struct distortions
{
  /** @name - input parameters initialized by user in input module (all other quantities are computed in this module,
   *   given these parameters and the content of the 'precision', 'background', 'thermodynamics' and
   *  'primordial' structures) */

  //@{

  int sd_branching_approx;                      /* Which approximation to use for the branching ratios? */

  int sd_PCA_size;

  FileName sd_detector_file_name;               /* Name of detector list file */

  DetectorName sd_detector_name;                /* Name of detector */
  double sd_detector_nu_min;                    /* Minimum frequency of chosen detector */
  double sd_detector_nu_max;                    /* Maximum frequency of chosen detector */
  double sd_detector_nu_delta;                  /* Bin size of chosen detector */
  int sd_detector_bin_number;
  double sd_detector_delta_Ic;

  int sd_reio_type;

  //@}


  /** @name - Public tables and parameters */

  //@{

  /* Precision parameters */
  double z_muy;
  double z_th;

  double z_min;                              /* Minimum redshift */
  double z_max;                              /* Maximum redshift */
  int z_size;                                /* Lenght of redshift array */
  double z_delta;                            /* Redshift intervals */
  double * z;                                /* z[index_z] = list of values */

  double * z_weights;

  /* Can be specified if no noisefile */
  double x_min;                              /* Minimum dimentionless frequency */
  double x_max;                              /* Maximum dimentionless frequency */
  double x_delta;                            /* dimentionless frequency intervals */
  /* Will always be specified */
  int x_size;                                /* Lenght of dimentionless frequency array */
  double * x;                                /* x[index_x] = list of values */

  double * x_weights;

  /* Unit conversions */
  double x_to_nu;                            /* Conversion factor nu[GHz] = x_to_nu * x */
  double DI_units;                           /* Conversion from unitless DI to DI[10^26 W m^-2 Hz^-1 sr^-1] */

  /* File names for the PCA */
  DetectorFileName sd_detector_noise_file;              /* Full path of detector noise file */
  DetectorFileName sd_PCA_file_generator;               /* Full path of PCA generator file */
  DetectorFileName sd_detector_list_file;               /* Full path of detector list file */


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

  /* Total distortion amplitude for residual distortions */
  double epsilon;

  /* Total heating function */
  double * dQrho_dz_tot;

  /* Total heating rate */
  double Drho_over_rho;

  /* Total spectral distortion */
  double * DI;                               /* DI[index_x] = list of values */

  /* Variables to read, allocate and interpolate external file branching_ratios_exact.dat */
  double * br_exact_z;
  int br_exact_Nz;

  double * f_g_exact;
  double * ddf_g_exact;
  double * f_y_exact;
  double * ddf_y_exact;
  double * f_mu_exact;
  double * ddf_mu_exact;

  double * E_vec;                            /* E_vec[index_e][index_z] with index_e=1-8 */
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

  double * S_vec;                            /* S_vec[index_s][index_x] with index_e=1-8 */
  double * ddS_vec;
  int S_vec_size;


  double * delta_Ic_array;                   /* delta_Ic[index_x] for detectors with given sensitivity in each bin */

  //@}


  /** @name - Flags and technical parameters */

  //@{

  int has_distortions;

  int user_defined_detector;
  int user_defined_name;

  int has_detector_file;

  int has_SZ_effect;

  int only_exotic; /**< flag specifying wether to only take exotic injection contributions */

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
                       struct thermo * pth,
                       struct perturbs * ppt,
                       struct primordial * ppm,
                       struct distortions * psd);

  int distortions_constants(struct precision* ppr,
                            struct background * pba,
                            struct thermo * pth,
                            struct distortions * psd);

  int distortions_free(struct distortions * psd);

  /* PCA decomposition (branching ratios and spectral shapes) for unknown detector */
  int distortions_generate_detector(struct precision * ppr,
                                    struct distortions * psd);

  int distortions_set_detector(struct precision * ppr,
                               struct distortions* psd);

  int distortions_read_detector_noisefile(struct precision * ppr,
                                          struct distortions * psd);

  /* Indices and lists */
  int distortions_indices(struct distortions * psd);

  int distortions_get_xz_lists(struct precision * ppr,
                               struct background* pba,
                               struct thermo* pth,
                               struct distortions* psd);

  /* The main computation methods */
  int distortions_compute_branching_ratios(struct precision * ppr,
                                           struct distortions* psd);

  int distortions_compute_heating_rate(struct background* pba,
                                       struct thermo * pth,
                                       struct perturbs * ppt,
                                       struct primordial * ppm,
                                       struct distortions * psd);

  int distortions_compute_spectral_shapes(struct precision * ppr,
                                          struct background * pba,
                                          struct thermo * pth,
                                          struct distortions * psd);

  /* Additional sources of distortions due to recombination and LSS formation */
  int distortions_add_effects_reio(struct background * pba,
                                   struct thermo * pth,
                                   struct distortions * psd,
                                   double T_e,
                                   double Dtau,
                                   double beta,
                                   double beta_z,
                                   double x,
                                   double * y_reio,
                                   double * DI);

  /* PCA decomposition (branching ratios and spectral shapes) for known detector */
  int distortions_read_br_data(struct precision * ppr,
                               struct distortions * psd);
  int distortions_spline_br_data(struct distortions* psd);
  int distortions_interpolate_br_data(struct distortions* psd,
                                      double z,
                                      double* f_g,
                                      double* f_y,
                                      double* f_mu,
                                      double* E,
                                      int * last_index);
  int distortions_free_br_data(struct distortions * psd);

  int distortions_read_sd_data(struct precision * ppr,
                               struct distortions * psd);
  int distortions_spline_sd_data(struct distortions* psd);
  int distortions_interpolate_sd_data(struct distortions* psd,
                                      double nu,
                                      double * G_T,
                                      double * Y_SZ,
                                      double * M_mu,
                                      double * S,
                                      int * index);
  int distortions_free_sd_data(struct distortions * psd);

  /* Output */
  int distortions_output_heat_titles(struct distortions * psd, char titles[_MAXTITLESTRINGLENGTH_]);
  int distortions_output_heat_data(struct distortions * psd,
                                   int number_of_titles,
                                   double * data);

  int distortions_output_sd_titles(struct distortions * psd, char titles[_MAXTITLESTRINGLENGTH_]);
  int distortions_output_sd_data(struct distortions * psd,
                                 int number_of_titles,
                                 double * data);

#ifdef __cplusplus
}
#endif

/**************************************************************/


#endif
/* @endcond */
