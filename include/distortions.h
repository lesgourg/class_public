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
#include "noninjection.h"

#define _MAX_DETECTOR_NAME_LENGTH_ 100
typedef char DetectorName[_MAX_DETECTOR_NAME_LENGTH_];
typedef char DetectorFileName[_FILENAMESIZE_+_MAX_DETECTOR_NAME_LENGTH_+256];

/** List of possible branching ratio approximations */

enum br_approx {bra_sharp_sharp,bra_sharp_soft,bra_soft_soft,bra_soft_soft_cons,bra_exact};

/** List of possible schemes to compute relativistic contribution from
    reionization and structure formatio */

enum reio_approx {sd_reio_Nozawa, sd_reio_Chluba};

/**
 * distorsions structure, containing all the distortion-related parameters and
 * evolution that other modules need to know.
 */

struct distortions
{
  /** @name - input parameters initialized by user in input module
   *   (all other quantities are computed in this module, given these
   *   parameters and the content of the 'precision', 'background',
   *   'thermodynamics' and 'primordial' structures) */

  //@{

  int sd_branching_approx;                      /**< Which approximation to use for the branching ratios? */

  int sd_PCA_size;                              /**< Number of PCA components for the calculation of residual distortions */

  DetectorFileName sd_detector_file_name;       /**< Name of detector list file */

  DetectorName sd_detector_name;                /**< Name of detector */
  double sd_detector_nu_min;                    /**< Minimum frequency of chosen detector */
  double sd_detector_nu_max;                    /**< Maximum frequency of chosen detector */
  double sd_detector_nu_delta;                  /**< Bin size of chosen detector */
  int sd_detector_bin_number;                   /**< Number of frequency bins of chosen detector */
  double sd_detector_delta_Ic;                  /**< Sensitivity of the chosen detector */

  enum reio_approx sd_reio_type;                /**< Calculation method for Sunyaev Zeldovich contributions from re-ionization */

  double sd_add_y;                              /**< Possible additional y contribution (manually) to the SD signal */
  double sd_add_mu;                             /**< Possible additional mu contribution (manually) to the SD signal */

  //@}

  /** @name - Public tables and parameters */

  //@{

  /* Parameters related to redshift (z) sampling */
  double z_muy;                              /**< Redshift of the transition of mu to y era */
  double z_th;                               /**< Redshift of the transition from thermal shift to mu era */

  double z_min;                              /**< Minimum redshift */
  double z_max;                              /**< Maximum redshift */
  int z_size;                                /**< Lenght of redshift array */
  double z_delta;                            /**< Redshift intervals */
  double * z;                                /**< Redshift list z[index_z] = list of values */

  double * z_weights;                        /**< Weights for integration over z */

  /* Can be specified if no noisefile */
  double x_min;                              /**< Minimum dimentionless frequency */
  double x_max;                              /**< Maximum dimentionless frequency */
  double x_delta;                            /**< dimentionless frequency intervals */

  /* Will always be specified */
  int x_size;                                /**< Lenght of dimentionless frequency array */
  double * x;                                /**< Dimensionless frequency x[index_x] = list of values */
  double * x_weights;                        /**< Weights for integration over x */

  /* Unit conversions */
  double x_to_nu;                            /**< Conversion factor nu[GHz] = x_to_nu * x */
  double DI_units;                           /**< Conversion from unitless DI to DI[10^26 W m^-2 Hz^-1 sr^-1] */

  /* File names for the PCA */
  DetectorFileName sd_detector_noise_file;              /**< Full path of detector noise file */
  DetectorFileName sd_PCA_file_generator;               /**< Full path of PCA generator file */
  DetectorFileName sd_detector_list_file;               /**< Full path of detector list file */


  /* Tables storing branching ratios, distortions amplitudes and spectral distoritons for all types of distortios */
  double ** br_table;              /**< Branching ratios br_table[index_type][index_z] */
  double * sd_parameter_table;     /**< Spectral Distortion parameters (g,mu,y,r) sd_parameter_table[index_type] */
  double ** sd_shape_table;        /**< Spectral Distortion shapes (G,M,Y,R) sd_shape_table[index_type][index_x] */
  double ** sd_table;              /**< Spectral Distortion Intensities (final deltaI seperated by component) sd_table[index_type][index_x] */

  /* indices for the type of distortion */
  int index_type_g;                /**< temperature shift/g type distortion */
  int index_type_mu;               /**< mu type distortion */
  int index_type_y;                /**< y type distortion */
  int index_type_PCA;              /**< PCA type distortion (first index) */
  int type_size;                   /**< Number of total components for the type array */

  /* Total distortion amplitude for residual distortions */
  double epsilon;

  /* Total heating function */
  double * dQrho_dz_tot;

  /* Total heating rate */
  double Drho_over_rho;

  /* Total spectral distortion */
  double * DI;                               /**< DI[index_x] = list of values */

  /* Variables to read, allocate and interpolate external file branching_ratios_exact.dat */
  double * br_exact_z;                       /**< Redshift array for reading from file br_exact_z[index_z] */
  int br_exact_Nz;                           /**< Number of redshift values for reading from file */

  double * f_g_exact;                        /**< temperature shift/g distortion branching ratio f_g_exact[index_z] */
  double * ddf_g_exact;                      /**< second derivative of the above ddf_g_exact[index_z] */
  double * f_y_exact;                        /**< y distortion branching ratio f_y_exact[index_z] */
  double * ddf_y_exact;                      /**< second derivative of the above ddf_y_exact[index_z] */
  double * f_mu_exact;                       /**< mu distortion shape branching ratio f_mu_exact[index_z] */
  double * ddf_mu_exact;                     /**< second derivative of the above ddf_mu_exact[index_z] */

  double * E_vec;                            /**< PCA component E branching ratio for reading from file E_vec[index_e*br_exact_Nz+index_z] with index_e=[1..8] */
  double * ddE_vec;                          /**< second derivative of the above ddE_vec[index_e*br_exact_Nz+index_z] */
  int E_vec_size;                            /**< number of PCA component E branching ratios */

  /* Variable to read, allocate and interpolate external file PCA_distortions_schape.dat */
  double * PCA_nu;                           /**< Frquency array for reading from file PCA_nu[index_nu] */
  int PCA_Nnu;                               /**< Number of frequency values for reading from file */

  double * PCA_G_T;                          /**< temperature shift/g distortion shape PCA_G_T[index_nu] */
  double * ddPCA_G_T;                        /**< second derivative of the above ddPCA_G_T[index_nu] */
  double * PCA_Y_SZ;                         /**< y distortion shape PCA_Y_SZ[index_nu] */
  double * ddPCA_Y_SZ;                       /**< second derivative of the above ddPCA_Y_SZ[index_nu] */
  double * PCA_M_mu;                         /**< mu distortion shape PCA_M_mu[index_nu] */
  double * ddPCA_M_mu;                       /**< second derivative of the above ddPCA_M_mu[index_nu] */

  double * S_vec;                            /**< PCA component S shape for reading from file S_vec[index_s*S_vec_size+index_x] with index_s=[1..8] */
  double * ddS_vec;                          /**< second derivative of the above ddS_vec[index_s*S_vec_size+index_x] */
  int S_vec_size;                            /**< number of PCA component S spectral shapes */


  double * delta_Ic_array;                   /**< delta_Ic[index_x] for detectors with given sensitivity in each bin */

  //@}


  /** @name - Flags and technical parameters */

  //@{

  int has_distortions;                      /**< do we need to compute spectral distortions? */

  int has_user_defined_detector;            /**< does the user specify their own detector? */
  int has_user_defined_name;                /**< does the user specify the name of their detector? */

  int has_detector_file;                    /**< do we have a file for the detector specification? */

  int has_SZ_effect;                        /**< do we include the SZ effect? */

  int include_only_exotic;                  /**< shall we only take exotic injection contributions? */
  int include_g_distortion;                 /**< shall we include the g distortion in the total distortion ?  */

  int has_noninjected;                      /**< do we have terms that are not injected (like dissipation of acoustic waves)? */

  struct noninjection ni;                   /**< noninjection file structure */

  short distortions_verbose;                /**< flag regulating the amount of information sent to standard output (none if set to zero) */

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
                       struct thermodynamics * pth,
                       struct perturbations * ppt,
                       struct primordial * ppm,
                       struct distortions * psd);

  int distortions_constants(struct precision* ppr,
                            struct background * pba,
                            struct thermodynamics * pth,
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
                               struct thermodynamics* pth,
                               struct distortions* psd);

  /* The main computation methods */
  int distortions_compute_branching_ratios(struct precision * ppr,
                                           struct distortions* psd);

  int distortions_compute_heating_rate(struct precision* ppr,
                                       struct background* pba,
                                       struct thermodynamics * pth,
                                       struct perturbations * ppt,
                                       struct primordial * ppm,
                                       struct distortions * psd);

  int distortions_compute_spectral_shapes(struct precision * ppr,
                                          struct background * pba,
                                          struct thermodynamics * pth,
                                          struct distortions * psd);

  /* Additional sources of distortions due to recombination and LSS formation */
  int distortions_add_effects_reio(struct background * pba,
                                   struct thermodynamics * pth,
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
