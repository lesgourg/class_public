#ifndef __NONINJECTION__
#define __NONINJECTION__

#include "common.h" //Use here ONLY the things required for defining the struct (i.e. common.h for the ErrorMsg)

struct noninjection{

  /** @name - Imported parameters */

  //@{
  /* Arrays related to wavenumbers */
  double k_min;
  double k_max;
  int k_size;
  double* k;
  double* k_weights;
  double* pk_primordial_k;

  /* Array related to WKB approximation for diss. of acc. waves */
  double* integrand_approx;

  /* Arrays related to redshift */
  double* z_table_coarse;
  int z_size_coarse;
  double logz_max;
  double * noninjection_table;
  double * ddnoninjection_table;

  int z_size;
  double* z_table;
  int last_index_z;

  /* Temporary quantities */
  // WKB approximation quantities
  double f_nu_wkb;
  double dkD_dz;
  double kD;
  // Fixed thermodynamic quantities
  double H0;
  double T_g0;
  double Omega0_b;
  double Omega0_cdm;
  double rho0_cdm;
  // Fixed thermodynamic quantities
  double fHe;
  double N_e0;
  // Varying background quantities
  double H;
  double a;
  double rho_g;
  // Varying thermodynamic quantities
  double heat_capacity;
  double nH;
  double T_b;
  double T_g;
  double x_e;

  //@}

  /** @name - Public tables and parameters */

  //@{

  double* photon_dep_table;  /* The table of energy depositions into the photon fluid */

  //@}

  /** @name - Flags and technical parameters */

  //@{

  /* Error message */
  ErrorMsg error_message;

  //@}

};

/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int noninjection_init(struct precision* ppr,
                        struct background* pba,
                        struct thermodynamics* pth,
                        struct perturbations* ppt,
                        struct primordial* ppm,
                        struct noninjection* pni);

  int noninjection_free(struct noninjection* pni);

  int noninjection_photon_heating_at_z(struct noninjection* pni,
                                       double z,
                                       double* heat);

  int noninjection_rate_adiabatic_cooling(struct noninjection * pni,
                                          double z,
                                          double * energy_rate);

  int noninjection_rate_acoustic_diss(struct noninjection * pni,
                                      double z,
                                      double * energy_rate);

  int noninjection_output_titles(struct noninjection * pni,
                                 char titles[_MAXTITLESTRINGLENGTH_]);

  int noninjection_output_data(struct noninjection * pni,
                               int number_of_titles,
                               double * data);

#ifdef __cplusplus
}
#endif

#endif
