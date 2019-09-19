#ifndef __HEATING__
#define __HEATING__

#include "common.h" //Use here ONLY the things required for defining the struct (i.e. common.h for the ErrorMsg)

/**
 * All heating parameters and evolution that other modules need to know.
 */
enum PBH_accretion_approx {spherical_accretion, disk_accretion};
enum f_eff_approx {f_eff_on_the_spot, f_eff_from_file};
enum chi_approx {chi_CK, chi_PF, chi_Galli_file, chi_Galli_analytic, chi_full_heating, chi_from_x_file, chi_from_z_file};

struct heating{

  /** @name - input parameters initialized by user in input module (all other quantities are computed in this module,
   *   given these parameters and the content of the 'precision', 'background', 'thermodynamics' and
   *  'primordial' structures) */

  //@{

  /* Exotic energy injection parameters */
  double DM_annihilation_efficiency;
  double DM_annihilation_cross_section;
  double DM_annihilation_mass;
  double DM_annihilation_fraction;
  double DM_annihilation_variation;
  double DM_annihilation_z;
  double DM_annihilation_zmax;
  double DM_annihilation_zmin;
  double DM_annihilation_f_halo;
  double DM_annihilation_z_halo;

  double DM_decay_fraction;
  double DM_decay_Gamma;

  double PBH_evaporation_fraction;
  double PBH_evaporation_mass;

  double PBH_accretion_fraction;
  double PBH_accretion_mass;
  double PBH_accretion_recipe;
  double PBH_accretion_relative_velocities;
  double PBH_accretion_ADAF_delta;
  double PBH_accretion_eigenvalue;

  /* Injection efficiency */
  int f_eff_type;
  FileName f_eff_file;

  /* Deposition function and injection efficiency */
  int chi_type;
  FileName chi_z_file;
  FileName chi_x_file;

  /* Approximation for energy injection of acoustic waves dissipation */
  int heating_rate_acoustic_diss_approx;

  //@}


  /** @name - Imported parameters */

  //@{

  /* Parameters from background structure */
  int Nz_size;
  double z_initial;
  double z_start_chi_approx;

  /* Parameters from background structure */
  /* Redshift independent, i.e. defined in heating_init */
  double H0;
  double T_g0;
  double Omega0_b;
  double Omega0_cdm;
  double rho0_cdm;
  double f_nu_wkb;
  /* Redshift dependent, i.e. defined in heating_calculate_at_z or heating_at_z_second_order */
  double H;
  double a;
  double t;
  double R;
  double rho_g;
  double rho_b;
  double rho_cdm;
  double rho_dcdm;

  /* Parameters from thermodynamics structure */
  /* Redshift independent, i.e. defined in heating_init */
  double Y_He;
  double fHe;
  double heat_capacity;
  double N_e0;
  double nH;
  /* Redshift dependent, i.e. defined in heating_calculate_at_z or heating_at_z_second_order */
  double T_b;
  double T_g;
  double x_e;
  double dkappa;
  double dkD_dz;
  double kD;

  /* Parameters from primordial structure */
  double k_max;
  double k_min;
  double k_size;
  double* k;
  double* pk_primordial_k;

  //@}

  /** @name - Public tables and parameters */

  //@{

  /* Redshift tables */
  double* z_table;
  int z_size;

  double tol_z_table;
  int filled_until_index_z;
  double filled_until_z;

  int last_index_z_feff;
  int last_index_z_chi;
  int last_index_z_inj;
  int last_index_z;

  int index_z_store;

  /* X_e table */
  int last_index_x_chi;

  /* PBH mass evolution table and PBH free parameters */
  double PBH_z_evaporation;
  double PBH_QCD_activation;
  double * PBH_table_z;
  double * PBH_table_mass;
  double * PBH_table_mass_dd;
  double * PBH_table_F;
  double * PBH_table_F_dd;
  int Nz_PBH;

  /* Energy injection table */
  double** injection_table;
  int index_inj_cool;
  int index_inj_diss;
  int index_inj_DM_ann;
  int index_inj_DM_dec;
  int index_inj_PBH_eva;
  int index_inj_PBH_acc;
  int index_inj_tot;
  int inj_size;                  /** All contributions + total */

  /* Injection efficiency table */
  double f_eff;
  int feff_z_size;
  double* feff_table;

  /* Deposition function tables */
  double* chiz_table;
  int chiz_size;
  double* chix_table;
  int chix_size;

  double** deposition_table; /* The table of energy depositions into the IGM of different deposition types */
  double* photon_dep_table;  /* The table of energy depositions into the photon fluid */
  double* chi;
  int index_dep_heat;
  int index_dep_ionH;
  int index_dep_ionHe;
  int index_dep_lya;
  int index_dep_lowE;
  int dep_size;

  /* Energy deposition vector */
  double* pvecdeposition;

   //@}

  /** @name - Flags and technical parameters */

  //@{

  /* Flags */
  int has_exotic_injection;

  int has_DM_ann;
  int has_DM_dec;
  int has_PBH_eva;
  int has_PBH_acc;

  int to_store;

  /* Book-keeping */
  int heating_verbose;

  ErrorMsg error_message;

};



/**************************************************************/

/* *
 * Putting this down here is important, because of the special nature of this module.
 * This allows the struct heating to already be defined and thus be a normal member
 * (as opposed to a pointer member) of the struct thermo in thermodynamics.h
 * */
struct background;
struct thermo;
struct perturbs;
struct primordial;

/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  /* Allocate, define indeces for and free heating tables */
  int heating_init(struct precision * ppr,
                   struct background* pba,
                   struct thermo* pth);

  int heating_indices(struct thermo* pth);

  int heating_free(struct thermo* pth);

  /* Main functions */
  int heating_calculate_at_z(struct background* pba,
                             struct thermo* pth,
                             double x,
                             double z,
                             double Tmat,
                             double* pvecback);

  int heating_energy_injection_at_z(struct heating* phe,
                                    double z,
                                    double* dEdz_inj);

  int heating_deposition_function_at_z(struct heating* phe,
                                       double x,
                                       double z);

  int heating_add_noninjected(struct background* pba,
                              struct thermo* pth,
                              struct perturbs* ppt,
                              struct primordial* ppm);

  int heating_photon_at_z(struct thermo* pth,
                          double z,
                          double* heat);

  int heating_baryon_at_z(struct thermo* pth,
                          double z);

  /* Heating functions */
  int heating_rate_adiabatic_cooling(struct heating * phe,
                                     double z,
                                     double * energy_rate);

  int heating_rate_acoustic_diss(struct heating * phe,
                                 double z,
                                 double * energy_rate);

  int heating_rate_DM_annihilation(struct heating * phe,
                                   double z,
                                   double * energy_rate);

  int heating_rate_DM_decay(struct heating * phe,
                            double z,
                            double * energy_rate);

  int heating_rate_PBH_evaporation_mass_evolution(struct background * pba,
                                                  struct heating * phe);

  int heating_rate_PBH_evaporation(struct heating * phe,
                                   double z,
                                   double * energy_rate);

  int heating_rate_PBH_accretion(struct heating * phe,
                                 double z,
                                 double * energy_rate);

  /* Injection efficiency */
  int heating_read_feff_from_file(struct precision* ppr,
                                   struct heating* phe,
                                   char* f_eff_file);

  /* Deposition function */
  int heating_read_chi_z_from_file(struct precision* ppr,
                                   struct heating* phe,
                                   char* chi_z_file);

  int heating_read_chi_x_from_file(struct precision* ppr,
                                   struct heating* phe,
                                   char* chi_x_file);

  int heating_output_titles(struct heating* phe,char* titles_heat);

  int heating_output_data(struct heating * phe,
                          int number_of_titles,
                          double * data);

#ifdef __cplusplus
}
#endif

#endif
