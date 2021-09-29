#ifndef __INJECTION__
#define __INJECTION__

#include "common.h" //Use here ONLY the things required for defining the struct (i.e. common.h for the ErrorMsg)

/**
 * All injection parameters and evolution that other modules need to know.
 */
enum PBH_accretion_approx {spherical_accretion, disk_accretion};
enum f_eff_approx {f_eff_on_the_spot, f_eff_from_file};
enum chi_approx {chi_CK, chi_PF, chi_Galli_file, chi_Galli_analytic, chi_full_heating, chi_from_x_file, chi_from_z_file};

struct injection{

  /** @name - input parameters initialized by user in input module (all other quantities are computed in this module,
   *   given these parameters and the content of the 'precision', 'background' and 'thermodynamics' structs) */

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
  enum PBH_accretion_approx PBH_accretion_recipe;
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

  //@}


  /** @name - Imported parameters */

  //@{

  /* Parameters from precision structure */
  int Nz_size;
  double z_initial;
  double z_start_chi_approx;

  /* Parameters from background structure */
  /* Redshift independent, i.e. defined in injection_init */
  double H0;
  double T_g0;
  double Omega0_b;
  double Omega0_cdm;
  double rho0_cdm;
  /* Redshift dependent, i.e. defined in injection_calculate_at_z */
  double H;
  double a;
  double t;
  double rho_g;
  double rho_b;
  double rho_cdm;

  /* Parameters from thermodynamics structure */
  /* Redshift independent, i.e. defined in injection_init */
  double fHe;
  double N_e0;
  /* Redshift dependent, i.e. defined in injection_calculate_at_z */
  double heat_capacity;
  double T_b;
  double x_e;
  double nH;

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

  ErrorMsg error_message;

  //@}
};




/**************************************************************/

/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  /* Allocate, define indeces for and free injection tables */
  int injection_init(struct precision * ppr,
                     struct background* pba,
                     struct thermodynamics* pth);

  int injection_indices(struct thermodynamics* pth);

  int injection_free(struct thermodynamics* pth);

  /* Main functions */
  int injection_calculate_at_z(struct background* pba,
                               struct thermodynamics* pth,
                               double x,
                               double z,
                               double Tmat,
                               double* pvecback);

  int injection_energy_injection_at_z(struct injection* phe,
                                      double z,
                                      double* dEdt_inj);

  int injection_deposition_function_at_z(struct injection* phe,
                                         double x,
                                         double z);

  int injection_deposition_at_z(struct thermodynamics* pth,
                                double z);

  /* injection functions */
  int injection_rate_DM_annihilation(struct injection * phe,
                                     double z,
                                     double * energy_rate);

  int injection_rate_DM_decay(struct injection * phe,
                              double z,
                              double * energy_rate);

  int injection_rate_PBH_evaporation_mass_evolution(struct background * pba,
                                                    struct injection * phe);

  int injection_rate_PBH_evaporation(struct injection * phe,
                                     double z,
                                     double * energy_rate);

  int injection_rate_PBH_accretion(struct injection * phe,
                                   double z,
                                   double * energy_rate);

  /* Injection efficiency */
  int injection_read_feff_from_file(struct precision* ppr,
                                     struct injection* phe,
                                     char* f_eff_file);

  /* Deposition function */
  int injection_read_chi_z_from_file(struct precision* ppr,
                                     struct injection* phe,
                                     char* chi_z_file);

  int injection_read_chi_x_from_file(struct precision* ppr,
                                     struct injection* phe,
                                     char* chi_x_file);

  int injection_output_titles(struct injection* phe,
                              char titles[_MAXTITLESTRINGLENGTH_]);

  int injection_output_data(struct injection * phe,
                            int number_of_titles,
                            double * data);

#ifdef __cplusplus
}
#endif

#endif
