#ifndef __HEATING__
#define __HEATING__

#include "common.h" //Use here ONLY the things required for defining the struct (i.e. common.h)

struct heating{

  /* Flags */
  int deposit_energy_as;
  int has_exotic_injection;
  int has_dcdm;
  
  int has_DM_ann;
  int has_DM_dec;
 
  int to_store;

  /* Redshift tables */
  double* z_table;
  int z_size;

  double tol_z_table;
  int filled_until_index_z_dep;
  double filled_until_z_dep;
  int filled_until_index_z_inj;
  double filled_until_z_inj;
  
  int last_index_z_dep;
  int last_index_z_inj;
  int last_index_z_feff;
  int last_index_z_chi;

  /* TODO */
  int last_index_chix;

  /* chi tables */
  double* chiz_table;
  int chiz_size;
  double* chix_table;
  int chix_size;

  /* f_eff table */
  int feff_z_size;
  double* feff_table;
    
  /* Parameters from background structure */
  int to_store;

  int has_DM_ann;
  int has_DM_dec;
  
  double* injection_table;
  int index_inj_DM_ann;
  int index_inj_DM_dec;
  int index_inj_BAO;
  int index_inj_tot;
  //int index_dep_lowE;
  int inj_size;       //All contributions + total

  /* Deposition table */
  double* chi_table;
  int chi_type;
  double* deposition_table;
  int index_dep_heat;
  int index_dep_ionH;
  int index_dep_ionHe;
  int index_dep_lya;
  int index_dep_lowE;
  int dep_size;

  /* Parameters from background structure */
  double H0;
  double rho_crit0;
  double nH0;
  double Omega0_b;
  double Omega0_cdm;
  double Omega0_dcdmdr
  double rho_cdm;
  double rho_dcdm;
  double t;
  double Gamma_dcdm;
  double T_b;
  double x_e;

  int last_index_bg;

  /* Heating parameters */
  short has_on_the_spot;         /**< flag to specify if we want to use the on-the-spot approximation **/
  double f_eff;

  double annihilation_efficiency;/**< parameter describing CDM annihilation (f <sigma*v> / m_cdm, see e.g. 0905.0003) */
  double annihilation_variation; /**< if this parameter is non-zero, the function F(z)=(f <sigma*v>/m_cdm)(z) will be a parabola in
                                      log-log scale between zmin and zmax, with a curvature given by annihlation_variation (must be
                                      negative), and with a maximum in zmax; it will be constant outside this range */
  double annihilation_z;         /**< if annihilation_variation is non-zero, this is the value of z at which the parameter annihilation is defined, i.e.
                                      F(annihilation_z)=annihilation */
  double annihilation_zmax;      /**< if annihilation_variation is non-zero, redshift above which annihilation rate is maximal */
  double annihilation_zmin;      /**< if annihilation_variation is non-zero, redshift below which annihilation rate is constant */
  double annihilation_f_halo;    /**< takes the contribution of DM annihilation in halos into account*/
  double annihilation_z_halo;    /**< characteristic redshift for DM annihilation in halos*/

  double decay;                  /**< parameter describing CDM decay (f/tau, see e.g. 1109.6322)*/
  double decay_fraction;
  
  /* Heat injection table */
  double* injection_table;
  int index_inj_DM_ann;
  int index_inj_DM_dec;
  int index_inj_BAO;
  int index_inj_tot;
  //int index_dep_lowE;
  int inj_size;                  /** All contributions + total */

  /* Energy deposition table */
  double* chi_table;
  int chi_type;

  double* pvecdeposition;

  double* deposition_table;
  int index_dep_heat;
  int index_dep_ionH;
  int index_dep_ionHe;
  int index_dep_lya;
  //int index_dep_lowE;
  int dep_size;

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
#include "perturbations.h"

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
  int heating_at_z(struct background* pba,
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

  /* Branching ratios into the different channels */
  int heating_read_chi_z_from_file(struct precision* ppr,
                                   struct heating* phe);

  int heating_read_chi_x_from_file(struct precision* ppr,
                                   struct heating* phe);
  
  /* Efficiency of energy deposition */
  int heating_read_feff_from_file(struct precision* ppr,
                                  struct heating* phe);
 
  /* Heating functions */
  int heating_DM_annihilation(struct heating * phe,
                              double z,
                              double * energy_rate);

  int heating_DM_decay(struct heating * phe,
                       double z,
                       double * energy_rate);

  int heating_add_second_order_terms(struct background* pba,
                                     struct thermo* pth,
                                     struct perturbs* ppt);

#ifdef __cplusplus
}
#endif

#endif
