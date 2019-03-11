#ifndef __HEATING__
#define __HEATING__

#include "common.h" //Use here ONLY the things required for defining the struct (i.e. common.h)
struct heating{

  /* Flags */
  int deposit_energy_as;
  int has_exotic_injection;
  int has_dcdm;

  double* z_table;
  int z_size;
  int last_index_z_dep;
  int last_index_z_inj;
  int last_index_z_feff;
  int last_index_chix;
  int last_index_chiz;

  int filled_until_index_z_dep;
  double filled_until_z_dep;
  int filled_until_index_z_inj;
  double filled_until_z_inj;

  double* pvecdeposition;
  double tol_z_table;

  double* chiz_table;
  int chiz_size;
  double* chix_table;
  int chix_size;

  int feff_z_size;
  double* feff_table;

  int to_store;

  int index_ht_CRR;
  int index_ht_BAO;
  int ht_size;

  int has_DM_ann;
  int has_DM_dec;
  int has_BH_evap;
  int has_BH_acc;

  double* injection_table;
  int index_inj_BH_evap;
  int index_inj_BH_acc;
  int index_inj_DM_ann;
  int index_inj_DM_dec;
  int index_inj_BAO;
  int index_inj_CRR;
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
  //int index_dep_lowE;
  int dep_size;


  /* Background stuff etc. */
  double H0;
  double rho_crit0;
  double Omega0_cdm;
  double rho_cdm;
  double rho_dcdm;
  double t;
  double Gamma_dcdm;
  int last_index_bg;


  double f_eff;

  /* Parameters */
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

  short has_on_the_spot;         /**< flag to specify if we want to use the on-the-spot approximation **/

  double decay;                  /**< parameter describing CDM decay (f/tau, see e.g. 1109.6322)*/
  double decay_fraction;

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

  /* Outward functions */
  int heating_init(struct precision * ppr, struct background* pba, struct thermo* pth);

  int heating_at_z(struct background* pba, struct thermo* pth, double x, double z, double* pvecback);

  int heating_free(struct thermo* pth);

  /* Own functions */
  int heating_indices(struct thermo* pth);

  int heating_deposition_function(struct heating* phe, double x, double z);

  int heating_energy_injection_at_z(struct heating* phe, double z, double* dEdz_inj);

  int heating_deposit_analytical_integral(struct background* pba, struct thermo* pth, double z, double* energy_rate);

  int heating_read_feff_from_file(struct precision* ppr, struct heating* phe);

  /* DM annihilation */
  int heating_DM_annihilation(struct heating * phe,
                              double z,
                              double * energy_rate);

  int heating_DM_decay(struct heating * phe,
                       double z,
                       double * energy_rate);
#ifdef __cplusplus
}
#endif

#endif
