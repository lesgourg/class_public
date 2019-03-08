#ifndef __HEATING__
#define __HEATING__

#include "common.h" //Use here ONLY the things required for defining the struct (i.e. common.h)
struct heating{

  /* Flags */
  int deposit_energy_as;


  double* z_table;
  int z_size;

  int index_ht_CRR;
  int index_ht_BAO;
  int ht_size;

  int has_DM_ann;
  int has_DM_dec;
  int has_pbh_evap;
  int has_pbh_acc;

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

  int heating_at_z(struct background* pba, struct thermo* pth, double z, double* dQdz, double* dxdz, double* pvecback);

  int heating_free(struct thermo* pth);


  /* Own functions */
  int heating_indices(struct thermo* pth);


  /* DM annihilation */
  int heating_DM_annihilation(struct heating * phe,
                            double z,
                            double * energy_rate);

#ifdef __cplusplus
}
#endif

#endif
