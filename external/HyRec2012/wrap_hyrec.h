#ifndef __WRAP_HYREC__
#define __WRAP_HYREC__


//Think about include guard
#ifndef TWOG_FILE
#include "hyrectools.h"
#include "helium.h"
#include "hydrogen.h"
#include "hyrec_params.h"
#endif

#include "common.h" //Use here ONLY the things required for defining the struct (i.e. common.h for the ErrorMsg)

struct thermohyrec{

  HRATEEFF* rate_table;
  TWO_PHOTON_PARAMS *twog_params;

  double ** Dfnu_hist;
  double ** Dfminus_hist;
  double *  Dfminus_Ly_hist[3];
  int N_TR; //Number of T_radiation values in precomputed files
  int N_TM; //Number of T_matter values in precomputed files
  int N_LY;
  int N_VIRT;
  short stage;

  double T_cmb;
  double nH0;
  double T0;

  double zstart;
  double zend;
  double dlna;
  long nz;


  double fsR;
  double meR;
  double fHe;

  int izH0;
  double zH0;

  double* xe_output;
  short to_store;

  int filled_until_index_z;
  double z_prev;
  double H_prev;
  double TR_prev;
  double TM_prev;
  double ion_prev;
  double exclya_prev;

  double xHeIII;
  double xHeIII_limit;

  double xHeII;
  double xH1s;

  double dxHIIdlna;
  double dxHIIdlna_prev;
  double dxHIIdlna_prev2;
  double dxHeIIdlna;
  double dxHeIIdlna_prev;
  double dxHeIIdlna_prev2;

  short saha_flag;

  FileName alpha_file;
  FileName rr_file;
  FileName twog_file;
  ErrorMsg error_message;

  int thermohyrec_verbose;
};

/**************************************************************/

/* *
 * Putting this down here is important, because of the special nature of this wrapper
 * */

struct thermo* pth;

/* Boilerplate for C++ */
#ifdef __cplusplus
extern "C" {
#endif

  int thermodynamics_hyrec_init(struct precision* ppr, double Nnow, double T_cmb, double fHe, double zstart_hyrec, struct thermohyrec* phy);

  int thermodynamics_hyrec_calculate_xe(struct thermo * pth, struct thermohyrec * phy,
                                        double z, double H_in, double T_b, double T_gamma,
                                        double* x_e, double* dxe_dlna);

  int thermodynamics_hyrec_get_xe(struct thermohyrec * phy, double z, double* x_e, double* dxdlna);

  int thermodynamics_hyrec_free(struct thermohyrec* phy);

#ifdef __cplusplus
}
#endif

/* Ionization energy of hydrogen HI */
#define _E_H_ion_   13.5984336478
/* Lyman-Alpha transition energy of hydrogen , approximately 3/4 of the ionization energy because (1s->2p transition and E~1/n^2) */
#define _E_H_lya_   10.1988356821

#endif
