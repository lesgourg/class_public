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

  double n_H_now;
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

  // Purely for my own convenience
  double* xe_output;


  int filled_until_index_z;
  double z_prev;
  double TR_prev;
  double TM_prev;

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

/* Boilerplate for C++ */
#ifdef __cplusplus
extern "C" {
#endif

  int thermodynamics_hyrec_init(struct precision* ppr, double Nnow, double T_cmb, double fHe, struct thermohyrec* phy);

  int thermodynamics_hyrec_get_xe(struct thermohyrec * phy,
                                  double z, double H, double T_b, double T_gamma,
                                  double* x_e, double* dxe_dlna, double energy_injection);

  int thermodynamics_hyrec_free(struct thermohyrec* phy);

#ifdef __cplusplus
}
#endif
#endif
