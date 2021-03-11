#ifndef __WRAP_HYREC__
#define __WRAP_HYREC__


//Think about include guard
#ifndef TWOG_FILE
#include "helium.h"
#include "hydrogen.h"
#include "history.h"
#endif

#include "common.h" //Use here ONLY the things required for defining the struct (i.e. common.h for the ErrorMsg)

struct thermohyrec{

  HYREC_DATA *data;

  double zstart;
  double zend;

  double xHeII_limit;

  ErrorMsg error_message;

  int thermohyrec_verbose;
};

/**************************************************************/

/* Boilerplate for C++ */
#ifdef __cplusplus
extern "C" {
#endif

  int thermodynamics_hyrec_init(struct precision* ppr, struct background * pba, struct thermodynamics * pth, double Nnow, double T_cmb, double fHe, double zstart_hyrec, struct thermohyrec* phy);

  int thermodynamics_hyrec_calculate_xe(struct thermodynamics * pth, struct thermohyrec * phy,
                                        double z, double H_in, double T_b, double T_gamma,
                                        double* x_e, double* dxe_dlna);

  int thermodynamics_hyrec_get_xe(struct thermohyrec * phy, double z, double* x_e, double* dxdlna);

  int thermodynamics_hyrec_free(struct thermohyrec* phy);

  int hyrec_dx_H_dz(struct thermodynamics* pth, struct thermohyrec* phy, double x_H, double x_He, double xe, double nH, double z, double Hz, double Tmat, double Trad, double *dx_H_dz);
  int hyrec_dx_He_dz(struct thermodynamics* pth, struct thermohyrec* phy, double x_H, double x_He, double xe, double nH, double z, double Hz, double Tmat, double Trad, double *dx_He_dz);

#ifdef __cplusplus
}
#endif

/* Ionization energy of hydrogen HI */
#define _E_H_ion_   13.5984336478
/* Lyman-Alpha transition energy of hydrogen , approximately 3/4 of the ionization energy because (1s->2p transition and E~1/n^2) */
#define _E_H_lya_   10.1988356821

#endif
