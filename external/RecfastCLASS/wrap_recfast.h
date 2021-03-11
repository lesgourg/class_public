#ifndef __RECFAST_CLASS__
#define __RECFAST_CLASS__

#include "common.h" //Use here ONLY the things required for defining the struct (i.e. common.h for the ErrorMsg)

// By default, recfast assumes photo-ionization coefficients depending on Tmat
// However, this is an approximation (see e.g. arXiV:1605.03928 page 10, arXiV:1503.04827 page 2, right column)
// We thus also give the user the ability to use the dependence on Tgamma = T_photon(z) = Tcmb * (1+z) instead
enum recfast_photoion_modes {recfast_photoion_Tmat, recfast_photoion_Trad};

struct thermorecfast {

  double CDB;     /**< defined as in RECFAST */
  double CR;      /**< defined as in RECFAST */
  double CK;      /**< defined as in RECFAST */
  double CL;      /**< defined as in RECFAST */
  double fHe;     /**< defined as in RECFAST */
  double YHe;
  double CDB_He;  /**< defined as in RECFAST */
  double CK_He;   /**< defined as in RECFAST */
  double CL_He;   /**< defined as in RECFAST */
  double CL_Het;
  double H_frac;  /**< defined as in RECFAST */
  double Tnow;    /**< defined as in RECFAST */
  double Nnow;    /**< defined as in RECFAST */
  double CDB_He2s2p;   /**< Bfact in RECFAST */
  double CB1;     /**< defined as in RECFAST */
  double CB1_He1; /**< defined as in RECFAST */
  double CB1_He2; /**< defined as in RECFAST */

  double H0;

  double AGauss1;
  double AGauss2;
  double zGauss1;
  double zGauss2;
  double wGauss1;
  double wGauss2;
  int Heswitch;
  int Hswitch;
  double x_H0_trigger2;
  double x_He0_trigger2;
  double z_switch_late;

  double x_He_trigger_small;
  double fudge_He;
  double fudge_H;
  double x_H_limit_KHe;
  double x_H_limit_CfHe_t;

  double max_exp_boltz;
  double Bfact;
  double CT;

  enum recfast_photoion_modes photoion_mode;

  ErrorMsg error_message;

};


/**************************************************************/

/* Boilerplate for C++ */
#ifdef __cplusplus
extern "C" {
#endif

  int recfast_init(struct precision* ppr,
                   struct background* pba,
                   struct thermodynamics * pth,
                   struct thermorecfast * precfast,
                   enum recfast_photoion_modes recfast_photoion_mode,
                   double fHe);

  int recfast_dx_H_dz(struct thermodynamics* pth, struct thermorecfast * pre,
                      double x_H, double x, double n,
                      double z, double Hz, double Tmat, double Trad,
                      double* dxH_dz);

  int recfast_dx_He_dz(struct thermodynamics* pth, struct thermorecfast * pre,
                       double x_He, double x, double x_H, double n,
                       double z, double Hz, double Tmat, double Trad,
                       double* dxHe_dz);
#ifdef __cplusplus
}
#endif

/**
 * @name Some specific constants needed by RECFAST:
 */

//@{

#define _Lambda_            8.2245809
#define _Lambda_He_         51.3
/* Ionization inv wavenumber of hydrogen in 1/m */
#define _L_H_ion_           1.096787737e7
/* Lyman-alpha transition inv wavenumber of hydrogen in 1/m, approximately 3/4 of the ionization because (1s->2p transition and E~1/lambda~1/n^2) */
#define _L_H_alpha_         8.225916453e6
/* Ionization inv wavenumber of helium HeI in eV */
#define _L_He1_ion_         1.98310772e7
/* Ionization inv wavenumber of helium HeII in eV */
#define _L_He2_ion_         4.389088863e7
/* Inv Wavenumber of 1s->2s transition of helium HeI in eV */
#define _L_He_2s_           1.66277434e7
/* Inv Wavenumber of 1s->2p transition of helium HeI in eV */
#define _L_He_2p_           1.71134891e7

#define _A2P_s_             1.798287e9     /*updated like in recfast 1.4*/
#define _A2P_t_             177.58e0       /*updated like in recfast 1.4*/
#define _L_He_2Pt_          1.690871466e7  /*updated like in recfast 1.4*/
#define _L_He_2St_          1.5985597526e7 /*updated like in recfast 1.4*/
#define _L_He2St_ion_       3.8454693845e6 /*updated like in recfast 1.4*/
#define _sigma_He_2Ps_      1.436289e-22   /*updated like in recfast 1.4*/
#define _sigma_He_2Pt_      1.484872e-22   /*updated like in recfast 1.4*/
//@}

/**
 * @name Some specific constants needed by recfast_derivs:
 */

//@{

#define _a_PPB_           4.309
#define _b_PPB_          -0.6166
#define _c_PPB_           0.6703
#define _d_PPB_           0.5300
#define _T_0_             pow(10.,0.477121)   /* from recfast 1.4 */
#define _a_VF_            pow(10.,-16.744)
#define _b_VF_            0.711
#define _T_1_             pow(10.,5.114)
#define _a_trip_          pow(10.,-16.306) /* from recfast 1.4 */
#define _b_trip_          0.761            /* from recfast 1.4 */

//@}

#endif
