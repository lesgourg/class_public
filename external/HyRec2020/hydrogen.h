/*******************************************************************************************************/
/*                 HYREC-2: Hydrogen and Helium Recombination Code                                      */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (2010-17)                                     */
/*             with contributions from Nanoom Lee (2020)                                                */
/*                                                                                                      */
/*         hydrogen.h: all functions related to Hydrogen recombination                                  */
/*                                                                                                      */
/*         Units used: cgs + eV (all temperatures in eV)                                                */
/*                                                                                                      */
/*         Revision history:                                                                            */
/*            - January 2020 : - Added new mode, SWIFT                                                  */
/*                             - Two timestep parameters for FULL mode and other modes.                 */
/*            - December 2014: - Accounts for additional energy injection                               */
/*            - May 2012:   - Using the photon distortion instead of absolute value of radiation field  */
/*                          - Accounting for explicit dependence on alpha and m_e                       */
/*                          - Some definitions moved to header file history.h                           */
/*            - January 2011: updated value of 2s--1s decay rate,                                       */
/*                            changed temperature range for effective rates                             */
/*            - Written November 2010.                                                                  */
/********************************************************************************************************/

#ifndef __HYDROGEN__
#define __HYDROGEN__

#include "energy_injection.h"

/* Definition of different recombination models  */

#define PEEBLES   0    /* Peebles's effective three-level atom */
#define RECFAST   1    /* Effective three-level atom for hydrogen with fudge factor F = 1.14 */
#define EMLA2s2p  2    /* Correct EMLA model, with standard decay rates from 2s and 2p only (accounts for nmax = infinity, l-resolved) */
#define FULL      3    /* All radiative transfer effects included. Additional switches in header file hydrogen.h */
#define SWIFT     4    /* Fast calculation with fitting function which is calculated based on FULL mode */


/* When the probability of being ionized from n=2 becomes lower than PION_MAX,
   switch off radiative transfer calculation as it becomes irrelevant */
#define PION_MAX  1e-2


/****** CONSTANTS IN CGS + EV UNIT SYSTEM *******/

#define EI   13.598286071938324        /* Hydrogen ionization energy in eV, reduced mass, no relativistic corrections */

/* Energy differences between excited levels of hydrogen -- used often */
#define E21  10.198714553953742
#define E31  12.087365397278509
#define E41  12.748393192442178
#define E32  1.8886508433247664
#define E42  2.5496786384884356

#define hPc       1.239841874331e-04   /* hc in eV cm */
#define mH        0.93878299831e9      /* Hydrogen atom mass in eV/c^2 */
#define kBoltz    8.617343e-5          /* Boltzmann constant in eV/K */
#define L2s1s     8.2206               /* 2s -> 1s two-photon decay rate in s^{-1} (Labzowsky et al 2005) */


/************* EFFECTIVE MULTI LEVEL ATOM *******************/

/*** Effective rate tables and associated parameters ***/

#define ALPHA_FILE  "Alpha_inf.dat"                 /* Effective recombination coefficients to 2s and 2p */
#define RR_FILE     "R_inf.dat"                     /* Effective transfer rate R_{2p,2s} */
#define TR_MIN 0.004                                /* Minimum Tr in eV */
#define TR_MAX 0.4                                  /* Maximum Tr in eV */
#define NTR    100                                  /* Number of Tr values */
#define T_RATIO_MIN 0.1                             /* T_RATIO is min(Tm/Tr, Tr/Tm) */
#define T_RATIO_MAX 1.0
#define NTM 40

/************* CORRECTION FUNCTION AND ITS FIRST DERIVATIVE FOR SWIFT MODE *****/
#define FIT_FILE  "fit_swift.dat"              /* Correction function and first derivative for SWIFT mode */
#define DKK_SIZE 265

/*** Tables and parameters for radiative transfer calculation ***/
#define TWOG_FILE "two_photon_tables.dat"

#define NSUBLYA  140
#define NSUBLYB  271
#define NVIRT    311
#define NDIFF    80

#define DLNA_HYREC   8.49e-5    /* Timestep for FULL mode. Maximum compatible with these tables is 8.49e-5 */
//#define DLNA_HYREC   2.e-6    /* Timestep used in FULL mode for SWIFT correction function calculation*/
#define DLNA_SWIFT   4.e-3      /* Timestep for any other mode.*/

#define SIZE_ErrorM   2048
#define SIZE_InputFile   512

/* Higher-resolution tables  */
/*
#define TWOG_FILE_CLASS "two_photon_tables_hires.dat"
#define NSUBLYA  408
#define NSUBLYB  1323
#define NVIRT    1493
#define NDIFF    300
#define DLNA    1.e-7
*/

/**** Structure containing all atomic data for hydrogen ****/

typedef struct {
  /* Tables of effective rates */
  double logTR_tab[NTR];
  double T_RATIO_tab[NTM];
  double **logAlpha_tab[4];
  double logR2p2s_tab[NTR];
  double DlogTR, DT_RATIO;

  /* Tables of 2-photon rates */
  double Eb_tab[NVIRT];       /* Energies of the virtual levels in eV */
  double A1s_tab[NVIRT];      /* 3*A2p1s*phi(E)*DE */
  double A2s_tab[NVIRT];      /* dLambda_2s/dE * DeltaE if E < Elya dK2s/dE * Delta E if E > Elya */
  double A3s3d_tab[NVIRT];    /* (dLambda_3s/dE + 5*dLambda_3d/dE) * Delta E for E < ELyb, Raman scattering rate for E > ELyb */
  double A4s4d_tab[NVIRT];    /* (dLambda_4s/dE + 5*dLambda_4d/dE) * Delta E */

} HYREC_ATOMIC;

typedef struct {
  double *swift_func[5];
} FIT_FUNC;


/**** Structure containing all radiative transfer tables ****/

typedef struct {

  double z0;               // first redshift at which radiation fields are stored
  long iz_rad_0;
  double **Dfminus_hist;
  double **Dfnu_hist;
  double **Dfminus_Ly_hist;

} RADIATION;


/* Structure for HYREC-2 internal parameters */

typedef struct {
  double h;                                         /* Hubble constant */
  double T0;                                        /* CMB temperature today in K*/
  double obh2, ocbh2, odeh2, okh2, orh2, onuh2;     /* density parameters */
  double w0, wa;                                    /* Dark energy equation of state parameters */
  double Neff;                                      /* total effective number of neutrinos (massive + massless) */
  double Nur;                                       /* number of massless neutrinos */
  double Nmnu;                                      /* number of massive neutrinos */
  double mnu[3];                                    /* neutrino masses */
  double fHe;                                       /* Helium fraction by number */
  double nH0;                                       /* density of hydrogen today in cm^{-3} [Changed from m^{-3} in February 2015] */
  double YHe;                                       /* Helium fraction */
  double fsR, meR;                                  /* fine-structure constant alpha/alpha(today)
                                                       and me/me(today) (Added April 2012)*/
  double dlna, nz;

  INJ_PARAMS *inj_params;                           /* Structure containing all Energy-injection parameters */

} REC_COSMOPARAMS;

/* Structure for HYREC-2 data */

typedef struct{
  HYREC_ATOMIC *atomic;
  REC_COSMOPARAMS *cosmo;
  double zmax;
  double zmin;
  long int Nz;
  double *xe_output;
  double *Tm_output;
  int error;
  int quasi_eq;
  char *error_message;
  char *path_to_hyrec;
  RADIATION *rad;
  FIT_FUNC *fit;
} HYREC_DATA;

/*********** EFFECTIVE 3-LEVEL A LA PEEBLES ***************/
double SAHA_FACT(double fsR, double meR);
double LYA_FACT(double fsR, double meR);
double L2s_rescaled(double fsR, double meR);
void rescale_T(double *T, double fsR, double meR);

double alphaB_PPB(double TM, double fsR, double meR);
double rec_TLA_dxHIIdlna(REC_COSMOPARAMS *cosmo, double xe, double xHII, double nH, double H, double TM, double TR, double Fudge);



void allocate_radiation(RADIATION *rad, long int Nz, int *error, char error_message[SIZE_ErrorM]);
void free_radiation(RADIATION *rad);


void allocate_and_read_atomic(HYREC_ATOMIC *atomic, int *error, char *path_to_hyrec, char error_message[SIZE_ErrorM]);
void free_atomic(HYREC_ATOMIC *atomic);
void allocate_and_read_fit(FIT_FUNC *fit, int *error, char *path_to_hyrec, char error_message[SIZE_ErrorM]);
void free_fit(FIT_FUNC *fit);
void interpolate_rates(double Alpha[2], double DAlpha[2], double Beta[2], double *R2p2s, double TR, double TM_TR,
                       HYREC_ATOMIC *atomic, double fsR, double meR, int *error, char error_message[SIZE_ErrorM]);
double rec_swift_hyrec_dxHIIdlna(HYREC_DATA *data, double xe, double xHII, double nH, double Hubble, double TM, double TR, double z);
double rec_HMLA_dxHIIdlna(HYREC_DATA *data, double xe, double xHII, double nH, double H, double TM, double TR);
void populate_Diffusion(double *Aup, double *Adn, double *A2p_up, double *A2p_dn,
                        double TM, double Eb_tab[NVIRT], double A1s_tab[NVIRT]);
void populateTS_2photon(double Trr[2][2], double *Trv[2], double *Tvr[2], double *Tvv[3],
                        double sr[2], double sv[NVIRT], double Dtau[NVIRT],
                        double xe, double xHII, double TM, double TR, double nH, double H, HYREC_ATOMIC *atomic,
                        double Dfplus[NVIRT], double Dfplus_Ly[],
                        double Alpha[2], double DAlpha[2], double Beta[2],
                        double fsR, double meR, double exclya, int *error, char error_message[SIZE_ErrorM]);
void solveTXeqB(double *diag, double *updiag, double *dndiag, double *X, double *B, unsigned N, int *error, char error_message[SIZE_ErrorM]);
void solve_real_virt(double xr[2], double xv[NVIRT], double Trr[2][2], double *Trv[2], double *Tvr[2],
                     double *Tvv[3], double sr[2], double sv[NVIRT], int *error, char error_message[SIZE_ErrorM]);
double interp_Dfnu(double x0, double dx, double *ytab, unsigned int Nx, double x);
void fplus_from_fminus(double fplus[NVIRT], double fplus_Ly[], double **Dfminus_hist, double **Dfminus_Ly_hist,
                       double TR, double zstart, unsigned iz, double z, double Eb_tab[NVIRT]);
double rec_HMLA_2photon_dxedlna(HYREC_DATA *data, double xe, double nH, double H, double TM, double TR, unsigned iz, double z);
double rec_dxHIIdlna(HYREC_DATA *data, int model, double xe, double xHII, double nH, double H, double TM, double TR,
                     unsigned iz, double z);


/************ SWITCHES FOR RADIATIVE TRANSFER. ALL SWITCHES SET TO 1 ARE THE DEFAULT MODEL  ************/

#define EFFECT_A    1    /* 2s-->1s stimulated two-photon decays and non-thermal absorptions */
#define EFFECT_B    1    /* Sub-Lyman alpha two-photon transitions 3s/3d<--> 1s and 4s/4d<-->1s */
#define EFFECT_C    1    /* Super-Lyman alpha two-photon transitions 3s/3d<--> 1s and 4s/4d<-->1s */
#define EFFECT_D    1    /* Raman scattering from 2s and 3s/3d */
#define DIFFUSION   1    /* Lyman alpha frequency diffusion */

#endif
