/********************************************************************************************************/
/*                 HYREC: Hydrogen and Helium Recombination Code                                        */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                                     */
/*                                                                                                      */
/*         hydrogen.h: all functions related to Hydrogen recombination                                  */
/*                                                                                                      */
/*         Units used: cgs + eV (all temperatures in eV)                                                */
/*                                                                                                      */
/*         Version: December 2014                                                                       */
/*                                                                                                      */
/*         Revision history:                                                                            */
/*            - Written November 2010.                                                                  */
/*            - January 2011: updated value of 2s--1s decay rate,                                       */
/*                            changed temperature range for effective rates                             */
/*            - May 2012:   - Using the photon distortion instead of absolute value of radiation field  */
/*                          - Accounting for explicit dependence on alpha and m_e                       */
/*                          - Some definitions moved to header file history.h                           */
/*            - December 2014: - Accounts for additional energy injection                               */
/********************************************************************************************************/
#include "energy_injection.h"
/* Definition of different recombination models  */

#define PEEBLES   0    /* Peebles's effective three-level atom */
#define RECFAST   1    /* Effective three-level atom for hydrogen with fudge factor F = 1.14 */
#define EMLA2s2p  2    /* Correct EMLA model, with standard decay rates from 2s and 2p only (accounts for nmax = infinity, l-resolved) */
#define FULL      3    /* All radiative transfer effects included. Additional switches in header file hydrogen.h */

/* When the probability of being ionized from n=2 becomes lower than PION_MAX,
   switch off radiative transfer calculation as it becomes irrelevant */
#define PION_MAX  1e-2


/****** CONSTANTS IN CGS + EV UNIT SYSTEM *******/

#define EI   13.598286071938324              /* Hydrogen ionization energy in eV, reduced mass, no relativistic corrections */

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

/* Structure for HyRec internal parameters */

typedef struct {
	double h;
  double T0;                                /* CMB temperature today in K*/
  double obh2, omh2, odeh2, okh2, orh2;     /* density parameters */
  double Nnueff;                            /* effective number of neutrinos */
  double fHe;                               /* Helium fraction by number */
  double nH0;                               /* density of hydrogen today in cm^{-3} (changed from m^{-3}) */
  double Y;
  double fsR, meR;              /* fine-structure constant alpha/alpha(today)
                                    and me/me(today) (Added April 2012)*/

  INJ_PARAMS *inj_params;     /* Structure containing all Energy-injection parameters */

} REC_COSMOPARAMS;
/*********** EFFECTIVE 3-LEVEL A LA PEEBLES ***************/
double SAHA_FACT(double fsR, double meR);
double LYA_FACT(double fsR, double meR);
double L2s_rescaled(double fsR, double meR);
void rescale_T(double *T, double fsR, double meR);

double alphaB_PPB(double TM, double fsR, double meR);
double rec_TLA_dxHIIdlna(double z, double xe, double xHII, double nH, double H, double TM, double TR,
			 double Fudge, double fsR, double meR, double dEdtdV, REC_COSMOPARAMS *params);


/************* EFFECTIVE MULTI LEVEL ATOM *******************/

/*** Effective rate tables and associated parameters ***/

#define ALPHA_FILE  "HyRec_2017/data/Alpha_inf.dat"     /* Effective recombination coefficients to 2s and 2p */
#define RR_FILE     "HyRec_2017/data/R_inf.dat"         /* Effective transfer rate R_{2p,2s} */
#define TR_MIN 0.004                         /* Minimum Tr in eV */
#define TR_MAX 0.4                           /* Maximum Tr in eV */
#define NTR    100                           /* Number of Tr values */
#define TM_TR_MIN 0.1                        /* Same thing for ratio Tm/Tr*/
#define TM_TR_MAX 1.0
#define NTM 40

/*** Tables and parameters for radiative transfer calculation ***/

#define TWOG_FILE "HyRec_2017/data/two_photon_tables.dat"
#define NSUBLYA  140
#define NSUBLYB  271
#define NVIRT    311
#define NDIFF    80

#define DLNA     8.49e-5    /* Timestep. Maximum compatible with these tables is 8.49e-5 */

/* Higher-resolution tables  */
/* #define TWOG_FILE "data/two_photon_tables_hires.dat"
/* #define NSUBLYA  408 */
/* #define NSUBLYB  1323 */
/* #define NVIRT    1493 */
/* #define NDIFF    300 */
/* #define DLNA    8.47e-5 */


/**** Structure containing all atomic data for hydrogen ****/

typedef struct {
  /* Tables of effective rates */
  double logTR_tab[NTR];
  double TM_TR_tab[NTM];
  double **logAlpha_tab[2];
  double logR2p2s_tab[NTR];
  double DlogTR, DTM_TR;

  /* Tables of 2-photon rates */
  double Eb_tab[NVIRT];       /* Energies of the virtual levels in eV */
  double A1s_tab[NVIRT];      /* 3*A2p1s*phi(E)*DE */
  double A2s_tab[NVIRT];      /* dLambda_2s/dE * DeltaE if E < Elya dK2s/dE * Delta E if E > Elya */
  double A3s3d_tab[NVIRT];    /* (dLambda_3s/dE + 5*dLambda_3d/dE) * Delta E for E < ELyb, Raman scattering rate for E > ELyb */
  double A4s4d_tab[NVIRT];    /* (dLambda_4s/dE + 5*dLambda_4d/dE) * Delta E */

} HYREC_ATOMIC;


/**** Structure containing all radiative transfer tables ****/

typedef struct {

  double z0;               // first redshift at which radiation fields are stored
  double **Dfminus_hist;
  double **Dfnu_hist;
  double **Dfminus_Ly_hist;

} RADIATION;

void allocate_radiation(RADIATION *rad, long int Nz);
void free_radiation(RADIATION *rad);


void allocate_and_read_atomic(HYREC_ATOMIC *atomic);
void free_atomic(HYREC_ATOMIC *atomic);
void interpolate_rates(double Alpha[2], double DAlpha[2], double Beta[2], double *R2p2s, double TR, double TM_TR,
		       HYREC_ATOMIC *atomic, double fsR, double meR, int *error);
double rec_HMLA_dxHIIdlna(double z, double xe, double xHII, double nH, double H, double TM, double TR,
		          HYREC_ATOMIC *atomic, double fsR, double meR, double dE_dtdV, int *error, REC_COSMOPARAMS *params);
void populate_Diffusion(double *Aup, double *Adn, double *A2p_up, double *A2p_dn,
                        double TM, double Eb_tab[NVIRT], double A1s_tab[NVIRT]);
void populateTS_2photon(double Trr[2][2], double *Trv[2], double *Tvr[2], double *Tvv[3],
                        double sr[2], double sv[NVIRT], double Dtau[NVIRT],double z,
                        double xe, double xHII, double TM, double TR, double nH, double H, HYREC_ATOMIC *atomic,
												REC_COSMOPARAMS *params, double Dfplus[NVIRT], double Dfplus_Ly[],
                        double Alpha[2], double DAlpha[2], double Beta[2],
                        double fsR, double meR, double dE_dtdV, int *error);
void solveTXeqB(double *diag, double *updiag, double *dndiag, double *X, double *B, unsigned N);
void solve_real_virt(double xr[2], double xv[NVIRT], double Trr[2][2], double *Trv[2], double *Tvr[2],
                     double *Tvv[3], double sr[2], double sv[NVIRT]);
double interp_Dfnu(double x0, double dx, double *ytab, unsigned int Nx, double x);
void fplus_from_fminus(double fplus[NVIRT], double fplus_Ly[], double **Dfminus_hist, double **Dfminus_Ly_hist,
                       double TR, double zstart, unsigned iz, double z, double Eb_tab[NVIRT]);
double rec_HMLA_2photon_dxedlna(double xe, double nH, double H, double TM, double TR,
                                HYREC_ATOMIC *atomic,REC_COSMOPARAMS *params,
                                double **Dfminus_hist, double **Dfminus_Ly_hist, double **Dfnu_hist,
                                double zstart, unsigned iz, double z, double fsR, double meR, double dE_dtdV, int *error);
double rec_dxHIIdlna(int model, double xe, double xHII, double nH, double H, double TM, double TR,
                     HYREC_ATOMIC *atomic, RADIATION *rad, unsigned iz, double z,
		     double fsR, double meR, double dEdtdV, int *error, REC_COSMOPARAMS *params);



/************ SWITCHES FOR RADIATIVE TRANSFER. ALL SWITCHES SET TO 1 ARE THE DEFAULT MODEL  ************/

#define EFFECT_A    1    /* 2s-->1s stimulated two-photon decays and non-thermal absorptions */
#define EFFECT_B    1    /* Sub-Lyman alpha two-photon transitions 3s/3d<--> 1s and 4s/4d<-->1s */
#define EFFECT_C    1    /* Super-Lyman alpha two-photon transitions 3s/3d<--> 1s and 4s/4d<-->1s */
#define EFFECT_D    1    /* Raman scattering from 2s and 3s/3d */
#define DIFFUSION   1    /* Lyman alpha frequency diffusion */
