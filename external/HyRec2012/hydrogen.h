/********************************************************************************************************/
/*                 HYREC: Hydrogen and Helium Recombination Code                                        */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                                     */
/*                                                                                                      */
/*         hydrogen.h: all functions related to Hydrogen recombination                                  */
/*                                                                                                      */
/*         Units used: cgs + eV (all temperatures in eV)                                                */
/*                                                                                                      */
/*         Version: May 2012                                                                            */
/*                                                                                                      */
/*         Revision history:                                                                            */
/*            - Written November 2010.                                                                  */
/*            - January 2011: updated value of 2s--1s decay rate,                                       */
/*                            changed temperature range for effective rates                             */
/*            - May 2012:   - Using the photon distortion instead of absolute value of radiation field  */
/*                          - Accounting for explicit dependence on alpha and m_e                       */
/*                          - Some definitions moved to header file history.h                           */                   
/********************************************************************************************************/ 

#include "hyrec_params.h"  /* Needed for definition of switches for derivatives, and numerical integration parameters */
                           /* The two-photon table and parameter definition were moved to history.h (May 2012) */


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



double square(double x);
double cube(double x);
 
/*********** PEEBLES + POST-SAHA + RECFAST ***************/ 
double SAHA_FACT(double fsR, double meR);
double LYA_FACT(double fsR, double meR);
double L2s_rescaled(double fsR, double meR);
void rescale_T(double *T, double fsR, double meR);

double xeSaha(double nH, double TR, double fsR, double meR);
double dxeSaha_dlna(double nH, double TR, double fsR, double meR);
double alphaB_PPB(double TM, double fsR, double meR);
double rec_TLA_dxedlna(double xe, double nH, double H, double TM, double TR, double F, double fsR, double meR);
double xe_PostSaha_TLA(double nH, double H, double T, double F, double *Dxe, double fsR, double meR);

/************* EFFECTIVE MULTI LEVEL ATOM *******************/

#define ALPHA_FILE  "Alpha_inf.dat"     /* Contains the effective recombination coefficients to 2s and 2p */
#define RR_FILE     "R_inf.dat"         /* Contains the effective transfer rate R_{2p,2s} */


/* Boundaries and number of elements of temperature tables */
#define TR_MIN 0.004            /* Tr parameters */
#define TR_MAX 0.4        
#define NTR    100                        
#define TM_TR_MIN 0.1           /* Tm/Tr parameters */
#define TM_TR_MAX 1.0
#define NTM 40                       

/**** Effective rate coefficients structure ****/
typedef struct {
  double *logTR_tab;
  double *TM_TR_tab;
  double **logAlpha_tab[2];
  double *logR2p2s_tab;
  double DlogTR, DTM_TR;
} HRATEEFF;

void read_rates(HRATEEFF *rate_table);
void interpolate_rates(double Alpha[2], double DAlpha[2], double Beta[2], double *R2p2s, double TR, double TM_TR, HRATEEFF *rate_table, double fsR, double meR);
double rec_HMLA_dxedlna(double xe, double nH, double H, double TM, double TR, HRATEEFF *rate_table, double fsR, double meR);
double xe_PostSaha_HMLA(double nH, double H, double T, HRATEEFF *rate_table, double *Dxe, double fsR, double meR);

/************ TWO-PHOTON PROCESSES AND DIFFUSION. ALL SWITCHES SET TO 1 ARE THE DEFAULT MODEL  ************/

#define EFFECT_A    1    /* 2s-->1s stimulated two-photon decays and non-thermal absorptions */
#define EFFECT_B    1    /* Sub-Lyman alpha two-photon transitions 3s/3d<--> 1s and 4s/4d<-->1s */
#define EFFECT_C    1    /* Super-Lyman alpha two-photon transitions 3s/3d<--> 1s and 4s/4d<-->1s */
#define EFFECT_D    1    /* Raman scattering from 2s and 3s/3d */
#define DIFFUSION   1    /* Lyman alpha frequency diffusion */

/**** Structure containing two-photon transition tables ****/
 
typedef struct {     
    double Eb_tab[NVIRT];       /* Energies of the virtual levels in eV */
    double A1s_tab[NVIRT];      /* 3*A2p1s*phi(E)*DE */ 
    double A2s_tab[NVIRT];      /* dLambda_2s/dE * DeltaE if E < Elya dK2s/dE * Delta E if E > Elya */
    double A3s3d_tab[NVIRT];    /* (dLambda_3s/dE + 5*dLambda_3d/dE) * Delta E for E < ELyb, Raman scattering rate for E > ELyb */
    double A4s4d_tab[NVIRT];    /* (dLambda_4s/dE + 5*dLambda_4d/dE) * Delta E */
}  TWO_PHOTON_PARAMS;


void read_twog_params(TWO_PHOTON_PARAMS *twog);
void populate_Diffusion(double *Aup, double *Adn, double *A2p_up, double *A2p_dn, 
                        double TM, double Eb_tab[NVIRT], double A1s_tab[NVIRT]);
void populateTS_2photon(double Trr[2][2], double *Trv[2], double *Tvr[2], double *Tvv[3], 
                        double sr[2], double sv[NVIRT], double Dtau[NVIRT],
                        double xe, double xHII, double TM, double TR, double nH, double H, HRATEEFF *rate_table,
                        TWO_PHOTON_PARAMS *twog, double Dfplus[NVIRT], double Dfplus_Ly[], 
                        double Alpha[2], double DAlpha[2], double Beta[2], double fsR, double meR);
void solveTXeqB(double *diag, double *updiag, double *dndiag, double *X, double *B, unsigned N);
void solve_real_virt(double xr[2], double xv[NVIRT], double Trr[2][2], double *Trv[2], double *Tvr[2], 
                     double *Tvv[3], double sr[2], double sv[NVIRT]);
double interp_Dfnu(double x0, double dx, double *ytab, unsigned int Nx, double x);
void fplus_from_fminus(double fplus[NVIRT], double fplus_Ly[], double **Dfminus_hist, double *Dfminus_Ly_hist[], 
                       double TR, double zstart, unsigned iz, double z, double Eb_tab[NVIRT]);
double rec_HMLA_2photon_dxedlna(double xe, double nH, double H, double TM, double TR,
                                HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog,
                                double **Dfminus_hist, double *Dfminus_Ly_hist[], double **Dfnu_hist,
                                double zstart, unsigned iz, double z, double fsR, double meR);
double rec_dxHIIdlna(int model, double xe, double xHII, double nH, double H, double TM, double TR, 
                     HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog, double **Dfminus_hist, double *Dfminus_Ly_hist[], 
                     double **Dfnu_hist, double zstart, unsigned iz, double z, double fsR, double meR);
