/*************************************************************************************************/
/*                 HYREC: Hydrogen and Helium Recombination Code                                 */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                              */
/*                                                                                               */
/*         hydrogen.c: all functions related to Hydrogen recombination                           */
/*                                                                                               */
/*         Units used: cgs + eV (all temperatures in eV)                                         */
/*                                                                                               */
/*         Version: January 2011                                                                 */
/*         Revision history:                                                                     */
/*            - written November 2010                                                            */
/*            - January 2011: - changed the post-Saha expansion to use the full derivative       */
/*                       (including two-photon processes and diffusion) rather than Peebles'ODE  */
/*	                      - post-Saha expansion can now pass the difference from Saha value  */
/*                                  to external routines                                         */
/*                            - differential 2s--1s rate is now systematically normalized to     */
/*                                 total 2s--1s rate that can be set by user in hydrogen.h       */
/*************************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "hyrectools.h"
#include "hydrogen.h"


/**************************************************************************************************
Case-B recombination coefficient, fit of Pequignot et al 1991, in cm^3 s^{-1}
***************************************************************************************************/

double alphaB_PPB(double TM) {
   double t4;
   t4 = TM/kBoltz/1e4;
   return 4.309e-13*pow(t4,-0.6166)/(1.+ 0.6703*pow(t4,0.5300));
}

/**************************************************************************************************
Peebles recombination rate
***************************************************************************************************/

double rec_HPeebles_dxedlna(double xe, double nH, double H, double TM, double TR, double energy_rate) {

  double RLya, alphaB, four_betaB, C;
  double chi_ion_H;

  RLya   = 4.662899067555897e15 * H / nH / (1.-xe);
  alphaB = alphaB_PPB(TM);
  four_betaB  = 3.016103031869581e21 *TR*sqrt(TR) *exp(-0.25*EI/TR) *alphaB;

  C = (3.*RLya + L2s1s)/(3.*RLya + L2s1s + four_betaB);

  // chi_ion_H = (1.-xe)/3.; // old approximation from Chen and Kamionkowski
  if (xe < 1.)
    chi_ion_H = 0.369202*pow(1.-pow(xe,0.463929),1.70237); // coefficient as revised by Galli et al. 2013 (in fact it is a fit by Vivian Poulin of columns 1 and 4 in Table V of Galli et al. 2013)
  else
    chi_ion_H = 0.;

  return (-nH*xe*xe*alphaB + four_betaB*(1.-xe)*exp(-E21/TR))*C/H
    +chi_ion_H/nH*energy_rate*(1./EI+(1.-C)/E21)/H;

}

/****************************************************************************************************
RecFast recombination rate (Seager et al 1999, 2000): effective three-level atom
with a fudge factor F = 1.14
****************************************************************************************************/

double rec_HRecFast_dxedlna(double xe, double nH, double H, double TM, double TR, double energy_rate) {

  double RLya, alphaB, four_betaB, C;
  double chi_ion_H;

  RLya   = 4.662899067555897e15 * H / nH / (1.-xe);
  alphaB = 1.14 * alphaB_PPB(TM);
  four_betaB  = 3.016103031869581e21 *TR*sqrt(TR) *exp(-0.25*EI/TR) *alphaB;

  C = (3.*RLya + L2s1s)/(3.*RLya + L2s1s + four_betaB);

  //chi_ion_H = (1.-xe)/3.; // old approximation from Chen and Kamionkowski
  if (xe < 1.)
    chi_ion_H = 0.369202*pow(1.-pow(xe,0.463929),1.70237); // coefficient as revised by Galli et al. 2013 (in fact it is a fit by Vivian Poulin of columns 1 and 4 in Table V of Galli et al. 2013)
  else
    chi_ion_H = 0.;

  return (-nH*xe*xe*alphaB + four_betaB*(1.-xe)*exp(-E21/TR))*C/H
    +chi_ion_H/nH*energy_rate*(1./EI+(1.-C)/E21)/H;

}

/**********************************************************************************************
Store tabulated temperatures and read effective MLA rates from files.
**********************************************************************************************/

void read_rates(HRATEEFF *rate_table){

   FILE *fA = fopen(ALPHA_FILE, "r");
   FILE *fR = fopen(RR_FILE, "r");

   unsigned i, j, l;
   int fscanf_result;

   maketab(log(TR_MIN), log(TR_MAX), NTR, rate_table->logTR_tab);
   maketab(TM_TR_MIN, TM_TR_MAX, NTM, rate_table->TM_TR_tab);
   rate_table->DlogTR = rate_table->logTR_tab[1] - rate_table->logTR_tab[0];
   rate_table->DTM_TR = rate_table->TM_TR_tab[1] - rate_table->TM_TR_tab[0];

   for (i = 0; i < NTR; i++) {
      for (j = 0; j < NTM; j++) {
	 for (l = 0; l <= 1; l++) {
           fscanf_result = fscanf(fA, "%le", &(rate_table->logAlpha_tab[l][j][i]));
           if(fscanf_result != 1){printf("Hyrec Warning :: Could not read log Alpha table (Alpha_inf.dat)");}
           rate_table->logAlpha_tab[l][j][i] = log(rate_table->logAlpha_tab[l][j][i]);
        }
      }

      fscanf_result = fscanf(fR, "%le", &(rate_table->logR2p2s_tab[i]));
      if(fscanf_result != 1){printf("Hyrec Warning :: Could not read rate table (R_inf.dat)");}
      rate_table->logR2p2s_tab[i] = log(rate_table->logR2p2s_tab[i]);

   }
   fclose(fA);
   fclose(fR);
}

/************************************************************************************************
Interpolation of tabulated effective rates
To be more (slightly) efficient, not using the external interpolation routine.
************************************************************************************************/

void interpolate_rates(double Alpha[2], double Beta[2], double *R2p2s, double TR, double TM_TR, HRATEEFF *rate_table) {
    double factor;
    unsigned l, k;
    long iTM, iTR;
    double frac1, frac2;
    double logTR;
    double coeff1[4], coeff2[4], temp[4];

    logTR = log(TR);

    /* Check if TM/TR is in the range tabulated */
    if (TM_TR < TM_TR_MIN || TM_TR > TM_TR_MAX) {
      fprintf(stderr, "Error: TM/TR-value is out of range in interpolate_rates.\n");
      exit(1);
    }

    /* Check if log(TR) is in the range tabulated */
    if (TR < TR_MIN || TR > TR_MAX) {
      fprintf(stderr,"Error: TR-value is out of range in interpolate_rates.\n");
      exit(1);
    }

    /* Identify location to interpolate in TM/TR */
    iTM = (long)floor((TM_TR - TM_TR_MIN)/rate_table->DTM_TR);
    if (iTM < 1) iTM = 1;
    if (iTM > NTM-3) iTM = NTM-3;
    frac1 = (TM_TR - TM_TR_MIN)/rate_table->DTM_TR - iTM;
    coeff1[0] = frac1*(frac1-1.)*(2.-frac1)/6.;
    coeff1[1] = (1.+frac1)*(1.-frac1)*(2.-frac1)/2.;
    coeff1[2] = (1.+frac1)*frac1*(2.-frac1)/2.;
    coeff1[3] = (1.+frac1)*frac1*(frac1-1.)/6.;

    /* Identify location to interpolate in log(TR) */
    iTR = (long)floor((logTR - log(TR_MIN))/rate_table->DlogTR);
    if (iTR < 1) iTR = 1;
    if (iTR > NTR-3) iTR = NTR-3;
    frac2 = (logTR - log(TR_MIN))/rate_table->DlogTR - iTR;
    coeff2[0] = frac2*(frac2-1.)*(2.-frac2)/6.;
    coeff2[1] = (1.+frac2)*(1.-frac2)*(2.-frac2)/2.;
    coeff2[2] = (1.+frac2)*frac2*(2.-frac2)/2.;
    coeff2[3] = (1.+frac2)*frac2*(frac2-1.)/6.;


    for (l = 0; l <= 1; l++) {
    /* effective recombination coefficient to each level */
       for (k = 0; k < 4; k++) {
           temp[k] = rate_table->logAlpha_tab[l][iTM-1+k][iTR-1]*coeff2[0]
                   + rate_table->logAlpha_tab[l][iTM-1+k][iTR]*coeff2[1]
                   + rate_table->logAlpha_tab[l][iTM-1+k][iTR+1]*coeff2[2]
                   + rate_table->logAlpha_tab[l][iTM-1+k][iTR+2]*coeff2[3];
       }

      Alpha[l] = exp(temp[0]*coeff1[0]+temp[1]*coeff1[1]
                    +temp[2]*coeff1[2]+temp[3]*coeff1[3]);

    /* Alpha evaluated at Tm = Tr, for use in detailed balance */
    Beta[l] = exp(rate_table->logAlpha_tab[l][NTM-1][iTR-1]*coeff2[0]
                 +rate_table->logAlpha_tab[l][NTM-1][iTR]*coeff2[1]
                 +rate_table->logAlpha_tab[l][NTM-1][iTR+1]*coeff2[2]
                 +rate_table->logAlpha_tab[l][NTM-1][iTR+2]*coeff2[3]);
    }

    /* Beta obtained by detailed balance */
    /* factor = pow(2.0 * M_PI * mue *TR / hPc / hPc, 1.5)) * exp(-0.25*EI/TR) */
    factor = 3.016103031869581e21 *  TR*sqrt(TR) * exp(-3.399571517984581/TR);
    Beta[0] *= factor;              /* 2s */
    Beta[1] *= factor/3.;           /* 2p */

    /* Effective 2p->2s rate */
    *R2p2s = exp(rate_table->logR2p2s_tab[iTR-1]*coeff2[0]
                +rate_table->logR2p2s_tab[iTR]*coeff2[1]
                +rate_table->logR2p2s_tab[iTR+1]*coeff2[2]
                +rate_table->logR2p2s_tab[iTR+2]*coeff2[3]);

}

/************************************************************************************************
Solves for the populations of the 2s and 2p states in steady-state, and returns dxe/dt.
Uses standard rate for 2s-->1s decay and Sobolev for Lyman alpha (no feedback)
Inputs: xe, nH in cm^{-3}, H in s^{-1}, TM, TR in eV. Output: dxe/dlna
************************************************************************************************/

double rec_HMLA_dxedlna(double xe, double nH, double Hubble, double TM, double TR, double energy_rate, HRATEEFF *rate_table){

   double Alpha[2];
   double Beta[2];
   double R2p2s;
   double matrix[2][2];
   double RHS[2];
   double det, RLya;
   double x2[2];
   double x1s_db;
   double C_2p;
   double chi_ion_H;

   interpolate_rates(Alpha, Beta, &R2p2s, TR, TM / TR, rate_table);

   x1s_db = (1.-xe)*exp(-E21/TR);

   /******** 2s *********/
   matrix[0][0] = Beta[0] + 3.*R2p2s + L2s1s;
   matrix[0][1] = -R2p2s;
   RHS[0]       = xe*xe*nH *Alpha[0] + L2s1s *x1s_db;

   /******** 2p *********/
   RLya = 4.662899067555897e15 *Hubble /nH/(1.-xe);   /*8 PI H/(3 nH x1s lambda_Lya^3) */

   matrix[1][1] = Beta[1] + R2p2s + RLya;
   matrix[1][0] = -3.*R2p2s;
   RHS[1]       = xe*xe*nH *Alpha[1] + 3.*RLya *x1s_db;

   /**** invert the 2 by 2 system matrix_ij * xj = RHSi *****/
   det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
   x2[0] = (matrix[1][1] * RHS[0] - matrix[0][1] * RHS[1])/det;
   x2[1] = (matrix[0][0] * RHS[1] - matrix[1][0] * RHS[0])/det;

   C_2p=(RLya+R2p2s*L2s1s/matrix[0][0])/(matrix[1][1]-R2p2s*3.*R2p2s/matrix[0][0]);

   //chi_ion_H = (1.-xe)/3.; // old approximation from Chen and Kamionkowski
    if (xe < 1.)
      chi_ion_H = 0.369202*pow(1.-pow(xe,463929),1.70237); // coefficient as revised by Galli et al. 2013 (in fact it is a fit by Vivian Poulin of columns 1 and 4 in Table V of Galli et al. 2013)
    else
      chi_ion_H = 0.;

   return  (x1s_db*(L2s1s + 3.*RLya) -x2[0]*L2s1s -x2[1]*RLya)/Hubble
     +chi_ion_H/nH*energy_rate*(1./EI+(1.-C_2p)/E21)/Hubble;

}

/*********************************************************************************************
Read two-photon rates from table
**********************************************************************************************/

void read_twog_params(TWO_PHOTON_PARAMS *twog){

   FILE *fA;
   unsigned b;
   double L2s1s_current;
   int fscanf_result;

   fA = fopen(TWOG_FILE, "r");

   for (b = 0; b < NVIRT; b++) {
      fscanf_result = 0;
      fscanf_result += fscanf(fA, "%le", &(twog->Eb_tab[b]));
      fscanf_result += fscanf(fA, "%le", &(twog->A1s_tab[b]));
      fscanf_result += fscanf(fA, "%le", &(twog->A2s_tab[b]));
      fscanf_result += fscanf(fA, "%le", &(twog->A3s3d_tab[b]));
      fscanf_result += fscanf(fA, "%le", &(twog->A4s4d_tab[b]));
      if(fscanf_result!=5){printf("Hyrec Warning :: Could not read Two Photon table (two_photon_tables.dat)");}
   }
   fclose(fA);

   /* Normalize 2s--1s differential decay rate to L2s1s (can be set by user in hydrogen.h) */
   L2s1s_current = 0.;
   for (b = 0; b < NSUBLYA; b++) L2s1s_current += twog->A2s_tab[b];
   for (b = 0; b < NSUBLYA; b++) twog->A2s_tab[b] *= L2s1s/L2s1s_current;


  /* Switches for the various effects considered in Hirata (2008) and diffusion:
      Effect A: correct 2s-->1s rate, with stimulated decays and absorptions of non-thermal photons
      Effect B: Sub-Lyman-alpha two-photon decays
      Effect C: Super-Lyman-alpha two-photon decays
      Effect D: Raman scattering */

   #if (EFFECT_A == 0)
     for (b = 0; b < NSUBLYA; b++) twog->A2s_tab[b] = 0;
   #endif
   #if (EFFECT_B == 0)
     for (b = 0; b < NSUBLYA; b++) twog->A3s3d_tab[b] = twog->A4s4d_tab[b] = 0;
   #endif
   #if (EFFECT_C == 0)
      for (b = NSUBLYA; b < NVIRT; b++) twog->A3s3d_tab[b] = twog->A4s4d_tab[b] = 0;
   #endif
   #if (EFFECT_D == 0)
      for (b = NSUBLYA; b < NVIRT; b++) twog->A2s_tab[b] = 0;
      for (b = NSUBLYB; b < NVIRT; b++) twog->A3s3d_tab[b] = 0;
   #endif
   #if (DIFFUSION == 0)
      for (b = 0; b < NVIRT; b++) twog->A1s_tab[b] = 0;
   #endif

}

/******************************************************************************************************
Square and cube, used often
******************************************************************************************************/

double square(double x) {
   return x*x;
}

double cube(double x) {
   return x*x*x;
}

/********************************************************************************************************
Compute the A_{b,b+/-1} "Einstein A-"coefficients between virtual states, due to diffusion
(In the notation of the paper, A_{b,b\pm1} = R_{b,b\pm1})
Aup[b] = A_{b, b+1}    Adn[b] = A_{b, b-1}
********************************************************************************************************/

void populate_Diffusion(double *Aup, double *Adn, double *A2p_up, double *A2p_dn,
                        double TM, double Eb_tab[NVIRT], double A1s_tab[NVIRT]) {

    unsigned b;
    double DE2;

    DE2 = E21*E21*2.*TM/mH;

   /****** RED WING ******/
    b = NSUBLYA - NDIFF/2;
    Aup[b] = DE2 / square(Eb_tab[b+1] - Eb_tab[b]) * A1s_tab[b];    /* A{0,1}. Assume A{0,-1} = 0 */

    for (b = NSUBLYA - NDIFF/2 + 1; b < NSUBLYA-1; b++) {
       Adn[b] = exp((Eb_tab[b] - Eb_tab[b-1])/TM) * Aup[b-1];                      /* Detailed balance */
       Aup[b] = (DE2 * A1s_tab[b] - square(Eb_tab[b] - Eb_tab[b-1]) * Adn[b])
                       /square(Eb_tab[b+1] - Eb_tab[b]);        /* Aup[b] , Adn[b] must add up to correct diffusion rate */
    }
    /* Last bin below Lyman alpha */
    b = NSUBLYA - 1;
    Adn[b] = exp((Eb_tab[b] - Eb_tab[b-1])/TM) * Aup[b-1];

    Aup[b] = (DE2 * A1s_tab[b] - square(Eb_tab[b] - Eb_tab[b-1]) * Adn[b])
                    /square(E21 - Eb_tab[b]);
    *A2p_dn = exp((E21 - Eb_tab[b])/TM)/3.* Aup[b];   /* 2p -> NSUBLYA-1 rate obtained by detailed balance */


    /****** BLUE WING ******/
    b = NSUBLYA + NDIFF/2 - 1;
    Adn[b] = DE2 / square(Eb_tab[b] - Eb_tab[b-1]) * A1s_tab[b];

    for (b = NSUBLYA + NDIFF/2 - 2; b > NSUBLYA; b--) {
       Aup[b] = exp((Eb_tab[b] - Eb_tab[b+1])/TM) * Adn[b+1];
       Adn[b] = (DE2 * A1s_tab[b] - square(Eb_tab[b+1] - Eb_tab[b]) * Aup[b])
                   /square(Eb_tab[b] - Eb_tab[b-1]);
    }
    /* First bin above Lyman alpha */
    b = NSUBLYA;
    Aup[b] = exp((Eb_tab[b] - Eb_tab[b+1])/TM) * Adn[b+1];

    Adn[b] = (DE2 * A1s_tab[b] - square(Eb_tab[b+1] - Eb_tab[b]) * Aup[b])
                 /square(Eb_tab[b] - E21);
    *A2p_up = exp((E21 - Eb_tab[b])/TM)/3. * Adn[b];   /* 2p -> NSUBLYA rate obtained by detailed balance */


}

/*********************************************************************************************************
 Populate the real-real, real-virtual, virtual-real and virtual-virtual T-matrices,
 as well as the source vectors sr, sv,
WITH DIFFUSION. Tvv[0][b] is the diagonal element Tbb, Tvv[1][b] = T{b,b-1} and Tvv[2][b] = T{b,b+1}
Also, computes and stores the optical depths Delta tau_b for future use
**********************************************************************************************************/

void populateTS_2photon(double Trr[2][2], double *Trv[2], double *Tvr[2], double *Tvv[3],
                        double sr[2], double sv[NVIRT], double Dtau[NVIRT],
                        double xe, double TM, double TR, double nH, double H, HRATEEFF *rate_table,
                        TWO_PHOTON_PARAMS *twog, double fplus[NVIRT], double fplus_Ly[],
                        double Alpha[2], double Beta[2], double z) {

   double R2p2s;
   unsigned b;
   double  RLya, Gammab, Pib, dbfact;

   double *Aup, *Adn;
   double A2p_up, A2p_dn;

   Aup = create_1D_array(NVIRT);
   Adn = create_1D_array(NVIRT);

   RLya = 4.662899067555897e15 *H /nH/(1.-xe);   /*8 PI H/(3 nH x1s lambda_Lya^3) */

   interpolate_rates(Alpha, Beta, &R2p2s, TR, TM / TR, rate_table);

   /****** 2s row and column ******/

   Trr[0][0] = Beta[0] + 3.*R2p2s
             + 3.* RLya * (1.664786871919931 *exp(-E32/TR)     /* Ly-beta escape */
	                          + 1.953125 *exp(-E42/TR));   /* Ly-gamma escape */

   Trr[0][1] = -R2p2s;
   sr[0]     = nH * Alpha[0] * xe*xe
             + 3.* RLya * (1.-xe) * (1.664786871919931 *fplus_Ly[1]
                                    + 1.953125 *exp(-E41/TR));

   #if (EFFECT_A == 0)
       /* Standard treatment of 2s-->1s two-photon decays */
       Trr[0][0] += L2s1s;
       sr[0]     += L2s1s * (1.-xe) * exp(-E21/TR);
   #endif


   /****** 2p row and column ******/

   Trr[1][1] = Beta[1] + R2p2s + RLya;
   Trr[1][0] = -3.*R2p2s;
   sr[1]     = nH * Alpha[1] * xe*xe + 3.*RLya * (1.-xe) * fplus_Ly[0];


   /***** Two-photon transitions: populating Trv, Tvr and updating Trr ******/

   for (b = 0; b < NVIRT; b++) {
       dbfact = exp((twog->Eb_tab[b] - E21)/TR);

       Trr[0][0] -= Tvr[0][b] = -twog->A2s_tab[b]/fabs(exp((twog->Eb_tab[b] - E21)/TR)-1.);
       Trv[0][b]  = Tvr[0][b] *dbfact;

       Trr[1][1] -= Tvr[1][b] = -exp(-E32/TR)/3. * twog->A3s3d_tab[b]/fabs(exp((twog->Eb_tab[b] - E31)/TR)-1.)
                                -exp(-E42/TR)/3. * twog->A4s4d_tab[b]/fabs(exp((twog->Eb_tab[b] - E41)/TR)-1.);
       Trv[1][b] = Tvr[1][b] *3.*dbfact;
   }

    /****** Tvv and sv. Accounting for DIFFUSION ******/

    populate_Diffusion(Aup, Adn, &A2p_up, &A2p_dn, TM, twog->Eb_tab, twog->A1s_tab);

    /* Updating Tvr, Trv, Trr for diffusion between line center ("2p state") and two neighboring bins */

    Trr[1][1] += (A2p_dn + A2p_up);

    for (b = 0; b < NVIRT; b++) {

      Gammab = -(Trv[0][b] + Trv[1][b]) + Aup[b] + Adn[b];    /* Inverse lifetime of virtual state b */

      /*** Diffusion region ***/
      if (  (b >= NSUBLYA - NDIFF/2 && b < NSUBLYA - 1)
          ||(b > NSUBLYA && b < NSUBLYA + NDIFF/2)) {
	  Tvv[1][b] = -Aup[b-1];
          Tvv[2][b] = -Adn[b+1];
      }
      /* Bins neighboring Lyman alpha. */
      if (b == NSUBLYA-1)  {
          Tvv[1][b] = -Aup[b-1];
          Tvv[2][b] = 0.;
          Tvr[1][b] -= A2p_dn;
          Trv[1][b] -= Aup[b];
      }
      if (b == NSUBLYA)  {
          Tvv[1][b] = 0.;
          Tvv[2][b] = -Adn[b+1];
          Tvr[1][b] -= A2p_up;
          Trv[1][b] -= Adn[b];
      }
      /*********************/

      Dtau[b] = Gammab * (1.-xe) * cube(hPc/twog->Eb_tab[b]) * nH /8. /M_PI /H;

       if (Dtau[b] > 1e-30) {
          Pib = (1.-exp(-Dtau[b]))/Dtau[b];
          Tvv[0][b] = Gammab/(1.-Pib);
          sv[b]  = Tvv[0][b] * (1.-xe) * fplus[b] * Pib;
       }
       else {  /* Nearly vanishing optical depth: free streaming */
          Tvv[0][b] = 1.;
          Trv[0][b] = Trv[1][b] = Tvr[0][b] = Tvr[1][b] = 0;
          sv[b] = (1.-xe) * fplus[b];
       }
   }

    free(Aup);
    free(Adn);
}

/*********************************************************************
Solves the linear system T*X = B, where T is a DIAGONALLY DOMINANT
tridiagonal matrix, and X and B are both one-column vectors of N elements.
diag[i] = T_{ii}, updiag[i] = T_{i,i+1}, dndiag[i] = T_{i,i-1}
IMPORTANT: This is NOT THE MOST GENERAL ALGORITHM. Only adapted for the
case we will consider, i.e. |T_{ii}| > |T_{i,i+1}| + |T_{i,i-1}|.
**********************************************************************/

void solveTXeqB(double *diag, double *updiag, double *dndiag,
                double *X, double *B, unsigned N){
     int i;
     double denom;
     double *alpha = create_1D_array(N);
     double *gamma = create_1D_array(N);  /* X[i] = gamma[i] - alpha[i] * X[i+1] */

     alpha[0] = updiag[0] / diag[0];
     gamma[0] = B[0] / diag[0];

     for (i = 1; i < N; i++) {
         denom = diag[i] - dndiag[i] * alpha[i-1];
         alpha[i] = updiag[i] / denom;
         gamma[i] = (B[i] - dndiag[i] * gamma[i-1]) / denom;
     }

     X[N-1] = gamma[N-1];
     for (i = N-2; i >= 0; i--) X[i] = gamma[i] - alpha[i] * X[i+1];

     free(alpha);
     free(gamma);
}

/**************************************************************************************************************
Solves for the populations of the real (2s, 2p) and virtual states
***************************************************************************************************************/

void solve_real_virt(double xr[2], double xv[NVIRT], double Trr[2][2], double *Trv[2], double *Tvr[2],
                     double *Tvv[3], double sr[2], double sv[NVIRT]){

   double *Tvv_inv_Tvr[2];
   double *Tvv_inv_sv;
   double Trr_new[2][2];
   double sr_new[2];
   unsigned i, j, b;
   double det;

   unsigned NSUBDIFF;
   NSUBDIFF = NSUBLYA - NDIFF/2;     /* lowest bin of the diffusion region */

   /** Allocate memory **/
   for (i = 0; i < 2; i++) Tvv_inv_Tvr[i] = create_1D_array(NVIRT);
   Tvv_inv_sv = create_1D_array(NVIRT);


   /*** Computing Tvv^{-1}.Tvr ***/

   for (i = 0; i < 2; i++) {
      for (b = 0; b < NSUBDIFF; b++)               Tvv_inv_Tvr[i][b] = Tvr[i][b]/Tvv[0][b];
      for (b = NSUBLYA + NDIFF/2; b < NVIRT; b++)  Tvv_inv_Tvr[i][b] = Tvr[i][b]/Tvv[0][b];
      solveTXeqB(Tvv[0]+NSUBDIFF, Tvv[2]+NSUBDIFF, Tvv[1]+NSUBDIFF, Tvv_inv_Tvr[i]+NSUBDIFF, Tvr[i]+NSUBDIFF, NDIFF);
   }

   /*** Trr_new = Trr - Trv.Tvv^{-1}.Tvr ***/
   for (i = 0; i < 2; i++) for (j = 0; j < 2; j++) {
          Trr_new[i][j] = Trr[i][j];
          for (b = 0; b < NVIRT; b++) Trr_new[i][j] -= Trv[i][b]*Tvv_inv_Tvr[j][b];
   }

   /*** Computing Tvv^{-1}.sv ***/
   for (b = 0; b < NSUBDIFF; b++)               Tvv_inv_sv[b] = sv[b]/Tvv[0][b];
   for (b = NSUBLYA + NDIFF/2; b < NVIRT; b++)  Tvv_inv_sv[b] = sv[b]/Tvv[0][b];
   solveTXeqB(Tvv[0]+NSUBDIFF, Tvv[2]+NSUBDIFF, Tvv[1]+NSUBDIFF, Tvv_inv_sv+NSUBDIFF, sv+NSUBDIFF, NDIFF);

   /*** sr_new = sr - Trv.Tvv^{-1}sv ***/
   for (i = 0; i < 2; i++) {
      sr_new[i] = sr[i];
      for (b = 0; b < NVIRT; b++) sr_new[i] -= Trv[i][b]*Tvv_inv_sv[b];
   }

   /*** Solve 2 by 2 system Trr_new.xr = sr_new ***/
   det = Trr_new[0][0] * Trr_new[1][1] - Trr_new[0][1] * Trr_new[1][0];
   xr[0] = (Trr_new[1][1] * sr_new[0] - Trr_new[0][1] * sr_new[1])/det;
   xr[1] = (Trr_new[0][0] * sr_new[1] - Trr_new[1][0] * sr_new[0])/det;

   /*** xv = Tvv^{-1}(sv - Tvr.xr) ***/
   for (b = 0; b < NVIRT; b++) xv[b] = Tvv_inv_sv[b] - Tvv_inv_Tvr[0][b]*xr[0] - Tvv_inv_Tvr[1][b]*xr[1];


   /** Free memory **/
   for (i = 0; i < 2; i++) free(Tvv_inv_Tvr[i]);
   free(Tvv_inv_sv);


}

/*************************************************************************************************************
Obtain fplus at each bin, given the history of fminus (simple free-streaming). iz is the current time step.
fminus[0..iz-1] is known.
Assume the Lyman lines are optically thick
logfminus_hist is a NVIRT by nz array of previous log(f(nu_b - epsilon)(z))
logfminus_Ly_hist is a 3 by nz array of previous log(f(nu - epsilon)) redward of Ly alpha, beta and gamma lines
*************************************************************************************************************/

void fplus_from_fminus(double fplus[NVIRT], double fplus_Ly[], double **logfminus_hist, double *logfminus_Ly_hist[],
                       double TR, double zstart, double dlna, unsigned iz, double z, double Eb_tab[NVIRT])
{
   unsigned b;
   double ainv, lna_start, zp1;

   zp1 = 1.+z;
   lna_start = -log(1.+zstart);

   /*** Bins below Lyman alpha ***/
   for (b = 0; b < NSUBLYA-1; b++) {
      ainv = zp1*Eb_tab[b+1]/Eb_tab[b];
      fplus[b] = exp(rec_interp1d(lna_start, dlna, logfminus_hist[b+1], iz, -log(ainv)));
   }

   /*** highest bin below Ly-alpha: feedback from optically thick Ly-alpha ***/
   b = NSUBLYA-1;
   ainv = zp1*E21/Eb_tab[b];
   fplus[b] = exp(rec_interp1d(lna_start, dlna, logfminus_Ly_hist[0], iz, -log(ainv)));

   /*** incoming photon occupation number at Lyman alpha ***/
   b = NSUBLYA;     /* next highest bin */
   ainv = zp1*Eb_tab[b]/E21;
   fplus_Ly[0] = exp(rec_interp1d(lna_start, dlna, logfminus_hist[b], iz, -log(ainv)));

   /*** Bins between Lyman alpha and beta ***/
   for (b = NSUBLYA; b < NSUBLYB-1; b++) {
     ainv = zp1*Eb_tab[b+1]/Eb_tab[b];
     fplus[b] = exp(rec_interp1d(lna_start, dlna, logfminus_hist[b+1], iz, -log(ainv)));
   }

   /*** highest bin below Ly-beta: feedback from Ly-beta ***/
   b = NSUBLYB-1;
   ainv = zp1*E31/Eb_tab[b];
   fplus[b] = exp(rec_interp1d(lna_start, dlna, logfminus_Ly_hist[1], iz, -log(ainv)));

   /*** incoming photon occupation number at Lyman beta ***/
   b = NSUBLYB;     /* next highest bin */
   ainv = zp1*Eb_tab[b]/E31;
   fplus_Ly[1] = exp(rec_interp1d(lna_start, dlna, logfminus_hist[b], iz, -log(ainv)));


   /*** Bins between Lyman beta and gamma ***/
   for (b = NSUBLYB; b < NVIRT-1; b++) {
     ainv = zp1*Eb_tab[b+1]/Eb_tab[b];
     fplus[b] = exp(rec_interp1d(lna_start, dlna, logfminus_hist[b+1], iz, -log(ainv)));
   }

   /*** highest energy bin: feedback from Ly-gamma ***/
   b = NVIRT-1;
   ainv = zp1*E41/Eb_tab[b];
   fplus[b] = exp(rec_interp1d(lna_start, dlna, logfminus_Ly_hist[2], iz, -log(ainv)));

}

/******************************************************************************************************************
d(x_e)/dt when including two-photon processes
This is independent of whether diffusion is present or not.
fminus[0..iz-1] is known. Update fminus[iz]
******************************************************************************************************************/

double rec_HMLA_2photon_dxedlna(double xe, double nH, double H, double TM, double TR,
                                HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog,
                                double zstart, double dlna, double **logfminus_hist, double *logfminus_Ly_hist[],
                                unsigned iz, double z, double energy_rate){

   double xr[2];
   double xv[NVIRT];
   double xedot, Pib, feq;
   double fplus[NVIRT], fplus_Ly[3];
   unsigned b, i;

   double Trr[2][2];
   double matrix[2][2];
   double *Trv[2];
   double *Tvr[2];
   double *Tvv[3];
   double sr[2];
   double sv[NVIRT];
   double Dtau[NVIRT];
   double Alpha[2], Beta[2];

   double RLya;
   double R2p2s;
   double C_2p;

   double chi_ion_H;

   for (i = 0; i < 2; i++) Trv[i] = create_1D_array(NVIRT);
   for (i = 0; i < 2; i++) Tvr[i] = create_1D_array(NVIRT);
   for (i = 0; i < 3; i++) Tvv[i] = create_1D_array(NVIRT);

   /* Redshift photon occupation number from previous times and higher energy bins */
   fplus_from_fminus(fplus, fplus_Ly, logfminus_hist, logfminus_Ly_hist, TR,
                     zstart, dlna, iz, z, twog->Eb_tab);

   /* Compute real-real, real-virtual and virtual-virtual transition rates */
   populateTS_2photon(Trr, Trv, Tvr, Tvv, sr, sv, Dtau, xe, TM, TR, nH, H, rate_table,
                      twog, fplus, fplus_Ly, Alpha, Beta, z);

   /* Solve for the population of the real and virtual states */
   solve_real_virt(xr, xv, Trr, Trv, Tvr, Tvv, sr, sv);


   /*************************************************************/

   /* Dark matter annihilation*/
   RLya = 4.662899067555897e15 *H /nH/(1.-xe);   /*8 PI H/(3 nH x1s lambda_Lya^3) */

   interpolate_rates(Alpha, Beta, &R2p2s, TR, TM / TR, rate_table);

   matrix[0][0] = Beta[0] + 3.*R2p2s + L2s1s;
   matrix[1][1] = Beta[1] + R2p2s + RLya;

   C_2p=(RLya+R2p2s*L2s1s/matrix[0][0])/(matrix[1][1]-R2p2s*3.*R2p2s/matrix[0][0]);


   /*************************************************************/

   //chi_ion_H = (1.-xe)/3.; // old approximation from Chen and Kamionkowski
   if (xe < 1.)
     chi_ion_H = 0.369202*pow(1.-pow(xe,0.463929),1.70237); // coefficient as revised by Galli et al. 2013 (in fact it is a fit by Vivian Poulin of columns 1 and 4 in Table V of Galli et al. 2013)
   else
     chi_ion_H = 0.;

   /* Obtain xe_dot */
   xedot = -nH*xe*xe*(Alpha[0]+Alpha[1]) + xr[0]*Beta[0] + xr[1]*Beta[1]
	+chi_ion_H/nH*energy_rate*(1./EI+(1.-C_2p)/E21);


   /* Update fminuses */

   for (b = 0; b < NVIRT; b++) {
     if (Dtau[b] != 0) {
         Pib = (1.-exp(-Dtau[b]))/Dtau[b];
         feq  = -xr[0]*Tvr[0][b] - xr[1]*Tvr[1][b];
         feq -= (b == 0       ?  xv[1]*Tvv[2][0]:
                   b == NVIRT-1 ?  xv[NVIRT-2]*Tvv[1][NVIRT-1]:
                   xv[b+1]*Tvv[2][b] + xv[b-1]*Tvv[1][b]);
         feq /= (1.-xe)*(1.-Pib)*Tvv[0][b];

         logfminus_hist[b][iz] = log(fplus[b] + (feq - fplus[b])*(1.-exp(-Dtau[b])));
     }
     else logfminus_hist[b][iz] = log(fplus[b]);
   }

   /* log(xnp/(3 x1s)) , n = 2, 3, 4, assuming 3p and 4p in Boltzmann eq. with 2s */
   logfminus_Ly_hist[0][iz] = log(xr[1]/3./(1.-xe));
   logfminus_Ly_hist[1][iz] = log(xr[0]/(1.-xe)) - E32/TR;
   logfminus_Ly_hist[2][iz] = log(xr[0]/(1.-xe)) - E42/TR;


   for (i = 0; i < 2; i++) free(Trv[i]);
   for (i = 0; i < 2; i++) free(Tvr[i]);
   for (i = 0; i < 3; i++) free(Tvv[i]);

   return xedot/H;
}

/******************************************************************************************
Post-Saha expansion using the full derivative with two-photon processes and diffusion
ADDED JANUARY 2011
*******************************************************************************************/

double xe_PostSahaH(double nH, double H, double T, HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog,
                           double zstart, double dlna, double **logfminus_hist, double *logfminus_Ly_hist[],
		    unsigned iz, double z, double *Dxe, int model, double energy_rate){

    double s, xeSaha, dxeSaha_dlna, Ddxedlna_Dxe;

    s = 3.016103031869581e21 *T*sqrt(T) *exp(-EI/T)/nH;                     /* xe^2/(1-xe) in Saha eq.*/
    xeSaha = 2./(1.+sqrt(1.+4./s));                                         /* Saha equilibrium */
    dxeSaha_dlna = -(EI/T - 1.5)/(2.*xeSaha + s)*xeSaha*xeSaha;             /* Analytic derivative of above expression */

    if (model == 0) {         /* Peebles model */
      Ddxedlna_Dxe = (rec_HPeebles_dxedlna(xeSaha+0.01*(1.-xeSaha), nH, H, T, T, energy_rate)
		      - rec_HPeebles_dxedlna(xeSaha-0.01*(1.-xeSaha), nH, H, T, T, energy_rate))/0.02/(1.-xeSaha);
    }
    else if (model == 1)  {   /* "Recfast" model */
        Ddxedlna_Dxe = (rec_HRecFast_dxedlna(xeSaha+0.01*(1.-xeSaha), nH, H, T, T, energy_rate)
		      - rec_HRecFast_dxedlna(xeSaha-0.01*(1.-xeSaha), nH, H, T, T, energy_rate))/0.02/(1.-xeSaha);
    }
    else if (model == 2) { /* EMLA model with 2s and 2p decays only, no radiative transfer */
        Ddxedlna_Dxe = (rec_HMLA_dxedlna(xeSaha+0.01*(1.-xeSaha), nH, H, T, T, energy_rate, rate_table)
		      - rec_HMLA_dxedlna(xeSaha-0.01*(1.-xeSaha), nH, H, T, T, energy_rate, rate_table))/0.02/(1.-xeSaha);
    }
    else {    /* Default mode, with radiative transfer */
        Ddxedlna_Dxe = (rec_HMLA_2photon_dxedlna(xeSaha+0.01*(1.-xeSaha),nH, H, T, T,
                         rate_table, twog, zstart, dlna, logfminus_hist, logfminus_Ly_hist, iz, z, energy_rate)
                      - rec_HMLA_2photon_dxedlna(xeSaha-0.01*(1.-xeSaha), nH, H, T, T,
	                 rate_table, twog, zstart, dlna, logfminus_hist, logfminus_Ly_hist, iz, z, energy_rate))/0.02/(1.-xeSaha);
    }

    *Dxe = dxeSaha_dlna/Ddxedlna_Dxe;

     /* Compute derivative again just so the fminuses are properly updated */
    //dxedlna = rec_HMLA_2photon_dxedlna(xeSaha + *Dxe, nH, H, T, T,rate_table, twog, zstart, dlna, logfminus_hist, logfminus_Ly_hist, iz, z, energy_rate);

    return xeSaha + *Dxe;

}

/*****************************************************************************************************
Update fminus(z) given the history of previous fminus's, assuming:
- thermal radiation field and Boltzmann equilibrium for the excited states if func_select = 0
- excited states are in Saha eq. with xe, and free-streaming for the radiation field between
  optically thick Lyman lines if func_select = 1
*****************************************************************************************************/

void update_fminus_Saha(double **logfminus_hist, double *logfminus_Ly_hist[],
                        double xe, double TR, double nH, TWO_PHOTON_PARAMS *twog,
			double zstart, double dlna, unsigned iz, double z, int func_select){

    double fplus[NVIRT];
    double fplus_Ly[2];
    unsigned b;
    double common_fact;

    if (func_select == 0) {
       for (b = 0; b < NVIRT; b++) logfminus_hist[b][iz] = -twog->Eb_tab[b]/TR;
       logfminus_Ly_hist[0][iz] = -E21/TR;
       logfminus_Ly_hist[1][iz] = -E31/TR;
       logfminus_Ly_hist[2][iz] = -E41/TR;
    }
    else {
        fplus_from_fminus(fplus, fplus_Ly, logfminus_hist, logfminus_Ly_hist, TR,
                          zstart, dlna, iz, z, twog->Eb_tab);

         for (b = 0; b < NVIRT; b++) logfminus_hist[b][iz] = log(fplus[b]);  /* free streaming */

         common_fact = log(nH*xe*xe/(1.-xe)/(3.016103031869581e21 *TR*sqrt(TR)));
         logfminus_Ly_hist[0][iz] = common_fact + EI/4./TR;         /* Saha equilibrium with the continuum */
         logfminus_Ly_hist[1][iz] = common_fact + EI/9./TR;
         logfminus_Ly_hist[2][iz] = common_fact + EI/16./TR;
    }
}

/***********************************************************************************************************/
