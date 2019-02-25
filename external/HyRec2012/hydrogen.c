/*************************************************************************************************/
/*                 HYREC: Hydrogen and Helium Recombination Code                                 */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                              */
/*                                                                                               */
/*         hydrogen.c: all functions related to Hydrogen recombination                           */
/*                                                                                               */
/*         Units used: cgs + eV (all temperatures in eV)                                         */
/*                                                                                               */
/*         Version: May 2012                                                                     */ 
/*                                                                                               */
/*         Revision history:                                                                     */
/*            - Written November 2010                                                            */
/*            - January 2011: - changed the post-Saha expansion to use the full derivative       */
/*                       (including two-photon processes and diffusion) rather than Peebles'ODE  */
/*	                      - post-Saha expansion can now pass the difference from Saha value  */
/*                                  to external routines                                         */              
/*                            - differential 2s--1s rate is now systematically normalized to     */
/*                                 total 2s--1s rate that can be set by user in hydrogen.h       */
/*            - May 2012: - Now solve for the photon distortion instead of absolute value        */
/*                               of radiation field (less prone to numerical errors)             */
/*                          - Improved the post-Saha expansion to properly account for           */
/*                               non-thermal distortions                                         */
/*                          - Added explicit dependence on fine-structure constant               */
/*                           (fsR = alpha/alpha0) and electron mass (meR = me/me0 ~ mue/mue0)    */
/*                          - Included dependence on xe and xHII in case HeII still present      */
/*************************************************************************************************/ 

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "hyrectools.h"
#include "hydrogen.h"


/***********************************************************************************************************
Some constants appropriately rescaled for different values of the fine-structure constant and electron mass 
***********************************************************************************************************/

double SAHA_FACT(double fsR, double meR) {
  return 3.016103031869581e21 * cube(fsR*meR);     /* (2 pi mu_e * EI/EI0)^(3/2)/h^3 in eV^(-3/2) cm^(-3), used for Saha equilibrium and detailed balance */
}

double LYA_FACT(double fsR, double meR) {
  return 4.662899067555897e15 * cube(fsR*fsR*meR); /* 8pi/(3 lambda_Lya^3) in cm^(-3), used as prefactor for Lyman-alpha Sobolev escape probability */
}

double L2s_rescaled(double fsR, double meR) {     /* Rescaled two-photon decay rate */ 
   return L2s1s * square(fsR*fsR*fsR*fsR) * meR;
}

/*****************************************************************************************************
Temperature rescaling given fine-structure constant and electron mass, so we can use today's EI value
*****************************************************************************************************/

void rescale_T(double *T, double fsR, double meR) {
   *T /= fsR*fsR*meR;
}

/**************************************************************************************************
Case-B recombination coefficient and photoionization rate, fit of Pequignot et al 1991, in cm^3 s^{-1}
INPUT TEMPERATURE ASSUMED TO BE ALREADY RESCALED FOR VALUES OF alpha_fs and me
***************************************************************************************************/

double alphaB_PPB(double TM, double fsR, double meR) {
   double t4;   
  
   t4 = TM/kBoltz/1e4;
   return square(fsR/meR) * 4.309e-13*pow(t4,-0.6166)/(1.+ 0.6703*pow(t4,0.5300));
}

/**************************************************************************************************
Effective three-level atom model with adjustable fudge factor.
Fudge = 1 is Peebles' model. Fudge = 1.14 is similar to RecFast (Seager et al. 1999, 2000).
Changes (May 2012):
- Correction: detailed balance implies that the photoionization rate is proportional to alpha_B(Tr) rather than alpha_B(Tm)
- Analytically subtract nearly-cancelling terms at high-z
- Explicit dependence on alpha_fs and m_e now accounted for
- Added explicit dependence on xHII, which is not necessarily equal to xe if Helium has not entirely recombined
***************************************************************************************************/

double rec_TLA_dxHIIdlna(double xe, double xHII, double nH, double H, double TM, double TR, double Fudge, 
                         double fsR, double meR) {

   double RLya, alphaB_TM, alphaB_TR, four_betaB, C, s, Dxe2, DalphaB;
  
   rescale_T(&TM, fsR, meR);
   rescale_T(&TR, fsR, meR);    

   RLya       = LYA_FACT(fsR, meR) * H / nH / (1.-xHII);
   alphaB_TM  = Fudge * alphaB_PPB(TM, fsR, meR);
   alphaB_TR  = Fudge * alphaB_PPB(TR, fsR, meR);
   
   four_betaB = SAHA_FACT(fsR, meR) *TR*sqrt(TR) *exp(-0.25*EI/TR) * alphaB_TR;
   C          = (3.*RLya + L2s_rescaled(fsR, meR))/(3.*RLya + L2s_rescaled(fsR, meR) + four_betaB);    /* Peebles' C-factor */ 
 
   s          = SAHA_FACT(fsR, meR) *TR*sqrt(TR) *exp(-EI/TR)/nH;
   Dxe2       = xe*xHII - s*(1.-xHII);    /* xe xp - xe xp[Saha eq with 1s] -- gives more compact expressions */
   DalphaB    = alphaB_TM - alphaB_TR; 

   return -nH*(s*(1.-xHII)*DalphaB + Dxe2*alphaB_TM)*C/H; 

}

/**********************************************************************************************
Store tabulated temperatures and read effective MLA rates from files.
**********************************************************************************************/

void read_rates(HRATEEFF *rate_table){

   FILE *fA = fopen(ALPHA_FILE, "r");
   FILE *fR = fopen(RR_FILE, "r"); 

   unsigned i, j, l;

   maketab(log(TR_MIN), log(TR_MAX), NTR, rate_table->logTR_tab);
   maketab(TM_TR_MIN, TM_TR_MAX, NTM, rate_table->TM_TR_tab);
   rate_table->DlogTR = rate_table->logTR_tab[1] - rate_table->logTR_tab[0];
   rate_table->DTM_TR = rate_table->TM_TR_tab[1] - rate_table->TM_TR_tab[0];  

   for (i = 0; i < NTR; i++) {
       for (j = 0; j < NTM; j++) for (l = 0; l <= 1; l++) {
           fscanf(fA, "%le", &(rate_table->logAlpha_tab[l][j][i]));
           rate_table->logAlpha_tab[l][j][i] = log(rate_table->logAlpha_tab[l][j][i]);
       }
      fscanf(fR, "%le", &(rate_table->logR2p2s_tab[i]));
      rate_table->logR2p2s_tab[i] = log(rate_table->logR2p2s_tab[i]);   
   }
   fclose(fA);
   fclose(fR);
}

/************************************************************************************************
Interpolation of tabulated effective rates
To be more (slightly) efficient, not using the external interpolation routine.
Gets the correct rates for given fine-structure constant and electron mass.
- Modified May 2012:   - Accounts for different alpha_fs and m_e
                       - Also returns DAlpha[2], table of ALpha(Tm, Tr) - ALpha(Tr, Tr) 
INPUT TEMPERATURE ASSUMED TO BE ALREADY RESCALED FOR VALUES OF alpha_fs and me
************************************************************************************************/

void interpolate_rates(double Alpha[2], double DAlpha[2], double Beta[2], double *R2p2s, double TR, double TM_TR, 
                       HRATEEFF *rate_table, double fsR, double meR) {   
    double factor;
    unsigned l, k;
    long iTM, iTR;
    double frac1, frac2;
    double logTR;
    double coeff1[4], coeff2[4], temp[4];

    logTR = log(TR);

    /* Check if TM/TR is in the range tabulated */
    if (TM_TR < TM_TR_MIN || TM_TR > TM_TR_MAX) {
      fprintf(stderr, "Error: TM/TR = %f is out of range in interpolate_rates.\n", TM_TR);
      exit(1);
    }

    /* Check if log(TR) is in the range tabulated */
    if (TR < TR_MIN || TR > TR_MAX) {
      fprintf(stderr, "Error: TR = %f is out of range in interpolate_rates.\n", TR);
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
 
      Alpha[l] = square(fsR/meR)* exp(temp[0]*coeff1[0]+temp[1]*coeff1[1]
                                      +temp[2]*coeff1[2]+temp[3]*coeff1[3]);

      /* Alpha evaluated at Tm = Tr, for use in detailed balance */
      Beta[l] = square(fsR/meR)* exp(rate_table->logAlpha_tab[l][NTM-1][iTR-1]*coeff2[0] 
                                  +rate_table->logAlpha_tab[l][NTM-1][iTR]*coeff2[1]
                                  +rate_table->logAlpha_tab[l][NTM-1][iTR+1]*coeff2[2] 
                                  +rate_table->logAlpha_tab[l][NTM-1][iTR+2]*coeff2[3]); 

      DAlpha[l] = Alpha[l] - Beta[l];  /* Alpha(Tm, Tr) - Alpha(Tr, Tr) */

    }

    /* Beta obtained by detailed balance */
    /* factor = pow(2.0 * M_PI * mue *TR / hPc / hPc, 1.5)) * exp(-0.25*EI/TR) */
    factor = SAHA_FACT(fsR, meR) * TR*sqrt(TR) * exp(-0.25*EI/TR);
    Beta[0] *= factor;              /* 2s */        
    Beta[1] *= (factor/3.);         /* 2p */

    /* Effective 2p->2s rate */
    *R2p2s = fsR*fsR*fsR*fsR*fsR*meR * 
             exp(rate_table->logR2p2s_tab[iTR-1]*coeff2[0] 
                +rate_table->logR2p2s_tab[iTR]*coeff2[1]
                +rate_table->logR2p2s_tab[iTR+1]*coeff2[2] 
                +rate_table->logR2p2s_tab[iTR+2]*coeff2[3]);
  
}

/************************************************************************************************
Solves for the populations of the 2s and 2p states in steady-state, and returns dxe/dlna.
Uses standard rate for 2s-->1s decay and Sobolev for Lyman alpha (no feedback), 
and fully accounts for virtually all transitions among excited states through effective rates.
Inputs: xe, nH in cm^{-3}, H in s^{-1}, TM, TR in eV. Output: dxe/dlna
Changes May 2012: 
- Now use the analytic expressions using the generalized Peebles' C factors 
[see companion paper, Eqs. (43)-(46)]. 
- Also explicitly subtract nearly cancelling terms at early times. 
- Account for explicit dependence on alpha_fs and m_e
- Added explicit dependence on xHII, which is not necessarily equal to xe if Helium has not entirely recombined
************************************************************************************************/

double rec_HMLA_dxHIIdlna(double xe, double xHII, double nH, double H, double TM, double TR, 
                          HRATEEFF *rate_table, double fsR, double meR){

   double Alpha[2], DAlpha[2], Beta[2], R2p2s, RLya;   
   double Gamma_2s, Gamma_2p, C2s, C2p, s, Dxe2;
   double ratio;

   ratio = TM/TR;
   rescale_T(&TR, fsR, meR);
   TM = ratio * TR;   /* This way ensure that TM<=TR is preserved */   
  
   interpolate_rates(Alpha, DAlpha, Beta, &R2p2s, TR, TM/TR, rate_table, fsR, meR);  

   RLya = LYA_FACT(fsR, meR) *H/nH/(1.-xHII);   /* 8 PI H/(3 nH x1s lambda_Lya^3) */  
   
   /* Effective inverse lifetimes of 2s and 2p states */
   Gamma_2s = Beta[0] + 3.*R2p2s + L2s_rescaled(fsR, meR);     
   Gamma_2p = Beta[1] + R2p2s + RLya;     

   /* Generalization of Peebles' C factor */
   C2s = (L2s_rescaled(fsR, meR) + 3.*R2p2s * RLya/Gamma_2p)/(Gamma_2s - 3.*R2p2s*R2p2s/Gamma_2p);
   C2p = (RLya  +   R2p2s * L2s_rescaled(fsR, meR)/Gamma_2s)/(Gamma_2p - 3.*R2p2s*R2p2s/Gamma_2s);
   
   s    = SAHA_FACT(fsR, meR) *TR*sqrt(TR) *exp(-EI/TR)/nH;
   Dxe2 = xe*xHII - s*(1.-xHII); /* xe^2 - xe^2[Saha eq with 1s] -- gives more compact expressions */

   return -nH/H *( (s*(1.-xHII)*DAlpha[0] + Alpha[0]*Dxe2)*C2s + (s*(1.-xHII)*DAlpha[1] + Alpha[1]*Dxe2)*C2p );
    
}

/*********************************************************************************************
Read two-photon rates from table and store in structure.
**********************************************************************************************/
  
void read_twog_params(TWO_PHOTON_PARAMS *twog){ 
   
   FILE *fA;
   unsigned b;
   double L2s1s_current, max_DLNA, DlnE;

   fA = fopen(TWOG_FILE, "r");  

   for (b = 0; b < NVIRT; b++) { 
      fscanf(fA, "%le", &(twog->Eb_tab[b]));
      fscanf(fA, "%le", &(twog->A1s_tab[b]));
      fscanf(fA, "%le", &(twog->A2s_tab[b]));
      fscanf(fA, "%le", &(twog->A3s3d_tab[b]));
      fscanf(fA, "%le", &(twog->A4s4d_tab[b]));
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

   /* Added May 2012: check right away that the chosen time-step is small enough given the tabulated spectrum */

   max_DLNA = 1.;
   for (b = 0; b < NSUBLYA-1; b++) {
     DlnE = log(twog->Eb_tab[b+1]/twog->Eb_tab[b]);
     max_DLNA = max_DLNA < DlnE ? max_DLNA : DlnE;
   }
   DlnE = log(E21/twog->Eb_tab[NSUBLYA-1]); 
   max_DLNA = max_DLNA < DlnE ? max_DLNA : DlnE;
   DlnE = log(twog->Eb_tab[NSUBLYA]/E21); 
   max_DLNA = max_DLNA < DlnE ? max_DLNA : DlnE;
   for (b = NSUBLYA; b < NSUBLYB-1; b++) {
     DlnE = log(twog->Eb_tab[b+1]/twog->Eb_tab[b]);
     max_DLNA = max_DLNA < DlnE ? max_DLNA : DlnE;
   } 
   DlnE = log(E31/twog->Eb_tab[NSUBLYB-1]); 
   max_DLNA = max_DLNA < DlnE ? max_DLNA : DlnE;
   DlnE = log(twog->Eb_tab[NSUBLYB]/E31); 
   max_DLNA = max_DLNA < DlnE ? max_DLNA : DlnE;
   for (b = NSUBLYB; b < NVIRT-1; b++) {
     DlnE = log(twog->Eb_tab[b+1]/twog->Eb_tab[b]);
     max_DLNA = max_DLNA < DlnE ? max_DLNA : DlnE;
   } 
   DlnE = log(E41/twog->Eb_tab[NVIRT-1]); 
   max_DLNA = max_DLNA < DlnE ? max_DLNA : DlnE; 

   if (max_DLNA < DLNA) {
       fprintf(stderr, "Error: the time-step used (DLNA = %E) is too large for the used frequency grid (max DLNA = %E). Change DLNA accordingly in history.h.\n", DLNA, max_DLNA);
       exit(1); 
   }
   
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
 as well as the source vectors sr, sv, given Dxe = xe - xe[Saha] as an input. 
WITH DIFFUSION. Tvv[0][b] is the diagonal element Tbb, Tvv[1][b] = T{b,b-1} and Tvv[2][b] = T{b,b+1} 
Also, computes and stores the optical depths Delta tau_b for future use

Modified May 2012: - now uses the photon distortion 
instead of absolute photon occupation number.
- Accounts for explicit dependence on alpha_fs and m_e 
INPUT TEMPERATURE ASSUMED TO BE ALREADY RESCALED FOR VALUES OF alpha_fs and me
- Added explicit dependence on xHII, which is not necessarily equal to xe if Helium has not entirely recombined
**********************************************************************************************************/

void populateTS_2photon(double Trr[2][2], double *Trv[2], double *Tvr[2], double *Tvv[3], 
                        double sr[2], double sv[NVIRT], double Dtau[NVIRT],
                        double xe, double xHII, double TM, double TR, double nH, double H, HRATEEFF *rate_table,
                        TWO_PHOTON_PARAMS *twog, double Dfplus[NVIRT], double Dfplus_Ly[], 
                        double Alpha[2], double DAlpha[2], double Beta[2], double fsR, double meR) {

   double R2p2s;
   unsigned b;
   double  RLya, Gammab, one_minus_Pib, dbfact, x1s, s, Dxe2;
  
   double *Aup, *Adn;
   double A2p_up, A2p_dn; 

   double rescale2g, rescalediff;

   /*** Added May 2012: rescalings for dependence on alpha and me ***/
   rescale2g   = square(fsR*fsR*fsR*fsR)*meR;  /* for two-photon rates */
   rescalediff = rescale2g * fsR*fsR*meR;      
         /* diffusion rate ~ two-photon rate * TM ~ two-photon rate * alpha^2 me * rescaled(TM) [which is assumed as an input] */   
   

   Aup = create_1D_array(NVIRT);  
   Adn = create_1D_array(NVIRT);

   x1s  = 1.-xHII;
   s    = SAHA_FACT(fsR, meR) *TR*sqrt(TR) *exp(-EI/TR)/nH;
   Dxe2 = xe*xHII - s*x1s;
  
   RLya = LYA_FACT(fsR, meR) *H /nH/x1s;   /*8 PI H/(3 nH x1s lambda_Lya^3) */ 
 
   interpolate_rates(Alpha, DAlpha, Beta, &R2p2s, TR, TM / TR, rate_table, fsR, meR); 
    

  /****** 2s row and column ******/

   Trr[0][0] = Beta[0] + 3.*R2p2s
             + 3.* RLya * (1.664786871919931 *exp(-E32/TR)   /* Ly-beta escape */              
	                     + 1.953125 *exp(-E42/TR));          /* Ly-gamma escape */
			                  
   Trr[0][1] = -R2p2s;
   sr[0]     = nH * (s*x1s*DAlpha[0] + Alpha[0]*Dxe2) + 3.* RLya * x1s * 1.664786871919931 *Dfplus_Ly[1]; 
                          
   #if (EFFECT_A == 0) /* Standard treatment of 2s-->1s two-photon decays */
       Trr[0][0] += L2s1s*rescale2g;  /* rescaled for alpha, me */
   #endif
  
     
   /****** 2p row and column ******/

   Trr[1][1] = Beta[1] + R2p2s + RLya;			   
   Trr[1][0] = -3.*R2p2s;
   sr[1]     = nH * (s*x1s*DAlpha[1] + Alpha[1]*Dxe2) + 3.*RLya * x1s * Dfplus_Ly[0]; 
                                         
    
   /***** Two-photon transitions: populating Trv, Tvr and updating Trr ******/

   for (b = 0; b < NVIRT; b++) {
       dbfact = exp((twog->Eb_tab[b] - E21)/TR);

       Trr[0][0] -= Tvr[0][b] = -rescale2g*twog->A2s_tab[b]/fabs(exp((twog->Eb_tab[b] - E21)/TR)-1.);
       Trv[0][b]  = Tvr[0][b] *dbfact;
  
       Trr[1][1] -= Tvr[1][b] = -exp(-E32/TR)/3. * rescale2g*twog->A3s3d_tab[b]/fabs(exp((twog->Eb_tab[b] - E31)/TR)-1.)
                                -exp(-E42/TR)/3. * rescale2g*twog->A4s4d_tab[b]/fabs(exp((twog->Eb_tab[b] - E41)/TR)-1.); 
       Trv[1][b] = Tvr[1][b] *3.*dbfact;        
   }

    /****** Tvv and sv. Accounting for DIFFUSION ******/

    populate_Diffusion(Aup, Adn, &A2p_up, &A2p_dn, TM, twog->Eb_tab, twog->A1s_tab); 
     
    /*** Added May 2012: rescale for dependence on alpha and me ***/
    A2p_up *= rescalediff;
    A2p_dn *= rescalediff; 
    for (b = 0; b < NVIRT; b++) {
          Aup[b] *= rescalediff;    
          Adn[b] *= rescalediff;
    }

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

      Dtau[b] = Gammab * x1s * cube(hPc/twog->Eb_tab[b]/fsR/fsR/meR) * nH /8. /M_PI /H; 
        /* Rescaled for alpha, me*/

      one_minus_Pib = Dtau[b] > 1e-6 ? 1.- (1.-exp(-Dtau[b]))/Dtau[b] : Dtau[b]/2. - square(Dtau[b])/6.;
      Tvv[0][b] = Dtau[b] > 0.? Gammab/one_minus_Pib : 2./(x1s * cube(hPc/twog->Eb_tab[b]/fsR/fsR/meR) * nH /8. /M_PI /H);  /* Added May 2012: proper limit Dtau->0 */
      sv[b]  = Tvv[0][b] * x1s * Dfplus[b] * (1.-one_minus_Pib);         
       
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
   for (b = 0; b < NVIRT; b++) {
     xv[b] = Tvv_inv_sv[b] - Tvv_inv_Tvr[0][b]*xr[0] - Tvv_inv_Tvr[1][b]*xr[1];   
   }
 
  /** Free memory **/
   for (i = 0; i < 2; i++) free(Tvv_inv_Tvr[i]); 
   free(Tvv_inv_sv);   
  
}

/*********************************************************************************************
Interpolation of the photon distortion used to get f+ from f- at a higher frequency bin and earlier time. 
Use a simple linear interpolation so the spectrum is always positive.
Added May 2012
*********************************************************************************************/

double interp_Dfnu(double lna_start, double dlna, double *ytab, unsigned int iz, double lna){  
   long ind;
   double frac;

   /* If iz = 0 or 1, radiation field at earlier times is still thermal. 
      Also thermal if iz > 1 and lna < lna_start. */
   if (iz == 0 || iz == 1 || lna < lna_start) return 0.;
   
   /* Check if in range */
   if (lna >= lna_start + dlna*(iz-1)) {
       fprintf(stderr, "Error in interp_Dfnu: lna-value out of range in interpolation routine\n"); 
       fprintf(stderr, "The time-step used is probably too large\n");
       exit(1);
    }
 
    /* If iz >= 2, do a linear interpolation so the spectrum is always positive */
    ind  = (long) floor((lna-lna_start)/dlna);
    frac = (lna-lna_start)/dlna - ind;
    
    return (1.-frac)*ytab[ind] + frac*ytab[ind+1];  

}  

/*************************************************************************************************************
Obtain fplus at each bin, given the history of fminus (simple free-streaming). iz is the current time step.
fminus[0..iz-1] is known.
Assume the Lyman lines are optically thick
Dfminus_hist is a NVIRT by nz array of previous Delta f(nu_b - epsilon)(z)
Dfminus_Ly_hist is a 3 by nz array of previous Deltaf(nu - epsilon) redward of Ly alpha, beta and gamma lines

Changed May 2012: Now using the interpolation function interp_Dfnu, which only interpolates 
                    over 2 nearest neighbors, which ensures that the distortion is always positive
*************************************************************************************************************/

void fplus_from_fminus(double Dfplus[NVIRT], double Dfplus_Ly[], double **Dfminus_hist, double *Dfminus_Ly_hist[], 
                       double TR, double zstart, unsigned iz, double z, double Eb_tab[NVIRT]) {
   unsigned b;
   double ainv, lna_start, zp1;
  
   zp1 = 1.+z;
   lna_start = -log(1.+zstart);

   /*** Bins below Lyman alpha ***/
   for (b = 0; b < NSUBLYA-1; b++) {
      ainv = zp1*Eb_tab[b+1]/Eb_tab[b]; 
      Dfplus[b] = interp_Dfnu(lna_start, DLNA, Dfminus_hist[b+1], iz, -log(ainv));
   }
 
   /*** highest bin below Ly-alpha: feedback from optically thick Ly-alpha ***/
   b = NSUBLYA-1; 
   ainv = zp1*E21/Eb_tab[b];
   Dfplus[b] = interp_Dfnu(lna_start, DLNA, Dfminus_Ly_hist[0], iz, -log(ainv));
      
   /*** incoming photon occupation number at Lyman alpha ***/ 
   b = NSUBLYA;     /* next highest bin */
   ainv = zp1*Eb_tab[b]/E21;
   Dfplus_Ly[0] = interp_Dfnu(lna_start, DLNA, Dfminus_hist[b], iz, -log(ainv)); 
   
   /*** Bins between Lyman alpha and beta ***/   
   for (b = NSUBLYA; b < NSUBLYB-1; b++) {
     ainv = zp1*Eb_tab[b+1]/Eb_tab[b]; 
     Dfplus[b] = interp_Dfnu(lna_start, DLNA, Dfminus_hist[b+1], iz, -log(ainv));
   }
   
   /*** highest bin below Ly-beta: feedback from Ly-beta ***/
   b = NSUBLYB-1; 
   ainv = zp1*E31/Eb_tab[b];
   Dfplus[b] = interp_Dfnu(lna_start, DLNA, Dfminus_Ly_hist[1], iz, -log(ainv));
   
   /*** incoming photon occupation number at Lyman beta ***/ 
   b = NSUBLYB;     /* next highest bin */
   ainv = zp1*Eb_tab[b]/E31;
   Dfplus_Ly[1] = interp_Dfnu(lna_start, DLNA, Dfminus_hist[b], iz, -log(ainv)); 
  
   /*** Bins between Lyman beta and gamma ***/   
   for (b = NSUBLYB; b < NVIRT-1; b++) {
     ainv = zp1*Eb_tab[b+1]/Eb_tab[b]; 
     Dfplus[b] = interp_Dfnu(lna_start, DLNA, Dfminus_hist[b+1], iz, -log(ainv));
   }

   /*** highest energy bin: feedback from Ly-gamma ***/
   b = NVIRT-1;
   ainv = zp1*E41/Eb_tab[b];
   Dfplus[b] = interp_Dfnu(lna_start, DLNA, Dfminus_Ly_hist[2], iz, -log(ainv));
             
}

/******************************************************************************************************************
dxe/dlna when including two-photon processes. 
Assume fminus[0..iz-1] is known. Update fminus[iz]

Modified May 2012: 
- now use the photon distortion instead of absolute photon occupation number
- Accounts for explicit dependence on alpha_fs and m_e
- Added Dfnu_hist as a variable. Will contain the *average* distortion within each bin
******************************************************************************************************************/

double rec_HMLA_2photon_dxHIIdlna(double xe, double xHII, double nH, double H, double TM, double TR, 
                                  HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog,  
                                  double **Dfminus_hist, double *Dfminus_Ly_hist[], double **Dfnu_hist,
                                  double zstart, unsigned iz, double z, double fsR, double meR){
   
   double xr[2], xv[NVIRT], Dfplus[NVIRT], Dfplus_Ly[2]; /* Assume incoming radiation blueward of Ly-gamma is Blackbody */
   double dxedlna, one_minus_Pib, one_minus_exptau, Dfeq, s, x1s, Dxe2;
   unsigned b, i;
   
   double Trr[2][2];
   double *Trv[2];
   double *Tvr[2];
   double *Tvv[3];
   double sr[2];
   double sv[NVIRT];
   double Dtau[NVIRT]; 
   double Alpha[2], DAlpha[2], Beta[2];
   double ratio;

   ratio = TM/TR;
   rescale_T(&TR, fsR, meR);
   TM = ratio * TR;   /* This way ensure that TM<=TR is preserved */   
     
   
   for (i = 0; i < 2; i++) Trv[i] = create_1D_array(NVIRT);
   for (i = 0; i < 2; i++) Tvr[i] = create_1D_array(NVIRT);
   for (i = 0; i < 3; i++) Tvv[i] = create_1D_array(NVIRT); 

   /* Redshift photon occupation number from previous times and higher energy bins */
   fplus_from_fminus(Dfplus, Dfplus_Ly, Dfminus_hist, Dfminus_Ly_hist, TR, zstart, iz, z, twog->Eb_tab);
  
   /* Compute real-real, real-virtual and virtual-virtual transition rates */
   populateTS_2photon(Trr, Trv, Tvr, Tvv, sr, sv, Dtau, xe, xHII, TM, TR, nH, H, rate_table, 
                      twog, Dfplus, Dfplus_Ly, Alpha, DAlpha, Beta, fsR, meR);
    
   /* Solve for the population of the real and virtual states 
      (in fact, for the difference xi - xi[eq with 1s]) */   
   solve_real_virt(xr, xv, Trr, Trv, Tvr, Tvv, sr, sv);

   /* Obtain xe_dot */
   x1s  = 1.-xHII;
   s    = SAHA_FACT(fsR, meR) *TR*sqrt(TR) *exp(-EI/TR)/nH;
   Dxe2 = xe*xHII - s*x1s; 

   dxedlna = -(nH *(s*x1s*DAlpha[0] + Alpha[0]*Dxe2) - xr[0]*Beta[0] 
              +nH *(s*x1s*DAlpha[1] + Alpha[1]*Dxe2) - xr[1]*Beta[1])/H; 
    
   /* Update fminuses */
 
   for (b = 0; b < NVIRT; b++) {
     if (Dtau[b] > 1e-30) {
         one_minus_Pib = Dtau[b] > 1e-6 ? 1.- (1.-exp(-Dtau[b]))/Dtau[b] : Dtau[b]/2. - square(Dtau[b])/6.;
         Dfeq  = -xr[0]*Tvr[0][b] - xr[1]*Tvr[1][b];
         Dfeq -=  (b == 0       ?  xv[1]*Tvv[2][0]:
                   b == NVIRT-1 ?  xv[NVIRT-2]*Tvv[1][NVIRT-1]:    
                   xv[b+1]*Tvv[2][b] + xv[b-1]*Tvv[1][b]);        
         Dfeq /= x1s*one_minus_Pib*Tvv[0][b];
         one_minus_exptau = Dtau[b] > 1e-6 ? 1.-exp(-Dtau[b]) : Dtau[b] - square(Dtau[b])/2.;               
                                
         Dfminus_hist[b][iz] = Dfplus[b] + (Dfeq - Dfplus[b])*one_minus_exptau;
     }
     else Dfminus_hist[b][iz] = Dfplus[b];
   }
    
   Dfminus_Ly_hist[0][iz] = xr[1]/3./x1s;    
   Dfminus_Ly_hist[1][iz] = xr[0]/x1s*exp(-E32/TR); 
   Dfminus_Ly_hist[2][iz] = xr[0]/x1s*exp(-E42/TR); 
   
   /* Average radiation field in each bin */
   for (b = 0; b < NVIRT; b++) Dfnu_hist[b][iz] = xv[b]/x1s;   

   for (i = 0; i < 2; i++) free(Trv[i]);
   for (i = 0; i < 2; i++) free(Tvr[i]);
   for (i = 0; i < 3; i++) free(Tvv[i]);

   return dxedlna;
 
}

/*****************************************************************************************
Single function that handles all cases depending on the switch. 
When the two-photon rate is required, makes a simple estimate of the probability of 
photoionization from n=2. If larger than PION_MAX (set in history.h), 
use full radiative transfer equations, otherwise just use the simple EMLA.
******************************************************************************************/

double rec_dxHIIdlna(int model, double xe, double xHII, double nH, double H, double TM, double TR, 
                     HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog, double **Dfminus_hist, double *Dfminus_Ly_hist[], 
                     double **Dfnu_hist, double zstart, unsigned iz, double z, double fsR, double meR){
   
    if      (model == PEEBLES)  return rec_TLA_dxHIIdlna(xe, xHII, nH, H, TM, TR, 1.00, fsR, meR);
    else if (model == RECFAST)  return rec_TLA_dxHIIdlna(xe, xHII, nH, H, TM, TR, 1.14, fsR, meR);
    else if (model == EMLA2s2p) return rec_HMLA_dxHIIdlna(xe, xHII, nH, H, TM, TR, rate_table, fsR, meR);
    else if (model == FULL)     return rec_HMLA_2photon_dxHIIdlna(xe, xHII, nH, H, TM, TR, rate_table, twog, Dfminus_hist,  
                                                                  Dfminus_Ly_hist, Dfnu_hist, zstart, iz, z, fsR, meR);
    else {
      fprintf(stderr, "Error in rec_dxedlna: model = %i is undefined.\n", model);
      exit(1);
    }
}

/*********************************************************************************************/
