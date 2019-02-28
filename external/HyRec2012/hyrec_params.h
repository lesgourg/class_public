/*************************************************************************************************/
/*                 HYREC: Hydrogen and Helium Recombination Code                                 */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                              */
/*                                                                                               */
/*         hyrec_params.h: switches and accuracy parameters                                      */
/*         Version: May 2012                                                                     */
/*************************************************************************************************/



#define PROMPT 1      /* Set to zero to suppress initial prompts */

/**** Switch to choose the physical model used for hydrogen ****/ 

/* definitions */
#define PEEBLES   0    /* Peebles effective three-level atom */
#define RECFAST   1    /* Effective three-level atom for hydrogen with fudge factor F = 1.14 */
#define EMLA2s2p  2    /* Correct EMLA model, with standard decay rates from 2s and 2p only (accounts for nmax = infinity, l-resolved) */
#define FULL      3    /* All radiative transfer effects included. Additional switches in header file hydrogen.h */

/** Here is the switch **/

#define MODEL FULL     /* default setting: FULL */


/*** Switch to print the spectrum (nu_Ly_alpha/2 to nu_Ly_gamma). 
If so, a smoother spectrum can be obtained by using the higher accuracy settings below (but standard parameters are fine too). ***/

#define PRINT_SPEC 0                                                              
#define SPEC_FILE "hyrec_spectrum.dat"   /*File to which the spectrum is printed */
#define NSPEC 8                         /* NSPEC spectra for z evenly spaced between ZMIN_SPEC and ZMAX_SPEC will be printed */
#define ZMIN_SPEC 900.
#define ZMAX_SPEC 1600.


/*********** Numerical integration parameters ************/

#define ZSTART  8000.             /* Starting redshift */
#define ZEND    0.                /* End redshift */

#define DLNA          8.49e-5    /* Time step Delta ln(a). Must be smaller than 8.5e-5 if using radiative transfer */
#define DXHEII_MAX    1e-5       /* If xHeII - xHeII(Saha) < DXEHII_MAX, use post-Saha expansion for Helium. Lower value = higher accuracy. */
#define DXHII_MAX     3e-4       /* If xHII - xHII(Saha) < DXHII_MAX, use post-Saha expansion for Hydrogen. Switch to ODE integration after that.
                                        IMPORTANT: do not set to a lower value unless using a smaller time-step */
#define XHEII_MIN     1e-6       /* Stop considering Helium recombination once xHeII < XHEII_MIN */
#define DLNT_MAX      5e-4       /* Use the steady-state approximation for Tm as long as 1-Tm/Tr < DLNT_MAX, then switch to ODE integration */
#define PION_MAX      1e-2       /* When the probability of being ionized from n=2 becomes lower than PION_MAX, switch off radiative transfer calculation as it becomes irrelevant */

/*** Tables and parameters for radiative transfer calculation ***/

#define TWOG_FILE "two_photon_tables.dat" /* Maximum DLNA = 8.49e-5 */
#define NSUBLYA  140
#define NSUBLYB  271
#define NVIRT    311
#define NDIFF    80


/**** Set of higher-accuracy parameters to use if getting the spectrum (make a smoother spectrum)****/

/* #define DLNA       2e-5 */
/* #define DXHEII_MAX 1e-6 */
/* #define DXHII_MAX  2e-5 */
/* #define XHEII_MIN  1e-6 */
/* #define DLNT_MAX   5e-5 */
/* #define PION_MAX   1e-4 */

/* /\* Higher-resolution tables *\/ */
/* #define TWOG_FILE "two_photon_tables_hires.dat"  /\* maximum dlna = 8.47e-5 *\/ */
/* #define NSUBLYA  408 */
/* #define NSUBLYB  1323 */
/* #define NVIRT    1493 */
/* #define NDIFF    300 */

