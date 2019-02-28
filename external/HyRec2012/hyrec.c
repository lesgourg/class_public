/*********************************************************************************************************/
/*                          HYREC: Hydrogen and Helium Recombination Code                                */
/*                     Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                          */
/*                                                                                                       */
/*         hyrec.c: main module                                                                          */
/*                                                                                                       */
/*         Version: May 2012                                                                             */
/*                                                                                                       */
/*         Revision history:                                                                             */
/*            - written November 2010                                                                    */
/*            - January 2011: changed various switches (notably for post-Saha expansions)                */
/*                             so that they remain valid for arbitrary cosmologies                       */
/*            - November 2011: extended integration down to z = 0 with Peeble's model for z < 20         */
/*                             changed dTm/dlna so it can be called at all times                         */
/*            - May 2012: - included explicit dependence on fine-structure constant and electron mass    */
/*                          - the user can now extract the Lyman-lines distortion (up to Ly-gamma)       */
/*********************************************************************************************************/ 

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "hyrectools.h"
#include "helium.h"
#include "hydrogen.h"
#include "history.h"
#include "hyrec_params.h"



int main(void) {

   REC_COSMOPARAMS param;
   HRATEEFF rate_table;
   TWO_PHOTON_PARAMS twog_params;
   double *xe_output, *Tm_output, **Dfnu_hist, *Dfminus_Ly_hist[3];
   long iz;

   double dz = -1;
   unsigned nz = 8001;
   double z, xe, Tm;

   #if PRINT_SPEC == 1
      FILE *fp;
      double z_spec[NSPEC];
      unsigned b;
      long izmax_spec;
      double prefact, Dfnu;
   #endif

   /* Build effective rate table */
   rate_table.logTR_tab = create_1D_array(NTR);
   rate_table.TM_TR_tab = create_1D_array(NTM);
   rate_table.logAlpha_tab[0] = create_2D_array(NTM, NTR);
   rate_table.logAlpha_tab[1] = create_2D_array(NTM, NTR);
   rate_table.logR2p2s_tab = create_1D_array(NTR);

   read_rates(&rate_table);
  
   /* Read two-photon rate tables */
   read_twog_params(&twog_params);

   /* Get cosmological parameters */
   rec_get_cosmoparam(stdin, stderr, &param);
   
   /* allocate memory for output (only nzrt for spectrum arrays) */
   xe_output          = create_1D_array(param.nz);
   Tm_output          = create_1D_array(param.nz);
   Dfnu_hist          = create_2D_array(NVIRT, param.nzrt);
   Dfminus_Ly_hist[0] = create_1D_array(param.nzrt);        /* Ly-alpha */
   Dfminus_Ly_hist[1] = create_1D_array(param.nzrt);        /* Ly-beta  */
   Dfminus_Ly_hist[2] = create_1D_array(param.nzrt);        /* Ly-gamma */

  

   /* Compute the recombination history */
   rec_build_history(&param, &rate_table, &twog_params, xe_output, Tm_output, Dfnu_hist, Dfminus_Ly_hist);
    
   
   /* Interpolate at the desired output redshifts */
   for(iz=0; iz<nz; iz++) {
       z = ZSTART + dz * iz;    /* print output every dz */
       xe = rec_interp1d(-log(1.+ZSTART), DLNA, xe_output, param.nz, -log(1.+z));
       Tm = rec_interp1d(-log(1.+ZSTART), DLNA, Tm_output, param.nz, -log(1.+z));
       printf("%7.2lf %15.15lf %15.13lf\n", z, xe, Tm/param.T0/(1.+z));
    }
  

    /***** Printing out the Lyman-lines spectral distortion *****/

   #if PRINT_SPEC == 1
   for (iz = 0; iz < NSPEC; iz++) z_spec[iz] = ZMIN_SPEC + (ZMAX_SPEC-ZMIN_SPEC)/(NSPEC-1) *iz;  /* Redshifts at which spectrum is printed */
        izmax_spec = (long) floor(2+log((1.+param.zH0)/(1.+ZMIN_SPEC))/DLNA);    /* Max index used for interpolation */
      
        fp = fopen(SPEC_FILE, "w");
        
        fprintf(fp, "Number of spectral distortion photons per hydrogen atom per log-frequency interval, (8 pi nu^3)/(c^3 nH) Delta f_nu"); 
        fprintf(fp, "nu/nuLya  z=");
        for (iz = 0; iz < NSPEC; iz++) fprintf(fp, " %E", z_spec[iz]);
        fprintf(fp, "\n");
      
        for (b = 0; b < NSUBLYA; b++) { /* Sub-Lyman alpha bins */
           fprintf(fp, "%E", twog_params.Eb_tab[b]/E21);
           for (iz = 0; iz < NSPEC; iz++) {
	     Dfnu = interp_Dfnu(-log(1.+param.zH0), DLNA, Dfnu_hist[b], izmax_spec, -log(1.+z_spec[iz]));
             prefact = 8.*M_PI*cube(twog_params.Eb_tab[b]/E21)/(param.nH0 * cube((1.+z_spec[iz]) * 1216e-10));
             fprintf(fp," %E", prefact*Dfnu);
           }
           fprintf(fp, "\n");
        }
        fprintf(fp, "%E", 1.);    /* Lyman-alpha */
        for (iz = 0; iz < NSPEC; iz++) {
    	    Dfnu = interp_Dfnu(-log(1.+param.zH0), DLNA, Dfminus_Ly_hist[0], izmax_spec, -log(1.+z_spec[iz]));
            prefact = 8.*M_PI/(param.nH0 * cube((1.+z_spec[iz]) * 1216e-10));
            fprintf(fp," %E", prefact*Dfnu);
        }
        fprintf(fp, "\n");
        for (; b < NSUBLYB; b++) { /* Bins between Ly-alpha and Ly-beta */
           fprintf(fp, "%E", twog_params.Eb_tab[b]/E21);
           for (iz = 0; iz < NSPEC; iz++) {
              Dfnu = interp_Dfnu(-log(1.+param.zH0), DLNA, Dfnu_hist[b], izmax_spec, -log(1.+z_spec[iz]));
              prefact = 8.*M_PI*cube(twog_params.Eb_tab[b]/E21)/(param.nH0 * cube((1.+z_spec[iz]) * 1216e-10));
              fprintf(fp," %E", prefact*Dfnu);                
           }
           fprintf(fp, "\n");
        }
        fprintf(fp, "%E", E31/E21);    /* Lyman-beta */
        for (iz = 0; iz < NSPEC; iz++) {
             Dfnu =  interp_Dfnu(-log(1.+param.zH0), DLNA, Dfminus_Ly_hist[1], izmax_spec, -log(1.+z_spec[iz])); 
             prefact = 8.*M_PI*cube(E31/E21)/(param.nH0 * cube((1.+z_spec[iz]) * 1216e-10));
             fprintf(fp," %E", prefact*Dfnu);
        }
        fprintf(fp, "\n");
        for (; b < NVIRT; b++) {   /* Bins between Ly-beta and Ly-gamma */
           fprintf(fp, "%E", twog_params.Eb_tab[b]/E21);
           for (iz = 0; iz < NSPEC; iz++) {
              Dfnu = interp_Dfnu(-log(1.+param.zH0), DLNA, Dfnu_hist[b], izmax_spec, -log(1.+z_spec[iz]));
              prefact = 8.*M_PI*cube(twog_params.Eb_tab[b]/E21)/(param.nH0 * cube((1.+z_spec[iz]) * 1216e-10));
              fprintf(fp," %E", prefact*Dfnu);
           }
           fprintf(fp, "\n");
        }
        fprintf(fp, "%E", E41/E21);    /* Lyman-gamma */
        for (iz = 0; iz < NSPEC; iz++) {
            Dfnu = interp_Dfnu(-log(1.+param.zH0), DLNA, Dfminus_Ly_hist[2], izmax_spec, -log(1.+z_spec[iz]));
            prefact = 8.*M_PI*cube(E41/E21)/(param.nH0 * cube((1.+z_spec[iz]) * 1216e-10));    
            fprintf(fp," %E", prefact*Dfnu);
        }
        fprintf(fp, "\n");

        fclose(fp);

    #endif
    /*** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** ***/  


   /* Cleanup */
    free(rate_table.logTR_tab);
    free(rate_table.TM_TR_tab);
    free_2D_array(rate_table.logAlpha_tab[0], NTM);
    free_2D_array(rate_table.logAlpha_tab[1], NTM);
    free(rate_table.logR2p2s_tab);
    free(xe_output);
    free(Tm_output);
    free_2D_array(Dfnu_hist, NVIRT);
    free(Dfminus_Ly_hist[0]);
    free(Dfminus_Ly_hist[1]);
    free(Dfminus_Ly_hist[2]);

  return(0);
}

/****************************************************************************************************/
