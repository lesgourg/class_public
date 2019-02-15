/*************************************************************************************************/
/*                 HYREC: Hydrogen and Helium Recombination Code                                 */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                              */
/*                                                                                               */
/*         hyrec.c: main module                                                                  */
/*                                                                                               */
/*         Version: November 2011                                                                 */
/*         Revision history:                                                                     */
/*            - written November 2010                                                            */
/*            - January 2011: changed various switches (notably for post-Saha expansions)        */
/*                             so that they remain valid for arbitrary cosmologies               */
/*            - November 2011: extended integration down to z = 0 with Peeble's model for z < 20   */
/*                             changed dTm/dlna so it can be used at all times                     */
/*************************************************************************************************/ 
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "hyrec.h"

int main(void) {

   REC_COSMOPARAMS param;
   HRATEEFF rate_table;
   TWO_PHOTON_PARAMS twog_params;
   double *xe_output, *Tm_output;
   long iz;

   double dz = -1;
   unsigned nz = 8001;
   double z, xe, Tm;

      
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

   /* Compute the recombination history */
   xe_output = (double*)malloc((size_t)(param.nz*sizeof(double)));
   Tm_output = (double*)malloc((size_t)(param.nz*sizeof(double)));
  
   rec_build_history(&param, &rate_table, &twog_params, xe_output, Tm_output);

   /* Interpolate at the desired output redshifts */
   for(iz=0; iz<nz; iz++) {
       z = param.zstart + dz * iz;    /* print output every dz */
       xe = rec_interp1d(-log(1.+param.zstart), param.dlna, xe_output, param.nz, -log(1.+z));
       Tm = rec_interp1d(-log(1.+param.zstart), param.dlna, Tm_output, param.nz, -log(1.+z));
       printf("%7.2lf %15.13lf %15.13lf\n", z, xe, Tm/param.T0/(1.+z));
   }
   
 
    /* Cleanup */
    free((char*)xe_output);
    free((char*)Tm_output);
    free(rate_table.logTR_tab);
    free(rate_table.TM_TR_tab);
    free_2D_array(rate_table.logAlpha_tab[0], NTM);
    free_2D_array(rate_table.logAlpha_tab[1], NTM);
    free(rate_table.logR2p2s_tab);

  return(0);
}
