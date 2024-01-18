/*************************************************************************************************/
/*                 HYREC-2: Hydrogen and Helium Recombination Code                               */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (2010-17)                              */
/*           with contributions from Nanoom Lee (2020)                                           */
/*                                                                                               */
/*         helium.h: all functions related to helium recombination (and Saha equilibria)         */
/*                                                                                               */
/*************************************************************************************************/

#ifndef __HELIUM__
#define __HELIUM__

#include "hydrogen.h"

#define SIZE_ErrorM   2048

double rec_xesaha_HeII_III(REC_COSMOPARAMS *cosmo, double z, double *xHeIII);
double rec_saha_xHeII(REC_COSMOPARAMS *cosmo, double z);
double rec_saha_xH1s(REC_COSMOPARAMS *cosmo, double z, double xHeII);
double rec_helium_dxHeIIdlna(HYREC_DATA *data, double z, double xH1s, double xHeII, double H);

#endif
