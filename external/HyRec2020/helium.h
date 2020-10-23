/*************************************************************************************************/
/*                 HYREC-2: Hydrogen and Helium Recombination Code                               */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (2010-17)                              */
/*           with contributions from Nanoom Lee (2020)                                           */
/*                                                                                               */
/*         helium.h: all functions related to helium recombination (and Saha equilibria)         */
/*                                                                                               */
/*************************************************************************************************/ 
#define SIZE_ErrorM   2048

double rec_xesaha_HeII_III(double nH0, double Tr0, double fHe, double z, double *xHeIII, double fsR, double meR);
double rec_saha_xHeII(double nH0, double Tr0, double fHe, double z, double fsR, double meR);
double rec_saha_xH1s(double xHeII, double nH0, double T0, double z, double fsR, double meR);
double rec_helium_dxHeIIdlna(double xH1s, double xHeII, double nH0, double Tr0, double fHe,
			     double H, double z, double fsR, double meR, int *error, char error_message[SIZE_ErrorM]);
double xe_PostSahaHe(double nH0, double Tr0, double fHe, double H, double z);

