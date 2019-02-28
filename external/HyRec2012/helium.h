/*************************************************************************************************/
/*                 HYREC: Hydrogen and Helium Recombination Code                                 */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                              */
/*                                                                                               */
/*         helium.h: all functions related to helium recombination (and Saha equilibria)         */
/*                                                                                               */
/*         Version: May 2012  (first released November 2010)                                     */
/*************************************************************************************************/ 

double rec_xesaha_HeII_III(double nH0, double Tr0, double fHe, double z, double *xHeIII, double fsR, double meR);
double rec_saha_xHeII(double nH0, double Tr0, double fHe, double z, double fsR, double meR);
double rec_saha_xH1s(double xHeII, double nH0, double T0, double z, double fsR, double meR);
double rec_helium_dxHeIIdlna(double xH1s, double xHeII, double nH0, double Tr0, double fHe, double H, double z,
                             double fsR, double meR);

