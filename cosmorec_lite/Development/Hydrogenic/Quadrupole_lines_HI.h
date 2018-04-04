//==================================================================================================
// Author: Jens Chluba 
//
// first implementation: May 2011
// last modification: May 2011
// 
// Routines to compute electic quadrupole transitions in hydrogenic 
// atoms. The formulae were taken from Grin et al. 2010 and Hey 2006.
//==================================================================================================

#ifndef QUADRUPOLE_LINES_HI_H
#define QUADRUPOLE_LINES_HI_H

//==================================================================================================
// quadrupole transition rates for hydrogen.
// (the mass of is assumed to be MH==mp, where mp is the proton mass)
// routines seems to work very well even for n>~500
//==================================================================================================
double A_E2_Quadrupole(int n, int l, int np, int lp);

//==================================================================================================
// quadrupole transition for hydrogenic atoms of 
// charge Z and mass Ma=Np*mp, where mp is the mass a proton. 
//==================================================================================================
double A_E2_Quadrupole(int Z, double Np, int n, int l, int np, int lp);

//==================================================================================================
// dipole and quadrupole transition rates for hydrogenic atoms of 
// charge Z and mass Ma=Np*mp, where mp is the mass a proton. 
// Here all rates out of the level are computed at once to save time.
//==================================================================================================
void A_E1_E2_Quadrupole(int Z, double Np, int n, int l, int np,
                        double &Am1, double &Ap1, 
                        double &Am2, double &Amp, double &Ap2);

#endif

//==================================================================================================
//==================================================================================================
