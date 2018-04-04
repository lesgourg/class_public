//====================================================================
// Author: Jens Chluba 
// last modification: Oct 2010
// 
// The routines are based on Storey & Hummer, 1991. Notation and 
// recursion relations can be found there. However, we modified these 
// slightly to achieve better stability.
//====================================================================

#ifndef OSCILLATOR_STRENGTH_H
#define OSCILLATOR_STRENGTH_H

//====================================================================
// oscillator strength and transition rate for hydrogen.
// (the mass of is assumed to be MH==mp, where mp is the proton mass)
// routines seems to work very well even for n>~500
//====================================================================
double A_SH(int &n, int &l, int &np, int &lp);
double f_SH(int &n, int &l, int &np, int &lp);

//====================================================================
// oscillator strength and transition rate for hydrogenic atoms of 
// charge Z and mass Ma=Np*mp, where mp is the mass a proton. 
// Drake 1996, Atomic, Molecular & Optical Physics
//====================================================================
double A_SH(int Z, double Np, int n, int l, int np, int lp);
double f_SH(int Z, double Np, int n, int l, int np, int lp);

#endif
