//=======================================================================================
// Author Jens Chluba May 2009
//=======================================================================================
// 30.07.2014: Xi_Data is passed on read only

#ifndef LOAD_POP_HI_H
#define LOAD_POP_HI_H

#include "Atom.h"
#include <string>

//=======================================================================================
// Load the data for Xi
//=======================================================================================
void Load_Xi_HI(int nS, double zs, double ze, Gas_of_Atoms &H_Atoms, string filename="");
void clear_Xi_HI();
void Set_Load_Populations_HI_verbosity(int v);

//=======================================================================================
// to directly pass on the solutions for the populations
//=======================================================================================
void compute_Xi_HI_splines(double zs, double ze, 
                           Gas_of_Atoms &H_Atoms, 
                           const vector<vector<double> > &X_Data);

//=======================================================================================
// functions for access
//=======================================================================================
double calc_HI_Xe(double z);
double calc_HI_rho(double z);
double calc_HI_X1s(double z);

//=======================================================================================
double calc_HI_Xnl(double z, int n, int l);
double calc_HI_Xi(double z, int i);

#endif

//=======================================================================================
//=======================================================================================
