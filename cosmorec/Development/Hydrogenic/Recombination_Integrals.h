//======================================================================================
// Author: Jens Chluba 
// last modification: Oct 2010
// 
// Purpose: carry out the integrals over photoionization cross sections
//======================================================================================

#ifndef RECOMBINATION_INTS_H
#define RECOMBINATION_INTS_H

#include <iostream>
#include <string>

using namespace std;

//======================================================================================
// routines for integration for photo-ionization & recombination rates
//======================================================================================
double Integrate_Ric_Rci_JC_Patterson(int choice, double lga, double lgb, double epsrel, double epsabs, double (*fptr)(double));  // f(x) in log scale
double Integrate_Ric_Rci_JC_Patterson(int choice, double lga, double lgb, double epsrel, double epsabs, double (*fptr)(double, void *p), void *p);  // f(x) in log scale

double Integrate_Ric_Rci_JC_Chebyshev(int choice, double xl, double xu, double epsrel, double epsabs, double (*fptr)(double), double xc=1.0); // f(x) in linear scale
double Integrate_Ric_Rci_JC_Chebyshev_log(int choice, double lgxl, double lgxu, double epsrel, double epsabs, double (*fptr)(double), double lgxc=1.0); // f(x) in log scale
double Integrate_Ric_Rci_JC_Chebyshev_log(int choice, string mess, double lgxl, double lgxu, double epsrel, double epsabs, double (*fptr)(double), double lgxc=1.0); // f(x) in log scale

#endif
