//==================================================================================================
// Author Jens Chluba
//
// first implementation: July 2010
// last modification   : July 2014
//==================================================================================================
// 30.07.2014: Xi_Data is passed on read only
// 04.06.2011: added nD-1s transitions to transfer problem
// 26.03.2011: added Ly-n integrals & electron scattering Kernel support

#ifndef SOLVE_PDES_H
#define SOLVE_PDES_H

#include "Cosmos.h"
#include "Atom.h"

//==================================================================================================
// switch on/off HI 2s-1s channel; this overwrites all other settings for 2s-1s channel
//==================================================================================================
void switch_off_1s_2s_correction();
void switch_on_1s_2s_correction();

//==================================================================================================
// this allows to switch the 1s-2s absorption correction off; in the radiative transfer 
// everything is included, but the integral is only taken over the emission term. If the HI 2s-1s 
// channel is switch off, then calling this function has no effect.
//==================================================================================================
void switch_off_1s_2s_absorption_correction();

//==================================================================================================
int compute_DPesc_with_diffusion_equation_effective(vector<double> &DF_vec_z, 
                                                    vector<vector<double> > &DF_Ly_n, 
                                                    vector<vector<double> > &DF_2_gamma_vec, 
                                                    vector<vector<double> > &DF_Raman_vec, 
                                                    vector<vector<double> > &DF_nD1s_vec, 
                                                    vector<double> &DI1_2s_vec, 
                                                    double zs, double ze, 
                                                    int nmax_2g_corrs, int nmax_R_corrs,
                                                    Cosmos &cos, Gas_of_Atoms &HIA,
                                                    const vector<vector<double> > &HI_Solution,
                                                    int it_num);

#endif

//==================================================================================================
//==================================================================================================
