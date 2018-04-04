//====================================================================================================================
// Authors Jens Chluba & Jeffrey Fung Feb-May 2011
//====================================================================================================================
// 30.07.2014: Xi_Data is passed on read only

#ifndef SOLVE_PDES_HEI_H
#define SOLVE_PDES_HEI_H

#include "Cosmos.h"
#include "Atom.h"
#include "HeI_Atom.h"

//====================================================================================================================
// switch on/off kernel approximation for electron scattering (default is off).
// comment: When the kernel correction is switched on the code takes rather long. The additional
// correction is small, so the default is recommended.
//====================================================================================================================
void switch_off_e_scat_kernel_HeI();
void switch_on_e_scat_kernel_HeI();

//====================================================================================================================
// switch on/off HeI singlet 2s-1s channel; this overwrites all other settings for 2s-1s channel
//====================================================================================================================
void switch_off_1s_2s_HeI_correction();
void switch_on_1s_2s_HeI_correction();

//====================================================================================================================
// this allows to switch the 1s-2s absorption correction off; in the radiative transfer
// everything is included, but the integral is only taken over the emission term. If the HeI 2s-1s 
// channel is switch off, then calling this function has no effect.
//====================================================================================================================
void switch_off_1s_2s_HeI_absorption_correction();

//====================================================================================================================
int compute_HeI_DPesc_with_diffusion_equation_effective(int nHeI_Trans,
                                                        vector<double> &DF_vec_z, 
                                                        vector<vector<double> > &DF_vec,
                                                        vector<double> &DI1_2s_vec,
                                                        vector<int> &HeI_level_indices,
                                                        double zs, double ze, 
                                                        Cosmos &cos, 
                                                        Gas_of_HeI_Atoms &HeIAtoms,
                                                        const vector<vector<double> > &HI_HeI_Solution,
                                                        int it_num, bool last_it);

#endif

//====================================================================================================================
//====================================================================================================================
