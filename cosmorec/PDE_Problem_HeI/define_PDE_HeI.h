//=====================================================================================
// Author Jeffrey Fung and Jens Chluba Feb/March 2011
//=====================================================================================

//=====================================================================================
//
// Purpose: Setup routine for PDE of form
//
// du/dt = A(x,t) d^2u/dx^2 + B(x,t) du/dx + C(x,t) u + D(x,t)
//
// for the helium diffusion problem. Electron scattering, line scattering can be 
// treated, as well as resonance emission and absorption using Lorentzian profiles
//
//=====================================================================================

#ifndef DEFINE_PDE_HEI_H
#define DEFINE_PDE_HEI_H

#include <vector>

#include "Cosmos.h"
#include "HeI_Atom.h"
#include "Voigtprofiles.h"

using namespace std;

//=====================================================================================
//
// setting different run-flags
//
//=====================================================================================
void HeI_switch_on_line_scattering();
void HeI_switch_off_line_scattering();

void HeI_switch_on_e_scattering();
void HeI_switch_off_e_scattering();

void HeI_switch_on_line_em_abs();
void HeI_switch_off_line_em_abs();

void HeI_switch_on_HI_abs();
void HeI_switch_off_HI_abs();

void HeI_switch_on_2s1s();
void HeI_switch_off_2s1s();

//=====================================================================================
struct PDE_solver_functions_HeI
{
    double nu2s1s;
    int emission_2s1s;
    
    double (* HI_Xe)(double z);
    double (* HI_rho)(double z);
    double (* HI_X1s)(double z);
    //
    double (* HeI_X1s)(double z);
    double (* HeI_Xi)(double z, int i);
    double (* HeI_Dnem_i)(double z, int i);
    double (* HeI_pd_i)(double z, int i);
    
    vector<double> exp_x;
    vector<double> sig_Lyc_x;
    vector<double> phi_2s1s_x;
    vector<double> Kernel_plus;
    vector<double> Kernel_minus;

    vector<int> res_index_of_lines;
    vector<vector<double> > Voigt_profiles;
    vector<vector<double> > dVoigt_profiles_dx;
    vector<double> Voigt_profiles_aV;
    vector<double> Voigt_profiles_Dnuj1s_nuj1s;
    vector<double> Voigt_profiles_nui1s;
    vector<double> Voigt_profiles_Ai1s;
    vector<double> Voigt_profiles_pd;
    vector<double> Voigt_profiles_w;
    vector<double> dum_factors;
};

//=====================================================================================
// for numerical purposes it can be beneficial to add a very small diffusion term, to
// damp small scale perturbations. The effect of electron scattering is very small on
// Xe and if one does not want e-scattering be important reducing its efficiency by ~100 
// will make it extremely insignificant but allows stabilizing the system. This can be 
// controlled with the e_sc_efficiency parameter. ==1 means full electron scattering.
// Default is 1.0 if this function is not called.
//=====================================================================================
void set_electron_scattering_efficiency_HeI(double val);

//=====================================================================================
// profile memory
//=====================================================================================
void init_profile_memory_HeI(int nxi, int nres, PDE_solver_functions_HeI &PDE_funcs);

//=====================================================================================
// Set PDE Data
//=====================================================================================
void Setup_PDE_Data_HeI(Cosmos &cosm, Gas_of_HeI_Atoms &HeIA, 
                        PDE_solver_functions_HeI &PDE_funcs);

//=====================================================================================
// PDE setup
//=====================================================================================
void def_PDE_HeI(double z, const vector<double> &xi, 
                 vector<double> &Ai, vector<double> &Bi, 
                 vector<double> &Ci, vector<double> &Di);

//=====================================================================================
// profile access (index i is the line index)
//=====================================================================================
double phi_i_pd_HeI(int i);
double phi_i_aV_HeI(int i);
double phi_i_HeI(int i, int k);
double phi_i_nui1s_HeI(int i);
double phi_i_Ai1s_HeI(int i);
double phi_i_Dnu_HeI(int i);
double phi_i_w_HeI(int i);
Voigtprofile_Dawson& get_HeI_profile(int line);

//=====================================================================================
// access two HeI Singlet 2s-1s emission profile (with stimulated emission)
//=====================================================================================
double phi_2s1s_HeI(int k);

//=====================================================================================
// PDE setup with electron scattering treated using kernel (added 24.02.2011 JC)
//=====================================================================================
void electron_kernel_corr_func(double z, 
                               const vector<double> &xi, 
                               const vector<double> &yi, 
                               vector<double> &corr_v);

//=====================================================================================
void compute_kernel_factors(int ix, double nu21, double Tm, 
                            const vector<double> &xarr, 
                            vector<double> &Kernel_p, 
                            vector<double> &Kernel_m, 
                            double eps);

double compute_kernel_Integral_polint(int ix, double nu21, double Tm, 
                                      const vector<double> &xi, 
                                      const vector<double> &yi);

#endif

//=====================================================================================
//=====================================================================================
