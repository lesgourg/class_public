//==================================================================================================
// Author Jens Chluba July 2010
// last modification: Aug 2012
//==================================================================================================

//==================================================================================================
//
// Purpose: Setup routine for PDE of form
//
// du/dt = A(x,t) d^2u/dx^2 + B(x,t) du/dx + C(x,t) u + D(x,t)
//
//==================================================================================================

#ifndef DEFINE_PDE_H
#define DEFINE_PDE_H

#include <vector>

#include "Cosmos.h"
#include "Atom.h"

using namespace std;

//==================================================================================================
//
// setting different run-flags
//
//==================================================================================================
void HI_switch_on_line_scattering();
void HI_switch_off_line_scattering();

void HI_switch_on_e_scattering();
void HI_switch_off_e_scattering();

void HI_switch_on_line_em_abs();
void HI_switch_off_line_em_abs();

void HI_switch_on_2s1s();
void HI_switch_off_2s1s();

//==================================================================================================
struct Raman_2g_profiles
{
    vector<double> ratio;
    vector<double> profile;
};

//==================================================================================================
struct PDE_solver_functions
{
    double (* HI_Xe)(double z);
    double (* HI_rho)(double z);
    double (* HI_X1s)(double z);
    //
    double (* HI_Xnl)(double z, int n, int l);
    double (* HI_Dnem_nl)(double z, int n, int l);
    double (* HI_pd_nl)(double z, int n, int l);
        
    vector<int> index_emission; // for given n this is the index for 
                                // two-photon emission & Raman-absorption
    int index_2;
    //
    vector<Raman_2g_profiles > Raman_ns_1s;
    vector<Raman_2g_profiles > Raman_nd_1s;
    //
    vector<Raman_2g_profiles > two_g_ns_1s;
    vector<Raman_2g_profiles > two_g_nd_1s;
    //
    vector<vector<double> > Voigt_profiles_A_npns;
    vector<vector<double> > Voigt_profiles_A_npnd;
    //
    vector<double> exp_x;
    vector<vector<double> > Voigt_profiles;
    vector<vector<double> > HP_Voigt_profiles;
    vector<double> Voigt_profiles_aV;
    vector<double> Voigt_profiles_Dnuj1s_nuj1s;
    vector<double> Voigt_profiles_pd;
    vector<double> dum_factors;
    // for Lyman-series 1+1 photon parts
    vector<double> Voigt_profiles_pd_eff;
    vector<double> Voigt_profiles_Dnem_eff;
};

//==================================================================================================
// for numerical purposes it can be beneficial to add a very small diffusion term, to
// damp small scale perturbations. The effect of electron scattering is very small on
// Xe and if one does not want e-scattering be important reducing its efficiency by ~100 
// will make it extremely insignificant but allows stabilizing the system. This can be 
// controlled with the e_sc_efficiency parameter. ==1 means full electron scattering.
// Default is 1.0 if this function is not called.
//==================================================================================================
void set_electron_scattering_efficiency(double val);

//==================================================================================================
// setting for two-gamma terms
//==================================================================================================
int  Get_nmax_two_g_correction();
void switch_off_two_g_corrections();
void switch_on_two_g_corrections(int tg_nmax, int nS);

//==================================================================================================
// setting for Raman terms
//==================================================================================================
int  Get_nmax_Raman_correction();
void switch_off_Raman_corrections();
void switch_on_Raman_corrections(int R_nmax, int nS);

//==================================================================================================
// profile memory
//==================================================================================================
void init_profile_memory(int nxi, int nres, PDE_solver_functions &PDE_funcs);

//==================================================================================================
// Set PDE Data
//==================================================================================================
void Setup_PDE_Data(Cosmos &cosm, Gas_of_Atoms &HIA, PDE_solver_functions &PDE_funcs);

//==================================================================================================
// PDE setup
//==================================================================================================
void def_PDE_Lyn_and_2s1s(double z, const vector<double> &xi, vector<double> &Ai, 
                          vector<double> &Bi, vector<double> &Ci, vector<double> &Di);

//==================================================================================================
// PDE setup with electron scattering treated using kernel (added 30.02.2011 JC)
//==================================================================================================
void electron_kernel_corr_func_HI(double z, 
                                  const vector<double> &xi, 
                                  const vector<double> &yi, 
                                  vector<double> &corr_v);

//==================================================================================================
// profile access
//==================================================================================================
double phi_ns1s_2g(int n, int k);
double phi_nd1s_2g(int n, int k);
//
double phi_ns1s_Raman(int n, int k);
double phi_nd1s_Raman(int n, int k);

double phi_Ly_n(int n, int k);
//
double phi_Ly_n_Dnu_nu(int n);
double phi_Ly_n_Dnu(int n);
//
double phi_Ly_n_aV(int n);
double phi_Ly_n_pd(int n);

#endif

//==================================================================================================
//==================================================================================================
