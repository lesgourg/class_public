//====================================================================================================================
// Author Jens Chluba Feb-May 2011
//====================================================================================================================

//====================================================================================================================
//
// Purpose: Setup routine for PDE of form
//
// du/dt = A(x,t) d^2u/dx^2 + B(x,t) du/dx + C(x,t) u + D(x,t)
//
// for the helium diffusion problem. Electron scattering, line scattering can be 
// treated, as well as resonance emission and absorption using Lorentzian profiles
//
// 22.04.2011: added simple support for openmp
//====================================================================================================================

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>

#include "physical_consts.h"
#include "routines.h"
#include "Cosmos.h"
#include "HeI_Atom.h"
#include "get_effective_rates.HeI.h"
#include "Voigtprofiles.h"

#include "define_PDE_HeI.h"
#include "Patterson.h"
#include "Compton_Kernel.h"

using namespace std;

//====================================================================================================================
//
// run-flags
//
//====================================================================================================================
struct define_PDE_HeI_runflags
{
    bool line_sc;         // resonance scattering on/off
    bool el_sc;           // electron scattering on/off
    bool line_em_abs;     // line emission & absorption on/off
    bool A2s1s_em_abs;    // HeI Singlet 2s-1s emission & absorption on/off
    bool HI_absorption;   // HI continuum absorption on/off
};

define_PDE_HeI_runflags def_HeI_run={1, 1, 1, 1, 1};  // default setting

//====================================================================================================================
void HeI_switch_on_line_scattering(){ def_HeI_run.line_sc=1; return; }
void HeI_switch_off_line_scattering(){ def_HeI_run.line_sc=0; return; }

void HeI_switch_on_e_scattering(){ def_HeI_run.el_sc=1; return; }
void HeI_switch_off_e_scattering(){ def_HeI_run.el_sc=0; return; }

void HeI_switch_on_line_em_abs(){ def_HeI_run.line_em_abs=1; return; }
void HeI_switch_off_line_em_abs(){ def_HeI_run.line_em_abs=0; return; }

void HeI_switch_on_HI_abs(){ def_HeI_run.HI_absorption=1; return; }
void HeI_switch_off_HI_abs(){ def_HeI_run.HI_absorption=0; return; }

void HeI_switch_on_2s1s(){ def_HeI_run.A2s1s_em_abs=1; return; }
void HeI_switch_off_2s1s(){ def_HeI_run.A2s1s_em_abs=0; return; }

//====================================================================================================================
// for numerical purposes it can be beneficial to add a very small diffusion term, to
// damp small scale perturbations. The effect of electron scattering (in the Fokker-Planck 
// approximation) on Xe is rather notable for Helium recombinationon. If one does not want 
// e-scattering be important reducing its efficiency by ~100-1000 will make it extremely 
// insignificant but allows stabilizing the system. This can be controlled with the e_sc_efficiency 
// parameter. ==1 means full electron scattering. Default is 1.0 if this function is not called.
//====================================================================================================================
double e_sc_efficiency_HeI=1.0;
void set_electron_scattering_efficiency_HeI(double val){ e_sc_efficiency_HeI=val; return; }

//====================================================================================================================
//
// PDE Data
//
//====================================================================================================================
struct PDE_Setup_Data_HeI
{
    Cosmos *cosm;
    Gas_of_HeI_Atoms *HeI;
    PDE_solver_functions_HeI *PDE_funcs;
};

static PDE_Setup_Data_HeI PDE_Setup_D_HeI;

//====================================================================================================================
// memory
//====================================================================================================================
void init_profile_memory_HeI(int nxi, vector<double> &VP)
{
    VP.clear();
    VP.resize(nxi); 
    return;
}

//====================================================================================================================
void init_profile_memory_HeI(int nxi, int nres, PDE_solver_functions_HeI &PDE_funcs)
{
    init_profile_memory_HeI(nxi, PDE_funcs.exp_x);
    init_profile_memory_HeI(nxi, PDE_funcs.sig_Lyc_x);
    init_profile_memory_HeI(nxi, PDE_funcs.phi_2s1s_x);
    init_profile_memory_HeI(nxi, PDE_funcs.Kernel_plus);
    init_profile_memory_HeI(nxi, PDE_funcs.Kernel_minus);
    
    //==============================================================
    // allocate memory for Voigt profiles 
    //==============================================================
    PDE_funcs.Voigt_profiles.clear();
    
    for(int k=0; k<nres; k++)
    {
        PDE_funcs.Voigt_profiles.push_back(PDE_funcs.exp_x);
        PDE_funcs.dVoigt_profiles_dx.push_back(PDE_funcs.exp_x);
    }   
    
    init_profile_memory_HeI(nres, PDE_funcs.Voigt_profiles_aV);
    init_profile_memory_HeI(nres, PDE_funcs.Voigt_profiles_Dnuj1s_nuj1s);
    init_profile_memory_HeI(nres, PDE_funcs.Voigt_profiles_pd);
    init_profile_memory_HeI(nres, PDE_funcs.Voigt_profiles_nui1s);
    init_profile_memory_HeI(nres, PDE_funcs.Voigt_profiles_Ai1s);
    init_profile_memory_HeI(nres, PDE_funcs.Voigt_profiles_w);
    init_profile_memory_HeI(nres, PDE_funcs.dum_factors);

    //==============================================================
    // Temperature independent values for profiles
    //==============================================================
    for(int k=0; k<(int)PDE_funcs.Voigt_profiles_aV.size(); k++)
    {
        Voigtprofile_Dawson p=get_HeI_profile(k);
        
        PDE_funcs.Voigt_profiles_nui1s[k]=p.Get_nu21();
        PDE_funcs.Voigt_profiles_Ai1s[k] =p.Get_A21();
        
        int HeI_i=get_HeI_index(PDE_funcs.res_index_of_lines[k]);
        PDE_funcs.Voigt_profiles_w[k]=PDE_Setup_D_HeI.HeI->Get_gw(HeI_i)/PDE_Setup_D_HeI.HeI->Get_gw(0);
    }
    
    
    return;
}

//====================================================================================================================
// Setup for PDE Data
//====================================================================================================================
void Setup_PDE_Data_HeI(Cosmos &cosm, Gas_of_HeI_Atoms &HeIA, PDE_solver_functions_HeI &PDE_funcs)
{
    PDE_Setup_D_HeI.cosm=&cosm;
    PDE_Setup_D_HeI.HeI=&HeIA;
    PDE_Setup_D_HeI.PDE_funcs=&PDE_funcs;
    
    return;
}

//====================================================================================================================
// scattering terms
//====================================================================================================================
void add_Diff_and_recoil_f_HeI(double D_x, double dlnD_x_dx, double xifac, double &Ai, double &Bi, double &Ci)
{
    //=============================================================================
    // Doppler broadening
    //=============================================================================
    Ai+=D_x;        
    Bi+=dlnD_x_dx;
    //=============================================================================
    
    //=============================================================================
    // recoil term
    //=============================================================================
    Bi+=D_x*xifac;
    Ci+=dlnD_x_dx*xifac;
    //=============================================================================

    return;
}

//====================================================================================================================
// electron scattering terms
//====================================================================================================================
void add_electron_scattering_HeI(double z, double nu21, double Te, double Ne, double Hz, const vector<double> &xi, 
                                 vector<double> &Ai, vector<double> &Bi, vector<double> &Ci)
{
    //----------------------------------------------------------
    // xifac= h nu_Lya/ k Te; important for recoil term
    //----------------------------------------------------------
    double xifac=const_h_kb*nu21/Te;
    double kappa_electron_sc=-const_sigT*Ne*(const_kB*Te/const_me_gr/const_cl)/Hz/(1.0+z);
    kappa_electron_sc*=e_sc_efficiency_HeI;
    
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=0; k<(int)xi.size(); k++) 
    {
        double x=xi[k];
        double D_x=kappa_electron_sc*x*x;
        double dlnD_x_dx=D_x*4.0/x;
        add_Diff_and_recoil_f_HeI(D_x, dlnD_x_dx, xifac, Ai[k], Bi[k], Ci[k]);
    }       
    
    return;
}

//====================================================================================================================
// resonance scattering and em/abs terms
//====================================================================================================================
double pd_n_func_HeI(int line, double z)
{ return PDE_Setup_D_HeI.PDE_funcs->HeI_pd_i(z, PDE_Setup_D_HeI.PDE_funcs->res_index_of_lines[line]); }

double Dnem_n_func_HeI(int line, double z)
{ return PDE_Setup_D_HeI.PDE_funcs->HeI_Dnem_i(z, PDE_Setup_D_HeI.PDE_funcs->res_index_of_lines[line]); }

//====================================================================================================================
// exponential-factor exp(-h nu/k Tg)
//====================================================================================================================
void update_exp_x_HeI(double Tg, double nu21, const vector<double> &xi)
{
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=0; k<(int)xi.size(); k++) 
        PDE_Setup_D_HeI.PDE_funcs->exp_x[k]=exp(-const_h_kb*nu21*xi[k]/Tg);    
    
    return;
}

//====================================================================================================================
// update Voigt profile
//====================================================================================================================
void update_given_profile_HeI(int k, double Tm, double psc, double nu21, 
                              Voigtprofile_Dawson &p, const vector<double> &xi)
{
    double aV=PDE_Setup_D_HeI.PDE_funcs->Voigt_profiles_aV[k]=p.aVoigt(Tm);
    PDE_Setup_D_HeI.PDE_funcs->Voigt_profiles_Dnuj1s_nuj1s[k]=p.DnuT_nu12(Tm);
    double A21=PDE_Setup_D_HeI.PDE_funcs->Voigt_profiles_Ai1s[k];
    //
    PDE_Setup_D_HeI.PDE_funcs->Voigt_profiles_aV[k]=aV*A21/p.Get_Gamma()/psc;
    
    //================================================================================================================
    // Ly-n profiles
    //================================================================================================================
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int j=0; j<(int)xi.size(); j++) 
    {
        double x=xi[j];
        double xD=p.nu2x(nu21*x, Tm);
        PDE_Setup_D_HeI.PDE_funcs->Voigt_profiles    [k][j]=p.phi(xD, aV);
        PDE_Setup_D_HeI.PDE_funcs->dVoigt_profiles_dx[k][j]=p.dphi_dx(xD, aV);
    }
    
    return;
}

//====================================================================================================================
Voigtprofile_Dawson& get_HeI_profile(int line)
{
    int i_HeI=get_HeI_index(PDE_Setup_D_HeI.PDE_funcs->res_index_of_lines[line]);
    int n=PDE_Setup_D_HeI.HeI->Get_n(i_HeI);
    int l=PDE_Setup_D_HeI.HeI->Get_l(i_HeI);
    int s=PDE_Setup_D_HeI.HeI->Get_S(i_HeI);
    
    if(l==1 && s==0) return PDE_Setup_D_HeI.HeI->nP_S_profile(n);
    else if(l==1 && s==1) return PDE_Setup_D_HeI.HeI->nP_T_profile(n);
    else if(l==2 && s==0) return PDE_Setup_D_HeI.HeI->nD_S_profile(n);
    
    return PDE_Setup_D_HeI.HeI->nP_S_profile(0);
}

void update_Voigt_profiles_HeI(double z, double nu21, const vector<double> &xi)
{
    double Tm=PDE_Setup_D_HeI.PDE_funcs->HI_rho(z)*PDE_Setup_D_HeI.cosm->TCMB(z);

    for(int k=0; k<(int)PDE_Setup_D_HeI.PDE_funcs->Voigt_profiles.size(); k++)
    {
        PDE_Setup_D_HeI.PDE_funcs->Voigt_profiles_pd[k]=pd_n_func_HeI(k, z);
        double psc_np=1.0-PDE_Setup_D_HeI.PDE_funcs->Voigt_profiles_pd[k];

        update_given_profile_HeI(k, Tm, psc_np, nu21, get_HeI_profile(k), xi);
    }
    
    return;
}

//====================================================================================================================
// HeI Singlet 2s-1s two photon decay profile function; normalised to 2
//====================================================================================================================
double phi_HeI_2s_decay_spectrum(double y)
{
	if(y<0.0000001 || y> 0.9999999) return 0.0;
    double w=y*(1.0-y);
    
    return 19.6039602*pow(w, 3)*( 1.742-7.2*w+12.8*w*w)/pow(w+0.03, 2);
}

void update_HeI_2s_1s_profile(double z, double nu21, const vector<double> &xi)
{
    double Tg=PDE_Setup_D_HeI.cosm->TCMB(z);
    double exp_2s1s=exp(-const_h_kb*PDE_Setup_D_HeI.PDE_funcs->nu2s1s/Tg);
    
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int j=0; j<PDE_Setup_D_HeI.PDE_funcs->emission_2s1s; j++) 
    {
        double y=xi[j]*nu21/PDE_Setup_D_HeI.PDE_funcs->nu2s1s;
        double phi=phi_HeI_2s_decay_spectrum(y);
        double npl_fac=1.0/(1.0-PDE_Setup_D_HeI.PDE_funcs->exp_x[j]);
        double nplp_fac=1.0/(1.0-exp_2s1s/PDE_Setup_D_HeI.PDE_funcs->exp_x[j]);
        
        PDE_Setup_D_HeI.PDE_funcs->phi_2s1s_x[j]=const_HeI_A2s_1s*phi*npl_fac*nplp_fac
                                                   /PDE_Setup_D_HeI.PDE_funcs->nu2s1s;
    }
    
    return;
}

//====================================================================================================================
//
// access two HeI Singlet 2s-1s emission profile (with stimulated emission)
//
//====================================================================================================================
double phi_2s1s_HeI(int k){ return PDE_Setup_D_HeI.PDE_funcs->phi_2s1s_x[k]; }

//====================================================================================================================
//
// access line properties (i==line index)
//
//====================================================================================================================
double phi_i_HeI(int i, int k){ return PDE_Setup_D_HeI.PDE_funcs->Voigt_profiles[i][k]; }
double dphi_i_dx_HeI(int i, int k){ return PDE_Setup_D_HeI.PDE_funcs->dVoigt_profiles_dx[i][k]; }
double phi_i_Dnu_nu_HeI(int i){ return PDE_Setup_D_HeI.PDE_funcs->Voigt_profiles_Dnuj1s_nuj1s[i]; }
double phi_i_aV_HeI(int i){ return PDE_Setup_D_HeI.PDE_funcs->Voigt_profiles_aV[i]; }
double phi_i_pd_HeI(int i){ return PDE_Setup_D_HeI.PDE_funcs->Voigt_profiles_pd[i]; }
double phi_i_nui1s_HeI(int i){ return PDE_Setup_D_HeI.PDE_funcs->Voigt_profiles_nui1s[i]; }
double phi_i_Ai1s_HeI(int i){ return PDE_Setup_D_HeI.PDE_funcs->Voigt_profiles_Ai1s[i]; }
double phi_i_Dnu_HeI(int i){ return phi_i_nui1s_HeI(i)*PDE_Setup_D_HeI.PDE_funcs->Voigt_profiles_Dnuj1s_nuj1s[i]; }
double phi_i_w_HeI(int i){ return PDE_Setup_D_HeI.PDE_funcs->Voigt_profiles_w[i]; }

//====================================================================================================================
//
// PDE setup functions
//
//====================================================================================================================
void add_Ly_n_sc_only_HeI(int line, double z, double nu21, double Tg, double Te, double Ne, double N1s, double Hz, 
                          const vector<double> &xi, vector<double> &Ai, vector<double> &Bi, vector<double> &Ci)
{
    //----------------------------------------------------------
    // xifac= h nu_Lya/ k Te; important for recoil term
    //----------------------------------------------------------
    double xifac=const_h_kb*nu21/Te;    
    double nun1=phi_i_nui1s_HeI(line);
    
    //----------------------------------------------------------
    // resonance scattering
    //----------------------------------------------------------
    double lambda=const_cl/nu21;         // lambdaj1*(nuj1/nu21)
    double DnuD_nu=phi_i_Dnu_nu_HeI(line);
    double pd_np=phi_i_pd_HeI(line);
    double psc_np=1.0-pd_np;
    
    double const_sigr=0.5*phi_i_w_HeI(line)*lambda*lambda*phi_i_Ai1s_HeI(line)/FOURPI/phi_i_Dnu_HeI(line);
    double kappa_Ly_n_sc =-psc_np*const_sigr*N1s*(const_kB*Te/const_mHe4_gr/const_cl)/Hz/(1.0+z);
    
    //----------------------------------------------------------
    // add terms
    //----------------------------------------------------------
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=0; k<(int)xi.size(); k++) 
    {
        double x=xi[k];
        
        //------------------------------------------------------
        // line-scattering terms 
        //------------------------------------------------------
        double phi_x=phi_i_HeI(line, k);
        double dphi_x_dx=dphi_i_dx_HeI(line, k) / DnuD_nu * nu21/nun1;
        double D_x=kappa_Ly_n_sc*phi_x;
        double dlnD_x_dx=kappa_Ly_n_sc*(phi_x*2.0/x+dphi_x_dx);
        
        add_Diff_and_recoil_f_HeI(D_x, dlnD_x_dx, xifac, Ai[k], Bi[k], Ci[k]);
    }

    return;
}

//====================================================================================================================
void add_Ly_n_em_only_HeI(int line, double z, double nu21, double Tg, double Te, double Ne, double N1s, double Hz, 
                          const vector<double> &xi, vector<double> &Ci, vector<double> &Di)
{   
    //----------------------------------------------------------
    // resonance
    //----------------------------------------------------------
    double lambda=const_cl/nu21;         // lambdaj1*(nuj1/nu21)
    double pd_np=phi_i_pd_HeI(line);
    double const_sigr=0.5*phi_i_w_HeI(line)*lambda*lambda*phi_i_Ai1s_HeI(line)/FOURPI/phi_i_Dnu_HeI(line);
    
    //----------------------------------------------------------
    // Ly-n emission/absorption term
    //----------------------------------------------------------
    double nun1=phi_i_nui1s_HeI(line);
    double x_c=const_h_kb*nun1/Tg;
    double Dnem =Dnem_n_func_HeI(line, z)*nu21;
    //
    double kappa_Ly_n_abs=-pd_np*const_sigr*N1s*const_cl/Hz/(1.0+z); 
    double exp_dum=exp(-x_c);

#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=0; k<(int)xi.size(); k++) 
    {
        double x=xi[k];
        
        //------------------------------------------------------
        // prepare for line-scattering and em/abs terms 
        //------------------------------------------------------
        double phi_x=phi_i_HeI(line, k);
        
        //------------------------------------------------------
        // Source-term from emission & absorption 
        //------------------------------------------------------
        double exp_fac=exp_dum*(1.0/PDE_Setup_D_HeI.PDE_funcs->exp_x[k]-1.0);
        double D_x=kappa_Ly_n_abs*phi_x/x/x;
        
        //------------------------------------------------------
        Ci[k]+=-D_x*exp_fac;
        Di[k]+= D_x*Dnem;
        //------------------------------------------------------
    }
        
    return;
}

//====================================================================================================================
void add_Ly_n_sc_em_full_HeI(int line, double z, double nu21, double Tg, double Te, double Ne, double N1s, double Hz, 
                             const vector<double> &xi, vector<double> &Ai, vector<double> &Bi,
                             vector<double> &Ci, vector<double> &Di)
{
    //----------------------------------------------------------
    // xifac= h nu_Lya/ k Te; important for recoil term
    //----------------------------------------------------------
    double xifac=const_h_kb*nu21/Te;    
    double nun1=phi_i_nui1s_HeI(line);
    
    //----------------------------------------------------------
    // resonance scattering
    //----------------------------------------------------------
    double lambda=const_cl/nu21;         // lambdaj1*(nuj1/nu21)
    double DnuD_nu=phi_i_Dnu_nu_HeI(line);
    double pd_np=phi_i_pd_HeI(line);
    double psc_np=1.0-pd_np;
    double const_sigr=0.5*phi_i_w_HeI(line)*lambda*lambda*phi_i_Ai1s_HeI(line)/FOURPI/phi_i_Dnu_HeI(line);
    double kappa_Ly_n_sc =-psc_np*const_sigr*N1s*(const_kB*Te/const_mHe4_gr/const_cl)/Hz/(1.0+z);
    
    //----------------------------------------------------------
    // Ly-n emission/absorption term
    //----------------------------------------------------------
    double x_c=const_h_kb*nun1/Tg;
    double Dnem =Dnem_n_func_HeI(line, z)*nu21;
    //
    double kappa_Ly_n_abs=-pd_np*const_sigr*N1s*const_cl/Hz/(1.0+z); 
    double exp_dum=exp(-x_c);
    
    //----------------------------------------------------------
    // add terms
    //----------------------------------------------------------
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=0; k<(int)xi.size(); k++) 
    {
        double x=xi[k];
        
        //------------------------------------------------------
        // line-scattering terms 
        //------------------------------------------------------
        double phi_x=phi_i_HeI(line, k);
        double dphi_x_dx=dphi_i_dx_HeI(line, k) / DnuD_nu * nu21/nun1;
        double D_x=kappa_Ly_n_sc*phi_x;
        double dlnD_x_dx=kappa_Ly_n_sc*(phi_x*2.0/x+dphi_x_dx);
        
        add_Diff_and_recoil_f_HeI(D_x, dlnD_x_dx, xifac, Ai[k], Bi[k], Ci[k]);

        //------------------------------------------------------
        // Source-term from emission & absorption 
        //------------------------------------------------------
        double exp_fac=exp_dum*(1.0/PDE_Setup_D_HeI.PDE_funcs->exp_x[k]-1.0);
        D_x=kappa_Ly_n_abs*phi_x/x/x;
        
        //------------------------------------------------------
        Ci[k]+=-D_x*exp_fac;
        Di[k]+= D_x*Dnem;
        //------------------------------------------------------
    }
    
    return;
}


//====================================================================================================================
//
// HI absorption term
//
//====================================================================================================================
void add_HI_absorption_term(double z, double nu21, double Tg, double N1s, double Hz, 
                            const vector<double> &xi, vector<double> &Ci)
{
    double kappa_abs=-const_cl*N1s/Hz/(1.0+z);
  
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=0; k<(int)xi.size(); k++) 
        Ci[k]+=-kappa_abs*PDE_Setup_D_HeI.PDE_funcs->sig_Lyc_x[k];
    
    return;
}

//====================================================================================================================
//
// HeI Singlet 2s-1s term
//
//====================================================================================================================

void add_HeI_2s_1s_terms(double z, double nu21, double Tg, double N1s, double Hz, 
                         const vector<double> &xi, vector<double> &Ci, vector<double> &Di)
{
	update_HeI_2s_1s_profile(z, nu21, xi);
    
    //----------------------------------------------------------
	// 2s-1s emission/absorption term
	//----------------------------------------------------------
	double sig_2s=pow(const_cl/nu21, 2)/2.0/FOURPI;
	double kappa_2s1s=-sig_2s*N1s*const_cl/Hz/(1.0+z);
	//
	double Dnem_2s=PDE_Setup_D_HeI.PDE_funcs->HeI_Dnem_i(z, 0)*nu21;
	double x_c=const_h_kb*PDE_Setup_D_HeI.PDE_funcs->nu2s1s/Tg;
	double exp_dum=exp(-x_c);
	
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
	for(int k=0; k<PDE_Setup_D_HeI.PDE_funcs->emission_2s1s; k++) 
	{
		double x=xi[k];
        double phi_em=phi_2s1s_HeI(k);
		double exp_fac_2s=exp_dum*(1.0/PDE_Setup_D_HeI.PDE_funcs->exp_x[k]-1.0);
		double DD=kappa_2s1s*phi_em/x/x;
		
		//------------------------------------------------------
		Ci[k]+=-DD*exp_fac_2s;
		Di[k]+= DD*Dnem_2s;
		//------------------------------------------------------
	}
	
	return;
}

//====================================================================================================================
//
// main wrapper
//
//====================================================================================================================
void def_PDE_HeI(double z, const vector<double> &xi, 
                 vector<double> &Ai, vector<double> &Bi, 
                 vector<double> &Ci, vector<double> &Di)
{
    int nmax=PDE_Setup_D_HeI.PDE_funcs->Voigt_profiles_aV.size();
    
    //----------------------------------------------------------
    // xi = nu/nu_Lya
    //----------------------------------------------------------
    double nu21=PDE_Setup_D_HeI.HeI->nP_S_profile(2).Get_nu21();
    
    //----------------------------------------------------------
    // time-dependent variables 
    //----------------------------------------------------------
    double Hz=PDE_Setup_D_HeI.cosm->H(z);
    double NH=PDE_Setup_D_HeI.cosm->NH(z);
    double Tg=PDE_Setup_D_HeI.cosm->TCMB(z);
    
    //----------------------------------------------------------
    // solutions from rec-code
    //----------------------------------------------------------
    double Ne=NH*PDE_Setup_D_HeI.PDE_funcs->HI_Xe(z);
    double N1s=NH*PDE_Setup_D_HeI.PDE_funcs->HI_X1s(z);
    double NHeI1s=NH*PDE_Setup_D_HeI.PDE_funcs->HeI_X1s(z);
    double rho=PDE_Setup_D_HeI.PDE_funcs->HI_rho(z);
    double Te=Tg*rho;
    
    //----------------------------------------------------------
    // exp(x) is at least needed in calling routine
    //----------------------------------------------------------
    update_exp_x_HeI(Tg, nu21, xi);

    //----------------------------------------------------------
    // reset things
    //----------------------------------------------------------
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=0; k<(int)xi.size(); k++)
    {
        Ai[k]=Ci[k]=Di[k]=0.0;
        Bi[k]=-xi[k]/(1.0+z);  // Hubble term
    }
    
    //----------------------------------------------------------
    // electron scattering terms first (smallest term)
    //----------------------------------------------------------
    if(def_HeI_run.el_sc) add_electron_scattering_HeI(z, nu21, Te, Ne, Hz, xi, Ai, Bi, Ci);
    
    //----------------------------------------------------------
    // normal Ly-n scatterring and emission/absorption
    //----------------------------------------------------------
    update_Voigt_profiles_HeI(z, nu21, xi);
    
    if(def_HeI_run.line_sc || def_HeI_run.line_em_abs)
    {
        bool local_line_em_abs=def_HeI_run.line_em_abs;
        bool local_line_sc=def_HeI_run.line_sc;
        
        for(int n=0; n<nmax; n++)
        {
            if(local_line_sc && local_line_em_abs)
                add_Ly_n_sc_em_full_HeI(n, z, nu21, Tg, Te, Ne, NHeI1s, Hz, xi, Ai, Bi, Ci, Di);
            else if(local_line_sc && !local_line_em_abs)
                add_Ly_n_sc_only_HeI(n, z, nu21, Tg, Te, Ne, NHeI1s, Hz, xi, Ai, Bi, Ci);
            else if(!local_line_sc && local_line_em_abs)
                add_Ly_n_em_only_HeI(n, z, nu21, Tg, Te, Ne, NHeI1s, Hz, xi, Ci, Di);
        }
    }
        
    //----------------------------------------------------------
    // account for HI absorption
    //----------------------------------------------------------
    if(def_HeI_run.HI_absorption) add_HI_absorption_term(z, nu21, Tg, N1s, Hz, xi, Ci);
            
    //----------------------------------------------------------
    // Source-term HeI Singlet 2s-1s emission and absorption
    //----------------------------------------------------------
    if(def_HeI_run.A2s1s_em_abs) add_HeI_2s_1s_terms(z, nu21, Tg, NHeI1s, Hz, xi, Ci, Di);
    
    return;
}

//====================================================================================================================
//
// for computations with kernel corrections
//
//====================================================================================================================
#include "./define_PDE_HeI.Kernel.cpp"

//====================================================================================================================
//====================================================================================================================
