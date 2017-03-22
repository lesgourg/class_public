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

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>

#include "physical_consts.h"
#include "routines.h"
#include "Cosmos.h"
#include "Atom.h"

#include "define_PDE.h"
#include "define_PDE_HeI.h"

using namespace std;

int define_PDE_verbosity=0;

//==================================================================================================
//
// run-flags
//
//==================================================================================================
struct define_PDE_HI_runflags
{
    bool line_sc;         // resonance scattering on/off
    bool el_sc;           // electron scattering on/off
    bool line_em_abs;     // line emission & absorption on/off
    bool A2s1s_em_abs;    // HI 2s-1s emission & absorption on/off
    bool nD_Quadrupole;   // HI nD-1s quadrupole lines
};

define_PDE_HI_runflags def_HI_run={1, 1, 1, 1, 0};  // default setting

//==================================================================================================
void HI_switch_on_line_scattering(){ def_HI_run.line_sc=1; return; }
void HI_switch_off_line_scattering(){ def_HI_run.line_sc=0; return; }

void HI_switch_on_e_scattering(){ def_HI_run.el_sc=1; return; }
void HI_switch_off_e_scattering(){ def_HI_run.el_sc=0; return; }

void HI_switch_on_line_em_abs(){ def_HI_run.line_em_abs=1; return; }
void HI_switch_off_line_em_abs(){ def_HI_run.line_em_abs=0; return; }

void HI_switch_on_2s1s(){ def_HI_run.A2s1s_em_abs=1; return; }
void HI_switch_off_2s1s(){ def_HI_run.A2s1s_em_abs=0; return; }

//==================================================================================================
//
// more delicate run-flags; These are just for detailed computations, and should
// normally not be changed.
//
//==================================================================================================
bool recoil_l=1;
bool diffusion_l=1;
//
bool recoil_e=1;
bool diffusion_e=1;
//==================================================================================================

//==================================================================================================
bool two_g_process=0;
int two_g_process_nmax=0;

//==================================================================================================
bool Raman_process=0;
int Raman_process_nmax=0;

//==================================================================================================
bool use_full_width=1;
//==================================================================================================

//==================================================================================================
// for numerical purposes it can be beneficial to add a very small diffusion term, to
// damp small scale perturbations. The effect of electron scattering is very small on
// Xe and if one does not want e-scattering be important reducing its efficiency by ~100 
// will make it extremely insignificant but allows stabilizing the system. This can be 
// controlled with the e_sc_efficiency parameter. ==1 means full electron scattering.
// Default is 1.0 if this function is not called.
//==================================================================================================
double e_sc_efficiency=1.0;
void set_electron_scattering_efficiency(double val){ e_sc_efficiency=val; return; }

//==================================================================================================
//
// setting for two-photon terms
//
//==================================================================================================
int  Get_nmax_two_g_correction(){ return two_g_process_nmax; }

void switch_off_two_g_corrections()
{
    two_g_process=0;
    two_g_process_nmax=0;
    
    if(define_PDE_verbosity>=1) 
        cout << " switch_off_two_g_corrections :: Switching two-gamma corrections off. " << endl;
    
    return;
}

void switch_on_two_g_corrections(int tg_nmax, int nS)
{
    two_g_process=1;
    two_g_process_nmax=(int)min(8, min(tg_nmax, nS));
    
    if(define_PDE_verbosity>=1)
    {
        cout << " switch_on_two_g_corrections :: Switching on two-gamma corrections. " << endl;
        cout << " switch_on_two_g_corrections :: nmax will be set to " 
             << two_g_process_nmax << endl;
    }
    
    return;
}



//==================================================================================================
//
// setting for Raman terms
//
//==================================================================================================
int  Get_nmax_Raman_correction(){ return Raman_process_nmax; }

void switch_off_Raman_corrections()
{
    Raman_process=0;
    Raman_process_nmax=0;
    
    if(define_PDE_verbosity>=1)
        cout << " switch_off_Raman_corrections :: Switching Raman corrections off. " << endl;
    
    return;
}

void switch_on_Raman_corrections(int R_nmax, int nS)
{
    Raman_process=1;
    Raman_process_nmax=(int)min(7, min(R_nmax, nS-1));
    
    if(define_PDE_verbosity>=1)
    {
        cout << " switch_on_Raman_corrections :: Switching on Raman corrections. " << endl;
        cout << " switch_on_Raman_corrections :: nmax will be set to " 
             << Raman_process_nmax << endl;
    }

    return;
}


//==================================================================================================
//
// PDE Data
//
//==================================================================================================
struct PDE_Setup_Data
{
    Cosmos *cosm;
    Gas_of_Atoms *HI;
    PDE_solver_functions *PDE_funcs;
};

static PDE_Setup_Data PDE_Setup_D;

//==================================================================================================
// memory
//==================================================================================================
void init_profile_memory(int nxi, vector<double> &VP)
{
    VP.clear();
    VP.resize(nxi); 
    return;
}

//==================================================================================================
void init_profile_memory(int nxi, Raman_2g_profiles &RGP)
{
    init_profile_memory(nxi, RGP.ratio);
    init_profile_memory(nxi, RGP.profile);
    return;
}

//==================================================================================================
void init_profile_memory(int nxi, int nres, PDE_solver_functions &PDE_funcs)
{
    init_profile_memory(nxi, PDE_funcs.exp_x);
    //
    if(nres>1)
    {
        PDE_funcs.Raman_ns_1s.clear(); PDE_funcs.Raman_ns_1s.resize(nres-1);
        PDE_funcs.Raman_nd_1s.clear(); PDE_funcs.Raman_nd_1s.resize(nres-1);
        //
        PDE_funcs.two_g_ns_1s.clear(); PDE_funcs.two_g_ns_1s.resize(nres-1);
        PDE_funcs.two_g_nd_1s.clear(); PDE_funcs.two_g_nd_1s.resize(nres-1);
        //
        PDE_funcs.Voigt_profiles_A_npns.clear();
        PDE_funcs.Voigt_profiles_A_npnd.clear();
        
        vector<double> dum_v;
        for(int k=0; k<nres-1; k++) // slight waste of memory, but that is ok...
        {
            init_profile_memory(nxi, PDE_funcs.Raman_ns_1s[k]);
            init_profile_memory(nxi, PDE_funcs.two_g_ns_1s[k]);
            init_profile_memory(nxi, PDE_funcs.Raman_nd_1s[k]);
            init_profile_memory(nxi, PDE_funcs.two_g_nd_1s[k]);
            //
            PDE_funcs.Voigt_profiles_A_npns.push_back(dum_v);
            PDE_funcs.Voigt_profiles_A_npnd.push_back(dum_v);
        }
        
        //==============================================================
        // allocate memory for Voigt profiles of Ly-n
        //==============================================================
        PDE_funcs.Voigt_profiles.clear();
        PDE_funcs.HP_Voigt_profiles.clear();
        
        for(int k=0; k<nres; k++)
        {
            PDE_funcs.Voigt_profiles.push_back(PDE_funcs.Raman_ns_1s[0].profile);
            PDE_funcs.HP_Voigt_profiles.push_back(PDE_funcs.Raman_ns_1s[0].profile);
        }
        
        init_profile_memory(nres, PDE_funcs.Voigt_profiles_aV);
        init_profile_memory(nres, PDE_funcs.Voigt_profiles_Dnuj1s_nuj1s);
        init_profile_memory(nres, PDE_funcs.Voigt_profiles_pd);
        init_profile_memory(nres, PDE_funcs.dum_factors);
        init_profile_memory(nres, PDE_funcs.Voigt_profiles_pd_eff);
        init_profile_memory(nres, PDE_funcs.Voigt_profiles_Dnem_eff);
        
        for(int k=0; k<nres-1; k++)
        {
            init_profile_memory(nres, PDE_funcs.Voigt_profiles_A_npns[k]);
            init_profile_memory(nres, PDE_funcs.Voigt_profiles_A_npnd[k]);
        }
    }
    else
    {
        PDE_funcs.two_g_ns_1s.clear(); PDE_funcs.two_g_ns_1s.resize(1);
        init_profile_memory(nxi, PDE_funcs.two_g_ns_1s[0]);
    }
    
    return;
}

//==================================================================================================
// Setup for PDE Data
//==================================================================================================
void Setup_PDE_Data(Cosmos &cosm, Gas_of_Atoms &HIA, PDE_solver_functions &PDE_funcs)
{
    PDE_Setup_D.cosm=&cosm;
    PDE_Setup_D.HI=&HIA;
    PDE_Setup_D.PDE_funcs=&PDE_funcs;
    
    return;
}

//==================================================================================================
// scattering terms
//==================================================================================================
void add_Diff_and_recoil_f(double D_x, double dlnD_x_dx, double xifac, 
                           double &Ai, double &Bi, double &Ci)
{
    Ai+=D_x;
    Bi+=dlnD_x_dx;
    Bi+=D_x*xifac;
    Ci+=dlnD_x_dx*xifac;
    return;
}

void add_Diff_f(double D_x, double dlnD_x_dx, double xifac, double &Ai, double &Bi, double &Ci)
{
    Ai+=D_x;
    Bi+=dlnD_x_dx;
    return;
}

void add_recoil_f(double D_x, double dlnD_x_dx, double xifac, double &Ai, double &Bi, double &Ci)
{
    Bi+=D_x*xifac;
    Ci+=dlnD_x_dx*xifac;
    return;
}

//==================================================================================================
// electron scattering terms
//==================================================================================================
void add_electron_scattering(double z, double nu21, double Te, double Ne, double Hz, 
                             const vector<double> &xi, vector<double> &Ai, 
                             vector<double> &Bi, vector<double> &Ci)
{
    //----------------------------------------------------------
    // xifac= h nu_Lya/ k Te; important for recoil term
    //----------------------------------------------------------
    double xifac=const_h_kb*nu21/Te;
    double kappa_electron_sc=-const_sigT*Ne*(const_kB*Te/const_me_gr/const_cl)/Hz/(1.0+z);
    kappa_electron_sc*=e_sc_efficiency;
    
    void (* add_ptr_f)(double, double, double, double &, double &, double &)=NULL;  
    if(diffusion_e && recoil_e) add_ptr_f=add_Diff_and_recoil_f;
    else if(diffusion_e && !recoil_e) add_ptr_f=add_Diff_f;
    else if(!diffusion_e && recoil_e) add_ptr_f=add_recoil_f;
    
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=0; k<(int)xi.size(); k++) 
    {
        double x=xi[k];
        double D_x=kappa_electron_sc*x*x;
        double dlnD_x_dx=D_x*4.0/x;
        add_ptr_f(D_x, dlnD_x_dx, xifac, Ai[k], Bi[k], Ci[k]);
    }       
    
    return;
}

//==================================================================================================
// Ly-n scattering and em/abs terms
//==================================================================================================
double pd_n_func(int n, double z)
{
    if(n<=10) return PDE_Setup_D.PDE_funcs->HI_pd_nl(z, n, 1);
    else{ cerr << " NOT DEFINED " << endl; exit(0); }
}

double Dnem_n_func(int n, double z)
{
    if(n<=10) return PDE_Setup_D.PDE_funcs->HI_Dnem_nl(z, n, 1);
    else{ cerr << " NOT DEFINED " << endl; exit(0); }
}

//==================================================================================================
// exponential-factor exp(-h\nu/k Tg)
//==================================================================================================
void update_exp_x(double Tg, double nu21, const vector<double> &xi)
{
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=0; k<(int)xi.size(); k++) 
        PDE_Setup_D.PDE_funcs->exp_x[k]=exp(-const_h_kb*nu21*xi[k]/Tg);
    
    return;
}

//==================================================================================================
// Ly-n Voigt profile
//==================================================================================================
void update_Voigt_profile_Ly_n(double z, double nu21, int n, const vector<double> &xi)
{
    double Tm=PDE_Setup_D.PDE_funcs->HI_rho(z)*PDE_Setup_D.cosm->TCMB(z);
    Voigtprofile_Dawson &Lyn=PDE_Setup_D.HI->HI_Lyn_profile(n);
    
    PDE_Setup_D.PDE_funcs->Voigt_profiles_pd[n-2]=pd_n_func(n, z);
    PDE_Setup_D.PDE_funcs->Voigt_profiles_Dnuj1s_nuj1s[n-2]=Lyn.DnuT_nu12(Tm);
    
    double psc_np=1.0-PDE_Setup_D.PDE_funcs->Voigt_profiles_pd[n-2];
    
    if(use_full_width) PDE_Setup_D.PDE_funcs->Voigt_profiles_aV[n-2]=Lyn.aVoigt(Tm)/Lyn.Get_Gamma()
                                                                    *Lyn.Get_A21()/psc_np;
    else PDE_Setup_D.PDE_funcs->Voigt_profiles_aV[n-2]=Lyn.aVoigt(Tm);
        
    double aV=PDE_Setup_D.PDE_funcs->Voigt_profiles_aV[n-2];
        
    //==============================================================================================
    // Ly-n profiles
    //==============================================================================================
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=0; k<(int)xi.size(); k++) 
    {
        double x=xi[k];
        double xD=Lyn.nu2x(nu21*x, Tm);
        PDE_Setup_D.PDE_funcs->Voigt_profiles[n-2][k]=Lyn.phi(xD, aV);
    }
    
    //==============================================================================================
    // corresponding low frequency profiles (only for two-photon profiles)
    //==============================================================================================
    if(two_g_process)
    {
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
        for(int k=0; k<(int)xi.size(); k++) 
        {
            double x=(1.0-1.0/n/n)/0.75-xi[k];
            double xD=Lyn.nu2x(nu21*x, Tm);
            PDE_Setup_D.PDE_funcs->HP_Voigt_profiles[n-2][k]=Lyn.phi(xD, aV);
        }
    }       
    
    return;
}

//==================================================================================================
// access to Lyman-series properties
//==================================================================================================
double phi_Ly_n(int n, int k){ return PDE_Setup_D.PDE_funcs->Voigt_profiles[n-2][k]; }
//
double phi_Ly_n_Dnu_nu(int n){ return PDE_Setup_D.PDE_funcs->Voigt_profiles_Dnuj1s_nuj1s[n-2]; }
double phi_Ly_n_Dnu(int n){ return phi_Ly_n_Dnu_nu(n)*PDE_Setup_D.HI->HI_Lyn_profile(n).Get_nu21();}
//
double phi_Ly_n_aV(int n){ return PDE_Setup_D.PDE_funcs->Voigt_profiles_aV[n-2]; }
double phi_Ly_n_pd(int n){ return PDE_Setup_D.PDE_funcs->Voigt_profiles_pd[n-2]; }
//
double phi_HP_n(int n, int k){ return PDE_Setup_D.PDE_funcs->HP_Voigt_profiles[n-2][k]; }

//==================================================================================================
//
// Raman profiles
//
//==================================================================================================
double Raman_emission_fac(double exp_x, double exp_dum)
{
    double npl_xp=exp_x/(exp_dum-exp_x);
    return npl_xp;
}

//==================================================================================================
void update_Raman_profiles_nsd_1s(double z, int ni, Raman_2g_profiles &RGP, 
                                  vector<double> &Voigt_profiles_A_np_nsd)
{
    double nun1=PDE_Setup_D.HI->Level(ni, 0).Get_Dnu_1s();
    double Tg=PDE_Setup_D.cosm->TCMB(z);
    double x_c=const_h_kb*nun1/Tg;
    double exp_dum=exp(-x_c);
    int klow=PDE_Setup_D.PDE_funcs->index_emission[ni-2];
    int mmax=Get_nmax_two_g_correction()-1; // always just include the resonances up 
                                            // ni_max for Raman-process
    if(Get_nmax_two_g_correction()==0) mmax=PDE_Setup_D.PDE_funcs->Voigt_profiles_aV.size();
        
    //-----------------------------------------------------------------------------
    // setup factors for profile conversion
    //-----------------------------------------------------------------------------
    for(int m=ni-2+1; m<mmax; m++)
        PDE_Setup_D.PDE_funcs->dum_factors[m]=( 1.0-phi_Ly_n_pd(m+2) )
                                                *Voigt_profiles_A_np_nsd[m]/phi_Ly_n_Dnu(m+2);
    
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=klow; k<(int)RGP.profile.size(); k++)
    {
        RGP.profile[k]=0.0;
        
        //-----------------------------------------------------------------------------
        // Ly-ni is not included in ratio
        //-----------------------------------------------------------------------------
        for(int m=ni-2+1; m<mmax; m++) 
            RGP.profile[k]+=phi_Ly_n(m+2, k)*PDE_Setup_D.PDE_funcs->dum_factors[m];
        
        RGP.profile[k]*=RGP.ratio[k];
        RGP.profile[k]*=Raman_emission_fac(PDE_Setup_D.PDE_funcs->exp_x[k], exp_dum);   
    }
    
    return;
}

//==================================================================================================
//
// 2 gamma profiles
//
//==================================================================================================
double two_g_emission_fac(double exp_mx, double exp_mx0)
{
    double npl_x =exp_mx/(1.0-exp_mx);
    double npl_xp=exp_mx0/(exp_mx-exp_mx0);
    return (1.0+npl_xp)*(1.0+npl_x);
}

//==================================================================================================
void update_two_g_profiles_nsd_1s(double z, int ni, Raman_2g_profiles &RGP, 
                                  vector<double> &Voigt_profiles_A_np_nsd)
{
    if(ni<3) return;
    
    double nun1=PDE_Setup_D.HI->Level(ni, 0).Get_Dnu_1s();
    double Tg=PDE_Setup_D.cosm->TCMB(z);
    double x_c=const_h_kb*nun1/Tg;
    double exp_dum=exp(-x_c);
    int kmax=PDE_Setup_D.PDE_funcs->index_emission[ni-2];
    
    //-----------------------------------------------------------------------------
    // setup factors for profile conversion
    //-----------------------------------------------------------------------------
    for(int m=0; m<ni-2; m++)
        PDE_Setup_D.PDE_funcs->dum_factors[m]=( 1.0-phi_Ly_n_pd(m+2) )
                                                *Voigt_profiles_A_np_nsd[m]/phi_Ly_n_Dnu(m+2);
    
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=0; k<kmax; k++)
    {
        RGP.profile[k]=0.0;
        
        //-----------------------------------------------------------------------------
        // Ly-ni is not included in ratio
        //-----------------------------------------------------------------------------
        for(int m=0; m<ni-2; m++)
        {
            double dum=phi_Ly_n(m+2, k);
            //
            if(k>=PDE_Setup_D.PDE_funcs->index_2) dum+=phi_HP_n(m+2, k);
            else dum+=phi_HP_n(m+2, PDE_Setup_D.PDE_funcs->index_2);
            
            RGP.profile[k]+=dum*PDE_Setup_D.PDE_funcs->dum_factors[m];
        }

        RGP.profile[k]*=RGP.ratio[k];
        RGP.profile[k]*=two_g_emission_fac(PDE_Setup_D.PDE_funcs->exp_x[k], exp_dum);   
    }
    
    return;
}

//==================================================================================================
void update_two_g_profiles_2s1s(double z, Raman_2g_profiles &RGP)
{
    double nun1=PDE_Setup_D.HI->Level(2, 0).Get_Dnu_1s();
    double Tg=PDE_Setup_D.cosm->TCMB(z);
    double x_c=const_h_kb*nun1/Tg;
    double exp_dum=exp(-x_c);
    int kmax=PDE_Setup_D.PDE_funcs->index_emission[0];
    
    for(int k=0; k<kmax; k++)
        RGP.profile[k]=RGP.ratio[k]*two_g_emission_fac(PDE_Setup_D.PDE_funcs->exp_x[k],exp_dum)/nun1;
    
    return;
}

//==================================================================================================
//
// access to profiles
//
//==================================================================================================
double phi_ns1s_2g(int n, int k)
{ return PDE_Setup_D.PDE_funcs->two_g_ns_1s[n-2].profile[k]; }

double phi_nd1s_2g(int n, int k)
{ return PDE_Setup_D.PDE_funcs->two_g_nd_1s[n-2].profile[k]; }

double phi_ns1s_Raman(int n, int k)
{ return PDE_Setup_D.PDE_funcs->Raman_ns_1s[n-2].profile[k]; }

double phi_nd1s_Raman(int n, int k)
{ return PDE_Setup_D.PDE_funcs->Raman_nd_1s[n-2].profile[k]; }

//==================================================================================================
//
// ns/nd-1s Raman/2g profile
//
//==================================================================================================
void update_Raman_profiles(double z)
{
    if(Raman_process_nmax>=2)
    { 
        update_Raman_profiles_nsd_1s(z, 2, PDE_Setup_D.PDE_funcs->Raman_ns_1s[0], 
                                     PDE_Setup_D.PDE_funcs->Voigt_profiles_A_npns[0]);
        
        if(Raman_process_nmax>=3)
        { 
            for(int m=1; m<(int)PDE_Setup_D.PDE_funcs->Raman_ns_1s.size()-1; m++)
            {
                update_Raman_profiles_nsd_1s(z, m+2, PDE_Setup_D.PDE_funcs->Raman_ns_1s[m], 
                                             PDE_Setup_D.PDE_funcs->Voigt_profiles_A_npns[m]);
                
                update_Raman_profiles_nsd_1s(z, m+2, PDE_Setup_D.PDE_funcs->Raman_nd_1s[m], 
                                             PDE_Setup_D.PDE_funcs->Voigt_profiles_A_npnd[m]);
            }
        }
    }
    
    return;
}

//==================================================================================================
void update_2gamma_profiles(double z)
{
    // this is always there if one entered the diffusion solver
    update_two_g_profiles_2s1s(z, PDE_Setup_D.PDE_funcs->two_g_ns_1s[0]);
        
    for(int m=1; m<(int)PDE_Setup_D.PDE_funcs->two_g_ns_1s.size(); m++)
    {
        update_two_g_profiles_nsd_1s(z, m+2, PDE_Setup_D.PDE_funcs->two_g_ns_1s[m], 
                                     PDE_Setup_D.PDE_funcs->Voigt_profiles_A_npns[m]);
        
        update_two_g_profiles_nsd_1s(z, m+2, PDE_Setup_D.PDE_funcs->two_g_nd_1s[m], 
                                     PDE_Setup_D.PDE_funcs->Voigt_profiles_A_npnd[m]);
    }
    
    return;
}   

//==================================================================================================
//
// correction to pd and Dnem; depends on two-photon and raman processes
//
//==================================================================================================
void Dp_Dnem_Raman_update(int ni, int li, int n, double N1s, double Nj, double Tg, 
                          double &DRp, double &DRm)
{
    double nuij=PDE_Setup_D.HI->Level(n, 1).Get_nu21(ni, li);
    double exp_x=exp(-const_h_kb*nuij/Tg);
    double npij=exp_x/(1.0-exp_x);
    double Aij=PDE_Setup_D.HI->Level(n, 1).Get_A21(ni, li);
    //
    DRm+=Aij*(1.0+npij);
    DRp+=3.0/(2.0*li+1.0)*Nj*Aij*npij;
    
    return;
}

void Dp_Dnem_2gamma_update(int ni, int li, int n, double N1s, double Nj, double Tg, 
                           double &DRp, double &DRm)
{
    double nuij=PDE_Setup_D.HI->Level(ni, li).Get_nu21(n, 1);
    double exp_x=exp(-const_h_kb*nuij/Tg);
    double npij=exp_x/(1.0-exp_x);
    double Aij=PDE_Setup_D.HI->Level(ni, li).Get_A21(n, 1);
    //
    DRm+=(2.0*li+1.0)/3.0*Aij*npij;
    DRp+=Aij*(1.0+npij)*Nj;
    
    return;
}

//==================================================================================================
//
// combined calls
//
//==================================================================================================
void Dp_Dnem_function(int n, double z, double Tg, double pd_np, double nem, double N1s, 
                      double &Dp, double &Dn)
{
    Dp=Dn=0.0;

    double DRp=0.0, DRm=0.0;
    double NH=PDE_Setup_D.cosm->NH(z);
    
    //----------------------------------------------------------
    // Raman-scatterings
    //----------------------------------------------------------
    if(Raman_process)
    {
        for(int ni=(int)min(7, n-1); ni>=3; ni--)
        {
            Dp_Dnem_Raman_update(ni, 0, n, N1s, PDE_Setup_D.PDE_funcs->HI_Xnl(z, ni, 0)*NH, 
                                 Tg, DRp, DRm);
            
            Dp_Dnem_Raman_update(ni, 2, n, N1s, PDE_Setup_D.PDE_funcs->HI_Xnl(z, ni, 2)*NH, 
                                 Tg, DRp, DRm);
        }   

        if(Raman_process_nmax>=2 && n>2) 
        	Dp_Dnem_Raman_update(2, 0, n, N1s, PDE_Setup_D.PDE_funcs->HI_Xnl(z, 2, 0)*NH, 
                                 Tg, DRp, DRm);
    }
    
    //----------------------------------------------------------
    // 2 gamma events
    //----------------------------------------------------------
    if(two_g_process)
    {
        double DRp2=0.0, DRm2=0.0;
        for(int ni=(int)min(8, two_g_process_nmax); ni>=n+1; ni--)
        {
            Dp_Dnem_2gamma_update(ni, 0, n, N1s, PDE_Setup_D.PDE_funcs->HI_Xnl(z, ni, 0)*NH, 
                                  Tg, DRp2, DRm2);
            
            Dp_Dnem_2gamma_update(ni, 2, n, N1s, PDE_Setup_D.PDE_funcs->HI_Xnl(z, ni, 2)*NH, 
                                  Tg, DRp2, DRm2);
        }   
        
        DRp+=DRp2; DRm+=DRm2;
    }
    
    //----------------------------------------------------------
    // return if no correction is made
    //----------------------------------------------------------
    if(DRm==0.0) return;
    
    //----------------------------------------------------------
    // corrected death probability and Dnem
    //----------------------------------------------------------
    double Rm_tot=PDE_Setup_D.HI->Level(n, 1).Get_A21(1, 0)/(1.0-pd_np); // R-,tot = A + R-
    //double Rm    =Rm_tot-PDE_Setup_D.HI->Level(n, 1).Get_A21(1, 0); // JC mod 13th Aug 2012
    double Rm    =PDE_Setup_D.HI->Level(n, 1).Get_A21(1, 0)*pd_np/(1.0-pd_np);
    //
    Dp=DRm/Rm_tot;
    //
    double Rp=3.0*nem*N1s*Rm;
    //double nem_t=(Rp-DRp)/(Rm-DRm)/3.0/N1s;                         // JC mod 13th Aug 2012
    //Dn=nem_t-nem;                                                   // JC mod 13th Aug 2012
    Dn=( (Rp-DRp) - 3.0*nem*N1s*(Rm-DRm) )/(Rm-DRm)/3.0/N1s;

    return;
}

//==================================================================================================
//
// PDE setup functions
//
//==================================================================================================
void add_Ly_n_sc_only(int n, double z, double nu21, double Tg, double Te, double Ne, double N1s, 
                      double Hz, const vector<double> &xi, 
                      vector<double> &Ai, vector<double> &Bi, vector<double> &Ci)
{
    //----------------------------------------------------------
    // xifac= h nu_Lya/ k Te; important for recoil term
    //----------------------------------------------------------
    double xifac=const_h_kb*nu21/Te;
    double nun1=PDE_Setup_D.HI->Level(n, 1).Get_Dnu_1s();
    
    //----------------------------------------------------------
    // Ly-n scattering
    //----------------------------------------------------------
    double lambda=const_cl/nu21; // lambdaj1*(nuj1/nu21)
    double aV=phi_Ly_n_aV(n);
    double DnuD_nu=phi_Ly_n_Dnu_nu(n);
    //
    double pd_np=phi_Ly_n_pd(n);
    double psc_np=1.0-pd_np;
    double const_sigr=1.5*lambda*lambda*(PDE_Setup_D.HI->HI_Lyn_profile(n).Get_A21()
                                         /FOURPI/phi_Ly_n_Dnu(n));
    //
    double kappa_Ly_n_sc =-psc_np*const_sigr*N1s*(const_kB*Te/const_mH_gr/const_cl)/Hz/(1.0+z);
    
    void (* add_ptr_f)(double, double, double, double &, double &, double &)=NULL;  
    if(diffusion_l && recoil_l) add_ptr_f=add_Diff_and_recoil_f;
    else if(diffusion_l && !recoil_l) add_ptr_f=add_Diff_f;
    else if(!diffusion_l && recoil_l) add_ptr_f=add_recoil_f;
    
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=0; k<(int)xi.size(); k++) 
    {
        double x=xi[k];
        
        //------------------------------------------------------
        // line-scattering terms 
        //------------------------------------------------------
        double xD   =PDE_Setup_D.HI->HI_Lyn_profile(n).nu2x(nu21*x, Te);
        double phi_x=phi_Ly_n(n, k);
        double dphi_x_dx=PDE_Setup_D.HI->HI_Lyn_profile(n).dphi_dx(xD, aV) / DnuD_nu * nu21/nun1;
        double D_x=kappa_Ly_n_sc*phi_x;
        double dlnD_x_dx=kappa_Ly_n_sc*(phi_x*2.0/x+dphi_x_dx);
        
        add_ptr_f(D_x, dlnD_x_dx, xifac, Ai[k], Bi[k], Ci[k]);
    }
    
    return;
}

//==================================================================================================
void add_Ly_n_em_only(int n, double z, double nu21, double Tg, double Te, double Ne, double N1s, 
                      double Hz, const vector<double> &xi, vector<double> &Ci, vector<double> &Di)
{   
    //----------------------------------------------------------
    // Ly-n scattering
    //----------------------------------------------------------
    double lambda=const_cl/nu21; // lambdaj1*(nuj1/nu21)
    double pd_np=phi_Ly_n_pd(n);
    double const_sigr=1.5*lambda*lambda*(PDE_Setup_D.HI->HI_Lyn_profile(n).Get_A21()
                                         /FOURPI/phi_Ly_n_Dnu(n));
    
    //----------------------------------------------------------
    // Ly-n emission/absorption term
    //----------------------------------------------------------
    double nun1=PDE_Setup_D.HI->Level(n, 1).Get_Dnu_1s();
    double x_c=const_h_kb*nun1/Tg;
    double Dnem =Dnem_n_func(n, z)*nu21;
    //
    double kappa_Ly_n_abs=-pd_np*const_sigr*N1s*const_cl/Hz/(1.0+z); 
    double exp_dum=exp(-x_c);

    //----------------------------------------------------------
    // correcting the death probability and Dnem
    //----------------------------------------------------------
    if(Raman_process || two_g_process) 
    {
        double Dp, Dn;
        Dp_Dnem_function(n, z, Tg, pd_np, Dnem/nu21+exp_dum, N1s, Dp, Dn);
        //
        kappa_Ly_n_abs*=(1.0-Dp/pd_np);
        pd_np*=(1.0-Dp/pd_np);
        Dnem+=Dn*nu21;
    }

    //----------------------------------------------------------
    // save info for computation of DPesc
    //----------------------------------------------------------
    PDE_Setup_D.PDE_funcs->Voigt_profiles_pd_eff[n-2]=pd_np;
    PDE_Setup_D.PDE_funcs->Voigt_profiles_Dnem_eff[n-2]=Dnem/nu21;
    
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=0; k<(int)xi.size(); k++) 
    {
        double x=xi[k];
        
        //------------------------------------------------------
        // prepare for line-scattering and em/abs terms 
        //------------------------------------------------------
        double phi_x=phi_Ly_n(n, k);

        //------------------------------------------------------
        // Doppler core emission/absorption; Lorentzian wings 
        // will not be included in the 1+1 terms when switched on
        //------------------------------------------------------
        //double xD=PDE_Setup_D.HI->HI_Lyn_profile(n).nu2x(nu21*x, Te); 
        //double phi_x=PDE_Setup_D.HI->HI_Lyn_profile(n).phi_DC(xD);
        
        //------------------------------------------------------
        // Source-term from emission & absorption 
        //------------------------------------------------------
        double exp_fac=exp_dum/PDE_Setup_D.PDE_funcs->exp_x[k];
        double D_x=kappa_Ly_n_abs*phi_x/x/x;
        
        //------------------------------------------------------
        Ci[k]+=-D_x*exp_fac;
        Di[k]+= D_x*Dnem;
        //------------------------------------------------------
    }
    
    return;
}

//==================================================================================================
void add_Ly_n_sc_em_full(int n, double z, double nu21, double Tg, double Te, double Ne, double N1s, 
                         double Hz, const vector<double> &xi, 
                         vector<double> &Ai, vector<double> &Bi, 
                         vector<double> &Ci, vector<double> &Di)
{
    //----------------------------------------------------------
    // xifac= h nu_Lya/ k Te; important for recoil term
    //----------------------------------------------------------
    double xifac=const_h_kb*nu21/Te;
    
    //----------------------------------------------------------
    // Ly-n scattering
    //----------------------------------------------------------
    double pd_np=phi_Ly_n_pd(n);
    double psc_np=1.0-pd_np;
    //
    double lambda=const_cl/nu21; // lambdaj1*(nuj1/nu21)
    double aV=phi_Ly_n_aV(n);
    double DnuD_nu=phi_Ly_n_Dnu_nu(n);
    //
    double const_sigr=1.5*lambda*lambda*(PDE_Setup_D.HI->HI_Lyn_profile(n).Get_A21()
                                         /FOURPI/phi_Ly_n_Dnu(n));
    
    double kappa_Ly_n_sc =-psc_np*const_sigr*N1s*(const_kB*Te/const_mH_gr/const_cl)/Hz/(1.0+z);
    
    void (* add_ptr_f)(double, double, double, double &, double &, double &)=NULL;  
    if(diffusion_l && recoil_l) add_ptr_f=add_Diff_and_recoil_f;
    else if(diffusion_l && !recoil_l) add_ptr_f=add_Diff_f;
    else if(!diffusion_l && recoil_l) add_ptr_f=add_recoil_f;
    
    //----------------------------------------------------------
    // Ly-n emission/absorption term
    //----------------------------------------------------------
    double nun1=PDE_Setup_D.HI->Level(n, 1).Get_Dnu_1s();
    double x_c=const_h_kb*nun1/Tg;
    double Dnem =Dnem_n_func(n, z)*nu21;
    //
    double kappa_Ly_n_abs=-pd_np*const_sigr*N1s*const_cl/Hz/(1.0+z); 
    double exp_dum=exp(-x_c);
    
    //----------------------------------------------------------
    // correcting the death probability and Dnem
    //----------------------------------------------------------
    if(Raman_process || two_g_process) 
    {
        double Dp, Dn;
        Dp_Dnem_function(n, z, Tg, pd_np, Dnem/nu21+exp_dum, N1s, Dp, Dn);
        //
        kappa_Ly_n_abs*=1.0-Dp/pd_np;
        pd_np*=1.0-Dp/pd_np;
        Dnem+=Dn*nu21;
    }
    
    //----------------------------------------------------------
    // save info for computation of DPesc
    //----------------------------------------------------------
    PDE_Setup_D.PDE_funcs->Voigt_profiles_pd_eff[n-2]=pd_np;
    PDE_Setup_D.PDE_funcs->Voigt_profiles_Dnem_eff[n-2]=Dnem/nu21;
    
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=0; k<(int)xi.size(); k++) 
    {
        double x=xi[k];
        
        //------------------------------------------------------
        // prepare for line-scattering and em/abs terms 
        //------------------------------------------------------
        double phi_x=phi_Ly_n(n, k);
        
        //------------------------------------------------------
        // line-scattering terms 
        //------------------------------------------------------
        double xD   =PDE_Setup_D.HI->HI_Lyn_profile(n).nu2x(nu21*x, Te);
        double dphi_x_dx=PDE_Setup_D.HI->HI_Lyn_profile(n).dphi_dx(xD, aV) / DnuD_nu * nu21/nun1;
        double D_x=kappa_Ly_n_sc*phi_x;
        double dlnD_x_dx=kappa_Ly_n_sc*(phi_x*2.0/x+dphi_x_dx);
        
        add_ptr_f(D_x, dlnD_x_dx, xifac, Ai[k], Bi[k], Ci[k]);
        
        //------------------------------------------------------
        // Source-term from emission & absorption 
        //------------------------------------------------------
        double exp_fac=exp_dum/PDE_Setup_D.PDE_funcs->exp_x[k];
        
        //------------------------------------------------------
        // Doppler core emission/absorption; Lorentzian wings 
        // will not be included in the 1+1 terms when switched on
        //------------------------------------------------------
        //if(n==2) phi_x=PDE_Setup_D.HI->HI_Lyn_profile(n).phi_DC(xD); 
        //phi_x=0.0;
        
        D_x=kappa_Ly_n_abs*phi_x/x/x;
        
        //------------------------------------------------------
        Ci[k]+=-D_x*exp_fac;
        Di[k]+= D_x*Dnem;
        //------------------------------------------------------
    }
    
    return;
}

//==================================================================================================
// 
// add emission, absorption and scattering from nD-1s quadrupole lines (added June 2011)
//
//==================================================================================================
void add_nD_1s_sc_em_full(int n, double z, double nu21, double Tg, double Te, double Ne, double N1s, 
                          double Hz, const vector<double> &xi, 
                          vector<double> &Ai, vector<double> &Bi, 
                          vector<double> &Ci, vector<double> &Di)
{
    Voigtprofile_Dawson *p=&PDE_Setup_D.HI->HI_nD1s_profile(n);

    //----------------------------------------------------------
    // nD-1s lines; scattering can be neglected
    //----------------------------------------------------------
    double lambda=const_cl/nu21; // lambdaj1*(nuj1/nu21)
    double aV=p->aVoigt(Te);
    double nun1=p->Get_nu21();
    double DnuD_nu=p->DnuT_nu12(Te);
    double DnuD=DnuD_nu*nun1;
    //
    double const_sigr=2.5*lambda*lambda*p->Get_A21()/FOURPI/DnuD;
    
    //----------------------------------------------------------
    // nD-1s emission/absorption term
    //----------------------------------------------------------
    double x_c=const_h_kb*nun1/Tg;
    double Dnem =PDE_Setup_D.PDE_funcs->HI_Dnem_nl(z, n, 2)*nu21;
    //
    double kappa_nD_1s_abs=-const_sigr*N1s*const_cl/Hz/(1.0+z); 
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
        double xD=p->nu2x(nu21*x, Te);
        double phi_x=p->phi(xD, aV);	
        
        //------------------------------------------------------
        // Source-term from emission & absorption 
        //------------------------------------------------------
        double exp_fac=exp_dum/PDE_Setup_D.PDE_funcs->exp_x[k];
        double D_x=kappa_nD_1s_abs*phi_x/x/x;
        
        //------------------------------------------------------
        Ci[k]+=-D_x*exp_fac;
        Di[k]+= D_x*Dnem;
        //------------------------------------------------------
    }
    
    return;
}

//==================================================================================================
//
// Raman-scattering terms
//
//==================================================================================================
void add_nsd_1s_Raman_terms(int n, int l, double (*phi_func)(int, int), double Dnem,
                            double z, double nu21, double Tg, double N1s, double Hz, 
                            const vector<double> &xi, vector<double> &Ci, vector<double> &Di)
{
    //----------------------------------------------------------
    // ns/nd-1s Raman term
    //----------------------------------------------------------    
    double nun1=PDE_Setup_D.HI->Level(n, 1).Get_Dnu_1s();
    //
    double sig_nsd=(2.0*l+1.0)/2.0/FOURPI*pow(const_cl/nu21, 2); 
    double kappa_nsd1s=-sig_nsd*N1s*const_cl/Hz/(1.0+z);
    //
    double Dnem_nsd=Dnem*nu21;
    double x_c=const_h_kb*nun1/Tg;
    double exp_dum=exp(-x_c);
    //
    int klow=PDE_Setup_D.PDE_funcs->index_emission[n-2];    

#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=klow; k<(int)xi.size(); k++) 
    {
        double x=xi[k];
        double phi_nsd_Raman_em=phi_func(n, k); // vacuum profile versus dy=dnu/nuj1
        //
        double npl_x =PDE_Setup_D.PDE_funcs->exp_x[k]/(1.0-PDE_Setup_D.PDE_funcs->exp_x[k]);
        //  
        double exp_fac_nsd=exp_dum/npl_x;
        double DD=kappa_nsd1s*phi_nsd_Raman_em/x/x;
        
        //------------------------------------------------------
        Ci[k]+=-DD*exp_fac_nsd;
        Di[k]+= DD*Dnem_nsd;
        //------------------------------------------------------
    }
    
    return;
}

//==================================================================================================
//
// 2gamma-terms
//
//==================================================================================================
void add_nsd_1s_2gamma_terms(int n, int l, double (*phi_func)(int, int), double Dnem,
                             double z, double nu21, double Tg, double N1s, double Hz, 
                             const vector<double> &xi, vector<double> &Ci, vector<double> &Di)
{
    //----------------------------------------------------------
    // ns/nd-1s 2 gamma term
    //----------------------------------------------------------    
    double nun1=PDE_Setup_D.HI->Level(n, l).Get_Dnu_1s();
    //
    double sig_nsd=(2.0*l+1.0)/2.0/FOURPI*pow(const_cl/nu21, 2); 
    double kappa_nsd1s=-sig_nsd*N1s*const_cl/Hz/(1.0+z);
    //
    double Dnem_nsd=Dnem*nu21;
    double x_c=const_h_kb*nun1/Tg;
    double exp_dum=exp(-x_c);
    //
    int kmax=PDE_Setup_D.PDE_funcs->index_emission[n-2];    
    
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=0; k<kmax; k++) 
    {
        double x=xi[k];
        double phi_nsd_2g_em=phi_func(n, k); // vacuum profile versus dy=dnu/nuj1
        //
        double npl_x=PDE_Setup_D.PDE_funcs->exp_x[k]/(1.0-PDE_Setup_D.PDE_funcs->exp_x[k]);
        //
        double exp_fac_nsd=exp_dum/npl_x;
        double DD=kappa_nsd1s*phi_nsd_2g_em/x/x;
        
        //------------------------------------------------------
        Ci[k]+=-DD*exp_fac_nsd;
        Di[k]+= DD*Dnem_nsd;
        //------------------------------------------------------
    }
    
    return;
}

//==================================================================================================
//
// simple output of profile functions
//
//==================================================================================================
void plot_Raman(const vector<double> &xi)
{
    if(Get_nmax_Raman_correction()<2) return;
    
    wait_f_r(Get_nmax_Raman_correction());

    string fname="./temp/Raman.dat";
    ofstream ofile(fname.c_str());
    ofile.precision(10);
    
    for(int k=0; k<(int)xi.size(); k++) 
    {
        ofile << xi[k] << " " << xi[k]-1.0 << " ";
        if(Get_nmax_Raman_correction()>=2) ofile << phi_ns1s_Raman(2, k) << " ";
        for(int m=1; m<=Get_nmax_Raman_correction()-2; m++) 
            ofile << phi_ns1s_Raman(m+2, k) << " " << phi_nd1s_Raman(m+2, k) << " ";        
        ofile << endl;  
    }
    
    ofile.close();
    
    wait_f_r(" Raman ");
    
    return;
}   

void plot_two_g(const vector<double> &xi, int nresmax)
{
    if(Get_nmax_two_g_correction()<2) return;

    wait_f_r(Get_nmax_two_g_correction());
    
    string fname="./temp/two_g.dat";
    ofstream ofile(fname.c_str());
    ofile.precision(10);
    
    for(int k=0; k<(int)xi.size(); k++) 
    {
        ofile << xi[k] << " " << phi_ns1s_2g(2, k) << " ";
        for(int m=1; m<=nresmax-2; m++) 
            ofile << phi_ns1s_2g(m+2, k) << " " << phi_nd1s_2g(m+2, k) << " ";
        ofile << endl;  
    }
    
    ofile.close();
    
    wait_f_r(" two-g ");
    
    return;
}   
        

//==================================================================================================
//
// main wrapper
//
//==================================================================================================
void def_PDE_Lyn_and_2s1s(double z, const vector<double> &xi, 
                          vector<double> &Ai, vector<double> &Bi, 
                          vector<double> &Ci, vector<double> &Di)
{
    int nmax=PDE_Setup_D.PDE_funcs->Voigt_profiles_aV.size();
    
    //----------------------------------------------------------
    // xi = nu/nu_Lya
    //----------------------------------------------------------
    double nu21=PDE_Setup_D.HI->Level(2, 1).Get_Dnu_1s();
    
    //----------------------------------------------------------
    // time-dependent variables 
    //----------------------------------------------------------
    double Hz=PDE_Setup_D.cosm->H(z);
    double NH=PDE_Setup_D.cosm->NH(z);
    double Tg=PDE_Setup_D.cosm->TCMB(z);
    
    //----------------------------------------------------------
    // solutions from rec-code
    //----------------------------------------------------------
    double Ne=NH*PDE_Setup_D.PDE_funcs->HI_Xe(z);
    double N1s=NH*PDE_Setup_D.PDE_funcs->HI_X1s(z);
    double rho=PDE_Setup_D.PDE_funcs->HI_rho(z);
    double Te=Tg*rho;
    
    //----------------------------------------------------------
    // exp(x) is at least needed in calling routine
    //----------------------------------------------------------
    update_exp_x(Tg, nu21, xi);
    
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
    if(def_HI_run.el_sc) add_electron_scattering(z, nu21, Te, Ne, Hz, xi, Ai, Bi, Ci);

    //----------------------------------------------------------
    // normal Ly-n scatterring and emission/absorption
    //----------------------------------------------------------
    if(def_HI_run.line_sc || def_HI_run.line_em_abs)
    {
        bool local_line_em_abs=def_HI_run.line_em_abs;
        bool local_line_sc=def_HI_run.line_sc;
        
        for(int n=2; n<=nmax; n++)
        {
            update_Voigt_profile_Ly_n(z, nu21, n, xi);
            
            if(local_line_sc && local_line_em_abs) 
                add_Ly_n_sc_em_full(n, z, nu21, Tg, Te, Ne, N1s, Hz, xi, Ai, Bi, Ci, Di);
            
            else if(local_line_sc && !local_line_em_abs) 
                add_Ly_n_sc_only(n, z, nu21, Tg, Te, Ne, N1s, Hz, xi, Ai, Bi, Ci);
            
            else if(!local_line_sc && local_line_em_abs) 
                add_Ly_n_em_only(n, z, nu21, Tg, Te, Ne, N1s, Hz, xi, Ci, Di);
        }
    }
    
    //----------------------------------------------------------
    // Raman-scattering
    //----------------------------------------------------------
    update_Raman_profiles(z);
    
    if(Raman_process)
    {
        if(Raman_process_nmax>=2) 
            add_nsd_1s_Raman_terms(2, 0, phi_ns1s_Raman, 
                                   PDE_Setup_D.PDE_funcs->HI_Dnem_nl(z, 2, 0), 
                                   z, nu21, Tg, N1s, Hz, xi, Ci, Di);

        for(int ni=3; ni<=(int)min(7, Raman_process_nmax); ni++)
        {
            add_nsd_1s_Raman_terms(ni, 0, phi_ns1s_Raman, 
                                   PDE_Setup_D.PDE_funcs->HI_Dnem_nl(z, ni, 0), 
                                   z, nu21, Tg, N1s, Hz, xi, Ci, Di);
            
            add_nsd_1s_Raman_terms(ni, 2, phi_nd1s_Raman, 
                                   PDE_Setup_D.PDE_funcs->HI_Dnem_nl(z, ni, 2), 
                                   z, nu21, Tg, N1s, Hz, xi, Ci, Di);
        }   
    }
    
//  if(z<=1100.0) plot_Raman(xi);

    //----------------------------------------------------------
    // two-photon emission/absorption
    //----------------------------------------------------------
    update_2gamma_profiles(z);
    
    if(two_g_process)
    {
        for(int ni=3; ni<=(int)min(8, two_g_process_nmax); ni++)
        {
            add_nsd_1s_2gamma_terms(ni, 0, phi_ns1s_2g, 
                                    PDE_Setup_D.PDE_funcs->HI_Dnem_nl(z, ni, 0), 
                                    z, nu21, Tg, N1s, Hz, xi, Ci, Di);
            
            add_nsd_1s_2gamma_terms(ni, 2, phi_nd1s_2g, 
                                    PDE_Setup_D.PDE_funcs->HI_Dnem_nl(z, ni, 2), 
                                    z, nu21, Tg, N1s, Hz, xi, Ci, Di);
        }
    }
    
//  if(z<=1100.0) plot_two_g(xi, nmax);
    
    //----------------------------------------------------------
    // Source-term 2s-1s emission and absorption
    //----------------------------------------------------------
    if(def_HI_run.A2s1s_em_abs) 
        add_nsd_1s_2gamma_terms(2, 0, phi_ns1s_2g, PDE_Setup_D.PDE_funcs->HI_Dnem_nl(z, 2, 0), 
                                z, nu21, Tg, N1s, Hz, xi, Ci, Di);
    
    //----------------------------------------------------------
    // include nD-1s quadrupole lines (added June 2011)
    //----------------------------------------------------------
    for(int n=3; def_HI_run.nD_Quadrupole && n<=(int)min(8, two_g_process_nmax); n++)
        add_nD_1s_sc_em_full(n, z, nu21, Tg, Te, Ne, N1s, Hz, xi, Ai, Bi, Ci, Di);

    return;
}


//==================================================================================================
//
// routines to account for corrections to electron scattering using Kernel
//
//==================================================================================================
void electron_kernel_corr_func_HI(double z, const vector<double> &xi, 
                                  const vector<double> &yi, 
                                  vector<double> &corr_v)
{
    //----------------------------------------------------------
    // xi = nu/nu_Lya
    //----------------------------------------------------------
    double nu21=PDE_Setup_D.HI->Level(2, 1).Get_Dnu_1s();
    
    //----------------------------------------------------------
    // time-dependent variables 
    //----------------------------------------------------------
    double Hz=PDE_Setup_D.cosm->H(z);
    double NH=PDE_Setup_D.cosm->NH(z);
    double Tg=PDE_Setup_D.cosm->TCMB(z);
    
    //----------------------------------------------------------
    // solutions from rec-code
    //----------------------------------------------------------
    double Ne=NH*PDE_Setup_D.PDE_funcs->HI_Xe(z);
    double rho=PDE_Setup_D.PDE_funcs->HI_rho(z);
    double Te=Tg*rho;    
    double kap_e=-const_sigT*Ne*(const_kB*Te/const_me_gr/const_cl)/Hz/(1.0+z);
    
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=1; k<(int)xi.size()-1; k++) 
        corr_v[k]=kap_e*compute_kernel_Integral_polint(k, nu21, Te, xi, yi);
    
    return;
}

//==================================================================================================
//==================================================================================================
