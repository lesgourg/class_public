//======================================================================================
// Author: Jens Chluba 
// last modification: Oct 2010
// 
// Purpose: compute the recombination & photoionzation rates of a hydrogenic atom in a
// blackbody ambient radiation field of temperature Tg. The electrons have temperature Te.
//======================================================================================

//======================================================================================
// class based on Storey & Hummer, 1991 
// these routines have been checked for n<=300 but should work until n~1000
//======================================================================================

#include <iostream>
#include <string>
#include <cmath>

#include "Photoionization_cross_section.h"
#include "Rec_Phot_BB.h"
#include "routines.h"
#include "physical_consts.h"
#include "Recombination_Integrals.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

using namespace std;

//======================================================================================
// parameters for the NAG integration routines
// reference model: double NAG_Rec_Phot_BB_epsInt=1.0e-10;
//======================================================================================
double NAG_Rec_Phot_BB_epsInt=2.0e-9;
double NAG_Rec_Phot_BB_absInt=1.0e-100;
double NAG_Rec_Phot_BB_absInt_fac=0.0e-10;
//
double NAG_xi_max_fac_SH=1.0e+8;      
double NAG_xi_max_fac_SH_QSP=1.0e+8;
// 
int NAG_np_spline=1024;
double NAG_xi_max_spline=5.0e+7;      // for qubic spline interpolant
double NAG_xi_max_spline_split=30.0;  // split point for spline. Half of the spline points will be there
double NAG_xi_max_limit=1.0e+8;       // 250 worked // limit on nu/nu0. Beyond gnl=0
//
double Rec_Phot_BB_rho;
double Rec_Phot_BB_T_g;
double Rec_Phot_BB_xi;
static Rec_Phot_BB_SH *ptr_to_Rec_Phot_BB_SH;

//======================================================================================
//
// Konstructors and Destructors
//
//======================================================================================
void Rec_Phot_BB_SH::init(int nv, int lv, int Z, double Np, int mflag)
{
    Photoionization_cross_section_SH::init(nv, lv, Z, Np, mflag);
    
    g_l=2.0*(2.0*lv+1.0);                           // statistical weight; 2 from e-spin 
    if(nv<=2) accRec=2;
    else accRec=0;

    xi_small=nu_sig_phot_ion_small(1.0e-30)/Get_nu_ionization();
    
    return;
}

Rec_Phot_BB_SH::Rec_Phot_BB_SH():Photoionization_cross_section_SH(){ }

Rec_Phot_BB_SH::Rec_Phot_BB_SH(int nv, int lv, int Z, double Np, int mflag):Photoionization_cross_section_SH()
{ init(nv, lv, Z, Np, mflag); }

double Rec_Phot_BB_SH::sig_phot_ion_lim(double nu)
{ 
  if(nu<Get_nu_ionization() || nu>=NAG_xi_max_limit*Get_nu_ionization()) return 0.0;
  return sig_phot_ion(nu); 
}

//======================================================================================
// integrations
//======================================================================================
double Rec_Phot_BB_SH::R_nl_c_Int(double T_g){return R_nl_c_JC(T_g);}                
//
double Rec_Phot_BB_SH::R_c_nl_Int(double T_g, double rho){return R_c_nl_JC(T_g, rho);}
//
double Rec_Phot_BB_SH::dR_c_nl_dTe_Int(double T_g, double rho){return dR_c_nl_dTe_JC(T_g, rho);}

//======================================================================================
// functions for integration routine 
//======================================================================================
double Rec_Phot_BB_SH_dR_nl_c(double x)
{ return x*x/(1.0-exp(-x))*exp(Rec_Phot_BB_xi-x)*ptr_to_Rec_Phot_BB_SH->sig_phot_ion_lim(x*Rec_Phot_BB_T_g/const_h_kb); }

double Rec_Phot_BB_SH_dR_nl_c_log(double lgx)
{   
    double x=exp(lgx);
    return x*Rec_Phot_BB_SH_dR_nl_c(x); 
}

double Rec_Phot_BB_SH_dR_c_nl(double x)
{ 
    // rho=TM/Tg
    double rho=Rec_Phot_BB_rho;    
    double xi=Rec_Phot_BB_xi;
    
    return x*x/(1.0-exp(-x))*exp((xi-x)/rho)*ptr_to_Rec_Phot_BB_SH->sig_phot_ion_lim(x*Rec_Phot_BB_T_g/const_h_kb); 
}

double Rec_Phot_BB_SH_dR_c_nl_log(double lgx)
{ 
    double x=exp(lgx);
    return x*Rec_Phot_BB_SH_dR_c_nl(x);
}

//======================================================================================
double Rec_Phot_BB_SH::R_nl_c_JC(double T_g)
{
    double xi=nu2xx(Get_nu_ionization(), T_g);
    // 08.10.2009: do not compute photonionization rate, when things are FAR TOO small
    if(xi>=400) return 1.0e-300;
    
    ptr_to_Rec_Phot_BB_SH=this;
    Rec_Phot_BB_T_g=T_g; 
    Rec_Phot_BB_xi=xi; 
    
    double r=0.0;
    double epsInt=NAG_Rec_Phot_BB_epsInt, epsabs=0.0;
    double ab=xi;
    double bb=xi*NAG_xi_max_fac_SH;
    //=============================================
    // exponential cut-off due to BB
    //=============================================
    if(xi>=10.0) bb=min(bb, xi-log(1.0e-30));
    if(xi< 10.0) bb=min(bb, xi+10.0-log(1.0e-30));
    //=============================================
    // is cross-section falling faster or planck?
    //=============================================
    bb=min(bb, xi*xi_small);
    
    r=Integrate_Ric_Rci_JC_Chebyshev_log(1, log(ab), log(bb), epsInt, epsabs, Rec_Phot_BB_SH_dR_nl_c_log);
    
    return 2.0*FOURPI*const_cl*pow(T_g/const_hcl_kb, 3)*r*exp(-xi);    
}

//======================================================================================
double Rec_Phot_BB_SH::R_c_nl_JC(double T_g, double rho)
{
    double xi=nu2xx(Get_nu_ionization(), T_g);
    //=============================================
    // rho=TM/Tg
    // changed from x to xe=x/rho
    //=============================================
    ptr_to_Rec_Phot_BB_SH=this;
    Rec_Phot_BB_rho=rho; 
    Rec_Phot_BB_T_g=T_g; 
    Rec_Phot_BB_xi=xi;
    
    double r=0.0;
    double epsInt=NAG_Rec_Phot_BB_epsInt, epsabs=0.0;
    double ab=xi;
    double bb=xi*NAG_xi_max_fac_SH;
    //--------------------------------------------------------------
    // recombination cross section exponentially cuts of for x-xi >> rho
    //--------------------------------------------------------------
    if(xi>=10.0) bb=min(bb, xi-rho*log(1.0e-30));
    if(xi< 10.0) bb=min(bb, xi+10.0-rho*log(1.0e-30));
    //=============================================
    // is cross-section falling faster or planck?
    //=============================================
    bb=min(bb, xi*xi_small);

    r=Integrate_Ric_Rci_JC_Patterson(2, log(ab), log(bb), epsInt, epsabs, Rec_Phot_BB_SH_dR_c_nl_log);
    
    return 2.0*FOURPI*const_cl*(g_l/2.0)*r*pow(const_kb_mec2/Get_mu_red()*T_g/rho/TWOPI, 1.5);
}

//======================================================================================
// Computing the derivatives with respect to Te explicitly 
//======================================================================================
double Rec_Phot_BB_SH_ddR_c_nl_dTe(double x)
{ 
    //=============================================
    // rho=TM/Tg
    //=============================================
    double rho=Rec_Phot_BB_rho;    
    double xi=Rec_Phot_BB_xi;
    
    return (x-xi-1.5*rho)*x*x/(1.0-exp(-x))*exp((xi-x)/rho)*ptr_to_Rec_Phot_BB_SH->sig_phot_ion_lim(x*Rec_Phot_BB_T_g/const_h_kb); 
}

double Rec_Phot_BB_SH_ddR_c_nl_dTe_log(double lgx)
{ 
    double x=exp(lgx);
    return x*Rec_Phot_BB_SH_ddR_c_nl_dTe(x);
}

//======================================================================================
double Rec_Phot_BB_SH::dR_c_nl_dTe_JC(double T_g, double rho)
{
    double xi=nu2xx(Get_nu_ionization(), T_g);
    //=============================================
    // rho=TM/Tg
    //=============================================
    ptr_to_Rec_Phot_BB_SH=this;
    Rec_Phot_BB_rho=rho; 
    Rec_Phot_BB_T_g=T_g; 
    Rec_Phot_BB_xi=xi;
    
    double r=0.0;
    double epsInt=NAG_Rec_Phot_BB_epsInt, epsabs=0.0;
    double ab=xi;
    double bb=xi*NAG_xi_max_fac_SH;
    //--------------------------------------------------------------
    // recombination cross section exponentially cuts of for x-xi >> rho
    //--------------------------------------------------------------
    if(xi>=10.0) bb=min(bb, xi-rho*log(1.0e-30));
    if(xi< 10.0) bb=min(bb, xi+10.0-rho*log(1.0e-30));
    //=============================================
    // is cross-section falling faster or planck?
    //=============================================
    bb=min(bb, xi*xi_small);
    
    r=Integrate_Ric_Rci_JC_Chebyshev_log(3, log(ab), log(bb), epsInt, epsabs, Rec_Phot_BB_SH_ddR_c_nl_dTe_log, log(xi+1.5*rho));

    return 2.0*FOURPI*const_cl*(g_l/2.0)*r*pow(const_kb_mec2/Get_mu_red()*T_g/rho/TWOPI, 1.5)/(rho*rho*T_g);
}


//======================================================================================
//
// Konstructors and Destructors
// class with qubic splines for Gaunt-factors
//
//======================================================================================
void Rec_Phot_BB_SH_QSP::init(int nv, int lv, int Z, double Np, int mflag)
{
    Photoionization_cross_section_SH::init(nv, lv, Z, Np, mflag);
    
    g_l=2.0*(2.0*lv+1.0);                           // statistical weight; 2 from e-spin 
    if(nv<=2) accRec=2;
    else accRec=0;
    
    nxi=NAG_np_spline;
    xi_small=nu_sig_phot_ion_small(1.0e-30)/Photoionization_cross_section_SH::Get_nu_ionization();  
    
    // don't set up any splines if this is not improving the performance
    if(lv>=nv-10)
    {
        spline_is_allocated=0;
    }
    else
    {
        acc=gsl_interp_accel_alloc();
        spline=gsl_spline_alloc(gsl_interp_cspline, nxi);
        spline_is_allocated=1;
        
        calc_coeff_for_spline();
    }
    
    return;
}

//======================================================================================
Rec_Phot_BB_SH_QSP::Rec_Phot_BB_SH_QSP():Photoionization_cross_section_SH()
{ spline_is_allocated=0; }

Rec_Phot_BB_SH_QSP::Rec_Phot_BB_SH_QSP(int nv, int lv, int Z, double Np, int mflag):Photoionization_cross_section_SH()
{
    spline_is_allocated=0;
    init(nv, lv, Z, Np, mflag);
}

Rec_Phot_BB_SH_QSP::~Rec_Phot_BB_SH_QSP()
{
    if(spline_is_allocated==1) 
    {   
        gsl_spline_free(spline);
        gsl_interp_accel_free(acc); 
    }       
}

//======================================================================================
// GSL interpolation
//======================================================================================
int Rec_Phot_BB_SH_QSP::calc_coeff_for_spline()
{
    double *xi=new double[nxi];
    double *lgxi=new double[nxi];
    double *gaunt=new double[nxi];
    
    if(mess_flag>=2)
        cout << " Rec_Phot_BB_SH_QSP::calc_coeff_for_spline:" 
        << "\n calculating coefficients for Gaunt-factor interpolation for level (n,l) " 
        << nn << " " << ll << endl;
    
    //==================================================================================
    // fill xi array in log-scale but in two pieces if necessary
    //==================================================================================
    if(min(xi_small, NAG_xi_max_spline)<=NAG_xi_max_spline_split) init_xarr(1.0, min(xi_small, NAG_xi_max_spline)*1.001, xi, nxi, 1, 0);
    else 
    {
        int npl=nxi/2, nup=nxi-npl; 
        //
        init_xarr(1.0, NAG_xi_max_spline_split, xi, npl, 1, 0); 
        //
        init_xarr(NAG_xi_max_spline_split, min(xi_small, NAG_xi_max_spline)*1.001, &xi[npl-1], nup+1, 1, 0);
    }
    
    //==================================================================================
    // fill gaunt array
    //==================================================================================
    double nu0=Get_nu_ionization();
    for(int i=0; i<nxi; i++){ lgxi[i]=log(xi[i]); gaunt[i]=log(Photoionization_cross_section_SH::g_phot_ion(xi[i]*nu0)); }
    if(mess_flag>=2) for(int i=0; i<nxi; i++) cout << xi[i]*nu0 << " " << gaunt[i] << " " << Photoionization_cross_section_SH::g_phot_ion(xi[i]*nu0) << endl;
    
    gsl_spline_init(spline, lgxi, gaunt, nxi);
    
    if(mess_flag>=2) cout << " done..." << endl << endl;
    
    //==================================================================================
    // clean up
    //==================================================================================
    delete[] xi;
    delete[] lgxi;
    delete[] gaunt;
    return 0;
}

double Rec_Phot_BB_SH_QSP::calc_gaunt_fac_spline(double xiv)
{
    //--------------------------------------------------------------------------------------
    // for very high shells one has to limit the maximal value of nu/nu0 to avoid divergence
    //--------------------------------------------------------------------------------------
    if(xiv<1.0 || xiv>=NAG_xi_max_limit) return 0.0;
    if(xiv>=min(xi_small, NAG_xi_max_spline)) return Photoionization_cross_section_SH::g_phot_ion(xiv*Get_nu_ionization());
    if(xiv<=1.0+1.0e-10) return Photoionization_cross_section_SH::gaunt_nuc();

    if(spline_is_allocated==1) return exp(gsl_spline_eval(spline, log(xiv), acc));

    return Photoionization_cross_section_SH::g_phot_ion(xiv*Get_nu_ionization());
}

double Rec_Phot_BB_SH_QSP::sig_phot_ion_lim(double nu){ return sig_phot_ion_Kramers(nu)*calc_gaunt_fac_spline(nu/Get_nu_ionization()); }
double Rec_Phot_BB_SH_QSP::g_phot_ion_lim(double nu){ return calc_gaunt_fac_spline(nu/Get_nu_ionization()); }

//======================================================================================
// integrations
//======================================================================================
double Rec_Phot_BB_SH_QSP::R_nl_c_Int(double T_g){return R_nl_c_JC(T_g);}                
//
double Rec_Phot_BB_SH_QSP::R_c_nl_Int(double T_g, double rho){return R_c_nl_JC(T_g, rho);}
//
double Rec_Phot_BB_SH_QSP::dR_c_nl_dTe_Int(double T_g, double rho){return dR_c_nl_dTe_JC(T_g, rho);}

//======================================================================================
// functions for integration routine 
//======================================================================================
struct Ric_vars_SH_QSP
{
    double T_g;
    double xi;
    Rec_Phot_BB_SH_QSP *SH_QSP_ptr;
};

double Rec_Phot_BB_SH_QSP_dR_nl_c_P(double x, void *p)
{ 
    Ric_vars_SH_QSP *V=((Ric_vars_SH_QSP *) p); 
    
    return x*x/(1.0-exp(-x))*exp(V->xi-x)*V->SH_QSP_ptr->sig_phot_ion_lim(x*V->T_g/const_h_kb); 
}

double Rec_Phot_BB_SH_QSP_dR_nl_c_log_P(double lgx, void *p)
{   
    double x=exp(lgx);
    return x*Rec_Phot_BB_SH_QSP_dR_nl_c_P(x, p); 
}

//======================================================================================
struct Rci_vars_SH_QSP
{
    double rho;
    double T_g;
    double xi;
    Rec_Phot_BB_SH_QSP *SH_QSP_ptr;
};

double Rec_Phot_BB_SH_QSP_dR_c_nl_P(double x, void *p)
{ 
    // rho=TM/Tg
    Rci_vars_SH_QSP *V=((Rci_vars_SH_QSP *) p); 
    
    return x*x/(1.0-exp(-x))*exp((V->xi-x)/V->rho)*V->SH_QSP_ptr->sig_phot_ion_lim(x*V->T_g/const_h_kb); 
}

double Rec_Phot_BB_SH_QSP_dR_c_nl_log_P(double lgx, void *p)
{ 
    double x=exp(lgx);
    return x*Rec_Phot_BB_SH_QSP_dR_c_nl_P(x, p);
}

//======================================================================================
double Rec_Phot_BB_SH_QSP::R_nl_c_JC(double T_g)
{
    double xi=nu2xx(Get_nu_ionization(), T_g);
    // 08.10.2009: do not compute photonionization rate, when things are FAR TOO small
    if(xi>=400) return 1.0e-300;
    //
    Ric_vars_SH_QSP Vars;
    //
    Vars.T_g=T_g; 
    Vars.xi=xi;
    Vars.SH_QSP_ptr=this;   
    void *params=&Vars;
    
    double r=0.0;
    double epsInt=NAG_Rec_Phot_BB_epsInt, epsabs=0.0;
    double ab=xi;
    double bb=xi*NAG_xi_max_fac_SH_QSP;
    //=============================================
    // exponential cut-off due to BB
    //=============================================
    if(xi>=10.0) bb=min(bb, xi-log(1.0e-30));
    if(xi< 10.0) bb=min(bb, xi+10.0-log(1.0e-30));
    //=============================================
    // is cross-section falling faster or planck?
    //=============================================
    bb=min(bb, xi*xi_small);
    
    r=Integrate_Ric_Rci_JC_Patterson(1, log(ab), log(bb), epsInt, epsabs, Rec_Phot_BB_SH_QSP_dR_nl_c_log_P, params);
    
    return 2.0*FOURPI*const_cl*pow(T_g/const_hcl_kb, 3)*r*exp(-xi);    
}

//======================================================================================
double Rec_Phot_BB_SH_QSP::R_c_nl_JC(double T_g, double rho)
{
    double xi=nu2xx(Get_nu_ionization(), T_g);
    //=============================================
    // rho=TM/Tg
    //=============================================
    Rci_vars_SH_QSP Vars;
    //
    Vars.rho=rho;
    Vars.T_g=T_g; 
    Vars.xi=xi;
    Vars.SH_QSP_ptr=this;   
    void *params=&Vars;
    
    double r=0.0;
    double epsInt=NAG_Rec_Phot_BB_epsInt, epsabs=0.0;
    double ab=xi;
    double bb=xi*NAG_xi_max_fac_SH_QSP;
    //--------------------------------------------------------------
    // recombination cross section exponentially cuts of for x-xi >> rho
    //--------------------------------------------------------------
    if(xi>=10.0) bb=min(bb, xi-rho*log(1.0e-30));
    if(xi< 10.0) bb=min(bb, xi+10.0-rho*log(1.0e-30));
    //=============================================
    // is cross-section falling faster or planck?
    //=============================================
    bb=min(bb, xi*xi_small);
    
    r=Integrate_Ric_Rci_JC_Patterson(2, log(ab), log(bb), epsInt, epsabs, Rec_Phot_BB_SH_QSP_dR_c_nl_log_P, params);
    
    return 2.0*FOURPI*const_cl*(g_l/2.0)*r*pow(const_kb_mec2/Get_mu_red()*T_g/rho/TWOPI, 1.5);
}

//======================================================================================
// Computing the derivatives with respect to Te explicitly 
//======================================================================================
struct dRci_dT_vars_SH_QSP
{
    double rho;
    double T_g;
    double xi;
    Rec_Phot_BB_SH_QSP *SH_QSP_ptr;
};

double Rec_Phot_BB_SH_QSP_ddR_c_nl_dTe_P(double x, void *p)
{ 
    //=============================================
    // rho=TM/Tg
    //=============================================
    dRci_dT_vars_SH_QSP *V=((dRci_dT_vars_SH_QSP *) p); 
    
    return (x-V->xi-1.5*V->rho)*x*x/(1.0-exp(-x))*exp((V->xi-x)/V->rho)*V->SH_QSP_ptr->sig_phot_ion_lim(x*V->T_g/const_h_kb); 
}

double Rec_Phot_BB_SH_QSP_ddR_c_nl_dTe_log_P(double lgx, void *p)
{ 
    double x=exp(lgx);
    return x*Rec_Phot_BB_SH_QSP_ddR_c_nl_dTe_P(x, p); 
}

//======================================================================================
double Rec_Phot_BB_SH_QSP::dR_c_nl_dTe_JC(double T_g, double rho)
{
    double xi=nu2xx(Get_nu_ionization(), T_g);
    //=============================================
    // rho=TM/Tg
    //=============================================
    dRci_dT_vars_SH_QSP Vars;
    //
    Vars.rho=rho;
    Vars.T_g=T_g; 
    Vars.xi=xi;
    Vars.SH_QSP_ptr=this;   
    void *params=&Vars;
    
    double r=0.0;
    double epsInt=NAG_Rec_Phot_BB_epsInt, epsabs=0.0;
    double ab=xi;
    double bb=xi*NAG_xi_max_fac_SH_QSP;
    //--------------------------------------------------------------
    // recombination cross section exponentially cuts of for x-xi >> rho
    //--------------------------------------------------------------
    if(xi>=10.0) bb=min(bb, xi-rho*log(1.0e-30));
    if(xi< 10.0) bb=min(bb, xi+10.0-rho*log(1.0e-30));
    //=============================================
    // is cross-section falling faster or planck?
    //=============================================
    bb=min(bb, xi*xi_small);
    
    r=Integrate_Ric_Rci_JC_Patterson(3, log(ab), log(bb), epsInt, epsabs, Rec_Phot_BB_SH_QSP_ddR_c_nl_dTe_log_P, params);
    
    return 2.0*FOURPI*const_cl*(g_l/2.0)*r*pow(const_kb_mec2/Get_mu_red()*T_g/rho/TWOPI, 1.5)/(rho*rho*T_g);
}

