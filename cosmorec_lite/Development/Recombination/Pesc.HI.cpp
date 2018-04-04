//============================================================================
// Author: Jens Chluba
// July 2007
// Purpose: compute the correction to the escape probability in the helium 
// lines caused by the absorption from hydrogen
//============================================================================
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>

#include "Pesc.HI.h"
#include "Photoionization_cross_section.h"
#include "Voigtprofiles.h"
#include "routines.h"
#include "Integration_routines.h"

using namespace std;

//============================================================================
// global data
//============================================================================
struct Var_Pesc_SH_Integral
{
    double Te;
    double eta_S;
    double eta_c;
    double a;
    double xD;
    double xDm;
    Voigtprofile_Dawson *P;
    Photoionization_cross_section_SH *L;
};

Var_Pesc_SH_Integral Pesc_SH_Data_I;

//============================================================================
// inner integral performed analytically (Rubino-Martin et al, 2008)
//============================================================================
double Inner_Int_appr(double Te, double eta_S, double eta_c, double xD, double a, double chimax,
                      Voigtprofile_Dawson &P, Photoionization_cross_section_SH &L)
{
    double phi=P.phi(xD, a);
    double nu=P.x2nu(xD, Te);
    double eta_cont=eta_c*L.sig_phot_ion(nu)/L.sig_phot_ion_nuc()*P.DnuT(Te)/nu/phi;
    double xi=eta_S/(eta_S+eta_cont);
    double ym=chimax-P.xi_Int(xD, a);
    double dtau_c=(eta_S+eta_cont)*ym; 
    double dtau_S=eta_S*ym;
    
    return phi*( 1.0-exp(-dtau_S)-xi*( 1.0-exp(-dtau_c) ));
}

//============================================================================
// the outer integral using symmetry
//============================================================================
double dPesc_appr_I_sym(double xD, void *p)
{ 
  return Inner_Int_appr(Pesc_SH_Data_I.Te, Pesc_SH_Data_I.eta_S, Pesc_SH_Data_I.eta_c, 
                        xD, Pesc_SH_Data_I.a, Pesc_SH_Data_I.xDm, *Pesc_SH_Data_I.P, *Pesc_SH_Data_I.L)
        +Inner_Int_appr(Pesc_SH_Data_I.Te, Pesc_SH_Data_I.eta_S, Pesc_SH_Data_I.eta_c, 
                        -xD, Pesc_SH_Data_I.a, Pesc_SH_Data_I.xDm, *Pesc_SH_Data_I.P, *Pesc_SH_Data_I.L); 
}

double DPesc_appr_I_sym(double Te, double eta_S, double eta_c, Voigtprofile_Dawson &P, 
                        Photoionization_cross_section_SH &L, double abs)
{
    Pesc_SH_Data_I.Te=Te;
    Pesc_SH_Data_I.eta_S=eta_S;
    Pesc_SH_Data_I.eta_c=eta_c;
    Pesc_SH_Data_I.P=&P;
    Pesc_SH_Data_I.L=&L;
    Pesc_SH_Data_I.a=P.aVoigt(Te);
    
    double xDl=0.0, xDu=1.0e+4;
    Pesc_SH_Data_I.xDm=P.xi_Int(xDu, Pesc_SH_Data_I.a);
    
    double r=0.0, r1;  
    double epsabs=abs, epsrel=1.0e-6;
    double a=xDl, b=xDu;
    void *params=NULL;
    
    a=xDl; b=xDl+1.0; r1=0.0;
    r1=Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dPesc_appr_I_sym, params);
    r+=r1;
//  cout << " A: " << r << endl;
    
    epsabs=fabs(r)*epsrel;
    a=b; b=xDl+2.0; r1=0.0;
    r1=Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dPesc_appr_I_sym, params);
    r+=r1;
//  cout << " B: " << r << endl;
    
    epsabs=fabs(r)*epsrel;
    a=b; b=xDl+4.0; r1=0.0;
    r1=Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dPesc_appr_I_sym, params);
    r+=r1;
//  cout << " C: " << r << endl;
    
    epsabs=fabs(r)*epsrel;
    a=b; b=xDl+10.0; r1=0.0;
    r1=Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dPesc_appr_I_sym, params);
    r+=r1;
//  cout << " D: " << r << endl;
    
    epsabs=fabs(r)*epsrel;
    a=b; b=xDl+30.0; r1=0.0; 
    r1=Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dPesc_appr_I_sym, params);
    r+=r1;
//  cout << " E: " << r << endl;

    epsabs=fabs(r)*epsrel;
    a=b; b=xDl+100.0; r1=0.0; 
    r1=Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dPesc_appr_I_sym, params);
    r+=r1;
//  cout << " F: " << r << endl;    
    
    epsabs=fabs(r)*epsrel;
    a=b; b=xDu; r1=0.0; 
    r1=Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dPesc_appr_I_sym, params);
    r+=r1;
//  cout << " G: " << r << endl;
    
    return r;
}

//============================================================================
// the outer integral using symmetry & log-coordinates
//============================================================================
double dPesc_appr_I_sym_lg(double lgxD, void *p)
{  
  double xD=exp(lgxD);
  return xD*Inner_Int_appr(Pesc_SH_Data_I.Te, Pesc_SH_Data_I.eta_S, Pesc_SH_Data_I.eta_c, 
               -xD, Pesc_SH_Data_I.a, Pesc_SH_Data_I.xDm, *Pesc_SH_Data_I.P, *Pesc_SH_Data_I.L)
        +xD*Inner_Int_appr(Pesc_SH_Data_I.Te, Pesc_SH_Data_I.eta_S, Pesc_SH_Data_I.eta_c, 
                xD, Pesc_SH_Data_I.a, Pesc_SH_Data_I.xDm, *Pesc_SH_Data_I.P, *Pesc_SH_Data_I.L); 
}

double DPesc_appr_I_sym_lg(double Te, double eta_S, double eta_c, Voigtprofile_Dawson &P, 
                           Photoionization_cross_section_SH &L, double abs)
{
    Pesc_SH_Data_I.Te=Te;
    Pesc_SH_Data_I.eta_S=eta_S;
    Pesc_SH_Data_I.eta_c=eta_c;
    Pesc_SH_Data_I.P=&P;
    Pesc_SH_Data_I.L=&L;
    Pesc_SH_Data_I.a=P.aVoigt(Te);
    
    double xDu=1.0e+4;
    Pesc_SH_Data_I.xDm=P.xi_Int(xDu, Pesc_SH_Data_I.a);
    
    double r=0.0;  
    double epsabs=abs, epsInt=1.0e-6;
    double a, b;
    
    a=log(1.0e-20); b=log(xDu);
    void *p=NULL;
    r=Integrate_using_Patterson_adaptive(a, b, epsInt, epsabs, dPesc_appr_I_sym_lg, p);
    
    return r;
}

