//====================================================================================================================
// Author Jens Chluba March 2011
//====================================================================================================================

//====================================================================================================================
// Purpose: Implementation of electron scattering kernel functions following
// Sazonov & Sunyaev, 2000, Apj 543, p. 28-55.
//====================================================================================================================

//--------------------------------------------------------------------------------------------------------------------
// Th_e == k Te / me c^2
//--------------------------------------------------------------------------------------------------------------------

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>

#include "Compton_Kernel.h"
#include "physical_consts.h"
#include "routines.h"

using namespace std;

//====================================================================================================================
//
// definitions according to Sazonov & Sunyaev, 2000, Apj 543, p. 28-55
//
//====================================================================================================================
inline double Kernel_p0(double delta)
{
    double F=exp(-delta*delta);
    double G=0.5*SQRT_PI*erfc_JC(fabs(delta));    
    double d2=delta*delta;
    return (0.55+(0.8+0.4*d2)*d2)*F-fabs(delta)*(1.5+(2.0+0.8*d2)*d2)*G;
}

//====================================================================================================================
//
// table of the kernel function p0 to avoid calling erfc(|delta|) and exp(-delta^2)
//
//====================================================================================================================
bool p0_Kernel_spline_was_used=0;
int p0_Kernel_spline_memindex_plus;
int p0_Kernel_spline_memindex_minus;
double p0_Kernel_spline_dmax=sqrt(-log(1.0e-8));

//====================================================================================================================
// interpolation function setup and call
//====================================================================================================================
double Kernel_p0_spline(double delta)
{
    if(!p0_Kernel_spline_was_used)
    {
        p0_Kernel_spline_was_used=1;

        int np=500;
        vector <double> xa(np), ya(np), xar(np);
        xa[0]=0.0;
        init_xarr(1.0e-8, p0_Kernel_spline_dmax, &xa[1], np-1, 1, 0);
        
        for(int k=0; k<np; k++)
        {
            double d=xa[k];
            ya[k]=log(Kernel_p0(d));
        }
        
        p0_Kernel_spline_memindex_plus=calc_spline_coeffies_JC(np, &xa[0], &ya[0], "Kernel_p0_spline :: + ");

        for(int k=0; k<np; k++)
        {
            xar[k]=-xa[np-1-k];
            double d=xar[k];
            ya[k]=log(Kernel_p0(d));
        }
        
        p0_Kernel_spline_memindex_minus=calc_spline_coeffies_JC(np, &xar[0], &ya[0], "Kernel_p0_spline :: - ");
    }
    
    if(delta>p0_Kernel_spline_dmax || delta<-p0_Kernel_spline_dmax) return Kernel_p0(delta);
    
    if(delta>=0) return exp(calc_spline_JC(delta, p0_Kernel_spline_memindex_plus));
    else return exp(calc_spline_JC(delta, p0_Kernel_spline_memindex_minus));

    return 0.0;
}

//====================================================================================================================
//
// lowest order kernel; no recoil, no Doppler boosting;
//
//====================================================================================================================
double P0_Kernel(double nu, double nup, double Th_e)
{
    double delta=Kernel_delta(nu, nup, Th_e);
    double p_0=Kernel_p0(delta);
    return sqrt(2.0/Th_e)/SQRT_PI*p_0/nu;
}

double P0_Kernel_spline(double nu, double nup, double Th_e)
{
    double delta=Kernel_delta(nu, nup, Th_e);
    double p_0=Kernel_p0_spline(delta);
    return sqrt(2.0/Th_e)/SQRT_PI*p_0/nu;
}

double norm_P0(double nu, double Th_e)
{ return 1.0+1.5*Th_e; }

double P0_Kernel_AliHaimoud(double nu, double nup, double Th_e)
{
    double P0=P0_Kernel(nu, nup, Th_e);
    double del=const_h_mec2*(nup-nu)/Th_e;
    double fac=2.0/( 1.0+pow(nu/nup, 2)*exp(del) );
    return P0*fac;
}

//====================================================================================================================
//
// 'Kompaneets' kernel; Doppler broadening, Doppler boosting, and electron recoil terms are included
//
//====================================================================================================================
double PK_Kernel(double nu, double nup, double Th_e)
{
    double delta=Kernel_delta(nu, nup, Th_e);
    double p_0=Kernel_p0(delta);
    return sqrt(2.0/Th_e)/SQRT_PI*p_0/nu*(1.0+sqrt(2.0*Th_e)*delta*(1.0-const_h_mec2*nu/Th_e));
}

double PK_Kernel_spline(double nu, double nup, double Th_e)
{
    double delta=Kernel_delta(nu, nup, Th_e);
    double p_0=Kernel_p0_spline(delta);
    return sqrt(2.0/Th_e)/SQRT_PI*p_0/nu*(1.0+sqrt(2.0*Th_e)*delta*(1.0-const_h_mec2*nu/Th_e));
}

double norm_PK(double nu, double Th_e)
{ return 1.0+2.5*Th_e-const_h_mec2*nu; }


//====================================================================================================================
//
// 'full' kernel according to Eq. (19) of Sazonov & Sunyaev, 2000
//
//====================================================================================================================
double P_Compton(double nu, double nup, double Th_e)
{
    double delta=Kernel_delta(nu, nup, Th_e), absd=fabs(delta);
    double omega=const_h_mec2*nu;
    double d2=delta*delta;
    double xe=omega/Th_e;
    
    double F=exp(-delta*delta);
    double G=0.5*SQRT_PI*erfc_JC(fabs(delta));   
    
    double p0=(0.55+(0.8+0.4*d2)*d2)*F-absd*(1.5+(2.0+0.8*d2)*d2)*G;
    double pt=( 
               (-1091.0/1120.0+(-507.0/560.0+(57.0/35.0+68.0/35.0*d2)*d2)*d2)*F
                         +absd*(9.0/4.0+(1.0-(26.0/5.0+136.0/35.0*d2)*d2)*d2)*G
              )*Th_e;
    
    double pr=( 
               (-23.0/280.0+(26.0/35.0+(34.0/35.0+16.0/35.0*d2)*d2)*d2)*F
                                  -absd*(2.0+(2.4+32.0/35.0*d2)*d2)*d2 *G
              )*xe*omega;
    
    double sqrt_T=sqrt(Th_e);
    //============================================================================
    // here the analytic renormalization, Eq. (26), was used
    //============================================================================
    double norm=SQRT_2/sqrt_T/SQRT_PI/nu/(1.0-2.0*omega+(-5.3-44.0/5.0*xe+63.0/20*xe*xe)*Th_e*Th_e);
    //double norm=sqrt_2/(sqrt_PI*sqrt_T*T_e*xe)*(1.0+(5.3+3.8*xe+2.05*xe*xe)*T_e*T_e);
    
    return norm*( (1.0+delta*(SQRT_2*(1.0-xe)*sqrt_T
                         -4.0*delta*omega
                         +2.0*SQRT_2*d2*(-2.0+xe*xe/3.0)*omega*sqrt_T)
                         )*p0
                 +(1.0+delta* SQRT_2*(1.0-xe)*sqrt_T)*pt
                 +(1.0+delta* SQRT_2*(3.0-xe)*sqrt_T)*pr
                 );
}

//====================================================================================================================
//====================================================================================================================
