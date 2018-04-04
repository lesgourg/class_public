//======================================================================
// Standards
//======================================================================
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <iomanip>
#include <cmath>
#include <limits.h>
#include <vector>

//======================================================================
// some definitions and useful stuff
//======================================================================
#include "physical_consts.h"
#include "routines.h"

//======================================================================
// Hydrogen & Helium Atom
//======================================================================
#include "Atom.h"
#include "HeI_Atom.h"

//======================================================================
// local header
//======================================================================
#include "ODE_effective.h"

using namespace std;


//======================================================================
// beginning of Namespace "ODE_HI_effective_funcs"                    //
//======================================================================

namespace ODE_HI_effective_funcs 
{

    //==================================================================
    // constants in the computation
    //==================================================================
    const double rho_g_fac=PI*PI/15.0*const_kB*pow(const_kB/(const_hbar*const_cl), 3)/pow(const_cl, 2)/const_me_gr;

    //==================================================================
    // simple functions
    //==================================================================
    double rho_g_1_cm3(double Tg){return ODE_HI_effective_funcs::rho_g_fac*pow(Tg, 4);}
    //
    double tau_S_Ly_n(const double &Aij, const double &lij, const double &N1s, 
                      const double &Nnp, const double &H_z, const double &w)
                     { return Aij*pow(lij, 3)/(8.0*PI*H_z)*(N1s*w - Nnp); }

    double dtau_S_Ly_n_dX1s(const double &Aij, const double &lij, const double &NH, 
                            const double &H_z, const double &w)
                           { return  Aij*pow(lij, 3)/(8.0*PI*H_z)*3.0*NH; }

    double dtau_S_Ly_n_dXnp(const double &Aij, const double &lij, const double &NH, 
                            const double &H_z)
                           { return -Aij*pow(lij, 3)/(8.0*PI*H_z)*NH;     }
    //
    double p_ij(const double &tau_S){ return (1.0-exp(-tau_S))/tau_S; }             
    double dp_ij_dtauS(const double &tau_S){ return (exp(-tau_S)*(1.0+tau_S)-1.0)/tau_S/tau_S; }                
}

//======================================================================
// End of Namespace "ODE_HI_effective_funcs"                          //
//======================================================================



//======================================================================
// beginning of Namespace "ODE_effective"                             //
//======================================================================

namespace ODE_effective 
{
    
    //==================================================================
    //
    // electron temperature equation (rho=Te/Tg)
    //
    //==================================================================
    void evaluate_TM(double z, double Xe, double fHe, double rho, double Tg, double Hz, double &drho_dt)
    {
        //--------------------------------------------------------------
        // Compton cooling
        //--------------------------------------------------------------
        drho_dt=8.0/3.0*const_sigT*const_cl*ODE_HI_effective_funcs::rho_g_1_cm3(Tg)*Xe/(1.0+Xe+fHe)*(1.0-rho);
        
        //--------------------------------------------------------------
        // adiabatic expansion
        //--------------------------------------------------------------
        drho_dt+=-Hz*rho;
        
        return;
    }
    
    //==================================================================
    void df_drho_evaluate_TM(double z, double Xe, double fHe, double Tg, double Hz, double &drho_dt)
    {
        drho_dt =-8.0/3.0*const_sigT*const_cl*ODE_HI_effective_funcs::rho_g_1_cm3(Tg)*Xe/(1.0+Xe+fHe);
        drho_dt+=-Hz;       
        return;
    }
    
    //==================================================================
    void df_dXHI1s_evaluate_TM(double z, double Xe, double fHe, double rho, 
                               double Tg, double Hz, double &drho_dt)
    {
        drho_dt=-8.0/3.0*const_sigT*const_cl*ODE_HI_effective_funcs::rho_g_1_cm3(Tg)
                        *(1.0+fHe)/(1.0+Xe+fHe)/(1.0+Xe+fHe)*(1.0-rho);
        return;
    }
    
    //==================================================================
    void df_dXHeI1s_evaluate_TM(double z, double Xe, double fHe, double rho, 
                                double Tg, double Hz, double &drho_dt)
    {
        df_dXHI1s_evaluate_TM(z, Xe, fHe, rho, Tg, Hz, drho_dt);
        return;
    }
    
    
    //==================================================================
    //
    // 2s - 1s two-photon decay rate
    //
    //==================================================================
    double Lambda_ind(double zz)
    {
        return 1.0;
        
        //==============================================================
        // this is the correction to the 2s-1s decay rate caused by 
        // stimulated emission (Chluba & Sunyaev 2006). Since this is 
        // included in the diffusion code this functions is not used 
        // normally, but can be activated optionally.
        //==============================================================
        double aa1=1.181e-3, aa2=2.177e-2, aa3=-1.958e-3;
        //double aa1=8.97397e-4, aa2=2.23429e-2, aa3=-2.20249e-3;
        double zhigh=2500.0, zlow=500.0, xi;
        
        if(zz > zhigh) xi=zhigh/1100.0;
        else if(zlow > zz) xi=0.0;
        else xi=zz/1100.0;
        
        return 1.0+(aa1+(aa2+aa3*xi)*xi)*xi;
    }

    //==================================================================
    void evaluate_2s_two_photon_decay(double z, double Tg, double X1s, double X2s, double nu21, 
                                      double A2s1s, double &dX1s_H_dt, double &dX2s_H_dt) 
    {
        double dum=A2s1s*Lambda_ind(z)*( X2s-X1s*exp(-const_h_kb*nu21/Tg) );
        // 1s and 2s contributions
        dX1s_H_dt+= dum;
        dX2s_H_dt+=-dum;
        
        return;
    }
    
    //==================================================================
    void df_dX1s_evaluate_2s_two_photon_decay(double z, double Tg, double nu21, double A2s1s, 
                                              double &dX1s_H_dt, double &dX2s_H_dt) 
    {
        double dum=-A2s1s*Lambda_ind(z)*exp(-const_h_kb*nu21/Tg);
        // 1s and 2s contributions
        dX1s_H_dt+= dum;
        dX2s_H_dt+=-dum;
        
        return;
    }
    
    //==================================================================
    void df_dX2s_evaluate_2s_two_photon_decay(double z, double Tg, double nu21, double A2s1s, 
                                              double &dX1s_H_dt, double &dX2s_H_dt) 
    {
        double dum=A2s1s*Lambda_ind(z);
        // 1s and 2s contributions
        dX1s_H_dt+= dum;
        dX2s_H_dt+=-dum;
        
        return;
    }
    
    //==================================================================
    //
    // nP - 1s channel 
    // 
    // 28.12.2010: added weight w as option for quadrupole lines; for 
    // normal Ly-n w=3 for quadrupole lines w=5
    //
    //==================================================================
    void evaluate_Ly_n_channel(double z, double Tg, double X1s, double Xnp, double NH, 
                               double H_z, double A21, double lambda21, double nu21, 
                               double &dX1s_H_dt, double &dXnp_H_dt, double w) 
    {
        double tau_S=ODE_HI_effective_funcs::tau_S_Ly_n(A21, lambda21, X1s*NH, Xnp*NH, H_z, w);
        double pij=ODE_HI_effective_funcs::p_ij(tau_S);
        double efac=exp(-const_h_kb*nu21/Tg);
        double dum=pij*A21/(1.0-efac)*( Xnp - w*X1s*efac ); 
        // 1s and nP contributions      
        dX1s_H_dt+= dum;
        dXnp_H_dt+=-dum;
        
        return;
    }
    
    void df_dX1s_evaluate_Ly_n_channel(double z, double Tg, double X1s, double Xnp, double NH, 
                                       double H_z, double A21, double lambda21, double nu21, 
                                       double &dX1s_H_dt, double &dXnp_H_dt, double w)
    {
        double tau_S=ODE_HI_effective_funcs::tau_S_Ly_n(A21, lambda21, X1s*NH, Xnp*NH, H_z, w);
        double dtau_S_dx=ODE_HI_effective_funcs::dtau_S_Ly_n_dX1s(A21, lambda21, NH, H_z, w);
        //
        double pij=ODE_HI_effective_funcs::p_ij(tau_S);
        double dpij_dtau=ODE_HI_effective_funcs::dp_ij_dtauS(tau_S);
        //
        double efac=exp(-const_h_kb*nu21/Tg);
        double dum=A21/(1.0-efac)*( ( Xnp - w*X1s*efac )*dpij_dtau*dtau_S_dx-pij*w*efac ); 
        // 1s and nP contributions      
        dX1s_H_dt+= dum;
        dXnp_H_dt+=-dum;
        
        return;
    }
    
    void df_dXnp_evaluate_Ly_n_channel(double z, double Tg, double X1s, double Xnp, double NH, 
                                       double H_z, double A21, double lambda21, double nu21, 
                                       double &dX1s_H_dt, double &dXnp_H_dt, double w)
    {
        double tau_S=ODE_HI_effective_funcs::tau_S_Ly_n(A21, lambda21, X1s*NH, Xnp*NH, H_z, w);
        double dtau_S_dx=ODE_HI_effective_funcs::dtau_S_Ly_n_dXnp(A21, lambda21, NH, H_z);
        //
        double pij=ODE_HI_effective_funcs::p_ij(tau_S);
        double dpij_dtau=ODE_HI_effective_funcs::dp_ij_dtauS(tau_S);
        //
        double efac=exp(-const_h_kb*nu21/Tg);
        double dum=A21/(1.0-efac)*( ( Xnp - w*X1s*efac )*dpij_dtau*dtau_S_dx+pij ); 
        // 1s and nP contributions      
        dX1s_H_dt+= dum;
        dXnp_H_dt+=-dum;
        
        return;
    }

    //==================================================================
    //
    // nD - 1s quadrupole channel (added May 2011)
    //
    //==================================================================
    void evaluate_nD_1s_Q_channel(double Xnp, double Xnd, double A21, 
                                  double &dXnp_H_dt, 
                                  double &dXnd_H_dt)
    {
//        return;
        double DR=A21*(Xnd-5.0/3.0*Xnp);
        
        dXnp_H_dt+= DR;
        dXnd_H_dt+=-DR;
        
        return;
    }
    
    void df_dXnp_evaluate_nD_1s_Q_channel(double Xnd, double A21, 
                                          double &dXnp_H_dt, 
                                          double &dXnd_H_dt)
    {
//        return;
        double DR=-A21*5.0/3.0;
        
        dXnp_H_dt+= DR;
        dXnd_H_dt+=-DR;
        
        return;
    }
    
    void df_dXnd_evaluate_nD_1s_Q_channel(double Xnp, double A21, 
                                          double &dXnp_H_dt, 
                                          double &dXnd_H_dt)
    {
//        return;
        double DR=A21;
        
        dXnp_H_dt+= DR;
        dXnd_H_dt+=-DR;
        
        return;
    }    
        
}   

//======================================================================
// End of Namespace "ODE_effective"                                   //
//======================================================================



//======================================================================
// beginning of Namespace "ODE_HI_effective"                          //
//======================================================================

namespace ODE_HI_effective {

    //==================================================================
    //
    // effective recombination and photonionization rates
    //
    //==================================================================
    void evaluate_effective_Rci_Ric_terms(double z, double Xe, double Np, Gas_of_Atoms &HIA, 
                                          vector<double> &Ai, vector<double> &Bi, 
                                          int (*get_HI_index)(int), double *dXi_H_dt)
    {
        int index_HI;
        for(int m=0; m<(int) Ai.size(); m++)
        {
            index_HI=get_HI_index(m);
            dXi_H_dt[index_HI]+=Bi[m]*(Ai[m]*Xe*Np-HIA.X(index_HI)); // here Ai=Ai_rec/Bi_phot
        }
        
        return;
    }
    
    //==================================================================
    // Bi[m] is independent of Te
    //==================================================================
    void evaluate_effective_Rci_Ric_terms_for_Te_derivative(double z, double Xe, double Nc, 
                                                            vector<double> &Ai, vector<double> &Bi, 
                                                            int (*get_HI_index)(int), double *dXi_H_dt)
    {
        int index_HI;
        for(int m=0; m<(int) Ai.size(); m++)
        {
            index_HI=get_HI_index(m);
            dXi_H_dt[index_HI]+=Bi[m]*Ai[m]*Xe*Nc;
        }
        
        return;
    }

    //==================================================================
    void df_dXHI1s_evaluate_effective_Rci_Ric_terms(double z, double Ne, double Nc, 
                                                    vector<double> &Ai, vector<double> &Bi, 
                                                    int (*get_HI_index)(int), double *dXi_H_dt)
    {
        // Xe= 1-X1s + (fHe-XHeI1s) & Nc=NH*(1-X1s) --> d (Xe*Nc) / dXH1s = -(Nc+Ne)
        int index_HI;
        for(int m=0; m<(int) Ai.size(); m++)
        {
            index_HI=get_HI_index(m);
            dXi_H_dt[index_HI]+=-Bi[m]*Ai[m]*(Nc+Ne);
        }
        
        return;
    }

    //==================================================================
    void df_dXHeI1s_evaluate_effective_Rci_Ric_terms(double z, double Nc, vector<double> &Ai, vector<double> &Bi, 
                                                     int (*get_HI_index)(int), double *dXi_H_dt)
    {
        // Xe= 1-X1s + (fHe-XHeI1s) & Nc=NH*(1-X1s) --> d (Xe*Nc) / dXHeI1s = -Nc
        int index_HI;
        for(int m=0; m<(int) Ai.size(); m++)
        {
            index_HI=get_HI_index(m);
            dXi_H_dt[index_HI]+=-Bi[m]*Ai[m]*Nc;
        }
        
        return;
    }
        
    //==================================================================
    void df_dXi_evaluate_effective_Rci_Ric_terms(int i, double z, vector<double> &Ai, vector<double> &Bi, 
                                                 int (*get_HI_index)(int), double *dXi_H_dt)
    {
        // i is the index of the resolved state; get_HI_index(i) cannot be XHI1s
        int index_HI=get_HI_index(i);
        dXi_H_dt[index_HI]+=-Bi[i];
        return;
    }
    
    
    //==================================================================
    //
    // effective Rij-rates
    //
    //==================================================================
    void evaluate_effective_Rij_terms(double z, double Tg, Gas_of_Atoms &HIA, vector<vector<double> > &Rij, int (*get_HI_index)(int), double *dXi_H_dt)
    {
        int index_HI_i, index_HI_j;
        double dum, xij;
        for(int i=0; i<(int)Rij.size(); i++)
        {
            index_HI_i=get_HI_index(i);

            for(int j=i+1; j<(int)Rij[i].size(); j++)
            {
                index_HI_j=get_HI_index(j);
                // for nu21 also the sign is important...
                xij=const_h_kb*HIA.Level(index_HI_i).nu_ul( HIA.Get_n_of_Level(index_HI_i), HIA.Get_n_of_Level(index_HI_j) )/Tg;
                dum=HIA.Level(index_HI_i).Get_gw()/HIA.Level(index_HI_j).Get_gw();
                //
                dum=Rij[i][j]*( HIA.X(index_HI_i) - dum*exp(-xij)*HIA.X(index_HI_j) );
                //
                dXi_H_dt[index_HI_i]-=dum;
                dXi_H_dt[index_HI_j]+=dum;
            }
        }
        
        return;
    }
    
    //==================================================================
    void df_dXi_evaluate_effective_Rij_terms(int i, double z, double Tg, Gas_of_Atoms &HIA, vector<vector<double> > &Rij, 
                                             int (*get_HI_index)(int), double *dXi_H_dt)
    {
        int index_HI_i=get_HI_index(i), index_HI_j;

        for(int j=i+1; j<(int)Rij[i].size(); j++)
        {
            index_HI_j=get_HI_index(j);
            dXi_H_dt[index_HI_i]-=Rij[i][j];
            dXi_H_dt[index_HI_j]+=Rij[i][j];
        }

        double dum, xij;
        index_HI_j=get_HI_index(i);
        for(int m=0; m<(int)Rij.size(); m++)
        {
            index_HI_i=get_HI_index(m);
            // for nu21 also the sign is important...
            xij=const_h_kb*HIA.Level(index_HI_i).nu_ul( HIA.Get_n_of_Level(index_HI_i), HIA.Get_n_of_Level(index_HI_j) )/Tg;
            dum=HIA.Level(index_HI_i).Get_gw()/HIA.Level(index_HI_j).Get_gw();
            dum=-Rij[m][i]*dum*exp(-xij);
            //
            dXi_H_dt[index_HI_i]-=dum;
            dXi_H_dt[index_HI_j]+=dum;
        }

        return;
    }
}

//======================================================================
// End of Namespace "ODE_HI_effective"                                //
//======================================================================




//======================================================================
// beginning of Namespace "ODE_HeI_effective"                         //
//======================================================================

namespace ODE_HeI_effective 
{
    
    //==================================================================
    // effective recombination and photonionization rates
    //==================================================================
    void evaluate_effective_Rci_Ric_terms(double z, double Xe, double Nc, Gas_of_HeI_Atoms &HeIA, vector<double> &Ai, 
                                          vector<double> &Bi, int (*get_HeI_index)(int), double *dXi_He_dt)
    {
        int index_HeI;
        for(int m=0; m<(int) Ai.size(); m++)
        {
            index_HeI=get_HeI_index(m);
            dXi_He_dt[index_HeI]+=Bi[m]*(Ai[m]*Xe*Nc-HeIA.Xi(index_HeI)); // here Ai=Ai_rec/Bi_phot
        }
        
        return;
    }
    
    //==================================================================
    void df_dXHI1s_evaluate_effective_Rci_Ric_terms(double z, double Nc, vector<double> &Ai, vector<double> &Bi, 
                                                    int (*get_HeI_index)(int), double *dXi_He_dt)
    {
        // Xe= 1-X1s + (fHe-XHeI1s) & Nc=NH*(fHe-XHeI1s) --> d (Xe*Nc) / dXHI1s = -Nc
        int index_HeI;
        for(int m=0; m<(int) Ai.size(); m++)
        {
            index_HeI=get_HeI_index(m);
            dXi_He_dt[index_HeI]+=-Bi[m]*Ai[m]*Nc; // here Ai=Ai_rec/Bi_phot
        }
        
        return;
    }

    //==================================================================
    void df_dXHeI1s_evaluate_effective_Rci_Ric_terms(double z, double Ne, double Nc, vector<double> &Ai, vector<double> &Bi, 
                                                     int (*get_HeI_index)(int), double *dXi_He_dt)
    {
        // Xe= 1-X1s + (fHe-XHeI1s) & Nc=NH*(fHe-XHeI1s) --> d (Xe*Nc) / dXHeI1s = -(Nc+Ne)
        int index_HeI;
        for(int m=0; m<(int) Ai.size(); m++)
        {
            index_HeI=get_HeI_index(m);
            dXi_He_dt[index_HeI]+=-Bi[m]*Ai[m]*(Nc+Ne); // here Ai=Ai_rec/Bi_phot
        }
        
        return;
    }
    
    //==================================================================
    void df_dXi_evaluate_effective_Rci_Ric_terms(int i, double z, vector<double> &Ai, vector<double> &Bi, 
                                                 int (*get_HeI_index)(int), double *dXi_He_dt)
    {
        // i is the index of the resolved state; get_HeI_index(i) cannot be XHeI1s
        int index_HeI=get_HeI_index(i);
        dXi_He_dt[index_HeI]+=-Bi[i]; 
        
        return;
    }
    
    
    //==================================================================
    // effective Rij-rates
    //==================================================================
    void evaluate_effective_Rij_terms(double z, double Tg, Gas_of_HeI_Atoms &HeIA, vector<vector<double> > &Rij, 
                                      int (*get_HeI_index)(int), double *dXi_He_dt)
    {
        int index_HeI_i, index_HeI_j;
        double dum, xij, nuij;
        for(int i=0; i<(int)Rij.size(); i++)
        {
            index_HeI_i=get_HeI_index(i);
            
            for(int j=i+1; j<(int)Rij[i].size(); j++)
            {
                index_HeI_j=get_HeI_index(j);
                // for nu21 also the sign is important...
                nuij=HeIA.Get_nu_ion(index_HeI_j)-HeIA.Get_nu_ion(index_HeI_i);
                xij=const_h_kb*nuij/Tg;
                dum=HeIA.Get_gw(index_HeI_i)/HeIA.Get_gw(index_HeI_j);              
                //
                dum=Rij[i][j]*( HeIA.Xi(index_HeI_i) - dum*exp(-xij)*HeIA.Xi(index_HeI_j) );
                //
                dXi_He_dt[index_HeI_i]-=dum;
                dXi_He_dt[index_HeI_j]+=dum;
            }
        }
        
        return;
    }
    
    //==================================================================
    void df_dXi_evaluate_effective_Rij_terms(int i, double z, double Tg, Gas_of_HeI_Atoms &HeIA, vector<vector<double> > &Rij, 
                                             int (*get_HeI_index)(int), double *dXi_He_dt)
    {
        int index_HeI_i=get_HeI_index(i), index_HeI_j;
        for(int j=i+1; j<(int)Rij[i].size(); j++)
        {
            index_HeI_j=get_HeI_index(j);
            dXi_He_dt[index_HeI_i]-=Rij[i][j];
            dXi_He_dt[index_HeI_j]+=Rij[i][j];
        }
        
        double dum, xij, nuij;
        index_HeI_j=get_HeI_index(i);
        for(int m=0; m<(int)Rij.size(); m++)
        {
            index_HeI_i=get_HeI_index(m);
            // for nu21 also the sign is important...
            nuij=HeIA.Get_nu_ion(index_HeI_j)-HeIA.Get_nu_ion(index_HeI_i);
            xij=const_h_kb*nuij/Tg;
            dum=HeIA.Get_gw(index_HeI_i)/HeIA.Get_gw(index_HeI_j);              
            dum=-Rij[m][i]*dum*exp(-xij);
            //
            dXi_He_dt[index_HeI_i]-=dum;
            dXi_He_dt[index_HeI_j]+=dum;
        }
        
        return;
    }
}

//======================================================================
// End of Namespace "ODE_HeI_effective"                               //
//======================================================================



