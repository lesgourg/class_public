//==========================================================================================
// Author Jens Chluba Aug/Sept 2010
// purpose: compute transitions rates for np, ns & nd-states of hydrogen
// last modification: Aug 2012
//==========================================================================================

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <vector>

#include "HI_Transition_Data.h"
#include "Oscillator_strength.h"
#include "physical_consts.h"

using namespace std;

//==========================================================================================
//==========================================================================================
const int nmax_setup=10;

namespace HI_Transition_Data 
{
    //======================================================================================
    // Gamma_n=A_tot/(4 pi) for the first nmax p-states (vacuum)
    //======================================================================================
    vector<double> Gamma_np_vals;
    
    //======================================================================================
    // get np vacuum life-time; if needed set up
    //======================================================================================
    double Get_Gamma_np(int n)
    {
        if(Gamma_np_vals.size()==0)
        {
            Gamma_np_vals.resize(nmax_setup+1);
            
            Gamma_np_vals[0]=0.0; // n=0 (n/a)
            Gamma_np_vals[1]=0.0; // n=1
            Gamma_np_vals[2]=A_SH(1, 1.0, 2, 1, 1, 0)/FOURPI; // n=2
            
            for(int d=3; d<=nmax_setup; d++)
            {
                double Atot=A_SH(1, 1.0, d, 1, 1, 0);
                Atot+=A_SH(1, 1.0, d, 1, 2, 0);
                
                for(int k=3; k<d; k++) 
                    Atot+=A_SH(1, 1.0, d, 1, k, 0)+A_SH(1, 1.0, d, 1, k, 2);
                
                Gamma_np_vals[d]=Atot/FOURPI;
            }
        }
        
        return Gamma_np_vals[n];
    }
    
    //======================================================================================
    // A_npks transition rates for the first nmax p-states (vacuum)
    // k==1: just Ly-n Anp->1s rates
    // k>1 : rates into (!) the np-state
    //       n<k --> Aks-->np
    //       n>k --> 3*Anp->ks
    //======================================================================================
    vector<vector<double> > A_npks_vals;
    
    double Get_A_npks(int k, int n)
    {
        if(A_npks_vals.size()==0)
        {
            vector<double> dum(nmax_setup+1, 0);
            A_npks_vals.push_back(dum); // n=0 (n/a)
            
            // Ly-n rates
            for(int nv=2; nv<=nmax_setup; nv++)
                dum[nv]=A_SH(1, 1.0, nv, 1, 1, 0);
            
            A_npks_vals.push_back(dum);

            for(int kv=2; kv<=nmax_setup; kv++)
            {
                dum[kv]=0.0; // no transition np->ns
                
                // transitions from higher levels
                for(int nv=1; nv<kv; nv++)
                    dum[nv]=A_SH(1, 1.0, kv, 0, nv, 1);
                
                // transitions from lower levels
                for(int nv=kv+1; nv<=nmax_setup; nv++)
                    dum[nv]=3.0*A_SH(1, 1.0, nv, 1, kv, 0);

                A_npks_vals.push_back(dum);
            }
        }
        
        return A_npks_vals[k][n];
    }

    //======================================================================================
    // access to Ly-n rates (vacuum)
    //======================================================================================
    double Get_A_np1s(int n){ return Get_A_npks(1, n); }

    //======================================================================================
    // A_npkd transition rates for the first nmax p-states (vacuum)
    // k>1 : rates into (!) the np-state
    //       n<k --> Akd-->np 
    //       n>k --> 3/5*Anp->kd
    //======================================================================================
    vector<vector<double> > A_npkd_vals;
    
    double Get_A_npkd(int k, int n)
    {
        if(A_npkd_vals.size()==0)
        {
            vector<double> dum(nmax_setup+1, 0);
            A_npkd_vals.push_back(dum); // n=0 (n/a)
            A_npkd_vals.push_back(dum); // n=1 (n/a)
            A_npkd_vals.push_back(dum); // n=2 (n/a)
            
            for(int kv=3; kv<=nmax_setup; kv++)
            {
                dum[kv]=0.0; // no transition np->nd
                
                // transitions from higher levels
                for(int nv=1; nv<kv; nv++)
                    dum[nv]=A_SH(1, 1.0, kv, 2, nv, 1);
                
                // transitions from lower levels
                for(int nv=kv+1; nv<=nmax_setup; nv++)
                    dum[nv]=3.0/5.0*A_SH(1, 1.0, nv, 1, kv, 2);
                
                A_npkd_vals.push_back(dum);
            }
        }

        return A_npkd_vals[k][n];
    } 
}

//==========================================================================================
//==========================================================================================
