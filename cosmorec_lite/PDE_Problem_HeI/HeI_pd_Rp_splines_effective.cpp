//===========================================================================================================
// Authors Jeffrey Fung & Jens Chluba Feb/March 2010
//
// comment: HeI populations have to be pre-loaded; Also the effective rates have to be initialized
//          Also note that here calc_Xi and all other splines are using the resolved level index 
//          instead of the HeI-level index!
//===========================================================================================================
// 30.07.2014: Xi_Data is passed on read only
// 22.07.2014: restart of spline setup if number of redshift points changes
// 22.07.2014: changed setup routine for splines (smaller buffer region)

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>

#include "HeI_pd_Rp_splines_effective.h"
#include "get_effective_rates.HeI.h"
#include "Load.Populations.HI.h"
#include "Load.Populations.HeI.h"
#include "HeI_Atom.h"
#include "Cosmos.h"
#include "physical_consts.h"
#include "routines.h"

using namespace std;

//===========================================================================================================
static bool HeI_pd_Rp_splines_effective_spline_memory_was_allocated=0;
int HeI_pd_Rp_splines_effective_verbose=1;

//===========================================================================================================
void Set_HeI_pd_Rp_splines_effective_verbosity(int v)
{ HeI_pd_Rp_splines_effective_verbose=v; return; }

//===========================================================================================================
//
// compute pd and Rp/Rm
//
//===========================================================================================================
void calc_HeI_Rp_Rm_pd_all(double z, Gas_of_HeI_Atoms &HeIA, double Tg, double Te, double Ne, 
                          double Xp, double NH, vector<double> &Rp_Rm, vector<double> &pd)
{
    int ntot=Rp_Rm.size();
    
    vector<double> Ai(ntot);
    vector<double> Bi(ntot), Bitot(ntot);
    vector<vector<double> > RijVec(ntot, vector<double>(ntot));

    //=======================================================================================================
    get_rates_all_HeI(Tg, Ai, Bi, Bitot, RijVec);

    //=======================================================================================================
    // Ai[nr] :: effective recombination rate to level nr
    // Bi[nr] :: effective photonionization rate from level nr
    // RijVec[i][j] :: effective transition rate to between level i and j
    // Bitot[nr] :: total rate for exciting level nr
    //=======================================================================================================
    
    double Rp, Rm;  
    for(int nr=0; nr<ntot; nr++)
    {
        int index_HeI_i=get_HeI_index(nr);
        
        //===================================================================================================
        // compute Rp & Rm
        //===================================================================================================
        Rp=Ne*Xp*Bi[nr]*Ai[nr];
        Rm=Bitot[nr];

        for(int j=0; j<nr; j++) Rp+=RijVec[j][nr]*calc_HeI_Xi(z, j);
        for(int j=nr+1; j<ntot; j++) Rp+=RijVec[j][nr]*calc_HeI_Xi(z, j);
 
        Rp_Rm[nr]=Rp*NH/Rm;
        
        if(index_HeI_i==1) pd[nr]=Rm/(const_HeI_A2s_1s+Rm);
        else pd[nr]=Rm/(HeIA.Get_A(HeIA.Get_n(index_HeI_i), HeIA.Get_l(index_HeI_i), 
                                   HeIA.Get_S(index_HeI_i), HeIA.Get_J(index_HeI_i), 1, 0, 0, 0)+Rm);
    }

    return;
}

//===========================================================================================================
//
// For precalulated death probabilities and emission rates
//
//===========================================================================================================
struct HeI_Rp_pd_eff_Data
{
    vector<double> Rp;
    int memindex;
};

vector<HeI_Rp_pd_eff_Data> HeI_pd_eff_Data;
vector<HeI_Rp_pd_eff_Data> HeI_Dnem_eff_Data;

//===========================================================================================================
void clear_HeI_pd_Rp_effective_splines()
{
    if(HeI_pd_Rp_splines_effective_spline_memory_was_allocated)
    {
        for(int k=0; k<(int)HeI_pd_eff_Data.size(); k++) 
            free_spline_JC(HeI_pd_eff_Data[k].memindex, "set_up_splines_for_HeI_pd_Rp_effective");
    }
    
    HeI_pd_Rp_splines_effective_spline_memory_was_allocated=0;

    return;
}

//===========================================================================================================
void clear_HeI_pd_Rp_effective_memory()
{
    clear_HeI_pd_Rp_effective_splines();
    
    HeI_pd_eff_Data.clear();
    HeI_Dnem_eff_Data.clear();
    
    return;
}

//===========================================================================================================
void set_up_splines_for_HeI_pd_Rp_effective(double zend, double zstart,
                                            Cosmos &cosmos, Gas_of_HeI_Atoms &HeA, 
                                            double *Rp_zarr, int nz)
{   
    //===========================================================================
    int ntot=get_number_of_resolved_levels_HeI();
        
    //===========================================================================
    // check if spline memory was already allocated for as many levels
    //===========================================================================
    if(HeI_pd_Rp_splines_effective_spline_memory_was_allocated)
    {
        if(ntot+1!=(int)HeI_pd_eff_Data.size())
        { 
            clear_HeI_pd_Rp_effective_splines();
            HeI_pd_eff_Data.resize(ntot+1);
            HeI_Dnem_eff_Data.resize(ntot+1);
        }            
    }
    else 
    {
        HeI_pd_eff_Data.resize(ntot+1);
        HeI_Dnem_eff_Data.resize(ntot+1);
    }   
    
    //===========================================================================
    // check if the number of redshift points has changed
    //===========================================================================
    if(nz!=(int)HeI_pd_eff_Data[0].Rp.size())
        for(int i=0; i<ntot+1; i++)
        { 
            HeI_pd_eff_Data[i].Rp.resize(nz); 
            HeI_Dnem_eff_Data[i].Rp.resize(nz); 
        }
    
    if(HeI_pd_Rp_splines_effective_verbose>=1) cout << endl;

    //===========================================================================
    // fill ALL y-arrays
    //===========================================================================
    double z, Tg, Te, NH;
    vector<double> Rp_Rm(ntot);
    vector<double> pd(ntot);
    vector<double> Dnem_arr(ntot);
    Transition_Data_HeI_A T;
    //
    for(int i=0; i<nz; i++)
    {
        z=Rp_zarr[i];
        Tg=cosmos.TCMB(z);
        Te=Tg*calc_HI_rho(z);
        NH=cosmos.NH(z);

        calc_HeI_Rp_Rm_pd_all(z, HeA, Tg, Te, NH*calc_HI_Xe(z), cosmos.fHe()-calc_HeI_X1s(z), NH, Rp_Rm, pd);

        //=====================================================================
        // common data
        //=====================================================================
        double x_c=const_h_kb*HeA.Sing.Level(2, 0).Get_Dnu_1s2()/Tg;
        double exp_21=exp(-x_c);
        double N1s=NH*calc_HeI_X1s(z);
        double N2s=NH*calc_HeI_Xi(z, 0);

        //=====================================================================
        // Dnem for singlet 2s-state (pd==1)
        //=====================================================================
        Dnem_arr[0]=N2s/N1s-exp_21;
            
        //=====================================================================
        // compute the nP - 1s and nD - 1s Rp_Rm using the net rate equation
        //=====================================================================
        for(int res_i=1; res_i<ntot; res_i++)
        {
            int HeI_i=get_HeI_index(res_i);
            double Ni=NH*calc_HeI_Xi(z, res_i);
            T=HeA.Get_Trans_Data(HeI_i, 1, 0, 0, 0);
            double wr=1.0*T.gwp/HeA.Get_gw(HeI_i);
            //
            x_c=const_h_kb*T.Dnu/Tg;
            exp_21=exp(-x_c);
            //
            double nL=wr*Ni/N1s;
            double tauS=T.A21*pow(T.lambda21, 3)/(8.0*PI*cosmos.H(z))*(N1s/wr - Ni);
            double PS=(1.0-exp(-tauS))/tauS;
            //
            Dnem_arr[res_i]=(nL-exp_21)*(1.0+(1.0/pd[res_i]-1.0)*PS);
        }
        
        //=====================================================================
        // copy all stuff to arrays for splines
        //=====================================================================
        for(int nr=0; nr<ntot; nr++)
        { 
            HeI_pd_eff_Data[nr].Rp[i]=log(pd[nr]); 
            HeI_Dnem_eff_Data[nr].Rp[i]=Dnem_arr[nr]; 
        }
        
        //=====================================================================
        // set the the entry with res-index=ntot to zero 
        // (needed for spline interpolation)
        //=====================================================================
        HeI_pd_eff_Data[ntot].Rp[i]=HeI_Dnem_eff_Data[ntot].Rp[i]=0.0;
        
        if(HeI_pd_Rp_splines_effective_verbose>=1) 
            cout << "\033[2K" << " z= " << z << endl << "\033[1A";
    }
    
    if(HeI_pd_Rp_splines_effective_verbose>=1) 
        cout << " set_up_splines_for_HeI_pd_Rp:: rates tabulated " << endl;
    
    //========================================================================
    // calc splines + allocalate memory (if necessary)
    //========================================================================
    if(HeI_pd_Rp_splines_effective_spline_memory_was_allocated)
    {
        for(int n=0; n<ntot+1; n++)
        {
            update_spline_coeffies_JC(HeI_pd_eff_Data[n].memindex, nz, Rp_zarr, &HeI_pd_eff_Data[n].Rp[0]);
            update_spline_coeffies_JC(HeI_Dnem_eff_Data[n].memindex, nz, Rp_zarr, &HeI_Dnem_eff_Data[n].Rp[0]);
        }
    }
    else 
    {
        for(int n=0; n<ntot+1; n++)
        {
            HeI_pd_eff_Data[n].memindex=calc_spline_coeffies_JC(nz, Rp_zarr, &HeI_pd_eff_Data[n].Rp[0], 
                                                               "set_up_splines_for_HeI_pd_Rp :: pd "
                                                               + int_to_string(n));

            HeI_Dnem_eff_Data[n].memindex=calc_spline_coeffies_JC(nz, Rp_zarr, &HeI_Dnem_eff_Data[n].Rp[0], 
                                                                 "set_up_splines_for_HeI_pd_Rp :: Dnem "
                                                                 + int_to_string(n));
        }
        
        HeI_pd_Rp_splines_effective_spline_memory_was_allocated=1;
    }
    
    if(HeI_pd_Rp_splines_effective_verbose>=1) 
        cout << " set_up_splines_for_HeI_pd_Rp:: splines done." << endl; 

    return;
}

//===========================================================================================================
void set_up_splines_for_HeI_pd_Rp_effective(double zend, double zstart, Cosmos &cosmos, Gas_of_HeI_Atoms &HeA)
{   
    int nz=2000;
    double *Rp_zarr=new double[nz];
    
    //===========================================================================
    // fill z-array (z0 < z1 < ... < zn)
    //===========================================================================
    init_xarr(zend/1.0001, zstart*1.0001, Rp_zarr, nz, 1, 0);
    
    //===========================================================================
    // compute splines
    //===========================================================================
    set_up_splines_for_HeI_pd_Rp_effective(zend, zstart, cosmos, HeA, Rp_zarr, nz);
    
    //========================================================================
    // clean up
    //========================================================================
    delete[] Rp_zarr;
    
    return;
}

//===========================================================================================================
void set_up_splines_for_HeI_pd_Rp_effective(double zend, double zstart, 
                                           Cosmos &cosmos, Gas_of_HeI_Atoms &HeA,
                                           const vector<vector<double> > &X_Data)
{   
    int nmax=X_Data.size();
    vector<double> Rp_zarr(nmax, 0.0);
    
    for(int k=0; k<nmax; k++)
    {
        double z=X_Data[k][0];
        if(z>=zend/1.0001 && z<=zstart*1.0001) Rp_zarr[nmax-1-k]=z;
    }
    
    for(int k=0; k<nmax; k++)
        if(Rp_zarr[0]==0.0) Rp_zarr.erase(Rp_zarr.begin());
        else break;
    
    int nz=Rp_zarr.size();
    
    //===========================================================================
    // fill z-array (z0 < z1 < ... < zn)
    //===========================================================================
    init_xarr(zend/1.0001, zstart*1.0001, &Rp_zarr[0], nz, 0, 0);

    //===========================================================================
    // compute splines
    //===========================================================================
    set_up_splines_for_HeI_pd_Rp_effective(zend, zstart, cosmos, HeA, &Rp_zarr[0], nz);
    
    //========================================================================
    // clean up
    //========================================================================
    Rp_zarr.clear();
    
    return;
}

//===========================================================================================================
//
// pd-functions
//
//===========================================================================================================
double calc_HeI_pd_i_splines_effective(double z, int ires)
{ return exp(calc_spline_JC(z, HeI_pd_eff_Data[ires].memindex)); }

//===========================================================================================================
//
// Dnem-functions
//
//===========================================================================================================
double calc_HeI_Dnem_i_splines_effective(double z, int ires)
{ return calc_spline_JC(z, HeI_Dnem_eff_Data[ires].memindex); }

//===========================================================================================================
//===========================================================================================================
