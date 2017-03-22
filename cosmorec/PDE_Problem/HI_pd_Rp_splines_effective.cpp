//===========================================================================================================
// Author Jens Chluba May 2009
// comment: HI populations have to be pre-loaded; Also the effective rates have to be initialized
//===========================================================================================================
// 30.07.2014: Xi_Data is passed on read only
// 22.07.2014: restart of spline setup if number of redshift points changes
// 22.07.2014: changed setup routine for splines (smaller buffer region)

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>

#include "HI_pd_Rp_splines_effective.h"
#include "get_effective_rates.HI.h"
#include "Load.Populations.HI.h"
#include "Atom.h"
#include "Cosmos.h"
#include "physical_consts.h"
#include "routines.h"

using namespace std;

//===========================================================================================================
static bool HI_pd_Rp_splines_effective_spline_memory_was_allocated=0;
int HI_pd_Rp_splines_effective_verbose=1;

//===========================================================================================================
void Set_HI_pd_Rp_splines_effective_verbosity(int v)
{ HI_pd_Rp_splines_effective_verbose=v; return; }

//===========================================================================================================
//
// compute pd and Rp/Rm
//
//===========================================================================================================
void calc_HI_Rp_Rm_pd_all(double z, Gas_of_Atoms &HIA, double Tg, double Te, double Ne, double Xp, double NH, 
                          vector<double> &Rp_Rm, vector<double> &pd)
{
    int ntot=Rp_Rm.size();
    
    vector<double> Ai(ntot);
    vector<double> Bi(ntot), Bitot(ntot);
    vector<vector<double> > RijVec(ntot, vector<double>(ntot));

    //=======================================================================================================
    get_rates_all(Tg, Te, Ai, Bi, Bitot, RijVec);
    
    //=======================================================================================================
    // Ai[nr] :: effective recombination rate to level nr
    // Bi[nr] :: effective photonionization rate from level nr
    // RijVec[i][j] :: effective transition rate to between level i and j
    // Bitot[nr] :: total rate for exciting level nr
    //=======================================================================================================
    
    double Rp, Rm;  
    for(int nr=0; nr<ntot; nr++)
    {
        int index_HI_i=get_HI_index(nr);
        
        //===================================================================================================
        // compute Rp & Rm
        //===================================================================================================
        Rp=Ne*Xp*Bi[nr]*Ai[nr];
        Rm=Bitot[nr];

        for(int j=0; j<nr; j++)
        {
            int index_HI_j=get_HI_index(j);
            Rp+=RijVec[j][nr]*calc_HI_Xi(z, index_HI_j);
        }
        
        for(int j=nr+1; j<ntot; j++)
        {
            int index_HI_j=get_HI_index(j);
            Rp+=RijVec[j][nr]*calc_HI_Xi(z, index_HI_j);
        }

        Rp_Rm[nr]=Rp*NH/Rm;
        
        if(index_HI_i==1) pd[nr]=Rm/(const_HI_A2s_1s+Rm);
        else pd[nr]=Rm/(HIA.Level(index_HI_i).Get_A21(1, 0)+Rm);
    }

    return;
}

//===========================================================================================================
//
// For precalulated death probabilities and emission rates
//
//===========================================================================================================
struct HI_Rp_pd_eff_Data
{
    vector<double> Rp;
    int memindex;
};

vector<HI_Rp_pd_eff_Data> HI_pd_eff_Data;
vector<HI_Rp_pd_eff_Data> HI_Dnem_eff_Data;

//===========================================================================================================
void clear_HI_pd_Rp_effective_splines()
{
    if(HI_pd_Rp_splines_effective_spline_memory_was_allocated)
    {
        for(int k=0; k<(int)HI_pd_eff_Data.size(); k++) 
            free_spline_JC(HI_pd_eff_Data[k].memindex, "set_up_splines_for_HI_pd_Rp_effective");
    }
    
    HI_pd_Rp_splines_effective_spline_memory_was_allocated=0;

    return;
}

//===========================================================================================================
void clear_HI_pd_Rp_effective_memory()
{
    clear_HI_pd_Rp_effective_splines();
    
    HI_pd_eff_Data.clear();
    HI_Dnem_eff_Data.clear();
    
    return;
}

//===========================================================================================================
void set_up_splines_for_HI_pd_Rp_effective(double zend, double zstart, 
                                           Cosmos &cosmos, Gas_of_Atoms &HA, 
                                           double *Rp_zarr, int nz)
{   
    //===========================================================================
    int ntot=get_number_of_resolved_levels();
        
    //===========================================================================
    // check if spline memory was already allocated for as many levels
    //===========================================================================
    if(HI_pd_Rp_splines_effective_spline_memory_was_allocated)
    {
        if(ntot+1!=(int)HI_pd_eff_Data.size())
        {
            clear_HI_pd_Rp_effective_splines();
            HI_pd_eff_Data.resize(ntot+1);
            HI_Dnem_eff_Data.resize(ntot+1);
        }
    }
    else 
    {
        HI_pd_eff_Data.resize(ntot+1);
        HI_Dnem_eff_Data.resize(ntot+1);
    }
    
    //===========================================================================
    // check if the number of redshift points has changed
    //===========================================================================
    if(nz!=(int)HI_pd_eff_Data[0].Rp.size())
        for(int i=0; i<ntot+1; i++)
        { 
            HI_pd_eff_Data[i].Rp.resize(nz); 
            HI_Dnem_eff_Data[i].Rp.resize(nz); 
        }
        
    if(HI_pd_Rp_splines_effective_verbose>=1) cout << endl;

    //===========================================================================
    // fill ALL y-arrays
    //===========================================================================
    double z, Tg, Te, NH;
    vector<double> Rp_Rm(ntot);
    vector<double> pd(ntot);
    vector<double> Dnem_arr(ntot);
    //
    for(int i=0; i<nz; i++)
    {
        z=Rp_zarr[i];
        Tg=cosmos.TCMB(z);
        Te=Tg*calc_HI_rho(z);
        NH=cosmos.NH(z);
        
        calc_HI_Rp_Rm_pd_all(z, HA, Tg, Te, NH*calc_HI_Xe(z), 1.0-calc_HI_X1s(z), NH, Rp_Rm, pd);

        //=====================================================================
        // common data
        //=====================================================================
        double nu21=HA.Level(2, 1).Get_Dnu_1s();
        double x_c=const_h_kb*nu21/Tg;
        double exp_21=exp(-x_c);
        double N1s=NH*calc_HI_Xi(z, 0);
        double N2s=NH*calc_HI_Xi(z, 1);
        double N2p=NH*calc_HI_Xi(z, 2);

        //=====================================================================
        // compute Dnem for 2s-state (pd==1)
        //=====================================================================
        int res_i=get_res_index(1);
        Dnem_arr[res_i]=N2s/N1s-exp_21;
        
        //=====================================================================
        // compute Dnem for 2p-state
        //=====================================================================
        res_i=get_res_index(2);
        double nL=N2p/3.0/N1s;      
        double tauS=HA.Level(2, 1).Get_A21(1, 0)*pow(HA.Level(2, 1).Get_lambda21(1, 0), 3)
                    /(8.0*PI*cosmos.H(z))*(N1s*3.0 - N2p);      
        double PS=(1.0-exp(-tauS))/tauS;

        Dnem_arr[res_i]=(nL-exp_21)*(1.0+(1.0/pd[res_i]-1.0)*PS);
        
        //=====================================================================
        // compute the nP - 1s Rp_Rm using the net rate equation
        //=====================================================================
        for(int nr=0; nr<ntot; nr++)
        {
            int HI_i=get_HI_index(nr);
            
            if(HA.Level(HI_i).Get_l()==1 && HA.Level(HI_i).Get_n()>2)
            {
                res_i=nr;
                double Ni=NH*calc_HI_Xi(z, HI_i);
                //
                nu21=HA.Level(HI_i).Get_Dnu_1s();
                x_c=const_h_kb*nu21/Tg;
                exp_21=exp(-x_c);
                //
                nL=Ni/3.0/N1s;
                tauS=HA.Level(HI_i).Get_A21(1, 0)*pow(HA.Level(HI_i).Get_lambda21(1, 0), 3)
                     /(8.0*PI*cosmos.H(z))*(N1s*3.0 - Ni);
                PS=(1.0-exp(-tauS))/tauS;
                //
                Dnem_arr[res_i]=(nL-exp_21)*(1.0+(1.0/pd[res_i]-1.0)*PS);
            }
            else if( (HA.Level(HI_i).Get_l()==0 || HA.Level(HI_i).Get_l()==2) && HA.Level(HI_i).Get_n()>2)
            {
                res_i=nr;
                double Ni=NH*calc_HI_Xi(z, HI_i);
                //
                nu21=HA.Level(HI_i).Get_Dnu_1s();
                x_c=const_h_kb*nu21/Tg;
                exp_21=exp(-x_c);
                //
                Dnem_arr[res_i]=Ni/N1s/(2.0*HA.Level(HI_i).Get_l()+1.0)-exp_21;
            }
        }
        
        //=====================================================================
        // copy all stuff to arrays for splines
        //=====================================================================
        for(int nr=0; nr<ntot; nr++)
        { 
            HI_pd_eff_Data[nr].Rp[i]=log(pd[nr]); 
            HI_Dnem_eff_Data[nr].Rp[i]=Dnem_arr[nr]; 
        }
        
        //=====================================================================
        // set the the entry with res-index=ntot to zero 
        // (needed for spline interpolation)
        //=====================================================================
        HI_pd_eff_Data[ntot].Rp[i]=HI_Dnem_eff_Data[ntot].Rp[i]=0.0;
        
        if(HI_pd_Rp_splines_effective_verbose>=1) 
            cout << "\033[2K" << " z= " << z << endl << "\033[1A";
    }
    
    if(HI_pd_Rp_splines_effective_verbose>=1) 
        cout << " set_up_splines_for_HI_pd_Rp:: rates tabulated " << endl;
    
    //========================================================================
    // calc splines + allocalate memory (if necessary)
    //========================================================================
    if(HI_pd_Rp_splines_effective_spline_memory_was_allocated)
    {
        for(int n=0; n<ntot+1; n++)
        {
            update_spline_coeffies_JC(HI_pd_eff_Data[n].memindex, nz, Rp_zarr, &HI_pd_eff_Data[n].Rp[0]);
            update_spline_coeffies_JC(HI_Dnem_eff_Data[n].memindex, nz, Rp_zarr, &HI_Dnem_eff_Data[n].Rp[0]);
        }
    }
    else 
    {
        for(int n=0; n<ntot+1; n++)
        {
            HI_pd_eff_Data[n].memindex=calc_spline_coeffies_JC(nz, Rp_zarr, &HI_pd_eff_Data[n].Rp[0], 
                                                               "set_up_splines_for_HI_pd_Rp :: pd "
                                                               + int_to_string(n));

            HI_Dnem_eff_Data[n].memindex=calc_spline_coeffies_JC(nz, Rp_zarr, &HI_Dnem_eff_Data[n].Rp[0], 
                                                                 "set_up_splines_for_HI_pd_Rp :: Dnem "
                                                                 + int_to_string(n));
        }
        
        HI_pd_Rp_splines_effective_spline_memory_was_allocated=1;
    }
    
    if(HI_pd_Rp_splines_effective_verbose>=1) 
        cout << " set_up_splines_for_HI_pd_Rp:: splines done." << endl; 
    
    return;
}

//===========================================================================================================
void set_up_splines_for_HI_pd_Rp_effective(double zend, double zstart, Cosmos &cosmos, Gas_of_Atoms &HA)
{   
    int nz=2000;
    double *Rp_zarr=new double[nz];

    //===========================================================================
    // fill z-array (z0 < z1 < ... < zn)
    //===========================================================================
    init_xarr(zend/1.0001, zstart*1.0001, Rp_zarr, nz, 0, 0);
    
    //===========================================================================
    // compute splines
    //===========================================================================
    set_up_splines_for_HI_pd_Rp_effective(zend, zstart, cosmos, HA, Rp_zarr, nz);
    
    //========================================================================
    // clean up
    //========================================================================
    delete[] Rp_zarr;
        
    return;
}

//===========================================================================================================
void set_up_splines_for_HI_pd_Rp_effective(double zend, double zstart, 
                                           Cosmos &cosmos, Gas_of_Atoms &HA,
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
    set_up_splines_for_HI_pd_Rp_effective(zend, zstart, cosmos, HA, &Rp_zarr[0], nz);
    
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
double calc_HI_pd_nl_splines_effective(double z, int n, int l)
{ return exp(calc_spline_JC(z, HI_pd_eff_Data[get_res_index(n*(n-1)/2+l)].memindex)); }

double calc_HI_pd_i_splines_effective(double z, int i)
{ return exp(calc_spline_JC(z, HI_pd_eff_Data[get_res_index(i)].memindex)); }

//===========================================================================================================
//
// Dnem-functions
//
//===========================================================================================================
double calc_HI_Dnem_nl_splines_effective(double z, int n, int l)
{ return calc_spline_JC(z, HI_Dnem_eff_Data[get_res_index(n*(n-1)/2+l)].memindex); }

double calc_HI_Dnem_i_splines_effective(double z, int i)
{ return calc_spline_JC(z, HI_Dnem_eff_Data[get_res_index(i)].memindex); }

