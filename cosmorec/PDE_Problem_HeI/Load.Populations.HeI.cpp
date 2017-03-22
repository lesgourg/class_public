//====================================================================================================================
// Authors Jeffrey Fung & Jens Chluba Feb/March 2010
//====================================================================================================================
// 30.07.2014: Xi_Data is passed on read only

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>

#include "Load.Populations.HeI.h"
#include "HeI_Atom.h"
#include "routines.h"

using namespace std;

static bool Load_Populations_HeI_spline_memory_was_allocated=0;
int Load_Populations_HeI_verbose=1;

//==================================================================================
//
// global variable that will contain the data for spline routine
//
//==================================================================================
struct spline_data_HeI
{
    int memindex;
};

vector<spline_data_HeI > HeI_Xi_spline;

//==================================================================================
void Set_Load_Populations_HeI_verbosity(int v)
{ Load_Populations_HeI_verbose=v; return; }

//==================================================================================
//
// access data (local)
//
//==================================================================================
double HeI_Get_z(int iz, const vector<vector<double> > &Xi_Data)
{ return Xi_Data[iz][0]; }

double HeI_Get_Xi(int iz, const vector<vector<double> > &Xi_Data, int i)
{ return Xi_Data[iz][i]; }

//==================================================================================
//
// functions
//
//==================================================================================
double calc_HeI_X1s(double z)
{ return exp(calc_spline_JC(z, HeI_Xi_spline[0].memindex)); }

double calc_HeI_Xi(double z, int ires)
{ return exp(calc_spline_JC(z, HeI_Xi_spline[ires+1].memindex)); }

//==================================================================================
//
// Load the data for Xi
//
//==================================================================================
void compute_Xi_HeI_splines(double zs, double ze, 
                            Gas_of_HeI_Atoms &HeI_Atoms,
                            const vector<vector<double> > &X_Data)
{
    //==============================================================================
    // allocate memory for splines
    //==============================================================================
    int nz=X_Data.size()-1;
    int nres=X_Data[0].size()-2; // number of resolved levels (1st is z, 2nd X1s)
    
    //==============================================================================
    // check if the number of levels or redshift-points has changes
    //==============================================================================
    if(Load_Populations_HeI_spline_memory_was_allocated)
        if((int)HeI_Xi_spline.size()!=nres) clear_Xi_HeI();    
    
    //==============================================================================
    // check if redshift ranges is withing bounds
    //==============================================================================
    if(zs>HeI_Get_z(0, X_Data) || ze<HeI_Get_z(nz-1, X_Data))
    { 
        cerr << " compute_Xi_HeI_splines: check redshift range " << zs << " " 
             << ze << " " << HeI_Get_z(0, X_Data) << " " 
             << HeI_Get_z(nz-1, X_Data) << endl; 
        exit(0); 
    }  
    
    double *za=new double[nz];
    double *ya=new double[nz];
    
    //==============================================================================
    // now compute all the spline arrays
    //==============================================================================
    for(int i=0; i<nz; i++) za[nz-1-i]=HeI_Get_z(i, X_Data);
    
    if(Load_Populations_HeI_spline_memory_was_allocated)
    {
        // HeI-X1s
        for(int k=0; k<nz; k++) ya[nz-1-k]=log(HeI_Get_Xi(k, X_Data, 1));
        update_spline_coeffies_JC(HeI_Xi_spline[0].memindex, nz, za, ya);

        // other resolved states
        for(int i=0; i<nres; i++)
        { 
            for(int k=0; k<nz; k++) ya[nz-1-k]=log(HeI_Get_Xi(k, X_Data, i+2));
            update_spline_coeffies_JC(HeI_Xi_spline[i+1].memindex, nz, za, ya);
        }
    }
    else
    {
        HeI_Xi_spline.clear();
        spline_data_HeI dummy; 
        
        // HeI-X1s
        for(int k=0; k<nz; k++) ya[nz-1-k]=log(HeI_Get_Xi(k, X_Data, 1));
        dummy.memindex=calc_spline_coeffies_JC(nz, za, ya, "compute_Xi_HeI_splines:: HeI-X1s");
        HeI_Xi_spline.push_back(dummy);

        for(int i=0; i<nres; i++)
        { 
            for(int k=0; k<nz; k++) ya[nz-1-k]=log(HeI_Get_Xi(k, X_Data, i+2));
            dummy.memindex=calc_spline_coeffies_JC(nz, za, ya, "compute_Xi_HeI_splines::HeI-res: "+int_to_string(i));
            HeI_Xi_spline.push_back(dummy);
        }
        
        Load_Populations_HeI_spline_memory_was_allocated=1;
    }
    
    //==============================================================================
    // clean up
    //==============================================================================
    delete [] za;
    delete [] ya;
    
    return;      
}

//==================================================================================
//
// Clearing the data for Xi
//
//==================================================================================
void clear_Xi_HeI()
{
    if(Load_Populations_HeI_spline_memory_was_allocated)
        for(int i=0; i<(int)HeI_Xi_spline.size(); i++)
            free_spline_JC(HeI_Xi_spline[i].memindex);
    
    HeI_Xi_spline.clear();
    
    Load_Populations_HeI_spline_memory_was_allocated=0;
    
    return;
}

//====================================================================================================================
//====================================================================================================================
