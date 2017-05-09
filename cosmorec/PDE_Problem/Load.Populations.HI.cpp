//==================================================================================
// Author Jens Chluba May 2009
//==================================================================================
// 30.07.2014: Xi_Data is passed on read only

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>

#include "Load.Populations.HI.h"
#include "Atom.h"
#include "routines.h"

using namespace std;

static bool Load_Populations_HI_spline_memory_was_allocated=0;
int Load_Populations_HI_verbose=1;

//==================================================================================
//
// global variable that will contain the data for spline routine
//
//==================================================================================
struct spline_data_HI
{
    int memindex;
};

vector<spline_data_HI > HI_Xi_spline;

//==================================================================================
void Set_Load_Populations_HI_verbosity(int v)
{ Load_Populations_HI_verbose=v; return; }

//==================================================================================
//
// to read the results of a previous calculation from a file
//
//==================================================================================
void read_HI_Population_from_file(ifstream &ifile, int nS,
                                  vector<vector<double> > &Xi_Data)
{
    //==============================================================================
    // vector<vector<double> > &Xi_Data will contain
    // z, Xe, X1s, X2s, X2p, ...., rho, Xp, XHeII
    //==============================================================================
    if(Load_Populations_HI_verbose>=1)
        cout << "\n read_HI_Population_from_file:"
             << " Reading HI data from previous calculation..." << endl;
    
    int neq=1+nS*(nS+1)/2+2+2;
    vector<double> vdum(neq);

    while(ifile)
    {
        for(int i=0; i<neq; i++) ifile >> vdum[i]; 
        Xi_Data.push_back(vdum);
    }
    
    if(Load_Populations_HI_verbose>=1)
        cout << " read_HI_Population_from_file: Number of redshift points: " 
             << Xi_Data.size()-1 << endl;
    
    return;
}

//==================================================================================
//
// access data (local)
//
//==================================================================================
double HI_Get_z(int iz, const vector<vector<double> > &Xi_Data)
{ return Xi_Data[iz][0]; }

double HI_Get_Xi(int iz, const vector<vector<double> > &Xi_Data, int i)
{ return Xi_Data[iz][i]; }

//==================================================================================
//
// functions
//
//==================================================================================
double calc_HI_Xe(double z)
{ return exp(calc_spline_JC(z, HI_Xi_spline[0].memindex)); }

double calc_HI_rho(double z)
{ return calc_spline_JC(z, HI_Xi_spline[HI_Xi_spline.size()-1].memindex)+1.0; }

double calc_HI_X1s(double z)
{ return exp(calc_spline_JC(z, HI_Xi_spline[1].memindex)); }

double calc_HI_Xnl(double z, int n, int l)
{ return exp(calc_spline_JC(z, HI_Xi_spline[n*(n-1)/2+l+1].memindex)); }

double calc_HI_Xi(double z, int i)
{ return exp(calc_spline_JC(z, HI_Xi_spline[i+1].memindex)); }

//==================================================================================
//
// Load the data for Xi
//
//==================================================================================
void compute_Xi_HI_splines(double zs, double ze, Gas_of_Atoms &H_Atoms, 
                           const vector<vector<double> > &X_Data)
{
    //==============================================================================
    // allocate memory for splines
    //==============================================================================
    int nz=X_Data.size()-1;
    
    //==============================================================================
    // check if the number of levels has changes
    //==============================================================================
    if(Load_Populations_HI_spline_memory_was_allocated)
        if(HI_Xi_spline.size()!=X_Data[0].size()-3) clear_Xi_HI();

    //==============================================================================
    // check if redshift ranges is withing bounds
    //==============================================================================
    if(zs>HI_Get_z(0, X_Data) || ze<HI_Get_z(nz-1, X_Data))
    { 
        cerr << " compute_Xi_HI_splines: check redshift range " << zs << " "
             << ze << " " << HI_Get_z(0, X_Data) << " " 
             << HI_Get_z(nz-1, X_Data) << endl; 
        exit(0); 
    }  
    
    double *za=new double[nz];
    double *ya=new double[nz];
    
    //==============================================================================
    // now compute all the spline arrays
    //==============================================================================
    for(int i=0; i<nz; i++) za[nz-1-i]=HI_Get_z(i, X_Data);
    
    if(Load_Populations_HI_spline_memory_was_allocated)
    {
        for(int i=1; i<(int)X_Data[0].size()-3; i++)
        { 
            for(int k=0; k<nz; k++) ya[nz-1-k]=log(HI_Get_Xi(k, X_Data, i));
            update_spline_coeffies_JC(HI_Xi_spline[i-1].memindex, nz, za, ya);
        }
        
        for(int k=0; k<nz; k++) ya[nz-1-k]=HI_Get_Xi(k, X_Data, X_Data[0].size()-3)-1.0;
        update_spline_coeffies_JC(HI_Xi_spline[HI_Xi_spline.size()-1].memindex, nz, za, ya);
    }
    else
    {
        HI_Xi_spline.clear();
        spline_data_HI dummy; 
        
        for(int i=1; i<(int)X_Data[0].size()-3; i++)
        { 
            for(int k=0; k<nz; k++) ya[nz-1-k]=log(HI_Get_Xi(k, X_Data, i));
            dummy.memindex=calc_spline_coeffies_JC(nz, za, ya, "compute_Xi_HI_splines::"+int_to_string(i));
            HI_Xi_spline.push_back(dummy);
        }
        
        for(int k=0; k<nz; k++) ya[nz-1-k]=HI_Get_Xi(k, X_Data, X_Data[0].size()-3)-1.0;
        dummy.memindex=calc_spline_coeffies_JC(nz, za, ya, 
                                               "compute_Xi_HI_splines::"+int_to_string(X_Data[0].size()-3));
        HI_Xi_spline.push_back(dummy);
        
        Load_Populations_HI_spline_memory_was_allocated=1;
    }
    
    //==============================================================================
    // clean up
    //==============================================================================
    delete [] za;
    delete [] ya;
    
    return;      
}

//==================================================================================
void Load_Xi_HI(int nS, double zs, double ze, Gas_of_Atoms &H_Atoms, string fname)
{
    if(Load_Populations_HI_verbose>=1)
        cout << "\n Loading solution from : " << fname << " || nShells=" << nS << endl;
    
    ifstream ifile(fname.c_str());
    if(ifile.is_open()!=1)
    { 
        cerr << " Error opening file." << fname << " Exiting." << endl; 
        exit(0); 
    }
    
    //==============================================================================
    // read populations
    //==============================================================================
    vector<vector<double> > X_Data;
    read_HI_Population_from_file(ifile, nS, X_Data);
    ifile.close();

    compute_Xi_HI_splines(zs, ze, H_Atoms, X_Data);
    
    X_Data.clear();
    
    return;      
}

//==================================================================================
//
// Clearing the data for Xi
//
//==================================================================================
void clear_Xi_HI()
{
    if(Load_Populations_HI_spline_memory_was_allocated)
        for(int i=0; i<(int)HI_Xi_spline.size(); i++) 
            free_spline_JC(HI_Xi_spline[i].memindex);
    
    HI_Xi_spline.clear();
    
    Load_Populations_HI_spline_memory_was_allocated=0;
    
    return;
}

//==================================================================================
//==================================================================================
