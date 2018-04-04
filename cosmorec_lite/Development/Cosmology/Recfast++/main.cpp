//========================================================================================
// Author: Jens Chluba
// Last modification: Oct 2010
// CITA, University of Toronto
// All rights reserved.
//========================================================================================

#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "recombination.Recfast.h"
#include "cosmology.Recfast.h"
#include "constants.Recfast.h"

using namespace std;

bool verbose=1;

//========================================================================================
// define default HI A2s1s two-photon decay rate. Added to allow changing it globally.
// (added 23.07.2014 JC)
//========================================================================================
double RECFAST_atomic_constants::RF_Lam2s1sH = 8.22458;  // Recfast++ default value

//========================================================================================
void show_cosmological_parameters(vector<double> &param)
{
    cout << "\n Recfast++ :: show_cosmological_parameters: " << endl;
    cout << " running Recfast++ history for the following parameters: " << endl;   
    
    cout << "\n zs: " << param[1] << "\t ze: " << param[2] << "\t nzpts: " << (int)param[0] << endl;
    cout << " Y_p: " << param[3] << "\t TCMB: " << param[4] << endl;
    cout << " OmegaM: " << param[5] << "\t OmegaB: " << param[6] 
         << " OmegaL: " << param[7] << "\t OmegaK: " << param[8] << endl;
    cout << " Omega_rel (photons & neutrinos): " << calc_Orel(param[4], param[10], param[9]) 
         << " Nnu= " << param[10] << endl;
    //
    cout << " Hubble constant in units of 100 km s-1 Mpc-1: " << param[9] << endl;
    cout << " Fudge factor for H recombination: " << param[11] << endl;

    if(param[12]!=0.0) 
        cout << " DM-annihilations with efficiency " << param[12] << " will be included " << endl;
    if(param[13]==1) 
        cout << " Will include recombination corrections according to Chluba & Thomas 2010 " << endl;
        
    cout << endl;
    
    return;
}

//========================================================================================
void read_parameter_file(string fname, vector<double> &params)
{
    if(verbose) 
        cout << " Recfast++ :: read_parameter_file : reading parameters from: " 
             << fname << endl;
    
    ifstream ifile(fname.c_str());
    if(!ifile)
    { 
        cerr << " Recfast++ :: read_parameter_file :" 
             << " error reading parameter-file. Exiting. " << endl; 
        exit(0); 
    }
    
    //====================================================================================
    ifile >> params[3];  // Yp
    ifile >> params[4];  // Temperature of CMB at z=0
    ifile >> params[5];  // Omega matter 
    ifile >> params[6];  // Omega Baryons 
    ifile >> params[7];  // Omega Lambda
    ifile >> params[8];  // Omega Curvature 
    ifile >> params[9];  // h100
    ifile >> params[10]; // effective number of neutrinos 
    ifile >> params[11]; // fudge-factor; normally F=1.14
    
    ifile >> params[12]; // fDM [eV/s] which gives annihilation efficiency; 
                         // typical value fDM=2.0e-24 eV/s (see Chluba 2010 for definitions)
    ifile >> params[13]; // switch on/off recombination corrections (Chluba & Thomas 2010)


    //====================================================================================
    // compute Omega_L according to Omega_K
    //====================================================================================
    if(params[7]<=0.0) 
        params[7]=1.0-params[5]-params[8]-calc_Orel(params[4], params[10], params[9]);
    
    ifile.close();
    
    if(verbose) show_cosmological_parameters(params);
    
    return;
}

//========================================================================================
void write_Recfast_output(string oname, vector<double> &zarr, vector<double> &Xe_H, 
                          vector<double> &Xe_He, vector<double> &Xe, vector<double> &TM)
{
    if(verbose) 
        cout << " Recfast++ :: write_Recfast_output : writing output into: " 
             << oname << endl;

    ofstream ofile(oname.c_str());
    ofile.precision(8);
    
    for(int i=0; i<(int)zarr.size(); i++)
        ofile << zarr[i] << " " << Xe[i] << " " << Xe_H[i] << " " 
              << Xe_He[i] << " " << TM[i] << endl;
        
    return;
}

//========================================================================================
// main program
//========================================================================================
int main(int narg, char *args[])
{
    //====================================================================================
    // System call to clear screen 
    //====================================================================================
    system("clear");
    
    if(narg < 2) 
    {
        cout << " Recfast-error: usage $ ./Recfast++ parameters.dat" << endl;
        exit(0);
    }
    
    string fname=args[narg-1];

    //====================================================================================
    // reading parameters
    //====================================================================================
    vector<double> params(14);
    int npz=10000;
    double zstart=1.0e+4, zend=0.001;
    
    params[0]=npz;
    params[1]=zstart;
    params[2]=zend;
    
    read_parameter_file(fname, params);
    
    //====================================================================================
    // running Recfast part
    //====================================================================================
    if(verbose) cout << " Recfast++ :: Entering computation. " << endl;

    vector<double> zarr(npz), Xe_H(npz), Xe_He(npz), Xe(npz), TM(npz);
    
    //====================================================================================
    // main code
    //====================================================================================
    Xe_frac(&params[0], &zarr[0], &Xe_H[0], &Xe_He[0], &Xe[0], &TM[0]);
    
    string oname="./output/Xe_Recfast";
    if(params[12]!=0.0) oname+=".DM_annihilations";
    if(params[13]==1) oname+=".Rec_corrs_CT2010";
    oname+=".dat";
    
    write_Recfast_output(oname, zarr, Xe_H, Xe_He, Xe, TM);

    if(verbose) cout << " Recfast++ :: Run completed. Exiting. " << endl;

    return 0;
}
