//========================================================================================
// Author: Jens Chluba
// Last modification: Oct 2010
// CITA, University of Toronto
// All rights reserved.
//========================================================================================

//========================================================================================
// This module allows to add the recombination corrections as described by 
// Chluba & Thomas 2010. Simple linear interpolation of DNe/Ne is used, as explained by
// Rubino-Martin et al. 2010.
//========================================================================================

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

#include "constants.Recfast.h"
#include "Rec_corrs_CT.Recfast.h"

using namespace std;
using namespace RECFAST_physical_constants; // defined in constants.Recfast.h

//========================================================================================
// data-file with correction
//========================================================================================
string Rec_corrs_CT_file=(string)RECFASTPPPATH+"./src/Data/DXe_Xe.CT2010.dat";

//========================================================================================
// some global variables for compution
//========================================================================================
int Rec_corrs_CT_start_index=0;
vector<vector<double> > Rec_corrs_CT_Data;

//========================================================================================
// reading the correction factor
//========================================================================================
void read_recombination_correction(int mflag)
{
    Rec_corrs_CT_start_index=0;
    if(Rec_corrs_CT_Data.size()>0) return;
    
    ifstream ifile(Rec_corrs_CT_file.c_str());
    if(!ifile)
    { 
        cerr << " Recfast++ :: read_recombination_correction:" 
             << " file with correction factor ( " << Rec_corrs_CT_file 
             << " ) not found. Exiting. " << endl; 
             exit(0); 
    }
    
    vector<double> dum(2);
    Rec_corrs_CT_Data.clear();
    
    if(mflag>=1) 
        cout << " Recfast++ :: read_recombination_correction:"
             << " reading recombination correction according to Chluba & Thomas 2010 " 
             << endl;
    
    do
    {
        ifile >> dum[0];
        ifile >> dum[1];
        
        if(ifile) Rec_corrs_CT_Data.push_back(dum);
    }
    while(ifile);

    return;
}

//========================================================================================
double interpolate_data(double z)
{
    //====================================================================================
    // find next k index with z>z(k)
    //====================================================================================
    if(z<=Rec_corrs_CT_Data[Rec_corrs_CT_start_index][0]) Rec_corrs_CT_start_index++;
    
    //====================================================================================
    // simple linear interpolation
    //====================================================================================
    double dz=z-Rec_corrs_CT_Data[Rec_corrs_CT_start_index-1][0];
    double df_dz=(Rec_corrs_CT_Data[Rec_corrs_CT_start_index][1]
                 -Rec_corrs_CT_Data[Rec_corrs_CT_start_index-1][1])
                /(Rec_corrs_CT_Data[Rec_corrs_CT_start_index][0]
                 -Rec_corrs_CT_Data[Rec_corrs_CT_start_index-1][0]);
    
    double fcorr=Rec_corrs_CT_Data[Rec_corrs_CT_start_index-1][1]+df_dz*dz; 
    return 1.0+fcorr;
}

//========================================================================================
double Recombination_correction_factor(double z)
{
    int nz=Rec_corrs_CT_Data.size();
    
    if(Rec_corrs_CT_Data[0][0]<=z) return 1.0;
    if(Rec_corrs_CT_Data[nz-1][0]>=z) return 1.0+Rec_corrs_CT_Data[nz-1][1];
    return interpolate_data(z);
}
