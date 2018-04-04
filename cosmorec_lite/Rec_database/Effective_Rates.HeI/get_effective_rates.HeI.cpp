//======================================================================
// get_effective_rates.HeI.cpp
//
// Created by Jens Chluba (Oct 2010).
// Copyright 2010 CITA (U of T). All rights reserved.
//
//======================================================================

//======================================================================
// Standards
//======================================================================
#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>
#include <iomanip>
#include <cmath>
#include <limits.h>
#include <vector>

//======================================================================
// Relevant lib files
//======================================================================
#include "HeI_Atom.h"
#include "physical_consts.h"
#include "routines.h"
#include "get_effective_rates.HeI.h"

//======================================================================
//======================================================================

struct res_state_Data_HeI
{
    int HeI_index;
    
    vector<double> BiVec;
    vector<double> BitotVec;
    vector<vector<double> > RijMatrix;
    vector<double> AiVec;
    
    // Tg related quantities
    double lgTgmin, lgTgmax; 
    int Npoints_Tg;
    
    // # of resolved states 
    int neqres;
    unsigned long row, col; // search indices; accelerate the process
    
    vector<double> lgTgarr;
};

vector<res_state_Data_HeI> res_state_D_HeI;
vector<int> resolved_level_map_HeI;

//======================================================================
//======================================================================

//======================================================================
// main functions
//======================================================================
void read_parameter_file(string fname, Gas_of_HeI_Atoms &HeIA, vector<int> &res_Level_indices)
{
    int nShells = HeIA.Get_nShells();
    
    res_Level_indices.clear();
    
    ifstream ifile(fname.c_str());
    if(!ifile.is_open()){ cerr << " read_parameter_file (HeI) :: File : " << fname << " does not exist.. exiting " << endl; exit(0); }
    
    int n, l, s, j;
    
    while(ifile)
    {
        ifile >> n;
        ifile >> l;
        ifile >> s;
        ifile >> j;
        
        if(n > nShells)
        {
            cout << " nShells is smaller than resolved states asked for" << endl;
            cout<< " Skipping levels with n> " << nShells << endl;
            break;
        }
        if(ifile)
        {   
            res_Level_indices.push_back(HeIA.Get_Level_index(n, l, s, j));
            cout << " resolved level # " << res_Level_indices.size()-1 << " :: " << n << " " << l << " s= " << s << " j= " << j << " HeI_index= " << res_Level_indices.back() << endl;
        }
    }
    
    ifile.close();
    return;
}

//======================================================================
void load_rates(string effective_rate_path, int nShell_max, int index_of_res_state, Gas_of_HeI_Atoms &HeIA, res_state_Data_HeI &Data_res)
{   
    Data_res.BiVec.clear();
    Data_res.BitotVec.clear();
    Data_res.AiVec.clear();
    Data_res.RijMatrix.clear();
    Data_res.lgTgarr.clear();
    
    // Open file to read 
    
    int nres = HeIA.Get_n(index_of_res_state);
    int lres = HeIA.Get_l(index_of_res_state);
    int sres = HeIA.Get_S(index_of_res_state);
    int jres = HeIA.Get_J(index_of_res_state);
    double dum;
    
    string filenamei = effective_rate_path + "HeI_Rates_n"+int_to_string(nres)+"_l"+int_to_string(lres)+"_S"
                        +int_to_string(sres)+"_J"+int_to_string(jres)+".nS_"+int_to_string(nShell_max)+".dat";
    
    ifstream  Ratesfile(filenamei.c_str());
    if(Ratesfile.is_open())
    {
        Ratesfile >> Data_res.Npoints_Tg;
        Ratesfile >> Data_res.neqres;

        for(int m=0; m<Data_res.Npoints_Tg; m++){ Ratesfile >> dum; Data_res.lgTgarr.push_back(dum); } 
        
        Data_res.lgTgmin=Data_res.lgTgarr.front();
        Data_res.lgTgmax=Data_res.lgTgarr.back();
        Data_res.row=Data_res.col=0;
    }
    else
    { cerr << " load_rates::File : " << filenamei << " does not exist.. exiting \n" << endl; exit(0); }
    
    // ---------------Data is stored as follows----------------------
    // lines 1 to 2 is for meta-data. Following that every line (row)  
    // corresponts to a particular Tg. Column is ordered as Bi, Rijs,
    // and ALL Ai. 
    // --------------------------------------------------------------
    vector<double> dumVecRij(Data_res.neqres);
    
    for(int i=0; i<Data_res.Npoints_Tg; i++){
        
        // Reading Bi
        Ratesfile >> dum;
        Data_res.BiVec.push_back(dum);
        
        // Reading Bi_tot
        Ratesfile >> dum;
        Data_res.BitotVec.push_back(dum);

        // Reading Rij 
        for(int j=0; j<Data_res.neqres; j++) Ratesfile >> dumVecRij[j];
        Data_res.RijMatrix.push_back(dumVecRij);
            
        // Reading Ai
        Ratesfile >> dum;
        Data_res.AiVec.push_back(dum);
    }
    
    Ratesfile.close();
        
    return;
}

void load_rates(string path, int nShell_max, Gas_of_HeI_Atoms &HeIA)
{
    string effective_rate_path=path;
    string fname=effective_rate_path + "res_state_list.dat";
    vector<int> res_Level_indices;
    res_state_D_HeI.clear();
    resolved_level_map_HeI.clear();
    
    cout << "\n load_rates::Loading effective Rate data (HeI) " << endl;
    
    // read configuration from setup-file
    read_parameter_file(fname, HeIA, res_Level_indices);
    
    res_state_Data_HeI dum;
    int nres=res_Level_indices.size();
    
    for(int i=0; i<(int)res_Level_indices.size(); i++)
    {
        dum.HeI_index=res_Level_indices[i];
        load_rates(effective_rate_path, nShell_max, dum.HeI_index, HeIA, dum);      
        res_state_D_HeI.push_back(dum);
    }

    // create level map
    int ntl=HeIA.Get_total_number_of_Levels();
    for(int lev=0; lev<ntl; lev++) resolved_level_map_HeI.push_back(nres);
    for(int i=0; i<(int)res_Level_indices.size(); i++) resolved_level_map_HeI[get_HeI_index(i)]=i;
        
    cout << " load_rates::finished. # of resolved states: " << res_state_D_HeI.size() 
         << " Tg-points: " << res_state_D_HeI[0].lgTgarr.size() << endl << endl;
    
    return;
}


//========================================================================================
// functions for mapping between HI-level and res-state indicies
//========================================================================================
int get_number_of_resolved_levels_HeI(){ return res_state_D_HeI.size(); }

int get_HeI_index(int index_of_res_state)
{ return res_state_D_HeI[index_of_res_state].HeI_index; }

//--------------------------------------------------------------------
// for resolved states this function returns the resolved level index;
// for unresolved states this function returns index=nres
//--------------------------------------------------------------------
int get_res_index_HeI(int index_of_HeI_state){ return resolved_level_map_HeI[index_of_HeI_state]; }

//========================================================================================
// compute effective rates
//========================================================================================
void get_rates_HeI(double Tg, vector<double> &Ai, vector<double> &Bi, vector<vector<double> > &RijVec)
{
    if(res_state_D_HeI[0].BiVec.size()==0){ cerr << " get_rates_HeI:: It seems like the transition data has not been loaded. Exiting... " << endl; exit(0); }
    if(Ai.size()!=res_state_D_HeI.size()){ cerr << " get_rates_HeI:: It seems like the Ai, Bi and RijVec vectors were not set. Exiting... " << endl; exit(0); } 
    
    double logTg  = log(Tg);
        
    // Make sure all points are in the table including the buffer for interpolation
    if( logTg < res_state_D_HeI[0].lgTgmin || logTg > res_state_D_HeI[0].lgTgmax )
    {
        cerr << " get_rates_HeI:: Tg = " << Tg << " is out of table range precomputed." << endl;    
        exit(0);
    }
    
    locate_JC(&res_state_D_HeI[0].lgTgarr[0], res_state_D_HeI[0].Npoints_Tg, logTg, &res_state_D_HeI[0].row); 

    int lx = (int) res_state_D_HeI[0].row;    // index of x coordinate of adjacent grid point to left of P

    // ---------------------------------------------
    // Boundary condition                           
    // ---------------------------------------------
    if(lx+3 > res_state_D_HeI[0].Npoints_Tg || lx<0)
    { cerr << " get_rates_HeI:: out of bounds for Tg = " << Tg << ". Exiting! " << endl; exit(0); } 
    
    double x = logTg;
    double x1 = res_state_D_HeI[0].lgTgarr[lx];
    double x2 = res_state_D_HeI[0].lgTgarr[lx+1];
    double x3 = res_state_D_HeI[0].lgTgarr[lx+2];
    double x4 = res_state_D_HeI[0].lgTgarr[lx+3];
    
    // ------------------------------------------- 
    // Coefficients for 1-D interpolations         
    // ------------------------------------------- 
    double a[4];
    //
    a[0]= (x - x2)*(x - x3)*(x - x4)/((x1 - x2)*(x1 - x3)*(x1 - x4));
    a[1]= (x - x1)*(x - x3)*(x - x4)/((x1 - x2)*(x2 - x3)*(x4 - x2));
    a[2]= (x - x1)*(x - x2)*(x - x4)/((x1 - x3)*(x2 - x3)*(x3 - x4));       
    a[3]= (x - x1)*(x - x2)*(x - x3)/((x3 - x4)*(x4 - x1)*(x2 - x4));
    
    double fx;  
    for(int m=0; m<(int)res_state_D_HeI.size(); m++)
    {
        // ------------------------------------- 
        // Simple 1-D interpolation to find Ai   
        // -------------------------------------        
        fx=0.0;
        for(int k=0; k<4; k++) fx+= a[k]*(res_state_D_HeI[m].AiVec[k+lx]-res_state_D_HeI[m].BiVec[k+lx]);
        
        Ai[m] = exp(fx);
        
        // ------------------------------------- 
        // Simple 1-D interpolation to find Bi   
        // -------------------------------------        
        fx=0.0;
        for(int k=0; k<4; k++) fx+= a[k]*res_state_D_HeI[m].BiVec[k+lx];
        
        Bi[m] = exp(fx);
        
        // ------------------------------------- 
        // Getting all Ri->j interpolated values 
        // ------------------------------------- 
        for(int i=0; i<=m; i++) RijVec[m][i] = 0.0; // for these the detailed balance relations inside the code will be used
        
        for(int i=m+1; i<res_state_D_HeI[m].neqres; i++)
        {
            fx=0.0;
            for(int k=0; k<4; k++) fx+= a[k]*res_state_D_HeI[m].RijMatrix[k+lx][i];
            
            RijVec[m][i] = exp(fx);
        }
    }
    
    return;
}

//========================================================================================
void get_rates_all_HeI(double Tg, vector<double> &Ai, vector<double> &Bi, vector<double> &Bitot, vector<vector<double> > &RijVec)
{
    if(res_state_D_HeI[0].BiVec.size()==0){ cerr << " get_rates_HeI:: It seems like the transition data has not been loaded. Exiting... " << endl; exit(0); }
    if(Ai.size()!=res_state_D_HeI.size()){ cerr << " get_rates_HeI:: It seems like the Ai, Bi and RijVec vectors were not set. Exiting... " << endl; exit(0); } 
    
    double logTg  = log(Tg);
    
    // Make sure all points are in the table including the buffer for interpolation
    if( logTg < res_state_D_HeI[0].lgTgmin || logTg > res_state_D_HeI[0].lgTgmax )
    {
        cerr << " get_rates_HeI:: Tg = " << Tg << " is out of table range precomputed." << endl;    
        exit(0);
    }
    
    locate_JC(&res_state_D_HeI[0].lgTgarr[0], res_state_D_HeI[0].Npoints_Tg, logTg, &res_state_D_HeI[0].row); 
    
    int lx = (int) res_state_D_HeI[0].row;    // index of x coordinate of adjacent grid point to left of P
    
    // ---------------------------------------------
    // Boundary condition                           
    // ---------------------------------------------
    if(lx+3 > res_state_D_HeI[0].Npoints_Tg || lx<0)
    { cerr << " get_rates_HeI:: out of bounds for Tg = " << Tg << ". Exiting! " << endl; exit(0); } 
    
    double x = logTg;
    double x1 = res_state_D_HeI[0].lgTgarr[lx];
    double x2 = res_state_D_HeI[0].lgTgarr[lx+1];
    double x3 = res_state_D_HeI[0].lgTgarr[lx+2];
    double x4 = res_state_D_HeI[0].lgTgarr[lx+3];
    
    // ------------------------------------------- 
    // Coefficients for 1-D interpolations         
    // ------------------------------------------- 
    double a[4];
    //
    a[0]= (x - x2)*(x - x3)*(x - x4)/((x1 - x2)*(x1 - x3)*(x1 - x4));
    a[1]= (x - x1)*(x - x3)*(x - x4)/((x1 - x2)*(x2 - x3)*(x4 - x2));
    a[2]= (x - x1)*(x - x2)*(x - x4)/((x1 - x3)*(x2 - x3)*(x3 - x4));       
    a[3]= (x - x1)*(x - x2)*(x - x3)/((x3 - x4)*(x4 - x1)*(x2 - x4));
    
    double fx;  
    for(int m=0; m<(int)res_state_D_HeI.size(); m++)
    {
        // --------------------------------------- 
        // Simple 1-D interpolation to find Ai   
        // --------------------------------------- 
        fx=0.0;
        for(int k=0; k<4; k++) fx+= a[k]*(res_state_D_HeI[m].AiVec[k+lx]-res_state_D_HeI[m].BiVec[k+lx]);
        
        Ai[m] = exp(fx);
        
        // --------------------------------------- 
        // Simple 1-D interpolation to find Bi   
        // --------------------------------------- 
        fx=0.0;
        for(int k=0; k<4; k++) fx+= a[k]*res_state_D_HeI[m].BiVec[k+lx];
        
        Bi[m] = exp(fx);
        
        // --------------------------------------- 
        // Simple 1-D interpolation to find Bitot   
        // ---------------------------------------      
        fx=0.0;
        for(int k=0; k<4; k++) fx+= a[k]*res_state_D_HeI[m].BitotVec[k+lx];
        
        Bitot[m] = exp(fx);
        
        // --------------------------------------- 
        // Getting all Ri->j interpolated values 
        // --------------------------------------- 
        for(int i=0; i<=m; i++) RijVec[m][i] = 0.0; // for these the detailed balance relations inside the code will be used
        
        for(int i=m+1; i<res_state_D_HeI[m].neqres; i++)
        {
            fx=0.0;
            for(int k=0; k<4; k++) fx+= a[k]*res_state_D_HeI[m].RijMatrix[k+lx][i];
            
            RijVec[m][i] = exp(fx);
        }
    }
    
    return;
}

//========================================================================================
// compute effective rate out of level
//========================================================================================
void get_Bitot_HeI(double Tg, vector<double> &Bitot)
{
    if(res_state_D_HeI[0].BiVec.size()==0){ cerr << " get_Bitot_HeI:: It seems like the transition data has not been loaded. Exiting... " << endl; exit(0); }
    if(Bitot.size()!=res_state_D_HeI.size()){ cerr << " get_Bitot_HeI:: It seems like the Bitot vector was not set. Exiting... " << endl; exit(0); }    
    
    double logTg  = log(Tg);
    
    // Make sure all points are in the table including the buffer for interpolation
    if( logTg < res_state_D_HeI[0].lgTgmin || logTg > res_state_D_HeI[0].lgTgmax )
    {
        cerr << " get_Bitot_HeI:: Tg = " << Tg << " is out of table range precomputed." << endl;    
        exit(0);
    }
    
    locate_JC(&res_state_D_HeI[0].lgTgarr[0], res_state_D_HeI[0].Npoints_Tg, logTg, &res_state_D_HeI[0].row); 
    
    int lx = (int) res_state_D_HeI[0].row;    // index of x coordinate of adjacent grid point to left of P
    
    // ---------------------------------------------
    // Boundary condition                           
    // ---------------------------------------------
    if(lx+3 > res_state_D_HeI[0].Npoints_Tg || lx<0)
    { cerr << " get_Bitot_HeI:: out of bounds for Tg = " << Tg << ". Exiting! " << endl; exit(0); } 
    
    double x = logTg;
    double x1 = res_state_D_HeI[0].lgTgarr[lx];
    double x2 = res_state_D_HeI[0].lgTgarr[lx+1];
    double x3 = res_state_D_HeI[0].lgTgarr[lx+2];
    double x4 = res_state_D_HeI[0].lgTgarr[lx+3];
    
    // ------------------------------------------- 
    // Coefficients for 1-D interpolations         
    // ------------------------------------------- 
    double a[4];
    //
    a[0]= (x - x2)*(x - x3)*(x - x4)/((x1 - x2)*(x1 - x3)*(x1 - x4));
    a[1]= (x - x1)*(x - x3)*(x - x4)/((x1 - x2)*(x2 - x3)*(x4 - x2));
    a[2]= (x - x1)*(x - x2)*(x - x4)/((x1 - x3)*(x2 - x3)*(x3 - x4));       
    a[3]= (x - x1)*(x - x2)*(x - x3)/((x3 - x4)*(x4 - x1)*(x2 - x4));
    
    double fx;  
    for(int m=0; m<(int)res_state_D_HeI.size(); m++)
    {
        // --------------------------------------- 
        // Simple 1-D interpolation to find Bitot   
        // ---------------------------------------      
        fx=0.0;
        for(int k=0; k<4; k++) fx+= a[k]*res_state_D_HeI[m].BitotVec[k+lx];
        
        Bitot[m] = exp(fx);
    }
    
    return;
}

//========================================================================================
// compute effective rate out of level
//========================================================================================
void get_Bitot_HeI(double Tg, double &Bitot, int ires)
{
    if(res_state_D_HeI[0].BiVec.size()==0){ cerr << " get_Bitot_HeI:: It seems like the transition data has not been loaded. Exiting... " << endl; exit(0); }
    if(ires>=(int)res_state_D_HeI.size() || ires < 0){ cerr << " get_Bitot_HeI:: It seems like you are asking for a non-existing level. Exiting... " << endl; exit(0); }    
    
    double logTg  = log(Tg);
    
    // Make sure all points are in the table including the buffer for interpolation
    if( logTg < res_state_D_HeI[0].lgTgmin || logTg > res_state_D_HeI[0].lgTgmax )
    {
        cerr << " get_Bitot_HeI:: Tg = " << Tg << " is out of table range precomputed." << endl;    
        exit(0);
    }
    
    locate_JC(&res_state_D_HeI[0].lgTgarr[0], res_state_D_HeI[0].Npoints_Tg, logTg, &res_state_D_HeI[0].row); 
    
    int lx = (int) res_state_D_HeI[0].row;    // index of x coordinate of adjacent grid point to left of P
    
    // ---------------------------------------------
    // Boundary condition                           
    // ---------------------------------------------
    if(lx+3 > res_state_D_HeI[0].Npoints_Tg || lx<0)
    { cerr << " get_Bitot_HeI:: out of bounds for Tg = " << Tg << ". Exiting! " << endl; exit(0); } 
    
    double x = logTg;
    double x1 = res_state_D_HeI[0].lgTgarr[lx];
    double x2 = res_state_D_HeI[0].lgTgarr[lx+1];
    double x3 = res_state_D_HeI[0].lgTgarr[lx+2];
    double x4 = res_state_D_HeI[0].lgTgarr[lx+3];
    
    // ------------------------------------------- 
    // Coefficients for 1-D interpolations         
    // ------------------------------------------- 
    double a[4];
    //
    a[0]= (x - x2)*(x - x3)*(x - x4)/((x1 - x2)*(x1 - x3)*(x1 - x4));
    a[1]= (x - x1)*(x - x3)*(x - x4)/((x1 - x2)*(x2 - x3)*(x4 - x2));
    a[2]= (x - x1)*(x - x2)*(x - x4)/((x1 - x3)*(x2 - x3)*(x3 - x4));       
    a[3]= (x - x1)*(x - x2)*(x - x3)/((x3 - x4)*(x4 - x1)*(x2 - x4));
    
    // --------------------------------------- 
    // Simple 1-D interpolation to find Bitot   
    // ---------------------------------------      
    double fx=0.0;
    for(int k=0; k<4; k++) fx+= a[k]*res_state_D_HeI[ires].BitotVec[k+lx];
    
    Bitot = exp(fx);
    
    return;
}

