//======================================================================
// get_effective_rates.HI.h
//
//
// Created by Rajat Thomas and Jens Chluba on 10-07-20.
// Copyright 2010 CITA ( U of T). All rights reserved.
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

//========================================================================================
// Relevant lib files
//========================================================================================
#include "Atom.h"
#include "physical_consts.h"
#include "routines.h"

#include "get_effective_rates.HI.h"

//========================================================================================
// global variables for computations
//========================================================================================
struct res_state_Data
{
    int HI_index;
    
    vector<double> BiVec;
    vector<double> BitotVec;
    vector<vector<double> > RijMatrix;
    vector<vector<double> > AiMatrix;
    
    /* Tg related quantities */
    double lgTgmin, lgTgmax; 
    int Npoints_Tg;
    
    /* Te related quantities */
    double lgrhomin, lgrhomax; 
    int Npoints_rho;
    
    /* # of resolved states */
    int neqres;
    unsigned long row, col; // search indices; accelerate the process
    
    vector<double> lgTgarr;
    vector<double> lgrhoarr;
};


vector<res_state_Data> res_state_D;
vector<int> resolved_level_map;

//========================================================================================
// loading the rates
//========================================================================================
void read_parameter_file(string fname, Gas_of_Atoms &HIA, vector<int> &res_Level_indices)
{
    int nShells = HIA.Get_nShells();
    
    res_Level_indices.clear();
    
    ifstream ifile(fname.c_str());
    if(!ifile.is_open()){ cerr << " read_parameter_file (HI) :: File : " << fname << " does not exist.. exiting " << endl; exit(0); }
    
    int n, l;
    
    while(ifile)
    {
        ifile >> n;
        ifile >> l;
        
        if(n > nShells)
        {
            cout << " nShells is smaller than resolved states asked for" << endl;
            cout<< " Skipping levels with n> " << nShells << endl;
            break;
        }
        if(ifile)
        {   
            res_Level_indices.push_back(HIA.Get_Level_index(n, l));
            cout << " resolved level # " << res_Level_indices.size()-1 << " :: " << n << " " << l << " HI_index= " << res_Level_indices.back() << endl;
        }
    }
    
    ifile.close();
    return;
}



void load_rates(string effective_rate_path, int nShell_max, int index_of_res_state, Gas_of_Atoms &HIA, res_state_Data &Data_res)
{   
    Data_res.BiVec.clear();
    Data_res.BitotVec.clear();
    Data_res.AiMatrix.clear();
    Data_res.RijMatrix.clear();
    Data_res.lgTgarr.clear();
    Data_res.lgrhoarr.clear();
    
    /* Open file to read */
    
    int nres = HIA.Get_n_of_Level(index_of_res_state);
    int lres = HIA.Get_l_of_Level(index_of_res_state);
    double dum;
    
    string filenamei = effective_rate_path + "Rates_n"+int_to_string(nres)+"_l"+int_to_string(lres)+".nS_"+int_to_string(nShell_max)+".dat";
    
    ifstream  Ratesfile(filenamei.c_str());
    if(Ratesfile.is_open())
    {
        Ratesfile >> Data_res.Npoints_Tg;
        Ratesfile >> Data_res.Npoints_rho;
        Ratesfile >> Data_res.neqres;

        for(int m=0; m<Data_res.Npoints_Tg; m++){ Ratesfile >> dum; Data_res.lgTgarr.push_back(dum); } 
        for(int m=0; m<Data_res.Npoints_rho; m++){ Ratesfile >> dum; Data_res.lgrhoarr.push_back(dum); }
        
        Data_res.lgTgmin=Data_res.lgTgarr.front();
        Data_res.lgTgmax=Data_res.lgTgarr.back();
        Data_res.lgrhomin=Data_res.lgrhoarr.front();
        Data_res.lgrhomax=Data_res.lgrhoarr.back();
        Data_res.row=Data_res.col=0;
    }
    else
    { cerr << " load_rates::File : " << filenamei << " does not exist.. exiting \n" << endl; exit(0); }
    
    /* ---------------Data is stores as follows----------------------
     lines 1 to 3 is for meta-data. Following that every line (row)  
     corresponts to a particular Tg. Column is ordered as Bi, Rijs,
     and ALL Ai for each Te value. 
     So the # cols = 1(Bi) + neqres(Rij) + Npoints_rho(Ai)
     # rows = 3(meta-data+ + Npoints_Tg
     --------------------------------------------------------------*/
    vector<double> dumVecRij(Data_res.neqres);
    vector<double> dumVecAi(Data_res.Npoints_rho);

    for(int i=0; i<Data_res.Npoints_Tg; i++){
        
        /* Reading Bi */
        Ratesfile >> dum;
        Data_res.BiVec.push_back(dum);
        
        /* Reading Bi_tot */
        Ratesfile >> dum;
        Data_res.BitotVec.push_back(dum);
        
        /* Reading Rij */
        for(int j=0; j<Data_res.neqres; j++) Ratesfile >> dumVecRij[j];
        Data_res.RijMatrix.push_back(dumVecRij);
            
        /* Reading Ai */
        for(int j=0; j<Data_res.Npoints_rho; j++) Ratesfile >> dumVecAi[j];
        Data_res.AiMatrix.push_back(dumVecAi);
    }
    
    Ratesfile.close();
    
    return;
}

void load_rates(string path, int nShell_max, Gas_of_Atoms &HIA)
{
    string effective_rate_path=path;
    string fname=effective_rate_path + "res_state_list.dat";
    vector<int> res_Level_indices;
    res_state_D.clear();
    resolved_level_map.clear();
    
    cout << "\n load_rates::Loading effective Rate data (HI) " << endl;
    
    // read configuration from setup-file
    read_parameter_file(fname, HIA, res_Level_indices);
    
    res_state_Data dum;
    int nres=res_Level_indices.size();
    for(int i=0; i<nres; i++)
    {
        dum.HI_index=res_Level_indices[i];
        load_rates(effective_rate_path, nShell_max, dum.HI_index, HIA, dum);        
        res_state_D.push_back(dum);
    }
    
    // create level map
    int ntl=HIA.Get_total_number_of_Levels();
    for(int lev=0; lev<ntl; lev++) resolved_level_map.push_back(nres);
    for(int i=0; i<(int)res_Level_indices.size(); i++) resolved_level_map[get_HI_index(i)]=i;

    cout << " load_rates::finished. # of resolved states: " << res_state_D.size() 
         << " Tg-points: " << res_state_D[0].lgTgarr.size() << endl << endl;
    
    return;
}

//========================================================================================
// functions for mapping between HI-level and res-state indicies
//========================================================================================
int get_number_of_resolved_levels(){ return res_state_D.size(); }

int get_HI_index(int index_of_res_state)
{ return res_state_D[index_of_res_state].HI_index; }

//--------------------------------------------------------------------
// for resolved states this function returns the resolved level index;
// for unresolved states this function returns index=nres
//--------------------------------------------------------------------
int get_res_index(int index_of_HI_state){ return resolved_level_map[index_of_HI_state]; }

//========================================================================================
// computations of rates
//========================================================================================
void get_rates(double Tg, double Te, vector<double> &Ai, vector<double> &Bi, vector<vector<double> > &RijVec)
{
    if(res_state_D[0].BiVec.size()==0){ cerr << " get_rates:: It seems like the transition data has not been loaded. Exiting... " << endl; exit(0); }
    if(Ai.size()!=res_state_D.size()){ cerr << " get_rates:: It seems like the Ai, Bi and RijVec vectors were not set. Exiting... " << endl; exit(0); } 
    
    double logTg  = log(Tg);
    double logrho  = log(Te/Tg);
        
    /* Make sure all points are in the table including the buffer for interpolation */
    if( logTg < res_state_D[0].lgTgmin || logTg > res_state_D[0].lgTgmax || logrho < res_state_D[0].lgrhomin || logrho > res_state_D[0].lgrhomax )
    {
        cerr << " get_rates:: The pair of (Tg,Te) = (" << Tg << "," << Te << ") == Te/Tg = " << Te/Tg << " is out of table range precomputed." << endl  ;   
        exit(0);
    }
    
    locate_JC(&res_state_D[0].lgTgarr[0], res_state_D[0].Npoints_Tg, logTg, &res_state_D[0].row); 
    locate_JC(&res_state_D[0].lgrhoarr[0], res_state_D[0].Npoints_rho, logrho, &res_state_D[0].col); 

    int lx = (int) res_state_D[0].row;    // index of x coordinate of adjacent grid point to left of P
    int ly = (int) res_state_D[0].col;    // index of y coordinate of adjacent grid point below P

    /* --------------------------------------------- */
    /* Boundary condition                            */
    /* --------------------------------------------- */
    if(lx+3 > res_state_D[0].Npoints_Tg || ly+3 > res_state_D[0].Npoints_rho || lx<0 || ly<0)
    { cerr << " get_rates:: out of bounds for Tg = " << Tg << " Te= " << Te << ". Exiting! " << endl; exit(0); } 
    
    double x = logTg;
    double y = logrho;
    double x1 = res_state_D[0].lgTgarr[lx];
    double x2 = res_state_D[0].lgTgarr[lx+1];
    double x3 = res_state_D[0].lgTgarr[lx+2];
    double x4 = res_state_D[0].lgTgarr[lx+3];
    double y1 = res_state_D[0].lgrhoarr[ly];
    double y2 = res_state_D[0].lgrhoarr[ly+1];
    double y3 = res_state_D[0].lgrhoarr[ly+2];
    double y4 = res_state_D[0].lgrhoarr[ly+3];
    
    /* ------------------------------------------- */
    /* Coefficients for 1-D interpolations         */
    /* ------------------------------------------- */
    double a[4], b[4];
    //
    a[0]= (x - x2)*(x - x3)*(x - x4)/((x1 - x2)*(x1 - x3)*(x1 - x4));
    a[1]= (x - x1)*(x - x3)*(x - x4)/((x1 - x2)*(x2 - x3)*(x4 - x2));
    a[2]= (x - x1)*(x - x2)*(x - x4)/((x1 - x3)*(x2 - x3)*(x3 - x4));       
    a[3]= (x - x1)*(x - x2)*(x - x3)/((x3 - x4)*(x4 - x1)*(x2 - x4));
    //
    b[0]= (y - y2)*(y - y3)*(y - y4)/((y1 - y2)*(y1 - y3)*(y1 - y4));
    b[1]= (y - y1)*(y - y3)*(y - y4)/((y1 - y2)*(y2 - y3)*(y4 - y2));
    b[2]= (y - y1)*(y - y2)*(y - y4)/((y1 - y3)*(y2 - y3)*(y3 - y4));       
    b[3]= (y - y1)*(y - y2)*(y - y3)/((y3 - y4)*(y4 - y1)*(y2 - y4));
    
    double fxy, fx; 
    for(int m=0; m<(int)res_state_D.size(); m++)
    {
        /* ------------------------------------- */
        /* The 2-D interpolated valued of Ai     */
        /* ------------------------------------- */
        fxy=0.0;
        for(int i=0; i<4; i++)
        {
            fx=0.0;
            for(int k=0; k<4; k++) fx+=a[k]*(res_state_D[m].AiMatrix[k+lx][i+ly]-res_state_D[m].BiVec[k+lx]); 
            fxy+=b[i]*fx;                                                          // interpolation in y 
        }
        
        Ai[m] = exp(fxy);

        /* ------------------------------------- */
        /* Simple 1-D interpolation to find Bi   */
        /* ------------------------------------- */     
        fx=0.0;
        for(int k=0; k<4; k++) fx+= a[k]*res_state_D[m].BiVec[k+lx];
        
        Bi[m] = exp(fx);

        /* ------------------------------------- */
        /* Getting all Ri->j interpolated values */
        /* ------------------------------------- */
        for(int i=0; i<=m; i++) RijVec[m][i] = 0.0; // for these the detailed balance relations inside the code will be used
        
        for(int i=m+1; i<res_state_D[m].neqres; i++)
        {
            fx=0.0;
            for(int k=0; k<4; k++) fx+= a[k]*res_state_D[m].RijMatrix[k+lx][i];
            
            RijVec[m][i] = exp(fx);
        }
    }
    
    return;
}

//========================================================================================
void get_rates_all(double Tg, double Te, vector<double> &Ai, vector<double> &Bi, vector<double> &Bitot, vector<vector<double> > &RijVec)
{
    if(res_state_D[0].BiVec.size()==0){ cerr << " get_rates_all:: It seems like the transition data has not been loaded. Exiting... " << endl; exit(0); }
    if(Ai.size()!=res_state_D.size()){ cerr << " get_rates_all:: It seems like the Ai, Bi and RijVec vectors were not set. Exiting... " << endl; exit(0); } 
    
    double logTg  = log(Tg);
    double logrho  = log(Te/Tg);
    
    /* Make sure all points are in the table including the buffer for interpolation */
    if( logTg < res_state_D[0].lgTgmin || logTg > res_state_D[0].lgTgmax || logrho < res_state_D[0].lgrhomin || logrho > res_state_D[0].lgrhomax )
    {
        cerr << " get_rates_all:: The pair of (Tg,Te) = (" << Tg << "," << Te << ") == Te/Tg = " << Te/Tg << " is out of table range precomputed." << endl  ;   
        exit(0);
    }
    
    locate_JC(&res_state_D[0].lgTgarr[0], res_state_D[0].Npoints_Tg, logTg, &res_state_D[0].row); 
    locate_JC(&res_state_D[0].lgrhoarr[0], res_state_D[0].Npoints_rho, logrho, &res_state_D[0].col); 
    
    int lx = (int) res_state_D[0].row;    // index of x coordinate of adjacent grid point to left of P
    int ly = (int) res_state_D[0].col;    // index of y coordinate of adjacent grid point below P
    
    /* --------------------------------------------- */
    /* Boundary condition                            */
    /* --------------------------------------------- */
    if(lx+3 > res_state_D[0].Npoints_Tg || ly+3 > res_state_D[0].Npoints_rho || lx<0 || ly<0)
    { cerr << " get_rates_all:: out of bounds for Tg = " << Tg << " Te= " << Te << ". Exiting! " << endl; exit(0); } 
    
    double x = logTg;
    double y = logrho;
    double x1 = res_state_D[0].lgTgarr[lx];
    double x2 = res_state_D[0].lgTgarr[lx+1];
    double x3 = res_state_D[0].lgTgarr[lx+2];
    double x4 = res_state_D[0].lgTgarr[lx+3];
    double y1 = res_state_D[0].lgrhoarr[ly];
    double y2 = res_state_D[0].lgrhoarr[ly+1];
    double y3 = res_state_D[0].lgrhoarr[ly+2];
    double y4 = res_state_D[0].lgrhoarr[ly+3];
    
    /* ------------------------------------------- */
    /* Coefficients for 1-D interpolations         */
    /* ------------------------------------------- */
    double a[4], b[4];
    //
    a[0]= (x - x2)*(x - x3)*(x - x4)/((x1 - x2)*(x1 - x3)*(x1 - x4));
    a[1]= (x - x1)*(x - x3)*(x - x4)/((x1 - x2)*(x2 - x3)*(x4 - x2));
    a[2]= (x - x1)*(x - x2)*(x - x4)/((x1 - x3)*(x2 - x3)*(x3 - x4));       
    a[3]= (x - x1)*(x - x2)*(x - x3)/((x3 - x4)*(x4 - x1)*(x2 - x4));
    //
    b[0]= (y - y2)*(y - y3)*(y - y4)/((y1 - y2)*(y1 - y3)*(y1 - y4));
    b[1]= (y - y1)*(y - y3)*(y - y4)/((y1 - y2)*(y2 - y3)*(y4 - y2));
    b[2]= (y - y1)*(y - y2)*(y - y4)/((y1 - y3)*(y2 - y3)*(y3 - y4));       
    b[3]= (y - y1)*(y - y2)*(y - y3)/((y3 - y4)*(y4 - y1)*(y2 - y4));
    
    double fxy, fx; 
    for(int m=0; m<(int)res_state_D.size(); m++)
    {
        /* ------------------------------------- */
        /* The 2-D interpolated valued of Ai     */
        /* ------------------------------------- */
        fxy=0.0;
        for(int i=0; i<4; i++)
        {
            fx=0.0;
            for(int k=0; k<4; k++) fx+=a[k]*(res_state_D[m].AiMatrix[k+lx][i+ly]-res_state_D[m].BiVec[k+lx]); 
            fxy+=b[i]*fx;                                                          // interpolation in y 
        }
        
        Ai[m] = exp(fxy);
        
        /* ------------------------------------- */
        /* Simple 1-D interpolation to find Bi   */
        /* ------------------------------------- */     
        fx=0.0;
        for(int k=0; k<4; k++) fx+= a[k]*res_state_D[m].BiVec[k+lx];
        
        Bi[m] = exp(fx);
        
        /* ------------------------------------- */
        /* Simple 1-D interpolation for Bitot    */
        /* ------------------------------------- */     
        fx=0.0;
        for(int k=0; k<4; k++) fx+= a[k]*res_state_D[m].BitotVec[k+lx];
        
        Bitot[m] = exp(fx);
        
        /* ------------------------------------- */
        /* Getting all Ri->j interpolated values */
        /* ------------------------------------- */
        for(int i=0; i<m; i++)
        {
            fx=0.0;
            for(int k=0; k<4; k++) fx+= a[k]*res_state_D[m].RijMatrix[k+lx][i];
            
            RijVec[m][i] = exp(fx);
        }
        
        RijVec[m][m] = 0.0;
        
        for(int i=m+1; i<res_state_D[m].neqres; i++)
        {
            fx=0.0;
            for(int k=0; k<4; k++) fx+= a[k]*res_state_D[m].RijMatrix[k+lx][i];
            
            RijVec[m][i] = exp(fx);
        }
    }
    
    return;
}

//========================================================================================
void get_Bitot(double Tg, vector<double> &Bitot)
{
    if(res_state_D[0].BiVec.size()==0){ cerr << " get_Bitot:: It seems like the transition data has not been loaded. Exiting... " << endl; exit(0); }
    if(Bitot.size()!=res_state_D.size()){ cerr << " get_Bitot:: It seems like the Bitot vector was not set. Exiting... " << endl; exit(0); }    
    
    double logTg  = log(Tg);
    
    // Make sure all points are in the table including the buffer for interpolation
    if( logTg < res_state_D[0].lgTgmin || logTg > res_state_D[0].lgTgmax )
    {
        cerr << " get_Bitot:: Tg = " << Tg << " is out of table range precomputed." << endl;    
        exit(0);
    }
    
    locate_JC(&res_state_D[0].lgTgarr[0], res_state_D[0].Npoints_Tg, logTg, &res_state_D[0].row); 
    
    int lx = (int) res_state_D[0].row;    // index of x coordinate of adjacent grid point to left of P
    
    // ---------------------------------------------
    // Boundary condition                           
    // ---------------------------------------------
    if(lx+3 > res_state_D[0].Npoints_Tg || lx<0)
    { cerr << " get_Bitot:: out of bounds for Tg = " << Tg << ". Exiting! " << endl; exit(0); } 
    
    double x = logTg;
    double x1 = res_state_D[0].lgTgarr[lx];
    double x2 = res_state_D[0].lgTgarr[lx+1];
    double x3 = res_state_D[0].lgTgarr[lx+2];
    double x4 = res_state_D[0].lgTgarr[lx+3];
    
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
    for(int m=0; m<(int)res_state_D.size(); m++)
    {
        // --------------------------------------- 
        // Simple 1-D interpolation to find Bitot   
        // ---------------------------------------      
        fx=0.0;
        for(int k=0; k<4; k++) fx+= a[k]*res_state_D[m].BitotVec[k+lx];
        
        Bitot[m] = exp(fx);
    }
    
    return;
}

//========================================================================================
void get_Bitot(double Tg, double &Bitot, int ires)
{
    if(res_state_D[0].BiVec.size()==0){ cerr << " get_Bitot:: It seems like the transition data has not been loaded. Exiting... " << endl; exit(0); }
    if(ires>=(int)res_state_D.size() || ires < 0){ cerr << " get_Bitot:: It seems like you are asking for a non-existing level. Exiting... " << endl; exit(0); }    
    
    double logTg  = log(Tg);
    
    // Make sure all points are in the table including the buffer for interpolation
    if( logTg < res_state_D[0].lgTgmin || logTg > res_state_D[0].lgTgmax )
    {
        cerr << " get_Bitot:: Tg = " << Tg << " is out of table range precomputed." << endl;    
        exit(0);
    }
    
    locate_JC(&res_state_D[0].lgTgarr[0], res_state_D[0].Npoints_Tg, logTg, &res_state_D[0].row); 
    
    int lx = (int) res_state_D[0].row;    // index of x coordinate of adjacent grid point to left of P
    
    // ---------------------------------------------
    // Boundary condition                           
    // ---------------------------------------------
    if(lx+3 > res_state_D[0].Npoints_Tg || lx<0)
    { cerr << " get_Bitot:: out of bounds for Tg = " << Tg << ". Exiting! " << endl; exit(0); } 
    
    double x = logTg;
    double x1 = res_state_D[0].lgTgarr[lx];
    double x2 = res_state_D[0].lgTgarr[lx+1];
    double x3 = res_state_D[0].lgTgarr[lx+2];
    double x4 = res_state_D[0].lgTgarr[lx+3];
    
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
    for(int k=0; k<4; k++) fx+= a[k]*res_state_D[ires].BitotVec[k+lx];
    
    Bitot = exp(fx);
    
    return;
}


