//========================================================================================
// Author: Jens Chluba
// Last modification: Oct 2010
// CITA, University of Toronto
// All rights reserved.
//========================================================================================

//========================================================================================
// In this module the solution to the RECFAST ODE-system is computed. 
// Different driver routines are provided. 
//========================================================================================

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "constants.Recfast.h"
#include "ODE_solver.Recfast.h"
#include "cosmology.Recfast.h"
#include "evalode.Recfast.h"
#include "recombination.Recfast.h"
#include "Rec_corrs_CT.Recfast.h"

using namespace std;
using namespace RECFAST_physical_constants; // defined in constants.Recfast.h
using namespace RECFAST_atomic_constants;   // defined in constants.Recfast.h

//========================================================================================
// Global Variables; defined in cosmology.Recfast.h
//========================================================================================
struct Input input;


//========================================================================================
// Saha formula for HeIII
//========================================================================================
double SahaBoltz_HeIII(double nH, double fHe, double T)
{
    double c1 = 2.0*RF_PI*RF_mElect*RF_kBoltz*T/(RF_hPlanck*RF_hPlanck);
    double g = pow(c1, 1.5)*exp(-RF_EionHeII/(RF_kBoltz*T))/nH, d=1.0+fHe-g;
    double Xe= 0.5*(d + sqrt(d*d + 4.0*g*(1.0+2.0*fHe)));
    
    return Xe;
}

//========================================================================================
// Saha formula for HeII
//========================================================================================
double SahaBoltz_HeII(double nH, double fHe, double T)
{
    double c1 = 2.0*RF_PI*RF_mElect*RF_kBoltz*T/(RF_hPlanck*RF_hPlanck);
    double g = 4.0*pow(c1, 1.5)*exp(-RF_EionHeII/(RF_kBoltz*T))/nH, d=1.0-g;
    double Xe= 0.5*(d + sqrt(d*d + 4.0*g*(1.0+fHe)));
    
    return Xe;
}


//========================================================================================
// functions to communicate with cosmology object
//========================================================================================
int Xe_frac(double params[14], double *zarr, double *Xe_H, 
            double *Xe_He, double *Xe, double *TM, int mflag)
{
  int flg;
  long int nzpts=(int)params[0];
  double *dXe=new double[nzpts];
  double *dX_H=new double[nzpts];
  
  flg=Xe_frac(params, zarr, Xe_H, Xe_He, Xe, dXe, dX_H, TM, mflag);
  
  delete[] dXe;
  delete[] dX_H;

  return flg;
}  

//========================================================================================
// functions to communicate with cosmology object
//========================================================================================
int Xe_frac(double params[14], double *zarr, double *Xe_H, double *Xe_He, 
            double *Xe, double *dXe, double *TM, int mflag)
{
  int flg;
  long int nzpts=(int)params[0];
  double *dX_H=new double[nzpts];
  
  flg=Xe_frac(params, zarr, Xe_H, Xe_He, Xe, dXe, dX_H, TM, mflag);
  
  delete[] dX_H;

  return flg;
}  


//========================================================================================
// main functions
//========================================================================================
void Set_input_variables(double params[14])
{
    double fHe= params[3]/(4.0*RF_fac_mHemH*(1.0-params[3]));

    input.zstart=params[1];
    input.zend=params[2];
    input.YP=params[3];
    input.To=params[4];
    input.OmegaM=params[5];
    input.OmegaB=params[6];
    input.OmegaL=params[7];
    input.OmegaK=params[8];
    input.h100 = params[9];
    input.Nnu = params[10];
    if(params[11]==0.0) input.F=1.14;
    else input.F=params[11];
    
    input.fDM=params[12];
    input.switch_on_recombination_corrs=(int)params[13];
    
    input.H0 = params[9]*100.0*1.0e+5/RF_Mpc;
    input.fHe=fHe;  
    
    return;
}

//========================================================================================
// compute ionization history like in Recfast-module
//========================================================================================
int Xe_frac(double params[14], double *zarr, double *Xe_H, double *Xe_He, double *Xe, 
            double *dXe, double *dX_H, double *TM, int mflag)
{
    //====================================================================================
    // params has to contain 11 parameters in the following sequence: 
    // nzpts, zs, ze, Y_p, T0, OmT, OmB, OmL, OmK, h100, fudge
    //------------------------------------------------------------------------------------
    // nzpts: number of redshift points;
    // zs, ze: starting and ending redshift; 
    // if the fudge factor ==0.0 --> will use standard value 1.14
    //====================================================================================
    long int i, j=0, nzpts=(int)params[0];
    double zc_HeIII=9000.0, zc_HeIIIS=4800.0, zc_HeII=3600.0;
    
    //====================================================================================
    // Set in starting values
    //====================================================================================
    Set_input_variables(params);
    double fHe=input.fHe;
    
    //====================================================================================
    // above zc_HeIII everything is completely ionized
    //====================================================================================
    if (input.zstart > zc_HeIII && input.zstart > input.zend) 
    {
        for (j=0, i=(int)input.zstart; i>=max(zc_HeIII, input.zend)+10; i-=10, j++) 
        {
            zarr[j] = i;
            Xe_H[j] = 1.0;
            Xe_He[j]= fHe;
            TM[j]   = input.To*(1.0+zarr[j]);
            Xe[j]   = 1.0+2.0*fHe;
            dXe[j]  = 0.0;
            dX_H[j] = 0.0;
        }
        input.zstart = max(zc_HeIII, input.zend);
    }
    
    //====================================================================================
    // now HeIII recombination starts and HeII <--> HeIII are in Saha equilibrium
    //====================================================================================
    if (input.zstart > zc_HeIIIS && input.zstart > input.zend) 
    {
        for (i=(int)input.zstart; i>=max(zc_HeIIIS, input.zend)+1; i--, j++) 
        {
            zarr[j] = i;
            Xe_H[j] = 1.0;
            Xe_He[j]= fHe;
            TM[j]   = input.To*(1.0+zarr[j]);
            Xe[j]   = SahaBoltz_HeIII(NH(zarr[j]), fHe, TM[j]);
            dXe[j]   = 0.0;
            dX_H[j]   = 0.0;
        }
        input.zstart = max(zc_HeIIIS, input.zend);
    }
    
    //====================================================================================
    // Only need to start integrating ODE before HeI recombination. Anything before that 
    // time is just a constant ionization fraction. xe = 1.0 + fHe 
    //====================================================================================  
    if (input.zstart > zc_HeII && input.zstart > input.zend) 
    {
        for (i=(int)input.zstart; i>=max(zc_HeII, input.zend)+10; i-=10, j++) 
        {
            zarr[j] = i;
            Xe_H[j] = 1.0;
            Xe_He[j]= fHe;
            TM[j]   = input.To*(1.0+zarr[j]);
            Xe[j]   = 1.0+fHe;
            dXe[j]   = 0.0;
            dX_H[j]   = 0.0;
        }
        input.zstart = max(zc_HeII, input.zend);
    }
    
    zarr[j] = input.zstart;
    Xe_H[j] = 1.0;
    Xe_He[j]= fHe;
    TM[j]   = input.To*(1.0+zarr[j]);
    Xe[j]   = 1.0+fHe;
    dXe[j]   = 0.0;
    dX_H[j]   = 0.0;
    
    //====================================================================================
    // reading recombination correction as given by Chluba & Thomas 2010
    //====================================================================================
    if(input.switch_on_recombination_corrs==1) read_recombination_correction(mflag);
    
    //====================================================================================
    // Now He I recombination starts
    //====================================================================================
    if (input.zstart > input.zend) 
    {
        double zs, ze;
        
        //================================================================================
        // initialise zarr (linear in z)
        //================================================================================
        double Dz=(input.zstart-input.zend)/(double)(nzpts-j-1);
        for(zarr[j]=input.zstart, i=j+1; i<nzpts; i++) zarr[i]=zarr[i-1]-Dz;
        zarr[nzpts-1]=0.0;
        
        //================================================================================
        // Set starting solution. Everything is ionized, 
        // and the matter temperature = the radiation temperature           
        // 0 == xHe = nHe/nHTot
        // 1 == xp = np/nHTot
        // 2 == TM (matter temperature)
        //================================================================================
        
        //================================================================================
        // set up Solution memory for ODE-solver
        //================================================================================
        ODE_solver_Solution Sz;
        zs=zarr[j];
        //
        Sz.z=zs;
        Sz.y.resize(3);
        Sz.y[0]=Xe_He[j];
        Sz.y[1]=Xe_H[j];
        Sz.y[2]=TM[j];
        //
        Sz.dy.resize(3);
        ODE_Solver_set_up_solution_and_memory(zs, Sz, fcn);     
        
        //================================================================================
        // solve for Xe(z)
        //================================================================================
        for(i=j+1; i<nzpts; i++)
        {
            ze=zarr[i];
            
            ODE_Solver_Solve_history(zs, ze, Sz);
            
            //============================================================================
            // save step
            //============================================================================
            zs=ze;
            Xe_He[i]= Sz.y[0];
            Xe_H[i] = Sz.y[1];
            TM[i]   = Sz.y[2];
            Xe[i]   = Xe_H[i]+Xe_He[i];
            if(input.switch_on_recombination_corrs==1) Xe[i]*=Recombination_correction_factor(ze);
            //
            dXe[i]  = Sz.dy[0]+Sz.dy[1];
            dX_H[i] = Sz.dy[1];
        }
    }
    
    return 0;
}


//========================================================================================
// to start computation with a given initial solution at low z
// zarr contains the points at which the solution should be obtained
// rescaling of f_Xe from initial condition is used; For this dXei is needed
//========================================================================================
int Xe_frac_rescaled(double params[14], const double *zarr, double *Xe_H, 
                     double *Xe_He, double *Xe, double *TM, 
                     double Xe_Hi, double Xe_Hei, double Xei, double TMi, double dXei,
                     int mflag)
{
    //====================================================================================
    // params has to contain 11 parameters in the following sequence: 
    // nzpts, zs, ze, Y_p, T0, OmT, OmB, OmL, OmK, h100, fudge
    //------------------------------------------------------------------------------------
    // nzpts: number of redshift points;
    // zs, ze: starting and ending redshift; 
    // if the fudge factor ==0.0 --> will use standard value 1.14
    //====================================================================================
    long int i, j=0, nzpts=(int)params[0];
    
    //====================================================================================
    // Set in starting values
    //====================================================================================
    Set_input_variables(params);

    Xe_H[j] = Xe_Hi;
    Xe_He[j]= Xe_Hei;
    TM[j]   = TMi;
    Xe[j]   = Xei;
    
    //====================================================================================
    // reading recombination correction as given by Chluba & Thomas 2010
    //====================================================================================
    if(input.switch_on_recombination_corrs==1) read_recombination_correction(mflag);
        
    //====================================================================================
    // everything starts with HeI recombination
    //====================================================================================
    if (input.zstart > input.zend) 
    {
        double zs, ze;
        
        //================================================================================
        // set up Solution memory for ODE-solver
        //================================================================================
        ODE_solver_Solution Sz;
        zs=zarr[j];
        //
        Sz.z=zs;
        Sz.y.resize(3);
        Sz.y[0]=Xe_He[j];
        Sz.y[1]=Xe_H[j];
        Sz.y[2]=TM[j];
        //
        Sz.dy.resize(3);
        ODE_Solver_set_up_solution_and_memory(zs, Sz, fcn_rescaled);        
        
        //================================================================================
        // determine rescaling factor
        //================================================================================
        double fcn_Xe[4];
        int neq=3;
        fcn(&neq, &zs, &Sz.y[0], fcn_Xe);
        double ff=dXei/fcn_Xe[1];
        set_rescaling(ff);
        
        if(mflag>=1) cout << " Xe_frac:: rescale factor= " << ff << endl;
        
        //================================================================================
        // solve for Xe(z)
        //================================================================================
        for(i=j+1; i<nzpts; i++)
        {
            ze=zarr[i];
            
            ODE_Solver_Solve_history(zs, ze, Sz);
            
            //============================================================================
            // save step
            //============================================================================
            zs=ze;
            Xe_He[i]= Sz.y[0];
            Xe_H[i] = Sz.y[1];
            TM[i]   = Sz.y[2];
            Xe[i]   = Xe_H[i]+Xe_He[i];
            if(input.switch_on_recombination_corrs==1) Xe[i]*=Recombination_correction_factor(ze);
        }
    }
    
    return 0;
}
