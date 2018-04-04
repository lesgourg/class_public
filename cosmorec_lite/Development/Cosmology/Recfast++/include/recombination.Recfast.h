//========================================================================================
// Author: Jens Chluba
// Last modification: Oct 2010
// CITA, University of Toronto
// All rights reserved.
//========================================================================================

//========================================================================================
// In this module the solution to the RECFAST ODE-system is computed. 
// Different driver routines are provided
//========================================================================================

#ifndef RECOMBINATION_RECFAST_H
#define RECOMBINATION_RECFAST_H

//========================================================================================
// here Xe=ne/nH !!!
// params contains parameters as defined in recombination.c 
//========================================================================================

//========================================================================================
// compute ionization history like in Recfast-module
//========================================================================================
int Xe_frac(double *params, double *zarr, double *Xe_H, 
            double *Xe_He, double *Xe, double *TM, int mflag=1);

int Xe_frac(double *params, double *zarr, double *Xe_H, 
            double *Xe_He, double *Xe, double *dXe, double *TM, int mflag=1);

int Xe_frac(double *params, double *zarr, double *Xe_H, double *Xe_He, 
            double *Xe, double *dXe, double *dX_H, double *TM, int mflag=1);

//========================================================================================
// to start computation with a given initial solution at low z
// zarr contains the points at which the solution should be obtained
// rescaling of f_Xe from initial condition is used; For this dXei is needed
//========================================================================================
int Xe_frac_rescaled(double *params, const double *zarr, double *Xe_H, 
                     double *Xe_He, double *Xe, double *TM, double Xe_Hi, 
                     double Xe_Hei, double Xei, double TMi, double dXei,
                     int mflag=1);

#endif 
