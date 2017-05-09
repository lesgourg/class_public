//========================================================================================
// Author: Jens Chluba
// Last modification: Oct 2010
// CITA, University of Toronto
// All rights reserved.
//========================================================================================

//========================================================================================
// This module contains the evaluation of the RECFAST ODE-system. The code is based on the
// C-version for RECFAST from Seager 1999
//========================================================================================

#include <cmath>
#include <iostream>
#include <cstdlib>

#include "cosmology.Recfast.h"
#include "constants.Recfast.h"
#include "evalode.Recfast.h"
#include "DM_annihilation.Recfast.h"

using namespace std;
using namespace RECFAST_physical_constants; // defined in constants.Recfast.h
using namespace RECFAST_atomic_constants;   // defined in constants.Recfast.h

//========================================================================================
// Global Variables; defined in cosmology.Recfast.h
//========================================================================================
extern struct Input input;


//========================================================================================
// Boltzmann factor
//========================================================================================
double Boltzmann(double gj, double gi, double E, double T )
{ return max(1.0e-300, ( (gj/gi) * exp(-E/(RF_kBoltz*T)) )); }

//========================================================================================
// Boltzmann factor for detailed balance relation
//========================================================================================
double SahaBoltz(double gi, double gc, double ne, double E_ion, double T)
{
    double c1;
    c1 = (RF_hPlanck/(2.0*RF_PI*RF_mElect)) * (RF_hPlanck/RF_kBoltz);
    return (ne * gi/(2.0*gc) * pow(c1, 1.5) * pow(T, -1.5) 
               * exp(min(680.0, E_ion/(RF_kBoltz*T))) );
}


//========================================================================================
// case B recombination rate for Hydrogen (in m^3 / s)
//========================================================================================
double alphaH_func(double TM)
{
  double a1 = 4.309, a2 = -0.6166, a3=0.6703, a4=0.5300, t=TM/1.0e+4;
  return a1*1.0e-19*pow(t, a2)/(1.0+a3*pow(t, a4)); 
}

//========================================================================================
// case B recombination rate for neutral Helium (in m^3 / s)
//========================================================================================
double alphaHe_func(double TM)
{
  double a1 = pow(10.0, -16.744), a2 = 0.711, T0=3.0, T1=pow(10.0, 5.114);
  double sqrt_TMT0=sqrt(TM/T0), sqrt_TMT1=sqrt(TM/T1); 
  return a1/(sqrt_TMT0*pow(1.0+sqrt_TMT0, 1.0-a2)*pow(1.0+sqrt_TMT1, 1.0+a2));
}


//========================================================================================
// evaluation of Recfast-system
//========================================================================================
void evaluate_Recfast_System(double z, double y[], double fvec[], int n)
{
    double xe, xp, xHep, nH1, nHe1;
    double TR, TM, Comp;
    double Hz;
    
    double nHTot;
    double KH, KHe, fHe;
    double alphaH, betaH;
    double alphaHe, betaHe;
    double BH, BHe, BHe2p2s;
    double SH, SHe;
    double CH, CHe;
    
    //====================================================================================
    // get Hubble & NH
    //====================================================================================
    Hz=H_z(z);
    nHTot = NH(z);
    fHe = input.fHe;
    
    //====================================================================================
    // set variables
    //====================================================================================
    xHep = y[0]; 
    xp   = y[1];
    TM   = y[2];
    xe   = xp+xHep;
    
    nH1 = (1.0-xp)*nHTot;
    nHe1 = nHTot*(fHe-xHep);
    
    //====================================================================================
    // radiation temperature at z
    //====================================================================================
    TR = TCMB(z);
    
    //====================================================================================
    // Compton-cooling term
    //====================================================================================
    Comp = 8.0*RF_sigmaT*RF_aRad*pow(TR, 4)/(3.0*Hz*(1.0+z)*RF_mElect*RF_cLight);
    
    //====================================================================================
    // Boltzmann and detailed balance factors
    //====================================================================================
    SHe = SahaBoltz(1.0, 2.0, 1.0, RF_EionHe2s, TM);
    BHe = Boltzmann(1.0, 1.0, RF_E2s1sHe, TM);
    //
    SH  = SahaBoltz(2.0, 1.0, 1.0, RF_EionH2s, TM);
    BH  = Boltzmann(1.0, 1.0, RF_E2s1sH,  TM);
    
    //====================================================================================
    // For HeI, the energy levels of 2s and 2p are different,
    // so that a boltzmann factor remains.
    //====================================================================================
    BHe2p2s = 1.0/Boltzmann(1.0, 1.0, RF_E2p2sHe, TM);
    
    //====================================================================================
    // H recombination coefficient + fugde-factor
    //====================================================================================
    alphaH = alphaH_func(TM);
    alphaH*= input.F;
    betaH  = alphaH/SH;
    
    //====================================================================================
    // He recombination coefficient 
    //====================================================================================
    alphaHe = alphaHe_func(TM);
    betaHe  = alphaHe/SHe;
    
    //====================================================================================
    // tau-Sobolev for helium and hydrogen
    //====================================================================================
    KHe = pow(RF_LyalphaHe, 3)/(Hz*8.0*RF_PI);
    KH  = pow(RF_LyalphaH , 3)/(Hz*8.0*RF_PI);
    
    //====================================================================================
    // Inhibition factors for HeI & HI
    //====================================================================================
    CHe = (1.0+KHe*RF_Lam2s1sHe*nHe1*BHe2p2s)/(1.0+KHe*(RF_Lam2s1sHe+betaHe)*nHe1*BHe2p2s);
    CH  = (1.0+KH *RF_Lam2s1sH *nH1         )/(1.0+KH *(RF_Lam2s1sH +betaH )*nH1         );
    
    //====================================================================================
    // derivatives; 1 == He, 2 == H,  3 == TM 
    //====================================================================================
    fvec[0] = (alphaHe*xe*xHep*nHTot - betaHe*(fHe-xHep)*BHe)*CHe/((1.0+z)*Hz);
    fvec[1] = (alphaH* xe*xp  *nHTot - betaH *(1.0-xp  )*BH )*CH /((1.0+z)*Hz);
    fvec[2] = Comp*xe/(1.0+xe+fHe)*(TM-TR) + 2.0*TM/(1.0+z);
    
    //====================================================================================
    // add terms due to DM annihilation
    //====================================================================================
    if(input.fDM!=0.0) 
        evaluate_DM_annihilation_terms(z, Hz, fHe, xp, xHep, CHe, CH, input.fDM, fvec);
        
    return;
}


//========================================================================================
// wrapper for ODE-solver
//========================================================================================
void fcn(int *neq, double *z, double *y, double *f)
{
    evaluate_Recfast_System(*z, y, f, *neq);
    return;
}

//========================================================================================
// wrapper for ODE-solver with rescaled derivative for Xe
// this is important to stitch the recfast solution to the output of the multi-level code
//========================================================================================
double fcn_rescaled_fac=1.0;
void set_rescaling(double ff){ fcn_rescaled_fac=ff; return; }

void fcn_rescaled(int *neq, double *z, double *y, double *f)
{
    evaluate_Recfast_System(*z, y, f, *neq);
    f[1]*=fcn_rescaled_fac;
    return;
}
