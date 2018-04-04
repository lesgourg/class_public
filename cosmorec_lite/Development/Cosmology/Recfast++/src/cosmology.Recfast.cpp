//========================================================================================
// Author: Jens Chluba
// first implementation: June 2010
// Last modification: June 2012
// CITA, University of Toronto
// All rights reserved.
//========================================================================================
// 08.06.2012: added option to use external hubble factor

#include <cmath>
#include <iostream>

#include "constants.Recfast.h"
#include "cosmology.Recfast.h"

using namespace std;
using namespace RECFAST_physical_constants; // defined in constants.Recfast.h

//========================================================================================
// Global Variables; defined in cosmology.Recfast.h
//========================================================================================
extern struct Input input;

//========================================================================================
// Hubble-function in 1/sec
//========================================================================================
double H_z_loc(double z)
{
  double Fnu, Zeq, z1=1.0+z;
  
  Fnu = input.Nnu*(7.0/8.0)*pow(4.0/11.0, 4.0/3.0); 
  Zeq = 3.0*pow(input.H0*RF_cLight, 2)
           /(8.0*RF_PI*RF_G*RF_aRad*(1.0+Fnu))
          /(pow(input.To,4))*input.OmegaM-1.0;

  return input.H0*sqrt(input.OmegaL+pow(z1, 2)
                       *(input.OmegaK+z1*input.OmegaM*(1.0+z1/(1.0+Zeq)) ) );
}

//========================================================================================
// allow setting Hubble function from outside of Recfast++ (added 08.06.2012)
//========================================================================================
double (*H_z_ptr)(double)=H_z_loc;
void set_H_pointer(double (*Hz_p)(double)){ H_z_ptr=Hz_p; return; }
void reset_H_pointer(){ H_z_ptr=H_z_loc; return; }

//========================================================================================
double H_z(double z)
{ return H_z_ptr(z); }

//========================================================================================
// hydrogen number density in m^-3
//========================================================================================
double NH(double z) 
{
  double mu_H=1.0/(1.0-input.YP);
  return 3.0*pow(input.H0, 2)*input.OmegaB/(8.0*RF_PI*RF_G*RF_mHatom*mu_H)*pow(1.0+z, 3);
}

//========================================================================================
// CMB temperature at z
//========================================================================================
double TCMB(double z){ return input.To*(1.0+z); }

//========================================================================================
// compute total contribution from relativistic particles (photon & neutrinos)
//========================================================================================
double calc_Orel(double TCMB0, double Nnu, double h100)
{ 
    double H0=h100*100.0*1.0e+5/RF_Mpc;
    double a=RF_aRad*pow(TCMB0, 4)/RF_cLight/RF_cLight; 
    double b=3.0*pow(H0, 2)/(8.0*RF_PI*RF_G);  
    return a/b*(1.0+Nnu*(7.0/8.0)*pow(4.0/11.0, 4.0/3.0));
}

//========================================================================================
//========================================================================================
