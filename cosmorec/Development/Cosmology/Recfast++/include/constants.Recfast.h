//========================================================================================================
// Author: Jens Chluba
// Last modification: Oct 2010
// CITA, University of Toronto
// All rights reserved.
//========================================================================================================
// 23.07.2014: added possiblity to explicitly set HI A2s1s rate

#ifndef RECFAST_CONSTANTS_H
#define RECFAST_CONSTANTS_H

//========================================================================================================
// Macro for value of PI
//========================================================================================================
#ifndef  RF_PI
#define  RF_PI 3.1415926535897932
#endif

namespace RECFAST_physical_constants 
{
    //====================================================================================================
    // Define physical constants. Here we use MKS definitions throughout.
    //====================================================================================================
    const double RF_cLight   =2.99792458E+08;      // Speed of light [m/s]      
    const double RF_hPlanck  =6.62606876E-34;      // Planck constant [Js]    
    const double RF_kBoltz   =1.3806503E-23;       // Boltzman constant [J/K] 
    const double RF_mElect   =9.10938188E-31;      // Electron mass [kg]        

    //====================================================================================================
    const double RF_amu     = 1.660538782e-24;     // gr  | NIST 2008, error 5.0e-8 
    const double RF_mH_gr   = 1.0078250321*RF_amu; // Hydrogen mass in gr  | NIST 2010

    //====================================================================================================
    //const double RF_mHatom   =1.67262158E-27;    // H atom mass   [kg] from c-version 
    //const double RF_mHatom   =1.673725E-27;      // H atom mass   [kg] from Recfast.for v1.2 
    //const double RF_mHatom   =1.673575E-27;      // H atom mass   [kg] from Recfast.for v1.3/1.4 
    const double RF_mHatom =RF_mH_gr*1.0e-3;       // H atom mass   [kg] latest value      

    //====================================================================================================
    const double RF_G      =6.67428e-8*1.0e-3;     // Gravitational constant [N m^2/kg^2] 
                                                   // new value from Scott/Wiki/NIST 2008

    //====================================================================================================
    const double RF_aRad   =4.0*5.670400e-8/RF_cLight; // radiation constant
    const double RF_sigmaT =6.6524616e-29;             // Thomson cross section in m^2 
    
    //====================================================================================================
    // 28.05.2008
    // This is an option which is included to account for the fact that
    // the helium mass is not 4*mH (Wong et al 2008). 
    //====================================================================================================
    const double RF_mHe4_mH  = 3.97152594;         // ratio computed from NIST values
    const double RF_fac_mHemH=RF_mHe4_mH/4.0;      // Wong, Moss, Scott, 2008 --> helium mass is not 4*mH 
    //const double RF_fac_mHemH=1.0;               // In this case helium mass is assumed to be 4*mH 

    //====================================================================================================  
    // Mega-parsec
    //====================================================================================================  
    const double RF_Mpc     =3.08568025e+24;       // June 2010; from web
}

//========================================================================================================
// Define atomic data for H and He
// This part is from the c-version of Recfast
//========================================================================================================
using namespace RECFAST_physical_constants;

namespace RECFAST_atomic_constants
{
    //====================================================================================================  
    // two-photon transition rates
    //====================================================================================================  
    // hydrogen
    //const double RF_Lam2s1sH=8.2206;  // H 2s-1s two photon rate in [s^-1] from Labzowsky 2005 
    //const double RF_Lam2s1sH=8.22458; // H 2s-1s two photon rate in [s^-1]
    extern double RF_Lam2s1sH;
    
    // helium
    //const double RF_Lam2s1sHe=52.607; // HeI 2s-1s two photon rate in [s^-1] from Labzowsky 2005 
    const double RF_Lam2s1sHe=51.3;     // HeI 2s-1s two photon rate in [s^-1] 
    
    //====================================================================================================  
    // transition wavelenght
    //====================================================================================================  
    // hydrogen
    const double RF_LyalphaH=1.215670e-07;  // H Lyman alpha wavelength in [m] 
    
    // helium
    const double RF_LyalphaHe=5.843344e-08; // HeI 2^1p - 1^1s wavelength in [m] 
    
    //====================================================================================================  
    // energy values from Recfast v1.4
    //====================================================================================================  
    // hydrogen
    const double RF_L_H_ion   =1.096787737e+7;   // level for H ion. (in m^-1)
    const double RF_L_H_alpha =8.225916453e+6;   // averaged over 2 levels
    
    // helium
    const double RF_L_He1_ion   =1.98310772e+7;  // from Drake (1993)
    const double RF_L_He2_ion   =4.389088863e+7; // from JPhysChemRefData (1987)
    const double RF_L_He_2s     =1.66277434e+7;  // from Drake (1993)
    const double RF_L_He_2p     =1.71134891e+7;  // from Drake (1993)
    
    //====================================================================================================  
    // hydrogen
    const double RF_EionH2s=(RF_L_H_ion-RF_L_H_alpha)*RF_cLight*RF_hPlanck; // H 2s ionization energy in [J] 
    const double RF_E2s1sH = RF_L_H_alpha*RF_cLight*RF_hPlanck;             // H 2s energy from 1s in [J]
    
    // helium
    const double RF_EionHe2s=6.363254e-19;  // HeI 2s ionization energy in [J]      --> nu = 9.60336246 10^14 Hz
    const double RF_EionHeII=8.7186944e-18; // HeII ionization energy in [J]        --> nu = 1.315817072 10^16 Hz
    
    const double RF_E2s1sHe=3.30301387e-18; // HeI 2s energy from 1s^2 in [J]       --> nu = 4.984877142 10^15 Hz
    const double RF_E2p2sHe=9.64908313e-20; // HeI 2p - 2s energy difference in [J] --> nu = 1.456230456 10^14 Hz 
}

#endif
