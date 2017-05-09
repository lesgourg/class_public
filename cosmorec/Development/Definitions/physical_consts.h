//-------------------------------------------------------------------------------------------------------
// Author Jens Chluba May 2008
// comment: I checked NIST and updated all the constants that I found
//-------------------------------------------------------------------------------------------------------
// 23.07.2014: added possiblity to explicitly set HI A2s1s rate

#ifndef PHYSICAL_CONSTS_H
#define PHYSICAL_CONSTS_H

#include <cmath>
#include "Definitions.h"

//-------------------------------------------------------------------------------------------------------
// atomic data for HI
//-------------------------------------------------------------------------------------------------------
//const double const_HI_A2s_1s = 8.22458;         // Hydrogen 2s decay rate (old Recfast value)
//const double const_HI_A2s_1s = 8.2206;          // Hydrogen 2s decay rate (Labzowsky et al. 2005)
extern double const_HI_A2s_1s;

//-------------------------------------------------------------------------------------------------------
// atomic data for DI
//-------------------------------------------------------------------------------------------------------
const double const_DI_A2s_1s = const_HI_A2s_1s;   // Deuterium 2s decay rate

//-------------------------------------------------------------------------------------------------------
// atomic data for Helium I
//-------------------------------------------------------------------------------------------------------
const double const_HeI_A2s_1s = 51.3 ;            // HeI 2s-1s two photon rate in [s^-1] ( Drake (1969) )
//const double const_HeI_A2s_1s = 50.94;          // HeI 2s-1s two photon rate in [s^-1] ( Eric Switzer / Drake 1986 ) --> small difference with normal value
const double const_HeI_A23s_1s = 4.02e-9;         // HeI 2^3 s-1^1 s two photon rate in [s^-1] Drake (1969)
const double const_HeI_A23s_1s_M1 = 1.272426e-4;  // HeI 2^3 s-1^1 s (M1-transition) in [s^-1] Lach & Pachucki (2001)

//-------------------------------------------------------------------------------------------------------
// atomic data for Helium II
//-------------------------------------------------------------------------------------------------------
const double const_HeII_A2s_1s = 526.5;           // Helium II 2s-1s two-photon decay rate

//-------------------------------------------------------------------------------------------------------
// most fundamental constants
//-------------------------------------------------------------------------------------------------------
const double const_alpha   = 1.0/137.035999679;   //       | NIST 2008, error 6.8e-10, old: 1.0/137.03599976; 
const double const_cl      = 2.99792458e+10;      // cm/sec| NIST 2008, exact, old: same
const double const_kB_J    = 1.3806504e-23;       // J/K   | NIST 2008, error 1.7e-6, old: 1.3806503e-23;
const double const_kB      = const_kB_J*1.0e+7;   // erg/K 
const double const_e       = 1.602176487e-12;     // erg/V | NIST 2008, error 2.5e-8, old: 1.60217733e-12; 

const double const_h       = 6.62606896e-27;      // erg sec | NIST 2008, error 5.0e-8, old: 6.62606876e-27;
const double const_hbar    = 1.054571628e-27;     // erg sec | NIST 2008, error 5.0e-8, old: h/(2.0*PI);

//-------------------------------------------------------------------------------------------------------
// graviational constant
//-------------------------------------------------------------------------------------------------------
const double const_G       = 6.67428e-8;          // cm^3 g^-1 s^-2 new value from Scott/Wiki/NIST 2008
//const double const_G       = 6.67259e-8;          // cm^3 g^-1 s^-2 old value from Recfast v1.2

const double const_sigB    = 5.670400e-5;         // erg s^-1 cm^-2 K^-4 | NIST 2008, error 7.0e-6, old: same

//-------------------------------------------------------------------------------------------------------
// Bohr Radius
//-------------------------------------------------------------------------------------------------------
const double const_a0      = 5.2917720859e-9;     // cm    | NIST 2008, error 6.8e-10, old: --

//-------------------------------------------------------------------------------------------------------
//  Rydberg constant for m --> infinity (hc R ~ e^4 me/(2 hbar^2))
//-------------------------------------------------------------------------------------------------------
const double const_Ry_inf_erg=2.17987197e-11;     // erg   | NIST 2008, error 5.0e-8, old: 2.1798719e-11;
const double const_Ry_inf_icm=1.0973731568527e+5; // cm^-1 | NIST 2008, error 6.6e-12, 
                                                  // old: 1.0973731568525e+5 from Wikipedia (value 2002)
//-------------------------------------------------------------------------------------------------------
// masses 
//-------------------------------------------------------------------------------------------------------
const double const_me_gr   = 9.10938215e-28;   // gr  | NIST 2008, error 5.0e-8, old: 9.10938188e-28;
const double const_me      = 510.998910;       // keV | NIST 2008, error 2.5e-8, old: 510.998902;
const double const_mp_gr   = 1.672621637e-24;  // gr  | NIST 2008, error 5.0e-8, old: 1.67262158e-24;
const double const_mp      = 938.272013e+3;    // keV | NIST 2008, error 2.5e-8, old: 938.271998e+3; 
const double const_mn_gr   = 1.674927211e-24;  // gr  | NIST 2008, error 5.0e-8, old: -- 
const double const_mn      = 939.565346e+3;    // keV | NIST 2008, error 2.5e-8, old: 939.565330e+3;

const double const_me_mp   = 5.4461702177e-4;  // --  | NIST 2008, error 4.3e-10, old: --
const double const_me_malp = 1.37093355570e-4; // --  | NIST 2008, error 4.2e-10, old: me_mp/4.0;

const double const_amu     = 1.660538782e-24;         // gr  | NIST 2008, error 5.0e-8, old: --

const double const_mH_gr   = 1.0078250321*const_amu;  // Hydrogen mass in gr  | NIST 2010
const double const_mD_gr   = 2.0141017780*const_amu;  // Deuterium mass in gr | NIST 2010 
const double const_mDeuteron_gr = 3.34358320e-24;     // Deuteron mass in gr | NIST 2008, error 5 x 10^-8     
const double const_mHe3_gr   = 3.0160293097*const_amu;// He3 mass in gr | NIST 2010
const double const_mHe4_gr   = 4.0026032497*const_amu;// He4 mass in gr | NIST 2010

const double const_Msol    = 1.9891e+33;       // gr  | source Weigert ?

//-------------------------------------------------------------------------------------------------------
// ionisation potential of hydrogenic atoms
//
// To get the real potential one has to use E=EH/(1+me/M), 
// where M is the mass of the nucleus
//-------------------------------------------------------------------------------------------------------
const double const_EH_inf      = 13.60569193;         // eV   | NIST 2008, error 2.5e-8, old: 13.60569816;
const double const_EH_inf_ergs = const_Ry_inf_erg;    // ergs 
const double const_EH_inf_Hz   = 3.289841960361e+15;  // Hz   | NIST 2008, error 6.6e-12, old: --

//-------------------------------------------------------------------------------------------------------
// other constants
//-------------------------------------------------------------------------------------------------------
const double const_R_me    = 2.8179402894e-13;        // cm | NIST 2008, error 2.1e-9, old: 2.817940285e-13;
const double const_lambdac = 2.4263102175e-10;        // cm | NIST 2008, error 1.4e-9, old: 2.426310215e-10;

const double const_aRad    = 4.0*const_sigB/const_cl;                             // erg cm^-3 K^-4
const double const_sigT    = 8.0*PI/3.0*const_R_me*const_R_me;                    // cm^2
const double const_PIe2_mec= 3.0/8.0*const_sigT*const_cl/const_R_me;              // cm^2/sec 

//-------------------------------------------------------------------------------------------------------
// alpha-particle to proton mass ratio
//-------------------------------------------------------------------------------------------------------
const double const_malpha_mp  = 3.97259968951;           // -- | NIST 2008, error 1e-10, old: 4.0;

//-------------------------------------------------------------------------------------------------------
// Helium to hydrogen mass
//-------------------------------------------------------------------------------------------------------
//const double const_mHe4_mH  = 3.971;                   // -- | according to Wong, Moss, Scott, 2008
//const double const_mHe4_mH  = 3.9715;                  // -- | according to Grin & Hirata 2009
const double const_mHe4_mH  = 3.97152594;                // ratio computed from NIST values

const double const_mHe3_mH  = 2.99261203;                // ratio computed from NIST values

//-------------------------------------------------------------------------------------------------------
// Deuteron to proton mass ratio
//-------------------------------------------------------------------------------------------------------
const double const_mDeuteron_mp  = 1.99900750108;        // -- | NIST 2008, error 1e-10;

//-------------------------------------------------------------------------------------------------------
// Deuterium to Hydrogen mass ratio
//-------------------------------------------------------------------------------------------------------
const double const_mD_mH  = 1.99846374;                  // ratio computed from NIST values

//-------------------------------------------------------------------------------------------------------
// 28.05.2008
// 
// This is an option which we included to account for the fact that
// the helium mass is not 4*mH (Wong et al 2008) However, we only
// changed those variables (in Cosmos and the Recombination routines)
// that are important for the recombination computations.
//
//-------------------------------------------------------------------------------------------------------
const double fac_mHemH=const_mHe4_mH/4.0;        // Wong, Moss, Scott, 2008 --> helium mass is not 4*mH 
//const double fac_mHemH=1.0;                    // In this case helium mass is assumed to be 4*mH 

//-------------------------------------------------------------------------------------------------------
// Mega-parsec
//-------------------------------------------------------------------------------------------------------
//const double const_Mpc     =3.0856775807e+24;    // cm | some book
const double const_Mpc     =3.08568025e+24;      // June 2010; from web

//-------------------------------------------------------------------------------------------------------
// for conversions
//-------------------------------------------------------------------------------------------------------
const double const_h_kb    = const_h/const_kB;                           // K sec
const double const_kb_mec2 = const_kB/const_me_gr/const_cl/const_cl;     // K^-1 sec
const double const_h_mec2 = const_h/const_me_gr/const_cl/const_cl;       // sec
const double const_hcl_kb  = const_h_kb*const_cl;                        // K cm

//-------------------------------------------------------------------------------------------------------
// bremsstrahlung emission coefficient
//-------------------------------------------------------------------------------------------------------
const double const_kappa_br= 1.0/(2.0*PI)/sqrt(6.0*PI)*pow(const_lambdac, 3.0)*const_alpha; // cm^3

//-------------------------------------------------------------------------------------------------------
// DC emission coefficient
//-------------------------------------------------------------------------------------------------------
const double const_kappa_dc= 4.0*const_alpha/(3.0*PI);

#endif
