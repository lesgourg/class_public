//========================================================================================
// Author: Jens Chluba (May 2003)
// Last modification: June 2012
// CITA, University of Toronto
// All rights reserved.
//========================================================================================

//========================================================================================
// 18.06.2012
//========================================================================================
// added possibility to set Hubble factor using external function pointer
//========================================================================================

//========================================================================================
// 08.06.2012
//========================================================================================
// added possibility to load Hubble factor from some external table
//========================================================================================

//========================================================================================
// 28.05.2008
//========================================================================================
// Set the variable fac_mHemH in "physical_consts.h" to take into
// account the fact that the helium mass is not 4*mH (Wong et al 2008)
// However, we only changed those variables that are important for the
// recombination computations.
//========================================================================================

#ifndef COSMOS_H
#define COSMOS_H

#include <iostream>
#include <cmath>
#include <string>
#include <vector>

#include "physical_consts.h"

using namespace std;

//========================================================================================
// first added 08.06.2012
//========================================================================================
class Hubble
{
    
private:
    
    int mem_index_Hz;
    double zmax, zmin;
    double (* ext_Hptr)(const double *);

    bool loaded;
    bool pointer;
    
public:

    Hubble();
    Hubble(const vector<double> &z, const vector<double> &Hz);
    Hubble(double (* Hptr)(const double *));
    ~Hubble();
    
    void init(const vector<double> &z, const vector<double> &Hz);
    void init(const double *z, const double *Hz, const int nz);
    void init(double (*Hptr)(const double *));
    void clear();
    
    bool check_limits(double z);
    
    double H(double z);
};


//========================================================================================
class Cosmos
{

private:
// --------------- Cosmology constants ------------

    double Nnu;                                 // number of neutrino species: 3 
    double g_rel;                               // 1.0+3.0*(7.0/8.0)*pow(4.0/11.0, 4.0/3.0)
    double zsRe;                                // start recombination

    //====================================================================================
    // initial data to be set 
    //====================================================================================
    double h100;                                // h=H0/100 km s^-1 Mpc^-1
    double T_CMB0;                              // Temperature of the CMB today in [K]
    double Y_p;                                 // Helium energy density 
    double O_T;                                 // total energy density 
    double O_mh2;                               // matter energy density 
    double O_bh2;                               // baryon energy density 
    double O_L;                                 // cosmological constant
    double fudge;                               // Recfast fudge factor
    double fDM_annihilation;
    int CR_correction; 
    
    //====================================================================================
    // dependent consts
    //====================================================================================
    double h2, H0;                              // Hubble constant today 
    double rho_c, rho_g, rho_g_gr;              // gr cm^-3  
    double O_Th2, O_kh2, O_Lh2, O_m, O_b, O_k;
    double T27, Theta_CMB0, Theta27;
    double O_g, O_gh2;                          // photon density
    double O_rel, O_relh2;                      // relativistic particle density 
    double z_eq;                                // Matter radiation equality (Wu 1995)
    double Nb0;                                 // Baryon number density today

    //====================================================================================
    // auxilliary variables
    //====================================================================================
    void set_initials(double Nnuval);
    double power;
    
    //====================================================================================
    // auxilliary functions
    //====================================================================================
    void calc_dependent_consts();

    //====================================================================================
    // setup for Xe & X1s etc
    //====================================================================================
    int recombination_history(int nzpts, double zs, double ze, double *zarr, double *Xe_H, 
                              double *Xe_He, double *Xe, double *TM);
    
    int recombination_history(int nzpts, double zs, double ze, double *zarr, 
                              double *Xe_H, double *Xe_He, double *Xe, 
                              double *dXe, double *dX_H, double *TM);

    //====================================================================================
    // spline setup
    //====================================================================================
    int calc_coeff_X_spline(int M, double zstart, double zend);
    
    int n_Xe;
    int memindex_Xe; 
    int memindex_dXe; 
    int memindex_XeH; 
    int memindex_dXeH; 
    int memindex_rho; 
    int memindex_XeHe; 
    bool splines_are_charged;
    bool spline_memory_allocated;

    void check_splines();
    double calc_spline(double z, int memindex);
    double calc_spline_log(double z, int memindex);
    //====================================================================================

    double RF_hPlanck, RF_kBoltz, RF_mElect; // local physical constants     
    
    //====================================================================================
    // to load an external hubble function from a table (08.06.2012)
    //====================================================================================
    Hubble loaded_Hz;
    
public:

    double zs, ze;

    //====================================================================================
    // Konstructors and Destructors
    //====================================================================================
    Cosmos();
    Cosmos(const double h100, const double T0, const double Yp, 
           const double densities[7], const double zstart, 
           const double zend, const int n_Xe_pts, 
           const double Nnu=3.04);
    
    ~Cosmos();

    void init(const double h100, const double T0, const double Yp, 
              const double densities[7], const double zstart, 
              const double zend, const int n_Xe_pts, 
              const double Nnu=3.04);
    
    void init_splines();
    
    //====================================================================================
    // to load an external hubble function from a table (08.06.2012)
    //====================================================================================
    void init_Hubble(const vector<double> &z, const vector<double> &Hz);
    void init_Hubble(const double *z, const double *Hz, const int nz);
    void init_Hubble(double (* Hptr)(const double *));
    void clear_Hubble(){ loaded_Hz.clear(); }
    
    //====================================================================================
    void display_cosmology();
    void dump_cosmology(string filename);

    double Get_Y_p() const { return Y_p;} 
    double Get_Nnu() const { return Nnu;}
    //
    double Get_O_bh2() const { return O_bh2;} 
    double Get_O_b() const   { return O_b;  } 
    //
    double Get_O_DMh2() const { return O_mh2-O_bh2;} 
    double Get_O_DM() const   { return O_m-O_b;    } 
    //
    double Get_O_mh2() const { return O_mh2;} 
    double Get_O_m() const   { return O_m;  } 
    //
    double Get_O_Lh2() const { return O_Lh2;  } 
    double Get_O_L() const   { return O_L;    } 
    // 
    double Get_O_g() const  { return O_g;   } 
    double Get_rho_c() const { return rho_c; } 
    double Get_rho_g() const { return rho_g; } 
    double Get_T27()   const { return T27;   } 
    double Get_T_CMB0() const{ return T_CMB0;} 
    double Get_power() const { return power; } 

    //====================================================================================
    // CMB temperature
    //====================================================================================
    double TCMB(double z)    {return T_CMB0*(1.0+z);} 
    double Ttoz(double T_g)  {return T_g/T_CMB0-1.0;} 

    //====================================================================================
    // CMB temperature Theta = k Tg / me c^2
    //====================================================================================
    double ThetaCMB(double z){return Theta_CMB0*(1.0+z);} 
    double Thetatoz(double Theta_g){return Theta_g/Theta_CMB0-1.0;} 

    //====================================================================================
    // conversion Theta <--> Tg
    //====================================================================================
    double TtoTheta(double T_g)    {return const_kb_mec2*T_g;} 
    double ThetatoT(double Theta_g){return Theta_g/const_kb_mec2;} 

    double calc_Orel(double TCMB0, double Nnu, double h100)
    { 
        double H0=100.0*h100*(1.0e+3*1.0e+2)/const_Mpc;
        double a=PI*PI/15.0*const_kB*TCMB0*pow(const_kB*TCMB0/(const_hbar*const_cl), 3)
                      /pow(const_cl, 2); 
        double b=3.0*pow(H0, 2)/(8.0*PI*const_G);  
        return a/b*(1.0+Nnu*(7.0/8.0)*pow(4.0/11.0, 4.0/3.0));
    }
    
    //====================================================================================
    // compute Recfast system with a given initial solution at low z
    //====================================================================================
    int recombine_using_Recfast_system(int nzpts, double zi, double ze, double fDM, 
                                       double Xe_Hi, double Xe_Hei, double Xei, double TMi, 
                                       double dXei, const double *zarr, double *Xe_H, 
                                       double *Xe_He, double *Xe, double *TM);

    //====================================================================================
    // Xe & X1s
    //====================================================================================
    double Xe(double z)                // this is Ne/Nb/[Ne(inf)/Nb(inf)] ~ Ne/Nb/(1-Yp/2)
    { 
        if(z>=zsRe) return 1.0; 
        // --- Xe (1-Y_p/2) =ne/nb ! but in recombination Xe=ne/nh 
        // --> factor of (1.0-Y_p)/(1.0-Y_p/2.0) 
        return (1.0-Y_p)/(1.0-Y_p/2.0*(2.0-1.0/fac_mHemH))*calc_spline_log(z, memindex_Xe);
    } 
    
    double Xe_Seager(double z)
    {                                              // xe as given by Seager 2000. Xe=ne/nh
        if(z>=zsRe) return (1.0-Y_p/2.0*(2.0-1.0/fac_mHemH))/(1.0-Y_p); 
        return calc_spline_log(z, memindex_Xe); 
    } 
    
    double dXe_dz_Seager(double z)
    {                                              // xe as given by Seager 2000. Xe=ne/nh
        if(z>=zsRe) return 0.0; 
        return calc_spline(z, memindex_dXe); 
    } 
    
    double dXe_dt_Seager(double z)
    {                                              // xe as given by Seager 2000. Xe=ne/nh
        if(z>=zsRe) return 0.0; 
        return -H(z)*(1.0+z)*calc_spline(z, memindex_dXe); 
    }
    
    //====================================================================================
    // X1s = NH1s/NH
    //====================================================================================
    double dX1s_dz(double z)
    { 
        if(z>=zsRe) return 0.0; 
        return -calc_spline(z, memindex_dXeH);
    } 
    
    double dX1s_dt(double z)
    { 
        if(z>=zsRe) return 0.0; 
        return +H(z)*(1.0+z)*calc_spline(z, memindex_dXeH);
    } 

    //====================================================================================
    // 30.06.2009
    //====================================================================================
    double SahaBoltz_HeII(double T)
    {
        // HeII ionization energy in [J] --> nu = 1.315817072 10^16 Hz
        //double RF_EionHeII =8.7186944e-18;        
        double RF_EionHeII =4.0*const_EH_inf_ergs*1.0e-7/(1.0+const_me_malp);  

        double c1 = 2.0*PI*RF_mElect*RF_kBoltz*T/(RF_hPlanck*RF_hPlanck);
        double NHm3=NH(Ttoz(T))*1.0e+6;
        double g = 4.0*pow(c1, 1.5)*exp(-RF_EionHeII/(RF_kBoltz*T))/NHm3, d=1.0-g;
        double Xe= 0.5*(d + sqrt(d*d + 4.0*g*(1.0+fHe())));
        
        return Xe;
    }
    
    double SahaBoltz_HeI1s(double T)
    {
        double RF_EionHeI  =3.93933e-18;         // HeI ionization energy in [J]       
        
        double c1 = 2.0*PI*RF_mElect*RF_kBoltz*T/(RF_hPlanck*RF_hPlanck);
        double NHm3=NH(Ttoz(T))*1.0e+6;
        double g = 4.0*pow(c1, 1.5)*exp(-RF_EionHeI/(RF_kBoltz*T))/NHm3;
        double Xe=Xe_Seager(Ttoz(T));
        double XHeI1s=Xe*(1.0+2.0*fHe()-Xe)/(2.0*Xe+g);
        
        return XHeI1s;
    }
    
    double XHe1s(double z)
    { 
        if(z>=zsRe) return 0.0; 
        if(z>=3500.0) return SahaBoltz_HeI1s(TCMB(z));
        return fHe()-calc_spline_log(z, memindex_XeHe); 
    }
    
    double SahaBoltz_HI1s(double T)
    {
        // HeI ionization energy in [J] 
        //double RF_EionHI   =2.17869e-18;               
        double RF_EionHI   =const_EH_inf_ergs*1.0e-7/(1.0+const_me_mp);
        
        double c1 = 2.0*PI*RF_mElect*RF_kBoltz*T/(RF_hPlanck*RF_hPlanck);
        double NHm3=NH(Ttoz(T))*1.0e+6;
        double g = pow(c1, 1.5)*exp(-RF_EionHI/(RF_kBoltz*T))/NHm3;
        double Xe=Xe_Seager(Ttoz(T));
        double X1s=Xe/(Xe+g);
        
        return X1s;
    }

    double X1s(double z)
    { 
        if(z>=zsRe) return 0.0; 
        if(z>=3500.0) return SahaBoltz_HI1s(TCMB(z));
        return 1.0-calc_spline_log(z, memindex_XeH);
    }

    //====================================================================================
    // (approximate) free proton fraction from Recfast
    //====================================================================================
    double Xp(double z){ return 1.0-X1s(z);}          
    
    //====================================================================================
    // Te_Tg from Recfast
    //====================================================================================
    double Te_Tg(double z){ if(z>=zsRe) return 1.0; return calc_spline(z, memindex_rho); } 
    double Te(double z){return TCMB(z)*Te_Tg(z);} 

    //====================================================================================
    // different number densities
    //====================================================================================
    double Nb(double z){return Nb0*pow(1.0+z,3);}                                   // Nb==rho_b/mH
    double Ne(double z){return Xe(z)*(1.0-Y_p/2.0*(2.0-1.0/fac_mHemH))*Nb(z);}
    double N1s(double z){return X1s(z)*NH(z);}
    double NH(double z){return (1.0-Y_p)*Nb(z); }                                   // number density of hydrogen nuclei 1/cm^3
    double Np(double z){return Xp(z)*NH(z); }                                       // number density of free protons    1/cm^3
    //
    double NHe(double z){return Y_p/(4.0*fac_mHemH)*Nb(z);}                         // number density of helium nuclei   1/cm^3
    double NHeI(double z);                                                          // number density of HeI    1/cm^3
    double NHeII(double z);                                                         // number density of HeII   1/cm^3
    double NHeIII(double z);                                                        // number density of HeIII  1/cm^3
    //
    double fHe(double Ypv){ return Ypv/(4.0*fac_mHemH)/(1.0-Ypv); }                 // nHe/nH
    double fHe()          { return fHe(Y_p); }                                      // nHe/nH
    //
    double rho_g_g_cm3(double zz){return rho_g*pow(1.0+zz, 4);}                     // photon energy density g/cm^3
    double rho_g_1_cm3(double zz){return rho_g_gr*pow(1.0+zz, 4);}                  // photon energy density 1/cm^3 (mec2)
    double n_gamma(double zz){return 2.0*1.20206*1.7595e+30*pow(ThetaCMB(zz), 3);}  // photon number density in 1/cm^3

    //====================================================================================
    // Hubble parameter 
    //====================================================================================
    double H(double z);
    
    //====================================================================================
    // different derivatives 
    //====================================================================================
    double dt_dz(double z){ return -1.0/H(z)/(1.0+z);}                        
    double dz_dt(double z){ return -H(z)*(1.0+z);}                        

    //====================================================================================
    // expansion time
    //====================================================================================
    double t_exp(double z){ return  1.0/H(z); } 
    
    //====================================================================================
    // cosmological time
    //====================================================================================
    double t(double z);                 

    //====================================================================================
    // times for energy exchange
    //====================================================================================
    double t_eg(double z){ return 1.0/(const_sigT*const_cl*Ne(z));}                   // 1/(sigT ne c) 
    double t_C(double z){ return t_eg(z);}                                            // 1/(sigT ne c) 
    double t_K(double z){return t_eg(z)/(4.0*ThetaCMB(z));}                           // 1/(4 Theta sigT ne c)
    double dtau_C_dz(double z){ return dt_dz(z)/t_eg(z);}                             // dtau_C/dz 
    double t_ep(double z);                                                            // Spitzer formula

    //====================================================================================
    // Integral^zh_z (1+z)^pow |dt_dz| dz 
    //====================================================================================
    double E_Int(double zh, double z, double pow);                                 
};

#endif
