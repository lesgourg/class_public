//========================================================================================
// Author: Jens Chluba (May 2003)
// Last modification: June 2012
// CITA, University of Toronto
// All rights reserved.
//========================================================================================

//========================================================================================
// 18.06.2012
//========================================================================================
// added possibility to set Hubble factor using external function
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

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cmath>

#include "Cosmos.h"
#include "cosmology.Recfast.h"
#include "recombination.Recfast.h"
#include "physical_consts.h"
#include "routines.h"
#include "Integration_routines.h"

using namespace std;

//========================================================================================
// pointer for GSL integrals (avoid type-casting)
//========================================================================================
Cosmos *GSL_cosmosptr=NULL;

//========================================================================================
// level of communication
//========================================================================================
int Cosmos_mess_flag=0;

//========================================================================================
//
// Class Hubble
//
//========================================================================================

//========================================================================================
// Konstructors and Destructors
//========================================================================================
Hubble :: Hubble(){ loaded=pointer=0; mem_index_Hz=-10; ext_Hptr=NULL; }

Hubble :: Hubble(const vector<double> &z, const vector<double> &Hz)
{ 
    loaded=pointer=0; mem_index_Hz=-10; ext_Hptr=NULL;
    init(z, Hz); 
}

Hubble :: Hubble(double (*Hptr)(const double *))
{
    loaded=pointer=1; mem_index_Hz=-10; ext_Hptr=Hptr;
}

Hubble :: ~Hubble(){ if(loaded) clear(); }

//========================================================================================
void Hubble :: init(const double *z, const double *Hz, const int nz)
{
    vector<double> lz(nz), lHz(nz);  
    
    // check ordering
    if(z[0]<z[1]) 
    {
        if(z[0]==0) lz[0]=log(1.0e-10);
        else lz[0]=log(z[0]);
        
        zmin=exp(lz[0]); 
        zmax=z[nz-1];
        
        for(int k=1; k<nz; k++) 
        {
            lz[k]=log(z[k]);
            lHz[k]=log(Hz[k]);            
        }
    }
    else 
    {
        if(z[nz-1]==0) lz[0]=log(1.0e-10);
        else lz[0]=log(z[nz-1]);
        
        zmin=exp(lz[0]); 
        zmax=z[0];
        
        for(int k=1; k<nz; k++) 
        {
            lz[k]=log(z[nz-1-k]);
            lHz[k]=log(Hz[nz-1-k]);            
        }
    }
    
    if(mem_index_Hz==-10) 
        mem_index_Hz=calc_spline_coeffies_JC(nz, &lz[0], &lHz[0], "Hubble:: Hz");
    
    else update_spline_coeffies_JC(mem_index_Hz, nz, &lz[0], &lHz[0], "Hubble:: Hz");
    
    loaded=1;
    lz.clear();
    lHz.clear();
    
    return;
}

void Hubble :: init(const vector<double> &z, const vector<double> &Hz)
{
    int nz=z.size();
    
    init(&z[0], &Hz[0], nz);
    
    return;
}

void Hubble :: init(double (*Hptr)(const double *))
{
    loaded=pointer=1; ext_Hptr=Hptr;
    return;
}

//========================================================================================
void Hubble :: clear()
{
    if(mem_index_Hz!=-10) free_spline_JC(mem_index_Hz);
    loaded=pointer=0; 
    mem_index_Hz=-10;
    ext_Hptr=NULL;
    
    return;
}

bool Hubble :: check_limits(double z)
{
    if(pointer) return 1;
    else if(loaded==1 && z>zmin && z<zmax) return 1;
    
    return 0;
}

//========================================================================================
double Hubble :: H(double z)
{
    if(pointer) return ext_Hptr(&z); 
    else if(loaded) return exp(calc_spline_JC(log(z), mem_index_Hz, " Hubble:: H(z)"));
        
    cerr << " Hubble:: Hubble was not initialized " << endl;
    
    return 0;
}
//========================================================================================



//========================================================================================
//
// Class Cosmos
//
//========================================================================================

//========================================================================================
// Konstructors and Destructors
//========================================================================================
void Cosmos :: init(const double h, const double T0, const double Yp, 
                    const double densities[7], 
                    const double zstart, const double zend, 
                    const int n_Xe_pts, const double Nnuval)
{
    //====================================================================================
    // The array densities[] has to contain the densities in the sequence:     
    // O_T, O_mh2, O_bh2, O_L                                                  
    // The array flag[] has to contain different flags in the sequence:        
    //====================================================================================
    h100=h; T_CMB0=T0; Y_p=Yp; Nnu=Nnuval;
    O_T=densities[0]; O_mh2=densities[1]; O_bh2=densities[2]; O_L=densities[3];
    O_T+=calc_Orel(T0, Nnu, h100);
    set_initials(Nnu);
    
    fudge=densities[4];
    if(fudge==0.0) fudge=1.14;
    
    fDM_annihilation=densities[5];
    CR_correction=(int)densities[6]; 

    //====================================================================================
    // for Xe-spline interpolations...
    //====================================================================================
    zs=zstart; ze=zend; 
    n_Xe=n_Xe_pts;

    calc_dependent_consts();
 
    splines_are_charged=0;
    
    //RF_hPlanck  =6.62606876E-34;      // Planck's constant [Js]    
    //RF_kBoltz   =1.3806503E-23;       // Boltzman's constant [J/K] 
    //RF_mElect   =9.10938188E-31;      // Electron mass [kg]  
    RF_hPlanck  =const_h*1.0e-7;      // Planck's constant [Js]    
    RF_kBoltz   =const_kB_J;          // Boltzman's constant [J/K] 
    RF_mElect   =const_me_gr*1.0e-3;  // Electron mass [kg]        
    
    return;
}

void Cosmos :: init_Hubble(const vector<double> &z, const vector<double> &Hz)
{
    loaded_Hz.init(z, Hz);
    
    return;
}

void Cosmos :: init_Hubble(const double *z, const double *Hz, const int nz)
{
    loaded_Hz.init(z, Hz, nz);
    
    return;
}

void Cosmos :: init_Hubble(double (* Hptr)(const double *))
{
    loaded_Hz.init(Hptr);
    
    return;
}

//========================================================================================
Cosmos :: Cosmos(){ spline_memory_allocated=0; }
Cosmos :: Cosmos(const double h, const double T0, const double Yp, 
                 const double densities[7], 
                 const double zstart, const double zend, 
                 const int n_Xe_pts, const double Nnu)
{ 
    spline_memory_allocated=0; 
    init(h, T0, Yp, densities, zstart, zend, n_Xe_pts, Nnu); 
}

Cosmos :: ~Cosmos(){}

//========================================================================================
// private functions
//========================================================================================
void Cosmos :: set_initials(double Nnuval)
{
    g_rel=1.0+Nnuval*(7.0/8.0)*pow(4.0/11.0, 4.0/3.0);
    zsRe=1.0e+4;

    return;
}

void Cosmos :: init_splines()
{
    calc_coeff_X_spline(n_Xe, zsRe, 0.0);
    splines_are_charged=1;

    return;
}

void Cosmos :: calc_dependent_consts()
{
    h2=h100*h100; H0=100.0*h100*(1.0e+3*1.0e+2)/const_Mpc;  // H = h *100 * km/s/Mpc 
                                                            //   = h *100 * 10^5  cm/s/Mpc  
    T27=T_CMB0/2.7; Theta_CMB0=TtoTheta(T_CMB0); Theta27=T27;

    rho_c=3.0*pow(H0, 2)/(8.0*PI*const_G);                  //  1.87882e-29 h2    g cm^-3  

    rho_g=PI*PI/15.0*const_kB*T_CMB0                        //  4.47378e-34 T27^4 g cm^-3
          *pow(const_kB*T_CMB0/(const_hbar*const_cl), 3)/pow(const_cl, 2);

    rho_g_gr=rho_g/const_me_gr;
    
    O_g=rho_g/rho_c; O_gh2=O_g*h2;
    O_rel=O_g*g_rel; O_relh2=O_rel*h2;
    z_eq=O_mh2/O_relh2 -1.0;                     //  4.199628e+4/g_rel*pow(T27,-4)*O_mh2;

    O_m=O_mh2/h2; O_b=O_bh2/h2;
    O_Th2=O_T*h2; O_Lh2=O_L*h2;
    O_k=1.0-O_T; O_kh2=O_k*h2;

    //====================================================================================
    // 28.05.2008
    // we changed the definition of Nb to Nb=rho_b/mH, where mH is the neutral HI mass
    //====================================================================================
    //Nb0=O_b*rho_c/const_mp_gr/XIFAC;           //  1.12328e-5 cm^-3
    Nb0=O_b*rho_c/const_mH_gr;                   //  1.123xxe-5 cm^-3
    return;
}

//========================================================================================
int Cosmos :: recombination_history(int nzpts, double zstart, double zend, double *zarr, 
                                    double *Xe_H, double *Xe_He, double *Xe, double *TM)
{
  double *dXe=new double[nzpts];
  double *dX_H=new double[nzpts];

  int flg=recombination_history(nzpts, zstart, zend, zarr, Xe_H, Xe_He, Xe, dXe, dX_H, TM);

  delete[] dXe;
  delete[] dX_H;
  
  return flg;
}

//========================================================================================
double Hfunction(double z){ return GSL_cosmosptr->H(z); }

int Cosmos :: recombination_history(int nzpts, double zstart, double zend, 
                                    double *zarr, double *Xe_H, 
                                    double *Xe_He, double *Xe, double *dXe, 
                                    double *dX_H, double *TM)
{
    double *param=new double[14];

    // set parameters for Xe_frac()
    param[0]=nzpts;
    param[1]=zstart;
    param[2]=zend;
    param[3]=Y_p;
    param[4]=T_CMB0;
    param[5]=O_m;
    param[6]=O_b;
    param[7]=O_L;
    param[8]=O_k;
    param[9]=h100;
    param[10]=Nnu;
    param[11]=fudge;
    param[12]=fDM_annihilation;
    param[13]=CR_correction;

    if(Cosmos_mess_flag>=2)
    {
        cout << "\n Cosmos :: recombination_history: " << endl;
        cout << " running recombination_history with the following parameters: " << endl;   
        
        cout << " zs: " << param[1] << "\t ze: " << param[2] 
             << "\t nzpts: " << (int)param[0] << endl;
        cout << " Y_p: " << param[3] << "\t TCMB: " << param[4] << endl;
        cout << " OmegaM: " << param[5] << "\t OmegaB: " << param[6] 
             << "\n OmegaL: " << param[7] << "\t OmegaK: " << param[8] << endl;
        cout << " Omega_gamma h2: " << O_gh2 << "\t Omega_rel h2: " << O_relh2 
             << " Nnu= " << Nnu << endl;
        cout << " Hubble constant in units of 100 km s-1 Mpc-1: " << param[9] << endl;
        cout << " Fudge factor for H recombination: " << param[11] << endl;
        cout << " DM-annihilation " << param[12] << endl;
        cout << " Chluba & Thomas correction function: " << param[13] << endl;
    }
    
    // call with hubble function of cosmology object!
    GSL_cosmosptr=this;
    set_H_pointer(Hfunction);
    
    Xe_frac(param, zarr, Xe_H, Xe_He, Xe, dXe, dX_H, TM, Cosmos_mess_flag);

    reset_H_pointer();
    GSL_cosmosptr=NULL;
    
    delete [] param;
    return 0;
}

//========================================================================================
int Cosmos :: calc_coeff_X_spline(int M, double zstart, double zend)
{
    double *zarray=new double[M];
    double *zGSL=new double[M];
    double *Xe_H=new double[M];
    double *dXe_H=new double[M];
    double *Xe_He=new double[M];
    double *Xe=new double[M];
    double *dXe=new double[M];
    double *YGSL=new double[M];
    double *TM=new double[M];
    
    if(Cosmos_mess_flag>=1)
        cout << "\n Cosmos::calc_coeff_X_spline:" 
             << "\n calculating coeffients for X spline interpolation..." << endl;
    
    //====================================================================================
    // here arrays are filled with values
    //====================================================================================
    recombination_history(M, zsRe, zend, zarray, Xe_H, Xe_He, Xe, dXe, dXe_H, TM);
    
    //====================================================================================
    // resort zarray and Xe array so that ordering is zarray[0]<zarray[1]< ...
    //====================================================================================
    for(int j=0, i=M-1; i>=0; i--, j++) 
    {
        zGSL[j]=zarray[i];
        YGSL[j]=log(Xe[i]);
    }
    if(!spline_memory_allocated) 
        memindex_Xe=calc_spline_coeffies_JC(M, zGSL, YGSL, 
                                            "Cosmos :: calc_coeff_X_spline: Xe");
    else 
        update_spline_coeffies_JC(memindex_Xe, M, zGSL, YGSL, 
                                  "Cosmos :: calc_coeff_X_spline: Xe");

    //====================================================================================
    // dXe
    //====================================================================================
    for(int j=0, i=M-1; i>=0; i--, j++) YGSL[j]=dXe[i];
    if(!spline_memory_allocated) 
        memindex_dXe=calc_spline_coeffies_JC(M, zGSL, YGSL, 
                                             "Cosmos :: calc_coeff_X_spline: dXe");
    else 
        update_spline_coeffies_JC(memindex_dXe, M, zGSL, YGSL, 
                                  "Cosmos :: calc_coeff_X_spline: dXe");

    //====================================================================================
    // X1s
    //====================================================================================
    for(int j=0, i=M-1; i>=0; i--, j++) YGSL[j]=log(Xe_H[i]);
    if(!spline_memory_allocated) 
        memindex_XeH=calc_spline_coeffies_JC(M, zGSL, YGSL, 
                                             "Cosmos :: calc_coeff_X_spline: XeH");
    else 
        update_spline_coeffies_JC(memindex_XeH, M, zGSL, YGSL, 
                                  "Cosmos :: calc_coeff_X_spline: XeH");
                                   
    //====================================================================================
    // dXe_H
    //====================================================================================
    for(int j=0, i=M-1; i>=0; i--, j++) YGSL[j]=dXe_H[i];
    if(!spline_memory_allocated) 
        memindex_dXeH=calc_spline_coeffies_JC(M, zGSL, YGSL, 
                                              "Cosmos :: calc_coeff_X_spline: dXeH");
    else 
        update_spline_coeffies_JC(memindex_dXeH, M, zGSL, YGSL, 
                                  "Cosmos :: calc_coeff_X_spline: dXeH");
                                   
    //====================================================================================
    // rho=Te/Tg 
    //====================================================================================
    for(int j=0, i=M-1; i>=0; i--, j++) YGSL[j]=TM[i]/TCMB(zarray[i]);
    if(!spline_memory_allocated) 
        memindex_rho=calc_spline_coeffies_JC(M, zGSL, YGSL, 
                                             "Cosmos :: calc_coeff_X_spline: rho");
    else 
        update_spline_coeffies_JC(memindex_rho, M, zGSL, YGSL, 
                                  "Cosmos :: calc_coeff_X_spline: rho");
                                   
    //====================================================================================
    // XHe1s
    //====================================================================================
    for(int j=0, i=M-1; i>=0; i--, j++) YGSL[j]=log(Xe_He[i]);
    if(!spline_memory_allocated) 
        memindex_XeHe=calc_spline_coeffies_JC(M, zGSL, YGSL, 
                                              "Cosmos :: calc_coeff_X_spline: Xe_He");
    else 
        update_spline_coeffies_JC(memindex_XeHe, M, zGSL, YGSL, 
                                  "Cosmos :: calc_coeff_X_spline: Xe_He");
    
    //====================================================================================
    // clean up
    //====================================================================================
    delete[] zarray;
    delete[] zGSL;
    delete[] Xe_H;
    delete[] dXe_H;
    delete[] Xe_He;
    delete[] Xe;
    delete[] dXe;
    delete[] YGSL;
    delete[] TM;
    
    spline_memory_allocated=1;

    return 0;
}

//========================================================================================
// local spline-functions
//========================================================================================
void Cosmos :: check_splines()
{   
    if(splines_are_charged) return; 
    else
    { 
        cerr << " Cosmos::check_splines: recombination splines are not armed."
             << " Exiting. " << endl; 
        exit(0); 
    }
}

double Cosmos :: calc_spline(double z, int memindex)    
{   
    check_splines(); 
    return calc_spline_JC(z, memindex, " Cosmos::calc_spline "); 
}
double Cosmos :: calc_spline_log(double z, int memindex)
{ 
    check_splines(); 
    return exp(calc_spline_JC(z, memindex, " Cosmos::calc_spline_log ")); 
}

//========================================================================================
// public functions
//========================================================================================
void Cosmos :: display_cosmology()
{
    cout << "\n*----------------------------------------------"
         << "--------------------------------------------*" << endl;
    cout << "* Cosmos :: These are the most important cosmological"
         << " parameters used in this computation: *" << endl;
    cout << "*------------------------------------------------"
         << "------------------------------------------*\n" << endl;

    cout << " zs: " << zs << "\t\t ze: " << ze << endl;
    cout << " Y_p: " << Y_p << "\t\t fHe: " << fHe() << "\t\t TCMB: " << T_CMB0 << endl;
    cout << " OmegaM   : " << O_m   << "\t OmegaB   : " << O_b 
         << "\t OmegaL   : " << O_L   << "\t OmegaK   : " << O_k << endl;
    cout << " OmegaM h2: " << O_mh2 << "\t OmegaB h2: " << O_bh2 
         << "\t OmegaL h2: " << O_Lh2 << "\t OmegaK h2: " << O_kh2 << endl;
    cout << " Omega_gamma h2: " << O_gh2 << "\t Omega_rel h2: " << O_relh2 
         << "\t Nnu= " << Nnu << endl;
    cout << " Hubble constant in units of 100 km s-1 Mpc-1: " << h100 
         << " 1Mpc= " << const_Mpc << " cm " << endl;
    cout << " radiation-matter equality redshift z_eq: " << z_eq << endl;

    cout << "\n*----------------------------------------------"
         << "--------------------------------------------*" << endl;
    cout << "*------------------------------------------------"
         << "------------------------------------------*\n" << endl;

    return;
}

//========================================================================================
void Cosmos :: dump_cosmology(string filename)
{
    ofstream f;
    f.open(filename.c_str());
    f.precision(8);
    
    f << " zs: " << zs << "\t ze: " << ze << endl;
    f << " Y_p: " << Y_p << "\t fHe: " << fHe() << "\t TCMB: " << T_CMB0 << endl;
    f << " OmegaM   : " << O_m   << "\t OmegaB   : " << O_b 
      << "\t OmegaL   : " << O_L   << "\t OmegaK   : " << O_k << endl;
    f << " OmegaM h2: " << O_mh2 << "\t OmegaB h2: " << O_bh2 
      << "\t OmegaL h2: " << O_Lh2 << "\t OmegaK h2: " << O_kh2 << endl;
    f << " Omega_gamma h2: " << O_gh2 << "\t Omega_rel h2: " << O_relh2 
      << "\t Nnu= " << Nnu << endl;
    f << " Hubble constant in units of 100 km s-1 Mpc-1: " << h100 
      << " 1Mpc= " << const_Mpc << " cm " << endl;
    f << " radiation-matter equality redshift z_eq: " << z_eq << endl;

    f.close();

    return;
}

//========================================================================================
// Hubble-function
//========================================================================================
double Cosmos :: H(double z)
{
    if(loaded_Hz.check_limits(z)) return loaded_Hz.H(z); 
    
    //====================================================================================
    // This is equivalent with Seager 2000 ApJSS 128
    // return H0*sqrt(O_L+zp1*zp1*(O_k+zp1*O_m*(1.0+zp1/(1.0+z_eq)))); 
    // (Seager 2000 ApJSS 128)
    //====================================================================================
    double zp1=1.0+z;
    return H0*sqrt(O_L+zp1*zp1*(O_k+zp1*(O_m+O_rel*zp1)));       
}

//========================================================================================
// cosmological time
//========================================================================================
double dtdz_GSL(double logz, void *p)
{ 
    double zz=exp(logz); 
    return -zz*GSL_cosmosptr->dt_dz(zz);
}

double Cosmos :: t(double z)
{
    GSL_cosmosptr=this;
    double epsInt=1.0e-12, epsabs=0.0;
    double a=log(z), b=log(zs), r=0.0;
    
    r=Integrate_gk_GSL(6, a, b, epsInt, epsabs, dtdz_GSL);  

    GSL_cosmosptr=NULL;
    return r;
}

//========================================================================================
// Spitzer formula
//========================================================================================
double Cosmos :: t_ep(double z)                                                
{ 
    return t_eg(z)*sqrt(PI/2.0)*const_mp/(20.0*const_me)*(1.0-Y_p/2.0)
                  /(1.0-Y_p)*pow(TCMB(z)*(1.0+const_me/const_mp), 1.5); 
}

//========================================================================================
// integrals over Hubble; Integral^zh_z (1+z)^pow |dt_dz| dz
//========================================================================================
double dE_GSL(double logz, void *p)
{ 
    double z1=exp(logz); 
    return -z1*pow(z1, GSL_cosmosptr->Get_power())*GSL_cosmosptr->dt_dz(z1-1.0);
}

double Cosmos :: E_Int(double zh, double z, double pow)
{
    GSL_cosmosptr=this;
    power=pow;
    double epsInt=1.0e-12, epsabs=0.0;
    double a=log(z+1.0), b=log(zh+1.0), r=0.0;

    r=Integrate_gk_GSL(6, a, b, epsInt, epsabs, dE_GSL);

    GSL_cosmosptr=NULL;
    return r;
}

//========================================================================================
// helium nuclei
//========================================================================================
double Cosmos :: NHeII(double z)
{
  // it is assumed that NHeIII << NHeI & NHeII
  if(z<=4000) return Ne(z)-Np(z);

  // it is assumed that NHeI << NHeII & NHeIII
  return 2.0*fHe()*NH(z)+Np(z)-Ne(z);
}

double Cosmos :: NHeIII(double z)
{
  if(z<=4000) return 0.0;

  // it is assumed that NHeI << NHeII & NHeIII
  return Ne(z)-Np(z)-fHe()*NH(z);
}

double Cosmos :: NHeI(double z)
{ return fHe()*NH(z)-NHeII(z)-NHeIII(z); }

//========================================================================================
// compute Recfast system with a given initial solution at low z
//========================================================================================
int Cosmos :: recombine_using_Recfast_system(int nzpts, double zi, double ze, double fDM, 
                                             double Xe_Hi, double Xe_Hei, 
                                             double Xei, double TMi, double dXei, 
                                             const double *zarr, double *Xe_H, 
                                             double *Xe_He, double *Xe, double *TM)
{
    double *param=new double[14];

    // set parameters for Xe_frac()
    param[0]=nzpts;
    param[1]=zi;
    param[2]=ze;
    param[3]=Y_p;
    param[4]=T_CMB0;
    param[5]=O_m;
    param[6]=O_b;
    param[7]=O_L;
    param[8]=O_k;
    param[9]=h100;
    param[10]=Nnu;
    param[11]=fudge;
    param[12]=fDM;
    param[13]=0;
    
    if(Cosmos_mess_flag>=2)
    {
        cout << "\n Cosmos :: recombination_history: " << endl;
        cout << " running recombination_history with the following parameters: " << endl;   
        
        cout << " zs: " << param[1] << "\t ze: " << param[2] 
             << "\t nzpts: " << (int)param[0] << endl;
        cout << " Y_p: " << param[3] << "\t TCMB: " << param[4] << endl;
        cout << " OmegaM: " << param[5] << "\t OmegaB: " << param[6] 
             << "\n OmegaL: " << param[7] << "\t OmegaK: " << param[8] << endl;
        cout << " Omega_gamma h2: " << O_gh2 << "\t Omega_rel h2: " << O_relh2 
             << " Nnu= " << Nnu << endl;
        cout << " Hubble constant in units of 100 km s-1 Mpc-1: " << param[9] << endl;
        cout << " Fudge factor for H recombination: " << param[11] << endl;
        
        cout << "\n and initial values " << Xe_Hi << " " << Xe_Hei << " " 
             << Xei << " " << TMi << endl << endl;
    }
    
    // call with hubble function of cosmology object!
    GSL_cosmosptr=this;
    set_H_pointer(Hfunction);
                  
    Xe_frac_rescaled(param, zarr, Xe_H, Xe_He, Xe, TM, Xe_Hi, 
                     Xe_Hei, Xei, TMi, dXei, Cosmos_mess_flag);

    reset_H_pointer();
    GSL_cosmosptr=NULL;
    
    delete [] param;
    return 0;
}
