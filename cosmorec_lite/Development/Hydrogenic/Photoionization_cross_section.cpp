//======================================================================================
// Author: Jens Chluba 
// last modification: Oct 2010
// purpose: compute photoionization cross sections for hydrogenic atoms
// 
// The routines are based on Storey & Hummer, 1991. Notation and 
// recursion relations can be found there. However, we modified these 
// slightly to achieve better stability.
//======================================================================================

#include <iostream>
#include <string>
#include <cmath>
#include <stdlib.h>

#include "Photoionization_cross_section.h"
#include "routines.h"
#include "physical_consts.h"

using namespace std;

//======================================================================================
// global variables
//======================================================================================
double Photoionization_cross_section_thres_E=1.0e-8;
double Photoionization_cross_section_thres_nu=1.0e-8;
double Photoionization_cross_section_xi_max_limit_SH=1.0e+8;    // limit on nu/nu0. Beyond gnl=0


//======================================================================================
// Storey & Hummer, 1991
//======================================================================================
// kappa=sqrt(E)
//
//======================================================================================
double PCS_l_Ct(double E, int l){ return sqrt(1.0+l*l*E); } // Eq. 2.29 * l
double PCS_l_C(int n, int l){ return sqrt(1.0*n*n-l*l)/n; } // Eq. 2.25 * l

//======================================================================================
// formula 2.31
//======================================================================================
double PCS_Rt_0(int np)
{
  double log10r=0.5*(log10(PI/32.0)-log10factorial(2*np-1))+(np+2.0)*log10(4.0*np)-2.0*np*log10(exp(1.0));
  return pow(10.0, log10r);
}

//===============================================================================================   
// formula 2.32 in rewritten form
//===============================================================================================   
double ln_f_Boardman_E_SH_root(int i, double E, double y)
{
    if(i<1) return 1.0;
    
    double ff=pow((1.0+i*i*E)/(4.0*i), 2), r=log(1.0/ff/ff);
    for(int l=1; l<=i; l++) r+=log( (1.0+E*l*l)/ff/( 1.0*l*(2*i-(l-1)) ) );
    return r+log(2.0*i);
}

double PCS_Rt_stab_root_np(int np, double E)
{
    double y=sqrt(E);
    double logP=0.5*( ln_f_Boardman_E_SH_root(np, E, y) - 4.0*atan(np*y)/y - log( 1.0-exp( -2.0*PI/y ) ) ) ;    
    return exp(logP/np);  
}

//===============================================================================================   
// working version
//===============================================================================================   
void PCS_R_np(int np, int lstop, double E, vector<double> &Rnp_lplus1, vector<double> &Rnp_lminus1)
{
    //===============================================================================================   
    // JC 06.01.2010:
    //===============================================================================================   
    // Improving this routine was crucial for n>~300; the main problem was that R[n-1, n] for large 
    // n is very close to ~0 for large energies. However, due to ~n mutiplications the value for 
    // R[0, 1] ~ 10^-15. Numerically this was very unstable.
    // The solution was to rescale all R[l, l+1] by R[n-1, n]^((l+1)/n). In this way every multiplication 
    // increases the values of R a bit and the initial condition is R'[n-1, n]=1, which is numerically 
    // stable again. A similar method should be used for Anl_n'l'
    //===============================================================================================   
    
    // Rnp_lplus1 is array of length np+1 which will contain all the R[lp, lp+1] afterward
    // Rnp_lminus1 is array of length np+1 which will contain all the R[lp, lp-1] afterward
    
    //===============================================================================================   
    // NOTE: for l >= 0 the full array will be updated
    //       for l < 0 it is assumed that only the value for l is needed 
    //       and no further call for other l's and fixed energy E will follow
    //===============================================================================================   
    
    //===============================================================================================   
    // these levels don't exist --> no integrals R(lp, l)
    //===============================================================================================   
    Rnp_lplus1[np]=0.0;
    Rnp_lminus1[np]=0.0; 
    
    int lp=np-1;     // lp denotes the value of l for the free electron!
    //===============================================================================================   
    // formula 2.32 modified to obtain more stable results
    //===============================================================================================   
    double f=PCS_Rt_stab_root_np(np, E);
    double f1=f;
    double f2=pow(f1, 2);
    //
    Rnp_lplus1[lp]=sqrt(PI/32.0);
    
    //===============================================================================================       
    // from recurson formula 2.24
    //===============================================================================================       
    Rnp_lminus1[lp]=f2*PCS_l_Ct(E, lp+1)*Rnp_lplus1[lp]/( (lp+1.0)*2.0*PCS_l_Ct(E, lp) );
    
    for(lp=np-2; lp>=lstop; lp--)
    {
        //===================================================
        // formula 2.23
        //===================================================
        Rnp_lplus1[lp] =((2.0*(lp+1.0)+1.0)*PCS_l_Ct(E, lp+2)*Rnp_lplus1 [lp+1]*f1
                                           +PCS_l_C (np,lp+2)*Rnp_lminus1[lp+2])
                                           /( (lp+2.0)*2.0*PCS_l_C(np, lp+1) );
        //===================================================
        // formula 2.24      
        //===================================================
        Rnp_lminus1[lp]=(                   PCS_l_Ct(E, lp+1)*Rnp_lplus1[lp]*f2 
                              +(2.0*lp+1.0)*PCS_l_C (np,lp+1)*Rnp_lminus1[lp+1]*f1 ) 
                                           /( (lp+1.0)*2.0*PCS_l_Ct(E, lp) );
        
   }

    //===============================================================================================       
    // multiply an additional factor of sqrt(1+n^2*E)/n
    // see Eq. 2.20
    //===============================================================================================       
    double fac=sqrt(1.0+np*np*E)/np;
    for(lp=lstop; lp<np; lp++){ Rnp_lplus1[lp]*=fac*pow(f, lp+1); Rnp_lminus1[lp]*=fac*pow(f, lp-1); }
    
    return;
}

double PCS_sig_1(int l, double Rp1)
{ 
    //lp=l+1 > l
    return (l+1.0)*pow(Rp1, 2)/(2.0*l+1.0); // cm^2
}

double PCS_sig_2(int l, double Rm1)
{ 
    //lp=l-1 < l
    return l*pow(Rm1, 2)/(2.0*l+1.0); // cm^2
}



//======================================================================================
//
// Konstructors and Destructors for class Photoionization_cross_section
//
//======================================================================================
void Photoionization_cross_section::Set_Ry(double NNp)
{
    if(NNp<=0){ cout << " Photoionization_cross_section::Set_Ry: please check Np: " << NNp << endl; exit(0); }

    Np=NNp;
    mu_red=1.0/(1.0+const_me_mp/NNp);
    
    //==================================================================================
    // Rydberg of the atom: R_inf/(1+me/M)
    //==================================================================================
    Ry   =const_EH_inf_ergs*mu_red;
    Ry_Hz=const_EH_inf_Hz*mu_red;
    Ry_eV=const_EH_inf*mu_red;
    
    return;
}

void Photoionization_cross_section::init(int nv, int lv, int ZZ, double Np, int mflag)
{
    nn=nv; ll=lv; 
    Z=ZZ;
    Set_Ry(Np);
    nc=0; 
    
    if(ll>nn-1)
    {
        cout << " Photoionization_cross_section:: This is not a proper electron orbit: (n, l): " 
        << nn << " " << ll << endl;
        exit(0);
    }
    
    mess_flag=mflag;
    nu_ionization=nu_n(nn);
    E_ionization=E_n(nn);
    E_ionization_eV=E_n_eV(nn);
    
    return;
}

Photoionization_cross_section::Photoionization_cross_section(int nv, int lv, int Z, double Np, int mflag)
{ init(nv, lv, Z, Np, mflag); }

//======================================================================================
double Photoionization_cross_section::nu_threshold()
{ return nu_ionization*(1.0+Photoionization_cross_section_thres_nu);}

//======================================================================================
// Kramers formula as given in Karzas & Latter Eq. (39) comment:
// sig_Kramers ~ Z^4 relative to hydrogen value. This just comes from
// the scaling of nu_ion ~ Z^2
//======================================================================================
double Photoionization_cross_section::sig_phot_ion_Kramers(double nu)
{ 
    if(nu<nu_ionization) return 0.0;
    return const_PIe2_mec*16.0/( 3.0*sqrt(3.0)*PI*nn )*pow(nu_ionization/nu, 3)/nu_ionization; 
}



//===================================================================================
//
// Constructor and destructor for class Photoionization_cross_section_SH
//
//===================================================================================
void Photoionization_cross_section_SH::init(int nv, int lv, int Z, double Np, int mflag)
{
    Photoionization_cross_section::init(nv, lv, Z, Np, mflag);
    nc=1000; 

    if(nn>nc && mess_flag>0)
    {
        cout << " %==============================================================%" << endl;
        cout << " Photoionization_cross_section_SH::init:\n" 
             << " using Kramers-approximation for level ( " << nv << ", " << lv << " )" << endl << endl;
    }
    
    sigma_nucval=sig_phot_ion(nu_threshold());
    gaunt_nucval=g_phot_ion(nu_threshold());
    
    xi_small=nu_sig_phot_ion_small(1.0e-30)/Get_nu_ionization();

    return;
}

Photoionization_cross_section_SH::Photoionization_cross_section_SH(int nv, int lv, int Z, double Np, int mflag)
{ init(nv, lv, Z, Np, mflag); }


//===================================================================================
double Photoionization_cross_section_SH::sig_phot_ion(const vector<double> &Rp1, const vector<double> &Rm1)
{ 
    if(ll<0 || ll>=nn) return 0.0;  // state does not exist
    //===============================================================================
    // n and E are fixed when one has obtained *Rp1 and *Rm1
    // FOURPI*const_alpha*pow(const_a0, 2)/3.0=8.5596557e-19 according to S&H 1991
    //===============================================================================
    double f_sig=FOURPI*const_alpha*pow(const_a0, 2)/3.0;
    
    return f_sig*pow(1.0*Z*mu_red, -2)*(PCS_sig_1(ll, Rp1[ll])+PCS_sig_2(ll, Rm1[ll]));
}

double Photoionization_cross_section_SH::sig_phot_ion(double nu)
{
    if(nn>nc) return sig_phot_ion_Kramers(nu);
    if(nu<nu_ionization || nu>Photoionization_cross_section_xi_max_limit_SH*nu_ionization) return 0.0;
    
    //===============================================================================
    // to avoid the singularity at threshold energy
    //===============================================================================
    if(nu<=nu_threshold()) nu=nu_threshold();
    
    Rp1.resize(nn+1);
    Rm1.resize(nn+1);
    
    PCS_R_np(nn, ll, nu2E_e(nu), Rp1, Rm1);
    
    double r=sig_phot_ion(Rp1, Rm1);
    
    Rp1.clear();
    Rm1.clear();
    
    return r;
}

double Photoionization_cross_section_SH::g_phot_ion(double nu)
{ return sig_phot_ion(nu)/sig_phot_ion_Kramers(nu); }

//===================================================================================
// functions to find where the cross section becomes smaller than eps*sig_nuc
//===================================================================================
struct Photoionization_cross_section_SH_nu_sig_phot_ion_small
{
    double val;
    Photoionization_cross_section_SH *SH_ptr;
};

double Photoionization_cross_section_SH_func_nu2sig_phot_ion_small(double *nu, void *p)
{ 
    Photoionization_cross_section_SH_nu_sig_phot_ion_small *V=((Photoionization_cross_section_SH_nu_sig_phot_ion_small *) p); 
    return ( *nu)*( *nu)*V->SH_ptr->sig_phot_ion( *nu)/V->val-1.0; 
}

double Photoionization_cross_section_SH::nu_sig_phot_ion_small(double eps)
{
    double nu1=Get_nu_ionization();
    double nu2=Get_nu_ionization()*Photoionization_cross_section_xi_max_limit_SH;

    Photoionization_cross_section_SH_nu_sig_phot_ion_small Vars;
    Vars.SH_ptr=this;
    Vars.val=eps*nu1*nu1*sig_phot_ion_nuc();
    void *params=&Vars;
    
    double ff=Photoionization_cross_section_SH_func_nu2sig_phot_ion_small(&nu1, params)*Photoionization_cross_section_SH_func_nu2sig_phot_ion_small(&nu2, params);
    
    if(ff<0.0) return find_root_brent(Photoionization_cross_section_SH_func_nu2sig_phot_ion_small, params, nu1, nu2, 1.0e-3);
    else return Get_nu_ionization()*Photoionization_cross_section_xi_max_limit_SH;
}


//===================================================================================
//
// Konstructors and Destructors for class Photoionization_Lyc
//
//===================================================================================
void Photoionization_Lyc::init(int mflag)
{
    mess_flag=mflag;
    nu_ionization=nu_n();
    E_ionization=E_n();
    E_ionization_eV=E_n_eV();
    
    return;
}

Photoionization_Lyc::Photoionization_Lyc(int mflag)
{ init(mflag); }

//===================================================================================
double Photoionization_Lyc::nu_threshold()
{ return nu_ionization*(1.0+Photoionization_cross_section_thres_nu); }

//===================================================================================
// Kramers formula as given in Karzas & Latter Eq. (39) 
//===================================================================================
double Photoionization_Lyc::sig_phot_ion_Kramers(double nu)
{ 
    if(nu<nu_ionization) return 0.0;
    return const_PIe2_mec*16.0/( 3.0*sqrt(3.0)*PI )*pow(nu_ionization/nu, 3)/nu_ionization; 
}

double Photoionization_Lyc::sig_phot_ion(double nu)
{ return sig_phot_ion_Kramers(nu)*g_phot_ion(nu); }

double Photoionization_Lyc::g_phot_ion(double nu)
{ 
    if(nu<nu_ionization || nu>Photoionization_cross_section_xi_max_limit_SH*nu_ionization) return 0.0;
    //===================================================================================
    // to avoid the singularity at threshold energy
    //===================================================================================
    if(nu<=nu_threshold()) nu=nu_threshold();
    
    double E=nu2E_e(nu), y=sqrt(E);
    
    return 8.0*sqrt(3)*PI/(1.0+E)*exp(-4.0*atan(y)/y)/( 1.0-exp( -min(2.0*PI/y, 500.0) ) );
}

