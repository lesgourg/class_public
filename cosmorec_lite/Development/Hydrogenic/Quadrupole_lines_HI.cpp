//==================================================================================================
// Author: Jens Chluba 
//
// first implementation: May 2011
// last modification: May 2011
// 
// Routines to compute electic quadrupole transitions in hydrogenic 
// atoms. The formulae were taken from Grin et al. 2010 and Hey 2006.
// Also the formulae of Storey & Hummer 1991 are used.
//==================================================================================================

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdlib.h>

#include "routines.h"
#include "physical_consts.h"
#include "Quadrupole_lines_HI.h"

using namespace std;

//==================================================================================================
// reduced angular tensor operator C^(2) (Eq. (49) of Grin et al. 2010) 
//==================================================================================================
double reduced_angular_tensor_C2(int l, int lp)
{
    double sig = (l%2 == 0 ? 1.0 : -1.0);
    return sig*4.0*sqrt((0.5+l)*(0.5+lp))*Wigner_3J_Symbol(l, 2, lp, 0, 0, 0);
}


//==================================================================================================
// Storey & Hummer, 1991
// Recursion relatations for dipole matrix elements
//==================================================================================================
double Ql_C(int n, int l){ return (l>=n ? 1.0e-300 : sqrt(1.0*n*n-l*l)/n); }


//==================================================================================================
// initial condition for the dipole matrix element
//==================================================================================================
double Dipole_R0(int n, int np)
{
    //=============================================================================
    // n has to be bigger than np. 
    //=============================================================================
    double log10r=(np+2.0)*log10(4.0*np*n)-log10(4.0)+0.5*(log10factorial(n+np)
                                                           -log10factorial(n-np-1)
                                                           -log10factorial(2*np-1))
                                          +log10(1.0*n-np)*(n-np-2.0)
                                          -log10(1.0*n+np)*(n+np+2.0);
    
    return pow(10.0, log10r);
}

//==================================================================================================
double Quadrupole_R0_Dl1(int n, int np)
{
    double log10r=log10(np)+(np+2.0)*log10(np*n)+np*log10(4.0)+0.5*(log10factorial(n+np)
                                                                   -log10factorial(n-np-1)
                                                                   -log10factorial(2*np-1))
                                                +log10(1.0*n-np)*(n-np-3.0)
                                                -log10(1.0*n+np)*(n+np+3.0)
                                                +log10(1.0*np*np*(np+1.0)-n*n*(np - 3.0));
    
    return pow(10.0, log10r);
}    

//==================================================================================================
double Quadrupole_R0_Dlm2(int n, int np)
{
    double log10r=log10(np)+(np+3.0)*log10(np*n)+(np+2)*log10(4.0)+0.5*(log10factorial(n+np+1)
                                                                       -log10factorial(n-np-2)
                                                                       -log10factorial(2*np-1))
                                                +log10(1.0*n-np)*(n-np-3.0)
                                                -log10(1.0*n+np)*(n+np+3.0);
    
    return -pow(10.0, log10r);
}    


//==================================================================================================
// dipole matrix elements for hydrogenic atoms 
//==================================================================================================
void R_E1_Dipole(int n, int l, int np, int lp,     
                 vector<double> &Xnp_lm1_E1, 
                 vector<double> &Xnp_lp1_E1)
{
    //----------------------------------------------------------------------------------------------
    // initial condition for (1)^R^{n np}_{np np-1}
    //----------------------------------------------------------------------------------------------
    Xnp_lm1_E1[np]=Dipole_R0(n, np);
    
    //----------------------------------------------------------------------------------------------
    // from recurson formulae from Chluba et al. 2011
    //----------------------------------------------------------------------------------------------
    for(int lt=np-1; lt>=(int)max(0, lp-1); lt--)
    {
        //------------------------------------------------------------------------------------------
        // dipole elements
        //------------------------------------------------------------------------------------------
        Xnp_lp1_E1[lt] =( (2.0*lt+3.0)*Ql_C(np, lt+2)*Xnp_lp1_E1[lt+1]
                                      +Ql_C(n , lt+2)*Xnp_lm1_E1[lt+2] )
                        /(2.0*(lt+2.0)*Ql_C(n , lt+1));
        
        //------------------------------------------------------------------------------------------
        // dipole elements
        //------------------------------------------------------------------------------------------
        if(lt-1>=0) 
            Xnp_lm1_E1[lt]=(               Ql_C(np, lt+1)*Xnp_lp1_E1[lt] 
                             +(2.0*lt+1.0)*Ql_C(n , lt+1)*Xnp_lm1_E1[lt+1] )
                            /((lt+1.0)*2.0*Ql_C(np, lt));
    }
    
    return;
}

//==================================================================================================
// quadrupole matrix elements for hydrogenic atoms (Dl==0)
//==================================================================================================
void R_E2_Quadrupole_Dl0(int n, int l, int np, int lp,     
                         const vector<double> &Xnp_lm1_E1, 
                         vector<double> &Xnp_lmp_E2)
{
    //----------------------------------------------------------------------------------------------
    // from recurson formulae from Chluba et al. 2011
    // these elements are always needed when quadrupole lines are computed.
    //----------------------------------------------------------------------------------------------
    for(int lt=np-1; lt>=(int)max(0, lp-1); lt--)
    {
        //------------------------------------------------------------------------------------------
        // diagobal quadrupole element
        //------------------------------------------------------------------------------------------
        Xnp_lmp_E2[lt]=( Ql_C(np,lt+1)*Xnp_lmp_E2[lt+1] 
                         +2.0*(lt+1.0)*Xnp_lm1_E1[lt+1] )/ Ql_C(n, lt+1);
    }
    
    return;
}

//==================================================================================================
// quadrupole matrix elements for hydrogenic atoms (Dl==-2)
//==================================================================================================
void R_E2_Quadrupole_Dlm2(int n, int l, int np, int lp,     
                           const vector<double> &Xnp_lm1_E1, 
                           const vector<double> &Xnp_lmp_E2, 
                           vector<double> &Xnp_lm2_E2)
{
    //----------------------------------------------------------------------------------------------
    // quadrupole element l-2
    //----------------------------------------------------------------------------------------------
    if(lp==l-2 && n-np>1)  Xnp_lm2_E2[np+1]=-2.0*np*Ql_C(n, np+1)/Ql_C(n, np)*Xnp_lmp_E2[np-1];
    //if(lp==l-2 && n-np>1)  Xnp_lm2_E2[np+1]=Quadrupole_R0_Dlm2(n, np);
    
    for(int lt=np; lt>=lp+2; lt--)
    {
        Xnp_lm2_E2[lt]=( (2.0*lt-1.0)     *Ql_C(n ,lt+1)*Xnp_lm2_E2[lt+1]
                        + 2.0             *Ql_C(np,lt  )*Xnp_lmp_E2[lt]
                        + 2.0*(2.0*lt-1.0)*(3.0*lt+2.0) *Xnp_lm1_E1[lt] )
                       /( (2.0*lt+1.0)*Ql_C(np, lt-1) );  
    }

    return;
}

//==================================================================================================
// quadrupole matrix elements for hydrogenic atoms (Dl==+2)
//==================================================================================================
void R_E2_Quadrupole_Dlp2(int n, int l, int np, int lp,     
                          const vector<double> &Xnp_lp1_E1,
                          const vector<double> &Xnp_lmp_E2,
                          vector<double> &Xnp_lp2_E2)
{
    //----------------------------------------------------------------------------------------------
    // quadrupole element l+2
    //----------------------------------------------------------------------------------------------
    for(int lt=np-3; lt>=lp-2; lt--)
    {
        Xnp_lp2_E2[lt]=( (2.0*lt+3.0)     *Ql_C(np,lt+3)*Xnp_lp2_E2[lt+1]
                        + 2.0             *Ql_C(n ,lt+2)*Xnp_lmp_E2[lt+2]
                        + 2.0*(2.0*lt+3.0)*(3.0*lt+8.0) *Xnp_lp1_E1[lt+1] )
                       /( (2.0*lt+5.0)*Ql_C(n, lt+1) );   
    }
    
    return;
}
//==================================================================================================


//==================================================================================================
// quadrupole transition rates for hydrogen.
// (the mass of is assumed to be MH==mp, where mp is the proton mass)
// routines seems to work very well even for n~500
//==================================================================================================
double A_E2_Quadrupole(int n, int l, int np, int lp)
{ return A_E2_Quadrupole(1, 1.0, n, l, np, lp); }


//==================================================================================================
// quadrupole transition rates for hydrogenic atoms of 
// charge Z and mass Ma=Np*mp, where mp is the mass a proton. 
//==================================================================================================
double A_E2_Quadrupole(int Z, double Np, int n, int l, int np, int lp)
{
    //----------------------------------------------------------------------------------------------
    // NOTE: in recursions np, lp are used for the lower states
    // quadrupole transitions have l-lp= +/- 2 and 0 
    //----------------------------------------------------------------------------------------------
    if(!(l-lp==2 || lp-l==2 || lp-l==0) || np>=n) return 0.0; 
    
    //----------------------------------------------------------------------------------------------
    // no levels with l>=n
    //----------------------------------------------------------------------------------------------
    if(l>=n || lp>=np) return 0.0;
    
    //----------------------------------------------------------------------------------------------
    // no levels with l<0
    //----------------------------------------------------------------------------------------------
    if(l<0 || lp<0) return 0.0;
    if(n<2 || np<1) return 0.0;
    
    //----------------------------------------------------------------------------------------------
    // arrays for (1)^X^{n l}_{np l-1} and (1)^X^{n l}_{np l+1}
    //----------------------------------------------------------------------------------------------
    vector<double> Xnp_lm1_E1(np+4, 0.0);
    vector<double> Xnp_lp1_E1(np+4, 0.0);
    
    //----------------------------------------------------------------------------------------------
    // arrays for (2)^X^{n l}_{np l}, (2)^X^{n l}_{np l-2} and (2)^X^{n l}_{np l+2}
    //----------------------------------------------------------------------------------------------
    vector<double> Xnp_lmp_E2(np+4, 0.0);
    vector<double> Xnp_lm2_E2(np+4, 0.0);
    vector<double> Xnp_lp2_E2(np+4, 0.0);

    //----------------------------------------------------------------------------------------------
    // compute dipole matrix elements
    //----------------------------------------------------------------------------------------------
    R_E1_Dipole(n, l, np, lp, Xnp_lm1_E1, Xnp_lp1_E1);

    //----------------------------------------------------------------------------------------------
    // compute quadrupole matrix elements (Dl==0)
    //----------------------------------------------------------------------------------------------
    R_E2_Quadrupole_Dl0(n, l, np, lp, Xnp_lm1_E1, Xnp_lmp_E2);
    
    //----------------------------------------------------------------------------------------------
    // compute quadrupole matrix elements (Dl==-2)
    //----------------------------------------------------------------------------------------------
    if(lp==l-2) R_E2_Quadrupole_Dlm2(n, l, np, lp, Xnp_lm1_E1, Xnp_lmp_E2, Xnp_lm2_E2);
    
    //----------------------------------------------------------------------------------------------
    // compute quadrupole matrix elements (Dl==+2)
    //----------------------------------------------------------------------------------------------
    if(lp==l+2) R_E2_Quadrupole_Dlp2(n, l, np, lp, Xnp_lp1_E1, Xnp_lmp_E2, Xnp_lp2_E2);
    
    //----------------------------------------------------------------------------------------------
    double mu_red=1.0/(1.0+const_me_mp/Np);
    double red_C2=reduced_angular_tensor_C2(l, lp);
    
    //----------------------------------------------------------------------------------------------
    // formula according to Grin et al. 2010 (with Z and reduced mass added):
    //
    // norm_fac == const_alpha * mu_red * pow(Z*Z*TWOPI*nu_21_Hyd, 5) 
    //                         * pow(const_a0/const_cl/Z, 4) / ( 60.0*(0.5+l) );
    //
    // This can be rewritten using definitions of Ryd and a0.
    //----------------------------------------------------------------------------------------------
    double norm_fac=(const_cl/const_a0) * mu_red * pow(const_alpha*Z, 6) * pow(1.0/np/np-1.0/n/n, 5) 
                    / ( 15.0*256.0*(0.5+l) ); 
/*    
    cout << n << " " << l << " " << np << " " << lp << endl;
    for(int lt=np+1; lt>=(int)max(0, lp-2); lt--)
        cout << lt << " E2_0 " << Xnp_lmp_E2[lt] << " E2_m2 " << Xnp_lm2_E2[lt] << " E2_p2 " << Xnp_lp2_E2[lt] << endl;
    wait_f_r();
*/    
    if(lp==l)        red_C2*=Xnp_lmp_E2[l];
    else if(lp==l-2) red_C2*=Xnp_lm2_E2[l];
    else if(lp==l+2) red_C2*=Xnp_lp2_E2[l];
    
    return red_C2*red_C2*norm_fac;
}

//==================================================================================================
// dipole and quadrupole transition rates for hydrogenic atoms of 
// charge Z and mass Ma=Np*mp, where mp is the mass a proton. 
// Here all rates out of the level are computed at once to save time.
//==================================================================================================
void A_E1_E2_Quadrupole(int Z, double Np, int n, int l, int np,
                        double &Am1, double &Ap1, 
                        double &Am2, double &Amp, double &Ap2)
{
    Am1=Ap1=Am2=Amp=Ap2=0.0;
    
    //----------------------------------------------------------------------------------------------
    // NOTE: in recursions np, lp are used for the lower states.
    // dipole transitions have l-lp= +/- 1  
    // quadrupole transitions have l-lp= +/- 2 and 0 
    //----------------------------------------------------------------------------------------------
    if(np>=n || l>=n || l<0 || n<2 || np<1) return; 
    
    //----------------------------------------------------------------------------------------------
    // arrays for (1)^X^{n l}_{np l-1} and (1)^X^{n l}_{np l+1}
    //----------------------------------------------------------------------------------------------
    vector<double> Xnp_lm1_E1(np+4, 0.0);
    vector<double> Xnp_lp1_E1(np+4, 0.0);
    
    //----------------------------------------------------------------------------------------------
    // arrays for (2)^X^{n l}_{np l}, (2)^X^{n l}_{np l-2} and (2)^X^{n l}_{np l+2}
    //----------------------------------------------------------------------------------------------
    vector<double> Xnp_lmp_E2(np+4, 0.0);
    vector<double> Xnp_lm2_E2(np+4, 0.0);
    vector<double> Xnp_lp2_E2(np+4, 0.0);
    
    //----------------------------------------------------------------------------------------------
    // compute dipole matrix elements
    //----------------------------------------------------------------------------------------------
    R_E1_Dipole(n, l, np, l-2, Xnp_lm1_E1, Xnp_lp1_E1);
    
    //----------------------------------------------------------------------------------------------
    // compute quadrupole matrix elements (Dl==0, -/+2)
    //----------------------------------------------------------------------------------------------
    R_E2_Quadrupole_Dl0 (n, l, np, l-1, Xnp_lm1_E1, Xnp_lmp_E2);
    R_E2_Quadrupole_Dlm2(n, l, np, l-2, Xnp_lm1_E1, Xnp_lmp_E2, Xnp_lm2_E2);
    R_E2_Quadrupole_Dlp2(n, l, np, l+2, Xnp_lp1_E1, Xnp_lmp_E2, Xnp_lp2_E2);
    
    //----------------------------------------------------------------------------------------------
    // Dipole transitions
    //----------------------------------------------------------------------------------------------
    double mu_red=1.0/(1.0+const_me_mp/Np);
    
    double norm_fac_D=(const_cl/const_a0) * mu_red * pow(const_alpha*Z, 4) 
                       * pow(1.0/np/np-1.0/n/n, 3) / ( 12.0*(0.5+l) ); 
    
    if(l-1>=0) Am1=norm_fac_D*l      *pow(Xnp_lm1_E1[l], 2);
    if(l+1<np) Ap1=norm_fac_D*(l+1.0)*pow(Xnp_lp1_E1[l], 2);

    //----------------------------------------------------------------------------------------------
    // formula according to Grin et al. 2010 (with Z and reduced mass added):
    //
    // norm_fac == const_alpha * mu_red * pow(Z*Z*TWOPI*nu_21_Hyd, 5) 
    //                         * pow(const_a0/const_cl/Z, 4) / ( 60.0*(0.5+l) );
    //
    // This can be rewritten using definitions of Ryd and a0.
    //----------------------------------------------------------------------------------------------
    double norm_fac_Q=(const_cl/const_a0) * mu_red * pow(const_alpha*Z, 6) 
                       * pow(1.0/np/np-1.0/n/n, 5) / ( 15.0*256.0*(0.5+l) ); 

    if(l-2>=0)       Am2=norm_fac_Q*pow(reduced_angular_tensor_C2(l, l-2)*Xnp_lm2_E2[l], 2);
    if(l!=0 && l<np) Amp=norm_fac_Q*pow(reduced_angular_tensor_C2(l, l)  *Xnp_lmp_E2[l], 2);
    if(l+2<np)       Ap2=norm_fac_Q*pow(reduced_angular_tensor_C2(l, l+2)*Xnp_lp2_E2[l], 2);
    
    return;
}

//==================================================================================================
//==================================================================================================
