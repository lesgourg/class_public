//=================================================================================
// Author: Jens Chluba 
// last modification: May 2011
// 
// The routines are based on Storey & Hummer, 1991. Notation and 
// recursion relations can be found there. However, we modified these 
// slightly to achieve better stability.
//=================================================================================

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdlib.h>

#include "routines.h"
#include "physical_consts.h"
#include "Oscillator_strength.h"

using namespace std;

//=================================================================================
// Storey & Hummer, 1991
// Recursion relatations for A(n,l-->n',l')
//=================================================================================
double l_C(int n, int l){ return sqrt(1.0*n*n-l*l)/n; }             // Eq. 2.25 * l

double R0_root(int n, int np)
{
    //=============================================================================
    // n has to be bigger than np
    //=============================================================================
    double log10r=(np+2.0)*log10(4.0*np*n)-log10(4.0)+0.5*(log10factorial(n+np)
                                                          -log10factorial(n-np-1)
                                                          -log10factorial(2*np-1))
                                                          +log10(1.0*n-np)*(n-np-2.0)
                                                          -log10(1.0*n+np)*(n+np+2.0);
    return pow(10.0, log10r/np);
}

//=================================================================================
// main function
//=================================================================================
double A_SH_func(int &n, int &l, int &np, int &lp)
{
    //=============================================================================
    // NOTE: in recursions np, lp are used for the lower states
    // dipole transitions l-lp= +/- 1 
    //=============================================================================
    if(abs(l-lp)!=1.0) return 0.0;
    if(n<=np) return 0.0;

    //=============================================================================
    // no levels with l>=n
    //=============================================================================
    if(l>=n || lp>=np) return 0.0;

    //=============================================================================
    // no levels with l<0
    //=============================================================================
    if(l<0 || lp<0) return 0.0;
    if(n<2 || np<1) return 0.0;
    
    vector<double> Rnp_lplus1(np+1);
    vector<double> Rnp_lminus1(np+1);
    
    //=============================================================================
    // these levels don't exist --> no integrals R(lp, l)
    //=============================================================================
    Rnp_lplus1[np]=0.0;
    Rnp_lminus1[np]=0.0; 
    
    int lt=np-1;
    //=============================================================================
    // formula 2.27
    //=============================================================================
    double f=R0_root(n, np);
    double f1=f;
    double f2=pow(f1, 2);
    //
    Rnp_lplus1[lt]=1.0;

    //=============================================================================
    // from recurson formula 2.24
    //=============================================================================
    Rnp_lminus1[lt]=f2*l_C(n, lt+1)*Rnp_lplus1[lt]/( (lt+1.0)*2.0*l_C(n, lt) );
    
    for(lt=np-2; lt>=lp; lt--)
    {
        //=========================================================================
        // formula 2.23
        //=========================================================================
        Rnp_lplus1[lt] =((2.0*lt+3.0)*l_C(n, lt+2)*Rnp_lplus1 [lt+1]*f1 
                                     +l_C(np,lt+2)*Rnp_lminus1[lt+2] )
                        /( (lt+2.0)*2.0*l_C(np, lt+1) );

        //=========================================================================
        // formula 2.24
        //=========================================================================
        Rnp_lminus1[lt]=(                l_C(n, lt+1)*Rnp_lplus1[lt]   *f2 
                           +(2.0*lt+1.0)*l_C(np,lt+1)*Rnp_lminus1[lt+1]*f1 )
                        /( (lt+1.0)*2.0*l_C(n, lt) );
    }
    
    //=============================================================================
    // see Eq. 2.21, hydrogen
    //=============================================================================
    double mu_red=1.0/(1.0+const_me_mp), Z=1.0; 
    double fac_A=pow(const_alpha, 4)*const_cl/6.0/const_a0;  // 2.6775015e+9 according to S&H 1991
    double fac=fac_A*mu_red*pow(Z, 4)*pow(1.0/np/np-1.0/n/n, 3);
    
    //=============================================================================
    //lp=l+1 > l
    //=============================================================================
    if(lp==l+1) return fac*(l+1.0)*pow(Rnp_lminus1[lp]*pow(f, lp-1), 2)/(2.0*l+1.0);   // 1/sec

    //=============================================================================
    //lp=l-1 < l
    //=============================================================================
    if(lp==l-1) return fac*l*pow(Rnp_lplus1[lp]*pow(f, lp+1), 2)/(2.0*l+1.0);          // 1/sec
    
    return 0.0;
}

//=================================================================================
double A_SH(int &n, int &l, int &np, int &lp)
{ 
    double A=A_SH_func(n, l, np, lp); 

    if(isnan_JC(A))
    {
        cout << " A_SH:: detected nan " << A << " || " << n << " " << l 
             << " <--> " << np << " " << lp << endl;
        wait_f_r();
    }
    
    return A;
}

//=================================================================================
double f_SH(int &n, int &l, int &np, int &lp)
{
    double A;
    if(n<np) A=A_SH(np, lp, n, l);
    else A=A_SH(n, l, np, lp);
    
    double lambda=fabs(const_cl/( const_EH_inf_Hz*(1.0/n/n-1.0/np/np)/(1.0+const_me_mp) ));
    double f=A*pow(lambda, 2)/(2.0*FOURPI*const_PIe2_mec*(1.0+const_me_mp));
    
    if(n<np) return f*(2.0*lp+1.0)/(2.0*l+1.0);
    else return f;
}

//=================================================================================
// oscillator strength and transition rates for 
// hydrogenic atoms of charge Z and mass Ma=Np*mp, where mp is the mass
// a proton. Drake 1996, Atomic, Molecular & Optical Physics
//=================================================================================
double A_SH(int Z, double Np, int n, int l, int np, int lp)
{ 
    double A=pow(1.0*Z, 4)*A_SH(n, l, np, lp)*(1.0+const_me_mp)/(1.0+const_me_mp/Np); 
    
    if(isnan_JC(A))
    {
        cout << " A_SH:: detected nan " << A << " || " << n << " " << l 
             << " <--> " << np << " " << lp << endl;
        wait_f_r();
    }
    
    return A;
}

double f_SH(int Z, double Np, int n, int l, int np, int lp){ return f_SH(n, l, np, lp); }



