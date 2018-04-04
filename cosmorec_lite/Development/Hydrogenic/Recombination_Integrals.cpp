//================================================================================================
// Author: Jens Chluba 
// last modification: Oct 2010
// 
// Purpose: carry out the integrals over photoionization cross sections
//================================================================================================

#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>

#include "routines.h"
#include "Definitions.h"
#include "Integration_routines.h"

using namespace std;

//================================================================================================
// photonionization rate. 
//================================================================================================

//================================================================================================
// Integral using Patterson formulae
//================================================================================================
int Integrate_Ric_Rci_JC_Patterson_refine(int &reclev, double lga, double lgb, double epsrel, double epsabs, double (*fptr)(double), double *r)
{
    int debug=0;
    *r=0.0;
    
    //-------------------------------------------------------------------
    // to avoid SEGfault for strange solution
    //-------------------------------------------------------------------
    reclev++;
    if(reclev>=100){ cout << " Integrate_Ric_Rci_JC_Patterson_refine::recursion_level has exceeded limit of 100. Exiting \n" << endl; exit(0); }
        
    if(debug>=1) cout << " entering subdivision " << exp(lga) << " " << exp(lgb)<< endl;
    int neval, ntot=0, ifail=0;
    double r1;
    
    int subint=(int)max(2, 4*(lgb-lga)/log(10.0)); // about 4 intervals per decade
    double a1, b1=lga, fa=(lgb-lga)/subint, err;
    
    for(int k=1; k<=subint; k++)
    {
        a1=b1; 
        b1=a1+fa; 
        r1=0.0;
        err=max(epsabs, fabs( *r)*epsrel);
        
        ifail=Integrate_using_Patterson(a1, b1, epsrel, err, fptr, &neval, &r1);
        ntot+=neval;
        
        //===========================================================================================
        // do subdivision (recursive calling)
        //===========================================================================================
        if(ifail!=0)
        { 
            neval=Integrate_Ric_Rci_JC_Patterson_refine(reclev, a1, b1, epsrel, err, fptr, &r1); 
            ifail=0; 
            ntot+=neval;
        }
        
        *r+=r1;

        if(debug>=2) cout << " refined Ric_Rci " << k << " " << r1 << " " << *r << " ratio= " << fabs(r1/ *r) 
                          << " ( " << a1 << " " << b1 << " ) xmax= " << lgb << " " << neval << " " << ifail << endl;
        if(fabs(r1)*(subint-k)<=err) break;
    }
    if(debug>=1) cout << " refined Ric_Rci= " << *r << " " << ntot << " " << subint << endl;
    if(debug>=2) wait_f_r();
    
    return ntot;
}   

//================================================================================================
double Integrate_Ric_Rci_JC_Patterson_part(int *stopped_before_end_of_interval, double lga, double lgb, double epsrel, double epsabs, double (*fptr)(double))
{
    if(lga>=lgb) return 0.0;
    
    int debug=0;
    int neval, ntot=0, ifail=0;
    double r=0.0, r1;
    
    int subint=(int)max(2, 4*(lgb-lga)/log(10.0)); // about 4 intervals per decade
    double a1, b1=lga, fa=(lgb-lga)/subint, err;
    
    int k;

    //-------------------------------------------------------------------
    // to avoid SEGfault for strange solution
    //-------------------------------------------------------------------
    int reclev=0;

    for(k=1; k<=subint; k++)
    {
        a1=b1;
        b1=a1+fa;
        r1=0.0;
        err=max(epsabs, fabs(r)*epsrel);
        
        ifail=Integrate_using_Patterson(a1, b1, epsrel, err, fptr, &neval, &r1);
        ntot+=neval;
        
        //===========================================================================================
        // do subdivision
        //===========================================================================================
        if(ifail!=0)
        { 
            neval=Integrate_Ric_Rci_JC_Patterson_refine(reclev, a1, b1, epsrel, err, fptr, &r1); 
            ifail=0; 
            ntot+=neval;
        }
        
        r+=r1;
        
        if(debug>=2) cout << " Ric_Rci " << k << " " << r1 << " " << r << " ratio= " << fabs(r1/r) 
                          << " ( " << a1 << " " << b1 << " ) xmax= " << lgb << " " << neval << " " << ifail << endl;
        if(fabs(r1)*(subint-k)<=err) break;
    }
    if(debug>=1) cout << " Ric_Rci= " << r << " " << ntot << " " << subint << " " << k << "\n out"<< endl;
    if(debug>=2) wait_f_r();
    
    if(k<subint) *stopped_before_end_of_interval=1;
        
    return r;
}

//================================================================================================
double Integrate_Ric_Rci_JC_Patterson(int choice, double lga, double lgb, double epsrel, double epsabs, double (*fptr)(double))
{   
    int debug=0;
    //===========================================================================================================
    // in both cases, Ric & Rci for x>>1 things cut off exponentially due to Planck & Maxwell Boltzmann.
    // Therefore one case but an upper limit on x. Alternatively, one can integrate the rates using Gauss-Laguere 
    // rules, but this means one has to transform from log(x) --> x (this does not converge yet...)
    //===========================================================================================================
    double xWien=4.0;
    double xmax=xWien-log(1.0e-30);
    double xl=exp(lga), xu=exp(lgb);
    //
    double a, b;
    double r=0.0, r1=0.0, errabs;
    int stopped_before_end_of_interval=0;
    
    //-------------------------------------------------------------------
    // Only integrate over Wien-tail and return
    //-------------------------------------------------------------------
    if(xl>xWien) 
    {
        a=lga; b=lgb;
        r=Integrate_Ric_Rci_JC_Patterson_part(&stopped_before_end_of_interval, a, b, epsrel, epsabs, fptr);

        return r;
    }
    
    //-------------------------------------------------------------------
    // Integrate both over Wien-tail and low frequency part. 
    // Start with low nu part and stop when contributions become small.
    //-------------------------------------------------------------------
    stopped_before_end_of_interval=0;

    a=lga; b=log(xWien);
    r=Integrate_Ric_Rci_JC_Patterson_part(&stopped_before_end_of_interval, a, b, epsrel, epsabs, fptr);
    // 
    if(stopped_before_end_of_interval==0) // only continue if the other part was not already successful
    {
        a=b; b=min(lgb, log(xmax));
        errabs=min(fabs(r)*epsrel, epsabs);

        if(debug>=1) cout << " Integrate_Ric_Rci_JC_Patterson: " << exp(a) << " " << exp(b) << " " << xu << endl;
        if(debug>=2) wait_f_r();
        
        r1=Integrate_Ric_Rci_JC_Patterson_part(&stopped_before_end_of_interval, a, b, epsrel, errabs, fptr);
        r+=r1;
    }

    return r;
}

//================================================================================================
// Integral using Patterson formulae with parameter
//================================================================================================
int Integrate_Ric_Rci_JC_Patterson_refine(int &reclev, double lga, double lgb, double epsrel, double epsabs, double (*fptr)(double, void *p), double *r, void *p)
{
    int debug=0;
    *r=0.0;
    
    //-------------------------------------------------------------------
    // to avoid SEGfault for strange solution
    //-------------------------------------------------------------------
    reclev++;
    if(reclev>=100){ cout << " Integrate_Ric_Rci_JC_Patterson_refine::recursion_level has exceeded limit of 100. Exiting \n" << endl; exit(0); }
    
    if(debug>=1) cout << " entering subdivision " << exp(lga) << " " << exp(lgb)<< endl;
    int neval, ntot=0, ifail=0;
    double r1;
    
    int subint=(int)max(2, 4*(lgb-lga)/log(10.0)); // about 4 intervals per decade
    double a1, b1=lga, fa=(lgb-lga)/subint, err;
    
    for(int k=1; k<=subint; k++)
    {
        a1=b1; 
        b1=a1+fa; 
        r1=0.0;
        err=max(epsabs, fabs( *r)*epsrel);
        
        ifail=Integrate_using_Patterson(a1, b1, epsrel, err, fptr, &neval, &r1, p);
        ntot+=neval;
        
        //===========================================================================================
        // do subdivision (recursive calling)
        //===========================================================================================
        if(ifail!=0)
        { 
            neval=Integrate_Ric_Rci_JC_Patterson_refine(reclev, a1, b1, epsrel, err, fptr, &r1, p); 
            ifail=0; 
            ntot+=neval;
        }
        
        *r+=r1;
        
        if(debug>=2) cout << " refined Ric_Rci " << k << " " << r1 << " " << *r << " ratio= " << fabs(r1/ *r) 
                          << " ( " << a1 << " " << b1 << " ) xmax= " << lgb << " " << neval << " " << ifail << endl;
        if(fabs(r1)*(subint-k)<=err) break;
    }
    if(debug>=1) cout << " refined Ric_Rci= " << *r << " " << ntot << " " << subint << endl;
    if(debug>=2) wait_f_r();
    
    return ntot;
}   

//================================================================================================
double Integrate_Ric_Rci_JC_Patterson_part(int *stopped_before_end_of_interval, double lga, double lgb, double epsrel, double epsabs, double (*fptr)(double, void *p), void *p)
{
    if(lga>=lgb) return 0.0;
    
    int debug=0;
    int neval, ntot=0, ifail=0;
    double r=0.0, r1;
    
    int subint=(int)max(2, 4*(lgb-lga)/log(10.0)); // about 4 intervals per decade
    double a1, b1=lga, fa=(lgb-lga)/subint, err;
    
    int k;
    
    //-------------------------------------------------------------------
    // to avoid SEGfault for strange solution
    //-------------------------------------------------------------------
    int reclev=0;
    
    for(k=1; k<=subint; k++)
    {
        a1=b1;
        b1=a1+fa;
        r1=0.0;
        err=max(epsabs, fabs(r)*epsrel);
        
        ifail=Integrate_using_Patterson(a1, b1, epsrel, err, fptr, &neval, &r1, p);
        ntot+=neval;
        
        //===========================================================================================
        // do subdivision
        //===========================================================================================
        if(ifail!=0)
        { 
            neval=Integrate_Ric_Rci_JC_Patterson_refine(reclev, a1, b1, epsrel, err, fptr, &r1, p); 
            ifail=0; 
            ntot+=neval;
        }
        
        r+=r1;
        
        if(debug>=2) cout << " Ric_Rci " << k << " " << r1 << " " << r << " ratio= " << fabs(r1/r) 
                          << " ( " << a1 << " " << b1 << " ) xmax= " << lgb << " " << neval << " " << ifail << endl;
        if(fabs(r1)*(subint-k)<=err) break;
    }
    if(debug>=1) cout << " Ric_Rci= " << r << " " << ntot << " " << subint << " " << k << "\n out"<< endl;
    if(debug>=2) wait_f_r();
    
    if(k<subint) *stopped_before_end_of_interval=1;
    
    return r;
}

//================================================================================================
double Integrate_Ric_Rci_JC_Patterson(int choice, double lga, double lgb, double epsrel, double epsabs, double (*fptr)(double, void *p), void *p)
{   
    int debug=0;
    //===========================================================================================================
    // in both cases, Ric & Rci for x>>1 things cut off exponentially due to Planck & Maxwell Boltzmann.
    // Therefore one case but an upper limit on x. Alternatively, one can integrate the rates using Gauss-Laguere 
    // rules, but this means one has to transform from log(x) --> x (this does not converge yet...)
    //===========================================================================================================
    double xWien=4.0;
    double xmax=xWien-log(1.0e-30);
    double xl=exp(lga), xu=exp(lgb);
    //
    double a, b;
    double r=0.0, r1=0.0, errabs;
    int stopped_before_end_of_interval=0;
    
    //-------------------------------------------------------------------
    // Only integrate over Wien-tail and return
    //-------------------------------------------------------------------
    if(xl>xWien) 
    {
        a=lga; b=lgb;
        r=Integrate_Ric_Rci_JC_Patterson_part(&stopped_before_end_of_interval, a, b, epsrel, epsabs, fptr, p);
        
        return r;
    }
    
    //-------------------------------------------------------------------
    // Integrate both over Wien-tail and low frequency part. 
    // Start with low nu part and stop when contributions become small.
    //-------------------------------------------------------------------
    stopped_before_end_of_interval=0;
    
    a=lga; b=log(xWien);
    r=Integrate_Ric_Rci_JC_Patterson_part(&stopped_before_end_of_interval, a, b, epsrel, epsabs, fptr, p);
    //  cout << r << endl; 
    // 
    if(stopped_before_end_of_interval==0) // only continue if the other part was not already successful
    {
        a=b; b=min(lgb, log(xmax));
        errabs=min(fabs(r)*epsrel, epsabs);
        
        if(debug>=1) cout << " Integrate_Ric_Rci_JC_Patterson: " << exp(a) << " " << exp(b) << " " << xu << endl;
        if(debug>=2) wait_f_r();
        
        r1=Integrate_Ric_Rci_JC_Patterson_part(&stopped_before_end_of_interval, a, b, epsrel, errabs, fptr, p);
        r+=r1;
    }
    
    return r;
}



//================================================================================================
// Integral using Chebyshev formulae
//================================================================================================
double Integrate_Ric_Rci_JC_Chebyshev(int choice, double xl, double xu, double epsrel, double epsabs, double (*fptr)(double), double xc)
{   
    int debug=0;
    //===========================================================================================================
    // in both cases, Ric & Rci for x>>1 things cut off exponentially due to Planck & Maxwell Boltzmann.
    // Therefore one can put an upper limit on x. 
    //===========================================================================================================
    double xWien=2.0;
    double xmax=xWien-log(1.0e-40);
    //
    double a, b, x0;
    double r=0.0;
    int neval=0;
    
    //======================================================
    // Use Gauss-Laguerre in Wien tail
    //======================================================
    if(xl>=xWien && choice!=3)
    {
        a=xl; b=min(xu, xmax); x0=min((3.0*a+b)/4.0, 1.2*a);
        
        Integrate_using_Chebyshev(a, x0, b, epsrel, epsabs, fptr, &neval, &r);

        if(debug>=2) cout << r << " " << neval << " " << a << " " << b << " " << x0 << " choice= " << choice << " " << endl;
        if(debug>=3) wait_f_r(" Wien");
        return r;
    }
    
    //======================================================
    // When interval is small, just do simple thing
    //======================================================
    if(xu-xl<=5)
    {
        a=xl; b=xu; x0=min((3.0*a+b)/8.0, 2.0*a);
        Integrate_using_Chebyshev(a, x0, b, epsrel, epsabs, fptr, &neval, &r);
        
        if(debug>=2) cout << r << " " << neval << " " << a << " " << b << " " << x0 << " choice= " << choice << " " << endl;
        if(debug>=3) wait_f_r(" short int ");
        return r;
    }
        
    //======================================================
    // case for dRci_dTe
    //======================================================
    if(choice==3) 
    {
        a=xl; b=min(xu, xmax); x0=min((3.0*a+b)/4.0, xc);
        Integrate_using_Chebyshev(a, x0, b, epsrel, epsabs, fptr, &neval, &r);

        if(debug>=2) cout << r << " " << neval << " " << a << " " << b << " " << x0 << " choice= " << choice << " " << endl;
        if(debug>=3) wait_f_r(" dRci_dTe ");
        return r;
    }
    
    //======================================================
    // other cases
    //======================================================    
    if(xu<=xWien)
    {
        a=xl; b=xu; x0=min((7.0*a+b)/8.0, 20.0*a);
        Integrate_using_Chebyshev(a, x0, b, epsrel, epsabs, fptr, &neval, &r);
        
        if(debug>=2) cout << r << " " << neval << " " << a << " " << b << " " << x0 << " choice= " << choice << " " << endl;
        if(debug>=3) wait_f_r(" no wien");
    }
    else 
    {
        a=xl; b=min(xu, xmax); x0=min((3.0*a+b)/4.0, 10.0*a);
        Integrate_using_Chebyshev(a, x0, b, epsrel, epsabs, fptr, &neval, &r);
        
        if(debug>=2) cout << r << " " << neval << " " << a << " " << b << " " << x0 << " choice= " << choice << " " << endl;
        if(debug>=3) wait_f_r(" both ");
    }
    
    return r;
}

//================================================================================================
// Integral using Chebyshev formulae in log space
//================================================================================================
inline void message_Chebyshev(string mess, double r, int neval, double a, double x0, double b, int choice, int debug)
{
    if(debug>=1)
    {
        cout << r << mess << " neval= " << neval << " " << a << " " << x0 << " " << b << " choice= " << choice << " " << endl;
        if(debug>=2) wait_f_r(mess);
    }
    
    return;
}

//================================================================================================
double Integrate_Ric_Rci_JC_Chebyshev_log(int choice, string mess, double lgxl, double lgxu, double epsrel, double epsabs, double (*fptr)(double), double lgxc)
{   
    int debug=0;
    //===========================================================================================================
    // in both cases, Ric & Rci for x>>1 things cut off exponentially due to Planck & Maxwell Boltzmann.
    // Therefore one can put an upper limit on x. 
    //===========================================================================================================
    double xWien=5.0;
    double xl=exp(lgxl), xu=exp(lgxu), xc=exp(lgxc);
    //
    double a, b, x0;
    double k1=7.0, k2=8.0;
    double r=0.0;
    int neval=0;
    
    //======================================================
    // Use Gauss-Laguerre in Wien tail
    //======================================================
    if(xl>=xWien && choice!=3)
    {
        a=xl; b=xu; x0=min((k1*a+b)/k2, 1.5*a);

        Integrate_using_Chebyshev(log(a), log(x0), log(b), epsrel, epsabs, fptr, &neval, &r, " log-Wien"+mess);
        
        message_Chebyshev(" log-Wien"+mess, r, neval, a, x0, b, choice, debug);

        return r;
    }
    
    //======================================================
    // When interval is small, just do simple thing
    //======================================================
    if(xu/xl<=5)
    {
        a=xl; b=xu; x0=min((k1*a+b)/k2, 1.2*a);
        Integrate_using_Chebyshev(log(a), log(x0), log(b), epsrel, epsabs, fptr, &neval, &r, " log-short"+mess+" "+int_to_string(choice));
        
        message_Chebyshev(" log-short"+mess, r, neval, a, x0, b, choice, debug);
        return r;
    }
    
    //======================================================
    // case for dRci_dTe
    //======================================================
    if(choice==3) 
    {
        if(xu>xc){ a=xl; b=xu; x0=min((k1*a+b)/k2, xc); }
        else{ a=xl; b=xu; x0=min((k1*a+b)/k2, 1.2*a); } 
        Integrate_using_Chebyshev(log(a), log(x0), log(b), epsrel, epsabs, fptr, &neval, &r, " log-dRci_dTe"+mess);

        message_Chebyshev(" log-dRci_dTe "+mess, r, neval, a, x0, b, choice, debug);
        return r;
    }
    
    //======================================================
    // other cases
    //======================================================    
    if(xu<=xWien)
    {
        double r1=0.0;
        int n1=0;
        a=xl; b=min(xl*5.0, xu); x0=min((k1*a+b)/k2, 2.0*a);
        Integrate_using_Chebyshev(log(a), log(x0), log(b), epsrel, epsabs, fptr, &neval, &r1, " log-no wien I"+mess);
        n1=neval; r+=r1;
        
        a=b; b=min(xl*10.0, xu); x0=min((k1*a+b)/k2, 2.0*a); r1=0.0; neval=0;
        Integrate_using_Chebyshev(log(a), log(x0), log(b), epsrel, max(epsabs, fabs(r)*epsrel), fptr, &neval, &r1, " log-no wien II"+mess);
        n1+=neval; r+=r1;
        
        a=b; b=xu; x0=min((k1*a+b)/k2, 3.0*a);
        Integrate_using_Chebyshev(log(a), log(x0), log(b), epsrel, max(epsabs, fabs(r)*epsrel), fptr, &neval, &r1, " log-no wien III "+mess);
        neval+=n1; r+=r1;

        message_Chebyshev(" log-no Wien I&II&III"+mess, r, neval, a, x0, b, choice, debug);
    }
    else 
    {
        double r1=0.0;
        int n1=0;
        a=xl; b=min(xl*20.0, xu); x0=min((k1*a+b)/k2, 3.0*a);
        Integrate_using_Chebyshev(log(a), log(x0), log(b), epsrel, epsabs, fptr, &neval, &r1, " log-both I"+mess);
        n1=neval; r+=r1;
        
        a=b; b=min(xWien, xu); x0=min((k1*a+b)/k2, 10.0*a); r1=0.0; neval=0;
        Integrate_using_Chebyshev(log(a), log(x0), log(b), epsrel, max(epsabs, fabs(r)*epsrel), fptr, &neval, &r1, " log-both II"+mess);
        n1+=neval; r+=r1;
        
        a=b; b=xu; x0=min((k1*a+b)/k2, 1.2*a); r1=0.0; neval=0;
        Integrate_using_Chebyshev(log(a), log(x0), log(b), epsrel, max(epsabs, fabs(r)*epsrel), fptr, &neval, &r1, " log-both III"+mess);
        neval+=n1; r+=r1;

        message_Chebyshev(" log-both I&II&III"+mess, r, neval, a, x0, b, choice, debug);
    }
    
    return r;
}

//================================================================================================
double Integrate_Ric_Rci_JC_Chebyshev_log(int choice, double lgxl, double lgxu, double epsrel, double epsabs, double (*fptr)(double), double lgxc)
{ return Integrate_Ric_Rci_JC_Chebyshev_log(choice, "", lgxl, lgxu, epsrel, epsabs, fptr, lgxc); }

