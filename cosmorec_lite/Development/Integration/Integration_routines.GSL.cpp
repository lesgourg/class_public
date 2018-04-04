//==================================================================================================
// Author Jens Chluba Jan 2010
//==================================================================================================
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>

#include "routines.h"
#include "Integration_routines.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

using namespace std;

//==================================================================================================
// routines for integration from GSL
//==================================================================================================
double Integrate_gk15_GSL(double a, double b, double epsrel, double epsabs,
                          double (*fptr)(double, void *), void *para)
{
    if(a>=b) return 0.0;
    
    int neval=5000;
    gsl_integration_workspace *w=gsl_integration_workspace_alloc(neval);
    gsl_function F;
    F.function = fptr;
    F.params=para;

    double r=0.0;
    double epsest;

    gsl_integration_qag (&F, a, b, epsabs, epsrel, neval, 1, w, &r, &epsest); 
    gsl_integration_workspace_free(w);
        
    return r;
}

double Integrate_gk31_GSL(double a, double b, double epsrel, double epsabs, 
                          double (*fptr)(double, void *), void *para)
{
    if(a>=b) return 0.0;
    
    int neval=5000;
    gsl_integration_workspace *w=gsl_integration_workspace_alloc(neval);
    gsl_function F;
    F.function = fptr;
    F.params=para;
    
    double r=0.0;
    double epsest;
    
    gsl_integration_qag (&F, a, b, epsabs, epsrel, neval, 3, w, &r, &epsest); 
    gsl_integration_workspace_free(w);
    
    return r;
}

double Integrate_gk61_GSL(double a, double b, double epsrel, double epsabs, 
                          double (*fptr)(double, void *), void *para)
{
    if(a>=b) return 0.0;

    int neval=5000;
    gsl_integration_workspace *w=gsl_integration_workspace_alloc(neval);
    gsl_function F;
    F.function = fptr;
    F.params=para;
    
    double r=0.0;
    double epsest;
    
    gsl_integration_qag (&F, a, b, epsabs, epsrel, neval, 6, w, &r, &epsest); 
    gsl_integration_workspace_free(w);
    
    return r;
}

double Integrate_gk_GSL(int key, double a, double b, double epsrel, double epsabs, 
                        double (*fptr)(double, void *), void *para)
{
    if(a>=b) return 0.0;
    
    int neval=5000;
    gsl_integration_workspace *w=gsl_integration_workspace_alloc(neval);
    gsl_function F;
    F.function = fptr;
    F.params=para;
    
    double r=0.0;
    double epsest;
    key=min(6, key);
    gsl_integration_qag (&F, a, b, epsabs, epsrel, neval, key, w, &r, &epsest); 
    gsl_integration_workspace_free(w);
    
    return r;
}

//==================================================================================================
//==================================================================================================

