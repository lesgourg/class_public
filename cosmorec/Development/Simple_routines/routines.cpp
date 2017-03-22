//======================================================================================
// Author: Jens Chluba 
// first implementation: Jan 2002
// last modification: June 2012
//
// Purpose: collection of several simple routines
//======================================================================================
// Jul 2012: fixed a memory problem with gsl-1.15 compatibility of splines.
// Jun 2012: fixed a memory issue with spline setup routines.
// Jan 2012: added simple routines to load tables of data and create splines
// Dec 2011: added routines for Gamma and incomplete Gamma functions

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>

#include "routines.h"
#include "Definitions.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_gamma.h>

//======================================================================================
// for xmgrace output
//======================================================================================
#ifdef GRACE_DEFINED
#include <grace_np.h>
#endif

using namespace std;

//======================================================================================
// Special functions
//======================================================================================

//======================================================================================
// Dawson integral; Based on Antia
// accuracy ~10^-15
//======================================================================================
/*  Coefficients of rational function approximations */

static double Dawson_a[9] = { 5.574810327568612837e-01,  1.480548760676749464e-01,
    2.473412335012745505e-02,  2.885568838534784757e-03,
    2.460310593071530085e-04,  1.551300911064284532e-05,
    7.087756590256507931e-07,  2.179833744510066470e-08,
    3.580284771465818339e-10};
static double Dawson_b[9] = { 1.000000000000000001e+00, -1.091856339098055081e-01,
    4.306752089643668630e-02, -1.497994816990499693e-03,
    3.339251349042084215e-04, -1.857560947421670758e-06,
    6.710614181693725933e-07,  4.960681954907270745e-09,
    2.335008164471329286e-10};

static double Dawson_a1[9] = { 6.262545365722985257e-01,  1.885254515728064716e-01,
    3.762507762042703303e-02,  4.987142830383369467e-03,
    6.094443149907482678e-04,  3.731642154042092311e-05,
    4.239542422752186968e-06,  5.982569991142279706e-08,
    1.174473438334967860e-08};
static double Dawson_b1[10] = { 1.000651029373811169e+00, -4.210501168492219455e-02,
    3.981798325257397695e-02,  1.025635538449943440e-03,
    4.099305572015143676e-04,  1.596933865686699852e-05,
    2.244388174606844525e-06,  3.090088866526039173e-08,
    5.894053468010250910e-09, -1.103344321928670385e-13};

static double Dawson_a2[8] = {-4.779664691344904517e+00,  3.275776248263644817e+00,
    -8.957008446633805629e-01,  1.364310927775567374e-01,
    -1.254364750946865384e-02,  7.203420550034988907e-04,
    -2.404381247095368252e-05,  3.964910408777002362e-07};
static double Dawson_b2[9] = {-1.813546033357404877e+00,  1.455784804742178438e+00,
    -4.178728305708361334e-01,  6.532985079094380749e-02,
    -6.100423692880221490e-03,  3.543095925141692234e-04,
    -1.192279451351772837e-05,  1.982456062768223261e-07,
    -2.989826117053882457e-16};

static double Dawson_a3[7] = {-7.277595909795730480e+01,  2.030133303927516449e+03,
    -2.623905050877405481e+04,  1.572452946972041743e+05,
    -3.959782020608236881e+05,  3.242235333091555331e+05,
    -3.695969139799184010e+04};
static double Dawson_b3[7] = { 4.999999999999999987e-01, -3.613797954897864312e+01,
    9.972476621892572722e+02, -1.263834541306106681e+04,
    7.275923849721134328e+04, -1.668382015827219241e+05,
    1.031530728542497468e+05};

double Dawson_Int(const double &x)
{
    double xa,y,fn,fd,f;
    
    xa=fabs(x);
    if(xa<2.5) {
        y=x*x;
        fn=(((((((Dawson_b[8]*y+Dawson_b[7])*y+Dawson_b[6])*y+Dawson_b[5])*y
               +Dawson_b[4])*y+Dawson_b[3])*y+Dawson_b[2])*y+Dawson_b[1])*y+Dawson_b[0];
        fd=((((((((Dawson_a[8]*y+Dawson_a[7])*y+Dawson_a[6])*y+Dawson_a[5])*y
                +Dawson_a[4])*y+Dawson_a[3])*y+Dawson_a[2])*y+Dawson_a[1])*y+Dawson_a[0])*y+1;
        f=x*fn/fd;
    }
    
    else if(xa<4.0) {
        y=x*x;
        fn=((((((((Dawson_b1[9]*y+Dawson_b1[8])*y+Dawson_b1[7])*y+Dawson_b1[6])*y
                +Dawson_b1[5])*y+Dawson_b1[4])*y+Dawson_b1[3])*y+Dawson_b1[2])*y
            +Dawson_b1[1])*y+Dawson_b1[0];
        fd=((((((((Dawson_a1[8]*y+Dawson_a1[7])*y+Dawson_a1[6])*y+Dawson_a1[5])*y
                +Dawson_a1[4])*y+Dawson_a1[3])*y+Dawson_a1[2])*y+Dawson_a1[1])*y
            +Dawson_a1[0])*y+1;
        f=x*fn/fd;
    }
    
    else if(xa<5.5) {
        y=x*x;
        fn=(((((((Dawson_b2[8]*y+Dawson_b2[7])*y+Dawson_b2[6])*y+Dawson_b2[5])*y
               +Dawson_b2[4])*y+Dawson_b2[3])*y+Dawson_b2[2])*y+Dawson_b2[1])*y+Dawson_b2[0];
        fd=(((((((Dawson_a2[7]*y+Dawson_a2[6])*y+Dawson_a2[5])*y+Dawson_a2[4])*y
               +Dawson_a2[3])*y+Dawson_a2[2])*y+Dawson_a2[1])*y+Dawson_a2[0])*y+1;
        f=x*fn/fd;
    }
    
    else {
        y=1./(x*x);
        fn=(((((Dawson_b3[6]*y+Dawson_b3[5])*y+Dawson_b3[4])*y+Dawson_b3[3])*y
             +Dawson_b3[2])*y+Dawson_b3[1])*y+Dawson_b3[0];
        fd=((((((Dawson_a3[6]*y+Dawson_a3[5])*y+Dawson_a3[4])*y+Dawson_a3[3])*y
              +Dawson_a3[2])*y+Dawson_a3[1])*y+Dawson_a3[0])*y+1;
        f=fn/(fd*x);
    }
    return f;
}

//======================================================================================
// To calculate the  Error function for a real argument; Based on Antia
// accuracy ~10^-15
//======================================================================================
/*  Coefficients of rational function approximations */
double erf_a[7] = {4.891837297874833514e-01,  1.110175596090841322e-01,
    1.526977188817169289e-02,  1.388143322498740953e-03,
    8.446154421878829637e-05,  3.239915935631695227e-06,
    6.200069065781009292e-08};
double erf_b[8] = {1.128379167095512575e+00,  1.758583405424390318e-01,
    5.411290829613803886e-02,  3.805760982281134630e-03,
    4.314620532020313106e-04,  1.423201530735308681e-05,
    6.188698151904667110e-07,  2.293915472550064153e-09};
double erf_a1[7] = {3.040844316651787557e+01,  3.358927547920980388e+02,
    1.703170048554673596e+03,  4.133195606956137156e+03,
    4.555611776312694034e+03,  1.935778559646575488e+03,
    2.051089518116697193e+02};
double erf_b1[8] = {5.641895835477562749e-01,  1.687403209467950089e+01,
    1.813522721872712655e+02,  8.779664433851135986e+02,
    1.965115658619443782e+03,  1.865558781108286245e+03,
    5.828865035607128761e+02,  2.558559157228883880e+01};

double erf_JC(const double &x)
{
    double fn,fd,f,y;
    
    if(fabs(x)<2.0) {
        /*  Use approximation for Erf(x)/x */
        y=x*x;
        fn=((((((erf_b[7]*y+erf_b[6])*y+erf_b[5])*y+erf_b[4])*y+erf_b[3])*y
             +erf_b[2])*y+erf_b[1])*y+erf_b[0];
        fd=((((((erf_a[6]*y+erf_a[5])*y+erf_a[4])*y+erf_a[3])*y+erf_a[2])*y
             +erf_a[1])*y+erf_a[0])*y+1;
        f=x*fn/fd;
    }
    
    else {
        /*  Use approximation for x*Exp(x^2)*Erfc(x) */
        y=1.0/(x*x);
        fn=((((((erf_b1[7]*y+erf_b1[6])*y+erf_b1[5])*y+erf_b1[4])*y+erf_b1[3])*y
             +erf_b1[2])*y+erf_b1[1])*y+erf_b1[0];
        fd=((((((erf_a1[6]*y+erf_a1[5])*y+erf_a1[4])*y+erf_a1[3])*y+erf_a1[2])*y
             +erf_a1[1])*y+erf_a1[0])*y+1;
        f=1.-exp(-x*x)*fn/(fd*fabs(x));
        if(x<0.0) f=-f;
    }
    return f;
}

//======================================================================================
// To calculate the complementary Error function for a real argument; Based on Antia
// accuracy ~10^-15
//======================================================================================
/*  Coefficients of rational function approximations */
static double erfc_a[7] = {4.891837297874833514e-01,  1.110175596090841322e-01,
    1.526977188817169289e-02,  1.388143322498740953e-03,
    8.446154421878829637e-05,  3.239915935631695227e-06,
    6.200069065781009292e-08};
static double erfc_b[8] = {1.128379167095512575e+00,  1.758583405424390318e-01,
    5.411290829613803886e-02,  3.805760982281134630e-03,
    4.314620532020313106e-04,  1.423201530735308681e-05,
    6.188698151904667110e-07,  2.293915472550064153e-09};
static double erfc_a1[7] = {3.040844316651787557e+01,  3.358927547920980388e+02,
    1.703170048554673596e+03,  4.133195606956137156e+03,
    4.555611776312694034e+03,  1.935778559646575488e+03,
    2.051089518116697193e+02};
static double erfc_b1[8] = {5.641895835477562749e-01,  1.687403209467950089e+01,
    1.813522721872712655e+02,  8.779664433851135986e+02,
    1.965115658619443782e+03,  1.865558781108286245e+03,
    5.828865035607128761e+02,  2.558559157228883880e+01};

double erfc_JC(const double &x)
{
    double fn,fd,f,y;
    
    if(fabs(x)<2.0) {
        /*  Use approximation for Erf(x)/x */
        y=x*x;
        fn=((((((erfc_b[7]*y+erfc_b[6])*y+erfc_b[5])*y+erfc_b[4])*y+erfc_b[3])*y
             +erfc_b[2])*y+erfc_b[1])*y+erfc_b[0];
        fd=((((((erfc_a[6]*y+erfc_a[5])*y+erfc_a[4])*y+erfc_a[3])*y+erfc_a[2])*y
             +erfc_a[1])*y+erfc_a[0])*y+1;
        f=1.-x*fn/fd;
    }
    
    else {
        /*  Use approximation for x*Exp(x^2)*Erfc(x) */
        y=1.0/(x*x);
        fn=((((((erfc_b1[7]*y+erfc_b1[6])*y+erfc_b1[5])*y+erfc_b1[4])*y+erfc_b1[3])*y
             +erfc_b1[2])*y+erfc_b1[1])*y+erfc_b1[0];
        fd=((((((erfc_a1[6]*y+erfc_a1[5])*y+erfc_a1[4])*y+erfc_a1[3])*y+erfc_a1[2])*y
             +erfc_a1[1])*y+erfc_a1[0])*y+1;
        f=exp(-x*x)*fn/(fd*fabs(x));
        if(x<0.0) f=2.-f;
    }
    return f;
}

//======================================================================================
double Gamma_JC(const double &x)                     // Gamma-function
{ return gsl_sf_gamma(min(x, 500.0)); }

//======================================================================================
double Gamma_JC(const double &x, const double &a)   // incomplete Gamma-function
{ return gsl_sf_gamma_inc(a, min(x, 500.0)); }

//======================================================================================
double max(double a, double b){ return a>b ? a : b;}
double min(double a, double b){ return a<b ? a : b;}

//======================================================================================
double dlog10factorial(double dn)
{ if(dn<=1.0) return 0.0; return log10(dn)+dlog10factorial(dn-1.0); }

double log10factorial_full(int n)
{ return dlog10factorial( (double) n); }

double dfactorial(double dn)
{ if(dn<=1.0) return 1.0; return dn*dfactorial(dn-1.0); }

double factorial(int n)
{ 
    if(n>=20) return sqrt(TWOPI*(n+1.0))*exp(-(n+1.0))*pow(n+1.0, n)*
        (1.0+8.333333333333e-2/(n+1.0)
         +3.472222222222e-3/(n+1.0)/(n+1.0)
         -2.681327160494e-3/(n+1.0)/(n+1.0)/(n+1.0)
         -2.294720936214e-4/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)
         +7.840392217201e-4/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)
         +6.972813758366e-5/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)
         -5.921664373537e-4/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0) 
         );
    
    return dfactorial( (double) n); 
}

double log10factorial(int n)
{ 
    if(n>=20) return 0.5*log10(TWOPI*(n+1.0))+n*log10(n+1.0)-(n+1.0)*0.43429448190325176
        +log10(
               1.0+8.333333333333e-2/(n+1.0)
               +3.472222222222e-3/(n+1.0)/(n+1.0)
               -2.681327160494e-3/(n+1.0)/(n+1.0)/(n+1.0)
               -2.294720936214e-4/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)
               +7.840392217201e-4/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)
               +6.972813758366e-5/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)
               -5.921664373537e-4/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0) 
               );
    
    return dlog10factorial( (double) n); 
}

double factorial_corrfac(int n)
{ 
    return      (1.0+8.333333333333e-2/(n+1.0)
                 +3.472222222222e-3/(n+1.0)/(n+1.0)
                 -2.681327160494e-3/(n+1.0)/(n+1.0)/(n+1.0)
                 -2.294720936214e-4/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)
                 +7.840392217201e-4/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)
                 +6.972813758366e-5/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)
                 -5.921664373537e-4/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0)/(n+1.0) 
                 );
}

//======================================================================================
// checking for nan
//======================================================================================
bool isnan_JC(double a)
{
#if defined(__APPLE__) 
    return !(a==a);
#else 
    return isnan(a);
#endif
}

//===================================================================================
//
// Spline interpolation using the GSL libray
//
//===================================================================================
vector<bool> allocvec;
vector<string> strvec;
vector<gsl_interp_accel *> accvec; 
vector<gsl_spline *> splinevec; 

int calc_spline_coeffies_JC(int nxi, const double *za, const double *ya, string variable)
{
    if(variable=="") strvec.push_back("variable name was not given");
    else strvec.push_back(variable);
    
    gsl_interp_accel *acc;
    accvec.push_back(acc);
    //
    gsl_spline *spline;
    splinevec.push_back(spline);
    //
    int element=splinevec.size()-1;
    
    accvec[element]=gsl_interp_accel_alloc();
    splinevec[element]=gsl_spline_alloc(gsl_interp_cspline, nxi);
    
    gsl_spline_init (splinevec[element], za, ya, nxi);
    
    allocvec.push_back(1);
    
    return element;
}

//===================================================================================
void update_spline_coeffies_JC(int memindex, int nxi, 
                               const double *za, const double *ya, 
                               string variable)
{
    if(allocvec[memindex])
    {
        if(variable!="") strvec[memindex]=variable;
        gsl_spline_free (splinevec[memindex]);
        gsl_interp_accel_free (accvec[memindex]);
        
        accvec[memindex]=gsl_interp_accel_alloc();
        splinevec[memindex]=gsl_spline_alloc(gsl_interp_cspline, nxi);  
        
        gsl_spline_init(splinevec[memindex], za, ya, nxi);
    }
    else
    { 
        cerr << " update_spline_coeffies_JC::This part of memory was not set... " << endl; 
        exit(0); 
    }

    return;
}

//===================================================================================
double calc_spline_JC(double x, int memindex, string mess)
{   
/*  
    cout.precision(16);
    cout << " x= "<< x << " memory index: " << memindex 
         << " called by: " << mess 
         << " name of variable: " << strvec[memindex] 
         << " xmin= " << splinevec[memindex]->interp->xmin 
         << " xmax= " << splinevec[memindex]->interp->xmax << endl;
*/
    // no extrapolation but this is needed to make GSL not complain if x~xmin or xmax
    if(x<=splinevec[memindex]->interp->xmin*(1.0+1.0e-14) && 
       x>=splinevec[memindex]->interp->xmin*(1.0-1.0e-14)) 
        
        x=splinevec[memindex]->interp->xmin*(1.0+1.0e-14);
    
    if(x>=splinevec[memindex]->interp->xmax*(1.0-1.0e-14) &&
       x<=splinevec[memindex]->interp->xmax*(1.0+1.0e-14)) 
        
        x=splinevec[memindex]->interp->xmax*(1.0-1.0e-14);
    
    return gsl_spline_eval(splinevec[memindex], x, accvec[memindex]); 
}

//===================================================================================
void free_spline_JC(int memindex, string mess)
{
    if( memindex < (int)allocvec.size() )
    {
        if(allocvec[memindex])
        {
            // cout << " free_spline_JC:: called by " << mess << endl;
            strvec[memindex]=" splines were deleted by " + mess;
            gsl_spline_free (splinevec[memindex]);
            gsl_interp_accel_free (accvec[memindex]);
        }
    
        allocvec[memindex]=0;
    }
    
    return;
}

void free_all_splines_JC()
{
    cout << "\n %------------------------------------------------% " << endl;
    cout << " % Clearing all GSL-splines (routines.cpp) " << endl;
    cout << " %------------------------------------------------%\n " << endl;
    
    for(int k=0; k<(int)splinevec.size(); k++)
    {
        if(allocvec[k])
        {
            gsl_spline_free (splinevec[k]);
            gsl_interp_accel_free (accvec[k]);
        }
    }
    
    splinevec.clear();
    accvec.clear();
    strvec.clear();
    allocvec.clear();
    
    return;
}

//===================================================================================
void show_spline_memory()
{
    for(int i=0; i<(int)strvec.size(); i++) 
        cout << " memory index: " << i << " name of variable: " << strvec[i] << endl;
    
    if(strvec.size()==0) cout << " show_spline_memory:: No splines allocated. " << endl;
    
    return;
}


//======================================================================================
// load tables of data (added Jan 2012)
//======================================================================================
void load_data_from_file(string fname, int cols, vector<int> &spline_mem_indices, 
                         bool logx, bool logy)
{
    cout << " load_data_from_file :: loading data from file " << fname << endl;
    cout << " logx= " << logx << " logy= " << logy << endl;
    
    ifstream ifile(fname.c_str());
    spline_mem_indices.resize(cols, -1);
    
    vector<double> vdum(cols);
    vector<vector<double> > Xi_Data;
    
    while(ifile)
    {
        for(int i=0; i<cols && ifile; i++) ifile >> vdum[i]; 
        if(ifile) Xi_Data.push_back(vdum);
    }
    
    cout << " load_data_from_file :: setting up splines " << endl;
    
    vector<double> xarr(Xi_Data.size());
    vector<double> yarr(Xi_Data.size());
    
    if(Xi_Data[0][0]<Xi_Data[1][0])
    {       
        
        if(logx) for(int i=0; i< (int)Xi_Data.size(); i++) xarr[i] = log(Xi_Data[i][0]);
        else for(int i=0; i< (int)Xi_Data.size(); i++) xarr[i] = Xi_Data[i][0];
        
        for(int c=1; c<cols; c++)
        {
            if(logy) for(int i=0; i<(int)Xi_Data.size(); i++) yarr[i] = log(Xi_Data[i][c]);
            else for(int i=0; i<(int)Xi_Data.size(); i++) yarr[i] = Xi_Data[i][c];
            
            spline_mem_indices[c-1]=calc_spline_coeffies_JC((int)Xi_Data.size(), 
                                                            &xarr[0], &yarr[0], 
                                                            " load_data_from_file ");
        }
    }
    else 
    {
        int nmax=Xi_Data.size()-1;
        
        if(logx) for(int i=0; i<(int)Xi_Data.size(); i++) xarr[nmax-i] = log(Xi_Data[i][0]);
        else for(int i=0; i< (int)Xi_Data.size(); i++) xarr[nmax-i] = Xi_Data[i][0];
        
        for(int c=1; c<cols; c++)
        {
            if(logy) for(int i=0; i<(int)Xi_Data.size(); i++) yarr[nmax-i] = log(Xi_Data[i][c]);
            else for(int i=0; i<(int)Xi_Data.size(); i++) yarr[nmax-i] = Xi_Data[i][c];
            
            spline_mem_indices[c-1]=calc_spline_coeffies_JC((int)Xi_Data.size(), 
                                                            &xarr[0], &yarr[0], 
                                                            " load_data_from_file ");
        }        
    }
    
    vdum.clear();
    Xi_Data.clear();
    xarr.clear();
    yarr.clear();
    
    return;
}

void load_data_from_file_loglin(string fname, int cols, vector<int> &spline_mem_indices)
{
    load_data_from_file(fname, cols, spline_mem_indices, 1, 0);
    return;
}

void load_data_from_file_loglog(string fname, int cols, vector<int> &spline_mem_indices)
{
    load_data_from_file(fname, cols, spline_mem_indices, 1, 1);
    return;
}

//======================================================================================
// npol-1 is degree of the interpolating polynomial
//======================================================================================
void polint_routines(const double *xa, const double *ya, int n, const double x, double &y, double &dy)
{
    int i,m,ns=0;
    double den,dif,dift,ho,hp,w;
    
    vector<double> c(n),d(n);
    dif=fabs(x-xa[0]);
    for (i=0;i<n;i++) {
        if ((dift=fabs(x-xa[i])) < dif) {
            ns=i;
            dif=dift;
        }
        c[i]=ya[i];
        d[i]=ya[i];
    }
    y=ya[ns--];
    for (m=1;m<n;m++) {
        for (i=0;i<n-m;i++) {
            ho=xa[i]-x;
            hp=xa[i+m]-x;
            w=c[i+1]-d[i];
            if ((den=ho-hp) == 0.0)
            { 
                cout << " polint error " << ho << " " << x << " " << hp << " " 
                     << i << " " << xa[i] << " " << ya[i] << " " 
                     << m << " " << xa[i+m] << " " << ya[i+m] 
                     << endl; 
                
                wait_f_r(); 
            }
            den=w/den;
            d[i]=hp*den;
            c[i]=ho*den;
        }
        y += (dy=(2*(ns+1) < (n-m) ? c[ns+1] : d[ns--]));
    }
}

//======================================================================================
void polint_JC(const double *xa, const double *ya, int na, const double x, 
               int &istart, int npol, double *y, double *dy)
{
    // npol-1 is degree of the interpolating polynomial
    long unsigned int j=istart;
    locate_JC(xa, na-1, x, &j);      // find index corresponding to x (start at i=istart)
    istart=j;
    
    if(istart==na-1){ *y=0.0; *dy=0.0; return; } // no point above x[j]! 

    int npol_loc=(int)min(na, npol); // can only use as many as na points for interpolation
    int nlow=(npol_loc-1)/2; //, nup=npol_loc-1-nlow; // nlow+nup+1 == npol_loc!!!
    int ks=(int)min(max(istart-nlow, 0), na-1-npol_loc);  // defining region arround index

    polint_routines(&xa[ks], &ya[ks], npol_loc, x, *y, *dy);
        
    return;
}

//======================================================================================
void polint_JC(const double *xa, const double *ya, int na, const double x, 
               int npol, double *y, double *dy)
{
    // npol-1 is degree of the interpolating polynomial
    int j=0;
    polint_JC(xa, ya, na, x, j, npol, y, dy);
    return;
}

//======================================================================================
// simple grid setup
//======================================================================================
void init_xarr(double x0, double xm, double *xarr, int npts, int method_flag, int mess_flg)
{
    if(method_flag!=0 && method_flag!=1 && method_flag!=2) 
    {
    cout << " Choose initialisation strategy for x-array" 
         << " (0 for linear or 1 for log or 2 for log-linear) !" << endl;
    cin  >> method_flag;
    }
    
    double dx=0.0;
    int i;
    
    // linear
    if(method_flag==0)
    {
      if(mess_flg==1) cout << "\n Initializing linear in x with " << npts 
                           << " points\n" << endl;
    dx=(xm-x0)/(npts-1);
    // filling array    
    for(xarr[0]=x0, i=1; i<npts; i++) xarr[i]=xarr[i-1]+dx;
    }
    
    // log     
    if(method_flag==1)
    {
    if(mess_flg==1) cout << "\n Initializing logarithmic in x with " << npts 
                         << " points\n" << endl;
    dx=pow(xm/x0, 1.0/(double)(npts-1));
    // filling array
    for(xarr[0]=x0, i=1; i<npts; i++) xarr[i]=xarr[i-1]*dx;
    }
    
    // log (1/4) - linear (3/4)   
    if(method_flag==2)
    {
    if(mess_flg==1) cout << "\n Initializing first 1/30 logarithmic the last 29/30 linear in x with " 
                 << npts << " points\n" << endl;
    init_xarr(x0, xm/30.0, xarr, npts/4, 1, 0);
    init_xarr(xm/30.0, xm, &xarr[npts/4-1], npts-(npts/4-1), 0, 0);
    }
    
    return;
}  

void init_xarr(double x0, double xm, double *xarr, int npts, int method_flag)
{ init_xarr(x0, xm, xarr, npts, method_flag, 1); return; }

//======================================================================================
void wait_f_r()
{
    char c;
    cout << " To go on press return! " << endl;
    cin.get(c);
    return;
}

void wait_f_r(string mess)
{
    char c;
    cout << mess << endl;
    cin.get(c);
    return;
}

void wait_f_r(int mess)
{
    char c;
    cout << mess << endl;
    cin.get(c);
    return;
}

void wait_f_r(double mess)
{
    char c;
    cout << mess << endl;
    cin.get(c);
    return;
}

//======================================================================================
void locate_JC(const double xx[], unsigned long n, double x, unsigned long *j)
{
    // returns index j with xx[j]<= x; outside the boundaries it 
    // will return the upper or lower index.
    unsigned long ju,jm,jl;
    int ascnd=(xx[n-1] > xx[0]);
    
    if(ascnd==1){
        if(x <= xx[0]){ *j=0; return; }
        if(x >= xx[n-1]){ *j=n-1; return; }
    }
    else{
        if(x >= xx[0]){ *j=0; return; }
        if(x <= xx[n-1]){ *j=n-1; return; } 
    }
    
    jl=0;
    ju=n;
    while (ju-jl > 1) {
        jm=(ju+jl) >> 1;
        if ( (x >= xx[jm]) == ascnd)
            jl=jm;
        else
            ju=jm;
    }
    *j=jl;
    
    return;
}

void hunt(const double xx[], unsigned long n, double x, unsigned long *jlo)
{
    unsigned long jm,jhi,inc;
    int ascnd;

    ascnd=(xx[n] > xx[1]);
    if (*jlo <= 0 || *jlo > n) {
        *jlo=0;
        jhi=n+1;
    } else {
        inc=1;
        if ( (x >= xx[*jlo]) == ascnd) {
            if (*jlo == n) return;
            jhi=(*jlo)+1;
            while ( (x >= xx[jhi]) == ascnd) {
                *jlo=jhi;
                inc += inc;
                jhi=(*jlo)+inc;
                if (jhi > n) {
                    jhi=n+1;
                    break;
                }
            }
        } else {
            if (*jlo == 1) {
                *jlo=0;
                return;
            }
            jhi=(*jlo)--;
            while ( (x < xx[*jlo]) == ascnd) {
                jhi=(*jlo);
                inc <<= 1;
                if (inc >= jhi) {
                    *jlo=0;
                    break;
                }
                else *jlo=jhi-inc;
            }
        }
    }
    if (x == xx[1]) {*jlo=1;  return;}
    while (jhi-(*jlo) != 1) {
        jm=(jhi+(*jlo)) >> 1;
        if ( (x > xx[jm]) == ascnd)
            *jlo=jm;
        else
            jhi=jm;
    }
}

//======================================================================================
// i/o modules
//======================================================================================
string int_to_string(int i)
{
    if(i==0) return "0";
    
    char   buf[16];
    int dig=10, num=0;
    string name="";
    
    if(i<0)
    {
        name="-";
        i=-i;
    }
    
    while(i>=dig) dig*=10;
    
    while(dig!=1)
    {
        dig/=10;
        num = (int)(i/dig);
        sprintf(buf, "%d", num); 
        name += *buf;
        i-=num*dig;
    }
    
    return name;
}

string int_to_string(int i, int ni)
{
    string name=int_to_string(i);
    int len=name.length();

    for(int l=len; l<ni; l++) name="0"+name;
    
    return name;
}

//======================================================================================
// root-finding methods
//======================================================================================
inline double SIGN(const double &a, const double &b)
{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

// bisection
// xacc is relative accuracy
double find_root(double (* func)(double *), double x1, double x2, double xacc)
{
    const int find_root_JMAX=40;
    int j;
    double dx,f,fmid,xmid,rtb;
    
    f=func(&x1);
    fmid=func(&x2);
    if (f*fmid >= 0.0)
    { cerr << " Root must be bracketed for bisection" << endl; return x2;}
    
    rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
    
    for (j=0;j<find_root_JMAX;j++) {
        dx*= 0.5;
        xmid=rtb+dx;
        fmid=func(&xmid);
        if (fmid <= 0.0) rtb=xmid;
        if (fabs(dx/xmid) < xacc || fmid == 0.0) return rtb;
    }
    
    cerr << "Too many bisections in find_root" << endl;
    return x2;
}

// Brent's method
// xacc is relative accuracy
double find_root_brent(double (* func)(double *), double x1, double x2, double xacc)
{
    const int find_root_brent_ITMAX=100;
    const double find_root_brent_EPS=1.0e-10;
    
    int iter;
    double a=x1,b=x2,c=x2,d=0.0,e=0.0 ,min1,min2;
    double fa=func(&a),fb=func(&b),fc,p,q,r,s,tol1,xm;
    
    if (fa*fb >= 0.0)
    { cerr << " Root must be bracketed for bisection" << endl; return x2; }
    
    fc=fb;
    for (iter=0;iter<find_root_brent_ITMAX;iter++) {
        if (fb*fc>=0.0) {
            c=a;
            fc=fa;
            e=d=b-a;
        }
        
        if (fabs(fc) < fabs(fb)) {
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        }
        
        tol1=2.0*find_root_brent_EPS*fabs(b)+fabs(xacc*b);
        xm=0.5*(c-b);
        
        if (fabs(xm) <= tol1 || fb == 0.0) return b;
        
        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
            s=fb/fa;
            if (a == c) {
                p=2.0*xm*s;
                q=1.0-s;
            } else {
                q=fa/fc;
                r=fb/fc;
                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                q=(q-1.0)*(r-1.0)*(s-1.0);
            }
            if (p > 0.0) q = -q;
            p=fabs(p);
            min1=3.0*xm*q-fabs(tol1*q);
            min2=fabs(e*q);
            if (2.0*p < (min1 < min2 ? min1 : min2)) {
                e=d;
                d=p/q;
            } else {
                d=xm;
                e=d;
            }
        } else {
            d=xm;
            e=d;
        }
        a=b;
        fa=fb;
        if (fabs(d) > tol1)
            b += d;
        else
            b += SIGN(tol1,xm);
        fb=func(&b);
    }
    
    cerr << "Maximum number of iterations exceeded in find_root_brent" << endl;
    return x2;
}

// Brent's method
// xacc is relative accuracy
double find_root_brent(double (* func)(double *, void *p), void *para, 
                       double x1, double x2, double xacc)
{
    const int find_root_brent_ITMAX=100;
    const double find_root_brent_EPS=1.0e-10;
    
    int iter;
    double a=x1,b=x2,c=x2,d=0.0,e=0.0 ,min1,min2;
    double fa=func(&a, para),fb=func(&b, para),fc,p,q,r,s,tol1,xm;
    
    if (fa*fb >= 0.0)
    { cerr << " Root must be bracketed for bisection" << endl; return x2; }
    
    fc=fb;
    for (iter=0;iter<find_root_brent_ITMAX;iter++) {
        if (fb*fc>=0.0) {
            c=a;
            fc=fa;
            e=d=b-a;
        }
        
        if (fabs(fc) < fabs(fb)) {
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        }
        
        tol1=2.0*find_root_brent_EPS*fabs(b)+fabs(xacc*b);
        xm=0.5*(c-b);
        
        if (fabs(xm) <= tol1 || fb == 0.0) return b;
        
        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
            s=fb/fa;
            if (a == c) {
                p=2.0*xm*s;
                q=1.0-s;
            } else {
                q=fa/fc;
                r=fb/fc;
                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                q=(q-1.0)*(r-1.0)*(s-1.0);
            }
            if (p > 0.0) q = -q;
            p=fabs(p);
            min1=3.0*xm*q-fabs(tol1*q);
            min2=fabs(e*q);
            if (2.0*p < (min1 < min2 ? min1 : min2)) {
                e=d;
                d=p/q;
            } else {
                d=xm;
                e=d;
            }
        } else {
            d=xm;
            e=d;
        }
        a=b;
        fa=fb;
        if (fabs(d) > tol1)
            b += d;
        else
            b += SIGN(tol1,xm);
        fb=func(&b, para);
    }
    
    cerr << "Maximum number of iterations exceeded in find_root_brent" << endl;
    return x2;
}

//======================================================================================
// for xmgrace output; added April 2011; does not work very well yet...
//======================================================================================
#ifdef GRACE_DEFINED
void plot_xy_function(const vector<double> &xarr, const vector<double> &yarr, int flg=0)
{
    if (GraceOpen(1028) == -1) 
	{
		cerr << " can't run grace " << endl; 
		exit(-1);
	}
	
	//===============================================================
	// Send some initialization commands to Grace 
	//===============================================================
	GracePrintf("s0 COLOR 1");
	GracePrintf("s0 LINESTYLE 1");
	GracePrintf("s0 on");
    
    if(flg==2) GracePrintf("XAXES SCALE LOGARITHMIC");
    else GracePrintf("XAXES SCALE NORMAL");
	if(flg==0) GracePrintf("YAXES SCALE NORMAL");
    else GracePrintf("YAXES SCALE LOGARITHMIC");
	
	//===============================================================
	// plot set
	//===============================================================
    for(int i=0; i<(int)xarr.size(); i++) 
    {
        GracePrintf("g0.s0 point %10.16f, %10.16f", xarr[i], yarr[i]);
    }
    
    GracePrintf("AUTOTICKS");
    GracePrintf("AUTOSCALE s0");
	GracePrintf("PRINT");	
	GracePrintf("redraw");
    
	//===============================================================
	// wait before cleaning up
	//===============================================================
    wait_f_r();
    GraceClose();
    
	return;
}
#else
void plot_xy_function(const vector<double> &xarr, const vector<double> &yarr, int flg=0)
{
    cerr << " plot_xy_function::Grace support was not activated. Exiting. " << endl;
    exit(0);
    return;
}
#endif

void plot_xy_function_linear(const vector<double> &xarr, const vector<double> &yarr)
{ plot_xy_function(xarr, yarr, 0); return; }

void plot_xy_function_log(const vector<double> &xarr, const vector<double> &yarr)
{ plot_xy_function(xarr, yarr, 1); return; }

void plot_xy_function_loglog(const vector<double> &xarr, const vector<double> &yarr)
{ plot_xy_function(xarr, yarr, 2); return; }

//======================================================================================
// Wigner 3J symbol (added 17.05.2011)
//======================================================================================
double Wigner_3J_Symbol(int j1, int j2, int j3, int m1, int m2, int m3)
{ return gsl_sf_coupling_3j (2*j1, 2*j2, 2*j3, 2*m1, 2*m2, 2*m3); }

//======================================================================================
// divide and conquer sum
//======================================================================================
double DC_sum(const double *yarr, int M)
{
    if(M<=4) 
    {
        double r=yarr[0];
        for(int i=1; i<M; i++) r+=yarr[i];
        return r;
    }

    int N=floor(M/2);
    return DC_sum(&yarr[0], N)+DC_sum(&yarr[N], M-N);
}

//======================================================================================
double DC_sum(const vector<double> &yarr)
{
    int M=yarr.size();
    
    if(M<=4) 
    {
        double r=yarr[0];
        for(int i=1; i<M; i++) r+=yarr[i];
        return r;
    }

    int N=floor(M/2);
    return DC_sum(&yarr[0], N)+DC_sum(&yarr[N], M-N);
}

//======================================================================================
// divide and conquer sum of product
//======================================================================================
double DC_sumprod(const double *yarr, const double *zarr, int M)
{
    if(M<=4) 
    {
        double r=yarr[0]*zarr[0];
        for(int i=1; i<M; i++) r+=yarr[i]*zarr[i];
        return r;
    }

    int N=floor(M/2);
    return DC_sumprod(&yarr[0], &zarr[0], N)+DC_sumprod(&yarr[N], &zarr[N], M-N);
}

//======================================================================================
double DC_sumprod(const vector<double> &yarr, const vector<double> &zarr)
{
    int M=yarr.size();
    
    if(M<=4) 
    {
        double r=yarr[0]*zarr[0];
        for(int i=1; i<M; i++) r+=yarr[i]*zarr[i];
        return r;
    }

    int N=floor(M/2);
    return DC_sumprod(&yarr[0], &zarr[0], N)+DC_sumprod(&yarr[N], &zarr[N], M-N);
}

//======================================================================================
//======================================================================================
