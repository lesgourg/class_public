//==========================================================================================
// Author: Jens Chluba 
// last modified: Oct 2010
//
// Purpose: computation of the Voigt-profile and integrals over it
//==========================================================================================
// note: x  = (nu-nu0)/ DnuT
//       a == Voigt parameter
//==========================================================================================

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cmath>

#include "Voigtprofiles.h"
#include "physical_consts.h"
#include "routines.h"

using namespace std;

bool Voigtprofiles_verbose=0;

//======================================================================================
// Class Voigtprofile_Base
//======================================================================================
Voigtprofile_Base :: Voigtprofile_Base()
{
    //==================================================================================
    // values for HI Lyman-alpha line as default
    //==================================================================================
    nu21=2.46604e+15;                                  // Hz
    lambda21=1.21568e-5;                               // cm
    A21=6.2649e+8;                                     // 1/s
    f12=0.416423;                                      // absorption oscillator strenght
    
    Gamma=A21;
    AM=1.0;
    
    k_mHc2=const_kB/const_mH_gr/const_cl/const_cl;
    
    if(Voigtprofiles_verbose) Voigtprofile_Base :: display_atomic_data();
}

Voigtprofile_Base :: ~Voigtprofile_Base(){}

void Voigtprofile_Base :: Set_atomic_data(double nu21v, double lambda21v, 
                                          double A21v, double f12v, double AMv)
{
    nu21=nu21v;
    lambda21=lambda21v;
    A21=A21v;
    f12=f12v;
    Gamma=A21;
    AM=AMv;
    
    if(Voigtprofiles_verbose)
    {
        cout << " Atomic data was reset " << endl; 
        Voigtprofile_Base :: display_atomic_data();
    }
    
    return;
}

void Voigtprofile_Base :: Set_atomic_data_Gamma(double nu21v, double lambda21v, 
                                                double A21v, double f12v, 
                                                double Gv, double AMv)
{
    nu21=nu21v;
    lambda21=lambda21v;
    A21=A21v;
    f12=f12v;
    Gamma=Gv;
    AM=AMv;
    
    if(Voigtprofiles_verbose)
    {
        cout << " Atomic data was reset " << endl; 
        Voigtprofile_Base :: display_atomic_data();
    }
    
    return;
}

void Voigtprofile_Base :: display_atomic_data() const
{
  cout << "\n %==============================================================%" << endl;
  cout << " % Voigtprofile_Base::display_atomic_data: " << endl;
  cout << " % " << nu21 << " " << lambda21 << endl;
  cout << " % " << A21 << " f12: " << f12 << " Gamma: " << Gamma << endl;
  cout << " % Atomic Mass in units of mp: " << AM << endl;
  cout << " %==============================================================%" << endl << endl;
  
  return;
}

//======================================================================================
// Voigt profile
//======================================================================================
double Voigtprofile_Base :: phi_DC(double x)
{ return max(exp(-x*x)/SQRT_PI, 1.0e-100); }

double Voigtprofile_Base :: phi_DC(double x, double a)
{ return phi_DC(x)*exp(a*a)*erfc_JC(a); }

double Voigtprofile_Base :: phi_wings(double x, double a)
{ 
  register double x2=x*x, d=a*a/x2; 
  register double K0=1.0+(1.5+(15.0/4.0+105.0/8.0/x2)/x2)/x2;
  register double K1=1.0+(5.0+105.0/4.0/x2)/x2;
  return a/(x2*PI)*(K0-d*K1); 
}

//======================================================================================
// 1. derivative of Voigt profile
//======================================================================================
double Voigtprofile_Base :: dphi_dx_DC(double x)
{ return max(-2.0*x*exp(-x*x)/SQRT_PI, 1.0e-100); }

double Voigtprofile_Base :: dphi_dx_DC(double x, double a)
{ return dphi_dx_DC(x)*exp(a*a)*erfc_JC(a); }

double Voigtprofile_Base :: dphi_dx_wings(double x, double a)
{ 
  register double x2=x*x, a2=a*a; 
  return -2.0*a/(x*x2*PI)*(1.0 + (2.0*(1.5-a2)+(15.0*(0.75-a2)+105.0*(0.5-a2)/x2)/x2)/x2 ); 
}

//======================================================================================
// Integral over wing region Int_-infinty^x phi_wings dx
//======================================================================================
double Voigtprofile_Base :: xi_Int_wings(double x, double a)
{ 
  register double x2=x*x, a2=a*a; 
  return -a/PI/x*(1.0 +((0.5-a2/3.0) + ((0.75-a2) + (1.875-15.0*a2/4.0)/x2)/x2)/x2 );
}



//======================================================================================
// Class Voigtprofile_Dawson
//======================================================================================
Voigtprofile_Dawson :: Voigtprofile_Dawson():Voigtprofile_Base(){}
Voigtprofile_Dawson :: ~Voigtprofile_Dawson(){}

//======================================================================================
// for fabs(x)> xwing use the wing approximation 
//======================================================================================
const double Voigtprofile_Dawson_x_wing=30.0;

//======================================================================================
// Dawson-integral exp(-x^2) * int_0^x  exp(y^2) dy
//======================================================================================
double Voigtprofile_Dawson :: Dawson(double x)
{ return Dawson_Int(x); }

void Voigtprofile_Dawson :: eval_Hn(double x, double *H)
{
  double x2=x*x;

  double F=Dawson(x);
  H[0]=exp(-x2);
  H[1]=2.0/SQRT_PI*(2.0*x*F-1.0);
  H[2]=H[0]*(1.0-2.0*x2);
  H[3]=4.0/3.0/SQRT_PI*( x2-1.0 + F*x*(3.0-2.0*x2) );
  H[4]=H[0]/6.0*(3.0+4.0*x2*(x2-3.0));
  H[5]=2.0/15.0/SQRT_PI*( -4.0+(9.0-2.0*x2)*x2+F*x*(15.0-4.0*x2*(5.0-x2)) );
  H[6]=H[0]/90.0*(15.0-2.0*(45.0-2.0*(15.0-2.0*x2)*x2)*x2 );
  H[7]=-2.0*(24.0-(87.0-(40.0-4.0*x2)*x2)*x2-x*(105.0-(210.0-(84.0-8.0*x2)*x2)*x2)*F)/315.0/SQRT_PI;
  H[8]=H[0]*(105.0 - (840.0-(840.0-(224.0-16.0*x2)*x2)*x2)*x2)/2520.0;
    
  return;
}

//======================================================================================
// Voigt-profile based on expansion in a (Mihalas) 
// for a<~ 0.01 these should be accurate <~ 1.0e-8
// the functions were optimized for a~0.001. There is a peak close 
// the center, but otherwise things go fine.
//======================================================================================
double Voigtprofile_Dawson :: phi(double x, double a)
{
  if(fabs(x)>=Voigtprofile_Dawson_x_wing) return phi_wings(x, a);
  double r=0.0, H[9];
  eval_Hn(x, H);
  r=H[0]+(H[1]+(H[2]+(H[3]+(H[4]+(H[5]+(H[6]+(H[7]+H[8]*a)*a)*a)*a)*a)*a)*a)*a;

  return r/SQRT_PI;
}

void Voigtprofile_Dawson :: eval_dHn(double x, double *dH)
{
  double x2=x*x;

  double F=Dawson(x);
  dH[0]=-2.0*x*exp(-x2);
  dH[1]= 4.0/SQRT_PI*(F+x-2.0*F*x2);
  dH[2]=dH[0]*(3.0-2.0*x2);
  dH[3]=4.0/3.0/SQRT_PI*( (5.0-2.0*x2)*x + F*(3.0+4.0*(x2-3.0)*x2) );
  dH[4]=dH[0]/6.0*(15.0-4.0*(5.0-x2)*x2);
  dH[5]=2.0/15.0/SQRT_PI*( (33.0+4.0*(x2-7.0)*x2)*x + F*( 15.0-(90.0-(60.0-8.0*x2)*x2)*x2) );
  dH[6]=dH[0]/90.0*(105.0-(210.0-(84.0-8.0*x2)*x2)*x2);
  dH[7]=2.0/315.0/SQRT_PI*( (279.0-(370.0-(108.0-8.0*x2)*x2)*x2)*x 
                        + F*(105.0-(840.0-(840.0-(224.0-16.0*x2)*x2)*x2)*x2) );
  dH[8]=dH[0]/315.0*(945.0/8.0-(315.0-(189.0-(36.0-2.0*x2)*x2)*x2)*x2);

  return;
}

//======================================================================================
// 1. derivative of Voigt-profile based on expansion in a (Mihalas) 
// for a<~ 0.01 these should be accurate <~ 1.0e-8
// the functions were optimized for a~0.001. There is a peak close 
// the center, but otherwise things go fine.
//======================================================================================
double Voigtprofile_Dawson :: dphi_dx(double x, double a)
{
  if(fabs(x)>=Voigtprofile_Dawson_x_wing) return dphi_dx_wings(x, a);
  double r=0.0, dH[9];
  eval_dHn(x, dH);
  r=dH[0]+(dH[1]+(dH[2]+(dH[3]+(dH[4]+(dH[5]+(dH[6]+(dH[7]+dH[8]*a)*a)*a)*a)*a)*a)*a)*a;

  return r/SQRT_PI;
}

//======================================================================================
// approximation for int_{-xmin}^xmax phi(y) dy
//======================================================================================
void Voigtprofile_Dawson :: eval_IHn(double t1, double t2, double *IH)
{
    double t12=t1*t1, t22=t2*t2;
    double F1=Dawson(t1), F2=Dawson(t2);
    double Exp1=exp(-t12), Exp2=exp(-t22);
    double E1=erf_JC(t1), E2=erf_JC(t2);
    
    IH[0]=(E2-E1)*SQRT_PI*0.5;
    IH[1]=2.0/SQRT_PI*(F1-F2);
    IH[2]=t2*Exp2-t1*Exp1;
    IH[3]=-2.0*( t2 - t1 + (1.0 - 2.0*t22)*F2 - (1.0 - 2.0*t12)*F1 )/3.0/SQRT_PI;
    IH[4]=( (3.0 - 2.0*t22)*t2*Exp2 - (3.0 - 2.0*t12)*t1*Exp1 )/6.0;
    IH[5]=( (5.0 - 2.0*t12)*t1 - (5.0 - 2.0*t22)*t2 
          + (3.0 + 4.0*t12*(-3.0 + t12))*F1 - (3.0 + 4.0*t22*(-3.0 + t22))*F2 )/15.0/SQRT_PI;
    IH[6]=( t2*(15.0 + 4*t22*(-5.0 + t22))*Exp2 - t1*(15.0 + 4*t12*(-5.0 + t12))*Exp1 )/90.0;
    IH[7]=( (33.0 - (28.0-4.0*t12)*t12)*t1
           -(33.0 - (28.0-4.0*t22)*t22)*t2 
           +(15.0 - (90.0-(60.0-8.0*t12)*t12)*t12)*F1
           -(15.0 - (90.0-(60.0-8.0*t22)*t22)*t22)*F2 )/315.0/SQRT_PI;
    IH[8]=( t2*( 105.0 - (210.0-(84.0-8.0*t22)*t22)*t22)*Exp2
           -t1*( 105.0 - (210.0-(84.0-8.0*t12)*t12)*t12)*Exp1)/2520.0;
    
    return;
}

//======================================================================================
double Voigtprofile_Dawson :: xi_Int_Dawson(double xmin, double xmax, double a)
{
    double r=0.0, IH[9];
    eval_IHn(xmin, xmax, IH);
    r=IH[0]+(IH[1]+(IH[2]+(IH[3]+(IH[4]+(IH[5]+(IH[6]+(IH[7]+IH[8]*a)*a)*a)*a)*a)*a)*a)*a;
    
    return r/SQRT_PI;
}

//======================================================================================
// int_{-xmin}^xmax phi(y) dy
//======================================================================================
double Voigtprofile_Dawson :: xi_Int(double xmin, double xmax, double a)
{
    if(xmin==xmax) return 0.0;

    double r=0.0, fac=1.0;
    // change order of integration
    if(xmin>xmax){ double dum=xmin; xmin=xmax; xmax=dum; fac=-1.0; }
    //
    double ab, bb;
    double r1=0.0;
    
    if(xmin<=0 && xmax<=0)
    { 
        // part from the wing region
        ab=min(xmin,-Voigtprofile_Dawson_x_wing); bb=min(xmax,-Voigtprofile_Dawson_x_wing); r1=0.0;   
        if(ab<bb) r1=xi_Int_wings(bb, a)-xi_Int_wings(ab, a);
        r+=r1;
    
        ab=max(xmin, -Voigtprofile_Dawson_x_wing); bb=max(xmax, -Voigtprofile_Dawson_x_wing); r1=0.0;   
        if(ab<bb) r1=xi_Int_Dawson(ab, bb, a);
        r+=r1;
    }
    else if(xmin>0 && xmax>0)
    { 
        // part from the wing region
        ab=min(-xmax,-Voigtprofile_Dawson_x_wing); bb=min(-xmin,-Voigtprofile_Dawson_x_wing); r1=0.0;   
        if(ab<bb) r1=xi_Int_wings(bb, a)-xi_Int_wings(ab, a);
        r+=r1;
        
        ab=max(-xmax, -Voigtprofile_Dawson_x_wing); bb=max(-xmin, -Voigtprofile_Dawson_x_wing); r1=0.0;       
        if(ab<bb) r1=xi_Int_Dawson(ab, bb, a);
        r+=r1;
    }
    else if(xmin<=0 && xmax>0)
    { 
        r1=0.0; 
        r1=1.0-xi_Int(xmax, xmax*1.0e+6, a);
        r+=r1;
        
        r1=0.0; 
        r1=-xi_Int(-1.0e+6, xmin, a);       
        r+=r1;
    }
    else cout << " Voigtprofile_Dawson :: xi_Int: this case was not expected..." << endl; 
    
    return r*fac;
}

//======================================================================================
// int_{-infinity}^x phi(y) dy
//======================================================================================
double Voigtprofile_Dawson :: xi_Int(double x, double a)
{ return xi_Int(-1.0e+4, x, a); }

//======================================================================================
// int_{-xc}^x phi(y) dy, with xc = -nu0/DnuD
//======================================================================================
double Voigtprofile_Dawson :: xi_Int_xc(double x, double a)
{ return xi_Int(-1.0/a2Dn_n0(a), x, a); }

//======================================================================================
// int_x^{infinity} phi(y) dy
//======================================================================================
double Voigtprofile_Dawson :: chi_Int(double x, double a)
{ return xi_Int(x, 1.0e+5, a); }

