//==================================================================================================
//
// Author Jens Chluba July 2010
// last modification:  Aug 2012
//
//==================================================================================================
// 27.06.2011: added PDE setup for integro-differential equation systems
// 24.02.2011: added PDE solver with Integral part as correction 
// 01.02.2011: added possibility to take integral part over +/-2 off-diags into account

#ifndef PDE_SOLVER_H
#define PDE_SOLVER_H

#include <vector>

using namespace std;

//==================================================================================================
struct Lagrange_interpolation_coefficients
{
    vector<vector<double> > dli_dx;
    vector<vector<double> > d2li_d2x;
    vector<double> Int_li;
};

//==================================================================================================
struct PDE_Stepper_Data
{
    vector<double> *ptr;    
    //
    vector<double> Ai;
    vector<double> Bi;
    vector<double> Ci;
    vector<double> Di;
    vector<double> Ii;
    // used as workspace
    vector<double> Ui;
    vector<double> Vi;
    vector<double> bi;
    vector<double> zi; 
    //
    vector<double> Aip;
    vector<double> Bip;
    vector<double> Cip;
    vector<double> Dip;
    vector<double> Iip;
    //
    vector<double> *Ai_ptr; 
    vector<double> *Bi_ptr; 
    vector<double> *Ci_ptr; 
    vector<double> *Di_ptr; 
    vector<double> *Ii_ptr; 
    //
    vector<double> *Aip_ptr;    
    vector<double> *Bip_ptr;    
    vector<double> *Cip_ptr;    
    vector<double> *Dip_ptr;    
    vector<double> *Iip_ptr;    
    //
    vector<vector<double> > lambda; // temporary data for banded solve
    vector<vector<double> > lambdap;

    Lagrange_interpolation_coefficients LG;
    
    //===========================================================
    // workspace for CODE solver
    //===========================================================
    vector<double> dyi;
    vector<double> ddyi;    
    vector<double> Gi;
    vector<double> Hi;    
    vector<double> Gip;
    vector<double> Hip; 
    //
    vector<double> *Hi_ptr; 
    vector<double> *Gi_ptr; 
    vector<double> *Hip_ptr; 
    vector<double> *Gip_ptr; 
};

//==================================================================================================
void init_PDE_Stepper_Data(PDE_Stepper_Data &PDE_D, int npts);
void reset_PDE_solver_variables();

//==================================================================================================
void setup_Lagrange_interpolation_coefficients_O1(PDE_Stepper_Data &PDE_D, vector<double> &xarr);
void setup_Lagrange_interpolation_coefficients_O2(PDE_Stepper_Data &PDE_D, vector<double> &xarr);

void Step_PDE_O1(double zs, double ze, const vector<double> &xi, vector<double> &yi, 
                 double yi_low, double yi_up, PDE_Stepper_Data &PDE_D,
                 void (*func)(double z, const vector<double> &xi, vector<double> &Ai, 
                              vector<double> &Bi, vector<double> &Ci, vector<double> &Di));

void Step_PDE_O2t(double theta, double zs, double ze, const vector<double> &xi, vector<double> &yi, 
                  double yi_low, double yi_up, PDE_Stepper_Data &PDE_D,
                  void (*func)(double z, const vector<double> &xi, vector<double> &Ai, 
                               vector<double> &Bi, vector<double> &Ci, vector<double> &Di),
                  bool iteration_mode=0);

//==================================================================================================
// PDE solver with Integral part as correction 
// added 24.02.2011
//==================================================================================================
void setup_Lagrange_interpolation_coefficients_O2_Int(PDE_Stepper_Data &PDE_D, 
                                                      vector<double> &xarr);

void Step_PDE_O2t_Int(double theta, int it_max, double zs, double ze, 
                      const vector<double> &xi, vector<double> &yi, 
                      double yi_low, double yi_up, PDE_Stepper_Data &PDE_D,
                      //
                      void (*func)(double z, const vector<double> &xi, 
                                   vector<double> &Ai, vector<double> &Bi, 
                                   vector<double> &Ci, vector<double> &Di),
                      //
                      void (*corr_func)(double z, 
                                        const vector<double> &xi, 
                                        const vector<double> &yi,
                                        vector<double> &corr_v));

//==================================================================================================
// PDE setup for coupled ODE problem (at the moment only one PDE)
// added 27.06.2011
//==================================================================================================
void setup_Lagrange_interpolation_coefficients_O2_CODE(PDE_Stepper_Data &PDE_D, 
                                                       vector<double> &xarr, 
                                                       bool setup_Integral);

void compute_derivatives(const vector<double> &yi, 
                         const vector<double> &Ai, const vector<double> &Bi, 
                         PDE_Stepper_Data &PDE_D, 
                         vector<double> &dy, vector<double> &ddy);

#endif

//==================================================================================================
//==================================================================================================

