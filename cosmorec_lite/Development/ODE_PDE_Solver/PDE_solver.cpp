//==================================================================================================
//
// Author Jens Chluba July 2010
// last modification:  Aug 2012
//
//==================================================================================================
// 27.06.2011: added PDE setup for integro-differential equation systems
// 24.02.2011: added PDE solver with Integral part as correction 
// 01.02.2011: added possibility to take integral part over +/-2 off-diags into account

//==================================================================================================
//
// Purpose: Solve linear (!) parabolic differential equation of the form
//
// du/dt = A(x,t) d^2u/dx^2 + B(x,t) du/dx + C(x,t) u + D(x,t)
//
// The mesh of grid-points given by xi[...] is not assumed to be uniform.
//
//==================================================================================================

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>

#include "PDE_solver.h"
#include "routines.h"

using namespace std;

int define_PDEsolver_verbosity=0;

//==================================================================================================
//
// global structure
//
//==================================================================================================
struct previous_solution_info
{
    vector<double> Fsol;
    bool is_set_up;
    double zsol;
};

bool PDE_solver_initial_call_of_solver=0;

//==================================================================================================
void init_PDE_Stepper_Data(PDE_Stepper_Data &PDE_D, int npts)
{
    PDE_D.ptr=NULL;
    //
    PDE_D.Ai.resize(npts);
    PDE_D.Bi.resize(npts);
    PDE_D.Ci.resize(npts);
    PDE_D.Di.resize(npts);
    PDE_D.Ii.resize(npts);
    //
    PDE_D.Ai_ptr=&PDE_D.Ai;
    PDE_D.Bi_ptr=&PDE_D.Bi;
    PDE_D.Ci_ptr=&PDE_D.Ci;
    PDE_D.Di_ptr=&PDE_D.Di;
    PDE_D.Ii_ptr=&PDE_D.Ii;
    //
    PDE_D.Ui.resize(npts);
    PDE_D.Vi.resize(npts);
    PDE_D.bi.resize(npts);
    PDE_D.zi.resize(npts);
    //
    PDE_D.Aip.resize(npts);
    PDE_D.Bip.resize(npts);
    PDE_D.Cip.resize(npts);
    PDE_D.Dip.resize(npts);
    PDE_D.Iip.resize(npts);
    //
    PDE_D.Aip_ptr=&PDE_D.Aip;
    PDE_D.Bip_ptr=&PDE_D.Bip;
    PDE_D.Cip_ptr=&PDE_D.Cip;
    PDE_D.Dip_ptr=&PDE_D.Dip;
    PDE_D.Iip_ptr=&PDE_D.Iip;
        
    PDE_D.lambda. resize(npts, vector<double> (5));
    PDE_D.lambdap.resize(npts, vector<double> (5));

    return;
}

//==================================================================================================
//
// reset solver
//
//==================================================================================================
void reset_PDE_solver_variables()
{
    PDE_solver_initial_call_of_solver=0;
    return;
}

//-------------------------------------------------------------------------------------
// 1. derivative
//-------------------------------------------------------------------------------------
double dli_dx_xj(int i, int j, double *xarr, int np)
{
    double dl_dx=0.0;
    
    if(i==j) 
    {
        for(int k=0; k<i; k++)    dl_dx+=1.0/(xarr[i]-xarr[k]);
        for(int k=i+1; k<np; k++) dl_dx+=1.0/(xarr[i]-xarr[k]);
    }
    else if(i<j) 
    {
        dl_dx=1.0/(xarr[i]-xarr[j]);
        for(int k=0; k<i; k++)    dl_dx*=(xarr[j]-xarr[k])/(xarr[i]-xarr[k]);
        for(int k=i+1; k<j; k++)  dl_dx*=(xarr[j]-xarr[k])/(xarr[i]-xarr[k]);
        for(int k=j+1; k<np; k++) dl_dx*=(xarr[j]-xarr[k])/(xarr[i]-xarr[k]);
    }
    else 
    {
        dl_dx=1.0/(xarr[i]-xarr[j]);
        for(int k=0; k<j; k++)    dl_dx*=(xarr[j]-xarr[k])/(xarr[i]-xarr[k]);
        for(int k=j+1; k<i; k++)  dl_dx*=(xarr[j]-xarr[k])/(xarr[i]-xarr[k]);
        for(int k=i+1; k<np; k++) dl_dx*=(xarr[j]-xarr[k])/(xarr[i]-xarr[k]);
    }
    
    return dl_dx;
}

//-------------------------------------------------------------------------------------
// 2. derivative
//-------------------------------------------------------------------------------------
double d2li_d2x_xj(int i, int j, double *xarr, int np)
{
    double d2l_d2x=0.0, dum;
    
    if(i==j) 
    {
        for(int k=0; k<i; k++)
        {
            dum=0.0;
            
            for(int m=0; m<k; m++) dum+=1.0/(xarr[i]-xarr[m]);
            for(int m=k+1; m<i; m++) dum+=1.0/(xarr[i]-xarr[m]);
            for(int m=i+1; m<np; m++) dum+=1.0/(xarr[i]-xarr[m]);
            
            d2l_d2x+=1.0/(xarr[i]-xarr[k])*dum;
        }
        
        for(int k=i+1; k<np; k++)
        {
            dum=0.0;
            
            for(int m=0; m<i; m++) dum+=1.0/(xarr[i]-xarr[m]);
            for(int m=i+1; m<k; m++) dum+=1.0/(xarr[i]-xarr[m]);
            for(int m=k+1; m<np; m++) dum+=1.0/(xarr[i]-xarr[m]);
            
            d2l_d2x+=1.0/(xarr[i]-xarr[k])*dum;
        }
    }
    else 
    {
        for(int m=0; m<np; m++) 
        {
            if(m!=i && m!=j)
            {
                dum=1.0;
                for(int k=0; k<np; k++) 
                    if(k!=i && k!=j && k!=m) dum*=(xarr[j]-xarr[k])/(xarr[i]-xarr[k]);
                d2l_d2x+=dum/(xarr[i]-xarr[m]);
            }
        }
        
        d2l_d2x*=2.0/(xarr[i]-xarr[j]);
    }
    
    return d2l_d2x;
}

//--------------------------------------------------------------------------------------------------
void print_mess_PDE(string mess)
{
    if(define_PDEsolver_verbosity>0) cout << mess << endl;   
    return;
}

//==================================================================================================
//
// 1. order scheme
//
//==================================================================================================

//-------------------------------------------------------------------------------------
// precompute coefficients
//-------------------------------------------------------------------------------------
void setup_Lagrange_interpolation_coefficients_O1(PDE_Stepper_Data &PDE_D, vector<double> &xarr)
{
    int np=xarr.size();
    
    print_mess_PDE("\n setup_Lagrange_interpolation_coefficients_O1:: setting up grid ");
    
    PDE_D.LG.dli_dx.clear();
    PDE_D.LG.d2li_d2x.clear();
    
    //====================================================================
    // allocate memory
    //====================================================================
    vector<double> dum(3, 0.0);
    for(int k=0; k<np; k++)
    {
        PDE_D.LG.dli_dx  .push_back(dum);
        PDE_D.LG.d2li_d2x.push_back(dum);
    }
    
    //====================================================================
    // all intereor points
    //====================================================================
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int i=1; i<np-1; i++)
        for(int j=0; j<3; j++)
        {
            PDE_D.LG.dli_dx  [i][j]=dli_dx_xj  (j, 1, &xarr[i-1], 3);
            PDE_D.LG.d2li_d2x[i][j]=d2li_d2x_xj(j, 1, &xarr[i-1], 3);
        }
    
    print_mess_PDE(" setup_Lagrange_interpolation_coefficients_O1:: done ");

    return;
}


//-------------------------------------------------------------------------------------
// compute lambda_i
//-------------------------------------------------------------------------------------
void compute_lambda_i_01(double dz, int i, vector<double> &Ai, 
                         vector<double> &Bi, vector<double> &Ci, 
                         PDE_Stepper_Data &PDE_D, double *lambda)
{       
    for(int j=0; j<3; j++) lambda[j]=PDE_D.LG.d2li_d2x[i][j]*Ai[i]+PDE_D.LG.dli_dx[i][j]*Bi[i];
    lambda[1]+=Ci[i];
    for(int j=0; j<3; j++) lambda[j]*=-dz;
    lambda[1]+=1.0;
}

//-------------------------------------------------------------------------------------
// compute solution
//-------------------------------------------------------------------------------------
void Step_PDE_O1(double zs, double ze, const vector<double> &xi, vector<double> &yi, 
                 double yi_low, double yi_up, PDE_Stepper_Data &PDE_D,
                 void (*func)(double z, const vector<double> &xi, vector<double> &Ai, 
                              vector<double> &Bi, vector<double> &Ci, vector<double> &Di))
{
    double dz=ze-zs;
    int npxi=xi.size();
    
    //=============================================
    // evaluate coefficients of PDE at future time
    //=============================================
    func(ze, xi, PDE_D.Ai, PDE_D.Bi, PDE_D.Ci, PDE_D.Di);
    
    //=============================================
    // lower & upper boundary condition
    //=============================================
    yi[0]=yi_low;
    yi[npxi-1]=yi_up;
    
    //=============================================
    // compute Ui* and bi*
    //=============================================
    double lambda[3], dum=0.0; lambda[2]=0.0;
    //=============================================
    PDE_D.Ui[0]=0.0;
    PDE_D.bi[0]=yi[0];  
    //=============================================
    for(int i=1; i<npxi-1; i++)
    {
        compute_lambda_i_01(dz, i, PDE_D.Ai, PDE_D.Bi, PDE_D.Ci, PDE_D, lambda);
        //
        dum=(lambda[1]-lambda[0]*PDE_D.Ui[i-1]);
        PDE_D.Ui[i]=lambda[2]/dum;
        PDE_D.bi[i]=(yi[i]+dz*PDE_D.Di[i]-lambda[0]*PDE_D.bi[i-1])/dum;
    }
    //=============================================
    PDE_D.bi[npxi-2]-=yi[npxi-1]*lambda[2]/dum;
    
    //=============================================
    // compute solution
    //=============================================
    yi[npxi-2]=PDE_D.bi[npxi-2];
    for(int i=npxi-3; i>0; i--) yi[i]=PDE_D.bi[i]-PDE_D.Ui[i]*yi[i+1];  
    
    return;
}


//==================================================================================================
//
// 2. order scheme (in spacial derivative)
//
//==================================================================================================

//-------------------------------------------------------------------------------------
// precompute coefficients
//-------------------------------------------------------------------------------------
void setup_Lagrange_interpolation_coefficients_O2(PDE_Stepper_Data &PDE_D, vector<double> &xarr)
{
    int np=xarr.size();

    print_mess_PDE("\n setup_Lagrange_interpolation_coefficients_O2:: setting up grid ");
    
    PDE_D.LG.dli_dx.clear();
    PDE_D.LG.d2li_d2x.clear();

    //====================================================================
    // allocate memory
    //====================================================================
    vector<double> dum(5, 0.0);
    for(int k=0; k<np; k++)
    {
        PDE_D.LG.dli_dx  .push_back(dum);
        PDE_D.LG.d2li_d2x.push_back(dum);
    }
    
    //====================================================================
    // 1. interior point (lower boundary not needed)
    //====================================================================
    int i=1;
    for(int j=0; j<5; j++)
    {
        PDE_D.LG.dli_dx  [i][j]=dli_dx_xj  (j, 1, &xarr[i-1], 5);
        PDE_D.LG.d2li_d2x[i][j]=d2li_d2x_xj(j, 1, &xarr[i-1], 5);
    }
    
    //====================================================================
    // last interior point (upper boundary not needed)
    //====================================================================
    i=np-2;
    for(int j=0; j<5; j++)
    {
        PDE_D.LG.dli_dx  [i][j]=dli_dx_xj  (j, 3, &xarr[i-3], 5);
        PDE_D.LG.d2li_d2x[i][j]=d2li_d2x_xj(j, 3, &xarr[i-3], 5);
    }

    //====================================================================
    // the rest
    //====================================================================
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(i=2; i<np-2; i++)
        for(int j=0; j<5; j++)
        {
            PDE_D.LG.dli_dx  [i][j]=dli_dx_xj  (j, 2, &xarr[i-2], 5);
            PDE_D.LG.d2li_d2x[i][j]=d2li_d2x_xj(j, 2, &xarr[i-2], 5);
        }

    print_mess_PDE(" setup_Lagrange_interpolation_coefficients_O2:: done ");

    return;
}
        

//-------------------------------------------------------------------------------------
// compute lambda_i
//-------------------------------------------------------------------------------------
void compute_lambda_i_02(const double dz, const int i, const int center, 
                         const vector<double> &Ai, const vector<double> &Bi, 
                         const vector<double> &Ci, 
                         PDE_Stepper_Data &PDE_D, vector<double> &lambda)
{       
    for(int j=0; j<5; j++) lambda[j]=PDE_D.LG.d2li_d2x[i][j]*Ai[i]+PDE_D.LG.dli_dx[i][j]*Bi[i];
    lambda[center]+=Ci[i];
    lambda[center]-=1.0/(dz+1.0e-100);
}

//-------------------------------------------------------------------------------------
// compute solution
//
// 27.04.2011: rearranged some of the evaluation to make things more compact
// 22.04.2011: Added simple openmp support for matrix element evaluations
// 28.02.2011: Found small bug in the matrix solving part. Affected the boundary.
// 24.02.2011: This routine was optimized to reduce the number of operations per call
//
//-------------------------------------------------------------------------------------
void Step_PDE_O2t(double theta, double zs, double ze, 
                  const vector<double> &xi, vector<double> &yi, 
                  double yi_low, double yi_up, PDE_Stepper_Data &PDE_D,
                  void (*func)(double z, const vector<double> &xi, 
                               vector<double> &Ai, vector<double> &Bi, 
                               vector<double> &Ci, vector<double> &Di),
                  bool iteration_mode)
{
    double kap=theta-1.0, eta=kap/theta;
    double dz=ze-zs;
    int npxi=xi.size();
    
    //=================================================================================
    // evaluate coefficients of PDE at future time
    //================================================================================= 
    if(!PDE_solver_initial_call_of_solver)
    {
        func(zs, xi, *PDE_D.Aip_ptr, *PDE_D.Bip_ptr, *PDE_D.Cip_ptr, *PDE_D.Dip_ptr);
        PDE_solver_initial_call_of_solver=1;
    }
    
    func(ze, xi, *PDE_D.Ai_ptr, *PDE_D.Bi_ptr, *PDE_D.Ci_ptr, *PDE_D.Di_ptr);
    
    //=================================================================================
    // compute all matrix elements
    //=================================================================================
    int i=1;
    compute_lambda_i_02(dz*theta, i, 1, *PDE_D.Ai_ptr, *PDE_D.Bi_ptr, 
                        *PDE_D.Ci_ptr , PDE_D, PDE_D.lambda [i]);
    
    compute_lambda_i_02(dz*kap  , i, 1, *PDE_D.Aip_ptr,*PDE_D.Bip_ptr,
                        *PDE_D.Cip_ptr, PDE_D, PDE_D.lambdap[i]); 
    
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(i=2; i<npxi-2; i++)
    {
        compute_lambda_i_02(dz*theta, i, 2, *PDE_D.Ai_ptr, *PDE_D.Bi_ptr, 
                            *PDE_D.Ci_ptr , PDE_D, PDE_D.lambda [i]);
        
        compute_lambda_i_02(dz*kap  , i, 2, *PDE_D.Aip_ptr,*PDE_D.Bip_ptr,
                            *PDE_D.Cip_ptr, PDE_D, PDE_D.lambdap[i]); 
    }
    
    i=npxi-2;
    compute_lambda_i_02(dz*theta, i, 3, *PDE_D.Ai_ptr, *PDE_D.Bi_ptr, 
                        *PDE_D.Ci_ptr , PDE_D, PDE_D.lambda [i]);
    
    compute_lambda_i_02(dz*kap  , i, 3, *PDE_D.Aip_ptr,*PDE_D.Bip_ptr,
                        *PDE_D.Cip_ptr, PDE_D, PDE_D.lambdap[i]); 

    //================================================================================= 
    // initialize all b-vectors with common parts
    //================================================================================= 
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(i=1; i<npxi-1; i++) PDE_D.bi[i]=eta*(*PDE_D.Dip_ptr)[i]-(*PDE_D.Di_ptr)[i];

    //================================================================================= 
    i=1;
    for(int m=0; m<5; m++) PDE_D.bi[i]+=PDE_D.lambdap[i][m]*eta*yi[i+m-1];
    PDE_D.bi[i]-=PDE_D.lambda[i][0]*yi_low; // JC changed sign 28.02.2011

    i=2;
    PDE_D.bi[i]-=PDE_D.lambda[i][0]*yi_low; // JC changed sign 28.02.2011

#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(i=2; i<npxi-2; i++) 
    {
        for(int m=0; m<5; m++) PDE_D.bi[i]+=PDE_D.lambdap[i][m]*eta*yi[i+m-2];
    }

    i=npxi-3;
    PDE_D.bi[i]-=PDE_D.lambda[i][4]*yi_up;
    
    i=npxi-2;
    for(int m=0; m<5; m++) PDE_D.bi[i]+=PDE_D.lambdap[i][m]*eta*yi[i+m-3];
    PDE_D.bi[i]-=PDE_D.lambda[i][4]*yi_up;
    //================================================================================= 
    
    //================================================================================= 
    // compute Ui* and bi*
    //================================================================================= 
    double dum1, dum2, dum3;
    double Zlow;
    
    //================================================================================= 
    // 1. intereor point
    //=================================================================================  
    i=1;
    //
    PDE_D.bi[i]/=PDE_D.lambda[i][1];
    PDE_D.Ui[i] =PDE_D.lambda[i][2]/PDE_D.lambda[i][1];
    PDE_D.Vi[i] =PDE_D.lambda[i][3]/PDE_D.lambda[i][1];
    Zlow        =PDE_D.lambda[i][4]/PDE_D.lambda[i][1];
    
    //================================================================================= 
    // 2. intereor point
    //=================================================================================  
    i=2;
    //
    dum2=PDE_D.lambda[i][2]-PDE_D.lambda[i][1]*PDE_D.Ui[i-1];
    //
    PDE_D.bi[i]-=PDE_D.lambda[i][1]*PDE_D.bi[i-1];
    PDE_D.bi[i]/=dum2;
    //
    PDE_D.Ui[i] =(PDE_D.lambda[i][3]-PDE_D.lambda[i][1]*PDE_D.Vi[i-1])/dum2;
    PDE_D.Vi[i] =(PDE_D.lambda[i][4]-PDE_D.lambda[i][1]*Zlow)/dum2;
    
    //================================================================================= 
    // 3. intereor point
    //=================================================================================  
    i=3;
    //
    dum1=PDE_D.lambda[i][1]-PDE_D.lambda[i][0]*PDE_D.Ui[i-2];
    dum2=PDE_D.lambda[i][2]-PDE_D.lambda[i][0]*PDE_D.Vi[i-2];
    //
    PDE_D.bi[i]-=PDE_D.lambda[i][0]*PDE_D.bi[i-2];
    PDE_D.Ui[i] =PDE_D.lambda[i][3]-PDE_D.lambda[i][0]*Zlow;
    //---------------------------------------------
    dum2-=dum1*PDE_D.Ui[i-1];
    //
    PDE_D.bi[i]-=dum1*PDE_D.bi[i-1];
    PDE_D.bi[i]/=dum2;
    //
    PDE_D.Ui[i]-=dum1*PDE_D.Vi[i-1];
    PDE_D.Ui[i]/=dum2;
    PDE_D.Vi[i] =PDE_D.lambda[i][4]/dum2;
    
    //================================================================================= 
    // intereor point 4...n-4
    //=================================================================================  
    for(i=4; i<npxi-3; i++)
    {
        dum1=PDE_D.lambda[i][1]-PDE_D.lambda[i][0]*PDE_D.Ui[i-2];
        dum2=PDE_D.lambda[i][2]-PDE_D.lambda[i][0]*PDE_D.Vi[i-2];
        //
        PDE_D.bi[i]-=PDE_D.lambda[i][0]*PDE_D.bi[i-2];
        PDE_D.Ui[i] =PDE_D.lambda[i][3];
        //---------------------------------------------
        dum2-=dum1*PDE_D.Ui[i-1];
        //
        PDE_D.bi[i]-=dum1*PDE_D.bi[i-1];
        PDE_D.bi[i]/=dum2;
        //
        PDE_D.Ui[i]-=dum1*PDE_D.Vi[i-1];
        PDE_D.Ui[i]/=dum2;
        PDE_D.Vi[i] =PDE_D.lambda[i][4]/dum2;
    }
    
    //================================================================================= 
    // intereor point n-3
    //=================================================================================  
    i=npxi-3;
    //
    dum1=PDE_D.lambda[i][1]-PDE_D.lambda[i][0]*PDE_D.Ui[i-2];
    dum2=PDE_D.lambda[i][2]-PDE_D.lambda[i][0]*PDE_D.Vi[i-2];
    //
    PDE_D.bi[i]-=PDE_D.lambda[i][0]*PDE_D.bi[i-2];
    PDE_D.Ui[i] =PDE_D.lambda[i][3];
    //---------------------------------------------
    dum2-=dum1*PDE_D.Ui[i-1];
    //
    PDE_D.bi[i]-=dum1*PDE_D.bi[i-1];
    PDE_D.bi[i]/=dum2;
    //
    PDE_D.Ui[i]-=dum1*PDE_D.Vi[i-1];
    PDE_D.Ui[i]/=dum2;
    PDE_D.Vi[i] =0.0;
    
    //================================================================================= 
    // intereor point n-2
    //=================================================================================  
    i=npxi-2;
    //
    dum1=PDE_D.lambda[i][1]-PDE_D.lambda[i][0]*PDE_D.Ui[i-3];
    dum2=PDE_D.lambda[i][2]-PDE_D.lambda[i][0]*PDE_D.Vi[i-3];
    dum3=PDE_D.lambda[i][3];
    //
    PDE_D.bi[i]-=PDE_D.lambda[i][0]*PDE_D.bi[i-3];
    //---------------------------------------------
    dum2-=dum1*PDE_D.Ui[i-2];
    dum3-=dum1*PDE_D.Vi[i-2];
    //
    PDE_D.bi[i]-=dum1*PDE_D.bi[i-2];
    //---------------------------------------------
    dum3-=dum2*PDE_D.Ui[i-1];
    //
    PDE_D.bi[i]-=dum2*PDE_D.bi[i-1];
    PDE_D.bi[i]/=dum3;
    PDE_D.Ui[i]=PDE_D.Vi[i]=0.0;
    
    //================================================================================= 
    // lower & upper boundary condition
    //================================================================================= 
    yi[0]=yi_low;
    yi[npxi-1]=yi_up;
    
    //================================================================================= 
    // compute solution
    //================================================================================= 
    yi[npxi-2]=PDE_D.bi[npxi-2];
    //
    for(i=npxi-3; i>0; i--) yi[i]=PDE_D.bi[i]-PDE_D.Ui[i]*yi[i+1]-PDE_D.Vi[i]*yi[i+2];
    //
    yi[1]-=Zlow*yi[4];  
    
    //================================================================================= 
    // copy old coefficients
    //================================================================================= 
    if(!iteration_mode)
    {
        PDE_D.ptr=PDE_D.Aip_ptr;
        PDE_D.Aip_ptr=PDE_D.Ai_ptr;
        PDE_D.Ai_ptr=PDE_D.ptr;
        //
        PDE_D.ptr=PDE_D.Bip_ptr;
        PDE_D.Bip_ptr=PDE_D.Bi_ptr;
        PDE_D.Bi_ptr=PDE_D.ptr;
        //
        PDE_D.ptr=PDE_D.Cip_ptr;
        PDE_D.Cip_ptr=PDE_D.Ci_ptr;
        PDE_D.Ci_ptr=PDE_D.ptr;
        //  
        PDE_D.ptr=PDE_D.Dip_ptr;
        PDE_D.Dip_ptr=PDE_D.Di_ptr;
        PDE_D.Di_ptr=PDE_D.ptr;
        //
        PDE_D.ptr=NULL;
    }
    
    return;
}




//==================================================================================================
//
// 2. order scheme with correction from integral
//
//==================================================================================================

//-------------------------------------------------------------------------------------
// precompute coefficients
//-------------------------------------------------------------------------------------
void setup_Lagrange_interpolation_coefficients_O2_Int(PDE_Stepper_Data &PDE_D, 
                                                      vector<double> &xarr)
{
    int np=xarr.size();
    
    print_mess_PDE("\n setup_Lagrange_interpolation_coefficients_O2_Int:: setting up grid ");
    
    PDE_D.LG.dli_dx.clear();
    PDE_D.LG.d2li_d2x.clear();
    
    //====================================================================
    // allocate memory
    //====================================================================
    vector<double> dum(5, 0.0);
    vector<double> dum2(np, 0.0);

    for(int k=0; k<np; k++)
    {
        PDE_D.LG.dli_dx  .push_back(dum);
        PDE_D.LG.d2li_d2x.push_back(dum);
    }
    
    //====================================================================
    // 1. interior point (lower boundary not needed)
    //====================================================================
    int i=1;
    for(int j=0; j<5; j++)
    {
        PDE_D.LG.dli_dx  [i][j]=dli_dx_xj  (j, 1, &xarr[i-1], 5);
        PDE_D.LG.d2li_d2x[i][j]=d2li_d2x_xj(j, 1, &xarr[i-1], 5);
    }
    
    //====================================================================
    // last interior point (upper boundary not needed)
    //====================================================================
    i=np-2;
    for(int j=0; j<5; j++)
    {
        PDE_D.LG.dli_dx  [i][j]=dli_dx_xj  (j, 3, &xarr[i-3], 5);
        PDE_D.LG.d2li_d2x[i][j]=d2li_d2x_xj(j, 3, &xarr[i-3], 5);
    }
    
    //====================================================================
    // the rest
    //====================================================================
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(i=2; i<np-2; i++)
        for(int j=0; j<5; j++)
        {
            PDE_D.LG.dli_dx  [i][j]=dli_dx_xj  (j, 2, &xarr[i-2], 5);
            PDE_D.LG.d2li_d2x[i][j]=d2li_d2x_xj(j, 2, &xarr[i-2], 5);
        }
    
    print_mess_PDE(" setup_Lagrange_interpolation_coefficients_O2_Int:: done ");
    
    return;
}

//==================================================================================================
//
// 2. order scheme with correction from integral; Iterative scheme
//
//==================================================================================================

void Step_PDE_O2t_Int(double theta, int it_max, double zs, double ze, 
                      const vector<double> &xi, vector<double> &yi, 
                      double yi_low, double yi_up, PDE_Stepper_Data &PDE_D,
                      //
                      void (*func)(double z, const vector<double> &xi, 
                                   vector<double> &Ai, vector<double> &Bi, 
                                   vector<double> &Ci, vector<double> &Di),
                      //
                      void (*corr_func)(double z, const vector<double> &xi, 
                                        const vector<double> &yi,
                                        vector<double> &corr_v)
                      )
{
    double kap=theta-1.0, eta=kap/theta;
    double dz=ze-zs;
    int npxi=xi.size();
    
    //=================================================================================
    // evaluate coefficients of PDE at future time
    //=================================================================================
    if(!PDE_solver_initial_call_of_solver)
    {
        func(zs, xi, *PDE_D.Aip_ptr, *PDE_D.Bip_ptr, *PDE_D.Cip_ptr, *PDE_D.Dip_ptr);
        PDE_solver_initial_call_of_solver=1;
        
        for(int i=0; i<npxi; i++) (*PDE_D.Iip_ptr)[i]=(*PDE_D.Ii_ptr)[i]=0.0;
    }
    
    func(ze, xi, *PDE_D.Ai_ptr, *PDE_D.Bi_ptr, *PDE_D.Ci_ptr, *PDE_D.Di_ptr);
    
    //=================================================================================
    // compute all matrix elements
    //=================================================================================
    int i=1;
    compute_lambda_i_02(dz*theta, i, 1, *PDE_D.Ai_ptr, *PDE_D.Bi_ptr, 
                        *PDE_D.Ci_ptr , PDE_D, PDE_D.lambda [i]);
    
    compute_lambda_i_02(dz*kap  , i, 1, *PDE_D.Aip_ptr,*PDE_D.Bip_ptr,
                        *PDE_D.Cip_ptr, PDE_D, PDE_D.lambdap[i]); 
    
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(i=2; i<npxi-2; i++)
    {
        compute_lambda_i_02(dz*theta, i, 2, *PDE_D.Ai_ptr, *PDE_D.Bi_ptr, 
                            *PDE_D.Ci_ptr , PDE_D, PDE_D.lambda [i]);
        
        compute_lambda_i_02(dz*kap  , i, 2, *PDE_D.Aip_ptr,*PDE_D.Bip_ptr,
                            *PDE_D.Cip_ptr, PDE_D, PDE_D.lambdap[i]); 
    }
    
    i=npxi-2;
    compute_lambda_i_02(dz*theta, i, 3, *PDE_D.Ai_ptr, *PDE_D.Bi_ptr, 
                        *PDE_D.Ci_ptr , PDE_D, PDE_D.lambda [i]);
    
    compute_lambda_i_02(dz*kap  , i, 3, *PDE_D.Aip_ptr,*PDE_D.Bip_ptr,
                        *PDE_D.Cip_ptr, PDE_D, PDE_D.lambdap[i]); 
    
    //=================================================================================
    // compute Ui* and bi*
    //=================================================================================
    double dum1, dum2, dum3;
    double Zlow;
    //=================================================================================
    for(int it=0; it<=it_max; it++)
    {        
        //==============================================================================
        // 1. intereor point
        //==============================================================================
        i=1;

        if(it==0) 
        {
            PDE_D.bi[i]=PDE_D.lambdap[i][0]*eta*yi[0];
            for(int m=1; m<5; m++) PDE_D.bi[i]+=PDE_D.lambdap[i][m]*eta*yi[i+m-1];
            PDE_D.bi[i]+=eta*(*PDE_D.Dip_ptr)[i]-(*PDE_D.Di_ptr)[i];
            PDE_D.bi[i]-=PDE_D.lambda[i][0]*yi_low; // JC changed sign 28.02.2011
            PDE_D.bi[i]+=eta*(*PDE_D.Iip_ptr)[i];
        }
        else PDE_D.bi[i]=-(*PDE_D.Ii_ptr)[i];        

        PDE_D.bi[i]/=PDE_D.lambda[i][1];
        PDE_D.Ui[i] =PDE_D.lambda[i][2]/PDE_D.lambda[i][1];
        PDE_D.Vi[i] =PDE_D.lambda[i][3]/PDE_D.lambda[i][1];
        Zlow        =PDE_D.lambda[i][4]/PDE_D.lambda[i][1];
        
        //==============================================================================
        // 2. intereor point
        //==============================================================================
        i=2;

        if(it==0) 
        {
            PDE_D.bi[i]=PDE_D.lambdap[i][0]*eta*yi[0];
            for(int m=1; m<5; m++) PDE_D.bi[i]+=PDE_D.lambdap[i][m]*eta*yi[i+m-2];
            PDE_D.bi[i]+=eta*(*PDE_D.Dip_ptr)[i]-(*PDE_D.Di_ptr)[i];
            PDE_D.bi[i]-=PDE_D.lambda[i][0]*yi_low; // JC changed sign 28.02.2011
            PDE_D.bi[i]+=eta*(*PDE_D.Iip_ptr)[i];
        }
        else PDE_D.bi[i]=-(*PDE_D.Ii_ptr)[i];        
        
        dum2=PDE_D.lambda[i][2]-PDE_D.lambda[i][1]*PDE_D.Ui[i-1];
        //
        PDE_D.bi[i]-=PDE_D.lambda[i][1]*PDE_D.bi[i-1];
        PDE_D.bi[i]/=dum2;
        PDE_D.Ui[i] =(PDE_D.lambda[i][3]-PDE_D.lambda[i][1]*PDE_D.Vi[i-1])/dum2;
        PDE_D.Vi[i] =(PDE_D.lambda[i][4]-PDE_D.lambda[i][1]*Zlow)/dum2;
        
        //==============================================================================
        // 3. intereor point
        //==============================================================================
        i=3;

        if(it==0) 
        {
            PDE_D.bi[i]=PDE_D.lambdap[i][0]*eta*yi[i-2];
            for(int m=1; m<5; m++) PDE_D.bi[i]+=PDE_D.lambdap[i][m]*eta*yi[i+m-2];
            PDE_D.bi[i]+=eta*(*PDE_D.Dip_ptr)[i]-(*PDE_D.Di_ptr)[i];
            PDE_D.bi[i]+=eta*(*PDE_D.Iip_ptr)[i];
        }
        else PDE_D.bi[i]=-(*PDE_D.Ii_ptr)[i];        
        
        dum1=PDE_D.lambda[i][1]-PDE_D.lambda[i][0]*PDE_D.Ui[i-2];
        dum2=PDE_D.lambda[i][2]-PDE_D.lambda[i][0]*PDE_D.Vi[i-2];
        //
        PDE_D.bi[i]-=PDE_D.lambda[i][0]*PDE_D.bi[i-2];
        PDE_D.Ui[i] =PDE_D.lambda[i][3]-PDE_D.lambda[i][0]*Zlow;
        //
        dum2-=dum1*PDE_D.Ui[i-1];
        //
        PDE_D.bi[i]-=dum1*PDE_D.bi[i-1];
        PDE_D.bi[i]/=dum2;
        PDE_D.Ui[i]-=dum1*PDE_D.Vi[i-1];
        PDE_D.Ui[i]/=dum2;
        PDE_D.Vi[i] =PDE_D.lambda[i][4]/dum2;
        
        //==============================================================================
        // intereor point 4...n-4
        //==============================================================================
        if(it==0) for(i=4; i<npxi-3; i++)
        {
            PDE_D.bi[i]=PDE_D.lambdap[i][0]*eta*yi[i-2];
            for(int m=1; m<5; m++) PDE_D.bi[i]+=PDE_D.lambdap[i][m]*eta*yi[i+m-2];
            PDE_D.bi[i]+=eta*(*PDE_D.Dip_ptr)[i]-(*PDE_D.Di_ptr)[i];
            PDE_D.bi[i]+=eta*(*PDE_D.Iip_ptr)[i];
            
            dum1=PDE_D.lambda[i][1]-PDE_D.lambda[i][0]*PDE_D.Ui[i-2];
            dum2=PDE_D.lambda[i][2]-PDE_D.lambda[i][0]*PDE_D.Vi[i-2];
            //
            PDE_D.bi[i]-=PDE_D.lambda[i][0]*PDE_D.bi[i-2];
            PDE_D.Ui[i] =PDE_D.lambda[i][3];
            //
            dum2-=dum1*PDE_D.Ui[i-1];
            //
            PDE_D.bi[i]-=dum1*PDE_D.bi[i-1];
            PDE_D.bi[i]/=dum2;
            PDE_D.Ui[i]-=dum1*PDE_D.Vi[i-1];
            PDE_D.Ui[i]/=dum2;
            PDE_D.Vi[i] =PDE_D.lambda[i][4]/dum2;
        }
        else for(i=4; i<npxi-3; i++)
        {
            PDE_D.bi[i]=-(*PDE_D.Ii_ptr)[i];
            
            dum1=PDE_D.lambda[i][1]-PDE_D.lambda[i][0]*PDE_D.Ui[i-2];
            dum2=PDE_D.lambda[i][2]-PDE_D.lambda[i][0]*PDE_D.Vi[i-2];
            //
            PDE_D.bi[i]-=PDE_D.lambda[i][0]*PDE_D.bi[i-2];
            PDE_D.Ui[i] =PDE_D.lambda[i][3];
            //
            dum2-=dum1*PDE_D.Ui[i-1];
            //
            PDE_D.bi[i]-=dum1*PDE_D.bi[i-1];
            PDE_D.bi[i]/=dum2;
            PDE_D.Ui[i]-=dum1*PDE_D.Vi[i-1];
            PDE_D.Ui[i]/=dum2;
            PDE_D.Vi[i] =PDE_D.lambda[i][4]/dum2;
        }        
        
        //==============================================================================
        // intereor point n-3
        //==============================================================================
        i=npxi-3;

        if(it==0) 
        {
            PDE_D.bi[i]=PDE_D.lambdap[i][0]*eta*yi[i-2];
            for(int m=1; m<5; m++) PDE_D.bi[i]+=PDE_D.lambdap[i][m]*eta*yi[i+m-2];
            PDE_D.bi[i]+=eta*(*PDE_D.Dip_ptr)[i]-(*PDE_D.Di_ptr)[i];
            PDE_D.bi[i]-=PDE_D.lambda[i][4]*yi_up;
            PDE_D.bi[i]+=eta*(*PDE_D.Iip_ptr)[i];
        }
        else PDE_D.bi[i]=-(*PDE_D.Ii_ptr)[i];        
        
        dum1=PDE_D.lambda[i][1]-PDE_D.lambda[i][0]*PDE_D.Ui[i-2];
        dum2=PDE_D.lambda[i][2]-PDE_D.lambda[i][0]*PDE_D.Vi[i-2];
        //
        PDE_D.bi[i]-=PDE_D.lambda[i][0]*PDE_D.bi[i-2];
        PDE_D.Ui[i] =PDE_D.lambda[i][3];
        PDE_D.Vi[i] =0.0;
        //
        dum2-=dum1*PDE_D.Ui[i-1];
        //
        PDE_D.bi[i]-=dum1*PDE_D.bi[i-1];
        PDE_D.bi[i]/=dum2;
        PDE_D.Ui[i]-=dum1*PDE_D.Vi[i-1];
        PDE_D.Ui[i]/=dum2;
        
        //==============================================================================
        // intereor point n-2
        //==============================================================================
        i=npxi-2;

        if(it==0) 
        {
            PDE_D.bi[i]=PDE_D.lambdap[i][0]*eta*yi[i-3];
            for(int m=1; m<5; m++) PDE_D.bi[i]+=PDE_D.lambdap[i][m]*eta*yi[i+m-3];
            PDE_D.bi[i]+=eta*(*PDE_D.Dip_ptr)[i]-(*PDE_D.Di_ptr)[i];
            PDE_D.bi[i]-=PDE_D.lambda[i][4]*yi_up;
            PDE_D.bi[i]+=eta*(*PDE_D.Iip_ptr)[i];
        }
        else PDE_D.bi[i]=-(*PDE_D.Ii_ptr)[i];        
        
        dum1=PDE_D.lambda[i][1]-PDE_D.lambda[i][0]*PDE_D.Ui[i-3];
        dum2=PDE_D.lambda[i][2]-PDE_D.lambda[i][0]*PDE_D.Vi[i-3];
        dum3=PDE_D.lambda[i][3];
        //
        PDE_D.bi[i]-=PDE_D.lambda[i][0]*PDE_D.bi[i-3];
        PDE_D.Ui[i]=PDE_D.Vi[i]=0.0;
        //
        dum2-=dum1*PDE_D.Ui[i-2];
        dum3-=dum1*PDE_D.Vi[i-2];
        //
        PDE_D.bi[i]-=dum1*PDE_D.bi[i-2];
        //
        dum3-=dum2*PDE_D.Ui[i-1];
        //
        PDE_D.bi[i]-=dum2*PDE_D.bi[i-1];
        PDE_D.bi[i]/=dum3;
        
        //==============================================================================
        // lower & upper boundary condition
        //==============================================================================
        yi[0]=yi_low;
        yi[npxi-1]=yi_up;
        
        //==============================================================================
        // compute solution
        //==============================================================================
        if(it==0){ PDE_D.zi[npxi-1]=yi[npxi-1]; PDE_D.zi[0]=yi[0]; } 
        else PDE_D.zi[npxi-1]=PDE_D.zi[0]=0.0;
        //
        PDE_D.zi[npxi-2]=PDE_D.bi[npxi-2];
        //
        for(i=npxi-3; i>0; i--)
            PDE_D.zi[i]=PDE_D.bi[i]-PDE_D.Ui[i]*PDE_D.zi[i+1]-PDE_D.Vi[i]*PDE_D.zi[i+2];
        //
        PDE_D.zi[1]-=Zlow*PDE_D.zi[4];  
        //
        if(it==0) for(i=1; i<npxi-1; i++) yi[i]=PDE_D.zi[i];
        else for(i=1; i<npxi-1; i++) yi[i]+=PDE_D.zi[i];
        
        //==============================================================================
        // compute correction from Integral part
        //==============================================================================
        corr_func(ze, xi, PDE_D.zi, *PDE_D.Ii_ptr);
        
        if(it==0) for(i=1; i<npxi-1; i++) (*PDE_D.Iip_ptr)[i]=(*PDE_D.Ii_ptr)[i];
        else for(i=1; i<npxi-1; i++) (*PDE_D.Iip_ptr)[i]+=(*PDE_D.Ii_ptr)[i];
    }
    
    //==============================================================================
    // swap old coefficients
    //==============================================================================
    PDE_D.ptr=PDE_D.Aip_ptr;
    PDE_D.Aip_ptr=PDE_D.Ai_ptr;
    PDE_D.Ai_ptr=PDE_D.ptr;
    //
    PDE_D.ptr=PDE_D.Bip_ptr;
    PDE_D.Bip_ptr=PDE_D.Bi_ptr;
    PDE_D.Bi_ptr=PDE_D.ptr;
    //
    PDE_D.ptr=PDE_D.Cip_ptr;
    PDE_D.Cip_ptr=PDE_D.Ci_ptr;
    PDE_D.Ci_ptr=PDE_D.ptr;
    //  
    PDE_D.ptr=PDE_D.Dip_ptr;
    PDE_D.Dip_ptr=PDE_D.Di_ptr;
    PDE_D.Di_ptr=PDE_D.ptr;
    //
    PDE_D.ptr=NULL;
    
    return;
}


//=====================================================================================
//
// 2. order scheme with correction from integral
//
//=====================================================================================

//-------------------------------------------------------------------------------------
// Integral aux functions
//-------------------------------------------------------------------------------------
double Int_fac(int j, double *xi)
{
    double dx=(xi[4] - xi[0]);
    double P=-1.0/60.0;
    for(int l=0  ; l<j; l++) P*=dx/(xi[l] - xi[j]);
    for(int l=j+1; l<5; l++) P*=dx/(xi[l] - xi[j]);
    return P;
}

//-------------------------------------------------------------------------------------
// 5 point formula for x0-x4
//-------------------------------------------------------------------------------------
double Int_li_xj_0(double *xi)
{ 
    // here j=0 (0, 1, 2, 3, 4)
    double dx1=xi[1]-xi[0];
    double dx2=xi[2]-xi[0];
    double dx3=xi[3]-xi[0];
    double dx4=xi[4]-xi[0];
    
    double r=(3.0*dx4-5.0*(dx1+dx2+dx3));
    r+=10.0*( dx2*dx3 + dx1*(dx2+dx3) )/dx4;
    r-=30.0*dx1*(dx2/dx4)*(dx3/dx4);     
    //
    r*=Int_fac(0, xi);
    
    return r;
}

//-------------------------------------------------------------------------------------
// 5 point formula for x0-x4
//-------------------------------------------------------------------------------------
double Int_li_xj_4(double *xi)
{ 
    // here j=4 (0, 1, 2, 3, 4)
    double dx1=xi[1]-xi[0];
    double dx2=xi[2]-xi[0];
    double dx3=xi[3]-xi[0];
    double dx4=xi[4]-xi[0];
    
    double r=12.0*dx4-15.0*(dx1+dx2+dx3);
    r+=20.0*( dx2*dx3 + dx1*(dx2+dx3) )/dx4;
    r-=30.0*dx1*(dx2/dx4)*(dx3/dx4);     
    //
    r*=-Int_fac(4, xi);
    
    return r;
}

//-------------------------------------------------------------------------------------
double Int_li_xj(int j, double *xi)
{ 
    //=================================================================
    if(j==0) return Int_li_xj_0(xi);
    if(j==4) return Int_li_xj_4(xi);
    //=================================================================
    
    //=================================================================
    // here j!=0 && j!=4 and the initial point is i=1, 2, 3
    //=================================================================
    double dx4=xi[4]-xi[0];
    double dx1=0.0, dx2=0.0;
    if(j==1){ dx1=xi[2]-xi[0];  dx2=xi[3]-xi[0]; }
    if(j==2){ dx1=xi[1]-xi[0];  dx2=xi[3]-xi[0]; }
    if(j==3){ dx1=xi[1]-xi[0];  dx2=xi[2]-xi[0]; }
    
    double r=3.0*dx4-5.0*(dx1+dx2);
    r+=10.0*dx1*(dx2/dx4);
    //
    r*=Int_fac(j, xi);
    
    return r;
}

//-------------------------------------------------------------------------------------
// precompute coefficients
//-------------------------------------------------------------------------------------
void setup_Lagrange_interpolation_coefficients_O2_CODE(PDE_Stepper_Data &PDE_D, 
                                                       vector<double> &xarr, 
                                                       bool setup_Integral)
{
    int np=xarr.size();
    
    print_mess_PDE("\n setup_Lagrange_interpolation_coefficients_O2_CODE:: setting up grid ");
    
    PDE_D.LG.dli_dx.clear();
    PDE_D.LG.d2li_d2x.clear();
    PDE_D.dyi.resize(np+3, 0.0);
    PDE_D.ddyi.resize(np+3, 0.0);
    PDE_D.Hi.resize(np+3, 0.0);
    PDE_D.Gi.resize(np+3, 0.0);
    PDE_D.Hip.resize(np+3, 0.0);
    PDE_D.Gip.resize(np+3, 0.0);
    //
    PDE_D.Hi_ptr=&PDE_D.Hi;
    PDE_D.Gi_ptr=&PDE_D.Gi;
    PDE_D.Hip_ptr=&PDE_D.Hip;
    PDE_D.Gip_ptr=&PDE_D.Gip;

    //====================================================================
    // allocate memory
    //====================================================================
    vector<double> dum(6, 0.0);
    
    for(int k=0; k<np; k++)
    {
        PDE_D.LG.dli_dx  .push_back(dum);
        PDE_D.LG.d2li_d2x.push_back(dum);
    }
    
    //====================================================================
    // lower boundary
    //====================================================================
    int i=0;
    for(int j=0; j<5; j++)
    {
        PDE_D.LG.dli_dx  [i][j]=dli_dx_xj  (j, 0, &xarr[0], 5);
        PDE_D.LG.d2li_d2x[i][j]=d2li_d2x_xj(j, 0, &xarr[0], 5);
    }
    
    //====================================================================
    // 1. interior point (lower boundary not needed)
    //====================================================================
    i=1;
    for(int j=0; j<5; j++)
    {
        PDE_D.LG.dli_dx  [i][j]=dli_dx_xj  (j, 1, &xarr[i-1], 5);
        PDE_D.LG.d2li_d2x[i][j]=d2li_d2x_xj(j, 1, &xarr[i-1], 5);
    }
    
    //====================================================================
    // last interior point
    //====================================================================
    i=np-2;
    for(int j=0; j<5; j++)
    {
        PDE_D.LG.dli_dx  [i][j]=dli_dx_xj  (j, 3, &xarr[i-3], 5);
        PDE_D.LG.d2li_d2x[i][j]=d2li_d2x_xj(j, 3, &xarr[i-3], 5);
    }
    
    //====================================================================
    // upper boundary
    //====================================================================
    i=np-1;
    for(int j=0; j<5; j++)
    {
        PDE_D.LG.dli_dx  [i][j]=dli_dx_xj  (j, 4, &xarr[i-4], 5);
        PDE_D.LG.d2li_d2x[i][j]=d2li_d2x_xj(j, 4, &xarr[i-4], 5);
    }
    
    //====================================================================
    // the rest
    //====================================================================
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(i=2; i<np-2; i++)
        for(int j=0; j<5; j++)
        {
            PDE_D.LG.dli_dx  [i][j]=dli_dx_xj  (j, 2, &xarr[i-2], 5);
            PDE_D.LG.d2li_d2x[i][j]=d2li_d2x_xj(j, 2, &xarr[i-2], 5);
        }
    
    //====================================================================
    // integral over frequency domain using 5-point formula
    //====================================================================
    if(setup_Integral)
    {
        PDE_D.LG.Int_li.resize(np, 0.0);
        
        //================================================================
        // internal points with 5-point formula
        //================================================================
        int ix;
        for(ix=1; ix<np-5; ix+=4) 
            for(int j=0; j<5; j++) PDE_D.LG.Int_li[ix+j]+=Int_li_xj(j, &xarr[ix]);

        //================================================================
        // finish off using simple trapeziodal rule
        //================================================================
        for(ix-=4; ix<np-2; ix++) PDE_D.LG.Int_li[ix]+=0.5*(xarr[ix+1]-xarr[ix]);
    }

    print_mess_PDE(" setup_Lagrange_interpolation_coefficients_O2_CODE:: done ");
    
    return;
}



//-------------------------------------------------------------------------------------
// compute dy and ddy
//-------------------------------------------------------------------------------------
void compute_derivatives(const vector<double> &yi, 
                         const vector<double> &Ai, const vector<double> &Bi, 
                         PDE_Stepper_Data &PDE_D, 
                         vector<double> &dy, vector<double> &ddy)
{       
    //==========================================================================
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int i=2; i<(int)yi.size()-2; i++)
    {
        dy[i]=ddy[i]=0.0;
        
        double dum1=PDE_D.LG.d2li_d2x[i][0]*Ai[i]*yi[i+0-2]
                   +PDE_D.LG.d2li_d2x[i][4]*Ai[i]*yi[i+4-2];
        double dum2=PDE_D.LG.d2li_d2x[i][1]*Ai[i]*yi[i+1-2]
                   +PDE_D.LG.d2li_d2x[i][3]*Ai[i]*yi[i+3-2];
        
        ddy[i]=PDE_D.LG.d2li_d2x[i][2]*Ai[i]*yi[i];
        ddy[i]+=dum2;
        ddy[i]+=dum1;

        dum1=PDE_D.LG.dli_dx[i][0]*Bi[i]*yi[i+0-2]
            +PDE_D.LG.dli_dx[i][4]*Bi[i]*yi[i+4-2];
        dum2=PDE_D.LG.dli_dx[i][1]*Bi[i]*yi[i+1-2]
            +PDE_D.LG.dli_dx[i][3]*Bi[i]*yi[i+3-2];
        
        dy[i]=PDE_D.LG.dli_dx[i][2]*Bi[i]*yi[i];
        dy[i]+=dum2;
        dy[i]+=dum1;
        
    }
    
    //==========================================================================
    // derivatives close to lower boundaries
    //==========================================================================
    dy[0]=ddy[0]=0.0;
    dy[1]=ddy[1]=0.0;
    
    //==========================================================================
    for(int j=0; j<5; j++)
    {
        ddy[0]+=PDE_D.LG.d2li_d2x[0][j]*Ai[1]*yi[j];
        dy [0]+=PDE_D.LG.dli_dx  [0][j]*Bi[1]*yi[j];

        ddy[1]+=PDE_D.LG.d2li_d2x[1][j]*Ai[1]*yi[j];
        dy [1]+=PDE_D.LG.dli_dx  [1][j]*Bi[1]*yi[j];
    }

    //==========================================================================
    // derivatives close to upper boundaries
    //==========================================================================
    int i=yi.size()-1;
    dy[i]=ddy[i]=0.0;
    dy[i-1]=ddy[i-1]=0.0;
    
    //==========================================================================
    for(int j=0; j<5; j++)
    {
        ddy[i-1]+=PDE_D.LG.d2li_d2x[i-1][j]*Ai[i-1]*yi[i-4+j];
        dy [i-1]+=PDE_D.LG.dli_dx  [i-1][j]*Bi[i-1]*yi[i-4+j];
        ddy[i]  +=PDE_D.LG.d2li_d2x[i][j]*Ai[i]*yi[i-4+j];
        dy [i]  +=PDE_D.LG.dli_dx  [i][j]*Bi[i]*yi[i-4+j];
    }
    //==========================================================================

    return;
}

//==================================================================================================
//==================================================================================================
