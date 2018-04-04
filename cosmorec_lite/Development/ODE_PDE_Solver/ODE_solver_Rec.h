//==================================================================================================
// ODE_solver routines 
//
// Purpose: solve a system of stiff ODEs with moderate dimension (n<~100). 
// 
// Basic Aspects: The solver is based on a 6th order Gears method (implicit && stiffly-stable) with 
// variable time-step. The linear algebra system is solved iteratively using a biconjugate gradient 
// method. A guess for the next solution is extrapolated from the previous solutions; 
//
// Author: Jens Chluba with contributions from Geoff Vasil at CITA
//
// First implementation: 12/01/2010
// Last modification   : 03/01/2012
//==================================================================================================
// 24/08/2012: improved some of the functions and tidied the code up
//
// 15/06/2011: set max order to 5; This was making the code a bit faster; 
//             tidied the code a bit;
//==================================================================================================

#ifndef _ODE_SOLVER_REC_H_
#define _ODE_SOLVER_REC_H_

#include <vector>
#include <string>

using namespace std;

namespace ODE_solver_Rec
{
    //==============================================================================================
    // Structure that contains the solution y and its derivative dy at redshift z
    //==============================================================================================
    struct ODE_solver_Solution
    {
        double z;
        vector<double> y;
        vector<double> dy;
    };
    
    //==============================================================================================
    struct ODE_solver_accuracies
    {
        vector<double> rel;
        vector<double> abs;
    };
    
    //==============================================================================================
    struct ODE_solver_matrix
    {
        vector<double> A;
        vector<int> row;
        vector<int> col;
        vector<int> diags;
    };

    //==============================================================================================
    struct ODE_Solver_data
    {
        void (*fcn_col_ptr)(int *neq, double *z, double *y, double *f, int col);
        void (*fcn_col_para_ptr)(int *neq, double *z, double *y, double *f, int col, void *p);
        bool num_jac;
        //
        ODE_solver_Solution Stemp;
        ODE_solver_Solution *Sptr;
        ODE_solver_Solution *Snewptr, *Sguessptr;
        ODE_solver_Solution *Snptr[6];
        //
        ODE_solver_matrix *Jacobian_ptr;
        //
        vector<double> F, dY, abs_vector, rel_vector;
        //
        double Dz_z_last;
        double zstart_solver;
        double tolSol;
        int order;
        int count;
        int Jac_is_set;
        // counting number of significant up and down-steps
        int n_up, n_down;
        double direction;
        // to pass on parameters for the setup
        void *p;
    };
    
    //==============================================================================================
    void allocate_memory(ODE_solver_Solution &Sz, ODE_solver_accuracies &tols, int neq); 
    
    void clear_memory(ODE_solver_Solution &Sz, ODE_solver_accuracies &tols);
    
    void clear_SolverData(ODE_Solver_data &SD);
    
    void clear_all_memory(ODE_solver_Solution &Sz, 
                          ODE_solver_accuracies &tols, 
                          ODE_Solver_data &SD);
    
    void set_verbosity(int verb);
    
    //==============================================================================================
    // ODE_Solver_set_up_solution_and_memory::
    //
    // This routine has to be called before the solver is called. In particular the memory is set 
    // and the function pointer for the ODE evaluation is given. 
    //
    // ODE_solver_Solution Sz.z should be set to the initial redshift;
    // ODE_solver_Solution rtol[] & atol[] should contain the relative error request;
    // ODE_solver_Solution Sz.y[] should contain the initial solution;
    // ODE_solver_Solution Sz.dy[] will be computed using fn_ptr;
    //
    // ODE_solver_accuracies tols.rel[] should contain the relative error request;
    // ODE_solver_accuracies tols.abs[] should contain the relative error request;
    //
    // fn_ptr declares the ODE system according to y'= f(z, y);
    //
    //==============================================================================================
    void ODE_Solver_set_up_solution_and_memory(ODE_solver_Solution &Sz, 
                                               ODE_solver_accuracies &tols, 
                                               ODE_Solver_data &ODE_Solver_info, 
                                               void (*fcn_ptr)(int *neq, double *z, 
                                                               double *y, double *f, 
                                                               int col), 
                                               bool num_jac=1);
    
    //==============================================================================================
    // to reset error requirements
    //==============================================================================================
    void ODE_Solver_set_errors(ODE_solver_accuracies &tols, ODE_Solver_data &ODE_Solver_info);
    
    //==============================================================================================
    // This routine has can only be called after the setup was done 
    // (ODE_Solver_set_up_solution_and_memory)
    //
    // After the cal ODE_solver_Solution Sz will contain the solution at time zend; initially Sz 
    // should contain the solution at zs. Previous solutions will also be stored by the routine.
    //
    //==============================================================================================
    int ODE_Solver_Solve_history(double zs, double zend,
                                 ODE_solver_Solution &Sz, 
                                 ODE_Solver_data &ODE_Solver_info);

}

#endif

//==================================================================================================
//==================================================================================================
