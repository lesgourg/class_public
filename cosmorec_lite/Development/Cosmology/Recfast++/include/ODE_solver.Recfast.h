//=====================================================================================================
// ODE_solver_Recfast routines 
//
// Purpose: solve the stiff ODE for the RECFAST problem. Since the dimension of the problem 
// is very low (i.e. 3 equations), not much effort was made in terms of optimization for the
// Linear-Algebra methods or memory consumption. The solver is 'infinitely' fast compared to 
// routines that might call RECFAST (e.g. CMBfast, CAMB, etc).
// 
// Basic Aspects: The solver is based on a 6th order Gears method (implicit && stiffly-stable) with 
// variable time-step. The linear algebra system is solved iteratively using a biconjugate gradient 
// method. A guess for the next solution is extrapolated from the previous solutions.
//
// A similar solver was used for the work of Chluba, Vasil & Dursi, 2010, MNRAS, 407, pp. 599-612.
//
// Author: Jens Chluba with contributions from Geoff Vasil at CITA
// Date: 12/01/2010
//=====================================================================================================

#ifndef _ODE_SOLVER_RECFAST_H_
#define _ODE_SOLVER_RECFAST_H_

#include <vector>

using namespace std;

//=====================================================================================================
// Structure that contains the solution y and its derivative dy at redshift z
//=====================================================================================================
struct ODE_solver_Solution
{
    double z;
    vector<double> y;
    vector<double> dy;
};

//=====================================================================================================
// This routine has to be called before the solver is called. In particular the memory is set and the 
// function pointer for the ODE evaluation is given. 
//
// ODE_solver_Solution Sz.z should be set to the initial redshift;
// ODE_solver_Solution Sz.y[] should contain the initial solution;
// ODE_solver_Solution Sz.dy[] will be computed using fn_ptr;
// fn_ptr declares the ODE system according to y'= f(z, y)
//
//=====================================================================================================
void ODE_Solver_set_up_solution_and_memory(double z, ODE_solver_Solution &Sz, 
                                           void (*fcn_ptr)(int *neq, double *z, double *y, double *f));

//=====================================================================================================
// This routine has can only be called after the setup was done (ODE_Solver_set_up_solution_and_memory)
//
// After the cal ODE_solver_Solution Sz will contain the solution at time zend; initially Sz should 
// contain the solution at zs. Previous solutions will also be stored by the routine.
//
//=====================================================================================================
int ODE_Solver_Solve_history(double zs, double zend, ODE_solver_Solution &Sz);

#endif
