//======================================================================================
// Author: Jens Chluba 
// first implementation: Jan 2002
// last modification: June 2012
//
// Purpose: collection of several simple routines
//======================================================================================
// Jun 2012: fixed a memory issue with spline setup routines.
// Jan 2012: added simple routines to load tables of data and create splines
// Dec 2011: added routines for Gamma and incomplete Gamma functions

#ifndef ROUTINES_H
#define ROUTINES_H

#include <iostream>
#include <string>
#include <vector>

using namespace std;

//======================================================================================
// Special functions
//======================================================================================
double Dawson_Int(const double &x); // Dawson integral; Based on Antia
double erf_JC(const double &x);     // Error function for a real argument
double erfc_JC(const double &x);    // complementary Error function for a real argument
double Gamma_JC(const double &x);                    // Gamma-function
double Gamma_JC(const double &x, const double &a);   // incomplete Gamma-function

double max(double a, double b);
double min(double a, double b);

double factorial(int n);
double log10factorial(int n);
double log10factorial_full(int n);
double factorial_corrfac(int n);

//======================================================================================
// checking for nan
//======================================================================================
bool isnan_JC(double a);

//======================================================================================
// routines for interpolation; based in GSL
//======================================================================================
int calc_spline_coeffies_JC(int nxi, const double *za, const double *ya, 
                            string variable="");

void update_spline_coeffies_JC(int memindex, int nxi, 
                               const double *za, const double *ya, 
                               string variable="");

double calc_spline_JC(double x, int memindex, string mess="");

void free_spline_JC(int memindex, string mess="");
void free_all_splines_JC();
void show_spline_memory();

//======================================================================================
// load tables of data (added Jan 2012)
//======================================================================================
void load_data_from_file(string fname, int cols, vector<int> &spline_mem_indices, 
                         bool logx, bool logy);

void load_data_from_file_loglin(string fname, int cols, vector<int> &spline_mem_indices);
void load_data_from_file_loglog(string fname, int cols, vector<int> &spline_mem_indices);


//======================================================================================
// npol-1 is degree of the interpolating polynomial
//======================================================================================
void polint_JC(const double *xa, const double *ya, int na, const double x, int npol, 
               double *y, double *dy);

void polint_JC(const double *xa, const double *ya, int na, const double x, int &istart, 
               int npol, double *y, double *dy);

//======================================================================================
// simple grid setup
//======================================================================================
void init_xarr(double x0, double xm, double *xarr, int npts, int method_flag, int mess_flg);
void init_xarr(double x0, double xm, double *xarr, int npts, int method_flag);

//======================================================================================
void wait_f_r();
void wait_f_r(string mess);
void wait_f_r(int num);
void wait_f_r(double num);

//======================================================================================
void locate_JC(const double xx[], unsigned long n, double x, unsigned long *j);
void hunt(const double xx[], unsigned long n, double x, unsigned long *jlo);

//======================================================================================
// i/o modules
//======================================================================================
string int_to_string(int i);
string int_to_string(int i, int ni);

//======================================================================================
// root-finding methods
//======================================================================================
double find_root(double (* func)(double *), double x1, double x2, double xacc);
double find_root_brent(double (* func)(double *), double x1, double x2, double xacc);
double find_root_brent(double (* func)(double *, void *p), void *p, 
                       double x1, double x2, double xacc);

//======================================================================================
// for xmgrace output
//======================================================================================
void plot_xy_function_linear(const vector<double> &xarr, const vector<double> &yarr);
void plot_xy_function_log   (const vector<double> &xarr, const vector<double> &yarr);
void plot_xy_function_loglog(const vector<double> &xarr, const vector<double> &yarr);

//======================================================================================
// Wigner 3J symbol (added 17.05.2011)
//======================================================================================
double Wigner_3J_Symbol(int j1, int j2, int j3, int m1, int m2, int m3);

//======================================================================================
// divide and conquer sum
//======================================================================================
double DC_sum(const double *yarr, int M);
double DC_sum(const vector<double> &yarr);

//======================================================================================
// divide and conquer sum of product
//======================================================================================
double DC_sumprod(const double *yarr, const double *zarr, int M);
double DC_sumprod(const vector<double> &yarr, const vector<double> &zarr);

#endif
