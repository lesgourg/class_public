#ifndef __EVO__
#define __EVO__
#include "common.h"
// #include "perturbations.h"
#include "sparse.h"
#define TINY 1e-50
/**************************************************************/

struct jacobian{
/*Stuff for normal method: */
	double **dfdy;
	double *jacvec; /*Stores experience gained from subsequent calls */
	double **LU;
	double *LUw;
	int *luidx;
	/*Sparse stuff:*/
	int use_sparse;
	int sparse_stuff_initialized;
	int max_nonzero;     /*Maximal number of non-zero entries to be considered sparse */
	int repeated_pattern;
	int trust_sparse; /* Number of times a pattern is repeated (actually included) before we trust it. */
	int has_grouping;
	int has_pattern;
	int new_jacobian; /* True if sp_ludcmp has not been run on the current jacobian. */
	int cnzmax;
	int *col_group; /* Column grouping. Groups go from 0 to max_group*/
	int *col_wi; /* Workarray for column grouping*/
	int max_group; /*Number of columngroups -1 */
	sp_mat *spJ; /* Stores the matrix we want to decompose */
	double *xjac; /*Stores the values of the sparse jacobian. (Same pattern as spJ) */
	sp_num *Numerical; /*Stores the LU decomposition.*/
	int *Cp; /* Stores the column pointers of the spJ+spJ' sparsity pattern. */
	int *Ci; /* Stores the row indices of the  spJ+spJ' sparsity pattern. */
};

struct numjac_workspace{
	/* Allocate vectors and matrices: */
	double *yscale;
	double *del;
	double * Difmax;
	double * absFdelRm;
	double * absFvalue;
	double * absFvalueRm;
	double * Fscale;
	double * ffdel;
	double * yydel;
	double * tmp;

	double **ydel_Fdel;

	int * logj;
	int * Rowmax;
};

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int initialize_jacobian(struct jacobian *jac, int neq, ErrorMsg error_message);
  int uninitialize_jacobian(struct jacobian *jac);
  int initialize_numjac_workspace(struct numjac_workspace * nj_ws,int neq, ErrorMsg error_message);
  int uninitialize_numjac_workspace(struct numjac_workspace * nj_ws);
  int calc_C(struct jacobian *jac);
  int interp_from_dif(double tinterp,double tnew,double *ynew,double h,double **dif,int k, double *yinterp,
		      double *ypinterp, double *yppinterp, int* index, int neq, int output);
  int new_linearisation(struct jacobian *jac,double hinvGak,int neq, ErrorMsg error_message);
  int adjust_stepsize(double **dif, double abshdivabshlast, int neq,int k);
  void eqvec(double *datavec,double *emptyvec, int n);
  int lubksb(double **a, int n, int *indx, double b[]);
  int ludcmp(double **a, int n, int *indx, double *d, double *vv);
  int fzero_Newton(int (*func)(double *x,
                               int x_size,
                               void *param,
                               double *F,
                               ErrorMsg error_message),
                   double *x_inout,
                   double *dxdF,
                   int x_size,
                   double tolx,
                   double tolF,
                   void *param,
                   int *fevals,
                   ErrorMsg error_message);

  int numjac(int (*derivs)(double x,double * y,double * dy,void * parameters_and_workspace,ErrorMsg error_message),
	     double t, double *y, double *fval, struct jacobian *jac, struct numjac_workspace *nj_ws,
	     double thresh, int neq, int *nfe,
	     void * parameters_and_workspace_for_derivs, ErrorMsg error_message);


int evolver_ndf15(
	int (*derivs)(double x,double * y,double * dy,
		void * parameters_and_workspace, ErrorMsg error_message),
	double x_ini,
	double x_final,
	double * y_inout,
 	int * used_in_output,
	int neq,
	void * parameters_and_workspace_for_derivs,
	double rtol,
	double minimum_variation,
	int (*timescale_and_approximation)(double x,
					   void * parameters_and_workspace,
					   double * timescales,
					   ErrorMsg error_message),
	double timestep_over_timescale,
	double * t_vec,
	int t_res,
	int (*output)(double x,double y[],double dy[],int index_x,void * parameters_and_workspace,
		ErrorMsg error_message),
	int (*print_variables)(double x, double y[], double dy[], void *parameters_and_workspace,
		ErrorMsg error_message),
	ErrorMsg error_message);


#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
