#ifndef __QSS__
#define __QSS__

#define _MIN_NUMBER_OF_LAGUERRE_POINTS_ 5

/******************************************/
/* Quadrature Sampling Strategy for CLASS */
/* 10/12 2010                             */
/* Thomas Tram                            */
/******************************************/
#include "common.h"

enum ncdm_quadrature_method {qm_auto, qm_Laguerre, qm_trapz_indefinite, qm_trapz};

/* Structures for QSS */

typedef struct adaptive_integration_tree_node{
  /* binary tree node: */
  double I;		/* Estimate of integral */
  double err;		/* Estimated error */
  double *x;		/* Pointer to the abscissas of node */
  double *w;		/* Pointer to the corresponding weights */
  int leaf_childs;/* Number of leafs under current node. 1 means that the node is a leaf. */
  /* Pointer to children: */
  struct  adaptive_integration_tree_node *left, *right;	/* Pointer to left child. */
} qss_node;

    /**
     * Boilerplate for C++
     */
#ifdef __cplusplus
    extern "C" {
#endif
      int get_qsampling(double *x,
			double *w,
			int *N,
			int N_max, double rtol,
			double *qvec,
			int qsiz,
			int (*test)(void * params_for_function, double q, double *psi),
			int (*function)(void * params_for_function, double q, double *f0),
			void * params_for_function,
			ErrorMsg errmsg);
       int get_qsampling_manual(double *x,
				double *w,
				int N,
				double qmax,
				enum ncdm_quadrature_method method,
				double *qvec,
				int qsiz,
				int (*function)(void * params_for_function, double q, double *f0),
				void * params_for_function,
				ErrorMsg errmsg);

      int sort_x_and_w(double *x, double *w, double *workx, double *workw, int startidx, int endidx);
      int get_leaf_x_and_w(qss_node *node, int *ind, double *x, double *w,int isindefinite);
      int reduce_tree(qss_node *node, int level);
      int burn_tree(qss_node *node);
      int leaf_count(qss_node *node);
      double get_integral(qss_node *node, int level);
      int gk_adapt(
		   qss_node **node,
		   int (*test)(void * params_for_function, double q, double *psi),
		   int (*function)(void * params_for_function, double q, double *f0),
		   void * params_for_function,
		   double tol,
		   int treemode,
		   double a,
		   double b,
		   int isindefinite,
		   ErrorMsg errmsg);
      int compute_Hermite(double *x, double *w, int N, int alpha, double *b, double *c);
      int compute_Laguerre(double *x, double *w, int N, double alpha, double *b, double *c, int totalweight);
      int gk_quad(int (*test)(void * params_for_function, double q, double *psi),
		  int (*function)(void * params_for_function, double q, double *f0),
		  void * params_for_function,
		  qss_node* node,
		  double a,
		  double b,
		  int isindefinite);
      double testfun(double x);

      int quadrature_gauss_legendre(
				    double *mu,
				    double *w8,
				    int n,
				    double tol,
				    ErrorMsg error_message);

      int quadrature_gauss_legendre_2D(
				       int n,
				       double * x,
				       double * y,
				       double * w,
				       ErrorMsg error_message);

      int quadrature_in_rectangle(
				  double xl,
				  double xr,
				  double yl,
				  double yr,
				  int *n,
				  double ** x,
				  double ** y,
				  double ** w,
				  ErrorMsg error_message);

      int cubature_order_eleven(
				double xl,
				double xr,
				double yl,
				double yr,
				double *x,
				double *y,
				double *w,
				ErrorMsg error_message);


#ifdef __cplusplus
    }
#endif


#endif
