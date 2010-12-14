#ifndef __QSS__
#define __QSS__

/******************************************/
/* Quadrature Sampling Strategy for CLASS */
/* 10/12 2010                             */
/* Thomas Tram                            */
/******************************************/
#include "common.h"

/* Structures for QSS */

typedef struct adaptive_integration_tree_node{
  /* binary tree node: */
  double I;		/* Estimate of integral */
  double err;		/* Estimated error */
  double *x;		/* Pointer to the abscissas of node */
  double *w;		/* Pointer to the coresponding weights */
  int leaf_childs;/* Number of leafes under current node. 1 means that the node is a leaf. */
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
		    int (*test)(void * params_for_function, double q, double *psi),
		    int (*function)(void * params_for_function, double q, double *f0),
		    void * params_for_function,
		    ErrorMsg errmsg);

  int get_leaf_x_and_w(qss_node *node, int *ind, double *x, double *w);
  int reduce_tree(qss_node *node, int level);
  int burn_tree(qss_node *node);
  int leaf_count(qss_node *node);
  double get_integral(qss_node *node, int level);
  qss_node* gk_adapt(
		     int (*test)(void * params_for_function, double q, double *psi),
		     int (*function)(void * params_for_function, double q, double *f0),
		     void * params_for_function,
		     double tol, 
		     int treemode, 
		     double a, 
		     double b);
  int compute_Laguerre(double *x, double *w, int N, double alpha, double *b, double *c);
  int gk_quad(int (*test)(void * params_for_function, double q, double *psi),
	      int (*function)(void * params_for_function, double q, double *f0),
	      void * params_for_function,
	      qss_node* node, 
	      double a, 
	      double b);
  double testfun(double x);




#ifdef __cplusplus
}
#endif


#endif
