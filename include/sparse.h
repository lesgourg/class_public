#ifndef __SPA__
#define __SPA__
/****************************************/
/* Sparse Matrix algorithms for CLASS   */
/* 15/11 2010                           */
/* Thomas Tram                          */
/****************************************/
#include "common.h"

/* Structures: */
typedef struct sparse_matrix{
	/* Sparse matrix in compressed column form: */
	int ncols;		/* Number of columns */
	int nrows;		/* Number of rows */
	int maxnz;		/* Maximum number of non-zero entries*/
	int *Ap;		/* Ap[0..ncols]. Ap[k+1]-Ap[k] is the number of entries in the k'th column. */
	int *Ai;		/* Ai[0..(maxnz-1)]. Contains the row indices of the entries. */
	double *Ax;		/* Ax[0..(maxnz-1)]. Contains the values of the entries. */
} sp_mat;

typedef struct sparse_numerical{
	/* Sparse LU decomposition along with enough information to do a fast refactorization: */
	int n;			/*Matrix assumed square, [nxn] */
	sp_mat *L;		/*L and U is the factors of the decomposed matrix.*/
	sp_mat *U;
	int **xi;		/*xi[k] points to a row of xi, which holds the topological ordered indices.*/
	int *topvec;	/*topvec[k] holds the first index in xi[k].*/
	int *pinv;		/*Inverse row permutation. */
	int *p;			/*Row permutation. */
	int *q;			/* Column permutation */
	int *wamd;		/* Work array for sp_amd */
	double *w;		/* Work array for sp_lu */
} sp_num;


/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif
/* Routines and macros: */
int sp_mat_alloc(sp_mat** A, int ncols, int nrows, int maxnz, ErrorMsg error_message);
int sp_mat_free(sp_mat *A);
int sp_num_alloc(sp_num** N, int n,ErrorMsg error_message);
int sp_num_free(sp_num *N);
int reachr(sp_mat *G, sp_mat *B,int k, int *xik,int *pinv);
void dfsr(int j, sp_mat *G, int *top, int *xik, int *pinv);
int sp_splsolve(sp_mat *G, sp_mat *B, int k, int*xik, int top, double *x, int *pinv);
int sp_ludcmp(sp_num *N, sp_mat *A, double pivtol);
int sp_lusolve(sp_num *N, double *b, double *x);
int sp_refactor(sp_num *N, sp_mat *A);
int column_grouping(sp_mat *G, int *col_g, int *col_wi);
int sp_amd(int *Cp, int *Ci, int n, int cnzmax, int *P, int *W);
int sp_wclear(int mark, int lemax, int *w, int n);
int sp_tdfs(int j, int k, int *head, const int *next, int *post, int *stack);


#define SPFLIP(i) (-(i)-2)
#define SPUNFLIP(i) (((i)<0) ? SPFLIP(i) : (i))
#define SPMARKED(w,j) (w[j] < 0)
#define SPMARK(w,j) {w[j] = SPFLIP(w[j]);}

#ifdef __cplusplus
}
#endif


#endif
