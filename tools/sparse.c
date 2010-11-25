/****************************************/
/* Sparse Matrix algorithms for CLASS   */
/* 15/11 2010                           */
/* Thomas Tram                          */
/****************************************/
/*	This module is used for solving sparse linear systems, arising
	when doing Newton iterations in evolver_ndf15 with a sparse Jacobian.
	The LU factorization is a left-looking algorithm based on the algorithms
	presented in "Direct Methods for Sparse Linear Systems", ISBN 978-0-898716-13-9.
	The primary modification is the ability to refactor a matrix, based on an
	earlier factorization. The routine column_grouping calculates the 'first fit'
	column grouping of a sparse matrix, which is used for evaluating the Jacobian
	with far fewer function evaluations.*/

#include "common.h"
#include "sparse.h"
int sp_mat_alloc(sp_mat** A, int ncols, int nrows, int maxnz, ErrorMsg error_message){
	int ncp =  ncols+1;
	class_alloc((*A),sizeof(sp_mat),error_message);
	class_alloc((*A)->Ax,maxnz*sizeof(double),error_message);
	class_alloc((*A)->Ai,maxnz*sizeof(int),error_message);
	class_alloc((*A)->Ap,(ncp*sizeof(int)),error_message);
	(*A)->ncols = ncols;
	(*A)->nrows = nrows;
	(*A)->maxnz = maxnz;
	return _SUCCESS_;
}

int sp_mat_free(sp_mat *A){
	free(A->Ax);
	free(A->Ai);
	free(A->Ap);
	free(A);
	return _SUCCESS_;
}

int sp_num_alloc(sp_num** N, int n, ErrorMsg error_message){
	int maxnz, k;
	class_alloc((*N),sizeof(sp_num),error_message);
	maxnz = n*(n+1);
	maxnz /=2;
	(*N)->n = n;
	class_call(sp_mat_alloc(&((*N)->L), n, n, maxnz, error_message),
		error_message,error_message);
	class_call(sp_mat_alloc(&((*N)->U), n, n, maxnz, error_message),
		error_message,error_message);
	class_alloc((*N)->xi,n*sizeof(int*),error_message); 
	/* I really want xi to be a vector of pointers to vectors. */
	class_alloc((*N)->xi[0],n*n*sizeof(int),error_message);
	for (k=1;k<n;k++)	(*N)->xi[k] = (*N)->xi[k-1]+n; 
	/*Assign pointers to rows.*/
	class_alloc((*N)->topvec,n*sizeof(int),error_message);
	class_alloc((*N)->pinv,n*sizeof(int),error_message);
	class_alloc((*N)->p,n*sizeof(int),error_message);
	class_alloc((*N)->w,n*sizeof(double),error_message); 
	return _SUCCESS_;
}

int sp_num_free(sp_num *N){
	sp_mat_free(N->L);
	sp_mat_free(N->U);
	free(N->xi[0]);
	free(N->xi);
	free(N->topvec);
	free(N->pinv);
	free(N->p);
	free(N->w);
	free(N);
	return _SUCCESS_;
}

int reachr(sp_mat *G, sp_mat *B,int k, int *xik,int *pinv){
	int p, n, top, *Bp, *Bi, *Gp;
	n=G->ncols; Bp = B->Ap; Bi = B->Ai;Gp = G->Ap;
	top = n;
	for (p=Bp[k];p<Bp[k+1];p++){ /* For each entry in the k'th column of B */
		if (!SPMARKED(Gp,Bi[p])){ /* If node is not marked... */
			dfsr(Bi[p],G,&top,xik,pinv); /* ...start a depth first search at this entry.*/
		}
	}
	for (p=top; p<n; p++) SPMARK(Gp, xik[p]);
	return top;
}

void dfsr(int j, sp_mat *G, int *top, int *xik, int *pinv){
	int i, p, p1, p2, jnew, *Gp = G->Ap, *Gi=G->Ai;
	jnew = pinv[j];
	SPMARK(Gp,j);
	if (jnew>=0){	/*We should consider the jnew column.*/
		p1 = SPUNFLIP(Gp[jnew]); /*Get true column pointers of neighbours:*/
		p2 = SPUNFLIP(Gp[jnew+1]);
		for(p=p1; p<p2; p++){ /* Iterate over neighboors */
			i = Gi[p];
			if (!SPMARKED(Gp,i)){ /* If any unmarked neighboors are found... */
				dfsr(i,G,top,xik,pinv); /*... do depth first search from that node. */
			}
		}
	}
	xik[--(*top)] = j; /*Put column value on stack. */
}

int sp_splsolve(sp_mat *G, sp_mat *B, int k, int*xik, int top, double *x, int *pinv){
	int j, J, p, q, px, n, *Gp, *Gi, *Bp, *Bi;
	double *Gx, *Bx;
	Gp = G->Ap; Gi = G->Ai; Gx = G->Ax;
	Bp = B->Ap; Bi = B->Ai; Bx = B->Ax;
	n = G->ncols;
	
	for (p=top; p<n; p++) x[xik[p]] = 0;
	for (p=Bp[k];p<Bp[k+1];p++) x[Bi[p]] = Bx[p];
	for (px=top; px<n; px++){
		j=xik[px];
		J = pinv[j];
		if (J<0) continue;
		x[j] /= Gx[Gp[J]];
		p = Gp[J]+1;
		q = Gp[J+1];
		for( ;p<q; p++){
			x[Gi[p]] -=Gx[p] * x[j];
		}
	}
	return _SUCCESS_;
}

int sp_ludcmp(sp_num *N, sp_mat *A, double pivtol){
	double pivot, *Lx, *Ux, *x, a, t;
	int *Lp, *Li, *Up, *Ui, *pinv, *pvec, n, ipiv, k, top, p, i, col, lnz, unz;
	n = A->ncols;
	Li = N->L->Ai; Lp = N->L->Ap; Lx = N->L->Ax;
	Ui = N->U->Ai; Up = N->U->Ap; Ux = N->U->Ax;
	lnz = 0; unz = 0;
	x = N->w; pinv = N->pinv; pvec = N->p;
	for (i=0; i<n; i++) x[i]=0;
	for (i=0; i<n; i++) pinv[i] = -1;
	for (k=0; k<=n; k++) Lp[k] = 0;
	
	for(k=0; k<n; k++){
		/* Triangular solve: */
		Lp[k] = lnz;
		Up[k] = unz;
		col = k;
		
		top = reachr(N->L, A, col, N->xi[k], pinv);
		N->topvec[k] = top;
		sp_splsolve(N->L, A, col, N->xi[k], top, x, pinv);
		/* Find pivot: */
		ipiv = -1;
		a = -1;
		for(p=top; p<n; p++){
			i = N->xi[k][p];
			if (pinv[i]<0){
				t = fabs(x[i]);
				if (t>a){
					a = t;
					ipiv = i;
				}
			}
			else{
				Ui[unz] = pinv[i];
				Ux[unz] = x[i];
				unz++;
			}
		}
		if ((ipiv == -1)||(a<=0)) return _FAILURE_;
		if ((pinv[col]<0) && (fabs(x[col])>=a*pivtol)) ipiv = col;
		/* Divide by pivot: */
		pivot = x[ipiv];
		Ui[unz] = k;
		Ux[unz] = pivot;
		unz++;
		pinv[ipiv] = k;
		pvec[k] = ipiv;
		Li[lnz] = ipiv;
		Lx[lnz] = 1.0;
		lnz++;
		for (p=top; p<n; p++){
			i = N->xi[k][p];
			if (pinv[i]<0){
				Li[lnz] = i;
				Lx[lnz] = x[i]/pivot;
				lnz++;
			}
			x[i] = 0;
		}
	}
	/* Finalize: */
	Lp[n] = lnz;
	Up[n] = unz;
	for(p=0; p<lnz; p++) Li[p] = pinv[Li[p]];
	return _SUCCESS_;
}

int sp_lusolve(sp_num *N, double *b, double *x){
	int p, j, n, *Ap, *Ai;
	double *Ax;
	n=N->n;
	/* permute b and initialize x:*/
	for (j=0; j<n; j++) x[N->pinv[j]] = b[j];
	/* lower solve: */
	Ap = N->L->Ap; Ai = N->L->Ai; Ax = N->L->Ax;
	for (j=0; j<n; j++){
		x[j] /=Ax[Ap[j]];
		for (p=Ap[j]+1; p<Ap[j+1]; p++){
			x[Ai[p]] -=Ax[p]*x[j];
		}
	}
	/* upper solve: */
	Ap = N->U->Ap; Ai = N->U->Ai; Ax = N->U->Ax;
	for (j=n-1; j>=0; j--){
		x[j] /=Ax[Ap[j+1]-1];
		for (p=Ap[j];p<Ap[j+1]-1; p++){
			x[Ai[p]] -= Ax[p]*x[j];
		}
	}
	return _SUCCESS_;
}

int sp_refactor(sp_num *N, sp_mat *A){
	double pivot, *Lx, *Ux, *x;
	int *Lp, *Li, *Up, *Ui, *pinv, *pvec, n, ipiv, k, top, p, i, col, lnz, unz;
	n = A->ncols;
	Li = N->L->Ai; Lp = N->L->Ap; Lx = N->L->Ax;
	Ui = N->U->Ai; Up = N->U->Ap; Ux = N->U->Ax;
	lnz = 0; unz = 0;
	x = N->w; pinv = N->pinv; pvec = N->p;
	for (i=0; i<n; i++) x[i]=0;
	for (k=0; k<=n; k++) Lp[k] = 0;
	for(k=0; k<n; k++){
		/* Triangular solve: */
		Lp[k] = lnz;
		Up[k] = unz;
		col = k;
		
		top = N->topvec[k];
		sp_splsolve(N->L, A, col, N->xi[k], top, x, pinv);
		/* Assign values to U and L: */
		ipiv = pvec[k];
		pivot = x[ipiv];
		Li[lnz] = ipiv;
		Lx[lnz] = 1;
		lnz++;
		for (p=top; p<n; p++){
			i = N->xi[k][p];
			if (pinv[i]<k){
				Ui[unz] = pinv[i];
				Ux[unz] = x[i];
				unz++;
			}
			if (pinv[i]>k){
				Li[lnz] = i;
				Lx[lnz] = x[i]/pivot;
				lnz++;
			}
			x[i] = 0;
		}
		Ui[unz] = k;
		Ux[unz] = pivot;
		unz++;
	}
	Lp[n] = lnz;
	Up[n] = unz;
	for(p=0; p<lnz; p++) Li[p] = pinv[Li[p]];
	return _SUCCESS_;
}

int column_grouping(sp_mat *G, int *col_g, int *filled){
	int curcol,testcol,groupnum,fitted,neq;
	int i, *Ap, *Ai;
	
	neq = G->ncols;	Ai = G->Ai; Ap = G->Ap;
	for(i=0;i<neq;i++) col_g[i]=-1;
  
	groupnum=-1;
	for(curcol=0;curcol<neq;curcol++){
		/*Loop through columns..*/
		if (col_g[curcol]==-1){
			/* If current column is not in a group,
			go to next(first) group and assign groupnum to current column: */
			groupnum++;
			col_g[curcol] = groupnum;
			/* put fillness vector equal to current column */
			for(i=0;i<neq;i++) filled[i]=0;
			for(i=Ap[curcol];i<Ap[curcol+1];i++) filled[Ai[i]] = 1; /* A bit convoluted, but it should do it..*/
			/* Try to fit any of the ungrouped columns into this group: */
			for (testcol=(curcol+1);testcol<neq;testcol++){
				if (col_g[testcol]==-1){
					/* loop over row numbers in testcol. If a row number is already filled, break and go to next column.*/
					fitted = 1;
					for(i=Ap[testcol];i<Ap[testcol+1];i++){
						if (filled[Ai[i]] == 1){
							/* We have hit an existing entry */
							fitted=0;
							break;
						}
					}
					if (fitted){
						/* Membership accepted..*/
						col_g[testcol] = groupnum;
						/* Add to filled.. */
						for(i=Ap[testcol];i<Ap[testcol+1];i++) filled[Ai[i]]=1;
					}
				}
			}
		}
	}
	return groupnum;
}
