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
	/* Has to be n+1 because sp_amd uses it for storage:*/
	class_alloc((*N)->q,(n+1)*sizeof(int),error_message);
	class_alloc((*N)->w,n*sizeof(double),error_message);
	class_alloc((*N)->wamd,(8*(n+1))*sizeof(int),error_message);
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
	free(N->q);
	free(N->w);
	free(N->wamd);
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
	int *Lp, *Li, *Up, *Ui, *pinv, *pvec, *q;
	int n, ipiv, k, top, p, i, col, lnz, unz;
	n = A->ncols; q = N->q;
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
		col = q ? (q[k]) : k;

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
	double *Ax, *w;
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
	if (N->q!=NULL){
		/* We must permute once more..*/
		w = N->w;
		for(j=0;j<n;j++) w[j] = x[j];
		for(j=0; j<n; j++) x[N->q[j]] = w[j];
	}
	return _SUCCESS_;
}

int sp_refactor(sp_num *N, sp_mat *A){
	double pivot, *Lx, *Ux, *x;
	int *Lp, *Li, *Up, *Ui, *pinv, *pvec, *q;
	int n, ipiv, k, top, p, i, col, lnz, unz;
	n = A->ncols;
	Li = N->L->Ai; Lp = N->L->Ap; Lx = N->L->Ax;
	Ui = N->U->Ai; Up = N->U->Ap; Ux = N->U->Ax;
	q = N->q;
	lnz = 0; unz = 0;
	x = N->w; pinv = N->pinv; pvec = N->p;
	for (i=0; i<n; i++) x[i]=0;
	for (k=0; k<=n; k++) Lp[k] = 0;
	for(k=0; k<n; k++){
		/* Triangular solve: */
		Lp[k] = lnz;
		Up[k] = unz;
		col = q ? (q[k]) : k;

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
  int testcol,groupnum,fitted,neq;
  int i, *Ap, *Ai,done;

  neq = G->ncols;
  Ai = G->Ai; Ap = G->Ap;
  for(i=0;i<neq;i++)
    col_g[i]=-1;

  for(groupnum=0; groupnum<neq; groupnum++){
    //Reset filled vector:
    for (i=0; i<neq; i++){
      filled[i] = 0;
    }
    //Try to assign remaining columns to current group:
    done = _TRUE_;
    for (testcol=0; testcol<neq; testcol++){
      if (col_g[testcol]!=-1){
	continue;
      }
      done = _FALSE_;
      //Is current testcol in conflict with the groups filling vector?
      fitted = _TRUE_;
      for (i=Ap[testcol]; i<Ap[testcol+1]; i++){
	if (filled[Ai[i]] != 0){
	  fitted = _FALSE_;
	  break;
	}
      }
      if (fitted == _TRUE_){
	//Test succesfull
	col_g[testcol] = groupnum;
	for (i=Ap[testcol]; i<Ap[testcol+1]; i++){
	  filled[Ai[i]] = 1;
	}
      }
    }
    if (done == _TRUE_){
      //No columns remaining.
      //No column belongs to current groupnum.
      break;
    }
  }
  return (groupnum-1);
}

int sp_amd(int *Cp, int *Ci, int n, int nzmax, int *P, int *W){
	int *last, *len, *nv, *next, *head, *elen, *degree, *w, *hhead;
	int d, dk, dext, lemax=0, e, elenk, eln, i, j, k, k1, k2, k3, jlast, ln;
	int dense, mindeg=0, nvi, nvj, nvk, mark, wnvi, ok, nel=0;
	int p, p1, p2, p3, p4, pj, pk, pk1, pk2, pn, q, cnz;
	unsigned int h;
	/*	I assume that the sparse matrix C is symmetrix (C = A + A' in our case)
		and that the diagonal elements has been removed. C must be large enough,
		C->max_nonzero >= (6/5)*(C->Ap[n]) + 2n. Work array W must be 8*(n+1).
	*/
	dense = MAX(16,10*sqrt((double) n));
	dense = MIN(n-2,dense);
	cnz = Cp[n];
	/*	Assign pointers to positions in work array:*/
	len = W;
	nv = W+(n+1);
	next = W+2*(n+1);
	head = W+3*(n+1);
	elen = W+4*(n+1);
	degree = W+5*(n+1);
	w = W+6*(n+1);
	hhead = W+7*(n+1);
	last = P;
	/* Initialise quotient graph: */
	for(k=0; k<n; k++) len[k] = Cp[k+1]-Cp[k];
	len[n] = 0;
	for(i=0; i<=n; i++){
		head[i] = -1;
		last[i] = -1;
		next[i] = -1;
		hhead[i] = -1;
		nv[i] = 1;
		w[i] = 1;
		elen[i] = 0;
		degree[i] = len[i];
	}
	mark = sp_wclear(0, 0, w, n);
	elen[n] = -2;
	Cp[n] = -1;
	w[n] = 0;
	/* Initialize degree lists */
	for(i=0; i<n; i++){
		d = degree[i];
		if(d==0){
			elen[i] = -2;
			nel++;
			Cp[i] = -1;
			w[i] = 0;
		}
		else if(d>dense){
			nv[i] = 0;
			elen[i] = -1;
			nel++;
			Cp[i] = SPFLIP(n);
			nv[n]++;
		}
		else{
			if(head[d]!=-1) last[head[d]] = i;
			next[i] = head[d];
			head[d] = i;
		}
	}
	while(nel<n){
		/* Select node ofminimum degree: */
		for(k=-1; (mindeg<n)&&((k=head[mindeg])==-1); mindeg++);
		if(next[k]!=-1) last[next[k]]=-1;
		head[mindeg] = next[k];
		elenk = elen[k];
		nvk = nv[k];
		nel +=nvk;
		/* Garbage collection */
		if((elenk>0)&&(cnz+mindeg >=nzmax)){
			for(j=0; j<n; j++){
				if((p=Cp[j])>=0){
					Cp[j] = Ci[p];
					Ci[p] = SPFLIP(j);
				}
			}
			for(q=0,p=0; p<cnz; ){
				if((j=SPFLIP(Ci[p++]))>=0){
					Ci[q] = Cp[j];
					Cp[j] = q++;
					for(k3=0; k3<len[j]-1; k3++) Ci[q++] = Ci[p++];
				}
			}
			cnz = q;
		}
		/* Construct new element: */
		dk = 0;
		nv[k] = -nvk;
		p = Cp[k];
		pk1 = (elenk==0) ? p : cnz;
		pk2 = pk1;
		for(k1=1; k1<=elenk + 1; k1++){
			if (k1>elenk){
				e=k;
				pj = p;
				ln = len[k] - elenk;
			}
			else{
				e = Ci[p++];
				pj = Cp[e];
				ln = len[e];
			}
			for(k2=1; k2<=ln; k2++){
				i=Ci[pj++];
				if((nvi=nv[i])<=0) continue;
				dk += nvi;
				nv[i] = -nvi;
				Ci[pk2++] = i;
				if (next[i]!=-1) last[next[i]] =last[i];
				if (last[i]!=-1){
					next[last[i]] = next[i];
				}
				else{
					head[degree[i]] = next[i];
				}
			}
			if(e!=k){
				Cp[e] = SPFLIP(k);
				w[e] = 0;
			}
		}
		if (elenk!=0) cnz = pk2;
		degree[k] = dk;
		Cp[k] = pk1;
		len[k] = pk2 - pk1;
		elen[k] = -2;
		/* Find set differences: */
		mark = sp_wclear(mark, lemax, w, n);
		for(pk=pk1; pk<pk2; pk++){
			i = Ci[pk];
			if((eln = elen[i]) <=0) continue;
			nvi = -nv[i];
			wnvi = mark - nvi;
			for(p=Cp[i]; p<=Cp[i] + eln -1; p++){
				e = Ci[p];
				if(w[e] >=mark){
					w[e] -= nvi;
				}
				else if(w[e]!=0){
					w[e] = degree[e] + wnvi;
				}
			}
		}
		/* Degree update:*/
		for(pk=pk1; pk<pk2; pk++){
			i=Ci[pk];
			p1 = Cp[i];
			p2 = p1 + elen[i] - 1;
			pn = p1;
			for(h=0, d=0, p=p1; p<=p2; p++){
				e = Ci[p];
				if(w[e]!=0){
					dext = w[e] - mark;
					if (dext>0){
						d+=dext;
						Ci[pn++] = e;
						h +=e;
					}
					else{
						Cp[e] = SPFLIP(k);
						w[e] = 0;
					}
				}
			}
			elen[i] = pn-p1+1;
			p3 = pn;
			p4 = p1 + len[i];
			for(p=p2+1; p<p4; p++){
				j=Ci[p];
				if((nvj = nv[j]) <= 0) continue;
				d += nvj;
				Ci[pn++] = j;
				h += j;
			}
			if (d==0){
				Cp[i] = SPFLIP(k);
				nvi = -nv[i];
				dk -= nvi;
				nvk += nvi;
				nel +=nvi;
				nv[i] = 0;
				elen[i] = -1;
			}
			else{
				degree[i] = MIN(degree[i],d);
				Ci[pn] = Ci[p3];
				Ci[p3] = Ci[p1];
				Ci[p1] = k;
				len[i] = pn -p1 +1;
				h %= n;
				next[i] = hhead[h];
				hhead[h] = i;
				last[i] = h;
			}
		}
		degree[k] = dk;
		lemax = MAX(lemax,dk);
		mark = sp_wclear(mark+lemax, lemax, w, n);
		/* Supernode detection: */
		for(pk=pk1; pk <pk2; pk++){
			i = Ci[pk];
			if(nv[i]>=0) continue;
			h = last[i];
			i = hhead[h];
			hhead[h] = -1;
			for( ;(i!=-1)&&(next[i]!=-1);i = next[i],mark++){
				ln = len[i];
				eln = elen[i];
				for(p=Cp[i] + 1; p<=Cp[i]+ln-1; p++) w[Ci[p]] = mark;
				jlast = i;
				for(j=next[i]; j!=-1; ){
					ok = (len[j]==ln)&&(elen[j] == eln);
					for(p=Cp[j] +1;ok&&p<=Cp[j]+ln-1; p++){
						if(w[Ci[p]]!=mark) ok=0;
					}
					if (ok){
						Cp[j] = SPFLIP(i);
						nv[i] +=nv[j];
						nv[j] = 0;
						elen[j] = -1;
						j = next[j];
						next[jlast] = j;
					}
					else{
						jlast =j;
						j = next[j];
					}
				}
			}
		}
		/* Finalize new element: */
		for(p=pk1, pk=pk1; pk<pk2; pk++){
			i = Ci[pk];
			if ((nvi = -nv[i]) <=0) continue;
			nv[i] = nvi;
			d = degree[i] +dk -nvi;
			d = MIN(d, (n-nel-nvi));
			if(head[d]!=-1) last[head[d]] = i;
			next[i] = head[d];
			last[i] = -1;
			head[d] = i;
			mindeg = MIN(mindeg,d);
			degree[i] = d;
			Ci[p++] = i;
		}
		nv[k] = nvk;
		if((len[k] = p-pk1)==0){
			Cp[k] = -1;
			w[k] = 0;
		}
		if(elenk!=0) cnz = p;
	}
	/* Postordering */
	for(i=0; i<n; i++) Cp[i] = SPFLIP(Cp[i]);
	for(j=0; j<=n; j++) head[j] = -1;
	for(j=n; j>=0; j--){
		if(nv[j] > 0) continue;
		next[j] = head[Cp[j]];
		head[Cp[j]] = j;
	}
	for(e=n; e>=0; e--){
		if(nv[e]<=0) continue;
		if(Cp[e] !=-1){
			next[e] = head[Cp[e]];
			head[Cp[e]] = e;
		}
	}
	for(k=0,i=0; i<=n; i++){
		if (Cp[i] == -1) k = sp_tdfs(i, k, head, next, P, w);
	}
	return _SUCCESS_;
}

int sp_wclear(int mark, int lemax, int *w, int n){
	int k;
	if (mark<2 || (mark+lemax<0)){
		for(k=0; k<n; k++) if(w[k]!=0) w[k] = 1;
		mark = 2;
	}
	return (mark);
}

int sp_tdfs(int j, int k, int *head, const int *next, int *post, int *stack){
	int i, p, top =0;
	stack[0] = j;
	while(top>=0){
		p = stack[top];
		i = head[p];
		if (i==-1){
			top--;
			post[k++] = p;
		}
		else{
			head[p] = next[i];
			stack[++top] = i;
		}
	}
	return (k);
}


