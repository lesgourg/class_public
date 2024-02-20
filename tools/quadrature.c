/******************************************/
/* Quadrature Sampling Strategy for CLASS */
/* 10/12 2010                             */
/* Thomas Tram                            */
/******************************************/
#include "quadrature.h"

int get_qsampling_manual(double *x,
			 double *w,
			 int N,
			 double qmax,
			 enum ncdm_quadrature_method method,
			 double *qvec,
			 int qsiz,
			 int (*function)(void * params_for_function, double q, double *f0),
			 void * params_for_function,
			 ErrorMsg errmsg) {

  double y, h, t;
  double *b, *c;
  int i;
  switch (method){ 
  case (qm_auto) :
    return _FAILURE_;
  case (qm_Laguerre) :
    /* Allocate storage for Laguerre coefficients: */
    class_alloc(b,N*sizeof(double),errmsg);
    class_alloc(c,N*sizeof(double),errmsg);
    compute_Laguerre(x,w,N,0.0,b,c,_TRUE_);
    for (i=0; i<N; i++){
      (*function)(params_for_function,x[i],&y);
      w[i] *= y;
    }
    free(b);
    free(c);
    return _SUCCESS_;
  case (qm_trapz) :
    for (i=0; i<N; i++){
      /** Note that we count q=0 as an extra point with weight 0 */
      h = qmax/N;
      x[i] = h + i*h;
      (*function)(params_for_function,x[i],&y);
      w[i] = y*h;
      if (i==N-1)
	w[i] *=0.5;
    }
    return _SUCCESS_;
  case (qm_trapz_indefinite) :
    /** We do the variable transformation q = 1/t-1. The trapezoidal rule is closed, but since the distribution function
	goes to zero in both limits, we can use an effectively N+2 rule simply by not using the exterior points. */
    for (i=0; i<N; i++){
      h = 1.0/(N+1.0);
      t = h + i*h;
      x[i] = 1.0/t-1.0;
      (*function)(params_for_function,x[i],&y);
      w[i] = y*h/t/t;
    }
    return _SUCCESS_;
  }
  return _SUCCESS_;
}
int get_qsampling(double *x,
		  double *w,
		  int *N,
		  int N_max,
		  double rtol,
		  double *qvec,
		  int qsiz,
		  int (*test)(void * params_for_function, double q, double *psi),
		  int (*function)(void * params_for_function, double q, double *f0),
		  void * params_for_function,
		  ErrorMsg errmsg) {

  /* This routine returns the fewest possible number of abscissas and weights under
     the requirement that a test function folded with the neutrino distribution function
     can be integrated to an accuracy of rtol. If the distribution function is Fermi-Dirac
     or close, a Laguerre quadrature formula is often the best choice.

     This function combines two completely different strategies: Adaptive Gauss-Kronrod
     quadrature and Laguerres quadrature formula. */

  int i, NL=2,NR,level,Nadapt=0,NLag,NLag_max,Nold=NL;
  int adapt_converging=_FALSE_,Laguerre_converging=_FALSE_,combined_converging=_FALSE_;
  double y,y1,y2,I,Igk,err,ILag,*b,*c;
  qss_node *root,*root_comb;
  double I_comb,I_atzero,I_atinf,I_comb2;
  int N_comb=0,N_comb_lag=16,N_comb_leg=4;
  double a_comb,b_comb;
  double q_leg[4],w_leg[4];
  double q_lag[N_comb_lag],w_lag[N_comb_lag];
  char method_chosen[40];
  double qmin=0., qmax=0., qmaxm1=0.;
  double *wcomb2=NULL,delq;
  double Itot=0.0;
  int zeroskip=0;

  /* Set roots and weights for Gauss-Legendre 4 point rule: */
  q_leg[3] = sqrt((3.0+2.0*sqrt(6.0/5.0))/7.0);
  q_leg[2] = sqrt((3.0-2.0*sqrt(6.0/5.0))/7.0);
  q_leg[1] = -q_leg[2];
  q_leg[0] = -q_leg[3];
  w_leg[3] = (18.0-sqrt(30.0))/36.0;
  w_leg[2] = (18.0+sqrt(30.0))/36.0;
  w_leg[1] = w_leg[2];
  w_leg[0] = w_leg[3];

  /* Allocate storage for Laguerre coefficients: */
  class_alloc(b,N_max*sizeof(double),errmsg);
  class_alloc(c,N_max*sizeof(double),errmsg);
  /* If a vector of q values has been passed, use it: */
  if ((qvec!=NULL)&&(qsiz>1)){
    qmin = qvec[0];
    qmax = qvec[qsiz-1];
    qmaxm1 = qvec[qsiz-2];
  }
  else{
    qvec = NULL;
  }

  /* First do the adaptive quadrature - this will also give the value of the integral: */
  gk_adapt(&root,(*test),(*function), params_for_function,
	   rtol*1e-4, 1, 0.0, 1.0, _TRUE_, errmsg);
  /* Do a leaf count: */
  leaf_count(root);
  /* I can get the integral now: */
  I = get_integral(root, 1);
  //printf("I = %le |%le, used points: %d\n",I,I-1.0,15*root->leaf_childs);
  Itot += I;

  /* Starting from the top, move down in levels until tolerance is met: */
  for(level=root->leaf_childs; level>=1; level--){
    Igk = get_integral(root,level);
    err = I-Igk;
    if (fabs(err/Itot)<rtol) break;
  }
  if (level>0){
    /* Reduce tree to the found level:*/
    reduce_tree(root,level);
    /* Count the new leafs: */
    leaf_count(root);
    /* I know know how many function evaluations is
       required by the adaptively sampled grid:*/
    Nadapt = 15*root->leaf_childs;
    /* The adaptive routine could not recieve required precision
       using less than the required maximal number of points.*/
    if (Nadapt <= N_max) adapt_converging = _TRUE_;
  }


  /* Combined adaptive quadrature and Laguerre rescaled quadrature: */
  if(qvec!=NULL){
    /* Evaluate [0;qmin] using 4 point Gauss-Legendre rule: */
    (*function)(params_for_function,qmin,&y2);
    for(i=0,I_atzero=0.0; i<N_comb_leg; i++){
      q_leg[i] = 0.5*qmin*(1.0+q_leg[i]);
      w_leg[i] = w_leg[i]*0.5*qmin*y2;
      (*test)(params_for_function,q_leg[i],&y);
      I_atzero +=w_leg[i]*y;
    }

    /* Find asymptotic extrapolation:*/
    (*function)(params_for_function,qmaxm1,&y1);
    (*function)(params_for_function,qmax,&y2);

    b_comb = (y1/y2-1.0)/(qmax-qmaxm1);
    b_comb = MAX(b_comb,1e-100);
    //c_comb = -b_comb*qmax;
    a_comb = y2*exp(b_comb*qmax);
    // printf("f(q) = %g*exp(-%g*q) \n",a_comb,b_comb);
    //(*function)(params_for_function,100,&y2);
    //printf("f(100) = %e ?= %e\n",y2,a_comb*exp(-b_comb*100));

    /* Evaluate tail using 6 point Laguerre: */
    compute_Laguerre(q_lag,w_lag,N_comb_lag,0.0,b,c,_TRUE_);
    for (i=0,I_atinf=0.0; i<N_comb_lag; i++){
      w_lag[i] *= exp(-q_lag[i]);
      q_lag[i] = qmax + q_lag[i]/b_comb;
      w_lag[i] = a_comb/b_comb*exp(-b_comb*qmax)*w_lag[i];
      (*test)(params_for_function,q_lag[i],&y);
      I_atinf +=w_lag[i]*y;
    }

    /* Do the adaptive quadrature - this will also give the main part of the integral: */
    gk_adapt(&root_comb,(*test),(*function), params_for_function,
	     rtol*1e-2, 1, qmin, qmax, _FALSE_, errmsg);
    /* Do a leaf count: */
    leaf_count(root_comb);
    /* Starting from the top, move down in levels until tolerance is met: */
    for(level=root_comb->leaf_childs; level>=1; level--){
      I_comb = get_integral(root_comb,level);
      //printf("%le + %le + %le = %le | %le\n",
      //     I_atzero,I_atinf,I_comb,I_comb+I_atinf+I_atzero,I_comb+I_atinf+I_atzero-1.0);
      I_comb +=(I_atinf+I_atzero);
      err = I-I_comb;
      if (fabs(err/Itot)<rtol) break;
    }
    /* Reduce tree to the found level:*/
    if (level>0){
      reduce_tree(root_comb,level);
      /* Count the new leafs: */
      leaf_count(root_comb);
      N_comb = 15*root_comb->leaf_childs+N_comb_lag+N_comb_leg;
      if (N_comb <= N_max) combined_converging = _TRUE_;
    }

    /* Do the second combined quadrature: Same as above, but with trapezoidal rule
       using the given q grid:  */
    class_alloc(wcomb2,qsiz*sizeof(double),errmsg);
    I_comb2 = 0.0;
    for(i=0; i<qsiz; i++){
      if(i==0){
	delq = qvec[1]-qvec[0];
      }
      else if(i==qsiz-1){
	delq = qvec[qsiz-1]-qvec[qsiz-2];
      }
      else{
	delq = qvec[i+1]-qvec[i-1];
      }
      (*function)(params_for_function,qvec[i],&y2);
      wcomb2[i] = 0.5*y2*delq;
      (*test)(params_for_function,qvec[i],&y);
      I_comb2 +=wcomb2[i]*y;
    }
    I_comb2 +=(I_atzero+I_atinf);
    err = I - I_comb2;
    //    if(fabs(err/Itot)<rtol) combined2_converging= _TRUE_;
    //printf("I_comb2 = %e, rerr = %e\n",I_comb2,fabs(err/I));
  }


  /* Search for the minimal Laguerre quadrature rule: */
  NLag_max = MIN(N_max,80);
  for (NLag=NL; NLag<=NLag_max; NLag = MIN(NLag_max,NLag+10)){
    /* Evaluate integral: */
    compute_Laguerre(x,w,NLag,0.0,b,c,_TRUE_);
    ILag = 0.0;
    for (i=0; i<NLag; i++){
      (*test)(params_for_function,x[i],&y);
      (*function)(params_for_function,x[i],&y2);
      w[i] *= y2;
      ILag += y*w[i];
    }
    err = I-ILag;
    //fprintf(stderr,"\n Computing Laguerre, N=%d, I=%g and err=%g.\n",NLag,ILag,err);
    if (fabs(err/I)<rtol){
      Laguerre_converging = _TRUE_;
      break;
    }
    Nold = NLag;
    if (Nold == NLag_max) break;
  }

  if (Laguerre_converging == _TRUE_){
    /* We must refine NLag: */
    NL = Nold;
    NR = NLag;
    while ((NR-NL)>1) {
      NLag = (NL+NR)/2;
      /* Evaluate integral: */
      compute_Laguerre(x,w,NLag,0.0,b,c,_TRUE_);
      ILag = 0.0;
      for (i=0; i<NLag; i++){
	(*test)(params_for_function,x[i],&y);
	(*function)(params_for_function,x[i],&y2);
	w[i] *= y2;
	ILag += y*w[i];
      }
      err = I-ILag;
      //fprintf(stderr,"\n NLag=%d, rerr=%g.\n",NLag,fabs(err/I));
      if (fabs(err/Itot)<rtol){
	NR = NLag;
      }
      else{
	NL = NLag;
      }
    }
  }

  /* Choose best method if both works: */
  *N = N_max;
  //Laguerre_converging = _FALSE_;
  if (adapt_converging == _TRUE_) {
    *N = Nadapt;
  }
  if (combined_converging == _TRUE_){
    if(N_comb <= *N){
      *N = N_comb;
      adapt_converging = _FALSE_;
    }
    else{
      combined_converging = _FALSE_;
    }
  }
  if (Laguerre_converging == _TRUE_){
    if (NLag <= *N){
      *N = NLag;
      combined_converging = _FALSE_;
      adapt_converging = _FALSE_;
    }
    else{
      Laguerre_converging = _FALSE_;
    }
  }
  //printf("N_adapt=%d, N_combined=%d at level=%d, Nlag=%d\n",Nadapt,N_comb,level,NLag);
  if (adapt_converging==_TRUE_){
    class_sprintf(method_chosen,"Adaptive Gauss-Kronrod Quadrature");
    /* Gather weights and xvalues from tree: */
    i = Nadapt-1;
    get_leaf_x_and_w(root,&i,x,w,_TRUE_);
  }
  else if (Laguerre_converging==_TRUE_){
    class_sprintf(method_chosen,"Gauss-Laguerre Quadrature");
    /* x and w is already populated in this case. */
  }
  else if (combined_converging == _TRUE_){
    class_sprintf(method_chosen,"Combined Quadrature");
    for(i=0; i<N_comb_leg; i++){
      x[i] = q_leg[i];
      w[i] = w_leg[i];
    }
    i = N_comb_leg;
    get_leaf_x_and_w(root_comb,&i,x,w,_FALSE_);
    //printf("from %d to %d\n",N_comb_leg,i);

    for(i=0; i<N_comb_lag; i++){
      x[N_comb-N_comb_lag+i] = q_lag[i];
      w[N_comb-N_comb_lag+i] = w_lag[i];
    }
  }
  else{
    /* Failed to converge! */
    class_stop(errmsg,
		"get_qsampling fails to obtain a relative tolerance of %g as required using atmost %d points. If the PSD is interpolated from a file, try increasing the resolution and the q-interval (qmin;qmax) if possible, or decrease tol_ncdm/tol_ncdm_bg. As a last resort, increase _QUADRATURE_MAX_/_QUADRATURE_MAX_BG_.",rtol,N_max);
  }
  /* Trim weights to avoid zero weights: */
  for(i=0,zeroskip=0; i<*N; i++){
    for( ;(i<*N)&&(w[i+zeroskip]==0.0); zeroskip++,(*N)--);
    x[i] = x[i+zeroskip];
    w[i] = w[i+zeroskip];
  }


  //printf("Chosen sampling: %s, with %d points.\n",method_chosen,*N);
  //for(i=0; i<*N; i++) printf("(q,w) = (%g,%g)\n",x[i],w[i]);
  /* Deallocate tree: */
  burn_tree(root);
  if(qvec!=NULL){
    burn_tree(root_comb);
    free(wcomb2);
  }

  free(b);
  free(c);

  return _SUCCESS_;
}

int sort_x_and_w(double *x, double *w, double *workx, double *workw, int startidx, int endidx){
  int i,top=endidx,bot=startidx;
  double pivot;
  /* End recursion if only one element left in array: */
  if ((endidx-startidx)<1){
    return _SUCCESS_;
  }
  else{
    /*Copy x and w to workarray: */
    for (i=startidx; i<=endidx; i++){
      workx[i] = x[i];
      workw[i] = w[i];
    }
    pivot = x[endidx];
    //printf("pivot chosen: x[%d] = %g\n",endidx,pivot);
    for (i=startidx; i<endidx; i++){
      if (workx[i]<=pivot){
	//printf("<--%g  ",workx[i]);
	x[bot] = workx[i];
	w[bot++] = workw[i];
      }
      else{
	//printf("  %g-->",workx[i]);
	x[top] = workx[i];
	w[top--] = workw[i];
      }
    }
    //printf("\n top=%d, bot=%d, left=%d, right=%d\n",top,bot,startidx,endidx);
    x[top] = pivot;
    w[top] = workw[endidx];
    /* Recursive call: */
    sort_x_and_w(x,w,workx,workw,startidx,bot-1);
    sort_x_and_w(x,w,workx,workw,top+1,endidx);
    return _SUCCESS_;
  }

}

int get_leaf_x_and_w(qss_node *node, int *ind, double *x, double *w,int isindefinite){
  /* x and w should be exactly 15*root_node->leafchilds, and a leaf count should have
     been performed. Or perhaps I just use the fact that a leaf won't have children.
     Nah, let me use the leaf-count then. */
  int k;
  if (node->leaf_childs==1){
    for(k=0;k<15;k++){
      x[*ind] = node->x[k];
      w[*ind] = node->w[k];
      if (isindefinite == _TRUE_){
	(*ind)--;
      }
      else{
	(*ind)++;
      }
    }
  }
  else{
    /* Do recursive call: */
    get_leaf_x_and_w(node->left,ind,x,w,isindefinite);
    get_leaf_x_and_w(node->right,ind,x,w,isindefinite);
  }
  return _SUCCESS_;
}

int reduce_tree(qss_node *node, int level){
  /* Reduce the tree to a given level. Make all nodes with
     node->leaf_childs==level into leafs.
     If we call reduce_tree(root,1), nothing happens.*/
  if(node->leaf_childs==level){
    burn_tree(node->left);
    burn_tree(node->right);
    node->left = NULL;
    node->right = NULL;
  }
  else if(node->leaf_childs>level){
    /* else try to see if children nodes can be simplified: */
    reduce_tree(node->left,level);
    reduce_tree(node->right,level);
  }
  /* If called on a node which has leaf_childs<level, it does nothing. */
  return _SUCCESS_;
}


int burn_tree(qss_node *node){
  /* Burn node and all subnodes. */
  /* Call burn_branch recursively on children nodes: */
  /* This node and all its subnodes */
  if (node!=NULL){
    if (node->left!=NULL) burn_tree(node->left);
    if (node->right!=NULL) burn_tree(node->right);

    if (node->x!=NULL) free(node->x);
    if (node->w!=NULL) free(node->w);
    free(node);
  }
  return _SUCCESS_;
}

int leaf_count(qss_node *node){
  /* Count the amount of leafs under a given node and write the number in the node. */
  /* We call recursively, until a node is a leaf - then we add the numbers on our
     way back:*/
  if (node->left!=NULL){
    /* This is not a leaf, do recursive call: */
    leaf_count(node->left);
    leaf_count(node->right);
    node->leaf_childs = node->left->leaf_childs + node->right->leaf_childs;
    return _SUCCESS_;
  }
  else{
    /* This is a leaf, by definition leaf_childs = 1: */
    node->leaf_childs = 1;
    return _SUCCESS_;
  }
}

double get_integral(qss_node *node, int level){
  /* Traverse the tree and return the estimate of the integral at a given level.
     level 1 is the best estimate. */
  double IL,IR;
  /* An updated leaf_count is assumed. */
  if (node->leaf_childs<=level){
    return node->I;
  }
  else{
    IL = get_integral(node->left, level);
    IR = get_integral(node->right, level);
    /* Combine the integrals: */
    return (IL+IR);
  }
}



int gk_adapt(
	     qss_node** node,
	     int (*test)(void * params_for_function, double q, double *psi),
	     int (*function)(void * params_for_function, double q, double *f0),
	     void * params_for_function,
	     double tol,
	     int treemode,
	     double a,
	     double b,
	     int isindefinite,
	     ErrorMsg errmsg){
  /* Do adaptive Gauss-Kronrod quadrature, while building the
     recurrence tree. If treemode!=0, store x-values and weights aswell.
     At first call, a and b should be 0 and 1 if isdefinite==_TRUE_. */
  double mid;
  /* Allocate current node: */
  class_alloc(*node,sizeof(qss_node),errmsg);
  if (treemode==0){
    (*node)->x = NULL;
    (*node)->w = NULL;
  }
  else{
    class_alloc((*node)->x,15*sizeof(double),errmsg);
    class_alloc((*node)->w,15*sizeof(double),errmsg);
  }
  (*node)->left = NULL; (*node)->right = NULL;

  gk_quad((*test), (*function), params_for_function, *node, a, b, isindefinite);
  if ((fabs((*node)->err/(*node)->I) < tol)||(tol>=1.0)){
    /* Stop recursion and return. tol>=1.0 in case of I=0 infinite recursion */
    return _SUCCESS_;
  }
  else{
    /* Call gk_adapt recursively on children:*/
    mid = 0.5*(a+b);
    //printf("<-%g,%g,%g,%g",mid,tol,(*node)->err,(*node)->I);
    gk_adapt(&((*node)->left),(*test),(*function), params_for_function, 1.5*tol,
	     treemode, a, mid, isindefinite, errmsg);
    //printf("%g->",mid);
    gk_adapt(&((*node)->right),(*test),(*function), params_for_function, 1.5*tol,
	     treemode, mid, b, isindefinite, errmsg);
    /* Update integral and error in this node and return: */
    /* Actually, it is more convenient just to keep the nodes own estimate of the
       integral for our purposes.
       (*node)->I = (*node)->left->I + (*node)->right->I;
       (*node)->err = sqrt(pow(node->left->err,2)+pow(node->right->err,2));
    */
    return _SUCCESS_;
  }
}

int compute_Hermite(double *x, double *w, int N, int alpha, double *b, double *c){
  int NLag,i;
  double alpha_Lag;

  NLag = N/2;
  /* In case N is uneven, zero the N'th weight:*/
  w[N-1] = 0.0;
  alpha_Lag = (alpha-1.0)/2.0;

  /* Compute the positive roots and weights (up to some simple manipulation): */
  compute_Laguerre(x+NLag,w+NLag,NLag,alpha_Lag,b,c,_FALSE_);

  /* Do manipulations:*/
  for(i=NLag; i<2*NLag; i++){
    x[i] = sqrt(x[i]);
    w[i] *=0.5;
  }

  /* Set the negative roots and weights:*/
  for(i=0; i<NLag; i++){
    x[i] = -x[2*NLag-i-1];
    w[i] = w[2*NLag-i-1];
    if (alpha%2!=0){
      w[i] = -w[i];
    }
  }
  return _SUCCESS_;
}


int compute_Laguerre(double *x, double *w, int N, double alpha, double *b, double *c,int totalweight){
  int i,j,iter,maxiter=10;
  double x0=0.,r1,r2,ratio,d,logprod,logcc;
  double p0,p1,p2,dp0,dp1,dp2;
  double eps=1e-14;
  /* Initialise recursion coefficients: */
  for(i=0; i<N; i++){
    b[i] = alpha + 2.0*i +1.0;
    c[i] = i*(alpha+i);
  }
  logprod = 0.0;
  for(i=1; i<N; i++) logprod +=log(c[i]);
  logcc = lgamma(alpha+1)+logprod;

  /* Loop over roots: */
  for (i=0; i<N; i++){
    /* Estimate root: */
    if (i==0) {
      x0 =(1.0+alpha)*(3.0+0.92*alpha)/( 1.0+2.4*N+1.8*alpha);
    }
    else if (i==1){
      x0 += (15.0+6.25*alpha)/( 1.0+0.9*alpha+2.5*N);
    }
    else{
      r1 = (1.0+2.55*(i-1))/( 1.9*(i-1));
      r2 = 1.26*(i-1)*alpha/(1.0+3.5*(i-1));
      ratio = (r1+r2)/(1.0+0.3*alpha);
      x0 += ratio*(x0-x[i-2]);
    }
    /* Refine root using Newtons method: */
    for(iter=1; iter<=maxiter; iter++){
      /* We need to find p2=L_N(x0), dp2=L'_N(x0) and
	 p1 = L_(N-1)(x0): */
      p1 = 1.0;
      dp1 = 0.0;
      p2 = x0 - alpha - 1.0;
      dp2 = 1.0;
      for (j=1; j<N; j++ ){
	p0 = p1;
	dp0 = dp1;
	p1 = p2;
	dp1 = dp2;
	p2  = (x0-b[j])*p1 - c[j]*p0;
	dp2 = (x0-b[j])*dp1 + p1 - c[j]*dp0;
      }
      /* New guess at root: */
      d = p2/dp2;
      x0 -= d;
      if (fabs(d)<=eps*(fabs(x0)+1.0)) break;
    }
    /* Okay, write root and weight: */
    x[i] = x0;

    if (totalweight == _TRUE_)
      w[i] = exp(x0+logcc-log(dp2*p1));
    else
       w[i] = exp(logcc-log(dp2*p1));
  }

  return _SUCCESS_;

}




int gk_quad(int (*test)(void * params_for_function, double q, double *psi),
	    int (*function)(void * params_for_function, double q, double *f0),
	    void * params_for_function,
	    qss_node* node,
	    double a,
	    double b,
	    int isindefinite){
  const double z_k[15]={-0.991455371120813,
			-0.949107912342759,
			-0.864864423359769,
			-0.741531185599394,
			-0.586087235467691,
			-0.405845151377397,
			-0.207784955007898,
			0.0,
			0.207784955007898,
			0.405845151377397,
			0.586087235467691,
			0.741531185599394,
			0.864864423359769,
			0.949107912342759,
			0.991455371120813};
  const double w_k[15]={0.022935322010529,
			0.063092092629979,
			0.104790010322250,
			0.140653259715525,
			0.169004726639267,
			0.190350578064785,
			0.204432940075298,
			0.209482141084728,
			0.204432940075298,
			0.190350578064785,
			0.169004726639267,
			0.140653259715525,
			0.104790010322250,
			0.063092092629979,
			0.022935322010529};
  const double w_g[7]={0.129484966168870,
		       0.279705391489277,
		       0.381830050505119,
		       0.417959183673469,
		       0.381830050505119,
		       0.279705391489277,
		       0.129484966168870};
  int i,j;
  double x,wg,wk,t,Ik,Ig,y,y2;

  /* 	Loop through abscissas, transform the interval and form the Kronrod
     15 point estimate of the integral.
     Every second time we update the Gauss 7 point quadrature estimate. */

  Ik=0.0;
  Ig=0.0;
  for (i=0;i<15;i++){
    /* Transform z into t in interval between a and b: */
    t = 0.5*(a*(1-z_k[i])+b*(1+z_k[i]));
    /* Modify weight such that it reflects the linear transformation above: */
    wk = 0.5*(b-a)*w_k[i];
    if (isindefinite==_TRUE_){
      /* Transform t into x in interval between 0 and inf: */
      x = 1.0/t-1.0;
      /* Modify weight accordingly: */
      wk = wk/(t*t);
    }
    else{
      x = t;
    }
    (*test)(params_for_function,x,&y);
    (*function)(params_for_function,x,&y2);
    wk *= y2;
    /* Update Kronrod integral: */
    Ik +=wk*y;
    /* If node->x and node->w is allocated, store values: */
    if (node->x!=NULL) node->x[i] = x;
    if (node->w!=NULL) node->w[i] = wk;
    /* If i is uneven, update Gauss integral: */
    if ((i%2)==1){
      j = (i-1)/2;
      /* Transform weight according to linear transformation: */
      wg = 0.5*(b-a)*w_g[j];
      if (isindefinite == _TRUE_){
        /* Transform weight according to non-linear transformation x = 1/t -1: */
        wg = wg/(t*t);
      }
      /* Update integral: */
      Ig +=wg*y*y2;
    }
  }
  node->err = pow(200*fabs(Ik-Ig),1.5);
  node->I = Ik;
  return _SUCCESS_;
}

/**
 * This routine computes the weights and abscissas of a Gauss-Legendre quadrature between -1 and 1
 *
 * @param mu     Input/output: Vector of cos(beta) values
 * @param w8     Input/output: Vector of quadrature weights
 * @param n      Input       : Number of quadrature points
 * @param tol    Input       : tolerance on each mu
 *
 * From Numerical recipes
 **/

int quadrature_gauss_legendre(
			      double *mu,
			      double *w8,
			      int n,
			      double tol,
			      ErrorMsg error_message) {

  int m,j,i,counter;
  double z1,z,pp,p3,p2,p1;

  m=(n+1)/2;
  for (i=1;i<=m;i++) {
    z=cos(_PI_*((double)i-0.25)/((double)n+0.5));
    counter=0;
    do {
      p1=1.0;
      p2=0.0;
      for (j=1;j<=n;j++) {
        p3=p2;
        p2=p1;
        p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp=n*(z*p1-p2)/(z*z-1.0);
      z1=z;
      z=z1-p1/pp;
      counter++;
      class_test(counter == _MAX_IT_,
		 error_message,
		 "maximum number of iteration reached: increase either _MAX_IT_ or tol\n");
    } while (fabs(z-z1) > tol);
    mu[i-1]=-z;
    mu[n-i]=z;
    w8[i-1]=2.0/((1.0-z*z)*pp*pp);
    w8[n-i]=w8[i-1];

  }
  return _SUCCESS_;
}

int quadrature_in_rectangle(
			    double xl,
			    double xr,
			    double yl,
			    double yr,
			    int *n,
			    double ** x,
			    double ** y,
			    double ** w,
			    ErrorMsg error_message) {


  double xl_tile,xr_tile,yl_tile,yr_tile;
  int N;

  N=24;
  xl_tile = xl;
  xr_tile = xr;
  yl_tile = yl;
  yr_tile = yr;

  *n = N;


  class_alloc(*x,sizeof(double)*N,error_message);
  class_alloc(*y,sizeof(double)*N,error_message);
  class_alloc(*w,sizeof(double)*N,error_message);
  class_call(cubature_order_eleven(xl_tile,
				   xr_tile,
				   yl_tile,
				   yr_tile,
				   *x+0,
				   *y+0,
				   *w+0,
				   error_message),
	     error_message,
	     error_message);


  return _SUCCESS_;
}



int cubature_order_eleven(
			  double xl,
			  double xr,
			  double yl,
			  double yr,
			  double *x,
			  double *y,
			  double *w,
			  ErrorMsg error_message){

  double wi[6]={0.48020763350723814563e-01,
		0.66071329164550595674e-01,
		0.97386777358668164196e-01,
		0.21173634999894860050e+00,
		0.22562606172886338740e+00,
		0.35115871839824543766e+00};
  double xi[6]={0.98263922354085547295e+00,
		0.82577583590296393730e+00,
		0.18858613871864195460e+00,
		0.81252054830481310049e+00,
		0.52532025036454776234e+00,
		0.41658071912022368274e-01};
  double yi[6]={0.69807610454956756478e+00,
		0.93948638281673690721e+00,
		0.95353952820153201585e+00,
		0.31562343291525419599e+00,
		0.71200191307533630655e+00,
		0.42484724884866925062e+00};

  int idx,i;
  double a1,a2,b1,b2;

  a1 = 2./(xr-xl);
  a2 = 2./(yr-yl);
  b1 = 1.-2*xr/(xr-xl);
  b2 = 1.-2*yr/(yr-yl);

  for (i=0,idx=0; i<6; i++,idx++){
    // Upper right corner:
    x[idx] = (xi[i]-b1)/a1;
    y[idx] = (yi[i]-b2)/a2;
    w[idx] = wi[i]/a1/a2;
  }
  for (i=0,idx=6; i<6; i++,idx++){
    // Upper left corner:
    x[idx] = (-yi[i]-b1)/a1;
    y[idx] = (xi[i]-b2)/a2;
    w[idx] = wi[i]/a1/a2;
  }
  for (i=0,idx=12; i<6; i++,idx++){
    // Lower left corner:
    x[idx] = (-xi[i]-b1)/a1;
    y[idx] = (-yi[i]-b2)/a2;
    w[idx] = wi[i]/a1/a2;
  }
  for (i=0,idx=18; i<6; i++,idx++){
    // Lower right corner:
    x[idx] = (yi[i]-b1)/a1;
    y[idx] = (-xi[i]-b2)/a2;
    w[idx] = wi[i]/a1/a2;
  }



  return _SUCCESS_;
}

