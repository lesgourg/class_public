/******************************************/
/* Quadrature Sampling Strategy for CLASS */
/* 10/12 2010                             */
/* Thomas Tram                            */
/******************************************/
#include "quadrature.h"

int get_qsampling(double *x,
		  double *w, 
		  int *N, 
		  int N_max, 
		  double rtol,
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
	
  int i, NL=5,NR,level,Nadapt,NLag,Nold=0;
  int adapt_converging=_TRUE_,Laguerre_converging=_TRUE_;
  double y,y2,I,Igk,err,ILag,*b,*c;
  qss_node* root;
  /* Allocate storage for Laguerre coefficients: */
  b = malloc(N_max*sizeof(double));
  c = malloc(N_max*sizeof(double));
  /* First do the adaptive quadrature - this will also give the value of the integral: */
  root = gk_adapt((*test),(*function), params_for_function, rtol*1e-2, 1, 0.0, 1.0);
  /* Do a leaf count: */
  leaf_count(root);
  /* I can get the integral now: */
  I = get_integral(root, 1);
  /* Starting from the top, move down in levels until tolerance is met: */
  for(level=root->leaf_childs; level>=1; level--){
    Igk = get_integral(root,level);
    err = I-Igk;
    if (fabs(err/Igk)<rtol) break;
  }
  /* Reduce tree to the found level:*/
  reduce_tree(root,level);
  /* Count the new leafs: */
  leaf_count(root);
  /* I know know how many function evaluations is 
     required by the adaptively sampled grid:*/
  Nadapt = 15*root->leaf_childs;
  /* The adaptive routine could not recieve required precision 
     using less than the required maximal number of points.*/
  if (Nadapt > N_max) adapt_converging = 0;
  /* Search for the minimal Laguerre quadrature rule: */
  for (NLag=NL; NLag<=N_max; NLag = min(N_max,NLag+15)){
    /* Evaluate integral: */
    compute_Laguerre(x,w,NLag,0.0,b,c);
    ILag = 0.0;
    for (i=0; i<NLag; i++){
      (*test)(params_for_function,x[i],&y);
      (*function)(params_for_function,x[i],&y2);
      w[i] *= y2;
      ILag += y*w[i];
    }
    err = I-ILag;
    //fprintf(stderr,"\n Computing Laguerre, N=%d, I=%g and err=%g.\n",NLag,ILag,err);
    if (fabs(err/I)<rtol) break;
    if (NLag == N_max){
      Laguerre_converging = _FALSE_;
      break;
    }
    Nold = NLag;
  }
  if (Laguerre_converging == _TRUE_){
    /* We must refine NLag: */
    NL = Nold;
    NR = NLag;
    while ((NR-NL)>1) {
      NLag = (NL+NR)/2;
      /* Evaluate integral: */
      compute_Laguerre(x,w,NLag,0.0,b,c);
      ILag = 0.0;
      for (i=0; i<NLag; i++){
	(*test)(params_for_function,x[i],&y);
	(*function)(params_for_function,x[i],&y2);
	w[i] *= y2;
	ILag += y*w[i];
      }
      err = I-ILag;
      //fprintf(stderr,"\n NLag=%d, rerr=%g.\n",NLag,fabs(err/I));
      if (fabs(err/I)<rtol){
	NR = NLag;
      }
      else{
	NL = NLag;
      }
    }
  }
	
  /* Choose best method if both works: */
  if ((adapt_converging==_TRUE_)&&(Laguerre_converging==_TRUE_)){
    if (Nadapt<=NLag){
      Laguerre_converging = _FALSE_;
    }
    else{
      adapt_converging = _FALSE_;
    }
  }
  if (adapt_converging==_TRUE_){
    /* Gather weights and xvalues from tree: */
    i = 0;
    get_leaf_x_and_w(root,&i,x,w);
    *N = Nadapt;
  }
  else if (Laguerre_converging==_TRUE_){
    /* x and w is already populated in this case. */
    *N = NLag;
  }
  else{
    /* Failed to converge! */
    /* status = _FAILURE_; */
    class_test(0==0,
	       errmsg,
	       "Fails to converge! Check that your distribution function is not weird, then increase _QUADRATURE_MAX_ or decrease tol_ncdm_integration");

  }
  /* Deallocate tree: */
  burn_tree(root);
  free(b);
  free(c);
  return _SUCCESS_;
}
	
	

int get_leaf_x_and_w(qss_node *node, int *ind, double *x, double *w){
  /* x and w should be exactly 15*root_node->leafchilds, and a leaf count should have
     been performed. Or perhaps I just use the fact that a leaf won't have children.
     Nah, let me use the leaf-count then. */
  int k;
  if (node->leaf_childs==1){
    for(k=0;k<15;k++){
      x[*ind] = node->x[k];
      w[*ind] = node->w[k];
      (*ind)++;
    }
  }
  else{
    /* Do recursive call: */
    get_leaf_x_and_w(node->left,ind,x,w);
    get_leaf_x_and_w(node->right,ind,x,w);
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
	


qss_node* gk_adapt(
		   int (*test)(void * params_for_function, double q, double *psi),
		   int (*function)(void * params_for_function, double q, double *f0),
		   void * params_for_function,
		   double tol, 
		   int treemode, 
		   double a, 
		   double b){
  /* Do adaptive Gauss-Kronrod quadrature, while building the
     recurrence tree. If treemode!=0, store x-values and weights aswell.
     At first call, a and b should be 0 and 1. */
  double mid;
  qss_node *node;
	
  /* Allocate current node: */
  node = malloc(sizeof(qss_node));
  if (treemode==0){
    node->x = NULL;
    node->w = NULL;
  }
  else{
    node->x = malloc(15*sizeof(double));
    node->w = malloc(15*sizeof(double));
  }
  node->left = NULL; node->right = NULL;
	
  gk_quad((*test), (*function), params_for_function, node, a, b);
  if (node->err/node->I < tol){
    /* Stop recursion and return pointer to node: */
    return node;
  }
  else{
    /* Call gk_adapt recursively on children:*/
    mid = 0.5*(a+b);
    node->left = gk_adapt((*test),(*function), params_for_function, 1.5*tol, treemode, a, mid);
    node->right = gk_adapt((*test),(*function), params_for_function, 1.5*tol, treemode, mid, b);
    /* Update integral and error in this node and return: */
    /* Actually, it is more convenient just to keep the nodes own estimate of the
       integral for our purposes.
       node->I = node->left->I + node->right->I;
       node->err = sqrt(pow(node->left->err,2)+pow(node->right->err,2));
    */
    return node;
  }
}
	
int compute_Laguerre(double *x, double *w, int N, double alpha, double *b, double *c){
  int i,j,iter,maxiter=10;
  double prod,cc,x0=0.,r1,r2,ratio,d,logprod,logcc;
  double p0,p1,p2,dp0,dp1,dp2;
  double eps=1e-14;
  /* Initialise recursion coefficients: */
  for(i=0; i<N; i++){
    b[i] = alpha + 2.0*i +1.0;
    c[i] = i*(alpha+i);
  }
  prod=1.0;
  logprod = 0.0;
  for(i=1; i<N; i++) logprod +=log(c[i]);
  prod = exp(logprod);
  logcc = lgamma(alpha+1)+logprod;
  cc = exp(logcc);
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
    /*		w[i] = (cc/dp2)/p1;		*/
    w[i] = exp(x0+logcc-log(dp2*p1));
    /*printf("\n logcc=%g, logdp2=%g, 
      logp1=%g.",logcc,log(dp2),log(p1));
      printf("\n cc=%g, dp2=%g, p1=%g, prod=%g,tg(1)=%g",
      cc,dp2,p1,prod,tgamma(1.0));*/
  }

  return _SUCCESS_;

}


int gk_quad(int (*test)(void * params_for_function, double q, double *psi),
	    int (*function)(void * params_for_function, double q, double *f0),
	    void * params_for_function,
	    qss_node* node, 
	    double a, 
	    double b){
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
    /* Transform t into x in interval between 0 and inf: */
    x = 1.0/t-1.0;
    /* Modify weight accordingly: */
    wk = wk/(t*t);
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
      /* Transform weight according to non-linear transformation x = 1/t -1: */
      wg = wg/(t*t);
      /* Update integral: */
      Ig +=wg*y*y2;
    }
  }
  node->err = pow(200*fabs(Ik-Ig),1.5);
  node->I = Ik;
  return _SUCCESS_;
}
	
