/** @file class.c 
 * Julien Lesgourgues, 20.04.2010    
 */
#include "class.h"
#include <math.h> 
#include <time.h>
#define TINy 1.0e-1
#define NMAX 100000
#define GET_PSUM \
  for (j=0;j<ndim;j++){\
    for (sum=0.0,i=0;i<mpts;i++) sum += p[i][j];\
    psum[j] = sum;}
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
#define _NPARAMS_ 5
#define l_max 2000

struct precision pr;        /* for precision parameters */
struct background ba;       /* for cosmological background */
struct thermo th;           /* for thermodynamics */
struct perturbs pt;         /* for source functions */
struct bessels bs;          /* for bessel functions */
struct transfers tr;        /* for transfer functions */
struct primordial pm;       /* for primordial spectra */
struct spectra sp;          /* for output spectra */
struct nonlinear nl;	      /* for non linear spectra */
struct lensing le;          /* for lensed spectra */
struct output op;           /* for output files */
struct file_content fc;

double ** cl_ref;
double ** noise;
ErrorMsg errmsg;



double get_chi2( double * param){

  double **cl;
  int l;
  cl = calloc((l_max+1),sizeof(double*));
  for (l=2; l <= l_max; l++) {
    cl[l] = calloc(4,sizeof(double));
  }
  double chi2;
  int i;
  /*fprintf(stderr,"YHe:%e, h:%e, omega_cdm:%e\n",param[0],param[1],param[2]);*/
  sprintf(fc.value[4],"%g",param[0]);
  /*fprintf(stderr,"fc value for h is %g\n",param[1]);*/
  sprintf(fc.value[5],"%g",param[1]);
  sprintf(fc.value[6],"%g",param[2]);
  sprintf(fc.value[8],"%g",param[3]);
  sprintf(fc.value[9],"%g",param[4]);
  if (input_init(&fc,&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&nl,&le,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }

  class_assuming_bessels_computed(&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&nl,&le,&op,cl,errmsg);
  chi2_planck(&sp,cl,cl_ref,noise,&chi2);
  
  for (l=2; l <= l_max; l++) {
    free(cl[l]);
  }
  free(cl);
  /*fprintf(stderr,"here, YHe = %s, h = %s, omega_c = %s, n_s = %s, A_s = %s, chi2 = %e\n",fc.value[4],fc.value[5],fc.value[6],fc.value[8],fc.value[9],chi2);*/

  /*fprintf(stderr,"YHe:%e, h:%e, omega_cdm:%e, chi2:%e\n",param[0],param[1],param[2],chi2);*/
  return chi2;
}

void amoeba(double **p,double y[],int ndim,double ftol, double (*funk)(double [])){
  
  double amotry(double **p,double y[], double psum[], int ndim, double (*funk)(double []),int ihi,double fac);
  int i,ihi,ilo,inhi,j,mpts=ndim+1;
  int k;
  double rtol,sum,swap,ysave,ytry,*psum;
  psum = calloc(ndim,sizeof(double));
  int nfunk=0;
  GET_PSUM
  for (;;) {
      ilo=0;
      /*First we must determine which point is the highest (worst), next-highest, and lowest (best), by looping over the points in the simplex.*/
      ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
      for (i=0;i<mpts;i++) {
	if (y[i] <= y[ilo]) 
	  ilo=i; 
	if (y[i] > y[ihi]) {
	  inhi=ihi;
	  ihi=i;
	} 
	else if (y[i] > y[inhi] && i != ihi) inhi=i;
      }
      rtol=fabs(y[ihi]-y[ilo]); //Compute the fractional range from highest to lowest and return if satisfactory.
      /*fprintf(stderr,"--> y[ihi]:%e y[ilo]:%e rtol %e\n\n",y[ihi],y[ilo],rtol);*/
	if (rtol < ftol) { //If returning, put best point and value in slot 1.
	  SWAP(y[0],y[ilo]);
	  for (i=0;i<ndim;i++) 
	    SWAP(p[0][i],p[ilo][i]);
	  break;
	}
      if (nfunk >= NMAX) fprintf(stderr,"NMAX exceeded\n");
      nfunk += 2;
      //Begin a new iteration. First extrapolate by a factor −1 through the face of the simplex across from the high point, i.e., reflect the simplex from the high point. 
    
      ytry=amotry(p,y,psum,ndim,funk,ihi,-1.0);
      if (ytry <= y[ilo]){
	//Gives a result better than the best point, so try an additional extrapolation by a factor 2.
	/*fprintf(stderr,"is better, going twice as far\n");*/
	ytry=amotry(p,y,psum,ndim,funk,ihi,2.0);}
      else if (ytry >= y[inhi]) {
	/*fprintf(stderr,"is worse, trying half\n");*/
	//The reflected point is worse than the second-highest, so look for an intermediate lower point, i.e., do a one-dimensional contraction.
	ysave=y[ihi];
	ytry=amotry(p,y,psum,ndim,funk,ihi,0.5);
	if (ytry >= ysave) { //Can’t seem to get rid of that high point. Better 
	  for (k=0;k<mpts;k++) { //contract around the lowest (best) point.
	    /*fprintf(stderr,"ilo:%d, k:%d\n",ilo,k);*/
	    if (k != ilo) {
	      for (j=0;j<ndim;j++){
		/*fprintf(stderr,"---- %g %g -- %g --\n",p[k][j],p[ilo][j],psum[0]);*/
		psum[j] = 0.5*(p[k][j]+p[ilo][j]);
		/*fprintf(stderr,"-SUM IS--- %g \n",sum);*/
		p[k][j] = psum[j];
		/*fprintf(stderr,"---- %g -----\n",p[k][j]);*/
	      }
	      y[k]=(*funk)(psum);
	    }
	  }
	  nfunk += ndim; 
	  GET_PSUM
	}
      } 
      else --(nfunk);
    } 
  free(psum);
}
  /*Multidimensional minimization of the function funk(x) where x[1..ndim] is a*/
  /*vector in ndim dimensions, by the downhill simplex method of Nelder and Mead.*/
  /*The matrix p[1..ndim+1] [1..ndim] is input. Its ndim+1 rows are*/
  /*ndim-dimensional vectors which are the vertices of the starting simplex.*/
  /*Also input is the vector y[1..ndim+1], whose components must be pre-*/
  /*initialized to the values of funk evaluated at the ndim+1 vertices (rows) of*/
  /*p; and ftol the fractional convergence tolerance to be achieved in the*/
  /*function value (n.b.!). On output, p and y will have been reset to ndim+1 new*/
  /*points all within ftol of a minimum function value, and nfunk gives the*/
  /*number of function evaluations taken.*/

/*Extrapolates by a factor fac through the face of the simplex across from the
 * high point, tries it, and replaces the high point if the new point is
 * better.*/
double amotry(double **p, double y[], double psum[], int ndim, double (*funk)(double []), int ihi, double fac)
{
  int j;
  double fac1,fac2,ytry,*ptry;
  ptry = malloc(ndim*sizeof(double));
  fac1=(1.0-fac)/(ndim*1.0);
  fac2=fac1-fac;
  for (j=0;j<ndim;j++){
    ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;}
  ytry=(*funk)(ptry); //Evaluate the function at the trial point.
  /*fprintf(stderr,"for YHe:%e, y:%e\n",ptry[0],ytry);*/
  if (ytry < y[ihi]) { //If it’s better than the highest, then replace the highest.
    y[ihi]=ytry;
    for (j=0;j<ndim;j++) {
      psum[j] += ptry[j]-p[ihi][j];
      p[ihi][j]=ptry[j]; }
  } 
  free(ptry); 
  return ytry;
}

int class_assuming_bessels_computed(
				    struct precision * ppr,
				    struct background * pba,
				    struct thermo * pth,
				    struct perturbs * ppt,
				    struct bessels * pbs,
				    struct transfers * ptr,
				    struct primordial * ppm,
				    struct spectra * psp,
				    struct nonlinear *pnl,
				    struct lensing *ple,
				    struct output * pop,
				    double ** cl,
				    ErrorMsg errmsg);

int chi2_simple(
		struct spectra * psp,
		double ** cl1,
		double ** cl2,
		double * chi2);

int chi2_planck(
		struct spectra * psp,
		double ** cl1,
		double ** cl2,
		double ** noise,
		double * chi2);

int noise_planck(
		 struct background * pba,
		 struct thermo * pth,
		 struct spectra * psp,
		 double ** noise
		 );

int main(int argc, char **argv) {


  int i,l,j,k;
  double storage;

  double parameter_initial,parameter_step;
  double * parameter;
  int param_num;

  double ** cl;
  double *chi2;
  double chi2_temp;
  double percentage,max_percentage;

  int ref_run;
  
  parser_init(&fc,11,errmsg);

  strcpy(fc.name[0],"output");
  strcpy(fc.value[0],"tCl,pCl");

  strcpy(fc.name[1],"l_max_scalars");
  sprintf(fc.value[1],"%d",l_max);
 
/*   strcpy(fc.name[2],"perturbations_verbose"); */
/*   sprintf(fc.value[2],"%d",2); */

/*******************************************************/

  strcpy(fc.name[2],"lensing");
  strcpy(fc.value[2],"no");

  strcpy(fc.name[3],"N_eff");
  strcpy(fc.name[4],"YHe");

  strcpy(fc.name[5],"h");
  strcpy(fc.name[6],"omega_cdm");
  strcpy(fc.name[7],"omega_b");
  strcpy(fc.name[8],"n_s");
  strcpy(fc.name[9],"A_s");
  strcpy(fc.name[10],"recombination");
  strcpy(fc.value[10],"RECFAST");


/*******************************************************/
  // Fixed parameter: Neff
  parameter_initial=3.046;
  parameter_step=0.02;

  param_num=41;
  ref_run=20;

/*******************************************************/
  // Varied parameters with Neff: omega_c and h
  double omega_cdm;
  double omega_b; 
  double h;
  double YHe;
  double A_s;
  double n_s;

  omega_cdm = 0.112;
  omega_b   = 0.027;
  h = 0.69;
  YHe = 0.25;
  n_s = 0.967;
  A_s = 2.3e-9;

/*******************************************************/
  // Minimisation of the chi2 over YHe
  // Quantities for amoebia

  double **p;
  double * starting_values;
  int ndim = _NPARAMS_;
  class_calloc(starting_values,ndim+1,sizeof(double),errmsg);
  class_calloc(p,ndim+1,sizeof(double*),errmsg);
  for (i=0;i<=ndim;i++){
    class_calloc(p[i],ndim,sizeof(double),errmsg);}
  

/*******************************************************/
  // defining starting values for YHe
  starting_values[0] = 0.245;
  starting_values[1] = 0.255;
/*******************************************************/

  class_calloc(cl,(l_max+1),sizeof(double*),errmsg);
  class_calloc(cl_ref,(l_max+1),sizeof(double*),errmsg);
  class_calloc(noise,(l_max+1),sizeof(double*),errmsg);
  for (l=2; l <= l_max; l++) {
    class_calloc(noise[l],4,sizeof(double),errmsg);
    class_calloc(cl[l],4,sizeof(double),errmsg);
    class_calloc(cl_ref[l],4,sizeof(double),errmsg);
  }
  class_calloc(chi2,ndim+1,sizeof(double),errmsg);

  class_alloc(parameter,param_num*sizeof(double),errmsg);

  double ftol = 0.1;
  int ihi;

  // Create values for parameter to vary
  for (i=0; i<param_num; i++) {

    parameter[i] = (parameter_initial - (param_num*1.0-1.)*parameter_step/2.) + i*parameter_step ;
  }

  // Create reference run
  sprintf(fc.value[3],"%g",parameter[ref_run]);
  sprintf(fc.value[4],"%g",YHe);
  sprintf(fc.value[5],"%g",h);
  sprintf(fc.value[6],"%g",omega_cdm);
  sprintf(fc.value[7],"%g",omega_b);
  sprintf(fc.value[8],"%g",n_s);
  sprintf(fc.value[9],"%g",A_s);

  if (input_init(&fc,&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&nl,&le,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }
  noise_planck(&ba,&th,&sp,noise);
  class_assuming_bessels_computed(&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&nl,&le,&op,cl_ref,errmsg);
  double alpha;
  double random_number;

  // Looping on all values of parameter
  for (i=0; i<param_num; i++) {

    /*alpha =(1.+0.2271*parameter[i])/(1.+0.2271*parameter[ref_run]);*/
    //initialisation of the minimization in 1d
    /*p[0][0] = starting_values[0];*/
    /*p[1][0] = starting_values[1];*/

    srand(time(NULL));
    sprintf(fc.value[3],"%g",parameter[i]);
    fprintf(stderr,"#run %d/%d with %s\n",i+1,param_num,fc.value[3]);
    for (j=0; j<=4; j++){
      random_number = random()*0.1/RAND_MAX; 
      //initialisation of the minimization in 5d
      // YHe
      p[0][0] = 0.245*(1.-random_number);
      p[1][0] = 0.255*(1.+random_number);
      p[2][0] = 0.245*(1.-random_number);
      p[3][0] = 0.245*(1.-random_number);
      p[4][0] = 0.245*(1.-random_number);
      p[5][0] = 0.245*(1.-random_number);

      // h
      p[0][1] = 0.68*(1.-random_number);
      p[1][1] = 0.68*(1.-random_number);
      p[2][1] = 0.70*(1.+random_number);
      p[3][1] = 0.68*(1.-random_number);
      p[4][1] = 0.68*(1.-random_number);
      p[5][1] = 0.68*(1.-random_number);

      // omega_cdm
      p[0][2] = 0.11*(1.-random_number);
      p[1][2] = 0.11*(1.-random_number);
      p[2][2] = 0.11*(1.-random_number);
      p[3][2] = 0.115*(1.+random_number);
      p[4][2] = 0.11*(1.-random_number);
      p[5][2] = 0.11*(1.-random_number);

      // ns
      p[0][3] = 0.965*(1.-random_number);
      p[1][3] = 0.965*(1.-random_number);
      p[2][3] = 0.965*(1.-random_number);
      p[3][3] = 0.965*(1.-random_number);
      p[4][3] = 0.97*(1.+random_number);
      p[5][3] = 0.965*(1.-random_number);

      // As
      p[0][4] = 2.25e-9*(1.-random_number);
      p[1][4] = 2.25e-9*(1.-random_number);
      p[2][4] = 2.25e-9*(1.-random_number);
      p[3][4] = 2.25e-9*(1.-random_number);
      p[4][4] = 2.25e-9*(1.-random_number);
      p[5][4] = 2.35e-9*(1.+random_number);

      /*sprintf(fc.value[5],"%g",h*sqrt(alpha));*/
      /*sprintf(fc.value[6],"%g",(omega_cdm+omega_b)*alpha-omega_b);*/
      /*sprintf(fc.value[7],"%g",omega_b);*/


      // Initialization of the simplex method, compute ndim+1 points.
      for (k=0; k<ndim+1; k++){
	sprintf(fc.value[4],"%g",p[k][0]);
	sprintf(fc.value[5],"%g",p[k][1]);
	sprintf(fc.value[6],"%g",p[k][2]);
	sprintf(fc.value[8],"%g",p[k][3]);
	sprintf(fc.value[9],"%g",p[k][4]);
	if (input_init(&fc,&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&nl,&le,&op,errmsg) == _FAILURE_) {
	  printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg);
	  return _FAILURE_;
	}

	class_assuming_bessels_computed(&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&nl,&le,&op,cl,errmsg);
	chi2_planck(&sp,cl,cl_ref,noise,&chi2[k]);

      }
      amoeba(p,chi2,ndim,ftol,get_chi2);
      if (j==0){
	storage = (chi2[0]+chi2[1]+chi2[2]+chi2[3]+chi2[4]+chi2[5])/6.;
      }
      fprintf(stderr,"run %d, chi2:%f\n",j,(chi2[0]+chi2[1]+chi2[2]+chi2[3]+chi2[4]+chi2[5])/6.);
      if ((chi2[0]+chi2[1]+chi2[2]+chi2[3]+chi2[4]+chi2[5])/6. < storage){
	storage = (chi2[0]+chi2[1]+chi2[2]+chi2[3]+chi2[4]+chi2[5])/6.;
      }
      if (j==4){
	fprintf(stderr,"final chi2: %f\n",storage);
      }
	// Initialization of the simplex method, compute ndim+1 points.
    }
    

    

    
    fprintf(stdout,"%e %e %e %e %e %e %e\n",parameter[i],(p[0][0]+p[1][0]+p[2][0]+p[3][0]+p[4][0]+p[5][0])/6.,
	(p[0][1]+p[1][1]+p[2][1]+p[3][1]+p[4][1]+p[5][1])/6.,
	(p[0][2]+p[1][2]+p[2][2]+p[3][2]+p[4][1]+p[5][1])/6.,
	(p[0][3]+p[1][3]+p[2][3]+p[3][3]+p[4][3]+p[5][3])/6.,
	(p[0][4]+p[1][4]+p[2][4]+p[3][4]+p[4][4]+p[5][4])/6.,
	storage);
  }
    /*fprintf(stdout,"%e %e %e %e %e\n",parameter[i],(p[0][0]+p[1][0]+p[2][0]+p[3][0])/4.,*/
	/*(p[0][1]+p[1][1]+p[2][1]+p[3][1])/4.,*/
	/*(p[0][2]+p[1][2]+p[2][2]+p[3][2])/4.,*/
	/*chi2[0]);*/
    /*fprintf(stdout,"%e %e %e\n",parameter[i],(p[0][0]+p[1][0])/2.,chi2[0]);*/
    return _SUCCESS_;

  }


int class_assuming_bessels_computed(
				    struct precision * ppr,
				    struct background * pba,
				    struct thermo * pth,
				    struct perturbs * ppt,
				    struct bessels * pbs,
				    struct transfers * ptr,
				    struct primordial * ppm,
				    struct spectra * psp,
				    struct nonlinear * pnl,
				    struct lensing * ple,
				    struct output * pop,
				    double ** cl,
				    ErrorMsg errmsg) {
  
  int l;
  double ** junk1;
  double ** junk2;

  if (background_init(ppr,pba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",pba->error_message);
    return _FAILURE_;
  }
  /*fprintf(stderr,"h = %e, omga_c = %e, omega_b = %e, Omega_Lambda = %e, Omega_ur = %e, z_eq = %e\n",pba->h,pba->Omega0_cdm*pba->h*pba->h,pba->Omega0_b*pba->h*pba->h,pba->Omega0_lambda,pba->Omega0_ur,(pba->Omega0_b+pba->Omega0_cdm)/(pba->Omega0_g+pba->Omega0_ur));*/
    
  if (thermodynamics_init(ppr,pba,pth) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",pth->error_message);
    return _FAILURE_;
  }

  if (perturb_init(ppr,pba,pth,ppt) == _FAILURE_) {
    printf("\n\nError in perturb_init \n=>%s\n",ppt->error_message);
    return _FAILURE_;
  }

  if (bessel_init(&pr,&bs) == _FAILURE_) {
    printf("\n\nError in bessel_init \n =>%s\n",bs.error_message);
    return _FAILURE_;
  }

  if (transfer_init(ppr,pba,pth,ppt,pbs,ptr) == _FAILURE_) {
    printf("\n\nError in transfer_init \n=>%s\n",ptr->error_message);
    return _FAILURE_;
  }

  if (primordial_init(ppr,ppt,ppm) == _FAILURE_) {
    printf("\n\nError in primordial_init \n=>%s\n",ppm->error_message);
    return _FAILURE_;
  }

  if (spectra_init(ppr,pba,ppt,ptr,ppm,psp) == _FAILURE_) {
    printf("\n\nError in spectra_init \n=>%s\n",psp->error_message);
    return _FAILURE_;
  }

  if (nonlinear_init(ppr,pba,pth,ppt,pbs,ptr,ppm,psp,pnl) == _FAILURE_) {
    printf("\n\nError in nonlinear_init \n=>%s\n",psp->error_message);
    return _FAILURE_;
  }

  if (lensing_init(ppr,ppt,psp,pnl,ple) == _FAILURE_) {
    printf("\n\nError in lensing_init \n=>%s\n",ple->error_message);
    return _FAILURE_;
  }

  for (l=2; l <= l_max; l++) {

    if (output_total_cl_at_l(psp,ple,pop,(double)l,cl[l]) == _FAILURE_) {
      printf("\n\nError in spectra_cl_at_l \n=>%s\n",psp->error_message);
      return _FAILURE_;}
    }
  

  /****** done ******/

  if (lensing_free(ple) == _FAILURE_) {
    printf("\n\nError in spectra_free \n=>%s\n",ple->error_message);
    return _FAILURE_;
  }

  if (nonlinear_free(pnl) == _FAILURE_) {
    printf("\n\nError in nonlinear_free \n=>%s\n",psp->error_message);
    return _FAILURE_;
  }
    
  if (spectra_free(psp) == _FAILURE_) {
    printf("\n\nError in spectra_free \n=>%s\n",psp->error_message);
    return _FAILURE_;
  }
    
  if (primordial_free(ppm) == _FAILURE_) {
    printf("\n\nError in primordial_free \n=>%s\n",ppm->error_message);
    return _FAILURE_;
  }
  
  if (transfer_free(ptr) == _FAILURE_) {
    printf("\n\nError in transfer_free \n=>%s\n",ptr->error_message);
    return _FAILURE_;
  }

  if (bessel_free(&bs) == _FAILURE_)  {
    printf("\n\nError in bessel_free \n=>%s\n",bs.error_message);
    return _FAILURE_;
  }

  if (perturb_free(ppt) == _FAILURE_) {
    printf("\n\nError in perturb_free \n=>%s\n",ppt->error_message);
    return _FAILURE_;
  }

  if (thermodynamics_free(pth) == _FAILURE_) {
    printf("\n\nError in thermodynamics_free \n=>%s\n",pth->error_message);
    return _FAILURE_;
  }

  if (background_free(pba) == _FAILURE_) {
    printf("\n\nError in background_free \n=>%s\n",pba->error_message);
    return _FAILURE_;
  }

  return _SUCCESS_;

}

int chi2_simple(
		struct spectra * psp,
		double ** cl1,
		double ** cl2,
		double * chi2) {

  int l;

  *chi2=0.;
  for (l=2; l <= l_max; l++) {
    *chi2 += 
      pow(((cl1[l][psp->index_ct_tt]-cl2[l][psp->index_ct_tt])
	   /cl2[l][psp->index_ct_tt]),2) +
      pow(((cl1[l][psp->index_ct_te]-cl2[l][psp->index_ct_te])
	   /cl2[l][psp->index_ct_te]),2) +
      pow(((cl1[l][psp->index_ct_ee]-cl2[l][psp->index_ct_ee])
	   /cl2[l][psp->index_ct_ee]),2);
  }

  return _SUCCESS_;

}

int chi2_planck(
		struct spectra * psp,
		double ** cl1, /* treated as 'theoretical spectrum' */
		double ** cl2, /* treated as 'observed spectrum' */
		double ** noise,
		double * chi2) {

  int l;
  double fsky=0.65;
  double clTT_th,clEE_th,clTE_th;
  double clTT_obs,clEE_obs,clTE_obs;
  double det_mixed,det_th,det_obs;
  double factor;

  *chi2=0.;
  for (l=2; l <= l_max; l++) {

/*     if (psp->ct_size == 1) { */
    if (0==1) {

      *chi2 += fsky*(2.*l+1.)*((cl2[l][0]+noise[l][0])/
			       (cl1[l][0]+noise[l][0])+
			       log((cl1[l][0]+noise[l][0])/(cl2[l][0]+noise[l][0]))-1.);

/*       *chi2 += fsky*(2.*l+1.)*((cl2[l][0])/ */
/* 			       (cl1[l][0])+ */
/* 			       log((cl1[l][0])/(cl2[l][0]))-1.); */

      
/*       factor=l*(l+1.)/2./_PI_; */
/*       printf("%d %e %e %e\n",l,factor*cl1[l][0],factor*cl2[l][0],factor*noise[l][0]); */

    }
    else {

      clTT_th = cl1[l][psp->index_ct_tt]+noise[l][psp->index_ct_tt];
      clEE_th = cl1[l][psp->index_ct_ee]+noise[l][psp->index_ct_ee];
      clTE_th = cl1[l][psp->index_ct_te];
      
      clTT_obs = cl2[l][psp->index_ct_tt]+noise[l][psp->index_ct_tt];
      clEE_obs = cl2[l][psp->index_ct_ee]+noise[l][psp->index_ct_ee];
      clTE_obs = cl2[l][psp->index_ct_te];

      det_mixed = 0.5*(clTT_th*clEE_obs+clTT_obs*clEE_th)-clTE_th*clTE_obs;

      det_th = clTT_th*clEE_th-clTE_th*clTE_th;

      det_obs = clTT_obs*clEE_obs-clTE_obs*clTE_obs;

      *chi2 += fsky*(2.*l+1.)*(2.*(det_mixed/det_th-1.)+log(det_th/det_obs));
    }
  }

  return _SUCCESS_;

}

int noise_planck(
		 struct background * pba,
		 struct thermo * pth,
		 struct spectra * psp,
		 double ** noise
		 ) {

  int num_channel,channel,l;
  double * theta, * deltaT, * deltaP; 

  num_channel=3;

  theta=malloc(num_channel*sizeof(double));
  deltaT=malloc(num_channel*sizeof(double));
  deltaP=malloc(num_channel*sizeof(double));
  
  /* 100 GHz*/
  theta[0]=9.5/3437.75;          /* converted from arcmin to radian */
  deltaT[0]=(6.8e-6/pba->T_cmb);  /* converted from K to dimensionless */
  deltaP[0]=(10.9e-6/pba->T_cmb); /* converted from K to dimensionless */
  
  /* 143 GHz*/
  theta[1]=7.1/3437.75;
  deltaT[1]=(6.0e-6/pba->T_cmb);
  deltaP[1]=(11.4e-6/pba->T_cmb);
  
  /* 217 GHz*/
  theta[2]=5.0/3437.75;
  deltaT[2]=(13.1e-6/pba->T_cmb);
  deltaP[2]=(26.7e-6/pba->T_cmb);

  for (l=2; l<=l_max; l++) {

    for (channel=0; channel<num_channel; channel++) {

      if (psp->has_tt) {
	noise[l][psp->index_ct_tt] +=
	  pow((exp(-l*(l+1.)*theta[channel]*theta[channel]/16./log(2.))
	       /theta[channel]/deltaT[channel]),2);
      }
      if (psp->has_ee) {
	noise[l][psp->index_ct_ee] +=
	  pow((exp(-l*(l+1.)*theta[channel]*theta[channel]/16./log(2.))
	       /theta[channel]/deltaP[channel]),2);
      }
    }

    if (psp->has_tt) {
      noise[l][psp->index_ct_tt] = 1./noise[l][psp->index_ct_tt];
    }
    if (psp->has_ee) {
    noise[l][psp->index_ct_ee] = 1./noise[l][psp->index_ct_ee];
    }
  }

  free(theta);
  free(deltaT);
  free(deltaP);

  return _SUCCESS_;
  
}

