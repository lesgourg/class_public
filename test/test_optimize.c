/** @file class.c 
 * Julien Lesgourgues, 20.04.2010    
 */
 
#include "class.h"

int class(
	  struct precision * ppr,
	  struct background * pba,
	  struct thermo * pth,
	  struct perturbs * ppt,
	  struct bessels * pbs,
	  struct transfers * ptr,
	  struct primordial * ppm,
	  struct spectra * psp,
	  struct nonlinear * pnl,
	  struct lensing *ple,
	  struct output * pop,
	  int l_max,
	  double ** cl,
	  ErrorMsg errmsg);

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
				    struct lensing *ple,
				    struct output * pop,
				    int l_max,
				    double ** cl,
				    ErrorMsg errmsg);

int chi2_simple(
		struct spectra * psp,
		double ** cl1,
		double ** cl2,
		int lmax,
		double * chi2);

int chi2_planck(
		struct spectra * psp,
		double ** cl1,
		double ** cl2,
		double ** nl,
		int lmax,
		double * chi2);

int noise_planck(
		 struct background * pba,
		 struct thermo * pth,
		 struct spectra * psp,
		 double ** nl,
		 int lmax
		 );

int main(int argc, char **argv) {

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions */
  struct bessels bs;          /* for bessel functions */
  struct transfers tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct spectra sp;          /* for output spectra */
  struct nonlinear nl;        /* for non-linear spectra */
  struct lensing le;          /* for lensed spectra */
  struct output op;           /* for output files */

  ErrorMsg errmsg;

  int i,l_max,l;

  double parameter_initial,parameter_logstep;
  double * parameter;
  int param_num;

  double *** cl;
  double ** noise;
  double chi2,chi2_bis;
  double percentage,max_percentage;
  int max_l;
  int ref_run;

  double * param;

  FILE * output;
  FILE * results;
  char filename[60];
  char results_name[60];
  char junk_string[120];
  int l_read;

/*******************************/

  param = &(pr.perturb_sampling_stepsize);
  //param = &(pr.radiation_streaming_trigger_tau_over_tau_k);
  //param = &(pr.radiation_streaming_trigger_tau_c_over_tau);

  parameter_initial=0.15;
  parameter_logstep=1.1;

  param_num=10;
  ref_run=-1;

  /* if ref_run<0, the reference is taken in the following external file: */

  sprintf(filename,"output/cl_ref_cl_lensed.dat");

/*******************************************************/

  if (input_init_from_arguments(argc,argv,&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&nl,&le,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }

  //l_max = pt.l_scalar_max;
  l_max=3000;

  class_alloc(parameter,param_num*sizeof(double),errmsg);

  class_calloc(noise,(l_max+1),sizeof(double*),errmsg);
  for (l=2; l <= l_max; l++) {
    class_calloc(noise[l],2,sizeof(double),errmsg);
  }


  if (ref_run >= param_num) {
    fprintf(stderr,"ref_run=%d out of allowed range\n",ref_run);
    return _FAILURE_;
  }

  else if (ref_run >= 0) {
    class_alloc(cl,param_num*sizeof(double**),errmsg);
    for (i=0; i<param_num; i++) {
      class_alloc(cl[i],(l_max+1)*sizeof(double*),errmsg);
      for (l=2; l <= l_max; l++) {
	class_alloc(cl[i][l],6*sizeof(double),errmsg);
      }
    }
  }

  else {
    class_alloc(cl,(param_num+1)*sizeof(double**),errmsg);
    for (i=0; i<(param_num+1); i++) {
      class_alloc(cl[i],(l_max+1)*sizeof(double*),errmsg);
      for (l=2; l <= l_max; l++) {
        class_alloc(cl[i][l],6*sizeof(double),errmsg);
      }
    }

    ref_run = param_num;

    /* read file and fill cl[param_num] */
    output = fopen(filename,"r");
    fprintf(stderr,"Read reference Cls in file %s\n",filename);
    fgets(junk_string,120,output);
    //    fprintf(stderr,"%s\n",junk_string);
    fgets(junk_string,120,output);
    //    fprintf(stderr,"%s\n",junk_string);
    fgets(junk_string,120,output);
    //    fprintf(stderr,"%s\n",junk_string);
    fgets(junk_string,120,output);
    //    fprintf(stderr,"%s\n",junk_string);
    float cl_read;
    for (l=2; l <= l_max; l++) {
      fscanf(output,"%d",&l_read);
      //fprintf(stderr,"%d",l_read);

      fscanf(output,"%e",&cl_read);
      //fprintf(stderr," %e",cl_read);
      cl[ref_run][l][0]=(double)cl_read*2.*_PI_/l/(l+1);
      //fprintf(stderr,"%d %e\n",l_read,cl[ref_run][l][sp.index_ct_tt]);

      fscanf(output,"%e",&cl_read);
      //fprintf(stderr," %e",cl_read);
      cl[ref_run][l][1]=(double)cl_read*2.*_PI_/l/(l+1);

      fscanf(output,"%e",&cl_read);
      cl[ref_run][l][2]=(double)cl_read*2.*_PI_/l/(l+1);
      //fprintf(stderr," %e\n",cl_read);

      fscanf(output,"%e",&cl_read);
      fscanf(output,"%e",&cl_read);
      fscanf(output,"%e",&cl_read);
      if (l_read != l) {
	printf("l_read != l: %d %d\n",l_read,l);
      }
    }

    fclose(output);

  }


  if (bessel_init(&pr,&bs) == _FAILURE_) {
    printf("\n\nError in bessel_init \n =>%s\n",bs.error_message);
    return _FAILURE_;
  }

  for (i=0; i<param_num; i++) {

    parameter[i] = parameter_initial * exp((double)i*log(parameter_logstep));
    
    *param = parameter[i];

    fprintf(stderr,"#run %d/%d with %g=%g\n",i+1,param_num,parameter[i],*param);

    //class(&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&nl,&le,&op,l_max,cl[i],errmsg);
    class_assuming_bessels_computed(&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&nl,&le,&op,l_max,cl[i],errmsg);
    
  }

  if (bessel_free(&bs) == _FAILURE_)  {
    printf("\n\nError in bessel_free \n=>%s\n",bs.error_message);
    return _FAILURE_;
  }

  noise_planck(&ba,&th,&sp,noise,l_max);

  for (i=0; i<param_num; i++) {
   
    chi2_planck(&sp,cl[i],cl[ref_run],noise,l_max,&chi2);

    if (chi2>10.) {
      fprintf(stderr,"parameter=%e BAD: chi2=%e \n",
	      parameter[i],chi2);
    }
    else {
      fprintf(stderr,"parameter=%e OK : chi2=%e \n",
	      parameter[i],chi2);
    }
  }

  return _SUCCESS_;

}


int class(
	  struct precision * ppr,
	  struct background * pba,
	  struct thermo * pth,
	  struct perturbs * ppt,
	  struct bessels * pbs,
	  struct transfers * ptr,
	  struct primordial * ppm,
	  struct spectra * psp,
	  struct nonlinear * pnl,
	  struct lensing *ple,
	  struct output * pop,
	  int l_max,
	  double ** cl,
	  ErrorMsg errmsg) {

  int l;
  double ** junk1;
  double ** junk2;

  if (background_init(ppr,pba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",pba->error_message);
    return _FAILURE_;
  }
    
  if (thermodynamics_init(ppr,pba,pth) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",pth->error_message);
    return _FAILURE_;
  }

  if (perturb_init(ppr,pba,pth,ppt) == _FAILURE_) {
    printf("\n\nError in perturb_init \n=>%s\n",ppt->error_message);
    return _FAILURE_;
  }

  if (bessel_init(ppr,pbs) == _FAILURE_) {
    printf("\n\nError in bessel_init \n =>%s\n",pbs->error_message);
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

  if (nonlinear_init(ppr,pba,pth,ppm,psp,pnl) == _FAILURE_) {
    printf("\n\nError in nonlinear_init \n=>%s\n",pnl->error_message);
    return _FAILURE_;
  }

  if (lensing_init(ppr,ppt,psp,pnl,ple) == _FAILURE_) {
    printf("\n\nError in lensing_init \n=>%s\n",ple->error_message);
    return _FAILURE_;
  }

  for (l=2; l <= l_max; l++) {

    if (output_total_cl_at_l(psp,ple,pop,(double)l,cl[l]) == _FAILURE_) {
      printf("\n\nError in spectra_cl_at_l \n=>%s\n",psp->error_message);
      return _FAILURE_;
    }

  }

  /****** done ******/

  if (lensing_free(ple) == _FAILURE_) {
    printf("\n\nError in spectra_free \n=>%s\n",ple->error_message);
    return _FAILURE_;
  }

  if (nonlinear_free(pnl) == _FAILURE_) {
    printf("\n\nError in nonlinear_free \n=>%s\n",pnl->error_message);
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

  if (bessel_free(pbs) == _FAILURE_)  {
    printf("\n\nError in bessel_free \n=>%s\n",pbs->error_message);
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
				    int l_max,
				    double ** cl,
				    ErrorMsg errmsg) {
  
  int l;
  double ** junk1;
  double ** junk2;

  if (background_init(ppr,pba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",pba->error_message);
    return _FAILURE_;
  }
    
  if (thermodynamics_init(ppr,pba,pth) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",pth->error_message);
    return _FAILURE_;
  }

  if (perturb_init(ppr,pba,pth,ppt) == _FAILURE_) {
    printf("\n\nError in perturb_init \n=>%s\n",ppt->error_message);
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

  if (nonlinear_init(ppr,pba,pth,ppm,psp,pnl) == _FAILURE_) {
    printf("\n\nError in nonlinear_init \n=>%s\n",pnl->error_message);
    return _FAILURE_;
  }

  if (lensing_init(ppr,ppt,psp,pnl,ple) == _FAILURE_) {
    printf("\n\nError in lensing_init \n=>%s\n",ple->error_message);
    return _FAILURE_;
  }

  for (l=2; l <= l_max; l++) {

    if (output_total_cl_at_l(psp,ple,pop,(double)l,cl[l]) == _FAILURE_) {
      printf("\n\nError in spectra_cl_at_l \n=>%s\n",psp->error_message);
      return _FAILURE_;
    }

    //fprintf(stderr,"%d %e %e %e\n",l,cl[l][0],cl[l][1],cl[l][2]);

  }

  /****** done ******/

  if (lensing_free(ple) == _FAILURE_) {
    printf("\n\nError in spectra_free \n=>%s\n",ple->error_message);
    return _FAILURE_;
  }

  if (nonlinear_free(pnl) == _FAILURE_) {
    printf("\n\nError in nonlinear_free \n=>%s\n",pnl->error_message);
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
		int lmax,
		double * chi2) {

  int l;

  *chi2=0.;
  for (l=2; l <= lmax; l++) {
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
		double ** nl,
		int lmax,
		double * chi2) {

  int l;
  double fsky=0.8;
  double clTT_th,clEE_th,clTE_th;
  double clTT_obs,clEE_obs,clTE_obs;
  double det_mixed,det_th,det_obs;
  double factor;

  *chi2=0.;
  for (l=2; l <= lmax; l++) {

/*     if (psp->ct_size == 1) { */
    if (0==1) {

      *chi2 += fsky*(2.*l+1.)*((cl2[l][0]+nl[l][0])/
			       (cl1[l][0]+nl[l][0])+
			       log((cl1[l][0]+nl[l][0])/(cl2[l][0]+nl[l][0]))-1.);

      //fprintf(stderr,"%d %e %e %e\n",l,cl1[l][0],cl1[l][1],cl1[l][2]);
      //fprintf(stderr,"%d %e %e %e\n",l,cl2[l][0],cl2[l][1],cl2[l][2]);
      //fprintf(stderr,"%d %e %e\n",l,nl[l][0],nl[l][1]);

/*       *chi2 += fsky*(2.*l+1.)*((cl2[l][0])/ */
/* 			       (cl1[l][0])+ */
/* 			       log((cl1[l][0])/(cl2[l][0]))-1.); */

      
/*       factor=l*(l+1.)/2./_PI_; */
/*       printf("%d %e %e %e\n",l,factor*cl1[l][0],factor*cl2[l][0],factor*nl[l][0]); */

    }
    else {

      clTT_th = cl1[l][psp->index_ct_tt]+nl[l][psp->index_ct_tt];
      clEE_th = cl1[l][psp->index_ct_ee]+nl[l][psp->index_ct_ee];
      clTE_th = cl1[l][psp->index_ct_te];
      
      clTT_obs = cl2[l][psp->index_ct_tt]+nl[l][psp->index_ct_tt];
      clEE_obs = cl2[l][psp->index_ct_ee]+nl[l][psp->index_ct_ee];
      clTE_obs = cl2[l][psp->index_ct_te];

      //printf("%e %e %e %e %e %e\n",clTT_th,clTT_obs,clEE_th,clEE_obs,clTE_th,clTE_obs);

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
		 double ** nl,
		 int lmax
		 ) {

  int num_channel,channel,l;
  double * theta, * deltaT, * deltaP; 

  num_channel=3;

  theta=malloc(num_channel*sizeof(double));
  deltaT=malloc(num_channel*sizeof(double));
  deltaP=malloc(num_channel*sizeof(double));
  
  /* 100 GHz*/
  theta[0]=9.5/3437.75;          /* converted from arcmin to radian */
  deltaT[0]=(6.8e-6/pba->Tcmb);  /* converted from K to dimensionless */
  deltaP[0]=(10.9e-6/pba->Tcmb); /* converted from K to dimensionless */
  
  /* 143 GHz*/
  theta[1]=7.1/3437.75;
  deltaT[1]=(6.0e-6/pba->Tcmb);
  deltaP[1]=(11.4e-6/pba->Tcmb);
  
  /* 217 GHz*/
  theta[2]=5.0/3437.75;
  deltaT[2]=(13.1e-6/pba->Tcmb);
  deltaP[2]=(26.7e-6/pba->Tcmb);

  for (l=2; l<=lmax; l++) {

    for (channel=0; channel<num_channel; channel++) {

      if (psp->has_tt) {
	nl[l][0] +=
	  pow((exp(-l*(l+1.)*theta[channel]*theta[channel]/16./log(2.))
	       /theta[channel]/deltaT[channel]),2);
      }
      if (psp->has_ee) {
	nl[l][1] +=
	  pow((exp(-l*(l+1.)*theta[channel]*theta[channel]/16./log(2.))
	       /theta[channel]/deltaP[channel]),2);
      }
    }

    if (psp->has_tt) {
      nl[l][0] = 1./nl[l][0];
    }
    if (psp->has_ee) {
    nl[l][1] = 1./nl[l][1];
    }

    //fprintf(stderr,"%d %e %e\n",l,nl[l][0],nl[l][1]);

  }

  free(theta);
  free(deltaT);
  free(deltaP);

  return _SUCCESS_;
  
}
