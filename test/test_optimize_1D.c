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
  struct lensing le;          /* for lensed spectra */
  struct output op;           /* for output files */
  struct spectra_nl nl;       /* for calculation of non-linear spectra */

  ErrorMsg errmsg;

  int i,l_max,l;

  struct file_content fc;

  double parameter_initial,parameter_logstep;
  double * parameter;
  int param_num;

  double *** cl;
  double ** noise;
  double chi2,chi2_bis;
  double percentage,max_percentage;
  int max_l;
  int ref_run;
  
  l_max=2500;

  parser_init(&fc,4,errmsg);

  strcpy(fc.name[0],"output");
  strcpy(fc.value[0],"tCl,pCl");

  strcpy(fc.name[1],"l_max_scalars");
  sprintf(fc.value[1],"%d",l_max);
 
/*   strcpy(fc.name[2],"perturbations_verbose"); */
/*   sprintf(fc.value[2],"%d",2); */

/*******************************************************/

  strcpy(fc.name[2],"lensing");
  strcpy(fc.value[2],"no");

  strcpy(fc.name[3],"transfer_cut_threshold_cl");

  parameter_initial=1.e-7;
  parameter_logstep=1.3;

  param_num=30;
  ref_run=0;

/*******************************************************/

  class_alloc(cl,param_num*sizeof(double**),errmsg);
  class_alloc(parameter,param_num*sizeof(double),errmsg);
  for (i=0; i<param_num; i++) {
    class_alloc(cl[i],(l_max+1)*sizeof(double*),errmsg);
    class_calloc(noise,(l_max+1),sizeof(double*),errmsg);
    for (l=2; l <= l_max; l++) {
      class_alloc(cl[i][l],3*sizeof(double),errmsg);
      class_calloc(noise[l],3,sizeof(double),errmsg);
    }
  }

  for (i=0; i<param_num; i++) {

    parameter[i] = parameter_initial * exp((double)i*log(parameter_logstep));

/*       parameter[i] = parameter_initial -i; */

/*     if (i==0) { */
/*       sprintf(fc.value[2],"%d",tc_osc); */
/*       strcpy(fc.name[3],"transfer_cut_threshold_osc"); */
/*       sprintf(fc.value[3],"%g",0.013); */
/*     } */
/*     else  { */
/*       sprintf(fc.value[2],"%d",tc_cl); */
/*       strcpy(fc.name[3],"transfer_cut_threshold_cl"); */
/*       sprintf(fc.value[3],"%g",8.e-7); */
/*     } */

/*     sprintf(fc.value[2],"%g",parameter[i]); */
      sprintf(fc.value[3],"%g",parameter[i]);
/*     sprintf(fc.value[3],"%d",(int)parameter[i]); */
 /*    sprintf(fc.value[2],"%d",1); */

    fprintf(stderr,"#run %d/%d with %s\n",i+1,param_num,fc.value[3]);

    if (input_init(&fc,&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&le,&op,&nl,errmsg) == _FAILURE_) {
      printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
      return _FAILURE_;
    }

    if (i==0) {
      
      if (bessel_init(&pr,&bs) == _FAILURE_) {
	printf("\n\nError in bessel_init \n =>%s\n",bs.error_message);
	return _FAILURE_;
      }
    }


 /*    class(&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&le,&op,l_max,cl[i],errmsg); */
    class_assuming_bessels_computed(&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&le,&op,l_max,cl[i],errmsg);
    

/*     for (l=2; l <= 2500; l++) */
/*       printf("%d %e %e %e\n",l, */
    /* 	     (l*(l+1)/2./_PI_)*cl[i][l][0],  */ /* here cl is the dimensionless C_l */ 
    /* 	     (l*(l+1)/2./_PI_)*cl[i][l][1],  */ /* multiply by pow((th.Tcmb*1.e6),2) for muK */ 
/* 	     (l*(l+1)/2./_PI_)*cl[i][l][2]); */
      
/*     printf("\n"); */
    
  }

  if (bessel_free(&bs) == _FAILURE_)  {
    printf("\n\nError in bessel_free \n=>%s\n",bs.error_message);
    return _FAILURE_;
  }

/*   for (i=0; i<param_num-1; i++) { */

/*     chi2=0; */
/*     max_percentage = 0.; */
/*     max_l = 0; */
/*     for (l=2; l <= 2500; l++) { */
/*       chi2 += pow(((cl[i][l][0]-cl[param_num-1][l][0])/cl[param_num-1][l][0]),2); */
/*       percentage = fabs(cl[i][l][0]/cl[param_num-1][l][0]-1.)*100.; */
/*       if (percentage > max_percentage) { */
/* 	max_percentage = percentage; */
/* 	max_l = l; */
/*       } */
/*     } */

/*     chi2_simple(cl[i],cl[param_num-1],2500,&chi2_bis); */

/*     fprintf(stderr,"parameter=%e chi2=%e (%e) l=%d percentage=%g\n", */
/* 	    parameter[i],chi2,chi2_bis,max_l,max_percentage); */
/*   } */

  noise_planck(&ba,&th,&sp,noise,l_max);

/*   for (l=2;l<=l_max;l++) { */
/*     printf("%d  %e  %e  %e  %e  %e  %e\n",l, */
/* 	   (l*(l+1)/2./_PI_)*cl[0][l][0]*pow(th.Tcmb*1.e6,2), */
/* 	   (l*(l+1)/2./_PI_)*cl[0][l][1]*pow(th.Tcmb*1.e6,2), */
/* 	   (l*(l+1)/2./_PI_)*cl[0][l][2]*pow(th.Tcmb*1.e6,2), */
/* 	   (l*(l+1)/2./_PI_)*noise[l][0]*pow(th.Tcmb*1.e6,2), */
/* 	   (l*(l+1)/2./_PI_)*noise[l][1]*pow(th.Tcmb*1.e6,2), */
/* 	   (l*(l+1)/2./_PI_)*noise[l][2]*pow(th.Tcmb*1.e6,2)); */
/*   } */

  for (i=0; i<param_num; i++) {
   
    chi2_planck(&sp,cl[i],cl[ref_run],noise,l_max,&chi2);

    if (chi2>0.1) {
      fprintf(stderr,"parameter=%e BAD: chi2=%2g \n",
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

  if (spectra_init(pba,ppt,ptr,ppm,psp) == _FAILURE_) {
    printf("\n\nError in spectra_init \n=>%s\n",psp->error_message);
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

  if (spectra_init(pba,ppt,ptr,ppm,psp) == _FAILURE_) {
    printf("\n\nError in spectra_init \n=>%s\n",psp->error_message);
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
	nl[l][psp->index_ct_tt] +=
	  pow((exp(-l*(l+1.)*theta[channel]*theta[channel]/16./log(2.))
	       /theta[channel]/deltaT[channel]),2);
      }
      if (psp->has_ee) {
	nl[l][psp->index_ct_ee] +=
	  pow((exp(-l*(l+1.)*theta[channel]*theta[channel]/16./log(2.))
	       /theta[channel]/deltaP[channel]),2);
      }
    }

    if (psp->has_tt) {
      nl[l][psp->index_ct_tt] = 1./nl[l][psp->index_ct_tt];
    }
    if (psp->has_ee) {
    nl[l][psp->index_ct_ee] = 1./nl[l][psp->index_ct_ee];
    }
  }

  free(theta);
  free(deltaT);
  free(deltaP);

  return _SUCCESS_;
  
}
