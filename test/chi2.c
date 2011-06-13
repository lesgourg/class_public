/** @file class.c 
 * Julien Lesgourgues, 20.04.2010    
 */

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "float.h"

#define _PI_ 3.14159

int chi2_simple(
		double ** cl1,
		double ** cl2,
		int lmax,
		double * chi2);

int chi2_planck(
		double ** cl1,
		double ** cl2,
		double ** nl,
		int lmax,
		double * chi2);

int noise_planck(
		 double ** nl,
		 int lmax
		 );

int main(int argc, char **argv) {

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
  char filename[200];
  char results_name[200];
  char junk_string[200];
  int l_read;

  /*******************************/

  l_max = 3000;
  
  cl=malloc(2*sizeof(double**));
  for (i=0; i<2; i++) {
    cl[i]=malloc((l_max+1)*sizeof(double*));
    for (l=2; l <= l_max; l++) {
      cl[i][l]=malloc(3*sizeof(double));
    }
  }
  
  noise=calloc((l_max+1),sizeof(double*));
  for (l=2; l <= l_max; l++) {
    noise[l]=calloc(3,sizeof(double));
  }

  sprintf(filename,"output/chi2pl0.1_cl_lensed.dat");

  /* read file and fill cl[param_num] */
  output = fopen(filename,"r");
  fprintf(stderr,"Read reference in file %s\n",filename);
  fgets(junk_string,200,output);
  //    fprintf(stderr,"%s\n",junk_string);
  fgets(junk_string,200,output);
  //    fprintf(stderr,"%s\n",junk_string);
  fgets(junk_string,200,output);
  //    fprintf(stderr,"%s\n",junk_string);
  fgets(junk_string,200,output);
  fgets(junk_string,200,output);
  fgets(junk_string,200,output);
  fgets(junk_string,200,output);
  //    fprintf(stderr,"%s\n",junk_string);
  //fgets(junk_string,200,output);
  //fgets(junk_string,200,output);
  //fgets(junk_string,200,output);
  float cl_read;
  for (l=2; l <= l_max; l++) {
    fscanf(output,"%d",&l_read);
    //fprintf(stderr,"%d\n",l_read);
    fscanf(output,"%e",&cl_read);
    //fprintf(stderr," %e\n",cl_read);
    cl[0][l][0]=(double)cl_read*2.*_PI_/l/(l+1);
    //      fprintf(stderr,"%d %e\n",l_read,cl[ref_run][l][sp.index_ct_tt]);
    fscanf(output,"%e",&cl_read);
    cl[0][l][1]=(double)cl_read*2.*_PI_/l/(l+1);
    fscanf(output,"%e",&cl_read);
    cl[0][l][2]=(double)cl_read*2.*_PI_/l/(l+1);
    //fprintf(stderr," %e %d\n",cl_read,sp.index_ct_te);
    fscanf(output,"%e",&cl_read);
    fscanf(output,"%e",&cl_read);
    fscanf(output,"%e",&cl_read);
    fscanf(output,"%e",&cl_read);
    if (l_read != l) {
      printf("l_read != l: %d %d\n",l_read,l);
    }
  }
  
  fclose(output);
    
  sprintf(filename,"output/cl_ref_cl_lensed.dat");
  
  /* read file and fill cl[param_num] */
  output = fopen(filename,"r");
  fprintf(stderr,"Read degraded precision in file %s\n",filename);
  fgets(junk_string,200,output);
  //fprintf(stderr,"%s\n",junk_string);
  fgets(junk_string,200,output);
  //    fprintf(stderr,"%s\n",junk_string);
  fgets(junk_string,200,output);
  //    fprintf(stderr,"%s\n",junk_string);
  fgets(junk_string,200,output);
  //    fprintf(stderr,"%s\n",junk_string);
  fgets(junk_string,200,output);
  fgets(junk_string,200,output);
  fgets(junk_string,200,output);
  //float cl_read;
  //fprintf(stderr,"get here\n");
  for (l=2; l <= l_max; l++) {
    fscanf(output,"%d",&l_read);
    //      fprintf(stderr,"%d",l_read);
    fscanf(output,"%e",&cl_read);
    //      fprintf(stderr," %e\n",cl_read);
    cl[1][l][0]=(double)cl_read*2.*_PI_/l/(l+1);
    //      fprintf(stderr,"%d %e\n",l_read,cl[ref_run][l][sp.index_ct_tt]);
    fscanf(output,"%e",&cl_read);
    cl[1][l][1]=(double)cl_read*2.*_PI_/l/(l+1);
    fscanf(output,"%e",&cl_read);
    cl[1][l][2]=(double)cl_read*2.*_PI_/l/(l+1);
    //fprintf(stderr," %e %d\n",cl_read,sp.index_ct_te);
    fscanf(output,"%e",&cl_read);
    fscanf(output,"%e",&cl_read);
    fscanf(output,"%e",&cl_read);
    fscanf(output,"%e",&cl_read);
    //fscanf(output,"%e",&cl_read);
    if (l_read != l) {
      printf("l_read != l: %d %d\n",l_read,l);
    }
  }
  
  fclose(output);

  fprintf(stderr,"get here\n");

  

  noise_planck(noise,l_max);
   
  fprintf(stderr,"get here\n");

  chi2_planck(cl[0],cl[1],noise,l_max,&chi2);
  
  fprintf(stderr,"get here\n");

  fprintf(stderr,"chi2=%2g \n",chi2);
  
  return 0;
  
}

int chi2_simple(
		double ** cl1,
		double ** cl2,
		int lmax,
		double * chi2) {

  int l;

  *chi2=0.;
  for (l=2; l <= lmax; l++) {
    *chi2 += 
      pow(((cl1[l][0]-cl2[l][0])
	   /cl2[l][0]),2) +
      pow(((cl1[l][2]-cl2[l][2])
	   /cl2[l][2]),2) +
      pow(((cl1[l][1]-cl2[l][1])
	   /cl2[l][1]),2);
  }

  return 0;

}

int chi2_planck(
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

      clTT_th = cl1[l][0]+nl[l][0];
      clEE_th = cl1[l][1]+nl[l][1];
      clTE_th = cl1[l][2];
      
      clTT_obs = cl2[l][0]+nl[l][0];
      clEE_obs = cl2[l][1]+nl[l][1];
      clTE_obs = cl2[l][2];

      //printf("%e %e %e %e %e %e\n",clTT_th,clTT_obs,clEE_th,clEE_obs,clTE_th,clTE_obs);

      det_mixed = 0.5*(clTT_th*clEE_obs+clTT_obs*clEE_th)-clTE_th*clTE_obs;

      det_th = clTT_th*clEE_th-clTE_th*clTE_th;

      det_obs = clTT_obs*clEE_obs-clTE_obs*clTE_obs;

      *chi2 += fsky*(2.*l+1.)*(2.*(det_mixed/det_th-1.)+log(det_th/det_obs));
    }
  }

  return 0;

}

int noise_planck(
		 double ** nl,
		 int lmax
		 ) {
  
  int num_channel,channel,l;
  double * theta, * deltaT, * deltaP; 
  double Tcmb = 2.726;

  num_channel=3;

  theta=malloc(num_channel*sizeof(double));
  deltaT=malloc(num_channel*sizeof(double));
  deltaP=malloc(num_channel*sizeof(double));
  
  /* 100 GHz*/
  theta[0]=9.5/3437.75;          /* converted from arcmin to radian */
  deltaT[0]=(6.8e-6/Tcmb);  /* converted from K to dimensionless */
  deltaP[0]=(10.9e-6/Tcmb); /* converted from K to dimensionless */
  
  /* 143 GHz*/
  theta[1]=7.1/3437.75;
  deltaT[1]=(6.0e-6/Tcmb);
  deltaP[1]=(11.4e-6/Tcmb);
  
  /* 217 GHz*/
  theta[2]=5.0/3437.75;
  deltaT[2]=(13.1e-6/Tcmb);
  deltaP[2]=(26.7e-6/Tcmb);

  for (l=2; l<=lmax; l++) {

    for (channel=0; channel<num_channel; channel++) {

      nl[l][0] +=
	pow((exp(-l*(l+1.)*theta[channel]*theta[channel]/16./log(2.))
	     /theta[channel]/deltaT[channel]),2);
      
      nl[l][1] +=
	pow((exp(-l*(l+1.)*theta[channel]*theta[channel]/16./log(2.))
	     /theta[channel]/deltaP[channel]),2);
    }

    nl[l][0] = 1./nl[l][0];
    nl[l][1] = 1./nl[l][1];
    
  }

  free(theta);
  free(deltaT);
  free(deltaP);

  return 0;
  
}
