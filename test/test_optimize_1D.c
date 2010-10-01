/** @file class.c 
 * Julien Lesgourgues, 20.04.2010    
 */
 
#include "class.h"

main(int argc, char **argv) {

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions */
  struct bessels bs;          /* for bessel functions */
  struct transfers tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct spectra sp;          /* for output spectra */
  struct output op;           /* for output files */
 
  ErrorMsg errmsg;

  int i,l_max,l;

  struct file_content fc;

  double parameter,parameter_initial,parameter_logstep;
  int param_num;

  double *** cl;

  parameter_initial=0.02;
  parameter_logstep=1.5;
  param_num=11;
  
  l_max=2500;

  class_alloc(cl,param_num*sizeof(double**),errmsg);
  for (i=0; i<param_num; i++) {
    class_alloc(cl[i],(l_max+1)*sizeof(double*),errmsg);
    for (l=2; l <= l_max; l++)
      class_alloc(cl[i][l],6*sizeof(double),errmsg);
  }

  parser_init(&fc,2,errmsg);

  strcpy(fc.name[0],"output");
  strcpy(fc.value[0],"tCl,pCl");

  strcpy(fc.name[1],"l_max_scalars");
  sprintf(fc.value[1],"%d",l_max);
 
/*   strcpy(fc.name[2],"back_integration_stepsize"); */

  parameter_initial=0.02;
  parameter_logstep=1.5;

  for (i=0; i<param_num; i++) {

    parameter = parameter_initial / exp((double)i*log(parameter_logstep));

/*     sprintf(fc.value[2],"%e",parameter); */
/*     printf("#run with %s\n",fc.value[2]); */

    if (input_init(&fc,&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&op,errmsg) == _FAILURE_) {
      printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
      return _FAILURE_;
    }

    class(&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&op,l_max,cl[i],errmsg);

    for (l=2; l <= 2500; l++)
      printf("%d %e %e %e\n",l,
	     (l*(l+1)/2./_PI_)*cl[i][l][0],  /* here cl is the dimensionless C_l */
	     (l*(l+1)/2./_PI_)*cl[i][l][1],  /* multiply by pow((th.Tcmb*1.e6),2) for muK */
	     (l*(l+1)/2./_PI_)*cl[i][l][2]);
      
    printf("\n");
    
    if (i>20) {
      chi2=0;
      for (l=2; l <= 2500; l++)
	chi2 += ((cl[i][l][0]-cl[0][l][0])/cl[0][l][0])**2+((cl[i][l][1]-cl[0][l][1])/cl[0][l][1])**2;
      printf("%g %g\n",parameter,chi2);
    }
     
  }
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

  for (l=2; l <= l_max; l++) {

    if (spectra_cl_at_l(psp,(double)l,cl[l],junk1,junk2) == _FAILURE_) {
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

}
