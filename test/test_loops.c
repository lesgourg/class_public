/** @file class.c 
 * Julien Lesgourgues, 17.04.2011    
 */
 
/* this main loops over several models. It illustrates how the code
   could be interfaced with a parameter extraction code. */

#include "class.h"

int main(int argc, char **argv) {

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions */
  struct bessels bs;          /* for bessel functions */
  struct transfers tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct spectra sp;          /* for output spectra */
  struct output op;           /* for output files */
  struct spectra_nl nl;       /* for calculation of non-linear spectra */
  ErrorMsg errmsg;            /* for error messages */

  int i; 
  int l_max;
  struct file_content fc;

  /* choose a value of l_max in C_l's */
  l_max=2400;
  
  /* all parameters for which we don't want to keep default values
     should be passed to the code through a file_content
     structure. Create such a structure with the size you need: 9 in
     this exemple */
  parser_init(&fc,9,errmsg);

  /* assign values to these 9 parameters. Some will be fixed, some
     will be varied in the loop. */
  strcpy(fc.name[0],"output");
  strcpy(fc.value[0],"tCl,pCl");
  
  strcpy(fc.name[1],"l_max_scalars");
  sprintf(fc.value[1],"%d",l_max);

  strcpy(fc.name[2],"T_cmb");
  sprintf(fc.value[2],"%e",2.726);

  strcpy(fc.name[3],"H0");
  sprintf(fc.value[3],"%e",72.);

  strcpy(fc.name[4],"omega_b");
  sprintf(fc.value[4],"%e",0.024);

  strcpy(fc.name[5],"omega_cdm");
  sprintf(fc.value[5],"%e",0.05);
  
  strcpy(fc.name[6],"z_reio");
  sprintf(fc.value[6],"%e",10.);

  strcpy(fc.name[7],"A_s_ad");
  sprintf(fc.value[7],"%e",2.3e-9);

  strcpy(fc.name[8],"n_s_ad");
  sprintf(fc.value[8],"%e",1.);

  /* now loop over values of some of these parameters: in this
     exemple, omega_b */
  for (i=0; i<=10; i++) {

    /* assign one value to omega_b */
    sprintf(fc.value[4],"%e",0.001+i*0.004);
    printf("#run with omega_b = %s\n",fc.value[4]);

    /* calls class (except the bessel module which has been called
       once and for all) and return the C_l's*/
    class_assuming_bessels_computed(&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&nl,&le,&op,l_max,cl[i],errmsg);
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
	  struct nonlinear * pnl,
	  struct lensing *ple,
	  struct output * pop,
	  int l_max,
	  double ** cl,
	  ErrorMsg errmsg) {

  int l;

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

  if (nonlinear_init(ppr,pba,pth,ppm,psp,pnl) == _FAILURE_) {
    printf("\n\nError in nonlinear_init \n=>%s\n",pnl->error_message);
    return _FAILURE_;
  }

  if (lensing_init(ppr,ppt,psp,ple) == _FAILURE_) {
    printf("\n\nError in lensing_init \n=>%s\n",ple->error_message);
    return _FAILURE_;
  }

  for (l=2; l <= l_max; l++) {

    if (output_total_cl_at_l(psp,ple,pop,(double)l,cl[l]) == _FAILURE_) {
      printf("\n\nError in spectra_cl_at_l \n=>%s\n",psp->error_message);
      return _FAILURE_;
    }

  }

  /****** all calculations done, now free the structures ******/

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

  if (nonlinear_init(ppr,pba,pth,ppm,psp,pnl) == _FAILURE_) {
    printf("\n\nError in nonlinear_init \n=>%s\n",pnl->error_message);
    return _FAILURE_;
  }

  if (lensing_init(ppr,ppt,psp,ple) == _FAILURE_) {
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

  /****** all calculations done, now free the structures ******/

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
