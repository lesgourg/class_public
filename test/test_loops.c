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
  struct spectra_nl nl;       /* for calculation of non-linear spectra */
 
  ErrorMsg errmsg;

  int i,l_max;

  struct file_content fc;

  l_max=2400;
  
  parser_init(&fc,9,errmsg);

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

  for (i=0; i<=10; i++) {

    sprintf(fc.value[4],"%e",0.001+i*0.004);
    printf("#run with omega_b = %s\n",fc.value[4]);

    class(&fc,&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&op,&nl,l_max,errmsg);

  }

}

int class(
	  struct file_content * pfc,
	  struct precision * ppr,
	  struct background * pba,
	  struct thermo * pth,
	  struct perturbs * ppt,
	  struct bessels * pbs,
	  struct transfers * ptr,
	  struct primordial * ppm,
	  struct spectra * psp,
	  struct output * pop,
	  struct spectra_nl * pnl,
	  int l_max,
	  ErrorMsg errmsg) {

  int l;
  double cl[6]; /* note: the actual size of this vector should be sp.ct_size, but 
		   this value is only computed in spectra_init(). Anyway 6 is the 
		   maximum value that sp.ct_size can take, including TT, EE, TE, BB, pp and Tp. */
  
  double ** junk1;
  double ** junk2;

  if (input_init(pfc,ppr,pba,pth,ppt,pbs,ptr,ppm,psp,pop,pnl,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }

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
      
    if (spectra_cl_at_l(psp,(double)l,cl,junk1,junk2) == _FAILURE_) {
      printf("\n\nError in spectra_cl_at_l \n=>%s\n",psp->error_message);
      return _FAILURE_;
    }
    
    printf("%d %e %e %e\n",l,
	   (l*(l+1)/2./_PI_)*cl[psp->index_ct_tt],  /* here cl is the dimensionless C_l */
	   (l*(l+1)/2./_PI_)*cl[psp->index_ct_ee],  /* multiply by pow((th.Tcmb*1.e6),2) for muK */
	   (l*(l+1)/2./_PI_)*cl[psp->index_ct_te]);
      
  }

  printf("\n");

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
