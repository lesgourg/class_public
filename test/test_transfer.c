/** @file class.c 
 * Julien Lesgourgues, 17.04.2011    
 */
 
/* this main only runs the modules up to the transfer one */

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
  struct nonlinear nl;        /* for non-linear spectra */
  struct lensing le;          /* for lensed spectra */
  struct output op;           /* for output files */
  ErrorMsg errmsg;            /* for error messages */

  if (input_init_from_arguments(argc, argv,&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&nl,&le,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }

  if (background_init(&pr,&ba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_init(&pr,&ba,&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (perturb_init(&pr,&ba,&th,&pt) == _FAILURE_) {
    printf("\n\nError in perturb_init \n=>%s\n",pt.error_message);
    return _FAILURE_;
  }

  if (bessel_init(&pr,&bs) == _FAILURE_) {
    printf("\n\nError in bessel_init \n =>%s\n",bs.error_message);
    return _FAILURE_;
  }

  if (transfer_init(&pr,&ba,&th,&pt,&bs,&tr) == _FAILURE_) {
    printf("\n\nError in transfer_init \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  /****** output the transfer functions ******/

  printf("Output of transfer functions (l, k, Delta)\n");

  /* 1) select the mode, initial condition, type and multipole of the
     function you want to plot: */

  int index_mode=pt.index_md_scalars;
  int index_ic  =pt.index_ic_ad;
  int index_type=tr.index_tt_lcmb;
  //int index_type=tr.index_tt_density+2;

  /* 2) here is an illustration of how to output the transfer
     functions at some (k,l)'s of your choice */
 
  /*
  int index_l = 0;
  double k=3.6e-4;
  double transfer;

  if (transfer_functions_at_k(&tr,
			      index_mode,
			      index_ic,
			      index_type,
			      index_l,
			      k,
			      &transfer
			      ) == _FAILURE_) {
    printf("\n\nError in transfer_function_at_k \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  printf("%d %e %e\n",tr.l[index_l],k,transfer);
   
  */
 
  /* 3) here you can output the full tabulated arrays for all k and l's*/

  int index_k;
  int index_l;
  double transfer;

  for (index_l=0; index_l<tr.l_size[index_mode]; index_l++) { 
    //for (index_l=20; index_l<21; index_l++) { 
    for (index_k=0; index_k<tr.k_size[index_mode]; index_k++) { 
      
      transfer=tr.transfer[index_mode]
	[((index_ic * tr.tt_size[index_mode] + index_type)
	  * tr.l_size[index_mode] + index_l)
	 * tr.k_size[index_mode] + index_k];
      
      if (transfer != 0.) {
	printf("%d %e %e\n",tr.l[index_l],tr.k[index_mode][index_k],transfer); 
      }
    }
    
    printf("\n\n");
    
  } 

  /****** all calculations done, now free the structures ******/

  if (transfer_free(&tr) == _FAILURE_) {
    printf("\n\nError in transfer_free \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  if (bessel_free(&bs) == _FAILURE_) {
    printf("\n\nError in bessel_free \n=>%s\n",bs.error_message);
    return _FAILURE_;
  }

  if (perturb_free(&pt) == _FAILURE_) {
    printf("\n\nError in perturb_free \n=>%s\n",pt.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_free(&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_free \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (background_free(&ba) == _FAILURE_) {
    printf("\n\nError in background_free \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  return _SUCCESS_;

}
