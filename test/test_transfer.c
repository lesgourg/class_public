/** @file class.c 
 * Julien Lesgourgues, 18.04.2010    
 */
 
#include "precision.h"
#include "background.h"
#include "thermodynamics.h"
#include "perturbations.h"
#include "bessel.h"
#include "transfer.h"
#include "primordial.h"
#include "spectra.h"

main() {

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions */
  struct bessels bs;          /* for bessel functions */
  struct transfers tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct spectra sp;          /* for output spectra */
  struct spectra op;          /* for output files */
 
  if (precision_init(&pr) == _FAILURE_) {
    printf("\n\nError running precision_init \n=>%s\n",pr.error_message); 
    return _FAILURE_;
  }

  if (input_init(&ba,&th,&pt,&bs,&tr,&pm,&sp,&op) == _FAILURE_) {
    printf("\n\nError running input_init"); 
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

  if (bessel_init(&pr,&ba,&pt,&bs) == _FAILURE_) {
    printf("\n\nError in bessel_init \n =>%s\n",bs.error_message);
    return _FAILURE_;
  }

  if (transfer_init(&pr,&ba,&th,&pt,&bs,&tr) == _FAILURE_) {
    printf("\n\nError in transfer_init \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  /****** output the transfer functions ******/

  printf("Output of transfer functions\n");

  int index_mode=pt.index_md_scalars;
  int index_ic  =pt.index_ic_ad;
  int index_type=pt.index_tp_t;
  int index_l=tr.l_size[index_mode]-5;
/*   int index_l = 20; */

  /* here you can output the transfer functions 
     at some k's of your choice */
 
/*   double k; */
/*   double transfer; */

/*   for (k=1.e-4; k<1.e0; k*=1.1) {  */

/*     if (transfer_functions_at_k( */
/* 				index_mode, */
/* 				index_ic, */
/* 				index_type, */
/* 				index_l, */
/* 				k, */
/* 				&transfer */
/* 				) == _FAILURE_) { */
/*       printf("\n\nError in transfer_function_at_k \n=>%s\n",tr.error_message); */
/*       return _FAILURE_;; */
/*     } */

/*     printf("%e %e\n",k,transfer);  */

/*   }  */

  /* here you can output the full tabulated arrays*/

  int index_k;
  double transfer;

  for (index_k=0; index_k<tr.k_size[index_mode]; index_k++) { 

    transfer=tr.transfer[index_mode]
      [((index_ic * pt.tp_size + index_type)
	* tr.l_size[index_mode] + index_l)
       * tr.k_size[index_mode] + index_k];
    
    printf("%d %e %e\n",tr.l[index_mode][index_l],tr.k[index_mode][index_k],transfer); 

  } 

  /************************************************************/

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

}
