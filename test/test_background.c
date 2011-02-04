/** @file class.c 
 * Julien Lesgourgues, 18.04.2010    
 */
 
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
  ErrorMsg errmsg;

  if (input_init_from_arguments(argc, argv,&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&nl,&le,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }

  if (background_init(&pr,&ba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  /****** here you could output the background evolution ******/

  int index_eta;

  for (index_eta=0; index_eta<ba.bt_size; index_eta++) {

/*     fprintf(stdout, */
/* 	    "%e %e %e %e\n", */
/* 	    ba.eta_table[index_eta], */
/* 	    1./(1+ba.z_table[index_eta]), */
/* 	    ba.background_table[index_eta*ba.bg_size+ba.index_bg_a], */
/* 	    ba.background_table[index_eta*ba.bg_size+ba.index_bg_H]); */

    fprintf(stdout,
	    "%e %e %e %e %e\n",
	    ba.eta_table[index_eta],
	    1./(1+ba.z_table[index_eta]),
	    ba.background_table[index_eta*ba.bg_size+ba.index_bg_rho_nur],
	    ba.background_table[index_eta*ba.bg_size+ba.index_bg_rho_ncdm1],
	    ba.background_table[index_eta*ba.bg_size+ba.index_bg_p_ncdm1]);

  }

  if (background_free(&ba) == _FAILURE_) {
    printf("\n\nError in background_free \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  return _SUCCESS_;

}
