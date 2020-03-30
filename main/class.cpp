/** @file class.c
 * Julien Lesgourgues, 17.04.2011
 */

#include "class.h"

int main(int argc, char **argv) {

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions */
  struct primordial pm;       /* for primordial spectra */
  struct nonlinear nl;        /* for non-linear spectra */
  struct transfers tr;        /* for transfer functions */
  struct spectra sp;          /* for output spectra */
  struct lensing le;          /* for lensed spectra */
  struct output op;           /* for output files */
  ErrorMsg errmsg;            /* for error messages */

  if (input_init_from_arguments(argc, argv,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,errmsg) == _FAILURE_) {
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

  if (primordial_init(&pr,&pt,&pm) == _FAILURE_) {
    printf("\n\nError in primordial_init \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }

  if (nonlinear_init(&pr,&ba,&th,&pt,&pm,&nl) == _FAILURE_) {
    printf("\n\nError in nonlinear_init \n=>%s\n",nl.error_message);
    return _FAILURE_;
  }

  if (transfer_init(&pr,&ba,&th,&pt,&nl,&tr) == _FAILURE_) {
    printf("\n\nError in transfer_init \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  if (spectra_init(&pr,&ba,&pt,&pm,&nl,&tr,&sp) == _FAILURE_) {
    printf("\n\nError in spectra_init \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }

  if (lensing_init(&pr,&pt,&sp,&nl,&le) == _FAILURE_) {
    printf("\n\nError in lensing_init \n=>%s\n",le.error_message);
    return _FAILURE_;
  }

  if (output_init(&ba,&th,&pt,&pm,&tr,&sp,&nl,&le,&op) == _FAILURE_) {
    printf("\n\nError in output_init \n=>%s\n",op.error_message);
    return _FAILURE_;
  }

  /****** all calculations done, now free the structures ******/

  if (lensing_free(&le) == _FAILURE_) {
    printf("\n\nError in lensing_free \n=>%s\n",le.error_message);
    return _FAILURE_;
  }

  if (spectra_free(&sp) == _FAILURE_) {
    printf("\n\nError in spectra_free \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }

  if (transfer_free(&tr) == _FAILURE_) {
    printf("\n\nError in transfer_free \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  if (nonlinear_free(&nl) == _FAILURE_) {
    printf("\n\nError in nonlinear_free \n=>%s\n",nl.error_message);
    return _FAILURE_;
  }

  if (primordial_free(&pm) == _FAILURE_) {
    printf("\n\nError in primordial_free \n=>%s\n",pm.error_message);
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
