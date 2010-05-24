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

  /****** here you could output the background evolution ******/

  if (background_free(&ba) == _FAILURE_) {
    printf("\n\nError in background_free \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

}
