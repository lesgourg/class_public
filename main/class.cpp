/** @file class.c
 * Julien Lesgourgues, 17.04.2011
 */

#include "class.h"
#include "output_module.h"
#include "cosmology.h"

int main(int argc, char **argv) {

  Input input;
  precision& pr = input.precision_;      /* for precision parameters */
  background& ba = input.background_;    /* for cosmological background */
  thermo& th = input.thermodynamics_;    /* for thermodynamics */
  perturbs& pt = input.perturbations_;   /* for source functions */
  primordial& pm = input.primordial_;    /* for primordial spectra */
  nonlinear& nl = input.nonlinear_;      /* for non-linear spectra */
  transfers& tr = input.transfers_;      /* for transfer functions */
  spectra& sp = input.spectra_;          /* for output spectra */
  lensing& le = input.lensing_;          /* for lensed spectra */
  output& op = input.output_;            /* for output files */
  ErrorMsg errmsg;                       /* for error messages */

  if (input_init_from_arguments(argc, argv,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg);
    return _FAILURE_;
  }

  if (background_init(&pr,&ba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  Cosmology cosmology = Cosmology(input);
  OutputModule output_module(input, cosmology.GetThermodynamicsModule(), cosmology.GetPerturbationsModule(), cosmology.GetPrimordialModule(), cosmology.GetNonlinearModule(), cosmology.GetSpectraModule(), cosmology.GetLensingModule());

  /****** all calculations done, now free the structures ******/

  if (background_free(&ba) == _FAILURE_) {
    printf("\n\nError in background_free \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  return _SUCCESS_;

}
