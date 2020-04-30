/** @file class.c
 * Julien Lesgourgues, 17.04.2011
 */

#include "class.h"
#include "cosmology.h"
#include "output_module.h"

int main(int argc, char **argv) {

  FileContent fc;
  ErrorMsg error_message;
  if (InputModule::file_content_from_arguments(argc, argv, fc, error_message) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n", error_message);
    return _FAILURE_;
  }

  Cosmology cosmology{fc};

  OutputModule output_module(cosmology.GetInputModule(), cosmology.GetBackgroundModule(), cosmology.GetThermodynamicsModule(), cosmology.GetPerturbationsModule(), cosmology.GetPrimordialModule(), cosmology.GetNonlinearModule(), cosmology.GetSpectraModule(), cosmology.GetLensingModule());

  return _SUCCESS_;

}
