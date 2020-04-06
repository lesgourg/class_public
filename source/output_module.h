#ifndef OUTPUT_MODULE_H
#define OUTPUT_MODULE_H

#include "base_module.h"
#include "input.h"
#include "lensing_module.h"
#include "spectra_module.h"

class OutputModule : public BaseModule {
public:
  OutputModule(const Input& input, const SpectraModule& spectra_module, const LensingModule& lensing_module);
private:
  int output_total_cl_at_l(int l, double* cl);
  int output_init();
  int output_cl();
  int output_pk(enum pk_outputs pk_output);
  int output_tk();
  int output_background();
  int output_thermodynamics();
  int output_perturbations();
  int output_primordial();
  int output_print_data(FILE *out, const char titles[_MAXTITLESTRINGLENGTH_], double *dataptr, int tau_size);
  int output_open_cl_file(FILE ** clfile, FileName filename, char * first_line, int lmax);
  int output_one_line_of_cl(FILE * clfile, double l, double * cl, int ct_size);
  int output_open_pk_file(FILE ** pkfile, FileName filename, char * first_line, double z);
  int output_one_line_of_pk(FILE * tkfile, double one_k, double one_pk);

  const SpectraModule& spectra_module_;
  const LensingModule& lensing_module_;
};

#endif //OUTPUT_MODULE_H
