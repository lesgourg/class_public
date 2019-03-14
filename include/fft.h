/**
 * Function definitions of the FFT tool, used by transfer.c and spectra.c
 * For more information, see fft.c
 *
 * @author Nils Sch.
 */

#ifndef FFT_DEFINED
#define FFT_DEFINED

struct FFT_plan{

  int N;
  double* cos_vals;
  double* sin_vals;

  double* temp_real;
  double* temp_imag;

};

struct DCT_plan{

  int N;
  struct FFT_plan* fft_plan;

  double* input_real;
  double* input_imag;
  double* output_real;
  double* output_imag;
};


void FFT_noplan(double* input_real, double* input_imag, double* output_real, double* output_imag, int N);
void FFT_planned(double* input_real, double* input_imag, double* output_real, double* output_imag, struct FFT_plan* plan);

void FFT_real_short_noplan(double* input_real_1, double* input_real_2, double* output_real_1,double* output_imag_1, double* output_real_2,double* output_imag_2, int N);
void FFT_real_noplan(double* input_real_1, double* input_real_2, double* output_real_1,double* output_imag_1, double* output_real_2,double* output_imag_2, int N);

void FFT_real_short_planned(double* input_real_1, double* input_real_2, double* output_real_1,double* output_imag_1, double* output_real_2,double* output_imag_2, struct FFT_plan* plan);
void FFT_real_planned(double* input_real_1, double* input_real_2, double* output_real_1,double* output_imag_1, double* output_real_2,double* output_imag_2, struct FFT_plan* plan);

void FFT_planner_init(int N, struct FFT_plan** planned);
void FFT_planner_free(struct FFT_plan** planned);


void DCT_noplan(double* input, double* output, int N);
void DCT_planned(double* input, double* output, struct DCT_plan* dct_plan);

void DCT_planner_init(int N, struct DCT_plan** planned);
void DCT_planner_free(struct DCT_plan** planned);
#endif
