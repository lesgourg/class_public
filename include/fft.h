/**
 * Function definitions of the FFT tool, used by transfer.c and spectra.c
 * For more information, see fft.c
 */

#ifndef FFT_DEFINED
#define FFT_DEFINED
/**
 * Compute FFT of two real inputs and store into two complex outputs of sizes N
 * Assumes arrays are allocated and of size N
 */
void FFT_real(double* input_real_1, double* input_real_2, double* output_real_1, double* output_imag_1, double* output_real_2, double* output_imag_2, int N);
/**
 * Compute FFT of two real inputs and store into two complex outputs of sizes N
 * Assumes arrays are allocated and of size N
 * Only gives up nonredundant part for real FFT which have c_(N-n)=c_n
 */
void FFT_real_short(double* input_real_1, double* input_real_2, double* output_real_1, double* output_imag_1, double* output_real_2, double* output_imag_2, int N);
/**
 * Compute FFT of single complex input and stores into single complex output of size N
 * Assumes arrays are allocated and of size N
 */
void FFT(double* input_real, double* input_imag, double* output_real, double* output_imag, int N, int stepsize);
#endif
