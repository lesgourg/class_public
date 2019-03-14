/** @file fft.c Documented fast fourier transform module
 *
 * Nils Schoenberg, 12.10.2017
 *
 * This module computes the fast fourier transform (FFT) of any function
 */
 #include "fft.h"
#include <math.h>
//Mathematical constants used in this file
 #define MATH_PI_2 6.2831853071795864769252867665590057683943387987502
 #define INV_SQRT_2 1/sqrt(2)
//Function implementations
 /**
 * Computes the fast fourier transform of some arbitrary input array input_real and input_double of sizes N
 * Returns the output by writing into output_real and output_imag
 * 
 * It is assumed that all the arrays are allocated and of size N
 * 
 * For recursion there is a stepsize parameter
 * If the full FFT should be calculated, set 
 * 
 *          ** stepsize = 1 ** 
 * 
 * */
void FFT(double* input_real, double* input_imag, double* output_real, double* output_imag, int N, int stepsize){
  //Larger base cases up to N == 8 have proven advantageous
  if (N == 8){
    //FFT N==4 EVEN
    //FFT(input_real, input_imag, output_real, output_imag, 4, 2 * stepsize);
    // i = 0
    double temp_even_real = input_real[0] + input_real[4 * stepsize];
    double temp_even_imag = input_imag[0] + input_imag[4 * stepsize];
    double temp_odd_real = input_real[2*stepsize] + input_real[6 * stepsize];
    double temp_odd_imag = input_imag[2*stepsize] + input_imag[6 * stepsize];
    output_real[0] = temp_even_real + temp_odd_real;
    output_imag[0] = temp_even_imag + temp_odd_imag;
    output_real[2] = temp_even_real - temp_odd_real;
    output_imag[2] = temp_even_imag - temp_odd_imag;
    // i = 1
    temp_even_real = input_real[0] - input_real[4 * stepsize];
    temp_even_imag = input_imag[0] - input_imag[4 * stepsize];
    temp_odd_real = input_real[2*stepsize] - input_real[6 * stepsize];
    temp_odd_imag = input_imag[2*stepsize] - input_imag[6 * stepsize];
    output_real[1] = temp_even_real + temp_odd_imag;
    output_imag[1] = temp_even_imag - temp_odd_real;
    output_real[3] = temp_even_real - temp_odd_imag;
    output_imag[3] = temp_even_imag + temp_odd_real;
    //FFT N==4 ODD
    //FFT(input_real + stepsize, input_imag + stepsize, output_real + 4, output_imag + 4, 4, 2 * stepsize);
    // i = 0
    temp_even_real = input_real[stepsize] + input_real[5 * stepsize];
    temp_even_imag = input_imag[stepsize] + input_imag[5 * stepsize];
    temp_odd_real = input_real[3*stepsize] + input_real[7 * stepsize];
    temp_odd_imag = input_imag[3*stepsize] + input_imag[7 * stepsize];
    output_real[4] = temp_even_real + temp_odd_real;
    output_imag[4] = temp_even_imag + temp_odd_imag;
    output_real[6] = temp_even_real - temp_odd_real;
    output_imag[6] = temp_even_imag - temp_odd_imag;
    // i = 1
    temp_even_real = input_real[stepsize] - input_real[5 * stepsize];
    temp_even_imag = input_imag[stepsize] - input_imag[5 * stepsize];
    temp_odd_real = input_real[3*stepsize] - input_real[7 * stepsize];
    temp_odd_imag = input_imag[3*stepsize] - input_imag[7 * stepsize];
    output_real[5] = temp_even_real + temp_odd_imag;
    output_imag[5] = temp_even_imag - temp_odd_real;
    output_real[7] = temp_even_real - temp_odd_imag;
    output_imag[7] = temp_even_imag + temp_odd_real;
    //FINAL FFT N==8
    //i==0
    temp_even_real = output_real[0];
    temp_even_imag = output_imag[0];
    temp_odd_real = output_real[4];
    temp_odd_imag = output_imag[4];
    output_real[0] = temp_even_real + temp_odd_real;
    output_imag[0] = temp_even_imag + temp_odd_imag;
    output_real[4] = temp_even_real - temp_odd_real;
    output_imag[4] = temp_even_imag - temp_odd_imag;
    //i==1
    temp_even_real = output_real[1];
    temp_even_imag = output_imag[1];
    temp_odd_real = output_real[5];
    temp_odd_imag = output_imag[5];
    output_real[1] = temp_even_real + INV_SQRT_2*temp_odd_real + INV_SQRT_2*temp_odd_imag;
    output_imag[1] = temp_even_imag + INV_SQRT_2*temp_odd_imag - INV_SQRT_2*temp_odd_real;
    output_real[5] = temp_even_real - INV_SQRT_2*temp_odd_real - INV_SQRT_2*temp_odd_imag;
    output_imag[5] = temp_even_imag + INV_SQRT_2*temp_odd_real - INV_SQRT_2*temp_odd_imag;
    //i==2
    temp_even_real = output_real[2];
    temp_even_imag = output_imag[2];
    temp_odd_real = output_real[6];
    temp_odd_imag = output_imag[6];
    output_real[2] = temp_even_real + temp_odd_imag;
    output_imag[2] = temp_even_imag - temp_odd_real;
    output_real[6] = temp_even_real - temp_odd_imag;
    output_imag[6] = temp_even_imag + temp_odd_real;
    //i==3
    temp_even_real = output_real[3];
    temp_even_imag = output_imag[3];
    temp_odd_real = output_real[7];
    temp_odd_imag = output_imag[7];
    output_real[3] = temp_even_real - INV_SQRT_2*temp_odd_real + INV_SQRT_2*temp_odd_imag;
    output_imag[3] = temp_even_imag - INV_SQRT_2*temp_odd_imag - INV_SQRT_2*temp_odd_real;
    output_real[7] = temp_even_real + INV_SQRT_2*temp_odd_real - INV_SQRT_2*temp_odd_imag;
    output_imag[7] = temp_even_imag + INV_SQRT_2*temp_odd_real + INV_SQRT_2*temp_odd_imag;
    return;
  }
  else if (N < 8){
    if (N == 4){
      // i = 0
      double temp_even_real = input_real[0] + input_real[2 * stepsize];
      double temp_even_imag = input_imag[0] + input_imag[2 * stepsize];
      double temp_odd_real = input_real[stepsize] + input_real[3 * stepsize];
      double temp_odd_imag = input_imag[stepsize] + input_imag[3 * stepsize];
      output_real[0] = temp_even_real + temp_odd_real;
      output_imag[0] = temp_even_imag + temp_odd_imag;
      output_real[2] = temp_even_real - temp_odd_real;
      output_imag[2] = temp_even_imag - temp_odd_imag;
      // i = 1
      temp_even_real = input_real[0] - input_real[2 * stepsize];
      temp_even_imag = input_imag[0] - input_imag[2 * stepsize];
      temp_odd_real = input_real[stepsize] - input_real[3 * stepsize];
      temp_odd_imag = input_imag[stepsize] - input_imag[3 * stepsize];
      output_real[1] = temp_even_real + temp_odd_imag;
      output_imag[1] = temp_even_imag - temp_odd_real;
      output_real[3] = temp_even_real - temp_odd_imag;
      output_imag[3] = temp_even_imag + temp_odd_real;
      return;
    }
    else
    if (N == 2){
      output_real[0] = input_real[0] + input_real[stepsize];
      output_real[1] = input_real[0] - input_real[stepsize];
      output_imag[0] = input_imag[0] + input_imag[stepsize];
      output_imag[1] = input_imag[0] - input_imag[stepsize];
      return;
    }
    else
    if (N == 1){
      output_real[0] = input_real[0];
      output_imag[0] = input_imag[0];
      return;
    }
  }
  else{
    //Use the butterfly algorithm to compute the fft of even and odd sets
    FFT(input_real, input_imag, output_real, output_imag, N / 2, 2 * stepsize);
    FFT(input_real + stepsize, input_imag + stepsize, output_real + N / 2, output_imag + N / 2, N / 2, 2 * stepsize);

        int i,j;
    //Reunite even and odd sets 
    for (i = 0 , j = N/2; i < N / 2; ++i,++j){
      double temp_even_real = output_real[i];
      double temp_even_imag = output_imag[i];
      double temp_odd_real = output_real[j];
      double temp_odd_imag = output_imag[j];
      //These twiddle factors cos_val and sin_val could be read from file instead
      //It will depend on many things whether or not that is advantagous
      // (e.g. on the number of frequqencies N)
      double cos_val = cos(i*MATH_PI_2 / ((double)(N)));
      double sin_val = sin(i*MATH_PI_2 / ((double)(N)));
      output_real[i] = temp_even_real + cos_val*temp_odd_real + sin_val*temp_odd_imag;
      output_imag[i] = temp_even_imag + cos_val*temp_odd_imag - sin_val*temp_odd_real;
      output_real[j] = temp_even_real - cos_val*temp_odd_real - sin_val*temp_odd_imag;
      output_imag[j] = temp_even_imag + sin_val*temp_odd_real - cos_val*temp_odd_imag;
    }
    return;
  }
}
 /**
 * Computes the fast fourier transform of some purely real inputs input_real_1 and input_real_2 of size N
 * Returns the output by writing into output_real_i and output_imag_i for input_real_i with i=1,2
 * 
 * It is assumed that all the arrays are allocated and of size N
 * 
 * Returns full output, arrays of size N
 * */
void FFT_real(double* input_real_1, double* input_real_2, double* output_real_1,double* output_imag_1, double* output_real_2,double* output_imag_2, int N){
  FFT(input_real_1, input_real_2, output_real_1, output_real_2, N,1);
  //output_real_1[0] remains the same
  double temp1,temp2;
  int i;
  output_imag_1[0] = 0.0;
  output_imag_2[0] = 0.0;
  output_imag_1[N/2] = 0.0;
  output_imag_2[N/2] = 0.0;
  for (i = 1; i < N/2; ++i){
    temp1 = output_real_1[i];
    temp2 = output_real_2[i];
    output_real_1[i] = 0.5*(temp1 + output_real_1[N - i]);
    output_real_2[i] = 0.5*(temp2 + output_real_2[N - i]);
    output_imag_1[i] = 0.5*(temp2 - output_real_2[N - i]);
    output_imag_2[i] = 0.5*(output_real_1[N - i] - temp1);
  }
  for (i = 0; i < N / 2; ++i){
    output_real_1[i + N / 2] = output_real_1[N / 2 - i];
    output_real_2[i + N / 2] = output_real_2[N / 2 - i];
    output_imag_1[i + N / 2] = - output_imag_1[N / 2 - i];
    output_imag_2[i + N / 2] = - output_imag_2[N / 2 - i];
  }
}
 /**
 * Computes the fast fourier transform of some purely real inputs input_real_1 and input_real_2 of size N
 * Returns the output by writing into output_real_i and output_imag_i for input_real_i with i=1,2
 * 
 * It is assumed that all the arrays are allocated and of size N
 * 
 * Only returns N/2 output arrays (still have to be allocated at size N)
 * 
 * For any real fourier transformation c_(-n) = c_n and thus c_(N-n) = c_n for finite fourier transformation of size N
 * */
void FFT_real_short(double* input_real_1, double* input_real_2, double* output_real_1, double* output_imag_1, double* output_real_2, double* output_imag_2, int N){
  //Only computes first N/2 elements, since others are related by symmetry
  FFT(input_real_1, input_real_2, output_real_1, output_real_2, N, 1);
  double temp1,temp2;
  int i;
  output_imag_1[0] = 0.0;
  output_imag_2[0] = 0.0;
  output_imag_1[N/2] = 0.0;
  output_imag_2[N/2] = 0.0;
  for (i = 1; i < N/2; ++i){
  temp1 = output_real_1[i];
  temp2 = output_real_2[i];
    output_real_1[i] = 0.5*(temp1 + output_real_1[N - i]);
    output_real_2[i] = 0.5*(temp2 + output_real_2[N - i]);
    output_imag_1[i] = 0.5*(temp2 - output_real_2[N - i]);
    output_imag_2[i] = 0.5*(output_real_1[N - i] - temp1);
  }
}

