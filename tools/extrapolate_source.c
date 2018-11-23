/***
 * Module with functions providing some extrapolations for the source functions
 * These functionbs are used by nonlinear.c
 * Written by Samuel Brieden, 2018
 */

#include "extrapolate_source.h"

/**
 * Fill the (previously allocated) array of extrapolated source values
 */
int extrapolate_source(
                       double* k_extrapolated,
                       int k_size_original,
                       int k_size_extrapolated,
                       double* source_array,
                       short extrapolation_method,
                       double* source_extrapolated,
                       double k_eq,
                       double h,
                       ErrorMsg errMsg
                       ) {
  double k_max;
  double source_max;
  int index_extr;
  int index_k_0;
  double k_0;

  /**
   * Copy existing wavenumbers k and sources
   * */
  class_protect_memcpy(source_extrapolated,source_array,k_size_original*sizeof(double));

  /**
   * Get last source and k, which are used in (almost) all methods
   * */
  k_max = k_extrapolated[k_size_original-1];
  source_max = source_array[k_size_original-1];

  index_k_0 = k_size_original-2;
  k_0 = k_extrapolated[index_k_0];

  for(index_extr=k_size_original;index_extr<k_size_extrapolated;++index_extr){

    switch(extrapolation_method){
      /**
       * Extrapolate by assuming the source to vanish
       * Has terrible discontinuity
       * */
    case extrapolation_zero:
      {
        source_extrapolated[index_extr]=0.0;
        break;
      }
      /**
       * Extrapolate starting from the maximum value, assuming  growth ~ ln(k)
       * Has a terrible bend in log slope, discontinuity only in derivative
       * */
    case extrapolation_only_max:
      {
        source_extrapolated[index_extr]=source_max*(log(k_extrapolated[index_extr])/log(k_max));
        break;
      }
      /**
       * Extrapolate starting from the maximum value, assuming  growth ~ ln(k)
       * Here we use k in h/Mpc instead of 1/Mpc as it is done in the CAMB implementation of HMcode
       * Has a terrible bend in log slope, discontinuity only in derivative
       * */
    case extrapolation_only_max_units:
      {
        source_extrapolated[index_extr]=source_max*(log(1/h*k_extrapolated[index_extr])/log(1/h*k_max));
        break;
      }
      /**
       * Extrapolate assuming source ~ ln(a*k) where a is obtained from the data at k_0
       * Mostly continuous derivative, quite good
       * */
    case extrapolation_max_scaled:
      {
        double scaled_factor = exp((source_array[index_k_0]*log(k_max)-source_array[k_size_original-1]*log(k_0))/(source_array[k_size_original-1]-source_array[index_k_0]));
        source_extrapolated[index_extr]=source_max*(log(scaled_factor*k_extrapolated[index_extr])/log(scaled_factor*k_max));
        break;
      }
      /**
       * Extrapolate assuming source ~ ln(e+a*k) where a is estimated like is done in original HMCode
       * */
    case extrapolation_hmcode:
      {
        //Extrapolation formula taken from original HMCode
        double scaled_factor = 1.8/(13.41*k_eq);
        source_extrapolated[index_extr]=source_max*(log(_E_+scaled_factor*k_extrapolated[index_extr])/log(_E_+scaled_factor*k_max));
        break;
      }
      /**
       * If the user has a complicated model and wants to interpolate differently,
       * they can define their interpolation here and switch to using it instead
       * */
    case extrapolation_user_defined:
      {
        class_stop(errMsg,"Method of source extrapolation 'user_defined' was not yet defined.");
        break;
      }
    }
    //end switch
  }
  //end while
  return _SUCCESS_;
}

/**
 * Define an extended array of log-spaced k values up to k_max_extrapolated
 */
int extrapolate_k(
                  double* k_array,
                  int k_size_original,
                  double* k_extrapolated,
                  double k_per_decade,
                  double k_max_extrapolated,
                  ErrorMsg errMsg
                  ){
  double k_max;
  int index_extr;
  int index_running_k;
  double running_k=0;

  /**
   * Copy existing wavenumbers k
   * */
  class_protect_memcpy(k_extrapolated,k_array,k_size_original*sizeof(double));

  /**
   * Get last k
   * */
  k_max = k_array[k_size_original-1];
  index_running_k=0;

  while(running_k<k_max_extrapolated && index_running_k < _MAX_NUM_EXTRAPOLATION_){
    index_extr = index_running_k+k_size_original;
    index_running_k++;
    running_k = k_max * pow(10,index_running_k/k_per_decade);
    k_extrapolated[index_extr] = running_k;
  }
  //end while
  class_test(index_running_k==_MAX_NUM_EXTRAPOLATION_,
             errMsg,
             "could not reach k_max_extrapolated = %.10e starting from k_max = %.10e with k_per_decade of %.10e in %i maximal steps (_MAX_NUM_INTERPOLATION_ in extrapolate_soruce.h)",
             k_max_extrapolated,k_max,k_per_decade,_MAX_NUM_EXTRAPOLATION_
             );
  return _SUCCESS_;
}

/**
 * Find the size of an extended array of log-spaced k values up to k_max_extrapolated
 */
int get_extrapolated_source_size(
                                 double k_per_decade,
                                 double k_max,
                                 double k_max_extrapolated,
                                 int k_size_original,
                                 int* size_extrapolated_source,
                                 ErrorMsg errMsg
                                 ){
  double running_k=0;
  int index_running_k;
  /**
   * The k sampling should be the same lines of code as in the extrapolate_source function
   * */
  index_running_k=0;
  while(running_k<k_max_extrapolated && index_running_k < _MAX_NUM_EXTRAPOLATION_){
    index_running_k++;
    running_k = k_max * pow(10,index_running_k/k_per_decade);
  }
  class_test(index_running_k==_MAX_NUM_EXTRAPOLATION_,
             errMsg,
             "could not reach k_max_extrapolated = %.10e starting from k_max = %.10e with k_per_decade of %.10e in %i maximal steps (_MAX_NUM_INTERPOLATION_ in extrapolate_soruce.h)",
             k_max_extrapolated,k_max,k_per_decade,_MAX_NUM_EXTRAPOLATION_
             );
  *size_extrapolated_source = k_size_original+index_running_k;
  return _SUCCESS_;
}
