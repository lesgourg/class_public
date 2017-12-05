/***
 * extrapolate_transfer provides extrapolations for the source functions
 * This file is used by both matter.c and nonlinear.c
 */
#include "extrapolate_source.h"
#define _E_ 2.718281828459045235360287471352662497757247093699959574966967627724076630353547594571382178525166427427466391932003059921817413596629043572900334295260595630738132328627943490763233829880753195251019011573834187930702154089149934884167509244761460668082264800168477411853742345442437107539077744992069551702761838606261331384583000752044933826560297606737113200709328709127443747047230696977209310141692836819025515108657463772111252389784425056953696
int extrapolate_source(
                struct background* pba,
                double* k_array,
                int k_size,
                double* source_array,
                short extrapolation_method,
                double* k_extrapolated,
                double* source_extrapolated,
                double k_per_decade,
                double k_max_extrapolated,
                ErrorMsg errMsg
                ){
  double k_max;
  double source_max;
  int index_extr;
  int index_running_k;
  double running_k;
  /**
   * Kopy existing wavenumbers k and sources
   * */
  class_protect_memcpy(k_extrapolated,k_array,k_size*sizeof(double));
  class_protect_memcpy(source_extrapolated,source_array,k_size*sizeof(double));
  /**
   * Get last source and k, which are used in (almost) all methods
   * */
  k_max = k_array[k_size-1];
  source_max = source_array[k_size-1];
  index_running_k=0;
  
  
  int index_k_0 = k_size-2;
  double k_0 = k_array[index_k_0];
  double k_eq = pba->a_eq*pba->H_eq;
        
  while(running_k<k_max_extrapolated && index_running_k < _MAX_NUM_EXTRAPOLATION_){
    index_extr = index_running_k+k_size;
    index_running_k++;
    running_k = k_max * pow(10,index_running_k/k_per_decade);
    k_extrapolated[index_extr] = running_k;
    switch(extrapolation_method){
      /**
       * Extrapolate by assuming the source to vanish
       * */
      case extrapolation_zero:
        source_extrapolated[index_extr]=0.0;
        break;
      /**
       * Extrapolate starting from the maximum value, assuming  growth ~ ln(k)
       * */
      case extrapolation_only_max:
        source_extrapolated[index_extr]=source_max*(log(k_extrapolated[index_extr])/log(k_max));
        break;
      /**
       * Extrapolate assuming source ~ ln(a*k) where a is obtained from the data at k_0
       * */
      case extrapolation_max_a:
        double scaled_factor = (source_array[index_k_0]*log(k_max)-source_array[k_size-1]*log(k_0))/(source_array[k_size-1]-source_array[index_k_0]);
        source_extrapolated[index_extr]=source_max*(log(scaled_factor*k_extrapolated[index_extr])/log(scaled_factor*k_max));
        break;
      /**
       * Extrapolate assuming source ~ ln(a*k) where a is estimated like is done in HMCode
       * */
      case hm_code:
        //Extrapolation formula taken from HM_Code
        double scaled_factor = 1.8/(13.41*k_eq);
        source_extrapolated[index_extr]=source_max*(log(_E_+scaled_factor*k_extrapolated[index_extr])/log(_E_+scaled_factor*k_max));
        break;
      /**
       * If the user has a complicated model and wants to interpolate differently, 
       * they can define their interpolation here and switch to using it instead
       * */
      case user_defined:
        class_stop(errMsg,"Method of source extrapolation 'user_defined' was not yet defined.");
        break;
    }
  }
  class_test(index_running_k==_MAX_NUM_EXTRAPOLATION_,
             errMsg
             "could not reach k_max_extrapolated = %.10e starting from k_max = %.10e with k_per_decade of %.10e in %i maximal steps (_MAX_NUM_INTERPOLATION_ in extrapolate_soruce.h)",
             k_max_extrapolated,k_max,k_per_decade,_MAX_NUM_EXTRAPOLATION_,
             );
  return _SUCCESS_;
}
int get_extrapolated_source_size(
          double k_per_decade,
          double k_max,
          double k_max_extrapolated,
          int* size_extrapolated_source,
          ErrorMsg errMsg
          ){
  double running_k;
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
             errMsg
             "could not reach k_max_extrapolated = %.10e starting from k_max = %.10e with k_per_decade of %.10e in %i maximal steps (_MAX_NUM_INTERPOLATION_ in extrapolate_soruce.h)",
             k_max_extrapolated,k_max,k_per_decade,_MAX_NUM_EXTRAPOLATION_,
             );
  *size_extrapolated_source = index_running_k+1;
  return _SUCCESS_;
}
