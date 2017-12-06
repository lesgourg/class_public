#ifndef __EXTRAPOLATE_SOURCE__
#define __EXTRAPOLATE_SOURCE__

#include "common.h"
#define extrapolation_zero 0
#define extrapolation_only_max 1
#define extrapolation_max_scaled 2
#define extrapolation_hmcode 3
#define extrapolation_user_defined 4
#define _MAX_NUM_EXTRAPOLATION_ 10000
int extrapolate_source(
        double* k_extrapolated,
        int k_size_original,
        int k_size_extrapolated,
        double* source_array,
        short extrapolation_method,
        double* source_extrapolated,
        double k_eq,
        ErrorMsg errMsg
        );
int extrapolate_k(
        double* k_array,
        int k_size_original,
        double* k_extrapolated,
        double k_per_decade,
        double k_max_extrapolated,
        ErrorMsg errMsg
        );
int get_extrapolated_source_size(
        double k_per_decade,
        double k_max,
        double k_max_extrapolated,
        int k_size_original,
        int* size_extrapolated_source,
        ErrorMsg errMsg
        );
#endif
