#ifndef __EVO__
#define __EVO__

#include "dei_rkck.h"

/**************************************************************/

/**
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int generic_evolver(int (*derivs)(double x, 
				    double * y, 
				    double * dy, 
				    void * parameters_and_workspace,
				    ErrorMsg error_message),
		      double * x,
		      double * y, 
		      int y_size,
		      void * parameters_and_workspace_for_derivs,
		      double tolerance, 
		      double minimum_variation,
		      double timescale,
		      int (*timescale_and_approximation)(double x, 
							 void * parameters_and_workspace,
							 double * timescales,
							 int * approximation_is_changing,
							 ErrorMsg error_message),
		      double timestep_over_timescale,
		      double * x_sampling,
		      int x_size,
		      int (*output)(double x,
				    double y[], 
				    int index_x,
				    void * parameters_and_workspace,
				    ErrorMsg error_message),
		      ErrorMsg error_message);

#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
