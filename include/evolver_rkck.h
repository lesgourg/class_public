#ifndef __EVO_RK__
#define __EVO_RK__

#include "dei_rkck.h"

/**************************************************************/

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int evolver_rk(int (*derivs)(double x,
				    double * y,
				    double * dy,
				    void * parameters_and_workspace,
				    ErrorMsg error_message),
		      double x_ini,
		      double x_end,
		      double * y,
		      int * used_in_output,
		      int y_size,
		      void * parameters_and_workspace_for_derivs,
		      double tolerance,
		      double minimum_variation,
		      int (*evaluate_timescale)(double x,
						void * parameters_and_workspace,
						double * timescale,
						ErrorMsg error_message),
		      double timestep_over_timescale,
		      double * x_sampling,
		      int x_size,
		      int (*output)(double x,
				    double y[],
				    double dy[],
				    int index_x,
				    void * parameters_and_workspace,
				    ErrorMsg error_message),
		      int (*print_variables)(double x,
					     double y[],
					     double dy[],
					     void * parameters_and_workspace,
					     ErrorMsg error_message),
		      ErrorMsg error_message);

#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
