#include "evolver_rkck.h"

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
		    ErrorMsg error_message) {

  int next_index_x;
  double x1,x2=0.,timestep,timescale;
  struct generic_integrator_workspace gi;
  double * dy;
  short call_output;

  class_test(x_ini > x_sampling[x_size-1],
	     error_message,
	     "called with x=%e, last x_sampling=%e",x_ini,x_sampling[x_size-1]);

  next_index_x=0;

  while (x_sampling[next_index_x] < x_ini) next_index_x++;

  class_call(initialize_generic_integrator(y_size, &gi),
	     gi.error_message,
	     error_message);

  class_alloc(dy,y_size*sizeof(double),error_message);

  x1=x_ini;

  call_output = _FALSE_;

  while ((x1 < x_end) && (next_index_x<x_size)) {

    class_call((*evaluate_timescale)(x1,
				     parameters_and_workspace_for_derivs,
				     &timescale,
				     error_message),
	       error_message,
	       error_message);

    timestep = timestep_over_timescale * timescale;

    class_test(fabs(timestep/x1) < minimum_variation,
	       error_message,
	       "integration step =%e < machine precision : leads either to numerical error or infinite loop",fabs(timestep/x1));

    if (x1 + 2.* timestep < x_sampling[next_index_x]) {
      x2 = x1 + timestep;
    }
    else {
      x2 = x_sampling[next_index_x];
      call_output = _TRUE_;
    }

    if (x2 > x_end) {
      x2 = x_end;
      call_output = _FALSE_;
    }

    if (print_variables != NULL) {

      if (x1 == x_ini) {

	class_call((*derivs)(x1,
			     y,
			     dy,
			     parameters_and_workspace_for_derivs,
			     error_message),
		   error_message,
		   error_message);
      }

      class_call((*print_variables)(x1,
				    y,
				    dy,
				    parameters_and_workspace_for_derivs,
				    error_message),
		 error_message,
		 error_message);
    }

    class_call(generic_integrator(derivs,
				  x1,
				  x2,
				  y,
				  parameters_and_workspace_for_derivs,
				  tolerance,
				  x1*minimum_variation,
				  &gi),
	       gi.error_message,
	       error_message);

    if (call_output == _TRUE_) {

      class_call((*derivs)(x2,
			   y,
			   dy,
			   parameters_and_workspace_for_derivs,
			   error_message),
		 error_message,
		 error_message);

      class_call((*output)(x2,
			   y,
			   dy,
			   next_index_x,
			   parameters_and_workspace_for_derivs,
			   error_message),
		 error_message,
		 error_message);

      call_output = _FALSE_;

      next_index_x++;

    }

    x1 = x2;

  }

  /* a last call is compulsory to ensure that all quantitites in
     y,dy,parameters_and_workspace_for_derivs are updated to the last
     point in the covered range */
  class_call((*derivs)(x1,
		       y,
		       dy,
		       parameters_and_workspace_for_derivs,
		       error_message),
	     error_message,
	     error_message);

  if (print_variables != NULL)
    class_call((*print_variables)(x1,
				  y,
				  dy,
				  parameters_and_workspace_for_derivs,
				  error_message),
	       error_message,
	       error_message);

  class_call(cleanup_generic_integrator(&gi),
	     gi.error_message,
	     error_message);

  free(dy);

  return _SUCCESS_;

}
