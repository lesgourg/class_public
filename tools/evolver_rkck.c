#include "evolver_rkck.h"

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
				  double dy[],
				  int index_x,
				  void * parameters_and_workspace,
				  ErrorMsg error_message),
		    ErrorMsg error_message) {

  int next_index_x;
  double x1,x2,timestep;
  int approximation_is_changing;
  struct generic_integrator_workspace gi;
  double * dy;

  x1=*x;

  class_test(x1 > x_sampling[x_size-1],
	     error_message,
	     "called with x=%e, last x_sampling=%e",x1,x_sampling[x_size-1]);
  
  next_index_x=0;
  while (x_sampling[next_index_x] < x1) next_index_x++; 

  class_call(initialize_generic_integrator(y_size, &gi),
	     gi.error_message,
	     error_message);
  
  class_alloc(dy,y_size*sizeof(double),error_message);

  timestep = timestep_over_timescale * timescale;
      
  class_test(timestep < minimum_variation,
	     error_message,
	     "integration step =%e < machine precision : leads either to numerical error or infinite loop",timestep);
  
  while (next_index_x < x_size) {

    while (x1 < x_sampling[next_index_x]) {
      
      if (x1 + 2.* timestep < x_sampling[next_index_x]) 
	x2 = x1 + timestep;
      else 
	x2 = x_sampling[next_index_x];
      
      class_call(generic_integrator(derivs,
				    x1,
				    x2,
				    y,
				    parameters_and_workspace_for_derivs,
				    tolerance,
				    minimum_variation,
				    &gi),
		 gi.error_message,
		 error_message);

      x1 = x2;

      class_call((*timescale_and_approximation)(x1, 
						parameters_and_workspace_for_derivs,
						&timescale,
						&approximation_is_changing,
						error_message),
		 error_message,
		 error_message);

      if (approximation_is_changing > 0) {
	
	class_test(approximation_is_changing > 1,
		   error_message,
		   "cannot change more than one approximation at a time");
	
	class_call(cleanup_generic_integrator(&gi),
		   gi.error_message,
		   error_message);
		  
	free(dy);
    
	*x = x1;

	return _SUCCESS_;
      
      }

      timestep = timestep_over_timescale * timescale;
      
    }
    
    class_call((*derivs)(x1,
			 y,
			 dy,
			 parameters_and_workspace_for_derivs,
			 error_message),
	       error_message,
	       error_message);
    
    class_call((*output)(x1,
			 y,
			 dy,
			 next_index_x,
			 parameters_and_workspace_for_derivs,
			 error_message),
	       error_message,
	       error_message);
    
    next_index_x++;
  
  }

  class_call(cleanup_generic_integrator(&gi),
	     gi.error_message,
	     error_message);
  
  free(dy);

  *x = x1;

  return _SUCCESS_;

}
