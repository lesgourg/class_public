/**
 * definitions for module thermodynamics.c
 */

#ifndef __TOOLS_ARRAYS__
#define __TOOLS_ARRAYS__

#include "common.h"

#define _SPLINE_NATURAL_ 0 /**< natural spline: ddy0=ddyn=0 */
#define _SPLINE_EST_DERIV_ 1 /**< spline with estimation of first derivative on both edges */

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int array_derive(
		   double * array,
		   int n_columns,
		   int n_lines,
		   int index_x,   /** from 0 to (n_columns-1) */
		   int index_y,
		   int index_dydx,
		   ErrorMsg errmsg);

  int array_derive_spline(
			  double * x_array,
			  int n_lines,
			  double * array,
			  double * array_splined,
			  int n_columns,
			  int index_y,
			  int index_dydx,
			  ErrorMsg errmsg);

  int array_derive_spline_table_line_to_line(
				       double * x_array,
				       int n_lines,
				       double * array,
				       int n_columns,
				       int index_y,
				       int index_ddy,
				       int index_dy,
				       ErrorMsg errmsg);

  int array_derive1_order2_table_line_to_line(
				       double * x_array,
				       int n_lines,
				       double * array,
				       int n_columns,
				       int index_y,
				       int index_dy,
				       ErrorMsg errmsg);

  int array_derive2_order2_table_line_to_line(
				       double * x_array,
				       int n_lines,
				       double * array,
				       int n_columns,
				       int index_y,
				       int index_dy,
				       int index_ddy,
				       ErrorMsg errmsg);

  int array_derive_two(
		   double * array,
		   int n_columns,
		   int n_lines,
		   int index_x,   /** from 0 to (n_columns-1) */
		   int index_y,
		   int index_dydx,
		   int index_ddydxdx,
		   ErrorMsg errmsg);



  int array_spline(
		   double * array,
		   int n_columns,
		   int n_lines,
		   int index_x,   /** from 0 to (n_columns-1) */
		   int index_y,
		   int index_ddydx2,
		   short spline_mode,
		   ErrorMsg errmsg);

  int array_spline_table_line_to_line(
				      double * x, /* vector of size x_size */
				      int x_size,
				      double * array,
				      int n_columns,
				      int index_y,
				      int index_ddydx2,
				      short spline_mode,
				      ErrorMsg errmsg);

  int array_spline_table_columns(
		       double * x,
		       int x_size,
		       double * y_array,
		       int y_size,
		       double * ddy_array,
		       short spline_mode,
		       ErrorMsg errmsg);

  int array_spline_table_lines(
		       double * x,
		       int x_size,
		       double * y_array,
		       int y_size,
		       double * ddy_array,
		       short spline_mode,
		       ErrorMsg errmsg
		       );

  int array_logspline_table_lines(
				  double * x,
				  int x_size,
				  double * y_array,
				  int y_size,
				  double * ddlny_array,
				  short spline_mode,
				  ErrorMsg errmsg
				  );

  int array_spline_table_one_column(
				    double * x, /* vector of size x_size */
				    int x_size,
				    double * y_array, /* array of size x_size*y_size with elements
							 y_array[index_y*x_size+index_x] */
				    int y_size,
				    int index_y,
				    double * ddy_array, /* array of size x_size*y_size */
				    short spline_mode,
				    ErrorMsg errmsg
				    );

  int array_logspline_table_one_column(
				    double * x, /* vector of size x_size */
				    int x_size,
				    int x_stop,
				    double * y_array, /* array of size x_size*y_size with elements
							 y_array[index_y*x_size+index_x] */
				    int y_size,
				    int index_y,
				    double * ddlogy_array, /* array of size x_size*y_size */
				    short spline_mode,
				    ErrorMsg errmsg
				    );

  int array_integrate_all_spline(
				 double * array,
				 int n_columns,
				 int n_lines,
				 int index_x,
				 int index_y,
				 int index_ddy,
				 double * result,
				 ErrorMsg errmsg
				 );

  int array_integrate(
		      double * array,
		      int n_columns,
		      int n_lines,
		      int index_x,   /** from 0 to (n_columns-1) */
		      int index_y,
		      int index_int_y_dx,
		      ErrorMsg errmsg);

  int array_integrate_all(
		      double * array,
		      int n_columns,
		      int n_lines,
		      int index_x,   /** from 0 to (n_columns-1) */
		      int index_y,
		      double * result);

  int array_integrate_ratio(
			    double * array,
			    int n_columns,
			    int n_lines,
			    int index_x,   /** from 0 to (n_columns-1) */
			    int index_y1,
			    int index_y2,
			    int index_int_y1_over_y2_dx,
			    ErrorMsg errmsg);

  int array_interpolate(
			double * array,
			int n_columns,
			int n_lines,
			int index_x,   /** from 0 to (n_columns-1) */
			double x,
			int * last_index,
			double * result,
			int result_size,
			ErrorMsg errmsg); /** from 1 to n_columns */

  int array_interpolate_spline(
			       double * x_array,
			       int n_lines,
			       double * array,
			       double * array_splined,
			       int n_columns,
			       double x,
			       int * last_index,
			       double * result,
			       int result_size, /** from 1 to n_columns */
			       ErrorMsg errmsg);

  int array_interpolate_logspline(
				  double * x_array,
				  int n_lines,
				  double * array,
				  double * array_logsplined,
				  int n_columns,
				  double x,
				  int * last_index,
				  double * result,
				  int result_size, /** from 1 to n_columns */
				  ErrorMsg errmsg);

  int array_interpolate_growing_closeby(
					double * array,
					int n_columns,
					int n_lines,
					int index_x,   /** from 0 to (n_columns-1) */
					double x,
					int * last_index,
					double * result,
					int result_size,
					ErrorMsg errmsg); /** from 1 to n_columns */

  int array_interpolate_spline_growing_closeby(
					       double * x_array,
					       int n_lines,
					       double * array,
					       double * array_splined,
					       int n_columns,
					       double x,
					       int * last_index,
					       double * result,
					       int result_size, /** from 1 to n_columns */
					       ErrorMsg errmsg);

  int array_interpolate_spline_growing_hunt(
					       double * x_array,
					       int n_lines,
					       double * array,
					       double * array_splined,
					       int n_columns,
					       double x,
					       int * last_index,
					       double * result,
					       int result_size, /** from 1 to n_columns */
					       ErrorMsg errmsg);

  int array_interpolate_spline_one_column(
					  double * x_array,
					  int x_size,
					  double * y_array, /* array of size x_size*y_size with elements
							       y_array[index_y*x_size+index_x] */
					  int y_size,
					  int index_y,
					  double * ddy_array, /* array of size x_size*y_size */
					  double x,   /* input */
					  double * y, /* output */
					  ErrorMsg errmsg
					  );

  int array_interpolate_extrapolate_spline_one_column(
					  double * x_array,
					  int x_size,
					  double * y_array, /* array of size x_size*y_size with elements
							       y_array[index_y*x_size+index_x] */
					  int y_size,
					  int index_y,
					  double * ddy_array, /* array of size x_size*y_size */
					  double x,   /* input */
					  double * y, /* output */
					  ErrorMsg errmsg
					  );

  int array_interpolate_extrapolate_logspline_loglinear_one_column(
								   double * x_array,
								   int x_size,
								   int x_stop,
								   double * y_array, /* array of size x_size*y_size with elements
											y_array[index_y*x_size+index_x] */
								   int y_size,
								   int index_y,
								   double * ddlogy_array, /* array of size x_size*y_size */
								   double x,   /* input */
								   double * y, /* output */
								   ErrorMsg errmsg
								   );

  int array_interpolate_extrapolate_logspline_one_column(
					  double * x_array,
					  int x_size,
					  int x_stop,
					  double * y_array, /* array of size x_size*y_size with elements
							       y_array[index_y*x_size+index_x] */
					  int y_size,
					  int index_y,
					  double * ddlogy_array, /* array of size x_size*y_size */
					  double x,   /* input */
					  double * y, /* output */
					  ErrorMsg errmsg
					  );

  /** interpolate to get y(x), when x_i and y_i are in two different arrays*/
  int array_interpolate_two(
			    double * array_x,
			    int n_columns_x,
			    int index_x,   /** from 0 to (n_columns_x-1) */
			    double * array_y,
			    int n_columns_y,
			    int n_lines,  /** must be the same for array_x and array_y */
			    double x,
			    double * result,
			    int result_size, /** from 1 to n_columns_y */
			    ErrorMsg errmsg);

  int array_interpolate_two_arrays_one_column(
					      double * array_x, /* assumed to be a vector (i.e. one column array) */
					      double * array_y,
					      int n_columns_y,
					      int index_y, /* between 0 and (n_columns_y-1) */
					      int n_lines,  /** must be the same for array_x and array_y */
					      double x,
					      double * result,
					      ErrorMsg errmsg);

  int array_interpolate_equal(
			    double * array,
			    int n_colums,
			    int n_lines,
			    double x,
			    double x_min,
			    double x_max,
			    double * result,
			    ErrorMsg errmsg);

  int array_interpolate_parabola(double x1,
				 double x2,
				 double x3,
				 double x,
				 double y1,
				 double y2,
				 double y3,
				 double * y,
				 double * dy,
				 double * ddy,
				 ErrorMsg errmsg);

  int array_smooth(double * array,
		   int n_columns,
		   int n_lines,
		   int index, /** from 0 to (n_columns-1) */
		   int radius,
		   ErrorMsg errmsg);


#ifdef __cplusplus
}
#endif

#endif
