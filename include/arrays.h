/**
 * definitions for module thermodynamics.c
 */

#ifndef __ARRAYS__
#define __ARRAYS__

#include "common.h"

#define _SPLINE_NATURAL_ 0 /**< natural spline: ddy0=ddyn=0 */
#define _SPLINE_EST_DERIV_ 1 /**< spline with estimation of first derivative on both edges */

#define GAUSS_EPSILON 1e-10 /**< Accuracy of determination of gaussian quadrature weights,absiccas */
//#define GAUSS_HERMITE_MAXITER 50 //Currently unused
#define gauss_type_trapezoid 0/**< Trapezoidal integration in [-1,1], unweighted integration, absicca spacing is constant, weights are constant */
#define gauss_type_legendre 1/**< Legendre integration in [-1,1], unweighted integration, absiccas are the roots of the Legendre polynomials, weights are Legendre derivatives */
#define gauss_type_legendre_half 2/**< Legendre integration in [-1,1], unweighted integration, absiccas are only the POSITIVE roots of the Legendre polynomials (half of the roots only), weights are corresponding Legendre derivatives*/
#define gauss_type_chebyshev_1 3/**< Chebyshev integration in [-1,1], weighted integration 1/sqrt(1-x*x), absiccas are the (cosine) roots of the Chebyshev polynomails, weights are constant */
#define gauss_type_chebyshev_2 4/**< Chebyshev integration in [-1,1], weighted integration sqrt(1-x*x), absiccas are the (cosine) roots of the Chebyshev polynomails, weights are sin*sin */
//#define gauss_type_hermite 5/**< Hermite integration in [-inf,inf], weighted integration e^(-x*x), absicccas are the roots of Hermite polynomials, weights are exponentially decreasing */ //Currently unused

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

  int array_spline_table_columns2(
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

int array_integrate_all_trapzd_or_spline(
		   double * array,
		   int n_columns,
		   int n_lines,
           int index_start_spline,
		   int index_x,   /** from 0 to (n_columns-1) */
		   int index_y,
		   int index_ddy,
		   double * result,
		   ErrorMsg errmsg);

  int array_integrate_spline_table_line_to_line(
						double * x_array,
						int n_lines,
						double * array,
						int n_columns,
						int index_y,
						int index_ddy,
						int index_inty,
						ErrorMsg errmsg);

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
			       double * __restrict__ x_array,
			       int n_lines,
			       double * __restrict__ array,
			       double * __restrict__ array_splined,
			       int n_columns,
			       double x,
			       int * __restrict__ last_index,
			       double * __restrict__ result,
			       int result_size, /** from 1 to n_columns */
			       ErrorMsg errmsg);

  int array_interpolate_linear(
			       double * x_array,
			       int n_lines,
			       double * array,
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
					ErrorMsg errmsg);

  int array_interpolate_one_growing_closeby(
                                            double * array,
                                            int n_columns,
                                            int n_lines,
                                            int index_x,   /** from 0 to (n_columns-1) */
                                            double x,
                                            int * last_index,
                                            int index_y,
                                            double * result,
                                            ErrorMsg errmsg);

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

  int array_interpolate_two_bis(
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

  int array_interpolate_cubic_equal(
				    double x0,
				    double dx,
				    double *yarray,
				    int Nx,
				  double x,
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

  int array_trapezoidal_weights(double * __restrict__ x,
                                int n,
                                double * __restrict__ w_trapz,
                                ErrorMsg errmsg);

  int array_trapezoidal_mweights(double * __restrict__ x,
                                int n,
                                double * __restrict__ w_trapz,
                                ErrorMsg errmsg);

  int array_trapezoidal_integral(double * __restrict__ integrand,
                                 int n,
                                 double * __restrict__ w_trapz,
                                 double * __restrict__ I,
                                 ErrorMsg errmsg);

  int array_trapezoidal_convolution(double * __restrict__ integrand1,
                                    double * __restrict__ integrand2,
                                    int n,
                                    double * __restrict__ w_trapz,
                                    double * __restrict__ I,
                                    ErrorMsg errmsg);

  int array_weights_gauss(double* xarray,
                          double* warray,
                          int N,
                          short gauss_type,
                          ErrorMsg err_msg);

  int array_weights_gauss_limits(double* xarray,
                                 double* warray,
                                 double xmin,
                                 double xmax,
                                 int N,
                                 short gauss_type,
                                 ErrorMsg err_msg);

  int array_weights_gauss_rescale_limits(double* xarray,
                                         double* warray,
                                         double* xarrayres,
                                         double* warrayres,
                                         double xmin,
                                         double xmax,
                                         int N,
                                         ErrorMsg err_msg);

#ifdef __cplusplus
}
#endif

#endif
