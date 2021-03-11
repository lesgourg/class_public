/******************************* HYRECTOOLS.H ********************************
Multidimensional array creation and freeing functions.
Also, function to make linear arrays and interpolation routines.
*****************************************************************************/

#ifndef __HYRECTOOLS__
#define __HYRECTOOLS__

#define SIZE_ErrorM   2048
double square(double x); 
double cube(double x);
double *create_1D_array(unsigned n1, int *error, char error_message[SIZE_ErrorM]);
double **create_2D_array(unsigned n1, unsigned n2, int *error, char error_message[SIZE_ErrorM]);
void free_2D_array(double **matrix, unsigned n1);
double ***create_3D_array(unsigned n1, unsigned n2, unsigned n3, int *error, char error_message[SIZE_ErrorM]);
void free_3D_array(double ***matrix, unsigned n1, unsigned n2);
void maketab(double xmin, double xmax, unsigned Nx, double *xtab);
double rec_interp1d(double x0, double dx, double *ytab, unsigned int Nx, double x, int *error, char error_message[SIZE_ErrorM]);
double rec_interp2d(double x10, double dx1, double x20, double dx2, double **ytab,
                    unsigned int Nx1, unsigned int Nx2, double x1, double x2, int *error, char error_message[SIZE_ErrorM]);
double rec_interpol_G(double x, double *xtab, double *ytab, unsigned int Nx, int *error, char error_message[SIZE_ErrorM]);

#endif
