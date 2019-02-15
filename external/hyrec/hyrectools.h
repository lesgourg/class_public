/*************************** HYRECTOOLS.H ********************************
Multidimensional array creation and freeing functions.
Also, function to make linear arrays and interpolation routines.
Version: January 2011 (identical to November 2010)
**************************************************************************/

double *create_1D_array(unsigned n1);
double **create_2D_array(unsigned n1, unsigned n2);
void free_2D_array(double **matrix, unsigned n1);
double ***create_3D_array(unsigned n1, unsigned n2, unsigned n3);
void free_3D_array(double ***matrix, unsigned n1, unsigned n2);
void maketab(double xmin, double xmax, unsigned Nx, double *xtab);
double rec_interp1d(double x0, double dx, double *ytab, unsigned int Nx, double x);
double rec_interp2d(double x10, double dx1, double x20, double dx2, double **ytab,
		    unsigned int Nx1, unsigned int Nx2, double x1, double x2);
