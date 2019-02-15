/*************************** HYRECTOOLS.C ********************************
Multidimensional array creation and freeing functions.
Also, function to make linear arrays and interpolation routines.

Version: January 2011 (unchanged from November 2010 version)
**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
      
#include "hyrectools.h"

/************************************************************************** 
Creates a [n1] array. 
***************************************************************************/

double *create_1D_array(unsigned n1){

   double *matrix = (double *) calloc(n1, sizeof(double));
   if (matrix == NULL) {
      fprintf(stderr, "memory issue in create_1D_array\n");
      exit(1);
    }
   return matrix;
}

/************************************************************************** 
Creates a [n1][n2] array. 
***************************************************************************/

double **create_2D_array(unsigned n1, unsigned n2){
  
 unsigned i;
   double **matrix = (double **) calloc(n1, sizeof(double *));
   if (matrix == NULL){
      fprintf(stderr, "memory issue in create_2D_array\n");
      exit(1);
    }
   for (i = 0; i < n1; i++)  matrix[i] = create_1D_array(n2);
   
   return matrix;
}

/******************************************************************************* 
Frees the memory of a [n1][] array. 
********************************************************************************/

void free_2D_array(double **matrix, unsigned n1){

    unsigned i;
    for (i = 0; i < n1; i++) free(matrix[i]);
   
    free(matrix);
}

/********************************************************************************* 
Creates a [n1][n2][n3] matrix. 
**********************************************************************************/

double ***create_3D_array(unsigned n1, unsigned n2, unsigned n3){

    unsigned i;
    double ***matrix = (double ***) calloc(n1, sizeof(double **));
    if (matrix == NULL) {
      fprintf(stderr, "memory issue in create_3D_array\n");
      exit(1);
    }
    for (i = 0; i < n1; i++)  matrix[i] = create_2D_array(n2, n3);
         
    return matrix;
}

/*********************************************************************************** 
Frees memory of a [n1][n2][] matrix
***********************************************************************************/

void free_3D_array(double ***matrix, unsigned n1, unsigned n2) {
   unsigned i;
   for (i = 0; i < n1; i++)  free_2D_array(matrix[i], n2);
   
   free(matrix);
}

/********************************************************************************************
Making an evenly spaced array.
Input: xmin, xmax, Nx, xtab.
xtab is changed, s.t. xtab[0] = xmin, xtab[Nx-1] = xmax, evenly spaced.
ATTENTION: there will be Nx points, and Nx-1 intervals.
**********************************************************************************************/

void maketab(double xmin, double xmax, unsigned Nx, double *xtab){
    unsigned i;
    double h = (xmax - xmin)/(Nx - 1.0);

    for (i = 0; i < Nx; i++) xtab[i] = xmin + i * h;
}

/************************************************************************************
 Interpolation routine for 1-D table.  Uses cubic interpolation assuming
 uniformly spaced x-values, x0 ... x0+(Nx-1)*dx.
*************************************************************************************/

double rec_interp1d(double x0, double dx, double *ytab, unsigned int Nx, double x) {

  long ix;
  double frac;

  /* Check if in range */
  if (dx > 0 && (x<x0 || x>x0+dx*(Nx-1))) {
    fprintf(stderr, "Error: rec_interp1d: x-value out of range in interpolation routine.\n");
    exit(1);
  }
  if (dx < 0 && (x>x0 || x<x0+dx*(Nx-1))) {
    fprintf(stderr, "Error: rec_interp1d: x-value out of range in interpolation routine.\n");
    exit(1);
  }

  /* Identify location to interpolate */
  ix = (long)floor((x-x0)/dx);
  if (ix<1) ix=1;
  if (ix>Nx-3) ix=Nx-3;
  frac = (x-x0)/dx-ix;
  ytab += ix-1;

  /* Return value */
  return(
    -ytab[0]*frac*(1.-frac)*(2.-frac)/6.
    +ytab[1]*(1.+frac)*(1.-frac)*(2.-frac)/2.
    +ytab[2]*(1.+frac)*frac*(2.-frac)/2.
    +ytab[3]*(1.+frac)*frac*(frac-1.)/6.
  );
}

/************************************************************************************
 Interpolation routine for 2-D table.  Uses bicubic interpolation assuming
 uniformly spaced x1 and x2-values.
*************************************************************************************/

double rec_interp2d(double x10, double dx1, double x20, double dx2, double **ytab,
  unsigned int Nx1, unsigned int Nx2, double x1, double x2) {

  int j;
  long ix1;
  double frac1;
  double temp[4];

  /* Check if in range in x1 */
  if (x1<x10 || x1>x10+dx1*(Nx1-1)) {
    fprintf(stderr, "Error: rec_interp2d: x1-value out of range in interpolation routine.\n");
    exit(1);
  }

  /* Identify location to interpolate in x1 */
  ix1 = (long)floor((x1-x10)/dx1);
  if (ix1<1) ix1=1;
  if (ix1>Nx1-3) ix1=Nx1-3;
  frac1 = (x1-x10)/dx1-ix1;
  ytab += ix1-1;  

  /* Get values interpolated over x2 at the 4 neighboring points in x1 */
  for(j=0;j<4;j++) temp[j] = rec_interp1d(x20,dx2,ytab[j],Nx2,x2);

  return(
    -temp[0]*frac1*(1.-frac1)*(2.-frac1)/6.
    +temp[1]*(1.+frac1)*(1.-frac1)*(2.-frac1)/2.
    +temp[2]*(1.+frac1)*frac1*(2.-frac1)/2.
    +temp[3]*(1.+frac1)*frac1*(frac1-1.)/6.
  );
}

/*********************************************************************************************/
