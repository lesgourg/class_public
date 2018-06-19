
#ifndef __TRIGONOMETRIC_INTEGRALS__
#define __TRIGONOMETRIC_INTEGRALS__

#include "common.h"


  int cosine_integral(
				 double x,
				 double *Ci,
         ErrorMsg error_message
				 );
	
  int sine_integral(
				 double x,
				 double *Si,
         ErrorMsg error_message
				 );
         
#endif
