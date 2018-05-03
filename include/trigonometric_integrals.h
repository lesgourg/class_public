
#ifndef __TRIGONOMETRIC_INTEGRALS__
#define __TRIGONOMETRIC_INTEGRALS__

#include "common.h"


  int nonlinear_hmcode_ci(
				 double x,
				 double *Ci,
         ErrorMsg error_message
				 );
	
  int nonlinear_hmcode_si(
				 double x,
				 double *Si,
         ErrorMsg error_message
				 );
         
#endif
