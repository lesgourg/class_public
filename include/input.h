/** @file input.h Documented includes for input module */

#include "background.h"
#include "thermodynamics.h"
#include "perturbations.h"
#include "bessel.h"
#include "transfer.h"
#include "primordial.h"
#include "spectra.h"
#include "output.h"

#ifndef __INPUT__
#define __INPUT__

/**************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int input_init(
		 struct background *pba,
		 struct thermo *pth,
		 struct perturbs *ppt,
		 struct bessels *pbs,
		 struct transfers *ptr,
		 struct primordial *ppm,
		 struct spectra *psp,
		 struct output *pop
		 );

#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
