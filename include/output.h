/** @file output.h Documented includes for output module */

#ifndef __OUTPUT__
#define __OUTPUT__

#include "spectra.h"

struct output {

  FileArg cls_ad;

  FileArg pk;
  
  double z_pk;

  /** @name - flag regulating the amount of information sent to standard output (none if set to zero) */

  //@{
  
  short output_verbose;

  //@}

  ErrorMsg error_message; /**< zone for writing error messages */
};

/*************************************************************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int output_init(
		  struct precision * ppr,
		  struct background * pba,
		  struct perturbs * ppt,
		  struct transfers * ptr,
		  struct spectra * psp,
		  struct output * pop
		  );

#ifdef __cplusplus
}
#endif

#endif
