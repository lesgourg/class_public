/** @file input.h Documented includes for input module */

#ifndef __INPUT__
#define __INPUT__

#include "class.h"

/**************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int input_init(
		 int argc, 
		 char **argv,
		 struct background *pba,
		 struct thermo *pth,
		 struct perturbs *ppt,
		 struct bessels *pbs,
		 struct transfers *ptr,
		 struct primordial *ppm,
		 struct spectra *psp,
		 struct output *pop,
		 ErrorMsg errmsg
		 );

  int input_init_params(
			struct file_content * pfc,
			struct background *pba,
			struct thermo *pth,
			struct perturbs *ppt,
			struct bessels * pbs,
			struct transfers *ptr,
			struct primordial *ppm,
			struct spectra *psp,
			struct output *pop,
			ErrorMsg errmsg
			);

  int input_init_default(
			   struct background *pba,
			   struct thermo *pth,
			   struct perturbs *ppt,
			   struct bessels * pbs,
			   struct transfers *ptr,
			   struct primordial *ppm,
			   struct spectra *psp,
			   struct output *pop);

int input_check_arguments_of_main(
				  int argc, 
				  char **argv, 
				  char * input,
				  char * precision,
				  ErrorMsg errmsg);

#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
