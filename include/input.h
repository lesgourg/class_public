/** @file input.h Documented includes for input module */

#ifndef __INPUT__
#define __INPUT__

#include "common.h"
#include "parser.h"
#include "background.h"
#include "thermodynamics.h"
#include "perturbations.h"
#include "bessel.h"
#include "transfer.h"
#include "primordial.h"
#include "spectra.h"
#include "output.h"

/* macro for opening file and returning error if it failed */
#define class_read_double(name,destination)				\
  do {									\
    flag1=parser_read_double(pfc,name,&param1,errmsg);			\
    if (flag1 == _SUCCESS_)						\
      destination = param1;						\
  } while(0);


#define class_read_int(name,destination)				\
  do {									\
    flag1=parser_read_int(pfc,name,&int1,errmsg);			\
    if (flag1 == _SUCCESS_)						\
      destination = int1;						\
  } while(0);

#define class_read_string(name,destination)				\
  do {									\
    flag1=parser_read_string(pfc,name,&string1,errmsg);			\
    if (flag1 == _SUCCESS_)						\
      strcpy(destination,string1);					\
  } while(0);


/**************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int input_init_from_arguments(
		 int argc, 
		 char **argv,
		 struct precision * ppr,
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

  int input_init(
		 struct file_content * pfc,
		 struct precision * ppr,
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

  int input_default_params(
			   struct background *pba,
			   struct thermo *pth,
			   struct perturbs *ppt,
			   struct bessels * pbs,
			   struct transfers *ptr,
			   struct primordial *ppm,
			   struct spectra *psp,
			   struct output *pop
			   );
  
  int input_default_precision(
			      struct precision * ppp
			      );

  int get_machine_precision(double * smallest_allowed_variation);

#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
