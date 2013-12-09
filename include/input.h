/** @file input.h Documented includes for input module */

#ifndef __INPUT__
#define __INPUT__

#include "common.h"
#include "parser.h"
#include "quadrature.h"
#include "background.h"
#include "thermodynamics.h"
#include "perturbations.h"
#include "transfer.h"
#include "primordial.h"
#include "spectra.h"
#include "nonlinear.h"
#include "lensing.h"
#include "output.h"

/* macro for reading parameter values with routines from the parser */
#define class_read_double(name,destination)				\
  do {									\
    class_call(parser_read_double(pfc,name,&param1,&flag1,errmsg),      \
	       errmsg,							\
	       errmsg);							\
    if (flag1 == _TRUE_)						\
      destination = param1;						\
  } while(0);


#define class_read_int(name,destination)				\
  do {									\
    class_call(parser_read_int(pfc,name,&int1,&flag1,errmsg),		\
 	       errmsg,							\
	       errmsg);							\
    if (flag1 == _TRUE_)						\
      destination = int1;						\
  } while(0);

#define class_read_string(name,destination)				\
  do {									\
    class_call(parser_read_string(pfc,name,&string1,&flag1,errmsg),	\
 	       errmsg,							\
	       errmsg);							\
    if (flag1 == _TRUE_)						\
      strcpy(destination,string1);					\
  } while(0);

#define class_read_double_one_of_two(name1,name2,destination)		\
  do {									\
    class_call(parser_read_double(pfc,name1,&param1,&flag1,errmsg),	\
	       errmsg,							\
	       errmsg);							\
    class_call(parser_read_double(pfc,name2,&param2,&flag2,errmsg),	\
	       errmsg,							\
	       errmsg);							\
    class_test((flag1 == _TRUE_) && (flag2 == _TRUE_),			\
	       errmsg,							\
	       "In input file, you can only enter one of %s, %s, choose one", \
	       name1,name2);						\
    if (flag1 == _TRUE_)						\
      destination = param1;						\
    if (flag2 == _TRUE_)						\
      destination = param2;						\
  } while(0);

#define class_at_least_two_of_three(a,b,c)		\
  ((a == _TRUE_) && (b == _TRUE_)) ||		\
  ((a == _TRUE_) && (c == _TRUE_)) ||		\
  ((b == _TRUE_) && (c == _TRUE_))

#define class_none_of_three(a,b,c)				\
  (a == _FALSE_) && (b == _FALSE_) && (c == _FALSE_)

/* macro for reading parameter values with routines from the parser */
#define class_read_list_of_doubles_or_default(name,destination,default,siz)	\
  do {									\
    class_call(parser_read_list_of_doubles(pfc,name,			\
	&entries_read,&(destination),&flag1,errmsg),			\
	       errmsg,							\
	       errmsg);							\
    if (flag1 == _TRUE_){						\
        class_test(entries_read != siz,errmsg,			\
             "Number of entries in %s, %d, does not match expected number, %d.", \
		name,entries_read,siz);				\
    }else{								\
	class_alloc(destination,siz*sizeof(double),errmsg);		\
	for(n=0; n<siz; n++) destination[n] = default;		\
    }									\
  } while(0);

#define class_read_list_of_integers_or_default(name,destination,default,siz) \
  do {									\
    class_call(parser_read_list_of_integers(pfc,name,			\
	&entries_read,&(destination),&flag1,errmsg),			\
	       errmsg,							\
	       errmsg);							\
    if (flag1 == _TRUE_){						\
        class_test(entries_read != siz,errmsg,			\
             "Number of entries in %s, %d, does not match expected number, %d.", \
		name,entries_read,siz);				\
    }else{								\
	class_alloc(destination,siz*sizeof(int),errmsg);		\
	for(n=0; n<siz; n++) destination[n] = default;		\
    }									\
  } while(0);

#define class_read_list_of_doubles(name,destination,siz)			\
  do {									\
    class_call(parser_read_list_of_doubles(pfc,name,			\
	&entries_read,&(destination),&flag1,errmsg),			\
	       errmsg,							\
	       errmsg);							\
    class_test(flag1 == _FALSE_,errmsg,					\
	"Entry %s is required but not found!",name)			\
        class_test(entries_read != siz,errmsg,			\
             "Number of entries in %s, %d, does not match expected number, %d.", \
		name,entries_read,siz);				\
  } while(0);

#define class_read_list_of_integers(name,destination,siz)			\
  do {									\
    class_call(parser_read_list_of_integers(pfc,name,			\
	&entries_read,&(destination),&flag1,errmsg),			\
	       errmsg,							\
	       errmsg);							\
    class_test(flag1 == _FALSE_,errmsg,					\
	"Entry %s is required but not found!",name)			\
        class_test(entries_read != siz,errmsg,			\
             "Number of entries in %s, %d, does not match expected number, %d.", \
		name,entries_read,siz);				\
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
		 struct transfers *ptr,
		 struct primordial *ppm,
		 struct spectra *psp,
		 struct nonlinear *pnl,
		 struct lensing *ple,
		 struct output *pop,
		 ErrorMsg errmsg
		 );

  int input_init(
		 struct file_content * pfc,
		 struct precision * ppr,
		 struct background *pba,
		 struct thermo *pth,
		 struct perturbs *ppt,
		 struct transfers *ptr,
		 struct primordial *ppm,
		 struct spectra *psp,
		 struct nonlinear *pnl,
		 struct lensing *ple,
		 struct output *pop,
		 ErrorMsg errmsg
		 );

  int input_default_params(
			   struct background *pba,
			   struct thermo *pth,
			   struct perturbs *ppt,
			   struct transfers *ptr,
			   struct primordial *ppm,
			   struct spectra *psp,
			   struct nonlinear *pnl,
			   struct lensing *ple,
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
