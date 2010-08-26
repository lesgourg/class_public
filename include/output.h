/** @file output.h Documented includes for output module */

#ifndef __OUTPUT__
#define __OUTPUT__

#include "spectra.h"

/**
 * Structure containing various informations on the output format, 
 * initialized in input_init() in the input module
 *
 */

struct output {

   /** @name - root for all file names */

  //@{

  FileArg root;

  //@}

  /** @name - number and value(s) of redshift at which P(k,z) should be written */

  //@{

  double z_pk_num;
  double * z_pk;

  //@}

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
		  struct background * pba,
		  struct perturbs * ppt,
		  struct spectra * psp,
		  struct output * pop
		  );

  int output_cl(
		struct perturbs * ppt,
		struct spectra * psp,
		struct output * pop
		);

  int output_pk(
		struct background * pba,
		struct perturbs * ppt,
		struct spectra * psp,
		struct output * pop
		);

  int output_open_cl_file(
			  struct spectra * psp,
			  struct output * pop,
			  FILE * * clfile,
			  FileArg filename,
			  char * first_line,
			  int lmax
			  );

  int output_one_line_of_cl(
			    FILE * clfile,
			    double l,
			    double * cl,
			    int ct_size
			    );

  int output_open_pk_file(
			  struct spectra * psp,
			  struct output * pop,
			  FILE * * clfile,
			  FileArg filename,
			  char * first_line,
			  double z
			  );

  int output_one_line_of_pk(
			    FILE * clfile,
			    double one_k,
			    double one_pk
			    );
#ifdef __cplusplus
}
#endif

#endif
