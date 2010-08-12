/** @file output.h Documented includes for output module */

#ifndef __OUTPUT__
#define __OUTPUT__

#include "spectra.h"

struct output {

  FileArg cl;
  FileArg cls;
  FileArg clt;
  FileArg cls_ad;
  FileArg cls_bi;
  FileArg cls_cdi;
  FileArg cls_nid;
  FileArg cls_niv;
  FileArg cls_ad_bi;
  FileArg cls_ad_cdi;
  FileArg cls_ad_nid;
  FileArg cls_ad_niv;
  FileArg cls_bi_cdi;
  FileArg cls_bi_nid;
  FileArg cls_bi_niv;
  FileArg cls_cdi_nid;
  FileArg cls_cdi_niv;
  FileArg cls_nid_niv;

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
  
#ifdef __cplusplus
}
#endif

#endif
