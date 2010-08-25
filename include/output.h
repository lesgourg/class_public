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
  FileArg pk_ad;
  FileArg pk_bi;
  FileArg pk_cdi;
  FileArg pk_nid;
  FileArg pk_niv;
  FileArg pk_ad_bi;
  FileArg pk_ad_cdi;
  FileArg pk_ad_nid;
  FileArg pk_ad_niv;
  FileArg pk_bi_cdi;
  FileArg pk_bi_nid;
  FileArg pk_bi_niv;
  FileArg pk_cdi_nid;
  FileArg pk_cdi_niv;
  FileArg pk_nid_niv;

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
