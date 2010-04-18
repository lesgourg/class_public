/** @file output.c Documented output module
 * Julien Lesgourgues, 18.04.2010    
 *
 * This module writes the output in files.
 */

#include "output.h"

int output_init(
		struct perturbs * ppt_input,
		struct transfers * ptr_input,
		struct spectra * psp_input,
		struct output * pop
		) {

  FILE * * out;
  int index_mode,index_ic,index_ct,l;
  double *cl;

  if (pop->output_verbose > 0)
    printf("Writing in output files \n");

  if (ptr_input->has_cls == _TRUE_) {

    for (index_mode = 0; index_mode < ppt_input->md_size; index_mode++) {

      out = malloc(ppt_input->ic_size[index_mode]*sizeof(FILE *));
      if (out == NULL) {
	sprintf(pop->error_message,"%s(L:%d) : Could not allocate out",__func__,__LINE__);
	return _FAILURE_;
      }

      cl = malloc(ppt_input->ic_size[index_mode]*sizeof(double));
      if (cl == NULL) {
	sprintf(pop->error_message,"%s(L:%d) : Could not allocate cl",__func__,__LINE__);
	return _FAILURE_;
      }

      for (index_ic = 0; index_ic < ppt_input->ic_size[index_mode]; index_ic++) {

	if ((ppt_input->has_scalars) && (index_mode == ppt_input->index_md_scalars)) {
	  if ((ppt_input->has_ad) && (index_ic == ppt_input->index_ic_ad)) {
	    out[index_ic]=fopen(pop->cls_ad,"w");
	    if (out[index_ic] == NULL) {
	      sprintf(pop->error_message,"%s(L:%d) : Could not open ou[index_ic]",__func__,__LINE__);
	      return _FAILURE_;
	    }
	  }
	}
      }

      for (l = 2; l <= psp_input->l[index_mode][psp_input->l_size[index_mode]-1]; l++) {  

	if (spectra_cl_at_l((double)l,index_mode,cl) == _FAILURE_) {
	  sprintf(pop->error_message,"%s(L:%d) : error in spectra_cl_at_l()\n=>%s",__func__,__LINE__,psp_input->error_message);
	  return _FAILURE_;
	}

	for (index_ic = 0; index_ic < ppt_input->ic_size[index_mode]; index_ic++) {

	  fprintf(out[index_ic],"%d",l);
	  for (index_ct=0; index_ct < psp_input->ct_size; index_ct++) {
	    fprintf(out[index_ic]," ",l);
	    fprintf(out[index_ic],"%e",l*(l+1)*cl[index_ic * psp_input->ct_size + index_ct]);
	  }
	  fprintf(out[index_ic],"\n");	
	  
	}

      }

      for (index_ic = 0; index_ic < ppt_input->ic_size[index_mode]; index_ic++) {

	fclose(out[index_ic]);

      }

      free(out);
      free(cl);

    }
  }

  return _SUCCESS_;

}
