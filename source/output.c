/** @file output.c Documented output module
 * Julien Lesgourgues, 18.04.2010    
 *
 * This module writes the output in files.
 */

#include "output.h"

int output_init(
		struct precision * ppr,
		struct background * pba,
		struct perturbs * ppt,
		struct transfers * ptr,
		struct spectra * psp,
		struct output * pop
		) {

  FILE * * out;
  FILE * outbis;
  int index_mode,index_ic,index_ct,l,index_k;
  double * cl_output;
  double * pk_output;

  if (ppt->tp_size == NULL) {
    if (pop->output_verbose > 0)
      printf("No spectra requested. Output module skipped.\n");
    return _SUCCESS_;
  }
  else {
    if (pop->output_verbose > 0)
      printf("Writing in output files \n");
  }

  /* deal with all C_l's */

  if (ptr->tt_size != NULL) {

    for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {

      class_alloc(out,ppt->ic_size[index_mode]*sizeof(FILE *),pop->error_message);

      class_alloc(cl_output,ppt->ic_size[index_mode]*psp->ct_size*sizeof(double),pop->error_message);

      for (index_ic = 0; index_ic < ppt->ic_size[index_mode]; index_ic++) {

	if ((ppt->has_scalars) && (index_mode == ppt->index_md_scalars)) {

	  if ((ppt->has_ad) && (index_ic == ppt->index_ic_ad)) {

	    class_open(out[index_ic],pop->cls_ad,"w",pop->error_message);

	    fprintf(out[index_ic],"# dimensionless [l(l+1)/2pi] C_l for scalar adiabatic mode\n");
	  }
	  else {
	    class_test(0==0,
		       pop->error_message,
		       "coding isocurvature modes not finished");
	  }	  
	}

	if ((ppt->has_tensors) && (index_mode == ppt->index_md_tensors)) {
	
	  class_open(out[index_ic],pop->clt,"w",pop->error_message);

	  fprintf(out[index_ic],"# dimensionless [l(l+1)/2pi] C_l for tensor mode\n");

	}

	fprintf(out[index_ic],"# number of values of l:\n");
	fprintf(out[index_ic],"%d\n",(int)(psp->l[index_mode][psp->l_size[index_mode]-1]-1));
        fprintf(out[index_ic],"#  l ");

	if (psp->has_tt == _TRUE_)
	  fprintf(out[index_ic],"TT           ");
	if (psp->has_ee == _TRUE_)
	  fprintf(out[index_ic],"EE           ");
	if (psp->has_te == _TRUE_)
	  fprintf(out[index_ic],"TE            ");
	if (psp->has_bb == _TRUE_)
	  fprintf(out[index_ic],"BB             ");
	if (psp->has_pp == _TRUE_)
	  fprintf(out[index_ic],"phiphi       ");
	if (psp->has_tp == _TRUE_)
	  fprintf(out[index_ic],"Tphi           ");

/* 	for (index_ct=0; index_ct < psp->ct_size[index_mode]; index_ct++) { */
/* 	  if ((ppt->has_cl_cmb_temperature == _TRUE_) && */
/* 	      (index_ct == psp->index_ct_tt)) fprintf(out[index_ic],"TT           "); */
/* 	  if ((ppt->has_cl_cmb_polarization == _TRUE_) && */
/* 	      (index_ct == psp->index_ct_ee)) fprintf(out[index_ic],"EE           "); */
/* 	  if ((ppt->has_cl_cmb_temperature == _TRUE_) && */
/* 	      (ppt->has_cl_cmb_polarization == _TRUE_) && */
/* 	      (index_ct == psp->index_ct_te)) fprintf(out[index_ic],"TE            "); */
/* 	  if ((ppt->has_scalars) && (index_mode == ppt->index_md_scalars)) { */
/* 	    if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && */
/* 		(index_ct == psp->index_ct_pp)) fprintf(out[index_ic],"phiphi       "); */
/* 	    if ((ppt->has_cl_cmb_temperature == _TRUE_) && */
/* 		(ppt->has_cl_cmb_lensing_potential == _TRUE_) && */
/* 		(index_ct == psp->index_ct_tp)) fprintf(out[index_ic],"Tphi           "); */
/* 	  } */
/* 	  if ((ppt->has_tensors) && (index_mode == ppt->index_md_tensors)) { */
/* 	    if ((ppt->has_cl_cmb_polarization == _TRUE_) && */
/* 		(index_ct == psp->index_ct_bb)) fprintf(out[index_ic],"BB           "); */
/* 	  } */
/* 	} */

	fprintf(out[index_ic],"\n");

      }

      for (l = 2; l <= psp->l[index_mode][psp->l_size[index_mode]-1]; l++) {  

	class_call(spectra_cl_at_l(psp,index_mode,(double)l,cl_output),
		   psp->error_message,
		   pop->error_message);

	for (index_ic = 0; index_ic < ppt->ic_size[index_mode]; index_ic++) {

	  fprintf(out[index_ic],"%4d",l);
	  for (index_ct=0; index_ct < psp->ct_size; index_ct++) {
	    fprintf(out[index_ic]," ",l);
	    fprintf(out[index_ic],"%e",l*(l+1)/2./_PI_*cl_output[index_ic * psp->ct_size + index_ct]);
	  }
	  fprintf(out[index_ic],"\n");	
	  
	}

      }

      for (index_ic = 0; index_ic < ppt->ic_size[index_mode]; index_ic++) {

	fclose(out[index_ic]);

      }

      free(out);
      free(cl_output);

    }
  }

  /* deal with all Fourier space spectra */

  if (ppt->has_pk_matter == _TRUE_) {
    
    class_test((pop->z_pk > psp->z_max_pk),
	       pop->error_message,
	       "P(k,z) computed up to z=%f but requested at z=%f. Must increase z_max_pk in precision file.",psp->z_max_pk,pop->z_pk);

    index_mode=ppt->index_md_scalars;

    /* if z_pk = 0, no interpolation needed, 
       just let pk_output point to the right address */
    if (pop->z_pk == 0) {
      pk_output=&(psp->pk[(psp->eta_size-1) * ppt->ic_size[index_mode] * psp->k_size]);
    }

    /* if 0 <= z_pk <= z_max_pk, interpolation needed, */
    else {
      class_alloc(pk_output,sizeof(double)*ppt->ic_size[index_mode]*psp->k_size,pop->error_message);
      class_call(spectra_pk_at_z(pba,psp,ppt->index_md_scalars,pop->z_pk,pk_output),
		 psp->error_message,
		 pop->error_message);
    }
      
    class_open(outbis,pop->pk,"w",pop->error_message);

    fprintf(outbis,"# Matter power spectrum P(k) at redshift z=%f\n",pop->z_pk);
    fprintf(outbis,"# Number of wavenumbers k:\n");
    fprintf(outbis,"%d\n",psp->k_size);
    fprintf(outbis,"# k (h/Mpc)  P (Mpc/h)^3:\n");
    /* if isocurvature modes present, should improve these preliminary lines */

    for (index_k=0; index_k<psp->k_size; index_k++) {

      fprintf(outbis,"%e",psp->k[index_k]/pba->h);

      for (index_ic = 0; index_ic < ppt->ic_size[index_mode]; index_ic++) {
	fprintf(outbis," %e",
		pow(pba->h,3)*pk_output[index_ic * psp->k_size + index_k]);
      }
      fprintf(outbis,"\n");
    }

    fclose(outbis);

    if (pop->z_pk != 0.)
      free(pk_output);

  }

  return _SUCCESS_;

}
