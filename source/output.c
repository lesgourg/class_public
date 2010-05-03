/** @file output.c Documented output module
 * Julien Lesgourgues, 18.04.2010    
 *
 * This module writes the output in files.
 */

#include "output.h"

int output_init(
		struct precision * ppr_input,
		struct background * pba_input,
		struct perturbs * ppt_input,
		struct transfers * ptr_input,
		struct spectra * psp_input,
		struct output * pop
		) {

  FILE * * out;
  FILE * outbis;
  int index_mode,index_ic,index_ct,l,index_k;
  double * cl_output;
  double * pk_output;

  if (pop->output_verbose > 0)
    printf("Writing in output files \n");

  if (ptr_input->tt_size >0) {

    for (index_mode = 0; index_mode < ppt_input->md_size; index_mode++) {

      out = malloc(ppt_input->ic_size[index_mode]*sizeof(FILE *));
      if (out == NULL) {
	sprintf(pop->error_message,"%s(L:%d) : Could not allocate out",__func__,__LINE__);
	return _FAILURE_;
      }

      cl_output= malloc(ppt_input->ic_size[index_mode]*psp_input->ct_size*sizeof(double));

      if (cl_output == NULL) {
	sprintf(pop->error_message,"%s(L:%d) : Could not allocate cl_output",__func__,__LINE__);
	return _FAILURE_;
      }

      for (index_ic = 0; index_ic < ppt_input->ic_size[index_mode]; index_ic++) {

	if ((ppt_input->has_scalars) && (index_mode == ppt_input->index_md_scalars)) {
	  if ((ppt_input->has_ad) && (index_ic == ppt_input->index_ic_ad)) {
	    out[index_ic]=fopen(pop->cls_ad,"w");
	    if (out[index_ic] == NULL) {
	      sprintf(pop->error_message,"%s(L:%d) : Could not open out[index_ic]",__func__,__LINE__);
	      return _FAILURE_;
	    }
	    fprintf(out[index_ic],"# dimensionless [l(l+1)/2pi] C_l for adiabatic mode\n");
	  }
	}
	
	fprintf(out[index_ic],"# number of values of l:\n");
	fprintf(out[index_ic],"%d\n",(int)(psp_input->l[index_mode][psp_input->l_size[index_mode]-1]-1));
        fprintf(out[index_ic],"#  l ");
	for (index_ct=0; index_ct < psp_input->ct_size; index_ct++) {
	  if ((ppt_input->has_cl_cmb_temperature == _TRUE_) &&
	      (index_ct == psp_input->index_ct_tt)) fprintf(out[index_ic],"TT           ");
	  if ((ppt_input->has_cl_cmb_polarization == _TRUE_) &&
	      (index_ct == psp_input->index_ct_ee)) fprintf(out[index_ic],"EE           ");
	  if ((ppt_input->has_cl_cmb_temperature == _TRUE_) &&
	      (ppt_input->has_cl_cmb_polarization == _TRUE_) &&
	      (index_ct == psp_input->index_ct_te)) fprintf(out[index_ic],"TE            ");
	  if ((ppt_input->has_cl_cmb_lensing_potential == _TRUE_) &&
	      (index_ct == psp_input->index_ct_pp)) fprintf(out[index_ic],"phiphi       ");
	  if ((ppt_input->has_cl_cmb_temperature == _TRUE_) &&
	      (ppt_input->has_cl_cmb_lensing_potential == _TRUE_) &&
	      (index_ct == psp_input->index_ct_tp)) fprintf(out[index_ic],"Tphi           ");
	}
	fprintf(out[index_ic],"\n");

      }

      for (l = 2; l <= psp_input->l[index_mode][psp_input->l_size[index_mode]-1]; l++) {  

	if (spectra_cl_at_l((double)l,index_mode,cl_output) == _FAILURE_) {
	  sprintf(pop->error_message,"%s(L:%d) : error in spectra_cl_at_l()\n=>%s",__func__,__LINE__,psp_input->error_message);
	  return _FAILURE_;
	}

	for (index_ic = 0; index_ic < ppt_input->ic_size[index_mode]; index_ic++) {

	  fprintf(out[index_ic],"%4d",l);
	  for (index_ct=0; index_ct < psp_input->ct_size; index_ct++) {
	    fprintf(out[index_ic]," ",l);
	    fprintf(out[index_ic],"%e",l*(l+1)/2./_PI_*cl_output[index_ic * psp_input->ct_size + index_ct]);
	  }
	  fprintf(out[index_ic],"\n");	
	  
	}

      }

      for (index_ic = 0; index_ic < ppt_input->ic_size[index_mode]; index_ic++) {

	fclose(out[index_ic]);

      }

      free(out);
      free(cl_output);

    }
  }

  if (ppt_input->has_pk_matter == _TRUE_) {
    
    if (pop->z_pk > ppr_input->z_max_pk) {
      sprintf(pop->error_message,"%s(L:%d) : P(k,z) computed up to z=%f but requested at z=%f. Must increase z_max_pk in precision file.",__func__,__LINE__,ppr_input->z_max_pk,pop->z_pk);
      return _FAILURE_;
    }

    index_mode=ppt_input->index_md_scalars;

    /* if z_pk = 0, no interpolation needed, 
       just let pk_output point to the right address */
    if (pop->z_pk == 0) {
      pk_output=&(psp_input->pk[(psp_input->eta_size-1) * ppt_input->ic_size[index_mode] * psp_input->k_size]);
    }

    /* if 0 <= z_pk <= z_max_pk, interpolation needed, */
    else {
      pk_output = malloc(sizeof(double)*ppt_input->ic_size[index_mode]*psp_input->k_size);
      if (pk_output == NULL) {
	sprintf(pop->error_message,"%s(L:%d) : Could not allocate pk_output",__func__,__LINE__);
	return _FAILURE_;
      }
      if (spectra_pk_at_z(pop->z_pk,pk_output) == _FAILURE_) {
	sprintf(pop->error_message,"%s(L:%d) : error in spectra_pk_at_z()\n=>%s",__func__,__LINE__,psp_input->error_message);
	return _FAILURE_;
      }
		      
    }
      
    outbis=fopen(pop->pk,"w");

    fprintf(outbis,"# Matter power spectrum P(k) at redshift z=%f\n",pop->z_pk);
    fprintf(outbis,"# Number of wavenumbers k:\n");
    fprintf(outbis,"%d\n",psp_input->k_size);
    fprintf(outbis,"# k (h/Mpc)  P (Mpc/h)^3:\n");
    /* if isocurvature modes present, should improve these preliminary lines */

    for (index_k=0; index_k<psp_input->k_size; index_k++) {

      fprintf(outbis,"%e",psp_input->k[index_k]/pba_input->h);

      for (index_ic = 0; index_ic < ppt_input->ic_size[index_mode]; index_ic++) {
	fprintf(outbis," %e",
		pow(pba_input->h,3)*pk_output[index_ic * psp_input->k_size + index_k]);
      }
      fprintf(outbis,"\n");
    }

    fclose(outbis);

    if (pop->z_pk != 0.)
      free(pk_output);

  }

  return _SUCCESS_;

}
