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
  int lmax;

  if (ppt->has_perturbations == _FALSE_) {
    if (pop->output_verbose > 0)
      printf("No spectra requested. Output module skipped.\n");
    return _SUCCESS_;
  }
  else {
    if (pop->output_verbose > 0)
      printf("Writing in output files \n");
  }

  /* deal with all C_l's */

  if (ppt->has_cls == _TRUE_) {

    /* loop over modes */

    for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {

      /* create bunch of files for this mode and all initial conditions */

      class_alloc(out,ppt->ic_size[index_mode]*sizeof(FILE *),pop->error_message);

      /* create temporary vector for writing the C_l's */

      class_alloc(cl_output,ppt->ic_size[index_mode]*psp->ct_size*sizeof(double),pop->error_message);

      /* give its name to each file (for each mode and initial conditions) */

      for (index_ic = 0; index_ic < ppt->ic_size[index_mode]; index_ic++) {

	if ((ppt->has_scalars) && (index_mode == ppt->index_md_scalars)) {

	  if ((ppt->has_ad) && (index_ic == ppt->index_ic_ad)) {

	    class_open(out[index_ic],pop->cls_ad,"w",pop->error_message);

	    fprintf(out[index_ic],"# dimensionless [l(l+1)/2pi] C_l for scalar adiabatic mode\n");
	  }

	  if ((ppt->has_bi) && (index_ic == ppt->index_ic_bi)) {

	    class_open(out[index_ic],pop->cls_bi,"w",pop->error_message);

	    fprintf(out[index_ic],"# dimensionless [l(l+1)/2pi] C_l for scalar baryon isocurvature mode\n");
	  }

	  if ((ppt->has_cdi) && (index_ic == ppt->index_ic_cdi)) {

	    class_open(out[index_ic],pop->cls_cdi,"w",pop->error_message);

	    fprintf(out[index_ic],"# dimensionless [l(l+1)/2pi] C_l for scalar CDM isocurvature mode\n");
	  }

	  if ((ppt->has_nid) && (index_ic == ppt->index_ic_nid)) {

	    class_open(out[index_ic],pop->cls_nid,"w",pop->error_message);

	    fprintf(out[index_ic],"# dimensionless [l(l+1)/2pi] C_l for scalar neutrino density isocurvature mode\n");
	  }

	  if ((ppt->has_niv) && (index_ic == ppt->index_ic_niv)) {

	    class_open(out[index_ic],pop->cls_niv,"w",pop->error_message);

	    fprintf(out[index_ic],"# dimensionless [l(l+1)/2pi] C_l for scalar neutrino velocity isocurvature mode\n");
	  }

	}

	if ((ppt->has_tensors) && (index_mode == ppt->index_md_tensors)) {
	
	  class_open(out[index_ic],pop->clt,"w",pop->error_message);

	  fprintf(out[index_ic],"# dimensionless [l(l+1)/2pi] C_l for tensor mode\n");

	}

	/* write headers */

	fprintf(out[index_ic],"# number of values of l:\n");
	fprintf(out[index_ic],"%d\n",(int)(psp->l_max[index_mode]-1));
        fprintf(out[index_ic],"#  l ");

	if (psp->has_tt == _TRUE_)
	  fprintf(out[index_ic],"TT           ");
	if (psp->has_ee == _TRUE_)
	  fprintf(out[index_ic],"EE           ");
	if (psp->has_te == _TRUE_)
	  fprintf(out[index_ic],"TE           "); 
	if (psp->has_bb == _TRUE_)
	  fprintf(out[index_ic],"BB           ");
	if (psp->has_pp == _TRUE_)
	  fprintf(out[index_ic],"phiphi       ");
	if (psp->has_tp == _TRUE_)
	  fprintf(out[index_ic],"Tphi         ");

	fprintf(out[index_ic],"\n");

      }

      /* loop over l */

      for (l = 2; l <= psp->l_max[index_mode]; l++) {  

	/* get results for one given mode (but all ic's and types) */

	class_call(spectra_cl_at_l(psp,index_mode,(double)l,cl_output),
		   psp->error_message,
		   pop->error_message);

	/* write result in corresponding file for each ic */

	for (index_ic = 0; index_ic < ppt->ic_size[index_mode]; index_ic++) {

	  fprintf(out[index_ic],"%4d",l);
	  for (index_ct=0; index_ct < psp->ct_size; index_ct++) {
	    fprintf(out[index_ic]," ",l);
	    fprintf(out[index_ic],"%e",l*(l+1)/2./_PI_*cl_output[index_ic * psp->ct_size + index_ct]);
	  }
	  fprintf(out[index_ic],"\n");	
	  
	}

      }

      /* close files */

      for (index_ic = 0; index_ic < ppt->ic_size[index_mode]; index_ic++) {

	fclose(out[index_ic]);

      }

      free(out);
      free(cl_output);

    }

    /* if more than one mode or more than one initial condition, output also the total C_l's summed over modes and ic's */

    if ((psp->md_size > 1) || ((psp->md_size == 1) && (psp->ic_size[0] > 1))) {

      class_alloc(cl_output,psp->ct_size*sizeof(double),pop->error_message);

      class_open(outbis,pop->cltot,"w",pop->error_message);

      lmax=0;
      for (index_mode=0; index_mode<psp->md_size; index_mode++) {
	lmax=max(lmax,psp->l_max[index_mode]);
      }

      fprintf(outbis,"# dimensionless [l(l+1)/2pi] C_l summed over modes and initial conditions \n");

      fprintf(outbis,"# number of values of l:\n");
      fprintf(outbis,"%d\n",(int)(lmax-1));
      fprintf(outbis,"#  l ");

      if (psp->has_tt == _TRUE_)
	fprintf(outbis,"TT           ");
      if (psp->has_ee == _TRUE_)
	fprintf(outbis,"EE           ");
      if (psp->has_te == _TRUE_)
	fprintf(outbis,"TE           "); 
      if (psp->has_bb == _TRUE_)
	fprintf(outbis,"BB           ");
      if (psp->has_pp == _TRUE_)
	fprintf(outbis,"phiphi       ");
      if (psp->has_tp == _TRUE_)
	fprintf(outbis,"Tphi         ");
      fprintf(outbis,"\n");

      for (l = 2; l <= lmax; l++) {  

	class_call(spectra_cl_tot_at_l(psp,(double)l,cl_output),
		   psp->error_message,
		   pop->error_message);

	fprintf(outbis,"%4d",l);
	for (index_ct=0; index_ct < psp->ct_size; index_ct++) {
	  fprintf(outbis," ");
	  fprintf(outbis,"%e",l*(l+1)/2./_PI_*cl_output[index_ct]);
	}
	fprintf(outbis,"\n");	
	  
      }

      fclose(outbis);

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
