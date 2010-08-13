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

  FILE * * * out_md_ic;  /* out_md_ic[index_mode][index_ic1_ic2] */
  FILE * * out_md;       /* out_md[index_mode] (cls summed eventually over ic's) */
  FILE * out;            /* (cls summed eventually over modes and ic's) */
  double * * cl_md_ic;   /* cl_md_ic[index_mode][index_ic1_ic2*psp->ct_size+index_ct] */
  double * * cl_md;      /* cl_md[index_mode][index_ct] */
  double * cl_tot;           /* cl_tot[index_ct] */

  FILE * outbis;

  int index_mode,index_ic1,index_ic2,index_ic1_ic2,index_ct,l,index_k;

  double * pk_output;

  FileArg file_name;
  char * first_line;

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

    /* first, allocate all arrays of files and cls */

    class_alloc(out_md_ic,
		psp->md_size*sizeof(FILE * *),
		pop->error_message);
        
    class_alloc(cl_md_ic,
		psp->md_size*sizeof(double *),
		pop->error_message);

    class_alloc(out_md,
		psp->md_size*sizeof(FILE *),
		pop->error_message);
    
    class_alloc(cl_md,
		psp->md_size*sizeof(double *),
		pop->error_message);

    for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {

      class_alloc(out_md_ic[index_mode],
		  psp->ic_ic_size[index_mode]*sizeof(FILE *),
		  pop->error_message);
      
    }

    /* second, open only the relevant files, and write a heading in each of them */

    class_call(output_open_cl_file(psp,
				   pop,
				   &out,
				   pop->cl,
				   "# dimensionless total [l(l+1)/2pi] C_l's",
				   psp->l_max_tot
				   ),
	       pop->error_message,
	       pop->error_message);
   
    class_alloc(cl_tot,
		psp->ct_size*sizeof(double),
		pop->error_message);

    if (ppt->md_size > 1) {

      for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {
	
	if ((ppt->has_scalars) && (index_mode == ppt->index_md_scalars)) {
	  
	  strcpy(file_name,pop->cls);
	  first_line="# dimensionless [l(l+1)/2pi] C_l's for scalar mode";

	}

	if ((ppt->has_tensors) && (index_mode == ppt->index_md_tensors)) {
	  	  
	  strcpy(file_name,pop->clt);
	  first_line="# dimensionless [l(l+1)/2pi] C_l's for tensor mode";

	}
	
	class_call(output_open_cl_file(psp,
				       pop,
				       &(out_md[index_mode]),
				       file_name,
				       first_line,
				       psp->l_max[index_mode]
				       ),
		   pop->error_message,
		   pop->error_message);
	
	class_alloc(cl_md[index_mode],
		    psp->ct_size*sizeof(double),
		    pop->error_message);

      }
    }

    for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {

      if (ppt->ic_size[index_mode] > 1) {

	for (index_ic1 = 0; index_ic1 < ppt->ic_size[index_mode]; index_ic1++) {
	  
	  for (index_ic2 = index_ic1; index_ic2 < ppt->ic_size[index_mode]; index_ic2++) {
	    
	    if ((ppt->has_scalars) && (index_mode == ppt->index_md_scalars)) {

	      if ((ppt->has_ad) && 
		  (index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_ad)) {

		strcpy(file_name,pop->cls_ad);
		first_line="# dimensionless [l(l+1)/2pi] C_l for scalar adiabatic (AD) mode";
	      }

	      if ((ppt->has_bi) && 
		  (index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_bi)) {

		strcpy(file_name,pop->cls_bi);
		first_line="# dimensionless [l(l+1)/2pi] C_l for scalar baryon isocurvature (BI) mode";
	      }
	      
	      if ((ppt->has_cdi) && 
		  (index_ic1 == ppt->index_ic_cdi) && (index_ic2 == ppt->index_ic_cdi)) {

		strcpy(file_name,pop->cls_cdi);
		first_line="# dimensionless [l(l+1)/2pi] C_l for scalar CDM isocurvature (CDI) mode";
	      }

	      if ((ppt->has_nid) && 
		  (index_ic1 == ppt->index_ic_nid) && (index_ic2 == ppt->index_ic_nid)) {

		strcpy(file_name,pop->cls_nid);
		first_line="# dimensionless [l(l+1)/2pi] C_l for scalar CDM isocurvature (CDI) mode";
	      }

	      if ((ppt->has_niv) && 
		  (index_ic1 == ppt->index_ic_niv) && (index_ic2 == ppt->index_ic_niv)) {

		strcpy(file_name,pop->cls_niv);
		first_line="# dimensionless [l(l+1)/2pi] C_l for scalar neutrino velocity isocurvature (NIV) mode";
	      }

	      if ((ppt->has_ad) && 
		  (ppt->has_bi) && (index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_bi)) {

		strcpy(file_name,pop->cls_ad_bi);
		first_line="# dimensionless [l(l+1)/2pi] C_l for scalar cross ADxBI mode";
	      }
	      
	      if ((ppt->has_ad) && (ppt->has_cdi) && 
		  (index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_cdi)) {

		strcpy(file_name,pop->cls_ad_cdi);
		first_line="# dimensionless [l(l+1)/2pi] C_l for scalar cross ADxCDI mode";
	      }

	      if ((ppt->has_ad) && (ppt->has_nid) && 
		  (index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_nid)) {

		strcpy(file_name,pop->cls_ad_nid);
		first_line="# dimensionless [l(l+1)/2pi] C_l for scalar cross ADxNID mode";
	      }

	      if ((ppt->has_ad) && (ppt->has_niv) && 
		  (index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_niv)) {

		strcpy(file_name,pop->cls_ad_niv);
		first_line="# dimensionless [l(l+1)/2pi] C_l for scalar cross ADxNIV mode";
	      }

	      if ((ppt->has_bi) && (ppt->has_cdi) && 
		  (index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_cdi)) {

		strcpy(file_name,pop->cls_bi_cdi);
		first_line="# dimensionless [l(l+1)/2pi] C_l for scalar cross BIxCDI mode";
	      }

	      if ((ppt->has_bi) && (ppt->has_nid) && 
		  (index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_nid)) {

		strcpy(file_name,pop->cls_bi_nid);
		first_line="# dimensionless [l(l+1)/2pi] C_l for scalar cross BIxNID mode";
	      }

	      if ((ppt->has_bi) && (ppt->has_niv) && 
		  (index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_niv)) {

		strcpy(file_name,pop->cls_bi_niv);
		first_line="# dimensionless [l(l+1)/2pi] C_l for scalar cross BIxNIV mode";
	      }

	      if ((ppt->has_cdi) && (ppt->has_nid) && 
		  (index_ic1 == ppt->index_ic_cdi) && (index_ic2 == ppt->index_ic_nid)) {

		strcpy(file_name,pop->cls_cdi_nid);
		first_line="# dimensionless [l(l+1)/2pi] C_l for scalar cross CDIxNID mode";
	      }

	      if ((ppt->has_cdi) && (ppt->has_niv) && 
		  (index_ic1 == ppt->index_ic_cdi) && (index_ic2 == ppt->index_ic_niv)) {

		strcpy(file_name,pop->cls_cdi_niv);
		first_line="# dimensionless [l(l+1)/2pi] C_l for scalar cross CDIxNIV mode";
	      }

	      if ((ppt->has_nid) && (ppt->has_niv) && 
		  (index_ic1 == ppt->index_ic_nid) && (index_ic2 == ppt->index_ic_niv)) {

		strcpy(file_name,pop->cls_nid_niv);
		first_line="# dimensionless [l(l+1)/2pi] C_l for scalar cross NIDxNIV mode";
	      }

	    }

	    if ((ppt->has_tensors) && (index_mode == ppt->index_md_tensors)) {

	      class_test(0==0,
			 pop->error_message,
			 "Seems that we have mixed initial conditions for tensors? Should not happen!\n");
	      
	    }

	    /* index value for the coefficients of the symmetric index_ic1*index_ic2 matrix; 
	       takes values between 0 and N(N+1)/2-1 with N=ppt->ic_size[index_mode] */
	    index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_mode]);

	    if ((index_ic1 == index_ic2) || (psp->is_non_zero[index_mode][index_ic1_ic2] == _TRUE_)) {

	      class_call(output_open_cl_file(psp,
					     pop,
					     &(out_md_ic[index_mode][index_ic1_ic2]),
					     file_name,
					     first_line,
					     psp->l_max[index_mode]
					     ),
			 pop->error_message,
			 pop->error_message);

	    }
	  }
	}

	class_alloc(cl_md_ic[index_mode],
			psp->ic_ic_size[index_mode]*psp->ct_size*sizeof(double),
			pop->error_message);
      }
    }

    for (l = 2; l <= psp->l_max_tot; l++) {  

      class_call(spectra_cl_at_l(psp,(double)l,cl_tot,cl_md,cl_md_ic),
		 psp->error_message,
		 pop->error_message);

      class_call(output_one_line_of_cl(out,l,cl_tot,psp->ct_size),
		 pop->error_message,
		 pop->error_message);

      if (ppt->md_size > 1) {
	for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {
	  if (l <= psp->l_max[index_mode]) {

	    class_call(output_one_line_of_cl(out_md[index_mode],l,cl_md[index_mode],psp->ct_size),
		       pop->error_message,
		       pop->error_message);
	  }
	}
      }
	  
      for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {
	if ((ppt->ic_size[index_mode] > 1) && (l <= psp->l_max[index_mode])) {
	  for (index_ic1 = 0; index_ic1 < ppt->ic_size[index_mode]; index_ic1++) {
	    for (index_ic2 = index_ic1; index_ic2 < ppt->ic_size[index_mode]; index_ic2++) {
	      /* index value for the coefficients of the symmetric index_ic1*index_ic2 matrix; 
		 takes values between 0 and N(N+1)/2-1 with N=ppt->ic_size[index_mode] */
	      index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_mode]);

	      if ((index_ic1 == index_ic2) || (psp->is_non_zero[index_mode][index_ic1_ic2] == _TRUE_)) {

		class_call(output_one_line_of_cl(out_md_ic[index_mode][index_ic1_ic2],l,&(cl_md_ic[index_mode][index_ic1_ic2*psp->ct_size]),psp->ct_size),
			   pop->error_message,
			   pop->error_message);
	      }
	    }
	  }
	}
      }
    }

    /* close files and free arrays of files and cls */

    for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {
      if (ppt->ic_size[index_mode] > 1) {
	for (index_ic1 = 0; index_ic1 < ppt->ic_size[index_mode]; index_ic1++) {
	  for (index_ic2 = index_ic1; index_ic2 < ppt->ic_size[index_mode]; index_ic2++) {
	    index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_mode]);
	    if ((index_ic1 == index_ic2) || (psp->is_non_zero[index_mode][index_ic1_ic2] == _TRUE_)) {
	      fclose(out_md_ic[index_mode][index_ic1_ic2]);
	    }
	  }
	}
	free(cl_md_ic[index_mode]);
      }
    }
    if (ppt->md_size > 1) {
      for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {
	fclose(out_md[index_mode]);
	free(cl_md[index_mode]);
      }
    }
    fclose(out);
    free(cl_tot);
    for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {
      free(out_md_ic[index_mode]);
    }
    free(out_md_ic);
    free(cl_md_ic);
    free(out_md);
    free(cl_md);

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

      for (index_ic1 = 0; index_ic1 < ppt->ic_size[index_mode]; index_ic1++) {
	fprintf(outbis," %e",
		pow(pba->h,3)*pk_output[index_ic1 * psp->k_size + index_k]);
      }
      fprintf(outbis,"\n");
    }

    fclose(outbis);

    if (pop->z_pk != 0.)
      free(pk_output);

  }

  return _SUCCESS_;

}

int output_open_cl_file(
			struct spectra * psp,
			struct output * pop,
			FILE * * clfile,
			FileArg filename,
			char * first_line,
			int lmax
			) {

  class_open(*clfile,filename,"w",pop->error_message);

  fprintf(*clfile,"%s\n",first_line); 
  fprintf(*clfile,"# for l=2 to %d, i.e. number of lines equal to\n",(int)lmax);
  fprintf(*clfile,"%d\n",(int)(lmax-1));
  fprintf(*clfile,"#  l ");

  if (psp->has_tt == _TRUE_)
    fprintf(*clfile,"TT           ");
  if (psp->has_ee == _TRUE_)
    fprintf(*clfile,"EE           ");
  if (psp->has_te == _TRUE_)
    fprintf(*clfile,"TE           "); 
  if (psp->has_bb == _TRUE_)
    fprintf(*clfile,"BB           ");
  if (psp->has_pp == _TRUE_)
    fprintf(*clfile,"phiphi       ");
  if (psp->has_tp == _TRUE_)
    fprintf(*clfile,"Tphi         ");
  fprintf(*clfile,"\n");

  return _SUCCESS_;

}

int output_one_line_of_cl(
			  FILE * clfile,
			  double l,
			  double * cl,
			  int ct_size
			  ) {
  int index_ct;

  fprintf(clfile,"%4d",(int)l);

  for (index_ct=0; index_ct < ct_size; index_ct++) {
    fprintf(clfile," ");
    fprintf(clfile,"%e",l*(l+1)/2./_PI_*cl[index_ct]);
  }
  fprintf(clfile,"\n");	
    
  return _SUCCESS_;
    
}
