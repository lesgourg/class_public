/** @file output.c Documented output module
 *
 * Julien Lesgourgues, 26.08.2010    
 *
 * This module writes the output in files.
 *
 * The following function can be called from other modules or from the main:
 *
 * -# output_init() (must be called after spectra_init())
 *
 * No memory needs to be deallocated after that, 
 * hence there is no output_free() routine like in other modules.
 */

#include "output.h"

/** 
 * This routine writes the output in files.
 *
 * 
 * @param pba Input: pointer to background structure (needed for calling spectra_pk_at_z())
 * @param ppt Input : pointer perturbation structure
 * @param psp Input : pointer to spectra structure
 * @param pop Input : pointer to output structure
 */

int output_init(
		struct background * pba,
		struct perturbs * ppt,
		struct spectra * psp,
		struct output * pop
		) {

  /** Summary: */

  /** - check that we really want to output at least one spectrum */

  if ((ppt->has_cls == _FALSE_) && (ppt->has_pk_matter == _FALSE_)) {
    if (pop->output_verbose > 0)
      printf("No spectra requested. Output module skipped.\n");
    return _SUCCESS_;
  }
  else {
    if (pop->output_verbose > 0)
      printf("Writing in output files \n");
  }

  /** - deal with all anisotropy power spectra C_l's */

  if (ppt->has_cls == _TRUE_) {

    class_call(output_cl(ppt,psp,pop),
	       pop->error_message,
	       pop->error_message);
  }

  /** - deal with all Fourier matter power spectra P(k)'s */

  if (ppt->has_pk_matter == _TRUE_) {

    class_call(output_pk(pba,ppt,psp,pop),
	       pop->error_message,
	       pop->error_message);
  }

  return _SUCCESS_;

}

/** 
 * This routines writes the output in files for anisotropy power spectra C_l's.
 *
 * @param ppt Input : pointer perturbation structure
 * @param psp Input : pointer to spectra structure
 * @param pop Input : pointer to output structure
 */

int output_cl(
	      struct perturbs * ppt,
	      struct spectra * psp,
	      struct output * pop
	      ) {

  /** Summary: */

  /** - define local variables */

  FILE *** out_md_ic; /* array of pointers to files with argument 
			 out_md_ic[index_mode][index_ic1_ic2] 
			 (will contain cl's for each mode and pairs of initial conditions) */

  FILE ** out_md;     /* array of pointers to files with argument 
			 out_md[index_mode] 
			 (will contain cl's for each mode, summed eventually over ic's) */

  FILE * out;         /* (will contain total cl's, summed eventually over modes and ic's) */

  double ** cl_md_ic; /* array with argument 
			 cl_md_ic[index_mode][index_ic1_ic2*psp->ct_size+index_ct] */

  double ** cl_md;    /* array with argument 
			 cl_md[index_mode][index_ct] */

  double * cl_tot;    /* array with argument 
			 cl_tot[index_ct] */

  int index_mode;
  int index_ic1,index_ic2,index_ic1_ic2;
  int l;

  FileName file_name;
  char first_line[_LINE_LENGTH_MAX_];

  /** - first, allocate all arrays of files and cls */

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

  /** - second, open only the relevant files, and write a heading in each of them */

  sprintf(file_name,"%s%s",pop->root,"cl.dat");

  class_call(output_open_cl_file(psp,
				 pop,
				 &out,
				 file_name,
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
	  
	sprintf(file_name,"%s%s",pop->root,"cls.dat");
	strcpy(first_line,"# dimensionless [l(l+1)/2pi] C_l's for scalar mode");

      }

      if ((ppt->has_tensors) && (index_mode == ppt->index_md_tensors)) {
	  	  
	sprintf(file_name,"%s%s",pop->root,"clt.dat");
	strcpy(first_line,"# dimensionless [l(l+1)/2pi] C_l's for tensor mode");

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

	      sprintf(file_name,"%s%s",pop->root,"cls_ad.dat");
	      strcpy(first_line,"# dimensionless [l(l+1)/2pi] C_l for scalar adiabatic (AD) mode");
	    }

	    if ((ppt->has_bi) && 
		(index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_bi)) {

	      sprintf(file_name,"%s%s",pop->root,"cls_bi.dat");
	      strcpy(first_line,"# dimensionless [l(l+1)/2pi] C_l for scalar baryon isocurvature (BI) mode");
	    }
	      
	    if ((ppt->has_cdi) && 
		(index_ic1 == ppt->index_ic_cdi) && (index_ic2 == ppt->index_ic_cdi)) {

	      sprintf(file_name,"%s%s",pop->root,"cls_cdi.dat");
	      strcpy(first_line,"# dimensionless [l(l+1)/2pi] C_l for scalar CDM isocurvature (CDI) mode");
	    }

	    if ((ppt->has_nid) && 
		(index_ic1 == ppt->index_ic_nid) && (index_ic2 == ppt->index_ic_nid)) {

	      sprintf(file_name,"%s%s",pop->root,"cls_nid.dat");
	      strcpy(first_line,"# dimensionless [l(l+1)/2pi] C_l for scalar neutrino density isocurvature (NID) mode");
	    }

	    if ((ppt->has_niv) && 
		(index_ic1 == ppt->index_ic_niv) && (index_ic2 == ppt->index_ic_niv)) {

	      sprintf(file_name,"%s%s",pop->root,"cls_niv.dat");
	      strcpy(first_line,"# dimensionless [l(l+1)/2pi] C_l for scalar neutrino velocity isocurvature (NIV) mode");
	    }

	    if ((ppt->has_ad) && 
		(ppt->has_bi) && (index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_bi)) {

	      sprintf(file_name,"%s%s",pop->root,"cls_ad_bi.dat");
	      strcpy(first_line,"# dimensionless [l(l+1)/2pi] C_l for scalar cross ADxBI mode");
	    }
	      
	    if ((ppt->has_ad) && (ppt->has_cdi) && 
		(index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_cdi)) {

	      sprintf(file_name,"%s%s",pop->root,"cls_ad_cdi.dat");
	      strcpy(first_line,"# dimensionless [l(l+1)/2pi] C_l for scalar cross ADxCDI mode");
	    }

	    if ((ppt->has_ad) && (ppt->has_nid) && 
		(index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_nid)) {

	      sprintf(file_name,"%s%s",pop->root,"cls_ad_nid.dat");
	      strcpy(first_line,"# dimensionless [l(l+1)/2pi] C_l for scalar cross ADxNID mode");
	    }

	    if ((ppt->has_ad) && (ppt->has_niv) && 
		(index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_niv)) {

	      sprintf(file_name,"%s%s",pop->root,"cls_ad_niv.dat");
	      strcpy(first_line,"# dimensionless [l(l+1)/2pi] C_l for scalar cross ADxNIV mode");
	    }

	    if ((ppt->has_bi) && (ppt->has_cdi) && 
		(index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_cdi)) {

	      sprintf(file_name,"%s%s",pop->root,"cls_bi_cdi.dat");
	      strcpy(first_line,"# dimensionless [l(l+1)/2pi] C_l for scalar cross BIxCDI mode");
	    }

	    if ((ppt->has_bi) && (ppt->has_nid) && 
		(index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_nid)) {

	      sprintf(file_name,"%s%s",pop->root,"cls_bi_nid.dat");
	      strcpy(first_line,"# dimensionless [l(l+1)/2pi] C_l for scalar cross BIxNID mode");
	    }

	    if ((ppt->has_bi) && (ppt->has_niv) && 
		(index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_niv)) {

	      sprintf(file_name,"%s%s",pop->root,"cls_bi_niv.dat");
	      strcpy(first_line,"# dimensionless [l(l+1)/2pi] C_l for scalar cross BIxNIV mode");
	    }

	    if ((ppt->has_cdi) && (ppt->has_nid) && 
		(index_ic1 == ppt->index_ic_cdi) && (index_ic2 == ppt->index_ic_nid)) {

	      sprintf(file_name,"%s%s",pop->root,"cls_cdi_nid.dat");
	      strcpy(first_line,"# dimensionless [l(l+1)/2pi] C_l for scalar cross CDIxNID mode");
	    }

	    if ((ppt->has_cdi) && (ppt->has_niv) && 
		(index_ic1 == ppt->index_ic_cdi) && (index_ic2 == ppt->index_ic_niv)) {

	      sprintf(file_name,"%s%s",pop->root,"cls_cdi_niv.dat");
	      strcpy(first_line,"# dimensionless [l(l+1)/2pi] C_l for scalar cross CDIxNIV mode");
	    }

	    if ((ppt->has_nid) && (ppt->has_niv) && 
		(index_ic1 == ppt->index_ic_nid) && (index_ic2 == ppt->index_ic_niv)) {

	      sprintf(file_name,"%s%s",pop->root,"cls_nid_niv.dat");
	      strcpy(first_line,"# dimensionless [l(l+1)/2pi] C_l for scalar cross NIDxNIV mode");
	    }

	  }

	  if ((ppt->has_tensors) && (index_mode == ppt->index_md_tensors)) {

	    class_test(0==0,
		       pop->error_message,
		       "Seems that we have mixed initial conditions for tensors? Should not happen!\n");
	      
	  }

	  index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_mode]);

	  if (psp->is_non_zero[index_mode][index_ic1_ic2] == _TRUE_) {

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

  /** - third, perfomr loop over l. For each multipole, get all C_l's
      by calling spectra_cl_at_l() and distribute the results to 
      relevant files */

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
	for (index_ic1_ic2 = 0; index_ic1_ic2 < psp->ic_ic_size[index_mode]; index_ic1_ic2++) {
	  if (psp->is_non_zero[index_mode][index_ic1_ic2] == _TRUE_) {

	    class_call(output_one_line_of_cl(out_md_ic[index_mode][index_ic1_ic2],l,&(cl_md_ic[index_mode][index_ic1_ic2*psp->ct_size]),psp->ct_size),
		       pop->error_message,
		       pop->error_message);
	  }
	}
      }
    }
  }

  /** - finally, close files and free arrays of files and cls */

  for (index_mode = 0; index_mode < ppt->md_size; index_mode++) {
    if (ppt->ic_size[index_mode] > 1) {
      for (index_ic1_ic2 = 0; index_ic1_ic2 < psp->ic_ic_size[index_mode]; index_ic1_ic2++) {
	if (psp->is_non_zero[index_mode][index_ic1_ic2] == _TRUE_) {
	  fclose(out_md_ic[index_mode][index_ic1_ic2]);
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

  return _SUCCESS_;

} 

/** 
 * This routines writes the output in files for Fourier matter power spectra P(k)'s.
 * 
 * @param pba Input: pointer to background structure (needed for calling spectra_pk_at_z())
 * @param ppt Input : pointer perturbation structure
 * @param psp Input : pointer to spectra structure
 * @param pop Input : pointer to output structure
 */

int output_pk(
	      struct background * pba,
	      struct perturbs * ppt,
	      struct spectra * psp,
	      struct output * pop
	      ) {

  /** Summary: */

  /** - define local variables */

  FILE ** out_ic; /* array of pointers to files with argument 
		     out_ic[index_ic1_ic2] 
		     (will contain P(k)'s for each pair of initial conditions) */
  
  FILE * out;     /* (will contain total P(k) summed eventually over initial conditions) */
  
  double * pk_ic;  /* array with argument 
		      pk_ic[index_ic1_ic2 * psp->ln_k_size + index_k] */

  double * pk_tot; /* array with argument 
		      pk_tot[index_k] */

  int index_mode;
  int index_ic1,index_ic2,index_ic1_ic2;
  int l;
  int index_k;
  int index_z;

  FileName file_name;
  FileName redshift_suffix;
  char first_line[_LINE_LENGTH_MAX_];
    
  index_mode=ppt->index_md_scalars;

  for (index_z = 0; index_z < pop->z_pk_num; index_z++) {

    /** - first, check that requested redshift z_pk is consistent */

    class_test((pop->z_pk[index_z] > psp->z_max_pk),
	       pop->error_message,
	       "P(k,z) computed up to z=%f but requested at z=%f. Must increase z_max_pk in precision file.",psp->z_max_pk,pop->z_pk[index_z]);

    if (pop->z_pk_num == 1) 
      sprintf(redshift_suffix,"");
    else
      sprintf(redshift_suffix,"z%d_",index_z+1);

    /** - second, open only the relevant files, and write a heading in each of them */
    
    sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk.dat");

    class_call(output_open_pk_file(psp,
				   pop,
				   &out,
				   file_name,
				   "",
				   pop->z_pk[index_z]
				   ),
	       pop->error_message,
	       pop->error_message);
   
    class_alloc(pk_tot,
		psp->ln_k_size*sizeof(double),
		pop->error_message);

    if (psp->ic_size[index_mode] > 1) {

      class_alloc(out_ic,
		  psp->ic_ic_size[index_mode]*sizeof(FILE *),
		  pop->error_message);

      class_alloc(pk_ic,
		  psp->ln_k_size*psp->ic_ic_size[index_mode]*sizeof(double),
		  pop->error_message);

      for (index_ic1 = 0; index_ic1 < ppt->ic_size[index_mode]; index_ic1++) {
	  
	for (index_ic2 = index_ic1; index_ic2 < ppt->ic_size[index_mode]; index_ic2++) {
	  
	  if ((ppt->has_ad) && 
	      (index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_ad)) {
	    
	    sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_ad.dat");
	    strcpy(first_line,"for adiabatic (AD) mode ");
	  }

	  if ((ppt->has_bi) && 
	      (index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_bi)) {
	    
	    sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_bi.dat");
	    strcpy(first_line,"for baryon isocurvature (BI) mode ");
	  }
	  
	  if ((ppt->has_cdi) && 
	      (index_ic1 == ppt->index_ic_cdi) && (index_ic2 == ppt->index_ic_cdi)) {
	    
	    sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_cdi.dat");
	    strcpy(first_line,"for CDM isocurvature (CDI) mode ");
	  }
	  
	  if ((ppt->has_nid) && 
	      (index_ic1 == ppt->index_ic_nid) && (index_ic2 == ppt->index_ic_nid)) {
	    
	    sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_nid.dat");
	    strcpy(first_line,"for neutrino density isocurvature (NID) mode ");
	  }
	  
	  if ((ppt->has_niv) && 
	      (index_ic1 == ppt->index_ic_niv) && (index_ic2 == ppt->index_ic_niv)) {
	    
	    sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_niv.dat");
	    strcpy(first_line,"for neutrino velocity isocurvature (NIV) mode ");
	  }
	  
	  if ((ppt->has_ad) && 
	      (ppt->has_bi) && (index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_bi)) {
	    
	    sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_ad_bi.dat");
	    strcpy(first_line,"for cross ADxBI mode ");
	  }
	  
	  if ((ppt->has_ad) && (ppt->has_cdi) && 
	      (index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_cdi)) {
	    
	    sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_ad_cdi.dat");
	    strcpy(first_line,"for cross ADxCDI mode ");
	  }
	  
	  if ((ppt->has_ad) && (ppt->has_nid) && 
	      (index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_nid)) {
	    
	    sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_ad_nid.dat");
	    strcpy(first_line,"for scalar cross ADxNID mode ");
	  }
	  
	  if ((ppt->has_ad) && (ppt->has_niv) && 
	      (index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_niv)) {
	    
	    sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_ad_niv.dat");
	    strcpy(first_line,"for cross ADxNIV mode ");
	  }
	  
	  if ((ppt->has_bi) && (ppt->has_cdi) && 
	      (index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_cdi)) {
	    
	    sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_bi_cdi.dat");
	    strcpy(first_line,"for cross BIxCDI mode ");
	  }
	  
	  if ((ppt->has_bi) && (ppt->has_nid) && 
	      (index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_nid)) {
	    
	    sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_bi_nid.dat");
	    strcpy(first_line,"for cross BIxNID mode ");
	  }
	  
	  if ((ppt->has_bi) && (ppt->has_niv) && 
	      (index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_niv)) {
	    
	    sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_bi_niv.dat");
	    strcpy(first_line,"for cross BIxNIV mode ");
	  }
	  
	  if ((ppt->has_cdi) && (ppt->has_nid) && 
	      (index_ic1 == ppt->index_ic_cdi) && (index_ic2 == ppt->index_ic_nid)) {
	    
	    sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_cdi_nid.dat");
	    strcpy(first_line,"for cross CDIxNID mode ");
	  }
	  
	  if ((ppt->has_cdi) && (ppt->has_niv) && 
	      (index_ic1 == ppt->index_ic_cdi) && (index_ic2 == ppt->index_ic_niv)) {
	    
	    sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_cdi_niv.dat");
	    strcpy(first_line,"for cross CDIxNIV mode ");
	  }
	  
	  if ((ppt->has_nid) && (ppt->has_niv) && 
	      (index_ic1 == ppt->index_ic_nid) && (index_ic2 == ppt->index_ic_niv)) {
	    
	    sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_nid_niv.dat");
	    strcpy(first_line,"for cross NIDxNIV mode ");
	  }

	  index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_mode]);

	  if (psp->is_non_zero[index_mode][index_ic1_ic2] == _TRUE_) {

	    class_call(output_open_pk_file(psp,
					   pop,
					   &(out_ic[index_ic1_ic2]),
					   file_name,
					   first_line,
					   pop->z_pk[index_z]
					   ),
		       pop->error_message,
		       pop->error_message);
	  }
	} 
      }
    }
  
    /** - third, compute P(k) for each k (if several ic's, compute it for each ic and compute also the total); if z_pk = 0, this i9s done by directly reading inside the pre-computed table; if not, this is done by interpolating the table at the correct value of eta. */
  
    /* if z_pk = 0, no interpolation needed */

    if (pop->z_pk[index_z] == 0.) {

      for (index_k=0; index_k<psp->ln_k_size; index_k++) {

	if (psp->ic_size[index_mode] == 1) {
	  pk_tot[index_k] = exp(psp->ln_pk[(psp->ln_eta_size-1) * psp->ln_k_size + index_k]);
	}
	else {
	  pk_tot[index_k] = 0.;
	  for (index_ic1=0; index_ic1 < psp->ic_size[index_mode]; index_ic1++) {
	    index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_mode]);
	    pk_ic[index_ic1_ic2 * psp->ln_k_size + index_k] = exp(psp->ln_pk[((psp->ln_eta_size-1) * psp->ic_ic_size[index_mode] + index_ic1_ic2) * psp->ln_k_size + index_k]);
	    pk_tot[index_k] += pk_ic[index_ic1_ic2 * psp->ln_k_size + index_k];
	  }
	  for (index_ic1=0; index_ic1 < psp->ic_size[index_mode]; index_ic1++) {
	    for (index_ic2 = index_ic1+1; index_ic2 < psp->ic_size[index_mode]; index_ic2++) {
	      pk_ic[index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_mode])* psp->ln_k_size + index_k] = 
		psp->ln_pk[index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_mode])* psp->ln_k_size + index_k]
		*sqrt(pk_ic[index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_mode])* psp->ln_k_size + index_k] *
		      pk_ic[index_symmetric_matrix(index_ic2,index_ic2,psp->ic_size[index_mode])* psp->ln_k_size + index_k]);
	      pk_tot[index_k] += 2.*pk_ic[index_ic1_ic2 * psp->ln_k_size + index_k];
	    }
	  }
	}
      }
    }
      
    /* if 0 <= z_pk <= z_max_pk, interpolation needed, */
    else {

      class_call(spectra_pk_at_z(pba,
				 psp,
				 linear,
				 pop->z_pk[index_z],
				 pk_tot,
				 pk_ic),
		 psp->error_message,
		 pop->error_message);
    }

    /** - fourth, write in files */

    for (index_k=0; index_k<psp->ln_k_size; index_k++) {
	
      class_call(output_one_line_of_pk(out,
				       exp(psp->ln_k[index_k]),
				       pk_tot[index_k]),
		 pop->error_message,
		 pop->error_message);

      if (psp->ic_size[index_mode] > 1) {
	  
	for (index_ic1_ic2 = 0; index_ic1_ic2 < psp->ic_ic_size[index_mode]; index_ic1_ic2++) {

	  if (psp->is_non_zero[index_mode][index_ic1_ic2] == _TRUE_) {

	    class_call(output_one_line_of_pk(out_ic[index_ic1_ic2],
					     exp(psp->ln_k[index_k]),
					     pk_ic[index_ic1_ic2*psp->ln_k_size+index_k]),
		       pop->error_message,
		       pop->error_message);
	  }
	}
      }
    }

    /** - fifth, free memory and close files */

    free(pk_tot);
    fclose(out);

    if (psp->ic_size[index_mode] > 1) {
      for (index_ic1_ic2 = 0; index_ic1_ic2 < psp->ic_ic_size[index_mode]; index_ic1_ic2++) {
	if (psp->is_non_zero[index_mode][index_ic1_ic2] == _TRUE_) {
	  fclose(out_ic[index_ic1_ic2]);
	}
      }
      free(out_ic);
      free(pk_ic);
    }
  
  }
  
  return _SUCCESS_;

}

/**
 * This routine opens one file where some C_l's will be written, and writes 
 * a heading with some general information concerning its content.
 *
 * @param psp        Input : pointer to spectra structure
 * @param pop        Input : pointer to output structure
 * @param clfile     Output: returned pointer to file pointer
 * @param filename   Input : name of the file
 * @param first_line Input : text describing the content (mode, initial condition..)
 * @param lmax       Input : last multipole in the file (the first one is assmued to be 2)
 * @return the error status
 */

int output_open_cl_file(
			struct spectra * psp,
			struct output * pop,
			FILE * * clfile,
			FileName filename,
			char * first_line,
			int lmax
			) {

  class_open(*clfile,filename,"w",pop->error_message);

  fprintf(*clfile,"%s\n",first_line); 
  fprintf(*clfile,"# for l=2 to %d,\n",lmax);
  fprintf(*clfile,"# i.e. number of multipoles equal to %d\n",lmax-1);
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

/**
 * This routine write one line with l and all C_l's for all types (TT, TE...)
 *
 * @param clfile  Input : file pointer
 * @param l       Input : multipole
 * @param cl      Input : C_l's for all types
 * @param ct_size Input : number of types
 * @return the error status
 */

int output_one_line_of_cl(
			  FILE * clfile,
			  double l,
			  double * cl, /* array with argument cl[index_ct] */
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

/**
 * This routine opens one file where some P(k)'s will be written, and writes 
 * a heading with some general information concerning its content.
 *
 * @param psp        Input : pointer to spectra structure
 * @param pop        Input : pointer to output structure
 * @param clfile     Output: returned pointer to file pointer
 * @param filename   Input : name of the file
 * @param first_line Input : text describing the content (initial conditions, ...)
 * @param z          Input : redshift of the output
 * @return the error status
 */

int output_open_pk_file(
			struct spectra * psp,
			struct output * pop,
			FILE * * clfile,
			FileName filename,
			char * first_line,
			double z
			) {

  class_open(*clfile,filename,"w",pop->error_message);

  fprintf(*clfile,"# Matter power spectrum P(k) %sat redshift z=%g\n",first_line,z); 
  fprintf(*clfile,"# for k=%g to %g h/Mpc,\n",exp(psp->ln_k[0]),exp(psp->ln_k[psp->ln_k_size-1]));
  fprintf(*clfile,"# number of wavenumbers equal to %d\n",psp->ln_k_size);
  fprintf(*clfile,"# k (h/Mpc)  P (Mpc/h)^3:\n");

  return _SUCCESS_;
}

/**
 * This routine write one line with k and P(k)
 *
 * @param clfile  Input : file pointer
 * @param one_k   Input : wavenumber
 * @param one_pk  Input : matter power sectrum
 * @return the error status
 */

int output_one_line_of_pk(
			  FILE * clfile,
			  double one_k,
			  double one_pk
			  ) {

  fprintf(clfile,"%e %e\n",one_k,one_pk);
    
  return _SUCCESS_;
    
}
