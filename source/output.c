/** @file output.c Documented output module
 *
 * Julien Lesgourgues, 26.08.2010
 *
 * This module writes the output in files.
 *
 * The following function can be called from other modules or from the main:
 *
 * -# output_init() (must be called after spectra_init())
 * -# output_total_cl_at_l() (can be called even before output_init())
 *
 * No memory needs to be deallocated after that,
 * hence there is no output_free() routine like in other modules.
 */

#include "output.h"

int output_total_cl_at_l(
                         struct spectra * psp,
                         struct lensing * ple,
                         struct output * pop,
                         int l,
                         double * cl
                         ){

  double ** cl_md_ic; /* array with argument
                         cl_md_ic[index_md][index_ic1_ic2*psp->ct_size+index_ct] */

  double ** cl_md;    /* array with argument
                         cl_md[index_md][index_ct] */

  int index_md;

  if (ple->has_lensed_cls == _TRUE_) {
    class_call(lensing_cl_at_l(ple,
                               l,
                               cl),
               ple->error_message,
               pop->error_message);
  }
  else {

    class_alloc(cl_md_ic,
                psp->md_size*sizeof(double *),
                pop->error_message);

    class_alloc(cl_md,
                psp->md_size*sizeof(double *),
                pop->error_message);

    for (index_md = 0; index_md < psp->md_size; index_md++) {

      if (psp->md_size > 1)

        class_alloc(cl_md[index_md],
                    psp->ct_size*sizeof(double),
                    ple->error_message);

      if (psp->ic_size[index_md] > 1)

        class_alloc(cl_md_ic[index_md],
                    psp->ic_ic_size[index_md]*psp->ct_size*sizeof(double),
                    ple->error_message);
    }

    class_call(spectra_cl_at_l(psp,
                               (double)l,
                               cl,
                               cl_md,
                               cl_md_ic),
               psp->error_message,
               pop->error_message);

    for (index_md = 0; index_md < psp->md_size; index_md++) {

      if (psp->md_size > 1)
        free(cl_md[index_md]);

      if (psp->ic_size[index_md] > 1)
        free(cl_md_ic[index_md]);

    }

    free(cl_md_ic);
    free(cl_md);

  }

  return _SUCCESS_;

}

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
                struct thermo * pth,
                struct perturbs * ppt,
                struct primordial * ppm,
                struct transfers * ptr,
                struct spectra * psp,
                struct nonlinear * pnl,
                struct lensing * ple,
                struct output * pop
                ) {

  /** Summary: */

  /** - check that we really want to output at least one file */

  if ((ppt->has_cls == _FALSE_) && (ppt->has_pk_matter == _FALSE_) && (ppt->has_density_transfers == _FALSE_) && (ppt->has_velocity_transfers == _FALSE_) && (pop->write_background == _FALSE_) && (pop->write_thermodynamics == _FALSE_) && (pop->write_primordial == _FALSE_)) {
    if (pop->output_verbose > 0)
      printf("No output files requested. Output module skipped.\n");
    return _SUCCESS_;
  }
  else {
    if (pop->output_verbose > 0)
      printf("Writing output files in %s... \n",pop->root);
  }

  /** - deal with all anisotropy power spectra C_l's */

  if (ppt->has_cls == _TRUE_) {

    class_call(output_cl(pba,ppt,psp,ple,pop),
               pop->error_message,
               pop->error_message);
  }

  /** - deal with all Fourier matter power spectra P(k)'s */

  if (ppt->has_pk_matter == _TRUE_) {

    class_call(output_pk(pba,ppt,psp,pop),
               pop->error_message,
               pop->error_message);

    if (pnl->method != nl_none) {
          class_call(output_pk_nl(pba,ppt,psp,pop),
                     pop->error_message,
                     pop->error_message);
    }
  }

  /** - deal with density and matter power spectra */

  if ((ppt->has_density_transfers == _TRUE_) || (ppt->has_velocity_transfers == _TRUE_)) {

    class_call(output_tk(pba,ppt,psp,pop),
               pop->error_message,
               pop->error_message);
  }

  /** - deal with background quantitites */

  if (pop->write_background == _TRUE_) {

    class_call(output_background(pba,pop),
               pop->error_message,
               pop->error_message);

  }

  /** - deal with thermodynamics quantitites */

  if (pop->write_thermodynamics == _TRUE_) {

    class_call(output_thermodynamics(pba,pth,pop),
               pop->error_message,
               pop->error_message);

  }

  /** - deal with primordial spectra */

  if (pop->write_primordial == _TRUE_) {

    class_call(output_primordial(ppt,ppm,pop),
               pop->error_message,
               pop->error_message);

  }

  return _SUCCESS_;

}

/**
 * This routines writes the output in files for anisotropy power spectra C_l's.
 *
 * @param pba Input: pointer to background structure (needed for T_cmb)
 * @param ppt Input : pointer perturbation structure
 * @param psp Input : pointer to spectra structure
 * @param ple Input : pointer to lensing structure
 * @param pop Input : pointer to output structure
 */

int output_cl(
              struct background * pba,
              struct perturbs * ppt,
              struct spectra * psp,
              struct lensing * ple,
              struct output * pop
              ) {

  /** Summary: */

  /** - define local variables */

  FILE *** out_md_ic; /* array of pointers to files with argument
                         out_md_ic[index_md][index_ic1_ic2]
                         (will contain cl's for each mode and pairs of initial conditions) */

  FILE ** out_md;     /* array of pointers to files with argument
                         out_md[index_md]
                         (will contain cl's for each mode, summed eventually over ic's) */

  FILE * out;         /* (will contain total cl's, summed eventually over modes and ic's) */

  FILE * out_lensed;         /* (will contain total lensed cl's) */

  double ** cl_md_ic; /* array with argument
                         cl_md_ic[index_md][index_ic1_ic2*psp->ct_size+index_ct] */

  double ** cl_md;    /* array with argument
                         cl_md[index_md][index_ct] */

  double * cl_tot;    /* array with argument
                         cl_tot[index_ct] */

  int index_md;
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

  for (index_md = 0; index_md < ppt->md_size; index_md++) {

    class_alloc(out_md_ic[index_md],
                psp->ic_ic_size[index_md]*sizeof(FILE *),
                pop->error_message);

  }

  /** - second, open only the relevant files, and write a heading in each of them */

  sprintf(file_name,"%s%s",pop->root,"cl.dat");

  class_call(output_open_cl_file(psp,
                                 pop,
                                 &out,
                                 file_name,
                                 "total [l(l+1)/2pi] C_l's",
                                 psp->l_max_tot
                                 ),
             pop->error_message,
             pop->error_message);

  class_alloc(cl_tot,
              psp->ct_size*sizeof(double),
              pop->error_message);


  if (ple->has_lensed_cls == _TRUE_) {

    sprintf(file_name,"%s%s",pop->root,"cl_lensed.dat");

    class_call(output_open_cl_file(psp,
                                   pop,
                                   &out_lensed,
                                   file_name,
                                   "total lensed [l(l+1)/2pi] C_l's",
                                   ple->l_lensed_max
                                   ),
               pop->error_message,
               pop->error_message);
  }

  if (ppt->md_size > 1) {

    for (index_md = 0; index_md < ppt->md_size; index_md++) {

      if (_scalars_) {

        sprintf(file_name,"%s%s",pop->root,"cls.dat");
        strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar mode");

      }

      if (_tensors_) {

        sprintf(file_name,"%s%s",pop->root,"clt.dat");
        strcpy(first_line,"[l(l+1)/2pi] C_l's for tensor mode");

      }

      class_call(output_open_cl_file(psp,
                                     pop,
                                     &(out_md[index_md]),
                                     file_name,
                                     first_line,
                                     psp->l_max[index_md]
                                     ),
                 pop->error_message,
                 pop->error_message);

      class_alloc(cl_md[index_md],
                  psp->ct_size*sizeof(double),
                  pop->error_message);

    }
  }

  for (index_md = 0; index_md < ppt->md_size; index_md++) {

    if (ppt->ic_size[index_md] > 1) {

      for (index_ic1 = 0; index_ic1 < ppt->ic_size[index_md]; index_ic1++) {

        for (index_ic2 = index_ic1; index_ic2 < ppt->ic_size[index_md]; index_ic2++) {

          if (_scalars_) {

            if ((ppt->has_ad == _TRUE_) &&
                (index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_ad)) {

              sprintf(file_name,"%s%s",pop->root,"cls_ad.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar adiabatic (AD) mode");
            }

            if ((ppt->has_bi == _TRUE_) &&
                (index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_bi)) {

              sprintf(file_name,"%s%s",pop->root,"cls_bi.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar baryon isocurvature (BI) mode");
            }

            if ((ppt->has_cdi == _TRUE_) &&
                (index_ic1 == ppt->index_ic_cdi) && (index_ic2 == ppt->index_ic_cdi)) {

              sprintf(file_name,"%s%s",pop->root,"cls_cdi.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar CDM isocurvature (CDI) mode");
            }

            if ((ppt->has_nid == _TRUE_) &&
                (index_ic1 == ppt->index_ic_nid) && (index_ic2 == ppt->index_ic_nid)) {

              sprintf(file_name,"%s%s",pop->root,"cls_nid.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar neutrino density isocurvature (NID) mode");
            }

            if ((ppt->has_niv == _TRUE_) &&
                (index_ic1 == ppt->index_ic_niv) && (index_ic2 == ppt->index_ic_niv)) {

              sprintf(file_name,"%s%s",pop->root,"cls_niv.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar neutrino velocity isocurvature (NIV) mode");
            }

            if ((ppt->has_ad == _TRUE_) &&
                (ppt->has_bi == _TRUE_) && (index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_bi)) {

              sprintf(file_name,"%s%s",pop->root,"cls_ad_bi.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar cross ADxBI mode");
            }

            if ((ppt->has_ad == _TRUE_) && (ppt->has_cdi == _TRUE_) &&
                (index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_cdi)) {

              sprintf(file_name,"%s%s",pop->root,"cls_ad_cdi.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar cross ADxCDI mode");
            }

            if ((ppt->has_ad == _TRUE_) && (ppt->has_nid == _TRUE_) &&
                (index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_nid)) {

              sprintf(file_name,"%s%s",pop->root,"cls_ad_nid.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar cross ADxNID mode");
            }

            if ((ppt->has_ad == _TRUE_) && (ppt->has_niv == _TRUE_) &&
                (index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_niv)) {

              sprintf(file_name,"%s%s",pop->root,"cls_ad_niv.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar cross ADxNIV mode");
            }

            if ((ppt->has_bi == _TRUE_) && (ppt->has_cdi == _TRUE_) &&
                (index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_cdi)) {

              sprintf(file_name,"%s%s",pop->root,"cls_bi_cdi.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar cross BIxCDI mode");
            }

            if ((ppt->has_bi == _TRUE_) && (ppt->has_nid == _TRUE_) &&
                (index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_nid)) {

              sprintf(file_name,"%s%s",pop->root,"cls_bi_nid.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar cross BIxNID mode");
            }

            if ((ppt->has_bi == _TRUE_) && (ppt->has_niv == _TRUE_) &&
                (index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_niv)) {

              sprintf(file_name,"%s%s",pop->root,"cls_bi_niv.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar cross BIxNIV mode");
            }

            if ((ppt->has_cdi == _TRUE_) && (ppt->has_nid == _TRUE_) &&
                (index_ic1 == ppt->index_ic_cdi) && (index_ic2 == ppt->index_ic_nid)) {

              sprintf(file_name,"%s%s",pop->root,"cls_cdi_nid.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar cross CDIxNID mode");
            }

            if ((ppt->has_cdi == _TRUE_) && (ppt->has_niv == _TRUE_) &&
                (index_ic1 == ppt->index_ic_cdi) && (index_ic2 == ppt->index_ic_niv)) {

              sprintf(file_name,"%s%s",pop->root,"cls_cdi_niv.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar cross CDIxNIV mode");
            }

            if ((ppt->has_nid == _TRUE_) && (ppt->has_niv == _TRUE_) &&
                (index_ic1 == ppt->index_ic_nid) && (index_ic2 == ppt->index_ic_niv)) {

              sprintf(file_name,"%s%s",pop->root,"cls_nid_niv.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar cross NIDxNIV mode");
            }

          }

          if (_tensors_) {

            class_test(0==0,
                       pop->error_message,
                       "Seems that we have mixed initial conditions for tensors? Should not happen!\n");

          }

          index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md]);

          if (psp->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {

            class_call(output_open_cl_file(psp,
                                           pop,
                                           &(out_md_ic[index_md][index_ic1_ic2]),
                                           file_name,
                                           first_line,
                                           psp->l_max[index_md]
                                           ),
                       pop->error_message,
                       pop->error_message);

          }
        }
      }

      class_alloc(cl_md_ic[index_md],
                  psp->ic_ic_size[index_md]*psp->ct_size*sizeof(double),
                  pop->error_message);
    }
  }

  /** - third, perform loop over l. For each multipole, get all C_l's
      by calling spectra_cl_at_l() and distribute the results to
      relevant files */

  for (l = 2; l <= psp->l_max_tot; l++) {

    class_call(spectra_cl_at_l(psp,(double)l,cl_tot,cl_md,cl_md_ic),
               psp->error_message,
               pop->error_message);

    class_call(output_one_line_of_cl(pba,psp,pop,out,(double)l,cl_tot,psp->ct_size),
               pop->error_message,
               pop->error_message);

    if ((ple->has_lensed_cls == _TRUE_) && (l<=ple->l_lensed_max)) {

      class_call(lensing_cl_at_l(ple,
                                 (double)l,
                                 cl_tot),
                 ple->error_message,
                 pop->error_message);

      class_call(output_one_line_of_cl(pba,psp,pop,out_lensed,l,cl_tot,psp->ct_size),
                 pop->error_message,
                 pop->error_message);
    }

    if (ppt->md_size > 1) {
      for (index_md = 0; index_md < ppt->md_size; index_md++) {
        if (l <= psp->l_max[index_md]) {

          class_call(output_one_line_of_cl(pba,psp,pop,out_md[index_md],l,cl_md[index_md],psp->ct_size),
                     pop->error_message,
                     pop->error_message);
        }
      }
    }

    for (index_md = 0; index_md < ppt->md_size; index_md++) {
      if ((ppt->ic_size[index_md] > 1) && (l <= psp->l_max[index_md])) {
        for (index_ic1_ic2 = 0; index_ic1_ic2 < psp->ic_ic_size[index_md]; index_ic1_ic2++) {
          if (psp->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {

            class_call(output_one_line_of_cl(pba,psp,pop,out_md_ic[index_md][index_ic1_ic2],l,&(cl_md_ic[index_md][index_ic1_ic2*psp->ct_size]),psp->ct_size),
                       pop->error_message,
                       pop->error_message);
          }
        }
      }
    }
  }

  /** - finally, close files and free arrays of files and cls */

  for (index_md = 0; index_md < ppt->md_size; index_md++) {
    if (ppt->ic_size[index_md] > 1) {
      for (index_ic1_ic2 = 0; index_ic1_ic2 < psp->ic_ic_size[index_md]; index_ic1_ic2++) {
        if (psp->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {
          fclose(out_md_ic[index_md][index_ic1_ic2]);
        }
      }
      free(cl_md_ic[index_md]);
    }
  }
  if (ppt->md_size > 1) {
    for (index_md = 0; index_md < ppt->md_size; index_md++) {
      fclose(out_md[index_md]);
      free(cl_md[index_md]);
    }
  }
  fclose(out);
  if (ple->has_lensed_cls == _TRUE_) {
    fclose(out_lensed);
  }
  free(cl_tot);
  for (index_md = 0; index_md < ppt->md_size; index_md++) {
    free(out_md_ic[index_md]);
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

  FILE ** out_ic=NULL; /* array of pointers to files with argument
                          out_ic[index_ic1_ic2]
                          (will contain P(k)'s for each pair of initial conditions) */

  FILE * out;     /* (will contain total P(k) summed eventually over initial conditions) */

  double * pk_ic=NULL;  /* array with argument
                           pk_ic[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2] */

  double * pk_tot; /* array with argument
                      pk_tot[index_k] */

  int index_md;
  int index_ic1,index_ic2;
  int index_ic1_ic2=0;
  int index_k;
  int index_z;

  FileName file_name;
  FileName redshift_suffix;
  char first_line[_LINE_LENGTH_MAX_];

  index_md=ppt->index_md_scalars;

  for (index_z = 0; index_z < pop->z_pk_num; index_z++) {

    /** - first, check that requested redshift z_pk is consistent */

    class_test((pop->z_pk[index_z] > psp->z_max_pk),
               pop->error_message,
               "P(k,z) computed up to z=%f but requested at z=%f. Must increase z_max_pk in precision file.",psp->z_max_pk,pop->z_pk[index_z]);

    if (pop->z_pk_num == 1)
      redshift_suffix[0]='\0';
    else
      sprintf(redshift_suffix,"z%d_",index_z+1);

    /** - second, open only the relevant files, and write a heading in each of them */

    sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk.dat");

    class_call(output_open_pk_file(pba,
                                   psp,
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

    if (psp->ic_size[index_md] > 1) {

      class_alloc(out_ic,
                  psp->ic_ic_size[index_md]*sizeof(FILE *),
                  pop->error_message);

      class_alloc(pk_ic,
                  psp->ln_k_size*psp->ic_ic_size[index_md]*sizeof(double),
                  pop->error_message);

      for (index_ic1 = 0; index_ic1 < ppt->ic_size[index_md]; index_ic1++) {

        for (index_ic2 = index_ic1; index_ic2 < ppt->ic_size[index_md]; index_ic2++) {

          if ((ppt->has_ad == _TRUE_) &&
              (index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_ad)) {

            sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_ad.dat");
            strcpy(first_line,"for adiabatic (AD) mode ");
          }

          if ((ppt->has_bi == _TRUE_) &&
              (index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_bi)) {

            sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_bi.dat");
            strcpy(first_line,"for baryon isocurvature (BI) mode ");
          }

          if ((ppt->has_cdi == _TRUE_) &&
              (index_ic1 == ppt->index_ic_cdi) && (index_ic2 == ppt->index_ic_cdi)) {

            sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_cdi.dat");
            strcpy(first_line,"for CDM isocurvature (CDI) mode ");
          }

          if ((ppt->has_nid == _TRUE_) &&
              (index_ic1 == ppt->index_ic_nid) && (index_ic2 == ppt->index_ic_nid)) {

            sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_nid.dat");
            strcpy(first_line,"for neutrino density isocurvature (NID) mode ");
          }

          if ((ppt->has_niv == _TRUE_) &&
              (index_ic1 == ppt->index_ic_niv) && (index_ic2 == ppt->index_ic_niv)) {

            sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_niv.dat");
            strcpy(first_line,"for neutrino velocity isocurvature (NIV) mode ");
          }

          if ((ppt->has_ad == _TRUE_) &&
              (ppt->has_bi == _TRUE_) && (index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_bi)) {

            sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_ad_bi.dat");
            strcpy(first_line,"for cross ADxBI mode ");
          }

          if ((ppt->has_ad == _TRUE_) && (ppt->has_cdi == _TRUE_) &&
              (index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_cdi)) {

            sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_ad_cdi.dat");
            strcpy(first_line,"for cross ADxCDI mode ");
          }

          if ((ppt->has_ad == _TRUE_) && (ppt->has_nid == _TRUE_) &&
              (index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_nid)) {

            sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_ad_nid.dat");
            strcpy(first_line,"for scalar cross ADxNID mode ");
          }

          if ((ppt->has_ad == _TRUE_) && (ppt->has_niv == _TRUE_) &&
              (index_ic1 == ppt->index_ic_ad) && (index_ic2 == ppt->index_ic_niv)) {

            sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_ad_niv.dat");
            strcpy(first_line,"for cross ADxNIV mode ");
          }

          if ((ppt->has_bi == _TRUE_) && (ppt->has_cdi == _TRUE_) &&
              (index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_cdi)) {

            sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_bi_cdi.dat");
            strcpy(first_line,"for cross BIxCDI mode ");
          }

          if ((ppt->has_bi == _TRUE_) && (ppt->has_nid == _TRUE_) &&
              (index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_nid)) {

            sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_bi_nid.dat");
            strcpy(first_line,"for cross BIxNID mode ");
          }

          if ((ppt->has_bi == _TRUE_) && (ppt->has_niv == _TRUE_) &&
              (index_ic1 == ppt->index_ic_bi) && (index_ic2 == ppt->index_ic_niv)) {

            sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_bi_niv.dat");
            strcpy(first_line,"for cross BIxNIV mode ");
          }

          if ((ppt->has_cdi == _TRUE_) && (ppt->has_nid == _TRUE_) &&
              (index_ic1 == ppt->index_ic_cdi) && (index_ic2 == ppt->index_ic_nid)) {

            sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_cdi_nid.dat");
            strcpy(first_line,"for cross CDIxNID mode ");
          }

          if ((ppt->has_cdi == _TRUE_) && (ppt->has_niv == _TRUE_) &&
              (index_ic1 == ppt->index_ic_cdi) && (index_ic2 == ppt->index_ic_niv)) {

            sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_cdi_niv.dat");
            strcpy(first_line,"for cross CDIxNIV mode ");
          }

          if ((ppt->has_nid == _TRUE_) && (ppt->has_niv == _TRUE_) &&
              (index_ic1 == ppt->index_ic_nid) && (index_ic2 == ppt->index_ic_niv)) {

            sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_nid_niv.dat");
            strcpy(first_line,"for cross NIDxNIV mode ");
          }

          index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md]);

          if (psp->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {

            class_call(output_open_pk_file(pba,
                                           psp,
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

    /** - third, compute P(k) for each k (if several ic's, compute it for each ic and compute also the total); if z_pk = 0, this is done by directly reading inside the pre-computed table; if not, this is done by interpolating the table at the correct value of tau. */

    /* if z_pk = 0, no interpolation needed */

    if (pop->z_pk[index_z] == 0.) {

      for (index_k=0; index_k<psp->ln_k_size; index_k++) {

        if (psp->ic_size[index_md] == 1) {
          pk_tot[index_k] = exp(psp->ln_pk[(psp->ln_tau_size-1) * psp->ln_k_size + index_k]);
        }
        else {
          pk_tot[index_k] = 0.;
          for (index_ic1=0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
            index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_md]);
            pk_ic[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2] = exp(psp->ln_pk[((psp->ln_tau_size-1) * psp->ln_k_size + index_k) * psp->ic_ic_size[index_md] + index_ic1_ic2]);
            pk_tot[index_k] += pk_ic[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2];
          }
          for (index_ic1=0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
            for (index_ic2 = index_ic1+1; index_ic2 < psp->ic_size[index_md]; index_ic2++) {
              pk_ic[index_k * psp->ic_ic_size[index_md] + index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md])] =
                psp->ln_pk[index_k * psp->ic_ic_size[index_md] + index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md])]
                *sqrt(pk_ic[index_k * psp->ic_ic_size[index_md] + index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_md])] *
                      pk_ic[index_k * psp->ic_ic_size[index_md] + index_symmetric_matrix(index_ic2,index_ic2,psp->ic_size[index_md])]);
              pk_tot[index_k] += 2.*pk_ic[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2];
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
                                       exp(psp->ln_k[index_k])/pba->h,
                                       pk_tot[index_k]*pow(pba->h,3)),
                 pop->error_message,
                 pop->error_message);

      if (psp->ic_size[index_md] > 1) {

        for (index_ic1_ic2 = 0; index_ic1_ic2 < psp->ic_ic_size[index_md]; index_ic1_ic2++) {

          if (psp->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {

            class_call(output_one_line_of_pk(out_ic[index_ic1_ic2],
                                             exp(psp->ln_k[index_k])/pba->h,
                                             pk_ic[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2]*pow(pba->h,3)),
                       pop->error_message,
                       pop->error_message);
          }
        }
      }
    }

    /** - fifth, free memory and close files */

    free(pk_tot);
    fclose(out);

    if (psp->ic_size[index_md] > 1) {
      for (index_ic1_ic2 = 0; index_ic1_ic2 < psp->ic_ic_size[index_md]; index_ic1_ic2++) {
        if (psp->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {
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
 * This routines writes the output in files for Fourier non-linear matter power spectra P(k)'s.
 *
 * @param pba Input: pointer to background structure (needed for calling spectra_pk_at_z())
 * @param ppt Input : pointer perturbation structure
 * @param psp Input : pointer to spectra structure
 * @param pop Input : pointer to output structure
 */

int output_pk_nl(
                 struct background * pba,
                 struct perturbs * ppt,
                 struct spectra * psp,
                 struct output * pop
                 ) {

  /** Summary: */

  /** - define local variables */

  FILE * out;     /* (will contain total P(k) summed eventually over initial conditions) */

  double * pk_tot; /* array with argument pk_tot[index_k] */

  int index_k;
  int index_z;

  FileName file_name;
  FileName redshift_suffix;

  for (index_z = 0; index_z < pop->z_pk_num; index_z++) {

    /** - first, check that requested redshift z_pk is consistent */

    class_test((pop->z_pk[index_z] > psp->z_max_pk),
               pop->error_message,
               "P(k,z) computed up to z=%f but requested at z=%f. Must increase z_max_pk in precision file.",psp->z_max_pk,pop->z_pk[index_z]);

    if (pop->z_pk_num == 1)
      redshift_suffix[0]='\0';
    else
      sprintf(redshift_suffix,"z%d_",index_z+1);

    /** - second, open only the relevant files, and write a heading in each of them */

    sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"pk_nl.dat");

    class_call(output_open_pk_file(pba,
                                   psp,
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

    /** - third, compute P(k) for each k (if several ic's, compute it for each ic and compute also the total); if z_pk = 0, this is done by directly reading inside the pre-computed table; if not, this is done by interpolating the table at the correct value of tau. */

    /* if z_pk = 0, no interpolation needed */

    if (pop->z_pk[index_z] == 0.) {

      for (index_k=0; index_k<psp->ln_k_size; index_k++) {

        pk_tot[index_k] = exp(psp->ln_pk_nl[(psp->ln_tau_size-1) * psp->ln_k_size + index_k]);

      }
    }

    /* if 0 <= z_pk <= z_max_pk, interpolation needed, */
    else {

      class_call(spectra_pk_nl_at_z(pba,
                                    psp,
                                    linear,
                                    pop->z_pk[index_z],
                                    pk_tot),
                 psp->error_message,
                 pop->error_message);
    }

    /** - fourth, write in files */

    for (index_k=0; index_k<psp->ln_k_size; index_k++) {

      class_call(output_one_line_of_pk(out,
                                       exp(psp->ln_k[index_k])/pba->h,
                                       pk_tot[index_k]*pow(pba->h,3)),
                 pop->error_message,
                 pop->error_message);

    }

    /** - fifth, free memory and close files */

    fclose(out);
    free(pk_tot);

  }

  return _SUCCESS_;

}

/**
 * This routines writes the output in files for matter transfer functions T_i(k)'s.
 *
 * @param pba Input: pointer to background structure (needed for calling spectra_pk_at_z())
 * @param ppt Input : pointer perturbation structure
 * @param psp Input : pointer to spectra structure
 * @param pop Input : pointer to output structure
 */

int output_tk(
              struct background * pba,
              struct perturbs * ppt,
              struct spectra * psp,
              struct output * pop
              ) {

  /** Summary: */

  /** - define local variables */

  FILE ** out_ic; /* array of pointers to files with argument
                     out_ic[index_ic]
                     (will contain T_i(k)'s for each initial conditions) */

  double * tk;  /* array with argument
                   pk_ic[(index_k * psp->ic_size[index_md] + index_ic)*psp->tr_size+index_tr] */

  double * tk_cmbfast = NULL; /* array with argument tk_cmbfast[index_tr] */


  int index_md;
  int index_ic;
  int index_k;
  int index_z;
  int index_tr;

  FileName file_name;
  FileName redshift_suffix;
  char first_line[_LINE_LENGTH_MAX_];

  index_md=ppt->index_md_scalars;

  if (pop->output_format == camb_format) {

    class_test(pba->N_ncdm>1,
               pop->error_message,
               "you wish to output the transfer functions in CMBFAST/CAMB format but you have more than one non-cold dark matter (ncdm) species. The two are not compatible (since CMBFAST/CAMB only have one ncdm species): switch to CLASS output format or keep only on ncdm species");

    class_test(ppt->has_velocity_transfers == _TRUE_,
               pop->error_message,
               "you wish to output the transfer functions in CMBFAST/CAMB format, but you requested velocity transfer functions. The two are not compatible (since CMBFAST/CAMB do not compute velocity transfer functions): switch to CLASS output format, or ask only for density transfer function");
  }

  for (index_z = 0; index_z < pop->z_pk_num; index_z++) {

    /** - first, check that requested redshift z_pk is consistent */

    class_test((pop->z_pk[index_z] > psp->z_max_pk),
               pop->error_message,
               "T_i(k,z) computed up to z=%f but requested at z=%f. Must increase z_max_pk in precision file.",psp->z_max_pk,pop->z_pk[index_z]);

    if (pop->z_pk_num == 1)
      redshift_suffix[0]='\0';
    else
      sprintf(redshift_suffix,"z%d_",index_z+1);

    /** - second, open only the relevant files, and write a heading in each of them */

    class_alloc(out_ic,
                psp->ic_size[index_md]*sizeof(FILE *),
                pop->error_message);

    class_alloc(tk,
                psp->ln_k_size*psp->ic_size[index_md]*psp->tr_size*sizeof(double),
                pop->error_message);

    for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++) {

      if ((ppt->has_ad == _TRUE_) && (index_ic == ppt->index_ic_ad)) {

        if (ppt->ic_size[index_md] == 1)
          sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"tk.dat");
        else
          sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"tk_ad.dat");
        strcpy(first_line,"for adiabatic (AD) mode (normalized to initial curvature=1) ");
      }

      if ((ppt->has_bi == _TRUE_) && (index_ic == ppt->index_ic_bi)) {

        sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"tk_bi.dat");
        strcpy(first_line,"for baryon isocurvature (BI) mode (normalized to initial entropy=1)");
      }

      if ((ppt->has_cdi == _TRUE_) && (index_ic == ppt->index_ic_cdi)) {

        sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"tk_cdi.dat");
        strcpy(first_line,"for CDM isocurvature (CDI) mode (normalized to initial entropy=1)");
      }

      if ((ppt->has_nid == _TRUE_) && (index_ic == ppt->index_ic_nid)) {

        sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"tk_nid.dat");
        strcpy(first_line,"for neutrino density isocurvature (NID) mode (normalized to initial entropy=1)");
      }

      if ((ppt->has_niv == _TRUE_) && (index_ic == ppt->index_ic_niv)) {

        sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"tk_niv.dat");
        strcpy(first_line,"for neutrino velocity isocurvature (NIV) mode (normalized to initial entropy=1)");
      }

      class_call(output_open_tk_file(pba,
                                     ppt,
                                     psp,
                                     pop,
                                     &(out_ic[index_ic]),
                                     file_name,
                                     first_line,
                                     pop->z_pk[index_z]
                                     ),
                 pop->error_message,
                 pop->error_message);

    }

    /** - third, compute T_i(k) for each k (if several ic's, compute it for each ic; if z_pk = 0, this is done by directly reading inside the pre-computed table; if not, this is done by interpolating the table at the correct value of tau. */

    /* if z_pk = 0, no interpolation needed */

    if (pop->z_pk[index_z] == 0.) {

      for (index_k=0; index_k<psp->ln_k_size; index_k++) {
        for (index_tr=0; index_tr<psp->tr_size; index_tr++) {
          for (index_ic=0; index_ic<psp->ic_size[index_md]; index_ic++) {
            tk[(index_k * psp->ic_size[index_md] + index_ic) * psp->tr_size + index_tr] = psp->matter_transfer[(((psp->ln_tau_size-1)*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + index_tr];
          }
        }
      }
    }

    /* if 0 <= z_pk <= z_max_pk, interpolation needed, */
    else {

      class_call(spectra_tk_at_z(pba,
                                 psp,
                                 pop->z_pk[index_z],
                                 tk),
                 psp->error_message,
                 pop->error_message);
    }

    /** - fourth, write in files */

    if (pop->output_format == camb_format)
      class_alloc(tk_cmbfast,
                  6*sizeof(double),
                  pop->error_message);

    for (index_k=0; index_k<psp->ln_k_size; index_k++) {
      for (index_ic = 0; index_ic < psp->ic_size[index_md]; index_ic++) {

        if (pop->output_format == class_format) {

          class_call(output_one_line_of_tk(out_ic[index_ic],
                                           exp(psp->ln_k[index_k])/pba->h,
                                           &(tk[(index_k * psp->ic_size[index_md] + index_ic) * psp->tr_size]),
                                           psp->tr_size),
                     pop->error_message,
                     pop->error_message);

        }
        else if (pop->output_format == camb_format) {

          /* rescale and reorder the matter transfer functions following the CMBFAST/CAMB convention */

          if (pba->has_cdm == _TRUE_)
            tk_cmbfast[0]=-tk[(index_k*psp->ic_size[index_md]+index_ic)*psp->tr_size+psp->index_tr_delta_cdm]/exp(2.*psp->ln_k[index_k]);
          else
            tk_cmbfast[0]= 0.;
          tk_cmbfast[1]=-tk[(index_k*psp->ic_size[index_md]+index_ic)*psp->tr_size+psp->index_tr_delta_b]/exp(2.*psp->ln_k[index_k]);
          tk_cmbfast[2]=-tk[(index_k*psp->ic_size[index_md]+index_ic)*psp->tr_size+psp->index_tr_delta_g]/exp(2.*psp->ln_k[index_k]);
          if (pba->has_ur == _TRUE_)
            tk_cmbfast[3]=-tk[(index_k*psp->ic_size[index_md]+index_ic)*psp->tr_size+psp->index_tr_delta_ur]/exp(2.*psp->ln_k[index_k]);
          else
            tk_cmbfast[3]=0.;
          if (pba->has_ncdm == _TRUE_)
            tk_cmbfast[4]=-tk[(index_k*psp->ic_size[index_md]+index_ic)*psp->tr_size+psp->index_tr_delta_ncdm1]/exp(2.*psp->ln_k[index_k]);
          else
            tk_cmbfast[4]=0.;
          tk_cmbfast[5]=-tk[(index_k*psp->ic_size[index_md]+index_ic)*psp->tr_size+psp->index_tr_delta_tot]/exp(2.*psp->ln_k[index_k]);

          class_call(output_one_line_of_tk(out_ic[index_ic],
                                           exp(psp->ln_k[index_k])/pba->h,
                                           tk_cmbfast,
                                           6),
                     pop->error_message,
                     pop->error_message);

        }
      }
    }

    /** - fifth, free memory and close files */

    if (pop->output_format == camb_format)
      free(tk_cmbfast);

    free(tk);

    for (index_ic = 0; index_ic < psp->ic_size[index_md]; index_ic++) {
      fclose(out_ic[index_ic]);
    }
    free(out_ic);

  }

  return _SUCCESS_;

}

int output_background(
                      struct background * pba,
                      struct output * pop
                      ) {

  FILE * out;
  FileName file_name;
  int index_eta;

  sprintf(file_name,"%s%s",pop->root,"background.dat");

  class_call(output_open_background_file(pba,
                                         pop,
                                         &out,
                                         file_name
                                         ),
             pop->error_message,
             pop->error_message);

  for (index_eta=0; index_eta<pba->bt_size; index_eta++) {

    class_call(output_one_line_of_background(pba,
                                             out,
                                             &(pba->background_table[index_eta*pba->bg_size])
                                             ),
               pop->error_message,
               pop->error_message);

  }

  fclose(out);

  return _SUCCESS_;

}

int output_thermodynamics(
                          struct background * pba,
                          struct thermo * pth,
                          struct output * pop
                      ) {

  FILE * out;
  FileName file_name;
  int index_z;
  double tau;

  sprintf(file_name,"%s%s",pop->root,"thermodynamics.dat");

  class_call(output_open_thermodynamics_file(
                                             pth,
                                             pop,
                                             &out,
                                             file_name
                                             ),
             pop->error_message,
             pop->error_message);

  for (index_z=0; index_z<pth->tt_size; index_z++) {

    class_call(background_tau_of_z(
                                   pba,
                                   pth->z_table[index_z],
                                   &tau
                                   ),
               pop->error_message,
               pop->error_message);

    class_call(output_one_line_of_thermodynamics(
                                                 pth,
                                                 out,
                                                 tau,
                                                 pth->z_table[index_z],
                                                 pth->thermodynamics_table+index_z*pth->th_size
                                                 ),
               pop->error_message,
               pop->error_message);

  }

  fclose(out);

  return _SUCCESS_;

}

int output_primordial(
                      struct perturbs * ppt,
                      struct primordial * ppm,
                      struct output * pop
                      ) {

  FILE * out;
  FileName file_name;
  int index_k;

  sprintf(file_name,"%s%s",pop->root,"primordial_Pk.dat");

  class_call(output_open_primordial_file(ppt,
                                         ppm,
                                         pop,
                                         &out,
                                         file_name
                                         ),
             pop->error_message,
             pop->error_message);

  for (index_k=0; index_k<ppm->lnk_size; index_k++) {

    class_call(output_one_line_of_primordial(ppt,
                                             ppm,
                                             out,
                                             index_k
                                             ),
               pop->error_message,
               pop->error_message);

  }

  fclose(out);

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

  int index_d1,index_d2;

  class_open(*clfile,filename,"w",pop->error_message);

  if (pop->write_header == _TRUE_) {

    if (pop->output_format == class_format) {
      fprintf(*clfile,"# dimensionless %s\n",first_line);
    }
    if (pop->output_format == camb_format) {
      fprintf(*clfile,"# %s (units: [microK]^2)\n",first_line);
    }

    fprintf(*clfile,"# for l=2 to %d, i.e. number of multipoles equal to %d\n",lmax,lmax-1);
    fprintf(*clfile,"#\n");

    if (pop->output_format == class_format) {
      fprintf(*clfile,"# -> if you prefer output in CAMB/HealPix/LensPix units/order, set 'format' to 'camb' in input file\n");
    }

    fprintf(*clfile,"# -> if you don't want to see such a header, set 'headers' to 'no' in input file\n");
    fprintf(*clfile,"#\n");
    fprintf(*clfile,"#  l ");

    if (pop->output_format == class_format) {
      if (psp->has_tt == _TRUE_)
        fprintf(*clfile,"TT               ");
      if (psp->has_ee == _TRUE_)
        fprintf(*clfile,"EE               ");
      if (psp->has_te == _TRUE_)
        fprintf(*clfile,"TE                ");
      if (psp->has_bb == _TRUE_)
        fprintf(*clfile,"BB               ");
      if (psp->has_pp == _TRUE_)
        fprintf(*clfile,"phiphi           ");
      if (psp->has_tp == _TRUE_)
        fprintf(*clfile,"Tphi             ");
      if (psp->has_ep == _TRUE_)
        fprintf(*clfile,"Ephi             ");
      if (psp->has_dd == _TRUE_)
        for (index_d1=0; index_d1<psp->d_size; index_d1++)
          for (index_d2=index_d1; index_d2<=MIN(index_d1+psp->non_diag,psp->d_size-1); index_d2++)
            fprintf(*clfile,"dens[%d]-dens[%d]  ",index_d1+1,index_d2+1);
      if (psp->has_td == _TRUE_)
        for (index_d1=0; index_d1<psp->d_size; index_d1++)
          fprintf(*clfile,"T-dens[%d]        ",index_d1+1);
      if (psp->has_pd == _TRUE_)
        for (index_d1=0; index_d1<psp->d_size; index_d1++)
          fprintf(*clfile,"phi-dens[%d]      ",index_d1+1);
      if (psp->has_ll == _TRUE_)
        for (index_d1=0; index_d1<psp->d_size; index_d1++)
          for (index_d2=index_d1; index_d2<=MIN(index_d1+psp->non_diag,psp->d_size-1); index_d2++)
            fprintf(*clfile,"lens[%d]-lens[%d]  ",index_d1+1,index_d2+1);
      if (psp->has_tl == _TRUE_)
        for (index_d1=0; index_d1<psp->d_size; index_d1++)
          fprintf(*clfile,"T-lens[%d]        ",index_d1+1);
      if (psp->has_dl == _TRUE_)
        for (index_d1=0; index_d1<psp->d_size; index_d1++)
          for (index_d2=index_d1; index_d2<=MIN(index_d1+psp->non_diag,psp->d_size-1); index_d2++)
            fprintf(*clfile,"dens[%d]-lens[%d]  ",index_d1+1,index_d2+1);
      fprintf(*clfile,"\n");
    }

    if (pop->output_format == camb_format) {
      if (psp->has_tt == _TRUE_)
        fprintf(*clfile,"TT               ");
      if (psp->has_ee == _TRUE_)
        fprintf(*clfile,"EE               ");
      if (psp->has_bb == _TRUE_)
        fprintf(*clfile,"BB               ");
      if (psp->has_te == _TRUE_)
        fprintf(*clfile,"TE                ");
      if (psp->has_pp == _TRUE_)
        fprintf(*clfile,"dd               ");
      if (psp->has_tp == _TRUE_)
        fprintf(*clfile,"dT               ");
      if (psp->has_ep == _TRUE_)
        fprintf(*clfile,"dE               ");
      if (psp->has_dd == _TRUE_)
        for (index_d1=0; index_d1<psp->d_size; index_d1++)
          for (index_d2=index_d1; index_d2<=MIN(index_d1+psp->non_diag,psp->d_size-1); index_d2++)
            fprintf(*clfile,"dens[%d]-dens[%d]  ",index_d1+1,index_d2+1);
      if (psp->has_td == _TRUE_)
        for (index_d1=0; index_d1<psp->d_size; index_d1++)
          fprintf(*clfile,"T-dens[%d]        ",index_d1+1);
      if (psp->has_pd == _TRUE_)
        for (index_d1=0; index_d1<psp->d_size; index_d1++)
          fprintf(*clfile,"phi-dens[%d]      ",index_d1+1);
      if (psp->has_ll == _TRUE_)
        for (index_d1=0; index_d1<psp->d_size; index_d1++)
          for (index_d2=index_d1; index_d2<=MIN(index_d1+psp->non_diag,psp->d_size-1); index_d2++)
            fprintf(*clfile,"lens[%d]-lens[%d]  ",index_d1+1,index_d2+1);
      if (psp->has_tl == _TRUE_)
        for (index_d1=0; index_d1<psp->d_size; index_d1++)
          fprintf(*clfile,"T-lens[%d]        ",index_d1+1);
      if (psp->has_dl == _TRUE_)
        for (index_d1=0; index_d1<psp->d_size; index_d1++)
          for (index_d2=index_d1; index_d2<=MIN(index_d1+psp->non_diag,psp->d_size-1); index_d2++)
            fprintf(*clfile,"dens[%d]-lens[%d]  ",index_d1+1,index_d2+1);
      fprintf(*clfile,"\n");

    }
  }

  return _SUCCESS_;

}

/**
 * This routine write one line with l and all C_l's for all types (TT, TE...)
 *
 * @param pba        Input: pointer to background structure (needed for T_cmb)
 * @param psp        Input : pointer to spectra structure
 * @param pop        Input : pointer to output structure
 * @param clfile  Input : file pointer
 * @param l       Input : multipole
 * @param cl      Input : C_l's for all types
 * @param ct_size Input : number of types
 * @return the error status
 */

int output_one_line_of_cl(
                          struct background * pba,
                          struct spectra * psp,
                          struct output * pop,
                          FILE * clfile,
                          double l,
                          double * cl, /* array with argument cl[index_ct] */
                          int ct_size
                          ) {
  int index_ct;
  double factor;

  factor = l*(l+1)/2./_PI_;

  fprintf(clfile,"%4d",(int)l);

  if (pop->output_format == class_format) {

    for (index_ct=0; index_ct < ct_size; index_ct++) {
      fprintf(clfile," %16.10e",factor*cl[index_ct]);
    }
    fprintf(clfile,"\n");
  }

  if (pop->output_format == camb_format) {

    if (psp->has_tt == _TRUE_)
      fprintf(clfile," %16.10e",factor*pow(pba->T_cmb*1.e6,2)*cl[psp->index_ct_tt]);
    if (psp->has_ee == _TRUE_)
      fprintf(clfile," %16.10e",factor*pow(pba->T_cmb*1.e6,2)*cl[psp->index_ct_ee]);
    if (psp->has_bb == _TRUE_)
      fprintf(clfile," %16.10e",factor*pow(pba->T_cmb*1.e6,2)*cl[psp->index_ct_bb]);
    if (psp->has_te == _TRUE_)
      fprintf(clfile," %16.10e",factor*pow(pba->T_cmb*1.e6,2)*cl[psp->index_ct_te]);
    if (psp->has_pp == _TRUE_)
      fprintf(clfile," %16.10e",l*(l+1)*factor*cl[psp->index_ct_pp]);
    if (psp->has_tp == _TRUE_)
      fprintf(clfile," %16.10e",sqrt(l*(l+1))*factor*pba->T_cmb*1.e6*cl[psp->index_ct_tp]);
    if (psp->has_ep == _TRUE_)
      fprintf(clfile," %16.10e",sqrt(l*(l+1))*factor*pba->T_cmb*1.e6*cl[psp->index_ct_ep]);
    if (psp->has_dd == _TRUE_)
      for (index_ct=0; index_ct<(psp->d_size*(psp->d_size+1) - (psp->d_size-psp->non_diag)*(psp->d_size-1-psp->non_diag))/2; index_ct++)
        fprintf(clfile," %16.10e",factor*cl[psp->index_ct_dd+index_ct]);
    if (psp->has_td == _TRUE_)
      for (index_ct=0; index_ct<psp->d_size; index_ct++)
        fprintf(clfile," %16.10e",factor*cl[psp->index_ct_td+index_ct]);
    if (psp->has_pd == _TRUE_)
      for (index_ct=0; index_ct<psp->d_size; index_ct++)
        fprintf(clfile," %16.10e",factor*cl[psp->index_ct_pd+index_ct]);
    if (psp->has_ll == _TRUE_)
      for (index_ct=0; index_ct<(psp->d_size*(psp->d_size+1) - (psp->d_size-psp->non_diag)*(psp->d_size-1-psp->non_diag))/2; index_ct++)
        fprintf(clfile," %16.10e",factor*cl[psp->index_ct_ll+index_ct]);
    if (psp->has_tl == _TRUE_)
      for (index_ct=0; index_ct<psp->d_size; index_ct++)
        fprintf(clfile," %16.10e",factor*cl[psp->index_ct_tl+index_ct]);
    if (psp->has_dl == _TRUE_)
      for (index_ct=0; index_ct<(psp->d_size*(psp->d_size+1) - (psp->d_size-psp->non_diag)*(psp->d_size-1-psp->non_diag))/2; index_ct++)
        fprintf(clfile," %16.10e",factor*cl[psp->index_ct_dl+index_ct]);
    fprintf(clfile,"\n");
  }

  return _SUCCESS_;

}

/**
 * This routine opens one file where some P(k)'s will be written, and writes
 * a heading with some general information concerning its content.
 *
 * @param pba        Input: pointer to background structure (needed for h)
 * @param psp        Input : pointer to spectra structure
 * @param pop        Input : pointer to output structure
 * @param tkfile     Output: returned pointer to file pointer
 * @param filename   Input : name of the file
 * @param first_line Input : text describing the content (initial conditions, ...)
 * @param z          Input : redshift of the output
 * @return the error status
 */

int output_open_pk_file(
                        struct background * pba,
                        struct spectra * psp,
                        struct output * pop,
                        FILE * * pkfile,
                        FileName filename,
                        char * first_line,
                        double z
                        ) {

  class_open(*pkfile,filename,"w",pop->error_message);

  if (pop->write_header == _TRUE_) {
    fprintf(*pkfile,"# Matter power spectrum P(k) %sat redshift z=%g\n",first_line,z);
    fprintf(*pkfile,"# for k=%g to %g h/Mpc,\n",
            exp(psp->ln_k[0])/pba->h,
            exp(psp->ln_k[psp->ln_k_size-1])/pba->h);
    fprintf(*pkfile,"# number of wavenumbers equal to %d\n",psp->ln_k_size);
    fprintf(*pkfile,"# k (h/Mpc)  P (Mpc/h)^3:\n");
  }

  return _SUCCESS_;
}

/**
 * This routine writes one line with k and P(k)
 *
 * @param tkfile  Input : file pointer
 * @param one_k   Input : wavenumber
 * @param one_pk  Input : matter power sectrum
 * @return the error status
 */

int output_one_line_of_pk(
                          FILE * pkfile,
                          double one_k,
                          double one_pk
                          ) {

  fprintf(pkfile,"%e %16.10e\n",one_k,one_pk);

  return _SUCCESS_;

}

/**
 * This routine opens one file where some T_i(k)'s will be written, and writes
 * a heading with some general information concerning its content.
 *
 * @param psp        Input : pointer to spectra structure
 * @param pop        Input : pointer to output structure
 * @param tkfile     Output: returned pointer to file pointer
 * @param filename   Input : name of the file
 * @param first_line Input : text describing the content (initial conditions, ...)
 * @param z          Input : redshift of the output
 * @return the error status
 */

int output_open_tk_file(
                        struct background * pba,
                        struct perturbs * ppt,
                        struct spectra * psp,
                        struct output * pop,
                        FILE * * tkfile,
                        FileName filename,
                        char * first_line,
                        double z
                        ) {

  int n_ncdm;

  class_open(*tkfile,filename,"w",pop->error_message);

  if (pop->write_header == _TRUE_) {

    if (pop->output_format == class_format) {

      fprintf(*tkfile,"# Transfer functions T_i(k) %sat redshift z=%g\n",first_line,z);
      fprintf(*tkfile,"# for k=%g to %g h/Mpc,\n",exp(psp->ln_k[0])/pba->h,exp(psp->ln_k[psp->ln_k_size-1])/pba->h);
      fprintf(*tkfile,"# number of wavenumbers equal to %d\n",psp->ln_k_size);
      if (ppt->has_density_transfers == _TRUE_) {
        fprintf(*tkfile,"# d_i   stands for (delta rho_i/rho_i)(k,z) with above normalization \n");
        fprintf(*tkfile,"# d_tot stands for (delta rho_tot/rho_tot)(k,z) with rho_Lambda NOT included in rho_tot\n");
        fprintf(*tkfile,"# (note that this differs from the transfer function output from CAMB/CMBFAST, which gives the same\n");
        fprintf(*tkfile,"#  quantities divided by -k^2 with k in Mpc^-1; use format=camb to match CAMB)\n");
      }
      if (ppt->has_velocity_transfers == _TRUE_) {
        fprintf(*tkfile,"# t_i   stands for theta_i(k,z) with above normalization \n");
        fprintf(*tkfile,"# t_tot stands for (sum_i [rho_i+p_i] theta_i)/(sum_i [rho_i+p_i]))(k,z)\n");
      }
      fprintf(*tkfile,"#\n");
      fprintf(*tkfile,"# k (h/Mpc)       ");
      if (ppt->has_density_transfers == _TRUE_) {
        fprintf(*tkfile,"d_g                ");
        fprintf(*tkfile,"d_b                ");
        if (pba->has_cdm == _TRUE_)
          fprintf(*tkfile,"d_cdm              ");
        if (pba->has_fld == _TRUE_)
          fprintf(*tkfile,"d_de               ");
        if (pba->has_ur == _TRUE_)
          fprintf(*tkfile,"d_ur               ");
        if (pba->has_ncdm == _TRUE_) {
          for (n_ncdm=0; n_ncdm < pba->N_ncdm; n_ncdm++) {
            fprintf(*tkfile,"d_ncdm[%d]          ",n_ncdm);
          }
        }
        fprintf(*tkfile,"d_tot              ");
      }
      if (ppt->has_velocity_transfers == _TRUE_) {
        fprintf(*tkfile,"t_g                ");
        fprintf(*tkfile,"t_b                ");
        if ((pba->has_cdm == _TRUE_) && (ppt->gauge != synchronous))
          fprintf(*tkfile,"t_cdm              ");
        if (pba->has_fld == _TRUE_)
          fprintf(*tkfile,"t_de               ");
        if (pba->has_ur == _TRUE_)
          fprintf(*tkfile,"t_ur               ");
        if (pba->has_ncdm == _TRUE_) {
          for (n_ncdm=0; n_ncdm < pba->N_ncdm; n_ncdm++) {
            fprintf(*tkfile,"t_ncdm[%d]          ",n_ncdm);
          }
        }
        fprintf(*tkfile,"t_tot              ");
      }
      fprintf(*tkfile,"\n");
    }

    else if (pop->output_format == camb_format) {

      fprintf(*tkfile,"# Rescaled matter transfer functions [-T_i(k)/k^2] %sat redshift z=%g\n",first_line,z);
      fprintf(*tkfile,"# for k=%g to %g h/Mpc,\n",exp(psp->ln_k[0])/pba->h,exp(psp->ln_k[psp->ln_k_size-1])/pba->h);
      fprintf(*tkfile,"# number of wavenumbers equal to %d\n",psp->ln_k_size);
      fprintf(*tkfile,"# T_i   stands for (delta rho_i/rho_i)(k,z) with above normalization \n");
      fprintf(*tkfile,"# The rescaling factor [-1/k^2] with k in 1/Mpc is here to match the CMBFAST/CAMB output convention\n");
      fprintf(*tkfile,"#\n");
      fprintf(*tkfile,"# k (h/Mpc)       ");
      fprintf(*tkfile,"-T_cdm/k2         ");
      fprintf(*tkfile,"-T_b/k2           ");
      fprintf(*tkfile,"-T_g/k2           ");
      fprintf(*tkfile,"-T_ur/k2          ");
      fprintf(*tkfile,"-T_ncdm/k2        ");
      fprintf(*tkfile,"-T_tot/k2         ");
      fprintf(*tkfile,"\n");

    }

  }

  return _SUCCESS_;
}

/**
 * This routine writes one line with k and T_i(k)'s
 *
 * @param tkfile  Input : file pointer
 * @param one_k   Input : wavenumber
 * @param tk      Input : vector of transfer functions tk[index_tr]
 * @param tr_size Input : number of transfer functions
 * @return the error status
 */

int output_one_line_of_tk(
                          FILE * tkfile,
                          double one_k,
                          double * tk,
                          int tr_size
                          ) {

  int index_tr;

  fprintf(tkfile,"%16.10e",one_k);

  for (index_tr=0; index_tr<tr_size; index_tr++)
    fprintf(tkfile,"  %16.10e",tk[index_tr]);

  fprintf(tkfile,"\n");

  return _SUCCESS_;

}

/**
 * This routine opens one file where some background quantitites will be written, and writes
 * a heading with some general information concerning its content.
 *
 * @param pba        Input: pointer to background structure
 * @param pop        Input : pointer to output structure
 * @param backfile   Output: returned pointer to file pointer
 * @param filename   Input : name of the file
 * @return the error status
 */

int output_open_background_file(
                                struct background * pba,
                                struct output * pop,
                                FILE * * backfile,
                                FileName filename
                                ) {

  int n;
  char tmp[30]; //A fixed number here is ok, since it should just correspond to the largest string which is printed to tmp.
  class_open(*backfile,filename,"w",pop->error_message);

  if (pop->write_header == _TRUE_) {
    fprintf(*backfile,"# Table of selected background quantitites\n");
    fprintf(*backfile,"# All densities are mutiplied by (8piG/3) (below, shortcut notation (.) for this factor) \n");
    /** Length of the columntitle should be less than _OUTPUTPRECISION_+6 to be indented correctly,
        but it can be as long as . */
    fprintf(*backfile,"#");
    class_fprintf_columntitle(*backfile,"z",_TRUE_);
    class_fprintf_columntitle(*backfile,"proper time [Gyr]",_TRUE_);
    class_fprintf_columntitle(*backfile,"conf. time * c [Mpc]",_TRUE_);
    class_fprintf_columntitle(*backfile,"H / c [1/Mpc]",_TRUE_);
    class_fprintf_columntitle(*backfile,"comov. dist. [Mpc]",_TRUE_);
    class_fprintf_columntitle(*backfile,"ang.diam.dist. [Mpc]",_TRUE_);
    class_fprintf_columntitle(*backfile,"lum. dist. [Mpc]",_TRUE_);
    class_fprintf_columntitle(*backfile,"comov.snd.hrz. [Mpc]",_TRUE_);
    class_fprintf_columntitle(*backfile,"(.)rho_g [Mpc^-2]",_TRUE_);
    class_fprintf_columntitle(*backfile,"(.)rho_b [Mpc^-2]",_TRUE_);
    class_fprintf_columntitle(*backfile,"(.)rho_cdm [Mpc^-2]",pba->has_cdm);
    if (pba->has_ncdm == _TRUE_){
      for (n=0; n<pba->N_ncdm; n++){
        sprintf(tmp,"(.)rho_ncdm[%d] [Mpc^-2]",n);
        class_fprintf_columntitle(*backfile,tmp,_TRUE_);
      }
    }
    class_fprintf_columntitle(*backfile,"(.)rho_lambda [Mpc^-2]",pba->has_lambda);
    class_fprintf_columntitle(*backfile,"(.)rho_fld [Mpc^-2]",pba->has_fld);
    class_fprintf_columntitle(*backfile,"(.)rho_ur [Mpc^-2]",pba->has_ur);
    class_fprintf_columntitle(*backfile,"(.)rho_crit[Mpc^-2]",_TRUE_);
    fprintf(*backfile,"\n");
  }

  return _SUCCESS_;
}

/**
 * This routine writes one line with background quantitites
 *
 * @param pba        Input: pointer to background structure
 * @param backfile   Input : file pointer
 * @param pvecback   Input : vector of background quantitites
 * @return the error status
 */

int output_one_line_of_background(
                                  struct background * pba,
                                  FILE * backfile,
                                  double * pvecback
                                  ) {

  int n;

  fprintf(backfile," ");
  class_fprintf_double(backfile,pba->a_today/pvecback[pba->index_bg_a]-1.,_TRUE_);
  class_fprintf_double(backfile,pvecback[pba->index_bg_time]/_Gyr_over_Mpc_,_TRUE_);
  class_fprintf_double(backfile,pba->conformal_age-pvecback[pba->index_bg_conf_distance],_TRUE_);
  class_fprintf_double(backfile,pvecback[pba->index_bg_H],_TRUE_);
  class_fprintf_double(backfile,pvecback[pba->index_bg_conf_distance],_TRUE_);
  class_fprintf_double(backfile,pvecback[pba->index_bg_ang_distance],_TRUE_);
  class_fprintf_double(backfile,pvecback[pba->index_bg_lum_distance],_TRUE_);
  class_fprintf_double(backfile,pvecback[pba->index_bg_rs],_TRUE_);
  class_fprintf_double(backfile,pvecback[pba->index_bg_rho_g],_TRUE_);
  class_fprintf_double(backfile,pvecback[pba->index_bg_rho_b],_TRUE_);
  class_fprintf_double(backfile,pvecback[pba->index_bg_rho_cdm],pba->has_cdm);
  if (pba->has_ncdm == _TRUE_){
    for (n=0; n<pba->N_ncdm; n++,_TRUE_)
      class_fprintf_double(backfile,pvecback[pba->index_bg_rho_ncdm1+n],_TRUE_);
  }
  class_fprintf_double(backfile,pvecback[pba->index_bg_rho_lambda],pba->has_lambda);
  class_fprintf_double(backfile,pvecback[pba->index_bg_rho_fld],pba->has_fld);
  class_fprintf_double(backfile,pvecback[pba->index_bg_rho_ur],pba->has_ur);
  class_fprintf_double(backfile,pvecback[pba->index_bg_rho_crit],_TRUE_);

  fprintf(backfile,"\n");

  return _SUCCESS_;

}

/**
 * This routine opens one file where some thermodynamics quantitites will be written, and writes
 * a heading with some general information concerning its content.
 *
 * @param pth        Input: pointer to thermodynamics structure
 * @param pop        Input : pointer to output structure
 * @param thermofile Output: returned pointer to file pointer
 * @param filename   Input : name of the file
 * @return the error status
 */

int output_open_thermodynamics_file(
                                    struct thermo * pth,
                                    struct output * pop,
                                    FILE ** thermofile,
                                    FileName filename
                                    ) {

  class_open(*thermofile,filename,"w",pop->error_message);

  if (pop->write_header == _TRUE_) {
    fprintf(*thermofile,"# Table of selected thermodynamics quantitites\n");
    fprintf(*thermofile,"# The following notation is used in column titles:\n");
    fprintf(*thermofile,"#    x_e = electron ionisation fraction\n");
    fprintf(*thermofile,"# -kappa = optical depth\n");
    fprintf(*thermofile,"# kappa' = Thomson scattering rate, prime denotes conformal time derivatives\n");
    fprintf(*thermofile,"#      g = kappa' e^-kappa = visibility function \n");
    fprintf(*thermofile,"#     Tb = baryon temperature \n");
    fprintf(*thermofile,"#  c_b^2 = baryon sound speed squared \n");
    fprintf(*thermofile,"#  tau_d = baryon drag optical depth \n");

    /** Length of the columntitle should be less than _OUTPUTPRECISION_+6 to be indented correctly,
        but it can be as long as _COLUMNTITLE_. */

    fprintf(*thermofile,"#");
    class_fprintf_columntitle(*thermofile,"z",_TRUE_);
    class_fprintf_columntitle(*thermofile,"conf. time [Mpc]",_TRUE_);
    class_fprintf_columntitle(*thermofile,"x_e",_TRUE_);
    class_fprintf_columntitle(*thermofile,"kappa' [Mpc^-1]",_TRUE_);
    //class_fprintf_columntitle(*thermofile,"kappa''",_TRUE_);
    //class_fprintf_columntitle(*thermofile,"kappa'''",_TRUE_);
    class_fprintf_columntitle(*thermofile,"exp(-kappa)",_TRUE_);
    class_fprintf_columntitle(*thermofile,"g [Mpc^-1]",_TRUE_);
    //class_fprintf_columntitle(*thermofile,"g'",_TRUE_);
    //class_fprintf_columntitle(*thermofile,"g''",_TRUE_);
    class_fprintf_columntitle(*thermofile,"Tb [K]",_TRUE_);
    class_fprintf_columntitle(*thermofile,"c_b^2",_TRUE_);
    class_fprintf_columntitle(*thermofile,"tau_d",_TRUE_);
    //class_fprintf_columntitle(*thermofile,"max. rate",_TRUE_);
    fprintf(*thermofile,"\n");
  }

  return _SUCCESS_;
}

/**
 * This routine writes one line with thermodynamics quantitites
 *
 * @param pth        Input : pointer to thermodynamics structure
 * @param thermofile Input : file pointer
 * @param tau        Input : conformal time
 * @param z          Input : redshift
 * @param pvecthermo Input : vector of thermodynamics quantitites
 * @return the error status
 */

int output_one_line_of_thermodynamics(
                                      struct thermo * pth,
                                      FILE * thermofile,
                                      double tau,
                                      double z,
                                      double * pvecthermo
                                      ) {

  fprintf(thermofile," ");
  class_fprintf_double(thermofile,z,_TRUE_);
  class_fprintf_double(thermofile,tau,_TRUE_);
  class_fprintf_double(thermofile,pvecthermo[pth->index_th_xe],_TRUE_);
  class_fprintf_double(thermofile,pvecthermo[pth->index_th_dkappa],_TRUE_);
  //class_fprintf_double(thermofile,pvecthermo[pth->index_th_ddkappa],_TRUE_);
  //class_fprintf_double(thermofile,pvecthermo[pth->index_th_dddkappa],_TRUE_);
  class_fprintf_double(thermofile,pvecthermo[pth->index_th_exp_m_kappa],_TRUE_);
  class_fprintf_double(thermofile,pvecthermo[pth->index_th_g],_TRUE_);
  //class_fprintf_double(thermofile,pvecthermo[pth->index_th_dg],_TRUE_);
  //class_fprintf_double(thermofile,pvecthermo[pth->index_th_ddg],_TRUE_);
  class_fprintf_double(thermofile,pvecthermo[pth->index_th_Tb],_TRUE_);
  class_fprintf_double(thermofile,pvecthermo[pth->index_th_cb2],_TRUE_);
  class_fprintf_double(thermofile,pvecthermo[pth->index_th_tau_d],_TRUE_);
  //class_fprintf_double(thermofile,pvecthermo[pth->index_th_rate],_TRUE_);

  fprintf(thermofile,"\n");

  return _SUCCESS_;

}

/**
 * This routine opens one file where the primordial spectrum/spectra will be written,
 * and writes a heading with some general information concerning its content.
 *
 * @param ppt        Input: pointer to perturbation structure
 * @param ppm        Input: pointer to primordial structure
 * @param pop        Input : pointer to output structure
 * @param outputfile Output: returned pointer to file pointer
 * @param filename   Input : name of the file
 * @return the error status
 */

int output_open_primordial_file(
                                struct perturbs * ppt,
                                struct primordial * ppm,
                                struct output * pop,
                                FILE * * outputfile,
                                FileName filename
                                ) {

  class_open(*outputfile,filename,"w",pop->error_message);

  if (pop->write_header == _TRUE_) {
    fprintf(*outputfile,"# Dimensionless primordial spectrum, equal to [k^3/2pi^2] P(k) \n");
    fprintf(*outputfile,"                k [1/Mpc]");
    fprintf(*outputfile,"              P_scalar(k)");
    if (ppt->has_tensors == _TRUE_)
      fprintf(*outputfile,"              P_tensor(k)");
    fprintf(*outputfile,"\n");
  }

  return _SUCCESS_;
}

/**
 * This routine writes one line with the primordial spectrum/spectra
 *
 * @param ppt        Input: pointer to perturbation structure
 * @param ppm        Input: pointer to primordial structure
 * @param outputfile   Input : file pointer
 * @param index_k    Input: index for wavenumebr in ppm structure
 * @return the error status
 */

int output_one_line_of_primordial(
                                  struct perturbs * ppt,
                                  struct primordial * ppm,
                                  FILE * outputfile,
                                  int index_k
                                  ) {

  fprintf(outputfile,"%25.12e",exp(ppm->lnk[index_k]));
  fprintf(outputfile,"%25.12e",exp(ppm->lnpk[ppt->index_md_scalars][index_k]));
  if (ppt->has_tensors == _TRUE_)
    fprintf(outputfile,"%25.12e",exp(ppm->lnpk[ppt->index_md_tensors][index_k]));
  fprintf(outputfile,"\n");

  return _SUCCESS_;

}
