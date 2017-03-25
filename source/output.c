/** @file output.c Documented output module
 *
 * Julien Lesgourgues, 26.08.2010
 *
 * This module writes the output in files.
 *
 * The following functions can be called from other modules or from the main:
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
 * @param pth Input: pointer to thermodynamics structure
 * @param ppt Input: pointer perturbation structure
 * @param ppm Input: pointer to primordial structure
 * @param ptr Input: pointer to transfer structure
 * @param psp Input: pointer to spectra structure
 * @param pnl Input: pointer to nonlinear structure
 * @param ple Input: pointer to lensing structure
 * @param pop Input: pointer to output structure
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

  /** - deal with all anisotropy power spectra \f$ C_l\f$'s */

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

  /** - deal with background quantities */

  if (pop->write_background == _TRUE_) {

    class_call(output_background(pba,pop),
               pop->error_message,
               pop->error_message);

  }

  /** - deal with thermodynamics quantities */

  if (pop->write_thermodynamics == _TRUE_) {

    class_call(output_thermodynamics(pba,pth,pop),
               pop->error_message,
               pop->error_message);

  }

  /** - deal with perturbation quantities */

  if (pop->write_perturbations == _TRUE_) {

    class_call(output_perturbations(pba,ppt,pop),
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
 * This routines writes the output in files for anisotropy power spectra \f$ C_l\f$'s.
 *
 * @param pba Input: pointer to background structure (needed for \f$ T_{cmb}\f$)
 * @param ppt Input: pointer perturbation structure
 * @param psp Input: pointer to spectra structure
 * @param ple Input: pointer to lensing structure
 * @param pop Input: pointer to output structure
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

  /** - first, allocate all arrays of files and \f$ C_l\f$'s */

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

  /** - third, perform loop over l. For each multipole, get all \f$ C_l\f$'s
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

  /** - finally, close files and free arrays of files and \f$ C_l\f$'s */

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
 * @param ppt Input: pointer perturbation structure
 * @param psp Input: pointer to spectra structure
 * @param pop Input: pointer to output structure
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

    /** - second, open only the relevant files and write a heading in each of them */

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
 * @param ppt Input: pointer perturbation structure
 * @param psp Input: pointer to spectra structure
 * @param pop Input: pointer to output structure
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
 * This routines writes the output in files for matter transfer functions \f$ T_i(k)\f$'s.
 *
 * @param pba Input: pointer to background structure (needed for calling spectra_pk_at_z())
 * @param ppt Input: pointer perturbation structure
 * @param psp Input: pointer to spectra structure
 * @param pop Input: pointer to output structure
 */

int output_tk(
              struct background * pba,
              struct perturbs * ppt,
              struct spectra * psp,
              struct output * pop
              ) {

  /** Summary: */

  /** - define local variables */
  char titles[_MAXTITLESTRINGLENGTH_]={0};
  double * data;
  int size_data, number_of_titles;

  FILE * tkfile;

  int index_md;
  int index_ic;
  int index_z;

  double z;

  FileName file_name;
  FileName redshift_suffix;
  char first_line[_LINE_LENGTH_MAX_];
  FileName ic_suffix;

  index_md=ppt->index_md_scalars;

  if (pop->output_format == camb_format) {

    class_test(pba->N_ncdm>1,
               pop->error_message,
               "you wish to output the transfer functions in CMBFAST/CAMB format but you have more than one non-cold dark matter (ncdm) species. The two are not compatible (since CMBFAST/CAMB only have one ncdm species): switch to CLASS output format or keep only on ncdm species");

    class_test(ppt->has_velocity_transfers == _TRUE_,
               pop->error_message,
               "you wish to output the transfer functions in CMBFAST/CAMB format, but you requested velocity transfer functions. The two are not compatible (since CMBFAST/CAMB do not compute velocity transfer functions): switch to CLASS output format, or ask only for density transfer function");
  }


  class_call(spectra_output_tk_titles(pba,ppt,pop->output_format,titles),
             pba->error_message,
             pop->error_message);
  number_of_titles = get_number_of_titles(titles);
  size_data = number_of_titles*psp->ln_k_size;

  class_alloc(data, sizeof(double)*psp->ic_size[index_md]*size_data, pop->error_message);

  for (index_z = 0; index_z < pop->z_pk_num; index_z++) {

    z = pop->z_pk[index_z];

    /** - first, check that requested redshift z_pk is consistent */

    class_test((pop->z_pk[index_z] > psp->z_max_pk),
               pop->error_message,
               "T_i(k,z) computed up to z=%f but requested at z=%f. Must increase z_max_pk in precision file.",psp->z_max_pk,pop->z_pk[index_z]);

    if (pop->z_pk_num == 1)
      redshift_suffix[0]='\0';
    else
      sprintf(redshift_suffix,"z%d_",index_z+1);

    /** - second, open only the relevant files, and write a heading in each of them */

    class_call(spectra_output_tk_data(pba,
                                      ppt,
                                      psp,
                                      pop->output_format,
                                      pop->z_pk[index_z],
                                      number_of_titles,
                                      data
                                      ),
               psp->error_message,
               pop->error_message);

    for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++) {

      class_call(spectra_firstline_and_ic_suffix(ppt, index_ic, first_line, ic_suffix),
                 pop->error_message, pop->error_message);

      if ((ppt->has_ad == _TRUE_) && (ppt->ic_size[index_md] == 1) )
        sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"tk.dat");
      else
        sprintf(file_name,"%s%s%s%s%s",pop->root,redshift_suffix,"tk_",ic_suffix,".dat");

      class_open(tkfile, file_name, "w", pop->error_message);

      if (pop->write_header == _TRUE_) {
        if (pop->output_format == class_format) {
          fprintf(tkfile,"# Transfer functions T_i(k) %sat redshift z=%g\n",first_line,z);
          fprintf(tkfile,"# for k=%g to %g h/Mpc,\n",exp(psp->ln_k[0])/pba->h,exp(psp->ln_k[psp->ln_k_size-1])/pba->h);
          fprintf(tkfile,"# number of wavenumbers equal to %d\n",psp->ln_k_size);
          if (ppt->has_density_transfers == _TRUE_) {
            fprintf(tkfile,"# d_i   stands for (delta rho_i/rho_i)(k,z) with above normalization \n");
            fprintf(tkfile,"# d_tot stands for (delta rho_tot/rho_tot)(k,z) with rho_Lambda NOT included in rho_tot\n");
            fprintf(tkfile,"# (note that this differs from the transfer function output from CAMB/CMBFAST, which gives the same\n");
            fprintf(tkfile,"#  quantities divided by -k^2 with k in Mpc^-1; use format=camb to match CAMB)\n");
          }
          if (ppt->has_velocity_transfers == _TRUE_) {
            fprintf(tkfile,"# t_i   stands for theta_i(k,z) with above normalization \n");
            fprintf(tkfile,"# t_tot stands for (sum_i [rho_i+p_i] theta_i)/(sum_i [rho_i+p_i]))(k,z)\n");
          }
          fprintf(tkfile,"#\n");
        }
        else if (pop->output_format == camb_format) {

          fprintf(tkfile,"# Rescaled matter transfer functions [-T_i(k)/k^2] %sat redshift z=%g\n",first_line,z);
          fprintf(tkfile,"# for k=%g to %g h/Mpc,\n",exp(psp->ln_k[0])/pba->h,exp(psp->ln_k[psp->ln_k_size-1])/pba->h);
          fprintf(tkfile,"# number of wavenumbers equal to %d\n",psp->ln_k_size);
          fprintf(tkfile,"# T_i   stands for (delta rho_i/rho_i)(k,z) with above normalization \n");
          fprintf(tkfile,"# The rescaling factor [-1/k^2] with k in 1/Mpc is here to match the CMBFAST/CAMB output convention\n");
          fprintf(tkfile,"#\n");
          fprintf(tkfile,"#");
          fprintf(tkfile,"\n");

        }
      }

      output_print_data(tkfile,
                        titles,
                        data+index_ic*size_data,
                        size_data);

      /** - free memory and close files */
      fclose(tkfile);

    }

  }

  free(data);

  return _SUCCESS_;

}

int output_background(
                      struct background * pba,
                      struct output * pop
                      ) {

  FILE * backfile;
  FileName file_name;

  char titles[_MAXTITLESTRINGLENGTH_]={0};
  double * data;
  int size_data, number_of_titles;

  class_call(background_output_titles(pba,titles),
             pba->error_message,
             pop->error_message);
  number_of_titles = get_number_of_titles(titles);
  size_data = number_of_titles*pba->bt_size;
  class_alloc(data,sizeof(double)*size_data,pop->error_message);
  class_call(background_output_data(pba,
                                    number_of_titles,
                                    data),
             pba->error_message,
             pop->error_message);

  sprintf(file_name,"%s%s",pop->root,"background.dat");
  class_open(backfile,file_name,"w",pop->error_message);

  if (pop->write_header == _TRUE_) {
    fprintf(backfile,"# Table of selected background quantities\n");
    fprintf(backfile,"# All densities are multiplied by (8piG/3) (below, shortcut notation (.) for this factor) \n");
    fprintf(backfile,"# Densities are in units [Mpc^-2] while all distances are in [Mpc]. \n");
    if (pba->has_scf == _TRUE_){
      fprintf(backfile,"# The units of phi, tau in the derivatives and the potential V are the following:\n");
      fprintf(backfile,"# --> phi is given in units of the reduced Planck mass m_Pl = (8 pi G)^(-1/2)\n");
      fprintf(backfile,"# --> tau in the derivative of V(phi) is given in units of Mpc.\n");
      fprintf(backfile,"# --> the potential V(phi) is given in units of m_Pl^2/Mpc^2.\n");
    }
  }

  output_print_data(backfile,
                    titles,
                    data,
                    size_data);

  free(data);
  fclose(backfile);

  return _SUCCESS_;

}

int output_thermodynamics(
                          struct background * pba,
                          struct thermo * pth,
                          struct output * pop
                      ) {

  FileName file_name;
  FILE * thermofile;
  char titles[_MAXTITLESTRINGLENGTH_]={0};
  double * data;
  int size_data, number_of_titles;

  class_call(thermodynamics_output_titles(pba,pth,titles),
             pth->error_message,
             pop->error_message);
  number_of_titles = get_number_of_titles(titles);
  size_data = number_of_titles*pth->tt_size;
  class_alloc(data,sizeof(double)*size_data,pop->error_message);
  class_call(thermodynamics_output_data(pba,
                                        pth,
                                        number_of_titles,
                                        data),
             pth->error_message,
             pop->error_message);

  sprintf(file_name,"%s%s",pop->root,"thermodynamics.dat");
  class_open(thermofile,file_name,"w",pop->error_message);

  if (pop->write_header == _TRUE_) {
    fprintf(thermofile,"# Table of selected thermodynamics quantities\n");
    fprintf(thermofile,"# The following notation is used in column titles:\n");
    fprintf(thermofile,"#    x_e = electron ionization fraction\n");
    fprintf(thermofile,"# -kappa = optical depth\n");
    fprintf(thermofile,"# kappa' = Thomson scattering rate, prime denotes conformal time derivatives\n");
    fprintf(thermofile,"#      g = kappa' e^-kappa = visibility function \n");
    fprintf(thermofile,"#     Tb = baryon temperature \n");
    fprintf(thermofile,"#  c_b^2 = baryon sound speed squared \n");
    fprintf(thermofile,"#  tau_d = baryon drag optical depth \n");
    if (pth->compute_damping_scale == _TRUE_) {
      fprintf(thermofile,"#  r_d = simplest analytic approximation to photon comoving damping scale \n");
    }
  }

  output_print_data(thermofile,
                    titles,
                    data,
                    size_data);

  free(data);
  fclose(thermofile);

  return _SUCCESS_;

}


int output_perturbations(
                         struct background * pba,
                         struct perturbs * ppt,
                         struct output * pop
                         ) {

  FILE * out;
  FileName file_name;
  int index_ikout, index_md;
  double k;

  for (index_ikout=0; index_ikout<ppt->k_output_values_num; index_ikout++){

    if (ppt->has_scalars == _TRUE_){
      index_md = ppt->index_md_scalars;
      k = ppt->k[index_md][ppt->index_k_output_values[index_md*ppt->k_output_values_num+index_ikout]];
      sprintf(file_name,"%s%s%d%s",pop->root,"perturbations_k",index_ikout,"_s.dat");
      class_open(out, file_name, "w", ppt->error_message);
      fprintf(out,"#scalar perturbations for mode k = %.*e Mpc^(-1)\n",_OUTPUTPRECISION_,k);
      output_print_data(out,
                        ppt->scalar_titles,
                        ppt->scalar_perturbations_data[index_ikout],
                        ppt->size_scalar_perturbation_data[index_ikout]);

      fclose(out);
    }
    if (ppt->has_vectors == _TRUE_){
      index_md = ppt->index_md_vectors;
      k = ppt->k[index_md][ppt->index_k_output_values[index_md*ppt->k_output_values_num+index_ikout]];
      sprintf(file_name,"%s%s%d%s",pop->root,"perturbations_k",index_ikout,"_v.dat");
      class_open(out, file_name, "w", ppt->error_message);
      fprintf(out,"#vector perturbations for mode k = %.*e Mpc^(-1)\n",_OUTPUTPRECISION_,k);
      output_print_data(out,
                        ppt->vector_titles,
                        ppt->vector_perturbations_data[index_ikout],
                        ppt->size_vector_perturbation_data[index_ikout]);

      fclose(out);
    }
    if (ppt->has_tensors == _TRUE_){
      index_md = ppt->index_md_tensors;
      k = ppt->k[index_md][ppt->index_k_output_values[index_md*ppt->k_output_values_num+index_ikout]];
      sprintf(file_name,"%s%s%d%s",pop->root,"perturbations_k",index_ikout,"_t.dat");
      class_open(out, file_name, "w", ppt->error_message);
      fprintf(out,"#tensor perturbations for mode k = %.*e Mpc^(-1)\n",_OUTPUTPRECISION_,k);
      output_print_data(out,
                        ppt->tensor_titles,
                        ppt->tensor_perturbations_data[index_ikout],
                        ppt->size_tensor_perturbation_data[index_ikout]);

      fclose(out);
    }


  }
  return _SUCCESS_;

}

int output_primordial(
                      struct perturbs * ppt,
                      struct primordial * ppm,
                      struct output * pop
                      ) {
  FileName file_name;
  FILE * out;
  char titles[_MAXTITLESTRINGLENGTH_]={0};
  double * data;
  int size_data, number_of_titles;

  sprintf(file_name,"%s%s",pop->root,"primordial_Pk.dat");

  class_call(primordial_output_titles(ppt,ppm,titles),
             ppm->error_message,
             pop->error_message);
  number_of_titles = get_number_of_titles(titles);
  size_data = number_of_titles*ppm->lnk_size;
  class_alloc(data,sizeof(double)*size_data,pop->error_message);
  class_call(primordial_output_data(ppt,
                                    ppm,
                                    number_of_titles,
                                    data),
             ppm->error_message,
             pop->error_message);

  class_open(out,file_name,"w",pop->error_message);
  if (pop->write_header == _TRUE_) {
    fprintf(out,"# Dimensionless primordial spectrum, equal to [k^3/2pi^2] P(k) \n");
  }

  output_print_data(out,
                    titles,
                    data,
                    size_data);

  free(data);
  fclose(out);

  return _SUCCESS_;
}


int output_print_data(FILE *out,
                      char titles[_MAXTITLESTRINGLENGTH_],
                      double *dataptr,
                      int size_dataptr){
  int colnum=1, number_of_titles;
  int index_title, index_tau;
  char thetitle[_MAXTITLESTRINGLENGTH_];
  char *pch;

  /** Summary*/

  /** - First we print the titles */
  fprintf(out,"#");

  strcpy(thetitle,titles);
  pch = strtok(thetitle,_DELIMITER_);
  while (pch != NULL){
    class_fprintf_columntitle(out, pch, _TRUE_, colnum);
    pch = strtok(NULL,_DELIMITER_);
  }
  fprintf(out,"\n");

  /** - Then we print the data */
  number_of_titles = colnum-1;
  if (number_of_titles>0){
    for (index_tau=0; index_tau<size_dataptr/number_of_titles; index_tau++){
      fprintf(out," ");
      for (index_title=0; index_title<number_of_titles; index_title++){
        class_fprintf_double(out, dataptr[index_tau*number_of_titles+index_title], _TRUE_);
      }
      fprintf(out,"\n");
    }
  }
  return _SUCCESS_;
}


/**
 * This routine opens one file where some \f$ C_l\f$'s will be written, and writes
 * a heading with some general information concerning its content.
 *
 * @param psp        Input: pointer to spectra structure
 * @param pop        Input: pointer to output structure
 * @param clfile     Output: returned pointer to file pointer
 * @param filename   Input: name of the file
 * @param first_line Input: text describing the content (mode, initial condition..)
 * @param lmax       Input: last multipole in the file (the first one is assumed to be 2)
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
  /** Summary */

  int index_d1,index_d2;
  int colnum = 1;
  char tmp[60]; //A fixed number here is ok, since it should just correspond to the largest string which is printed to tmp.

  class_open(*clfile,filename,"w",pop->error_message);

  if (pop->write_header == _TRUE_) {

    /** - First we deal with the entries that are dependent of format type */

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

    if (psp->has_pp == _TRUE_) {
      if (pop->output_format == class_format) {
        fprintf(*clfile,"# -> for CMB lensing (phi), these are C_l^phi-phi for the lensing potential.\n");
      }
      if (pop->output_format == camb_format) {
        fprintf(*clfile,"# -> for CMB lensing (d), these are C_l^dd for the deflection field.\n");
      }
    }

    if (psp->has_ll == _TRUE_) {
      fprintf(*clfile,"# -> for galaxy lensing (lens[i]), these are C_l^phi-phi for the lensing potential.\n");
    }

    if (psp->has_pp == _TRUE_ || psp->has_ll == _TRUE_) {
      fprintf(*clfile,"#    Remember the conversion factors:\n");
      fprintf(*clfile,"#    C_l^dd (deflection) = l(l+1) C_l^phi-phi\n");
      fprintf(*clfile,"#    C_l^gg (shear/convergence) = 1/4 (l(l+1))^2 C_l^phi-phi\n");
    }

    fprintf(*clfile,"#\n");

    if (0==1){
      fprintf(*clfile,"#");
      class_fprintf_columntitle(*clfile,"l",_TRUE_,colnum);
    }
    else{
      fprintf(*clfile,"# 1:l ");
      colnum++;
    }
    if (pop->output_format == class_format) {
      class_fprintf_columntitle(*clfile,"TT",psp->has_tt,colnum);
      class_fprintf_columntitle(*clfile,"EE",psp->has_ee,colnum);
      class_fprintf_columntitle(*clfile,"TE",psp->has_te,colnum);
      class_fprintf_columntitle(*clfile,"BB",psp->has_bb,colnum);
      class_fprintf_columntitle(*clfile,"phiphi",psp->has_pp,colnum);
      class_fprintf_columntitle(*clfile,"TPhi",psp->has_tp,colnum);
      class_fprintf_columntitle(*clfile,"Ephi",psp->has_ep,colnum);
    }
    else if (pop->output_format == camb_format) {
      class_fprintf_columntitle(*clfile,"TT",psp->has_tt,colnum);
      class_fprintf_columntitle(*clfile,"EE",psp->has_ee,colnum);
      class_fprintf_columntitle(*clfile,"BB",psp->has_bb,colnum);
      class_fprintf_columntitle(*clfile,"TE",psp->has_te,colnum);
      class_fprintf_columntitle(*clfile,"dd",psp->has_pp,colnum);
      class_fprintf_columntitle(*clfile,"dT",psp->has_tp,colnum);
      class_fprintf_columntitle(*clfile,"dE",psp->has_ep,colnum);
    }

    /** - Next deal with entries that are independent of format type */

    if (psp->has_dd == _TRUE_){
      for (index_d1=0; index_d1<psp->d_size; index_d1++){
        for (index_d2=index_d1; index_d2<=MIN(index_d1+psp->non_diag,psp->d_size-1); index_d2++){
          sprintf(tmp,"dens[%d]-dens[%d]",index_d1+1,index_d2+1);
          class_fprintf_columntitle(*clfile,tmp,_TRUE_,colnum);
        }
      }
    }
    if (psp->has_td == _TRUE_){
      for (index_d1=0; index_d1<psp->d_size; index_d1++){
        sprintf(tmp,"T-dens[%d]",index_d1+1);
        class_fprintf_columntitle(*clfile,tmp,_TRUE_,colnum);
      }
    }
    if (psp->has_pd == _TRUE_){
      for (index_d1=0; index_d1<psp->d_size; index_d1++){
        sprintf(tmp,"phi-dens[%d]",index_d1+1);
        class_fprintf_columntitle(*clfile,tmp,_TRUE_,colnum);
      }
    }
    if (psp->has_ll == _TRUE_){
      for (index_d1=0; index_d1<psp->d_size; index_d1++){
        for (index_d2=index_d1; index_d2<=MIN(index_d1+psp->non_diag,psp->d_size-1); index_d2++){
          sprintf(tmp,"lens[%d]-lens[%d]",index_d1+1,index_d2+1);
          class_fprintf_columntitle(*clfile,tmp,_TRUE_,colnum);
        }
      }
    }
    if (psp->has_tl == _TRUE_){
      for (index_d1=0; index_d1<psp->d_size; index_d1++){
        sprintf(tmp,"T-lens[%d]",index_d1+1);
        class_fprintf_columntitle(*clfile,tmp,_TRUE_,colnum);
      }
    }
    if (psp->has_dl == _TRUE_){
      for (index_d1=0; index_d1<psp->d_size; index_d1++){
        for (index_d2=MAX(index_d1-psp->non_diag,0); index_d2<=MIN(index_d1+psp->non_diag,psp->d_size-1); index_d2++) {
          sprintf(tmp,"dens[%d]-lens[%d]",index_d1+1,index_d2+1);
          class_fprintf_columntitle(*clfile,tmp,_TRUE_,colnum);
        }
      }
    }
    fprintf(*clfile,"\n");
  }

  return _SUCCESS_;

}

/**
 * This routine write one line with l and all \f$ C_l\f$'s for all types (TT, TE...)
 *
 * @param pba        Input: pointer to background structure (needed for \f$ T_{cmb}\f$)
 * @param psp        Input: pointer to spectra structure
 * @param pop        Input: pointer to output structure
 * @param clfile  Input: file pointer
 * @param l       Input: multipole
 * @param cl      Input: \f$ C_l\f$'s for all types
 * @param ct_size Input: number of types
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
  int index_ct, index_ct_rest;
  double factor;

  factor = l*(l+1)/2./_PI_;

  fprintf(clfile," ");

  if (0==1){
    class_fprintf_int(clfile, (int)l, _TRUE_);
  }
  else{
    fprintf(clfile,"%4d ",(int)l);
  }

  if (pop->output_format == class_format) {

    for (index_ct=0; index_ct < ct_size; index_ct++) {
      class_fprintf_double(clfile, factor*cl[index_ct], _TRUE_);
    }
    fprintf(clfile,"\n");
  }

  if (pop->output_format == camb_format) {
    class_fprintf_double(clfile, factor*pow(pba->T_cmb*1.e6,2)*cl[psp->index_ct_tt], psp->has_tt);
    class_fprintf_double(clfile, factor*pow(pba->T_cmb*1.e6,2)*cl[psp->index_ct_ee], psp->has_ee);
    class_fprintf_double(clfile, factor*pow(pba->T_cmb*1.e6,2)*cl[psp->index_ct_bb], psp->has_bb);
    class_fprintf_double(clfile, factor*pow(pba->T_cmb*1.e6,2)*cl[psp->index_ct_te], psp->has_te);
    class_fprintf_double(clfile, l*(l+1)*factor*cl[psp->index_ct_pp], psp->has_pp);
    class_fprintf_double(clfile, sqrt(l*(l+1))*factor*pba->T_cmb*1.e6*cl[psp->index_ct_tp], psp->has_tp);
    class_fprintf_double(clfile, sqrt(l*(l+1))*factor*pba->T_cmb*1.e6*cl[psp->index_ct_ep], psp->has_ep);
    index_ct_rest = 0;
    if (psp->has_tt == _TRUE_)
      index_ct_rest++;
    if (psp->has_ee == _TRUE_)
      index_ct_rest++;
    if (psp->has_bb == _TRUE_)
      index_ct_rest++;
    if (psp->has_te == _TRUE_)
      index_ct_rest++;
    if (psp->has_pp == _TRUE_)
      index_ct_rest++;
    if (psp->has_tp == _TRUE_)
      index_ct_rest++;
    if (psp->has_ep == _TRUE_)
      index_ct_rest++;
    /* Now print the remaining (if any) entries:*/
    for (index_ct=index_ct_rest; index_ct < ct_size; index_ct++) {
      class_fprintf_double(clfile, factor*cl[index_ct], _TRUE_);
    }

    fprintf(clfile,"\n");

  }
  return _SUCCESS_;

}

/**
 * This routine opens one file where some P(k)'s will be written, and writes
 * a heading with some general information concerning its content.
 *
 * @param pba        Input: pointer to background structure (needed for h)
 * @param psp        Input: pointer to spectra structure
 * @param pop        Input: pointer to output structure
 * @param pkfile     Output: returned pointer to file pointer
 * @param filename   Input: name of the file
 * @param first_line Input: text describing the content (initial conditions, ...)
 * @param z          Input: redshift of the output
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

  int colnum = 1;
  class_open(*pkfile,filename,"w",pop->error_message);

  if (pop->write_header == _TRUE_) {
    fprintf(*pkfile,"# Matter power spectrum P(k) %sat redshift z=%g\n",first_line,z);
    fprintf(*pkfile,"# for k=%g to %g h/Mpc,\n",
            exp(psp->ln_k[0])/pba->h,
            exp(psp->ln_k[psp->ln_k_size-1])/pba->h);
    fprintf(*pkfile,"# number of wavenumbers equal to %d\n",psp->ln_k_size);

    fprintf(*pkfile,"#");
    class_fprintf_columntitle(*pkfile,"k (h/Mpc)",_TRUE_,colnum);
    class_fprintf_columntitle(*pkfile,"P (Mpc/h)^3",_TRUE_,colnum);

    fprintf(*pkfile,"\n");
  }

  return _SUCCESS_;
}

/**
 * This routine writes one line with k and P(k)
 *
 * @param pkfile  Input: file pointer
 * @param one_k   Input: wavenumber
 * @param one_pk  Input: matter power spectrum
 * @return the error status
 */

int output_one_line_of_pk(
                          FILE * pkfile,
                          double one_k,
                          double one_pk
                          ) {

  fprintf(pkfile," ");
  class_fprintf_double(pkfile,one_k,_TRUE_);
  class_fprintf_double(pkfile,one_pk,_TRUE_);
  fprintf(pkfile,"\n");

  return _SUCCESS_;

}
