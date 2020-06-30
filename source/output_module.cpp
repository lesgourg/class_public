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

#include "background_module.h"
#include "thermodynamics_module.h"
#include "perturbations_module.h"
#include "primordial_module.h"
#include "nonlinear_module.h"
#include "lensing_module.h"
#include "spectra_module.h"
#include "output_module.h"

OutputModule::OutputModule(InputModulePtr input_module, BackgroundModulePtr background_module, ThermodynamicsModulePtr thermodynamics_module, PerturbationsModulePtr perturbations_module, PrimordialModulePtr primordial_module, NonlinearModulePtr nonlinear_module, SpectraModulePtr spectra_module, LensingModulePtr lensing_module)
: BaseModule(std::move(input_module))
, background_module_(std::move(background_module))
, thermodynamics_module_(std::move(thermodynamics_module))
, perturbations_module_(std::move(perturbations_module))
, primordial_module_(std::move(primordial_module))
, nonlinear_module_(std::move(nonlinear_module))
, spectra_module_(std::move(spectra_module))
, lensing_module_(std::move(lensing_module)) {

  if (output_init() != _SUCCESS_) {
    throw std::runtime_error(error_message_);
  }
}


int OutputModule::output_total_cl_at_l(
                         int l,
                         double * cl
                         ){

  double ** cl_md_ic; /* array with argument
                         cl_md_ic[index_md][index_ic1_ic2*spectra_module_->ct_size_+index_ct] */

  double ** cl_md;    /* array with argument
                         cl_md[index_md][index_ct] */

  int index_md;

  if (ple->has_lensed_cls == _TRUE_) {
    class_call(lensing_module_->lensing_cl_at_l(l, cl),
               lensing_module_->error_message_,
               error_message_);
  }
  else {

    class_alloc(cl_md_ic,
                spectra_module_->md_size_*sizeof(double *),
                error_message_);

    class_alloc(cl_md,
                spectra_module_->md_size_*sizeof(double *),
                error_message_);

    for (index_md = 0; index_md < spectra_module_->md_size_; index_md++) {

      if (spectra_module_->md_size_ > 1)

        class_alloc(cl_md[index_md],
                    spectra_module_->ct_size_*sizeof(double),
                    error_message_);

      if (spectra_module_->ic_size_[index_md] > 1)

        class_alloc(cl_md_ic[index_md],
                    spectra_module_->ic_ic_size_[index_md]*spectra_module_->ct_size_*sizeof(double),
                    error_message_);
    }

    class_call(spectra_module_->spectra_cl_at_l((double)l, cl, cl_md, cl_md_ic),
               psp->error_message,
               error_message_);

    for (index_md = 0; index_md < spectra_module_->md_size_; index_md++) {

      if (spectra_module_->md_size_ > 1)
        free(cl_md[index_md]);

      if (spectra_module_->ic_size_[index_md] > 1)
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

int OutputModule::output_init() {

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

    class_call(output_cl(),
               error_message_,
               error_message_);
  }

  /** - deal with all Fourier matter power spectra P(k)'s */

  if (ppt->has_pk_matter == _TRUE_) {

    class_call(output_pk(pk_linear),
               error_message_,
               error_message_);

    if (pnl->method != nl_none) {

      class_call(output_pk(pk_nonlinear),
                 error_message_,
                 error_message_);

    }
  }

  /** - deal with density and matter power spectra */

  if ((ppt->has_density_transfers == _TRUE_) || (ppt->has_velocity_transfers == _TRUE_)) {

    class_call(output_tk(),
               error_message_,
               error_message_);

  }

  /** - deal with background quantities */

  if (pop->write_background == _TRUE_) {

    class_call(output_background(),
               error_message_,
               error_message_);

  }

  /** - deal with thermodynamics quantities */

  if (pop->write_thermodynamics == _TRUE_) {

    class_call(output_thermodynamics(),
               error_message_,
               error_message_);

  }

  /** - deal with perturbation quantities */

  if (pop->write_perturbations == _TRUE_) {

    class_call(output_perturbations(),
               error_message_,
               error_message_);

  }

  /** - deal with primordial spectra */

  if (pop->write_primordial == _TRUE_) {

    class_call(output_primordial(),
               error_message_,
               error_message_);

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

int OutputModule::output_cl() {

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
                         cl_md_ic[index_md][index_ic1_ic2*spectra_module_->ct_size_+index_ct] */

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
              spectra_module_->md_size_*sizeof(FILE * *),
              error_message_);

  class_alloc(cl_md_ic,
              spectra_module_->md_size_*sizeof(double *),
              error_message_);

  class_alloc(out_md,
              spectra_module_->md_size_*sizeof(FILE *),
              error_message_);

  class_alloc(cl_md,
              spectra_module_->md_size_*sizeof(double *),
              error_message_);

  for (index_md = 0; index_md < perturbations_module_->md_size_; index_md++) {

    class_alloc(out_md_ic[index_md],
                spectra_module_->ic_ic_size_[index_md]*sizeof(FILE *),
                error_message_);

  }

  /** - second, open only the relevant files, and write a heading in each of them */

  sprintf(file_name,"%s%s",pop->root,"cl.dat");

  class_call(output_open_cl_file(&out,
                                 file_name,
                                 "total [l(l+1)/2pi] C_l's",
                                 spectra_module_->l_max_tot_
                                 ),
             error_message_,
             error_message_);

  class_alloc(cl_tot,
              spectra_module_->ct_size_*sizeof(double),
              error_message_);


  if (ple->has_lensed_cls == _TRUE_) {

    sprintf(file_name,"%s%s",pop->root,"cl_lensed.dat");

    class_call(output_open_cl_file(&out_lensed,
                                   file_name,
                                   "total lensed [l(l+1)/2pi] C_l's",
                                   lensing_module_->l_lensed_max_
                                   ),
               error_message_,
               error_message_);
  }

  if (perturbations_module_->md_size_ > 1) {

    for (index_md = 0; index_md < perturbations_module_->md_size_; index_md++) {

      if (_scalarsEXT_) {

        sprintf(file_name,"%s%s",pop->root,"cls.dat");
        strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar mode");

      }

      if (_tensorsEXT_) {

        sprintf(file_name,"%s%s",pop->root,"clt.dat");
        strcpy(first_line,"[l(l+1)/2pi] C_l's for tensor mode");

      }

      class_call(output_open_cl_file(&(out_md[index_md]),
                                     file_name,
                                     first_line,
                                     spectra_module_->l_max_[index_md]
                                     ),
                 error_message_,
                 error_message_);

      class_alloc(cl_md[index_md],
                  spectra_module_->ct_size_*sizeof(double),
                  error_message_);

    }
  }

  for (index_md = 0; index_md < perturbations_module_->md_size_; index_md++) {

    if (perturbations_module_->ic_size_[index_md] > 1) {

      for (index_ic1 = 0; index_ic1 < perturbations_module_->ic_size_[index_md]; index_ic1++) {

        for (index_ic2 = index_ic1; index_ic2 < perturbations_module_->ic_size_[index_md]; index_ic2++) {

          if (_scalarsEXT_) {

            if ((ppt->has_ad == _TRUE_) &&
                (index_ic1 == perturbations_module_->index_ic_ad_) && (index_ic2 == perturbations_module_->index_ic_ad_)) {

              sprintf(file_name,"%s%s",pop->root,"cls_ad.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar adiabatic (AD) mode");
            }

            if ((ppt->has_bi == _TRUE_) &&
                (index_ic1 == perturbations_module_->index_ic_bi_) && (index_ic2 == perturbations_module_->index_ic_bi_)) {

              sprintf(file_name,"%s%s",pop->root,"cls_bi.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar baryon isocurvature (BI) mode");
            }

            if ((ppt->has_cdi == _TRUE_) &&
                (index_ic1 == perturbations_module_->index_ic_cdi_) && (index_ic2 == perturbations_module_->index_ic_cdi_)) {

              sprintf(file_name,"%s%s",pop->root,"cls_cdi.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar CDM isocurvature (CDI) mode");
            }

            if ((ppt->has_nid == _TRUE_) &&
                (index_ic1 == perturbations_module_->index_ic_nid_) && (index_ic2 == perturbations_module_->index_ic_nid_)) {

              sprintf(file_name,"%s%s",pop->root,"cls_nid.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar neutrino density isocurvature (NID) mode");
            }

            if ((ppt->has_niv == _TRUE_) &&
                (index_ic1 == perturbations_module_->index_ic_niv_) && (index_ic2 == perturbations_module_->index_ic_niv_)) {

              sprintf(file_name,"%s%s",pop->root,"cls_niv.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar neutrino velocity isocurvature (NIV) mode");
            }

            if ((ppt->has_ad == _TRUE_) &&
                (ppt->has_bi == _TRUE_) && (index_ic1 == perturbations_module_->index_ic_ad_) && (index_ic2 == perturbations_module_->index_ic_bi_)) {

              sprintf(file_name,"%s%s",pop->root,"cls_ad_bi.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar cross ADxBI mode");
            }

            if ((ppt->has_ad == _TRUE_) && (ppt->has_cdi == _TRUE_) &&
                (index_ic1 == perturbations_module_->index_ic_ad_) && (index_ic2 == perturbations_module_->index_ic_cdi_)) {

              sprintf(file_name,"%s%s",pop->root,"cls_ad_cdi.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar cross ADxCDI mode");
            }

            if ((ppt->has_ad == _TRUE_) && (ppt->has_nid == _TRUE_) &&
                (index_ic1 == perturbations_module_->index_ic_ad_) && (index_ic2 == perturbations_module_->index_ic_nid_)) {

              sprintf(file_name,"%s%s",pop->root,"cls_ad_nid.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar cross ADxNID mode");
            }

            if ((ppt->has_ad == _TRUE_) && (ppt->has_niv == _TRUE_) &&
                (index_ic1 == perturbations_module_->index_ic_ad_) && (index_ic2 == perturbations_module_->index_ic_niv_)) {

              sprintf(file_name,"%s%s",pop->root,"cls_ad_niv.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar cross ADxNIV mode");
            }

            if ((ppt->has_bi == _TRUE_) && (ppt->has_cdi == _TRUE_) &&
                (index_ic1 == perturbations_module_->index_ic_bi_) && (index_ic2 == perturbations_module_->index_ic_cdi_)) {

              sprintf(file_name,"%s%s",pop->root,"cls_bi_cdi.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar cross BIxCDI mode");
            }

            if ((ppt->has_bi == _TRUE_) && (ppt->has_nid == _TRUE_) &&
                (index_ic1 == perturbations_module_->index_ic_bi_) && (index_ic2 == perturbations_module_->index_ic_nid_)) {

              sprintf(file_name,"%s%s",pop->root,"cls_bi_nid.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar cross BIxNID mode");
            }

            if ((ppt->has_bi == _TRUE_) && (ppt->has_niv == _TRUE_) &&
                (index_ic1 == perturbations_module_->index_ic_bi_) && (index_ic2 == perturbations_module_->index_ic_niv_)) {

              sprintf(file_name,"%s%s",pop->root,"cls_bi_niv.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar cross BIxNIV mode");
            }

            if ((ppt->has_cdi == _TRUE_) && (ppt->has_nid == _TRUE_) &&
                (index_ic1 == perturbations_module_->index_ic_cdi_) && (index_ic2 == perturbations_module_->index_ic_nid_)) {

              sprintf(file_name,"%s%s",pop->root,"cls_cdi_nid.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar cross CDIxNID mode");
            }

            if ((ppt->has_cdi == _TRUE_) && (ppt->has_niv == _TRUE_) &&
                (index_ic1 == perturbations_module_->index_ic_cdi_) && (index_ic2 == perturbations_module_->index_ic_niv_)) {

              sprintf(file_name,"%s%s",pop->root,"cls_cdi_niv.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar cross CDIxNIV mode");
            }

            if ((ppt->has_nid == _TRUE_) && (ppt->has_niv == _TRUE_) &&
                (index_ic1 == perturbations_module_->index_ic_nid_) && (index_ic2 == perturbations_module_->index_ic_niv_)) {

              sprintf(file_name,"%s%s",pop->root,"cls_nid_niv.dat");
              strcpy(first_line,"[l(l+1)/2pi] C_l's for scalar cross NIDxNIV mode");
            }

          }

          if (_tensorsEXT_) {

            class_test(0==0,
                       error_message_,
                       "Seems that we have mixed initial conditions for tensors? Should not happen!\n");

          }

          index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,spectra_module_->ic_size_[index_md]);

          if (spectra_module_->is_non_zero_[index_md][index_ic1_ic2] == _TRUE_) {

            class_call(output_open_cl_file(&(out_md_ic[index_md][index_ic1_ic2]),
                                           file_name,
                                           first_line,
                                           spectra_module_->l_max_[index_md]
                                           ),
                       error_message_,
                       error_message_);

          }
        }
      }

      class_alloc(cl_md_ic[index_md],
                  spectra_module_->ic_ic_size_[index_md]*spectra_module_->ct_size_*sizeof(double),
                  error_message_);
    }
  }

  /** - third, perform loop over l. For each multipole, get all \f$ C_l\f$'s
      by calling spectra_cl_at_l() and distribute the results to
      relevant files */

  for (l = 2; l <= spectra_module_->l_max_tot_; l++) {

    class_call(spectra_module_->spectra_cl_at_l((double)l, cl_tot, cl_md, cl_md_ic),
               psp->error_message,
               error_message_);

    class_call(output_one_line_of_cl(out, (double)l, cl_tot, spectra_module_->ct_size_),
               error_message_,
               error_message_);

    if ((ple->has_lensed_cls == _TRUE_) && (l <= lensing_module_->l_lensed_max_)) {

      class_call(lensing_module_->lensing_cl_at_l((double)l, cl_tot),
                 lensing_module_->error_message_,
                 error_message_);

      class_call(output_one_line_of_cl(out_lensed, l, cl_tot, spectra_module_->ct_size_),
                 error_message_,
                 error_message_);
    }

    if (perturbations_module_->md_size_ > 1) {
      for (index_md = 0; index_md < perturbations_module_->md_size_; index_md++) {
        if (l <= spectra_module_->l_max_[index_md]) {

          class_call(output_one_line_of_cl(out_md[index_md], l, cl_md[index_md], spectra_module_->ct_size_),
                     error_message_,
                     error_message_);
        }
      }
    }

    for (index_md = 0; index_md < perturbations_module_->md_size_; index_md++) {
      if ((perturbations_module_->ic_size_[index_md] > 1) && (l <= spectra_module_->l_max_[index_md])) {
        for (index_ic1_ic2 = 0; index_ic1_ic2 < spectra_module_->ic_ic_size_[index_md]; index_ic1_ic2++) {
          if (spectra_module_->is_non_zero_[index_md][index_ic1_ic2] == _TRUE_) {

            class_call(output_one_line_of_cl(out_md_ic[index_md][index_ic1_ic2], l, &(cl_md_ic[index_md][index_ic1_ic2*spectra_module_->ct_size_]), spectra_module_->ct_size_),
                       error_message_,
                       error_message_);
          }
        }
      }
    }
  }

  /** - finally, close files and free arrays of files and \f$ C_l\f$'s */

  for (index_md = 0; index_md < perturbations_module_->md_size_; index_md++) {
    if (perturbations_module_->ic_size_[index_md] > 1) {
      for (index_ic1_ic2 = 0; index_ic1_ic2 < spectra_module_->ic_ic_size_[index_md]; index_ic1_ic2++) {
        if (spectra_module_->is_non_zero_[index_md][index_ic1_ic2] == _TRUE_) {
          fclose(out_md_ic[index_md][index_ic1_ic2]);
        }
      }
      free(cl_md_ic[index_md]);
    }
  }
  if (perturbations_module_->md_size_ > 1) {
    for (index_md = 0; index_md < perturbations_module_->md_size_; index_md++) {
      fclose(out_md[index_md]);
      free(cl_md[index_md]);
    }
  }
  fclose(out);
  if (ple->has_lensed_cls == _TRUE_) {
    fclose(out_lensed);
  }
  free(cl_tot);
  for (index_md = 0; index_md < perturbations_module_->md_size_; index_md++) {
    free(out_md_ic[index_md]);
  }
  free(out_md_ic);
  free(cl_md_ic);
  free(out_md);
  free(cl_md);

  return _SUCCESS_;

}

/**
 * This routines writes the output in files for Fourier matter power spectra P(k)'s
 * (linear or non-linear)
 *
 * @param pba       Input: pointer to background structure (needed for calling spectra_pk_at_z())
 * @param ppt       Input: pointer perturbation structure
 * @param pnl       Input: pointer to nonlinear structure
 * @param pop       Input: pointer to output structure
 * @param pk_output Input: pk_linear or pk_nonlinear
 */

int OutputModule::output_pk(enum pk_outputs pk_output) {

  /** Summary: */

  /** - define local variables */

  FILE ** out_pk_ic = NULL;  /* out_pk_ic[index_ic1_ic2] is a pointer to a file with P(k) for each pair of ic */
  FILE * out_pk;             /* out_pk[index_pk] is a pointer to a file with total P(k) summed over ic */

  double * ln_pk_ic = NULL;  /* array ln_pk_ic[index_k * nonlinear_module_->ic_ic_size_ + index_ic1_ic2] */
  double * ln_pk;            /* array ln_pk[index_k] */

  int index_ic1,index_ic2;
  int index_ic1_ic2=0;
  int index_k;
  int index_z;
  int index_pk;

  FileName file_name;

  char redshift_suffix[7]; // 7 is enough to write "z%d_" as long as there are at most 10'000 bins
  char type_suffix[9];     // 6 is enough to write "pk_cb_nl" plus closing character \0
  char first_line[_LINE_LENGTH_MAX_];
  short do_ic = _FALSE_;

  /** - preliminary: check whether we need to output the decomposition into contributions from each initial condition */

  if ((pk_output == pk_linear) && (nonlinear_module_->ic_size_ > 1))
    do_ic = _TRUE_;

  /** - allocate arrays to store the P(k) */

  class_alloc(ln_pk,
              nonlinear_module_->k_size_*sizeof(double),
              error_message_);

  if (do_ic == _TRUE_) {

    class_alloc(ln_pk_ic,
                nonlinear_module_->k_size_*nonlinear_module_->ic_ic_size_*sizeof(double),
                error_message_);

    /** - allocate pointer to output files */

    class_alloc(out_pk_ic,
                nonlinear_module_->ic_ic_size_*sizeof(FILE *),
                error_message_);
  }

  /** - loop over pk type (_cb, _m) */

  for (index_pk = 0; index_pk < nonlinear_module_->pk_size_; index_pk++) {

    if ((nonlinear_module_->has_pk_m_ == _TRUE_) && (index_pk == nonlinear_module_->index_pk_m_)) {
      if (pk_output == pk_linear)
        sprintf(type_suffix,"pk");
      else
        sprintf(type_suffix,"pk_nl");
    }
    if ((nonlinear_module_->has_pk_cb_ == _TRUE_) && (index_pk == nonlinear_module_->index_pk_cb_)) {
      if (pk_output == pk_linear)
        sprintf(type_suffix,"pk_cb");
      else
        sprintf(type_suffix,"pk_cb_nl");
    }

    /** - loop over z */

    for (index_z = 0; index_z < pop->z_pk_num; index_z++) {

      /** - first, check that requested redshift z_pk is consistent */

      class_test((pop->z_pk[index_z] > ppt->z_max_pk),
                 error_message_,
                 "P(k,z) computed up to z=%f but requested at z=%f. Must increase z_max_pk in precision file.",ppt->z_max_pk,pop->z_pk[index_z]);

      if (pop->z_pk_num == 1)
        redshift_suffix[0]='\0';
      else
        sprintf(redshift_suffix,"z%d_",index_z+1);

      /** - second, open only the relevant files and write a header in each of them */

      sprintf(file_name,"%s%s%s%s",pop->root,redshift_suffix,type_suffix,".dat");

      class_call(output_open_pk_file(&out_pk,
                                     file_name,
                                     "",
                                     pop->z_pk[index_z]
                                     ),
                 error_message_,
                 error_message_);

      if (do_ic == _TRUE_) {

        for (index_ic1 = 0; index_ic1 < nonlinear_module_->ic_size_; index_ic1++) {

          for (index_ic2 = index_ic1; index_ic2 < nonlinear_module_->ic_size_; index_ic2++) {

            if ((ppt->has_ad == _TRUE_) && (index_ic1 == perturbations_module_->index_ic_ad_) && (index_ic2 == perturbations_module_->index_ic_ad_)) {
              sprintf(file_name,"%s%s%s%s",pop->root,redshift_suffix,type_suffix,"_ad.dat");
              strcpy(first_line,"for adiabatic (AD) mode ");
            }

            if ((ppt->has_bi == _TRUE_) && (index_ic1 == perturbations_module_->index_ic_bi_) && (index_ic2 == perturbations_module_->index_ic_bi_)) {
              sprintf(file_name,"%s%s%s%s",pop->root,redshift_suffix,type_suffix,"_bi.dat");
              strcpy(first_line,"for baryon isocurvature (BI) mode ");
            }

            if ((ppt->has_cdi == _TRUE_) && (index_ic1 == perturbations_module_->index_ic_cdi_) && (index_ic2 == perturbations_module_->index_ic_cdi_)) {
              sprintf(file_name,"%s%s%s%s",pop->root,redshift_suffix,type_suffix,"_cdi.dat");
              strcpy(first_line,"for CDM isocurvature (CDI) mode ");
            }

            if ((ppt->has_nid == _TRUE_) && (index_ic1 == perturbations_module_->index_ic_nid_) && (index_ic2 == perturbations_module_->index_ic_nid_)) {
              sprintf(file_name,"%s%s%s%s",pop->root,redshift_suffix,type_suffix,"_nid.dat");
              strcpy(first_line,"for neutrino density isocurvature (NID) mode ");
            }

            if ((ppt->has_niv == _TRUE_) && (index_ic1 == perturbations_module_->index_ic_niv_) && (index_ic2 == perturbations_module_->index_ic_niv_)) {
              sprintf(file_name,"%s%s%s%s",pop->root,redshift_suffix,type_suffix,"_niv.dat");
              strcpy(first_line,"for neutrino velocity isocurvature (NIV) mode ");
            }

            if ((ppt->has_ad == _TRUE_) && (ppt->has_bi == _TRUE_) && (index_ic1 == perturbations_module_->index_ic_ad_) && (index_ic2 == perturbations_module_->index_ic_bi_)) {
              sprintf(file_name,"%s%s%s%s",pop->root,redshift_suffix,type_suffix,"_ad_bi.dat");
              strcpy(first_line,"for cross ADxBI mode ");
            }

            if ((ppt->has_ad == _TRUE_) && (ppt->has_cdi == _TRUE_) && (index_ic1 == perturbations_module_->index_ic_ad_) && (index_ic2 == perturbations_module_->index_ic_cdi_)) {
              sprintf(file_name,"%s%s%s%s",pop->root,redshift_suffix,type_suffix,"_ad_cdi.dat");
              strcpy(first_line,"for cross ADxCDI mode ");
            }

            if ((ppt->has_ad == _TRUE_) && (ppt->has_nid == _TRUE_) && (index_ic1 == perturbations_module_->index_ic_ad_) && (index_ic2 == perturbations_module_->index_ic_nid_)) {
              sprintf(file_name,"%s%s%s%s",pop->root,redshift_suffix,type_suffix,"_ad_nid.dat");
              strcpy(first_line,"for scalar cross ADxNID mode ");
            }

            if ((ppt->has_ad == _TRUE_) && (ppt->has_niv == _TRUE_) && (index_ic1 == perturbations_module_->index_ic_ad_) && (index_ic2 == perturbations_module_->index_ic_niv_)) {
              sprintf(file_name,"%s%s%s%s",pop->root,redshift_suffix,type_suffix,"_ad_niv.dat");
              strcpy(first_line,"for cross ADxNIV mode ");
            }

            if ((ppt->has_bi == _TRUE_) && (ppt->has_cdi == _TRUE_) && (index_ic1 == perturbations_module_->index_ic_bi_) && (index_ic2 == perturbations_module_->index_ic_cdi_)) {
              sprintf(file_name,"%s%s%s%s",pop->root,redshift_suffix,type_suffix,"_bi_cdi.dat");
              strcpy(first_line,"for cross BIxCDI mode ");
            }

            if ((ppt->has_bi == _TRUE_) && (ppt->has_nid == _TRUE_) && (index_ic1 == perturbations_module_->index_ic_bi_) && (index_ic2 == perturbations_module_->index_ic_nid_)) {
              sprintf(file_name,"%s%s%s%s",pop->root,redshift_suffix,type_suffix,"_bi_nid.dat");
              strcpy(first_line,"for cross BIxNID mode ");
            }

            if ((ppt->has_bi == _TRUE_) && (ppt->has_niv == _TRUE_) && (index_ic1 == perturbations_module_->index_ic_bi_) && (index_ic2 == perturbations_module_->index_ic_niv_)) {
              sprintf(file_name,"%s%s%s%s",pop->root,redshift_suffix,type_suffix,"_bi_niv.dat");
              strcpy(first_line,"for cross BIxNIV mode ");
            }

            if ((ppt->has_cdi == _TRUE_) && (ppt->has_nid == _TRUE_) && (index_ic1 == perturbations_module_->index_ic_cdi_) && (index_ic2 == perturbations_module_->index_ic_nid_)) {
              sprintf(file_name,"%s%s%s%s",pop->root,redshift_suffix,type_suffix,"_cdi_nid.dat");
              strcpy(first_line,"for cross CDIxNID mode ");
            }

            if ((ppt->has_cdi == _TRUE_) && (ppt->has_niv == _TRUE_) && (index_ic1 == perturbations_module_->index_ic_cdi_) && (index_ic2 == perturbations_module_->index_ic_niv_)) {
              sprintf(file_name,"%s%s%s%s",pop->root,redshift_suffix,type_suffix,"_cdi_niv.dat");
              strcpy(first_line,"for cross CDIxNIV mode ");
            }

            if ((ppt->has_nid == _TRUE_) && (ppt->has_niv == _TRUE_) && (index_ic1 == perturbations_module_->index_ic_nid_) && (index_ic2 == perturbations_module_->index_ic_niv_)) {
              sprintf(file_name,"%s%s%s%s",pop->root,redshift_suffix,type_suffix,"_nid_niv.dat");
              strcpy(first_line,"for cross NIDxNIV mode ");
            }

            index_ic1_ic2 = index_symmetric_matrix(index_ic1, index_ic2, nonlinear_module_->ic_size_);

            if (nonlinear_module_->is_non_zero_[index_ic1_ic2] == _TRUE_) {

              class_call(output_open_pk_file(&(out_pk_ic[index_ic1_ic2]),
                                             file_name,
                                             first_line,
                                             pop->z_pk[index_z]
                                             ),
                         error_message_,
                         error_message_);
            }
          }
        }
      }

      /** - third, compute P(k) for each k */

      class_call(nonlinear_module_->nonlinear_pk_at_z(logarithmic, pk_output, pop->z_pk[index_z], index_pk, ln_pk, ln_pk_ic),
                 nonlinear_module_->error_message_,
                 error_message_);

      /** - fourth, write in files */

      for (index_k = 0; index_k < nonlinear_module_->k_size_; index_k++) {

        class_call(output_one_line_of_pk(out_pk,
                                         exp(nonlinear_module_->ln_k_[index_k])/pba->h,
                                         exp(ln_pk[index_k])*pow(pba->h,3)
                                         ),
                   error_message_,
                   error_message_);

        if (do_ic == _TRUE_) {

          for (index_ic1_ic2 = 0; index_ic1_ic2 < nonlinear_module_->ic_ic_size_; index_ic1_ic2++) {

            if (nonlinear_module_->is_non_zero_[index_ic1_ic2] == _TRUE_) {

              class_call(output_one_line_of_pk(out_pk_ic[index_ic1_ic2],
                                               exp(nonlinear_module_->ln_k_[index_k])/pba->h,
                                               exp(ln_pk_ic[index_k*nonlinear_module_->ic_ic_size_ + index_ic1_ic2])*pow(pba->h, 3)),
                         error_message_,
                         error_message_);
            }
          }
        }
      } /* end loop over k */

      /** - fifth, close files */

      fclose(out_pk);

      if (do_ic == _TRUE_) {
        for (index_ic1_ic2 = 0; index_ic1_ic2 < nonlinear_module_->ic_ic_size_; index_ic1_ic2++) {
          if (nonlinear_module_->is_non_zero_[index_ic1_ic2] == _TRUE_) {
            fclose(out_pk_ic[index_ic1_ic2]);
          }
        }
      }

    } /* end loop over index_z */

  } /* end loop over index_pk */

  /* free arrays and pointers */
  free(ln_pk);
  if (pk_output == pk_linear) {
    free(ln_pk_ic);
    free(out_pk_ic);
  }

  return _SUCCESS_;
}

/**
 * This routines writes the output in files for matter transfer functions \f$ T_i(k)\f$'s.
 *
 * @param pba Input: pointer to background structure (needed for calling spectra_pk_at_z())
 * @param ppt Input: pointer perturbation structure
 * @param pop Input: pointer to output structure
 */

int OutputModule::output_tk() {

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
  char redshift_suffix[7]; // 7 is enough to write "z%d_" as long as there are at most 10'000 bins
  char first_line[_LINE_LENGTH_MAX_];
  char ic_suffix[4];   // 4 is enough to write "ad", "bi", "cdi", "nid", "niv", ...


  index_md = perturbations_module_->index_md_scalars_;

  if (pop->output_format == camb_format) {

    class_test(pba->N_ncdm>1,
               error_message_,
               "you wish to output the transfer functions in CMBFAST/CAMB format but you have more than one non-cold dark matter (ncdm) species. The two are not compatible (since CMBFAST/CAMB only have one ncdm species): switch to CLASS output format or keep only on ncdm species");

    class_test(ppt->has_velocity_transfers == _TRUE_,
               error_message_,
               "you wish to output the transfer functions in CMBFAST/CAMB format, but you requested velocity transfer functions. The two are not compatible (since CMBFAST/CAMB do not compute velocity transfer functions): switch to CLASS output format, or ask only for density transfer function");
  }


  class_call(perturbations_module_->perturb_output_titles(pop->output_format, titles),
             perturbations_module_->error_message_,
             error_message_);
  number_of_titles = get_number_of_titles(titles);
  size_data = number_of_titles*perturbations_module_->k_size_[index_md];

  class_alloc(data, sizeof(double)*perturbations_module_->ic_size_[index_md]*size_data, error_message_);

  for (index_z = 0; index_z < pop->z_pk_num; index_z++) {

    z = pop->z_pk[index_z];

    /** - first, check that requested redshift z_pk is consistent */

    class_test((pop->z_pk[index_z] > ppt->z_max_pk),
               error_message_,
               "T_i(k,z) computed up to z=%f but requested at z=%f. Must increase z_max_pk in precision file.",ppt->z_max_pk,pop->z_pk[index_z]);

    if (pop->z_pk_num == 1)
      redshift_suffix[0]='\0';
    else
      sprintf(redshift_suffix,"z%d_",index_z+1);

    /** - second, open only the relevant files, and write a heading in each of them */

    class_call(perturbations_module_->perturb_output_data(pop->output_format, pop->z_pk[index_z], number_of_titles, data),
               perturbations_module_->error_message_,
               error_message_);

    for (index_ic = 0; index_ic < perturbations_module_->ic_size_[index_md]; index_ic++) {

      class_call(perturbations_module_->perturb_output_firstline_and_ic_suffix(index_ic, first_line, ic_suffix),
                 perturbations_module_->error_message_,
                 error_message_);

      if ((ppt->has_ad == _TRUE_) && (perturbations_module_->ic_size_[index_md] == 1))
        sprintf(file_name,"%s%s%s",pop->root,redshift_suffix,"tk.dat");
      else
        sprintf(file_name,"%s%s%s%s%s",pop->root,redshift_suffix,"tk_",ic_suffix,".dat");

      class_open(tkfile, file_name, "w", error_message_);

      if (pop->write_header == _TRUE_) {
        if (pop->output_format == class_format) {
          fprintf(tkfile,"# Transfer functions T_i(k) %sat redshift z=%g\n",first_line,z);
          fprintf(tkfile, "# for k=%g to %g h/Mpc,\n", perturbations_module_->k_[index_md][0]/pba->h, perturbations_module_->k_[index_md][perturbations_module_->k_size_[index_md] - 1]/pba->h);
          fprintf(tkfile, "# number of wavenumbers equal to %d\n", perturbations_module_->k_size_[index_md]);
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
          fprintf(tkfile, "# for k=%g to %g h/Mpc,\n", perturbations_module_->k_[index_md][0]/pba->h, perturbations_module_->k_[index_md][perturbations_module_->k_size_[index_md] - 1]/pba->h);
          fprintf(tkfile, "# number of wavenumbers equal to %d\n", perturbations_module_->k_size_[index_md]);
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

int OutputModule::output_background() {

  FILE * backfile;
  FileName file_name;

  char titles[_MAXTITLESTRINGLENGTH_]={0};
  double * data;
  int size_data, number_of_titles;

  class_call(background_module_->background_output_titles(titles),
             background_module_->error_message_,
             error_message_);
  number_of_titles = get_number_of_titles(titles);
  size_data = number_of_titles*background_module_->bt_size_;
  class_alloc(data, sizeof(double)*size_data, error_message_);
  class_call(background_module_->background_output_data(number_of_titles, data),
             background_module_->error_message_,
             error_message_);

  sprintf(file_name,"%s%s",pop->root,"background.dat");
  class_open(backfile, file_name, "w", error_message_);

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

int OutputModule::output_thermodynamics() {

  FileName file_name;
  FILE * thermofile;
  char titles[_MAXTITLESTRINGLENGTH_]={0};
  double * data;
  int size_data, number_of_titles;

  class_call(thermodynamics_module_->thermodynamics_output_titles(titles),
             thermodynamics_module_->error_message_,
             error_message_);
  number_of_titles = get_number_of_titles(titles);
  size_data = number_of_titles*thermodynamics_module_->tt_size_;
  class_alloc(data, sizeof(double)*size_data, error_message_);
  class_call(thermodynamics_module_->thermodynamics_output_data(number_of_titles, data),
             thermodynamics_module_->error_message_,
             error_message_);

  sprintf(file_name,"%s%s",pop->root,"thermodynamics.dat");
  class_open(thermofile, file_name, "w", error_message_);

  if (pop->write_header == _TRUE_) {
    fprintf(thermofile,"# Table of selected thermodynamics quantities\n");
    fprintf(thermofile,"# The following notation is used in column titles:\n");
    fprintf(thermofile,"#         x_e = electron ionization fraction\n");
    fprintf(thermofile,"#      -kappa = optical depth\n");
    fprintf(thermofile,"#      kappa' = Thomson scattering rate, prime denotes conformal time derivatives\n");
    fprintf(thermofile,"#           g = kappa' e^-kappa = visibility function \n");
    fprintf(thermofile,"#          Tb = baryon temperature \n");
    fprintf(thermofile,"#         w_b = baryon equation of state parameter \n");
    fprintf(thermofile,"#       c_b^2 = baryon sound speed squared \n");
    fprintf(thermofile,"#       tau_d = baryon drag optical depth \n");
    if (pth->compute_damping_scale == _TRUE_)
      fprintf(thermofile,"#         r_d = approximate comoving value of photon damping scale \n");
    if(pba->has_idm_dr == _TRUE_) {
      fprintf(thermofile,"#  dmu_idm_dr = scattering rate of idr with idm_dr (i.e. idr opacity to idm_dr scattering) (units 1/Mpc)\n");
      fprintf(thermofile,"# ddmu_idm_dr = derivative of this rate\n");
      fprintf(thermofile,"#  tau_idm_dr = optical depth of idm_dr (due to interactions with idr) \n");
      fprintf(thermofile,"#     tau_idr = optical depth of idr (due to self-interactions) \n");
      fprintf(thermofile,"#    g_idm_dr = visibility function of idm_idr \n");
      fprintf(thermofile,"#  c_idm_dr^2 = interacting dark matter squared sound speed \n");
      fprintf(thermofile,"#    T_idm_dr = temperature of DM interacting with DR \n");
      fprintf(thermofile,"#     dmu_idr = idr self-interaction rate \n");
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


int OutputModule::output_perturbations() {

  FILE * out;
  FileName file_name;
  int index_ikout, index_md;
  double k;

  for (index_ikout=0; index_ikout<ppt->k_output_values_num; index_ikout++){

    if (ppt->has_scalars == _TRUE_){
      index_md = perturbations_module_->index_md_scalars_;
      k = perturbations_module_->k_[index_md][perturbations_module_->index_k_output_values_[index_md*ppt->k_output_values_num + index_ikout]];
      sprintf(file_name,"%s%s%d%s",pop->root,"perturbations_k",index_ikout,"_s.dat");
      class_open(out, file_name, "w", error_message_);
      fprintf(out,"#scalar perturbations for mode k = %.*e Mpc^(-1)\n",_OUTPUTPRECISION_,k);
      output_print_data(out,
                        perturbations_module_->scalar_titles_,
                        perturbations_module_->scalar_perturbations_data_[index_ikout],
                        perturbations_module_->size_scalar_perturbation_data_[index_ikout]);

      fclose(out);
    }
    if (ppt->has_vectors == _TRUE_){
      index_md = perturbations_module_->index_md_vectors_;
      k = perturbations_module_->k_[index_md][perturbations_module_->index_k_output_values_[index_md*ppt->k_output_values_num + index_ikout]];
      sprintf(file_name,"%s%s%d%s",pop->root,"perturbations_k",index_ikout,"_v.dat");
      class_open(out, file_name, "w", error_message_);
      fprintf(out,"#vector perturbations for mode k = %.*e Mpc^(-1)\n",_OUTPUTPRECISION_,k);
      output_print_data(out,
                        perturbations_module_->vector_titles_,
                        perturbations_module_->vector_perturbations_data_[index_ikout],
                        perturbations_module_->size_vector_perturbation_data_[index_ikout]);

      fclose(out);
    }
    if (ppt->has_tensors == _TRUE_){
      index_md = perturbations_module_->index_md_tensors_;
      k = perturbations_module_->k_[index_md][perturbations_module_->index_k_output_values_[index_md*ppt->k_output_values_num + index_ikout]];
      sprintf(file_name,"%s%s%d%s",pop->root,"perturbations_k",index_ikout,"_t.dat");
      class_open(out, file_name, "w", error_message_);
      fprintf(out,"#tensor perturbations for mode k = %.*e Mpc^(-1)\n",_OUTPUTPRECISION_,k);
      output_print_data(out,
                        perturbations_module_->tensor_titles_,
                        perturbations_module_->tensor_perturbations_data_[index_ikout],
                        perturbations_module_->size_tensor_perturbation_data_[index_ikout]);

      fclose(out);
    }


  }
  return _SUCCESS_;

}

int OutputModule::output_primordial() {
  FileName file_name;
  FILE * out;
  char titles[_MAXTITLESTRINGLENGTH_]={0};
  double * data;
  int size_data, number_of_titles;

  sprintf(file_name,"%s%s",pop->root,"primordial_Pk.dat");

  class_call(primordial_module_->primordial_output_titles(titles),
             primordial_module_->error_message_,
             error_message_);
  number_of_titles = get_number_of_titles(titles);
  size_data = number_of_titles*primordial_module_->lnk_size_;
  class_alloc(data,sizeof(double)*size_data,error_message_);
  class_call(primordial_module_->primordial_output_data(number_of_titles, data),
             primordial_module_->error_message_,
             error_message_);

  class_open(out, file_name, "w", error_message_);
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


int OutputModule::output_print_data(FILE *out,
                      const char titles[_MAXTITLESTRINGLENGTH_],
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

int OutputModule::output_open_cl_file(
                        FILE * * clfile,
                        FileName filename,
                        char * first_line,
                        int lmax
                        ) {
  /** Summary */

  int index_d1,index_d2;
  int colnum = 1;
  char tmp[60]; //A fixed number here is ok, since it should just correspond to the largest string which is printed to tmp.

  class_open(*clfile, filename, "w", error_message_);

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

    if (spectra_module_->has_pp_ == _TRUE_) {
      if (pop->output_format == class_format) {
        fprintf(*clfile,"# -> for CMB lensing (phi), these are C_l^phi-phi for the lensing potential.\n");
      }
      if (pop->output_format == camb_format) {
        fprintf(*clfile,"# -> for CMB lensing (d), these are C_l^dd for the deflection field.\n");
      }
    }

    if (spectra_module_->has_ll_ == _TRUE_) {
      fprintf(*clfile,"# -> for galaxy lensing (lens[i]), these are C_l^phi-phi for the lensing potential.\n");
    }

    if (spectra_module_->has_pp_ == _TRUE_ || spectra_module_->has_ll_ == _TRUE_) {
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
      class_fprintf_columntitle(*clfile, "TT",     spectra_module_->has_tt_, colnum);
      class_fprintf_columntitle(*clfile, "EE",     spectra_module_->has_ee_, colnum);
      class_fprintf_columntitle(*clfile, "TE",     spectra_module_->has_te_, colnum);
      class_fprintf_columntitle(*clfile, "BB",     spectra_module_->has_bb_, colnum);
      class_fprintf_columntitle(*clfile, "phiphi", spectra_module_->has_pp_, colnum);
      class_fprintf_columntitle(*clfile, "TPhi",   spectra_module_->has_tp_, colnum);
      class_fprintf_columntitle(*clfile, "Ephi",   spectra_module_->has_ep_, colnum);
    }
    else if (pop->output_format == camb_format) {
      class_fprintf_columntitle(*clfile, "TT", spectra_module_->has_tt_, colnum);
      class_fprintf_columntitle(*clfile, "EE", spectra_module_->has_ee_, colnum);
      class_fprintf_columntitle(*clfile, "BB", spectra_module_->has_bb_, colnum);
      class_fprintf_columntitle(*clfile, "TE", spectra_module_->has_te_, colnum);
      class_fprintf_columntitle(*clfile, "dd", spectra_module_->has_pp_, colnum);
      class_fprintf_columntitle(*clfile, "dT", spectra_module_->has_tp_, colnum);
      class_fprintf_columntitle(*clfile, "dE", spectra_module_->has_ep_, colnum);
    }

    /** - Next deal with entries that are independent of format type */

    if (spectra_module_->has_dd_ == _TRUE_){
      for (index_d1=0; index_d1 < spectra_module_->d_size_; index_d1++){
        for (index_d2 = index_d1; index_d2 <= MIN(index_d1 + psp->non_diag,spectra_module_->d_size_ - 1); index_d2++){
          sprintf(tmp,"dens[%d]-dens[%d]",index_d1+1,index_d2+1);
          class_fprintf_columntitle(*clfile,tmp,_TRUE_,colnum);
        }
      }
    }
    if (spectra_module_->has_td_ == _TRUE_){
      for (index_d1 = 0; index_d1<spectra_module_->d_size_; index_d1++){
        sprintf(tmp,"T-dens[%d]",index_d1+1);
        class_fprintf_columntitle(*clfile,tmp,_TRUE_,colnum);
      }
    }
    if (spectra_module_->has_pd_ == _TRUE_){
      for (index_d1 = 0; index_d1<spectra_module_->d_size_; index_d1++){
        sprintf(tmp,"phi-dens[%d]",index_d1+1);
        class_fprintf_columntitle(*clfile,tmp,_TRUE_,colnum);
      }
    }
    if (spectra_module_->has_ll_ == _TRUE_){
      for (index_d1 = 0; index_d1 < spectra_module_->d_size_; index_d1++){
        for (index_d2 = index_d1; index_d2 <= MIN(index_d1 + psp->non_diag, spectra_module_->d_size_ - 1); index_d2++){
          sprintf(tmp,"lens[%d]-lens[%d]",index_d1+1,index_d2+1);
          class_fprintf_columntitle(*clfile,tmp,_TRUE_,colnum);
        }
      }
    }
    if (spectra_module_->has_tl_ == _TRUE_){
      for (index_d1 = 0; index_d1 < spectra_module_->d_size_; index_d1++){
        sprintf(tmp,"T-lens[%d]",index_d1+1);
        class_fprintf_columntitle(*clfile,tmp,_TRUE_,colnum);
      }
    }
    if (spectra_module_->has_dl_ == _TRUE_){
      for (index_d1 = 0; index_d1 < spectra_module_->d_size_; index_d1++){
        for (index_d2 = MAX(index_d1-psp->non_diag, 0); index_d2 <= MIN(index_d1 + psp->non_diag, spectra_module_->d_size_ - 1); index_d2++) {
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

int OutputModule::output_one_line_of_cl(
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
    class_fprintf_double(clfile, factor*pow(pba->T_cmb*1.e6,2)*cl[spectra_module_->index_ct_tt_], spectra_module_->has_tt_);
    class_fprintf_double(clfile, factor*pow(pba->T_cmb*1.e6,2)*cl[spectra_module_->index_ct_ee_], spectra_module_->has_ee_);
    class_fprintf_double(clfile, factor*pow(pba->T_cmb*1.e6,2)*cl[spectra_module_->index_ct_bb_], spectra_module_->has_bb_);
    class_fprintf_double(clfile, factor*pow(pba->T_cmb*1.e6,2)*cl[spectra_module_->index_ct_te_], spectra_module_->has_te_);
    class_fprintf_double(clfile, l*(l+1)*factor*cl[spectra_module_->index_ct_pp_], spectra_module_->has_pp_);
    class_fprintf_double(clfile, sqrt(l*(l+1))*factor*pba->T_cmb*1.e6*cl[spectra_module_->index_ct_tp_], spectra_module_->has_tp_);
    class_fprintf_double(clfile, sqrt(l*(l+1))*factor*pba->T_cmb*1.e6*cl[spectra_module_->index_ct_ep_], spectra_module_->has_ep_);
    index_ct_rest = 0;
    if (spectra_module_->has_tt_ == _TRUE_)
      index_ct_rest++;
    if (spectra_module_->has_ee_ == _TRUE_)
      index_ct_rest++;
    if (spectra_module_->has_bb_ == _TRUE_)
      index_ct_rest++;
    if (spectra_module_->has_te_ == _TRUE_)
      index_ct_rest++;
    if (spectra_module_->has_pp_ == _TRUE_)
      index_ct_rest++;
    if (spectra_module_->has_tp_ == _TRUE_)
      index_ct_rest++;
    if (spectra_module_->has_ep_ == _TRUE_)
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
 * @param pnl        Input: pointer to nonlinear structure
 * @param pop        Input: pointer to output structure
 * @param pkfile     Output: returned pointer to file pointer
 * @param filename   Input: name of the file
 * @param first_line Input: text describing the content (initial conditions, ...)
 * @param z          Input: redshift of the output
 * @return the error status
 */

int OutputModule::output_open_pk_file(
                        FILE * * pkfile,
                        FileName filename,
                        char * first_line,
                        double z
                        ) {

  int colnum = 1;
  class_open(*pkfile, filename, "w", error_message_);

  if (pop->write_header == _TRUE_) {
    fprintf(*pkfile, "# Matter power spectrum P(k) %sat redshift z=%g\n", first_line, z);
    fprintf(*pkfile, "# for k=%g to %g h/Mpc,\n",
            exp(nonlinear_module_->ln_k_[0])/pba->h,
            exp(nonlinear_module_->ln_k_[nonlinear_module_->k_size_ - 1])/pba->h);
    fprintf(*pkfile, "# number of wavenumbers equal to %d\n", nonlinear_module_->k_size_);

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

int OutputModule::output_one_line_of_pk(
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
