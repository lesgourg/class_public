/** @file test_spectra.c
 *
 * Julien Lesgourgues, 26.08.2010
 *
 * main intended for computing power spectra, not using the output module.
 *
 */

#include "class.h"

int main(int argc, char **argv) {

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions */
  struct primordial pm;       /* for primordial spectra */
  struct nonlinear nl;        /* for non-linear spectra */
  struct transfers tr;        /* for transfer functions */
  struct spectra sp;          /* for output spectra */
  struct lensing le;          /* for lensed spectra */
  struct output op;           /* for output files */
  ErrorMsg errmsg;            /* for error messages */

  if (input_init_from_arguments(argc, argv,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg);
    return _FAILURE_;
  }

  if (background_init(&pr,&ba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_init(&pr,&ba,&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (perturb_init(&pr,&ba,&th,&pt) == _FAILURE_) {
    printf("\n\nError in perturb_init \n=>%s\n",pt.error_message);
    return _FAILURE_;
  }

  if (primordial_init(&pr,&pt,&pm) == _FAILURE_) {
    printf("\n\nError in primordial_init \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }

  if (nonlinear_init(&pr,&ba,&th,&pt,&pm,&nl) == _FAILURE_) {
    printf("\n\nError in nonlinear_init \n=>%s\n",nl.error_message);
    return _FAILURE_;
  }

  if (transfer_init(&pr,&ba,&th,&pt,&nl,&tr) == _FAILURE_) {
    printf("\n\nError in transfer_init \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  if (spectra_init(&pr,&ba,&pt,&pm,&nl,&tr,&sp) == _FAILURE_) {
    printf("\n\nError in spectra_init \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }

  /************************/

  //double z=2.;
  double z=0.;

  double * out_pk_l;
  double * out_pk_ic_l;
  int index_k;
  int index_ic_ic;

  /*
  out_pk_l = calloc(nl.k_size,sizeof(double));
  out_pk_ic_l = calloc(nl.k_size*nl.ic_ic_size,sizeof(double));

  if(nonlinear_pk_linear_at_z(&ba,
                              &nl,
                              linear,
                              z,
                              0,
                              out_pk_l,
                              out_pk_ic_l) == _FAILURE_) {
    printf("\n\nError in spectra_init \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }

  for (index_k=0; index_k<nl.k_size; index_k++) {
    fprintf(stdout,"%e   ",out_pk_l[index_k]);
    if (nl.ic_size > 1) {
      for (index_ic_ic = 0; index_ic_ic < nl.ic_ic_size; index_ic_ic++) {
        fprintf(stdout,"%e   ",out_pk_ic_l[index_k*nl.ic_ic_size+index_ic_ic]);
      }
    }
    fprintf(stdout,"\n");
  }
  */

  /************************/

  double k=0.02;
  //double k=10.;

  double * pk;
  double * pk_ic;
  double * pk_cb;
  double * pk_cb_ic;

  pk=calloc(1,sizeof(double));
  pk_cb=calloc(1,sizeof(double));
  if (nl.ic_size > 1) {
    pk_ic=calloc(nl.ic_ic_size,sizeof(double));
    pk_cb_ic=calloc(nl.ic_ic_size,sizeof(double));
  }

  if (spectra_pk_at_k_and_z_new(&ba,
                                &pm,
                                &sp,
                                k,
                                z,
                                pk,
                                pk_ic,
                                pk_cb,
                                pk_cb_ic) == _FAILURE_) {
    printf("\n\nError in spectra_init \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }

  if (nl.has_pk_m == _TRUE_) fprintf(stdout,"old P_m =%e\n",*pk);
  if (nl.has_pk_cb == _TRUE_) fprintf(stdout,"old P_cb=%e\n",*pk_cb);
  if (nl.ic_size > 1) {
    for (index_ic_ic = 0; index_ic_ic < nl.ic_ic_size; index_ic_ic++) {
      if (nl.has_pk_m == _TRUE_) fprintf(stdout,"old P_m [ic_ic=%d]=%e\n",index_ic_ic,pk_ic[index_ic_ic]);
      if (nl.has_pk_cb == _TRUE_) fprintf(stdout,"old P_cb[ic_ic=%d]=%e\n",index_ic_ic,pk_cb_ic[index_ic_ic]);
    }
  }

  if (nl.has_pk_m == _TRUE_) *pk=0.;
  if (nl.has_pk_cb == _TRUE_) *pk_cb=0.;
  if (nl.ic_size > 1) {
    for (index_ic_ic = 0; index_ic_ic < nl.ic_ic_size; index_ic_ic++) {
      if (nl.has_pk_m == _TRUE_) pk_ic[index_ic_ic]=0.;
      if (nl.has_pk_cb == _TRUE_) pk_cb_ic[index_ic_ic]=0.;
    }
  }

  if (spectra_pk_at_k_and_z_new(&ba,
                                &pm,
                                &sp,
                                k,
                                z,
                                pk,
                                pk_ic,
                                pk_cb,
                                pk_cb_ic) == _FAILURE_) {
    printf("\n\nError in spectra_init \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }

  if (nl.has_pk_m == _TRUE_) fprintf(stdout,"new P_m =%e\n",*pk);
  if (nl.has_pk_cb == _TRUE_) fprintf(stdout,"new P_cb=%e\n",*pk_cb);
  if (nl.ic_size > 1) {
    for (index_ic_ic = 0; index_ic_ic < nl.ic_ic_size; index_ic_ic++) {
      if (nl.has_pk_m == _TRUE_) fprintf(stdout,"new P_m [ic_ic=%d]=%e\n",index_ic_ic,pk_ic[index_ic_ic]);
      if (nl.has_pk_cb == _TRUE_) fprintf(stdout,"new P_cb[ic_ic=%d]=%e\n",index_ic_ic,pk_cb_ic[index_ic_ic]);
    }
  }


  //FILE * output;
  //int index_mode=0;

  /****** output Cls ******/

  /*
  int index_ic1_ic2=0;
  int index_ct=0;
  int index_l;

  if (pt.has_cmb == _TRUE_) {

    output=fopen("output/testing_cls.dat","w");

    for (index_l=0; index_l < sp.l_size[index_mode]; index_l++)
      fprintf(output,"%g %g\n",
	      sp.l[index_mode][index_l],
	      sp.cl[index_mode][(index_l * sp.ic_ic_size[index_mode] + index_ic1_ic2) * sp.ct_size + index_ct]);

    fclose(output);

  }
  */

  /****** output P(k) ******/

  /*
  int index_eta = sp.ln_eta_size-1;
  int index_k;
  double pk;
  double junk;

  if (pt.has_pk_matter == _TRUE_) {

    output=fopen("output/testing_pks.dat","w");

    for (index_k=0; index_k < sp.ln_k_size; index_k++)
      fprintf(output,"%g %g\n",
	      sp.ln_k[index_k],
	      sp.ln_pk[(index_eta * sp.ln_k_size + index_k) * sp.ic_ic_size[index_mode] + index_ic1_ic2]);

    fclose(output);

  }
  */

  /****** output T_i(k) ******/
  /*
  int index_ic=0;
  int index_tr;
  double * tk;
  double * tkk;

  if (pt.has_matter_transfers == _TRUE_) {

    output=fopen("output/testing_tks.dat","w");

    for (index_k=0; index_k < sp.ln_k_size; index_k++) {
      fprintf(output,"%g",sp.ln_k[index_k]);
      for (index_tr=0; index_tr < sp.tr_size; index_tr++) {
	fprintf(output,"  %g",
		sp.matter_transfer[((index_eta * sp.ln_k_size + index_k) * sp.ic_size[index_mode] + index_ic) * sp.tr_size + index_tr]);
      }
      fprintf(output,"\n");
    }
  }

  */

  /****************************/

  if (spectra_free(&sp) == _FAILURE_) {
    printf("\n\nError in spectra_free \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }

  if (transfer_free(&tr) == _FAILURE_) {
    printf("\n\nError in transfer_free \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  if (nonlinear_free(&nl) == _FAILURE_) {
    printf("\n\nError in nonlinear_free \n=>%s\n",nl.error_message);
    return _FAILURE_;
  }

  if (primordial_free(&pm) == _FAILURE_) {
    printf("\n\nError in primordial_free \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }

  if (perturb_free(&pt) == _FAILURE_) {
    printf("\n\nError in perturb_free \n=>%s\n",pt.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_free(&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_free \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (background_free(&ba) == _FAILURE_) {
    printf("\n\nError in background_free \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  return _SUCCESS_;

}
