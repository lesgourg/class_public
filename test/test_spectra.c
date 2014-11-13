/** @file test_spectra.c
 * 
 * Julien Lesgourgues, 26.08.2010
 *
 * main intended for computing power spectra, not using the output module.
 *     
 */
 
#include "class.h"

main(int argc, char **argv) {

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions */
  struct bessels bs;          /* for bessel functions */
  struct transfers tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct spectra sp;          /* for output spectra */
  struct output op;           /* for output files */
  struct spectra_nl nl;       /* for calculation of non-linear spectra */

  ErrorMsg errmsg;

  if (input_init_from_arguments(argc, argv,&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&op,&nl,errmsg) == _FAILURE_) {
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

  if (bessel_init(&pr,&bs) == _FAILURE_) {
    printf("\n\nError in bessel_init \n =>%s\n",bs.error_message);
    return _FAILURE_;
  }

  if (transfer_init(&pr,&ba,&th,&pt,&bs,&tr) == _FAILURE_) {
    printf("\n\nError in transfer_init \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  if (primordial_init(&pr,&pt,&pm) == _FAILURE_) {
    printf("\n\nError in primordial_init \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }

  if (spectra_init(&ba,&pt,&tr,&pm,&sp) == _FAILURE_) {
    printf("\n\nError in spectra_init \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }

  FILE * output;
  int index_mode=0;

  /****** output Cls ******/

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

  /****** output P(k) ******/

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

  /****** output T_i(k) ******/

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

  

  /****************************/

  if (primordial_free(&pm) == _FAILURE_) {
    printf("\n\nError in primordial_free \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }

  if (transfer_free(&tr) == _FAILURE_) {
    printf("\n\nError in transfer_free \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  if (bessel_free(&bs) == _FAILURE_)  {
    printf("\n\nError in bessel_free \n=>%s\n",bs.error_message);
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

}
