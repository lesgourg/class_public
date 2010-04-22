/** @file class.c 
 * Julien Lesgourgues, 18.04.2010    
 */
 
#include "precision.h"
#include "background.h"
#include "thermodynamics.h"
#include "perturbations.h"
#include "bessel.h"
#include "transfer.h"
#include "primordial.h"
#include "spectra.h"

main() {

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions */
  struct bessels bs;          /* for bessel functions */
  struct transfers tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct spectra sp;          /* for output spectra */
  struct spectra op;          /* for output files */
 
  if (precision_init(&pr) == _FAILURE_) {
    printf("\n\nError running precision_init \n=>%s\n",pr.error_message); 
    return _FAILURE_;
  }

  if (input_init(&ba,&th,&pt,&bs,&tr,&pm,&sp,&op) == _FAILURE_) {
    printf("\n\nError running input_init"); 
    return _FAILURE_;
  }

  if (background_init(&pr,&ba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_init(&ba,&pr,&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (perturb_init(&ba,&th,&pr,&pt) == _FAILURE_) {
    printf("\n\nError in perturb_init \n=>%s\n",pt.error_message);
    return _FAILURE_;
  }

  /****** here you can output the source functions ******/

  FILE * output;
  int index_k,index_eta;
  int index_mode=pt.index_md_scalars;
  int index_type=pt.index_tp_l;
  int index_ic=pt.index_ic_ad;

  output=fopen("output/source.dat","w");

  for (index_k=0; index_k < pt.k_size[index_mode]; index_k++) {
    for (index_eta=0; index_eta < pt.eta_size; index_eta++) { 

      fprintf(output,"%e %e %e\n",
	      pt.eta_sampling[index_eta],
	      pt.k[index_mode][index_k],
	      pt.sources[index_mode]
	      [index_ic * pt.tp_size + index_type]
	      [index_eta * pt.k_size[index_mode] + index_k]
	      );
    }
    fprintf(output,"\n");
  }

  fclose(output);

  /******************************************************/

  if (perturb_free() == _FAILURE_) {
    printf("\n\nError in perturb_free \n=>%s\n",pt.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_free() == _FAILURE_) {
    printf("\n\nError in thermodynamics_free \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (background_free() == _FAILURE_) {
    printf("\n\nError in background_free \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

}
