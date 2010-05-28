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
#include "output.h"

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

  if (thermodynamics_init(&pr,&ba,&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  /********************************************/
  /***** output thermodynamics quantities *****/
  /********************************************/
  
  int i;

  printf("#1: redshift z\n");
  printf("#2: electron ionization fraction x_e\n");
  printf("#3: exponential of optical depth e^-kappa\n");
  printf("#4: Thomson scattering rate kappa'\n");
  printf("#5: visibility function g=kappa' e^-kappa \n");
  for (i=0; i < th.tt_size; i++)
    printf("%e %e %e %e %e\n",
	   th.z_table[i],
	   th.thermodynamics_table[i*th.th_size+th.index_th_xe],
	   th.thermodynamics_table[i*th.th_size+th.index_th_exp_m_kappa],
	   th.thermodynamics_table[i*th.th_size+th.index_th_dkappa],
	   th.thermodynamics_table[i*th.th_size+th.index_th_g]);

  /********************************************/

  if (thermodynamics_free(&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_free \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (background_free(&ba) == _FAILURE_) {
    printf("\n\nError in background_free \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

}
