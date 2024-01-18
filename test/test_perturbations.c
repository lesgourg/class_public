/** @file class.c
 * Julien Lesgourgues, 17.04.2011
 */

/* this main runs only the background, thermodynamics and perturbation part */

#include "class.h"

int main(int argc, char **argv) {

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermodynamics th;           /* for thermodynamics */
  struct perturbations pt;         /* for source functions */
  struct transfer tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct harmonic hr;          /* for output spectra */
  struct fourier fo;        /* for non-linear spectra */
  struct lensing le;          /* for lensed spectra */
  struct distortions sd;      /* for spectral distortions */
  struct output op;           /* for output files */
  ErrorMsg errmsg;            /* for error messages */

  if (input_init(argc, argv,&pr,&ba,&th,&pt,&tr,&pm,&hr,&fo,&le,&sd,&op,errmsg) == _FAILURE_) {
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

  if (perturbations_init(&pr,&ba,&th,&pt) == _FAILURE_) {
    printf("\n\nError in perturbations_init \n=>%s\n",pt.error_message);
    return _FAILURE_;
  }

  if (pt.has_perturbations == _TRUE_) {

    /*********************************************************************/
    /*  here you can output the source function S(k,tau) of your choice  */
    /*********************************************************************/

    FILE * output;
    int index_k,index_tau;

    /* choose a mode (scalar, tensor, ...) */
    int index_md=pt.index_md_scalars;

    /* choose a type (temperature, polarization, grav. pot., ...) */
    int index_type=pt.index_tp_t0;

    /* choose an initial condition (ad, bi, cdi, nid, niv, ...) */
    int index_ic=pt.index_ic_ad;

    output=fopen("output/source.dat","w");
    fprintf(output,"#   k       tau       S\n");

    for (index_k=0; index_k < pt.k_size[index_md]; index_k++) {
      for (index_tau=0; index_tau < pt.tau_size; index_tau++) {

        fprintf(output,"%e %e %e\n",
                pt.k[index_md][index_k],
                pt.tau_sampling[index_tau],
                pt.sources[index_md]
                [index_ic * pt.tp_size[index_md] + index_type]
                [index_tau * pt.k_size[index_md] + index_k]
                );
      }
      fprintf(output,"\n");
    }

    fclose(output);

  }

  /****** all calculations done, now free the structures ******/

  if (perturbations_free(&pt) == _FAILURE_) {
    printf("\n\nError in perturbations_free \n=>%s\n",pt.error_message);
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
