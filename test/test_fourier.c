/** @file class.c
 * Julien Lesgourgues, 07.03.2014
 */

/* this main only runs the modules up to the nonlinear one */

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

  if (primordial_init(&pr,&pt,&pm) == _FAILURE_) {
    printf("\n\nError in primordial_init \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }

  if (fourier_init(&pr,&ba,&th,&pt,&pm,&fo) == _FAILURE_) {
    printf("\n\nError in fourier_init \n=>%s\n",fo.error_message);
    return _FAILURE_;
  }

  /****** output the transfer functions ******/

  double z,k_nl,k_nl_cb,k;
  FILE * output;
  double * pvecback;
  int index_tau,index_k;
  int junk;
  double r_nl;

  if (fo.method == nl_halofit) {

    printf("Non-linear scale k_NL found by halofit:\n");

    z=0.;
    if (fourier_k_nl_at_z(&ba,&fo,z,&k_nl,&k_nl_cb) == _FAILURE_) {
      printf("\n\nError in fourier_k_nl_at_z \n=>%s\n",fo.error_message);
      return _FAILURE_;
    }
    printf("  z=%f   k_nl=%e\n",z,k_nl);

    z=0.5;
    if (fourier_k_nl_at_z(&ba,&fo,z,&k_nl,&k_nl_cb) == _FAILURE_) {
      printf("\n\nError in fourier_k_nl_at_z \n=>%s\n",fo.error_message);
      return _FAILURE_;
    }
    printf("  z=%f   k_nl=%e\n",z,k_nl);

    z=1.0;
    if (fourier_k_nl_at_z(&ba,&fo,z,&k_nl,&k_nl_cb) == _FAILURE_) {
      printf("\n\nError in fourier_k_nl_at_z \n=>%s\n",fo.error_message);
      return _FAILURE_;
    }
    printf("  z=%f   k_nl=%e\n",z,k_nl);

    printf("Non-linear correction factor r_nl=sqrt(P_nl/P_l) written in file with columns (z, k, r_nl\n");

    output=fopen("output/r_fo.dat","w");

    pvecback=malloc(ba.bg_size_short*sizeof(double));

    for (index_tau=0; index_tau<fo.tau_size; index_tau++) {

      if (background_at_tau(&ba,
                            fo.tau[index_tau],
                            short_info,
                            inter_normal,
                            &junk,
                            pvecback) == _FAILURE_) {
        printf("\n\nError in background_at_tau \n=>%s\n",ba.error_message);
        return _FAILURE_;
      }

      z=1./pvecback[ba.index_bg_a]-1.;

      for (index_k=0; index_k<fo.k_size; index_k++) {

        k=fo.k[index_k];
        r_nl = fo.nl_corr_density[fo.index_pk_m][index_tau * fo.k_size + index_k];

        fprintf(output,"%e  %e  %e\n",z,k,r_nl);
      }
      fprintf(output,"\n\n");
    }

    fclose(output);
    free(pvecback);
  }

  /****** all calculations done, now free the structures ******/

  if (fourier_free(&fo) == _FAILURE_) {
    printf("\n\nError in fourier_free \n=>%s\n",fo.error_message);
    return _FAILURE_;
  }

  if (primordial_free(&pm) == _FAILURE_) {
    printf("\n\nError in primordial_free \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }

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
