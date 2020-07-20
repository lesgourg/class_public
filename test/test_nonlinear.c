/** @file class.c
 * Julien Lesgourgues, 07.03.2014
 */

/* this main only runs the modules up to the nonlinear one */

#include "class.h"

int main(int argc, char **argv) {

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions */
  struct transfers tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct spectra sp;          /* for output spectra */
  struct nonlinear nl;        /* for non-linear spectra */
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

  /****** output the transfer functions ******/

  double z,k_nl,k_nl_cb,k;
  FILE * output;
  double * pvecback;
  int index_tau,index_k;
  int junk;
  double r_nl;

  if (nl.method == nl_halofit) {

    printf("Non-linear scale k_NL found by halofit:\n");

    z=0.;
    if (nonlinear_k_nl_at_z(&ba,&nl,z,&k_nl,&k_nl_cb) == _FAILURE_) {
      printf("\n\nError in nonlinear_k_nl_at_z \n=>%s\n",nl.error_message);
      return _FAILURE_;
    }
    printf("  z=%f   k_nl=%e\n",z,k_nl);

    z=0.5;
    if (nonlinear_k_nl_at_z(&ba,&nl,z,&k_nl,&k_nl_cb) == _FAILURE_) {
      printf("\n\nError in nonlinear_k_nl_at_z \n=>%s\n",nl.error_message);
      return _FAILURE_;
    }
    printf("  z=%f   k_nl=%e\n",z,k_nl);

    z=1.0;
    if (nonlinear_k_nl_at_z(&ba,&nl,z,&k_nl,&k_nl_cb) == _FAILURE_) {
      printf("\n\nError in nonlinear_k_nl_at_z \n=>%s\n",nl.error_message);
      return _FAILURE_;
    }
    printf("  z=%f   k_nl=%e\n",z,k_nl);

    printf("Non-linear correction factor r_nl=sqrt(P_nl/P_l) written in file with columns (z, k, r_nl\n");

    output=fopen("output/r_nl.dat","w");

    pvecback=malloc(ba.bg_size_short*sizeof(double));

    for (index_tau=0; index_tau<nl.tau_size; index_tau++) {

      if (background_at_tau(&ba,
                            nl.tau[index_tau],
                            ba.short_info,
                            ba.inter_normal,
                            &junk,
                            pvecback) == _FAILURE_) {
        printf("\n\nError in background_at_tau \n=>%s\n",ba.error_message);
        return _FAILURE_;
      }

      z=ba.a_today/pvecback[ba.index_bg_a]-1.;

      for (index_k=0; index_k<nl.k_size; index_k++) {

        k=nl.k[index_k];
        r_nl = nl.nl_corr_density[nl.index_pk_m][index_tau * nl.k_size + index_k];

        fprintf(output,"%e  %e  %e\n",z,k,r_nl);
      }
      fprintf(output,"\n\n");
    }

    fclose(output);
    free(pvecback);
  }

  /****** all calculations done, now free the structures ******/

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
