/** @file class.c
 * Julien Lesgourgues, 17.04.2011
 */

/* this main only runs the modules up to the transfer one */

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

  if (transfer_init(&pr,&ba,&th,&pt,&fo,&tr) == _FAILURE_) {
    printf("\n\nError in transfer_init \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  /****** output the transfer functions ******/

  printf("Output of transfer functions (l, q, k, nu, Delta)\n");
  printf("(in flat space, q=k and nu=inf) \n");

  /* 1) select the mode, initial condition, type and multipole of the
     function you want to plot: */

  int index_mode=pt.index_md_scalars;
  int index_ic  =pt.index_ic_ad;
  int index_type=tr.index_tt_t0;

  /* 2) here is an illustration of how to output the transfer
     functions at some (k,l)'s of your choice */

  /*
  int index_l = 0;
  double q=3.6e-4;
  double transfer;

  if (transfer_functions_at_q(&tr,
                              index_mode,
                              index_ic,
                              index_type,
                              index_l,
                              q,
                              &transfer
                              ) == _FAILURE_) {
    printf("\n\nError in transfer_function_at_k \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  printf("%d %e %e\n",tr.l[index_l],q,transfer);
  */

  /* 3) here you can output the full tabulated arrays for all k and l's*/

  int index_q;
  int index_l;
  double transfer;
  FILE * output;

  output=fopen("output/test.trsf","w");

  for (index_l=0; index_l<tr.l_size[index_mode]; index_l++) {
    for (index_q=0; index_q<tr.q_size; index_q++) {

          /* use this to plot a single type : */

          transfer = tr.transfer[index_mode]
            [((index_ic * tr.tt_size[index_mode] + index_type)
              * tr.l_size[index_mode] + index_l)
             * tr.q_size + index_q];

          /* or use this to plot the full temperature transfer function: */
          /*
          transfer =
            tr.transfer[index_mode][((index_ic * tr.tt_size[index_mode] + tr.index_tt_t0) * tr.l_size[index_mode] + index_l) * tr.q_size + index_q] +
            tr.transfer[index_mode][((index_ic * tr.tt_size[index_mode] + tr.index_tt_t1) * tr.l_size[index_mode] + index_l) * tr.q_size + index_q] +
            tr.transfer[index_mode][((index_ic * tr.tt_size[index_mode] + tr.index_tt_t2) * tr.l_size[index_mode] + index_l) * tr.q_size + index_q];
          */

          if (transfer != 0.) {
            fprintf(output,"%d %e %e %e %e\n",
                    tr.l[index_l],
                    tr.q[index_q],
                    tr.k[index_mode][index_q],
                    tr.q[index_q]/sqrt(ba.sgnK*ba.K),
                    transfer);
          }
        }

        fprintf(output,"\n\n");
        //}
  }

  fclose(output);

  /****** all calculations done, now free the structures ******/

  if (transfer_free(&tr) == _FAILURE_) {
    printf("\n\nError in transfer_free \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

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
