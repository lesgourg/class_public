/** @file cl.c Documented nonlinear module
 *
 * Benjamin Audren and Julien Lesgourgues, 21.12.2010    
 *
 */

#include "nonlinear.h"

int nonlinear_init(
		   struct precision *ppr,
		   struct background *pba,
		   struct thermo *pth,
		   struct primordial *ppm,
		   struct spectra *psp,
		   struct nonlinear *pnl
		   ) {

  int index_eta,index_k;

  class_test((pnl->method < nl_none) || (pnl->method > nl_trg),
	     pnl->error_message,
	     "Your non-linear method variable is set to %d, our of the range defined in nonlinear.h",pnl->method);

  if (pnl->method == nl_none) {
    if (pnl->nonlinear_verbose > 0)
      printf("No non-linear spectra requested. Nonlinear module skipped.\n");
    return _SUCCESS_;
  }
  else {
    if (pnl->nonlinear_verbose > 0)
      printf("Compute non-linear matter power spectrum using trg module.\n");

    struct spectra_nl trg;

    if (pnl->method == nl_trg_linear)
      trg.mode = 0;
    if (pnl->method == nl_trg_one_loop)
      trg.mode = 1;
    if (pnl->method == nl_trg)
      trg.mode = 2;

    trg.k_max = exp(psp->ln_k[psp->ln_k_size-1]) * pba->h - 1.;

    trg.spectra_nl_verbose = pnl->nonlinear_verbose;

    class_call(trg_init(ppr,pba,pth,ppm,psp,&trg),
	       trg.error_message,
	       pnl->error_message);

    /* copy non-linear spectrum in pnl */

    pnl->eta_size = trg.eta_size;
    pnl->k_size = trg.k_size;

    class_calloc(pnl->p_dd,
		 pnl->k_size*pnl->eta_size,
		 sizeof(double),
		 pnl->error_message);
    class_calloc(pnl->p_dt,
		 pnl->k_size*pnl->eta_size,
		 sizeof(double),
		 pnl->error_message);
    class_calloc(pnl->p_tt,
		 pnl->k_size*pnl->eta_size,
		 sizeof(double),
		 pnl->error_message);

    for (index_eta=0; index_eta<pnl->eta_size; index_eta++) {
      for (index_k=0; index_k<pnl->k_size; index_k++) {

	pnl->p_dd[index_eta*pnl->k_size+index_k]=trg.p_11_nl[index_eta*pnl->k_size+index_k];
	pnl->p_dt[index_eta*pnl->k_size+index_k]=trg.p_12_nl[index_eta*pnl->k_size+index_k];
	pnl->p_tt[index_eta*pnl->k_size+index_k]=trg.p_22_nl[index_eta*pnl->k_size+index_k];

      }
    }

    class_call(trg_free(&trg),
	       trg.error_message,
	       pnl->error_message);
  }
  
  return _SUCCESS_;
}

int nonlinear_free(
		   struct nonlinear *pnl
		   ) {

  if (pnl->method > nl_none) {
    free(pnl->p_dd);
    free(pnl->p_dt);
    free(pnl->p_tt);
  }

  return _SUCCESS_;

}
