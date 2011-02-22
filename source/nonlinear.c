/** @file cl.c Documented nonlinear module
 *
 * Benjamin Audren and Julien Lesgourgues, 21.12.2010    
 *
 */

#include "nonlinear.h"

int nonlinear_pk_at_z(
		      struct nonlinear * pnl,
		      double z,
		      double * pz_density,
		      double * pz_velocity,
		      double * pz_cross
		      ) {

  int last_index;

  class_call(array_interpolate_spline(pnl->z,
				      pnl->z_size,
				      pnl->p_density,
				      pnl->ddp_density,
				      pnl->k_size,
				      z,
				      &last_index,
				      pz_density,
				      pnl->k_size,
				      pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);

  class_call(array_interpolate_spline(pnl->z,
				      pnl->z_size,
				      pnl->p_velocity,
				      pnl->ddp_velocity,
				      pnl->k_size,
				      z,
				      &last_index,
				      pz_velocity,
				      pnl->k_size,
				      pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);

  class_call(array_interpolate_spline(pnl->z,
				      pnl->z_size,
				      pnl->p_cross,
				      pnl->ddp_cross,
				      pnl->k_size,
				      z,
				      &last_index,
				      pz_cross,
				      pnl->k_size,
				      pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);

  return _SUCCESS_;
}

int nonlinear_pk_at_k_and_z(
			    struct nonlinear * pnl,
			    double k,
			    double z,
			    double * pk_density,
			    double * pk_velocity,
			    double * pk_cross
			    ) {
  
  double * pz_density;
  double * pz_velocity;
  double * pz_cross;
  double * ddpz_density;
  double * ddpz_velocity;
  double * ddpz_cross;
  int last_index;

  class_alloc(pz_density,pnl->k_size*sizeof(double),pnl->error_message);
  class_alloc(pz_velocity,pnl->k_size*sizeof(double),pnl->error_message);
  class_alloc(pz_cross,pnl->k_size*sizeof(double),pnl->error_message);
  class_alloc(ddpz_density,pnl->k_size*sizeof(double),pnl->error_message);
  class_alloc(ddpz_velocity,pnl->k_size*sizeof(double),pnl->error_message);
  class_alloc(ddpz_cross,pnl->k_size*sizeof(double),pnl->error_message);

  class_call(nonlinear_pk_at_z(pnl,z,pz_density,pz_velocity,pz_cross),
	     pnl->error_message,
	     pnl->error_message);

  class_call(array_spline_table_lines(pnl->k,
				      pnl->k_size,
				      pz_density,
				      1,
				      ddpz_density,
				      _SPLINE_NATURAL_,
				      pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);
      
  class_call(array_interpolate_spline(pnl->k,
				      pnl->k_size,
				      pz_density,
				      ddpz_density,
				      1,
				      k,
				      &last_index,
				      pk_density,
				      1,
				      pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);

  class_call(array_spline_table_lines(pnl->k,
				      pnl->k_size,
				      pz_velocity,
				      1,
				      ddpz_velocity,
				      _SPLINE_NATURAL_,
				      pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);
      
  class_call(array_interpolate_spline(pnl->k,
				      pnl->k_size,
				      pz_velocity,
				      ddpz_velocity,
				      1,
				      k,
				      &last_index,
				      pk_velocity,
				      1,
				      pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);

  class_call(array_spline_table_lines(pnl->k,
				      pnl->k_size,
				      pz_cross,
				      1,
				      ddpz_cross,
				      _SPLINE_NATURAL_,
				      pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);
      
  class_call(array_interpolate_spline(pnl->k,
				      pnl->k_size,
				      pz_cross,
				      ddpz_cross,
				      1,
				      k,
				      &last_index,
				      pk_cross,
				      1,
				      pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);

  free(pz_density);
  free(pz_velocity);
  free(pz_cross);
  free(ddpz_density);
  free(ddpz_velocity);
  free(ddpz_cross);

  return _SUCCESS_;
}

int nonlinear_init(
		   struct precision *ppr,
		   struct background *pba,
		   struct thermo *pth,
		   struct primordial *ppm,
		   struct spectra *psp,
		   struct nonlinear *pnl
		   ) {

  int index_z,index_k;

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
      printf("Computing non-linear matter power spectrum using trg module\n");

    struct spectra_nl trg;

    if (pnl->method == nl_trg_linear)
      trg.mode = 0;
    if (pnl->method == nl_trg_one_loop)
      trg.mode = 1;
    if (pnl->method == nl_trg)
      trg.mode = 2;

    trg.k_max = exp(psp->ln_k[psp->ln_k_size-1]) * pba->h - 1.;

    trg.double_escape = ppr->double_escape;
    trg.has_bc_spectrum = ppr->has_bc_spectrum;
    trg.z_ini = ppr->z_ini;
    trg.eta_size = ppr->eta_size;
    trg.k_L = ppr->k_L;
    trg.k_min = ppr->k_min;
    trg.logstepx_min = ppr->logstepx_min;

    trg.spectra_nl_verbose = pnl->nonlinear_verbose;

    class_call(trg_init(ppr,pba,pth,ppm,psp,&trg),
	       trg.error_message,
	       pnl->error_message);

      fprintf(stderr," -> done with trg_init\n");

    /* copy non-linear spectrum in pnl */

    pnl->z_size = trg.eta_size;
    pnl->k_size = trg.k_size;

    class_calloc(pnl->k,
		 pnl->k_size,
		 sizeof(double),
		 pnl->error_message);

    class_calloc(pnl->z,
		 pnl->z_size,
		 sizeof(double),
		 pnl->error_message);

    class_calloc(pnl->p_density,
		 pnl->k_size*pnl->z_size,
		 sizeof(double),
		 pnl->error_message);
    class_calloc(pnl->p_cross,
		 pnl->k_size*pnl->z_size,
		 sizeof(double),
		 pnl->error_message);
    class_calloc(pnl->p_velocity,
		 pnl->k_size*pnl->z_size,
		 sizeof(double),
		 pnl->error_message);

    class_calloc(pnl->ddp_density,
		 pnl->k_size*pnl->z_size,
		 sizeof(double),
		 pnl->error_message);
    class_calloc(pnl->ddp_cross,
		 pnl->k_size*pnl->z_size,
		 sizeof(double),
		 pnl->error_message);
    class_calloc(pnl->ddp_velocity,
		 pnl->k_size*pnl->z_size,
		 sizeof(double),
		 pnl->error_message);

    for (index_k=0; index_k<pnl->k_size; index_k++) {

      pnl->k[index_k] = trg.k[index_k];

    }

    for (index_z=0; index_z<pnl->z_size; index_z++) {
      
      pnl->z[index_z] = trg.z[index_z];

      for (index_k=0; index_k<pnl->k_size; index_k++) {

	pnl->p_density[index_z*pnl->k_size+index_k]=trg.p_11_nl[index_z*pnl->k_size+index_k];
	pnl->p_cross[index_z*pnl->k_size+index_k]=trg.p_12_nl[index_z*pnl->k_size+index_k];
	pnl->p_velocity[index_z*pnl->k_size+index_k]=trg.p_22_nl[index_z*pnl->k_size+index_k];

      }
    }

    class_call(trg_free(&trg),
	       trg.error_message,
	       pnl->error_message);

    class_call(array_spline_table_lines(pnl->z,
					pnl->z_size,
					pnl->p_density,
					pnl->k_size,
					pnl->ddp_density,
					_SPLINE_EST_DERIV_,
					pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);

    class_call(array_spline_table_lines(pnl->z,
					pnl->z_size,
					pnl->p_cross,
					pnl->k_size,
					pnl->ddp_cross,
					_SPLINE_EST_DERIV_,
					pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);

    class_call(array_spline_table_lines(pnl->z,
					pnl->z_size,
					pnl->p_velocity,
					pnl->k_size,
					pnl->ddp_velocity,
					_SPLINE_EST_DERIV_,
					pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);

  }
  
  return _SUCCESS_;
}

int nonlinear_free(
		   struct nonlinear *pnl
		   ) {

  if (pnl->method > nl_none) {
    free(pnl->k);
    free(pnl->z);
    free(pnl->p_density);
    free(pnl->p_cross);
    free(pnl->p_velocity);
    free(pnl->ddp_density);
    free(pnl->ddp_cross);
    free(pnl->ddp_velocity);
  }

  return _SUCCESS_;

}
