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
		      double * pz_cross,
		      int * k_size_at_z
		      ) {

  int last_index;
  int index_z;
  int index_k;

  class_test(pnl->method == nl_none,
	     pnl->error_message,
	     "No non-linear spectra requested. You cannot call the function non_linear_pk_at_z()");

  if (pnl->z_size == 1) {

    class_test(z !=0,
	       pnl->error_message,
	       "if you want the non-linear pk at z=%g, you must set the input parameter z_max_pk to at least this value",z);
   
    for (index_k=0; index_k<pnl->k_size[0]; index_k++) {

      pz_density[index_k] = pnl->p_density[index_k];
      if ((pnl->method >= nl_trg_linear) && (pnl->method <= nl_trg)) {
	pz_velocity[index_k] = pnl->p_velocity[index_k];
	pz_cross[index_k] = pnl->p_cross[index_k];
      }
    }
  }

  else {

    class_call(array_interpolate_spline(pnl->z,
					pnl->z_size,
					pnl->p_density,
					pnl->ddp_density,
					pnl->k_size[0],
					z,
					&last_index,
					pz_density,
					pnl->k_size[0],
					pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    
    if ((pnl->method >= nl_trg_linear) && (pnl->method <= nl_trg)) {
      
      class_call(array_interpolate_spline(pnl->z,
					  pnl->z_size,
					  pnl->p_velocity,
					  pnl->ddp_velocity,
					  pnl->k_size[0],
					  z,
					  &last_index,
					  pz_velocity,
					  pnl->k_size[0],
					  pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);
      
      class_call(array_interpolate_spline(pnl->z,
					  pnl->z_size,
					  pnl->p_cross,
					  pnl->ddp_cross,
					  pnl->k_size[0],
					  z,
					  &last_index,
					  pz_cross,
					  pnl->k_size[0],
					  pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);
      
    }
  }

  for (index_z=0; pnl->z[index_z] > z; index_z++);
  * k_size_at_z = pnl->k_size[index_z];

  return _SUCCESS_;
}

int nonlinear_pk_at_k_and_z(
			    struct nonlinear * pnl,
			    double k,
			    double z,
			    double * pk_density,
			    double * pk_velocity,
			    double * pk_cross,
			    int * k_size_at_z
			    ) {
  
  double * pz_density;
  double * pz_velocity=NULL;
  double * pz_cross=NULL;
  double * ddpz_density;
  double * ddpz_velocity=NULL;
  double * ddpz_cross=NULL;
  int last_index;

  class_test(pnl->method == nl_none,
	     pnl->error_message,
	     "No non-linear spectra requested. You cannot call the function non_linear_pk_at_z()");

  class_alloc(pz_density,pnl->k_size[0]*sizeof(double),pnl->error_message);
  class_alloc(ddpz_density,pnl->k_size[0]*sizeof(double),pnl->error_message);

  if ((pnl->method >= nl_trg_linear) && (pnl->method <= nl_trg)) {

    class_alloc(pz_velocity,pnl->k_size[0]*sizeof(double),pnl->error_message);
    class_alloc(pz_cross,pnl->k_size[0]*sizeof(double),pnl->error_message);
    class_alloc(ddpz_velocity,pnl->k_size[0]*sizeof(double),pnl->error_message);
    class_alloc(ddpz_cross,pnl->k_size[0]*sizeof(double),pnl->error_message);

  }

  class_call(nonlinear_pk_at_z(pnl,z,pz_density,pz_velocity,pz_cross,k_size_at_z),
	     pnl->error_message,
	     pnl->error_message);

  class_call(array_spline_table_lines(pnl->k,
				      *k_size_at_z,
				      pz_density,
				      1,
				      ddpz_density,
				      _SPLINE_NATURAL_,
				      pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);
      
  class_call(array_interpolate_spline(pnl->k,
				      *k_size_at_z,
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

  if ((pnl->method >= nl_trg_linear) && (pnl->method <= nl_trg)) {

    class_call(array_spline_table_lines(pnl->k,
					*k_size_at_z,
					pz_velocity,
					1,
					ddpz_velocity,
					_SPLINE_NATURAL_,
					pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    
    class_call(array_interpolate_spline(pnl->k,
					*k_size_at_z,
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
					*k_size_at_z,
					pz_cross,
					1,
					ddpz_cross,
					_SPLINE_NATURAL_,
					pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    
    class_call(array_interpolate_spline(pnl->k,
					*k_size_at_z,
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

    free(pz_velocity);
    free(pz_cross);
    free(ddpz_velocity);
    free(ddpz_cross);
  }

  free(pz_density);
  free(ddpz_density);

  return _SUCCESS_;
}

int nonlinear_k_nl_at_z(
			struct nonlinear * pnl,
			double z,
			double * k_nl
			) {

  class_test(z>pnl->z[0],
	     pnl->error_message,
	     "you want to interpolate at z=%g greater than zmax=%g\n",z,pnl->z[0]);
	     
  if (pnl->z_size == 1) {
    *k_nl = pnl->k_nl[0];
  }
  else {
    class_call(array_interpolate_two(pnl->z,
				     1,
				     0,
				     pnl->k_nl,
				     1,
				     pnl->z_size,
				     z,
				     k_nl,
				     1,
				     pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
  }
  
  return _SUCCESS_;
}

int nonlinear_init(
		   struct precision *ppr,
		   struct background *pba,
		   struct thermo *pth,
		   struct perturbs *ppt,
		   struct bessels * pbs,
		   struct transfers * ptr,
		   struct primordial *ppm,
		   struct spectra *psp,
		   struct nonlinear *pnl
		   ) {

  int index_z,index_k;
  int last_density;
  int last_cross;
  int last_velocity;
  double * pk_ic=NULL;  /* array with argument 
		      pk_ic[index_k * psp->ic_ic_size[index_mode] + index_ic1_ic2] */

  double * pk_tot; /* array with argument 
		      pk_tot[index_k] */
  int index_ncdm;

  /** (a) if non non-linear corrections requested */

  if (pnl->method == nl_none) {
    if (pnl->nonlinear_verbose > 0)
      printf("No non-linear spectra requested. Nonlinear module skipped.\n");
  }

  /** (b) for HALOFIT non-linear spectrum */

  else if (pnl->method == nl_halofit) {
    if (pnl->nonlinear_verbose > 0)
      printf("Computing non-linear matter power spectrum with Halofit (including update by Bird et al 2011)\n");

    if (pba->has_ncdm) {
      for (index_ncdm=0;index_ncdm < pba->N_ncdm; index_ncdm++){
	if (pba->m_ncdm_in_eV[index_ncdm] >  _M_EV_TOO_BIG_FOR_HALOFIT_)
	  fprintf(stdout,"Warning: Halofit is proved to work for CDM, and also with a small HDM component thanks to Bird et al.'s update. But it sounds like you are running with a WDM component of mass %f eV, which makes the use of Halofit suspicious.\n",pba->m_ncdm_in_eV[index_ncdm]);
      }
    }

    /** - define values of z */

    pnl->z_size = (int)(psp->z_max_pk/ppr->halofit_dz)+1;

    class_alloc(pnl->z,
		pnl->z_size*sizeof(double),
		pnl->error_message);

    class_alloc(pnl->k_nl,
		pnl->z_size*sizeof(double),
		pnl->error_message);

    if (pnl->z_size == 1) {
      pnl->z[0]=0.;
    }
    else {
      for (index_z=0; index_z < pnl->z_size; index_z++) {
	pnl->z[index_z]=(double)(pnl->z_size-1-index_z)/(double)(pnl->z_size-1)*psp->z_max_pk;  /* z stored in decreasing order */
      }
    }

    /** - define values of k */

    class_alloc(pnl->k_size,
		pnl->z_size*sizeof(int),
		pnl->error_message);
    
    for (index_z=0; index_z < pnl->z_size; index_z++) {
      pnl->k_size[index_z] = psp->ln_k_size;
    }

    class_alloc(pnl->k,
		pnl->k_size[0]*sizeof(double),
		pnl->error_message);

    for (index_k=0; index_k < pnl->k_size[0]; index_k++) {
      pnl->k[index_k] = exp(psp->ln_k[index_k]);
    }

    /** - allocate p_density (= pk_nonlinear) and fill it with linear power spectrum */

    class_alloc(pnl->p_density,
		pnl->k_size[0]*pnl->z_size*sizeof(double),
		pnl->error_message);

    class_alloc(pk_tot,
		psp->ln_k_size*sizeof(double),
		pnl->error_message);

    if (psp->ic_size[psp->index_md_scalars] > 1) {

      class_alloc(pk_ic,
		  psp->ln_k_size*psp->ic_ic_size[psp->index_md_scalars]*sizeof(double),
		  pnl->error_message);

    }

    for (index_z=0; index_z < pnl->z_size; index_z++) {

      class_call(spectra_pk_at_z(pba,
				 psp,
				 linear,
				 pnl->z[index_z],
				 pk_tot,
				 pk_ic),
		 psp->error_message,
		 pnl->error_message);

      for (index_k=0; index_k < pnl->k_size[index_z]; index_k++) {
	
	pnl->p_density[index_z*pnl->k_size[index_z]+index_k] = pk_tot[index_k];

      }      
    }

    free(pk_tot);
    if (psp->ic_size[psp->index_md_scalars] > 1) {
      free(pk_ic);
    }

    /** - apply non-linear corrections */

    class_call(nonlinear_halofit(ppr,pba,ppm,psp,pnl),
	       pnl->error_message,
	       pnl->error_message);

    /** - take second derivative w.r.t z in view of spline inteprolation */

    if (pnl->z_size > 1) {

      class_alloc(pnl->ddp_density,
		  pnl->k_size[0]*pnl->z_size*sizeof(double),
		  pnl->error_message);

      class_call(array_spline_table_lines(pnl->z,
					  pnl->z_size,
					  pnl->p_density,
					  pnl->k_size[0],
					  pnl->ddp_density,
					  _SPLINE_EST_DERIV_,
					  pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);
      
    }
  }

  /** (c) for TRG non-linear spectrum */

  else if ((pnl->method >= nl_trg_linear) && (pnl->method <= nl_trg)) {
    if (pnl->nonlinear_verbose > 0)
      printf("Computing non-linear matter power spectrum with trg module\n");

    struct spectra_nl trg;

    if (pnl->method == nl_trg_linear)
      trg.mode = 0;
    if (pnl->method == nl_trg_one_loop)
      trg.mode = 1;
    if (pnl->method == nl_trg)
      trg.mode = 2;

    trg.k_max = exp(psp->ln_k[psp->ln_k_size-1]) * pba->h - 1.;

    trg.double_escape = ppr->double_escape;
    trg.z_ini = ppr->z_ini;
    trg.eta_size = ppr->eta_size;
    trg.k_L = ppr->k_L;
    trg.k_min = ppr->k_min;
    trg.logstepx_min = ppr->logstepx_min;
    trg.k_growth_factor = ppr->k_growth_factor;
    trg.ic = pnl->ic;

    trg.spectra_nl_verbose = pnl->nonlinear_verbose;

    class_call(trg_init(ppr,pba,pth,ppm,psp,&trg),
	       trg.error_message,
	       pnl->error_message);

    /* copy non-linear spectrum in pnl */
    
    pnl->z_size = (trg.eta_size-1)/2+1;

    class_calloc(pnl->k_size,
		 pnl->z_size,
		 sizeof(int),
		 pnl->error_message);

    for (index_z=0; index_z < pnl->z_size; index_z++) {
      pnl->k_size[index_z] = trg.k_size-4*trg.double_escape*index_z;
    }

    class_calloc(pnl->k,
		 pnl->k_size[0],
		 sizeof(double),
		 pnl->error_message);

    class_calloc(pnl->z,
		 pnl->z_size,
		 sizeof(double),
		 pnl->error_message);

    class_calloc(pnl->p_density,
		 pnl->k_size[0]*pnl->z_size,
		 sizeof(double),
		 pnl->error_message);
    class_calloc(pnl->p_cross,
		 pnl->k_size[0]*pnl->z_size,
		 sizeof(double),
		 pnl->error_message);
    class_calloc(pnl->p_velocity,
		 pnl->k_size[0]*pnl->z_size,
		 sizeof(double),
		 pnl->error_message);

    class_calloc(pnl->ddp_density,
		 pnl->k_size[0]*pnl->z_size,
		 sizeof(double),
		 pnl->error_message);
    class_calloc(pnl->ddp_cross,
		 pnl->k_size[0]*pnl->z_size,
		 sizeof(double),
		 pnl->error_message);
    class_calloc(pnl->ddp_velocity,
		 pnl->k_size[0]*pnl->z_size,
		 sizeof(double),
		 pnl->error_message);

    for (index_k=0; index_k<pnl->k_size[0]; index_k++) {

      pnl->k[index_k] = trg.k[index_k];

    }

    for (index_z=0; index_z<pnl->z_size; index_z++) {
      
      pnl->z[index_z] = trg.z[2*index_z];

      for (index_k=0; index_k<pnl->k_size[0]; index_k++) {

	pnl->p_density[index_z*pnl->k_size[0]+index_k]=trg.p_11_nl[2*index_z*pnl->k_size[0]+index_k];
	pnl->p_cross[index_z*pnl->k_size[0]+index_k]=trg.p_12_nl[2*index_z*pnl->k_size[0]+index_k];
	pnl->p_velocity[index_z*pnl->k_size[0]+index_k]=trg.p_22_nl[2*index_z*pnl->k_size[0]+index_k];

      }
    }

    /* for non-computed values: instead oif leaving zeros, leave last
       computed value: the result is more smooth and will not fool the
       interpolation routine; but these values are not outputed at the
       end. In order to have an even better intrpolation, would be
       better to extrapolate with a parabola rather than a
       constant. */
 
    for (index_k=0; index_k<pnl->k_size[0]; index_k++) {

      last_density = pnl->p_density[index_k];
      last_cross = pnl->p_cross[index_k];
      last_velocity = pnl->p_velocity[index_k];

      for (index_z=0; index_z<pnl->z_size; index_z++) {

	if (pnl->p_density[index_z*pnl->k_size[0]+index_k] == 0.)
	  pnl->p_density[index_z*pnl->k_size[0]+index_k] = last_density;
	else
	  last_density = pnl->p_density[index_z*pnl->k_size[0]+index_k];

	if (pnl->p_cross[index_z*pnl->k_size[0]+index_k] == 0.)
	  pnl->p_cross[index_z*pnl->k_size[0]+index_k] = last_cross;
	else
	  last_cross = pnl->p_cross[index_z*pnl->k_size[0]+index_k];
	
	if (pnl->p_velocity[index_z*pnl->k_size[0]+index_k] == 0.)
	  pnl->p_velocity[index_z*pnl->k_size[0]+index_k] = last_velocity;
	else
	  last_velocity = pnl->p_velocity[index_z*pnl->k_size[0]+index_k];
      }
    }


    class_call(trg_free(&trg),
	       trg.error_message,
	       pnl->error_message);

    class_call(array_spline_table_lines(pnl->z,
					pnl->z_size,
					pnl->p_density,
					pnl->k_size[0],
					pnl->ddp_density,
					_SPLINE_EST_DERIV_,
					pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);

    class_call(array_spline_table_lines(pnl->z,
					pnl->z_size,
					pnl->p_cross,
					pnl->k_size[0],
					pnl->ddp_cross,
					_SPLINE_EST_DERIV_,
					pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);

    class_call(array_spline_table_lines(pnl->z,
					pnl->z_size,
					pnl->p_velocity,
					pnl->k_size[0],
					pnl->ddp_velocity,
					_SPLINE_EST_DERIV_,
					pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);

  }

  else {
    class_stop(pnl->error_message,
		"Your non-linear method variable is set to %d, out of the range defined in nonlinear.h",pnl->method);
  }
  
  return _SUCCESS_;
}

int nonlinear_free(
		   struct nonlinear *pnl
		   ) {

  if (pnl->method > nl_none) {

    free(pnl->k);
    free(pnl->k_size);
    free(pnl->z);
    free(pnl->k_nl);
    free(pnl->p_density);
    if (pnl->z_size > 1) {
      free(pnl->ddp_density);
    }

    if ((pnl->method >= nl_trg_linear) && (pnl->method <= nl_trg)) {
      free(pnl->p_cross);
      free(pnl->p_velocity);
      free(pnl->ddp_cross);
      free(pnl->ddp_velocity);
    }
  }

  return _SUCCESS_;

}

int nonlinear_halofit(
		      struct precision *ppr,
		      struct background *pba,
		      struct primordial *ppm,
		      struct spectra *psp,
		      struct nonlinear *pnl
		      ) {

  double Omega_m,Omega_v,fnu,Omega0_m;

  /** determine non linear ratios (from pk) **/

  int index_z,index_k;
  double pk_lin,pk_quasi,pk_halo,rk;
  double sigma,rknl,rneff,rncur,d1,d2;
  double diff,xlogr1,xlogr2,rmid;

  double extragam,gam,a,b,c,xmu,xnu,alpha,beta,f1,f2,f3;
  double pk_linaa;
  double y;
  double f1a,f2a,f3a,f1b,f2b,f3b,frac;

  double * pvecback;

  int last_index;
  int counter;
  double sum1,sum2,sum3;
  double anorm;
  double tau;

  double *integrand_array;

  double x2,R;

  class_alloc(pvecback,pba->bg_size*sizeof(double),pnl->error_message);

  Omega0_m = (pba->Omega0_cdm + pba->Omega0_b + pba->Omega0_ncdm_tot);
  fnu      = pba->Omega0_ncdm_tot/(pba->Omega0_b + pba->Omega0_cdm+pba->Omega0_ncdm_tot);
  anorm    = 1./(2*pow(_PI_,2));

  class_alloc(integrand_array,pnl->k_size[0]*7*sizeof(double),pnl->error_message);

  for (index_k=0; index_k < pnl->k_size[0]; index_k++) {
    integrand_array[index_k*7 + 0] = pnl->k[index_k];
  }

  for (index_z=0; index_z < pnl->z_size; index_z++) {

    class_call(background_tau_of_z(pba,pnl->z[index_z],&tau),pba->error_message,pnl->error_message);
    class_call(background_at_tau(pba,tau,pba->long_info,pba->inter_normal,&last_index,pvecback),
	pba->error_message,
	pnl->error_message);

    Omega_m = pvecback[pba->index_bg_Omega_m];
    Omega_v = 1.-pvecback[pba->index_bg_Omega_m]-pvecback[pba->index_bg_Omega_r];

    /* minimum value of R such that the integral giving sigma_R is converged */
    R=sqrt(-log(ppr->halofit_sigma_precision))/pnl->k[pnl->k_size[0]-1];

    /* corresponding value of sigma_R */
    sum1=0.;
    for (index_k=0; index_k < pnl->k_size[index_z]; index_k++) {
      x2 = pow(pnl->k[index_k]*R,2);
      integrand_array[index_k*7 + 1] = pnl->p_density[index_z*pnl->k_size[index_z]+index_k]*pow(pnl->k[index_k],2)*anorm*exp(-x2);
    }
    /* fill in second derivatives */
    class_call(array_spline(integrand_array,7,pnl->k_size[index_z],0,1,2,_SPLINE_EST_DERIV_,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    /* integrate */
    class_call(array_integrate_all_spline(integrand_array,7,pnl->k_size[index_z],0,1,2,&sum1,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    sigma  = sqrt(sum1);

    class_test(sigma < 1.,
	       pnl->error_message,
	       "Your k_max=%g 1/Mpc is too small for Halofit to find the non-linearity scale z_nl at z=%g. Increase input parameter P_k_max_h/Mpc or P_k_max_1/Mpc",
	       pnl->k[pnl->k_size[0]-1],
	       pnl->z[index_z]);  

    xlogr1 = log(R)/log(10.);

    /* maximum value of R in the bisection algorithm leading to the determination of R_nl */
    R=1./ppr->halofit_min_k_nonlinear;
    
    /* corresponding value of sigma_R */
    sum1=0.;
    for (index_k=0; index_k < pnl->k_size[index_z]; index_k++) {
      x2 = pow(pnl->k[index_k]*R,2);
      integrand_array[index_k*7 + 1] = pnl->p_density[index_z*pnl->k_size[index_z]+index_k]*pow(pnl->k[index_k],2)*anorm*exp(-x2);
    }
    /* fill in second derivatives */
    class_call(array_spline(integrand_array,7,pnl->k_size[index_z],0,1,2,_SPLINE_EST_DERIV_,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    /* integrate */
    class_call(array_integrate_all_spline(integrand_array,7,pnl->k_size[index_z],0,1,2,&sum1,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    sigma  = sqrt(sum1);

    class_test(sigma > 1.,
	       pnl->error_message,
	       "Your input value for the precision parameter halofit_min_k_nonlinear=%e is too large, the non-linear wavenumber k_nl must be smaller than that",
	       ppr->halofit_min_k_nonlinear);

    xlogr2 = log(R)/log(10.);

    counter = 0;
    do {
      rmid = pow(10,(xlogr2+xlogr1)/2.0);
      counter ++;

      /* in original halofit, this is the function wint() */
      sum1=0.;
      sum2=0.;
      sum3=0.;

      for (index_k=0; index_k < pnl->k_size[index_z]; index_k++) {
	x2 = pnl->k[index_k]*pnl->k[index_k]*rmid*rmid;
	integrand_array[index_k*7 + 1] = pnl->p_density[index_z*pnl->k_size[index_z]+index_k]*pow(pnl->k[index_k],2)*anorm*exp(-x2);
	integrand_array[index_k*7 + 3] = pnl->p_density[index_z*pnl->k_size[index_z]+index_k]*pow(pnl->k[index_k],2)*anorm*2.*x2*exp(-x2);
	integrand_array[index_k*7 + 5] = pnl->p_density[index_z*pnl->k_size[index_z]+index_k]*pow(pnl->k[index_k],2)*anorm*4.*x2*(1.-x2)*exp(-x2);
      }
      /* fill in second derivatives */
      class_call(array_spline(integrand_array,7,pnl->k_size[index_z],0,1,2,_SPLINE_EST_DERIV_,pnl->error_message),
	  pnl->error_message,
	  pnl->error_message);
      class_call(array_spline(integrand_array,7,pnl->k_size[index_z],0,3,4,_SPLINE_EST_DERIV_,pnl->error_message),
	  pnl->error_message,
	  pnl->error_message);
      class_call(array_spline(integrand_array,7,pnl->k_size[index_z],0,5,6,_SPLINE_EST_DERIV_,pnl->error_message),
	  pnl->error_message,
	  pnl->error_message);

      /* integrate */
      class_call(array_integrate_all_spline(integrand_array,7,pnl->k_size[index_z],0,1,2,&sum1,pnl->error_message),
	  pnl->error_message,
	  pnl->error_message);
      class_call(array_integrate_all_spline(integrand_array,7,pnl->k_size[index_z],0,3,4,&sum2,pnl->error_message),
	  pnl->error_message,
	  pnl->error_message);
      class_call(array_integrate_all_spline(integrand_array,7,pnl->k_size[index_z],0,5,6,&sum3,pnl->error_message),
	  pnl->error_message,
	  pnl->error_message);

      sigma  = sqrt(sum1);
      d1 = -sum2/sum1;
      d2 = -sum2*sum2/sum1/sum1 - sum3/sum1;
      /* in original halofit, this is the end of the function wint() */
 
      diff = sigma - 1.0;
      /*fprintf(stderr,"xlogr1 = %g, xlogr2 = %g, rmid = %g, diff: =%g, abs(diff) = %g\n",xlogr1,xlogr2,log10(rmid),diff,fabs(diff));*/
      if (diff>0.001){
	xlogr1=log10(rmid);
	/*fprintf(stderr,"going up  , new xlogr1=%g\n",xlogr1);*/
      }
      else if (diff < -0.001) {
	xlogr2 = log10(rmid);
	/*fprintf(stderr,"going down, new xlogr2=%g\n",xlogr2);*/
      }

      class_test(counter > _MAX_IT_,
		 pnl->error_message,
		 "could not converge within maximum allowed number of iterations");

    } while (fabs(diff) > 0.001);

    rknl  = 1./rmid;
    rneff = -3.-d1;
    rncur = -d2;

    pnl->k_nl[index_z] = rknl;

    for (index_k = 0; index_k < pnl->k_size[index_z]; index_k++){
      
      rk = pnl->k[index_k];
      
      if (rk > ppr->halofit_min_k_nonlinear) {
	
	pk_lin = pnl->p_density[index_z*pnl->k_size[index_z]+index_k]*pow(pnl->k[index_k],3)*anorm;
	
	/* in original halofit, this is the beginning of the function halofit() */
	
	/*SPB11: Standard halofit underestimates the power on the smallest
	   * scales by a factor of two. Add an extra correction from the
	   * simulations in Bird, Viel,Haehnelt 2011 which partially accounts for
	   * this.*/
	
	extragam = 0.3159 -0.0765*rneff -0.8350*rncur;
	gam=extragam+0.86485+0.2989*rneff+0.1631*rncur;
	
	a=1.4861+1.83693*rneff+1.67618*rneff*rneff+0.7940*rneff*rneff*rneff+0.1670756*rneff*rneff*rneff*rneff-0.620695*rncur;
	a=pow(10,a);
	b=pow(10,(0.9463+0.9466*rneff+0.3084*rneff*rneff-0.940*rncur));
	c=pow(10,(-0.2807+0.6669*rneff+0.3214*rneff*rneff-0.0793*rncur));
	xmu   = pow(10,(-3.54419+0.19086*rneff));
	xnu   = pow(10,(0.95897+1.2857*rneff));
	alpha = 1.38848+0.3701*rneff-0.1452*rneff*rneff;
	beta  = 0.8291+0.9854*rneff+0.3400*pow(rneff,2)+fnu*(-6.4868+1.4373*pow(rneff,2));
	if(fabs(1-Omega_m)>0.01) { /*then omega evolution */
	  f1a=pow(Omega_m,(-0.0732));
	  f2a=pow(Omega_m,(-0.1423));
	  f3a=pow(Omega_m,(0.0725));
	  f1b=pow(Omega_m,(-0.0307));
	  f2b=pow(Omega_m,(-0.0585));
	  f3b=pow(Omega_m,(0.0743));     
	  frac=Omega_v/(1.-Omega_m); 
	  f1=frac*f1b + (1-frac)*f1a;
	  f2=frac*f2b + (1-frac)*f2a;
	  f3=frac*f3b + (1-frac)*f3a;
	}
	else {
	  f1=1.;
	  f2=1.;
	  f3=1.;
	}
	
	y=(rk/rknl);
	pk_halo = a*pow(y,f1*3.)/(1.+b*pow(y,f2)+pow(f3*c*y,3.-gam));
	pk_halo=pk_halo/(1+xmu*pow(y,-1)+xnu*pow(y,-2))*(1+fnu*(2.080-12.39*(Omega0_m-0.3))/(1+1.201e-03*pow(y,3)));
	pk_linaa=pk_lin*(1+fnu*26.29*pow(rk,2)/(1+1.5*pow(rk,2)));
	pk_quasi=pk_lin*pow((1+pk_linaa),beta)/(1+pk_linaa*alpha)*exp(-y/4.0-pow(y,2)/8.0);
	
	pnl->p_density[index_z*pnl->k_size[index_z]+index_k] = (pk_halo+pk_quasi)/pow(pnl->k[index_k],3)/anorm;
	
	/* in original halofit, this is the end of the function halofit() */
      }
    }
  }						
  free(pvecback);
  free(integrand_array);
  return _SUCCESS_;
}
