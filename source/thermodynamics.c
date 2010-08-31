/** @file thermodynamics.c Documented thermodynamics module
 * Julien Lesgourgues, 18.04.2010    
 *
 * Deals with the thermodynamical evolution.
 * This module has two purposes: 
 *
 * - at the beginning, to initialize the thermodynamics, i.e. to
 integrate the thermodynamical equations, and store all
 thermodynamical quantities as a function of redshift inside an
 interpolation table. The current version is based on RECFAST.
 *
 * - at any time in the code, to evaluate any thermodynamical quantity
 for a given redshft value (by interpolating within the interpolation
 table).
 *
 * Hence the following functions can be called from other modules:
 *
 * -# thermodynamics_init() at the beginning (but after background_init()) 
 * -# thermodynamics_at_z() at any later time
 * -# thermodynamics_free() at the end, when no more calls to thermodynamics_at_z() are needed
 */

#include "thermodynamics.h"

/** 
 * Thermodynamics quantities at given redshift z. 
 *
 * Evaluates all thermodynamics quantities at a given value of
 * the redshift by reading the pre-computed table ant interpolating.
 * This function can be called from whatever module at whatever time,
 * provided that thermodynamics_init() has been called before, and 
 * thermodynamics_free() has not been called yet.
 *
 * @param z Input: redshift
 * @param intermode Input: interpolation mode (normal or growing_closeby)
 * @param last_index Input/Ouput: index of the previous/current point in the interpolation array (input only for closeby mode, output for both) 
 * @param pvecthermo_local Output: vector (assumed to be already allocated) of thermodynamics quantities
 * @return the error status
 */
int thermodynamics_at_z(
			struct background * pba,
			struct thermo * pth,
			double z,
			enum interpolation_mode intermode,
			int * last_index,
			double * pvecback,
			double * pvecthermo
			) {

  /** Summary: */

  /** - define local variables */
  int i;
  double x0;

  /** - chech that z is in the pre-computed range; if \f$ z>z_{max}
      \f$, extrapolate */

  /* deal with case z < z_min=0, unphysical */
  class_test(z < pth->z_table[0],
	     pth->error_message,"out of range: z < z_min=%e",pth->z_table[0]);

  /* deal with case z > z_max, not stored inside interpolation table,
     but doable with simple analytic approximations */
  if (z >= pth->z_table[pth->tt_size-1]) {

    x0= pth->thermodynamics_table[(pth->tt_size-1)*pth->th_size+pth->index_th_xe];
    
    /* ionization fraction */
    pvecthermo[pth->index_th_xe] = x0;

    /* Calculate Tb */
    pvecthermo[pth->index_th_Tb] = pth->Tcmb*(1.+z);

    /* Calculate cb2 (cb2 = (k_B/mu) Tb (1-1/3 dlnTb/dlna) = (k_B/mu) Tb (1+1/3 (1+z) dlnTb/dz)) */
    /* note that m_H / mu = 1 + (m_H/m_He-1) Y_p + x_e (1-Y_p) */
    pvecthermo[pth->index_th_cb2] = _k_B_ / ( _c_ * _c_ * _m_H_ ) * (1. + (1./_not4_ - 1.) * pth->YHe + x0 * (1.-pth->YHe)) * pth->Tcmb * (1.+z) * 4. / 3.;

    /* Calculate dkappa/deta (dkappa/deta = a n_e x_e sigma_T = a^{-2} n_e(today) x_e sigma_T in units of 1/Mpc) */
    pvecthermo[pth->index_th_dkappa] = (1.+z) * (1.+z) * pth->n_e * x0 * _sigma_ * _Mpc_over_m_;


    /* Calculate dz/deta = -H with background_functions() */
    class_call(background_functions(pba,1./(1.+z),short_info,pvecback),
	       pba->error_message,
	       pth->error_message);

    /* Calculate d2kappa/deta2 = dz/deta d/dz[dkappa/deta] */
    pvecthermo[pth->index_th_ddkappa] = -pvecback[pba->index_bg_H] * 2. / (1.+z) * pvecthermo[pth->index_th_dkappa];

    /* Calculate d3kappa/deta3 = dz/deta d/dz[d2kappa/deta2] */
    /*     pvecthermo[pth->index_th_ddkappa] = (pvecback[pba->index_bg_Hdot] / pvecback[pba->index_bg_H] - */
    /*        pvecback[pba->index_bg_H]	/ (1.+z)) * pvecthermo[pth->index_th_ddkappa]; */

    /* visibility = g = (dkappa/deta) * exp(- kappa) */
    pvecthermo[pth->index_th_g]=0.;

  }

  /** - otherwise, just interpolate in table with array_interpolate() */
  else {

    if (intermode == normal) {

      class_call(array_interpolate_spline(
					  pth->z_table,
					  pth->tt_size,
					  pth->thermodynamics_table,
					  pth->d2thermodynamics_dz2_table,
					  pth->th_size,
					  z,
					  last_index,
					  pvecthermo,
					  pth->th_size,
					  pth->error_message),
		 pth->error_message,
		 pth->error_message);

    }

    if (intermode == closeby) {

      class_call(array_interpolate_spline_growing_closeby(
							  pth->z_table,
							  pth->tt_size,
							  pth->thermodynamics_table,
							  pth->d2thermodynamics_dz2_table,
							  pth->th_size,
							  z,
							  last_index,
							  pvecthermo,
							  pth->th_size,
							  pth->error_message),
		 pth->error_message,
		 pth->error_message);
    }

  }

  return _SUCCESS_;
}

/** 
 * Initialize the thermo structure, including thermodynamics interpolation table.
 * 
 * Reads the cosmological parameters and initialize all fields in the
 * structure thermo, in particular:
 *
 * - initializes all indices in the thermodynamics vectors with thermodynamics_indices()
 *
 * - integrates the thermodynamics using thermodynamics_recombination_with_recfast(), and initializes the thermodynamics interpolation table using array_derive(), thermodynamics_cure_discontinuity() and array_integrate_ratio()().
 *
 * This function shall be called at the beginning of each run, but
 * only after background_init(). It allocates memory spaces which
 * should be freed later with thermodynamics_free().
 *
 * @param pba_input Input : Initialized background structure
 * @param ppr_input Input : Parameters describing how the computation is to be performed
 * @param pth_output Output : Initialized thermodynamics structure
 * @return the error status
 */
int thermodynamics_init(
			struct precision * ppr,
			struct background * pba,
			struct thermo * pth
			) {

  /** Summary: */

  /** - define local variables */

  /* running index over vector of thermodynamics variables */
  int i, index_th, index_re;
  /* visibility function */
  double g, g_previous, eta_visibility_max;
  /* for calling background_at_eta() */
  int last_index_back;
  double * eta_table; /**< list of eta values associated with z values in pth->z_table */

  struct recombination reco;
  struct reionization reio;
  struct recombination * preco;
  struct reionization * preio;

  double * pvecback;

  preco=&reco;
  preio=&reio;

  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);

  if (pth->thermodynamics_verbose > 0)
    printf("Computing thermodynamics\n");

  /* Tcmb in K */
  class_test((pth->Tcmb < _TCMB_SMALL_)||(pth->Tcmb > _TCMB_BIG_),
	     pth->error_message,
	     "Tcmb=%g out of bounds (%g<Tcmb<%g)",pth->Tcmb,_TCMB_SMALL_,_TCMB_BIG_);

  class_test(fabs(pba->Omega0_g/((4.*_sigma_B_/_c_*pow(pth->Tcmb,4.)) / (3.*_c_*_c_*1.e10*pba->h*pba->h/_Mpc_over_m_/_Mpc_over_m_/8./_PI_/_G_))-1.) > ppr->smallest_allowed_variation,
	     pth->error_message,
	     "inconsistency between photon temperature and density (fixed by stefan-Boltzmann law): you have Tcmb=%f K, but Omega0_g=%e and omega0_g=%e",
	     pth->Tcmb,
	     pba->Omega0_g,
	     pba->Omega0_g*pba->h*pba->h);

  /* Y_He */
  class_test((pth->YHe < _YHE_SMALL_)||(pth->YHe > _YHE_BIG_),
	     pth->error_message,
	     "Y_He=%g out of bounds (%g<Y_He<%g)",pth->YHe,_YHE_SMALL_,_YHE_BIG_);

  /* tests in order to prevent segmentation fault in the following */
  class_test(_not4_ == 0.,
	     pth->error_message,
	     "stop to avoid division by zero");
  class_test(pth->YHe == 1.,
	     pth->error_message,
	     "stop to avoid division by zero");

  /** - assign values to all indices in the structures with thermodynamics_indices()*/
  class_call(thermodynamics_indices(pth,preco,preio),
	     pth->error_message,
	     pth->error_message);

  /** - solve recombination and store values of \f$ z, x_e, d \kappa / d \eta, T_b, c_b^2 \f $ with thermodynamics_recombination() */
  class_call(thermodynamics_recombination(ppr,pba,pth,preco,pvecback),
	     pth->error_message,
	     pth->error_message);

  /** - solve reionization and store values of \f$ z, x_e, d \kappa / d z \f $ with thermodynamics_recombination()*/
  if (pth->reio_parametrization != reio_none) {
    class_call(thermodynamics_reionization(ppr,pba,pth,preco,preio,pvecback),
	       pth->error_message,
	       pth->error_message);
  }
  else {
    preio->rt_size=0;
    preio->index_reco_when_reio_start=-1;
  }

  /*   FILE * output; */

  /*   output=fopen("test_output/reco","w"); */
  /*   for (i=0;i < preco->rt_size;i++) { */
  /*     fprintf(output,"%e %e\n", */
  /* 	    preco->recombination_table[i*preco->re_size+preco->index_re_z], */
  /* 	    preco->recombination_table[i*preco->re_size+preco->index_re_xe]); */
  /*   } */
  /*   fclose(output); */
  /*   output=fopen("test_output/reio","w"); */
  /*   for (i=0;i < preio->rt_size;i++) { */
  /*     fprintf(output,"%e %e\n", */
  /* 	    preio->reionization_table[i*preio->re_size+preio->index_re_z], */
  /* 	    preio->reionization_table[i*preio->re_size+preio->index_re_xe]); */
  /*   } */
  /*   fclose(output); */


  /** - allocate memory for thermodynamics interpolation tables */

  /* first, a little check that the two tables match each other and can be merged */
  if (pth->reio_parametrization != reio_none) {
    class_test(preco->recombination_table[preio->index_reco_when_reio_start*preco->re_size+preco->index_re_z] !=
	       preio->reionization_table[(preio->rt_size -1)*preio->re_size+preio->index_re_z],
	       pth->error_message,
	       "mismatch which should never happen");
  }

  /* number of redshift in full table = number in reco + number in reio - overlap */
  pth->tt_size = ppr->recfast_Nz0 + preio->rt_size - preio->index_reco_when_reio_start - 1;
  class_alloc(pth->z_table,pth->tt_size*sizeof(double),pth->error_message);
  class_alloc(eta_table,pth->tt_size*sizeof(double),pth->error_message);
  class_alloc(pth->thermodynamics_table,pth->th_size*pth->tt_size*sizeof(double),pth->error_message);
  class_alloc(pth->d2thermodynamics_dz2_table,pth->th_size*pth->tt_size*sizeof(double),pth->error_message);  

  for (i=0; i < preio->rt_size; i++) {
    pth->z_table[i]=
      preio->reionization_table[i*preio->re_size+preio->index_re_z];
    pth->thermodynamics_table[i*pth->th_size+pth->index_th_xe]=
      preio->reionization_table[i*preio->re_size+preio->index_re_xe];
    pth->thermodynamics_table[i*pth->th_size+pth->index_th_dkappa]=
      preio->reionization_table[i*preio->re_size+preio->index_re_dkappadeta];
    pth->thermodynamics_table[i*pth->th_size+pth->index_th_Tb]=
      preio->reionization_table[i*preio->re_size+preio->index_re_Tb];
    pth->thermodynamics_table[i*pth->th_size+pth->index_th_cb2]=
      preio->reionization_table[i*preio->re_size+preio->index_re_cb2];
  }
  for (i=0; i < ppr->recfast_Nz0 - preio->index_reco_when_reio_start - 1; i++) {
    index_th=i+preio->rt_size;
    index_re=i+preio->index_reco_when_reio_start+1;
    pth->z_table[index_th]=
      preco->recombination_table[index_re*preco->re_size+preco->index_re_z];
    pth->thermodynamics_table[index_th*pth->th_size+pth->index_th_xe]=
      preco->recombination_table[index_re*preco->re_size+preco->index_re_xe];
    pth->thermodynamics_table[index_th*pth->th_size+pth->index_th_dkappa]=
      preco->recombination_table[index_re*preco->re_size+preco->index_re_dkappadeta];
    pth->thermodynamics_table[index_th*pth->th_size+pth->index_th_Tb]=
      preco->recombination_table[index_re*preco->re_size+preco->index_re_Tb];
    pth->thermodynamics_table[index_th*pth->th_size+pth->index_th_cb2]=
      preco->recombination_table[index_re*preco->re_size+preco->index_re_cb2];
  }

  free(preco->recombination_table);

  if (pth->reio_parametrization != reio_none)
    free(preio->reionization_table);

  /** -compute table of corresponding conformal times */

  for (i=0; i < pth->tt_size; i++) {
    class_call(background_eta_of_z(pba,pth->z_table[i],eta_table+i),
	       pba->error_message,
	       pth->error_message);
  }

  /** - fill missing columns (quantities not computed previously but
      related); eventually, reduce numerical errors with
      thermodynamics_cure_discontinuity()
  */

  /* -> second derivative with respect to eta of dkappa */
  /** - fill tables of second derivatives (in view of spline interpolation) */
  class_call(array_spline_table_line_to_line(eta_table,
					     pth->tt_size,
					     pth->thermodynamics_table,
					     pth->th_size,
					     pth->index_th_dkappa,
					     pth->index_th_dddkappa,
					     _SPLINE_EST_DERIV_,
					     pth->error_message),
	     pth->error_message,
	     pth->error_message);

  /* -> first derivative with respect to eta of dkappa */
  class_call(array_derive_spline_table_line_to_line(eta_table,
						    pth->tt_size,
						    pth->thermodynamics_table,
						    pth->th_size,
						    pth->index_th_dkappa,
						    pth->index_th_dddkappa,
						    pth->index_th_ddkappa,
						    pth->error_message),
	     pth->error_message,
	     pth->error_message);

  /* -> compute -kappa = [int_{eta_today}^{eta} deta dkappa/deta], store temporarily in column "g" */ 
  class_call(array_integrate_spline_table_line_to_line(eta_table,
						       pth->tt_size,
						       pth->thermodynamics_table,
						       pth->th_size,
						       pth->index_th_dkappa,
						       pth->index_th_dddkappa,
						       pth->index_th_g,
						       pth->error_message),
	     pth->error_message,
	     pth->error_message);

  free(eta_table);

  /** 
   * - compute visibility : \f$ g= (d \kappa/d \eta) e^{- \kappa} \f$; compute also redshift when (in order of growing \f$ z \f$): 
   *      -# \f$ g \f$ cannot be approximated by zero anymore (time for switching on source function)
   *      -# \f$ g \f$ is maximum (recombination time)
   *      -# \f$ g \f$ falls again below some threshold (can be used for switching on free-streaming approximation, but not in this version)
   * 
   */

  pth->z_visibility_start_sources=0.;
  pth->z_visibility_free_streaming=0.;
  pth->z_visibility_max=0.;
  g = 0.;

  /* loop on z (decreasing z, increasing time) */
  for (i=pth->tt_size-1; i>=0; i--) {

    g_previous = g;
    g = pth->thermodynamics_table[i*pth->th_size+pth->index_th_dkappa] *
      exp(pth->thermodynamics_table[i*pth->th_size+pth->index_th_g]);

    if (g > 0.) {
      if (pth->z_visibility_start_sources == 0.) {
	if (g >= ppr->visibility_threshold_start_sources)
	  pth->z_visibility_start_sources=pth->z_table[i];
      }
      else {
	if (pth->z_visibility_max ==0.) 
	  if (g < g_previous)
	    pth->z_visibility_max=pth->z_table[i+1];
	  else {
	    if (pth->z_visibility_free_streaming == 0.)
	      if ((g < g_previous) && (g <= ppr->visibility_threshold_free_streaming))
		pth->z_visibility_free_streaming=pth->z_table[i];
	  }
      }
    }
    
    /* exp(-kappa) */    
    pth->thermodynamics_table[i*pth->th_size+pth->index_th_exp_m_kappa] = 
      exp(pth->thermodynamics_table[i*pth->th_size+pth->index_th_g]);

    /* compute g' (the plus sign of the second term is correct, see def of -kappa in thermodynamics module!) */
    pth->thermodynamics_table[i*pth->th_size+pth->index_th_dg] = 
      (pth->thermodynamics_table[i*pth->th_size+pth->index_th_ddkappa] +
       pth->thermodynamics_table[i*pth->th_size+pth->index_th_dkappa] *
       pth->thermodynamics_table[i*pth->th_size+pth->index_th_dkappa]) *
      exp(pth->thermodynamics_table[i*pth->th_size+pth->index_th_g]);
    
    /* g */
    pth->thermodynamics_table[i*pth->th_size+pth->index_th_g] = g;

    class_test(pth->thermodynamics_table[i*pth->th_size+pth->index_th_dkappa] == 0.,
	       pth->error_message,
	       "variation rate diverges");

    pth->thermodynamics_table[i*pth->th_size+pth->index_th_rate] =
      sqrt(pow(pth->thermodynamics_table[i*pth->th_size+pth->index_th_dkappa],2.)
	   +pow(pth->thermodynamics_table[i*pth->th_size+pth->index_th_ddkappa]/ 
		pth->thermodynamics_table[i*pth->th_size+pth->index_th_dkappa],2.)
	   +fabs(pth->thermodynamics_table[i*pth->th_size+pth->index_th_dddkappa]/ 
		 pth->thermodynamics_table[i*pth->th_size+pth->index_th_dkappa]));

  }

  /* -> smooth the rate (detaile of smoothing unimportant: only the order of magnitude of the rate matters) */ 
  class_call(array_smooth(pth->thermodynamics_table,
			  pth->th_size,
			  pth->tt_size,
			  pth->index_th_rate,
			  ppr->thermo_rate_smoothing_radius,
			  pth->error_message),
	     pth->error_message,
	     pth->error_message);

  /* check consistency of these values */

  class_test(pth->z_visibility_start_sources == 0.,
	     pth->error_message,
	     "source functions cannot be sampled, probably because ppr->visibility_threshold_start_sources=%e is too large and exceeds maximum of g",ppr->visibility_threshold_start_sources);

  class_test(pth->z_visibility_max == 0.,
	     pth->error_message,
	     "maximum of visibility function is today, hence recombination time could not be found");

  class_test(pth->z_visibility_free_streaming > pth->z_visibility_max,
	     pth->error_message,
	     "pth->z_visibility_free_streaming=%e should never be larger than pth->z_visibility_max=%e",pth->z_visibility_free_streaming,pth->z_visibility_max);

  /** - find conformal recombination time using background_eta_of_z() **/
  class_call(background_eta_of_z(pba,pth->z_visibility_max,&eta_visibility_max),
	     pba->error_message,
	     pth->error_message);

  pth->eta_rec = eta_visibility_max;
  /*   printf("eta_rec=%e\n",eta_visibility_max); */
  /*   printf("eta_0=%e\n",pba->conformal_age-eta_visibility_max); */
  /*   printf("eta_0-eta_rec=%e\n",pba->conformal_age-eta_visibility_max); */

  class_call(background_at_eta(pba,pth->eta_rec, long_info, normal, &last_index_back, pvecback),
	     pba->error_message,
	     pth->error_message);

  pth->rs_rec=pvecback[pba->index_bg_rs];

  /** - fill tables of second derivatives with respect to z (in view of spline interpolation) */
  class_call(array_spline_table_lines(pth->z_table,
				      pth->tt_size,
				      pth->thermodynamics_table,
				      pth->th_size,
				      pth->d2thermodynamics_dz2_table,
				      _SPLINE_EST_DERIV_,
				      pth->error_message),
	     pth->error_message,
	     pth->error_message);

  if (pth->thermodynamics_verbose > 0) {
    printf(" -> recombination at z = %f\n",pth->z_visibility_max);
    if (pth->reio_parametrization != reio_none) {
      if (pth->reio_z_or_tau==reio_tau)
	printf(" -> reionization  at z = %f\n",pth->z_reio);
      if (pth->reio_z_or_tau==reio_z)
	printf(" -> reionization with optical depth = %f\n",pth->tau_reio);
    }
  }

  free(pvecback);

  return _SUCCESS_;
}

/**
 * Free all memory space allocated by thermodynamics_init().
 * 
 * To be called at the end of each run, only when no further calls to
 * thermodynamics_at_z() are needed.
 *
 * @return the error status
 */
int thermodynamics_free(
			struct thermo * pth
			) {

  free(pth->z_table);
  free(pth->thermodynamics_table);
  free(pth->d2thermodynamics_dz2_table);

  return _SUCCESS_;
}

/**
 * Assign value to each relevant index in vectors of thermodynamical quantities.
 *
 * Called once by thermodynamics_init().
 *
 * @param Input/Output: pointer to thermo structure
 * @param Input/Output: pointer to recombination structure
 * @param Input/Output: pointer to reionization structure 
 * @return the error status
 */
int thermodynamics_indices(
			   struct thermo * pth,
			   struct recombination * preco,
			   struct reionization * preio
			   ) {

  /** Summary: */

  /** - define local variables */

  /* a running index for the vector of thermodynamics quantities */
  int index;

  /** - intialization of all indices and flags in thermo structure */
  index = 0;
  
  pth->index_th_xe = index;
  index++;
  pth->index_th_dkappa = index;
  index++;
  pth->index_th_ddkappa = index;
  index++;
  pth->index_th_dddkappa = index;
  index++;
  pth->index_th_exp_m_kappa = index;
  index++;
  pth->index_th_g = index;
  index++;
  pth->index_th_dg = index;
  index++;
  pth->index_th_Tb = index;
  index++;
  pth->index_th_cb2 = index;
  index++;
  pth->index_th_rate = index;
  index++;

  /* end of indices */
  pth->th_size = index;

  /** - intialization of all indices and flags in recombination structure */
  index = 0;

  preco->index_re_z = index;
  index++;
  preco->index_re_xe = index;
  index++;
  preco->index_re_dkappadeta = index;
  index++;
  preco->index_re_Tb = index;
  index++;
  preco->index_re_cb2 = index;
  index++;

  /* end of indices */
  preco->re_size = index;

  /** - intialization of all indices and flags in reionization structure */
  index = 0;

  preio->index_re_z = index;
  index++;
  preio->index_re_xe = index;
  index++;
  preio->index_re_Tb = index;
  index++;
  preio->index_re_cb2 = index;
  index++;
  preio->index_re_dkappadeta = index;
  index++;
  preio->index_re_dkappadz = index;
  index++;
  preio->index_re_d3kappadz3 = index;
  index++;

  /* end of indices */
  preio->re_size = index;

  /* same with parameters of the function x_e(z) */
  
  /* case wher x_e(z) taken like in CAMB (other cases can be added) */ 
  if (pth->reio_parametrization == reio_camb) {

    index=0;

    preio->index_reio_redshift = index;
    index++;
    preio->index_reio_start = index;
    index++;
    preio->index_reio_xe_before = index;
    index++;
    preio->index_reio_xe_after = index;
    index++;
    preio->index_reio_exponent = index;
    index++;
    preio->index_reio_width = index;
    index++;
    preio->index_helium_fullreio_fraction = index;
    index++;
    preio->index_helium_fullreio_redshift = index;
    index++;
    preio->index_helium_fullreio_width = index;
    index++;

    preio->reio_num_params = index;

  }

  return _SUCCESS_;
}

/**
 * This subroutine contains the reionization function \f$ X_e(z) \f$
 *
 * @param z Input: redshift
 * @param Input: pointer to reionization structure, containing the parameters of the function \f$ X_e(z) \f$
 * @param xe Ouput: \f$ X_e(z) \f$
 */
int thermodynamics_reionization_function(
					 double z,
					 struct thermo * pth,
					 struct reionization * preio,
					 double * xe
					 ) {

  double argument;

  /* ionization function like in CAMB */
  if (pth->reio_parametrization == reio_camb) {

    if (z > preio->reionization_parameters[preio->index_reio_start]) {
      *xe = preio->reionization_parameters[preio->index_reio_xe_before];
    }
    else {
      argument = (pow((1.+preio->reionization_parameters[preio->index_reio_redshift]),
		      preio->reionization_parameters[preio->index_reio_exponent]) 
		  - pow((1.+z),preio->reionization_parameters[preio->index_reio_exponent]))
	/(preio->reionization_parameters[preio->index_reio_exponent] 
	  /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
	  *pow((1.+preio->reionization_parameters[preio->index_reio_redshift]),
	       (preio->reionization_parameters[preio->index_reio_exponent]-1.)))
	/preio->reionization_parameters[preio->index_reio_width]; 
      /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */

      *xe = (preio->reionization_parameters[preio->index_reio_xe_after]
	     -preio->reionization_parameters[preio->index_reio_xe_before])
	*(tanh(argument)+1.)/2.
	+preio->reionization_parameters[preio->index_reio_xe_before];

    }
    
    argument = (preio->reionization_parameters[preio->index_helium_fullreio_redshift] - z)
      /preio->reionization_parameters[preio->index_helium_fullreio_width];
    /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
    *xe += preio->reionization_parameters[preio->index_helium_fullreio_fraction] 
      * (tanh(argument)+1.)/2.;

    return _SUCCESS_;

  }

  class_test(0 == 0,
	     pth->error_message,
	     "value of reio_parametrization=%d unclear",pth->reio_parametrization);
}


int thermodynamics_get_xe_before_reionization(
					      struct precision * ppr,
					      struct thermo * pth,
					      struct recombination * preco,
					      double z,
					      double * xe
					      ) {

  int i;

  i=0;
  while (preco->recombination_table[i*preco->re_size+preco->index_re_z] < z) {
    i++;
    class_test(i == ppr->recfast_Nz0,
	       pth->error_message,
	       "z = %e > largest redshift in thermodynamics table \n",ppr->reionization_z_start_max);
  }
  
  *xe = preco->recombination_table[i*preco->re_size+preco->index_re_xe];

  return _SUCCESS_;

}
  

/** 
 * Add reionization to recombination history.
 *
 * Called once by thermodynamics_init().
 * 
 * @param Input/Output: pointer to recombination structure
 * @param Input/Output: pointer to reionization structure 
 * @return the error status
 */
int thermodynamics_reionization(
				struct precision * ppr,
				struct background * pba,
				struct thermo * pth,
				struct recombination * preco,
				struct reionization * preio,
				double * pvecback
				) {

  int i,counter;
  double z_sup,z_mid,z_inf;
  double tau_sup,tau_mid,tau_inf;

  /** (a) if reionization implemented like in CAMB */

  if (pth->reio_parametrization == reio_camb) {

    /** - allocate the parameters of the function \f$ X_e(z) \f$ */

    class_alloc(preio->reionization_parameters,preio->reio_num_params*sizeof(double),pth->error_message); 
    
    /** - set values of these parameters, excepted those depending on the reionization redshift */

    /* preio->reionization_parameters[preio->index_reio_redshift]: reionization redshift to be found later given the optical depth */
    /* preio->reionization_parameters[preio->index_reio_start]: starting redshift to be found later given the optical depth */
    /* preio->reionization_parameters[preio->index_reio_xe_before]: initial X_e to be found later given the optical depth */
    preio->reionization_parameters[preio->index_reio_xe_after] = 1. + pth->YHe/(_not4_*(1.-pth->YHe));    /* xe_after_reio: H + singly ionized He (note: segmentation fault impossible, checked before that denominator is non-zero) */
    preio->reionization_parameters[preio->index_reio_exponent] = ppr->reionization_exponent; /* reio_exponent */
    preio->reionization_parameters[preio->index_reio_width] = ppr->reionization_width;    /* reio_width */
    preio->reionization_parameters[preio->index_helium_fullreio_fraction] = pth->YHe/(_not4_*(1.-pth->YHe)); /* helium_fullreio_fraction (note: segmentation fault impossible, checked before that denominator is non-zero) */
    preio->reionization_parameters[preio->index_helium_fullreio_redshift] = ppr->helium_fullreio_redshift; /* helium_fullreio_redshift */
    preio->reionization_parameters[preio->index_helium_fullreio_width] = ppr->helium_fullreio_width;    /* helium_fullreio_width */

    class_test(preio->reionization_parameters[preio->index_reio_exponent]==0,
	       pth->error_message,
	       "stop to avoid division by zero");

    class_test(preio->reionization_parameters[preio->index_reio_width]==0,
	       pth->error_message,
	       "stop to avoid division by zero");

    class_test(preio->reionization_parameters[preio->index_helium_fullreio_width]==0,
	       pth->error_message,
	       "stop to avoid division by zero");

    /** - if reionization redshift given as an input, initialize the remaining values and fill reionization table*/

    if (pth->reio_z_or_tau == reio_z) {
      
      /* reionization redshift */
      preio->reionization_parameters[preio->index_reio_redshift] = pth->z_reio; 
      /* infer starting redshift */
      preio->reionization_parameters[preio->index_reio_start] = preio->reionization_parameters[preio->index_reio_redshift]+ppr->reionization_start_factor*ppr->reionization_width;
      class_test(preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max,
		 pth->error_message,
		 "starting redshift for reionization > reionization_z_start_max = %e\n",ppr->reionization_z_start_max);

      /* infer xe_before_reio */
      class_call(thermodynamics_get_xe_before_reionization(ppr,
							   pth,
							   preco,
							   preio->reionization_parameters[preio->index_reio_redshift],
							   &(preio->reionization_parameters[preio->index_reio_xe_before])),
		 pth->error_message,
		 pth->error_message);

      /* fill reionization table */
      class_call(thermodynamics_reionization_discretize(ppr,pba,pth,preco,preio,pvecback),
		 pth->error_message,
		 pth->error_message);

      pth->tau_reio=preio->reionization_optical_depth;

    }

    /** - if reionization optical depth given as an input, find reionization redshift by dichotomy and initialize the remaining values */

    if (pth->reio_z_or_tau == reio_tau) {

      /* upper value */

      z_sup = ppr->reionization_z_start_max-ppr->reionization_start_factor*ppr->reionization_width;
      class_test(z_sup < 0.,
		 pth->error_message,
		 "parameters are such that reionization cannot take place before today while starting after z_start_max; need to increase z_start_max");

      /* reionization redshift */
      preio->reionization_parameters[preio->index_reio_redshift] = z_sup; 
      /* infer starting redshift */
      preio->reionization_parameters[preio->index_reio_start] = ppr->reionization_z_start_max;
      /* infer xe_before_reio */
      class_call(thermodynamics_get_xe_before_reionization(ppr,
							   pth,
							   preco,
							   preio->reionization_parameters[preio->index_reio_redshift],
							   &(preio->reionization_parameters[preio->index_reio_xe_before])),
		 pth->error_message,
		 pth->error_message);

      /* fill reionization table */
      class_call(thermodynamics_reionization_discretize(ppr,pba,pth,preco,preio,pvecback),
		 pth->error_message,
		 pth->error_message);

      tau_sup=preio->reionization_optical_depth;

      class_test(tau_sup < pth->tau_reio,
		 pth->error_message,
		 "parameters are such that reionization cannot start after z_start_max");

      /* lower value */

      z_inf = 0.;
      tau_inf = 0.;      

      /* try intermediate values */
    
      counter=0;
      while ((tau_sup-tau_inf) > pth->tau_reio * ppr->reionization_optical_depth_tol) {
	z_mid=0.5*(z_sup+z_inf);
      
	/* reionization redshift */
	preio->reionization_parameters[preio->index_reio_redshift] = z_mid;
	/* infer starting redshift */
	preio->reionization_parameters[preio->index_reio_start] = preio->reionization_parameters[preio->index_reio_redshift]+ppr->reionization_start_factor*ppr->reionization_width;
	class_test(preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max,
		   pth->error_message,
		   "starting redshift for reionization > reionization_z_start_max = %e",ppr->reionization_z_start_max);

	/* infer xe_before_reio */
	class_call(thermodynamics_get_xe_before_reionization(ppr,
							     pth,
							     preco,
							     preio->reionization_parameters[preio->index_reio_redshift],
							     &(preio->reionization_parameters[preio->index_reio_xe_before])),
		   pth->error_message,
		   pth->error_message);

	/* clean and fill reionization table */
	free(preio->reionization_table);
	class_call(thermodynamics_reionization_discretize(ppr,pba,pth,preco,preio,pvecback),
		   pth->error_message,
		   pth->error_message);

	tau_mid=preio->reionization_optical_depth;
	
	/* trial */

	if (tau_mid > pth->tau_reio) {
	  z_sup=z_mid;
	  tau_sup=tau_mid;
	}
	else {
	  z_inf=z_mid;
	  tau_inf=tau_mid;
	}

	counter++;
	class_test(counter > _MAX_IT_,
		   pth->error_message,
		   "while searching for reionization_optical_depth, maximum number of iterations exceeded");
      }

      /* store z reionization in thermodynamics structure */
      pth->z_reio=preio->reionization_parameters[preio->index_reio_redshift];

    }

    free(preio->reionization_parameters);

    return _SUCCESS_;

  }

  class_test(0 == 0,
	     pth->error_message,
	     "value of reio_z_or_tau=%d unclear",pth->reio_z_or_tau);

}

int thermodynamics_reionization_discretize(
					   struct precision * ppr,
					   struct background * pba,
					   struct thermo * pth,
					   struct recombination * preco,
					   struct reionization * preio,
					   double * pvecback
					   ) {

  /* a growing table (since the number of redshift steps is not known a priori) */
  growTable gTable;
  /* needed for growing table */
  double * pData;
  /* needed for growing table */
  void * memcopy_result;
  /* current vector of values related to reionization */
  double * reio_vector;
  /* running index inside thermodynamics table */
  int i;
  int number_of_redshifts;
  /* values of z, dz, X_e */
  double dz,dz_max;
  double z,z_next;
  double xe,xe_next;
  double dkappadz,dkappadz_next;
  double Tb,Tba2,Yp;
  double dkappadeta,dkappadeta_next;

  Yp = pth->YHe;

  /** (a) allocate vector of values related to reionization */
  class_alloc(reio_vector,preio->re_size*sizeof(double),pth->error_message);

  /** (b) create a growTable with gt_init() */
  class_call(gt_init(&gTable),
	     gTable.error_message,
	     pth->error_message);

  /** (c) first line is taken from thermodynamics table, just before reionization starts */

  /** - look where to start in current thermodynamics table */
  i=0;
  while (preco->recombination_table[i*preco->re_size+preco->index_re_z] < preio->reionization_parameters[preio->index_reio_start]) {
    i++;
    class_test(i == ppr->recfast_Nz0,
	       pth->error_message,
	       "reionization_z_start_max = %e > largest redshift in thermodynamics table",ppr->reionization_z_start_max);
  }

  /** - get redshift */
  z=preco->recombination_table[i*preco->re_size+preco->index_re_z];
  reio_vector[preio->index_re_z]=z;
  preio->index_reco_when_reio_start=i;

  /** - get \f$ X_e \f$ */
  class_call(thermodynamics_reionization_function(z,pth,preio,&xe),
	     pth->error_message,
	     pth->error_message);

  reio_vector[preio->index_re_xe] = xe;

  /** - get \f$ d kappa / d z = (d kappa / d eta) * (d eta / d z) = - (d kappa / d eta) / H \f$ */
  class_call(background_functions(pba,1./(1.+z),short_info,pvecback),
	     pba->error_message,
	     pth->error_message);

  reio_vector[preio->index_re_dkappadeta] = (1.+z) * (1.+z) * pth->n_e * xe * _sigma_ * _Mpc_over_m_;

  class_test(pvecback[pba->index_bg_H] == 0.,
	     pth->error_message,
	     "stop to avoid division by zero");

  reio_vector[preio->index_re_dkappadz] = reio_vector[preio->index_re_dkappadeta] / pvecback[pba->index_bg_H];

  dkappadz = reio_vector[preio->index_re_dkappadz];
  dkappadeta = reio_vector[preio->index_re_dkappadeta];

  /** - get baryon temperature **/
  Tb = preco->recombination_table[i*preco->re_size+preco->index_re_Tb];
  reio_vector[preio->index_re_Tb] = Tb;

  /** - after recombination, Tb scales like (1+z)**2. Compute constant factor Tb/(1+z)**2. */
  Tba2 = Tb/(1+z)/(1+z);

  /** - get baryon sound speed */
  reio_vector[preio->index_re_cb2] = 5./3. * _k_B_ / ( _c_ * _c_ * _m_H_ ) * (1. + (1./_not4_ - 1.) * Yp + xe * (1.-Yp)) * Tb;

  /** - store these values in growing table */
  class_call(gt_add(&gTable,_GT_END_,(void *) reio_vector,sizeof(double)*(preio->re_size)),
	     gTable.error_message,
	     pth->error_message);

  number_of_redshifts=1;

  /** (d) set the maximum step value (equal to the step in thermodynamics table) */
  dz_max=preco->recombination_table[i*preco->re_size+preco->index_re_z]
    -preco->recombination_table[(i-1)*preco->re_size+preco->index_re_z];

  /** (e) loop over redshift values */
  while (z > 0.) {

    /** - try default step */
    dz = dz_max;
    z_next=z-dz;
    if (z_next < 0.) z_next=0.;

    class_call(thermodynamics_reionization_function(z_next,pth,preio,&xe_next),
	       pth->error_message,
	       pth->error_message);

    class_call(background_functions(pba,1./(1.+z_next),short_info,pvecback),
	       pba->error_message,
	       pth->error_message);

    class_test(pvecback[pba->index_bg_H] == 0.,
	       pth->error_message,
	       "stop to avoid division by zero");

    dkappadz_next= (1.+z_next) * (1.+z_next) * pth->n_e * xe_next * _sigma_ * _Mpc_over_m_ / pvecback[pba->index_bg_H];

    dkappadeta_next= (1.+z_next) * (1.+z_next) * pth->n_e * xe_next * _sigma_ * _Mpc_over_m_;

    class_test((dkappadz == 0.) || (dkappadeta == 0.),
	       pth->error_message,
	       "dkappadz=%e, dkappadeta=%e, stop to avoid division by zero",dkappadz,dkappadeta);

    /** - reduce step if necessary */
    while (((fabs(dkappadz_next-dkappadz)/dkappadz) > ppr->reionization_sampling) || 
	   ((fabs(dkappadeta_next-dkappadeta)/dkappadeta) > ppr->reionization_sampling)) {

      dz*=0.9;

      class_test(dz < ppr->smallest_allowed_variation,
		 pth->error_message,
		 "integration step =%e < machine precision : leads either to numerical error or infinite loop",dz);

      z_next=z-dz;

      class_call(thermodynamics_reionization_function(z_next,pth,preio,&xe_next),
		 pth->error_message,
		 pth->error_message);

      class_call(background_functions(pba,1./(1.+z_next),short_info,pvecback),
		 pba->error_message,
		 pth->error_message);

      class_test(pvecback[pba->index_bg_H] == 0.,
		 pth->error_message,
		 "stop to avoid division by zero");

      dkappadz_next= (1.+z_next) * (1.+z_next) * pth->n_e * xe_next * _sigma_ * _Mpc_over_m_ / pvecback[pba->index_bg_H];

      dkappadeta_next= (1.+z_next) * (1.+z_next) * pth->n_e * xe_next * _sigma_ * _Mpc_over_m_;
    }

    /** - get \f$ z, X_e, d kappa / d z \f$ and store in growing table */
    z=z_next;
    xe=xe_next;
    dkappadz=dkappadz_next;
    dkappadeta= dkappadeta_next;

    class_test((dkappadz == 0.) || (dkappadeta == 0.),
	       pth->error_message,
	       "dkappadz=%e, dkappadeta=%e, stop to avoid division by zero",dkappadz,dkappadeta);

    reio_vector[preio->index_re_z] = z;   
    reio_vector[preio->index_re_xe] = xe;
    reio_vector[preio->index_re_dkappadz] = dkappadz;
    reio_vector[preio->index_re_dkappadeta] = dkappadz * pvecback[pba->index_bg_H];

    /** - get baryon temperature **/
    Tb = Tba2*(1+z)*(1+z);
    reio_vector[preio->index_re_Tb] = Tb;
    
    /** - get baryon sound speed */
    reio_vector[preio->index_re_cb2] = 5./3. * _k_B_ / ( _c_ * _c_ * _m_H_ ) * (1. + (1./_not4_ - 1.) * Yp + xe * (1.-Yp)) * Tb;

    class_call(gt_add(&gTable,_GT_END_,(void *) reio_vector,sizeof(double)*(preio->re_size)),
	       gTable.error_message,
	       pth->error_message);

    number_of_redshifts++;
  }
  
  /** - allocate reionization_table with correct size */
  class_alloc(preio->reionization_table,preio->re_size*number_of_redshifts*sizeof(double),pth->error_message);

  preio->rt_size=number_of_redshifts;

  /** - retrieve data stored in the growTable with gt_getPtr() */
  class_call(gt_getPtr(&gTable,(void**)&pData),
	     gTable.error_message,
	     pth->error_message);

  /** -- copy growTable to reionization_temporary_table (invert order of lines, so that redshift is growing, like in recombination table) */
  for (i=0; i < preio->rt_size; i++) {
    memcopy_result = memcpy(preio->reionization_table+i*preio->re_size,pData+(preio->rt_size-i-1)*preio->re_size,preio->re_size*sizeof(double));
    class_test(memcopy_result != preio->reionization_table+i*preio->re_size,
	       pth->error_message,
	       "cannot copy data back to reionization_temporary_table");

  }

  /** - free the growTable with gt_free() , free vector of reionization variables */
  class_call(gt_free(&gTable),
	     gTable.error_message,
	     pth->error_message);
  
  free(reio_vector);

  /** - spline \f$ d tau / dz \f$ with respect to z in view of integrating for optical depth */
  class_call(array_spline(preio->reionization_table,
			  preio->re_size,
			  preio->rt_size,
			  preio->index_re_z,
			  preio->index_re_dkappadz,
			  preio->index_re_d3kappadz3,
			  _SPLINE_EST_DERIV_,
			  pth->error_message),
	     pth->error_message,
	     pth->error_message);
  
  /** - integrate for optical depth */
  class_call(array_integrate_all_spline(preio->reionization_table,
					preio->re_size,
					preio->rt_size,
					preio->index_re_z,
					preio->index_re_dkappadz,
					preio->index_re_d3kappadz3,
					&(preio->reionization_optical_depth),
					pth->error_message),
	     pth->error_message,
	     pth->error_message);

  return _SUCCESS_;

}

/** 
 * Integrate thermodynamics with RECFAST.
 *
 * Integrate thermodynamics with RECFAST, allocate and fill the part
 * of the thermodynamics interpolation table (the rest is filled in
 * thermodynamics_init()). Called once by
 * thermodynamics_init().
 *
 *******************************************************************************
 * RECFAST is an integrator for Cosmic Recombination of Hydrogen and Helium,   *
 * developed by Douglas Scott (dscott@astro.ubc.ca)                            *
 * based on calculations in the paper Seager, Sasselov & Scott                 *
 * (ApJ, 523, L1, 1999).                                                       *
 * and "fudge" updates in Wong, Moss & Scott (2008).                           *
 *                                                                             *
 * Permission to use, copy, modify and distribute without fee or royalty at    *
 * any tier, this software and its documentation, for any purpose and without  *
 * fee or royalty is hereby granted, provided that you agree to comply with    *
 * the following copyright notice and statements, including the disclaimer,    *
 * and that the same appear on ALL copies of the software and documentation,   *
 * including modifications that you make for internal use or for distribution: *
 *                                                                             *
 * Copyright 1999-2010 by University of British Columbia.  All rights reserved.*
 *                                                                             *
 * THIS SOFTWARE IS PROVIDED "AS IS", AND U.B.C. MAKES NO                      *
 * REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.                          *
 * BY WAY OF EXAMPLE, BUT NOT LIMITATION,                                      *
 * U.B.C. MAKES NO REPRESENTATIONS OR WARRANTIES OF                            *
 * MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT               *
 * THE USE OF THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE         *
 * ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.            *
 *******************************************************************************
 *
 * Version 1.5: includes extra fitting function from
 *              Rubino-Martin et al. arXiv:0910.4383v1 [astro-ph.CO]
 *
 * @param Input/Ouput: pointer to recombination structure
 * @return the error status
 */
int thermodynamics_recombination(
				 struct precision * ppr,
				 struct background * pba,
				 struct thermo * pth,
				 struct recombination * preco,
				 double * pvecback
				 ) {

  /* vector of variables to be integrated: xH, xHe, Tmat */
  double y[3],dy[3];

  double OmegaB,Yp,zinitial,zfinal,x_H0,x_He0;
  double z,mu_H,n,Lalpha,Lalpha_He,DeltaB,DeltaB_He,mu_T;
  double x0,w0,w1,Lw0,Lw1,hW;
  double zstart,zend,rhs,Trad,Tmat;
  int i,Nz;

  /* contains all quantities relevant for the integration algorithm */
  struct generic_integrator_workspace gi;
  /* contains all fixed parameters which should be passed to thermodynamics_derivs_with_recfast */
  struct thermodynamics_parameters_and_workspace tpaw;
  
  /** Summary: */

  /** - allocate memory for thermodynamics interpolation tables (size known in advance) */
  preco->rt_size = ppr->recfast_Nz0;
  class_alloc(preco->recombination_table,preco->re_size*preco->rt_size*sizeof(double),pth->error_message);

  /** - initialize generic integrator with initialize_generic_integrator() */
  class_call(initialize_generic_integrator(_RECFAST_INTEG_SIZE_, &gi),
	     gi.error_message,
	     pth->error_message);
  
  /** - read a few precision/cosmological parameters */

  /* Nz */
  Nz=ppr->recfast_Nz0;

  /* Ho : must be h in 100km/s/Mpc * bigH */
  preco->H0 = pba->H0 * _bigH_;

  /* Omega_b */
  OmegaB = pba->Omega0_b;
  /*printf("Omega_b = %f \n",OmegaB);*/

  /* Yp */
  Yp = pth->YHe;
  /*printf("Y_He = %f \n",Yp);*/

  /* Tnow */
  preco->Tnow = pth->Tcmb;

  /* z_initial and z_final */
  zinitial=ppr->recfast_z_initial;
  zfinal=ppr->recfast_z_final;

  /* H_frac */ 
  preco->H_frac = ppr->recfast_H_frac;

  /* H fudging */
  class_test((ppr->recfast_Hswitch != _TRUE_) && (ppr->recfast_Hswitch != _FALSE_),
	     pth->error_message,
	     "RECFAST error: unknown H fudging scheme");
  preco->fu = ppr->recfast_fudge_H;
  if (ppr->recfast_Hswitch == _TRUE_) 
    preco->fu += ppr->recfast_delta_fudge_H;

  /* He fudging */
  class_test((ppr->recfast_Heswitch < 0) || (ppr->recfast_Heswitch > 6),
	     pth->error_message,
	     "RECFAST error: unknown He fudging scheme");

  /* related quantities */ 
  z=zinitial;
  mu_H = 1./(1.-Yp);
  mu_T = _not4_ /(_not4_ - (_not4_-1.)*Yp); /* recfast 1.4*/
  preco->fHe = Yp/(_not4_ *(1.-Yp)); /* recfast 1.4 */
  preco->Nnow = 3.*preco->H0*preco->H0*OmegaB/(8.*_PI_*_G_*mu_H*_m_H_);
  pth->n_e = preco->Nnow;
  /*  printf("Nnow= %e\n",Nnow); */

  /* quantities related to constants defined in thermodynamics.h */
  n = preco->Nnow * pow((1.+z),3);
  Lalpha = 1./_L_H_alpha_;
  Lalpha_He = 1./_L_He_2p_;
  DeltaB = _h_P_*_c_*(_L_H_ion_-_L_H_alpha_);
  preco->CDB = DeltaB/_k_B_;
  DeltaB_He = _h_P_*_c_*(_L_He1_ion_-_L_He_2s_);
  preco->CDB_He = DeltaB_He/_k_B_;
  preco->CB1 = _h_P_*_c_*_L_H_ion_/_k_B_;
  preco->CB1_He1 = _h_P_*_c_*_L_He1_ion_/_k_B_;
  preco->CB1_He2 = _h_P_*_c_*_L_He2_ion_/_k_B_;
  preco->CR = 2.*_PI_*(_m_e_/_h_P_)*(_k_B_/_h_P_);
  preco->CK = pow(Lalpha,3)/(8.*_PI_);
  preco->CK_He = pow(Lalpha_He,3)/(8.*_PI_);
  preco->CL = _c_*_h_P_/(_k_B_*Lalpha);
  preco->CL_He = _c_*_h_P_/(_k_B_/_L_He_2s_);
  preco->CT = (8./3.)*(_sigma_/(_m_e_*_c_))*_a_;
  preco->Bfact = _h_P_*_c_*(_L_He_2p_-_L_He_2s_)/_k_B_;

  /* C1P3P = _C2p1P_-_C2p3P_; */
  /* cc3P1P = _C2p3P_/_C2p1P_;  */
  /* cccP = pow(cc3P1P,3); */
  /* hck=_h_P_*_c_/_k_B_;  */

  tpaw.pba = pba;
  tpaw.ppr = ppr;
  tpaw.preco = preco;
  tpaw.pvecback = pvecback;

  class_test(zinitial < 8000.,
	     pth->error_message,
	     "z_initial=%f<8000, should get recfast initial conditions from get_init()",zinitial);

  /* uncomment this part only if the initial redshift is not larger than 8000 */
  /*   if (get_init(z,y)== _FAILURE_) { */
  /*     sprintf(pth->error_message,"%s(L:%d): error in calling get_init()",__func__,__LINE__); */
  /*     return _FAILURE_; */
  /*   } */

  /** - impose initial conditions */
  y[0] = 1.;
  y[1] = 1.;
  x0 = 1.+2.*preco->fHe;
  y[2] = preco->Tnow*(1.+z);

  w0=1./ sqrt(1. + zinitial);
  w1=1./ sqrt(1. + zfinal);
  Lw0 = log(w0);
  Lw1 = log(w1);
  hW=(Lw1-Lw0)/(1.*Nz);

  /** - loop over redshift steps Nz : perform each step with
      generic_integrator(), store the results in the table using thermodynamics_derivs_with_recfast()*/

  for(i=0; i <Nz; i++) {

    zstart = zinitial  + (zfinal-zinitial)*i/(1.*Nz);
    zend   = zinitial  + (zfinal-zinitial)*(i+1)/(1.*Nz);
    z = zend;
    
    if (zend > 8000.) {
      x_H0 = 1.;
      x_He0 = 1.;
      x0 = 1.+2.*preco->fHe;
      y[0] = x_H0;
      y[1] = x_He0;
      y[2] = preco->Tnow*(1.+z);
    }
    else 
      if (z > 5000.) {
	x_H0 = 1.;
	x_He0 = 1.;
	rhs = exp( 1.5*log(preco->CR*preco->Tnow/(1.+z)) - preco->CB1_He2/(preco->Tnow*(1.+z)) ) / preco->Nnow;
	rhs = rhs*1.;
	x0 = 0.5*(sqrt(pow((rhs-1.-preco->fHe),2) + 4.*(1.+2.*preco->fHe)*rhs) - (rhs-1.-preco->fHe));
	y[0] = x_H0;
	y[1] = x_He0;
	y[2] = preco->Tnow*(1.+z);
      }
      else 
	if (z > 3500.) {
	  x_H0 = 1.;
	  x_He0 = 1.;
	  x0 = x_H0 + preco->fHe*x_He0;
	  y[0] = x_H0;
	  y[1] = x_He0;
	  y[2] = preco->Tnow*(1.+z);
	}
	else 
	  if (y[1] > ppr->recfast_x_He0_trigger) {
	    x_H0 = 1.;
	    rhs = exp(1.5*log(preco->CR*preco->Tnow/(1.+z)) - preco->CB1_He1/(preco->Tnow*(1.+z)))/preco->Nnow;
	    rhs = rhs*4.;
	    x_He0 = 0.5*(sqrt(pow((rhs-1.),2) + 4.*(1.+preco->fHe)*rhs )- (rhs-1.));
	    x0 = x_He0;
	    x_He0 = (x0-1.)/preco->fHe;
	    y[0] = x_H0;
	    y[1] = x_He0;
	    y[2] = preco->Tnow*(1.+z);
	  }
	  else 
	    if (y[0] > ppr->recfast_x_H0_trigger) {
	      rhs = exp(1.5*log(preco->CR*preco->Tnow/(1.+z)) - preco->CB1/(preco->Tnow*(1.+z)))/preco->Nnow;
	      x_H0 = 0.5*(sqrt(pow(rhs,2)+4.*rhs) - rhs);

	      class_call(generic_integrator(thermodynamics_derivs_with_recfast,
					    zstart,
					    zend,
					    y,
					    &tpaw,
					    ppr->tol_thermo_integration,
					    ppr->smallest_allowed_variation,
					    &gi),
			 gi.error_message,
			 pth->error_message);

	      y[0] = x_H0;
	      x0 = y[0] + preco->fHe*y[1];

	      

	    }      
	    else {

	      class_call(generic_integrator(thermodynamics_derivs_with_recfast,
					    zstart,
					    zend,
					    y,
					    &tpaw,
					    ppr->tol_thermo_integration,
					    ppr->smallest_allowed_variation,
					    &gi),
			 gi.error_message,
			 pth->error_message);

	      x0 = y[0] + preco->fHe*y[1];

	    }

    /* store the results in the table */
    /* results are obtained in order of decreasing z, and stored in order of growing z */

    /* -> redshift */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_z)=zend;

    /* -> ionization fraction */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_xe)=x0;

    /* -> Tb */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_Tb)=y[2];

    /* -> get dTb/dz=dy[2] */
    class_call(thermodynamics_derivs_with_recfast(zend, y, dy, &tpaw,pth->error_message),
	       pth->error_message,
	       pth->error_message);

    /* -> cb2 = (k_B/mu) Tb (1-1/3 dlnTb/dlna) = (k_B/mu) Tb (1+1/3 (1+z) dlnTb/dz) */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_cb2)
      = _k_B_ / ( _c_ * _c_ * _m_H_ ) * (1. + (1./_not4_ - 1.) * Yp + x0 * (1.-Yp)) * y[2] * (1. + (1.+zend) * dy[2] / y[2] / 3.);

    /* -> dkappa/deta = a n_e x_e sigma_T = a^{-2} n_e(today) x_e sigma_T (in units of 1/Mpc) */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_dkappadeta)
      = (1.+zend) * (1.+zend) * preco->Nnow * x0 * _sigma_ * _Mpc_over_m_;

/*     fprintf(stdout,"%g %g %g %g %g\n", */
/* 	    zend, */
/* 	    y[0], */
/* 	    y[1], */
/* 	    y[2], */
/* 	    x0); */
    
  }

/*   class_test(0==0,pth->error_message,"stop here for testing") */


  /* RECFAST is done */

  /** - cleanup generic integrator with cleanup_generic_integrator() */

  class_call(cleanup_generic_integrator(&gi),
	     gi.error_message,
	     pth->error_message);
  
  return _SUCCESS_;
}

/**
 * Smoothing routine for reducing numerical errors arising from
 * discontinuities in derivatives.
 *
 * Called by thermodynamics_init() when computing the derivative of the
 * Thomson scattering rate.
 *
 * @param index_th index of the column which needs smoothing in the thermodynamics table 
 @return the error status
*/
/* int thermodynamics_cure_discontinuity(int index_th) { */

/*   int i, i_start, i_stop; */

/*   double z_previous,z_current; */
/*   double z_left,z_right; */
/*   double y_left,y_right; */
/*   double weigth; */

/*   i_start = 2; */
/*   i_stop = pth->tt_size - 2; */

/*   z_current = pth->z_table[i_start-1]; */

/*   for (i=i_start; i<= i_stop; i++) { */
    
/*     z_previous=z_current; */
/*     z_current=pth->z_table[i-1]; */

/*     if ( ((z_current > 5000.) && (z_previous <= 5000.)) ||  */
/* 	 ((z_current > 8000.) && (z_previous <= 8000.)) ||  */
/* 	 ((z_current > 3500.) && (z_previous <= 3500.))  */
/* 	 ){ */

/*       y_left = *(pth->thermodynamics_table+(i-2)*pth->th_size+index_th); */
/*       y_right = *(pth->thermodynamics_table+(i+1)*pth->th_size+index_th); */
	
/*       z_left = pth->z_table[i-2]; */
/*       z_right = pth->z_table[i+1]; */
	
/*       weigth = (z_previous - z_left)/(z_right - z_left); */
	
/*       *(pth->thermodynamics_table+(i-1)*pth->th_size+index_th) = y_left  */
/* 	+ (y_right - y_left) * weigth;  */
	
/*       weigth = (z_current - z_left)/(z_right - z_left); */
	
/*       *(pth->thermodynamics_table+i*pth->th_size+index_th) = y_left  */
/* 	+ (y_right - y_left) * weigth;  */

/*     } */
/*   } */
  
/*   return _SUCCESS_; */
  
/* } */
  


/**
 * Subroutine evaluating the derivative with respect to redshift of thermodynamical quantities (from RECFAST version 1.4). 
 *
 * Computes derivatives of the three variables to integrate: 
 * \f$ d x_H / dz, d x_{He} / dz, d T_{mat} / dz \f$.
 * 
 * This is one of the few functions in the code which are passed to the generic_integrator() routine. 
 * Since generic_integrator() should work with functions passed from various modules, the format of the arguments
 * is a bit special:
 * - fixed parameters and workspaces are passed through a generic pointer. Here, this pointer contains the 
 *   precision, background and recombination structures, plus a background vector, but generic_integrator() 
 *   doesn't know its fine structure.
 * - the error management is a bit special: errors are not written as usual to pth->error_message, but to a generic 
 *   error_message passed in the list of arguments.
 *
 * @param z Input : redshift
 * @param y Input : vector of variable
 * @param dy Output : its derivative (already allocated)
 * @param fixed_parameters Input: pointer to fixed parameters (e.g. indices)
 * @param error_message Output : error message
 */
int thermodynamics_derivs_with_recfast(
				       double z,
				       double * y,
				       double * dy,
				       void * parameters_and_workspace,
				       ErrorMsg error_message
				       ) {

  double x,n,n_He,Trad,Tmat,x_H,x_He,Hz,dHdz,epsilon;
  double Rup,Rdown,K,K_He,Rup_He,Rdown_He,He_Boltz;
  double timeTh,timeH;
  double sq_0,sq_1,a_PPB,b_PPB,c_PPB,d_PPB;

  double cc3P,tau;
  int L,L2,LL;

  /* new in recfast 1.4: */
  double Rdown_trip,Rup_trip,tauHe_s,pHe_s,Doppler,gamma_2Ps,pb,qb,AHcon;
  double tauHe_t,pHe_t,CfHe_t,CL_PSt,gamma_2Pt;
  int Heflag;

  struct thermodynamics_parameters_and_workspace * ptpaw;
  struct precision * ppr;
  struct background * pba;
  struct recombination * preco;
  double * pvecback;

  ptpaw = parameters_and_workspace;
  ppr = ptpaw->ppr;
  pba = ptpaw->pba;
  preco = ptpaw->preco;
  pvecback = ptpaw->pvecback;

  x_H = y[0];
  x_He = y[1];
  x = x_H + preco->fHe * x_He;
  Tmat = y[2];

  n = preco->Nnow * pow((1.+z),3);
  n_He = preco->fHe * preco->Nnow * pow((1.+z),3);
  Trad = preco->Tnow * (1.+z);

  class_call(background_functions(pba,1./(1.+z),short_info,pvecback),
	     pba->error_message,
	     error_message);
  
  Hz=pvecback[pba->index_bg_H]/_Mpc_in_sec_;
  dHdz=-pvecback[pba->index_bg_H_prime]/pvecback[pba->index_bg_H]/pba->a_today;

  Rdown=1.e-19*_a_PPB_*pow((Tmat/1.e4),_b_PPB_)/(1.+_c_PPB_*pow((Tmat/1.e4),_d_PPB_));
  Rup = Rdown * pow((preco->CR*Tmat),1.5)*exp(-preco->CDB/Tmat);

  sq_0 = sqrt(Tmat/_T_0_);
  sq_1 = sqrt(Tmat/_T_1_);
  Rdown_He = _a_VF_/(sq_0 * pow((1.+sq_0),(1.-_b_VF_)) * pow((1. + sq_1),(1. + _b_VF_)));
  Rup_He = Rdown_He*pow((preco->CR*Tmat),1.5)*exp(-preco->CDB_He/Tmat);
  Rup_He = 4.*Rup_He;
  K = preco->CK/Hz;

  /* following is from recfast 1.5 */

  if (ppr->recfast_Hswitch == _TRUE_ )
    K *= 1.
      + ppr->recfast_AGauss1*exp(-pow((log(1.+z)-ppr->recfast_zGauss1)/ppr->recfast_wGauss1,2)) 
      + ppr->recfast_AGauss2*exp(-pow((log(1.+z)-ppr->recfast_zGauss2)/ppr->recfast_wGauss2,2));

  /* end of new recfast 1.5 piece */

  /* following is from recfast 1.4 */

  Rdown_trip = _a_trip_/(sq_0*pow((1.+sq_0),(1.-_b_trip_)) * pow((1.+sq_1),(1.+_b_trip_)));
  Rup_trip = Rdown_trip*exp(-_h_P_*_c_*_L_He2St_ion_/(_k_B_*Tmat))*pow(preco->CR*Tmat,1.5)*4./3.;

  if ((x_He < 5.e-9) || (x_He > 0.980)) 
    Heflag = 0;
  else
    Heflag = ppr->recfast_Heswitch;

  if (Heflag == 0)
    K_He = preco->CK_He/Hz;
  else {
    tauHe_s = _A2P_s_*preco->CK_He*3.*n_He*(1.-x_He)/Hz;
    pHe_s = (1.-exp(-tauHe_s))/tauHe_s;
    K_He = 1./(_A2P_s_*pHe_s*3.*n_He*(1.-x_He));

    /*    if (((Heflag == 2) || (Heflag >= 5)) && (x_H < 0.99999)) { */
    if (((Heflag == 2) || (Heflag >= 5)) && (x_H < 0.9999999)) { /* threshold changed by Antony Lewis in 2008 to get smoother Helium */

      Doppler = 2.*_k_B_*Tmat/(_m_H_*_not4_*_c_*_c_);
      Doppler = _c_*_L_He_2p_*sqrt(Doppler);
      gamma_2Ps = 3.*_A2P_s_*preco->fHe*(1.-x_He)*_c_*_c_
	/(sqrt(_PI_)*_sigma_He_2Ps_*8.*_PI_*Doppler*(1.-x_H))
	/pow(_c_*_L_He_2p_,2.);
      pb = 0.36;
      qb = ppr->recfast_fudge_He;
      AHcon = _A2P_s_/(1.+pb*pow(gamma_2Ps,qb));
      K_He=1./((_A2P_s_*pHe_s+AHcon)*3.*n_He*(1.-x_He));
    }

    if (Heflag >= 3) {
      tauHe_t = _A2P_t_*n_He*(1.-x_He)*3./(8.*_PI_*Hz*pow(_L_He_2Pt_,3.));
      pHe_t = (1. - exp(-tauHe_t))/tauHe_t;
      CL_PSt = _h_P_*_c_*(_L_He_2Pt_ - _L_He_2St_)/_k_B_;
      if ((Heflag == 3) || (Heflag == 5) || (x_H >= 0.99999)) {
	CfHe_t = _A2P_t_*pHe_t*exp(-CL_PSt/Tmat);
	CfHe_t = CfHe_t/(Rup_trip+CfHe_t);
      }
      else {
	Doppler = 2.*_k_B_*Tmat/(_m_H_*_not4_*_c_*_c_);
	Doppler = _c_*_L_He_2Pt_*sqrt(Doppler);
	gamma_2Pt = 3.*_A2P_t_*preco->fHe*(1.-x_He)*_c_*_c_
	  /(sqrt(_PI_)*_sigma_He_2Pt_*8.*_PI_*Doppler*(1.-x_H))
	  /pow(_c_*_L_He_2Pt_,2.);
	pb = 0.66;
	qb = 0.9;
	AHcon = _A2P_t_/(1.+pb*pow(gamma_2Pt,qb))/3.;
	CfHe_t = (_A2P_t_*pHe_t+AHcon)*exp(-CL_PSt/Tmat);
	CfHe_t = CfHe_t/(Rup_trip+CfHe_t);
      }
    }
  }

  /* end of new recfast 1.4 piece */

  timeTh=(1./(preco->CT*pow(Trad,4.)))*(1.+x+preco->fHe)/x;
  timeH=2./(3.*preco->H0*pow(1.+z,1.5)); 

  if (x_H > 0.99)
    dy[0] = 0.;
  else {
    if (x_H > 0.985) {
      dy[0] = (x*x_H*n*Rdown - Rup*(1.-x_H)*exp(-preco->CL/Tmat)) /(Hz*(1.+z));

    /*   fprintf(stderr,"%e %e %e %e %e %e %e %e %e %e %e\n",x,x_H,n,Rdown,Rup,preco->CL,Tmat, */
/* 	      x*x_H*n*Rdown, */
/* 	      Rup*(1.-x_H)*exp(-preco->CL/Tmat), */
/* 	      x*x_H*n*Rdown-Rup*(1.-x_H)*exp(-preco->CL/Tmat), */
/* 	      dy[0]); */

/*       class_test(0==0,error_message,""); */
    }
    else {
      dy[0] = ((x*x_H*n*Rdown - Rup*(1.-x_H)*exp(-preco->CL/Tmat)) *(1. + K*_Lambda_*n*(1.-x_H))) /(Hz*(1.+z)*(1./preco->fu+K*_Lambda_*n*(1.-x)/preco->fu +K*Rup*n*(1.-x)));
    }
  }

  if (x_He < 1.e-15) 
    dy[1]=0.;
  else {

    if (preco->Bfact/Tmat < 680.) 
      He_Boltz=exp(preco->Bfact/Tmat);
    else 
      He_Boltz=exp(680.);

    dy[1] = ((x*x_He*n*Rdown_He - Rup_He*(1.-x_He)*exp(-preco->CL_He/Tmat)) *(1. + K_He*_Lambda_He_*n_He*(1.-x_He)*He_Boltz)) /(Hz*(1+z) * (1. + K_He*(_Lambda_He_+Rup_He)*n_He*(1.-x_He)*He_Boltz));

    /* following is from recfast 1.4 */

    if (Heflag >= 3)
      dy[1] = dy[1] + 
	(x*x_He*n*Rdown_trip
	 - (1.-x_He)*3.*Rup_trip*exp(-_h_P_*_c_*_L_He_2St_/(_k_B_*Tmat)))
	*CfHe_t/(Hz*(1.+z));

    /* end of new recfast 1.4 piece */
  }

  if (timeTh < preco->H_frac*timeH) {
    dy[2]=Tmat/(1.+z);
    /* v 1.5: like in camb, add here a smoothing term as suggested by Adam Moss */
/*     epsilon = Hz * (1.+x+preco->fHe) / (preco->CT*pow(Trad,3)*x); */
/*     dy[2] = preco->Tnow + epsilon*((1.+preco->fHe)/(1.+preco->fHe+x))*((y[0]+preco->fHe*y[1])/x)  */
/*       - epsilon* dHdz/Hz + 3.*epsilon/(1.+z); */
    
  }
  else
    dy[2]= preco->CT * pow(Trad,4) * x / (1.+x+preco->fHe) * (Tmat-Trad) / (Hz*(1.+z)) + 2.*Tmat/(1.+z);

  return _SUCCESS_;
}      

/* int get_init( */
/* 	     double z,	     /\**< Input  : redshift *\/ */
/* 	     double * y      /\**< Ouput : x_H0, x_He0, T *\/ */
/* 	     ) { */
  
/*   double rhs,x_H0,x_He0,x0; */
  
/*   if(z > 8000.) { */
/*     x_H0 = 1.; */
/*     x_He0 = 1.; */
/*     x0 = 1.+2.*preco->fHe; */
/*   } */
/*   else { */
/*     if (z > 3500.) { */
/*       x_H0 = 1.; */
/*       x_He0 = 1.; */
/*       rhs = exp( 1.5 * log(preco->CR*preco->Tnow/(1.+z)) - preco->CB1_He2/(preco->Tnow*(1.+z)) ) / preco->Nnow; */
/*       rhs = rhs*1.;  /\*ratio of g's is 1 for He++ <-> He+*\/ */
/*       x0 = 0.5 * ( sqrt( pow((rhs-1.-preco->fHe),2) + 4.*(1.+2.*preco->fHe)*rhs) - (rhs-1.-preco->fHe) ); */
/*     } */
/*     else { */
/*       if(z > 2000.) { */
/*         x_H0 = 1.; */
/* 	rhs = exp( 1.5 * log(preco->CR*preco->Tnow/(1.+z)) - preco->CB1_He1/(preco->Tnow*(1.+z)) ) / preco->Nnow; */
/*         rhs = rhs*4.;  /\*ratio of g's is 4 for He+ <-> He0*\/ */
/* 	x_He0 = 0.5  * ( sqrt( pow((rhs-1.),2) + 4.*(1.+preco->fHe)*rhs )- (rhs-1.)); */
/* 	x0 = x_He0; */
/* 	x_He0 = (x0 - 1.)/preco->fHe; */
/*       } */
/*       else { */
/* 	rhs = exp( 1.5 * log(preco->CR*preco->Tnow/(1.+z)) - preco->CB1/(preco->Tnow*(1.+z)) ) / preco->Nnow; */
/* 	x_H0 = 0.5 * (sqrt( pow(rhs,2)+4.*rhs ) - rhs ); */
/* 	x_He0 = 0.; */
/* 	x0 = x_H0; */
/*       } */
/*     } */
/*   } */

/*   y[0]=x_H0; */
/*   y[1]=x_He0; */
/*   y[2]=preco->Tnow*(1.+z); */

/*   return _SUCCESS_; */
/* } */
