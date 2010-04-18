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

/** @name - structures used within the thermodynamics module: */

//@{

struct background * pba; /**< a cosmo structure pointer for internal use in the thermodynamics module */
struct precision * ppr; /**< a precision_params structure pointer for internal use in the thermodynamics module */
struct thermo * pth; /**< a thermo structure pointer for internal use in the thermodynamics module */

double * pvecback_th; /**< vector of background quantities, used
			 throughout the thermodynamics module.Use a
			 global variable in order to avoid
			 reallocating it many times. */

/** @name - recfast variables: */

//@{
double CDB,CR,CK,CL,CT,fHe,CDB_He,CK_He,CL_He,fu,H_frac,Tnow,Nnow,Bfact;
double C1P3P,cc3P1P,cccP,hck,CB1,CB1_He1,CB1_He2,H0; 
double kb_over_mH; /* k_B/mH */

//@}

/** @name - miscellaneous: */

//@{

char * errmsg; /**< error management pointer */
char Transmit_Error_Message[2048]; /**< contains error message */

//@}

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
			double z,
			enum interpolation_mode intermode,
			int * last_index,
			double * pvecthermo_local
			) {

  /** Summary: */

  /** - define local variables */
  int i;
  double x0;

  /** - chech that z is in the pre-computed range; if \f$ z>z_{max}
      \f$, extrapolate */

  /* deal with case z < z_min=0, unphysical */
  if (z < pth->z_table[0]) {
    sprintf(pth->error_message,"%s(L:%d) : z < z_min=%d",__func__,__LINE__,pth->z_table[0]);
    return _FAILURE_;
  }

  /* deal with case z > z_max, not stored inside interpolation table,
     but doable with simple analytic approximations */
  if (z >= pth->z_table[pth->tt_size-1]) {

    x0= pth->thermodynamics_table[(pth->tt_size-1)*pth->th_size+pth->index_th_xe];
    
    /* ionization fraction */
    pvecthermo_local[pth->index_th_xe] = x0;

    /* Calculate Tb */
    pvecthermo_local[pth->index_th_Tb] = Tnow*(1.+z);

    /* Calculate cb2 (cb2 = (k_B/mu) Tb (1-1/3 dlnTb/dlna) = (k_B/mu) Tb (1+1/3 (1+z) dlnTb/dz)) */
    /* note that m_H / mu = 1 + (m_H/m_He-1) Y_p + x_e (1-Y_p) */
    pvecthermo_local[pth->index_th_cb2] = kb_over_mH * (1. + (1./_not4_ - 1.) * pth->YHe + x0 * (1.-pth->YHe)) * Tnow * (1.+z) * 4. / 3.;

    /* Calculate dkappa/deta (dkappa/deta = a n_e x_e sigma_T = a^{-2} n_e(today) x_e sigma_T in units of 1/Mpc) */
    pvecthermo_local[pth->index_th_dkappa] = (1.+z) * (1.+z) * Nnow * x0 * _sigma_ * _Mpc_over_m_;


    /* Calculate dz/deta = -H with background_functions_of_a() */
    if (background_functions_of_a(1./(1.+z),short_info,pvecback_th) == _FAILURE_) {
      sprintf(pth->error_message,"%s(L:%d) : Error calling background_functions_of_a \n=>%s",__func__,__LINE__,pba->error_message);
      return _FAILURE_;
    }

    /* Calculate d2kappa/deta2 = dz/deta d/dz[dkappa/deta] */
    pvecthermo_local[pth->index_th_ddkappa] = -pvecback_th[pba->index_bg_H] * 2. / (1.+z) * pvecthermo_local[pth->index_th_dkappa];

    /* Calculate d3kappa/deta3 = dz/deta d/dz[d2kappa/deta2] */
    /*     pvecthermo_local[pth->index_th_ddkappa] = (pvecback_th[pba->index_bg_Hdot] / pvecback_th[pba->index_bg_H] - */
    /*        pvecback_th[pba->index_bg_H]	/ (1.+z)) * pvecthermo_local[pth->index_th_ddkappa]; */

    /* visibility = g = (dkappa/deta) * exp(- kappa) */
    pvecthermo_local[pth->index_th_g]=0.;

  }

  /** - otherwise, just interpolate in table with array_interpolate() */
  else {

    if (intermode == normal) {

      if (array_interpolate_spline(
				   pth->z_table,
				   pth->tt_size,
				   pth->thermodynamics_table,
				   pth->d2thermodynamics_dz2_table,
				   pth->th_size,
				   z,
				   last_index,
				   pvecthermo_local,
				   pth->th_size,
				   errmsg)== _FAILURE_) {
	sprintf(pth->error_message,"%s(L:%d) : error in array_interpolate_spline() \n=>%s",__func__,__LINE__,errmsg);
	return _FAILURE_;
      }

    }

    if (intermode == closeby) {

      if (array_interpolate_spline_growing_closeby(
						   pth->z_table,
						   pth->tt_size,
						   pth->thermodynamics_table,
						   pth->d2thermodynamics_dz2_table,
						   pth->th_size,
						   z,
						   last_index,
						   pvecthermo_local,
						   pth->th_size,
						   errmsg)== _FAILURE_) {
	sprintf(pth->error_message,"%s(L:%d) : error in array_interpolate_spline_growing_closeby() \n=>%s",__func__,__LINE__,errmsg);
	return _FAILURE_;
      }

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
			struct background * pba_input,
			struct precision * ppr_input,
			struct thermo * pth_output
			) {

  /** Summary: */

  /** - define local variables */

  /* running index over vector of thermodynamics variables */
  int i, index_th, index_re;
  /* visibility function */
  double g, g_previous, eta_visibility_max;
  /* for calling background_at_eta() */
  int last_index_back;
  double * pvecback_th_long;
  double * eta_table; /**< list of eta values associated with z values in pth->z_table */

  struct recombination reco;
  struct reionization reio;
  struct recombination * preco;
  struct reionization * preio;

  /** - identify the cosmo, precision and thermo structures pba, ppr, pth (used throughout thermodynamics.c as global variables) to the input/output structures of this function (pba and ppr are already filled, pth will be filled by this function) */
  pba = pba_input;
  ppr = ppr_input;
  pth = pth_output;

  preco=&reco;
  preio=&reio;

  if (pth->thermodynamics_verbose > 0)
    printf("Computing thermodynamics\n");

  /* Tcmb in K */
  if ((pth->Tcmb < _TCMB_SMALL_)||(pth->Tcmb > _TCMB_BIG_)) {
    sprintf(pth->error_message,"%s(L:%d): Tcmb out of bounds (got %g, expected %g<Tcmb<%g)",__func__,__LINE__,pth->Tcmb,_TCMB_SMALL_,_TCMB_BIG_);
    return _FAILURE_;
  }

  /* Y_He */
  if ((pth->YHe < _YHE_SMALL_)||(pth->YHe > _YHE_BIG_)) {
    sprintf(pth->error_message,"%s(L:%d): Y_He out of bounds (got %g, expected %g<Y_He<%g)",__func__,__LINE__,pth->YHe,_YHE_SMALL_,_YHE_BIG_);
    return _FAILURE_;
  }

  /* tests in order to prevent segmentaltion fault in the following */
  if (_not4_ == 0.) {
    sprintf(pth->error_message,"%s(L:%d) : you have _not4_=0, stop to avoid division by zero",__func__,__LINE__);
    return _FAILURE_;
  }
  if (pth->YHe == 1.) {
    sprintf(pth->error_message,"%s(L:%d) : you have YHe=1, stop to avoid division by zero",__func__,__LINE__);
    return _FAILURE_;
  }

  /** - assign values to all indices in the structures with thermodynamics_indices()*/
  if (thermodynamics_indices(pth,preco,preio) == _FAILURE_) {
    sprintf(Transmit_Error_Message,"%s",pth->error_message);
    sprintf(pth->error_message,"%s(L:%d) : error in thermodynamics_indices()\n=>%s",__func__,__LINE__,Transmit_Error_Message);
    return _FAILURE_;
  }

  /** - allocate memory for the background vector pvecback_th (global variable in thermodynamics.c, used as long as thermodynamics_free() is not called) */
  pvecback_th = malloc(pba->bg_size_short*sizeof(double));
  if (pvecback_th==NULL) {
    sprintf(pth->error_message,"%s(L:%d): Cannot allocate pvecback_th \n",__func__,__LINE__);
    return _FAILURE_;
  }

  /** - solve recombination and store values of \f$ z, x_e, d \kappa / d \eta, T_b, c_b^2 \f $ with thermodynamics_recombination() */
  if (thermodynamics_recombination(preco) == _FAILURE_) {
    sprintf(Transmit_Error_Message,"%s",pth->error_message);
    sprintf(pth->error_message,"%s(L:%d) : error in thermodynamics_recombination()\n=>%s",__func__,__LINE__,Transmit_Error_Message);
    return _FAILURE_;
  }

  /** - solve reionization and store values of \f$ z, x_e, d \kappa / d z \f $ with thermodynamics_recombination()*/
  if (pth->reio_parametrization != reio_none) {
    if (thermodynamics_reionization(preco,preio) == _FAILURE_) {
      sprintf(Transmit_Error_Message,"%s",pth->error_message);
      sprintf(pth->error_message,"%s(L:%d) : error in thermodynamics_reionization()\n=>%s",__func__,__LINE__,Transmit_Error_Message);
      return _FAILURE_;
    }
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
    if (preco->recombination_table
	[preio->index_reco_when_reio_start*preco->re_size+preco->index_re_z] !=
	preio->reionization_table[(preio->rt_size -1)*preio->re_size+preio->index_re_z]) {
      sprintf(pth->error_message,"%s(L:%d) : mismatch which should never happen \n",__func__,__LINE__);
      return _FAILURE_;
    }
  }

  /* number of redshift in full table = number in reco + number in reio - overlap */
  pth->tt_size = ppr->recfast_Nz0 + preio->rt_size - preio->index_reco_when_reio_start - 1;
  pth->z_table = malloc(pth->tt_size*sizeof(double));
  if (pth->z_table==NULL) {
    sprintf(pth->error_message,"%s(L:%d): Cannot allocate pth->z_table \n",__func__,__LINE__);
    return _FAILURE_;
  }
  eta_table = malloc(pth->tt_size*sizeof(double));
  if (eta_table==NULL) {
    sprintf(pth->error_message,"%s(L:%d): Cannot allocate eta_table \n",__func__,__LINE__);
    return _FAILURE_;
  }
  pth->thermodynamics_table = malloc(pth->th_size*pth->tt_size*sizeof(double));
  if (pth->thermodynamics_table==NULL) {
    sprintf(pth->error_message,"%s(L:%d): Cannot allocate pth->thermodynamics_table \n",__func__,__LINE__);
    return _FAILURE_;
  }
  pth->d2thermodynamics_dz2_table = malloc(pth->th_size*pth->tt_size*sizeof(double));  
  if (pth->d2thermodynamics_dz2_table==NULL) {
    sprintf(pth->error_message,"%s(L:%d): Cannot allocate pth->d2thermodynamics_dz2_table \n",__func__,__LINE__);
    return _FAILURE_;
  }

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
    if (background_eta_of_z(pth->z_table[i],eta_table+i) == _FAILURE_){
      sprintf(pth->error_message,"%s(L:%d) : error in background_eta_of_z()\n=>%s",__func__,__LINE__,pba->error_message);
      return _FAILURE_;
    }
  }

  /** - fill missing columns (quantities not computed previously but
      related); eventually, reduce numerical errors with
      thermodynamics_cure_discontinuity()
  */

  /* -> second derivative with respect to eta of dkappa */
  /** - fill tables of second derivatives (in view of spline interpolation) */
  if (array_spline_table_line_to_line(eta_table,
				      pth->tt_size,
				      pth->thermodynamics_table,
				      pth->th_size,
				      pth->index_th_dkappa,
				      pth->index_th_dddkappa,
				      _SPLINE_EST_DERIV_,
				      errmsg) == _FAILURE_) {
    sprintf(pth->error_message,"%s(L:%d) : error in array_spline_table_line_to_line()\n=>%s",__func__,__LINE__,errmsg);
    return _FAILURE_;
  }

  /* -> first derivative with respect to eta of dkappa */
  if (array_derive_spline_table_line_to_line(eta_table,
					     pth->tt_size,
					     pth->thermodynamics_table,
					     pth->th_size,
					     pth->index_th_dkappa,
					     pth->index_th_dddkappa,
					     pth->index_th_ddkappa,
					     errmsg) == _FAILURE_) {
    sprintf(pth->error_message,"%s(L:%d) : error in array_derive_spline_table_line_to_line()\n=>%s",__func__,__LINE__,errmsg);
    return _FAILURE_;
  }

  /* -> compute -kappa = [int_{eta_today}^{eta} deta dkappa/deta], store temporarily in column "g" */ 
  if (array_integrate_spline_table_line_to_line(eta_table,
						pth->tt_size,
						pth->thermodynamics_table,
						pth->th_size,
						pth->index_th_dkappa,
						pth->index_th_dddkappa,
						pth->index_th_g,
						errmsg) == _FAILURE_) {
    sprintf(pth->error_message,"%s(L:%d) : error in array_integrate_spline_table_line_to_line()\n=>%s",__func__,__LINE__,errmsg);
    return _FAILURE_;
  }

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

    /* maximum variation rate */
    if (pth->thermodynamics_table[i*pth->th_size+pth->index_th_dkappa] != 0.) {
      pth->thermodynamics_table[i*pth->th_size+pth->index_th_rate] =
	sqrt(pow(pth->thermodynamics_table[i*pth->th_size+pth->index_th_dkappa],2.)
	     +pow(pth->thermodynamics_table[i*pth->th_size+pth->index_th_ddkappa]/ 
		  pth->thermodynamics_table[i*pth->th_size+pth->index_th_dkappa],2.)
	     +fabs(pth->thermodynamics_table[i*pth->th_size+pth->index_th_dddkappa]/ 
		   pth->thermodynamics_table[i*pth->th_size+pth->index_th_dkappa]));
    }
    else {
      sprintf(pth->error_message,"%s(L:%d) : [dkappa/deta vanishes], variation rate diverges",__func__,__LINE__);
      return _FAILURE_;
    }

  }

  /* -> compute smooth the rate (detaile dof smoothing unimportant: only the order of magnitude of the rate matters) */ 
  if (array_smooth(pth->thermodynamics_table,
		   pth->th_size,
		   pth->tt_size,
		   pth->index_th_rate,
		   ppr->thermo_rate_smoothing_radius,
		   errmsg) == _FAILURE_) {
    sprintf(pth->error_message,"%s(L:%d) : error in array_smooth()\n=>%s",__func__,__LINE__,errmsg);
    return _FAILURE_;
  }

  /* check consistency of these values */

  if (pth->z_visibility_start_sources == 0.) {
    sprintf(pth->error_message,"%s(L:%d) : Source functions cannot be sampled, probably because ppr->visibility_threshold_start_sources=%e is too large and exceeds maximum of g",__func__,__LINE__);
    return _FAILURE_;
  }

  if (pth->z_visibility_max == 0.) {
    sprintf(pth->error_message,"%s(L:%d) : Maximum of visibility function, hence recombination time could not be found",__func__,__LINE__);
    return _FAILURE_;
  }

  if (pth->z_visibility_free_streaming > pth->z_visibility_max) {
    sprintf(pth->error_message,"%s(L:%d) : pth->z_visibility_free_streaming=%e should never be larger than pth->z_visibility_max=%e \n=>%s",__func__,__LINE__,pth->z_visibility_free_streaming,pth->z_visibility_max);
    return _FAILURE_;
  }

  /** - find conformal recombination time using background_eta_of_z() **/
  if (background_eta_of_z(pth->z_visibility_max,&eta_visibility_max)
      == _FAILURE_) {
    sprintf(pth->error_message,"%s(L:%d) : error in background_eta_of_z() \n=>%s",__func__,__LINE__,pba->error_message);
    return _FAILURE_;
  }

  pth->eta_rec = eta_visibility_max;
  /*   printf("eta_rec=%e\n",eta_visibility_max); */
  /*   printf("eta_0=%e\n",pba->conformal_age-eta_visibility_max); */
  /*   printf("eta_0-eta_rec=%e\n",pba->conformal_age-eta_visibility_max); */

  /** - find sound horizon at recombination using background_at_eta() */
  pvecback_th_long = malloc(pba->bg_size*sizeof(double));
  if (pvecback_th_long==NULL) {
    sprintf(pth->error_message,"%s(L:%d): Cannot allocate pvecback_th_long \n",__func__,__LINE__);
    return _FAILURE_;
  }
  if (background_at_eta(pth->eta_rec, long_info, normal, &last_index_back, pvecback_th_long) == _FAILURE_) {
    sprintf(pth->error_message,"%s(L:%d) : error in background_at_eta()\n=>%s",__func__,__LINE__,pba->error_message);
    return _FAILURE_;
  }  
  pth->rs_rec=pvecback_th_long[pba->index_bg_rs];
  free(pvecback_th_long);


  /** - fill tables of second derivatives with respect to z (in view of spline interpolation) */
  if (array_spline_table_lines(pth->z_table,
			       pth->tt_size,
			       pth->thermodynamics_table,
			       pth->th_size,
			       pth->d2thermodynamics_dz2_table,
			       _SPLINE_EST_DERIV_,
			       errmsg) == _FAILURE_) {
    sprintf(pth->error_message,"%s(L:%d) : error in array_spline_table_lines()\n=>%s",__func__,__LINE__,errmsg);
    return _FAILURE_;
  }

  if (pth->thermodynamics_verbose > 0) {
    printf(" -> recombination at z = %f\n",pth->z_visibility_max);
    if (pth->reio_parametrization != reio_none) {
      if (pth->reio_z_or_tau==reio_tau)
	printf(" -> reionization  at z = %f\n",pth->z_reio);
      if (pth->reio_z_or_tau==reio_z)
	printf(" -> reionization with optical depth = %f\n",pth->tau_reio);
    }
  }

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
int thermodynamics_free() {

  free(pth->z_table);
  free(pth->thermodynamics_table);
  free(pth->d2thermodynamics_dz2_table);
  free(pvecback_th);

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
int thermodynamics_indices(struct thermo * pthermo,
			   struct recombination * preco,
			   struct reionization * preio) {

  /** Summary: */

  /** - define local variables */

  /* a running index for the vector of thermodynamics quantities */
  int index;

  /** - intialization of all indices and flags in thermo structure */
  index = 0;
  
  pthermo->index_th_xe = index;
  index++;
  pthermo->index_th_dkappa = index;
  index++;
  pthermo->index_th_ddkappa = index;
  index++;
  pthermo->index_th_dddkappa = index;
  index++;
  pthermo->index_th_exp_m_kappa = index;
  index++;
  pthermo->index_th_g = index;
  index++;
  pthermo->index_th_dg = index;
  index++;
  pthermo->index_th_Tb = index;
  index++;
  pthermo->index_th_cb2 = index;
  index++;
  pthermo->index_th_rate = index;
  index++;

  /* end of indices */
  pthermo->th_size = index;

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
					 struct reionization * preio,
					 double * xe
					 ) {

  double argument;

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

  sprintf(pth->error_message,"%s(L:%d) : value of reio_parametrization=%d unclear",__func__,__LINE__,pth->reio_parametrization);
  return _FAILURE_;

}


int thermodynamics_get_xe_before_reionization(double z,
					      struct recombination * preco,
					      double * xe) {
  int i;

  i=0;
  while (preco->recombination_table[i*preco->re_size+preco->index_re_z] < z) {
    i++;
    if (i == ppr->recfast_Nz0) {
      sprintf(pth->error_message,"%s(L:%d) : z = %e > largest redshift in thermodynamics table \n",__func__,__LINE__,ppr->reionization_z_start_max);
      return _FAILURE_;
    }
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
int thermodynamics_reionization(struct recombination * preco,
				struct reionization * preio) {

  int i,counter;
  double z_sup,z_mid,z_inf;
  double tau_sup,tau_mid,tau_inf;

  /** (a) if reionization implemented like in CAMB */

  if (pth->reio_parametrization == reio_camb) {

    /** - allocate the parameters of the function \f$ X_e(z) \f$ */

    preio->reionization_parameters=malloc(preio->reio_num_params*sizeof(double)); 

    if (preio->reionization_parameters==NULL) {
      sprintf(pth->error_message,"%s(L:%d): Cannot allocate reionization_parameters \n",__func__,__LINE__);
      return _FAILURE_;
    }
    
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

    if (preio->reionization_parameters[preio->index_reio_exponent]==0) {
      sprintf(pth->error_message,"%s(L:%d) : reio_exponent is null, stop to avoid division by zero",__func__,__LINE__);
      return _FAILURE_;
    }
    if (preio->reionization_parameters[preio->index_reio_width]==0) {
      sprintf(pth->error_message,"%s(L:%d) : reio_width is null, stop to avoid division by zero",__func__,__LINE__);
      return _FAILURE_;
    }
    if (preio->reionization_parameters[preio->index_helium_fullreio_width]==0) {
      sprintf(pth->error_message,"%s(L:%d) : fullreio_width is null, stop to avoid division by zero",__func__,__LINE__);
      return _FAILURE_;
    }

    /** - if reionization redshift given as an input, initialize the remaining values and fill reionization table*/

    if (pth->reio_z_or_tau == reio_z) {
      
      /* reionization redshift */
      preio->reionization_parameters[preio->index_reio_redshift] = pth->z_reio; 
      /* infer starting redshift */
      preio->reionization_parameters[preio->index_reio_start] = preio->reionization_parameters[preio->index_reio_redshift]+ppr->reionization_start_factor*ppr->reionization_width;
      if (preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max) {
	sprintf(pth->error_message,"%s(L:%d) : starting redshift for reionization > reionization_z_start_max = %e\n",__func__,__LINE__,ppr->reionization_z_start_max);
	return _FAILURE_;
      }
      /* infer xe_before_reio */
      if(thermodynamics_get_xe_before_reionization(preio->reionization_parameters[preio->index_reio_redshift],
						   preco,
						   &(preio->reionization_parameters[preio->index_reio_xe_before])) == _FAILURE_) {
	sprintf(Transmit_Error_Message,"%s",pth->error_message);
	sprintf(pth->error_message,"%s(L:%d) : error in thermodynamics_get_xe_before_reionization()\n=>%s",__func__,__LINE__,Transmit_Error_Message);
	return _FAILURE_;
      }
      /* fill reionization table */
      if (thermodynamics_reionization_discretize(preco,preio) == _FAILURE_) {
	sprintf(Transmit_Error_Message,"%s",pth->error_message);
	sprintf(pth->error_message,"%s(L:%d) : error in thermodynamics_reionization_discretize()\n=>%s",__func__,__LINE__,Transmit_Error_Message);
	return _FAILURE_;
      }

      pth->tau_reio=preio->reionization_optical_depth;

    }

    /** - if reionization optical depth given as an input, find reionization redshift by dichotomy and initialize the remaining values */

    if (pth->reio_z_or_tau == reio_tau) {

      /* upper value */

      z_sup = ppr->reionization_z_start_max-ppr->reionization_start_factor*ppr->reionization_width;
      if (z_sup < 0.) {
	sprintf(pth->error_message,"%s(L:%d) : Parameters are such that reionization cannot take place before today while starting after z_start_max; need to increase z_start_max",__func__,__LINE__);
	return _FAILURE_;
      }

      /* reionization redshift */
      preio->reionization_parameters[preio->index_reio_redshift] = z_sup; 
      /* infer starting redshift */
      preio->reionization_parameters[preio->index_reio_start] = ppr->reionization_z_start_max;
      /* infer xe_before_reio */
      if(thermodynamics_get_xe_before_reionization(preio->reionization_parameters[preio->index_reio_redshift],
						   preco,
						   &(preio->reionization_parameters[preio->index_reio_xe_before])) == _FAILURE_) {
	sprintf(Transmit_Error_Message,"%s",pth->error_message);
	sprintf(pth->error_message,"%s(L:%d) : error in thermodynamics_get_xe_before_reionization()\n=>%s",__func__,__LINE__,Transmit_Error_Message);
	return _FAILURE_;
      }
      /* fill reionization table */
      if (thermodynamics_reionization_discretize(preco,preio) == _FAILURE_) {
	sprintf(Transmit_Error_Message,"%s",pth->error_message);
	sprintf(pth->error_message,"%s(L:%d) : error in thermodynamics_reionization_discretize()\n=>%s",__func__,__LINE__,Transmit_Error_Message);
	return _FAILURE_;
      }
      tau_sup=preio->reionization_optical_depth;
      if (tau_sup < pth->tau_reio) {
	sprintf(pth->error_message,"%s(L:%d) : parameters are such that reionization cannot start after z_start_max",__func__,__LINE__,pth->error_message);
	return _FAILURE_;
      }

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
	if (preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max) {
	  sprintf(pth->error_message,"%s(L:%d) : starting redshift for reionization > reionization_z_start_max = %e\n",__func__,__LINE__,ppr->reionization_z_start_max);
	  return _FAILURE_;
	}
	/* infer xe_before_reio */
	if(thermodynamics_get_xe_before_reionization(preio->reionization_parameters[preio->index_reio_redshift],
						     preco,
						     &(preio->reionization_parameters[preio->index_reio_xe_before])) == _FAILURE_) {
	  sprintf(pth->error_message,"%s(L:%d) : reionization_z_start_max = %e > largest z in thermodynamics table \n",__func__,__LINE__,ppr->reionization_z_start_max);
	  return _FAILURE_;
	}
	/* clean and fill reionization table */
	free(preio->reionization_table);
	if (thermodynamics_reionization_discretize(preco,preio) == _FAILURE_) {
	  sprintf(Transmit_Error_Message,"%s",pth->error_message);
	  sprintf(pth->error_message,"%s(L:%d) : error in thermodynamics_reionization_discretize()\n=>%s",__func__,__LINE__,Transmit_Error_Message);
	  return _FAILURE_;
	}
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
	if (counter > _MAX_IT_) {
	  sprintf(pth->error_message,"%s(L:%d) : while searching for reionization_optical_depth, maximum number of iterations exceeded",__func__,__LINE__);
	  return _FAILURE_;
	}
	
      }

      /* store z reionization in thermodynamics structure */
      pth->z_reio=preio->reionization_parameters[preio->index_reio_redshift];

    }

    free(preio->reionization_parameters);

    return _SUCCESS_;

  }

  sprintf(pth->error_message,"%s(L:%d) : value of reio_z_or_tau=%d unclear",__func__,__LINE__,pth->reio_z_or_tau);
  return _FAILURE_;

}

int thermodynamics_reionization_discretize(
					   struct recombination * preco,
					   struct reionization * preio
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
  reio_vector=malloc((preio->re_size)*sizeof(double));
  if (reio_vector==NULL) {
    sprintf(pth->error_message,"%s(L:%d): Cannot allocate reio_vector \n",__func__,__LINE__);
    return _FAILURE_;
  }


  /** (b) create a growTable with gt_init() */
  if (gt_init(&gTable) == _FAILURE_) {
    sprintf(pth->error_message,"%s(L:%d) : error in gt_init() \n=>%s",__func__,__LINE__,gTable.error_message);
  }

  /** (c) first line is taken from thermodynamics table, just before reionization starts */

  /** - look where to start in current thermodynamics table */
  i=0;
  while (preco->recombination_table[i*preco->re_size+preco->index_re_z] < preio->reionization_parameters[preio->index_reio_start]) {
    i++;
    if (i == ppr->recfast_Nz0) {
      sprintf(pth->error_message,"%s(L:%d) : reionization_z_start_max = %e > largest redshift in thermodynamics table",__func__,__LINE__,ppr->reionization_z_start_max);
      return _FAILURE_;
    }
  }

  /** - get redshift */
  z=preco->recombination_table[i*preco->re_size+preco->index_re_z];
  reio_vector[preio->index_re_z]=z;
  preio->index_reco_when_reio_start=i;

  /** - get \f$ X_e \f$ */
  if (thermodynamics_reionization_function(z,preio,&xe) == _FAILURE_) {
    sprintf(Transmit_Error_Message,"%s",pth->error_message);
    sprintf(pth->error_message,"%s(L:%d) : error in thermodynamics_reionization_function()\n=>%s",__func__,__LINE__,Transmit_Error_Message);
    return _FAILURE_;
  }
  reio_vector[preio->index_re_xe] = xe;

  /** - get \f$ d kappa / d z = (d kappa / d eta) * (d eta / d z) = - (d kappa / d eta) / H \f$ */
  if (background_functions_of_a(1./(1.+z),short_info,pvecback_th) == _FAILURE_){
    sprintf(pth->error_message,"%s(L:%d) : error in background_functions_of_a()\n=>%s",__func__,__LINE__,pba->error_message);
    return _FAILURE_;
  }
  reio_vector[preio->index_re_dkappadeta] = (1.+z) * (1.+z) * Nnow * xe * _sigma_ * _Mpc_over_m_;
  if (pvecback_th[pba->index_bg_H] != 0.) {
    reio_vector[preio->index_re_dkappadz] = reio_vector[preio->index_re_dkappadeta] / pvecback_th[pba->index_bg_H];
  }
  else {
    sprintf(pth->error_message,"%s(L:%d) : H is null, stop to avoid division by zero",__func__,__LINE__);
    return _FAILURE_;
  }
  dkappadz = reio_vector[preio->index_re_dkappadz];
  dkappadeta = reio_vector[preio->index_re_dkappadeta];

  /** - get baryon temperature **/
  Tb = preco->recombination_table[i*preco->re_size+preco->index_re_Tb];
  reio_vector[preio->index_re_Tb] = Tb;

  /** - after recombination, Tb scales like (1+z)**2. Compute constant factor Tb/(1+z)**2. */
  Tba2 = Tb/(1+z)/(1+z);

  /** - get baryon sound speed */
  reio_vector[preio->index_re_cb2] = 5./3. * kb_over_mH * (1. + (1./_not4_ - 1.) * Yp + xe * (1.-Yp)) * Tb;

  /** - store these values in growing table */
  if (gt_add(&gTable,_GT_END_,(void *) reio_vector,sizeof(double)*(preio->re_size)) == _FAILURE_) { 
    sprintf(pth->error_message,"%s(L:%d) : error in gt_add()\n=>%s",__func__,__LINE__,gTable.error_message);
    return _FAILURE_;
  }
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

    if (thermodynamics_reionization_function(z_next,preio,&xe_next) == _FAILURE_) {
      sprintf(Transmit_Error_Message,"%s",pth->error_message);
      sprintf(pth->error_message,"%s(L:%d) : error in thermodynamics_reionization_function()\n=>%s",__func__,__LINE__,Transmit_Error_Message);
      return _FAILURE_;
    }  
    if (background_functions_of_a(1./(1.+z_next),short_info,pvecback_th) == _FAILURE_){
      sprintf(pth->error_message,"%s(L:%d) : error in background_functions_of_a()\n=>%s",__func__,__LINE__,pba->error_message);
      return _FAILURE_;
    }
    if (pvecback_th[pba->index_bg_H] != 0.) {
      dkappadz_next= (1.+z_next) * (1.+z_next) * Nnow * xe_next * _sigma_ * _Mpc_over_m_ / pvecback_th[pba->index_bg_H];
    }
    else {
      sprintf(pth->error_message,"%s(L:%d) : H is null, stop to avoid division by zero",__func__,__LINE__);
      return _FAILURE_;
    }
    dkappadeta_next= (1.+z_next) * (1.+z_next) * Nnow * xe_next * _sigma_ * _Mpc_over_m_;

    if ((dkappadz == 0.) && (dkappadeta == 0.)) {
      sprintf(pth->error_message,"%s(L:%d) : dkappadz=%e, dkappadeta=%e, stop to avoid division by zero",__func__,__LINE__,dkappadz,dkappadeta);
      return _FAILURE_;
    }

    /** - reduce step if necessary */
    while (((fabs(dkappadz_next-dkappadz)/dkappadz) > ppr->reionization_sampling) || 
	   ((fabs(dkappadeta_next-dkappadeta)/dkappadeta) > ppr->reionization_sampling)) {
      dz*=0.9;
      if (dz < ppr->smallest_allowed_variation) {
	sprintf(pth->error_message,"%s(L:%d) : error : integration step =%e < machine precision : leads either to numerical error or infinite loop\n",__func__,__LINE__,dz);
	return _FAILURE_;
      }
      z_next=z-dz;
      if (thermodynamics_reionization_function(z_next,preio,&xe_next) == _FAILURE_) {
	sprintf(Transmit_Error_Message,"%s",pth->error_message);
	sprintf(pth->error_message,"%s(L:%d) : error in thermodynamics_reionization_function()\n=>%s",__func__,__LINE__,Transmit_Error_Message);
	return _FAILURE_;
      }
      if (background_functions_of_a(1./(1.+z_next),short_info,pvecback_th) == _FAILURE_){
	sprintf(pth->error_message,"%s(L:%d) : error in background_functions_of_a()\n=>%s",__func__,__LINE__,pba->error_message);
	return _FAILURE_;
      }
      if (pvecback_th[pba->index_bg_H] != 0.) {
	dkappadz_next= (1.+z_next) * (1.+z_next) * Nnow * xe_next * _sigma_ * _Mpc_over_m_ / pvecback_th[pba->index_bg_H];
      }
      else {
	sprintf(pth->error_message,"%s(L:%d) : H is null, stop to avoid division by zero",__func__,__LINE__);
	return _FAILURE_;
      }
      dkappadeta_next= (1.+z_next) * (1.+z_next) * Nnow * xe_next * _sigma_ * _Mpc_over_m_;
    }

    /** - get \f$ z, X_e, d kappa / d z \f$ and store in growing table */
    z=z_next;
    xe=xe_next;
    dkappadz=dkappadz_next;
    dkappadeta= dkappadeta_next;

    if ((dkappadz == 0.) && (dkappadeta == 0.)) {
      sprintf(pth->error_message,"%s(L:%d) : dkappadz=%e, dkappadeta=%e, stop to avoid division by zero",__func__,__LINE__,dkappadz,dkappadeta);
      return _FAILURE_;
    }

    reio_vector[preio->index_re_z] = z;   
    reio_vector[preio->index_re_xe] = xe;
    reio_vector[preio->index_re_dkappadz] = dkappadz;
    reio_vector[preio->index_re_dkappadeta] = dkappadz * pvecback_th[pba->index_bg_H];

    /** - get baryon temperature **/
    Tb = Tba2*(1+z)*(1+z);
    reio_vector[preio->index_re_Tb] = Tb;
    
    /** - get baryon sound speed */
    reio_vector[preio->index_re_cb2] = 5./3. * kb_over_mH * (1. + (1./_not4_ - 1.) * Yp + xe * (1.-Yp)) * Tb;

    if (gt_add(&gTable,_GT_END_,(void *) reio_vector,sizeof(double)*(preio->re_size)) == _FAILURE_) { 
      sprintf(pth->error_message,"%s(L:%d) : error in gt_add() \n=>%s",__func__,__LINE__,gTable.error_message);
      return _FAILURE_;
    }
    number_of_redshifts++;
  }
  
  /** - allocate reionization_table with correct size */
  preio->reionization_table=malloc(preio->re_size*number_of_redshifts*sizeof(double));
  if (preio->reionization_table==NULL) {
    sprintf(pth->error_message,"%s(L:%d): Cannot allocate preio->reionization_table \n",__func__,__LINE__);
    return _FAILURE_;
  }
  preio->rt_size=number_of_redshifts;

  /** - retrieve data stored in the growTable with gt_getPtr() */
  if (gt_getPtr(&gTable,(void**)&pData) == _FAILURE_) {
    sprintf(pth->error_message,"%s(L:%d) : error in gt_getPtr()\n=>%s",__func__,__LINE__,gTable.error_message);
    return _FAILURE_;
  }

  /** -- copy growTable to reionization_temporary_table (invert order of lines, so that redshift is growing, like in recombination table) */
  for (i=0; i < preio->rt_size; i++) {
    memcopy_result = memcpy(preio->reionization_table+i*preio->re_size,pData+(preio->rt_size-i-1)*preio->re_size,preio->re_size*sizeof(double));
    if (memcopy_result != preio->reionization_table+i*preio->re_size) {
      sprintf(pth->error_message,"%s(L:%d) : Cannot copy data back to reionization_temporary_table",__func__,__LINE__);
      return _FAILURE_;
    }
  }

  /** - free the growTable with gt_free() , free vector of reionization variables */
  gt_free(&gTable);
  free(reio_vector);

  /** - spline \f$ d tau / dz \f$ with respect to z in view of integrating for optical depth */
  if (array_spline(preio->reionization_table,
		   preio->re_size,
		   preio->rt_size,
		   preio->index_re_z,
		   preio->index_re_dkappadz,
		   preio->index_re_d3kappadz3,
		   _SPLINE_EST_DERIV_,
		   errmsg) == _FAILURE_) {
    sprintf(pth->error_message,"%s(L:%d) : error in array_spline()\n=>%s",__func__,__LINE__,errmsg);
    return _FAILURE_;
  }
  
  /** - integrate for optical depth */
  if (array_integrate_all_spline(preio->reionization_table,
				 preio->re_size,
				 preio->rt_size,
				 preio->index_re_z,
				 preio->index_re_dkappadz,
				 preio->index_re_d3kappadz3,
				 &(preio->reionization_optical_depth),
				 errmsg) == _FAILURE_) {
    sprintf(pth->error_message,"%s(L:%d) : error in array_integrate_all_spline()\n=>%s",__func__,__LINE__,errmsg);
    return _FAILURE_;
  }

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
 * @param Input/Ouput: pointer to recombination structure
 * @return the error status
 */
int thermodynamics_recombination(struct recombination * preco) {

  /* vector of variables to be integrated: xH, xHe, Tmat */
  double y[3],dy[3];

  double OmegaB,Yp,zinitial,zfinal,x_H0,x_He0;
  double z,mu_H,n,Lalpha,Lalpha_He,DeltaB,DeltaB_He,mu_T;
  double x0,w0,w1,Lw0,Lw1,hW;
  double zstart,zend,rhs,Trad,Tmat;
  int i,Nz;
  Generic_integrator_struct generic_integrator_in;
  
  /** Summary: */

  /** - allocate memory for thermodynamics interpolation tables (size known in advance) */
  preco->rt_size = ppr->recfast_Nz0;
  preco->recombination_table = malloc(preco->re_size*preco->rt_size*sizeof(double));
  if (preco->recombination_table==NULL) {
    sprintf(pth->error_message,"%s(L:%d): Cannot allocate recombination_table",__func__,__LINE__);
    return _FAILURE_;
  }

  /** - initialize generic integrator with initialize_generic_integrator() */
  if (initialize_generic_integrator(_RECFAST_INTEG_SIZE_, &generic_integrator_in) == _FAILURE_){
    sprintf(pth->error_message,"%s(L:%d) : error in initialize_generic_integrator() \n=>%s",__func__,__LINE__,generic_integrator_in.error_message);
    return _FAILURE_;
  }
  
  /** - read a few precision/cosmological parameters */

  /* Nz */
  Nz=ppr->recfast_Nz0;

  /* Ho : must be h in 100km/s/Mpc * bigH */
  H0 = pba->H0 * _bigH_;
  /* printf("h=H0/(100km/s/Mpc)= %f, H0= %e SI \n",pba->.H0*2999.7,H0);*/

  /* Omega_b */
  OmegaB = pba->Omega0_b;
  /*printf("Omega_b = %f \n",OmegaB);*/

  /* Yp */
  Yp = pth->YHe;
  /*printf("Y_He = %f \n",Yp);*/

  /* Tnow */
  Tnow = pth->Tcmb;
  /*printf("T_cmb = %f \n",Tnow);*/

  /* z_initial and z_final */
  zinitial=ppr->recfast_z_initial;
  zfinal=ppr->recfast_z_final;

  /* H_frac */ 
  H_frac = ppr->recfast_H_frac;

  /* fudge */
  fu = ppr->recfast_fudge;
 
  /* related quantities */ 
  z=zinitial;
  mu_H = 1./(1.-Yp);
  mu_T = _not4_ /(_not4_ - (_not4_-1.)*Yp); /* recfast 1.4*/
  kb_over_mH= _k_B_ / ( _C_ * _C_ * _m_H_ );
  fHe = Yp/(_not4_ *(1.-Yp)); /* recfast 1.4 */
  Nnow = 3.*H0*H0*OmegaB/(8.*_PI_*_G_*mu_H*_m_H_);
  /*  printf("Nnow= %e\n",Nnow); */

  /* quantities related to constants defined in thermodynamics.h */
  n = Nnow * pow((1.+z),3);
  Lalpha = 1./_L_H_alpha_;
  Lalpha_He = 1./_L_He_2p_;
  DeltaB = _h_P_*_C_*(_L_H_ion_-_L_H_alpha_);
  CDB = DeltaB/_k_B_;
  DeltaB_He = _h_P_*_C_*(_L_He1_ion_-_L_He_2s_);
  CDB_He = DeltaB_He/_k_B_;
  CB1 = _h_P_*_C_*_L_H_ion_/_k_B_;
  CB1_He1 = _h_P_*_C_*_L_He1_ion_/_k_B_;
  CB1_He2 = _h_P_*_C_*_L_He2_ion_/_k_B_;
  CR = 2.*_PI_*(_m_e_/_h_P_)*(_k_B_/_h_P_);
  CK = pow(Lalpha,3)/(8.*_PI_);
  CK_He = pow(Lalpha_He,3)/(8.*_PI_);
  CL = _C_*_h_P_/(_k_B_*Lalpha);
  CL_He = _C_*_h_P_/(_k_B_/_L_He_2s_);
  CT = (8./3.)*(_sigma_/(_m_e_*_C_))*_a_;
  Bfact = _h_P_*_C_*(_L_He_2p_-_L_He_2s_)/_k_B_;

  C1P3P = _C2p1P_-_C2p3P_;
  cc3P1P = _C2p3P_/_C2p1P_; 
  cccP=pow(cc3P1P,3);
  hck=_h_P_*_C_/_k_B_; 

  if (zinitial < 8000.) {
    sprintf(pth->error_message,"%s(L:%d): z_initial=%f<8000, should get recfast initial conditions from get_init()",__func__,__LINE__,zinitial);
    return _FAILURE_;
  }
  /* uncomment this part only if the initial redshift is not larger than 8000 */
  /*   if (get_init(z,y)== _FAILURE_) { */
  /*     sprintf(pth->error_message,"%s(L:%d): error in calling get_init()",__func__,__LINE__); */
  /*     return _FAILURE_; */
  /*   } */

  /** - impose initial conditions */
  y[0] = 1.;
  y[1] = 1.;
  x0 = 1.+2.*fHe;
  y[2] = Tnow*(1.+z);

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
      x0 = 1.+2.*fHe;
      y[0] = x_H0;
      y[1] = x_He0;
      y[2] = Tnow*(1.+z);
    }
    else 
      if (z > 5000.) {
	x_H0 = 1.;
	x_He0 = 1.;
	rhs = exp( 1.5*log(CR*Tnow/(1.+z)) - CB1_He2/(Tnow*(1.+z)) ) / Nnow;
	rhs = rhs*1.;
	x0 = 0.5*(sqrt(pow((rhs-1.-fHe),2) + 4.*(1.+2.*fHe)*rhs) - (rhs-1.-fHe));
	y[0] = x_H0;
	y[1] = x_He0;
	y[2] = Tnow*(1.+z);
      }
      else 
	if (z > 3500.) {
	  x_H0 = 1.;
	  x_He0 = 1.;
	  x0 = x_H0 + fHe*x_He0;
	  y[0] = x_H0;
	  y[1] = x_He0;
	  y[2] = Tnow*(1.+z);
	}
	else 
	  if (y[1] > ppr->recfast_x_He0_trigger) {
	    x_H0 = 1.;
	    rhs = exp(1.5*log(CR*Tnow/(1.+z)) - CB1_He1/(Tnow*(1.+z)))/Nnow;
	    rhs = rhs*4.;
	    x_He0 = 0.5*(sqrt(pow((rhs-1.),2) + 4.*(1.+fHe)*rhs )- (rhs-1.));
	    x0 = x_He0;
	    x_He0 = (x0-1.)/fHe;
	    y[0] = x_H0;
	    y[1] = x_He0;
	    y[2] = Tnow*(1.+z);
	  }
	  else 
	    if (y[0] > ppr->recfast_x_H0_trigger) {
	      rhs = exp(1.5*log(CR*Tnow/(1.+z)) - CB1/(Tnow*(1.+z)))/Nnow;
	      x_H0 = 0.5*(sqrt(pow(rhs,2)+4.*rhs) - rhs);

	      if (generic_integrator(thermodynamics_derivs_with_recfast,
				     zstart,
				     zend,
				     y,
				     ppr->tol_thermo_integration,
				     &generic_integrator_in)== _FAILURE_) {
		sprintf(Transmit_Error_Message,"%s",pth->error_message);
		sprintf(pth->error_message,"%s(L:%d) : error in generic_integrator()\n=>%s\n=>%s",__func__,__LINE__,generic_integrator_in.error_message,Transmit_Error_Message);
		return _FAILURE_;
	      }

	      y[0] = x_H0;
	      x0 = y[0] + fHe*y[1];

	      

	    }      
	    else {

	      if (generic_integrator(thermodynamics_derivs_with_recfast,
				     zstart,
				     zend,
				     y,
				     ppr->tol_thermo_integration,
				     &generic_integrator_in)== _FAILURE_) {
		sprintf(Transmit_Error_Message,"%s",pth->error_message);
		sprintf(pth->error_message,"%s(L:%d) : error in generic_integrator()\n=>%s\n=>%s",__func__,__LINE__,generic_integrator_in.error_message,Transmit_Error_Message);
		return _FAILURE_;
	      }

	      x0 = y[0] + fHe*y[1];

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
    sprintf(pth->error_message,"");                             // initialise the error_message to "".
    thermodynamics_derivs_with_recfast(zend, y, dy);            // run a 'void' function -> we don't know if a problem occured in this function.
    if(strcmp(pth->error_message,"")){                          // check if the error_message is still "" to look if the function worked corectly.
      sprintf(Transmit_Error_Message,"%s",pth->error_message);
      sprintf(pth->error_message,"%s(L:%d) : error in thermodynamics_derivs_with_recfast()\n=>%s",__func__,__LINE__,Transmit_Error_Message);
      return _FAILURE_;
    }

    /* -> cb2 = (k_B/mu) Tb (1-1/3 dlnTb/dlna) = (k_B/mu) Tb (1+1/3 (1+z) dlnTb/dz) */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_cb2)
      = kb_over_mH * (1. + (1./_not4_ - 1.) * Yp + x0 * (1.-Yp)) * y[2] * (1. + (1.+zend) * dy[2] / y[2] / 3.);

    /* -> dkappa/deta = a n_e x_e sigma_T = a^{-2} n_e(today) x_e sigma_T (in units of 1/Mpc) */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_dkappadeta)
      = (1.+zend) * (1.+zend) * Nnow * x0 * _sigma_ * _Mpc_over_m_;
    
  }

  /* RECFAST is done */

  /** - cleanup generic integrator with cleanup_generic_integrator() */

  if (cleanup_generic_integrator(&generic_integrator_in) == _FAILURE_){
    sprintf(pth->error_message,"%s(L:%d) : error in cleanup_generic_integrator()\n=>%s",__func__,__LINE__,generic_integrator_in.error_message);
    return _FAILURE_;
  }
  
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
 * \f$ d x_H / dz, d x_{He} / dz, d T_{mat} / dz \f$. Passed as an argument to
 * the generic_integrator() function.
 *
 * @param z Input : redshift
 * @param y Input : vector of variable
 * @param dy Output : its derivative (already allocated)
 */
void thermodynamics_derivs_with_recfast(
					double z,
					double * y,
					double * dy
					) {

  double x,n,n_He,Trad,Tmat,x_H,x_He,Hz;
  double Rup,Rdown,K,K_He,Rup_He,Rdown_He,He_Boltz;
  double timeTh,timeH;
  double sq_0,sq_1,a_PPB,b_PPB,c_PPB,d_PPB;

  double cc3P,tau;
  int L,L2,LL;

  /* new in recfast 1.4: */
  double Rdown_trip,Rup_trip,tauHe_s,pHe_s,Doppler,gamma_2Ps,pb,qb,AHcon;
  double tauHe_t,pHe_t,CfHe_t,CL_PSt,gamma_2Pt;
  short Heflag;

  x_H = y[0];
  x_He = y[1];
  x = x_H + fHe * x_He;
  Tmat = y[2];

  n = Nnow * pow((1.+z),3);
  n_He = fHe * Nnow * pow((1.+z),3);
  Trad = Tnow * (1.+z);

  if (background_functions_of_a(1./(1.+z),short_info,pvecback_th) == _FAILURE_){
    sprintf(pth->error_message,"%s(L:%d) : error in background_functions_of_a()\n=>%s",__func__,__LINE__,pba->error_message);
    return;
  }
  
  Hz=pvecback_th[pba->index_bg_H]/_Mpc_in_sec_;

  Rdown=1.e-19*_a_PPB_*pow((Tmat/1.e4),_b_PPB_)/(1.+_c_PPB_*pow((Tmat/1.e4),_d_PPB_));
  Rup = Rdown * pow((CR*Tmat),1.5)*exp(-CDB/Tmat);

  sq_0 = sqrt(Tmat/_T_0_);
  sq_1 = sqrt(Tmat/_T_1_);
  Rdown_He = _a_VF_/(sq_0 * pow((1.+sq_0),(1.-_b_VF_)) * pow((1. + sq_1),(1. + _b_VF_)));
  Rup_He = Rdown_He*pow((CR*Tmat),1.5)*exp(-CDB_He/Tmat);
  Rup_He = 4.*Rup_He;
  K = CK/Hz;

  /* following is from recfast 1.4 */

  Rdown_trip = _a_trip_/(sq_0*pow((1.+sq_0),(1.-_b_trip_)) * pow((1.+sq_1),(1.+_b_trip_)));
  Rup_trip = Rdown_trip*exp(-_h_P_*_C_*_L_He2St_ion_/(_k_B_*Tmat))*pow(CR*Tmat,1.5)*4./3.;

  if ((x_He < 5.e-9) || (x_He > 0.980)) 
    Heflag = 0;
  else
    Heflag = ppr->recfast_Heswitch;

  if (Heflag == 0)
    K_He = CK_He/Hz;
  else {
    tauHe_s = _A2P_s_*CK_He*3.*n_He*(1.-x_He)/Hz;
    pHe_s = (1.-exp(-tauHe_s))/tauHe_s;
    K_He = 1./(_A2P_s_*pHe_s*3.*n_He*(1.-x_He));

    if (((Heflag == 2) || (Heflag >= 5)) && (x_H < 0.99999)) {
      Doppler = 2.*_k_B_*Tmat/(_m_H_*_not4_*_C_*_C_);
      Doppler = _C_*_L_He_2p_*sqrt(Doppler);
      gamma_2Ps = 3.*_A2P_s_*fHe*(1.-x_He)*_C_*_C_
	/(sqrt(_PI_)*_sigma_He_2Ps_*8.*_PI_*Doppler*(1.-x_H))
	/pow(_C_*_L_He_2p_,2.);
      pb = 0.36;
      qb = ppr->recfast_fudge_He;
      AHcon = _A2P_s_/(1.+pb*pow(gamma_2Ps,qb));
      K_He=1./((_A2P_s_*pHe_s+AHcon)*3.*n_He*(1.-x_He));
    }

    if (Heflag >= 3) {
      tauHe_t = _A2P_t_*n_He*(1.-x_He)*3./(8.*_PI_*Hz*pow(_L_He_2Pt_,3.));
      pHe_t = (1. - exp(-tauHe_t))/tauHe_t;
      CL_PSt = _h_P_*_C_*(_L_He_2Pt_ - _L_He_2St_)/_k_B_;
      if ((Heflag == 3) || (Heflag == 5) || (x_H >= 0.99999)) {
	CfHe_t = _A2P_t_*pHe_t*exp(-CL_PSt/Tmat);
	CfHe_t = CfHe_t/(Rup_trip+CfHe_t);
      }
      else {
	Doppler = 2.*_k_B_*Tmat/(_m_H_*_not4_*_C_*_C_);
	Doppler = _C_*_L_He_2Pt_*sqrt(Doppler);
	gamma_2Pt = 3.*_A2P_t_*fHe*(1.-x_He)*_C_*_C_
	  /(sqrt(_PI_)*_sigma_He_2Pt_*8.*_PI_*Doppler*(1.-x_H))
	  /pow(_C_*_L_He_2Pt_,2.);
	pb = 0.66;
	qb = 0.9;
	AHcon = _A2P_t_/(1.+pb*pow(gamma_2Pt,qb))/3.;
	CfHe_t = (_A2P_t_*pHe_t+AHcon)*exp(-CL_PSt/Tmat);
	CfHe_t = CfHe_t/(Rup_trip+CfHe_t);
      }
    }
  }

  /* end of new recfast 1.4 piece */

  timeTh=(1./(CT*pow(Trad,4.)))*(1.+x+fHe)/x;
  timeH=2./(3.*H0*pow(1.+z,1.5)); 

  if (x_H > 0.99)
    dy[0] = 0.;
  else {
    if (x_H > 0.985) {
      dy[0] = (x*x_H*n*Rdown - Rup*(1.-x_H)*exp(-CL/Tmat)) /(Hz*(1.+z));
    }
    else {
      dy[0] = ((x*x_H*n*Rdown - Rup*(1.-x_H)*exp(-CL/Tmat)) *(1. + K*_Lambda_*n*(1.-x_H))) /(Hz*(1.+z)*(1./fu+K*_Lambda_*n*(1.-x)/fu +K*Rup*n*(1.-x)));
    }
  }

  if (x_He < 1.e-15) 
    dy[1]=0.;
  else {

    if (Bfact/Tmat < 680.) 
      He_Boltz=exp(Bfact/Tmat);
    else 
      He_Boltz=exp(680.);

    dy[1] = ((x*x_He*n*Rdown_He - Rup_He*(1.-x_He)*exp(-CL_He/Tmat)) *(1. + K_He*_Lambda_He_*n_He*(1.-x_He)*He_Boltz)) /(Hz*(1+z) * (1. + K_He*(_Lambda_He_+Rup_He)*n_He*(1.-x_He)*He_Boltz));

    if (Heflag == 3)
      dy[1] = dy[1] + 
	(x*x_He*n*Rdown_trip
	 - (1.-x_He)*3.*Rup_trip*exp(-_h_P_*_C_*_L_He_2St_/(_k_B_*Tmat)))
	*CfHe_t/(Hz*(1.+z));
	
  }

  if (timeTh < H_frac*timeH)
    dy[2]=Tmat/(1.+z);
  else
    dy[2]= CT * pow(Trad,4) * x / (1.+x+fHe) * (Tmat-Trad) / (Hz*(1.+z)) + 2.*Tmat/(1.+z);
  return;
}      

/* int get_init( */
/* 	     double z,	     /\**< Input  : redshift *\/ */
/* 	     double * y      /\**< Ouput : x_H0, x_He0, T *\/ */
/* 	     ) { */
  
/*   double rhs,x_H0,x_He0,x0; */
  
/*   if(z > 8000.) { */
/*     x_H0 = 1.; */
/*     x_He0 = 1.; */
/*     x0 = 1.+2.*fHe; */
/*   } */
/*   else { */
/*     if (z > 3500.) { */
/*       x_H0 = 1.; */
/*       x_He0 = 1.; */
/*       rhs = exp( 1.5 * log(CR*Tnow/(1.+z)) - CB1_He2/(Tnow*(1.+z)) ) / Nnow; */
/*       rhs = rhs*1.;  /\*ratio of g's is 1 for He++ <-> He+*\/ */
/*       x0 = 0.5 * ( sqrt( pow((rhs-1.-fHe),2) + 4.*(1.+2.*fHe)*rhs) - (rhs-1.-fHe) ); */
/*     } */
/*     else { */
/*       if(z > 2000.) { */
/*         x_H0 = 1.; */
/* 	rhs = exp( 1.5 * log(CR*Tnow/(1.+z)) - CB1_He1/(Tnow*(1.+z)) ) / Nnow; */
/*         rhs = rhs*4.;  /\*ratio of g's is 4 for He+ <-> He0*\/ */
/* 	x_He0 = 0.5  * ( sqrt( pow((rhs-1.),2) + 4.*(1.+fHe)*rhs )- (rhs-1.)); */
/* 	x0 = x_He0; */
/* 	x_He0 = (x0 - 1.)/fHe; */
/*       } */
/*       else { */
/* 	rhs = exp( 1.5 * log(CR*Tnow/(1.+z)) - CB1/(Tnow*(1.+z)) ) / Nnow; */
/* 	x_H0 = 0.5 * (sqrt( pow(rhs,2)+4.*rhs ) - rhs ); */
/* 	x_He0 = 0.; */
/* 	x0 = x_H0; */
/*       } */
/*     } */
/*   } */

/*   y[0]=x_H0; */
/*   y[1]=x_He0; */
/*   y[2]=Tnow*(1.+z); */

/*   return _SUCCESS_; */
/* } */
