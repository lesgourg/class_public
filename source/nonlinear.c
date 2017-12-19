/** @file nonlinear.c Documented nonlinear module
 *
 * Julien Lesgourgues, 6.03.2014
 *
 * New module replacing an older one present up to version 2.0 The new
 * module is located in a better place in the main, allowing it to
 * compute non-linear correction to \f$ C_l\f$'s and not just \f$ P(k)\f$. It will
 * also be easier to generalize to new methods.  The old implementation
 * of one-loop calculations and TRG calculations has been dropped from
 * this version, they can still be found in older versions.
 *
 */

#include "nonlinear.h"
#include "extrapolate_source.h"

int nonlinear_k_nl_at_z(
                        struct background *pba,
                        struct nonlinear * pnl,
                        double z,
                        double * k_nl
                        ) {

  double tau;

  class_call(background_tau_of_z(pba,
                                 z,
                                 &tau),
             pba->error_message,
             pnl->error_message);

  if (pnl->tau_size == 1) {
    *k_nl = pnl->k_nl[0];
  }
  else {
    class_call(array_interpolate_two(pnl->tau,
                                     1,
                                     0,
                                     pnl->k_nl,
                                     1,
                                     pnl->tau_size,
                                     tau,
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
                   struct primordial *ppm,
                   struct nonlinear *pnl
                   ) {

  int index_ncdm;
  int index_k;
  int index_tau;
  int size_extrapolated_source;
  double *pk_l;
  double *pk_nl;
  double *lnk_l;
  double *lnpk_l;
  double *ddlnpk_l;
  short print_warning=_FALSE_;
  double * pvecback;
  int last_index;
  double a,z;
  enum halofit_statement halofit_found_k_max;

  /** Summary
   *
   * (a) First deal with the case where non non-linear corrections requested */

  if (pnl->method == nl_none) {
    if (pnl->nonlinear_verbose > 0)
      printf("No non-linear spectra requested. Nonlinear module skipped.\n");
  }

  /** (b) Compute for HALOFIT non-linear spectrum */

  else if ((pnl->method == nl_halofit) || ((pnl->method == nl_HMcode))) {
    if ((pnl->nonlinear_verbose > 0) && (pnl->method == nl_halofit))
      printf("Computing non-linear matter power spectrum with Halofit (including update Takahashi et al. 2012 and Bird 2014)\n");
	
	if ((pnl->nonlinear_verbose > 0) && (pnl->method == nl_HMcode))
      printf("Computing non-linear matter power spectrum with HMcode \n");
	
    if (pba->has_ncdm) {
      for (index_ncdm=0;index_ncdm < pba->N_ncdm; index_ncdm++){
        if (pba->m_ncdm_in_eV[index_ncdm] >  _M_EV_TOO_BIG_FOR_HALOFIT_)
          fprintf(stdout,"Warning: Halofit is proved to work for CDM, and also with a small HDM component thanks to Bird et al.'s update. But it sounds like you are running with a WDM component of mass %f eV, which makes the use of Halofit suspicious.\n",pba->m_ncdm_in_eV[index_ncdm]);
      }
    }

    /** - copy list of (k,tau) from perturbation module */

    pnl->k_size = ppt->k_size[ppt->index_md_scalars];
    class_alloc(pnl->k,pnl->k_size*sizeof(double),pnl->error_message);
    for (index_k=0; index_k<pnl->k_size; index_k++)
      pnl->k[index_k] = ppt->k[ppt->index_md_scalars][index_k];

    pnl->tau_size = ppt->tau_size;
    class_alloc(pnl->tau,pnl->tau_size*sizeof(double),pnl->error_message);
    for (index_tau=0; index_tau<pnl->tau_size; index_tau++)
      pnl->tau[index_tau] = ppt->tau_sampling[index_tau];

    class_alloc(pnl->nl_corr_density,pnl->tau_size*pnl->k_size*sizeof(double),pnl->error_message);
    class_alloc(pnl->k_nl,pnl->tau_size*sizeof(double),pnl->error_message);

    class_alloc(pk_l,pnl->k_size*sizeof(double),pnl->error_message);
    class_alloc(pk_nl,pnl->k_size*sizeof(double),pnl->error_message);

    class_alloc(lnk_l,pnl->k_size*sizeof(double),pnl->error_message);
    class_alloc(lnpk_l,pnl->k_size*sizeof(double),pnl->error_message);
    class_alloc(ddlnpk_l,pnl->k_size*sizeof(double),pnl->error_message);


		if (pnl->method == nl_HMcode){
			class_call(get_extrapolated_source_size(
																						ppr->k_per_decade_for_pk,
																						pnl->k[pnl->k_size-1],
																						ppr->hmcode_max_k_extra,
																						pnl->k_size,
																						&size_extrapolated_source,
																						pnl->error_message),
														pnl->error_message,
														pnl->error_message);
			pnl->k_size_extra = size_extrapolated_source;
   
			class_alloc(pnl->k_extra,pnl->k_size_extra*sizeof(double),pnl->error_message);
		
			class_realloc(pk_l,pk_l,pnl->k_size_extra*sizeof(double),pnl->error_message);
   
		
		
			class_call(extrapolate_k(
														 pnl->k,
														 pnl->k_size,
														 pnl->k_extra,
														 ppr->k_per_decade_for_pk,
														 ppr->hmcode_max_k_extra,
														 pnl->error_message),
										pnl->error_message,
										pnl->error_message);
		}

    /** - loop over time */

    for (index_tau = pnl->tau_size-1; index_tau>=0; index_tau--) {

      /* get P_L(k) at this time */
			class_call(nonlinear_pk_l(pba,ppt,ppm,pnl,index_tau,pk_l,lnk_l,lnpk_l,ddlnpk_l),
                 pnl->error_message,
                 pnl->error_message);

       /* get P_NL(k) at this time */
      if (pnl->method == nl_halofit) {
				if (print_warning == _FALSE_) {
					class_call(nonlinear_halofit(ppr,
																			pba,
																			ppm,
																			pnl,
																			pnl->tau[index_tau],
																			pk_l,
																			pk_nl,
																			lnk_l,
																			lnpk_l,
																			ddlnpk_l,
																			&(pnl->k_nl[index_tau]),
																			&halofit_found_k_max),
										pnl->error_message,
										pnl->error_message);

					if (halofit_found_k_max == ok) {

						// for debugging:
						/*
							for (index_k=0; index_k<pnl->k_size; index_k++) {
							fprintf(stdout,"%e  %e  %e\n",pnl->k[index_k],pk_l[index_k],pk_nl[index_k]);
							}
							fprintf(stdout,"\n\n");
						*/

						for (index_k=0; index_k<pnl->k_size; index_k++) {
							pnl->nl_corr_density[index_tau * pnl->k_size + index_k] = sqrt(pk_nl[index_k]/pk_l[index_k]);
						}
					}
					else {
						/* when Halofit found k_max too small, use 1 as the
							non-linear correction for this redshift/time, store the
							last index which worked, and print a warning. */
						print_warning = _TRUE_;
						pnl->index_tau_min_nl = index_tau+1;
						for (index_k=0; index_k<pnl->k_size; index_k++) {
							pnl->nl_corr_density[index_tau * pnl->k_size + index_k] = 1.;
						}
						if (pnl->nonlinear_verbose > 0) {
							class_alloc(pvecback,pba->bg_size*sizeof(double),pnl->error_message);
							class_call(background_at_tau(pba,pnl->tau[index_tau],pba->short_info,pba->inter_normal,&last_index,pvecback),
												pba->error_message,
												pnl->error_message);
							a = pvecback[pba->index_bg_a];
							z = pba->a_today/a-1.;
							fprintf(stdout,
											" -> [WARNING:] Halofit non-linear corrections could not be computed at redshift z=%5.2f and higher.\n    This is because k_max is too small for Halofit to be able to compute the scale k_NL at this redshift.\n    If non-linear corrections at such high redshift really matter for you,\n    just try to increase one of the parameters P_k_max_h/Mpc or P_k_max_1/Mpc or halofit_min_k_max (the code will take the max of these parameters) until reaching desired z.\n",
											z);
							free(pvecback);
						}
					}
				}
				else {
					/* if Halofit found k_max too small at a previous
							time/redhsift, use 1 as the non-linear correction for all
							higher redshifts/earlier times. */
					for (index_k=0; index_k<pnl->k_size; index_k++) {
						pnl->nl_corr_density[index_tau * pnl->k_size + index_k] = 1.;
					}
				}
			}
			else if (pnl->method == nl_HMcode) {
				 /* only fill sigma and grow table once (at redshift 0) */
				int i;
				if (index_tau==pnl->tau_size-1){
					class_call(nonlinear_fill_growtab(ppr,pba,pnl), 
						pnl->error_message, pnl->error_message);	      
					class_call(nonlinear_fill_sigtab(ppr,pba,ppt,ppm,pnl,pk_l), 
						pnl->error_message, pnl->error_message);
					/*if	(index_tau == pnl->tau_size-1) {
						fprintf(stdout, "i,  R         sigma\n");
						for (i=0;i<64;i++){
							fprintf(stdout, "%d, %e, %e\n",i, pnl->rtab[i]*pba->h, pnl->stab[i]);
						}
					}*/       
        } 
				/** Compute prediction using HMcode at conformal time tau */
				class_call(nonlinear_HMcode(ppr,
																		pba,
																		ppt,
																		ppm,
																		pnl,
																		pnl->tau[index_tau],
																		pk_l,
																		pk_nl,
																		lnk_l,
																		lnpk_l,
																		ddlnpk_l,
																		&(pnl->k_nl[index_tau])),
															pnl->error_message,
															pnl->error_message);
				for (index_k=0; index_k<pnl->k_size; index_k++) {
					pnl->nl_corr_density[index_tau * pnl->k_size + index_k] = sqrt(pk_nl[index_k]/pk_l[index_k]);
				}
      }
    }
		// for debugging
    /*
    for (index_tau = pnl->tau_size-1; index_tau>=0; index_tau--) {
      for (index_k=0; index_k<pnl->k_size; index_k++) {
        fprintf(stdout,"%e  %e\n",pnl->k[index_k],pnl->nl_corr_density[index_tau * pnl->k_size + index_k]);
      }
      fprintf(stdout,"\n\n");
    }
    */

    free(pk_l);
    free(pk_nl);
    free(lnk_l);
    free(lnpk_l);
    free(ddlnpk_l);
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

    if (pnl->method == nl_halofit) {
      free(pnl->k);
      free(pnl->tau);
      free(pnl->nl_corr_density);
      free(pnl->k_nl);
    }
		else if (pnl->method == nl_HMcode){
			free(pnl->k);
      free(pnl->tau);
      free(pnl->nl_corr_density);
      free(pnl->k_nl);
      free(pnl->k_extra);
			free(pnl->rtab);
			free(pnl->stab);
			free(pnl->ddstab);
			free(pnl->sigtab);
			free(pnl->growtable);
			free(pnl->tautable);
			free(pnl->ztable);
		}
  }

  if (pnl->has_pk_eq == _TRUE_) {
    free(pnl->eq_tau);
    free(pnl->eq_w_and_Omega);
    free(pnl->eq_ddw_and_ddOmega);
  }

  return _SUCCESS_;

}

int nonlinear_pk_l(
                   struct background *pba,
                   struct perturbs *ppt,
                   struct primordial *ppm,
                   struct nonlinear *pnl,
                   int index_tau,
                   double *pk_l,
                   double *lnk,
                   double *lnpk,
                   double *ddlnpk) {

  int index_md;
  int index_k;
  int index_ic;
  int index_ic1,index_ic2,index_ic1_ic2;
  double * primordial_pk;
  double source_ic1,source_ic2;
  double * source_ic_extra;

  index_md = ppt->index_md_scalars;

  class_alloc(primordial_pk,ppm->ic_ic_size[index_md]*sizeof(double),pnl->error_message);
	
	if (pnl->method == nl_HMcode){
		
		class_alloc(source_ic_extra,ppm->ic_size[index_md]*pnl->k_size_extra*sizeof(double),pnl->error_message);
		
		for (index_ic=0; index_ic<ppm->ic_size[index_md]; index_ic++){

			class_call(extrapolate_source(
																		pnl->k_extra,
																		pnl->k_size,
																		pnl->k_size_extra,
																		ppt->sources[index_md][index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_m]+index_tau * pnl->k_size,
																		extrapolation_only_max_units,
																		source_ic_extra+index_ic*pnl->k_size_extra,
																		pba->a_eq*pba->H_eq,
																		pba->h,																		
																		pnl->error_message), 
															pnl->error_message,
															pnl->error_message)
			
		}

		for (index_k=0; index_k<pnl->k_size_extra; index_k++) {
			//fprintf(stdout, "%e, %e\n", pnl->k_extra[index_k], source_ic_extra[index_k]);
    
			class_call(primordial_spectrum_at_k(ppm,
                                        index_md,
                                        linear,
                                        pnl->k_extra[index_k],
                                        primordial_pk),
               ppm->error_message,
               pnl->error_message);

			pk_l[index_k] = 0;
			
				// part diagonal in initial conditions 
			for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++) {

					index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic1,ppm->ic_size[index_md]);

			
			pk_l[index_k] += 2.*_PI_*_PI_/pow(pnl->k_extra[index_k],3)\
					*source_ic_extra[index_ic1*pnl->k_size_extra+index_k]*source_ic_extra[index_ic1*pnl->k_size_extra+index_k]\
					*primordial_pk[index_ic1_ic2];
			}
		
			// part non-diagonal in initial conditions 
			for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++) {
				for (index_ic2 = index_ic1+1; index_ic2 < ppm->ic_size[index_md]; index_ic2++) {

					index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,ppm->ic_size[index_md]);

					if (ppm->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {
					
						pk_l[index_k] += 2.*2.*_PI_*_PI_/pow(pnl->k_extra[index_k],3)
							*source_ic_extra[index_ic1*pnl->k_size_extra+index_k]*source_ic_extra[index_ic2*pnl->k_size_extra+index_k]
							*primordial_pk[index_ic1_ic2]; // extra 2 factor (to include the symmetric term ic2,ic1)					
					}
        }
      }
			if (index_k<pnl->k_size) {
				lnk[index_k] = log(pnl->k[index_k]);
				lnpk[index_k] = log(pk_l[index_k]);
			}
			// for debugging:
			//if (index_tau == pnl->tau_size-1) fprintf(stdout, "%e %e\n", pnl->k_extra[index_k], pk_l[index_k]);
    
		}
		free(source_ic_extra);
		
		
	} else  { 
			for (index_k=0; index_k<pnl->k_size; index_k++) {

					class_call(primordial_spectrum_at_k(ppm,
                                        index_md,
                                        linear,
                                        pnl->k[index_k],
                                        primordial_pk),
               ppm->error_message,
               pnl->error_message);

				pk_l[index_k] = 0;
				

				// part diagonal in initial conditions 
				for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++) {

					index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic1,ppm->ic_size[index_md]);

					source_ic1 = ppt->sources[index_md]
						[index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_m]
						[index_tau * ppt->k_size[index_md] + index_k];

					pk_l[index_k] += 2.*_PI_*_PI_/pow(pnl->k[index_k],3)
						*source_ic1*source_ic1
						*primordial_pk[index_ic1_ic2];
				}

			// part non-diagonal in initial conditions 
				for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++) {
					for (index_ic2 = index_ic1+1; index_ic2 < ppm->ic_size[index_md]; index_ic2++) {

						index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,ppm->ic_size[index_md]);

						if (ppm->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {

							source_ic1 = ppt->sources[index_md]
								[index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_m]
								[index_tau * ppt->k_size[index_md] + index_k];

							source_ic2 = ppt->sources[index_md]
								[index_ic2 * ppt->tp_size[index_md] + ppt->index_tp_delta_m]
								[index_tau * ppt->k_size[index_md] + index_k];

							pk_l[index_k] += 2.*2.*_PI_*_PI_/pow(pnl->k[index_k],3)
								*source_ic1*source_ic2
								*primordial_pk[index_ic1_ic2]; // extra 2 factor (to include the symmetric term ic2,ic1)					
				        						
					}
        }
      }

			lnk[index_k] = log(pnl->k[index_k]);
			lnpk[index_k] = log(pk_l[index_k]);
		}
	}
	
  class_call(array_spline_table_columns(lnk,
                                        pnl->k_size,
                                        lnpk,
                                        1,
                                        ddlnpk,
                                        _SPLINE_NATURAL_,
                                        pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  free(primordial_pk);


  return _SUCCESS_;

}

int nonlinear_halofit(
                      struct precision *ppr,
                      struct background *pba,
                      struct primordial *ppm,
                      struct nonlinear *pnl,
                      double tau,
                      double *pk_l,
                      double *pk_nl,
                      double *lnk_l,
                      double *lnpk_l,
                      double *ddlnpk_l,
                      double *k_nl,
                      enum halofit_statement * halofit_found_k_max
                      ) {

  double Omega_m,Omega_v,fnu,Omega0_m, w0, dw_over_da_fld, integral_fld;

  /** Determine non linear ratios (from pk) **/

  int index_k;
  double pk_lin,pk_quasi,pk_halo,rk;
  double sigma,rknl,rneff,rncur,d1,d2;
  double diff,xlogr1,xlogr2,rmid;

  double gam,a,b,c,xmu,xnu,alpha,beta,f1,f2,f3;
  double pk_linaa;
  double y;
  double f1a,f2a,f3a,f1b,f2b,f3b,frac;

  double * pvecback;

  int last_index=0;
  int counter;
  double sum1,sum2,sum3;
  double anorm;

  double *integrand_array;
  int integrand_size;
  int index_ia_k;
  int index_ia_pk;
  int index_ia_sum;
  int index_ia_ddsum;
  /*
  int index_ia_sum2;
  int index_ia_ddsum2;
  int index_ia_sum3;
  int index_ia_ddsum3;
  */
  int ia_size;
  int index_ia;

  double k_integrand;
  double lnpk_integrand;

  double R;

  double * w_and_Omega;

  class_alloc(pvecback,pba->bg_size*sizeof(double),pnl->error_message);

  Omega0_m = (pba->Omega0_cdm + pba->Omega0_b + pba->Omega0_ncdm_tot + pba->Omega0_dcdm);
  fnu      = pba->Omega0_ncdm_tot/Omega0_m;

  if (pnl->has_pk_eq == _FALSE_) {

    /* default method to compute w0 = w_fld today, Omega_m(tau) and Omega_v=Omega_DE(tau),
       all required by HALFIT fitting formulas */

    class_call(background_w_fld(pba,pba->a_today,&w0,&dw_over_da_fld,&integral_fld), pba->error_message, pnl->error_message);

    class_call(background_at_tau(pba,tau,pba->long_info,pba->inter_normal,&last_index,pvecback),
               pba->error_message,
               pnl->error_message);

    Omega_m = pvecback[pba->index_bg_Omega_m];
    Omega_v = 1.-pvecback[pba->index_bg_Omega_m]-pvecback[pba->index_bg_Omega_r];

  }
  else {

    /* alternative method called PK-equal, described in 0810.0190 and
                      1601.0723, extending the range of validity of
                      HALOFIT from constant w to w0-wa models. In that
                      case, some effective values of w0(tau_i) and
                      Omega_m(tau_i) have been pre-computed in the input
                      module, and we just ned to interpolate within
                      tabulated arrays, to get them at the current tau
                      value. */

    class_alloc(w_and_Omega,pnl->eq_size*sizeof(double),pnl->error_message);

    class_call(array_interpolate_spline(
                                        pnl->eq_tau,
                                        pnl->eq_tau_size,
                                        pnl->eq_w_and_Omega,
                                        pnl->eq_ddw_and_ddOmega,
                                        pnl->eq_size,
                                        tau,
                                        &last_index,
                                        w_and_Omega,
                                        pnl->eq_size,
                                        pnl->error_message),
               pnl->error_message,
               pnl->error_message);

    w0 = w_and_Omega[pnl->index_eq_w];
    Omega_m = w_and_Omega[pnl->index_eq_Omega_m];
    Omega_v = 1.-Omega_m;

    free(w_and_Omega);
  }

  anorm    = 1./(2*pow(_PI_,2));

  /*      Until the 17.02.2015 the values of k used for integrating sigma(R) quantities needed by Halofit where the same as in the perturbation module.
          Since then, we sample these integrals on more values, in order to get more precise integrals (thanks Matteo Zennaro for noticing the need for this).

     We create a temporary integrand_array which columns will be:
     - k in 1/Mpc
     - just linear P(k) in Mpc**3
     - 1/(2(pi**2)) P(k) k**2 exp(-(kR)**2) or 1/(2(pi**2)) P(k) k**2 2 (kR) exp(-(kR)**2) or 1/(2(pi**2)) P(k) k**2 4 (kR)(1-kR) exp(-(kR)**2)
     - second derivative of previous line with spline
  */

  index_ia=0;
  class_define_index(index_ia_k,     _TRUE_,index_ia,1);
  class_define_index(index_ia_pk,    _TRUE_,index_ia,1);
  class_define_index(index_ia_sum,   _TRUE_,index_ia,1);
  class_define_index(index_ia_ddsum, _TRUE_,index_ia,1);
  ia_size = index_ia;

  integrand_size=(int)(log(pnl->k[pnl->k_size-1]/pnl->k[0])/log(10.)*ppr->halofit_k_per_decade)+1;

  class_alloc(integrand_array,integrand_size*ia_size*sizeof(double),pnl->error_message);

  //fprintf(stderr,"Omega_m=%e,  fnu=%e\n",Omega0_m,fnu);

  /* we fill integrand_array with values of k and P(k) using interpolation */

  last_index=0;

  for (index_k=0; index_k < integrand_size; index_k++) {

    k_integrand=pnl->k[0]*pow(10.,index_k/ppr->halofit_k_per_decade);

    class_call(array_interpolate_spline(lnk_l,
                                        pnl->k_size,
                                        lnpk_l,
                                        ddlnpk_l,
                                        1,
                                        log(k_integrand),
                                        &last_index,
                                        &lnpk_integrand,
                                        1,
                                        pnl->error_message),
               pnl->error_message,
               pnl->error_message);

    integrand_array[index_k*ia_size + index_ia_k] = k_integrand;
    integrand_array[index_k*ia_size + index_ia_pk] = exp(lnpk_integrand);

  }

  class_call(background_at_tau(pba,tau,pba->long_info,pba->inter_normal,&last_index,pvecback),
             pba->error_message,
             pnl->error_message);

  Omega_m = pvecback[pba->index_bg_Omega_m];
  Omega_v = 1.-pvecback[pba->index_bg_Omega_m]-pvecback[pba->index_bg_Omega_r];

  // for debugging:
  //printf("Call Halofit at z=%e\n",pba->a_today/pvecback[pba->index_bg_a]-1.);

  /* minimum value of R such that the integral giving sigma_R is
     converged.  The parameter halofit_sigma_precision should be
     understood as follows: we trust our calculation of sigma(R) as
     long as the integral reaches a value k_max such that the factor
     exp(-(Rk_max)**2) is already as low as halofit_sigma_precisio,
     shoing that the integreal is converged.  In practise this
     condition is tested only for R_max, the highest value of R in our
     bisection algorithm. Hence a smaller value of
     halofit_sigma_precision will lead to a more precise halofit
     result at the *highest* redshift at which halofit can make
     computations, at the expense of requiring a larger k_max; but
     this parameter is not relevant for the precision on P_nl(k,z) at
     other redshifts, so there is normally no need to change i
   */

  R=sqrt(-log(ppr->halofit_sigma_precision))/integrand_array[(integrand_size-1)*ia_size + index_ia_k];

  class_call(nonlinear_halofit_integrate(
                                         pnl,
                                         integrand_array,
                                         integrand_size,
                                         ia_size,
                                         index_ia_k,
                                         index_ia_pk,
                                         index_ia_sum,
                                         index_ia_ddsum,
                                         R,
                                         halofit_integral_one,
                                         &sum1
                                         ),
             pnl->error_message,
             pnl->error_message);

  sigma  = sqrt(sum1);

  /* the following error should not stop the code: it will arrive
     inevitably at some large redshift, and then the code should not
     stop, but just give up computing P_NL(k,z). This is why we have a
     special error handling here (using class_test_except and free()
     commands to avoid memory leaks, and calling this whole function
     not through a class_call) */

  /*
  class_test_except(sigma < 1.,
                    pnl->error_message,
                    free(pvecback);free(integrand_array),
                    "Your k_max=%g 1/Mpc is too small for Halofit to find the non-linearity scale z_nl at z=%g. Increase input parameter P_k_max_h/Mpc or P_k_max_1/Mpc",
                    pnl->k[pnl->k_size-1],
                    pba->a_today/pvecback[pba->index_bg_a]-1.);
  */

  if (sigma < 1.) {
    * halofit_found_k_max = too_small;
    free(pvecback);
    free(integrand_array);
    return _SUCCESS_;
  }
  else {
    * halofit_found_k_max = ok;
  }

  xlogr1 = log(R)/log(10.);

  /* maximum value of R in the bisection algorithm leading to the
     determination of R_nl.  For this value we can make a
     conservaitive guess: 1/halofit_min_k_nonlinear, where
     halofit_min_k_nonlinear is the minimum value of k at which we ask
     halofit to give us an estimate of P_nl(k,z). By assumption we
     treat all smaller k's as linear, so we know that
     sigma(1/halofit_min_k_nonlinear) must be <<1 (and if it is not
     the test below will alert us) */

  R=1./ppr->halofit_min_k_nonlinear;

  /* corresponding value of sigma_R */
  class_call(nonlinear_halofit_integrate(
                                         pnl,
                                         integrand_array,
                                         integrand_size,
                                         ia_size,
                                         index_ia_k,
                                         index_ia_pk,
                                         index_ia_sum,
                                         index_ia_ddsum,
                                         R,
                                         halofit_integral_one,
                                         &sum1
                                         ),
             pnl->error_message,
             pnl->error_message);

  sigma  = sqrt(sum1);

  class_test(sigma > 1.,
             pnl->error_message,
             "Your input value for the precision parameter halofit_min_k_nonlinear=%e is too large, such that sigma(R=1/halofit_min_k_nonlinear)=% > 1. For self-consistency, it should have been <1. Decrease halofit_min_k_nonlinear",
             ppr->halofit_min_k_nonlinear,sigma);

  xlogr2 = log(R)/log(10.);

  counter = 0;
  do {
    rmid = pow(10,(xlogr2+xlogr1)/2.0);
    counter ++;

    class_call(nonlinear_halofit_integrate(
                                           pnl,
                                           integrand_array,
                                           integrand_size,
                                           ia_size,
                                           index_ia_k,
                                           index_ia_pk,
                                           index_ia_sum,
                                           index_ia_ddsum,
                                           rmid,
                                           halofit_integral_one,
                                           &sum1
                                           ),
             pnl->error_message,
             pnl->error_message);

    sigma  = sqrt(sum1);

    diff = sigma - 1.0;

    if (diff > ppr->halofit_tol_sigma){
      xlogr1=log10(rmid);
    }
    else if (diff < -ppr->halofit_tol_sigma) {
      xlogr2 = log10(rmid);
    }

    /* The first version of this test woukld let the code continue: */
    /*
    class_test_except(counter > _MAX_IT_,
                      pnl->error_message,
                      free(pvecback);free(integrand_array),
                      "could not converge within maximum allowed number of iterations");
    */
    /* ... but in this situation it sounds better to make it stop and return an error! */
    class_test(counter > _MAX_IT_,
               pnl->error_message,
               "could not converge within maximum allowed number of iterations");

  } while (fabs(diff) > ppr->halofit_tol_sigma);

  /* evaluate all the other integrals at R=rmid */

  class_call(nonlinear_halofit_integrate(
                                         pnl,
                                         integrand_array,
                                         integrand_size,
                                         ia_size,
                                         index_ia_k,
                                         index_ia_pk,
                                         index_ia_sum,
                                         index_ia_ddsum,
                                         rmid,
                                         halofit_integral_two,
                                         &sum2
                                         ),
             pnl->error_message,
             pnl->error_message);

  class_call(nonlinear_halofit_integrate(
                                         pnl,
                                         integrand_array,
                                         integrand_size,
                                         ia_size,
                                         index_ia_k,
                                         index_ia_pk,
                                         index_ia_sum,
                                         index_ia_ddsum,
                                         rmid,
                                         halofit_integral_three,
                                         &sum3
                                         ),
             pnl->error_message,
             pnl->error_message);

  sigma  = sqrt(sum1);
  d1 = -sum2/sum1;
  d2 = -sum2*sum2/sum1/sum1 - sum3/sum1;

  rknl  = 1./rmid;
  rneff = -3.-d1;
  rncur = -d2;

  *k_nl = rknl;

  //fprintf(stderr,"Here\n");

  for (index_k = 0; index_k < pnl->k_size; index_k++){

    rk = pnl->k[index_k];

    if (rk > ppr->halofit_min_k_nonlinear) {

      pk_lin = pk_l[index_k]*pow(pnl->k[index_k],3)*anorm;

      /* in original halofit, this is the beginning of the function halofit() */

      /*SPB11: Standard halofit underestimates the power on the smallest
       * scales by a factor of two. Add an extra correction from the
       * simulations in Bird, Viel,Haehnelt 2011 which partially accounts for
       * this.*/
      /*SPB14: This version of halofit is an updated version of the fit to the massive neutrinos
       * based on the results of Takahashi 2012, (arXiv:1208.2701).
       */
      gam=0.1971-0.0843*rneff+0.8460*rncur;
      a=1.5222+2.8553*rneff+2.3706*rneff*rneff+0.9903*rneff*rneff*rneff+ 0.2250*rneff*rneff*rneff*rneff-0.6038*rncur+0.1749*Omega_v*(1.+w0);
      a=pow(10,a);
      b=pow(10, (-0.5642+0.5864*rneff+0.5716*rneff*rneff-1.5474*rncur+0.2279*Omega_v*(1.+w0)));
      c=pow(10, 0.3698+2.0404*rneff+0.8161*rneff*rneff+0.5869*rncur);
      xmu=0.;
      xnu=pow(10,5.2105+3.6902*rneff);
      alpha=fabs(6.0835+1.3373*rneff-0.1959*rneff*rneff-5.5274*rncur);
      beta=2.0379-0.7354*rneff+0.3157*pow(rneff,2)+1.2490*pow(rneff,3)+0.3980*pow(rneff,4)-0.1682*rncur + fnu*(1.081 + 0.395*pow(rneff,2));

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
      pk_halo=pk_halo/(1+xmu*pow(y,-1)+xnu*pow(y,-2))*(1+fnu*(0.977-18.015*(Omega0_m-0.3)));
      // rk is in 1/Mpc, 47.48and 1.5 in Mpc**-2, so we need an h**2 here (Credits Antonio J. Cuesta)
      pk_linaa=pk_lin*(1+fnu*47.48*pow(rk/pba->h,2)/(1+1.5*pow(rk/pba->h,2)));
      pk_quasi=pk_lin*pow((1+pk_linaa),beta)/(1+pk_linaa*alpha)*exp(-y/4.0-pow(y,2)/8.0);

      pk_nl[index_k] = (pk_halo+pk_quasi)/pow(pnl->k[index_k],3)/anorm;

      /* in original halofit, this is the end of the function halofit() */
    }
    else {
      pk_nl[index_k] = pk_l[index_k];
    }
  }

  free(pvecback);
  free(integrand_array);
  return _SUCCESS_;
}

/* in original halofit, this is equivalent to the function wint() */
int nonlinear_halofit_integrate(
                                struct nonlinear *pnl,
                                double * integrand_array,
                                int integrand_size,
                                int ia_size,
                                int index_ia_k,
                                int index_ia_pk,
                                int index_ia_sum,
                                int index_ia_ddsum,
                                double R,
                                enum halofit_integral_type type,
                                double * sum
                                ) {

  double k,pk,x2,integrand;
  int index_k;
  double anorm = 1./(2*pow(_PI_,2));

    for (index_k=0; index_k < integrand_size; index_k++) {
      k = integrand_array[index_k*ia_size + index_ia_k];
      pk = integrand_array[index_k*ia_size + index_ia_pk];
      x2 = k*k*R*R;

      integrand = pk*k*k*anorm*exp(-x2);
      if (type == halofit_integral_two) integrand *= 2.*x2;
      if (type == halofit_integral_three) integrand *= 4.*x2*(1.-x2);

      integrand_array[index_k*ia_size + index_ia_sum] = integrand;
    }

    /* fill in second derivatives */
    class_call(array_spline(integrand_array,
                            ia_size,
                            integrand_size,
                            index_ia_k,
                            index_ia_sum,
                            index_ia_ddsum,
                            _SPLINE_NATURAL_,
                            pnl->error_message),
               pnl->error_message,
               pnl->error_message);

    /* integrate */
    class_call(array_integrate_all_spline(integrand_array,
                                          ia_size,
                                          integrand_size,
                                          index_ia_k,
                                          index_ia_sum,
                                          index_ia_ddsum,
                                          sum,
                                          pnl->error_message),
               pnl->error_message,
               pnl->error_message);

  return _SUCCESS_;
}
/* #define _E_ exp(1)

int nonlinear_sigma(
                  struct precision * ppr,
                  struct background * pba,
                  struct perturbs * ppt,
                  struct primordial * ppm,
                  struct nonlinear * pnl,
                  double R,
                  double *pk_l,
                  double *lnk_l,                 
                  double * sigma
                  ) {
 
  double * array_for_sigma;
  double * k_array;
  double * pk_l_new;
  double lnk_min, lnk_max;
  int k_size_new;
  int index_num;
  int index_k;
  int index_y;
  int index_ddy;
  int i;

  double t,k,W,x;
  
  // extrapolate k array up to k_max=1e4 h/Mpc  
  
  lnk_max = log(1.e8*pba->h);
  lnk_min = log(ppt->k[ppt->index_md_scalars][pnl->k_size-1]);//log(ppt->k_max_for_pk);
  k_size_new = pnl->k_size+(lnk_max-lnk_min)*(ppr->k_per_decade_for_pk);//2.302585e-01);
  //fprintf(stdout, "%d, %d, %e, %e\n", pnl->k_size, k_size_new, ppt->k_max_for_pk, ppr->k_per_decade_for_pk);
  class_alloc(k_array,
              k_size_new*sizeof(double),
              pnl->error_message); 
  
  class_alloc(pk_l_new,
              k_size_new*sizeof(double),
              pnl->error_message); 
  
  for (i=0;i<k_size_new;i++){
		if (i<pnl->k_size-1){
			k_array[i] = exp(lnk_l[i]);
			pk_l_new[i] = pk_l[i];
		}else{
			k_array[i] = exp(lnk_min+(lnk_max-lnk_min)*(i-pnl->k_size)/(k_size_new-1.-pnl->k_size));
			pk_l_new[i] = pk_l[pnl->k_size-2]*pow(log(k_array[i])/log(k_array[pnl->k_size-2]), 2)*pow((k_array[i])/(k_array[pnl->k_size-2]), ppm->n_s-4.); 
		}
		//fprintf(stdout, "%e %e\n", k_array[i], pk_l_new[i]);
	}
  //fprintf(stdout, "%e, %e\n", k_array[pnl->k_size-1], pk_l_new[pnl->k_size-1]);            
  //fprintf(stdout, "%e, %e\n", k_array[pnl->k_size], pk_l_new[pnl->k_size]);            
  //fprintf(stdout, "%e, %e\n", k_array[pnl->k_size+1], pk_l_new[pnl->k_size+1]);           
  
  //_E_+*1.8/(13.41*0.012)
  i=0;
  index_k=i;
  i++;
  index_y=i;
  i++;
  index_ddy=i;
  i++;
  index_num=i;

  class_alloc(array_for_sigma,
              k_size_new*index_num*sizeof(double),
              pnl->error_message);
            
  for (i=0;i<k_size_new;i++) {
    //k=exp(lnk_l[i]);
    k=k_array[i];
    if (i == (k_size_new-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
    x=k*R;
    if (x<0.01) {
		W = 1.-(pow(x, 2.)/10.);
    }
    else {
		W=3./x/x/x*(sin(x)-x*cos(x));
    }
    //if (R=1e3/pba->h) fprintf(stdout, "%e\n", W);
    array_for_sigma[i*index_num+index_k]=k;
    array_for_sigma[i*index_num+index_y]=k*k*pk_l_new[i]*W*W;
  }



  class_call(array_spline(array_for_sigma,
                          index_num,
                          k_size_new, //pnl->k_size,
                          index_k,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  class_call(array_integrate_all_spline(array_for_sigma,
                                        index_num,
                                        k_size_new, //pnl->k_size,
                                        index_k,
                                        index_y,
                                        index_ddy,
                                        sigma,
                                        pnl->error_message),
             pnl->error_message,
             pnl->error_message);
	for (i=0;i<k_size_new;i++) {
		//fprintf(stdout, "%e, %e, %e\n",array_for_sigma[i*index_num+index_k], array_for_sigma[i*index_num+index_y], array_for_sigma[i*index_num+index_ddy]);
	}	
	
	free(k_array);
	free(pk_l_new);
  free(array_for_sigma);

  *sigma = sqrt(*sigma/(2.*_PI_*_PI_));
	//fprintf(stdout, "%e\n", *sigma);
  return _SUCCESS_;  

}*/
/*
int nonlinear_sigma(
                  struct precision * ppr,
                  struct background * pba,
                  struct perturbs * ppt,
                  struct primordial * ppm,
                  struct nonlinear * pnl,
                  double R,
                  double *pk_l,                
                  double * sigma
                  ) {
  double pk;

  double * array_for_sigma;
  int index_num;
  int index_k;
  int index_y;
  int index_ddy;
  int i;

  double k,W,x;
  
  i=0;
  index_k=i;
  i++;
  index_y=i;
  i++;
  index_ddy=i;
  i++;
  index_num=i;

  class_alloc(array_for_sigma,
              pnl->k_size_extra*index_num*sizeof(double),
              pnl->error_message);
            
  for (i=0;i<pnl->k_size_extra;i++) {
    k=pnl->k_extra[i];
    //t = 1./(1.+k);
    //if (i>0) fprintf(stdout, "%e, %e, %e, %e\n", k, pk_l[i], t, pnl->k_extra[i]-pnl->k_extra[i-1]);
    if (i == (pnl->k_size_extra-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
    x=k*R;
	if (x<0.01) {
		W = 1.-(pow(x, 2.)/10.);
	}	
	else {
		 W = 3./x/x/x*(sin(x)-x*cos(x));
    }
    array_for_sigma[i*index_num+index_k] = k;
    array_for_sigma[i*index_num+index_y] = k*k*pk_l[i]*W*W;
  }

  class_call(array_spline(array_for_sigma,
                          index_num,
                          pnl->k_size_extra,
                          index_k,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  class_call(array_integrate_all_trapzd_or_spline(array_for_sigma,
                                        index_num,
                                        pnl->k_size_extra,
                                        pnl->k_size_extra-1,
                                        index_k,
                                        index_y,
                                        index_ddy,
                                        sigma,
                                        pnl->error_message),
             pnl->error_message,
             pnl->error_message);

//	for (i=0;i<pnl->k_size_extra;i++) {
//		fprintf(stdout, "%e, %e, %e, %e\n", *sigma, array_for_sigma[i*index_num+index_k], array_for_sigma[i*index_num+index_y], array_for_sigma[i*index_num+index_ddy]);    		
//	}

  free(array_for_sigma);

  *sigma = sqrt(*sigma/(2.*_PI_*_PI_));
  fprintf(stdout, "%e\n", *sigma); 

  return _SUCCESS_;

  

}*/


int nonlinear_sigma(
                  struct precision * ppr,
                  struct background * pba,
                  struct perturbs * ppt,
                  struct primordial * ppm,
                  struct nonlinear * pnl,
                  double R,
                  double *pk_l,                
                  double * sigma
                  ) {
  double pk;

  double * array_for_sigma;
  int index_num;
  int index_k;
  int index_y;
  int index_ddy;
  int i;

  double k,W,x,t;
  
  i=0;
  index_k=i;
  i++;
  index_y=i;
  i++;
  index_ddy=i;
  i++;
  index_num=i;

  class_alloc(array_for_sigma,
              pnl->k_size_extra*index_num*sizeof(double),
              pnl->error_message);
            
  for (i=pnl->k_size_extra-1;i>=0;i--) {
    k=pnl->k_extra[i];
    t = 1./(1.+k);
    //if (i>0) fprintf(stdout, "%e, %e, %e, %e\n", k, pk_l[i], t, pnl->k_extra[i]-pnl->k_extra[i-1]);
    if (i == (pnl->k_size_extra-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
    x=k*R;
	if (x<0.01) {
		W = 1.-(pow(x, 2.)/10.);
	}	
	else {
		 W = 3./x/x/x*(sin(x)-x*cos(x));
    }
    array_for_sigma[(pnl->k_size_extra-1-i)*index_num+index_k] = t;
    array_for_sigma[(pnl->k_size_extra-1-i)*index_num+index_y] = k*k*k*pk_l[i]*W*W/(t*(1.-t));
  }

  class_call(array_spline(array_for_sigma,
                          index_num,
                          pnl->k_size_extra,
                          index_k,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  class_call(array_integrate_all_trapzd_or_spline(array_for_sigma,
                                        index_num,
                                        pnl->k_size_extra,
                                        0, //pnl->k_size_extra-1,
                                        index_k,
                                        index_y,
                                        index_ddy,
                                        sigma,
                                        pnl->error_message),
             pnl->error_message,
             pnl->error_message);

//	for (i=0;i<pnl->k_size_extra;i++) {
//		fprintf(stdout, "%e, %e, %e, %e\n", *sigma, array_for_sigma[i*index_num+index_k], array_for_sigma[i*index_num+index_y], array_for_sigma[i*index_num+index_ddy]);    		
//	}

  free(array_for_sigma);

  *sigma = sqrt(*sigma/(2.*_PI_*_PI_));
  //fprintf(stdout, "%e\n", *sigma); 

  return _SUCCESS_;

  

}

/*
int nonlinear_sigma(
                  struct precision * ppr,
                  struct background * pba,
                  struct perturbs * ppt,
                  struct primordial * ppm,
                  struct nonlinear * pnl,
                  double R,
                  double *pk_l,                
                  double * sigma
                  ) {
  double pk;

  double * array_for_sigma;
  int index_num;
  int index_k;
  int index_y;
  int index_ddy;
  int i, j, n;

  double k,W,product,t;
  
  double a, b;
  double x, dx;
  double f1, f2, fx;
  double sum_n, sum_2n, sum_new, sum_old
  int jmin, jmax
  
  i=0;
  index_k=i;
  i++;
  index_y=i;
  i++;
  index_ddy=i;
  i++;
  index_num=i;

  class_alloc(array_for_sigma,
              pnl->k_size_extra*index_num*sizeof(double),
              pnl->error_message);
            
  for (i=0;i<pnl->k_size_extra;i++) {
    
    
    k=pnl->k_extra[i];
    t = 1./(1.+k);
    //if (i>0) fprintf(stdout, "%e, %e, %e, %e\n", k, pk_l[i], t, pnl->k_extra[i]-pnl->k_extra[i-1]);
    if (i == (pnl->k_size_extra-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
    product=k*R;
	if (product<0.01) {
		W = 1.-(pow(product, 2.)/10.);
	}	
	else {
		 W = 3./product/product/product*(sin(product)-product*cos(product));
    }
    array_for_sigma[i*index_num+index_k]=t;
    array_for_sigma[i*index_num+index_y]=pow(k, 3)*pk_l[i]*W*W/(t*(1.-t));
  }

  class_call(array_spline(array_for_sigma,
                          index_num,
                          pnl->k_size_extra,
                          index_k,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          pnl->error_message),
             pnl->error_message,
             pnl->error_message);

	jmin=5;
	jmax=30;
	a=0.;
	b=1.;
	
	if (a==b) {
		sigint0=0.;
	} 
	else {
		sum_2n=0.d0

       for (j=1;jmax;j++) {
          
          //Note, you need this to be 1+2**n for some integer n
          //j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+pow(2,(j-1));

          //Calculate the dx interval for this value of 'n'
          dx=(b-a)/(n-1.);

          if (j==1) {
             
             //The first go is just the trapezium of the end points
             f1=array_for_sigma[0];
             f2=array_for_sigma[pnl->k_size_extra];
             sum_2n=0.5*(f1+f2)*dx;
             sum_new=sum_2n;
             
          }
          else {

             //Loop over only new even points to add these to the integral
             for (i=2;n;i=i+2) {
                x=a+(b-a)*(i-1.)/(n-1.)
                fx=sigma_integrand_transformed(x,r,f0_rapid,z,cosm)
                sum_2n=sum_2n+fx
             }
             //END DO

             //Now create the total using the old and new parts
             sum_2n=sum_n/2.d0+sum_2n*dx

             //Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.d0*sum_2n-sum_n)/3.d0 !This is Simpson's rule and cancels error
             ELSE
                STOP 'SIGINT0: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.d0+sum_new/sum_old)<acc)) THEN
             //jmin avoids spurious early convergence
             sigint0=REAL(sum_new)
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'SIGINT0: Integration timed out'
          ELSE
             //Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.d0
          END IF

       }//END DO
	}
	// end if



	for (i=0;i<pnl->k_size_extra;i++) {
		fprintf(stdout, "%e, %e, %e, %e\n", *sigma, array_for_sigma[i*index_num+index_k], array_for_sigma[i*index_num+index_y], array_for_sigma[i*index_num+index_ddy]);    		
	}

  free(array_for_sigma);

  *sigma = sqrt(*sigma/(2.*_PI_*_PI_));
  fprintf(stdout, "%e\n", *sigma); 

  return _SUCCESS_;

  

}
*/
/*
int nonlinear_sigma_disp(
                  struct precision * ppr,
                  struct background * pba,
                  struct perturbs * ppt,
                  struct primordial * ppm,
                  struct nonlinear * pnl,
                  double R,
                  double *pk_l,
                  double *lnk_l,                 
                  double * sigma_disp
                  ) {
  
  double * array_for_sigma_disp;
  double * k_array;
  double * pk_l_new;
  double lnk_min, lnk_max;
  int k_size_new;
  int index_num;
  int index_k;
  int index_y;
  int index_ddy;
  int i;

  double t,k,W,x;
  
  // extrapolate k array up to k_max=1e4 h/Mpc  
  
  lnk_max = log(1.e4)*pba->h;
  lnk_min = log(ppt->k_max_for_pk);
  k_size_new = pnl->k_size+(lnk_max-lnk_min)*(ppr->k_per_decade_for_pk);//2.302585e-01);
  //fprintf(stdout, "%d, %d, %e, %e\n", pnl->k_size, k_size_new, ppt->k_max_for_pk, ppr->k_per_decade_for_pk);
  class_alloc(k_array,
              k_size_new*sizeof(double),
              pnl->error_message); 
  
  class_alloc(pk_l_new,
              k_size_new*sizeof(double),
              pnl->error_message); 
  
  for (i=0;i<k_size_new;i++){
		if (i<pnl->k_size){
			k_array[i] = exp(lnk_l[i]);
			pk_l_new[i] = pk_l[i];
		}else{
			k_array[i] = exp(lnk_min+(lnk_max-lnk_min)*(i-pnl->k_size)/(k_size_new-1.-pnl->k_size));
			pk_l_new[i] = pk_l[pnl->k_size-1]*pow(log(k_array[i])/log(k_array[pnl->k_size-1]), 2)*pow((k_array[i])/(k_array[pnl->k_size-1]), ppm->n_s-4.); 
		}
		//fprintf(stdout, "%e, %e\n", k_array[i], pk_l_new[i]);
	}
  //fprintf(stdout, "%e, %e\n", k_array[pnl->k_size-1], pk_l_new[pnl->k_size-1]);            
  //fprintf(stdout, "%e, %e\n", k_array[pnl->k_size], pk_l_new[pnl->k_size]);            
  //fprintf(stdout, "%e, %e\n", k_array[pnl->k_size+1], pk_l_new[pnl->k_size+1]);           
      
  
  i=0;
  index_k=i;
  i++;
  index_y=i;
  i++;
  index_ddy=i;
  i++;
  index_num=i;

  class_alloc(array_for_sigma_disp,
              k_size_new*index_num*sizeof(double),
              pnl->error_message);
            
  for (i=0;i<k_size_new;i++) {
    //k=exp(lnk_l[i]);
    k=k_array[i];
    if (i == (k_size_new-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
    x=k*R;
    if (x<0.01) {
		W = 1.-(pow(x, 2.)/10.);
    }
    else {
		W=3./x/x/x*(sin(x)-x*cos(x));
    }
    array_for_sigma_disp[i*index_num+index_k]=k;
    array_for_sigma_disp[i*index_num+index_y]=k*k*pk_l_new[i]*W*W;
  }



  class_call(array_spline(array_for_sigma_disp,
                          index_num,
                          k_size_new, //pnl->k_size,
                          index_k,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  class_call(array_integrate_all_spline(array_for_sigma_disp,
                                        index_num,
                                        k_size_new, //pnl->k_size,
                                        index_k,
                                        index_y,
                                        index_ddy,
                                        sigma_disp,
                                        pnl->error_message),
             pnl->error_message,
             pnl->error_message);
	for (i=0;i<k_size_new;i++) {
		//fprintf(stdout, "%e, %e, %e\n",array_for_sigma_disp[i*index_num+index_k], array_for_sigma_disp[i*index_num+index_y], array_for_sigma_disp[i*index_num+index_ddy]);
	}	
	
	free(k_array);
	free(pk_l_new);
  free(array_for_sigma_disp);

  *sigma_disp = sqrt(*sigma_disp/(2.*_PI_*_PI_)/3); // unit: [Mpc]
	//fprintf(stdout, "%e\n", *sigma_disp);
  return _SUCCESS_;  

}*/

int nonlinear_sigma_disp(
                  struct precision * ppr,
                  struct background * pba,
                  struct perturbs * ppt,
                  struct primordial * ppm,
                  struct nonlinear * pnl,
                  double R,
                  double *pk_l,               
                  double * sigma_disp
                  ) {
  double pk;

  double * array_for_sigma_disp;
  int index_num;
  int index_k;
  int index_y;
  int index_ddy;
  int i;

  double k,W,x,t;
  
  i=0;
  index_k=i;
  i++;
  index_y=i;
  i++;
  index_ddy=i;
  i++;
  index_num=i;

  class_alloc(array_for_sigma_disp,
              pnl->k_size_extra*index_num*sizeof(double),
              pnl->error_message);
            
  for (i=0;i<pnl->k_size_extra;i++) {
    k=pnl->k[i];
    if (i == (pnl->k_size_extra-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
    x=k*R;
	if (x<0.01) {
		W = 1.-(pow(x, 2.)/10.);
	}	
	else {
		 W = 3./x/x/x*(sin(x)-x*cos(x));
    }
    array_for_sigma_disp[i*index_num+index_k]=k;
    array_for_sigma_disp[i*index_num+index_y]=pk_l[i]*W*W;
  }
/*
for (i=pnl->k_size_extra-1;i>=0;i--) {
    k=pnl->k_extra[i];
    t = 1./(1.+k);
    //if (i>0) fprintf(stdout, "%e, %e, %e, %e\n", k, pk_l[i], t, pnl->k_extra[i]-pnl->k_extra[i-1]);
    if (i == (pnl->k_size_extra-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
    x=k*R;
	if (x<0.01) {
		W = 1.-(pow(x, 2.)/10.);
	}	
	else {
		 W = 3./x/x/x*(sin(x)-x*cos(x));
    }
    array_for_sigma_disp[(pnl->k_size_extra-1-i)*index_num+index_k] = t;
    array_for_sigma_disp[(pnl->k_size_extra-1-i)*index_num+index_y] = k*pk_l[i]*W*W/(t*(1.-t));
  }*/



  class_call(array_spline(array_for_sigma_disp,
                          index_num,
                          pnl->k_size,
                          index_k,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  class_call(array_integrate_all_spline(array_for_sigma_disp,
                                        index_num,
                                        pnl->k_size,
                                        index_k,
                                        index_y,
                                        index_ddy,
                                        sigma_disp,
                                        pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  free(array_for_sigma_disp);

  *sigma_disp = sqrt(*sigma_disp/(2.*_PI_*_PI_)/3); // unit: [Mpc]

  return _SUCCESS_;

  

}

int nonlinear_fill_sigtab(
              struct precision *ppr,
						  struct background * pba,
						  struct perturbs * ppt,
						  struct primordial * ppm,
						  struct nonlinear * pnl,
						  double *pk_l             
						  ) {
	
  //double rmin;
  //double rmax;
  //int nsig;
  double r;
  double rmin, rmax;
  double sig;
  int i, index_r, index_sig, index_ddsig, index_n, nsig;
	
	rmin = ppr->rmin_for_sigtab/pba->h;
  rmax = ppr->rmax_for_sigtab/pba->h;
  nsig = ppr->n_hmcode_tables;
  
  i=0;
  index_r=i;
  i++;
  index_sig=i;
  i++;
  index_ddsig=i;
  i++;
  index_n=i;
  
  class_alloc((pnl->sigtab),(nsig*index_n*sizeof(double)),pnl->error_message);
  class_alloc(pnl->rtab,nsig*sizeof(double),pnl->error_message);
  class_alloc(pnl->stab,nsig*sizeof(double),pnl->error_message);
  class_alloc(pnl->ddstab,nsig*sizeof(double),pnl->error_message);
  
  for (i=0;i<nsig;i++){
	r=exp(log(rmin)+log(rmax/rmin)*i/(nsig-1));
    class_call(nonlinear_sigma(ppr,pba,ppt,ppm,pnl,r,pk_l,&sig), 
      pnl->error_message, pnl->error_message); 
    pnl->sigtab[i*index_n+index_r]=r;
    pnl->sigtab[i*index_n+index_sig]=sig;
  }
  
  class_call(array_spline(pnl->sigtab,
						  index_n,
						  nsig,
						  index_r,
						  index_sig,
						  index_ddsig,
						  _SPLINE_EST_DERIV_,
						  pnl->error_message), 
					pnl->error_message,
					pnl->error_message); 
  for (i=0;i<nsig;i++){
	pnl->rtab[i] = pnl->sigtab[i*index_n+index_r];
	pnl->stab[i] = pnl->sigtab[i*index_n+index_sig];
	pnl->ddstab[i] = pnl->sigtab[i*index_n+index_ddsig];  
	//if (i==0) fprintf(stdout, "sigma max = %e\n", pnl->stab[i]);  
	//if (i==nsig-1) fprintf(stdout, "sigma min = %e\n", pnl->stab[i]);  

	//fprintf(stdout, "%e, %e, %e\n",pnl->sigtab[i*index_n+index_r], pnl->sigtab[i*index_n+index_sig], pnl->sigtab[i*index_n+index_ddsig]);
  }      
  return _SUCCESS_; 
}	

int nonlinear_fill_growtab(
              struct precision *ppr,
						  struct background * pba,
						  struct nonlinear * pnl            
						  ){
	
	double z, ainit, amax, scalefactor, tau_growth;
  int i, index_z, index_tau_growth, index_growth, index_gcol, last_index, ng;  
  double * pvecback;
  
  ng = ppr->n_hmcode_tables;
  ainit = ppr->ainit_for_growtab;
  amax = ppr->amax_for_growtab;
  
  i=0;
  index_z = i;
  i++;
  index_tau_growth = i;
  i++;
  index_growth = i;
  i++;
  index_gcol = i;
  
	last_index = 0;
	
  class_alloc(pnl->growtable,ng*sizeof(double),pnl->error_message);
  class_alloc(pnl->ztable,ng*sizeof(double),pnl->error_message);
  class_alloc(pnl->tautable,ng*sizeof(double),pnl->error_message);
  class_alloc(pvecback,pba->bg_size*sizeof(double),pnl->error_message);
  
  for (i=0;i<ng;i++){
		scalefactor = ainit+(amax-ainit)*(i)/(ng-1);
		z = 1./scalefactor-1.;
		
		pnl->ztable[i] = z;
		
		class_call(background_tau_of_z(
                        pba,
                        z,
                        &tau_growth
                        ),
					pnl->error_message, pnl->error_message);
					
		pnl->tautable[i] = tau_growth;		
			
		class_call(background_at_tau(pba,tau_growth,pba->long_info,pba->inter_normal,&last_index,pvecback),
             pba->error_message,
             pnl->error_message);
             
    pnl->growtable[i] = pvecback[pba->index_bg_D];   
    //fprintf(stdout, "%e, %e %e\n", z, tau_growth, pnl->growtable[i]);      			
	}								

	free(pvecback);

}

/* in original HMcode, this is equivalent to the function Ci(x) */
int nonlinear_ci(
								 double x,
								 double *Ci
								 ){
	
	double x2, y, f, g, ci8;
	double em_const = 0.577215664901532861e0;
	
	if (fabs(x)<=4.){

       x2=x*x;

       ci8=em_const+log(x)+x2*(-0.25e0+x2*(7.51851524438898291e-3+x2*(-1.27528342240267686e-4
            +x2*(1.05297363846239184e-6+x2*(-4.68889508144848019e-9+x2*(1.06480802891189243e-11
            +x2*(-9.93728488857585407e-15)))))))/ (1.+x2*(1.1592605689110735e-2+
            x2*(6.72126800814254432e-5+x2*(2.55533277086129636e-7+x2*(6.97071295760958946e-10+
            x2*(1.38536352772778619e-12+x2*(1.89106054713059759e-15+x2*(1.39759616731376855e-18))))))));

       *Ci=ci8;
	}
    else {

       y=1./(x*x); 

       f = (1.e0 + y*(7.44437068161936700618e2 + y*(1.96396372895146869801e5 + 
            y*(2.37750310125431834034e7 +y*(1.43073403821274636888e9 + y*(4.33736238870432522765e10
            + y*(6.40533830574022022911e11 + y*(4.20968180571076940208e12 + y*(1.00795182980368574617e13
            + y*(4.94816688199951963482e12 +y*(-4.94701168645415959931e11)))))))))))/
            (x*(1. +y*(7.46437068161927678031e2 +y*(1.97865247031583951450e5 +
            y*(2.41535670165126845144e7 + y*(1.47478952192985464958e9 + 
            y*(4.58595115847765779830e10 +y*(7.08501308149515401563e11 + y*(5.06084464593475076774e12 
            + y*(1.43468549171581016479e13 + y*(1.11535493509914254097e13)))))))))));   

       g = y*(1.e0 + y*(8.1359520115168615e2 + y*(2.35239181626478200e5 + y*(3.12557570795778731e7
            + y*(2.06297595146763354e9 + y*(6.83052205423625007e10 +
            y*(1.09049528450362786e12 + y*(7.57664583257834349e12 +
            y*(1.81004487464664575e13 + y*(6.43291613143049485e12 +y*(-1.36517137670871689e12)))))))))))
            / (1. + y*(8.19595201151451564e2 +y*(2.40036752835578777e5 +
            y*(3.26026661647090822e7 + y*(2.23355543278099360e9 + y*(7.87465017341829930e10 
            + y*(1.39866710696414565e12 + y*(1.17164723371736605e13 + y*(4.01839087307656620e13 +y*(3.99653257887490811e13))))))))));
       *Ci=f*sin(x)-g*cos(x);
	}
																
	return _SUCCESS_;
}

/* in original HMcode, this is equivalent to the function Si(x) */
int nonlinear_si(
								 double x,
								 double *Si
								 ){
	
	double x2, y, f, g, si8;
	double pi8=3.1415926535897932384626433;
	
	if (fabs(x)<=4.){

       x2=x*x;

       si8 = x*(1.e0+x2*(-4.54393409816329991e-2+x2*(1.15457225751016682e-3
            +x2*(-1.41018536821330254e-5+x2*(9.43280809438713025e-8+x2*(-3.53201978997168357e-10
            +x2*(7.08240282274875911e-13+x2*(-6.05338212010422477e-16))))))))/ 
            (1.+x2*(1.01162145739225565e-2 +x2*(4.99175116169755106e-5+
            x2*(1.55654986308745614e-7+x2*(3.28067571055789734e-10+x2*(4.5049097575386581e-13
            +x2*(3.21107051193712168e-16)))))));

       *Si=si8;
	}
    else {

       y=1./(x*x);

       f = (1.e0 + y*(7.44437068161936700618e2 + y*(1.96396372895146869801e5 +
            y*(2.37750310125431834034e7 +y*(1.43073403821274636888e9 + y*(4.33736238870432522765e10 
            + y*(6.40533830574022022911e11 + y*(4.20968180571076940208e12 + 
            y*(1.00795182980368574617e13 + y*(4.94816688199951963482e12 +
            y*(-4.94701168645415959931e11)))))))))))/ (x*(1. +y*(7.46437068161927678031e2 +
            y*(1.97865247031583951450e5 +y*(2.41535670165126845144e7 + 
            y*(1.47478952192985464958e9 + y*(4.58595115847765779830e10 +
            y*(7.08501308149515401563e11 + y*(5.06084464593475076774e12 + 
            y*(1.43468549171581016479e13 + y*(1.11535493509914254097e13)))))))))));


       g = y*(1.e0 + y*(8.1359520115168615e2 + y*(2.35239181626478200e5 + 
            y*(3.12557570795778731e7 + y*(2.06297595146763354e9 + y*(6.83052205423625007e10 +
            y*(1.09049528450362786e12 + y*(7.57664583257834349e12 +y*(1.81004487464664575e13 +
            y*(6.43291613143049485e12 +y*(-1.36517137670871689e12)))))))))))/
            (1. + y*(8.19595201151451564e2 +y*(2.40036752835578777e5 + y*(3.26026661647090822e7 
            + y*(2.23355543278099360e9 + y*(7.87465017341829930e10 + y*(1.39866710696414565e12 
            + y*(1.17164723371736605e13 + y*(4.01839087307656620e13 +y*(3.99653257887490811e13))))))))));

       *Si=pi8/2.-f*cos(x)-g*sin(x);

	}
																
	return _SUCCESS_;
}

/* in original HMcode, this is equivalent to the function winnfw(k,rv,c) */
int nonlinear_window_nfw(
						  					 struct nonlinear * pnl,
						  					 double k,
						  					 double rv,
						  					 double c,
												 double *window_nfw
												 ){	
	double si1, si2, ci1, ci2, ks;
	double p1, p2, p3;									 
																
	ks = k*rv/c;
	
	class_call(nonlinear_si(
													ks*(1.+c),
													&si2
													),
					pnl->error_message, pnl->error_message);	
	
	class_call(nonlinear_si(
													ks,
													&si1
													),
					pnl->error_message, pnl->error_message);				
	
	class_call(nonlinear_ci(
													ks*(1.+c),
													&ci2
													),
					pnl->error_message, pnl->error_message);	
	
	class_call(nonlinear_ci(
													ks,
													&ci1
													),
					pnl->error_message, pnl->error_message);		

  p1=cos(ks)*(ci2-ci1);
  p2=sin(ks)*(si2-si1);
  p3=sin(ks*c)/(ks*(1.+c));

  *window_nfw=p1+p2-p3;
  *window_nfw=*window_nfw/(log(1.+c)-c/(1.+c));
	
	return _SUCCESS_;
}

  




int nonlinear_HMcode(
                      struct precision *ppr,
                      struct background *pba,
                      struct perturbs *ppt,
                      struct primordial *ppm,
                      struct nonlinear *pnl,
                      double tau,
                      double *pk_l,
                      double *pk_nl,
                      double *lnk_l,
                      double *lnpk_l,
                      double *ddlnpk_l,
                      double *k_nl
                      ) {

  double Omega_m,Omega_v,fnu,Omega0_m, w0, dw_over_da_fld, integral_fld;	
  double z_at_tau;
  double mmin, mmax, m, r, nu, nu_nl, z_form, g_form;
  int n, i, ng, nsig;
  int index_k, index_lnsig, index_dlnsig, index_ncol;
  int last_index;  
  double pk_lin, pk_2h, pk_1h;
  double sigma_disp, sigma_disp100, sigma8, sig, sigmatest, sigf, r_nl, ln_r_nl, n_eff, radius;
  double * rtab;
  double * sigmatab;
  double * ddy_sigtab;
  double P;
  double * pvecback;
  double growth;
  
  double anorm;
  
  double * conc;
  double * mass;
  double * sigma_r;
	double * sigmaf_r;
	double * ln_sigma_squared;
	double * ln_dsigma_dr;
  double * r_virial;
  double * r_real;
  double * ln_r_real;
  double * nu_arr;
  
  double delta_c, Delta_v;
  double fraction;
  double Abary, eta, eta0;
  double fdamp;
  double alpha;
  double dlnsigdlnR;
  
  double window_nfw;
  
  /* include precision parameters that control the number of entries in the growth and sigma tables */
  ng = ppr->n_hmcode_tables;
  nsig = ppr->n_hmcode_tables;
  
  /* Call all the relevant background parameters at this tau */
  class_alloc(pvecback,pba->bg_size*sizeof(double),pnl->error_message);
	last_index=0;
  Omega0_m = (pba->Omega0_cdm + pba->Omega0_b + pba->Omega0_ncdm_tot + pba->Omega0_dcdm);
  class_call(background_w_fld(pba,pba->a_today,&w0,&dw_over_da_fld,&integral_fld), pba->error_message, pnl->error_message);

  fnu      = pba->Omega0_ncdm_tot/Omega0_m;
  anorm    = 1./(2*pow(_PI_,2));
  
  class_call(background_at_tau(pba,tau,pba->long_info,pba->inter_normal,&last_index,pvecback),
             pba->error_message,
             pnl->error_message);
	
  Omega_m = pvecback[pba->index_bg_Omega_m];
  Omega_v = 1.-pvecback[pba->index_bg_Omega_m]-pvecback[pba->index_bg_Omega_r];
  growth = pvecback[pba->index_bg_D];
  z_at_tau = 1./pvecback[pba->index_bg_a]-1.;
 
  free(pvecback);
    
  /** Get sigma(R=8), sigma_disp(R=0), sigma_disp(R=100)  */

  class_call(nonlinear_sigma(ppr,pba,ppt,ppm,pnl,8./pba->h,pk_l,&sigma8), 
			pnl->error_message, pnl->error_message);
	class_call(nonlinear_sigma(ppr,pba,ppt,ppm,pnl,1.e-4/pba->h,pk_l,&sigmatest), 
			pnl->error_message, pnl->error_message);
			//fprintf(stdout, "%e", sigmatest);		
  class_call(nonlinear_sigma_disp(ppr,pba,ppt,ppm,pnl,0.,pk_l,&sigma_disp), 
			pnl->error_message, pnl->error_message);
  class_call(nonlinear_sigma_disp(ppr,pba,ppt,ppm,pnl,100./pba->h,pk_l,&sigma_disp100), 
			pnl->error_message, pnl->error_message);			
  
   /** Initialisation steps for the 1-Halo Power Integral */
  mmin=ppr->mmin_for_p1h_integral/pba->h; //Minimum mass for integration; (unit conversion from  m[Msun/h] to m[Msun]  )
  mmax=ppr->mmax_for_p1h_integral/pba->h; //Maximum mass for integration;
  n=ppr->nsteps_for_p1h_integral; //Number of points for integration;

  class_alloc(mass,n*sizeof(double),pnl->error_message);
  class_alloc(r_real,n*sizeof(double),pnl->error_message);
  class_alloc(r_virial,n*sizeof(double),pnl->error_message); 
  class_alloc(sigma_r,n*sizeof(double),pnl->error_message);
  class_alloc(sigmaf_r,n*sizeof(double),pnl->error_message);
  class_alloc(nu_arr,n*sizeof(double),pnl->error_message);
  
  /* Linear theory density perturbation threshold for spherical collapse */
  // 
  delta_c=1.59+0.0314*log(sigma8);
  delta_c=delta_c*(1.+0.0123*log10(Omega_m));
  Delta_v=418.*pow(Omega_m, -0.352);
  fraction = pow(0.01, 1./3.);
	
	/* Fill the arrays for the P1H Integral */
  for (i=0;i<n;i++){
	  m = exp(log(mmin)+log(mmax/mmin)*(i)/(n-1));
	  r = pow((3.*m/(4.*_PI_*(2.775e11*pow(pba->h, 2))*Omega_m)), (1./3.)); // 2.775e11=rho_c/h^2 [Msun/Mpc^3]
	  mass[i] = m;
	  r_real[i] = r;
	  r_virial[i] = r_real[i]/pow(Delta_v, 1./3.);
	  class_call(array_interpolate_spline(pnl->rtab,
										  nsig,
										  pnl->stab,
										  pnl->ddstab,
										  1,
										  r,
										  &last_index,
										  &sig,
										  1,
										  pnl->error_message),
					pnl->error_message, pnl->error_message);
	  class_call(array_interpolate_spline(pnl->rtab,
										  nsig,
										  pnl->stab,
										  pnl->ddstab,
										  1,
										  r*fraction,
										  &last_index,
										  &sigf,
										  1,
										  pnl->error_message),
					pnl->error_message, pnl->error_message);
	  sig=sig*growth;
		sigf=sigf*growth;
		nu=delta_c/sig;
	  sigma_r[i] = sig;
	  sigmaf_r[i] = sigf;
	  nu_arr[i] = nu;	
	  //if (tau==pba->conformal_age) fprintf(stdout, "%e, ", r_real[i]);				
	  //if (tau==pba->conformal_age) fprintf(stdout, "%e\n", r*fraction);				
  }
  /* //this is the computational exhaustive way
  for (i=0;i<n;i++){	  
	  class_call(nonlinear_sigma(ppr,pba,ppt,ppm,pnl,r/pba->h,pk_l,lnk_l,&sig), 
			pnl->error_message, pnl->error_message); //too expensive computationally, lut needed
	  class_call(nonlinear_sigma(ppr,pba,ppt,ppm,pnl,fraction*r/pba->h,pk_l,lnk_l,&sigf), 
			pnl->error_message, pnl->error_message); //too expensive computationally, lut needed
	  nu=delta_c/sig;
	  class_call(array_interpolate(pnl->sigtab,2,nsig,i,r_real[i],last_index,&sig,n,pnl->error_message),
		pnl->error_message, pnl->error_message)
	  sigma_r[i] = sig;
	  nu_arr[i] = nu;
	  
      if (tau > 14164){
	    fprintf(stdout, "%e\n", nu_arr[i]);
	    fprintf(stdout, "%e\n", sigma_disp);
	    fprintf(stdout, "%e\n", sigma_disp100);	
      }
 	  
  }*/
  
  // find nonlinear scales k_nl and r_nl:
  nu_nl = 1.;
  class_call(array_interpolate_two_arrays_one_column(
																			nu_arr,
																			r_real,
																			1,
																			0,
																			n,
																			nu_nl,
																			&r_nl,
																			pnl->error_message),
					pnl->error_message, pnl->error_message);
	*k_nl = 1./r_nl;
	ln_r_nl = log(r_nl);
  
  //Calculate neff and alpha:
  i=0;
  index_lnsig = i;
  i++;
  index_dlnsig = i;
  i++;
  index_ncol = i;
    
  class_alloc(ln_sigma_squared,n*index_ncol*sizeof(double),pnl->error_message);
  class_alloc(ln_r_real,n*sizeof(double),pnl->error_message);
  class_alloc(ln_dsigma_dr,n*sizeof(double),pnl->error_message);
  
  for (i=0;i<n;i++){
		ln_sigma_squared[i*index_ncol+index_lnsig] = log(pow(sigma_r[i], 2.));
		ln_r_real[i] = log(r_real[i]);
	} 
  
  class_call(array_derive1_order2_table_line_to_line(
				       ln_r_real,
				       n,
				       ln_sigma_squared,
				       index_ncol,
				       index_lnsig,
				       index_dlnsig,
				       pnl->error_message),
				       pnl->error_message, pnl->error_message);
	for (i=0;i<n;i++){
		ln_dsigma_dr[i] = ln_sigma_squared[i*index_ncol+index_dlnsig];
	} 

	class_call(array_interpolate_two_arrays_one_column(
										nu_arr,
										ln_dsigma_dr,
										1,
										0,
										n,
										nu_nl,
										&dlnsigdlnR,
										pnl->error_message),
					pnl->error_message, pnl->error_message);						       
/*	for (i=0;i<n;i++){
		if (tau==pba->conformal_age) fprintf(stdout, "%e, %e, %e\n",r_real[i], sigma_r[i], ln_sigma_squared[i*index_ncol+index_dlnsig]);
	}
  fprintf(stdout, "%e, %e\n",log(r_nl), dlnsigdlnR);*/
  n_eff = -3.- dlnsigdlnR;
  alpha = 3.24*pow(1.85, n_eff);
  free(ln_sigma_squared);
  free(ln_dsigma_dr);
  
  // Calculate eta:
	Abary = 3.13;
	//eta0 = 0.98-0.12*Abary;
  eta0 = 0.603;
  eta = eta0 - 0.3*sigma8;
  
  //Calculate concentration (conc): 
	class_alloc(conc,n*sizeof(double),pnl->error_message);
  /*
  if (tau==pba->conformal_age) fprintf(stdout, "tau = %f; z = %f\n", tau, z_at_tau);
  if (tau==pba->conformal_age) fprintf(stdout, "--------------------------------------------------------------------\n");
  if (tau==pba->conformal_age) fprintf(stdout, "M,		rf,		s(fM),		zf,		gf,		c \n");
  if (tau==pba->conformal_age) fprintf(stdout, "--------------------------------------------------------------------\n");
*/
  for (i=0;i<n;i++){
		g_form = delta_c*growth/sigmaf_r[i];
		if (g_form > 1.) g_form = 1.;
		class_call(array_interpolate_two_arrays_one_column(
																pnl->growtable,
																pnl->ztable,
																1,
																0,
																ng,
																g_form,
																&z_form,
																pnl->error_message),
					pnl->error_message, pnl->error_message);	
		if (z_form < z_at_tau){
			conc[i] = Abary;
		} else {
			conc[i] = Abary*(1.+z_form)/(1.+z_at_tau);
	  } 	
		//if (tau==pba->conformal_age) fprintf(stdout, "%e, %e, %e, %e, %e, %e\n",mass[i], r_real[i]*fraction*pba->h, sigmaf_r[i], z_form, g_form, conc[i]);
	}
  //if (tau==pba->conformal_age) fprintf(stdout, "#################################################################\n");
  
  /*//prints out k versus power spectrum
  for (i=0;i<pnl->k_size;i++){
		if (tau==pba->conformal_age) fprintf(stdout, "%e, %e\n",exp(lnk_l[i])/pba->h, anorm*pk_l[i]*exp(3.*lnk_l[i]));
	}*/
  
  /* //prints out r-sigma table
	for (i=0;i<nsig*3;i=i+3){
		if (tau==pba->conformal_age) fprintf(stdout, "%e %e\n",pnl->rtab[i], pnl->sigtab[i+1]);
	}*/  
	
	/** Now compute the nonlinear correction */
	
	double  g, fac, et, ks, wk;
	double  * integrand; 
	double gst, p, q, A; //Parameters relevant for Sheth-Tormen Mass function
	
	
	ks=0.584*pow(sigma_disp, -1.);   //Damping wavenumber of the 1-halo term at very large scales;
	fdamp = 0.0095*pow(sigma_disp100*pba->h, 1.37); //Damping factor for 2-halo term 
	if (fdamp<1.e-3) fdamp=1.e-3;
  if (fdamp>0.99)  fdamp=0.99;
	
	p=0.3;
	q=0.707;
	A=0.21616;
	
	int index_nu;
  int index_y;
  int index_ddy;
  
  i=0;
  index_nu=i;
  i++;
  index_y=i;
  i++;
  index_ddy=i;
  i++;
  index_ncol=i;
	
	for (index_k = 0; index_k < pnl->k_size; index_k++){
		
		class_alloc(integrand,n*index_ncol*sizeof(double),pnl->error_message);
		
		pk_lin = pk_l[index_k]*pow(pnl->k[index_k],3)*anorm; //convert P_k in Delta_k^2
		
		for (i=0; i<n; i++){ //Calculates the integrand for the ph1 integral at all nu values
			class_call(nonlinear_window_nfw(
																	pnl,
																	pow(nu_arr[i], eta)*pnl->k[index_k],
																	r_virial[i],
																	conc[i],
																  &window_nfw),
					pnl->error_message, pnl->error_message);	
			
			/* This is the function gst(nu) in HMcode*/
			gst=A*(1.+(pow(q*nu_arr[i]*nu_arr[i], -p)))*exp(-q*nu_arr[i]*nu_arr[i]/2.);
			
			integrand[i*index_ncol+index_nu] = nu_arr[i];
			
			integrand[i*index_ncol+index_y] = mass[i]*gst*pow(window_nfw, 2.);
			
		} 
		class_call(array_spline(integrand,
                            index_ncol,
                            n,
                            index_nu,
                            index_y,
                            index_ddy,
                            _SPLINE_EST_DERIV_,
                            pnl->error_message),
             pnl->error_message,
             pnl->error_message);
		
		class_call(array_integrate_all_trapzd_or_spline(
																				integrand,
                                        index_ncol,
                                        n, 
                                        n-1, //0 or n-1
                                        index_nu,
                                        index_y,
                                        index_ddy,
                                        &pk_1h,
                                        pnl->error_message),
             pnl->error_message,
             pnl->error_message);
    
    if (pow(pnl->k[index_k]/ks, 2)>7.) fac = 0.;
    else fac=exp(-pow((pnl->k[index_k]/ks), 2.));
    
    pk_1h = pk_1h*anorm*pow(pnl->k[index_k],3)*(1.-fac)/(2.775e11*pow(pba->h, 2)*Omega_m);  // dimensionless power
		
		if (fdamp==0){
			pk_2h=pk_lin;
		}else{	
			pk_2h=pk_lin*(1.-fdamp*pow(tanh(pnl->k[index_k]*sigma_disp/sqrt(fdamp)), 2.)); //dimensionless power
		}
		if (pk_2h<0.) pk_2h=0.;	
		pk_nl[index_k] = pow((pow(pk_1h, alpha) + pow(pk_2h, alpha)), (1./alpha))/pow(pnl->k[index_k],3)/anorm; //converted back to P_k
		//if (tau==pba->conformal_age) fprintf(stdout, "%e %e %e %e %e\n", pnl->k[index_k], pk_lin, pk_1h, pk_2h, pk_nl[index_k]*pow(pnl->k[index_k],3)*anorm);			
		//if (tau==pba->conformal_age) fprintf(stdout, "%e, %e\n",pnl->k[index_k]/pba->h, pk_nl[index_k]);
		
		free(integrand);
	}
	
	  // print values for comparison with Fortran code 
  
  if (pnl->nonlinear_verbose > 0 && tau==pba->conformal_age){
		fprintf(stdout, " -> Parameters at redshift z = 0:\n"); 
		fprintf(stdout, "    sigd [Mpc/h]:	%e\n", sigma_disp*pba->h);
		fprintf(stdout, "    sigd100 [Mpc/h]:    %e\n", sigma_disp100*pba->h);
		fprintf(stdout, "    sigma8:		%e\n", sigma8);
		fprintf(stdout, "    nu min:		%e\n", nu_arr[0]);
		fprintf(stdout, "    nu max:		%e\n", nu_arr[n-1]);
		fprintf(stdout, "    r_v min [Mpc/h]:    %e\n", r_virial[0]*pba->h);
		fprintf(stdout, "    r_v max [Mpc/h]:    %e\n", r_virial[n-1]*pba->h);				
		fprintf(stdout, "    r_nl [Mpc/h]:	%e\n", r_nl*pba->h);
		fprintf(stdout, "    k_nl [h/Mpc]:	%e\n", *k_nl/pba->h);
		fprintf(stdout, "    neff:		%e\n", n_eff);
		fprintf(stdout, "    c min:		%e\n", conc[n-1]);
		fprintf(stdout, "    c max:		%e\n", conc[0]);			
		fprintf(stdout, "    Dv:			%e\n", Delta_v);
		fprintf(stdout, "    dc:			%e\n", delta_c);	
		fprintf(stdout, "    eta:		%e\n", eta);	
		fprintf(stdout, "    k*:			%e\n", ks);
		fprintf(stdout, "    Abary:		%e\n", Abary);			
		fprintf(stdout, "    fdamp:		%e\n", fdamp);		
		fprintf(stdout, "    alpha:		%e\n", alpha);		
  } 
	
  free(conc);	
  free(mass);
  free(r_real);
  free(r_virial);
  free(sigma_r);
  free(sigmaf_r);
  free(nu_arr);

  return _SUCCESS_;
}
  





