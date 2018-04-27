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
                        double * k_nl,
                        double * k_nl_cb
                        ) {

  double tau;
  int index_pk;

  class_call(background_tau_of_z(pba,
                                 z,
                                 &tau),
             pba->error_message,
             pnl->error_message);

  if (pnl->tau_size == 1) {
    *k_nl = pnl->k_nl[pnl->index_pk_m][0];
  }
  else {
    class_call(array_interpolate_two(pnl->tau,
                                     1,
                                     0,
                                     pnl->k_nl[pnl->index_pk_m],
                                     1,
                                     pnl->tau_size,
                                     tau,
                                     k_nl,
                                     1,
                                     pnl->error_message),
               pnl->error_message,
               pnl->error_message);
  }



 if (pba->has_ncdm){
 
 if (pnl->tau_size == 1) {
    *k_nl_cb = pnl->k_nl[pnl->index_pk_cb][0];
  }
  else {
    class_call(array_interpolate_two(pnl->tau,
                                     1,
                                     0,
                                     pnl->k_nl[pnl->index_pk_cb],
                                     1,
                                     pnl->tau_size,
                                     tau,
                                     k_nl_cb,
                                     1,
                                     pnl->error_message),
               pnl->error_message,
               pnl->error_message);
  }

  }
  else{
    *k_nl_cb = 1.e30;
  }

  return _SUCCESS_;
}


int nonlinear_hmcode_sigmaR_at_z(
                        struct precision *ppr,
                        struct background *pba,
                        struct nonlinear * pnl,
                        double R,
                        double z,
                        double * sigma_R
                        ) {
  
  int index_r, index_tau;
  double tau, sigma;
  double *sigma_array;
  
  class_alloc(sigma_array,ppr->n_hmcode_tables*sizeof(double),pnl->error_message);

  class_call(background_tau_of_z(pba,
                                 z,
                                 &tau),
             pba->error_message,
             pnl->error_message);
  
  for (index_r=0;index_r<ppr->n_hmcode_tables;index_r++){
    
    class_call(array_interpolate_two(pnl->tau,
                                     1,
                                     0,
                                     pnl->stab+index_r*pnl->tau_size,
                                     1,
                                     pnl->tau_size,
                                     tau,
                                     &sigma,
                                     1,
                                     pnl->error_message),
               pnl->error_message,
               pnl->error_message);
    
    sigma_array[index_r] = sigma;
  }            
  
  class_call(array_interpolate_two(pnl->rtab,
                                     1,
                                     0,
                                     sigma_array,
                                     1,
                                     ppr->n_hmcode_tables,
                                     R,
                                     sigma_R,
                                     1,
                                     pnl->error_message),
               pnl->error_message,
               pnl->error_message);
 /*
  for (index_r=0;index_r<ppr->n_hmcode_tables;index_r++){
    fprintf(stdout, "%e %e %e\n", pnl->rtab[index_r], sigma_array[index_r], pnl->stab[index_r*(pnl->tau_size)+pnl->tau_size-1]);
  }
*/
  free(sigma_array);
  
  return _SUCCESS_;  
}

int nonlinear_hmcode_sigma8_at_z(
                        struct background *pba,
                        struct nonlinear * pnl,
                        double z,
                        double * sigma_8
                        ) {

  double tau;

  class_call(background_tau_of_z(pba,
                                 z,
                                 &tau),
             pba->error_message,
             pnl->error_message);
             
  if (pnl->tau_size == 1) {
    *sigma_8 = pnl->sigma_8[0];
  } 
  else {
    class_call(array_interpolate_two(pnl->tau,
                                     1,
                                     0,
                                     pnl->sigma_8,
                                     1,
                                     pnl->tau_size,
                                     tau,
                                     sigma_8,
                                     1,
                                     pnl->error_message),
               pnl->error_message,
               pnl->error_message);
  }
  
  
  return _SUCCESS_;  
}

int nonlinear_hmcode_sigmadisp_at_z(
                        struct background *pba,
                        struct nonlinear * pnl,
                        double z,
                        double * sigma_disp
                        ) {

  double tau;

  class_call(background_tau_of_z(pba,
                                 z,
                                 &tau),
             pba->error_message,
             pnl->error_message);
             
  if (pnl->tau_size == 1) {
    *sigma_disp = pnl->sigma_disp[0];
  } 
  else {
    class_call(array_interpolate_two(pnl->tau,
                                     1,
                                     0,
                                     pnl->sigma_disp,
                                     1,
                                     pnl->tau_size,
                                     tau,
                                     sigma_disp,
                                     1,
                                     pnl->error_message),
               pnl->error_message,
               pnl->error_message);
  }
  
  
  return _SUCCESS_;  
}

int nonlinear_hmcode_sigmadisp100_at_z(
                        struct background *pba,
                        struct nonlinear * pnl,
                        double z,
                        double * sigma_disp_100
                        ) {

  double tau;

  class_call(background_tau_of_z(pba,
                                 z,
                                 &tau),
             pba->error_message,
             pnl->error_message);
             
  if (pnl->tau_size == 1) {
    *sigma_disp_100 = pnl->sigma_disp_100[0];
  } 
  else {
    class_call(array_interpolate_two(pnl->tau,
                                     1,
                                     0,
                                     pnl->sigma_disp_100,
                                     1,
                                     pnl->tau_size,
                                     tau,
                                     sigma_disp_100,
                                     1,
                                     pnl->error_message),
               pnl->error_message,
               pnl->error_message);
  }
  
  
  return _SUCCESS_;  
}

int nonlinear_hmcode_sigmaprime_at_z(
                        struct background *pba,
                        struct nonlinear * pnl,
                        double z,
                        double * sigma_prime
                        ) {

  double tau;

  class_call(background_tau_of_z(pba,
                                 z,
                                 &tau),
             pba->error_message,
             pnl->error_message);
             
  if (pnl->tau_size == 1) {
    *sigma_prime = pnl->sigma_prime[0];
  } 
  else {
    class_call(array_interpolate_two(pnl->tau,
                                     1,
                                     0,
                                     pnl->sigma_prime,
                                     1,
                                     pnl->tau_size,
                                     tau,
                                     sigma_prime,
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

  int index_pk;
  double **pk_l;
  double **pk_nl;
  double **lnk_l;
  double **lnpk_l;
  double **ddlnpk_l;

  short print_warning=_FALSE_;
  double * pvecback;
  int last_index;
  double a,z;
  short halofit_found_k_max;
  int pk_type;

  /** Summary
   *
   * (a) First deal with the case where non non-linear corrections requested */

  if (pnl->method == nl_none) {
    if (pnl->nonlinear_verbose > 0)
      printf("No non-linear spectra requested. Nonlinear module skipped.\n");
  }

  /** (b) Compute for HALOFIT or HMcode non-linear spectrum */

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

    index_pk = 0;
    class_define_index(pnl->index_pk_m,  _TRUE_, index_pk,1);
    class_define_index(pnl->index_pk_cb,  pba->has_ncdm, index_pk,1);
    pnl->pk_size = index_pk;
    //printf("pk_size=%d, index_pk_m=%d, index_pk_cb=%d\n",pnl->pk_size,pnl->index_pk_m,pnl->index_pk_cb);

    /** - copy list of (k,tau) from perturbation module */

    pnl->k_size = ppt->k_size[ppt->index_md_scalars];
    class_alloc(pnl->k,pnl->k_size*sizeof(double),pnl->error_message);
    for (index_k=0; index_k<pnl->k_size; index_k++)
      pnl->k[index_k] = ppt->k[ppt->index_md_scalars][index_k];

    pnl->tau_size = ppt->tau_size;
    class_alloc(pnl->tau,pnl->tau_size*sizeof(double),pnl->error_message);
    for (index_tau=0; index_tau<pnl->tau_size; index_tau++)
      pnl->tau[index_tau] = ppt->tau_sampling[index_tau];

    if (pnl->method == nl_HMcode){
      class_alloc(pnl->sigma_8,pnl->tau_size*sizeof(double),pnl->error_message);
      class_alloc(pnl->sigma_disp,pnl->tau_size*sizeof(double),pnl->error_message);    
      class_alloc(pnl->sigma_disp_100,pnl->tau_size*sizeof(double),pnl->error_message);    
      class_alloc(pnl->sigma_prime,pnl->tau_size*sizeof(double),pnl->error_message);
    }

    class_alloc(pnl->nl_corr_density,
                pnl->pk_size*sizeof(double *),
                pnl->error_message);

    class_alloc(pnl->k_nl,
                pnl->pk_size*sizeof(double *),
                pnl->error_message);

    class_alloc(pk_l,
                pnl->pk_size*sizeof(double *),
                pnl->error_message);

    class_alloc(pk_nl,
                pnl->pk_size*sizeof(double *),
                pnl->error_message);

    class_alloc(lnk_l,
                pnl->pk_size*sizeof(double *),
                pnl->error_message);

    class_alloc(lnpk_l,
                pnl->pk_size*sizeof(double *),
                pnl->error_message);

    class_alloc(ddlnpk_l,
                pnl->pk_size*sizeof(double *),
                pnl->error_message);

    for (index_pk=0; index_pk<pnl->pk_size; index_pk++){
    
    class_alloc(pnl->nl_corr_density[index_pk],pnl->tau_size*pnl->k_size*sizeof(double),pnl->error_message);
    class_alloc(pnl->k_nl[index_pk],pnl->tau_size*sizeof(double),pnl->error_message);

    class_alloc(pk_l[index_pk],pnl->k_size*sizeof(double),pnl->error_message);
    class_alloc(pk_nl[index_pk],pnl->k_size*sizeof(double),pnl->error_message);

    class_alloc(lnk_l[index_pk],pnl->k_size*sizeof(double),pnl->error_message);//this is not really necessary
    class_alloc(lnpk_l[index_pk],pnl->k_size*sizeof(double),pnl->error_message);
    class_alloc(ddlnpk_l[index_pk],pnl->k_size*sizeof(double),pnl->error_message);
  
    }

    if (pnl->method == nl_HMcode){
    /** initialise the source extrapolation */
      class_call(get_extrapolated_source_size(ppr->k_per_decade_for_pk,
                                                pnl->k[pnl->k_size-1], 
                                                ppr->hmcode_max_k_extra,
                                                pnl->k_size,
                                                &size_extrapolated_source,
                                                pnl->error_message),
          pnl->error_message,
          pnl->error_message);
      pnl->k_size_extra = size_extrapolated_source;
   
      class_alloc(pnl->k_extra,pnl->k_size_extra*sizeof(double),pnl->error_message);
      for (index_pk=0; index_pk<pnl->pk_size; index_pk++) {
        class_realloc(pk_l[index_pk],pk_l[index_pk],pnl->k_size_extra*sizeof(double),pnl->error_message);   
        class_realloc(lnk_l[index_pk],lnk_l[index_pk],pnl->k_size_extra*sizeof(double),pnl->error_message);
        class_realloc(lnpk_l[index_pk],lnpk_l[index_pk],pnl->k_size_extra*sizeof(double),pnl->error_message);
        class_realloc(ddlnpk_l[index_pk],ddlnpk_l[index_pk],pnl->k_size_extra*sizeof(double),pnl->error_message);
      }
      class_call(extrapolate_k(
						 pnl->k,
						 pnl->k_size,
					         pnl->k_extra,
						 ppr->k_per_decade_for_pk,
						 ppr->hmcode_max_k_extra,
					         pnl->error_message),
				   pnl->error_message,
				   pnl->error_message);
      
      /** fill table with scale independent growth factor */
	    class_call(nonlinear_hmcode_fill_growtab(ppr,pba,pnl), 
				pnl->error_message, pnl->error_message);	
       
     /** Set the baryonic feedback parameters according to the chosen feedback models */
      if (pnl->feedback == emu_dmonly){
        pnl->eta_0 = 0.603;
        pnl->c_min = 3.13;
      }
      if (pnl->feedback == owls_dmonly){
        pnl->eta_0 = 0.64;
        pnl->c_min = 3.43;
      }			
      if (pnl->feedback == owls_ref){
        pnl->eta_0 = 0.68;
        pnl->c_min = 3.91;
      }		
      if (pnl->feedback == owls_agn){
        pnl->eta_0 = 0.76;
        pnl->c_min = 2.32;
      }		
      if (pnl->feedback == owls_dblim){
        pnl->eta_0 = 0.70;
        pnl->c_min = 3.01;
      }
    }
    
/* //This is not necessary anymore, since the order of index_pk does not play a role  
  for (pk_type=pnl->pk_size-1; pk_type>=0; pk_type--) {

      if(pk_type == pnl->index_pk_m){
        index_pk = pnl->index_pk_m;
      }
      else if((pba->has_ncdm)&&(pk_type == pnl->index_pk_cb)){
        index_pk = pnl->index_pk_cb;
      }
      else {
        class_stop(pnl->error_message,"P(k) is set neither to total matter nor to cold dark matter + baryons, pk_type=%d \n",pk_type);
      } 
*/
    
    pnl->index_tau_min_nl = 0;  
    
    /** - loop over time */

    for (index_tau = pnl->tau_size-1; index_tau>=0; index_tau--) {
      
      //clock_t begin = clock();
      
      for (pk_type=0; pk_type<pnl->pk_size; pk_type++) {
        
        if(pk_type == 0) {
          if(pba->has_ncdm) {
            index_pk=pnl->index_pk_cb;
          }
          else {
            index_pk = pnl->index_pk_m;
          }
        }
        else if(pk_type == 1) {
          if(pba->has_ncdm){
            index_pk = pnl->index_pk_m;
          }
          else {
            class_stop(pnl->error_message,"looks like pk_size=2 even if you do not have any massive neutrinos");
          }
        }
        else {
         class_stop(pnl->error_message,"P(k) is set neither to total matter nor to cold dark matter + baryons, pk_type=%d \n",pk_type);
        }
        
      /* get P_L(k) at this time */
      class_call(nonlinear_pk_l(pba,ppt,ppm,pnl,index_pk,index_tau,pk_l[index_pk],lnk_l[index_pk],lnpk_l[index_pk],ddlnpk_l[index_pk]),
                 pnl->error_message,
                 pnl->error_message);

      // get P_NL(k) at this time with Halofit 
      if (pnl->method == nl_halofit) {
        if (print_warning == _FALSE_) {
    
            class_call(nonlinear_halofit(
                 ppr,
                 pba,
                 ppt,
                 ppm,
                 pnl,
                 index_pk,
                 pnl->tau[index_tau],
                 pk_l[index_pk],
                 pk_nl[index_pk],
                 lnk_l[index_pk],
                 lnpk_l[index_pk],
                 ddlnpk_l[index_pk],
                 &(pnl->k_nl[index_pk][index_tau]),																		      
                 &halofit_found_k_max),
                 pnl->error_message,
                 pnl->error_message);
          
          if (halofit_found_k_max == _TRUE_) {

						// for debugging:
            /*
							for (index_k=0; index_k<pnl->k_size; index_k++) {
								fprintf(stdout,"%e  %e  %e\n",pnl->k[index_k],pk_l[index_k],pk_nl[index_k]);
							 }
								fprintf(stdout,"\n\n");
            */
              for (index_k=0; index_k<pnl->k_size; index_k++) {
                pnl->nl_corr_density[index_pk][index_tau * pnl->k_size + index_k] = sqrt(pk_nl[index_pk][index_k]/pk_l[index_pk][index_k]);
              } 
          }
          else {
          /* when Halofit found k_max is false, use 1 as the
            non-linear correction for this redshift/time, store the
            last index which worked, and print a warning. */
            print_warning = _TRUE_;
            pnl->index_tau_min_nl = index_tau+1;
              for (index_k=0; index_k<pnl->k_size; index_k++) {
                pnl->nl_corr_density[index_pk][index_tau * pnl->k_size + index_k] = 1.;
              }
            if (pnl->nonlinear_verbose > 0) {
              class_alloc(pvecback,pba->bg_size*sizeof(double),pnl->error_message);
              class_call(background_at_tau(pba,pnl->tau[index_tau],pba->short_info,pba->inter_normal,&last_index,pvecback),
                pba->error_message,
                pnl->error_message);
              a = pvecback[pba->index_bg_a];
              z = pba->a_today/a-1.;
              fprintf(stdout,
											" -> [WARNING:] Halofit non-linear corrections could not be computed at redshift z=%5.2f and higher.\n    This is because k_max is too small for Halofit to be able to compute the scale k_NL at this redshift.\n    If non-linear corrections at such high redshift really matter for you,\n    just try to increase one of the parameters P_k_max_h/Mpc or P_k_max_1/Mpc or halofit_min_k_max (the code will take the max of these parameters) until reaching desired z.\n",z);
              free(pvecback);
            }
          }
        }
        else {
           /* if Halofit found k_max too small at a previous
              time/redhsift, use 1 as the non-linear correction for all
              higher redshifts/earlier times. */
            for (index_k=0; index_k<pnl->k_size; index_k++) {
              pnl->nl_corr_density[index_pk][index_tau * pnl->k_size + index_k] = 1.;
            }
        
        }
      }
      // get P_NL(k) at this time with HMcode
      else if (pnl->method == nl_HMcode) {				
        if (print_warning == _FALSE_) {
          if (pk_type==0){
            class_call(nonlinear_hmcode_fill_sigtab(ppr,pba,ppt,ppm,pnl,index_tau,lnk_l[index_pk],lnpk_l[index_pk],ddlnpk_l[index_pk]), 
              pnl->error_message, pnl->error_message);
          }
							/*if	(index_tau == pnl->tau_size-1) {
								fprintf(stdout, "i,  R         sigma\n");
								for (i=0;i<64;i++){
									fprintf(stdout, "%d, %e, %e\n",i, pnl->rtab[i*pnl->tau_size+index_tau]*pba->h, pnl->stab[i*pnl->tau_size+index_tau]);
								}
							} */    
            class_call(nonlinear_hmcode(ppr,
                       pba,
                       ppt,
                       ppm,
                       pnl,
                       index_pk,
                       index_tau,
                       pnl->tau[index_tau],
                       pk_l[index_pk],
                       pk_nl[index_pk],
                       lnk_l,
                       lnpk_l,
                       ddlnpk_l,
                       &(pnl->k_nl[index_pk][index_tau]),
                       &halofit_found_k_max),
              pnl->error_message,
              pnl->error_message);
          
          if (halofit_found_k_max == _TRUE_) {

						// for debugging:
						/*
							for (index_k=0; index_k<pnl->k_size_extra; index_k++) {
							if (index_tau == pnl->tau_size-1) fprintf(stdout,"%e  %e  %e\n",pnl->k_extra[index_k],pk_l[index_k],pk_l_cb[index_k]);
							}
							if (index_tau == pnl->tau_size-1) fprintf(stdout,"\n\n");
						*/
            
              for (index_k=0; index_k<pnl->k_size; index_k++) {
                pnl->nl_corr_density[index_pk][index_tau * pnl->k_size + index_k] = sqrt(pk_nl[index_pk][index_k]/pk_l[index_pk][index_k]);
              }
            
          }
          else {
						/* when HMcode found k_max is false, use 1 as the
							non-linear correction for this redshift/time, store the
							last index which worked, and print a warning. */
            print_warning = _TRUE_;
            pnl->index_tau_min_nl = index_tau+1;
            
              for (index_k=0; index_k<pnl->k_size; index_k++) {
                pnl->nl_corr_density[index_pk][index_tau * pnl->k_size + index_k] = 1.;
              }
            
            if (pnl->nonlinear_verbose > 0) {
              class_alloc(pvecback,pba->bg_size*sizeof(double),pnl->error_message);
              class_call(background_at_tau(pba,pnl->tau[index_tau],pba->short_info,pba->inter_normal,&last_index,pvecback),
                pba->error_message,
                pnl->error_message);
              a = pvecback[pba->index_bg_a];
              z = pba->a_today/a-1.;
              fprintf(stdout,
											" -> [WARNING:] HMcode non-linear corrections could not be computed at redshift z=%5.2f and higher.\n    This is because k_max is too small for HMcode to be able to compute the scale k_NL at this redshift.\n    If non-linear corrections at such high redshift really matter for you,\n    just try to increase one of the parameters P_k_max_h/Mpc or P_k_max_1/Mpc or hmcode_min_k_max (the code will take the max of these parameters) until reaching desired z.\n",z);
              free(pvecback);
            }
          }
        }
        else {
					/* if HMcode found k_max too small at a previous
							time/redhsift, use 1 as the non-linear correction for all
							higher redshifts/earlier times. */
     
            for (index_k=0; index_k<pnl->k_size; index_k++) {
              pnl->nl_corr_density[index_pk][index_tau * pnl->k_size + index_k] = 1.;
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

    } //end loop over pk_type
    
    //show the time spent for each tau:
    //clock_t end = clock();
    //double time_spent = ((double)(end - begin))/CLOCKS_PER_SEC;
    //fprintf(stdout, "tau = %e, time spent: %e s\n", pnl->tau[index_tau], time_spent);    
    
    
    } //end loop over index_tau

    for (index_pk=0; index_pk<pnl->pk_size; index_pk++){
      free(pk_l[index_pk]);
      free(pk_nl[index_pk]);
      free(lnk_l[index_pk]);
      free(lnpk_l[index_pk]);
      free(ddlnpk_l[index_pk]);    
    }
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
  /* //time spent per tau:
  clock_t end = clock();
  double time_spent = ((double)(end - begin))/CLOCKS_PER_SEC;
  fprintf(stdout, "time spent %e\n", time_spent);
  */
  return _SUCCESS_;
}

int nonlinear_free(
                   struct nonlinear *pnl
                   ) {
  int index_pk;
  
  if (pnl->method > nl_none) {

    if (pnl->method == nl_halofit) {
      free(pnl->k);
      free(pnl->tau);
      for(index_pk=0;index_pk<pnl->pk_size;++index_pk){
        free(pnl->nl_corr_density[index_pk]);
        free(pnl->k_nl[index_pk]);
      }
      free(pnl->nl_corr_density);
      free(pnl->k_nl);
    }
		else if (pnl->method == nl_HMcode){
			free(pnl->k);
      free(pnl->tau);
      for(index_pk=0;index_pk<pnl->pk_size;++index_pk){
        free(pnl->nl_corr_density[index_pk]);
        free(pnl->k_nl[index_pk]);
      }
      free(pnl->nl_corr_density);
      free(pnl->k_nl);
      free(pnl->k_extra);
			free(pnl->rtab);
			free(pnl->stab);
			free(pnl->ddstab);
			free(pnl->growtable);
			free(pnl->tautable);
			free(pnl->ztable);
      free(pnl->sigma_8);
      free(pnl->sigma_disp);
      free(pnl->sigma_disp_100);
      free(pnl->sigma_prime);
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
                   int index_pk,
                   int index_tau,
                   double *pk_l,
                   double *lnk,
                   double *lnpk,
                   double *ddlnpk) {

  int index_md;
  int index_k;
  int index_ic;
  int index_delta;
  int index_ic1,index_ic2,index_ic1_ic2;
  int k_size; 
  double * primordial_pk;
  double source_ic1,source_ic2;
  double * source_ic_extra;
//  double * source_ic_extra_cb;

  index_md = ppt->index_md_scalars;
  
  // Initialize first, then assign correct value
  index_delta = ppt->index_tp_delta_m;
  if(index_pk == pnl->index_pk_m){
    index_delta = ppt->index_tp_delta_m;
  }
    
  else if((pba->has_ncdm)&&(index_pk == pnl->index_pk_cb)){
    index_delta = ppt->index_tp_delta_cb;
  }
  else {
    class_stop(pnl->error_message,"P(k) is set neither to total matter nor to cold dark matter + baryons, index_pk=%d \n",index_pk);
  }

  class_alloc(primordial_pk,ppm->ic_ic_size[index_md]*sizeof(double),pnl->error_message);
		
	if (pnl->method == nl_HMcode){
		
		class_alloc(source_ic_extra,ppm->ic_size[index_md]*pnl->k_size_extra*sizeof(double),pnl->error_message);
    
		//fprintf(stdout, "%d\n", ppt->pk_only_cdm_bar);
		for (index_ic=0; index_ic<ppm->ic_size[index_md]; index_ic++){

			class_call(extrapolate_source(pnl->k_extra,
						      pnl->k_size,
						      pnl->k_size_extra,
						      ppt->sources[index_md][index_ic * ppt->tp_size[index_md] + index_delta]+index_tau * pnl->k_size,
						      extrapolation_only_max_units,
						      source_ic_extra+index_ic*pnl->k_size_extra,
						      pba->a_eq*pba->H_eq,
						      pba->h,																		                           pnl->error_message), 
				pnl->error_message,
				pnl->error_message)
		
    }
		
		for (index_k=0; index_k<pnl->k_size_extra; index_k++) {
			//fprintf(stdout, "%e %e\n", pnl->k_extra[index_k], source_ic_extra[index_k]);
			
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
			
      lnk[index_k] = log(pnl->k_extra[index_k]);
      lnpk[index_k] = log(pk_l[index_k]);
		}

		free(source_ic_extra);
		
    class_call(array_spline_table_columns(lnk,
                                          pnl->k_size_extra,
                                          lnpk,
                                          1,
                                          ddlnpk,
                                          _SPLINE_NATURAL_,
                                          pnl->error_message),
             pnl->error_message,
             pnl->error_message);
       
    free(primordial_pk);
        
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
            [index_ic1 * ppt->tp_size[index_md] + index_delta]
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
            [index_ic1 * ppt->tp_size[index_md] + index_delta]
            [index_tau * ppt->k_size[index_md] + index_k];

          source_ic2 = ppt->sources[index_md]
            [index_ic2 * ppt->tp_size[index_md] + index_delta]
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

//??? this array_spline table columns has to be replaced with another function
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
  }

  return _SUCCESS_;

}

int nonlinear_halofit(
                      struct precision *ppr,
                      struct background *pba,
                      struct perturbs *ppt,
                      struct primordial *ppm,
                      struct nonlinear *pnl,
                      int index_pk,
                      double tau,
                      double *pk_l,
                      double *pk_nl,
                      double *lnk_l,
                      double *lnpk_l,
                      double *ddlnpk_l,
                      double *k_nl,
                      short * halofit_found_k_max
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

  //Initialize first, then assign correct value
  fnu = pba->Omega0_ncdm_tot/Omega0_m;
  if (index_pk == pnl->index_pk_m){
    fnu = pba->Omega0_ncdm_tot/Omega0_m;
  }
  else if((pba->has_ncdm)&&(index_pk == pnl->index_pk_cb)){
    fnu = 0.;
  }
  else {
    class_stop(pnl->error_message,"P(k) is set neither to total matter nor to cold dark matter + baryons, index_pk=%d \n",index_pk);
  } 

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
    * halofit_found_k_max = _FALSE_;
    free(pvecback);
    free(integrand_array);
    return _SUCCESS_;
  }
  else {
    * halofit_found_k_max = _TRUE_;
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

/* //sigma Function with integration over k
int nonlinear_hmcode_sigma(
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

/** Calculates the sigma integral for a given scale R*/

int nonlinear_hmcode_sigma(
                  struct precision * ppr,
                  struct background * pba,
                  struct perturbs * ppt,
                  struct primordial * ppm,
                  struct nonlinear * pnl,
                  double R,
                  double *lnk_l,
                  double *lnpk_l,
                  double *ddlnpk_l,                                  
                  double * sigma
                  ) {
  double pk, lnpk;

  double * array_for_sigma;
  int index_num;
  int index_k;
  int index_sigma;
  int index_ddsigma;
  int i;
  int integrand_size;
  int last_index=0;

  double k,W,x,t;
  
  i=0;
  index_k=i;
  i++;
  index_sigma=i;
  i++;
  index_ddsigma=i;
  i++;
  index_num=i;
  
  integrand_size=(int)(log(pnl->k_extra[pnl->k_size_extra-1]/pnl->k_extra[0])/log(10.)*ppr->hmcode_k_per_decade)+1;
  class_alloc(array_for_sigma,
              integrand_size*index_num*sizeof(double),
              pnl->error_message);
          
  for (i=integrand_size-1;i>=0;i--) {
    k=pnl->k_extra[0]*pow(10.,i/ppr->hmcode_k_per_decade);
    t = 1./(1.+k);
    if (i == (integrand_size-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
    x=k*R;
    if (x<0.01) {
      W = 1.-(pow(x, 2.)/10.);
    }	
    else {
      W = 3./x/x/x*(sin(x)-x*cos(x));
    }
    
    class_call(array_interpolate_spline(lnk_l,
                                        pnl->k_size_extra,
                                        lnpk_l,
                                        ddlnpk_l,
                                        1,
                                        log(k),
                                        &last_index,
                                        &lnpk,
                                        1,
                                        pnl->error_message),
               pnl->error_message,
               pnl->error_message);
    
    pk = exp(lnpk);
    
    array_for_sigma[(integrand_size-1-i)*index_num+index_k] = t;
    array_for_sigma[(integrand_size-1-i)*index_num+index_sigma] = k*k*k*pk*W*W/(t*(1.-t));
    //if (i<pnl->k_size && R==ppr->rmin_for_sigtab/pba->h) fprintf(stdout, "%e %e\n", k, array_for_sigma[(integrand_size-1-i)*index_num+index_sigma]); 
  }

  class_call(array_spline(array_for_sigma,
                          index_num,
                          integrand_size,
                          index_k,
                          index_sigma,
                          index_ddsigma,
                          _SPLINE_EST_DERIV_,
                          pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  class_call(array_integrate_all_trapzd_or_spline(array_for_sigma,
                                        index_num,
                                        integrand_size,
                                        0, //integrand_size-1,
                                        index_k,
                                        index_sigma,
                                        index_ddsigma,
                                        sigma,
                                        pnl->error_message),
             pnl->error_message,
             pnl->error_message);

	//for (i=0;i<pnl->k_size;i++) {
		//fprintf(stdout, "%e %e %e\n", pnl->k[i], array_for_sigma[i*index_num+index_sigma], array_for_sigma[i*index_num+index_ddsigma]);    		
//	}

  free(array_for_sigma);

  *sigma = sqrt(*sigma/(2.*_PI_*_PI_));
  //fprintf(stdout, "%e\n", *sigma); 

  return _SUCCESS_;

  

}

/** Calculates the d(sigma)/dR integral for a given scale R*/

int nonlinear_hmcode_sigma_prime(
                  struct precision * ppr,
                  struct background * pba,
                  struct perturbs * ppt,
                  struct primordial * ppm,
                  struct nonlinear * pnl,
                  double R,
                  double *lnk_l,
                  double *lnpk_l,
                  double *ddlnpk_l,                                  
                  double * sigma_prime
                  ) {
  double pk, lnpk;

  double * array_for_sigma_prime;
  int index_num;
  int index_k;
  int index_sigma_prime;
  int index_ddsigma_prime;
  int integrand_size;
  int last_index=0;
  int i;

  double k,W,W_prime,x,t;
  
  i=0;
  index_k=i;
  i++;
  index_sigma_prime=i;
  i++;
  index_ddsigma_prime=i;
  i++;
  index_num=i;
  
  integrand_size=(int)(log(pnl->k_extra[pnl->k_size_extra-1]/pnl->k_extra[0])/log(10.)*ppr->hmcode_k_per_decade)+1;
  class_alloc(array_for_sigma_prime,
              integrand_size*index_num*sizeof(double),
              pnl->error_message);
            
  for (i=integrand_size-1;i>=0;i--) {
    k=pnl->k_extra[0]*pow(10.,i/ppr->hmcode_k_per_decade);
    t = 1./(1.+k);
    if (i == (integrand_size-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
    x=k*R;
	if (x<0.01) {
		W = 1.-(x*x/10.);
		W_prime = -0.2*x;
		
	}	
	else {
		 W = 3./x/x/x*(sin(x)-x*cos(x));
		 W_prime=3./x/x*sin(x)-9./x/x/x/x*(sin(x)-x*cos(x));

    }
    
    class_call(array_interpolate_spline(lnk_l,
                                        pnl->k_size_extra,
                                        lnpk_l,
                                        ddlnpk_l,
                                        1,
                                        log(k),
                                        &last_index,
                                        &lnpk,
                                        1,
                                        pnl->error_message),
               pnl->error_message,
               pnl->error_message);
    
    pk = exp(lnpk);
    
    array_for_sigma_prime[(integrand_size-1-i)*index_num+index_k] = t;
    array_for_sigma_prime[(integrand_size-1-i)*index_num+index_sigma_prime] = k*k*k*pk*2.*k*W*W_prime/(t*(1.-t));
  }

  class_call(array_spline(array_for_sigma_prime,
                          index_num,
                          integrand_size,
                          index_k,
                          index_sigma_prime,
                          index_ddsigma_prime,
                          _SPLINE_EST_DERIV_,
                          pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  class_call(array_integrate_all_trapzd_or_spline(array_for_sigma_prime,
                                        index_num,
                                        integrand_size,
                                        0, //integrand_size-1,
                                        index_k,
                                        index_sigma_prime,
                                        index_ddsigma_prime,
                                        sigma_prime,
                                        pnl->error_message),
             pnl->error_message,
             pnl->error_message);

//	for (i=0;i<integrand_size;i++) {
//		fprintf(stdout, "%e, %e, %e, %e\n", *sigma_prime, array_for_sigma_prime[i*index_num+index_k], array_for_sigma_prime[i*index_num+index_sigma_prime], array_for_sigma_prime[i*index_num+index_ddsigma_prime]);    		
//	}

  free(array_for_sigma_prime);

  *sigma_prime = *sigma_prime/(2.*_PI_*_PI_);
  //fprintf(stdout, "%e\n", *sigma_prime); 

  return _SUCCESS_;

  

}

/** Calculates the sigma_velocitydispersion integral for a given scale R*/

int nonlinear_hmcode_sigma_disp(
                  struct precision * ppr,
                  struct background * pba,
                  struct perturbs * ppt,
                  struct primordial * ppm,
                  struct nonlinear * pnl,
                  double R,
                  double *lnk_l,
                  double *lnpk_l,
                  double *ddlnpk_l,               
                  double * sigma_disp
                  ) {
  double pk, lnpk;

  double * array_for_sigma_disp;
  int index_num;
  int index_k;
  int index_y;
  int index_ddy;
  int integrand_size;
  int last_index=0;
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
  
  integrand_size=(int)(log(pnl->k_extra[pnl->k_size_extra-1]/pnl->k_extra[0])/log(10.)*ppr->hmcode_k_per_decade)+1;
  class_alloc(array_for_sigma_disp,
              integrand_size*index_num*sizeof(double),
              pnl->error_message);
            
  for (i=0;i<integrand_size;i++) {
    k=pnl->k_extra[0]*pow(10.,i/ppr->hmcode_k_per_decade);
    if (i == (integrand_size-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
    x=k*R;
	if (x<0.01) {
		W = 1.-(pow(x, 2.)/10.);
	}	
	else {
		 W = 3./x/x/x*(sin(x)-x*cos(x));
    }
    
    class_call(array_interpolate_spline(lnk_l,
                                        pnl->k_size_extra,
                                        lnpk_l,
                                        ddlnpk_l,
                                        1,
                                        log(k),
                                        &last_index,
                                        &lnpk,
                                        1,
                                        pnl->error_message),
               pnl->error_message,
               pnl->error_message);
    
    pk = exp(lnpk);
   
    
    array_for_sigma_disp[i*index_num+index_k]=k;
    array_for_sigma_disp[i*index_num+index_y]=pk*W*W;
  }

  class_call(array_spline(array_for_sigma_disp,
                          index_num,
                          integrand_size,
                          index_k,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  class_call(array_integrate_all_spline(array_for_sigma_disp,
                                        index_num,
                                        integrand_size,
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

/** Function that fills pnl->rtab, pnl->stab and pnl->ddstab with (r, sigma, ddsigma)
 * logarithmically spaced in r.
 * Called by nonlinear_init at for all tau to account for scale-dependant growth
 * before nonlinear_hmcode is called */

int nonlinear_hmcode_fill_sigtab(
              struct precision * ppr,
						  struct background * pba,
						  struct perturbs * ppt,
						  struct primordial * ppm,
						  struct nonlinear * pnl,
              int index_tau,
              double *lnk_l,
						  double *lnpk_l,
              double *ddlnpk_l             
						  ) {
	
  //double rmin;
  //double rmax;
  //int nsig;
  double r;
  double rmin, rmax;
  double sig;
  double * sigtab;
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
  
  class_alloc((sigtab),(nsig*index_n*sizeof(double)),pnl->error_message);
  
  if (index_tau == pnl->tau_size-1){
    class_alloc(pnl->rtab,nsig*sizeof(double),pnl->error_message);
    class_alloc(pnl->stab,nsig*pnl->tau_size*sizeof(double),pnl->error_message);
    class_alloc(pnl->ddstab,nsig*pnl->tau_size*sizeof(double),pnl->error_message);
  }
  
  for (i=0;i<nsig;i++){
    r=exp(log(rmin)+log(rmax/rmin)*i/(nsig-1));
    class_call(nonlinear_hmcode_sigma(ppr,pba,ppt,ppm,pnl,r,lnk_l,lnpk_l,ddlnpk_l,&sig), 
      pnl->error_message, pnl->error_message); 
    sigtab[i*index_n+index_r]=r;
    sigtab[i*index_n+index_sig]=sig;
  }
  
  class_call(array_spline(sigtab,
						  index_n,
						  nsig,
						  index_r,
						  index_sig,
						  index_ddsig,
						  _SPLINE_EST_DERIV_,
						  pnl->error_message), 
					pnl->error_message,
					pnl->error_message); 
  if (index_tau == pnl->tau_size-1){        
    for (i=0;i<nsig;i++){
      pnl->rtab[i] = sigtab[i*index_n+index_r];
      pnl->stab[i*pnl->tau_size+index_tau] = sigtab[i*index_n+index_sig];
      pnl->ddstab[i*pnl->tau_size+index_tau] = sigtab[i*index_n+index_ddsig];  
      //if (i==0) fprintf(stdout, "sigma max = %e\n", pnl->stab[i]);  
      //if (i==nsig-1) fprintf(stdout, "sigma min = %e\n", pnl->stab[i]);  
      //fprintf(stdout, "%e, %e, %e\n",sigtab[i*index_n+index_r], sigtab[i*index_n+index_sig], sigtab[i*index_n+index_ddsig]);
    }
  } 
  else{
    for (i=0;i<nsig;i++){  
      pnl->stab[i*pnl->tau_size+index_tau] = sigtab[i*index_n+index_sig];
      pnl->ddstab[i*pnl->tau_size+index_tau] = sigtab[i*index_n+index_ddsig];
    }
  }        
  
  free(sigtab); 
  
  return _SUCCESS_; 
}	


/** Function that fills pnl->tautable and pnl->growtable with (tau, D(tau))
 * linearly spaced in scalefactor a.
 * Called by nonlinear_init at z=0 before nonlinear_hmcode is called  */

int nonlinear_hmcode_fill_growtab(
              struct precision * ppr,
						  struct background * pba,
						  struct nonlinear * pnl            
						  ){
	
	double z, ainit, amax, scalefactor, tau_growth, f_class, f_camb, gamma, growth;
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
		scalefactor = log(ainit)+(log(amax)-log(ainit))*(i)/(ng-1);
		z = 1./exp(scalefactor)-1.;
		
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
    //fprintf(stdout, "%e %e\n", exp(scalefactor), pnl->growtable[i]/exp(scalefactor));      			
	}								

	free(pvecback);
	  
	return _SUCCESS_; 
}


int nonlinear_hmcode_growint(
              struct precision * ppr,
						  struct background * pba,
						  struct nonlinear * pnl,            
						  double a,
              double w0,
              double wa,
						  double * growth
						  ){
	
	double w1, z, ainit, amax, scalefactor, tau_growth, f_class, f_camb, gamma, Omega_m, Omega0_m, Omega0_v, Omega0_k, Hubble2, X_de;
  int i, index_a, index_growth, index_ddgrowth, index_gcol, last_index, ng;  
  double * pvecback;
  double * integrand;
  
  ng = 2048;
  ainit = a;
  amax = 1.;
  
  Omega0_m = (pba->Omega0_cdm + pba->Omega0_b + pba->Omega0_ncdm_tot + pba->Omega0_dcdm);
  Omega0_v = 1. - (Omega0_m + pba->Omega0_g + pba->Omega0_ur);
  Omega0_k = 1. - (Omega0_m + Omega0_v + pba->Omega0_g + pba->Omega0_ur);
  
  w1 = w0; 
  
  i=0;
  index_a = i;
  i++;
  index_growth = i;
  i++;
  index_ddgrowth = i;
  i++;
  index_gcol = i;
  
	last_index = 0;
	
  class_alloc(integrand,ng*index_gcol*sizeof(double),pnl->error_message);
  class_alloc(pvecback,pba->bg_size*sizeof(double),pnl->error_message);
  
  if (ainit == amax) {
			*growth = 1.;
	}	
	else {
  
		for (i=0;i<ng;i++){
			scalefactor = ainit+(amax-ainit)*(i)/(ng-1);
			z = 1./scalefactor-1.;
      X_de = pow(scalefactor, -3.*(1.+w0+wa))*exp(-3.*wa*(1.-scalefactor));
      Hubble2 = (Omega0_m*pow((1.+z), 3.) + Omega0_k*pow((1.+z), 2.) + Omega0_v*X_de);
      Omega_m = (Omega0_m*pow((1.+z), 3.))/Hubble2;
 
			if (w1 == -1.){
				gamma = 0.55;
			} 
			else if (w1 < -1.){
				gamma = 0.55+0.02*(1+w1);
			}
			else {
				gamma = 0.55+0.05*(1+w1);
			}
    
			integrand[i*index_gcol+index_a] = scalefactor;
			integrand[i*index_gcol+index_growth]= -pow(Omega_m, gamma)/scalefactor;
			
			class_call(array_spline(integrand,
                          index_gcol,
                          ng,
                          index_a,
                          index_growth,
                          index_ddgrowth,
                          _SPLINE_EST_DERIV_,
                          pnl->error_message),
             pnl->error_message,
             pnl->error_message);

			class_call(array_integrate_all_trapzd_or_spline(integrand,
                                        index_gcol,
                                        ng,
                                        0, //ng-1,
                                        index_a,
                                        index_growth,
                                        index_ddgrowth,
                                        growth,
                                        pnl->error_message),
             pnl->error_message,
             pnl->error_message);
    
				*growth = exp(*growth);
		}								
	}
  //fprintf(stdout, "%e %e \n", a, *growth);
	free(pvecback);
	free(integrand);
	
  return _SUCCESS_; 	
}

/** this is the fourier transform of the NFW density profile
 ** it depends on k, the virial radius rv and the concentration c
 ** of a halo  */
int nonlinear_hmcode_window_nfw(
						  					 struct nonlinear * pnl,
						  					 double k,
						  					 double rv,
						  					 double c,
												 double *window_nfw
												 ){	
	double si1, si2, ci1, ci2, ks;
	double p1, p2, p3;									 
																
	ks = k*rv/c;
	
	class_call(nonlinear_hmcode_si(
													ks*(1.+c),
													&si2
													),
					pnl->error_message, pnl->error_message);	
	
	class_call(nonlinear_hmcode_si(
													ks,
													&si1
													),
					pnl->error_message, pnl->error_message);				
	
	class_call(nonlinear_hmcode_ci(
													ks*(1.+c),
													&ci2
													),
					pnl->error_message, pnl->error_message);	
	
	class_call(nonlinear_hmcode_ci(
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

/** This is the Sheth-Tormen halo mass function (1999, MNRAS, 308, 119) */
int nonlinear_hmcode_halomassfunction(
                                      double nu, 
                                      double *hmf
                                      ){
                                        
  double p, q, A;
	
  p=0.3;
	q=0.707;
	A=0.21616;
  
  *hmf=A*(1.+(pow(q*nu*nu, -p)))*exp(-q*nu*nu/2.);
  
  return _SUCCESS_;
}


/** Computes the nonlinear correction on the linear power spectrum via 
 * the method presented in Mead et al. 1505.07833 */

int nonlinear_hmcode(
                      struct precision *ppr,
                      struct background *pba,
                      struct perturbs *ppt,
                      struct primordial *ppm,
                      struct nonlinear *pnl,
                      int index_pk,
                      int index_tau,
                      double tau,
                      double *pk_l,                      
                      double *pk_nl,
                      double **lnk_l,
                      double **lnpk_l,
                      double **ddlnpk_l,
                      double *k_nl,
                      short * halofit_found_k_max                        
                      ) {
  
  /* integers */
  int n, i, ng, nsig;
  int index_k, index_lnsig, index_dlnsig, index_ncol;
  int last_index=0;  
  int index_ncdm;
  int index_pk_cb;
  int counter, index_nl;
  
	int index_nu;
  int index_y;
  int index_ddy;
  
  /* Background parameters */
  double Omega_m,Omega_v,fnu,Omega0_m, w0, dw_over_da_fld, integral_fld;	
  double z_at_tau;
  double rho_crit_today_in_msun_mpc3, rho_crit_today_in_msun_mpc3_class;
  double growth;  
  double anorm;
  
  /* temporary numbers */ 
  double m, r, nu, sig, sigf;
	double diff, rmid, r1, r2;
  
  /* HMcode parameters */
  double mmin, mmax, nu_min;
  
  double sigma_disp, sigma_disp100, sigma8;
  double delta_c, Delta_v;
  double fraction;
  
  double sigma_nl, nu_nl, r_nl;
  double sigma_prime;
  double dlnsigdlnR;  
  double n_eff; 
  double alpha;

  double z_form, g_form;
  double z_infinity, tau_infinity;
  double g_lcdm, g_wcdm;
  double de_correction;
  
  double eta;
  double gst, window_nfw;
  double fac, k_star, fdamp;
  double pk_lin, pk_2h, pk_1h;
  
  /* data fields */ 
  double * rtab;
  double * stab;
  double * ddstab;
  
  double * pvecback;  
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
  
  double * p1h_integrand; 
  
  
  /** include precision parameters that control the number of entries in the growth and sigma tables */
  ng = ppr->n_hmcode_tables;
  nsig = ppr->n_hmcode_tables;

  
  /** Call all the relevant background parameters at this tau */
  class_alloc(pvecback,pba->bg_size*sizeof(double),pnl->error_message);

  Omega0_m = (pba->Omega0_cdm + pba->Omega0_b + pba->Omega0_ncdm_tot + pba->Omega0_dcdm);
  class_call(background_w_fld(pba,pba->a_today,&w0,&dw_over_da_fld,&integral_fld), pba->error_message, pnl->error_message);

  fnu      = pba->Omega0_ncdm_tot/Omega0_m;
  anorm    = 1./(2*pow(_PI_,2));  
  
  class_call(background_w_fld(pba,pba->a_today,&w0,&dw_over_da_fld,&integral_fld), pba->error_message, pnl->error_message);
  
  class_call(background_at_tau(pba,tau,pba->long_info,pba->inter_normal,&last_index,pvecback),
             pba->error_message,
             pnl->error_message);

  Omega_m = pvecback[pba->index_bg_Omega_m];
  Omega_v = 1.-pvecback[pba->index_bg_Omega_m]-pvecback[pba->index_bg_Omega_r];


  growth = pvecback[pba->index_bg_D];

  z_at_tau = 1./pvecback[pba->index_bg_a]-1.;
  
  /* The number below does the unit conversion to solar masses over Mpc^3: 2.77474589e11=rho_c/h^2 [Msun/Mpc^3] with rho_c = 8*pi*G/3*c^2 */
  rho_crit_today_in_msun_mpc3 = 2.77474589e11*pow(pba->h, 2); 
  
  free(pvecback);
  
  /** Calculate the Dark Energy correction: */  
  if (pba->has_fld==_TRUE_){
    if (index_tau == pnl->tau_size-1){         
      class_call(nonlinear_hmcode_growint(ppr,pba,pnl,1./(1.+pnl->z_infinity),-1.,0.,&g_lcdm ),
        pnl->error_message, pnl->error_message);
      class_call(nonlinear_hmcode_growint(ppr,pba,pnl,1./(1.+pnl->z_infinity),w0,dw_over_da_fld*(-1),&g_wcdm ),
        pnl->error_message, pnl->error_message);
      pnl->dark_energy_correction = pow(g_wcdm/g_lcdm, 1.5);
      //fprintf(stdout, "%e %e %e %e\n", dw_over_da_fld, g_wcdm, g_lcdm, pnl->dark_energy_correction); 
    }
  }
  else {
    pnl->dark_energy_correction = 1.;
  }
  /** Test whether pk_cb has to be taken into account (only if we have massive neutrinos)*/
  if (pba->has_ncdm==_TRUE_){
    index_pk_cb = pnl->index_pk_cb;
  }
  else {
    index_pk_cb = index_pk;
  }
  
  //double sigma_test;  
  /** Get sigma(R=8 Mpc/h), sigma_disp(R=0), sigma_disp(R=100 Mpc/h) and write the into pnl structure */
  class_call(nonlinear_hmcode_sigma(ppr,pba,ppt,ppm,pnl,8./pba->h,lnk_l[index_pk],lnpk_l[index_pk],ddlnpk_l[index_pk],&sigma8), 
			pnl->error_message, pnl->error_message);	
  class_call(nonlinear_hmcode_sigma_disp(ppr,pba,ppt,ppm,pnl,0.,lnk_l[index_pk],lnpk_l[index_pk],ddlnpk_l[index_pk],&sigma_disp), 
			pnl->error_message, pnl->error_message);
  class_call(nonlinear_hmcode_sigma_disp(ppr,pba,ppt,ppm,pnl,100./pba->h,lnk_l[index_pk],lnpk_l[index_pk],ddlnpk_l[index_pk],&sigma_disp100), 
			pnl->error_message, pnl->error_message);			
  /*class_call(nonlinear_hmcode_sigma(ppr,pba,ppt,ppm,pnl,4.,lnk_l[index_pk],lnpk_l[index_pk],ddlnpk_l[index_pk],&sigma_test), 
			pnl->error_message, pnl->error_message);
  fprintf(stdout, "z: %e, sigma: %e\n", z_at_tau, sigma_test); */   
  pnl->sigma_8[index_tau] = sigma8;
  pnl->sigma_disp[index_tau] = sigma_disp;
  pnl->sigma_disp_100[index_tau] = sigma_disp100;
  
  //fprintf(stdout, "%e %e\n", z_at_tau, sigma8);
  
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
  
  class_alloc(stab,ppr->n_hmcode_tables*sizeof(double),pnl->error_message);
  class_alloc(ddstab,ppr->n_hmcode_tables*sizeof(double),pnl->error_message);
  
  for (i=0;i<ppr->n_hmcode_tables;i++){
    stab[i] = pnl->stab[i*pnl->tau_size+index_tau];
    ddstab[i] = pnl->ddstab[i*pnl->tau_size+index_tau];
  }  
  
  // Linear theory density perturbation threshold for spherical collapse
  delta_c = 1.59+0.0314*log(sigma8); //Mead et al. (2015; arXiv 1505.07833)
  delta_c = delta_c*(1.+0.0123*log10(Omega_m)); //Nakamura & Suto (1997) fitting formula for LCDM models
  delta_c = delta_c*(1.+0.262*fnu); //Mead et al. (2016; arXiv 1602.02154) neutrino addition
  
  // virialized overdensity
  Delta_v=418.*pow(Omega_m, -0.352); //Mead et al. (2015; arXiv 1505.07833)          
  Delta_v=Delta_v*(1.+0.916*fnu); //Mead et al. (2016; arXiv 1602.02154) neutrino addition
  
  // mass or radius fraction respectively
  fraction = pow(0.01, 1./3.);
	
	/* Fill the arrays needed for the P1H Integral: mass, r_real, r_virial, nu_arr, sigma_r, sigmaf_r
   * The P1H Integral is an integral over nu=delta_c/sigma(M), where M is connected to R via R=(3M)/(4*pi*rho_m).
   * The Integrand is M*Window^2{nu(M)*k, Rv(M), c(M)}*f(nu) with the window being the fouriertransformed
   * NFW profile, Rv = R/Delta_v^(1/3) and Sheth-Thormen halo mass function f.
   * The halo concentration-mass-relation c(M) will be found later.  */
  for (i=0;i<n;i++){
	  m = exp(log(mmin)+log(mmax/mmin)*(i)/(n-1));
	  r = pow((3.*m/(4.*_PI_*rho_crit_today_in_msun_mpc3*Omega0_m)), (1./3.)); 
	  mass[i] = m;
	  r_real[i] = r;
	  r_virial[i] = r_real[i]/pow(Delta_v, 1./3.);
	  class_call(array_interpolate_spline(pnl->rtab,
										  nsig,
										  stab,
										  ddstab,
										  1,
										  r,
										  &last_index,
										  &sig,
										  1,
										  pnl->error_message),
					pnl->error_message, pnl->error_message);
	  class_call(array_interpolate_spline(pnl->rtab,
										  nsig,
										  stab,
										  ddstab,
										  1,
										  r*fraction,
										  &last_index,
										  &sigf,
										  1,
										  pnl->error_message),
					pnl->error_message, pnl->error_message);

		nu=delta_c/sig;
	  sigma_r[i] = sig;
	  sigmaf_r[i] = sigf;
	  nu_arr[i] = nu;				
  }
  
  free(stab);
  free(ddstab);
  
  /** find nonlinear scales k_nl and r_nl and the effective spectral index n_eff */
  nu_nl = 1.;
  nu_min = nu_arr[0];
	//fprintf(stdout, "z: %e, numin: %e\n", z_at_tau, nu_min);
  // stop calculating the nonlinear correction if the nonlinear scale is not reached in the table:
	if (nu_min > nu_nl) {
    class_stop(pnl->error_message,"Your input of the mass range is such that the nonlinear scale cannot be found at this redshift %5.2. Try to decrease the parameter mmin_for_p1h_integral\n",z_at_tau);
  }

  
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

	class_call(array_search_bisect(n,nu_arr,nu_nl,&index_nl,pnl->error_message), pnl->error_message, pnl->error_message);
	
	r1 = r_real[index_nl-1];
	r2 = r_real[index_nl+2];
  
  // for debugging: (if it happens that r_nl is not between r1 and r2, which should never be the case) 
  //fprintf(stdout, "%e %e %e %e\n", r1, nu_arr[index_nl-1], r2, nu_arr[index_nl+2]);
	
  /* do bisectional iteration between r1 and r2 to find the precise value of r_nl */
  counter = 0;
  do {
    r_nl = (r1+r2)/2.;
    counter ++;
	  
		class_call(nonlinear_hmcode_sigma(ppr,pba,ppt,ppm,pnl,r_nl,lnk_l[index_pk_cb],lnpk_l[index_pk_cb],ddlnpk_l[index_pk_cb],&sigma_nl),
			pnl->error_message, pnl->error_message);		

    diff = sigma_nl - delta_c;

    if (diff > ppr->halofit_tol_sigma){
      r1=r_nl;
    }
    else if (diff < -ppr->halofit_tol_sigma) {
      r2 = r_nl;
    }

    class_test(counter > _MAX_IT_,
               pnl->error_message,
               "could not converge within maximum allowed number of iterations");

  } while (fabs(diff) > ppr->halofit_tol_sigma);
	
  if (pnl->nonlinear_verbose>5){
    fprintf(stdout, "number of iterations for r_nl at z = %e: %d\n", z_at_tau, counter);
  }
	*k_nl = 1./r_nl;
  
  if (*k_nl > pnl->k[pnl->k_size-1]) {
    * halofit_found_k_max = _FALSE_;
    free(mass);
    free(r_real);
    free(r_virial);
    free(sigma_r);
    free(sigmaf_r);
    free(nu_arr);
    return _SUCCESS_;
  }
	else {
		* halofit_found_k_max = _TRUE_;
	}
  
  /* call sigma_prime function at r_nl to find the effective spectral index n_eff */
	class_call(nonlinear_hmcode_sigma_prime(ppr,pba,ppt,ppm,pnl,r_nl,lnk_l[index_pk_cb],lnpk_l[index_pk_cb],ddlnpk_l[index_pk_cb],&sigma_prime),
			pnl->error_message, pnl->error_message);
	dlnsigdlnR = r_nl*pow(sigma_nl, -2)*sigma_prime;
  n_eff = -3.- dlnsigdlnR;
  alpha = 3.24*pow(1.85, n_eff);
  
  pnl->sigma_prime[index_tau] = sigma_prime;
  
  /** Calculate halo concentration-mass relation conc(mass) */ 
	class_alloc(conc,n*sizeof(double),pnl->error_message);

  if (tau==pba->conformal_age && pnl->nonlinear_verbose > 9) {
    fprintf(stdout, "#################################################################\n");
    fprintf(stdout, "tau = %f; z = %f\n", tau, z_at_tau);
    fprintf(stdout, "--------------------------------------------------------------------\n");
    fprintf(stdout, "M,		rf,		s(fM),		zf,		gf,		c \n");
    fprintf(stdout, "--------------------------------------------------------------------\n");
  }
    
  for (i=0;i<n;i++){
		//find growth rate at formation
    g_form = delta_c*growth/sigmaf_r[i];
    if (g_form > 1.) g_form = 1.;
		
    // 
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
			conc[i] = pnl->c_min;
		} else {
			conc[i] = pnl->c_min*(1.+z_form)/(1.+z_at_tau)*pnl->dark_energy_correction;
	  } 	
		if (tau==pba->conformal_age && pnl->nonlinear_verbose > 9) fprintf(stdout, "%e %e %e %e %e %e\n",mass[i], r_real[i]*fraction*pba->h, sigmaf_r[i], z_form, g_form, conc[i]);
	}
  if (tau==pba->conformal_age && pnl->nonlinear_verbose > 9) fprintf(stdout, "#################################################################\n");
  
  
  
  
	/** Now compute the nonlinear correction */
  eta = pnl->eta_0 - 0.3*sigma8; //halo bloating parameter
  k_star=0.584/sigma_disp;   //Damping wavenumber of the 1-halo term at very large scales;
	fdamp = 0.0095*pow(sigma_disp100*pba->h, 1.37); //Damping factor for 2-halo term 
	if (fdamp<1.e-3) fdamp=1.e-3;
  if (fdamp>0.99)  fdamp=0.99;
  
  double nu_cut = 10.;
  int index_cut;
  class_call(array_search_bisect(n,nu_arr,nu_cut,&index_cut,pnl->error_message), pnl->error_message, pnl->error_message);
  
  i=0;
  index_nu=i;
  i++;
  index_y=i;
  i++;
  index_ddy=i;
  i++;
  index_ncol=i;
	
	for (index_k = 0; index_k < pnl->k_size; index_k++){
		
		class_alloc(p1h_integrand,index_cut*index_ncol*sizeof(double),pnl->error_message);
		
		pk_lin = pk_l[index_k]*pow(pnl->k[index_k],3)*anorm; //convert P_k to Delta_k^2
		
		for (i=0; i<index_cut; i++){ //Calculates the integrand for the ph1 integral at all nu values
			//get the nu^eta-value of the window
      class_call(nonlinear_hmcode_window_nfw(
																	pnl,
																	pow(nu_arr[i], eta)*pnl->k[index_k],
																	r_virial[i],
																	conc[i],
																  &window_nfw),
					pnl->error_message, pnl->error_message);	
			//get the value of the halo mass function
      class_call(nonlinear_hmcode_halomassfunction(
																	nu_arr[i],
																  &gst),
					pnl->error_message, pnl->error_message);	

			p1h_integrand[i*index_ncol+index_nu] = nu_arr[i];
			
			p1h_integrand[i*index_ncol+index_y] = mass[i]*gst*pow(window_nfw, 2.);
      //if ((tau==pba->conformal_age) && (index_k == 0)) {
        //fprintf(stdout, "%d %e %e\n", index_cut, p1h_integrand[i*index_ncol+index_nu], p1h_integrand[i*index_ncol+index_y]);
			//}
		} 
		class_call(array_spline(p1h_integrand,
                            index_ncol,
                            index_cut,
                            index_nu,
                            index_y,
                            index_ddy,
                            _SPLINE_EST_DERIV_,
                            pnl->error_message),
             pnl->error_message,
             pnl->error_message);
		
		class_call(array_integrate_all_trapzd_or_spline(
																				p1h_integrand,
                                        index_ncol,
                                        index_cut, 
                                        index_cut-1, //0 or n-1
                                        index_nu,
                                        index_y,
                                        index_ddy,
                                        &pk_1h,
                                        pnl->error_message),
             pnl->error_message,
             pnl->error_message);
    

    if (pow(pnl->k[index_k]/k_star, 2)>7.){
      fac = 0.;     //prevents problems if (k/k*)^2 is large
    }
    else{
      fac=exp(-pow((pnl->k[index_k]/k_star), 2.));
    }
    /*if (tau==pba->conformal_age){
      fprintf(stdout, "%e %e %e\n",pnl->k[index_k]/pba->h, fac, pk_1h);
    }*/
    pk_1h = pk_1h*anorm*pow(pnl->k[index_k],3)*(1.-fac)/(rho_crit_today_in_msun_mpc3*Omega0_m);  // dimensionless power
		
		if (fdamp==0){
			pk_2h=pk_lin;
		}else{	
			pk_2h=pk_lin*(1.-fdamp*pow(tanh(pnl->k[index_k]*sigma_disp/sqrt(fdamp)), 2.)); //dimensionless power
		}
		if (pk_2h<0.) pk_2h=0.;	
		pk_nl[index_k] = pow((pow(pk_1h, alpha) + pow(pk_2h, alpha)), (1./alpha))/pow(pnl->k[index_k],3)/anorm; //converted back to P_k
		
    // for debugging:
    //if (tau==pba->conformal_age) fprintf(stdout, "%e %e %e %e %e\n", pnl->k[index_k], pk_lin, pk_1h, pk_2h, pk_nl[index_k]*pow(pnl->k[index_k],3)*anorm);			
		//if (tau==pba->conformal_age) fprintf(stdout, "%e, %e\n",pnl->k[index_k]/pba->h, pk_nl[index_k]);
		
		free(p1h_integrand);
	}
	
  // print parameter values  
  if ((pnl->nonlinear_verbose > 1 && tau==pba->conformal_age) || pnl->nonlinear_verbose > 5){
		fprintf(stdout, " -> Parameters at redshift z = %e:\n", z_at_tau); 
		fprintf(stdout, "    fnu:		%e\n", fnu);
		fprintf(stdout, "    sigd [Mpc/h]:	%e\n", sigma_disp*pba->h);
		fprintf(stdout, "    sigd100 [Mpc/h]:    %e\n", sigma_disp100*pba->h);
		fprintf(stdout, "    sigma8:		%e\n", sigma8);
		fprintf(stdout, "    nu min:		%e\n", nu_arr[0]);
		fprintf(stdout, "    nu max:		%e\n", nu_arr[n-1]);
		fprintf(stdout, "    r_v min [Mpc/h]:    %e\n", r_virial[0]*pba->h);
		fprintf(stdout, "    r_v max [Mpc/h]:    %e\n", r_virial[n-1]*pba->h);				
		fprintf(stdout, "    r_nl [Mpc/h]:	%e\n", r_nl*pba->h);
		fprintf(stdout, "    k_nl [h/Mpc]:	%e\n", *k_nl/pba->h);
		fprintf(stdout, "    sigma_nl:		%e\n", sigma_nl/delta_c);
		fprintf(stdout, "    neff:		%e\n", n_eff);
		fprintf(stdout, "    c min:		%e\n", conc[n-1]);
		fprintf(stdout, "    c max:		%e\n", conc[0]);			
		fprintf(stdout, "    Dv:			%e\n", Delta_v);
		fprintf(stdout, "    dc:			%e\n", delta_c);	
		fprintf(stdout, "    eta:		%e\n", eta);	
		fprintf(stdout, "    k*:			%e\n", k_star/pba->h);
		fprintf(stdout, "    Abary:		%e\n", pnl->c_min);			
		fprintf(stdout, "    fdamp:		%e\n", fdamp);		
		fprintf(stdout, "    alpha:		%e\n", alpha);
		fprintf(stdout, "    ksize, kmin, kmax:   %d, %e, %e\n", pnl->k_size, pnl->k[0]/pba->h, pnl->k[pnl->k_size-1]/pba->h);	
		
    // for debugging:
    //for (i=0;i<ppr->n_hmcode_tables;i++){
			//fprintf(stdout, "%d %e %e %e\n",i+1, pnl->rtab[i*pnl->tau_size+index_tau]*pba->h, pnl->stab[i*pnl->tau_size+index_tau], pnl->ddstab[i*pnl->tau_size+index_tau]);
			//fprintf(stdout, "%e %e\n",1./(1.+pnl->ztable[i]), pnl->growtable[i]);
  	//}
  	
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
  





