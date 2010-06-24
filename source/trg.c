/** @file trg.c Document Time Renormalization Group module
 *  Benjamin Audren, 24.06.2010
 *
 *  Calculates the non linear matter spectra P_k
 *
 **/

#include "precision.h"
#include "background.h"
#include "thermodynamics.h"
#include "perturbations.h"
#include "bessel.h"
#include "transfer.h"
#include "primordial.h"
#include "spectra.h"
#include "trg.h"
#include "time.h"
#ifdef _OPENMP
#include "omp.h"
#endif


/********************
 *
 * Definition of gamma functions, carrying all non-linear information
 *
 ********************/

int trg_gamma_121(
		  double  k, 
		  double  q, 
		  double  p, 
		  double  *result
		  ){
  *result =  (1./ (4*q*q)) * (- p*p + k*k + q*q);
  return _SUCCESS_;
}



int trg_gamma_222(
		  double k, 
		  double p, 
		  double q,
		  double * result,
		  char * errmsg
		  ){
  *result=k*k/(4*q*q*p*p) * (k*k - p*p - q*q);
  return _SUCCESS_;
}


/********************
 *
 * Calls the spectrum computed by spectra.c at any time eta, and
 * multiplies it by the exponential to get to the same unit as the
 * rest of the code.  
 * It serves at the initialisation for both 1-loop
 * et TRG, and is after only used by the 1-loop computation to input
 * in the A integrals.
 *
 ********************/


int trg_p11_at_k(
		 struct background * pba,
		 struct primordial * ppm,
		 struct spectra * psp,
		 struct spectra_nl * pnl,
		 int index_eta,
		 int index_ic,
		 double k,
		 double *result
		 ){

  class_call(spectra_pk_at_k_and_z(pba,ppm,psp,0,index_ic,k,pnl->z[index_eta],result),
	     psp->error_message,
	     pnl->error_message);

  *result=exp(-2*pnl->eta[index_eta])* *result;
  
  return _SUCCESS_;
}

/********************
 *
 * As soon as one leaves the initial condition, one needs to
 * interpolate pk_nl, p_12_nl and p_22_nl from the already computed
 * non linear spectrum. This is done by the following two commands
 *
 ********************/

/*====>*/

   /********************
    *
    * Fill the second derivative table of p_ab up to the index_eta
    * subscript
    *
    ********************/

int trg_ddp_ab(
	       struct spectra_nl * pnl,
	       double *p_ab,
	       int index_eta,
	       double *result,
	       char *errmsg
	       ){
  class_call(array_spline_table_one_column(pnl->k,pnl->k_size,p_ab,pnl->eta_size,index_eta,result,_SPLINE_EST_DERIV_,errmsg),
	     errmsg,
	     errmsg);
  return _SUCCESS_;

}

    /********************
     *
     * Interpolate any structure 
     *
     ********************/

int trg_p_ab_at_any_k(
		      struct spectra_nl * pnl,
		      double *p_ab, /*filled until [any_index_k+k_size*index_eta]*/
		      double *ddp_ab,/*filled until the same index*/
		      int index_eta,
		      double any_k,
		      double *result, /* will get a value of p_ab at any k different from k[index_k]*/
		      char *errmsg
		      ){

  class_call(array_interpolate_spline_one_column(pnl->k,pnl->k_size,p_ab,pnl->eta_size,index_eta,ddp_ab,any_k,result,errmsg),
	     errmsg,
	     pnl->error_message);

  return _SUCCESS_;
}

/********************
 *
 * Argument of the integral of all A's, B's, given at any time eta,
 * usually written as F( k,x+y/sqrt(2),x-y/sqrt(2) ) + F( y <-> -y) p
 * stands for x+y/sqrt(2) (plus) and m for x-y/sqrt(2) (minus).  This
 * function is used when mode is set to 2, to perform full TRG
 * computation.
 *
 ********************/

 

int trg_A_arg(
	      struct spectra_nl * pnl,
	      enum name_A name, 
	      double k, 
	      double p, 
	      double m, 
	      int index_eta,
	      int index_k, /**< used only for testing, could be supressed */
	      double * result, 
	      char * errmsg){

  double p_11k,p_11m,p_11p;
  double p_12k,p_12m,p_12p;
  double p_22k,p_22m,p_22p;

  double gamma1_kpm,gamma1_pkm,gamma1_mkp,gamma1_mpk,gamma1_pmk,gamma1_kmp;
  double gamma2_kpm,gamma2_pkm,gamma2_mkp,gamma2_mpk,gamma2_pmk,gamma2_kmp;

  /**********
   * Under and above the two chosen cut-off k_max and k_min, force the
   * result to 0.
   * m can only be under k_min, while p might be under
   * k_min and above k_max.
   **********/

  if (m < pnl->k[0]) {
    *result=0.;
    return _SUCCESS_;
  }


  if ((p<pnl->k[0]) || (p>pnl->k[pnl->k_size-1])) {
    *result=0.;
    return _SUCCESS_;
  }

  /**********
   * Every function needed for the argument of a certain integral is
   * only called inside the switch case instance, in order to lower the
   * execution time - this way the lecture of this function is somehow
   * harder, but the argument is called in the very end of each case
   * bloc. Moreover, the call of functions is in this order : P_ab, the
   * gamma2, then finally gamma1.
   **********/

  switch(name){
  case _A0_:
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,k,&p_22k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,p,&p_22p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,m,&p_22m,errmsg),
	       errmsg,
	       pnl->error_message);
      
    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,errmsg),
	       errmsg,
	       pnl->error_message);

    trg_gamma_121(k,p,m,&gamma1_kpm);
    trg_gamma_121(k,m,p,&gamma1_kmp);
   
   
    *result= gamma1_kpm *( gamma2_kpm*p_22p*p_22m + gamma2_pmk*p_22m*p_22k + 
			   gamma2_mkp*p_22k*p_22p )
      + gamma1_kmp *( gamma2_kmp*p_22m*p_22p + gamma2_mpk*p_22p*p_22k + 
		      gamma2_pkm*p_22k*p_22m );

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  case _A11_:
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,k,&p_12k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,p,&p_12p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,m,&p_12m,errmsg),
	       errmsg,
	       pnl->error_message);

    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,k,&p_22k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,p,&p_22p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,m,&p_22m,errmsg),
	       errmsg,
	       pnl->error_message);

    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,errmsg),
	       errmsg,
	       pnl->error_message);

    trg_gamma_121(k,p,m,&gamma1_kpm);
    trg_gamma_121(k,m,p,&gamma1_kmp);

       


    *result= gamma1_kpm *( gamma1_kpm*p_22p*p_12m + gamma1_kmp*p_12p*p_22m + 
			   gamma2_pmk*p_22m*p_12k + gamma2_mkp*p_12k*p_22p)
      + gamma1_kmp *( gamma1_kmp*p_22m*p_12p + gamma1_kpm*p_12m*p_22p + 
		      gamma2_mpk*p_22p*p_12k + gamma2_pkm*p_12k*p_22m);

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  case _A12_:
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,k,&p_12k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,p,&p_12p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,m,&p_12m,errmsg),
	       errmsg,
	       pnl->error_message);

    
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,k,&p_22k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,p,&p_22p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,m,&p_22m,errmsg),
	       errmsg,
	       pnl->error_message);
    
   
    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,errmsg),
	       errmsg,
	       pnl->error_message);
  
    
    trg_gamma_121(k,p,m,&gamma1_kpm);
    trg_gamma_121(k,m,p,&gamma1_kmp);
    trg_gamma_121(p,k,m,&gamma1_pkm);
    trg_gamma_121(p,m,k,&gamma1_pmk);
    trg_gamma_121(m,p,k,&gamma1_mpk);
    trg_gamma_121(m,k,p,&gamma1_mkp);

    *result= gamma1_kpm *( gamma2_kpm*p_12p*p_22m + gamma1_pmk*p_22m*p_12k +
			   gamma1_pkm*p_12m*p_22k + gamma2_mkp*p_22k*p_12p)
      + gamma1_kmp *( gamma2_kmp*p_12m*p_22p + gamma1_mpk*p_22p*p_12k +
		      gamma1_mkp*p_12p*p_22k + gamma2_pkm*p_22k*p_12m);

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  case _A13_:


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,k,&p_12k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,p,&p_12p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,m,&p_12m,errmsg),
	       errmsg,
	       pnl->error_message);

    
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,k,&p_22k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,p,&p_22p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,m,&p_22m,errmsg),
	       errmsg,
	       pnl->error_message);
    
   
    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,errmsg),
	       errmsg,
	       pnl->error_message);
  
    
    trg_gamma_121(k,p,m,&gamma1_kpm);
    trg_gamma_121(k,m,p,&gamma1_kmp);
    trg_gamma_121(p,k,m,&gamma1_pkm);
    trg_gamma_121(p,m,k,&gamma1_pmk);
    trg_gamma_121(m,p,k,&gamma1_mpk);
    trg_gamma_121(m,k,p,&gamma1_mkp);

    *result= gamma1_kpm *( gamma2_kpm*p_22p*p_12m + gamma2_pmk*p_12m*p_22k +
			   gamma1_mkp*p_22k*p_12p + gamma1_mpk*p_12k*p_22p)
      + gamma1_kmp *( gamma2_kmp*p_22m*p_12p + gamma2_mpk*p_12p*p_22k +
		      gamma1_pkm*p_22k*p_12m + gamma1_pmk*p_12k*p_22m);

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  case _A21_:
    

    class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,p,&p_11p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,m,&p_11m,errmsg),
	       errmsg,
	       pnl->error_message);
    

    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,k,&p_12k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,p,&p_12p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,m,&p_12m,errmsg),
	       errmsg,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,k,&p_22k,errmsg),
	       errmsg,
	       pnl->error_message);

    
    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,errmsg),
	       errmsg,
	       pnl->error_message);
    
    
    trg_gamma_121(k,p,m,&gamma1_kpm);
    trg_gamma_121(k,m,p,&gamma1_kmp);
    trg_gamma_121(p,k,m,&gamma1_pkm);
    trg_gamma_121(p,m,k,&gamma1_pmk);
    trg_gamma_121(m,p,k,&gamma1_mpk);
    trg_gamma_121(m,k,p,&gamma1_mkp);

    *result= gamma1_kpm *( gamma2_kpm*p_12p*p_12m + gamma1_pmk*p_12m*p_12k +
			   gamma1_pkm*p_11m*p_22k + gamma1_mkp*p_22k*p_11p +
			   gamma1_mpk*p_12k*p_12p)
      + gamma1_kmp *( gamma2_kmp*p_12m*p_12p + gamma1_mpk*p_12p*p_12k +
		      gamma1_mkp*p_11p*p_22k + gamma1_pkm*p_22k*p_11m +
		      gamma1_pmk*p_12k*p_12m);

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  case _A22_:
    class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,k,&p_11k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,p,&p_11p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,m,&p_11m,errmsg),
	       errmsg,
	       pnl->error_message);
    

    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,k,&p_12k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,p,&p_12p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,m,&p_12m,errmsg),
	       errmsg,
	       pnl->error_message);

    
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,p,&p_22p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,m,&p_22m,errmsg),
	       errmsg,
	       pnl->error_message);
      
   
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,errmsg),
	       errmsg,
	       pnl->error_message);
  
    
    trg_gamma_121(k,p,m,&gamma1_kpm);
    trg_gamma_121(k,m,p,&gamma1_kmp);
    trg_gamma_121(p,k,m,&gamma1_pkm);
    trg_gamma_121(p,m,k,&gamma1_pmk);
    trg_gamma_121(m,p,k,&gamma1_mpk);
    trg_gamma_121(m,k,p,&gamma1_mkp);


    *result= gamma1_kpm *( gamma1_kpm*p_22p*p_11m + gamma1_kmp*p_12p*p_12m +
			   gamma2_pmk*p_12m*p_12k + gamma1_mkp*p_12k*p_12p +
			   gamma1_mpk*p_11k*p_22p)
      + gamma1_kmp *( gamma1_kmp*p_22m*p_11p + gamma1_kpm*p_12m*p_12p +
		      gamma2_mpk*p_12p*p_12k + gamma1_pkm*p_12k*p_12m +
		      gamma1_pmk*p_11k*p_22m);

    *result *= m*p/2./pow(2*_PI_,3); 

    return _SUCCESS_;
    break;

    /****************************************/

  case _A23_:
    class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,k,&p_11k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,p,&p_11p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,m,&p_11m,errmsg),
	       errmsg,
	       pnl->error_message);
    

    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,k,&p_12k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,p,&p_12p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,m,&p_12m,errmsg),
	       errmsg,
	       pnl->error_message);

    
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,p,&p_22p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,m,&p_22m,errmsg),
	       errmsg,
	       pnl->error_message);
      
   
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,errmsg),
	       errmsg,
	       pnl->error_message);
  
    
    trg_gamma_121(k,p,m,&gamma1_kpm);
    trg_gamma_121(k,m,p,&gamma1_kmp);
    trg_gamma_121(p,k,m,&gamma1_pkm);
    trg_gamma_121(p,m,k,&gamma1_pmk);
    trg_gamma_121(m,p,k,&gamma1_mpk);
    trg_gamma_121(m,k,p,&gamma1_mkp);


    *result= gamma1_kpm *( gamma1_kpm*p_12p*p_12m + gamma1_kmp*p_11p*p_22m +
			   gamma1_pmk*p_22m*p_11k + gamma1_pkm*p_12m*p_12k +
			   gamma2_mkp*p_12k*p_12p)
      + gamma1_kmp *( gamma1_kmp*p_12m*p_12p + gamma1_kpm*p_11m*p_22p +
		      gamma1_mpk*p_22p*p_11k + gamma1_mkp*p_12p*p_12k +
		      gamma2_pkm*p_12k*p_12m);

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  case _A3_:
    class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,k,&p_11k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,p,&p_11p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,m,&p_11m,errmsg),
	       errmsg,
	       pnl->error_message);
    

    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,k,&p_12k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,p,&p_12p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,m,&p_12m,errmsg),
	       errmsg,
	       pnl->error_message);

    
    trg_gamma_121(k,p,m,&gamma1_kpm);
    trg_gamma_121(k,m,p,&gamma1_kmp);
    trg_gamma_121(p,k,m,&gamma1_pkm);
    trg_gamma_121(p,m,k,&gamma1_pmk);
    trg_gamma_121(m,p,k,&gamma1_mpk);
    trg_gamma_121(m,k,p,&gamma1_mkp);


    *result= (gamma1_kpm+gamma1_kmp) *(
				       gamma1_kpm*p_12p*p_11m + gamma1_kmp*p_11p*p_12m +
				       gamma1_pmk*p_12m*p_11k + gamma1_pkm*p_11m*p_12k +
				       gamma1_mkp*p_12k*p_11p + gamma1_mpk*p_11k*p_12p);

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  case _B0_:
    class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,k,&p_11k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,p,&p_11p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,m,&p_11m,errmsg),
	       errmsg,
	       pnl->error_message);
    

    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,k,&p_12k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,p,&p_12p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,m,&p_12m,errmsg),
	       errmsg,
	       pnl->error_message);

    
    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,errmsg),
	       errmsg,
	       pnl->error_message);
    
    
    trg_gamma_121(k,p,m,&gamma1_kpm);
    trg_gamma_121(k,m,p,&gamma1_kmp);
    trg_gamma_121(p,k,m,&gamma1_pkm);
    trg_gamma_121(p,m,k,&gamma1_pmk);
    trg_gamma_121(m,p,k,&gamma1_mpk);
    trg_gamma_121(m,k,p,&gamma1_mkp);

    *result= (gamma2_kpm+gamma2_kmp) *(
				       gamma1_kpm*p_12p*p_11m + gamma1_kmp*p_11p*p_12m +
				       gamma1_pmk*p_12m*p_11k + gamma1_pkm*p_11m*p_12k +
				       gamma1_mkp*p_12k*p_11p + gamma1_mpk*p_11k*p_12p);

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  case _B11_:
    class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,k,&p_11k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,p,&p_11p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,m,&p_11m,errmsg),
	       errmsg,
	       pnl->error_message);
    

    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,k,&p_12k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,p,&p_12p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,m,&p_12m,errmsg),
	       errmsg,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,k,&p_22k,errmsg),
	       errmsg,
	       pnl->error_message);

    
    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,errmsg),
	       errmsg,
	       pnl->error_message);
    

    trg_gamma_121(p,k,m,&gamma1_pkm);
    trg_gamma_121(p,m,k,&gamma1_pmk);
    trg_gamma_121(m,p,k,&gamma1_mpk);
    trg_gamma_121(m,k,p,&gamma1_mkp);


    *result= gamma2_kpm *( gamma2_kpm*p_12p*p_12m + gamma1_pmk*p_12m*p_12k +
			   gamma1_pkm*p_11m*p_22k + gamma1_mkp*p_22k*p_11p +
			   gamma1_mpk*p_12k*p_12p)
      + gamma2_kmp *( gamma2_kmp*p_12m*p_12p + gamma1_mpk*p_12p*p_12k +
		      gamma1_mkp*p_11p*p_22k + gamma1_pkm*p_22k*p_11m +
		      gamma1_pmk*p_12k*p_12m);

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  case _B12_:


    class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,k,&p_11k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,p,&p_11p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,m,&p_11m,errmsg),
	       errmsg,
	       pnl->error_message);
    

    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,k,&p_12k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,p,&p_12p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,m,&p_12m,errmsg),
	       errmsg,
	       pnl->error_message);

    
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,p,&p_22p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,m,&p_22m,errmsg),
	       errmsg,
	       pnl->error_message);
    
   
    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,errmsg),
	       errmsg,
	       pnl->error_message);
       
    
    trg_gamma_121(k,p,m,&gamma1_kpm);
    trg_gamma_121(k,m,p,&gamma1_kmp);
    trg_gamma_121(p,k,m,&gamma1_pkm);
    trg_gamma_121(p,m,k,&gamma1_pmk);
    trg_gamma_121(m,p,k,&gamma1_mpk);
    trg_gamma_121(m,k,p,&gamma1_mkp);

    *result= gamma2_kpm *( gamma1_kpm*p_22p*p_11m + gamma1_kmp*p_12p*p_12m +
			   gamma2_pmk*p_12m*p_12k + gamma1_mkp*p_12k*p_12p +
			   gamma1_mpk*p_11k*p_22p)
      + gamma2_kmp *( gamma1_kmp*p_22m*p_11p + gamma1_kpm*p_12m*p_12p +
		      gamma2_mpk*p_12p*p_12k + gamma1_pkm*p_12k*p_12m +
		      gamma1_pmk*p_11k*p_22m);

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  case _B21_:
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,k,&p_12k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,p,&p_12p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,m,&p_12m,errmsg),
	       errmsg,
	       pnl->error_message);

    
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,k,&p_22k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,p,&p_22p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,m,&p_22m,errmsg),
	       errmsg,
	       pnl->error_message);
    
   
    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,errmsg),
	       errmsg,
	       pnl->error_message);
  
    
    trg_gamma_121(k,p,m,&gamma1_kpm);
    trg_gamma_121(k,m,p,&gamma1_kmp);


    *result= gamma2_kpm *( gamma1_kpm*p_22p*p_12m + gamma1_kmp*p_12p*p_22m + 
			   gamma2_pmk*p_22m*p_12k + gamma2_mkp*p_12k*p_22p)
      + gamma2_kmp *( gamma1_kmp*p_22m*p_12p + gamma1_kpm*p_12m*p_22p + 
		      gamma2_mpk*p_22p*p_12k + gamma2_pkm*p_12k*p_22m);

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  case _B22_:
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,k,&p_12k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,p,&p_12p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,m,&p_12m,errmsg),
	       errmsg,
	       pnl->error_message);

    
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,k,&p_22k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,p,&p_22p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,m,&p_22m,errmsg),
	       errmsg,
	       pnl->error_message);
    
   
    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,errmsg),
	       errmsg,
	       pnl->error_message);
  
    
    trg_gamma_121(p,k,m,&gamma1_pkm);
    trg_gamma_121(p,m,k,&gamma1_pmk);
    trg_gamma_121(m,p,k,&gamma1_mpk);
    trg_gamma_121(m,k,p,&gamma1_mkp);

    *result= gamma2_kpm *( gamma2_kpm*p_12p*p_22m + gamma1_pmk*p_22m*p_12k +
			   gamma1_pkm*p_12m*p_22k + gamma2_mkp*p_22k*p_12p)
      + gamma2_kmp *( gamma2_kmp*p_12m*p_22p + gamma1_mpk*p_22p*p_12k +
		      gamma1_mkp*p_12p*p_22k + gamma2_pkm*p_22k*p_12m);

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  case _B3_:
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,k,&p_22k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,p,&p_22p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,m,&p_22m,errmsg),
	       errmsg,
	       pnl->error_message);
    
   
    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,errmsg),
	       errmsg,
	       pnl->error_message);
  
    
    *result= gamma2_kpm *( gamma2_kpm*p_22p*p_22m + gamma2_pmk*p_22m*p_22k + 
			   gamma2_mkp*p_22k*p_22p )
      + gamma2_kmp *( gamma2_kmp*p_22m*p_22p + gamma2_mpk*p_22p*p_22k + 
		      gamma2_pkm*p_22k*p_22m );

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  default:
    sprintf(pnl->error_message,"%s(L:%d): non valid argument in integrals A of I, %d is not defined\n",__func__,__LINE__,name);
    return _FAILURE_;
    break;
  }

}

/******************
 *
 * Argument of the integral A's when mode is set to 1, i.e. for a
 * one-loop computation. The only difference with the above function
 * is that every non linear spectrum involved in the calculation is
 * replaced by the linear one extracted from the CLASS code.
 *
 ******************/


int trg_A_arg_one_loop(
		       struct spectra_nl * pnl,
		       enum name_A name, 
		       double k, 
		       double p, 
		       double m, 
		       int index_eta,
		       int index_k, /**< used only for testing, could be supressed */
		       double * result, 
		       char * errmsg){

  double p_11k,p_11m,p_11p;
  double p_12k,p_12m,p_12p;
  double p_22k,p_22m,p_22p;

  double gamma1_kpm,gamma1_pkm,gamma1_mkp,gamma1_mpk,gamma1_pmk,gamma1_kmp;
  double gamma2_kpm,gamma2_pkm,gamma2_mkp,gamma2_mpk,gamma2_pmk,gamma2_kmp;

  /**********
   * Under and above the two chosen cut-off k_max and k_min, force the
   * result to 0.
   * m can only be under k_min, while p might be under
   * k_min and above k_max.
   **********/

  if (m < pnl->k[0]) {
    *result=0.;
    return _SUCCESS_;
  }

  if ((p<pnl->k[0]) || (p>pnl->k[pnl->k_size-1])) {
    *result=0.;
    return _SUCCESS_;
  }

  /**********
   * Every function needed for the argument of a certain integral is
   * only called inside the switch case instance, in order to lower the
   * execution time - this way the lecture of this function is somehow
   * harder, but the argument is called in the very end of each case
   * bloc. Moreover, the call of functions is in this order : P_ab, the
   * gamma2, then finally gamma1.
   **********/

  switch(name){
  case _A0_:
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,k,&p_22k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,p,&p_22p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,m,&p_22m,errmsg),
	       errmsg,
	       pnl->error_message);
      
    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,errmsg),
	       errmsg,
	       pnl->error_message);

    trg_gamma_121(k,p,m,&gamma1_kpm);
    trg_gamma_121(k,m,p,&gamma1_kmp);
   
   
    *result= gamma1_kpm *( gamma2_kpm*p_22p*p_22m + gamma2_pmk*p_22m*p_22k + 
			   gamma2_mkp*p_22k*p_22p )
      + gamma1_kmp *( gamma2_kmp*p_22m*p_22p + gamma2_mpk*p_22p*p_22k + 
		      gamma2_pkm*p_22k*p_22m );

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  case _A11_:
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,k,&p_12k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,p,&p_12p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,m,&p_12m,errmsg),
	       errmsg,
	       pnl->error_message);

    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,k,&p_22k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,p,&p_22p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,m,&p_22m,errmsg),
	       errmsg,
	       pnl->error_message);

    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,errmsg),
	       errmsg,
	       pnl->error_message);

    trg_gamma_121(k,p,m,&gamma1_kpm);
    trg_gamma_121(k,m,p,&gamma1_kmp);

       


    *result= gamma1_kpm *( gamma1_kpm*p_22p*p_12m + gamma1_kmp*p_12p*p_22m + 
			   gamma2_pmk*p_22m*p_12k + gamma2_mkp*p_12k*p_22p)
      + gamma1_kmp *( gamma1_kmp*p_22m*p_12p + gamma1_kpm*p_12m*p_22p + 
		      gamma2_mpk*p_22p*p_12k + gamma2_pkm*p_12k*p_22m);

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  case _A12_:
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,k,&p_12k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,p,&p_12p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,m,&p_12m,errmsg),
	       errmsg,
	       pnl->error_message);

    
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,k,&p_22k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,p,&p_22p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,m,&p_22m,errmsg),
	       errmsg,
	       pnl->error_message);
    
   
    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,errmsg),
	       errmsg,
	       pnl->error_message);
  
    
    trg_gamma_121(k,p,m,&gamma1_kpm);
    trg_gamma_121(k,m,p,&gamma1_kmp);
    trg_gamma_121(p,k,m,&gamma1_pkm);
    trg_gamma_121(p,m,k,&gamma1_pmk);
    trg_gamma_121(m,p,k,&gamma1_mpk);
    trg_gamma_121(m,k,p,&gamma1_mkp);

    *result= gamma1_kpm *( gamma2_kpm*p_12p*p_22m + gamma1_pmk*p_22m*p_12k +
			   gamma1_pkm*p_12m*p_22k + gamma2_mkp*p_22k*p_12p)
      + gamma1_kmp *( gamma2_kmp*p_12m*p_22p + gamma1_mpk*p_22p*p_12k +
		      gamma1_mkp*p_12p*p_22k + gamma2_pkm*p_22k*p_12m);

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  case _A13_:


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,k,&p_12k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,p,&p_12p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,m,&p_12m,errmsg),
	       errmsg,
	       pnl->error_message);

    
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,k,&p_22k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,p,&p_22p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,m,&p_22m,errmsg),
	       errmsg,
	       pnl->error_message);
    
   
    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,errmsg),
	       errmsg,
	       pnl->error_message);
  
    
    trg_gamma_121(k,p,m,&gamma1_kpm);
    trg_gamma_121(k,m,p,&gamma1_kmp);
    trg_gamma_121(p,k,m,&gamma1_pkm);
    trg_gamma_121(p,m,k,&gamma1_pmk);
    trg_gamma_121(m,p,k,&gamma1_mpk);
    trg_gamma_121(m,k,p,&gamma1_mkp);

    *result= gamma1_kpm *( gamma2_kpm*p_22p*p_12m + gamma2_pmk*p_12m*p_22k +
			   gamma1_mkp*p_22k*p_12p + gamma1_mpk*p_12k*p_22p)
      + gamma1_kmp *( gamma2_kmp*p_22m*p_12p + gamma2_mpk*p_12p*p_22k +
		      gamma1_pkm*p_22k*p_12m + gamma1_pmk*p_12k*p_22m);

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  case _A21_:
    

    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,index_eta,p,&p_11p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,index_eta,m,&p_11m,errmsg),
	       errmsg,
	       pnl->error_message);
    

    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,k,&p_12k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,p,&p_12p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,m,&p_12m,errmsg),
	       errmsg,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,k,&p_22k,errmsg),
	       errmsg,
	       pnl->error_message);

    
    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,errmsg),
	       errmsg,
	       pnl->error_message);
    
    
    trg_gamma_121(k,p,m,&gamma1_kpm);
    trg_gamma_121(k,m,p,&gamma1_kmp);
    trg_gamma_121(p,k,m,&gamma1_pkm);
    trg_gamma_121(p,m,k,&gamma1_pmk);
    trg_gamma_121(m,p,k,&gamma1_mpk);
    trg_gamma_121(m,k,p,&gamma1_mkp);

    *result= gamma1_kpm *( gamma2_kpm*p_12p*p_12m + gamma1_pmk*p_12m*p_12k +
			   gamma1_pkm*p_11m*p_22k + gamma1_mkp*p_22k*p_11p +
			   gamma1_mpk*p_12k*p_12p)
      + gamma1_kmp *( gamma2_kmp*p_12m*p_12p + gamma1_mpk*p_12p*p_12k +
		      gamma1_mkp*p_11p*p_22k + gamma1_pkm*p_22k*p_11m +
		      gamma1_pmk*p_12k*p_12m);

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  case _A22_:
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,index_eta,k,&p_11k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,index_eta,p,&p_11p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,index_eta,m,&p_11m,errmsg),
	       errmsg,
	       pnl->error_message);
    

    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,k,&p_12k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,p,&p_12p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,m,&p_12m,errmsg),
	       errmsg,
	       pnl->error_message);

    
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,p,&p_22p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,m,&p_22m,errmsg),
	       errmsg,
	       pnl->error_message);
      
   
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,errmsg),
	       errmsg,
	       pnl->error_message);
  
    
    trg_gamma_121(k,p,m,&gamma1_kpm);
    trg_gamma_121(k,m,p,&gamma1_kmp);
    trg_gamma_121(p,k,m,&gamma1_pkm);
    trg_gamma_121(p,m,k,&gamma1_pmk);
    trg_gamma_121(m,p,k,&gamma1_mpk);
    trg_gamma_121(m,k,p,&gamma1_mkp);


    *result= gamma1_kpm *( gamma1_kpm*p_22p*p_11m + gamma1_kmp*p_12p*p_12m +
			   gamma2_pmk*p_12m*p_12k + gamma1_mkp*p_12k*p_12p +
			   gamma1_mpk*p_11k*p_22p)
      + gamma1_kmp *( gamma1_kmp*p_22m*p_11p + gamma1_kpm*p_12m*p_12p +
		      gamma2_mpk*p_12p*p_12k + gamma1_pkm*p_12k*p_12m +
		      gamma1_pmk*p_11k*p_22m);

    *result *= m*p/2./pow(2*_PI_,3);    

    return _SUCCESS_;
    break;

    /****************************************/

  case _A23_:
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,index_eta,k,&p_11k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,index_eta,p,&p_11p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,index_eta,m,&p_11m,errmsg),
	       errmsg,
	       pnl->error_message);
    

    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,k,&p_12k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,p,&p_12p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,m,&p_12m,errmsg),
	       errmsg,
	       pnl->error_message);

    
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,p,&p_22p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,m,&p_22m,errmsg),
	       errmsg,
	       pnl->error_message);
      
   
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,errmsg),
	       errmsg,
	       pnl->error_message);
  
    
    trg_gamma_121(k,p,m,&gamma1_kpm);
    trg_gamma_121(k,m,p,&gamma1_kmp);
    trg_gamma_121(p,k,m,&gamma1_pkm);
    trg_gamma_121(p,m,k,&gamma1_pmk);
    trg_gamma_121(m,p,k,&gamma1_mpk);
    trg_gamma_121(m,k,p,&gamma1_mkp);


    *result= gamma1_kpm *( gamma1_kpm*p_12p*p_12m + gamma1_kmp*p_11p*p_22m +
			   gamma1_pmk*p_22m*p_11k + gamma1_pkm*p_12m*p_12k +
			   gamma2_mkp*p_12k*p_12p)
      + gamma1_kmp *( gamma1_kmp*p_12m*p_12p + gamma1_kpm*p_11m*p_22p +
		      gamma1_mpk*p_22p*p_11k + gamma1_mkp*p_12p*p_12k +
		      gamma2_pkm*p_12k*p_12m);

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  case _A3_:
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,index_eta,k,&p_11k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,index_eta,p,&p_11p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,index_eta,m,&p_11m,errmsg),
	       errmsg,
	       pnl->error_message);
    

    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,k,&p_12k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,p,&p_12p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,m,&p_12m,errmsg),
	       errmsg,
	       pnl->error_message);

    
    trg_gamma_121(k,p,m,&gamma1_kpm);
    trg_gamma_121(k,m,p,&gamma1_kmp);
    trg_gamma_121(p,k,m,&gamma1_pkm);
    trg_gamma_121(p,m,k,&gamma1_pmk);
    trg_gamma_121(m,p,k,&gamma1_mpk);
    trg_gamma_121(m,k,p,&gamma1_mkp);


    *result= (gamma1_kpm+gamma1_kmp) *(
				       gamma1_kpm*p_12p*p_11m + gamma1_kmp*p_11p*p_12m +
				       gamma1_pmk*p_12m*p_11k + gamma1_pkm*p_11m*p_12k +
				       gamma1_mkp*p_12k*p_11p + gamma1_mpk*p_11k*p_12p);

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  case _B0_:
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,index_eta,k,&p_11k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,index_eta,p,&p_11p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,index_eta,m,&p_11m,errmsg),
	       errmsg,
	       pnl->error_message);
    

    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,k,&p_12k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,p,&p_12p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,m,&p_12m,errmsg),
	       errmsg,
	       pnl->error_message);

    
    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,errmsg),
	       errmsg,
	       pnl->error_message);
    
    
    trg_gamma_121(k,p,m,&gamma1_kpm);
    trg_gamma_121(k,m,p,&gamma1_kmp);
    trg_gamma_121(p,k,m,&gamma1_pkm);
    trg_gamma_121(p,m,k,&gamma1_pmk);
    trg_gamma_121(m,p,k,&gamma1_mpk);
    trg_gamma_121(m,k,p,&gamma1_mkp);

    *result= (gamma2_kpm+gamma2_kmp) *(
				       gamma1_kpm*p_12p*p_11m + gamma1_kmp*p_11p*p_12m +
				       gamma1_pmk*p_12m*p_11k + gamma1_pkm*p_11m*p_12k +
				       gamma1_mkp*p_12k*p_11p + gamma1_mpk*p_11k*p_12p);

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  case _B11_:
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,index_eta,k,&p_11k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,index_eta,p,&p_11p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,index_eta,m,&p_11m,errmsg),
	       errmsg,
	       pnl->error_message);
    

    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,k,&p_12k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,p,&p_12p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,m,&p_12m,errmsg),
	       errmsg,
	       pnl->error_message);

    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,k,&p_22k,errmsg),
	       errmsg,
	       pnl->error_message);
    
    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,errmsg),
	       errmsg,
	       pnl->error_message);
    

    trg_gamma_121(p,k,m,&gamma1_pkm);
    trg_gamma_121(p,m,k,&gamma1_pmk);
    trg_gamma_121(m,p,k,&gamma1_mpk);
    trg_gamma_121(m,k,p,&gamma1_mkp);


    *result= gamma2_kpm *( gamma2_kpm*p_12p*p_12m + gamma1_pmk*p_12m*p_12k +
			   gamma1_pkm*p_11m*p_22k + gamma1_mkp*p_22k*p_11p +
			   gamma1_mpk*p_12k*p_12p)
      + gamma2_kmp *( gamma2_kmp*p_12m*p_12p + gamma1_mpk*p_12p*p_12k +
		      gamma1_mkp*p_11p*p_22k + gamma1_pkm*p_22k*p_11m +
		      gamma1_pmk*p_12k*p_12m);

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  case _B12_:


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,index_eta,k,&p_11k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,index_eta,p,&p_11p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,index_eta,m,&p_11m,errmsg),
	       errmsg,
	       pnl->error_message);
    

    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,k,&p_12k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,p,&p_12p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,m,&p_12m,errmsg),
	       errmsg,
	       pnl->error_message);

    
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,p,&p_22p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,m,&p_22m,errmsg),
	       errmsg,
	       pnl->error_message);
    
   
    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,errmsg),
	       errmsg,
	       pnl->error_message);
       
    
    trg_gamma_121(k,p,m,&gamma1_kpm);
    trg_gamma_121(k,m,p,&gamma1_kmp);
    trg_gamma_121(p,k,m,&gamma1_pkm);
    trg_gamma_121(p,m,k,&gamma1_pmk);
    trg_gamma_121(m,p,k,&gamma1_mpk);
    trg_gamma_121(m,k,p,&gamma1_mkp);

    *result= gamma2_kpm *( gamma1_kpm*p_22p*p_11m + gamma1_kmp*p_12p*p_12m +
			   gamma2_pmk*p_12m*p_12k + gamma1_mkp*p_12k*p_12p +
			   gamma1_mpk*p_11k*p_22p)
      + gamma2_kmp *( gamma1_kmp*p_22m*p_11p + gamma1_kpm*p_12m*p_12p +
		      gamma2_mpk*p_12p*p_12k + gamma1_pkm*p_12k*p_12m +
		      gamma1_pmk*p_11k*p_22m);

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  case _B21_:
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,k,&p_12k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,p,&p_12p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,m,&p_12m,errmsg),
	       errmsg,
	       pnl->error_message);

    
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,k,&p_22k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,p,&p_22p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,m,&p_22m,errmsg),
	       errmsg,
	       pnl->error_message);
    
   
    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,errmsg),
	       errmsg,
	       pnl->error_message);
  
    
    trg_gamma_121(k,p,m,&gamma1_kpm);
    trg_gamma_121(k,m,p,&gamma1_kmp);


    *result= gamma2_kpm *( gamma1_kpm*p_22p*p_12m + gamma1_kmp*p_12p*p_22m + 
			   gamma2_pmk*p_22m*p_12k + gamma2_mkp*p_12k*p_22p)
      + gamma2_kmp *( gamma1_kmp*p_22m*p_12p + gamma1_kpm*p_12m*p_22p + 
		      gamma2_mpk*p_22p*p_12k + gamma2_pkm*p_12k*p_22m);

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  case _B22_:
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,k,&p_12k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,p,&p_12p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,index_eta,m,&p_12m,errmsg),
	       errmsg,
	       pnl->error_message);

    
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,k,&p_22k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,p,&p_22p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,m,&p_22m,errmsg),
	       errmsg,
	       pnl->error_message);
    
   
    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,errmsg),
	       errmsg,
	       pnl->error_message);
  
    
    trg_gamma_121(p,k,m,&gamma1_pkm);
    trg_gamma_121(p,m,k,&gamma1_pmk);
    trg_gamma_121(m,p,k,&gamma1_mpk);
    trg_gamma_121(m,k,p,&gamma1_mkp);

    *result= gamma2_kpm *( gamma2_kpm*p_12p*p_22m + gamma1_pmk*p_22m*p_12k +
			   gamma1_pkm*p_12m*p_22k + gamma2_mkp*p_22k*p_12p)
      + gamma2_kmp *( gamma2_kmp*p_12m*p_22p + gamma1_mpk*p_22p*p_12k +
		      gamma1_mkp*p_12p*p_22k + gamma2_pkm*p_22k*p_12m);

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  case _B3_:
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,k,&p_22k,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,p,&p_22p,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,index_eta,m,&p_22m,errmsg),
	       errmsg,
	       pnl->error_message);
    
   
    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,errmsg),
	       errmsg,
	       pnl->error_message);
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,errmsg),
	       errmsg,
	       pnl->error_message);
  
    
    *result= gamma2_kpm *( gamma2_kpm*p_22p*p_22m + gamma2_pmk*p_22m*p_22k + 
			   gamma2_mkp*p_22k*p_22p )
      + gamma2_kmp *( gamma2_kmp*p_22m*p_22p + gamma2_mpk*p_22p*p_22k + 
		      gamma2_pkm*p_22k*p_22m );

    *result *= m*p/2./pow(2*_PI_,3);

    return _SUCCESS_;
    break;

    /****************************************/

  default:
    sprintf(pnl->error_message,"%s(L:%d): non valid argument in integrals A of I, %d is not defined\n",__func__,__LINE__,name);
    return _FAILURE_;
    break;
  }

}




/********************
 *
 * Integration with the Simpson's method (over y, and simple trapeze for x)
 *
 ********************/

int trg_integrate_xy_at_eta(
			    struct background * pba,
			    struct primordial * ppm,
			    struct spectra * psp,
			    struct spectra_nl * pnl,
			    enum name_A name,
			    int index_eta,
			    double * result,
			    char *errmsg
			    ){
  
  int index_k;
  double k,k_min,k_max;

  int x_size,index_x;
  double x;
  double * xx;
  double logstepx;
    
  int y_size,index_y;
  double y;
  double * yy;
  double logstepy;

  double * h_up;
  double * h_do;
  double * v_le;
  double * v_ri;
  
  int il;

  double * partial_sum;
  double * partial_area;
  double sum,area,max;
  double increment_sum,increment_area;

  /**********
   * The following approximation consists in taking into account the
   * fact that for small k's, the influence of non linearity is
   * small. Hence, one could take the linear evolution for P's for
   * sufficiently small k's. Here, arbitrarily, one takes values of
   * index_k before 10 to follow linear regime, i.e. putting to zero
   * their A's.
   **********/

  /**********
   * This approximation is here to calculate more precisely the first
   * k points on x-y domains. The other solution (cleaner, but more
   * time consuming) would consist in modifying the lower bound in the
   * table_array routines by chosing an arbitrarily small cut-off for
   * P's (well under the k_min).
   **********/

  k_min=pnl->k[0];
  k_max=pnl->k[pnl->k_size-1];

  for(index_k=0; index_k<pnl->k_size; index_k++){

    k=pnl->k[index_k];

    logstepx=min(1.1,1+0.01/pow(k,1)); /***> PRECISION PARAMETER <***/

    if(logstepx<1.0035) 
      logstepx=1.0037;                 /***> PRECISION PARAMETER <***/

    logstepy=logstepx;

    /**********
     * Linear theory is applied below the scale defined by k_L, chosen
     * in trg.c
     **********/

    if(index_k<pnl->index_k_L){

      result[index_k+pnl->k_size*index_eta]=0.;

    }


     /***                                                          ***/
    /**** TO CHECK : DO WITHOUT : MAYBE COMPENSATIONS THERE ALSO ! ****/
     /***                                                          ***/ 



    /**********
     * Near the highest k's computed, the necessary cut-off in the
     * computation of the integral introduces bad values : the
     * different constituant do not compensate anymore.
     **********/

    else if(index_k>pnl->k_size-5){
    
      result[index_k+pnl->k_size*index_eta]=
    	((result[pnl->k_size-5+pnl->k_size*index_eta]-result[pnl->k_size-6+pnl->k_size*index_eta])/
    	 (pnl->k[pnl->k_size-5]-pnl->k[pnl->k_size-6])) * (pnl->k[index_k]-pnl->k[pnl->k_size-5]) +
    	result[pnl->k_size-5+pnl->k_size*index_eta];
      
    }

    /**********
     * Normal computation for all the remaining points
     **********/
    
    else{


      x_size = (int)(log(2.*k_max/k)/log(logstepx)) + 2;

      class_calloc(xx,x_size,sizeof(double),pnl->error_message);

      class_calloc(h_up,x_size,sizeof(double),pnl->error_message);
      class_calloc(h_do,x_size,sizeof(double),pnl->error_message);

      index_x = 0; 

      do {
	
	class_test(index_x >= x_size,
		   pnl->error_message," ");

	xx[index_x] = k/sqrt(2.)*pow(logstepx,index_x); /* set value of x */

	if (xx[index_x] >= k_max*sqrt(2.)) xx[index_x] = k_max*sqrt(2.); /* correct the last value */

	index_x ++;

      } while (xx[index_x-1] < k_max*sqrt(2.));

      x_size = index_x;

      y_size = (int)(log(2.)/log(logstepy)) + 2;

      class_calloc(yy,y_size,sizeof(double),pnl->error_message);

      class_calloc(v_le,y_size,sizeof(double),pnl->error_message);
      class_calloc(v_ri,y_size,sizeof(double),pnl->error_message);
 
      class_calloc(partial_sum,y_size-1,sizeof(double),pnl->error_message);
      class_calloc(partial_area,y_size-1,sizeof(double),pnl->error_message);

      index_y = 0;

      do {

	class_test(index_y >= y_size,
		   pnl->error_message," ");

	yy[index_y] = k*sqrt(2.) - k/sqrt(2.)*pow(logstepy,index_y); /* set value of y */

	if (yy[index_y] < 0.) yy[index_y] = 0.; /* correct the last value */

	index_y ++;
	
      } while (yy[index_y-1] > 0.);

      y_size = index_y;

      h_do[0]=0.;
      v_ri[0]=h_do[0];

      for (index_x=1; index_x < x_size; index_x ++) {

	x=xx[index_x];
	y=yy[0];

	if (x <= sqrt(2.)*k_max-y) {
	  if(pnl->mode==1){
	    class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&h_do[index_x],errmsg),
		       errmsg,
		       pnl->error_message);}
	  else if(pnl->mode==2){
	    class_call(trg_A_arg(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&h_do[index_x],errmsg),
		       errmsg,
		       pnl->error_message);}
	}
	else {
	  h_do[index_x]=0.;
	}
	
      }

      for (index_y=1; index_y < y_size; index_y ++) {

	x=xx[0];
	y=yy[index_y];

	if(pnl->mode==1){
	  class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&v_ri[index_y],errmsg),
		     errmsg,
		     pnl->error_message);}

	else if(pnl->mode==2){
	  class_call(trg_A_arg(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&v_ri[index_y],errmsg),
		     errmsg,
		     pnl->error_message);}
      }

      sum = 0.;
      area = 0.;
      max = 0.;

      /********************* loop over L-shaped regions **********************/

      for (il=0; il < y_size-1; il++) {

	/* move previous bottom-line to up-line, and previous right-line to left-line 
	   (remember that some point may have not been calculated, 
	   in this case the end of the lines is filled with zeros */

	for (index_x=il; index_x < x_size; index_x ++)
	  h_up[index_x] = h_do[index_x];

	for (index_y=il; index_y < y_size; index_y ++)
	  v_le[index_y] = v_ri[index_y];

	/* one new point on the diagonal, integral of cell on diagonal */

	x=xx[il+1];
	y=yy[il+1];

	if(pnl->mode==1){
	  class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&h_do[il+1],errmsg),
		     errmsg,
		     pnl->error_message);}
	else if(pnl->mode==2){
	  class_call(trg_A_arg(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&h_do[il+1],errmsg),
		     errmsg,
		     pnl->error_message);}

	v_ri[il+1]= h_do[il+1];

	increment_sum = (xx[il+1]-xx[il])*(yy[il]-yy[il+1])*0.25*(h_up[il]+h_up[il+1]+v_le[il+1]+v_ri[il+1]);
	increment_area = (xx[il+1]-xx[il])*(yy[il]-yy[il+1]);

	partial_sum[il] = increment_sum;
	partial_area[il] = increment_area;

	/**************** new points on the horizontal and diagonal, within the square ******************/

	for (index_x=il+1; index_x < y_size-1; index_x ++) {

	  x=xx[index_x+1];
	  y=yy[il];
	    
	  if (x <= sqrt(2)*k_max-y) {
	    if(pnl->mode==1){
	      class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&h_up[index_x+1],errmsg),
			 errmsg,
			 pnl->error_message);}
	    else if(pnl->mode==2){
	      class_call(trg_A_arg(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&h_up[index_x+1],errmsg),
			 errmsg,
			 pnl->error_message);}
	  } 
	  else {
	    h_up[index_x+1]=0.;
	  }
	    
	  

	  /* the point h_do[index_x+1] is new; calculate */

	  x=xx[index_x+1];
	  y=yy[il+1];

	  if (x <= sqrt(2)*k_max-y) {
	    if(pnl->mode==1){
	      class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&h_do[index_x+1],errmsg),
			 errmsg,
			 pnl->error_message);}
	    else if(pnl->mode==2){
	      class_call(trg_A_arg(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&h_do[index_x+1],errmsg),
			 errmsg,
			 pnl->error_message);}
	  } 
	  else {
	    h_do[index_x+1]=0.;
	  }

	  /* the point v_le[index_x+1] may have not been calculated at the previous stage; check and calculate */

	  x=xx[il];
	  y=yy[index_x+1];
	    
	  if(pnl->mode==1){
	    class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&v_le[index_x+1],errmsg),
		       errmsg,
		       pnl->error_message);}
	  else if(pnl->mode==2){
	    class_call(trg_A_arg(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&v_le[index_x+1],errmsg),
		       errmsg,
		       pnl->error_message);}
	  

	  /* the point v_ri[index_x+1] is new; calculate */

	  x=xx[il+1];
	  y=yy[index_x+1];

	  if(pnl->mode==1){
	    class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&v_ri[index_x+1],errmsg),
		       errmsg,
		       pnl->error_message);}
	  else if(pnl->mode==2){
	    class_call(trg_A_arg(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&v_ri[index_x+1],errmsg),
		       errmsg,
		       pnl->error_message);}

	  /* now integrate on the two new cells */

	  increment_sum = (xx[il+1]-xx[il])*(yy[index_x]-yy[index_x+1])*0.25*
	    (v_le[index_x]+v_le[index_x+1]+v_ri[index_x]+v_ri[index_x+1])
	    + (xx[index_x+1]-xx[index_x])*(yy[il]-yy[il+1])*0.25*
	    (h_up[index_x]+h_up[index_x+1]+h_do[index_x]+h_do[index_x+1]);

	  increment_area = (xx[il+1]-xx[il])*(yy[index_x]-yy[index_x+1]) 
	    + (xx[index_x+1]-xx[index_x])*(yy[il]-yy[il+1]);

	  partial_sum[il] += increment_sum;
	  partial_area[il] += increment_area;

	}

	/***************** new points on the horizontal, beyond the square *******************/

	for (index_x=y_size-1; index_x < x_size-1; index_x ++) {

	  /* the point h_up[index_x+1] may have not been calculated at the previous stage; check and calculate */
	  
	  x=xx[index_x+1];
	  y=yy[il];
	      
	  if (x <= sqrt(2)*k_max-y) {
	    if(pnl->mode==1){
	      class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&h_up[index_x+1],errmsg),
			 errmsg,
			 pnl->error_message);}
	    else if(pnl->mode==2){
	      class_call(trg_A_arg(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&h_up[index_x+1],errmsg),
			 errmsg,
			 pnl->error_message);}
	  } 
	  else {
	    h_up[index_x+1]=0.;
	  }
	  
	  /* the point h_do[index_x+1] is new; calculate */

	  x=xx[index_x+1];
	  y=yy[il+1];

	  if (x <= sqrt(2)*k_max-y) {
	    if(pnl->mode==1){
	      class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&h_do[index_x+1],errmsg),
			 errmsg,
			 pnl->error_message);}
	    else if(pnl->mode==2){
	      class_call(trg_A_arg(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&h_do[index_x+1],errmsg),
			 errmsg,
			 pnl->error_message);}
	  } 
	  else {
	    h_do[index_x+1]=0.;
	  }

	  /* now integrate on the new cell */
	  
	  increment_sum = (xx[index_x+1]-xx[index_x])*(yy[il]-yy[il+1])*0.25*
	    (h_up[index_x]+h_up[index_x+1]+h_do[index_x]+h_do[index_x+1]);
	  
	  increment_area = (xx[index_x+1]-xx[index_x])*(yy[il]-yy[il+1]);
	  
	  partial_sum[il] += increment_sum;
	  partial_area[il] += increment_area;
	  
	}
	

	/* update the total sum with the new L-shaped region */

	sum += partial_sum[il];
	area += partial_area[il];
	
      }

      result[index_k+pnl->k_size*index_eta]=sum;

      free(xx);
      free(h_up);
      free(h_do);
      free(yy);
      free(v_le);
      free(v_ri);
      free(partial_sum);
      free(partial_area);
    
    }

  }

  return _SUCCESS_;
  
}

/*******************************
 *
 * Summary : taking the power spectrum of matter, matter-velocity and
 *           velocity, at a chosen time a_ini when matter dominated
 *           already but before all modes become non linear, it
 *           computes the fully non-linear evolution of these 3
 *           quantities for a range of eta=log(a/a_ini) up to today.
 *
 ******************************/

/******************************
 *
 *	       It requires the definition of many quantities (14+14)
 *	       built on the Bispectra of matter and velocity and on
 *	       the spectra as well, in order to put the equations in a
 *	       more numerically friendly shape, as explained in annex
 *	       B of the Pietroni article (astro-ph 0806.0971v3).
 *
 *********************************/

int trg_init (
	      struct precision  * ppr, /* input */
	      struct background * pba, /**< structure containing all background evolution */
	      struct perturbs   * ppt, 
	      struct primordial * ppm,
	      struct spectra    * psp, /**< structure containing list of k, z and P(k,z) values) */
	      struct spectra_nl * pnl /* output */ 
	      ) {

  FILE *nl_spectra;

  time_t time_1,time_2;

  int index_k;
  int index_eta;

  int index_ic=ppt->index_ic; /* for adiabatic modes = 0 */

  double * temp_k;
  double temp; 

  double k_max;
  double k_L;
  double k_min;
  
  int index_name,name_size;
  double ** AA; /* Definition of all the A's integral in one vector of
		   matrix */

  /*
    Variables of calculation, k and eta, and initialization
  */

  double a_ini;
  double z_ini;
  double eta_max;
  double * Omega_m, * H, *H_prime;

  double eta_step; /* step for integrations */

  double exp_eta;
  double fourpi_over_k;

  int index;  
  int index_plus; /* intermediate variables for time integration */
  
  double Omega_11,Omega_12; /* Definition of the matrix Omega that
			       mixes terms together, two are time
			       indepedant and two dependant --- SO
			       FAR ---*/
  double *Omega_21,*Omega_22;

  double *a0,*a11,*a12,*a13,*a21,*a22,*a23,*a3;
  double *b0,*b11,*b12,*b21,*b22,*b3; /* Definition of the variables
					 1ijk,lmn and 2 ijk,lmn that
					 replace the Bispectra
					 variables*/

  double * pvecback_nl;

  /* This code can be optionally compiled with the openmp option for parallel computation.
     Inside parallel regions, the use of the command "return" is forbidden.
     For error management, instead of "return _FAILURE_", we will set the variable below
     to "abort = _TRUE_". This will lead to a "return _FAILURE_" jus after leaving the 
     parallel region. */
  int abort;

#ifdef _OPENMP
  /* instrumentation times */
  double tstart, tstop;
#endif

  pnl->spectra_nl_verbose=1; /* 0 gives no information at all, 1 the
				minimum amount of it, 2 gives a
				percentage evolution of time
				integrator while 3 gives also the time
				spend in every parallel region */

  pnl->mode=2; /* 0 is linear evolution, 1 one loop and 2 full trg */

  if (pnl->spectra_nl_verbose > 0)
    printf("Computing non-linear spectra with TRG method\n");

  class_alloc(pvecback_nl,pba->bg_size*sizeof(double),pnl->error_message);

  /* Define initial eta and redshift */
  a_ini = 2e-2;   /* gives a starting redshift of 49 */
  pnl->z_ini = pba->a_today/a_ini - 1.;


  /* Define eta_max, where eta=log(a/a_ini) */
  eta_max = log(pba->a_today/a_ini);


  /* Define size and step for integration in eta */
  pnl->eta_size = 75; /***> PRECISION PARAMETER <***/
  pnl->eta_step = (eta_max)/(pnl->eta_size-1);
  eta_step = pnl->eta_step;

  /* at any time, a = a_ini * exp(eta) and z=exp(eta)-1*/

  
  /* Define k parameters */ 

  k_max = pnl->k_max;  /**< PRECISION PARAMETER */ /* set in test_trg.c */
  k_L   = 1.e-3;        /**< PRECISION PARAMETER */  /* scale above which we take linear theory */
  k_min = 1.e-4;      /**< PRECISION PARAMETER */  /* minimun cut-off scale  */

  /* Defining the values of k with three set of log steps */
  index_k=0;
  class_calloc(temp_k,2000,sizeof(double),pnl->error_message);
  temp_k[0]=k_min;

  while( temp_k[index_k] < 0.01 ) {
    index_k++;
    temp_k[index_k]=k_min*pow(1.2,index_k); 
  } /* To have roughly 10 points per decade until 0.01 h/Mpc */
  temp=index_k;
  temp_k[index_k]=0.01;
  while( temp_k[index_k] < 50 ) {
    index_k++;
    temp_k[index_k]=0.01*pow(1.02,(index_k-temp));
  } /* To have roughly 100 points per decade until 50 h/Mpc */
  temp=index_k;
  temp_k[index_k]=50;
  while( temp_k[index_k] < k_max ) {
    index_k++;
    temp_k[index_k]=50*pow(1.12,(index_k-temp));
  } /* To have roughly 30 points per decade afterwards */
  temp_k[index_k]=k_max;

  pnl->k_size=index_k+1;

  class_calloc(pnl->k,pnl->k_size,sizeof(double),pnl->error_message);

  /* then fill in the values of k, while determining the index below
     which we use only the linear theory */

  temp=0;

  for(index_k=0; index_k<pnl->k_size; index_k++){
    pnl->k[index_k]=temp_k[index_k];
    if( (pnl->k[index_k] > k_L*pba->h) && temp==0){
      pnl->index_k_L=index_k;
      temp++;
    }
  }
 
  free(temp_k);
 
  /* fill array of eta values, and pick two values z_1 and z_2
     resp. near 2.33 and 1.00 (for output) */

  double z_1,z_2;

  z_1=2.33;
  z_2=1.00;

  int index_eta_1,index_eta_2;

  pnl->eta = calloc(pnl->eta_size,sizeof(double));
  pnl->z   = calloc(pnl->eta_size,sizeof(double));

  int temp_z;

  temp_z=0;

  /* Fill in values of eta and z, while determining the precise index
     of z_1 and z_2 */

  for (index_eta=0; index_eta<pnl->eta_size; index_eta++) {
    pnl->eta[index_eta] = index_eta*eta_step;
    pnl->z[index_eta]   = exp(-pnl->eta[index_eta])*(pba->a_today/a_ini)-1;
    if(pnl->z[index_eta]<z_1 && temp_z==0){
      index_eta_1=index_eta;
      temp_z=1;
    }
    if(pnl->z[index_eta]<z_2 && temp_z==1){
      index_eta_2=index_eta;
      temp_z=2;
    }
    if(pnl->z[index_eta]<0) pnl->z[index_eta]=0;
  }

  if (pnl->spectra_nl_verbose > 0)
    printf(" -> starting calculation at redshift z = %2.2f\n",pnl->z[0]);

  /* definition of background values for each eta */

  class_alloc(Omega_m,pnl->eta_size*sizeof(double),pnl->error_message);
  class_alloc(H,pnl->eta_size*sizeof(double),pnl->error_message);
  class_alloc(H_prime,pnl->eta_size*sizeof(double),pnl->error_message);

  for (index_eta=0; index_eta<pnl->eta_size; index_eta++){
    class_call(background_functions_of_a(
					 pba,
					 a_ini*exp(pnl->eta[index_eta]),
					 long_info,
					 pvecback_nl
					 ),
	       pba->error_message,
	       pnl->error_message);
    
    Omega_m[index_eta] = pvecback_nl[pba->index_bg_Omega_b];
    if (pba->has_cdm == _TRUE_) {
      Omega_m[index_eta] += pvecback_nl[pba->index_bg_Omega_cdm];
    }
    
    H[index_eta] = pvecback_nl[pba->index_bg_H] * a_ini * exp(pnl->eta[index_eta]);
    H_prime[index_eta] =H[index_eta]*(1 + pvecback_nl[pba->index_bg_H_prime] / a_ini * exp(-pnl->eta[index_eta])/pvecback_nl[pba->index_bg_H]/pvecback_nl[pba->index_bg_H]);
    
  }

  /* now Omega_m, H are known at each eta */

  /* Definition of the matrix elements Omega_11,Omega_12, Omega_22,
     Omega_12 for each eta, to modify in order to take into account
     more physics with a k dependence, for instance */

  Omega_11 = 1.;
  Omega_12 = -1.;
  
  class_alloc(Omega_21,pnl->eta_size*sizeof(double),pnl->error_message);
  class_alloc(Omega_22,pnl->eta_size*sizeof(double),pnl->error_message);

  for(index_eta=0; index_eta<pnl->eta_size; index_eta++) {
      Omega_21[index_eta] = -3./2 * Omega_m[index_eta];
      Omega_22[index_eta] = 2 + H_prime[index_eta]/H[index_eta];
    }

  /* Definition of P_ab=<\phi_a \phi_b> and initialization at
     eta[0]=0, k[0]. At the first time, no previous non-gaussianity is
     considered, so the spectra are all linear. Hence P_ab=p_ab usual,
     without time factors. Moreover, P_ab=P_cd for all a b c d, since
     in the linear theory, whatever domination era it is, this
     equality holds.

     Convention is taken so that all the upcoming matrices are lines
     at constant eta and column at constant k. So the first line will
     constitute of the initial condition given by precedent parts of
     the code, and the time evolution will fill down the rest step by
     step.

     As a result of the definition in only one line (to avoid
     splitting the lines of the matrices everywhere physically),
     proper call for all these matrices is A[index_k +
     k_size*index_eta] instead of the more intuitive version
     A[index_k][index_eta].*/

  class_calloc(pnl->pk_nl,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(pnl->p_12_nl,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(pnl->p_22_nl,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);

  class_calloc(pnl->p_11,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(pnl->p_12,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(pnl->p_22,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);


  for(index_k=0; index_k<pnl->k_size; index_k++){

    class_call(trg_p11_at_k(pba,ppm,psp,pnl,0,index_ic,pnl->k[index_k],&pnl->pk_nl[index_k]),
	       pnl->error_message,
	       pnl->error_message);

    pnl->p_11[index_k]=pnl->pk_nl[index_k];
    
    class_call(trg_p11_at_k(pba,ppm,psp,pnl,0,index_ic,pnl->k[index_k],&pnl->p_12_nl[index_k]),
	       pnl->error_message,
	       pnl->error_message);
    
    pnl->p_12[index_k]=pnl->p_12_nl[index_k];

    class_call(trg_p11_at_k(pba,ppm,psp,pnl,0,index_ic,pnl->k[index_k],&pnl->p_22_nl[index_k]),
	       pnl->error_message,
	       pnl->error_message);

    pnl->p_22[index_k]=pnl->p_22_nl[index_k];

  }


  /* Initialisation and definition of second derivatives */

  class_calloc(pnl->ddpk_nl,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(pnl->ddp_12_nl,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(pnl->ddp_22_nl,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);

  class_calloc(pnl->ddp_11,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(pnl->ddp_12,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(pnl->ddp_22,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  
  class_call(trg_ddp_ab(pnl,pnl->pk_nl,0,pnl->ddpk_nl,pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);

  class_call(trg_ddp_ab(pnl,pnl->p_12_nl,0,pnl->ddp_12_nl,pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);

  class_call(trg_ddp_ab(pnl,pnl->p_22_nl,0,pnl->ddp_22_nl,pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);


  class_call(trg_ddp_ab(pnl,pnl->p_11,0,pnl->ddp_11,pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);

  class_call(trg_ddp_ab(pnl,pnl->p_12,0,pnl->ddp_12,pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);

  class_call(trg_ddp_ab(pnl,pnl->p_22,0,pnl->ddp_22,pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);


  /* Definition of I_{121,lmn} (here a0, a11,...) and I_{222,lmn}
     (here b0,b11,...) and initialization (directly with calloc
     for assuming no initial non gaussianity in the spectrum) */

  /**********
   * The convention for naming the different I_{121,lmn} in a's
   * (resp. I_{222,lmn} in b's) is the following : 
   * --> the numbers after the a (resp. b) indicate first the 
   *     number of 1 (resp. 2) in the last three indices
   * --> the possibly second number indicates, when such an
   *     information is needed, the position of the lonely index
   **********/

  /**********
   * As an example worth more than arid explanations, I_{121,212} will
   * be called a12. Indeed, the three first indices are 121, so it is an
   * 'a'. Then there is only one '1' in the last three indices, so
   * there is a '1', and finally the only '1' in these three indices
   * is in the 2nd position, so a '2'. a12. On the other hand,
   * I_{222,212} will be called b22, because there a two '2', and the
   * lonely '1' in the last indices is at the 2nd position, so '2'
   * again. Nothing cryptic !
   **********/

  class_calloc(a0 ,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(a11,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(a12,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(a13,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(a21,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(a22,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(a23,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(a3 ,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);

  class_calloc(b0 ,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(b11,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(b12,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(b21,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(b22,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(b3 ,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message); 

  /* Definition of the A's that appear in time-evolution of
     the 1's and 2's */

  name_size = _B3_+1;
  class_calloc(AA,name_size,sizeof(double*),pnl->error_message);
  for (index_name=0; index_name<name_size; index_name++) 
    class_calloc(AA[index_name],pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);


  /********************
   *
   * Integration at eta=0 with trg_integrate_at_xy().
   *
   ********************/
  double temp1;

  if (pnl->spectra_nl_verbose > 1)
    printf(" -> initialisation\n");
  
  if(pnl->mode > 0){

    /* initialize error management flag */
    abort = _FALSE_;

    /*** beginning of parallel region ***/

#pragma omp parallel							\
  shared(name_size,abort,pba,ppm,psp,pnl,AA)				\
  private(tstart,index_name,tstop)
    
    {

#ifdef _OPENMP
      tstart = omp_get_wtime();
#endif

#pragma omp for schedule (dynamic)
      for (index_name=0; index_name<name_size; index_name++) {

#pragma omp flush(abort)

	class_call_parallel(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,index_name,0,AA[index_name],pnl->error_message),
			    pnl->error_message,
			    pnl->error_message);
      }

#ifdef _OPENMP
	tstop = omp_get_wtime();
	if (pnl->spectra_nl_verbose > 2)
	  printf("In %s: time spent in parallel region (loop over names) = %e s for thread %d\n",
		 __func__,tstop-tstart,omp_get_thread_num());
#endif

    } /* end of parallel region */

    if (abort == _TRUE_) return _FAILURE_;

  }

  class_open(nl_spectra,"output/nl_ini.dat","wr",pnl->error_message);

  for(index_k=0; index_k<pnl->k_size; index_k++){ 
    fprintf(nl_spectra,"%e ",pnl->k[index_k]/pba->h);
    for (index_name=0; index_name<name_size; index_name++)
      fprintf(nl_spectra,"%e ",AA[index_name][index_k]);
    fprintf(nl_spectra,"%e %e %e\n",
	    pnl->pk_nl[index_k]*pow(pba->h,3),
	    pnl->p_12_nl[index_k]*pow(pba->h,3),
	    pnl->p_22_nl[index_k]*pow(pba->h,3));
  }
  
  fclose(nl_spectra);


  /********************
   * Calculate the time evolution with a very simple integrator
   ********************/

  class_open(nl_spectra,"output/nl_spectra.dat","wr",pnl->error_message)

    if (pnl->spectra_nl_verbose > 1){
      printf(" -> progression:\n");}

  time_1=time(NULL);
 
  for (index_eta=1; index_eta<pnl->eta_size; index_eta++){
    exp_eta=exp(pnl->eta[index_eta-1]);

    for (index_k=0; index_k<pnl->k_size; index_k++){

      /**********
       * Computation for every k of the function at the next time step
       **********/

      fourpi_over_k=4*_PI_/(pnl->k[index_k]);
      index = index_k+pnl->k_size*(index_eta-1);
      index_plus = index_k+pnl->k_size*index_eta;

     

      pnl->pk_nl[index_plus]= eta_step *(
					 -2*Omega_11*pnl->pk_nl[index]
					 -2*Omega_12*pnl->p_12_nl[index]
					 +exp_eta*4*fourpi_over_k*a22[index] )
	+ pnl->pk_nl[index];

      pnl->p_22_nl[index_plus] = eta_step *(
					    -2*Omega_22[index_eta-1]*pnl->p_22_nl[index]
					    -2*Omega_21[index_eta-1]*pnl->p_12_nl[index]
					    +exp_eta*2*fourpi_over_k*b3[index] )
	+ pnl->p_22_nl[index];

      pnl->p_12_nl[index_plus] = eta_step *(
					    -pnl->p_12_nl[index]*(Omega_11+Omega_22[index_eta-1])
					    -Omega_12*pnl->p_22_nl[index]
					    -Omega_21[index_eta-1]*pnl->pk_nl[index]
					    +exp_eta*fourpi_over_k*(2*a13[index]+b21[index]))
	+ pnl->p_12_nl[index];


      a0[index_plus]         = eta_step *(
					  -Omega_21[index_eta-1]*(a11[index]+a12[index]+a13[index])
					  -3*Omega_22[index_eta-1]*a0[index]
					  +2*exp_eta*AA[_A0_][index])
	+ a0[index];

      a11[index_plus]        = eta_step *(
					  -a11[index]*(2*Omega_22[index_eta-1]+Omega_11)
					  -Omega_12*a0[index]
					  -Omega_21[index_eta-1]*(a22[index]+a23[index])
					  +2*exp_eta*AA[_A11_][index])
	+ a11[index];

      a12[index_plus]        = eta_step *(
					  -a12[index]*(2*Omega_22[index_eta-1]+Omega_11)
					  -Omega_21[index_eta-1]*(a23[index]+a21[index])
					  -Omega_12*a0[index]
					  +2*exp_eta*AA[_A12_][index])
	+ a12[index];

      a13[index_plus]        = eta_step *(
					  -a13[index]*(2*Omega_22[index_eta-1]+Omega_11)
					  -Omega_12*a0[index]
					  -Omega_21[index_eta-1]*(a22[index]+a21[index])
					  +2*exp_eta*AA[_A13_][index])
	+ a13[index];

      a21[index_plus]        = eta_step *(
					  -a21[index]*(2*Omega_11+Omega_22[index_eta-1])
					  -Omega_12*(a12[index]+a13[index])
					  -Omega_21[index_eta-1]*a3[index]
					  +2*exp_eta*AA[_A21_][index])
	+ a21[index];

      a22[index_plus]        = eta_step *(
					  -a22[index]*(2*Omega_11+Omega_22[index_eta-1])
					  -Omega_12*(a13[index]+a11[index])
					  -Omega_21[index_eta-1]*a3[index]
					  +2*exp_eta*AA[_A22_][index])
	+ a22[index];

      a23[index_plus]        = eta_step *(
					  -a23[index]*(2*Omega_11+Omega_22[index_eta-1])
					  -Omega_12*(a12[index]+a11[index])
					  -Omega_21[index_eta-1]*a3[index]
					  +2*exp_eta*AA[_A23_][index])
	+ a23[index];

      a3[index_plus]         = eta_step *(
					  -a3[index]*3*Omega_11
					  -Omega_12*(a22[index]+a21[index]+a23[index])
					  +2*exp_eta*AA[_A3_][index])
	+ a3[index];

      b0[index_plus]         = eta_step *(
					  -3*b0[index]*Omega_11
					  -Omega_12*(b11[index]+2*b12[index])
					  +2*exp_eta*AA[_B0_][index])
	+ b0[index];

      b11[index_plus]        = eta_step *(
					  -b11[index]*(2*Omega_11+Omega_22[index_eta-1])
					  -2*Omega_12*b22[index]
					  -Omega_21[index_eta-1]*b0[index]
					  +2*exp_eta*AA[_B11_][index])
	+ b11[index];

      b12[index_plus]        = eta_step *(
					  -b12[index]*(2*Omega_11+Omega_22[index_eta-1])
					  -Omega_12*(b22[index]+b21[index])
					  -Omega_21[index_eta-1]*b0[index]
					  +2*exp_eta*AA[_B12_][index])
	+ b12[index];

      b21[index_plus]        = eta_step *(
					  -b21[index]*(2*Omega_22[index_eta-1]+Omega_11)
					  -2*Omega_21[index_eta-1]*b12[index]
					  -Omega_12*b3[index]
					  +2*exp_eta*AA[_B21_][index])
	+ b21[index];

      b22[index_plus]        = eta_step *(
					  -b22[index]*(2*Omega_22[index_eta-1]+Omega_11)
					  -Omega_21[index_eta-1]*(b12[index]+b11[index])
					  -Omega_12*b3[index]
					  +2*exp_eta*AA[_B22_][index])
	+ b22[index];

      b3[index_plus]         = eta_step *(
					  -3*Omega_22[index_eta-1]*b3[index]
					  -Omega_21[index_eta-1]*(b21[index]+2*b22[index])
					  +2*exp_eta*AA[_B3_][index])
	+ b3[index];
	
    }


    /**********************************************************************/


    /**********
     * Update of second derivatives for interpolation
     **********/
    if(pnl->mode==2){
      class_call(trg_ddp_ab(pnl,pnl->pk_nl,index_eta,pnl->ddpk_nl,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);

      class_call(trg_ddp_ab(pnl,pnl->p_12_nl,index_eta,pnl->ddp_12_nl,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);

      class_call(trg_ddp_ab(pnl,pnl->p_22_nl,index_eta,pnl->ddp_22_nl,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);
    }

    else if(pnl->mode==1){

      for(index_k=0; index_k<pnl->k_size; index_k++){

	class_call(trg_p11_at_k(pba,ppm,psp,pnl,index_eta,0,pnl->k[index_k],&pnl->p_11[index_k+pnl->k_size*index_eta]),
		   pnl->error_message,
		   pnl->error_message);

	class_call(trg_p11_at_k(pba,ppm,psp,pnl,index_eta,0,pnl->k[index_k],&pnl->p_12[index_k+pnl->k_size*index_eta]),
		   pnl->error_message,
		   pnl->error_message);

	class_call(trg_p11_at_k(pba,ppm,psp,pnl,index_eta,0,pnl->k[index_k],&pnl->p_22[index_k+pnl->k_size*index_eta]),
		   pnl->error_message,
		   pnl->error_message);

      }

      class_call(trg_ddp_ab(pnl,pnl->p_11,index_eta,pnl->ddp_11,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);

      class_call(trg_ddp_ab(pnl,pnl->p_12,index_eta,pnl->ddp_12,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);

      class_call(trg_ddp_ab(pnl,pnl->p_22,index_eta,pnl->ddp_22,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);
    }
      


    /**********
     * Update of A's and B's function at the new index_eta which needs
     * the interpolation of P's
     **********/
    
    if(pnl->mode>0){

    /* initialize error management flag */
    abort = _FALSE_;

    /*** beginning of parallel region ***/

#pragma omp parallel							\
  shared(name_size,abort,pba,ppm,psp,pnl,index_eta,AA)			\
  private(tstart,index_name,tstop)
    
    {

#ifdef _OPENMP
      tstart = omp_get_wtime();
#endif

#pragma omp for schedule (dynamic)
      for (index_name=0; index_name<name_size; index_name++) {

#pragma omp flush(abort)

	class_call_parallel(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,index_name,index_eta,AA[index_name],pnl->error_message),
		   pnl->error_message,
		   pnl->error_message);
      }


#ifdef _OPENMP
	tstop = omp_get_wtime();
	if (pnl->spectra_nl_verbose > 2)
	  printf("In %s: time spent in parallel region (loop over names) = %e s for thread %d\n",
		 __func__,tstop-tstart,omp_get_thread_num());
#endif

    } /* end of parallel region */

    if (abort == _TRUE_) return _FAILURE_;

    }

    if(pnl->spectra_nl_verbose > 1)
      printf("    %2.1f%% done\n",100.*index_eta/(pnl->eta_size-1.));
    
    
    time_2=time(NULL);
    if(pnl->spectra_nl_verbose > 0){
      if(index_eta==9){
	printf("  elapsed time after ten loops : %2.f s\n",difftime(time_2, time_1));
	printf("  estimated remaining : %3.1f min\n",difftime(time_2,time_1)*(pnl->eta_size-2)/60/10);
      }
    }
    class_test( isnan(pnl->pk_nl[50+pnl->k_size*index_eta]) !=0 ,
	       pnl->error_message,"Time integrator produces nan !!\n");
	       
	       
  }

  if(pnl->spectra_nl_verbose > 0)
    printf("Done in %2.f min !\n",difftime(time_2,time_1)/60);

  /* Printing results in a file */

  fprintf(nl_spectra,"## k (h/Mpc) P_NL ([Mpc/h]^3) P_L ([Mpc/h]^3) \n");
 
  for(index_eta=4; index_eta<pnl->eta_size; index_eta+=5) {

    fprintf(nl_spectra,"##at z =%e\n",pnl->z[index_eta]);

    for(index_k=0; index_k<pnl->k_size; index_k++){ 

      class_call(
		 spectra_pk_at_k_and_z(pba,ppm,psp,0,index_ic,pnl->k[index_k],pnl->z[index_eta],&temp1),
		 psp->error_message,
		 pnl->error_message);
      
      fprintf(nl_spectra,"%e\t%e\t%e\t",pnl->k[index_k]/pba->h,pow(pba->h,3)*pnl->pk_nl[index_k+(pnl->k_size*(index_eta))]*exp(pnl->eta[index_eta]*2),pow(pba->h,3)*temp1);
      
      for (index_name=0; index_name<name_size; index_name++)
	fprintf(nl_spectra,"%e\t",AA[index_name][index_k+pnl->k_size*index_eta]);
      fprintf(nl_spectra,"\n");
    }
    fprintf(nl_spectra,"\n\n");
  }
  
  fclose(nl_spectra);
  
  


  for (index_name=0; index_name<name_size; index_name++) 
    free(AA[index_name]);
  free(AA);
  
  free(b3);
  free(b22);
  free(b21);
  free(b12);
  free(b11);
  free(b0);
  
  free(a3);
  free(a23);
  free(a22);
  free(a21);
  free(a13);
  free(a12);
  free(a11);
  free(a0);
  
  free(Omega_21);
  free(Omega_22);

  return _SUCCESS_;
   
}

int trg_free(
	     struct spectra_nl * pnl
	     ){
  
  free(pnl->k);
  free(pnl->pk_nl);
  free(pnl->p_12_nl);
  free(pnl->p_22_nl);
  free(pnl->ddpk_nl);
  free(pnl->ddp_12);
  free(pnl->ddp_22);
  free(pnl->eta);

  return _SUCCESS_;  
}
