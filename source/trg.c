/** @file trg.c Document Time Renormalization Group module
 *  Benjamin Audren, 03.05.2010
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
#include <time.h>

ErrorMsg Transmit_Error_Message; /**< contains error message */

double * pvecback_nl;


/* if possible, no other global variable here */

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

int trg_p12_ini(
		struct background * pba, /* all struct are input here */
		struct primordial * ppm,
		struct spectra * psp,
		struct spectra_nl * pnl,
		int index_ic,
		double * result
		){

  double temp1,temp2,temp3;
  int index_k;

  for(index_k=0; index_k<pnl->k_size; index_k++){

    class_call(spectra_pk_at_k_and_z(pba,ppm,psp,0,index_ic,pnl->k[index_k],pnl->z[0],&temp1),
	       psp->error_message,
	       pnl->error_message);
    
    class_call(spectra_pk_at_k_and_z(pba,ppm,psp,0,index_ic,pnl->k[index_k],pnl->z[1],&temp2),
	       psp->error_message,
	       pnl->error_message);

    temp3   = (sqrt(temp2) - sqrt(temp1))/(pnl->z[1]-pnl->z[0]);
    result[index_k] = - (1+pnl->z[0]) * sqrt(temp1) * temp3 ;
  }
  return _SUCCESS_;
}




int trg_p22_ini(
		struct background * pba,
		struct primordial * ppm,
		struct spectra * psp,
		struct spectra_nl * pnl,
		int index_ic,
		double * result
		){
  double temp1,temp2,temp3;
  int index_k;

  for(index_k=0; index_k<pnl->k_size; index_k++){

    class_call(spectra_pk_at_k_and_z(pba,ppm,psp,0,index_ic,pnl->k[index_k],pnl->z[0],&temp1),
	       psp->error_message,
	       pnl->error_message);

    class_call(spectra_pk_at_k_and_z(pba,ppm,psp,0,index_ic,pnl->k[index_k],pnl->z[1],&temp2),
	       psp->error_message,
	       pnl->error_message);
        
    temp3   = (sqrt(temp2) - sqrt(temp1))/ (pnl->z[0] - pnl->z[1]);
    result[index_k] = (1+pnl->z[0])*(1+pnl->z[0])*temp3*temp3;
  }
  return _SUCCESS_;
}


/* As soon as one leaves the initial condition, one needs to
   interpolate pk from the already computed non linear spectrum. This
   is done by the following command */

int trg_pk_nl_ini(
		  struct background * pba,
		  struct primordial * ppm,
		  struct spectra * psp,
		  struct spectra_nl * pnl,
		  int index_ic,
		  double *result
		  ){

  int index_k;

  for(index_k=0; index_k<pnl->k_size; index_k++){
    class_call(spectra_pk_at_k_and_z(pba,ppm,psp,0,index_ic,pnl->k[index_k],pnl->z_ini,&result[index_k]),
	       psp->error_message,
	       pnl->error_message);
  }
  return _SUCCESS_;
}

 /********************
  * Fill the second derivative table of p_ab up to the index_eta
  * subscript
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
   * Argument of the integral of all A's, B's, given at any time eta,
   * usually written as F( k,x+y/sqrt(2),x-y/sqrt(2) ) + F( y <-> -y)
   ********************/

 

  int trg_A_arg(
		struct spectra_nl * pnl,
		enum name_A name, 
		double k, 
		double p, 
		double m, 
		int index_eta,
		double * result, 
		char * errmsg){

    double p_11k,p_11m,p_11p;
    double p_12k,p_12m,p_12p;
    double p_22k,p_22m,p_22p;

    double gamma1_kpm,gamma1_pkm,gamma1_mkp,gamma1_mpk,gamma1_pmk,gamma1_kmp;
    double gamma2_kpm,gamma2_pkm,gamma2_mkp,gamma2_mpk,gamma2_pmk,gamma2_kmp;

/********************
 * Every function needed for the argument of a certain integral is
 * only called inside the switch case instance, in order to lower the
 * execution time - this way the lecture of this function is somehow
 * harder, but the argument is called in the very end of each case
 * bloc. Moreover, the call of functions is in this order : P_ab, the
 * gamma2, then finally gamma1.
 ********************/

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
			     gamma2_pmk*p_22m*p_22k + gamma2_mkp*p_12k*p_22p)
	     + gamma1_kmp *( gamma1_kmp*p_22m*p_12p + gamma1_kpm*p_12m*p_22p + 
			     gamma2_mpk*p_22p*p_22k + gamma2_pkm*p_12k*p_22m);
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

       *result= gamma1_kpm *( gamma2_kpm*p_22p*p_22m + gamma1_pmk*p_22m*p_12k +
			      gamma1_pkm*p_12m*p_22k + gamma2_mkp*p_12k*p_22p)
	      + gamma1_kmp *( gamma2_kmp*p_22m*p_22p + gamma1_mpk*p_22p*p_12k +
			      gamma1_mkp*p_12p*p_22m + gamma2_pkm*p_12k*p_22m);
      return _SUCCESS_;
      break;

/****************************************/

    case _A21_:
      class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,k,&p_11k,errmsg),
		 errmsg,
		 pnl->error_message);
      class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,p,&p_11p,errmsg),
		 errmsg,
		 pnl->error_message);
      class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,m,&p_11m,errmsg),
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

      *result= gamma1_kpm *( gamma2_kpm*p_12p*p_12m + gamma1_pmk*p_12m*p_11k +
			     gamma1_pkm*p_11m*p_12k + gamma1_mkp*p_12k*p_11p +
			     gamma1_mpk*p_11k*p_12p)
	     + gamma1_kmp *( gamma2_kmp*p_12m*p_12p + gamma1_mpk*p_12p*p_11k +
			     gamma1_mkp*p_11p*p_12k + gamma1_pkm*p_12k*p_11m +
			     gamma1_pmk*p_11k*p_12m);
      return _SUCCESS_;
      break;

/****************************************/

    case _A22_:
      class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,p,&p_11p,errmsg),
		 errmsg,
		 pnl->error_message);
      class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,m,&p_11m,errmsg),
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
			     gamma2_pmk*p_22m*p_12k + gamma1_mkp*p_22k*p_11p +
			     gamma1_mpk*p_12k*p_12p)
	     + gamma1_kmp *( gamma1_kmp*p_22m*p_11p + gamma1_kpm*p_12m*p_12p +
			     gamma2_mpk*p_22p*p_12k + gamma1_pkm*p_22k*p_11m +
			     gamma1_pmk*p_12k*p_12m);
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
    

      trg_gamma_121(p,k,m,&gamma1_pkm);
      trg_gamma_121(p,m,k,&gamma1_pmk);
      trg_gamma_121(m,p,k,&gamma1_mpk);
      trg_gamma_121(m,k,p,&gamma1_mkp);


      *result= gamma2_kpm *( gamma2_kpm*p_12p*p_12m + gamma1_pmk*p_12m*p_11k +
			     gamma1_pkm*p_11m*p_12k + gamma1_mkp*p_12k*p_11p +
			     gamma1_mpk*p_11k*p_12p)
	     + gamma2_kmp *( gamma2_kmp*p_12m*p_12p + gamma1_mpk*p_12p*p_11k +
			     gamma1_mkp*p_11p*p_12k + gamma1_pkm*p_12k*p_11m +
			     gamma1_pmk*p_11k*p_12m);
      return _SUCCESS_;
      break;

/****************************************/

    case _B12_:
      class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,p,&p_11p,errmsg),
		 errmsg,
		 pnl->error_message);
      class_call(trg_p_ab_at_any_k(pnl,pnl->pk_nl,pnl->ddpk_nl,index_eta,m,&p_11m,errmsg),
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
			      gamma2_pmk*p_22m*p_12k + gamma1_mkp*p_22k*p_11p +
			      gamma1_mpk*p_12k*p_12p)
	      + gamma2_kmp *( gamma1_kmp*p_22m*p_11p + gamma1_kpm*p_12m*p_12p +
			      gamma2_mpk*p_22p*p_12k + gamma1_pkm*p_22k*p_11m +
			      gamma1_pmk*p_12k*p_12m);
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
			     gamma2_pmk*p_22m*p_22k + gamma2_mkp*p_12k*p_22p)
	     + gamma2_kmp *( gamma1_kmp*p_22m*p_12p + gamma1_kpm*p_12m*p_22p + 
			     gamma2_mpk*p_22p*p_22k + gamma2_pkm*p_12k*p_22m);
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
      class_call(trg_gamma_222(m,k,p,&gamma2_mkp,errmsg),
		 errmsg,
		 pnl->error_message);
  
    
      trg_gamma_121(p,k,m,&gamma1_pkm);
      trg_gamma_121(p,m,k,&gamma1_pmk);
      trg_gamma_121(m,p,k,&gamma1_mpk);
      trg_gamma_121(m,k,p,&gamma1_mkp);

      *result= gamma2_kpm *( gamma2_kpm*p_22p*p_22m + gamma1_pmk*p_22m*p_12k +
			     gamma1_pkm*p_12m*p_22k + gamma2_mkp*p_12k*p_22p)
	     + gamma2_kmp *( gamma2_kmp*p_22m*p_22p + gamma1_mpk*p_22p*p_12k +
			     gamma1_mkp*p_12p*p_22m + gamma2_pkm*p_12k*p_22m);
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
      return _SUCCESS_;
      break;

/****************************************/

      default:
	sprintf(pnl->error_message,"%s(L:%d): non valid argument in integrals A of I, %d is not defined\n",__func__,__LINE__,name);
	return _FAILURE_;
	break;
    }

  }

int trg_integrate_xy_at_eta(
			    struct spectra_nl * pnl,
			    enum name_A name,
			    int index_eta,
			    int k_size, /* functions of k for
					   integration are defined
					   from [0+k_size*index_eta]
					   to
					   [k_size-1+k_size*index_eta]*/
			    double * result,
			    char *errmsg
			    ){
  
  int index_k;
  double k_max=pnl->k[pnl->k_size-1]; /* choice of maximum k when integrating over x*/
  double *x,*y;
  double x_calc,y_calc;
  double logstep;
  double temp;
  int index_y,index_x;
  int x_size,y_size;
  
  double arg_plus,arg_minus;
  double x_min,x_max;
  double y_min,y_max;
  
  double *sum_y;
  double diff_sum_y;

  double epsilon; /* cut-off scale to avoid numerical divergence*/
  double epsilon_prime; /* smallest change in integrands that we want
			   to consider, one should test its influence
			   on the final result, if any*/

  /* The following approximation consists in taking into account the
     fact that for small k's, the influence of non linearity is
     small. Hence, one could take the linear evolution for P's for
     sufficiently small k's. Here, arbitrarily, one takes values of
     index_k before 10 to follow linear regime, i.e. putting to zero
     their A's. 

     This approximation is here to calculate more precisely the first
     k points on x-y domains. The other solution (cleaner, but more
     time consuming) would consist in modifying the lower bound in the
     table_array routines by chosing an arbitrarily small cut-off for
     P's (well under the k_min).

  */

  epsilon_prime=1e-6;

  epsilon=pnl->k[0];

  diff_sum_y=1;

  for(index_k=0; index_k<pnl->k_size; index_k++){
    
    if(index_k<pnl->index_k_1){
      result[index_k+k_size*index_eta]=0;
    }
    
    else{
      x_min=pnl->k[index_k]/sqrt(2);
      
      logstep=1+2*epsilon/pnl->k[index_k];
      
      /*checking number of integration steps for x */

      x_calc=x_min;
      
      index_x=0;
      while(x_calc<k_max*sqrt(2)){
	index_x++;
	x_calc=x_min*pow(logstep,index_x);
      }
      
      x_size=index_x+1;
      
      /*recording of x values */
      
      x=calloc(x_size,sizeof(double));
      x[0]=x_min;
      for(index_x=1; index_x<x_size-1; index_x++){
	x[index_x]=x[index_x-1]*logstep;
      }
      x[x_size-1]=sqrt(2)*k_max;
      
      /*checking number of integration steps for y*/
      
      index_y=0;
      y_calc=x_min;
      while(y_calc<(2*x_min)){
	index_y++;
	y_calc=x_min*pow(logstep,index_y);
      }
      
      y_size=index_y+1;
      
      /*recording of y values*/
      
      y=calloc(y_size,sizeof(double));
      y[0]=pnl->k[index_k]/sqrt(2);
      y_calc=y[0];
      for(index_y=1; index_y<y_size-1; index_y++){
	y_calc=y[0]*pow(logstep,index_y);
	y[index_y]=2*y[0]-y_calc;
      }
      y[y_size-1]=0;
      
      /* Loops */
      
      sum_y=calloc(x_size,sizeof(double));
          
	/* beware the sign change to compensate integration
	   boundaries */
	
	/* For each x, perform integration on y */
      
	/* First, the case index_x==0 , the integration starts at index_y=1*/
      
      for(index_x=0; index_x<x_size; index_x++){
	if(index_x==0){
	  
	  if((x[0]-y[1])<epsilon*sqrt(2)){
	    class_call(trg_A_arg(pnl,name,pnl->k[index_k],(x[0]+y[1])/sqrt(2),epsilon,index_eta,&temp,errmsg),
		       errmsg,
		       pnl->error_message);
	  }
	  else {
	    class_call(trg_A_arg(pnl,name,pnl->k[index_k],(x[0]+y[1])/sqrt(2),(x[0]-y[1])/sqrt(2),index_eta,&temp,errmsg),
		       errmsg,
		       pnl->error_message);
	  }
	  arg_plus=-(x[0]*x[0]-y[1]*y[1])/4 * temp;
	  
	  for(index_y=2; index_y<y_size; index_y++){
	    arg_minus=arg_plus;
	    if((x[0]+y[index_y])>k_max*sqrt(2)){
	      temp=0;
	    }
	    else {
	      class_call(trg_A_arg(pnl,name,pnl->k[index_k],(x[0]+y[index_y])/sqrt(2),(x[0]-y[index_y])/sqrt(2),index_eta,&temp,errmsg),
			 errmsg,
			 pnl->error_message);
	      
	      arg_plus=-(x[0]*x[0]-y[index_y]*y[index_y])/4 * temp;
	      sum_y[0]+=(y[index_y]-y[index_y-1])*0.5*(arg_plus+arg_minus);
	    }
	  }
	}
	
	/* Now the rest of the cases */
	
	else if(fabs(diff_sum_y)<epsilon_prime){
	  sum_y[index_x]=0;
	  diff_sum_y=0;
	  
	}
	else{
	  if((x[index_x]+y[0])>k_max*sqrt(2)){
	    temp=0;
	  }
	  else if((x[index_x]-y[0])<epsilon*sqrt(2)){
	    class_call(trg_A_arg(pnl,name,pnl->k[index_k],(x[index_x]+y[0])/sqrt(2),epsilon,index_eta,&temp,errmsg),
		       errmsg,
		       pnl->error_message);
	    
	  }
	  else{
	    class_call(trg_A_arg(pnl,name,pnl->k[index_k],(x[index_x]+y[0])/sqrt(2),(x[index_x]-y[0])/sqrt(2),index_eta,&temp,errmsg),
		       errmsg,
		       pnl->error_message);
	  }
	  arg_plus=-(x[index_x]*x[index_x]-y[0]*y[0])/4*temp;
	  
	  for(index_y=1; index_y<y_size; index_y++){
	    arg_minus=arg_plus;
	    if((x[index_x]+y[index_y])>k_max*sqrt(2)){
	      temp=0;
	    }
	    else if((x[index_x]-y[index_y])<epsilon*sqrt(2)){
	      class_call(trg_A_arg(pnl,name,pnl->k[index_k],(x[index_x]+y[index_y])/sqrt(2),epsilon,index_eta,&temp,errmsg),
			 errmsg,
			 pnl->error_message);
	    }
	    else{
	      class_call(trg_A_arg(pnl,name,pnl->k[index_k],(x[index_x]+y[index_y])/sqrt(2),(x[index_x]-y[index_y])/sqrt(2),index_eta,&temp,errmsg),
			 errmsg,
			 pnl->error_message);
	    }
	    arg_plus=-(x[index_x]*x[index_x]-y[index_y]*y[index_y])/4*temp;
	    sum_y[index_x]+=(y[index_y]-y[index_y-1])*0.5*(arg_plus+arg_minus);
	  }
	  diff_sum_y=sum_y[index_x]-sum_y[index_x-1];
	}
      }
      for(index_x=1; index_x<x_size; index_x++){
	  result[index_k+k_size*index_eta]+=(x[index_x]-x[index_x-1])*0.5*(sum_y[index_x]+sum_y[index_x-1]); /* integration over x */
	}
	
	free(sum_y);
	free(x);
	free(y);
    }
  }
  return _SUCCESS_;
}





   /*******************************
   *
   * Summary : taking the power spectrum of matter, matter-velocity
   *           and velocity, at a chosen time a_ini when matter dominated
   *           already but before all modes become non linear, it computes the
   *           fully non-linear evolution of these 3 quantities for a range of
   *           eta=log(a/a_ini) up to today.
   *
   *	       It requires the definition of many quantities (12+12) built
   *	       on the Bispectra of matter and velocity and on the spectra
   *           as well, in order to put the equations in a more
   *           numerically friendly shape, as explained in annex B of the
   *	       Pietroni article (astro-ph 0806.0971v3).
   *
   *********************************/

int trg_init (
	      struct precision * ppr, /* input */
	      struct background * pba, /**< structure containing all background evolution */
	      struct primordial * ppm,
	      struct spectra * psp, /**< structure containing list of k, z and P(k,z) values) */
	      struct spectra_nl * pnl /* output */ 
	      ) {

  FILE *nl_spectra;

  time_t time_1,time_2;

  int index_k;
  int index_eta;

  int index_ic=0; /*or ppt->index_ic; adiabatic=0 */

  double temp_k;
  

  /*
    Variables of calculation, k and eta, and initialization
  */

  double a_ini;
  double z_ini;
  double eta_max;
  double * Omega_m, * H, *H_prime;

  double eta_step; /* step for integrations */

  double exp_eta;
  double pi_over_k;
  int index;  
  int index_plus; /* intermediate variables for time integration */
  
  double Omega_11,Omega_12; /* Definition of the matrix Omega that mixes terms together*/
  double *Omega_21,*Omega_22;

  double *a0,*a11,*a12,*a21,*a22,*a3;
  double *b0,*b11,*b12,*b21,*b22,*b3; /* Definition of the variables
					 1ijk,lmn and 2 ijk,lmn that
					 replace the Bispectra
					 variables*/

  double *A0,*A11,*A12,*A21,*A22,*A3;
  double *B0,*B11,*B12,*B21,*B22,*B3; /* Definition of the variables
					 Aijk,lmn */

  pnl->spectra_nl_verbose=1;

  pvecback_nl = calloc(pba->bg_size,sizeof(double));

  if (pnl->spectra_nl_verbose > 0)
    printf("Computing non-linear spectra with TRG method\n");

  /* define initial eta and redshift */
  a_ini = 2e-2;
  pnl->z_ini = pba->a_today/a_ini - 1.;

  /* define eta_max, where eta=log(a/a_ini) */
  eta_max = log(pba->a_today/a_ini);

  /* define size and step for integration in eta */
  pnl->eta_size = 200;
  pnl->eta_step = (eta_max)/(pnl->eta_size-1);
  eta_step = pnl->eta_step;

  /* at any time, a = a_ini * exp(eta) and z=exp(eta)-1*/
  
  /* fill array of k values  */

  /* first define the total length, to reach the k_max bound
     for the integrator */

  double k_max;
  double k_min;
  double k_1;
  double k_2;

  k_1=1e-3;

  k_max=70;
  k_min=7e-5;

  index_k=0;
  temp_k=k_min;

  while(temp_k<k_max) {
    index_k++;
    temp_k=k_min*pow(1.2,index_k); /* in order to have roughly 10 points
				   per decade */
  }

  pnl->k_size=index_k+1;

  pnl->k = calloc(pnl->k_size,sizeof(double));

  /* then fill in the values of k, while determining the actual upper
     bound k_2 for calculation so as to avoid integration
     troubles. This k_2 is set arbitrary to 1.5 h/Mpc to reproduce
     existing data. */

  temp_k=0;

  for(index_k=0; index_k<pnl->k_size-1; index_k++){
    pnl->k[index_k]=k_min*pow(1.2,index_k);
    if( (pnl->k[index_k] > k_1*pba->h) && temp_k==0){
      pnl->index_k_1=index_k;
      temp_k++;
    }
    if( (pnl->k[index_k] > 1.5*pba->h) && temp_k==1 ){
      pnl->k_size_eff=index_k+1;
      k_2=pnl->k[index_k];
      temp_k++;
    }
    printf("%d %e\n",index_k,pnl->k[index_k]);
  }
  pnl->k[pnl->k_size-1]=k_max;



  /* fill array of eta values, and pick two values z_1 and z_2
     resp. near 2.33 and 1.00 */

  double z_1,z_2;

  z_1=2.33;
  z_2=1.00;

  int index_eta_1,index_eta_2;

  pnl->eta = calloc(pnl->eta_size,sizeof(double));
  pnl->z   = calloc(pnl->eta_size,sizeof(double));

  int temp_z;

  temp_z=0;

  printf("\n\n");

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
    printf("%d %e %e\n",index_eta,pnl->eta[index_eta],pnl->z[index_eta]);
  }

  if (pnl->spectra_nl_verbose > 0)
    printf(" -> starting calculation at redshift z = %2.2f\n",pnl->z[0]);

  /* definition of background values for each eta */

  Omega_m = calloc(pnl->eta_size,sizeof(double));
  H = calloc(pnl->eta_size,sizeof(double));
  H_prime = calloc(pnl->eta_size,sizeof(double));

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
  
  Omega_21 = calloc(pnl->eta_size,sizeof(double));
  Omega_22 = calloc(pnl->eta_size,sizeof(double));

  for(index_eta=0; index_eta<pnl->eta_size; index_eta++) {
    Omega_21[index_eta] = -3./2 * Omega_m[index_eta];
    Omega_22[index_eta] = 2 + H_prime[index_eta]/H[index_eta];
  } 

  /* Definition of P_11=pk_nl, P_12=<delta theta> and P_22=<theta
     theta>, and initialization at eta[0]=0, k[0] 

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

  pnl->pk_nl =  calloc(pnl->k_size*pnl->eta_size,sizeof(double));
  pnl->p_12  =  calloc(pnl->k_size*pnl->eta_size,sizeof(double));
  pnl->p_22  =  calloc(pnl->k_size*pnl->eta_size,sizeof(double));


  if(pnl->pk_nl != NULL && pnl->p_12 != NULL && pnl->p_22 != NULL){
    class_call(trg_pk_nl_ini(pba,ppm,psp,pnl,index_ic,pnl->pk_nl),
	       pnl->error_message,
	       pnl->error_message);
    
    class_call(trg_p12_ini(pba,ppm,psp,pnl,index_ic,pnl->p_12),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_p22_ini(pba,ppm,psp,pnl,index_ic,pnl->p_22),
	       pnl->error_message,
	       pnl->error_message);

  }
  else {
    sprintf(Transmit_Error_Message,"%s(L:%d): could not allocate memory !",__func__,__LINE__);
    sprintf(pnl->error_message,"%s",Transmit_Error_Message);
    return _FAILURE_;
  }

  /* Initialisation and definition of second derivatives */

  pnl->ddpk_nl=calloc(pnl->k_size*pnl->eta_size,sizeof(double));
  pnl->ddp_12=calloc(pnl->k_size*pnl->eta_size,sizeof(double));
  pnl->ddp_22=calloc(pnl->k_size*pnl->eta_size,sizeof(double));
  
  if(pnl->pk_nl != NULL && pnl->p_12 != NULL && pnl->p_22 != NULL){
    class_call(trg_ddp_ab(pnl,pnl->pk_nl,0,pnl->ddpk_nl,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_ddp_ab(pnl,pnl->p_12,0,pnl->ddp_12,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_ddp_ab(pnl,pnl->p_22,0,pnl->ddp_22,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
  }
  else {
    sprintf(Transmit_Error_Message,"%s(L:%d): could not allocate memory !",__func__,__LINE__);
    sprintf(pnl->error_message,"%s",Transmit_Error_Message);
    return _FAILURE_;
  }

  /* Definition of 1_0, 1_11,(here a0, a11,...) etc, and 2_0, 2_11,
     (here b0,b11,...) etc.. and initialization (directly with calloc
     for assuming no initial non gaussianity in the spectrum) */

  a0  = calloc(pnl->k_size*pnl->eta_size,sizeof(double));
  a11 = calloc(pnl->k_size*pnl->eta_size,sizeof(double));
  a12 = calloc(pnl->k_size*pnl->eta_size,sizeof(double));
  a21 = calloc(pnl->k_size*pnl->eta_size,sizeof(double));
  a22 = calloc(pnl->k_size*pnl->eta_size,sizeof(double));
  a3  = calloc(pnl->k_size*pnl->eta_size,sizeof(double));

  b0  = calloc(pnl->k_size*pnl->eta_size,sizeof(double));
  b11 = calloc(pnl->k_size*pnl->eta_size,sizeof(double));
  b12 = calloc(pnl->k_size*pnl->eta_size,sizeof(double));
  b21 = calloc(pnl->k_size*pnl->eta_size,sizeof(double));
  b22 = calloc(pnl->k_size*pnl->eta_size,sizeof(double));
  b3  = calloc(pnl->k_size*pnl->eta_size,sizeof(double)); 

  /* Definition of the A_121,122's that appear in time-evolution of
     the 1's and 2's, and initialization with pk_nl, p_12,
     p_22. Convention is taken for A_121,any = A any and A_222,any = B any */

  A0  = calloc(pnl->k_size*pnl->eta_size,sizeof(double));
  A11 = calloc(pnl->k_size*pnl->eta_size,sizeof(double));
  A12 = calloc(pnl->k_size*pnl->eta_size,sizeof(double));
  A21 = calloc(pnl->k_size*pnl->eta_size,sizeof(double));
  A22 = calloc(pnl->k_size*pnl->eta_size,sizeof(double));
  A3  = calloc(pnl->k_size*pnl->eta_size,sizeof(double));

  B0  = calloc(pnl->k_size*pnl->eta_size,sizeof(double));
  B11 = calloc(pnl->k_size*pnl->eta_size,sizeof(double));
  B12 = calloc(pnl->k_size*pnl->eta_size,sizeof(double));
  B21 = calloc(pnl->k_size*pnl->eta_size,sizeof(double));
  B22 = calloc(pnl->k_size*pnl->eta_size,sizeof(double));
  B3  = calloc(pnl->k_size*pnl->eta_size,sizeof(double));

  /********************
   *
   * Integration at eta=0 with n_xy number of steps for integration
   * over x and n_xy over y.
   *
   ********************/
  double temp1,temp2,temp3;

  if (pnl->spectra_nl_verbose > 0)
    printf(" -> initialisation\n");
  
  class_call(trg_integrate_xy_at_eta(pnl,_A0_,0,pnl->k_size,A0,pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);
				     
  class_call(trg_integrate_xy_at_eta(pnl,_A11_,0,pnl->k_size,A11,pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);
  
  class_call(trg_integrate_xy_at_eta(pnl,_A12_,0,pnl->k_size,A12,pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);

  class_call(trg_integrate_xy_at_eta(pnl,_A21_,0,pnl->k_size,A21,pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);

  class_call(trg_integrate_xy_at_eta(pnl,_A22_,0,pnl->k_size,A22,pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);

  class_call(trg_integrate_xy_at_eta(pnl,_A3_,0,pnl->k_size,A3,pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);

  class_call(trg_integrate_xy_at_eta(pnl,_B0_,0,pnl->k_size,B0,pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);
  
  class_call(trg_integrate_xy_at_eta(pnl,_B11_,0,pnl->k_size,B11,pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);

  class_call(trg_integrate_xy_at_eta(pnl,_B12_,0,pnl->k_size,B12,pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);

  class_call(trg_integrate_xy_at_eta(pnl,_B21_,0,pnl->k_size,B21,pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);

  class_call(trg_integrate_xy_at_eta(pnl,_B22_,0,pnl->k_size,B22,pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);

  class_call(trg_integrate_xy_at_eta(pnl,_B3_,0,pnl->k_size,B3,pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);



  /********************
   * Now we calculate the time evolution with a very simple integrator
   ********************/

  if((nl_spectra = fopen("nl_spectra.dat","wr")) == NULL)
    printf("could not create nl_spectra.dat !!\n");

  if (pnl->spectra_nl_verbose > 0){
    printf(" -> progression:\n");}

  time_1=time(NULL);
 
  for (index_eta=1; index_eta<pnl->eta_size; index_eta++){
    exp_eta=exp(pnl->eta[index_eta-1]);

    for (index_k=0; index_k<pnl->k_size; index_k++){

      /**********
       * Computation for every k of the function at the next time step
       **********/

      pi_over_k=4*_PI_/(pnl->k[index_k]);
      index = index_k+pnl->k_size*(index_eta-1);
      index_plus = index_k+pnl->k_size*index_eta;

     

      pnl->pk_nl[index_plus]= eta_step *(
					 -2*Omega_11*pnl->pk_nl[index]
					 -2*Omega_12*pnl->p_12[index]
					 +exp_eta*4*pi_over_k*a22[index] )
	                    + pnl->pk_nl[index];

      pnl->p_22[index_plus] = eta_step *(
					 -2*Omega_22[index_eta-1]*pnl->p_22[index]
					 -2*Omega_21[index_eta-1]*pnl->p_12[index]
					 +exp_eta*2*pi_over_k*b3[index] )
	                    + pnl->p_22[index];

      pnl->p_12[index_plus] = eta_step *(
					 -pnl->p_12[index]*(Omega_11+Omega_22[index_eta-1])
					 -Omega_12*pnl->p_22[index]
					 -Omega_21[index_eta-1]*pnl->pk_nl[index]
					 +exp_eta*pi_over_k*(2*a22[index]+b21[index]))
	                     + pnl->p_12[index];


      a0[index_plus]         = eta_step *(
					  -Omega_21[index_eta-1]*(a11[index]+2*a12[index])
					  -3*Omega_22[index_eta-1]*a0[index]
					  +2*exp_eta*A0[index])
	                     + a0[index];

      a11[index_plus]        = eta_step *(
					  -a11[index]*(2*Omega_22[index_eta-1]+Omega_11)
					  -Omega_12*a0[index]
					  -2*Omega_21[index_eta-1]*a22[index]
					  +2*exp_eta*A11[index])
	                     + a11[index];

      a12[index_plus]        = eta_step *(
					  -a12[index]*(2*Omega_22[index_eta-1]+Omega_11)
					  -Omega_21[index_eta-1]*(a22[index]+a21[index])
					  -Omega_12*a0[index]
					  +2*exp_eta*A12[index])
	                     + a12[index];

      a21[index_plus]        = eta_step *(
					  -a21[index]*(2*Omega_11+Omega_22[index_eta-1])
					  -2*Omega_12*a12[index]
					  -Omega_21[index_eta-1]*a3[index]
					  +2*exp_eta*A21[index])
	                     + a21[index];

      a22[index_plus]        = eta_step *(
					  -a22[index]*(2*Omega_11+Omega_22[index_eta-1])
					  -Omega_12*(a12[index]+a11[index])
					  -Omega_21[index_eta-1]*a3[index]
					  +2*exp_eta*A22[index])
	                     + a22[index];

      a3[index_plus]         = eta_step *(
					  -a3[index]*3*Omega_11
					  -Omega_12*(2*a22[index]+a21[index])
					  +2*exp_eta*A3[index])
	                     + a3[index];

      b0[index_plus]         = eta_step *(
					  -3*b0[index]*Omega_11
					  -Omega_12*(b11[index]+2*b12[index])
					  +2*exp_eta*B0[index])
 	                     + b0[index];

      b11[index_plus]        = eta_step *(
					  -b11[index]*(2*Omega_11+Omega_22[index_eta-1])
					  -2*Omega_12*b22[index]
					  -Omega_21[index_eta-1]*b0[index]
					  +2*exp_eta*B11[index])
	                     + b11[index];

      b12[index_plus]        = eta_step *(
					  -b12[index]*(2*Omega_11+Omega_22[index_eta-1])
					  -Omega_12*(b22[index]+b21[index])
					  -Omega_21[index_eta-1]*b0[index]
					  +2*exp_eta*B12[index])
	                     + b12[index];

      b21[index_plus]        = eta_step *(
					  -b21[index]*(2*Omega_22[index_eta-1]+Omega_11)
					  -2*Omega_21[index_eta-1]*b12[index]
					  -Omega_12*b3[index]
					  +2*exp_eta*B21[index])
	                     + b21[index];

      b22[index_plus]        = eta_step *(
					  -b22[index]*(2*Omega_22[index_eta-1]+Omega_11)
					  -Omega_21[index_eta-1]*(b12[index]+b11[index])
					  -Omega_12*b3[index]
					  +2*exp_eta*B22[index])
	                     + b22[index];

      b3[index_plus]         = eta_step *(
					  -3*Omega_22[index_eta-1]*b3[index]
					  -Omega_21[index_eta-1]*(b21[index]+2*b22[index])
					  +2*exp_eta*B3[index])
	                     + b3[index];
	
    }

    /**********
     * Update of second derivatives for interpolation
     **********/

    class_call(trg_ddp_ab(pnl,pnl->pk_nl,index_eta,pnl->ddpk_nl,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_ddp_ab(pnl,pnl->p_12,index_eta,pnl->ddp_12,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_ddp_ab(pnl,pnl->p_22,index_eta,pnl->ddp_22,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);


    /**********
     * Update of A's and B's function at the new index_eta
     **********/
    
    class_call(trg_integrate_xy_at_eta(pnl,_A0_,index_eta,pnl->k_size,A0,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_integrate_xy_at_eta(pnl,_A11_,index_eta,pnl->k_size,A11,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_integrate_xy_at_eta(pnl,_A12_,index_eta,pnl->k_size,A12,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_integrate_xy_at_eta(pnl,_A21_,index_eta,pnl->k_size,A21,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_integrate_xy_at_eta(pnl,_A22_,index_eta,pnl->k_size,A22,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_integrate_xy_at_eta(pnl,_A3_,index_eta,pnl->k_size,A3,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_integrate_xy_at_eta(pnl,_B0_,index_eta,pnl->k_size,B0,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_integrate_xy_at_eta(pnl,_B11_,index_eta,pnl->k_size,B11,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_integrate_xy_at_eta(pnl,_B12_,index_eta,pnl->k_size,B12,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_integrate_xy_at_eta(pnl,_B21_,index_eta,pnl->k_size,B21,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_integrate_xy_at_eta(pnl,_B22_,index_eta,pnl->k_size,B22,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_integrate_xy_at_eta(pnl,_B3_,index_eta,pnl->k_size,B3,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    
    printf("    %2.1f%% done\n",100.*index_eta/(pnl->eta_size-1.));

    time_2=time(NULL);
    /*    if(index_eta==9){
      printf("elapsed time after ten loops : %2.f s\n",difftime(time_2, time_1));
      printf("estimated remaining : %3.1f min\n",difftime(time_2,time_1)*(pnl->eta_size-2)/60/10);
    }
    */
    if(isnan(pnl->pk_nl[50+pnl->k_size*index_eta])!=0){
      printf("ca marche pas !!!nan!!!\n");
    }
  }

  printf("Done !\n");

  /* printing in a file */

  fprintf(nl_spectra,"##at z_1=%e\n ",z_1);

  for(index_k=0; index_k<pnl->k_size; index_k++){
    class_call(
	       spectra_pk_at_k_and_z(pba,ppm,psp,0,index_ic,pnl->k[index_k],pnl->z[index_eta_1],&temp1),
	       psp->error_message,
	       pnl->error_message);

    fprintf(nl_spectra,"%e\t%e\t%e\n",pnl->k[index_k]/pba->h,pow(pba->h,3)*pnl->pk_nl[index_k+(pnl->k_size*(index_eta_1))]*exp(pnl->eta[index_eta_1]*2),pow(pba->h,3)*temp1);
  }

  fprintf(nl_spectra,"\n\n##at z_2=%e\n ",z_2);

  for(index_k=0; index_k<pnl->k_size; index_k++){
    class_call(
	       spectra_pk_at_k_and_z(pba,ppm,psp,0,index_ic,pnl->k[index_k],pnl->z[index_eta_2],&temp1),
	       psp->error_message,
	       pnl->error_message);

    fprintf(nl_spectra,"%e\t%e\t%e\n",pnl->k[index_k]/pba->h,pow(pba->h,3)*pnl->pk_nl[index_k+(pnl->k_size*(index_eta_2))]*exp(pnl->eta[index_eta_2]*2),pow(pba->h,3)*temp1);
  }
  /*
  fprintf(nl_spectra,"\n\n##for %d values of k\n##z is equal to %f\n## k\tpk_nl\tpk_lin at last eta(today)\n",pnl->k_size_eff,pnl->z[pnl->eta_size-1]);
  
  for(index_k=0; index_k<pnl->k_size_eff; index_k++){
    class_call(
	       spectra_pk_at_k_and_z(pba,ppm,psp,0,index_ic,pnl->k[index_k],pnl->z[pnl->eta_size-1],&temp1),
	       psp->error_message,
	       pnl->error_message);

    fprintf(nl_spectra,"%e\t%e\t%e\n",pnl->k[index_k]/pba->h,pow(pba->h,3)*pnl->pk_nl[index_k+(pnl->k_size*(pnl->eta_size-1))]*exp(pnl->eta[pnl->eta_size-1]*2),pow(pba->h,3)*temp1);
  }
  */

  fclose(nl_spectra);
  
  free(B3);
  free(B22);
  free(B21);
  free(B12);
  free(B11);
  free(B0);

  free(A3);
  free(A22);
  free(A21);
  free(A12);
  free(A11);
  free(A0);

  free(b3);
  free(b22);
  free(b21);
  free(b12);
  free(b11);
  free(b0);

  free(a3);
  free(a22);
  free(a21);
  free(a12);
  free(a11);
  free(a0);

  return _SUCCESS_;
   
}

int trg_free(
	     struct spectra_nl * pnl
	     ){
  
  free(pnl->k);
  free(pnl->pk_nl);
  free(pnl->p_12);
  free(pnl->p_22);
  free(pnl->ddpk_nl);
  free(pnl->ddp_12);
  free(pnl->ddp_22);
  free(pnl->eta);

  return _SUCCESS_;  
}
