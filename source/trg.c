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
#include "time.h"

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

/* calculate spectrum p_12 at index_eta and k */

int trg_p12_at_k(
		 struct background * pba, /* all struct are input here */
		 struct primordial * ppm,
		 struct spectra * psp,
		 struct spectra_nl * pnl,
		 int index_eta,
		 int index_ic,
		 double k,
		 double * result
		 ){

  double temp1,temp2,temp3;
  
  class_call(spectra_pk_at_k_and_z(pba,ppm,psp,0,index_ic,k,pnl->z[index_eta],&temp1),
	     psp->error_message,
	     pnl->error_message);
    
  class_call(spectra_pk_at_k_and_z(pba,ppm,psp,0,index_ic,k,pnl->z[index_eta+1],&temp2),
	     psp->error_message,
	     pnl->error_message);

  temp3   = (sqrt(temp2) - sqrt(temp1))/(pnl->z[index_eta+1]-pnl->z[index_eta]);
  *result = - exp(-2*pnl->eta[index_eta])*(1+pnl->z[index_eta]) * sqrt(temp1) * temp3 ;
  
  return _SUCCESS_;
}

/* calculate spectrum p_22 at index_eta and k */


int trg_p22_at_k(
		 struct background * pba,
		 struct primordial * ppm,
		 struct spectra * psp,
		 struct spectra_nl * pnl,
		 int index_eta,
		 int index_ic,
		 double k,
		 double * result
		 ){
  double temp1,temp2,temp3;
  
  class_call(spectra_pk_at_k_and_z(pba,ppm,psp,0,index_ic,k,pnl->z[index_eta],&temp1),
	     psp->error_message,
	     pnl->error_message);
  
  class_call(spectra_pk_at_k_and_z(pba,ppm,psp,0,index_ic,k,pnl->z[index_eta+1],&temp2),
	     psp->error_message,
	     pnl->error_message);
        
  temp3   = (sqrt(temp2) - sqrt(temp1))/ (pnl->z[index_eta+1] - pnl->z[index_eta]);
  *result = exp(-2*pnl->eta[index_eta])*(1+pnl->z[index_eta])*(1+pnl->z[index_eta])*temp3*temp3;
  
  return _SUCCESS_;
}


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

/* As soon as one leaves the initial condition, one needs to
   interpolate pk from the already computed non linear spectrum. This
   is done by the following two commands */

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

/**** interpolate any structure *****/

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
 * p stands for x+y/sqrt(2) (plus) and m for x-y/sqrt(2) (minus)
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

  /* our scheme is such that the first value of (x-y)/sqrt(2)=m is k_min. However
     there might be a small rounding error. After checking that this error is at most of 0.1per cent,
     impose manually that the minimum (x-y)/sqrt(2) is exactly epsilon */

  class_test(m < 0.999*pnl->k[0],
	     pnl->error_message,
	     "error in cut-off implementation (k=%e, (x-y)/sqrt(2)=%e x=%e y=%e k_min=%e)",
	     k,m,(p+m)/sqrt(2),(p-m)/sqrt(2),pnl->k[0]);
  if (m<pnl->k[0]) m=pnl->k[0];

  /* set argument to zero when p is not between the two cut-off values */

  if ((p<pnl->k[0]) || (p>pnl->k[pnl->k_size-1])) {
    *result=0.;
    return _SUCCESS_;
  }

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

    /* if (index_k==103) */
    if (index_k==30 || index_k==pnl->k_size-1)
      printf("%e %e %e %e %e %e\n",k,p,m,(p+m)/sqrt(2.),(p-m)/sqrt(2.),p*m*(*result)); 
    

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

    if (index_k==30 || index_k==pnl->k_size-1)
      printf("%e %e %e %e %e %e\n",k,p,m,(p+m)/sqrt(2.),(p-m)/sqrt(2.),p*m*(*result)); 

    return _SUCCESS_;
    break;

    /****************************************/

  default:
    sprintf(pnl->error_message,"%s(L:%d): non valid argument in integrals A of I, %d is not defined\n",__func__,__LINE__,name);
    return _FAILURE_;
    break;
  }

}


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

  /* our scheme is such that the first value of (x-y)/sqrt(2)=m is k_min. However
     there might be a small rounding error. After checking that this error is at most of 0.1per cent,
     impose manually that the minimum (x-y)/sqrt(2) is exactly epsilon */

  /*   class_test(m < 0.999*pnl->k[0], */
  /* 	     pnl->error_message, */
  /* 	     "error in cut-off implementation (k=%e, (x-y)/sqrt(2)=%e x=%e y=%e k_min=%e)", */
  /* 	     k,m,(p+m)/sqrt(2),(p-m)/sqrt(2),pnl->k[0]); */
  /*  if (m<pnl->k[0]) m=pnl->k[0]; */

  if (m < pnl->k[0]) {
    *result=0.;
    return _SUCCESS_;
  }


  /* set argument to zero when p is not between the two cut-off values */

  if ((p<pnl->k[0]) || (p>pnl->k[pnl->k_size-1])) {
    *result=0.;
    return _SUCCESS_;
  }

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

    *result *= m*p/2.;

    /* printf("%e %e %e %e\n",k,(p+m)/sqrt(2.),(p-m)/sqrt(2.),(*result)); */

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

    *result *= m*p/2.;

    /* printf("%e %e %e %e\n",k,(p+m)/sqrt(2.),(p-m)/sqrt(2.),(*result)); */

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

    *result *= m*p/2.;

    /* printf("%e %e %e %e\n",k,(p+m)/sqrt(2.),(p-m)/sqrt(2.),(*result)); */

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

    *result *= m*p/2.;

    /* printf("%e %e %e %e\n",k,(p+m)/sqrt(2.),(p-m)/sqrt(2.),(*result)); */

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

    *result *= m*p/2.;

    /* printf("%e %e %e %e\n",k,(p+m)/sqrt(2.),(p-m)/sqrt(2.),(*result)); */

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

    /* if (index_k==103) */
    /*     if (index_k==30 || index_k==pnl->k_size-1) */
    /*        printf("%e %e %e %e %e %e\n",k,p,m,(p+m)/sqrt(2.),(p-m)/sqrt(2.),p*m*(*result)); */ 

    *result *= m*p/2.;    

    /* printf("%e %e %e %e\n",k,(p+m)/sqrt(2.),(p-m)/sqrt(2.),(*result)); */

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

    *result *= m*p/2.;

/*     printf("%e %e %e %e\n",k,(p+m)/sqrt(2.),(p-m)/sqrt(2.),(*result)); */

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

    *result *= m*p/2.;

    /* printf("%e %e %e %e\n",k,(p+m)/sqrt(2.),(p-m)/sqrt(2.),(*result)); */

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

    *result *= m*p/2.;

    /* printf("%e %e %e %e\n",k,(p+m)/sqrt(2.),(p-m)/sqrt(2.),(*result)); */

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

    *result *= m*p/2.;

    /* printf("%e %e %e %e\n",k,(p+m)/sqrt(2.),(p-m)/sqrt(2.),(*result)); */

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

    *result *= m*p/2.;

    /* printf("%e %e %e %e\n",k,(p+m)/sqrt(2.),(p-m)/sqrt(2.),(*result)); */

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

    *result *= m*p/2.;

    /* printf("%e %e %e %e\n",k,(p+m)/sqrt(2.),(p-m)/sqrt(2.),(*result)); */

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

    *result *= m*p/2.;

    /* printf("%e %e %e %e\n",k,(p+m)/sqrt(2.),(p-m)/sqrt(2.),(*result)); */

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

    *result *= m*p/2.;


    /*     if (index_k==30 || index_k==pnl->k_size-1) */
    /*        printf("%e %e %e %e %e %e\n",k,p,m,(p+m)/sqrt(2.),(p-m)/sqrt(2.),p*m*(*result));  */

    /* printf("%e %e %e %e\n",k,(p+m)/sqrt(2.),(p-m)/sqrt(2.),(*result)); */

    return _SUCCESS_;
    break;

    /****************************************/

  default:
    sprintf(pnl->error_message,"%s(L:%d): non valid argument in integrals A of I, %d is not defined\n",__func__,__LINE__,name);
    return _FAILURE_;
    break;
  }

}


/****** this is only for testing and plotting ****************/
int trg_show_borders(
		     struct background * pba,
		     struct primordial * ppm,
		     struct spectra * psp,
		     struct spectra_nl * pnl,
		     enum name_A name,
		     int index_eta,
		     double * result,
		     char *errmsg
		     ){

  double k,x,y,temp;

  k=5.;

  x=k/sqrt(2.);

  for (y=0.; y <= k/sqrt(2.); y += k/sqrt(2.)/1000.) {

    if (((x+y)/sqrt(2) >= pnl->k[0]) && ((x-y)/sqrt(2) >= pnl->k[0])) {
  
      class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2),(x-y)/sqrt(2),index_eta,0,&temp,errmsg),
		 errmsg,
		 pnl->error_message);

      printf("%e %e %e\n",y,temp*(x*x-y*y)/4.,0.);

    }

  }

  printf("\n\n");

  y=k/sqrt(2.);

  for (x=k/sqrt(2.); x<= k*sqrt(2.); x += k/sqrt(2.)/1000.) {

    if (((x+y)/sqrt(2) >= pnl->k[0]) && ((x-y)/sqrt(2) >= pnl->k[0])) {

      class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2),(x-y)/sqrt(2),index_eta,0,&temp,errmsg),
		 errmsg,
		 pnl->error_message);
      
      printf("%e %e %e\n",sqrt(2)*k-x,0.,-temp*(x*x-y*y)/4.);

    }

  }

  printf("\n\n");

  return _SUCCESS_;


}

/***************
 *
 * Integration with the Simpson's method (over y, and simple trapeze for x)
 *
 ***************/

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
  
  int il,index_stop;

  double * partial_sum;
  double * partial_area;
  double sum,area,max;
  double increment_sum,increment_area;

  double local_average_value,previous_average_value,total_average_value;

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

  k_min=pnl->k[0];
  k_max=pnl->k[pnl->k_size-1];

  for(index_k=0; index_k<pnl->k_size; index_k++){
  /* for(index_k=50; index_k<51; index_k++){ */

    k=pnl->k[index_k];

    /*    logstepx=min(1.1,1+0.01/pow(k,2)); */

    logstepx=min(1.1,1+0.01/pow(k,1));
   
    logstepy=logstepx;

    if(index_k<pnl->index_k_L){ /*size under which we pick linear
				  theory, index_k_L defined in
				  trg_init */

      result[index_k+pnl->k_size*index_eta]=0.;
    }

    /* extrapolate high k's behaviour */

    /*    else if(index_k>pnl->k_size-4){ */
    
    /*       result[index_k+pnl->k_size*index_eta]= */
    /* 	((result[pnl->k_size-7+pnl->k_size*index_eta]-result[pnl->k_size-8+pnl->k_size*index_eta])/ */
    /* 	 (pnl->k[pnl->k_size-7]-pnl->k[pnl->k_size-8])) * (pnl->k[index_k]-pnl->k[pnl->k_size-7]) + */
    /* 	result[pnl->k_size-7+pnl->k_size*index_eta]; */
	
    /*     } */
    
    /*     else{  */


    x_size = (int)(log(2.*k_max/k)/log(logstepx)) + 2;

    class_calloc(xx,x_size,sizeof(double),pnl->error_message);

    class_calloc(h_up,x_size,sizeof(double),pnl->error_message);
    class_calloc(h_do,x_size,sizeof(double),pnl->error_message);

    /*       class_calloc(sum_y,x_size,sizeof(double),pnl->error_message); */

    index_x = 0; 

    do {
	
      class_test(index_x >= x_size,
		 pnl->error_message," ");

      xx[index_x] = k/sqrt(2.)*pow(logstepx,index_x); /* set value of x */

      if (xx[index_x] >= k_max*sqrt(2.)) xx[index_x] = k_max*sqrt(2.); /* correct the last value */

      index_x ++;

    } while (xx[index_x-1] < k_max*sqrt(2.));

    if (x_size != index_x) printf("x: %d %d\n", x_size,index_x);

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

    if (y_size != index_y) printf("y: %d %d\n", y_size,index_y);

    y_size = index_y;

    /*       printf("integrate for name=%d, index_k=%d\n",(int)name,index_k); */


    /* compute first h and v lines */

    h_do[0]=0.;
    v_ri[0]=h_do[0];

    for (index_x=1; index_x < x_size; index_x ++) {

      x=xx[index_x];
      y=yy[0];

      if (x <= sqrt(2.)*k_max-y) {
	class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&h_do[index_x],errmsg),
		   errmsg,
		   pnl->error_message);
      }
      else {
	h_do[index_x]=0.;
      }
	
    }

    for (index_y=1; index_y < y_size; index_y ++) {

      x=xx[0];
      y=yy[index_y];

      class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&v_ri[index_y],errmsg),
		 errmsg,
		 pnl->error_message);

    }

    sum = 0.;
    area = 0.;
    max = 0.;

    /********************* loop over L-shaped regions **********************/

    for (il=0; il < y_size-1; il++) {

      /* intialize stop index */ 

      index_stop = x_size-1;

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
      class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&h_do[il+1],errmsg),
		 errmsg,
		 pnl->error_message);
      v_ri[il+1]= h_do[il+1];

      increment_sum = (xx[il+1]-xx[il])*(yy[il]-yy[il+1])*0.25*(h_up[il]+h_up[il+1]+v_le[il+1]+v_ri[il+1]);
      increment_area = (xx[il+1]-xx[il])*(yy[il]-yy[il+1]);

      partial_sum[il] = increment_sum;
      partial_area[il] = increment_area;

      /**************** new points on the horizontal and diagonal, within the square ******************/

      for (index_x=il+1; index_x < y_size-1; index_x ++) {

	/* the point h_up[index_x+1] may have not been calculated at the previous stage; check and calculate */

	if (h_up[index_x+1] == 0.) {

	  x=xx[index_x+1];
	  y=yy[il];
	    
	  if (x <= sqrt(2)*k_max-y) {
	    class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&h_up[index_x+1],errmsg),
		       errmsg,
		       pnl->error_message);
	  } 
	  else {
	    h_up[index_x+1]=0.;
	  }
	    
	}

	/* the point h_do[index_x+1] is new; calculate */

	x=xx[index_x+1];
	y=yy[il+1];

	if (x <= sqrt(2)*k_max-y) {
	  class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&h_do[index_x+1],errmsg),
		     errmsg,
		     pnl->error_message);
	} 
	else {
	  h_do[index_x+1]=0.;
	}

	/* the point v_le[index_x+1] may have not been calculated at the previous stage; check and calculate */

	if (v_le[index_x+1] == 0.) {
	    
	  x=xx[il];
	  y=yy[index_x+1];

	  class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&v_le[index_x+1],errmsg),
		     errmsg,
		     pnl->error_message);
	    
	}

	/* the point v_ri[index_x+1] is new; calculate */

	x=xx[il+1];
	y=yy[index_x+1];

	class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&v_ri[index_x+1],errmsg),
		   errmsg,
		   pnl->error_message);

	/* now integrate on the two new cells */

	increment_sum = (xx[il+1]-xx[il])*(yy[index_x]-yy[index_x+1])*0.25*
	  (v_le[index_x]+v_le[index_x+1]+v_ri[index_x]+v_ri[index_x+1])
	  + (xx[index_x+1]-xx[index_x])*(yy[il]-yy[il+1])*0.25*
	  (h_up[index_x]+h_up[index_x+1]+h_do[index_x]+h_do[index_x+1]);

	increment_area = (xx[il+1]-xx[il])*(yy[index_x]-yy[index_x+1]) 
	  + (xx[index_x+1]-xx[index_x])*(yy[il]-yy[il+1]);

	partial_sum[il] += increment_sum;
	partial_area[il] += increment_area;

	/* if (fabs(increment_sum/increment_area/((sum+partial_sum[il])/(area+partial_area[il]))) < _STOP_INT_) { */
	/* if (fabs(increment_sum/increment_area/(sum+partial_sum[il])) < _STOP_INT_) { */
/* 	  index_stop = index_x+1; /\* will remember where we stoped *\/ */
/* 	  index_x = y_size-1;     /\* to exit this loop *\/ */
/* 	} */

      }

      /***************** new points on the horizontal, beyond the square *******************/

      if (index_stop == x_size-1) {
	for (index_x=y_size-1; index_x < x_size-1; index_x ++) {

	  /* the point h_up[index_x+1] may have not been calculated at the previous stage; check and calculate */

	  if (h_up[index_x+1] == 0.) {
	      
	    x=xx[index_x+1];
	    y=yy[il];
	      
	    if (x <= sqrt(2)*k_max-y) {
	      class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&h_up[index_x+1],errmsg),
			 errmsg,
			 pnl->error_message);
	    } 
	    else {
	      h_up[index_x+1]=0.;
	    }
	    
	  }

	  /* the point h_do[index_x+1] is new; calculate */

	  x=xx[index_x+1];
	  y=yy[il+1];

	  if (x <= sqrt(2)*k_max-y) {
	    class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,index_k,&h_do[index_x+1],errmsg),
		       errmsg,
		       pnl->error_message);
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

	  /* if (fabs(increment_sum/increment_area/((sum+partial_sum[il])/(area+partial_area[il]))) < _STOP_INT_) { */
	  /* if (fabs(increment_sum/(sum+partial_sum[il])) < _STOP_INT_) { */
/* 	    index_stop = index_x+1; /\* will remember where we stoped *\/ */
/* 	    index_x = x_size;     /\* to exit this loop *\/ */
/* 	  } */

	}
      }

      /*************** for non-computed points fill new line/column with zeros **********/

      /* for (index_x = index_stop+1; index_x < x_size; index_x++) */
/* 	h_do[index_x] = 0.; */
      
/*       for (index_y = index_stop+1; index_y < y_size; index_y++) */
/* 	v_ri[index_y] = 0.; */
      
      /* update the total sum with the new L-shaped region */

      sum += partial_sum[il];
      area += partial_area[il];
	
      /* if (il > 0) { */

/* 	local_average_value = partial_sum[il]/partial_area[il]; */
/* 	previous_average_value = partial_sum[il-1]/partial_area[il-1]; */
/* 	total_average_value = sum/area; */

/* 	/\* if (fabs(local_average_value/total_average_value) < _STOP_INT_) *\/ */

/* 	/\* if (fabs(local_average_value/total_average_value) < 1e-2/\\*_STOP_INT_*\\/) *\/ */
/* /\* 	  /\\* 	if ((fabs(local_average_value/total_average_value) < _STOP_INT_) && *\\/ *\/ */
/* /\* 	  /\\* 	    (fabs((local_average_value- previous_average_value) *\\/ *\/ */
/* /\* 	  /\\* 		  /(xx[il]-xx[il-1])/total_average_value) < _STOP_INT_)) *\\/ *\/ */
/* /\* 	  il=y_size-1; *\/ */
	
/* 	/\******** new stopping scheme *******\/ */

/* 	/\* if(fabs(partial_sum[il])>max) *\/ */
/* /\* 	  max=fabs(partial_sum[il]); *\/ */
	
/* /\* 	if(fabs((partial_sum[il]-partial_sum[il-1])/(xx[il]-xx[il-1])) < _STOP_INT_ && *\/ */
/* /\* 	   fabs(partial_sum[il]/max) < _STOP_INT_) *\/ */
/* /\* 	  il=y_size-1; *\/ */
/*       } */
    }

    /*       printf("finished index_k=%d, name=%d\n",index_k,(int)name);  */

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
  double fourpi_over_k;
  int index;  
  int index_plus; /* intermediate variables for time integration */
  
  double Omega_11,Omega_12; /* Definition of the matrix Omega that
			       mixes terms together, two are k
			       indepedant and two dependant*/
  double *Omega_21,*Omega_22;

  double *a0,*a11,*a12,*a13,*a21,*a22,*a23,*a3;
  double *b0,*b11,*b12,*b21,*b22,*b3; /* Definition of the variables
					 1ijk,lmn and 2 ijk,lmn that
					 replace the Bispectra
					 variables*/

  double *A0,*A11,*A12,*A13,*A21,*A22,*A23,*A3;
  double *B0,*B11,*B12,*B21,*B22,*B3; /* Definition of the variables
					 Aijk,lmn */

  double * pvecback_nl;

  pnl->spectra_nl_verbose=0;
  pnl->mode=1; /* 0 is linear evolution, 1 one loop and 2 full trg */

  if (pnl->spectra_nl_verbose > 0)
    printf("Computing non-linear spectra with TRG method\n");

  class_alloc(pvecback_nl,pba->bg_size*sizeof(double),pnl->error_message);

  /* define initial eta and redshift */
  a_ini = 2e-2;   /* gives a starting redshift of 49 */
  pnl->z_ini = pba->a_today/a_ini - 1.;

  /* define eta_max, where eta=log(a/a_ini) */
  eta_max = log(pba->a_today/a_ini);

  /* define size and step for integration in eta */
  pnl->eta_size = 100; /* to calculate fast */
  pnl->eta_step = (eta_max)/(pnl->eta_size-1);
  eta_step = pnl->eta_step;

  /* at any time, a = a_ini * exp(eta) and z=exp(eta)-1*/
  
  /* fill array of k values  */

  /* first define the total length, to reach the k_max bound
     for the integrator */
  double k_max;
  double k_L;
  double k_min;
  /*   double k_2; */
  double logstepk;

  k_max=20;    /**< PRECISION PARAMETER */ /* not above k_scalar_max*h in test_trg */
  k_L=1e-3;        /**< PRECISION PARAMETER */
  k_min=1e-4;      /**< PRECISION PARAMETER */

  logstepk=1.2;    /**< PRECISION PARAMETER */ /* with 1.01 we map
						  finely the
						  oscillations,
						  though very much
						  time consuming */

  /* find total number of k values in the module */
  index_k=0;
  temp_k=k_min;
  while(temp_k<k_max) {
    index_k++;
    temp_k=k_min*pow(logstepk,index_k);
  }
  pnl->k_size=index_k+1;

  class_alloc(pnl->k,pnl->k_size*sizeof(double),pnl->error_message);

  /* then fill in the values of k, while determining the index below
     which we use only the linear theory */

  temp_k=0;

  for(index_k=0; index_k<pnl->k_size-1; index_k++){
    pnl->k[index_k]=k_min*pow(logstepk,index_k);
    if( (pnl->k[index_k] > k_L*pba->h) && temp_k==0){
      pnl->index_k_L=index_k;
      temp_k++;
    }
    /*     if( (pnl->k[index_k] > 1.5*pba->h) && temp_k==1 ){ */
    /*       pnl->k_size_eff=index_k+1; */
    /*       k_2=pnl->k[index_k]; */
    /*       temp_k++; */
    /*     } */
  }
  pnl->k[pnl->k_size-1]=k_max;

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

  
  
  /*   return _SUCCESS_; */
  
  /* Definition of 1_0, 1_11,(here a0, a11,...) etc, and 2_0, 2_11,
     (here b0,b11,...) etc.. and initialization (directly with calloc
     for assuming no initial non gaussianity in the spectrum) */

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

  /* Definition of the A_121,122's that appear in time-evolution of
     the 1's and 2's, and initialization with pk_nl, p_12,
     p_22. Convention is taken for A_121,any = A any and A_222,any = B any */

  class_calloc(A0 ,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(A11,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(A12,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(A13,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(A21,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(A22,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(A23,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(A3 ,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);

  class_calloc(B0 ,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(B11,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(B12,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(B21,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(B22,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(B3 ,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);

  /********************
   *
   * Integration at eta=0 with n_xy number of steps for integration
   * over x and n_xy over y.
   *
   ********************/
  double temp1,temp2,temp3;

  /**** testing zone ***/
  /*   class_call(trg_show_borders(pba,ppm,psp,pnl,_A0_,0,A0,pnl->error_message), */
  /* 	     pnl->error_message, */
  /* 	     pnl->error_message); */
  /*   class_call(trg_show_borders(pba,ppm,psp,pnl,_A11_,0,A0,pnl->error_message), */
  /* 	     pnl->error_message, */
  /* 	     pnl->error_message); */
  /*   class_call(trg_show_borders(pba,ppm,psp,pnl,_A12_,0,A0,pnl->error_message), */
  /* 	     pnl->error_message, */
  /* 	     pnl->error_message); */
  /*   class_call(trg_show_borders(pba,ppm,psp,pnl,_A13_,0,A0,pnl->error_message), */
  /* 	     pnl->error_message, */
  /* 	     pnl->error_message); */
  /*   class_call(trg_show_borders(pba,ppm,psp,pnl,_A21_,0,A0,pnl->error_message), */
  /* 	     pnl->error_message, */
  /* 	     pnl->error_message); */
  /*   class_call(trg_show_borders(pba,ppm,psp,pnl,_A22_,0,A0,pnl->error_message), */
  /* 	     pnl->error_message, */
  /* 	     pnl->error_message); */
  /*   class_call(trg_show_borders(pba,ppm,psp,pnl,_A23_,0,A0,pnl->error_message), */
  /* 	     pnl->error_message, */
  /* 	     pnl->error_message); */
  /*   class_call(trg_show_borders(pba,ppm,psp,pnl,_A3_,0,A0,pnl->error_message), */
  /* 	     pnl->error_message, */
  /* 	     pnl->error_message); */
  /*   class_call(trg_show_borders(pba,ppm,psp,pnl,_B0_,0,A0,pnl->error_message), */
  /* 	     pnl->error_message, */
  /* 	     pnl->error_message); */
  /*   class_call(trg_show_borders(pba,ppm,psp,pnl,_B11_,0,A0,pnl->error_message), */
  /* 	     pnl->error_message, */
  /* 	     pnl->error_message); */
  /*   class_call(trg_show_borders(pba,ppm,psp,pnl,_B12_,0,A0,pnl->error_message), */
  /* 	     pnl->error_message, */
  /* 	     pnl->error_message); */
  /*   class_call(trg_show_borders(pba,ppm,psp,pnl,_B21_,0,A0,pnl->error_message), */
  /* 	     pnl->error_message, */
  /* 	     pnl->error_message); */
  /*   class_call(trg_show_borders(pba,ppm,psp,pnl,_B22_,0,A0,pnl->error_message), */
  /* 	     pnl->error_message, */
  /* 	     pnl->error_message); */
  /*   class_call(trg_show_borders(pba,ppm,psp,pnl,_B3_,0,A0,pnl->error_message), */
  /* 	     pnl->error_message, */
  /* 	     pnl->error_message); */
  /*   return _SUCCESS_; */


  if (pnl->spectra_nl_verbose > 0)
    printf(" -> initialisation\n");
  
  if(pnl->mode > 0){

    class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_A0_,0,A0,pnl->error_message),
      	       pnl->error_message,
      	       pnl->error_message);
    
/*     printf("\n\n"); */

    class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_A11_,0,A11,pnl->error_message),
      	       pnl->error_message,
      	       pnl->error_message);
    
/*      printf("\n\n");  */

    class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_A12_,0,A12,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    

/*      printf("\n\n");  */

    class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_A13_,0,A13,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    
/*      printf("\n\n");  */

    class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_A21_,0,A21,pnl->error_message),
      	       pnl->error_message,
      	       pnl->error_message);
    
/*      printf("\n\n");  */

    class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_A22_,0,A22,pnl->error_message),
      	       pnl->error_message,
      	       pnl->error_message);
    
/*      printf("\n\n");  */

    class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_A23_,0,A23,pnl->error_message),
      	       pnl->error_message,
      	       pnl->error_message);
    
/*      printf("\n\n");  */

    class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_A3_,0,A3,pnl->error_message),
      	       pnl->error_message,
      	       pnl->error_message);
    
/*      printf("\n\n");  */

    class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_B0_,0,B0,pnl->error_message),
      	       pnl->error_message,
      	       pnl->error_message);
    
/*      printf("\n\n");  */

    class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_B11_,0,B11,pnl->error_message),
      	       pnl->error_message,
      	       pnl->error_message);
    
/*      printf("\n\n");  */

    class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_B12_,0,B12,pnl->error_message),
      	       pnl->error_message,
      	       pnl->error_message);

/*      printf("\n\n");  */

    class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_B21_,0,B21,pnl->error_message),
      	       pnl->error_message,
      	       pnl->error_message);
    
     printf("\n\n"); 

    class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_B22_,0,B22,pnl->error_message),
      	       pnl->error_message,
      	       pnl->error_message);

/*      printf("\n\n");  */

    class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_B3_,0,B3,pnl->error_message),
      	       pnl->error_message,
      	       pnl->error_message);
  }

  class_open(nl_spectra,"output/nl_ini.dat","wr",pnl->error_message);

  for(index_k=0; index_k<pnl->k_size; index_k++){ fprintf(nl_spectra,
							  "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
							  pnl->k[index_k],
							  A0[index_k],A11[index_k],
							  A12[index_k],A13[index_k],A21[index_k],
							  A22[index_k],A23[index_k],A3[index_k],
							  B0[index_k],B11[index_k],
							  B12[index_k],B21[index_k],
							  B22[index_k],B3[index_k],
							  pnl->pk_nl[index_k],pnl->p_12_nl[index_k],
							  pnl->p_22_nl[index_k]);}
  
  fclose(nl_spectra);

  return _SUCCESS_;


  /********************
   * Now we calculate the time evolution with a very simple integrator
   ********************/

  class_open(nl_spectra,"output/nl_spectra.dat","wr",pnl->error_message)

    if (pnl->spectra_nl_verbose > 0){
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
					  +2*exp_eta*A0[index])
	+ a0[index];

      a11[index_plus]        = eta_step *(
					  -a11[index]*(2*Omega_22[index_eta-1]+Omega_11)
					  -Omega_12*a0[index]
					  -Omega_21[index_eta-1]*(a22[index]+a23[index])
					  +2*exp_eta*A11[index])
	+ a11[index];

      a12[index_plus]        = eta_step *(
					  -a12[index]*(2*Omega_22[index_eta-1]+Omega_11)
					  -Omega_21[index_eta-1]*(a23[index]+a21[index])
					  -Omega_12*a0[index]
					  +2*exp_eta*A12[index])
	+ a12[index];

      a13[index_plus]        = eta_step *(
					  -a13[index]*(2*Omega_22[index_eta-1]+Omega_11)
					  -Omega_12*a0[index]
					  -Omega_21[index_eta-1]*(a22[index]+a21[index])
					  +2*exp_eta*A13[index])
	+ a13[index];

      a21[index_plus]        = eta_step *(
					  -a21[index]*(2*Omega_11+Omega_22[index_eta-1])
					  -Omega_12*(a12[index]+a13[index])
					  -Omega_21[index_eta-1]*a3[index]
					  +2*exp_eta*A21[index])
	+ a21[index];

      a22[index_plus]        = eta_step *(
					  -a22[index]*(2*Omega_11+Omega_22[index_eta-1])
					  -Omega_12*(a13[index]+a11[index])
					  -Omega_21[index_eta-1]*a3[index]
					  +2*exp_eta*A22[index])
	+ a22[index];

      a23[index_plus]        = eta_step *(
					  -a23[index]*(2*Omega_11+Omega_22[index_eta-1])
					  -Omega_12*(a12[index]+a11[index])
					  -Omega_21[index_eta-1]*a3[index]
					  +2*exp_eta*A23[index])
	+ a23[index];

      a3[index_plus]         = eta_step *(
					  -a3[index]*3*Omega_11
					  -Omega_12*(a22[index]+a21[index]+a23[index])
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

	class_call(trg_p12_at_k(pba,ppm,psp,pnl,index_eta,0,pnl->k[index_k],&pnl->p_12[index_k+pnl->k_size*index_eta]),
		   pnl->error_message,
		   pnl->error_message);

	class_call(trg_p22_at_k(pba,ppm,psp,pnl,index_eta,0,pnl->k[index_k],&pnl->p_22[index_k+pnl->k_size*index_eta]),
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
     * Update of A's and B's function at the new index_eta which needs the interpolation of P's
     **********/
    
    if(pnl->mode>0){
      class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_A0_,index_eta,A0,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);
      
      class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_A11_,index_eta,A11,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);
      
      class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_A12_,index_eta,A12,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);
      
      class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_A13_,index_eta,A13,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);

      class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_A21_,index_eta,A21,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);

      class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_A22_,index_eta,A22,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);
    
      class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_A23_,index_eta,A23,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);

      class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_A3_,index_eta,A3,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);

      class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_B0_,index_eta,B0,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);

      class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_B11_,index_eta,B11,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);

      class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_B12_,index_eta,B12,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);

      class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_B21_,index_eta,B21,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);
      
      class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_B22_,index_eta,B22,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);

      class_call(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,_B3_,index_eta,B3,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);
    }
    
    printf("    %2.1f%% done\n",100.*index_eta/(pnl->eta_size-1.));

    /* to plot some values, directly on screen in case of divergence and the program does not execute to the end */

    /* if(index_eta==2 || index_eta==20 || index_eta==30){  for(index_k=0; index_k<pnl->k_size; index_k++){ printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", */
    /* 									       pnl->k[index_k], */
    /* 									       A0[index_k+pnl->k_size*index_eta],A11[index_k+pnl->k_size*index_eta], */
    /* 									       A12[index_k+pnl->k_size*index_eta],A13[index_k+pnl->k_size*index_eta],A21[index_k+pnl->k_size*index_eta], */
    /* 									       A22[index_k+pnl->k_size*index_eta],A23[index_k+pnl->k_size*index_eta],A3[index_k+pnl->k_size*index_eta], */
    /* 									       B0[index_k+pnl->k_size*index_eta],B11[index_k+pnl->k_size*index_eta], */
    /* 									       B12[index_k+pnl->k_size*index_eta],B21[index_k+pnl->k_size*index_eta], */
    /* 									       B22[index_k+pnl->k_size*index_eta],B3[index_k+pnl->k_size*index_eta], */
    /* 									       pnl->pk_nl[index_k+pnl->k_size*index_eta],pnl->p_12_nl[index_k+pnl->k_size*index_eta], */
    /* 														pnl->p_22_nl[index_k+pnl->k_size*index_eta]);} */
    /*     } */
    
    
    time_2=time(NULL);
    if(index_eta==9){
      printf("elapsed time after ten loops : %2.f s\n",difftime(time_2, time_1));
      printf("estimated remaining : %3.1f min\n",difftime(time_2,time_1)*(pnl->eta_size-2)/60/10);
    }
    if(isnan(pnl->pk_nl[50+pnl->k_size*index_eta])!=0){
      printf("ca marche pas !!!nan!!!\n");
    }
  }

  printf("Done in %2.f min !\n",difftime(time_2,time_1)/60);

  /* printing on screen */

  fprintf(nl_spectra,"## k (h/Mpc) P_NL ([Mpc/h]^3) P_L ([Mpc/h]^3) \n");
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
  
  fprintf(nl_spectra,"\n\n##for %d values of k\n##z is equal to %f\n## k\tpk_nl\tpk_lin at last eta(today)\n",pnl->k_size,pnl->z[pnl->eta_size-1]);

  for(index_k=0; index_k<pnl->k_size; index_k++){
    class_call(
	       spectra_pk_at_k_and_z(pba,ppm,psp,0,index_ic,pnl->k[index_k],pnl->z[pnl->eta_size-1],&temp1),
	       psp->error_message,
	       pnl->error_message);

    fprintf(nl_spectra,"%e\t%e\t%e\n",pnl->k[index_k]/pba->h,pow(pba->h,3)*pnl->pk_nl[index_k+(pnl->k_size*(pnl->eta_size-1))]*exp(pnl->eta[pnl->eta_size-1]*2),pow(pba->h,3)*temp1);
  }

  fclose(nl_spectra);
  
  free(B3);
  free(B22);
  free(B21);
  free(B12);
  free(B11);
  free(B0);

  free(A3);
  free(A23);
  free(A22);
  free(A21);
  free(A13);
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
