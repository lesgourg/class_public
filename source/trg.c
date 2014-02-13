/** @file trg.c Documented Time Renormalization Group module
 *
 * Benjamin Audren, 06.06.2011
 *
 * Computes the non linear matter spectra P_k with Time
 * Renormalization Group method or one-loop given the linear power
 * spectrum at z_ini. No initial non-gaussianity assumed. The matter
 * density perturbation is here understood as the barycenter of the
 * baryon density and cdm density.  Hereafter, 'matter' is used in
 * this sense.  To select the mode of execution, please change the
 * value of the "non linear" parameter in your .ini file, respectively
 * to 'trg' or 'one-loop'. For testing purpose, you might want to use
 * the value 'test-linear', which will then output the linear spectra
 * out of the equations.
 *
 * The logic is the following :
 *
 * - in a first step, recover all relevant information from all other modules.
 *
 * - in a second step, initialize all relevant quantities, i.e. the matter power
 *   spectra and the I's variables (that contain the matter bispectra).
 *
 * - in a third step, compute the A variables (complicated functions of the matter
 *   power spectra), at first time.
 *
 * - the final step consists in a loop, where the matter power spectra
 *   and I's are computed step by step in time, by integrating the
 *   system of equation with the predictor corrector algorithm
 *   described in the released paper. With the new spectra, the code
 *   computes the new A's, and continues the loop.
 *
 * - the three quantities pnl->p_11_nl, pnl->p_12_nl, pnl->p_22_nl
 *   corresponding to, resp, the density, velocity and cross
 *   correlation are written after the call of trg_init() in the
 *   nonlinear.c module, ready for plotting or further use.
 **/

#include "background.h"
#include "thermodynamics.h"
#include "perturbations.h"
#include "transfer.h"
#include "primordial.h"
#include "spectra.h"
#include "trg.h"
#include "time.h"
#ifdef _OPENMP
#include "omp.h"
#endif


/**
 * Gamma function 121, for any set of wavevectors (k,p,q)
 *
 * Encodes the non linear behaviour of the evolution
 *
 * @param k       Input: wave-vector
 * @param p       Input: wave-vector
 * @param q       Input: wave-vector
 * @param result  Output: computed the rational interaction function
 */

int trg_gamma_121(
		  double  k,
		  double  p,
		  double  q,
		  double  *result
		  ){
  *result =  (1./ (4*p*p)) * (- q*q + k*k + p*p);
  return _SUCCESS_;
}

/**
 * Gamma function 222, for any set of wavevectors (k,p,q)
 *
 * Encodes the non linear behaviour of the evolution
 *
 * @param k       Input: wave-vector
 * @param p       Input: wave-vector
 * @param q       Input: wave-vector
 * @param result  Output: computed the rational interaction function
 */


int trg_gamma_222(
		  double k,
		  double p,
		  double q,
		  double * result,
		  char   * errmsg
		  ){

  class_test(pow(k,2)==pow(p,2)+pow(q,2),
	     errmsg,
	     "Function trg_gamma_222 is dividing by zero, check boundaries in trg_integrate_");

  *result=k*k/(4*q*q*p*p) * (k*k - p*p - q*q);
  return _SUCCESS_;
}

/**
 * Interpolates-extrapolates p_ab for any k different from pnl->k
 *
 * It is just a renaming of the already existing structure
 * array_interpolate_extrapolate_logspline_loglinear_one_column for simple
 * convenient reasons : this routine is called many times in the trg_arg_*
 * functions, with many redundancy in the arguments.
 *
 * @param pnl       Input: pointer to spectra_nl structure
 * @param p_ab      Input: table of values containing any spectrum, already allocated, and supposed filled until [any_index_k+pnl->k_size*index_eta]
 * @param ddp_ab    Input: table of values containing the second derivative of any spectrum, already allocated, and supposed filled until [any_index_k+pnl->k_size*index_eta]
 * @param index_eta Input: index of the time evolution
 * @param any_k     Input: desired value of the wavenumber to calculate the spectrum. If it is already part of pnl->k, then the answer is given. If any_k < k_max, the result is interpolated. If any_k > k_max, the result is extrapolated
 * @param result    Output: P_ab at the given any_k for the desired eta (pnl->eta[index_eta])
 * @return the error status
 */

int trg_p_ab_at_any_k(
		      struct spectra_nl * pnl,
		      double *p_ab,
		      double *ddp_ab,
		      int    index_eta,
		      double any_k,
		      double *result
		      ){

  class_call(array_interpolate_extrapolate_logspline_loglinear_one_column(pnl->k,
									  pnl->k_size,
									  pnl->k_size-pnl->double_escape*2*index_eta,
									  p_ab,
									  pnl->eta_size,
									  index_eta,
									  ddp_ab,
									  any_k,
									  result,
									  pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);

  return _SUCCESS_;
}

/**
 * Argument of the integral defining AA[name] for the TRG method
 *
 * Returns the argument (symmetrized in p,m) of the desired integral,
 * called by the 'name' argument, at the desired wave-number values,
 * and time pnl->eta[index_eta]. The names 'p' and 'm' stands resp. for
 * 'plus' and 'minus', standing themselves resp. for
 * '(x + y)/sqrt(2)' and '(x - y)/sqrt(2)'
 *
 * This trg version uses the full non-linear spectrum computed at each time
 *
 * @param pnl       Input: pointer to spectra_nl structure
 * @param name      Input: select the desired A function to get the argument from, from _A0_ to _B3_
 * @param k         Input: wave-vector
 * @param p         Input: wave-vector corresponding to (x+y)/sqrt(2)
 * @param m         Input: wave-vector corresponding to (x-y)/sqrt(2)
 * @param index_eta Input: index of time evolution
 * @param result    Output: argument of the integral for the given wave numbers and at the given time
 * @return the error status
 */



int trg_A_arg_trg(
	      struct spectra_nl * pnl,
	      enum name_A name,
	      double k,
	      double p,
	      double m,
	      int index_eta,
	      double * result
	      ){

  /** - define local variables */

  /* shorter names for the spectra */

  double p_11k,p_11m,p_11p;
  double p_12k,p_12m,p_12p;
  double p_22k,p_22m,p_22p;

  /* shorter names for the interaction functions */

  double gamma1_kpm,gamma1_pkm,gamma1_mkp,gamma1_mpk,gamma1_pmk,gamma1_kmp;
  double gamma2_kpm,gamma2_pkm,gamma2_mkp,gamma2_mpk,gamma2_pmk,gamma2_kmp;

  /** - take care of the cases where one of the argument is below k_min */

  if (m<pnl->k[0] || p<pnl->k[0]){
    *result=0.;
    return _SUCCESS_;
  }

  /** - for each name (_A0_ to _B3_), assign all and only the relevant parameters
      to the local variables, and return the correct argument. */

  /* each case is separated for computation time reason, though it makes the function
     harder to read. For precision : each 'case' statement starts by calling P_11, P_12, P_22
     then gamma_222, then gamma_121. It then computes the function, and finally multiplies
     by the Fourier factor and the p.m factor (x^2-y^2) */

  switch(name){
  case _A0_:
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,k,&p_22k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,p,&p_22p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,m,&p_22m),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,pnl->error_message),
	       pnl->error_message,
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
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,k,&p_12k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,p,&p_12p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,m,&p_12m),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,k,&p_22k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,p,&p_22p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,m,&p_22m),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,pnl->error_message),
	       pnl->error_message,
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
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,k,&p_12k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,p,&p_12p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,m,&p_12m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,k,&p_22k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,p,&p_22p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,m,&p_22m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,pnl->error_message),
	       pnl->error_message,
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


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,k,&p_12k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,p,&p_12p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,m,&p_12m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,k,&p_22k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,p,&p_22p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,m,&p_22m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,pnl->error_message),
	       pnl->error_message,
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


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11_nl,pnl->ddp_11_nl,index_eta,p,&p_11p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11_nl,pnl->ddp_11_nl,index_eta,m,&p_11m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,k,&p_12k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,p,&p_12p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,m,&p_12m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,k,&p_22k),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,pnl->error_message),
	       pnl->error_message,
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
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11_nl,pnl->ddp_11_nl,index_eta,k,&p_11k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11_nl,pnl->ddp_11_nl,index_eta,p,&p_11p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11_nl,pnl->ddp_11_nl,index_eta,m,&p_11m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,k,&p_12k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,p,&p_12p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,m,&p_12m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,p,&p_22p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,m,&p_22m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,pnl->error_message),
	       pnl->error_message,
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
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11_nl,pnl->ddp_11_nl,index_eta,k,&p_11k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11_nl,pnl->ddp_11_nl,index_eta,p,&p_11p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11_nl,pnl->ddp_11_nl,index_eta,m,&p_11m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,k,&p_12k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,p,&p_12p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,m,&p_12m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,p,&p_22p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,m,&p_22m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,pnl->error_message),
	       pnl->error_message,
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
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11_nl,pnl->ddp_11_nl,index_eta,k,&p_11k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11_nl,pnl->ddp_11_nl,index_eta,p,&p_11p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11_nl,pnl->ddp_11_nl,index_eta,m,&p_11m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,k,&p_12k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,p,&p_12p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,m,&p_12m),
	       pnl->error_message,
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
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11_nl,pnl->ddp_11_nl,index_eta,k,&p_11k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11_nl,pnl->ddp_11_nl,index_eta,p,&p_11p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11_nl,pnl->ddp_11_nl,index_eta,m,&p_11m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,k,&p_12k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,p,&p_12p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,m,&p_12m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,pnl->error_message),
	       pnl->error_message,
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
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11_nl,pnl->ddp_11_nl,index_eta,k,&p_11k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11_nl,pnl->ddp_11_nl,index_eta,p,&p_11p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11_nl,pnl->ddp_11_nl,index_eta,m,&p_11m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,k,&p_12k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,p,&p_12p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,m,&p_12m),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,k,&p_22k),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,pnl->error_message),
	       pnl->error_message,
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


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11_nl,pnl->ddp_11_nl,index_eta,k,&p_11k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11_nl,pnl->ddp_11_nl,index_eta,p,&p_11p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11_nl,pnl->ddp_11_nl,index_eta,m,&p_11m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,k,&p_12k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,p,&p_12p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,m,&p_12m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,p,&p_22p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,m,&p_22m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,pnl->error_message),
	       pnl->error_message,
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
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,k,&p_12k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,p,&p_12p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,m,&p_12m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,k,&p_22k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,p,&p_22p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,m,&p_22m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,pnl->error_message),
	       pnl->error_message,
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
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,k,&p_12k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,p,&p_12p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12_nl,pnl->ddp_12_nl,index_eta,m,&p_12m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,k,&p_22k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,p,&p_22p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,m,&p_22m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,pnl->error_message),
	       pnl->error_message,
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
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,k,&p_22k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,p,&p_22p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22_nl,pnl->ddp_22_nl,index_eta,m,&p_22m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,pnl->error_message),
	       pnl->error_message,
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

/**
 * Argument of the integral defining AA[name] for the TRG method
 *
 * Returns the argument (symmetrized in p,m) of the desired integral,
 * called by the 'name' argument, at the desired wave-number values,
 * and time pnl->eta[index_eta]. The names 'p' and 'm' stands resp. for
 * 'plus' and 'minus', standing themselves resp. for
 * '(x + y)/sqrt(2)' and '(x - y)/sqrt(2)'
 *
 * This one-loop version uses the linear spectrum computed initially, and
 * hence do not need any index_eta. It is indeed called only once, before
 * the time evolution loop, and is copied from time step to time step to
 * keep the same structure than for the trg computation.
 *
 * @param pnl       Input: pointer to spectra_nl structure
 * @param name      Input: select the desired A function to get the argument from, from _A0_ to _B3_
 * @param k         Input: wave-vector
 * @param p         Input: wave-vector corresponding to (x+y)/sqrt(2)
 * @param m         Input: wave-vector corresponding to (x-y)/sqrt(2)
 * @param result    Output: argument of the integral for the given wave numbers and at the given time
 * @return the error status
 */

int trg_A_arg_one_loop(
		       struct spectra_nl * pnl,
		       enum name_A name,
		       double k,
		       double p,
		       double m,
		       double * result
		       ){

  /** - define local variables */

  /* shorter names for the spectra */

  double p_11k,p_11m,p_11p;
  double p_12k,p_12m,p_12p;
  double p_22k,p_22m,p_22p;

  /* shorter names for the interaction functions */

  double gamma1_kpm,gamma1_pkm,gamma1_mkp,gamma1_mpk,gamma1_pmk,gamma1_kmp;
  double gamma2_kpm,gamma2_pkm,gamma2_mkp,gamma2_mpk,gamma2_pmk,gamma2_kmp;

  /** - take care of the cases where one of the argument is below k_min */

  if (m<pnl->k[0] || p<pnl->k[0]){
    *result=0.;
    return _SUCCESS_;
  }

  /** - for each name (_A0_ to _B3_), assign all and only the relevant parameters
      to the local variables, and return the correct argument. */

  /* each case is separated for computation time reason, though it makes the function
     harder to read. For precision : each 'case' statement starts by calling P_11, P_12, P_22
     then gamma_222, then gamma_121. It then computes the function, and finally multiplies
     by the Fourier factor and the p.m factor (x^2-y^2) */


  switch(name){
  case _A0_:
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,k,&p_22k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,p,&p_22p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,m,&p_22m),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,pnl->error_message),
	       pnl->error_message,
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
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,k,&p_12k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,p,&p_12p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,m,&p_12m),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,k,&p_22k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,p,&p_22p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,m,&p_22m),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,pnl->error_message),
	       pnl->error_message,
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
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,k,&p_12k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,p,&p_12p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,m,&p_12m),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,k,&p_22k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,p,&p_22p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,m,&p_22m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,pnl->error_message),
	       pnl->error_message,
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


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,k,&p_12k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,p,&p_12p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,m,&p_12m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,k,&p_22k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,p,&p_22p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,m,&p_22m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,pnl->error_message),
	       pnl->error_message,
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


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,p,&p_11p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,m,&p_11m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,k,&p_12k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,p,&p_12p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,m,&p_12m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,k,&p_22k),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,pnl->error_message),
	       pnl->error_message,
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
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,k,&p_11k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,p,&p_11p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,m,&p_11m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,k,&p_12k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,p,&p_12p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,m,&p_12m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,p,&p_22p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,m,&p_22m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,pnl->error_message),
	       pnl->error_message,
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
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,k,&p_11k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,p,&p_11p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,m,&p_11m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,k,&p_12k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,p,&p_12p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,m,&p_12m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,p,&p_22p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,m,&p_22m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,pnl->error_message),
	       pnl->error_message,
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
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,k,&p_11k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,p,&p_11p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,m,&p_11m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,k,&p_12k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,p,&p_12p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,m,&p_12m),
	       pnl->error_message,
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
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,k,&p_11k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,p,&p_11p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,m,&p_11m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,k,&p_12k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,p,&p_12p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,m,&p_12m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,pnl->error_message),
	       pnl->error_message,
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
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,k,&p_11k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,p,&p_11p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,m,&p_11m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,k,&p_12k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,p,&p_12p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,m,&p_12m),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,k,&p_22k),
	       pnl->error_message,
	       pnl->error_message);

    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,pnl->error_message),
	       pnl->error_message,
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


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,k,&p_11k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,p,&p_11p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,m,&p_11m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,k,&p_12k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,p,&p_12p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,m,&p_12m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,p,&p_22p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,m,&p_22m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,pnl->error_message),
	       pnl->error_message,
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
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,k,&p_12k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,p,&p_12p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,m,&p_12m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,k,&p_22k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,p,&p_22p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,m,&p_22m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,pnl->error_message),
	       pnl->error_message,
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
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,k,&p_12k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,p,&p_12p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,m,&p_12m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,k,&p_22k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,p,&p_22p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,m,&p_22m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,pnl->error_message),
	       pnl->error_message,
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
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,k,&p_22k),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,p,&p_22p),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,m,&p_22m),
	       pnl->error_message,
	       pnl->error_message);


    class_call(trg_gamma_222(k,p,m,&gamma2_kpm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(k,m,p,&gamma2_kmp,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,k,m,&gamma2_pkm,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(p,m,k,&gamma2_pmk,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,p,k,&gamma2_mpk,pnl->error_message),
	       pnl->error_message,
	       pnl->error_message);
    class_call(trg_gamma_222(m,k,p,&gamma2_mkp,pnl->error_message),
	       pnl->error_message,
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




/**
 * Performs the integration of the arguments given by the trg_A_arg functions
 *
 * The integration pattern used here is rather complicated. It is described step
 * by step inside the function, but here is summarized the ideas. For each pnl->k,
 *
 * - in a first step, determine the boundaries of the integration on x and y.
 *
 * - in a second step, perform the integration over the computed domain. To smooth out
 *   the function, the mean value of the argument on a square is integrated, and not
 *   simply the argument. Cases are thus defined on the grid domain.
 *
 * - to avoid computing several times the function over the same point, local variables
 *   are used, and the computation takes place on two adjacent horizontal lines (h_up and h_do(wn))
 *   and two adjacent vertical lines (h_le(ft) and h_ri(ght)). Every step, the lines are
 *   switched one step right and down, moving all computed values from h_do and l_ri to h_up and l_le,
 *   leaving the rest to be computed.
 *
 * - this integration scheme consists thus in integrating over L-shaped region, starting
 *   from the upper left corner. Since the region is not exactly squared, the last lower-right
 *   corner is integrated line per line at the end. This scheme has been designed to
 *   better follow the compensations in the computation:
 *     - the majority of the features of the argument takes place in the first "L" region, with
 *     opposite signs on the horizontal and the vertical
 *     - though, the integration has to be performed over all the region because small contributions
 *     continue to add even in the lower-right corner
 *
 * - a large k_max is defined in this function to replace the 'infinity' upper limit of the integral
 *   over y. The value has been defined to get the convergence of the integral,
 *   while keeping the computing time reasonnable.
 *
 * - finally, a 'double escape' scheme has been implemented to get rid of numerical divergences.
 *   It concretely consists in forgetting a fixed number of points (pnl->double_escape) at each half-step
 *   in time, meaning every time new P's are calculated and every time new A's are calculated.
 *
 *
 * @param pba       Input: pointer to background structure
 * @param ppm       Input: pointer to primordial structure
 * @param psp       Input: pointer to spectra structure
 * @param pnl       Input: pointer to spectra_nl structure
 * @param name      Input: select the desired A function to compute, from _A0_ to _B3_
 * @param index_eta Input: index of time evolution
 * @param result    Output: returns the computed integral corresponding to the whole function A[name]
 * @return the error status
 */

int trg_integrate_xy_at_eta(
			    struct background * pba,
			    struct primordial * ppm,
			    struct spectra * psp,
			    struct spectra_nl * pnl,
			    enum name_A name,
			    int index_eta,
			    double * result
			    ){

 /** - define local variables */

  int index_k;
  double k;
  double k_max;

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
  double sum,area;
  double increment_sum,increment_area;

  /** - set a value for the upper bound in x, hard-coded. */

  k_max=3000.;

  /** - enter the loop over k-values still under the 'double_escape'
        condition */

  for(index_k=0; index_k<pnl->k_size-pnl->double_escape*(2*index_eta+1); index_k++){

    k=pnl->k[index_k];

    /** - define the steps on x and y */

    /* for large k's, one must reduce the steps to be accurate
       enough. These are precision parameters hard-coded, maybe there
       is some execution time to save here, but not much. */

    logstepx=MIN(1.07,1+0.01/pow(k,1));

    /* however, a pnl->logstepx_min is defined to restrain the
       computing time. This is the key element controling the speed of
       the computation. For information, with eta_size=100,
       pnl->logstepx_min = 1.07 : 70 minutes, pnl->logstepx_min =
       1.02: 4 hours. Can be changed in trg.ini */

    if(logstepx< pnl->logstepx_min)  logstepx= pnl->logstepx_min;

    logstepy=logstepx;

    /** A) deal with the case k < k_linear. Can be changed in trg.ini */

    /* for sufficiently small k's, the modes are supposed to stay
       linear at all time during the evolution, their A function is
       thus equal to zero */
    if(index_k<pnl->index_k_L){

      result[index_k+pnl->k_size*index_eta]=0.;

    }

    /** B) deal with the case k > k_linear : full calculation */

    else{

       /* define the integration domain over x */

      x_size = (int)(log(2.*k_max/k)/log(logstepx)) + 2;

      class_calloc(xx,x_size,sizeof(double),pnl->error_message);

      class_calloc(h_up,x_size,sizeof(double),pnl->error_message);
      class_calloc(h_do,x_size,sizeof(double),pnl->error_message);

      index_x = 0;

      /* affect the values of x */

      do {

	class_test(index_x >= x_size,
		   pnl->error_message," ");

	xx[index_x] = k/sqrt(2.)*pow(logstepx,index_x);

	if (xx[index_x] >= k_max*sqrt(2.)) xx[index_x] = k_max*sqrt(2.); /* correct the last value */

	index_x ++;

      } while (xx[index_x-1] < k_max*sqrt(2.));

      x_size = index_x; /* to enforce it and avoid problems */



      /* define the integration domain over y */

      y_size = (int)(log(2.)/log(logstepy)) + 2;

      class_calloc(yy,y_size,sizeof(double),pnl->error_message);

      class_calloc(v_le,y_size,sizeof(double),pnl->error_message);
      class_calloc(v_ri,y_size,sizeof(double),pnl->error_message);

      class_calloc(partial_sum,y_size-1,sizeof(double),pnl->error_message);
      class_calloc(partial_area,y_size-1,sizeof(double),pnl->error_message);

      index_y = 0;

      /* affect the values of y */

      do {

	class_test(index_y >= y_size,
		   pnl->error_message," ");

	yy[index_y] = k*sqrt(2.) - k/sqrt(2.)*pow(logstepy,index_y);

	if (yy[index_y] < 0.) yy[index_y] = 0.; /* correct the last value */

	index_y ++;

      } while (yy[index_y-1] > 0.);

      y_size = index_y; /* to enforce it and avoid problems */

      /*
       *compute first h and v lines
       */

      h_do[0]=0.;
      v_ri[0]=h_do[0];

      for (index_x=1; index_x < x_size; index_x ++) {

	x=xx[index_x];
	y=yy[0];

	if (x <= sqrt(2.)*k_max) {
	  if(pnl->mode==1){
	    class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),&h_do[index_x]),
		       pnl->error_message,
		       pnl->error_message);
	  }
	  if(pnl->mode==2){
	   class_call(trg_A_arg_trg(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,&h_do[index_x]),
		       pnl->error_message,
		       pnl->error_message);
	  }


	}
	else {
	  h_do[index_x]=0.;
	}

      }

      for (index_y=1; index_y < y_size; index_y ++) {

	x=xx[0];
	y=yy[index_y];

	if(pnl->mode==1){
	  class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),&v_ri[index_y]),
		     pnl->error_message,
		     pnl->error_message);
	}
	if(pnl->mode==2){
	  class_call(trg_A_arg_trg(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,&v_ri[index_y]),
		     pnl->error_message,
		     pnl->error_message);
	}

      }

      sum = 0.;
      area = 0.;

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
	  class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),&h_do[il+1]),
		     pnl->error_message,
		     pnl->error_message);
	}
	if(pnl->mode==2){
	  class_call(trg_A_arg_trg(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,&h_do[il+1]),
		     pnl->error_message,
		     pnl->error_message);
	}

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

	    if (x <= sqrt(2)*k_max) {
	      if(pnl->mode==1){
		class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),&h_up[index_x+1]),
			   pnl->error_message,
			   pnl->error_message);
	      }
	      if(pnl->mode==2){
		class_call(trg_A_arg_trg(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,&h_up[index_x+1]),
			   pnl->error_message,
			   pnl->error_message);
	      }
	    }
	    else {
	      h_up[index_x+1]=0.;
	    }

	  }

	  /* the point h_do[index_x+1] is new; calculate */

	  x=xx[index_x+1];
	  y=yy[il+1];

	  if (x <= sqrt(2)*k_max) {
	    if(pnl->mode==1){
	      class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),&h_do[index_x+1]),
			 pnl->error_message,
			 pnl->error_message);
	    }
	    if(pnl->mode==2){
	      class_call(trg_A_arg_trg(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,&h_do[index_x+1]),
			 pnl->error_message,
			 pnl->error_message);
	    }
	  }
	  else {
	    h_do[index_x+1]=0.;
	  }

	  /* the point v_le[index_x+1] may have not been calculated at the previous stage; check and calculate */

	  if (v_le[index_x+1] == 0.) {

	    x=xx[il];
	    y=yy[index_x+1];

	    if(pnl->mode==1){
	      class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),&v_le[index_x+1]),
			 pnl->error_message,
			 pnl->error_message);
	    }
	    if(pnl->mode==2){
	      class_call(trg_A_arg_trg(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,&v_le[index_x+1]),
			 pnl->error_message,
			 pnl->error_message);
	    }

	  }

	  /* the point v_ri[index_x+1] is new; calculate */

	  x=xx[il+1];
	  y=yy[index_x+1];

	  if(pnl->mode==1){
	    class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),&v_ri[index_x+1]),
		       pnl->error_message,
		       pnl->error_message);
	  }
	  if(pnl->mode==2){
	    class_call(trg_A_arg_trg(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,&v_ri[index_x+1]),
		       pnl->error_message,
		       pnl->error_message);
	  }
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

	  if (h_up[index_x+1] == 0.) {

	    x=xx[index_x+1];
	    y=yy[il];

	    if (x <= sqrt(2)*k_max) {
	      if(pnl->mode==1){
		class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),&h_up[index_x+1]),
		    pnl->error_message,
		    pnl->error_message);
	      }
	      if(pnl->mode==2){
		class_call(trg_A_arg_trg(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,&h_up[index_x+1]),
		    pnl->error_message,
		    pnl->error_message);
	      }
	    }
	    else {
	      h_up[index_x+1]=0.;
	    }

	  }

	  /* the point h_do[index_x+1] is new; calculate */

	  x=xx[index_x+1];
	  y=yy[il+1];

	  if (x <= sqrt(2)*k_max) {
	    if(pnl->mode==1){
	      class_call(trg_A_arg_one_loop(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),&h_do[index_x+1]),
		  pnl->error_message,
		  pnl->error_message);
	    }
	    if(pnl->mode==2){
	      class_call(trg_A_arg_trg(pnl,name,k,(x+y)/sqrt(2.),(x-y)/sqrt(2.),index_eta,&h_do[index_x+1]),
		  pnl->error_message,
		  pnl->error_message);
	    }
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

/**
 * Logarithmic step in k for the low k values
 *
 * It simply affects the correct logstep corresponding to a certain
 * value of k.  A hyperbolic tangent has been selected to get a smooth
 * tranisition between the two ends of the spectrum, in order to avoid
 * numerical troubles in case of a brutal change in step.
 *
 * The values are adjusted to fit the default values of steps in eta
 * and of double escape parameters. Change at your own risk !
 *
 * @param k 		Input:  current value of the wavenumber
 * @param logstep	Output: logstep associated with the desired wavenumber
 */
int trg_logstep1_k (
		   struct precision * ppr,
		   double k,
		   double * logstep) {

 *logstep = ppr->logstepk1 - ppr->logstepk2*tanh((k-ppr->logstepk4)*ppr->logstepk3);
  return _SUCCESS_;
}

/**
 * Logarithmic step in k for the high k values
 *
 * It simply affects the correct logstep corresponding to a certain
 * value of k.  A hyperbolic tangent has been selected to get a smooth
 * tranisition between the two ends of the spectrum, in order to avoid
 * numerical troubles in case of a brutal change in step.
 *
 * The values are adjusted to fit the default values of steps in eta
 * and of double escape parameters. Change at your own risk !
 *
 * @param k 		Input:  current value of the wavenumber
 * @param logstep	Output: logstep associated with the desired wavenumber
 */

int trg_logstep2_k (
		    struct precision * ppr,
		    double k,
		    double * logstep) {
  *logstep = ppr->logstepk5 - ppr->logstepk6*tanh((k-ppr->logstepk8)*ppr->logstepk7);
  return _SUCCESS_;
}

/**
 * Computes the non-linear spectra for matter (cold dark matter and
 * baryons, so far) and fills in the spectra_nl structure with pk_nl,
 * p_12_nl, p_22_nl.
 *
 *
 * Since this sub-routine uses a different k table, one should also
 * uses pnl->k to output the results.
 *
 * The program can be told to compute the non linear spectra with
 * (resp.) the Time Renormalization Group method, the one-loop method,
 * or the simple linear method (for testing purposes), by letting the
 * non-linearity mode option in a .ini file be (resp.) 2, 1 or 0.
 *
 * The starting redshift of the non-linear computation can be selected
 * via the .ini file as well (CHECK TROUBLES IF THIS VALUE IS CHANGED
 * AND NOT Z_MAX_PK).  It is set by default to z=35.
 *
 * It is possible (so far) to compute exactly the cdm + baryon
 * spectra, or the approximated cdm + baryon via the option b+c
 * spectrum (YET TO BE IMPLEMENTED)
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure (provides H, Omega_m at redshift of interest)
 * @param pth Input: pointer to thermo structure (provides the baryon sound speed)
 * @param ppm Input: pointer to primordial structure
 * @param psp Input: pointer to spectra structure
 * @param pnl Output: pointer to initialized spectra_nl structure
 * @return the error status
 */

int trg_init (
	      struct precision * ppr,
	      struct background * pba,
	      struct thermo *pth,
	      struct primordial * ppm,
	      struct spectra * psp,
	      struct spectra_nl * pnl
	      ) {

  /** Summary: 	the code takes as an input the matter power spectrum at desired starting redshift,
   * 		and computes its evolution afterwards up to z=0.
   *
   * 		For each step in time, it first computes the non-linear quantities called AA in this
   * 		module, a bunch of 16 non equivalent terms. This is the time-expensive part of the code.
   * 		Then it computes the spectra at the following step.
   *
   * 		To avoid numerical divergences, a method called 'double escape' is applied. It consists
   * 		in droping every half step some points (defined with the double_escape parameter). This
   * 		way, every starting numerical instabily may be taken care of). As it is, the values of
   * 		this parameter, along with eta_step and logstep_k are fixed and balanced to produce an
   * 		output at z=0 with no divergences. Changes are at your own risks.*/

  /** Variables for time control of the computation */

  time_t time_1,time_2;

  /** Wave-number (k) quantities */

  int index_k;

  double * temp_k;
  double logstepk;
  double fourpi_over_k;

  /** Time (t or eta) quantities */

  int index_eta;

  int index;
  int index_plus;

  double a_ini;
  double eta_max;
  double exp_eta;

  /** Variables of calculation */

  double *p_11_linear;
  double *p_12_linear;
  double *p_22_linear;

  double time_step;
  int index_int;

  /** Background quantities */

  double * pvecback_nl;
  double * Omega_m, * H, *H_prime;
  double * growth_factor;
  double * corrected_growth_factor;

  int a,b;
  double gprime;
  double aprime;

  /**
   * Definition of the matrix Omega that mixes terms together,
   * two are k indepedant and two dependant
   */

  double Omega_11,Omega_12;
  double *Omega_21,*Omega_22;

  /**
   * Definition of the variables 1 ijk,lmn and 2 ijk,lmn that
   * replace the Bispectra variables, and of the A's.
   */

  double *a0,*a11,*a12,*a13,*a21,*a22,*a23,*a3;
  double *b0,*b11,*b12,*b21,*b22,*b3;

  int index_name,name_size;
  double ** AA;

  /** Additionnal temporary variables */

  double temp;
  double * junk;
  double pk,pk_ini;

  double dtau;
  double dz_p;
  double delta_minus,delta_plus,d_delta_m_over_dz;

  /** Meaningless initialisations to remove warnings (some variables are only defined in conditional loops) **/

  logstepk=0;
  class_calloc(junk,2,sizeof(double),pnl->error_message); // This one would be to use several different IC in spectra_pk_at_k_and_z


  /**
   * This code can be optionally compiled with the openmp option for parallel computation.
   *
   * Inside parallel regions, the use of the command "return" is forbidden.
   *
   * For error management, instead of "return _FAILURE_", we will set the variable below
   * to "abort = _TRUE_". This will lead to a "return _FAILURE_" jus after leaving the
   * parallel region.
   */
  int abort;

#ifdef _OPENMP
  /** instrumentation times */
  double tstart, tstop;
#endif

  if (pnl->spectra_nl_verbose > 0){
    if(pnl->mode==0)
      printf(" -> testing trg.c module with linear computation\n");
    if(pnl->mode==1)
      printf(" -> using one-loop method\n");
    if(pnl->mode==2)
      printf(" -> using TRG method\n");
  }

  class_test(pba->has_ncdm == _TRUE_,
	     pnl->error_message,
	     "module trg.c not yet generalized to cases with non-cold relics");

  class_test(pba->has_fld == _TRUE_,
	     pnl->error_message,
	     "module trg.c not yet generalized to cases with dark energy");

  class_calloc(pvecback_nl,pba->bg_size,sizeof(double),pnl->error_message);

  /** define initial eta, redshift and size factor */

  a_ini = pba->a_today/(pnl->z_ini + 1.);

  /** define eta_max, where eta=log(a/a_ini) */

  eta_max = log(pba->a_today/a_ini);

  /** define size and step for integration in eta, at any time now, a = a_ini * exp(eta) and z=exp(eta)-1 */

  pnl->eta_step = (eta_max)/(pnl->eta_size-1);

  /** find total number of k values in the module */

  index_k=0;
  class_calloc(temp_k,20000,sizeof(double),pnl->error_message);
  temp_k[0]=pnl->k_min;

  while(temp_k[index_k]<1){
    class_test(index_k>=20000,pnl->error_message,"Change initial size of temp_k\n");
    class_call(trg_logstep1_k(ppr,temp_k[index_k],&logstepk),
	       pnl->error_message,
	       pnl->error_message);
    index_k++;
    temp_k[index_k]=temp_k[index_k-1]*logstepk;
  }

  temp_k[index_k]=temp_k[index_k-1]*logstepk;

  while(temp_k[index_k]<pnl->k_max){
    class_test(index_k>=20000,pnl->error_message,"Change initial size of temp_k\n");
    class_call(trg_logstep2_k(ppr,temp_k[index_k],&logstepk),
	       pnl->error_message,
	       pnl->error_message);
    index_k++;
    temp_k[index_k]=temp_k[index_k-1]*logstepk;
  }

  pnl->k_size=index_k;

  /** Define the values of k, and find the index separating the linear approximation from the non-linear one */

  class_calloc(pnl->k,pnl->k_size,sizeof(double),pnl->error_message);

  temp=0;

  for(index_k=0; index_k<pnl->k_size; index_k++){
    pnl->k[index_k]=temp_k[index_k];
    if( (pnl->k[index_k] > pnl->k_L*pba->h) && temp==0){
      pnl->index_k_L=index_k;
      temp++;
    }
  }

  free(temp_k);

  /** Define time-related values. First check wheter eta_size is odd:
      otherwise the pattern of integration with the predictor
      corrector algorithm is not well behaved */

  class_test((pnl->eta_size==pnl->eta_size % 2),
	     pnl->error_message,
	     "Please use an odd number for eta_size in trg.c module\n");

  class_calloc(pnl->eta,pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(pnl->z,  pnl->eta_size,sizeof(double),pnl->error_message);

  for(index_eta=0;index_eta<pnl->eta_size; index_eta++){
    pnl->eta[index_eta]=index_eta*pnl->eta_step;
    pnl->z[index_eta]= exp(-pnl->eta[index_eta])*(pba->a_today/a_ini)-1.;
    if(index_eta==pnl->eta_size-1){
      pnl->z[index_eta]=0;
    }
  }

  class_test((pnl->z[0]>psp->z_max_pk),
	     pnl->error_message,
	     "\nStarting redshift for non-linear is %2.1f > maximum computed redshift %2.1f\nYou probably want to increase z_max_pk\n",
	     pnl->z[0],psp->z_max_pk);

  if (pnl->spectra_nl_verbose > 1)
    printf(" -> starting calculation at redshift z = %2.2f\n",pnl->z[0]);

  /** Definition of background values */

  class_calloc(Omega_m,pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(H,      pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(H_prime,pnl->eta_size,sizeof(double),pnl->error_message);

  for (index_eta=0; index_eta < pnl->eta_size; index_eta++){
    class_call(background_functions(
				    pba,
				    a_ini*exp(pnl->eta[index_eta]),
				    pba->long_info,
				    pvecback_nl
				    ),
	       pba->error_message,
	       pnl->error_message);

    Omega_m[index_eta] = pvecback_nl[pba->index_bg_Omega_m];

    H[index_eta] = pvecback_nl[pba->index_bg_H] * a_ini * exp(pnl->eta[index_eta]);
    H_prime[index_eta] =H[index_eta]*(1 + pvecback_nl[pba->index_bg_H_prime] / a_ini * exp(-pnl->eta[index_eta])/pvecback_nl[pba->index_bg_H]/pvecback_nl[pba->index_bg_H]);

  }

  free(pvecback_nl);

  /** Definition of the matrix elements Omega_11,Omega_12, Omega_22,Omega_12 for each k and eta */

  Omega_11 = 1.;
  Omega_12 = -1.;

  class_calloc(Omega_21, pnl->eta_size * pnl->k_size,sizeof(double),pnl->error_message);
  class_calloc(Omega_22, pnl->eta_size * pnl->k_size,sizeof(double),pnl->error_message);

  for(index_eta=0; index_eta < pnl->eta_size; index_eta++) {
    for(index_k=0; index_k < pnl->k_size; index_k++) {

      index=index_k+pnl->k_size*index_eta;

      Omega_21[index] = -3./2 * Omega_m[index_eta];
      Omega_22[index] = 2 + H_prime[index_eta]/H[index_eta];
    }
  }

  /** Implementing the E1L method (see paper). We recover the computed
      growth_factor, normalized to unity at starting redshift, and for
      a k_growth_factor chosen in the .ini. The dependence in
      k_growth_factor has been tested to be negligible. */

  class_calloc(growth_factor,pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(corrected_growth_factor,pnl->eta_size,sizeof(double),pnl->error_message);

  class_call(spectra_pk_at_k_and_z(pba,ppm,psp,pnl->k_growth_factor,pnl->z[0],&pk_ini,junk),
	     psp->error_message,
	     pnl->error_message);

  growth_factor[0] = 1.;
  corrected_growth_factor[0] = 1.;

  for(index_eta=1; index_eta < pnl->eta_size; index_eta++) {

    class_call(spectra_pk_at_k_and_z(pba,ppm,psp,pnl->k_growth_factor,pnl->z[index_eta],&pk,junk),
               psp->error_message,
               pnl->error_message);

    /* Putting this following value to 1 for each eta would recover the U1L method */
      growth_factor[index_eta]=sqrt(pk/pk_ini)*(1+pnl->z[index_eta])/(1+pnl->z[0]);

      gprime=(growth_factor[index_eta]-growth_factor[index_eta-1]);
      aprime=-1./pow(1+pnl->z[index_eta],2)*(pnl->z[index_eta]-pnl->z[index_eta-1]);

      corrected_growth_factor[index_eta]=growth_factor[index_eta]+1./(1+pnl->z[index_eta])*gprime/aprime;
  }


  /**
   * Definition of P_11, P_12 and P_22, the two points correlators
   * of the density/density (resp. density/velocity and velocity/velocity),
   * initialized at eta=0 (by default, z=35).
   *
   * The p_ij_nl version are keeping track of the non linear spectra. The p_ij
   * are only here for the one-loop computation. The p_ij_linear, not stored
   * in the nl_spectra module, keep trace of the linear evolution of this code.
   *
   * With the default parameter choice, the relative difference between the linear
   * evolution in this part and the more carefull time integration from CLASS
   * is no more than 1.6% for the peak height in the spectrum, at z=0.
   */

  class_calloc(pnl->p_11_nl,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(pnl->p_12_nl,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(pnl->p_22_nl,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);

  class_calloc(pnl->p_11,pnl->k_size,sizeof(double),pnl->error_message);
  class_calloc(pnl->p_12,pnl->k_size,sizeof(double),pnl->error_message);
  class_calloc(pnl->p_22,pnl->k_size,sizeof(double),pnl->error_message);

  class_calloc(p_11_linear,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(p_12_linear,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(p_22_linear,pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);

  /* There is the possibility to add a cutoff to take into account
   * exponential suppression of the power spectrum at high k. By
   * default, however, set on 1. (no cutoff).
   */

  /* To determine precisely the velocity at initial time, the code
     evaluates the time derivative of the density (relation true for
     linear perturbation theory) using the following local variables.
   */

  dtau=0.0001; /*conformal age of the universe is 14000 */
  dz_p=-dtau*pba->a_today*H[0]/a_ini;

  for(index_k=0; index_k<pnl->k_size; index_k++){

    class_call(spectra_pk_at_k_and_z(pba,ppm,psp,pnl->k[index_k],pnl->z[0],&pnl->p_11_nl[index_k],junk),
               psp->error_message,
               pnl->error_message);

    class_call(spectra_pk_at_k_and_z(pba,ppm,psp,pnl->k[index_k],pnl->z[0]+dz_p,&delta_plus,junk),
               psp->error_message,
               pnl->error_message);

    class_call(spectra_pk_at_k_and_z(pba,ppm,psp,pnl->k[index_k],pnl->z[0]-dz_p,&delta_minus,junk),
               psp->error_message,
               pnl->error_message);


    pnl->p_11[index_k]    = pnl->p_11_nl[index_k];
    p_11_linear[index_k]  = pnl->p_11_nl[index_k];

    d_delta_m_over_dz=(sqrt(delta_plus)-sqrt(delta_minus))/(2*dz_p);

    pnl->p_12_nl[index_k] = -sqrt(pnl->p_11_nl[index_k])*d_delta_m_over_dz*pba->a_today/a_ini;
    pnl->p_12[index_k]    = pnl->p_12_nl[index_k];
    p_12_linear[index_k]  = pnl->p_12_nl[index_k];

    pnl->p_22_nl[index_k] = pow(d_delta_m_over_dz*pba->a_today/a_ini,2);
    pnl->p_22[index_k]    = pnl->p_22_nl[index_k];
    p_22_linear[index_k]  = pnl->p_22_nl[index_k];
  }


  free(junk);
  free(Omega_m);
  free(H_prime);

  /** Initialisation and definition of second derivatives */

  class_calloc(pnl->ddp_11_nl,  pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(pnl->ddp_12_nl,  pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(pnl->ddp_22_nl,  pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);

  class_calloc(pnl->ddp_11,     pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(pnl->ddp_12,     pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);
  class_calloc(pnl->ddp_22,     pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);

  class_call(array_logspline_table_one_column(pnl->k,pnl->k_size,pnl->k_size,
					      pnl->p_11_nl,
					      pnl->eta_size,0,
					      pnl->ddp_11_nl,
					      _SPLINE_NATURAL_,pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);

  class_call(array_logspline_table_one_column(pnl->k,pnl->k_size,pnl->k_size,
					      pnl->p_12_nl,
					      pnl->eta_size,0,
					      pnl->ddp_12_nl,
					      _SPLINE_NATURAL_,pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);

  class_call(array_logspline_table_one_column(pnl->k,pnl->k_size,pnl->k_size,
					      pnl->p_22_nl,
					      pnl->eta_size,0,
					      pnl->ddp_22_nl,
					      _SPLINE_NATURAL_,pnl->error_message),
	     pnl->error_message,
	     pnl->error_message);

  for (index_k = 0; index_k < pnl->k_size ; index_k++) {
    pnl->ddp_11[index_k] = pnl->ddp_11_nl[index_k];
    pnl->ddp_12[index_k] = pnl->ddp_12_nl[index_k];
    pnl->ddp_22[index_k] = pnl->ddp_22_nl[index_k];
  }
   /* Definition of 1_0, 1_11,(here a0, a11,...) etc, and 2_0, 2_11,
    * (here b0,b11,...) etc.. and initialization (directly with calloc
    * for assuming no initial non gaussianity in the spectrum) The
    * convention for 1_0, 1_11, 1_22 is : I_121_222, I_121_122,
    * I_121_121, in more details:
    * 0 : there is no 1 in the second row of indices,
    * 11: there is one 1 in the second row of indices, and it is at the first place,
    * 22: there are two 1's and the only 2 is at the second place.

    * The convention for 2_0, etc.. are reversed, i.e. 2_11 means
    * I_222_211. Though not very clear at first, I found no clearer
    * and faster way to designate all this cumbersome indices, and in
    * the end, it is not so bad.
    */

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
   * the 1's and 2's, and initialization with pk_nl, p_12,
   * p_22. Convention is taken for A_121,any = AA[first indices] and A_222,any = AA[last indices]
   */

  name_size = _B3_+1;
  class_calloc(AA,name_size,sizeof(double*),pnl->error_message);
  for (index_name=0; index_name<name_size; index_name++)
    class_calloc(AA[index_name],pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);

  /** First step in time, computing the non-linear quantities AA */

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

	class_call_parallel(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,index_name,0,AA[index_name]),
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

  free(pnl->p_11);
  free(pnl->p_12);
  free(pnl->p_22);

  /** Correct setup of initial conditions
      There are two choices, the default is setting-up with 1-loop initial bispectrum,
      with the parameter trg_inicond = pt,
      or with gaussian initial condtion, with trg_inicond = lin */

  if(pnl->ic==nl_pt){

    for (index_k=0; index_k<pnl->k_size; index_k++){
      a0[index_k] = 2.*AA[_A0_][index_k];
      a11[index_k]= 2.*AA[_A11_][index_k];
      a12[index_k]= 2.*AA[_A12_][index_k];
      a13[index_k]= 2.*AA[_A13_][index_k];
      a21[index_k]= 2.*AA[_A21_][index_k];
      a22[index_k]= 2.*AA[_A22_][index_k];
      a23[index_k]= 2.*AA[_A23_][index_k];
      a3[index_k] = 2.*AA[_A3_][index_k];

      b0[index_k] = 2.*AA[_B0_][index_k];
      b11[index_k]= 2.*AA[_B11_][index_k];
      b12[index_k]= 2.*AA[_B12_][index_k];
      b21[index_k]= 2.*AA[_B21_][index_k];
      b22[index_k]= 2.*AA[_B22_][index_k];
      b3[index_k] = 2.*AA[_B3_][index_k];
    }
  }

  /** Now we calculate the time evolution with a predictor corrector
      algorithm. */

  /* At each step, first the new power spectra are computed on a half
   * time step (and pnl->double_escape points are dropped) then the AA
   * functions are updated (for 1-loop method, they are just copied
   * from last step) (and pnl->double_escape points are dropped).
   * Then, the new power spectra are evaluated in the whole time step
   * using quantities computed on the half time step, and the new AA
   * functions are updated. The pnl->double_escape procedure also
   * takes place during this time.
   */

  if (pnl->spectra_nl_verbose > 1){
    printf(" -> progression:\n");}

  time_1=time(NULL);
  time_2=time(NULL);

  for (index_eta=2; index_eta<pnl->eta_size; index_eta+=2){
    exp_eta=exp(pnl->eta[index_eta-2]);
    time_step=pnl->eta[index_eta-1]-pnl->eta[index_eta-2];

    if(pnl->k_size-pnl->double_escape*2*(index_eta)/1 < 2){
      printf("  --> Wrong choice of double escape and eta_size parameters\n  --> Stopped at %d, z=%2.2e Try Again\n",index_eta, pnl->z[index_eta]);
      return _FAILURE_;
    }

    for (index_k=0; index_k<pnl->k_size-pnl->double_escape*2*(index_eta-1)/1; index_k++){

      /* Some useful intermediate variables */

      fourpi_over_k=4*_PI_/(pnl->k[index_k]);
      index = index_k+pnl->k_size*(index_eta-2);
      index_plus = index_k+pnl->k_size*(index_eta-1);

     /* For each k, compute the power spectra and bispectra (through ai) at the new time */

      pnl->p_11_nl[index_plus]= (time_step) *(
					 -2*Omega_11*pnl->p_11_nl[index]
					 -2*Omega_12*pnl->p_12_nl[index]
					 +exp_eta*4*fourpi_over_k*a22[index] )
	+ pnl->p_11_nl[index];

      pnl->p_22_nl[index_plus] = (time_step) *(
					    -2*Omega_22[index]*pnl->p_22_nl[index]
					    -2*Omega_21[index]*pnl->p_12_nl[index]
					    +exp_eta*2*fourpi_over_k*b3[index] )
	+ pnl->p_22_nl[index];

      pnl->p_12_nl[index_plus] = (time_step) *(
					    -pnl->p_12_nl[index]*(Omega_11+Omega_22[index])
					    -Omega_12*pnl->p_22_nl[index]
					    -Omega_21[index]*pnl->p_11_nl[index]
					    +exp_eta*fourpi_over_k*(2*a13[index]+b21[index]))
	+ pnl->p_12_nl[index];


      a0[index_plus]         = (time_step) *(
					  -Omega_21[index]*(a11[index]+a12[index]+a13[index])
					  -3*Omega_22[index]*a0[index]
					  +2*exp_eta*AA[_A0_][index])
	+ a0[index];

      a11[index_plus]        = (time_step) *(
					  -a11[index]*(2*Omega_22[index]+Omega_11)
					  -Omega_12*a0[index]
					  -Omega_21[index]*(a22[index]+a23[index])
					  +2*exp_eta*AA[_A11_][index])
	+ a11[index];

      a12[index_plus]        = (time_step) *(
					  -a12[index]*(2*Omega_22[index]+Omega_11)
					  -Omega_21[index]*(a23[index]+a21[index])
					  -Omega_12*a0[index]
					  +2*exp_eta*AA[_A12_][index])
	+ a12[index];

      a13[index_plus]        = (time_step) *(
					  -a13[index]*(2*Omega_22[index]+Omega_11)
					  -Omega_12*a0[index]
					  -Omega_21[index]*(a22[index]+a21[index])
					  +2*exp_eta*AA[_A13_][index])
	+ a13[index];

      a21[index_plus]        = (time_step) *(
					  -a21[index]*(2*Omega_11+Omega_22[index])
					  -Omega_12*(a12[index]+a13[index])
					  -Omega_21[index]*a3[index]
					  +2*exp_eta*AA[_A21_][index])
	+ a21[index];

      a22[index_plus]        = (time_step) *(
					  -a22[index]*(2*Omega_11+Omega_22[index])
					  -Omega_12*(a13[index]+a11[index])
					  -Omega_21[index]*a3[index]
					  +2*exp_eta*AA[_A22_][index])
	+ a22[index];

      a23[index_plus]        = (time_step) *(
					  -a23[index]*(2*Omega_11+Omega_22[index])
					  -Omega_12*(a12[index]+a11[index])
					  -Omega_21[index]*a3[index]
					  +2*exp_eta*AA[_A23_][index])
	+ a23[index];

      a3[index_plus]         = (time_step) *(
					  -a3[index]*3*Omega_11
					  -Omega_12*(a22[index]+a21[index]+a23[index])
					  +2*exp_eta*AA[_A3_][index])
	+ a3[index];

      b0[index_plus]         = (time_step) *(
					  -3*b0[index]*Omega_11
					  -Omega_12*(b11[index]+2*b12[index])
					  +2*exp_eta*AA[_B0_][index])
	+ b0[index];

      b11[index_plus]        = (time_step) *(
					  -b11[index]*(2*Omega_11+Omega_22[index])
					  -2*Omega_12*b22[index]
					  -Omega_21[index]*b0[index]
					  +2*exp_eta*AA[_B11_][index])
	+ b11[index];

      b12[index_plus]        = (time_step) *(
					  -b12[index]*(2*Omega_11+Omega_22[index])
					  -Omega_12*(b22[index]+b21[index])
					  -Omega_21[index]*b0[index]
					  +2*exp_eta*AA[_B12_][index])
	+ b12[index];

      b21[index_plus]        = (time_step) *(
					  -b21[index]*(2*Omega_22[index]+Omega_11)
					  -2*Omega_21[index]*b12[index]
					  -Omega_12*b3[index]
					  +2*exp_eta*AA[_B21_][index])
	+ b21[index];

      b22[index_plus]        = (time_step) *(
					  -b22[index]*(2*Omega_22[index]+Omega_11)
					  -Omega_21[index]*(b12[index]+b11[index])
					  -Omega_12*b3[index]
					  +2*exp_eta*AA[_B22_][index])
	+ b22[index];

      b3[index_plus]         = (time_step) *(
					  -3*Omega_22[index]*b3[index]
					  -Omega_21[index]*(b21[index]+2*b22[index])
					  +2*exp_eta*AA[_B3_][index])
	+ b3[index];

      /* The linear quantities are followed through this simplified integrator, for reference */

      p_11_linear[index_plus]= (time_step) *(
					  -2*Omega_11*p_11_linear[index]
					  -2*Omega_12*p_12_linear[index]
					  )
	+ p_11_linear[index];

      p_22_linear[index_plus] = (time_step) *(
					   -2*Omega_22[index]*p_22_linear[index]
					   -2*Omega_21[index]*p_12_linear[index]
					   )
	+ p_22_linear[index];

      p_12_linear[index_plus] = (time_step) *(
					   -p_12_linear[index]*(Omega_11+Omega_22[index])
					   -Omega_12*p_22_linear[index]
					   -Omega_21[index]*p_11_linear[index]
					   )
	+ p_12_linear[index];
    }

    /** Update of second derivatives for interpolation, only in TRG mode */

    if(pnl->mode==2){
      class_call(array_logspline_table_one_column(pnl->k,pnl->k_size,pnl->k_size-pnl->double_escape*2*(index_eta)/1,
						  pnl->p_11_nl,
						  pnl->eta_size,index_eta-1,
						  pnl->ddp_11_nl,
						  _SPLINE_NATURAL_,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);

      class_call(array_logspline_table_one_column(pnl->k,pnl->k_size,pnl->k_size-pnl->double_escape*2*(index_eta)/1,
						  pnl->p_12_nl,
						  pnl->eta_size,index_eta-1,
						  pnl->ddp_12_nl,
						  _SPLINE_NATURAL_,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);

      class_call(array_logspline_table_one_column(pnl->k,pnl->k_size,pnl->k_size-pnl->double_escape*2*(index_eta)/1,
						  pnl->p_22_nl,
						  pnl->eta_size,index_eta-1,
						  pnl->ddp_22_nl,
						  _SPLINE_NATURAL_,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);
    }

    /** Update of AA's at the new time
     */

    /* For Evolved-1-Loop, copy the previous values, taking into
     * account the correct powers of growth factor g for density, and
     * f=g+ag'/a' for velocity perturbation) */

    if(pnl->mode==1){
      for (index_name=0; index_name<name_size; index_name++){
	if(index_name==0)
	  a=0;
	else if(index_name==1||index_name==2||index_name==3||index_name==11||index_name==12)
	  a=1;
	else if(index_name==4||index_name==5||index_name==6||index_name==9||index_name==10)
	  a=2;
	else if(index_name==7||index_name==8)
	  a=3;
	else
	  a=4;
	b=4-a;
	for(index_k=0; index_k<pnl->k_size-pnl->double_escape*2*(index_eta)/1; index_k++){
	  AA[index_name][index_k+pnl->k_size*(index_eta-1)]=pow(growth_factor[index_eta-1],a)*pow(corrected_growth_factor[index_eta-1],b)*AA[index_name][index_k];}
      }
    }

    /* For TRG, simply recomputes the AA values entirely */

    else if(pnl->mode == 2){

      /* initialize error management flag */
      abort = _FALSE_;

      /*** beginning of parallel region ***/

#pragma omp parallel							\
      shared(name_size,abort,pba,ppm,psp,pnl,index_eta,AA)		\
      private(tstart,index_name,tstop)

      {

#ifdef _OPENMP
	tstart = omp_get_wtime();
#endif

#pragma omp for schedule (dynamic)
	for (index_name=0; index_name<name_size; index_name++) {

#pragma omp flush(abort)

	  class_call_parallel(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,index_name,index_eta-1,AA[index_name]),
	      pnl->error_message,
	      pnl->error_message);
	}


#ifdef _OPENMP
	tstop = omp_get_wtime();
	if ((pnl->spectra_nl_verbose > 2) && (pnl->mode > 1))
	  printf("In %s: time spent in parallel region (loop over names) = %e s for thread %d\n",
	      __func__,tstop-tstart,omp_get_thread_num());
#endif

      } /* end of parallel region */

      if (abort == _TRUE_) return _FAILURE_;

    }

    //CORRECTOR
    exp_eta=exp(pnl->eta[index_eta-1]);
    time_step=pnl->eta[index_eta]-pnl->eta[index_eta-2];

    for (index_k=0; index_k<pnl->k_size-pnl->double_escape*2*(index_eta)/1; index_k++){

      /* Some useful intermediate variables */

      fourpi_over_k=4*_PI_/(pnl->k[index_k]);
      index_int= index_k+pnl->k_size*(index_eta-2);
      index = index_k+pnl->k_size*(index_eta-1);
      index_plus = index_k+pnl->k_size*index_eta;

     /* For each k, compute the power spectra and bispectra (through ai) at the new time */

      pnl->p_11_nl[index_plus]= (time_step) *(
					 -2*Omega_11*pnl->p_11_nl[index]
					 -2*Omega_12*pnl->p_12_nl[index]
					 +exp_eta*4*fourpi_over_k*a22[index] )
	+ pnl->p_11_nl[index_int];

      pnl->p_22_nl[index_plus] = (time_step) *(
					    -2*Omega_22[index]*pnl->p_22_nl[index]
					    -2*Omega_21[index]*pnl->p_12_nl[index]
					    +exp_eta*2*fourpi_over_k*b3[index] )
	+ pnl->p_22_nl[index_int];

      pnl->p_12_nl[index_plus] = (time_step) *(
					    -pnl->p_12_nl[index]*(Omega_11+Omega_22[index])
					    -Omega_12*pnl->p_22_nl[index]
					    -Omega_21[index]*pnl->p_11_nl[index]
					    +exp_eta*fourpi_over_k*(2*a13[index]+b21[index]))
	+ pnl->p_12_nl[index_int];


      a0[index_plus]         = (time_step) *(
					  -Omega_21[index]*(a11[index]+a12[index]+a13[index])
					  -3*Omega_22[index]*a0[index]
					  +2*exp_eta*AA[_A0_][index])
	+ a0[index_int];

      a11[index_plus]        = (time_step) *(
					  -a11[index]*(2*Omega_22[index]+Omega_11)
					  -Omega_12*a0[index]
					  -Omega_21[index]*(a22[index]+a23[index])
					  +2*exp_eta*AA[_A11_][index])
	+ a11[index_int];

      a12[index_plus]        = (time_step) *(
					  -a12[index]*(2*Omega_22[index]+Omega_11)
					  -Omega_21[index]*(a23[index]+a21[index])
					  -Omega_12*a0[index]
					  +2*exp_eta*AA[_A12_][index])
	+ a12[index_int];

      a13[index_plus]        = (time_step) *(
					  -a13[index]*(2*Omega_22[index]+Omega_11)
					  -Omega_12*a0[index]
					  -Omega_21[index]*(a22[index]+a21[index])
					  +2*exp_eta*AA[_A13_][index])
	+ a13[index_int];

      a21[index_plus]        = (time_step) *(
					  -a21[index]*(2*Omega_11+Omega_22[index])
					  -Omega_12*(a12[index]+a13[index])
					  -Omega_21[index]*a3[index]
					  +2*exp_eta*AA[_A21_][index])
	+ a21[index_int];

      a22[index_plus]        = (time_step) *(
					  -a22[index]*(2*Omega_11+Omega_22[index])
					  -Omega_12*(a13[index]+a11[index])
					  -Omega_21[index]*a3[index]
					  +2*exp_eta*AA[_A22_][index])
	+ a22[index_int];

      a23[index_plus]        = (time_step) *(
					  -a23[index]*(2*Omega_11+Omega_22[index])
					  -Omega_12*(a12[index]+a11[index])
					  -Omega_21[index]*a3[index]
					  +2*exp_eta*AA[_A23_][index])
	+ a23[index_int];

      a3[index_plus]         = (time_step) *(
					  -a3[index]*3*Omega_11
					  -Omega_12*(a22[index]+a21[index]+a23[index])
					  +2*exp_eta*AA[_A3_][index])
	+ a3[index_int];

      b0[index_plus]         = (time_step) *(
					  -3*b0[index]*Omega_11
					  -Omega_12*(b11[index]+2*b12[index])
					  +2*exp_eta*AA[_B0_][index])
	+ b0[index_int];

      b11[index_plus]        = (time_step) *(
					  -b11[index]*(2*Omega_11+Omega_22[index])
					  -2*Omega_12*b22[index]
					  -Omega_21[index]*b0[index]
					  +2*exp_eta*AA[_B11_][index])
	+ b11[index_int];

      b12[index_plus]        = (time_step) *(
					  -b12[index]*(2*Omega_11+Omega_22[index])
					  -Omega_12*(b22[index]+b21[index])
					  -Omega_21[index]*b0[index]
					  +2*exp_eta*AA[_B12_][index])
	+ b12[index_int];

      b21[index_plus]        = (time_step) *(
					  -b21[index]*(2*Omega_22[index]+Omega_11)
					  -2*Omega_21[index]*b12[index]
					  -Omega_12*b3[index]
					  +2*exp_eta*AA[_B21_][index])
	+ b21[index_int];

      b22[index_plus]        = (time_step) *(
					  -b22[index]*(2*Omega_22[index]+Omega_11)
					  -Omega_21[index]*(b12[index]+b11[index])
					  -Omega_12*b3[index]
					  +2*exp_eta*AA[_B22_][index])
	+ b22[index_int];

      b3[index_plus]         = (time_step) *(
					  -3*Omega_22[index]*b3[index]
					  -Omega_21[index]*(b21[index]+2*b22[index])
					  +2*exp_eta*AA[_B3_][index])
	+ b3[index_int];

      /* The linear quantities are followed through this simplified integrator, for reference */

      p_11_linear[index_plus]= (time_step) *(
					  -2*Omega_11*p_11_linear[index]
					  -2*Omega_12*p_12_linear[index]
					  )
	+ p_11_linear[index_int];

      p_22_linear[index_plus] = (time_step) *(
					   -2*Omega_22[index]*p_22_linear[index]
					   -2*Omega_21[index]*p_12_linear[index]
					   )
	+ p_22_linear[index_int];

      p_12_linear[index_plus] = (time_step) *(
					   -p_12_linear[index]*(Omega_11+Omega_22[index])
					   -Omega_12*p_22_linear[index]
					   -Omega_21[index]*p_11_linear[index]
					   )
	+ p_12_linear[index_int];
    }

    /** Update of second derivatives for interpolation, only in TRG mode */

    if(pnl->mode==2){
      class_call(array_logspline_table_one_column(pnl->k,pnl->k_size,pnl->k_size-pnl->double_escape*2*(index_eta)/1,
						  pnl->p_11_nl,
						  pnl->eta_size,index_eta,
						  pnl->ddp_11_nl,
						  _SPLINE_NATURAL_,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);

      class_call(array_logspline_table_one_column(pnl->k,pnl->k_size,pnl->k_size-pnl->double_escape*2*(index_eta)/1,
						  pnl->p_12_nl,
						  pnl->eta_size,index_eta,
						  pnl->ddp_12_nl,
						  _SPLINE_NATURAL_,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);

      class_call(array_logspline_table_one_column(pnl->k,pnl->k_size,pnl->k_size-pnl->double_escape*2*(index_eta)/1,
						  pnl->p_22_nl,
						  pnl->eta_size,index_eta,
						  pnl->ddp_22_nl,
						  _SPLINE_NATURAL_,pnl->error_message),
		 pnl->error_message,
		 pnl->error_message);
    }

    /** Update of AA's at the new time (for 1-loop, copy the previous
	values, taking into account the growth factor) */

    /* For Evolved-1-Loop, copy the previous values, taking into
     * account the correct powers of growth factor g for density, and
     * f=g+ag'/a' for velocity perturbation) */

    if(pnl->mode==1){
      for (index_name=0; index_name<name_size; index_name++){
	if(index_name==0)
	  a=0;
	else if(index_name==1||index_name==2||index_name==3||index_name==11||index_name==12)
	  a=1;
	else if(index_name==4||index_name==5||index_name==6||index_name==9||index_name==10)
	  a=2;
	else if(index_name==7||index_name==8)
	  a=3;
	else
	  a=4;
	b=4-a;
	for(index_k=0; index_k<pnl->k_size-pnl->double_escape*2*(index_eta)/1; index_k++){
	  AA[index_name][index_k+pnl->k_size*index_eta]=pow(growth_factor[index_eta-1],a)*pow(corrected_growth_factor[index_eta-1],b)*AA[index_name][index_k];}
      }
    }

    else if(pnl->mode == 2){

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

	  class_call_parallel(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,index_name,index_eta,AA[index_name]),
	      pnl->error_message,
	      pnl->error_message);
	}


#ifdef _OPENMP
	tstop = omp_get_wtime();
	if ((pnl->spectra_nl_verbose > 2) && (pnl->mode > 1))
	  printf("In %s: time spent in parallel region (loop over names) = %e s for thread %d\n",
	      __func__,tstop-tstart,omp_get_thread_num());
#endif

      } /* end of parallel region */

      if (abort == _TRUE_) return _FAILURE_;

    }

    if(pnl->spectra_nl_verbose>1)
      printf("    %2.1f%% done\n",100.*index_eta/(pnl->eta_size-1.));

    time_2=time(NULL);
    if(pnl->spectra_nl_verbose>0){
      if(index_eta==9 && pnl->mode > 1){
	printf("     elapsed time after ten loops : %2.f s\n",difftime(time_2, time_1));
	printf("     max estimated remaining : %3.1f min\n",difftime(time_2,time_1)*(pnl->eta_size-2)/60/10);
      }
    }

    class_test(isnan(pnl->p_11_nl[50+pnl->k_size*index_eta])!=0,pnl->error_message,"Code returns nan!\n");

  }

  if(pnl->spectra_nl_verbose>1) printf("Done in %2.f min\n",difftime(time_2,time_1)/60);

  /** End of the computation, beginning of cleaning */

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

  free(p_11_linear);
  free(p_12_linear);
  free(p_22_linear);

  free(Omega_21);
  free(Omega_22);

  free(growth_factor);
  free(corrected_growth_factor);

  /** Pk_nl values are transformed back into real power spectrum
      values, from the phi doublet notation. */

  for(index_eta=0; index_eta < pnl->eta_size; index_eta++) {
    for(index_k=0; index_k < pnl->k_size-pnl->double_escape*index_eta; index_k++) {

      pnl->p_11_nl[index_k+pnl->k_size*index_eta]*=exp(-log( (pnl->z[index_eta]+1.) * a_ini / pba->a_today )*2);
      pnl->p_12_nl[index_k+pnl->k_size*index_eta]*=exp(-log( (pnl->z[index_eta]+1.) * a_ini / pba->a_today )*2)*sqrt(H[index_eta]);
      pnl->p_22_nl[index_k+pnl->k_size*index_eta]*=exp(-log( (pnl->z[index_eta]+1.) * a_ini / pba->a_today )*2)*H[index_eta];

    }
  }

  free(H);

  return _SUCCESS_;

}

/**
 * This routine frees all the memory space allocated by trg_init().
 *
 * To be called at the end of each run, when no further use of the nl_functions are needed.
 *
 * @param pnl Input : pointer to spectra_nl structure (which fields must be freed).
 * @return the error status
 */

int trg_free(
	     struct spectra_nl * pnl
	     ){

  free(pnl->k);
  free(pnl->p_11_nl);
  free(pnl->p_12_nl);
  free(pnl->p_22_nl);
  free(pnl->ddp_11_nl);
  free(pnl->ddp_12_nl);
  free(pnl->ddp_22_nl);
  free(pnl->ddp_11);
  free(pnl->ddp_12);
  free(pnl->ddp_22);
  free(pnl->eta);
  free(pnl->z);

  return _SUCCESS_;
}
