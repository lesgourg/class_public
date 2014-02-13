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
#include "bessel.h"
#include "transfer.h"
#include "primordial.h"
#include "spectra.h"
#include "trg_bispectrum.h"
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



int trg_G_terms(
    struct spectra_nl * pnl,
    enum name_B name,
    double k,
    double q,
    double p,
    int index_eta,
    double * result
    ){

  /** - define local variables */

  /* shorter names for the spectra */

  double p_11k,p_11q,p_11p;
  double p_12k,p_12q,p_12p;
  double p_22k,p_22q,p_22p;

  /*shorter names for the interaction functions*/

  double gamma1_kqp,gamma1_kpq,gamma1_qpk,gamma1_qkp,gamma1_pkq,gamma1_pqk;
  double gamma2_kqp,gamma2_qpk,gamma2_pkq;

  switch(name){
    case _111_:
      class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,k,&p_11k),
	  pnl->error_message,
	  pnl->error_message);
      class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,q,&p_11q),
	  pnl->error_message,
	  pnl->error_message);
      class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,p,&p_11p),
	  pnl->error_message,
	  pnl->error_message);

      class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,k,&p_12k),
	  pnl->error_message,
	  pnl->error_message);
      class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,q,&p_12q),
	  pnl->error_message,
	  pnl->error_message);
      class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,p,&p_12p),
	  pnl->error_message,
	  pnl->error_message);

      trg_gamma_121(k,q,p,&gamma1_kqp);
      trg_gamma_121(k,p,q,&gamma1_kpq);
      trg_gamma_121(q,p,k,&gamma1_qpk);
      trg_gamma_121(q,k,p,&gamma1_qkp);
      trg_gamma_121(p,k,q,&gamma1_pkq);
      trg_gamma_121(p,q,k,&gamma1_pqk);

      *result = gamma1_kqp*p_12q*p_11p + gamma1_kpq*p_11q*p_12p +
	gamma1_qpk*p_12p*p_11k + gamma1_qkp*p_11p*p_12k +
	gamma1_pkq*p_12k*p_11q + gamma1_pqk*p_11k*p_12q;

      return _SUCCESS_;
      break;

    /****************************************/

    case _121_:
      class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,k,&p_11k),
	  pnl->error_message,
	  pnl->error_message);
      class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,p,&p_11p),
	  pnl->error_message,
	  pnl->error_message);

      class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,k,&p_12k),
	  pnl->error_message,
	  pnl->error_message);
      class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,q,&p_12q),
	  pnl->error_message,
	  pnl->error_message);
      class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,p,&p_12p),
	  pnl->error_message,
	  pnl->error_message);

      class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,q,&p_22q),
	  pnl->error_message,
	  pnl->error_message);


      class_call(trg_gamma_222(q,p,k,&gamma2_qpk,pnl->error_message),
	  pnl->error_message,
	  pnl->error_message);

      trg_gamma_121(k,q,p,&gamma1_kqp);
      trg_gamma_121(k,p,q,&gamma1_kpq);
      trg_gamma_121(p,k,q,&gamma1_pkq);
      trg_gamma_121(p,q,k,&gamma1_pqk);

      *result = gamma1_kqp*p_22q*p_11p + gamma1_kpq*p_12q*p_12p +
	gamma2_qpk*p_12p*p_12k +
	gamma1_pkq*p_12k*p_12q + gamma1_pqk*p_11k*p_22q;

      return _SUCCESS_;
      break;

    /****************************************/

    case _211_:
      class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,p,&p_11p),
	  pnl->error_message,
	  pnl->error_message);
      class_call(trg_p_ab_at_any_k(pnl,pnl->p_11,pnl->ddp_11,0,q,&p_11q),
	  pnl->error_message,
	  pnl->error_message);

      class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,k,&p_12k),
	  pnl->error_message,
	  pnl->error_message);
      class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,q,&p_12q),
	  pnl->error_message,
	  pnl->error_message);
      class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,p,&p_12p),
	  pnl->error_message,
	  pnl->error_message);

      class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,k,&p_22k),
	  pnl->error_message,
	  pnl->error_message);


      class_call(trg_gamma_222(k,q,p,&gamma2_kqp,pnl->error_message),
	  pnl->error_message,
	  pnl->error_message);

      trg_gamma_121(k,q,p,&gamma1_qpk);
      trg_gamma_121(k,p,q,&gamma1_qkp);
      trg_gamma_121(p,k,q,&gamma1_pkq);
      trg_gamma_121(p,q,k,&gamma1_pqk);

      *result = gamma2_kqp*p_12q*p_12p +
       	gamma1_qpk*p_12p*p_12k + gamma1_qkp*p_11p*p_22k +
	gamma1_pkq*p_22k*p_11q + gamma1_pqk*p_12k*p_12q;

      return _SUCCESS_;
      break;

    /****************************************/
    case _122_:

      class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,k,&p_12k),
	  pnl->error_message,
	  pnl->error_message);
      class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,q,&p_12q),
	  pnl->error_message,
	  pnl->error_message);
      class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,p,&p_12p),
	  pnl->error_message,
	  pnl->error_message);

      class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,q,&p_22q),
	  pnl->error_message,
	  pnl->error_message);
      class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,p,&p_22p),
	  pnl->error_message,
	  pnl->error_message);


      class_call(trg_gamma_222(q,p,k,&gamma2_qpk,pnl->error_message),
	  pnl->error_message,
	  pnl->error_message);

      class_call(trg_gamma_222(p,k,q,&gamma2_pkq,pnl->error_message),
	  pnl->error_message,
	  pnl->error_message);

      trg_gamma_121(k,q,p,&gamma1_kqp);
      trg_gamma_121(k,p,q,&gamma1_kpq);

      *result = gamma1_kqp*p_22q*p_12p + gamma1_kpq*p_12q*p_22p +
	gamma2_qpk*p_22p*p_12k + gamma2_pkq*p_12k*p_22q;

      return _SUCCESS_;
      break;

    /****************************************/

    case _221_:

      class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,k,&p_12k),
	  pnl->error_message,
	  pnl->error_message);
      class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,q,&p_12q),
	  pnl->error_message,
	  pnl->error_message);
      class_call(trg_p_ab_at_any_k(pnl,pnl->p_12,pnl->ddp_12,0,p,&p_12p),
	  pnl->error_message,
	  pnl->error_message);

      class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,q,&p_22q),
	  pnl->error_message,
	  pnl->error_message);
      class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,p,&p_22p),
	  pnl->error_message,
	  pnl->error_message);


      class_call(trg_gamma_222(q,p,k,&gamma2_qpk,pnl->error_message),
	  pnl->error_message,
	  pnl->error_message);

      class_call(trg_gamma_222(p,k,q,&gamma2_pkq,pnl->error_message),
	  pnl->error_message,
	  pnl->error_message);

      trg_gamma_121(k,q,p,&gamma1_kqp);
      trg_gamma_121(k,p,q,&gamma1_kpq);

      *result = gamma1_kqp*p_22q*p_12p + gamma1_kpq*p_12q*p_22p +
	gamma2_qpk*p_22p*p_12k + gamma2_pkq*p_12k*p_22q;

      return _SUCCESS_;
      break;

    /****************************************/
    case _222_:

      class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,k,&p_22k),
	  pnl->error_message,
	  pnl->error_message);
      class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,q,&p_22q),
	  pnl->error_message,
	  pnl->error_message);
      class_call(trg_p_ab_at_any_k(pnl,pnl->p_22,pnl->ddp_22,0,p,&p_22p),
	  pnl->error_message,
	  pnl->error_message);


      class_call(trg_gamma_222(k,q,p,&gamma2_kqp,pnl->error_message),
	  pnl->error_message,
	  pnl->error_message);
      class_call(trg_gamma_222(q,p,k,&gamma2_qpk,pnl->error_message),
	  pnl->error_message,
	  pnl->error_message);
      class_call(trg_gamma_222(p,k,q,&gamma2_pkq,pnl->error_message),
	  pnl->error_message,
	  pnl->error_message);

      *result = gamma2_kqp*p_22q*p_22p + gamma2_qpk*p_22p*p_22k + gamma2_pkq*p_22k*p_22q;

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

int trg_I_arg(
    struct spectra_nl * pnl,
    enum name_I name,
    double k,
    double p,
    double m,
    int i,
    int index_k,
    int index_eta,
    double * result
    ){

  /** - define local variables */

  /* shorter names for the interaction functions */

  double gamma1_kpm,gamma1_pkm,gamma1_mkp,gamma1_mpk,gamma1_pmk,gamma1_kmp;
  double gamma2_kpm,gamma2_pkm,gamma2_mkp,gamma2_mpk,gamma2_pmk,gamma2_kmp;

  /*[>* - take care of the cases where one of the argument is below k_min <]*/

  /*if (m<pnl->k[0] || p<pnl->k[0]){*/
    /**result=0.;*/
    /*return _SUCCESS_;*/
  /*}*/

  /** - for each name (_A0_ to _B3_), assign all and only the relevant parameters
      to the local variables, and return the correct argument. */

  /* each case is separated for computation time reason, though it makes the function
     harder to read. For precision : each 'case' statement starts by calling P_11, P_12, P_22
     then gamma_222, then gamma_121. It then computes the function, and finally multiplies
     by the Fourier factor and the p.m factor (x^2-y^2) */


  switch(name){
    case _11_:

      trg_gamma_121(k,p,m,&gamma1_kpm);
      trg_gamma_121(k,m,p,&gamma1_kmp);

      *result=2.*gamma1_kpm*pnl->b_121_nl[index_k+pnl->k_size*(index_eta+pnl->eta_size*i)]+
	2.*gamma1_kmp*pnl->b_121_nl[index_k+pnl->k_size*(index_eta+pnl->eta_size*2*i)];
      return _SUCCESS_;
      break;

      /****************************************/

    case _12_:
      class_call(trg_gamma_222(k,p,m,&gamma2_kpm,pnl->error_message),
	  pnl->error_message,
	  pnl->error_message);

      class_call(trg_gamma_222(k,m,p,&gamma2_kmp,pnl->error_message),
	  pnl->error_message,
	  pnl->error_message);

      trg_gamma_121(k,p,m,&gamma1_kpm);
      trg_gamma_121(k,m,p,&gamma1_kmp);

      *result=gamma1_kpm*pnl->b_221_nl[index_k+pnl->k_size*(index_eta+pnl->eta_size*i)]+
	gamma1_kmp*pnl->b_221_nl[index_k+pnl->k_size*(index_eta+pnl->eta_size*2*i)]+
	gamma2_kpm*pnl->b_122_nl[index_k+pnl->k_size*(index_eta+pnl->eta_size*i)]+
	gamma2_kmp*pnl->b_122_nl[index_k+pnl->k_size*(index_eta+pnl->eta_size*2*i)];

      return _SUCCESS_;
      break;

      /****************************************/

    case _22_:

      class_call(trg_gamma_222(k,p,m,&gamma2_kpm,pnl->error_message),
	  pnl->error_message,
	  pnl->error_message);

      class_call(trg_gamma_222(k,m,p,&gamma2_kmp,pnl->error_message),
	  pnl->error_message,
	  pnl->error_message);
      *result=gamma2_kpm*pnl->b_222_nl[index_k+pnl->k_size*(index_eta+pnl->eta_size*i)]+
	gamma2_kmp*pnl->b_222_nl[index_k+pnl->k_size*(index_eta+pnl->eta_size*2*i)];
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
			    enum name_I name,
			    int index_eta,
			    double * result
			    ){

 /** - define local variables */

  int index_k;
  double k;
  double k_max;

  double arg;

  double xl,xr,yl,yr;
  double p,m;

  double *w;
  double *x;
  double *y;

  int n;
  int i;
  /** - set a value for the upper bound in x, hard-coded. */

  k_max=3000.;

  /** - enter the loop over k-values still under the 'double_escape'
        condition */

  for(index_k=0; index_k<pnl->k_size-pnl->double_escape*(2*index_eta+1); index_k++){
    k=pnl->k[index_k];

    /** determine boundary of the integral to perform over xy */

    xl=k/sqrt(2.);
    xr=k_max;
    yl=0.;
    yr=k/sqrt(2.);

    /** implement quadrature integral. The result will be the weighted sum of the bispectra
     * evaluated in carefully chosen points */

    n=24;

    class_call(quadrature_in_rectangle(xl,xr,yl,yr,&n,&x,&y,&w,pnl->error_message),
	pnl->error_message,
	pnl->error_message);

    for(i=0;i<n;i++){

      p=(x[i]+y[i])/sqrt(2.);
      m=(x[i]-y[i])/sqrt(2.);

      class_call(trg_I_arg(pnl,name,k,p,m,i,index_k,index_eta,&arg),
	  pnl->error_message,
	  pnl->error_message);

      *result+=((x[i])*(x[i])-y[i]*y[i])/2.*arg*w[i];
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
   * 		For each step in time, it first computes the non-linear quantities called I in this
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
  int temp_index_k;

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

  double cutoff;

  double time_step;
  int index_int;

  /** Background quantities */

  double * pvecback_nl;
  double * Omega_m, * H, *H_prime;
  double * growth_factor;
  double * corrected_growth_factor;

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

  int index_name,name_size;
  double **G;

  /** Additionnal temporary variables */

  double temp;
  double * junk;
  double pk,pk_ini;

  double dtau;
  double dz_p;
  double delta_minus,delta_plus,d_delta_m_over_dz;

  int i,n;
  n=24;

  int index_b,index_b_plus;

  double *x,*y,*w;
  double xl,xr,yl,yr,k_max;
  double pl,pr,ql,qr;
  k_max=3000;
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
  temp_index_k=index_k;

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
				    long_info,
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

  /* Implementation of interpolation for Omega matrices elements.*/

  double *ddOmega_21,*ddOmega_22;

  class_calloc(ddOmega_21, pnl->eta_size*pnl->k_size,sizeof(double),pnl->error_message);
  class_calloc(ddOmega_22, pnl->eta_size*pnl->k_size,sizeof(double),pnl->error_message);



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

    /* There is the possibility to add a cutoff to take into account
     * exponential suppression of the power spectrum at high k. By
     * default, however, set on 1. (no cutoff).
     */

    cutoff=1.;

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

  /**
   * Definition of xy size for quadrature integration is given
   * in input file by the parameter n. If it is incorrect, it will
   * changed (always to the closest lower value). This way, the
   * initialization can take place now.
   */

  pnl->xy_size=24; // would be xy_size=ppr->n;

  /**
   * Definition of B_111,B_121,B_122,B_222, the four independent three points
   * correlators of the density/veloctiy fields, initialized at eta=0 (by default,
   * z=35)
   */

  class_calloc(pnl->b_111_nl,pnl->k_size*pnl->eta_size*pnl->xy_size*2,sizeof(double),pnl->error_message);
  class_calloc(pnl->b_121_nl,pnl->k_size*pnl->eta_size*pnl->xy_size*2,sizeof(double),pnl->error_message);
  class_calloc(pnl->b_211_nl,pnl->k_size*pnl->eta_size*pnl->xy_size*2,sizeof(double),pnl->error_message);
  class_calloc(pnl->b_122_nl,pnl->k_size*pnl->eta_size*pnl->xy_size*2,sizeof(double),pnl->error_message);
  class_calloc(pnl->b_221_nl,pnl->k_size*pnl->eta_size*pnl->xy_size*2,sizeof(double),pnl->error_message);
  class_calloc(pnl->b_222_nl,pnl->k_size*pnl->eta_size*pnl->xy_size*2,sizeof(double),pnl->error_message);

  class_calloc(pnl->b_111,pnl->k_size*pnl->xy_size*2,sizeof(double),pnl->error_message);
  class_calloc(pnl->b_121,pnl->k_size*pnl->xy_size*2,sizeof(double),pnl->error_message);
  class_calloc(pnl->b_122,pnl->k_size*pnl->xy_size*2,sizeof(double),pnl->error_message);
  class_calloc(pnl->b_222,pnl->k_size*pnl->xy_size*2,sizeof(double),pnl->error_message);

  name_size = _222_+1;
  class_calloc(G,name_size,sizeof(double*),pnl->error_message);
  for (index_name=0; index_name<name_size; index_name++)
    class_calloc(G[index_name],pnl->k_size*pnl->eta_size*pnl->xy_size*2,sizeof(double),pnl->error_message);

  int name_size_I;
  double **I;
  name_size_I = _22_+1;
  class_calloc(I,name_size_I,sizeof(double*),pnl->error_message);
  for (index_name=0; index_name<name_size_I; index_name++)
    class_calloc(I[index_name],pnl->k_size*pnl->eta_size,sizeof(double),pnl->error_message);

  /** Initialization of bispectra
   */

  //TO DO


  /** First step in time, computing the non-linear integrals over bispectrum I*/

  if (pnl->spectra_nl_verbose > 1)
    printf(" -> initialisation\n");

  if(pnl->mode > 0){

    /* initialize error management flag */
    abort = _FALSE_;

    /*** beginning of parallel region ***/

#pragma omp parallel							\
  shared(name_size_I,abort,pba,ppm,psp,pnl,I)				\
  private(tstart,index_name,tstop)

    {

#ifdef _OPENMP
      tstart = omp_get_wtime();
#endif

#pragma omp for schedule (dynamic)
      for (index_name=0; index_name<name_size_I; index_name++) {

#pragma omp flush(abort)

	class_call_parallel(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,index_name,0,I[index_name]),
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

  /** Now we calculate the time evolution with a predictor corrector
      algorithm. */

  /* At each step, first the new power spectra are computed on a half
   * time step (and pnl->double_escape points are dropped) then the I
   * functions are updated (for 1-loop method, they are just copied
   * from last step) (and pnl->double_escape points are dropped).
   * Then, the new power spectra are evaluated in the whole time step
   * using quantities computed on the half time step, and the new I
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

    double p,q;
    double omega_21_p,omega_21_q,omega_22_p,omega_22_q;
    int index_2b;

    for (index_k=0; index_k<pnl->k_size-pnl->double_escape*2*(index_eta-1)/1; index_k++){

      /* Some useful intermediate variables */

      fourpi_over_k=4*_PI_/(pnl->k[index_k]);
      index = index_k+pnl->k_size*(index_eta-2);
      index_plus = index_k+pnl->k_size*(index_eta-1);

     /* For each k, compute the power spectra and bispectra (through ai) at the new time */

      pnl->p_11_nl[index_plus]= (time_step) *(
					 -2*Omega_11*pnl->p_11_nl[index]
					 -2*Omega_12*pnl->p_12_nl[index]
					 +exp_eta*2*fourpi_over_k*(I[_11_][0]))
	+ pnl->p_11_nl[index];

      pnl->p_22_nl[index_plus] = (time_step) *(
					    -2*Omega_22[index]*pnl->p_22_nl[index]
					    -2*Omega_21[index]*pnl->p_12_nl[index]
					    +exp_eta*2*fourpi_over_k*(I[_22_][0]))
	+ pnl->p_22_nl[index];

      pnl->p_12_nl[index_plus] = (time_step) *(
					    -pnl->p_12_nl[index]*(Omega_11+Omega_22[index])
					    -Omega_12*pnl->p_22_nl[index]
					    -Omega_21[index]*pnl->p_11_nl[index]
					    +exp_eta*fourpi_over_k*I[_12_][0])
	+ pnl->p_12_nl[index];


      for(i=0;i<2*n;i++){

	index_b = index_k+pnl->k_size*(index_eta-2+pnl->eta_size*i);
	if(i<n){
	  index_2b = index_k+pnl->k_size*(index_eta-2+pnl->eta_size*(i+n));
	}
	else
	  index_2b = index_k+pnl->k_size*(index_eta-2+pnl->eta_size*(i-n));
	index_b_plus = index_k+pnl->k_size*(index_eta-1+pnl->eta_size*i);
	printf("%d %d %d %d\n",pnl->k_size*pnl->eta_size*pnl->xy_size*2,index_b,index_2b,index_b_plus);
	xl=pnl->k[index_k]/sqrt(2.);
	xr=k_max;
	yl=0.;
	yr=pnl->k[index_k]/sqrt(2.);

	pl=(xl-yl)/sqrt(2.);
	pr=(xr-yr)/sqrt(2.);
	ql=(xl+yl)/sqrt(2.);
	qr=(xr+yr)/sqrt(2.);

	class_call(quadrature_in_rectangle(pl,pr,ql,qr,&n,&p,&q,&w,pnl->error_message),
	    pnl->error_message,
	    pnl->error_message);

	for(index_name=0;index_name<name_size;index_name++){
	  class_call(trg_G_terms(pnl,index_name,pnl->k[index_k],q,p,index_eta,G[index_name]),
	      pnl->error_message,
	      pnl->error_message);
	}

	class_call(array_interpolate_extrapolate_spline_one_column(pnl->k,pnl->k_size,Omega_21,pnl->eta_size,
	      index_eta,ddOmega_21,p,&omega_21_p,pnl->error_message),
	    pnl->error_message,
	    pnl->error_message);

	class_call(array_interpolate_extrapolate_spline_one_column(pnl->k,pnl->k_size,Omega_21,pnl->eta_size,
	      index_eta,ddOmega_21,q,&omega_21_q,pnl->error_message),
	    pnl->error_message,
	    pnl->error_message);

	class_call(array_interpolate_extrapolate_spline_one_column(pnl->k,pnl->k_size,Omega_22,pnl->eta_size,
	      index_eta,ddOmega_22,p,&omega_22_p,pnl->error_message),
	    pnl->error_message,
	    pnl->error_message);

	class_call(array_interpolate_extrapolate_spline_one_column(pnl->k,pnl->k_size,Omega_22,pnl->eta_size,
	      index_eta,ddOmega_22,q,&omega_22_q,pnl->error_message),
	    pnl->error_message,
	    pnl->error_message);

	pnl->b_111_nl[index_b_plus] = (time_step) *(
	    -pnl->b_111_nl[index_b]*(3.*Omega_11 ) - pnl->b_211_nl[index_b]*Omega_12 - (pnl->b_121_nl[index_b]+pnl->b_121_nl[index_2b])*Omega_12 + 2*exp_eta*G[_111_][index_b]) + pnl->b_111_nl[index_b];

	pnl->b_121_nl[index_b_plus] = (time_step) *( -pnl->b_121_nl[index_b]*(Omega_11+omega_22_q+omega_22_p)-pnl->b_222_nl[index_b]*Omega_12
	    -pnl->b_121_nl[index_b]*omega_21_p - pnl->b_121_nl[index_2b]*omega_21_q + 2.*exp_eta*G[_121_][index_b]) + pnl->b_121_nl[index_b];

	pnl->b_211_nl[index_b_plus] = (time_step) *( -pnl->b_211_nl[index_b]*(Omega_22[index] + 2.*Omega_11)-pnl->b_111_nl[index_b]*Omega_21[index]
	    -pnl->b_221_nl[index_b]*Omega_12 - pnl->b_221_nl[index_2b]*Omega_12 + 2.*exp_eta*G[_211_][index_b]) + pnl->b_211_nl[index_b];

	pnl->b_122_nl[index_b_plus] = (time_step) *(
	    -pnl->b_122_nl[index_b]*(Omega_11+omega_22_q+omega_22_p) - pnl->b_222_nl[index_b]*Omega_12-pnl->b_121_nl[index_b]*omega_21_p
	    -pnl->b_121_nl[index_2b]*omega_21_q + 2.*exp_eta*G[_122_][index_b]) + pnl->b_122_nl[index_b];

	pnl->b_221_nl[index_b_plus] = (time_step) *(
	    -pnl->b_221_nl[index_b]*(Omega_22[index]+omega_22_q+Omega_12) - pnl->b_121_nl[index_b]*Omega_21[index]-pnl->b_211_nl[index_b]*omega_21_q
	    -pnl->b_222_nl[index_b]*Omega_12 + 2.*exp_eta*G[_221_][index_b]) + pnl->b_221_nl[index_b];

	pnl->b_222_nl[index_b_plus] = (time_step) *(
	    -pnl->b_222_nl[index_b]*(Omega_21[index]+omega_22_q+omega_22_p)-pnl->b_122_nl[index_b]*Omega_21[index]-pnl->b_221_nl[index_2b]*omega_21_q
	    -pnl->b_221_nl[index_b]*omega_21_p + 2.*exp_eta*G[_222_][index_b]) + pnl->b_222_nl[index_b];

      }

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

    /*if(pnl->mode==1){*/
      /*for (index_name=0; index_name<name_size; index_name++){*/
	/*if(index_name==0)*/
	  /*a=0;*/
	/*else if(index_name==1||index_name==2||index_name==3||index_name==11||index_name==12)*/
	  /*a=1;*/
	/*else if(index_name==4||index_name==5||index_name==6||index_name==9||index_name==10)*/
	  /*a=2;*/
	/*else if(index_name==7||index_name==8)*/
	  /*a=3;*/
	/*else*/
	  /*a=4;*/
	/*b=4-a;*/
	/*for(index_k=0; index_k<pnl->k_size-pnl->double_escape*2*(index_eta)/1; index_k++){*/
	  /*AA[index_name][index_k+pnl->k_size*(index_eta-1)]=pow(growth_factor[index_eta-1],a)*pow(corrected_growth_factor[index_eta-1],b)*AA[index_name][index_k];}*/
      /*}*/
    /*}*/

    /* For TRG, simply recomputes the AA values entirely */

    /*else if(pnl->mode == 2){*/

      /*[> initialize error management flag <]*/
      /*abort = _FALSE_;*/

      /*[>** beginning of parallel region **<]*/

/*#pragma omp parallel							\*/
      /*shared(name_size,abort,pba,ppm,psp,pnl,index_eta,AA)		\*/
      /*private(tstart,index_name,tstop)*/

      /*{*/

/*#ifdef _OPENMP*/
	/*tstart = omp_get_wtime();*/
/*#endif*/

/*#pragma omp for schedule (dynamic)*/
	/*for (index_name=0; index_name<name_size; index_name++) {*/

/*#pragma omp flush(abort)*/

	  /*class_call_parallel(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,index_name,index_eta-1,AA[index_name]),*/
	      /*pnl->error_message,*/
	      /*pnl->error_message);*/
	/*}*/


/*#ifdef _OPENMP*/
	/*tstop = omp_get_wtime();*/
	/*if ((pnl->spectra_nl_verbose > 2) && (pnl->mode > 1))*/
	  /*printf("In %s: time spent in parallel region (loop over names) = %e s for thread %d\n",*/
	      /*__func__,tstop-tstart,omp_get_thread_num());*/
/*#endif*/

      /*} [> end of parallel region <]*/

      /*if (abort == _TRUE_) return _FAILURE_;*/

    /*}*/

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
					 +exp_eta*4*fourpi_over_k)
	+ pnl->p_11_nl[index_int];

      pnl->p_22_nl[index_plus] = (time_step) *(
					    -2*Omega_22[index]*pnl->p_22_nl[index]
					    -2*Omega_21[index]*pnl->p_12_nl[index]
					    +exp_eta*2*fourpi_over_k )
	+ pnl->p_22_nl[index_int];

      pnl->p_12_nl[index_plus] = (time_step) *(
					    -pnl->p_12_nl[index]*(Omega_11+Omega_22[index])
					    -Omega_12*pnl->p_22_nl[index]
					    -Omega_21[index]*pnl->p_11_nl[index]
					    +exp_eta*fourpi_over_k/*(2*a13[index]+b21[index])*/)
	+ pnl->p_12_nl[index_int];


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

    /*if(pnl->mode==1){*/
      /*for (index_name=0; index_name<name_size; index_name++){*/
	/*if(index_name==0)*/
	  /*a=0;*/
	/*else if(index_name==1||index_name==2||index_name==3||index_name==11||index_name==12)*/
	  /*a=1;*/
	/*else if(index_name==4||index_name==5||index_name==6||index_name==9||index_name==10)*/
	  /*a=2;*/
	/*else if(index_name==7||index_name==8)*/
	  /*a=3;*/
	/*else*/
	  /*a=4;*/
	/*b=4-a;*/
	/*for(index_k=0; index_k<pnl->k_size-pnl->double_escape*2*(index_eta)/1; index_k++){*/
	  /*AA[index_name][index_k+pnl->k_size*index_eta]=pow(growth_factor[index_eta-1],a)*pow(corrected_growth_factor[index_eta-1],b)*AA[index_name][index_k];}*/
      /*}*/
    /*}*/

    /*else if(pnl->mode == 2){*/

      /*[> initialize error management flag <]*/
      /*abort = _FALSE_;*/

      /*[>** beginning of parallel region **<]*/

/*#pragma omp parallel							\*/
      /*shared(name_size,abort,pba,ppm,psp,pnl,index_eta,AA)			\*/
      /*private(tstart,index_name,tstop)*/

      /*{*/

/*#ifdef _OPENMP*/
	/*tstart = omp_get_wtime();*/
/*#endif*/

/*#pragma omp for schedule (dynamic)*/
	/*for (index_name=0; index_name<name_size; index_name++) {*/

/*#pragma omp flush(abort)*/

	  /*class_call_parallel(trg_integrate_xy_at_eta(pba,ppm,psp,pnl,index_name,index_eta,AA[index_name]),*/
	      /*pnl->error_message,*/
	      /*pnl->error_message);*/
	/*}*/


/*#ifdef _OPENMP*/
	/*tstop = omp_get_wtime();*/
	/*if ((pnl->spectra_nl_verbose > 2) && (pnl->mode > 1))*/
	  /*printf("In %s: time spent in parallel region (loop over names) = %e s for thread %d\n",*/
	      /*__func__,tstop-tstart,omp_get_thread_num());*/
/*#endif*/

      /*} [> end of parallel region <]*/

      /*if (abort == _TRUE_) return _FAILURE_;*/

    /*}*/

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

  free(p_11_linear);
  free(p_12_linear);
  free(p_22_linear);

  free(Omega_21);
  free(Omega_22);

  free(ddOmega_21);
  free(ddOmega_22);

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
