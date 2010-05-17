#include "spectra.h"

#ifndef __SPECTRA_NL__
#define __SPECTRA_NL__

struct spectra_nl {

  double * k;
  double * eta;
  double * pk_nl;
  double * p_12;
  double * p_22;
  double * ddpk_nl;
  double * ddp_12;
  double * ddp_22;
  double * z;

  double k_max;

  double z_ini;
  double eta_step;
  int eta_size;
  int k_size;

  short spectra_nl_verbose;  /**< from 0 to 1: amount of information written in standard output */

  ErrorMsg error_message; /**< zone for writing error messages */
};

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  enum name_A{
    A0,
    A11,
    A12,
    A21,
    A22,
    A3
  };

  int trg_gamma_121(
		    double k, 
		    double q, 
		    double p, 
		    double * result
		    );
  
  int trg_gamma_222(
		    double k, 
		    double p, 
		    double q,
		    double * result,
		    char * errmsg
		    );

  int trg_p12_ini(
		  double *k,
		  int index_ic,
		  double *H,
		  double * result,
		  char * errmsg				
		  );

  int trg_p22_ini(
		  double *k,
		  int index_ic,
		  double *H,
		  double * result,
		  char * errmsg				
		  );
  
  int trg_pk_nl_ini(
		    double *k,
		    int index_ic,
		    double * result,
		    char * errmsg				
		    );
  
  int trg_ddp_ab(
		 double *p_ab,
		 int index_eta,
		 double *result,
		 char * errmsg
		 );
		 

  int trg_p_ab_at_any_k(
			double * p_ab,
			double * ddp_ab,
			int index_eta,
			double any_k,
			double * result,
			char *errmsg
			);


  /********************
   * Argument of the integral of all A's, B's, given at any time eta,
   * usually written as F( k,x+y/sqrt(2),x-y/sqrt(2) ) + F( y <-> -y)
   *
   *
   ********************/

  int trg_A_arg(
		enum name_A name, 
		double k, 
		double p, 
		double m, 
		int index_eta,
		double * result, 
		char * errmsg);



  /**************
   * Function that performs the integration with a simple trapeze
   * method over x and y with number of steps n_xy of any
   * A0..., B0..., at any chosen time eta=eta[index_eta]
   **************/

  int trg_integrate_xy_at_eta(
			      enum name_A name,
			      int index_eta,
			      int n_xy,
			      int k_size, /* functions of k for
					     integration are defined
					     from [0] to [k_size-1]*/
			      double * result,
			      char * errmsg
			      );


  int trg_init(
	       struct precision *ppr_input,
	       struct background *pba_input,
	       struct primordial *ppm_input,
	       struct spectra *psp_input,
	       struct spectra_nl *pnl_output
	       );

  int trg_free();
      
#ifdef __cplusplus
}
#endif

/**** below some define commands ***/

#endif
