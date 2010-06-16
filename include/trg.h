#include "spectra.h"

#ifndef __SPECTRA_NL__
#define __SPECTRA_NL__

struct spectra_nl {

  double * k;
  int k_size; /**< total number of k values */
  int index_k_L;  /**< for index_k=0, ...,(index_k_L-1), use linear theory only */

  double * eta;
  int eta_size;
  double eta_step;

  double * z;
  double z_ini;

  double * pk_nl;
  double * p_12_nl;
  double * p_22_nl;

  double * p_11;
  double * p_12;
  double * p_22;

  double * ddpk_nl;
  double * ddp_12_nl;
  double * ddp_22_nl;

  double * ddp_11;
  double * ddp_12;
  double * ddp_22;

  short spectra_nl_verbose;  /**< from 0 to 1: amount of information written in standard output */
  short mode; /**< from 0 to 2: 0 being linear theory, 1 for one loop and 2 for full trg calculation*/

  ErrorMsg error_message; /**< zone for writing error messages */
};

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  enum name_A{
    _A0_,
    _A11_,
    _A12_,
    _A13_,
    _A21_,
    _A22_,
    _A23_,
    _A3_,
    _B0_,
    _B11_,
    _B12_,
    _B21_,
    _B22_,
    _B3_
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

  int trg_p12_at_k(
	      struct background * pba,
	      struct primordial * ppm,
	      struct spectra * psp,
	      struct spectra_nl * pnl,
	      int index_eta,
	      int index_ic,
	      double k,
	      double * result
	      );

  int trg_p22_at_k(
	      struct background * pba,
	      struct primordial * ppm,
	      struct spectra * psp,
	      struct spectra_nl * pnl,
	      int index_eta,
	      int index_ic,
	      double k,
	      double * result
	      );
  
  int trg_p11_at_k(
	      struct background * pba,
	      struct primordial * ppm,
	      struct spectra * psp,
	      struct spectra_nl * pnl,
	      int index_eta,
	      int index_ic,
	      double k,
	      double *result				
	      );
  
  int trg_ddp_ab(
		 struct spectra_nl * pnl,
		 double *p_ab,
		 int index_eta,
		 double *result,
		 char * errmsg
		 );
		 

  int trg_p_ab_at_any_k(
			struct spectra_nl * pnl,
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
		struct spectra_nl *pnl,
		enum name_A name, 
		double k, 
		double p, 
		double m, 
		int index_eta,
		int index_k, /**< used only for testing, could be supressed */
		double * result, 
		char * errmsg);

  int trg_A_arg_one_loop(
			 struct spectra_nl * pnl,
			 enum name_A name, 
			 double k, 
			 double p, 
			 double m, 
			 int index_eta,
			 int index_k, /**< used only for testing, could be supressed */
			 double * result, 
			 char * errmsg);



  /**************
   * Function that performs the integration with a simple trapeze
   * method over x and y with number of steps n_xy of any
   * A0..., B0..., at any chosen time eta=eta[index_eta]
   **************/

  int trg_integrate_xy_at_eta(
			      struct background * pba,
			      struct primordial * ppm,
			      struct spectra * psp,
			      struct spectra_nl * pnl,
			      enum name_A name,
			      int index_eta,
			      double * result,
			      char * errmsg
			      );


  int trg_init(
	       struct precision *ppr,
	       struct background *pba,
	       struct primordial *ppm,
	       struct spectra *psp,
	       struct spectra_nl *pnl
	       );

  int trg_free(
	       struct spectra_nl *pnl
	       );
      
#ifdef __cplusplus
}
#endif

/**** below some define commands ***/

#define _STOP_INT_ 1.e-3

#endif
