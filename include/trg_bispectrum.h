/** @file trg.h Document includes for trg module */

#include "spectra.h"

#ifndef __SPECTRA_NL__
#define __SPECTRA_NL__

/**
 * Structure containing all the non-linear information.
 *
 * Once initialised by spectra_init(), contains a table for all two points correlation functions
 * and for all the ai,bj functions (containing the three points correlation functions), for each
 * time and wave-number.
 */
struct spectra_nl {

  /** @name - input parameters initialized by user in input module
      (all other quantities are computed in this module , given these parameters
      and the content of the 'precision', 'background', 'thermo',
      'primordial' and 'spectra' structure */

  //@{

  short mode; 		 	/**< from 0 to 2: 0 being linear theory, 1 for one loop and 2 for full trg calculation*/
  int double_escape;		/**< Usual value is 2 (which means the code drops 2 points every half-step), you might want to increase it a bit if you have a finner k-grid */
  double z_ini;			/**< Starting value of the redshift for the non-linear computation */
  int eta_size;			/**< Number of steps in time, to change in agreement with the precision, and the escape parameter */
  double k_min;			/**< Minimal k (in h/Mpc ? CHECK) where the spectra is computed. k_max is infered from the spectra module */
  double k_L;			/**< k value (in h/MPC ? CHECK) under which we consider the spectra to stay linear for all computed redshifts */
  double logstepx_min;		/**< Set the precision of the xy integrator */

  //@}

  /** @name - table of k values, and related quantities */

  //@{

  double * k;		/**< table containing the values of k used in this module */
  double k_max;		/**< maximum value of k where the spectra is computed, infered from spectra module */
  int k_size; 		/**< total number of k values */
  int xy_size; 		/**< total number of (x,y) values on the integral domain */
  int index_k_L;  	/**< for index_k=0, ...,(index_k_L-1), use linear theory only */
  double k_growth_factor;  /**< value used to define linear growth factor */
  //@}

  /** @name - tables of time values, and related quantities */

  //@{

  double * eta;		/**< table containing eta values used in this module(new time variable defined as log(a/a_ini)) */
  double * z;		/**< table containing z   values used in this module */
  double eta_step;	/**< Mean step in eta (time variable), the actual time steps are defined in trg_init() for a better precision */

  //@}

  /** @name - tables of non-linear spectra and their derivatives*/

  //@{

  double * p_11_nl;
  double * p_12_nl;
  double * p_22_nl;

  double * p_11;
  double * p_12;
  double * p_22;

  double * ddp_11_nl;
  double * ddp_12_nl;
  double * ddp_22_nl;

  double * ddp_11;
  double * ddp_12;
  double * ddp_22;

  //@}

  /** @name - tables of non-linear bispectra*/

  //@{

  double * b_111_nl;
  double * b_121_nl;
  double * b_122_nl;
  double * b_221_nl;
  double * b_211_nl;
  double * b_222_nl;

  double * b_111;
  double * b_121;
  double * b_122;
  double * b_221;
  double * b_211;
  double * b_222;

  //@}

  /** @name - technical parameters */

  //@{

  short spectra_nl_verbose;  	/**< from 0: amount of information written in standard output */

  ErrorMsg error_message; 	/**< zone for writing error messages */

  //@}
};

/********************************************************************************/

/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  /**
   * List of names for the different bispectra terms.
   */

  enum name_B{
    _111_,
    _121_,
    _211_,
    _122_,
    _221_,
    _222_
  };

  enum name_I{
    _11_,
    _12_,
    _22_
  };

  int trg_gamma_121(
		    double k,
		    double p,
		    double q,
		    double * result
		    );

  int trg_gamma_222(
		    double k,
		    double p,
		    double q,
		    double * result,
		    char * errmsg
		    );

  int trg_p_ab_at_any_k(
			struct spectra_nl * pnl,
			double * p_ab,
			double * ddp_ab,
			int index_eta,
			double any_k,
			double * result
			);

  int trg_G_arg(
		struct spectra_nl *pnl,
		enum name_B name,
		double k,
		double p,
		double m,
		int index_eta,
		double * result
		);

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
      );

  int trg_integrate_xy_at_eta(
			      struct background * pba,
			      struct primordial * ppm,
			      struct spectra * psp,
			      struct spectra_nl * pnl,
			      enum name_I name,
			      int index_eta,
			      double * result
			      );

  int trg_logstep1_k(
		    struct precision * ppr,
		    double k,
		    double * logstep);

  int trg_logstep2_k(
		    struct precision * ppr,
		    double k,
		    double * logstep);

  int trg_init(
	       struct precision *ppr,
	       struct background *pba,
	       struct thermo *pth,
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

#endif
