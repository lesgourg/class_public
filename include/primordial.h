/** @file primordial.h Documented includes for primordial module */

#ifndef __PRIMORDIAL__
#define __PRIMORDIAL__

#include "transfer.h"

enum primordial_spectrum_type {
  analytic_Pk,
  inflation_V,
  inflation_H
};

/**
 * All tables of transfer functions \f$ \Delta_l^{X} (k) \f$.
 *
 * Once initialized by transfer_init(), contains all tables of
 * transfer functions used for interpolation in other modules, for all
 * requested modes (scalar/vector/tensor), initial conditions, type
 * (temperature, polarization, etc), l and k.
 */
struct primordial {

  int md_size; /**< number of modes included in computation */
  int * ic_size;       /**< for a given mode, ic_size[index_mode] = number of initial conditions included in computation */
  int * ic_ic_size;    /**< number of pairs of (index_ic1, index_ic2) with index_ic2 >= index_ic1; this number is just ic_size[index_mode](ic_size[index_mode]+1)/2  */

  enum primordial_spectrum_type primordial_spec_type;

  double k_pivot; /* pivot scale in Mpc-1 */

  double A_s;  /* usual scalar amplitude = curvature power spectrum at pivot scale */
  double n_s;  /* usual scalar tilt = [curvature power spectrum tilt at pivot scale -1] */
  double alpha_s; /* usual scalar running */

  double r;    /* usual tensor to scalar ratio of power spectra, r=A_T/A_S=P_h/P_R */
  double n_t;  /* usual tensor tilt = [GW power spectrum tilt at pivot scale] */
  double alpha_t; /* usual tensor running */

  double f_bi;  /* baryon isocurvature (BI) entropy-to-curvature ratio S_bi/R */
  double n_bi;  /* BI tilt */
  double alpha_bi; /* BI running */

  double f_cdi;  /* CDM isocurvature (CDI) entropy-to-curvature ratio S_cdi/R */
  double n_cdi;  /* CDI tilt */
  double alpha_cdi; /* CDI running */

  double f_nid;  /* neutrino density isocurvature (NID) entropy-to-curvature ratio S_nid/R */
  double n_nid;  /* NID tilt */
  double alpha_nid; /* NID running */

  double f_niv;  /* neutrino velocity isocurvature (NIV) entropy-to-curvature ratio S_niv/R */
  double n_niv;  /* NIV tilt */
  double alpha_niv; /* NIV running */

  double c_ad_bi; /* ADxBI cross-correlation at pivot scale, from -1 to 1 */
  double n_ad_bi; /* ADxBI cross-correlation tilt */
  double alpha_ad_bi; /* ADxBI cross-correlation running */

  double c_ad_cdi; /* ADxCDI cross-correlation at pivot scale, from -1 to 1 */
  double n_ad_cdi; /* ADxCDI cross-correlation tilt */
  double alpha_ad_cdi; /* ADxCDI cross-correlation running */

  double c_ad_nid; /* ADxNID cross-correlation at pivot scale, from -1 to 1 */
  double n_ad_nid; /* ADxNID cross-correlation tilt */
  double alpha_ad_nid; /* ADxNID cross-correlation running */

  double c_ad_niv; /* ADxNIV cross-correlation at pivot scale, from -1 to 1 */
  double n_ad_niv; /* ADxNIV cross-correlation tilt */
  double alpha_ad_niv; /* ADxNIV cross-correlation running */

  double c_bi_cdi; /* BIxCDI cross-correlation at pivot scale, from -1 to 1 */
  double n_bi_cdi; /* BIxCDI cross-correlation tilt */
  double alpha_bi_cdi; /* BIxCDI cross-correlation running */

  double c_bi_nid; /* BIxNIV cross-correlation at pivot scale, from -1 to 1 */
  double n_bi_nid; /* BIxNIV cross-correlation tilt */
  double alpha_bi_nid; /* BIxNIV cross-correlation running */

  double c_bi_niv; /* BIxNIV cross-correlation at pivot scale, from -1 to 1 */
  double n_bi_niv; /* BIxNIV cross-correlation tilt */
  double alpha_bi_niv; /* BIxNIV cross-correlation running */

  double c_cdi_nid; /* CDIxNID cross-correlation at pivot scale, from -1 to 1 */
  double n_cdi_nid; /* CDIxNID cross-correlation tilt */
  double alpha_cdi_nid; /* CDIxNID cross-correlation running */

  double c_cdi_niv; /* CDIxNIV cross-correlation at pivot scale, from -1 to 1 */
  double n_cdi_niv; /* CDIxNIV cross-correlation tilt */
  double alpha_cdi_niv; /* CDIxNIV cross-correlation running */

  double c_nid_niv; /* NIDxNIV cross-correlation at pivot scale, from -1 to 1 */
  double n_nid_niv; /* NIDxNIV cross-correlation tilt */
  double alpha_nid_niv; /* NIDxNIV cross-correlation running */

  /* above parameters are stored more conveniently in symmetric matrices */
  double ** amplitude; /* amplitude[index_mode][index_ic1_ic2] */
  double ** tilt;      /* tilt[index_mode][index_ic1_ic2] */
  double ** running;   /* running[index_mode][index_ic1_ic2] */
  
  double * lnk; /* list of ln(k) values lnk[index_k] */
  int lnk_size; /* number of ln(k) values */

  double ** lnpk; /* depends on indices index_mode, index_ic1, index_ic2, index_k as:
		     lnpk[index_mode][index_k*ppm->ic_ic_size[index_mode]+index_ic1_ic2]
		     where index_ic1_ic2 labels ordered pairs (index_ic1, index_ic2) (since 
                     the primordial spectrum is symmetric in (index_ic1, index_ic2)).
		     - for diagonal elements (index_ic1 = index_ic2) this arrays contains
                       ln[P(k)] where P(k) is positive by construction.
		     - for non-diagonal elements this arrays contains the k-dependent 
		       cosine of the correlation angle, namely
		       P(k)_(index_ic1, index_ic2)/sqrt[P(k)_index_ic1 P(k)_index_ic2]
		     This choice is convenient since the sign of the non-diagonal cross-correlation 
		     is arbitrary. For fully correlated or anti-correlated initial conditions,
		     this non-diagonal element is independent on k, and equal to +1 or -1.
		  */
  double ** ddlnpk; /* second derivative of above array, for spline interpolation. So: 
		       - for index_ic1 = index_ic), we spline ln[P(k)] vs. ln(k), which is
		         good since this function is usually smooth.
		       - for non-diagonal coefficients, we spline  
		         P(k)_(index_ic1, index_ic2)/sqrt[P(k)_index_ic1 P(k)_index_ic2]
			 vs. ln(k), which is fine since this quantity is often assumed to be
                         constant (e.g for fully correlated/anticorrelated initial conditions)
			 or nearly constant, and with arbitrary sign.
		    */
  short ** is_non_zero; /* is_non_zero[index_mode][index_ic1_ic2] set to false if pair
			    (index_ic1, index_ic2) is uncorrelated 
			    (ensures more precision and saves time with respect to the option
			    of simply setting P(k)_(index_ic1, index_ic2) to zero) */

  /** @name - flag regulating the amount of information sent to standard output (none if set to zero) */

  //@{

  short primordial_verbose;

  //@}

  ErrorMsg error_message; /**< zone for writing error messages */

};

/*************************************************************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int primordial_spectrum_at_k(
			       struct primordial * ppm,
			       int index_mode,
			       double k,
			       double * pk
			       );

  int primordial_init(
		      struct precision  * ppr,
		      struct perturbs   * ppt,
		      struct primordial * ppm
		      );
    
  int primordial_free(
		      struct primordial * ppm
		      );
    
  int primordial_indices(
			 struct perturbs   * ppt,
			 struct primordial * ppm
			 );

  int primordial_analytic_spectrum_init(
					struct perturbs   * ppt,
					struct primordial * ppm
					);

  int primordial_analytic_spectrum(
				   struct primordial * ppm,
				   int index_mode,
				   int index_ic1_ic2,
				   double k,
				   double * pk
				   );

  int primordial_get_lnk_list(
			      struct primordial * ppm,
			      double kmin,
			      double kmax,
			      double k_per_decade
			      );

#ifdef __cplusplus
}
#endif

/**  
 * @name Constants
 */

//@{


//@}

#endif
