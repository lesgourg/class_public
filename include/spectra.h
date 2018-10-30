/** @file spectra.h Documented includes for spectra module */

#ifndef __SPECTRA__
#define __SPECTRA__

#include "transfer.h"

/**
 * Structure containing everything about anisotropy and Fourier power spectra that other modules need to know.
 *
 * Once initialized by spectra_init(), contains a table of all
 * \f$ C_l\f$'s and P(k) as a function of multipole/wavenumber,
 * mode (scalar/tensor...), type (for \f$ C_l\f$'s: TT, TE...),
 * and pairs of initial conditions (adiabatic, isocurvatures...).
 */

struct spectra {

  /** @name - input parameters initialized by user in input module
      (all other quantities are computed in this module, given these parameters
      and the content of the 'background', 'perturbs', 'transfers' and
      'primordial' structures) */

  //@{

  double z_max_pk;  /**< maximum value of z at which matter spectrum P(k,z) will be evaluated; keep fixed to zero if P(k) only needed today */


  int non_diag; /**< sets the number of cross-correlation spectra
                   that you want to calculate: 0 means only
                   auto-correlation, 1 means only adjacent bins,
                   and number of bins minus one means all
                   correlations */

  //@}

  /** @name - information on number of modes and pairs of initial conditions */

  //@{

  int md_size;           /**< number of modes (scalar, tensor, ...) included in computation */
  int index_md_scalars; /**< index for scalar modes */

  int * ic_size;         /**< for a given mode, ic_size[index_md] = number of initial conditions included in computation */
  int * ic_ic_size;      /**< for a given mode, ic_ic_size[index_md] = number of pairs of (index_ic1, index_ic2) with index_ic2 >= index_ic1; this number is just N(N+1)/2  where N = ic_size[index_md] */
  short ** is_non_zero; /**< for a given mode, is_non_zero[index_md][index_ic1_ic2] is set to true if the pair of initial conditions (index_ic1, index_ic2) are statistically correlated, or to false if they are uncorrelated */

  //@}

  /** @name - information on number of type of C_l's (TT, TE...) */

  //@{

  int has_tt; /**< do we want \f$ C_l^{TT}\f$? (T = temperature) */
  int has_ee; /**< do we want \f$ C_l^{EE}\f$? (E = E-polarization) */
  int has_te; /**< do we want \f$ C_l^{TE}\f$? */
  int has_bb; /**< do we want \f$ C_l^{BB}\f$? (B = B-polarization) */
  int has_pp; /**< do we want \f$ C_l^{\phi\phi}\f$? (\f$ \phi \f$ = CMB lensing potential) */
  int has_tp; /**< do we want \f$ C_l^{T\phi}\f$? */
  int has_ep; /**< do we want \f$ C_l^{E\phi}\f$? */
  int has_dd; /**< do we want \f$ C_l^{dd}\f$? (d = density) */
  int has_td; /**< do we want \f$ C_l^{Td}\f$? */
  int has_pd; /**< do we want \f$ C_l^{\phi d}\f$? */
  int has_ll; /**< do we want \f$ C_l^{ll}\f$? (l = galaxy lensing potential) */
  int has_tl; /**< do we want \f$ C_l^{Tl}\f$? */
  int has_dl; /**< do we want \f$ C_l^{dl}\f$? */

  int index_ct_tt; /**< index for type \f$ C_l^{TT} \f$*/
  int index_ct_ee; /**< index for type \f$ C_l^{EE} \f$*/
  int index_ct_te; /**< index for type \f$ C_l^{TE} \f$*/
  int index_ct_bb; /**< index for type \f$ C_l^{BB} \f$*/
  int index_ct_pp; /**< index for type \f$ C_l^{\phi\phi} \f$*/
  int index_ct_tp; /**< index for type \f$ C_l^{T\phi} \f$*/
  int index_ct_ep; /**< index for type \f$ C_l^{E\phi} \f$*/
  int index_ct_dd; /**< first index for type \f$ C_l^{dd} \f$((d_size*d_size-(d_size-non_diag)*(d_size-non_diag-1)/2) values) */
  int index_ct_td; /**< first index for type \f$ C_l^{Td} \f$(d_size values) */
  int index_ct_pd; /**< first index for type \f$ C_l^{pd} \f$(d_size values) */
  int index_ct_ll; /**< first index for type \f$ C_l^{ll} \f$((d_size*d_size-(d_size-non_diag)*(d_size-non_diag-1)/2) values) */
  int index_ct_tl; /**< first index for type \f$ C_l^{Tl} \f$(d_size values) */
  int index_ct_dl; /**< first index for type \f$ C_l^{dl} \f$(d_size values) */

  int d_size;      /**< number of bins for which density Cl's are computed */

  int ct_size; /**< number of \f$ C_l \f$ types requested */

  //@}

  /** @name - table of pre-computed C_l values, and related quantities */

  //@{

  int * l_size;   /**< number of multipole values for each requested mode, l_size[index_md] */

  int l_size_max; /**< greatest of all l_size[index_md] */

  double * l;    /**< list of multipole values l[index_l] */


  int ** l_max_ct;    /**< last multipole (given as an input) at which
                         we want to output \f$ C_l\f$'s for a given mode and type;
                         l[index_md][l_size[index_md]-1] can be larger
                         than l_max[index_md], in order to ensure a
                         better interpolation with no boundary effects */

  int * l_max;    /**< last multipole (given as an input) at which
                     we want to output \f$ C_l\f$'s for a given mode (maximized over types);
                     l[index_md][l_size[index_md]-1] can be larger
                     than l_max[index_md], in order to ensure a
                     better interpolation with no boundary effects */

  int l_max_tot; /**< last multipole (given as an input) at which
                    we want to output \f$ C_l\f$'s (maximized over modes and types);
                    l[index_md][l_size[index_md]-1] can be larger
                    than l_max[index_md], in order to ensure a
                    better interpolation with no boundary effects */

  double ** cl;   /**< table of anisotropy spectra for each mode, multipole, pair of initial conditions and types, cl[index_md][(index_l * psp->ic_ic_size[index_md] + index_ic1_ic2) * psp->ct_size + index_ct] */
  double ** ddcl; /**< second derivatives of previous table with respect to l, in view of spline interpolation */

  double alpha_II_2_20;	/**< parameter describing adiabatic versus isocurvature contribution in mutipole range [2,20] (see Planck parameter papers) */
  double alpha_RI_2_20;	/**< parameter describing adiabatic versus isocurvature contribution in mutipole range [2,20] (see Planck parameter papers) */
  double alpha_RR_2_20;	/**< parameter describing adiabatic versus isocurvature contribution in mutipole range [2,20] (see Planck parameter papers) */

  double alpha_II_21_200; /**< parameter describing adiabatic versus isocurvature contribution in mutipole range [21,200] (see Planck parameter papers) */
  double alpha_RI_21_200; /**< parameter describing adiabatic versus isocurvature contribution in mutipole range [21,200] (see Planck parameter papers) */
  double alpha_RR_21_200; /**< parameter describing adiabatic versus isocurvature contribution in mutipole range [21,200] (see Planck parameter papers) */

  double alpha_II_201_2500; /**< parameter describing adiabatic versus isocurvature contribution in mutipole range [201,2500] (see Planck parameter papers) */
  double alpha_RI_201_2500; /**< parameter describing adiabatic versus isocurvature contribution in mutipole range [201,2500] (see Planck parameter papers) */
  double alpha_RR_201_2500; /**< parameter describing adiabatic versus isocurvature contribution in mutipole range [201,2500] (see Planck parameter papers) */

  double alpha_II_2_2500; /**< parameter describing adiabatic versus isocurvature contribution in mutipole range [2,2500] (see Planck parameter papers) */
  double alpha_RI_2_2500; /**< parameter describing adiabatic versus isocurvature contribution in mutipole range [2,2500] (see Planck parameter papers) */
  double alpha_RR_2_2500; /**< parameter describing adiabatic versus isocurvature contribution in mutipole range [2,2500] (see Planck parameter papers) */

  double alpha_kp; /**< parameter describing adiabatic versus isocurvature contribution at pivot scale (see Planck parameter papers) */
  double alpha_k1; /**< parameter describing adiabatic versus isocurvature contribution at scale k1 (see Planck parameter papers) */
  double alpha_k2; /**< parameter describing adiabatic versus isocurvature contribution at scale k2 (see Planck parameter papers) */

  //@}

  /** @name - table of pre-computed matter power spectrum P(k) values, and related quantities */

  //@{

  int ln_k_size;    /**< number ln(k) values */
  double * ln_k;    /**< list of ln(k) values ln_k[index_k] */

  int ln_tau_size;  /**< number of ln(tau) values, for the matter
                       power spectrum and the matter transfer
                       functions, (only one if z_max_pk = 0) */

  double * ln_tau;  /**< list of ln(tau) values ln_tau[index_tau], for
                       the matter power spectrum and the matter
                       transfer functions, in growing order. So
                       exp(ln_tau[0]) is the earliest time
                       (i.e. highest redshift), while
                       exp(ln_tau[ln_tau_size-1]) is today (i.e
                       z=0). */

  double * ln_pk;   /**< Matter power spectrum.
                       depends on indices index_ic1_ic2, index_k, index_tau as:
                       ln_pk[(index_tau * psp->k_size + index_k)* psp->ic_ic_size[index_md] + index_ic1_ic2]
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

  double * ddln_pk; /**< second derivative of above array with respect to log(tau), for spline interpolation. So:
                       - for index_ic1 = index_ic, we spline ln[P(k)] vs. ln(k), which is
                       good since this function is usually smooth.
                       - for non-diagonal coefficients, we spline
                       P(k)_(index_ic1, index_ic2)/sqrt[P(k)_index_ic1 P(k)_index_ic2]
                       vs. ln(k), which is fine since this quantity is often assumed to be
                       constant (e.g for fully correlated/anticorrelated initial conditions)
                       or nearly constant, and with arbitrary sign.
                    */

  double sigma8;    /**< sigma8 parameter */

  double sigma8_cb; /**< if ncdm present: contribution to sigma8 from only baryons and cdm */

  double * ln_pk_l;   /**q< Total linear matter power spectrum, just
                           depending on indices index_k, index_tau as:
                           ln_pk[index_tau * psp->k_size + index_k]
                           Range of k and tau value identical to
                           ln_pk array. */

  double * ddln_pk_l; /**< second derivative of above array with respect to log(tau), for spline interpolation. */

  int ln_tau_nl_size;  /**< number of ln(tau) values for non-linear
                          spectrum (possibly smaller than ln_tau_size,
                          because the non-linear spectrum is stored
                          only in the time/redhsift range where the
                          non-linear corrections were really computed,
                          to avoid dealing with discontinuities in
                          the spline interpolation) */

  double * ln_tau_nl;  /**< list of ln(tau) values
                          ln_tau_nl[index_tau], for the non-linear
                          power spectrum, in growing order. So
                          exp(ln_tau_nl[0]) is the earliest time
                          (i.e. highest redshift), while
                          exp(ln_tau_nl[ln_tau_nl_size-1]) is today
                          (i.e z=0). */

  double * ln_pk_nl;   /**< Non-linear matter power spectrum.
                          depends on indices index_k, index_tau as:
                          ln_pk_nl[index_tau * psp->k_size + index_k] */
  double * ddln_pk_nl; /**< second derivative of above array with respect to log(tau), for spline interpolation. */

  double * ln_pk_cb;           /**< same as ln_pk for baryon+cdm component only */
  double * ddln_pk_cb;         /**< same as ddln_pk for baryon+cdm component only */

  double * ln_pk_cb_l;         /**< same as ln_pk_l for baryon+cdm component only */
  double * ddln_pk_cb_l;       /**< same as ddln_pk_l for baryon+cdm component only */

  double * ln_pk_cb_nl;        /**< same as ln_pk_nl for baryon+cdm component only */
  double * ddln_pk_cb_nl;      /**< same as ddln_pk_nl for baryon+cdm component only */

  int index_tr_delta_g;        /**< index of gamma density transfer function */
  int index_tr_delta_b;        /**< index of baryon density transfer function */
  int index_tr_delta_cdm;      /**< index of cold dark matter density transfer function */
  int index_tr_delta_dcdm;     /**< index of decaying cold dark matter density transfer function */
  int index_tr_delta_scf;      /**< index of scalar field phi transfer function */
  int index_tr_delta_fld;      /**< index of dark energy fluid density transfer function */
  int index_tr_delta_ur;       /**< index of ultra-relativistic neutrinos/relics density transfer function */
  int index_tr_delta_dr;       /**< index of decay radiation density transfer function */
  int index_tr_delta_ncdm1;    /**< index of first species of non-cold dark matter (massive neutrinos, ...) density transfer function */
  int index_tr_delta_tot;      /**< index of total matter density transfer function */
  int index_tr_theta_g;        /**< index of gamma velocity transfer function */
  int index_tr_theta_b;        /**< index of baryon velocity transfer function */
  int index_tr_theta_cdm;      /**< index of cold dark matter velocity transfer function */
  int index_tr_theta_dcdm;     /**< index of decaying cold dark matter velocity transfer function */
  int index_tr_theta_scf;      /**< index of derivative of scalar field phi transfer function */
  int index_tr_theta_fld;      /**< index of dark energy fluid velocity transfer function */
  int index_tr_theta_ur;       /**< index of ultra-relativistic neutrinos/relics velocity transfer function */
  int index_tr_theta_dr;       /**< index of decay radiation velocity transfer function */
  int index_tr_theta_ncdm1;    /**< index of first species of non-cold dark matter (massive neutrinos, ...) velocity transfer function */
  int index_tr_theta_tot;      /**< index of total matter velocity transfer function */
  int index_tr_phi;            /**< index of Bardeen potential phi */
  int index_tr_psi;            /**< index of Bardeen potential psi */
  int index_tr_phi_prime;      /**< index of derivative of Bardeen potential phi */
  int index_tr_h;              /**< index of synchronous gauge metric perturbation h */
  int index_tr_h_prime;        /**< index of synchronous gauge metric perturbation h' */
  int index_tr_eta;            /**< index of synchronous gauge metric perturbation eta */
  int index_tr_eta_prime;      /**< index of synchronous gauge metric perturbation eta' */
  int tr_size;                 /**< total number of species in transfer functions */

  double * matter_transfer;   /**< Matter transfer functions.
                                 Depends on indices index_md,index_tau,index_ic,index_k, index_tr as:
                                 matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + index_tr]
                              */
  double * ddmatter_transfer; /**< second derivative of above array with respect to log(tau), for spline interpolation. */

  /* double * LddCl; /\**< density Cl's in the Limber plus thin shell approximation (then, there are no non-diagonal correlations between various shells of different redshifts); depends on index_tau,index_l as: LddCl[index_tau*psp->psp->l_size[psp->index_md_scalars]+index_l] *\/ */

  /* double * LTdCl; /\**< cross (temperature * density) Cl's in the Limber plus thin shell approximation; depends on index_tau,index_l as: LTdCl[index_tau*psp->psp->l_size[psp->index_md_scalars]+index_l] *\/ */

  //@}

  /** @name - technical parameters */

  //@{

  short spectra_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */

  //@}
};

/*************************************************************************************************************/
/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int spectra_bandpower(
                        struct spectra * psp,
                        int l1,
                        int l2,
                        double * TT_II,
                        double * TT_RI,
                        double * TT_RR
                        );

  int spectra_cl_at_l(
                      struct spectra * psp,
                      double l,
                      double * cl,
                      double ** cl_md,
                      double ** cl_md_ic
                      );

  int spectra_pk_at_z(
                      struct background * pba,
                      struct spectra * psp,
                      enum linear_or_logarithmic mode,
                      double z,
                      double * output_tot,
                      double * output_ic,
                      double * output_cb_tot,
                      double * output_cb_ic
                      );

  int spectra_pk_at_k_and_z(
                            struct background * pba,
                            struct primordial * ppm,
                            struct spectra * psp,
                            double k,
                            double z,
                            double * pk,
                            double * pk_ic,
                            double * pk_cb,
                            double * pk_cb_ic
                            );

  int spectra_pk_nl_at_z(
                         struct background * pba,
                         struct spectra * psp,
                         enum linear_or_logarithmic mode,
                         double z,
                         double * output_tot,
                         double * output_cb_tot
                         );

  int spectra_pk_nl_at_k_and_z(
                               struct background * pba,
                               struct primordial * ppm,
                               struct spectra * psp,
                               double k,
                               double z,
                               double * pk_tot,
                               double * pk_cb_tot
                               );

  int spectra_tk_at_z(
                      struct background * pba,
                      struct spectra * psp,
                      double z,
                      double * output
                      );

  int spectra_tk_at_k_and_z(
                            struct background * pba,
                            struct spectra * psp,
                            double k,
                            double z,
                            double * output
                            );

  int spectra_init(
                   struct precision * ppr,
                   struct background * pba,
                   struct perturbs * ppt,
                   struct primordial * ppm,
                   struct nonlinear *pnl,
                   struct transfers * ptr,
                   struct spectra * psp
                   );

  int spectra_free(
                   struct spectra * psp
                   );

  int spectra_indices(
                      struct background * pba,
                      struct perturbs * ppt,
                      struct transfers * ptr,
                      struct primordial * ppm,
                      struct spectra * psp
                      );

  int spectra_cls(
                  struct background * pba,
                  struct perturbs * ppt,
                  struct transfers * ptr,
                  struct primordial * ppm,
                  struct spectra * psp
                  );

  int spectra_compute_cl(
                         struct background * pba,
                         struct perturbs * ppt,
                         struct transfers * ptr,
                         struct primordial * ppm,
                         struct spectra * psp,
                         int index_md,
                         int index_ic1,
                         int index_ic2,
                         int index_l,
                         int cl_integrand_num_columns,
                         double * cl_integrand,
                         double * primordial_pk,
                         double * transfer_ic1,
                         double * transfer_ic2
                         );

  int spectra_k_and_tau(
                        struct background * pba,
                        struct perturbs * ppt,
                        struct nonlinear *pnl,
                        struct spectra * psp
                        );

  int spectra_pk(
                 struct background * pba,
                 struct perturbs * ppt,
                 struct primordial * ppm,
                 struct nonlinear *pnl,
                 struct spectra * psp
                 );

  int spectra_sigma(
                    struct background * pba,
                    struct primordial * ppm,
                    struct spectra * psp,
                    double R,
                    double z,
                    double *sigma
                    );

  int spectra_sigma_cb(
                    struct background * pba,
                    struct primordial * ppm,
                    struct spectra * psp,
                    double R,
                    double z,
                    double *sigma_cb
                    );

  int spectra_matter_transfers(
                               struct background * pba,
                               struct perturbs * ppt,
                               struct spectra * psp
                               );

  int spectra_output_tk_titles(struct background *pba,
                               struct perturbs *ppt,
                               enum file_format output_format,
                               char titles[_MAXTITLESTRINGLENGTH_]
                               );

  int spectra_output_tk_data(
                             struct background * pba,
                             struct perturbs * ppt,
                             struct spectra * psp,
                             enum file_format output_format,
                             double z,
                             int number_of_titles,
                             double *data
                             );

  int spectra_firstline_and_ic_suffix(struct perturbs *ppt,
                                     int index_ic,
                                     char first_line[_LINE_LENGTH_MAX_],
                                     FileName ic_suffix);

  int spectra_fast_pk_at_kvec_and_zvec(
				       struct background * pba,
				       struct spectra * psp,
				       double * kvec,
				       int kvec_size,
				       double * zvec,
				       int zvec_size,
				       double * pk_tot_out, /* (must be already allocated with kvec_size*zvec_size) */
                       double * pk_cb_tot_out,
				       int nonlinear);

#ifdef __cplusplus
}
#endif

#endif
/* @endcond */
