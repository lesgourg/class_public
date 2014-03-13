/** @file spectra.h Documented includes for spectra module */

#ifndef __SPECTRA__
#define __SPECTRA__

#include "transfer.h"

/**
 * Structure containing everything about anisotropy and Fourier power spectra that other modules need to know.
 *
 * Once initialized by spectra_init(), contains a table of all
 * C_l's and P(k) as a function of multipole/wavenumber,
 * mode (scalar/tensor...), type (for C_l's: TT, TE...),
 * and pairs of initial conditions (adiabatic, isocurvatures...).
 */

struct spectra {

  /** @name - input parameters initialized by user in input module
      (all other quantitites are computed in this module, given these parameters
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

  int has_tt; /**< do we want C_l^TT ? (T = temperature) */
  int has_ee; /**< do we want C_l^EE ? (E = E-polarization) */
  int has_te; /**< do we want C_l^TE ? */
  int has_bb; /**< do we want C_l^BB ? (B = B-polarization) */
  int has_pp; /**< do we want C_l^phi-phi ? (phi = CMB lensing potential) */
  int has_tp; /**< do we want C_l^T-phi ? */
  int has_ep; /**< do we want C_l^E-phi ? */
  int has_dd; /**< do we want C_l^dd ? (d = density) */
  int has_td; /**< do we want C_l^Td ? */
  int has_pd; /**< do we want C_l^phi-d ? */
  int has_ll; /**< do we want C_l^l-l ? (l = galaxy lensing potential) */
  int has_tl; /**< do we want C_l^T-l ? */
  int has_dl; /**< do we want C_l^d-l ? */

  int index_ct_tt; /**< index for type C_l^TT */
  int index_ct_ee; /**< index for type C_l^EE */
  int index_ct_te; /**< index for type C_l^TE */
  int index_ct_bb; /**< index for type C_l^BB */
  int index_ct_pp; /**< index for type C_l^phi-phi */
  int index_ct_tp; /**< index for type C_l^T-phi */
  int index_ct_ep; /**< index for type C_l^E-phi */
  int index_ct_dd; /**< first index for type C_l^dd ((d_size*d_size-(d_size-non_diag)*(d_size-non_diag-1)/2) values) */
  int index_ct_td; /**< first index for type C_l^Td (d_size values) */
  int index_ct_pd; /**< first index for type C_l^pd (d_size values) */
  int index_ct_ll; /**< first index for type C_l^ll ((d_size*d_size-(d_size-non_diag)*(d_size-non_diag-1)/2) values) */
  int index_ct_tl; /**< first index for type C_l^Tl (d_size values) */
  int index_ct_dl; /**< first index for type C_l^dl (d_size values) */

  int d_size;

  int ct_size; /**< number of C_l types requested */

  //@}

  /** @name - table of pre-computed C_l values, and related quantitites */

  //@{

  int * l_size;   /**< number of multipole values for each requested mode, l_size[index_md] */

  int l_size_max; /**< greatest of all l_size[index_md] */

  double * l;    /**< list of multipole values l[index_l] */


  int ** l_max_ct;    /**< last multipole (given as an input) at which
                         we want to output C_ls for a given mode and type;
                         l[index_md][l_size[index_md]-1] can be larger
                         than l_max[index_md], in order to ensure a
                         better interpolation with no boundary effects */

  int * l_max;    /**< last multipole (given as an input) at which
                     we want to output C_ls for a given mode (maximized over types);
                     l[index_md][l_size[index_md]-1] can be larger
                     than l_max[index_md], in order to ensure a
                     better interpolation with no boundary effects */

  int l_max_tot; /**< last multipole (given as an input) at which
                    we want to output C_ls (maximized over modes and types);
                    l[index_md][l_size[index_md]-1] can be larger
                    than l_max[index_md], in order to ensure a
                    better interpolation with no boundary effects */

  double ** cl;   /**< table of anisotropy spectra for each mode, multipole, pair of initial conditions and types, cl[index_md][(index_l * psp->ic_ic_size[index_md] + index_ic1_ic2) * psp->ct_size + index_ct] */
  double ** ddcl; /**< second derivatives of previous table with respect to l, in view of spline interpolation */

  double alpha_II_2_20;
  double alpha_RI_2_20;
  double alpha_RR_2_20;

  double alpha_II_21_200;
  double alpha_RI_21_200;
  double alpha_RR_21_200;

  double alpha_II_201_2500;
  double alpha_RI_201_2500;
  double alpha_RR_201_2500;

  double alpha_II_2_2500;
  double alpha_RI_2_2500;
  double alpha_RR_2_2500;

  double alpha_kp;
  double alpha_k1;
  double alpha_k2;

  //@}

  /** @name - table of pre-computed matter power spectrum P(k) values, and related quantitites */

  //@{

  int ln_k_size;    /**< number ln(k) values */
  double * ln_k;    /**< list of ln(k) values ln_k[index_k] */

  int ln_tau_size;  /**< number ln(tau) values (only one if z_max_pk = 0) */
  double * ln_tau;  /**< list of ln(tau) values ln_tau[index_tau] */

  double * ln_pk;   /**< Matter power spectrum.
                       depends on indices index_md, index_ic1, index_ic2, index_k, index_tau as:
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

  double * ln_pk_nl;   /**< Non-linear matter power spectrum.
                          depends on indices index_k, index_tau as:
                          ln_pk_nl[index_tau * psp->k_size + index_k]
                    */
  double * ddln_pk_nl; /**< second derivative of above array with respect to log(tau), for spline interpolation. */

  int index_tr_delta_g;        /**< index of gamma density transfer function */
  int index_tr_delta_b;        /**< index of baryon density transfer function */
  int index_tr_delta_cdm;      /**< index of cold dark matter density transfer function */
  int index_tr_delta_fld;      /**< index of dark energy fluid density transfer function */
  int index_tr_delta_ur;       /**< index of ultra-relativistic neutrinos/relics density transfer function */
  int index_tr_delta_ncdm1;    /**< index of first species of non-cold dark matter (massive neutrinos, ...) density transfer function */
  int index_tr_delta_tot;      /**< index of total matter density transfer function */
  int index_tr_theta_g;        /**< index of gamma velocity transfer function */
  int index_tr_theta_b;        /**< index of baryon velocity transfer function */
  int index_tr_theta_cdm;      /**< index of cold dark matter velocity transfer function */
  int index_tr_theta_fld;      /**< index of dark energy fluid velocity transfer function */
  int index_tr_theta_ur;       /**< index of ultra-relativistic neutrinos/relics velocity transfer function */
  int index_tr_theta_ncdm1;    /**< index of first species of non-cold dark matter (massive neutrinos, ...) velocity transfer function */
  int index_tr_theta_tot;      /**< index of total matter velocity transfer function */
  int tr_size;                 /**< total number of species in transfer functions */

  double * matter_transfer;   /**< Matter transfer functions.
                                 Depends on indices index_md,index_tau,index_ic,index_k, index_tr as:
                                 matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + index_tr]
                              */
  double * ddmatter_transfer; /**< second derivative of above array with respect to log(tau), for spline interpolation. */

  /* double * LddCl; /\**< density Cl's in the Limber plus thin shell approximation (then, there are no non-diagonal correlations betzeen various shells of different redshifts); depends on index_tau,index_l as: LddCl[index_tau*psp->psp->l_size[psp->index_md_scalars]+index_l] *\/ */

  /* double * LTdCl; /\**< cross (temperature * density) Cl's in the Limber plus thin shell approximation; depends on index_tau,index_l as: LTdCl[index_tau*psp->psp->l_size[psp->index_md_scalars]+index_l] *\/ */

  //@}

  /** @name - technical parameters */

  //@{

  short spectra_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */

  //@}
};

/*************************************************************************************************************/

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
                      double * output_ic
                      );

  int spectra_pk_at_k_and_z(
                            struct background * pba,
                            struct primordial * ppm,
                            struct spectra * psp,
                            double k,
                            double z,
                            double * pk,
                            double * pk_ic
                            );

  int spectra_pk_nl_at_z(
                         struct background * pba,
                         struct spectra * psp,
                         enum linear_or_logarithmic mode,
                         double z,
                         double * output_tot
                         );

  int spectra_pk_nl_at_k_and_z(
                               struct background * pba,
                               struct primordial * ppm,
                               struct spectra * psp,
                               double k,
                               double z,
                               double * pk_tot
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

  int spectra_matter_transfers(
                               struct background * pba,
                               struct perturbs * ppt,
                               struct spectra * psp
                               );

#ifdef __cplusplus
}
#endif

#endif
