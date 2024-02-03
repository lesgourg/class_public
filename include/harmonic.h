/** @file harmonic.h Documented includes for harmonic module */

#ifndef __HARMONIC__
#define __HARMONIC__

#include "transfer.h"

/**
 * Structure containing everything about anisotropy and Fourier power spectra that other modules need to know.
 *
 * Once initialized by harmonic_init(), contains a table of all
 * \f$ C_l\f$'s and P(k) as a function of multipole/wavenumber,
 * mode (scalar/tensor...), type (for \f$ C_l\f$'s: TT, TE...),
 * and pairs of initial conditions (adiabatic, isocurvatures...).
 */

struct harmonic {

  /** @name - input parameters initialized by user in input module
      (all other quantities are computed in this module, given these parameters
      and the content of the 'background', 'perturbs', 'transfers' and
      'primordial' structures) */

  //@{

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

  double ** cl;   /**< table of anisotropy spectra for each mode, multipole, pair of initial conditions and types, cl[index_md][(index_l * phr->ic_ic_size[index_md] + index_ic1_ic2) * phr->ct_size + index_ct] */
  double ** ddcl; /**< second derivatives of previous table with respect to l, in view of spline interpolation */

  //@}

  /** @name - technical parameters */

  //@{

  struct fourier * pfo; /**< a pointer to the fourier structure is
                           stored in the harmonic structure. This odd,
                           unusual and unelegant feature has been
                           introduced in v2.8 in order to keep in use
                           some deprecated functions harmonic_pk_...()
                           that are now pointing at new function
                           fourier_pk_...(). In the future, if the
                           deprecated functions are removed, it will
                           be possible to remove also this pointer. */

  short harmonic_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

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

  /* external functions (meant to be called from other modules) */

  int harmonic_cl_at_l(
                       struct harmonic * phr,
                       double l,
                       double * cl,
                       double ** cl_md,
                       double ** cl_md_ic
                       );

  /* internal functions */

  int harmonic_init(
                    struct precision * ppr,
                    struct background * pba,
                    struct perturbations * ppt,
                    struct primordial * ppm,
                    struct fourier *pfo,
                    struct transfer * ptr,
                    struct harmonic * phr
                    );

  int harmonic_free(
                    struct harmonic * phr
                    );

  int harmonic_indices(
                       struct background * pba,
                       struct perturbations * ppt,
                       struct transfer * ptr,
                       struct primordial * ppm,
                       struct harmonic * phr
                       );

  int harmonic_cls(
                   struct precision * ppr,
                   struct background * pba,
                   struct perturbations * ppt,
                   struct transfer * ptr,
                   struct primordial * ppm,
                   struct harmonic * phr
                   );

  int harmonic_compute_cl(
                          struct precision * ppr,
                          struct background * pba,
                          struct perturbations * ppt,
                          struct transfer * ptr,
                          struct primordial * ppm,
                          struct harmonic * phr,
                          int index_md,
                          int index_ic1,
                          int index_ic2,
                          int index_l,
                          int cl_integrand_num_columns,
                          double * cl_integrand,
                          double * cl_integrand_limber,
                          double * primordial_pk,
                          double * transfer_ic1,
                          double * transfer_ic2
                          );

  int harmonic_k_and_tau(
                         struct background * pba,
                         struct perturbations * ppt,
                         struct fourier *pfo,
                         struct harmonic * phr
                         );

  /* deprecated functions (since v2.8) */

  int harmonic_pk_at_z(
                       struct background * pba,
                       struct harmonic * phr,
                       enum linear_or_logarithmic mode,
                       double z,
                       double * output_tot,
                       double * output_ic,
                       double * output_cb_tot,
                       double * output_cb_ic
                       );

  int harmonic_pk_at_k_and_z(
                             struct background * pba,
                             struct primordial * ppm,
                             struct harmonic * phr,
                             double k,
                             double z,
                             double * pk,
                             double * pk_ic,
                             double * pk_cb,
                             double * pk_cb_ic
                             );

  int harmonic_pk_nl_at_z(
                          struct background * pba,
                          struct harmonic * phr,
                          enum linear_or_logarithmic mode,
                          double z,
                          double * output_tot,
                          double * output_cb_tot
                          );

  int harmonic_pk_nl_at_k_and_z(
                                struct background * pba,
                                struct primordial * ppm,
                                struct harmonic * phr,
                                double k,
                                double z,
                                double * pk_tot,
                                double * pk_cb_tot
                                );

  int harmonic_fast_pk_at_kvec_and_zvec(
                                        struct background * pba,
                                        struct harmonic * phr,
                                        double * kvec,
                                        int kvec_size,
                                        double * zvec,
                                        int zvec_size,
                                        double * pk_tot_out,
                                        double * pk_cb_tot_out,
                                        int nonlinear);

  int harmonic_sigma(
                     struct background * pba,
                     struct primordial * ppm,
                     struct harmonic * phr,
                     double R,
                     double z,
                     double *sigma
                     );

  int harmonic_sigma_cb(
                        struct background * pba,
                        struct primordial * ppm,
                        struct harmonic * phr,
                        double R,
                        double z,
                        double *sigma_cb
                        );

  /* deprecated functions (since v2.1) */

  int harmonic_tk_at_z(
                       struct background * pba,
                       struct harmonic * phr,
                       double z,
                       double * output
                       );

  int harmonic_tk_at_k_and_z(
                             struct background * pba,
                             struct harmonic * phr,
                             double k,
                             double z,
                             double * output
                             );

  /* end deprecated functions */

#ifdef __cplusplus
}
#endif

#endif
/* @endcond */
