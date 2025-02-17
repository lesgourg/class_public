/** @file hmcode.h Documented includes for HMcode module */

#include "primordial.h"

#ifndef __HMCODE__
#define __HMCODE__

struct hmcode_growth {
  double * a_table;
  int a_size;
  double * growth_table;
  double * normgrowth_table;
  double * ddnormgrowth_table;
  double * dda_table;

  int gt_size;
  int index_gt_g;
  int index_gt_ddg;
  int index_gt_intg;
  int index_gt_ddintg;
  int index_gt_Omnorad;
  int index_gt_ddOmnorad;

  double Omega_cdm;
  double Omega_b;
  double Omega_v;
  double Omega_w;
  double w0;
  double wa;
  double Omega_m;
  double Omega_nu;
  double a_nu;
  double Tcmb;

  double om_m;
  double fnu;

  double a_ini;
  double a_final;

  double smallest_allowed_variation;

  int index_fg_g;
  int index_fg_dg;
  int fg_size;

  int index_norad_H2;
  int index_norad_AH;
  int index_norad_Om;
  int norad_size;
  double * pvecnorad;

  double * P_nw;
  double lnk_norm;

  ErrorMsg error_message; 	/**< zone for writing error messages */

};

/**
 * Structure containing variables used only internally in fourier module by various functions.
 *
 */

struct hmcode_workspace {

  /** @name - quantitites used by HMcode */

  //@{

  double * rtab; /** List of R values */
  double * stab; /** List of Sigma Values */
  double * ddstab; /** Splined sigma */

  double * growtable;
  double * ztable;
  double * tautable;

  double ** sigma_8;
  double ** sigma_disp;
  double ** sigma_disp_100;
  double ** sigma_prime;

  double dark_energy_correction; /** this is the ratio [g_wcdm(z_infinity)/g_lcdm(z_infinity)]^1.5
                                  * (power comes from Dolag et al. (2004) correction)
                                  * it is 1, if has_fld == _FALSE_ */
  double* lnk_wiggle;
  double* pk_wiggle;
  double* ddpk_wiggle;

  struct hmcode_growth* phg;
  //@}

};

/********************************************************************************/

/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int hmcode(
             struct precision *ppr,
             struct background *pba,
             struct perturbations *ppt,
             struct primordial *ppm,
             struct fourier *pfo,
             int index_pk,
             int index_tau,
             double tau,
             double *pk_nl,
             double **lnpk_l,
             double **ddlnpk_l,
             double *k_nl,
             short *halofit_found_k_max,
             struct hmcode_workspace *phw
             );

  int hmcode_compute(
                     struct precision *ppr,
                     struct background *pba,
                     struct perturbations *ppt,
                     struct primordial *ppm,
                     struct fourier *pfo,
                     int index_pk,
                     int index_tau,
                     double tau,
                     double *pk_nl,
                     double **lnpk_l,
                     double **ddlnpk_l,
                     double *k_nl,
                     short *halofit_found_k_max,
                     struct hmcode_workspace *phw
                     );

  int hmcode_workspace_init(
                            struct precision *ppr,
                            struct background *pba,
                            struct fourier *pfo,
                            struct hmcode_workspace *phw
                            );

  int hmcode_workspace_free(
                            struct fourier *pfo,
                            struct hmcode_workspace *phw
                            );

  int hmcode_dark_energy_correction(
                                    struct precision *ppr,
                                    struct background *pba,
                                    struct fourier *pfo,
                                    struct hmcode_workspace *phw
                                    );

  int hmcode_baryonic_feedback(
                               struct fourier *pfo
                               );


  int hmcode_fill_sigtab(
                         struct precision *ppr,
                         struct background *pba,
                         struct perturbations *ppt,
                         struct primordial *ppm,
                         struct fourier *pfo,
                         int index_tau,
                         int index_pk,
                         double **lnpk_l,
                         double **ddlnpk_l,
                         struct hmcode_workspace *phw
                         );

  int hmcode_fill_growtab(
                          struct precision *ppr,
                          struct background *pba,
                          struct fourier *pfo,
                          struct hmcode_workspace *phw
                          );

  int hmcode_growint(
                     struct precision *ppr,
                     struct background *pba,
                     struct fourier *pfo,
                     double a,
                     double w,
                     double wa,
                     double *growth
                     );

  int hmcode_window_nfw(
                        struct fourier *pfo,
                        double k,
                        double rv,
                        double c,
                        double *window_nfw
                        );

  int hmcode_halomassfunction(
                              double nu,
                              double *hmf
                              );

  int hmcode_sigma8_at_z(
                         struct background *pba,
                         struct fourier *pfo,
                         double z,
                         double *sigma_8,
                         double *sigma_8_cb,
                         struct hmcode_workspace *phw
                         );

  int hmcode_sigmadisp_at_z(
                            struct background *pba,
                            struct fourier * pfo,
                            double z,
                            double * sigma_disp,
                            double * sigma_disp_cb,
                            struct hmcode_workspace * phw
                            );

  int hmcode_sigmadisp100_at_z(
                               struct background *pba,
                               struct fourier *pfo,
                               double z,
                               double *sigma_disp_100,
                               double *sigma_disp_100_cb,
                               struct hmcode_workspace *phw
                               );

  int hmcode_sigmaprime_at_z(
                             struct background *pba,
                             struct fourier *pfo,
                             double z,
                             double *sigma_prime,
                             double *sigma_prime_cb,
                             struct hmcode_workspace *phw
                             );

  int hmcode_nowiggle_init(
                           struct fourier *pfo,
                           double **lnpk_l,
                           double **ddlnpk_l,
                           int index_pk,
                           struct hmcode_workspace *phw);

  int hmcode_wnw_split(
                       struct precision *ppr,
                       struct background *pba,
                       struct primordial * ppm,
                       struct fourier *pfo
                       );

  int hmcode_nowiggle(
                      struct fourier *pfo,
                      double * lnpk_l,
                      double * ddlnpk_l,
                      int nw_size,
                      double * ln_k_nw,
                      double * pk_nw,
                      double * pk_w
                      );

  int hmcode_eisenstein_hu(
                           struct precision *ppr,
                           struct background *pba,
                           struct primordial * ppm,
                           struct fourier *pfo,
                           double * ln_k,
                           int k_size,
                           double * ln_pk_nowiggle
                           );

  int hmcode_noradiation_growth_init(
                                     struct precision *ppr,
                                     struct background *pba,
                                     struct fourier *pfo,
                                     struct hmcode_workspace *phw
                                     );

  int hmcode_noradiation_growth_free(
                                     struct fourier *pfo,
                                     struct hmcode_workspace *phw
                                     );

  int hmcode_noradiation_growth_compute(
                                        struct fourier *pfo,
                                        struct hmcode_workspace *phw
                                        );

  int hmcode_growth_derivs(
                           double a,
                           double *y,
                           double *dy,
                           void *parameters_and_workspace,
                           ErrorMsg error_message
                           );

  int hmcode_growth_sources(
                            double a,
                            double *y,
                            double *dy,
                            int index_a,
                            void *parameters_and_workspace,
                            ErrorMsg error_message
                            );

  int hmcode_norad(
                   struct hmcode_growth *phg,
                   double a
                   );

  int hmcode_cbnu_ratio(
                        double k,
                        double z,
                        double fnu,
                        double omega_m,
                        double Tcmb,
                        double growth,
                        double *cbnu_ratio
                        );

#ifdef __cplusplus
}
#endif

#endif
/* @endcond */
