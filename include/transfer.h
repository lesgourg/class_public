/** @file transfer.h Documented includes for transfer module. */

#ifndef __TRANSFER__
#define __TRANSFER__

#include "fourier.h"
#include "hyperspherical.h"
#include <sys/shm.h>
#include <sys/stat.h>
#include "errno.h"

/* macro: test if index_tt is in the range between index and index+num, while the flag is true */
#define _index_tt_in_range_(index,num,flag) (flag == _TRUE_) && (index_tt >= index) && (index_tt < index+num)
/* macro: test if index_tt corresponds to an integrated nCl/sCl contribution */
#define _integrated_ncl_ (_index_tt_in_range_(ptr->index_tt_lensing, ppt->selection_num, ppt->has_cl_lensing_potential)) || \
  (_index_tt_in_range_(ptr->index_tt_nc_lens, ppt->selection_num, ppt->has_nc_lens)) || \
    (_index_tt_in_range_(ptr->index_tt_nc_g4,   ppt->selection_num, ppt->has_nc_gr)) || \
    (_index_tt_in_range_(ptr->index_tt_nc_g5,   ppt->selection_num, ppt->has_nc_gr))
/* macro: test if index_tt corresponds to an non-integrated nCl/sCl contribution */
#define _nonintegrated_ncl_ (_index_tt_in_range_(ptr->index_tt_density, ppt->selection_num, ppt->has_nc_density)) || \
  (_index_tt_in_range_(ptr->index_tt_rsd,     ppt->selection_num, ppt->has_nc_rsd)) || \
    (_index_tt_in_range_(ptr->index_tt_d0,      ppt->selection_num, ppt->has_nc_rsd)) || \
    (_index_tt_in_range_(ptr->index_tt_d1,      ppt->selection_num, ppt->has_nc_rsd)) || \
    (_index_tt_in_range_(ptr->index_tt_nc_g1,   ppt->selection_num, ppt->has_nc_gr))  || \
    (_index_tt_in_range_(ptr->index_tt_nc_g2,   ppt->selection_num, ppt->has_nc_gr))  || \
    (_index_tt_in_range_(ptr->index_tt_nc_g3,   ppt->selection_num, ppt->has_nc_gr))
/* macro: bin number associated to particular redshift bin and selection function for non-integrated contributions*/
#define _get_bin_nonintegrated_ncl_(index_tt)                           \
  if (_index_tt_in_range_(ptr->index_tt_density, ppt->selection_num, ppt->has_nc_density)) \
    bin = index_tt - ptr->index_tt_density;                             \
  if (_index_tt_in_range_(ptr->index_tt_rsd,     ppt->selection_num, ppt->has_nc_rsd)) \
    bin = index_tt - ptr->index_tt_rsd;                                 \
  if (_index_tt_in_range_(ptr->index_tt_d0,      ppt->selection_num, ppt->has_nc_rsd)) \
    bin = index_tt - ptr->index_tt_d0;                                  \
  if (_index_tt_in_range_(ptr->index_tt_d1,      ppt->selection_num, ppt->has_nc_rsd)) \
    bin = index_tt - ptr->index_tt_d1;                                  \
  if (_index_tt_in_range_(ptr->index_tt_nc_g1,   ppt->selection_num, ppt->has_nc_gr)) \
    bin = index_tt - ptr->index_tt_nc_g1;                               \
  if (_index_tt_in_range_(ptr->index_tt_nc_g2,   ppt->selection_num, ppt->has_nc_gr)) \
    bin = index_tt - ptr->index_tt_nc_g2;                               \
  if (_index_tt_in_range_(ptr->index_tt_nc_g3,   ppt->selection_num, ppt->has_nc_gr)) \
    bin = index_tt - ptr->index_tt_nc_g3;
/* macro: bin number associated to particular redshift bin and selection function for integrated contributions*/
#define _get_bin_integrated_ncl_(index_tt)                              \
  if (_index_tt_in_range_(ptr->index_tt_lensing, ppt->selection_num, ppt->has_cl_lensing_potential)) \
    bin = index_tt - ptr->index_tt_lensing;                             \
  if (_index_tt_in_range_(ptr->index_tt_nc_lens, ppt->selection_num, ppt->has_nc_lens)) \
    bin = index_tt - ptr->index_tt_nc_lens;                             \
  if (_index_tt_in_range_(ptr->index_tt_nc_g4,   ppt->selection_num, ppt->has_nc_gr)) \
    bin = index_tt - ptr->index_tt_nc_g4;                               \
  if (_index_tt_in_range_(ptr->index_tt_nc_g5,   ppt->selection_num, ppt->has_nc_gr)) \
    bin = index_tt - ptr->index_tt_nc_g5;
/**
 * Structure containing everything about transfer functions in
 * harmonic space \f$ \Delta_l^{X} (q) \f$ that other modules need to
 * know.
 *
 * Once initialized by transfer_init(), contains all tables of
 * transfer functions used for interpolation in other modules, for all
 * requested modes (scalar/vector/tensor), initial conditions, types
 * (temperature, polarization, etc), multipoles l, and wavenumbers q.
 *
 * Wavenumbers are called q in this module and k in the perturbation
 * module. In flat universes k=q. In non-flat universes q and k differ
 * through q2 = k2 + K(1+m), where m=0,1,2 for scalar, vector,
 * tensor. q should be used throughout the transfer module, except
 * when interpolating or manipulating the source functions S(k,tau)
 * calculated in the perturbation module: for a given value of q, this
 * should be done at the corresponding k(q).
 *
 * The content of this structure is entirely computed in this module,
 * given the content of the 'precision', 'bessels', 'background',
 * 'thermodynamics' and 'perturbation' structures.
 */

struct transfer {

  /** @name - input parameters initialized by user in input module
   *  (all other quantities are computed in this module, given these
   *  parameters and the content of previous structures) */

  //@{

  double lcmb_rescale; /**< normally set to one, can be used
                          exceptionally to rescale by hand the CMB
                          lensing potential */
  double lcmb_tilt;    /**< normally set to zero, can be used
                          exceptionally to tilt by hand the CMB
                          lensing potential */
  double lcmb_pivot;   /**< if lcmb_tilt non-zero, corresponding pivot
                          scale */

  double selection_bias[_SELECTION_NUM_MAX_];               /**< light-to-mass bias in the transfer function of density number count */
  double selection_magnification_bias[_SELECTION_NUM_MAX_]; /**< magnification bias in the transfer function of density number count */

  short has_nz_file;     /**< Has dN/dz (selection function) input file? */
  short has_nz_analytic; /**< Use analytic form for dN/dz (selection function) distribution? */
  FileName nz_file_name; /**< dN/dz (selection function) input file name */
  int nz_size;           /**< number of redshift values in input tabulated selection function */
  double * nz_z;         /**< redshift values in input tabulated selection function */
  double * nz_nz;        /**< input tabulated values of selection function */
  double * nz_ddnz;      /**< second derivatives in splined selection function*/

  short has_nz_evo_file;      /**< Has dN/dz (evolution function) input file? */
  short has_nz_evo_analytic;  /**< Use analytic form for dN/dz (evolution function) distribution? */
  FileName nz_evo_file_name;  /**< dN/dz (evolution function) input file name */
  int nz_evo_size;            /**< number of redshift values in input tabulated evolution function */
  double * nz_evo_z;          /**< redshift values in input tabulated evolution function */
  double * nz_evo_nz;         /**< input tabulated values of evolution function */
  double * nz_evo_dlog_nz;    /**< log of tabulated values of evolution function */
  double * nz_evo_dd_dlog_nz; /**< second derivatives in splined log of evolution function */

  //@}

  /** @name - flag stating whether we need transfer functions at all */

  //@{

  short has_cls; /**< copy of same flag in perturbation structure */

  //@}

  /** @name - number of modes and transfer function types */

  //@{

  int md_size;       /**< number of modes included in computation */

  int index_tt_t0;      /**< index for transfer type = temperature (j=0 term) */
  int index_tt_t1;      /**< index for transfer type = temperature (j=1 term) */
  int index_tt_t2;      /**< index for transfer type = temperature (j=2 term) */
  int index_tt_e;       /**< index for transfer type = E-polarization */
  int index_tt_b;       /**< index for transfer type = B-polarization */
  int index_tt_lcmb;    /**< index for transfer type = CMB lensing */
  int index_tt_density; /**< index for first bin of transfer type = matter density */
  int index_tt_lensing; /**< index for first bin of transfer type = galaxy lensing */

  int index_tt_rsd;     /**< index for first bin of transfer type = redshift space distortion of number count */
  int index_tt_d0;      /**< index for first bin of transfer type = doppler effect for of number count (j=0 term) */
  int index_tt_d1;      /**< index for first bin of transfer type = doppler effect for of number count (j=1 term) */
  int index_tt_nc_lens; /**< index for first bin of transfer type = lensing for of number count */
  int index_tt_nc_g1;   /**< index for first bin of transfer type = gravity term G1 for of number count */
  int index_tt_nc_g2;   /**< index for first bin of transfer type = gravity term G2 for of number count */
  int index_tt_nc_g3;   /**< index for first bin of transfer type = gravity term G3 for of number count */
  int index_tt_nc_g4;   /**< index for first bin of transfer type = gravity term G3 for of number count */
  int index_tt_nc_g5;   /**< index for first bin of transfer type = gravity term G3 for of number count */

  int * tt_size;     /**< number of requested transfer types tt_size[index_md] for each mode */

  //@}

  /** @name - number and list of multipoles */

  //@{

  int ** l_size_tt;  /**< number of multipole values for which we effectively compute the transfer function,l_size_tt[index_md][index_tt] */

  int * l_size;   /**< number of multipole values for each requested mode, l_size[index_md] */

  int l_size_max; /**< greatest of all l_size[index_md] */

  int * l;        /**< list of multipole values l[index_l] */

  //int * l_size_bessel; /**< for each wavenumber, maximum value of l at which bessel functions must be evaluated */

  double angular_rescaling; /**< correction between l and k space due to curvature (= comoving angular diameter distance to recombination / comoving radius to recombination) */

  //@}

  /** @name - number and list of wavenumbers */

  //@{

  size_t q_size; /**< number of wavenumber values corresponding to k up to k_max_cl */

  double * q;  /**< list of wavenumber values, q[index_q] */

  double ** k; /**< list of wavenumber values for each requested mode, k[index_md][index_q]. In flat universes k=q. In non-flat universes q and k differ through q2 = k2 + K(1+m), where m=0,1,2 for scalar, vector, tensor. q should be used throughout the transfer module, excepted when interpolating or manipulating the source functions S(k,tau): for a given value of q this should be done in k(q). */

  int index_q_flat_approximation; /**< index of the first q value using the flat rescaling approximation */

  short do_lcmb_full_limber; /**< in this particular run, will we use the full Limber scheme? */

  size_t q_size_limber; /**< number of wavenumber values corresponding to k up to k_max */

  double * q_limber;  /**< list of wavenumber values used in full limber scheme, q_limber[index_q] */

  double ** k_limber; /**< list of wavenumber values used in full limber scheme */

  //@}

  /** @name - transfer functions */

  //@{

  double ** transfer; /**< table of transfer functions for each mode, initial condition, type, multipole and wavenumber, with argument transfer[index_md][((index_ic * ptr->tt_size[index_md] + index_tt) * ptr->l_size[index_md] + index_l) * ptr->q_size + index_q] */

  double ** transfer_limber; /**< table of transfer functions used in full limber scheme */

  //@}

  /** @name - technical parameters */

  //@{

  short transfer_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */

  //@}
};

/**
 * Structure containing all the quantities that each thread needs to
 * know for computing transfer functions (but that can be forgotten
 * once the transfer functions are known, otherwise they would be
 * stored in the transfer module)
 */

struct transfer_workspace {

  /** @name - quantities related to Bessel functions */

  //@{

  HyperInterpStruct HIS; /**< structure containing all hyperspherical bessel functions (flat case) or all hyperspherical bessel functions for a given value of beta=q/sqrt(|K|) (non-flat case). HIS = Hyperspherical Interpolation Structure. */

  int HIS_allocated; /**< flag specifying whether the previous structure has been allocated */

  HyperInterpStruct * pBIS;  /**< pointer to structure containing all the spherical bessel functions of the flat case (used even in the non-flat case, for approximation schemes). pBIS = pointer to Bessel Interpolation Structure. */

  int l_size;        /**< number of l values */

  //@}

  /** @name - quantities related to the integrand of the transfer functions (most of them are arrays of time) */

  //@{

  int tau_size;                  /**< number of discrete time values for a given type */
  int tau_size_max;              /**< maximum number of discrete time values for all types */
  double * interpolated_sources; /**< interpolated_sources[index_tau]:
                                    sources interpolated from the
                                    perturbation module at the right
                                    value of k */
  double * sources;              /**< sources[index_tau]: sources
                                    used in transfer module, possibly
                                    differing from those in the
                                    perturbation module by some
                                    resampling or rescaling */
  double * tau0_minus_tau;       /**< tau0_minus_tau[index_tau]: values of (tau0 - tau) */
  double * w_trapz;              /**< w_trapz[index_tau]: values of weights in trapezoidal integration (related to time steps) */
  double * chi;                  /**< chi[index_tau]: value of argument of bessel
                                    function: k(tau0-tau) (flat case)
                                    or sqrt(|K|)(tau0-tau) (non-flat
                                    case) */
  double * cscKgen;              /**< cscKgen[index_tau]: useful trigonometric function */
  double * cotKgen;              /**< cotKgen[index_tau]: useful trigonometric function */

  //@}

  /** @name - parameters defining the spatial curvature (copied from background structure) */

  //@{

  double K; /**< curvature parameter (see background module for details) */
  int sgnK; /**< 0 (flat), 1 (positive curvature, spherical, closed), -1 (negative curvature, hyperbolic, open) */

  //@}

  double tau0_minus_tau_cut; /**< critical value of (tau0-tau) in time cut approximation for the wavenumber at hand */
  short neglect_late_source; /**< flag stating whether we use the time cut approximation for the wavenumber at hand */
};

/**
 * enumeration of possible source types. This looks redundant with
 * respect to the definition of indices index_tt_... This definition is however
 * convenient and time-saving: it allows to use a "case" statement in
 * transfer_radial_function()
 */

typedef enum {SCALAR_TEMPERATURE_0,
              SCALAR_TEMPERATURE_1,
              SCALAR_TEMPERATURE_2,
              SCALAR_POLARISATION_E,
              VECTOR_TEMPERATURE_1,
              VECTOR_TEMPERATURE_2,
              VECTOR_POLARISATION_E,
              VECTOR_POLARISATION_B,
              TENSOR_TEMPERATURE_2,
              TENSOR_POLARISATION_E,
              TENSOR_POLARISATION_B,
              NC_RSD} radial_function_type;

enum Hermite_Interpolation_Order {HERMITE3, HERMITE4, HERMITE6};

/*************************************************************************************************************/
/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int transfer_functions_at_q(
                              struct transfer * ptr,
                              int index_md,
                              int index_ic,
                              int index_type,
                              int index_l,
                              double q,
                              double * ptransfer_local
                              );

  int transfer_init(
                    struct precision * ppr,
                    struct background * pba,
                    struct thermodynamics * pth,
                    struct perturbations * ppt,
                    struct fourier * pfo,
                    struct transfer * ptr
                    );

  int transfer_free(
                    struct transfer * ptr
                    );

  int transfer_indices(
                       struct precision * ppr,
                       struct perturbations * ppt,
                       struct transfer * ptr,
                       double q_period,
                       double K,
                       int sgnK
                       );

  int transfer_perturbation_copy_sources_and_nl_corrections(
                                                            struct perturbations * ppt,
                                                            struct fourier * pfo,
                                                            struct transfer * ptr,
                                                            double *** sources
                                                            );

  int transfer_perturbation_source_spline(
                                          struct perturbations * ppt,
                                          struct transfer * ptr,
                                          double *** sources,
                                          double *** sources_spline
                                          );

  int transfer_perturbation_sources_free(
                                         struct perturbations * ppt,
                                         struct fourier * pfo,
                                         struct transfer * ptr,
                                         double *** sources
                                         );

  int transfer_perturbation_sources_spline_free(
                                                struct perturbations * ppt,
                                                struct transfer * ptr,
                                                double *** sources_spline
                                                );

  int transfer_get_l_list(
                          struct precision * ppr,
                          struct perturbations * ppt,
                          struct transfer * ptr
                          );

  int transfer_get_q_list(
                          struct precision * ppr,
                          struct perturbations * ppt,
                          struct transfer * ptr,
                          double q_period,
                          double K,
                          int sgnK
                          );

  int transfer_get_q_limber_list(
                                 struct precision * ppr,
                                 struct perturbations * ppt,
                                 struct transfer * ptr,
                                 double K,
                                 int sgnK
                                 );

  int transfer_get_k_list(
                          struct perturbations * ppt,
                          struct transfer * ptr,
                          double K
                          );

  int transfer_get_source_correspondence(
                                         struct perturbations * ppt,
                                         struct transfer * ptr,
                                         int ** tp_of_tt
                                         );

  int transfer_free_source_correspondence(
                                          struct transfer * ptr,
                                          int ** tp_of_tt
                                          );

  int transfer_source_tau_size_max(
                                   struct precision * ppr,
                                   struct background * pba,
                                   struct perturbations * ppt,
                                   struct transfer * ptr,
                                   double tau_rec,
                                   double tau0,
                                   int * tau_size_max
                                   );

  int transfer_source_tau_size(
                               struct precision * ppr,
                               struct background * pba,
                               struct perturbations * ppt,
                               struct transfer * ptr,
                               double tau_rec,
                               double tau0,
                               int index_md,
                               int index_tt,
                               int * tau_size
                               );

  int transfer_compute_for_each_q(
                                  struct precision * ppr,
                                  struct background * pba,
                                  struct perturbations * ppt,
                                  struct transfer * ptr,
                                  int ** tp_of_tt,
                                  int index_q,
                                  int tau_size_max,
                                  double tau_rec,
                                  double *** sources,
                                  double *** sources_spline,
                                  double * window,
                                  struct transfer_workspace * ptw,
                                  short use_full_limber
                                  );

  int transfer_radial_coordinates(
                                  struct transfer * ptr,
                                  struct transfer_workspace * ptw,
                                  int index_md,
                                  int index_q
                                  );

  int transfer_interpolate_sources(
                                   struct perturbations * ppt,
                                   struct transfer * ptr,
                                   double k,
                                   int index_md,
                                   int index_ic,
                                   int index_type,
                                   double * sources,
                                   double * source_spline,
                                   double * interpolated_sources
                                   );

  int transfer_sources(
                       struct precision * ppr,
                       struct background * pba,
                       struct perturbations * ppt,
                       struct transfer * ptr,
                       double * interpolated_sources,
                       double tau_rec,
                       double k,
                       int index_md,
                       int index_tt,
                       double * sources,
                       double * window,
                       int tau_size_max,
                       double * tau0_minus_tau,
                       double * delta_tau,
                       int * tau_size_out
                       );

  int transfer_selection_function(
                                  struct precision * ppr,
                                  struct perturbations * ppt,
                                  struct transfer * ptr,
                                  int bin,
                                  double z,
                                  double * selection);

  int transfer_dNdz_analytic(
                             struct transfer * ptr,
                             double z,
                             double * dNdz,
                             double * dln_dNdz_dz);

  int transfer_selection_sampling(
                                  struct precision * ppr,
                                  struct background * pba,
                                  struct perturbations * ppt,
                                  struct transfer * ptr,
                                  int bin,
                                  double * tau0_minus_tau,
                                  int tau_size);

  int transfer_lensing_sampling(
                                struct precision * ppr,
                                struct background * pba,
                                struct perturbations * ppt,
                                struct transfer * ptr,
                                int bin,
                                double tau0,
                                double * tau0_minus_tau,
                                int tau_size);

  int transfer_source_resample(
                               struct precision * ppr,
                               struct background * pba,
                               struct perturbations * ppt,
                               struct transfer * ptr,
                               int bin,
                               double * tau0_minus_tau,
                               int tau_size,
                               int index_md,
                               double tau0,
                               double * interpolated_sources,
                               double * sources);

  int transfer_selection_times(
                               struct precision * ppr,
                               struct background * pba,
                               struct perturbations * ppt,
                               struct transfer * ptr,
                               int bin,
                               double * tau_min,
                               double * tau_mean,
                               double * tau_max);

  int transfer_selection_compute(
                                 struct precision * ppr,
                                 struct background * pba,
                                 struct perturbations * ppt,
                                 struct transfer * ptr,
                                 double * selection,
                                 double * tau0_minus_tau,
                                 double * delta_tau,
                                 int tau_size,
                                 double * pvecback,
                                 double tau0,
                                 int bin);

  int transfer_compute_for_each_l(
                                  struct transfer_workspace * ptw,
                                  struct precision * ppr,
                                  struct perturbations * ppt,
                                  struct transfer * ptr,
                                  int index_q,
                                  int index_md,
                                  int index_ic,
                                  int index_tt,
                                  int index_l,
                                  double l,
                                  double q_max_bessel,
                                  radial_function_type radial_type,
                                  short use_full_limber
                                  );

  int transfer_use_limber(
                          struct precision * ppr,
                          struct perturbations * ppt,
                          struct transfer * ptr,
                          double q_max_bessel,
                          int index_md,
                          int index_tt,
                          double q,
                          double l,
                          short * use_limber
                          );

  int transfer_integrate(
                         struct perturbations * ppt,
                         struct transfer * ptr,
                         struct transfer_workspace *ptw,
                         int index_q,
                         int index_md,
                         int index_tt,
                         double l,
                         int index_l,
                         double q,
                         radial_function_type radial_type,
                         double * trsf
                         );

  int transfer_limber(
                      struct transfer * ptr,
                      struct transfer_workspace * ptw,
                      int index_md,
                      int index_q,
                      double l,
                      double q,
                      radial_function_type radial_type,
                      double * trsf
                      );

  int transfer_limber_interpolate(
                                  struct transfer * ptr,
                                  double * tau0_minus_tau,
                                  double * sources,
                                  int tau_size,
                                  double tau0_minus_tau_limber,
                                  double * S
                                  );

  int transfer_limber2(
                       int tau_size,
                       struct transfer * ptr,
                       int index_md,
                       int index_q,
                       double l,
                       double q,
                       double * tau0_minus_tau,
                       double * sources,
                       radial_function_type radial_type,
                       double * trsf
                       );

  int transfer_can_be_neglected(
                                struct precision * ppr,
                                struct perturbations * ppt,
                                struct transfer * ptr,
                                int index_md,
                                int index_ic,
                                int index_tt,
                                double ra_rec,
                                double q,
                                double l,
                                short * neglect
                                );

  int transfer_late_source_can_be_neglected(
                                            struct precision * ppr,
                                            struct perturbations * ppt,
                                            struct transfer * ptr,
                                            int index_md,
                                            int index_tt,
                                            double l,
                                            short * neglect);

  int transfer_select_radial_function(
                                      struct perturbations * ppt,
                                      struct transfer * ptr,
                                      int index_md,
                                      int index_tt,
                                      radial_function_type *radial_type
                                      );

  int transfer_radial_function(
                               struct transfer_workspace * ptw,
                               struct perturbations * ppt,
                               struct transfer * ptr,
                               double k,
                               int index_q,
                               int index_l,
                               int x_size,
                               double * radial_function,
                               radial_function_type radial_type
                               );

  int transfer_init_HIS_from_bessel(
                                    struct transfer * ptr,
                                    HyperInterpStruct *pHIS
                                    );

  int transfer_global_selection_read(
                                     struct transfer * ptr
                                     );

  int transfer_workspace_init(
                              struct transfer * ptr,
                              struct precision * ppr,
                              struct transfer_workspace **ptw,
                              int perturbations_tau_size,
                              int tau_size_max,
                              double K,
                              int sgnK,
                              double tau0_minus_tau_cut,
                              HyperInterpStruct * pBIS
                              );

  int transfer_workspace_free(
                              struct transfer * ptr,
                              struct transfer_workspace *ptw
                              );

  int transfer_update_HIS(
                          struct precision * ppr,
                          struct transfer * ptr,
                          struct transfer_workspace * ptw,
                          int index_q,
                          double tau0
                          );

  int transfer_get_lmax(int (*get_xmin_generic)(int sgnK,
                                                int l,
                                                double nu,
                                                double xtol,
                                                double phiminabs,
                                                double *x_nonzero,
                                                int *fevals),
                        int sgnK,
                        double nu,
                        int *lvec,
                        int lsize,
                        double phiminabs,
                        double xmax,
                        double xtol,
                        int *index_l_left,
                        int *index_l_right,
                        ErrorMsg error_message);

  int transfer_precompute_selection(
                                    struct precision * ppr,
                                    struct background * pba,
                                    struct perturbations * ppt,
                                    struct transfer * ptr,
                                    double tau_rec,
                                    int tau_size_max,
                                    double ** window
                                    );

  int transfer_f_evo(
                     struct background* pba,
                     struct transfer * ptr,
                     double* pvecback,
                     int last_index,
                     double cotKgen,
                     double* f_evo
                     );

#ifdef __cplusplus
}
#endif

#endif
/* @endcond */
