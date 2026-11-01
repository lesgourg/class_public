/** @file transfer.h Documented includes for transfer module. */

#ifndef __TRANSFER__
#define __TRANSFER__

#include "fourier.h"
#include "hyperspherical.h"
#include "errno.h"

/** maximum number of characters in a line of the file containing tabulated selection functions (selection=tabulated) */
#define _SELECTION_FILE_LINE_LENGTH_MAX_ 8192

/* macro: test if index_tt is in the range between index and index+num, while the flag is true */
#define _index_tt_in_range_(index,num,flag) (flag == _TRUE_) && (index_tt >= index) && (index_tt < index+num)

/**
 * enumeration of possible radial function types. This looks redundant with
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

/**
 * enumeration of the possible window functions \f$ W^X(\tau) \f$
 * multiplying the perturbation source \f$ S^{Y(X)}(k,\tau) \f$ in the
 * line-of-sight integral of a transfer type X
 */

enum transfer_window_type {
  window_none,           /**< W=1 (CMB temperature and polarisation types) */
  window_cmb_lensing,    /**< CMB lensing kernel with support [tau_rec, tau0] */
  window_selection,      /**< selection function of a redshift bin times a function of the background, with support [tau_min, tau_max] of the bin (local number count contributions) */
  window_integrated      /**< selection function of a redshift bin convolved with a kernel, with support [tau_min, tau0] (galaxy lensing and number count contributions integrated along the line of sight) */
};

/**
 * enumeration of the two ways of sampling the line-of-sight integrals in time
 */

enum transfer_sampling_type {
  sampling_from_perturbations, /**< the times of the perturbation module (ppt->tau_sampling) that are larger than tau_min */
  sampling_uniform,            /**< tau_size evenly spaced times between tau_min and tau_max (just tau_min if tau_size=1) */
  sampling_logarithmic         /**< the times of the perturbation module larger than tau_min, as long as their step is smaller than (tau0-tau)*log_step; then times with logarithmic steps Delta tau = (tau0-tau)*log_step until tau_max; then tau0 */
};

/**
 * enumeration of the possible integration methods for computing a transfer function
 * at a given wavenumber and multipole
 */

enum transfer_integration_method {
  transfer_set_zero,             /**< the transfer function vanishes or can be neglected */
  transfer_limber_approximation, /**< Limber approximation */
  transfer_line_of_sight         /**< line-of-sight integral of the transfer source times a radial function */
};

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

  FileName selection_file_name;    /**< for selection=tabulated: name of the file containing the redshift z in the first column and the selection function dN/dz of each bin in the next columns */
  int selection_file_z_size;       /**< for selection=tabulated: number of redshift values in this file */
  double * selection_file_z;       /**< for selection=tabulated: redshift values */
  double * selection_file_dNdz;    /**< for selection=tabulated: selection function of each bin, selection_file_dNdz[bin*selection_file_z_size+index_z] */
  double * selection_file_dddNdz;  /**< for selection=tabulated: second derivatives of the previous array with respect to z, for spline interpolation */
  double selection_z_min[_SELECTION_NUM_MAX_]; /**< for selection=tabulated: smallest redshift at which the selection function of each bin is non-zero (set in the input module) */
  double selection_z_max[_SELECTION_NUM_MAX_]; /**< for selection=tabulated: largest redshift at which the selection function of each bin is non-zero (set in the input module) */

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
  int index_tt_t2;      /**< index for transfer type = temperature (j=2 term) for scalar, vector and tensor modes */
  int index_tt_e;       /**< index for transfer type = E-polarization (for scalar, vector and tensor modes) */
  int index_tt_b;       /**< index for transfer type = B-polarization for tensor modes (m=2) */
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

  int index_tt_vector_t1;    /**< index for transfer type = temperature (j=1, m=1 term) for vector mode */
  int index_tt_vector_b;     /**< index for transfer type = B-polarization for vector modes (m=1) */

  int * tt_size;     /**< number of requested transfer types tt_size[index_md] for each mode */

  //@}

  /** @name - properties of each transfer type X, defined once for all in transfer_type_properties() */

  //@{

  int ** tp_of_tt;                     /**< index of the perturbation source Y(X) used by each transfer type, tp_of_tt[index_md][index_tt] */
  int ** bin_of_tt;                    /**< redshift bin of the number count and galaxy lensing types (-1 for CMB types), bin_of_tt[index_md][index_tt] */
  radial_function_type ** radial_type; /**< type of radial function (Bessel function or derivatives) in the line-of-sight integral, radial_type[index_md][index_tt] */
  enum transfer_window_type ** window_type; /**< type of window function W(tau) multiplying the perturbation source, window_type[index_md][index_tt] */
  double ** k_power;                   /**< the transfer source W(tau) S(k,tau) is further multiplied by (k/k_pivot)^k_power, k_power[index_md][index_tt] */
  double ** k_pivot;                   /**< pivot wavenumber of the previous factor, k_pivot[index_md][index_tt] */
  double ** l_limber;                  /**< the Limber approximation is used for l >= l_limber[index_md][index_tt] (_HUGE_ if never) */
  double ** delta_k_neglect;           /**< the transfer function is set to zero for l < (q-delta_k_neglect)*r_a(tau_rec), delta_k_neglect[index_md][index_tt] (_HUGE_ if never) */
  short ** has_late_source_cut;        /**< can the source be neglected at late times (tau>tau_cut) for l > ppr->transfer_neglect_late_source? has_late_source_cut[index_md][index_tt] */
  short ** has_full_limber;            /**< is the type also computed in the full Limber q-sampling scheme (on the grid q_limber)? has_full_limber[index_md][index_tt] (only for the CMB lensing potential of scalar modes, when ptr->do_lcmb_full_limber is true) */

  //@}

  /** @name - time samplings of the line-of-sight integrals, defined once for all in transfer_time_samplings() */

  //@{

  int * ts_size;                /**< number of distinct time sampling schemes for each mode, ts_size[index_md]; several types can share a time sampling */
  int ** ts_of_tt;              /**< index of the time sampling used by each type, ts_of_tt[index_md][index_tt] */
  int ** tau_size_ts;           /**< number of sampled times, tau_size_ts[index_md][index_ts] */
  int tau_size_max;             /**< largest value in tau_size_ts */
  double *** tau0_minus_tau;    /**< sampled values of (tau0-tau) in decreasing order, tau0_minus_tau[index_md][index_ts][index_tau] */
  double *** w_trapz;           /**< corresponding weights of the trapezoidal integration rule, w_trapz[index_md][index_ts][index_tau] */
  int *** index_tau_pt;         /**< for each sampled time, index of the time in ppt->tau_sampling just before it (or equal to it): used for resampling perturbation sources by linear interpolation, index_tau_pt[index_md][index_ts][index_tau] */
  double *** weight_tau_pt;     /**< corresponding weight of the next time in ppt->tau_sampling in the linear interpolation, weight_tau_pt[index_md][index_ts][index_tau] */
  double ** q_max_sampled;      /**< largest wavenumber q for which the Bessel functions are sampled with at least ppr->selection_sampling_bessel_min points per radian, q_max_sampled[index_md][index_ts] (_HUGE_ for samplings that are not evenly spaced, which are not designed for this criterion) */
  double *** window;            /**< window function W(tau) by which the perturbation source is multiplied, sampled like the type, window[index_md][index_tt][index_tau] */

  //@}

  /** @name - conformal times delimiting the selection function of each redshift bin, computed once for all in transfer_selection_times() */

  //@{

  double selection_tau_min[_SELECTION_NUM_MAX_];  /**< smallest time at which the selection function of each bin is non-zero */
  double selection_tau_mean[_SELECTION_NUM_MAX_]; /**< time corresponding to the mean redshift of each bin */
  double selection_tau_max[_SELECTION_NUM_MAX_];  /**< largest time at which the selection function of each bin is non-zero */

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

  short do_lcmb_full_limber; /**< in this particular run, will we use the full Limber q-sampling scheme? */

  size_t q_size_limber; /**< number of wavenumber values corresponding to k up to k_max */

  double * q_limber;  /**< list of wavenumber values used in full Limber q-sampling scheme, q_limber[index_q] */

  double ** k_limber; /**< list of wavenumber values used in full Limber q-sampling scheme */

  //@}

  /** @name - transfer functions */

  //@{

  double ** transfer; /**< table of transfer functions for each mode, initial condition, type, multipole and wavenumber, with argument transfer[index_md][((index_ic * ptr->tt_size[index_md] + index_tt) * ptr->l_size[index_md] + index_l) * ptr->q_size + index_q] */

  double ** transfer_limber; /**< table of transfer functions used in full Limber q-sampling scheme */

  //@}

  /** @name - technical parameters */

  //@{

  short transfer_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */

  short is_allocated; /**< flag is set to true if allocated */

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

  int ic_size_max;               /**< maximum number of initial conditions over all modes */
  int tp_size_max;               /**< maximum number of perturbation source types over all modes */
  int tt_per_ts_max;             /**< maximum number of transfer types sharing a time sampling */

  double ** interpolated_sources; /**< interpolated_sources[index_ic*tp_size_max+index_tp][index_tau]:
                                     perturbation sources interpolated at the current value of k, on the time
                                     sampling of the perturbation module (only allocated for sources needed in this module) */
  short * is_interpolated;       /**< is_interpolated[index_ic*tp_size_max+index_tp]: has the previous source been
                                    interpolated at the current wavenumber and mode? */
  double * sources;              /**< sources[(index_ic*tt_per_ts_max+index_tt_of_ts)*tau_size_max+index_tau]: transfer sources
                                    W(tau) S(k,tau) (k/k_pivot)^k_power of the types using the current time sampling scheme (index_tt_of_ts runs over these types) */
  double * chi;                  /**< chi[index_tau]: value of argument of bessel
                                    function: k(tau0-tau) (flat case)
                                    or sqrt(|K|)(tau0-tau) (non-flat
                                    case) */
  double * cscKgen;              /**< cscKgen[index_tau]: useful trigonometric function */
  double * cotKgen;              /**< cotKgen[index_tau]: useful trigonometric function */

  //@}

  /** @name - quantities related to the current multipole l, shared by all types using the current time sampling */

  //@{

  enum transfer_integration_method * method; /**< method[index_ic*tt_per_ts_max+index_tt_of_ts]: method for computing the transfer function of each type using the current time sampling scheme */
  int * index_tau_max;           /**< index_tau_max[index_ic*tt_per_ts_max+index_tt_of_ts]: for transfer_line_of_sight, index of the last sampled time contributing to the line-of-sight integral */
  short * bessel_truncation;     /**< bessel_truncation[index_ic*tt_per_ts_max+index_tt_of_ts]: for transfer_line_of_sight, is the integral truncated by the lower cut of the Bessel function (at tau0_minus_tau_min_bessel)? */

  int bessel_size;               /**< number of times of the current time sampling (starting from the first one) at which Bessel functions have been interpolated */
  double rescale_argument;       /**< rescaling of the argument of the Bessel functions in the flat rescaling approximation (curved space) */
  double * Phi;                  /**< Bessel function at the radial coordinates of the current time sampling, in reverse order (i.e. order of increasing argument): Phi[bessel_size-1-index_tau] */
  double * dPhi;                 /**< first derivative of Bessel function, same order */
  double * d2Phi;                /**< second derivative of Bessel function, same order */
  double * chireverse;           /**< argument of the Bessel function, same order */
  double * rescale_function;     /**< rescaling function of the flat rescaling approximation in curved space, same order */
  double * radial_function;      /**< radial function of a given type, in the order of the time sampling */

  //@}

  /** @name - parameters defining the spatial curvature (copied from background structure) */

  //@{

  double K; /**< curvature parameter (see background module for details) */
  int sgnK; /**< 0 (flat), 1 (positive curvature, spherical, closed), -1 (negative curvature, hyperbolic, open) */

  //@}

  double tau0_minus_tau_cut; /**< critical value of (tau0-tau) in time cut approximation */
};

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

  int transfer_type_properties(
                               struct precision * ppr,
                               struct perturbations * ppt,
                               struct transfer * ptr
                               );

  int transfer_type_sampling(
                             struct precision * ppr,
                             struct perturbations * ppt,
                             struct transfer * ptr,
                             double tau_rec,
                             double tau0,
                             int index_md,
                             int index_tt,
                             enum transfer_sampling_type * sampling_type,
                             double * tau_min,
                             double * tau_max,
                             double * log_step,
                             int * tau_size
                             );

  int transfer_uniform_sampling(
                                double tau0,
                                double tau_min,
                                double tau_max,
                                int tau_size,
                                double * tau0_minus_tau
                                );

  int transfer_logarithmic_sampling(
                                    struct perturbations * ppt,
                                    double tau0,
                                    double tau_min,
                                    double tau_max,
                                    double log_step,
                                    double * tau0_minus_tau,
                                    int * tau_size,
                                    ErrorMsg error_message
                                    );

  int transfer_time_samplings(
                              struct precision * ppr,
                              struct background * pba,
                              struct perturbations * ppt,
                              struct transfer * ptr,
                              double tau_rec
                              );

  int transfer_compute_for_each_q(
                                  struct precision * ppr,
                                  struct background * pba,
                                  struct perturbations * ppt,
                                  struct transfer * ptr,
                                  int index_q,
                                  double tau_rec,
                                  double *** sources,
                                  double *** sources_spline,
                                  struct transfer_workspace * ptw,
                                  short use_full_limber
                                  );

  int transfer_radial_coordinates(
                                  struct transfer_workspace * ptw,
                                  double * tau0_minus_tau,
                                  int tau_size,
                                  double k
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
                       struct transfer * ptr,
                       double * interpolated_sources,
                       double k,
                       int index_md,
                       int index_tt,
                       double * sources
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

  int transfer_selection_times(
                               struct precision * ppr,
                               struct background * pba,
                               struct perturbations * ppt,
                               struct transfer * ptr
                               );

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

  int transfer_select_method(
                             struct precision * ppr,
                             struct transfer * ptr,
                             struct transfer_workspace * ptw,
                             int index_q,
                             int index_md,
                             int index_tt,
                             int index_l,
                             double l,
                             double q,
                             double ra_rec,
                             double q_max_bessel,
                             double tau0_minus_tau_min_bessel,
                             double * tau0_minus_tau,
                             int tau_size,
                             double * sources,
                             short use_full_limber,
                             enum transfer_integration_method * method,
                             int * index_tau_max,
                             short * bessel_truncation
                             );

  int transfer_use_limber(
                          struct transfer * ptr,
                          double q_max_bessel,
                          int index_md,
                          int index_tt,
                          double q,
                          double l,
                          short * use_limber
                          );

  int transfer_bessel_cut(
                          struct transfer * ptr,
                          struct transfer_workspace * ptw,
                          int index_q,
                          int index_l,
                          double l,
                          double k,
                          double * tau0_minus_tau_min_bessel
                          );

  int transfer_integrate(
                         struct transfer * ptr,
                         struct transfer_workspace *ptw,
                         int index_md,
                         int index_tt,
                         int index_l,
                         double k,
                         double * tau0_minus_tau,
                         double * w_trapz,
                         double * sources,
                         int index_tau_max,
                         short bessel_truncation,
                         double tau0_minus_tau_min_bessel,
                         double * trsf
                         );

  int transfer_limber(
                      struct transfer * ptr,
                      struct transfer_workspace * ptw,
                      double l,
                      double q,
                      double * tau0_minus_tau,
                      int tau_size,
                      double * sources,
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
                                struct transfer * ptr,
                                int index_md,
                                int index_tt,
                                double ra_rec,
                                double q,
                                double l,
                                short * neglect
                                );

  int transfer_late_source_can_be_neglected(
                                            struct precision * ppr,
                                            struct transfer * ptr,
                                            int index_md,
                                            int index_tt,
                                            double l,
                                            short * neglect);

  int transfer_radial_function_needs(
                                     radial_function_type radial_type,
                                     short * need_Phi,
                                     short * need_dPhi,
                                     short * need_d2Phi
                                     );

  int transfer_bessel_interpolate(
                                  struct transfer_workspace * ptw,
                                  struct transfer * ptr,
                                  int index_q,
                                  int index_l,
                                  int x_size,
                                  short need_Phi,
                                  short need_dPhi,
                                  short need_d2Phi
                                  );

  int transfer_radial_function(
                               struct transfer_workspace * ptw,
                               struct transfer * ptr,
                               double k,
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

  int transfer_selection_file_read(
                                   char * file_name,
                                   int * z_size_out,
                                   int * bin_num_out,
                                   double ** z_out,
                                   double ** dNdz_out,
                                   double * z_mean,
                                   double * z_min,
                                   double * z_max,
                                   ErrorMsg error_message
                                   );

  int transfer_tabulated_selection_read(
                                        struct perturbations * ppt,
                                        struct transfer * ptr
                                        );

  int transfer_workspace_init(
                              struct transfer * ptr,
                              struct perturbations * ppt,
                              struct transfer_workspace *ptw,
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

  int transfer_precompute_windows(
                                  struct precision * ppr,
                                  struct background * pba,
                                  struct perturbations * ppt,
                                  struct transfer * ptr,
                                  double tau_rec
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
