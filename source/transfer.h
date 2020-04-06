/** @file transfer.h Documented includes for transfer module. */

#ifndef __TRANSFER__
#define __TRANSFER__

#include "nonlinear.h"
#include "hyperspherical.h"

/* macro: test if index_tt is in the range between index and index+num, while the flag is true */
#define _index_tt_in_range_(index,num,flag) (flag == _TRUE_) && (index_tt >= index) && (index_tt < index+num)
/* macro: test if index_tt corresponds to an integrated nCl/sCl contribution */
#define _integrated_ncl_ (_index_tt_in_range_(index_tt_lensing_, ppt->selection_num, ppt->has_cl_lensing_potential)) || \
          (_index_tt_in_range_(index_tt_nc_lens_, ppt->selection_num, ppt->has_nc_lens)) || \
          (_index_tt_in_range_(index_tt_nc_g4_,   ppt->selection_num, ppt->has_nc_gr)) || \
          (_index_tt_in_range_(index_tt_nc_g5_,   ppt->selection_num, ppt->has_nc_gr))
/* macro: test if index_tt corresponds to an non-integrated nCl/sCl contribution */
#define _nonintegrated_ncl_ (_index_tt_in_range_(index_tt_density_, ppt->selection_num, ppt->has_nc_density)) || \
          (_index_tt_in_range_(index_tt_rsd_,     ppt->selection_num, ppt->has_nc_rsd)) || \
          (_index_tt_in_range_(index_tt_d0_,      ppt->selection_num, ppt->has_nc_rsd)) || \
          (_index_tt_in_range_(index_tt_d1_,      ppt->selection_num, ppt->has_nc_rsd)) || \
          (_index_tt_in_range_(index_tt_nc_g1_,   ppt->selection_num, ppt->has_nc_gr))  || \
          (_index_tt_in_range_(index_tt_nc_g2_,   ppt->selection_num, ppt->has_nc_gr))  || \
          (_index_tt_in_range_(index_tt_nc_g3_,   ppt->selection_num, ppt->has_nc_gr))
/* macro: bin number associated to particular redshift bin and selection function for non-integrated contributions*/
#define _get_bin_nonintegrated_ncl_(index_tt)                                                      \
      if (_index_tt_in_range_(index_tt_density_, ppt->selection_num, ppt->has_nc_density))     \
        bin = index_tt - index_tt_density_;                                                    \
      if (_index_tt_in_range_(index_tt_rsd_,     ppt->selection_num, ppt->has_nc_rsd))         \
        bin = index_tt - index_tt_rsd_;                                                        \
      if (_index_tt_in_range_(index_tt_d0_,      ppt->selection_num, ppt->has_nc_rsd))         \
        bin = index_tt - index_tt_d0_;                                                         \
      if (_index_tt_in_range_(index_tt_d1_,      ppt->selection_num, ppt->has_nc_rsd))         \
        bin = index_tt - index_tt_d1_;                                                         \
      if (_index_tt_in_range_(index_tt_nc_g1_,   ppt->selection_num, ppt->has_nc_gr))          \
        bin = index_tt - index_tt_nc_g1_;                                                      \
      if (_index_tt_in_range_(index_tt_nc_g2_,   ppt->selection_num, ppt->has_nc_gr))          \
        bin = index_tt - index_tt_nc_g2_;                                                      \
      if (_index_tt_in_range_(index_tt_nc_g3_,   ppt->selection_num, ppt->has_nc_gr))          \
        bin = index_tt - index_tt_nc_g3_;
/* macro: bin number associated to particular redshift bin and selection function for integrated contributions*/
#define _get_bin_integrated_ncl_(index_tt)                                                               \
      if (_index_tt_in_range_(index_tt_lensing_, ppt->selection_num, ppt->has_cl_lensing_potential)) \
        bin = index_tt - index_tt_lensing_;                                                          \
      if (_index_tt_in_range_(index_tt_nc_lens_, ppt->selection_num, ppt->has_nc_lens))              \
        bin = index_tt - index_tt_nc_lens_;                                                          \
      if (_index_tt_in_range_(index_tt_nc_g4_,   ppt->selection_num, ppt->has_nc_gr))                \
        bin = index_tt - index_tt_nc_g4_;                                                            \
      if (_index_tt_in_range_(index_tt_nc_g5_,   ppt->selection_num, ppt->has_nc_gr))                \
        bin = index_tt - index_tt_nc_g5_;
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

struct transfers {

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

  short has_nz_evo_file;      /**< Has dN/dz (evolution function) input file? */
  short has_nz_evo_analytic;  /**< Use analytic form for dN/dz (evolution function) distribution? */
  FileName nz_evo_file_name;  /**< dN/dz (evolution function) input file name */

  //@}



  /** @name - technical parameters */

  //@{

  short initialise_HIS_cache; /**< only true if we are using CLASS for setting up a cache of HIS structures */

  short transfer_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

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

#endif
