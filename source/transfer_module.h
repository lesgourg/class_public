#ifndef TRANSFER_MODULE_H
#define TRANSFER_MODULE_H

#include "input_module.h"
#include "base_module.h"

class TransferModule : public BaseModule {
public:
  TransferModule(InputModulePtr input_module, BackgroundModulePtr background_module, ThermodynamicsModulePtr thermodynamics_module, PerturbationsModulePtr perturbations_module, NonlinearModulePtr nonlinear_module);
  ~TransferModule();

  /** @name - number of modes and transfer function types */
  //@{

  int index_tt_t0_;      /**< index for transfer type = temperature (j=0 term) */
  int index_tt_t1_;      /**< index for transfer type = temperature (j=1 term) */
  int index_tt_t2_;      /**< index for transfer type = temperature (j=2 term) */
  int index_tt_e_;       /**< index for transfer type = E-polarization */
  int index_tt_b_;       /**< index for transfer type = B-polarization */
  int index_tt_lcmb_;    /**< index for transfer type = CMB lensing */
  int index_tt_density_; /**< index for first bin of transfer type = matter density */
  int index_tt_lensing_; /**< index for first bin of transfer type = galaxy lensing */

  int index_tt_rsd_;     /**< index for first bin of transfer type = redshift space distortion of number count */
  int index_tt_d0_;      /**< index for first bin of transfer type = doppler effect for of number count (j=0 term) */
  int index_tt_d1_;      /**< index for first bin of transfer type = doppler effect for of number count (j=1 term) */
  int index_tt_nc_lens_; /**< index for first bin of transfer type = lensing for of number count */
  int index_tt_nc_g1_;   /**< index for first bin of transfer type = gravity term G1 for of number count */
  int index_tt_nc_g2_;   /**< index for first bin of transfer type = gravity term G2 for of number count */
  int index_tt_nc_g3_;   /**< index for first bin of transfer type = gravity term G3 for of number count */
  int index_tt_nc_g4_;   /**< index for first bin of transfer type = gravity term G3 for of number count */
  int index_tt_nc_g5_;   /**< index for first bin of transfer type = gravity term G3 for of number count */

  int * tt_size_;     /**< number of requested transfer types tt_size[index_md] for each mode */
  //@}

  /** @name - number and list of multipoles */
  //@{
  int l_size_max_; /**< greatest of all l_size[index_md] */
  int ** l_size_tt_;  /**< number of multipole values for which we effectively compute the transfer function,l_size_tt[index_md][index_tt] */
  int * l_size_;   /**< number of multipole values for each requested mode, l_size[index_md] */
  int * l_;        /**< list of multipole values l[index_l] */
  //@}

  /** @name - number and list of wavenumbers */
  //@{
  size_t q_size_; /**< number of wavenumber values */
  double * q_;  /**< list of wavenumber values, q[index_q] */
  double ** k_; /**< list of wavenumber values for each requested mode, k[index_md][index_q]. In flat universes k=q. In non-flat universes q and k differ through q2 = k2 + K(1+m), where m=0,1,2 for scalar, vector, tensor. q should be used throughout the transfer module, excepted when interpolating or manipulating the source functions S(k,tau): for a given value of q this should be done in k(q). */
  int index_q_flat_approximation_; /**< index of the first q value using the flat rescaling approximation */
  //@}
  /** @name - transfer functions */
  //@{
  double ** transfer_; /**< table of transfer functions for each mode, initial condition, type, multipole and wavenumber, with argument transfer[index_md][((index_ic * transfer_module_->tt_size_[index_md] + index_tt) * transfer_module_->l_size_[index_md] + index_l) * transfer_module_->q_size_ + index_q] */
  //@}



private:
  int transfer_functions_at_q(int index_md, int index_ic, int index_type, int index_l, double q, double * ptransfer_local);
  int transfer_init();
  int transfer_free();
  int transfer_indices_of_transfers(double q_period, double K, int sgnK);
  int transfer_perturbation_copy_sources_and_nl_corrections(double *** sources);
  int transfer_perturbation_source_spline(double *** sources, double *** sources_spline);
  int transfer_perturbation_sources_free(double *** sources);
  int transfer_perturbation_sources_spline_free(double *** sources_spline);
  int transfer_get_l_list();
  int transfer_get_q_list(double q_period, double K, int sgnK);
  int transfer_get_q_list_v1(double q_period, double K, int sgnK);
  int transfer_get_k_list(double K);
  int transfer_get_source_correspondence(int ** tp_of_tt);
  int transfer_free_source_correspondence(int ** tp_of_tt);
  int transfer_source_tau_size_max(double tau_rec, double tau0, int * tau_size_max);
  int transfer_source_tau_size(double tau_rec, double tau0, int index_md, int index_tt, int * tau_size);
  int transfer_compute_for_each_q(int ** tp_of_tt, int index_q, int tau_size_max, double tau_rec, double *** sources, double *** sources_spline, double * window, struct transfer_workspace * ptw);
  int transfer_radial_coordinates(struct transfer_workspace * ptw, int index_md, int index_q);
  int transfer_interpolate_sources(int index_q, int index_md, int index_ic, int index_type, double * sources, double * source_spline, double * interpolated_sources);
  int transfer_sources(double * interpolated_sources, double tau_rec, int index_q, int index_md, int index_tt, double * sources, double * window, int tau_size_max, double * tau0_minus_tau, double * delta_tau, int * tau_size_out);
  int transfer_selection_function(int bin, double z, double * selection);
  int transfer_dNdz_analytic(double z, double * dNdz, double * dln_dNdz_dz);
  int transfer_selection_sampling(int bin, double * tau0_minus_tau, int tau_size);
  int transfer_lensing_sampling(int bin, double tau0, double * tau0_minus_tau, int tau_size);
  int transfer_source_resample(int bin, double * tau0_minus_tau, int tau_size, int index_md, double tau0, double * interpolated_sources, double * sources);
  int transfer_selection_times(int bin, double * tau_min, double * tau_mean, double * tau_max);
  int transfer_selection_compute(double * selection, double * tau0_minus_tau, double * delta_tau, int tau_size, double * pvecback, double tau0, int bin);
  int transfer_compute_for_each_l(struct transfer_workspace * ptw, int index_q, int index_md, int index_ic, int index_tt, int index_l, double l, double q_max_bessel, radial_function_type radial_type);
  int transfer_use_limber(double q_max_bessel, int index_md, int index_tt, double q, double l, short * use_limber);
  int transfer_integrate(struct transfer_workspace *ptw, int index_q, int index_md, int index_tt, double l, int index_l, double q, radial_function_type radial_type, double * trsf);
  int transfer_limber(struct transfer_workspace * ptw, int index_md, int index_q, double l, double q, radial_function_type radial_type, double * trsf);
  int transfer_limber_interpolate(double * tau0_minus_tau, double * sources, int tau_size, double tau0_minus_tau_limber, double * S);
  int transfer_limber2(int tau_size, int index_md, int index_q, double l, double q, double * tau0_minus_tau, double * sources, radial_function_type radial_type, double * trsf);
  int transfer_can_be_neglected(int index_md, int index_ic, int index_tt, double ra_rec, double q, double l, short * neglect);
  int transfer_late_source_can_be_neglected(int index_md, int index_tt, double l, short * neglect);
  int transfer_select_radial_function(int index_md, int index_tt, radial_function_type *radial_type);
  int transfer_radial_function(struct transfer_workspace * ptw, double k, int index_q, int index_l, int x_size, double * radial_function, radial_function_type radial_type);
  int transfer_init_HIS_from_bessel(HyperInterpStruct *pHIS);
  int transfer_global_selection_read();
  int transfer_workspace_init(struct transfer_workspace **ptw, int perturb_tau_size, int tau_size_max, double K, int sgnK, double tau0_minus_tau_cut, HyperInterpStruct * pBIS);
  int transfer_workspace_free(struct transfer_workspace *ptw);
  int transfer_update_HIS(struct transfer_workspace * ptw, int index_q, double tau0);
  int transfer_get_lmax(int (*get_xmin_generic)(int sgnK, int l, double nu, double xtol, double phiminabs, double *x_nonzero, int *fevals),
                        int sgnK, double nu, int * lvec, int lsize, double phiminabs, double xmax, double xtol, int * index_l_left, int * index_l_right, ErrorMsg error_message);
  int transfer_precompute_selection(double tau_rec, int tau_size_max, double ** window);
  int transfer_f_evo(double * pvecback, int last_index, double cotKgen, double * f_evo);

  BackgroundModulePtr background_module_;
  ThermodynamicsModulePtr thermodynamics_module_;
  PerturbationsModulePtr perturbations_module_;
  NonlinearModulePtr nonlinear_module_;

  short has_cls_; /**< copy of same flag in perturbation structure */
  int md_size_;       /**< number of modes included in computation */

  int nz_size_;           /**< number of redshift values in input tabulated selection function */
  double * nz_z_;         /**< redshift values in input tabulated selection function */
  double * nz_nz_;        /**< input tabulated values of selection function */
  double * nz_ddnz_;      /**< second derivatives in splined selection function*/

  int nz_evo_size_;            /**< number of redshift values in input tabulated evolution function */
  double * nz_evo_z_;          /**< redshift values in input tabulated evolution function */
  double * nz_evo_nz_;         /**< input tabulated values of evolution function */
  double * nz_evo_dlog_nz_;    /**< log of tabulated values of evolution function */
  double * nz_evo_dd_dlog_nz_; /**< second derivatives in splined log of evolution function */

};

#endif //TRANSFER_MODULE_H
