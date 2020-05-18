#ifndef SPECTRA_MODULE_H
#define SPECTRA_MODULE_H

#include "input_module.h"
#include "base_module.h"

#include <map>
#include <string>
#include <vector>

class SpectraModule : public BaseModule {
public:
  SpectraModule(InputModulePtr input_module, PerturbationsModulePtr perturbations_module, PrimordialModulePtr primordial_module_, NonlinearModulePtr nonlinear_module, TransferModulePtr transfer_module);
  ~SpectraModule();
  int spectra_cl_at_l(double l, double * cl, double ** cl_md, double ** cl_md_ic) const;
  std::map<std::string, int> cl_output_index_map() const;
  std::map<std::string, std::vector<double>> cl_output(int lmax) const;
  void cl_output_no_copy(int lmax, std::vector<double*>& output_pointers) const;

  /** @name - information on number of modes and pairs of initial conditions */

  //@{
  int md_size_;           /**< number of modes (scalar, tensor, ...) included in computation */
  int * ic_size_;         /**< for a given mode, ic_size[index_md] = number of initial conditions included in computation */
  int * ic_ic_size_;      /**< for a given mode, ic_ic_size[index_md] = number of pairs of (index_ic1, index_ic2) with index_ic2 >= index_ic1; this number is just N(N+1)/2  where N = ic_size[index_md] */
  short ** is_non_zero_; /**< for a given mode, is_non_zero[index_md][index_ic1_ic2] is set to true if the pair of initial conditions (index_ic1, index_ic2) are statistically correlated, or to false if they are uncorrelated */
  //@}

  /** @name - information on number of type of C_l's (TT, TE...) */

  //@{

  int has_tt_; /**< do we want \f$ C_l^{TT}\f$? (T = temperature) */
  int has_ee_; /**< do we want \f$ C_l^{EE}\f$? (E = E-polarization) */
  int has_te_; /**< do we want \f$ C_l^{TE}\f$? */
  int has_bb_; /**< do we want \f$ C_l^{BB}\f$? (B = B-polarization) */
  int has_pp_; /**< do we want \f$ C_l^{\phi\phi}\f$? (\f$ \phi \f$ = CMB lensing potential) */
  int has_tp_; /**< do we want \f$ C_l^{T\phi}\f$? */
  int has_ep_; /**< do we want \f$ C_l^{E\phi}\f$? */
  int has_dd_; /**< do we want \f$ C_l^{dd}\f$? (d = density) */
  int has_td_; /**< do we want \f$ C_l^{Td}\f$? */
  int has_pd_; /**< do we want \f$ C_l^{\phi d}\f$? */
  int has_ll_; /**< do we want \f$ C_l^{ll}\f$? (l = galaxy lensing potential) */
  int has_tl_; /**< do we want \f$ C_l^{Tl}\f$? */
  int has_dl_; /**< do we want \f$ C_l^{dl}\f$? */

  int index_ct_tt_; /**< index for type \f$ C_l^{TT} \f$*/
  int index_ct_ee_; /**< index for type \f$ C_l^{EE} \f$*/
  int index_ct_te_; /**< index for type \f$ C_l^{TE} \f$*/
  int index_ct_bb_; /**< index for type \f$ C_l^{BB} \f$*/
  int index_ct_pp_; /**< index for type \f$ C_l^{\phi\phi} \f$*/
  int index_ct_tp_; /**< index for type \f$ C_l^{T\phi} \f$*/
  int index_ct_ep_; /**< index for type \f$ C_l^{E\phi} \f$*/
  int index_ct_dd_; /**< first index for type \f$ C_l^{dd} \f$((d_size*d_size-(d_size-non_diag)*(d_size-non_diag-1)/2) values) */
  int index_ct_td_; /**< first index for type \f$ C_l^{Td} \f$(d_size values) */
  int index_ct_pd_; /**< first index for type \f$ C_l^{pd} \f$(d_size values) */
  int index_ct_ll_; /**< first index for type \f$ C_l^{ll} \f$((d_size*d_size-(d_size-non_diag)*(d_size-non_diag-1)/2) values) */
  int index_ct_tl_; /**< first index for type \f$ C_l^{Tl} \f$(d_size values) */
  int index_ct_dl_; /**< first index for type \f$ C_l^{dl} \f$(d_size values) */

  int ct_size_; /**< number of \f$ C_l \f$ types requested */
  int d_size_;      /**< number of bins for which density Cl's are computed */
  //@}

  int l_size_max_; /**< greatest of all l_size[index_md] */
  double * l_;    /**< list of multipole values l[index_l] */
  int ** l_max_ct_;    /**< last multipole (given as an input) at which
                         we want to output \f$ C_l\f$'s for a given mode and type;
                         l[index_md][l_size[index_md]-1] can be larger
                         than l_max[index_md], in order to ensure a
                         better interpolation with no boundary effects */
  int * l_max_;    /**< last multipole (given as an input) at which
                     we want to output \f$ C_l\f$'s for a given mode (maximized over types);
                     l[index_md][l_size[index_md]-1] can be larger
                     than l_max[index_md], in order to ensure a
                     better interpolation with no boundary effects */
  int l_max_tot_; /**< last multipole (given as an input) at which
                    we want to output \f$ C_l\f$'s (maximized over modes and types);
                    l[index_md][l_size[index_md]-1] can be larger
                    than l_max[index_md], in order to ensure a
                    better interpolation with no boundary effects */


private:
  int spectra_init();
  int spectra_free();
  int spectra_indices();
  int spectra_cls();
  int spectra_compute_cl(int index_md, int index_ic1, int index_ic2, int index_l, int cl_integrand_num_columns, double * cl_integrand, double * primordial_pk, double * transfer_ic1, double * transfer_ic2);
  int spectra_k_and_tau();
  /* deprecated functions (since v2.8) */
  int spectra_pk_at_z(enum linear_or_logarithmic mode, double z, double * output_tot, double * output_ic, double * output_cb_tot, double * output_cb_ic);
  int spectra_pk_at_k_and_z(double k, double z, double * pk, double * pk_ic, double * pk_cb, double * pk_cb_ic);
  int spectra_pk_nl_at_z(enum linear_or_logarithmic mode, double z, double * output_tot, double * output_cb_tot);
  int spectra_pk_nl_at_k_and_z(double k, double z, double * pk_tot, double * pk_cb_tot);
  int spectra_fast_pk_at_kvec_and_zvec(double * kvec, int kvec_size, double * zvec, int zvec_size, double * pk_tot_out, double * pk_cb_tot_out, int nonlinear);
  int spectra_sigma(double R, double z, double * sigma);
  int spectra_sigma_cb(double R, double z, double * sigma_cb);
  /* deprecated functions (since v2.1) */
  int spectra_tk_at_z(double z, double * output);
  int spectra_tk_at_k_and_z(double k, double z, double * output);

  PerturbationsModulePtr perturbations_module_;
  PrimordialModulePtr primordial_module_;
  NonlinearModulePtr nonlinear_module_;
  TransferModulePtr transfer_module_;

  int index_md_scalars_; /**< index for scalar modes */

  /** @name - table of pre-computed C_l values, and related quantities */

  //@{

  int * l_size_;   /**< number of multipole values for each requested mode, l_size[index_md] */


  double ** cl_;   /**< table of anisotropy spectra for each mode, multipole, pair of initial conditions and types, cl[index_md][(index_l * psp->ic_ic_size[index_md] + index_ic1_ic2) * psp->ct_size + index_ct] */
  double ** ddcl_; /**< second derivatives of previous table with respect to l, in view of spline interpolation */

  //@}


};

#endif //SPECTRA_MODULE_H
