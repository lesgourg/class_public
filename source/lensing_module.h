#ifndef LENSING_MODULE_H
#define LENSING_MODULE_H

#include "input_module.h"
#include "base_module.h"

#include <map>
#include <string>
#include <vector>

class LensingModule : public BaseModule {
public:
  LensingModule(InputModulePtr input_module, SpectraModulePtr spectra_module);
  ~LensingModule();
  std::map<std::string, std::vector<double>> cl_output(int lmax) const;
  int lensing_cl_at_l(int l, double * cl_lensed) const;

  int l_unlensed_max_;    /**< last multipole in all calculations (same as in spectra module)*/
  int l_lensed_max_;    /**< last multipole at which lensed spectra are computed */

private:
  int lensing_init();
  int lensing_free();
  int lensing_indices();
  int lensing_lensed_cl_tt(double *ksi, double **d00, double *w8, int nmu);
  int lensing_lensed_cl_te(double *ksiX, double **d20, double *w8, int nmu);
  int lensing_lensed_cl_ee_bb(double *ksip, double *ksim, double **d22, double **d2m2, double *w8, int nmu);
  int lensing_addback_cl_tt(double *cl_tt);
  int lensing_addback_cl_te(double *cl_te);
  int lensing_addback_cl_ee_bb(double *cl_ee, double *cl_bb);
  int lensing_X000(double * mu, int num_mu, int lmax, double * sigma2, double ** X000);
  int lensing_Xp000(double * mu, int num_mu, int lmax, double * sigma2, double ** Xp000);
  int lensing_X220(double * mu, int num_mu, int lmax, double * sigma2, double ** X220);
  int lensing_X022(double * mu, int num_mu, int lmax, double * sigma2, double ** X022);
  int lensing_Xp022(double * mu, int num_mu, int lmax, double * sigma2, double ** Xp022);
  int lensing_X121(double * mu, int num_mu, int lmax, double * sigma2, double ** X121);
  int lensing_X132(double * mu, int num_mu, int lmax, double * sigma2, double ** X132);
  int lensing_X242(double * mu, int num_mu, int lmax, double * sigma2, double ** X242);
  int lensing_d00(double * mu, int num_mu, int lmax, double ** d00);
  int lensing_d11(double * mu, int num_mu, int lmax, double ** d11);
  int lensing_d1m1(double * mu, int num_mu, int lmax, double ** d1m1);
  int lensing_d2m2(double * mu, int num_mu, int lmax, double ** d2m2);
  int lensing_d22(double * mu, int num_mu, int lmax, double ** d22);
  int lensing_d20(double * mu, int num_mu, int lmax, double ** d20);
  int lensing_d31(double * mu, int num_mu, int lmax, double ** d3m1);
  int lensing_d3m1(double * mu, int num_mu, int lmax, double ** d3m1);
  int lensing_d3m3(double * mu, int num_mu, int lmax, double ** d3m3);
  int lensing_d40(double * mu, int num_mu, int lmax, double ** d40);
  int lensing_d4m2(double * mu, int num_mu, int lmax, double ** d4m2);
  int lensing_d4m4(double * mu, int num_mu, int lmax, double ** d4m4);

  /** @name - information on number of type of C_l's (TT, TE...) */

  //@{

  int has_tt_; /**< do we want lensed \f$ C_l^{TT}\f$? (T = temperature) */
  int has_ee_; /**< do we want lensed \f$ C_l^{EE}\f$? (E = E-polarization) */
  int has_te_; /**< do we want lensed \f$ C_l^{TE}\f$? */
  int has_bb_; /**< do we want \f$ C_l^{BB}\f$? (B = B-polarization) */
  int has_pp_; /**< do we want \f$ C_l^{\phi\phi}\f$? (\f$ \phi \f$ = CMB lensing potential) */
  int has_tp_; /**< do we want \f$ C_l^{T\phi}\f$? */
  int has_dd_; /**< do we want \f$ C_l^{dd}\f$? (d = matter density) */
  int has_td_; /**< do we want \f$ C_l^{Td}\f$? */
  int has_ll_; /**< do we want \f$ C_l^{ll}\f$? (l = lensing potential) */
  int has_tl_; /**< do we want \f$ C_l^{Tl}\f$? */

  int index_lt_tt_; /**< index for type \f$ C_l^{TT} \f$*/
  int index_lt_ee_; /**< index for type \f$ C_l^{EE} \f$*/
  int index_lt_te_; /**< index for type \f$ C_l^{TE} \f$*/
  int index_lt_bb_; /**< index for type \f$ C_l^{BB} \f$*/
  int index_lt_pp_; /**< index for type \f$ C_l^{\phi\phi} \f$*/
  int index_lt_tp_; /**< index for type \f$ C_l^{T\phi} \f$*/
  int index_lt_dd_; /**< index for type \f$ C_l^{dd} \f$*/
  int index_lt_td_; /**< index for type \f$ C_l^{Td} \f$*/
  int index_lt_ll_; /**< index for type \f$ C_l^{dd} \f$*/
  int index_lt_tl_; /**< index for type \f$ C_l^{Td} \f$*/

  int lt_size_; /**< number of \f$ C_l\f$ types requested */

  //@}

  /** @name - table of pre-computed C_l values, and related quantities */

  //@{


  /* interpolable version: */

  int l_size_;       /**< number of l values */

  int * l_max_lt_;    /**< last multipole (given as an input) at which
        we want to output \f$ C_l \f$'s for a given mode and type */

  double * l_;       /**< table of multipole values l[index_l] */
  double * cl_lens_; /**< table of anisotropy spectra for each
         multipole and types,
         cl[index_l * ple->lt_size + index_lt] */

  double * ddcl_lens_; /**< second derivatives for interpolation */

  //@}
  SpectraModulePtr spectra_module_;

};

#endif //LENSING_MODULE_H
