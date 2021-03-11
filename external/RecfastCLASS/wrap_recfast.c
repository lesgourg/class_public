#include "wrap_recfast.h"
#include "thermodynamics.h"


/**
  *********************************************************************************************************************************
  * RECFAST is an integrator for Cosmic Recombination of Hydrogen and Helium, developed by Douglas Scott (dscott@astro.ubc.ca)
  * based on calculations in the paper Seager, Sasselov & Scott (ApJ, 523, L1, 1999) and "fudge" updates in Wong, Moss & Scott (2008).
  *
  * Permission to use, copy, modify and distribute without fee or royalty at any tier, this software and its documentation, for any
  * purpose and without fee or royalty is hereby granted, provided that you agree to comply with the following copyright notice and
  * statements, including the disclaimer, and that the same appear on ALL copies of the software and documentation,
  * including modifications that you make for internal use or for distribution:
  *
  * Copyright 1999-2010 by University of British Columbia.  All rights reserved.
  *
  * THIS SOFTWARE IS PROVIDED "AS IS", AND U.B.C. MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
  * BY WAY OF EXAMPLE, BUT NOT LIMITATION, U.B.C. MAKES NO REPRESENTATIONS OR WARRANTIES OF MERCHANTABILITY
  * OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE
  * ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.
  *********************************************************************************************************************************
  *
  * Version 1.5: includes extra fitting function from Rubino-Martin et al. arXiv:0910.4383v1 [astro-ph.CO]
  */

int recfast_init(struct precision* ppr,
                 struct background* pba,
                 struct thermodynamics* pth,
                 struct thermorecfast * pre,
                 enum recfast_photoion_modes recfast_photoion_mode,
                 double fHe) {

  /** Define local quantities */
  double Lalpha,Lalpha_He;

  /** - Import some thermodynamical quantities */
  pre->fHe = fHe;
  pre->photoion_mode = recfast_photoion_mode;

  /** - read a few precision/cosmological parameters */
  pre->AGauss1 = ppr->recfast_AGauss1;
  pre->AGauss2 = ppr->recfast_AGauss2;
  pre->zGauss1 = ppr->recfast_zGauss1;
  pre->zGauss2 = ppr->recfast_zGauss2;
  pre->wGauss1 = ppr->recfast_wGauss1;
  pre->wGauss2 = ppr->recfast_wGauss2;
  pre->Hswitch = ppr->recfast_Hswitch;
  pre->Heswitch = ppr->recfast_Heswitch;
  pre->x_H0_trigger2 = ppr->recfast_x_H0_trigger2;
  pre->x_He0_trigger2 = ppr->recfast_x_He0_trigger2;
  pre->fudge_H = ppr->recfast_fudge_H;
  pre->fudge_He = ppr->recfast_fudge_He;
  pre->x_H_limit_KHe = 0.9999999; /* threshold changed by Antony Lewis in 2008 to get smoother Helium */
  pre->x_H_limit_CfHe_t = 0.99999;
  pre->max_exp_boltz = 680.; /* The actual threshold is more around 709.8, but this leaves some wiggle-room */
  pre->x_He_trigger_small = 5.0e-9;

  pre->z_switch_late = ppr->recfast_z_switch_late;

  /** - Adjust fudging factors if needed */
  if (pre->Hswitch == _TRUE_){
    pre->fudge_H += ppr->recfast_delta_fudge_H;
  }


  /** - Assign quantities that have to be calculated first */
  /* These factors use inverse wavenumbers _L_ since they are most accurately measured */
  /* Lyman alpha wavelengths */
  Lalpha = 1./_L_H_alpha_;
  Lalpha_He = 1./_L_He_2p_;
  /* Ionization-lya temperature differences */
  pre->CDB = _h_P_*_c_*(_L_H_ion_-_L_H_alpha_)/_k_B_;
  pre->CDB_He = _h_P_*_c_*(_L_He1_ion_-_L_He_2s_)/_k_B_;
  /* Ionization temperatures */
  pre->CB1 = _h_P_*_c_*_L_H_ion_/_k_B_;           // equivalent to ptw->const_Tion_H
  pre->CB1_He1 = _h_P_*_c_*_L_He1_ion_/_k_B_;     // equivalent to ptw->const_Tion_HeI
  pre->CB1_He2 = _h_P_*_c_*_L_He2_ion_/_k_B_;     // equivalent to ptw->const_Tion_HeII
  /* Constants defined for the Peeble's factors  */
  pre->CR = 2.*_PI_*(_m_e_/_h_P_)*(_k_B_/_h_P_);  // equivalent to ptw->const_NR_numberdens
  pre->CK = pow(Lalpha,3)/(8.*_PI_);
  pre->CK_He = pow(Lalpha_He,3)/(8.*_PI_);
  /* Lyman alpha temperature */
  pre->CL = _c_*_h_P_/(_k_B_*Lalpha);
  pre->CL_He = _c_*_h_P_/(_k_B_/_L_He_2s_);
  pre->CL_Het = _h_P_*_c_*_L_He_2St_/_k_B_;
  /* Helium 2s->2p transition temperature*/
  pre->CDB_He2s2p = _h_P_*_c_*(_L_He_2p_-_L_He_2s_)/_k_B_;

  /** - Test schemes */
  /* He fudging */
  class_test((ppr->recfast_Heswitch < 0) || (ppr->recfast_Heswitch > 6),
             pre->error_message,
             "RECFAST error: unknown He fudging scheme");
  /* H fudging */
  class_test((ppr->recfast_Hswitch != _TRUE_) && (ppr->recfast_Hswitch != _FALSE_),
             pre->error_message,
             "RECFAST error: unknown H fudging scheme");

  return _SUCCESS_;

}

/**
 * The hydrogen switches are
 * - 0 => Only take normal K_H
 * - 1 => Add simple corections to K_H from gaussian fit
 * */
int recfast_dx_H_dz(struct thermodynamics* pth, struct thermorecfast * pre, double x_H, double x, double nH,
                    double z, double Hz, double Tmat, double Trad,
                    double* dxH_dz) {

  /** Define local variables */
  struct injection* pin = &(pth->in);
  /* new in recfast 1.4: */
  double Rup,Rdown,K,C,C_nofudge;
  double ion_H,ion_He,ion_lya;

  /** - Get necessary coefficients */
  Rdown = 1.e-19*_a_PPB_*pow((Tmat/1.e4),_b_PPB_)/(1.+_c_PPB_*pow((Tmat/1.e4),_d_PPB_));
  switch (pre->photoion_mode){
    case recfast_photoion_Tmat:
      Rup = Rdown * pow((pre->CR*Tmat),1.5)*exp(-pre->CDB/Tmat);
    break;
    default:
    case recfast_photoion_Trad:
      Rup = 1.e-19*_a_PPB_*pow((Trad/1.e4),_b_PPB_)/(1.+_c_PPB_*pow((Trad/1.e4),_d_PPB_)) * pow((pre->CR*Trad),1.5)*exp(-pre->CDB/Trad);
    break;
  }
  K = pre->CK/Hz;

  /* following is from recfast 1.5 */
  /** - Adjust the K constant with double gaussian fit */
  if (pre->Hswitch == _TRUE_ ){
    K *= 1.
      + pre->AGauss1*exp(-pow((log(1.+z)-pre->zGauss1)/pre->wGauss1,2))
      + pre->AGauss2*exp(-pow((log(1.+z)-pre->zGauss2)/pre->wGauss2,2));
  }
  /* end of new recfast 1.5 piece */

  /** - Calculate Peebles' coefficient */
  /* Peebles' coefficient (approximated as one when the Hydrogen
   * ionization fraction is very close to one) */
  if (x_H < pre->x_H0_trigger2 || z < pre->z_switch_late) {
    C = (1. + K*_Lambda_*nH*(1.-x_H))/(1./pre->fudge_H+K*_Lambda_*nH*(1.-x_H)/pre->fudge_H +K*Rup*nH*(1.-x_H));
    C_nofudge = (1. + K*_Lambda_*nH*(1.-x_H))/(1.+K*_Lambda_*nH*(1.-x_H) + K*Rup*nH*(1.-x_H));
  }
  else {
    C = 1.;
    C_nofudge = 1.;
  }

  /** - Evolve system by fudged Peebles' equation, use fudged Peebles' coefficient C */
  *dxH_dz = (x*x_H*nH*Rdown - Rup*(1.-x_H)*exp(-pre->CL/Tmat)) * C / (Hz*(1.+z));

  /** - Energy injection */
  if (pth->has_exotic_injection == _TRUE_) {
    ion_H = pin->pvecdeposition[pin->index_dep_ionH];
    ion_He = pin->pvecdeposition[pin->index_dep_ionHe];
    ion_lya = pin->pvecdeposition[pin->index_dep_lya];

    *dxH_dz += -1./nH*((ion_H+ion_He)/(_E_H_ion_*_eV_)+ion_lya*(1.-C_nofudge)/(_E_H_lya_*_eV_))/(Hz*(1.+z));
  }

  return _SUCCESS_;
}

/**
 * The helium switches are
 * - 0 => Only take normal K_He
 * - 1 => Add simple corections to K_He
 * - 2 => Add simple + doppler corrections to K_He
 * - 3 => Add simple + triple corrections to K_He, but not doppler
 * - 4 => Add simple + triple + triple doppler corrections to K_He, but not normal doppler
 * - 5 => Add simple + triple + doppler corrections to K_He, but not triple doppler
 * - 6 => Add simple + triple + doppler + triple doppler corrections to K_He
 * */
int recfast_dx_He_dz(struct thermodynamics* pth, struct thermorecfast * pre, double x_He, double x, double x_H, double nH,
                     double z, double Hz, double Tmat, double Trad,
                     double* dxHe_dz) {

  /** Define local variables */
  double Rdown_trip,Rup_trip,tauHe_s,pHe_s,tauHe_t,pHe_t,CL_PSt;
  double Doppler,gamma_2Ps,gamma_2Pt,pb,qb,AHcon;
  double sq_0,sq_1;
  double K_He,Rup_He,Rdown_He,He_Boltz;
  double CfHe_t=0.;
  double C_He;
  double n_He;
  int Heflag;

  /* This is just to prevent the compiler complaining:
     Technically Rdown_trip/Rup_trip should always be fine */
  Rdown_trip = 0.; Rup_trip = 0.;

  /** - Local variables and coefficients */
  n_He = pre->fHe * nH;
  sq_0 = sqrt(Tmat/_T_0_);
  sq_1 = sqrt(Tmat/_T_1_);

  Rdown_He = _a_VF_/(sq_0 * pow((1.+sq_0),(1.-_b_VF_)) * pow((1. + sq_1),(1. + _b_VF_)));
  switch (pre->photoion_mode){
    case recfast_photoion_Tmat:
      Rup_He = 4.*Rdown_He*pow((pre->CR*Tmat),1.5)*exp(-pre->CDB_He/Tmat);
      break;
    case recfast_photoion_Trad:
    default:
      Rup_He = 4.*_a_VF_/(sqrt(Trad/_T_0_) * pow((1.+sqrt(Trad/_T_0_)),(1.-_b_VF_)) * pow((1. + sqrt(Trad/_T_1_)),(1. + _b_VF_))) * pow((pre->CR*Trad),1.5)*exp(-pre->CDB_He/Trad);
      break;
  }

  /** - The K_He is calculated up to the required accuracy  */
  if ((x_He < pre->x_He_trigger_small) || (x_He > pre->x_He0_trigger2)){
    Heflag = 0;
  }
  else if(z<pre->z_switch_late){
    Heflag = 2;
  }
  else{
    Heflag = pre->Heswitch;
  }

  /* Simplest case : as defined in the equation */
  if (Heflag == 0){
    K_He = pre->CK_He/Hz;
  }
  /* More difficult cases : Take into account additional contributions */
  else {
    tauHe_s = _A2P_s_*pre->CK_He*3.*n_He*(1.-x_He)/Hz;
    pHe_s = (1.-exp(-tauHe_s))/tauHe_s;
    K_He = 1./(_A2P_s_*pHe_s*3.*n_He*(1.-x_He));

    /* In doppler mode (2) or in all mode (>=5), take doppler corrections */
    if (((Heflag == 2) || (Heflag >= 5)) && (x_H < pre->x_H_limit_KHe)) {

      Doppler = 2.*_k_B_*Tmat/(_m_H_*_not4_*_c_*_c_);
      Doppler = _c_*_L_He_2p_*sqrt(Doppler);
      gamma_2Ps = 3.*_A2P_s_*pre->fHe*(1.-x_He)*_c_*_c_
          /(sqrt(_PI_)*_sigma_He_2Ps_*8.*_PI_*Doppler*(1.-x_H))
          /pow(_c_*_L_He_2p_,2);
      pb = 0.36;
      qb = pre->fudge_He;
      AHcon = _A2P_s_/(1.+pb*pow(gamma_2Ps,qb));
      K_He=1./((_A2P_s_*pHe_s+AHcon)*3.*n_He*(1.-x_He));
    }

    /* In modes where triple He is added (>=3), calculate the triple He CfHe_t */
    if (Heflag >= 3) {
      Rdown_trip = _a_trip_/(sq_0*pow((1.+sq_0),(1.-_b_trip_)) * pow((1.+sq_1),(1.+_b_trip_)));
      switch (pre->photoion_mode){
        case recfast_photoion_Tmat:
          Rup_trip = Rdown_trip*exp(-_h_P_*_c_*_L_He2St_ion_/(_k_B_*Tmat))*pow(pre->CR*Tmat,1.5)*4./3.;
          break;
        case recfast_photoion_Trad:
        default:
          Rup_trip = _a_trip_/(sqrt(Trad/_T_0_)*pow((1.+sqrt(Trad/_T_0_)),(1.-_b_trip_)) * pow((1.+sqrt(Trad/_T_1_)),(1.+_b_trip_))) *exp(-_h_P_*_c_*_L_He2St_ion_/(_k_B_*Trad))*pow(pre->CR*Trad,1.5)*4./3.;
          break;
      }
      tauHe_t = _A2P_t_*n_He*(1.-x_He)*3./(8.*_PI_*Hz*pow(_L_He_2Pt_,3));
      pHe_t = (1. - exp(-tauHe_t))/tauHe_t;
      CL_PSt = _h_P_*_c_*(_L_He_2Pt_ - _L_He_2St_)/_k_B_;
      /* In triple He default mode, or mode 5, take simple term */
      if ((Heflag == 3) || (Heflag == 5) || (x_H >= pre->x_H_limit_CfHe_t)) {
        CfHe_t = _A2P_t_*pHe_t*exp(-CL_PSt/Tmat);
        CfHe_t = CfHe_t/(Rup_trip+CfHe_t);
      }
      /* In triple He doppler mode (4) or all mode (>=6), take doppler corrections */
      else {
        Doppler = 2.*_k_B_*Tmat/(_m_H_*_not4_*_c_*_c_);
        Doppler = _c_*_L_He_2Pt_*sqrt(Doppler);
        gamma_2Pt = 3.*_A2P_t_*pre->fHe*(1.-x_He)*_c_*_c_
          /(sqrt(_PI_)*_sigma_He_2Pt_*8.*_PI_*Doppler*(1.-x_H))
          /pow(_c_*_L_He_2Pt_,2);
        pb = 0.66;
        qb = 0.9;
        AHcon = _A2P_t_/(1.+pb*pow(gamma_2Pt,qb))/3.;
        CfHe_t = (_A2P_t_*pHe_t+AHcon)*exp(-CL_PSt/Tmat);
        CfHe_t = CfHe_t/(Rup_trip+CfHe_t);
      }
    }
  }

  /** - Final helium equations, again fudged Peebles' equations are used */
  if (x_He < 1.e-15){
    *dxHe_dz=0.;
  }
  else {

    /* Calculate first the boltzmann factor (limited for numerical reasons) */
    if (pre->CDB_He2s2p/Tmat < pre->max_exp_boltz){
      He_Boltz=exp(pre->CDB_He2s2p/Tmat);
    }
    else{
      He_Boltz=exp(pre->max_exp_boltz);
    }

    C_He = (1. + K_He*_Lambda_He_*n_He*(1.-x_He)*He_Boltz)/(1. + K_He*(_Lambda_He_+Rup_He)*n_He*(1.-x_He)*He_Boltz);

    /** Final He quations by Peebles with K_He*/
    *dxHe_dz = (x*x_He*nH*Rdown_He - Rup_He*(1.-x_He)*exp(-pre->CL_He/Tmat)) * C_He / (Hz*(1+z));

    /* following is from recfast 1.4 (now reordered) */
    /* this correction is not self-consistent when there is energy injection from dark matter, and leads to nan's at small redshift (unimportant when reionization takes over before that redshift) */
    /* Corrections due to triple helium, only when Heflag >=3 */
    if (Heflag >= 3){
      /** CfHe_t correction */
      *dxHe_dz += (x*x_He*nH*Rdown_trip - (1.-x_He)*3.*Rup_trip*exp(-pre->CL_Het/Tmat)) *CfHe_t/(Hz*(1.+z));
    }
    /* end of new recfast 1.4 piece */
  }

  /** - No He Energy injection */

  return _SUCCESS_;
}
