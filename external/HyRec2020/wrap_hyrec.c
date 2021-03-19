/**
 * Small wrapper file for the hyrec code to be used in thermodynamics
 * Nils Schoeneberg Feb 2019
 * */

#include "common.h"
#include "thermodynamics.h"
#include "wrap_hyrec.h"


/**
 * Initialize the thermohyrec structure, and in particular HyRec 2020.
 *
 * @param ppr   Input: pointer to precision structure
 * @param pba   Input: pointer to background structure
 * @param pth   Input: pointer to thermodynamics structure
 * @param phy Input/Output: pointer to thermohyrec structure
 * @return the error status
 */
int thermodynamics_hyrec_init(struct precision* ppr, struct background * pba, struct thermodynamics * pth, double Nnow, double T_cmb, double fHe, double zstart_hyrec, struct thermohyrec* phy){

  /** Summary: */

  if(phy->thermohyrec_verbose > 0){
    printf(" -> Using the hyrec wrapper programmed by Nils Sch. (Oct2020)\n");
    printf("    implements HyRec2 version Oct 2020 by Yacine Ali-Haimoud, Chris Hirata, and Nanoom Lee\n");
  }

  /** - allocate necessary data structures */
  class_alloc(phy->data,
              sizeof(HYREC_DATA),
              phy->error_message);

  phy->zend = 0.;
  phy->zstart = zstart_hyrec;

  if(phy->thermohyrec_verbose > 1){
    printf("    Starting HyRec at z = %.10e until z = %.10e\n",phy->zstart, phy->zend);
  }

  /** - allocate hyrec internally */
  phy->data->path_to_hyrec = ppr->hyrec_path;
  hyrec_allocate(phy->data, phy->zstart, phy->zend);
  /* Error during allocation */
  if(phy->data->error != 0){
    class_call_message(phy->error_message,"hyrec_allocate",phy->data->error_message);
    return _FAILURE_;
  }

  /** - set cosmological parameters for hyrec */
  phy->data->cosmo->T0 = T_cmb;
  phy->data->cosmo->obh2 = pba->Omega0_b*pba->h*pba->h;
  phy->data->cosmo->ocbh2 = (pba->Omega0_b+pba->Omega0_cdm)*pba->h*pba->h;

  phy->data->cosmo->YHe = pth->YHe;
  phy->data->cosmo->Neff = pba->Neff;
  phy->data->cosmo->fHe = fHe; /* abundance of helium relative to hydrogen by number */
  phy->data->cosmo->fsR = 1.;
  phy->data->cosmo->meR = 1.;
  phy->data->cosmo->nH0 = Nnow*1e-6;

  /** - set other parameters for hyrec */
  /* XEII_MIN = 1e-6 defined in history.h
     HYREC-2 calculates Helium recombinations only until xHeII ~ XEII_MIN */
  phy->xHeII_limit = XHEII_MIN;

  return _SUCCESS_;
}


/**
 * Free all memory space allocated by thermodynamics_hyrec_init
 *
 * @param pth Input/Output: pointer to thermohyrec structure (to be freed)
 * @return the error status
 */
int thermodynamics_hyrec_free(struct thermohyrec* phy){

  /* We just need to free hyrec (without error management) */
  hyrec_free(phy->data);
  free(phy->data);

  return _SUCCESS_;
}

/**
 * Calculate the derivative of the hydrogen HII ionization fraction
 *
 * @param pth   Input: pointer to thermodynamics structure
 * @param phy   Input: pointer to thermohyrec structure
 * @param x_H   Input: hydrogen HII ionization fraction
 * @param x_He  Input: helium HeIII ionization fraction
 * @param nH    Input: comoving total number of hydrogen atoms
 * @param z     Input: current cosmological redshift
 * @param Hz    Input: current value of hubble parameter in 1/s
 * @param Tmat  Input: temperature of baryons in Kelvin
 * @param Trad  Input: temperature of photons in Kelvin
 * @param dx_H_dz        Output: change in ionization fraction of hydrogen HII
 * @return the error status
 */
int hyrec_dx_H_dz(struct thermodynamics* pth, struct thermohyrec* phy, double x_H, double x_He, double xe, double nH, double z, double Hz, double Tmat, double Trad, double *dx_H_dz) {

  /** Summary: */

  /** - define local variables */
  struct injection* pin = &(pth->in);
  long iz = 0;
  int model;

  /** - do tests */
  class_test(MODEL == FULL, phy->error_message, "FULL mode is currently not allowed for HyRec-2");

  /** - assign variables */
  if (pth->has_exotic_injection == _TRUE_) {
    phy->data->cosmo->inj_params->ion = pin->pvecdeposition[pin->index_dep_ionH]/nH/(_E_H_ion_*_eV_);
    phy->data->cosmo->inj_params->exclya = pin->pvecdeposition[pin->index_dep_lya]/nH/(_E_H_lya_*_eV_);
  }
  else{
    phy->data->cosmo->inj_params->ion = 0.;
    phy->data->cosmo->inj_params->exclya = 0.;
  }

  double Trad_phys = Trad*kBoltz/phy->data->cosmo->fsR/phy->data->cosmo->fsR/phy->data->cosmo->meR;

  if (Trad_phys <= TR_MIN || Tmat/Trad <= T_RATIO_MIN) { model = PEEBLES; }
  else { model = MODEL; }

  /** - convert to correct units, and retrieve derivative */
  *dx_H_dz = -1./(1.+z)* rec_dxHIIdlna(phy->data, model, xe, x_H, nH*1e-6, Hz, Tmat*kBoltz, Trad*kBoltz, iz, z);

  /** - do error management */
  if(phy->data->error != 0){
    class_call_message(phy->error_message,"rec_dxHIIdlna",phy->data->error_message);
    return _FAILURE_;
  }

  return _SUCCESS_;
}

/**
 * Calculate the derivative of the helium HeIII ionization fraction
 *
 * @param pth   Input: pointer to thermodynamics structure
 * @param phy   Input: pointer to thermohyrec structure
 * @param x_H   Input: hydrogen HII ionization fraction
 * @param x_He  Input: helium HeIII ionization fraction
 * @param xe    Input: sum total ionization fraction
 * @param nH    Input: comoving total number of hydrogen atoms
 * @param z     Input: current cosmological redshift
 * @param Hz    Input: current value of hubble parameter in 1/s
 * @param Tmat  Input: temperature of baryons in Kelvin
 * @param Trad  Input: temperature of photons in Kelvin
 * @param dx_He_dz       Output: change in ionization fraction of helium HeIII
 * @return the error status
 */
int hyrec_dx_He_dz(struct thermodynamics* pth, struct thermohyrec* phy, double x_H, double x_He, double xe, double nH, double z, double Hz, double Tmat, double Trad, double *dx_He_dz) {

  /** Summary: */

  /** - define local variables */
  struct injection* pin = &(pth->in);
  double xHeII = x_He*phy->data->cosmo->fHe;    // Different definitions between CLASS and HYREC-2

  /** - do tests */
  class_test(MODEL == FULL, phy->error_message, "FULL mode is currently not allowed for HyRec-2");

  /** - assign variables */
  if (pth->has_exotic_injection == _TRUE_) {
    phy->data->cosmo->inj_params->ion = pin->pvecdeposition[pin->index_dep_ionHe]/nH/(_E_H_ion_*_eV_);
    phy->data->cosmo->inj_params->exclya = pin->pvecdeposition[pin->index_dep_lya]/nH/(_E_H_lya_*_eV_);
  }
  else{
    phy->data->cosmo->inj_params->ion = 0.;
    phy->data->cosmo->inj_params->exclya = 0.;
  }

  if (xHeII<phy->xHeII_limit) {
    /** - don't evolve He below the limit */
    *dx_He_dz=0;
  }
  else {

    /** - convert to correct units, and retrieve derivative */
    *dx_He_dz = -1./(1.+z)* rec_helium_dxHeIIdlna(phy->data, z, 1.-x_H, xHeII, Hz) / phy->data->cosmo->fHe;

    /** - do error management */
    if(phy->data->error != 0){
      class_call_message(phy->error_message,"rec_helium_dxHeIIdlna",phy->data->error_message);
      return _FAILURE_;
    }
  }

  return _SUCCESS_;
}
