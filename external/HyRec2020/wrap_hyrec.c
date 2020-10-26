/**
 * Small wrapper file for the hyrec code to be used in thermodynamics
 * Nils Schoeneberg Feb 2019
 * */

#include "common.h"
#include "thermodynamics.h"
#include "wrap_hyrec.h"

/* Boundaries and number of elements of temperature tables */
#define TR_MIN 0.004            /* Tr parameters */
#define TR_MAX 0.4
#define NTR    100
#define TM_TR_MIN 0.1           /* Tm/Tr parameters */
#define TM_TR_MAX 1.0
#define NTM 40

/* Forward-declare these just for convenience, so they are not in the .h file */
#define _HYREC_N_EXTRAPOLATION_ 30
int thermodynamics_hyrec_init(struct precision* ppr, struct background * pba, struct thermo * pth, double Nnow, double T_cmb, double fHe, double zstart_hyrec, struct thermohyrec* phy){

  if(phy->thermohyrec_verbose > 0){
    printf(" -> Using the hyrec wrapper programmed by Nils Sch. (Feb2019)\n");
    printf("    implements HyRec2 version Oct 2020 by Yacine Ali-Haimoud, Chris Hirata, and Nanoom Lee\n");
  }
  int i;
  double dN_safety;
  double Nur;
  class_alloc(phy->data,
              sizeof(HYREC_DATA),
              phy->error_message);
  
  HYREC_DATA *rec_data = phy->data;
 
  phy->N_LY = 3;
  phy->N_VIRT = NVIRT;
  phy->N_TM = NTM;
  phy->N_TR = NTR;

  phy->nH0 = Nnow*1e-6;
  phy->T_cmb = T_cmb;
  phy->T0 = phy->T_cmb;
  phy->fHe = fHe;

  phy->zend = 0.;
  if (MODEL == 4) phy->dlna = DLNA_SWIFT;
  else phy->dlna = DLNA_HYREC;
  dN_safety = 2*_HYREC_N_EXTRAPOLATION_; //Has to be >= _HYREC_N_EXTRAPOLATION_
  phy->zstart = (1.+zstart_hyrec)*exp(phy->dlna*dN_safety)-1.; //Make sure hyrec filled some values before first call from class (at z=8050)

  phy->nz = (long) floor(2.+log((1.+phy->zstart)/(1.+phy->zend))/phy->dlna);

  if(phy->thermohyrec_verbose > 1){
    printf("    Starting HyRec at z = %.10e until z = %.10e with %ld points\n",phy->zstart, phy->zend, phy->nz);
  }
  phy->data->path_to_hyrec = "external/HyRec2020/"; 
  hyrec_allocate(phy->data, phy->zstart, phy->zend);
  /* Error during allocation */
  if(phy->data->error != 0){
    class_call_message(phy->error_message,"hyrec_allocate",phy->data->error_message);
    return _FAILURE_;
  }

  phy->data->cosmo->Nmnu = pba->N_ncdm;
  if (pba->N_ncdm != 0) {
      for (i=0;i<pba->N_ncdm;i++) phy->data->cosmo->mnu[i] = pba->m_ncdm_in_eV[i];
  }
  phy->data->cosmo->T0 = phy->T0;
  phy->data->cosmo->obh2 = pba->Omega0_b*pba->h*pba->h;
  phy->data->cosmo->ocbh2 = (pba->Omega0_b+pba->Omega0_cdm)*pba->h*pba->h;
  phy->data->cosmo->onuh2 = pba->Omega0_ncdm_tot*pba->h*pba->h;
  
  phy->data->cosmo->YHe = pth->YHe;
  phy->data->cosmo->Nnueff = pba->Omega0_ur/pba->Omega0_g *8./7./pow(4./11.,4./3.);
  phy->data->cosmo->fHe = phy->fHe;              /* abundance of helium by number */
  phy->data->cosmo->fsR = 1.;
  phy->data->cosmo->meR = 1.;
  phy->data->cosmo->nH0 = phy->nH0;

  /* Mutiple massive neutrinos*/
  for (i=0;i<phy->data->cosmo->Nmnu;i++) {
	  phy->data->cosmo->mnu[i] = pba->m_ncdm_in_eV[i];
  }

  /* history of x_e */
  phy->xe_output    = create_1D_array(phy->nz, &phy->data->error, phy->data->error_message);

  /* history of photon occupation numbers */
  
  //Ratios of fine structure constant and electron mass at recombination compared to now
  phy->fsR = 1.;
  phy->meR = 1.;
  phy->xHeIII = phy->fHe;   /* Delta_xe = xHeIII here */
  phy->xHeII = 0.; //Uninitialized

  phy->izH0 = (long) floor(1 + log(_k_B_/_eV_*phy->T_cmb/phy->fsR/phy->fsR/phy->meR*(1.+phy->zstart)/(TR_MAX*0.99))/phy->dlna);
  //When T_gamma reaches TR_MAX (just approximately!), 1% error margin
  phy->zH0  = (1.+phy->zstart)*exp(-phy->izH0 * phy->dlna) - 1.;
  phy->data->rad->z0 = phy->zH0;
  
  phy->stage = 0;
  phy->saha_flag = 1;

  phy->filled_until_index_z = 0;

  phy->TR_prev = phy->T_cmb*(1+phy->zstart);
  phy->TM_prev = phy->TR_prev;
  phy->z_prev = (1+phy->zstart);
  phy->xHeIII_limit = 1e-8;

  return _SUCCESS_;
}

int thermodynamics_hyrec_free(struct thermohyrec* phy){
  
  hyrec_free (phy->data);

  return _SUCCESS_;
}

int hyrec_dx_H_dz(struct thermohyrec* phy, double x_H, double x_He, double xe, double nH, double z, double Hz, double Tmat, double Trad,
                  double dEdtdV_ion, double dEdtdV_lya, char* error_message, double *dx_H_dz) {
  long iz;
  double temp;
  int model;
  nH *= 1e-6;
  Tmat *= kBoltz;
  Trad *= kBoltz;
  if (MODEL == FULL) return 1; // or return error. FULL mode is not allowed in CLASS for now
  else iz = 0; // iz is only for FULL mode
  
  phy->data->cosmo->inj_params->ion = dEdtdV_ion;
  phy->data->cosmo->inj_params->exclya = dEdtdV_lya;

  double Trad_phys = Trad/phy->data->cosmo->fsR/phy->data->cosmo->fsR/phy->data->cosmo->meR;

  if (Trad_phys <= TR_MIN || Tmat/Trad <= TM_TR_MIN) { model = PEEBLES; }
  else { model = MODEL; }

  *dx_H_dz = -1./(1.+z)* rec_dxHIIdlna(phy->data, model, xe, x_H, nH, Hz, Tmat, Trad, iz, z);

  if(phy->data->error != 0){
    class_call_message(error_message,"rec_dxHIIdlna",phy->data->error_message);
    return _FAILURE_;
  }

  return _SUCCESS_;
}
int hyrec_dx_He_dz(struct thermohyrec* phy, double x_H, double x_He, double xe, double nH, double z, double Hz, double Tmat, double Trad,
                  double dEdtdV_ion, double dEdtdV_lya, char* error_message, double *dx_He_dz) {
  long iz;
  double xHeII = x_He*phy->data->cosmo->fHe;    // Different definitions between CLASS and HYREC-2
  double xH1;
  double temp;
 
  nH *= 1e-6;
  Tmat *= kBoltz;
  Trad *= kBoltz;
  
  if (MODEL == FULL) return 1; // or return error. FULL mode is not allowed in CLASS for now
  else iz = 0; // iz is only for FULL mode
  
  phy->data->cosmo->inj_params->ion = dEdtdV_ion;
  phy->data->cosmo->inj_params->exclya = dEdtdV_lya;
  /* XEII_MIN = 1e-6 defined in history.h
     HYREC-2 calculates Helium recombinations until xHeII ~ 1e-6 */
  if (xHeII<XHEII_MIN) {
    *dx_He_dz=0;
  }
  else {
    xH1 = rec_saha_xH1s(phy->data->cosmo, z, xHeII);
    *dx_He_dz = -1./(1.+z)* rec_helium_dxHeIIdlna(phy->data, z, xH1, xHeII, Hz) / phy->data->cosmo->fHe;
    if(phy->data->error != 0){
      class_call_message(error_message,"rec_helium_dxHeIIdlna",phy->data->error_message);
      return _FAILURE_;
    }
  }

  return _SUCCESS_;
}
