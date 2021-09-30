/******************************************************************************************************/
/*                           HYREC-2: Hydrogen and Helium Recombination Code                          */
/*                      Written by Yacine Ali-Haimoud and Chris Hirata (2010-17)                      */
/*                         with contributions from Nanoom Lee (2020)                                  */
/*                                                                                                    */
/*         history.c: functions for numerical integration of the recombination history                */
/*                                                                                                    */
/*                                                                                                    */
/*         Revision history:                                                                          */
/*             - January 2020 : - added massive neutrino part in Hubble rate                          */
/*                              - added new mode, SWIFT                                               */
/*                              - separated DLNA in two cases, MODEL = FULL or else                   */
/*                              - added error_massage control                                         */
/*             -         2015: - added DM annihilation and 21 cm routines.                            */
/*                             - changed cosmological parameters input form                           */
/*                             - possibility to change fine structure constant/ electron mass         */
/*                             - nH0 now in cm^-3 (instead of m^-3 which was only used in helium.c)   */
/*             - October 2012: - added some wrapper functions for running CAMB with HyRec             */
/*                              (courtesy of Antony Lewis)                                            */
/*                              - possibility to change fine structure constant/ electron mass        */
/*             - May 2012:  - added explicit dependence on fine structure constant and electron mass  */
/*                          - modified call of rec_build_history                                      */
/*                            and improved numerical radiative transfer equations                     */
/*                            so the Lyman-lines spectrum can be extracted                            */
/*                          - split some functions for more clarity                                   */
/*             - November 2011: - extended integration down to z = 0 with Peeble's model for z < 20   */
/*                              - changed dTm/dlna so it can be called at all times                   */
/*             - January 2011: changed various switches (notably for post-Saha expansions)            */
/*                             so that they remain valid for arbitrary cosmologies                    */
/*             - written November 2010                                                                */
/******************************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <libgen.h>

#include "history.h"
#include "helium.h"

/*************************************************************************************
Hubble expansion rate in sec^-1.
*************************************************************************************/

#ifdef CAMB

/* Use the Hubble rate from CAMB */

extern double exported_dtauda(double *);

double rec_HubbleRate(REC_COSMOPARAMS *cosmo, double z) {
  double a;

  a = 1./(1.+z);
  /* conversion from d tau/ da in Mpc to H(z) in 1/s */
  return 1./(a*a)/exported_dtauda(&a) /3.085678e22 * 2.99792458e8;
}

HYREC_DATA rec_data;
double logstart, dlna;
long int Nz;
int firstTime=0;

void hyrec_init() {
  /* Allocate spaces for hyrec data */
  double zmax = 8000.;
  double zmin = 0.;
  char *buffer = (char *) malloc (4096);
  getcwd (buffer, 4096);
  chdir(HYRECPATH);
  rec_data.path_to_hyrec = "";
  hyrec_allocate(&rec_data, zmax, zmin);
  chdir(buffer);
  free(buffer);
}

void rec_build_history_camb_(const double* OmegaC, const double* OmegaB, const double* h0inp, const double* tcmb,
                             const double* yp, const double* num_nu, double *xe, double *Tm, long int* nz) {

  double zmax = 8000.;
  double zmin = 0.;
  double h = *h0inp/100.;
  double h2 = h*h;
  char sub_message[SIZE_ErrorM];

  /* To load tables only once */
  if (firstTime==0){
    hyrec_init();
    logstart = -log(1.+zmax);
    Nz = rec_data.Nz;
    firstTime =1;
  }

  if (*nz != Nz) {
    rec_data.error = 1;
    sprintf(sub_message, "  called from rec_build_history_camb_\n");
    strcat(sub_message, "  Wrong number of redshifts is given from CAMB\n");
    strcat(sub_message, "  For the default SWIFT model, Nz should be 2248 in hyrec.f90\n");
    strcat(sub_message, "  For the rest of models in HYREC-2, Nz should be 105859 in hyrec.f90\n");
    strcat(sub_message, "  Check the line 18 of history.h and line 22 of hyrec.f90\n");
    strcat(rec_data.error_message, sub_message);
    printf("\n%s\n",rec_data.error_message);
    hyrec_free(&rec_data);
    exit(1);
   }

  /* Defining HYREC-2 parameters from CAMB */
  /* Note that there are some parameters which don't have to be defined here
     since Hubble rate is directly given from CAMB: orh2, okh2, odeh2, onuh2, w0, wa. */
  rec_data.cosmo->h = h;
  rec_data.cosmo->T0 = *tcmb;
  rec_data.cosmo->obh2 = *OmegaB * h2;
  rec_data.cosmo->ocbh2 = (*OmegaB + *OmegaC) * h2;
  rec_data.cosmo->YHe = *yp;
  rec_data.cosmo->Neff = *num_nu;  /* effective number of neutrino species */
  rec_data.cosmo->fsR = rec_data.cosmo->meR = 1.;   /* Default: today's values */
  rec_data.cosmo->nH0 = 11.223846333047e-6*rec_data.cosmo->obh2*(1.-rec_data.cosmo->YHe);  // number density of hudrogen today in cm-3
  rec_data.cosmo->fHe = rec_data.cosmo->YHe/(1.-rec_data.cosmo->YHe)/3.97153;              // abundance of helium by number
  if (MODEL == SWIFT) rec_data.cosmo->dlna = DLNA_SWIFT;
  else rec_data.cosmo->dlna = DLNA_HYREC;
  dlna = rec_data.cosmo->dlna;

  /* It seems there are no parameters related to DM annihilation in CAMB
     So, following parameters (inj_params) are meaningless here          */
  rec_data.cosmo->inj_params->pann = 0.;
  rec_data.cosmo->inj_params->pann_halo = 0.;
  rec_data.cosmo->inj_params->ann_z = 1.;
  rec_data.cosmo->inj_params->ann_zmax = 1.;
  rec_data.cosmo->inj_params->ann_zmin = 1.;
  rec_data.cosmo->inj_params->ann_var = 1.;
  rec_data.cosmo->inj_params->ann_z_halo = 1.;
  rec_data.cosmo->inj_params->on_the_spot = 1;
  rec_data.cosmo->inj_params->decay = 0.;

  /* Primodial black hole parameters */
  rec_data.cosmo->inj_params->Mpbh = 1.;
  rec_data.cosmo->inj_params->fpbh = 0.;

  rec_data.cosmo->inj_params->odmh2 = *OmegaC * h2;

  /* This is a flag to use Hubble rate defined history.c
     This flag is for the HYREC-2 implementation in CLASS    */
  double Hubble_flag[1];
  Hubble_flag[0] = -1.;

  rec_data.xe_output = xe; rec_data.Tm_output = Tm;

  rec_build_history(&rec_data, MODEL, Hubble_flag);
  if (rec_data.error == 1) {
    printf("\n%s\n",rec_data.error_message);
    exit(1);
  }
}


double hyrec_xe_(double* a, double *xe){
  double loga = log(*a);
  int error;   /* error and error_message are meaningless here */
  char *error_message;
  error_message = malloc(SIZE_ErrorM);
  if (loga < logstart) return xe[0];
  return rec_interp1d(logstart, dlna, xe, Nz, loga, &error, error_message);
}

double hyrec_tm_(double* a, double *tm){
  double loga = log(*a);
  int error;  /* error and error_message are meaningless here */
  char *error_message;
  if (loga < logstart) return tm[0];
  return rec_interp1d(logstart, dlna, tm, Nz, loga, &error, error_message);
}

#else

/* Hyrec Hubble expansion rate. */

double rec_HubbleRate(REC_COSMOPARAMS *cosmo, double z) {
  double a = 1./(1.+z), y, Tnu0;
  int i;
  cosmo->onuh2 = 0;
  Tnu0 = cosmo->T0 * pow(4./11.,1./3.) ;
  if (cosmo->Nmnu == 0) cosmo->onuh2 = 0;
  else{
    /* Massive neutrino energy density using fitting function */
    for (i=0;i<cosmo->Nmnu;i++) {
      y = cosmo->mnu[i]/kBoltz/(Tnu0/a);
      cosmo->onuh2 = cosmo->onuh2 + (1.+0.317322*0.0457584*pow(y,3.47446+1.)
          + 2.05298*0.0457584*pow(y,3.47446-1.))/(1.+0.0457584*pow(y,3.47446))*(3.45e-8*Tnu0*Tnu0*Tnu0*Tnu0)*5.6822*2;
    }
  }

  /* Total density parameter, including curvature */
  double rho = cosmo->ocbh2 /a/a/a     /* Matter (baryon + CDM) */
             + cosmo->okh2 /a/a       /* Curvature */
             + cosmo->odeh2 * pow(1./a, 3*(1+cosmo->w0)) * exp(3*cosmo->wa* (log(1./a)-1.+a))  /* Dark energy */
             + cosmo->orh2 /a/a/a/a   /* Radiation (photons + massless neutrinos) */
             + cosmo->onuh2 /a/a/a/a; /* Massive neutrinos */
  /* Conversion to Hubble rate in sec-1 */

  return( 3.2407792896393e-18 * sqrt(rho) );
}

#endif

/*************************************************************************************************
Cosmological parameters Input/Output
*************************************************************************************************/

void rec_get_cosmoparam(FILE *fin, FILE *fout, REC_COSMOPARAMS *param) {
  double Omega_b, Omega_cb, Omega_k, y, Tnu0;
  int i;
  if (fout!=NULL) fprintf(fout, "Enter Hubble parameter (h) : \n");
  if(fscanf(fin, "%lg", &(param->h))!=1){if (fout!=NULL) fprintf(fout, "Error in rec_get_cosmoparam when reading parameter 'h'\n");exit(1);};
  if (fout!=NULL) fprintf(fout, "Enter CMB temperature today [Kelvin]: \n");
  if(fscanf(fin, "%lg", &(param->T0))!=1){if (fout!=NULL) fprintf(fout, "Error in rec_get_cosmoparam when reading parameter 'CMB temperature today [Kelvin]'\n");exit(1);};

  if (fout!=NULL) fprintf(fout, "Enter baryon density, Omega_b: \n");
  if(fscanf(fin, "%lg", &(Omega_b))!=1){if (fout!=NULL) fprintf(fout, "Error in rec_get_cosmoparam when reading parameter 'Omega_b'\n");exit(1);};
  if (fout!=NULL) fprintf(fout, "Enter matter (CDM+baryons) density, Omega_cb: \n");
  if(fscanf(fin, "%lg", &(Omega_cb))!=1){if (fout!=NULL) fprintf(fout, "Error in rec_get_cosmoparam when reading parameter 'Omega_cb'\n");exit(1);};
  if (fout!=NULL) fprintf(fout, "Enter curvature, Omega_k: \n");
  if(fscanf(fin, "%lg", &(Omega_k))!=1){if (fout!=NULL) fprintf(fout, "Error in rec_get_cosmoparam when reading parameter 'Omega_k'\n");exit(1);};
  if (fout!=NULL) fprintf(fout, "Enter dark energy equation of state parameters, w wa: ");
  if(fscanf(fin, "%lg %lg", &(param->w0), &(param->wa))!=2){if (fout!=NULL) fprintf(fout, "Error in rec_get_cosmoparam when reading parameters 'w wa'\n");exit(1);};

  if (fout!=NULL) fprintf(fout, "Enter number of massive neutrino species, Nmnu: \n");
  if(fscanf(fin, "%lg", &(param->Nmnu))!=1){if (fout!=NULL) fprintf(fout, "Error in rec_get_cosmoparam when reading parameter 'Nmnu'\n");exit(1);};
  if (fout!=NULL) fprintf(fout, "Enter mass of neutrino1, mnu: \n");
  if(fscanf(fin, "%lg", &(param->mnu[0]))!=1){if (fout!=NULL) fprintf(fout, "Error in rec_get_cosmoparam when reading parameter 'mnu1'\n");exit(1);};
  if (fout!=NULL) fprintf(fout, "Enter mass of neutrino2, mnu: \n");
  if(fscanf(fin, "%lg", &(param->mnu[1]))!=1){if (fout!=NULL) fprintf(fout, "Error in rec_get_cosmoparam when reading parameter 'mnu2'\n");exit(1);};
  if (fout!=NULL) fprintf(fout, "Enter mass of neutrino3, mnu: \n");
  if(fscanf(fin, "%lg", &(param->mnu[2]))!=1){if (fout!=NULL) fprintf(fout, "Error in rec_get_cosmoparam when reading parameter 'mnu3'\n");exit(1);};

  if (fout!=NULL) fprintf(fout, "Enter primordial helium mass fraction, Y: \n");
  if(fscanf(fin, "%lg", &(param->YHe))!=1){if (fout!=NULL) fprintf(fout, "Error in rec_get_cosmoparam when reading parameter 'Y'\n");exit(1);};
  if (fout!=NULL) fprintf(fout, "Enter total effective number of neutrino species, N_eff: \n");
  if(fscanf(fin, "%lg", &(param->Neff))!=1){if (fout!=NULL) fprintf(fout, "Error in rec_get_cosmoparam when reading parameter 'N_eff'\n");exit(1);};

  if (fout!=NULL) fprintf(fout, "ratio of fine structure constant at recombination to today's value, fsR: \n");
  if(fscanf(fin, "%lg", &(param->fsR))!=1){if (fout!=NULL) fprintf(fout, "Error in rec_get_cosmoparam when reading parameter 'fsR'\n");exit(1);};
  if (fout!=NULL) fprintf(fout, "ratio of electron mass at recombination to today's value, meR: \n");
  if(fscanf(fin, "%lg", &(param->meR))!=1){if (fout!=NULL) fprintf(fout, "Error in rec_get_cosmoparam when reading parameter 'meR'\n");exit(1);};

  param->Nur = param->Neff - param->Nmnu;
  param->orh2  = 4.48162687719e-7 *param->T0*param->T0*param->T0*param->T0 *(1. + 0.227107317660239 *param->Nur);
  param->ocbh2  = Omega_cb *param->h*param->h;
  param->obh2  = Omega_b *param->h*param->h;
  param->okh2  = Omega_k *param->h*param->h;

  Tnu0 = param->T0 * pow(4./11.,1./3.);
  param->onuh2 = 0.;
  if (param->Nmnu != 0){
    for (i=0;i<param->Nmnu;i++){
      y = param->mnu[i]/kBoltz/Tnu0;
      param->onuh2 = param->onuh2 + (1.+0.317322*0.0457584*pow(y,3.47446+1.)
                    + 2.05298*0.0457584*pow(y,3.47446-1.))/(1.+0.0457584*pow(y,3.47446))*(3.450e-8 *Tnu0*Tnu0*Tnu0*Tnu0)*5.6822*2;
    }
  }

  param->odeh2 = (1. -Omega_cb -Omega_k - param->orh2/param->h/param->h- param->onuh2/param->h/param->h)*param->h*param->h;
  param->nH0 = 11.223846333047e-6*param->obh2*(1.-param->YHe);  // number density of hudrogen today in cm-3
  param->fHe = param->YHe/(1-param->YHe)/3.97153;              // abundance of helium by number

  if (fout!=NULL) fprintf(fout, "dark matter annihilation parameter, in cm^3/s/GeV, pann: \n");
  if(fscanf(fin, "%lg", &(param->inj_params->pann))!=1){if (fout!=NULL) fprintf(fout, "Error in rec_get_cosmoparam when reading parameter 'pann'\n");exit(1);};
  if (fout!=NULL) fprintf(fout, "pann_halo: \n");
  if(fscanf(fin, "%lg", &(param->inj_params->pann_halo))!=1){if (fout!=NULL) fprintf(fout, "Error in rec_get_cosmoparam when reading parameter 'pann_halo'\n");exit(1);};
  if (fout!=NULL) fprintf(fout, "ann_z: \n");
  if(fscanf(fin, "%lg", &(param->inj_params->ann_z))!=1){if (fout!=NULL) fprintf(fout, "Error in rec_get_cosmoparam when reading parameter 'ann_z'\n");exit(1);};
  if (fout!=NULL) fprintf(fout, "ann_zmax: \n");
  if(fscanf(fin, "%lg", &(param->inj_params->ann_zmax))!=1){if (fout!=NULL) fprintf(fout, "Error in rec_get_cosmoparam when reading parameter 'ann_zmax'\n");exit(1);};
  if (fout!=NULL) fprintf(fout, "ann_zmin: \n");
  if(fscanf(fin, "%lg", &(param->inj_params->ann_zmin))!=1){if (fout!=NULL) fprintf(fout, "Error in rec_get_cosmoparam when reading parameter 'ann_zmin'\n");exit(1);};
  if (fout!=NULL) fprintf(fout, "ann_var: \n");
  if(fscanf(fin, "%lg", &(param->inj_params->ann_var))!=1){if (fout!=NULL) fprintf(fout, "Error in rec_get_cosmoparam when reading parameter 'ann_var'\n");exit(1);};
  if (fout!=NULL) fprintf(fout, "ann_z_halo: \n");
  if(fscanf(fin, "%lg", &(param->inj_params->ann_z_halo))!=1){if (fout!=NULL) fprintf(fout, "Error in rec_get_cosmoparam when reading parameter 'ann_z_halo'\n");exit(1);};
  if (fout!=NULL) fprintf(fout, "on_the_spot: \n");
  if(fscanf(fin, "%d", &(param->inj_params->on_the_spot))!=1){if (fout!=NULL) fprintf(fout, "Error in rec_get_cosmoparam when reading parameter 'on_the_spot'\n");exit(1);};
  if (fout!=NULL) fprintf(fout, "decay: \n");
  if(fscanf(fin, "%lg", &(param->inj_params->decay))!=1){if (fout!=NULL) fprintf(fout, "Error in rec_get_cosmoparam when reading parameter 'decay'\n");exit(1);};

  if (fout!=NULL) fprintf(fout, "Mpbh: \n");
  if(fscanf(fin, "%lg", &(param->inj_params->Mpbh))!=1){if (fout!=NULL) fprintf(fout, "Error in rec_get_cosmoparam when reading parameter 'Mpbh'\n");exit(1);};
  if (fout!=NULL) fprintf(fout, "fpbh: \n");
  if(fscanf(fin, "%lg", &(param->inj_params->fpbh))!=1){if (fout!=NULL) fprintf(fout, "Error in rec_get_cosmoparam when reading parameter 'fpbh'\n");exit(1);};

  param->inj_params->odmh2      = param->ocbh2 - param->obh2;

  if (MODEL == SWIFT) param->dlna = DLNA_SWIFT;
  else param->dlna = DLNA_HYREC;

  if (fout!=NULL) fprintf(fout, "\n");
}



/*****************************************************************************************
Matter temperature -- 1st order steady state, from Hirata 2008.
The input and output temperatures are in KELVIN.
Added December 2014: possibility for additional energy deposition dEdtdV in eV/s/cm^3.
******************************************************************************************/

double rec_Tmss(double z, double xe, REC_COSMOPARAMS *cosmo, double dEdtdV, double H) {

  double fsR = cosmo->fsR;
  double meR = cosmo->meR;
  double Tr  = cosmo->T0 *(1.+z);
  double nH  = cosmo->nH0 *cube(1.+z);

  /* Coefficient = 8 sigma_T a_r / (3 m_e c) */
  /* Here Tr, Tm are the actual (not rescaled) temperatures */
  double coeff  = fsR*fsR/meR/meR/meR*4.91466895548409e-22*Tr*Tr*Tr*Tr*xe/(1.+xe+cosmo->fHe)/H;
  double Tm = Tr/(1.+1./coeff)
              + (1.+2.*xe)/3.*dEdtdV/kBoltz /(1.5 *nH*(1.+xe+cosmo->fHe))/H /(1.+coeff);

  return Tm;
}

/******************************************************************************************
Matter temperature evolution derivative. Input and output temperatures are in KELVIN.
Added May 2012: when Tm = Tr, return -Tr (needed by CLASS)
Corrected June 2016: in the presence of heating, Tm can cross Tr, so added a criterion
for returning -Tr when quasi-steady state. Note: this is not very "clean", there should
be some flag for quasi-steady-state, will eventually fix.
Added December 2014: possibility of additional energy deposition dEdtdV in eV/s/cm^3.
******************************************************************************************/

double rec_dTmdlna(double z, double xe, double Tm, REC_COSMOPARAMS *cosmo, double dEdtdV, double H) {
  double fsR = cosmo->fsR;
  double meR = cosmo->meR;
  double Tr  = cosmo->T0 *(1.+z);
  double nH  = cosmo->nH0 *cube(1.+z);
  return ( (Tr/Tm-1.<1e-10  && Tr > 3000.)  ? -Tr :
          -2.*Tm + fsR*fsR/meR/meR/meR*4.91466895548409e-22*Tr*Tr*Tr*Tr*xe/(1.+xe+cosmo->fHe)*(Tr-Tm)/H
          + (1.+2.*xe)/3. *dEdtdV /kBoltz /(1.5 *nH*(1.+xe+cosmo->fHe))/H);
   /* Coefficient = 8 sigma_T a_r / (3 m_e c) */
   /* Here Tr, Tm are the actual (not rescaled) temperatures */
}

double Tm_implicit(double z, double xe, double Tm, REC_COSMOPARAMS *cosmo, double dEdtdV, double H, double DLNA) {
  double fsR = cosmo->fsR;
  double meR = cosmo->meR;
  double Tr  = cosmo->T0 *(1.+z);
  double nH  = cosmo->nH0 *cube(1.+z);
  double gamma = fsR*fsR/meR/meR/meR*4.91466895548409e-22*Tr*Tr*Tr*Tr*xe/(1.+xe+cosmo->fHe)/H;
  return (Tm + DLNA*(gamma*Tr + (1.+2.*xe)/3. *dEdtdV /kBoltz /(1.5 *nH*(1.+xe+cosmo->fHe))/H) )/ ( 1.+(2.+gamma)*DLNA );
}


/************************************************************************
Explicit integrator.
Inputs: deriv: derivative at current time step
        deriv_prev: derivatives at two previous timesteps
        (deriv_prev[0] = previous time, deriv_prev[1] = 2 timesteps ago)
Output: change in x per dt
*************************************************************************/
double hyrec_integrator(double deriv, double deriv_prev[2], double z) {
  double result;

  if (MODEL == FULL) result = 1.25 * deriv - 0.25 *deriv_prev[1];
  else result = 23./12.*deriv -16./12. * deriv_prev[0] + 5./12. *deriv_prev[1];

  // update derivatives
  deriv_prev[1] = deriv_prev[0];
  deriv_prev[0] = deriv;

  return result;
}

/**********************************************************************************************
Second-order explicit integrator for HeII abundance during HeII->HeI recombination.
Assumes hydrogen ionization is in Saha equilibrium (free electrons being provided by both H and He).
If xHeII is close enough to Saha equilibrium do a post-Saha expansion.
dxHeIIdlna_prev[0] = derivative at previous timestep
dxHeIIdlna_prev[1] = derivative 2 timesteps prior
z = input z
H = Hubble parameter(z)
***********************************************************************************************/

void rec_get_xe_next1_He(HYREC_DATA *data, double z_in, double *xHeII, double dxHeIIdlna_prev[2],
                         double *hubble_array, int flag) {

  REC_COSMOPARAMS *cosmo = data->cosmo;
  double Nz = data->Nz, dz = data->zmax/Nz;
  int *error = &data->error;

  double xH1, xH1_p, xH1_m, xHeIISaha, dxHeIISaha_dlna, DdxHeIIdlna_Dxe, dxHeIIdlna, z_out, Dxe, DLNA;
  char sub_message[SIZE_ErrorM];
  double H;
  if (flag==10) DLNA = cosmo->dlna;
  else DLNA = cosmo->dlna/10.;

  if (hubble_array[0]==-1.) H  = rec_HubbleRate(cosmo, z_in);
  else H = rec_interp1d(.0, dz, hubble_array, Nz, z_in, error, sub_message);

  xH1        = rec_saha_xH1s(cosmo, z_in, *xHeII);
  dxHeIIdlna = rec_helium_dxHeIIdlna(data, z_in, xH1, *xHeII, H);
  /* Post-Saha approximation during the early phase of HeII->HeI recombination */
  if (data->quasi_eq == 1) {
    z_out     = (1.+z_in)*exp(-DLNA)-1.;
    xHeIISaha = rec_saha_xHeII(cosmo, z_out);

    dxHeIISaha_dlna  = (1.+z_out)*(rec_saha_xHeII(cosmo, z_out-0.5)
                                  -rec_saha_xHeII(cosmo, z_out+0.5));

    Dxe    = 0.01*(cosmo->fHe - xHeIISaha);
    xH1_p = rec_saha_xH1s(cosmo, z_out, xHeIISaha+Dxe);
    xH1_m = rec_saha_xH1s(cosmo, z_out, xHeIISaha-Dxe);

    if (hubble_array[0]==-1.) H  = rec_HubbleRate(cosmo, z_out);
    else H = rec_interp1d(.0, dz, hubble_array, Nz, z_out, error, sub_message);

    DdxHeIIdlna_Dxe  = (rec_helium_dxHeIIdlna(data, z_out, xH1_p, xHeIISaha+Dxe, H)
                       -rec_helium_dxHeIIdlna(data, z_out, xH1_m, xHeIISaha-Dxe, H))/2./Dxe;

    *xHeII = xHeIISaha + dxHeIISaha_dlna/DdxHeIIdlna_Dxe;
    /* Check that the post-Saha approximation is still valid. If not, switch it off for future iterations */
    if (fabs(*xHeII - xHeIISaha) > DXHEII_MAX) data->quasi_eq = 0;
  }

  /* Otherwise integrate ODE */
  else *xHeII += DLNA * hyrec_integrator(dxHeIIdlna, dxHeIIdlna_prev,z_in);

  if (*error == 1) {
    sprintf(sub_message, "  called from rec_get_xe_next1_He at z = %f\n", z_in);
    strcat(data->error_message, sub_message);
  }
}


/***************************************************************************************************
Quasi-equilibrium solution for the hydrogen neutral density at z
Used at early stages when the system is stiff. After that switch to explicit second-order solver.
input:  xH1 at previous time step
output: xH1 at z.
iz_rad is the index for the radiation field at z
***************************************************************************************************/

void rec_xH1_stiff(HYREC_DATA *data, int model, double z, double xHeII, double *xH1, unsigned iz_rad, double H){

  REC_COSMOPARAMS *cosmo = data->cosmo;
  int *error = &data->error;

  double ainv, xH1sSaha, xHIISaha, dxH1sSaha_dlna, dxH1sdlna_Saha, DdxH1sdlna_DxH1s, T, nH, Dxe;
  int model_stiff;	// To use EMLA2p2s model for PostSaha in FULL and SWIFT mode
  char sub_message[SIZE_ErrorM];

  // Set model for rec_xH1_stiff. FULL mode uses EMLA2s2p for stiff.
  // SWIFT mode uses EMLA2s2p for stiff when z is not in the range of redshifts for SWIFT fitting function.

  T = kBoltz*cosmo->T0 * (ainv=1.+z);

  if (model == FULL) {
    model_stiff = EMLA2s2p;
  }
  else{
    model_stiff = model;
    if (model == SWIFT){
      if (T/kBoltz/cosmo->fsR/cosmo->fsR/cosmo->meR > data->fit->swift_func[0][DKK_SIZE-1] ) model_stiff = EMLA2s2p;
    }
  }

  nH = cosmo->nH0 *cube(1.+z);
  xH1sSaha = rec_saha_xH1s(cosmo, z, xHeII);
  xHIISaha = 1.-xH1sSaha;

  dxH1sSaha_dlna = (1.+z)*(rec_saha_xH1s(cosmo, z-0.5, xHeII)
                          -rec_saha_xH1s(cosmo, z+0.5, xHeII));
                  /* (partial xHII)/(partial lna). Use xH1s = 1-xHII for better accuracy. */
  if (xHeII != 0.) {
    dxH1sSaha_dlna = (1.+z)*(rec_saha_xH1s(cosmo, z-0.5, xHeII)
                            -rec_saha_xH1s(cosmo, z+0.5, xHeII));
    Dxe = 0.01*(cosmo->fHe - xHeII);
    dxH1sSaha_dlna += (rec_saha_xH1s(cosmo, z, xHeII+Dxe)
                      -rec_saha_xH1s(cosmo, z, xHeII-Dxe))/2./Dxe
                      *rec_helium_dxHeIIdlna(data, z, xH1sSaha, xHeII, H);
                   /* (partial xHII)/(partial xHeII).dxHeII/dlna */
  }

  dxH1sdlna_Saha = -rec_dxHIIdlna(data, model_stiff, xHIISaha + xHeII, xHIISaha, nH, H, T, T, iz_rad, z);
  Dxe            = 0.01*xH1sSaha;
  DdxH1sdlna_DxH1s = (rec_dxHIIdlna(data, model_stiff, xHIISaha+Dxe + xHeII, xHIISaha+Dxe, nH, H, T, T, iz_rad, z)
                       -rec_dxHIIdlna(data, model_stiff, xHIISaha-Dxe + xHeII, xHIISaha-Dxe, nH, H, T, T, iz_rad, z))/2./Dxe;
  *xH1 = xH1sSaha + (dxH1sSaha_dlna - dxH1sdlna_Saha)/DdxH1sdlna_DxH1s;


  if (fabs(*xH1 - xH1sSaha) > DXHII_MAX)  data->quasi_eq = 0;
  //if (z<1700.) *stiff = 0;   /* Used when calculating the correction function for SWIFT mode. */

  /* Update photon population when MODEL = FULL */
  if (model == FULL) rec_dxHIIdlna(data, model, xHeII + 1.-*xH1, 1.-*xH1, nH, H, T, T, iz_rad, z);

  if (*xH1 < 0. || *xH1 != *xH1) {
    sprintf(sub_message, "xH1 < 0 in rec_xH1_stiff: at z = %f, xH1 = %E\n", z, *xH1);
    strcat(data->error_message, sub_message);
    *error = 1;
  }

  if (*error == 1) {
    sprintf(sub_message, "  called from rec_xH1_stiff at z = %f\n", z);
    strcat(data->error_message, sub_message);
    return;
  }
  return;
}



/*****************************************************************************************************
Second-order integrator (unless stiff system) used to evolve simultaneously HI and HeII.
When hydrogen is close to Saha equilibrium but there is still a significant amount of HeII,
use a post-Saha expansion for hydrogen. The other way around is not treated (simply integrate H and He until
there is almost no HeII left, then integrate H only)
******************************************************************************************************/

void get_rec_next2_HHe(HYREC_DATA *data, int model, double z_in, long iz, double *xH1, double *xHeII,
                       double dxHIIdlna_prev[2], double dxHeIIdlna_prev[2], double H) {

  REC_COSMOPARAMS *cosmo = data->cosmo;
  int *error = &data->error;
  double Tm = data->Tm_output[iz-1];
  long iz_rad = iz-1-data->rad->iz_rad_0;

  double dxHeIIdlna, dxHIIdlna=0., z_out, xe;
  double nH, TR, DLNA;
  char sub_message[SIZE_ErrorM];
  DLNA = cosmo->dlna;
  xe = *xHeII + 1.- (*xH1);
  nH = cosmo->nH0 *cube(1.+z_in);
  TR = kBoltz * cosmo->T0 *(1.+z_in);
  /* Evolve HeII by solving ODE */
  dxHeIIdlna  = rec_helium_dxHeIIdlna(data, z_in, *xH1, *xHeII, H);
  *xHeII     += DLNA * hyrec_integrator(dxHeIIdlna, dxHeIIdlna_prev, z_in);

  /* Compute dxHII/dlna. This also correctly updates the radiation field at z_in,
     which is required even when using the stiff approximation */

  /* If system is stiff use the quasi-equilibrium solution */
  if(data->quasi_eq == 1){
    z_out = (1.+z_in)*exp(-DLNA)-1.;
    rec_xH1_stiff(data, model, z_out, *xHeII, xH1, iz_rad+1, H);
    hyrec_integrator(dxHIIdlna, dxHIIdlna_prev, z_in);
  }

    /* Otherwise use second-order explicit solver */
  else {
    dxHIIdlna = rec_dxHIIdlna(data, model, xe, 1.-(*xH1), nH, H, kBoltz*Tm, TR, iz_rad, z_in);
    *xH1 -= DLNA * hyrec_integrator(dxHIIdlna, dxHIIdlna_prev, z_in);
  }

    /* Checking for errors */
  if (*xH1 < 0. || *xH1 > 1. || *xH1 != *xH1) {
    sprintf(sub_message, "xH1 < 0 or xH1 > 1 in get_rec_next2_HHe: at z = %f, xH1 = %E\n", z_out, *xH1);
    strcat(data->error_message, sub_message);
    *error = 1;
  }
  if (*error == 1) {
    sprintf(sub_message, "  called from get_rec_next2_HHe at z = %f\n", z_out);
    strcat(data->error_message, sub_message);
    return;
  }
}

/*********************************************************************************************************
Second-order integrator (unless stiff system) to evolve hydrogen only, assuming helium has entirely recombined.
Tm is given as an input (to avoid computing it twice) and fixed to quasi-equilibrium value with Tr.
Input : xe [at z_in]
Output: xe [at next timestep]
**********************************************************************************************************/

void rec_get_xe_next1_H(HYREC_DATA *data, int model, double z_in, long iz, double xe_in, double Tm_in,
                        double *xe_out, double *Tm_out, double dxedlna_prev[2], double H, int flag) {

  REC_COSMOPARAMS *cosmo = data->cosmo;
  int *error = &data->error;
  long iz_rad = iz-1-data->rad->iz_rad_0;

  double dxedlna, z_out;
  double nH, TR, xH1, dEdtdV, DLNA;
  char sub_message[SIZE_ErrorM];
  if (flag==10) DLNA = cosmo->dlna;
  else DLNA = cosmo->dlna/10.;

  z_out   = (1.+z_in)*exp(-DLNA)-1.;
  nH = cosmo->nH0 *cube(1.+z_in);
  TR = kBoltz *cosmo->T0 *(1.+z_in);

  dEdtdV = cosmo->inj_params->ion*3*nH*EI / (1-xe_in);
  /* Compute dxHII/dlna. This also correctly updates the radiation field at z_in,
     which is required even when using the stiff approximation */

  /* If system is stiff use the quasi-equilibrium solution */
  if (data->quasi_eq == 1) {
    xH1 = 1.-xe_in;
    rec_xH1_stiff(data, model, z_out, 0, &xH1, iz_rad+1, H);
    *xe_out = 1.-xH1;
    dxedlna_prev[1] = dxedlna_prev[0];
    dxedlna_prev[0] = (*xe_out-xe_in)/DLNA;
  }
    /* Otherwise use second-order explicit solver */
  else {
    dxedlna = rec_dxHIIdlna(data, model, xe_in, xe_in, nH, H, kBoltz*Tm_in, TR, iz_rad, z_in);
    *xe_out = xe_in + DLNA * hyrec_integrator(dxedlna, dxedlna_prev, z_in);
  }

  /* Quasi-steady state solution for Tm */
  *Tm_out = rec_Tmss(z_out, *xe_out, cosmo, dEdtdV, H);

  // Test that the outcome is sensible
  if (*xe_out > 1. || *xe_out < 0. || *xe_out != *xe_out) {
    sprintf(sub_message, "xe > 1 or xe < 0 in get_rec_next1_H at z = %E, xe = %E\n", z_out, *xe_out);
    strcat(data->error_message, sub_message);
    *error = 1;
  }
  if (*error == 1) {
    sprintf(sub_message, "  called from get_rec_next1_H at z = %f\n", z_out);
    strcat(data->error_message, sub_message);
    return;
  }
}


/**********************************************************************************************
Second-order integrator for evolving xe and Tm simultaneously
Used for Hydrogen recombination only
May 2012: added a switch so Peebles model can be used at low redshift.
September 2016: added dEdtdV_dep, the *deposited* energy
(as opposed to dE_dtdV which is the instantaenously injected energy)
July 2020: added implicit integration method
***********************************************************************************************/

void rec_get_xe_next2_HTm(HYREC_DATA *data, int model, double z_in, long iz, double dxedlna_prev[2],
                          double dTmdlna_prev[2], double H, double z_out, double H_next) {

  REC_COSMOPARAMS *cosmo = data->cosmo;
  double xe_in = data->xe_output[iz-1];
  double Tm_in = data->Tm_output[iz-1];
  int *error = &data->error;
  long iz_rad = iz-1-data->rad->iz_rad_0;

  double dxedlna, dTmdlna, nH, TR, dEdtdV, DLNA;
  char sub_message[SIZE_ErrorM];
  DLNA = cosmo->dlna;

  nH = cosmo->nH0 *cube(1.+z_in);
  TR = kBoltz *cosmo->T0 *(1.+z_in);
  dEdtdV = cosmo->inj_params->ion*3*nH*EI / (1-xe_in);
  /*For low redshifts (z < 20 or so) use Peeble's model.
    The precise model does not metter much here as
    1) the free electron fraction is basically zero (~1e-4) in any case and
    2) the universe is going to be reionized around that epoch */

  if (TR/cosmo->fsR/cosmo->fsR/cosmo->meR <= TR_MIN
     || kBoltz*Tm_in/TR <= T_RATIO_MIN) model = PEEBLES;
    dxedlna = rec_dxHIIdlna(data, model, xe_in, xe_in, nH, H, kBoltz*Tm_in, TR, iz_rad, z_in);

  dTmdlna = rec_dTmdlna(z_in, xe_in, Tm_in, cosmo, dEdtdV, H);

  data->Tm_evolve_implicit = 1;
  if (fabs(1-dTmdlna_prev[0]/dTmdlna)<DTM_DIFF_MAX) data->Tm_evolve_implicit = 0;
  if (data->Tm_evolve_implicit == 1){
    dTmdlna_prev[0] = dTmdlna;
    dTmdlna_prev[1] = dTmdlna_prev[0];
    data->xe_output[iz] = xe_in + DLNA *hyrec_integrator(dxedlna, dxedlna_prev, z_in);
    data->Tm_output[iz] = Tm_implicit(z_out, data->xe_output[iz], Tm_in, cosmo, dEdtdV, H_next, DLNA);
  }
  else {
    data->xe_output[iz] = xe_in + DLNA *hyrec_integrator(dxedlna, dxedlna_prev, z_in);
    data->Tm_output[iz] = Tm_in + DLNA *hyrec_integrator(dTmdlna, dTmdlna_prev, z_in);
  }

  if (*error == 1) {
    sprintf(sub_message, "  called from rec_get_xe_next2_HTm at z = %f\n", z_in);
    strcat(data->error_message, sub_message);
    return;
  }
}


/****************************************************************************************************
Builds a recombination history with a given model
(model = PEEBLES, RECFAST, EMLA2s2p, FULL, or SWIFT)
Added May 2012: The radiation field was added as an input so it can be extracted if desired
Replaced condition iz < param->nz by equivalent z >= 0 (so param.nz not needed anymore)
Added September 2016: follow the deposited energy dEdtdV_dep
Added July 2020: Implicit integration for Tm
				 Adopting a smaller time step when diff.eq is stiff in SWIFT mode
****************************************************************************************************/

char* rec_build_history(HYREC_DATA *data, int model, double *hubble_array){

  REC_COSMOPARAMS *cosmo = data->cosmo;
  int *error = &data->error;

  double *xe_output = data->xe_output, *Tm_output = data->Tm_output;
  int Nz = data->Nz;
  double zstart = data->zmax, zend = data->zmin;
  long iz;
  double dxHIIdlna_prev[2], dTmdlna_prev[2], dxHeIIdlna_prev[2];
  double dxHIIdlna_prev_sub[2], dxHeIIdlna_prev_sub[2], xHeII_prev[4];
  double z, dz, DLNA, Delta_xe, xHeII, xH1, dEdtdV_dep, nH, H;
  double *ion = &cosmo->inj_params->ion;
  double *exclya = &cosmo->inj_params->exclya;
  double z_out=0., H_next;
  int flag=10;
  double xe_i, Tm_i;

  DLNA = cosmo->dlna;
  dz = zstart/Nz;
  // Index at which we start integrating Hydrogen recombination, and corresponding redshift
  data->rad->iz_rad_0  = (long) floor(1 + log(kBoltz*cosmo->T0/square(cosmo->fsR)/cosmo->meR*(1.+zstart)/TR_MAX)/DLNA);
  data->rad->z0   = (1.+zstart)*exp(-data->rad->iz_rad_0 * DLNA) - 1.;

  z = zstart;
  /********* He III -> II Saha phase. Tm = Tr. Stop when xHeIII = 1e-8 *********/
  Delta_xe = cosmo->fHe;   /* Delta_xe = xHeIII here */
  for(iz = 0; z >= 0. && Delta_xe > 1e-8; iz++) {
    z = (1.+zstart)*exp(-cosmo->dlna*iz) - 1.;
    xe_output[iz] = rec_xesaha_HeII_III(cosmo, z, &Delta_xe);
    Tm_output[iz] = cosmo->T0 * (1.+z);
  }

  /******** He II -> I recombination.
    Hydrogen in Saha equilibrium.
    Tm fixed to steady state.
    Neglect any energy injection.
    Integrate until TR is low enough that can start integrating hydrogen recombination
    (this occurs at index izH0 computed in rec_get_cosmoparam).
    Start with quasi-equilibrium approximation.
  ********/

  dxHeIIdlna_prev[1] = dxHeIIdlna_prev[0] = 0.;

  xHeII    = rec_saha_xHeII(cosmo, z);
  data->quasi_eq = 1;                    /* Start with post-saha expansion */

  dxHeIIdlna_prev_sub[1] = dxHeIIdlna_prev[1];
  dxHeIIdlna_prev_sub[0] = dxHeIIdlna_prev[0];
  xHeII_prev[3] = xHeII;
  xHeII_prev[2] = xHeII;
  xHeII_prev[1] = xHeII;
  xHeII_prev[0] = xHeII;

  data->loop_after_quasi=1;
  for(; iz <= data->rad->iz_rad_0; iz++) {

    if (model == SWIFT && data->quasi_eq == 0 && data->loop_after_quasi == 1){
      xe_i = xe_output[iz-1]; Tm_i = Tm_output[iz-1];
      for (flag=0;flag<10;flag++) {
        rec_get_xe_next1_He(data, z, &xHeII, dxHeIIdlna_prev_sub, hubble_array, flag);
        z  = (1.+zstart)*exp(-DLNA*(iz-1+(flag+1)/10.)) - 1.;
        xH1 = rec_saha_xH1s(cosmo, z, xHeII);
        xe_i = 1.-xH1 + xHeII;

        if (hubble_array[0]==-1.) H  = rec_HubbleRate(cosmo, z);
        else H = rec_interp1d(.0, dz, hubble_array, Nz, z, error, data->error_message);
        Tm_i = rec_Tmss(z, xe_i, cosmo, 0., H);
      }
      xe_output[iz] = xe_i; Tm_output[iz] = Tm_i;

      xHeII_prev[3] = xHeII_prev[2];
      xHeII_prev[2] = xHeII_prev[1];
      xHeII_prev[1] = xHeII_prev[0];
      xHeII_prev[0] = xHeII;
      dxHeIIdlna_prev[1] = (xHeII_prev[1] - xHeII_prev[3])/2./DLNA;
      dxHeIIdlna_prev[0] = (xHeII_prev[0] - xHeII_prev[2])/2./DLNA;

      if (fabs(1-dxHeIIdlna_prev[1]/dxHeIIdlna_prev[0])<DXHEII_DIFF_MAX) data->loop_after_quasi = 0;
    }
    else{
      rec_get_xe_next1_He(data, z, &xHeII, dxHeIIdlna_prev, hubble_array, flag);
      z             = (1.+zstart)*exp(-DLNA*iz) - 1.;
      xH1           = rec_saha_xH1s(cosmo, z, xHeII);
      xe_output[iz] = 1.-xH1 + xHeII;

      if (hubble_array[0]==-1.) H  = rec_HubbleRate(cosmo, z);
      else H = rec_interp1d(.0, dz, hubble_array, Nz, z, error, data->error_message);

      Tm_output[iz] = rec_Tmss(z, xe_output[iz], cosmo, 0., H);
    }

    if (*error == 1) return data->error_message;
  }

  /******** H II -> I and He II -> I simultaneous recombination (rarely needed but just in case)
        Tm fixed to steady state.
        Integrate H and He simultaneously until xHeII < XHEII_MIN
        Start with post-saha expansion for hydrogen
        Now account for possible energy injection.
        Solve for dEdtdV_dep simultaneously;
  ********/

  dxHIIdlna_prev[1] = (xe_output[iz-2] - xe_output[iz-4])/2./DLNA - dxHeIIdlna_prev[1];
  dxHIIdlna_prev[0] = (xe_output[iz-1] - xe_output[iz-3])/2./DLNA - dxHeIIdlna_prev[0];
  data->quasi_eq          = 1;

  // Initialize energy *deposition*
  dEdtdV_dep = 0.;
  nH = cosmo->nH0*cube(1.+z);

  if (hubble_array[0]==-1.) H  = rec_HubbleRate(cosmo, z);
  else H = rec_interp1d(.0, dz, hubble_array, Nz, z, error, data->error_message);

  update_dEdtdV_dep(z, DLNA, xe_output[iz-1], Tm_output[iz-1], nH, H, cosmo->inj_params, &dEdtdV_dep);
  *ion = dEdtdV_dep/3. /nH *xH1 /EI;
  *exclya = *ion /0.75;

  for(; z >= 0. && xHeII > XHEII_MIN; iz++) {

    get_rec_next2_HHe(data, model, z, iz, &xH1, &xHeII, dxHIIdlna_prev, dxHeIIdlna_prev, H);
    xe_output[iz] = 1.-xH1 + xHeII;

    z  = (1.+zstart)*exp(-DLNA*iz) - 1.;

    if (hubble_array[0]==-1.) H  = rec_HubbleRate(cosmo, z);
    else H = rec_interp1d(.0, dz, hubble_array, Nz, z, error, data->error_message);

    nH = cosmo->nH0*cube(1.+z);
    Tm_output[iz] = rec_Tmss(z, xe_output[iz], cosmo, dEdtdV_dep, H);
    update_dEdtdV_dep(z, DLNA, xe_output[iz], Tm_output[iz], nH, H, cosmo->inj_params, &dEdtdV_dep);
    *ion = dEdtdV_dep/3. /nH *xH1 /EI;
    *exclya = *ion /0.75;

    if (*error == 1) return data->error_message;
  }

  /******** H recombination. Helium assumed entirely neutral.
        Tm fixed to steady-state until its relative difference from Tr is DLNT_MAX
  ********/
  dxHIIdlna_prev_sub[1] = dxHIIdlna_prev[1];
  dxHIIdlna_prev_sub[0] = dxHIIdlna_prev[0];
  data->loop_after_quasi = 1;
  for (; z >= 0. && fabs(1.-Tm_output[iz-1]/cosmo->T0/(1.+z)) < DLNT_MAX; iz++) {

    if (model == SWIFT && data->loop_after_quasi == 1){
      xe_i = xe_output[iz-1]; Tm_i = Tm_output[iz-1];
      for (flag=0;flag<10;flag++) {
        rec_get_xe_next1_H(data, model, z, iz, xe_i, Tm_i, &xe_i, &Tm_i, dxHIIdlna_prev_sub, H, flag);
        z  = (1.+zstart)*exp(-DLNA*(iz-1+(flag+1)/10.)) - 1.;
        nH = cosmo->nH0*cube(1.+z);

        if (hubble_array[0]==-1.) H  = rec_HubbleRate(cosmo, z);
        else H = rec_interp1d(.0, dz, hubble_array, Nz, z, error, data->error_message);

        update_dEdtdV_dep(z, DLNA/10., xe_i, Tm_i, nH, H, cosmo->inj_params, &dEdtdV_dep);
        *ion = dEdtdV_dep/3. /nH *(1.-xe_i) /EI;
        *exclya = *ion /0.75;
      }
      xe_output[iz] = xe_i; Tm_output[iz] = Tm_i;
      dxHIIdlna_prev[1] = (xe_output[iz-1] - xe_output[iz-3])/2./DLNA;
      dxHIIdlna_prev[0] = (xe_output[iz] - xe_output[iz-2])/2./DLNA;

      if (fabs(1-dxHIIdlna_prev[1]/dxHIIdlna_prev[0])<DXHII_DIFF_MAX) data->loop_after_quasi=0;
    }

    else{
      rec_get_xe_next1_H(data, model, z, iz, xe_output[iz-1], Tm_output[iz-1], xe_output+iz, Tm_output+iz,
                         dxHIIdlna_prev, H, flag);
      if (data->quasi_eq == 1){
          dxHIIdlna_prev[1] = (xe_output[iz-2] - xe_output[iz-4])/2./DLNA;
          dxHIIdlna_prev[0] = (xe_output[iz-1] - xe_output[iz-3])/2./DLNA;
      }
      z  = (1.+zstart)*exp(-DLNA*iz) - 1.;
      z_out  = (1.+zstart)*exp(-DLNA*(iz+1)) - 1.;
      nH = cosmo->nH0*cube(1.+z);

      if (hubble_array[0]==-1.) H  = rec_HubbleRate(cosmo, z);
      else H = rec_interp1d(.0, dz, hubble_array, Nz, z, error, data->error_message);

      update_dEdtdV_dep(z, DLNA, xe_output[iz], Tm_output[iz], nH, H, cosmo->inj_params, &dEdtdV_dep);
      *ion = dEdtdV_dep/3. /nH *(1.-xe_output[iz]) /EI;
      *exclya = *ion /0.75;
    }
    if (*error == 1) return data->error_message;
  }

  /******** Evolve xe and Tm simultaneously until z = zend
      Note that the radiative transfer calculation is switched off automatically in the functions
      rec_get_xe_next1_H and rec_get_xe_next2_HTm when it is no longer relevant.
  ********/

  dTmdlna_prev[1] = (Tm_output[iz-2] - Tm_output[iz-4])/2./DLNA;
  dTmdlna_prev[0] = (Tm_output[iz-1] - Tm_output[iz-3])/2./DLNA;

  if (hubble_array[0]==-1.) {
    H  = rec_HubbleRate(cosmo, z);
    H_next  = rec_HubbleRate(cosmo, z_out);
  }
  else {
    H = rec_interp1d(.0, dz, hubble_array, Nz, z, error, data->error_message);
    H_next = rec_interp1d(.0, dz, hubble_array, Nz, z_out, error, data->error_message);
  }

  for(; z > zend; iz++) {
    rec_get_xe_next2_HTm(data, model, z, iz, dxHIIdlna_prev, dTmdlna_prev, H, z_out, H_next);
    z  = (1.+zstart)*exp(-DLNA*iz) - 1.;
    z_out  = (1.+zstart)*exp(-DLNA*(iz+1)) - 1.;
    if (z < zend) z=0.;

    if (hubble_array[0]==-1.) {
      H  = rec_HubbleRate(cosmo, z);
      H_next  = rec_HubbleRate(cosmo, z_out);
    }
    else {
      H = rec_interp1d(.0, dz, hubble_array, Nz, z, error, data->error_message);
      if (z_out > zend) H_next = rec_interp1d(.0, dz, hubble_array, Nz, z_out, error, data->error_message);
      else H_next = H;
    }
    nH = cosmo->nH0*cube(1.+z);
    update_dEdtdV_dep(z, DLNA, xe_output[iz], Tm_output[iz], nH, H, cosmo->inj_params, &dEdtdV_dep);
    *ion = dEdtdV_dep/3. /nH *(1.-xe_output[iz]) /EI;
    *exclya = *ion /0.75;

    if (*error == 1) return data->error_message;
  }
  return data->error_message;
}


/***********************************************************
Function to allocate and initialize HYREC-2 internal tables
Note that path_to_hyrec in HYREC_DATA should be defined first
before calling hyrec_allocate().
***********************************************************/

void hyrec_allocate(HYREC_DATA *data, double zmax, double zmin) {
  double DLNA;
  if (MODEL == SWIFT) DLNA = DLNA_SWIFT;
  else DLNA = DLNA_HYREC;

  data->error = 0;
  data->error_message=malloc(SIZE_ErrorM);
  sprintf(data->error_message, "\n**** ERROR HAS OCCURRED in HYREC-2 ****\n");

  data->zmax = (zmax > 3000.? zmax : 3000.);
  data->zmin = zmin;

  data->atomic = (HYREC_ATOMIC *) malloc(sizeof(HYREC_ATOMIC));
  allocate_and_read_atomic(data->atomic, &data->error, data->path_to_hyrec, data->error_message);

  data->fit = (FIT_FUNC *) malloc(sizeof(FIT_FUNC));
  allocate_and_read_fit(data->fit, &data->error, data->path_to_hyrec, data->error_message);

  data->cosmo  = (REC_COSMOPARAMS *) malloc(sizeof(REC_COSMOPARAMS));
  data->cosmo->inj_params = (INJ_PARAMS *)  malloc(sizeof(INJ_PARAMS));

  data->Nz = (long int) (log((1.+zmax)/(1.+zmin))/DLNA) + 2;
  data->rad = (RADIATION *) malloc(sizeof(RADIATION));

  // For now assume that radiation field never needed over more than 1 decade in redshift
  // (typically from z ~ 1700 to 800 for recombination history)
  // Will have to adapt for outputting radiation fields at lower z
  if (MODEL == FULL)  allocate_radiation(data->rad, (long int) (log(10.)/DLNA), &data->error, data->error_message);

  data->xe_output = create_1D_array(data->Nz, &data->error, data->error_message);
  data->Tm_output = create_1D_array(data->Nz, &data->error, data->error_message);
}


void hyrec_free(HYREC_DATA *data) {
  free_atomic(data->atomic);
  free(data->cosmo->inj_params);
  free(data->cosmo);
  free(data->atomic);
  free(data->xe_output);
  free(data->Tm_output);
  free(data->error_message);
  if (MODEL == FULL) free_radiation(data->rad);
  free(data->rad);
  free_fit(data->fit);
  free(data->fit);
}

/******************************************************************
Compute a recombination history given input cosmological parameters
********************************************************************/

void hyrec_compute(HYREC_DATA *data, int model){
  /* if Hubble_flag[0]=-1., HYREC-2 uses its own Hubble rate.
  When HYREC-2 is included in CLASS, the Hubble rate should be given from CLASS.
  In that case, rec_build_history() would be used directly not this function, hyrec_compute(). */
  double Hubble_flag[1];
  Hubble_flag[0] = -1.;

  rec_build_history(data, model, Hubble_flag);
}

/*****
     Once HYREC_DATA outputs are computed, obtain xe(z) and Tm(z) by interpolation
*****/

double hyrec_xe(double z, HYREC_DATA *data) {
  if (z > data->zmax) return data->xe_output[0];
  if (z < data->zmin) {
    fprintf(stderr, "\033[1m\033[31m error\033[22;30m in hyrec_xe: requesting x_e at z = %f ", z);
    fprintf(stderr, "lower than zmin\n");
    exit(1);
  }
  double DLNA = data->cosmo->dlna;
  return rec_interp1d(-log(1.+data->zmax), DLNA, data->xe_output, data->Nz, -log(1.+z), &data->error, data->error_message);
}

double hyrec_Tm(double z, HYREC_DATA *data) {
  if(z > data->zmax) return data->cosmo->T0*(1.+z);
  if (z < data->zmin) {
    fprintf(stderr, "\033[1m\033[31m error\033[22;30m in hyrec_Tm: requesting x_e at z = %f ", z);
    fprintf(stderr, "lower than zmin\n");
    exit(1);
  }
  double DLNA = data->cosmo->dlna;
  return rec_interp1d(-log(1.+data->zmax), DLNA, data->Tm_output, data->Nz, -log(1.+z), &data->error, data->error_message);
}
