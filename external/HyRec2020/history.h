/*************************************************************************************************/
/*                 HYREC-2: Hydrogen and Helium Recombination Code                               */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (2010-17)                              */
/*            with contributions from Nanoom Lee (2020)                                          */
/*                                                                                               */
/*         history.h: functions for recombination history                                        */
/*                                                                                               */
/*************************************************************************************************/


#include "hyrectools.h"
#include "hydrogen.h"
#include "helium.h"
#include "energy_injection.h"

#define HYREC_VERSION "2020"

#define MODEL SWIFT					/* SWIFT is the default model. Four more models can be used (PEEBLES, RECFAST, EMLA2s2p, FULL). */
									/* Each model is defined in hydrogen.h */

#define SIZE_ErrorM      2048

#define DXHEII_MAX       1e-5       /* If xHeII - xHeII(Saha) < DXEHII_MAX, use post-Saha expansion for Helium. Lower value = higher accuracy. */

#define DXHII_MAX_FULL   3e-4      
#define DXHII_MAX        3e-4       /* If xHII - xHII(Saha) < DXHII_MAX, use post-Saha expansion for Hydrogen. Switch to ODE integration after that.
                                    IMPORTANT: do not set to a lower value unless using a smaller time-step */

#define XHEII_MIN        1e-6       /* Stop considering Helium recombination once xHeII < XHEII_MIN */
//#define XHEII_MIN      1e-10      /* Used when calculating correction function in SWIFT mode */ 

#define DLNT_MAX         5e-4       /* Use the steady-state approximation for Tm as long as 1-Tm/Tr < DLNT_MAX, then switch to ODE integration */


/* Structure for HyRec internal parameters */ 

typedef struct {
  double h;                                 /* Hubble constant */
  double T0;                                /* CMB temperature today in K*/
  double obh2, ocbh2, odeh2, okh2, orh2, onuh2;     /* density parameters */ 
  double Nmnu, Nnueff;                             /* effective number of neutrinos */
  double mnu[3];                                   /* neutrino masses */
  double fHe;                               /* Helium fraction by number */
  double nH0;                               /* density of hydrogen today in cm^{-3} [Changed from m^{-3} in February 2015] */ 
  double YHe;                               /* Helium fraction */ 
  double fsR, meR;              /* fine-structure constant alpha/alpha(today) 
                                    and me/me(today) (Added April 2012)*/
  double zstart, zend, dlna, nz; 

  INJ_PARAMS *inj_params;     /* Structure containing all Energy-injection parameters */

} REC_COSMOPARAMS;

void rec_get_cosmoparam(FILE *fin, FILE *fout, REC_COSMOPARAMS *param, int *error, char error_message[SIZE_ErrorM]);

double rec_HubbleRate(REC_COSMOPARAMS *param, double z, int *error, char error_message[SIZE_ErrorM]);

double rec_Tmss(double z, double xe, REC_COSMOPARAMS *cosmo, double dEdtdV, double H);

double rec_dTmdlna(double z, double xe, double Tm, REC_COSMOPARAMS *cosmo, double dEdtdV, double H);

double Tm_implicit(double z, double xe, double Tm, REC_COSMOPARAMS *cosmo, double dEdtdV, double H, double DLNA);

void rec_get_xe_next1_He(REC_COSMOPARAMS *param, double z_in, double *xHeII,
			 double *dxHeIIdlna_prev, int *post_saha, double *hubble_array, int Nz, double dz, int *error, char error_message[SIZE_ErrorM], int flag);

void rec_xH1_stiff(int model, REC_COSMOPARAMS *param, double z, double xHeII, double *xH1, 
		   HYREC_ATOMIC *atomic, RADIATION *rad, FIT_FUNC *fit, unsigned iz_rad,
		   double ion, double exclya, int *stiff, int *error, char error_message[SIZE_ErrorM], double H);

void get_rec_next2_HHe(int model, REC_COSMOPARAMS *param, double z_in, double Tm,
                       double *xH1, double *xHeII, HYREC_ATOMIC *atomic, RADIATION *rad, FIT_FUNC *fit, unsigned iz_rad,
		       double dxHIIdlna_prev[2], double dxHeIIdlna_prev[2], double ion, double exclya, int *stiff, int *error, char error_message[SIZE_ErrorM], double H);

void rec_get_xe_next1_H(int model ,REC_COSMOPARAMS *param, double z_in, double xe_in, double Tm_in,
			double *xe_out, double *Tm_out, HYREC_ATOMIC *atomic, RADIATION *rad, FIT_FUNC *fit, unsigned iz_rad,
			double dxedlna_prev[2], double ion, double exclya, int *stiff, int *error, char error_message[SIZE_ErrorM], double H, int flag);

void rec_get_xe_next2_HTm(int model, REC_COSMOPARAMS *param,
			  double z_in, double xe_in, double Tm_in, double *xe_out, double *Tm_out,
                          HYREC_ATOMIC *atomic, RADIATION *rad, FIT_FUNC *fit, unsigned iz_rad,
			  double dxedlna_prev[2], double dTmdlna_prev[2], double ion, double exclya, int *error, char error_message[SIZE_ErrorM], double H, double z_out, double H_next);

char* rec_build_history(int model, double zstart, double zend,
		       REC_COSMOPARAMS *param, HYREC_ATOMIC *atomic, RADIATION *rad, FIT_FUNC *fit,
		       double *xe_output, double *Tm_output, double *hubble_array, int Nz, int *error, char error_message[]);


typedef struct{
  HYREC_ATOMIC *atomic;
  REC_COSMOPARAMS *cosmo;
  double zmax;
  double zmin;
  long int Nz;
  double *xe_output;
  double *Tm_output;
  int error;
  char *error_message;
  char *path_to_hyrec;
  RADIATION *rad; 
  FIT_FUNC *fit;
} HYREC_DATA;


void hyrec_allocate(HYREC_DATA *data, double zmax, double zmin);
void hyrec_free(HYREC_DATA *data);
void hyrec_compute(HYREC_DATA *data, int model);
double hyrec_xe(double z, HYREC_DATA *data);
double hyrec_Tm(double z, HYREC_DATA *data);
