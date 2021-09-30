/*************************************************************************************************/
/*                 HYREC-2: Hydrogen and Helium Recombination Code                               */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (2010-17)                              */
/*            with contributions from Nanoom Lee (2020)                                          */
/*                                                                                               */
/*         history.h: functions for recombination history                                        */
/*                                                                                               */
/*************************************************************************************************/

#ifndef __HISTORY__
#define __HISTORY__

#include "hyrectools.h"
#include "hydrogen.h"

#define HYREC_VERSION "2020"

#define MODEL SWIFT	                /* SWIFT is the default model. Four more models can be used (PEEBLES, RECFAST, EMLA2s2p, FULL). */
                                    /* Each model is defined in hydrogen.h */

/* !!!!!  Do NOT change any numbers below unless you know what's going on with each parameter exactly !!!!! */

#define SIZE_ErrorM      2048

#define DXHEII_MAX       1e-5       /* If xHeII - xHeII(Saha) < DXEHII_MAX, use post-Saha expansion for Helium.*/
#define DXHEII_DIFF_MAX  5e-2       /* If |1-dxHeIIdlna_prev/dxHeIIdlna| > DXHEII_DIFF_MAX, do loop with 10 times smaller time step */

#define DXHII_MAX        3e-4       /* If xHII - xHII(Saha) < DXHII_MAX, use post-Saha expansion for Hydrogen. Switch to ODE integration after that.
                                    IMPORTANT: do not set to a lower value unless using a smaller time-step */
#define DXHII_DIFF_MAX   5e-2       /* If |1-dxHIIdlna_prev/dxHIIdlna| > DXHII_DIFF_MAX, do loop with 10 times smaller time step */

#define XHEII_MIN        1e-6       /* Stop considering Helium recombination once xHeII < XHEII_MIN */
//#define XHEII_MIN      1e-10      /* Used when calculating correction function in SWIFT mode */

#define DLNT_MAX         5e-4       /* Use the steady-state approximation for Tm as long as 1-Tm/Tr < DLNT_MAX, then switch to ODE integration */
#define DTM_DIFF_MAX     5e-2       /* If |1-dTmdlna_prev/dTmdlna| > DTM_DIFF_MAX, evole Tm with implicit method */

void rec_get_cosmoparam(FILE *fin, FILE *fout, REC_COSMOPARAMS *param);

double rec_HubbleRate(REC_COSMOPARAMS *param, double z);

double rec_Tmss(double z, double xe, REC_COSMOPARAMS *cosmo, double dEdtdV, double H);

double rec_dTmdlna(double z, double xe, double Tm, REC_COSMOPARAMS *cosmo, double dEdtdV, double H);

double Tm_implicit(double z, double xe, double Tm, REC_COSMOPARAMS *cosmo, double dEdtdV, double H, double DLNA);

void rec_get_xe_next1_He(HYREC_DATA *data, double z_in, double *xHeII, double dxHeIIdlna_prev[2],
                         double *hubble_array, int flag);

void rec_xH1_stiff(HYREC_DATA *data, int model, double z, double xHeII, double *xH1, unsigned iz_rad, double H);

void get_rec_next2_HHe(HYREC_DATA *data, int model, double z_in, long iz, double *xH1, double *xHeII,
                       double dxHIIdlna_prev[2], double dxHeIIdlna_prev[2], double H);

void rec_get_xe_next1_H(HYREC_DATA *data, int model, double z_in, long iz, double xe_in, double Tm_in,
                        double *xe_out, double *Tm_out, double dxedlna_prev[2], double H, int flag);

void rec_get_xe_next2_HTm(HYREC_DATA *data, int model, double z_in, long iz, double dxedlna_prev[2],
                          double dTmdlna_prev[2], double H, double z_out, double H_next);

char* rec_build_history(HYREC_DATA *data, int model, double *hubble_array);

void hyrec_allocate(HYREC_DATA *data, double zmax, double zmin);
void hyrec_free(HYREC_DATA *data);
void hyrec_compute(HYREC_DATA *data, int model);
double hyrec_xe(double z, HYREC_DATA *data);
double hyrec_Tm(double z, HYREC_DATA *data);

#endif
