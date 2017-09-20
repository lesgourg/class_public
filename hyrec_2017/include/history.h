/*************************************************************************************************/
/*                 HYREC: Hydrogen and Helium Recombination Code                                 */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                              */
/*                                                                                               */
/*         history.h: functions for recombination history                                        */
/*                                                                                               */
/*************************************************************************************************/

#include "hyrectools.h"
#include "hydrogen.h"
#include "helium.h"



double rec_HubbleRate(REC_COSMOPARAMS *param, double z);

double rec_Tmss(double z, double xe, REC_COSMOPARAMS *cosmo, double dEdtdV);

double rec_dTmdlna(double z, double xe, double Tm, REC_COSMOPARAMS *cosmo, double dEdtdV);

void rec_get_xe_next1_He(REC_COSMOPARAMS *param, double z_in, double *xHeII,
			 double *dxHeIIdlna_prev, int *post_saha);

void rec_xH1_stiff(int model, REC_COSMOPARAMS *param, double z, double xHeII, double *xH1,
		   HYREC_ATOMIC *atomic, RADIATION *rad, unsigned iz_rad,
		   double dEdtdV, int *stiff, int *error);

void get_rec_next2_HHe(int model, REC_COSMOPARAMS *param, double z_in, double Tm,
                       double *xH1, double *xHeII, HYREC_ATOMIC *atomic, RADIATION *rad, unsigned iz_rad,
		       double dxHIIdlna_prev[2], double dxHeIIdlna_prev[2], double dEdtdV, int *stiff, int *error);

void rec_get_xe_next1_H(int model, REC_COSMOPARAMS *param, double z_in, double xe_in, double Tm_in,
			double *xe_out, double *Tm_out, HYREC_ATOMIC *atomic, RADIATION *rad, unsigned iz_rad,
			double dxedlna_prev[2], double dEdtdV, int *stiff, int *error);

void rec_get_xe_next2_HTm(int model, REC_COSMOPARAMS *param,
			  double z_in, double xe_in, double Tm_in, double *xe_out, double *Tm_out,
                          HYREC_ATOMIC *atomic, RADIATION *rad, unsigned iz_rad,
			  double dxedlna_prev[2], double dTmdlna_prev[2], double dEdtdV, int *error);

void rec_build_history(int model, double zstart, double zend,
		       REC_COSMOPARAMS *param, HYREC_ATOMIC *atomic, RADIATION *rad,
		       double *xe_output, double *Tm_output);


typedef struct{
  HYREC_ATOMIC *atomic;
  REC_COSMOPARAMS *cosmo;
  double zmax;
  double zmin;
  long int Nz;
  double *xe_output;
  double *Tm_output;
  RADIATION *rad;
} HYREC_DATA;
