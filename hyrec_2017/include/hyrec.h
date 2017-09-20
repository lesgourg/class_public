#define HYREC_VERSION "2017"

#include "history.h"


void hyrec_allocate(HYREC_DATA *data, double zmax, double zmin);

void hyrec_free(HYREC_DATA *data);

void hyrec_compute_CLASS(HYREC_DATA *data, int model);

void hyrec_compute(HYREC_DATA *data, int model,
		   double h, double T0, double Omega_b, double Omega_m, double Omega_k, double YHe, double Nnueff,
		   double alphaR, double meR, double pann, double pann_halo, double ann_z, double ann_zmax,
		   double ann_zmin, double ann_var, double ann_z_halo, double Mpbh, double fpbh, int feedback_pbh, int on_the_spot);

double hyrec_xe(double z, HYREC_DATA *rec_data);

double hyrec_Tm(double z, HYREC_DATA *rec_data);

double hyrec_dTmdlna(double z, HYREC_DATA *rec_data);
