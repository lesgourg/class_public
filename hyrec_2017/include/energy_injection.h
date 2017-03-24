/* Structure with all energy injection parameters */
/* If adding a new energy injection process 
   make sure to add relevant parameters here */

typedef struct {

  double odmh2;                 /* Omega_dm h^2 */
  
  double pann, pann_halo;       /* DM annihilation parameter in the smooth background and in haloes */
                                /* Units of pann and pann_halo are cm^3/s/GeV */

  double ann_z, ann_zmax, ann_zmin, ann_var; /* Parameters for the variation of pann(z) */
  double ann_z_halo;                         /* Characteristic redshift for annihilation in haloes */
    
  double Mpbh, fpbh;           /* Mass and fraction of DM made of primordial black holes */
  int coll_ion;                /* If 1: assume gas gest collisionally ionized. 
                                  If 0: assume gas gets photoionized by PBH radiation */     

  int on_the_spot;            /* if set to 1 assume energy deposition rate = injection rate */
                              /* Otherwise solves for deposition given injection with simple recipe */
  
} INJ_PARAMS;

double beta_pbh(double Mpbh, double z, double xe, double Teff);
double gamma_pbh(double Mpbh, double z, double xe, double Teff);
double lambda_pbh(double Mpbh, double z, double xe, double Teff);
double Mdot_pbh(double Mpbh, double z, double xe, double Teff);
double TS_over_me_pbh(double Mpbh, double z, double xe, double Teff, int coll_ion);
double eps_over_mdot_pbh(double Mpbh, double z, double xe, double Teff, int coll_ion);
double L_pbh(double Mpbh, double z, double xe, double Teff, int coll_ion);
double vbc_rms_func(double z);
double L_pbh_av(double Mpbh, double z, double xe, double Tgas, int coll_ion);
double dEdtdV_pbh(double fpbh, double Mpbh, double z, double xe, double Tgas, int coll_ion);
double dEdtdV_inj(double z, double xe, double Tgas, INJ_PARAMS *params);
void update_dEdtdV_dep(double z_out, double dlna, double xe, double Tgas,
		       double nH, double H, INJ_PARAMS *params, double *dEdtdV_dep);
