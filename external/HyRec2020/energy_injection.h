/* Structure with all energy injection parameters */
/* If adding a new energy injection process 
   make sure to add relevant parameters here */

#ifndef __ENERGY_INJECTION__
#define __ENERGY_INJECTION__

typedef struct {

  double odmh2;                 /* Omega_dm h^2 */
  
  double pann, pann_halo;       /* DM annihilation parameter in the smooth background and in haloes */
                                /* Units of pann and pann_halo are cm^3/s/GeV */

  double ann_z, ann_zmax, ann_zmin, ann_var; /* Parameters for the variation of pann(z) */
  double ann_z_halo;                         /* Characteristic redshift for annihilation in haloes */
  
  double decay;

  double Mpbh, fpbh;           /* Mass and fraction of DM made of primordial black holes */

  int on_the_spot;            /* if set to 1 assume energy deposition rate = injection rate */
                              /* Otherwise solves for deposition given injection with simple recipe */

  double ion, exclya;
  
} INJ_PARAMS;

void update_dEdtdV_dep(double z_out, double dlna, double xe, double Tgas,
                       double nH, double H, INJ_PARAMS *params, double *dEdtdV_dep);

#endif
