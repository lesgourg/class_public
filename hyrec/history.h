/*************************************************************************************************/
/*                 HYREC: Hydrogen and Helium Recombination Code                                 */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                              */
/*                                                                                               */
/*         history.h: functions for recombination history                                        */
/*                                                                                               */
/*************************************************************************************************/

#define PROMPT 1      /* Set to zero to suppress initial prompts */

/**** Switch to choose the physical model used for hydrogen ****/ 

/* definitions*/
#define PEEBLES   0    /* Peebles effective three-level atom */
#define RECFAST   1    /* Effective three-level atom for hydrogen with fudge factor F = 1.14 */
#define EMLA2s2p  2    /* Correct EMLA model, with standard decay rates from 2s and 2p only */
#define FULL      3    /* All radiative transfer effects included. Additional switches in header file hydrogen.h */

/** here is the switch **/
#define MODEL RECFAST     /* default setting: FULL */

/***** Switches for derivative d(xe)/dt *****/

#define FUNC_HEI     1
#define FUNC_H2G     2
#define FUNC_HMLA    3
#define FUNC_PEEBLES 4

/**** Cosmological parameters. 
      Include information on starting and ending redshit and timestep  ****/

typedef struct {
   double T0;                   /* CMB temperature today in K*/
   double obh2, omh2, okh2;     /* cosmological parameters */
   double odeh2, w0, wa;        /* dark energy parameters */
   double Y;                    /* primordial helium abundance */
   double Nnueff;               /* effective number of neutrinos */

   /* Secondary parameters, to avoid recalculating every time */
   double nH0;                  /* density of hydrogen today in m^{-3} */  
   double fHe;                  /* Helium fraction by number */

   double zstart, zend, dlna;   /* initial and final redshift and step size in log a */
   long nz;                     /* total number of redshift steps */

   /** parameters for energy injection */

   double annihilation; /** parameter describing CDM annihilation (f <sigma*v> / m_cdm, see e.g. 0905.0003) */
  
   short has_on_the_spot; /** do we want to use the on-the-spot approximation? */

   double decay; /** parameter descibing CDM decay (f/tau, see e.g. 1109.6322)*/

   double annihilation_variation; /** if this parameter is non-zero,
				     the function F(z)=(f <sigma*v> /
				     m_cdm)(z) will be a parabola in
				     log-log scale between zmin and
				     zmax, with a curvature given by
				     annihlation_variation (must ne
				     negative), and with a maximum in
				     zmax; it will be constant outside
				     this range */

   double annihilation_z; /** if annihilation_variation is non-zero,
			     this is the value of z at which the
			     parameter annihilation is defined, i.e.
			     F(annihilation_z)=annihilation */
  
   double annihilation_zmax; /** if annihilation_variation is non-zero,
				redhsift above which annihilation rate
				is maximal */

   double annihilation_zmin; /** if annihilation_variation is non-zero,
				redhsift below which annihilation rate
				is constant */

   double annihilation_f_halo; /* takes the contribution of DM annihilation in halos into account*/
   double annihilation_z_halo; /*characteristic redshift for DM annihilation in halos*/

} REC_COSMOPARAMS;

void rec_get_cosmoparam(FILE *fin, FILE *fout, REC_COSMOPARAMS *param);
double rec_HubbleConstant(REC_COSMOPARAMS *param, double z);
double rec_Tmss(double xe, double Tr, double H, double fHe, double nH, double energy_rate);
double rec_dTmdlna(double xe, double Tm, double Tr, double H, double fHe , double nH, double energy_rate);
void rec_get_xe_next1(REC_COSMOPARAMS *param, double z1, double xe_in, double *xe_out,
                      HRATEEFF *rate_table, int func_select, unsigned iz, TWO_PHOTON_PARAMS *twog_params,
		      double **logfminus_hist, double *logfminus_Ly_hist[], 
                      double *z_prev, double *dxedlna_prev, double *z_prev2, double *dxedlna_prev2);
void rec_get_xe_next2(REC_COSMOPARAMS *param, double z1, double xe_in, double Tm_in, double *xe_out, double *Tm_out,
                      HRATEEFF *rate_table, int func_select, unsigned iz, TWO_PHOTON_PARAMS *twog_params,
		      double **logfminus_hist, double *logfminus_Ly_hist[], 
                      double *z_prev, double *dxedlna_prev, double *dTmdlna_prev, 
                      double *z_prev2, double *dxedlna_prev2, double *dTmdlna_prev2);
void rec_build_history(REC_COSMOPARAMS *param, HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog_params,
                       double *xe_output, double *Tm_output);

double energy_injection_rate(REC_COSMOPARAMS *param, double z);
