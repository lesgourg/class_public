/*************************************************************************************************/
/*                 HYREC: Hydrogen and Helium Recombination Code                                 */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                              */
/*                                                                                               */
/*         hydrogen.h: all functions related to Hydrogen recombination                           */
/*                                                                                               */
/*         Units used: cgs + eV (all temperatures in eV)                                         */
/*         Version: Januray 2011    (updated value of 2s--1s decay rate,                         */
/*                                  changed temperature range for effective rates)               */
/*************************************************************************************************/
/****** CONSTANTS IN CGS + EV UNIT SYSTEM *******/

#define EI   13.598286071938324              /* Hydrogen ionization energy in eV, reduced mass, no relativistic corrections */

/* Energy differences between excited levels of hydrogen -- used often */
#define E21  10.198714553953742
#define E31  12.087365397278509
#define E41  12.748393192442178
#define E32  1.8886508433247664
#define E42  2.5496786384884356

#define hPc       1.239841874331e-04   /* hc in eV cm */
#define mH        0.93878299831e9      /* Hydrogen atom mass in eV/c^2 */
#define kBoltz    8.617343e-5          /* Boltzmann constant in eV/K */
#define L2s1s     8.2206               /* 2s -> 1s two-photon decay rate in s^{-1} (Labzowsky et al 2005) */



double square(double x);
double cube(double x);

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
   double annihilation_boost_factor;
   double annihilation_m_DM;
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
   int annihil_coef_num_lines;
   double *annihil_coef_heat;
   double *annihil_coef_ionH;
   double *annihil_coef_ionHe;
   double *annihil_coef_lya;
   double *annihil_coef_lowE;
   double *annihil_coef_xe;
   double *annihil_coef_dd_heat;
   double *annihil_coef_dd_ionH;
   double *annihil_coef_dd_ionHe;
   double *annihil_coef_dd_lya;
   double *annihil_coef_dd_lowE;
   double chi_heat;
   double chi_lya;
   double chi_ionH;
   double chi_ionHe;
   double chi_lowE;

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
   double * annihil_z;
   double * annihil_f_eff;
   double * annihil_dd_f_eff;
   short  energy_repart_functions; /**< energy repartition functions */
   double energy_deposition_treatment;
   double f_eff;
   int annihil_f_eff_num_lines;
   int reio_parametrization; /*Do we want the reio by stars based on SFR modeling ? 0 = no, 1 = yes*/
   int star_heating_parametrization; /*Do we want heating by stars based on SFR modeling ? 0 = no, 1 = yes*/

   /** A few parameters if the scheme reio_stars_sfr_source_term is chosen */

   double f_esc; /**< fraction of photons produced by stellar populations that escape to ionize the IGM */
   double Zeta_ion; /**< Lyman continuum photon production efficiency of the stellar population */
   double Log10_Zeta_ion;/**< The log10 of former parameter. */
   double fx; /**< X-ray efficiency fudge factor of photons responsible for heating the medium. */
   double Ex; /**< Associated normalization from Pober et al. 1503.00045. */
   double ap;   /**<  a few parameters entering the fit of the star formation rate (SFR), introduced in Madau & Dickinson, Ann.Rev.Astron.Astrophys. 52 (2014) 415-486, updated in Robertson & al. 1502.02024.*/
   double bp;
   double cp;
   double dp;
   double z_start_reio_stars; /**< Controls the beginning of star reionisation, the SFR experiences is put to 0 above this value. */


   double Omega0_g;
   double odcdmh2,ocdmh2;
   double Omega0_b;
   double Omega0_cdm;
   double Omega0_dcdm;
   double Omega0_lambda;
   double Gamma_dcdm;
   double H0;


} REC_COSMOPARAMS;


/*********** PEEBLES + POST-SAHA + RECFAST ***************/

double alphaB_PPB(double TM);
double rec_HPeebles_dxedlna(double xe, double z, double nH, double H, double TM, double TR, double energy_rate, REC_COSMOPARAMS *param);
double rec_HRecFast_dxedlna(double xe, double z, double nH, double H, double TM, double TR, double energy_rate, REC_COSMOPARAMS *param);

/************* EFFECTIVE MULTI LEVEL ATOM *******************/

#define ALPHA_FILE  "Alpha_inf.dat"     /* Contains the effective recombination coefficients to 2s and 2p */
#define RR_FILE     "R_inf.dat"         /* Contains the effective transfer rate R_{2p,2s} */


/* Boundaries and number of elements of temperature tables */
#define TR_MIN 0.004            /* Tr parameters */
#define TR_MAX 0.4
#define NTR    100
#define TM_TR_MIN 0.1           /* Tm/Tr parameters */
#define TM_TR_MAX 1.0
#define NTM 40

/* Effective rate coefficients structure */
typedef struct {
  double *logTR_tab;
  double *TM_TR_tab;
  double **logAlpha_tab[2];
  double *logR2p2s_tab;
  double DlogTR, DTM_TR;
}
HRATEEFF;

void read_rates(HRATEEFF *rate_table);
void interpolate_rates(double Alpha[2], double Beta[2], double *R2p2s, double TR, double TM_TR, HRATEEFF *rate_table);
double rec_HMLA_dxedlna(double xe, double z, double nH, double Hubble, double TM, double TR, double energy_rate, REC_COSMOPARAMS *param, HRATEEFF *rate_table);

/************ TWO-PHOTON PROCESSES AND DIFFUSION  ************/

#define EFFECT_A    1    /* 2s-->1s stimulated two-photon decays and non-thermal absorptions */
#define EFFECT_B    1    /* Sub-Lyman alpha two-photon transitions 3s/3d<--> 1s and 4s/4d<-->1s */
#define EFFECT_C    1    /* Super-Lyman alpha two-photon transitions 3s/3d<--> 1s and 4s/4d<-->1s */
#define EFFECT_D    1    /* Raman scattering from 2s and 3s/3d */
#define DIFFUSION   1    /* Lyman alpha frequency diffusion */


#define NSUBLYA  140
#define NSUBLYB  271
#define NVIRT    311
#define NDIFF    80
#define TWOG_FILE "two_photon_tables.dat"
/* maximum dlna = 8.5e-5 */



typedef struct {
    double Eb_tab[NVIRT];       /* Energies of the virtual levels in eV */
    double A1s_tab[NVIRT];      /* 3*A2p1s*phi(E)*DE */
    double A2s_tab[NVIRT];      /* dLambda_2s/dE * DeltaE if E < Elya dK2s/dE * Delta E if E > Elya */
    double A3s3d_tab[NVIRT];    /* (dLambda_3s/dE + 5*dLambda_3d/dE) * Delta E for E < ELyb, Raman scattering rate for E > ELyb */
    double A4s4d_tab[NVIRT];    /* (dLambda_4s/dE + 5*dLambda_4d/dE) * Delta E */
}  TWO_PHOTON_PARAMS;


void read_twog_params(TWO_PHOTON_PARAMS *twog);
void populate_Diffusion(double *Aup, double *Adn, double *A2p_up, double *A2p_dn,
                        double TM, double Eb_tab[NVIRT], double A1s_tab[NVIRT]);
void populateTS_2photon(double Trr[2][2], double *Trv[2], double *Tvr[2], double *Tvv[3],
                        double sr[2], double sv[NVIRT], double Dtau[NVIRT],
                        double xe, double TM, double TR, double nH, double H, HRATEEFF *rate_table,
                        TWO_PHOTON_PARAMS *twog, double fplus[NVIRT], double fplus_Ly[],
                        double Alpha[], double Beta[], double z);
void solveTXeqB(double *diag, double *updiag, double *dndiag, double *X, double *B, unsigned N);
void solve_real_virt(double xr[2], double xv[NVIRT], double Trr[2][2], double *Trv[2], double *Tvr[2],
                     double *Tvv[3], double sr[2], double sv[NVIRT]);
void fplus_from_fminus(double fplus[NVIRT], double fplus_Ly[], double **logfminus_hist, double *logfminus_Ly_hist[],
                       double TR, double zstart, double dlna, unsigned iz, double z, double Eb_tab[NVIRT]);
double rec_HMLA_2photon_dxedlna(double xe, double nH, double H, double TM, double TR,
                                HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog,
                                double zstart, double dlna, double **logfminus_hist, double *logfminus_Ly_hist[], unsigned iz, double z,
								double energy_rate, REC_COSMOPARAMS *param);
double xe_PostSahaH(double nH, double H, double T, HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog,
                    double zstart, double dlna, double **logfminus_hist, double *logfminus_Ly_hist[],
                    unsigned iz, double z, double *Dxe, int model, double energy_rate, REC_COSMOPARAMS *param);
void update_fminus_Saha(double **logfminus_hist, double *logfminus_Ly_hist[],
                        double xe, double TR, double nH, TWO_PHOTON_PARAMS *twog,
			double zstart, double dlna, unsigned iz, double z, int func_select);

int hyrec_annihilation_coefficients_interpolate(REC_COSMOPARAMS *param,
                                                        double xe,
                                                        double * chi_heat,
                                                        double * chi_lya,
                                                        double * chi_ionH,
                                                        double * chi_ionHe,
                                                        double * chi_lowE
                                                      );
