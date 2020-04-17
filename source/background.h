/** @file background.h Documented includes for background module */

#ifndef __BACKGROUND__
#define __BACKGROUND__

#include "common.h"
#include "quadrature.h"
#include "growTable.h"
#include "arrays.h"
#include "dei_rkck.h"
#include "parser.h"

//The name for this macro can be at most 30 characters total
#define _class_print_species_(name,type) \
printf("-> %-30s Omega = %-15g , omega = %-15g\n",name,pba->Omega0_##type,pba->Omega0_##type*pba->h*pba->h);

/** list of possible types of spatial curvature */

enum spatial_curvature {flat,open,closed};

/** list of possible parametrisations of the DE equation of state */

enum equation_of_state {CLP,EDE};

/**
 * All background parameters and evolution that other modules need to know.
 *
 * Once initialized by the backgound_init(), contains all necessary
 * information on the background evolution (except thermodynamics),
 * and in particular, a table of all background quantities as a
 * function of time and scale factor, used for interpolation in other
 * modules.
 */

struct background
{
  /** @name - input parameters initialized by user in input module
   *  (all other quantities are computed in this module, given these parameters
   *   and the content of the 'precision' structure)
   *
   * The background cosmological parameters listed here form a parameter
   * basis which is directly usable by the background module. Nothing
   * prevents from defining the input cosmological parameters
   * differently, and to pre-process them into this format, using the input
   * module (this might require iterative calls of background_init()
   * e.g. for dark energy or decaying dark matter). */

  //@{

  double H0; /**< \f$ H_0 \f$: Hubble parameter (in fact, [\f$H_0/c\f$]) in \f$ Mpc^{-1} \f$ */

  double Omega0_g; /**< \f$ \Omega_{0 \gamma} \f$: photons */

  double T_cmb; /**< \f$ T_{cmb} \f$: current CMB temperature in Kelvins */

  double Omega0_b; /**< \f$ \Omega_{0 b} \f$: baryons */

  double Omega0_cdm; /**< \f$ \Omega_{0 cdm} \f$: cold dark matter */

  double Omega0_lambda; /**< \f$ \Omega_{0_\Lambda} \f$: cosmological constant */

  double Omega0_fld; /**< \f$ \Omega_{0 de} \f$: fluid */

  enum equation_of_state fluid_equation_of_state; /**< parametrisation scheme for fluid equation of state */

  double w0_fld; /**< \f$ w0_{DE} \f$: current fluid equation of state parameter */
  double wa_fld; /**< \f$ wa_{DE} \f$: fluid equation of state parameter derivative */
  double Omega_EDE; /**< \f$ wa_{DE} \f$: Early Dark Energy density parameter */

  double cs2_fld; /**< \f$ c^2_{s~DE} \f$: sound speed of the fluid
		     in the frame comoving with the fluid (so, this is
		     not [delta p/delta rho] in the synchronous or
		     newtonian gauge!) */

  short use_ppf; /**< flag switching on PPF perturbation equations
                    instead of true fluid equations for
                    perturbations. It could have been defined inside
                    perturbation structure, but we leave it here in
                    such way to have all fld parameters grouped. */

  double c_gamma_over_c_fld; /**< ppf parameter defined in eq. (16) of 0808.3125 [astro-ph] */

  double Omega0_ur; /**< \f$ \Omega_{0 \nu r} \f$: ultra-relativistic neutrinos */

  double Omega0_idr; /**< \f$ \Omega_{0 idr} \f$: interacting dark radiation */
  double T_idr;      /**< \f$ T_{idr} \f$: current temperature of interacting dark radiation in Kelvins */

  double Omega0_idm_dr; /**< \f$ \Omega_{0 idm_dr} \f$: dark matter interacting with dark radiation */

  double Omega0_dcdmdr; /**< \f$ \Omega_{0 dcdm}+\Omega_{0 dr} \f$: decaying cold dark matter (dcdm) decaying to dark radiation (dr) */

  double Gamma_dcdm; /**< \f$ \Gamma_{dcdm} \f$: decay constant for decaying cold dark matter */

  double Omega_ini_dcdm;    /**< \f$ \Omega_{ini,dcdm} \f$: rescaled initial value for dcdm density (see 1407.2418 for definitions) */

  double Omega0_scf;        /**< \f$ \Omega_{0 scf} \f$: scalar field */
  short attractor_ic_scf;   /**< whether the scalar field has attractor initial conditions */
  double phi_ini_scf;       /**< \f$ \phi(t_0) \f$: scalar field initial value */
  double phi_prime_ini_scf; /**< \f$ d\phi(t_0)/d\tau \f$: scalar field initial derivative wrt conformal time */
  double * scf_parameters;  /**< list of parameters describing the scalar field potential */
  int scf_parameters_size;  /**< size of scf_parameters */
  int scf_tuning_index;     /**< index in scf_parameters used for tuning */
  //double scf_lambda; /**< \f$ \lambda \f$ : scalar field exponential potential slope */
  //double scf_alpha;  /**< \f$ \alpha \f$ : Albrecht-Skordis polynomial slope */
  //double scf_B; /**< \f$ \alpha \f$ : Albrecht-Skordis field shift */
  //double scf_A; /**< \f$ \alpha \f$ : Albrecht-Skordis offset */

  double Omega0_k; /**< \f$ \Omega_{0_k} \f$: curvature contribution */

  int N_ncdm;                            /**< Number of distinguishable ncdm species */
  double * M_ncdm;                       /**< vector of masses of non-cold relic:
                                             dimensionless ratios m_ncdm/T_ncdm */
  double * Omega0_ncdm, Omega0_ncdm_tot; /**< Omega0_ncdm for each species and for the total Omega0_ncdm */
  double * deg_ncdm, deg_ncdm_default;   /**< vector of degeneracy parameters in factor
                                             of p-s-d: 1 for one family of neutrinos
                                             (= one neutrino plus its anti-neutrino,
                                             total g*=1+1=2, so deg = 0.5 g*); and its
					     default value */

  /* the following parameters help to define the analytical ncdm phase space distributions (p-s-d) */
  double * T_ncdm,T_ncdm_default;       /**< list of 1st parameters in
					     p-s-d of non-cold relics:
					     relative temperature
					     T_ncdm1/T_gamma; and its
					     default value */
  double * ksi_ncdm, ksi_ncdm_default;  /**< list of 2nd parameters in
					     p-s-d of non-cold relics:
					     relative chemical potential
					     ksi_ncdm1/T_ncdm1; and its
					     default value */
  double * ncdm_psd_parameters;         /**< list of parameters for specifying/modifying
                                             ncdm p.s.d.'s, to be customized for given model
                                             (could be e.g. mixing angles) */
  /* end of parameters for analytical ncdm p-s-d */

  /* the following parameters help to define tabulated ncdm p-s-d passed in file */
  int * got_files;                      /**< list of flags for each species, set to true if
					     p-s-d is passed through file */
  char * ncdm_psd_files;                /**< list of filenames for tabulated p-s-d */
  /* end of parameters for tabulated ncdm p-s-d */

  //@}

  /** @name - related parameters */

  //@{

  double h; /**< reduced Hubble parameter */
  double K; /**< \f$ K \f$: Curvature parameter \f$ K=-\Omega0_k*a_{today}^2*H_0^2\f$; */
  int sgnK; /**< K/|K|: -1, 0 or 1 */
  double * m_ncdm_in_eV; /**< list of ncdm masses in eV (inferred from M_ncdm and other parameters above) */

  //@}

  /** @name - other background parameters */

  //@{

  double a_today; /**< scale factor today (arbitrary and irrelevant for most purposes) */

  //@}

  /** @name - flags describing the absence or presence of cosmological
      ingredients
      *
      * having one of these flag set to zero allows to skip the
      * corresponding contributions, instead of adding null contributions.
      */


  //@{

  short has_cdm;       /**< presence of cold dark matter? */
  short has_dcdm;      /**< presence of decaying cold dark matter? */
  short has_dr;        /**< presence of relativistic decay radiation? */
  short has_scf;       /**< presence of a scalar field? */
  short has_ncdm;      /**< presence of non-cold dark matter? */
  short has_lambda;    /**< presence of cosmological constant? */
  short has_fld;       /**< presence of fluid with constant w and cs2? */
  short has_ur;        /**< presence of ultra-relativistic neutrinos/relics? */
  short has_idr;       /**< presence of interacting dark radiation? */
  short has_idm_dr;    /**< presence of dark matter interacting with dark radiation? */
  short has_curvature; /**< presence of global spatial curvature? */

  //@}

  /**
   *@name - arrays related to sampling and integration of ncdm phase space distributions
   */


  //@{
  int * ncdm_quadrature_strategy; /**< Vector of integers according to quadrature strategy. */
  int * ncdm_input_q_size; /**< Vector of numbers of q bins */
  double * ncdm_qmax;   /**< Vector of maximum value of q */
  double ** q_ncdm_bg;  /**< Pointers to vectors of background sampling in q */
  double ** w_ncdm_bg;  /**< Pointers to vectors of corresponding quadrature weights w */
  double ** q_ncdm;     /**< Pointers to vectors of perturbation sampling in q */
  double ** w_ncdm;     /**< Pointers to vectors of corresponding quadrature weights w */
  double ** dlnf0_dlnq_ncdm; /**< Pointers to vectors of logarithmic derivatives of p-s-d */
  int * q_size_ncdm_bg; /**< Size of the q_ncdm_bg arrays */
  int * q_size_ncdm;    /**< Size of the q_ncdm arrays */
  double * factor_ncdm; /**< List of normalization factors for calculating energy density etc.*/

  //@}

  /**
   *@name - some flags needed for calling background functions
   */

  //@{

  short short_info;  /**< flag for calling background_at_eta and return little information */
  short normal_info; /**< flag for calling background_at_eta and return medium information */
  short long_info;   /**< flag for calling background_at_eta and return all information */

  short inter_normal;  /**< flag for calling background_at_eta and find position in interpolation table normally */
  short inter_closeby; /**< flag for calling background_at_eta and find position in interpolation table starting from previous position in previous call */

  //@}

  /** @name - technical parameters */

  //@{

  short shooting_failed;  /**< flag is set to true if shooting failed. */

  ErrorMsg shooting_error; /**< Error message from shooting failed. */

  short background_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */

  //@}
};

/**
 * temporary parameters and workspace passed to phase space distribution function
 */

struct background_parameters_for_distributions {

  /* structures containing fixed input parameters (indices, ...) */
  struct background * pba;

  /* Additional parameters */

  /* Index of current distribution function */
  int n_ncdm;

  /* Used for interpolating in file of tabulated p-s-d: */
  int tablesize;
  double *q;
  double *f0;
  double *d2f0;
  int last_index;

};
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif
  int background_free_input(
                            struct background *pba
                            );
  int background_ncdm_distribution(
				  void *pba,
				  double q,
				  double * f0
				  );

  int background_ncdm_test_function(
				     void *pba,
				     double q,
				     double * test
				     );

  int background_ncdm_init(
			    struct precision *ppr,
			    struct background *pba
			    );


  int background_ncdm_momenta(
                             double * qvec,
                             double * wvec,
                             int qsize,
                             double M,
                             double factor,
                             double z,
                             double * n,
		             double * rho,
                             double * p,
                             double * drho_dM,
			     double * pseudo_p
                             );

  int background_ncdm_M_from_Omega(
				    struct precision *ppr,
				    struct background *pba,
					int species
				    );
#ifdef __cplusplus
}
#endif
/**
 * @name Some conversion factors and fundamental constants needed by background module:
 */

//@{

#define _Mpc_over_m_ 3.085677581282e22  /**< conversion factor from meters to megaparsecs */
/* remark: CAMB uses 3.085678e22: good to know if you want to compare  with high accuracy */

#define _Gyr_over_Mpc_ 3.06601394e2 /**< conversion factor from megaparsecs to gigayears
				         (c=1 units, Julian years of 365.25 days) */
#define _c_ 2.99792458e8            /**< c in m/s */
#define _G_ 6.67428e-11             /**< Newton constant in m^3/Kg/s^2 */
#define _eV_ 1.602176487e-19        /**< 1 eV expressed in J */

/* parameters entering in Stefan-Boltzmann constant sigma_B */
#define _k_B_ 1.3806504e-23
#define _h_P_ 6.62606896e-34
/* remark: sigma_B = 2 pi^5 k_B^4 / (15h^3c^2) = 5.670400e-8
                   = Stefan-Boltzmann constant in W/m^2/K^4 = Kg/K^4/s^3 */

//@}

/**
 * @name Some limits on possible background parameters
 */

//@{

#define _h_BIG_ 1.5            /**< maximal \f$ h \f$ */
#define _h_SMALL_ 0.3         /**< minimal \f$ h \f$ */
#define _omegab_BIG_ 0.039    /**< maximal \f$ omega_b \f$ */
#define _omegab_SMALL_ 0.005  /**< minimal \f$ omega_b \f$ */

//@}

/**
 * @name Some limits imposed in other parts of the module:
 */

//@{

#define _SCALE_BACK_ 0.1  /**< logarithmic step used when searching
			     for an initial scale factor at which ncdm
			     are still relativistic */

#define _PSD_DERIVATIVE_EXP_MIN_ -30 /**< for ncdm, for accurate computation of dlnf0/dlnq, q step is varied in range specified by these parameters */
#define _PSD_DERIVATIVE_EXP_MAX_ 2  /**< for ncdm, for accurate computation of dlnf0/dlnq, q step is varied in range specified by these parameters */

#define _zeta3_ 1.2020569031595942853997381615114499907649862923404988817922 /**< for quandrature test function */
#define _zeta5_ 1.0369277551433699263313654864570341680570809195019128119741 /**< for quandrature test function */

//@}


#endif
/* @endcond */
