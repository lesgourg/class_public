#include "macros_precision.h"

/*
 * Background Quantities
 */

/**
 * Default initial value of scale factor used in the integration of background quantities.
 * For models like ncdm, the code may decide to start the integration earlier.
 */
class_precision_parameter(a_ini_over_a_today_default,double,1.e-14)
/**
 * Number of background integration steps that are stored in the output vector
 */
class_precision_parameter(background_Nloga,int,40000)
/**
 * Evolver to be used for thermodynamics (rk, ndf15)
 */
class_type_parameter(background_evolver,int,enum evolver_type,ndf15)
/**
 * Tolerance of the background integration, giving the allowed relative integration error.
 * (used by both evolvers)
 */
class_precision_parameter(tol_background_integration,double,1.e-10)
/**
 * Only relevant for rk evolver: the default integration step is given
 * by this number multiplied by the timescale defined in
 * background_timescale (given by the sampling step)
 */
class_precision_parameter(background_integration_stepsize,double,0.5)
/**
 * Tolerance of the deviation of \f$ \Omega_r \f$ from 1 for which to start integration:
 * The starting point of integration will be chosen,
 * such that the Omega of radiation at that point is close to 1 within tolerance.
 * (Class starts background integration during complete radiation domination)
 */
class_precision_parameter(tol_initial_Omega_r,double,1.e-4)
/**
 * Tolerance of relative deviation of the used non-cold dark matter mass compared to that which would give the correct density.
 * The dark matter mass is estimated from the dark matter density using a Newton-Method.
 * In the nonrelativistic limit, this could be estimated using M=density/number density
 */
class_precision_parameter(tol_M_ncdm,double,1.e-7)
/**
 * Tolerance on the relative precision of the integration over
 * non-cold dark matter phase-space distributions.
 */
class_precision_parameter(tol_ncdm,double,1.e-3)
/**
 * Tolerance on the relative precision of the integration over
 * non-cold dark matter phase-space distributions in the synchronous gauge.
 */
class_precision_parameter(tol_ncdm_synchronous,double,1.e-3)
/**
 * Tolerance on the relative precision of the integration over
 * non-cold dark matter phase-space distributions in the newtonian gauge.
 */
class_precision_parameter(tol_ncdm_newtonian,double,1.e-5)
/**
 * Tolerance on the relative precision of the integration over
 * non-cold dark matter phase-space distributions during the background evolution.
 */
class_precision_parameter(tol_ncdm_bg,double,1.e-5)
/**
 * Tolerance on the initial deviation of non-cold dark matter from being fully relativistic.
 * Using w = pressure/density, this quantifies the maximum deviation from 1/3. (for relativistic species)
 */
class_precision_parameter(tol_ncdm_initial_w,double,1.e-3)
/**
 * Tolerance on the deviation of the conformal time of equality from the true value in 1/Mpc.
 */
class_precision_parameter(tol_tau_eq,double,1.e-6)

/**
 * Minimum amount of cdm to allow calculations in synchronous gauge comoving with cdm.
 */
class_precision_parameter(Omega0_cdm_min_synchronous,double,1.e-10)

/**
 * Absolute tolerance of root x during shooting (only 2D case)
 */
class_precision_parameter(tol_shooting_deltax,double,1.e-4)
/**
 * Absolute tolerance of function value F during shooting (only 2D case)
 */
class_precision_parameter(tol_shooting_deltaF,double,1.e-6)
/**
 * Relative tolerance of root x during shooting (only 1D case)
 */
class_precision_parameter(tol_shooting_deltax_rel,double,1.e-5)
/**
 * Tolerance on input of various fractions (e.g. f_idm)
 */
class_precision_parameter(tol_fraction_accuracy,double,1.e-10)
/**
 * Threshold value of M_ncdm=T_ncdm/m_ncdm above wich a species is
 * considered a "non-free-streaming" when comuting the parameter
 * Omega0_nfsm, relevant for HyRec and non-linear correction
 * algorithms
 */
class_precision_parameter(M_nfsm_threshold,double,1.e4)
/*
 * Currently unused parameter.
 */

//class_precision_parameter(safe_phi_scf,double,0.0)
/**
 * Big Bang Nucleosynthesis file path. The file specifies the predictions for
 * \f$ Y_\mathrm{He} \f$ for given \f$ \omega_b \f$ and \f$ N_\mathrm{eff} \f$.
 */
class_string_parameter(sBBN_file,"/external/bbn/sBBN_2025.dat","sBBN file")

/*
 *  Thermodynamical quantities
 */

/**
 * The initial z for the calculation of the recombination history
 */
class_precision_parameter(thermo_z_initial,double,5.e6)
/**
 * The initial z for the calculation of the recombination history in
 * presence of idm (unless the later is tightly-coupled at this
 * redshift)
 */
class_precision_parameter(thermo_z_initial_if_idm,double,1.e9)
/**
 * The switch z for the recfast calculation towards linear sampling
 */
class_precision_parameter(thermo_z_linear,double,1.e4)
/**
 * Number of recfast integration steps (linear sampling, intermdiate times between z_linear and reionization)
 */
class_precision_parameter(thermo_Nz_lin,int,20000)
/**
 * Number of recfast integration steps (logarithmnic sampling. early times between z-initial and z_linear)
 */
class_precision_parameter(thermo_Nz_log,int,5000)
/**
 * Evolver to be used for thermodynamics (rk, ndf15)
 */
class_type_parameter(thermo_evolver,int,enum evolver_type,ndf15)
/**
 * Tolerance of the relative value of integral during thermodynamical integration
 * (used by both evolvers)
 */
class_precision_parameter(tol_thermo_integration,double,1.0e-6)
/**
 * Only relevant for rk evolver: the default integration step is given
 * by this number multiplied by the timescale defined in
 * thermodynamics_timescale (given by the sampling step)
 */
class_precision_parameter(thermo_integration_stepsize,double,0.1)
/**
 * Smoothing in redshift of the variation rate of \f$ \exp(-\kappa) \f$, g, and \f$ \frac{dg}{d\tau} \f$ that is used as a timescale afterwards
 */
class_precision_parameter(thermo_rate_smoothing_radius,int,50)
/**
 * Redshift at which CLASS starts to test for too early re-ionization and/or incomplete recombination.
 */
class_precision_parameter(z_end_reco_test,double,500.)
/**
 * Number of sampling points in the case of primordial black holes in ln(1+z)
 */
class_precision_parameter(primordial_black_hole_Nz,int,75000)
/**
 * Number of sampling points in the case of the coarse sampling in noninjection.c in ln(1+z)
 */
class_precision_parameter(noninjection_Nz_log,int,1000)
/**
 * Number of discrete wavenumbers for dissipation of acoustic waves (found to be giving a reasonable precision)
 */
class_precision_parameter(noninjection_Nk_acc_diss,int,500)
class_precision_parameter(k_min_acc_diss,double,0.12) /**< Minimum wavenumber for dissipation of acoustic waves */
class_precision_parameter(k_max_acc_diss,double,1.e6) /**< Maximum wavenumber for dissipation of acoustic waves */
class_precision_parameter(z_wkb_acc_diss,double,1.e6) /**< Redshift of the WKB approximation for diss. of acoustic waves */

/*
 * Recfast 1.4/1.5 parameters
 */

class_precision_parameter(recfast_Heswitch,int,6)       /**< from recfast 1.4, specifies how accurate the Helium recombination should be handled */
class_precision_parameter(recfast_fudge_He,double,0.86) /**< from recfast 1.4, fugde factor for Peeble's equation coefficient of Helium */
class_precision_parameter(recfast_Hswitch,int,_TRUE_)   /**< from recfast 1.5, specifies how accurate the Hydrogen recombination should be handled */
class_precision_parameter(recfast_fudge_H,double,1.14)  /**< from recfast 1.4, fudge factor for Peeble's equation coeffient of Hydrogen */
class_precision_parameter(recfast_delta_fudge_H,double,-0.015) /**< from recfast 1.5.2, increasing Hydrogen fudge factor if Hswitch is enabled */
class_precision_parameter(recfast_AGauss1,double,-0.14) /**< from recfast 1.5, Gaussian Peeble prefactor fit, amplitude */
class_precision_parameter(recfast_AGauss2,double,0.079) /**< from recfast 1.5.2, Gaussian Peeble prefactor fit, amplitude */
class_precision_parameter(recfast_zGauss1,double,7.28)  /**< from recfast 1.5, Gaussian Peeble prefactor fit, center */
class_precision_parameter(recfast_zGauss2,double,6.73)  /**< from recfast 1.5.2, Gaussian Peeble prefactor fit, center */
class_precision_parameter(recfast_wGauss1,double,0.18)  /**< from recfast 1.5, Gaussian Peeble prefactor fit, width */
class_precision_parameter(recfast_wGauss2,double,0.33)  /**< from recfast 1.5, Gaussian Peeble prefactor fit, width */

class_precision_parameter(recfast_z_He_1,double,8000.0) /**< from recfast 1.4, Starting value of Helium recombination 1 */
class_precision_parameter(recfast_delta_z_He_1,double,50.0) /**< Smoothing factor for recombination approximation switching, found to be OK on 3.09.10 */
class_precision_parameter(recfast_z_He_2,double,4500.0) /**< from recfast 1.4, Ending value of Helium recombination 1, changed on 28.10.20 from 5000 to 4500 */
class_precision_parameter(recfast_delta_z_He_2,double,100.0)/**< Smoothing factor for recombination approximation switching, found to be OK on 3.09.10 */
class_precision_parameter(recfast_z_He_3,double,3500.0) /**< from recfast 1.4, Starting value of Helium recombination 2 */
class_precision_parameter(recfast_delta_z_He_3,double,50.0) /**< Smoothing factor for recombination approximation switching, found to be OK on 3.09.10 */

class_precision_parameter(recfast_z_early_H_recombination,double,2870.) /**< from class 3.0, redshift at beginning of early H-recombination (analytic approximation possible), replaces condition  */
class_precision_parameter(recfast_delta_z_early_H_recombination,double,50.) /**< from class 3.0, smoothing radius delta z for beginning of early H-recombination period  */

class_precision_parameter(recfast_z_full_H_recombination,double,1600.)  /**< from class 3.0, redshift at beignning of full H recombination (always use full equations), replaces condition x_H <  recfast_x_H0_trigger */
class_precision_parameter(recfast_delta_z_full_H_recombination,double,50.)  /**< from class 3.0, smoothing radius delta z for full H-recombination period  */

class_precision_parameter(recfast_delta_z_reio,double,2.)  /**< from class 3.0, smoothing radius delta z for reionization period  */

class_precision_parameter(recfast_x_He0_trigger,double,0.995) /**< Switch for Helium full calculation during reco, raised from 0.99 to 0.995 for smoother Helium */
class_precision_parameter(recfast_x_He0_trigger2,double,0.995)     /**< Switch for Helium full calculation during reco, for changing Helium flag, raised from 0.985 to same as previous one for smoother Helium */
class_precision_parameter(recfast_x_He0_trigger_delta,double,0.05) /**< Smoothing factor for recombination approximation switching, found to be OK on 3.09.10 */
class_precision_parameter(recfast_x_H0_trigger,double,0.995)       /**< Switch for Hydrogen full calculation during reco, raised from 0.99 to 0.995 for smoother Hydrogen */
class_precision_parameter(recfast_x_H0_trigger2,double,0.995)      /**< Switch for Hydrogen full calculation during reco, for changing Hydrogen flag, raised from 0.98 to same as previous one for smoother Hydrogen */
class_precision_parameter(recfast_x_H0_trigger_delta,double,0.05)  /**< Smoothing factor for recombination approximation switching, found to be OK on 3.09.10 */

//class_precision_parameter(recfast_H_frac,double,1.0e-3)  /**< from recfast 1.4, specifies the time at which the temperature evolution is calculated by the more precise equation, not used currently */
/**
 * This is an important flag for energy injections! It also modifies wether recfast will switch approximation schemes or not.
 */
class_precision_parameter(recfast_z_switch_late,double,800.)

/*
 * Hyrec Parameters
 */

class_string_parameter(hyrec_path,"/external/HyRec2020/","hyrec_path") /**< Path to hyrec */

/*
 * Reionization parameters
 */

class_precision_parameter(reionization_z_start_max,double,50.0) /**< Maximum starting value in z for reionization */
class_precision_parameter(reionization_sampling,double,1.5e-2)  /**< Minimum sampling density in z during reionization */
class_precision_parameter(reionization_optical_depth_tol,double,1.0e-4) /**< Relative tolerance on finding the user-given optical depth of reionization given a certain redshift of reionization */
class_precision_parameter(reionization_start_factor,double,8.0) /**< Searching optical depth corresponding to the redshift is started from an initial offset beyond z_reionization_start, multiplied by reionization_width */

/*
 * Heating parameters
 */

class_string_parameter(chi_z_Galli,"/external/heating/Galli_et_al_2013.dat","Galli_file") /**< File containing the chi approximation according to Galli et al 2013 */
class_precision_parameter(z_start_chi_approx,double,2.0e3) /**< Switching redshift from full heating to chosen approx for deposition function */

/*
 * Perturbation parameters
 */

class_precision_parameter(k_min_tau0,double,0.1) /**< number defining k_min for the computation of Cl's and P(k)'s (dimensionless): (k_min tau_0), usually chosen much smaller than one */

class_precision_parameter(k_max_tau0_over_l_max,double,1.8) /**< number defining k_max for the computation of Cl's (dimensionless): (k_max tau_0)/l_max, usually chosen around two. In v3.2.2: lowered from 2.4 to 1.8, because the high value 2.4 was needed to keep CMB lensing accurate enough. With the new Limber scheme, this will be the case anyway, and k_max can be lowered in other observables in order to speed up the code. */
class_precision_parameter(k_step_sub,double,0.05) /**< step in k space, in units of one period of acoustic oscillation at decoupling, for scales inside sound horizon at decoupling */
class_precision_parameter(k_step_super,double,0.002) /**< step in k space, in units of one period of acoustic oscillation at decoupling, for scales above sound horizon at decoupling */
class_precision_parameter(k_step_transition,double,0.2) /**< dimensionless number regulating the transition from 'sub' steps to 'super' steps. Decrease for more precision. */
class_precision_parameter(k_step_super_reduction,double,0.1) /**< the step k_step_super is reduced by this amount in the k-->0 limit (below scale of Hubble and/or curvature radius) */

class_precision_parameter(k_per_decade_for_pk,double,10.0) /**< if values needed between kmax inferred from k_oscillations and k_kmax_for_pk, this gives the number of k per decade outside the BAO region*/

class_precision_parameter(idmdr_boost_k_per_decade_for_pk,double,1.0) /**< boost factor for the case of DAO in idm-idr models */

class_precision_parameter(k_per_decade_for_bao,double,70.0) /**< if values needed between kmax inferred from k_oscillations and k_kmax_for_pk, this gives the number of k per decade inside the BAO region (for finer sampling)*/

class_precision_parameter(k_bao_center,double,3.0) /**< in ln(k) space, the central value of the BAO region where sampling is finer is defined as k_rec times this number (recommended: 3, i.e. finest sampling near 3rd BAO peak) */

class_precision_parameter(k_bao_width,double,4.0) /**< in ln(k) space, width of the BAO region where sampling is finer: this number gives roughly the number of BAO oscillations well resolved on both sides of the central value (recommended: 4, i.e. finest sampling from before first up to 3+4=7th peak) */

class_precision_parameter(start_small_k_at_tau_c_over_tau_h,double,0.0015) /**< largest wavelengths start being sampled when universe is sufficiently opaque. This is quantified in terms of the ratio of thermo to hubble time scales, \f$ \tau_c/\tau_H \f$. Start when start_largek_at_tau_c_over_tau_h equals this ratio. Decrease this value to start integrating the wavenumbers earlier in time. */

class_precision_parameter(start_large_k_at_tau_h_over_tau_k,double,0.07)  /**< largest wavelengths start being sampled when mode is sufficiently outside Hubble scale. This is quantified in terms of the ratio of hubble time scale to wavenumber time scale, \f$ \tau_h/\tau_k \f$ which is roughly equal to (k*tau). Start when this ratio equals start_large_k_at_tau_k_over_tau_h. Decrease this value to start integrating the wavenumbers earlier in time. */

/**
 * when to switch off tight-coupling approximation: first condition:
 * \f$ \tau_c/\tau_H \f$ > tight_coupling_trigger_tau_c_over_tau_h.
 * Decrease this value to switch off earlier in time.  If this
 * number is larger than start_sources_at_tau_c_over_tau_h, the code
 * returns an error, because the source computation requires
 * tight-coupling to be switched off.
 */
class_precision_parameter(tight_coupling_trigger_tau_c_over_tau_h,double,0.015)

/**
 * when to switch off tight-coupling approximation:
 * second condition: \f$ \tau_c/\tau_k \equiv k \tau_c \f$ <
 * tight_coupling_trigger_tau_c_over_tau_k.
 * Decrease this value to switch off earlier in time.
 */
class_precision_parameter(tight_coupling_trigger_tau_c_over_tau_k,double,0.01)

/**
 * when to switch off tight-coupling approximation:
 * third condition: for the case of idm with photons.
 */
class_precision_parameter(tight_coupling_trigger_tau_c_over_tau_dmu_idm_g, double, 0.01);

/**
 * when to switch off tight-coupling approximation:
 * fourth condition: for the case of idm with baryons.
 */
class_precision_parameter(tight_coupling_trigger_tau_c_over_tau_R_idm_b, double, 0.01)

class_precision_parameter(start_sources_at_tau_c_over_tau_h,double,0.008) /**< sources start being sampled when universe is sufficiently opaque. This is quantified in terms of the ratio of thermo to hubble time scales, \f$ \tau_c/\tau_H \f$. Start when start_sources_at_tau_c_over_tau_h equals this ratio. Decrease this value to start sampling the sources earlier in time. */

class_precision_parameter(tight_coupling_approximation,int,(int)compromise_CLASS) /**< method for tight coupling approximation */

class_precision_parameter(idm_dr_tight_coupling_trigger_tau_c_over_tau_k,double,0.01)  /**< when to switch off the dark-tight-coupling approximation, first condition (see normal tca for full definition) */
class_precision_parameter(idm_dr_tight_coupling_trigger_tau_c_over_tau_h,double,0.015) /**< when to switch off the dark-tight-coupling approximation, second condition (see normal tca for full definition) */

class_precision_parameter(l_max_g,int,12)     /**< number of momenta in Boltzmann hierarchy for photon temperature (scalar), at least 4 */
class_precision_parameter(l_max_pol_g,int,10) /**< number of momenta in Boltzmann hierarchy for photon polarization (scalar), at least 4 */
class_precision_parameter(l_max_dr,int,17)   /**< number of momenta in Boltzmann hierarchy for decay radiation, at least 4 */
class_precision_parameter(l_max_ur,int,17)   /**< number of momenta in Boltzmann hierarchy for relativistic neutrino/relics (scalar), at least 4 */
class_precision_parameter(l_max_idr,int,17)   /**< number of momenta in Boltzmann hierarchy for interacting dark radiation */
class_precision_parameter(l_max_ncdm,int,17)   /**< number of momenta in Boltzmann hierarchy for relativistic neutrino/relics (scalar), at least 4 */
class_precision_parameter(l_max_g_ten,int,5)     /**< number of momenta in Boltzmann hierarchy for photon temperature (tensor), at least 4 */
class_precision_parameter(l_max_pol_g_ten,int,5) /**< number of momenta in Boltzmann hierarchy for photon polarization (tensor), at least 4 */

class_precision_parameter(curvature_ini,double,1.0)     /**< initial condition for curvature for adiabatic */
class_precision_parameter(entropy_ini,double,1.0) /**< initial condition for entropy perturbation for isocurvature */
class_precision_parameter(gw_ini,double,1.0)      /**< initial condition for tensor metric perturbation h */

/**
 * default step \f$ d \tau \f$ in perturbation integration, in units of the timescale involved in the equations (usually, the min of \f$ 1/k \f$, \f$ 1/aH \f$, \f$ 1/\dot{\kappa} \f$)
 */
class_precision_parameter(perturbations_integration_stepsize,double,0.5)
/**
 * default step \f$ d \tau \f$ for sampling the source function, in units of the timescale involved in the sources: \f$ (\dot{\kappa}- \ddot{\kappa}/\dot{\kappa})^{-1} \f$
 */
class_precision_parameter(perturbations_sampling_stepsize,double,0.1)
/**
 * added in v 3.2.2: age fraction (between 0 and 1 ) such that, when
 * tau > conformal_age * age_fraction, the time sampling of sources is
 * twice finer, in order to boost the accuracy of the lensing
 * line-of-sight integrals (for l < l_switch_limber) without changing
 * that of unlensed CMB observables. Setting to 1.0 disables this
 * functionality.
*/
class_precision_parameter(perturbations_sampling_boost_above_age_fraction, double, 0.9)
/**
 * control parameter for the precision of the perturbation integration,
 * IMPORTANT FOR SETTING THE STEPSIZE OF NDF15
 */
class_precision_parameter(tol_perturbations_integration,double,1.0e-5)
/**
 * cutoff relevant for controlling stiffness in the PPF scheme. It is
 * neccessary for the Runge-Kutta evolver, but not for ndf15. However,
 * the approximation is excellent for a cutoff value of 1000, so we
 * leave it on for both evolvers. (CAMB uses a cutoff value of 30.)
 */
class_precision_parameter(c_gamma_k_H_square_max,double,1.0e3)
/**
 * precision with which the code should determine (by bisection) the
 * times at which sources start being sampled, and at which
 * approximations must be switched on/off (units of Mpc)
 */
class_precision_parameter(tol_tau_approx,double,1.0e-10)
/**
 * method for switching off photon perturbations
 */
class_precision_parameter(radiation_streaming_approximation,int,rsa_MD_with_reio)
/**
 * when to switch off photon perturbations, ie when to switch
 * on photon free-streaming approximation (keep density and thtau, set
 * shear and higher momenta to zero):
 * first condition: \f$ k \tau \f$ > radiation_streaming_trigger_tau_h_over_tau_k
 */
class_precision_parameter(radiation_streaming_trigger_tau_over_tau_k,double,45.0)
/**
 * when to switch off photon perturbations, ie when to switch
 * on photon free-streaming approximation (keep density and theta, set
 * shear and higher momenta to zero):
 * second condition:
 */
class_precision_parameter(radiation_streaming_trigger_tau_c_over_tau,double,5.0)

class_precision_parameter(idr_streaming_approximation,int,rsa_idr_none) /**< method for dark radiation free-streaming approximation */
class_precision_parameter(idr_streaming_trigger_tau_over_tau_k,double,50.0) /**< when to switch on dark radiation (idr) free-streaming approximation, first condition */
class_precision_parameter(idr_streaming_trigger_tau_c_over_tau,double,10.0) /**< when to switch on dark radiation (idr) free-streaming approximation, second condition */

class_precision_parameter(ur_fluid_approximation,int,ufa_CLASS) /**< method for ultra relativistic fluid approximation */
/**
 * when to switch off ur (massless neutrinos / ultra-relativistic
 * relics) fluid approximation
 */
class_precision_parameter(ur_fluid_trigger_tau_over_tau_k,double,30.0)
class_precision_parameter(ncdm_fluid_approximation,int,ncdmfa_CLASS) /**< method for non-cold dark matter fluid approximation */
/**
 * when to switch off ncdm (massive neutrinos / non-cold
 * relics) fluid approximation
 */
class_precision_parameter(ncdm_fluid_trigger_tau_over_tau_k,double,31.0)
/**
 * whether CMB source functions can be approximated as zero when
 * visibility function g(tau) is tiny
 */
class_precision_parameter(neglect_CMB_sources_below_visibility,double,1.0e-3)
/**
 * The type of evolver to use: options are ndf15 or rk
 */
class_type_parameter(evolver,int,enum evolver_type,ndf15)

/*
 * Primordial parameters
 */

class_precision_parameter(k_per_decade_primordial,double,10.0) /**< logarithmic sampling for primordial spectra (number of points per decade in k space) */

class_precision_parameter(primordial_inflation_ratio_min,double,100.0) /**< for each k, start following wavenumber when aH = k/primordial_inflation_ratio_min */
class_precision_parameter(primordial_inflation_ratio_max,double,1.0/50.0) /**< for each k, stop following wavenumber, at the latest, when aH = k/primordial_inflation_ratio_max */
class_precision_parameter(primordial_inflation_phi_ini_maxit,int,10000)      /**< maximum number of iteration when searching a suitable initial field value phi_ini (value reached when no long-enough slow-roll period before the pivot scale) */
class_precision_parameter(primordial_inflation_pt_stepsize,double,0.01)     /**< controls the integration timestep for inflaton perturbations */
class_precision_parameter(primordial_inflation_bg_stepsize,double,0.005)     /**< controls the integration timestep for inflaton background */
class_precision_parameter(primordial_inflation_tol_integration,double,1.0e-3) /**< controls the precision of the ODE integration during inflation */
class_precision_parameter(primordial_inflation_attractor_precision_pivot,double,0.001)   /**< targeted precision when searching attractor solution near phi_pivot */
class_precision_parameter(primordial_inflation_attractor_precision_initial,double,0.1) /**< targeted precision when searching attractor solution near phi_ini */
class_precision_parameter(primordial_inflation_attractor_maxit,int,10) /**< maximum number of iteration when searching attractor solution */
class_precision_parameter(primordial_inflation_tol_curvature,double,1.0e-3) /**< for each k, stop following wavenumber, at the latest, when curvature perturbation R is stable up to to this tolerance */
class_precision_parameter(primordial_inflation_aH_ini_target,double,0.9) /**< control the step size in the search for a suitable initial field value */
class_precision_parameter(primordial_inflation_end_dphi,double,1.0e-10) /**< first bracketing width, when trying to bracket the value phi_end at which inflation ends naturally */
class_precision_parameter(primordial_inflation_end_logstep,double,10.0) /**< logarithmic step for updating the bracketing width, when trying to bracket the value phi_end at which inflation ends naturally */
class_precision_parameter(primordial_inflation_small_epsilon,double,0.1) /**< value of slow-roll parameter epsilon used to define a field value phi_end close to the end of inflation (doesn't need to be exactly at the end): epsilon(phi_end)=small_epsilon (should be smaller than one) */
class_precision_parameter(primordial_inflation_small_epsilon_tol,double,0.01) /**< tolerance in the search for phi_end */
class_precision_parameter(primordial_inflation_extra_efolds,double,2.0) /**< a small number of efolds, irrelevant at the end, used in the search for the pivot scale (backward from the end of inflation) */

/*
 * Transfer function parameters
 */

class_precision_parameter(l_linstep,int,40) /**< factor for logarithmic spacing of values of l over which bessel and transfer functions are sampled */

class_precision_parameter(l_logstep,double,1.12) /**< maximum spacing of values of l over which Bessel and transfer functions are sampled (so, spacing becomes linear instead of logarithmic at some point) */

class_precision_parameter(hyper_x_min,double,1.0e-5)  /**< flat case: lower bound on the smallest value of x at which we sample \f$ \Phi_l^{\nu}(x)\f$ or \f$ j_l(x)\f$ */
class_precision_parameter(hyper_sampling_flat,double,8.0)  /**< flat case: number of sampled points x per approximate wavelength \f$ 2\pi \f$, should remain >7.5 */
class_precision_parameter(hyper_sampling_curved_low_nu,double,7.0)  /**< open/closed cases: number of sampled points x per approximate wavelength \f$ 2\pi/\nu\f$, when \f$ \nu \f$ smaller than hyper_nu_sampling_step */
class_precision_parameter(hyper_sampling_curved_high_nu,double,3.0) /**< open/closed cases: number of sampled points x per approximate wavelength \f$ 2\pi/\nu\f$, when \f$ \nu \f$ greater than hyper_nu_sampling_step */
class_precision_parameter(hyper_nu_sampling_step,double,1000.0)  /**< open/closed cases: value of nu at which sampling changes  */
class_precision_parameter(hyper_phi_min_abs,double,1.0e-10)  /**< small value of Bessel function used in calculation of first point x (\f$ \Phi_l^{\nu}(x) \f$ equals hyper_phi_min_abs) */
class_precision_parameter(hyper_x_tol,double,1.0e-4)  /**< tolerance parameter used to determine first value of x */
class_precision_parameter(hyper_flat_approximation_nu,double,4000.0)  /**< value of nu below which the flat approximation is used to compute Bessel function */

class_precision_parameter(q_linstep,double,0.45)         /**< asymptotic linear sampling step in q
                               space, in units of \f$ 2\pi/r_a(\tau_rec) \f$
                               (comoving angular diameter distance to
                               recombination), very important for CMB */

class_precision_parameter(q_logstep_spline,double,170.0) /**< initial logarithmic sampling step in q
                                space, in units of \f$ 2\pi/r_a(\tau_{rec})\f$
                                (comoving angular diameter distance to
                                recombination), very important for CMB and LSS */

class_precision_parameter(q_logstep_open,double,6.0)   /**< in open models, the value of
                                q_logstep_spline must be decreased
                                according to curvature. Increasing
                                this number will make the calculation
                                more accurate for large positive
                                Omega_k */

class_precision_parameter(q_logstep_trapzd,double,20.0) /**< initial logarithmic sampling step in q
                                space, in units of \f$ 2\pi/r_a(\tau_{rec}) \f$
                                (comoving angular diameter distance to
                                recombination), in the case of small
                                q's in the closed case, for which one
                                must used trapezoidal integration
                                instead of spline (the number of q's
                                for which this is the case decreases
                                with curvature and vanishes in the
                                flat limit) */

class_precision_parameter(q_numstep_transition,double,250.0) /**< number of steps for the transition
                                 from q_logstep_trapzd steps to
                                 q_logstep_spline steps (transition
                                 must be smooth for spline) */

class_precision_parameter(q_logstep_limber,double,1.025) /**< new in v3.2.2: in the new 'full limber' scheme, logarithmic step for the k-grid (and q-grid) */
class_precision_parameter(k_max_limber_over_l_max_scalars,double,0.001) /**< new in v3.2.2: in the new 'full limber' scheme, the integral runs up to k_max = l_max_scalars times this parameter (units of 1/Mpc) */

class_precision_parameter(transfer_neglect_delta_k_S_t0,double,0.15) /**< for temperature source function T0 of scalar mode, range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero */
class_precision_parameter(transfer_neglect_delta_k_S_t1,double,0.04) /**< same for temperature source function T1 of scalar mode */
class_precision_parameter(transfer_neglect_delta_k_S_t2,double,0.15) /**< same for temperature source function T2 of scalar mode */
class_precision_parameter(transfer_neglect_delta_k_S_e,double,0.11)  /**< same for polarization source function E of scalar mode */
class_precision_parameter(transfer_neglect_delta_k_V_t1,double,1.0) /**< same for temperature source function T1 of vector mode */
class_precision_parameter(transfer_neglect_delta_k_V_t2,double,1.0) /**< same for temperature source function T2 of vector mode */
class_precision_parameter(transfer_neglect_delta_k_V_e,double,1.0)  /**< same for polarization source function E of vector mode */
class_precision_parameter(transfer_neglect_delta_k_V_b,double,1.0)  /**< same for polarization source function B of vector mode */
class_precision_parameter(transfer_neglect_delta_k_T_t2,double,0.2) /**< same for temperature source function T2 of tensor mode */
class_precision_parameter(transfer_neglect_delta_k_T_e,double,0.25)  /**< same for polarization source function E of tensor mode */
class_precision_parameter(transfer_neglect_delta_k_T_b,double,0.1)  /**< same for polarization source function B of tensor mode */

class_precision_parameter(transfer_neglect_late_source,double,400.0)  /**< value of l below which the CMB source functions can be neglected at late time, excepted when there is a Late ISW contribution */

class_precision_parameter(l_switch_limber,double,10.) /**< when to use the Limber approximation for project gravitational potential cl's */
// For density Cl, we recommend not to use the Limber approximation
// at all, and hence to put here a very large number (e.g. 10000); but
// if you have wide and smooth selection functions you may wish to
// use it; then 100 might be OK
class_precision_parameter(l_switch_limber_for_nc_local_over_z,double,100.0) /**< when to use the Limber approximation for local number count contributions to cl's (relative to central redshift of each bin) */
// For terms integrated along the line-of-sight involving spherical
// Bessel functions (but not their derivatives), Limber
// approximation works well. High precision can be reached with 2000
// only. But if you have wide and smooth selection functions you may
// reduce to e.g. 30.
class_precision_parameter(l_switch_limber_for_nc_los_over_z,double,30.0) /**< when to use the Limber approximation for number count contributions to cl's integrated along the line-of-sight (relative to central redshift of each bin) */

class_precision_parameter(selection_cut_at_sigma,double,5.0)/**< in sigma units, where to cut gaussian selection functions */
class_precision_parameter(selection_sampling,double,50.0) /**< controls sampling of integral over time when selection functions vary quicker than Bessel functions. Increase for better sampling. */
class_precision_parameter(selection_sampling_bessel,double,20.0)/**< controls sampling of integral over time when selection functions vary slower than Bessel functions. Increase for better sampling. IMPORTANT for lensed contributions. */
class_precision_parameter(selection_sampling_bessel_los,double,ppr->selection_sampling_bessel)/**< controls sampling of integral over time when selection functions vary slower than Bessel functions. This parameter is specific to number counts contributions to Cl integrated along the line of sight. Increase for better sampling */
class_precision_parameter(selection_tophat_edge,double,0.1) /**< controls how smooth are the edge of top-hat window function (<<1 for very sharp, 0.1 for sharp) */

/*
 * Fourier module precision parameters
 * */

class_precision_parameter(sigma_k_per_decade,double,80.) /**< logarithmic stepsize controlling the precision of integrals for sigma(R,k) and similar quantitites */

class_precision_parameter(nonlinear_min_k_max,double,5.0) /**< when
                               using an algorithm to compute nonlinear
                               corrections, like halofit or hmcode,
                               k_max must be at least equal to this
                               value. Calculations are done internally
                               until this k_max, but the P(k,z) output
                               is still controlled by P_k_max_1/Mpc or
                               P_k_max_h/Mpc even if they are
                               smaller */

class_precision_parameter(k_max_for_pk_sigma8_min,double,10.) /**< minimal k_max for computation of sigma8 */
class_precision_parameter(k_max_for_pk_sigma8_max,double,100.) /**< maximal k_max for computation of sigma8 */

/** parameters relevant for HALOFIT computation */

class_precision_parameter(halofit_min_k_nonlinear,double,1.0e-4)/**< value of k in 1/Mpc below which non-linear corrections will be neglected */

class_precision_parameter(halofit_k_per_decade,double,80.0) /**< halofit needs to evalute integrals
                               (linear power spectrum times some
                               kernels). They are sampled using
                               this logarithmic step size. */

class_precision_parameter(halofit_sigma_precision,double,0.05) /**< a smaller value will lead to a
                               more precise halofit result at the *highest*
                               redshift at which halofit can make computations,
                               at the expense of requiring a larger k_max; but
                               this parameter is not relevant for the
                               precision on P_nl(k,z) at other redshifts, so
                               there is normally no need to change it */

class_precision_parameter(halofit_tol_sigma,double,1.0e-6) /**< tolerance required on sigma(R) when
                               matching the condition sigma(R_nl)=1,
                               whcih defines the wavenumber of
                               non-linearity, k_nl=1./R_nl */

class_precision_parameter(pk_eq_z_max,double,5.0)  /**< Maximum z for the pk_eq method */
class_precision_parameter(pk_eq_Nzlog,int,10)      /**< Number of logarithmically spaced redshift values for the pk_eq method */
class_precision_parameter(pk_eq_tol,double,1.0e-7) /**< Tolerance on the pk_eq method for finding the pk */

/** Parameters relevant for HMcode computation */

class_precision_parameter(hmcode_max_k_extra,double,1.e6) /**< parameter specifying the maximum k value for
                                                             the extrapolation of the linear power spectrum
                                                             (needed for the sigma computation) */

class_precision_parameter(hmcode_tol_sigma,double,1.e-6) /**< tolerance required on sigma(R) when matching the
                                                            condition sigma(R_nl)=1, which defines the wavenumber
                                                            of non-linearity, k_nl=1./R_nl */

/**
 * parameters controlling stepsize and min/max r & a values for
 * sigma(r) & grow table
 */
class_precision_parameter(n_hmcode_tables,int,64)
class_precision_parameter(rmin_for_sigtab,double,1.e-5)
class_precision_parameter(rmax_for_sigtab,double,1.e3)
class_precision_parameter(ainit_for_growtab,double,1.e-3)
class_precision_parameter(amax_for_growtab,double,1.)

/**
 * parameters controlling stepsize and min/max halomass values for the
 * 1-halo-power integral
 */
class_precision_parameter(nsteps_for_p1h_integral,int,256)
class_precision_parameter(mmin_for_p1h_integral,double,1.e3)
class_precision_parameter(mmax_for_p1h_integral,double,1.e18)


/*
 * Lensing precision parameters
 */

class_precision_parameter(accurate_lensing,int,_FALSE_) /**< switch between Gauss-Legendre quadrature integration and simple quadrature on a subdomain of angles */
class_precision_parameter(num_mu_minus_lmax,int,70) /**< difference between num_mu and l_max, increase for more precision */
class_precision_parameter(delta_l_max,int,500)/**< difference between l_max in unlensed and lensed spectra */
class_precision_parameter(tol_gauss_legendre,double,ppr->smallest_allowed_variation) /**< tolerance with which quadrature points are found: must be very small for an accurate integration (if not entered manually, set automatically to match machine precision) */

/*
 * Spectral distortions precision parameters
 */

class_precision_parameter(sd_z_min,double,1.02e3)
class_precision_parameter(sd_z_max,double,5.0e6)
class_precision_parameter(sd_z_size,int,400)

class_precision_parameter(sd_x_min,double,1.0e-2)
class_precision_parameter(sd_x_max,double,5.0e1)
class_precision_parameter(sd_x_size,int,500)

/**
 * Tolerance on the deviation of the distortions detector quality
 */
class_precision_parameter(tol_sd_detector,double,1.e-5)

class_string_parameter(sd_external_path,"/external/distortions","sd_external_path")


#undef class_precision_parameter
#undef class_string_parameter
#undef class_type_parameter
