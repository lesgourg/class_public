/**
 * Implementation of HMcode (a la 2015 and 2020)
 * Nils Schoeneberg, Julien Lesgourgues, Samuel Brieden
 * */

#include "hmcode.h"
#include "fourier.h"
#include "parallel.h"
#include "sys/time.h"

/**
 * Computes the nonlinear correction on the linear power spectrum with
 * HMcode 2015 (Mead et al. 1505.07833), 2016 (Mead et al. 1602.02154)
 * or 2020 (Mead et al. 2009.01858), wrapper for 'hmcode_compute'
 * method.
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param ppt        Input: pointer to perturbation structure
 * @param ppm        Input: pointer to primordial structure
 * @param pfo        Input: pointer to fourier structure
 * @param index_pk   Input: index of the pk type, either index_m or index_cb
 * @param index_tau  Input: index of tau, at which to compute the nl correction
 * @param tau        Input: tau, at which to compute the nl correction
 * @param pk_nl      Output:nonlinear power spectrum
 * @param lnpk_l     Input: logarithm of the linear power spectrum for both index_m and index_cb
 * @param ddlnpk_l   Input: spline of the logarithm of the linear power spectrum for both index_m and index_cb
 * @param k_nl       Output: nonlinear scale for index_m and index_cb
 * @param nl_corr_not_computable_at_this_k Ouput: was the computation doable?
 * @param phw        Input/Output: pointer to hmcode workspace
 * @return the error status
 */

int hmcode(
           struct precision *ppr,
           struct background *pba,
           struct perturbations *ppt,
           struct primordial *ppm,
           struct fourier *pfo,
           int index_pk,
           int index_tau,
           double tau,
           double *pk_nl,
           double **lnpk_l,
           double **ddlnpk_l,
           double *k_nl,
           short *nl_corr_not_computable_at_this_k,
           struct hmcode_workspace * phw
           ) {

  double *pk_nl_baseline;
  double *pk_nl_denominator;
  int index_k;

  /* The baryonic contribution is calculated as the ratio of three different
     simpler models, which are computed below
     (baryonic/unfitted_nobaryons * nobaryons) */
  if (pfo->hm_version == hmcode_version_2020_baryonic) {

    class_call(hmcode_compute(ppr, pba, ppt, ppm, pfo,
                              index_pk, index_tau, tau,
                              pk_nl, lnpk_l, ddlnpk_l,
                              k_nl, nl_corr_not_computable_at_this_k,
                              phw),
               pfo->error_message,
               pfo->error_message);

    class_alloc(pk_nl_baseline,
                pfo->k_size*sizeof(double),
                pfo->error_message);

    pfo->hm_version = hmcode_version_2020;

    class_call(hmcode_compute(ppr, pba, ppt, ppm, pfo,
                              index_pk, index_tau, tau,
                              pk_nl_baseline, lnpk_l, ddlnpk_l,
                              k_nl, nl_corr_not_computable_at_this_k,
                              phw),
               pfo->error_message,
               pfo->error_message);

    pfo->hm_version = hmcode_version_2020_unfitted;

    class_alloc(pk_nl_denominator,
                pfo->k_size*sizeof(double),
                pfo->error_message);

    class_call(hmcode_compute(ppr, pba, ppt, ppm, pfo,
                              index_pk, index_tau, tau,
                              pk_nl_denominator, lnpk_l, ddlnpk_l,
                              k_nl, nl_corr_not_computable_at_this_k,
                              phw),
               pfo->error_message,
               pfo->error_message);

    pfo->hm_version = hmcode_version_2020_baryonic;

    for(index_k=0;index_k<pfo->k_size;++index_k){
      pk_nl[index_k] = pk_nl[index_k] * pk_nl_baseline[index_k] / pk_nl_denominator[index_k];
    }

    free(pk_nl_baseline);
    free(pk_nl_denominator);
  }
  else { //In this case we can do the normal analysis
    class_call(hmcode_compute(ppr, pba, ppt, ppm, pfo,
                              index_pk, index_tau, tau,
                              pk_nl, lnpk_l, ddlnpk_l,
                              k_nl, nl_corr_not_computable_at_this_k,
                              phw),
               pfo->error_message,
               pfo->error_message);
  }

  return _SUCCESS_;
}

/**
 * Computes the nonlinear correction on the linear power spectrum via
 * HMcode 2015 (Mead et al. 1505.07833), 2016 (Mead et al. 1602.02154)
 * or 2020 (Mead et al. 2009.01858)
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param ppt        Input: pointer to perturbation structure
 * @param ppm        Input: pointer to primordial structure
 * @param pfo        Input: pointer to fourier structure
 * @param index_pk   Input: index of the pk type, either index_m or index_cb
 * @param index_tau  Input: index of tau, at which to compute the nl correction
 * @param tau        Input: tau, at which to compute the nl correction
 * @param pk_nl      Output:nonlinear power spectrum
 * @param lnpk_l     Input: logarithm of the linear power spectrum for both index_m and index_cb
 * @param ddlnpk_l   Input: spline of the logarithm of the linear power spectrum for both index_m and index_cb
 * @param k_nl       Output: nonlinear scale for index_m and index_cb
 * @param nl_corr_not_computable_at_this_k Ouput: was the computation doable?
 * @param phw        Input/Output: pointer to hmcode workspace
 * @return the error status
 */

int hmcode_compute(
                   struct precision *ppr,
                   struct background *pba,
                   struct perturbations *ppt,
                   struct primordial *ppm,
                   struct fourier *pfo,
                   int index_pk,
                   int index_tau,
                   double tau,
                   double *pk_nl,
                   double **lnpk_l,
                   double **ddlnpk_l,
                   double *k_nl,
                   short *nl_corr_not_computable_at_this_k,
                   struct hmcode_workspace *phw
                   ) {

  /* integers */
  int index_mass, i, ng, nsig;
  int index_k, index_ncol;
  int last_index=0;
  int index_pk_cb;
  int counter, index_nl;

  int index_nu, index_cut;
  int index_y;
  int index_ddy;

  /* Background parameters */
  double Omega_m,fnu,Omega0_m;
  double z_at_tau;
  double rho_crit_today_in_msun_mpc3;
  double growth;
  double anorm;
  double scale_indep_growth;
  double accumulated_growth;
  double Om_norad;

  /* temporary numbers */
  struct hmcode_growth* phg;
  double m, r, nu, sig, sigf;
  double diff, r1, r2, nu1, nu2;
  double a,b,h;
  double *head;
  double dv_x,dv_y;

  /* HMcode parameters */
  double mmin, mmax, nu_min;

  double sigma_disp, sigma_disp100, sigma8, sigma8_cb;
  double delta_c, Delta_v;
  double fraction;

  double sigma_nl, nu_nl, r_nl;
  double sigma_prime;
  double dlnsigdlnR;
  double n_eff;
  double alpha;

  double z_form, g_form;

  double eta;
  double Abary;
  double nu_cut;
  double k_star, fdamp;

  double kdamp;
  double DEcorr;
  double w0, dw_over_da_fld, integral_fld;
  double g_wcdm, g_lcdm;

  /* data fields */
  double *pvecback;
  double *conc;
  double *mass;
  double *sigma_r;
  double *sigmaf_r;
  double *r_virial;
  double *r_real;
  double *nu_arr;

  double *p1h_integrand;

  struct timeval begin, end;
  long seconds;
  long microseconds;
  double elapsed;
  int num_threads=1;
  double a_form;

  phg = phw->phg;

  /** include precision parameters that control the number of entries in the growth and sigma tables */
  ng = ppr->n_hmcode_tables;
  nsig = ppr->n_hmcode_tables;

  /** Compute background quantitites today */

  Omega0_m = pba->Omega0_m;
  fnu      = pba->Omega0_ncdm_tot/Omega0_m;

  /** If index_pk_cb, choose Omega0_cb as the matter density parameter.
   * If index_pk_m, choose Omega0_cbn as the matter density parameter. */
  if (index_pk==pfo->index_pk_cluster){
    Omega0_m = Omega0_m - pba->Omega0_ncdm_tot;
  }

  anorm    = 1./(2*pow(_PI_,2));

  /** Call all the relevant background parameters at this tau */
  class_alloc(pvecback,pba->bg_size*sizeof(double),pfo->error_message);

  class_call(background_at_tau(pba,tau,long_info,inter_normal,&last_index,pvecback),
             pba->error_message,
             pfo->error_message);

  Omega_m = pvecback[pba->index_bg_Omega_m];//TBC (i.e. check if for P_cb here we should use Omega_cb) -- here the total time varying Omega_m is used for delta_c and for Delta_v according to the Mead fit of the Massara simulations.

  /*growth = pvecback[pba->index_bg_D];*/ // This would be the simplest approach, but sadly HMcode and CLASS differ a bit in their growth definitions
  class_call(array_spline_hunt(phg->a_table,
                               phg->a_size,
                               pvecback[pba->index_bg_a],
                               &last_index,
                               &h,
                               &a,
                               &b,
                               pfo->error_message),
             pfo->error_message,
             pfo->error_message);

  growth = array_spline_eval(phg->normgrowth_table,phg->ddnormgrowth_table,last_index,last_index+1,h,a,b); //HM2020

  z_at_tau = 1./pvecback[pba->index_bg_a]-1.;

  /* The number below is the critical density today, rho_c = 3 * H0^2 / 8*pi*G, in units of M_sun over Mpc^3 */
  rho_crit_today_in_msun_mpc3 = 3.*pow(1.e5*pba->h, 2)/8./_PI_/_G_*_Mpc_over_m_/_M_SUN_;

  free(pvecback);

  /* Read additional quantities for HMcode 2020 if required, and setup the 2020 nowiggle workspace */
  if (pfo->hm_version == hmcode_version_2020 || pfo->hm_version == hmcode_version_2020_baryonic || pfo->hm_version == hmcode_version_2020_unfitted){

    class_call(hmcode_nowiggle_init(pfo,lnpk_l,ddlnpk_l,index_pk,phw),
               pfo->error_message,
               pfo->error_message);

    head = phg->growth_table+last_index*phg->gt_size;
    scale_indep_growth = array_spline_eval(head+phg->index_gt_g,head+phg->index_gt_ddg,0,phg->gt_size,h,a,b);
    accumulated_growth = array_spline_eval(head+phg->index_gt_intg,head+phg->index_gt_ddintg,0,phg->gt_size,h,a,b);
    Om_norad = array_spline_eval(head+phg->index_gt_Omnorad,head+phg->index_gt_ddOmnorad,0,phg->gt_size,h,a,b);
  }
  else{
    scale_indep_growth = 0;
    accumulated_growth = 0;
    Om_norad = 0;
  }

  /** Test whether pk_cb has to be taken into account (only if we have massive neutrinos)*/
  if (pba->has_ncdm==_TRUE_){
    index_pk_cb = pfo->index_pk_cb;
  }
  else {
    index_pk_cb = index_pk;
  }


  /** Get sigma(R=8 Mpc/h), sigma_disp(R=0), sigma_disp(R=100 Mpc/h) and write them into pfo structure */

  class_call(fourier_sigmas(pfo,
                            8./pba->h,
                            lnpk_l[index_pk],ddlnpk_l[index_pk],
                            pfo->k_size_extra,
                            ppr->sigma_k_per_decade,
                            out_sigma,
                            &sigma8),
             pfo->error_message,
             pfo->error_message);

  class_call(fourier_sigmas(pfo,
                            8./pba->h,
                            lnpk_l[index_pk_cb],ddlnpk_l[index_pk_cb],
                            pfo->k_size_extra,
                            ppr->sigma_k_per_decade,
                            out_sigma,
                            &sigma8_cb),
             pfo->error_message,
             pfo->error_message);

  class_call(fourier_sigmas(pfo,
                            0.,
                            lnpk_l[index_pk],ddlnpk_l[index_pk],
                            pfo->k_size_extra,
                            ppr->sigma_k_per_decade,
                            out_sigma_disp,
                            &sigma_disp),
             pfo->error_message,
             pfo->error_message);

  class_call(fourier_sigmas(pfo,
                            100./pba->h,
                            lnpk_l[index_pk],ddlnpk_l[index_pk],
                            pfo->k_size_extra,
                            ppr->sigma_k_per_decade,
                            out_sigma_disp,
                            &sigma_disp100),
             pfo->error_message,
             pfo->error_message);

  phw->sigma_8[index_pk][index_tau] = sigma8;
  phw->sigma_disp[index_pk][index_tau] = sigma_disp;
  phw->sigma_disp_100[index_pk][index_tau] = sigma_disp100;

  /** Initialisation steps for the 1-Halo Power Integral */
  mmin=ppr->mmin_for_p1h_integral/pba->h; //Minimum mass for integration; (unit conversion from  m[Msun/h] to m[Msun]  )
  mmax=ppr->mmax_for_p1h_integral/pba->h; //Maximum mass for integration;

  class_alloc(mass,ppr->nsteps_for_p1h_integral*sizeof(double),pfo->error_message);
  class_alloc(r_real,ppr->nsteps_for_p1h_integral*sizeof(double),pfo->error_message);
  class_alloc(r_virial,ppr->nsteps_for_p1h_integral*sizeof(double),pfo->error_message);
  class_alloc(sigma_r,ppr->nsteps_for_p1h_integral*sizeof(double),pfo->error_message);
  class_alloc(sigmaf_r,ppr->nsteps_for_p1h_integral*sizeof(double),pfo->error_message);
  class_alloc(nu_arr,ppr->nsteps_for_p1h_integral*sizeof(double),pfo->error_message);

  // Linear theory density perturbation threshold for spherical collapse
  // The Einstein de-Sitter value would be delta_c = (3./20.)*pow(12.*_PI_,2./3.) ~1.686
  // Also virialized overdensity
  switch (pfo->hm_version ){
  case hmcode_version_2015:
    delta_c = 1.59+0.0314*log(sigma8); //Mead et al. (2015; arXiv 1505.07833)
    delta_c = delta_c*(1.+0.012299*log10(Omega_m)); //Nakamura & Suto (1997) fitting formula for LCDM models (as in Mead 2016) //BUGFIX :: changed 0.0123 to 0.012299
    delta_c = delta_c*(1.+0.262*fnu); //Mead et al. (2016; arXiv 1602.02154) neutrino addition

    Delta_v=418.*pow(Omega_m, -0.352); //Mead et al. (2015; arXiv 1505.07833)
    Delta_v=Delta_v*(1.+0.916*fnu); //Mead et al. (2016; arXiv 1602.02154) neutrino addition
    break;
  case hmcode_version_2020_baryonic:
  case hmcode_version_2020_unfitted:
  case hmcode_version_2020:
    dv_x = scale_indep_growth; //phg->index_gt_g is already g(a)/a
    dv_y = accumulated_growth * (1.+z_at_tau);

    delta_c = 1.;
    delta_c = delta_c + (-0.0069 - 0.0208 * (1.-dv_x) + 0.0312 * (1.-dv_x) * (1.-dv_x) + 0.0021 * (1.-dv_y))*pow(log10(Om_norad),1);
    delta_c = delta_c + (+0.0001 - 0.0647 * (1.-dv_x) - 0.0417 * (1.-dv_x) * (1.-dv_x) + 0.0646 * (1.-dv_y))*pow(log10(Om_norad),0);

    delta_c = delta_c*((3./20.)*pow(12.*_PI_,2./3.))*(1.-0.041*fnu);

    Delta_v = 1.;
    Delta_v = Delta_v + (-0.79 - 10.17 * (1.-dv_x) + 2.51 * (1.-dv_x) * (1.-dv_x) +  6.51 * (1.-dv_y))*pow(log10(Om_norad),1);
    Delta_v = Delta_v + (-1.89 +  0.38 * (1.-dv_x) + 18.8 * (1.-dv_x) * (1.-dv_x) - 15.87 * (1.-dv_y))*pow(log10(Om_norad),2);

    Delta_v = Delta_v*(18.*_PI_*_PI_)*(1.+0.763*fnu);
    break;
    //  delta_c = 1.686;
    //  Delta_v = 200.;
    //break;
  }

  // mass or radius fraction respectively
  fraction = pow(0.01, 1./3.);

  /* Fill the arrays needed for the P1H Integral: mass, r_real, r_virial, nu_arr, sigma_r, sigmaf_r
   * The P1H Integral is an integral over nu=delta_c/sigma(M), where M is connected to R via R=(3M)/(4*pi*rho_m).
   * The Integrand is M*Window^2{nu(M)*k, Rv(M), c(M)}*f(nu) with the window being the fouriertransformed
   * NFW profile, Rv = R/Delta_v^(1/3) and Sheth-Thormen halo mass function f.
   * The halo concentration-mass-relation c(M) will be found later.  */

  for (index_mass=0;index_mass<ppr->nsteps_for_p1h_integral;index_mass++){

    m = exp(log(mmin)+log(mmax/mmin)*(index_mass)/(ppr->nsteps_for_p1h_integral-1));
    r = pow((3.*m/(4.*_PI_*rho_crit_today_in_msun_mpc3*Omega0_m)), (1./3.));
    mass[index_mass] = m;
    r_real[index_mass] = r;
    r_virial[index_mass] = r_real[index_mass]/pow(Delta_v, 1./3.);

    class_call(array_interpolate_spline(phw->rtab,
                                        nsig,
                                        phw->stab,
                                        phw->ddstab,
                                        1,
                                        r,
                                        &last_index,
                                        &sig,
                                        1,
                                        pfo->error_message),
               pfo->error_message, pfo->error_message);

    class_call(array_interpolate_spline(phw->rtab,
                                        nsig,
                                        phw->stab,
                                        phw->ddstab,
                                        1,
                                        r*fraction,
                                        &last_index,
                                        &sigf,
                                        1,
                                        pfo->error_message),
               pfo->error_message, pfo->error_message);

    nu=delta_c/sig;
    sigma_r[index_mass] = sig;
    sigmaf_r[index_mass] = sigf;
    nu_arr[index_mass] = nu;
  }

  /** find nonlinear scales k_nl and r_nl and the effective spectral index n_eff */
  nu_nl = 1.;
  nu_min = nu_arr[0];

  /* stop calculating the nonlinear correction if the nonlinear scale is not reached in the table: */
  if (nu_min > nu_nl) {
    if (pfo->fourier_verbose>0) fprintf(stdout, " -> [WARNING:] the minimum mass in the mass-table is too large to find the nonlinear scale at this redshift.\n   Decrease mmin_for_p1h_integral\n");
    * nl_corr_not_computable_at_this_k = _TRUE_;
    free(mass);
    free(r_real);
    free(r_virial);
    free(sigma_r);
    free(sigmaf_r);
    free(nu_arr);
    return _SUCCESS_;
  }

  /* make a first guess for the nonlinear scale */
  class_call(array_interpolate_two_arrays_one_column(
                                                     nu_arr,
                                                     r_real,
                                                     1,
                                                     0,
                                                     ppr->nsteps_for_p1h_integral,
                                                     nu_nl,
                                                     &r_nl,
                                                     pfo->error_message),
             pfo->error_message, pfo->error_message);

  class_call(array_search_bisect(ppr->nsteps_for_p1h_integral,nu_arr,nu_nl,&index_nl,pfo->error_message), pfo->error_message, pfo->error_message);

  r1 = r_real[index_nl-1];
  r2 = r_real[index_nl+2];
  nu1 = nu_arr[index_nl-1];
  nu2 = nu_arr[index_nl+2];

  /* As an initial guess, we adopt the linear interpolation */
  r_nl = (nu_nl-nu_arr[index_nl])/(nu_arr[index_nl+1]-nu_arr[index_nl])*(r_real[index_nl+1]-r_real[index_nl]) + r_real[index_nl];
  // for debugging: (if it happens that r_nl is not between r1 and r2, which should never be the case)
  //fprintf(stdout, "r1=%e nu1=%e r2=%e nu2=%e\n", r1, nu_arr[index_nl-1], r2, nu_arr[index_nl+2]);


  // do iteration between r1 and r2 to find the precise value of r_nl --> TODO :: speed up by semi-newtonian iterations
  counter = 0;
  do {

    class_call(fourier_sigmas(pfo,
                              r_nl,
                              lnpk_l[index_pk_cb],ddlnpk_l[index_pk_cb],
                              pfo->k_size_extra,
                              ppr->sigma_k_per_decade,
                              out_sigma,
                              &sigma_nl),
               pfo->error_message, pfo->error_message);

    diff = sigma_nl - delta_c;

    if (diff > ppr->hmcode_tol_sigma){
      r1 = r_nl;
    }
    else if (diff < -ppr->hmcode_tol_sigma) {
      r2 = r_nl;
    }
    r_nl = (r1+r2)/2.;
    counter ++;

    class_test(counter > _MAX_IT_,
               pfo->error_message,
               "could not converge within maximum allowed number of iterations");

  } while (fabs(diff) > ppr->hmcode_tol_sigma);

  if (pfo->fourier_verbose>5){
    fprintf(stdout, "number of iterations for r_nl at z = %e: %d\n", z_at_tau, counter);
  }
  *k_nl = 1./r_nl;

  if (*k_nl > pfo->k[pfo->k_size-1]) {
    * nl_corr_not_computable_at_this_k = _TRUE_;
    free(mass);
    free(r_real);
    free(r_virial);
    free(sigma_r);
    free(sigmaf_r);
    free(nu_arr);
    return _SUCCESS_;
  }
  else {
    * nl_corr_not_computable_at_this_k = _FALSE_;
  }

  /* call sigma_prime function at r_nl to find the effective spectral index n_eff */

  class_call(fourier_sigmas(pfo,
                            r_nl,
                            lnpk_l[index_pk_cb],ddlnpk_l[index_pk_cb],
                            pfo->k_size_extra,
                            ppr->sigma_k_per_decade,
                            out_sigma_prime,
                            &sigma_prime),
             pfo->error_message,
             pfo->error_message);

  dlnsigdlnR = r_nl*pow(sigma_nl, -2)*sigma_prime;
  n_eff = -3.- dlnsigdlnR;

  phw->sigma_prime[index_pk][index_tau] = sigma_prime;

  /** Calculate halo concentration-mass relation conc(mass) (Bullock et al. 2001) */
  class_alloc(conc,ppr->nsteps_for_p1h_integral*sizeof(double),pfo->error_message);

  switch (pfo->hm_version ){
  case hmcode_version_2015:
    DEcorr = phw->dark_energy_correction;
    break;
  case hmcode_version_2020:
    //  DEcorr = phw->dark_energy_correction;
    //break;
  case hmcode_version_2020_baryonic:
  case hmcode_version_2020_unfitted:
    /* Dark Energy Dolag correction */
    class_call(hmcode_growint(ppr,pba,pfo,1./(1.+z_at_tau),-1.,0.,&g_lcdm),
               pfo->error_message, pfo->error_message);
    //There are different ways to define this, but the uncommented one is the one consistent with HMcode 2020
    /*class_call(background_w_fld(pba,1.,&w0,&dw_over_da_fld,&integral_fld),
      pba->error_message,
      pfo->error_message);

      class_call(hmcode_growint(ppr,pba,pfo,1./(1.+z_at_tau),w0,dw_over_da_fld*(-1.),&g_wcdm),
      pfo->error_message,
      pfo->error_message);
      DEcorr = phw->dark_energy_correction*g_lcdm/g_wcdm;*/
    DEcorr = phw->dark_energy_correction*g_lcdm/growth;
    break;
  }
  switch (pfo->hm_version ){
  case hmcode_version_2015:
    Abary = pfo->c_min;
    alpha = 3.24 * pow(1.85, n_eff);
    break;
  case hmcode_version_2020:
    Abary = ( 5.1958429 * pow(sigma8_cb,0.) + 0. ) * pow(10,z_at_tau*0.);
    alpha = 1.8751465*pow(1.6029913,n_eff);
    break;
  case hmcode_version_2020_baryonic:
    Abary = ( (-0.496*(pfo->log10T_heat_hmcode-7.8)+3.44) * pow(sigma8_cb,0.) + 0.) * pow(10,z_at_tau*(-0.0371*(pfo->log10T_heat_hmcode-7.8)-0.0671));
    alpha = 1.0;
    break;
  case hmcode_version_2020_unfitted:
    Abary = 4.;
    alpha = 1.;
    break;
  }

  /* Fill concentration array */
  for (index_mass=0;index_mass<ppr->nsteps_for_p1h_integral;index_mass++){
    //find growth rate at formation
    g_form = delta_c*growth/sigmaf_r[index_mass];
    if (g_form > 1.) g_form = 1.;

    /*class_call(array_interpolate_two_arrays_one_column(
      phw->growtable,
      phw->ztable,
      1,
      0,
      ng,
      g_form,
      &z_form,
      pfo->error_message),
      pfo->error_message, pfo->error_message);*/

    class_call(array_spline_hunt(phg->normgrowth_table,
                                 phg->a_size,
                                 g_form,
                                 &last_index,
                                 &h,
                                 &a,
                                 &b,
                                 pfo->error_message),
               pfo->error_message,
               pfo->error_message);

    a_form = array_spline_eval(phg->a_table,phg->dda_table,last_index,last_index+1,h,a,b);
    z_form = 1./a_form-1.; //HM2020

    if (z_form < z_at_tau){
      conc[index_mass] = Abary;
    } else {
      conc[index_mass] = Abary*(1.+z_form)/(1.+z_at_tau);
    }
    if(z_at_tau < pfo->z_infinity){
      conc[index_mass] = conc[index_mass]*DEcorr;
    }
  }

  /* HMcode parameters for non-linear corerction */
  switch (pfo->hm_version ){
  case hmcode_version_2015:
    k_star=0.584/sigma_disp;   // Damping wavenumber of the 1-halo term at very large scales;
    eta = pfo->eta_0 - 0.3*sigma8; // halo bloating parameter
    fdamp = 0.0095*pow(sigma_disp100*pba->h, 1.37);
    break;
  case hmcode_version_2020:
    //In CAMB, these are cb sigma8
    eta = 0.1281210*pow(sigma8_cb, -0.3643963);
    k_star=0.0561778*pow(sigma8_cb,-1.0131066)*pba->h;
    fdamp = 0.2695822 * pow(sigma8_cb, 0.9403087);
    kdamp = 0.0569871*pow(sigma8_cb, -1.0890162) *pba->h;
    break;
  case hmcode_version_2020_baryonic:
  case hmcode_version_2020_unfitted:
    eta = 0.;
    k_star = 0.;
    kdamp = 0.;
    fdamp = 0.;
    break;
  }

  // Damping factor for 2-halo term
  if (fdamp<1.e-3) fdamp=1.e-3;
  if (fdamp>0.99)  fdamp=0.99;
  if ( pfo->hm_version == hmcode_version_2020_baryonic || pfo->hm_version == hmcode_version_2020_unfitted ){
    fdamp = 0.;
  }

  /* the 1h integral contains the halo mass function proportional to exp(-nu^2).
   * To save time, the integration loop cuts, when nu exceeds a large value,
   * where the integrand is 0 anyhow. This cut index is found here. */
  nu_cut = 10.;
  // for now we disable the automatic search for nu_cut, since it creates tiny deviations from HMcode 2020 // TBC: why?
  //if (nu_cut < nu_arr[ppr->nsteps_for_p1h_integral-1]){
  if (nu_cut < nu_arr[ppr->nsteps_for_p1h_integral-1] && _FALSE_){
    class_call(array_search_bisect(ppr->nsteps_for_p1h_integral,nu_arr,nu_cut,&index_cut,pfo->error_message), pfo->error_message, pfo->error_message);
  }
  else {
    index_cut = ppr->nsteps_for_p1h_integral;
  }

  i=0;
  index_nu=i;
  i++;
  index_y=i;
  i++;
  index_ddy=i;
  i++;
  index_ncol=i;

  class_setup_parallel();

  if (pfo->fourier_verbose>2) {
    num_threads = task_system.GetNumThreads();
    gettimeofday(&begin, 0);
  }

  for (index_k = 0; index_k < pfo->k_size; index_k++){

    class_run_parallel( \
                       =,

                       double declare_list_of_variables_inside_parallel_region(k,logk,pk_lin,window_nfw,gst,pk_1h,fac,pk_2h,pk_wig,fac_dewiggle);
                       double declare_list_of_variables_inside_parallel_region(sbar_bar,sbarz_bar,mbar_bar,mbarz_bar,mb,sb,ratio_m,fb,fc,fs);
                       double DMONLY_halo_mass_fraction;
                       int index_mass_p;
                       double * p1h_integrand;
                       int last_index_p;

                       k = pfo->k[index_k];
                       logk = log(k);

                       class_alloc(p1h_integrand,index_cut*index_ncol*sizeof(double),pfo->error_message);

                       pk_lin = exp(lnpk_l[index_pk][index_k])*pow(k,3)*anorm; //convert P_k to Delta_k^2

                       for (index_mass_p=0; index_mass_p<index_cut; index_mass_p++){ //Calculates the integrand for the ph1 integral at all nu values
                         //get the nu^eta-value of the window
                         class_call(hmcode_window_nfw(pfo,
                                                      pow(nu_arr[index_mass_p], eta)*k,
                                                      r_virial[index_mass_p],
                                                      conc[index_mass_p], //(rv/rs)
                                                      &window_nfw),
                                    pfo->error_message, pfo->error_message);

                         //get the value of the halo mass function
                         //equivalent to g_nu function (or g_ST)
                         class_call(hmcode_halomassfunction(nu_arr[index_mass_p],
                                                            &gst),
                                    pfo->error_message, pfo->error_message);

                         p1h_integrand[index_mass_p*index_ncol+index_nu] = nu_arr[index_mass_p];

                         switch (pfo->hm_version ){
                         case hmcode_version_2015:
                           break;
                         case hmcode_version_2020_unfitted:
                         case hmcode_version_2020:
                           window_nfw *= (1.-fnu);
                           break;
                         case hmcode_version_2020_baryonic:
                           sbar_bar = -0.0030*(pfo->log10T_heat_hmcode-7.8)+0.0201;
                           sbarz_bar = 0.0224*(pfo->log10T_heat_hmcode-7.8)+0.409;
                           sb = MIN(sbar_bar*pow(10,z_at_tau*sbarz_bar),pba->Omega0_b/pba->Omega0_m);

                           mbar_bar = pow(10.,1.81*(pfo->log10T_heat_hmcode-7.8)+13.87);
                           mbarz_bar = 0.195*(pfo->log10T_heat_hmcode-7.8)-0.108;

                           mb = mbar_bar*pow(10,mbarz_bar*z_at_tau);
                           ratio_m = pow(pba->h*mass[index_mass_p]/mb,2.); //hmod%nbar
                           fb = pba->Omega0_b/pba->Omega0_m;
                           fc = pba->Omega0_cdm/pba->Omega0_m;
                           fs = sb;
                           DMONLY_halo_mass_fraction = fc+(fb-fs)*ratio_m/(1.+ratio_m);

                           window_nfw = window_nfw*DMONLY_halo_mass_fraction; //Account for gas expulsion

                           window_nfw = window_nfw+sb; //Account for star formation
                           break;
                         }

                         p1h_integrand[index_mass_p*index_ncol+index_y] = mass[index_mass_p]*gst*pow(window_nfw, 2.);
                       }

                       class_call(array_spline(p1h_integrand,
                                               index_ncol,
                                               index_cut,
                                               index_nu,
                                               index_y,
                                               index_ddy,
                                               _SPLINE_EST_DERIV_,
                                               pfo->error_message),
                                  pfo->error_message,
                                  pfo->error_message);

                       class_call(array_integrate_all_trapzd_or_spline(
                                                                       p1h_integrand,
                                                                       index_ncol,
                                                                       index_cut,
                                                                       index_cut-1, //ranges from 0 to n-1
                                                                       index_nu,
                                                                       index_y,
                                                                       index_ddy,
                                                                       &pk_1h,
                                                                       pfo->error_message),
                                  pfo->error_message,
                                  pfo->error_message);

                       switch (pfo->hm_version ){
                       case hmcode_version_2015:
                         if (pow(k/k_star, 2)>7.){
                           fac = 1.;     //prevents problems if (k/k*)^2 is large
                         }
                         else{
                           fac = 1.-exp(-pow((k/k_star), 2.));
                         }
                         break;
                       case hmcode_version_2020:
                         if ( k_star == 0 ){
                           fac = 1.;
                         } else {
                           fac = pow(k/k_star,4)/(1.+pow(k/k_star,4));
                         }
                         break;
                       case hmcode_version_2020_unfitted:
                       case hmcode_version_2020_baryonic:
                         fac = 1.;
                         break;
                       }

                       pk_1h = pk_1h*anorm*pow(k,3)*fac/(rho_crit_today_in_msun_mpc3*Omega0_m);  // dimensionless power

                       switch (pfo->hm_version ){
                       case hmcode_version_2015:
                         pk_2h = pk_lin;
                         if(fdamp>0){
                           pk_2h = pk_2h*(1.-fdamp*pow(tanh(k*sigma_disp/sqrt(fdamp)), 2.));
                         }
                         break;
                       case hmcode_version_2020:
                       case hmcode_version_2020_unfitted:
                       case hmcode_version_2020_baryonic:
                         if(logk<phw->lnk_wiggle[0] || logk>phw->lnk_wiggle[pfo->nk_wiggle-1]){
                           pk_wig = 0.;
                         }
                         else{
                           class_call(array_interpolate_spline(phw->lnk_wiggle,
                                                               pfo->nk_wiggle,
                                                               phw->pk_wiggle,
                                                               phw->ddpk_wiggle,
                                                               1,
                                                               log(k),
                                                               &last_index_p,
                                                               &pk_wig,
                                                               1,
                                                               pfo->error_message),
                                      pfo->error_message,
                                      pfo->error_message);
                         }
                         fac_dewiggle = exp(-(k*k*sigma_disp*sigma_disp));
                         pk_2h = pk_lin + (fac_dewiggle-1.)*pk_wig;

                         if(fdamp>0){
                           pk_2h = pk_2h*(1.-fdamp*(pow(k/kdamp,2.8534197))/(1.+pow(k/kdamp,2.8534197)));
                         }
                         break;
                       }

                       //p_hm function: Combine pk_1h and pk_2h
                       class_test(pk_2h < 0. || pk_1h < 0.,pfo->error_message,"The 2 halo or 1 halo term is negative for HMcode, and the 'safe_negative' option is disabled. Aborting");
                       pk_nl[index_k] = pow((pow(pk_1h, alpha) + pow(pk_2h, alpha)), (1./alpha))/pow(k,3)/anorm;
                       free(p1h_integrand);

                       return _SUCCESS_;
                        );
  }

  class_finish_parallel();

  if (pfo->fourier_verbose>2) {
    gettimeofday(&end, 0);
    seconds = end.tv_sec - begin.tv_sec;
    microseconds = end.tv_usec - begin.tv_usec;
    elapsed = seconds + microseconds*1e-6;

    fprintf(stderr,"In %s: time spent in parallel region (loop over k's at tau[%d]= %f) = %.6f s using %d threads\n",
            __func__,
            index_tau,
            tau,
            elapsed,
            num_threads
            );
  }

  // print parameter values
  if ((pfo->fourier_verbose > 2 && tau==pba->conformal_age) || pfo->fourier_verbose > 3){
    fprintf(stdout, " -> Parameters at redshift z = %e:\n", z_at_tau);
    fprintf(stdout, "    fnu:		%e\n", fnu);
    fprintf(stdout, "    sigd [Mpc/h]:	%e\n", sigma_disp*pba->h);
    //fprintf(stdout, "    sigd100 [Mpc/h]:    %e\n", sigma_disp100*pba->h);
    fprintf(stdout, "    sigma8:		%e\n", sigma8);
    fprintf(stdout, "    nu min:		%e\n", nu_arr[0]);
    fprintf(stdout, "    nu max:		%e\n", nu_arr[ppr->nsteps_for_p1h_integral-1]);
    fprintf(stdout, "    r_v min [Mpc/h]:    %e\n", r_virial[0]*pba->h);
    fprintf(stdout, "    r_v max [Mpc/h]:    %e\n", r_virial[ppr->nsteps_for_p1h_integral-1]*pba->h);
    fprintf(stdout, "    r_nl [Mpc/h]:	%e\n", r_nl*pba->h);
    fprintf(stdout, "    k_nl [h/Mpc]:	%e\n", *k_nl/pba->h);
    fprintf(stdout, "    sigma_nl:		%e\n", sigma_nl/delta_c);
    fprintf(stdout, "    neff:		%e\n", n_eff);
    fprintf(stdout, "    c min:		%e\n", conc[ppr->nsteps_for_p1h_integral-1]);
    fprintf(stdout, "    c max:		%e\n", conc[0]);
    fprintf(stdout, "    Dv:			%e\n", Delta_v);
    fprintf(stdout, "    dc:			%e\n", delta_c);
    fprintf(stdout, "    eta:		%e\n", eta);
    fprintf(stdout, "    k*:			%e\n", k_star/pba->h);
    fprintf(stdout, "    Abary:		%e\n", Abary);
    fprintf(stdout, "    fdamp:		%e\n", fdamp);
    fprintf(stdout, "    alpha:		%e\n", alpha);
    fprintf(stdout, "    ksize, kmin, kmax:   %d, %e, %e\n", pfo->k_size, pfo->k[0]/pba->h, pfo->k[pfo->k_size-1]/pba->h);

  }

  free(conc);
  free(mass);
  free(r_real);
  free(r_virial);
  free(sigma_r);
  free(sigmaf_r);
  free(nu_arr);

  return _SUCCESS_;
}

/**
 * allocate and fill arrays of hmcode workspace
 *
 * @param ppr         Input: pointer to precision structure
 * @param pba         Input: pointer to background structure
 * @param pfo         Input: pointer to fourier structure
 * @param phw         Output: pointer to hmcode workspace
 * @return the error status
 */

int hmcode_workspace_init(
                          struct precision *ppr,
                          struct background *pba,
                          struct fourier *pfo,
                          struct hmcode_workspace *phw
                          ){

  int ng;
  int index_pk;

  /** - allocate arrays of the hmcode workspace */

  class_alloc(phw->rtab,ppr->n_hmcode_tables*sizeof(double),pfo->error_message);
  class_alloc(phw->stab,ppr->n_hmcode_tables*sizeof(double),pfo->error_message);
  class_alloc(phw->ddstab,ppr->n_hmcode_tables*sizeof(double),pfo->error_message);

  ng = ppr->n_hmcode_tables;

  class_alloc(phw->growtable,ng*sizeof(double),pfo->error_message);
  class_alloc(phw->ztable,ng*sizeof(double),pfo->error_message);
  class_alloc(phw->tautable,ng*sizeof(double),pfo->error_message);

  class_alloc(phw->sigma_8,pfo->pk_size*sizeof(double *),pfo->error_message);
  class_alloc(phw->sigma_disp,pfo->pk_size*sizeof(double *),pfo->error_message);
  class_alloc(phw->sigma_disp_100,pfo->pk_size*sizeof(double *),pfo->error_message);
  class_alloc(phw->sigma_prime,pfo->pk_size*sizeof(double *),pfo->error_message);

  for (index_pk=0; index_pk<pfo->pk_size; index_pk++){
    class_alloc(phw->sigma_8[index_pk],pfo->tau_size*sizeof(double),pfo->error_message);
    class_alloc(phw->sigma_disp[index_pk],pfo->tau_size*sizeof(double),pfo->error_message);
    class_alloc(phw->sigma_disp_100[index_pk],pfo->tau_size*sizeof(double),pfo->error_message);
    class_alloc(phw->sigma_prime[index_pk],pfo->tau_size*sizeof(double),pfo->error_message);
  }

  /** - fill table with scale independent growth factor */

  class_call(hmcode_fill_growtab(ppr,pba,pfo,phw),
             pfo->error_message,
             pfo->error_message);

  class_alloc(phw->pk_wiggle, pfo->nk_wiggle*sizeof(double), pfo->error_message);
  class_alloc(phw->ddpk_wiggle, pfo->nk_wiggle*sizeof(double), pfo->error_message);
  class_alloc(phw->lnk_wiggle, pfo->nk_wiggle*sizeof(double), pfo->error_message);

  return _SUCCESS_;
}

/**
 * Deallocate arrays in the hmcode workspace
 *
 * @param pfo Input: pointer to fourier structure
 * @param phw Input: pointer to hmcode workspace
 * @return the error status
 */

int hmcode_workspace_free(
                          struct fourier *pfo,
                          struct hmcode_workspace *phw
                          ) {

  int index_pk;

  free(phw->rtab);
  free(phw->stab);
  free(phw->ddstab);

  free(phw->growtable);
  free(phw->ztable);
  free(phw->tautable);

  for (index_pk=0; index_pk<pfo->pk_size; index_pk++){
    free(phw->sigma_8[index_pk]);
    free(phw->sigma_disp[index_pk]);
    free(phw->sigma_disp_100[index_pk]);
    free(phw->sigma_prime[index_pk]);
  }

  free(phw->sigma_8);
  free(phw->sigma_disp);
  free(phw->sigma_disp_100);
  free(phw->sigma_prime);

  free(phw->pk_wiggle);
  free(phw->ddpk_wiggle);
  free(phw->lnk_wiggle);

  return _SUCCESS_;
}

/**
 * Set the HMcode dark energy correction (if w is not -1)
 *
 * @param ppr         Input: pointer to precision structure
 * @param pba         Input: pointer to background structure
 * @param pfo         Input: pointer to fourier structure
 * @param phw         Output: pointer to hmcode workspace
 * @return the error status
 */

int hmcode_dark_energy_correction(
                                  struct precision *ppr,
                                  struct background *pba,
                                  struct fourier *pfo,
                                  struct hmcode_workspace *phw
                                  ) {

  int last_index;
  double *pvecback;
  double tau_growth;
  double g_lcdm,g_wcdm;
  double w0,dw_over_da_fld,integral_fld;
  double h,a,b;
  double growth;

  /** - if there is dynamical Dark Energy (w is not -1) modeled as a fluid */

  if (pba->has_fld==_TRUE_){

    class_alloc(pvecback,pba->bg_size*sizeof(double),pfo->error_message);

    class_call(background_tau_of_z(pba,pfo->z_infinity,&tau_growth),
               pba->error_message,
               pfo->error_message);

    class_call(background_at_tau(pba,tau_growth,long_info,inter_normal,&last_index,pvecback),
               pba->error_message,
               pfo->error_message);

    class_call(background_w_fld(pba,1.,&w0,&dw_over_da_fld,&integral_fld),
               pba->error_message,
               pfo->error_message);

    class_call(hmcode_growint(ppr,pba,pfo,1./(1.+pfo->z_infinity),-1.,0.,&g_lcdm),
               pfo->error_message, pfo->error_message);

    class_call(hmcode_growint(ppr,pba,pfo,1./(1.+pfo->z_infinity),w0,dw_over_da_fld*(-1.),&g_wcdm),
               pfo->error_message,
               pfo->error_message);

    free(pvecback);

    switch( pfo->hm_version ){
    case hmcode_version_2015:
      phw->dark_energy_correction = pow(g_wcdm/g_lcdm, 1.5);
      break;
    case hmcode_version_2020:
    case hmcode_version_2020_unfitted:
    case hmcode_version_2020_baryonic:
      //phw->dark_energy_correction = g_wcdm/g_lcdm; //Here, we are missing the redshift dependence, which we include only within the respective function
      //hmcode_norad(phw->phg,1./(pfo->z_infinity+1.));
      class_call(array_spline_hunt(phw->phg->a_table,
                                   phw->phg->a_size,
                                   1./(pfo->z_infinity+1.),
                                   &last_index,
                                   &h,
                                   &a,
                                   &b,
                                   pfo->error_message),
                 pfo->error_message,
                 pfo->error_message);

      growth = array_spline_eval(phw->phg->normgrowth_table,phw->phg->ddnormgrowth_table,last_index,last_index+1,h,a,b);
      phw->dark_energy_correction = growth/g_lcdm;
      break;
    }
  }

  /** - otherwise, we assume no dynamical Dark Energy (w is -1) */

  else {
    phw->dark_energy_correction = 1.;
  }

  return _SUCCESS_;
}

/**
 * Set the HMcode baryonic feedback parameters according to the chosen feedback model
 *
 * @param pfo   Output: pointer to fourier structure
 * @return the error status
 */

int hmcode_baryonic_feedback(
                             struct fourier *pfo
                             ) {

  switch (pfo->feedback) {

  case hmcode_emu_dmonly:
    {
      pfo->eta_0 = 0.603;
      pfo->c_min = 3.13;
      break;
    }

  case hmcode_owls_dmonly:
    {
      pfo->eta_0 = 0.64;
      pfo->c_min = 3.43;
      break;
    }

  case hmcode_owls_ref:
    {
      pfo->eta_0 = 0.68;
      pfo->c_min = 3.91;
      break;
    }

  case hmcode_owls_agn:
    {
      pfo->eta_0 = 0.76;
      pfo->c_min = 2.32;
      break;
    }

  case hmcode_owls_dblim:
    {
      pfo->eta_0 = 0.70;
      pfo->c_min = 3.01;
      break;
    }

  case hmcode_user_defined:
    {
      /* eta_0 and c_min already passed in input */
      break;
    }
  }
  return _SUCCESS_;
}

/**
 * Function that fills phw->rtab, phw->stab and phw->ddstab with (r,
 * sigma, ddsigma) logarithmically spaced in r.  Called by
 * fourier_init at for all tau to account for scale-dependant growth
 * before hmcode is called
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param ppt        Input: pointer to perturbation structure
 * @param ppm        Input: pointer to primordial structure
 * @param pfo        Input: pointer to fourier structure
 * @param index_tau  Input: index of tau, at which to compute the nl correction
 * @param index_pk   Input: index of the pk type, either index_m or index_cb
 * @param lnpk_l     Input: logarithm of the linear power spectrum for either index_m or index_cb
 * @param ddlnpk_l   Input: spline of the logarithm of the linear power spectrum for either index_m or index_cb
 * @param phw        Output: pointer to hmcode workspace
 * @return the error status
 */

int hmcode_fill_sigtab(
                       struct precision *ppr,
                       struct background *pba,
                       struct perturbations *ppt,
                       struct primordial *ppm,
                       struct fourier *pfo,
                       int index_tau,
                       int index_pk,
                       double **lnpk_l,
                       double **ddlnpk_l,
                       struct hmcode_workspace *phw
                       ) {

  double r;
  double rmin, rmax;
  double sig;
  double *sigtab;
  double *lpk, *ddlpk;
  int i, index_r, index_sig, index_ddsig, index_n, nsig;
  double cbcorr, z, growth;
  int last_index = 0;
  double a,b,h;
  double *head;
  double cbnu_ratio;
  int index_k;
  int index_pk_default;

  struct hmcode_growth *phg;

  phg = phw->phg;

  rmin = ppr->rmin_for_sigtab/pba->h;
  rmax = ppr->rmax_for_sigtab/pba->h;
  nsig = ppr->n_hmcode_tables;

  i=0;
  index_r=i;
  i++;
  index_sig=i;
  i++;
  index_ddsig=i;
  i++;
  index_n=i;

  class_alloc((sigtab),(nsig*index_n*sizeof(double)),pfo->error_message);
  class_alloc(lpk, pfo->k_size_extra*sizeof(double), pfo->error_message);
  class_alloc(ddlpk, pfo->k_size_extra*sizeof(double), pfo->error_message);
  for (index_k=0; index_k<pfo->k_size_extra; index_k++){
    cbcorr = 0.;
    if ((pfo->has_pk_m == _TRUE_ && pfo->has_pk_cb == _TRUE_ && index_pk == pfo->index_pk_m)
        && (pfo->hm_version == hmcode_version_2020 || pfo->hm_version == hmcode_version_2020_baryonic)){

      class_call(background_z_of_tau(pba,pfo->tau[index_tau],&z),
                 pba->error_message,
                 phg->error_message);

      class_call(array_spline_hunt(phg->a_table,
                                   phg->a_size,
                                   1./(z+1.),
                                   &last_index,
                                   &h,
                                   &a,
                                   &b,
                                   pfo->error_message),
                 pfo->error_message,
                 pfo->error_message);

      head = phg->growth_table+last_index*phg->gt_size;
      growth = array_spline_eval(head+phg->index_gt_g,head+phg->index_gt_ddg,0,phg->gt_size,h,a,b);
      // growth = g(a)/a = 1 in the EdS case, otherwise it is found in the growth table (without radiation, though)
      class_call(hmcode_cbnu_ratio(pfo->k[index_k],
                                   z,
                                   phg->fnu,
                                   phg->om_m,
                                   phg->Tcmb,
                                   growth/(1.+z),
                                   &cbnu_ratio),
                 pfo->error_message,
                 pfo->error_message);

      cbcorr = 2. * log(cbnu_ratio);
      lpk[index_k] = lnpk_l[index_pk][index_k]+cbcorr;
    }
    else{
      if (pfo->has_pk_cb == _TRUE_)
        index_pk_default = pfo->index_pk_cb;
      else
        index_pk_default = pfo->index_pk_m;
      lpk[index_k] = lnpk_l[index_pk_default][index_k];//Always take the P(k)_cb if we are a) in the cb case or b) not in HMcode 2020
    }
  }

  class_call(array_spline_table_columns(
                                        pfo->ln_k,
                                        pfo->k_size_extra,
                                        lpk,
                                        1,
                                        ddlpk,
                                        _SPLINE_EST_DERIV_,
                                        pfo->error_message),
             pfo->error_message,
             pfo->error_message);

  for (i=0;i<nsig;i++){
    r=exp(log(rmin)+log(rmax/rmin)*i/(nsig-1));

    class_call(fourier_sigmas(pfo,
                              r,
                              lpk,
                              ddlpk,
                              pfo->k_size_extra,
                              ppr->sigma_k_per_decade,
                              out_sigma,
                              &sig),
               pfo->error_message,
               pfo->error_message);

    sigtab[i*index_n+index_r]=r;
    sigtab[i*index_n+index_sig]=sig;
  }

  class_call(array_spline(sigtab,
						  index_n,
						  nsig,
						  index_r,
						  index_sig,
						  index_ddsig,
						  _SPLINE_EST_DERIV_,
						  pfo->error_message),
             pfo->error_message,
             pfo->error_message);

  if (index_tau == pfo->tau_size-1){
    for (i=0;i<nsig;i++){
      phw->rtab[i] = sigtab[i*index_n+index_r];
      phw->stab[i] = sigtab[i*index_n+index_sig];
      phw->ddstab[i] = sigtab[i*index_n+index_ddsig];
    }
  }
  else{
    for (i=0;i<nsig;i++){
      phw->stab[i] = sigtab[i*index_n+index_sig];
      phw->ddstab[i] = sigtab[i*index_n+index_ddsig];
    }
  }

  free(lpk);
  free(ddlpk);
  free(sigtab);

  return _SUCCESS_;
}

/**
 * Function that fills phw->tautable and phw->growtable with (tau, D(tau))
 * linearly spaced in scalefactor a.
 * Called by fourier_init at before the loop over tau
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure (will provide the scale independent growth factor)
 * @param pfo Input/Output: pointer to fourier structure
 * @param phw Output: pointer to hmcode workspace
 * @return the error status
 */

int hmcode_fill_growtab(
                        struct precision *ppr,
                        struct background *pba,
                        struct fourier *pfo,
                        struct hmcode_workspace *phw
                        ){

  double z, ainit, amax, scalefactor, tau_growth;
  int index_scalefactor, last_index, ng;
  double *pvecback;

  ng = ppr->n_hmcode_tables;
  ainit = ppr->ainit_for_growtab;
  amax = ppr->amax_for_growtab;

  last_index = 0;

  class_alloc(pvecback,pba->bg_size*sizeof(double),pfo->error_message);

  for (index_scalefactor=0;index_scalefactor<ng;index_scalefactor++){
    scalefactor = ainit+(amax-ainit)*(index_scalefactor)/(ng-1);
    z = 1./scalefactor-1.;

    phw->ztable[index_scalefactor] = z;

    class_call(background_tau_of_z(pba, z, &tau_growth),
               pba->error_message,
               pfo->error_message);

    phw->tautable[index_scalefactor] = tau_growth;

    class_call(background_at_tau(pba,tau_growth,long_info,inter_normal,&last_index,pvecback),
               pba->error_message,
               pfo->error_message);

    phw->growtable[index_scalefactor] = pvecback[pba->index_bg_D];

  }

  free(pvecback);

  return _SUCCESS_;
}

/**
 * This function finds the scale independent growth factor by
 * integrating the approximate relation d(lnD)/d(lna) =
 * Omega_m(z)^gamma by Linder & Cahn 2007
 *
 * @param ppr    Input: pointer to precision structure
 * @param pba    Input: pointer to background structure
 * @param pfo    Input: pointer to fourier structure
 * @param a      Input: scalefactor
 * @param w0     Input: dark energy equation of state today
 * @param wa     Input: dark energy equation of state varying with a: w=w0+(1-a)wa
 * @param growth Output: scale independent growth factor at a
 * @return the error status
 */

int hmcode_growint(
                   struct precision *ppr,
                   struct background *pba,
                   struct fourier *pfo,
                   double a,
                   double w0,
                   double wa,
                   double *growth
                   ){

  double z, ainit, amax, scalefactor, gamma, X_de, Hubble2, Omega_m;
  int i, index_scalefactor, index_a, index_growth, index_ddgrowth, index_gcol, ng; // index_scalefactor is a running index while index_a is a column index
  double *pvecback;
  double *integrand;

  ng = 1024; // number of growth values (stepsize of the integral), should not be hardcoded and replaced by a precision parameter
  ainit = a;
  amax = 1.;

  i=0;
  index_a = i;
  i++;
  index_growth = i;
  i++;
  index_ddgrowth = i;
  i++;
  index_gcol = i;

  class_alloc(integrand,ng*index_gcol*sizeof(double),pfo->error_message);
  class_alloc(pvecback,pba->bg_size*sizeof(double),pfo->error_message);

  if (ainit == amax) {
    *growth = 1.;
  }
  else {

    for (index_scalefactor=0;index_scalefactor<ng;index_scalefactor++){

      scalefactor = ainit+(amax-ainit)*(index_scalefactor)/(ng-1);
      z = 1./scalefactor-1.;

      /* This will compute Omega_m(z) for the input values of w0 and wa, to let the user compare the wCDM and LCDM cases. This is why we cannot extract Omega_m(z) fromn the background module in this place. */
      X_de = pow(scalefactor, -3.*(1.+w0+wa))*exp(-3.*wa*(1.-scalefactor));
      Hubble2 = (pba->Omega0_m*pow((1.+z), 3.) + pba->Omega0_k*pow((1.+z), 2.) + pba->Omega0_de*X_de);
      Omega_m = (pba->Omega0_m*pow((1.+z), 3.))/Hubble2;
      /* Samuel brieden: TBC: check that the matching between the
         background quantity and this fitting formula improves by
         using Omega_cb (as it is done in background). Carefull:
         Hubble remains with Omega0_m */

      if (w0 == -1.){
        gamma = 0.55;
      }
      else if (w0 < -1.){
        gamma = 0.55+0.02*(1+w0);
      }
      else {
        gamma = 0.55+0.05*(1+w0);
      }
      integrand[index_scalefactor*index_gcol+index_a] = scalefactor;
      integrand[index_scalefactor*index_gcol+index_growth]= -pow(Omega_m, gamma)/scalefactor;
    }

    class_call(array_spline(integrand,
                            index_gcol,
                            ng,
                            index_a,
                            index_growth,
                            index_ddgrowth,
                            _SPLINE_EST_DERIV_,
                            pfo->error_message),
               pfo->error_message,
               pfo->error_message);

    class_call(array_integrate_all_trapzd_or_spline(integrand,
                                                    index_gcol,
                                                    ng,
                                                    0, //ng-1,
                                                    index_a,
                                                    index_growth,
                                                    index_ddgrowth,
                                                    growth,
                                                    pfo->error_message),
               pfo->error_message,
               pfo->error_message);

    *growth = exp(*growth);

  }
  //fprintf(stdout, "%e %e \n", a, *growth);
  free(pvecback);
  free(integrand);

  return _SUCCESS_;
}

/**
 * This is the fourier transform of the NFW density profile.
 *
 * @param pfo Input: pointer to fourier structure
 * @param k   Input: wave vector
 * @param rv  Input: virial radius
 * @param c   Input: concentration = rv/rs (with scale radius rs)
 * @param window_nfw Output: Window Function of the NFW profile
 * @return the error status
 */

int hmcode_window_nfw(
                      struct fourier *pfo,
                      double k,
                      double rv,
                      double c,
                      double *window_nfw
                      ){

  double si1, si2, ci1, ci2, ks;
  double p1, p2, p3;

  ks = k*rv/c;

  class_call(sine_integral(ks*(1.+c),
                           &si2,
                           pfo->error_message),
             pfo->error_message,
             pfo->error_message);

  class_call(sine_integral(ks,
                           &si1,
                           pfo->error_message),
             pfo->error_message,
             pfo->error_message);

  class_call(cosine_integral(ks*(1.+c),
                             &ci2,
                             pfo->error_message),
             pfo->error_message,
             pfo->error_message);

  class_call(cosine_integral(ks,
                             &ci1,
                             pfo->error_message),
             pfo->error_message,
             pfo->error_message);

  p1=cos(ks)*(ci2-ci1);
  p2=sin(ks)*(si2-si1);
  p3=sin(ks*c)/(ks*(1.+c));

  *window_nfw=p1+p2-p3;
  *window_nfw=*window_nfw/(log(1.+c)-c/(1.+c));

  return _SUCCESS_;
}

/**
 * This is the Sheth-Tormen halo mass function (1999, MNRAS, 308, 119)
 *
 * @param nu   Input: the \f$ \nu \f$ parameter that depends on the halo mass via \f$ \nu(M) = \delta_c/\sigma(M) \f$
 * @param hmf  Output: Value of the halo mass function at this \f$ \nu \f$
 * @return the error status
 */

int hmcode_halomassfunction(
                            double nu,
                            double *hmf
                            ){

  double p, q, A;

  p=0.3;
  q=0.707;
  A=0.21616;

  *hmf=A*(1.+(pow(q*nu*nu, -p)))*exp(-q*nu*nu/2.);

  return _SUCCESS_;
}

/**
 * Compute sigma8(z)
 *
 * @param pba        Input: pointer to background structure
 * @param pfo        Input: pointer to fourier structure
 * @param z          Input: redshift
 * @param sigma_8    Output: sigma8(z)
 * @param sigma_8_cb Output: sigma8_cb(z)
 * @param phw        Output: pointer to hmcode workspace
 * @return the error status
 */

int hmcode_sigma8_at_z(
                       struct background *pba,
                       struct fourier *pfo,
                       double z,
                       double *sigma_8,
                       double *sigma_8_cb,
                       struct hmcode_workspace *phw
                       ) {

  double tau;

  class_call(background_tau_of_z(pba,
                                 z,
                                 &tau),
             pba->error_message,
             pfo->error_message);

  if (pfo->tau_size == 1) {
    *sigma_8 = phw->sigma_8[pfo->index_pk_m][0];
  }
  else {
    class_call(array_interpolate_two(pfo->tau,
                                     1,
                                     0,
                                     phw->sigma_8[pfo->index_pk_m],
                                     1,
                                     pfo->tau_size,
                                     tau,
                                     sigma_8,
                                     1,
                                     pfo->error_message),
               pfo->error_message,
               pfo->error_message);
  }


  if (pba->has_ncdm == _TRUE_){

    if (pfo->tau_size == 1) {
      *sigma_8_cb = phw->sigma_8[pfo->index_pk_cb][0];
    }
    else {
      class_call(array_interpolate_two(pfo->tau,
                                       1,
                                       0,
                                       phw->sigma_8[pfo->index_pk_cb],
                                       1,
                                       pfo->tau_size,
                                       tau,
                                       sigma_8_cb,
                                       1,
                                       pfo->error_message),
                 pfo->error_message,
                 pfo->error_message);
    }

  }
  else{
    *sigma_8_cb = *sigma_8;
  }

  return _SUCCESS_;
}

/**
 * Compute sigmadisp(z)
 *
 * @param pba           Input: pointer to background structure
 * @param pfo           Input: pointer to fourier structure
 * @param z             Input: redshift
 * @param sigma_disp    Output: sigmadisp(z)
 * @param sigma_disp_cb Output: sigmadisp_cb(z)
 * @param phw           Output: pointer to hmcode workspace
 * @return the error status
 */

int hmcode_sigmadisp_at_z(
                          struct background *pba,
                          struct fourier *pfo,
                          double z,
                          double *sigma_disp,
                          double *sigma_disp_cb,
                          struct hmcode_workspace *phw
                          ) {

  double tau;

  class_call(background_tau_of_z(pba,
                                 z,
                                 &tau),
             pba->error_message,
             pfo->error_message);

  if (pfo->tau_size == 1) {
    *sigma_disp = phw->sigma_disp[pfo->index_pk_m][0];
  }
  else {
    class_call(array_interpolate_two(pfo->tau,
                                     1,
                                     0,
                                     phw->sigma_disp[pfo->index_pk_m],
                                     1,
                                     pfo->tau_size,
                                     tau,
                                     sigma_disp,
                                     1,
                                     pfo->error_message),
               pfo->error_message,
               pfo->error_message);
  }

  if (pba->has_ncdm == _TRUE_){

    if (pfo->tau_size == 1) {
      *sigma_disp_cb = phw->sigma_disp[pfo->index_pk_cb][0];
    }
    else {
      class_call(array_interpolate_two(pfo->tau,
                                       1,
                                       0,
                                       phw->sigma_disp[pfo->index_pk_cb],
                                       1,
                                       pfo->tau_size,
                                       tau,
                                       sigma_disp_cb,
                                       1,
                                       pfo->error_message),
                 pfo->error_message,
                 pfo->error_message);
    }

  }
  else{
    *sigma_disp_cb = *sigma_disp;
  }

  return _SUCCESS_;
}

/**
 * Compute sigmadisp100(z)
 *
 * @param pba               Input: pointer to background structure
 * @param pfo               Input: pointer to fourier structure
 * @param z                 Input: redshift
 * @param sigma_disp_100    Output: sigmadisp100(z)
 * @param sigma_disp_100_cb Output: sigmadisp100_cb(z)
 * @param phw               Output: pointer to hmcode workspace
 * @return the error status
 */

int hmcode_sigmadisp100_at_z(
                             struct background *pba,
                             struct fourier *pfo,
                             double z,
                             double *sigma_disp_100,
                             double *sigma_disp_100_cb,
                             struct hmcode_workspace *phw
                             ) {

  double tau;

  class_call(background_tau_of_z(pba,
                                 z,
                                 &tau),
             pba->error_message,
             pfo->error_message);

  if (pfo->tau_size == 1) {
    *sigma_disp_100 = phw->sigma_disp_100[pfo->index_pk_m][0];
  }
  else {
    class_call(array_interpolate_two(pfo->tau,
                                     1,
                                     0,
                                     phw->sigma_disp_100[pfo->index_pk_m],
                                     1,
                                     pfo->tau_size,
                                     tau,
                                     sigma_disp_100,
                                     1,
                                     pfo->error_message),
               pfo->error_message,
               pfo->error_message);
  }

  if (pba->has_ncdm == _TRUE_){

    if (pfo->tau_size == 1) {
      *sigma_disp_100_cb = phw->sigma_disp_100[pfo->index_pk_cb][0];
    }
    else {
      class_call(array_interpolate_two(pfo->tau,
                                       1,
                                       0,
                                       phw->sigma_disp_100[pfo->index_pk_cb],
                                       1,
                                       pfo->tau_size,
                                       tau,
                                       sigma_disp_100_cb,
                                       1,
                                       pfo->error_message),
                 pfo->error_message,
                 pfo->error_message);
    }

  }
  else{
    *sigma_disp_100_cb = *sigma_disp_100;
  }

  return _SUCCESS_;
}

/**
 * Compute sigma'(z)
 *
 * @param pba            Input: pointer to background structure
 * @param pfo            Input: pointer to fourier structure
 * @param z              Input: redshift
 * @param sigma_prime    Output: sigma'(z)
 * @param sigma_prime_cb Output: sigma'_cb(z)
 * @param phw            Output: pointer to hmcode workspace
 * @return the error status
 */

int hmcode_sigmaprime_at_z(
                           struct background *pba,
                           struct fourier *pfo,
                           double z,
                           double *sigma_prime,
                           double *sigma_prime_cb,
                           struct hmcode_workspace *phw
                           ) {

  double tau;

  class_call(background_tau_of_z(pba,
                                 z,
                                 &tau),
             pba->error_message,
             pfo->error_message);

  if (pfo->tau_size == 1) {
    *sigma_prime = phw->sigma_prime[pfo->index_pk_m][0];
  }
  else {
    class_call(array_interpolate_two(pfo->tau,
                                     1,
                                     0,
                                     phw->sigma_prime[pfo->index_pk_m],
                                     1,
                                     pfo->tau_size,
                                     tau,
                                     sigma_prime,
                                     1,
                                     pfo->error_message),
               pfo->error_message,
               pfo->error_message);
  }

  if (pba->has_ncdm == _TRUE_){

    if (pfo->tau_size == 1) {
      *sigma_prime_cb = phw->sigma_prime[pfo->index_pk_cb][0];
    }
    else {
      class_call(array_interpolate_two(pfo->tau,
                                       1,
                                       0,
                                       phw->sigma_prime[pfo->index_pk_cb],
                                       1,
                                       pfo->tau_size,
                                       tau,
                                       sigma_prime_cb,
                                       1,
                                       pfo->error_message),
                 pfo->error_message,
                 pfo->error_message);
    }
  }
  else{
    *sigma_prime_cb = *sigma_prime;
  }

  return _SUCCESS_;
}

/**
 * Compute and store the dimensionless wiggle spectrum {cal P}_wiggle
 *
 * @param pfo            Input: pointer to fourier structure
 * @param lnpk_l         Input: linear power spectrum
 * @param ddlnpk_l       Input: its second derivatives
 * @param index_pk       Input: index of the pk type, either index_m or index_cb
 * @param phw            Output: pointer to hmcode workspace including the wiggle spectrum
 * @return the error status
 */

int hmcode_nowiggle_init(
                         struct fourier *pfo,
                         double **lnpk_l,
                         double **ddlnpk_l,
                         int index_pk,
                         struct hmcode_workspace *phw
                         ) {

  int i;
  double loganorm = log(1/(2.*_PI_*_PI_));
  double * pk_nw;

  /* Temporary array to store P_nowiggle (will not be used here) */

  class_alloc(pk_nw, pfo->nk_wiggle*sizeof(double), pfo->error_message);

  /* Get P_wiggle = P_l - P_nowiggle */

  class_call(hmcode_nowiggle(pfo,
                             lnpk_l[index_pk],
                             ddlnpk_l[index_pk],
                             pfo->nk_wiggle,
                             phw->lnk_wiggle,
                             pk_nw,
                             phw->pk_wiggle),
             pfo->error_message,
             pfo->error_message);

  /* Transform P_wiggle into a dimensionless linear spectrum {cal
     P}_wiggle = k^3/(2pi^2) P_l, where {cal P}_l would scale like
     k^(4+n_s) for k<k_eq. Store the result in phw->pk_wiggle. */

  for(i=0;i<pfo->nk_wiggle;++i){
    phw->pk_wiggle[i] = phw->pk_wiggle[i]*exp(3.*phw->lnk_wiggle[i]+loganorm);
  }

  /* Spline phw->pk_wiggle w.r.t. ln(k), to be able to
     interpolate later on. */

  class_call(array_spline_table_columns(phw->lnk_wiggle,
                                        pfo->nk_wiggle,
                                        phw->pk_wiggle,
                                        1,
                                        phw->ddpk_wiggle,
                                        _SPLINE_EST_DERIV_,
                                        pfo->error_message),
             pfo->error_message,
             pfo->error_message);

  free(pk_nw);

  return _SUCCESS_;
}

/**
 * Compute the decomposition of the linear power spectrum into a
 * wiggly and a non-wiggly part. Store the results in the fourier
 * structure.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param ppm Input: pointer to primordial structure
 * @param pfo Input/Output: pointer to fourier structure
 * @return the error status
 */

int hmcode_wnw_split(
                      struct precision *ppr,
                      struct background *pba,
                      struct primordial * ppm,
                      struct fourier *pfo
                      ) {

  int index_pk, index_k, index_tau;
  int last_index=0;
  double ln_k;

  double * lnpk_l;
  double * ddlnpk_l;

  double * ln_k_nw;
  double * pk_nw;
  double * pk_w;

  double * lnpk_nw;
  double * ddlnpk_nw;


  /** - single index for which we want to compute and store the dewiggled spectrum (refers to _cb or _m) */

  index_pk = *(pfo->pk_l_nw_index);

  /** - fill the nowiggle array with the original linear spectrum.
      Later, we will overwrite with the de-wiggled values */

  memcpy(pfo->ln_pk_l_nw_extra, pfo->ln_pk_l_extra[index_pk],
         pfo->ln_tau_size * pfo->k_size_extra * sizeof(double));

  /** - allocate temporary arrays */

  class_alloc(lnpk_l,pfo->k_size_extra*sizeof(double),pfo->error_message);
  class_alloc(ddlnpk_l,pfo->k_size_extra*sizeof(double),pfo->error_message);

  class_alloc(ln_k_nw, pfo->nk_wiggle*sizeof(double), pfo->error_message);
  class_alloc(pk_nw, pfo->nk_wiggle*sizeof(double), pfo->error_message);
  class_alloc(pk_w, pfo->nk_wiggle*sizeof(double), pfo->error_message);

  class_alloc(lnpk_nw, pfo->nk_wiggle*sizeof(double), pfo->error_message);
  class_alloc(ddlnpk_nw, pfo->nk_wiggle*sizeof(double), pfo->error_message);

  /** - loop over times at which we are storing the power spectra */

  for (index_tau=0; index_tau<pfo->ln_tau_size; index_tau++){

    /** - create a temporary copy of ln(P_l) at this time */
    memcpy(lnpk_l, &(pfo->ln_pk_l_extra[index_pk][index_tau*pfo->k_size_extra]),
           pfo->k_size_extra * sizeof(double));

    /** - spline it with respect to k, needed for interpolation inside hmcode_nowiggle() */

    class_call(array_spline_table_columns(pfo->ln_k,
                                          pfo->k_size_extra,
                                          lnpk_l,
                                          1,
                                          ddlnpk_l,
                                          _SPLINE_NATURAL_,
                                          pfo->error_message),
               pfo->error_message,
               pfo->error_message);

    /* - compute P_nowiggle and P_wiggle (here the later is not used,
         only the former). These spectra are sampled over an array
         ln_k_nw determined by the function hmcode_nowiggle() */

    class_call(hmcode_nowiggle(pfo,
                               lnpk_l,
                               ddlnpk_l,
                               pfo->nk_wiggle,
                               ln_k_nw,
                               pk_nw,
                               pk_w),
               pfo->error_message,
               pfo->error_message);

    /* - compute ln(P_nowiggle) */

    for (index_k=0; index_k<pfo->nk_wiggle; index_k++) {
      lnpk_nw[index_k] = log(pk_nw[index_k]);
    }

    /** - spline ln(P_nowiggle) with respect to ln(k) */

    class_call(array_spline_table_columns(ln_k_nw,
                                          pfo->nk_wiggle,
                                          lnpk_nw,
                                          1,
                                          ddlnpk_nw,
                                          _SPLINE_NATURAL_,
                                          pfo->error_message),
               pfo->error_message,
               pfo->error_message);

    /** - interpolate ln(P_nowiggle) at the ln(k) values
          used in the fourier structure */

    for (index_k=0; index_k<pfo->k_size_extra; index_k++) {

      ln_k = pfo->ln_k[index_k];

      if ((ln_k>ln_k_nw[0]) && (ln_k<ln_k_nw[pfo->nk_wiggle-1])) {

        class_call(array_interpolate_spline(ln_k_nw,
                                            pfo->nk_wiggle,
                                            lnpk_nw,
                                            ddlnpk_nw,
                                            1,
                                            ln_k,
                                            &last_index,
                                            &(pfo->ln_pk_l_nw_extra[index_tau * pfo->k_size_extra + index_k]),
                                            1,
                                            pfo->error_message),
                   pfo->error_message,
                   pfo->error_message);
      }
    }
  }

  /** - free temporary arrays */

  free(lnpk_l);
  free(ddlnpk_l);

  free(ln_k_nw);
  free(pk_nw);
  free(pk_w);

  free(lnpk_nw);
  free(ddlnpk_nw);

  return _SUCCESS_;
}

/**
 * Dewiggle the linear spectrum
 *
 * @param pfo            Input: pointer to fourier structure
 * @param lnpk_l         Input: linear power spectrum
 * @param ddlnpk_l       Input: its second derivatives
 * @param index_pk       Input: index of the pk type, either index_m or index_cb
 * @param nw_size        Input: number of sampled values of the output spectra
 * @param ln_k_nw        Ouput: ln(k) array of size nw_size
 * @param pk_nw          Output: P_nowiggle array of size nw_size
 * @param pk_w           Output: P_wiggle array of size nw_size
 * @return the error status
 */

int hmcode_nowiggle(
                    struct fourier *pfo,
                    double * lnpk_l,
                    double * ddlnpk_l,
                    int nw_size,
                    double * ln_k_nw,
                    double * pk_nw,
                    double * pk_w
                    ) {

  // JL: The following parameters define the smoothing algorithm.
  //     In principle, they should be precision parameters rather than hard-coded.
  //     However, the goal is to match exactly other versions of HMcode, so we leave these fixed.
  double wiggle_sigma = 0.25;
  double kmin_wiggle = 5e-3; // Mpc/h
  double logkmin_wiggle = log(kmin_wiggle);
  double kmax_wiggle = 5.; // Mpc/h
  double logkmax_wiggle = log(kmax_wiggle);

  int i;
  int last_index = 0;
  double * lnpk_lin;
  double * lnpk_analytic_nowiggle;
  double * pk_ratio;

  class_alloc(lnpk_lin, nw_size*sizeof(double), pfo->error_message);
  class_alloc(lnpk_analytic_nowiggle, nw_size*sizeof(double), pfo->error_message);
  class_alloc(pk_ratio, nw_size*sizeof(double), pfo->error_message);

  for(i=0; i<nw_size; ++i){

    /* Define a new sampling in k-space, from logkmin_wiggle to
       logkmax_wiggle with nw_size steps */

    ln_k_nw[i] = logkmin_wiggle+i*(logkmax_wiggle-logkmin_wiggle)/(nw_size-1);

    /* Interpolate the linear spectrum to fill lnpk_lin */

    class_call(array_interpolate_spline(pfo->ln_k,
                                        pfo->k_size_extra,
                                        lnpk_l,
                                        ddlnpk_l,
                                        1,
                                        ln_k_nw[i],
                                        &last_index,
                                        &(lnpk_lin[i]),
                                        1,
                                        pfo->error_message),
               pfo->error_message,
               pfo->error_message);

    /* Interpolate the smooth analytic spectrum to populate
       lnpk_analytic_nowiggle. This spectrum is only used to reduce the
       dynamical range before smoothing. */

    class_call(array_interpolate_spline(pfo->ln_k,
                                        pfo->k_size_extra,
                                        pfo->ln_pk_l_an_extra,
                                        pfo->ddln_pk_l_an_extra,
                                        1,
                                        ln_k_nw[i],
                                        &last_index,
                                        &(lnpk_analytic_nowiggle[i]),
                                        1,
                                        pfo->error_message),
               pfo->error_message,
               pfo->error_message);

    /* Compute the ratio P_l / P_analytic_nowiggle */

    pk_ratio[i] = exp(lnpk_lin[i]-lnpk_analytic_nowiggle[i]);

  }

  /* Perform a Gaussian smoothing of this ratio, the result is stored in pk_nw */

  class_call(array_smooth_Gaussian(ln_k_nw,
                                   pk_ratio,
                                   pk_nw,
                                   nw_size,
                                   wiggle_sigma,
                                   pfo->error_message),
             pfo->error_message,
             pfo->error_message);

  /* Multiply the smoothed ratio by the previously used analytic
     approximation to get the actual nowiggle spectrum P_nowiggle, and
     compute the difference P_wiggle = P_l - P_nowiggle */

  for(i=0; i<nw_size; ++i){

    pk_nw[i] = pk_nw[i]*exp(lnpk_analytic_nowiggle[i]);
    pk_w[i] = exp(lnpk_lin[i])-pk_nw[i];
  }

  free(lnpk_lin);
  free(lnpk_analytic_nowiggle);
  free(pk_ratio);

  return _SUCCESS_;
}

/**
 * No-wiggle spectrum form Eisenstein & Hu
 *
 * @param pfo            Input: pointer to fourier structure
 * @param phw            Input: pointer to hmcode workspace
 * @param pk_norm        Output: nowiggle power spectrum
 * @return the error status
 */

int hmcode_eisenstein_hu(
                         struct precision *ppr,
                         struct background *pba,
                         struct primordial * ppm,
                         struct fourier *pfo,
                         double * ln_k,
                         int k_size,
                         double * ln_pk_nowiggle
                         ){

  int index_k;
  double h,wm,wb,rb;
  double s, alpha, Gamma, q,L,C,Tk_nw,Pk_nw;
  double k;
  double ns,Om,OL,g,coeff,norm;

  h = pba->h;
  wm = (1.-pba->Omega0_lambda-pba->Omega0_fld)*h*h;
  wb = pba->Omega0_b*h*h;
  rb = wb/wm;

  s = 44.5*log(9.83/wm)/sqrt(1.+10.*pow(wb,0.75));
  alpha = 1.-0.328*log(431.*wm)*rb+0.38*log(22.3*wm)*rb*rb;

  // ns was hard-coded as ns=0.95 in old HMcode implementation (this has no
  // significant impact since the analytic power spectrum is used only to
  // reduce the dynamical range before interpolation!)
  //ns = 0.95;
  ns = ppm->n_s;

  // added by JL for an approximate but realistic normalisation of the spectrum:
  // - factor g(Omega_M) defined and approximated like in Kofman et al. 1993 (valid in a flat Universe)
  // - factor norm refering to eq. (6.34, first line) in Neutrino Cosmology book (CUP 2013)
  //   (the original HMcode was using norm=1, which makes zero difference for the calculation of P_NL)
  Om = (1.-pba->Omega0_lambda-pba->Omega0_fld);
  OL = pba->Omega0_lambda+pba->Omega0_fld;
  g = pow(Om,0.2) / (1.+0.003*pow(OL/Om,4./3.));
  coeff = g/Om/pba->H0/pba->H0;
  // if {\cal P}_R = A_s (k/kp)^(ns-1), then [{\cal P}_R / k^(ns-1)] = A_s / kp^(ns-1)
  norm = pow(coeff,2)*8*_PI_*_PI_/25.*ppm->A_s/pow(ppm->k_pivot,ns-1.);

  for(index_k=0; index_k<k_size; index_k++){
    k = exp(ln_k[index_k]);
    Gamma = (wm/h)*(alpha+(1.-alpha)/(1.+pow(0.43*k*s,4)));
    q = k/h*(pba->T_cmb/2.7)*(pba->T_cmb/2.7)/Gamma;
    L = log(2.*exp(1.)+1.8*q);
    C = 14.2+731./(1.+62.5*q);
    Tk_nw = L/(L+C*q*q);
    // According to eq. (6.34, first line) in Neutrino Cosmology book (CUP 2013),
    // P(k) = norm * k^ns * T(k) where T(k) is noramlised to 1 for small k
    ln_pk_nowiggle[index_k] = log(norm*pow(k,ns)*pow(Tk_nw,2));
  }

  return _SUCCESS_;
}

/**
 * Initialise workspace to compute the growth factor. The growth factor is
 * computed internally within HMcode by integrating some simplified
 * perturbation equations. This is redundant with things done by CLASS
 * in the background module, and it could be by-passed at some point.
 * But here we do it exactly like in the original HMcode code, in
 * order to find exactly the same P_NL/P_L ratio.
 *
 * @param ppr            Input: pointer to precision structure
 * @param pba            Input: pointer to background structure
 * @param pfo            Input: pointer to fourier structure
 * @param phw            Output: pointer to hmcode workspace
 * @return the error status
 */

int hmcode_noradiation_growth_init(
                                   struct precision *ppr,
                                   struct background *pba,
                                   struct fourier *pfo,
                                   struct hmcode_workspace *phw
                                   ){

  struct hmcode_growth *phg;
  int index_fg,index_gt, index_norad;
  int i;
  double a_ini;
  double a_final;
  double lna_ini;
  double lna_final;
  double neutrino_constant;
  double Omega_nu_rad;

  class_alloc(phg,sizeof(struct hmcode_growth),pfo->error_message);
  phw->phg = phg;

  index_fg = 0;
  class_define_index(phg->index_fg_g,_TRUE_,index_fg,1);
  class_define_index(phg->index_fg_dg,_TRUE_,index_fg,1);

  phg->fg_size = index_fg;

  index_gt = 0;
  class_define_index(phg->index_gt_g,_TRUE_,index_gt,1); //Actually g(a)/a
  class_define_index(phg->index_gt_ddg,_TRUE_,index_gt,1);
  class_define_index(phg->index_gt_intg,_TRUE_,index_gt,1);
  class_define_index(phg->index_gt_ddintg,_TRUE_,index_gt,1);
  class_define_index(phg->index_gt_Omnorad,_TRUE_,index_gt,1);
  class_define_index(phg->index_gt_ddOmnorad,_TRUE_,index_gt,1);

  phg->gt_size = index_gt;

  index_norad = 0;
  class_define_index(phg->index_norad_H2,_TRUE_,index_norad,1);
  class_define_index(phg->index_norad_AH,_TRUE_,index_norad,1);
  class_define_index(phg->index_norad_Om,_TRUE_,index_norad,1);

  phg->norad_size = index_norad;

  a_ini = 1e-4;
  a_final = 1;
  lna_ini = log(a_ini);
  lna_final = log(a_final);

  phg->a_size = 128;
  class_alloc(phg->a_table,phg->a_size*sizeof(double),pfo->error_message);
  class_alloc(phg->growth_table,phg->a_size*phg->gt_size*sizeof(double),pfo->error_message);
  class_alloc(phg->normgrowth_table,phg->a_size*sizeof(double),pfo->error_message);
  class_alloc(phg->ddnormgrowth_table,phg->a_size*sizeof(double),pfo->error_message);
  class_alloc(phg->dda_table,phg->a_size*sizeof(double),pfo->error_message);
  class_alloc(phg->pvecnorad,phg->norad_size*sizeof(double),pfo->error_message);

  for(i=0;i<phg->a_size;++i){
    phg->a_table[i] = a_ini*exp(i*(lna_final-lna_ini)/(phg->a_size-1));
  }

  phg->a_ini = a_ini;
  phg->a_final = a_final;
  phg->a_table[0]=a_ini;
  phg->a_table[phg->a_size-1]=a_final;
  phg->smallest_allowed_variation = ppr->smallest_allowed_variation;

  phg->Omega_cdm = pba->Omega0_m-pba->Omega0_b;//pba->Omega0_cdm;
  phg->Omega_b = pba->Omega0_b;
  phg->Omega_v = pba->Omega0_lambda;
  Omega_nu_rad = pba->Omega0_g*pba->Neff*7./8.*pow(4./11.,4./3.);

  phg->Omega_nu = Omega_nu_rad;
  phg->a_nu = 1;
  // These would be the expressions for mnu>0, but to emulate the default HMcode, we use the expressions above
  //phg->a_nu = pba->T_cmb*pow(4./11.,1./3.) /(mnu/pba->Neff)*(1.38064852e-23/1.60218e-19);
  //phg->Omega_nu = Omega_nu_rad*(pow(1.+pow(0.3173*1./phg->a_nu,1.83),1./1.83));

  phg->Omega_w = pba->Omega0_fld;
  phg->w0 = pba->w0_fld;
  phg->wa = pba->wa_fld;
  phg->Omega_m = 1.-phg->Omega_v-phg->Omega_w;
  phg->Tcmb = pba->T_cmb;
  phg->om_m = pba->Omega0_m*pba->h*pba->h;
  phg->fnu = pba->Omega0_ncdm_tot/pba->Omega0_m;
  class_test(pba->Omega0_scf>0,
             pfo->error_message,
             "Cannot use scalar field scf with HMcode (not yet coded)");

  class_call(hmcode_noradiation_growth_compute(pfo, phw),
             phw->phg->error_message,
             pfo->error_message);

  return _SUCCESS_;
}

/**
 * Free workspace to compute the growth factor. The growth factor is
 * computed internally within HMcode by integrating some simplified
 * perturbation equations. This is redundant with things done by CLASS
 * in the background module, and it could be by-passed at some point.
 * But here we do it exactly like in the original HMcode code, in
 * order to find exactly the same P_NL/P_L ratio.
 *
 * @param pfo            Input: pointer to fourier structure
 * @param phw            Output: pointer to hmcode workspace
 * @return the error status
 */

int hmcode_noradiation_growth_free(
                                   struct fourier *pfo,
                                   struct hmcode_workspace *phw
                                   ){

  struct hmcode_growth *phg;

  phg = phw->phg;
  free(phg->a_table);
  free(phg->growth_table);
  free(phg->normgrowth_table);
  free(phg->ddnormgrowth_table);
  free(phg->dda_table);
  free(phg->pvecnorad);
  free(phg);

  return _SUCCESS_;
}

/**
 * Compute the growth factor. The growth factor is computed internally
 * within HMcode by integrating some simplified perturbation
 * equations. This is redundant with things done by CLASS in the
 * background module, and it could be by-passed at some point. But
 * here we do it exactly like in the original HMcode code, in order to
 * find exactly the same P_NL/P_L ratio.
 *
 * @param pfo            Input: pointer to fourier structure
 * @param phw            Output: pointer to hmcode workspace
 * @return the error status
 */

int hmcode_noradiation_growth_compute(
                                      struct fourier *pfo,
                                      struct hmcode_workspace *phw
                                      ){

  struct hmcode_growth *phg;
  int i;
  double *pvec_growth;
  int *used_in_output;
  //int (*generic_evolver)() = evolver_ndf15;
  double f;
  /* This could be moved to the precision stucture, but to be sure to
     get the same as the original HMcode, we hard-code it to the same
     value: */
  double tol_hmcode_growth_integration = 1e-10;

  /* evolvers */
  extern int evolver_rk(EVOLVER_PROTOTYPE);
  extern int evolver_ndf15(EVOLVER_PROTOTYPE);
  int (*generic_evolver)(EVOLVER_PROTOTYPE) = evolver_ndf15;

  phg = phw->phg;

  class_alloc(pvec_growth,phg->fg_size*sizeof(double),phg->error_message);
  class_alloc(used_in_output,phg->fg_size*sizeof(int),phg->error_message);

  for(i=0;i<phg->fg_size;++i){
    used_in_output[i] = _TRUE_;
  }

  class_call(hmcode_norad(phg, phg->a_ini),
             phg->error_message,
             phg->error_message);

  f = 1.-phg->pvecnorad[phg->index_norad_Om];

  pvec_growth[phg->index_fg_g] = pow(phg->a_ini,1.-3.*f/5.);
  pvec_growth[phg->index_fg_dg] = (1.-3.*f/5.)*pow(phg->a_ini,-3.*f/5.);

  /** - perform the integration */
  class_call(generic_evolver(hmcode_growth_derivs,
                             phg->a_ini,
                             phg->a_final,
                             pvec_growth,
                             used_in_output,
                             phg->fg_size,
                             phg,
                             tol_hmcode_growth_integration,
                             phg->smallest_allowed_variation,
                             NULL,
                             1.,
                             phg->a_table,
                             phg->a_size,
                             hmcode_growth_sources,
                             NULL,
                             phg->error_message),
             phg->error_message,
             phg->error_message);

  class_call(array_spline_table_line_to_line(
                                             phg->a_table,
                                             phg->a_size,
                                             phg->growth_table,
                                             phg->gt_size,
                                             phg->index_gt_g,
                                             phg->index_gt_ddg,
                                             _SPLINE_EST_DERIV_,
                                             phg->error_message),
             phg->error_message,
             phg->error_message);

  class_call(array_integrate_spline_table_line_to_line(
                                                       phg->a_table,
                                                       phg->a_size,
                                                       phg->growth_table,
                                                       phg->gt_size,
                                                       phg->index_gt_g,
                                                       phg->index_gt_ddg,
                                                       phg->index_gt_intg,
                                                       phg->error_message),
             phg->error_message,
             phg->error_message);

  class_call(array_spline_table_line_to_line(
                                             phg->a_table,
                                             phg->a_size,
                                             phg->growth_table,
                                             phg->gt_size,
                                             phg->index_gt_intg,
                                             phg->index_gt_ddintg,
                                             _SPLINE_EST_DERIV_,
                                             phg->error_message),
             phg->error_message,
             phg->error_message);

  class_call(array_spline_table_line_to_line(
                                             phg->a_table,
                                             phg->a_size,
                                             phg->growth_table,
                                             phg->gt_size,
                                             phg->index_gt_Omnorad,
                                             phg->index_gt_ddOmnorad,
                                             _SPLINE_EST_DERIV_,
                                             phg->error_message),
             phg->error_message,
             phg->error_message);

  for(i=0;i<phg->a_size;++i){
    phg->normgrowth_table[i] /= phg->normgrowth_table[phg->a_size-1];
    //printf("%i - %.10e %.10e \n",i,phg->a_table[i],phg->normgrowth_table[i]);
  }

  class_call(array_spline_table_columns(
                                        phg->a_table,
                                        phg->a_size,
                                        phg->normgrowth_table,
                                        1,
                                        phg->ddnormgrowth_table,
                                        _SPLINE_EST_DERIV_,
                                        phg->error_message),
             phg->error_message,
             phg->error_message);

  class_call(array_spline_table_columns(
                                        phg->normgrowth_table,
                                        phg->a_size,
                                        phg->a_table,
                                        1,
                                        phg->dda_table,
                                        _SPLINE_EST_DERIV_,
                                        phg->error_message),
             phg->error_message,
             phg->error_message);

  for(i=0;i<phg->a_size;++i){
    //printf("%i - %.10e %.10e %.10e \n",i,phg->a_table[i],phg->ddnormgrowth_table[i],phg->dda_table[i]);
  }

  free(pvec_growth);
  free(used_in_output);

  return _SUCCESS_;
}

/**
 * System of equations used to compute the growth factor. The growth
 * factor is computed internally within HMcode by integrating some
 * simplified perturbation equations. This is redundant with things
 * done by CLASS in the background module, and it could be by-passed
 * at some point. But here we do it exactly like in the original
 * HMcode code, in order to find exactly the same P_NL/P_L ratio.
 *
 * @param a                        Input: scale factor
 * @param y                        Input: vector of variables in the differential system
 * @param dy                       Ouput: vector of derivative of these variables
 * @param parameters_and_workspace Input: quantities needed to write the differential system
 * @return the error status
 */

int hmcode_growth_derivs(
                         double a,
                         double *y,
                         double *dy,
                         void *parameters_and_workspace,
                         ErrorMsg error_message
                         ){

  struct hmcode_growth *phg;
  double Om_norad;
  double AH_norad, H2_norad;

  //phg = parameters_and_workspace;
  phg = (struct hmcode_growth *)parameters_and_workspace;

  class_call(hmcode_norad(phg, a),
             phg->error_message,
             phg->error_message);

  H2_norad = phg->pvecnorad[phg->index_norad_H2];
  AH_norad = phg->pvecnorad[phg->index_norad_AH];
  Om_norad = phg->pvecnorad[phg->index_norad_Om];

  dy[phg->index_fg_g] = y[phg->index_fg_dg];
  dy[phg->index_fg_dg] = 1.5*Om_norad*y[phg->index_fg_g]/(a*a) - (2.+AH_norad/H2_norad)*y[phg->index_fg_dg]/a;

  return _SUCCESS_;
}

/**
 * Quantities that should be stored in memory over the integration of
 * the differential system (in this case, just the growth factor). The
 * growth factor is computed internally within HMcode by integrating
 * some simplified perturbation equations. This is redundant with
 * things done by CLASS in the background module, and it could be
 * by-passed at some point. But here we do it exactly like in the
 * original HMcode code, in order to find exactly the same P_NL/P_L
 * ratio.
 *
 * @param a                        Input: scale factor
 * @param y                        Input: vector of variables in the differential system
 * @param dy                       Ouput: vector of derivative of these variables
 * @param index_a                  Input:index of the scale factor in the y vector
 * @param parameters_and_workspace Input: quantities needed to write the differential system
 * @return the error status
 */

int hmcode_growth_sources(
                          double a,
                          double *y,
                          double *dy,
                          int index_a,
                          void *parameters_and_workspace,
                          ErrorMsg error_message
                          ) {

  struct hmcode_growth *phg;

  //phg = parameters_and_workspace;
  phg = (struct hmcode_growth *)parameters_and_workspace;

  phg->a_table[index_a] = a;
  phg->growth_table[index_a*phg->gt_size+phg->index_gt_g] = y[phg->index_fg_g]/a; //It is actually g(a)/a
  phg->normgrowth_table[index_a] = y[phg->index_fg_g];

  class_call(hmcode_norad(phg, a),
             phg->error_message,
             phg->error_message);

  phg->growth_table[index_a*phg->gt_size+phg->index_gt_Omnorad] = phg->pvecnorad[phg->index_norad_Om];

  return _SUCCESS_;
}

/**
 * The system of equations used to compute the growth factor in HMcode
 * needs to know H(a) and Omega_m(a) in an approximate universe with
 * no radiation contribution.
 *
 * @param phw        Input/Output: pointer to hmcode groth structure
 * @param a          Input: scale factor
 * @return the error status
 */

int hmcode_norad(
                 struct hmcode_growth *phg,
                 double a
                 ){

  double AHnr_nu, H2nr_nu;
  double w_v = -1.;
  double w_de = phg->w0+(1.-a)*phg->wa;
  double X_de = pow(a, -3.*(1.+phg->w0+phg->wa))*exp(-3.*phg->wa*(1.-a));
  double Om_norad;
  double AH_norad, H2_norad;
  double a3;

  a3 = a*a*a;

  if( a > phg->a_nu ){
    AHnr_nu = phg->Omega_nu/a3;
    H2nr_nu = phg->Omega_nu/a3;
  }
  else{
    H2nr_nu = 0.;
    AHnr_nu = 0.;
  }

  AH_norad = -(phg->Omega_cdm/a3+phg->Omega_b/a3+AHnr_nu+phg->Omega_v*(1.+3.*w_v) + phg->Omega_w*(1.+3.*w_de)*X_de)/2.;
  H2_norad = (phg->Omega_cdm/a3 + phg->Omega_b/a3 + H2nr_nu + phg->Omega_v + phg->Omega_w*X_de + (1.-phg->Omega_m-phg->Omega_v-phg->Omega_w)/a/a);

  H2_norad+=phg->Omega_nu*(pow(1.+pow(0.3173*a/phg->a_nu,1.83),1./1.83)-1.)/a3/a;//Correction from inconsistent modeling between X_nu and H2_norad in HMcode

  if( a > phg->a_nu ){
    Om_norad = phg->Omega_m/a3/H2_norad;
  }
  else{
    Om_norad = (phg->Omega_cdm + phg->Omega_b)/a3/H2_norad;
  }
  phg->pvecnorad[phg->index_norad_H2] = H2_norad;
  phg->pvecnorad[phg->index_norad_AH] = AH_norad;
  phg->pvecnorad[phg->index_norad_Om] = Om_norad;

  return _SUCCESS_;
}

/**
 * Approximate ratio of growth factors with/without neutrinos
 *
 * @param phw        Input/Output: pointer to hmcode groth structure
 * @param a          Input: scale factor
 * @return the error status
 */

int hmcode_cbnu_ratio(
                      double k,
                      double z,
                      double fnu,
                      double om_m,
                      double Tcmb,
                      double growth,
                      double *cbnu_ratio
                      ){

  double pcb, BigT, zeq, D, q, Nnu, yfs, Dcb, Dcbnu;

  pcb=(5.-sqrt(1+24*(1.-fnu)))/4;

  BigT=Tcmb/2.7;

  zeq=(2.5e4)*om_m*pow(BigT,-4.);
  //IF (EdS_Tcold_growth) THEN
  //    D=(1.+zeq)/(1.+z) ! EdS solution
  //ELSE
  //    D=(1.+zeq)*ungrow(z, cosm) ! General solution
  //END IF
  D=(1.+zeq)*growth;//(1.+zeq)/(1.+z);

  q=k*BigT*BigT/om_m;

  Nnu = 1;
  yfs=17.2*fnu*(1.+0.488*pow(fnu,-7./6.))*pow(Nnu*q/fnu,2);

  Dcb=pow(1.+pow(D/(1.+yfs),0.7),pcb/0.7);
  Dcbnu=pow(pow(1-fnu,0.7/pcb)+pow(D/(1.+yfs),0.7),pcb/0.7);

  *cbnu_ratio =  Dcb/Dcbnu;

  return _SUCCESS_;
}
