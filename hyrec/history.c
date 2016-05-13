/*************************************************************************************************/
/*                 HYREC: Hydrogen and Helium Recombination Code                                 */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                              */
/*                                                                                               */
/*         history.c: compute reionization history                                               */
/*                                                                                               */
/*         Version: January 2011                                                                 */
/*         Revision history:                                                                     */
/*            - written November 2010                                                            */
/*            - January 2011: changed various switches (notably for post-Saha expansions)        */
/*                             so that they remain valid for arbitrary cosmologies               */
/*************************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "hyrec.h"


/*************************************************************************************************
Cosmological parameters Input/Output
*************************************************************************************************/

void rec_get_cosmoparam(FILE *fin, FILE *fout, REC_COSMOPARAMS *param) {

  /* Cosmology */
  if (fout!=NULL && PROMPT==1) fprintf(fout, "Enter CMB temperature today [Kelvin]: ");
  fscanf(fin, "%lg", &(param->T0));
  if (fout!=NULL && PROMPT==1) fprintf(fout, "Enter baryon density, omega_bh2: ");
  fscanf(fin, "%lg", &(param->obh2));
  if (fout!=NULL && PROMPT==1) fprintf(fout, "Enter total matter (CDM+baryons) density, omega_mh2: ");
  fscanf(fin, "%lg", &(param->omh2));
  if (fout!=NULL && PROMPT==1) fprintf(fout, "Enter curvature, omega_kh2: ");
  fscanf(fin, "%lg", &(param->okh2));
  if (fout!=NULL && PROMPT==1) fprintf(fout, "Enter dark energy density, omega_deh2: ");
  fscanf(fin, "%lg", &(param->odeh2));
  if (fout!=NULL && PROMPT==1) fprintf(fout, "Enter dark energy equation of state parameters, w wa: ");
  fscanf(fin, "%lg %lg", &(param->w0), &(param->wa));
  if (fout!=NULL && PROMPT==1) fprintf(fout, "Enter primordial helium mass fraction, Y: ");
  fscanf(fin, "%lg", &(param->Y));
  if (fout!=NULL && PROMPT==1) fprintf(fout, "Enter effective number of neutrino species, N_nu_eff: ");
  fscanf(fin, "%lg", &(param->Nnueff));

  param->nH0 = 11.223846333047*param->obh2*(1.-param->Y);  /* number density of hudrogen today in m-3 */
  param->fHe = param->Y/(1-param->Y)/3.97153;              /* abundance of helium by number */


  /* Redshift range */
  param->zstart = 8000.;
  param->zend = 0.;
  param->dlna = 8.49e-5;
  param->nz = (long) floor(2+log((1.+param->zstart)/(1.+param->zend))/param->dlna);

  if (fout!=NULL && PROMPT==1) fprintf(fout, "\n");
}

/*************************************************************************************
Hubble expansion parameter in sec^-1
*************************************************************************************/

double rec_HubbleConstant(REC_COSMOPARAMS *param, double z) {

   int j;
   double rho, rho_i[5], ainv;
   double ogh2;

   ainv = 1.+z; /* inverse scale factor */

   /* Matter */
   rho_i[0] = param->omh2 *ainv*ainv*ainv;

   /* Curvature */
   rho_i[1] = param->okh2 *ainv*ainv;

   /* Dark energy */
   rho_i[2] = param->odeh2 * pow(ainv, 3*(1+param->w0)) * exp(3*param->wa* (log(ainv)-1.+1./ainv));

   /* Get radiation density.
    * Coefficient based on 1 AU = 1.49597870691e11 m (JPL SSD) */
   ogh2 = 4.48162687719e-7 * param->T0 * param->T0 * param->T0 * param->T0;
   rho_i[3] = ogh2 *ainv*ainv*ainv*ainv;

   /* Neutrino density -- scaled from photons assuming lower temperature */
   rho_i[4] = 0.227107317660239 * rho_i[3] * param->Nnueff;

   /* Total density, including curvature contribution */
   rho = 0.; for(j=0; j<5; j++) rho += rho_i[j];

   /* Conversion to Hubble rate */
   return( 3.2407792896393e-18 * sqrt(rho) );
}

/*****************************************************************************************
Matter temperature -- 1st order steady state, from Hirata 2008
******************************************************************************************/

double rec_Tmss(double xe, double Tr, double H, double fHe, double nH, double energy_rate,REC_COSMOPARAMS *param) {

  int num_lines=10;
  int last_index;
  double chi_heat;
  ErrorMsg error_message;




  //chi_heat = (1.+2.*xe)/3.; // old approximation from Chen and Kamionkowski

  if (xe < 1.) {
    //Coefficient as revised by Galli et al. 2013 (in fact it is an interpolation by Vivian Poulin of columns 1 and 2 in Table V of Galli et al. 2013)
    array_interpolate_spline(param->annihil_coef_xe,
                            param->annihil_coef_num_lines,
                            param->annihil_coef_heat,
                            param->annihil_coef_dd_heat,
                            1,
                            xe,
                            &last_index,
                            &(chi_heat),
                            1,
                            error_message);
    // chi_heat = (1.+2.*xe)/3.; // old approximation from Chen and Kamionkowski
    if (chi_heat > 1.) chi_heat = 1.;
    // fprintf(stdout,"%e   %e     \n", xe,chi_heat);
    // chi_heat = 0.996857*(1.-pow(1.-pow(xe,0.300134),1.51035));   // coefficient as revised by Galli et al. 2013 (in fact it is a fit by Vivian Poulin of columns 1 and 2 in Table V of Galli et al. 2013)
  }
  else
    chi_heat = 1.;

  return Tr/(1.+H/4.91466895548409e-22/Tr/Tr/Tr/Tr*(1.+xe+fHe)/xe)
    +2./3./kBoltz*chi_heat/nH*energy_rate/(4.91466895548409e-22*pow(Tr,4)*xe);

  /* Coefficient = 8 sigma_T a_r / (3 m_e c) */
}

/******************************************************************************************
Matter temperature evolution derivative
******************************************************************************************/

double rec_dTmdlna(double xe, double Tm, double Tr, double H, double fHe, double nH, double energy_rate,REC_COSMOPARAMS *param) {
  int last_index;
  double chi_heat;
  ErrorMsg error_message;

  if (xe < 1.) {
    //Coefficient as revised by Galli et al. 2013 (in fact it is an interpolation by Vivian Poulin of columns 1 and 2 in Table V of Galli et al. 2013)
    array_interpolate_spline(param->annihil_coef_xe,
                                        param->annihil_coef_num_lines,
                                        param->annihil_coef_heat,
                                        param->annihil_coef_dd_heat,
                                        1,
                                        xe,
                                        &last_index,
                                        &(chi_heat),
                                        1,
                                        error_message);
    // chi_heat = (1.+2.*xe)/3.; // old approximation from Chen and Kamionkowski
    // chi_heat = 0.996857*(1.-pow(1.-pow(xe,0.300134),1.51035));// coefficient as revised by Galli et al. 2013 (in fact it is a fit by Vivian Poulin of columns 1 and 2 in Table V of Galli et al. 2013)
    if (chi_heat > 1.) chi_heat = 1.;
    // fprintf(stdout,"%e   %e     \n", xe,chi_heat);
  }
  else
    chi_heat = 1.;

  return -2.*Tm + 4.91466895548409e-22*Tr*Tr*Tr*Tr*xe/(1.+xe+fHe)*(Tr-Tm)/H
    +2./3./kBoltz*chi_heat/nH*energy_rate/(1.+xe+fHe)/H;
}

/**********************************************************************************************
Second order integrator using derivative from previous time steps
Evolves xe only, assumes Tm is given by the steady-state solution
***********************************************************************************************/

void rec_get_xe_next1(REC_COSMOPARAMS *param, double z1, double xe_in, double *xe_out,
                      HRATEEFF *rate_table, int func_select, unsigned iz, TWO_PHOTON_PARAMS *twog_params,
		      double **logfminus_hist, double *logfminus_Ly_hist[],
                      double *z_prev, double *dxedlna_prev, double *z_prev2, double *dxedlna_prev2) {

  double dxedlna, Tr, nH, ainv, H, Tm;

    Tr = param->T0 * (ainv=1.+z1);
    nH = param->nH0 * ainv*ainv*ainv;
    H = rec_HubbleConstant(param, z1);
    Tm = rec_Tmss(xe_in, Tr, H, param->fHe, nH*1e-6, energy_injection_rate(param,z1), param);

    #if (MODEL == PEEBLES)
        dxedlna = func_select==FUNC_HEI  ? rec_helium_dxedt(xe_in, param->nH0, param->T0, param->fHe, H, z1)/H:
                                           rec_HPeebles_dxedlna(xe_in, nH*1e-6, H, Tm*kBoltz, Tr*kBoltz, energy_injection_rate(param,z1), param);
    #elif (MODEL == RECFAST)
        dxedlna = func_select==FUNC_HEI  ? rec_helium_dxedt(xe_in, param->nH0, param->T0, param->fHe, H, z1)/H:
                                           rec_HRecFast_dxedlna(xe_in, nH*1e-6, H, Tm*kBoltz, Tr*kBoltz, energy_injection_rate(param,z1), param);
    #elif (MODEL == EMLA2s2p)
        dxedlna = func_select==FUNC_HEI  ? rec_helium_dxedt(xe_in, param->nH0, param->T0, param->fHe, H, z1)/H:
                                           rec_HMLA_dxedlna(xe_in, nH*1e-6, H, Tm*kBoltz, Tr*kBoltz, energy_injection_rate(param,z1), param, rate_table);
    #else
        dxedlna = func_select==FUNC_HEI  ? rec_helium_dxedt(xe_in, param->nH0, param->T0, param->fHe, H, z1)/H:
	  func_select==FUNC_H2G  ? rec_HMLA_2photon_dxedlna(xe_in, nH*1e-6, H, Tm*kBoltz, Tr*kBoltz, rate_table, twog_params,
							    param->zstart, param->dlna, logfminus_hist, logfminus_Ly_hist, iz, z1, energy_injection_rate(param,z1), param)
	  :rec_HMLA_dxedlna(xe_in, nH*1e-6, H, Tm*kBoltz, Tr*kBoltz, energy_injection_rate(param,z1), param, rate_table);
    #endif

    *xe_out = xe_in + param->dlna * (1.25 * dxedlna - 0.25 * (*dxedlna_prev2));

    *z_prev2       = *z_prev;
    *dxedlna_prev2 = *dxedlna_prev;
    *z_prev        = z1;
    *dxedlna_prev  = dxedlna;

}

/**********************************************************************************************
Second order integrator using derivative from previous time steps
Evolves xe and Tm simultaneously
***********************************************************************************************/

void rec_get_xe_next2(REC_COSMOPARAMS *param, double z1, double xe_in, double Tm_in, double *xe_out, double *Tm_out,
                      HRATEEFF *rate_table, int func_select, unsigned iz, TWO_PHOTON_PARAMS *twog_params,
		      double **logfminus_hist, double *logfminus_Ly_hist[],
                      double *z_prev, double *dxedlna_prev, double *dTmdlna_prev,
                      double *z_prev2, double *dxedlna_prev2, double *dTmdlna_prev2) {

    double dxedlna, dTmdlna, Tr, nH, ainv, H;

    Tr = param->T0 * (ainv=1.+z1);
    nH = param->nH0 * ainv*ainv*ainv;
    double f_esc=0.2;
    double Zeta_ion=pow(10,53.14);
    // double rho_sfr = 0.01376*pow(ainv,3.26)/(1+pow((ainv)/2.59,5.68))*ainv*ainv*ainv;//Comoving to physical
    double rho_sfr = 0.01376*pow(ainv,3.26)/(1+pow((ainv)/2.59,5.68))*ainv*ainv*ainv*exp(-z1/12);//Comoving to physical
    // rho_sfr =0;
    double dNion_over_dt;
    double erg_to_ev = 6.24150913*pow(10,11);
    double E_x = 3.4*pow(10,40)*erg_to_ev;
    // double rho_sfr = (0.01/pow(8,-3.6))*pow(1+z1,-3.6);
    double stars_xe;
    double MPCcube_to_mcube=2.93799895*pow(10,67);
    double MPCcube_to_cmcube=2.93799895*pow(10,73);
    dNion_over_dt = f_esc*Zeta_ion*rho_sfr;
    double f_abs =1;
    double f_X =0.2;
    H = rec_HubbleConstant(param, z1);
    double L_x = E_x  * f_X* rho_sfr/(3*MPCcube_to_mcube*kBoltz*nH*H*(1.+xe_in+param->fHe));
    // fprintf(stdout,"z= %e, L_x = %e ",z1,L_x);

    // if(z1>25)stars_xe = 0;
    // else
    stars_xe = dNion_over_dt/(MPCcube_to_mcube)/(H*nH);
    // stars_xe = 0;

    #if (MODEL == PEEBLES)
         dxedlna = func_select==FUNC_HEI  ? rec_helium_dxedt(xe_in, param->nH0, param->T0, param->fHe, H, z1)/H:
                                            rec_HPeebles_dxedlna(xe_in, nH*1e-6, H, Tm_in*kBoltz, Tr*kBoltz, energy_injection_rate(param,z1), param);
    #elif (MODEL == RECFAST)
         dxedlna = func_select==FUNC_HEI  ? rec_helium_dxedt(xe_in, param->nH0, param->T0, param->fHe, H, z1)/H:
                                            rec_HRecFast_dxedlna(xe_in, nH*1e-6, H, Tm_in*kBoltz, Tr*kBoltz, energy_injection_rate(param,z1), param);
    #elif (MODEL == EMLA2s2p)
         dxedlna = func_select==FUNC_HEI  ? rec_helium_dxedt(xe_in, param->nH0, param->T0, param->fHe, H, z1)/H:
                                            rec_HMLA_dxedlna(xe_in, nH*1e-6, H, Tm_in*kBoltz, Tr*kBoltz, energy_injection_rate(param,z1), param, rate_table);
    #else
         dxedlna = func_select==FUNC_HEI  ? rec_helium_dxedt(xe_in, param->nH0, param->T0, param->fHe, H, z1)/H:
                   func_select==FUNC_H2G  ? rec_HMLA_2photon_dxedlna(xe_in, nH*1e-6, H, Tm_in*kBoltz, Tr*kBoltz, rate_table, twog_params,
                                                                     param->zstart, param->dlna, logfminus_hist, logfminus_Ly_hist, iz, z1,
																	 energy_injection_rate(param,z1), param):
	           func_select==FUNC_HMLA ? rec_HMLA_dxedlna(xe_in, nH*1e-6, H, Tm_in*kBoltz, Tr*kBoltz, energy_injection_rate(param,z1), param, rate_table)
                                           :rec_HPeebles_dxedlna(xe_in, nH*1e-6, H, Tm_in*kBoltz, Tr*kBoltz, energy_injection_rate(param,z1), param); /* used for z < 20 only */
    #endif

	 dTmdlna = rec_dTmdlna(xe_in, Tm_in, Tr, H, param->fHe, nH*1e-6, energy_injection_rate(param,z1), param);
  //

  //
  if(param->reio_parametrization==1){
    dTmdlna+=L_x*(1+2*xe_in)/3.;
    // dTmdlna+=L_x*(1+2*xe_in)/3.*10;
    dxedlna+=stars_xe*((1-xe_in)/3)*2*3;
    // fprintf(stdout, "Computing star reionisation, rho_sfr = %e,z1 = %e\n",rho_sfr,z1 );
    /*******************Helium**********************/
    dxedlna+=stars_xe*param->fHe*(1+tanh((6-z1)/0.5));
    if(z1<6)dxedlna+=stars_xe*param->fHe*(1+tanh((3.5-z1)/0.5));
    /***********************************************/
  }


  //  fprintf(stdout, "%e\n",stars_xe);
    *xe_out = xe_in + param->dlna * (1.25 * (dxedlna) - 0.25 * (*dxedlna_prev2));
    /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
    *Tm_out = Tm_in + param->dlna * (1.25 * dTmdlna - 0.25 * (*dTmdlna_prev2));
    //  fprintf(stdout, "%e\n",*Tm_out);

    *z_prev2       = *z_prev;
    *dxedlna_prev2 = *dxedlna_prev;
    *dTmdlna_prev2 = *dTmdlna_prev;
    *z_prev        = z1;
    *dxedlna_prev  = dxedlna;
    *dTmdlna_prev  = dTmdlna;


}

/****************************************************************************************************
Builds a recombination history
****************************************************************************************************/

void rec_build_history(REC_COSMOPARAMS *param, HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog_params,
                       double *xe_output, double *Tm_output) {


   long iz;
   double **logfminus_hist;
   double *logfminus_Ly_hist[3];
   double H, z, z_prev, dxedlna_prev, z_prev2, dxedlna_prev2, dTmdlna_prev, dTmdlna_prev2;
   double Delta_xe;

   /* history of photon occupation numbers */
   logfminus_hist = create_2D_array(NVIRT, param->nz);
   logfminus_Ly_hist[0] = create_1D_array(param->nz);   /* Ly-alpha */
   logfminus_Ly_hist[1] = create_1D_array(param->nz);   /* Ly-beta  */
   logfminus_Ly_hist[2] = create_1D_array(param->nz);   /* Ly-gamma */


   z = param->zstart;


   /********* He II + III Saha phase *********/
   Delta_xe = 1.;   /* Delta_xe = xHeIII */

   for(iz=0; iz<param->nz && Delta_xe > 1e-9; iz++) {
      z = (1.+param->zstart)*exp(-param->dlna*iz) - 1.;
      xe_output[iz] = rec_sahaHeII(param->nH0,param->T0,param->fHe,z, &Delta_xe);
      Tm_output[iz] = param->T0 * (1.+z);
   }

   /******* He I + II post-Saha phase *********/
   Delta_xe = 0.;     /* Delta_xe = xe - xe(Saha) */

   for (; iz<param->nz && Delta_xe < 5e-4; iz++) {
      z = (1.+param->zstart)*exp(-param->dlna*iz) - 1.;
      xe_output[iz] = xe_PostSahaHe(param->nH0,param->T0,param->fHe, rec_HubbleConstant(param,z), z, &Delta_xe);
      Tm_output[iz] = param->T0 * (1.+z);
   }

   /****** Segment where we follow the helium recombination evolution, Tm fixed to steady state *******/

   z_prev2 = (1.+param->zstart)*exp(-param->dlna*(iz-3)) - 1.;
   dxedlna_prev2 = (xe_output[iz-2] - xe_output[iz-4])/2./param->dlna;

   z_prev = (1.+param->zstart)*exp(-param->dlna*(iz-2)) - 1.;
   dxedlna_prev = (xe_output[iz-1] - xe_output[iz-3])/2./param->dlna;

   Delta_xe = 1.;  /* Difference between xe and H-Saha value */

   for(; iz<param->nz && (Delta_xe > 1e-4 || z > 1650.); iz++) {

      rec_get_xe_next1(param, z, xe_output[iz-1], xe_output+iz, rate_table, FUNC_HEI, iz-1, twog_params,
		     logfminus_hist, logfminus_Ly_hist, &z_prev, &dxedlna_prev, &z_prev2, &dxedlna_prev2);

      z = (1.+param->zstart)*exp(-param->dlna*iz) - 1.;
      Tm_output[iz] = rec_Tmss(xe_output[iz], param->T0*(1.+z), rec_HubbleConstant(param, z), param->fHe, param->nH0*cube(1.+z), energy_injection_rate(param,z),param);

      /* Starting to populate the photon occupation number with thermal values */
      update_fminus_Saha(logfminus_hist, logfminus_Ly_hist, xe_output[iz], param->T0*(1.+z)*kBoltz,
                         param->nH0*cube(1.+z)*1e-6, twog_params, param->zstart, param->dlna, iz, z, 0);

      Delta_xe = fabs(xe_output[iz]- rec_saha_xe_H(param->nH0, param->T0, z));
    }


   /******* Hydrogen post-Saha equilibrium phase *********/
   Delta_xe = 0.;  /*Difference between xe and Saha value */

   for(; iz<param->nz && Delta_xe < 5e-5; iz++) {
      z = (1.+param->zstart)*exp(-param->dlna*iz) - 1.;
      H = rec_HubbleConstant(param,z);
      xe_output[iz] =  xe_PostSahaH(param->nH0*cube(1.+z)*1e-6, H, kBoltz*param->T0*(1.+z), rate_table, twog_params,
				    param->zstart, param->dlna, logfminus_hist, logfminus_Ly_hist, iz, z, &Delta_xe, MODEL, energy_injection_rate(param,z),param);
      Tm_output[iz] = rec_Tmss(xe_output[iz], param->T0*(1.+z), H, param->fHe, param->nH0*cube(1.+z), energy_injection_rate(param,z),param);
    }

    /******* Segment where we follow the hydrogen recombination evolution with two-photon processes
             Tm fixed to steady state ******/

    z_prev2 = (1.+param->zstart)*exp(-param->dlna*(iz-3)) - 1.;
    dxedlna_prev2 = (xe_output[iz-2] - xe_output[iz-4])/2./param->dlna;

    z_prev = (1.+param->zstart)*exp(-param->dlna*(iz-2)) - 1.;
    dxedlna_prev = (xe_output[iz-1] - xe_output[iz-3])/2./param->dlna;

    for(; iz<param->nz && 1.-Tm_output[iz-1]/param->T0/(1.+z) < 5e-4 && z > 700.; iz++) {

       rec_get_xe_next1(param, z, xe_output[iz-1], xe_output+iz, rate_table, FUNC_H2G, iz-1, twog_params,
               	      logfminus_hist, logfminus_Ly_hist, &z_prev, &dxedlna_prev, &z_prev2, &dxedlna_prev2);
       z = (1.+param->zstart)*exp(-param->dlna*iz) - 1.;
       Tm_output[iz] = rec_Tmss(xe_output[iz], param->T0*(1.+z), rec_HubbleConstant(param, z), param->fHe, param->nH0*cube(1.+z)*1e-6, energy_injection_rate(param,z),param);
    }

   /******* Segment where we follow the hydrogen recombination evolution with two-photon processes
            AND Tm evolution ******/

    dTmdlna_prev2 = rec_dTmdlna(xe_output[iz-3], Tm_output[iz-3], param->T0*(1.+z_prev2),
                                rec_HubbleConstant(param, z_prev2), param->fHe, param->nH0*cube(1.+z_prev2), energy_injection_rate(param,z_prev2),param);
    dTmdlna_prev  = rec_dTmdlna(xe_output[iz-2], Tm_output[iz-2], param->T0*(1.+z_prev),
                                rec_HubbleConstant(param, z_prev), param->fHe, param->nH0*cube(1.+z_prev), energy_injection_rate(param,z_prev),param);

    for(; iz<param->nz && z > 700.; iz++) {

        rec_get_xe_next2(param, z, xe_output[iz-1], Tm_output[iz-1], xe_output+iz, Tm_output+iz, rate_table, FUNC_H2G,
                        iz-1, twog_params, logfminus_hist, logfminus_Ly_hist, &z_prev, &dxedlna_prev, &dTmdlna_prev,
                        &z_prev2, &dxedlna_prev2, &dTmdlna_prev2);
        z = (1.+param->zstart)*exp(-param->dlna*iz) - 1.;
     }

    /***** Segment where we follow Tm as well as xe *****/
    /* Radiative transfer effects switched off here */
    for(; iz<param->nz && z>20.; iz++) {

        rec_get_xe_next2(param, z, xe_output[iz-1], Tm_output[iz-1], xe_output+iz, Tm_output+iz, rate_table, FUNC_HMLA,
                        iz-1, twog_params, logfminus_hist, logfminus_Ly_hist, &z_prev, &dxedlna_prev, &dTmdlna_prev,
                        &z_prev2, &dxedlna_prev2, &dTmdlna_prev2);
        z = (1.+param->zstart)*exp(-param->dlna*iz) - 1.;
    }

    /*** For z < 20 use Peeble's model. The precise model does not metter much here as
            1) the free electron fraction is basically zero (~1e-4) in any case and
            2) the universe is going to be reionized around that epoch
         Tm is still evolved explicitly ***/
    for(; iz<param->nz; iz++) {

        rec_get_xe_next2(param, z, xe_output[iz-1], Tm_output[iz-1], xe_output+iz, Tm_output+iz, rate_table, FUNC_PEEBLES,
                        iz-1, twog_params, logfminus_hist, logfminus_Ly_hist, &z_prev, &dxedlna_prev, &dTmdlna_prev,
                        &z_prev2, &dxedlna_prev2, &dTmdlna_prev2);
        z = (1.+param->zstart)*exp(-param->dlna*iz) - 1.;
    }

    /* Cleanup */
    free_2D_array(logfminus_hist, NVIRT);
    free(logfminus_Ly_hist[0]);
    free(logfminus_Ly_hist[1]);
    free(logfminus_Ly_hist[2]);

}
/* Old version, corrected by Vivian Poulin. Kept for comparaison */
// double onthespot_injection_rate(REC_COSMOPARAMS *param,
// 			     double z) {
//
//   double annihilation_at_z;
//   double rho_cdm_today;
//   double u_min;
//   double erfc;
//
//   /*redshift-dependent annihilation parameter*/
//
//   if (z>param->annihilation_zmax) {
//
//     annihilation_at_z = param->annihilation*
//       exp(-param->annihilation_variation*pow(log((param->annihilation_z+1.)/(param->annihilation_zmax+1.)),2));
//   }
//   else if (z>param->annihilation_zmin) {
//
//     annihilation_at_z = param->annihilation*
//       exp(param->annihilation_variation*(-pow(log((param->annihilation_z+1.)/(param->annihilation_zmax+1.)),2)
// 					 +pow(log((z+1.)/(param->annihilation_zmax+1.)),2)));
//   }
//   else {
//
//     annihilation_at_z = param->annihilation*
//       exp(param->annihilation_variation*(-pow(log((param->annihilation_z+1.)/(param->annihilation_zmax+1.)),2)
// 					 +pow(log((param->annihilation_zmin+1.)/(param->annihilation_zmax+1.)),2)));
//   }
//
//   rho_cdm_today = param->omh2*1.44729366e-9; /* energy density in Kg/m^3 */
//
//   u_min = (1+z)/(1+param->annihilation_z_halo);
//
//   erfc = pow(1.+0.278393*u_min+0.230389*u_min*u_min+0.000972*u_min*u_min*u_min+0.078108*u_min*u_min*u_min*u_min,-4);
//   // fprintf(stdout,"param->annihilation = %e\n",param->annihilation);
//
//   return (pow(rho_cdm_today,2)/2.99792458e8/2.99792458e8*pow((1.+z),3)*
//     (pow((1.+z),3)*annihilation_at_z+param->annihilation_f_halo*erfc)
//     +rho_cdm_today*pow((1+z),3)*param->decay)/1.e6/1.60217653e-19;
//
//   /* energy density rate in eV/cm^3/s (remember that annihilation_at_z is in m^3/s/Kg and decay in s^-1) */
//   /* note that the injection rate used by recfast, defined in therodynamics.c, is in J/m^3/s. Here we multiplied by 1/1.e6/1.60217653e-19 to convert to eV and cm. */
//
// }
/**********************************************************************************************/

/*************************New version, corrected by Vivian Poulin******************************/
double onthespot_injection_rate(REC_COSMOPARAMS *param,
			     double z) {

  double annihilation_at_z;
  double rho_ini_dcdm;
  double result_integrale;
  double rho_cdm_today;
  double _Mpc_over_m_;
  double energy_rate;
  rho_ini_dcdm = param->odcdmh2*1.44729366e-9; /* energy density in Kg/m^3 */
  rho_cdm_today = param->ocdmh2*1.44729366e-9; /* energy density in Kg/m^3 */
  _Mpc_over_m_ = 3.085677581282*pow(10,22);
    if(param->decay>0){

       result_integrale = exp(-param->Gamma_dcdm*2*((param->Omega0_b+param->Omega0_cdm)*pow(param->Omega0_g+(param->Omega0_b+param->Omega0_cdm)/(1+z),0.5)
      +2*pow(param->Omega0_g,1.5)*(1+z)-2*param->Omega0_g*pow((1+z)*(param->Omega0_g*(1+z)+(param->Omega0_b+param->Omega0_cdm)),0.5))/(3*pow((param->Omega0_b+param->Omega0_cdm),2)*(1+z)*param->H0));

      if(rho_ini_dcdm!=0)energy_rate=rho_ini_dcdm*pow((1+z),3)*param->decay*result_integrale*(param->Gamma_dcdm*2.99792458e8/_Mpc_over_m_)/1.e6/1.60217653e-19;
      else energy_rate=rho_cdm_today*pow((1+z),3)*param->decay*result_integrale*(param->Gamma_dcdm*2.99792458e8/_Mpc_over_m_)/1.e6/1.60217653e-19;
      // if(z>0)fprintf(stdout, "z = %e energy_rate = %e\n", z, energy_rate*1.e6*1.60217653e-19);
    }

if(param->annihilation>0.)energy_rate =  (pow(rho_cdm_today,2)/2.99792458e8/2.99792458e8*pow((1.+z),6)*param->annihilation)/1.e6/1.60217653e-19;

    /* Old version, kept for comparaison */
  // return (pow(rho_cdm_today,2)/2.99792458e8/2.99792458e8*pow((1.+z),6)*
  //   param->annihilation
  //   +rho_cdm_today*pow((1+z),3)*param->decay)/1.e6/1.60217653e-19;
  // fprintf(stdout,"z = %e, energy_rate = %e\n",z,energy_rate);

return energy_rate;
  /* energy density rate in eV/cm^3/s (remember that annihilation_at_z is in m^3/s/Kg and decay in s^-1) */
  /* note that the injection rate used by recfast, defined in therodynamics.c, is in J/m^3/s. Here we multiplied by 1/1.e6/1.60217653e-19 to convert to eV and cm. */

}
double beyond_onthespot_injection_rate( REC_COSMOPARAMS *param,
                                     double z) {

  double rho_cdm_today,rho_ini_dcdm,_Mpc_over_m_;
  _Mpc_over_m_ = 3.085677581282*pow(10,22);

  double sigma_thermal = 3*pow(10,-32); // Sigma_v in m^3/s
  double conversion = 1.8*pow(10,-27); // Conversion GeV => Kg
  double f_halos;
  int last_index;
  ErrorMsg error_message;
  double energy_rate;
  double zp,dz;
  double integrand,first_integrand;
  double factor;
  double Boost_factor;
  /*redshift-dependent annihilation parameter*/


  rho_ini_dcdm = param->odcdmh2*1.44729366e-9; /* energy density in Kg/m^3 */
  rho_cdm_today = param->omh2*1.44729366e-9; /* energy density in Kg/m^3 */
  if(rho_ini_dcdm==0)rho_ini_dcdm = rho_cdm_today;
  if (param->annihilation_f_halo > 0. && param->annihilation_z_halo > 0.){
  array_interpolate_spline(param->annihil_z,
                                      param->annihil_f_halos_num_lines,
                                      param->annihil_f_halos,
                                      param->annihil_dd_f_halos,
                                      1,
                                      z,
                                      &last_index,
                                      &(f_halos),
                                      1,
                                      error_message);
  Boost_factor = param->annihilation_f_halo*erfc((1+z)/(1+param->annihilation_z_halo))/pow(1+z,3);
  energy_rate = pow(rho_cdm_today,2)/2.99792458e8/2.99792458e8*pow((1+z),6)*(1+Boost_factor)*param->annihilation*f_halos/1.e6/1.60217653e-19;
  // fprintf(stdout, "%e %e\n",z,Boost_factor);
//
  /* energy density rate in eV/cm^3/s (remember that sigma_thermal/(preco->annihilation_m_DM*conversion) is in m^3/s/Kg) */
  // fprintf(stdout,"%e   %e   %e    %e\n", z,f_halos,energy_rate,1+Boost_factor);
  }
  else{
    array_interpolate_spline(param->annihil_z,
                                        param->annihil_f_halos_num_lines,
                                        param->annihil_f_halos,
                                        param->annihil_dd_f_halos,
                                        1,
                                        z,
                                        &last_index,
                                        &(f_halos),
                                        1,
                                        error_message);
    // if(f_halos<0)f_halos*=-1;
    // fprintf(stdout, "f_halos %e\n", f_halos );

    if(param->annihilation>0)energy_rate = pow(rho_cdm_today,2)/2.99792458e8/2.99792458e8*pow((1+z),6)*param->annihilation*f_halos/1.e6/1.60217653e-19;
    else if(param->decay>0) {
      energy_rate = (rho_ini_dcdm*pow((1+z),3)*param->decay*f_halos*
                                          (param->Gamma_dcdm*2.99792458e8/_Mpc_over_m_))/1.e6/1.60217653e-19;
    // fprintf(stdout,"fhalos = %e, z = %e, energy_rate = %e\n",f_halos,z,energy_rate);
  }

}
    /* energy density rate in eV/cm^3/s (remember that sigma_thermal/(preco->annihilation_m_DM*convers



  // else{
  //
  //       /* factor = c sigma_T n_H(0) / H(0) (dimensionless) */
  //       factor = 2.99792458e8 * 6.6524616e-29 * param->nH0 / (3.2407792896393e-18 * sqrt(param->omh2));
  //
  //       /* integral over z'(=zp) with step dz */
  //       dz=1.;
  //
  //       /* first point in trapezoidal integral */
  //       zp = z;
  //       first_integrand = factor*pow(1+z,8)/pow(1+zp,7.5)*exp(2./3.*factor*(pow(1+z,1.5)-pow(1+zp,1.5)))*onthespot_injection_rate(param,zp); // beware: versions before 2.4.3, there were rwrong exponents: 6 and 5.5 instead of 8 and 7.5
  //       energy_rate = 0.5*dz*first_integrand;
  //
  //       /* other points in trapezoidal integral */
  //       do {
  //
  // 	zp += dz;
  // 	integrand = factor*pow(1+z,8)/pow(1+zp,7.5)*exp(2./3.*factor*(pow(1+z,1.5)-pow(1+zp,1.5)))*onthespot_injection_rate(param,zp); // beware: versions before 2.4.3, there were rwrong exponents: 6 and 5.5 instead of 8 and 7.5
  // 	energy_rate += dz*integrand;
  // 	//moment += dz*integrand*(zp-z);
  //
  //       } while (integrand/first_integrand > 0.02);
  //
  // }
  return energy_rate;

}
double energy_injection_rate(REC_COSMOPARAMS *param,
					double z) {

double result;
double energy_rate;
double zp,dz;
double integrand,first_integrand;
double factor;
  if (param->annihilation > 0. || param->decay > 0.) {

    if (param->has_on_the_spot == 0) {
          //
          //     /* factor = c sigma_T n_H(0) / H(0) (dimensionless) */
          //     factor = 2.99792458e8 * 6.6524616e-29 * param->nH0 / (3.2407792896393e-18 * sqrt(param->omh2));
          //
          //     /* integral over z'(=zp) with step dz */
          //     dz=1.;
          //
          //     /* first point in trapezoidal integral */
          //     zp = z;
          //     first_integrand = factor*pow(1+z,8)/pow(1+zp,7.5)*exp(2./3.*factor*(pow(1+z,1.5)-pow(1+zp,1.5)))*onthespot_injection_rate(param,zp); // beware: versions before 2.4.3, there were rwrong exponents: 6 and 5.5 instead of 8 and 7.5
          //     result = 0.5*dz*first_integrand;
          //
          //     /* other points in trapezoidal integral */
          //     do {
          //
        	// zp += dz;
        	// integrand = factor*pow(1+z,8)/pow(1+zp,7.5)*exp(2./3.*factor*(pow(1+z,1.5)-pow(1+zp,1.5)))*onthespot_injection_rate(param,zp); // beware: versions before 2.4.3, there were rwrong exponents: 6 and 5.5 instead of 8 and 7.5
        	// result += dz*integrand;
        	// //moment += dz*integrand*(zp-z);
          //
          //     } while (integrand/first_integrand > 0.02);
      result = beyond_onthespot_injection_rate(param,z);
      /* test lines for printing energy rate rescaled by (1=z)^6 in J/m^3/s w/o approximation */
      /*  fprintf(stdout,"%e  %e  %e\n",
      1.+z,result/pow(1.+z,6)*1.602176487e-19*1.e6,
      onthespot_injection_rate(param,z)/pow(1.+z,6)*1.602176487e-19*1.e6);
      */
    }
    else {
      result = onthespot_injection_rate(param,z);
    }

    /* effective energy density rate in eV/cm^3/s  */
    return result;
  }
  else {
    return 0.;
  }

}
