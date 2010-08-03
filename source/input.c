/** @file input.c 
 * Julien Lesgourgues, 18.04.2010    
 */

#include "input.h" 

/* Use this routine to extract initial parameters 
   from files 'xxx.ini' and/or 'xxx.pre'. They can be the arguments of the main() routine.

   If class is embedded into another code, you will probably prefer to call directly 
   input_init() in order to pass input parameters through a 'file_content' structure.
 */

int input_init_from_arguments(
	       int argc, 
	       char **argv,
	       struct precision * ppr,
	       struct background *pba,
	       struct thermo *pth,
	       struct perturbs *ppt,
	       struct bessels * pbs,
	       struct transfers *ptr,
	       struct primordial *ppm,
	       struct spectra *psp,
	       struct output *pop,
	       ErrorMsg errmsg
	       ) {

  struct file_content fc;
  struct file_content fc_input;
  struct file_content fc_precision;

  char input_file[_ARGUMENT_LENGTH_MAX_];
  char precision_file[_ARGUMENT_LENGTH_MAX_];

  int i;
  char extension[5];

  /* Initialize some values. If no arguments are passed, they will remain null and
     inform init_params() that all parameters take default values. */

  fc.size = 0;
  fc_input.size = 0;
  fc_precision.size = 0;
  input_file[0]='\0';
  precision_file[0]='\0';

  /* If some arguments are passed, identify eventual 'xxx.ini' and 'xxx.pre' files, and store their name. */

  if (argc > 1) {
    for (i=1; i<argc; i++) {
      strncpy(extension,(argv[i]+strlen(argv[i])-4),4);
      extension[4]='\0';
      if (strcmp(extension,".ini") == 0) {
	class_test(input_file[0] != '\0',
		   errmsg,
		   "You have passed more than one input file with extension '.ini', choose one.");
	strcpy(input_file,argv[i]);
      }
      if (strcmp(extension,".pre") == 0) {
	class_test(precision_file[0] != '\0',
		   errmsg,
		   "You have passed more than one precision with extension '.pre', choose one.");
	strcpy(precision_file,argv[i]);
      }
    }
  }
  
  /* if there is an 'xxx.ini' file, read it and store its content. */

  if (input_file[0]!='\0')
    
    class_call(parser_read_file(input_file,&fc_input,errmsg),
		 errmsg,
		 errmsg);

  /* if there is an 'xxx.pre' file, read it and store its content. */

  if (precision_file[0]!='\0')
    
    class_call(parser_read_file(precision_file,&fc_precision,errmsg),
		 errmsg,
		 errmsg);

  /* if files were read, merge their contents in a single 'file_content' structure. */

  if ((input_file[0]!='\0') || (precision_file[0]!='\0'))

    class_call(parser_cat(&fc_input,&fc_precision,&fc,errmsg),
	       errmsg,
	       errmsg);

  class_call(parser_free(&fc_input),errmsg,errmsg);
  class_call(parser_free(&fc_precision),errmsg,errmsg);
  
  /* now, initialize all parameters given the input 'file_content' structure. 
     If its size is null, all parameters take their default values. */

  class_call(input_init(&fc,
			ppr,
			pba,
			pth,
			ppt,
			pbs,
			ptr,
			ppm,
			psp,
			pop,
			errmsg),
	     errmsg,
	     errmsg);
  
  class_call(parser_free(&fc),errmsg,errmsg);

  return _SUCCESS_;
}

/* Initialize all parameters given the input 'file_content' structure. 
   If its size is null, all parameters take their default values. */

int input_init(
	       struct file_content * pfc,
	       struct precision * ppr,
	       struct background *pba,
	       struct thermo *pth,
	       struct perturbs *ppt,
	       struct bessels * pbs,
	       struct transfers *ptr,
	       struct primordial *ppm,
	       struct spectra *psp,
	       struct output *pop,
	       ErrorMsg errmsg
	       ) {

  int flag1,flag2,flag3;
  double param1,param2,param3;
  int int1;
  double Omega_tot;
  char string1[_ARGUMENT_LENGTH_MAX_];
  char string2[_LINE_LENGTH_MAX_];
  int i;
  FILE * param_output;

  class_call(input_default_params(pba,
				  pth,
				  ppt,
				  pbs,
				  ptr,
				  ppm,
				  psp,
				  pop),
	     errmsg,
	     errmsg);

  class_call(input_default_precision(ppr),
	     errmsg,
	     errmsg);

  if (pfc->size == 0) 
    return _SUCCESS_;

  /* h (dimensionless) and H0 in Mpc^{-1} = h / 2999.7 */
  flag1=parser_read_double(pfc,"H0",&param1,errmsg);
  flag2=parser_read_double(pfc,"h",&param2,errmsg);
  class_test((flag1 == _SUCCESS_) && (flag2 == _SUCCESS_),
	     errmsg,
	     "In input file, you cannot enter both h and H0, choose one");
  if (flag1 == _SUCCESS_) {
    pba->H0 = param1 * 1.e3 / _c_;
    pba->h = param1 / 100.;
    }
  if (flag2 == _SUCCESS_) {
    pba->H0 = param2 *  1.e5 / _c_;
    pba->h = param2;
  }

  /* Omega_0_g (photons) and Tcmb */
  flag1=parser_read_double(pfc,"T_cmb",&param1,errmsg);
  flag2=parser_read_double(pfc,"Omega_g",&param2,errmsg);
  flag3=parser_read_double(pfc,"omega_g",&param3,errmsg);
  class_test(((flag1 == _SUCCESS_) && (flag2 == _SUCCESS_)) || ((flag1 == _SUCCESS_) && (flag3 == _SUCCESS_)) || ((flag2 == _SUCCESS_) && (flag3 == _SUCCESS_)),
	     errmsg,
	     "In input file, you can only enter one of Tcmb, Omega_g or omega_g, choose one");

  if ((flag1 == _FAILURE_) && (flag2 == _FAILURE_) && (flag3 == _FAILURE_)) {
    pba->Omega0_g = (4.*_sigma_B_/_c_*pow(pth->Tcmb,4.)) / (3.*_c_*_c_*1.e10*pba->h*pba->h/_Mpc_over_m_/_Mpc_over_m_/8./_PI_/_G_);
  }
  else {

    if (flag1 == _SUCCESS_) {
      /* Omega0_g = rho_g / rho_c0, each of them expressed in Kg/m/s^2 */
      /* rho_g = (4 sigma_B / c) T^4 */
      /* rho_c0 = 3 c^2 H0^2 / (8 pi G) */ 
      pba->Omega0_g = (4.*_sigma_B_/_c_*pow(param1,4.)) / (3.*_c_*_c_*1.e10*pba->h*pba->h/_Mpc_over_m_/_Mpc_over_m_/8./_PI_/_G_);
      pth->Tcmb=param1;
    }

    if (flag2 == _SUCCESS_) {
      pba->Omega0_g = param2;
      pth->Tcmb=pow(pba->Omega0_g * (3.*_c_*_c_*1.e10*pba->h*pba->h/_Mpc_over_m_/_Mpc_over_m_/8./_PI_/_G_) / (4.*_sigma_B_/_c_),0.25);
    }

    if (flag3 == _SUCCESS_) {
      pba->Omega0_g = param3/pba->h/pba->h;
      pth->Tcmb = pow(pba->Omega0_g * (3.*_c_*_c_*1.e10*pba->h*pba->h/_Mpc_over_m_/_Mpc_over_m_/8./_PI_/_G_) / (4.*_sigma_B_/_c_),0.25);
    }
  }

  Omega_tot = pba->Omega0_g;

  /* Omega_0_b (baryons) */
  flag1=parser_read_double(pfc,"Omega_b",&param1,errmsg);
  flag2=parser_read_double(pfc,"omega_b",&param2,errmsg);
  class_test(((flag1 == _SUCCESS_) && (flag2 == _SUCCESS_)),
	     errmsg,
	     "In input file, you can only enter one of Omega_b or omega_b, choose one");
  if (flag1 == _SUCCESS_)
    pba->Omega0_b = param1;
  if (flag2 == _SUCCESS_)
    pba->Omega0_b = param2/pba->h/pba->h;

  Omega_tot += pba->Omega0_b;

  /* Omega_0_nur (ultra-relativistic species / massless neutrino) */
  flag1=parser_read_double(pfc,"N_eff",&param1,errmsg);
  flag2=parser_read_double(pfc,"Omega_nur",&param2,errmsg);
  flag3=parser_read_double(pfc,"omega_nur",&param3,errmsg);
  class_test(((flag1 == _SUCCESS_) && (flag2 == _SUCCESS_)) || ((flag1 == _SUCCESS_) && (flag3 == _SUCCESS_)) || ((flag2 == _SUCCESS_) && (flag3 == _SUCCESS_)),
	     errmsg,
	     "In input file, you can only enter one of N_eff, Omega_nur or omega_nur, choose one");

  if ((flag1 == _FAILURE_) && (flag2 == _FAILURE_) && (flag3 == _FAILURE_)) {
    pba->Omega0_nur = 3.04*7./8.*pow(4./11.,4./3.)*pba->Omega0_g;
  }
  else {

    if (flag1 == _SUCCESS_) {
      pba->Omega0_nur = param1*7./8.*pow(4./11.,4./3.)*pba->Omega0_g;
    }
    if (flag2 == _SUCCESS_) {
      pba->Omega0_nur = param2;
    }
    if (flag3 == _SUCCESS_) {
      pba->Omega0_nur = param3/pba->h/pba->h;
    }
  }

  Omega_tot += pba->Omega0_nur;

  /* Omega_0_cdm (CDM) */
  flag1=parser_read_double(pfc,"Omega_cdm",&param1,errmsg);
  flag2=parser_read_double(pfc,"omega_cdm",&param2,errmsg);
  class_test(((flag1 == _SUCCESS_) && (flag2 == _SUCCESS_)),
	     errmsg,
	     "In input file, you can only enter one of Omega_cdm or omega_cdm, choose one");
  if (flag1 == _SUCCESS_)
    pba->Omega0_cdm = param1;
  if (flag2 == _SUCCESS_)
    pba->Omega0_cdm = param2/pba->h/pba->h;

  Omega_tot += pba->Omega0_cdm;

  /* Omega_0_k (curvature) */
  flag1=parser_read_double(pfc,"Omega_k",&param1,errmsg);
  if (flag1 == _SUCCESS_)
    pba->Omega0_k = param1;

  /* Omega_0_lambda (cosmological constant), Omega0_de (dark energy fluid) */
  flag1=parser_read_double(pfc,"Omega_Lambda",&param1,errmsg);
  flag2=parser_read_double(pfc,"Omega_de",&param2,errmsg);
  class_test((flag1 == _SUCCESS_) && (flag2 == _SUCCESS_),
	     errmsg,
	     "In input file, you can enter only two out of Omega_Lambda, Omega_de, Omega_k, the third one is inferred");

  if ((flag1 == _FAILURE_) && (flag2 == _FAILURE_)) {	
    pba->Omega0_lambda = 1.+pba->Omega0_k-pba->Omega0_g-pba->Omega0_nur-pba->Omega0_b-pba->Omega0_cdm;
  }
  else {
    if (flag1 == _SUCCESS_) {
      pba->Omega0_lambda= param1;
      pba->Omega0_de = 1. + pba->Omega0_k - param1 - Omega_tot;
    }
    if (flag2 == _SUCCESS_) {
      pba->Omega0_lambda= 1. + pba->Omega0_k - param2 - Omega_tot;
      pba->Omega0_de = param2;
    }
  }

  class_test(pba->Omega0_de != 0.,
	     errmsg,
	     "Dark energy fluid not tested yet");
  
  class_test(pba->Omega0_k != 0.,
	     errmsg,
	     "Open/close case not written yet");

  /* scale factor today (arbitrary) */
  flag1=parser_read_double(pfc,"a_today",&param1,errmsg);
  if (flag1 == _SUCCESS_)
    pba->a_today = param1;

  /** - assign values to thermodynamics cosmological parameters */

  /* scale factor today (arbitrary) */
  flag1=parser_read_double(pfc,"YHe",&param1,errmsg);
  if (flag1 == _SUCCESS_)
    pth->YHe = param1;

  /* reionization parametrization */
  flag1=parser_read_string(pfc,"reio_parametrization",&string1,errmsg);
  if (flag1 == _SUCCESS_) {
    flag2=_FALSE_;
    if (strcmp(string1,"reio_none") == 0) {
      pth->reio_parametrization=reio_none;
      flag2=_TRUE_;
      printf("Warning: you are computing a model without reionization (why not...)\n");
    }
    if (strcmp(string1,"reio_camb") == 0) {
      pth->reio_parametrization=reio_camb;
      flag2=_TRUE_;
    }
    class_test(flag2==_FALSE_,
	       errmsg,
	       "could not identify reionization_parametrization value, check that it is one of 'reio_none', 'reio_camb', ...");
  }

  /* reionization parameters if reio_parametrization=reio_camb */
  if (pth->reio_parametrization == reio_camb) {
    flag1=parser_read_double(pfc,"z_reio",&param1,errmsg);
    flag2=parser_read_double(pfc,"tau_reio",&param2,errmsg);
    class_test(((flag1 == _SUCCESS_) && (flag2 == _SUCCESS_)),
	       errmsg,
	       "In input file, you can only enter one of z_reio or tau_reio, choose one");
    if (flag1 == _SUCCESS_) {
      pth->z_reio=param1;
      pth->reio_z_or_tau=reio_z;
    }
    if (flag2 == _SUCCESS_) {
      pth->tau_reio=param2;
      pth->reio_z_or_tau=reio_tau;
    }
  }

  /** - define which perturbations and sources should be computed, and down to which scale */

  flag1=parser_read_string(pfc,"output",&string1,errmsg);
  if (flag1 == _SUCCESS_) {

    if ((strstr(string1,"tCl") != NULL) || (strstr(string1,"TCl") != NULL) || (strstr(string1,"TCL") != NULL))
      ppt->has_cl_cmb_temperature=_TRUE_;  

    if ((strstr(string1,"pCl") != NULL) || (strstr(string1,"PCl") != NULL) || (strstr(string1,"PCL") != NULL))
      ppt->has_cl_cmb_polarization=_TRUE_;  
    
    if ((strstr(string1,"lCl") != NULL) || (strstr(string1,"LCl") != NULL) || (strstr(string1,"LCL") != NULL))
      ppt->has_cl_cmb_lensing_potential=_TRUE_;

    if ((strstr(string1,"mPk") != NULL) || (strstr(string1,"MPk") != NULL) || (strstr(string1,"MPK") != NULL))
      ppt->has_pk_matter=_TRUE_; 
  }

  if ((ppt->has_cl_cmb_temperature == _TRUE_) ||
      (ppt->has_cl_cmb_polarization == _TRUE_) ||
      (ppt->has_cl_cmb_lensing_potential == _TRUE_) ||
      (ppt->has_pk_matter == _TRUE_)) {

    flag1=parser_read_string(pfc,"modes",&string1,errmsg);
    if (flag1 == _SUCCESS_) {

      /* if no modes are specified, the default is has_scalars=_TRUE_; 
	 but if they are specified we should reset has_scalars to _FALSE_ before reading */
      ppt->has_scalars=_FALSE_;

      if ((strstr(string1,"s") != NULL) || (strstr(string1,"S") != NULL))
	ppt->has_scalars=_TRUE_; 

      if ((strstr(string1,"v") != NULL) || (strstr(string1,"V") != NULL)) {
	ppt->has_vectors=_TRUE_;  
      }
      
      if ((strstr(string1,"t") != NULL) || (strstr(string1,"T") != NULL)) {
	ppt->has_tensors=_TRUE_;
      }

      class_test(ppt->has_scalars==_FALSE_ && ppt->has_vectors ==_FALSE_ && ppt->has_tensors ==_FALSE_,
		 errmsg,	       
		 "You wrote: modes=%s. Could not identify any of the modes ('s', 'v', 't') in such input",string1);
    }

    if (ppt->has_scalars == _TRUE_) {

      flag1=parser_read_string(pfc,"ic",&string1,errmsg);
      if (flag1 == _SUCCESS_) {

	/* if no initial conditions are specified, the default is has_ad=_TRUE_; 
	   but if they are specified we should reset has_ad to _FALSE_ before reading */
	ppt->has_ad=_FALSE_;

	if ((strstr(string1,"ad") != NULL) || (strstr(string1,"AD") != NULL))
	  ppt->has_ad=_TRUE_; 
	
	if ((strstr(string1,"bi") != NULL) || (strstr(string1,"BI") != NULL))
	  ppt->has_bi=_TRUE_; 
	
	if ((strstr(string1,"cdi") != NULL) || (strstr(string1,"CDI") != NULL))
	  ppt->has_cdi=_TRUE_; 
	
	if ((strstr(string1,"nid") != NULL) || (strstr(string1,"NID") != NULL))
	  ppt->has_nid=_TRUE_; 
	
	if ((strstr(string1,"niv") != NULL) || (strstr(string1,"NIV") != NULL))
	  ppt->has_niv=_TRUE_; 
      
	class_test(ppt->has_ad==_FALSE_ && ppt->has_bi ==_FALSE_ && ppt->has_cdi ==_FALSE_ && ppt->has_nid ==_FALSE_ && ppt->has_niv ==_FALSE_,
		   errmsg,	       
		   "You wrote: ic=%s. Could not identify any of the initial conditions ('ad', 'bi', 'cdi', 'nid', 'niv') in such input",string1);
	
      }
    }

    else {

      class_test(ppt->has_cl_cmb_lensing_potential == _TRUE_,
		 errmsg,
		 "Inconsistency: you want C_l's for cmb lensing potential, but no scalar modes\n");

      class_test(ppt->has_pk_matter == _TRUE_,
		 errmsg,
		 "Inconsistency: you want P(k) of matter, but no scalar modes\n");

    }

  }

  /** - define the primordial spectrum */

  flag1=parser_read_string(pfc,"P_k_ini type",&string1,errmsg);
  if (flag1 == _SUCCESS_) {
    flag2=_FALSE_;
    if (strcmp(string1,"analytic_Pk") == 0) {
      ppm->primordial_spec_type = analytic_Pk;
      flag2=_TRUE_;
    }
    class_test(flag2==_FALSE_,
	       errmsg,
	       "could not identify primordial spectrum type, check that it is one of 'analytic_pk', ...");
  }

  if (ppm->primordial_spec_type == analytic_Pk) {

    class_read_double("k_pivot",ppm->k_pivot);

    if (ppt->has_scalars == _TRUE_) {
      
      if (ppt->has_ad == _TRUE_) {

	class_read_double("A_s_ad",ppm->A_s_ad);
	class_read_double("n_s_ad",ppm->n_s_ad);
	class_read_double("alpha_s_ad",ppm->alpha_s_ad);

      }

    }

    if (ppt->has_tensors == _TRUE_) {
    
      	class_read_double("A_t",ppm->A_t);
	class_read_double("n_t",ppm->n_t);
	class_read_double("alpha_t",ppm->alpha_t);

    }

  }

  /** - parameters for output spectra */

  flag1=parser_read_string(pfc,"root",&string1,errmsg);
  if (flag1 == _SUCCESS_) {
    sprintf(pop->cls_ad,"%s%s",string1,"cls_ad.dat");
    sprintf(pop->cls_bi,"%s%s",string1,"cls_bi.dat");
    sprintf(pop->cls_cdi,"%s%s",string1,"cls_cdi.dat");
    sprintf(pop->cls_nid,"%s%s",string1,"cls_nid.dat");
    sprintf(pop->cls_niv,"%s%s",string1,"cls_niv.dat");
    sprintf(pop->clt,"%s%s",string1,"clt.dat");
    sprintf(pop->cltot,"%s%s",string1,"cltot.dat");
    sprintf(pop->pk,"%s%s",string1,"pk.dat");
  }

  pbs->l_max=0;

  if (ppt->has_scalars == _TRUE_) {
    class_read_double("l_max_scalars",ptr->l_scalar_max);
    pbs->l_max=ptr->l_scalar_max;
  }

  if (ppt->has_tensors == _TRUE_) {   
    class_read_double("l_max_tensors",ptr->l_tensor_max);
    pbs->l_max=max(pbs->l_max,ptr->l_tensor_max);
  }

  class_read_double("z_pk",pop->z_pk);

  class_read_double("P_k_max",ppt->k_scalar_kmax_for_pk);

  flag1=parser_read_double(pfc,"z_max_pk",&param1,errmsg);
  if (flag1==_SUCCESS_) {
    psp->z_max_pk = param1;
  }
  else {
    psp->z_max_pk = pop->z_pk;
  }

  /** - amount of information sent to standard output (none if all set to zero) */

  class_read_int("background_verbose",
		 pba->background_verbose);

  class_read_int("thermodynamics_verbose",
		 pth->thermodynamics_verbose);

  class_read_int("perturbations_verbose",
		 ppt->perturbations_verbose);

  class_read_int("bessels_verbose",
		 pbs->bessels_verbose);

  class_read_int("transfer_verbose",
		 ptr->transfer_verbose);

  class_read_int("primordial_verbose",
		 ppm->primordial_verbose);

  class_read_int("spectra_verbose",
		 psp->spectra_verbose);

  class_read_int("output_verbose",
		 pop->output_verbose);

  /** Precision parameters */

  /** - parameters related to the background */

  class_read_double("a_ini_over_a_today_default",ppr->a_ini_over_a_today_default);
  class_read_double("back_integration_stepsize",ppr->back_integration_stepsize);
  class_read_double("tol_background_integration",ppr->tol_background_integration);

  /** - parameters related to the thermodynamics */

  class_read_double("recfast_z_initial",ppr->recfast_z_initial);
  class_read_double("recfast_z_final",ppr->recfast_z_final);
  class_read_double("recfast_H_frac",ppr->recfast_H_frac);
  class_read_double("recfast_x_H0_trigger",ppr->recfast_x_H0_trigger);
  class_read_double("recfast_x_He0_trigger",ppr->recfast_x_He0_trigger);
  class_read_double("recfast_fudge",ppr->recfast_fudge);
  class_read_double("recfast_fudge_He",ppr->recfast_fudge_He);
  class_read_int("recfast_Heswitch",ppr->recfast_Heswitch);
  class_read_int("recfast_Nz0",ppr->recfast_Nz0);
  class_read_double("tol_thermo_integration",ppr->tol_thermo_integration);
  class_read_double("visibility_threshold_start_sources",ppr->visibility_threshold_start_sources);
  class_read_double("visibility_threshold_free_streaming",ppr->visibility_threshold_free_streaming);
  class_read_double("reionization_z_start_max",ppr->reionization_z_start_max);
  class_read_double("reionization_sampling",ppr->reionization_sampling);
  class_read_double("reionization_optical_depth_tol",ppr->reionization_optical_depth_tol);
  class_read_double("reionization_exponent",ppr->reionization_exponent);
  class_read_double("reionization_width",ppr->reionization_width);
  class_read_double("reionization_start_factor",ppr->reionization_start_factor);
  class_read_double("helium_fullreio_redshift",ppr->helium_fullreio_redshift);
  class_read_double("helium_fullreio_width",ppr->helium_fullreio_width);
  class_read_int("thermo_rate_smoothing_radius",ppr->thermo_rate_smoothing_radius);

  /** - parameters related to the perturbations */

  class_read_int("gauge",ppr->gauge);
  class_read_double("k_scalar_min",ppr->k_scalar_min);
  class_read_double("k_scalar_oscillations",ppr->k_scalar_oscillations);
  class_read_double("k_scalar_step_sub",ppr->k_scalar_step_sub);
  class_read_double("k_scalar_step_super",ppr->k_scalar_step_super);
  class_read_double("k_scalar_step_transition",ppr->k_scalar_step_transition);
  class_read_double("k_scalar_k_per_decade_for_pk",ppr->k_scalar_k_per_decade_for_pk);
  class_read_double("k_tensor_min",ppr->k_tensor_min);
  class_read_double("k_tensor_oscillations",ppr->k_tensor_oscillations);
  class_read_double("k_tensor_step_sub",ppr->k_tensor_step_sub);
  class_read_double("k_tensor_step_super",ppr->k_tensor_step_super);
  class_read_double("k_tensor_step_transition",ppr->k_tensor_step_transition);
  class_read_double("k_eta_min",ppr->k_eta_min);
  class_read_double("eta_min_over_sampling_min",ppr->eta_min_over_sampling_min);
  class_read_double("k_eta_ma",ppr->k_eta_max);
  class_read_int("l_max_g",ppr->l_max_g);
  class_read_int("l_max_pol_g",ppr->l_max_pol_g);
  class_read_int("l_max_nur",ppr->l_max_nur);
  class_read_int("l_max_g_ten",ppr->l_max_g_ten);
  class_read_int("l_max_pol_g_ten",ppr->l_max_pol_g_ten);
  class_read_double("phi_ini",ppr->phi_ini);
  class_read_double("entropy_ini",ppr->entropy_ini);
  class_read_double("gw_ini",ppr->gw_ini);
  class_read_double("perturb_integration_stepsize",ppr->perturb_integration_stepsize);
  class_read_double("tol_perturb_integration",ppr->tol_perturb_integration);
  class_read_double("perturb_sampling_stepsize",ppr->perturb_sampling_stepsize);
  class_read_double("tight_coupling_trigger_eta_g_over_eta_h",ppr->tight_coupling_trigger_eta_g_over_eta_h);
  class_read_double("tight_coupling_trigger_eta_g_over_eta_",ppr->tight_coupling_trigger_eta_g_over_eta_k);
  class_read_double("rad_pert_trigger_k_over_aH",ppr->rad_pert_trigger_k_over_aH);
  class_read_double("rad_pert_trigger_Omega_r",ppr->rad_pert_trigger_Omega_r);

  /** - parameter related to the Bessel functions */

  class_read_double("l_logstep",ppr->l_logstep);
  class_read_int("l_linstep",ppr->l_linstep);
  class_read_double("bessel_scalar_x_step",ppr->bessel_scalar_x_step);
  class_read_double("bessel_scalar_j_cut",ppr->bessel_scalar_j_cut);
  class_read_int("bessel_always_recompute",ppr->bessel_always_recompute);

  /** - parameter related to the primordial spectra */

  class_read_double("k_per_decade_primordial",ppr->k_per_decade_primordial);

  /** - parameter related to the transfer functions */

  class_read_double("k_step_trans_scalars",ppr->k_step_trans_scalars);
  class_read_double("k_step_trans_tensors",ppr->k_step_trans_tensors);
  class_read_double("transfer_cut",ppr->transfer_cut);
  class_read_double("transfer_cut_threshold_osc",ppr->transfer_cut_threshold_osc);
  class_read_double("transfer_cut_threshold_cl",ppr->transfer_cut_threshold_cl);

  flag1=parser_read_string(pfc,"parameters",&string1,errmsg);
  if (flag1 == _SUCCESS_) {
    class_open(param_output,string1,"w",errmsg);
    fprintf(param_output,"# List of input/precision parameters actually read\n");
    fprintf(param_output,"# (all other parameters set to default values)\n");
    fprintf(param_output,"#\n");
    fprintf(param_output,"# This file, written by CLASS, can be used as the input file\n");
    fprintf(param_output,"# of another run provided that it has a suffix xxx.ini or xxx.pre;\n");
    fprintf(param_output,"# in the next run you can change this name with the\n");
    fprintf(param_output,"# parameters = ... entry of your input file.\n");
    fprintf(param_output,"#\n");
    for (i=0; i<pfc->size; i++) {
      if (pfc->read[i] == _TRUE_)
	fprintf(param_output,"%s = %s\n",pfc->name[i],pfc->value[i]);
    }
    fprintf(param_output,"#\n");
  }

  return _SUCCESS_;

}

int input_default_params(
			 struct background *pba,
			 struct thermo *pth,
			 struct perturbs *ppt,
			 struct bessels * pbs,
			 struct transfers *ptr,
			 struct primordial *ppm,
			 struct spectra *psp,
			 struct output *pop
			 ) {

  ErrorMsg errmsg;
         
  pba->h = 0.7;
  pba->H0 = pba->h * 1.e5 / _c_;
  pth->Tcmb = 2.726;
  pba->Omega0_g = (4.*_sigma_B_/_c_*pow(pth->Tcmb,4.)) / (3.*_c_*_c_*1.e10*pba->h*pba->h/_Mpc_over_m_/_Mpc_over_m_/8./_PI_/_G_);
  pba->Omega0_nur = 3.04*7./8.*pow(4./11.,4./3.)*pba->Omega0_g;
  pba->Omega0_b = 0.05;
  pba->Omega0_cdm = 0.25;    
  pba->Omega0_k = 0.;
  pba->Omega0_lambda = 1.+pba->Omega0_k-pba->Omega0_g-pba->Omega0_nur-pba->Omega0_b-pba->Omega0_cdm;
  pba->Omega0_de = 0.;     
  pba->a_today = 1.;       
  pba->w_de=-1.;
  pba->cs2_de=1.;

  /* pth->Tcmb already fixed above */
  pth->YHe=0.25;            
  pth->reio_parametrization=reio_camb;
  pth->reio_z_or_tau=reio_z;
  pth->z_reio=10.;
  pth->tau_reio=0.08;

  ppt->has_cl_cmb_temperature = _FALSE_;
  ppt->has_cl_cmb_polarization = _FALSE_;
  ppt->has_cl_cmb_lensing_potential = _FALSE_;
  ppt->has_pk_matter = _FALSE_;

  ppt->has_ad=_TRUE_;  
  ppt->has_bi=_FALSE_;
  ppt->has_cdi=_FALSE_;
  ppt->has_nid=_FALSE_;
  ppt->has_niv=_FALSE_;

  ppt->has_scalars=_TRUE_;  
  ppt->has_vectors=_FALSE_;
  ppt->has_tensors=_FALSE_;  

  ppm->primordial_spec_type = analytic_Pk;
  ppm->k_pivot = 0.05;
  ppm->A_s_ad = 2.3e-9;
  ppm->n_s_ad = 1.;
  ppm->alpha_s_ad = 0.;
  ppm->A_t = 2.3e-9;
  ppm->n_t = 0.;
  ppm->alpha_t = 0.;

  sprintf(pop->cls_ad,"output/cls.dat");
  sprintf(pop->clt,"output/clt.dat");
  pbs->l_max=2500;
  ptr->l_scalar_max=2500;
  ptr->l_tensor_max=500;

  sprintf(pop->pk,"output/pk.dat");
  pop->z_pk = 0.;  
  psp->z_max_pk = pop->z_pk;
  
  pba->background_verbose = 0;
  pth->thermodynamics_verbose = 0;
  ppt->perturbations_verbose = 0;
  pbs->bessels_verbose = 0;
  ptr->transfer_verbose = 0;
  ppm->primordial_verbose = 0;
  psp->spectra_verbose = 0;
  pop->output_verbose = 0;

  return _SUCCESS_;

}

/** 
 * Initialize the precision parameter structure. 
 * 
 * @param ppr Input/Ouput: a precision_params structure pointer  
 * @return the error status
 *
 * All precision parameters used in the other moduels are assigned
 * here a default value, namely:
 */
int input_default_precision ( struct precision * ppr ) {

  /** Summary: */

  /**
   * - parameters related to the background
   */

  ppr->a_ini_over_a_today_default = 1.e-13;  /* 1.e-7 unless needs large k_max in P(k) */
  ppr->back_integration_stepsize = 2.e-2;   /* 0.02 */
  ppr->tol_background_integration = 1.e-3;  /* 0.002 */

  /**
   * - parameters related to the thermodynamics
   */

  ppr->recfast_z_initial=1.e4;
  ppr->recfast_z_final=0.;
  ppr->recfast_H_frac=1.e-3; /* from recfast */      
  ppr->recfast_x_H0_trigger=0.99; /* from recfast */   
  ppr->recfast_x_He0_trigger=0.99; /* from recfast */  
  ppr->recfast_fudge=1.14; /* from recfast */
  ppr->recfast_fudge_He=0.86; /* from recfast 1.4 */
  ppr->recfast_Heswitch=6.; /* from recfast 1.4 */ 
  ppr->recfast_Nz0=10000; /* smaller than 6000 gives bug in transfer, need to check why */
  ppr->tol_thermo_integration=1.e-3; /* optimized 9/09/08  */

  ppr->visibility_threshold_start_sources=3.5e-7; /* 3.5e-7 optimized 9/09/08  */
  ppr->visibility_threshold_free_streaming=1.e-5;

  ppr->reionization_z_start_max = 50.;
  ppr->reionization_sampling=1.e-2; /*1.e-2*/
  ppr->reionization_optical_depth_tol=1.e-2;
  ppr->reionization_exponent=1.5;
  ppr->reionization_width=0.5;
  ppr->reionization_start_factor=8.;
  ppr->helium_fullreio_redshift=3.5;
  ppr->helium_fullreio_width=0.5;

  ppr->thermo_rate_smoothing_radius=10;

  /**
   * - parameters related to the perturbations
   */

  ppr->gauge=synchronous;

  ppr->k_scalar_min=0.3; /* 0.3 -> 0.1 */
  ppr->k_scalar_oscillations=7.;  
  ppr->k_scalar_step_sub=0.1;  /* 0.02 -> 0.005 */
  ppr->k_scalar_step_super=0.005;  /* 0.01 -> 0.005 */
  ppr->k_scalar_step_transition=0.4;

  ppr->k_scalar_k_per_decade_for_pk=10.;

  ppr->k_tensor_min=0.1; /* 0.3 -> 0.1 */
  ppr->k_tensor_oscillations=3.5;  
  ppr->k_tensor_step_sub=0.01;  /* 0.02 -> 0.005 */
  ppr->k_tensor_step_super=0.0002;  /* 0.01 -> 0.005 */
  ppr->k_tensor_step_transition=0.2;

  ppr->k_eta_min=1.e-1; /* 4.5e-6 optimized 9/09/08  */
  ppr->eta_min_over_sampling_min=0.5;
  ppr->k_eta_max=10.; /* 600 */

  ppr->l_max_g=10; /* optimized 9/09/08  */
  ppr->l_max_pol_g=10; /* optimized 9/09/08  */
  ppr->l_max_nur=25;
  ppr->l_max_g_ten=5;
  ppr->l_max_pol_g_ten=5;

  ppr->phi_ini=1.;
  ppr->entropy_ini=1.;
  ppr->gw_ini=1.;

  ppr->perturb_integration_stepsize=0.5; /* 0.5 */ 
  ppr->tol_perturb_integration=1.e-3; 
  ppr->perturb_sampling_stepsize=0.1; /* 0.1 */

  ppr->tight_coupling_trigger_eta_g_over_eta_h=0.006; /* 0.006 */
  ppr->tight_coupling_trigger_eta_g_over_eta_k=1.5e-2; /*1.5e-2*/

  ppr->rad_pert_trigger_k_over_aH = 40.; /* 40 */
  ppr->rad_pert_trigger_Omega_r = 0.1; /* 0.1 */

  /**
   * - parameter related to the Bessel functions
   */

  ppr->l_logstep=1.2 /* 1.4*/;
  ppr->l_linstep=50;

  ppr->bessel_scalar_x_step=0.1; /* 1. 1.27 optimized 9/09/08 */
  ppr->bessel_scalar_j_cut=1.e-5; /* 8.1e-5 optimized 9/09/08 */
  ppr->bessel_always_recompute=_TRUE_;

  /**
   * - parameter related to the primordial spectra
   */

  ppr->k_per_decade_primordial = 10.; 

  /**
   * - parameter related to the transfer functions
   */
  
  ppr->k_step_trans_scalars=0.15; /* 0.1 sampling step in k space, in units of 2pi/(eta_0-eta_rec), which is the typical period of oscillations of |Delta_l(k)|^2 */
  ppr->k_step_trans_tensors=0.15;

  ppr->transfer_cut=tc_cl;
  ppr->transfer_cut_threshold_osc=0.01; /* 0.01 */
  ppr->transfer_cut_threshold_cl=2.e-6; /* 2.e-6 */

  /**
   * - automatic estimate of machine precision
   */

  get_machine_precision(&(ppr->smallest_allowed_variation));

  class_test(ppr->smallest_allowed_variation < 0,
	     ppr->error_message,
	     "smallest_allowed_variation = %e < 0",ppr->smallest_allowed_variation);

  return _SUCCESS_;

}

/** 
 * Computes automatically the machine precision. 
 *
 * @param smallest_allowed_variation a pointer to the smallest allowed variation
 *
 * Returns the smallest
 * allowed variation (minimum epsilon * _TOLVAR_)
 */
int get_machine_precision(double * smallest_allowed_variation) {
  double one, meps, sum;
  
  one = 1.0;
  meps = 1.0;
  do {
    meps /= 2.0;
    sum = one + meps;
  } while (sum != one);
  meps *= 2.0;
  
  *smallest_allowed_variation = meps * _TOLVAR_;

  return _SUCCESS_;

}
