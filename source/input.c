/** @file input.c 
 * Julien Lesgourgues, 18.04.2010    
 */

#include "input.h" 

/* If class is executed in a terminal, use this routine to extract initial parameters 
   from the arguments of the main.
   If instead class is embedded into another code, use directly input_init_params() to
   pass input parameters through a 'file_content' structure.
 */

int input_init(
	       int argc, 
	       char **argv,
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
  char input_file[_ARGUMENT_LENGTH_MAX_];
  char precision_file[_ARGUMENT_LENGTH_MAX_];

  /* if input_init is called with arg=0 or 1, set null input file names in order to read parameters from code */
  if (argc < 2) {
    input_file[0]=='\0';
    precision_file[0]=='\0';
  }
  else {
    class_call(input_check_arguments_of_main(argc, argv, input_file,precision_file,errmsg),
	       errmsg,
	       errmsg);
  }

  if (input_file[0]=='\0') {

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

  }
  else {

    class_call(parser_read_file(input_file,&fc,errmsg),
	       errmsg,
	       errmsg);

    class_call(input_init_params(&fc,
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
  }

  return _SUCCESS_;
}

  int input_init_params(
			struct file_content * pfc,
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
    
    if ((strstr(string1,"lCl") != NULL) || (strstr(string1,"LCl") != NULL) || (strstr(string1,"LCL") == NULL))
      ppt->has_cl_cmb_lensing_potential=_TRUE_;  

    if ((strstr(string1,"mPk") != NULL) || (strstr(string1,"MPk") != NULL) || (strstr(string1,"MPK") == NULL))
      ppt->has_pk_matter=_TRUE_;  
  }

  if ((ppt->has_cl_cmb_temperature == _TRUE_) ||
      (ppt->has_cl_cmb_polarization == _TRUE_) ||
      (ppt->has_cl_cmb_lensing_potential == _TRUE_) ||
      (ppt->has_pk_matter == _TRUE_)) {

    flag1=parser_read_string(pfc,"modes",&string1,errmsg);
    if (flag1 == _SUCCESS_) {
    
      if ((strstr(string1,"s") != NULL) || (strstr(string1,"S") != NULL))
	ppt->has_scalars=_TRUE_;  
      
      if ((strstr(string1,"v") != NULL) || (strstr(string1,"V") != NULL))
	ppt->has_vectors=_TRUE_;  
      
      if ((strstr(string1,"t") != NULL) || (strstr(string1,"T") != NULL))
	ppt->has_tensors=_TRUE_;  
      
      class_test(ppt->has_scalars==_FALSE_ && ppt->has_vectors ==_FALSE_ && ppt->has_tensors ==_FALSE_,
		 errmsg,	       
		 "You wrote: modes=%s. Could not identify any of the modes ('s', 'v', 't') in such input",string1);
    }

    if (ppt->has_scalars == _TRUE_) {

      flag1=parser_read_string(pfc,"ic",&string1,errmsg);
      if (flag1 == _SUCCESS_) {

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

    class_read_double("A_s_ad",ppm->A_s_ad);

    class_read_double("n_s_ad",ppm->n_s_ad);

    class_read_double("alpha_s_ad",ppm->alpha_s_ad);

    class_read_double("k_pivot",ppm->k_pivot);

  }

  /** - parameters for output spectra */

  class_read_string("cls",pop->cls_ad);

  class_read_double("l_max",pbs->l_max);
  ptr->l_scalar_max=pbs->l_max;
  ptr->l_tensor_max=pbs->l_max;

  class_read_string("pk",pop->pk);

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
  ppm->A_s_ad = 2.3e-9;
  ppm->n_s_ad = 1.;
  ppm->alpha_s_ad = 0.;
  ppm->k_pivot = 0.05;

  sprintf(pop->cls_ad,"output/cls.dat");
  pbs->l_max=2500;
  ptr->l_scalar_max=2500;
  ptr->l_tensor_max=2500;

  sprintf(pop->pk,"output/pk.dat");
  pop->z_pk = 0.;  
  psp->z_max_pk = pop->z_pk;
  
  pth->thermodynamics_verbose = 0;
  ppt->perturbations_verbose = 0;
  pbs->bessels_verbose = 0;
  ptr->transfer_verbose = 0;
  ppm->primordial_verbose = 0;
  psp->spectra_verbose = 0;
  pop->output_verbose = 0;

  return _SUCCESS_;

}

int input_check_arguments_of_main(
				  int argc, 
				  char **argv, 
				  char * input,
				  char * precision,
				  ErrorMsg errmsg) {

  int i;
  char extension[5];

  if (argc == 1) {
    input[0]='\0';
    precision[0]='\0';
    return _SUCCESS_;
  }

  input[0]='\0';
  precision[0]='\0';
  for (i=1; i<argc; i++) {
    strncpy(extension,(argv[i]+strlen(argv[i])-4),4);
    extension[4]='\0';
    if (strcmp(extension,".ini") == 0) {
      class_test(input[0] != '\0',
		 errmsg,
		 "You have passed more than one input file with extension '.ini', choose one.");
      strcpy(input,argv[i]);
    }
    if (strcmp(extension,".pre") == 0) {
      class_test(precision[0] != '\0',
		 errmsg,
		 "You have passed more than one precision with extension '.pre', choose one.");
      strcpy(precision,argv[i]);
    }
  }

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
