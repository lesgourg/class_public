/** @file input.c 
 * Julien Lesgourgues, 18.04.2010    
 */

#include "input.h" 

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

  int flag1,flag2,flag3;
  double param1,param2,param3;
  struct file_content fc;
  char input_file[_ARGUMENT_LENGTH_MAX_];
  double Omega_tot;
  char string1[_ARGUMENT_LENGTH_MAX_];
  char string2[_LINE_LENGTH_MAX_];

  /** - assign values to background cosmological parameters */
  /* the following parameters must be assigned:

     pba->H0
     pba->Omega0_g
     pba->Omega0_nur
     pba->Omega0_b
     pba->Omega0_cdm     (optional; 0. if not passed)
     pba->Omega0_lambda  (optional; 0. if not passed)
     pba->Omega0_de      (optional; 0. if not passed)
     pba->a_today        (optional; 1. if not passed)
     pth->Tcmb
     pth->YHe            (optional; 0.25 if not passed)
     pth->reio_parametrization
     pth->reio_z_or_tau
     pth->z_reio
     pth->tau_reio
  */

  class_call(input_check_arguments_of_main(argc, argv, input_file,errmsg),
	     errmsg,
	     errmsg);

  class_call(parser_read_file(input_file,&fc,errmsg),
	     errmsg,
	     errmsg);

  /* h (dimensionless) and H0 in Mpc^{-1} = h / 2999.7 */
  flag1=parser_read_double(&fc,"H0",&param1,errmsg);
  flag2=parser_read_double(&fc,"h",&param2,errmsg);
  class_test((flag1 == _FAILURE_) && (flag2 == _FAILURE_),
	     errmsg,
	     "In input file, enter either h (dimensionless) or H0 (in km/s/Mpc)");
  class_test((flag1 == _SUCCESS_) && (flag2 == _SUCCESS_),
	     errmsg,
	     "In input file, you cannot enter both h and H0, choose one");
  if (flag1 == _SUCCESS_) {
    pba->H0 = param1 * 1.e3 / _c_;
    pba->h = param1 / 100.;
    }
  else {
    pba->H0 = param2 *  1.e5 / _c_;
    pba->h = param2;
  }


  /* Omega_0_g (photons) and Tcmb */
  flag1=parser_read_double(&fc,"T_cmb",&param1,errmsg);
  flag2=parser_read_double(&fc,"Omega_g",&param2,errmsg);
  flag3=parser_read_double(&fc,"omega_g",&param3,errmsg);
  if((flag1 == _FAILURE_) && (flag2 == _FAILURE_) && (flag3 == _FAILURE_)) {
    pth->Tcmb=2.726;
    pba->Omega0_g = (4.*_sigma_B_/_c_*pow(pth->Tcmb,4.)) / (3.*_c_*_c_*1.e10*pba->h*pba->h/_Mpc_over_m_/_Mpc_over_m_/8./_PI_/_G_);
  }
  class_test(((flag1 == _SUCCESS_) && (flag2 == _SUCCESS_)) || ((flag1 == _SUCCESS_) && (flag3 == _SUCCESS_)) || ((flag2 == _SUCCESS_) && (flag3 == _SUCCESS_)),
	     errmsg,
	     "In input file, you can only enter one of Tcmb, Omega_g or omega_g, choose one");
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

  Omega_tot = pba->Omega0_g;

  /* Omega_0_b (baryons) */
  flag1=parser_read_double(&fc,"Omega_b",&param1,errmsg);
  flag2=parser_read_double(&fc,"omega_b",&param2,errmsg);
  class_test((flag1 == _FAILURE_) && (flag2 == _FAILURE_),
	     errmsg,
	     "In input file, enter either Omega_b or omega_b (to infer baryon density)");
  class_test(((flag1 == _SUCCESS_) && (flag2 == _SUCCESS_)),
	     errmsg,
	     "In input file, you can only enter one of Omega_b or omega_b, choose one");
  if (flag1 == _SUCCESS_)
    pba->Omega0_b = param1;
  else
    pba->Omega0_b = param2/pba->h/pba->h;

  Omega_tot += pba->Omega0_b;

  /* Omega_0_nur (ultra-relativistic species / massless neutrino) */
  flag1=parser_read_double(&fc,"N_eff",&param1,errmsg);
  flag2=parser_read_double(&fc,"Omega_nur",&param2,errmsg);
  flag3=parser_read_double(&fc,"omega_nur",&param3,errmsg);
  if((flag1 == _FAILURE_) && (flag2 == _FAILURE_) && (flag3 == _FAILURE_)) {
    pba->Omega0_nur = 3.04*7./8.*pow(4./11.,4./3.)*pba->Omega0_g;
  }
  class_test(((flag1 == _SUCCESS_) && (flag2 == _SUCCESS_)) || ((flag1 == _SUCCESS_) && (flag3 == _SUCCESS_)) || ((flag2 == _SUCCESS_) && (flag3 == _SUCCESS_)),
	     errmsg,
	     "In input file, you can only enter one of N_eff, Omega_nur or omega_nur, choose one");
  if (flag1 == _SUCCESS_) {
    pba->Omega0_nur = param1*7./8.*pow(4./11.,4./3.)*pba->Omega0_g;
  }
  if (flag2 == _SUCCESS_) {
    pba->Omega0_nur = param2;
  }
  if (flag3 == _SUCCESS_) {
    pba->Omega0_nur = param3/pba->h/pba->h;
  }

  Omega_tot += pba->Omega0_nur;

  /* Omega_0_cdm (CDM) */
  flag1=parser_read_double(&fc,"Omega_cdm",&param1,errmsg);
  flag2=parser_read_double(&fc,"omega_cdm",&param2,errmsg);
  if ((flag1 == _FAILURE_) && (flag2 == _FAILURE_)) {
    printf("Warning: you are computing a model without Cold Dark Matter (why not...)\n");
    pba->Omega0_cdm = 0.;
  }
  class_test(((flag1 == _SUCCESS_) && (flag2 == _SUCCESS_)),
	     errmsg,
	     "In input file, you can only enter one of Omega_cdm or omega_cdm, choose one");
  if (flag1 == _SUCCESS_)
    pba->Omega0_cdm = param1;
  else
    pba->Omega0_cdm = param2/pba->h/pba->h;

  Omega_tot += pba->Omega0_cdm;

  /* Omega_0_lambda (cosmological constant), Omega0_de (dark energy fluid), Omega0_k (curvature) */
  flag1=parser_read_double(&fc,"Omega_Lambda",&param1,errmsg);
  flag2=parser_read_double(&fc,"Omega_de",&param2,errmsg);
  flag3=parser_read_double(&fc,"Omega_k",&param3,errmsg);
  class_test(((flag1 == _FAILURE_) && (flag2 == _FAILURE_)) || ((flag1 == _FAILURE_) && (flag3 == _FAILURE_)) || ((flag2 == _FAILURE_) && (flag3 == _FAILURE_)),
	     errmsg,
	     "In input file, enter two parameters our of Omega_Lambda, Omega_de, Omega_k (to infer the third one)");
  class_test(((flag1 == _SUCCESS_) && (flag2 == _SUCCESS_) && (flag3 == _SUCCESS_)),
	     errmsg,
	     "In input file, you can enter only two out of Omega_Lambda, Omega_de, Omega_k, the third one is inferred");
  if (flag1 == _FAILURE_) {
    pba->Omega0_lambda= 1. + param3 - param2 - Omega_tot;
    pba->Omega0_de = param2;
    pba->Omega0_k = param3;
  }
  if (flag2 == _FAILURE_) {
    pba->Omega0_lambda= param1;
    pba->Omega0_de = 1. + param3 - param1 - Omega_tot;
    pba->Omega0_k = param3;
  }
  if (flag3 == _FAILURE_) {
    pba->Omega0_lambda= param1;
    pba->Omega0_de = param2;
    pba->Omega0_k = param1 + param2 + Omega_tot - 1.;
  }

/*   printf("h=%f\n",pba->h); */
/*   printf("Omega_g=%e\n",pba->Omega0_g); */
/*   printf("Omega_b=%f\n",pba->Omega0_b); */
/*   printf("Omega_nur=%e\n",pba->Omega0_nur); */
/*   printf("Omega_nur/Omega_g=%e\n",pba->Omega0_nur/pba->Omega0_g); */
/*   printf("Omega_cdm=%f\n",pba->Omega0_cdm); */
/*   printf("Omega_lambda=%f\n",pba->Omega0_lambda); */
/*   printf("Omega_de=%f\n",pba->Omega0_de); */
/*   printf("Omega_k=%f\n",pba->Omega0_k); */

  class_test(pba->Omega0_de != 0.,
	     errmsg,
	     "Dark energy fluid not tested yet");
  
  class_test(pba->Omega0_k != 0.,
	     errmsg,
	     "Open/close case not written yet");

  /* scale factor today (arbitrary) */
  flag1=parser_read_double(&fc,"a_today",&param1,errmsg);
  if (flag1 == _FAILURE_)
    pba->a_today = 1.;
  else 
    pba->a_today = param1;

  /** - assign values to thermodynamics cosmological parameters */

  /* scale factor today (arbitrary) */
  flag1=parser_read_double(&fc,"YHe",&param1,errmsg);
  if (flag1 == _FAILURE_)
    pth->YHe = 0.25;
  else 
    pth->YHe = param1;

  /* reionization parametrization */
  flag1=parser_read_string(&fc,"reio_parametrization",&string1,errmsg);
  if (flag1 == _FAILURE_) {
    pth->reio_parametrization=reio_none;
    printf("Warning: you are computing a model without reionization (why not...)\n");
  }
  else {
    printf("read:%s;\n",string1);
    flag2=_FALSE_;
    if (strcmp(string1,"reio_none") == 0) {
      pth->reio_parametrization=reio_none;
      printf("Warning: you are computing a model without reionization (why not...)\n");
    }
    else {
      if (strcmp(string1,"reio_camb") == 0)
	pth->reio_parametrization=reio_camb;
      else {
	flag2=_TRUE_;
      }
    }
    class_test(flag2==_TRUE_,
	       errmsg,
	       "could not identify reionization_parametrization value, check that it is one of 'reio_none', 'reio_camb', ...");
  }

  /* reionization parameters if reio_parametrization=reio_camb */
  if (pth->reio_parametrization == reio_camb) {
    flag1=parser_read_double(&fc,"z_reio",&param1,errmsg);
    flag2=parser_read_double(&fc,"tau_reio",&param2,errmsg);
    class_test((flag1 == _FAILURE_) && (flag2 == _FAILURE_),
	       errmsg,
	       "Since you have set reionization parameterization to reio_camb, enter one of z_reio or tau_reio");
    class_test(((flag1 == _SUCCESS_) && (flag2 == _SUCCESS_)),
	       errmsg,
	       "In input file, you can only enter one of z_reio or tau_reio, choose one");
    if (flag1 == _SUCCESS_) {
      pth->z_reio=param1;
      pth->reio_z_or_tau=reio_z;
    }
    else {
      pth->tau_reio=param2;
      pth->reio_z_or_tau=reio_tau;
    }
  }
  
  /** - define which perturbations and sources should be computed, and down to which scale */

  flag1=parser_read_string(&fc,"output",&string1,errmsg);
  if (flag1 == _FAILURE_) {
    ppt->has_cl_cmb_temperature = _FALSE_;
    ppt->has_cl_cmb_polarization = _FALSE_;
    ppt->has_cl_cmb_lensing_potential = _FALSE_;
    ppt->has_pk_matter = _FALSE_;
  }
  else {

    if ((strstr(string1,"tCl") == NULL) && (strstr(string1,"TCl") == NULL) && (strstr(string1,"TCL") == NULL))
      ppt->has_cl_cmb_temperature=_FALSE_;
    else
      ppt->has_cl_cmb_temperature=_TRUE_;  

    if ((strstr(string1,"pCl") == NULL) && (strstr(string1,"PCl") == NULL) && (strstr(string1,"PCL") == NULL))
      ppt->has_cl_cmb_polarization=_FALSE_;
    else
      ppt->has_cl_cmb_polarization=_TRUE_;  
    
    if ((strstr(string1,"lCl") == NULL) && (strstr(string1,"LCl") == NULL) && (strstr(string1,"LCL") == NULL))
      ppt->has_cl_cmb_lensing_potential=_FALSE_;
    else
      ppt->has_cl_cmb_lensing_potential=_TRUE_;  

    if ((strstr(string1,"mPk") == NULL) && (strstr(string1,"MPk") == NULL) && (strstr(string1,"MPK") == NULL))
      ppt->has_pk_matter=_FALSE_;
    else
      ppt->has_pk_matter=_TRUE_;  

  }

  flag1=parser_read_string(&fc,"modes",&string1,errmsg);
  if (flag1 == _FAILURE_) {
    ppt->has_scalars=_TRUE_;  
    ppt->has_vectors=_FALSE_;
    ppt->has_tensors=_FALSE_;
  }
  else {
    
    if ((strstr(string1,"s") == NULL) && (strstr(string1,"S") == NULL))
      ppt->has_scalars=_FALSE_;
    else
      ppt->has_scalars=_TRUE_;  

    if ((strstr(string1,"v") == NULL) && (strstr(string1,"V") == NULL))
      ppt->has_vectors=_FALSE_;
    else
      ppt->has_vectors=_TRUE_;  

    if ((strstr(string1,"t") == NULL) && (strstr(string1,"T") == NULL))
      ppt->has_tensors=_FALSE_;
    else
      ppt->has_tensors=_TRUE_;  

    class_test(ppt->has_scalars==_FALSE_ && ppt->has_vectors ==_FALSE_ && ppt->has_tensors ==_FALSE_,
	       errmsg,	       
               "You wrote: modes=%s. Could not identify any of the modes ('s', 'v', 't') in such input",string1);
  }

  if (ppt->has_scalars == _TRUE_) {
    flag1=parser_read_string(&fc,"ic",&string1,errmsg);
    if (flag1 == _FAILURE_) {
      ppt->has_ad=_TRUE_;  
      ppt->has_bi=_FALSE_;
      ppt->has_cdi=_FALSE_;
      ppt->has_nid=_FALSE_;
      ppt->has_niv=_FALSE_;
    }
    else {
      
      if ((strstr(string1,"ad") == NULL) && (strstr(string1,"AD") == NULL))
	ppt->has_ad=_FALSE_;
      else
	ppt->has_ad=_TRUE_; 
      
      if ((strstr(string1,"bi") == NULL) && (strstr(string1,"BI") == NULL))
	ppt->has_bi=_FALSE_;
      else
	ppt->has_bi=_TRUE_; 
      
      if ((strstr(string1,"cdi") == NULL) && (strstr(string1,"CDI") == NULL))
	ppt->has_cdi=_FALSE_;
      else
	ppt->has_cdi=_TRUE_; 
      
      if ((strstr(string1,"nid") == NULL) && (strstr(string1,"NID") == NULL))
	ppt->has_nid=_FALSE_;
      else
	ppt->has_nid=_TRUE_; 
      
      if ((strstr(string1,"niv") == NULL) && (strstr(string1,"NIV") == NULL))
	ppt->has_niv=_FALSE_;
      else
	ppt->has_niv=_TRUE_; 
      
      class_test(ppt->has_ad==_FALSE_ && ppt->has_bi ==_FALSE_ && ppt->has_cdi ==_FALSE_ && ppt->has_nid ==_FALSE_ && ppt->has_niv ==_FALSE_,
		 errmsg,	       
		 "You wrote: ic=%s. Could not identify any of the initial conditions ('ad', 'bi', 'cdi', 'nid', 'niv') in such input",string1);

    }
  }

  /** - define the primordial spectrum */

  flag1=parser_read_string(&fc,"P_k_ini type",&string1,errmsg);
  if (flag1 == _FAILURE_) {
    ppm->primordial_spec_type = analytic_Pk;
  }
  else {
    flag2=_FALSE_;
    if (strcmp(string1,"analytic_Pk") == 0) {
      ppm->primordial_spec_type = analytic_Pk;
    }
    else {
    	flag2=_TRUE_;
    }
    class_test(flag2==_TRUE_,
	       errmsg,
	       "could not identify primordial spectrum type, check that it is one of 'analytic_pk', ...");
  }

  if (ppm->primordial_spec_type == analytic_Pk) {

    flag1=parser_read_double(&fc,"A_s_ad",&param1,errmsg);
    if (flag1==_FAILURE_)
      ppm->A_s_ad = 1.;
    else
      ppm->A_s_ad = param1;

    flag1=parser_read_double(&fc,"n_s_ad",&param1,errmsg);
    if (flag1==_FAILURE_)
      ppm->n_s_ad = 1.;
    else
      ppm->n_s_ad = param1;

    flag1=parser_read_double(&fc,"alpha_s_ad",&param1,errmsg);
    if (flag1==_FAILURE_)
      ppm->alpha_s_ad = 0.;
    else
      ppm->alpha_s_ad = param1;

    flag1=parser_read_double(&fc,"k_pivot",&param1,errmsg);
    if (flag1==_FAILURE_)
      ppm->k_pivot = 0.05;
    else
      ppm->k_pivot = param1;

  }

  /** - define up to which value of z P(k) should be computed */

  psp->z_max_pk = 0.;

  /** - prefix for name of output files (blank if no output files needed);
   number of redshift for output P(k) (including 0 and z_max_pk if non-zero) */

  pop->cls_ad = "output/cls.dat";
  pop->pk = "output/pk.dat";
  pop->z_pk = 0.;

  /** - amount of information sent to standard output (none if all set to zero) */

  pba->background_verbose = 1;
  pth->thermodynamics_verbose = 1;
  ppt->perturbations_verbose = 2;
  pbs->bessels_verbose = 2;
  ptr->transfer_verbose = 2;
  ppm->primordial_verbose = 1;
  psp->spectra_verbose = 2;
  pop->output_verbose = 1;

  return _SUCCESS_;

}

int input_check_arguments_of_main(
				  int argc, 
				  char **argv, 
				  char * input,
				  ErrorMsg errmsg) {

  int i;
  char extension[5];

  class_test(argc == 1,
	     errmsg,
	     "No input file xxx.ini, run with e.g. \n >./class params.ini");

  input[0]='\0';
  for (i=1; i<argc; i++) {
    strncpy(extension,(argv[i]+strlen(argv[i])-4),4);
    extension[4]='\0';
    if (strcmp(extension,".ini") == 0) {
      class_test(input[0] != '\0',
		 errmsg,
		 "You have passed more than one input file xxx.ini. Choose one.");
      strcpy(input,argv[i]);
    }
  }
  class_test(input[0] == '\0',
	     errmsg,
	     "No input file xxx.ini, run with e.g. \n >./class params.ini");

  return _SUCCESS_;

}
