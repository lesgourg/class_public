#include "non_cold_dark_matter.h"

std::shared_ptr<NonColdDarkMatter> NonColdDarkMatter::Create(FileContent* pfc, const NcdmSettings& ncdm_settings) {
  int flag1;
  int N_ncdm;
  ErrorMsg error_message;
  parser_read_int(pfc, "N_ncdm", &N_ncdm, &flag1, error_message);
  if ((flag1 == _TRUE_) && (N_ncdm > 0)) {
    try {
      return std::shared_ptr<NonColdDarkMatter>(new NonColdDarkMatter(pfc, ncdm_settings));
    }
    catch (...) {
      return nullptr;
    }
  }
  else {
    return nullptr;
  }
}

NonColdDarkMatter::NonColdDarkMatter(FileContent* pfc, const NcdmSettings& ncdm_settings) {
  rho_nu_rel_ = 56.0/45.0*pow(_PI_, 6)*pow(4.0/11.0, 4.0/3.0)*_G_/pow(_h_P_, 3)/pow(_c_, 7)*pow(_Mpc_over_m_, 2)*pow(ncdm_settings.T_cmb*_k_B_, 4);
  int return_value = background_ncdm_init(pfc, ncdm_settings);
  ThrowInvalidArgumentIf(return_value == _FAILURE_, error_message_);
}

NonColdDarkMatter::~NonColdDarkMatter() {
  SafeFree(M_ncdm_);
  SafeFree(Omega0_ncdm_);
  SafeFree(omega0_ncdm_);
  SafeFree(m_ncdm_in_eV_);
  SafeFree(deg_ncdm_);
  SafeFree(T_ncdm_);
  SafeFree(ksi_ncdm_);
  SafeFree(ncdm_psd_parameters_);
  SafeFree(got_files_);
  SafeFree(ncdm_psd_files_);

  SafeFree(ncdm_quadrature_strategy_);
  SafeFree(ncdm_input_q_size_);
  SafeFree(q_size_ncdm_bg_);
  SafeFree(q_size_ncdm_);
  SafeFree(ncdm_qmax_);
  SafeFree(factor_ncdm_);

  SafeFree(q_ncdm_);
  SafeFree(w_ncdm_);
  SafeFree(dlnf0_dlnq_ncdm_);
  SafeFree(q_ncdm_bg_);
  SafeFree(w_ncdm_bg_);
}


/**
 * This is the routine where the distribution function f0(q) of each
 * ncdm species is specified (it is the only place to modify if you
 * need a partlar f0(q))
 *
 * @param pbadist Input:  structure containing all parameters defining f0(q)
 * @param q       Input:  momentum
 * @param f0      Output: phase-space distribution
 */

int NonColdDarkMatter::background_ncdm_distribution(void* pbadist, double q, double* f0) {
  background_parameters_for_distributions* pbadist_local = static_cast<background_parameters_for_distributions*>(pbadist);
  int n_ncdm = pbadist_local->n_ncdm;   /* extract index of ncdm species under consideration */
  double ksi = pbadist_local->ncdm->ksi_ncdm_[n_ncdm];      /* extract chemical potential */
  double qlast,dqlast,f0last,df0last;
  double* param = pbadist_local->ncdm->ncdm_psd_parameters_;
  int lastidx;
  /* Variables corresponding to entries in param: */
  //double square_s12,square_s23,square_s13;
  //double mixing_matrix[3][3];
  //int i;

  /** - shall we interpolate in file, or shall we use analytical formula below? */

  /** - a) deal first with the case of interpolating in files */
  if (pbadist_local->ncdm->got_files_[n_ncdm] == _TRUE_) {

    lastidx = pbadist_local->tablesize-1;
    if(q<pbadist_local->q[0]){
      //Handle q->0 case:
      *f0 = pbadist_local->f0[0];
    }
    else if(q>pbadist_local->q[lastidx]){
      //Handle q>qmax case (ensure continuous and derivable function with Boltzmann tail):
      qlast=pbadist_local->q[lastidx];
      f0last=pbadist_local->f0[lastidx];
      dqlast = qlast - pbadist_local->q[lastidx - 1];
      df0last = f0last - pbadist_local->f0[lastidx - 1];

      *f0 = f0last*exp(-(qlast - q)*df0last/f0last/dqlast);
    }
    else{
      //Do interpolation:
      class_call(array_interpolate_spline(
                                          pbadist_local->q,
                                          pbadist_local->tablesize,
                                          pbadist_local->f0,
                                          pbadist_local->d2f0,
                                          1,
                                          q,
                                          &pbadist_local->last_index,
                                          f0,
                                          1,
                                          pbadist_local->ncdm->error_message_),
                 pbadist_local->ncdm->error_message_, pbadist_local->ncdm->error_message_);
    }
  }

  /** - b) deal now with case of reading analytical function */
  else{
    /**
       Next enter your analytic expression(s) for the p.s.d.'s. If
       you need different p.s.d.'s for different species, put each
       p.s.d inside a condition, like for instance: if (n_ncdm==2)
       {*f0=...}.  Remember that n_ncdm = 0 refers to the first
       species.
    */

    /**************************************************/
    /*    FERMI-DIRAC INCLUDING CHEMICAL POTENTIALS   */
    /**************************************************/

    *f0 = 1.0/pow(2*_PI_, 3)*(1./(exp(q - ksi) + 1.) + 1./(exp(q + ksi) + 1.));

    /**************************************************/

    /** This form is only appropriate for approximate studies, since in
        reality the chemical potentials are associated with flavor
        eigenstates, not mass eigenstates. It is easy to take this into
        account by introducing the mixing angles. In the later part
        (not read by the code) we illustrate how to do this. */

    if (_FALSE_) {

      /* We must use the list of extra parameters read in input, stored in the
         ncdm_psd_parameter list, extracted above from the structure
         and now called param[..] */

      /* check that this list has been read */
      class_test(param == NULL,
                 pbadist_local->ncdm->error_message_,
                 "Analytic expression wants to use 'ncdm_psd_parameters', but they have not been entered!");

      /* extract values from the list (in this example, mixing angles) */
      double square_s12 = param[0];
      double square_s23 = param[1];
      double square_s13 = param[2];

      /* infer mixing matrix */
      double mixing_matrix[3][3];
      int i;

      mixing_matrix[0][0] = pow(fabs(sqrt((1 - square_s12)*(1 - square_s13))), 2);
      mixing_matrix[0][1] = pow(fabs(sqrt(square_s12*(1 - square_s13))), 2);
      mixing_matrix[0][2] = fabs(square_s13);
      mixing_matrix[1][0] = pow(fabs(sqrt((1 - square_s12)*square_s13*square_s23) + sqrt(square_s12*(1 - square_s23))), 2);
      mixing_matrix[1][1] = pow(fabs(sqrt(square_s12*square_s23*square_s13) - sqrt((1 - square_s12)*(1 - square_s23))), 2);
      mixing_matrix[1][2] = pow(fabs(sqrt(square_s23*(1 - square_s13))), 2);
      mixing_matrix[2][0] = pow(fabs(sqrt(square_s12*square_s23) - sqrt((1 - square_s12)*square_s13*(1 - square_s23))), 2);
      mixing_matrix[2][1] = pow(sqrt((1 - square_s12)*square_s23) + sqrt(square_s12*square_s13*(1 - square_s23)), 2);
      mixing_matrix[2][2] = pow(fabs(sqrt((1 - square_s13)*(1 - square_s23))), 2);

      /* loop over flavor eigenstates and compute psd of mass eigenstates */
      *f0=0.0;
      for(i = 0; i < 3; i++){
        *f0 += mixing_matrix[i][n_ncdm]*1.0/pow(2*_PI_, 3)*(1./(exp(q - pbadist_local->ncdm->ksi_ncdm_[i]) + 1.) + 1./(exp(q + pbadist_local->ncdm->ksi_ncdm_[i]) + 1.));
      }
    } /* end of region not used, but shown as an example */
  }

  return _SUCCESS_;
}

/**
 * This function is only used for the purpose of finding optimal
 * quadrature weights. The logic is: if we can accurately convolve
 * f0(q) with this function, then we can convolve it accurately with
 * any other relevant function.
 *
 * @param pbadist Input:  structure containing all background parameters
 * @param q       Input:  momentum
 * @param test    Output: value of the test function test(q)
 */

int NonColdDarkMatter::background_ncdm_test_function(
                                  void* pbadist,
                                  double q,
                                  double* test
                                  ) {

  double c = 2.0/(3.0*_zeta3_);
  double d = 120.0/(7.0*pow(_PI_, 4));
  double e = 2.0/(45.0*_zeta5_);

  /** Using a + bq creates problems for otherwise acceptable distributions
      which diverges as \f$ 1/r \f$ or \f$ 1/r^2 \f$ for \f$ r\to 0 \f$*/
  *test = pow(2.0*_PI_, 3)/6.0*(c*q*q - d*q*q*q - e*q*q*q*q);

  return _SUCCESS_;
}

/**
 * This function finds optimal quadrature weights for each ncdm
 * species
 *
 * @param ppr Input: precision structure
 * @param pba Input/Output: background structure
 */

int NonColdDarkMatter::background_ncdm_init(FileContent* pfc, const NcdmSettings& ncdm_settings) {
  int n;
  int int1;
  int flag1;
  int entries_read;
  char* errmsg = error_message_;
  class_read_int("N_ncdm", N_ncdm_);

  /* Quadrature modes, 0 is qm_auto. */
  class_read_list_of_integers_or_default("Quadrature strategy", ncdm_quadrature_strategy_, 0, N_ncdm_);
  /* Number of momentum bins */
  class_read_list_of_integers_or_default("Number of momentum bins", ncdm_input_q_size_, -1, N_ncdm_);

  /* qmax, if relevant */
  class_read_list_of_doubles_or_default("Maximum q", ncdm_qmax_, 15, N_ncdm_);

  /* Read temperatures: */
  class_read_list_of_doubles_or_default("T_ncdm", T_ncdm_, T_ncdm_default_, N_ncdm_);

  /* Read chemical potentials: */
  class_read_list_of_doubles_or_default("ksi_ncdm", ksi_ncdm_, ksi_ncdm_default, N_ncdm_);

  /* Read degeneracy of each ncdm species: */
  class_read_list_of_doubles_or_default("deg_ncdm", deg_ncdm_, deg_ncdm_default_, N_ncdm_);

  /* Read mass of each ncdm species: */
  class_read_list_of_doubles_or_default("m_ncdm", m_ncdm_in_eV_, 0.0, N_ncdm_);

  /* Read Omega of each ncdm species: */
  class_read_list_of_doubles_or_default("Omega_ncdm", Omega0_ncdm_, 0.0, N_ncdm_);

  /* Read omega of each ncdm species: */
  class_read_list_of_doubles_or_default("omega_ncdm", omega0_ncdm_, 0.0, N_ncdm_);

  /* Check for duplicate Omega/omega entries, missing mass definition and
     update Omega0_ncdm_:*/
  for (int n = 0; n < N_ncdm_; n++){
    if (omega0_ncdm_[n] != 0.0){
      class_test(Omega0_ncdm_[n] != 0,errmsg,
                 "Nonzero values for both Omega and omega for ncdm species %d are specified!", n);
      Omega0_ncdm_[n] = omega0_ncdm_[n]/ncdm_settings.h/ncdm_settings.h;
    }
    else {
      omega0_ncdm_[n] = Omega0_ncdm_[n]*ncdm_settings.h*ncdm_settings.h;
    }
    if ((Omega0_ncdm_[n] == 0.0) && (m_ncdm_in_eV_[n] == 0.0)) {
      /* this is the right place for passing the default value of
         the mass (all parameters must have a default value; most of
         them are defined in input_default_params{}, but the ncdm mass
         is a bit special and there is no better place for setting its
         default value). We put an arbitrary value m << 10^-3 eV,
         i.e. the ultra-relativistic limit.*/
      m_ncdm_in_eV_[n] = 1.e-5;
    }
  }

  /* Check if filenames for interpolation tables are given: */
  class_read_list_of_integers_or_default("use_ncdm_psd_files", got_files_, _FALSE_, N_ncdm_);

  if (flag1 == _TRUE_){
    int fileentries = 0;
    for (n = 0; n < N_ncdm_; n++){
      if (got_files_[n] == _TRUE_) fileentries++;
    }

    if (fileentries > 0) {
      /* Okay, read filenames.. */
      int flag2;
      class_call(parser_read_list_of_strings(pfc, "ncdm_psd_filenames", &entries_read, &ncdm_psd_files_, &flag2, errmsg),
                 errmsg,
                 errmsg);
      class_test(flag2 == _FALSE_, errmsg, "Input use_ncdm_files is found, but no filenames found!");
      class_test(entries_read != fileentries,
                 errmsg,
                 "Number of filenames found, %d, does not match number of _TRUE_ values in use_ncdm_files, %d",
                 entries_read,
                 fileentries);
    }
  }
  /* Read (optional) p.s.d.-parameters:*/
  parser_read_list_of_doubles(pfc,
                              "ncdm_psd_parameters",
                              &entries_read,
                              &ncdm_psd_parameters_,
                              &flag1,
                              errmsg);


  int index_q,tolexp;
  double f0m2,f0m1,f0,f0p1,f0p2,q,df0dq;

  /* Allocate pointer arrays: */
  class_alloc(q_ncdm_, sizeof(double*)*N_ncdm_, error_message_);
  class_alloc(w_ncdm_, sizeof(double*)*N_ncdm_, error_message_);
  class_alloc(q_ncdm_bg_, sizeof(double*)*N_ncdm_, error_message_);
  class_alloc(w_ncdm_bg_, sizeof(double*)*N_ncdm_, error_message_);
  class_alloc(dlnf0_dlnq_ncdm_, sizeof(double*)*N_ncdm_, error_message_);

  for (int n = 0; n < N_ncdm_; ++n) {
    q_ncdm_[n] = nullptr;
    w_ncdm_[n] = nullptr;
    q_ncdm_bg_[n] = nullptr;
    w_ncdm_bg_[n] = nullptr;
    dlnf0_dlnq_ncdm_[n] = nullptr;
  }

  /* Allocate pointers: */
  class_alloc(q_size_ncdm_, sizeof(int)*N_ncdm_, error_message_);
  class_alloc(q_size_ncdm_bg_, sizeof(int)*N_ncdm_, error_message_);
  class_alloc(factor_ncdm_, sizeof(double)*N_ncdm_, error_message_);
  class_alloc(M_ncdm_, sizeof(double)*N_ncdm_, error_message_);

  int filenum = 0;
  for(int k = 0; k < N_ncdm_; k++){
    background_parameters_for_distributions pbadist(this, k);
    pbadist.n_ncdm = k;
    pbadist.q = NULL;
    pbadist.tablesize = 0;
    /*Do we need to read in a file to interpolate the distribution function? */
    if ((got_files_ != NULL) && (got_files_[k] == _TRUE_)) {
      FILE* psdfile = fopen(ncdm_psd_files_ + filenum*_ARGUMENT_LENGTH_MAX_, "r");
      class_test(psdfile == NULL, error_message_, "Could not open file %s!", ncdm_psd_files_ + filenum*_ARGUMENT_LENGTH_MAX_);
      // Find size of table:
      double tmp1;
      double tmp2;
      int status = 0;
      int row;
      for (row = 0; status == 2; row++) {
        status = fscanf(psdfile, "%lf %lf", &tmp1, &tmp2);
      }
      rewind(psdfile);
      pbadist.tablesize = row - 1;

      /*Allocate room for interpolation table: */
      class_alloc(pbadist.q, sizeof(double)*pbadist.tablesize, error_message_);
      class_alloc(pbadist.f0, sizeof(double)*pbadist.tablesize, error_message_);
      class_alloc(pbadist.d2f0, sizeof(double)*pbadist.tablesize, error_message_);
      for (row = 0; row < pbadist.tablesize; row++){
        status = fscanf(psdfile, "%lf %lf", &pbadist.q[row], &pbadist.f0[row]);
      }
      fclose(psdfile);
      /* Call spline interpolation: */
      class_call(array_spline_table_lines(pbadist.q,
                                          pbadist.tablesize,
                                          pbadist.f0,
                                          1,
                                          pbadist.d2f0,
                                          _SPLINE_EST_DERIV_,
                                          error_message_),
                 error_message_,
                 error_message_);
      filenum++;
    }

    /* Handle perturbation qsampling: */
    if (ncdm_quadrature_strategy_[k] == qm_auto) {
      /** Automatic q-sampling for this species */
      class_alloc(q_ncdm_[k], _QUADRATURE_MAX_*sizeof(double), error_message_);
      class_alloc(w_ncdm_[k], _QUADRATURE_MAX_*sizeof(double), error_message_);

      class_call(get_qsampling(q_ncdm_[k],
                               w_ncdm_[k],
                               &(q_size_ncdm_[k]),
                               _QUADRATURE_MAX_,
                               ncdm_settings.tol_ncdm,
                               pbadist.q,
                               pbadist.tablesize,
                               background_ncdm_test_function,
                               background_ncdm_distribution,
                               &pbadist,
                               error_message_),
                 error_message_,
                 error_message_);
      q_ncdm_[k] = (double*)realloc(q_ncdm_[k], q_size_ncdm_[k]*sizeof(double));
      w_ncdm_[k] = (double*)realloc(w_ncdm_[k], q_size_ncdm_[k]*sizeof(double));

      /* Handle background q_sampling: */
      class_alloc(q_ncdm_bg_[k], _QUADRATURE_MAX_BG_*sizeof(double), error_message_);
      class_alloc(w_ncdm_bg_[k], _QUADRATURE_MAX_BG_*sizeof(double), error_message_);

      class_call(get_qsampling(q_ncdm_bg_[k],
                               w_ncdm_bg_[k],
                               &(q_size_ncdm_bg_[k]),
                               _QUADRATURE_MAX_BG_,
                               ncdm_settings.tol_ncdm_bg,
                               pbadist.q,
                               pbadist.tablesize,
                               background_ncdm_test_function,
                               background_ncdm_distribution,
                               &pbadist,
                               error_message_),
                 error_message_,
                 error_message_);


      q_ncdm_bg_[k] = (double*)realloc(q_ncdm_bg_[k], q_size_ncdm_bg_[k]*sizeof(double));
      w_ncdm_bg_[k] = (double*)realloc(w_ncdm_bg_[k], q_size_ncdm_bg_[k]*sizeof(double));
    }
    else{
      /** Manual q-sampling for this species. Same sampling used for both perturbation and background sampling, since this will usually be a high precision setting anyway */
      q_size_ncdm_bg_[k] = ncdm_input_q_size_[k];
      q_size_ncdm_[k] = ncdm_input_q_size_[k];
      class_alloc(q_ncdm_bg_[k], q_size_ncdm_bg_[k]*sizeof(double), error_message_);
      class_alloc(w_ncdm_bg_[k], q_size_ncdm_bg_[k]*sizeof(double), error_message_);
      class_alloc(q_ncdm_[k], q_size_ncdm_[k]*sizeof(double), error_message_);
      class_alloc(w_ncdm_[k], q_size_ncdm_[k]*sizeof(double), error_message_);
      class_call(get_qsampling_manual(q_ncdm_[k],
                                      w_ncdm_[k],
                                      q_size_ncdm_[k],
                                      ncdm_qmax_[k],
                                      (enum ncdm_quadrature_method)ncdm_quadrature_strategy_[k],
                                      pbadist.q,
                                      pbadist.tablesize,
                                      background_ncdm_distribution,
                                      &pbadist,
                                      error_message_),
                 error_message_,
                 error_message_);
      for (index_q = 0; index_q < q_size_ncdm_[k]; index_q++) {
        q_ncdm_bg_[k][index_q] = q_ncdm_[k][index_q];
        w_ncdm_bg_[k][index_q] = w_ncdm_[k][index_q];
      }

    }

    class_alloc(dlnf0_dlnq_ncdm_[k],
                q_size_ncdm_[k]*sizeof(double),
                error_message_);


    for (index_q = 0; index_q < q_size_ncdm_[k]; index_q++) {
      q = q_ncdm_[k][index_q];
      class_call(background_ncdm_distribution(&pbadist, q, &f0),
                 error_message_, error_message_);

      //Loop to find appropriate dq:
      double dq = 1.;
      for(tolexp = _PSD_DERIVATIVE_EXP_MIN_; tolexp<_PSD_DERIVATIVE_EXP_MAX_; tolexp++){

        if (index_q == 0){
          dq = MIN((0.5 - DBL_EPSILON)*q, 2*exp(tolexp)*(q_ncdm_[k][index_q + 1] - q));
        }
        else if (index_q == q_size_ncdm_[k] - 1){
          dq = exp(tolexp)*2.0*(q_ncdm_[k][index_q] - q_ncdm_[k][index_q - 1]);
        }
        else{
          dq = exp(tolexp)*(q_ncdm_[k][index_q + 1] - q_ncdm_[k][index_q - 1]);
        }

        class_call(background_ncdm_distribution(&pbadist, q - 2*dq, &f0m2),
                   error_message_, error_message_);
        class_call(background_ncdm_distribution(&pbadist, q + 2*dq, &f0p2),
                   error_message_, error_message_);

        if (fabs((f0p2 - f0m2)/f0) > sqrt(DBL_EPSILON)) break;
      }

      class_call(background_ncdm_distribution(&pbadist, q - dq, &f0m1),
                 error_message_, error_message_);
      class_call(background_ncdm_distribution(&pbadist, q + dq, &f0p1),
                 error_message_, error_message_);
      //5 point estimate of the derivative:
      df0dq = (+f0m2 - 8*f0m1 + 8*f0p1 - f0p2)/12.0/dq;
      //Avoid underflow in extreme tail:
      if (fabs(f0) == 0.)
        dlnf0_dlnq_ncdm_[k][index_q] = -q; /* valid for whatever f0 with exponential tail in exp(-q) */
      else
        dlnf0_dlnq_ncdm_[k][index_q] = q/f0*df0dq;
    }

    factor_ncdm_[k] = deg_ncdm_[k]*4*_PI_*pow(ncdm_settings.T_cmb*T_ncdm_[k]*_k_B_, 4)*8*_PI_*_G_
      /3./pow(_h_P_/2./_PI_, 3)/pow(_c_, 7)*_Mpc_over_m_*_Mpc_over_m_;

    /* If allocated, deallocate interpolation table:  */
  }

  /* We must calculate M from omega or vice versa if one of them is missing.
     If both are present, we must update the degeneracy parameter to
     reflect the implicit normalization of the distribution function.*/
  double H0 = ncdm_settings.h*1.e5/_c_;
  for (n=0; n < N_ncdm_; n++){
    if (m_ncdm_in_eV_[n] != 0.0){
      /* Case of only mass or mass and Omega/omega: */
      M_ncdm_[n] = m_ncdm_in_eV_[n]/_k_B_*_eV_/T_ncdm_[n]/ncdm_settings.T_cmb;
      double rho_ncdm;
      class_call(background_ncdm_momenta(n, 0., NULL, &rho_ncdm, NULL, NULL, NULL),
                 error_message_,
                 errmsg);
      if (Omega0_ncdm_[n] == 0.0){
        Omega0_ncdm_[n] = rho_ncdm/H0/H0;
        omega0_ncdm_[n] = Omega0_ncdm_[n]*ncdm_settings.h*ncdm_settings.h;
      }
      else{
        double fnu_factor = (H0*H0*Omega0_ncdm_[n]/rho_ncdm);
        factor_ncdm_[n] *= fnu_factor;
        // dlnf0dlnq is already computed, but it is independent of any normalization of f0.
        deg_ncdm_[n] *= fnu_factor;
      }
    }
    else {
      /* Case of only Omega/omega: */
      M_ncdm_[n] = background_ncdm_M_from_Omega(n, H0, Omega0_ncdm_[n], ncdm_settings.tol_M_ncdm);
      m_ncdm_in_eV_[n] = _k_B_/_eV_*T_ncdm_[n]*M_ncdm_[n]*ncdm_settings.T_cmb;
    }
  }
  return _SUCCESS_;
}

/**
 * For a given ncdm species: given the quadrature weights, the mass
 * and the redshift, find background quantities by a quick weighted
 * sum over.  Input parameters passed as NULL pointers are not
 * evaluated for speed-up
 *
 * @param qvec     Input: sampled momenta
 * @param wvec     Input: quadrature weights
 * @param qsize    Input: number of momenta/weights
 * @param M        Input: mass
 * @param factor   Input: normalization factor for the p.s.d.
 * @param z        Input: redshift
 * @param n        Output: number density
 * @param rho      Output: energy density
 * @param p        Output: pressure
 * @param drho_dM  Output: derivative used in next function
 * @param pseudo_p Output: pseudo-pressure used in perturbation module for fluid approx
 *
 */

int NonColdDarkMatter::background_ncdm_momenta(int n_ncdm, double z, double* n, double* rho, double* p, double* drho_dM, double* pseudo_p) const {
  return background_ncdm_momenta_mass(n_ncdm, M_ncdm_[n_ncdm], z, n, rho, p, drho_dM, pseudo_p);
}

int NonColdDarkMatter::background_ncdm_momenta_mass(int n_ncdm, double M, double z, double* n, double* rho, double* p, double* drho_dM, double* pseudo_p) const {

  /** Summary: */

  /** - rescale normalization at given redshift */
  double factor2 = factor_ncdm_[n_ncdm]*pow(1 + z, 4);
  double* qvec = q_ncdm_bg_[n_ncdm];
  double* wvec = w_ncdm_bg_[n_ncdm];
  double qsize = q_size_ncdm_bg_[n_ncdm];
  /** - initialize quantities */
  if (n != NULL) *n = 0.;
  if (rho != NULL) *rho = 0.;
  if (p != NULL) *p = 0.;
  if (drho_dM != NULL) *drho_dM = 0.;
  if (pseudo_p != NULL) *pseudo_p = 0.;

  /** - loop over momenta */
  for (int index_q = 0; index_q < qsize; index_q++) {

    /* squared momentum */
    double q2 = qvec[index_q]*qvec[index_q];

    /* energy */
    double epsilon = sqrt(q2 + M*M/(1. + z)/(1. + z));

    /* integrand of the various quantities */
    if (n != NULL) *n += q2*wvec[index_q];
    if (rho != NULL) *rho += q2*epsilon*wvec[index_q];
    if (p != NULL) *p += q2*q2/3./epsilon*wvec[index_q];
    if (drho_dM != NULL) *drho_dM += q2*M/(1. + z)/(1. + z)/epsilon*wvec[index_q];
    if (pseudo_p != NULL) *pseudo_p += pow(q2/epsilon, 3)/3.0*wvec[index_q];
  }

  /** - adjust normalization */
  if (n != NULL) *n *= factor2/(1. + z);
  if (rho != NULL) *rho *= factor2;
  if (p != NULL) *p *= factor2;
  if (drho_dM != NULL) *drho_dM *= factor2;
  if (pseudo_p != NULL) *pseudo_p *= factor2;

  return _SUCCESS_;
}

/**
 * When the user passed the density fraction Omega_ncdm or
 * omega_ncdm in input but not the mass, infer the mass with Newton iteration method.
 *
 * @param ppr    Input: precision structure
 * @param pba    Input/Output: background structure
 * @param n_ncdm Input: index of ncdm species
 */

double NonColdDarkMatter::background_ncdm_M_from_Omega(int n_ncdm, double H0, double Omega0, double tol_M_ncdm) {
  int maxiter = 50;
  double rho0 = H0*H0*Omega0; /*Remember that rho is defined such that H^2=sum(rho_i) */
  double M = 0;
  double rho;
  double n;
  background_ncdm_momenta_mass(n_ncdm, M, 0., &n, &rho, NULL, NULL, NULL);

  /* Is the value of Omega less than a massless species?*/
  class_test(rho0 < rho, error_message_,
             "The value of Omega for the %dth species, %g, is less than for a massless species! It should be at least %g. Check your input.",
             n_ncdm, Omega0, Omega0*rho/rho0);

  /* In the strict NR limit we have rho = n*(M) today, giving a zeroth order guess: */
  M = rho0/n; /* This is our guess for M. */
  for (int iter = 1; iter <= maxiter; iter++){

    /* Newton iteration. First get relevant quantities at M: */
    double drhodM;
    background_ncdm_momenta_mass(n_ncdm, M, 0., NULL, &rho, NULL, &drhodM, NULL);

    double deltaM = (rho0 - rho)/drhodM; /* By definition of the derivative */
    if ((M + deltaM) < 0.0) {
      deltaM = -M/2.0; /* Avoid overshooting to negative M value. */
    }
    M += deltaM; /* Update value of M.. */
    if (fabs(deltaM/M) < tol_M_ncdm){
      /* Accuracy reached.. */
      return M;
    }
  }
  ThrowInvalidArgument("Newton iteration could not converge on a mass for some reason.");
  return 0.;
}

double NonColdDarkMatter::GetOmega0() const {
  double res = 0;
  for (int i = 0; i < N_ncdm_; ++i) {
    res += Omega0_ncdm_[i];
  }
  return res;
}

void NonColdDarkMatter::SafeFree(double* pointer) {
  if (pointer) {
    free(pointer);
  }
}

void NonColdDarkMatter::SafeFree(int* pointer) {
  if (pointer) {
    free(pointer);
  }
}

void NonColdDarkMatter::SafeFree(char* pointer) {
  if (pointer) {
    free(pointer);
  }
}

void NonColdDarkMatter::SafeFree(double** double_pointer) {
  if (!double_pointer) {
    return;
  }
  for (int n = 0; n < N_ncdm_; ++n) {
    SafeFree(double_pointer[n]);
  }
  free(double_pointer);
}

void NonColdDarkMatter::PrintNeffInfo() const {
  int filenum = 0;
  for (int n_ncdm = 0; n_ncdm < N_ncdm_; ++n_ncdm) {
    /* inform if p-s-d read in files */
    if (got_files_[n_ncdm] == _TRUE_) {
      printf(" -> ncdm species i=%d read from file %s\n", n_ncdm + 1, ncdm_psd_files_ + filenum*_ARGUMENT_LENGTH_MAX_);
      filenum++;
    }
    double rho_ncdm_rel;
    background_ncdm_momenta_mass(n_ncdm, 0., 0., NULL, &rho_ncdm_rel, NULL, NULL, NULL);

    /* inform user of the contribution of each species to
       radiation density (in relativistic limit): should be
       between 1.01 and 1.02 for each active neutrino species;
       evaluated as rho_ncdm/rho_nu_rel where rho_nu_rel is the
       density of one neutrino in the instantaneous decoupling
       limit, i.e. assuming T_nu=(4/11)^1/3 T_gamma (this comes
       from the definition of N_eff) */

    printf(" -> ncdm species i=%d sampled with %d (resp. %d) points for purpose of background (resp. perturbation) integration. In the relativistic limit it gives Delta N_eff = %g\n",
           n_ncdm + 1,
           q_size_ncdm_bg_[n_ncdm],
           q_size_ncdm_[n_ncdm],
           rho_ncdm_rel/rho_nu_rel_);
  }
}

void NonColdDarkMatter::PrintMassInfo() const {
  for (int n_ncdm = 0; n_ncdm < N_ncdm_; n_ncdm++) {
    printf(" -> non-cold dark matter species with i=%d has m_i = %e eV (so m_i / omega_i =%e eV)\n",
           n_ncdm + 1,
           m_ncdm_in_eV_[n_ncdm],
           m_ncdm_in_eV_[n_ncdm]*deg_ncdm_[n_ncdm]/omega0_ncdm_[n_ncdm]);
  }
}

void NonColdDarkMatter::PrintOmegaInfo() const {
  for (int index_ncdm = 0; index_ncdm < N_ncdm_; ++index_ncdm){
    printf("-> %-26s%-4d Omega = %-15g , omega = %-15g\n",
           "Neutrino Species Nr.",
           index_ncdm + 1,
           Omega0_ncdm_[index_ncdm],
           omega0_ncdm_[index_ncdm]);
  }
}

double NonColdDarkMatter::GetNeff(double z) const {
  double Neff = 0;
  for (int n_ncdm = 0; n_ncdm < N_ncdm_; ++n_ncdm) {
    /* call this function to get rho_ncdm */
    double rho_ncdm_rel;
    background_ncdm_momenta_mass(n_ncdm, 0., z, NULL, &rho_ncdm_rel, NULL, NULL, NULL);
    Neff += rho_ncdm_rel/rho_nu_rel_;
  }
  return Neff;
}

double NonColdDarkMatter::GetMassInElectronvolt(int n_ncdm) const {
  if ((n_ncdm < 0) || (n_ncdm >= N_ncdm_)) {
    return 0.;
  }
  else {
    return m_ncdm_in_eV_[n_ncdm];
  }
}
