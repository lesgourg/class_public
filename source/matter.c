/** @file New way of calculating angular power spectra
 *
 * Nils Schöneberg, 16.10.2017
 *
 */

#include "matter.h"
#include "hypergeom.h"
#include <time.h>

#define CHUNK_SIZE 4
#define WARN_LARGE_NU_IMAG 40.0
#define MATTER_REWRITE_PRINTING _FALSE_

#define MATTER_VERBOSITY_TIMING 1
#define MATTER_VERBOSITY_INDICES 5
#define MATTER_VERBOSITY_FUNCTIONS 3
#define MATTER_VERBOSITY_PARAMETERS 2
#define MATTER_VERBOSITY_RANGES 5
#define MATTER_VERBOSITY_BESSEL 5
#define MATTER_VERBOSITY_CLCALCULATION 4
#define MATTER_VERBOSITY_CLCALCULATION_PARTIAL 5
#define MATTER_VERBOSITY_CLRESULTS 2

#define matter_is_index(index_from,index_to,condition)              \
((condition) && ((index_from)==(index_to)))

#define matter_is_integrated(index_radtp)                           \
((pma->has_lensing_terms && index_radtp == pma->radtp_nclens)       \
|| (pma->has_gravitational_terms &&                                 \
(index_radtp == pma->radtp_g4 || index_radtp == pma->radtp_g5) )    \
|| (pma->has_cl_shear && index_radtp == pma->radtp_shlens)          \
)

#define matter_get_t(index_t)                                       \
if(integrate_logarithmically && !pma->uses_limber_approximation){   \
  double y = pma->t_spline_sampling[(index_t)]*(y_max-y_min)+y_min; \
  t = 1.0-exp(-y);                                                  \
}else if(!pma->uses_limber_approximation){                          \
  t = pma->t_spline_sampling[(index_t)]*(t_max-t_min)+t_min;        \
}                                                                   \
else{                                                               \
  t = 1.0;                                                          \
}
#define matter_get_t_orig(index_t)                                  \
if(integrate_logarithmically){                                      \
  double y = pma->t_sampling[index_t]*(y_max-y_min)+y_min;          \
  t = 1.0-exp(-y);                                                  \
}else{                                                              \
  t = pma->t_sampling[index_t]*(t_max-t_min)+t_min;                 \
}
#define matter_get_tij_limits(index_wd1,index_wd2)                    \
t_min = MAX((pma->tau0-pma->tw_max[index_wd1])/(pma->tau0-pma->tw_min[index_wd2]),0.0+pma->bi_maximal_t_offset);\
t_max = MIN((pma->tau0-pma->tw_min[index_wd1])/(pma->tau0-pma->tw_max[index_wd2]),1.0-pma->bi_maximal_t_offset);
#define matter_get_t_limits(index_wd1,index_wd2)                    \
matter_get_tij_limits(index_wd2,index_wd1)                          \
if(t_min>t_max){matter_get_tij_limits(index_wd1,index_wd2);}        \
else{                                                               \
  double temp_t_min = t_min; double temp_t_max=t_max;               \
  matter_get_tij_limits(index_wd1,index_wd2);                       \
  t_min = MIN(temp_t_min,t_min); t_max = MAX(temp_t_max,t_max);     \
}
/* macro for re-allocating memory, returning error if it failed */
#define class_realloc_parallel(pointer, newname, size, error_message_output)  {                                  \
  if(abort==_FALSE_){                                                                                            \
    pointer=realloc(newname,size);                                                                               \
    if (pointer == NULL) {                                                                                       \
      int size_int;                                                                                              \
      size_int = size;                                                                                           \
      class_alloc_message(error_message_output,#pointer, size_int);                                              \
      abort=_TRUE_;                                                                                              \
    }                                                                                                            \
  }                                                                                                              \
}
/**
 * Anisotropy matter power spectra \f$ C_l\f$'s for all types, windows and initial conditions.
 * The mode is always scalar.
 *
 * This routine evaluates all the \f$C_l\f$'s at a given value of l by
 * interpolating in the pre-computed table. When relevant, it also
 * sums over all initial conditions.
 *
 * This function can be called from whatever module at whatever time, provided that
 * matter_init() has been called before, and matter_free() has not
 * been called yet.
 *
 * @param pma        Input: pointer to matter structure (containing pre-computed table)
 * @param l          Input: multipole number
 * @param cl_tot     Output: total \f$C_l\f$'s for all types (d_n d_m, L_n L_m)
 * @param cl_ic      Output: \f$C_l\f$'s for all types (d_n d_m, L_n L_m) decomposed by pairs of initial conditions (adiabatic, isocurvatures) when relevant
 * @return the error status
 */

int matter_cl_at_l(
                    struct matters* pma,
                    double l,
                    double * cl_tot,    // array with argument cl_tot[index_cltp*pma->num_window_grid+index_wd_grid] (must be already allocated)
                    double ** cl_ic      // array with argument cl_ic[index_ic1_ic2][index_cltp*pma->num_window_grid+index_wd_grid] (must be already allocated for a given mode only if several ic's)
                    ) {
  class_test(pma->has_cls && pma->l_size<=0,pma->error_message,"Matter was never calculated. Cannot obtain Cl's");

  /**
   * Initialize local variables
   * */

  int last_index;
  int index_ic1,index_ic2,index_ic1_ic2;
  int index_cltp;
  int index_wd_grid;
  last_index = 0;

  /**
   * (a) treat case in which there is only one initial condition.
   * Then, only cl_tot needs to be filled.
   *
   * We set those Cl's above l_max to 0
   * */

  if (pma->ic_size <= 1) {
    if ((int)l <= pma->l_sampling[pma->l_size-1]) {
      for (index_cltp=0; index_cltp<pma->cltp_size; index_cltp++){
        /**
         * Interpolate at given value of l for all cltps
         * */
        class_call(matter_interpolate_spline_growing_hunt(
                                            pma->l_sampling,
                                            pma->l_size,
                                            pma->cl[index_cltp],
                                            pma->ddcl[index_cltp],
                                            //(pma->num_windows*(pma->num_windows+1))/2,
                                            pma->num_window_grid,
                                            l,
                                            &last_index,
                                            cl_tot+index_cltp*pma->num_window_grid,
                                            _FALSE_,
                                            pma->error_message),
                   pma->error_message,
                   pma->error_message);
      }
      //End cltp
    }
    else {
      /** Set zero elements to zero */
      for (index_cltp=0; index_cltp<pma->cltp_size; index_cltp++){
        for(index_wd_grid=0; index_wd_grid < pma->num_window_grid; ++ index_wd_grid){
          cl_tot[index_cltp]=0.;
        }
        //End wd grid
      }
      //end cltp
    }
    //Ifend l<lmax
  }
  /**
   * (b) treat case in which there are several initial condition.
   *  Fill cl_ic and sum it to get cl_tot.
   *
   * We set those Cl's above l_max to 0
   * */
  else{
    for (index_cltp=0; index_cltp<pma->cltp_size; index_cltp++)
      cl_tot[index_cltp]=0.;
    for (index_ic1 = 0; index_ic1 < pma->ic_size; index_ic1++) {
      for (index_ic2 = index_ic1; index_ic2 < pma->ic_size; index_ic2++) {
        index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,pma->ic_size);
        if (((int)l <= pma->l_sampling[pma->l_size-1]) &&
            (pma->is_non_zero[index_ic1_ic2] == _TRUE_)) {

          for (index_cltp=0; index_cltp<pma->cltp_size; index_cltp++){
            /**
             * Interpolate at given value of l and ic combination for all cltps
             * */
            class_call(matter_interpolate_spline_growing_hunt(
                                                pma->l_sampling,
                                                pma->l_size,
                                                pma->cl[index_ic1_ic2*pma->cltp_size+index_cltp],
                                                pma->ddcl[index_ic1_ic2*pma->cltp_size+index_cltp],
                                                //(pma->num_windows*(pma->num_windows+1))/2,
                                                pma->num_window_grid,
                                                l,
                                                &last_index,
                                                cl_ic[index_ic1_ic2]+index_cltp*pma->num_window_grid,
                                                _FALSE_,
                                                pma->error_message),
                       pma->error_message,
                       pma->error_message);
          }
          //End cltp
        }
        else {
          for (index_cltp=0; index_cltp<pma->cltp_size; index_cltp++){
            for(index_wd_grid=0;index_wd_grid<pma->num_window_grid;++index_wd_grid){
              cl_ic[index_ic1_ic2][index_cltp*pma->num_window_grid+index_wd_grid]=0.;
            }
            //End wd grid
          }
          //End cltp
        }
        //Ifend l<lmax and nonzero

        /**
         * We compute the total Cl's by summing over the different initial conditions
         * */
        for (index_cltp=0; index_cltp<pma->cltp_size; index_cltp++) {
          for(index_wd_grid=0;index_wd_grid<pma->num_window_grid;++index_wd_grid){
            if (index_ic1 == index_ic2){
              cl_tot[index_cltp*pma->num_window_grid+index_wd_grid]+=cl_ic[index_ic1_ic2][index_cltp*pma->num_window_grid+index_wd_grid];
            }else{
              cl_tot[index_cltp*pma->num_window_grid+index_wd_grid]+=2.*cl_ic[index_ic1_ic2][index_cltp*pma->num_window_grid+index_wd_grid];
            }
            //Ifend ic1==ic2
          }
          //End wd grid
        }
        //End cltp
      }
      //End ic2
    }
    //End ic1
  }
  //Ifend only single ic
  return _SUCCESS_;
}
/**
 * Initiate the matter construct
 *
 */
int matter_init(
                struct precision * ppr,
                struct background * pba,
                struct thermo * pth,
                struct perturbs * ppt,
                struct primordial * ppm,
                struct nonlinear * pnl,
                struct matters * pma
             ){
  double start_time_omp = omp_get_wtime();
  clock_t start_time = clock();
  double clocksecs;
  if (pma->has_cls == _FALSE_){
    /** In this case, the user requested deactivation of this module */
    return _SUCCESS_;
  }
  if (ppt->has_cls == _FALSE_ || ppt->has_scalars == _FALSE_){
    if(ppt->has_scalars == _FALSE_){
      if (pma->matter_verbose > 0)
      printf("No scalar modes requested. Matter module skipped.\n");
    }
    else{
      if (pma->matter_verbose > 0)
        printf("No Cl spectra requested. Matter module skipped.\n");
    }
    return _SUCCESS_;
  }
  else if(!(ppt->has_cl_number_count||ppt->has_cl_lensing_potential)){
    if (pma->matter_verbose > 0)
      printf("No NumberCount/Shear Cl spectra requested. Matter module skipped.\n");
    return _SUCCESS_;
  }
  else if (pma->matter_verbose > 0){
    fprintf(stdout,"Computing matter spectra\n");
  }
  if(pma->matter_verbose > MATTER_VERBOSITY_PARAMETERS){
    printf(" -> Verbosity set to %i \n",pma->matter_verbose);
  }
  //class_test(ppt->selection == dirac,
  //           pma->error_message,
  //           "Dirac Windows currently not yet implemented -- You can use a very thin gaussian instead.");
  class_test(ppt->selection == dirac && ppt->has_nc_rsd,
             pma->error_message,
             "Including redshift space distortions for dirac functions not yet implemented.");
  class_test(ppt->selection == tophat && ppt->has_nc_rsd,
             pma->error_message,
            "Including redshift space distortions for tophat functions not yet implemented.");
  /**
   * Parameters taken from other structs
   * These should be assigned BEFORE any call to matter_obtain_indices
   * */
  pma->ic_size = ppt->ic_size[ppt->index_md_scalars];
  pma->ic_ic_size = (pma->ic_size*(pma->ic_size+1))/2;
  pma->num_windows = ppt->selection_num;
  pma->num_window_grid = (pma->num_windows*(pma->num_windows+1))/2;
  pma->has_cls = ppt->has_cls;

  /**
   * Flags for execution and approximation flags
   * */

  pma->uses_density_splitting = _FALSE_;
  pma->uses_all_l_sampling = _FALSE_;
  pma->uses_lensing_reduction = _TRUE_;
  pma->uses_rsd_combination = _TRUE_;
  pma->uses_limber_approximation = _FALSE_;
  pma->uses_relative_factors = _FALSE_;

  pma->uses_integration = matter_integrate_tw_t;


  pma->uses_bessel_recursion = _TRUE_;
  //TODO :: do automatic assignment except for when manual requested, then class_test combinations

  class_call(matter_obtain_indices(ppm,ppt,pma),
             pma->error_message,
             pma->error_message);
  class_call(matter_obtain_bi_indices(pma),
             pma->error_message,
             pma->error_message);

  //ALWAYS REMEMBER TO MAKE RECLASS FOR BIAS PROBLEMS
  //ANY CHANGE TO matter.h REQUIRES RECLASS (independently of matter.c changes)

  class_alloc(pma->short_pvecback,
              pba->bg_size_short*sizeof(double),
              pma->error_message);

  pma->h = pba->h;
  pma->tau0 = pba->conformal_age;

  pma->k_max_extr = ppr->matter_k_max_extrapolation;

  pma->ptw_size = ppr->matter_window_preparation_size;//800;//1200;//800;
  //IMPORTANT TO BE HIGH ptw size
  pma->ptw_integrated_size = ppr->matter_integrated_window_preparation_size;//1600;//1600;//160;
  //IMPORTANT TO BE HIGH ptwInt size
  //Tw size: 30 for lens and dens
  //IntTw size: 75 for lens and dens
  //t spline size: 50 for lens, 20 for dens

  pma->bi_wanted_samples = ppr->matter_bi_sample_size+1;
  pma->bessel_recursion_t_size = ppr->matter_bi_sample_size;

  pma->tau_size = MIN(2*(int)(0.5*ppt->tau_size+0.25),pma->tau_size_max);

  //TODO :: FIX CAN ONLY BE POWER OF 2 (?)
  pma->size_fft_result = pma->size_fft_cutoff;
  pma->bi_maximal_t_offset = ppr->matter_t_offset;

  //Offset for keeping numerical instabilities for
  //log(x) small for x->0 (e.g. exp(log(x1)-log(x2) ) = 0 for x1 approx = 0, not NaN
  //Also used for setting exp(log(x))!=x to a value that does not over/underflow the bounds

  pma->small_log_offset = ppr->matter_chi_offset;

  pma->k_weight_k_max = ppr->matter_k_weight_kmax;
  pma->k_weight_k_min = ppr->matter_k_weight_kmin;
  pma->k_weight_mode = ppr->matter_k_weight_mode;

  if(pma->matter_verbose > MATTER_VERBOSITY_PARAMETERS){
    printf("Parameter options as follows : \n");
    printf(" -> Parameter %s has value %s \n","uses integration",(pma->uses_integration==matter_integrate_tw_t?"tw_t":"tw_logt"));
    printf(" -> Parameter %s has value %s \n","uses seperability",(pma->uses_separability?"TRUE":"FALSE"));
    printf(" -> Parameter %s has value %s \n","allow extrapolation",(pma->allow_extrapolation?"TRUE":"FALSE"));
    printf("Integration method : %s \n\n",
         (pma->uses_integration==matter_integrate_tw_t?"tw t":"tw logt"));
    printf("-> Number of windows = %i\n",pma->num_windows);
    printf("-> tau0 = %.20e \n",pma->tau0);
    printf("-> h = %f \n",pma->h);
  }

  /**
   * Test wether the defined combinations would give a valid caclulation
   * */
  //TODO :: class_test if parameters are all correct (i.e. >0 )

  class_test(pma->size_fft_cutoff>pma->size_fft_result,
             pma->error_message,
             "the coefficient cutoff size (%i) has to be smaller or equal to the result size (%i)",
             pma->size_fft_cutoff,pma->size_fft_result);
  /*class_test(pnl->method != nl_none && pma->uses_separability,
             pma->error_message,
             "Cannot assume separability and nonlinear methods at the same time.");
  */
  class_test(pma->tau_size%2!=0,
             pma->error_message,
             "The tau_size parameter currently has to be a multiple of 2");

  /**
   * Obtain samplings in k and tau
   * */
  class_alloc(pma->logk_sampling,
              pma->size_fft_input*sizeof(double),
              pma->error_message);
  class_alloc(pma->k_sampling,
              pma->size_fft_input*sizeof(double),
              pma->error_message);
  class_call(matter_obtain_k_sampling(ppt,pma),
             pma->error_message,
             pma->error_message);

  class_call(matter_obtain_tau_sampling(ppr,pba,ppt,pma),
             pma->error_message,
             pma->error_message);

  /**
   * Obtain primordial spectrum and sources
   * */
  clock_t point_1_time = clock();
  double ** prim_spec;
  class_alloc(prim_spec,
              pma->ic_ic_size*sizeof(double*),
              pma->error_message);
  class_call(matter_obtain_primordial_spectrum(ppt,ppm,pma,prim_spec),
             pma->error_message,
             pma->error_message);

  double ** sources;
  class_alloc(sources,
              pma->ic_size*pma->stp_size*sizeof(double*),
              pma->error_message);
  class_call(matter_obtain_perturbation_sources(pba,ppt,pnl,pma,sources),
             pma->error_message,
             pma->error_message);

  /**
   * Obtain the growth factor from the sources
   * */
  clock_t point_2_time = clock();
  double* k_weights;
  class_alloc(k_weights,
              pma->k_size*sizeof(double),
              pma->error_message);
  class_call(matter_obtain_growth_factor_k_weights(
                                                   pma,
                                                   ppt,
                                                   k_weights
                                                   ),
             pma->error_message,
             pma->error_message);
  class_alloc(pma->growth_factor_tau,
              pma->ic_size*pma->stp_size*sizeof(double*),
              pma->error_message);
  class_call(matter_obtain_growth_factor(
                                         pma,
                                         sources,
                                         k_weights
                                         ),
            pma->error_message,
            pma->error_message);
  if(pma->uses_relative_factors){
    class_call(matter_obtain_relative_factor(
                                             pma,
                                             sources,
                                             k_weights
                                             ),
              pma->error_message,
              pma->error_message);
  }
  free(k_weights);
  class_alloc(pma->ddgrowth_factor_tau,
              pma->ic_size*pma->stp_size*sizeof(double*),
              pma->error_message);
  class_call(matter_spline_growth_factor(
                                         pma
                                        ),
            pma->error_message,
            pma->error_message);
  clock_t point_3_time = clock();
  /**
   * Extrapolate the sources
   * */
  double* perturbed_k_sampling;
  perturbed_k_sampling = ppt->k[ppt->index_md_scalars];
  if(pma->allow_extrapolation){

    double ** extrapolated_sources;
    class_alloc(extrapolated_sources,
                pma->ic_size*pma->stp_size*sizeof(double*),
                pma->error_message);

    class_call(matter_extrapolate_sources(pba,
                                          ppr,
                                          ppt,
                                          pma,
                                          &perturbed_k_sampling,
                                          sources,
                                          extrapolated_sources,
                                          pma->extrapolation_type),
               pma->error_message,
               pma->error_message);

    class_call(matter_free_perturbation_sources(ppt,pnl,pma,sources),
              pma->error_message,
              pma->error_message);
    /**
     * Here we use a special trick:
     * We allow the now empty sources pointer to point to the filled
     * extrapolated_sources array.
     *
     * Thus, the rest of the program can proceed in exactly the same
     * way as before, just that now sources points to a slightly
     * longer array.
     *
     * The deallocation of the extrapolated_sources array in physical memory
     * is being done by the deallocation of the sources pointer.
     *
     * The updating to the new k_size is already done in the extrapolation function.
     * */

     sources = extrapolated_sources;
  }

  /**
   * Now we sample the sources (extrapolated or not)
   *  in the desired k and tau sampling
   * */
  double** sampled_sources;
  class_alloc(sampled_sources,
              pma->ic_size*pma->stp_size*sizeof(double*),
              pma->error_message);
  class_call(matter_sample_sources(ppt,pma,sources,sampled_sources,perturbed_k_sampling),
             pma->error_message,
             pma->error_message);
  /**
   * We free those sources that are not required anymore
   * */
  class_call(matter_free_perturbation_sources(ppt,pnl,pma,sources),
             pma->error_message,
             pma->error_message);
  if(pma->allow_extrapolation){
    free(perturbed_k_sampling);
  }
  clock_t point_4_time = clock();

  /**
   * All sources are obtained,
   *  we can proceed with calculating the FFT in
   *  logarithmic k space
   * */
  double** fft_coeff_real;
  double** fft_coeff_imag;
  class_alloc(fft_coeff_real,
              pma->ic_ic_size*pma->stp_grid_size*sizeof(double*),
              pma->error_message);
  class_alloc(fft_coeff_imag,
              pma->ic_ic_size*pma->stp_grid_size*sizeof(double*),
              pma->error_message);
  class_call(matter_obtain_coeff_sampling(pma),
             pma->error_message,
             pma->error_message);
  class_call(matter_FFTlog_perturbation_sources_parallel(pba,ppt,ppm,pma,sampled_sources,prim_spec,fft_coeff_real,fft_coeff_imag),
             pma->error_message,
             pma->error_message);
  if((!pma->uses_separability)){
    /**
     * This function !replaces! the fft coefficients with their (nearly) constant counterparts
     * The originals can be obtained by re-multiplying with the growth factors
     * */
    class_call(matter_obtain_nonseparability(
                                          pma,
                                          fft_coeff_real,
                                          fft_coeff_imag),
               pma->error_message,
               pma->error_message);

  }
  /**
   * After the FFT, the initial sources are no longer needed,
   *  the information is now carried by the coefficients
   *  ( and/or growth factors )
   * */
  class_call(matter_free_perturbation_sources(ppt,pnl,pma,sampled_sources),
             pma->error_message,
             pma->error_message);
  class_call(matter_free_primordial(ppm,pma,prim_spec),
             pma->error_message,
             pma->error_message);

  class_call(matter_obtain_l_sampling(ppr,ppt,pth,pma),
             pma->error_message,
             pma->error_message);

  clock_t point_5_time = clock();

  /**
   * There are two big ways of obtaining the bessel integrals
   *
   * 1) Use the recursion relation of the bessel integrals
   *   This proves to be really fast, and surprisingly even more accurate
   *
   * 2) Using the direct representations through taylor series
   *   This older method proves to become unreliable due to
   *   floating point arithmetics, especially around
   *   high imaginary parts in nu, large l, and t of around 0.9-0.99
   *
   * We first obtain the integrals for a grid of t values,
   *  after which we spline them for exactly those t values that we require evaluation at
   *
   * It is done this way because the points at which evaluation is required
   *  changes depending on the l and nu
   * */
  double bessel_start_omp = omp_get_wtime();
  if(!pma->uses_limber_approximation){
    /**
     * Obtain the bessel integrals
     * */
    if(pma->uses_bessel_recursion){
      class_call(matter_obtain_bessel_recursion_parallel(pma),
                 pma->error_message,
                 pma->error_message);
    }else{
      class_call(matter_obtain_bessel_integrals(pma),
                 pma->error_message,
                 pma->error_message);
    }
    /**
     * Spline bessel integrals
     * */
    class_alloc(pma->ddbi_real,
                pma->tilt_grid_size*sizeof(double**),
                pma->error_message);
    class_alloc(pma->ddbi_imag,
                pma->tilt_grid_size*sizeof(double**),
                pma->error_message);
    if(pma->uses_bessel_recursion){
      class_call(matter_spline_bessel_integrals_recursion(pma),
               pma->error_message,
               pma->error_message);
    }
    else{
      class_call(matter_spline_bessel_integrals(pma),
                 pma->error_message,
                 pma->error_message);
    }
  }
  //Ifend obtain bessel integrals
  double bessel_end_omp = omp_get_wtime();
  clock_t point_6_time = clock();

  /**
   * Finally, we want to prepare the window functions.
   *  The original window functions given by dN/dz,
   *  and then the window functions multiplied with factors corresponding
   *  to the specific types, and integrating
   * */
  class_call(matter_obtain_time_sampling(ppr,ppt,pba,pma),
             pma->error_message,
             pma->error_message);
  class_call(matter_obtain_prepare_windows_parallel(ppr,ppt,pba,pma),
             pma->error_message,
             pma->error_message);
  clock_t point_7_time = clock();

  /**
   * Now we resample the growth factor for later use
   * */
  class_alloc(pma->growth_factor,
              pma->ic_size*pma->stp_size*sizeof(double*),
              pma->error_message);
  class_call(matter_resample_growth_factor(pma),
             pma->error_message,
             pma->error_message);

  class_call(matter_obtain_t_sampling(pma),
             pma->error_message,
             pma->error_message);
  clock_t point_8_time = clock();
  double point_8_time_omp = omp_get_wtime();
  /**
   * Now we have truly assembled all ingredients
   *  to integrate the final C_l's
   * - We have the window functions (including growth factors)
   * - We have the nu
   * - We have the FFT coefficients
   * - We have the bessel integrals
   *
   * Thus we are finally able to obtain the C_l
   * */
  class_call(matter_integrate_cl(ppr,
                                 pba,
                                 ppt,
                                 pma,
                                 fft_coeff_real,
                                 fft_coeff_imag),
             pma->error_message,
             pma->error_message);

  double point_9_time_omp = omp_get_wtime();
  clock_t point_9_time = clock();

  /**
   * Now we can free our fft coefficients
   *  and the window functions
   * */
  class_call(matter_free_fft(pma,
                             fft_coeff_real,
                             fft_coeff_imag),
               pma->error_message,
               pma->error_message);
  class_call(matter_free_prepare_window(pma),
             pma->error_message,
             pma->error_message);
  /**
   * Finally we spline the Cl's to interpolate
   *  for all l's
   * */
  class_call(matter_spline_cls(pma),
             pma->error_message,
             pma->error_message);


  clock_t point_10_time = clock();
  double end_time_omp = omp_get_wtime();

  if(pma->matter_verbose > MATTER_VERBOSITY_PARAMETERS){
    printf("\n\n\n PARAMETER SUMMARY \n\n ");
    printf("Calculationary parameters \n");
    printf(" -> tilt == %.10e \n",pma->bias);
    printf(" -> nu_imag max == %.10e \n",pma->nu_imag[pma->size_fft_result-1]);
    printf(" -> nu_imag step == %.10e \n",pma->nu_imag[1]);
    printf(" -> k_max == %.10e \n",pma->k_sampling[pma->size_fft_input-1]);
    printf(" -> k_min == %.10e \n",pma->k_sampling[0]);
    printf(" -> delta log(k) == %.10e \n \n\n",pma->deltalogk);

    printf("Parameter counts \n");
    printf(" -> Number of types %i \n",pma->stp_size);
    printf(" -> Number of radials %i \n",pma->radtp_size_total);
    printf(" -> Number of bessel integrals %i \n",pma->bitp_size);
    printf(" -> Number of tilts %i \n",pma->tilt_size);

    printf("Parameter options as follows : \n");
    printf(" -> Parameter '%s' has value %s \n","has_cls",(pma->has_cls?"TRUE":"FALSE"));
    printf(" -> Parameter '%s' has value %s \n","uses integration",(pma->uses_integration==matter_integrate_tw_t?"tw_t":"tw_logt"));
    printf(" -> Parameter '%s' has value %s \n","uses seperability",(pma->uses_separability?"TRUE":"FALSE"));
    printf(" -> Parameter '%s' has value %s \n","allow extrapolation",(pma->allow_extrapolation?"TRUE":"FALSE"));
    printf(" -> Parameter '%s' has value %s \n","uses density spltting",(pma->uses_density_splitting?"TRUE":"FALSE"));
    printf(" -> Parameter '%s' has value %s \n","uses bessel recursion",(pma->uses_bessel_recursion?"TRUE":"FALSE"));
    printf(" -> Parameter '%s' has value %s \n","uses intxi_interpolation",(pma->uses_intxi_interpolation?"TRUE":"FALSE"));
    printf(" -> Parameter '%s' has value %s \n","uses intxi_logarithmic",(pma->uses_intxi_logarithmic?"TRUE":"FALSE"));
    printf(" -> Parameter '%s' has value %s \n","uses intxi_symmetric",(pma->uses_intxi_symmetrized?"TRUE":"FALSE"));
    printf(" -> Parameter '%s' has value %s \n","uses intxi_asymptotic",(pma->uses_intxi_asymptotic?"TRUE":"FALSE"));
    //printf(" -> Parameter '%s' has value %s \n","uses analytic bessel",(pma->uses_bessel_analytic_integration?"TRUE":"FALSE"));
    printf(" -> Parameter '%s' has value %s \n","uses all ell",(pma->uses_all_l_sampling?"TRUE":"FALSE"));

    printf(" -> Parameter '%s' has value %s \n","uses RSD combination",(pma->uses_rsd_combination?"TRUE":"FALSE"));
    printf(" -> Parameter '%s' has value %s \n","uses relative factors",(pma->uses_relative_factors?"TRUE":"FALSE"));
    printf(" -> Parameter '%s' has value %s \n","uses limber approximation",(pma->uses_limber_approximation?"TRUE":"FALSE"));
    printf(" -> Parameter '%s' has value %i \n","nondiagonals",pma->non_diag);

    printf(" -> Parameter '%s' has value %i \n","tw size",pma->tw_size);
    printf(" -> Parameter '%s' has value %i \n","tw integrated size",pma->integrated_tw_size);
    printf(" -> Parameter '%s' has value %i \n","t size",pma->t_size);
    printf(" -> Parameter '%s' has value %i \n","max coeff",pma->size_fft_cutoff);

  }
  if(pma->matter_verbose > MATTER_VERBOSITY_TIMING){
    printf("\n\n\n TIMING SUMMARY (for %10d spectra, %5d l values) \n\n",((pma->non_diag+1)*(2*pma->num_windows-pma->non_diag))/2,pma->l_size);//,pma->num_window_grid,pma->l_size);
    clocksecs = ((double)(point_1_time-start_time))/(double)CLOCKS_PER_SEC;
    if(pma->matter_verbose>0){
      printf("Init took                        %15f seconds \n",clocksecs);
    }
    clocksecs = ((double)(point_2_time-point_1_time))/(double)CLOCKS_PER_SEC;
    if(pma->matter_verbose>0){
      printf("Obtaining Sources took           %15f seconds \n",clocksecs);
    }
    clocksecs = ((double)(point_3_time-point_2_time))/(double)CLOCKS_PER_SEC;
    if(pma->matter_verbose>0){
      printf("Obtaining Growth Factor took     %15f seconds \n",clocksecs);
    }
    clocksecs = ((double)(point_4_time-point_3_time))/(double)CLOCKS_PER_SEC;
    if(pma->matter_verbose>0){
      printf("Extrapolate Sources took         %15f seconds \n",clocksecs);
    }
    clocksecs = ((double)(point_5_time-point_4_time))/(double)CLOCKS_PER_SEC;
    if(pma->matter_verbose>0){
      printf("FFT took                         %15f seconds \n",clocksecs);
    }
    clocksecs = ((double)(point_6_time-point_5_time))/(double)CLOCKS_PER_SEC;
    if(pma->matter_verbose>0){
      printf("Obtain Bessel took               %15f seconds, %15f CPU time \n",bessel_end_omp-bessel_start_omp,clocksecs);
    }
    clocksecs = ((double)(point_7_time-point_6_time))/(double)CLOCKS_PER_SEC;
    if(pma->matter_verbose>0){
      printf("Obtain Windows took              %15f seconds \n",clocksecs);
    }
    clocksecs = ((double)(point_8_time-point_7_time))/(double)CLOCKS_PER_SEC;
    if(pma->matter_verbose>0){
      printf("Resample Growth Factor took      %15f seconds \n",clocksecs);
    }
    /*clocksecs = ((double)(point_9_time-point_8_time))/(double)CLOCKS_PER_SEC;
    if(pma->matter_verbose>0){
      printf("Integrating took                 %15f seconds (%10f per window) \n",clocksecs,clocksecs/pma->num_window_grid);
    }*/
    clocksecs = ((double)(point_9_time-point_8_time))/(double)CLOCKS_PER_SEC;
    if(pma->matter_verbose>0){
      printf("Integrating took                 %15f seconds (%10f per window), %15f CPU time \n",point_9_time_omp-point_8_time_omp,(point_9_time_omp-point_8_time_omp)/(((pma->non_diag+1)*(2*pma->num_windows-pma->non_diag))/2),clocksecs);
    }
    clocksecs = ((double)(point_10_time-point_9_time))/(double)CLOCKS_PER_SEC;
    if(pma->matter_verbose>0){
      printf("Cl spline and print took         %15f seconds \n",clocksecs);
    }
    clocksecs = ((double)(point_10_time-start_time))/(double)CLOCKS_PER_SEC;
    if(pma->matter_verbose>0){
      printf("All points excluding output took %15f seconds, %15f CPU time \n",end_time_omp-start_time_omp,clocksecs);
    }
  }
  clock_t end_time = clock();
  clocksecs = ((double)(end_time-start_time))/(double)CLOCKS_PER_SEC;
  if(pma->matter_verbose>MATTER_VERBOSITY_TIMING){
    printf("Matter took %f seconds \n",clocksecs);
  }
  return _SUCCESS_;
}











int matter_free(
                struct matters * pma
              ) {
  int i,j;
  if(pma->has_cls){
    if(pma->matter_verbose>2){
      printf("Matter starts to free \n");
    }
    if(pma->matter_verbose>2){
      printf("Deleting fft stuff \n");
    }
    free(pma->logk_sampling);
    free(pma->k_sampling);
    free(pma->tau_sampling);

    if(pma->matter_verbose>2){
      printf("Deleting window stuff \n");
    }
    free(pma->tw_sampling);
    free(pma->tw_weights);
    free(pma->integrated_tw_sampling);
    free(pma->integrated_tw_weights);
    free(pma->exp_integrated_tw_sampling);

    free(pma->tw_max);
    free(pma->tw_min);

    free(pma->t_sampling);
    free(pma->t_weights);

    if(pma->uses_intxi_interpolation){
      free(pma->t_spline_sampling);
    }
    for(i=0;i<pma->ic_size*pma->stp_size;++i){
      free(pma->growth_factor_tau[i]);
      free(pma->growth_factor[i]);
      free(pma->ddgrowth_factor_tau[i]);
    }
    free(pma->growth_factor_tau);
    free(pma->growth_factor);
    free(pma->ddgrowth_factor_tau);


    if(pma->matter_verbose>2){
      printf("Deleting stuff of method %s \n",(pma->uses_integration==matter_integrate_tw_t?"tw t":"tw logt"));
    }


    if(pma->matter_verbose>2){
      printf("Deleting bessel stuff \n");
    }
    free(pma->nu_real);
    free(pma->nu_imag);
    if(!pma->uses_limber_approximation){
      for(j=0;j<pma->tilt_grid_size;++j){
        if(pma->uses_bessel_recursion){
          for(i=0;i<pma->l_size_recursion*pma->size_fft_cutoff;++i){
            int delta = (i/pma->size_fft_cutoff);
            int index = delta*pma->size_fft_result+(i-pma->size_fft_cutoff*delta);
            //printf("Freeing %4d %4d %4d \n",j,delta,i-pma->size_fft_cutoff*delta);
            free(pma->bi_real[j][index]);
            free(pma->bi_imag[j][index]);
            free(pma->bi_sampling[j][index]);
            free(pma->ddbi_real[j][index]);
            free(pma->ddbi_imag[j][index]);
          }
        }
        else{
          for(i=0;i<pma->size_fft_result*pma->l_size;++i){
            free(pma->bi_real[j][i]);
            free(pma->bi_imag[j][i]);
            free(pma->bi_sampling[j][i]);
            free(pma->ddbi_real[j][i]);
            free(pma->ddbi_imag[j][i]);
          }

        }
        free(pma->bi_real[j]);
        free(pma->bi_imag[j]);
        free(pma->bi_sampling[j]);
        free(pma->ddbi_real[j]);
        free(pma->ddbi_imag[j]);
        free(pma->bi_size[j]);
        free(pma->bi_max[j]);
      }
      free(pma->bi_real);
      free(pma->bi_imag);
      free(pma->ddbi_real);
      free(pma->ddbi_imag);
      free(pma->bi_sampling);
      free(pma->bi_size);
      free(pma->bi_max);
    }

    if(pma->matter_verbose>2){
      printf("Deleting general stuff \n");
    }
    free(pma->l_sampling);
    free(pma->short_pvecback);
    free(pma->index_perturb_tp_of_stp);
    free(pma->index_stp_of_radtp);
    if(pma->has_cltp_nc){
      if(pma->has_bitp_normal){
        free(pma->radtps_of_bitp[pma->cltp_index_nc*pma->bitp_size+pma->bitp_index_normal]);
      }
      if(pma->has_bitp_nu_reduced){
        free(pma->radtps_of_bitp[pma->cltp_index_nc*pma->bitp_size+pma->bitp_index_nu_reduced]);
      }
      if(pma->has_bitp_lfactor){
        free(pma->radtps_of_bitp[pma->cltp_index_nc*pma->bitp_size+pma->bitp_index_lfactor]);
      }
    }
    if(pma->has_cltp_sh){
      if(pma->has_bitp_lfactor){
        free(pma->radtps_of_bitp[pma->cltp_index_sh*pma->bitp_size+pma->bitp_index_lfactor]);
      }
    }
    free(pma->radtps_of_bitp);
    free(pma->radtp_of_bitp_size);

    free(pma->is_non_zero);
    for(i=0;i<pma->ic_ic_size;++i){
      for(j=0;j<pma->cltp_size;++j){
        free(pma->cl[i*pma->cltp_size+j]);
        free(pma->ddcl[i*pma->cltp_size+j]);
      }
    }
    free(pma->cl);
    free(pma->ddcl);
  }
  return _SUCCESS_;
}
int matter_free_primordial(
                           struct primordial* ppm,
                           struct matters * pma,
                           double ** prim_spec
                           ) {
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS) {
    printf("Method :: Deleting primordial \n");
  }
  int index_ic1_ic2;
  for (index_ic1_ic2 = 0; index_ic1_ic2 < pma->ic_ic_size; index_ic1_ic2++) {
    free(prim_spec[index_ic1_ic2]);
  }
  free(prim_spec);
  return _SUCCESS_;
}
int matter_free_fft(
                   struct matters * pma,
                   double ** fft_real,
                   double ** fft_imag
                   ) {
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS) {
    printf("Method :: Deleting fft \n");
  }
  int index_ic1_ic2;
  int index_stp1_stp2;
  for(index_ic1_ic2=0;index_ic1_ic2<pma->ic_ic_size;++index_ic1_ic2){
    for(index_stp1_stp2=0;index_stp1_stp2<pma->stp_grid_size;++index_stp1_stp2){
      free(fft_real[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2]);
      free(fft_imag[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2]);
    }
  }
  free(fft_real);
  free(fft_imag);
  return _SUCCESS_;
}
int matter_free_perturbation_sources(
                                     struct perturbs * ppt,
                                     struct nonlinear * pnl,
                                     struct matters * pma,
                                     double ** sources
                                     ) {
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS) {
    printf("Method :: Deleting perturbation sources \n");
  }
  int index_ic;
  int index_stp;
  for (index_ic = 0; index_ic < pma->ic_size; index_ic++) {
    for (index_stp = 0; index_stp < pma->stp_size; index_stp++) {
      free(sources[index_ic * pma->stp_size + index_stp]);
    }
  }
  free(sources);
  return _SUCCESS_;
}
int matter_free_prepare_window(
                               struct matters* pma
                              ){
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS){
    printf("Method :: Freeing windows including growth factors etc. \n");
  }
  int index_delete;
  for(index_delete=0;index_delete<pma->radtp_size_total*pma->ic_size;++index_delete){
    free(pma->ptw_window[index_delete]);
    free(pma->ptw_dwindow[index_delete]);
    free(pma->ptw_ddwindow[index_delete]);
  }
  free(pma->ptw_window);
  free(pma->ptw_dwindow);
  free(pma->ptw_ddwindow);

  free(pma->ptw_sampling);
  free(pma->ptw_weights);
  free(pma->ptw_orig_window);

  free(pma->ptw_integrated_sampling);
  free(pma->ptw_integrated_weights);
  return _SUCCESS_;
}






int matter_spline_cls(
                      struct matters* pma
                      ){
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS){
    printf("Method :: Splining the final Cl's \n");
  }
  int index_cltp;
  int index_ic_ic;
  int index_wd1_wd2;
  class_alloc(pma->ddcl,
              pma->ic_ic_size*pma->cltp_size*sizeof(double*),
              pma->error_message);

  for(index_cltp=0;index_cltp<pma->cltp_size;++index_cltp){
    for(index_ic_ic=0;index_ic_ic<pma->ic_ic_size;++index_ic_ic){
      class_alloc(pma->ddcl[index_ic_ic*pma->cltp_size+index_cltp],
                  ((pma->num_windows+1)*pma->num_windows)/2*pma->l_size*sizeof(double),
                  pma->error_message);
      for(index_wd1_wd2=0;index_wd1_wd2<((pma->num_windows+1)*pma->num_windows)/2;++index_wd1_wd2){
        array_spline_table_columns(
                 pma->l_sampling,
                 pma->l_size,
                 pma->cl[index_ic_ic*pma->cltp_size+index_cltp]+index_wd1_wd2*pma->l_size, // array of size x_size*y_size with elements
                 1,
                 pma->ddcl[index_ic_ic*pma->cltp_size+index_cltp]+index_wd1_wd2*pma->l_size,
                 _SPLINE_EST_DERIV_,
                 pma->error_message
                 );
      }
      //End wd
    }
    //End ic_ic
  }
  //End cl tp
  return _SUCCESS_;
}








int matter_obtain_coeff_sampling(
                               struct matters * pma
                              ){
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS){
    printf("Method :: Obtaining FFT coefficient sampling\n");
  }
  int index_coeff;
  int index_tilt1,index_tilt2,index_tilt1_tilt2;
  double current_tilt_offset;
  class_alloc(pma->nu_real,
              pma->tilt_grid_size*sizeof(double),
              pma->error_message);
  /**
   * The real part of the coefficient
   *  depends only on the tilt
   * We have to iterate thus through every
   *  possible combination of tilts.
   * */
  for(index_tilt1=0;index_tilt1<pma->tilt_size;++index_tilt1){
    for(index_tilt2=index_tilt1;index_tilt2<pma->tilt_size;++index_tilt2){
      index_tilt1_tilt2 = index_symmetric_matrix(index_tilt1,index_tilt2,pma->tilt_size);
      if(matter_is_index(index_tilt1,pma->tilt_index_normal,pma->has_tilt_normal)){
        if(matter_is_index(index_tilt2,pma->tilt_index_normal,pma->has_tilt_normal)){
          current_tilt_offset = 0.0;
        }
        else if(matter_is_index(index_tilt2,pma->tilt_index_reduced,pma->has_tilt_reduced)){
          current_tilt_offset = 2.0;
        }
        else{
          class_stop(pma->error_message,"Tilt index %i not recognized",index_tilt2);
        }
      }
      else if(matter_is_index(index_tilt1,pma->tilt_index_reduced,pma->has_tilt_reduced)){
        if(matter_is_index(index_tilt2,pma->tilt_index_normal,pma->has_tilt_normal)){
          current_tilt_offset = 2.0;
        }
        else if(matter_is_index(index_tilt2,pma->tilt_index_reduced,pma->has_tilt_reduced)){
          current_tilt_offset = 4.0;
        }
        else{
          class_stop(pma->error_message,"Tilt index %i not recognized",index_tilt2);
        }
      }
      else{
        class_stop(pma->error_message,"Tilt index %i not recognized",index_tilt1);
      }
      if(pma->matter_verbose > MATTER_VERBOSITY_RANGES){
        printf(" -> Found nu_real = %f (offset from bias : %f) \n",pma->bias-current_tilt_offset,current_tilt_offset);
      }
      pma->nu_real[index_tilt1_tilt2]=pma->bias-current_tilt_offset;
    }
    //End tilt2
  }
  //End tilt1
  /**
   * The imaginary part of the coefficients
   *   depends only on the coefficient index
   * */
  class_alloc(pma->nu_imag,
              pma->size_fft_result*sizeof(double),
              pma->error_message);
  for(index_coeff=0;index_coeff<pma->size_fft_result;++index_coeff){
    /**
     * The factor of (N-1)/N might at first seem confusing, but it is necessary and mathematically correct:
     *
     * For any FFT, we want factors of exp(2*pi*i*m*n)
     *  In our case, the FFT goes over log(k),
     *  which was sampled as k_m = k_0 * exp(m/(N-1)*dkap)
     * (where dkap = log(k_max)-log(k_min)
     *  Notice the N-1. This is included to let k_(N-1) = k_max
     *
     * However, this (N-1) factor is not the one required by the FFT exponential
     *  The (N-1)/N is a sort of correction for this fact
     *
     * For this, let us calculate k_m*k^(nu_imag_n)
     *  k_m = k_0 * exp(m/(N-1)*dkap)
     *  k^(nu_imag_n) = exp(2*pi*n/dkap *(N-1)/N)
     *
     * =>
     *  k_m k^(nu_imag_n) = k_0*exp(2*pi*m*n/N)
     * Which is exactly of the form we want
     *
     * It was very important here, that the N-1 factor should cancel,
     *  which is only possible if we include this correction factor here
     * */
    pma->nu_imag[index_coeff]=_TWOPI_*(((double)(index_coeff))/(pma->deltalogk))*((double)pma->size_fft_input-1)/((double)pma->size_fft_input);
    /**
     * The algorithm used to obtain the bessel integrals is not always stable,
     *  especially not when the imaginary part of the coefficients gets too large
     * We want to warn the user about this fact
     * */
    if(!pma->uses_bessel_recursion && pma->nu_imag[index_coeff]> WARN_LARGE_NU_IMAG){
      printf("PROBLEM :: nu_imag > %f , for index %i, DLOGK = %f , nu = %f \n \
        The imaginary part of the coefficients was larger than \n \
        what is currently supported without bessel recursion. \n \
        Consider switching to bessel recursion, \n \
        increasing your delta-logk (e.g. by using source extrapolation), \n \
        or increasing the warning parameter 'WARN_LARGE_NU_IMAG' \n",
        WARN_LARGE_NU_IMAG,index_coeff,pma->deltalogk,pma->nu_imag[index_coeff]);
    }
    //Ifend user warning
  }
  //End coeff
  return _SUCCESS_;
}
int matter_obtain_tau_sampling(
                               struct precision* ppr,
                               struct background* pba,
                               struct perturbs * ppt,
                               struct matters * pma
                              ){
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS){
    printf("Method :: Obtaining conformal time sampling \n");
  }
  /**
   * Define variables to be used later
   * */
  double z_max = 0.0;
  double tau_min;
  int index_wd;
  int index_tau;
  int bin;
  /**
   *
   * Find the start in tau (maximal z) of the windows.
   * We always end at tau0.
   *
   * Afterwards generate new tau values with given sampling accuracy.
   * index_tau_perturbs tracks the tau value of perturbations.c
   * index_tau_matter tracks the tau value used in this matter.c
   * index_tau_perturbs_of_tau_matter gives the relation between those two.
   *
   * */
  for(index_wd=0;index_wd<pma->num_windows;++index_wd){
    bin = index_wd;
    if (ppt->selection==gaussian) {
      z_max = MAX(z_max,ppt->selection_mean[bin]+ppt->selection_width[bin]*ppr->selection_cut_at_sigma);
    }
    if (ppt->selection==tophat) {
      z_max = MAX(z_max,ppt->selection_mean[bin]+(1.+ppr->selection_cut_at_sigma*ppr->selection_tophat_edge)*ppt->selection_width[bin]);
    }
    if (ppt->selection==dirac) {
      z_max = MAX(z_max,ppt->selection_mean[bin]+pma->small_log_offset);
    }
  }
  class_call(background_tau_of_z(pba,
                                 z_max,
                                 &tau_min),
             pba->error_message,
             pma->error_message);
  /**
   * The minimum tau has been found
   * */
  if(pma->matter_verbose>MATTER_VERBOSITY_RANGES){
    printf(" -> Found minimum tau for windows : %.10e \n",tau_min);
  }
  /**
   * Find the first tau_perturbs, corresponding to the tau_min
   * */
  pma->index_tau_perturbs_beginning = 0;
  class_call(matter_spline_prepare_hunt(ppt->tau_sampling,
                             ppt->tau_size,
                             tau_min,
                             &(pma->index_tau_perturbs_beginning),
                             pma->error_message),
             pma->error_message,
             pma->error_message);
  if(pma->matter_verbose>MATTER_VERBOSITY_RANGES){
    printf(" -> Corresponding tau perturbation : %.10e (index %i) \n",ppt->tau_sampling[pma->index_tau_perturbs_beginning],pma->index_tau_perturbs_beginning);
  }
  /**
   * With this information we can figure out nthe tau size in perturbs.
   * Using that in conjuction with the sampling factor, we can also figure out the tau size in matter
   * After obtaining the matter tau size, we can find the matter sampling,
   * keeping track of the index relations.
   * */
  if(pma->matter_verbose>MATTER_VERBOSITY_RANGES){
    printf(" -> The tau size of matter was chosen to be : %i values ( perturbation size = %i ) \n",
          pma->tau_size,ppt->tau_size);
  }
  class_alloc(pma->tau_sampling,
              pma->tau_size*sizeof(double),
              pma->error_message);

  pma->tau_sampling[0] = tau_min;
  for(index_tau=1;index_tau<pma->tau_size-1;++index_tau){
    pma->tau_sampling[index_tau] = exp(
                                        log(tau_min)+//285-pma->small_log_offset)+
                                        (log(pma->tau0)-log(tau_min))
                                          *((double)index_tau)/((double)(pma->tau_size-1))
                                      );
  }
  pma->tau_sampling[pma->tau_size-1]=pma->tau0;
  /**
   * We want to make sure that mistakes in exp(log(x))=x due to floating point errors do not arise
   * This is only important at the very edge, where otherwise the array might be outside of splining range
   * */
  pma->tau_sampling[0] = tau_min;
  pma->tau_sampling[pma->tau_size-1] = pma->tau0;
  pma->tau_grid_size = pma->tau_size*pma->tau_size;

  return _SUCCESS_;
}
int matter_obtain_k_sampling(
                               struct perturbs* ppt,
                               struct matters * pma
                              ){
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS){
    printf("Method :: Obtaining k sampling\n");
  }
  double k_min,k_max_total;
  int index_coeff;

  /**
   * We obtain k_min and k_max (storing k_max for possible later extrapolation)
   * */
  k_min = ppt->k_min;
  pma->k_max = ppt->k[ppt->index_md_scalars][ppt->k_size[ppt->index_md_scalars]-1];
  k_max_total = pma->k_max;
  if(pma->matter_verbose > MATTER_VERBOSITY_RANGES){
    printf(" -> Calculated k_max of %.10e \n",pma->k_max);
  }
  if(pma->allow_extrapolation){
    printf(" -> Given extrapolation k_max of %.10e \n",pma->k_max_extr);
    if(pma->k_max_extr>pma->k_max){
      k_max_total = pma->k_max_extr;
    }
    else{
      printf("WARNING :: k_max of extrapolation is smaller than k_max_scalars, skipping extrapolation \n");
      pma->allow_extrapolation = _FALSE_;
    }
  }
  pma->k_size = ppt->k_size[ppt->index_md_scalars];

  class_test(k_min<=0,
             pma->error_message,
             "The calculated k_min %.10e was smaller or equal to 0, and thus the log(k) could not be calculated.",
             k_min);
  class_test(k_max_total<=0,
             pma->error_message,
             "The calculated k_max %.10e was smaller or equal to 0, and thus the log(k) could not be calculated.",
             k_max_total);

  /**
   * We want to smaple k logarithmically,
   *  also the difference in log(k) is important for the calculation of the imaginary part of the frequency nu
   * */
  pma->logmink = log(k_min);
  pma->deltalogk = log(k_max_total)-log(k_min);

  /**
   * We make sure the first and last k are set exactly, since exp(log(x))=x is sometimes not true for floats
   * */
  pma->logk_sampling[0] = pma->logmink;
  pma->k_sampling[0] = ppt->k_min;
  for(index_coeff=1;index_coeff<pma->size_fft_input-1;++index_coeff){
    pma->logk_sampling[index_coeff] = (pma->logmink)+((double)(index_coeff)/(double)(pma->size_fft_input-1))*(pma->deltalogk);
    pma->k_sampling[index_coeff] = exp(pma->logk_sampling[index_coeff]);
  }
  pma->k_sampling[pma->size_fft_input-1] = pma->k_max_extr;
  pma->logk_sampling[pma->size_fft_input-1] = log(pma->k_max_extr);

  if(pma->matter_verbose > MATTER_VERBOSITY_RANGES){
    printf("Wavenumbers k go from %.10e to %.10e with delta(log(k)) = %.10e \n",pma->k_sampling[0],pma->k_sampling[pma->size_fft_input-1],pma->deltalogk);
  }
  return _SUCCESS_;
}
int matter_obtain_l_sampling(
                            struct precision* ppr,
                            struct perturbs* ppt,
                            struct thermo* pth,
                            struct matters * pma
                            ){
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS){
    printf("Method :: Obtaining l sampling\n");
  }
  int index_l;
  int current_l;
  int increment;
  int l_min = 2;
  int l_max = ppt->l_lss_max;

  //The smallest stepsize is 1, so we can safely assume the maximum size being l_max
  class_alloc(pma->l_sampling,
              l_max*sizeof(double),
              pma->error_message);
  if(!pma->uses_all_l_sampling){
    /**
     * This is the normal logarithmic sampling that you will also
     *   see in other parts of class, like e.g. the spectra module
     *
     *
     * We start from l = 2 and increase it with a logarithmic step
     * */

    index_l = 0;
    current_l = l_min;
    increment = MAX((int)(current_l * (pow(ppr->l_logstep,pth->angular_rescaling)-1.)),1);
    pma->l_sampling[index_l]=current_l;

    while (((current_l+increment) < l_max) &&
           (increment < ppr->l_linstep*pth->angular_rescaling)) {

      index_l ++;
      current_l += increment;
      pma->l_sampling[index_l]=current_l;

      increment = MAX((int)(current_l * (pow(ppr->l_logstep,pth->angular_rescaling)-1.)),1);

    }

    /**
     * When the logarithmic step becomes larger than some linear step,
     * stick to this linear step until we reach l_max
     * */

    increment = MAX((int)(ppr->l_linstep*pth->angular_rescaling+0.5),1);

    while ((current_l+increment) <= l_max) {

      index_l ++;
      current_l += increment;
      pma->l_sampling[index_l]=current_l;

    }

    /**
     * The last value has to be set to exactly l_max
     *  (Otherwise there would be out-of-bounds problems with splining)
     * Of course we only need to add the additonal
     *  value of l_max, if we don't already hit it
     *  by accident
     * */

    if (current_l != l_max) {

      index_l ++;
      current_l = l_max;
      pma->l_sampling[index_l]=current_l;

    }

    pma->l_size = index_l+1;
    class_realloc(pma->l_sampling,
                  pma->l_sampling,
                  (index_l+1)*sizeof(double),
                  pma->error_message);
  }
  else{
    /**
     * The l_size is l_max-1,
     *  thus index_l is smaller or equal to l_max-2
     *  thus index_l+2 is smaller or equal to l_max,
     * as desired
     * */
    pma->l_size = l_max-1;
    for(index_l=0;index_l<pma->l_size;++index_l){
      pma->l_sampling[index_l]=index_l+2;
    }
    class_realloc(pma->l_sampling,
                  pma->l_sampling,
                  pma->l_size*sizeof(double),
                  pma->error_message);
  }
  return _SUCCESS_;
}







int matter_obtain_t_sampling(struct matters* pma){
/**
   * We want to do a normal integration,
   *  if we do not use limber approximation
   * Otherwise we only have the value t=1
   * */
  if(!pma->uses_limber_approximation){
    class_alloc(pma->t_sampling,
                pma->t_size*sizeof(double),
                pma->error_message);
    class_alloc(pma->t_weights,
                pma->t_size*sizeof(double),
                pma->error_message);
    class_call(array_weights_gauss_limits(
                                    pma->t_sampling,
                                    pma->t_weights,
                                    0.0,//1e-8 for trapz integration
                                    1.0,
                                    pma->t_size,
                                    gauss_type_legendre,//gauss_type_legendre_half,//gauss_type_trapezoid,//gauss_type_legendre,
                                    pma->error_message),
              pma->error_message,
              pma->error_message);
  }
  else{
    class_alloc(pma->t_sampling,
                sizeof(double),
                pma->error_message);
    class_alloc(pma->t_weights,
                sizeof(double),
                pma->error_message);
    pma->t_size=1;
    pma->t_sampling[0]=1.0;
    pma->t_weights[0]=1.0;
  }
  /**
   * If we want to obtain f_n^{ij}(t) for a subset of
   *  t values and spline it for all others,
   *  we can create a seperate sampling (t_spline_sampling)
   * Since we never integrate over that sampling,
   * we ingore the weights returned from this method
   *
   * Otherwise t spline and t are the same
   * */
  if(pma->uses_intxi_interpolation){
    double* ignore;
    class_alloc(ignore,
                pma->t_spline_size*sizeof(double),
                pma->error_message);
    class_alloc(pma->t_spline_sampling,
                pma->t_spline_size*sizeof(double),
                pma->error_message);
    class_call(array_weights_gauss_limits(
                                    pma->t_spline_sampling,
                                    ignore,
                                    0.0,//+1e-8 ? for trapz integration
                                    1.0,
                                    pma->t_spline_size,
                                    gauss_type_trapezoid,//gauss_type_trapezoid,//gauss_type_legendre,
                                    pma->error_message),
              pma->error_message,
              pma->error_message);
    double a = 0.75;//0.75;lens//1.0;//0.5;
    int index_t;
    for(index_t=0;index_t<pma->t_spline_size;++index_t){
      double t = pma->t_spline_sampling[index_t];
      pma->t_spline_sampling[index_t]=(1.-a)*t+a*(2.*t-t*t);
    }
    free(ignore);
  }
  else{
    pma->t_spline_size = pma->t_size;
    pma->t_spline_sampling = pma->t_sampling;
  }
  //Ifend xi interpolation
  return _SUCCESS_;
}
int matter_obtain_time_sampling(
                          struct precision* ppr,
                          struct perturbs* ppt,
                          struct background* pba,
                          struct matters* pma){
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS){
    printf("Method :: Obtaining time sampling \n");
  }
  int index_wd;
  double tau_min,tau_max;
  double zmin=0,zmax=0;
  /**
   * Allocate time samplings and weights
   *
   * Normal time sampling:
   *  tw = Tau of Window
   *  The tau sampling that is going to be used finally
   * Additionally:
   *  integrated_tw for integrated windows
   *
   * Since derivatives of the window functions are going to be
   * possibly calculated, we might want a higher time sampling accuracy
   * for just that task:
   *  ptw = Preparation Tau of Window
   *  The tau sampling that is used in the preparation
   *  (integration/differentiation etc.) of the window functions
   * Additionally:
   *  ptw_integrated for integrated windows
   *
   * */
  if(ppt->selection==dirac){
    pma->ptw_size=3;
    pma->tw_size=3;
  }
  class_alloc(pma->ptw_sampling,
              pma->ptw_size*pma->num_windows*sizeof(double),
              pma->error_message);
  class_alloc(pma->ptw_weights,
              pma->ptw_size*pma->num_windows*sizeof(double),
              pma->error_message);
  class_alloc(pma->tw_sampling,
              pma->tw_size*pma->num_windows*sizeof(double),
              pma->error_message);
  class_alloc(pma->tw_weights,
              pma->tw_size*pma->num_windows*sizeof(double),
              pma->error_message);
  class_alloc(pma->integrated_tw_sampling,
              pma->integrated_tw_size*pma->num_windows*sizeof(double),
              pma->error_message);
  class_alloc(pma->integrated_tw_weights,
              pma->integrated_tw_size*pma->num_windows*sizeof(double),
              pma->error_message);
  class_alloc(pma->tw_max,
              pma->num_windows*sizeof(double),
              pma->error_message);
  class_alloc(pma->tw_min,
              pma->num_windows*sizeof(double),
              pma->error_message);
  class_alloc(pma->ptw_integrated_sampling,
              pma->ptw_integrated_size*pma->num_windows*sizeof(double),
              pma->error_message);
  class_alloc(pma->ptw_integrated_weights,
              pma->ptw_integrated_size*pma->num_windows*sizeof(double),
              pma->error_message);
  class_alloc(pma->exp_integrated_tw_sampling,
              pma->integrated_tw_size*pma->num_windows*sizeof(double),
              pma->error_message);
  /**
   * For every window we first find the range of the window function
   * And we store this range in tw_min and tw_max
   * */
  for(index_wd=0;index_wd<pma->num_windows;++index_wd){
    /**
     * Find tau_min
     * */
    if (ppt->selection==gaussian) {
      zmin = ppt->selection_mean[index_wd]+ppt->selection_width[index_wd]*ppr->selection_cut_at_sigma;
    }
    if (ppt->selection==tophat) {
      zmin = ppt->selection_mean[index_wd]+(1.+ppr->selection_cut_at_sigma*ppr->selection_tophat_edge)*ppt->selection_width[index_wd];
    }
    if (ppt->selection==dirac) {
      zmin = ppt->selection_mean[index_wd]+pma->small_log_offset;
    }

    class_call(background_tau_of_z(pba,
                                   zmin,
                                   &tau_min),
               pba->error_message,
               pma->error_message);
    pma->tw_min[index_wd]=tau_min;

    /**
     * Find tau_max (and keep it smaller than tau0)
     * */
    if (ppt->selection==gaussian) {
      zmax = ppt->selection_mean[index_wd]-ppt->selection_width[index_wd]*ppr->selection_cut_at_sigma;
    }
    if (ppt->selection==tophat) {
      zmax = ppt->selection_mean[index_wd]-(1.+ppr->selection_cut_at_sigma*ppr->selection_tophat_edge)*ppt->selection_width[index_wd];
    }
    if (ppt->selection==dirac) {
      zmax = ppt->selection_mean[index_wd]-pma->small_log_offset;
    }

    /**
     * If the window function vanishes at exactly 0,
     * we need to make sure that log(chi_min) = log(0) doesn't give a
     * numerical divergence. Hence this small offset.
     *
     * If the window function does not vanish at 0, and reaches instead
     * into the regime z<0, a mathematically real but unphysical divergence
     * appears in our equations.
     *
     * (Reminder: t = chi_2/chi_1, chi=chi_1
     *  The limit t = fixed, and chi->0
     * is neither the same as chi_1 -> 0 with chi_2 = fixed
     *        nor the same as chi_2 -> 0 with chi_1 = fixed
     *
     * Instead, the function chi^(-nu) I_ell(nu,t) has a different limit
     * depending on t fixed or t loose.
     * The reason for this is the unphysicalness of the transformation
     * u = k chi_1 as the variable of integration for I_ell for chi_1->0
     * or similarly for ut = k chi_2 for chi_2->0
     *
     * If t is not held fixed, either chi_1 -> 0 corresponds to t->inf
     * leading to the function vanishing, or chi_2 ->0 corresponds to t->0
     * leading also to the function vanishing. Only the t = fixed case
     * is problematic
     *
     * We will investigate in the future different approaches that
     * do not substitute the variables as mentioned.
     * */
    if(zmax ==0){zmax+=pma->small_log_offset;}
    class_test(zmax<0,pma->error_message,
               "\n\
                The window function cannot reach z=0, otherwise numerical \n\
                divergencies (which are unphysical) will appear. This is \n\
                a matter of problem reformulation, which is currently \n\
                being worked on to be solved. For now, please choose \n\
                a window function which vanishes at z=0, \n\
                and make sure that this function calculates the correct \n\
                maximal time value.");

    class_call(background_tau_of_z(pba,
                                   zmax,
                                   &tau_max),
               pba->error_message,
               pma->error_message);

    tau_max = MIN(tau_max,pma->tau0);
    pma->tw_max[index_wd]=tau_max;

    if(pma->matter_verbose > MATTER_VERBOSITY_RANGES){
      printf(" -> Obtained limits %f to %f for window %i, \n \
              with mean = (z: %f, chi: %f) and width (z: %f, chi: %f) \n",tau_min,tau_max,index_wd,ppt->selection_mean[index_wd],pma->tau0-(tau_max+tau_min)*0.5,ppt->selection_width[index_wd],(tau_max-tau_min)*0.5/ppr->selection_cut_at_sigma);
    }

    /**
     * Now we compute weights for integration over
     * the window times (tw, tau window)
     *
     * The same we do for the integrated windows
     * */
    class_call(array_weights_gauss_limits(
                                            pma->tw_sampling+index_wd*pma->tw_size,
                                            pma->tw_weights+index_wd*pma->tw_size,
                                            tau_min,
                                            tau_max,
                                            pma->tw_size,
                                            gauss_type_trapezoid,
                                            pma->error_message
                                            ),
               pma->error_message,
               pma->error_message);

    if(pma->uses_intxi_logarithmic){
      double xmin = pma->small_log_offset;
      double xmax = pma->tau0-tau_min; //Is always <tau0
      if(pma->uses_lensing_reduction){
        double nu_reduction = 0.5*pma->bias-3;//pma->bias-4.0;
        class_test(nu_reduction > 0 ,
                   pma->error_message,
                   "The lensing reduction method is not applicable for such a large bias");
        /**
         * Assume fully reduced tilt for the lensing spectra
         * Then calculate how far their reach is.
         * The windows decay as xi^nu_real = xi^nu_reduction
         *
         * Then (xi/xi_min)^nu_reduction != window_epsilon
         * */
        double window_epsilon = 1e-6;
        double chi_epsilon = (pma->tau0-0.5*(tau_min+tau_max))*pow(window_epsilon,-1./nu_reduction);//*(1-(pma->tau0-tau_max)/(pma->tau0-tau_min));
        xmin = MAX(chi_epsilon,pma->small_log_offset);
      }
      class_call(array_weights_gauss_limits(
                                              pma->integrated_tw_sampling+index_wd*pma->integrated_tw_size,
                                              pma->integrated_tw_weights+index_wd*pma->integrated_tw_size,
                                              log(xmin),
                                              log(xmax),
                                              pma->integrated_tw_size,
                                              gauss_type_trapezoid,
                                              pma->error_message
                                              ),
                   pma->error_message,
                   pma->error_message);
      int index_tw;
      pma->exp_integrated_tw_sampling[index_wd*pma->integrated_tw_size] = xmin;
      for(index_tw=1;index_tw<pma->integrated_tw_size-1;++index_tw){
        /**
         * Take the exponential of the logarithmic sampling (obtaining the original sampling),
         *  and make sure it is smaller than tau0 even when rounding errors occur.
         * */
        pma->exp_integrated_tw_sampling[index_wd*pma->integrated_tw_size+index_tw] = exp(pma->integrated_tw_sampling[index_wd*pma->integrated_tw_size+index_tw]);
      }
      pma->exp_integrated_tw_sampling[index_wd*pma->integrated_tw_size+pma->integrated_tw_size-1] = xmax;
    }
    else{
      class_call(array_weights_gauss_limits(
                                        pma->integrated_tw_sampling+index_wd*pma->integrated_tw_size,
                                        pma->integrated_tw_weights+index_wd*pma->integrated_tw_size,
                                        tau_min,
                                        pma->tau0,//285-pma->small_log_offset,
                                        pma->integrated_tw_size,
                                        gauss_type_trapezoid,
                                        pma->error_message
                                        ),
           pma->error_message,
           pma->error_message);
    }
    /**
     * Now we compute weights for integration over
     * the window times (tw, tau window)
     *
     * The same we do for the integrated windows
     * */
    class_call(array_weights_gauss_limits(
                                            pma->ptw_sampling+index_wd*pma->ptw_size,
                                            pma->ptw_weights+index_wd*pma->ptw_size,
                                            tau_min,
                                            tau_max,
                                            pma->ptw_size,
                                            gauss_type_trapezoid,
                                            pma->error_message
                                            ),
               pma->error_message,
               pma->error_message);
    class_call(array_weights_gauss_limits(
                                            pma->ptw_integrated_sampling+index_wd*pma->ptw_integrated_size,
                                            pma->ptw_integrated_weights+index_wd*pma->ptw_integrated_size,
                                            tau_min,
                                            pma->tau0,//285-pma->small_log_offset,
                                            pma->ptw_integrated_size,
                                            gauss_type_trapezoid,
                                            pma->error_message
                                            ),
               pma->error_message,
               pma->error_message);
  }
  return _SUCCESS_;
}
int matter_obtain_prepare_windows(
                          struct precision* ppr,
                          struct perturbs* ppt,
                          struct background* pba,
                          struct matters* pma
                          ){
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS){
    printf("Method :: Prepare windows including growth factor, derivatives, etc. \n");
  }
  /**
   * Initialize variables
   *  for the first part of this function
   * This part is essentially the same as in
   *  matter_obtain_windows
   * */
  int index_wd;
  int index_ptw;
  int index_ic;
  int index_radtp;
  int index_stp;
  int index_ptw_integrated;
  double window_at_z;
  double z;

  double dNdz;
  double dln_dNdz_dz;
  double prefactor;
  double* f_evo;

  /**
   * Allocate corresponding arrays
   * */
  class_alloc(f_evo,
              pma->num_windows*pma->ptw_size*sizeof(double),
              pma->error_message);
  class_alloc(pma->ptw_orig_window,
              pma->ptw_size*pma->num_windows*sizeof(double),
              pma->error_message);
  class_alloc(pma->ptw_window,
              pma->radtp_size_total*pma->ic_size*sizeof(double*),
              pma->error_message);
  class_alloc(pma->ptw_dwindow,
              pma->radtp_size_total*pma->ic_size*sizeof(double*),
              pma->error_message);
  class_alloc(pma->ptw_ddwindow,
              pma->radtp_size_total*pma->ic_size*sizeof(double*),
              pma->error_message);
  for(index_ic=0;index_ic<pma->ic_size;++index_ic){
    for(index_radtp=0;index_radtp<pma->radtp_size_total;++index_radtp){
      int size_local = pma->ptw_size;
      /**
       * If we want to integrate the window function,
       *  we have to use the integrated size
       * */
       //???
      if(matter_is_integrated(index_radtp)){
        size_local = MAX(size_local,pma->ptw_integrated_size);
      }
      class_alloc(pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp],
                  size_local*pma->num_windows*sizeof(double),
                  pma->error_message);
      class_alloc(pma->ptw_dwindow[index_ic*pma->radtp_size_total+index_radtp],
                  size_local*pma->num_windows*sizeof(double),
                  pma->error_message);
      class_alloc(pma->ptw_ddwindow[index_ic*pma->radtp_size_total+index_radtp],
                  size_local*pma->num_windows*sizeof(double),
                  pma->error_message);
    }
    //End radtp
  }
  //End ic

  /**
   * We set the variable last_index to a reasonable value
   *  using inter_normal,
   * so that later we can use inter_closeby without initial
   *  slowdown
   * This variable we save for the second part of the function
   * */
  int last_index = 0;
  int last_index_initial = 0;
  if(pma->ptw_size>0){
    class_call(background_at_tau(pba,
                                 pma->tau0,
                                 pba->short_info,
                                 pba->inter_normal,
                                 &last_index_initial,
                                 pma->short_pvecback),
                   pba->error_message,
                   pma->error_message);
    last_index = last_index_initial;
  }
  /**
   * Obtain the original window functions
   *  (The other ones are going to be multiplied with the
   *   growth factors, derived, etc.)
   * */
  for(index_wd=0;index_wd<pma->num_windows;++index_wd){
    for(index_ptw=0;index_ptw<pma->ptw_size;index_ptw++){
      /**
       * We find the z corresponding to the tau,
       *  then evaluate W(z),
       *  which is normalized according to
       *  int W(z) dz = 1
       * We then multiply it with H(z) to find
       *  W(tau) , which is normalized according to
       *  int W(tau) dtau = 1
       * Since this is not necessarily exactly 1,
       *  we resum the window values such that the integral is exactly 1
       *
       * Without re-normalization the error would still
       *  be tiny. Since this takes up really no time at all,
       *  we still do this tiny error correction
       * */
      class_call(background_at_tau(pba,
                                   pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw],
                                   pba->short_info,
                                   pba->inter_closeby,
                                   &last_index,
                                   pma->short_pvecback),
                 pba->error_message,
                 pma->error_message);

      z = pba->a_today/pma->short_pvecback[pba->index_bg_a]-1.;

      class_call(matter_window_function(
                                        ppr,
                                        ppt,
                                        pma,
                                        index_wd,
                                        z,
                                        &window_at_z),
                 pma->error_message,
                 pma->error_message);

      pma->ptw_orig_window[index_wd*pma->ptw_size+index_ptw] = window_at_z*pma->short_pvecback[pba->index_bg_H]*pba->a_today;
    }
    /**
     * Here we re-normalize the window, even though we will find norm ~= 1.0
     *  This statement has a tiny error
     * */
    double norm = 0.0;
    for (index_ptw = 0; index_ptw < pma->ptw_size; index_ptw++) {
     norm+=pma->ptw_orig_window[index_wd*pma->ptw_size+index_ptw]*pma->ptw_weights[index_wd*pma->ptw_size+index_ptw];
    }
    /**
     * The final step of re-normalization is setting
     *  W(tau) = H(z)*W(z)/ (int H(z) W(z) dtau)
     * It has the property that
     *  int W(tau) dtau = 1
     *  (which can be seen obviously)
     * */
    for (index_ptw = 0; index_ptw < pma->ptw_size; index_ptw++) {
      pma->ptw_orig_window[index_wd*pma->ptw_size+index_ptw]/=norm;
    }
  }
  /**
   * Initialize variables for second part of
   *  this function:
   * Now we want to multiply each window function
   *  by the growth factor, derivatives,
   * and if we have defined any type combinations,
   *  we want to combine the window functions
   *  of those type combinations
   * */
  int inf = 0;
  double a=0.0,b=0.0,h=0.0;
  double window_val;
  double tau;
  double* rsd_combined_windows   = NULL;
  double* rsd_combined_dwindows  = NULL;
  double* rsd_combined_ddwindows = NULL;
  double prefactor_dop1;
  double prefactor_dop2;
  double prefactor_rsd;
  /**
   * If we have doppler or gravitational terms,
   *  we are going to need the f_evo
   * */
  if(pma->has_doppler_terms || pma->has_gravitational_terms){
    /**
     * If there is a way to obtain it, obtain it here
     * */
    if((pma->has_nz_evo_file == _TRUE_) || (pma->has_nz_evo_analytic == _TRUE_)){
      /**
       * For every window and window sampling,
       *  iterate through all possible time sampling steps,
       *  get the background at this step,
       *  and finally
       * Here we can reuse the 'reasonable' index we found above
       * */
      last_index = last_index_initial;
      for(index_wd=0;index_wd<pma->num_windows;++index_wd){
        for(index_ptw=0;index_ptw<pma->ptw_size;++index_ptw){
          tau = pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw];
          class_call(background_at_tau(pba,
                                   tau,
                                   pba->short_info,
                                   pba->inter_closeby,
                                   &last_index,
                                   pma->short_pvecback),
                 pba->error_message,
                 pma->error_message);
          /**
           * First include the analytic terms
           * */
          f_evo[index_wd*pma->ptw_size+index_ptw] = 2./pma->short_pvecback[pba->index_bg_H]/pma->short_pvecback[pba->index_bg_a]/(pma->tau0-tau)
            + pma->short_pvecback[pba->index_bg_H_prime]/pma->short_pvecback[pba->index_bg_H]/pma->short_pvecback[pba->index_bg_H]/pma->short_pvecback[pba->index_bg_a];

          z = pba->a_today/pma->short_pvecback[pba->index_bg_a]-1.;

          /**
           * Then add the non-analytic terms
           * */
          if (pma->has_nz_evo_file ==_TRUE_) {

            class_test((z<pma->nz_evo_z[0]) || (z>pma->nz_evo_z[pma->nz_evo_size-1]),
                       pma->error_message,
                       "Your input file for the selection function only covers the redshift range [%f : %f]. However, your input for the selection function requires z=%f",
                       pma->nz_evo_z[0],
                       pma->nz_evo_z[pma->nz_evo_size-1],
                       z);


            class_call(array_interpolate_spline(
                                                pma->nz_evo_z,
                                                pma->nz_evo_size,
                                                pma->nz_evo_dlog_nz,
                                                pma->nz_evo_dd_dlog_nz,
                                                1,
                                                z,
                                                &last_index,
                                                &dln_dNdz_dz,
                                                1,
                                                pma->error_message),
                       pma->error_message,
                       pma->error_message);

          }
          else {
            /**
             * Or calculate them over some analytic approximation
             * */
            class_call(matter_dNdz_analytic(pma,
                                            z,
                                            &dNdz,
                                            &dln_dNdz_dz),
                       pma->error_message,
                       pma->error_message);

          }
          f_evo[index_wd*pma->ptw_size+index_ptw] -= dln_dNdz_dz/pma->short_pvecback[pba->index_bg_a];
        }
        //End ptw
      }
      //End wd
    }
    else {
      /**
       * If there is no way to obtain f_evo, set it to 0
       * */
      for(index_wd=0;index_wd<pma->num_windows;++index_wd){
        for(index_ptw=0;index_ptw<pma->ptw_size;++index_ptw){
          f_evo[index_wd*pma->ptw_size+index_ptw] = 0.;
        }
        //End ptw
      }
      //End wd
    }
    //Ifend obtainability
  }
  else {
    /**
     * If f_evo is just not required at all, also set it to 0
     * */
    for(index_wd=0;index_wd<pma->num_windows;++index_wd){
      for(index_ptw=0;index_ptw<pma->ptw_size;++index_ptw){
        f_evo[index_wd*pma->ptw_size+index_ptw] = 0.;
      }
      //End ptw
    }
    //End wd
  }
  //Ifend f_evo required
  /**
   * Now after aquiring the f_evo,
   *  we can do the full methodology
   * */
  for(index_ic=0;index_ic<pma->ic_size;++index_ic){
    for(index_radtp=0;index_radtp<pma->radtp_size_total;++index_radtp){
      index_stp = pma->index_stp_of_radtp[index_radtp];
      if(pma->matter_verbose > 4){
        printf(" -> Obtaining window at radial type %i with stp %i \n",index_radtp,index_stp);
      }
      /**
       * If the window is integrated, treat is specially
       * */
      if(
        matter_is_integrated(index_radtp)
      ){
        for(index_wd=0;index_wd<pma->num_windows;++index_wd){
          class_call(matter_integrate_window_function(pba,
                                           pma,
                                           pma->ptw_orig_window+index_wd*pma->ptw_size,
                                           pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_integrated_size,
                                           pma->ptw_integrated_sampling+index_wd*pma->ptw_integrated_size,
                                           pma->ptw_sampling+index_wd*pma->ptw_size,
                                           pma->ptw_weights+index_wd*pma->ptw_size,
                                           pma->ptw_integrated_size,
                                           pma->ptw_size,
                                           index_wd,
                                           index_radtp,
                                           f_evo
                                          ),
                      pma->error_message,
                      pma->error_message);
          /**
           * Multiply the final function with the growth factor
           * */
          class_call(matter_spline_prepare_hunt(pma->tau_sampling,
                                                pma->tau_size,
                                                pma->ptw_integrated_sampling[index_wd*pma->ptw_integrated_size],
                                                &inf,
                                                pma->error_message),
                     pma->error_message,
                     pma->error_message);
          for(index_ptw_integrated=0;index_ptw_integrated<pma->ptw_integrated_size;++index_ptw_integrated){
            class_call(matter_spline_hunt(
                                  pma->tau_sampling,
                                  pma->tau_size,
                                  pma->ptw_integrated_sampling[index_wd*pma->ptw_integrated_size+index_ptw_integrated],
                                  &inf,
                                  &h,
                                  &a,
                                  &b,
                                  pma->error_message
                                 ),
                pma->error_message,
                pma->error_message);

            pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_integrated_size+index_ptw_integrated]*=
            a * pma->growth_factor_tau[index_ic*pma->stp_size+index_stp][inf] +
            b * pma->growth_factor_tau[index_ic*pma->stp_size+index_stp][inf+1] +
            ((a*a*a-a)* pma->ddgrowth_factor_tau[index_ic*pma->stp_size+index_stp][inf] +
             (b*b*b-b)* pma->ddgrowth_factor_tau[index_ic*pma->stp_size+index_stp][inf+1])*h*h/6.;
          }
          //End ptw_integrated
        }
        //End wd
      }else if(matter_is_index(index_radtp,pma->radtp_rsd_combined,!pma->uses_relative_factors && pma->has_redshift_space_distortion && pma->uses_rsd_combination)){
        /**
         * If we have the combined type of redshift space distortion terms,
         *  we want to add up the different combinations of redshift space distortion
         * */
        last_index = last_index_initial;
        class_alloc(rsd_combined_windows,
                    3*pma->ptw_size*sizeof(double),
                    pma->error_message);
        class_alloc(rsd_combined_dwindows,
                    2*pma->ptw_size*sizeof(double),
                    pma->error_message);
        class_alloc(rsd_combined_ddwindows,
                    pma->ptw_size*sizeof(double),
                    pma->error_message);
        for(index_wd=0;index_wd<pma->num_windows;++index_wd){
          class_call(matter_spline_prepare_hunt(pma->tau_sampling,
                                     pma->tau_size,
                                     pma->ptw_sampling[index_wd*pma->ptw_size],
                                     &inf,
                                     pma->error_message),
                     pma->error_message,
                     pma->error_message);
          for(index_ptw=0;index_ptw<pma->ptw_size;++index_ptw){
            /**
             * We get the window at the ptw sampling value
             * Also the background quantities and such
             * */
            window_val = pma->ptw_orig_window[index_wd*pma->ptw_size+index_ptw];
            tau = pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw];
            class_call(background_at_tau(pba,
                               pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw],
                               pba->short_info,
                               pba->inter_normal,
                               &last_index,
                               pma->short_pvecback),
                pma->error_message,
                pma->error_message);
            class_call(matter_spline_hunt(
                                  pma->tau_sampling,
                                  pma->tau_size,
                                  pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw],
                                  &inf,
                                  &h,
                                  &a,
                                  &b,
                                  pma->error_message
                                 ),
                pma->error_message,
                pma->error_message);
            window_val *=
            a * pma->growth_factor_tau[index_ic*pma->stp_size+index_stp][inf] +
            b * pma->growth_factor_tau[index_ic*pma->stp_size+index_stp][inf+1] +
            ((a*a*a-a)* pma->ddgrowth_factor_tau[index_ic*pma->stp_size+index_stp][inf] +
             (b*b*b-b)* pma->ddgrowth_factor_tau[index_ic*pma->stp_size+index_stp][inf+1])*h*h/6.;
            /**
             * The difference in windows is given by their relative prefactors
             *  (And taking more or less derivatives)
             * */
            prefactor_dop1 = (
                    1.0 +
                    pma->short_pvecback[pba->index_bg_H_prime]/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]*pma->short_pvecback[pba->index_bg_H]) +
                    (2.0-5.0*pma->selection_magnification_bias[index_wd])/(pma->tau0-tau)/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]) +
                    5*pma->selection_magnification_bias[index_wd] -
                    f_evo[index_wd*pma->ptw_size+index_ptw]
                  );
            prefactor_dop2 = (
                  (f_evo[index_wd*pma->ptw_size+index_ptw]-3.0)*pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]
                  );
            prefactor_rsd = (
                    1.0/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H])
                  );
            //Ordered by amount of derivatives
            rsd_combined_windows[index_ptw] = prefactor_rsd*window_val;
            rsd_combined_windows[index_ptw+pma->ptw_size] = prefactor_dop1*window_val;
            rsd_combined_windows[index_ptw+2*pma->ptw_size] = prefactor_dop2*window_val;
          }
          //End ptw
          /**
           * Here we take the derivatives of those window functions,
           *  we will require these later
           * */
          class_call(matter_derive(
                        pma->ptw_sampling+index_wd*pma->ptw_size,
                        rsd_combined_windows,
                        pma->ptw_size,
                        rsd_combined_dwindows,
                        pma->error_message
                     ),
                     pma->error_message,
                     pma->error_message);
          class_call(matter_derive(
                        pma->ptw_sampling+index_wd*pma->ptw_size,
                        rsd_combined_windows+pma->ptw_size,
                        pma->ptw_size,
                        rsd_combined_dwindows+pma->ptw_size,
                        pma->error_message
                     ),
                     pma->error_message,
                     pma->error_message);
          class_call(matter_derive(
                        pma->ptw_sampling+index_wd*pma->ptw_size,
                        rsd_combined_dwindows,
                        pma->ptw_size,
                        rsd_combined_ddwindows,
                        pma->error_message
                     ),
                     pma->error_message,
                     pma->error_message);
          /**
           * Now we add the correct derivatives together to get the final
           *  prefactor for the redshift space distortion source
           * ( S_theta(k,tau) )
           * */
          for(index_ptw=0;index_ptw<pma->ptw_size;++index_ptw){
            pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_size+index_ptw]
              = rsd_combined_windows[index_ptw+2*pma->ptw_size] + rsd_combined_dwindows[index_ptw+pma->ptw_size] + rsd_combined_ddwindows[index_ptw];
          }
          //End tw
          /**
           * Finally, we calculate the derivatives of this final function
           * */
          class_call(matter_derive(
                        pma->ptw_sampling+index_wd*pma->ptw_size,
                        pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_size,
                        pma->ptw_size,
                        pma->ptw_dwindow[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_size,
                        pma->error_message
                     ),
                     pma->error_message,
                     pma->error_message);
          class_call(matter_derive(
                        pma->ptw_sampling+index_wd*pma->ptw_size,
                        pma->ptw_dwindow[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_size,
                        pma->ptw_size,
                        pma->ptw_ddwindow[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_size,
                        pma->error_message
                     ),
                     pma->error_message,
                     pma->error_message);
        }
        //End wds
        free(rsd_combined_windows);
        free(rsd_combined_dwindows);
        free(rsd_combined_ddwindows);
      }
      else if(matter_is_index(index_radtp,pma->radtp_combined,pma->uses_relative_factors)){
        /**
         * For the total combinations,
         *  we add up all different types with their relative factors
         * Here, we first allocate the required temporary window functions
         * */
        double* combined_windows;
        int index_stp_combined;
        double* dens1_window = NULL;
        double* dens1_dwindow = NULL;
        double* dens1_ddwindow = NULL;
        class_alloc(combined_windows,
                    pma->ptw_size*sizeof(double),
                    pma->error_message);
        if(pma->has_redshift_space_distortion){
          class_alloc(rsd_combined_windows,
                      3*pma->ptw_size*sizeof(double),
                      pma->error_message);
          class_alloc(rsd_combined_dwindows,
                      2*pma->ptw_size*sizeof(double),
                      pma->error_message);
          class_alloc(rsd_combined_ddwindows,
                      pma->ptw_size*sizeof(double),
                      pma->error_message);
        }
        if(pma->has_stp_delta_m && pma->uses_density_splitting){
          class_alloc(dens1_window,
                      pma->ptw_size*sizeof(double),
                      pma->error_message);
          class_alloc(dens1_dwindow,
                      pma->ptw_size*sizeof(double),
                      pma->error_message);
          class_alloc(dens1_ddwindow,
                      pma->ptw_size*sizeof(double),
                      pma->error_message);
        }
        for(index_wd=0;index_wd<pma->num_windows;++index_wd){
          class_call(matter_spline_prepare_hunt(pma->tau_sampling,
                                     pma->tau_size,
                                     pma->ptw_sampling[index_wd*pma->ptw_size],
                                     &inf,
                                     pma->error_message),
                     pma->error_message,
                     pma->error_message);
          for(index_ptw=0;index_ptw<pma->ptw_size;++index_ptw){
            /**
             * At each point we get the growth factor and background values
             * */
            window_val = pma->ptw_orig_window[index_wd*pma->ptw_size+index_ptw];
            tau = pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw];
            class_call(background_at_tau(pba,
                               pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw],
                               pba->short_info,
                               pba->inter_normal,
                               &last_index,
                               pma->short_pvecback),
                pma->error_message,
                pma->error_message);
            class_call(matter_spline_hunt(
                                  pma->tau_sampling,
                                  pma->tau_size,
                                  pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw],
                                  &inf,
                                  &h,
                                  &a,
                                  &b,
                                  pma->error_message
                                 ),
                pma->error_message,
                pma->error_message);
            /**
             * The total window starts off as 0
             * */
            combined_windows[index_ptw]=0.0;
            for(index_stp_combined=0;index_stp_combined<pma->stp_size;++index_stp_combined){
              double window_val_cur = window_val;
              /**
               * We only want those types which are required,
               *  phi+psi only appears in integrated windows
               * Otherwise we multiply with the growth factor
               * */
              if(matter_is_index(index_stp_combined,pma->stp_index_phi_plus_psi,pma->has_stp_phi_plus_psi)){continue;}
              window_val_cur *=
              a * pma->growth_factor_tau[index_ic*pma->stp_size+index_stp_combined][inf] +
              b * pma->growth_factor_tau[index_ic*pma->stp_size+index_stp_combined][inf+1] +
              ((a*a*a-a)* pma->ddgrowth_factor_tau[index_ic*pma->stp_size+index_stp_combined][inf] +
               (b*b*b-b)* pma->ddgrowth_factor_tau[index_ic*pma->stp_size+index_stp_combined][inf+1])*h*h/6.;
              /**
               * If we have redshift space distortion,
               *  we save the result in a temporary array
               * This we derive later
               * */
              if(matter_is_index(index_stp_combined,pma->stp_index_theta_m,pma->has_redshift_space_distortion)){
                prefactor_dop1 = pma->relative_factors[index_stp_combined]*(
                        1.0 +
                        pma->short_pvecback[pba->index_bg_H_prime]/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]*pma->short_pvecback[pba->index_bg_H]) +
                        (2.0-5.0*pma->selection_magnification_bias[index_wd])/(pma->tau0-tau)/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]) +
                        5*pma->selection_magnification_bias[index_wd] -
                        f_evo[index_wd*pma->ptw_size+index_ptw]
                      );
                prefactor_dop2 = pma->relative_factors[index_stp_combined]*(
                      (f_evo[index_wd*pma->ptw_size+index_ptw]-3.0)*pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]
                      );
                prefactor_rsd = pma->relative_factors[index_stp_combined]*(
                        1.0/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H])
                      );
                rsd_combined_windows[index_ptw] = prefactor_rsd*window_val_cur;
                rsd_combined_windows[index_ptw+pma->ptw_size] = prefactor_dop1*window_val_cur;
                rsd_combined_windows[index_ptw+2*pma->ptw_size] = prefactor_dop2*window_val_cur;
              }
              /**
               * If we have density, we also save in temporary arrays,
               * if we use density splitting,
               * since then we also need to take derivatives
               * */
              else if(matter_is_index(index_stp_combined,pma->stp_index_delta_m,pma->has_stp_delta_m && pma->uses_density_splitting)){
                double prefactor_dens = pma->relative_factors[index_stp_combined]*pma->selection_bias[index_wd];
                dens1_window[index_ptw] = prefactor_dens*window_val_cur;
                class_test(pma->selection_bias[index_wd]!=1.0,
                            pma->error_message,
                            "Bias was not 1.0");
                class_test(pma->relative_factors[index_stp_combined]!=1.0,
                             pma->error_message,
                             "Relative factor was not 1.0");
              }
              /**
               * Otherwise we can directly add the contribution
               * */
              else{
                double prefactor = 0;
                if(matter_is_index(index_stp_combined,pma->stp_index_psi,pma->has_gravitational_terms)){
                  prefactor = pma->relative_factors[index_stp_combined]*(
                    2.0 + pma->short_pvecback[pba->index_bg_H_prime]/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]*pma->short_pvecback[pba->index_bg_H])
                    + (2.0-5.0*pma->selection_magnification_bias[index_wd])/pma->short_pvecback[pba->index_bg_a]/pma->short_pvecback[pba->index_bg_H]/(pma->tau0-tau)
                    +5.0*pma->selection_magnification_bias[index_wd]-f_evo[index_wd*pma->ptw_size+index_ptw]
                  );
                }
                else if(matter_is_index(index_stp_combined,pma->stp_index_phi,pma->has_gravitational_terms)){
                  prefactor = pma->relative_factors[index_stp_combined]*(-2.0+5.0*pma->selection_magnification_bias[index_wd]);
                }
                else if(matter_is_index(index_stp_combined,pma->stp_index_phi_prime,pma->has_gravitational_terms)){
                  prefactor = pma->relative_factors[index_stp_combined]/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]);
                }
                else{
                  class_stop(pma->error_message,
                             "Unrecognized radial type (in combined) %i",index_stp_combined);
                }
                combined_windows[index_ptw] += prefactor*window_val_cur;
              }
              //Ifend different types
            }
            //End stp combined
          }
          //End ptw
          /**
           * If we have RSD or density splitting and
           *  need to calculate further derivatives,
           *  do it here, and add them to the combined window
           * */
          if(pma->has_redshift_space_distortion){
            class_call(matter_derive(
                          pma->ptw_sampling+index_wd*pma->ptw_size,
                          rsd_combined_windows,
                          pma->ptw_size,
                          rsd_combined_dwindows,
                          pma->error_message
                       ),
                       pma->error_message,
                       pma->error_message);
            class_call(matter_derive(
                          pma->ptw_sampling+index_wd*pma->ptw_size,
                          rsd_combined_windows+pma->ptw_size,
                          pma->ptw_size,
                          rsd_combined_dwindows+pma->ptw_size,
                          pma->error_message
                       ),
                       pma->error_message,
                       pma->error_message);
            class_call(matter_derive(
                          pma->ptw_sampling+index_wd*pma->ptw_size,
                          rsd_combined_dwindows,
                          pma->ptw_size,
                          rsd_combined_ddwindows,
                          pma->error_message
                       ),
                       pma->error_message,
                       pma->error_message);
            for(index_ptw=0;index_ptw<pma->ptw_size;++index_ptw){
              combined_windows[index_ptw]
                += (rsd_combined_windows[index_ptw+2*pma->ptw_size] + rsd_combined_dwindows[index_ptw+pma->ptw_size] + rsd_combined_ddwindows[index_ptw]);
            }
            //End ptw
          }
          //Ifend RSD derivatives
          if(pma->has_stp_delta_m && pma->uses_density_splitting){
            class_call(matter_derive(
                          pma->ptw_sampling+index_wd*pma->ptw_size,
                          dens1_window,
                          pma->ptw_size,
                          dens1_dwindow,
                          pma->error_message
                       ),
                       pma->error_message,
                       pma->error_message);
            class_call(matter_derive(
                          pma->ptw_sampling+index_wd*pma->ptw_size,
                          dens1_dwindow,
                          pma->ptw_size,
                          dens1_ddwindow,
                          pma->error_message
                       ),
                       pma->error_message,
                       pma->error_message);
            for(index_ptw=0;index_ptw<pma->ptw_size;++index_ptw){
              double tau = pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw];
              combined_windows[index_ptw]
                += (-dens1_ddwindow[index_ptw]  -2.0/(pma->tau0-tau)*dens1_dwindow[index_ptw] -2.0/(pma->tau0-tau)/(pma->tau0-tau)* dens1_window[index_ptw]);
            }
            //End ptw
          }
          //Ifend Density derivatives
          /**
           * Finally we calculate the window functions and the derivatives
           *  we possibly require later
           * */
          for(index_ptw=0;index_ptw<pma->ptw_size;++index_ptw){
            pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_size+index_ptw] = combined_windows[index_ptw];
          }
          class_call(matter_derive(
                        pma->ptw_sampling+index_wd*pma->ptw_size,
                        pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_size,
                        pma->ptw_size,
                        pma->ptw_dwindow[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_size,
                        pma->error_message
                     ),
                     pma->error_message,
                     pma->error_message);
          class_call(matter_derive(
                          pma->ptw_sampling+index_wd*pma->ptw_size,
                          pma->ptw_dwindow[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_size,
                          pma->ptw_size,
                          pma->ptw_ddwindow[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_size,
                          pma->error_message
                       ),
                       pma->error_message,
                       pma->error_message);
        }
        //End wds
        /**
         * Now delete temporary arrays
         * */
        if(pma->has_redshift_space_distortion){
          free(rsd_combined_windows);
          free(rsd_combined_dwindows);
          free(rsd_combined_ddwindows);
        }
        if(pma->has_stp_delta_m && pma->uses_density_splitting){
          free(dens1_window);
        }
        free(combined_windows);
      }
      else{
        for(index_wd=0;index_wd<pma->num_windows;++index_wd){
          class_call(matter_spline_prepare_hunt(pma->tau_sampling,
                                     pma->tau_size,
                                     pma->ptw_sampling[index_wd*pma->ptw_size],
                                     &inf,
                                     pma->error_message),
                     pma->error_message,
                     pma->error_message);
          for(index_ptw=0;index_ptw<pma->ptw_size;++index_ptw){
            /**
             * At each point get the growth factor,
             *  background values and finally multiply with prefactors
             * */
            window_val = pma->ptw_orig_window[index_wd*pma->ptw_size+index_ptw];
            tau = pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw];
            class_call(background_at_tau(pba,
                               pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw],
                               pba->short_info,
                               pba->inter_normal,
                               &last_index,
                               pma->short_pvecback),
                pma->error_message,
                pma->error_message);
            class_call(matter_spline_hunt(
                                  pma->tau_sampling,
                                  pma->tau_size,
                                  pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw],
                                  &inf,
                                  &h,
                                  &a,
                                  &b,
                                  pma->error_message
                                 ),
                pma->error_message,
                pma->error_message);
            window_val *=
            a * pma->growth_factor_tau[index_ic*pma->stp_size+index_stp][inf] +
            b * pma->growth_factor_tau[index_ic*pma->stp_size+index_stp][inf+1] +
            ((a*a*a-a)* pma->ddgrowth_factor_tau[index_ic*pma->stp_size+index_stp][inf] +
             (b*b*b-b)* pma->ddgrowth_factor_tau[index_ic*pma->stp_size+index_stp][inf+1])*h*h/6.;
            /**
             * Now we get the correct prefactor for any type
             * */
            if(matter_is_index(index_radtp,pma->radtp_dens,pma->has_stp_delta_m && (!pma->uses_density_splitting))){
              prefactor = pma->selection_bias[index_wd];
              class_test(pma->selection_bias[index_wd]!=1.0,
                        pma->error_message,
                        "Bias was not 1.0");
            }
            else if(matter_is_index(index_radtp,pma->radtp_dens1,pma->has_stp_delta_m && pma->uses_density_splitting)){
              prefactor = pma->selection_bias[index_wd];
              class_test(pma->selection_bias[index_wd]!=1.0,
                        pma->error_message,
                        "Bias was not 1.0");
            }
            else if(matter_is_index(index_radtp,pma->radtp_dens2,pma->has_stp_delta_m && pma->uses_density_splitting)){
              prefactor = pma->selection_bias[index_wd]/(pma->tau0-pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw])/(pma->tau0-pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw]);
              class_test(pma->selection_bias[index_wd]!=1.0,
                        pma->error_message,
                        "Bias was not 1.0");
            }
            else if(matter_is_index(index_radtp,pma->radtp_dop1,pma->has_doppler_terms && (!pma->uses_rsd_combination))){
              prefactor = (
                    1.0 +
                    pma->short_pvecback[pba->index_bg_H_prime]/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]*pma->short_pvecback[pba->index_bg_H]) +
                    (2.0-5.0*pma->selection_magnification_bias[index_wd])/(pma->tau0-tau)/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]) +
                    5*pma->selection_magnification_bias[index_wd] -
                    f_evo[index_wd*pma->ptw_size+index_ptw]
                  );
            }
            else if(matter_is_index(index_radtp,pma->radtp_dop2,pma->has_doppler_terms && (!pma->uses_rsd_combination))){
              prefactor = (
                    (f_evo[index_wd*pma->ptw_size+index_ptw]-3.0)*pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]
                  );
            }
            else if(matter_is_index(index_radtp,pma->radtp_rsd,pma->has_redshift_space_distortion && (!pma->uses_rsd_combination))){
              prefactor = (
                    1.0/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H])
                  );
            }
            else if(matter_is_index(index_radtp,pma->radtp_g1,pma->has_gravitational_terms)){
              prefactor = (
                2.0 + pma->short_pvecback[pba->index_bg_H_prime]/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]*pma->short_pvecback[pba->index_bg_H])
                + (2.0-5.0*pma->selection_magnification_bias[index_wd])/pma->short_pvecback[pba->index_bg_a]/pma->short_pvecback[pba->index_bg_H]/(pma->tau0-tau)
                +5.0*pma->selection_magnification_bias[index_wd]-f_evo[index_wd*pma->ptw_size+index_ptw]
              );
            }
            else if(matter_is_index(index_radtp,pma->radtp_g2,pma->has_gravitational_terms)){
              prefactor = (-2.0+5.0*pma->selection_magnification_bias[index_wd]);
            }
            else if(matter_is_index(index_radtp,pma->radtp_g3,pma->has_gravitational_terms)){
              prefactor = 1.0/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]);
            }
            else{
              class_stop(pma->error_message,
                         "Unrecognized radial type");
            }
            //Ifend which prefactor
            pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_size+index_ptw] = prefactor*window_val;
          }
          //End ptw
          /**
           * Now we derive the window function.
           *  These values we need (e.g. when density splitting or RSD (without combination)
           * */
          class_call(matter_derive(
                        pma->ptw_sampling+index_wd*pma->ptw_size,
                        pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_size,
                        pma->ptw_size,
                        pma->ptw_dwindow[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_size,
                        pma->error_message
                     ),
                     pma->error_message,
                     pma->error_message);
          class_call(matter_derive(
                        pma->ptw_sampling+index_wd*pma->ptw_size,
                        pma->ptw_dwindow[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_size,
                        pma->ptw_size,
                        pma->ptw_ddwindow[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_size,
                        pma->error_message
                     ),
                     pma->error_message,
                     pma->error_message);
        }
        //End wds
      }
      //Ifend has window combinations
    }
    //End radtp
  }
  //End ic
  free(f_evo);
  return _SUCCESS_;
}

int matter_obtain_prepare_windows_parallel(
                          struct precision* ppr,
                          struct perturbs* ppt,
                          struct background* pba,
                          struct matters* pma
                          ){
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS){
    printf("Method :: Prepare windows including growth factor, derivatives, etc. (parallel) \n");
  }
  /**
   * Initialize variables
   *  for the first part of this function
   * This part is essentially the same as in
   *  matter_obtain_windows
   * */
  int index_wd;
  int index_ptw;
  int index_ic;
  int index_radtp;
  int index_stp;
  int index_ptw_integrated;
  double window_at_z;
  double z;

  double dNdz;
  double dln_dNdz_dz;
  double prefactor;
  double* f_evo;

  /**
   * Allocate corresponding arrays
   * */
  class_alloc(f_evo,
              pma->num_windows*pma->ptw_size*sizeof(double),
              pma->error_message);
  class_alloc(pma->ptw_orig_window,
              pma->ptw_size*pma->num_windows*sizeof(double),
              pma->error_message);
  class_alloc(pma->ptw_window,
              pma->radtp_size_total*pma->ic_size*sizeof(double*),
              pma->error_message);
  class_alloc(pma->ptw_dwindow,
              pma->radtp_size_total*pma->ic_size*sizeof(double*),
              pma->error_message);
  class_alloc(pma->ptw_ddwindow,
              pma->radtp_size_total*pma->ic_size*sizeof(double*),
              pma->error_message);
  for(index_ic=0;index_ic<pma->ic_size;++index_ic){
    for(index_radtp=0;index_radtp<pma->radtp_size_total;++index_radtp){
      int size_local = pma->ptw_size;
      /**
       * If we want to integrate the window function,
       *  we have to use the integrated size
       * */
       //???
      if(matter_is_integrated(index_radtp)){
        size_local = MAX(size_local,pma->ptw_integrated_size);
      }
      class_alloc(pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp],
                  size_local*pma->num_windows*sizeof(double),
                  pma->error_message);
      class_alloc(pma->ptw_dwindow[index_ic*pma->radtp_size_total+index_radtp],
                  size_local*pma->num_windows*sizeof(double),
                  pma->error_message);
      class_alloc(pma->ptw_ddwindow[index_ic*pma->radtp_size_total+index_radtp],
                  size_local*pma->num_windows*sizeof(double),
                  pma->error_message);
    }
    //End radtp
  }
  //End ic

  /**
   * We set the variable last_index to a reasonable value
   *  using inter_normal,
   * so that later we can use inter_closeby without initial
   *  slowdown
   * This variable we save for the second part of the function
   * */
  int last_index = 0;
  int last_index_initial = 0;
  if(pma->ptw_size>0){
    class_call(background_at_tau(pba,
                                 pma->tau0,
                                 pba->short_info,
                                 pba->inter_normal,
                                 &last_index_initial,
                                 pma->short_pvecback),
                   pba->error_message,
                   pma->error_message);
    last_index = last_index_initial;
  }
  /**
   * Obtain the original window functions
   *  (The other ones are going to be multiplied with the
   *   growth factors, derived, etc.)
   * */
  for(index_wd=0;index_wd<pma->num_windows;++index_wd){
    for(index_ptw=0;index_ptw<pma->ptw_size;index_ptw++){
      /**
       * We find the z corresponding to the tau,
       *  then evaluate W(z),
       *  which is normalized according to
       *  int W(z) dz = 1
       * We then multiply it with H(z) to find
       *  W(tau) , which is normalized according to
       *  int W(tau) dtau = 1
       * Since this is not necessarily exactly 1,
       *  we resum the window values such that the integral is exactly 1
       *
       * Without re-normalization the error would still
       *  be tiny. Since this takes up really no time at all,
       *  we still do this tiny error correction
       * */
      class_call(background_at_tau(pba,
                                   pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw],
                                   pba->short_info,
                                   pba->inter_closeby,
                                   &last_index,
                                   pma->short_pvecback),
                 pba->error_message,
                 pma->error_message);

      z = pba->a_today/pma->short_pvecback[pba->index_bg_a]-1.;

      class_call(matter_window_function(
                                        ppr,
                                        ppt,
                                        pma,
                                        index_wd,
                                        z,
                                        &window_at_z),
                 pma->error_message,
                 pma->error_message);

      pma->ptw_orig_window[index_wd*pma->ptw_size+index_ptw] = window_at_z*pma->short_pvecback[pba->index_bg_H]*pba->a_today;
    }
    /**
     * Here we re-normalize the window, even though we will find norm ~= 1.0
     *  This statement has a tiny error
     * */
    double norm = 0.0;
    for (index_ptw = 0; index_ptw < pma->ptw_size; index_ptw++) {
     norm+=pma->ptw_orig_window[index_wd*pma->ptw_size+index_ptw]*pma->ptw_weights[index_wd*pma->ptw_size+index_ptw];
    }
    /**
     * The final step of re-normalization is setting
     *  W(tau) = H(z)*W(z)/ (int H(z) W(z) dtau)
     * It has the property that
     *  int W(tau) dtau = 1
     *  (which can be seen obviously)
     * */
    for (index_ptw = 0; index_ptw < pma->ptw_size; index_ptw++) {
      pma->ptw_orig_window[index_wd*pma->ptw_size+index_ptw]/=norm;
    }
  }
  /**
   * Initialize variables for second part of
   *  this function:
   * Now we want to multiply each window function
   *  by the growth factor, derivatives,
   * and if we have defined any type combinations,
   *  we want to combine the window functions
   *  of those type combinations
   * */
  int inf = 0;
  double a=0.0,b=0.0,h=0.0;
  double window_val;
  double tau;
  double prefactor_dop1;
  double prefactor_dop2;
  double prefactor_rsd;
  /**
   * If we have doppler or gravitational terms,
   *  we are going to need the f_evo
   * */
  if(pma->has_doppler_terms || pma->has_gravitational_terms){
    /**
     * If there is a way to obtain it, obtain it here
     * */
    if((pma->has_nz_evo_file == _TRUE_) || (pma->has_nz_evo_analytic == _TRUE_)){
      /**
       * For every window and window sampling,
       *  iterate through all possible time sampling steps,
       *  get the background at this step,
       *  and finally
       * Here we can reuse the 'reasonable' index we found above
       * */
      last_index = last_index_initial;
      for(index_wd=0;index_wd<pma->num_windows;++index_wd){
        for(index_ptw=0;index_ptw<pma->ptw_size;++index_ptw){
          tau = pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw];
          class_call(background_at_tau(pba,
                                   tau,
                                   pba->short_info,
                                   pba->inter_closeby,
                                   &last_index,
                                   pma->short_pvecback),
                 pba->error_message,
                 pma->error_message);
          /**
           * First include the analytic terms
           * */
          f_evo[index_wd*pma->ptw_size+index_ptw] = 2./pma->short_pvecback[pba->index_bg_H]/pma->short_pvecback[pba->index_bg_a]/(pma->tau0-tau)
            + pma->short_pvecback[pba->index_bg_H_prime]/pma->short_pvecback[pba->index_bg_H]/pma->short_pvecback[pba->index_bg_H]/pma->short_pvecback[pba->index_bg_a];

          z = pba->a_today/pma->short_pvecback[pba->index_bg_a]-1.;

          /**
           * Then add the non-analytic terms
           * */
          if (pma->has_nz_evo_file ==_TRUE_) {

            class_test((z<pma->nz_evo_z[0]) || (z>pma->nz_evo_z[pma->nz_evo_size-1]),
                       pma->error_message,
                       "Your input file for the selection function only covers the redshift range [%f : %f]. However, your input for the selection function requires z=%f",
                       pma->nz_evo_z[0],
                       pma->nz_evo_z[pma->nz_evo_size-1],
                       z);


            class_call(array_interpolate_spline(
                                                pma->nz_evo_z,
                                                pma->nz_evo_size,
                                                pma->nz_evo_dlog_nz,
                                                pma->nz_evo_dd_dlog_nz,
                                                1,
                                                z,
                                                &last_index,
                                                &dln_dNdz_dz,
                                                1,
                                                pma->error_message),
                       pma->error_message,
                       pma->error_message);

          }
          else {
            /**
             * Or calculate them over some analytic approximation
             * */
            class_call(matter_dNdz_analytic(pma,
                                            z,
                                            &dNdz,
                                            &dln_dNdz_dz),
                       pma->error_message,
                       pma->error_message);

          }
          f_evo[index_wd*pma->ptw_size+index_ptw] -= dln_dNdz_dz/pma->short_pvecback[pba->index_bg_a];
        }
        //End ptw
      }
      //End wd
    }
    else {
      /**
       * If there is no way to obtain f_evo, set it to 0
       * */
      for(index_wd=0;index_wd<pma->num_windows;++index_wd){
        for(index_ptw=0;index_ptw<pma->ptw_size;++index_ptw){
          f_evo[index_wd*pma->ptw_size+index_ptw] = 0.;
        }
        //End ptw
      }
      //End wd
    }
    //Ifend obtainability
  }
  else {
    /**
     * If f_evo is just not required at all, also set it to 0
     * */
    for(index_wd=0;index_wd<pma->num_windows;++index_wd){
      for(index_ptw=0;index_ptw<pma->ptw_size;++index_ptw){
        f_evo[index_wd*pma->ptw_size+index_ptw] = 0.;
      }
      //End ptw
    }
    //End wd
  }
  //Ifend f_evo required
  /**
   * Now after aquiring the f_evo,
   *  we can do the full methodology
   * */
  double** rsd_combined_windows   = NULL;
  double** rsd_combined_dwindows  = NULL;
  double** rsd_combined_ddwindows = NULL;
  int N_threads;
  //#ifdef _USE_OMP_
  N_threads = omp_get_max_threads();
  //#else
  //N_threads = 1;
  //#endif
  int n_thread;
  if((!pma->uses_relative_factors && pma->has_redshift_space_distortion && pma->uses_rsd_combination)
  || (pma->has_redshift_space_distortion && pma->uses_relative_factors)){
    class_alloc(rsd_combined_windows,
                N_threads*sizeof(double*),
                pma->error_message);
    class_alloc(rsd_combined_dwindows,
                N_threads*sizeof(double*),
                pma->error_message);
    class_alloc(rsd_combined_ddwindows,
                N_threads*sizeof(double*),
                pma->error_message);
    for(n_thread=0;n_thread<N_threads;++n_thread){
      class_alloc(rsd_combined_windows[n_thread],
                  3*pma->ptw_size*sizeof(double),
                  pma->error_message);
      class_alloc(rsd_combined_dwindows[n_thread],
                  2*pma->ptw_size*sizeof(double),
                  pma->error_message);
      class_alloc(rsd_combined_ddwindows[n_thread],
                  pma->ptw_size*sizeof(double),
                  pma->error_message);
    }
  }
  double** combined_windows = NULL;
  double** dens1_window = NULL;
  double** dens1_dwindow = NULL;
  double** dens1_ddwindow = NULL;
  if(pma->uses_relative_factors){
    class_alloc(combined_windows,
                N_threads*sizeof(double*),
                pma->error_message);
    for(n_thread=0;n_thread<N_threads;++n_thread){
      class_alloc(combined_windows[n_thread],
                  pma->ptw_size,
                  pma->error_message);
    }
    if(pma->has_stp_delta_m && pma->uses_density_splitting){
      class_alloc(dens1_window,
                  N_threads*sizeof(double*),
                  pma->error_message);
      class_alloc(dens1_dwindow,
                  N_threads*sizeof(double*),
                  pma->error_message);
      class_alloc(dens1_ddwindow,
                  N_threads*sizeof(double*),
                  pma->error_message);
      for(n_thread=0;n_thread<N_threads;++n_thread){
        class_alloc(dens1_window[n_thread],
                    pma->ptw_size*sizeof(double),
                    pma->error_message);
        class_alloc(dens1_dwindow[n_thread],
                    pma->ptw_size*sizeof(double),
                    pma->error_message);
        class_alloc(dens1_ddwindow[n_thread],
                    pma->ptw_size*sizeof(double),
                    pma->error_message);
      }
    }
  }
  for(index_ic=0;index_ic<pma->ic_size;++index_ic){
    for(index_radtp=0;index_radtp<pma->radtp_size_total;++index_radtp){
      index_stp = pma->index_stp_of_radtp[index_radtp];
      if(pma->matter_verbose > 4){
        printf(" -> Obtaining window at radial type %i with stp %i \n",index_radtp,index_stp);
      }
      /**
       * If the window is integrated, treat is specially
       * */
      if(
        matter_is_integrated(index_radtp)
      ){
        int abort = _FALSE_;
        #pragma omp parallel for private(index_wd,index_ptw,index_ptw_integrated,window_at_z,z,prefactor,tau,a,b,h,inf,prefactor_dop1,prefactor_dop2,prefactor_rsd,window_val) firstprivate(index_ic,index_radtp,pma,index_stp,f_evo,rsd_combined_windows,rsd_combined_dwindows,rsd_combined_ddwindows,combined_windows,dens1_window,dens1_dwindow,dens1_ddwindow)
        for(index_wd=0;index_wd<pma->num_windows;++index_wd){
          class_call_parallel(matter_integrate_window_function(pba,
                                           pma,
                                           pma->ptw_orig_window+index_wd*pma->ptw_size,
                                           pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_integrated_size,
                                           pma->ptw_integrated_sampling+index_wd*pma->ptw_integrated_size,
                                           pma->ptw_sampling+index_wd*pma->ptw_size,
                                           pma->ptw_weights+index_wd*pma->ptw_size,
                                           pma->ptw_integrated_size,
                                           pma->ptw_size,
                                           index_wd,
                                           index_radtp,
                                           f_evo
                                          ),
                      pma->error_message,
                      pma->error_message);
          /**
           * Multiply the final function with the growth factor
           * */
          class_call_parallel(matter_spline_prepare_hunt(pma->tau_sampling,
                                                pma->tau_size,
                                                pma->ptw_integrated_sampling[index_wd*pma->ptw_integrated_size],
                                                &inf,
                                                pma->error_message),
                     pma->error_message,
                     pma->error_message);
          for(index_ptw_integrated=0;index_ptw_integrated<pma->ptw_integrated_size;++index_ptw_integrated){
            class_call_parallel(matter_spline_hunt(
                                  pma->tau_sampling,
                                  pma->tau_size,
                                  pma->ptw_integrated_sampling[index_wd*pma->ptw_integrated_size+index_ptw_integrated],
                                  &inf,
                                  &h,
                                  &a,
                                  &b,
                                  pma->error_message
                                 ),
                pma->error_message,
                pma->error_message);

            pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_integrated_size+index_ptw_integrated]*=
            a * pma->growth_factor_tau[index_ic*pma->stp_size+index_stp][inf] +
            b * pma->growth_factor_tau[index_ic*pma->stp_size+index_stp][inf+1] +
            ((a*a*a-a)* pma->ddgrowth_factor_tau[index_ic*pma->stp_size+index_stp][inf] +
             (b*b*b-b)* pma->ddgrowth_factor_tau[index_ic*pma->stp_size+index_stp][inf+1])*h*h/6.;
          }
          //End ptw_integrated
        }
        //End wd
        if(abort==_TRUE_){return _FAILURE_;}
      }else if(matter_is_index(index_radtp,pma->radtp_rsd_combined,!pma->uses_relative_factors && pma->has_redshift_space_distortion && pma->uses_rsd_combination)){
        /**
         * If we have the combined type of redshift space distortion terms,
         *  we want to add up the different combinations of redshift space distortion
         * */
        last_index = last_index_initial;
        int abort = _FALSE_;
        #pragma omp parallel for private(index_wd,index_ptw,index_ptw_integrated,window_at_z,z,prefactor,tau,a,b,h,inf,prefactor_dop1,prefactor_dop2,prefactor_rsd,window_val) firstprivate(index_ic,index_radtp,pma,index_stp,f_evo,rsd_combined_windows,rsd_combined_dwindows,rsd_combined_ddwindows,combined_windows,dens1_window,dens1_dwindow,dens1_ddwindow)
        for(index_wd=0;index_wd<pma->num_windows;++index_wd){
          int tid = omp_get_thread_num();
          class_call_parallel(matter_spline_prepare_hunt(pma->tau_sampling,
                                     pma->tau_size,
                                     pma->ptw_sampling[index_wd*pma->ptw_size],
                                     &inf,
                                     pma->error_message),
                     pma->error_message,
                     pma->error_message);
          for(index_ptw=0;index_ptw<pma->ptw_size;++index_ptw){
            /**
             * We get the window at the ptw sampling value
             * Also the background quantities and such
             * */
            window_val = pma->ptw_orig_window[index_wd*pma->ptw_size+index_ptw];
            tau = pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw];
            class_call_parallel(background_at_tau(pba,
                               pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw],
                               pba->short_info,
                               pba->inter_normal,
                               &last_index,
                               pma->short_pvecback),
                pma->error_message,
                pma->error_message);
            class_call_parallel(matter_spline_hunt(
                                  pma->tau_sampling,
                                  pma->tau_size,
                                  pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw],
                                  &inf,
                                  &h,
                                  &a,
                                  &b,
                                  pma->error_message
                                 ),
                pma->error_message,
                pma->error_message);
            window_val *=
            a * pma->growth_factor_tau[index_ic*pma->stp_size+index_stp][inf] +
            b * pma->growth_factor_tau[index_ic*pma->stp_size+index_stp][inf+1] +
            ((a*a*a-a)* pma->ddgrowth_factor_tau[index_ic*pma->stp_size+index_stp][inf] +
             (b*b*b-b)* pma->ddgrowth_factor_tau[index_ic*pma->stp_size+index_stp][inf+1])*h*h/6.;
            /**
             * The difference in windows is given by their relative prefactors
             *  (And taking more or less derivatives)
             * */
            prefactor_dop1 = (
                    1.0 +
                    pma->short_pvecback[pba->index_bg_H_prime]/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]*pma->short_pvecback[pba->index_bg_H]) +
                    (2.0-5.0*pma->selection_magnification_bias[index_wd])/(pma->tau0-tau)/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]) +
                    5*pma->selection_magnification_bias[index_wd] -
                    f_evo[index_wd*pma->ptw_size+index_ptw]
                  );
            prefactor_dop2 = (
                  (f_evo[index_wd*pma->ptw_size+index_ptw]-3.0)*pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]
                  );
            prefactor_rsd = (
                    1.0/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H])
                  );
            //Ordered by amount of derivatives
            rsd_combined_windows[tid][index_ptw] = prefactor_rsd*window_val;
            rsd_combined_windows[tid][index_ptw+pma->ptw_size] = prefactor_dop1*window_val;
            rsd_combined_windows[tid][index_ptw+2*pma->ptw_size] = prefactor_dop2*window_val;
          }
          //End ptw
          /**
           * Here we take the derivatives of those window functions,
           *  we will require these later
           * */
          class_call_parallel(matter_derive(
                        pma->ptw_sampling+index_wd*pma->ptw_size,
                        rsd_combined_windows[tid],
                        pma->ptw_size,
                        rsd_combined_dwindows[tid],
                        pma->error_message
                     ),
                     pma->error_message,
                     pma->error_message);
          class_call_parallel(matter_derive(
                        pma->ptw_sampling+index_wd*pma->ptw_size,
                        rsd_combined_windows[tid]+pma->ptw_size,
                        pma->ptw_size,
                        rsd_combined_dwindows[tid]+pma->ptw_size,
                        pma->error_message
                     ),
                     pma->error_message,
                     pma->error_message);
          class_call_parallel(matter_derive(
                        pma->ptw_sampling+index_wd*pma->ptw_size,
                        rsd_combined_dwindows[tid],
                        pma->ptw_size,
                        rsd_combined_ddwindows[tid],
                        pma->error_message
                     ),
                     pma->error_message,
                     pma->error_message);
          /**
           * Now we add the correct derivatives together to get the final
           *  prefactor for the redshift space distortion source
           * ( S_theta(k,tau) )
           * */
          for(index_ptw=0;index_ptw<pma->ptw_size;++index_ptw){
            pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_size+index_ptw]
              = rsd_combined_windows[tid][index_ptw+2*pma->ptw_size] + rsd_combined_dwindows[tid][index_ptw+pma->ptw_size] + rsd_combined_ddwindows[tid][index_ptw];
          }
          //End tw
          /**
           * Finally, we calculate the derivatives of this final function
           * */
          class_call_parallel(matter_derive(
                        pma->ptw_sampling+index_wd*pma->ptw_size,
                        pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_size,
                        pma->ptw_size,
                        pma->ptw_dwindow[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_size,
                        pma->error_message
                     ),
                     pma->error_message,
                     pma->error_message);
          class_call_parallel(matter_derive(
                        pma->ptw_sampling+index_wd*pma->ptw_size,
                        pma->ptw_dwindow[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_size,
                        pma->ptw_size,
                        pma->ptw_ddwindow[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_size,
                        pma->error_message
                     ),
                     pma->error_message,
                     pma->error_message);
        }
        //End wds
        if(abort==_TRUE_){return _FAILURE_;}
      }
      else if(matter_is_index(index_radtp,pma->radtp_combined,pma->uses_relative_factors)){
        /**
         * For the total combinations,
         *  we add up all different types with their relative factors
         * Here, we first allocate the required temporary window functions
         * */
        int index_stp_combined;
        int abort = _FALSE_;
        #pragma omp parallel for private(index_wd,index_ptw,index_ptw_integrated,window_at_z,z,prefactor,tau,a,b,h,inf,prefactor_dop1,prefactor_dop2,prefactor_rsd,window_val) firstprivate(index_ic,index_radtp,pma,index_stp,f_evo,rsd_combined_windows,rsd_combined_dwindows,rsd_combined_ddwindows,combined_windows,dens1_window,dens1_dwindow,dens1_ddwindow)
        for(index_wd=0;index_wd<pma->num_windows;++index_wd){
          int tid = omp_get_thread_num();
          class_call_parallel(matter_spline_prepare_hunt(pma->tau_sampling,
                                     pma->tau_size,
                                     pma->ptw_sampling[index_wd*pma->ptw_size],
                                     &inf,
                                     pma->error_message),
                     pma->error_message,
                     pma->error_message);
          for(index_ptw=0;index_ptw<pma->ptw_size;++index_ptw){
            /**
             * At each point we get the growth factor and background values
             * */
            window_val = pma->ptw_orig_window[index_wd*pma->ptw_size+index_ptw];
            tau = pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw];
            class_call_parallel(background_at_tau(pba,
                               pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw],
                               pba->short_info,
                               pba->inter_normal,
                               &last_index,
                               pma->short_pvecback),
                pma->error_message,
                pma->error_message);
            class_call_parallel(matter_spline_hunt(
                                  pma->tau_sampling,
                                  pma->tau_size,
                                  pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw],
                                  &inf,
                                  &h,
                                  &a,
                                  &b,
                                  pma->error_message
                                 ),
                pma->error_message,
                pma->error_message);
            /**
             * The total window starts off as 0
             * */
            combined_windows[tid][index_ptw]=0.0;
            for(index_stp_combined=0;index_stp_combined<pma->stp_size;++index_stp_combined){
              double window_val_cur = window_val;
              /**
               * We only want those types which are required,
               *  phi+psi only appears in integrated windows
               * Otherwise we multiply with the growth factor
               * */
              if(matter_is_index(index_stp_combined,pma->stp_index_phi_plus_psi,pma->has_stp_phi_plus_psi)){continue;}
              window_val_cur *=
              a * pma->growth_factor_tau[index_ic*pma->stp_size+index_stp_combined][inf] +
              b * pma->growth_factor_tau[index_ic*pma->stp_size+index_stp_combined][inf+1] +
              ((a*a*a-a)* pma->ddgrowth_factor_tau[index_ic*pma->stp_size+index_stp_combined][inf] +
               (b*b*b-b)* pma->ddgrowth_factor_tau[index_ic*pma->stp_size+index_stp_combined][inf+1])*h*h/6.;
              /**
               * If we have redshift space distortion,
               *  we save the result in a temporary array
               * This we derive later
               * */
              if(matter_is_index(index_stp_combined,pma->stp_index_theta_m,pma->has_redshift_space_distortion)){
                prefactor_dop1 = pma->relative_factors[index_stp_combined]*(
                        1.0 +
                        pma->short_pvecback[pba->index_bg_H_prime]/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]*pma->short_pvecback[pba->index_bg_H]) +
                        (2.0-5.0*pma->selection_magnification_bias[index_wd])/(pma->tau0-tau)/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]) +
                        5*pma->selection_magnification_bias[index_wd] -
                        f_evo[index_wd*pma->ptw_size+index_ptw]
                      );
                prefactor_dop2 = pma->relative_factors[index_stp_combined]*(
                      (f_evo[index_wd*pma->ptw_size+index_ptw]-3.0)*pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]
                      );
                prefactor_rsd = pma->relative_factors[index_stp_combined]*(
                        1.0/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H])
                      );
                rsd_combined_windows[tid][index_ptw] = prefactor_rsd*window_val_cur;
                rsd_combined_windows[tid][index_ptw+pma->ptw_size] = prefactor_dop1*window_val_cur;
                rsd_combined_windows[tid][index_ptw+2*pma->ptw_size] = prefactor_dop2*window_val_cur;
              }
              /**
               * If we have density, we also save in temporary arrays,
               * if we use density splitting,
               * since then we also need to take derivatives
               * */
              else if(matter_is_index(index_stp_combined,pma->stp_index_delta_m,pma->has_stp_delta_m && pma->uses_density_splitting)){
                double prefactor_dens = pma->relative_factors[index_stp_combined]*pma->selection_bias[index_wd];
                dens1_window[tid][index_ptw] = prefactor_dens*window_val_cur;
              }
              /**
               * Otherwise we can directly add the contribution
               * */
              else{
                double prefactor = 0;
                if(matter_is_index(index_stp_combined,pma->stp_index_psi,pma->has_gravitational_terms)){
                  prefactor = pma->relative_factors[index_stp_combined]*(
                    2.0 + pma->short_pvecback[pba->index_bg_H_prime]/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]*pma->short_pvecback[pba->index_bg_H])
                    + (2.0-5.0*pma->selection_magnification_bias[index_wd])/pma->short_pvecback[pba->index_bg_a]/pma->short_pvecback[pba->index_bg_H]/(pma->tau0-tau)
                    +5.0*pma->selection_magnification_bias[index_wd]-f_evo[index_wd*pma->ptw_size+index_ptw]
                  );
                }
                else if(matter_is_index(index_stp_combined,pma->stp_index_phi,pma->has_gravitational_terms)){
                  prefactor = pma->relative_factors[index_stp_combined]*(-2.0+5.0*pma->selection_magnification_bias[index_wd]);
                }
                else if(matter_is_index(index_stp_combined,pma->stp_index_phi_prime,pma->has_gravitational_terms)){
                  prefactor = pma->relative_factors[index_stp_combined]/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]);
                }
                else{
                  class_test_parallel(_FALSE_,pma->error_message,
                             "Unrecognized radial type (in combined) %i",index_stp_combined);
                }
                combined_windows[tid][index_ptw] += prefactor*window_val_cur;
              }
              //Ifend different types
            }
            //End stp combined
          }
          //End ptw
          /**
           * If we have RSD or density splitting and
           *  need to calculate further derivatives,
           *  do it here, and add them to the combined window
           * */
          if(pma->has_redshift_space_distortion){
            class_call_parallel(matter_derive(
                          pma->ptw_sampling+index_wd*pma->ptw_size,
                          rsd_combined_windows[tid],
                          pma->ptw_size,
                          rsd_combined_dwindows[tid],
                          pma->error_message
                       ),
                       pma->error_message,
                       pma->error_message);
            class_call_parallel(matter_derive(
                          pma->ptw_sampling+index_wd*pma->ptw_size,
                          rsd_combined_windows[tid]+pma->ptw_size,
                          pma->ptw_size,
                          rsd_combined_dwindows[tid]+pma->ptw_size,
                          pma->error_message
                       ),
                       pma->error_message,
                       pma->error_message);
            class_call_parallel(matter_derive(
                          pma->ptw_sampling+index_wd*pma->ptw_size,
                          rsd_combined_dwindows[tid],
                          pma->ptw_size,
                          rsd_combined_ddwindows[tid],
                          pma->error_message
                       ),
                       pma->error_message,
                       pma->error_message);
            for(index_ptw=0;index_ptw<pma->ptw_size;++index_ptw){
              combined_windows[tid][index_ptw]
                += (rsd_combined_windows[tid][index_ptw+2*pma->ptw_size] + rsd_combined_dwindows[tid][index_ptw+pma->ptw_size] + rsd_combined_ddwindows[tid][index_ptw]);
            }
            //End ptw
          }
          //Ifend RSD derivatives
          if(pma->has_stp_delta_m && pma->uses_density_splitting){
            class_call_parallel(matter_derive(
                          pma->ptw_sampling+index_wd*pma->ptw_size,
                          dens1_window[tid],
                          pma->ptw_size,
                          dens1_dwindow[tid],
                          pma->error_message
                       ),
                       pma->error_message,
                       pma->error_message);
            class_call_parallel(matter_derive(
                          pma->ptw_sampling+index_wd*pma->ptw_size,
                          dens1_dwindow[tid],
                          pma->ptw_size,
                          dens1_ddwindow[tid],
                          pma->error_message
                       ),
                       pma->error_message,
                       pma->error_message);
            for(index_ptw=0;index_ptw<pma->ptw_size;++index_ptw){
              double tau = pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw];
              combined_windows[tid][index_ptw]
                += (-dens1_ddwindow[tid][index_ptw]  -2.0/(pma->tau0-tau)*dens1_dwindow[tid][index_ptw] -2.0/(pma->tau0-tau)/(pma->tau0-tau)* dens1_window[tid][index_ptw]);
            }
            //End ptw
          }
          //Ifend Density derivatives
          /**
           * Finally we calculate the window functions and the derivatives
           *  we possibly require later
           * */
          for(index_ptw=0;index_ptw<pma->ptw_size;++index_ptw){
            pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_size+index_ptw] = combined_windows[tid][index_ptw];
          }
          class_call_parallel(matter_derive(
                        pma->ptw_sampling+index_wd*pma->ptw_size,
                        pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_size,
                        pma->ptw_size,
                        pma->ptw_dwindow[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_size,
                        pma->error_message
                     ),
                     pma->error_message,
                     pma->error_message);
          class_call_parallel(matter_derive(
                          pma->ptw_sampling+index_wd*pma->ptw_size,
                          pma->ptw_dwindow[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_size,
                          pma->ptw_size,
                          pma->ptw_ddwindow[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_size,
                          pma->error_message
                       ),
                       pma->error_message,
                       pma->error_message);
        }
        //End wds
        if(abort==_TRUE_){return _FAILURE_;}
      }
      else{
        int abort = _FALSE_;
        #pragma omp parallel for private(index_wd,index_ptw,index_ptw_integrated,window_at_z,z,prefactor,tau,a,b,h,inf,prefactor_dop1,prefactor_dop2,prefactor_rsd,window_val) firstprivate(index_ic,index_radtp,pma,index_stp,f_evo,rsd_combined_windows,rsd_combined_dwindows,rsd_combined_ddwindows,combined_windows,dens1_window,dens1_dwindow,dens1_ddwindow)
        for(index_wd=0;index_wd<pma->num_windows;++index_wd){
          prefactor = 0.0;
          class_call_parallel(matter_spline_prepare_hunt(pma->tau_sampling,
                                     pma->tau_size,
                                     pma->ptw_sampling[index_wd*pma->ptw_size],
                                     &inf,
                                     pma->error_message),
                     pma->error_message,
                     pma->error_message);
          for(index_ptw=0;index_ptw<pma->ptw_size;++index_ptw){
            /**
             * At each point get the growth factor,
             *  background values and finally multiply with prefactors
             * */
            window_val = pma->ptw_orig_window[index_wd*pma->ptw_size+index_ptw];
            tau = pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw];
            class_call_parallel(background_at_tau(pba,
                               pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw],
                               pba->short_info,
                               pba->inter_normal,
                               &last_index,
                               pma->short_pvecback),
                pma->error_message,
                pma->error_message);
            class_call_parallel(matter_spline_hunt(
                                  pma->tau_sampling,
                                  pma->tau_size,
                                  pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw],
                                  &inf,
                                  &h,
                                  &a,
                                  &b,
                                  pma->error_message
                                 ),
                pma->error_message,
                pma->error_message);
            window_val *=
            a * pma->growth_factor_tau[index_ic*pma->stp_size+index_stp][inf] +
            b * pma->growth_factor_tau[index_ic*pma->stp_size+index_stp][inf+1] +
            ((a*a*a-a)* pma->ddgrowth_factor_tau[index_ic*pma->stp_size+index_stp][inf] +
             (b*b*b-b)* pma->ddgrowth_factor_tau[index_ic*pma->stp_size+index_stp][inf+1])*h*h/6.;
            /**
             * Now we get the correct prefactor for any type
             * */
            if(matter_is_index(index_radtp,pma->radtp_dens,pma->has_stp_delta_m && (!pma->uses_density_splitting))){
              prefactor = pma->selection_bias[index_wd];
            }
            else if(matter_is_index(index_radtp,pma->radtp_dens1,pma->has_stp_delta_m && pma->uses_density_splitting)){
              prefactor = pma->selection_bias[index_wd];
            }
            else if(matter_is_index(index_radtp,pma->radtp_dens2,pma->has_stp_delta_m && pma->uses_density_splitting)){
              prefactor = pma->selection_bias[index_wd]/(pma->tau0-pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw])/(pma->tau0-pma->ptw_sampling[index_wd*pma->ptw_size+index_ptw]);
            }
            else if(matter_is_index(index_radtp,pma->radtp_dop1,pma->has_doppler_terms && (!pma->uses_rsd_combination))){
              prefactor = (
                    1.0 +
                    pma->short_pvecback[pba->index_bg_H_prime]/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]*pma->short_pvecback[pba->index_bg_H]) +
                    (2.0-5.0*pma->selection_magnification_bias[index_wd])/(pma->tau0-tau)/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]) +
                    5*pma->selection_magnification_bias[index_wd] -
                    f_evo[index_wd*pma->ptw_size+index_ptw]
                  );
            }
            else if(matter_is_index(index_radtp,pma->radtp_dop2,pma->has_doppler_terms && (!pma->uses_rsd_combination))){
              prefactor = (
                    (f_evo[index_wd*pma->ptw_size+index_ptw]-3.0)*pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]
                  );
            }
            else if(matter_is_index(index_radtp,pma->radtp_rsd,pma->has_redshift_space_distortion && (!pma->uses_rsd_combination))){
              prefactor = (
                    1.0/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H])
                  );
            }
            else if(matter_is_index(index_radtp,pma->radtp_g1,pma->has_gravitational_terms)){
              prefactor = (
                2.0 + pma->short_pvecback[pba->index_bg_H_prime]/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]*pma->short_pvecback[pba->index_bg_H])
                + (2.0-5.0*pma->selection_magnification_bias[index_wd])/pma->short_pvecback[pba->index_bg_a]/pma->short_pvecback[pba->index_bg_H]/(pma->tau0-tau)
                +5.0*pma->selection_magnification_bias[index_wd]-f_evo[index_wd*pma->ptw_size+index_ptw]
              );
            }
            else if(matter_is_index(index_radtp,pma->radtp_g2,pma->has_gravitational_terms)){
              prefactor = (-2.0+5.0*pma->selection_magnification_bias[index_wd]);
            }
            else if(matter_is_index(index_radtp,pma->radtp_g3,pma->has_gravitational_terms)){
              prefactor = 1.0/(pma->short_pvecback[pba->index_bg_a]*pma->short_pvecback[pba->index_bg_H]);
            }
            else{
              class_test_parallel(_FALSE_,
                         pma->error_message,
                         "Unrecognized radial type");
            }
            //Ifend which prefactor
            pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_size+index_ptw] = prefactor*window_val;
          }
          //End ptw
          /**
           * Now we derive the window function.
           *  These values we need (e.g. when density splitting or RSD (without combination)
           * */
          class_call_parallel(matter_derive(
                        pma->ptw_sampling+index_wd*pma->ptw_size,
                        pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_size,
                        pma->ptw_size,
                        pma->ptw_dwindow[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_size,
                        pma->error_message
                     ),
                     pma->error_message,
                     pma->error_message);
          class_call_parallel(matter_derive(
                        pma->ptw_sampling+index_wd*pma->ptw_size,
                        pma->ptw_dwindow[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_size,
                        pma->ptw_size,
                        pma->ptw_ddwindow[index_ic*pma->radtp_size_total+index_radtp]+index_wd*pma->ptw_size,
                        pma->error_message
                     ),
                     pma->error_message,
                     pma->error_message);
        }
        //End wds
        if(abort==_TRUE_){return _FAILURE_;}
      }
      //Ifend has window combinations
    }
    //End radtp
  }
  //End ic
  free(f_evo);
  if((!pma->uses_relative_factors && pma->has_redshift_space_distortion && pma->uses_rsd_combination)
  || (pma->has_redshift_space_distortion && pma->uses_relative_factors)){
    for(n_thread=0;n_thread<N_threads;++n_thread){
      free(rsd_combined_windows[n_thread]);
      free(rsd_combined_dwindows[n_thread]);
      free(rsd_combined_ddwindows[n_thread]);
    }
    free(rsd_combined_windows);
    free(rsd_combined_dwindows);
    free(rsd_combined_ddwindows);
  }
  if(pma->uses_relative_factors){
    for(n_thread=0;n_thread<N_threads;++n_thread){
      free(combined_windows[n_thread]);
      if(pma->has_stp_delta_m && pma->uses_density_splitting){
        free(dens1_window[n_thread]);
        free(dens1_dwindow[n_thread]);
        free(dens1_ddwindow[n_thread]);
      }
    }
    if(pma->has_stp_delta_m && pma->uses_density_splitting){
      free(dens1_window);
      free(dens1_dwindow);
      free(dens1_ddwindow);
    }
    free(combined_windows);
  }
  return _SUCCESS_;
}
int matter_get_prepared_window_at(
                          struct matters* pma,
                          double tau,
                          int index_ic,
                          int index_radtp,
                          int index_wd,
                          int* last,
                          int derivative_type,
                          double* win_val
                          ){
  double a,b,h;
  if(derivative_type>=0){
    if(tau<pma->ptw_sampling[index_wd*pma->ptw_size] || tau>pma->ptw_sampling[(index_wd+1)*pma->ptw_size-1]){
      *win_val = 0.0;
      return _SUCCESS_;
    }

    class_call(matter_spline_hunt(pma->ptw_sampling+index_wd*pma->ptw_size,
                                  pma->ptw_size,
                                  tau,
                                  last,
                                  &h,
                                  &a,
                                  &b,
                                  pma->error_message),
                 pma->error_message,
                 pma->error_message);
    int last_index = *last;

    if(derivative_type==0){
      *win_val = a*pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_size+last_index]
                +b*pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_size+last_index+1];
    }
    if(derivative_type==1){
      *win_val = a*pma->ptw_dwindow[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_size+last_index]
                +b*pma->ptw_dwindow[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_size+last_index+1];
    }
    if(derivative_type==2){
      *win_val = a*pma->ptw_ddwindow[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_size+last_index]
                +b*pma->ptw_ddwindow[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_size+last_index+1];
    }
    if(derivative_type==3){
      class_stop(pma->error_message,"Currently this is not implemented ...");
      /* *win_val = a*(
                    -pma->ptw_ddwindow[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_size+last_index]
                    -2.0/(pma->tau0-tau)*pma->ptw_dwindow[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_size+last_index]
                    +(l*(l+1.0)-2.0)*pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_size+last_index]/(pma->tau0-tau)/(pma->tau0-tau)
                    )
                +b*(
                    -pma->ptw_ddwindow[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_size+last_index+1]
                    -2.0/(pma->tau0-tau)*pma->ptw_dwindow[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_size+last_index+1]
                    +(l*(l+1.0)-2.0)*pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_size+last_index+1]/(pma->tau0-tau)/(pma->tau0-tau)
                  );*/
    }
    if(derivative_type==4){
      *win_val = a*(
                    -pma->ptw_ddwindow[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_size+last_index]
                    -2.0/(pma->tau0-tau)*pma->ptw_dwindow[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_size+last_index]
                    -2.0*pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_size+last_index]/(pma->tau0-tau)/(pma->tau0-tau)
                    )
                +b*(
                    -pma->ptw_ddwindow[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_size+last_index+1]
                    -2.0/(pma->tau0-tau)*pma->ptw_dwindow[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_size+last_index+1]
                    -2.0*pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_size+last_index+1]/(pma->tau0-tau)/(pma->tau0-tau)
                    );
    }
  }
  else if(derivative_type==-1){
     if(tau<pma->ptw_integrated_sampling[index_wd*pma->ptw_integrated_size] ||
        tau>pma->ptw_integrated_sampling[index_wd*pma->ptw_integrated_size+pma->ptw_integrated_size-1]){
      *win_val = 0.0;
      return _SUCCESS_;
    }

    class_call(matter_spline_hunt(pma->ptw_integrated_sampling+index_wd*pma->ptw_integrated_size,
                                  pma->ptw_integrated_size,
                                  tau,
                                  last,
                                  &h,
                                  &a,
                                  &b,
                                  pma->error_message),
                 pma->error_message,
                 pma->error_message);
    int last_index = *last;

    *win_val = a*pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_integrated_size+last_index]
              +b*pma->ptw_window[index_ic*pma->radtp_size_total+index_radtp][index_wd*pma->ptw_integrated_size+last_index+1];
  }
  return _SUCCESS_;
}


int matter_window_function(
                          struct precision * ppr,
                          struct perturbs * ppt,
                          struct matters* pma,
                          int bin,
                          double z,
                          double * window_at_z) {

  double x;
  double dNdz;
  double dln_dNdz_dz;
  int last_index;

  // trivial dirac case
  if (ppt->selection==dirac) {

    *window_at_z=1.;

    return _SUCCESS_;
  }

  // difference between z and the bin center (we can take the absolute
  // value as long as all selection functions are symmetric around x=0)
  x=fabs(z-ppt->selection_mean[bin]);

  /** gaussian case */
  if (ppt->selection==gaussian) {

    *window_at_z = exp(-0.5*pow(x/ppt->selection_width[bin],2))
      /ppt->selection_width[bin]/sqrt(2.*_PI_);
    if ((pma->has_nz_file == _TRUE_) || (pma->has_nz_analytic == _TRUE_)) {

      if (pma->has_nz_file == _TRUE_) {

        class_test((z<pma->nz_z[0]) || (z>pma->nz_z[pma->nz_size-1]),
                   pma->error_message,
                   "Your input file for the selection function only covers the redshift range [%f : %f]. However, your input for the selection function requires z=%f",
                   pma->nz_z[0],
                   pma->nz_z[pma->nz_size-1],
                   z);

        class_call(array_interpolate_spline(
                                            pma->nz_z,
                                            pma->nz_size,
                                            pma->nz_nz,
                                            pma->nz_ddnz,
                                            1,
                                            z,
                                            &last_index,
                                            &dNdz,
                                            1,
                                            pma->error_message),
                   pma->error_message,
                   pma->error_message);
      }
      else {

        class_call(matter_dNdz_analytic(pma,
                                          z,
                                          &dNdz,
                                          &dln_dNdz_dz),
                   pma->error_message,
                   pma->error_message);
      }

      *window_at_z *= dNdz;
    }

    return _SUCCESS_;
  }
  /** top-hat case, with smoothed edges. The problem with sharp edges
     is that the final result will be affected by random
     noise. Indeed, the values of k at which the transfer functions
     Delta_l(k) are sampled will never coincide with the actual edges
     of the true transfer function (computed with or even without the
     Limber approximation). Hence the integral Cl=\int dk
     Delta_l(k)**2 (...) will be imprecise and will fluctuate randomly
     with the resolution along k. With smooth edges, the problem is
     solved, and the final Cls become mildly dependent on the
     resolution along k. */

  if (ppt->selection==tophat) {

    /* selection function, centered on z=mean (i.e. on x=0), equal to
       one around x=0, with tanh step centered on x=width, of width
       delta x = 0.1*width
    */
    *window_at_z=(1.-tanh((x-ppt->selection_width[bin])/(ppr->selection_tophat_edge*ppt->selection_width[bin])))/2.;

    if ((pma->has_nz_file == _TRUE_) || (pma->has_nz_analytic == _TRUE_)) {

      if (pma->has_nz_file == _TRUE_) {

        class_call(array_interpolate_spline(
                                            pma->nz_z,
                                            pma->nz_size,
                                            pma->nz_nz,
                                            pma->nz_ddnz,
                                            1,
                                            z,
                                            &last_index,
                                            &dNdz,
                                            1,
                                            pma->error_message),
                   pma->error_message,
                   pma->error_message);
      }
      else {

        class_call(matter_dNdz_analytic(pma,
                                          z,
                                          &dNdz,
                                          &dln_dNdz_dz),
                   pma->error_message,
                   pma->error_message);
      }

      *window_at_z *= dNdz;
    }

    return _SUCCESS_;
  }

  // get here only if selection type was not recognized
  class_stop(pma->error_message,
             "choice of window function not recognized.");

  return _SUCCESS_;
}
int matter_dNdz_analytic(
                           struct matters * pma,
                           double z,
                           double * dNdz,
                           double * dln_dNdz_dz) {

  /**
   * You can implement your favorite analytic ansatz for the selection
   *  function here.
   *
   * A typical function for a photometric sample:
   *  dN/dz = (z/z0)^alpha exp[-(z/z0)^beta].
   *  Then: dln(dN/dz)/dz = (alpha - beta*(z/z0)^beta)/z.
   *
   * In principle, one is free to use a different
   *  ansatz for the selection function and the evolution function.
   *
   * Since the selection function uses only dN/dz, while the
   *  evolution uses only dln(dN/dz)/dz, it is possible to use
   *  different functions for dN/dz and dln(dN/dz)/dz
   * */

  double z0,alpha,beta;

  z0 = 0.55;
  alpha = 2.0;
  beta = 1.5;

  *dNdz = pow(z/z0,alpha) * exp(-pow(z/z0,beta));

  *dln_dNdz_dz = (alpha - pow(z/z0,beta)*beta)/z;

  return _SUCCESS_;

}
int matter_obtain_primordial_spectrum(
                                      struct perturbs * ppt,
                                      struct primordial * ppm,
                                      struct matters * pma,
                                      double ** prim_spec
                                      ){
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS) {
    printf("Method :: Obtaining primordial spectrum\n");
  }
  /**
   * Prepare initial variables
   * */
  double * temp;
  int index_coeff;
  int index_ic1_ic2;
  int index_md = ppt->index_md_scalars;
  if(pma->matter_verbose > MATTER_VERBOSITY_RANGES) {
    printf(" -> Allocating mode %i (scalars) with ic size %i \n",index_md,ppm->ic_ic_size[index_md]);
  }
  /**
   * Allocate temporary arrays and primordial spectrum
   * */
  class_alloc(temp,
              pma->ic_ic_size*sizeof(double),
              pma->error_message);
  for(index_ic1_ic2=0;index_ic1_ic2<pma->ic_ic_size;++index_ic1_ic2){
    class_alloc(prim_spec[index_ic1_ic2],
                pma->size_fft_input*sizeof(double),
                pma->error_message);
  }
  //End ic1_ic2
  /**
   * We can now get the primordial spectrum from the primordial module
   * */
  for(index_coeff=0;index_coeff<pma->size_fft_input;++index_coeff){
    class_call(primordial_spectrum_at_k(ppm,index_md,linear,pma->k_sampling[index_coeff],temp),
               ppm->error_message,
               pma->error_message);
    for(index_ic1_ic2=0;index_ic1_ic2<ppm->ic_ic_size[index_md];++index_ic1_ic2){
        prim_spec[index_ic1_ic2][index_coeff] = temp[index_ic1_ic2];
    }
    //End ic1_ic2
  }
  //End coeff
  free(temp);
  return _SUCCESS_;
}
int matter_obtain_indices(
                          struct primordial* ppm,
                          struct perturbs* ppt,
                          struct matters* pma
                          ){
  if(pma->matter_verbose> MATTER_VERBOSITY_FUNCTIONS){
    printf("Method :: Obtaining indices (source types and perturbation types) \n");
  }
  /**
   * First, copy flags and indices from other structs
   * */
  int index_ic1_ic2;
  class_alloc(pma->is_non_zero,
              pma->ic_ic_size*sizeof(int),
              pma->error_message);

  for(index_ic1_ic2=0;index_ic1_ic2<pma->ic_ic_size;++index_ic1_ic2){
    pma->is_non_zero[index_ic1_ic2] = ppm->is_non_zero[ppt->index_md_scalars][index_ic1_ic2];
  }

  /**
   * Setting flags of which types are included
   * */
  pma->has_cltp_nc = ppt->has_cl_number_count;
  pma->has_stp_delta_m = ppt->has_nc_density;
  pma->has_redshift_space_distortion = ppt->has_nc_rsd;
  pma->has_lensing_terms = ppt->has_nc_lens;
  pma->has_gravitational_terms = ppt->has_nc_gr;
  pma->has_doppler_terms = pma->has_redshift_space_distortion;

  pma->has_cltp_sh = ppt->has_cl_lensing_potential;
  pma->has_cl_shear = pma->has_cltp_sh;


  pma->has_integrated_windows = (pma->has_cl_shear || pma->has_lensing_terms || pma->has_gravitational_terms);
  pma->has_unintegrated_windows = (pma->has_cltp_nc && (pma->has_stp_delta_m || pma->has_redshift_space_distortion || pma->has_gravitational_terms || pma->has_doppler_terms));
  pma->has_stp_phi_plus_psi = (pma->has_lensing_terms || pma->has_gravitational_terms || pma->has_cl_shear);

  pma->has_window_differentiation = (pma->has_redshift_space_distortion || (pma->uses_density_splitting && pma->has_stp_delta_m));

  /**
   * Defining indices of source types
   *
   * Here we explicitly give the correspondence between the types
   *  used in the perturbations module and the matter module
   * */
  int stp_size_counter = 0;
  int radtp_size_counter = 0;
  class_define_index(pma->stp_index_delta_m,      pma->has_stp_delta_m,stp_size_counter,                                   1);
  class_define_index(pma->stp_index_theta_m,      pma->has_redshift_space_distortion,stp_size_counter,                     1);
  class_define_index(pma->stp_index_phi,          pma->has_gravitational_terms,stp_size_counter,                           1);
  class_define_index(pma->stp_index_phi_prime,    pma->has_gravitational_terms,stp_size_counter,                           1);
  class_define_index(pma->stp_index_psi,          pma->has_gravitational_terms,stp_size_counter,                           1);
  class_define_index(pma->stp_index_phi_plus_psi, pma->has_stp_phi_plus_psi,stp_size_counter,                              1);

  /**
   * Defining indices of radial types or radtps
   *
   * These correspond to a combination of
   *   1) the source type
   *   2) the bessel function type
   *
   * When different relations are used,
   *  it is possible to combine some of these
   *  or instead to split them up
   * */
  if(pma->uses_relative_factors){
    class_define_index(pma->radtp_combined,       pma->uses_relative_factors,radtp_size_counter,                           1);
    if(pma->uses_density_splitting){
      class_define_index(pma->radtp_dens2,        pma->has_stp_delta_m,radtp_size_counter,                                 1);
    }
    else{
      class_define_index(pma->radtp_dens,         pma->has_stp_delta_m,radtp_size_counter,                                 1);
    }
  }else{
    if(pma->uses_density_splitting){
      class_define_index(pma->radtp_dens1,        pma->has_stp_delta_m,radtp_size_counter,                                 1);
      class_define_index(pma->radtp_dens2,        pma->has_stp_delta_m,radtp_size_counter,                                 1);
    }else{
      class_define_index(pma->radtp_dens,         pma->has_stp_delta_m,radtp_size_counter,                                 1);
    }
    if(pma->uses_rsd_combination){
      class_define_index(pma->radtp_rsd_combined, pma->has_redshift_space_distortion,radtp_size_counter,                   1);
    }
    else{
      class_define_index(pma->radtp_rsd,          pma->has_redshift_space_distortion,radtp_size_counter,                   1);
      class_define_index(pma->radtp_dop1,         pma->has_redshift_space_distortion,radtp_size_counter,                   1);
      class_define_index(pma->radtp_dop2,         pma->has_redshift_space_distortion,radtp_size_counter,                   1);
    }
    class_define_index(pma->radtp_g1,             pma->has_gravitational_terms,radtp_size_counter,                         1);
    class_define_index(pma->radtp_g2,             pma->has_gravitational_terms,radtp_size_counter,                         1);
    class_define_index(pma->radtp_g3,             pma->has_gravitational_terms,radtp_size_counter,                         1);
  }
  class_define_index(pma->radtp_nclens,           pma->has_lensing_terms,radtp_size_counter,                               1);
  class_define_index(pma->radtp_shlens,           pma->has_cl_shear,radtp_size_counter,                                    1);
  class_define_index(pma->radtp_g4,               pma->has_gravitational_terms,radtp_size_counter,                         1);
  class_define_index(pma->radtp_g5,               pma->has_gravitational_terms,radtp_size_counter,                         1);


  pma->stp_size = stp_size_counter;
  pma->stp_grid_size = pma->stp_size*pma->stp_size;
  pma->radtp_size_total = radtp_size_counter;
  pma->radtp_grid_size = pma->radtp_size_total*pma->radtp_size_total;

  /**
   * Define how the perturbation types, source types, and radial types
   *  depend on each other
   * This has to be done after the size of each is known,
   *  which is why it is done separately
   * Also print out which types were found
   * */
  class_alloc(pma->index_perturb_tp_of_stp,
              pma->stp_size*sizeof(int),
              pma->error_message);
  class_alloc(pma->index_stp_of_radtp,
              pma->radtp_size_total*sizeof(int),
              pma->error_message);
  if(pma->has_stp_delta_m){
    if(ppt->has_source_delta_cb){
      pma->index_perturb_tp_of_stp[pma->stp_index_delta_m] = ppt->index_tp_delta_cb;
    }
    else{
      pma->index_perturb_tp_of_stp[pma->stp_index_delta_m] = ppt->index_tp_delta_m;
    }
    if(pma->matter_verbose > MATTER_VERBOSITY_INDICES){
      printf(" -> Found delta_m at %i \n",pma->stp_index_delta_m);
    }
    if(pma->uses_density_splitting){
      pma->index_stp_of_radtp[pma->radtp_dens1] = pma->stp_index_delta_m;
      if(pma->matter_verbose > MATTER_VERBOSITY_INDICES){
        printf(" -> Corrsepondence dens1 to delta_m at %i radtp to %i stp \n",pma->radtp_dens1,pma->stp_index_delta_m);
      }
      pma->index_stp_of_radtp[pma->radtp_dens2] = pma->stp_index_delta_m;
      if(pma->matter_verbose > MATTER_VERBOSITY_INDICES){
        printf(" -> Corrsepondence dens2 to delta_m at %i radtp to %i stp \n",pma->radtp_dens2,pma->stp_index_delta_m);
      }
    }
    else{
      pma->index_stp_of_radtp[pma->radtp_dens] = pma->stp_index_delta_m;
      if(pma->matter_verbose > MATTER_VERBOSITY_INDICES){
        printf(" -> Corrsepondence dens to delta_m at %i radtp to %i stp \n",pma->radtp_dens,pma->stp_index_delta_m);
      }
    }
  }
  if(pma->has_redshift_space_distortion){
    pma->index_perturb_tp_of_stp[pma->stp_index_theta_m] = ppt->index_tp_theta_m;
    if(pma->matter_verbose > MATTER_VERBOSITY_INDICES){
      printf(" -> Found theta_m at %i \n",pma->stp_index_theta_m);
    }
    if(pma->uses_rsd_combination){
      pma->index_stp_of_radtp[pma->radtp_rsd_combined] = pma->stp_index_theta_m;
      if(pma->matter_verbose > MATTER_VERBOSITY_INDICES){
        printf(" -> Corrsepondence rsd_combined to theta_m at %i radtp to %i stp \n",pma->radtp_rsd_combined,pma->stp_index_theta_m);
      }
    }
    else{
      pma->index_stp_of_radtp[pma->radtp_dop1] = pma->stp_index_theta_m;
      pma->index_stp_of_radtp[pma->radtp_dop2] = pma->stp_index_theta_m;
      pma->index_stp_of_radtp[pma->radtp_rsd] = pma->stp_index_theta_m;
      if(pma->matter_verbose > MATTER_VERBOSITY_INDICES){
        printf(" -> Corrsepondence dop1 to theta_m at %i radtp to %i stp \n",pma->radtp_dop1,pma->stp_index_theta_m);
        printf(" -> Corrsepondence dop2 to theta_m at %i radtp to %i stp \n",pma->radtp_dop2,pma->stp_index_theta_m);
        printf(" -> Corrsepondence rsd to theta_m at %i radtp to %i stp \n",pma->radtp_rsd,pma->stp_index_theta_m);
      }
    }
  }

  if(pma->has_gravitational_terms){
    pma->index_perturb_tp_of_stp[pma->stp_index_phi] = ppt->index_tp_phi;
    pma->index_perturb_tp_of_stp[pma->stp_index_phi_prime] = ppt->index_tp_phi_prime;
    pma->index_perturb_tp_of_stp[pma->stp_index_psi] = ppt->index_tp_psi;
    pma->index_stp_of_radtp[pma->radtp_g1] = pma->stp_index_psi;
    pma->index_stp_of_radtp[pma->radtp_g2] = pma->stp_index_phi;
    pma->index_stp_of_radtp[pma->radtp_g3] = pma->stp_index_phi_prime;
    pma->index_stp_of_radtp[pma->radtp_g5] = pma->stp_index_phi_prime;
    if(pma->matter_verbose > MATTER_VERBOSITY_INDICES){
      printf(" -> Found phi at %i \n",pma->stp_index_phi);
      printf(" -> Found phi' at %i \n",pma->stp_index_phi_prime);
      printf(" -> Found psi at %i \n",pma->stp_index_psi);
      printf(" -> Corrsepondence g1 to psi at %i radtp to %i stp \n",pma->radtp_g1,pma->stp_index_psi);
      printf(" -> Corrsepondence g2 to phi at %i radtp to %i stp \n",pma->radtp_g2,pma->stp_index_phi);
      printf(" -> Corrsepondence g3 to phi' at %i radtp to %i stp \n",pma->radtp_g3,pma->stp_index_phi_prime);
      printf(" -> Corrsepondence g5 to phi' at %i radtp to %i stp \n",pma->radtp_g5,pma->stp_index_phi_prime);
    }
  }
  if(pma->has_stp_phi_plus_psi){
    pma->index_perturb_tp_of_stp[pma->stp_index_phi_plus_psi] = ppt->index_tp_phi_plus_psi;
    if(pma->matter_verbose > MATTER_VERBOSITY_INDICES){
      printf(" -> Found phi+psi at %i ( perturb %i ) \n",pma->stp_index_phi_plus_psi,ppt->index_tp_phi_plus_psi);
    }
    if(pma->has_lensing_terms){
      pma->index_stp_of_radtp[pma->radtp_nclens] = pma->stp_index_phi_plus_psi;
      if(pma->matter_verbose > MATTER_VERBOSITY_INDICES){
        printf(" -> Corrsepondence nc lensing to phi_plus_psi at %i radtp to %i stp \n",pma->radtp_nclens,pma->stp_index_phi_plus_psi);
      }
    }
    if(pma->has_cl_shear){
      pma->index_stp_of_radtp[pma->radtp_shlens] = pma->stp_index_phi_plus_psi;
      if(pma->matter_verbose > MATTER_VERBOSITY_INDICES){
        printf(" -> Corrsepondence shear lensing to phi_plus_psi at %i radtp to %i stp \n",pma->radtp_shlens,pma->stp_index_phi_plus_psi);
      }
    }
    if(pma->has_gravitational_terms){
      pma->index_stp_of_radtp[pma->radtp_g4] = pma->stp_index_phi_plus_psi;
      if(pma->matter_verbose > MATTER_VERBOSITY_INDICES){
        printf(" -> Corrsepondence g4 to phi_plus_psi at %i radtp to %i stp \n",pma->radtp_g4,pma->stp_index_phi_plus_psi);
      }
    }
  }
  if(pma->matter_verbose > MATTER_VERBOSITY_INDICES){
    printf("Total number of source types for matter : %i \n",pma->stp_size);
    printf("Total number of radial types for matter : %i \n",pma->radtp_size_total);
  }
  /**
   * If we try to reexpress the c_n^XY through each other,
   *  using a relative factor,
   * we here select what function we derive the factor relatively to
   * */
  if(pma->uses_relative_factors){
    if(pma->has_stp_delta_m){
      pma->index_relative_stp = pma->stp_index_delta_m;
    }
    else if(pma->has_redshift_space_distortion){
      pma->index_relative_stp = pma->stp_index_theta_m;
    }
    else if(pma->has_lensing_terms || pma->has_cl_shear){
      pma->index_relative_stp = pma->stp_index_phi_plus_psi;
    }
    else if(pma->has_gravitational_terms){
      pma->index_relative_stp = pma->stp_index_phi_plus_psi;
    }
    else{
      class_stop(pma->error_message,
                 "No type to take relative factors of");
    }
    pma->index_stp_of_radtp[pma->radtp_combined] = pma->index_relative_stp;
  }
  /**
   * Finally we want to count the total number of
   * Cl - types like number-count Cl's (nCl/dCl) or shear Cl's (sCl)
   * Currently only nCl and sCl are supported
   * */
  int cltp_size_counter = 0;
  class_define_index(pma->cltp_index_nc,pma->has_cltp_nc,cltp_size_counter,                 1);
  class_define_index(pma->cltp_index_sh,pma->has_cltp_sh,cltp_size_counter,                 1);
  pma->cltp_size = cltp_size_counter;

  if(pma->matter_verbose>MATTER_VERBOSITY_INDICES){
    printf(" -> Found requested Cl's :\n");
    if(pma->has_cltp_nc){
      printf("    -> Number Count Cl's \n");
    }
    if(pma->has_cltp_sh){
      printf("    -> Shear Cl's \n");
    }
  }
  return _SUCCESS_;
}
int matter_obtain_bi_indices(
                              struct matters* pma
                              ){
  if(pma->matter_verbose >MATTER_VERBOSITY_FUNCTIONS ){
    printf("Method :: Obtaining bessel integral and tilt indices \n");
  }
  /**
   * Whereas before we found how the sources relate to each other
   *  (Cosmological part)
   * now we find how the bessel integrals and tilts are related
   *  (Geometrical part)
   * Of course sometimes we don't need to calculate all possible
   *  tilts of the bessel integrals, and this depends on not only
   *  our options, but also if the corresponding radial types are
   *  defined or not
   *
   * */
  int bi_size_counter   = 0;
  int tilt_size_counter = 0;
  int bi_normal_size_counter =0;
  int bi_reduced_size_counter=0;
  int bi_lfactor_size_nc_counter=0;
  int bi_lfactor_size_sh_counter=0;

  /**
   * Of course checking how many functions and which are really required
   *  does get a bit tedious at times
   * */
  pma->has_bitp_normal = _FALSE_;
  pma->has_bitp_nu_reduced = _FALSE_;
  pma->has_bitp_lfactor = _FALSE_;
  if(pma->uses_relative_factors){
    if(pma->has_stp_delta_m || pma->has_redshift_space_distortion || pma->has_gravitational_terms || pma->has_lensing_terms){
      bi_reduced_size_counter=1;
    }
  }
  if(pma->has_stp_delta_m){
    if(pma->uses_density_splitting){
      pma->has_bitp_nu_reduced = _TRUE_;
      pma->has_bitp_lfactor = _TRUE_;
      if(!pma->uses_relative_factors){
        bi_reduced_size_counter++;
      }
      bi_lfactor_size_nc_counter++;
    }
    else{
      pma->has_bitp_normal = _TRUE_;
      bi_normal_size_counter++;
    }
  }
  if(pma->has_redshift_space_distortion){
    if(!pma->uses_relative_factors){
      if(pma->uses_rsd_combination){
        bi_reduced_size_counter++;
      }
      else{
        bi_reduced_size_counter+=3;
      }
    }
    pma->has_bitp_nu_reduced = _TRUE_;
  }
  if(pma->has_lensing_terms){
    bi_lfactor_size_nc_counter++;
    pma->has_bitp_lfactor = _TRUE_;
  }
  if(pma->has_cl_shear){
    bi_lfactor_size_sh_counter++;
    pma->has_bitp_lfactor = _TRUE_;
  }
  if(pma->has_gravitational_terms){
    if(!pma->uses_relative_factors){
      bi_reduced_size_counter+=5;
    }
    else{
      bi_reduced_size_counter+=2;
    }
    pma->has_bitp_nu_reduced= _TRUE_;
  }
  pma->has_tilt_normal = pma->has_bitp_normal;
  pma->has_tilt_reduced = (pma->has_bitp_lfactor || pma->has_bitp_nu_reduced);

  /**
   * Once we have find which tilts should exist
   *  and which bessel integral types should exist,
   *  we now define corresponding indices
   * */
  // Normal bessel functions
  class_define_index(pma->bitp_index_normal, pma->has_bitp_normal,bi_size_counter,           1);
  class_define_index(pma->tilt_index_normal, pma->has_bitp_normal,tilt_size_counter,         1);
  // Nu-2 and Nu-4 bessel functions
  class_define_index(pma->bitp_index_nu_reduced, pma->has_bitp_nu_reduced ,bi_size_counter,  1);
  class_define_index(pma->tilt_index_reduced, pma->has_bitp_nu_reduced || pma->has_bitp_lfactor,tilt_size_counter,    1);
  // l(l+1) prefactor of bessel functions
  // (does not introduce new tilt)
  class_define_index(pma->bitp_index_lfactor, pma->has_bitp_lfactor,bi_size_counter,         1);

  pma->bitp_size = bi_size_counter;
  pma->tilt_size = tilt_size_counter;
  pma->tilt_grid_size = (pma->tilt_size*(pma->tilt_size+1))/2;


  /**
   * Now we once again build an analogy table,
   *  saying for each bessel integral types
   *  which radial types can be found.
   *
   * We also define a small macro that takes care of
   *  checking whether all desired indices have been correctly assigned
   * This is mostly a macro checking whether or not the function is
   *  written correctly is changed
   * */
  if(pma->matter_verbose>MATTER_VERBOSITY_INDICES){
    printf(" -> Analysis of radial and bessel integral type structure : \n");
    printf(" -> Found number of tilts %i (symmetric grid %i) \n", pma->tilt_size,pma->tilt_grid_size);
    printf(" -> Found number of bessel integral types %i \n",pma->bitp_size);
    printf("    -> Bitp normal is %sfound (%4d indices)\n",(pma->has_bitp_normal?"":"not "),bi_normal_size_counter);
    printf("    -> Bitp reduced is %sfound (%4d indices)\n",(pma->has_bitp_nu_reduced?"":"not "),bi_reduced_size_counter);
    printf("    -> Bitp lfactor is %sfound (%4d indices(nc only), %4d indices(sh only))\n",(pma->has_bitp_lfactor?"":"not "),bi_lfactor_size_nc_counter,bi_lfactor_size_sh_counter);
  }
  class_alloc(pma->radtps_of_bitp,
              pma->bitp_size*pma->cltp_size*sizeof(double*),
              pma->error_message);
  class_alloc(pma->radtp_of_bitp_size,
              pma->bitp_size*pma->cltp_size*sizeof(double),
              pma->error_message);

  #define matter_index_correspondence(store_array,condition,index_in_array_counter,original_index) \
  if((condition)){                                                                                 \
    (store_array)[--(index_in_array_counter)] = (original_index);                                  \
  }
  if(pma->has_bitp_normal){
    if(pma->has_cltp_nc){
      class_alloc(pma->radtps_of_bitp[pma->cltp_index_nc*pma->bitp_size+pma->bitp_index_normal],
                  bi_normal_size_counter*sizeof(double),
                  pma->error_message);
      pma->radtp_of_bitp_size[pma->cltp_index_nc*pma->bitp_size+pma->bitp_index_normal] = bi_normal_size_counter;
      matter_index_correspondence(pma->radtps_of_bitp[pma->cltp_index_nc*pma->bitp_size+pma->bitp_index_normal],(pma->has_stp_delta_m  && (!pma->uses_density_splitting)),bi_normal_size_counter,pma->radtp_dens)
      class_test(bi_normal_size_counter!=0,
                 pma->error_message,
                 "Number of radial types for bessel integral type 'normal' do not match up.");
    }
    if(pma->has_cltp_sh){
      pma->radtp_of_bitp_size[pma->cltp_index_sh*pma->bitp_size+pma->bitp_index_normal] = 0;
    }
  }
  if(pma->has_bitp_nu_reduced){
    if(pma->has_cltp_nc){
      class_alloc(pma->radtps_of_bitp[pma->cltp_index_nc*pma->bitp_size+pma->bitp_index_nu_reduced],
                  bi_reduced_size_counter*sizeof(double),
                  pma->error_message);
      pma->radtp_of_bitp_size[pma->cltp_index_nc*pma->bitp_size+pma->bitp_index_nu_reduced] = bi_reduced_size_counter;
      if(pma->uses_relative_factors){
        matter_index_correspondence(pma->radtps_of_bitp[pma->cltp_index_nc*pma->bitp_size+pma->bitp_index_nu_reduced],pma->uses_relative_factors,bi_reduced_size_counter,pma->radtp_combined)
        matter_index_correspondence(pma->radtps_of_bitp[pma->cltp_index_nc*pma->bitp_size+pma->bitp_index_nu_reduced],pma->has_gravitational_terms,bi_reduced_size_counter,pma->radtp_g4)
        matter_index_correspondence(pma->radtps_of_bitp[pma->cltp_index_nc*pma->bitp_size+pma->bitp_index_nu_reduced],pma->has_gravitational_terms,bi_reduced_size_counter,pma->radtp_g5)
      }
      else{
        if(pma->uses_rsd_combination){
          matter_index_correspondence(pma->radtps_of_bitp[pma->cltp_index_nc*pma->bitp_size+pma->bitp_index_nu_reduced],pma->has_redshift_space_distortion,bi_reduced_size_counter,pma->radtp_rsd_combined)
        }
        else{
          matter_index_correspondence(pma->radtps_of_bitp[pma->cltp_index_nc*pma->bitp_size+pma->bitp_index_nu_reduced],pma->has_redshift_space_distortion,bi_reduced_size_counter,pma->radtp_rsd)
          matter_index_correspondence(pma->radtps_of_bitp[pma->cltp_index_nc*pma->bitp_size+pma->bitp_index_nu_reduced],pma->has_redshift_space_distortion,bi_reduced_size_counter,pma->radtp_dop1)
          matter_index_correspondence(pma->radtps_of_bitp[pma->cltp_index_nc*pma->bitp_size+pma->bitp_index_nu_reduced],pma->has_redshift_space_distortion,bi_reduced_size_counter,pma->radtp_dop2)
        }
        matter_index_correspondence(pma->radtps_of_bitp[pma->cltp_index_nc*pma->bitp_size+pma->bitp_index_nu_reduced],(pma->has_stp_delta_m && pma->uses_density_splitting),bi_reduced_size_counter,pma->radtp_dens1)
        matter_index_correspondence(pma->radtps_of_bitp[pma->cltp_index_nc*pma->bitp_size+pma->bitp_index_nu_reduced],pma->has_gravitational_terms,bi_reduced_size_counter,pma->radtp_g1)
        matter_index_correspondence(pma->radtps_of_bitp[pma->cltp_index_nc*pma->bitp_size+pma->bitp_index_nu_reduced],pma->has_gravitational_terms,bi_reduced_size_counter,pma->radtp_g2)
        matter_index_correspondence(pma->radtps_of_bitp[pma->cltp_index_nc*pma->bitp_size+pma->bitp_index_nu_reduced],pma->has_gravitational_terms,bi_reduced_size_counter,pma->radtp_g3)
        matter_index_correspondence(pma->radtps_of_bitp[pma->cltp_index_nc*pma->bitp_size+pma->bitp_index_nu_reduced],pma->has_gravitational_terms,bi_reduced_size_counter,pma->radtp_g4)
        matter_index_correspondence(pma->radtps_of_bitp[pma->cltp_index_nc*pma->bitp_size+pma->bitp_index_nu_reduced],pma->has_gravitational_terms,bi_reduced_size_counter,pma->radtp_g5)
      }
      class_test(bi_reduced_size_counter!=0,
                 pma->error_message,
                 "Number of radial types for bessel integral type 'reduced' do not match up.");
    }
    if(pma->has_cltp_sh){
      pma->radtp_of_bitp_size[pma->cltp_index_sh*pma->bitp_size+pma->bitp_index_nu_reduced] = 0;
    }
  }
  if(pma->has_bitp_lfactor){
    if(pma->has_cltp_nc){
      class_alloc(pma->radtps_of_bitp[pma->cltp_index_nc*pma->bitp_size+pma->bitp_index_lfactor],
                  bi_lfactor_size_nc_counter*sizeof(double),
                  pma->error_message);
      pma->radtp_of_bitp_size[pma->cltp_index_nc*pma->bitp_size+pma->bitp_index_lfactor] = bi_lfactor_size_nc_counter;
      matter_index_correspondence(pma->radtps_of_bitp[pma->cltp_index_nc*pma->bitp_size+pma->bitp_index_lfactor],pma->has_lensing_terms,bi_lfactor_size_nc_counter,pma->radtp_nclens)
      matter_index_correspondence(pma->radtps_of_bitp[pma->cltp_index_nc*pma->bitp_size+pma->bitp_index_lfactor],(pma->has_stp_delta_m && pma->uses_density_splitting),bi_lfactor_size_nc_counter,pma->radtp_dens2)
      class_test(bi_lfactor_size_nc_counter!=0,
                 pma->error_message,
                 "Number of radial types for bessel integral type 'lfactor' do not match up.");
    }
    if(pma->has_cltp_sh){
      class_alloc(pma->radtps_of_bitp[pma->cltp_index_sh*pma->bitp_size+pma->bitp_index_lfactor],
                  bi_lfactor_size_sh_counter*sizeof(double),
                  pma->error_message);
      pma->radtp_of_bitp_size[pma->cltp_index_sh*pma->bitp_size+pma->bitp_index_lfactor] = bi_lfactor_size_sh_counter;
      matter_index_correspondence(pma->radtps_of_bitp[pma->cltp_index_sh*pma->bitp_size+pma->bitp_index_lfactor],pma->has_cl_shear,bi_lfactor_size_sh_counter,pma->radtp_shlens)
      class_test(bi_lfactor_size_sh_counter!=0,
                 pma->error_message,
                 "Number of radial types for bessel integral type 'lfactor' do not match up.");
    }
  }
  /**
   * Finally we can print what types were found
   *   and what bessel integral types they correspond to
   * */
  if(pma->matter_verbose > MATTER_VERBOSITY_INDICES){
    int i,j,k;
    for(k=0;k<pma->cltp_size;++k){
      if( (!matter_is_index(k,pma->cltp_index_nc,pma->has_cltp_nc))
         &&(!matter_is_index(k,pma->cltp_index_sh,pma->has_cltp_sh))
        ){continue;}
      printf(" -> Searching for correspondences at cltp %i \n",k);
      for(i=0;i<pma->bitp_size;++i){
        for(j=0;j<pma->radtp_of_bitp_size[k*pma->bitp_size+i];++j){
          if(pma->uses_relative_factors){
            if(matter_is_index(pma->radtps_of_bitp[k*pma->bitp_size+i][j],pma->radtp_combined,pma->uses_relative_factors)){
              if(matter_is_index(i,pma->bitp_index_nu_reduced,pma->has_bitp_nu_reduced)){
                printf("    -> Found in bitp 'reduced' index 'combined' \n");
              }
            }
          }
          else{
            if(matter_is_index(pma->radtps_of_bitp[k*pma->bitp_size+i][j],pma->radtp_dens,pma->has_stp_delta_m && (!pma->uses_density_splitting))){
              if(matter_is_index(i,pma->bitp_index_normal,pma->has_bitp_normal)){
                printf("    -> Found in bitp 'normal' index 'density' \n");
              }
              else if(matter_is_index(i,pma->bitp_index_nu_reduced,pma->has_bitp_nu_reduced)){
                printf("    -> Found in bitp 'reduced' index 'density' \n");
              }
            }
            else if(matter_is_index(pma->radtps_of_bitp[k*pma->bitp_size+i][j],pma->radtp_dens1,pma->has_stp_delta_m && pma->uses_density_splitting)){
              if(matter_is_index(i,pma->bitp_index_nu_reduced,pma->has_bitp_nu_reduced)){
                printf("    -> Found in bitp 'reduced' index 'density (part 1)' \n");
              }
            }
            else if(matter_is_index(pma->radtps_of_bitp[k*pma->bitp_size+i][j],pma->radtp_dens2,pma->has_stp_delta_m && pma->uses_density_splitting)){
              if(matter_is_index(i,pma->bitp_index_lfactor,pma->has_bitp_lfactor)){
                printf("    -> Found in bitp 'l(l+1) factor' index 'density (part 2)' \n");
              }
            }
            else if(matter_is_index(pma->radtps_of_bitp[k*pma->bitp_size+i][j],pma->radtp_dop1,pma->has_redshift_space_distortion && (!pma->uses_rsd_combination))){
              if(matter_is_index(i,pma->bitp_index_nu_reduced,pma->has_bitp_nu_reduced)){
                printf("    -> Found in bitp 'reduced' index 'doppler 1' \n");
              }
            }
            else if(matter_is_index(pma->radtps_of_bitp[k*pma->bitp_size+i][j],pma->radtp_dop2,pma->has_redshift_space_distortion && (!pma->uses_rsd_combination))){
              if(matter_is_index(i,pma->bitp_index_nu_reduced,pma->has_bitp_nu_reduced)){
                printf("    -> Found in bitp 'reduced' index 'doppler 2' \n");
              }
            }
            else if(matter_is_index(pma->radtps_of_bitp[k*pma->bitp_size+i][j],pma->radtp_rsd,pma->has_redshift_space_distortion && (!pma->uses_rsd_combination))){
              if(matter_is_index(i,pma->bitp_index_nu_reduced,pma->has_bitp_nu_reduced)){
                printf("    -> Found in bitp 'reduced' index 'RSD (dominant term)' \n");
              }
            }
            else if(matter_is_index(pma->radtps_of_bitp[k*pma->bitp_size+i][j],pma->radtp_rsd_combined,pma->has_redshift_space_distortion && (pma->uses_rsd_combination))){
              if(matter_is_index(i,pma->bitp_index_nu_reduced,pma->has_bitp_nu_reduced)){
                printf("    -> Found in bitp 'reduced' index 'RSD (all combined)' \n");
              }
            }
            else if(matter_is_index(pma->radtps_of_bitp[k*pma->bitp_size+i][j],pma->radtp_g1,pma->has_gravitational_terms)){
              if(matter_is_index(i,pma->bitp_index_nu_reduced,pma->has_bitp_nu_reduced)){
                printf("    -> Found in bitp 'reduced' index 'gr 1' \n");
              }
            }
            else if(matter_is_index(pma->radtps_of_bitp[k*pma->bitp_size+i][j],pma->radtp_g2,pma->has_gravitational_terms)){
              if(matter_is_index(i,pma->bitp_index_nu_reduced,pma->has_bitp_nu_reduced)){
                printf("    -> Found in bitp 'reduced' index 'gr 2' \n");
              }
            }
            else if(matter_is_index(pma->radtps_of_bitp[k*pma->bitp_size+i][j],pma->radtp_g3,pma->has_gravitational_terms)){
              if(matter_is_index(i,pma->bitp_index_nu_reduced,pma->has_bitp_nu_reduced)){
                printf("    -> Found in bitp 'reduced' index 'gr 3' \n");
              }
            }
            else if(matter_is_index(pma->radtps_of_bitp[k*pma->bitp_size+i][j],pma->radtp_g4,pma->has_gravitational_terms)){
              if(matter_is_index(i,pma->bitp_index_nu_reduced,pma->has_bitp_nu_reduced)){
                printf("    -> Found in bitp 'reduced' index 'gr 4' \n");
              }
            }
            else if(matter_is_index(pma->radtps_of_bitp[k*pma->bitp_size+i][j],pma->radtp_g5,pma->has_gravitational_terms)){
              if(matter_is_index(i,pma->bitp_index_nu_reduced,pma->has_bitp_nu_reduced)){
                printf("    -> Found in bitp 'reduced' index 'gr 5' \n");
              }
            }
            else if(matter_is_index(pma->radtps_of_bitp[k*pma->bitp_size+i][j],pma->radtp_nclens,pma->has_cltp_nc && pma->has_lensing_terms)){
              if(matter_is_index(i,pma->bitp_index_lfactor,pma->has_bitp_lfactor)){
                printf("    -> Found in bitp 'l(l+1) factor' index 'number count lensing' \n");
              }
            }
            else if(matter_is_index(pma->radtps_of_bitp[k*pma->bitp_size+i][j],pma->radtp_shlens,pma->has_cl_shear && pma->has_cltp_sh)){
              if(matter_is_index(i,pma->bitp_index_lfactor,pma->has_bitp_lfactor)){
                printf("    -> Found in bitp 'l(l+1) factor' index 'shear lensing' \n");
              }
            }
            //Ifend select radtp
          }
          //Ifend uses relative factors
        }
        //End radtps
      }
      //End bitp
    }
    //End cltp
  }
  return _SUCCESS_;
}
int matter_obtain_bessel_integrals(
                                  struct matters * pma
                                  ){
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS){
    printf("Method :: Obtaining bessel integrals directly (no recursion) \n");
  }
  clock_t start_bessel = clock();
  /**
   * Declare variables to be used later
   * */
  int index_tilt1,index_tilt2,index_tilt1_tilt2;
  int index_coeff;int index_l;
  int index_t;
  int current_bi_size;
  double t_val = 0.0;
  double y_val;
  double y_max;

  double t_min_estimate;
  double y_min_guess;

  double res_real,res_imag;

  /**
   * We sample in between 0.0 and 1.0-BI_SAMPLING_EPSILON
   *
   * The maximum number of samples is how many evaluations
   *  are possible at most.
   * */
  double BI_SAMPLING_EPSILON = 1e-8;
  int MAX_SAMPLES = 1000000+1;
  y_max = -log(BI_SAMPLING_EPSILON);
  /**
   * We allocate the bessel integral arrays,
   *  where we store the final results
   * */
  class_alloc(pma->bi_real,
              pma->tilt_grid_size*sizeof(double**),
              pma->error_message);
  class_alloc(pma->bi_imag,
              pma->tilt_grid_size*sizeof(double**),
              pma->error_message);
  class_alloc(pma->bi_size,
              pma->tilt_grid_size*sizeof(int*),
              pma->error_message);
  class_alloc(pma->bi_sampling,
              pma->tilt_grid_size*sizeof(double**),
              pma->error_message);
  class_alloc(pma->bi_max,
              pma->tilt_grid_size*sizeof(double*),
              pma->error_message);
  for(index_tilt1=0;index_tilt1<pma->tilt_size;++index_tilt1){
    for(index_tilt2=index_tilt1;index_tilt2<pma->tilt_size;++index_tilt2){
      index_tilt1_tilt2 = index_symmetric_matrix(index_tilt1,index_tilt2,pma->tilt_size);
      class_alloc(pma->bi_real[index_tilt1_tilt2],
                  pma->size_fft_result*pma->l_size*sizeof(double*),
                  pma->error_message);
      class_alloc(pma->bi_imag[index_tilt1_tilt2],
                  pma->size_fft_result*pma->l_size*sizeof(double*),
                  pma->error_message);
      class_alloc(pma->bi_size[index_tilt1_tilt2],
                  pma->size_fft_result*pma->l_size*sizeof(int),
                  pma->error_message);
      class_alloc(pma->bi_sampling[index_tilt1_tilt2],
                  pma->size_fft_result*pma->l_size*sizeof(double*),
                  pma->error_message);
      class_alloc(pma->bi_max[index_tilt1_tilt2],
                  pma->size_fft_result*pma->l_size*sizeof(double),
                  pma->error_message);
      if(pma->matter_verbose>MATTER_VERBOSITY_BESSEL){
        printf(" -> Obtaining bessel integrals for tilt %f \n",pma->bias-pma->nu_real[index_tilt1_tilt2]);
      }
      for(index_l=0;index_l<pma->l_size;++index_l){
        if(pma->matter_verbose >MATTER_VERBOSITY_BESSEL){
          printf(" -> Obtaining bessel integrals of l = %4d (index %3d/%3d)\n",
                 (int)pma->l_sampling[index_l],index_l,pma->l_size-1
                );
        }
        for(index_coeff=0;index_coeff<pma->size_fft_cutoff;++index_coeff){
          /**
           * This estimation simply serves to find the minimum t
           *  This we do simply to sample with approximately the same samples,
           *  no matter which t we are using
           * This is actually preferable to the recursion type of evaluation,
           *  but for that we sadly do require the same t point for all l,
           *  and thus this nice feature is sadly not possible there
           * */
          class_call(matter_estimate_t_max_bessel(
                                                  pma,
                                                  pma->l_sampling[index_l],
                                                  pma->nu_imag[index_coeff],
                                                  &t_min_estimate),
                    pma->error_message,
                    pma->error_message);
          /**
           *
           * We calculate the minimum y of our guessed t
           * First we find the sampling step size for a given wanted number of samples.
           *  This also allows us to allocate the maximum size of bessel integrals,
           *  that would be used if there was no t_min
           * Of course we need to cap this maximum number off at some point,
           *  so we do not run into the maximum integer limit
           * However, the real number of samples should still be around the wanted one
           *
           * The real size is found later through the criterion
           *  |I_l(nu,t)| < eps*|I_l(nu,1)|
           * */
          y_min_guess = -log(1.0-t_min_estimate);
          current_bi_size= MIN((int)(pma->bi_wanted_samples*y_max/(y_max-y_min_guess))+1,MAX_SAMPLES);
          pma->bi_size[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff] = current_bi_size;
          class_test(current_bi_size==0,
                     pma->error_message,
                     "the bessel-integral array size was found to be 0. Either this is a bug, or the parameter bi_wanted_samples is too low."
                    );
          /**
           * Allocate the arrays of the maximum possible size
           * */
          class_alloc(pma->bi_real[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                  current_bi_size*sizeof(double),
                  pma->error_message);
          class_alloc(pma->bi_imag[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                  current_bi_size*sizeof(double),
                  pma->error_message);
          class_alloc(pma->bi_sampling[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                  current_bi_size*sizeof(double),
                  pma->error_message);

          /**
           * First we calculate I_l(nu,1) and |I_l(nu,1)|
           * */
          double first_val_real; double first_val_imag;
          bessel_integral_transform(
                          pma->l_sampling[index_l],
                          pma->nu_real[index_tilt1_tilt2],
                          pma->nu_imag[index_coeff]+pma->bessel_imag_offset,
                          1.0,
                          &first_val_real,
                          &first_val_imag
                          );
          double first_val_abs_squared = first_val_real*first_val_real+first_val_imag*first_val_imag;
          /**
           * This point is also our first sampling point
           * */
          pma->bi_sampling[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff][0]= 0.0;
          pma->bi_real[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff][0] = first_val_real;
          pma->bi_imag[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff][0] = first_val_imag;

          /**
           * We iterate through all possible (!) values of t,
           *  and keep those that have an I_l that is large enough
           *  (these t are probably (but not necessarily) larger than t_min_estimate)
           * */
          for(index_t=1;index_t<current_bi_size;++index_t){
            y_val = y_max-y_max*((double)(index_t-1))/((double)(current_bi_size-2));
            t_val = 1.0-exp(-y_val);
            pma->bi_sampling[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff][index_t]= 1.0 - t_val;
            /**
             * The hypergeometric function has problems if nu is close to 0 in absolute value
             * Therefore, we add some infinitessimal 1e-5 to nu_imag to stabilize the function
             *  (The reason is that the arguments degenerate to integer offset (l integer),
             *  and this is where the gamma is near to a pole (divergence)
             * At the exact divergence, we could calculate it with a finite sum with no problem,
             *  but close to it, (say 1e-10 , 1e-20) , this would be a huge problem
             *
             * Of course it would be theoretically possible to do these cases as well
             *  by a suitable transformation, but this would require a lot more checking,
             *  slowing the program further down
             * */
            bessel_integral_transform(
                                      pma->l_sampling[index_l],
                                      pma->nu_real[index_tilt1_tilt2],
                                      pma->nu_imag[index_coeff]+pma->bessel_imag_offset,
                                      t_val,
                                      &res_real,
                                      &res_imag
                                      );
            /**
             * At some point we will hopefully reach |I_l(nu,t)|< eps*|I_l(nu,1)|,
             *  thus finding the ACTUAL size (which should be of order bi_wanted_samples)
             *  If we reach below t_val==0.0, some error has happened
             *  This should never happen.
             * */
            if(res_real*res_real+res_imag*res_imag<BESSEL_EPSILON*first_val_abs_squared || t_val<0.0){
                class_test(t_val<0.0,
                           pma->error_message,
                           "Return due to t<0 in bessel_integral. This is a bug in the code."
                          );
                /**
                 * Now we can store the maximum value of t and array size,
                 *  and finally reallocate the arrays
                 * */
                pma->bi_real[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff][index_t] = res_real;
                pma->bi_imag[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff][index_t] = res_imag;
                pma->bi_size[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff] = index_t+1;
                pma->bi_max[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff] = 1.0 - t_val;
                class_realloc(pma->bi_real[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                              pma->bi_real[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                              pma->bi_size[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff]*sizeof(double),
                              pma->error_message);
                class_realloc(pma->bi_imag[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                              pma->bi_imag[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                              pma->bi_size[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff]*sizeof(double),
                              pma->error_message);
                class_realloc(pma->bi_sampling[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                              pma->bi_sampling[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                              pma->bi_size[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff]*sizeof(double),
                              pma->error_message);
                /**
                 * If we have found |I_l(nu,t)|<eps*|I_l(nu,1)| for this l value,
                 *  we can ignore all t smaller than this one (because there the I_l(nu,t) is smaller in absolute value
                 * (It goes to 0.0 for t->0.0)
                 * */
                break;
            }
            pma->bi_real[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff][index_t] = res_real;
            pma->bi_imag[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff][index_t] = res_imag;

          }
          //End t
          /**
           * If the I_l(nu,t) never statisfied the criterion |I_l(nu,t)|<eps*|I_l(nu,1)|
           *  this is a mistake
           * The last y_val is always 0.0, which corresponds to t = 1.0-exp(-0.0) = 0.0,
           *  so here the function has to be exactly equal to 0.0,
           *  for all l>=2 , which are considered here.
           * Thus finding a mode where the criterion was not fulfilled at any point is
           *  clearly a bug. (Or l<1 was somehow enabled, which should NOT be true)
           * */
          class_test(index_t==current_bi_size,
                     pma->error_message,
                     "The sampling in t-space for the I_l(nu,t) did not stop within the required parameter range ( t>=0.0 ) \n \
                      This is a bug. Either l<1 was activated, or I_l(nu,0) =/= 0, both of which should not happen. \n \
                      Here are a few more parameters that could be useful for debugging : \n \
                      index_t = %i , t =%f , l = %i , nu= %.10e+%.10ej \n \
                      Predicted number of required t steps %i, stepsize %.10e, up to 1.0 \n",
                      index_t,t_val,(int)pma->l_sampling[index_l],pma->nu_real[index_tilt1_tilt2],pma->nu_imag[index_coeff],
                      current_bi_size,1.0/current_bi_size
                    );
          //End if
        }
        //End coeff
      }
      //End l
    }
    //End tilt2
  }
  //End tilt1
  clock_t end_bessel = clock();
  if(pma->matter_verbose > MATTER_VERBOSITY_TIMING ){
    printf(" -> Obtaining Bessel Integrals took %f seconds \n",((double)(end_bessel-start_bessel))/CLOCKS_PER_SEC);
  }
  return _SUCCESS_;
}
//TODO :: add final one with t=0, => Il = 0
int matter_obtain_bessel_recursion_parallel(struct matters* pma){
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS ){
    printf("Method :: Obtain bessel integrals from recursion \n");
  }
  clock_t start_bessel = clock();
  double start_bessel_omp = omp_get_wtime();
  int index_tilt1,index_tilt2,index_tilt1_tilt2;
  int index_coeff;
  int index_l;
  int index_t;
  double y_max;
  double y_min;
  int bessel_recursion_l_size;

  /**
   * The maximum l of recursion that we need to reach
   *  is simply given by the last l value in l sampling
   * */
  double bessel_recursion_l_max = pma->l_sampling[pma->l_size-1];
  /**
   * There is a semi-analytical formula for calculating how
   *  many l values are required to reach a given accuracy
   *  of the bessel integrals.
   * This analytic formula simply explodes for t->1,
   *  which we 'catch' by putting an arbitrary, but very very large
   *  number here
   * It can happen that actually more than this number of l
   *  would be required, giving us higher errors
   * We assume however, that this case is VERY rare
   * */
  int back_complicated_max_size = 25000;
  /**
   * One of the more complicated methods can actually
   * recognize its own failure
   *  (It starts a backward recursion, that can be connected
   *   to the analytic limit for l=0)
   * This is the allowed deviation from the analytic limit
   * */
  double BI_ALLOWED_ERROR = 1e-6;
  /**
   * Some of the simpler recursion techniques become
   *  unsafe for small imaginary parts
   * Here we explicitly switch those off
   * */
  double NU_IMAG_BI_RECURSION_SWITCH = -1.0;//20.0;//20.0;
  /**
   * We want as many samples as the user requested, plus two additional ones:
   *  t = 1.0 and t = 0.0
   * */
  int bi_recursion_t_size = pma->bessel_recursion_t_size+2;
  /**
   * We are going to use a formula of the type
   *  1/(1-t), and rounding that to an integer
   * Of course not all doubles are representable as an integer,
   *  so we are going to do a simple cutoff procedure:
   * We are going to use
   *  1/(1-t+inv_maximum_representable_integer)
   * which can be at most
   *  1/inv_maximum_representable_integer
   *
   * Thus, if our maximum representable integer is
   *  2*10^9, we would set the flag to 0.5e-9
   *
   * However, just to be a bit more careful, we
   *  choose 1e-5
   * */
  double inv_maximum_representable_integer = 1e-5;
  double nu_min_forward_real = 1.5;
  double l_max_forward_real = 1200;
  double l_max_backward_real = 5000;
  /**
   * Now we can set the y_max and y_min variables
   *
   * The cases t=1.0 and t=0.0 are handled completely seperately anyway,
   * so we just need reasonable limits that are not 'too' far off,
   *
   * Here, we choose the immediate values from
   * (BI_MIN_T up to 1.0-BI_SAMPLING_EPSILON)
   * which is inside the interval
   * (0.0,1.0)
   * */
  double BI_SAMPLING_EPSILON = 1e-8;//1e-8
  double BI_MIN_T = 1e-6;//1e-6
  y_max = -log(BI_SAMPLING_EPSILON);
  y_min = -log(1.0-BI_MIN_T);

  /**
   * If our l is sampled logarithmically,
   * the same can not be done for the recursion relations
   * (which require every single l)
   *
   * As such, we can have a different size here
   * (bessel_recursion_l_size instead of pma->l_size)
   * */
  bessel_recursion_l_size = bessel_recursion_l_max+1;
  pma->l_size_recursion = bessel_recursion_l_size;

  /**
   * It is possible to restart the recursion relations
   * for different l, thus reducing the peak error
   *
   * In practice, it slows the program down by more than it
   *  can reduce the error, since the starting points themselves
   *  are error-prone
   *
   * As such, they are currently disabled (very high) with the option
   *  to enable them again if desired
   * */
  int bessel_recursion_backward_min_l_step_low_nu = 500000;
  int bessel_recursion_backward_min_l_step_high_nu;
  int bessel_recursion_backward_max_l_step = 1000000;
  int bessel_recursion_forward_min_l_step_low_nu = 400000;
  int bessel_recursion_forward_min_l_step_high_nu;
  int bessel_recursion_forward_max_l_step = 200000;

  if(pma->has_tilt_normal){
    bessel_recursion_forward_min_l_step_high_nu = 60000;
    bessel_recursion_backward_min_l_step_high_nu = 500000;
  }
  else{
    bessel_recursion_forward_min_l_step_high_nu = 10000;
    bessel_recursion_backward_min_l_step_high_nu = 500000;
  }

  /**
   * Allocate the arrays in which we want to store the final bessel integrals
   *  (Bessel Integrals are shortened to BI)
   * */
  class_alloc(pma->bi_real,
              pma->tilt_grid_size*sizeof(double**),
              pma->error_message);
  class_alloc(pma->bi_imag,
              pma->tilt_grid_size*sizeof(double**),
              pma->error_message);
  class_alloc(pma->bi_size,
              pma->tilt_grid_size*sizeof(int*),
              pma->error_message);
  class_alloc(pma->bi_sampling,
              pma->tilt_grid_size*sizeof(double**),
              pma->error_message);
  class_alloc(pma->bi_max,
              pma->tilt_grid_size*sizeof(double*),
              pma->error_message);
  /**
   * Define and allocate temporary arrays,
   *  which are required for the evaluations
   * */
  for(index_tilt1=0;index_tilt1<pma->tilt_size;++index_tilt1){
    for(index_tilt2=index_tilt1;index_tilt2<pma->tilt_size;++index_tilt2){
      index_tilt1_tilt2 = index_symmetric_matrix(index_tilt1,index_tilt2,pma->tilt_size);

      class_alloc(pma->bi_real[index_tilt1_tilt2],
                  bessel_recursion_l_size*pma->size_fft_result*sizeof(double*),
                  pma->error_message);
      class_alloc(pma->bi_imag[index_tilt1_tilt2],
                  bessel_recursion_l_size*pma->size_fft_result*sizeof(double*),
                  pma->error_message);
      class_alloc(pma->bi_sampling[index_tilt1_tilt2],
                  bessel_recursion_l_size*pma->size_fft_result*sizeof(double*),
                  pma->error_message);
      class_alloc(pma->bi_size[index_tilt1_tilt2],
                  bessel_recursion_l_size*pma->size_fft_result*sizeof(int),
                  pma->error_message);
      class_alloc(pma->bi_max[index_tilt1_tilt2],
                  bessel_recursion_l_size*pma->size_fft_result*sizeof(double),
                  pma->error_message);
      /**
       * Allocate the real and imaginary arrays
       * Also allocate the sampling array
       * */
      for(index_l=0;index_l<=bessel_recursion_l_max;++index_l){
        for(index_coeff=0;index_coeff<pma->size_fft_cutoff;++index_coeff){
          class_alloc(pma->bi_real[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                      bi_recursion_t_size*sizeof(double),
                      pma->error_message);
          class_alloc(pma->bi_imag[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                      bi_recursion_t_size*sizeof(double),
                      pma->error_message);
          class_alloc(pma->bi_sampling[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                      bi_recursion_t_size*sizeof(double),
                      pma->error_message);
        }
      }
      if(pma->matter_verbose > MATTER_VERBOSITY_BESSEL){
        printf(" -> Obtaining recursion starting bessel integrals for tilt %f \n",pma->bias-pma->nu_real[index_tilt1_tilt2]);
      }
      int abort = _FALSE_;
      #pragma omp parallel private(index_l,index_t,index_coeff) firstprivate(y_min,y_max,pma,bi_recursion_t_size,bessel_recursion_l_size)
      {
        double* max_t;
        double* abi_real;
        double* abi_imag;
        double* initial_abs;
        class_alloc_parallel(max_t,
                    bessel_recursion_l_size*sizeof(double),
                    pma->error_message);
        class_alloc_parallel(abi_real,
                    (bessel_recursion_l_size+back_complicated_max_size)*sizeof(double),
                    pma->error_message);
        class_alloc_parallel(abi_imag,
                    (bessel_recursion_l_size+back_complicated_max_size)*sizeof(double),
                    pma->error_message);
        class_alloc_parallel(initial_abs,
                    bessel_recursion_l_size*sizeof(double),
                    pma->error_message);

        double back_simple_time = 0;
        double for_simple_time = 0;
        double complex_time = 0;
        double inverse_time = 0;
        double taylor_time = 0;
        #pragma omp for schedule(dynamic,(pma->size_fft_result>=CHUNK_SIZE*omp_get_max_threads()?CHUNK_SIZE:1))
        for(index_coeff=0;index_coeff<pma->size_fft_cutoff;++index_coeff){

          if(pma->matter_verbose > MATTER_VERBOSITY_BESSEL){
            if(!MATTER_REWRITE_PRINTING){
              printf(" -> Obtaining bessel from nu[%3d/%3d] = %.10e+%.10ej \n",index_coeff,pma->size_fft_result-1,pma->nu_real[index_tilt1_tilt2],pma->nu_imag[index_coeff]);
            }
            else{
              printf("\r -> Obtaining bessel from nu[%3d/%3d] = %.10e+%.10ej ",index_coeff,pma->size_fft_result-1,pma->nu_real[index_tilt1_tilt2],pma->nu_imag[index_coeff]);
              fflush(stdout);
            }
          }

          /**
           * Set the nu_real and nu_imag parameters
           *  (as shorthands for quicker writing)
           * */
          double nu_real = pma->nu_real[index_tilt1_tilt2];
          double nu_imag = pma->nu_imag[index_coeff];
          double nu_fraction = pma->nu_imag[index_coeff]/pma->nu_imag[pma->size_fft_result-1];

          /**
           * Obtain the initial bessel integrals and enter them already into the final array
           * */
          bessel_integral_recursion_initial_abs(bessel_recursion_l_max,nu_real,nu_imag,abi_real,abi_imag,initial_abs);

          for(index_l=0;index_l<=bessel_recursion_l_max;++index_l){
            pma->bi_sampling[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff][0]= 0.0;
            pma->bi_real[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff][0] = abi_real[index_l];
            pma->bi_imag[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff][0] = abi_imag[index_l];
          }
          /**
           * Set the minimum t for which |I_l(nu,t)|<eps |I_l(nu,1)|
           *  conservatively to 0, and correct later
           * */

          memset(max_t,0,bessel_recursion_l_size*sizeof(double));
          /**
           * This flag is going to tell us which is the current
           * maximum l for any given t
           * that does not statisfy the exit condition:
           *  |I_l(nu,t)|<eps |I_l(nu,1)|
           *
           * The second term is given by initial_abs
           * */

          int l_max_cur = bessel_recursion_l_max;
          clock_t TOT = 0;
          clock_t COPY = 0;
          clock_t T = clock();
          /**
           * This flag keeps track of overflows happening during the summation
           * of the hypergeometric functions. If no more overflows occur
           *  ( this flag being set to _FALSE_ )
           * then a simplified version of the summation can be used,
           * which does not check for further overflows.
           *
           * Always initialize as _TRUE_ !
           * */
          short overflow_flag = _TRUE_;
          for(index_t=1;index_t<bi_recursion_t_size;++index_t){
            /**
             * Obtain the t at which we want to sample
             * */
            double y = y_max-(y_max-y_min)*sqrt(((double)(index_t-1))/((double)(bi_recursion_t_size-2)));
            double t = 1.0-exp(-y);

            /**
             * Here are some semi-analytical determinations
             *  of the range to which each method extends
             *
             * Sadly, these only apply for a bias around ~1.9
             *  and for the maximum error of 1e-6
             *
             * Feel free to improve these formuli
             *
             * Forward simple:
             *  Limited to t close to 1,
             *  the deviation of which we fit,
             *  and capped off with a maximum value (that depends on nu_imag)
             *
             * Backward simple:
             *  Limited to t close to 1 (surprisingly),
             *  but with a much broader range
             * Update: Using the overflow-safe version, which requires the
             *  overflow_flag to keep track of overflows during calculation
             *
             * Self inverse taylor:
             *  Limited to very very close to 1,
             *  (otherwise too slow)
             *
             * Complex:
             *  Decides on forward or backward recursion
             *  using matrix invertibility criteria,
             *  but requires sufficient l to shave off initial errors in the
             *  starting conditions
             * */
            double forward_simple_factor_high_real = (nu_imag*nu_imag/10000.0+nu_imag/3.3+4.0);
            double forward_simple_factor_low_real = nu_imag/4.0*0.95;
            double forward_simple_factor = (forward_simple_factor_high_real-forward_simple_factor_low_real)*(nu_real+2.1)/4.0+forward_simple_factor_low_real;

            double forward_l_max_const = 3000.0+nu_imag*40.0;
            int l_max_forward_simple = (int)MIN((forward_simple_factor/(1-t+inv_maximum_representable_integer)),forward_l_max_const);

            double backward_simple_factor = 10.0+nu_imag/5.0;
            int l_max_backward_simple = (int)(backward_simple_factor/(1-t+inv_maximum_representable_integer));

            double self_inverse_taylor_factor = 15.0;
            int l_max_self_inverse_taylor = (int)(self_inverse_taylor_factor/(1-t+inv_maximum_representable_integer));

            int delta_l_required = (int)MIN((12./(1-t)),back_complicated_max_size);

            double backward_simple_lfactor = MAX(1.5-nu_imag/15,0.0);
            if(t<T_MIN_TAYLOR){
              //printf("USING TAYLOR \n\n");
              clock_t func_start = clock();
              double func_start_t = omp_get_wtime();
              bessel_integral_recursion_taylor(l_max_cur,nu_real,nu_imag,t,max_t,initial_abs,abi_real,abi_imag);
              taylor_time += omp_get_wtime()-func_start_t;
              TOT += clock()-func_start;
            }
            /**
             * We want the hypergeometric series to have only very vew terms
             * The arguments are approximately l^2/4 * (1-z)^2/(1+z)^2 << 1
             *  (neglecting any nu dependence)
             * Using z = t^2 = (1-eps)^2 ~ 1-2eps
             *  we find (1-z)^2/(1+z)^2 ~ eps^2
             * Then we find (l*eps/2)^2 << 1
             *  As such we get eps << 2/l => 1-t = alpha*2/l with alpha<<1
             *  We find quick convergence for 1-t = 2*alpha/l < 2*T_MIN_INVERSE_TAYLOR/l_max
             * ( l < l_max , alpha = T_MIN_INVERSE_TAYLOR)
             * */
            else if(l_max_self_inverse_taylor >l_max_cur && t>1.-T_MIN_INVERSE_TAYLOR){
              //printf("USING SELF \n\n");
              clock_t func_start = clock();
              double func_start_t = omp_get_wtime();
              bessel_integral_recursion_inverse_self(l_max_cur,nu_real,nu_imag,t,abi_real,abi_imag,max_t,initial_abs,pma->error_message);
              inverse_time += omp_get_wtime()-func_start_t;
              TOT += clock()-func_start;
            }
            else if(
              ( nu_imag>NU_IMAG_BI_RECURSION_SWITCH &&
              l_max_forward_simple > l_max_cur )
              || (nu_real>nu_min_forward_real && l_max_cur<l_max_forward_real)
            ){
              //printf("USING FORWARD SIMPLE \n\n");
              clock_t func_start = clock();
              double func_start_t = omp_get_wtime();
              class_call_parallel(bessel_integral_recursion_forward_simple(
                                                       l_max_cur,
                                                       nu_real,
                                                       nu_imag,
                                                       t,
                                                       abi_real,
                                                       abi_imag,
                                                       max_t,
                                                       initial_abs,
                                                       pma->error_message),
                         pma->error_message,
                         pma->error_message);
              for_simple_time += omp_get_wtime()-func_start_t;
              TOT += clock()-func_start;
            }
            else if(
              (nu_imag>NU_IMAG_BI_RECURSION_SWITCH &&
              l_max_backward_simple > l_max_cur) ||
              (l_max_cur<l_max_backward_real)
            ){
              //printf("USING BACKWARD SIMPLE \n\n");
              clock_t func_start = clock();
              double func_start_t = omp_get_wtime();
              class_call_parallel(bessel_integral_recursion_backward_simple_safe(
                                                       l_max_cur,
                                                       (1.1+backward_simple_lfactor)*l_max_cur,
                                                       nu_real,
                                                       nu_imag,
                                                       t,
                                                       abi_real,
                                                       abi_imag,
                                                       max_t,
                                                       initial_abs,
                                                       &overflow_flag,
                                                       pma->error_message),
                         pma->error_message,
                         pma->error_message);
              TOT += clock()-func_start;
              back_simple_time += omp_get_wtime()-func_start_t;
            }
            else{
                //printf("USING COMPLICATED \n\n");
                clock_t func_start = clock();
                double func_start_t = omp_get_wtime();
                class_call_parallel(bessel_integral_recursion_complicated(l_max_cur,
                                                                 l_max_cur+delta_l_required-1,
                                                                 nu_real,
                                                                 nu_imag,
                                                                 t,
                                                                 BI_ALLOWED_ERROR,
                                                                 abi_real,
                                                                 abi_imag,
                                                                 max_t,
                                                                 initial_abs,
                                                                 pma->error_message),
                           pma->error_message,
                           pma->error_message);
              TOT += clock()-func_start;
              complex_time += omp_get_wtime()-func_start_t;
            }
            /**
             * After obtaining the corresponding bessel integrals through recursion
             *  for this particular value for t,
             *  we need to store them within the storage arrays
             * */
            clock_t copy_start = clock();
            for(index_l=0;index_l<=l_max_cur;++index_l){
              pma->bi_sampling[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff][index_t]= 1.0-t;
              pma->bi_real[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff][index_t] = abi_real[index_l];
              pma->bi_imag[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff][index_t] = abi_imag[index_l];
              if(max_t[index_l]>=t){
                /**
                 * If the condition |I_l(nu,t)|<eps*|I_l(nu,1)|
                 * is fulfilled, we do not need to evaluate this mode
                 * for any smaller values of t
                 *
                 * Thus this mode is 'exiting' our evaluation range
                 *
                 * We can also then reallocate the arrays reaching to exactly this point
                 * And also set the size and maximum evaluable point before 0
                 * */
                pma->bi_size[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff] = index_t+1;
                pma->bi_max[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff] = 1.0 - t;

                //TODO :: figure out why parallel reallocation leads to segmentation faults
                /*class_realloc_parallel(pma->bi_real[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                              pma->bi_real[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                              pma->bi_size[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff]*sizeof(double),
                              pma->error_message);
                class_realloc_parallel(pma->bi_imag[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                              pma->bi_imag[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                              pma->bi_size[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff]*sizeof(double),
                              pma->error_message);
                class_realloc_parallel(pma->bi_sampling[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                              pma->bi_sampling[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                              pma->bi_size[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff]*sizeof(double),
                              pma->error_message);*/
              }
              //Ifend exit
            }
            COPY += clock()-copy_start;
            //End l
            /**
             * If a mode has 'exited',
             *  we need to check what the current
             *  maximum mode to still be evolved is
             * */
            index_l=l_max_cur;
            while(index_l>0 && max_t[index_l]>=t){
              index_l--;
              l_max_cur--;
            }
          }
          //End t
          T = clock()-T;

          /**
           * The analytic solutions at l=0 and l=1
           *  continue to t=1 without reaching their criterion
           * (The l=0 goes towards a constant, the l=1 decreases too slowly to exit in most cases)
           * */
          pma->bi_size[index_tilt1_tilt2][0*pma->size_fft_result+index_coeff] = bi_recursion_t_size;
          pma->bi_max[index_tilt1_tilt2][0*pma->size_fft_result+index_coeff] = 0.0;
          pma->bi_size[index_tilt1_tilt2][1*pma->size_fft_result+index_coeff] = bi_recursion_t_size;
          pma->bi_max[index_tilt1_tilt2][1*pma->size_fft_result+index_coeff] = 0.0;

        }
        //End coeff
        if(MATTER_REWRITE_PRINTING){
          printf("\r                                                                             \n");
        }
        /*
        printf("FORWARD = %.10e, (per coeff = %.10e) \n",for_simple_time,for_simple_time/(pma->size_fft_cutoff));
        printf("BACKWARD = %.10e, (per coeff = %.10e) \n",back_simple_time,back_simple_time/(pma->size_fft_cutoff));
        printf("COMPLEX = %.10e, (per coeff = %.10e) \n",complex_time,complex_time/(pma->size_fft_cutoff));
        printf("INVERSE = %.10e, (per coeff = %.10e) \n",inverse_time,inverse_time/(pma->size_fft_cutoff));
        printf("TAYLOR = %.10e, (per coeff = %.10e) \n",taylor_time,taylor_time/(pma->size_fft_cutoff));
        */

        //printf("TIME AVG = %.10e secs \n",(double)(timecount)/CLOCKS_PER_SEC/pma->size_fft_cutoff);
        free(initial_abs);
        free(max_t);
        free(abi_real);
        free(abi_imag);
      }
      //End parallel
      if (abort == _TRUE_) return _FAILURE_;
    }
    //End tilt2
  }
  //End tilt1
  /**
   * Delete temporary arrays
   * */
  clock_t end_bessel = clock();
  double end_bessel_omp = omp_get_wtime();
  if(pma->matter_verbose > MATTER_VERBOSITY_TIMING) {
    printf(" -> Obtaining (recursion) Bessel Integrals took %f seconds \n",((double)(end_bessel-start_bessel))/CLOCKS_PER_SEC);
    printf(" -> Obtaining (recursion) Bessel Integrals took %f REAL seconds \n",end_bessel_omp-start_bessel_omp);
  }
  return _SUCCESS_;
}
int matter_obtain_growth_factor_k_weights(
                                struct matters* pma,
                                struct perturbs* ppt,
                                double* k_weights
                                ){
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS) {
    printf("Method :: Obtain k weights for growth factors \n");
  }
  /**
   * Define variables to be used later
   * */
  int index_k;
  double log_k_max;
  double log_k_min;
  double k_weight_sum=0.0;
  /**
   * Now we want to make sure that the k weight min and max are properly defined
   * */
  class_test(pma->k_weight_k_min<=0,
             pma->error_message,
             "k_weight_k_min has to be positive");
  class_test(pma->k_weight_k_max<=pma->k_weight_k_min,
             pma->error_message,
             "k_weight_k_max has to be bigger than k_weight_k_min"
            );
  log_k_max = log(pma->k_weight_k_max);
  log_k_min = log(pma->k_weight_k_min);
  /**
   * For either mode, we obtain the weights, and sum them up to normalize them
   *
   * matter_k_weights_gaussian: Gaussian window with 3 sigma intervals down to k_min and up to k_max
   *
   * matter_k_weights_step: Step-function, equal weights between k_min and k_max
   *
   * */
  if(pma->k_weight_mode == matter_k_weights_gaussian){
    double log_k_sigma = (log_k_max-log_k_min)/6.0;
    double log_k_mid = 0.5*(log_k_max+log_k_min);
    for(index_k = 0; index_k < pma->k_size; ++index_k){
      double log_k = log(ppt->k[ppt->index_md_scalars][index_k]);
      k_weights[index_k] = exp(-0.5*(log_k-log_k_mid)*(log_k-log_k_mid)/(log_k_sigma*log_k_sigma));
      k_weight_sum+=k_weights[index_k];
    }
    //End k
  }
  else if(pma->k_weight_mode == matter_k_weights_step){
    for(index_k = 0; index_k < pma->k_size; ++index_k){
      double k=ppt->k[ppt->index_md_scalars][index_k];
      if(k<pma->k_weight_k_max && k>pma->k_weight_k_min){
        k_weights[index_k] = 1.0;
      }
      else{
        k_weights[index_k] = 0.0;
      }
      k_weight_sum+=k_weights[index_k];
    }
    //End k
  }
  else{
    class_stop(pma->error_message,
               "k_weight_mode parameter is none of the selectable options"
               );
  }
  //Ifend weight mode selection
  /**
   * Now we normalize the weights
   * */
  for(index_k = 0; index_k < pma->k_size; ++index_k){
    k_weights[index_k]/=k_weight_sum;
  }
  return _SUCCESS_;
}
int matter_obtain_growth_factor(
                                struct matters* pma,
                                double ** sources,
                                double * k_weights
                                ){
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS){
    printf("Method :: Obtaining growth factor \n");
  }
  /**
   * Loop over all sources and obtain growth factor as
   *  D(tau) = Sum ( S_X(tau,k)/S_X(tau0,k) w_k )
   *
   * */
  int index_ic,index_stp,index_k,index_tau;
  for (index_ic = 0; index_ic < pma->ic_size; index_ic++) {
    for (index_stp = 0; index_stp < pma->stp_size; index_stp++) {
      class_alloc(pma->growth_factor_tau[index_ic * pma->stp_size + index_stp],
                  pma->tau_size*sizeof(double),
                  pma->error_message);
      for(index_tau = 0; index_tau < pma->tau_size;++index_tau){
        pma->growth_factor_tau[index_ic * pma->stp_size + index_stp][index_tau] = 0.0;
        for(index_k = 0; index_k < pma->k_size; ++index_k){
          pma->growth_factor_tau[index_ic * pma->stp_size + index_stp]
                  [index_tau] +=
          k_weights[index_k] *
          sources[index_ic * pma->stp_size + index_stp]
                  [index_tau*pma->k_size+index_k]
          /sources[index_ic * pma->stp_size + index_stp]
                  [(pma->tau_size-1)*pma->k_size+index_k];
        }
        //End coeff
      }
      //End tau
    }
    //End stp
  }
  //End ic
  return _SUCCESS_;
}
int matter_obtain_relative_factor(
                                  struct matters* pma,
                                  double** sources,
                                  double* k_weights
                                 ){
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS){
    printf("Method :: Obtaining relative factor\n");
  }
  int index_ic,index_stp,index_k;
  int index_stp_relative = pma->index_relative_stp;
  class_alloc(pma->relative_factors,
              pma->stp_size*sizeof(double),
              pma->error_message);
  /**
   * Here we obtain the relative factors between the different
   *  sources using lambda = Sum(S_X(tau0,k)/S_Y(tau0,k)),
   *  Here Y is the type we want to compute the relative factor to
   * */
  for (index_ic = 0; index_ic < pma->ic_size; index_ic++) {
    for(index_stp = 0; index_stp<pma->stp_size;++index_stp){
      if(index_stp==index_stp_relative){
        pma->relative_factors[index_ic*pma->stp_size+index_stp]=1.0;
      }
      else{
        pma->relative_factors[index_ic*pma->stp_size+index_stp]=0.0;
        for(index_k = 0; index_k < pma->k_size; ++index_k){
          pma->relative_factors[index_ic * pma->stp_size + index_stp]+=
          k_weights[index_k] *
          sources[index_ic * pma->stp_size + index_stp]
                  [(pma->tau_size-1)*pma->k_size+index_k]
          /sources[index_ic * pma->stp_size + index_stp_relative]
                  [(pma->tau_size-1)*pma->k_size+index_k];
        }
        //End k
      }
      //End iff relative
    }
    //End stp
  }
  //End ic
  return _SUCCESS_;
}
int matter_obtain_nonseparability(
                               struct matters* pma,
                               double ** fft_coeff_real,
                               double ** fft_coeff_imag
                              ){
  if(pma->matter_verbose>MATTER_VERBOSITY_FUNCTIONS){
    printf("Method :: Obtain non-scale-invariance factor\n");
  }
  double** fft_coeff_factor_real;
  double** fft_coeff_factor_imag;
  class_alloc(fft_coeff_factor_real,
              pma->ic_ic_size*pma->stp_grid_size*sizeof(double*),
              pma->error_message);
  class_alloc(fft_coeff_factor_imag,
              pma->ic_ic_size*pma->stp_grid_size*sizeof(double*),
              pma->error_message);
  int index_ic1,index_ic2,index_ic1_ic2;
  int index_tau1,index_tau2,index_tau1_tau2,index_coeff;
  int index_stp1,index_stp2,index_stp1_stp2;
  for (index_ic1 = 0; index_ic1 < pma->ic_size; index_ic1++) {
    for (index_ic2 = index_ic1; index_ic2 < pma->ic_size; index_ic2++) {
      index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,pma->ic_size);
      for (index_stp1 = 0; index_stp1 < pma->stp_size; index_stp1++) {
        for (index_stp2 = 0; index_stp2 < pma->stp_size; index_stp2++) {
          index_stp1_stp2 = index_stp2*pma->stp_size+index_stp1;
          class_alloc(fft_coeff_factor_real[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2],
                pma->size_fft_input*(pma->tau_grid_size+1)*sizeof(double),
                pma->error_message);
          class_alloc(fft_coeff_factor_imag[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2],
                pma->size_fft_input*(pma->tau_grid_size+1)*sizeof(double),
                pma->error_message);
          for(index_tau1=0;index_tau1<pma->tau_size;++index_tau1){
            for(index_tau2=0;index_tau2<pma->tau_size;++index_tau2){
              index_tau1_tau2 = index_tau2*pma->tau_size+index_tau1;
              for(index_coeff = 0; index_coeff < pma->size_fft_result;++index_coeff){
                fft_coeff_factor_real[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2][index_tau1_tau2*pma->size_fft_input+index_coeff] =
                  fft_coeff_real[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2][index_tau1_tau2*pma->size_fft_input+index_coeff]/
                  (pma->growth_factor_tau[index_ic1*pma->stp_size+index_stp1][index_tau1]*pma->growth_factor_tau[index_ic2*pma->stp_size+index_stp2][index_tau2]);
                fft_coeff_factor_imag[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2][index_tau1_tau2*pma->size_fft_input+index_coeff] =
                  fft_coeff_imag[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2][index_tau1_tau2*pma->size_fft_input+index_coeff]/
                  (pma->growth_factor_tau[index_ic1*pma->stp_size+index_stp1][index_tau1]*pma->growth_factor_tau[index_ic2*pma->stp_size+index_stp2][index_tau2]);
              }
              //End index_coeff
            }
            //End tau 1
          }
          //End tau 2
          free(fft_coeff_real[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2]);
          free(fft_coeff_imag[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2]);
          fft_coeff_real[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2] = fft_coeff_factor_real[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2];
          fft_coeff_imag[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2] = fft_coeff_factor_imag[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2];
        }
        //End stp2
      }
      //End stp1
    }
    //End ic2
  }
  //End ic1
  free(fft_coeff_factor_real);
  free(fft_coeff_factor_imag);
  return _SUCCESS_;
}
int matter_obtain_perturbation_sources(
                                      struct background* pba,
                                      struct perturbs * ppt,
                                      struct nonlinear * pnl,
                                      struct matters * pma,
                                      double ** sources
                                      ) {
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS){
    printf("Method :: Obtain perturbation sources\n");
  }
  /**
   * Define indices to use
   * */
  int index_ic;
  int index_tp;
  int index_stp;
  int index_k;
  int index_tau_matter;
  int index_tau_perturbs;
  int last_index_tau;
  double tau_beginning;
  double tau_fraction;
  int index_md = ppt->index_md_scalars;

  /**
   * Allocate temporary source arrays used for splining
   * */
  double * source_temp;
  double * ddsource_temp;
  class_alloc(source_temp,
              ppt->tau_size*pma->k_size*sizeof(double),
              pma->error_message);
  class_alloc(ddsource_temp,
              ppt->tau_size*pma->k_size*sizeof(double),
              pma->error_message);
  if(pma->matter_verbose > MATTER_VERBOSITY_RANGES){
    printf(" -> Allocated perturbation sources with size %i*%i = %i \n",pma->ic_ic_size,pma->stp_size,pma->ic_size*pma->stp_size);
  }
  for (index_ic = 0; index_ic < pma->ic_size; index_ic++) {
    for (index_stp = 0; index_stp < pma->stp_size; index_stp++) {
      index_tp = pma->index_perturb_tp_of_stp[index_stp];
      class_alloc(sources[index_ic * pma->stp_size + index_stp],
                  pma->k_size*pma->tau_size*sizeof(double),
                  pma->error_message);
      /**
       * If we want to use nonlinear spectra,
       *  we need to respect those when copying the sources
       * */
      if (pnl->method != nl_none) {
        for(index_tau_perturbs=0;index_tau_perturbs<ppt->tau_size;++index_tau_perturbs){
          for (index_k=0; index_k<ppt->k_size[index_md]; index_k++) {
            if(pba->has_ncdm){
              //Here we trust nonlinear to set index_pk_cb correctly
              source_temp[index_k*ppt->tau_size+index_tau_perturbs] =
              ppt->sources[index_md][index_ic * ppt->tp_size[index_md] + index_tp]
                [index_tau_perturbs * ppt->k_size[index_md] + index_k]
              * pnl->nl_corr_density[pnl->index_pk_cb][index_tau_perturbs * ppt->k_size[index_md] + index_k];
            }
            else{
              source_temp[index_k*ppt->tau_size+index_tau_perturbs] =
              ppt->sources[index_md][index_ic * ppt->tp_size[index_md] + index_tp]
                [index_tau_perturbs * ppt->k_size[index_md] + index_k]
              * pnl->nl_corr_density[pnl->index_pk_m][index_tau_perturbs * ppt->k_size[index_md] + index_k];
            }
          }
          //End k
        }
        //End tau perturbs
      }
      /**
       * Otherwise, we simply copy the perturbation sources
       * */
      else {
        for(index_tau_perturbs=0;index_tau_perturbs<ppt->tau_size;++index_tau_perturbs){
          for (index_k=0; index_k<ppt->k_size[index_md]; index_k++) {
            source_temp[index_k*ppt->tau_size+index_tau_perturbs] =
            ppt->sources[index_md][index_ic * ppt->tp_size[index_md] + index_tp]
              [index_tau_perturbs * ppt->k_size[index_md] + index_k];
          }
          //End k
        }
        //End tau perturbs
      }
      //Ifend nonlinear check

      /**
       * We want to work with the sources
       *  phi k^2 and psi k^2, not with the original sources
       * This allows us to use the same bias for all sources,
       *  which simplifies a lot of things
       * */
      if(
          matter_is_index(index_stp, pma->stp_index_phi_plus_psi, pma->has_stp_phi_plus_psi)
          ||
          matter_is_index(index_stp,pma->stp_index_phi,pma->has_gravitational_terms)
          ||
          matter_is_index(index_stp,pma->stp_index_phi_plus_psi,pma->has_gravitational_terms)
          ||
          matter_is_index(index_stp,pma->stp_index_phi_prime,pma->has_gravitational_terms)
          ||
          matter_is_index(index_stp,pma->stp_index_psi,pma->has_gravitational_terms)
        ){
        for(index_k=0;index_k<ppt->k_size[index_md];++index_k){
          for(index_tau_perturbs=0;index_tau_perturbs<ppt->tau_size;++index_tau_perturbs){
            source_temp[index_k*ppt->tau_size+index_tau_perturbs] *= ppt->k[index_md][index_k]*ppt->k[index_md][index_k];
          }
        }
      }
      /**
       * Now we spline all sources so we can interpolate them at any k we want
       * */
      class_call(array_spline_table_columns(
                                 ppt->tau_sampling,
                                 ppt->tau_size,
                                 source_temp,
                                 ppt->k_size[index_md],
                                 ddsource_temp,
                                 _SPLINE_EST_DERIV_,
                                 pma->error_message
                                ),
                 pma->error_message,
                 pma->error_message);
      /**
       * Now we spline the sources for every tau we need to know them on
       * */
      for(index_tau_matter=0;index_tau_matter<pma->tau_size;++index_tau_matter){
        /**
         * This is a rough estimation of the index in the tau array of the perturbation structure
         * */
        tau_beginning = ppt->tau_sampling[pma->index_tau_perturbs_beginning];
        tau_fraction = (log(pma->tau_sampling[index_tau_matter])-log(tau_beginning))/(log(pma->tau0)-log(tau_beginning));
        last_index_tau = pma->index_tau_perturbs_beginning+(ppt->tau_size-1-pma->index_tau_perturbs_beginning)*tau_fraction;

        class_call(matter_interpolate_spline_growing_hunt_transposed(
                                               ppt->tau_sampling,
                                               ppt->tau_size,
                                               source_temp,
                                               ddsource_temp,
                                               pma->k_size,
                                               pma->tau_sampling[index_tau_matter],
                                               &last_index_tau,
                                               sources[index_ic*pma->stp_size+index_stp]+index_tau_matter*pma->k_size,
                                               _FALSE_,
                                               pma->error_message
                                               ),
                   pma->error_message,
                   pma->error_message);
      }
      //End tau matter
    }
    //End stp
  }
  //End ic
  free(source_temp);
  free(ddsource_temp);
  return _SUCCESS_;
}

int matter_extrapolate_sources(
                               struct background* pba,
                               struct precision* ppr,
                               struct perturbs* ppt,
                               struct matters * pma,
                               double** k_extrapolated,
                               double** sources,
                               double** extrapolated_sources,
                               short extrapolation_type
                               ){
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS){
    printf("Method :: Source extrapolation \n");
  }
  int index_ic;
  int index_stp;
  int index_tau;
  int k_size_extr;
  *k_extrapolated=NULL;

  if(pma->matter_verbose > MATTER_VERBOSITY_RANGES){
    printf(" -> Getting extrapolated source size : ");
  }
  class_call(get_extrapolated_source_size(
                                          ppr->k_per_decade_for_pk,
                                          pma->k_max,
                                          pma->k_max_extr,
                                          ppt->k_size[ppt->index_md_scalars],
                                          &k_size_extr,
                                          pma->error_message),
             pma->error_message,
             pma->error_message);
  if(pma->matter_verbose > MATTER_VERBOSITY_RANGES){
    printf(" %i \n",k_size_extr);
  }
  class_alloc(*k_extrapolated,
              k_size_extr*sizeof(double),
              pma->error_message);
  class_call(extrapolate_k(ppt->k[ppt->index_md_scalars],
                           pma->k_size,
                           *k_extrapolated,
                           ppr->k_per_decade_for_pk,
                           pma->k_max_extr,
                           pma->error_message),
             pma->error_message,
             pma->error_message);
  for(index_ic=0;index_ic<pma->ic_size;++index_ic){
    for(index_stp=0;index_stp<pma->stp_size;++index_stp){
      class_alloc(extrapolated_sources[index_ic*pma->stp_size+index_stp],
                  pma->tau_size*k_size_extr*sizeof(double),
                  pma->error_message);
      for(index_tau=0;index_tau<pma->tau_size;++index_tau){
        class_call(extrapolate_source(*k_extrapolated,
                                      pma->k_size,
                                      k_size_extr,
                                      sources[index_ic*pma->stp_size+index_stp]+index_tau*ppt->k_size[ppt->index_md_scalars],
                                      extrapolation_type,
                                      extrapolated_sources[index_ic*pma->stp_size+index_stp]+index_tau*k_size_extr,
                                      pba->a_eq*pba->H_eq,
                                      pma->h,
                                      pma->error_message),
                   pma->error_message,
                   pma->error_message);
      }
      //end tau
    }
    //end stp
  }
  //end ic
  pma->k_size = k_size_extr;
  return _SUCCESS_;
}
int matter_spline_bessel_integrals(
                                   struct matters * pma
                                  ) {
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS){
    printf("Method :: Spline bessel integrals\n");
  }
  clock_t spline_start = clock();
  double spline_start_omp = omp_get_wtime();
  /**
   * Initialize and allocate initial arrays
   * */
  int index_coeff;
  int index_l;
  int index_tilt1,index_tilt2,index_tilt1_tilt2;
  for(index_tilt1=0;index_tilt1<pma->tilt_size;++index_tilt1){
    for(index_tilt2=index_tilt1;index_tilt2<pma->tilt_size;++index_tilt2){
      index_tilt1_tilt2 = index_symmetric_matrix(index_tilt1,index_tilt2,pma->tilt_size);
      class_alloc(pma->ddbi_real[index_tilt1_tilt2],
                  pma->l_size*pma->size_fft_result*sizeof(double*),
                  pma->error_message);
      class_alloc(pma->ddbi_imag[index_tilt1_tilt2],
                  pma->l_size*pma->size_fft_result*sizeof(double*),
                  pma->error_message);
      for(index_l=0;index_l<pma->l_size;++index_l){
        for(index_coeff=0;index_coeff<pma->size_fft_cutoff;++index_coeff){
          class_alloc(pma->ddbi_real[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                      pma->bi_size[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff]*sizeof(double),
                      pma->error_message);
          class_alloc(pma->ddbi_imag[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                      pma->bi_size[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff]*sizeof(double),
                      pma->error_message);
          class_call(array_spline_table_columns(pma->bi_sampling[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                                                pma->bi_size[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                                                pma->bi_real[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                                                1,
                                                pma->ddbi_real[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                                                _SPLINE_EST_DERIV_,
                                                pma->error_message),
                    pma->error_message,
                    pma->error_message);
          class_call(array_spline_table_columns(pma->bi_sampling[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                                                pma->bi_size[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                                                pma->bi_imag[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                                                1,
                                                pma->ddbi_imag[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                                                _SPLINE_EST_DERIV_,
                                                pma->error_message),
                    pma->error_message,
                    pma->error_message);
        }
        //End coeff
      }
      //End l
    }
    //End tilt2
  }
  //End tilt1
  clock_t spline_end = clock();
  double spline_end_omp = omp_get_wtime();
  if(pma->matter_verbose > MATTER_VERBOSITY_TIMING){
    printf(" -> Splining bessel integrals took %f CPU  seconds \n",((double)(spline_end-spline_start))/CLOCKS_PER_SEC);
    printf(" -> Splining bessel integrals took %f REAL seconds \n",spline_end_omp-spline_start_omp);
  }
  return _SUCCESS_;
}
int matter_spline_bessel_integrals_recursion(
                                  struct matters * pma
                                  ) {
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS){
    printf("Method :: Spline (recursion) bessel integrals\n");
  }
  clock_t spline_start = clock();
  double spline_start_omp = omp_get_wtime();
  /**
   * Define indices and allocate arrays
   * */
  int index_coeff;
  int index_l;
  int index_tilt1,index_tilt2,index_tilt1_tilt2;
  for(index_tilt1=0;index_tilt1<pma->tilt_size;++index_tilt1){
    for(index_tilt2=index_tilt1;index_tilt2<pma->tilt_size;++index_tilt2){
      index_tilt1_tilt2 = index_symmetric_matrix(index_tilt1,index_tilt2,pma->tilt_size);
      class_alloc(pma->ddbi_real[index_tilt1_tilt2],
                  pma->l_size_recursion*pma->size_fft_result*sizeof(double*),
                  pma->error_message);
      class_alloc(pma->ddbi_imag[index_tilt1_tilt2],
                  pma->l_size_recursion*pma->size_fft_result*sizeof(double*),
                  pma->error_message);
      for(index_l=0;index_l<pma->l_size_recursion;++index_l){
        for(index_coeff=0;index_coeff<pma->size_fft_cutoff;++index_coeff){
          class_alloc(pma->ddbi_real[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                      pma->bi_size[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff]*sizeof(double),
                      pma->error_message);
          class_alloc(pma->ddbi_imag[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                      pma->bi_size[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff]*sizeof(double),
                      pma->error_message);
        }
      }
      int abort = _FALSE_;
      #pragma omp parallel private(index_l,index_coeff) firstprivate(pma)
      {
      #pragma omp for
      for(index_l=0;index_l<pma->l_size_recursion;++index_l){
        for(index_coeff=0;index_coeff<pma->size_fft_cutoff;++index_coeff){
          class_call_parallel(array_spline_table_columns(pma->bi_sampling[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                                                pma->bi_size[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                                                pma->bi_real[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                                                1,
                                                pma->ddbi_real[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                                                _SPLINE_EST_DERIV_,
                                                pma->error_message),
                    pma->error_message,
                    pma->error_message);
          class_call_parallel(array_spline_table_columns(pma->bi_sampling[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                                                pma->bi_size[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                                                pma->bi_imag[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                                                1,
                                                pma->ddbi_imag[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                                                _SPLINE_EST_DERIV_,
                                                pma->error_message),
                    pma->error_message,
                    pma->error_message);
        }
        //End coeff
      }
      //End l
      }
      if(abort == _TRUE_) {return _FAILURE_;}
    }
    //End tilt2
  }
  //End tilt1
  clock_t spline_end = clock();
  double spline_end_omp = omp_get_wtime();
  if(pma->matter_verbose > MATTER_VERBOSITY_TIMING ){
    printf(" -> Splining bessel integrals (recursion) took %f CPU  seconds \n",((double)(spline_end-spline_start))/CLOCKS_PER_SEC);
    printf(" -> Splining bessel integrals (recursion) took %f REAL seconds \n",spline_end_omp-spline_start_omp);
  }
  return _SUCCESS_;
}
 /**
  * interpolate to get y_i(x), x is arranged in growing order, and the point x is presumably close (but maybe not so close) to the previous point x from the last call of this function.
  *
  */
int matter_interpolate_spline_growing_hunt(
					     double * x_array,
					     int x_size,
					     double * array,
					     double * array_splined,
					     int y_size,
					     double x,
					     int * last_index,
					     double * result,
               int IS_PRINT,
					     ErrorMsg errmsg
               ) {

  int inf,sup,mid,index_y,inc;
  double h,a,b;
  //Deal with arrays of different size
  if(*last_index>=x_size) *last_index=x_size-1;
  if(*last_index<0) *last_index=0;

  inc=1;
  if (x >= x_array[*last_index]) {
    if (x > x_array[x_size-1]) {
      sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,
	      x,x_array[x_size-1]);
      return _FAILURE_;
    }
    /* try closest neighboor upward */
    inf = *last_index;
    sup = inf + inc;
    if (x > x_array[sup]) {
      /* hunt upward */
      while (x > x_array[sup]) {
        inf = sup;
        inc += 1;
        sup += inc;
        if (sup > x_size-1) {
          sup = x_size-1;
        }
      }
      /* bisect */
      while (sup-inf > 1) {
        mid=(int)(0.5*(inf+sup));
        if (x < x_array[mid]) {sup=mid;}
        else {inf=mid;}
      }
    }
  }
  else {
    if (x < x_array[0]) {
      sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,
	      x,x_array[0]);
      return _FAILURE_;
    }
    /* try closest neighboor downward */
    sup = *last_index;
    inf = sup - inc;
    if (x < x_array[inf]) {
      /* hunt downward */
      while (x < x_array[inf]) {
        sup = inf;
        inc += 1;
        inf -= inc;
        if (inf < 0) {
          inf = 0;
        }
      }
      /* bisect */
      while (sup-inf > 1) {
        mid=(int)(0.5*(inf+sup));
        if (x < x_array[mid]) {sup=mid;}
        else {inf=mid;}
      }
    }
  }

  *last_index = inf;

  h = x_array[sup] - x_array[inf];
  b = (x-x_array[inf])/h;
  a = 1.0-b;
  for (index_y=0; index_y<y_size; index_y++){
    *(result+index_y) =
      a * *(array+inf+x_size*index_y) +
      b * *(array+sup+x_size*index_y) +
      ((a*a*a-a)* *(array_splined+inf+x_size*index_y) +
       (b*b*b-b)* *(array_splined+sup+x_size*index_y))*h*h/6.;
  }
  return _SUCCESS_;
}
int matter_interpolate_spline_growing_hunt_transposed(
					     double * x_array,
					     int x_size,
					     double * array,
					     double * array_splined,
					     int y_size,
					     double x,
					     int * last_index,
					     double * result,
               int IS_PRINT,
					     ErrorMsg errmsg
               ) {
  double h,a,b;
  int index_y;
  int inf,sup;
  matter_spline_hunt(x_array,x_size,x,last_index,&h,&a,&b,errmsg);
  inf = *last_index;
  sup = inf+1;
  for (index_y=0; index_y<y_size; index_y++){
    *(result+index_y) =
      a * *(array+inf+x_size*index_y) +
      b * *(array+sup+x_size*index_y) +
      ((a*a*a-a)* *(array_splined+inf+x_size*index_y) +
       (b*b*b-b)* *(array_splined+sup+x_size*index_y))*h*h/6.;
  }
  return _SUCCESS_;
}
int matter_sample_sources(
                        struct perturbs * ppt,
                        struct matters* pma,
                        double** source,
                        double** sampled_source,
                        double* perturbed_k_sampling
                        ) {
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS ){
    printf("Method :: Sample sources\n");
  }
  int index_ic;
  int index_stp;
  int index_coeff;
  int last_index=0;
  double* ddsource;
  class_alloc(ddsource,
              pma->k_size*pma->tau_size*sizeof(double),
              pma->error_message);
  for(index_ic=0;index_ic<pma->ic_size;++index_ic){
    for(index_stp=0;index_stp<pma->stp_size;++index_stp){
      class_alloc(sampled_source[index_ic*pma->stp_size+index_stp],
            pma->size_fft_input*pma->tau_size*sizeof(double),
            pma->error_message);
      /**
       * Spline the sources for all values of tau at once
       * */
      class_call(array_spline_table_columns(perturbed_k_sampling,
                                            pma->k_size,
                                            source[index_ic * pma->stp_size + index_stp],
                                            pma->tau_size,
                                            ddsource,
                                            _SPLINE_EST_DERIV_,
                                            pma->error_message),
                 pma->error_message,
                 pma->error_message);
      for(index_coeff=0;index_coeff<pma->size_fft_input;++index_coeff){
        /**
         * At every value of k interpolate for all values of tau at once
         * */
        class_call(matter_interpolate_spline_growing_hunt(
                                      perturbed_k_sampling,
                                      pma->k_size,
                                      source[index_ic * pma->stp_size + index_stp],
                                      ddsource,
                                      pma->tau_size,
                                      pma->k_sampling[index_coeff],
                                      &last_index,
                                      sampled_source[index_ic * pma->stp_size + index_stp]+index_coeff*pma->tau_size,
                                      0,
                                      pma->error_message
                                      ),
                  pma->error_message,
                  pma->error_message);
      }
      //End coeff
    }
    //End stp
  }
  //End ic
  /**
   * Sampled_source is now filled, thus we don't need the second derivatives anymore
   * */
  free(ddsource);
  return _SUCCESS_;
}
int matter_FFTlog_perturbation_sources_parallel(
                                      struct background * pba,
                                      struct perturbs * ppt,
                                      struct primordial * ppm,
                                      struct matters * pma,
                                      double ** sampled_sources,
                                      double ** prim_spec,
                                      double ** fft_coeff_real,
                                      double ** fft_coeff_imag
                                      ) {
  int N_threads;
  //#ifdef _USE_OMP_
  N_threads = omp_get_max_threads();
  //#else
  //N_threads = 1;
  //#endif
  int n_thread;
  double** integrand1;
  double** integrand2;
  class_alloc(integrand1,
              N_threads*sizeof(double*),
              pma->error_message);
  class_alloc(integrand2,
              N_threads*sizeof(double*),
              pma->error_message);
  for(n_thread=0;n_thread<N_threads;++n_thread){
    class_alloc(integrand1[n_thread],
                pma->size_fft_input*sizeof(double),
                pma->error_message);
    class_alloc(integrand2[n_thread],
                pma->size_fft_input*sizeof(double),
                pma->error_message);
  }
  if(pma->uses_separability){
    if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS){
      printf("Method :: FFT only today\n");
    }
    int index_stp1,index_stp2,index_stp1_stp2;
    int index_ic1,index_ic2,index_ic1_ic2;
    int index_coeff;
    int index_tau1,index_tau2,index_tau1_tau2;
    /**
     * First we allocate the integrands and the fft coefficients
     * */
    for (index_ic1 = 0; index_ic1 < pma->ic_size; index_ic1++) {
      for (index_ic2 = index_ic1; index_ic2 < pma->ic_size; index_ic2++) {
        index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,pma->ic_size);
        for (index_stp1_stp2 = 0; index_stp1_stp2 < pma->stp_grid_size; index_stp1_stp2++) {
          class_alloc(fft_coeff_real[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2],
                2*pma->size_fft_input*sizeof(double),
                pma->error_message);
          class_alloc(fft_coeff_imag[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2],
                2*pma->size_fft_input*sizeof(double),
                pma->error_message);
        }
        //End stp combination
      }
      //End ic2
    }
    //End ic1
    /**
     * Now we iterate over all combinations of sources and do the FFT transformation
     * */
    for (index_ic1 = 0; index_ic1 < pma->ic_size; index_ic1++) {
      for (index_ic2 = index_ic1; index_ic2 < pma->ic_size; index_ic2++) {
        index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,pma->ic_size);
        int abort = _FALSE_;
        #pragma omp parallel for collapse(2) private(index_stp1,index_stp2,index_coeff,index_stp1_stp2,index_tau1,index_tau2,index_tau1_tau2) firstprivate(fft_coeff_real,fft_coeff_imag)
        for (index_stp1 = 0; index_stp1 < pma->stp_size; index_stp1++){
          for (index_stp2 = 0; index_stp2 < pma->stp_size; index_stp2++){
            int tid = omp_get_thread_num();
            index_stp1_stp2 = index_stp1*pma->stp_size+index_stp2;
            index_tau1=pma->tau_size-1;
            index_tau2=pma->tau_size-1;
            index_tau1_tau2 = index_tau1*pma->tau_size+index_tau2;
            /**
             * There is a neat trick with FFT transformations for real inputs,
             * we can do two transformations at once.
             *  In this case, we should simply ignore the second part,
             *  since there is only one tau value (tau0),
             * It was easier to simply set the second integrand to 0 than rewriting
             *  the entire FFT functionality
             * */
            for(index_coeff=0;index_coeff<pma->size_fft_input;++index_coeff){
              integrand1[tid][index_coeff] = sampled_sources[index_ic1*pma->stp_size+index_stp1][index_coeff*pma->tau_size+index_tau1]
                                        *sampled_sources[index_ic2*pma->stp_size+index_stp2][index_coeff*pma->tau_size+index_tau2]
                                        *prim_spec[index_ic1_ic2][index_coeff]
                                        *pow(pma->k_sampling[index_coeff]/pma->k_sampling[0],-pma->bias);
            }
            for(index_coeff=0;index_coeff<pma->size_fft_input;++index_coeff){
              integrand2[tid][index_coeff] = 0.0;
            }
            FFT_real_short(integrand1[tid],
                           integrand2[tid],
                           fft_coeff_real[index_ic1_ic2*pma->stp_size+index_stp1_stp2],
                           fft_coeff_imag[index_ic1_ic2*pma->stp_size+index_stp1_stp2],
                           fft_coeff_real[index_ic1_ic2*pma->stp_size+index_stp1_stp2]+1*pma->size_fft_input,
                           fft_coeff_imag[index_ic1_ic2*pma->stp_size+index_stp1_stp2]+1*pma->size_fft_input,
                           pma->size_fft_input);
            /**
             * The coefficients we have calculated are not yet the final coefficients
             *  For these, we have to multiply with k0^(nu_imag)
             * */
            for(index_coeff=0;index_coeff<pma->size_fft_result;++index_coeff){
              double exp_factor = pow(pma->k_sampling[0],-pma->bias);
              double phase = -log(pma->k_sampling[0])*pma->nu_imag[index_coeff];
              double exp_real = exp_factor*cos(phase);
              double exp_imag = exp_factor*sin(phase);
              double coeff_real = fft_coeff_real[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2][index_coeff];
              double coeff_imag = fft_coeff_imag[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2][index_coeff];
              double newcoeff_real = coeff_real*exp_real-coeff_imag*exp_imag;
              double newcoeff_imag = coeff_real*exp_imag+coeff_imag*exp_real;
              fft_coeff_real[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2][index_coeff] = newcoeff_real/pma->size_fft_input;
              fft_coeff_imag[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2][index_coeff] = newcoeff_imag/pma->size_fft_input;
            }
            //end coeffs
          }
          //End stp2
        }
        //End stp1
        if(abort == _TRUE_){return _FAILURE_;}
      }
      //end ic2
    }
    //end ic1
    if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS){
      printf(" -> Returning from FFTlog... \n");
    }
  }else{
    if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS){
      printf("Method :: FFT for every tau combination\n");
    }
    int index_stp1,index_stp2,index_stp1_stp2;
    int index_ic1,index_ic2,index_ic1_ic2;
    int index_tau1=0,index_tau2=0,index_tau1_tau2=0;
    int index_coeff;
    /**
     * Allocate the integrands and the coefficient arrays
     * */
    for (index_ic1 = 0; index_ic1 < pma->ic_size; index_ic1++) {
      for (index_ic2 = index_ic1; index_ic2 < pma->ic_size; index_ic2++) {
        index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,pma->ic_size);
        for (index_stp1_stp2 = 0; index_stp1_stp2 < pma->stp_grid_size; index_stp1_stp2++) {
          class_alloc(fft_coeff_real[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2],
                pma->size_fft_input*(pma->tau_grid_size+1)*sizeof(double),
                pma->error_message);
          class_alloc(fft_coeff_imag[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2],
                pma->size_fft_input*(pma->tau_grid_size+1)*sizeof(double),
                pma->error_message);
        }
        //End stp combination
      }
      //End ic2
    }
    //End ic1
    /**
     * Now we iterate over all combinations of sources and do the FFT transformation
     * */
    for (index_ic1 = 0; index_ic1 < pma->ic_size; index_ic1++) {
      for (index_ic2 = index_ic1; index_ic2 < pma->ic_size; index_ic2++) {
        index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,pma->ic_size);
        if(pma->matter_verbose > MATTER_VERBOSITY_RANGES){
          printf(" -> At combination of ics %i , %i \n",index_ic1,index_ic2);
        }
        for (index_stp1 = 0; index_stp1 < pma->stp_size; index_stp1++){
          for (index_stp2 = 0; index_stp2 < pma->stp_size; index_stp2++){
            index_stp1_stp2 = index_stp1*pma->stp_size+index_stp2;
            if(pma->matter_verbose > MATTER_VERBOSITY_RANGES){
              printf(" -> At combination of sources %i , %i \n",index_stp1,index_stp2);
            }
            int abort = _FALSE_;
            #pragma omp parallel for collapse(2) private(index_tau1,index_tau2,index_tau1_tau2,index_coeff) firstprivate(fft_coeff_real,fft_coeff_imag)
            for(index_tau1=0;index_tau1<pma->tau_size;++index_tau1){
              for(index_tau2=0;index_tau2<pma->tau_size;index_tau2+=2){
                int tid = omp_get_thread_num();
                index_tau1_tau2 = index_tau1*pma->tau_size+index_tau2;
                /**
                 * There is a neat trick with FFT transformations for real inputs,
                 * we can do two transformations at once.
                 * */
                for(index_coeff=0;index_coeff<pma->size_fft_input;++index_coeff){
                  integrand1[tid][index_coeff] = sampled_sources[index_ic1*pma->stp_size+index_stp1][index_coeff*pma->tau_size+index_tau1]
                                            *sampled_sources[index_ic2*pma->stp_size+index_stp2][index_coeff*pma->tau_size+index_tau2]
                                            *prim_spec[index_ic1_ic2][index_coeff]
                                            *pow(pma->k_sampling[index_coeff]/pma->k_sampling[0],-pma->bias);
                }
                for(index_coeff=0;index_coeff<pma->size_fft_input;++index_coeff){
                  integrand2[tid][index_coeff] = sampled_sources[index_ic1*pma->stp_size+index_stp1][index_coeff*pma->tau_size+index_tau1]
                                            *sampled_sources[index_ic2*pma->stp_size+index_stp2][index_coeff*pma->tau_size+index_tau2+1]
                                            *prim_spec[index_ic1_ic2][index_coeff]
                                            *pow(pma->k_sampling[index_coeff]/pma->k_sampling[0],-pma->bias);
                }
                FFT_real_short(integrand1[tid],
                               integrand2[tid],
                               fft_coeff_real[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2]+index_tau1_tau2*pma->size_fft_input,
                               fft_coeff_imag[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2]+index_tau1_tau2*pma->size_fft_input,
                               fft_coeff_real[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2]+(index_tau1_tau2+1)*pma->size_fft_input,
                               fft_coeff_imag[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2]+(index_tau1_tau2+1)*pma->size_fft_input,
                               pma->size_fft_input);
                /**
                 * The coefficients we have calculated are not yet the final coefficients
                 *  For these, we have to multiply with k0^(nu_imag)
                 * */
                for(index_coeff=0;index_coeff<pma->size_fft_result;++index_coeff){
                  double exp_factor = pow(pma->k_sampling[0],-pma->bias);
                  double phase = -log(pma->k_sampling[0])*pma->nu_imag[index_coeff];
                  double exp_real = exp_factor*cos(phase);
                  double exp_imag = exp_factor*sin(phase);
                  double coeff_real = fft_coeff_real[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2][index_tau1_tau2*pma->size_fft_input+index_coeff];
                  double coeff_imag = fft_coeff_imag[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2][index_tau1_tau2*pma->size_fft_input+index_coeff];
                  double newcoeff_real = coeff_real*exp_real-coeff_imag*exp_imag;
                  double newcoeff_imag = coeff_real*exp_imag+coeff_imag*exp_real;
                  fft_coeff_real[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2][index_tau1_tau2*pma->size_fft_input+index_coeff] = newcoeff_real/pma->size_fft_input;
                  fft_coeff_imag[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2][index_tau1_tau2*pma->size_fft_input+index_coeff] = newcoeff_imag/pma->size_fft_input;
                  coeff_real = fft_coeff_real[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2][(index_tau1_tau2+1)*pma->size_fft_input+index_coeff];
                  coeff_imag = fft_coeff_imag[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2][(index_tau1_tau2+1)*pma->size_fft_input+index_coeff];
                  newcoeff_real = coeff_real*exp_real-coeff_imag*exp_imag;
                  newcoeff_imag = coeff_real*exp_imag+coeff_imag*exp_real;
                  fft_coeff_real[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2][(index_tau1_tau2+1)*pma->size_fft_input+index_coeff] = newcoeff_real/pma->size_fft_input;
                  fft_coeff_imag[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2][(index_tau1_tau2+1)*pma->size_fft_input+index_coeff] = newcoeff_imag/pma->size_fft_input;
                }
                //end coeffs
              }
              //End tau2
            }
            //End tau1
            if(abort==_TRUE_){return _FAILURE_;}
          }
          //End stp1
        }
        //End stp2
      }
      //End ic2
    }
    //End ic1

    if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS ){
      printf(" -> Returning from FFTlog... \n");
    }

  }
  for(n_thread=0;n_thread<N_threads;++n_thread){
    free(integrand1[n_thread]);
    free(integrand2[n_thread]);
  }
  free(integrand1);
  free(integrand2);

  return _SUCCESS_;
}
int matter_integrate_cl(
                        struct precision* ppr,
                        struct background* pba,
                        struct perturbs * ppt,
                        struct matters* pma,
                        double ** fft_coeff_real,
                        double ** fft_coeff_imag){
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS){
    printf("Method :: Integrate Cl's\n");
  }
  /**
   * Initialize and allocate local variables
   * */
  int index_l;
  int index_tilt_grid;
  int index_coeff;
  int index_tw;
  int index_wd1,index_wd2,index_wd1_wd2;
  int index_ic1,index_ic2,index_ic1_ic2;
  int index_cltp;

  int index_radtp1_radtp2;
  double*** window_fft_real;
  double*** window_fft_imag;
  double*** window_bessel_real;
  double*** window_bessel_imag;

  /**
   *
   * Now we allocate the fft coefficient arrays
   *  we use to interpolate into
   * We also allocate the bessel integral arrays
   *  we use to interpolate into
   *
   * */
  int N_threads;
  //#ifdef _USE_OMP_
  N_threads = omp_get_max_threads();
  //#else
  //N_threads = 1;
  //#endif
  class_alloc(window_fft_real,
              N_threads*sizeof(double**),
              pma->error_message);
  class_alloc(window_fft_imag,
              N_threads*sizeof(double**),
              pma->error_message);
  int tw_max_size = MAX(pma->integrated_tw_size,pma->tw_size);
  if(pma->uses_integration == matter_integrate_tw_t || pma->uses_integration == matter_integrate_tw_logt){
    /**
     * First allocalte fft coefficient arrays
     * */
    int n_thread;
    for(n_thread=0;n_thread<N_threads;++n_thread){
      class_alloc(window_fft_real[n_thread],
                  2*tw_max_size*sizeof(double*),
                  pma->error_message);
      class_alloc(window_fft_imag[n_thread],
                  2*tw_max_size*sizeof(double*),
                  pma->error_message);

      if(!pma->uses_separability){
        for(index_tw=0;index_tw<2*tw_max_size;++index_tw){
            class_alloc(window_fft_real[n_thread][index_tw],
                    pma->size_fft_result*sizeof(double),
                    pma->error_message);
            class_alloc(window_fft_imag[n_thread][index_tw],
                    pma->size_fft_result*sizeof(double),
                    pma->error_message);
        }
        //End index_tw
      }
      //Ifend separability
    }
    /**
     * Now allocate bessel arrays
     * */
    class_alloc(window_bessel_real,
              pma->l_size*sizeof(double**),
              pma->error_message);
    class_alloc(window_bessel_imag,
                pma->l_size*sizeof(double**),
                pma->error_message);
    for(index_l=0;index_l<pma->l_size;++index_l){
      class_alloc(window_bessel_real[index_l],
                  pma->tilt_grid_size*pma->size_fft_result*sizeof(double*),
                  pma->error_message);
      class_alloc(window_bessel_imag[index_l],
                  pma->tilt_grid_size*pma->size_fft_result*sizeof(double*),
                  pma->error_message);
      for(index_tilt_grid=0;index_tilt_grid<pma->tilt_grid_size;++index_tilt_grid){
        for(index_coeff=0;index_coeff<pma->size_fft_result;++index_coeff){
          //HERE :: TODO :: check that 2:1 is okay instead of 4:1
          class_alloc(window_bessel_real[index_l][index_tilt_grid*pma->size_fft_result+index_coeff],
                    (pma->has_integrated_windows?2:1)*pma->t_size*sizeof(double),
                    pma->error_message);
          class_alloc(window_bessel_imag[index_l][index_tilt_grid*pma->size_fft_result+index_coeff],
                    (pma->has_integrated_windows?2:1)*pma->t_size*sizeof(double),
                    pma->error_message);
        }
        //End coeff
      }
      //End tilt grid
    }
    //End l
  }
  //Ifend integration type
  else{
    class_stop(pma->error_message,
               "tau integration type not recognized. (Neither tw_t, nor tw_logt) ");
    window_fft_real = NULL;
    window_fft_imag = NULL;
    window_bessel_real = NULL;
    window_bessel_imag = NULL;
  }
  //Ifend known integration method
  if(pma->has_integrated_windows && !pma->uses_limber_approximation){
    int abort = _FALSE_;
    #pragma omp parallel private(index_l) firstprivate(pma,pba,window_bessel_real,window_bessel_imag)
    {
      #pragma omp for
      for(index_l=0;index_l<pma->l_size;++index_l){
        int index_l_eval = (pma->uses_bessel_recursion?pma->l_sampling[index_l]:index_l);
        class_call_parallel(matter_get_bessel_fort_parallel_integrated(
                        pba,
                        pma,
                        0,
                        0,
                        index_l_eval,
                        window_bessel_real[index_l],
                        window_bessel_imag[index_l],
                        (pma->uses_integration == matter_integrate_tw_logt)
                        ),
                pma->error_message,
                pma->error_message);
        //Ifend limber
      }
      //End l
    }
    if(abort == _TRUE_) {return _FAILURE_;}
  }
  if(pma->uses_limber_approximation){
    int abort = _FALSE_;
    #pragma omp parallel private(index_l) firstprivate(pma,pba,window_bessel_real,window_bessel_imag)
    {
      #pragma omp for
      for(index_l=0;index_l<pma->l_size;++index_l){
        class_call_parallel(matter_get_bessel_limber(pma,
                                                  index_l,
                                                  window_bessel_real[index_l],
                                                  window_bessel_imag[index_l]
                                                  ),
                         pma->error_message,
                         pma->error_message);
      }
      //End l
    }
    if(abort == _TRUE_) {return _FAILURE_;}
  }
  /**
   * Finally, we allocate the Cl's array
   *  and iterate through all initial conditions and window functions
   * */
  class_alloc(pma->cl,
              pma->ic_ic_size*pma->cltp_size*sizeof(double*),
              pma->error_message);
  for(index_cltp=0;index_cltp<pma->cltp_size;++index_cltp){
    for(index_ic1=0;index_ic1<pma->ic_size;++index_ic1){
      for(index_ic2=index_ic1;index_ic2<pma->ic_size;++index_ic2){
        index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,pma->ic_size);
        class_alloc(pma->cl[index_ic1_ic2*pma->cltp_size+index_cltp],
                    //(pma->num_windows*(pma->num_windows+1))/2*pma->l_size*sizeof(double),
                    pma->num_window_grid*pma->l_size*sizeof(double),
                    pma->error_message);
        int win_counter = 0;
        for(index_wd1=0;index_wd1<pma->num_windows;++index_wd1){
          for(index_wd2=index_wd1;index_wd2<MIN(pma->num_windows,index_wd1+pma->non_diag+1);++index_wd2){
            index_wd1_wd2 = index_symmetric_matrix(index_wd1,index_wd2,pma->num_windows);
            win_counter++;
            if(pma->matter_verbose > MATTER_VERBOSITY_CLCALCULATION){
            //  printf(" -> Calculating for Window Combination %5d,%5d (%10d/%10d) \n",index_wd1,index_wd2,index_wd1_wd2,(pma->num_windows*(pma->num_windows+1))/2)
              printf(" -> Calculating for Window Combination %5d,%5d (%10d/%10d) \n",index_wd1,index_wd2,win_counter,((pma->non_diag+1)*(2*pma->num_windows-pma->non_diag))/2);
            }
            /**
             * Now allocate local arrays to store the function f_n^{ij}(t) in
             *
             * These can theoretically become quite big,
             *   so we allocate a single one for every window and ic combination
             * */
            class_alloc(pma->intxi_real,
                        pma->radtp_grid_size*sizeof(double*),
                        pma->error_message);
            class_alloc(pma->intxi_imag,
                        pma->radtp_grid_size*sizeof(double*),
                        pma->error_message);
            for(index_radtp1_radtp2=0;index_radtp1_radtp2<pma->radtp_grid_size;++index_radtp1_radtp2){
              class_alloc(pma->intxi_real[index_radtp1_radtp2],
                          pma->size_fft_result*pma->t_size*sizeof(double),
                          pma->error_message);
              class_alloc(pma->intxi_imag[index_radtp1_radtp2],
                          pma->size_fft_result*pma->t_size*sizeof(double),
                          pma->error_message);
            }
            /**
             * If we desire interpolation, those arrays also have to be allocated
             * */
            if(pma->uses_intxi_interpolation){
              class_alloc(pma->intxi_spline_real,
                          pma->radtp_grid_size*sizeof(double*),
                          pma->error_message);
              class_alloc(pma->intxi_spline_imag,
                          pma->radtp_grid_size*sizeof(double*),
                          pma->error_message);
              class_alloc(pma->ddintxi_spline_real,
                          pma->radtp_grid_size*sizeof(double*),
                          pma->error_message);
              class_alloc(pma->ddintxi_spline_imag,
                          pma->radtp_grid_size*sizeof(double*),
                          pma->error_message);
              for(index_radtp1_radtp2=0;index_radtp1_radtp2<pma->radtp_grid_size;++index_radtp1_radtp2){
                class_alloc(pma->intxi_spline_real[index_radtp1_radtp2],
                            pma->size_fft_result*pma->t_spline_size*sizeof(double*),
                            pma->error_message);
                class_alloc(pma->intxi_spline_imag[index_radtp1_radtp2],
                            pma->size_fft_result*pma->t_spline_size*sizeof(double*),
                            pma->error_message);
                class_alloc(pma->ddintxi_spline_real[index_radtp1_radtp2],
                            pma->size_fft_result*pma->t_spline_size*sizeof(double*),
                            pma->error_message);
                class_alloc(pma->ddintxi_spline_imag[index_radtp1_radtp2],
                            pma->size_fft_result*pma->t_spline_size*sizeof(double*),
                            pma->error_message);
              }
            }
            /**
             * Now integrate the actual integral
             * */
            if(pma->uses_integration == matter_integrate_tw_t  || pma->uses_integration == matter_integrate_tw_logt){
              class_call(matter_integrate_for_each_ttau_parallel_chi_pre(ppr,
                                                        pba,
                                                        ppt,
                                                        pma,
                                                        fft_coeff_real,
                                                        fft_coeff_imag,
                                                        index_ic1,
                                                        index_ic2,
                                                        index_ic1_ic2,
                                                        index_wd1,
                                                        index_wd2,
                                                        index_cltp,
                                                        window_fft_real,
                                                        window_fft_imag,
                                                        window_bessel_real,
                                                        window_bessel_imag,
                                                        (pma->uses_integration == matter_integrate_tw_logt)
                                                        ),
                         pma->error_message,
                         pma->error_message);
            }
            else{
              class_stop(pma->error_message,
                         "Unrecognized integration method \n");
            }
            /**
             * Finally delete the intxi storage arrays again
             * */
            if(pma->uses_intxi_interpolation){
              for(index_radtp1_radtp2=0;index_radtp1_radtp2<pma->radtp_grid_size;++index_radtp1_radtp2){
                free(pma->intxi_spline_real[index_radtp1_radtp2]);
                free(pma->intxi_spline_imag[index_radtp1_radtp2]);
                free(pma->ddintxi_spline_real[index_radtp1_radtp2]);
                free(pma->ddintxi_spline_imag[index_radtp1_radtp2]);
              }
              free(pma->intxi_spline_real);
              free(pma->intxi_spline_imag);
              free(pma->ddintxi_spline_real);
              free(pma->ddintxi_spline_imag);
            }
            for(index_radtp1_radtp2=0;index_radtp1_radtp2<pma->radtp_grid_size;++index_radtp1_radtp2){
              free(pma->intxi_real[index_radtp1_radtp2]);
              free(pma->intxi_imag[index_radtp1_radtp2]);
            }
            free(pma->intxi_real);
            free(pma->intxi_imag);
          }
          //End wd2
        }
        //End wd1
      }
      //End ic2
    }
    //End ic1
  }
  //End cltp
  /**
   * Finally print the obtained results
   * */
  if(pma->matter_verbose > MATTER_VERBOSITY_CLRESULTS){
    /**
     * Print the direct C_l's
     * */
    printf("RESULTS C_l = \n\n");
    printf(" -> l sampling : \n");
    for(index_l=0;index_l<pma->l_size;++index_l){
      printf("%i,",(int)pma->l_sampling[index_l]);
    }
    printf("\n");
    for(index_cltp=0;index_cltp<pma->cltp_size;++index_cltp){
      printf(" -> At cltp %4d \n",index_cltp);
      for(index_ic1_ic2=0;index_ic1_ic2<pma->ic_ic_size;++index_ic1_ic2){
        printf("    ->At icic %4d \n",index_ic1_ic2);
        for(index_wd1=0;index_wd1<pma->num_windows;++index_wd1){
          for(index_wd2=index_wd1;index_wd2<MIN(pma->num_windows,index_wd1+pma->non_diag+1);++index_wd2){
            index_wd1_wd2 = index_symmetric_matrix(index_wd1,index_wd2,pma->num_windows);
            printf(" -> At win (%4d,%4d) ... \n",index_wd1,index_wd2);
            printf("%.10e",
              pma->cl[index_ic1_ic2*pma->cltp_size+index_cltp][index_wd1_wd2*pma->l_size+0]
            );
            for(index_l=1;index_l<pma->l_size;++index_l){
              printf(",%.10e",
                pma->cl[index_ic1_ic2*pma->cltp_size+index_cltp][index_wd1_wd2*pma->l_size+index_l]
              );
            }
            printf("\n");
            //End l
          }
          //End wd2
        }
        //End wd1
      }
      //End icgrid
    }
    //End cltp
    printf("\n");
    /**
     * Now also print the l(l+1)/2pi C_l's
     * */
    printf("RESULTS l(l+1)/2pi C_l = \n\n");
    printf(" -> l sampling : \n");
    for(index_l=0;index_l<pma->l_size;++index_l){
      printf("%i,", (int)pma->l_sampling[index_l]);
    }
    printf("\n");
    for(index_cltp=0;index_cltp<pma->cltp_size;++index_cltp){
      printf(" -> At cltp %4d \n",index_cltp);
      for(index_ic1_ic2=0;index_ic1_ic2<pma->ic_ic_size;++index_ic1_ic2){
        printf("   -> At icic %4d \n",index_ic1_ic2);
        for(index_wd1=0;index_wd1<pma->num_windows;++index_wd1){
          for(index_wd2=index_wd1;index_wd2<MIN(pma->num_windows,index_wd1+pma->non_diag+1);++index_wd2){
            index_wd1_wd2 = index_symmetric_matrix(index_wd1,index_wd2,pma->num_windows);
            printf(" -> At win (%4d,%4d) ... \n",index_wd1,index_wd2);
            printf("%.10e",
              pma->l_sampling[0]*(pma->l_sampling[0]+1.0)/(_TWOPI_)*pma->cl[index_ic1_ic2*pma->cltp_size+index_cltp][index_wd1_wd2*pma->l_size+0]
            );
            for(index_l=1;index_l<pma->l_size;++index_l){
              printf(",%.10e",
                pma->l_sampling[index_l]*(pma->l_sampling[index_l]+1.0)/(_TWOPI_)*pma->cl[index_ic1_ic2*pma->cltp_size+index_cltp][index_wd1_wd2*pma->l_size+index_l]
              );
            }
            printf("\n");
            //End l
          }
          //End wd2
        }
        //End wd1
      }
      //End icgrid
    }
    //End cltp
    printf("\n\n");
  }
  //Ifend Cl printing
  /**
   * Finalyl delete also the temporary arrays
   *  for the coefficients and the bessel functions
   * */
  int n_thread;
  for(n_thread=0;n_thread<N_threads;++n_thread){
    if(!pma->uses_separability){
      for(index_tw=0;index_tw<2*tw_max_size;++index_tw){
        free(window_fft_real[n_thread][index_tw]);
        free(window_fft_imag[n_thread][index_tw]);
      }
      //End tw
    }
    //Ifend sep
    free(window_fft_real[n_thread]);
    free(window_fft_imag[n_thread]);
  }
  //End n_thread
  for(index_l=0;index_l<pma->l_size;++index_l){
    for(index_tilt_grid=0;index_tilt_grid<pma->tilt_grid_size;++index_tilt_grid){
      for(index_coeff=0;index_coeff<pma->size_fft_result;++index_coeff){
        free(window_bessel_real[index_l][index_tilt_grid*pma->size_fft_result+index_coeff]);
        free(window_bessel_imag[index_l][index_tilt_grid*pma->size_fft_result+index_coeff]);
      }
      //End coeff
    }
    //End tilt grid
    free(window_bessel_real[index_l]);
    free(window_bessel_imag[index_l]);
  }
  //End l
  free(window_bessel_real);
  free(window_bessel_imag);
  free(window_fft_real);
  free(window_fft_imag);
  return _SUCCESS_;
}


int matter_integrate_window_function(
                                    struct background* pba,
                                    struct matters* pma,
                                    double* window,
                                    double* integrated_window,
                                    double* integrated_sampling,
                                    double* oldtw_sampling,
                                    double* oldtw_weights,
                                    int integrated_size,
                                    int oldtw_size,
                                    int index_wd,
                                    int radtp,
                                    double* f_evo
                                    ){
  /**
   * A simple routine to integrate the window function
   * */
  int index_integrated,index_oldtw;
  double newsum;
  double val;
  double* bg_factor = NULL;
  /**
   * First we set the last_index to some reasonable value
   * */
  int last_index = 0;
  class_call(background_at_tau(pba,
                               oldtw_sampling[0],
                               pba->short_info,
                               pba->inter_normal,
                               &last_index,
                               pma->short_pvecback),
               pba->error_message,
               pma->error_message);
  /**
   * If there is some background values involved,
   *  get them now before the double loop
   * */
  if(matter_is_index(radtp,pma->radtp_g5,pma->has_gravitational_terms)){
    class_alloc(bg_factor,
                oldtw_size*sizeof(double),
                pma->error_message);
    for(index_oldtw=0;index_oldtw<oldtw_size;++index_oldtw){
      class_call(background_at_tau(pba,
                                   oldtw_sampling[index_oldtw],
                                   pba->short_info,
                                   pba->inter_closeby, //pba->inter_closeby
                                   &last_index,
                                   pma->short_pvecback),
               pba->error_message,
               pma->error_message);
      bg_factor[index_oldtw] = 2.0*(1.0 + pma->short_pvecback[pba->index_bg_H_prime]
                    /pma->short_pvecback[pba->index_bg_a] /pma->short_pvecback[pba->index_bg_H] /pma->short_pvecback[pba->index_bg_H]
                    + (2.0 - 5.0 *pma->selection_magnification_bias[index_wd])
                    /(pma->tau0-oldtw_sampling[index_oldtw])/pma->short_pvecback[pba->index_bg_a] /pma->short_pvecback[pba->index_bg_H]
                    + 5.0 * pma->selection_magnification_bias[index_wd]
                    - f_evo[index_wd*oldtw_size+index_oldtw]);
    }
  }
  for(index_integrated=0;index_integrated<integrated_size;++index_integrated){
    /**
     * Now we integrate the window function,
     *  where we specifically ignore the last value
     * (xi = 0, divergent as 1/xi^2, but regularized by a factor of xi^2)
     * */
    if(index_integrated == integrated_size-1){
      newsum = 0.0;
    }
    else{
      newsum = 0.0;
      for(index_oldtw=0;index_oldtw<oldtw_size;++index_oldtw){
        //TODO :: see where val = 0.0
        if(oldtw_sampling[index_oldtw]>=integrated_sampling[index_integrated]){
          val=0.0;
        }
        else{
          if(fabs(integrated_sampling[index_integrated]-pma->tau0)<2*pma->small_log_offset){
            val=0.0;
          }
          else{
            val = oldtw_weights[index_oldtw]*window[index_oldtw];
            if(matter_is_index(radtp,pma->radtp_nclens,pma->has_lensing_terms)){
              val*=(1.0-2.5*pma->selection_magnification_bias[index_wd])*(integrated_sampling[index_integrated]-oldtw_sampling[index_oldtw])/((pma->tau0-oldtw_sampling[index_oldtw])*(pma->tau0-integrated_sampling[index_integrated]));
            }
            else if(matter_is_index(radtp,pma->radtp_g4,pma->has_gravitational_terms)){
              val*=(2.-5.*pma->selection_magnification_bias[index_wd])/(pma->tau0-oldtw_sampling[index_oldtw]);
            }
            else if(matter_is_index(radtp,pma->radtp_g5,pma->has_gravitational_terms)){
              /**
               * Additional factor of two, since we do the approximation
               * (phi+psi)' = 2 phi'
               * */
              val*=bg_factor[index_oldtw];
            }
            else if(matter_is_index(radtp,pma->radtp_shlens,pma->has_cl_shear)){
              val*=(integrated_sampling[index_integrated]-oldtw_sampling[index_oldtw])/((pma->tau0-oldtw_sampling[index_oldtw])*(pma->tau0-integrated_sampling[index_integrated]));
            }
            else{
              class_stop(pma->error_message,
                         "unrecognized integrated (lensing/gravitational) radial type");
            }
            //Ifend select radial type
          }
          //Ifend too close to tau0
        }
        //Ifend tau'>tau
        newsum+=val;
      }
      //End integration value
    }
    //Ifend not-regularized value
    integrated_window[index_integrated] = newsum;
  }
  //End final integration limit
  if(matter_is_index(radtp,pma->radtp_g5,pma->has_gravitational_terms)){
    free(bg_factor);
  }
  return _SUCCESS_;
}

int matter_get_half_integrand(
                             struct background* pba,
                             struct perturbs* ppt,
                             struct matters* pma,
                             double t,
                             int index_ic1,
                             int index_ic2,
                             int index_radtp1,
                             int index_radtp2,
                             int index_tilt1_tilt2,
                             int index_ic1_ic2,
                             int index_stp1_stp2,
                             int index_cltp,
                             int index_wd1,
                             int index_wd2,
                             double* integrand_real,
                             double* integrand_imag,
                             double** fft_coeff_real,
                             double** fft_coeff_imag,
                             double** wint_fft_real,
                             double** wint_fft_imag,
                             double* tw_local_sampling,
                             int tw_local_size
                             ){
  int index_tw_local;
  double* fft_real = fft_coeff_real[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2];
  double* fft_imag = fft_coeff_imag[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2];
  int is_integrated = matter_is_integrated(index_radtp1) && matter_is_integrated(index_radtp2);

  double window0_val,window1_val;
  int inf0=0,inf1=0;
  double x0,x1;

  int fft_index00,fft_index01,fft_index10,fft_index11;
  int index_coeff;
  short x1flag;
  double h0,a0,b0;
  double h1,a1,b1;
  double *fft_00_ptr,*fft_01_ptr,*fft_10_ptr,*fft_11_ptr;

  class_call(matter_spline_prepare_hunt(
                                pma->tau_sampling,
                                pma->tau_size,
                                tw_local_sampling[index_wd1*tw_local_size],
                                &inf0,
                                pma->error_message
                               ),
            pma->error_message,
            pma->error_message);
  class_call(matter_spline_prepare_hunt(
                                pma->tau_sampling,
                                pma->tau_size,
                                pma->tau0*(1-t)+t*tw_local_sampling[index_wd1*tw_local_size],
                                &inf1,
                                pma->error_message
                               ),
            pma->error_message,
            pma->error_message);
  for(index_tw_local=0;index_tw_local<tw_local_size;++index_tw_local){

    x1flag = _TRUE_;

    if(pma->uses_intxi_logarithmic && is_integrated){
      x0 = pma->tau0-pma->exp_integrated_tw_sampling[index_wd1*tw_local_size+index_tw_local];//exp(tw_local_sampling[index_wd1*tw_local_size+index_tw_local]);
      x1 = pma->tau0-t*pma->exp_integrated_tw_sampling[index_wd1*tw_local_size+index_tw_local];// - pma->small_log_offset;
    }
    else{
      x0 = tw_local_sampling[index_wd1*tw_local_size+index_tw_local];
      x1 = pma->tau0-t*(pma->tau0-x0);// - pma->small_log_offset;
    }
    class_test(x0>pma->tau0,
               pma->error_message,
               "with x0 = %.10e , t = %.10e ,tau0 = %.10e , x0 = %.10e",x0,t,pma->tau0,x0);
    class_test(x1>pma->tau0,
               pma->error_message,
               "with x0 = %.10e , t = %.10e ,tau0 = %.10e , x1 = %.10e",x0,t,pma->tau0,x1);
    if(
      (x1>pma->tw_max[index_wd2] && (!(pma->has_integrated_windows && matter_is_integrated(index_radtp2))))
      ||x1<pma->tw_min[index_wd2]){
      //The point x1 is outside of the window w2
      x1flag=_FALSE_;
    }
    class_call(matter_spline_hunt(
                                  pma->tau_sampling,
                                  pma->tau_size,
                                  x0,//285+pma->small_log_offset,
                                  &inf0,
                                  &h0,
                                  &a0,
                                  &b0,
                                  pma->error_message
                                 ),
               pma->error_message,
               pma->error_message);
    if(x1flag==_TRUE_){
     class_call(matter_spline_hunt(
                                  pma->tau_sampling,
                                  pma->tau_size,
                                  x1,//285+pma->small_log_offset,
                                  &inf1,
                                  &h1,
                                  &a1,
                                  &b1,
                                  pma->error_message
                                 ),
               pma->error_message,
               pma->error_message);
    }
    if(!pma->uses_separability){
      if(x1flag==_TRUE_){
        fft_index00 = (inf0)*pma->tau_size+(inf1);
        fft_index01 = (inf0)*pma->tau_size+(inf1+1);
        fft_index10 = (inf0+1)*pma->tau_size+(inf1);
        fft_index11 = (inf0+1)*pma->tau_size+(inf1+1);
        fft_00_ptr = fft_real+fft_index00*pma->size_fft_input;
        fft_01_ptr = fft_real+fft_index01*pma->size_fft_input;
        fft_10_ptr = fft_real+fft_index10*pma->size_fft_input;
        fft_11_ptr = fft_real+fft_index11*pma->size_fft_input;
        for(index_coeff=0;index_coeff<pma->size_fft_result;++index_coeff){
          wint_fft_real[index_tw_local][index_coeff] = a0*a1*fft_00_ptr[index_coeff]+a0*b1*fft_01_ptr[index_coeff]+b0*a1*fft_10_ptr[index_coeff]+b0*b1*fft_11_ptr[index_coeff];
        }
        fft_00_ptr = fft_imag+fft_index00*pma->size_fft_input;
        fft_01_ptr = fft_imag+fft_index01*pma->size_fft_input;
        fft_10_ptr = fft_imag+fft_index10*pma->size_fft_input;
        fft_11_ptr = fft_imag+fft_index11*pma->size_fft_input;
        for(index_coeff=0;index_coeff<pma->size_fft_result;++index_coeff){
          wint_fft_imag[index_tw_local][index_coeff] = a0*(a1*fft_00_ptr[index_coeff]+b1*fft_01_ptr[index_coeff])+b0*(a1*fft_10_ptr[index_coeff]+b1*fft_11_ptr[index_coeff]);
        }
      }
      else{
        memset(wint_fft_real[index_tw_local],0,pma->size_fft_result*sizeof(double));
        memset(wint_fft_imag[index_tw_local],0,pma->size_fft_result*sizeof(double));
      }
      //End if x1
    }else{
      wint_fft_real[index_tw_local]=fft_real;
      wint_fft_imag[index_tw_local]=fft_imag;
      wint_fft_real[tw_local_size+index_tw_local]=fft_real;
      wint_fft_imag[tw_local_size+index_tw_local]=fft_imag;
    }
  }
  //End tw
  int last_index0,last_index1;
  class_call(matter_spline_prepare_hunt(
                              tw_local_sampling+index_wd1*tw_local_size,
                              tw_local_size,
                              tw_local_sampling[index_wd1*tw_local_size],
                              &last_index0,
                              pma->error_message),
            pma->error_message,
            pma->error_message);
  class_call(matter_spline_prepare_hunt(
                              tw_local_sampling+index_wd1*tw_local_size,
                              tw_local_size,
                              pma->tau0*(1-t)+t*tw_local_sampling[index_wd1*tw_local_size],
                              &last_index1,
                              pma->error_message),
            pma->error_message,
            pma->error_message);

  int derivative_type1 = 0;
  int derivative_type2 = 0;
  class_call(matter_get_derivative_type(pma,
                             &derivative_type1,
                             &derivative_type2,
                             index_radtp1,
                             index_radtp2),
             pma->error_message,
             pma->error_message);
  for(index_tw_local=0;index_tw_local<tw_local_size;++index_tw_local){
    if(is_integrated && pma->uses_intxi_logarithmic){
      x0 = pma->tau0-pma->exp_integrated_tw_sampling[index_wd1*tw_local_size+index_tw_local];//exp(tw_local_sampling[index_wd1*tw_local_size+index_tw_local]);
    }
    else{
      x0 = tw_local_sampling[index_wd1*tw_local_size+index_tw_local];
    }
    x1 = pma->tau0-t*(pma->tau0-x0);
    class_call(matter_get_prepared_window_at(pma,
                                             x0,
                                             index_ic1,
                                             index_radtp1,
                                             index_wd1,
                                             &last_index0,
                                             derivative_type1,
                                             &window0_val),
              pma->error_message,
              pma->error_message);
    class_call(matter_get_prepared_window_at(pma,
                                             x1,
                                             index_ic2,
                                             index_radtp2,
                                             index_wd2,
                                             &last_index1,
                                             derivative_type2,
                                             &window1_val),
              pma->error_message,
              pma->error_message);
    double wwval = window0_val*window1_val;
    for(index_coeff=0;index_coeff<pma->size_fft_result;++index_coeff){
      integrand_real[index_coeff*tw_local_size+index_tw_local] = wwval*wint_fft_real[index_tw_local][index_coeff];
      integrand_imag[index_coeff*tw_local_size+index_tw_local] = wwval*wint_fft_imag[index_tw_local][index_coeff];
    }
    //End coeff
  }
  //End tw
  return _SUCCESS_;
}
int matter_get_ttau_integrand(
                             struct background* pba,
                             struct perturbs* ppt,
                             struct matters* pma,
                             double t,
                             int index_ic1,
                             int index_ic2,
                             int index_radtp1,
                             int index_radtp2,
                             int index_tilt1_tilt2,
                             int index_ic1_ic2,
                             int index_stp1_stp2,
                             int index_cltp,
                             int index_wd1,
                             int index_wd2,
                             double* integrand_real,
                             double* integrand_imag,
                             double** fft_coeff_real,
                             double** fft_coeff_imag,
                             double** wint_fft_real,
                             double** wint_fft_imag,
                             double* tw_local_sampling,
                             int tw_local_size
                             ){
  int index_tw_local;
  double* fft_real = fft_coeff_real[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2];
  double* fft_imag = fft_coeff_imag[index_ic1_ic2*pma->stp_grid_size+index_stp1_stp2];
  int is_integrated = matter_is_integrated(index_radtp1) && matter_is_integrated(index_radtp2);
  double logt = log(t);
  double window0_val,window1_val,window2_val;
  int inf0=0,inf1=0,inf2=0;
  double x0,x1,x2;
  double exp_factor = exp(logt*(pma->nu_real[index_tilt1_tilt2]-2.0));

  int fft_index00,fft_index01,fft_index10,fft_index11;
  int index_coeff;
  short x1flag,x2flag;
  double h0,a0,b0;
  double h1,a1,b1;
  double h2,a2,b2;
  double* cos_val;
  double* sin_val;
  //Theoretically, these could be allocated outside of the t loop,
  //but currently their allocation is not time consuming at all.
  class_alloc(cos_val,
              pma->size_fft_cutoff*sizeof(double),
              pma->error_message);
  class_alloc(sin_val,
              pma->size_fft_cutoff*sizeof(double),
              pma->error_message);
  for(index_coeff=0;index_coeff<pma->size_fft_cutoff;++index_coeff){
    double phase = logt*pma->nu_imag[index_coeff];
    cos_val[index_coeff] = cos(phase);
    sin_val[index_coeff] = sin(phase);
  }
  double *fft_00_ptr,*fft_01_ptr,*fft_10_ptr,*fft_11_ptr;
  class_call(matter_spline_prepare_hunt(
                                pma->tau_sampling,
                                pma->tau_size,
                                tw_local_sampling[index_wd1*tw_local_size],
                                &inf0,
                                pma->error_message
                               ),
            pma->error_message,
            pma->error_message);
  class_call(matter_spline_prepare_hunt(
                                pma->tau_sampling,
                                pma->tau_size,
                                pma->tau0*(1-t)+t*tw_local_sampling[index_wd1*tw_local_size],
                                &inf1,
                                pma->error_message
                               ),
            pma->error_message,
            pma->error_message);
  class_call(matter_spline_prepare_hunt(
                                pma->tau_sampling,
                                pma->tau_size,
                                pma->tau0*(1-1.0/t)+(1.0/t)*tw_local_sampling[index_wd1*tw_local_size],
                                &inf2,
                                pma->error_message
                               ),
            pma->error_message,
            pma->error_message);

  for(index_tw_local=0;index_tw_local<tw_local_size;++index_tw_local){

    x1flag = _TRUE_;
    x2flag = _TRUE_;

    if(pma->uses_intxi_logarithmic && is_integrated){
      x0 = pma->tau0-pma->exp_integrated_tw_sampling[index_wd1*tw_local_size+index_tw_local];
    }else{
      x0 = tw_local_sampling[index_wd1*tw_local_size+index_tw_local];
    }
    x1 = pma->tau0-t*(pma->tau0-x0);// - pma->small_log_offset;
    x2 = pma->tau0-(1.0/t)*(pma->tau0-x0);// - pma->small_log_offset;
    class_test(x0>pma->tau0,
               pma->error_message,
               "with x0 = %.10e , t = %.10e ,tau0 = %.10e , x0 = %.10e",x0,t,pma->tau0,x0);
    class_test(x1>pma->tau0,
               pma->error_message,
               "with x0 = %.10e , t = %.10e ,tau0 = %.10e , x1 = %.10e",x0,t,pma->tau0,x1);
    class_test(x2>pma->tau0,
               pma->error_message,
               "with x0 = %.10e , t = %.10e ,tau0 = %.10e , x2 = %.10e",x0,t,pma->tau0,x2);
    if(
      (x1>pma->tw_max[index_wd2] && (!(pma->has_integrated_windows && matter_is_integrated(index_radtp2))))
      ||x1<pma->tw_min[index_wd2]){
      //The point x1 is outside of the window w2
      x1flag=_FALSE_;
    }
    if(
      (x2>pma->tw_max[index_wd2] && (!(pma->has_integrated_windows && matter_is_integrated(index_radtp2))))
      ||x2<pma->tw_min[index_wd2]){
      //The point x2 is outside of the window w2
      x2flag=_FALSE_;
    }
    class_call(matter_spline_hunt(
                                  pma->tau_sampling,
                                  pma->tau_size,
                                  x0,
                                  &inf0,
                                  &h0,
                                  &a0,
                                  &b0,
                                  pma->error_message
                                 ),
               pma->error_message,
               pma->error_message);
    if(x1flag==_TRUE_){
       class_call(matter_spline_hunt(
                                    pma->tau_sampling,
                                    pma->tau_size,
                                    x1,
                                    &inf1,
                                    &h1,
                                    &a1,
                                    &b1,
                                    pma->error_message
                                   ),
                 pma->error_message,
                 pma->error_message);
    }
    if(x2flag==_TRUE_){
        class_call(matter_spline_hunt(
                                      pma->tau_sampling,
                                      pma->tau_size,
                                      x2,
                                      &inf2,
                                      &h2,
                                      &a2,
                                      &b2,
                                      pma->error_message
                                     ),
                   pma->error_message,
                   pma->error_message);
    }
    if(!pma->uses_separability){
      if(x1flag==_TRUE_){
        fft_index00 = (inf0)*pma->tau_size+(inf1);
        fft_index01 = (inf0)*pma->tau_size+(inf1+1);
        fft_index10 = (inf0+1)*pma->tau_size+(inf1);
        fft_index11 = (inf0+1)*pma->tau_size+(inf1+1);
        fft_00_ptr = fft_real+fft_index00*pma->size_fft_input;
        fft_01_ptr = fft_real+fft_index01*pma->size_fft_input;
        fft_10_ptr = fft_real+fft_index10*pma->size_fft_input;
        fft_11_ptr = fft_real+fft_index11*pma->size_fft_input;
        for(index_coeff=0;index_coeff<pma->size_fft_result;++index_coeff){
          wint_fft_real[index_tw_local][index_coeff] = a0*a1*fft_00_ptr[index_coeff]+a0*b1*fft_01_ptr[index_coeff]+b0*a1*fft_10_ptr[index_coeff]+b0*b1*fft_11_ptr[index_coeff];
        }
        fft_00_ptr = fft_imag+fft_index00*pma->size_fft_input;
        fft_01_ptr = fft_imag+fft_index01*pma->size_fft_input;
        fft_10_ptr = fft_imag+fft_index10*pma->size_fft_input;
        fft_11_ptr = fft_imag+fft_index11*pma->size_fft_input;
        for(index_coeff=0;index_coeff<pma->size_fft_result;++index_coeff){
          wint_fft_imag[index_tw_local][index_coeff] = a0*(a1*fft_00_ptr[index_coeff]+b1*fft_01_ptr[index_coeff])+b0*(a1*fft_10_ptr[index_coeff]+b1*fft_11_ptr[index_coeff]);
        }
      }
      else{
        memset(wint_fft_real[index_tw_local],0,pma->size_fft_result*sizeof(double));
        memset(wint_fft_imag[index_tw_local],0,pma->size_fft_result*sizeof(double));
      }
      //End if x1
      if(x2flag==_TRUE_){
        fft_index00 = (inf0)*pma->tau_size+(inf2);
        fft_index01 = (inf0)*pma->tau_size+(inf2+1);
        fft_index10 = (inf0+1)*pma->tau_size+(inf2);
        fft_index11 = (inf0+1)*pma->tau_size+(inf2+1);
        fft_00_ptr = fft_real+fft_index00*pma->size_fft_input;
        fft_01_ptr = fft_real+fft_index01*pma->size_fft_input;
        fft_10_ptr = fft_real+fft_index10*pma->size_fft_input;
        fft_11_ptr = fft_real+fft_index11*pma->size_fft_input;
        for(index_coeff=0;index_coeff<pma->size_fft_result;++index_coeff){
          wint_fft_real[index_tw_local+tw_local_size][index_coeff] = a0*a2*fft_00_ptr[index_coeff]+a0*b2*fft_01_ptr[index_coeff]+b0*a2*fft_10_ptr[index_coeff]+b0*b2*fft_11_ptr[index_coeff];
        }
        fft_00_ptr = fft_imag+fft_index00*pma->size_fft_input;
        fft_01_ptr = fft_imag+fft_index01*pma->size_fft_input;
        fft_10_ptr = fft_imag+fft_index10*pma->size_fft_input;
        fft_11_ptr = fft_imag+fft_index11*pma->size_fft_input;
        for(index_coeff=0;index_coeff<pma->size_fft_result;++index_coeff){
          wint_fft_imag[index_tw_local+tw_local_size][index_coeff] = a0*(a2*fft_00_ptr[index_coeff]+b2*fft_01_ptr[index_coeff])+b0*(a2*fft_10_ptr[index_coeff]+b2*fft_11_ptr[index_coeff]);
        }
      }
      else{
        memset(wint_fft_real[index_tw_local+tw_local_size],0,pma->size_fft_result*sizeof(double));
        memset(wint_fft_imag[index_tw_local+tw_local_size],0,pma->size_fft_result*sizeof(double));
      }
      //End if x2
    }else{
      wint_fft_real[index_tw_local]=fft_real;
      wint_fft_imag[index_tw_local]=fft_imag;
      wint_fft_real[tw_local_size+index_tw_local]=fft_real;
      wint_fft_imag[tw_local_size+index_tw_local]=fft_imag;
    }
  }
  //End tw
  int last_index0,last_index1,last_index2;
  class_call(matter_spline_prepare_hunt(
                              tw_local_sampling+index_wd1*tw_local_size,
                              tw_local_size,
                              tw_local_sampling[index_wd1*tw_local_size],
                              &last_index0,
                              pma->error_message),
            pma->error_message,
            pma->error_message);
  class_call(matter_spline_prepare_hunt(
                              tw_local_sampling+index_wd1*tw_local_size,
                              tw_local_size,
                              pma->tau0*(1-t)+t*tw_local_sampling[index_wd1*tw_local_size],
                              &last_index1,
                              pma->error_message),
            pma->error_message,
            pma->error_message);
  class_call(matter_spline_prepare_hunt(
                              tw_local_sampling+index_wd1*tw_local_size,
                              tw_local_size,
                              pma->tau0*(1-1.0/t)+(1.0/t)*tw_local_sampling[index_wd1*tw_local_size],
                              &last_index2,
                              pma->error_message),
             pma->error_message,
             pma->error_message);

  int derivative_type1 = 0;
  int derivative_type2 = 0;
  class_call(matter_get_derivative_type(pma,
                             &derivative_type1,
                             &derivative_type2,
                             index_radtp1,
                             index_radtp2),
             pma->error_message,
             pma->error_message);
  for(index_tw_local=0;index_tw_local<tw_local_size;++index_tw_local){
    if(pma->uses_intxi_logarithmic && is_integrated){
      x0 = pma->tau0-pma->exp_integrated_tw_sampling[index_wd1*tw_local_size+index_tw_local];
    }else{
      x0 = tw_local_sampling[index_wd1*tw_local_size+index_tw_local];
    }
    x1 = pma->tau0*(1-t)+t*x0;
    x2 = pma->tau0*(1-1.0/t)+(1.0/t)*x0;
    class_call(matter_get_prepared_window_at(pma,
                                             x0,
                                             index_ic1,
                                             index_radtp1,
                                             index_wd1,
                                             &last_index0,
                                             derivative_type1,
                                             &window0_val),
              pma->error_message,
              pma->error_message);
    class_call(matter_get_prepared_window_at(pma,
                                             x1,
                                             index_ic2,
                                             index_radtp2,
                                             index_wd2,
                                             &last_index1,
                                             derivative_type2,
                                             &window1_val),
              pma->error_message,
              pma->error_message);
    class_call(matter_get_prepared_window_at(pma,
                                             x2,
                                             index_ic2,
                                             index_radtp2,
                                             index_wd2,
                                             &last_index2,
                                             derivative_type2,
                                             &window2_val),
              pma->error_message,
              pma->error_message);
    double temp_first,temp_second;
    temp_first = window0_val*window1_val;
    temp_second = exp_factor*window0_val*window2_val;
    for(index_coeff=0;index_coeff<pma->size_fft_result;++index_coeff){//was cutoff?
      integrand_real[index_coeff*tw_local_size+index_tw_local] = temp_first*wint_fft_real[index_tw_local][index_coeff]+temp_second*(wint_fft_real[index_tw_local+tw_local_size][index_coeff]*cos_val[index_coeff]-wint_fft_imag[index_tw_local+tw_local_size][index_coeff]*sin_val[index_coeff]);
      integrand_imag[index_coeff*tw_local_size+index_tw_local] = temp_first*wint_fft_imag[index_tw_local][index_coeff]+temp_second*(wint_fft_imag[index_tw_local+tw_local_size][index_coeff]*cos_val[index_coeff]+wint_fft_real[index_tw_local+tw_local_size][index_coeff]*sin_val[index_coeff]);
    }
    //End coeff
  }
  //End tw
  free(cos_val);
  free(sin_val);
  return _SUCCESS_;
}
int matter_asymptote(struct precision* ppr, struct matters* pma,double t, int index_wd1, int index_wd2,double* result){
  double x1 = pma->tau0-0.5*(pma->tw_max[index_wd1]+pma->tw_min[index_wd1]);
  double x2 = pma->tau0-0.5*(pma->tw_max[index_wd2]+pma->tw_min[index_wd2]);

  double sigma1 = 0.5*(pma->tw_max[index_wd1]-pma->tw_min[index_wd1])/ppr->selection_cut_at_sigma;
  double sigma2 = 0.5*(pma->tw_max[index_wd2]-pma->tw_min[index_wd2])/ppr->selection_cut_at_sigma;
  *result = exp(-0.5*(x1*t-x2)*(x1*t-x2)/(sigma1*sigma1*t*t+sigma2*sigma2))+exp(-0.5*(x2*t-x1)*(x2*t-x1)/(sigma2*sigma2*t*t+sigma1*sigma1));
  return _SUCCESS_;
}
int matter_precompute_chit_factors(struct matters* pma,
                                   double* tw_local_sampling,
                                   int tw_local_size,
                                   int index_tilt1_tilt2,
                                   int index_wd,
                                   short is_log,
                                   double* pref_real,
                                   double* pref_imag){
  int index_tw_local;
  double logxi;
  double exp_factor,phase;
  int index_coeff;
  double offset = 1.0;
  if(is_log == _TRUE_){
    for(index_tw_local=0;index_tw_local<tw_local_size;++index_tw_local){
      logxi = tw_local_sampling[index_wd*tw_local_size+index_tw_local];
      exp_factor = exp(logxi*(2.0-pma->nu_real[index_tilt1_tilt2]));
      //exp_factor = exp(logxi*(2.0-pma->bias));
      for(index_coeff=0;index_coeff<pma->size_fft_cutoff;++index_coeff){
        phase = -logxi*pma->nu_imag[index_coeff];
        pref_real[index_tw_local*pma->size_fft_result+index_coeff] = exp_factor*cos(phase);
        pref_imag[index_tw_local*pma->size_fft_result+index_coeff] = exp_factor*sin(phase);
      }
    }
    //End tw local
  }else{
    for(index_tw_local=0;index_tw_local<tw_local_size;++index_tw_local){
      logxi = log((pma->tau0-tw_local_sampling[index_wd*tw_local_size+index_tw_local]));//285+pma->small_log_offset));
      exp_factor = exp(logxi*(offset-pma->nu_real[index_tilt1_tilt2]));//was 1.0 instead of offset
      //exp_factor = exp(logxi*(1.0-pma->bias));
      for(index_coeff=0;index_coeff<pma->size_fft_cutoff;++index_coeff){
        phase = -logxi*pma->nu_imag[index_coeff];
        pref_real[index_tw_local*pma->size_fft_result+index_coeff] = exp_factor*cos(phase);
        pref_imag[index_tw_local*pma->size_fft_result+index_coeff] = exp_factor*sin(phase);
      }
    }
    //End tw local
  }
  return _SUCCESS_;
}
int matter_integrate_cosmological_function(
                                           struct precision* ppr,
                                           struct background* pba,
                                           struct perturbs* ppt,
                                           struct matters* pma,
                                           int index_ic1,
                                           int index_ic2,
                                           int index_radtp1,
                                           int index_radtp2,
                                           int index_stp1,
                                           int index_stp2,
                                           int index_tilt1_tilt2,
                                           int index_ic1_ic2,
                                           int index_stp1_stp2,
                                           int index_cltp,
                                           int index_wd1,
                                           int index_wd2,
                                           double** integrand_real,
                                           double** integrand_imag,
                                           double** fft_coeff_real,
                                           double** fft_coeff_imag,
                                           double*** win_fft_real,
                                           double*** win_fft_imag,
                                           double* tw_local_sampling,
                                           double* tw_local_weights,
                                           int tw_local_size,
                                           double* chi_pref_real,
                                           double* chi_pref_imag,
                                           int is_integrated1,
                                           int is_integrated2,
                                           int integrate_logarithmically,
                                           double t_min,
                                           double t_max,
                                           int tw_max_size
                                          ){
  int abort = _FALSE_;
  int index_spl;
  double intxi_local_real,intxi_local_imag;
  int t_size_local = (pma->uses_intxi_interpolation?pma->t_spline_size:pma->t_size);
  double t;
  double* int_real;
  double* int_imag;
  double** window_fft_real;
  double** window_fft_imag;
  double y_min,y_max;
  //int integrated_t_offset;
  y_min = -log(1-t_min);
  y_max = -log(1-t_max);
  //integrated_t_offset = pma->t_size*((is_integrated1==1 || is_integrated2 == 1)?1:0);//pma->t_size*(is_integrated1*2+is_integrated2);
  class_call(matter_precompute_chit_factors(pma,
                                            tw_local_sampling,
                                            tw_local_size,
                                            index_tilt1_tilt2,
                                            index_wd1,
                                            pma->uses_intxi_logarithmic  && is_integrated1 == 1 && is_integrated2 == 1 ,
                                            chi_pref_real,
                                            chi_pref_imag),
             pma->error_message,
             pma->error_message);
  if(pma->uses_intxi_symmetrized && is_integrated1 == 1 && is_integrated2 == 1  ){
    class_call(matter_precompute_chit_factors(pma,
                                              tw_local_sampling,
                                              tw_local_size,
                                              index_tilt1_tilt2,
                                              index_wd2,
                                              pma->uses_intxi_logarithmic && is_integrated1 == 1 && is_integrated2 == 1 ,
                                              chi_pref_real+tw_max_size*pma->size_fft_result,
                                              chi_pref_imag+tw_max_size*pma->size_fft_result),
             pma->error_message,
             pma->error_message);
  }
  int index_t,index_coeff,index_tw_local;
  #pragma omp parallel private(index_t,index_coeff,index_tw_local,t,int_real,int_imag,window_fft_real,window_fft_imag) firstprivate(pma,pba,ppt,is_integrated1,is_integrated2)
  {
    int tid = omp_get_thread_num();
    int_real = integrand_real[tid];
    int_imag = integrand_imag[tid];
    window_fft_real = win_fft_real[tid];
    window_fft_imag = win_fft_imag[tid];
  /**
   * Now obtain the f_n^{ij}(t) function
   * */
  if(pma->uses_intxi_logarithmic && pma->uses_intxi_symmetrized && is_integrated1 == 1 && is_integrated2 == 1 ){
    /**
     * Now obtain the f_n^{ij}(t) function
     * */
    int index_stp2_stp1 = index_stp2*pma->stp_size+index_stp1;
    int index_ic2_ic1 = index_ic2*pma->ic_size+index_ic1;
    #pragma omp for
    for(index_t=0;index_t<t_size_local;++index_t){
      matter_get_t(index_t)
      class_call_parallel(matter_get_half_integrand(pba,
                                ppt,
                                pma,
                                t,
                                index_ic1,
                                index_ic2,
                                index_radtp1,
                                index_radtp2,
                                index_tilt1_tilt2,
                                index_ic1_ic2,
                                index_stp1_stp2,
                                index_cltp,
                                index_wd1,
                                index_wd2,
                                int_real,
                                int_imag,
                                fft_coeff_real,
                                fft_coeff_imag,
                                window_fft_real,
                                window_fft_imag,
                                tw_local_sampling,
                                tw_local_size
                                ),
                  pma->error_message,
                  pma->error_message);
      class_call_parallel(matter_get_half_integrand(pba,
                                ppt,
                                pma,
                                t,
                                index_ic2,
                                index_ic1,
                                index_radtp2,
                                index_radtp1,
                                index_tilt1_tilt2,
                                index_ic2_ic1,
                                index_stp2_stp1,
                                index_cltp,
                                index_wd2,
                                index_wd1,
                                int_real+tw_max_size*pma->size_fft_result,
                                int_imag+tw_max_size*pma->size_fft_result,
                                fft_coeff_real,
                                fft_coeff_imag,
                                window_fft_real,
                                window_fft_imag,
                                tw_local_sampling,
                                tw_local_size
                                ),
                  pma->error_message,
                  pma->error_message);
      for(index_coeff=0;index_coeff<pma->size_fft_cutoff;++index_coeff){
        double sum_real =0.0;
        double sum_imag =0.0;
        for(index_tw_local=0;index_tw_local<tw_local_size;++index_tw_local){
          sum_real+=tw_local_weights[index_wd1*tw_local_size+index_tw_local]*(chi_pref_real[index_tw_local*pma->size_fft_result+index_coeff]*int_real[index_coeff*tw_local_size+index_tw_local]-chi_pref_imag[index_tw_local*pma->size_fft_result+index_coeff]*int_imag[index_coeff*tw_local_size+index_tw_local]);
          sum_imag+=tw_local_weights[index_wd1*tw_local_size+index_tw_local]*(chi_pref_real[index_tw_local*pma->size_fft_result+index_coeff]*int_imag[index_coeff*tw_local_size+index_tw_local]+chi_pref_imag[index_tw_local*pma->size_fft_result+index_coeff]*int_real[index_coeff*tw_local_size+index_tw_local]);
        }
        for(index_tw_local=0;index_tw_local<tw_local_size;++index_tw_local){
          sum_real+=tw_local_weights[index_wd2*tw_local_size+index_tw_local]*(chi_pref_real[tw_max_size*pma->size_fft_result+index_tw_local*pma->size_fft_result+index_coeff]*int_real[index_coeff*tw_local_size+index_tw_local+tw_max_size*pma->size_fft_result]-chi_pref_imag[tw_max_size*pma->size_fft_result+index_tw_local*pma->size_fft_result+index_coeff]*int_imag[index_coeff*tw_local_size+index_tw_local+tw_max_size*pma->size_fft_result]);
          sum_imag+=tw_local_weights[index_wd2*tw_local_size+index_tw_local]*(chi_pref_real[tw_max_size*pma->size_fft_result+index_tw_local*pma->size_fft_result+index_coeff]*int_imag[index_coeff*tw_local_size+index_tw_local+tw_max_size*pma->size_fft_result]+chi_pref_imag[tw_max_size*pma->size_fft_result+index_tw_local*pma->size_fft_result+index_coeff]*int_real[index_coeff*tw_local_size+index_tw_local+tw_max_size*pma->size_fft_result]);
        }
        //End tw integration
        if(pma->uses_intxi_interpolation){
          pma->intxi_spline_real[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_spline_size+index_t] = sum_real;
          pma->intxi_spline_imag[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_spline_size+index_t] = sum_imag;
        }
        else{
          pma->intxi_real[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_size+index_t] = sum_real;
          pma->intxi_imag[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_size+index_t] = sum_imag;
        }
        //Ifend isinterpolated
      }
      //End coeff
    }
    //End t
  }
  else if(pma->uses_intxi_logarithmic && is_integrated1==1 && is_integrated2 == 1){
    #pragma omp for
    for(index_t=0;index_t<t_size_local;++index_t){
      matter_get_t(index_t)
      class_test_parallel(t==0,
           pma->error_message,
           "stop to avoid division by zero or logarithm of zero");
      class_call_parallel(matter_get_ttau_integrand(pba,
                                ppt,
                                pma,
                                t,
                                index_ic1,
                                index_ic2,
                                index_radtp1,
                                index_radtp2,
                                index_tilt1_tilt2,
                                index_ic1_ic2,
                                index_stp1_stp2,
                                index_cltp,
                                index_wd1,
                                index_wd2,
                                int_real,
                                int_imag,
                                fft_coeff_real,
                                fft_coeff_imag,
                                window_fft_real,
                                window_fft_imag,
                                tw_local_sampling,
                                tw_local_size
                                ),
                  pma->error_message,
                  pma->error_message);
      for(index_coeff=0;index_coeff<pma->size_fft_cutoff;++index_coeff){
        double sum_real =0.0;
        double sum_imag =0.0;
        for(index_tw_local=0;index_tw_local<tw_local_size;++index_tw_local){
          sum_real+=tw_local_weights[index_wd1*tw_local_size+index_tw_local]*(chi_pref_real[index_tw_local*pma->size_fft_result+index_coeff]*int_real[index_coeff*tw_local_size+index_tw_local]-chi_pref_imag[index_tw_local*pma->size_fft_result+index_coeff]*int_imag[index_coeff*tw_local_size+index_tw_local]);
          sum_imag+=tw_local_weights[index_wd1*tw_local_size+index_tw_local]*(chi_pref_real[index_tw_local*pma->size_fft_result+index_coeff]*int_imag[index_coeff*tw_local_size+index_tw_local]+chi_pref_imag[index_tw_local*pma->size_fft_result+index_coeff]*int_real[index_coeff*tw_local_size+index_tw_local]);
        }
        //End tw integration
        if(pma->uses_intxi_interpolation){
          pma->intxi_spline_real[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_spline_size+index_t] = sum_real;
          pma->intxi_spline_imag[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_spline_size+index_t] = sum_imag;
        }
        else{
          pma->intxi_real[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_size+index_t] = sum_real;
          pma->intxi_imag[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_size+index_t] = sum_imag;
        }
        //Ifend isinterpolated
      }
      //End coeff
    }
    //End t
  }
  else if(pma->uses_intxi_symmetrized && (is_integrated1==1 && is_integrated2==1)){
    /**
     * Now obtain the f_n^{ij}(t) function
     * */
    int index_stp2_stp1 = index_stp2*pma->stp_size+index_stp1;
    int index_ic2_ic1 = index_ic2*pma->stp_size+index_ic1;
    #pragma omp for
    for(index_t=0;index_t<t_size_local;++index_t){
      matter_get_t(index_t)
      class_test_parallel(t==0,
           pma->error_message,
           "stop to avoid division by zero or logarithm of zero");
      class_call_parallel(matter_get_half_integrand(pba,
                                ppt,
                                pma,
                                t,
                                index_ic1,
                                index_ic2,
                                index_radtp1,
                                index_radtp2,
                                index_tilt1_tilt2,
                                index_ic1_ic2,
                                index_stp1_stp2,
                                index_cltp,
                                index_wd1,
                                index_wd2,
                                int_real,
                                int_imag,
                                fft_coeff_real,
                                fft_coeff_imag,
                                window_fft_real,
                                window_fft_imag,
                                tw_local_sampling,
                                tw_local_size
                                ),
                  pma->error_message,
                  pma->error_message);
      class_call_parallel(matter_get_half_integrand(pba,
                                ppt,
                                pma,
                                t,
                                index_ic2,
                                index_ic1,
                                index_radtp2,
                                index_radtp1,
                                index_tilt1_tilt2,
                                index_ic2_ic1,
                                index_stp2_stp1,
                                index_cltp,
                                index_wd2,
                                index_wd1,
                                int_real+tw_max_size*pma->size_fft_result,
                                int_imag+tw_max_size*pma->size_fft_result,
                                fft_coeff_real,
                                fft_coeff_imag,
                                window_fft_real,
                                window_fft_imag,
                                tw_local_sampling,
                                tw_local_size
                                ),
                  pma->error_message,
                  pma->error_message);
      for(index_coeff=0;index_coeff<pma->size_fft_cutoff;++index_coeff){
        double sum_real =0.0;
        double sum_imag =0.0;
        for(index_tw_local=0;index_tw_local<tw_local_size;++index_tw_local){
          sum_real+=tw_local_weights[index_wd1*tw_local_size+index_tw_local]*(chi_pref_real[index_tw_local*pma->size_fft_result+index_coeff]*int_real[index_coeff*tw_local_size+index_tw_local]-chi_pref_imag[index_tw_local*pma->size_fft_result+index_coeff]*int_imag[index_coeff*tw_local_size+index_tw_local]);
          sum_imag+=tw_local_weights[index_wd1*tw_local_size+index_tw_local]*(chi_pref_real[index_tw_local*pma->size_fft_result+index_coeff]*int_imag[index_coeff*tw_local_size+index_tw_local]+chi_pref_imag[index_tw_local*pma->size_fft_result+index_coeff]*int_real[index_coeff*tw_local_size+index_tw_local]);
        }
        for(index_tw_local=0;index_tw_local<tw_local_size;++index_tw_local){
          sum_real+=tw_local_weights[index_wd2*tw_local_size+index_tw_local]*(chi_pref_real[tw_max_size*pma->size_fft_result+index_tw_local*pma->size_fft_result+index_coeff]*int_real[index_coeff*tw_local_size+index_tw_local+tw_max_size*pma->size_fft_result]-chi_pref_imag[tw_max_size*pma->size_fft_result+index_tw_local*pma->size_fft_result+index_coeff]*int_imag[index_coeff*tw_local_size+index_tw_local+tw_max_size*pma->size_fft_result]);
          sum_imag+=tw_local_weights[index_wd2*tw_local_size+index_tw_local]*(chi_pref_real[tw_max_size*pma->size_fft_result+index_tw_local*pma->size_fft_result+index_coeff]*int_imag[index_coeff*tw_local_size+index_tw_local+tw_max_size*pma->size_fft_result]+chi_pref_imag[tw_max_size*pma->size_fft_result+index_tw_local*pma->size_fft_result+index_coeff]*int_real[index_coeff*tw_local_size+index_tw_local+tw_max_size*pma->size_fft_result]);
        }
        //End tw integration
        if(pma->uses_intxi_interpolation){
          pma->intxi_spline_real[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_spline_size+index_t] = sum_real;
          pma->intxi_spline_imag[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_spline_size+index_t] = sum_imag;
        }
        else{
          pma->intxi_real[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_size+index_t] = sum_real;
          pma->intxi_imag[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_size+index_t] = sum_imag;

        }
        //Ifend isinterpolated
      }
      //End coeff
    }
    //End t
  }
  else{
    #pragma omp for
    for(index_t=0;index_t<t_size_local;++index_t){
      matter_get_t(index_t)
      class_test_parallel(t==0,
           pma->error_message,
           "stop to avoid division by zero or logarithm of zero");
      class_call_parallel(matter_get_ttau_integrand(pba,
                                ppt,
                                pma,
                                t,
                                index_ic1,
                                index_ic2,
                                index_radtp1,
                                index_radtp2,
                                index_tilt1_tilt2,
                                index_ic1_ic2,
                                index_stp1_stp2,
                                index_cltp,
                                index_wd1,
                                index_wd2,
                                int_real,
                                int_imag,
                                fft_coeff_real,
                                fft_coeff_imag,
                                window_fft_real,
                                window_fft_imag,
                                tw_local_sampling,
                                tw_local_size
                                ),
                  pma->error_message,
                  pma->error_message);
      for(index_coeff=0;index_coeff<pma->size_fft_cutoff;++index_coeff){
        double sum_real =0.0;
        double sum_imag =0.0;
        for(index_tw_local=0;index_tw_local<tw_local_size;++index_tw_local){
          sum_real+=tw_local_weights[index_wd1*tw_local_size+index_tw_local]*(chi_pref_real[index_tw_local*pma->size_fft_result+index_coeff]*int_real[index_coeff*tw_local_size+index_tw_local]-chi_pref_imag[index_tw_local*pma->size_fft_result+index_coeff]*int_imag[index_coeff*tw_local_size+index_tw_local]);
          sum_imag+=tw_local_weights[index_wd1*tw_local_size+index_tw_local]*(chi_pref_real[index_tw_local*pma->size_fft_result+index_coeff]*int_imag[index_coeff*tw_local_size+index_tw_local]+chi_pref_imag[index_tw_local*pma->size_fft_result+index_coeff]*int_real[index_coeff*tw_local_size+index_tw_local]);
        }
        //End tw integration
        if(pma->uses_intxi_interpolation){
          pma->intxi_spline_real[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_spline_size+index_t] = sum_real;
          pma->intxi_spline_imag[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_spline_size+index_t] = sum_imag;
        }
        else{
          pma->intxi_real[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_size+index_t] = sum_real;
          pma->intxi_imag[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_size+index_t] = sum_imag;
        }
        //Ifend isinterpolated
      }
      //End coeff
    }
    //End t
  }
  //Ifend
  }
  if(abort == _TRUE_){return _FAILURE_;}
  //End parallel
  if(pma->uses_intxi_interpolation){
    if(pma->uses_intxi_asymptotic && !(matter_is_integrated(index_radtp1) || matter_is_integrated(index_radtp2))){
      for(index_t=0;index_t<t_size_local;++index_t){
        matter_get_t(index_t)
        class_test(t==0,
           pma->error_message,
           "stop to avoid division by zero or logarithm of zero");
        double temp;
        class_call(matter_asymptote(ppr, pma, t, index_wd1,index_wd2, &temp),pma->error_message,pma->error_message);
        for(index_coeff=0;index_coeff<pma->size_fft_cutoff;++index_coeff){
          pma->intxi_spline_real[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_spline_size+index_t] /= temp;
        }
        //End coeff
      }
      //End t
    }
    //Ifend asymptotic
    /**
     * If we want spline interpolation,
     *  we first have to calculate the splines,
     *  and then we interpolate said splines
     * */
    for(index_coeff=0;index_coeff<pma->size_fft_cutoff;++index_coeff){
      class_call(array_spline_table_columns(pma->t_spline_sampling,
                                            pma->t_spline_size,
                                            pma->intxi_spline_real[index_radtp1*pma->radtp_size_total+index_radtp2]+index_coeff*pma->t_spline_size,
                                            1,
                                            pma->ddintxi_spline_real[index_radtp1*pma->radtp_size_total+index_radtp2]+index_coeff*pma->t_spline_size,
                                            _SPLINE_EST_DERIV_,
                                            pma->error_message),
           pma->error_message,
           pma->error_message);
      class_call(array_spline_table_columns(pma->t_spline_sampling,
                                            pma->t_spline_size,
                                            pma->intxi_spline_imag[index_radtp1*pma->radtp_size_total+index_radtp2]+index_coeff*pma->t_spline_size,
                                            1,
                                            pma->ddintxi_spline_imag[index_radtp1*pma->radtp_size_total+index_radtp2]+index_coeff*pma->t_spline_size,
                                            _SPLINE_EST_DERIV_,
                                            pma->error_message),
            pma->error_message,
            pma->error_message);
      index_spl = 0;
      for(index_t=0;index_t<pma->t_size;++index_t){
        //TODO :: was without orig before
        matter_get_t_orig(index_t);
        double a,b,h;
        class_call(matter_spline_hunt(pma->t_spline_sampling,
                           pma->t_spline_size,
                           pma->t_sampling[index_t],
                           &index_spl,
                           &h,
                           &a,
                           &b,
                           pma->error_message),
                  pma->error_message,
                  pma->error_message);
        intxi_local_real = b*pma->intxi_spline_real[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_spline_size+index_spl+1]
                          +a*pma->intxi_spline_real[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_spline_size+index_spl]
                    +(
                      (b*b*b-b)*pma->ddintxi_spline_real[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_spline_size+index_spl+1]+
                      (a*a*a-a)*pma->ddintxi_spline_real[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_spline_size+index_spl]
                     )*h*h/6.;
        intxi_local_imag = b*pma->intxi_spline_imag[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_spline_size+index_spl+1]
                          +a*pma->intxi_spline_imag[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_spline_size+index_spl]
                    +(
                      (b*b*b-b)*pma->ddintxi_spline_imag[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_spline_size+index_spl+1]+
                      (a*a*a-a)*pma->ddintxi_spline_imag[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_spline_size+index_spl]
                     )*h*h/6.;
        pma->intxi_real[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_size+index_t] = intxi_local_real;
        pma->intxi_imag[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_size+index_t] = intxi_local_imag;
        if(pma->uses_intxi_asymptotic && !(matter_is_integrated(index_radtp1) || matter_is_integrated(index_radtp2))){
          double temp;
          //TODO :: can be optimized (no coeff dependence)
          class_call(matter_asymptote(ppr, pma, t, index_wd1,index_wd2, &temp),pma->error_message,pma->error_message);
          pma->intxi_real[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_size+index_t] *= temp;
        }
        //Ifend asymptotic
      }
      //End t
    }
    //End coeff
  }
  //Ifend interpolation
  return _SUCCESS_;
}
int matter_integrate_for_each_ttau_parallel_chi_pre(
                        struct precision* ppr,
                        struct background* pba,
                        struct perturbs * ppt,
                        struct matters* pma,
                        double ** fft_coeff_real,
                        double ** fft_coeff_imag,
                        int index_ic1,
                        int index_ic2,
                        int index_ic1_ic2,
                        int index_wd1,
                        int index_wd2,
                        int index_cltp,
                        double*** window_fft_real,
                        double*** window_fft_imag,
                        double*** window_bessel_real,
                        double*** window_bessel_imag,
                        short integrate_logarithmically
                        ){
  /**
   * Define and allocate local variables
   * */
  int index_wd1_wd2;
  int index_coeff;
  int index_radtp1,index_radtp2;
  int index_radtp_of_bitp1,index_radtp_of_bitp2;
  int index_bitp1,index_bitp2;
  int index_tilt1,index_tilt2,index_tilt1_tilt2;
  int thread_id;
  double sum_temp = 0.0;
  int index_l;
  double t_min,t_max,y_min,y_max;
  int index_t;
  double t;
  double intxi_local_real,intxi_local_imag;
  int integrated_t_offset;
  index_wd1_wd2 = index_symmetric_matrix(index_wd1,index_wd2,pma->num_windows);

  int print_total_index = 0;
  int tw_max_size = 0;
  if(pma->has_unintegrated_windows){
    tw_max_size = MAX(tw_max_size,pma->tw_size);
  }
  if(pma->has_integrated_windows){
    tw_max_size = MAX(tw_max_size,pma->integrated_tw_size);
  }
  double* tw_local_sampling;
  double* tw_local_weights;
  int tw_local_size;
  int is_integrated1,is_integrated2;
  double* sum_l;
  class_alloc(sum_l,
              pma->l_size*sizeof(double),
              pma->error_message);
  memset(sum_l,0,pma->l_size*sizeof(double));
  double** integrand_real = NULL;
  double** integrand_imag = NULL;
  int N_threads = 1;
  //#ifdef _USE_OMP_
  N_threads = omp_get_max_threads();
  //#endif
  class_alloc(integrand_real,
              N_threads*sizeof(double*),
              pma->error_message);
  class_alloc(integrand_imag,
              N_threads*sizeof(double*),
              pma->error_message);
  for(thread_id=0;thread_id<N_threads;++thread_id){
    if(pma->uses_intxi_symmetrized){
      class_alloc(integrand_real[thread_id],
                  2*tw_max_size*pma->size_fft_result*sizeof(double),
                  pma->error_message);
      class_alloc(integrand_imag[thread_id],
                  2*tw_max_size*pma->size_fft_result*sizeof(double),
                  pma->error_message);
    }
    else{
      class_alloc(integrand_real[thread_id],
                  tw_max_size*pma->size_fft_result*sizeof(double),
                  pma->error_message);
      class_alloc(integrand_imag[thread_id],
                  tw_max_size*pma->size_fft_result*sizeof(double),
                  pma->error_message);
    }
  }
  double* chi_pref_real;
  double* chi_pref_imag;
  if(pma->uses_intxi_symmetrized){
    class_alloc(chi_pref_real,
                2*tw_max_size*pma->size_fft_result*sizeof(double),
                pma->error_message);
    class_alloc(chi_pref_imag,
                2*tw_max_size*pma->size_fft_result*sizeof(double),
                pma->error_message);
  }
  else{
    class_alloc(chi_pref_real,
                tw_max_size*pma->size_fft_result*sizeof(double),
                pma->error_message);
    class_alloc(chi_pref_imag,
                tw_max_size*pma->size_fft_result*sizeof(double),
                pma->error_message);
  }
  /**
   * Get the necessary bessel integrals
   *   and store in pre-prepared arrays
   * */
  if(pma->has_unintegrated_windows){
    int abort = _FALSE_;
    #pragma omp parallel private(index_l) firstprivate(pma,pba,window_bessel_real,window_bessel_imag)
    {
      #pragma omp for
      for(index_l=0;index_l<pma->l_size;++index_l){
        if(!pma->uses_limber_approximation){
          int index_l_eval = (pma->uses_bessel_recursion?pma->l_sampling[index_l]:index_l);
          class_call_parallel(matter_get_bessel_fort_parallel(
                          pba,
                          pma,
                          index_wd1,
                          index_wd2,
                          index_l_eval,
                          window_bessel_real[index_l],
                          window_bessel_imag[index_l],
                          integrate_logarithmically
                          ),
                  pma->error_message,
                  pma->error_message);
        }
        //Ifend limber
      }
      //End l
    }
    if(abort == _TRUE_) {return _FAILURE_;}
  }

  /**
   * Now iterate through all bessel integral types
   *  and radial types
   * */
  short switched_flag = _FALSE_;
  short type_doubling = _FALSE_;
  for(index_bitp1 = 0; index_bitp1< pma->bitp_size;++index_bitp1){
    if(pma->has_bitp_normal && index_bitp1 == pma->bitp_index_normal){
      index_tilt1 = pma->tilt_index_normal;
    }
    else{
      index_tilt1 = pma->tilt_index_reduced;
    }
    for(index_bitp2 = 0; index_bitp2< pma->bitp_size;++index_bitp2){
      if(pma->has_bitp_normal && index_bitp2 == pma->bitp_index_normal){
        index_tilt2 = pma->tilt_index_normal;
      }
      else{
        index_tilt2 = pma->tilt_index_reduced;
      }
      index_tilt1_tilt2 = index_symmetric_matrix(index_tilt1,index_tilt2,pma->tilt_size);
      for(index_radtp_of_bitp1 =0; index_radtp_of_bitp1 < pma->radtp_of_bitp_size[index_cltp*pma->bitp_size+index_bitp1];++index_radtp_of_bitp1){
        for(index_radtp_of_bitp2 =0; index_radtp_of_bitp2 < pma->radtp_of_bitp_size[index_cltp*pma->bitp_size+index_bitp2]; ++index_radtp_of_bitp2){
          index_radtp1 = pma->radtps_of_bitp[index_cltp*pma->bitp_size+index_bitp1][index_radtp_of_bitp1];
          index_radtp2 = pma->radtps_of_bitp[index_cltp*pma->bitp_size+index_bitp2][index_radtp_of_bitp2];
          if((index_wd1==index_wd2) && (index_ic1 == index_ic2)){
            if(index_radtp2 > index_radtp1){type_doubling = _FALSE_;continue;}
            else if(index_radtp2 < index_radtp1){type_doubling = _TRUE_;}
            else{type_doubling = _FALSE_;}
          }else{type_doubling = _FALSE_;}
          if(matter_is_integrated(index_radtp1) && !(matter_is_integrated(index_radtp2))){
            int temp = index_wd2;
            index_wd2 = index_wd1;
            index_wd1 = temp;
            temp = index_radtp1;
            index_radtp1 = index_radtp2;
            index_radtp2 = temp;
            temp = index_ic1;
            index_ic2 = index_ic1;
            index_ic1 = temp;
            switched_flag = _TRUE_;
          }
          int index_stp1 = pma->index_stp_of_radtp[index_radtp1];
          int index_stp2 = pma->index_stp_of_radtp[index_radtp2];
          int index_stp1_stp2 = index_stp1*pma->stp_size+index_stp2;
          if(pma->uses_density_splitting && pma->has_stp_delta_m && pma->uses_limber_approximation && (index_radtp1 == pma->radtp_dens1 || index_radtp2 == pma->radtp_dens1)){
            continue;
          }

          print_total_index++;
          if(pma->matter_verbose > MATTER_VERBOSITY_CLCALCULATION && !MATTER_REWRITE_PRINTING){
            printf(" -> BI types [%1d,%1d] (sizes [%2d,%2d]), RAD types [%2d,%2d] (Total %3d/%3d)",index_bitp1,index_bitp2,pma->radtp_of_bitp_size[index_bitp1],pma->radtp_of_bitp_size[index_bitp2],index_radtp1,index_radtp2,print_total_index,pma->radtp_size_total*pma->radtp_size_total);
          }
          if(pma->matter_verbose > MATTER_VERBOSITY_CLCALCULATION && MATTER_REWRITE_PRINTING){
            printf("\r -> BI types [%1d,%1d] (sizes [%2d,%2d]), RAD types [%2d,%2d] (Total %3d/%3d)",index_bitp1,index_bitp2,pma->radtp_of_bitp_size[index_bitp1],pma->radtp_of_bitp_size[index_bitp2],index_radtp1,index_radtp2,print_total_index,pma->radtp_size_total*pma->radtp_size_total);
            fflush(stdout);
          }
          /**
           * First define correct t sampling
           * */
          matter_get_t_limits(index_wd1,index_wd2)
          tw_local_weights = pma->tw_weights;
          tw_local_size = pma->tw_size;
          tw_local_sampling = pma->tw_sampling;
          if(
            matter_is_integrated(index_radtp1)
          ){
            t_min = 0.0+pma->bi_maximal_t_offset;//186+pma->bi_maximal_t_offset;
            t_max = 1.0-pma->bi_maximal_t_offset;//186-pma->bi_maximal_t_offset found important;
            tw_local_weights = pma->integrated_tw_weights;
            tw_local_size = pma->integrated_tw_size;
            tw_local_sampling = pma->integrated_tw_sampling;
            is_integrated1=1;
          }
          else{
            is_integrated1=0;
          }
          if(
            matter_is_integrated(index_radtp2)
          ){
            t_min = 0.0+pma->bi_maximal_t_offset;//186+pma->bi_maximal_t_offset;
            t_max = 1.0-pma->bi_maximal_t_offset;//186-pma->bi_maximal_t_offset found important;
            is_integrated2=1;
          }
          else{
            is_integrated2=0;
          }
          if(pma->matter_verbose > MATTER_VERBOSITY_CLCALCULATION && pma->matter_verbose > MATTER_VERBOSITY_RANGES && !MATTER_REWRITE_PRINTING){
            printf(" -> t range from %.10e to %.10e \n",t_min,t_max);
          }
          //TODO :: reformulate
          class_test(t_min >= t_max,
                     pma->error_message,
                     "Adjust matter_t_offset \n");
          y_min = -log(1-t_min);
          y_max = -log(1-t_max);
          integrated_t_offset = pma->t_size*((is_integrated1==1 || is_integrated2 == 1)?1:0);//pma->t_size*(is_integrated1*2+is_integrated2);

          class_call(matter_integrate_cosmological_function(ppr,
                                                            pba,
                                                            ppt,
                                                            pma,
                                                            index_ic1,
                                                            index_ic2,
                                                            index_radtp1,
                                                            index_radtp2,
                                                            index_stp1,
                                                            index_stp2,
                                                            index_tilt1_tilt2,
                                                            index_ic1_ic2,
                                                            index_stp1_stp2,
                                                            index_cltp,
                                                            index_wd1,
                                                            index_wd2,
                                                            integrand_real,
                                                            integrand_imag,
                                                            fft_coeff_real,
                                                            fft_coeff_imag,
                                                            window_fft_real,
                                                            window_fft_imag,
                                                            tw_local_sampling,
                                                            tw_local_weights,
                                                            tw_local_size,
                                                            chi_pref_real,
                                                            chi_pref_imag,
                                                            is_integrated1,
                                                            is_integrated2,
                                                            integrate_logarithmically,
                                                            t_min,
                                                            t_max,
                                                            tw_max_size),
                       pma->error_message,
                       pma->error_message);
          /**
           * We have now obtained the f_n^{ij}(t)
           *
           * Now we can advance to integrating the final Cl's
           * */
          for(index_l=0;index_l<pma->l_size;++index_l){
            double sum_t = 0.0;
            for(index_t=0;index_t<pma->t_size;++index_t){
              //TODO :: was without orig before
              matter_get_t_orig(index_t)
              sum_temp = 0.0;
              double bes_local_real,bes_local_imag;

              /**
               * The reverse FFT has
               *
               *  index 0 with a factor of 1
               *  index N with a factor of 1
               *  index i with a factor of 2,
               *   since here there would theoretically be
               *   both index i and N-i,
               *   but we relate N-i to i by symmetry
               *   of having a real integrand
               *
               * */
              intxi_local_real = pma->intxi_real[index_radtp1*pma->radtp_size_total+index_radtp2][0*pma->t_size+index_t];
              intxi_local_imag = pma->intxi_imag[index_radtp1*pma->radtp_size_total+index_radtp2][0*pma->t_size+index_t];
              bes_local_real = window_bessel_real[index_l][index_tilt1_tilt2*pma->size_fft_result+0][index_t+integrated_t_offset];
              bes_local_imag = window_bessel_imag[index_l][index_tilt1_tilt2*pma->size_fft_result+0][index_t+integrated_t_offset];
              sum_temp +=intxi_local_real*bes_local_real-intxi_local_imag*bes_local_imag;

              if(pma->size_fft_cutoff==pma->size_fft_result){
                intxi_local_real = pma->intxi_real[index_radtp1*pma->radtp_size_total+index_radtp2][(pma->size_fft_result-1)*pma->t_size+index_t];
                intxi_local_imag = pma->intxi_imag[index_radtp1*pma->radtp_size_total+index_radtp2][(pma->size_fft_result-1)*pma->t_size+index_t];
                bes_local_real = window_bessel_real[index_l][index_tilt1_tilt2*pma->size_fft_result+(pma->size_fft_result-1)][index_t+integrated_t_offset];
                bes_local_imag = window_bessel_imag[index_l][index_tilt1_tilt2*pma->size_fft_result+(pma->size_fft_result-1)][index_t+integrated_t_offset];
                sum_temp +=intxi_local_real*bes_local_real-intxi_local_imag*bes_local_imag;
              }
              for(index_coeff=1;index_coeff<pma->size_fft_cutoff-1;++index_coeff){
                intxi_local_real = pma->intxi_real[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_size+index_t];
                intxi_local_imag = pma->intxi_imag[index_radtp1*pma->radtp_size_total+index_radtp2][index_coeff*pma->t_size+index_t];
                bes_local_real = window_bessel_real[index_l][index_tilt1_tilt2*pma->size_fft_result+index_coeff][index_t+integrated_t_offset];
                bes_local_imag = window_bessel_imag[index_l][index_tilt1_tilt2*pma->size_fft_result+index_coeff][index_t+integrated_t_offset];
                sum_temp +=2.0*(intxi_local_real*bes_local_real-intxi_local_imag*bes_local_imag);
              }
              if(integrate_logarithmically && !pma->uses_limber_approximation){
                sum_t+=(1-t)*pma->t_weights[index_t]*sum_temp;
              }
              else{
                sum_t+=pma->t_weights[index_t]*sum_temp;
              }
            }
            /**
             * Since we only integrate over t_min to t_max,
             *  we have to rescale the final value obtained
             * (The t_weights are always scaled as 0 to 1)
             * */
            if(!integrate_logarithmically && !pma->uses_limber_approximation){
              sum_t*=(t_max-t_min);
            }
            else if(!pma->uses_limber_approximation){
              sum_t*=(y_max-y_min);
            }
            if(pma->has_bitp_lfactor && index_bitp2 == pma->bitp_index_lfactor){
              sum_t*=pma->l_sampling[index_l]*(pma->l_sampling[index_l]+1.0);
            }
            if(pma->has_bitp_lfactor && index_bitp1 == pma->bitp_index_lfactor){
              sum_t*=pma->l_sampling[index_l]*(pma->l_sampling[index_l]+1.0);
            }
            if(pma->matter_verbose > MATTER_VERBOSITY_CLCALCULATION_PARTIAL && !MATTER_REWRITE_PRINTING){
              printf("(l:%i (bitp1 =%i ,bitp2 =%i, radtp1 =%i , radtp2 = %i) = %.10e (total = %.10e) \n",(int)pma->l_sampling[index_l],index_bitp1,index_bitp2,index_radtp1,index_radtp2,sum_t,sum_l[index_l]);
            }
            if(type_doubling == _TRUE_){sum_t*=2.;}
            sum_l[index_l] += sum_t;
          }
          //End l
          if(switched_flag==_TRUE_){
            int temp = index_wd2;
            index_wd2 = index_wd1;
            index_wd1 = temp;
            temp = index_radtp1;
            index_radtp1 = index_radtp2;
            index_radtp2 = temp;
            switched_flag = _FALSE_;
          }
        }
        //End radtp2
      }
      //End radtp1
    }
    //End bitp2
  }
  //End bitp1
  /**
   * Print final results
   * */
  if(MATTER_REWRITE_PRINTING && pma->matter_verbose > MATTER_VERBOSITY_CLCALCULATION){
    printf("\r -> Output for current Window Combination :                                 \n");
  }
  for(index_l=0;index_l<pma->l_size;++index_l){
    if(pma->matter_verbose > MATTER_VERBOSITY_CLCALCULATION){
      printf("(l:%i) = %.10e \n",(int)pma->l_sampling[index_l],sum_l[index_l]);
    }
    pma->cl[index_ic1_ic2*pma->cltp_size+index_cltp][index_wd1_wd2*pma->l_size+index_l] = sum_l[index_l];
  }
  for(thread_id=0;thread_id<N_threads;++thread_id){
    free(integrand_real[thread_id]);
    free(integrand_imag[thread_id]);
  }
  free(integrand_real);
  free(integrand_imag);
  free(chi_pref_real);
  free(chi_pref_imag);
  free(sum_l);
  return _SUCCESS_;
}
int matter_get_bessel_limber(
                            struct matters* pma,
                            int index_l,
                            double** window_bessel_real,
                            double** window_bessel_imag
                            ){
  int index_tilt1_tilt2;
  int index_coeff;
  double l0 = sqrt(pma->l_sampling[index_l]*(pma->l_sampling[index_l]+1.));
  for(index_tilt1_tilt2=0;index_tilt1_tilt2<pma->tilt_grid_size;++index_tilt1_tilt2){
    for(index_coeff=0;index_coeff<pma->size_fft_cutoff;++index_coeff){
      double exp_factor = exp(log(l0)*(pma->nu_real[index_tilt1_tilt2]-3.));
      double phase = log(l0)*pma->nu_imag[index_coeff];
      /**
       * Since at t=1, we cannot split integration into 1/t and t,
       *  only half of the single point t=1 contributes
       * The theoretical factor would be (2pi^2),
       *  but we get aformentioned additional factor of 1/2.
       * */
      window_bessel_real[index_tilt1_tilt2*pma->size_fft_result+index_coeff][0] = _PI_*_PI_*exp_factor*cos(phase);
      window_bessel_imag[index_tilt1_tilt2*pma->size_fft_result+index_coeff][0] = _PI_*_PI_*exp_factor*sin(phase);
      if(pma->has_integrated_windows){
        window_bessel_real[index_tilt1_tilt2*pma->size_fft_result+index_coeff][0+pma->t_size] = _PI_*_PI_*exp_factor*cos(phase);
        window_bessel_imag[index_tilt1_tilt2*pma->size_fft_result+index_coeff][0+pma->t_size] = _PI_*_PI_*exp_factor*sin(phase);
      }
    }
    //End coeff
  }
  //End tilt grid
  return _SUCCESS_;
}
int matter_get_derivative_type(
                               struct matters* pma,
                               int* derivative_type1,
                               int* derivative_type2,
                               int index_radtp1,
                               int index_radtp2
                              ){
  if(!pma->uses_relative_factors){
    if(pma->has_stp_delta_m && pma->uses_density_splitting){
      if(index_radtp1 == pma->radtp_dens1){
        *derivative_type1 = 4;
      }
      if(index_radtp2 == pma->radtp_dens1){
        *derivative_type2 = 4;
      }
    }
    if(pma->has_redshift_space_distortion && (!pma->uses_rsd_combination)){
      if(index_radtp1 == pma->radtp_dop1){
        *derivative_type1 = 1;
      }
      if(index_radtp1 == pma->radtp_rsd){
        *derivative_type1 = 2;
      }
      if(index_radtp2 == pma->radtp_dop1){
        *derivative_type2 = 1;
      }
      else if(index_radtp2 == pma->radtp_rsd){
        *derivative_type2 = 2;
      }
    }
  }
  if(pma->has_lensing_terms){
    if(index_radtp1 == pma->radtp_nclens){
      *derivative_type1 = -1;
    }
    if(index_radtp2 == pma->radtp_nclens){
      *derivative_type2 = -1;
    }
  }
  if(pma->has_cl_shear){
    if(index_radtp1 == pma->radtp_shlens){
      *derivative_type1 = -1;
    }
    if(index_radtp2 == pma->radtp_shlens){
      *derivative_type2 = -1;
    }
  }
  if(pma->has_gravitational_terms){
    if(index_radtp1 == pma->radtp_g4 || index_radtp1 == pma->radtp_g5){
      *derivative_type1 = -1;
    }
    if(index_radtp2 == pma->radtp_g4 || index_radtp2 == pma->radtp_g5){
      *derivative_type2 = -1;
    }
  }
  return _SUCCESS_;
}
int matter_get_bessel_fort_parallel(
                      struct background* pba,
                      struct matters* pma,
                      int index_wd1,
                      int index_wd2,
                      int index_l,
                      double** window_bessel_real,
                      double** window_bessel_imag,
                      short integrate_logarithmically
                      ){
  int index_coeff;
  double res_real,res_imag;
  int last_t = 0;
  double t_min,t_max;
  int index_tilt1,index_tilt2,index_tilt1_tilt2;
  matter_get_t_limits(index_wd1,index_wd2)
  double y_min,y_max;
  y_min = -log(1-t_min);
  y_max = -log(1-t_max);
  int index_t;
  double t,eval_pt=0.0;
  double a,b,h;
  double* bi_real_i;
  double* bi_imag_i;
  double* ddbi_real_i;
  double* ddbi_imag_i;
  //TODO :: check viability of loop reordering
  if(pma->has_unintegrated_windows){
    for(index_tilt1=0;index_tilt1<pma->tilt_size;++index_tilt1){
      for(index_tilt2=index_tilt1;index_tilt2<pma->tilt_size;++index_tilt2){
        index_tilt1_tilt2 = index_symmetric_matrix(index_tilt1,index_tilt2,pma->tilt_size);

        for(index_coeff=0;index_coeff<pma->size_fft_cutoff;++index_coeff){
          last_t = pma->bi_size[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff]-2;
          bi_real_i = pma->bi_real[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff];
          bi_imag_i = pma->bi_imag[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff];
          ddbi_real_i = pma->ddbi_real[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff];
          ddbi_imag_i = pma->ddbi_imag[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff];
          // HERE :: last_t t change
          last_t = 0;
          for(index_t=0;index_t<pma->t_size;++index_t){
            matter_get_t_orig(index_t);
            eval_pt = 1.0-t;
            if(eval_pt>pma->bi_max[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff]){
              window_bessel_real[index_tilt1_tilt2*pma->size_fft_result+index_coeff][index_t] = 0.0;
              window_bessel_imag[index_tilt1_tilt2*pma->size_fft_result+index_coeff][index_t] = 0.0;
              continue;
            }
            class_call(matter_spline_hunt(pma->bi_sampling[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                                          pma->bi_size[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                                          eval_pt,
                                          &last_t,
                                          &h,
                                          &a,
                                          &b,
                                          pma->error_message),
                         pma->error_message,
                         pma->error_message);
            res_real =a * bi_real_i[last_t] + b * bi_real_i[last_t+1] + ((a*a*a-a)* ddbi_real_i[last_t] + (b*b*b-b)* ddbi_real_i[last_t+1])*h*h/6.;
            res_imag =a * bi_imag_i[last_t] + b * bi_imag_i[last_t+1] + ((a*a*a-a)* ddbi_imag_i[last_t] + (b*b*b-b)* ddbi_imag_i[last_t+1])*h*h/6.;
            window_bessel_real[index_tilt1_tilt2*pma->size_fft_result+index_coeff][index_t] = res_real;
            window_bessel_imag[index_tilt1_tilt2*pma->size_fft_result+index_coeff][index_t] = res_imag;
          }
          //End t
        }
        //End coeff
      }
      //End tilt2
    }
    //End tilt1
  }
  //End iff
  return _SUCCESS_;
}
int matter_get_bessel_fort_parallel_integrated(
                      struct background* pba,
                      struct matters* pma,
                      int index_wd1,
                      int index_wd2,
                      int index_l,
                      double** window_bessel_real,
                      double** window_bessel_imag,
                      short integrate_logarithmically
                      ){
  int index_coeff;
  double res_real,res_imag;
  int last_t = 0;
  double t_min = 0.0,t_max = 1.0;
  int index_tilt1,index_tilt2,index_tilt1_tilt2;
  double a,b,h;
  double* bi_real_i;
  double* bi_imag_i;
  double* ddbi_real_i;
  double* ddbi_imag_i;
  if(pma->has_integrated_windows){
    class_test(!pma->has_tilt_reduced,
               pma->error_message,
               "Has integrated windows, but not reduced tilt. This is a bug.");
    for(index_tilt1=0;index_tilt1<pma->tilt_size;++index_tilt1){
      for(index_tilt2=index_tilt1;index_tilt2<pma->tilt_size;++index_tilt2){
        index_tilt1_tilt2 = index_symmetric_matrix(index_tilt1,index_tilt2,pma->tilt_size);
        t_max = 1.0-pma->bi_maximal_t_offset;//186-pma->bi_maximal_t_offset found important;
        t_min = 0.0+pma->bi_maximal_t_offset;//186+pma->bi_maximal_t_offset;
        double y_min,y_max;
        y_min = -log(1-t_min);
        y_max = -log(1-t_max);//+pma->bi_maximal_t_offset);
        int index_t;
        double t,eval_pt=0.0;
        for(index_coeff=0;index_coeff<pma->size_fft_cutoff;++index_coeff){
          last_t = pma->bi_size[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff]-2;
          bi_real_i = pma->bi_real[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff];
          bi_imag_i = pma->bi_imag[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff];
          ddbi_real_i = pma->ddbi_real[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff];
          ddbi_imag_i = pma->ddbi_imag[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff];
          // HERE :: last_t t change
          last_t = 0;
          for(index_t=0;index_t<pma->t_size;++index_t){
            matter_get_t_orig(index_t);
            eval_pt = 1.0-t;
            if(eval_pt>pma->bi_max[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff]){
              window_bessel_real[index_tilt1_tilt2*pma->size_fft_result+index_coeff][index_t+pma->t_size] = 0.0;
              window_bessel_imag[index_tilt1_tilt2*pma->size_fft_result+index_coeff][index_t+pma->t_size] = 0.0;
              continue;
            }
            class_call(matter_spline_hunt(pma->bi_sampling[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                                          pma->bi_size[index_tilt1_tilt2][index_l*pma->size_fft_result+index_coeff],
                                          eval_pt,
                                          &last_t,
                                          &h,
                                          &a,
                                          &b,
                                          pma->error_message),
                         pma->error_message,
                         pma->error_message);
            res_real =a * bi_real_i[last_t] + b * bi_real_i[last_t+1] + ((a*a*a-a)* ddbi_real_i[last_t] + (b*b*b-b)* ddbi_real_i[last_t+1])*h*h/6.;
            res_imag =a * bi_imag_i[last_t] + b * bi_imag_i[last_t+1] + ((a*a*a-a)* ddbi_imag_i[last_t] + (b*b*b-b)* ddbi_imag_i[last_t+1])*h*h/6.;
            window_bessel_real[index_tilt1_tilt2*pma->size_fft_result+index_coeff][index_t+pma->t_size] = res_real;
            window_bessel_imag[index_tilt1_tilt2*pma->size_fft_result+index_coeff][index_t+pma->t_size] = res_imag;
          }
          //End t
        }
        //End coeff
      }
      //End tilt2
    }
    //End tilt1
  }
  //End iff
  return _SUCCESS_;
}

int matter_spline_prepare_hunt(
  double* x_array,
  int x_size,
  double x,
  int* last,
  ErrorMsg errmsg){
  int inf,sup,mid;

  inf=0;
  sup=x_size-1;

  if (x_array[inf] < x_array[sup]){

    if (x < x_array[inf]) {
      *last = inf;
      return _SUCCESS_;
    }

    if (x > x_array[sup]) {
      *last = sup;
      return _SUCCESS_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x < x_array[mid]) {sup=mid;}
      else {inf=mid;}

    }

  }

  else {

    if (x < x_array[sup]) {
      *last = sup;
      return _SUCCESS_;
    }

    if (x > x_array[inf]) {
      *last = inf;
      return _SUCCESS_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x > x_array[mid]) {sup=mid;}
      else {inf=mid;}

    }

  }

  *last = inf;
  return _SUCCESS_;
}
int matter_spline_hunt(
  double* x_array,
  int x_size,
  double x,
  int* last,
  double* h,
  double* a,
  double* b,
  ErrorMsg errmsg
  ){
  int last_index = *last;
  //Old :: last_index = x_size-1 (and if >=x_size) , but here inf = last_index,sup=last_index+1 => x_size => eval of array AT x_size, which is disallowed
  if(last_index>=x_size-1){last_index=x_size-2;}
  if(last_index<0){last_index=0;}
  int inf,sup,mid,inc;
  inc=1;
  if (x >= x_array[last_index]) {
    if (x > x_array[x_size-1]) {
      sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,
        x,x_array[x_size-1]);
      return _FAILURE_;
    }
    /* try closest neighboor upward */
    inf = last_index;
    sup = inf + inc;
    if (x > x_array[sup]) {
      /* hunt upward */
      while (x > x_array[sup]) {
        inf = sup;
        inc += 1;
        sup += inc;
        if (sup > x_size-1) {
          sup = x_size-1;
        }
      }
      /* bisect */
      while (sup-inf > 1) {
        mid=(int)(0.5*(inf+sup));
        if (x < x_array[mid]) {sup=mid;}
        else {inf=mid;}
      }
    }
  }
  else {
    if (x < x_array[0]) {
      sprintf(errmsg,"%s(L:%d) : x=%.20e < x_min=%.20e",__func__,__LINE__,
        x,x_array[0]);
      return _FAILURE_;
    }
    /* try closest neighboor downward */
    sup = last_index;
    inf = sup - inc;
    if (x < x_array[inf]) {
      /* hunt downward */
      while (x < x_array[inf]) {
        sup = inf;
        inc += 1;
        inf -= inc;
        if (inf < 0) {
          inf = 0;
        }
      }
      /* bisect */
      while (sup-inf > 1) {
        mid=(int)(0.5*(inf+sup));
        if (x < x_array[mid]) {sup=mid;}
        else {inf=mid;}
      }
    }
  }
  last_index = inf;
  *last = last_index;
  *h = x_array[sup] - x_array[inf];
  *b = (x-x_array[inf])/(*h);
  *a = 1.0-(*b);
  return _SUCCESS_;
}
int matter_estimate_t_max_bessel(
                                 struct matters * pma,
                                 double l,
                                 double nu_imag,
                                 double* t_max_estimate
                                 ){
  double nu_imag_max = pma->nu_imag[pma->size_fft_result-1];
  if(l>40){
    //Exponential fit with two saturation curves to the t_max
    //Using additional l/(c+l) offsets with c constant
    //Done for nu=0, and nu=nu_max and interpolated between those two
    //(guess_off is basically a multiple of nu_max)
    double expval1 = (1.0-exp(-0.08405923801793776*l));
    double expval2 = (1.0-exp(-0.03269491513404876*l));
    double guess_loff_1 = 16.552260860083162;
    double guess_loff_2 = 86.60472131391394;
    double guess_off = 1.0/(2*fabs(nu_imag_max));//1.0/72.0;
    *t_max_estimate = (
                      expval1 * l / (guess_loff_1+l) +
                      guess_off*(
                                 -(expval1 * l / (guess_loff_1+l))+
                                 (expval2 * l / (guess_loff_2+l))
                               )*fabs(nu_imag)
                      );
  }
  else{
    //Parabolic fit with pts (0,0),(lowl,val_lowl),(highl,val_highl)
    //for nmin = nu_min (=0) and nmax = nu_max
    //Then interpolated linearly between both parabolas
    double l_offset = pma->l_sampling[0];
    double val_highl_nmin = 0.75;//0.835;
    double val_lowl_nmin = 0.28;//0.620;
    double val_highl_nmax = 0.44;//0.595;
    double val_lowl_nmax = 0.05;//0.270;
    double highl = 37-l_offset; //60;
    double lowl = 8-l_offset;//22;
    double slope_nmin = ((val_highl_nmin/highl-val_lowl_nmin/lowl)/(highl-lowl));
    double offset_nmin = val_lowl_nmin/lowl;
    double slope_nmax = ((val_highl_nmax/highl-val_lowl_nmax/lowl)/(highl-lowl));
    double offset_nmax = val_lowl_nmax/lowl;
    double val_at_nmin = (slope_nmin*(l-l_offset-lowl)+offset_nmin)*(l-l_offset);
    double val_at_nmax = (slope_nmax*(l-l_offset-lowl)+offset_nmax)*(l-l_offset);
    double interpol_val = val_at_nmin*(1.0-nu_imag/nu_imag_max) + val_at_nmax*nu_imag/nu_imag_max;
    *t_max_estimate = interpol_val;
  }
  return _SUCCESS_;
}
int matter_resample_growth_factor(
                                  struct matters * pma
                                  ){
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS){
    printf("Method :: Resample the growth factor\n");
  }
  int index_ic,index_stp,index_wd,index_tw;
  int last_index = 0;
  for (index_ic = 0; index_ic < pma->ic_size; index_ic++) {
    for (index_stp = 0; index_stp < pma->stp_size; index_stp++) {
      class_alloc(pma->growth_factor[index_ic * pma->stp_size + index_stp],
                  pma->num_windows*pma->tw_size*sizeof(double),
                  pma->error_message);
      for( index_wd = 0; index_wd < pma->num_windows ; ++index_wd){
        for( index_tw=0;index_tw<pma->tw_size;++index_tw){
          class_call(matter_interpolate_spline_growing_hunt(
                          pma->tau_sampling,
                          pma->tau_size,
                          pma->growth_factor_tau[index_ic*pma->stp_size + index_stp],
                          pma->ddgrowth_factor_tau[index_ic*pma->stp_size + index_stp],
                          1,
                          pma->tw_sampling[index_wd*pma->tw_size+index_tw],
                          &last_index,
                          &pma->growth_factor[index_ic * pma->stp_size + index_stp][index_wd*pma->tw_size+index_tw],
                          _FALSE_,
                          pma->error_message),
                     pma->error_message,
                     pma->error_message);
        }
        //End tw
      }
      //End wd
    }
    //End stp
  }
  //End ic
  return _SUCCESS_;
}
int matter_spline_growth_factor(
                  struct matters* pma
                  ){
  if(pma->matter_verbose > MATTER_VERBOSITY_FUNCTIONS){
    printf("Method :: Splining growth factor\n");
  }
  int index_ic,index_stp;
  for (index_ic = 0; index_ic < pma->ic_size; index_ic++) {
    for (index_stp = 0; index_stp < pma->stp_size; index_stp++) {
      class_alloc(pma->ddgrowth_factor_tau[index_ic * pma->stp_size + index_stp],
                  pma->tau_size*sizeof(double),
                  pma->error_message);
      array_spline_table_columns(pma->tau_sampling,
                                 pma->tau_size,
                                 pma->growth_factor_tau[index_ic * pma->stp_size + index_stp],
                                 1,
                                 pma->ddgrowth_factor_tau[index_ic * pma->stp_size + index_stp],
                                 _SPLINE_EST_DERIV_,
                                 pma->error_message);
    }
    //End stp
  }
  //End ic
  return _SUCCESS_;
}

int matter_derive(
                  double* x_array,
                  double* array,
                  int x_length,
                  double* dy_array,
                  ErrorMsg errmsg
                  ){
  int i;
  double dxp,dxm,dyp,dym;

  if (x_length < 2) {
    sprintf(errmsg,"%s(L:%d) routine called with n_lines=%d, should be at least 2",__func__,__LINE__,x_length);
    return _FAILURE_;
  }

  dxp = x_array[2] - x_array[1];
  dxm = x_array[0] - x_array[1];
  dyp = array[2] - array[1];
  dym = array[0] - array[1];

  if ((dxp*dxm*(dxm-dxp)) == 0.) {
    sprintf(errmsg,"%s(L:%d) stop to avoid division by zero",__func__,__LINE__);
    return _FAILURE_;
  }

  dy_array[1] = (dyp*dxm*dxm-dym*dxp*dxp)/(dxp*dxm*(dxm-dxp));

  dy_array[0] = dy_array[1] + dxm * 2.*(dyp*dxm-dym*dxp)/(dxp*dxm*(dxp-dxm));

  for (i=2; i<x_length-1; i++) {

    dxp = x_array[i+1] - x_array[i];
    dxm = x_array[i-1] - x_array[i];
    dyp = array[i+1] - array[i];
    dym = array[i-1] - array[i];

    if ((dxp*dxm*(dxm-dxp)) == 0.) {
      sprintf(errmsg,"%s(L:%d) stop to avoid division by zero",__func__,__LINE__);
      return _FAILURE_;
    }

    dy_array[i] = (dyp*dxm*dxm-dym*dxp*dxp)/(dxp*dxm*(dxm-dxp));

  }

  dy_array[x_length-1] = dy_array[x_length-2] + (x_array[x_length-1] - x_array[x_length-2]) * 2.*(dyp*dxm-dym*dxp)/(dxp*dxm*(dxp-dxm));

  return _SUCCESS_;
}
