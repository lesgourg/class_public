/** @file Temporary addition file for new transfer functionality
 *
 * Nils Schöneberg, 18.10.2017
 *
 */
 #ifndef __MATTER__
 #define __MATTER__
 #include "nonlinear.h"
 #include "fft.h"
 #include "extrapolate_source.h"
enum matter_integration_method {matter_integrate_tw_t,matter_integrate_tw_logt};
enum matter_k_weight_method {matter_k_weights_gaussian,matter_k_weights_step};
struct matters{
  short matter_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */
  ErrorMsg error_message; /**< zone for writing error messages */
  
  /**
   * Use/Has/Allow flags, general
   * */
  short has_cls;
  short uses_separability;
  
  short uses_intxi_interpolation;
  short uses_intxi_symmetrized;
  short uses_intxi_logarithmic;
  short uses_intxi_asymptotic;
  
  short uses_rsd_combination;
  short uses_density_splitting;
  short has_integrated_windows;
  short has_unintegrated_windows;
  short has_window_differentiation;
  short uses_bessel_recursion;
  short uses_limber_approximation;
  short uses_relative_factors;
  short uses_all_l_sampling;
  short uses_bessel_analytic_integration;
  short uses_lensing_reduction;
  short uses_filon_clenshaw_curtis;
  
  int* is_non_zero;
  int non_diag;
  enum matter_integration_method uses_integration;
  /**
   * Tau sampling and tau, other general quantities
   * */
  double* tau_sampling; // of size tau_size
  double tau0;
  double h;
  
  int tau_size_max; //The limit to growing tau_size
  int tau_size; //The number of samples in tau space
  int tau_grid_size; // = tau_size*(tau_size+1)/2;
  int index_tau_perturbs_beginning; //Index of start of tau_matter sampling in tau_perturbs
  
  
  /**
   * ICs
   * */
  int ic_size;
  int ic_ic_size;
  
  /**
   * Type indices and type 'has' conditions
   * */
  
  int stp_index_delta_m;
  short has_stp_delta_m;
  int stp_index_theta_m;
  short has_doppler_terms;
  short has_redshift_space_distortion;
  short has_gravitational_terms;
  int stp_index_phi;
  int stp_index_phi_prime;
  short has_lensing_terms;
  short has_cl_shear;
  short has_stp_phi_plus_psi;
  int stp_index_phi_plus_psi;
  int stp_index_psi;
  
  int stp_size;
  int stp_grid_size;
  
  int radtp_dens;
  int radtp_dens1;
  int radtp_dens2;
  int radtp_nclens;
  int radtp_dop1;
  int radtp_dop2;
  int radtp_rsd;
  int radtp_rsd_combined;
  int radtp_g1;
  int radtp_g2;
  int radtp_g3;
  int radtp_g4;
  int radtp_g5;
  int radtp_combined;
  
  int radtp_shlens;
  
  int* radtp_of_bitp_size;
  int** radtps_of_bitp;
  int radtp_size_total;
  int radtp_grid_size;
  
  int bitp_index_normal;
  short has_bitp_normal;
  int bitp_index_lfactor;
  short has_bitp_lfactor;
  int bitp_index_nu_reduced;
  short has_bitp_nu_reduced;
  
  short has_tilt_normal;
  int tilt_index_normal;
  short has_tilt_reduced;
  int tilt_index_reduced;
  
  int bitp_size;
  int tilt_size;
  int tilt_grid_size;
  
  int* index_perturb_tp_of_stp;
  int* index_stp_of_radtp;
  
  /**
   * Radial types (bessel or bessel derivatives) and
   * Cl types 
   * */
  int cltp_size;
  short has_cltp_nc;
  int cltp_index_nc;
  short has_cltp_sh;
  int cltp_index_sh;
  
  
  /**
   * Window max and mins
   * */
  double z_max;
  double tau_max;
  double tau_min; //Minimum tau given by normal window functions
  double* tw_max;
  double* tw_min;
  
  
  /**
   * FFT size
   * */
  int size_fft_input;
  int size_fft_result;
  int size_fft_cutoff;
  
  /**
   * Relative factors
   * */
  double* relative_factors;
  int index_relative_stp;
  
  /**
   * Growth factor
   * */
  double** growth_factor; //Sampled in tw
  double ** growth_factor_tau; // Sampled in tau
  double ** ddgrowth_factor_tau; //Sampled in tau
  double k_weight_k_max;
  double k_weight_k_min;
  enum matter_k_weight_method k_weight_mode; 
  
  /**
   * k sampling, normal and extrapolated
   * */
  short allow_extrapolation;
  short extrapolation_type;
  
  double* k_sampling;
  double* logk_sampling; // Of size size_fft_coeff due to the nature of fft
  double logmink;
  double deltalogk;
  
  double k_max;
  double k_max_extr;
  
  int k_size;
  
  /**
   * Window sampling and windows
   * */
  
  double small_log_offset;
  int num_windows;
  int num_window_grid;
  int tw_size; //Length in steps of tau-indices of the window
  double* tw_sampling; //[index_wd*tw_size+index_tw]The tau values sampled for the windows
  double* tw_weights; //Trapezoidal weights for tau window integration
  
  double* exp_integrated_tw_sampling; //exp(integrated_tw_sampling), ONLY in case of logarithmic sampling
  double* integrated_tw_sampling;
  double* integrated_tw_weights;
  int integrated_tw_size;
  
  int ptw_size;
  double* ptw_sampling;
  double* ptw_weights;
  double* ptw_orig_window;
  
  int ptw_integrated_size;
  double* ptw_integrated_sampling;
  double* ptw_integrated_weights;
  
  double** ptw_window;
  double** ptw_dwindow;
  double** ptw_ddwindow;
  
  /**
   * Intxi
   * */
  
  
  double** intxi_real;
  double** intxi_imag;
  
  double** intxi_spline_real;
  double** intxi_spline_imag;
  double** ddintxi_spline_real;
  double** ddintxi_spline_imag;
  
  /**
   * Short pvecback
   * */
  double* short_pvecback;//Used to obtain a(tau) to get z(tau)
  
  
  /**
   * T sampling
   * */
  int t_size;
  double* t_sampling;
  double* t_weights;
  
  int t_spline_size;
  double* t_spline_sampling;
  
  /**
   * Hypergeometric/Bessel integrals coefficients
   * */
  double bias;
  double* nu_real;
  double* nu_imag;
  double bessel_imag_offset;
  
  /**
   * Hypergeometric/Bessel arrays
   * */
  int bi_wanted_samples;
  double*** bi_real; //Bessel integral real 
  double*** bi_imag; //Bessel integral imag
  double*** ddbi_real;
  double*** ddbi_imag;
  int** bi_size;
  double*** bi_sampling;
  //[l*size_fft_coeff+ nu], [t]
  double** bi_max;
  double bi_maximal_t_offset;
  double bessel_recursion_t_size;
  
  /**
   * Experimental
   * */
  int bi_global_size;
  double* bi_global_sampling;
  double* chi_global_sampling;
  int chi_global_size;
  double*** cldd;
  double*** dl;
  /**
   * l sampling
   * */
  int l_size;
  int l_size_recursion;
  double* l_sampling;
  double* l_sampling_recursion;
  
  
  /**
   * Cl's
   * */
  double ** cl;
  double ** ddcl;
  
  
  /**
   * Selection/Window values, and other adjustments
   * */
  double selection_bias[_SELECTION_NUM_MAX_];               /**< light-to-mass bias in the transfer function of density number count */
  double selection_magnification_bias[_SELECTION_NUM_MAX_]; /**< magnification bias in the transfer function of density number count */
  short has_nz_file;     /**< Has dN/dz (selection function) input file? */
  short has_nz_analytic; /**< Use analytic form for dN/dz (selection function) distribution? */
  FileName nz_file_name; /**< dN/dz (selection function) input file name */
  int nz_size;           /**< number of redshift values in input tabulated selection function */
  double * nz_z;         /**< redshift values in input tabulated selection function */
  double * nz_nz;        /**< input tabulated values of selection function */
  double * nz_ddnz;      /**< second derivatives in splined selection function*/
  short has_nz_evo_file;      /**< Has dN/dz (evolution function) input file? */
  short has_nz_evo_analytic;  /**< Use analytic form for dN/dz (evolution function) distribution? */
  FileName nz_evo_file_name;  /**< dN/dz (evolution function) input file name */
  int nz_evo_size;            /**< number of redshift values in input tabulated evolution function */
  double * nz_evo_z;          /**< redshift values in input tabulated evolution function */
  double * nz_evo_nz;         /**< input tabulated values of evolution function */
  double * nz_evo_dlog_nz;    /**< log of tabulated values of evolution function */
  double * nz_evo_dd_dlog_nz; /**< second derivatives in splined log of evolution function */
};
 struct matters_workspace{};
#ifdef __cplusplus
extern "C" {
#endif
  int matter_init(
                  struct precision * ppr,
                  struct background * pba,
                  struct thermo * pth,
                  struct perturbs * ppt,
                  struct primordial * ppm,
                  struct nonlinear * pnl,
                  struct matters * pma
                  );
  int matter_free(
                  struct matters * pls
                  );
  int matter_free_prepare_window(
                  struct matters* pma
                  );
  int matter_obtain_tau_sampling(
                  struct precision* ppr,
                  struct background* pba,
                  struct perturbs * ppt,
                  struct matters * pma
                  );
  int matter_obtain_k_sampling(
                  struct perturbs* ppt,
                  struct matters* pma
                  );
  int matter_obtain_l_sampling(
                  struct precision* ppr,
                  struct perturbs* ppt,
                  struct thermo* pth,
                  struct matters * pma
                  );
  int matter_obtain_indices(
                  struct primordial* ppm,
                  struct perturbs* ppt,
                  struct matters* pma
                  );
  int matter_obtain_coeff_sampling(
                  struct matters * pma
                  );
  int matter_obtain_bessel_integrals(
                  struct matters * pma
                  );   
  int matter_obtain_perturbation_sources(
                  struct background * pba,
                  struct perturbs * ppt,
                  struct nonlinear * pnl,
                  struct matters * pma,
                  double ** sources
                  );
  int matter_obtain_primordial_spectrum(
                  struct perturbs* ppt,
                  struct primordial* ppm,
                  struct matters* pma,
                  double ** prim_spec
                  );
  int matter_obtain_growth_factor(
                  struct matters* pma,
                  double ** sources,
                  double * k_weights
                  );
  int matter_obtain_relative_factor(
                  struct matters* pma,
                  double ** sources,
                  double * k_weights
                  );
  int matter_obtain_growth_factor_k_weights(
                  struct matters* pma,
                  struct perturbs* ppt,
                  double * k_weights
                  );
  int matter_obtain_prepare_windows(
                  struct precision* ppr,
                  struct perturbs* ppt,
                  struct background* pba,
                  struct matters* pma
                  );
  int matter_get_prepared_window_at(
                  struct matters* pma,
                  double tau,
                  int index_ic,
                  int index_radtp,
                  int index_wd,
                  int* last,
                  int derivative_type,
                  double* win_val
                  );
  int matter_obtain_bi_indices(
                  struct matters* pma
                  );
  int matter_spline_growth_factor(
                  struct matters* pma
                  );
  int matter_obtain_nonseparability(
                  struct matters* pma,
                  double ** fft_coeff_real,
                  double ** fft_coeff_imag
                  );
  int matter_spline_cls(
                  struct matters* pma
                  );
  int matter_free_perturbation_sources(
                  struct perturbs * ppt,
                  struct nonlinear * pnl,
                  struct matters * pma,
                  double ** sources
                  ); 
  int matter_free_primordial(
                  struct primordial * ppm,
                  struct matters * pma,
                  double ** prim_spec
                  );
  int matter_free_fft(
                  struct matters * pma,
                  double ** fft_real,
                  double ** fft_imag
                  );
  int matter_FFTlog_perturbation_sources(
                  struct background * pba,
                  struct perturbs * ppt,
                  struct primordial * ppm,
                  struct matters * pma,
                  double ** sampled_sources,
                  double ** prim_spec,
                  double ** fft_coeff_real,
                  double ** fft_coeff_imag
                  );
  int matter_sample_sources(
                  struct perturbs* ppt,
                  struct matters* pma,
                  double** source,
                  double** sampled_source,
                  double* perturbed_k_sampling
                  );
  int matter_extrapolate_sources(
                  struct background* pba,
                  struct precision* ppr,
                  struct perturbs* ppt,
                  struct matters* pma,
                  double** k_extrapolated,
                  double** sources,
                  double** extrapolated_sources,
                  short extrapolation_type
                  );
  int matter_integrate_cl(
                  struct precision* ppr,
                  struct background* pba,
                  struct perturbs * ppt,
                  struct matters* pma,
                  double ** fft_coeff_real,
                  double ** fft_coeff_imag
                  );
  int matter_spline_bessel_integrals_recursion(
                  struct matters * pma
                  );
  int matter_spline_bessel_integrals(
                  struct matters * pma
                  );
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
                  );             
  int matter_interpolate_spline_growing_hunt(
                  double * x_array,
                  int n_lines,
                  double * array,
                  double * array_splined,
                  int n_columns,
                  double x,
                  int * last_index,
                  double * result,
                  int IS_PRINT,
                  ErrorMsg errmsg
                  );
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
                  );
  int matter_estimate_t_max_bessel(
                  struct matters * pma,
                  double l, 
                  double nu_imag,
                  double* t_max_estimate
                  );
  int matter_window_function(
                  struct precision * ppr,
                  struct perturbs * ppt,
                  struct matters* pma,
                  int bin,
                  double z,
                  double * selection);
  int matter_dNdz_analytic(
                  struct matters * pma,
                  double z,
                  double * dNdz,
                  double * dln_dNdz_dz);
  int array_integrate_gauss(
                  double* xarray, 
                  double* warray, 
                  int N,
                  short gauss_type,
                  ErrorMsg err_msg
                  );
  int array_integrate_gauss_limits(
                  double* xarray, 
                  double* warray,
                  double xmin,
                  double xmax, 
                  int N,
                  short gauss_type,
                  ErrorMsg err_msg
                  );
  int matter_derive(
                  double* x_array,
                  double* array,
                  int x_length,
                  double* dx_array,
                  ErrorMsg errmsg
                  );
  int matter_spline_prepare_hunt(
                  double* x_array,
                  int x_size,
                  double x,
                  int* last,
                  ErrorMsg errmsg
                  );
  int matter_resample_growth_factor(
                  struct matters * pma
                  );
  int matter_cl_at_l(
                  struct matters* pma,
                  double l,
                  double * cl_tot,    /* array with argument cl_tot[index_cltp*pma->num_window_grid+index_wd_grid] (must be already allocated) */
                  double ** cl_ic /* array with argument cl_ic[index_ic1_ic2][index_cltp*pma->num_window_grid+index_wd_grid]  (must be already allocated only if several ic's) */
                  );
  int matter_get_bessel_limber(
                  struct matters* pma,
                  int index_l,
                  double** window_bessel_real,
                  double** window_bessel_imag
                  );
  int matter_get_derivative_type(
                  struct matters* pma,
                  int* derivative_type1,
                  int* derivative_type2,
                  int index_radtp1,
                  int index_radtp2
                  );
  int matter_precompute_chit_factors(
                  struct matters* pma,
                  double* tw_local_sampling,
                  int tw_local_size,
                  int index_tilt1_tilt2,
                  int index_wd,
                  short integrate_logarithmically,
                  double* pref_real,
                  double* pref_imag);
  int matter_asymptote(
                  struct precision* ppr, 
                  struct matters* pma,
                  double t, 
                  int index_wd1, 
                  int index_wd2,
                  double* result
                  );
  int matter_spline_hunt(
                  double* x_array,
                  int x_size,
                  double x,
                  int* last,
                  double* h,
                  double* a,
                  double* b,
                  ErrorMsg errmsg
                  );
  int matter_obtain_t_sampling(
                  struct matters* pma
                  );
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
                  );
  int matter_get_bessel_fort_parallel(
                  struct background* pba,
                  struct matters* pma,
                  int index_wd1,
                  int index_wd2,
                  int index_l,
                  double** window_bessel_real,
                  double** window_bessel_imag,
                  short integrate_logarithmically
                  );
  int matter_obtain_time_sampling(
                  struct precision* ppr,
                  struct perturbs* ppt,
                  struct background* pba,
                  struct matters* pma);
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
                  );
  int matter_obtain_bessel_recursion_parallel(struct matters* pma);
  int matter_get_bessel_fort_parallel_integrated(
                  struct background* pba,
                  struct matters* pma,
                  int index_wd1,
                  int index_wd2,
                  int index_l,
                  double** window_bessel_real,
                  double** window_bessel_imag,
                  short integrate_logarithmically
                  );
  int array_integrate_gauss_rescale_limits(double* xarray,double* warray,double* xarrayres,double* warrayres,double xmin,double xmax,int N);
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
                  );
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
                  );
  int matter_obtain_prepare_windows_parallel(
                  struct precision* ppr,
                  struct perturbs* ppt,
                  struct background* pba,
                  struct matters* pma
                  );
  int matter_FFTlog_perturbation_sources_parallel(
                  struct background * pba,
                  struct perturbs * ppt,
                  struct primordial * ppm,
                  struct matters * pma,
                  double ** sampled_sources,
                  double ** prim_spec,
                  double ** fft_coeff_real,
                  double ** fft_coeff_imag
                  );
#ifdef __cplusplus
}
#endif
 #endif
