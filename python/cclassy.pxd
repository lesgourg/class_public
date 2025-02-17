# Bunch of declarations from C to python. The idea here is to define
# only the quantities that will be used, for input, output or
# intermediate manipulation, by the python wrapper. For instance, in
# the precision structure, the only item used here are the error
# message and one parameter used for an error message. That is why
# nothing more is defined from this structure. The rest is internal in
# Class.  If, for whatever reason, you need an other, existing
# parameter from Class, remember to add it inside this cdef.

cdef extern from "class.h":

    cdef char[10] _VERSION_

    ctypedef char FileArg[40]

    ctypedef char* ErrorMsg

    ctypedef char FileName[256]

    cdef enum interpolation_method:
        inter_normal
        inter_growing_closeby

    cdef enum vecback_format:
        short_info
        normal_info
        long_info

    cdef enum linear_or_logarithmic:
        linear
        logarithmic

    cdef enum file_format:
        class_format
        camb_format

    cdef enum non_linear_method:
        nl_none
        nl_halofit
        nl_HMcode

    cdef enum pk_outputs:
        pk_linear
        pk_nonlinear
        pk_numerical_nowiggle
        pk_analytic_nowiggle

    cdef enum out_sigmas:
        out_sigma
        out_sigma_prime
        out_sigma_disp

    cdef struct precision:
        double nonlinear_min_k_max
        ErrorMsg error_message

    cdef struct background:
        short is_allocated
        ErrorMsg error_message
        int bg_size
        int index_bg_ang_distance
        int index_bg_lum_distance
        int index_bg_conf_distance
        int index_bg_a
        int index_bg_H
        int index_bg_D
        int index_bg_f
        int index_bg_Omega_m
        int index_bg_rho_b
        int index_bg_rho_cdm
        int index_bg_rho_ncdm1
        int index_bg_rho_crit
        short has_cdm
        short  has_ncdm
        int N_ncdm
        double T_cmb
        double h
        double H0
        double age
        double conformal_age
        double K
        double * m_ncdm_in_eV
        double Neff
        double Omega0_g
        double Omega0_b
        double Omega0_idr
        double T_idr
        double Omega0_cdm
        double Omega0_idm
        double Omega0_dcdm
        double Omega0_ncdm_tot
        double Omega0_lambda
        double Omega0_fld
        double w0_fld
        double wa_fld
        double cs2_fld
        double Omega0_ur
        double Omega0_dcdmdr
        double Omega0_dr
        double Omega0_scf
        double Omega0_k
        int bt_size
        double Omega0_m
        double Omega0_r
        double Omega0_de
        double a_eq
        double H_eq
        double z_eq
        double tau_eq

    cdef struct thermodynamics:
        short is_allocated
        ErrorMsg error_message
        int th_size
        int index_th_xe
        int index_th_Tb
        double tau_reio
        double z_reio
        double z_rec
        double tau_rec
        double rs_rec
        double ds_rec
        double da_rec
        double z_star
        double tau_star
        double rs_star
        double ds_star
        double ra_star
        double da_star
        double rd_star
        double z_d
        double tau_d
        double ds_d
        double rs_d
        double conf_time_reio
        double YHe
        double n_e
        double a_idm_dr
        double b_idr
        double nindex_idm_dr
        double m_idm
        double cross_idm_g
        double u_idm_g
        double cross_idm_b
        double n_index_idm_b
        int tt_size

    cdef struct perturbations:
        short is_allocated
        ErrorMsg error_message
        short has_scalars
        short has_vectors
        short has_tensors

        short has_density_transfers
        short has_velocity_transfers

        int has_pk_matter
        int l_lss_max

        int store_perturbations
        int k_output_values_num
        double k_output_values[30]
        double k_max_for_pk
        int index_k_output_values[30]
        char scalar_titles[8000]
        char vector_titles[8000]
        char tensor_titles[8000]
        int number_of_scalar_titles
        int number_of_vector_titles
        int number_of_tensor_titles
        int index_md_scalars

        double * scalar_perturbations_data[30]
        double * vector_perturbations_data[30]
        double * tensor_perturbations_data[30]
        int size_scalar_perturbation_data[30]
        int size_vector_perturbation_data[30]
        int size_tensor_perturbation_data[30]

        double * alpha_idm_dr
        double * beta_idr

        # add source functions for comparison
        short has_source_t
        short has_source_p
        short has_source_delta_m
        short has_source_delta_cb
        short has_source_delta_tot
        short has_source_delta_g
        short has_source_delta_b
        short has_source_delta_cdm
        short has_source_delta_idm
        short has_source_delta_idr
        short has_source_delta_dcdm
        short has_source_delta_fld
        short has_source_delta_scf
        short has_source_delta_dr
        short has_source_delta_ur
        short has_source_delta_ncdm
        short has_source_theta_m
        short has_source_theta_cb
        short has_source_theta_tot
        short has_source_theta_g
        short has_source_theta_b
        short has_source_theta_cdm
        short has_source_theta_idm
        short has_source_theta_idr
        short has_source_theta_dcdm
        short has_source_theta_fld
        short has_source_theta_scf
        short has_source_theta_dr
        short has_source_theta_ur
        short has_source_theta_ncdm
        short has_source_phi
        short has_source_phi_prime
        short has_source_phi_plus_psi
        short has_source_psi
        short has_source_h
        short has_source_h_prime
        short has_source_eta
        short has_source_eta_prime
        short has_source_H_T_Nb_prime
        short has_source_k2gamma_Nb

        int index_tp_t0
        int index_tp_t1
        int index_tp_t2
        int index_tp_p
        int index_tp_delta_m
        int index_tp_delta_cb
        int index_tp_delta_tot
        int index_tp_delta_g
        int index_tp_delta_b
        int index_tp_delta_cdm
        int index_tp_delta_idm
        int index_tp_delta_dcdm
        int index_tp_delta_fld
        int index_tp_delta_scf
        int index_tp_delta_dr
        int index_tp_delta_ur
        int index_tp_delta_idr
        int index_tp_delta_ncdm1

        int index_tp_theta_m
        int index_tp_theta_cb
        int index_tp_theta_tot
        int index_tp_theta_g
        int index_tp_theta_b
        int index_tp_theta_cdm
        int index_tp_theta_dcdm
        int index_tp_theta_fld
        int index_tp_theta_scf
        int index_tp_theta_ur
        int index_tp_theta_idr
        int index_tp_theta_idm
        int index_tp_theta_dr
        int index_tp_theta_ncdm1

        int index_tp_phi
        int index_tp_phi_prime
        int index_tp_phi_plus_psi
        int index_tp_psi
        int index_tp_h
        int index_tp_h_prime
        int index_tp_eta
        int index_tp_eta_prime
        int index_tp_H_T_Nb_prime
        int index_tp_k2gamma_Nb


        double *** sources
        double * tau_sampling
        int tau_size
        int k_size_pk
        int * k_size
        double ** k
        int * ic_size
        int index_ic_ad
        int md_size
        int * tp_size
        double * ln_tau
        int ln_tau_size

    cdef struct transfer:
        short is_allocated
        ErrorMsg error_message

    cdef struct primordial:
        short is_allocated
        ErrorMsg error_message
        double k_pivot
        double A_s
        double n_s
        double alpha_s
        double beta_s
        double r
        double n_t
        double alpha_t
        double V0
        double V1
        double V2
        double V3
        double V4
        double f_cdi
        double n_cdi
        double c_ad_cdi
        double n_ad_cdi
        double f_nid
        double n_nid
        double c_ad_nid
        double n_ad_nid
        double f_niv
        double n_niv
        double c_ad_niv
        double n_ad_niv
        double phi_min
        double phi_max
        int lnk_size

    cdef struct harmonic:
        short is_allocated
        ErrorMsg error_message
        int has_tt
        int has_te
        int has_ee
        int has_bb
        int has_pp
        int has_tp
        int has_dd
        int has_td
        int has_ll
        int has_dl
        int has_tl
        int l_max_tot
        int ** l_max_ct
        int ct_size
        int * ic_size
        int * ic_ic_size
        int md_size
        int d_size
        int non_diag
        int index_ct_tt
        int index_ct_te
        int index_ct_ee
        int index_ct_bb
        int index_ct_pp
        int index_ct_tp
        int index_ct_dd
        int index_ct_td
        int index_ct_pd
        int index_ct_ll
        int index_ct_dl
        int index_ct_tl
        int * l_size
        int index_md_scalars

    cdef struct output:
        ErrorMsg error_message

    cdef struct distortions:
        short is_allocated
        double * sd_parameter_table
        int index_type_g
        int index_type_mu
        int index_type_y
        int index_type_PCA
        int type_size
        double * DI
        double * x
        double DI_units
        double x_to_nu
        int has_distortions
        int x_size
        ErrorMsg error_message

    cdef struct lensing:
        short is_allocated
        int has_tt
        int has_ee
        int has_te
        int has_bb
        int has_pp
        int has_tp
        int has_dd
        int has_td
        int has_ll
        int has_dl
        int has_tl
        int index_lt_tt
        int index_lt_te
        int index_lt_ee
        int index_lt_bb
        int index_lt_pp
        int index_lt_tp
        int index_lt_dd
        int index_lt_td
        int index_lt_ll
        int index_lt_dl
        int index_lt_tl
        int * l_max_lt
        int lt_size
        int has_lensed_cls
        int l_lensed_max
        int l_unlensed_max
        ErrorMsg error_message

    cdef struct fourier:
        short is_allocated
        short has_pk_matter
        short has_pk_numerical_nowiggle
        short has_pk_analytic_nowiggle
        int method
        int ic_size
        int ic_ic_size
        int k_size
        int k_size_pk
        int ln_tau_size
        int index_ln_tau_pk
        int tau_size
        int index_tau_min_nl
        double * k
        double * ln_tau
        double * tau
        double ** ln_pk_l
        double ** ln_pk_nl
        double * sigma8
        int has_pk_m
        int has_pk_cb
        int index_pk_m
        int index_pk_cb
        int index_pk_total
        int index_pk_cluster
        ErrorMsg error_message

    cdef struct file_content:
        char * filename
        int size
        FileArg * name
        FileArg * value
        short * read

    void lensing_free(void*)
    void harmonic_free(void*)
    void transfer_free(void*)
    void primordial_free(void*)
    void perturbations_free(void*)
    void thermodynamics_free(void*)
    void background_free(void*)
    void fourier_free(void*)
    void distortions_free(void*)

    cdef int _FAILURE_
    cdef int _FALSE_
    cdef int _TRUE_

    cdef double _Mpc_over_m_
    cdef double _c_
    cdef double _G_
    cdef double _eV_

    int input_read_from_file(void*, void*, void*, void*, void*, void*, void*, void*, void*,
        void*, void*, void*, char*)
    int background_init(void*,void*)
    int thermodynamics_init(void*,void*,void*)
    int perturbations_init(void*,void*,void*,void*)
    int primordial_init(void*,void*,void*)
    int fourier_init(void*,void*,void*,void*,void*,void*)
    int transfer_init(void*,void*,void*,void*,void*,void*)
    int harmonic_init(void*,void*,void*,void*,void*,void*,void*)
    int lensing_init(void*,void*,void*,void*,void*)
    int distortions_init(void*,void*,void*,void*,void*,void*)

    int background_tau_of_z(void* pba, double z,double* tau)
    int background_z_of_tau(void* pba, double tau,double* z)
    int background_at_z(void* pba, double z, int return_format, int inter_mode, int * last_index, double *pvecback)
    int background_at_tau(void* pba, double tau, int return_format, int inter_mode, int * last_index, double *pvecback)
    int background_output_titles(void * pba, char titles[8000])
    int background_output_data(void *pba, int number_of_titles, double *data)

    int thermodynamics_at_z(void * pba, void * pth, double z, int inter_mode, int * last_index, double *pvecback, double *pvecthermo)
    int thermodynamics_output_titles(void * pba, void *pth, char titles[8000])
    int thermodynamics_output_data(void *pba, void *pth, int number_of_titles, double *data)

    int perturbations_output_data_at_z(void *pba,void *ppt, file_format output_format, double z, int number_of_titles, double *data)
    int perturbations_output_data_at_index_tau(void *pba,void *ppt, file_format output_format, int ondex_tau, int number_of_titles, double *data)
    int perturbations_output_data(void *pba,void *ppt, file_format output_format, double * tkfull, int number_of_titles, double *data)
    int perturbations_output_firstline_and_ic_suffix(void *ppt, int index_ic, char first_line[1024], FileName ic_suffix)
    int perturbations_output_titles(void *pba, void *ppt,  file_format output_format, char titles[8000])

    int primordial_output_titles(void * ppt, void *ppm, char titles[8000])
    int primordial_output_data(void *ppt, void *ppm, int number_of_titles, double *data)

    int harmonic_cl_at_l(void* phr,double l,double * cl,double * * cl_md,double * * cl_md_ic)
    int lensing_cl_at_l(void * ple,int l,double * cl_lensed)

    int harmonic_pk_at_z(
        void * pba,
        void * phr,
        int mode,
        double z,
        double * output_tot,
        double * output_ic,
        double * output_cb_tot,
        double * output_cb_ic
        )
    int fourier_pk_at_z(
        void * pba,
        void *pfo,
        int mode,
        int pk_output,
        double z,
        int index_pk,
        double * out_pk,
        double * out_pk_ic
        )

    int harmonic_pk_at_k_and_z(
        void* pba,
        void * ppm,
        void * phr,
        double k,
        double z,
        double * pk,
        double * pk_ic,
        double * pk_cb,
        double * pk_cb_ic)

    int harmonic_pk_nl_at_k_and_z(
        void* pba,
        void * ppm,
        void * phr,
        double k,
        double z,
        double * pk,
        double * pk_cb)

    int harmonic_pk_nl_at_z(
        void * pba,
        void * phr,
        int mode,
        double z,
        double * output_tot,
        double * output_cb_tot)

    int fourier_pk_at_k_and_z(
        void * pba,
        void * ppm,
        void * pfo,
        int pk_output,
        double k,
        double z,
        int index_pk,
        double * out_pk,
        double * out_pk_ic)

    int fourier_pk_tilt_at_k_and_z(
        void * pba,
        void * ppm,
        void * pfo,
        int pk_output,
        double k,
        double z,
        int index_pk,
        double * pk_tilt)

    int fourier_sigmas_at_z(
        void * ppr,
        void * pba,
        void * pfo,
        double R,
        double z,
        int index_pk,
        int sigma_output,
        double * result)

    int fourier_pks_at_kvec_and_zvec(
        void * pba,
        void * pfo,
        int pk_output,
        double * kvec,
        int kvec_size,
        double * zvec,
        int zvec_size,
        double * out_pk,
        double * out_pk_cb)

    int fourier_k_nl_at_z(void* pba, void* pfo, double z, double* k_nl, double* k_nl_cb)

    int harmonic_firstline_and_ic_suffix(void *ppt, int index_ic, char first_line[1024], FileName ic_suffix)

    int harmonic_fast_pk_at_kvec_and_zvec(
                  void * pba,
                  void * phr,
                  double * kvec,
                  int kvec_size,
                  double * zvec,
                  int zvec_size,
                  double * pk_tot_out,
                  double * pk_cb_tot_out,
                  int nonlinear)
