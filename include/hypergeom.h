/**
 * Definitions of the basic bessel integrals I_l(nu,t) = 4*pi integral 0 to inf of v^(nu-1) j_l(v) j_l(vt)
 * Implemented as a simple version taken from the corresponding paper
 * and as a self-derived version (transform)
 */
 #ifndef _HYPERGEOM_DEFINED
 #define _HYPERGEOM_DEFINED

#define T_MIN_TAYLOR 1e-2
#define T_MIN_INVERSE_TAYLOR 1e-4
//#define T_MIN_INVERSE_TAYLOR 1e-10
//#define BESSEL_EPSILON 1e-8
#define BESSEL_EPSILON 1e-10
#define T_SWITCH_FORWARD_RECURSION_SIMPLE 0.9
//0.9999
//0.94
//void bessel_integral(double l, double nu_real, double nu_imag, double t, double* res_real, double* res_imag);
void bessel_integral_transform(double l, double nu_real, double nu_imag, double t, double* res_real, double* res_imag);
//void bessel_integral_l0(double nu_real, double nu_imag, double t, double* res_real, double* res_imag);
//void bessel_integral_l1(double nu_real, double nu_imag, double t, double* res_real, double* res_imag);
void bessel_integral_recursion_initial_abs(int l_max,double nu_real,double nu_imag,double* abi_real,double* abi_imag,double* initial_abs);
int bessel_integral_recursion_complicated(int l_max,int l_recursion_max,double nu_real, double nu_imag, double t,double bi_allowed_error,double* bi_real,double* bi_imag,double* max_t,double* initial_abs,ErrorMsg errmsg);
void bessel_integral_recursion_taylor(int l_max,double nu_real,double nu_imag,double t,double* max_t,double* initial_abs,double* bi_real,double* bi_imag);
int bessel_integral_recursion_backward_simple(
                                             int l_max,
                                             int l_recursion_max, 
                                             double nu_real,
                                             double nu_imag,
                                             double t,
                                             double nu_fraction,
                                             int bessel_recursion_backward_min_l_step_high_nu,
                                             int bessel_recursion_backward_min_l_step_low_nu,
                                             int bessel_recursion_backward_max_l_step,
                                             double* abi_real,
                                             double* abi_imag,
                                             double* max_t,
                                             double* initial_abs,
                                             ErrorMsg errmsg
                                            );
int bessel_integral_recursion_backward_simple_safe(
                                             int l_max,
                                             int l_recursion_max, 
                                             double nu_real,
                                             double nu_imag,
                                             double t,
                                             double nu_fraction,
                                             int bessel_recursion_backward_min_l_step_high_nu,
                                             int bessel_recursion_backward_min_l_step_low_nu,
                                             int bessel_recursion_backward_max_l_step,
                                             double* abi_real,
                                             double* abi_imag,
                                             double* max_t,
                                             double* initial_abs,
                                             short* overflow_flag,
                                             ErrorMsg errmsg
                                            );
int bessel_integral_recursion_forward_simple(
                                             int l_max, 
                                             double nu_real,
                                             double nu_imag,
                                             double t,
                                             double nu_fraction,
                                             int bessel_recursion_forward_min_l_step_higj_nu,
                                             int bessel_recursion_forward_min_l_step_low_nu,
                                             int bessel_recursion_forward_max_l_step,
                                             double* abi_real,
                                             double* abi_imag,
                                             double* max_t,
                                             double* initial_abs,
                                             ErrorMsg errmsg
                                             );
int bessel_integral_recursion_inverse_self(int l_max,
                                           double nu_real,
                                           double nu_imag,
                                           double t,
                                           double* abi_real,
                                           double* abi_imag,
                                           double* max_t,
                                           double* initial_abs,
                                           ErrorMsg errmsg);
 #endif
